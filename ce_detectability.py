#!/usr/bin/env python3
"""
Ce detectability in Gaia RVS spectra.

Given a star's (Teff, logg, [Fe/H], [Ce/Fe]) and an assumed S/N, report whether
the Ce II feature at 851.1 nm is detectable, following Vitali et al.,
"Detectability limits of Ce abundances in Gaia RVS".

Two criteria, reported SEPARATELY as well as combined:

  crb_pass   Cramer-Rao bound:  sigma_CRB = 1/(SNR*sqrt(sum_i g_i^2)) < 0.1 dex
             -> the spectrum carries enough Fisher information about [Ce/Fe]
  ew_pass    Cayrel (1988):     EW >= 3*sigma_W,  sigma_W = (1.5/SNR)*sqrt(FWHM*dx)
             -> the line is individually visible above the noise
  detectable crb_pass AND ew_pass

They are kept apart on purpose. The CRB adds signal in quadrature across pixels,
so a line can be "statistically detectable" while its core sits well below the
per-pixel noise -- a regime very sensitive to continuum placement and blending.
The Cayrel test is the reality check.

Both limits scale as 1/SNR, so any S/N is obtained by rescaling the tabulated
S/N = 50 reference.

Usage
-----
    python3 ce_detectability.py --star 4500 1.8 -0.2 0.4 --snr 50
    python3 ce_detectability.py --input stars.csv --output flagged.csv --snr 50

    from ce_detectability import CeDetectability
    ce  = CeDetectability()
    out = ce.classify(df.teff, df.logg, df.feh, df.cefe, snr=df.snr)
    out[['crb_pass', 'ew_pass', 'detectable']]

Note on [alpha/Fe]
------------------
The grid behind this tool fixes [alpha/Fe] to the standard MARCS
metallicity-dependent relation (+0.4 for [Fe/H] <= -1, falling linearly to 0.0
at [Fe/H] >= 0), as described in Appendix A of the paper. It is therefore NOT a
free axis and cannot be supplied: every answer assumes that relation.
"""

import os
import argparse
import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

HERE = os.path.dirname(os.path.abspath(__file__))

TABLE_SNR    = 50      # S/N at which the table was computed
TARGET_SIGMA = 0.1     # dex, CRB precision threshold
CE_FWHM_NM   = 0.0936  # nm, measured median Ce line FWHM (paper rounds to 0.094)
RVS_DX_NM    = 0.03    # nm, RVS pixel size
N_SIGMA      = 3       # Cayrel detection significance

# grid steps, so "nearest" is measured in comparable units on every axis
STEP = {'teff': 25.0, 'logg': 0.2, 'feh': 0.1, 'cefe': 0.1}
AXES = ('teff', 'logg', 'feh', 'cefe')

TABLE_CANDIDATES = [
    os.environ.get('CE_DETECTABILITY_TABLE'),
    os.path.join(HERE, 'data', 'ce_detectability_table.csv.gz'),
    os.path.join(HERE, 'ce_detectability_table.csv.gz'),
]


def marcs_alpha(feh):
    """The [alpha/Fe] implicitly assumed at a given [Fe/H] (MARCS standard)."""
    return np.round(np.clip(-0.4 * np.asarray(feh, float), 0.0, 0.4), 1)


def cayrel_sigma_W(snr, fwhm_nm=CE_FWHM_NM, dx_nm=RVS_DX_NM):
    """Cayrel (1988) uncertainty on the equivalent width, in nm."""
    return (1.5 / snr) * np.sqrt(fwhm_nm * dx_nm)


def _find_table(csv=None):
    for c in ([csv] if csv else []) + TABLE_CANDIDATES:
        if c and os.path.exists(c):
            return c
    raise FileNotFoundError(
        "No detectability table found. Expected one of:\n  " +
        "\n  ".join(str(c) for c in TABLE_CANDIDATES if c) +
        "\nPass one with --table, or set $CE_DETECTABILITY_TABLE.")


class CeDetectability:
    """Nearest-node lookup of the Ce detectability criteria."""

    def __init__(self, csv=None, target_sigma=TARGET_SIGMA, n_sigma=N_SIGMA,
                 fwhm=CE_FWHM_NM, dx=RVS_DX_NM):
        # The thresholds are the user's to choose. The table stores the measured
        # quantities (gradient, EW) plus the limits at the published settings;
        # changing a threshold only rescales those limits, so no re-run of the
        # spectra is needed. Changing the WINDOW or the continuum definition does
        # need a re-run -- that is what compute_gradients.py is for.
        self.target_sigma = target_sigma
        self.n_sigma = n_sigma
        # sigma_W scales as sqrt(FWHM*dx); the table was built at the defaults
        self.ew_scale = ((n_sigma / N_SIGMA) *
                         np.sqrt((fwhm * dx) / (CE_FWHM_NM * RVS_DX_NM)))
        self.csv = _find_table(csv)
        t = pd.read_csv(self.csv)
        self.table = t.rename(columns={'cefe_mid': 'cefe'}).reset_index(drop=True)
        pts = np.column_stack([self.table[k].values / STEP[k] for k in AXES])
        self._tree = cKDTree(pts)
        self.ranges = {k: (self.table[k].min(), self.table[k].max()) for k in AXES}

    def classify(self, teff, logg, feh, cefe, snr=TABLE_SNR):
        """One row per input star; see README for the column meanings."""
        q = pd.DataFrame({'teff': np.atleast_1d(teff).astype(float),
                          'logg': np.atleast_1d(logg).astype(float),
                          'feh':  np.atleast_1d(feh).astype(float),
                          'cefe': np.atleast_1d(cefe).astype(float)})
        snr = np.broadcast_to(np.atleast_1d(np.asarray(snr, float)), (len(q),)).copy()

        pts = np.column_stack([q[k].values / STEP[k] for k in AXES])
        dist, idx = self._tree.query(pts, k=1)
        m = self.table.iloc[idx].reset_index(drop=True)

        scale   = TABLE_SNR / snr               # both limits scale as 1/SNR
        sigma   = m.sigma_cefe.values * scale
        ewlim   = m.ew_limit.values * scale * self.ew_scale
        sigma_W = ewlim / self.n_sigma

        in_range = np.ones(len(q), bool)
        for k, (lo, hi) in self.ranges.items():
            in_range &= (q[k].values >= lo) & (q[k].values <= hi)

        out = pd.DataFrame({
            'teff_in': q.teff, 'logg_in': q.logg,
            'feh_in': q.feh, 'cefe_in': q.cefe, 'snr': snr,
            'alpha_assumed': marcs_alpha(q.feh.values),
            'teff_grid': m.teff.values, 'logg_grid': m.logg.values,
            'feh_grid': m.feh.values, 'cefe_grid': m.cefe.values,
            'grid_distance': dist,
            'in_grid_range': in_range,
            'sigma_cefe': sigma,
            'grad_rss': m.grad_rss.values,
            'ew_ce': m.ew_ce.values,
            'ew_limit': ewlim,
            'ew_snr': m.ew_ce.values / sigma_W,
            'delta_cefe': m.delta_cefe.values,
        })
        out['crb_pass'] = out.sigma_cefe < self.target_sigma
        out['ew_pass'] = out.ew_ce >= out.ew_limit
        out['detectable'] = out.crb_pass & out.ew_pass
        return out


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--star', nargs=4, type=float,
                    metavar=('TEFF', 'LOGG', 'FEH', 'CEFE'))
    ap.add_argument('--input', help='CSV with columns teff, logg, feh, cefe (+ optional snr)')
    ap.add_argument('--output', help='where to write the flagged catalogue')
    ap.add_argument('--snr', type=float, default=TABLE_SNR)
    ap.add_argument('--table', help='path to the detectability table')

    # the thresholds are choices; the defaults are the paper's
    t = ap.add_argument_group('thresholds (defaults = paper values)')
    t.add_argument('--target-sigma', type=float, default=TARGET_SIGMA,
                   help=f'CRB precision threshold, dex (default {TARGET_SIGMA})')
    t.add_argument('--n-sigma', type=float, default=N_SIGMA,
                   help=f'EW detection significance (default {N_SIGMA})')
    t.add_argument('--fwhm', type=float, default=CE_FWHM_NM,
                   help=f'line FWHM, nm (default {CE_FWHM_NM})')
    t.add_argument('--dx', type=float, default=RVS_DX_NM,
                   help=f'instrument pixel, nm (default {RVS_DX_NM})')
    args = ap.parse_args()

    if not (args.star or args.input):
        ap.error('give either --star or --input')

    ce = CeDetectability(csv=args.table, target_sigma=args.target_sigma,
                         n_sigma=args.n_sigma, fwhm=args.fwhm, dx=args.dx)
    print(f"table: {ce.csv}")
    print(f"grid : {len(ce.table)} entries, "
          f"{ce.table.groupby(['teff','logg','feh']).ngroups} (Teff, logg, [Fe/H]) nodes")

    if args.star:
        r = ce.classify(*args.star, snr=args.snr).iloc[0]
        print(f"\nrequested : Teff={r.teff_in:.0f}  logg={r.logg_in:.2f}  "
              f"[Fe/H]={r.feh_in:+.2f}  [Ce/Fe]={r.cefe_in:+.2f}  S/N={r.snr:.0f}")
        print(f"matched   : Teff={r.teff_grid:.0f}  logg={r.logg_grid:.2f}  "
              f"[Fe/H]={r.feh_grid:+.2f}  [Ce/Fe]={r.cefe_grid:+.2f}  "
              f"(distance {r.grid_distance:.2f} cells"
              f"{'' if r.in_grid_range else ', OUTSIDE grid range'})")
        print(f"            [alpha/Fe] = {r.alpha_assumed:+.1f} (MARCS standard, fixed by the grid)")
        print(f"CRB       : sigma = {r.sigma_cefe:.3f} dex  -> "
              f"{'PASS' if r.crb_pass else 'FAIL'} (< {args.target_sigma} dex)")
        print(f"EW        : {r.ew_snr:.1f} sigma            -> "
              f"{'PASS' if r.ew_pass else 'FAIL'} (>= {args.n_sigma:g} sigma)")
        print(f"==> {'DETECTABLE' if r.detectable else 'NOT detectable'}")

    if args.input:
        src = pd.read_csv(args.input)
        snr = src['snr'].values if 'snr' in src.columns else args.snr
        res = ce.classify(src.teff, src.logg, src.feh, src.cefe, snr=snr)
        res = pd.concat([src.reset_index(drop=True),
                         res.drop(columns=['teff_in', 'logg_in', 'feh_in', 'cefe_in'])],
                        axis=1)
        out = args.output or args.input.replace('.csv', '_cedetect.csv')
        res.to_csv(out, index=False, float_format='%.6g')
        print(f"\n{len(res)} stars -> {out}")
        print(f"  CRB only (Fisher) pass : {res.crb_pass.sum():6d} ({100*res.crb_pass.mean():.1f}%)")
        print(f"  EW  only (Cayrel) pass : {res.ew_pass.sum():6d} ({100*res.ew_pass.mean():.1f}%)")
        print(f"  detectable (CRB + EW)  : {res.detectable.sum():6d} ({100*res.detectable.mean():.1f}%)")
        if (~res.in_grid_range).any():
            print(f"  WARNING: {int((~res.in_grid_range).sum())} stars outside the grid "
                  f"range (snapped to the nearest edge node)")


if __name__ == '__main__':
    main()
