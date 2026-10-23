#!/usr/bin/env python3
"""
Build the detectability table from synthetic spectra.

This is the step that produced data/ce_detectability_table.csv.gz. It is included
so the method is reproducible and auditable; the full 1.6 M-spectrum grid is too
large to distribute, so the shipped table is its output.

Input, either:
  * a directory of synth_teff*_logg*_feh*_cefe*.txt files (the raw grid), or
  * a FITS file written by condense_grid.py (the grid reduced to the Ce window)

For each (Teff, logg, [Fe/H]) node it steps through [Ce/Fe] and computes, between
consecutive abundances:

    g_i      = (f_i(X+D) - f_i(X)) / D                    per-pixel gradient
    sigma    = 1 / (SNR * sqrt(sum_i g_i^2))              Cramer-Rao bound
    EW       = integral of the line against a local linear continuum
    sigma_W  = (1.5/SNR) * sqrt(FWHM * dx)                Cayrel (1988)

and flags crb_pass (sigma < 0.1 dex), ew_pass (EW >= 3 sigma_W) and their AND.

    python3 compute_gradients.py /path/to/synthetic_spectra2 -o table.csv
    python3 compute_gradients.py ce_grid_example.fits -o example_table.csv
"""

import os
import re
import glob
import argparse
import itertools
import numpy as np
import pandas as pd

CE_WAVE_MIN  = 851.05     # nm, gradient window (paper Appendix B)
CE_WAVE_MAX  = 851.25
ASSUMED_SNR  = 50
TARGET_SIGMA = 0.1        # dex
CE_FWHM_NM   = float(np.mean([0.0863, 0.0895, 0.0988, 0.0998]))  # = 0.0936 nm
                          # measured median Ce-line FWHM (the paper rounds it
                          # to 0.094; 0.0936 is the value used for the table)
RVS_DX_NM    = 0.03       # nm, RVS pixel size -- NOT the fine synthetic sampling
N_SIGMA      = 3
N_EDGE_FRAC  = 0.025      # fraction of the window used per side as continuum

WAVE0_SYNTH = 851.0
DX_SYNTH    = 0.001
N_HEADER    = 1
PAD         = 3

FNAME_RE = re.compile(
    r'teff(\d+)_logg([+-]?\d+\.?\d*)_feh([+-]?\d+\.?\d*)_cefe([+-]?\d+\.?\d*)\.txt$')

EW_LIMIT = N_SIGMA * (1.5 / ASSUMED_SNR) * np.sqrt(CE_FWHM_NM * RVS_DX_NM)


def measure_ew(wave, flux, n_edge):
    """EW (nm) against a local linear continuum taken from the window edges, so a
    sloped pseudo-continuum or a blended neighbour is not counted as Ce."""
    flux = np.atleast_2d(flux)
    cl = np.median(flux[..., :n_edge], axis=-1)
    cr = np.median(flux[..., -n_edge:], axis=-1)
    cont = cl[..., None] + (cr - cl)[..., None] * np.linspace(0, 1, flux.shape[-1])
    cont = np.where(cont <= 0, 1.0, cont)
    return np.trapezoid(np.clip(1.0 - flux / cont, 0.0, None), wave, axis=-1)


def load_from_fits(path):
    from astropy.io import fits
    with fits.open(path) as h:
        flux = np.array(h[0].data.T, float)
        w0, dw = np.asarray(h[1].data, float).ravel()
        t = h[2].data
        par = pd.DataFrame({'teff': t['teff'].astype(float),
                            'logg': np.round(t['logg'].astype(float), 2),
                            'feh':  np.round(t['feh'].astype(float), 2),
                            'cefe': np.round(t['cefe'].astype(float), 2)})
    return par, flux, w0 + dw * np.arange(flux.shape[1])


def load_from_txt(directory):
    files = sorted(glob.glob(os.path.join(directory, 'synth_*.txt')))
    if not files:
        raise SystemExit(f"no synth_*.txt in {directory}")
    lo = N_HEADER + max(int((CE_WAVE_MIN - WAVE0_SYNTH) / DX_SYNTH) - PAD, 0)
    hi = N_HEADER + int((CE_WAVE_MAX - WAVE0_SYNTH) / DX_SYNTH) + PAD

    rows, spectra, wave = [], [], None
    for k, p in enumerate(files):
        m = FNAME_RE.search(os.path.basename(p))
        if m is None:
            continue
        with open(p) as fh:
            chunk = list(itertools.islice(fh, lo, hi + 1))
        try:
            arr = np.array([ln.split('\t')[:2] for ln in chunk if ln.strip()], float)
        except ValueError:
            continue
        msk = (arr[:, 0] >= CE_WAVE_MIN) & (arr[:, 0] <= CE_WAVE_MAX)
        if wave is None:
            wave = arr[msk, 0]
        if msk.sum() != len(wave):
            continue
        spectra.append(arr[msk, 1])
        rows.append((float(m.group(1)), round(float(m.group(2)), 2),
                     round(float(m.group(3)), 2), round(float(m.group(4)), 2)))
        if (k + 1) % 100000 == 0:
            print(f"  {k+1}/{len(files)}", flush=True)
    par = pd.DataFrame(rows, columns=['teff', 'logg', 'feh', 'cefe'])
    return par, np.array(spectra), wave


def compute(par, flux, wave):
    n_edge = max(1, int(round(N_EDGE_FRAC * len(wave))))
    print(f"window {wave[0]:.3f}-{wave[-1]:.3f} nm, {len(wave)} pixels, "
          f"n_edge={n_edge}")

    par = par.assign(_node=pd.factorize(pd.Series(list(zip(par.teff, par.logg, par.feh))))[0])
    order = np.lexsort((par.cefe.values, par._node.values))
    par, flux = par.iloc[order].reset_index(drop=True), flux[order]

    same = par._node.values[1:] == par._node.values[:-1]
    lo = np.where(same)[0]
    hi = lo + 1

    delta = par.cefe.values[hi] - par.cefe.values[lo]
    grad = (flux[hi] - flux[lo]) / delta[:, None]
    grad_rss = np.sqrt(np.sum(grad ** 2, axis=1))
    with np.errstate(divide='ignore'):
        sigma = 1.0 / (ASSUMED_SNR * grad_rss)
    ew = measure_ew(wave, flux[hi], n_edge)

    out = pd.DataFrame({
        'teff': par.teff.values[hi], 'logg': par.logg.values[hi],
        'feh': par.feh.values[hi],
        'cefe_mid': np.round((par.cefe.values[hi] + par.cefe.values[lo]) / 2, 3),
        'delta_cefe': np.round(delta, 3),
        'sigma_cefe': sigma, 'grad_rss': grad_rss,
        'max_gradient': np.max(np.abs(grad), axis=1),
        'ew_ce': ew, 'ew_limit': EW_LIMIT,
    })
    out['crb_pass'] = out.sigma_cefe < TARGET_SIGMA
    out['ew_pass'] = out.ew_ce >= out.ew_limit
    out['is_useful'] = out.crb_pass & out.ew_pass
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('input', help='spectra directory, or a condense_grid.py FITS')
    ap.add_argument('-o', '--output', default='ce_detectability_table.csv')
    args = ap.parse_args()

    if os.path.isdir(args.input):
        par, flux, wave = load_from_txt(args.input)
    else:
        par, flux, wave = load_from_fits(args.input)
    print(f"{len(par)} spectra, "
          f"{par.groupby(['teff','logg','feh']).ngroups} (Teff, logg, [Fe/H]) nodes")

    tab = compute(par, flux, wave)
    comp = dict(compression='gzip') if args.output.endswith('.gz') else {}
    tab.to_csv(args.output, index=False, float_format='%.6g', **comp)
    print(f"wrote {args.output}: {len(tab)} rows, "
          f"detectable {tab.is_useful.sum()} ({100*tab.is_useful.mean():.1f}%)")
    print(f"  EW only  : {100*tab.ew_pass.mean():.1f}%")
    print(f"  CRB only : {100*tab.crb_pass.mean():.1f}%")


if __name__ == '__main__':
    main()
