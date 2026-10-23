#!/usr/bin/env python3
"""
Condense the synthetic Ce grid to the analysis window.

Reads  synth_teff*_logg*_feh*_cefe*.txt  (1 header line, tab-separated
wavelength/flux, fixed 0.001 nm grid starting at 851.0 nm) and writes ONE FITS
file holding only the Ce window, in the same compact layout as the externally
provided Ultimate_Ce_0.fits:

    HDU0 : float32 image, shape (npix, nspectra)   -- normalised flux
    HDU1 : float64 [[lambda_0], [dlambda]]         -- nm
    HDU2 : binary table                            -- teff, logg, feh, cefe

Nothing the EW / gradient / CRB analysis uses is discarded: it only ever reads
this window. 37 GB of text becomes ~1.3 GB of float32.

Run it on the machine holding the spectra:

    python3 condense_grid.py /path/to/synthetic_spectra2 -o ce_grid_window.fits

Options:
    --wmin / --wmax   window in nm (default 851.05 - 851.25, as in the paper)
    --dtype           float32 (default) or float16 to halve the size again
    --limit N         only process N files, for a quick test
"""

import os
import re
import sys
import glob
import argparse
import itertools
import numpy as np
from astropy.io import fits

WAVE0_SYNTH = 851.0     # nm, wavelength of the first data row
DX_SYNTH    = 0.001     # nm per row
N_HEADER    = 1         # header lines to skip
PAD         = 3         # safety rows either side when slicing

FNAME_RE = re.compile(
    r'teff(\d+)_logg([+-]?\d+\.?\d*)_feh([+-]?\d+\.?\d*)_cefe([+-]?\d+\.?\d*)\.txt$')


def line_range(wmin, wmax):
    lo = int((wmin - WAVE0_SYNTH) / DX_SYNTH) - PAD
    hi = int((wmax - WAVE0_SYNTH) / DX_SYNTH) + PAD
    return N_HEADER + max(lo, 0), N_HEADER + hi


def read_window(path, lo, hi, wmin, wmax):
    """Return (wave, flux) for the window, or None if the file is unusable."""
    with open(path) as fh:
        chunk = list(itertools.islice(fh, lo, hi + 1))
    if not chunk:
        return None
    try:
        arr = np.array([ln.split('\t')[:2] for ln in chunk if ln.strip()], float)
    except ValueError:
        return None
    w, f = arr[:, 0], arr[:, 1]
    if w.min() > wmin + 0.01 or w.max() < wmax - 0.01:
        return None                      # fixed-grid assumption broken
    m = (w >= wmin) & (w <= wmax)
    return w[m], f[m]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('synth_dir')
    ap.add_argument('-o', '--output', default='ce_grid_window.fits')
    ap.add_argument('--wmin', type=float, default=851.05)
    ap.add_argument('--wmax', type=float, default=851.25)
    ap.add_argument('--dtype', default='float32', choices=['float32', 'float16'])
    ap.add_argument('--limit', type=int, default=None,
                    help='only the first N files, for a smoke test')
    ap.add_argument('--sample-nodes', type=int, default=None,
                    help='keep only N (Teff, logg, [Fe/H]) nodes, spread evenly '
                         'over the grid, but ALL their [Ce/Fe] values. Use this '
                         'to build a small example grid for a repository: the '
                         'gradient is a finite difference in [Ce/Fe], so a node '
                         'must keep its full [Ce/Fe] sequence to stay usable.')
    ap.add_argument('--seed', type=int, default=0)
    args = ap.parse_args()

    files = sorted(glob.glob(os.path.join(args.synth_dir, 'synth_*.txt')))
    if not files:
        sys.exit(f"no synth_*.txt found in {args.synth_dir}")
    print(f"{len(files)} files found")

    if args.sample_nodes:
        # group by (teff, logg, feh) and keep whole nodes, so every retained
        # node still has its full [Ce/Fe] sequence and the gradient is defined
        nodes = {}
        for f in files:
            m = FNAME_RE.search(os.path.basename(f))
            if m:
                nodes.setdefault((m.group(1), m.group(2), m.group(3)), []).append(f)
        keys = sorted(nodes)
        if args.sample_nodes < len(keys):
            idx = np.linspace(0, len(keys) - 1, args.sample_nodes).round().astype(int)
            keys = [keys[i] for i in sorted(set(idx))]
        files = sorted(f for k in keys for f in nodes[k])
        print(f"sampled {len(keys)} nodes -> {len(files)} spectra "
              f"({len(files)/max(len(keys),1):.0f} [Ce/Fe] values per node)")

    if args.limit:
        files = files[:args.limit]

    lo, hi = line_range(args.wmin, args.wmax)

    # first good file fixes the wavelength grid and the array width
    wave = None
    for p in files:
        r = read_window(p, lo, hi, args.wmin, args.wmax)
        if r is not None:
            wave = r[0]
            break
    if wave is None:
        sys.exit("could not read a usable spectrum")
    npix = len(wave)
    print(f"window {wave[0]:.3f}-{wave[-1]:.3f} nm, {npix} pixels")
    print(f"estimated output: {len(files)*npix*np.dtype(args.dtype).itemsize/1e9:.2f} GB")

    flux = np.zeros((len(files), npix), dtype=args.dtype)
    teff = np.zeros(len(files), 'float32')
    logg = np.zeros(len(files), 'float32')
    feh  = np.zeros(len(files), 'float32')
    cefe = np.zeros(len(files), 'float32')

    bad, n = [], 0
    for i, p in enumerate(files):
        m = FNAME_RE.search(os.path.basename(p))
        r = read_window(p, lo, hi, args.wmin, args.wmax)
        if m is None or r is None or len(r[1]) != npix:
            bad.append(p)
            continue
        flux[n] = r[1]
        teff[n], logg[n], feh[n], cefe[n] = (float(m.group(1)), float(m.group(2)),
                                             float(m.group(3)), float(m.group(4)))
        n += 1
        if (i + 1) % 50000 == 0:
            print(f"  {i+1}/{len(files)}", flush=True)

    flux, teff, logg, feh, cefe = flux[:n], teff[:n], logg[:n], feh[:n], cefe[:n]
    print(f"kept {n} spectra, skipped {len(bad)}")
    if bad[:5]:
        print("  first skipped:", *[os.path.basename(b) for b in bad[:5]], sep='\n    ')

    hdul = fits.HDUList([
        fits.PrimaryHDU(flux.T.astype(args.dtype)),          # (npix, nspectra)
        fits.ImageHDU(np.array([[wave[0]], [DX_SYNTH]], 'float64')),
        fits.BinTableHDU.from_columns([
            fits.Column(name='teff', format='E', array=teff),
            fits.Column(name='logg', format='E', array=logg),
            fits.Column(name='feh',  format='E', array=feh),
            fits.Column(name='cefe', format='E', array=cefe)]),
    ])
    hdul[0].header['WMIN'] = (float(wave[0]), 'nm, first pixel')
    hdul[0].header['DWAVE'] = (DX_SYNTH, 'nm per pixel')
    hdul[0].header['NSPEC'] = (n, 'number of spectra')
    hdul.writeto(args.output, overwrite=True)
    print(f"wrote {args.output} ({os.path.getsize(args.output)/1e9:.2f} GB)")


if __name__ == '__main__':
    main()
