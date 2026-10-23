#!/usr/bin/env python3
"""
Ce detectability maps in the Kiel diagram.

Reproduces Fig. 2 and Fig. C.1 of Vitali et al., "Detectability limits of Ce
abundances in Gaia RVS": the Ce spectral sensitivity across (Teff, logg), masked
to where the combined CRB + EW criterion is satisfied, with contours at several
S/N.

    python3 plot_detectability.py                 # both figures, default S/N set
    python3 plot_detectability.py --snr 25 50 200
    python3 plot_detectability.py --outdir figures
"""

import os
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from scipy.ndimage import gaussian_filter

from ce_detectability import CeDetectability, TABLE_SNR, TARGET_SIGMA

SNR_STYLE = {25: ('darkorange', '-'), 50: ('black', '--'), 200: ('hotpink', '-')}
DEFAULT_SNR = [50, 25, 200]


def select_cefe(df, cefe, atol=0.03):
    mids = np.unique(df.cefe.values)
    return df[np.isclose(df.cefe, mids[np.argmin(np.abs(mids - cefe))], atol=atol)]


def build_grid(df_bin, snr_list, base_snr=TABLE_SNR):
    """(Teff, logg) maps of the gradient and of detectability at each S/N."""
    grad = df_bin.groupby(['logg', 'teff'])['max_gradient'].mean().unstack()
    logg_vals, teff_vals = grad.index.values, grad.columns.values
    masks = {}
    for snr in snr_list:
        d = df_bin
        det = ((d.sigma_cefe * (base_snr / snr) < TARGET_SIGMA) &
               (d.ew_ce >= d.ew_limit * (base_snr / snr)))
        masks[snr] = (d.assign(det=det).groupby(['logg', 'teff'])['det'].max()
                       .unstack().reindex(index=logg_vals, columns=teff_vals).values)
    return teff_vals, logg_vals, grad.values, masks


def draw_panel(ax, teff, logg, grad, masks, vmin, vmax, base_snr, smooth=1.0):
    det = np.nan_to_num(masks[base_snr], nan=0.0).astype(bool)
    Z = np.where(det, grad, np.nan)
    if smooth:
        Z = np.where(det, gaussian_filter(np.nan_to_num(Z), sigma=smooth), np.nan)
    im = ax.imshow(Z, aspect='auto', origin='upper',
                   extent=[teff.min(), teff.max(), logg.max(), logg.min()],
                   cmap='Blues', vmin=vmin, vmax=vmax)
    ax.invert_xaxis()
    TT, GG = np.meshgrid(teff, logg)
    for snr, M in masks.items():
        M = np.nan_to_num(M, nan=0.0).astype(float)
        if M.sum() == 0:
            continue
        colour, ls = SNR_STYLE.get(snr, ('grey', '-'))
        ax.contour(TT, GG, M, levels=[0.5], colors=colour, linestyles=ls, linewidths=1.6)
    return im


def figure(df, cefe_values, feh_bins, feh_labels, snr_list, ylim, fname, figsize):
    nrow, ncol = len(cefe_values), len(feh_bins)
    panels, labels = [], []
    for cf in cefe_values:
        sel = select_cefe(df, cf)
        labels.append(np.unique(sel.cefe.values)[0])
        row = []
        for i, (lo, hi) in enumerate(feh_bins):
            last = (i == ncol - 1)
            b = sel[(sel.feh >= lo) & ((sel.feh <= hi) if last else (sel.feh < hi))]
            row.append(None if b.empty else build_grid(b, snr_list))
        panels.append(row)

    allg = np.concatenate([p[2].ravel() for r in panels for p in r if p is not None])
    allg = allg[np.isfinite(allg)]
    vmin, vmax = np.quantile(allg, [0.01, 0.98])

    fig = plt.figure(figsize=figsize)
    gs = GridSpec(nrow, ncol + 1, hspace=0.04, wspace=0.05,
                  width_ratios=[1] * ncol + [0.05])
    axes = [[fig.add_subplot(gs[r, c]) for c in range(ncol)] for r in range(nrow)]
    cax = fig.add_subplot(gs[:, ncol])

    im = None
    for r in range(nrow):
        for c in range(ncol):
            p = panels[r][c]
            if p is None:
                axes[r][c].set_visible(False)
                continue
            im = draw_panel(axes[r][c], *p, vmin=vmin, vmax=vmax, base_snr=TABLE_SNR)

    for r in range(nrow):
        for c in range(ncol):
            ax = axes[r][c]
            ax.set_ylim(*ylim)
            ax.tick_params(direction='in', top=True, right=True)
            ax.set_xticklabels([]) if r != nrow - 1 else ax.set_xlabel(r'$T_{\rm eff}$ [K]')
            ax.set_yticklabels([]) if c != 0 else ax.set_ylabel(r'$\log g$')
            if r == 0:
                ax.set_title(feh_labels[c], fontsize=10)
            if c == ncol - 1:
                ax.text(0.97, 0.06, f'[Ce/Fe] = {labels[r]:+.2f}', transform=ax.transAxes,
                        ha='right', va='bottom', fontsize=9,
                        bbox=dict(fc='white', ec='0.6', alpha=0.85, pad=2))

    cb = fig.colorbar(im, cax=cax)
    cb.set_label(r'max $|\partial f/\partial$[Ce/Fe]$|$', fontsize=10)
    axes[0][0].legend(handles=[Line2D([], [], color=SNR_STYLE[s][0],
                                      ls=SNR_STYLE[s][1], lw=1.6, label=f'S/N = {s}')
                               for s in snr_list if s in SNR_STYLE],
                      fontsize=8, loc='upper left', framealpha=0.85)
    fig.savefig(fname, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print('wrote', fname)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--snr', type=float, nargs='+', default=DEFAULT_SNR)
    ap.add_argument('--outdir', default='figures')
    ap.add_argument('--table')
    args = ap.parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    df = CeDetectability(csv=args.table).table
    snr_list = [int(s) for s in args.snr]

    # Fig. 2 analogue: two [Ce/Fe], two metallicity bins
    figure(df, [0.0, 0.5], [(-0.5, 0.0), (0.0, 0.5)],
           [r'$-0.5 \leq$ [Fe/H] $< 0.0$', r'$0.0 \leq$ [Fe/H] $\leq 0.5$'],
           snr_list, (3.8, 0), os.path.join(args.outdir, 'fig2_detectability.png'), (7, 5))

    # Fig. C.1 analogue: full [Ce/Fe] range, three metallicity bins
    figure(df, [-0.45, 0.0, 0.5, 0.95],
           [(-1.5, -0.5), (-0.5, 0.0), (0.0, 0.5)],
           [r'$-1.5 \leq$ [Fe/H] $< -0.5$', r'$-0.5 \leq$ [Fe/H] $< 0.0$',
            r'$0.0 \leq$ [Fe/H] $\leq 0.5$'],
           snr_list, (5.0, 0), os.path.join(args.outdir, 'figC1_detectability.png'), (10, 8))


if __name__ == '__main__':
    main()
