# Ce detectability in *Gaia* RVS spectra

Is the cerium abundance of a star actually measurable from its *Gaia* RVS
spectrum?

This repository answers that for any `(Teff, log g, [Fe/H], [Ce/Fe])` and
signal-to-noise ratio, using the criterion of

> Vitali et al., *Detectability limits of Ce abundances in Gaia RVS*

The Ce II triplet near 851.1 nm is weak and unresolved at the RVS resolving
power (R ≈ 11 500), so a reported [Ce/Fe] is not automatically a measurement.
The tool tells you whether the feature carries usable information for a given
star, and — importantly — reports the two ingredients of that judgement
separately.

![Ce detectability in the Kiel diagram](figures/fig2_detectability.png)

## Install

```bash
git clone https://github.com/<user>/ce-rvs-detectability.git
cd ce-rvs-detectability
pip install -r requirements.txt
```

Needs `numpy`, `pandas`, `scipy` (and `matplotlib`, `astropy` for the figures
and the grid tools). No installation step beyond that — the scripts run in place.

## Use

A single star:

```bash
python3 ce_detectability.py --star 4500 1.8 -0.2 0.4 --snr 50
```

```
requested : Teff=4500  logg=1.80  [Fe/H]=-0.20  [Ce/Fe]=+0.40  S/N=50
matched   : Teff=4500  logg=1.80  [Fe/H]=-0.20  [Ce/Fe]=+0.45  (distance 0.50 cells)
            [alpha/Fe] = +0.1 (MARCS standard, fixed by the grid)
CRB       : sigma = 0.072 dex  -> PASS (< 0.1 dex)
EW        : 1.2 sigma            -> FAIL (>= 3 sigma)
==> NOT detectable
```

Arguments to `--star` are in the order `Teff log g [Fe/H] [Ce/Fe]`.

A catalogue — needs columns `teff, logg, feh, cefe`, optionally a per-star `snr`:

```bash
python3 ce_detectability.py --input example_stars.csv --output flagged.csv --snr 50
```

| name | snr | grid_distance | sigma_cefe | ew_snr | crb_pass | ew_pass | detectable |
|---|---|---|---|---|---|---|---|
| cool_giant | 50 | 0.50 | 0.022 | 4.96 | True | True | **True** |
| mild_giant | 50 | 0.50 | 0.072 | 1.22 | True | False | False |
| metal_poor_giant | 50 | 0.71 | 0.050 | 1.39 | True | False | False |
| solar_dwarf | 50 | 0.50 | 1.107 | 0.01 | False | False | False |
| high_snr_giant | 200 | 0.50 | 0.018 | 4.89 | True | True | **True** |

From Python:

```python
from ce_detectability import CeDetectability
ce  = CeDetectability()
out = ce.classify(df.teff, df.logg, df.feh, df.cefe, snr=df.snr)
out[['crb_pass', 'ew_pass', 'detectable']]
```

## The criterion

Two conditions, evaluated at the requested S/N:

| flag | test | question it answers |
|---|---|---|
| `crb_pass` | σ<sub>CRB</sub> < 0.1 dex | does the spectrum carry enough **information** to pin [Ce/Fe]? |
| `ew_pass` | EW ≥ 3 σ<sub>W</sub> | is the **line visible** above the noise? |
| `detectable` | both | the paper's criterion |

with

```
sigma_CRB = 1 / (SNR * sqrt(sum_i g_i^2)),   g_i = df_i/d[Ce/Fe]
sigma_W   = (1.5 / SNR) * sqrt(FWHM * dx),   FWHM = 0.0936 nm, dx = 0.03 nm
```

σ<sub>CRB</sub> is the Cramér–Rao bound (Cramér 1945; Rao 1945), the best precision
any unbiased estimator could reach. σ<sub>W</sub> is the Cayrel (1988) equivalent-width
uncertainty. The equivalent width is integrated against a local linear continuum
taken from the edges of the 851.05–851.25 nm window, so a sloped pseudo-continuum
or a blended neighbour is not counted as Ce absorption.

**Why the two flags are reported separately.** The Cramér–Rao bound sums signal
in quadrature across pixels, so a line can be "statistically detectable" while
its core sits well below the per-pixel noise — a regime highly sensitive to
continuum placement and blending. The Cayrel test is the reality check. Stars
with `crb_pass=True, ew_pass=False` are exactly that marginal population
(`mild_giant` above), and deserve separate treatment rather than being merged
into the non-detections.

Both limits scale as 1/SNR, so any S/N is obtained by rescaling the tabulated
S/N = 50 reference.

## Output columns

| column | meaning |
|---|---|
| `crb_pass`, `ew_pass`, `detectable` | the three flags above |
| `sigma_cefe` | Cramér–Rao bound on [Ce/Fe] at your S/N, dex |
| `ew_ce`, `ew_limit` | measured EW and the 3σ threshold, nm |
| `ew_snr` | detection significance, in σ (the test is ≥ 3) |
| `grad_rss` | √(Σ g<sub>i</sub>²), the Ce sensitivity of the spectrum |
| `*_grid` | the grid node the star was matched to |
| `grid_distance` | distance to that node, in grid cells |
| `alpha_assumed` | the [α/Fe] implied by [Fe/H] (see below) |
| `in_grid_range` | False if the star lies outside the grid's coverage |
| `delta_cefe` | [Ce/Fe] step of the finite difference (normally 0.1) |

## Two things to watch

**`grid_distance`** is measured in grid cells (Teff/25 K, log g/0.2 dex,
[Fe/H]/0.1 dex, [Ce/Fe]/0.1 dex). 0 means the star sits on a node, up to ~0.5
that it sits between nodes. The criterion is a **step function** at the boundary,
so a star with `grid_distance` ≳ 0.5 landing near the contour is genuinely
ambiguous — flag it rather than quoting a hard yes/no. `in_grid_range` catches
stars outside the grid entirely; those are snapped to the nearest edge node and
should not be trusted.

**[α/Fe] is not a free parameter.** The synthesis fixed it to the standard MARCS
metallicity-dependent relation (+0.4 for [Fe/H] ≤ −1, falling linearly to 0.0 at
[Fe/H] ≥ 0), so it is determined by [Fe/H] and cannot be supplied. The tool
reports the assumed value in `alpha_assumed`. A star whose true [α/Fe] departs
from that relation is not covered: α-enhancement deepens the blends that set the
local pseudo-continuum and suppresses the measured Ce EW, so the boundary shifts
slightly. The effect is subdominant to continuum-placement error over a realistic
thin/thick-disc spread, but it is not zero.

## Grid coverage

| axis | range | step |
|---|---|---|
| Teff | 3700 – 6025 K | 25 K |
| log g | 0.0 – 4.8 | 0.2 |
| [Fe/H] | −2.0 – +0.5 | 0.1 |
| [Ce/Fe] | −0.5 – +2.0 | 0.1 |

61 100 (Teff, log g, [Fe/H]) nodes, 1 527 482 table entries.

Gradients are finite differences between consecutive [Ce/Fe] values, so the
tabulated [Ce/Fe] are **midpoints** (±0.05, ±0.15, …). A request for
[Ce/Fe] = 0.0 matches the −0.05 or +0.05 midpoint; the figures are labelled with
the midpoint actually used. In 16 of the 1 527 482 entries a synthesis is
missing, so the difference spans 0.2 dex instead of 0.1 — see `delta_cefe`.

## How the table was produced

`data/ce_detectability_table.csv.gz` is derived from a grid of **1 588 582
synthetic spectra** computed with Turbospectrum (Alvarez & Plez 1998; Plez 2012),
MARCS model atmospheres (Gustafsson et al. 2008), the Gaia-ESO v6 line list
(Heiter et al. 2021) with hyperfine structure and isotopic splitting, and the
solar scale of Grevesse et al. (2007), convolved to R = 11 500 and sampled at
0.001 nm. Microturbulence was 1 km/s for dwarfs and 2 km/s for giants; no
rotation or macroturbulence. Full details in Appendix A of the paper.

That grid is ~37 GB and is not distributed here. The two scripts that produced
the table from it are included, so the method is reproducible and auditable:

```bash
# 1. reduce the raw spectra to the Ce window (37 GB -> ~1.3 GB)
python3 condense_grid.py /path/to/synthetic_spectra -o ce_grid_window.fits

# 2. compute gradients, CRB and EW -> the shipped table
python3 compute_gradients.py ce_grid_window.fits -o ce_detectability_table.csv
```

`compute_gradients.py` also accepts the raw `synth_*.txt` directory directly.
`condense_grid.py --sample-nodes N` extracts a small example grid, keeping whole
(Teff, log g, [Fe/H]) nodes so every retained node keeps its full [Ce/Fe]
sequence and the gradient stays defined.

## Figures

```bash
python3 plot_detectability.py
```

writes `figures/fig2_detectability.png` and `figures/figC1_detectability.png`,
the detectability maps across the Kiel diagram for several metallicity and
[Ce/Fe] bins, with contours at S/N = 25, 50 and 200.

## Caveats

The boundary assumes a perfectly known continuum and exact atmospheric
parameters. It neglects continuum-placement bias (which dominates above ~1%
error), uncertainties in Teff and log g, and is derived from synthetic spectra
rather than from the output of any particular pipeline. Treat it as a guide to
where Ce abundances *can* be trusted, not as a per-star error bar.

## Citation

If you use this, please cite the paper:

```bibtex
@article{vitali_ce_rvs,
  author  = {Vitali, S. and others},
  title   = {Detectability limits of Ce abundances in Gaia RVS},
  journal = {Astronomy \& Astrophysics},
  year    = {2026}
}
```

## Licence

MIT — see `LICENSE`.
