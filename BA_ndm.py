#!/usr/bin/env python3
"""
BA_ndm.py - Brain age from normative models (normative deviation mapping, NDM).

Prototype of a likelihood-based brain age estimator as an alternative to the
GPR-based BrainAGE (BA_gpr_ui.m).  Every feature gets its own normative model

    y_j ~ N(mu_j(age, sex, site, covariates), sigma_j(age, sex, covariates)^2)

with natural cubic splines of age for mu and log(sigma), fitted by maximum
likelihood (the location-scale RS algorithm of ComBatLS, ported from
combat_family.py in the ComCat repository).  Brain age is the age at which the
subject's data are most likely,

    age_hat = argmax_a  log N(z(a); 0, R) - sum_j log sigma_j(a),

where z(a) are the z-scores at candidate age a and R is the correlation
between features, estimated from the training z-scores.

Two kinds of models are fitted to the same data:

  normative model  a location-scale model for every voxel/vertex, which gives
                   the z-maps of test subjects (--zmaps; --normative-only fits
                   and applies only this model).  It never uses PCA.
  NDM brain age    the likelihood above.  By default (--pca 100) its features
                   are the leading principal component scores of the training
                   data (of the whole brain and of each region), so that R can
                   be estimated in full.  --pca, --rank, --psi-min and the age
                   grid only affect the brain age and its deviation scores.

With --pca 0 every voxel/vertex is a feature of the brain age and R is a
low-rank plus diagonal model; the voxel/vertex-wise normative model is then
shared with the z-maps and fitted only once.  The residual correlation of
smoothed voxels is spread over hundreds of dimensions, so this model is
overconfident (standard errors several times too small) and less accurate.
Treating features as independent (--pca 0 --rank 0) overcounts the evidence of
correlated features and is worse still.

The estimate is unbiased conditional on age without a trend correction, comes
with a per-subject standard error (Laplace approximation) and gives regional
brain ages when the likelihood is restricted to the features of a region (the
lobe atlas used by BA_gpr_ui.m with D.parcellation = 1).  Several models
(tissue, resolution, smoothing, surface measure) are combined by a weighted
average with weights that sum to one, which keeps the ensemble unbiased
conditional on age.

Besides brain age, every subject gets a non-aging deviation: the Mahalanobis
distance of its z-scores at the estimated brain age, i.e. the atypicality that
an older or younger brain does not explain (normal score of a chi-square
distance; NDM.Deviation, also per lobe), and the total deviation at
chronological age (NDM.Deviation_age).

Covariates such as image quality measures (IQMs) can enter the mean and the
log SD of all normative models (--train-cov, --test-cov, --cov-mean, --cov-sd,
--cov-df).  One table with a row per subject serves all models.  The models
condition on every subject's own covariates, so z-maps and brain ages are
adjusted for them without a reference value.  Covariates that are constant or
collinear with the rest of the design (e.g. constant within each site) are not
used.  In tests, an IQM covariate removed its effect from the z-maps but left
the brain age about as accurate as without it, and a covariate that correlates
strongly with age made the brain age less accurate.

Non-Gaussian data can be warped per voxel/vertex (sinh-arcsinh, fitted jointly
with a voxel-wise normative model, as in warped Bayesian linear regression):
--warp-zmaps for the z-maps, where it calibrates the tails (recommended), and
--warp for the brain age models, where it made no difference in tests.

Fitted models can be saved with --save-model (<name>.mat: MATLAB struct
NDMmodel, <name>.npz: NumPy archive) and applied to new data with --model
instead of --train, without the training data.  The model file contains a JSON
description (settings, training sample, covariates and the feature space of
every model), which is also written to <name>.json and used to check that new
data are compatible.

For comparison, a Python replica of the GPR BrainAGE (BA_gpr.m, linear kernel,
PCA, linear trend correction with trend_method = 1) is run on the same folds.
It needs the training data and is skipped with --model.

Input files are the mat-files written by BA_data2mat.m (Y, age, male and, for
surface data, ind; MATLAB v5 or v7.3).  Files joined with '+' are concatenated
and their position serves as site.  Subjects with non-finite age or age <= 0
(or with missing covariates) are excluded from training and evaluation.

Examples
--------
10-fold cross-validation on NKIe with 4 models and lobe-wise brain age:

    python BA_ndm.py --train s4rp1_4mm_NKIe1239_CAT12.9.mat \\
        s4rp1_8mm_NKIe1239_CAT12.9.mat s4rp2_4mm_NKIe1239_CAT12.9.mat \\
        s4rp2_8mm_NKIe1239_CAT12.9.mat --kfold 10 --parcellation --out NKIe_ndm

Train on one sample and predict another (models are matched by position).
The control subjects given by --adjust (1-based, as D.ind_adjust) correct the
estimates for the test site and estimate the ensemble weights.  By default
(--correction offset) their median BrainAGE is subtracted from all subjects,
without a trend correction:

    python BA_ndm.py --train s4rp1_8mm_A_CAT12.9.mat s4rp2_8mm_A_CAT12.9.mat \\
        --test s4rp1_8mm_B_CAT12.9.mat s4rp2_8mm_B_CAT12.9.mat --adjust 1:108

With --correction agefree, the normative models are adapted to the test site
without the controls' ages (iteratively at their estimated brain ages); the ages
are then only used for the final offset.

Fit the voxel-wise normative model with two IQMs in the mean, save it, and
compute z-maps of new data later without the training data:

    python BA_ndm.py --normative-only --train s4rp1_4mm_A_CAT12.9.mat \\
        --train-cov IQM_A.csv --cov-mean SIQR,res_ECR --save-model A_norm.mat
    python BA_ndm.py --normative-only --model A_norm.mat --parcellation \\
        --test s4rp1_4mm_B_CAT12.9.mat --test-cov IQM_B.csv --adjust 1:108 --out B

Brain age models are saved and applied in the same way without
--normative-only (with --zmaps, the voxel-wise models are included).

Results are saved as <out>.mat (readable by MATLAB) and <out>.csv, and z-maps
as <out>_zmaps_<model>.mat.  With --normative-only, only the z-maps are saved,
plus the lobe-wise mean z in <out>.csv with --parcellation.

Requirements: numpy, scipy, h5py (v7.3 files), nibabel (surface parcellation)
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import os
import re
import time
import warnings
from dataclasses import asdict, dataclass, replace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# regions of the lobe atlas that are not represented on the surface
_SURF_EXCLUDE = (5, 15)

# format of saved models (save_models, load_models)
MODEL_FORMAT, MODEL_VERSION = 'BA_ndm model', 1

_NOTES = set()


def _note(msg):
    """Print a message only once."""
    if msg not in _NOTES:
        _NOTES.add(msg)
        print(msg, flush=True)


# ---------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------

@dataclass
class Data:
    """Data of one model (segmentation/resolution/smoothing/surface measure)."""
    Y: np.ndarray                  # (n_subjects, n_features) float32
    age: np.ndarray                # (n_subjects,)
    male: np.ndarray               # (n_subjects,)
    site: np.ndarray               # (n_subjects,) index of '+'-joined file
    name: str
    ind: np.ndarray | None = None  # 1-based vertex indices of surface data
    is_surf: bool = False
    res: str | None = None         # resampling of volume data, e.g. '8'
    has_male: bool = True          # False if the mat-file contains no sex
    dim: np.ndarray | None = None  # volume dimensions stored in the mat-file
    C: np.ndarray | None = None    # (n_subjects, n_covariates) covariates
    cov_names: tuple = ()          # names of the columns of C

    @property
    def n(self):
        return self.Y.shape[0]

    def subset(self, idx):
        return replace(self, Y=self.Y[idx], age=self.age[idx], male=self.male[idx],
                       site=self.site[idx], C=None if self.C is None else self.C[idx])


@dataclass(frozen=True)
class Covariates:
    """Covariates of the normative models.

    names : columns of Data.C, in this order
    mean  : covariates in mu; sd: covariates in log sigma (subsets of names)
    df    : degrees of freedom of their natural splines (1 = linear)
    """
    names: tuple = ()
    mean: tuple = ()
    sd: tuple = ()
    df: int = 1

    @classmethod
    def from_state(cls, s):
        return cls(tuple(s['names']), tuple(s['mean']), tuple(s['sd']), int(s['df']))


def _read_mat(path):
    """Read Y, age, male, ind and dim of a BA_data2mat mat-file (v5 or v7.3)."""
    keys = ('Y', 'age', 'male', 'ind', 'dim')
    try:
        from scipy.io import loadmat
        m = loadmat(path, variable_names=list(keys))
        out = {k: (m[k] if k in m and m[k].size else None) for k in keys}
    except NotImplementedError:    # v7.3 is HDF5 and stores arrays transposed
        import h5py
        out = {}
        with h5py.File(path, 'r') as h:
            for k in keys:
                if k not in h or h[k].attrs.get('MATLAB_empty', 0):
                    out[k] = None
                else:
                    out[k] = h[k][()].T
    if out['Y'] is None or out['age'] is None:
        raise ValueError(f"{path} must contain Y and age.")
    return out


def load_data(spec: str) -> Data:
    """Load a mat-file, or several files joined with '+' (one site per file)."""
    parts = spec.split('+')
    Ys, ages, males, sites, ind, dim = [], [], [], [], None, None
    has_male = True
    for s, path in enumerate(parts):
        m = _read_mat(path)
        age = np.asarray(m['age'], dtype=np.float64).ravel()
        Y = m['Y']
        if Y.shape[0] != age.size and Y.shape[1] == age.size:
            Y = Y.T
        if Y.shape[0] != age.size:
            raise ValueError(f"{path}: Y {Y.shape} does not match {age.size} ages.")
        has_male &= m['male'] is not None
        male = (np.zeros_like(age) if m['male'] is None
                else np.asarray(m['male'], dtype=np.float64).ravel())
        if m['ind'] is not None:
            if ind is not None and not np.array_equal(ind, m['ind'].ravel()):
                raise ValueError(f"{path}: surface index differs between '+'-joined files.")
            ind = m['ind'].ravel()
        if m['dim'] is not None:
            dm = np.asarray(m['dim'], dtype=np.float64).ravel()
            if dim is not None and not np.array_equal(dim, dm):
                raise ValueError(f"{path}: volume dimensions differ between '+'-joined files.")
            dim = dm
        Ys.append(np.asarray(Y, dtype=np.float32))
        ages.append(age)
        males.append(male)
        sites.append(np.full(age.size, s))
    base = os.path.basename(parts[0])
    res = re.search(r'_(\d+)mm_', base)
    return Data(np.vstack(Ys), np.concatenate(ages), np.concatenate(males),
                np.concatenate(sites), '+'.join(os.path.basename(p) for p in parts),
                ind, 'mesh' in base, res.group(1) if res else None, has_male, dim)


def read_table(spec):
    """Covariate table(s) with one row per subject: (names, values, numeric).

    Text files separated by commas, semicolons, tabs or white space.  A first row
    with a non-numeric cell is the header; without header the columns are named
    by their 1-based number.  Empty cells and NaN/NA are missing values; other
    non-numeric cells (e.g. subject IDs) mark their column as non-numeric.
    Files joined with '+' (e.g. one per '+'-joined sample) are concatenated and
    must have the same columns.
    """
    def num(c):
        try:
            return float(c)
        except ValueError:
            return None

    missing = ('', 'nan', 'na')
    names, blocks, numeric = None, [], None
    for path in spec.split('+'):
        with open(path) as f:
            lines = [ln.rstrip('\r\n') for ln in f]
        lines = [ln for ln in lines if ln.strip() and not ln.lstrip().startswith('#')]
        if not lines:
            raise ValueError(f"{path}: empty table.")
        sep = next((s for s in (',', ';', '\t') if s in lines[0]), None)
        rows = [[c.strip().strip('"\'') for c in (ln.split(sep) if sep else ln.split())]
                for ln in lines]
        header = any(num(c) is None and c.lower() not in missing for c in rows[0])
        nm = rows[0] if header else [str(j + 1) for j in range(len(rows[0]))]
        body = rows[1:] if header else rows
        vals = np.full((len(body), len(nm)), np.nan)
        isnum = np.ones(len(nm), dtype=bool)
        for i, row in enumerate(body):
            if len(row) != len(nm):
                raise ValueError(f"{path}: line {i + 1 + header} has {len(row)} columns "
                                 f"instead of {len(nm)}.")
            for j, c in enumerate(row):
                v = num(c)
                if v is not None:
                    vals[i, j] = v
                elif c.lower() not in missing:
                    isnum[j] = False
        if names is not None and nm != names:
            raise ValueError(f"{path}: the columns differ from those of the first table.")
        names = nm
        numeric = isnum if numeric is None else numeric & isnum
        blocks.append(vals)
    return names, np.vstack(blocks), numeric


def _parse_columns(spec, names, numeric):
    """Column names of a comma-separated list of names or 1-based numbers."""
    out = []
    for tok in (t.strip() for t in spec.split(',')):
        if not tok:
            continue
        if tok in names:
            nm = tok
        elif tok.isdigit() and 1 <= int(tok) <= len(names):
            nm = names[int(tok) - 1]
        else:
            raise ValueError(f"Covariate {tok!r} is not a column of the table "
                             f"({', '.join(names)}).")
        if not numeric[names.index(nm)]:
            raise ValueError(f"Covariate {nm!r} has non-numeric values.")
        if nm not in out:
            out.append(nm)
    return out


def covariates_from_table(names, values, numeric, mean=None, sd=None, df=1):
    """Covariates of a table (read_table).  mean and sd are comma-separated column
    names or numbers as for --cov-mean and --cov-sd; by default all numeric
    columns with values enter the mean and none the log SD ('none': none)."""
    def cols(spec, default):
        if spec is None:
            return default
        return [] if spec.strip().lower() == 'none' else _parse_columns(spec, names, numeric)

    m = cols(mean, [nm for j, nm in enumerate(names)
                    if numeric[j] and np.any(np.isfinite(values[:, j]))])
    s = cols(sd, [])
    if not m and not s:
        raise ValueError("No covariates selected.")
    if df < 1:
        raise ValueError("The spline df of the covariates must be at least 1.")
    return Covariates(tuple(dict.fromkeys(m + s)), tuple(m), tuple(s), int(df))


def covariate_matrix(names, values, cov, n, label, required=None):
    """The columns cov.names of a covariate table (n, n_covariates).  Columns that
    are not required (default: all are) may be missing and are then NaN."""
    required = cov.names if required is None else required
    missing = [nm for nm in required if nm not in names]
    if missing:
        raise ValueError(f"{label}: covariate(s) {', '.join(missing)} not in the table "
                         f"({', '.join(names)}).")
    if len(values) != n:
        raise ValueError(f"{label}: {len(values)} rows for {n} subjects.")
    return np.column_stack([values[:, names.index(nm)] if nm in names else np.full(n, np.nan)
                            for nm in cov.names])


def _region_names(atlas_dir):
    names = {}
    path = os.path.join(atlas_dir, 'Brain_Lobes.csv')
    if os.path.isfile(path):
        with open(path) as f:
            for line in f.read().splitlines()[1:]:
                rid, _, nm = line.partition(';')
                if rid.strip():
                    names[int(rid)] = nm.strip()
    return names


def lobe_atlas(d: Data, atlas_dir=HERE):
    """Lobe label of every feature, the regions used, and their names.

    Same atlas files and conventions as BA_gpr_ui.m (D.parcellation = 1).
    """
    names = _region_names(atlas_dir)
    if d.is_surf:
        import nibabel.freesurfer as fs
        if d.ind is None:
            raise ValueError(f"No index 'ind' of surface values found in {d.name}. "
                             "Please re-create the data using BA_data2mat.")
        # merged 32k meshes are ordered lh first, then rh
        lab = np.concatenate([
            fs.read_annot(os.path.join(atlas_dir, f'{h}.Brain_Lobes.annot'),
                          orig_ids=True)[0] for h in ('lh', 'rh')])
        atlas = lab[d.ind.astype(int) - 1]
        exclude = _SURF_EXCLUDE
    else:
        from scipy.io import loadmat
        if d.res is None:
            raise ValueError(f"Cannot infer resolution from {d.name}.")
        atlas = loadmat(os.path.join(atlas_dir, f'Brain_Lobes_{d.res}mm.mat'))['atlas'].ravel()
        exclude = ()
    if atlas.size != d.Y.shape[1]:
        raise ValueError(f"Atlas has {atlas.size} entries but {d.name} has "
                         f"{d.Y.shape[1]} features.")
    regions = [int(r) for r in np.unique(atlas[atlas > 0]) if r not in exclude]
    if not regions:
        raise ValueError("No regions found in lobe atlas.")
    return atlas.astype(int), regions, [names.get(r, str(r)) for r in regions]


# ---------------------------------------------------------------------------
# Normative model
# ---------------------------------------------------------------------------

class NaturalSpline:
    """Natural cubic spline basis (of age or a covariate) without intercept (The Elements of 
    Statistical Learning eq. 5.4-5.5).

    Knots at quantiles of the training ages; the basis is linear beyond the
    boundary knots, so the model extrapolates linearly.  df=1 is linear.
    """

    def __init__(self, x, df):
        x = np.asarray(x, dtype=np.float64)
        knots = np.unique(np.quantile(x, np.linspace(0, 1, max(df, 1) + 1)))
        self.lo, self.hi = knots[0], knots[-1]
        self.knots = (knots - self.lo) / (self.hi - self.lo)
        self.df = len(self.knots) - 1

    def __call__(self, x):
        t = (np.asarray(x, dtype=np.float64) - self.lo) / (self.hi - self.lo)
        k = self.knots
        if len(k) <= 2:
            return t[:, None]

        def d(i):
            return (np.maximum(t - k[i], 0) ** 3
                    - np.maximum(t - k[-1], 0) ** 3) / (k[-1] - k[i])
        d_last = d(len(k) - 2)
        return np.column_stack([t] + [d(i) - d_last for i in range(len(k) - 2)])

    def get_state(self):
        return dict(lo=float(self.lo), hi=float(self.hi),
                    knots=np.asarray(self.knots, dtype=np.float64))

    @classmethod
    def from_state(cls, s):
        self = cls.__new__(cls)
        self.lo, self.hi = float(s['lo']), float(s['hi'])
        self.knots = np.asarray(s['knots'], dtype=np.float64).ravel()
        self.df = len(self.knots) - 1
        return self


def _deviance(r, eta):
    """-2 log-likelihood of N(mu, exp(eta)^2) given residuals r = y - mu."""
    return np.sum(np.log(2 * np.pi) + 2 * eta + r ** 2 * np.exp(-2 * eta), axis=0)


def fit_location_scale(Y, X, W, max_iter=2000, tol=1e-6, max_elements=4_000_000, init=None):
    """ML fit of y_j ~ N(X beta_j, exp(W theta_j)^2) for every column j of Y.

    Port of _fit_location_scale() in combat_family.py (ComCat) for data of
    shape (n_subjects, n_features).  RS algorithm as in gamlss (family NO,
    identity link for mu, log link for sigma): alternate a weighted
    least-squares update of beta with a Fisher-scoring update of theta,
    halving the theta step if the deviance increases.  W[:, 0] must be the
    intercept.  Features are processed in chunks to bound memory.  init =
    (beta, theta) starts from a previous solution instead of OLS.

    Returns beta (kx, p), theta (kw, p), converged (p,)
    """
    n, p = Y.shape
    kx, kw = X.shape[1], W.shape[1]
    XX = (X[:, :, None] * X[:, None, :]).reshape(n, kx * kx)
    X_pinv = np.linalg.pinv(X)
    W_pinv = np.linalg.pinv(W)

    beta = np.empty((kx, p))
    theta = np.empty((kw, p))
    converged = np.zeros(p, dtype=bool)

    step = max(1, max_elements // n)
    for start in range(0, p, step):
        y_all = np.asarray(Y[:, start:start + step], dtype=np.float64)

        if init is not None:
            b_all = np.array(init[0][:, start:start + step], dtype=np.float64)
            t_all = np.array(init[1][:, start:start + step], dtype=np.float64)
        else:                                      # start: OLS mean, constant sigma
            b_all = X_pinv @ y_all
            r = y_all - X @ b_all
            t_all = np.zeros((kw, y_all.shape[1]))
            t_all[0] = 0.5 * np.log(np.mean(r ** 2, axis=0))

        active = np.arange(y_all.shape[1])
        for _ in range(max_iter):
            y, b, t = y_all[:, active], b_all[:, active], t_all[:, active]
            eta = W @ t

            # mu step: weighted least squares with weights 1/sigma^2
            w = np.exp(-2 * eta)
            A = (XX.T @ w).T.reshape(-1, kx, kx)
            rhs = (X.T @ (w * y)).T[..., None]
            b_new = np.linalg.solve(A, rhs)[..., 0].T
            d_mu = np.max(np.abs(X @ (b_new - b)) * np.exp(-eta), axis=0)
            r = y - X @ b_new

            # sigma step: Fisher scoring on log(sigma)
            dev_old = _deviance(r, eta)
            d_t = W_pinv @ (0.5 * (r ** 2 * w - 1))
            t_new = t + d_t
            eta_new = W @ t_new
            dev_new = _deviance(r, eta_new)
            for _ in range(20):                    # step halving (gamlss autostep)
                worse = dev_new > dev_old
                if not np.any(worse):
                    break
                d_t[:, worse] *= 0.5
                t_new[:, worse] = t[:, worse] + d_t[:, worse]
                eta_new[:, worse] = W @ t_new[:, worse]
                dev_new[worse] = _deviance(r[:, worse], eta_new[:, worse])
            d_theta = np.max(np.abs(t_new - t), axis=0)

            b_all[:, active] = b_new
            t_all[:, active] = t_new
            done = (d_mu < tol) & (d_theta < tol)
            converged[start + active[done]] = True
            active = active[~done]
            if active.size == 0:
                break

        beta[:, start:start + step] = b_all
        theta[:, start:start + step] = t_all

    return beta, theta, converged


def _rank(M, rtol=1e-8):
    """Numerical rank of M relative to its largest singular value."""
    s = np.linalg.svd(M, compute_uv=False)
    return int(np.sum(s > rtol * s[0])) if s.size and s[0] > 0 else 0


class NormativeModel:
    """Location-scale normative model for every feature.

    mu        = 1 + ns(age, df_mu) + male + site + covariates   (site: fixed effects)
    log sigma = 1 + ns(age, df_sigma) + male + covariates

    Covariates (cov: columns of C) enter as natural splines with cov.df degrees
    of freedom (1 = linear), centred at their training median, so that the model
    without covariates (C=None) refers to median covariates.  Covariates that are
    constant or collinear with the rest of the design are not used.  Predictions
    condition on every subject's own covariates (missing values: median) and use
    the size-weighted mean of the training sites (as ComBat's stand_mean),
    optionally adapted to a new site by adapt().
    """

    def __init__(self, df_mu=5, df_sigma=3, cov=None, max_iter=2000, tol=1e-6):
        self.df_mu, self.df_sigma = df_mu, df_sigma
        self.cov = cov or Covariates()
        self.max_iter, self.tol = max_iter, tol

    def _cov_columns(self, C, n, names):
        """Spline bases (n, k) of the covariates names, centred at their median."""
        cols = []
        for nm in names:
            b, med = self.cov_basis[nm], self.cov_median[nm]
            if C is None:
                cols.append(np.zeros((n, b.df)))
                continue
            x = np.asarray(C[:, self.cov.names.index(nm)], dtype=np.float64)
            x = np.where(np.isfinite(x), x, med)
            cols.append(b(x) - b(np.array([med])))
        return np.hstack(cols) if cols else np.zeros((n, 0))

    def _design(self, age, male, site=None, C=None):
        n = len(age)
        one = np.ones((n, 1))
        X = [one, self.basis_mu(age)]
        W = [one, self.basis_sigma(age)]
        if self.use_male:
            X.append(np.asarray(male, dtype=np.float64)[:, None])
            W.append(np.asarray(male, dtype=np.float64)[:, None])
        if len(self.site_levels) > 1:
            if site is None:     # reference: size-weighted mean of the sites
                X.append(np.repeat(self.site_props[None, 1:], n, axis=0))
            else:
                codes = np.searchsorted(self.site_levels, site)
                X.append((codes[:, None] == np.arange(1, len(self.site_levels))).astype(float))
        X.append(self._cov_columns(C, n, self.cov_mean))      # covariates are the last
        W.append(self._cov_columns(C, n, self.cov_sd))        # columns (cov_effects)
        return np.hstack(X), np.hstack(W)

    def _setup(self, age, male, site, C=None):
        """Bases of age and covariates, sites and sex coding of the training data;
        drops constant and collinear covariates; returns site."""
        n = len(age)
        site = np.zeros(n, dtype=int) if site is None else np.asarray(site)
        self.site_levels, counts = np.unique(site, return_counts=True)
        self.site_props = counts / n
        self.use_male = np.unique(male).size > 1
        self.basis_mu = NaturalSpline(age, self.df_mu)
        self.basis_sigma = NaturalSpline(age, self.df_sigma)
        self.cov_basis, self.cov_median = {}, {}
        self.cov_mean, self.cov_sd = [], []
        if not self.cov.names:
            return site
        if C is None or C.shape[1] != len(self.cov.names):
            raise ValueError(f"The normative model needs the covariates {', '.join(self.cov.names)}.")
        for j, nm in enumerate(self.cov.names):
            x = np.asarray(C[:, j], dtype=np.float64)
            if not np.all(np.isfinite(x)):
                raise ValueError(f"Covariate {nm} has missing values in the training data.")
            if x.max() > x.min():
                self.cov_basis[nm] = NaturalSpline(x, self.cov.df)
                self.cov_median[nm] = float(np.median(x))
            else:
                _note(f"Covariate {nm} is constant in the training data and not used.")
        X, W = self._design(age, male, site)
        for names, M, used, part in ((self.cov.mean, X, self.cov_mean, 'mean'),
                                     (self.cov.sd, W, self.cov_sd, 'log SD')):
            r = _rank(M)
            for nm in names:
                if nm not in self.cov_basis:
                    continue
                cols = self._cov_columns(C, n, [nm])
                Mn = np.hstack([M, cols])
                rn = _rank(Mn)
                if rn == r + cols.shape[1]:
                    M, r = Mn, rn
                    used.append(nm)
                else:
                    _note(f"Covariate {nm} is collinear with the design of the {part} "
                          "and not used there.")
        return site

    def fit(self, Y, age, male, site=None, verbose=False, C=None):
        site = self._setup(age, male, site, C)

        finite = np.all(np.isfinite(Y), axis=0)
        sd = np.zeros(Y.shape[1])
        sd[finite] = np.std(Y[:, finite], axis=0)
        self.valid = np.flatnonzero(finite & (sd > 0))
        self.n_features = Y.shape[1]

        X, W = self._design(age, male, site, C)
        t0 = time.time()
        self.beta, self.theta, conv = fit_location_scale(
            Y[:, self.valid], X, W, self.max_iter, self.tol)
        if verbose:
            print(f"    normative model: {self.valid.size} features, "
                  f"{int(np.sum(~conv))} not converged ({time.time() - t0:.1f}s)")
        self.offset = np.zeros(self.valid.size)
        self.scale = np.ones(self.valid.size)
        return self

    def mu_sigma(self, age, male, site=None, cols=None, C=None):
        """mu and sigma (n, p_valid) at the given ages and covariates C; cols selects
        valid features."""
        cols = slice(None) if cols is None else cols
        X, W = self._design(np.atleast_1d(age), np.atleast_1d(male), site, C)
        sd = np.exp(W @ self.theta[:, cols])
        mu = X @ self.beta[:, cols] + self.offset[cols] * sd
        return mu, sd * self.scale[cols]

    def base_mu_sigma(self, age, male, cols=None):
        """Unadapted mu and sigma (n, p_valid) at the reference site and median
        covariates (male may be a scalar)."""
        cols = slice(None) if cols is None else cols
        age = np.atleast_1d(age)
        X, W = self._design(age, np.broadcast_to(np.atleast_1d(male), age.shape))
        return X @ self.beta[:, cols], np.exp(W @ self.theta[:, cols])

    def cov_effects(self, C, cols=None):
        """Effects of the covariates C on mu and on log sigma (n, p_valid) relative
        to median covariates; None if the model has no such covariates."""
        cols = slice(None) if cols is None else cols
        Xc = self._cov_columns(C, len(C), self.cov_mean)
        Wc = self._cov_columns(C, len(C), self.cov_sd)
        m = Xc @ self.beta[len(self.beta) - Xc.shape[1]:, cols] if Xc.shape[1] else None
        lw = Wc @ self.theta[len(self.theta) - Wc.shape[1]:, cols] if Wc.shape[1] else None
        return m, lw

    def zscores(self, Y, age, male, site=None, C=None):
        """Normative deviation maps (n, p_valid) at the given ages."""
        mu, sd = self.mu_sigma(age, male, site, C=C)
        return (Y[:, self.valid] - mu) / sd

    def adapt(self, Y, age, male, C=None):
        """Adapt to a new site using its control subjects.

        Removes the mean and scales the SD of the controls' z-scores to 0 and 1
        (the ComBat location/scale step without empirical Bayes).
        """
        self.offset[:] = 0
        self.scale[:] = 1
        z = self.zscores(Y, age, male, C=C)
        self.offset = z.mean(axis=0)
        self.scale = z.std(axis=0, ddof=1)
        self.scale[~(self.scale > 0)] = 1

    def copy(self):
        """Copy that shares the fitted parameters but has its own site adaptation."""
        new = copy.copy(self)
        new.offset, new.scale = self.offset.copy(), self.scale.copy()
        return new

    def get_state(self):
        used = list(self.cov_basis)
        return dict(
            df_mu=self.df_mu, df_sigma=self.df_sigma, cov=asdict(self.cov),
            use_male=bool(self.use_male), site_levels=np.asarray(self.site_levels),
            site_props=np.asarray(self.site_props, dtype=np.float64),
            basis_mu=self.basis_mu.get_state(), basis_sigma=self.basis_sigma.get_state(),
            cov_used=used, cov_basis=[self.cov_basis[nm].get_state() for nm in used],
            cov_median=[self.cov_median[nm] for nm in used],
            cov_mean=list(self.cov_mean), cov_sd=list(self.cov_sd),
            n_features=int(self.n_features), valid=self.valid, beta=self.beta,
            theta=self.theta, offset=self.offset, scale=self.scale)

    @classmethod
    def from_state(cls, s):
        self = cls(s['df_mu'], s['df_sigma'], Covariates.from_state(s['cov']))
        self.use_male = bool(s['use_male'])
        self.site_levels = np.asarray(s['site_levels']).ravel()
        self.site_props = np.asarray(s['site_props'], dtype=np.float64).ravel()
        self.basis_mu = NaturalSpline.from_state(s['basis_mu'])
        self.basis_sigma = NaturalSpline.from_state(s['basis_sigma'])
        self.cov_basis = {nm: NaturalSpline.from_state(b)
                          for nm, b in zip(s['cov_used'], s['cov_basis'])}
        self.cov_median = {nm: float(v) for nm, v in zip(s['cov_used'], s['cov_median'])}
        self.cov_mean, self.cov_sd = list(s['cov_mean']), list(s['cov_sd'])
        self.n_features = int(s['n_features'])
        self.valid = np.asarray(s['valid']).astype(np.int64).ravel()
        self.beta, self.theta = np.asarray(s['beta']), np.asarray(s['theta'])
        self.offset = np.array(s['offset'], dtype=np.float64).ravel()
        self.scale = np.array(s['scale'], dtype=np.float64).ravel()
        return self


def _affine_profile(w, mu, sd):
    """Maximum likelihood a, b (per column) of w ~ N(a + b mu, (b sd)^2)."""
    q = 1 / sd ** 2
    Q = q.sum(axis=0)
    wbar, mbar = (q * w).sum(axis=0) / Q, (q * mu).sum(axis=0) / Q
    dw = w - wbar
    Sww = (q * dw ** 2).sum(axis=0)
    Swm = (q * dw * (mu - mbar)).sum(axis=0)
    c = (Swm + np.sqrt(Swm ** 2 + 4 * w.shape[0] * Sww)) / (2 * Sww)   # c = 1 / b
    return (c * wbar - mbar) / c, 1 / c


def _warp_loglik(x, mu, sd, eps, delta, prior):
    """Log-likelihood per column of x warped by w = sinh(delta asinh(x) - eps):
    sum_i log N(w_i; mu_i, sd_i) + log dw/dx_i without the terms that are
    constant in the warp (-log sd, constants), plus the prior."""
    t = delta * np.arcsinh(x) - eps
    ll = -0.5 * ((np.sinh(t) - mu) / sd) ** 2 + np.log(np.cosh(t))
    return ll.sum(axis=0) + x.shape[0] * np.log(delta) - 0.5 * (eps ** 2 + np.log(delta) ** 2) / prior


def fit_warp(x, mu, sd, eps, delta, n_iter=20, prior=9.0, tol=1e-6):
    """Update the sinh-arcsinh warp w = sinh(delta asinh(x) - eps) of every column of x.

    x (n, p) standardized data; mu, sd (n, p) location and scale of the warped
    data.  Maximizes sum_i [log N(w_i; mu_i, sd_i) + log dw/dx_i] with a weak
    normal prior (variance prior) on eps and log(delta) centred on the identity
    warp.  Damped Newton steps in (eps, eta = log delta) with step halving, so the
    objective never decreases.  Returns eps, delta (p,).
    """
    n = x.shape[0]
    u = np.arcsinh(x)
    eps = np.array(eps, dtype=np.float64)
    eta = np.log(np.asarray(delta, dtype=np.float64))
    ll = _warp_loglik(x, mu, sd, eps, np.exp(eta), prior)
    active = np.arange(x.shape[1])
    for _ in range(n_iter):
        if not active.size:
            break
        ua, m, s2 = u[:, active], mu[:, active], sd[:, active] ** 2
        e, de = eps[active], np.exp(eta[active])
        t = de * ua - e
        c, w = np.cosh(t), np.sinh(t)
        th, sech2 = np.tanh(t), 1 / c ** 2
        r = (w - m) / s2
        c2s = c ** 2 / s2
        g_e = np.sum(r * c - th, axis=0) - e / prior
        g_d = np.sum(-r * ua * c + th * ua, axis=0) + n / de
        h_ee = np.sum(-c2s - r * w + sech2, axis=0) - 1 / prior
        h_ed = np.sum(ua * c2s + r * w * ua - ua * sech2, axis=0)
        h_dd = np.sum(ua ** 2 * (-c2s - r * w + sech2), axis=0) - n / de ** 2
        
        # derivatives with respect to (eps, eta = log delta)
        g_h = de * g_d - eta[active] / prior
        h_eh = de * h_ed
        h_hh = de ** 2 * h_dd + de * g_d - 1 / prior
        
        # make the Hessian negative definite (Levenberg damping), then Newton step
        lam_max = 0.5 * (h_ee + h_hh) + np.sqrt(0.25 * (h_ee - h_hh) ** 2 + h_eh ** 2)
        damp = np.maximum(0, lam_max + 1e-6 * n)
        h11, h12, h22 = h_ee - damp, h_eh, h_hh - damp
        det = h11 * h22 - h12 * h12
        step_e = -(h22 * g_e - h12 * g_h) / det
        step_h = -(-h12 * g_e + h11 * g_h) / det
        acc = np.zeros(active.size, dtype=bool)
        new_e, new_h = e.copy(), eta[active].copy()
        f = np.ones(active.size)
        for _ in range(15):                        # step halving
            todo = ~acc
            if not todo.any():
                break
            ce = np.clip(e[todo] + f[todo] * step_e[todo], -5, 5)
            ch = np.clip(eta[active][todo] + f[todo] * step_h[todo], np.log(0.1), np.log(10))
            cols = active[todo]
            lln = _warp_loglik(x[:, cols], mu[:, cols], sd[:, cols], ce, np.exp(ch), prior)
            better = lln >= ll[cols]
            idx = np.flatnonzero(todo)
            new_e[idx[better]], new_h[idx[better]] = ce[better], ch[better]
            ll[cols[better]] = lln[better]
            acc[idx[better]] = True
            f[idx[~better]] *= 0.5
        change = np.maximum(np.abs(new_e - e), np.abs(new_h - eta[active]))
        eps[active], eta[active] = new_e, new_h
        active = active[acc & (change > tol)]
    return eps, np.exp(eta)


class Warp:
    """Sinh-arcsinh warp of every feature, fitted jointly with the normative model.

    As in warped Bayesian linear regression, every feature is standardized and
    warped by w = sinh(delta asinh(x) - eps), and the warped data follow the
    location-scale model of NormativeModel.  The parameters are estimated by
    block-coordinate ascent on the joint likelihood, per feature until the
    log-likelihood gains less than tol per subject (at most max_outer rounds):

    1. Newton steps for (eps, delta) given location and scale (fit_warp)
    2. the location and scale are carried over to the new warp by the best
       affine rescaling a + b mu, b sd (closed form), then updated by a few
       warm-started RS iterations (fit_location_scale).

    Cheap rounds converge much faster than refitting the normative model to
    convergence each time, and the carry-over removes most of the coupling
    between the warp and the overall location and scale.  (Optimizing the warp
    over the affine rescaling within step 1 instead converges to poorer
    stationary points.)  The warp does not depend on age, so its Jacobian is
    constant in the brain age likelihood; the warped data simply replace the raw
    data in all models.
    """

    def __init__(self, df_mu=5, df_sigma=3, max_outer=40, tol=1e-5, prior=9.0,
                 warp_iter=2, rs_iter=3, chunk=4000, cov=None):
        self.df_mu, self.df_sigma = df_mu, df_sigma
        self.max_outer, self.tol, self.prior = max_outer, tol, prior
        self.warp_iter, self.rs_iter, self.chunk = warp_iter, rs_iter, chunk
        self.cov = cov or Covariates()

    def fit(self, Y, age, male, site=None, verbose=False, C=None):
        finite = np.all(np.isfinite(Y), axis=0)
        sd = np.zeros(Y.shape[1])
        sd[finite] = np.std(Y[:, finite], axis=0)
        self.valid = np.flatnonzero(finite & (sd > 0))
        self.center = np.mean(Y[:, self.valid], axis=0, dtype=np.float64)
        self.sd = sd[self.valid]
        self.eps = np.zeros(self.valid.size)
        self.delta = np.ones(self.valid.size)
        design = NormativeModel(self.df_mu, self.df_sigma, self.cov)
        X, W = design._design(age, male, design._setup(age, male, site, C), C)
        n = len(age)
        ll_start, ll_end, rounds = 0.0, 0.0, 0
        t0 = time.time()
        for c0 in range(0, self.valid.size, self.chunk):
            cols = np.arange(c0, min(c0 + self.chunk, self.valid.size))
            x = self._standardize(Y, cols)
            beta, theta, _ = fit_location_scale(x, X, W)      # identity warp
            e, d = self.eps[cols], self.delta[cols]
            ll = _warp_loglik(x, X @ beta, np.exp(W @ theta), e, d, self.prior) \
                - np.sum(W @ theta, axis=0)
            ll_start += ll.sum()
            active = np.arange(cols.size)
            for it in range(self.max_outer):
                xa, ba, ta = x[:, active], beta[:, active], theta[:, active]
                mu, s = X @ ba, np.exp(W @ ta)
                ea, da = fit_warp(xa, mu, s, e[active], d[active], self.warp_iter, self.prior)
                w = np.sinh(da * np.arcsinh(xa) - ea)
                shift, fac = _affine_profile(w, mu, s)
                ba = ba * fac                          # column 0 of both designs is the intercept
                ba[0] += shift
                ta = ta.copy()
                ta[0] += np.log(fac)
                ba, ta, _ = fit_location_scale(w, X, W, max_iter=self.rs_iter, init=(ba, ta))
                new = _warp_loglik(xa, X @ ba, np.exp(W @ ta), ea, da, self.prior) \
                    - np.sum(W @ ta, axis=0)
                keep = new >= ll[active]               # accept only improvements
                idx = active[keep]
                e[idx], d[idx] = ea[keep], da[keep]
                beta[:, idx], theta[:, idx] = ba[:, keep], ta[:, keep]
                gain = np.where(keep, new - ll[active], 0.0)
                ll[idx] = new[keep]
                active = active[gain > self.tol * n]
                if not active.size:
                    break
            rounds = max(rounds, it + 1)
            self.eps[cols], self.delta[cols] = e, d
            ll_end += ll.sum()
        self.loglik = (ll_start, ll_end)
        if verbose:
            print(f"    warp: {self.valid.size} features, log-likelihood {ll_start:.6g} -> {ll_end:.6g}, "
                  f"up to {rounds} rounds ({time.time() - t0:.0f}s)")
        return self

    def _standardize(self, Y, cols):
        v = self.valid[cols]
        return (np.asarray(Y[:, v], dtype=np.float64) - self.center[cols]) / self.sd[cols]

    def transform(self, Y):
        """Warped copy of Y (float32); features without a warp are left unchanged."""
        out = np.array(Y, dtype=np.float32)
        for c0 in range(0, self.valid.size, self.chunk):
            cols = slice(c0, c0 + self.chunk)
            x = self._standardize(Y, cols)
            out[:, self.valid[cols]] = np.sinh(self.delta[cols] * np.arcsinh(x) - self.eps[cols])
        return out

    def get_state(self):
        return dict(df_mu=self.df_mu, df_sigma=self.df_sigma, cov=asdict(self.cov),
                    chunk=int(self.chunk), valid=self.valid, center=self.center,
                    sd=self.sd, eps=self.eps, delta=self.delta)

    @classmethod
    def from_state(cls, s):
        self = cls(s['df_mu'], s['df_sigma'], chunk=int(s['chunk']),
                   cov=Covariates.from_state(s['cov']))
        self.valid = np.asarray(s['valid']).astype(np.int64).ravel()
        for key in ('center', 'sd', 'eps', 'delta'):
            setattr(self, key, np.asarray(s[key], dtype=np.float64).ravel())
        return self


class ResidualCorrelation:
    """Low-rank plus diagonal correlation R = V diag(lam) V' + diag(psi) of z-scores.

    V, lam: leading eigenvectors/values of the training z-score correlation;
    psi = max(1 - communality, psi_min).  zR^-1z is evaluated with the Woodbury
    identity as sum(D z^2) - ||z A||^2 with D = 1/psi and A (p, rank).
    rank = 0 treats the features as independent.
    """

    def __init__(self, rank=20, psi_min=0.05):
        self.rank, self.psi_min = rank, psi_min

    def fit(self, Z):
        n, p = Z.shape
        k = min(self.rank, n - 1, p)
        if k <= 0:
            self.D = np.ones(p)
            self.A = np.zeros((p, 0))
            return self
        G = Z @ Z.T / n                            # eigenvectors via the n x n Gram matrix
        lam, U = np.linalg.eigh(G)
        lam, U = lam[::-1][:k], U[:, ::-1][:, :k]
        V = (Z.T @ U) / np.sqrt(n * lam)           # (p, k), orthonormal columns
        psi = np.maximum(1 - (V ** 2) @ lam, self.psi_min)
        self.D = 1 / psi
        M = np.diag(1 / lam) + V.T @ (self.D[:, None] * V)
        L = np.linalg.cholesky(np.linalg.inv(M))
        self.A = (self.D[:, None] * V) @ L
        return self

    def get_state(self):
        return dict(rank=int(self.rank), psi_min=float(self.psi_min), D=self.D, A=self.A)

    @classmethod
    def from_state(cls, s):
        self = cls(int(s['rank']), float(s['psi_min']))
        self.D, self.A = np.asarray(s['D'], dtype=np.float64), np.asarray(s['A'], dtype=np.float64)
        return self


def _loglik_grid(model, corr, cols, Y, male_val, grid, C=None, max_bytes=2e8):
    """log p(y | a) (n, G) on the age grid for subjects of one sex.

    The z-scores are z(a) = u (y~ alpha(a) - beta(a)) - kappa, with the data y~
    without the covariate effects on mu, u = exp(-covariate effects on log sigma),
    alpha = 1/sigma(a) and beta = mu(a)/sigma(a) of the model at median covariates
    and kappa = offset/scale of a site adaptation.  sum_j D_j z_j^2 and the
    Woodbury term then expand into matrix products, so that all subjects and a
    chunk of grid ages are handled by a few GEMMs.  The covariate part of
    -sum_j log sigma_j is constant in a and omitted.
    """
    n, p = Y.shape
    k = corr.A.shape[1]
    Yt = np.asarray(Y, dtype=np.float64)
    U = None
    if C is not None:
        m, lw = model.cov_effects(C, cols)
        if m is not None:
            Yt = Yt - m
        if lw is not None:
            U = np.exp(-lw)
    r = model.scale[cols]
    kappa = model.offset[cols] / r
    mref = model.base_mu_sigma(np.median(grid), male_val, cols)[0][0]
    Yt = Yt - mref                                 # centring limits cancellation
    if U is None:
        Y2 = Yt ** 2
    else:
        U2 = U ** 2
        UY, U2Y = U * Yt, U2 * Yt
        U2Y2 = U2Y * Yt
    nb = 1 if U is None else 2                     # arrays of size (chunk, p, k)
    gc = int(max(1, min(len(grid), max_bytes // (8 * p * max(k, 1) * nb))))
    ll = np.empty((n, len(grid)))
    for g0 in range(0, len(grid), gc):
        ag = grid[g0:g0 + gc]
        mu0, s0 = model.base_mu_sigma(ag, male_val, cols)
        sd = r * s0
        alpha = 1 / sd
        beta = (mu0 - mref) * alpha
        Da = corr.D * alpha
        if U is None:
            beta = beta + kappa
            sumsq = (Y2 @ (Da * alpha).T - 2 * (Yt @ (Da * beta).T)
                     + np.sum(corr.D * beta ** 2, axis=1))
        else:
            sumsq = (U2Y2 @ (Da * alpha).T - 2 * (U2Y @ (Da * beta).T)
                     + U2 @ (corr.D * beta ** 2).T - 2 * (UY @ (Da * kappa).T)
                     + 2 * (U @ (corr.D * beta * kappa).T) + np.sum(corr.D * kappa ** 2))
        part = -0.5 * sumsq - np.sum(np.log(sd), axis=1)
        if k:
            B = alpha[:, :, None] * corr.A[None]      # (gc, p, k)
            Yb = Yt if U is None else UY
            P = (Yb @ B.transpose(1, 0, 2).reshape(p, -1)).reshape(n, len(ag), k)
            if U is None:
                P -= beta @ corr.A
            else:
                B = beta[:, :, None] * corr.A[None]
                P -= (U @ B.transpose(1, 0, 2).reshape(p, -1)).reshape(n, len(ag), k)
                P -= kappa @ corr.A
            part += 0.5 * np.sum(P ** 2, axis=2)
        ll[:, g0:g0 + gc] = part
    return ll


def _argmax_refine(ll, grid):
    """Grid argmax refined by a parabola; Laplace SD from its curvature."""
    g = np.argmax(ll, axis=1)
    h = grid[1] - grid[0]
    age = grid[g].astype(np.float64)
    sd = np.full(len(g), np.nan)
    at_bound = (g == 0) | (g == len(grid) - 1)
    i = np.flatnonzero(~at_bound)
    lm, l0, lp = ll[i, g[i] - 1], ll[i, g[i]], ll[i, g[i] + 1]
    curv = lm - 2 * l0 + lp
    ok = curv < 0
    age[i[ok]] += 0.5 * (lm[ok] - lp[ok]) / curv[ok] * h
    sd[i[ok]] = h / np.sqrt(-curv[ok])
    return age, sd, at_bound


def _mahalanobis(model, corr, cols, Y, age, male, C=None):
    """Squared Mahalanobis distance z' R^-1 z of the z-scores at the given ages (n,)."""
    mu, sd = model.mu_sigma(age, male, cols=cols, C=C)
    z = (np.asarray(Y, dtype=np.float64) - mu) / sd
    d2 = z ** 2 @ corr.D
    if corr.A.shape[1]:
        d2 -= np.sum((z @ corr.A) ** 2, axis=1)
    return d2


def deviation_z(d2, k):
    """Normal score of a chi-square(k) distributed distance (Wilson-Hilferty)."""
    return ((np.maximum(d2, 0) / k) ** (1 / 3) - (1 - 2 / (9 * k))) / np.sqrt(2 / (9 * k))


def _gram(A, B, chunk=8192):
    """A @ B.T in float64, accumulated over feature chunks of float32 inputs."""
    G = np.zeros((A.shape[0], B.shape[0]))
    for c in range(0, A.shape[1], chunk):
        G += np.asarray(A[:, c:c + chunk], np.float64) @ np.asarray(B[:, c:c + chunk], np.float64).T
    return G


def _pca(X, n_comp, chunk=8192):
    """Centre (p,) and leading principal axes (k, p) of X (n, p) via the Gram matrix."""
    center = X.mean(axis=0, dtype=np.float64)
    Xc = X - center.astype(X.dtype)
    lam, U = np.linalg.eigh(_gram(Xc, Xc))
    k = min(n_comp, len(lam) - 1)
    lam, U = lam[::-1][:k], U[:, ::-1][:, :k]
    Vt = np.empty((k, X.shape[1]))
    for c in range(0, X.shape[1], chunk):
        Vt[:, c:c + chunk] = U.T @ np.asarray(Xc[:, c:c + chunk], np.float64)
    return center, Vt / np.sqrt(lam)[:, None]


class _Part:
    """Features of the whole brain or of one region with their normative model.

    cols  : data columns used
    mcols : columns of model.beta/theta that belong to this part
    Voxel-wise mode shares one model between all parts; PCA mode fits a model
    to the principal component scores of each part (center, Vt).
    """

    def __init__(self, name, cols, mcols, model, corr, center=None, Vt=None):
        self.name, self.cols, self.mcols = name, cols, mcols
        self.model, self.corr = model, corr
        self.center, self.Vt = center, Vt

    def model_input(self, Y):
        """Data in the form the normative model was fitted to."""
        if self.Vt is None:
            return Y
        return (np.asarray(Y[:, self.cols], np.float64) - self.center) @ self.Vt.T

    def features(self, Y):
        """Features aligned with mcols."""
        return self.model_input(Y)[:, self.model.valid[self.mcols]]

    def get_state(self, model_index):
        return dict(name=self.name, cols=np.asarray(self.cols, dtype=np.int64),
                    mcols=np.asarray(self.mcols, dtype=np.int64), model=int(model_index),
                    corr=self.corr.get_state(), center=self.center, Vt=self.Vt)

    @classmethod
    def from_state(cls, s, models):
        return cls(s['name'], np.asarray(s['cols']).astype(np.int64).ravel(),
                   np.asarray(s['mcols']).astype(np.int64).ravel(), models[s['model']],
                   ResidualCorrelation.from_state(s['corr']), s['center'], s['Vt'])


class NDMBrainAge:
    """Normative model + residual correlation -> likelihood-based brain age.

    pca > 0 : normative models for the leading pca principal component scores
              of the whole brain (and of each region) with a full residual
              correlation between the scores.  PCA is only used here, for the
              brain age and its deviation scores, never for the voxel/vertex-wise
              normative model of the z-maps (VoxelModel).
    pca = 0 : normative model for every voxel/vertex with a low-rank (rank)
              plus diagonal residual correlation.  Its correlation model misses
              most of the dependence between voxels, so the standard errors
              are much too small and the estimates less accurate.  This model
              can be shared with the z-maps (VoxelModel.from_ndm).
    warp    : warp every voxel/vertex with a sinh-arcsinh function fitted jointly
              with a voxel-wise normative model (Warp) before all other steps
    cov     : covariates of all normative models (Covariates, columns of Data.C)
    """

    def __init__(self, df_mu=5, df_sigma=3, pca=100, rank=20, psi_min=0.01,
                 grid_step=0.25, grid_margin=5.0, parcellation=False,
                 atlas_dir=HERE, warp=False, cov=None, verbose=False):
        self.df_mu, self.df_sigma = df_mu, df_sigma
        self.pca, self.rank, self.psi_min = pca, rank, psi_min
        self.grid_step, self.grid_margin = grid_step, grid_margin
        self.parcellation, self.atlas_dir = parcellation, atlas_dir
        self.warp, self.cov, self.verbose = warp, cov or Covariates(), verbose
        self.warper = None

    def _prep(self, d: Data):
        """Data with the training warp applied (unchanged without warp)."""
        return d if self.warper is None else replace(d, Y=self.warper.transform(d.Y))

    def _normative(self, Y, d):
        return NormativeModel(self.df_mu, self.df_sigma, self.cov).fit(
            Y, d.age, d.male, d.site, self.verbose, C=d.C)

    def _fit_pca_part(self, d, name, cols):
        X = d.Y[:, cols]
        keep = np.all(np.isfinite(X), axis=0) & (np.std(X, axis=0) > 0)
        cols = cols[keep]
        center, Vt = _pca(X[:, keep], self.pca)
        part = _Part(name, cols, None, None, None, center, Vt)
        F = part.model_input(d.Y)
        part.model = self._normative(F, d)
        part.mcols = np.arange(part.model.valid.size)
        Z = part.model.zscores(F, d.age, d.male, d.site, C=d.C)
        part.corr = ResidualCorrelation(Z.shape[1], self.psi_min).fit(Z)
        return part

    def fit(self, d: Data):
        if self.warp:
            self.warper = Warp(self.df_mu, self.df_sigma, cov=self.cov).fit(
                d.Y, d.age, d.male, d.site, self.verbose, C=d.C)
            d = self._prep(d)
        groups = [('global', None, np.arange(d.Y.shape[1]))]
        if self.parcellation:
            atlas, regions, names = lobe_atlas(d, self.atlas_dir)
            groups += [(nm, r, np.flatnonzero(atlas == r)) for r, nm in zip(regions, names)]

        self.parts, self.regions, self.region_names = [], [], []
        if not self.pca:
            model = self._normative(d.Y, d)
            Z = model.zscores(d.Y, d.age, d.male, d.site, C=d.C)
        for nm, r, cols in groups:
            if self.pca:
                part = self._fit_pca_part(d, nm, cols)
            else:
                sel = np.flatnonzero(np.isin(model.valid, cols))
                if not sel.size:
                    continue
                part = _Part(nm, model.valid[sel], sel, model,
                             ResidualCorrelation(self.rank, self.psi_min).fit(Z[:, sel]))
            self.parts.append(part)
            if r is not None:
                self.regions.append(r)
                self.region_names.append(nm)

        lo = max(d.age.min() - self.grid_margin, 0.0)
        self.grid = np.arange(lo, d.age.max() + self.grid_margin + self.grid_step / 2,
                              self.grid_step)
        return self

    def _models(self):
        """One part per distinct normative model (voxel-wise mode shares one model)."""
        models = {}
        for part in self.parts:
            models.setdefault(id(part.model), part)
        return list(models.values())

    def copy(self):
        """Copy with its own site adaptation (the fitted parameters are shared)."""
        new = copy.copy(self)
        models, new.parts = {}, []
        for part in self.parts:
            q = copy.copy(part)
            q.model = models.setdefault(id(part.model), part.model.copy())
            new.parts.append(q)
        return new

    def adapt(self, d: Data):
        """Adapt the normative models to a new site using its controls d."""
        d = self._prep(d)
        for part in self._models():
            part.model.adapt(part.model_input(d.Y), d.age, d.male, C=d.C)
        return self

    def adapt_agefree(self, d: Data, n_iter=20, tol=0.01):
        """Adapt the normative models to a new site without the controls' ages.

        Alternates between estimating the brain ages of the controls d with the
        current adaptation and adapting location and scale of their z-scores at
        these estimated ages.  A site effect along the aging direction cannot be
        separated from a brain age offset without ages; it remains as an offset
        of the estimates, which correct_age(..., 'offset') removes.
        """
        d = self._prep(d)
        for part in self._models():
            F, X = part.model_input(d.Y), part.features(d.Y)
            prev = None
            for _ in range(n_iter):
                age = self._estimate(part, X, d.male, d.C)[0]
                part.model.adapt(F, age, d.male, C=d.C)
                if prev is not None and np.max(np.abs(age - prev)) < tol:
                    break
                prev = age
        return self

    def _estimate(self, part, X, male, C=None, chunk=256):
        """Maximum likelihood age, Laplace SD and boundary flag for the features X
        of a part (C: covariates of the subjects)."""
        n = len(X)
        age, sd, bound = np.full(n, np.nan), np.full(n, np.nan), np.zeros(n, dtype=bool)
        for male_val in np.unique(male):
            idx = np.flatnonzero(male == male_val)
            for c0 in range(0, idx.size, chunk):
                sub = idx[c0:c0 + chunk]
                ll = _loglik_grid(part.model, part.corr, part.mcols, X[sub], male_val,
                                  self.grid, None if C is None else C[sub])
                age[sub], sd[sub], bound[sub] = _argmax_refine(ll, self.grid)
        return age, sd, bound

    def predict(self, d: Data, chunk=256):
        """Brain age and non-aging deviation of every subject.

        Returns a dict with (n,) values for the whole brain and (n, n_regions)
        values with the prefix 'regional' ('regional' itself is the brain age):
        'age'           maximum likelihood (brain) age
        'sd'            its Laplace standard error
        'at_bound'      estimate at the boundary of the age grid
        'deviation'     non-aging deviation: normal score of the Mahalanobis distance
                        of the subject's z-scores at its brain age, i.e. the
                        atypicality that an older or younger brain does not explain
        'deviation_age' the same at chronological age (total deviation)
        """
        d = self._prep(d)
        n, n_parts = d.n, len(self.parts)
        res = {k: np.full((n, n_parts), np.nan) for k in ('age', 'sd', 'deviation', 'deviation_age')}
        res['at_bound'] = np.zeros((n, n_parts), dtype=bool)
        ok = np.flatnonzero(np.isfinite(d.age) & (d.age > 0))
        Cs = (lambda idx: None) if d.C is None else (lambda idx: d.C[idx])
        for q, part in enumerate(self.parts):
            X = part.features(d.Y)
            if not np.all(np.isfinite(X)):
                raise ValueError(f"{d.name}: non-finite values in test data.")
            age, res['sd'][:, q], res['at_bound'][:, q] = self._estimate(part, X, d.male, d.C, chunk)
            res['age'][:, q] = age
            k = part.mcols.size
            for c0 in range(0, n, chunk):
                sub = np.arange(c0, min(c0 + chunk, n))
                res['deviation'][sub, q] = deviation_z(_mahalanobis(
                    part.model, part.corr, part.mcols, X[sub], age[sub], d.male[sub], Cs(sub)), k)
            for c0 in range(0, ok.size, chunk):
                sub = ok[c0:c0 + chunk]
                res['deviation_age'][sub, q] = deviation_z(_mahalanobis(
                    part.model, part.corr, part.mcols, X[sub], d.age[sub], d.male[sub], Cs(sub)), k)
        out = {key: v[:, 0] for key, v in res.items()}
        out['regional'] = res['age'][:, 1:]
        out.update({f'regional_{key}': res[key][:, 1:] for key in ('sd', 'deviation', 'deviation_age')})
        return out

    def get_state(self):
        parts = self._models()
        index = {id(p.model): i for i, p in enumerate(parts)}
        return dict(
            df_mu=self.df_mu, df_sigma=self.df_sigma, pca=self.pca, rank=self.rank,
            psi_min=self.psi_min, grid_step=self.grid_step, grid_margin=self.grid_margin,
            parcellation=bool(self.parcellation), warp=bool(self.warp), cov=asdict(self.cov),
            grid=self.grid, regions=[int(r) for r in self.regions],
            region_names=list(self.region_names),
            warper=None if self.warper is None else self.warper.get_state(),
            normative=[p.model.get_state() for p in parts],
            parts=[p.get_state(index[id(p.model)]) for p in self.parts])

    @classmethod
    def from_state(cls, s, atlas_dir=HERE):
        self = cls(s['df_mu'], s['df_sigma'], s['pca'], s['rank'], s['psi_min'],
                   s['grid_step'], s['grid_margin'], s['parcellation'], atlas_dir, s['warp'],
                   Covariates.from_state(s['cov']))
        self.warper = None if s['warper'] is None else Warp.from_state(s['warper'])
        models = [NormativeModel.from_state(m) for m in s['normative']]
        self.parts = [_Part.from_state(p, models) for p in s['parts']]
        self.regions = [int(r) for r in s['regions']]
        self.region_names = list(s['region_names'])
        self.grid = np.asarray(s['grid'], dtype=np.float64).ravel()
        return self


class VoxelModel:
    """Voxel/vertex-wise normative model: the pure normative model, without PCA.

    A location-scale model for every voxel/vertex (age, sex, site and covariates
    as in NormativeModel), optionally after a sinh-arcsinh warp (warp), which
    gives z-maps.  from_ndm() shares the model of an NDMBrainAge with pca=0,
    which is then fitted only once.
    """

    def __init__(self, df_mu=5, df_sigma=3, warp=False, cov=None, verbose=False):
        self.df_mu, self.df_sigma = df_mu, df_sigma
        self.warp, self.cov, self.verbose = warp, cov or Covariates(), verbose
        self.model = self.warper = None
        self.shared = False          # model and warp belong to an NDMBrainAge

    def fit(self, d: Data, warper=None):
        """Fit to the training data d; warper: a fitted Warp of d to reuse."""
        if self.warp:
            self.warper = warper or Warp(self.df_mu, self.df_sigma, cov=self.cov).fit(
                d.Y, d.age, d.male, d.site, self.verbose, C=d.C)
            d = replace(d, Y=self.warper.transform(d.Y))
        self.model = NormativeModel(self.df_mu, self.df_sigma, self.cov).fit(
            d.Y, d.age, d.male, d.site, self.verbose, C=d.C)
        return self

    @classmethod
    def from_ndm(cls, est):
        """The voxel/vertex-wise normative model of an NDMBrainAge with pca=0."""
        if est.pca:
            raise ValueError("Only an NDMBrainAge with pca=0 has a voxel/vertex-wise model.")
        self = cls(est.df_mu, est.df_sigma, est.warp, est.cov, est.verbose)
        self.model, self.warper, self.shared = est.parts[0].model, est.warper, True
        return self

    def zmaps(self, d: Data, ctrl=None, ctrl_age=None, parcellation=False,
              atlas_dir=HERE, chunk=256):
        """z-maps of the subjects d at chronological age.

        Given the indices ctrl of control subjects of d (and their ages ctrl_age,
        by default their chronological ages), a copy of the model is first
        adapted to the site of d.  Returns Z (n, n_features) as float32, with NaN
        for features without a model and for subjects without a valid age, and,
        with parcellation, the lobe-wise mean z (n, n_regions), the region labels
        and their names.
        """
        Y = d.Y if self.warper is None else self.warper.transform(d.Y)
        Cs = (lambda idx: None) if d.C is None else (lambda idx: d.C[idx])
        model = self.model.copy()
        if ctrl is not None:
            model.adapt(Y[ctrl], d.age[ctrl] if ctrl_age is None else ctrl_age,
                        d.male[ctrl], C=Cs(ctrl))
        Z = np.full(d.Y.shape, np.nan, dtype=np.float32)
        ok = np.flatnonzero(np.isfinite(d.age) & (d.age > 0))
        for c0 in range(0, ok.size, chunk):
            sub = ok[c0:c0 + chunk]
            Z[np.ix_(sub, model.valid)] = model.zscores(Y[sub], d.age[sub], d.male[sub],
                                                        C=Cs(sub))
        out = dict(Z=Z)
        if parcellation:
            atlas, regions, names = lobe_atlas(d, atlas_dir)
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', RuntimeWarning)
                out['regional_z'] = np.column_stack([np.nanmean(Z[:, atlas == r], axis=1)
                                                     for r in regions])
            out['regions'] = np.array(regions)
            out['region_names'] = np.array(names, dtype=object)
        return out

    def get_state(self):
        return dict(df_mu=self.df_mu, df_sigma=self.df_sigma, warp=bool(self.warp),
                    cov=asdict(self.cov), shared=bool(self.shared),
                    model=None if self.shared else self.model.get_state(),
                    warper=None if self.shared or self.warper is None else self.warper.get_state())

    @classmethod
    def from_state(cls, s, ndm=None):
        self = cls(s['df_mu'], s['df_sigma'], s['warp'], Covariates.from_state(s['cov']))
        if s['shared']:
            if ndm is None:
                raise ValueError("The voxel-wise model refers to a missing brain age model.")
            self.model, self.warper, self.shared = ndm.parts[0].model, ndm.warper, True
        else:
            self.model = NormativeModel.from_state(s['model'])
            self.warper = None if s['warper'] is None else Warp.from_state(s['warper'])
        return self


# ---------------------------------------------------------------------------
# Saving and loading models
# ---------------------------------------------------------------------------

@dataclass
class ModelSet:
    """Models fitted to one input (segmentation/resolution/smoothing/surface
    measure); info describes its feature space and is checked against new data."""
    info: dict
    brainage: NDMBrainAge | None = None
    voxel: VoxelModel | None = None


def _model_key(name):
    """Model part of a file name, e.g. 's4rp1_4mm' of 's4rp1_4mm_IXI547_CAT12.9.mat'."""
    base = os.path.basename(name.split('+')[0])
    m = re.match(r'(.*?_\d+mm)_', base)
    return m.group(1) if m else base.split('_')[0]


def _ind_hash(ind):
    if ind is None:
        return None
    return hashlib.sha1(np.asarray(ind, dtype=np.int64).ravel().tobytes()).hexdigest()


def data_info(d: Data):
    """Feature space of the data d, checked against new data by check_compatible()."""
    return dict(name=d.name, key=_model_key(d.name), n_features=int(d.Y.shape[1]),
                is_surf=bool(d.is_surf), res=d.res,
                dim=None if d.dim is None else [int(x) for x in d.dim],
                n_vertices=None if d.ind is None else int(np.size(d.ind)),
                ind_sha1=_ind_hash(d.ind))


def check_compatible(info, d: Data):
    """Raise an error if the data d do not match the feature space info of a model."""
    err = []
    if d.Y.shape[1] != info['n_features']:
        err.append(f"{d.Y.shape[1]} features instead of {info['n_features']}")
    if bool(d.is_surf) != info['is_surf']:
        err.append("surface data" if d.is_surf else "volume data")
    if info['res'] and d.res and d.res != info['res']:
        err.append(f"resolution {d.res} mm instead of {info['res']} mm")
    if info['dim'] and d.dim is not None and [int(x) for x in d.dim] != info['dim']:
        err.append(f"dimensions {[int(x) for x in d.dim]} instead of {info['dim']}")
    if info['ind_sha1'] and _ind_hash(d.ind) != info['ind_sha1']:
        err.append("other vertex indices (ind)")
    if err:
        raise ValueError(f"{d.name} does not match the model of {info['name']}: "
                         + '; '.join(err) + '.')
    if _model_key(d.name) != info['key']:
        _note(f"Note: {d.name} is applied to the model of {info['name']}.")


def _uses_sex(s: ModelSet):
    """True if the models of s use sex."""
    m = s.brainage.parts[0].model if s.brainage is not None else s.voxel.model
    return bool(m.use_male)


def _used_covariates(sets):
    """Covariates used by the normative models of sets: (mean, log SD) name lists."""
    mean, sd = {}, {}
    for s in sets:
        models = [p.model for p in s.brainage._models()] if s.brainage is not None else []
        models += [s.voxel.model] if s.voxel is not None else []
        for m in models:
            mean.update(dict.fromkeys(m.cov_mean))
            sd.update(dict.fromkeys(m.cov_sd))
    return list(mean), list(sd)


def _git_version():
    """Short git commit of BA_ndm.py, if available."""
    try:
        import subprocess
        out = subprocess.run(['git', '-C', HERE, 'rev-parse', '--short', 'HEAD'],
                             capture_output=True, text=True, timeout=5)
        return out.stdout.strip() or None
    except Exception:
        return None


def model_description(train, kw, contents, warp_zmaps=False, age_range=(0, np.inf),
                      sets=None):
    """Description of the models sets fitted to train (list of Data with the
    training subjects used): settings, training sample and covariates (requested
    and, with sets, used)."""
    d0 = train[0]
    cov = kw.get('cov') or Covariates()
    sites = d0.name.split('+')
    settings = {k: v for k, v in kw.items() if k not in ('cov', 'atlas_dir', 'verbose')}
    settings.update(warp_zmaps=bool(warp_zmaps),
                    age_range=[float(x) if np.isfinite(x) else None for x in age_range])
    desc = dict(format=MODEL_FORMAT, version=MODEL_VERSION,
                created=time.strftime('%Y-%m-%d %H:%M:%S'), ba_ndm=_git_version(),
                contents=contents, settings=settings,
                training=dict(files=[d.name for d in train], n=int(d0.n),
                              age_min=float(np.min(d0.age)), age_max=float(np.max(d0.age)),
                              sites=sites,
                              n_per_site=[int(np.sum(d0.site == s)) for s in range(len(sites))],
                              both_sexes=bool(np.unique(d0.male).size > 1)),
                covariates=None)
    if cov.names:
        desc['covariates'] = dict(asdict(cov), min=[float(x) for x in np.min(d0.C, axis=0)],
                                  max=[float(x) for x in np.max(d0.C, axis=0)])
        if sets is not None:
            desc['covariates']['used_mean'], desc['covariates']['used_sd'] = _used_covariates(sets)
    return desc


def _to_tree(obj, path, arrays):
    """JSON tree of a nested state; arrays are replaced by references and collected."""
    if isinstance(obj, np.ndarray):
        key = '/'.join(map(str, path))
        arrays[key] = obj
        return {'__array__': key, 'shape': list(obj.shape), 'dtype': obj.dtype.str}
    if isinstance(obj, dict):
        return {str(k): _to_tree(v, path + [k], arrays) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_to_tree(v, path + [i], arrays) for i, v in enumerate(obj)]
    if isinstance(obj, np.generic):
        return obj.item()
    return obj


def _from_tree(tree, get):
    """Nested state of a JSON tree; get(key) returns the array of a reference."""
    if isinstance(tree, dict):
        if '__array__' in tree:
            a = np.asarray(get(tree['__array__']))
            return a.astype(np.dtype(tree['dtype']), copy=False).reshape(tree['shape'])
        return {k: _from_tree(v, get) for k, v in tree.items()}
    if isinstance(tree, list):
        return [_from_tree(v, get) for v in tree]
    return tree


def _to_mat(obj):
    """Nested state for savemat: dicts -> structs, lists -> cell arrays or numeric
    vectors, None -> []."""
    if isinstance(obj, dict):
        return {k: _to_mat(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        if obj and all(isinstance(v, (int, float, np.number)) and not isinstance(v, bool)
                       for v in obj):
            return np.asarray(obj, dtype=np.float64)
        cell = np.empty(len(obj), dtype=object)
        for i, v in enumerate(obj):
            cell[i] = _to_mat(v)
        return cell
    return np.zeros((0, 0)) if obj is None else obj


def _mat_get(x, key):
    """Array at the path key ('a/0/b') of a struct read with struct_as_record=False."""
    for p in key.split('/'):
        if p.isdigit():                            # element of a cell array
            if isinstance(x, np.ndarray):
                x = x.ravel()[int(p)]
        else:
            while isinstance(x, np.ndarray) and x.dtype == object:
                x = x.ravel()[0]
            x = getattr(x, p)
    return x


def save_models(path, sets, description):
    """Save fitted models (list of ModelSet) with their description.

    <name>.mat : MATLAB struct NDMmodel (the description and, for every input,
                 info, brainage and voxel) and the JSON string meta_json
    <name>.npz : NumPy archive with the arrays and the JSON string __meta__
    The description and the feature space of every input are also written to
    <name>.json; load_models() uses the copy in the model file.
    """
    ext = os.path.splitext(path)[1].lower()
    if ext not in ('.mat', '.npz'):
        raise ValueError(f"{path}: the model file needs the extension .mat or .npz.")
    state = dict(description, models=[dict(
        info=s.info, brainage=None if s.brainage is None else s.brainage.get_state(),
        voxel=None if s.voxel is None else s.voxel.get_state()) for s in sets])
    arrays = {}
    meta = json.dumps(_to_tree(state, [], arrays))
    if ext == '.npz':
        np.savez_compressed(path, __meta__=np.array(meta), **arrays)
    else:
        from scipy.io import savemat
        savemat(path, {'NDMmodel': _to_mat(state), 'meta_json': meta},
                do_compression=True, long_field_names=True)
    info = dict(description, models=[dict(s.info, brainage=s.brainage is not None,
                                          voxel=s.voxel is not None) for s in sets])
    json_path = os.path.splitext(path)[0] + '.json'
    with open(json_path, 'w') as f:
        json.dump(_to_tree(info, [], {}), f, indent=2)
    print(f"Saved models to {path} and their description to {json_path}")


def load_models(path, atlas_dir=HERE):
    """Models saved by save_models(): (list of ModelSet, description)."""
    ext = os.path.splitext(path)[1].lower()
    if ext == '.npz':
        with np.load(path, allow_pickle=False) as f:
            if '__meta__' not in f.files:
                raise ValueError(f"{path} is not a BA_ndm model file.")
            state = _from_tree(json.loads(str(f['__meta__'])), lambda key: f[key])
    elif ext == '.mat':
        from scipy.io import loadmat
        m = loadmat(path, squeeze_me=False, struct_as_record=False)
        if 'meta_json' not in m or 'NDMmodel' not in m:
            raise ValueError(f"{path} is not a BA_ndm model file.")
        tree = json.loads(str(np.asarray(m['meta_json']).ravel()[0]))
        state = _from_tree(tree, lambda key: _mat_get(m['NDMmodel'], key))
    else:
        raise ValueError(f"{path}: the model file needs the extension .mat or .npz.")
    if state.get('format') != MODEL_FORMAT:
        raise ValueError(f"{path} is not a BA_ndm model file.")
    if state.get('version', 0) > MODEL_VERSION:
        raise ValueError(f"{path} was saved by a newer version of BA_ndm.py.")
    sets = []
    for s in state.pop('models'):
        est = None if s['brainage'] is None else NDMBrainAge.from_state(s['brainage'], atlas_dir)
        vm = None if s['voxel'] is None else VoxelModel.from_state(s['voxel'], est)
        sets.append(ModelSet(s['info'], est, vm))
    return sets, state


def _describe(desc, path):
    """One-line summary of a model description."""
    tr = desc['training']
    s = (f"{path}: {desc['contents']}, fitted {desc['created']} to {tr['n']} subjects "
         f"aged {tr['age_min']:.1f}-{tr['age_max']:.1f} ({', '.join(tr['sites'])})")
    cov = desc.get('covariates')
    if cov:
        s += (f"; covariates in the mean: {', '.join(cov.get('used_mean', cov['mean'])) or '-'}, "
              f"in the log SD: {', '.join(cov.get('used_sd', cov['sd'])) or '-'}")
    return s


# ---------------------------------------------------------------------------
# GPR BrainAGE replica (BA_gpr.m / BA_gpr_core.m) and trend correction
# ---------------------------------------------------------------------------

def gpr_baseline(Y_train, age_train, Y_test, mean_hyp=100.0, lik_hyp=-1.0):
    """Predicted age of the GPR BrainAGE (linear kernel, PCA, no dropout).

    Scaling to 0..1, PCA with n-1 components, rescaling of the scores to 0..1
    and a GP with linear covariance and constant prior mean as in BA_gpr.m.
    PCA signs differ from MATLAB's svd, which changes the rescaled scores
    slightly, so results are close to but not identical with BA_gpr_ui.m.
    """
    mn, mx = float(Y_train.min()), float(Y_train.max())
    center = Y_train.mean(axis=0, dtype=np.float64).astype(np.float32)
    Xtr = (Y_train - center) / np.float32(mx - mn)
    Xte = (Y_test - center) / np.float32(mx - mn)
    lam, U = np.linalg.eigh(_gram(Xtr, Xtr))
    k = min(Xtr.shape) - 1
    lam, U = lam[::-1][:k], U[:, ::-1][:, :k]
    Mtr = U * np.sqrt(lam)                         # training scores (U S)
    Mte = _gram(Xte, Xtr) @ U / np.sqrt(lam)
    lo, hi = Mtr.min(), Mtr.max()
    Mtr, Mte = (Mtr - lo) / (hi - lo), (Mte - lo) / (hi - lo)
    K = Mtr @ Mtr.T
    alpha = np.linalg.solve(K + np.exp(2 * lik_hyp) * np.eye(len(K)),
                            np.asarray(age_train, float) - mean_hyp)
    return mean_hyp + Mte @ Mtr.T @ alpha


def correct_age(pred, age, ctrl, method):
    """Predicted age (n,) or (n, k) corrected with the control subjects ctrl (indices).

    'offset' subtracts the median BrainAGE of the controls (the offset-only
    correction of BA_gpr_ui.m), 'trend' subtracts a linear trend of BrainAGE on
    age estimated on the controls (trend_method = 1, trend_degree = 1),
    'none' returns the estimate unchanged.  Columns are corrected separately.
    """
    P = np.array(pred, dtype=np.float64).reshape(len(age), -1)
    ba = P - age[:, None]
    if method == 'offset':
        P -= np.nanmedian(ba[ctrl], axis=0)
    elif method == 'trend':
        for j in range(P.shape[1]):
            P[:, j] -= np.polyval(np.polyfit(age[ctrl], ba[ctrl, j], 1), age)
    elif method != 'none':
        raise ValueError(f"Unknown correction {method!r}.")
    return P.reshape(np.shape(pred))


def trend_correct(pred, age, ctrl):
    """BrainAGE with linear trend correction estimated on the controls ctrl."""
    return correct_age(pred, age, ctrl, 'trend') - age


# ---------------------------------------------------------------------------
# Ensemble and evaluation
# ---------------------------------------------------------------------------

def ensemble_weights(pred, age, method='gls'):
    """Weights summing to one that combine the age estimates of several models.

    pred : (n, m) estimates for control subjects
    'mean' equal weights, 'mae' 1/MAE^2 (BA_gpr_ui.m D.ensemble = 5),
    'gls'  minimum mean squared error under the sum-to-one constraint,
           w ~ C^-1 1 with C = E[e e'] of the errors e = pred - age.
    """
    m = pred.shape[1]
    E = pred - age[:, None]
    if m == 1 or method == 'mean':
        w = np.ones(m)
    elif method == 'mae':
        w = 1 / np.mean(np.abs(E), axis=0) ** 2
    elif method == 'gls':
        w = np.linalg.solve(E.T @ E / len(E), np.ones(m))
    else:
        raise ValueError(f"Unknown ensemble method {method!r}.")
    return w / w.sum()


def metrics(ba, age, sd=None):
    """Summary of BrainAGE values (predicted - chronological age) of controls."""
    out = dict(MAE=np.mean(np.abs(ba)), RMSE=np.sqrt(np.mean(ba ** 2)),
               r_pred_age=np.corrcoef(ba + age, age)[0, 1],
               r_BA_age=np.corrcoef(ba, age)[0, 1], mean_BA=np.mean(ba))
    if sd is not None:
        ok = np.isfinite(sd)
        out['cover95'] = np.mean(np.abs(ba[ok]) < 1.96 * sd[ok]) if ok.any() else np.nan
    return out


def _print_table(rows):
    cols = ['MAE', 'RMSE', 'r_pred_age', 'r_BA_age', 'mean_BA', 'cover95']
    w = max(len(r[0]) for r in rows) + 2
    print(f"{'':{w}s}" + ''.join(f"{c:>11s}" for c in cols))
    for name, m in rows:
        print(f"{name:{w}s}" + ''.join(
            f"{m[c]:11.3f}" if c in m and np.isfinite(m[c]) else f"{'':11s}" for c in cols))


def _parse_index(spec, n):
    """MATLAB-style 1-based index list '1:108,150,160:170' -> 0-based array."""
    idx = []
    for part in spec.split(','):
        a, _, b = part.partition(':')
        idx.extend(range(int(a) - 1, int(b) if b else int(a)))
    idx = np.array(idx)
    if idx.min() < 0 or idx.max() >= n:
        raise ValueError(f"Index {spec} out of range 1..{n}.")
    return idx


# ---------------------------------------------------------------------------
# Drivers
# ---------------------------------------------------------------------------

def _valid_train(d: Data, age_range):
    """Subjects usable for training: valid age within age_range, no missing covariates."""
    ok = np.isfinite(d.age) & (d.age > 0) & (d.age >= age_range[0]) & (d.age <= age_range[1])
    if d.C is not None:
        ok &= np.all(np.isfinite(d.C), axis=1)
    return ok


def _check_same_subjects(datas):
    for d in datas[1:]:
        if d.n != datas[0].n or not np.allclose(d.age, datas[0].age, equal_nan=True):
            raise ValueError(f"{d.name} and {datas[0].name} do not contain the same subjects.")


def _combine(res, ages, ctrl, method):
    """Ensemble of the global and regional estimates of all models."""
    pred = np.column_stack([r['age'] for r in res])
    w = ensemble_weights(pred[ctrl], ages[ctrl], method)
    out = dict(weights=w, age=pred @ w)
    regions = sorted({r_ for r in res for r_ in r['regions']})
    reg = np.full((len(ages), len(regions)), np.nan)
    for j, r_ in enumerate(regions):
        cols = [r['regional'][:, r['regions'].index(r_)] for r in res if r_ in r['regions']]
        P = np.column_stack(cols)
        reg[:, j] = P @ ensemble_weights(P[ctrl], ages[ctrl], method)
    out['regional'], out['regions'] = reg, regions
    return out


def cross_validate(datas, kfold=10, seed=0, age_range=(0, np.inf), gpr=True,
                   ensemble='gls', **kw):
    """k-fold cross-validation of the NDM brain age (and GPR baseline)."""
    _check_same_subjects(datas)
    d0 = datas[0]
    ok = _valid_train(d0, age_range)
    if np.sum(~ok):
        print(f"{int(np.sum(~ok))} subject(s) excluded (invalid age, outside age range or missing covariates).")
        
    # age-stratified folds with random tie-breaking
    rng = np.random.default_rng(seed)
    idx = np.flatnonzero(ok)
    order = idx[np.lexsort((rng.random(idx.size), d0.age[idx]))]
    fold = np.full(d0.n, -1)
    fold[order] = np.arange(order.size) % kfold

    n = d0.n
    res = []
    for d in datas:
        res.append(dict(name=d.name, age=np.full(n, np.nan), sd=np.full(n, np.nan),
                        deviation=np.full(n, np.nan), deviation_age=np.full(n, np.nan),
                        at_bound=np.zeros(n, bool), regional=None, regions=[],
                        region_names=[], gpr=np.full(n, np.nan)))
    for f in range(kfold):
        tr, te = np.flatnonzero(fold != f) , np.flatnonzero(fold == f)
        tr = tr[ok[tr]]
        print(f"fold {f + 1}/{kfold}: {tr.size} training, {te.size} test subjects", flush=True)
        for d, r in zip(datas, res):
            t0 = time.time()
            est = NDMBrainAge(**kw).fit(d.subset(tr))
            out = est.predict(d.subset(te))
            for key in ('age', 'sd', 'at_bound', 'deviation', 'deviation_age'):
                r[key][te] = out[key]
            if r['regional'] is None:
                for key in ('regional', 'regional_deviation', 'regional_deviation_age'):
                    r[key] = np.full((n, len(est.regions)), np.nan)
                r['regions'], r['region_names'] = est.regions, est.region_names
            for key in ('regional', 'regional_deviation', 'regional_deviation_age'):
                r[key][te] = out[key]
            if gpr:
                r['gpr'][te] = gpr_baseline(d.Y[tr], d.age[tr], d.Y[te])
            print(f"  {d.name}: {time.time() - t0:.1f}s", flush=True)

    if gpr:
        for r in res:
            r['gpr_ba'] = trend_correct(r['gpr'], d0.age, np.flatnonzero(ok))
    return _summarize(datas, res, d0.age, ok, ensemble, gpr, fold=fold)


def fit_models(train, kw, brainage=True, voxel=False, warp_zmaps=False):
    """Fit models to every training input (list of Data with the subjects to use).

    brainage   : NDMBrainAge(**kw) for the brain age
    voxel      : VoxelModel for the z-maps (warped with warp_zmaps); with pca=0 and
                 the same warp, the model of the brain age is shared (fitted once)
    Returns a list of ModelSet.
    """
    sets = []
    for d in train:
        t0 = time.time()
        est = NDMBrainAge(**kw).fit(d) if brainage else None
        vm = None
        if voxel:
            if est is not None and not est.pca and est.warp == warp_zmaps:
                vm = VoxelModel.from_ndm(est)
            else:
                vm = VoxelModel(kw.get('df_mu', 5), kw.get('df_sigma', 3), warp_zmaps,
                                kw.get('cov')).fit(
                    d, warper=est.warper if est is not None and warp_zmaps else None)
        sets.append(ModelSet(data_info(d), est, vm))
        print(f"  {d.name}: fitted in {time.time() - t0:.1f}s", flush=True)
    return sets


def apply_models(sets, test, adjust=None, correction='offset', ensemble='gls',
                 zmaps_out=False, normative_only=False, parcellation=False,
                 train=None, gpr=True):
    """Apply fitted models (list of ModelSet) to the test inputs in the same order.

    adjust         : indices of the control subjects of the test sample
    correction     : use of the controls for the test site.  Brain age: 'offset'
                     subtracts their median BrainAGE, 'adapt' adapts the normative
                     models (location and scale of the z-scores), 'agefree' adapts
                     them without the controls' ages (at their estimated brain ages)
                     and then subtracts their median BrainAGE, 'trend' removes a
                     linear trend of BrainAGE on age, 'none' only estimates the
                     ensemble weights.  The GPR baseline is corrected in the same way
                     ('offset' for 'adapt' and 'agefree').  The z-maps are adapted
                     with the controls unless 'none' (at their estimated brain ages
                     for 'agefree').
    zmaps_out      : also return voxel/vertex-wise z-maps at chronological age
    normative_only : only the z-maps of the voxel/vertex-wise models
    parcellation   : lobe-wise mean z of the z-maps (also if the brain age models
                     were fitted with parcellation)
    train          : training data of the models, needed for the GPR baseline
    """
    _check_same_subjects(test)
    if len(sets) != len(test):
        raise ValueError(f"{len(test)} test input(s) for {len(sets)} model(s).")
    for s, d in zip(sets, test):
        check_compatible(s.info, d)
        if (normative_only or zmaps_out) and s.voxel is None:
            raise ValueError(f"No voxel/vertex-wise normative model for {d.name}.")
        if not normative_only and s.brainage is None:
            raise ValueError(f"No brain age model for {d.name}.")
    n = test[0].n
    age = test[0].age
    ok = np.isfinite(age) & (age > 0)
    if adjust is None:
        correction = 'none'
        ctrl = np.flatnonzero(ok)
    else:
        ctrl = np.intersect1d(adjust, np.flatnonzero(ok))
    zc = None if correction == 'none' else ctrl
    cov_names = test[0].cov_names

    if normative_only:
        if correction == 'agefree':
            raise ValueError("The correction 'agefree' needs brain age models.")
        zm = []
        for s, d in zip(sets, test):
            t0 = time.time()
            z = s.voxel.zmaps(d, zc, parcellation=parcellation)
            zm.append(dict(z, model=d.name, ind=d.ind, age=age, male=d.male))
            print(f"  {d.name}: {time.time() - t0:.1f}s", flush=True)
        if zc is None:
            print("No adaptation to the test site: z-maps at the reference site of the "
                  "training data.")
        else:
            print(f"z-maps adapted to the test site with {zc.size} controls.")
        out = dict(age=age, male=test[0].male, models=np.array([d.name for d in test], dtype=object),
                   correction=correction, ind_control=ctrl + 1, normative_only=True, zmaps=zm)
        if cov_names:
            out['covariates'] = np.array(cov_names, dtype=object)
        return out

    post = correction if correction in ('offset', 'trend') else 'none'
    if correction == 'agefree':
        post = 'offset'
    gpr_post = 'offset' if correction in ('adapt', 'agefree') else post
    gpr = gpr and train is not None
    res = []
    for i, (s, dte) in enumerate(zip(sets, test)):
        t0 = time.time()
        est = s.brainage.copy()
        if correction == 'adapt':
            est.adapt(dte.subset(ctrl))
        elif correction == 'agefree':
            est.adapt_agefree(dte.subset(ctrl))
        out = est.predict(dte)
        r = dict(name=dte.name, age=correct_age(out['age'], age, ctrl, post),
                 sd=out['sd'], at_bound=out['at_bound'],
                 regional=correct_age(out['regional'], age, ctrl, post), regions=est.regions,
                 region_names=est.region_names, gpr=np.full(n, np.nan),
                 offset=np.nanmedian(out['age'][ctrl] - age[ctrl]))
        for key in ('deviation', 'deviation_age', 'regional_deviation', 'regional_deviation_age'):
            r[key] = out[key]
        if gpr:
            r['gpr'] = gpr_baseline(train[i].Y, train[i].age, dte.Y)
            r['gpr_ba'] = correct_age(r['gpr'], age, ctrl, gpr_post) - age
        if zmaps_out:
            r['zmaps'] = s.voxel.zmaps(dte, zc, out['age'][zc] if correction == 'agefree' else None,
                                       parcellation or est.parcellation, est.atlas_dir)
        res.append(r)
        print(f"  {dte.name}: {time.time() - t0:.1f}s", flush=True)
    if adjust is None:
        print("No --adjust given: no correction for the test site and ensemble weights "
              "estimated on all test subjects.")
    else:
        print(f"Correction for the test site: {correction} ({ctrl.size} controls)")
    gpr_label = {'offset': 'offset corr.', 'trend': 'trend corr.', 'none': 'uncorrected'}[gpr_post]
    out = _summarize(test, res, age, ok, ensemble, gpr, ctrl=ctrl, gpr_label=gpr_label)
    out['correction'] = correction
    out['offset'] = np.array([r['offset'] for r in res])  # median BrainAGE of controls before correction
    if cov_names:
        out['covariates'] = np.array(cov_names, dtype=object)
    if zmaps_out:
        out['zmaps'] = [dict(r['zmaps'], model=r['name'], ind=d.ind, age=age, male=d.male)
                        for r, d in zip(res, test)]
    return out


def train_test(train, test, adjust=None, correction='offset', age_range=(0, np.inf),
               gpr=True, ensemble='gls', zmaps_out=False, warp_zmaps=False, **kw):
    """Train on one sample and predict another (fit_models and apply_models).

    Arguments as for apply_models; kw are passed to NDMBrainAge.  Without sex in
    the test data, sex is dropped from all models.
    """
    if not all(d.has_male for d in test):
        print("Test data contain no sex: sex is not used in the normative models.")
        train = [replace(d, male=np.zeros_like(d.male)) for d in train]
        test = [replace(d, male=np.zeros_like(d.male)) for d in test]
    sel = np.flatnonzero(_valid_train(train[0], age_range))
    train = [d.subset(sel) for d in train]
    sets = fit_models(train, kw, True, zmaps_out, warp_zmaps)
    return apply_models(sets, test, adjust, correction, ensemble, zmaps_out, train=train, gpr=gpr)


def _check_covariates(C, desc, used):
    """Notes on missing test covariates and on values outside the training range
    (used: names of the covariates used by the models)."""
    cov = desc['covariates']
    j_used = [j for j, nm in enumerate(cov['names']) if nm in used]
    n_miss = int(np.sum(~np.all(np.isfinite(C[:, j_used]), axis=1)))
    if n_miss:
        print(f"{n_miss} test subject(s) with missing covariates: the training median is used.")
    for j in j_used:
        nm, lo, hi = cov['names'][j], cov['min'][j], cov['max'][j]
        n_out = int(np.sum((C[:, j] < lo) | (C[:, j] > hi)))
        if n_out:
            print(f"Covariate {nm}: {n_out} test subject(s) outside the training range "
                  f"{lo:.4g} to {hi:.4g} (extrapolated).")


def _check_ages(age, desc):
    """Note on test subjects outside the training age range."""
    lo, hi = desc['training']['age_min'], desc['training']['age_max']
    n_out = int(np.sum((age < lo) | (age > hi)))
    if n_out:
        print(f"{n_out} test subject(s) outside the training age range {lo:.1f} to {hi:.1f}: "
              "the normative models extrapolate linearly.")


def _summarize(datas, res, age, ok, ensemble, gpr, fold=None, ctrl=None, gpr_label='trend corr.'):
    ctrl = np.flatnonzero(ok) if ctrl is None else np.intersect1d(ctrl, np.flatnonzero(ok))
    rows = []
    for r in res:
        rows.append((f"NDM {r['name']}", metrics(r['age'][ctrl] - age[ctrl], age[ctrl],
                                                 r['sd'][ctrl])))
        nb = int(np.sum(r['at_bound'][ctrl]))
        if nb:
            print(f"{r['name']}: {nb} estimate(s) at the boundary of the age grid")
    ens = {m: _combine(res, age, ctrl, m) for m in ('mean', 'mae', 'gls')}
    for m in ('mean', 'mae', 'gls'):
        rows.append((f"NDM ensemble ({m}) w={np.round(ens[m]['weights'], 2)}",
                     metrics(ens[m]['age'][ctrl] - age[ctrl], age[ctrl])))
    if gpr:
        for r in res:
            rows.append((f"GPR {r['name']} ({gpr_label})",
                         metrics(r['gpr_ba'][ctrl], age[ctrl])))
                         
        # ensemble of the trend-corrected GPR predictions (as BA_gpr_ui.m)
        P = np.column_stack([r['gpr_ba'] + age for r in res])
        for m in ('mae', 'gls'):
            w = ensemble_weights(P[ctrl], age[ctrl], m)
            rows.append((f"GPR ensemble ({m}) w={np.round(w, 2)}",
                         metrics(P[ctrl] @ w - age[ctrl], age[ctrl])))
        gpr_ens = P @ ensemble_weights(P[ctrl], age[ctrl], ensemble) - age
    print()
    _print_table(rows)
    if gpr:
        ba_ndm = ens[ensemble]['age'] - age
        print(f"\ncorr(NDM ensemble BrainAGE, GPR ensemble BrainAGE) in controls: "
              f"{np.corrcoef(ba_ndm[ctrl], gpr_ens[ctrl])[0, 1]:.3f}")

    e = ens[ensemble]
    reg_names = _region_names(HERE)
    if e['regions']:
        print("\nRegional NDM BrainAGE (ensemble), controls: mean / MAE / r(BA, age)")
        for j, rid in enumerate(e['regions']):
            ba = e['regional'][ctrl, j] - age[ctrl]
            print(f"  {reg_names.get(rid, str(rid)):28s} {np.mean(ba):7.2f} "
                  f"{np.mean(np.abs(ba)):7.2f} {np.corrcoef(ba, age[ctrl])[0, 1]:+7.2f}")

    out = dict(
        age=age, male=datas[0].male,
        models=np.array([r['name'] for r in res], dtype=object),
        PredictedAge=np.column_stack([r['age'] for r in res]),
        PredictedAge_sd=np.column_stack([r['sd'] for r in res]),
        at_bound=np.column_stack([r['at_bound'] for r in res]),
        PredictedAge_ensemble=e['age'], BrainAGE_ensemble=e['age'] - age,
        weights=e['weights'], ensemble_method=ensemble,
        regions=np.array(e['regions']),
        region_names=np.array([reg_names.get(r_, str(r_)) for r_ in e['regions']], dtype=object),
        BrainAGE_regional_ensemble=e['regional'] - age[:, None] if e['regions'] else np.zeros((len(age), 0)),
        ind_control=ctrl + 1)
    out['BrainAGE'] = out['PredictedAge'] - age[:, None]
    
    # non-aging deviation (at brain age) and total deviation (at chronological age),
    # normal scores; the ensemble is the mean over models
    for key, name in (('deviation', 'Deviation'), ('deviation_age', 'Deviation_age')):
        out[name] = np.column_stack([r[key] for r in res])
        out[name + '_ensemble'] = np.mean(out[name], axis=1)
    dev = out['Deviation_ensemble']
    print(f"\nNon-aging deviation (normal score, ensemble) in controls: "
          f"mean {np.mean(dev[ctrl]):+.2f}, SD {np.std(dev[ctrl]):.2f}")
    if e['regions']:
        n_reg, n_mod = len(e['regions']), len(res)
        R = np.full((len(age), n_reg, n_mod), np.nan)
        Dr = {key: np.full((len(age), n_reg, n_mod), np.nan)
              for key in ('regional_deviation', 'regional_deviation_age')}
        for m, r in enumerate(res):
            for j, rid in enumerate(e['regions']):
                if rid in r['regions']:
                    col = r['regions'].index(rid)
                    R[:, j, m] = r['regional'][:, col] - age
                    for key in Dr:
                        Dr[key][:, j, m] = r[key][:, col]
        out['BrainAGE_regional'] = R
        out['Deviation_regional'] = Dr['regional_deviation']
        out['Deviation_age_regional'] = Dr['regional_deviation_age']
        with np.errstate(all='ignore'):
            out['Deviation_regional_ensemble'] = np.nanmean(Dr['regional_deviation'], axis=2)
    if gpr:
        out['PredictedAge_GPR'] = np.column_stack([r['gpr'] for r in res])
        out['BrainAGE_GPR'] = np.column_stack([r['gpr_ba'] for r in res])
        out['BrainAGE_GPR_ensemble'] = gpr_ens
    if fold is not None:
        out['fold'] = fold + 1
    return out


def save_results(out, prefix):
    """Save the results.

    Brain age: <prefix>.mat (struct NDM) and <prefix>.csv.  z-maps, if computed:
    one <prefix>_zmaps_<model>.mat per model (struct NDMzmap: Z (n, n_features) in
    the feature order of the input Y, regional_z, regions, region_names, age, male,
    model, ind).  With --normative-only, only the z-maps are saved and, with
    regional z, their lobe-wise means as <prefix>.csv.
    """
    from scipy.io import savemat
    out = dict(out)
    zm = out.pop('zmaps', None)
    names = [os.path.splitext(m)[0] for m in out['models']]
    for z, nm in zip(zm or [], names):
        z = {k: (np.zeros(0) if v is None else v) for k, v in z.items()}
        savemat(f'{prefix}_zmaps_{nm}.mat', {'NDMzmap': z}, do_compression=True)
    zfiles = f"{len(zm)} z-map file(s) {prefix}_zmaps_*.mat" if zm else ''
    if out.get('normative_only'):
        saved = [zfiles]
        if zm and 'regional_z' in zm[0]:
            cols = [f"z_{nm}_{str(r).replace(' ', '_')}"
                    for nm, z in zip(names, zm) for r in z['region_names']]
            with open(prefix + '.csv', 'w') as f:
                f.write(','.join(['subject', 'age', 'male'] + cols) + '\n')
                for i in range(len(out['age'])):
                    vals = [out['age'][i], out['male'][i]] + [v for z in zm for v in z['regional_z'][i]]
                    f.write(','.join([str(i + 1)] + [f'{v:.4f}' for v in vals]) + '\n')
            saved.insert(0, prefix + '.csv')
        print("\nSaved " + ' and '.join(saved))
        return
    savemat(prefix + '.mat', {'NDM': out}, do_compression=True)
    cols = ['age', 'male', 'BrainAGE_ensemble', 'Deviation_ensemble', 'Deviation_age_ensemble']
    if 'BrainAGE_GPR_ensemble' in out:
        cols.append('BrainAGE_GPR_ensemble')
    with open(prefix + '.csv', 'w') as f:
        f.write(','.join(['subject'] + cols + [f'BrainAGE_{m}' for m in names]
                         + [f'SD_{m}' for m in names] + [f'Deviation_{m}' for m in names]) + '\n')
        for i in range(len(out['age'])):
            vals = ([out[c][i] for c in cols] + list(out['BrainAGE'][i])
                    + list(out['PredictedAge_sd'][i]) + list(out['Deviation'][i]))
            f.write(','.join([str(i + 1)] + [f'{v:.4f}' for v in vals]) + '\n')
    print(f"\nSaved {prefix}.mat and {prefix}.csv" + (f" and {zfiles}" if zm else ''))


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Brain age and voxel/vertex-wise normative models (NDM prototype).",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    p.add_argument('--train', nargs='+', help="training mat-file per model ('+' joins sites)")
    p.add_argument('--model', help="saved models (--save-model) to apply to --test instead "
                                   "of fitting models to --train")
    p.add_argument('--save-model', help="save the fitted models to <name>.mat (MATLAB struct "
                                        "NDMmodel) or <name>.npz (NumPy), with a JSON "
                                        "description in <name>.json")
    p.add_argument('--test', nargs='+', help="test mat-file per model (same order as --train "
                                             "or the saved models)")
    p.add_argument('--adjust', help="1-based indices of test controls, e.g. 1:108")
    p.add_argument('--correction', choices=['offset', 'adapt', 'agefree', 'trend', 'none'],
                   default='offset',
                   help="use of the --adjust controls: offset = subtract their median "
                        "BrainAGE (default), adapt = adapt the normative models to the "
                        "test site, agefree = adapt them without the controls' ages, then "
                        "offset, trend = linear trend of BrainAGE on age, none.  z-maps "
                        "are adapted to the test site unless none")
    p.add_argument('--test-male', help="text file with the sex of the test subjects "
                                       "(1 = male), if the mat-files contain none")
    p.add_argument('--normative-only', action='store_true',
                   help="only the voxel/vertex-wise normative models (no PCA, no brain age, "
                        "no GPR): save them (--save-model) and/or z-maps of --test")
    p.add_argument('--zmaps', action='store_true',
                   help="also fit the voxel/vertex-wise normative models and save z-maps of "
                        "the test subjects (per voxel/vertex, independent of --pca)")
    p.add_argument('--warp', action='store_true',
                   help="warp every voxel/vertex for the brain age models (sinh-arcsinh, "
                        "fitted jointly with a voxel-wise normative model)")
    p.add_argument('--warp-zmaps', action='store_true',
                   help="warp the data of the voxel/vertex-wise normative models (recommended: "
                        "calibrated tails of the z-maps); independent of --warp, but shares "
                        "its fit")
    p.add_argument('--train-cov', help="covariate table of the training subjects (e.g. IQMs): "
                                       "one row per subject in the order of the mat-files, "
                                       "used for all models; header with column names "
                                       "optional; '+' joins the tables of '+'-joined samples")
    p.add_argument('--test-cov', help="covariate table of the test subjects (same columns)")
    p.add_argument('--cov-mean', help="covariates in the mean: comma-separated column names "
                                      "or 1-based numbers (default: all numeric columns; "
                                      "none = none)")
    p.add_argument('--cov-sd', help="covariates in the log SD (as --cov-mean; default: none)")
    p.add_argument('--cov-df', type=int, default=1,
                   help="spline df of the covariates (default 1 = linear)")
    p.add_argument('--kfold', type=int,
                   help="folds of the cross-validation of the brain age if no --test "
                        "(default 10; with --save-model 0 = no cross-validation)")
    p.add_argument('--parcellation', action='store_true',
                   help="lobe-wise brain age and lobe-wise mean z of the z-maps")
    p.add_argument('--pca', type=int, default=100,
                   help="NDM brain age only: normative models for this many PCA scores per "
                        "model/region (default 100); 0 = every voxel/vertex.  The voxel/"
                        "vertex-wise normative model of the z-maps never uses PCA")
    p.add_argument('--rank', type=int, default=20,
                   help="NDM brain age with --pca 0: rank of the residual correlation "
                        "(default 20)")
    p.add_argument('--psi-min', type=float, default=0.01,
                   help="NDM brain age: floor of the unique variance of the residual correlation")
    p.add_argument('--df-mu', type=int, default=5, help="spline df of age for mu")
    p.add_argument('--df-sigma', type=int, default=3, help="spline df of age for sigma")
    p.add_argument('--grid-step', type=float, default=0.25,
                   help="NDM brain age: age grid step [years]")
    p.add_argument('--grid-margin', type=float, default=5.0,
                   help="NDM brain age: extend age grid beyond training range [years]")
    p.add_argument('--age-range', nargs=2, type=float, default=[0, np.inf])
    p.add_argument('--ensemble', choices=['gls', 'mae', 'mean'], default='gls')
    p.add_argument('--no-gpr', action='store_true', help="skip the GPR baseline")
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--out', default='BA_ndm_results', help="output prefix")
    a = p.parse_args(argv)

    if bool(a.train) == bool(a.model):
        p.error("give either --train or --model.")
    if a.model and not a.test:
        p.error("--model needs --test.")
    if a.model and (a.save_model or a.train_cov or a.cov_mean or a.cov_sd):
        p.error("--save-model, --train-cov, --cov-mean and --cov-sd need --train "
                "(saved models only need --test-cov).")
    if (a.cov_mean or a.cov_sd) and not a.train_cov:
        p.error("--cov-mean and --cov-sd need --train-cov.")
    if a.save_model and os.path.splitext(a.save_model)[1].lower() not in ('.mat', '.npz'):
        p.error("--save-model needs the extension .mat or .npz.")
    kfold = a.kfold if a.kfold is not None else (0 if a.save_model else 10)
    brainage, voxel = not a.normative_only, a.normative_only or a.zmaps
    if a.normative_only:
        if not (a.test or a.save_model):
            p.error("--normative-only needs --test and/or --save-model.")
        if a.correction == 'agefree':
            p.error("--correction agefree needs brain age models (not with --normative-only).")
    elif not (a.test or a.save_model or kfold):
        p.error("nothing to do: --kfold 0 without --test or --save-model.")
    if a.warp_zmaps and not voxel:
        print("--warp-zmaps has no effect without --zmaps or --normative-only.")
    t0 = time.time()

    test = None
    if a.test:
        test = [load_data(s) for s in a.test]
        _check_same_subjects(test)
        if a.test_male:
            male = np.loadtxt(a.test_male, dtype=np.float64).ravel()
            if male.size != test[0].n:
                p.error(f"--test-male has {male.size} values for {test[0].n} test subjects.")
            test = [replace(d, male=male, has_male=True) for d in test]

    if a.train:
        if test is not None and len(test) != len(a.train):
            p.error("--test needs one file per --train file.")
        train = [load_data(s) for s in a.train]
        for d in train:
            print(f"{d.name}: {d.n} subjects, {d.Y.shape[1]} features"
                  f"{' (surface)' if d.is_surf else ''}")
        if test is not None and not all(d.has_male for d in test):
            print("Test data contain no sex: sex is not used in the normative models.")
            train = [replace(d, male=np.zeros_like(d.male)) for d in train]
            test = [replace(d, male=np.zeros_like(d.male)) for d in test]
        cov = None
        if a.train_cov:
            _check_same_subjects(train)
            names, values, numeric = read_table(a.train_cov)
            cov = covariates_from_table(names, values, numeric, a.cov_mean, a.cov_sd, a.cov_df)
            C = covariate_matrix(names, values, cov, train[0].n, a.train_cov)
            train = [replace(d, C=C, cov_names=cov.names) for d in train]
            print(f"Covariates ({'linear' if cov.df == 1 else f'natural splines, df {cov.df}'}"
                  f") in the mean: {', '.join(cov.mean) or '-'}; in the log SD: "
                  f"{', '.join(cov.sd) or '-'}")
        kw = dict(df_mu=a.df_mu, df_sigma=a.df_sigma, pca=a.pca, rank=a.rank,
                  psi_min=a.psi_min, grid_step=a.grid_step, grid_margin=a.grid_margin,
                  parcellation=a.parcellation, warp=a.warp, cov=cov)
        if test is not None or a.save_model:
            ok = _valid_train(train[0], a.age_range)
            if not np.all(ok):
                print(f"{int(np.sum(~ok))} training subject(s) excluded (invalid age, outside "
                      "age range or missing covariates).")
            fit_train = [d.subset(np.flatnonzero(ok)) for d in train]
            contents = ('voxel/vertex-wise normative models' if not brainage else
                        'brain age and voxel/vertex-wise normative models' if voxel else
                        'brain age models')
            print(f"Fitting {contents} to {fit_train[0].n} subjects", flush=True)
            sets = fit_models(fit_train, kw, brainage, voxel, a.warp_zmaps)
            desc = model_description(fit_train, kw, contents, a.warp_zmaps, a.age_range, sets)
            if a.save_model:
                save_models(a.save_model, sets, desc)
        if test is None:
            if brainage and kfold:
                if a.zmaps:
                    print("z-maps are only computed for --test.")
                out = cross_validate(train, kfold, a.seed, a.age_range, not a.no_gpr,
                                     a.ensemble, **kw)
                save_results(out, a.out)
            print(f"Total time {time.time() - t0:.0f}s")
            return
    else:
        sets, desc = load_models(a.model)
        print(_describe(desc, a.model))
        cov = Covariates.from_state(desc['covariates']) if desc.get('covariates') else None
        if len(sets) != len(test):
            p.error(f"--test needs one file per saved model ({len(sets)}).")
        if brainage and any(s.brainage is None for s in sets):
            p.error(f"{a.model} contains no brain age models; use --normative-only.")
        if voxel and any(s.voxel is None for s in sets):
            p.error(f"{a.model} contains no voxel/vertex-wise normative models (they are "
                    "saved with --zmaps or --normative-only).")
        if not all(d.has_male for d in test) and any(_uses_sex(s) for s in sets):
            p.error("The saved models use sex, but the test data contain none: give --test-male.")
        if brainage and not a.no_gpr:
            print("GPR baseline skipped: it needs the training data.")
        fit_train = None

    used = set().union(*_used_covariates(sets)) if cov is not None else set()
    if used:
        if not a.test_cov:
            p.error(f"The models use the covariates {', '.join(sorted(used))}: give --test-cov.")
        names, values, _ = read_table(a.test_cov)
        C = covariate_matrix(names, values, cov, test[0].n, a.test_cov,
                             [nm for nm in cov.names if nm in used])
        _check_covariates(C, desc, used)
        test = [replace(d, C=C, cov_names=cov.names) for d in test]
    elif a.test_cov:
        print("--test-cov is ignored: the models use no covariates.")
    _check_ages(test[0].age, desc)
    adjust = _parse_index(a.adjust, test[0].n) if a.adjust else None
    out = apply_models(sets, test, adjust, a.correction, a.ensemble, a.zmaps,
                       a.normative_only, a.parcellation, fit_train, not a.no_gpr)
    save_results(out, a.out)
    print(f"Total time {time.time() - t0:.0f}s")


if __name__ == '__main__':
    main()
