#!/usr/bin/env python3
"""
BA_ndm.py - Brain age from normative models (normative deviation mapping, NDM).

Prototype of a likelihood-based brain age estimator as an alternative to the
GPR-based BrainAGE (BA_gpr_ui.m).  Every feature gets its own normative model

    y_j ~ N(mu_j(age, sex, site), sigma_j(age, sex)^2)

with natural cubic splines of age for mu and log(sigma), fitted by maximum
likelihood (the location-scale RS algorithm of ComBatLS, ported from
combat_family.py in the ComCat repository).  Brain age is the age at which the
subject's data are most likely,

    age_hat = argmax_a  log N(z(a); 0, R) - sum_j log sigma_j(a),

where z(a) are the z-scores at candidate age a and R is the correlation
between features, estimated from the training z-scores.

By default (--pca 100) the features are the leading principal component
scores of the training data, so that R can be estimated in full.  With --pca 0
every voxel/vertex is a feature and R is a low-rank plus diagonal model; the
residual correlation of smoothed voxels is spread over hundreds of dimensions,
so this model is overconfident (standard errors several times too small) and
less accurate.  Treating features as independent (--pca 0 --rank 0)
overcounts the evidence of correlated features and is worse still.

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
chronological age (NDM.Deviation_age).  With --zmaps, voxel/vertex-wise z-maps
of the test subjects are saved as well.

For comparison, a Python replica of the GPR BrainAGE (BA_gpr.m, linear kernel,
PCA, linear trend correction with trend_method = 1) is run on the same folds.

Input files are the mat-files written by BA_data2mat.m (Y, age, male and, for
surface data, ind; MATLAB v5 or v7.3).  Files joined with '+' are concatenated
and their position serves as site.  Subjects with non-finite age or age <= 0
are excluded from training and evaluation.

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

Results are saved as <out>.mat (readable by MATLAB) and <out>.csv, and with
--zmaps as <out>_zmaps_<model>.mat.

Requirements: numpy, scipy, h5py (v7.3 files), nibabel (surface parcellation)
"""

from __future__ import annotations

import argparse
import os
import re
import time
from dataclasses import dataclass, replace

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))

# regions of the lobe atlas that are not represented on the surface
_SURF_EXCLUDE = (5, 15)


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

    @property
    def n(self):
        return self.Y.shape[0]

    def subset(self, idx):
        return Data(self.Y[idx], self.age[idx], self.male[idx], self.site[idx],
                    self.name, self.ind, self.is_surf, self.res, self.has_male)


def _read_mat(path):
    """Read Y, age, male and ind of a BA_data2mat mat-file (v5 or v7.3)."""
    keys = ('Y', 'age', 'male', 'ind')
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
    Ys, ages, males, sites, ind = [], [], [], [], None
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
        Ys.append(np.asarray(Y, dtype=np.float32))
        ages.append(age)
        males.append(male)
        sites.append(np.full(age.size, s))
    base = os.path.basename(parts[0])
    res = re.search(r'_(\d+)mm_', base)
    return Data(np.vstack(Ys), np.concatenate(ages), np.concatenate(males),
                np.concatenate(sites), '+'.join(os.path.basename(p) for p in parts),
                ind, 'mesh' in base, res.group(1) if res else None, has_male)


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
    """Natural cubic spline basis of age without intercept (ESL eq. 5.4-5.5).

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


def _deviance(r, eta):
    """-2 log-likelihood of N(mu, exp(eta)^2) given residuals r = y - mu."""
    return np.sum(np.log(2 * np.pi) + 2 * eta + r ** 2 * np.exp(-2 * eta), axis=0)


def fit_location_scale(Y, X, W, max_iter=2000, tol=1e-6, max_elements=4_000_000):
    """ML fit of y_j ~ N(X beta_j, exp(W theta_j)^2) for every column j of Y.

    Port of _fit_location_scale() in combat_family.py (ComCat) for data of
    shape (n_subjects, n_features).  RS algorithm as in gamlss (family NO,
    identity link for mu, log link for sigma): alternate a weighted
    least-squares update of beta with a Fisher-scoring update of theta,
    halving the theta step if the deviance increases.  W[:, 0] must be the
    intercept.  Features are processed in chunks to bound memory.

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

        # start: OLS mean, constant sigma
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


class NormativeModel:
    """Location-scale normative model for every feature.

    mu    = 1 + ns(age, df_mu) + male + site   (site: fixed effects)
    log sigma = 1 + ns(age, df_sigma) + male

    Predictions use the size-weighted mean of the training sites (as ComBat's
    stand_mean), optionally adapted to a new site by adapt().
    """

    def __init__(self, df_mu=5, df_sigma=3, max_iter=2000, tol=1e-6):
        self.df_mu, self.df_sigma = df_mu, df_sigma
        self.max_iter, self.tol = max_iter, tol

    def _design(self, age, male, site=None):
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
        return np.hstack(X), np.hstack(W)

    def fit(self, Y, age, male, site=None, verbose=False):
        n = len(age)
        site = np.zeros(n, dtype=int) if site is None else np.asarray(site)
        self.site_levels, counts = np.unique(site, return_counts=True)
        self.site_props = counts / n
        self.use_male = np.unique(male).size > 1
        self.basis_mu = NaturalSpline(age, self.df_mu)
        self.basis_sigma = NaturalSpline(age, self.df_sigma)

        finite = np.all(np.isfinite(Y), axis=0)
        sd = np.zeros(Y.shape[1])
        sd[finite] = np.std(Y[:, finite], axis=0)
        self.valid = np.flatnonzero(finite & (sd > 0))
        self.n_features = Y.shape[1]

        X, W = self._design(age, male, site)
        t0 = time.time()
        self.beta, self.theta, conv = fit_location_scale(
            Y[:, self.valid], X, W, self.max_iter, self.tol)
        if verbose:
            print(f"    normative model: {self.valid.size} features, "
                  f"{int(np.sum(~conv))} not converged ({time.time() - t0:.1f}s)")
        self.offset = np.zeros(self.valid.size)
        self.scale = np.ones(self.valid.size)
        return self

    def mu_sigma(self, age, male, site=None, cols=None):
        """mu and sigma (n, p_valid) at the given ages; cols selects valid features."""
        cols = slice(None) if cols is None else cols
        X, W = self._design(np.atleast_1d(age), np.atleast_1d(male), site)
        sd = np.exp(W @ self.theta[:, cols])
        mu = X @ self.beta[:, cols] + self.offset[cols] * sd
        return mu, sd * self.scale[cols]

    def zscores(self, Y, age, male, site=None):
        """Normative deviation maps (n, p_valid) at the given ages."""
        mu, sd = self.mu_sigma(age, male, site)
        return (Y[:, self.valid] - mu) / sd

    def adapt(self, Y, age, male):
        """Adapt to a new site using its control subjects.

        Removes the mean and scales the SD of the controls' z-scores to 0 and 1
        (the ComBat location/scale step without empirical Bayes).
        """
        self.offset[:] = 0
        self.scale[:] = 1
        z = self.zscores(Y, age, male)
        self.offset = z.mean(axis=0)
        self.scale = z.std(axis=0, ddof=1)
        self.scale[~(self.scale > 0)] = 1


def _warp_loglik(x, mu, sd, eps, delta, prior):
    """Log-likelihood per column of x warped by w = sinh(delta asinh(x) - eps):
    sum_i log N(w_i; mu_i, sd_i) + log dw/dx_i (without constants), plus the prior."""
    t = delta * np.arcsinh(x) - eps
    ll = -0.5 * ((np.sinh(t) - mu) / sd) ** 2 + np.log(np.cosh(t))
    return ll.sum(axis=0) + x.shape[0] * np.log(delta) - 0.5 * (eps ** 2 + np.log(delta) ** 2) / prior


def fit_warp(x, mu, sd, eps, delta, n_iter=20, prior=9.0, tol=1e-6):
    """Update the sinh-arcsinh warp w = sinh(delta asinh(x) - eps) of every column of x.

    x (n, p) standardized data; mu, sd (n, p) location and scale of the warped
    data.  Maximizes sum_i [log N(w_i; mu_i, sd_i) + log dw/dx_i] with a weak
    normal prior (variance prior) on eps and log(delta) centred on the identity
    warp.  Damped Newton steps in (eps, log delta) with step halving, so the
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
        # parameters (eps, eta = log delta)
        g_h = de * g_d - eta[active] / prior
        h_eh = de * h_ed
        h_hh = de ** 2 * h_dd + de * g_d - 1 / prior
        # make the Hessian negative definite (Levenberg damping), then Newton step
        lam_max = 0.5 * (h_ee + h_hh) + np.sqrt(0.25 * (h_ee - h_hh) ** 2 + h_eh ** 2)
        damp = np.maximum(0, lam_max + 1e-6 * n)
        a, b, d = h_ee - damp, h_eh, h_hh - damp
        det = a * d - b * b
        step_e = -(d * g_e - b * g_h) / det
        step_h = -(-b * g_e + a * g_h) / det
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
    coordinate ascent on the joint likelihood: the normative model is fitted to
    the warped data, then (eps, delta) are updated given its location and scale
    (fit_warp), n_outer times.  The warp does not depend on age, so its Jacobian
    is constant in the brain age likelihood; the warped data simply replace the
    raw data in all models.
    """

    def __init__(self, df_mu=5, df_sigma=3, n_outer=4, prior=9.0, chunk=4000):
        self.df_mu, self.df_sigma = df_mu, df_sigma
        self.n_outer, self.prior, self.chunk = n_outer, prior, chunk

    def fit(self, Y, age, male, site=None, verbose=False):
        finite = np.all(np.isfinite(Y), axis=0)
        sd = np.zeros(Y.shape[1])
        sd[finite] = np.std(Y[:, finite], axis=0)
        self.valid = np.flatnonzero(finite & (sd > 0))
        self.center = np.mean(Y[:, self.valid], axis=0, dtype=np.float64)
        self.sd = sd[self.valid]
        self.eps = np.zeros(self.valid.size)
        self.delta = np.ones(self.valid.size)
        self.loglik = []
        for _ in range(self.n_outer):
            model = NormativeModel(self.df_mu, self.df_sigma).fit(
                self.transform(Y), age, male, site)
            if not np.array_equal(model.valid, self.valid):
                raise ValueError("Warped data have constant or non-finite features.")
            X, W = model._design(age, male, site)
            total = 0.0
            for c0 in range(0, self.valid.size, self.chunk):
                cols = slice(c0, c0 + self.chunk)
                x = self._standardize(Y, cols)
                mu, sd = X @ model.beta[:, cols], np.exp(W @ model.theta[:, cols])
                self.eps[cols], self.delta[cols] = fit_warp(
                    x, mu, sd, self.eps[cols], self.delta[cols], prior=self.prior)
                total += np.sum(_warp_loglik(x, mu, sd, self.eps[cols], self.delta[cols], self.prior)
                                - np.sum(np.log(sd), axis=0))
            self.loglik.append(total)
            if verbose:
                print(f"    warp: log-likelihood {total:.6g}")
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


def _loglik_grid(model, corr, cols, Y, male_val, grid, max_bytes=2e8):
    """log p(y | a) (n, G) on the age grid for subjects of one sex.

    Expands sum_j D_j z_j^2 and the Woodbury term into matrix products so that
    all subjects and a chunk of grid ages are handled by one GEMM.
    """
    n, p = Y.shape
    k = corr.A.shape[1]
    mref = model.mu_sigma(np.median(grid), male_val, cols=cols)[0][0]
    Yt = np.asarray(Y, dtype=np.float64) - mref    # centring limits cancellation
    Y2 = Yt ** 2
    gc = int(max(1, min(len(grid), max_bytes // (8 * p * max(k, 1)))))
    ll = np.empty((n, len(grid)))
    for g0 in range(0, len(grid), gc):
        ag = grid[g0:g0 + gc]
        mu, sd = model.mu_sigma(ag, np.full(len(ag), male_val), cols=cols)
        mu -= mref
        s = 1 / sd
        Dq = corr.D * s ** 2
        sumsq = Y2 @ Dq.T - 2 * (Yt @ (Dq * mu).T) + np.sum(Dq * mu ** 2, axis=1)
        part = -0.5 * sumsq - np.sum(np.log(sd), axis=1)
        if k:
            B = corr.A[None] * s[:, :, None]      # (gc, p, k)
            c = np.einsum('gp,gpk->gk', mu, B)
            P = (Yt @ B.transpose(1, 0, 2).reshape(p, -1)).reshape(n, len(ag), k)
            part += 0.5 * np.sum((P - c) ** 2, axis=2)
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


def _mahalanobis(model, corr, cols, Y, age, male):
    """Squared Mahalanobis distance z' R^-1 z of the z-scores at the given ages (n,)."""
    mu, sd = model.mu_sigma(age, male, cols=cols)
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


class NDMBrainAge:
    """Normative model + residual correlation -> likelihood-based brain age.

    pca > 0 : normative models for the leading pca principal component scores
              of the whole brain (and of each region) with a full residual
              correlation between the scores.
    pca = 0 : normative model for every voxel/vertex with a low-rank (rank)
              plus diagonal residual correlation.  Its correlation model misses
              most of the dependence between voxels, so the standard errors
              are much too small and the estimates less accurate.
    warp    : warp every voxel/vertex with a sinh-arcsinh function fitted jointly
              with a voxel-wise normative model (Warp) before all other steps
    """

    def __init__(self, df_mu=5, df_sigma=3, pca=100, rank=20, psi_min=0.01,
                 grid_step=0.25, grid_margin=5.0, parcellation=False,
                 atlas_dir=HERE, warp=False, verbose=False):
        self.df_mu, self.df_sigma = df_mu, df_sigma
        self.pca, self.rank, self.psi_min = pca, rank, psi_min
        self.grid_step, self.grid_margin = grid_step, grid_margin
        self.parcellation, self.atlas_dir = parcellation, atlas_dir
        self.warp, self.verbose = warp, verbose
        self.warper = None

    def _prep(self, d: Data):
        """Data with the training warp applied (unchanged without warp)."""
        return d if self.warper is None else replace(d, Y=self.warper.transform(d.Y))

    def _normative(self, Y, d):
        return NormativeModel(self.df_mu, self.df_sigma).fit(
            Y, d.age, d.male, d.site, self.verbose)

    def _fit_pca_part(self, d, name, cols):
        X = d.Y[:, cols]
        keep = np.all(np.isfinite(X), axis=0) & (np.std(X, axis=0) > 0)
        cols = cols[keep]
        center, Vt = _pca(X[:, keep], self.pca)
        part = _Part(name, cols, None, None, None, center, Vt)
        F = part.model_input(d.Y)
        part.model = self._normative(F, d)
        part.mcols = np.arange(part.model.valid.size)
        Z = part.model.zscores(F, d.age, d.male, d.site)
        part.corr = ResidualCorrelation(Z.shape[1], self.psi_min).fit(Z)
        return part

    def fit(self, d: Data):
        if self.warp:
            self.warper = Warp(self.df_mu, self.df_sigma).fit(d.Y, d.age, d.male, d.site, self.verbose)
            d = self._prep(d)
        groups = [('global', None, np.arange(d.Y.shape[1]))]
        if self.parcellation:
            atlas, regions, names = lobe_atlas(d, self.atlas_dir)
            groups += [(nm, r, np.flatnonzero(atlas == r)) for r, nm in zip(regions, names)]

        self.parts, self.regions, self.region_names = [], [], []
        if not self.pca:
            model = self._normative(d.Y, d)
            Z = model.zscores(d.Y, d.age, d.male, d.site)
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

    def adapt(self, d: Data):
        """Adapt the normative models to a new site using its controls d."""
        d = self._prep(d)
        for part in self._models():
            part.model.adapt(part.model_input(d.Y), d.age, d.male)
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
                age = self._estimate(part, X, d.male)[0]
                part.model.adapt(F, age, d.male)
                if prev is not None and np.max(np.abs(age - prev)) < tol:
                    break
                prev = age
        return self

    def _estimate(self, part, X, male, chunk=256):
        """Maximum likelihood age, Laplace SD and boundary flag for the features X of a part."""
        n = len(X)
        age, sd, bound = np.full(n, np.nan), np.full(n, np.nan), np.zeros(n, dtype=bool)
        for male_val in np.unique(male):
            idx = np.flatnonzero(male == male_val)
            for c0 in range(0, idx.size, chunk):
                sub = idx[c0:c0 + chunk]
                ll = _loglik_grid(part.model, part.corr, part.mcols, X[sub], male_val, self.grid)
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
        for q, part in enumerate(self.parts):
            X = part.features(d.Y)
            if not np.all(np.isfinite(X)):
                raise ValueError(f"{d.name}: non-finite values in test data.")
            age, res['sd'][:, q], res['at_bound'][:, q] = self._estimate(part, X, d.male, chunk)
            res['age'][:, q] = age
            k = part.mcols.size
            for c0 in range(0, n, chunk):
                sub = np.arange(c0, min(c0 + chunk, n))
                res['deviation'][sub, q] = deviation_z(_mahalanobis(
                    part.model, part.corr, part.mcols, X[sub], age[sub], d.male[sub]), k)
            for c0 in range(0, ok.size, chunk):
                sub = ok[c0:c0 + chunk]
                res['deviation_age'][sub, q] = deviation_z(_mahalanobis(
                    part.model, part.corr, part.mcols, X[sub], d.age[sub], d.male[sub]), k)
        out = {key: v[:, 0] for key, v in res.items()}
        out['regional'] = res['age'][:, 1:]
        out.update({f'regional_{key}': res[key][:, 1:] for key in ('sd', 'deviation', 'deviation_age')})
        return out


def zmaps(dtr: Data, dte: Data, ctrl=None, ctrl_age=None, df_mu=5, df_sigma=3,
          parcellation=False, atlas_dir=HERE, chunk=256, warper=None):
    """Voxel/vertex-wise normative z-maps of the test subjects at chronological age.

    A location-scale model is fitted to every voxel/vertex of the training data
    (after the sinh-arcsinh warp of a fitted Warp, if given).  Given the indices
    ctrl of control subjects of the test data (and their ages ctrl_age, by
    default their chronological ages), the model is first adapted to the test
    site.  Returns Z (n, n_features) as float32, with NaN for features without a
    model and for subjects without a valid age, and, with parcellation, the
    lobe-wise mean z (n, n_regions) and the region labels.
    """
    if warper is not None:
        dtr, dte = replace(dtr, Y=warper.transform(dtr.Y)), replace(dte, Y=warper.transform(dte.Y))
    model = NormativeModel(df_mu, df_sigma).fit(dtr.Y, dtr.age, dtr.male, dtr.site)
    if ctrl is not None:
        model.adapt(dte.Y[ctrl], dte.age[ctrl] if ctrl_age is None else ctrl_age, dte.male[ctrl])
    Z = np.full(dte.Y.shape, np.nan, dtype=np.float32)
    ok = np.flatnonzero(np.isfinite(dte.age) & (dte.age > 0))
    for c0 in range(0, ok.size, chunk):
        sub = ok[c0:c0 + chunk]
        Z[np.ix_(sub, model.valid)] = model.zscores(dte.Y[sub], dte.age[sub], dte.male[sub])
    out = dict(Z=Z)
    if parcellation:
        atlas, regions, names = lobe_atlas(dte, atlas_dir)
        out['regional_z'] = np.column_stack([np.nanmean(Z[:, atlas == r], axis=1) for r in regions])
        out['regions'] = np.array(regions)
    return out


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

def _valid_age(d: Data, age_range):
    return np.isfinite(d.age) & (d.age > 0) & (d.age >= age_range[0]) & (d.age <= age_range[1])


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
    ok = _valid_age(d0, age_range)
    if np.sum(~ok):
        print(f"{int(np.sum(~ok))} subject(s) excluded (invalid age or outside age range).")
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


def train_test(train, test, adjust=None, correction='offset', age_range=(0, np.inf),
               gpr=True, ensemble='gls', zmaps_out=False, **kw):
    """Train on one sample and predict another.

    adjust     : indices of the control subjects of the test sample
    correction : use of the controls for the test site: 'offset' subtracts their
                 median BrainAGE, 'adapt' adapts the normative models (location
                 and scale of the z-scores), 'agefree' adapts them without the
                 controls' ages (at their estimated brain ages) and then subtracts
                 their median BrainAGE, 'trend' removes a linear trend of BrainAGE
                 on age, 'none' only estimates the ensemble weights.  The GPR
                 baseline is corrected in the same way ('offset' for 'adapt' and
                 'agefree').
    zmaps_out  : also return voxel/vertex-wise z-maps at chronological age, adapted
                 to the test site with the controls (at their estimated brain ages
                 for 'agefree')
    """
    _check_same_subjects(test)
    n = test[0].n
    age = test[0].age
    ok = np.isfinite(age) & (age > 0)
    if adjust is None:
        correction = 'none'
        ctrl = np.flatnonzero(ok)
    else:
        ctrl = np.intersect1d(adjust, np.flatnonzero(ok))
    if not all(d.has_male for d in test):
        print("Test data contain no sex: sex is not used in the normative models.")
        train = [replace(d, male=np.zeros_like(d.male)) for d in train]
        test = [replace(d, male=np.zeros_like(d.male)) for d in test]
    post = correction if correction in ('offset', 'trend') else 'none'
    if correction == 'agefree':
        post = 'offset'
    gpr_post = 'offset' if correction in ('adapt', 'agefree') else post
    res = []
    for dtr, dte in zip(train, test):
        sel = np.flatnonzero(_valid_age(dtr, age_range))
        t0 = time.time()
        est = NDMBrainAge(**kw).fit(dtr.subset(sel))
        if correction == 'adapt':
            est.adapt(dte.subset(adjust))
        elif correction == 'agefree':
            est.adapt_agefree(dte.subset(adjust))
        out = est.predict(dte)
        r = dict(name=dte.name, age=correct_age(out['age'], age, ctrl, post),
                 sd=out['sd'], at_bound=out['at_bound'],
                 regional=correct_age(out['regional'], age, ctrl, post), regions=est.regions,
                 region_names=est.region_names, gpr=np.full(n, np.nan),
                 offset=np.nanmedian(out['age'][ctrl] - age[ctrl]))
        for key in ('deviation', 'deviation_age', 'regional_deviation', 'regional_deviation_age'):
            r[key] = out[key]
        if gpr:
            r['gpr'] = gpr_baseline(dtr.Y[sel], dtr.age[sel], dte.Y)
            r['gpr_ba'] = correct_age(r['gpr'], age, ctrl, gpr_post) - age
        if zmaps_out:
            zc = None if correction == 'none' else adjust
            r['zmaps'] = zmaps(dtr.subset(sel), dte, zc,
                               out['age'][adjust] if correction == 'agefree' else None,
                               est.df_mu, est.df_sigma, est.parcellation, est.atlas_dir,
                               warper=est.warper)
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
    if zmaps_out:
        out['zmaps'] = [dict(r['zmaps'], model=r['name'], ind=d.ind, age=age)
                        for r, d in zip(res, test)]
    return out


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
    """Save <prefix>.mat (struct NDM), <prefix>.csv and, if computed, one
    <prefix>_zmaps_<model>.mat per model (struct NDMzmap: Z (n, n_features) in
    the feature order of the input Y, regional_z, regions, age, model, ind)."""
    from scipy.io import savemat
    out = dict(out)
    zm = out.pop('zmaps', None)
    savemat(prefix + '.mat', {'NDM': out}, do_compression=True)
    names = [os.path.splitext(m)[0] for m in out['models']]
    for z, nm in zip(zm or [], names):
        z = {k: (np.zeros(0) if v is None else v) for k, v in z.items()}
        savemat(f'{prefix}_zmaps_{nm}.mat', {'NDMzmap': z}, do_compression=True)
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
    print(f"\nSaved {prefix}.mat and {prefix}.csv"
          + (f" and {len(zm)} z-map file(s) {prefix}_zmaps_*.mat" if zm else ''))


def main(argv=None):
    p = argparse.ArgumentParser(
        description="Brain age from normative models (NDM prototype).",
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    p.add_argument('--train', nargs='+', required=True,
                   help="training mat-file per model ('+' joins sites)")
    p.add_argument('--test', nargs='+', help="test mat-file per model (same order as --train)")
    p.add_argument('--adjust', help="1-based indices of test controls, e.g. 1:108")
    p.add_argument('--correction', choices=['offset', 'adapt', 'agefree', 'trend', 'none'],
                   default='offset',
                   help="use of the --adjust controls: offset = subtract their median "
                        "BrainAGE (default), adapt = adapt the normative models to the "
                        "test site, agefree = adapt them without the controls' ages, then "
                        "offset, trend = linear trend of BrainAGE on age, none")
    p.add_argument('--test-male', help="text file with the sex of the test subjects "
                                       "(1 = male), if the mat-files contain none")
    p.add_argument('--zmaps', action='store_true',
                   help="also save voxel/vertex-wise z-maps of the test subjects")
    p.add_argument('--warp', action='store_true',
                   help="warp every voxel/vertex (sinh-arcsinh, fitted jointly with a "
                        "voxel-wise normative model) for non-Gaussian data")
    p.add_argument('--kfold', type=int, default=10, help="folds if no --test (default 10)")
    p.add_argument('--parcellation', action='store_true', help="lobe-wise brain age")
    p.add_argument('--pca', type=int, default=100,
                   help="normative models for this many PCA scores per model/region "
                        "(default 100); 0 = every voxel/vertex")
    p.add_argument('--rank', type=int, default=20,
                   help="rank of residual correlation for --pca 0 (default 20)")
    p.add_argument('--psi-min', type=float, default=0.01, help="floor of unique variance")
    p.add_argument('--df-mu', type=int, default=5, help="spline df of age for mu")
    p.add_argument('--df-sigma', type=int, default=3, help="spline df of age for sigma")
    p.add_argument('--grid-step', type=float, default=0.25, help="age grid step [years]")
    p.add_argument('--grid-margin', type=float, default=5.0,
                   help="extend age grid beyond training range [years]")
    p.add_argument('--age-range', nargs=2, type=float, default=[0, np.inf])
    p.add_argument('--ensemble', choices=['gls', 'mae', 'mean'], default='gls')
    p.add_argument('--no-gpr', action='store_true', help="skip the GPR baseline")
    p.add_argument('--seed', type=int, default=0)
    p.add_argument('--out', default='BA_ndm_results', help="output prefix")
    a = p.parse_args(argv)

    kw = dict(df_mu=a.df_mu, df_sigma=a.df_sigma, pca=a.pca, rank=a.rank, psi_min=a.psi_min,
              grid_step=a.grid_step, grid_margin=a.grid_margin,
              parcellation=a.parcellation, warp=a.warp)
    t0 = time.time()
    train = [load_data(s) for s in a.train]
    for d in train:
        print(f"{d.name}: {d.n} subjects, {d.Y.shape[1]} features"
              f"{' (surface)' if d.is_surf else ''}")
    if a.test:
        if len(a.test) != len(a.train):
            p.error("--test needs one file per --train file.")
        test = [load_data(s) for s in a.test]
        if a.test_male:
            male = np.loadtxt(a.test_male, dtype=np.float64).ravel()
            if male.size != test[0].n:
                p.error(f"--test-male has {male.size} values for {test[0].n} test subjects.")
            test = [replace(d, male=male, has_male=True) for d in test]
        adjust = _parse_index(a.adjust, test[0].n) if a.adjust else None
        out = train_test(train, test, adjust, a.correction, a.age_range, not a.no_gpr,
                         a.ensemble, a.zmaps, **kw)
    else:
        if a.zmaps:
            print("--zmaps is only available with --test.")
        out = cross_validate(train, a.kfold, a.seed, a.age_range, not a.no_gpr,
                             a.ensemble, **kw)
    save_results(out, a.out)
    print(f"Total time {time.time() - t0:.0f}s")


if __name__ == '__main__':
    main()
