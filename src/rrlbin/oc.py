"""Level-2 O-C analysis of season-delay series, shared by the real-data pipeline, the simulations and re-analyses.

A series is a dict of arrays over seasons:
  t [d], tau [d], err [d], band (0 = OGLE I, 1 = MACHO B, 2 = MACHO R),
  alpha, alpha_err (amplitude factors), coh, coh_err (harmonic coherence, fundamental minus higher harmonics [d])
and P (pulsation period [d]).

Model for the delays (days), with Gaussian errors of variance err^2 + s^2 (white jitter s profiled):
  H0: c0 + c1 x + c2 x^2 + sum_b Delta_b [band == b]            (x = (t - t_ref)/1000 d; b = MACHO bands)
  H1: H0 + sum_{h=1}^{n_harm} [A_h sin(2 pi h t / P_orb) + B_h cos(2 pi h t / P_orb)]
Optional Gaussian priors on the band offsets Delta_b (from the measured band-lag relations) enter as pseudo-observations
that do not receive the jitter term. D = 2 max_P [lnL1(P) - lnL0].
"""
from __future__ import annotations

import numpy as np

from .ltte import _gls_lnl

DAY = 86400.0
YR = 365.25
BANDS_EXTRA = (1, 2)            # bands with a free offset relative to OGLE I


# ------------------------------------------------------------------ unwrapping
def robust_unwrap(t, tau, band, P, n_iter=4, half_window=2):
    """Remove whole-cycle ambiguities within each band's series.

    First a sequential unwrap in time order; then, iteratively, each season is moved by the integer number of cycles
    that brings it closest to a local linear prediction from its neighbours (up to `half_window` seasons on each side,
    itself excluded). This repairs single-season cycle slips that the sequential unwrap propagates to all later
    seasons. Returns the unwrapped delays (same order as the input)."""
    t, tau, band = map(np.asarray, (t, tau, band))
    out = tau.astype(float).copy()
    for b in np.unique(band):
        idx = np.flatnonzero(band == b)
        idx = idx[np.argsort(t[idx])]
        y = out[idx]
        y = y[0] + np.unwrap(2 * np.pi * (y - y[0]) / P) * P / (2 * np.pi)
        if idx.size >= 4:
            for _ in range(n_iter):
                changed = False
                for j in range(idx.size):
                    nb = [k for k in range(max(0, j - half_window), min(idx.size, j + half_window + 1)) if k != j]
                    if len(nb) < 2:
                        continue
                    c = np.polyfit(t[idx[nb]], y[nb], 1)
                    pred = np.polyval(c, t[idx[j]])
                    k = np.round((pred - y[j]) / P)
                    if k != 0:
                        y[j] += k * P
                        changed = True
                if not changed:
                    break
        out[idx] = y
    return out


# ------------------------------------------------------------------ design and likelihood
def design(t, band, t_ref, P_orb=None, n_harm=1):
    x = (np.asarray(t) - t_ref) / 1000.0
    cols = [x ** 2, x, np.ones_like(x)]
    names = ["c2", "c1", "c0"]
    for b in BANDS_EXTRA:
        if np.any(band == b):
            cols.append((band == b).astype(float))
            names.append(f"off{b}")
    if P_orb is not None:
        for h in range(1, n_harm + 1):
            w = 2 * np.pi * h * np.asarray(t) / P_orb
            cols += [np.sin(w), np.cos(w)]
            names += [f"s{h}", f"c{h}"]
    return np.column_stack(cols), names


def _augment(X, y, names, priors):
    """Append prior pseudo-rows for offset parameters: (row, target, sd)."""
    rows, tgt, sd = [], [], []
    for b, (mu, s) in (priors or {}).items():
        key = f"off{b}"
        if key in names:
            r = np.zeros(X.shape[1])
            r[names.index(key)] = 1.0
            rows.append(r)
            tgt.append(mu)
            sd.append(s)
    return (np.array(rows).reshape(-1, X.shape[1]), np.array(tgt), np.array(sd))


def profile_lnl(X, y, err, names, s_grid, priors=None):
    """max over white jitter s of the Gaussian lnL of a WLS fit, including prior pseudo-rows (no jitter on them).
    Returns (lnL, s2, beta, cov_beta)."""
    Xp, yp, sp = _augment(X, y, names, priors)
    best = (-np.inf, 0.0, None, None)
    for s2 in s_grid ** 2:
        v = np.r_[err ** 2 + s2, sp ** 2]
        A = np.vstack([X, Xp])
        yy = np.r_[y, yp]
        w = 1 / np.sqrt(v)
        beta, *_ = np.linalg.lstsq(A * w[:, None], yy * w, rcond=None)
        r = yy - A @ beta
        lnl = -0.5 * np.sum(r ** 2 / v + np.log(2 * np.pi * v))
        if lnl > best[0]:
            cov = np.linalg.pinv((A * w[:, None]).T @ (A * w[:, None]))
            best = (lnl, s2, beta, cov)
    return best


def period_grid(baseline, p_min=400.0, p_max_factor=2.0, oversample=5):
    """Trial orbital periods uniform in frequency, 1/(p_max_factor * baseline) ... 1/p_min.
    p_min = 400 d: a season mean over ~240 d keeps |sinc(pi 240/P)| = 0.50 of the amplitude at 400 d (0.66 at 500 d);
    the irregular placement of seasons within the year (rms 40 d) largely breaks the annual alias (scripts/alias_test.py)."""
    f_min, f_max = 1 / (p_max_factor * baseline), 1 / p_min
    n = int(np.ceil((f_max - f_min) * baseline * oversample)) + 1
    return 1 / np.linspace(f_max, f_min, n)


def default_s_grid(err):
    return np.r_[0.0, np.geomspace(0.1, 30, 30) * np.median(err)]


def orbit_search(t, tau, err, band, periods, n_harm=1, priors=None, s_grid=None):
    """D(P) for H1 vs H0; returns dict(D, P_best, Dp, beta, names, s0, s1, lnl0)."""
    t_ref = np.average(t, weights=err ** -2)
    s_grid = default_s_grid(err) if s_grid is None else s_grid
    X0, n0 = design(t, band, t_ref)
    l0, s20, _, _ = profile_lnl(X0, tau, err, n0, s_grid, priors)
    Dp = np.empty(len(periods))
    best = (-np.inf, None)
    for k, P in enumerate(periods):
        X1, n1 = design(t, band, t_ref, P, n_harm)
        l1, s21, b, _ = profile_lnl(X1, tau, err, n1, s_grid, priors)
        Dp[k] = 2 * (l1 - l0)
        if l1 > best[0]:
            best = (l1, (b, n1, s21))
    k = int(np.argmax(Dp))
    b, n1, s21 = best[1]
    return dict(D=float(Dp[k]), P_best=float(periods[k]), Dp=Dp, beta=b, names=n1, s0=float(np.sqrt(s20)),
                s1=float(np.sqrt(s21)), t_ref=t_ref)


def orbit_amplitude(beta, names, P_orb, n_harm):
    """Semi-amplitude of the fitted orbital curve (half its peak-to-peak) and of the fundamental harmonic [d]."""
    tt = np.linspace(0, P_orb, 400)
    y = np.zeros_like(tt)
    for h in range(1, n_harm + 1):
        y += beta[names.index(f"s{h}")] * np.sin(2 * np.pi * h * tt / P_orb) + beta[names.index(f"c{h}")] * np.cos(2 * np.pi * h * tt / P_orb)
    return 0.5 * np.ptp(y), float(np.hypot(beta[names.index("s1")], beta[names.index("c1")]))


def alias_dD(Dp, periods, P_best):
    """D(P_best) minus the largest D within +-10% of the annual aliases 1/|1/P_best -+ 1/yr| (small = ambiguous)."""
    al = [1 / abs(1 / P_best - 1 / YR), 1 / (1 / P_best + 1 / YR)]
    m = np.zeros(periods.size, bool)
    for a in al:
        m |= np.abs(periods / a - 1) < 0.10
    return float(Dp.max() - Dp[m].max()) if m.any() else np.inf


def red_null(t, tau, err, band, priors=None, n=16):
    """ML white jitter s and random-walk strength q under H0 (grids). Returns (s [d], rms random walk over baseline [d])."""
    t_ref = np.average(t, weights=err ** -2)
    X, names = design(t, band, t_ref)
    Xp, yp, sp = _augment(X, tau, names, priors)
    e0, T = np.median(err), np.ptp(t)
    W = np.minimum.outer(t - t.min(), t - t.min())
    best = (-np.inf, 0, 0)
    for s in np.r_[0.0, np.geomspace(0.1, 30, n) * e0]:
        for q in np.r_[0.0, (np.geomspace(0.1, 30, n) * e0) ** 2 / T]:
            C = np.diag(err ** 2 + s ** 2) + q * W
            if len(yp):
                C = np.block([[C, np.zeros((C.shape[0], len(yp)))], [np.zeros((len(yp), C.shape[0])), np.diag(sp ** 2)]])
            l, _ = _gls_lnl(np.vstack([X, Xp]), np.r_[tau, yp], C)
            if l > best[0]:
                best = (l, s, q)
    return best[1], float(np.sqrt(best[2] * T))


def predictive_score(t, tau, err, band, test, P_orb, n_harm=1, priors=None):
    """lnL_pred(test | orbit at P_orb fitted on ~test) - lnL_pred(test | no orbit); parameter covariance included."""
    tr = ~test
    t_ref = np.average(t[tr], weights=err[tr] ** -2)
    s_grid = default_s_grid(err[tr])
    out = []
    for P in (P_orb, None):
        Xa, na = design(t, band, t_ref, P, n_harm)
        # columns must be estimable from the training set
        l, s2, b, cov = profile_lnl(Xa[tr], tau[tr], err[tr], na, s_grid, priors)
        B = Xa[test]
        C = np.diag(err[test] ** 2 + s2) + B @ cov @ B.T
        r = tau[test] - B @ b
        L = np.linalg.cholesky(C)
        z = np.linalg.solve(L, r)
        out.append(float(-0.5 * z @ z - np.sum(np.log(np.diag(L))) - 0.5 * r.size * np.log(2 * np.pi)))
    return out[0] - out[1]


def chi2nu_const(x, xe):
    """chi^2_nu of values about their weighted mean (NaN if fewer than 2)."""
    x, xe = np.asarray(x, float), np.asarray(xe, float)
    ok = np.isfinite(x) & np.isfinite(xe) & (xe > 0)
    if ok.sum() < 2:
        return np.nan
    w = xe[ok] ** -2
    return float(np.sum(w * (x[ok] - np.sum(w * x[ok]) / w.sum()) ** 2) / (ok.sum() - 1))


def chi2nu_const_by_band(x, xe, band):
    """Pooled chi^2_nu over bands, each band about its own weighted mean."""
    tot, dof = 0.0, 0
    for b in np.unique(band):
        m = band == b
        c = chi2nu_const(x[m], xe[m])
        if np.isfinite(c):
            tot += c * (m.sum() - 1)
            dof += m.sum() - 1
    return tot / dof if dof else np.nan


# ------------------------------------------------------------------ everything for one star
def oc_stats(series: dict, priors=None, n_harm_detect=1, p_min=400.0, with_red=True) -> dict:
    """All Level-2 statistics for one star. `priors`: {band: (mean [d], sd [d])} for the band offsets.

    Detection uses n_harm_detect harmonics (1 = circular); the 2-harmonic fit (moderate eccentricity) is also reported.
    """
    t, tau, err, band = (np.asarray(series[k], float) for k in ("t", "tau", "err", "band"))
    band = band.astype(int)
    P = float(series["P"])
    o = np.argsort(t)
    t, tau, err, band = t[o], tau[o], err[o], band[o]
    tau = robust_unwrap(t, tau, band, P)
    periods = period_grid(np.ptp(t), p_min=p_min)
    r1 = orbit_search(t, tau, err, band, periods, 1, priors)
    r2 = orbit_search(t, tau, err, band, periods, 2, priors)
    rd = r1 if n_harm_detect == 1 else r2
    A1, _ = orbit_amplitude(r1["beta"], r1["names"], r1["P_best"], 1)
    A2, A2f = orbit_amplitude(r2["beta"], r2["names"], r2["P_best"], 2)
    out = dict(D=rd["D"], P_best=rd["P_best"], D_circ=r1["D"], P_circ=r1["P_best"], amp_circ_s=A1 * DAY,
               D_2h=r2["D"], P_2h=r2["P_best"], amp_2h_s=A2 * DAY, dD_harm2=r2["D"] - r1["D"],
               amp_s=(A1 if n_harm_detect == 1 else A2) * DAY, jit0_s=rd["s0"] * DAY, jit1_s=rd["s1"] * DAY,
               alias_dD=alias_dD(rd["Dp"], periods, rd["P_best"]), n_season=int(t.size), baseline=float(np.ptp(t)),
               err_med_s=float(np.median(err) * DAY), has_M=bool((band > 0).any()), has_R=bool((band == 2).any()))
    for key in ("alpha", "coh"):
        if key in series and series[key] is not None:
            x, xe = np.asarray(series[key], float)[o], np.asarray(series[key + "_err"], float)[o]
            out[f"{key}_chi2nu"] = chi2nu_const_by_band(x, xe, band)
    if with_red:
        s, rw = red_null(t, tau, err, band, priors)
        out.update(s_red_s=s * DAY, rw_rms_s=rw * DAY)
    # predictive test: hold out the MACHO years before OGLE-II (HJD' < 450). Each MACHO band's offset must be tied to
    # OGLE in the training set: either OGLE-II seasons overlapping that band (both in 1997-2000) or a prior on the offset.
    test = (band > 0) & (t < 450)
    ok = test.sum() >= 3
    ogle_overlap = np.any((band == 0) & (t < 2000))
    for b in BANDS_EXTRA:
        if np.any(band == b):
            tied = (ogle_overlap and ((band == b) & ~test).sum() >= 1) or (priors is not None and b in priors)
            ok &= tied
    ok &= (~test).sum() >= 8
    if ok:
        # the orbital period is searched on the TRAINING seasons only (no information from the held-out seasons)
        tr = ~test
        rtr = orbit_search(t[tr], tau[tr], err[tr], band[tr], period_grid(np.ptp(t[tr]), p_min=p_min), n_harm_detect, priors)
        out["P_train"], out["D_train"] = rtr["P_best"], rtr["D"]
        out["pred_score"] = predictive_score(t, tau, err, band, test, rtr["P_best"], n_harm_detect, priors)
    out["tau_unwrapped"] = tau
    return out


# ------------------------------------------------------------------ common mode over the population
def residuals_h0(series, priors=None):
    """Residuals of a star's delays about its own H0 fit (quadratic + band offsets), white jitter profiled."""
    t, tau, err, band = (np.asarray(series[k], float) for k in ("t", "tau", "err", "band"))
    band = band.astype(int)
    tau = robust_unwrap(t, tau, band, float(series["P"]))
    X, names = design(t, band, np.average(t, weights=err ** -2))
    _, _, b, _ = profile_lnl(X, tau, err, names, default_s_grid(err), priors)
    return tau - X @ b


def common_mode(series_list, labels_fn, n_iter=5, min_stars=100, tol=1 / DAY, select=None):
    """Iterative per-(band, year) clipped mean of H0 residuals over the population, subtracted until converged.

    labels_fn(t) -> integer year label. `select`: boolean per series (well-behaved stars used for the estimate).
    Returns (dict {(band, year): correction [d]}, list of max |update| per iteration)."""
    cm, hist = {}, []
    sel = np.ones(len(series_list), bool) if select is None else np.asarray(select)
    for _ in range(n_iter):
        rows = []
        for s, use in zip(series_list, sel):
            if not use or len(s["t"]) < 6:
                continue
            s2 = dict(s)
            s2["tau"] = np.asarray(s["tau"]) - np.array([cm.get((int(b), int(y)), 0.0)
                                                        for b, y in zip(s["band"], labels_fn(np.asarray(s["t"])))])
            r = residuals_h0(s2)
            rows += list(zip(np.asarray(s["band"]).astype(int), labels_fn(np.asarray(s["t"])), r))
        if not rows:
            break
        b, y, r = map(np.asarray, zip(*rows))
        upd = {}
        for key in set(zip(b, y)):
            m = (b == key[0]) & (y == key[1])
            if m.sum() >= min_stars:
                x = r[m]
                med, mad = np.median(x), 1.4826 * np.median(np.abs(x - np.median(x)))
                x = x[np.abs(x - med) < 3 * mad] if mad > 0 else x
                upd[(int(key[0]), int(key[1]))] = float(np.mean(x))   # 3-sigma clipped mean (linear -> converges)
        old = dict(cm)
        for k, v in upd.items():
            cm[k] = cm.get(k, 0.0) + v
        cm = _project_out_degenerate(cm)
        hist.append(max(abs(cm[k] - old.get(k, 0.0)) for k in cm) if cm else 0.0)   # net change after projection
        if hist[-1] < tol:
            break
    return cm, hist


def _project_out_degenerate(cm):
    """Remove from the common-mode table the directions every star's H0 absorbs: a quadratic in time common to all bands
    and a constant per MACHO band. Such components cannot be measured (and would drift under iteration)."""
    if not cm:
        return cm
    keys = sorted(cm)
    b = np.array([k[0] for k in keys])
    y = np.array([k[1] for k in keys], float)
    v = np.array([cm[k] for k in keys])
    yc = (y - y.mean()) / max(y.std(), 1.0)
    cols = [np.ones_like(yc), yc, yc ** 2] + [(b == bb).astype(float) for bb in BANDS_EXTRA if np.any(b == bb)]
    X = np.column_stack(cols)
    coef, *_ = np.linalg.lstsq(X, v, rcond=None)
    v = v - X @ coef
    return dict(zip(keys, v))


def apply_common_mode(series, cm, labels_fn):
    s = dict(series)
    s["tau"] = np.asarray(series["tau"]) - np.array([cm.get((int(b), int(y)), 0.0)
                                                    for b, y in zip(series["band"], labels_fn(np.asarray(series["t"])))])
    return s
