"""Pulsation timing: Fourier template and per-season time delays ("O-C") of RR Lyrae light curves.

Model for one star and one band (magnitudes):

    m(t) = z_s + sum_{k=1}^{K} [a_k cos(2 pi k phi) + b_k sin(2 pi k phi)],
    phi  = (t - tau_j - T0) / P,

where s labels the survey segment (separate zero point), j the observing season, and tau_j is the
time delay of season j. A single tau_j shifts all harmonics coherently, which is the signature of a
light-travel-time effect (LTTE) or a pure phase (period) change. It is not the signature of
Blazhko modulation, which changes the harmonic amplitudes and relative phases.

Conventions: t in days (HJD - 2450000). tau > 0 means the pulsation arrives late (the star is farther
away than at the reference). The Fisher error of tau_j is [sum_i (dm/dt)_i^2 / sigma_i^2]^(-1/2),
scaled by sqrt(chi2_nu) when that exceeds 1.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

TWO_PI = 2.0 * np.pi


# ---------------------------------------------------------------- Fourier basis
def fourier_design(phi: np.ndarray, K: int) -> np.ndarray:
    """Design matrix [cos 2pi phi, sin 2pi phi, ..., cos 2pi K phi, sin 2pi K phi] (no constant), shape (N, 2K)."""
    k = np.arange(1, K + 1)
    arg = TWO_PI * np.outer(phi, k)
    X = np.empty((phi.size, 2 * K))
    X[:, 0::2] = np.cos(arg)
    X[:, 1::2] = np.sin(arg)
    return X


def fourier_eval(coef: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """Template value (without zero point) at phase phi; coef = [a1, b1, a2, b2, ...]."""
    return fourier_design(phi, coef.size // 2) @ coef


def fourier_dphi(coef: np.ndarray, phi: np.ndarray) -> np.ndarray:
    """d(template)/d(phi)."""
    K = coef.size // 2
    k = np.arange(1, K + 1)
    arg = TWO_PI * np.outer(phi, k)
    a, b = coef[0::2], coef[1::2]
    return TWO_PI * (-np.sin(arg) * k) @ a + TWO_PI * (np.cos(arg) * k) @ b


def harmonic_amp_phase(coef: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(A_k, phi_k) with a_k cos x + b_k sin x = A_k cos(x - phi_k)."""
    a, b = coef[0::2], coef[1::2]
    return np.hypot(a, b), np.arctan2(b, a)


# ---------------------------------------------------------------- seasons
def season_labels(t: np.ndarray, gap: float = 60.0) -> np.ndarray:
    """Integer season label per epoch: a new season starts after a gap > `gap` days. t need not be sorted."""
    order = np.argsort(t)
    brk = np.concatenate([[0], np.cumsum(np.diff(t[order]) > gap)])
    lab = np.empty(t.size, dtype=int)
    lab[order] = brk
    return lab


def year_labels(t: np.ndarray, phase0: float = 245.0) -> np.ndarray:
    """Season = observing year, with the boundary at HJD' mod 365.25 = phase0 (the middle of the OGLE LMC seasonal gap,
    measured from OGLE-III/IV epochs). Needed for MACHO, which observed the LMC nearly year-round."""
    return np.floor((np.asarray(t) - phase0) / 365.25).astype(int)


# ---------------------------------------------------------------- fits
@dataclass
class TimingFit:
    """Result of fit_timing for one star/band."""
    P: float
    T0: float
    coef: np.ndarray            # template Fourier coefficients (2K)
    zp: np.ndarray              # zero point per segment
    seg_names: np.ndarray
    season: np.ndarray          # unique season labels
    t_season: np.ndarray        # Fisher-weighted epoch of each season's delay
    tau: np.ndarray             # delay per season [d]
    tau_err: np.ndarray         # scaled Fisher error [d]
    alpha: np.ndarray           # amplitude scale per season (1 = template amplitude)
    alpha_err: np.ndarray
    dm: np.ndarray              # mean-magnitude offset per season
    chi2nu_season: np.ndarray
    n_season: np.ndarray
    seg_season: np.ndarray      # segment of each season
    mask: np.ndarray            # epochs kept after clipping
    chi2nu: float
    coh: np.ndarray = None      # per season: delay of the fundamental minus delay of the higher harmonics [d]
    coh_err: np.ndarray = None
    unstable: np.ndarray = None  # labels of seasons dropped by the stability check


def _template_fit(t, m, w, seg_idx, nseg, tau_ep, P, T0, K):
    """Weighted LSQ for zero points (per segment) + Fourier coefficients, at fixed per-epoch delays."""
    phi = (t - tau_ep - T0) / P
    Z = np.zeros((t.size, nseg))
    Z[np.arange(t.size), seg_idx] = 1.0
    X = np.hstack([Z, fourier_design(phi, K)])
    sw = np.sqrt(w)
    beta, *_ = np.linalg.lstsq(X * sw[:, None], m * sw, rcond=None)
    return beta[:nseg], beta[nseg:]


def _season_shift(t, m, w, zp_ep, coef, P, T0, tau0, n_iter=10, grid=True):
    """Gauss-Newton fit of (tau, dm, alpha) for one season at a fixed template:
    m = zp + dm + alpha * T((t - tau - T0)/P). alpha (amplitude scale) is a Blazhko/blending diagnostic.

    Returns dict: tau, dm, alpha, var_tau, var_alpha, chi2, t_eff (Fisher-weighted epoch of the delay)."""
    if grid:  # coarse global search over one full cycle, to avoid wrong local minima
        taus = tau0 + P * np.linspace(-0.5, 0.5, 101)[:-1]
        k = np.arange(1, coef.size // 2 + 1)
        r0 = m - zp_ep
        chi = np.empty(taus.size)
        for i0 in range(0, taus.size, 25):          # vectorized over blocks of grid points (memory ~ 25 x n x K)
            arg = TWO_PI * ((t[None, :] - taus[i0:i0 + 25, None] - T0) / P)[..., None] * k
            T = np.cos(arg) @ coef[0::2] + np.sin(arg) @ coef[1::2]
            chi[i0:i0 + 25] = ((r0[None, :] - T) ** 2 * w[None, :]).sum(1)
        tau0 = taus[int(np.argmin(chi))]
    p = np.array([tau0, 0.0, 1.0])

    def jac(p):
        phi = (t - p[0] - T0) / P
        T = fourier_eval(coef, phi)
        r = m - zp_ep - p[1] - p[2] * T
        J = np.column_stack([-p[2] * fourier_dphi(coef, phi) / P, np.ones_like(t), T])
        return r, J

    for _ in range(n_iter):
        r, J = jac(p)
        step = np.linalg.solve(J.T @ (J * w[:, None]), J.T @ (w * r))
        p = p + step
        if abs(step[0]) < 1e-8:
            break
    r, J = jac(p)
    cov = np.linalg.inv(J.T @ (J * w[:, None]))
    wt = w * J[:, 0] ** 2
    return dict(tau=p[0], dm=p[1], alpha=p[2], var_tau=cov[0, 0], var_alpha=cov[2, 2],
                chi2=float(np.sum(w * r ** 2)), t_eff=float(np.sum(wt * t) / np.sum(wt)))


def _season_coherence(t, m, w, zp_ep, coef, P, T0, tau, dm, alpha, n_iter=10):
    """Harmonic-coherence test for one season: fit separate delays to the fundamental (k = 1) and to the higher
    harmonics (k >= 2) of the template,

        m = zp + dm + alpha [T_1((t - tau_1 - T0)/P) + T_h((t - tau_h - T0)/P)],

    starting from the common delay. A light-travel delay shifts every harmonic by the same TIME, so tau_1 - tau_h = 0;
    changes of the light-curve shape (Blazhko modulation, changing blends) generally make it non-zero.
    Returns (tau_1 - tau_h, its error) in days."""
    c1 = coef.copy()
    c1[2:] = 0.0
    ch = coef.copy()
    ch[:2] = 0.0
    p = np.array([tau, tau, dm, alpha])
    for _ in range(n_iter):
        f1, fh = (t - p[0] - T0) / P, (t - p[1] - T0) / P
        T1, Th = fourier_eval(c1, f1), fourier_eval(ch, fh)
        r = m - zp_ep - p[2] - p[3] * (T1 + Th)
        J = np.column_stack([-p[3] * fourier_dphi(c1, f1) / P, -p[3] * fourier_dphi(ch, fh) / P, np.ones_like(t), T1 + Th])
        H = J.T @ (J * w[:, None])
        step = np.linalg.solve(H, J.T @ (w * r))
        p = p + step
        if np.max(np.abs(step[:2])) < 1e-8:
            break
    cov = np.linalg.inv(H)
    return p[0] - p[1], float(np.sqrt(cov[0, 0] + cov[1, 1] - 2 * cov[0, 1]))


def _season_minima(t, m, w, zp_ep, coef, P, T0, n_grid=100, n_keep=3):
    """Distinct local minima of one season's delay chi^2 (dm, alpha free) at a fixed template: coarse grid over one cycle,
    the deepest `n_keep` grid minima refined by Gauss-Newton; returns the _season_shift results sorted by chi^2 (duplicates
    within 0.02 P merged)."""
    taus = P * np.linspace(-0.5, 0.5, n_grid + 1)[:-1]
    k = np.arange(1, coef.size // 2 + 1)
    chi = np.empty(n_grid)
    for i0 in range(0, n_grid, 25):
        arg = TWO_PI * ((t[None, :] - taus[i0:i0 + 25, None] - T0) / P)[..., None] * k
        T = np.cos(arg) @ coef[0::2] + np.sin(arg) @ coef[1::2]
        # profile dm and alpha analytically per grid point (weighted LSQ of r = dm + alpha T)
        r = (m - zp_ep)[None, :]
        sw, swT, swTT = w.sum(), (w * T).sum(1), (w * T * T).sum(1)
        swr, swTr = (w * r).sum(1), (w * T * r).sum(1)
        det = sw * swTT - swT ** 2
        dm = (swTT * swr - swT * swTr) / det
        al = (sw * swTr - swT * swr) / det
        chi[i0:i0 + 25] = (w * (r - dm[:, None] - al[:, None] * T) ** 2).sum(1)
    loc = [i for i in range(n_grid) if chi[i] <= chi[i - 1] and chi[i] <= chi[(i + 1) % n_grid]]
    loc = sorted(loc, key=lambda i: chi[i])[:n_keep]
    sols = []
    for i in loc:
        r = _season_shift(t, m, w, zp_ep, coef, P, T0, taus[i], grid=False)
        if 0.3 < r["alpha"] < 3 and not any(abs((r["tau"] - q["tau"] + P / 2) % P - P / 2) < 0.02 * P for q in sols):
            sols.append(r)
    return sorted(sols, key=lambda q: q["chi2"]) or [_season_shift(t, m, w, zp_ep, coef, P, T0, 0.0, grid=True)]


def fit_timing(t, m, err, seg, P, T0, K=8, gap=60.0, n_outer=50, clip=4.0, min_season=15, tol=1e-7, labels=None,
               gross=8.0, grid_every=3, stability=True, tol_cycle=1e-5) -> TimingFit:
    """Joint template + per-season delay fit for one star in one band.

    Parameters
    ----------
    t, m, err : arrays
        Epochs [HJD-2450000], magnitudes, errors.
    seg : array of str
        Survey segment of each epoch (separate zero point per segment, shared template shape).
    P, T0 : float
        Reference ephemeris (period [d], epoch [HJD-2450000]); delays are relative to it.
    K : int
        Fourier order.
    gap : float
        Season-break gap [d] (ignored if `labels` is given).
    labels : array of int, optional
        Explicit season label per epoch (e.g. year_labels(t)).
    n_outer : int
        Maximum template <-> delay iterations.
    tol : float
        Convergence: stop when the mask is unchanged and max |change of tau| < tol [d].
    tol_cycle : float
        After 10 iterations, also stop when max |change of tau| < tol_cycle [d] (0.86 s) although the clipping mask still
        flips (a limit cycle of points at the threshold; found on the 2010-2026 light curves, amplitude ~0.6 s).
    clip : float
        Sigma-clipping threshold on normalized residuals (after each outer iteration).
    min_season : int
        Minimum epochs for a season to get its own delay (smaller seasons are dropped).
    gross : float
        Clipping threshold of the FIRST iteration (gross outliers only): the first template is blurred by the unaligned delays,
        so a tight clip there removes rising-branch points (the most timing information) and can lock a season into a
        wrong minimum. Later iterations use `clip`. Clipped points can re-enter at every iteration.
    grid_every : int
        Coarse grid search over one cycle for every season in the first `grid_every` iterations (v3: only the first; 0/False
        = v3 behaviour): escapes wrong minima reached in the first iteration (blurred template, outliers not yet clipped).
        The stability check repeats the grid search once more at the final template.
    stability : bool
        After convergence, each season (only > `gross` sigma points removed, final template) is searched for distinct chi^2
        minima (_season_minima). The season is UNSTABLE and DROPPED (labels returned in `unstable`) if (a) the fitted delay
        is > 3 sigma from the deepest minimum (clipping- or start-dependent solution), or (b) the two deepest minima are
        > 3 sigma apart but differ by Delta chi^2 < 9 (after scaling by chi2_nu): an ambiguous season (found in sparse
        seasons of the 2010-2026 light curves).

    Each season also gets a mean-magnitude offset (absorbs slow photometric drifts and blending changes).
    """
    t, m, err = map(np.asarray, (t, m, err))
    seg = np.asarray(seg)
    seg_names, seg_idx = np.unique(seg, return_inverse=True)
    lab = season_labels(t, gap) if labels is None else np.asarray(labels)
    seasons, cnt = np.unique(lab, return_counts=True)
    keep = np.isin(lab, seasons[cnt >= min_season])
    w = 1.0 / err ** 2
    mask = keep.copy()
    tau_s = {s: 0.0 for s in seasons[cnt >= min_season]}
    for it in range(n_outer):
        tau_prev = dict(tau_s)
        tau_ep = np.array([tau_s.get(s, 0.0) for s in lab])
        zp, coef = _template_fit(t[mask], m[mask], w[mask], seg_idx[mask], seg_names.size, tau_ep[mask], P, T0, K)
        res = {}
        for s in tau_s:
            sel = mask & (lab == s)
            if sel.sum() < min_season:
                continue
            res[s] = _season_shift(t[sel], m[sel], w[sel], zp[seg_idx[sel]], coef, P, T0, tau_s[s], grid=(it == 0 or it < int(grid_every)))
            tau_s[s] = res[s]["tau"]
        # sigma clipping on the full model
        tau_ep = np.array([tau_s.get(s, 0.0) for s in lab])
        dm_ep = np.array([res[s]["dm"] if s in res else 0.0 for s in lab])
        al_ep = np.array([res[s]["alpha"] if s in res else 1.0 for s in lab])
        model = zp[seg_idx] + dm_ep + al_ep * fourier_eval(coef, (t - tau_ep - T0) / P)
        z = (m - model) / err
        cl = max(clip, gross) if it == 0 else clip
        new = keep & np.isin(lab, list(res)) & (np.abs(z) < cl * max(1.0, 1.4826 * np.median(np.abs(z[mask]))))
        dtau = max(abs(tau_s[s] - tau_prev[s]) for s in res)
        if it > 0 and np.array_equal(new, mask) and dtau < tol:
            break
        if it >= 10 and dtau < tol_cycle:    # clipping-set limit cycle (points flipping at the threshold): delays settled
            mask = new
            break
        mask = new

    unstable = []
    if stability:
        sc_all = max(1.0, 1.4826 * np.median(np.abs(z[mask])))
        for s in list(res):
            sel = keep & (lab == s) & (np.abs(z) < gross * sc_all)
            if sel.sum() < min_season:
                continue
            sols = _season_minima(t[sel], m[sel], w[sel], zp[seg_idx[sel]], coef, P, T0)
            n8 = int(sel.sum())
            best = sols[0]
            scl = max(best["chi2"] / max(n8 - 3, 1), 1.0)
            s_best = np.sqrt(best["var_tau"] * scl)
            wrap = lambda x: (x + P / 2) % P - P / 2
            ambiguous = (len(sols) > 1 and (sols[1]["chi2"] - best["chi2"]) / scl < 9.0
                         and abs(wrap(sols[1]["tau"] - best["tau"])) > 3 * s_best)
            off_best = abs(wrap(res[s]["tau"] - best["tau"])) > 3 * s_best
            if ambiguous or off_best:
                unstable.append(s)
                del res[s]
                mask &= lab != s
    ss = np.array(sorted(res))
    n = np.array([int((mask & (lab == s)).sum()) for s in ss])
    chi2 = np.array([res[s]["chi2"] for s in ss])
    chi2nu = chi2 / np.maximum(n - 3, 1)
    scale = np.maximum(chi2nu, 1.0)
    seg_season = np.array([seg_names[np.bincount(seg_idx[mask & (lab == s)]).argmax()] for s in ss])
    chi_tot = chi2.sum() / max(mask.sum() - 2 * K - seg_names.size - 3 * ss.size, 1)
    g = lambda k: np.array([res[s][k] for s in ss])
    # harmonic coherence per season (at the final template); error scaled like tau's
    coh = np.array([_season_coherence(t[mask & (lab == s)], m[mask & (lab == s)], w[mask & (lab == s)],
                                      zp[seg_idx[mask & (lab == s)]], coef, P, T0, res[s]["tau"], res[s]["dm"],
                                      res[s]["alpha"]) for s in ss]).reshape(-1, 2)
    # Gauge: a constant delay is degenerate with the template phase. Fix it so that the fundamental harmonic of the
    # template has phase 0: T(phi) = A1 cos(2 pi phi') + ..., phi' = (t - tau' - T0)/P, tau' = tau + ph1 P / (2 pi).
    # Delays are then physical (time of the fundamental's maximum light relative to T0) and comparable across bands.
    _, ph = harmonic_amp_phase(coef)
    shift = ph[0] * P / TWO_PI
    k = np.arange(1, coef.size // 2 + 1)
    c, sn = np.cos(k * ph[0]), np.sin(k * ph[0])
    a, b = coef[0::2].copy(), coef[1::2].copy()
    coef = coef.copy()
    coef[0::2], coef[1::2] = a * c + b * sn, b * c - a * sn
    tau_g = g("tau") + shift
    tau_g = tau_g - P * np.round((np.median(tau_g) - 0.0) / P)   # bring the series near zero (whole cycles)
    return TimingFit(P=P, T0=T0, coef=coef, zp=zp, seg_names=seg_names, season=ss,
                     t_season=g("t_eff"), tau=tau_g, tau_err=np.sqrt(g("var_tau") * scale),
                     alpha=g("alpha"), alpha_err=np.sqrt(g("var_alpha") * scale), dm=g("dm"),
                     chi2nu_season=chi2nu, n_season=n, seg_season=seg_season, mask=mask, chi2nu=float(chi_tot),
                     coh=coh[:, 0], coh_err=coh[:, 1] * np.sqrt(scale), unstable=np.array(unstable, int))


def unwrap_delays(tau: np.ndarray, P: float) -> np.ndarray:
    """Remove whole-cycle jumps between consecutive seasons (delays are defined modulo P)."""
    return tau[0] + np.unwrap(TWO_PI * (tau - tau[0]) / P) * P / TWO_PI


def delays_fixed_template(t, m, err, labels, coef, P, T0, min_season=8, clip=5.0, n_clip=3, zp=None, full=False):
    """Per-season delays of a sparse light curve against a FIXED template shape (e.g. another band's well-measured template,
    for Gaia epoch photometry): each season fits (tau, dm, alpha) via _season_shift with a common zero point removed.
    zp: fixed zero point (e.g. the segment zero point of a fit_timing result); default: estimated from the medians.
    Returns (t_eff, tau, tau_err) per season (errors scaled by sqrt(max(chi2_nu, 1))); with full=True also
    (alpha, alpha_err, n_used, chi2nu, label)."""
    t, m, err, labels = map(np.asarray, (t, m, err, labels))
    if zp is None:
        zp = np.median(m) - np.median(fourier_eval(coef, (t - T0) / P))
    w = 1 / err ** 2
    out = []
    for s in np.unique(labels):
        sel = labels == s
        if sel.sum() < min_season:
            continue
        keep = sel.copy()
        for _ in range(n_clip):
            r = _season_shift(t[keep], m[keep], w[keep], np.full(keep.sum(), zp), coef, P, T0, 0.0, grid=True)
            mod = zp + r["dm"] + r["alpha"] * fourier_eval(coef, (t - r["tau"] - T0) / P)
            z = (m - mod) / err
            new = sel & (np.abs(z) < clip * max(1.0, 1.4826 * np.median(np.abs(z[keep]))))
            if np.array_equal(new, keep):
                break
            keep = new
        chi2nu = r["chi2"] / max(keep.sum() - 3, 1)
        sc = max(chi2nu, 1.0)
        out.append((r["t_eff"], r["tau"], np.sqrt(r["var_tau"] * sc), r["alpha"], np.sqrt(r["var_alpha"] * sc),
                    int(keep.sum()), chi2nu, s))
    nret = 8 if full else 3
    return tuple(np.array(x) for x in list(zip(*out))[:nret]) if out else (np.empty(0),) * nret
