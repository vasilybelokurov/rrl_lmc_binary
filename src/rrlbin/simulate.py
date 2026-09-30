"""Simulated RR Lyrae light curves on real cadences, with timing nuisances and LTTE orbits.

m(t) = zp_seg + A(t) * T((t - tau(t) - T0) / P) + noise,

where T is a Fourier template (from a fit to the real star), A(t) an amplitude modulation and
tau(t) the total delay: quadratic ephemeris (Pdot) + Blazhko phase modulation + abrupt period change
+ random-walk phase + LTTE. All delays in days.
"""
from __future__ import annotations

import numpy as np

from .ltte import ltte_delay
from .timing import fourier_eval


def delay_pdot(t, P, pdot_d_per_Myr, t_ref):
    """O-C of a constant period change, tau = (1/2) (Pdot/P) (t - t_ref)^2 (a later arrival as P grows). Pdot in d/Myr."""
    pdot = pdot_d_per_Myr / 3.6525e8
    return 0.5 * pdot / P * (t - t_ref) ** 2


def delay_jump(t, dP_over_P, t_break):
    """Abrupt period change dP/P at t_break: O-C grows linearly after the break."""
    return dP_over_P * np.clip(t - t_break, 0, None)


def delay_random_walk(t, rms_over_baseline, rng):
    """Brownian phase wander, normalized so that its rms deviation from a straight line over the data is ~rms."""
    ts = np.sort(np.unique(t))
    steps = rng.normal(0, 1, ts.size) * np.sqrt(np.diff(np.r_[ts[0], ts]))
    w = np.cumsum(steps)
    w = w - np.polyval(np.polyfit(ts, w, 1), ts)
    w *= rms_over_baseline / max(w.std(), 1e-12)
    return np.interp(t, ts, w)


def simulate_lc(t, err, seg, coef, zp_by_seg, P, T0, rng, noise_scale=1.0,
                pdot=0.0, t_ref=5000.0, blazhko=None, jump=None, rw_rms=0.0, ltte=None):
    """Simulate magnitudes at the real epochs t.

    blazhko : dict(P_B, eps_A, eps_phi [cycles], psi_A, psi_phi) or None
    jump    : dict(dP_over_P, t_break) or None
    ltte    : dict(P_orb, amp [d], e, omega, t_peri) or None
    """
    tau = delay_pdot(t, P, pdot, t_ref)
    A = np.ones_like(t)
    if blazhko is not None:
        b = blazhko
        ph = 2 * np.pi * t / b["P_B"]
        A = A + b["eps_A"] * np.sin(ph + b["psi_A"])
        tau = tau + b["eps_phi"] * P * np.sin(ph + b["psi_phi"])
    if jump is not None:
        tau = tau + delay_jump(t, jump["dP_over_P"], jump["t_break"])
    if rw_rms > 0:
        tau = tau + delay_random_walk(t, rw_rms, rng)
    if ltte is not None:
        tau = tau + ltte_delay(t, **ltte)
    zp = np.array([zp_by_seg[s] for s in seg])
    m = zp + A * fourier_eval(coef, (t - tau - T0) / P)
    return m + rng.normal(0, 1, t.size) * err * noise_scale, tau
