"""Mock-data illustrations for the write-up (docs/writeup). All synthetic; outputs PNG in plots/writeup/.

  mock_orbit_delay.png     orbit of the RR Lyrae around the barycentre and the resulting light-travel delay tau(t)
  mock_lightcurve.png      the same pulsation template observed at two orbital phases: a pure time shift
  mock_oc_zoo.png          O-C diagrams: LTTE only; LTTE + Pdot; and the nuisances (Blazhko, jump, random walk)
  mock_measurement.png     OGLE-like cadence: per-season delays measured by the pipeline vs the true tau(t), and the fit
  mock_two_bands.png       two surveys/bands: templates differ in shape and phase -> constant offset between delay series

The template is the fitted OGLE I template of a typical RRab (Fourier coefficients hard-coded below); the orbit is
M1 = 0.65 Msun, M2 = 0.6 Msun, P_orb = 3000 d, e = 0.3, i = 90 deg unless stated.

Usage
-----
    python scripts/plot_writeup_mocks.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from rrlbin.ltte import a1sini_over_c, kepler_E, ltte_delay, oc_search, period_grid  # noqa: E402
from rrlbin.simulate import delay_jump, delay_random_walk  # noqa: E402
from rrlbin.timing import fit_timing, fourier_eval, year_labels  # noqa: E402

OUT = Path("plots/writeup")
DAY = 86400.0
P, T0 = 0.5712, 0.0
A = np.array([0.155, 0.068, 0.049, 0.024, 0.010, 0.008, 0.006, 0.004])
PH = np.array([0.0, 2.3, 4.5, 0.6, 2.9, 5.1, 1.3, 3.6])
COEF = np.empty(16)
COEF[0::2], COEF[1::2] = A * np.cos(PH), A * np.sin(PH)
M1, M2, PORB, ECC, OMEGA = 0.65, 0.6, 3000.0, 0.3, 1.0
AMP = a1sini_over_c(PORB, M1, M2) / DAY                    # [d]


def orbit_xy(t):
    Mn = 2 * np.pi * t / PORB
    E = kepler_E(np.mod(Mn, 2 * np.pi), ECC)
    nu = 2 * np.arctan2(np.sqrt(1 + ECC) * np.sin(E / 2), np.sqrt(1 - ECC) * np.cos(E / 2))
    r = (1 - ECC ** 2) / (1 + ECC * np.cos(nu))
    return r * np.cos(nu + OMEGA), r * np.sin(nu + OMEGA)


def fig_orbit():
    t = np.linspace(0, 2 * PORB, 1000)
    x, z = orbit_xy(t)
    tau = ltte_delay(t, PORB, AMP * DAY, ECC, OMEGA, 0.0)
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
    a1 = AMP * DAY * 299792.458 / 1.496e8
    ax[0].plot(x * a1, z * a1, "k-", lw=1)
    ax[0].plot(0, 0, "k+", ms=12)
    for tt, col in [(0.15 * PORB, "C3"), (0.62 * PORB, "C0")]:
        xx, zz = orbit_xy(np.array([tt]))
        ax[0].plot(xx * a1, zz * a1, "o", color=col, ms=9)
    ax[0].annotate("to observer", xy=(0, -1.6 * a1), xytext=(0, -0.9 * a1), ha="center",
                   arrowprops=dict(arrowstyle="->"))
    ax[0].set(aspect="equal", xlabel="x [AU]", ylabel="z (along line of sight, away from us) [AU]",
              title=f"RR Lyrae orbit about the barycentre\n(M1={M1}, M2={M2} Msun, P={PORB:.0f} d, e={ECC})")
    ax[1].plot(t, tau, "k-")
    for tt, col in [(0.15 * PORB, "C3"), (0.62 * PORB, "C0")]:
        ax[1].plot(tt, ltte_delay(np.array([tt]), PORB, AMP * DAY, ECC, OMEGA, 0.0), "o", color=col, ms=9)
    ax[1].axhline(0, c="0.6", lw=0.8)
    ax[1].set(xlabel="time [d]", ylabel="light-travel delay τ(t) = z/c [s]",
              title=f"delay: a1 sin i / c = {AMP * DAY:.0f} s (= {AMP * DAY / 60:.1f} min)")
    fig.tight_layout()
    fig.savefig(OUT / "mock_orbit_delay.png", dpi=130)


def fig_lightcurve():
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.3))
    t = np.linspace(0, 2 * P, 2000)
    for d, col, lab in [(+AMP, "C3", f"far side: τ = +{AMP * DAY:.0f} s"), (-AMP, "C0", f"near side: τ = −{AMP * DAY:.0f} s")]:
        ax[0].plot(t * 24, fourier_eval(COEF, (t - d - T0) / P), color=col, label=lab)
        ax[1].plot(t * 24 * 60, fourier_eval(COEF, (t - d - T0) / P), color=col)
    ax[0].invert_yaxis()
    ax[0].set(xlabel="time [h]", ylabel="Δ I [mag]", title="the same pulsation seen at two orbital phases")
    ax[0].legend(fontsize=9)
    x = np.linspace(0, 1, 2000)
    xr = x[np.argmin(np.gradient(fourier_eval(COEF, x)))]
    ax[1].set(xlim=((xr * P - 0.05) * 1440, (xr * P + 0.05) * 1440), xlabel="time [min]", ylabel="Δ I [mag]",
              title=f"zoom on the rising branch: shift = 2 × {AMP * DAY / 60:.1f} min")
    ax[1].invert_yaxis()
    fig.tight_layout()
    fig.savefig(OUT / "mock_lightcurve.png", dpi=130)


def fig_oc_zoo():
    rng = np.random.default_rng(3)
    t = np.linspace(-1000, 7500, 400)
    tl = ltte_delay(t, PORB, AMP * DAY, ECC, OMEGA, 0.0)
    pdot = 0.5 * (0.15 / 3.6525e8) / P * (t - 3000) ** 2 * DAY
    bl = 0.02 * P * np.sin(2 * np.pi * t / 1800 + 1) * DAY
    jp = delay_jump(t, 1.5e-5, 3500) * DAY
    rw = delay_random_walk(t, 900 / DAY, rng) * DAY
    panels = [("LTTE only (binary)", tl), ("LTTE + constant Ṗ (quadratic)", tl + pdot), ("constant Ṗ only", pdot),
              ("long-period Blazhko phase modulation", bl), ("abrupt period change (ΔP/P = 1.5e-5)", jp - jp.mean()),
              ("random-walk phase (900 s rms)", rw)]
    fig, axs = plt.subplots(2, 3, figsize=(15, 7), sharex=True)
    for ax, (tt, y) in zip(axs.ravel(), panels):
        ax.plot(t, y, "k-")
        ax.axhline(0, c="0.6", lw=0.8)
        ax.set(title=tt, ylabel="O−C [s]")
    for ax in axs[1]:
        ax.set_xlabel("HJD − 2450000")
    fig.suptitle("What an O−C diagram can look like (mock, noiseless)", y=1.0)
    fig.tight_layout()
    fig.savefig(OUT / "mock_oc_zoo.png", dpi=130)


def cadence(rng, years=range(2001, 2017), n=(40, 110)):
    ts = []
    for y in years:
        t0 = 245 + 365.25 * (y - 1998) + 60
        ts.append(t0 + np.sort(rng.uniform(0, 240, rng.integers(*n))))
    return np.concatenate(ts)


def fig_measurement():
    rng = np.random.default_rng(5)
    t = cadence(rng)
    tau_true = lambda x: ltte_delay(x, PORB, AMP * DAY, ECC, OMEGA, 0.0) / DAY
    m = 18.8 + fourier_eval(COEF, (t - tau_true(t) - T0) / P) + rng.normal(0, 0.06, t.size)
    e = np.full(t.size, 0.06)
    f = fit_timing(t, m, e, np.full(t.size, "O4"), P, T0, K=8, labels=year_labels(t))
    tau = f.tau - np.average(f.tau, weights=f.tau_err ** -2)
    tt = np.linspace(t.min(), t.max(), 800)
    tru = tau_true(tt)
    tru = tru - np.average(tau_true(f.t_season), weights=f.tau_err ** -2)
    r = oc_search(f.t_season, tau, f.tau_err, period_grid(np.ptp(f.t_season)))
    fig, ax = plt.subplots(1, 3, figsize=(16, 4.3), gridspec_kw=dict(width_ratios=[1.2, 1.6, 1]))
    sel = year_labels(t) == year_labels(t)[0]
    ax[0].plot(np.mod((t[sel] - T0) / P, 1), m[sel], "k.", ms=3)
    ax[0].invert_yaxis()
    ax[0].set(xlabel="phase (linear ephemeris)", ylabel="I [mag]", title=f"one season: {sel.sum()} epochs, σ = 0.06 mag")
    ax[1].plot(tt, tru * DAY, "k-", lw=1, label="true τ(t) (mean removed)")
    ax[1].errorbar(f.t_season, tau * DAY, f.tau_err * DAY, fmt="o", color="C0", label="measured per-season delay")
    ax[1].set(xlabel="HJD − 2450000", ylabel="O−C [s]", title="pipeline measurement on OGLE-like cadence")
    ax[1].legend(fontsize=8)
    ax[2].plot(period_grid(np.ptp(f.t_season)), r["Dp"], "k-")
    ax[2].axvline(PORB, c="C3", ls=":", label="true P_orb")
    ax[2].axhline(40, c="C3", ls="--", lw=0.8)
    ax[2].set(xscale="log", xlabel="trial P_orb [d]", ylabel="D(P)", title=f"orbit search: D = {r['D']:.0f} at P = {r['P_best']:.0f} d")
    ax[2].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(OUT / "mock_measurement.png", dpi=130)
    return float(np.median(f.tau_err) * DAY), r


def fig_two_bands():
    x = np.linspace(0, 1, 800)
    coefB = COEF.copy()
    k = np.arange(1, 9)
    lag = +0.034 * 2 * np.pi                                     # rotating phases by +k*lag moves the maximum EARLIER by 0.034 cycles
    a, b = coefB[0::2] * 1.9, coefB[1::2] * 1.9
    c, s = np.cos(k * lag), np.sin(k * lag)
    coefB[0::2], coefB[1::2] = a * c + b * s, b * c - a * s
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.2))
    ax[0].plot(x, fourier_eval(COEF, x), "C0", label="OGLE I template")
    ax[0].plot(x, fourier_eval(coefB, x), "C1", label="MACHO B template (larger amplitude, earlier maximum)")
    ax[0].invert_yaxis()
    ax[0].set(xlabel="pulsation phase", ylabel="Δ mag", title="same star, two bands")
    ax[0].legend(fontsize=8)
    rng = np.random.default_rng(1)
    t = np.linspace(-1000, 7400, 23)
    tr = ltte_delay(t, PORB, AMP * DAY, ECC, OMEGA, 0.0)
    mac = t < 1600
    ogl = t > 400
    ax[1].plot(t[ogl], tr[ogl] + rng.normal(0, 150, ogl.sum()), "o", color="C0", label="OGLE I delays")
    ax[1].plot(t[mac], tr[mac] - 0.034 * P * DAY + rng.normal(0, 150, mac.sum()), "s", color="C1",
               label=f"MACHO B delays (offset −0.034 P = {-0.034 * P * DAY:.0f} s)")
    td = np.linspace(-1000, 7400, 600)
    ax[1].plot(td, ltte_delay(td, PORB, AMP * DAY, ECC, OMEGA, 0.0), "k-", lw=0.8, label="true τ(t)")
    ax[1].axvspan(400, 1600, color="0.9", zorder=0, label="overlap years fix the offset")
    ax[1].set(xlabel="HJD − 2450000", ylabel="O−C [s]", title="delay series: a constant band offset")
    ax[1].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(OUT / "mock_two_bands.png", dpi=130)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    fig_orbit()
    fig_lightcurve()
    fig_oc_zoo()
    err, r = fig_measurement()
    fig_two_bands()
    print(f"amp = {AMP * DAY:.0f} s; mock measurement: median season error {err:.0f} s, D = {r['D']:.1f}, P_best = {r['P_best']:.0f} d")


if __name__ == "__main__":
    main()
