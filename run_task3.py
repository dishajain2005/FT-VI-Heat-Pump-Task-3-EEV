"""
run_task3.py -- STUDENT 3: all EEV studies, CSV results, figures and the hand-over file.
Usage:  python run_task3.py      (needs CoolProp for report numbers; same backend logic as run_all.py)
Outputs go to  results_task3/  and  figures_task3/
"""
import os, csv, json, math
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.optimize import brentq

from config import BRIEF, ASSUME, SENS_AMBIENTS
from co2_props import get_props
import baseline_cycle as bc
import flash_tank_vi as fv
import eev
import ftvi_eev as fe

props = get_props()
HERE = os.path.dirname(os.path.abspath(__file__))
RES, FIG = os.path.join(HERE, "results_task3"), os.path.join(HERE, "figures_task3")
os.makedirs(RES, exist_ok=True); os.makedirs(FIG, exist_ok=True)

# ----------------------------- plot style ------------------------------------------------
SURF, INK, INK2, MUTED, GRID, AXIS = "#fcfcfb", "#0b0b0b", "#52514e", "#898781", "#e1e0d9", "#c3c2b7"
CAT = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#6250d6", "#e34948"]
BLUE_ORD = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]     # ordinal ramp, light -> dark
plt.rcParams.update({
    "font.size": 9, "figure.dpi": 160, "savefig.dpi": 200,
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": AXIS, "axes.labelcolor": INK2, "text.color": INK,
    "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelcolor": INK2, "ytick.labelcolor": INK2,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 1.8,
    "legend.frameon": False, "axes.titlesize": 9.5, "axes.titleweight": "bold",
})


def save_csv(name, rows, keys=None):
    keys = keys or [k for k in rows[0].keys() if not isinstance(rows[0][k], (dict, list))]
    with open(os.path.join(RES, name), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow({k: (round(v, 6) if isinstance(v, float) else v) for k, v in r.items() if k in keys})


def limit_line(ax, y, text, where="right", va="bottom"):
    ax.axhline(y, color=MUTED, lw=1.0, ls=(0, (4, 3)), zorder=1)
    x = 0.99 if where == "right" else 0.01
    ax.text(x, y, text, transform=ax.get_yaxis_transform(), ha=where, va=va, color=INK2, fontsize=7.5)


def end_label(ax, x, y, text, color=INK2, dy=0):
    ax.annotate(text, (x, y), xytext=(4, dy), textcoords="offset points", va="center", fontsize=7.5, color=color)


def fig_save(fig, name):
    fig.savefig(os.path.join(FIG, name), bbox_inches="tight")
    plt.close(fig)


R = {"backend": props.name}
PHI = np.round(np.linspace(0, 1, 51), 4)
PFT_SET = [3.5, 4.0, 4.5, 5.0, 6.0]
COLS = ["phi", "P_FT", "P_port", "T_evap", "N", "x", "dP_eev", "choked", "m_suc", "m_valve", "m_inj", "m_tot",
        "xi", "xi_max", "limited_by", "phi_sat", "Q_h", "W", "COP_h", "T_dis", "T3", "x9",
        "energy_balance_err", "tank_energy_residual", "feasible"]

# ======================================================================================
# 0. Design point, valve size, validation
# ======================================================================================
V = fe.V_disp(); CdA = fe.CdA_design()
R["V_disp_cc"] = V * 1e6
R["CdA_mm2"] = CdA * 1e6
R["A_max_mm2"] = CdA / eev.EEV["Cd"] * 1e6
R["d_eq_mm"] = math.sqrt(4 * R["A_max_mm2"] / math.pi)
R["EEV"] = {k: v for k, v in eev.EEV.items()}
R["R_crit"] = eev.R_CRIT
A = bc.baseline(V_disp=V)
R["baseline"] = {k: A[k] for k in ("Q_h", "W", "COP_h", "T_dis", "m_dot", "P_evap")}
D = fe.ftvi_eev()                                          # brief's design point: phi = 50 %, P_FT = 4.5 MPa
R["design"] = {k: D[k] for k in COLS}

val = []
for p in (3.5, 4.0, 4.5, 5.0, 6.0):
    s2 = fv.ftvi(V, P_inj=p)
    s3 = fe.ftvi_eev(phi=1.0, P_FT=p, P_port=p - 1e-7, CdA_max=1.0)
    val.append(dict(P_FT=p, **{f"S2_{k}": s2[k] for k in ("xi", "Q_h", "W", "COP_h", "T_dis")},
                    **{f"S3_{k}": s3[k] for k in ("xi", "Q_h", "W", "COP_h", "T_dis")},
                    max_rel_diff=max(abs(s3[k] / s2[k] - 1) for k in ("xi", "Q_h", "W", "COP_h", "T_dis"))))
save_csv("T3_validation_vs_student2.csv", val)
R["validation_max_rel_diff"] = max(v["max_rel_diff"] for v in val)
Z = fe.ftvi_eev(phi=0.0)
R["closed_vs_baseline"] = dict(Q_h=(Z["Q_h"], A["Q_h"]), W=(Z["W"], A["W"]), COP=(Z["COP_h"], A["COP_h"]),
                               T_dis=(Z["T_dis"], A["T_dis"]))

# design-point state table (phi = 50 %)
st = [("1", "Suction", D["P_evap"], D["T1"], D["h1"], "SH vapour"),
      ("2", "End of stage 1", D["P_port"], D["T2"], D["h2"], "SH vapour"),
      ("7", "Flash-tank vapour", D["P_FT"], D["T_FT"], D["h7"], "x = 1"),
      ("7a", "After injection EEV", D["P_port"], D["T7a"], D["h7a"],
       "x = %.3f" % ((D["h7a"] - props.hf(D["P_port"])) / (props.hg(D["P_port"]) - props.hf(D["P_port"])))
       if D["h7a"] < props.hg(D["P_port"]) else "SH vapour"),
      ("3", "After mixing", D["P_port"], D["T3"], D["h3"], "SH vapour"),
      ("4", "Discharge", D["P_dis"], D["T4"], D["h4"], "SC fluid"),
      ("5", "Gas-cooler outlet", D["P_dis"], BRIEF["T_gc_out"], D["h5"], "SC fluid"),
      ("6", "EEV-1 outlet", D["P_FT"], D["T_FT"], D["h6"], f"x = {D['x']:.3f}"),
      ("8", "Tank outlet to EEV-2", D["P_FT"], D["T_FT"], D["h8"], f"x = {D['x8']:.3f}"),
      ("9", "Evaporator inlet", D["P_evap"], props.Tsat(D["P_evap"]), D["h9"], f"x = {D['x9']:.3f}")]
with open(os.path.join(RES, "T3_state_table_phi50.csv"), "w", newline="") as f:
    w = csv.writer(f); w.writerow(["state", "description", "P_MPa", "T_C", "h_kJkg", "phase"])
    for r in st: w.writerow([r[0], r[1], round(r[2], 4), round(r[3], 2), round(r[4], 2), r[5]])
R["state_table"] = st

# ======================================================================================
# 1. EEV opening sweep at several flash-tank pressures (-10 C ambient, 4000 rpm)
# ======================================================================================
S1 = {p: [fe.ftvi_eev(phi=f, P_FT=p) for f in PHI] for p in PFT_SET}
save_csv("T3_phi_sweep_PFT.csv", [r for p in PFT_SET for r in S1[p]], COLS)
R["S1_summary"] = {str(p): dict(phi_sat=S1[p][0]["phi_sat"], xi_max=S1[p][0]["xi_max"],
                                Q0=S1[p][0]["Q_h"], Q1=S1[p][-1]["Q_h"], COP0=S1[p][0]["COP_h"], COP1=S1[p][-1]["COP_h"],
                                T0=S1[p][0]["T_dis"], T1=S1[p][-1]["T_dis"], xi1=S1[p][-1]["xi"],
                                choked=S1[p][-1]["choked"], dP=S1[p][-1]["dP_eev"]) for p in PFT_SET}
R["S1_table_45"] = [r for r in S1[4.5] if abs(r["phi"] * 10 - round(r["phi"] * 10)) < 1e-9]

# ======================================================================================
# 2. Ambient sensitivity (P_FT = 4.5 MPa, 4000 rpm)
# ======================================================================================
S2 = {}
for Ta in SENS_AMBIENTS:
    S2[Ta] = [fe.ftvi_eev(phi=f, T_evap=bc.evap_temperature(Ta)) for f in PHI]
    for r in S2[Ta]: r["T_amb"] = Ta
save_csv("T3_phi_sweep_ambient.csv", [r for Ta in SENS_AMBIENTS for r in S2[Ta]], ["T_amb"] + COLS)
R["S2_summary"] = {str(Ta): dict(P_port=S2[Ta][0]["P_port"], phi_sat=S2[Ta][0]["phi_sat"],
                                 Q0=S2[Ta][0]["Q_h"], Q50=S2[Ta][25]["Q_h"], Q1=S2[Ta][-1]["Q_h"],
                                 COP0=S2[Ta][0]["COP_h"], COP50=S2[Ta][25]["COP_h"], COP1=S2[Ta][-1]["COP_h"],
                                 T0=S2[Ta][0]["T_dis"], T50=S2[Ta][25]["T_dis"], T1=S2[Ta][-1]["T_dis"],
                                 load=bc.cabin_load(Ta), choked=S2[Ta][-1]["choked"]) for Ta in SENS_AMBIENTS}

# ======================================================================================
# 3. Compressor speed: 4000 rpm vs the speed at which the fully open EEV just meets 5 kW
# ======================================================================================
N5 = brentq(lambda n: fe.ftvi_eev(phi=1.0, N=n)["Q_h"] - BRIEF["Q_load"], BRIEF["N_min"], BRIEF["N_max"])
R["N_5kW_full_open"] = N5
SPEEDS = [BRIEF["N_min"], round(N5), BRIEF["N"]]
S3 = {n: [fe.ftvi_eev(phi=f, N=n) for f in PHI] for n in SPEEDS}
save_csv("T3_phi_sweep_speed.csv", [r for n in SPEEDS for r in S3[n]], COLS)
R["S3_summary"] = {}
for n in SPEEDS:
    feas = [r["phi"] for r in S3[n] if r["Q_h"] >= BRIEF["Q_load"]]
    R["S3_summary"][str(n)] = dict(phi_sat=S3[n][0]["phi_sat"], phi_min_5kW=(min(feas) if feas else None),
                                   Q1=S3[n][-1]["Q_h"], COP1=S3[n][-1]["COP_h"], T1=S3[n][-1]["T_dis"],
                                   xi1=S3[n][-1]["xi"])

# ======================================================================================
# 4. Saturation opening phi_sat vs P_FT (and ambient)
# ======================================================================================
PF = np.round(np.arange(3.5, 6.0001, 0.05), 3)
S4 = []
for Ta in SENS_AMBIENTS:
    for p in PF:
        try:
            r = fe.ftvi_eev(phi=1.0, P_FT=p, T_evap=bc.evap_temperature(Ta))
        except ValueError:
            continue
        S4.append(dict(T_amb=Ta, P_FT=p, P_port=r["P_port"], phi_sat=r["phi_sat"], xi_max=r["xi_max"],
                       xi_full_open=r["xi"], m_valve_full=r["m_valve"], m_inj_max=r["xi_max"] * r["m_suc"],
                       choked=r["choked"]))
save_csv("T3_phi_sat.csv", S4)

# ======================================================================================
# 5. Map over (P_FT, phi) at the design ambient and speed -- input for Student 5
# ======================================================================================
PFm = np.round(np.arange(3.5, 6.0001, 0.1), 3); PHm = np.round(np.linspace(0, 1, 21), 3)
M = {k: np.full((len(PHm), len(PFm)), np.nan) for k in ("COP_h", "Q_h", "T_dis", "xi")}
rows = []
for j, p in enumerate(PFm):
    for i, f in enumerate(PHm):
        r = fe.ftvi_eev(phi=f, P_FT=p)
        for k in M: M[k][i, j] = r[k]
        rows.append(r)
save_csv("T3_map_PFT_phi.csv", rows, COLS)
iopt = np.unravel_index(np.nanargmax(M["COP_h"]), M["COP_h"].shape)
R["map_best_COP"] = dict(P_FT=float(PFm[iopt[1]]), phi=float(PHm[iopt[0]]), COP=float(M["COP_h"][iopt]),
                         Q=float(M["Q_h"][iopt]), T=float(M["T_dis"][iopt]))
iq = np.unravel_index(np.nanargmax(M["Q_h"]), M["Q_h"].shape)
R["map_best_Q"] = dict(P_FT=float(PFm[iq[1]]), phi=float(PHm[iq[0]]), COP=float(M["COP_h"][iq]),
                       Q=float(M["Q_h"][iq]), T=float(M["T_dis"][iq]))

# ======================================================================================
# 6. Sensitivity to the two main Task-3 assumptions
# ======================================================================================
S6 = []
sens_cases = [("CdA x 0.8", dict(CdA_max=0.8 * CdA)), ("CdA nominal", dict(CdA_max=CdA)),
              ("CdA x 1.2", dict(CdA_max=1.2 * CdA)),
              ("r_port 1.4", dict(r_port=1.4)), ("r_port 1.5", dict(r_port=1.5)), ("r_port 1.6", dict(r_port=1.6))]
S6d = {}
for name, kw in sens_cases:
    S6d[name] = [fe.ftvi_eev(phi=f, **kw) for f in PHI]
    for r in S6d[name]:
        r2 = dict(r); r2["case"] = name; S6.append(r2)
save_csv("T3_sensitivity.csv", S6, ["case"] + COLS)
R["S6_summary"] = {n: dict(phi_sat=S6d[n][0]["phi_sat"], Q50=S6d[n][25]["Q_h"], COP50=S6d[n][25]["COP_h"],
                           T50=S6d[n][25]["T_dis"], Q1=S6d[n][-1]["Q_h"], COP1=S6d[n][-1]["COP_h"],
                           T1=S6d[n][-1]["T_dis"], P_port=S6d[n][0]["P_port"]) for n in S6d}

# ======================================================================================
# 7. Hand-over for Student 5: function spec + test cases
# ======================================================================================
tests = []
for (Ta, p, f, n) in [(-10, 4.5, 0.5, 4000), (-10, 4.5, 1.0, 4000), (-10, 3.5, 0.5, 4000), (-10, 6.0, 0.3, 4000),
                      (-5, 4.5, 0.5, 4000), (-15, 6.0, 1.0, 4000), (-10, 4.5, 0.5, 2500), (-10, 5.0, 0.8, 8000)]:
    r = fe.ftvi_eev(phi=f, P_FT=p, T_evap=bc.evap_temperature(Ta), N=n)
    tests.append(dict(T_amb=Ta, P_FT=p, phi=f, N=n, P_port=r["P_port"], m_suc=r["m_suc"], x=r["x"],
                      m_inj=r["m_inj"], xi=r["xi"], limited_by=r["limited_by"], choked=r["choked"],
                      Q_h=r["Q_h"], W=r["W"], COP_h=r["COP_h"], T_dis=r["T_dis"]))
save_csv("T3_handover_testcases.csv", tests)
handover = dict(
    function="eev.m_inj(phi, P_FT, P_port, m_suc, x, CdA_max)  -> dict(m_inj [kg/s], xi, limited_by, choked, ...)",
    cycle="ftvi_eev.ftvi_eev(phi, P_FT, T_evap, N, P_dis, ...) -> full cycle dict (Q_h, W, COP_h, T_dis, xi, ...)",
    units="P [MPa], T [degC], h [kJ/kg], m [kg/s], phi 0..1, areas [m2]",
    constants=dict(CdA_max_m2=CdA, Cd=eev.EEV["Cd"], A_max_m2=CdA / eev.EEV["Cd"], r_port=eev.EEV["r_port"],
                   k=eev.EEV["k"], R_crit=eev.R_CRIT, V_disp_m3=V),
    rules=["P_port = r_port * P_evap", "no injection if P_FT <= P_port",
           "m_inj <= x/(1-x) * m_suc (no liquid injection)",
           "un-injected flash vapour leaves with the liquid to EEV-2: h8 = (m_tot h5 - m_inj h7)/m_suc"],
    decision_variable_ranges=dict(P_FT_MPa=[3.5, 6.0], phi=[0.0, 1.0], N_rpm=[BRIEF["N_min"], BRIEF["N_max"]]),
    test_cases=tests,
)
with open(os.path.join(RES, "T3_handover_student5.json"), "w") as f:
    json.dump(handover, f, indent=2, default=float)

# ======================================================================================
# FIGURES
# ======================================================================================
pft_col = dict(zip(PFT_SET, BLUE_ORD))

# --- Fig 1: EEV flow characteristic --------------------------------------------------
fig, ax = plt.subplots(figsize=(6.4, 3.8))
for p in PFT_SET:
    xs = [r["phi"] * 100 for r in S1[p]]
    ax.plot(xs, [r["m_inj"] * 1e3 for r in S1[p]], color=pft_col[p], lw=1.8, zorder=3)
    ax.plot(xs, [min(r["m_valve"], 0.03) * 1e3 for r in S1[p]], color=pft_col[p], lw=0.9, ls=(0, (2, 2)), zorder=2)
    end_label(ax, 100, S1[p][-1]["m_inj"] * 1e3, f"{p:.1f} MPa", INK2, dy={4.0: -6, 5.0: 6}.get(p, 0))
    ps = S1[p][0]["phi_sat"]
    if ps <= 1:
        ax.plot(ps * 100, S1[p][0]["xi_max"] * S1[p][0]["m_suc"] * 1e3, "o", ms=5, color=pft_col[p],
                mec=SURF, mew=1.2, zorder=4)
ax.set_xlim(0, 100); ax.set_ylim(0, 15)
ax.set_xlabel("EEV opening, φ (%)"); ax.set_ylabel("Injection flow, ṁ$_{inj}$ (g/s)")
ax.set_title("Injection flow vs EEV opening, by flash-tank pressure P$_{FT}$", loc="left")
ax.text(0.01, 0.98, "Solid: actual injection   Dotted: what the valve alone could pass   ● tank runs out of vapour (φ$_{sat}$)\n"
        "−10 °C ambient · 4000 rpm · P$_{port}$ = %.2f MPa" % S1[4.5][0]["P_port"],
        transform=ax.transAxes, va="top", fontsize=7.5, color=INK2)
fig_save(fig, "fig_T3_1_eev_flow.png")

# --- Fig 2: performance vs opening (4 panels) ----------------------------------------
fig, axs = plt.subplots(2, 2, figsize=(7.6, 5.6), sharex=True)
panels = [("xi", "Injection ratio, ξ = ṁ$_{inj}$/ṁ$_{suc}$", 1), ("Q_h", "Heating capacity, Q$_H$ (kW)", 1),
          ("COP_h", "Heating COP", 1), ("T_dis", "Discharge temperature, T$_{dis}$ (°C)", 1)]
for ax, (k, lab, _) in zip(axs.flat, panels):
    for p in PFT_SET:
        ax.plot([r["phi"] * 100 for r in S1[p]], [r[k] for r in S1[p]], color=pft_col[p])
    ax.set_title(lab, loc="left"); ax.set_xlim(0, 100)
axs[0, 1].axhline(BRIEF["Q_load"], color=MUTED, lw=1, ls=(0, (4, 3)))
axs[0, 1].text(99, BRIEF["Q_load"], "5 kW load", ha="right", va="bottom", fontsize=7.5, color=INK2)
axs[1, 0].axhline(A["COP_h"], color=MUTED, lw=1, ls=(0, (4, 3)))
axs[1, 0].text(99, A["COP_h"], "baseline (Case A) %.2f" % A["COP_h"], ha="right", va="top", fontsize=7.5, color=INK2)
axs[1, 1].axhline(A["T_dis"], color=MUTED, lw=1, ls=(0, (4, 3)))
axs[1, 1].text(99, A["T_dis"], "baseline %.1f °C" % A["T_dis"], ha="right", va="bottom", fontsize=7.5, color=INK2)
for ax in axs[1]: ax.set_xlabel("EEV opening, φ (%)")
handles = [plt.Line2D([], [], color=pft_col[p], lw=2) for p in PFT_SET]
fig.legend(handles, [f"P$_{{FT}}$ = {p:.1f} MPa" for p in PFT_SET], loc="upper center", ncol=5,
           bbox_to_anchor=(0.5, 1.02), fontsize=8)
fig.text(0.5, -0.01, "−10 °C ambient · 4000 rpm · P$_{dis}$ = 9.5 MPa · injection port at %.2f MPa" % S1[4.5][0]["P_port"],
         ha="center", fontsize=7.5, color=INK2)
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig_save(fig, "fig_T3_2_performance_vs_opening.png")

# --- Fig 3: ambient sensitivity -------------------------------------------------------
amb_col = {-5.0: CAT[0], -10.0: CAT[1], -15.0: CAT[2]}
fig, axs = plt.subplots(1, 3, figsize=(8.6, 3.0))
for ax, (k, lab) in zip(axs, [("Q_h", "Q$_H$ (kW)"), ("COP_h", "Heating COP"), ("T_dis", "T$_{dis}$ (°C)")]):
    for Ta in SENS_AMBIENTS:
        ax.plot([r["phi"] * 100 for r in S2[Ta]], [r[k] for r in S2[Ta]], color=amb_col[Ta])
    ax.set_title(lab, loc="left"); ax.set_xlabel("EEV opening, φ (%)"); ax.set_xlim(0, 100)
axs[0].axhline(BRIEF["Q_load"], color=MUTED, lw=1, ls=(0, (4, 3)))
axs[0].text(99, BRIEF["Q_load"], "5 kW", ha="right", va="bottom", fontsize=7.5, color=INK2)
axs[2].axhline(BRIEF["T_dis_max"], color=MUTED, lw=1, ls=(0, (4, 3)))
axs[2].text(99, BRIEF["T_dis_max"], "130 °C limit", ha="right", va="top", fontsize=7.5, color=INK2)
handles = [plt.Line2D([], [], color=amb_col[t], lw=2) for t in SENS_AMBIENTS]
fig.legend(handles, [f"{t:.0f} °C ambient" for t in SENS_AMBIENTS], loc="upper center", ncol=3,
           bbox_to_anchor=(0.5, 1.06), fontsize=8)
fig.tight_layout(rect=(0, 0, 1, 0.93))
fig_save(fig, "fig_T3_3_ambient.png")

# --- Fig 4: compressor speed ----------------------------------------------------------
sp_col = {SPEEDS[0]: CAT[0], SPEEDS[1]: CAT[1], SPEEDS[2]: CAT[2]}
fig, axs = plt.subplots(1, 3, figsize=(8.6, 3.0))
for ax, (k, lab) in zip(axs, [("xi", "Injection ratio ξ"), ("Q_h", "Q$_H$ (kW)"), ("T_dis", "T$_{dis}$ (°C)")]):
    for n in SPEEDS:
        ax.plot([r["phi"] * 100 for r in S3[n]], [r[k] for r in S3[n]], color=sp_col[n])
    ax.set_title(lab, loc="left"); ax.set_xlabel("EEV opening, φ (%)"); ax.set_xlim(0, 100)
axs[1].axhline(BRIEF["Q_load"], color=MUTED, lw=1, ls=(0, (4, 3)))
axs[1].text(99, BRIEF["Q_load"] + 0.08, "5 kW load", ha="right", va="bottom", fontsize=7.5, color=INK2)
handles = [plt.Line2D([], [], color=sp_col[n], lw=2) for n in SPEEDS]
fig.legend(handles, [f"{n:.0f} rpm" for n in SPEEDS], loc="upper center", ncol=3, bbox_to_anchor=(0.5, 1.06), fontsize=8)
fig.tight_layout(rect=(0, 0, 1, 0.93))
fig_save(fig, "fig_T3_4_speed.png")

# --- Fig 5: saturation opening --------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(8.0, 3.2))
for Ta in SENS_AMBIENTS:
    rr = [r for r in S4 if r["T_amb"] == Ta]
    axs[0].plot([r["P_FT"] for r in rr], [min(r["phi_sat"], 2.0) * 100 for r in rr], color=amb_col[Ta])
    axs[1].plot([r["P_FT"] for r in rr], [r["xi_full_open"] for r in rr], color=amb_col[Ta])
    axs[1].plot([r["P_FT"] for r in rr], [r["xi_max"] for r in rr], color=amb_col[Ta], lw=0.9, ls=(0, (2, 2)))
axs[0].axhline(100, color=MUTED, lw=1, ls=(0, (4, 3)))
axs[0].text(5.98, 100, "fully open", ha="right", va="bottom", fontsize=7.5, color=INK2)
axs[0].set_ylim(0, 200); axs[0].set_xlim(3.5, 6.0)
axs[0].set_title("Saturation opening φ$_{sat}$ (%)", loc="left")
axs[0].text(0.02, 0.03, "Above 100 %: the valve can never empty the tank (axis capped at 200 %)", transform=axs[0].transAxes,
            fontsize=7.5, color=INK2)
axs[1].set_title("Injection ratio: fully open (solid), vapour available (dotted, same at all ambients)", loc="left", fontsize=8.5)
axs[1].set_xlim(3.5, 6.0)
for ax in axs: ax.set_xlabel("Flash-tank pressure, P$_{FT}$ (MPa)")
handles = [plt.Line2D([], [], color=amb_col[t], lw=2) for t in SENS_AMBIENTS]
fig.legend(handles, [f"{t:.0f} °C ambient" for t in SENS_AMBIENTS], loc="upper center", ncol=3,
           bbox_to_anchor=(0.5, 1.06), fontsize=8)
fig.tight_layout(rect=(0, 0, 1, 0.93))
fig_save(fig, "fig_T3_5_phi_sat.png")

# --- Fig 6: maps ----------------------------------------------------------------------
fig, axs = plt.subplots(1, 3, figsize=(9.6, 3.2), sharey=True)
cmap = matplotlib.colors.LinearSegmentedColormap.from_list("b", ["#cde2fb", "#6da7ec", "#2a78d6", "#104281", "#0d366b"])
for ax, (k, lab, fmt) in zip(axs, [("COP_h", "Heating COP", "%.2f"), ("Q_h", "Q$_H$ (kW)", "%.1f"),
                                   ("T_dis", "T$_{dis}$ (°C)", "%.0f")]):
    Zm = M[k]
    cm = cmap if k != "T_dis" else cmap.reversed()
    cf = ax.contourf(PFm, PHm * 100, Zm, levels=12, cmap=cm)
    cs = ax.contour(PFm, PHm * 100, Zm, levels=6, colors=SURF, linewidths=0.6)
    ax.clabel(cs, fmt=fmt, fontsize=7, colors="white" if True else INK)
    if k == "Q_h":
        c5 = ax.contour(PFm, PHm * 100, Zm, levels=[BRIEF["Q_load"]], colors=INK, linewidths=1.2, linestyles="--")
        ax.clabel(c5, fmt="5 kW", fontsize=7.5, colors=INK)
    ax.set_title(lab, loc="left"); ax.set_xlabel("P$_{FT}$ (MPa)"); ax.grid(False)
    ax.plot(4.5, 50, marker="o", ms=6, color=CAT[1], mec="white", mew=1.2)
axs[0].set_ylabel("EEV opening, φ (%)")
axs[0].plot(R["map_best_COP"]["P_FT"], R["map_best_COP"]["phi"] * 100, marker="*", ms=10, color=CAT[3], mec="white", mew=0.8)
fig.text(0.5, -0.03, "● brief's design point (4.5 MPa, 50 %)   ★ highest COP on this map · −10 °C ambient · 4000 rpm",
         ha="center", fontsize=7.5, color=INK2)
fig.tight_layout()
fig_save(fig, "fig_T3_6_maps.png")

# --- Fig 7: sensitivity ---------------------------------------------------------------
fig, axs = plt.subplots(1, 2, figsize=(8.0, 3.0), sharey=True)
for ax, names, title in [(axs[0], ["CdA x 0.8", "CdA nominal", "CdA x 1.2"], "Valve size C$_d$A$_{max}$ ±20 %"),
                         (axs[1], ["r_port 1.4", "r_port 1.5", "r_port 1.6"], "Port pressure ratio r$_{port}$")]:
    for c, n in zip([CAT[0], INK2, CAT[1]], names):
        lab = n.replace("CdA", "C$_d$A").replace("r_port", "r$_{port}$ =")
        if n in ("r_port 1.4", "r_port 1.5", "r_port 1.6"):
            lab += f"  (P$_{{port}}$ {S6d[n][0]['P_port']:.2f} MPa)"
        ax.plot([r["phi"] * 100 for r in S6d[n]], [r["Q_h"] for r in S6d[n]], color=c, label=lab,
                lw=1.4 if "nominal" in n or n == "r_port 1.5" else 1.8)
    ax.set_title(title, loc="left"); ax.set_xlabel("EEV opening, φ (%)"); ax.set_xlim(0, 100)
    ax.legend(fontsize=7.5, loc="lower right")
axs[0].set_ylabel("Q$_H$ (kW)")
fig.tight_layout()
fig_save(fig, "fig_T3_7_sensitivity.png")

# --- Fig 8: P-h diagram for three openings --------------------------------------------
fig, ax = plt.subplots(figsize=(7.4, 4.4))
Ts = np.linspace(-25, 30.97, 300)
dome_h = [props.hf(props.Psat(t)) for t in Ts] + [props.hg(props.Psat(t)) for t in Ts[::-1]]
dome_P = [props.Psat(t) for t in Ts] + [props.Psat(t) for t in Ts[::-1]]
ax.plot(dome_h, dome_P, color=AXIS, lw=1.2)
ax.text(props.hf(props.Psat(-5)) - 5, props.Psat(-5), "sat. liquid", ha="right", fontsize=7, color=MUTED)
ax.text(props.hg(props.Psat(-8)) - 5, props.Psat(-8), "sat. vapour", ha="right", fontsize=7, color=MUTED)
for f, c in [(0.0, CAT[0]), (0.5, CAT[1]), (1.0, CAT[2])]:
    r = fe.ftvi_eev(phi=f)
    Pe, Pp, Pf, Pd = r["P_evap"], r["P_port"], r["P_FT"], r["P_dis"]
    # compression path: 1 -> 2 -> 3 -> 4 ; high side 4 -> 5 ; EEV-1 5 -> 6 ; tank 6 -> 8 ; EEV-2 8 -> 9 ; evap 9 -> 1
    path_h = [r["h1"], r["h2"], r["h3"], r["h4"], r["h5"], r["h6"], r["h8"], r["h9"], r["h1"]]
    path_P = [Pe, Pp, Pp, Pd, Pd, Pf, Pf, Pe, Pe]
    ax.plot(path_h, path_P, color=c, lw=1.6, marker="o", ms=3, label=f"φ = {f*100:.0f} %  (COP {r['COP_h']:.2f}, T$_{{dis}}$ {r['T_dis']:.0f} °C)")
    if f > 0:
        ax.plot([r["h7"], r["h7a"], r["h3"]], [Pf, Pp, Pp], color=c, lw=1.0, ls=(0, (2, 2)))
ax.set_yscale("log"); ax.set_ylim(1.8, 11.5); ax.set_xlabel("Specific enthalpy, h (kJ/kg)"); ax.set_ylabel("Pressure, P (MPa)")
ax.set_yticks([2, 3, 4, 5, 6, 8, 10]); ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_title("P–h diagram at three EEV openings (P$_{FT}$ = 4.5 MPa, −10 °C, 4000 rpm)", loc="left")
ax.legend(fontsize=7.5, loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3, columnspacing=1.2, handlelength=1.5)
ax.text(0.5, -0.30, "Dotted: injected vapour 7 → 7a (injection EEV) → mixing at P$_{port}$",
        transform=ax.transAxes, ha="center", fontsize=7.5, color=INK2)
fig_save(fig, "fig_T3_8_ph.png")

with open(os.path.join(RES, "T3_summary.json"), "w") as f:
    json.dump(R, f, indent=2, default=lambda o: float(o) if isinstance(o, (np.floating, np.integer)) else str(o))

# ------------------------------- console summary ----------------------------------------
print("Backend:", props.name)
print("Cd*A_max = %.3f mm2 -> A_max %.3f mm2 (d %.2f mm)" % (R["CdA_mm2"], R["A_max_mm2"], R["d_eq_mm"]))
print("Validation vs Student 2 (ideal-port limit): max rel. diff = %.2e" % R["validation_max_rel_diff"])
print("Design point (phi 50%%, 4.5 MPa): xi=%.3f Q=%.2f kW W=%.2f kW COP=%.3f Tdis=%.1f C, limited by %s" % (
    D["xi"], D["Q_h"], D["W"], D["COP_h"], D["T_dis"], D["limited_by"]))
for p, s in R["S1_summary"].items(): print(" P_FT", p, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in s.items()})
for t, s in R["S2_summary"].items(): print(" T_amb", t, {k: (round(v, 3) if isinstance(v, float) else v) for k, v in s.items()})
print(" N for 5 kW fully open: %.0f rpm" % N5)
for n, s in R["S3_summary"].items(): print(" N", n, s)
print(" Map best COP", R["map_best_COP"]); print(" Map best Q", R["map_best_Q"])
for n, s in R["S6_summary"].items(): print(" ", n, {k: round(v, 3) for k, v in s.items()})
