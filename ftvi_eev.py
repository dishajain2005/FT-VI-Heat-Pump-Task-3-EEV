"""
ftvi_eev.py  --  STUDENT 3 : FT-VI cycle with injection-EEV flow control
=========================================================================
Extends Student 2's ftvi() (flash_tank_vi.py) without changing it.

State points (Student 2's numbering + 7a)
  1   suction                              (P_evap, T_evap + SH)
  2   end of stage 1                       (P_port)
  7   flash-tank vapour, sat. vapour       (P_FT)
  7a  after injection EEV (isenthalpic)    (P_port, h7a = h7)   superheated vapour
  3   after mixing                         (P_port, h3 = (m_suc h2 + m_inj h7)/m_tot)
  4   discharge                            (P_dis)
  5   gas-cooler outlet                    (P_dis, T_gc_out)
  6   EEV-1 outlet = tank inlet            (P_FT, h6 = h5), quality x
  8   tank outlet to EEV-2                 (P_FT, h8 = (m_tot h5 - m_inj h7)/m_suc)
                                           = sat. liquid only when ALL flash vapour is injected
  9   EEV-2 outlet = evaporator inlet      (P_evap, h9 = h8)

Flows: m_suc fixed by the compressor (swept volume); m_inj from eev.m_inj(); m_tot = m_suc + m_inj.
"""
from config import BRIEF, ASSUME
from co2_props import get_props
from baseline_cycle import vol_eff, calibrate_displacement, evap_temperature
import eev

props = get_props()


# ----------------------------------------------------------------------------------------
# Compressor (swappable: Student 4 can replace this one function with the same signature)
# ----------------------------------------------------------------------------------------
def compressor(P_evap, h1, s1, P_port, P_dis, m_suc, m_inj, h_inj, eta_is=BRIEF["eta_is"]):
    """Two-stage compression with mid-compression vapour injection at P_port.
    Returns dict with h2, T2, h3, T3, s3, h4, T4 (= T_dis), W [kW], W1, W2."""
    h2 = h1 + (props.h_ps(P_port, s1) - h1) / eta_is
    m_tot = m_suc + m_inj
    h3 = (m_suc * h2 + m_inj * h_inj) / m_tot
    T3 = props.T_ph(P_port, h3)
    s3 = props.s(P_port, T3)
    h4 = h3 + (props.h_ps(P_dis, s3) - h3) / eta_is
    T4 = props.T_ph(P_dis, h4)
    W1, W2 = m_suc * (h2 - h1), m_tot * (h4 - h3)
    return dict(h2=h2, T2=props.T_ph(P_port, h2), h3=h3, T3=T3, s3=s3, h4=h4, T4=T4,
                W1=W1, W2=W2, W=W1 + W2)


# ----------------------------------------------------------------------------------------
# Valve sizing at the design point
# ----------------------------------------------------------------------------------------
_V = {}


def V_disp():
    if "V" not in _V:
        _V["V"] = calibrate_displacement()       # Student 1's calibration (Case A = 5.0 kW)
    return _V["V"]


def size_valve(margin=None, r_port=None, P_dis=BRIEF["P_dis"], T_gc_out=BRIEF["T_gc_out"]):
    """Cd*A_max [m2] such that phi = 1 passes `margin` x the available vapour at the design point."""
    margin = eev.EEV["size_margin"] if margin is None else margin
    d = eev.EEV["design"]
    T_ev = evap_temperature(d["T_amb"])
    P_e = props.Psat(T_ev)
    T1 = T_ev + ASSUME["superheat"]
    rho1 = props.rho(P_e, T1)
    m_suc = rho1 * V_disp() * d["N"] / 60.0 * vol_eff(P_dis, P_e)
    x = eev.flash_quality(d["P_FT"], props.h(P_dis, T_gc_out))
    m_max = x / (1 - x) * m_suc
    m_unit, _, _ = eev.valve_flow(1.0, d["P_FT"], eev.port_pressure(P_e, r_port), 1.0)
    return margin * m_max / m_unit


_CDA = {}


def CdA_design():
    if "v" not in _CDA:
        _CDA["v"] = size_valve()
    return _CDA["v"]


# ----------------------------------------------------------------------------------------
# Cycle
# ----------------------------------------------------------------------------------------
def ftvi_eev(phi=BRIEF["eev_open"] / 100.0, P_FT=BRIEF["P_inj"], T_evap=BRIEF["T_evap"],
             N=BRIEF["N"], P_dis=BRIEF["P_dis"], T_gc_out=BRIEF["T_gc_out"],
             eta_is=BRIEF["eta_is"], SH=ASSUME["superheat"], CdA_max=None, r_port=None,
             P_port=None):
    """Steady-state FT-VI cycle with the injection EEV at opening phi (0..1).
    P_port can be forced (e.g. P_port -> P_FT reproduces Student 2's ideal-port model)."""
    CdA_max = CdA_design() if CdA_max is None else CdA_max
    P_e = props.Psat(T_evap)
    P_port = eev.port_pressure(P_e, r_port) if P_port is None else P_port
    if not (P_e < P_port < P_dis and P_e < P_FT < P_dis):
        raise ValueError("Constraint violated: need P_evap < P_port, P_FT < P_dis")
    T1 = T_evap + SH
    h1, s1, rho1 = props.h(P_e, T1), props.s(P_e, T1), props.rho(P_e, T1)
    eta_v = vol_eff(P_dis, P_e)
    m_suc = rho1 * V_disp() * N / 60.0 * eta_v

    h5 = props.h(P_dis, T_gc_out)
    x = eev.flash_quality(P_FT, h5)
    h7 = props.hg(P_FT)
    inj = eev.m_inj(phi, P_FT, P_port, m_suc, x, CdA_max)
    m_i = inj["m_inj"]
    m_tot = m_suc + m_i

    c = compressor(P_e, h1, s1, P_port, P_dis, m_suc, m_i, h7, eta_is)

    h8 = (m_tot * h5 - m_i * h7) / m_suc
    hf_e, hg_e = props.hf(P_e), props.hg(P_e)
    x9 = (h8 - hf_e) / (hg_e - hf_e)
    x8 = (h8 - props.hf(P_FT)) / (h7 - props.hf(P_FT))

    Q_h = m_tot * (c["h4"] - h5)
    Q_e = m_suc * (h1 - h8)
    W = c["W"]
    tank_res = m_tot * h5 - (m_i * h7 + m_suc * h8)          # adiabatic tank energy residual
    return dict(
        phi=phi, P_FT=P_FT, P_port=P_port, T_evap=T_evap, P_evap=P_e, P_dis=P_dis, N=N,
        T_FT=props.Tsat(P_FT), T_port_sat=props.Tsat(P_port) if P_port < 7.37 else float("nan"),
        eta_v=eta_v, x=x,
        h1=h1, T1=T1, h2=c["h2"], T2=c["T2"], h3=c["h3"], T3=c["T3"], h4=c["h4"], T4=c["T4"],
        h5=h5, h6=h5, h7=h7, h7a=h7, T7a=props.T_ph(P_port, h7), h8=h8, x8=x8, h9=h8, x9=x9,
        m_suc=m_suc, m_inj=m_i, m_tot=m_tot, xi=inj["xi"], xi_max=inj["m_inj_max"] / m_suc,
        m_valve=inj["m_valve"], limited_by=inj["limited_by"], choked=inj["choked"], dP_eev=inj["dP"],
        phi_sat=eev.phi_saturation(P_FT, P_port, m_suc, x, CdA_max),
        Q_h=Q_h, Q_e=Q_e, W=W, W1=c["W1"], W2=c["W2"], COP_h=Q_h / W, T_dis=c["T4"],
        mass_err=m_tot - (m_suc + m_i), tank_energy_residual=tank_res,
        energy_balance_err=Q_h - Q_e - W,
        feasible=(c["T4"] <= BRIEF["T_dis_max"]) and (P_dis <= BRIEF["P_max"]) and (Q_h >= BRIEF["Q_load"]),
    )


if __name__ == "__main__":
    print("Backend:", props.name)
    print("Cd*A_max = %.4f mm2  (A_max = %.3f mm2, d = %.2f mm at Cd = %.2f)" % (
        CdA_design() * 1e6, CdA_design() / eev.EEV["Cd"] * 1e6,
        (4 * CdA_design() / eev.EEV["Cd"] / 3.14159265) ** 0.5 * 1e3, eev.EEV["Cd"]))
    for p in (0.0, 0.25, 0.5, 0.75, 1.0):
        r = ftvi_eev(phi=p)
        print("phi=%.2f xi=%.3f Q=%.3f W=%.3f COP=%.3f Tdis=%.1f lim=%s  EB=%.1e" % (
            p, r["xi"], r["Q_h"], r["W"], r["COP_h"], r["T_dis"], r["limited_by"], r["energy_balance_err"]))
