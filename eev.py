"""
eev.py  --  STUDENT 3 : injection EEV and flow control
=======================================================
The injection EEV sits on the vapour line between the flash tank (state 7, sat. vapour at P_FT)
and the compressor injection port (pressure P_port).  Units as in co2_props:
P [MPa], T [degC], h [kJ/kg], rho [kg/m3], m [kg/s], areas [m2].

Model (stated assumptions, see report):
  1. Port pressure      P_port = r_port * P_evap          (fixed port position in the compression)
  2. Valve flow         m_valve = Cd * phi * A_max * Y * sqrt(2 * rho_v * dP)
                        Y  = 1 - 0.41 * dP / (k * P_FT)   (compressible expansion factor, thin orifice)
                        dP = min(P_FT - P_port, (1 - r_crit) * P_FT)   (choked-flow limit)
  3. Vapour cap         m_inj = min(m_valve, m_inj_max),  m_inj_max = x/(1-x) * m_suc
                        (the tank cannot supply more vapour than flashes in EEV-1; no liquid injection)
  4. Sizing rule        fully open valve passes SIZE_MARGIN x the available vapour at the design point
"""
import math
from config import BRIEF, ASSUME
from co2_props import get_props

props = get_props()

# ----------------------------- Task 3 assumptions (stated in report) --------------------
EEV = dict(
    r_port=1.5,         # -    port pressure / suction pressure (fixed port location)
    Cd=0.70,            # -    discharge coefficient of the EEV orifice
    k=1.30,             # -    isentropic exponent of CO2 vapour (choke + expansion factor)
    size_margin=1.20,   # -    fully-open valve passes 1.2 x available vapour at design point
    design=dict(T_amb=-10.0, P_FT=BRIEF["P_inj"], N=BRIEF["N"]),
)
R_CRIT = (2.0 / (EEV["k"] + 1.0)) ** (EEV["k"] / (EEV["k"] - 1.0))   # critical pressure ratio (~0.546)


def port_pressure(P_evap, r_port=None):
    """Compressor injection-port pressure [MPa]."""
    return (EEV["r_port"] if r_port is None else r_port) * P_evap


def valve_flow(phi, P_FT, P_port, CdA_max, k=None):
    """Mass flow [kg/s] the EEV can pass at opening phi (0..1), ignoring the vapour cap.
    Returns (m_valve, dP_eff [MPa], choked [bool])."""
    k = EEV["k"] if k is None else k
    phi = min(max(phi, 0.0), 1.0)
    dP = P_FT - P_port
    if dP <= 0.0 or phi == 0.0:
        return 0.0, max(dP, 0.0), False
    dP_crit = (1.0 - R_CRIT) * P_FT
    choked = dP > dP_crit
    dP_eff = min(dP, dP_crit)
    Y = 1.0 - 0.41 * dP_eff / (k * P_FT)
    rho_v = props.rhog(P_FT)
    m = CdA_max * phi * Y * math.sqrt(2.0 * rho_v * dP_eff * 1e6)
    return m, dP_eff, choked


def flash_quality(P_FT, h5):
    hf, hg = props.hf(P_FT), props.hg(P_FT)
    return (h5 - hf) / (hg - hf)


def m_inj(phi, P_FT, P_port, m_suc, x, CdA_max):
    """HAND-OVER FUNCTION for Student 5.
    phi      EEV opening 0..1          P_FT    flash-tank pressure [MPa]
    P_port   injection-port pressure [MPa]   m_suc  suction (= tank-liquid) flow [kg/s]
    x        quality entering the flash tank  CdA_max  Cd*A_max of the valve [m2]
    Returns dict: m_inj [kg/s], xi = m_inj/m_suc, m_valve, m_inj_max, limited_by ('valve'|'vapour'|'none'),
    choked flag."""
    m_valve, dP, choked = valve_flow(phi, P_FT, P_port, CdA_max)
    m_max = x / (1.0 - x) * m_suc
    m = min(m_valve, m_max)
    if m <= 0.0:
        lim = "none"
    else:
        lim = "vapour" if m_valve >= m_max else "valve"
    return dict(m_inj=m, xi=m / m_suc, m_valve=m_valve, m_inj_max=m_max,
                limited_by=lim, choked=choked, dP=dP)


def phi_saturation(P_FT, P_port, m_suc, x, CdA_max):
    """Smallest opening at which the tank's vapour is exhausted (opening further does nothing).
    Returns >1 (i.e. not reached) if the fully open valve cannot pass all the vapour."""
    m_full, _, _ = valve_flow(1.0, P_FT, P_port, CdA_max)
    if m_full <= 0.0:
        return float("inf")
    return (x / (1.0 - x) * m_suc) / m_full
