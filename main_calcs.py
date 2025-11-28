"""
main_calcs.py

Compute choked mass flow rate through an intake restrictor using:
    m_dot = C_d * A_t * P_t * sqrt(gamma / (R * T_t)) *
            ((gamma + 1) / 2) ** (-(gamma + 1) / (2 * (gamma - 1)))

Assumptions:
- Gas is dry air, modeled as ideal gas.
- Steady, 1-D, adiabatic, isentropic flow up to the throat.
- Flow is choked at the throat (Mach = 1 at minimum area) - we can go back and look at this
- Upstream “reservoir” is the airbox / ambient just before the restrictor:
    -> Flow speed is low, so stagnation ≈ static - we can revisit this as well
    -> P_t ≈ P_ambient,  T_t ≈ T_ambient.
- Discharge coefficient C_d can be modified as needed
"""

import math
from CoolProp.CoolProp import PropsSI

# Import your unit conversion helpers
from unit_conversions import *


def compute_restrictor_choked_mdot(
    d_mm: float,
    Cd: float,
    P_amb_kpa: float,
    T_amb_C: float,
    fluid: str = "Air",
) -> dict:
    """
    Compute choked mass flow rate through a circular restrictor.

    Inputs:
        d_mm       : Restrictor throat diameter [mm]
        Cd         : Discharge coefficient [-]
        P_amb_kpa  : Upstream ambient (reservoir) pressure [kPa]
        T_amb_C    : Upstream ambient (reservoir) temperature [°C]
        fluid      : Working fluid name for CoolProp (default: 'Air')

    Returns:
        dict with:
            - d_m          : diameter [m]
            - A_t          : throat area [m^2]
            - P_t          : upstream total pressure [Pa]
            - T_t          : upstream total temperature [K]
            - gamma        : cp/cv [-]
            - R            : specific gas constant [J/(kg*K)]
            - K_gamma      : isentropic factor [-]
            - prop_factor  : sqrt(gamma / (R*T_t)) [- / sqrt(K)]
            - m_dot_kg_s   : choked mass flow [kg/s]
            - m_dot_lb_hr  : choked mass flow [lb/hr]
    """

    # ------------------------------------------------------------------
    # 1) Convert user inputs to SI units
    #    Assumption: upstream reservoir is airbox / ambient before restrictor,
    #    so P_t ≈ P_ambient, T_t ≈ T_ambient.
    # ------------------------------------------------------------------
    d_m = mm_m(d_mm)                 # mm -> m
    P_t = kpa_pa(P_amb_kpa)          # kPa -> Pa
    T_t = c_k(T_amb_C)               # °C -> K

    # ------------------------------------------------------------------
    # 2) Geometry: throat area
    #    Equation: A_t = pi * d^2 / 4
    # ------------------------------------------------------------------
    A_t = math.pi * (d_m ** 2) / 4.0

    # ------------------------------------------------------------------
    # 3) Get gas properties from CoolProp at upstream state
    #    - cp = Cpmass(P_t, T_t)
    #    - cv = Cvmass(P_t, T_t)
    #    - gamma = cp / cv
    #    - R = cp - cv   (ideal-gas relation)
    # ------------------------------------------------------------------
    cp = PropsSI("Cpmass", "P", P_t, "T", T_t, fluid)
    cv = PropsSI("Cvmass", "P", P_t, "T", T_t, fluid)
    gamma = cp / cv
    R = cp - cv

    # ------------------------------------------------------------------
    # 4) Isentropic factor K(gamma)
    #    From NASA choked-flow equation:
    #    K(gamma) = ((gamma + 1)/2) ^ (-(gamma + 1) / (2 * (gamma - 1)))
    # ------------------------------------------------------------------
    K_gamma = ((gamma + 1.0) / 2.0) ** (-(gamma + 1.0) / (2.0 * (gamma - 1.0)))

    # ------------------------------------------------------------------
    # 5) Property factor: sqrt(gamma / (R * T_t))
    #    From m_dot = C_d * A_t * P_t * sqrt(gamma / (R*T_t)) * K(gamma)
    # ------------------------------------------------------------------
    prop_factor = math.sqrt(gamma / (R * T_t))

    # ------------------------------------------------------------------
    # 6) Choked mass flow rate
    #    Final NASA form:
    #    m_dot = C_d * A_t * P_t * sqrt(gamma / (R*T_t)) * K(gamma)
    # ------------------------------------------------------------------
    m_dot_kg_s = Cd * A_t * P_t * prop_factor * K_gamma

    # Optional: also return in lb/hr for engine world convenience
    m_dot_lb_hr = kgs_lbhr(m_dot_kg_s)

    return {
        "d_m": d_m,
        "A_t": A_t,
        "P_t": P_t,
        "T_t": T_t,
        "gamma": gamma,
        "R": R,
        "K_gamma": K_gamma,
        "prop_factor": prop_factor,
        "m_dot_kg_s": m_dot_kg_s,
        "m_dot_lb_hr": m_dot_lb_hr,
    }


if __name__ == "__main__":
    result = compute_restrictor_choked_mdot(
        d_mm=20.0,
        Cd=0.90,
        P_amb_kpa=101.325,
        T_amb_C=20.0,
    )

    print("=== Restrictor Choked Flow Calculation ===")
    print(f"Diameter (m)          : {result['d_m']:.6f}")
    print(f"Area (m^2)            : {result['A_t']:.6e}")
    print(f"P_t (Pa)              : {result['P_t']:.1f}")
    print(f"T_t (K)               : {result['T_t']:.2f}")
    print(f"gamma                 : {result['gamma']:.5f}")
    print(f"R (J/kg-K)            : {result['R']:.2f}")
    print(f"K(gamma)              : {result['K_gamma']:.6f}")
    print(f"prop_factor           : {result['prop_factor']:.6e}")
    print(f"m_dot (kg/s)          : {result['m_dot_kg_s']:.5f}")
    print(f"m_dot (lb/hr)         : {result['m_dot_lb_hr']:.1f}")
