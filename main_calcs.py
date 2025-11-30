"""
main_calcs.py

Compute choked mass flow rate through an intake restrictor using:
    m_dot = C_d * A_t * P_t * sqrt(gamma / (R * T_t)) *
            ((gamma + 1) / 2) ** (-(gamma + 1) / (2 * (gamma - 1)))

Assumptions:
- Gas is dry air, modeled as ideal gas.
- Steady, 1-D, adiabatic, isentropic flow up to the throat.
- Flow is choked at the throat (Mach = 1 at minimum area) 
    - Note that we are not actually choked at throat
    - Separate calcs exist for true throat velocity and can be recalc'd here for better estimation
- Upstream “reservoir” is the airbox / ambient just before the restrictor:
    -> Flow speed is low, so stagnation ≈ static - we can revisit this as well to account for forward velocity
    -> P_t ≈ P_ambient,  T_t ≈ T_ambient.
- Discharge coefficient C_d can be modified as needed
"""

import math
import numpy as np
import matplotlib.pyplot as plt
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
    #    Assumption: upstream reservoir is ambient before restrictor,
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
    #    - R = cp - cv 
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

    # Optional: also return in lb/hr for engine spec comparison
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

# P-Star Calculation (Critical pressure at throat)

# ----------------------------------------------------------------------
# Critical pressure at the throat (p-star)
#
# Assumptions:
# - Same gas model as main calculation (dry air, ideal gas behavior).
# - Stagnation/total pressure upstream is P_t.
# - Flow is isentropic up to the throat.
#
# Equation (isentropic relation for M = 1 at throat):
#   p_star = P_t / ((gamma + 1)/2) ** (gamma / (gamma - 1))
#
# This is the static pressure at the throat when the flow is just choked.
# ----------------------------------------------------------------------
def compute_p_star(P_t: float, gamma: float) -> float:
    """
    Compute critical (sonic) pressure p* at the throat.

    Inputs:
        P_t   : upstream total (stagnation) pressure [Pa]
        gamma : heat capacity ratio cp/cv [-]

    Returns:
        p_star : critical static pressure at throat [Pa]
    """
    p_star = P_t / (((gamma + 1.0) / 2.0) ** (gamma / (gamma - 1.0)))
    return p_star

# Engine flow rate calculation

# ----------------------------------------------------------------------
# Engine mass flow rate (referenced to ambient conditions)
#
# Assumptions:
# - 4-stroke engine: one intake stroke every 2 revolutions.
# - Volumetric efficiency VE is defined w.r.t. ambient conditions.
# - Ambient state (P_amb_kpa, T_amb_C) is the same upstream
#   ambient state used for the restrictor calculation.
#
# Equations:
#   Vd_m3      = displacement_cc * 1e-6
#   N_int      = RPM / (2 * 60)                [intake events per second]
#   Vdot_eng   = Vd_m3 * VE * N_int            [m^3/s] - from 
#   rho_amb    = rho(P_amb, T_amb)             [CoolProp]
#   m_dot_eng  = Vdot_eng * rho_amb            [kg/s]
# ----------------------------------------------------------------------
def compute_engine_mass_flow_ambient(
    displacement_cc: float,
    rpm: float,
    VE: float,
    P_amb_kpa: float,
    T_amb_C: float,
    fluid: str = "Air",
) -> dict:
    """
    Engine volumetric and mass flow rate referenced to ambient conditions.

    Inputs:
        displacement_cc : total engine displacement [cc]
        rpm             : engine speed [rev/min]
        VE              : volumetric efficiency [0–1]
        P_amb_kpa       : ambient pressure [kPa]
        T_amb_C         : ambient temperature [°C]
        fluid           : working fluid for CoolProp (default 'Air')

    Returns:
        dict with:
            - Vd_m3        : displacement [m^3]
            - N_int        : intake events per second [-]
            - Vdot_m3_s    : engine volumetric flow [m^3/s]
            - rho_amb      : ambient density [kg/m^3]
            - m_dot_kg_s   : engine mass flow [kg/s]
            - m_dot_lb_hr  : engine mass flow [lb/hr]
    """
    # Displacement in m^3
    Vd_m3 = displacement_cc * 1e-6

    # 4-stroke: one intake event every 2 revs
    N_int = rpm / (2.0 * 60.0)

    # Volumetric flow based on displacement and VE
    Vdot_m3_s = Vd_m3 * VE * N_int

    # Ambient state - can adjust this to P and T conditions of the intake manifold with dyno data, coolprop will calc the density
    P_amb_Pa = kpa_pa(P_amb_kpa)
    T_amb_K = c_k(T_amb_C)

    # Ambient density from CoolProp
    rho_amb = PropsSI("Dmass", "P", P_amb_Pa, "T", T_amb_K, fluid) # change dat shit here

    # Mass flow
    m_dot_kg_s = Vdot_m3_s * rho_amb
    m_dot_lb_hr = kgs_lbhr(m_dot_kg_s)

    return {
        "Vd_m3": Vd_m3,
        "N_int": N_int,
        "Vdot_m3_s": Vdot_m3_s,
        "rho_amb": rho_amb,
        "m_dot_kg_s": m_dot_kg_s,
        "m_dot_lb_hr": m_dot_lb_hr,
    }

def compute_engine_mass_flow_map(
    displacement_cc: float,
    rpm: float,
    VE: float,
    MAP_Pa: float,
    T_amb_C: float,
    fluid: str = "Air",
) -> dict:
    """
    Engine mass flow rate using MAP (manifold absolute pressure) conditions.
    More accurate than ambient-referenced calculation.

    Inputs:
        displacement_cc : total engine displacement [cc]
        rpm             : engine speed [rev/min]
        VE              : volumetric efficiency [0–1]
        MAP_Pa          : manifold absolute pressure [Pa]
        T_amb_C         : ambient temperature [°C] (assuming isothermal)
        fluid           : working fluid for CoolProp (default 'Air')

    Returns:
        dict with:
            - Vd_m3        : displacement [m^3]
            - N_int        : intake events per second [-]
            - Vdot_m3_s    : engine volumetric flow [m^3/s]
            - rho_map      : density at MAP conditions [kg/m^3]
            - m_dot_kg_s   : engine mass flow [kg/s]
            - m_dot_lb_hr  : engine mass flow [lb/hr]
    """
    # Displacement in m^3
    Vd_m3 = displacement_cc * 1e-6

    # 4-stroke: one intake event every 2 revs
    N_int = rpm / (2.0 * 60.0)

    # Volumetric flow based on displacement and VE
    Vdot_m3_s = Vd_m3 * VE * N_int

    # Density at MAP conditions
    T_amb_K = c_k(T_amb_C)
    rho_map = PropsSI("Dmass", "P", MAP_Pa, "T", T_amb_K, fluid)

    # Mass flow
    m_dot_kg_s = Vdot_m3_s * rho_map
    m_dot_lb_hr = kgs_lbhr(m_dot_kg_s)

    return {
        "Vd_m3": Vd_m3,
        "N_int": N_int,
        "Vdot_m3_s": Vdot_m3_s,
        "rho_map": rho_map,
        "m_dot_kg_s": m_dot_kg_s,
        "m_dot_lb_hr": m_dot_lb_hr,
    }

# Pumping power calculation

# ----------------------------------------------------------------------
# Pumping power through intake (using measured MAP)
#
# Assumptions:
# - Upstream ambient pressure is P_t (ambient before restrictor).
# - Downstream pressure for the intake system is MAP (measured).
# - Volumetric flow rate through the intake equals engine volumetric flow
#   referenced to ambient (Vdot_m3_s from compute_engine_mass_flow_ambient).
#
# Equation:
#   Δp      = P_t - P_MAP        [Pa]
#   P_pump  = Δp * Vdot_eng      [W = J/s]
# ----------------------------------------------------------------------
def compute_pumping_power(
    P_t_Pa: float,
    MAP_Pa: float,
    Vdot_m3_s: float,
) -> dict:
    """
    Compute realistic intake pumping power.

    Inputs:
        P_t_Pa     : upstream total / ambient pressure [Pa]
        MAP_Pa     : manifold absolute pressure [Pa]
        Vdot_m3_s  : engine volumetric flow rate [m^3/s]

    Returns:
        dict with:
            - delta_p_Pa : pressure drop across intake [Pa]
            - P_pump_W   : pumping power [W]
            - P_pump_kW  : pumping power [kW]
            - P_pump_hp  : pumping power [hp]
    """
    delta_p_Pa = P_t_Pa - MAP_Pa
    P_pump_W = delta_p_Pa * Vdot_m3_s
    P_pump_kW = P_pump_W / 1000.0
    P_pump_hp = P_pump_W / 745.7  # 1 hp ≈ 745.7 W

    return {
        "delta_p_Pa": delta_p_Pa,
        "P_pump_W": P_pump_W,
        "P_pump_kW": P_pump_kW,
        "P_pump_hp": P_pump_hp,
    }

# Engine total power calculation

# ----------------------------------------------------------------------
# Engine power estimate from air mass flow
#
# Workflow:
# 1) Engine air mass flow (kg/s) - from engine mass flow calculation
# 2) Assume stoichiometric AFR
# 3) Compute fuel mass flow: m_dot_fuel = m_dot_air / AFR.
# 4) Convert to power: P_chem = m_dot_fuel * LHV.
# 5) Apply overall efficiency to get brake power: P_brake = eta * P_chem.
#
# Units:
# - m_dot_air_kg_s in kg/s
# - LHV in J/kg
# - P_chem, P_brake in W; then converted to kW and hp.
# ----------------------------------------------------------------------
def compute_engine_power_from_airflow(
    m_dot_air_kg_s: float,
    AFR: float = 14.7,
    LHV_J_per_kg: float = 43e6,
    eta: float = 0.30,
) -> dict:
    """
    Estimate engine brake power from air mass flow.

    Inputs:
        m_dot_air_kg_s : air mass flow [kg/s]
        AFR            : air-fuel ratio (mass basis) [-]
        LHV_J_per_kg   : fuel lower heating value [J/kg]
        eta            : overall brake efficiency [-]

    Returns:
        dict with:
            - m_dot_fuel_kg_s : fuel mass flow [kg/s]
            - P_chem_W        : chemical power [W]
            - P_chem_kW       : chemical power [kW]
            - P_brake_W       : brake power [W]
            - P_brake_kW      : brake power [kW]
            - P_brake_hp      : brake power [hp]
    """
    # Fuel mass flow at given AFR
    m_dot_fuel_kg_s = m_dot_air_kg_s / AFR

    # Chemical power from fuel
    P_chem_W = m_dot_fuel_kg_s * LHV_J_per_kg
    P_chem_kW = P_chem_W / 1000.0

    # Brake power with efficiency
    P_brake_W = eta * P_chem_W
    P_brake_kW = P_brake_W / 1000.0
    P_brake_hp = P_brake_W / 745.7  # W -> hp

    return {
        "m_dot_fuel_kg_s": m_dot_fuel_kg_s,
        "P_chem_W": P_chem_W,
        "P_chem_kW": P_chem_kW,
        "P_brake_W": P_brake_W,
        "P_brake_kW": P_brake_kW,
        "P_brake_hp": P_brake_hp,
    }


def plot_restrictor_vs_engine_flow(
    d_mm: float,
    Cd: float,
    P_amb_kpa: float,
    T_amb_C: float,
    displacement_cc: float,
    VE: float,
    MAP_Pa: float,
    rpm_min: float = 2000.0,
    rpm_max: float = 14000.0,
    rpm_step: float = 100.0,
    fluid: str = "Air",
    save_plot: bool = False,
    filename: str = "restrictor_vs_engine_flow.png",
) -> dict:
    """
    Plot restrictor choked flow (horizontal line) vs engine mass flow vs RPM.
    Shows where they intersect (where restrictor becomes limiting).
    Uses MAP for engine flow calculation.

    Inputs:
        d_mm           : Restrictor throat diameter [mm]
        Cd             : Discharge coefficient [-]
        P_amb_kpa      : Ambient pressure [kPa]
        T_amb_C        : Ambient temperature [°C]
        displacement_cc: Engine displacement [cc]
        VE             : Volumetric efficiency [0-1]
        MAP_Pa         : Manifold absolute pressure [Pa] (constant for all RPM)
        rpm_min        : Minimum RPM for plot [rev/min]
        rpm_max        : Maximum RPM for plot [rev/min]
        rpm_step       : RPM step size [rev/min]
        fluid          : Working fluid for CoolProp (default 'Air')
        save_plot      : Whether to save the plot to file
        filename       : Filename if saving plot

    Returns:
        dict with:
            - restrictor_mdot_kg_s : restrictor choked flow [kg/s]
            - restrictor_mdot_lbhr : restrictor choked flow [lb/hr]
            - intersection_rpm     : RPM where lines intersect [rev/min] (None if no intersection)
            - fig                  : matplotlib figure object
            - ax                   : matplotlib axes object
    """
    # Calculate restrictor choked flow (constant)
    restrictor_result = compute_restrictor_choked_mdot(
        d_mm=d_mm,
        Cd=Cd,
        P_amb_kpa=P_amb_kpa,
        T_amb_C=T_amb_C,
        fluid=fluid,
    )
    restrictor_mdot_kg_s = restrictor_result["m_dot_kg_s"]
    restrictor_mdot_lbhr = restrictor_result["m_dot_lb_hr"]

    # Generate RPM range
    rpm_array = np.arange(rpm_min, rpm_max + rpm_step, rpm_step)
    engine_mdot_kg_s = np.zeros_like(rpm_array)
    engine_mdot_lbhr = np.zeros_like(rpm_array)

    # Calculate engine mass flow for each RPM using MAP
    for i, rpm in enumerate(rpm_array):
        eng_result = compute_engine_mass_flow_map(
            displacement_cc=displacement_cc,
            rpm=rpm,
            VE=VE,
            MAP_Pa=MAP_Pa,
            T_amb_C=T_amb_C,
            fluid=fluid,
        )
        engine_mdot_kg_s[i] = eng_result["m_dot_kg_s"]
        engine_mdot_lbhr[i] = eng_result["m_dot_lb_hr"]

    # Find intersection point (where engine flow >= restrictor flow)
    intersection_rpm = None
    intersection_idx = None
    for i in range(len(rpm_array)):
        if engine_mdot_kg_s[i] >= restrictor_mdot_kg_s:
            intersection_rpm = rpm_array[i]
            intersection_idx = i
            break

    # Create the plot
    fig, ax = plt.subplots(figsize=(10, 6))

    # Plot engine mass flow vs RPM
    ax.plot(rpm_array, engine_mdot_kg_s, 'b-', linewidth=2, label='Engine Mass Flow')

    # Plot restrictor choked flow (horizontal line)
    ax.axhline(y=restrictor_mdot_kg_s, color='r', linestyle='--', linewidth=2, 
               label=f'Restrictor Choked Flow ({restrictor_mdot_kg_s:.4f} kg/s)')

    # Mark intersection point if it exists
    if intersection_rpm is not None:
        ax.plot(intersection_rpm, restrictor_mdot_kg_s, 'ro', markersize=10, 
                label=f'Intersection at {intersection_rpm:.0f} RPM', zorder=5)
        ax.axvline(x=intersection_rpm, color='gray', linestyle=':', linewidth=1.5, alpha=0.7)
        ax.text(intersection_rpm, restrictor_mdot_kg_s * 1.1, 
                f'{intersection_rpm:.0f} RPM', 
                ha='center', va='bottom', fontsize=10, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.5', facecolor='yellow', alpha=0.7))

    # Formatting
    ax.set_xlabel('Engine RPM [rev/min]', fontsize=12, fontweight='bold')
    ax.set_ylabel('Mass Flow Rate [kg/s]', fontsize=12, fontweight='bold')
    ax.set_title('Restrictor Choked Flow vs Engine Mass Flow', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3, linestyle='--')
    ax.legend(loc='best', fontsize=10)

    # Add secondary y-axis for lb/hr
    ax2 = ax.twinx()
    ax2.set_ylabel('Mass Flow Rate [lb/hr]', fontsize=12, fontweight='bold')
    
    # Set secondary axis limits to match primary (convert kg/s to lb/hr)
    y_min_kg_s = min(np.min(engine_mdot_kg_s), restrictor_mdot_kg_s * 0.9)
    y_max_kg_s = max(np.max(engine_mdot_kg_s), restrictor_mdot_kg_s * 1.1)
    y_min_lbhr = kgs_lbhr(y_min_kg_s)
    y_max_lbhr = kgs_lbhr(y_max_kg_s)
    ax2.set_ylim([y_min_lbhr, y_max_lbhr])
    ax2.tick_params(axis='y', labelsize=10)

    # Add text box with key parameters
    MAP_kpa = MAP_Pa / 1000.0
    textstr = f'Restrictor: d={d_mm:.1f} mm, Cd={Cd:.2f}\n'
    textstr += f'Engine: {displacement_cc:.0f} cc, VE={VE:.2f}\n'
    textstr += f'Ambient: P={P_amb_kpa:.1f} kPa, T={T_amb_C:.1f}°C\n'
    textstr += f'MAP: {MAP_kpa:.2f} kPa (constant)'
    if intersection_rpm is not None:
        textstr += f'\n\nRestrictor limits at: {intersection_rpm:.0f} RPM'
    else:
        textstr += f'\n\nRestrictor not limiting in this RPM range'
    
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    ax.text(0.02, 0.98, textstr, transform=ax.transAxes, fontsize=9,
            verticalalignment='top', bbox=props)

    plt.tight_layout()

    if save_plot:
        plt.savefig(filename, dpi=300, bbox_inches='tight')
        print(f"Plot saved to {filename}")

    return {
        "restrictor_mdot_kg_s": restrictor_mdot_kg_s,
        "restrictor_mdot_lbhr": restrictor_mdot_lbhr,
        "intersection_rpm": intersection_rpm,
        "fig": fig,
        "ax": ax,
    }


if __name__ == "__main__":
    # inputs
    d_mm = 20.0
    Cd = 0.90
    P_amb_kpa = 101.325
    T_amb_C = 20.0

    # engine inputs
    displacement_cc = 599.0
    rpm = 12000.0
    VE = 0.95

    # Manifold absolute pressure (from dyno data @ 12000 RPM)
    MAP_psia = 11.7 # flow bench data from Cobb said we choked when TMAP dropped to 10.9 PSI, 
    # a value never recorded on dyno
    MAP_Pa = psi_pa(MAP_psia)

    # Restrictor choked-flow calculation
    restrictor = compute_restrictor_choked_mdot(
        d_mm=d_mm,
        Cd=Cd,
        P_amb_kpa=P_amb_kpa, # ambient before restrictor
        T_amb_C=T_amb_C, # ambient before restrictor
    )
    p_star = compute_p_star(restrictor["P_t"], restrictor["gamma"])

    # Engine mass flow using MAP
    engine = compute_engine_mass_flow_map(
        displacement_cc=displacement_cc,
        rpm=rpm,
        VE=VE,
        MAP_Pa=MAP_Pa, # manifold absolute pressure as measured on dyno
        T_amb_C=T_amb_C, # assuming ambient, didn't record intake manifold temperature
    )

    # Determine limiting mass flow rate
    mfr_restrictor = restrictor["m_dot_kg_s"]
    mfr_engine = engine["m_dot_kg_s"]
    mfr_limiting = min(mfr_restrictor, mfr_engine)

    # Calculate ratio and choking check
    ratio = mfr_engine / mfr_restrictor
    is_choking = ratio >= 1.0

    # Pumping power calculations
    pump_map = compute_pumping_power(
        P_t_Pa=restrictor["P_t"],
        MAP_Pa=MAP_Pa,
        Vdot_m3_s=engine["Vdot_m3_s"],
    )

    pump_pstar = compute_pumping_power(
        P_t_Pa=restrictor["P_t"],
        MAP_Pa=p_star,
        Vdot_m3_s=engine["Vdot_m3_s"],
    )

    # Power estimate from limiting MFR
    power = compute_engine_power_from_airflow(
        m_dot_air_kg_s=mfr_limiting,
        AFR=14.7,
        LHV_J_per_kg=43e6,
        eta=0.30,
    )

    # Net power calculations (subtract pumping losses)
    P_net_map_hp = power["P_brake_hp"] - pump_map["P_pump_hp"]
    P_net_pstar_hp = power["P_brake_hp"] - pump_pstar["P_pump_hp"]

    # Output
    print("=" * 50)
    print(f"Restrictor MFR (kg/s)     : {mfr_restrictor:.5f}")
    print(f"Engine MFR (kg/s)          : {mfr_engine:.5f}")
    print(f"p* (kPa)                   : {p_star/1000:.2f}")
    print(f"Engine / Restrictor Ratio  : {ratio:.3f}")
    print(f"Restrictor Choking Check   : {'YES' if is_choking else 'NO'}")
    print()
    print(f"Pumping Power (MAP) (hp)  : {pump_map['P_pump_hp']:.3f}")
    print(f"Pumping Power (p*) (hp)    : {pump_pstar['P_pump_hp']:.3f}")
    print()
    print(f"Engine Power (hp)          : {power['P_brake_hp']:.2f}")
    print(f"Net Power (MAP approach)  : {P_net_map_hp:.2f} hp")
    print(f"Net Power (p* approach)    : {P_net_pstar_hp:.2f} hp")
    print("=" * 50)

    # Plot restrictor vs engine flow
    print("\nGenerating plot...")
    plot_result = plot_restrictor_vs_engine_flow(
        d_mm=d_mm,
        Cd=Cd,
        P_amb_kpa=P_amb_kpa,
        T_amb_C=T_amb_C,
        displacement_cc=displacement_cc,
        VE=VE,
        MAP_Pa=MAP_Pa,
        rpm_min=2000.0,
        rpm_max=14000.0,
        rpm_step=100.0,
        save_plot=True,
    )
    
    if plot_result["intersection_rpm"] is not None:
        print(f"Restrictor becomes limiting at: {plot_result['intersection_rpm']:.0f} RPM")
    
    plt.show()






