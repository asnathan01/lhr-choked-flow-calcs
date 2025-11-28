import pint

from pint import UnitRegistry
ureg = UnitRegistry()

# --------------- MASS CONVERSIONS ---------------
def lb_kg(value):       return (value * ureg.pound).to(ureg.kilogram).magnitude
def kg_lb(value):       return (value * ureg.kilogram).to(ureg.pound).magnitude
def oz_g(value):        return (value * ureg.ounce).to(ureg.gram).magnitude
def g_oz(value):        return (value * ureg.gram).to(ureg.ounce).magnitude
def tonUS_kg(value):    return (value * ureg.short_ton).to(ureg.kilogram).magnitude
def kg_tonUS(value):    return (value * ureg.kilogram).to(ureg.short_ton).magnitude
def tonmetric_kg(value):return (value * ureg.metric_ton).to(ureg.kilogram).magnitude
def kg_tonmetric(value):return (value * ureg.kilogram).to(ureg.metric_ton).magnitude
def mg_g(value):        return (value * ureg.milligram).to(ureg.gram).magnitude
def g_mg(value):        return (value * ureg.gram).to(ureg.milligram).magnitude
def lb_oz(value):       return (value * ureg.pound).to(ureg.ounce).magnitude
def oz_lb(value):       return (value * ureg.ounce).to(ureg.pound).magnitude
def slug_kg(value):     return (value * ureg.slug).to(ureg.kilogram).magnitude
def kg_slug(value):     return (value * ureg.kilogram).to(ureg.slug).magnitude


# --------------- TEMPERATURE CONVERSIONS (VERSION-SAFE) ---------------
def c_f(value):         return ureg.Quantity(value, ureg.degC).to(ureg.degF).magnitude
def f_c(value):         return ureg.Quantity(value, ureg.degF).to(ureg.degC).magnitude
def c_k(value):         return ureg.Quantity(value, ureg.degC).to(ureg.kelvin).magnitude
def k_c(value):         return ureg.Quantity(value, ureg.kelvin).to(ureg.degC).magnitude
def f_k(value):         return ureg.Quantity(value, ureg.degF).to(ureg.kelvin).magnitude
def k_f(value):         return ureg.Quantity(value, ureg.kelvin).to(ureg.degF).magnitude
def r_f(value):         return ureg.Quantity(value, ureg.degR).to(ureg.degF).magnitude
def f_r(value):         return ureg.Quantity(value, ureg.degF).to(ureg.degR).magnitude
def r_k(value):         return ureg.Quantity(value, ureg.degR).to(ureg.kelvin).magnitude
def k_r(value):         return ureg.Quantity(value, ureg.kelvin).to(ureg.degR).magnitude


# --------------- LENGTH CONVERSIONS ---------------
def in_mm(value):       return (value * ureg.inch).to(ureg.millimeter).magnitude
def mm_in(value):       return (value * ureg.millimeter).to(ureg.inch).magnitude
def in_cm(value):       return (value * ureg.inch).to(ureg.centimeter).magnitude
def cm_in(value):       return (value * ureg.centimeter).to(ureg.inch).magnitude
def ft_m(value):        return (value * ureg.foot).to(ureg.meter).magnitude
def m_ft(value):        return (value * ureg.meter).to(ureg.foot).magnitude
def yd_m(value):        return (value * ureg.yard).to(ureg.meter).magnitude
def m_yd(value):        return (value * ureg.meter).to(ureg.yard).magnitude
def mi_km(value):       return (value * ureg.mile).to(ureg.kilometer).magnitude
def km_mi(value):       return (value * ureg.kilometer).to(ureg.mile).magnitude
def mm_m(value):        return (value * ureg.millimeter).to(ureg.meter).magnitude
def m_mm(value):        return (value * ureg.meter).to(ureg.millimeter).magnitude

# --------------- MASS FLOW RATE CONVERSIONS ---------------
def lbhr_kgs(value):    return (value * ureg.pound / ureg.hour).to(ureg.kilogram / ureg.second).magnitude
def kgs_lbhr(value):    return (value * ureg.kilogram / ureg.second).to(ureg.pound / ureg.hour).magnitude
def gs_kgs(value):      return (value * ureg.gram / ureg.second).to(ureg.kilogram / ureg.second).magnitude
def kgs_gs(value):      return (value * ureg.kilogram / ureg.second).to(ureg.gram / ureg.second).magnitude
def kghr_kgs(value):    return (value * ureg.kilogram / ureg.hour).to(ureg.kilogram / ureg.second).magnitude
def kgs_kghr(value):    return (value * ureg.kilogram / ureg.second).to(ureg.kilogram / ureg.hour).magnitude

# --------------- VOLUMETRIC FLOW RATE CONVERSIONS ---------------
def cms_ls(value):      return (value * ureg.cubic_meter / ureg.second).to(ureg.liter / ureg.second).magnitude
def ls_cms(value):      return (value * ureg.liter / ureg.second).to(ureg.cubic_meter / ureg.second).magnitude
def lpm_gpm(value):     return (value * ureg.liter / ureg.minute).to(ureg.gallon / ureg.minute).magnitude
def gpm_lpm(value):     return (value * ureg.gallon / ureg.minute).to(ureg.liter / ureg.minute).magnitude
def cfm_cmh(value):     return (value * ureg.cubic_foot / ureg.minute).to(ureg.cubic_meter / ureg.hour).magnitude
def cmh_cfm(value):     return (value * ureg.cubic_meter / ureg.hour).to(ureg.cubic_foot / ureg.minute).magnitude
def cfs_cms(value):     return (value * ureg.cubic_foot / ureg.second).to(ureg.cubic_meter / ureg.second).magnitude
def cms_cfs(value):     return (value * ureg.cubic_meter / ureg.second).to(ureg.cubic_foot / ureg.second).magnitude

# --------------- DENSITY CONVERSIONS ---------------
def kg_m3_slug_ft3(value): return (value * ureg.kilogram / ureg.meter**3).to(ureg.slug / ureg.foot**3).magnitude
def slug_ft3_kg_m3(value): return (value * ureg.slug / ureg.foot**3).to(ureg.kilogram / ureg.meter**3).magnitude
def kg_m3_lb_ft3(value):   return (value * ureg.kilogram / ureg.meter**3).to(ureg.pound / ureg.foot**3).magnitude
def lb_ft3_kg_m3(value):   return (value * ureg.pound / ureg.foot**3).to(ureg.kilogram / ureg.meter**3).magnitude

# --------------- PRESSURE CONVERSIONS ---------------
def pa_bar(value):      return (value * ureg.pascal).to(ureg.bar).magnitude
def bar_pa(value):      return (value * ureg.bar).to(ureg.pascal).magnitude
def pa_kpa(value):      return (value * ureg.pascal).to(ureg.kilopascal).magnitude
def kpa_pa(value):      return (value * ureg.kilopascal).to(ureg.pascal).magnitude
def pa_psi(value):      return (value * ureg.pascal).to(ureg.psi).magnitude
def psi_pa(value):      return (value * ureg.psi).to(ureg.pascal).magnitude
def atm_pa(value):      return (value * ureg.atm).to(ureg.pascal).magnitude
def pa_atm(value):      return (value * ureg.pascal).to(ureg.atm).magnitude
def atm_psi(value):     return (value * ureg.atm).to(ureg.psi).magnitude
def psi_atm(value):     return (value * ureg.psi).to(ureg.atm).magnitude
def bar_psi(value):     return (value * ureg.bar).to(ureg.psi).magnitude
def psi_bar(value):     return (value * ureg.psi).to(ureg.bar).magnitude
def torr_pa(value):     return (value * ureg.torr).to(ureg.pascal).magnitude
def pa_torr(value):     return (value * ureg.pascal).to(ureg.torr).magnitude
def inhg_pa(value):     return (value * ureg.inch_Hg).to(ureg.pascal).magnitude
def pa_inhg(value):     return (value * ureg.pascal).to(ureg.inch_Hg).magnitude
def mH2O_kpa(value):    return (value * ureg.meter * ureg.water).to(ureg.kilopascal).magnitude
def kpa_mH2O(value):    return (value * ureg.kilopascal).to(ureg.meter * ureg.water).magnitude

# --------------- POWER CONVERSIONS ---------------
def hp_kw(value):       return (value * ureg.horsepower).to(ureg.kilowatt).magnitude
def kw_hp(value):       return (value * ureg.kilowatt).to(ureg.horsepower).magnitude
def hp_w(value):        return (value * ureg.horsepower).to(ureg.watt).magnitude
def w_hp(value):        return (value * ureg.watt).to(ureg.horsepower).magnitude
def kw_w(value):        return (value * ureg.kilowatt).to(ureg.watt).magnitude
def w_kw(value):        return (value * ureg.watt).to(ureg.kilowatt).magnitude
def mw_kw(value):       return (value * ureg.megawatt).to(ureg.kilowatt).magnitude
def kw_mw(value):       return (value * ureg.kilowatt).to(ureg.megawatt).magnitude
def mw_w(value):        return (value * ureg.megawatt).to(ureg.watt).magnitude
def w_mw(value):        return (value * ureg.watt).to(ureg.megawatt).magnitude
def btuhr_w(value):     return (value * ureg.BTU / ureg.hour).to(ureg.watt).magnitude
def w_btuhr(value):     return (value * ureg.watt).to(ureg.BTU / ureg.hour).magnitude
def btuhr_kw(value):    return (value * ureg.BTU / ureg.hour).to(ureg.kilowatt).magnitude
def kw_btuhr(value):    return (value * ureg.kilowatt).to(ureg.BTU / ureg.hour).magnitude
def btuhr_hp(value):    return (value * ureg.BTU / ureg.hour).to(ureg.horsepower).magnitude
def hp_btuhr(value):    return (value * ureg.horsepower).to(ureg.BTU / ureg.hour).magnitude