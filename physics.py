import numpy as np
from numba import njit

# Constantes termodinâmicas
L_o = 333.6 * 1000  # Calor latente de fusão da água (J/kg) [cite: 304]
T_o = 273.15        # Zero Celsius em Kelvin [cite: 334]

@njit(cache=True)
def calc_properties(T_celsius, xw, xp, xf, xc, xa, T_f):
    """
    Calcula as propriedades termofísicas baseadas em Choi & Okos (1986) [cite: 271]
    e Chen (1985) para a mudança de fase[cite: 337, 338].
    """
    # 1. Calor Específico (C) acima do congelamento (J/kg.K) [cite: 297, 329]
    cp_w = (4.1762 - 9.0864e-5*T_celsius + 5.4731e-6*T_celsius**2) * 1000
    cp_p = (2.0082 + 1.2089e-3*T_celsius - 1.3129e-6*T_celsius**2) * 1000
    cp_f = (1.9842 + 1.4733e-3*T_celsius - 4.8008e-6*T_celsius**2) * 1000
    cp_c = (1.5488 + 1.9625e-3*T_celsius - 5.9399e-6*T_celsius**2) * 1000
    cp_a = (1.0926 + 1.8896e-3*T_celsius - 3.6817e-6*T_celsius**2) * 1000
    
    c_u = cp_w*xw + cp_p*xp + cp_f*xf + cp_c*xc + cp_a*xa

    # 2. Fração de sólidos e água ligada [cite: 306, 331]
    x_s = xp + xf + xc + xa
    x_b = 0.4 * xp 
    
    # 3. Calor específico aparente (C_app) para região congelada (Modelo de Chen) 
    # Usado para lidar com o calor latente sendo liberado gradativamente
    if T_celsius < T_f:
        # Fator de conversão e proteção contra divisão por zero
        T_safe = min(T_celsius, -0.1)
        c_app = (1.55 + 1.26*x_s - ((xw - x_b) * (L_o/1000) * T_f) / (T_safe**2)) * 1000
    else:
        c_app = c_u

    # 4. Condutividade Térmica Média (k) (W/m.K) [cite: 297]
    k_w = 0.57109 + 1.7625e-3*T_celsius - 6.7036e-6*T_celsius**2
    k_p = 0.17881 + 1.1958e-3*T_celsius - 2.7178e-6*T_celsius**2
    k_f = 0.18071 - 2.7604e-3*T_celsius - 1.7749e-7*T_celsius**2
    k_c = 0.20141 + 1.3874e-3*T_celsius - 4.3312e-6*T_celsius**2
    k_a = 0.32962 + 1.4011e-3*T_celsius - 2.9069e-6*T_celsius**2
    k_ice = 2.2196 - 6.2489e-3*T_celsius + 1.0154e-4*T_celsius**2

    # Se estiver congelando, parte da água vira gelo [cite: 319]
    x_ice = 0.0
    if T_celsius < T_f:
        x_ice = (xw - x_b) * (1 - (T_f / min(T_celsius, -0.1)))
    
    x_w_unfrozen = max(0.0, xw - x_ice)
    
    # Modelo paralelo simplificado de condutividade térmica [cite: 376]
    k_eff = k_w*x_w_unfrozen + k_ice*x_ice + k_p*xp + k_f*xf + k_c*xc + k_a*xa

    # 5. Densidade Efetiva (rho) (kg/m3) [cite: 297]
    rho_w = 997.18 + 3.1439e-3*T_celsius - 3.7574e-3*T_celsius**2
    rho_p = 1329.9 - 0.5184*T_celsius
    rho_f = 925.59 - 0.41757*T_celsius
    rho_c = 1599.1 - 0.31046*T_celsius
    rho_a = 2423.8 - 0.28063*T_celsius
    rho_ice = 916.89 - 0.13071*T_celsius

    # Equação da densidade por conservação de volume [cite: 322]
    # Prevenindo divisão por zero
    den_sum = (x_w_unfrozen/rho_w + x_ice/rho_ice + xp/rho_p + 
               xf/rho_f + xc/rho_c + xa/rho_a)
    
    if den_sum > 0:
        rho_eff = 1.0 / den_sum
    else:
        rho_eff = 1000.0

    return k_eff, rho_eff, c_app