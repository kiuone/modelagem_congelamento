# ================================================FILE: physics.py
from numba import njit

# Constantes Termodinâmicas
L_o = 333600.0  # Calor latente de fusão da água pura (J/kg)

@njit(cache=True)
def calc_properties(T_celsius: float, xw: float, xp: float, xf: float, xc: float, xa: float, T_f: float) -> tuple[float, float, float]:
    """
    Calcula as propriedades termofísicas efetivas (condutividade, densidade e calor específico)
    da matriz alimentar em uma temperatura específica.
    
    A função particiona a água entre ligada e livre. Na região de congelamento (T < T_f),
    aplica o modelo empírico de Chen para determinar a fração de gelo e o calor específico aparente.
    Na região descongelada, aplica os polinômios de Choi & Okos (1986) utilizando o Método de Horner
    para otimização de processamento.
    
    Args:
        T_celsius: Temperatura atual do nó em °C.
        xw, xp, xf, xc, xa: Frações mássicas de água, proteína, gordura, carboidrato e cinzas.
        T_f: Temperatura inicial de congelamento do produto em °C.
        
    Returns:
        tuple: (k_eff, rho_eff, c_app) correspondendo a Condutividade (W/m.K), Densidade (kg/m3) 
               e Calor Específico Aparente (J/kg.K).
    """
    
    # Fração de sólidos totais
    x_s = xp + xf + xc + xa
    
    # Fração de água ligada à proteína (estimativa de Schwartzberg/Chen para carnes)
    x_b = 0.4 * xp 
    
    if T_celsius < T_f:
        # Prevenção de divisão por zero caso a temperatura se aproxime assintoticamente de 0
        T_safe = min(T_celsius, -0.1)
        
        # Modelo de Chen: Calor específico aparente (incorporando o calor latente de fusão)
        c_app = (1.55 + 1.26 * x_s - ((xw - x_b) * (L_o / 1000.0) * T_f) / (T_safe**2)) * 1000.0
        
        # Modelo de Chen: Fração de água convertida em gelo
        x_ice = (xw - x_b) * (1.0 - (T_f / T_safe))
    else:
        # Cálculo de calor sensível por Choi & Okos otimizado via Método de Horner
        cp_w = (4.1762 + T_celsius * (-9.0864e-5 + T_celsius * 5.4731e-6)) * 1000.0
        cp_p = (2.0082 + T_celsius * (1.2089e-3  - T_celsius * 1.3129e-6)) * 1000.0
        cp_f = (1.9842 + T_celsius * (1.4733e-3  - T_celsius * 4.8008e-6)) * 1000.0
        cp_c = (1.5488 + T_celsius * (1.9625e-3  - T_celsius * 5.9399e-6)) * 1000.0
        cp_a = (1.0926 + T_celsius * (1.8896e-3  - T_celsius * 3.6817e-6)) * 1000.0
        
        # Calor específico da mistura descongelada
        c_app = cp_w*xw + cp_p*xp + cp_f*xf + cp_c*xc + cp_a*xa
        x_ice = 0.0

    # Fração de água remanescente em estado líquido
    x_w_unfrozen = max(0.0, xw - x_ice)
    
    # Condutividade térmica dos componentes individuais (Choi & Okos via Horner)
    k_w = 0.57109 + T_celsius * (1.7625e-3 - T_celsius * 6.7036e-6)
    k_p = 0.17881 + T_celsius * (1.1958e-3 - T_celsius * 2.7178e-6)
    k_f = 0.18071 - T_celsius * (2.7604e-3 + T_celsius * 1.7749e-7)
    k_c = 0.20141 + T_celsius * (1.3874e-3 - T_celsius * 4.3312e-6)
    k_a = 0.32962 + T_celsius * (1.4011e-3 - T_celsius * 2.9069e-6)
    k_ice = 2.2196 - T_celsius * (6.2489e-3 - T_celsius * 1.0154e-4)
    
    # Modelo estrutural paralelo para condutividade térmica efetiva
    k_eff = k_w*x_w_unfrozen + k_ice*x_ice + k_p*xp + k_f*xf + k_c*xc + k_a*xa

    # Densidade dos componentes individuais (Choi & Okos via Horner)
    rho_w = 997.18 + T_celsius * (3.1439e-3 - T_celsius * 3.7574e-3)
    rho_p = 1329.9 - 0.5184 * T_celsius
    rho_f = 925.59 - 0.41757 * T_celsius
    rho_c = 1599.1 - 0.31046 * T_celsius
    rho_a = 2423.8 - 0.28063 * T_celsius
    rho_ice = 916.89 - 0.13071 * T_celsius

    # Densidade efetiva calculada por conservação de volume
    den_sum = (x_w_unfrozen/rho_w + x_ice/rho_ice + xp/rho_p + xf/rho_f + xc/rho_c + xa/rho_a)
    rho_eff = 1.0 / den_sum if den_sum > 0.0 else 1000.0

    return k_eff, rho_eff, c_app
