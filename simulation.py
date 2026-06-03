# ================================================FILE: simulation.py
import numpy as np
import polars as pl
from numba import njit
from physics import calc_properties

@njit(cache=True)
def _core_sim_loop(steps, N, dt, dx, T, T_inf, U_top, U_bottom, record_interval, 
                   xw, xp, xf, xc, xa, T_f):
    """
    Motor numérico explícito em diferenças finitas (1D) compilado via Numba.
    Resolve a equação de condução de calor transiente com propriedades variáveis.
    """
    
    # Pré-alocação de vetores de armazenamento para evitar overhead de memória
    num_records = (steps // record_interval) + 1
    time_records = np.zeros(num_records)
    top_records = np.zeros(num_records)
    center_records = np.zeros(num_records)
    bottom_records = np.zeros(num_records)
    hotspot_records = np.zeros(num_records)
    
    T_new = np.copy(T)
    record_idx = 0
    
    for step in range(steps):
        # 1. Avaliação dos nós internos (Condução pura)
        # O laço não afeta as fronteiras (índices 0 e N-1)
        for i in range(1, N - 1):
            k, rho, cp = calc_properties(T[i], xw, xp, xf, xc, xa, T_f)
            alpha = k / (rho * cp)
            T_new[i] = T[i] + alpha * dt * (T[i+1] - 2*T[i] + T[i-1]) / (dx**2)
            
        # 2. Avaliação da fronteira superior (Convecção + Filme Plástico)
        # Utiliza-se o fator 2 para compensar o meio-volume de controle da malha de contorno
        k_top, rho_top, cp_top = calc_properties(T[0], xw, xp, xf, xc, xa, T_f)
        T_new[0] = T[0] + (2.0 * dt / (rho_top * cp_top * dx)) * (
            U_top * (T_inf - T[0]) - k_top * (T[0] - T[1]) / dx
        )
        
        # 3. Avaliação da fronteira inferior (Convecção + Filme + Papelão + Esteira)
        k_bot, rho_bot, cp_bot = calc_properties(T[N-1], xw, xp, xf, xc, xa, T_f)
        T_new[N-1] = T[N-1] + (2.0 * dt / (rho_bot * cp_bot * dx)) * (
            U_bottom * (T_inf - T[N-1]) - k_bot * (T[N-1] - T[N-2]) / dx
        )
        
        # Atualização simultânea da malha para o próximo passo de tempo
        for i in range(N):
            T[i] = T_new[i]
        
        # Extração de dados apenas nos intervalos configurados
        if step % record_interval == 0:
            time_records[record_idx] = (step * dt) / 3600.0 # Conversão para Horas
            top_records[record_idx] = T[0]
            center_records[record_idx] = T[N // 2]
            bottom_records[record_idx] = T[N - 1]
            # O ponto mais quente acompanha o deslocamento do centro térmico devido à assimetria
            hotspot_records[record_idx] = np.max(T) 
            record_idx += 1
            
    return time_records[:record_idx], top_records[:record_idx], center_records[:record_idx], bottom_records[:record_idx], hotspot_records[:record_idx]

def run_simulation(tunnel_params, box_params, comp_params, sim_params):
    """
    Prepara as condições de contorno, resistências térmicas e inicializa o laço compilado.
    Retorna um DataFrame estruturado no Polars.
    """
    H = box_params['altura'] / 1000.0
    N = int(sim_params['nodes'])
    dx = H / (N - 1)
    dt = float(sim_params['dt'])
    
    time_total = sim_params['total_time_min'] * 60.0
    steps = int(time_total / dt)
    record_interval = max(1, int(60 / dt)) # Amostragem a cada 1 minuto de simulação
    
    h_tunel = float(tunnel_params['h'])
    
    # Cálculo das Resistências Térmicas (R = espessura / condutividade)
    R_plastico = (box_params['f_espessura'] / 1000.0) / box_params['f_k']
    R_papelao = (box_params['p_espessura'] / 1000.0) / box_params['p_k']
    
    # Coeficiente Global de Transferência de Calor (U) - Face Superior (Destampada)
    U_top = 1.0 / ((1.0 / h_tunel) + R_plastico)
    
    # Coeficiente Global de Transferência de Calor (U) - Face Inferior
    # Fator de 0.6 estima perda de área útil de convecção devido ao contato físico com a esteira
    U_bottom = (1.0 / ((1.0 / h_tunel) + R_plastico + R_papelao)) * 0.6
    
    T_inicial = float(comp_params['t_initial'])
    T = np.full(N, T_inicial, dtype=np.float64)
    T_inf = float(tunnel_params['temp'])
    
    # Chamada da função compilada em C++
    t_rec, top_rec, center_rec, bot_rec, hot_rec = _core_sim_loop(
        steps, N, dt, dx, T, T_inf, U_top, U_bottom, record_interval,
        float(comp_params['xw']), float(comp_params['xp']), float(comp_params['xf']),
        float(comp_params['xc']), float(comp_params['xa']), float(comp_params['T_f'])
    )
    
    # Estruturação colunar de alta performance
    return pl.DataFrame({
        "Tempo (h)": t_rec,
        "Superfície Topo": top_rec,
        "Centro Geométrico": center_rec,
        "Superfície Fundo": bot_rec,
        "Ponto Mais Quente": hot_rec
    })
