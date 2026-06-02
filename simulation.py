import numpy as np
import polars as pl
from numba import njit
from physics import calc_properties

# Este é o motor isolado que vai rodar em velocidade C++ (cache=True)
@njit(cache=True)
def _core_sim_loop(steps, N, dt, dx, T, T_inf, U_top, U_bottom, record_interval, xw, xp, xf, xc, xa, T_f):
    # Pré-alocando as matrizes (Numba adora matrizes de tamanho fixo)
    num_records = (steps // record_interval) + 1
    time_records = np.zeros(num_records)
    top_records = np.zeros(num_records)
    center_records = np.zeros(num_records)
    bottom_records = np.zeros(num_records)

    T_new = np.copy(T)
    k_arr = np.zeros(N)
    rho_arr = np.zeros(N)
    cp_arr = np.zeros(N)

    record_idx = 0

    # Este é o gargalo que agora rodará em milissegundos
    for step in range(steps):
        # Atualiza as propriedades térmicas para cada nó
        for i in range(N):
            k_arr[i], rho_arr[i], cp_arr[i] = calc_properties(T[i], xw, xp, xf, xc, xa, T_f)
            
        # Calcula a condução térmica interna
        for i in range(1, N-1):
            k_right = (k_arr[i] + k_arr[i+1]) / 2.0
            k_left = (k_arr[i] + k_arr[i-1]) / 2.0
            term_right = k_right * (T[i+1] - T[i])
            term_left = k_left * (T[i] - T[i-1])
            T_new[i] = T[i] + (dt / (rho_arr[i] * cp_arr[i] * dx**2)) * (term_right - term_left)

        # Contorno Superior (Topo)
        q_top = U_top * (T_inf - T[0])
        T_new[0] = T[0] + (dt / (rho_arr[0] * cp_arr[0] * dx)) * (q_top + k_arr[0] * (T[1] - T[0]) / dx)

        # Contorno Inferior (Fundo)
        q_bottom = U_bottom * (T_inf - T[N-1])
        T_new[N-1] = T[N-1] + (dt / (rho_arr[N-1] * cp_arr[N-1] * dx)) * (q_bottom + k_arr[N-1] * (T[N-2] - T[N-1]) / dx)
        
        # Atualiza os nós para o próximo ciclo
        for i in range(N):
            T[i] = T_new[i]
        
        # Grava os dados apenas no intervalo configurado para economizar memória
        if step % record_interval == 0:
            time_records[record_idx] = (step * dt) / 60.0
            top_records[record_idx] = T[0]
            center_records[record_idx] = T[N // 2]
            bottom_records[record_idx] = T[N-1]
            record_idx += 1
            
    # O Numba corta os arrays para o tamanho exato preenchido e devolve para o Python
    return time_records[:record_idx], top_records[:record_idx], center_records[:record_idx], bottom_records[:record_idx]

def run_simulation(tunnel_params, box_params, comp_params, sim_params):
    # Setup dos parâmetros da malha geométrica
    H = box_params['height'] / 1000.0
    N = int(sim_params['nodes'])
    dx = H / (N - 1)
    dt = float(sim_params['dt'])
    time_total = sim_params['total_time_min'] * 60
    
    T_inf = float(tunnel_params['temp'])
    h_conv = float(tunnel_params['h'])
    
    # Cálculo das resistências da embalagem
    R_top = box_params['plastic_thickness'] / 1000.0 / box_params['plastic_k']
    U_top = 1.0 / ((1.0 / h_conv) + R_top)
    
    R_bottom = (box_params['carton_thickness'] / 1000.0 / box_params['carton_k'] + 
                box_params['plastic_thickness'] / 1000.0 / box_params['plastic_k'])
    U_bottom = 1.0 / ((1.0 / h_conv) + R_bottom)

    T = np.full(N, comp_params['t_initial'], dtype=np.float64)
    record_interval = max(1, int(60 / dt))
    steps = int(time_total / dt)

    # Extraindo as composições para enviar limpas para o Numba
    xw = float(comp_params['water'])
    xp = float(comp_params['protein'])
    xf = float(comp_params['fat'])
    xc = float(comp_params['carbohydrate'])
    xa = float(comp_params['ash'])
    T_f = float(comp_params['t_f'])

    # --- CHAMA O MOTOR NUMBA ---
    # Aqui a mágica de velocidade acontece
    time_res, top_res, center_res, bottom_res = _core_sim_loop(
        steps, N, dt, dx, T, T_inf, U_top, U_bottom, record_interval,
        xw, xp, xf, xc, xa, T_f
    )

    # O Polars assume o controle e organiza a tabela
    df_results = pl.DataFrame({
        "Tempo (min)": time_res,
        "Superfície Topo (°C)": top_res,
        "Centro Térmico (°C)": center_res,
        "Superfície Fundo (°C)": bottom_res
    })
    
    return df_results