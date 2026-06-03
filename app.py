# ================================================FILE: app.py
import streamlit as st
import plotly.express as px
from simulation import run_simulation

# Configuração global da página
st.set_page_config(page_title="Simulador de Congelamento CVale", layout="wide")
st.title("🧊 Simulador de Congelamento Industrial")
st.markdown("Modelo 1D por Diferenças Finitas - Otimizado para Caixas Destampadas")

# ==============================================================================
# BARRA LATERAL: ENTRADA DE PARÂMETROS
# ==============================================================================
st.sidebar.header("Parâmetros de Controle")

# --- 1. TÚNEL DE CONGELAMENTO ---
with st.sidebar.expander("🌬️ Condições do Túnel", expanded=True):
    tunel_opcoes = {
        "Blast Freezer (Alta Vel.)": {"temp": -33.0, "h": 120.0},
        "Túnel Convencional": {"temp": -28.0, "h": 50.0},
        "Câmara Estática (S/ Vento)": {"temp": -22.0, "h": 15.0},
        "Customizado": {"temp": -30.0, "h": 60.0}
    }
    tunel_selecionado = st.selectbox("Perfil Operacional", list(tunel_opcoes.keys()))
    
    # Variáveis de estado do túnel
    t_temp = st.number_input("Temperatura do Ar (°C)", value=tunel_opcoes[tunel_selecionado]["temp"], step=1.0)
    t_h = st.number_input("Coef. Convectivo 'h' (W/m².K)", value=tunel_opcoes[tunel_selecionado]["h"], step=5.0, min_value=5.0, max_value=300.0)
    st.caption("*Nota: Valores de 'h' abaixo de 20 representam ar quase parado. Valores acima de 60 representam ventilação forçada extrema.*")

# --- 2. GEOMETRIA E EMBALAGEM ---
with st.sidebar.expander("📦 Propriedades da Caixa", expanded=False):
    # Dimensões globais para avaliação de dimensionalidade (Fator E) e malha 1D
    c_altura = st.number_input("Altura do Bloco de Carne (mm)", value=118.8, step=1.0)
    c_largura = st.number_input("Largura da Caixa (mm)", value=385.6, step=1.0)
    c_comprimento = st.number_input("Comprimento da Caixa (mm)", value=575.6, step=1.0)
    
    st.markdown("---")
    st.markdown("**Isolamento Inferior e Lateral (Papelão)**")
    p_espessura = st.number_input("Espessura Papelão (mm)", value=3.8, step=0.1)
    p_k = st.number_input("Condutividade Papelão (W/m.K)", value=0.064, step=0.01)
    
    st.markdown("**Isolamento de Contato (Filme Plástico)**")
    f_espessura = st.number_input("Espessura Plástico (mm)", value=0.1, step=0.05)
    f_k = st.number_input("Condutividade Plástico (W/m.K)", value=0.33, step=0.05)

# --- 3. MATRIZ BIOLÓGICA ---
with st.sidebar.expander("🍗 Composição do Produto", expanded=False):
    st.markdown("Frações Mássicas (A soma ideal é 1.0)")
    comp_w = st.slider("Água", 0.0, 1.0, 0.748)
    comp_p = st.slider("Proteína", 0.0, 1.0, 0.224)
    comp_f = st.slider("Gordura", 0.0, 1.0, 0.011)
    comp_c = st.slider("Carboidrato", 0.0, 1.0, 0.000)
    comp_a = st.slider("Cinzas", 0.0, 1.0, 0.017)
    
    soma_comp = comp_w + comp_p + comp_f + comp_c + comp_a
    
    st.markdown("---")
    t_initial = st.number_input("Temperatura de Entrada (°C)", value=5.0, step=1.0)
    T_f = st.number_input("Ponto de Congelamento Inicial (°C)", value=-2.0, step=0.1)

# --- 4. PARÂMETROS NUMÉRICOS ---
with st.sidebar.expander("⚙️ Motor Numérico (Diferenças Finitas)", expanded=False):
    sim_time = st.number_input("Tempo de Simulação (horas)", value=14.0, step=1.0)
    sim_nodes = st.number_input("Número de Nós (Fatias na Altura)", min_value=20, value=150, step=10)
    sim_dt = st.number_input("Passo de Tempo 'dt' (segundos)", value=0.01, step=0.01, format="%.3f")


# ==============================================================================
# PAINEL DE DIAGNÓSTICO E INDICADORES (PRÉ-SIMULAÇÃO)
# ==============================================================================
st.markdown("### 📊 Diagnóstico Termodinâmico e Numérico")

# 1. Validador de Matriz Biológica
if abs(soma_comp - 1.0) > 0.01:
    st.warning(f"⚠️ **Atenção:** A soma da composição centesimal é {soma_comp:.3f}. Para precisão termodinâmica máxima, ajuste os sliders para que a soma seja exata a 1.0.")

col1, col2, col3, col4 = st.columns(4)

# 2. Critério de Estabilidade Numérica (Número de Fourier)
dx_m = (c_altura / 1000.0) / (sim_nodes - 1)
alpha_max = 1.4e-6 # Difusividade aproximada da carne congelada
Fo = (alpha_max * sim_dt) / (dx_m ** 2)

with col1:
    st.metric("Número de Fourier (Fo)", f"{Fo:.4f}")
    if Fo <= 0.5:
        st.success("Malha Estável")
    else:
        st.error("Malha Instável (Reduza o dt!)")

# 3. Fator de Dimensionalidade (E) - Efeito das Bordas
E_factor = 1.0 + (c_altura / c_largura)**2 + (c_altura / c_comprimento)**2
with col2:
    st.metric("Fator de Forma (E)", f"{E_factor:.3f}")
    st.caption(f"Aceleração 3D estimada em +{(E_factor - 1) * 100:.1f}%")

# 4. Números de Biot (Top e Bottom)
# Espelha o cálculo de coeficientes globais (U) presentes no motor de simulação
R_plastico = (f_espessura / 1000.0) / f_k
R_papelao = (p_espessura / 1000.0) / p_k

U_top_UI = 1.0 / ((1.0 / t_h) + R_plastico)
U_bottom_UI = (1.0 / ((1.0 / t_h) + R_plastico + R_papelao)) * 0.6

k_frango_frozen = 1.35 # W/m.K médio pós-transição de fase
L_c = (c_altura / 1000.0) / 2.0 # Comprimento característico (meia espessura)

Bi_top = (U_top_UI * L_c) / k_frango_frozen
Bi_bottom = (U_bottom_UI * L_c) / k_frango_frozen

with col3:
    st.metric("Biot (Topo/Plástico)", f"{Bi_top:.2f}")
    if Bi_top > 10: st.caption("🔴 Domínio: Carne")
    elif Bi_top < 0.1: st.caption("🟡 Domínio: Túnel")
    else: st.caption("🟢 Regime Misto")

with col4:
    st.metric("Biot (Fundo/Papelão)", f"{Bi_bottom:.2f}")
    if Bi_bottom > 10: st.caption("🔴 Domínio: Carne")
    elif Bi_bottom < 0.1: st.caption("🟡 Domínio: Túnel")
    else: st.caption("🟢 Regime Misto")

st.markdown("---")

# ==============================================================================
# EXECUÇÃO DO MOTOR NUMÉRICO
# ==============================================================================
if st.button("▶️ Iniciar Integração Numérica", type="primary", use_container_width=True):
    
    # Trava de Segurança Matemática
    if Fo > 0.5:
        st.error("🚨 Simulação Abortada. O Critério de Fourier foi violado. Reduza o 'dt' ou diminua a quantidade de Nós.")
    else:
        # Empacotamento de Dicionários (Garante alinhamento exato com simulation.py)
        tunnel_params = {'temp': t_temp, 'h': t_h}
        box_params = {
            'altura': c_altura, 
            'p_espessura': p_espessura, 'p_k': p_k,
            'f_espessura': f_espessura, 'f_k': f_k
        }
        comp_params = {
            'xw': comp_w, 'xp': comp_p, 'xf': comp_f, 
            'xc': comp_c, 'xa': comp_a,
            'T_f': T_f, 't_initial': t_initial
        }
        sim_params = {'nodes': sim_nodes, 'dt': sim_dt, 'total_time_min': sim_time * 60}
        
        with st.spinner("Processando Matrizes Transientes via Numba/C++..."):
            try:
                # Dispara a simulação (Integração com simulation.py e physics.py)
                df_results = run_simulation(tunnel_params, box_params, comp_params, sim_params)
                
                # Transformação Melt para plotagem unificada no Plotly
                df_plot = df_results.to_pandas().melt(
                    id_vars=["Tempo (h)"], 
                    value_vars=["Superfície Topo", "Centro Geométrico", "Superfície Fundo", "Ponto Mais Quente"],
                    var_name="Sensor Virtual", 
                    value_name="Temperatura (°C)"
                )
                
                # Renderização Gráfica
                fig = px.line(df_plot, x="Tempo (h)", y="Temperatura (°C)", color="Sensor Virtual",
                              title="Perfil Térmico Assimétrico - Histórico de Temperaturas")
                
                # Linhas de Referência Físicas
                fig.add_hline(y=t_temp, line_dash="dash", line_color="red", opacity=0.5, annotation_text="Ar do Túnel")
                fig.add_hline(y=T_f, line_dash="dot", line_color="blue", opacity=0.5, annotation_text="Início Congelamento")
                fig.add_hline(y=-18.0, line_dash="dash", line_color="green", opacity=0.5, annotation_text="Meta Seg. Alimentar (-18°C)")
                
                fig.update_layout(hovermode="x unified", legend_title_text='Localização na Caixa')
                st.plotly_chart(fig, use_container_width=True)
                
            except Exception as e:
                import traceback
                st.error(f"Erro Crítico durante o cálculo da malha: {e}")
                st.code(traceback.format_exc())
