# ================================================FILE: app.py
import streamlit as st
import plotly.express as px
import polars as pl
from simulation import run_simulation

# Configuração global da página
st.set_page_config(page_title="Simulador de Congelamento", layout="wide")
st.title("Simulador de Congelamento Industrial")
st.markdown("Modelo 1D por Diferenças Finitas - Otimizado para Caixas Destampadas")

# ==============================================================================
# BARRA LATERAL: ENTRADA DE PARÂMETROS
# ==============================================================================
st.sidebar.header("Parâmetros de Controle")

# --- 1. TÚNEL DE CONGELAMENTO ---
with st.sidebar.expander("🌬️ Condições do Túnel", expanded=True):
    tunel_opcoes = {
        "Customizado (6 m/s)": {"temp": -31.0, "h": 60.0},
        "Blast Freezer (Alta Vel.)": {"temp": -33.0, "h": 80.0},
        "Túnel Convencional": {"temp": -30.0, "h": 40.0},
        "Câmara Estática (S/ Vento)": {"temp": -22.0, "h": 15.0}
    }
    tunel_selecionado = st.selectbox("Perfil Operacional", list(tunel_opcoes.keys()))
    
    # Variáveis de estado do túnel
    t_temp = st.number_input("Temperatura do Ar (°C)", value=tunel_opcoes[tunel_selecionado]["temp"], step=1.0)
    t_h = st.number_input("Coef. Convectivo 'h' (W/m².K)", value=tunel_opcoes[tunel_selecionado]["h"], step=5.0, min_value=5.0, max_value=300.0)
    # st.caption("*Nota: 6 m/s equivale a aproximadamente h = 60 W/m²K pela fórmula empírica (12 + 8v).*")

# --- 2. GEOMETRIA E EMBALAGEM ---
with st.sidebar.expander("📦 Propriedades da Caixa", expanded=False):
    # Dimensões globais para avaliação de dimensionalidade (Fator E) e malha 1D
    c_altura = st.number_input("Altura do Bloco de Carne (mm)", value=118.8, step=1.0)
    c_largura = st.number_input("Largura da Caixa (mm)", value=385.6, step=1.0)
    c_comprimento = st.number_input("Comprimento da Caixa (mm)", value=575.6, step=1.0)
    
    st.markdown("---")
    st.markdown("**Isolamento Inferior e Lateral (Papelão)**")
    p_espessura = st.number_input("Espessura Papelão (mm)", value=3.8, step=0.1)
    p_k = st.number_input("Condutividade Papelão (W/m.K)", value=0.12, step=0.01)
    
    st.markdown("**Isolamento de Contato (Filme Plástico)**")
    f_espessura = st.number_input("Espessura Plástico (mm)", value=0.1, step=0.05)
    f_k = st.number_input("Condutividade Plástico (W/m.K)", value=0.33, step=0.05)

# --- 3. MATRIZ BIOLÓGICA ---
with st.sidebar.expander("🍗 Composição do Produto", expanded=False):
    st.markdown("Frações Mássicas (A soma ideal é 1.0)")
    comp_w = st.slider("Água", 0.0, 1.0, 0.669)
    comp_p = st.slider("Proteína", 0.0, 1.0, 0.301)
    comp_f = st.slider("Gordura", 0.0, 1.0, 0.020)
    comp_c = st.slider("Carboidrato", 0.0, 1.0, 0.000)
    comp_a = st.slider("Cinzas", 0.0, 1.0, 0.010)
    
    soma_comp = comp_w + comp_p + comp_f + comp_c + comp_a
    
    st.markdown("---")
    t_initial = st.number_input("Temperatura de Entrada (°C)", value=3.1, step=1.0)
    T_f = st.number_input("Ponto de Congelamento Inicial (°C)", value=-2.0, step=0.1)

# --- 4. PARÂMETROS NUMÉRICOS ---
with st.sidebar.expander("⚙️ Motor Numérico (Diferenças Finitas)", expanded=False):
    sim_time = st.number_input("Tempo de Simulação (horas)", value=24.0, step=1.0)
    sim_nodes = st.number_input("Número de Nós (Fatias na Altura)", min_value=5, value=120, step=10)
    sim_dt = st.number_input("Passo de Tempo 'dt' (segundos)", value=0.01, step=0.01, format="%.3f")


# ==============================================================================
# PAINEL DE DIAGNÓSTICO E INDICADORES (PRÉ-SIMULAÇÃO)
# ==============================================================================
st.markdown("### 📊 Diagnóstico Termodinâmico e Numérico")

if abs(soma_comp - 1.0) > 0.01:
    st.warning(f"⚠️ **Atenção:** A soma da composição centesimal é {soma_comp:.3f}. Para precisão termodinâmica máxima, ajuste os sliders para que a soma seja exata a 1.0.")

col1, col2, col3, col4 = st.columns(4)

dx_m = (c_altura / 1000.0) / (sim_nodes - 1)
alpha_max = 1.4e-6 
Fo = (alpha_max * sim_dt) / (dx_m ** 2)

with col1:
    st.metric("Número de Fourier (Fo)", f"{Fo:.4f}")
    if Fo <= 0.5:
        st.success("Malha Estável")
    else:
        st.error("Malha Instável (Reduza o dt!)")

E_factor = 1.0 + (c_altura / c_largura)**2 + (c_altura / c_comprimento)**2
with col2:
    st.metric("Fator de Forma (E)", f"{E_factor:.3f}")
    st.caption(f"Aceleração 3D estimada em +{(E_factor - 1) * 100:.1f}%")

R_plastico = (f_espessura / 1000.0) / f_k
R_papelao = (p_espessura / 1000.0) / p_k

U_top_UI = 1.0 / ((1.0 / t_h) + R_plastico)
U_bottom_UI = (1.0 / ((1.0 / t_h) + R_plastico + R_papelao)) * 0.6

k_frango_frozen = 1.35 
L_c = (c_altura / 1000.0) / 2.0 

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
if st.button("▶️ Iniciar Integração Numérica", type="primary", width='stretch'):
    
    if Fo > 0.5:
        st.error("🚨 Simulação Abortada. O Critério de Fourier foi violado. Reduza o 'dt' ou diminua a quantidade de Nós.")
    else:
        tunnel_params = {'temp': t_temp, 'h': t_h}
        box_params = {
            'altura': c_altura, 
            'p_espessura': p_espessura, 'p_k': p_k,
            'f_espessura': f_espessura, 'f_k': f_k
        }
        comp_params = {
            'xw': comp_w, 'xp': comp_p, 'xf': comp_f, 
            'xc': comp_c, 'xa': comp_a,
            'T_f': T_f, 't_initial': t_initial}
        sim_params = {'nodes': sim_nodes, 'dt': sim_dt, 'total_time_min': sim_time * 60}
        
        with st.spinner("Processando Matrizes Transientes via Numba/C++..."):
            try:
                # Dispara a simulação
                df_results = run_simulation(tunnel_params, box_params, comp_params, sim_params)
                
                # Divisão da tela: Gráfico na esquerda (75%) e Relatório na direita (25%)
                col_grafico, col_relatorio = st.columns([3, 1])
                
                # --- LADO ESQUERDO: GRÁFICO ---
                with col_grafico:
                    df_plot = df_results.to_pandas().melt(
                        id_vars=["Tempo (h)"], 
                        value_vars=["Superfície Topo", "Centro Geométrico", "Superfície Fundo", "Ponto Mais Quente"],
                        var_name="Sensor Virtual", 
                        value_name="Temperatura (°C)")
                    
                    fig = px.line(df_plot, x="Tempo (h)", y="Temperatura (°C)", color="Sensor Virtual",
                                  title="Perfil Térmico Assimétrico - Histórico de Temperaturas")
                    
                    fig.add_hline(y=t_temp, line_dash="dash", line_color="red", opacity=0.5, annotation_text="Ar do Túnel")
                    fig.add_hline(y=T_f, line_dash="dot", line_color="blue", opacity=0.5, annotation_text="Início Congelamento")
                    fig.add_hline(y=-18.0, line_dash="dash", line_color="green", opacity=0.5, annotation_text="(-18°C)")
                    fig.update_layout(hovermode="closest", legend_title_text='Localização')
                    st.plotly_chart(fig, width='stretch')
                
                # --- LADO DIREITO: PAINEL DE TEMPO ATÉ -18°C ---
                with col_relatorio:
                    st.markdown("### 🎯 Meta -18°C")
                    st.markdown("Tempo necessário para cada ponto atingir o alvo térmico:")
                    
                    sensores = ["Superfície Topo", "Superfície Fundo", "Centro Geométrico", "Ponto Mais Quente"]
                    
                    for sensor in sensores:
                        # Filtra no Polars onde a temperatura caiu de -18.0
                        df_filtro = df_results.filter(pl.col(sensor) <= -18.0)
                        
                        if len(df_filtro) > 0:
                            # Pega a primeira linha que cruzou -18°C
                            tempo_h = df_filtro.get_column("Tempo (h)")[0]
                            tempo_m = int(tempo_h * 60)
                            st.info(f"**{sensor}**\n\n⏱️ {tempo_h:.2f} h ({tempo_m} min)")
                        else:
                            st.warning(f"**{sensor}**\n\n⏳ Não atingiu o alvo.")
                            
            except Exception as e:
                st.error(f"Erro Crítico durante o cálculo da malha: {e}")