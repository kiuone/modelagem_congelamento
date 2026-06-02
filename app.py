import streamlit as st
import plotly.express as px
from simulation import run_simulation

st.set_page_config(page_title="Simulador de Túnel de Congelamento", layout="wide")

st.title("Simulador de congelamento industrial")
st.markdown("Modelo 1D por Diferenças Finitas - Considerando Assimetria de Embalagens")

# --- BARRA LATERAL PARA CONTROLES ---
st.sidebar.header("Parâmetros do Processo")

# 1. Configurações do Túnel
with st.sidebar.expander("🌬️ Propriedades do Túnel", expanded=True):
    tunel_opcoes = {
        "Túnel 1": {"temp": -40.0, "h": 30.0},
        "Túnel 2": {"temp": -33.0, "h": 30.0},
        "Túnel 3": {"temp": -33.0, "h": 30.0}
    }
    tunel_selecionado = st.selectbox("Selecione o Túnel", list(tunel_opcoes.keys()))
    
    t_temp = st.number_input("Temperatura do Ar (°C)", value=tunel_opcoes[tunel_selecionado]["temp"], step=1.0)
    t_h = st.number_input("Coeficiente Convectivo (W/m².K)", value=tunel_opcoes[tunel_selecionado]["h"], step=1.0)

# 2. Configurações da Caixa e Embalagem
with st.sidebar.expander("📦 Propriedades da Caixa", expanded=False):
    # Dicionário com as regras de negócio e dimensões reais de cada modelo
    modelos_caixa = {
        "Fundo Baixo": {"altura": 118.8, "espessura": 3.8},
        "Fundo Alto": {"altura": 150.0, "espessura": 3.8},
        "Branca": {"altura": 100.0, "espessura": 4.0},      
        "Customizada": {"altura": 118.8, "espessura": 3.8}  }
    
    caixa_selecionada = st.selectbox("Selecione o Modelo da Caixa", list(modelos_caixa.keys()))
    
    # Se for customizada, permite editar. Se for fixa, apenas exibe e trava os valores.
    if caixa_selecionada == "Customizada":
        c_altura = st.number_input("Altura do Produto na Caixa (mm)", value=modelos_caixa[caixa_selecionada]["altura"], step=1.0)
        p_espessura = st.number_input("Espessura Papelão (mm)", value=modelos_caixa[caixa_selecionada]["espessura"], step=0.1)
    else:
        c_altura = modelos_caixa[caixa_selecionada]["altura"]
        p_espessura = modelos_caixa[caixa_selecionada]["espessura"]
        
        # Mostra as medidas travadas embaixo do seletor
        st.markdown(f"**Medidas Ativas para {caixa_selecionada}:**")
        st.caption(f"Altura do bloco de carne: {c_altura} mm")
        st.caption(f"Espessura do papelão: {p_espessura} mm")

    st.markdown("---")
    st.markdown("**Configurações Adicionais de Condutividade**")
    p_k = st.number_input("Condutividade Papelão/plástico branco (W/m.K)", value=0.064, step=0.01)
    
    st.markdown("**Filme Plástico (Interfolhado)**")
    f_espessura = st.number_input("Espessura Plástico (mm)", value=0.1, step=0.05)
    f_k = st.number_input("Condutividade Plástico (W/m.K)", value=0.33, step=0.05)

# 3. Composição do Produto
with st.sidebar.expander("🍗 Matriz Biológica (Produto)", expanded=False):
    st.markdown("Composição baseada em dados típicos do peito de frango (soma deve ser 1.0). Altere conforme o tipo de corte.")
    comp_w = st.slider("Água (Fração)", 0.0, 1.0, 0.748)
    comp_p = st.slider("Proteína (Fração)", 0.0, 1.0, 0.224)
    comp_f = st.slider("Gordura (Fração)", 0.0, 1.0, 0.011)
    comp_c = st.slider("Carboidrato (Fração)", 0.0, 1.0, 0.0)
    comp_a = st.slider("Cinzas (Fração)", 0.0, 1.0, 0.017)
    
    t_initial = st.number_input("Temperatura de Entrada (°C)", value=4.0, step=1.0)
    t_f = st.number_input("Ponto de Congelamento Inicial (°C)", value=-2.8, step=0.1)

# 4. Parâmetros Numéricos
with st.sidebar.expander("⚙️ Motor de Simulação", expanded=False):
    sim_time = st.number_input("Tempo Total Simulador (horas)", value=24.0, step=1.0)
    sim_nodes = st.number_input("Número de Nós da Malha", min_value=5.0, value=6.0, step=1.0)
    sim_dt = st.number_input("Passo de Tempo dt (s)", value=0.1, step=0.1)

    # MONITOR DE ESTABILIDADE (CRITÉRIO DE FOURIER)
    # dx é a espessura da fatia em metros (Altura do produto / (Nós - 1))
    dx_m = (c_altura / 1000.0) / (sim_nodes - 1) 
    alpha_max = 1.4e-6  # Difusividade térmica max aproximada (m²/s)
    
    # Cálculo do Número de Fourier atual
    Fo = (alpha_max * sim_dt) / (dx_m ** 2)
    
    st.markdown("---")
    st.markdown("**Monitor de Estabilidade (Fourier)**")
    
    # Lógica de exibição visual baseada no limite de 0.5
    if Fo <= 0.5:
        st.success(f"✅ **Fo = {Fo:.3f}** (Estável)")
        st.caption("*O limite máximo matemático é de 0.5, mas quanto menor o Fo, melhor.*")
    else:
        st.error(f"🚨 **Fo = {Fo:.3f}** (Instável!)")
        st.caption("*O modelo vai divergir. Reduza o 'dt' ou o 'Número de Nós' até o Fo ficar <= 0.5.*")

# --- EXECUÇÃO E GRÁFICOS ---
if st.sidebar.button("▶️ Rodar Simulação Térmica", width='stretch'):
    
    # Empacotando parâmetros
    tunnel_params = {'temp': t_temp, 'h': t_h}
    box_params = {
        'height': c_altura, 
        'carton_thickness': p_espessura, 'carton_k': p_k,
        'plastic_thickness': f_espessura, 'plastic_k': f_k
    }
    comp_params = {
        'water': comp_w, 'protein': comp_p, 'fat': comp_f, 
        'carbohydrate': comp_c, 'ash': comp_a,
        't_initial': t_initial, 't_f': t_f
    }
    sim_params = {'nodes': sim_nodes, 'dt': sim_dt, 'total_time_min': sim_time * 60}
    
    # Soma de segurança
    total_comp = comp_w + comp_p + comp_f + comp_c + comp_a
    if abs(total_comp - 1.0) > 0.01:
        st.error(f"Atenção: A soma das frações de composição é {total_comp:.2f}. O ideal é 1.0.")
    else:
        with st.spinner("Integrando modelo numérico e atualizando propriedades termofísicas..."):
            try:
                # O retorno é um polars DataFrame
                df_results = run_simulation(tunnel_params, box_params, comp_params, sim_params)
                
                # Conversão para pandas para injeção no Plotly express
                df_plot = df_results.to_pandas().melt(
                    id_vars=["Tempo (min)"], 
                    value_vars=["Superfície Topo (°C)", "Centro Térmico (°C)", "Superfície Fundo (°C)"],
                    var_name="Localização", 
                    value_name="Temperatura (°C)"
                )
                
                fig = px.line(
                    df_plot, x="Tempo (min)", y="Temperatura (°C)", color="Localização",
                    title="Curvas de Congelamento do Produto na Caixa"
                )
                # Linha da temperatura do túnel para referência
                fig.add_hline(y=t_temp, line_dash="dash", line_color="red", annotation_text="Temp. Túnel")
                # Linha de congelamento
                fig.add_hline(y=t_f, line_dash="dot", line_color="blue", annotation_text="Ponto de Transição")
                # Linhas de referência de temperatura (-12°C e -18°C)
                fig.add_hline(y=-12, line_dash="dash", line_color="orange", annotation_text="-12°C")
                fig.add_hline(y=-18, line_dash="dash", line_color="green", annotation_text="-18°C")
                
                st.plotly_chart(fig, width='stretch')
                
                st.success("Simulação concluída com sucesso! Os resultados convergem as resistências térmicas isoladas e o calor latente.")
                
            except Exception as e:
                st.error(f"Erro na simulação: {e}")
else:
    st.info("Ajuste os parâmetros na barra lateral e clique em 'Rodar Simulação Térmica' para ver a curva.")