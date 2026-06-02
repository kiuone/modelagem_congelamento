# modelagem congelamento
# 🧮 Fundamentos Matemáticos do Simulador Térmico

Este documento detalha o motor termodinâmico do projeto. A simulação da curva de congelação baseia-se num modelo de **Transferência de Calor Unidimensional (1D)**, resolvido numericamente através do **Método das Diferenças Finitas (MDF)** numa formulação explícita.

O objetivo do algoritmo é prever o arrefecimento desde a superfície até ao centro térmico do produto, lidando com o calor sensível, a extração do calor latente e a transição de fase gradual.

---

## 1. Discretização Espacial e Temporal

O domínio físico (a altura total do bloco de produto, $H$) é fatiado num número definido de nós ($N$). A distância entre cada nó é o passo espacial $\Delta x$:

$$\Delta x = \frac{H}{N - 1}$$

A simulação avança no tempo em pequenos incrementos, definidos pelo passo de tempo $\Delta t$ (em segundos). A cada incremento, a temperatura de cada fatia é recalculada.

---

## 2. Propriedades Termofísicas Variáveis (Choi & Okos)

A carne não possui propriedades térmicas constantes. À medida que a temperatura desce, a fração de água congela, alterando drasticamente a condutividade e a densidade. O modelo calcula propriedades dinâmicas a cada passo de tempo com base na composição centesimal (fração de água $X_w$, proteína $X_p$, gordura $X_f$, hidratos de carbono $X_c$ e cinzas $X_a$).

A propriedade efetiva global do nó é a soma ponderada das frações:

**Condutividade Térmica Efetiva ($k$):**
$$k_{eff} = k_w X_{w_{liq}} + k_{gelo} X_{gelo} + k_p X_p + k_f X_f + k_c X_c + k_a X_a$$

**Densidade Efetiva ($\rho$):**
Calculada pela conservação de volume:
$$\rho_{eff} = \frac{1}{\frac{X_{w_{liq}}}{\rho_w} + \frac{X_{gelo}}{\rho_{gelo}} + \frac{X_p}{\rho_p} + \frac{X_f}{\rho_f} + \frac{X_c}{\rho_c} + \frac{X_a}{\rho_a}}$$

---

## 3. Modelação da Mudança de Fase (Modelo de Chen)

O maior desafio matemático na congelação alimentar é o **calor latente**. O produto não congela a uma temperatura fixa (como a água pura a 0°C), mas sim de forma gradual a partir do seu ponto de congelação inicial ($T_f$). 

Para evitar descontinuidades matemáticas, o simulador utiliza o conceito de **Calor Específico Aparente ($c_{app}$)**, que engloba tanto o calor sensível como a energia latente de cristalização da água. 

Abaixo do ponto de congelação ($T < T_f$), aplica-se a equação de Chen:

$$c_{app} = \left( 1.55 + 1.26 X_s - \frac{(X_w - X_b) \cdot L_o \cdot T_f}{T^2} \right) \cdot 1000$$

Onde:
* $X_s$: Fração de sólidos totais.
* $X_b$: Fração de água fortemente ligada (que não congela).
* $L_o$: Calor latente de fusão da água ($333.6 \text{ kJ/kg}$).
* $T$: Temperatura atual do nó em °C.

---

## 4. O Motor de Diferenças Finitas (MDF)

A equação geral da condução de calor (Lei de Fourier) é discretizada para atualizar o estado de cada nó entre o instante atual ($n$) e o instante futuro ($n+1$).

### A. Nós Internos (Condução Pura)
O balanço de energia num nó interno ($i$) considera o fluxo de calor proveniente dos nós adjacentes ($i-1$ e $i+1$):

$$T_i^{n+1} = T_i^n + \frac{\Delta t}{\rho_i \cdot c_{app,i} \cdot \Delta x^2} \left[ k_{dir}(T_{i+1}^n - T_i^n) - k_{esq}(T_i^n - T_{i-1}^n) \right]$$

Sendo $k_{dir}$ e $k_{esq}$ as condutividades médias harmónicas nas interfaces entre os nós.

### B. Nós de Fronteira (Convecção + Condução com Assimetria)
As superfícies enfrentam resistências térmicas diferentes (ex: o topo protegido por um filme plástico fino, enquanto a base está protegida pelo filme plástico interfolhado mais a espessura da caixa de papelão).

O coeficiente global de transferência de calor ($U$) é calculado para cada lado:
$$U = \frac{1}{\frac{1}{h_{conv}} + \sum \frac{e_{emb}}{k_{emb}}}$$

O balanço no nó da superfície topo ($i = 0$), exposto ao ar frio do túnel ($T_{\infty}$), é governado por:
$$q_{topo} = U_{topo} \cdot (T_{\infty} - T_0^n)$$
$$T_0^{n+1} = T_0^n + \frac{\Delta t}{\rho_0 \cdot c_{app,0} \cdot \Delta x} \left[ q_{topo} + k_0 \frac{T_1^n - T_0^n}{\Delta x} \right]$$

*(A mesma lógica análoga é aplicada ao nó do fundo da embalagem $i = N-1$ com $U_{base}$)*.

---

## 5. Critério de Estabilidade Numérica (Número de Fourier)

Por se tratar de um método numérico explícito, o passo de tempo ($\Delta t$) não pode ser arbitrariamente grande, caso contrário a solução matemática diverge (instabilidade térmica). O algoritmo monitoriza continuamente o **Número de Fourier ($Fo$)**:

$$Fo = \frac{\alpha_{max} \cdot \Delta t}{\Delta x^2} \le 0.5$$

Onde $\alpha_{max}$ é a máxima difusividade térmica ($\frac{k}{\rho \cdot c_p}$) alcançada pela matriz biológica na fase congelada. O sistema está desenhado para interromper a execução e gerar um alerta caso o critério $Fo > 0.5$ seja violado pelas configurações introduzidas pelo utilizador.