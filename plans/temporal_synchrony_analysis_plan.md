# Plano de Implementação Atualizado: Análise de Sincronicidade e Dinâmica (Nature/Science)
## Consolidação da Estratégia de Análise Temporal

**Data:** 2026-01-12
**Status:** Atualizado para incluir Decomposição STL, Análise de Lag e Visualização com 3 Eixos
**Script Alvo:** `#######FINAL_CODES/TIME_SERIES_SYNCHRONY_ANALYSIS.py`

---

## 1. Visão Geral da Atualização

Este plano consolida a estratégia para elevar a análise de sincronicidade ao nível "Nature/Science". A principal mudança visual é a adoção de **3 eixos Y independentes** para as séries temporais, evitando o "achatamento" visual causado por escalas díspares (SST ~25-28°C, DLI ~10-40, Chl-a ~0.5-5.0). A análise estatística será reforçada com **Decomposição STL** para isolar pulsos de sedimentos.

---

## 2. Nova Estrutura de Visualização (3 Eixos Y)

Para evitar que a Clorofila (valores baixos) fique invisível ao lado do DLI (valores altos), usaremos a técnica de "Parasite Axes" do Matplotlib.

### Snippet de Implementação (Para o Script)

```python
def make_patch_spines_invisible(ax):
    ax.set_frame_on(True)
    ax.patch.set_visible(False)
    for sp in ax.spines.values():
        sp.set_visible(False)

def plot_three_axis_timeseries(ax, dates, sst, dli, chl, title):
    """
    Plota SST, DLI e Chl-a em 3 eixos Y independentes com cores distintas.
    """
    # Eixo 1: SST (Esquerda - Vermelho)
    ax.plot(dates, sst, color=COLOR_SST, linewidth=1.5, label='SST', zorder=2)
    ax.set_ylabel('SST (°C)', color=COLOR_SST, fontsize=12, fontweight='bold')
    ax.tick_params(axis='y', labelcolor=COLOR_SST)
    ax.yaxis.label.set_color(COLOR_SST)
    
    # Eixo 2: DLI (Direita 1 - Azul)
    ax2 = ax.twinx()
    ax2.plot(dates, dli, color=COLOR_DLI, linewidth=1.5, label='DLI', zorder=3)
    ax2.set_ylabel('DLI (mol m⁻² d⁻¹)', color=COLOR_DLI, fontsize=12, fontweight='bold')
    ax2.tick_params(axis='y', labelcolor=COLOR_DLI)
    ax2.spines['right'].set_color(COLOR_DLI)
    
    # Eixo 3: Chl-a (Direita 2 - Verde - Deslocado)
    ax3 = ax.twinx()
    # Deslocar eixo para a direita
    ax3.spines["right"].set_position(("axes", 1.15))
    make_patch_spines_invisible(ax3)
    ax3.spines["right"].set_visible(True)
    
    ax3.plot(dates, chl, color=COLOR_CHL, linewidth=1.5, label='Chl-a', zorder=4)
    ax3.set_ylabel('Chl-a (mg m⁻³)', color=COLOR_CHL, fontsize=12, fontweight='bold')
    ax3.tick_params(axis='y', labelcolor=COLOR_CHL)
    ax3.spines['right'].set_color(COLOR_CHL)

    ax.set_title(title, fontsize=14, fontweight='bold', loc='left', pad=15)
    return [ax, ax2, ax3]
```

---

## 3. Análise de Sincronicidade de Pulsos (STL Decomposition)

Para provar que a correlação não é apenas sazonal (inverno vs verão), isolaremos os "pulsos" (resíduos).

### Snippet de Implementação

```python
from statsmodels.tsa.seasonal import seasonal_decompose

def calculate_pulse_synchrony(sst_series, dli_series, chl_series, period=30):
    """
    Decompõe as séries e calcula correlação dos RESÍDUOS (Pulsos).
    """
    # Decompor (Model='additive' geralmente seguro para logs ou dados normalizados, 
    # mas 'multiplicative' pode ser melhor para Chl-a. Vamos usar additive após z-score)
    
    # Garantir sem NaNs e frequência definida (necessário para statsmodels)
    df = pd.DataFrame({'SST': sst_series, 'DLI': dli_series, 'CHL': chl_series}).dropna()
    
    # Decomposição
    res_sst = seasonal_decompose(df['SST'], period=period, extrapolate_trend='freq').resid
    res_dli = seasonal_decompose(df['DLI'], period=period, extrapolate_trend='freq').resid
    res_chl = seasonal_decompose(df['CHL'], period=period, extrapolate_trend='freq').resid
    
    # Correlações de Pulso
    pulse_corr_dli_chl = res_dli.corr(res_chl) # Esperado negativo (muita luz = menos chl? ou oposto?)
    # Nota: Em turbidez, sedimento bloqueia luz (DLI cai) e traz nutrientes (Chl sobe) -> Correlação Negativa
    
    pulse_corr_sst_chl = res_sst.corr(res_chl)
    
    return {
        'Pulse_Corr_DLI_CHL': pulse_corr_dli_chl,
        'Pulse_Corr_SST_CHL': pulse_corr_sst_chl,
        'Resid_DLI': res_dli,
        'Resid_CHL': res_chl
    }
```

---

## 4. Atualização da Figura Principal (Layout 2x3)

A figura deve ser expandida para acomodar a nova profundidade de análise.

| Coluna A (Inner Arc) | Coluna B (Outer Arc) | Coluna C (Estatística) |
| :--- | :--- | :--- |
| **Painel A:** Time Series (3 Eixos)<br>*SST, DLI, Chl (Raw + Smooth)* | **Painel B:** Time Series (3 Eixos)<br>*SST, DLI, Chl (Raw + Smooth)* | **Painel C:** Pearson Geral<br>*(Boxplot Inner vs Outer)* |
| **Painel D:** Cross-Correlation (Lag)<br>*Lag 0 a ±15 dias* | **Painel E:** Cross-Correlation (Lag)<br>*Lag 0 a ±15 dias* | **Painel F:** **Pulse Synchrony (STL)**<br>*(Boxplot dos Resíduos)* |

### Detalhes de Implementação da Figura

1.  **GridSpec:** `gs = fig.add_gridspec(2, 3, width_ratios=[1.5, 1.5, 1])` para dar mais espaço às séries temporais.
2.  **Painéis D e E (CCF):** Usar `plt.xcorr` ou cálculo manual de correlação cruzada para mostrar se o pico de correlação está em Lag 0 (imediato) ou Lag positivo (atrasado).
    *   *Hipótese:* Inner Arc tem pico alto e Lag curto. Outer Arc tem pico baixo e difuso.

---

## 5. Roteiro Consolidado de Execução

1.  **Setup & Load:**
    *   Carregar sites e dados de satélite (código existente).
    *   Agregar por arco (Inner/Outer).

2.  **Processamento Avançado (Novo):**
    *   Iterar sobre cada SITE.
    *   Calcular `Pearson_Geral` (código existente).
    *   Executar `seasonal_decompose` para obter Resíduos.
    *   Calcular `Pearson_Pulsos` (Resíduo DLI vs Resíduo Chl).
    *   Armazenar métricas em DataFrame mestre.

3.  **Visualização (Novo):**
    *   Gerar **Figura Principal** usando a função `plot_three_axis_timeseries` para os Painéis A e B.
    *   Gerar plots de CCF para Painéis D e E.
    *   Gerar Boxplots comparativos para Painéis C e F.

4.  **Estatística:**
    *   Comparar `Pearson_Pulsos` (Inner vs Outer) via Mann-Whitney U.

---

## 6. Próximos Passos (Para o Agente)

Ao receber o comando para executar, o agente deverá:
1.  Atualizar o script `TIME_SERIES_SYNCHRONY_ANALYSIS.py` com as novas funções de plotagem (3 eixos) e análise STL.
2.  Garantir que as bibliotecas `statsmodels` estejam disponíveis (ou adicionar bloco try/except com fallback).
3.  Executar o script e validar as figuras geradas.

Este plano garante que a complexidade visual (escalas diferentes) seja resolvida e que a complexidade científica (sazonalidade vs pulsos) seja abordada diretamente.
