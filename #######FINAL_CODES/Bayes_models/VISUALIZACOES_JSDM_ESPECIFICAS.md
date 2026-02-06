# MASTER Viz Pipeline v3 - Visualizações Específicas por Tipo de Modelo

## Figuras Geradas por Tipo de Modelo

### ZOIB (Abundância) e Gaussian (Health/RGR)

| # | Figura | Descrição | Específico |
|---|--------|-----------|------------|
| 1 | `FIGURE_1_Forest_Plot.png` | Efeitos fixos com intervalos de credibilidade (half-eye) | Genérico |
| 2 | `FIGURE_2_Marginal_Effects.png` | Ribbon plots com interações ARCH quando aplicável | Genérico |
| 3 | `FIGURE_3_Validation_PPC.png` | Posterior predictive checks (dens_overlay + zero check para ZOIB) | Genérico |
| 4 | `FIGURE_4_Diagnostics.png` | Trace plots e autocorrelação | Genérico |

### JSDM (Comunidade) - Adicional

| # | Figura | Descrição | Específico JSDM |
|---|--------|-----------|----------------|
| 5 | `FIGURE_5_Species_Forest_Plot.png` | Forest plots por espécie com cores | ✓ |
| 6 | `FIGURE_6_Species_Comparison.png` | Violin plots comparando respostas entre espécies | ✓ |
| 7 | `FIGURE_7_PPC_by_Species.png` | PPC para cada espécie individualmente | ✓ |

---

## Visualizações JSDM Explicadas

### FIGURE 5: Species Forest Plot (JSDM)

**Propósito:** Mostrar efeitos fixos separados por espécie

```r
# Estrutura do plot:
# - Cada espécie tem sua própria cor (paleta Okabe-Ito)
# - Efeitos fixos plotados com half-eye (median + ICs)
# - Cores diferenciam M. hispida de outras espécies

Exemplo visual:
┌────────────────────────────────────────────────────────┐
│ M. hispida     ▓─────●────▓ (efeito positivo)           │
│ Favia gravida   ▓──────●────▓ (efeito neutro)           │
│ Porites spp.    ▓──────●────▓ (efeito negativo)        │
└────────────────────────────────────────────────────────┘
             -0.5   0   0.5   1.0   (effect size)
```

**Interpretação:** Permite identificar quais espécies são mais sensíveis aos gradientes ambientais.

---

### FIGURE 6: Species Comparison (JSDM)

**Propósito:** Comparar distribuições posteriores de efeitos entre espécies

```r
# Estrutura do plot:
# - Violin plots para cada espécie
# - Ponto preto = mediana
# - Linha tracejada em x=0

Exemplo visual:
┌────────────────────────────────────────────────────────┐
│ Porites spp.    ╭──────╮                            │
│ Favia gravida   ╰─────╯                            │
│ M. hispida      ╭──────╮         (distribuição         │
│                  ╰─────╯          diferente)           │
└────────────────────────────────────────────────────────┘
    -0.5   0   0.5   1.0   (effect size)
```

**Interpretação:** Sobreposição de distribuições mostra respostas diferenciais entre espécies.

---

### FIGURE 7: PPC by Species (JSDM)

**Propósito:** Validar ajuste do modelo para cada espécie individualmente

```r
# Estrutura do plot:
# - PPC density overlay para cada espécie (até 4)
# - Linha preta = dados observados
# - Área azul = predições do modelo

Exemplo visual:
┌────────────────────────────────────────────────────────┐
│ M. hispida                                                │
│ ╭───────────────────────────────────────╮               │
│ │ ▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓▓ (obs)         │
│ ╰───────────────────────────────────────╯               │
│            ╭──────╮ (simulated)                          │
│            ╰─────╯                                       │
├────────────────────────────────────────────────────────┤
│ Favia gravida                                             │
│ ...                                                      │
└────────────────────────────────────────────────────────┘
```

**Interpretação:** Verifica se o modelo captura a distribuição de cada espécie adequadamente.

---

## Comparação com Script JSDM Original

| Visualização | Script JSDM Original | Pipeline v3 | Status |
|---------------|---------------------|--------------|--------|
| Forest plot por espécie | ✓ | ✓ | Mantido |
| Species comparison (violin) | ✓ | ✓ | Mantido |
| Conditional effects | ✓ | ✓ | Incluído em Figure 2 |
| PPC por espécie | ✓ | ✓ | Mantido |
| Diagnostics (trace/ACF) | - | ✓ | Adicionado |

---

## Resumo: Por que Visualizações JSDM são Únicas

### Modelos de Espécie Única (ZOIB, Gaussian)
- **Uma resposta** (COVER_PROP, HEALTH_PC1, RGR)
- **Foco:** Como preditores afetam essa resposta
- **Visualização:** Efeitos de preditores na resposta única

### JSDM (Comunidade)
- **Múltiplas respostas** (proporções de várias espécies)
- **Foco:** Como preditores afetam **composição** e cada espécie responde
- **Visualizações adicionais necessárias:**
  1. **Separação por espécie** - Cores diferentes para cada espécie nos plots
  2. **Comparação interespecífica** - Violin plots mostrando distribuições lado a lado
  3. **Validação por espécie** - PPC individual para verificar ajuste

---

## Exemplo de Saída

### Modelo ZOIB/Gaussian (4 figuras)
```
WINNER_zoib_CV_02_Conservative/
├── FIGURE_1_Forest_Plot.png
├── FIGURE_2_Marginal_Effects.png
├── FIGURE_3_Validation_PPC.png
└── FIGURE_4_Diagnostics.png
```

### Modelo JSDM (7 figuras)
```
WINNER_jsdm_beta_CV_30/
├── FIGURE_1_Forest_Plot.png          (efeitos fixos gerais)
├── FIGURE_2_Marginal_Effects.png     (ribbon plots por preditor)
├── FIGURE_3_Validation_PPC.png       (validação global)
├── FIGURE_4_Diagnostics.png          (trace + ACF)
├── FIGURE_5_Species_Forest_Plot.png   (✓ JSDM-específico)
├── FIGURE_6_Species_Comparison.png     (✓ JSDM-específico)
└── FIGURE_7_PPC_by_Species.png        (✓ JSDM-específico)
```

---

## Formatação Consistente

Todas as figuras seguem o mesmo padrão profissional:
- **Tema:** `theme_publication()` (Nature/Science)
- **Cores:** Okabe-Ito (colorblind-friendly)
- **DPI:** 300
- **Fundo:** Branco
- **Fontes:** Sans-serif, tamanhos padronizados

---

## Script Atualizado

Arquivo: `#######FINAL_CODES/new_glm_approach/05_MASTER_Viz_Pipeline_v3.R`

**Funções JSDM adicionadas:**
- `generate_jsdm_species_forest_plot()`
- `generate_jsdm_species_response_comparison()`
- `generate_jsdm_ppc_by_species()`

**Integração:** Detecta automaticamente se o modelo é DIRICHLET e gera as figuras específicas automaticamente.
