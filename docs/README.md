# Documentação do Projeto - Correções JSDM

Esta pasta contém a documentação completa das correções implementadas para os modelos JSDM (Joint Species Distribution Models) com distribuição Dirichlet.

## 📚 Documentos Disponíveis

### Guia Principal
- **[JSDM_Corrections_Guide.md](JSDM_Corrections_Guide.md)** - Guia completo com:
  - Contexto dos problemas
  - Soluções implementadas
  - Passo a passo para execução
  - Explicação técnica detalhada
  - Solução de problemas

## 🚀 Quick Start

### Opção 1: Pipeline Completo (Recomendado)
```r
source("#######FINAL_CODES/MASTER_Viz_Pipeline_v5_JSDM_FIXED.R")
```

### Opção 2: Demonstração
```r
source("#######FINAL_CODES/JSDM_Demo_and_Validation.R")
```

### Opção 3: Funções Específicas
```r
source("#######FINAL_CODES/JSDM_Correction_Module.R")
model <- readRDS("caminho/para/WINNER_JSDM_full_CV_ALL.rds")
p1 <- generate_jsdm_forest_plot_errorbar_CORRECTED(model, "JSDM")
p2 <- generate_jsdm_pdp_community_overlaid_CORRECTED(model, "JSDM")
```

## 📁 Localização dos Scripts

Todos os scripts de correção estão em:
```
#######FINAL_CODES/
├── JSDM_Fix_Mhispida_and_PDPs.R
├── JSDM_Correction_Module.R
├── MASTER_Viz_Pipeline_v5_JSDM_FIXED.R
└── JSDM_Demo_and_Validation.R
```

## ✅ Resumo das Correções

| Problema | Solução | Resultado |
|----------|---------|-----------|
| M. hispida não aparece no Forest Plot | Adicionar entradas dummy com effect = 0 | M. hispida visível como referência |
| PDPs mostram linha única | Usar `posterior_epred()` diretamente | 5 curvas separadas por espécie |

## 📞 Suporte

Para dúvidas técnicas, consulte o [JSDM_Corrections_Guide.md](JSDM_Corrections_Guide.md) completo.
