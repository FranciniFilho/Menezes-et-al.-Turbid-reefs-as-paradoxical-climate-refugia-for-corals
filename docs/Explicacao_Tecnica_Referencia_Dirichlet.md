# Por que M. hispida não tem coeficientes estimados?

## Explicação Técnica Completa

### A Pergunta
> "Por que M. hispida é a categoria de referência no modelo Dirichlet e não tem coeficientes estimados (log-ratio relativo a si mesmo = 0)?"

### Resposta Curta
Porque o **modelo Dirichlet é uma composição** (proporções que somam 1), e para estimar efeitos relativos entre categorias, precisamos de uma **categoria de referência** (baseline) contra a qual todas as outras são comparadas. Matematicamente, o log-ratio de algo relativo a si mesmo é sempre **log(1) = 0**.

---

## 1. O Conceito de Composição

### Dados Composicionais
No seu projeto, cada observação é uma **composição**:

```
MUSSISMILIA + TURF + CCA + CYANO + MACROALGAE = 100% (ou 1.0)
```

Exemplo de uma observação:
```
M. hispida:  0.25 (25%)
Turf:        0.30 (30%)
CCA:         0.20 (20%)
Cyano:       0.15 (15%)
Macroalgae:  0.10 (10%)
            ─────
Total:       1.00 (100%)
```

### O Problema da Dependência
Se eu sei que 4 categorias somam 70%, a 5ª **necessariamente** é 30%. Isso cria **dependência** entre as variáveis.

**Implicação:** Não podemos modelar cada proporção independentemente, pois elas estão matematicamente vinculadas.

---

## 2. Parametrização do Modelo Dirichlet

### A Transformação Log-Ratio

O modelo Dirichlet no `brms` usa uma **parametrização de log-ratio** (razão logarítmica) para lidar com composições:

```
Para cada observação i e categoria k (exceto a referência):

log(α_ik / α_i1) = η_ik

Onde:
• α_ik = parâmetro de concentração da categoria k
• α_i1 = parâmetro de concentração da categoria de referência (M. hispida)
• η_ik = preditor linear para categoria k
```

### Por que dividir pela referência?

Pense em **razões** (ratios):

| Comparação | Razão | Log-ratio |
|------------|-------|-----------|
| TURF / M. hispida | 1.2 | log(1.2) = 0.18 |
| CCA / M. hispida | 0.8 | log(0.8) = -0.22 |
| **M. hispida / M. hispida** | **1.0** | **log(1.0) = 0** |

**Quando compar algo consigo mesmo, a razão é 1, e log(1) = 0.**

---

## 3. O Problema de Identificabilidade

### O que é Identificabilidade?

Um modelo é **identificável** quando existe uma única solução (único conjunto de parâmetros) que maximiza a verossimilhança.

### Por que precisamos de uma referência?

Imagine que tentamos estimar efeitos para TODAS as 5 espécies sem referência:

```
log(α_TURF) = β₀_TURF + β₁_TURF·PC1 + ...
log(α_CCA)  = β₀_CCA  + β₁_CCA·PC1 + ...
...
log(α_MUSS) = β₀_MUSS + β₁_MUSS·PC1 + ...
```

**O problema:** Se eu adicionar 1 a TODOS os β₀, as razões **não mudam**:

```
(α_TURF · e¹) / (α_MUSS · e¹) = α_TURF / α_MUSS
```

Isso significa que existem **infinitas** combinações de parâmetros que produzem o mesmo ajuste!

### A Solução: Fixar uma categoria

Ao definir **M. hispida como referência**:

```
log(α_MUSS / α_MUSS) = log(1) = 0
```

Fixamos o sistema. Agora todos os outros parâmetros são **relativos a M. hispida**:

```
log(α_TURF / α_MUSS) = β₀_TURF + β₁_TURF·PC1 + ...
log(α_CCA / α_MUSS)  = β₀_CCA  + β₁_CCA·PC1 + ...
...
```

**Agora o modelo é identificável!**

---

## 4. Analogia Intuitiva

### Analogia: Alturas Relativas

Imagine que você quer medir as alturas de 5 pessoas **relativas entre si**:

**Abordagem 1: Sem referência**
- Pessoa A: 170 cm
- Pessoa B: 180 cm
- Pessoa C: 160 cm
- ...

Mas espere - **eu poderia adicionar 10 cm a TODOS** e as diferenças relativas seriam as mesmas:
- A: 180, B: 190, C: 170...
- Diferença A-B ainda é 10 cm!

**Sem uma referência fixa, as alturas absolutas são arbitrárias.**

**Abordagem 2: Com referência**
- Definir Pessoa A como referência (altura = 0)
- Pessoa B: +10 cm (relativo a A)
- Pessoa C: -10 cm (relativo a A)

**Agora as medidas são únicas e interpretáveis!**

### Analogia: Probabilidades Logit

No modelo logístico binário, temos:

```
log(p / (1-p)) = β₀ + β₁·X
```

Aqui, **(1-p) é a referência implícita**. Não estimamos coeficientes para "não-evento" - é o baseline contra o qual o evento é comparado.

O modelo Dirichlet é a **extensão multivariada** dessa ideia.

---

## 5. Por que M. hispida foi escolhida como referência?

### A Resposta está na Formulação

```r
Formula: cbind(MUSSISMILIA_prop, TURF_prop, CCA_prop, CYANO_prop, MACROALGAE_prop) ~ ...
```

No `brms`, para família Dirichlet:
- **A primeira coluna no cbind() é automaticamente a referência**
- Isso é uma convenção do pacote

### Você poderia mudar?

Sim! Se a fórmula fosse:

```r
cbind(TURF_prop, MUSSISMILIA_prop, CCA_prop, CYANO_prop, MACROALGAE_prop) ~ ...
```

Agora **TURF seria a referência**, e todos os coeficientes seriam:
- `muMUSSISMILIAprop_*`
- `muCCAprop_*`
- etc.

**Mas não teríamos `muTURFprop_*`**

### Qual a melhor escolha?

A escolha da referência é **estatisticamente arbitrária** (não afeta o ajuste), mas **biologicamente importante**:

✅ **Boa escolha (M. hispida):**
- É a espécie de interesse principal do estudo
- Interpretação natural: "Efeito do ambiente nas outras espécies RELATIVO à dominância de M. hispida"
- Facilita a narrativa científica

---

## 6. Interpretação dos Coeficientes

### O que significa cada coeficiente?

Dado o coeficiente `muTURFprop_PC1MAGNITUDE = -0.45`:

```
Significado: Um aumento de 1 unidade em PC1MAGNITUDE está associado a
uma mudança de -0.45 no log-ratio de TURF relativo a M. hispida.

Matematicamente:
log(Proporção_TURF / Proporção_MUSS) diminui em 0.45

Interpretação biológica:
Quando PC1MAGNITUDE aumenta, TURF diminui RELATIVAMENTE a M. hispida.
(ou seja, M. hispida "ganha" espaço do TURF nesse gradiente)
```

### Efeito positivo vs negativo

| Coeficiente | Interpretação |
|-------------|---------------|
| Positivo | A espécie AUMENTA relativo a M. hispida |
| Negativo | A espécie DIMINUI relativo a M. hispida |
| Zero | A espécie se mantém estável relativo a M. hispida |

---

## 7. Visualização: Por que M. hispida deve aparecer em 0

### No Forest Plot

```
│ MUSSISMILIA (ref) ♦─────── 0.00  ← Referência fixa
│ TURF: PC1         ●─────── -0.45  ← 45% menor que M. hispida
│ CCA: PC1          ●──────── 0.18  ← 18% maior que M. hispida
└───────────────────────────────→
        -0.5   0   0.5
```

**M. hispida em 0 é o ponto de ancoragem** contra o qual julgamos todas as outras espécies.

### Nos PDPs

Quando calculamos as **proporções preditas** (não os log-ratios), M. hispida aparece como uma curva normal porque usamos:

```
Proporção_MUSS = α_MUSS / (α_MUSS + α_TURF + α_CCA + α_CYANO + α_MACRO)
```

Onde `α_MUSS` é derivada do fato de que `log(α_MUSS/α_MUSS) = 0`.

---

## 8. Resumo Final

### Pergunta: "Por que M. hispida não tem coeficientes?"

**Resposta em 3 pontos:**

1. **Necessidade matemática:** Modelos de composição precisam de uma referência para serem identificáveis

2. **Definição da referência:** A primeira coluna no `cbind()` é automaticamente a referência no `brms`

3. **Log-ratio zero:** Por definição, `log(α_MUSS/α_MUSS) = log(1) = 0`, então não há parâmetro a estimar

### Analogia Final

> **Pergunta:** "Por que não temos coeficiente para M. hispida?"
>
> **Resposta:** "Por a mesma razão que não medimos a distância de alguém relativa a si mesmo - é sempre zero. Medimos as distâncias dos OUTROS relativo a essa pessoa."

---

## Referências Técnicas

1. **Aitchison, J. (1986)**. *The Statistical Analysis of Compositional Data*. Chapman & Hall.
   - Fundamentos teóricos de análise de dados composicionais

2. **Bürkner, P. C. (2017)**. brms: An R Package for Bayesian Multilevel Models.
   - Documentação da parametrização Dirichlet no brms

3. **Tsagris, M., & Stewart, C. (2018)**. A Review of Flexible Transformations for Modeling Compositional Data.
   - Discussão sobre transformações log-ratio

---

## Verificação Prática

No R, você pode verificar isso:

```r
# Ver os coeficientes do modelo
summary(model)$fixed

# Notar que todos começam com:
# - muTURFprop_*
# - muCCAprop_*
# - muCYANOprop_*
# - muMACROALGAEprop_*
# Mas NÃO existe muMUSSISMILIAprop_*

# Isso confirma que MUSSISMILIA é a referência!
```
