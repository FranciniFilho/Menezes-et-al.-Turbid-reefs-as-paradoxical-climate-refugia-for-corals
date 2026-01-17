# Especificação Técnica: Sistema de Atualização Robusta da Base de Dados MODIS/CRW

> [!IMPORTANT]
> Este documento serve como o "Modelo Mestre" para a implementação do script `Smart_MODIS_database_update_v2.py`. Ele detalha a lógica exata, tratamento de erros e arquitetura que deve ser seguida.

## 1. Visão Geral do Sistema

O objetivo é criar um script de manutenção de dados que seja **à prova de falhas**. Ele não deve apenas baixar arquivos, mas garantir que *o que está no disco é válido*.

### Arquitetura de 3 Fases
O script deve executar sequencialmente:
1.  **Auditoria (Local)**: Varre o disco, valida integridade real (abre o arquivo) e lista arquivos corrompidos.
2.  **Mapeamento (Rede)**: Verifica no servidor (NASA/NOAA) o que *deveria* existir versus o que temos.
3.  **Sincronização (Download)**: Baixa o necessário usando **Escritas Atômicas** e **Concorrência**.

---

## 2. Detalhamento Técnico das Fases

### Fase 1: Auditoria de Integridade Profunda
Não confie apenas na extensão `.nc`. O script deve verificar se o arquivo é um NetCDF válido.

**Requisitos da Função `verify_integrity(filepath)`:**
1.  **Verificação de Tamanho**: Se `size < 1KB`, marcar como `CORRUPTED_EMPTY`.
2.  **Verificação de Estrutura**: Tentar abrir com `xarray.open_dataset(..., engine='netcdf4')`.
    *   Se falhar com erro de I/O ou HDF5: marcar como `CORRUPTED_HDF_ERROR`.
3.  **Verificação de Conteúdo**:
    *   Verificar se a lista de variáveis (`ds.data_vars`) não está vazia.
    *   Verificar se existe a chave de tempo/data esperada.
4.  **Ação**: Arquivos corrompidos devem ser movidos para uma pasta `_QUARANTINE` ou listados para deleção, nunca deletados silenciosamente sem log.

### Fase 2: Mapeamento de Lacunas (Gap Analysis)
O script deve cruzar o **Inventário Local Válido** (da Fase 1) com o **Inventário Remoto**.

**Lógica de Detecção:**
*   **Para MODIS**:
    *   Verificar URL padrão (Final).
    *   Se 404, verificar URL NRT (Near Real-Time).
    *   Se ambos 404, marcar como `MISSING_ON_SERVER` (provável gap do satélite).
*   **Para CRW (Coral Reef Watch)**:
    *   Os arquivos são organizados por ano. O script deve iterar pelos anos relevantes.
    *   Deve respeitar a estrutura anual de diretórios.

**Otimização de Performance:**
*   Usar `concurrent.futures.ThreadPoolExecutor` para fazer requisições `HEAD` em paralelo.
*   Isso acelera a verificação de milhares de dias de horas para minutos.

### Fase 3: Download Inteligente e Atômico
Esta é a parte mais crítica para evitar arquivos corrompidos no futuro.

**Padrão de Escrita Atômica (Atomic Write Pattern):**
1.  Identificar URL alvo e arquivo de destino final (ex: `AQUA_MODIS...nc`).
2.  Baixar para arquivo temporário: `AQUA_MODIS...nc.tmp`.
3.  **Após download completo**:
    *   Verificar tamanho do arquivo `.tmp`.
    *   (Opcional) Verificar se abre com `xarray`.
4.  Somente se válido: `os.replace('...nc.tmp', '...nc')`.
    *   Isso garante que se o script for interrompido no meio, não sobra um arquivo `.nc` corrompido pela metade.

**Controle de Concorrência:**
*   Usar `ThreadPoolExecutor` com `max_workers=4` (evitar bloqueio pelo servidor da NASA).
*   Implementar `Backoff Exponencial` em caso de erro 5xx ou Timeout.

---

## 3. Estrutura de Configuration (Schema)

O script deve possuir um dicionário de configuração centralizado e fácil de editar.

```python
DATA_CONFIG = {
    "MODIS_CHL": {
        "enabled": True,
        "local_dir": r"H:\remote sensing\MODIS_DATA_FULL",
        "url_template": "https://oceandata.sci.gsfc.nasa.gov/cgi/getfile/AQUA_MODIS.{date}.L3m.DAY.CHL.chlor_a.4km.nc",
        "filename_pattern": r"AQUA_MODIS\.(\d{8})\.L3m\.DAY\.CHL\.chlor_a\.4km(\.NRT)?\.nc",
        "var_name_check": "chlor_a",
        "start_date": "2002-07-04",
        "end_date": "CURRENT", # Palavra chave para 'hoje'
        "auth_required": True,
        "nrt_fallback": True
    },
    # ... outros produtos (SST, PAR, KD490, CRW_DHW, CRW_SST) seguindo o mesmo schema
}
```

---

## 4. Prompt para Geração do Código (Meta-Prompt)

*Este bloco abaixo é destinado a ser copiado e colado para um LLM gerar o código final.*

```text
Atue como um Engenheiro de Dados Sênior Especialista em Python.
Escreva um script Python completo chamado "Smart_MODIS_Database_Update.py" que implemente as seguintes especificações rigorosamente:

1.  **Bibliotecas**: Use `requests`, `xarray`, `netCDF4`, `pandas`, `concurrent.futures`, `pathlib`, `logging`.
2.  **Classe `IntegrityChecker`**:
    - Método `check_file(path)`: Retorna bool. Abre o NetCDF, verifica se não está vazio e se tem variáveis.
    - Se arquivo corrupto: Logar erro e adicionar à lista de "To Download".
3.  **Classe `GapManager`**:
    - Gera lista de datas esperadas (pd.date_range).
    - Compara com arquivos locais válidos.
    - Usa `ThreadPoolExecutor` para fazer HEAD requests nas URLs da NASA/NOAA para verificar disponibilidade remota (lidar com 404).
    - Suporta lógica de fallback: Se URL principal falhar, tentar URL ".NRT.nc" (apenas para MODIS).
4.  **Classe `AtomicDownloader`**:
    - Baixa para `filename.nc.tmp`.
    - Só renomeia para `filename.nc` após download 100% concluído E validação básica de integridade.
    - Usa barra de progresso (`tqdm`) global ou por arquivo.
    - Retry logic: 3 tentativas com wait time exponencial.
5.  **Logging**:
    - Criar `audit.log` (erros de integridade encontrada).
    - Criar `download.log` (operações de download).
    - Printar resumo estatístico ao final (Total, Válidos, Corrompidos, Baixados, Falhas).

**Contexto dos Dados**:
- MODIS Aqua (CHL, SST, KD490, PAR) via oceanColor (requer Token NASA via Header).
- CRW (DHW, SST) via NOAA Star (diretórios anuais).

**Atenção**: O código deve ser robusto a falhas de rede e interrupções. Use Type Hinting.
```

## 5. Próximos Passos para Implementação

1.  Gerar o script usando o prompt acima.
2.  Configurar o `AUTH_TOKEN` correto.
3.  Executar em modo `--dry-run` (implementar flag) para ver o log de auditoria antes de baixar.
