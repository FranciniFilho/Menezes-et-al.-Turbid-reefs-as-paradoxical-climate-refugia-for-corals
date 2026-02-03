
import pandas as pd
import sys

# Path to the dataset (hardcoded for the CV_02 example, which represents structure for all)
DATA_PATH = r"c:\Users\rbfra\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######FINAL\#######FINAL_RESULTS\output_DADOS_FINAIS_PARA_MODELAGEM_cv_2\dados_finais_para_modelagem_com_ARCH.csv"

def analyze_dataset():
    print("--- INICIANDO ANÁLISE DE VIABILIDADE DO DATASET ---")
    try:
        # Load data
        df = pd.read_csv(DATA_PATH, sep=';', decimal=',')
        print(f"Dataset carregado com sucesso. Total de linhas: {len(df)}")
        print(f"Colunas: {list(df.columns)}")
        
        # 1. SAMPLE SIZES GLOBAL
        print("\n--- 1. TAMANHO AMOSTRAL (N) ---")
        print(f"Total N = {len(df)}")
        
        # 2. FACTOR LEVELS
        print("\n--- 2. CONTAGEM POR FATORES ---")
        print("\nHAB:")
        print(df['HAB'].value_counts())
        
        print("\nARCH:")
        print(df['ARCH'].value_counts())
        
        print("\nSITE (Total de Sites Únicos):")
        unique_sites = df['SITE'].unique()
        print(f"N_SITES = {len(unique_sites)}")
        print(unique_sites)
        
        # 3. INTERACTION BALANCE (HAB x ARCH)
        print("\n--- 3. BALANCEAMENTO ARCH x HAB ---")
        crosstab = pd.crosstab(df['ARCH'], df['HAB'])
        print(crosstab)
        
        # Check potential issue with interaction term
        if (crosstab == 0).any().any():
            print("\n⚠️ ALERTA: Existem combinações vazias de ARCH x HAB!")
        else:
            print("\n✓ Todas as combinações de ARCH x HAB têm dados.")

        # 4. REPLICATES PER SITE
        print("\n--- 4. RÉPLICAS POR SITE ---")
        site_counts = df['SITE'].value_counts()
        print(site_counts)
        print(f"\nMínimo de réplicas: {site_counts.min()}")
        print(f"Máximo de réplicas: {site_counts.max()}")
        print(f"Média de réplicas: {site_counts.mean():.2f}")
        
        if site_counts.min() < 3:
             print("\n⚠️ ALERTA: Alguns sites têm muito poucas réplicas (<3)!")
        
        # 5. RESPONSE VARIABLES CHECK
        print("\n--- 5. VARIÁVEIS RESPOSTA (NAs) ---")
        for var in ['RGR', 'HEALTH_PC1', 'HEALTH_PC2']:
            if var in df.columns:
                na_count = df[var].isna().sum()
                print(f"{var}: {na_count} NAs")
            else:
                print(f"{var}: COLUNA NÃO ENCONTRADA")

        # 6. HABITAT MERGE SIMULATION
        print("\n--- 6. SIMULAÇÃO DE MERGE DE HABITAT (RR + TP) ---")
        df['HAB_MERGED'] = df['HAB'].replace({'RR': 'RR_TP', 'TP': 'RR_TP'})
        print(df['HAB_MERGED'].value_counts())
        
        crosstab_merged = pd.crosstab(df['ARCH'], df['HAB_MERGED'])
        print("\nCruzamento ARCH x HAB_MERGED:")
        print(crosstab_merged)

    except Exception as e:
        print(f"Erro ao analisar o dataset: {e}")

if __name__ == "__main__":
    analyze_dataset()
