import os
import shutil

BASE = r"C:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######FINAL"
FINAL_DATA = os.path.join(BASE, "#######FINAL_DATA")
FINAL_RESULTS = os.path.join(BASE, "#######FINAL_RESULTS")
OLD_DATA = r"C:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######22.04.23\DATA"
CEBIMAR = r"C:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\########CEBIMAR\####PROJETOS\#####Coral trade offs"

# Step 1: Create subdirectory structure
subdirs = [
    "00_sites_metadata",
    "01_raw_biological",
    "02_PCA_environmental/CV_02",
    "02_PCA_environmental/CV_30",
    "02_PCA_environmental/CV_ALL",
    "03_modeling_data/CV_02",
    "03_modeling_data/CV_30",
    "03_modeling_data/CV_ALL",
    "04_health_growth_PCA",
    "05_benthic_cover",
]

for d in subdirs:
    path = os.path.join(FINAL_DATA, d)
    os.makedirs(path, exist_ok=True)
    print(f"[DIR] Created: {d}")

# Step 2: Copy site metadata files
copies = []

# From CEBIMAR
for f in ["sites_list_full.csv", "sites_list_full_clean.csv"]:
    src = os.path.join(CEBIMAR, f)
    dst = os.path.join(FINAL_DATA, "00_sites_metadata", f)
    if os.path.exists(src) and not os.path.exists(dst):
        shutil.copy2(src, dst)
        copies.append(f"[COPY] {f} -> 00_sites_metadata/")
    elif os.path.exists(dst):
        copies.append(f"[SKIP] {f} already exists in 00_sites_metadata/")
    else:
        copies.append(f"[MISS] {src} not found!")

# From old project DATA
for f in ["sites_list_full_clean_PLEST.csv"]:
    src = os.path.join(OLD_DATA, f)
    dst = os.path.join(FINAL_DATA, "00_sites_metadata", f)
    if os.path.exists(src) and not os.path.exists(dst):
        shutil.copy2(src, dst)
        copies.append(f"[COPY] {f} -> 00_sites_metadata/")
    elif os.path.exists(dst):
        copies.append(f"[SKIP] {f} already exists in 00_sites_metadata/")
    else:
        copies.append(f"[MISS] {src} not found!")

# Step 3: Copy raw biological data
src = os.path.join(OLD_DATA, "BENTHOS TEMPORAL CLEAN_v2.csv")
dst = os.path.join(FINAL_DATA, "01_raw_biological", "BENTHOS_TEMPORAL_CLEAN_v2.csv")
if os.path.exists(src) and not os.path.exists(dst):
    shutil.copy2(src, dst)
    copies.append(f"[COPY] BENTHOS TEMPORAL CLEAN_v2.csv -> 01_raw_biological/")
elif os.path.exists(dst):
    copies.append(f"[SKIP] BENTHOS_TEMPORAL_CLEAN_v2.csv already exists")
else:
    copies.append(f"[MISS] {src} not found!")

# Also copy existing FINAL_DATA files to proper subdirs
for f in ["Vitality and size_new.csv"]:
    src = os.path.join(FINAL_DATA, f)
    dst = os.path.join(FINAL_DATA, "01_raw_biological", f)
    if os.path.exists(src) and not os.path.exists(dst):
        shutil.copy2(src, dst)
        copies.append(f"[COPY] {f} -> 01_raw_biological/")

for f in ["######Vitality and size_new.xlsx"]:
    src = os.path.join(FINAL_DATA, f)
    dst = os.path.join(FINAL_DATA, "01_raw_biological", "Vitality_and_size_new.xlsx")
    if os.path.exists(src) and not os.path.exists(dst):
        shutil.copy2(src, dst)
        copies.append(f"[COPY] {f} -> 01_raw_biological/Vitality_and_size_new.xlsx")

# Copy sites_list_full.csv from root to 00_sites_metadata
src = os.path.join(FINAL_DATA, "sites_list_full.csv")
dst = os.path.join(FINAL_DATA, "00_sites_metadata", "sites_list_full.csv")
if os.path.exists(src) and not os.path.exists(dst):
    shutil.copy2(src, dst)
    copies.append(f"[COPY] sites_list_full.csv (root) -> 00_sites_metadata/")

# Step 4: Copy PCA environmental data
cv_map = {
    "CV_02": "#####output_local_PCA_CV_2_FINAL",
    "CV_30": "#####output_local_PCA_CV_30_FINAL",
    "CV_ALL": "#####output_local_PCA_CV_all_FINAL",
}

for cv_label, cv_dir in cv_map.items():
    src_dir = os.path.join(FINAL_RESULTS, cv_dir)
    dst_dir = os.path.join(FINAL_DATA, "02_PCA_environmental", cv_label)
    if os.path.exists(src_dir):
        for f in os.listdir(src_dir):
            if f.endswith(('.csv', '.xlsx')):
                src = os.path.join(src_dir, f)
                dst = os.path.join(dst_dir, f)
                if not os.path.exists(dst):
                    shutil.copy2(src, dst)
                    copies.append(f"[COPY] {f} -> 02_PCA_environmental/{cv_label}/")
                else:
                    copies.append(f"[SKIP] {f} already in {cv_label}/")
    else:
        copies.append(f"[MISS] Directory {cv_dir} not found!")

# Step 5: Copy modeling data
modeling_map = {
    "CV_02": "output_DADOS_FINAIS_PARA_MODELAGEM_cv_2",
    "CV_30": "output_DADOS_FINAIS_PARA_MODELAGEM_cv_30",
    "CV_ALL": "output_DADOS_FINAIS_PARA_MODELAGEM_cv_all",
}

for cv_label, cv_dir in modeling_map.items():
    src_dir = os.path.join(FINAL_RESULTS, cv_dir)
    dst_dir = os.path.join(FINAL_DATA, "03_modeling_data", cv_label)
    if os.path.exists(src_dir):
        for f in os.listdir(src_dir):
            if f.endswith(('.csv', '.xlsx')):
                src = os.path.join(src_dir, f)
                dst = os.path.join(dst_dir, f)
                if not os.path.exists(dst):
                    shutil.copy2(src, dst)
                    copies.append(f"[COPY] {f} -> 03_modeling_data/{cv_label}/")
                else:
                    copies.append(f"[SKIP] {f} already in {cv_label}/")
    else:
        copies.append(f"[MISS] Directory {cv_dir} not found!")

# Step 6: Copy health/growth PCA data
health_dir = os.path.join(FINAL_RESULTS, "#output_ANALISE_BIOLOGICA_boxplot_PCAnew_YEAR_RE")
dst_dir = os.path.join(FINAL_DATA, "04_health_growth_PCA")
if os.path.exists(health_dir):
    for f in os.listdir(health_dir):
        if f.endswith(('.csv', '.xlsx')):
            src = os.path.join(health_dir, f)
            dst = os.path.join(dst_dir, f)
            if not os.path.exists(dst):
                shutil.copy2(src, dst)
                copies.append(f"[COPY] {f} -> 04_health_growth_PCA/")
            else:
                copies.append(f"[SKIP] {f} already in 04_health_growth_PCA/")
else:
    copies.append(f"[MISS] Health PCA directory not found!")

# Step 7: Copy benthic cover data
benthic_dir = os.path.join(FINAL_RESULTS, "%%BOX_PLOTS_benthic_cover")
dst_dir = os.path.join(FINAL_DATA, "05_benthic_cover")
if os.path.exists(benthic_dir):
    for f in os.listdir(benthic_dir):
        if f.endswith(('.csv', '.xlsx')):
            src = os.path.join(benthic_dir, f)
            dst = os.path.join(dst_dir, f)
            if not os.path.exists(dst):
                shutil.copy2(src, dst)
                copies.append(f"[COPY] {f} -> 05_benthic_cover/")
            else:
                copies.append(f"[SKIP] {f} already in 05_benthic_cover/")
else:
    copies.append(f"[MISS] Benthic cover directory not found!")

# Print summary
print("\n=== COPY OPERATIONS SUMMARY ===")
for c in copies:
    print(c)

# Verify final structure
print("\n=== FINAL STRUCTURE ===")
for root, dirs, files in os.walk(FINAL_DATA):
    level = root.replace(FINAL_DATA, '').count(os.sep)
    indent = ' ' * 2 * level
    basename = os.path.basename(root) if level > 0 else "#######FINAL_DATA"
    print(f"{indent}{basename}/")
    subindent = ' ' * 2 * (level + 1)
    for file in sorted(files):
        size = os.path.getsize(os.path.join(root, file))
        print(f"{subindent}{file} ({size:,} bytes)")
