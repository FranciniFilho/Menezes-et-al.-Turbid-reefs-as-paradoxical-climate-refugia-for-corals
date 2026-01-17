import os
import sys

# Set working directory
script_dir = r'C:\Users\rbfra\ONEDRIVE_NOVO\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######FINAL\#######FINAL_CODES'
os.chdir(script_dir)
sys.path.insert(0, script_dir)

# Read the original script
with open('TIME_SERIES_VARIABILITY_ANALYSIS_GLM.py', 'r', encoding='utf-8') as f:
    script_content = f.read()

# Execute the script
exec(compile(script_content, 'TIME_SERIES_VARIABILITY_ANALYSIS_GLM.py', 'exec'))
