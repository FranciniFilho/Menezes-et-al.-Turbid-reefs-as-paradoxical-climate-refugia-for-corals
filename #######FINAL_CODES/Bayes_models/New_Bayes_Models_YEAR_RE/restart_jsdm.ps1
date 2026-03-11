# JSDM Restart Script - Runs with fixed prior code
$ErrorActionPreference = "Continue"

cd "C:\Users\rbfra\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######FINAL\#######FINAL_CODES\Bayes_models\New_Bayes_Models_YEAR_RE"

$LogFile = "C:\Users\rbfra\OneDrive\New_Bayes_Models_Output_PRIOR_SENSITIVITY_FULLGRID_v1\jsdm_fixed_v2.log"

Write-Host "Starting JSDM with fixed code..."
Write-Host "Log: $LogFile"

& Rscript "03_JSDM_YEAR_RE_Full_LOO_Selection.R" 2>&1 | Out-File -FilePath $LogFile -Encoding UTF8

Write-Host "JSDM script completed."
