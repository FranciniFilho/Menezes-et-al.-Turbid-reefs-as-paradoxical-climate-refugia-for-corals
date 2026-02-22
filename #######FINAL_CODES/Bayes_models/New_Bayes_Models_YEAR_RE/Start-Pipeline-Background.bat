@echo off
chcp 65001 >nul
echo ╔════════════════════════════════════════════════════════════╗
echo ║  INICIANDO PIPELINE BAYESIANO EM BACKGROUND                ║
echo ╚════════════════════════════════════════════════════════════╝
echo.
echo Diretório: %CD%
echo Log: C:\Users\rbfra\OneDrive\New_Bayes_Models_Output_PRIOR_SENSITIVITY_FULLGRID_v1\pipeline_background.log
echo.

:: Criar diretório de log se não existir
if not exist "C:\Users\rbfra\OneDrive\New_Bayes_Models_Output_PRIOR_SENSITIVITY_FULLGRID_v1" mkdir "C:\Users\rbfra\OneDrive\New_Bayes_Models_Output_PRIOR_SENSITIVITY_FULLGRID_v1"

:: Iniciar em background usando start /B
start /B "Pipeline Bayesian" cmd /c "Rscript 00_RUN_ALL_YEAR_RE_MODELS.R prior_sens_fullgrid WeaklyInformative 2^>^&1 ^| tee C:\Users\rbfra\OneDrive\New_Bayes_Models_Output_PRIOR_SENSITIVITY_FULLGRID_v1\pipeline_background.log"

echo ✓ Pipeline iniciado em background!
echo.
echo Para monitorar o progresso, use um dos comandos abaixo:
echo.
echo 1. Ver log em tempo real (novo PowerShell):
echo    Get-Content 'C:\Users\rbfra\OneDrive\New_Bayes_Models_Output_PRIOR_SENSITIVITY_FULLGRID_v1\pipeline_background.log' -Wait -Tail 30
echo.
echo 2. Ver últimas 50 linhas:
echo    Get-Content 'C:\Users\rbfra\OneDrive\New_Bayes_Models_Output_PRIOR_SENSITIVITY_FULLGRID_v1\pipeline_background.log' -Tail 50
echo.
echo 3. Contar arquivos gerados:
echo    (Get-ChildItem -Path 'C:\Users\rbfra\OneDrive\New_Bayes_Models_Output_PRIOR_SENSITIVITY_FULLGRID_v1' -Recurse -Filter '*.rds').Count
echo.
pause
