# Script PowerShell para executar pipeline Bayesian em background
# Autor: Claude
# Data: 2025-02-21

param(
    [string]$Action = "start",
    [string]$LogFile = "C:\Users\rbfra\OneDrive\New_Bayes_Models_Output_PRIOR_SENSITIVITY_FULLGRID_v1\pipeline_background.log"
)

$WorkingDir = "C:\Users\rbfra\OneDrive\########PUBLICACOES\############Menezes et al. Mus his distribution and abundance Abrolhos\######FINAL\#######FINAL_CODES\Bayes_models\New_Bayes_Models_YEAR_RE"
$RScript = "C:\Program Files\R\R-4.3.1\bin\Rscript.exe"

function Start-Pipeline {
    Write-Host "╔════════════════════════════════════════════════════════════╗" -ForegroundColor Green
    Write-Host "║  INICIANDO PIPELINE BAYESIANO EM BACKGROUND                ║" -ForegroundColor Green
    Write-Host "╚════════════════════════════════════════════════════════════╝" -ForegroundColor Green
    Write-Host ""
    Write-Host "Working Directory: $WorkingDir"
    Write-Host "Log File: $LogFile"
    Write-Host ""
    
    # Criar diretório do log se não existir
    $LogDir = Split-Path $LogFile -Parent
    if (!(Test-Path $LogDir)) {
        New-Item -ItemType Directory -Path $LogDir -Force | Out-Null
    }
    
    # Iniciar job em background
    $Job = Start-Job -ScriptBlock {
        param($workingDir, $logFile, $rscript)
        Set-Location $workingDir
        & $rscript "00_RUN_ALL_YEAR_RE_MODELS.R" "prior_sens_fullgrid" "WeaklyInformative" 2>&1 | Tee-Object -FilePath $logFile
    } -ArgumentList $WorkingDir, $LogFile, $RScript
    
    Write-Host "✓ Pipeline iniciado em background!" -ForegroundColor Green
    Write-Host "Job ID: $($Job.Id)"
    Write-Host ""
    Write-Host "Comandos para monitoramento:" -ForegroundColor Yellow
    Write-Host "  Get-Job                    - Lista todos os jobs"
    Write-Host "  Get-Job -Id $($Job.Id) | Format-List   - Detalhes do job"
    Write-Host "  Receive-Job -Id $($Job.Id) -Keep      - Ver output (mantém no buffer)"
    Write-Host "  Stop-Job -Id $($Job.Id)               - Parar o job"
    Write-Host ""
    Write-Host "Para ver log em tempo real:" -ForegroundColor Cyan
    Write-Host "  Get-Content '$LogFile' -Wait -Tail 20"
    Write-Host ""
    
    return $Job.Id
}

function Show-Status {
    Write-Host "╔════════════════════════════════════════════════════════════╗" -ForegroundColor Blue
    Write-Host "║  STATUS DO PIPELINE                                        ║" -ForegroundColor Blue
    Write-Host "╚════════════════════════════════════════════════════════════╝" -ForegroundColor Blue
    Write-Host ""
    
    $Jobs = Get-Job | Where-Object { $_.Name -like "Job*" }
    
    if ($Jobs.Count -eq 0) {
        Write-Host "Nenhum job em execução." -ForegroundColor Yellow
    } else {
        foreach ($Job in $Jobs) {
            $StatusColor = switch ($Job.State) {
                "Running" { "Green" }
                "Completed" { "Green" }
                "Failed" { "Red" }
                "Stopped" { "Yellow" }
                default { "White" }
            }
            
            Write-Host "Job ID: $($Job.Id)" -NoNewline
            Write-Host " | Estado: " -NoNewline
            Write-Host "$($Job.State)" -ForegroundColor $StatusColor -NoNewline
            Write-Host " | Comando: $($Job.Command.Substring(0, [Math]::Min(50, $Job.Command.Length)))..."
        }
    }
    Write-Host ""
    
    # Verificar arquivos gerados
    $OutputDir = "C:\Users\rbfra\OneDrive\New_Bayes_Models_Output_PRIOR_SENSITIVITY_FULLGRID_v1"
    if (Test-Path $OutputDir) {
        $RdsCount = (Get-ChildItem -Path $OutputDir -Recurse -Filter "*.rds" -ErrorAction SilentlyContinue).Count
        $CsvCount = (Get-ChildItem -Path $OutputDir -Recurse -Filter "*.csv" -ErrorAction SilentlyContinue).Count
        
        Write-Host "Progresso dos resultados:" -ForegroundColor Blue
        Write-Host "  Arquivos .rds gerados: $RdsCount / 210 (esperados)"
        Write-Host "  Arquivos .csv gerados: $CsvCount"
        Write-Host ""
    }
}

function Show-Log {
    param([int]$Lines = 50)
    
    if (Test-Path $LogFile) {
        Write-Host "╔════════════════════════════════════════════════════════════╗" -ForegroundColor Magenta
        Write-Host "║  ÚLTIMAS $Lines LINHAS DO LOG                              ║" -ForegroundColor Magenta
        Write-Host "╚════════════════════════════════════════════════════════════╝" -ForegroundColor Magenta
        Write-Host ""
        Get-Content $LogFile -Tail $Lines
        Write-Host ""
    } else {
        Write-Host "Arquivo de log não encontrado: $LogFile" -ForegroundColor Red
    }
}

# Execução principal
switch ($Action.ToLower()) {
    "start" { Start-Pipeline }
    "status" { Show-Status }
    "log" { Show-Log -Lines 50 }
    "taillog" { 
        Write-Host "Pressione Ctrl+C para parar de monitorar..." -ForegroundColor Yellow
        Get-Content $LogFile -Wait -Tail 20 
    }
    default { 
        Write-Host "Uso: .\Run-Pipeline-Background.ps1 -Action [start|status|log|taillog]" -ForegroundColor Yellow
        Write-Host ""
        Write-Host "Exemplos:" -ForegroundColor Cyan
        Write-Host "  .\Run-Pipeline-Background.ps1 -Action start     # Inicia o pipeline"
        Write-Host "  .\Run-Pipeline-Background.ps1 -Action status    # Ver status"
        Write-Host "  .\Run-Pipeline-Background.ps1 -Action log       # Ver últimas 50 linhas"
        Write-Host "  .\Run-Pipeline-Background.ps1 -Action taillog   # Monitorar em tempo real"
    }
}
