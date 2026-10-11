param([string]$Python = 'python', [switch]$SkipLong, [switch]$SkipDocs)
$ErrorActionPreference = 'Stop'
$Study = Join-Path $PSScriptRoot 'studies\nhanes3_testretest'
if (Test-Path (Join-Path $Study 'analysis_v1')) { throw 'analysis_v1 already exists: run in a fresh copy of this folder' }
$env:PYTHONUTF8 = '1'; $env:PYTHONDONTWRITEBYTECODE = '1'
Set-Location (Join-Path $PSScriptRoot 'scripts')
function Run($s) { Write-Host "== $s"; & $Python -u $s; if ($LASTEXITCODE -ne 0) { throw "$s failed" } }
Run nh3_testretest_stage12.py
Copy-Item (Join-Path $Study 'records\*') (Join-Path $Study 'analysis_v1') -Force
Run nh3_stage3.py
Run nh3_p8_independent.py
Run nh3_stage4.py
Run nh3_posthoc_v1.py
Run nh3_posthoc_v2.py
Run nh3_s9_correction.py
Run nh3_table1.py
Run nh3_reclassification_lookup.py
Run nh3_posthoc_repeat_strategies.py
Run nh3_posthoc_cvd75_curve.py
Run nh3_ext2001_labs.py
Run nh3_ext2021_burden.py
Run nh3_ext2_analyses.py
Run nh3_ext3_statin.py
Run nh3_ext4_30yr.py
if (-not $SkipLong) {
    Run nh3_ext5_stakes.py
    Run nh3_ext2_e4_uniform.py
    Run nh3_ext6_review.py
}
Run nh3_letter_numbers.py
Run nh3_figures.py
Run nh3_fig_e10.py
Run nh3_fig_e11.py
Run nh3_figures_v7.py
Run export_figures_ajpc.py
Run nh3_fig_ci_ajpc.py
Run build_evidence_ledger.py
if (-not $SkipDocs) {
    Run build_nh3_manuscript_v7.py
    Run build_nh3_supplement_v7.py
    Run build_nh3_manuscript_ajpc.py
}
Write-Host 'done'
