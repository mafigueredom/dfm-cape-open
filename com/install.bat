@echo off
REM Register the DFM CAPE-OPEN 1.1 unit (64-bit). Run elevated on Windows.
setlocal
set DIR=%~dp0
set REGASM=%WINDIR%\Microsoft.NET\Framework64\v4.0.30319\regasm.exe
if not exist "%REGASM%" (
  echo regasm.exe not found. Install .NET Framework 4.8, then rebuild.
  exit /b 1
)
set DLL=%DIR%DfmMethanation\bin\Release\net48\DfmMethanation.dll
if not exist "%DLL%" (
  echo DLL not found: %DLL%
  echo From the cape-open-aspen folder run:  build.bat
  exit /b 1
)
"%REGASM%" /codebase /tlb "%DLL%"
if errorlevel 1 exit /b 1
echo Registered %DLL%
echo Palette name: DFM Methanation Reactor
echo ProgID: PhD.DFM.Methanation.1
echo Engine image: dfm-methanation-cape:v1
echo Restart Aspen Plus so it re-reads the CAPE-OPEN category.
endlocal
