@echo off
REM Build the DFM CAPE engine image and the 64-bit COM DLL.
REM Run from this folder. Docker Desktop must be running.
setlocal
cd /d "%~dp0"

docker build -f docker\Dockerfile -t dfm-methanation-cape:v1 .
if errorlevel 1 exit /b 1

dotnet build -c Release com\DfmMethanation\DfmMethanation.csproj
if errorlevel 1 exit /b 1

echo.
echo Image: dfm-methanation-cape:v1
echo DLL:   com\DfmMethanation\bin\Release\net48\DfmMethanation.dll
echo Next, from an elevated command prompt:  com\install.bat
endlocal
