@echo off
setlocal
set DIR=%~dp0
set REGASM=%WINDIR%\Microsoft.NET\Framework64\v4.0.30319\regasm.exe
set DLL=%DIR%DfmMethanation\bin\Release\net48\DfmMethanation.dll
if not exist "%DLL%" (
  echo DLL not found: %DLL%
  exit /b 1
)
"%REGASM%" /unregister "%DLL%"
endlocal
