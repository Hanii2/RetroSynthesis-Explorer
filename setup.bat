@echo off
cd /d "%~dp0"
if exist .venv (
  echo Purana .venv mila. Agar pehle error aaya tha to ise delete karo aur dobara chalao.
  pause
  exit /b 1
)
set PYV=
py -3.12 --version >nul 2>&1 && set PYV=3.12
if "%PYV%"=="" py -3.11 --version >nul 2>&1 && set PYV=3.11
if "%PYV%"=="" py -3.10 --version >nul 2>&1 && set PYV=3.10
if "%PYV%"=="" (
  echo Python 3.10-3.12 nahi mila, jo Python pehle se hai wohi use kar raha hoon.
  python -m venv .venv
  goto :install
)
py -%PYV% -m venv .venv
:install
call .venv\Scripts\activate
python -m pip install --upgrade pip
pip install -r requirements.txt
echo.
echo Setup khatam. Ab run_app.bat par double-click karo.
pause
