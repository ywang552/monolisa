@echo off
setlocal enabledelayedexpansion

set USERNAME=ywang552
set PASSWORD=NOwjyq25!!
set ADDRESS=its-zest-login3
set LOCAL_PATH=C:\Users\ywang552\Downloads\strap-main\monolisa\src
set REMOTE_PATH=/home/ywang552/monolisa
set PLINK_PATH=C:\path\to\plink.exe

echo Starting file transfer...

REM Iterate through all files in LOCAL_PATH
for /r "%LOCAL_PATH%" %%F in (*) do (
    REM Debug: Print local file path
    echo Local File Path: %%F

    REM Extract the relative path of the file
    set FILE=%%F
    set RELATIVE_PATH=!FILE:%LOCAL_PATH%=!

    REM Remove leading backslash
    set RELATIVE_PATH=!RELATIVE_PATH:~1!

    REM Construct the remote file path
    set REMOTE_FILE=%REMOTE_PATH%/!RELATIVE_PATH!

    REM Debug: Print remote file path
    echo Remote File Path: !REMOTE_FILE!

    REM Transfer the file
    %PLINK_PATH% -batch -l %USERNAME% -pw %PASSWORD% "cat > !REMOTE_FILE!" < "%%F"
    
    if %ERRORLEVEL% neq 0 (
        echo Failed to transfer: %%F
    )
)

echo File transfer process completed.
pause
