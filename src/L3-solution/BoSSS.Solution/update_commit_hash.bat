@echo off
setlocal enabledelayedexpansion



:: Check if git is available
git --version 2>nul || (
    echo. > commit_hash.txt
    exit /b 0
)

::
:: Note: all the following is only to make sure that the `commit_hash.txt` is only updated if the commit hash is different;
:: otherwise, the incremental build will not work
::



if exist commit_hash.txt (
    set /p CURRENT_HASH=<commit_hash.txt
) else (
    set CURRENT_HASH=
)

echo "current hash:" %CURRENT_HASH%



FOR /F "tokens=* USEBACKQ" %%F IN (`git rev-parse HEAD`) DO (
    SET NEW_HASH=%%F
)



::echo %NEW_HASH%
::echo "equality " "!CURRENT_HASH!" equ "!NEW_HASH!"

if not %CURRENT_HASH% equ %NEW_HASH% (
    echo "writing new"
    echo !NEW_HASH! > commit_hash.txt
)
