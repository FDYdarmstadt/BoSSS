@echo off
setlocal EnableDelayedExpansion

"%BOSSS_INSTALL%\bin\Release\BoSSSpad.exe" --batch GetSomeInfo.bws
set "TARGET=E:\bosss_db_performance\sessions"
set "CUBEPATH=ilPSP.Cube_new"

::for /f %%I in (projectinfo.txt) do mkdir ilPSP.Cube_new\%%I

for /d %%D in (%TARGET%\*) do (
	set "sessionname=%%~nD"
	set "sessionpath=%%~fD"
	echo these are the bins of %%~nD
	for /f "tokens=1,2 delims=:" %%G in (sessioninfo.txt) do (
		if %%~nD == %%G (
		set counter=0
		for /f %%F in ('dir /b %%~fD ^| find "profiling_bin" ') do (
			copy "%%~fD\%%F" "."
			set /A counter=!counter!+1
		)
		echo !counter!
		(
		start cmd /c %CUBEPATH%\ilPSP.Cube.exe
		)| pause
		rename "calc.p!counter!.r1" "calc.%%H.r1"
		)
	)
)
