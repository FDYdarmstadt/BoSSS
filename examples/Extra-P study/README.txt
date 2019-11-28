preparation:
1. 	Build the project ~\internal\src\experimental\ilPSP.Cube\ilPSP.Cube.csproj in Release Mode
2. 	make folder "ilPSP.Cube_new" and copy the binaries there
3. 	Get the cubew.dll from https://www.scalasca.org/scalasca/software/cube-3.x/download.html
	and also put it into folder "ilPSP.Cube_new"
4. 	done.

execution:
1. 	type in cmd (%1: path to session): generating_CubeFiles %1
	for example: generating_CubeFiles: generating_CubeFiles.bat V:\SyncHHLR\DB_Cube_2
2.	The batch skript will execute GetSomeInfo.bws, which will generate session.txt,
	which specifies the sessions to copy the profiling_bin.* from. GetSomeInfo.bws can also be executed 
	manually. Do not forget to comment out the .bws execution in generating_CubeFiles.bat
3.	done.
	
