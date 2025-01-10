dotnet publish ../../src/L4-application/BoSSSpad/BoSSSpad.csproj -c Release -v q -o . > build.out   
dotnet publish ../../src/L4-application/StokesHelical_Ak/StokesHelical_Ak.csproj -c Release -v q -o . > build2.out 
echo $?