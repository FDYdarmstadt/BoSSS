dotnet publish ../../src/L4-application/BoSSSpad/BoSSSpad.csproj -c Release -v q -o .

jupyter.exe nbconvert --execute ./Part0_PrintingNip_Setup.ipynb --to html --output-dir ./PrintingNip/Output

jupyter.exe nbconvert --execute ./Part1_PrintingNip_Correlation_Run.ipynb --to html --output-dir ./PrintingNip/Output
jupyter.exe nbconvert --execute ./Part1_PrintingNip_Correlation_Evaluate.ipynb --to html --output-dir ./PrintingNip/Output

jupyter.exe nbconvert --execute ./Part2_PrintingNip_ConstantStagnationPoint_Run.ipynb --to html --output-dir ./PrintingNip/Output
jupyter.exe nbconvert --execute ./Part2_PrintingNip_ConstantStagnationPoint_Evaluate.ipynb --to html --output-dir ./PrintingNip/Output

jupyter.exe nbconvert --execute ./Part3_PrintingNip_SimulateExperiment_Run.ipynb --to html --output-dir ./PrintingNip/Output
jupyter.exe nbconvert --execute ./Part3_PrintingNip_SimulateExperiment_Evaluate.ipynb --to html --output-dir ./PrintingNip/Output
jupyter.exe nbconvert --execute ./Part4_PrintingNip_SimulateUXMap_Run.ipynb --to html --output-dir ./PrintingNip/Output

jupyter.exe nbconvert --execute ./Part6_PrintingNip_TearDown.ipynb --to html --output-dir ./PrintingNip/Output