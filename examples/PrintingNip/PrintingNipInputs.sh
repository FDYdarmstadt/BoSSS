#!/bin/bash

# change to script path
path=$(dirname "$0")
cd $path
cd ./PrintingNip

# clear file
echo -n > ./Output/PrintingNipInputs.tex

echo "\newpage" >> ./Output/PrintingNipInputs.tex
echo "\section{Part1}" >> ./Output/PrintingNipInputs.tex
echo "Parameter study with varying nip width \$\varepsilon\$, velocity boundary condition \$V_W\$ and pressure boundary condition \$\Delta p\$. The predicate H0, V0, P0 means that the nip wdith, velocity or pressure are kept constant for that plot. In case of the boundary conditions these are zero. The second predicate H, V, P denotes the variable, which is constant in each plot line. \newline" >> ./Output/PrintingNipInputs.tex
echo "\newpage" >> ./Output/PrintingNipInputs.tex

# scan dir
for filepath in ./Figures/Part1/*; do
    filename=$(basename -- "$filepath")
    stem="${filename%.*}"
    extension="${filename##*.}"
    
    # append import statement
    if [ $extension == "tex" ]; then
        echo "\begin{figure}[H]" >> ./Output/PrintingNipInputs.tex
        echo "\centering" >> ./Output/PrintingNipInputs.tex
        echo "\caption{\detokenize{${stem}}}" >> ./Output/PrintingNipInputs.tex
        echo "\resizebox{!}{0.4\textheight}{\input{./PrintingNip/Figures/Part1/${stem}.tex}}" >> ./Output/PrintingNipInputs.tex
        echo "\end{figure}" >> ./Output/PrintingNipInputs.tex
    fi
done

echo "\newpage" >> ./Output/PrintingNipInputs.tex
echo "\section{Part2}" >> ./Output/PrintingNipInputs.tex
echo "Here the pressure boundary condition is tuned, such that the stagnation point is kept at constant \$10mm\$. The plot lines show results for a constant velocity boundary condition.\newline" >> ./Output/PrintingNipInputs.tex
echo "\newpage" >> ./Output/PrintingNipInputs.tex

for filepath in ./Figures/Part2/*; do
    filename=$(basename -- "$filepath")
    stem="${filename%.*}"
    extension="${filename##*.}"
    
    # append import statement
    if [ $extension == "tex" ]; then
        echo "\begin{figure}[H]" >> ./Output/PrintingNipInputs.tex
        echo "\centering" >> ./Output/PrintingNipInputs.tex
        echo "\caption{\detokenize{${stem}}}" >> ./Output/PrintingNipInputs.tex
        echo "\resizebox{!}{0.4\textheight}{\input{./PrintingNip/Figures/Part2/${stem}.tex}}" >> ./Output/PrintingNipInputs.tex
        echo "\end{figure}" >> ./Output/PrintingNipInputs.tex
    fi
done

echo "\newpage" >> ./Output/PrintingNipInputs.tex
echo "\section{Part3}" >> ./Output/PrintingNipInputs.tex
echo "Here the pressure boundary condition is tuned, such that the stagnation point is kept at the experimentally calculated point for a certain velocity. In the experiment the geometry of the printing cylinders can be altered. In the simulation this is only considered through the stagnation point position. In plots with predicate H0 the nip width is not changed \$\varepsilon = 10^{-5}\$. In plots with R0 the Raster size is at constant \$80 lines / cm\$. \newline" >> ./Output/PrintingNipInputs.tex
echo "\newpage" >> ./Output/PrintingNipInputs.tex


for filepath in ./Figures/Part3/*; do
    filename=$(basename -- "$filepath")
    stem="${filename%.*}"
    extension="${filename##*.}"
    
    # append import statement
    if [ $extension == "tex" ]; then
        echo "\begin{figure}[H]" >> ./Output/PrintingNipInputs.tex
        echo "\centering" >> ./Output/PrintingNipInputs.tex
        echo "\caption{\detokenize{${stem}}}" >> ./Output/PrintingNipInputs.tex
        echo "\resizebox{!}{0.4\textheight}{\input{./PrintingNip/Figures/Part3/${stem}.tex}}" >> ./Output/PrintingNipInputs.tex
        echo "\end{figure}" >> ./Output/PrintingNipInputs.tex
    fi
done
echo "\newpage" >> ./Output/PrintingNipInputs.tex
echo "\section{Testplots}" >> ./Output/PrintingNipInputs.tex
echo "\newpage" >> ./Output/PrintingNipInputs.tex

for filepath in ./Figures/Test/*; do
    filename=$(basename -- "$filepath")
    stem="${filename%.*}"
    extension="${filename##*.}"
    
    # append import statement
    if [ $extension == "tex" ]; then
        echo "\begin{figure}[H]" >> ./Output/PrintingNipInputs.tex
        echo "\centering" >> ./Output/PrintingNipInputs.tex
        echo "\caption{\detokenize{${stem}}}" >> ./Output/PrintingNipInputs.tex
        echo "\resizebox{!}{0.4\textheight}{\input{./PrintingNip/Figures/Test/${stem}.tex}}" >> ./Output/PrintingNipInputs.tex
        echo "\end{figure}" >> ./Output/PrintingNipInputs.tex
    fi
done