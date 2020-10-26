cd ../../src/L4-application/TutorialTests/bin/Release
.\TutorialTests.exe
cd ../../../../../doc/handbook/
pdflatex BoSSShandbook.tex
biber BoSSShandbook.tex
pdflatex BoSSShandbook.tex
pdflatex BoSSShandbook.tex