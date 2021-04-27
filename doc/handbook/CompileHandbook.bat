cd ../../src/L4-application/TutorialTests/bin/Release
.\TutorialTests.exe
cd ../../../../../doc/handbook/
pdflatex -halt-on-error BoSSShandbook.tex
biber BoSSShandbook.tex
pdflatex -halt-on-error BoSSShandbook.tex
pdflatex -halt-on-error BoSSShandbook.tex
pdflatex -halt-on-error BoSSShandbook.tex
