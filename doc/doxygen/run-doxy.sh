rm -r ./latex
rm -r ./html
doxygen -b - < Doxyfile > doxyout.txt 2> doxyerr.txt
#doxygen Doxyfile 
#RES=0
#grep 'error: Problems running latex' doxyout.txt
#TESTRES=$?
