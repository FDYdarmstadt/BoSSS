mkdir public
copy "src\L4-application\PublicTestRunner\bin\Release\net6.0\*.html" ".\public\"
cd public    
"C:\Program Files\Git\bin\bash.exe" ../downloadValidationTests.sh