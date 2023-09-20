pwd 

echo 'some info: '
which bash
bash --version
echo '===================='

# using aliases for wget/find here; PATH for our gitlab runner is messed up,
# seems to use awkward versions of this
myWget='/c/ProgramData/chocolatey/bin/wget'
myFind='/c/Program Files/Git/usr/bin/find.exe'

# download tests run at FDY-internal HPC cluster

downloadFromJenkins() {
    $myWget $1
    unzip archive.zip
    cd archive
    ls
    FILES=$("$myFind" . -name "$2")
    echo "${FILES}"
    echo 'copy'
    cp $FILES ./..
    cd ..
    rm -r archive*
}

# next ...
downloadFromJenkins http://130.83.248.207:8080/view/BoSSS%20Master%20Pipelines/job/BoSSS-Stage4-ValidationTest-atFDYcluster/lastSuccessfulBuild/artifact/*zip*/archive.zip '*.html'
downloadFromJenkins http://130.83.248.207:8080/view/BoSSS%20Master%20Pipelines/job/BoSSS-Stage4-ValidationTest-atLichtenberg/lastSuccessfulBuild/artifact/*zip*/archive.zip '*.html'
downloadFromJenkins http://130.83.248.207:8080/view/BoSSS%20Master%20Pipelines/job/BoSSS-Stage1-Installer-Win/lastSuccessfulBuild/artifact/*zip*/archive.zip '*.exe'
downloadFromJenkins http://130.83.248.207:8080/view/BoSSS%20Master%20Pipelines/job/BoSSS-Stage1-Installer-Linux/lastSuccessfulBuild/artifact/*zip*/archive.zip '*.run'

cp `ls BoSSS-setup-*.exe` BoSSS-setup-latest.exe
cp `ls BoSSS-setup-*.run` BoSSS-setup-latest.run
