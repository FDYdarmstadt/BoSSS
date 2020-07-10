#!/usr/bin/bash
#
cd $BOSSS_INSTALL/public/src/L3-solution/BoSSS.Solution/

/usr/bin/sed -i -e "s:AssemblyInformationalVersion(\".*\"):AssemblyInformationalVersion(\"$(/usr/bin/git rev-parse HEAD)\"):g" ./Properties/AssemblyInfo.cs
