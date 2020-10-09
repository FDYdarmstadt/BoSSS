#!/usr/bin/bash
#
cd $BOSSS_INSTALL/public/src/L3-solution/BoSSS.Solution/

sed -i -e "s:AssemblyInformationalVersion(\".*\"):AssemblyInformationalVersion(\"\"):g" ./Properties/AssemblyInfo.cs
