#!/bin/bash
export WORKINGDIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" >/dev/null 2>&1 && pwd )"
echo $WORKINGDIR
export BOSSS_INSTALL=$WORKINGDIR
printf "BOSSS_INSTALL set to $BOSSS_INSTALL\n"
printf "To permanantly set this add \e[2mexport BOSSS_INSTALL=$BOSSS_INSTALL\e[0m to your .bashrc\n"
