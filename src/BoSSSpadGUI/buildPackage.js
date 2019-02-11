const packager = require('electron-packager');
var path = require('path');

//Win32_x64 Package
packager(
{ 
    "dir" : "./",
    "arch" : "x64",
    "platform" : "win32",
    "out" : "./InnoSetup/bin/",
    "overwrite" : "true",
    "prune" : "true",
    "ignore" : [new RegExp('InnoSetup'), new RegExp('src/cs/bin/Debug$')]
})

