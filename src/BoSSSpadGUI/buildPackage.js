const packager = require('electron-packager')
packager(
{
    "dir" : "./",
    "arch" : "x64",
    "platform" : "win32",
    "out" : "./InnoSetup/bin/",
    "overwrite" : "true"
}
)