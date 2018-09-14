const packager = require('electron-packager');
var path = require('path');
var copydir = require('copy-dir');

//Win32_x64 Package
packager(
{ 
    "dir" : "./",
    "arch" : "x64",
    "platform" : "win32",
    "out" : "./InnoSetup/bin/",
    "overwrite" : "true",
    "prune" : "true",
    "ignore" : [new RegExp('InnoSetup'), new RegExp('src/cs$')]
}).then(copyForWindows)

function copyForWindows()
{
    copydir.sync('./src/cs/bin/Release/', 
        './InnoSetup/bin/BoSSSpad-win32-x64/', 
        function(stat, filepath, filename)
        {
            if(stat === 'file' && path.extname(filepath) === '.dll') 
            {
                return true;
            }
            return false;
        }
    );
}
