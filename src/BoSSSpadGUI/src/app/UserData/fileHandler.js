const fs = require('fs').promises;
const fs1 = require('fs');

class FileHandler{

    constructor(dataPath){
        this.dataPath = dataPath;
    }

    async save(string){
        try{
            await this.writeFile(this.dataPath, string,'utf8');
        }catch(error){
            console.log(error);
        }
    }

    writeFile(path, data, options){
        return new Promise(function(resolve, reject){
            fs1.writeFile(path, data, options, function(error){
                if(error){
                    reject();
                } 
                else{
                    resolve(fs1);
                }
            });
        });
    }

    open(path, flags){
        return new Promise(function(resolve, reject){
            fs1.open(path, flags, function(error, fd){
                if(error){
                    reject();
                } 
                else{
                    resolve(fs1);
                }
            });
        });
    }

    async load(){
        try{
            var fileHandle = await fs.open(this.dataPath,  'a+');
            var string = fileHandle.readFile('utf8');
        }
        catch(error){
            console.log(error);
        }
        return string;
    }
}

module.exports = FileHandler;