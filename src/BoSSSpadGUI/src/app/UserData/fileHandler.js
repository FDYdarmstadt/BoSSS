const fs = require('fs');

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
            fs.writeFile(path, data, options, function(error){
                if(error){
                    reject(error);
                } 
                else{
                    resolve(fs);
                }
            });
        });
    }

    async load(){
        try{
            var string = await this.readFile(this.dataPath);
        }
        catch(error){
            console.log(error);
        }
        return string;
    }

    readFile(path){
        return new Promise(function(resolve, reject){
            fs.readFile(path, {encoding: 'utf8', flag: 'a+'}, function(error, data){
                if(error){
                    reject(error);
                } 
                else{
                    resolve(data);
                }
            });
        });
    }
}

module.exports = FileHandler;