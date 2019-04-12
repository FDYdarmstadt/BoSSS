const fs = require('fs').promises;

class FileHandler{

    constructor(dataPath){
        this.dataPath = dataPath;
    }

    async save(string){
        try{
            await fs.writeFile(this.dataPath, string,'utf8');
        }catch(error){
            console.log(error);
        }
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