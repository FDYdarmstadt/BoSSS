const FileHandler = require('./fileHandler.js');
//const convert = require('xml-js');

class File{

    static async initialize(path){
        var fileHandler, text;
        fileHandler = new FileHandler(path);
        text = await fileHandler.load();
        return new File(fileHandler, text);
    }

    constructor(fileHandler, text) { 
        this.fileHandler = fileHandler;
        this.object = this.createFile(text);
    }

    createFile(text){
        var object;
        if(text == ""){
            object = null;
        }
        else{
            object = File.objectFromText(text);
        }
        return object;
    }

    static objectFromText(text){
        var file;
        file = JSON.parse(text);
        return file;
    }

    async save(){
        var text;
        text = JSON.stringify(this.object);
        await this.fileHandler.save(text);
    }

    getObject(){
        return this.object;
    }

    setObject(object){
        this.object = object;
    }
}

module.exports = File;