const FileHandler = require('./fileHandler.js');
const convert = require('xml-js');

class File{

    static async initialize(xmlPath){
        var fileHandler, xmlText;
        fileHandler = new FileHandler(xmlPath);
        xmlText = await fileHandler.load();
        return new File(fileHandler, xmlText);
    }

    constructor(fileHandler, xmlText) { 
        this.fileHandler = fileHandler;
        this.object;
        
        this.object = this.createFile(xmlText);
    }

    createFile(xmlText){
        var object;
        if(xmlText == ""){
            object = null;
        }
        else{
            object = File.objectFromText(xmlText);
        }
        return object;
    }

    static objectFromText(xmlText){
        var file, settings;
        settings = File.getConverterSettings();
        file = convert.xml2js(xmlText, settings);
        return file;
    }

    static getConverterSettings(){
        return {compact: true, spaces: 4};
    }

    getObject(){
        return this.object;
    }

    setObject(object){
        this.object = object;
    }

    async save(){
        var xmlText, settings;
        settings = File.getConverterSettings();
        xmlText = convert.js2xml(this.object, settings);
        await this.fileHandler.save(xmlText);
    }
}

module.exports = File;