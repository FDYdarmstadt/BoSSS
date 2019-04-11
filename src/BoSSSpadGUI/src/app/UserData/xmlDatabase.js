const Database = require('./database.js')
const convert = require('xml-js');

class XMLFile{
    constructor(database, xmlData) { 
        this.xmlData = xmlData;
        this.database = database;
    }

    static async build(xmlPath){
        var database = new Database(xmlPath);
        var xmlText = await database.load();
        var xmlData = convert.xml2js(xmlText);
        return new XMLDatabase(database, xmlData);
    }

    getXML(){
        return xmlData;
    }

    setXML(xml){
        this.xmlData = xml;
    }

    async saveXML(){
        var xmlText = convert.js2xml(this.xmlData);
        await this.database.save(xmlText);
    }

    async test(){
        try{
            console.log(this.xmlData);
        }
        catch(e){
            console.log("Exception: " + e);
        }
    }
}

module.exports = XMLFile;