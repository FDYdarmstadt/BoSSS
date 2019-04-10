const fs = require('fs');

class UserData{

    constructor(dataPath){
        this.dataPath = dataPath;
    }

    save(){

    }

    load(){
        var load = new Promise(function(){
            fs.open(dataPath,'w+',);
        });
        
    }
}

const convert = require('xml-js');

class XMLData extends UserData{
    constructor(xmlPath) { 
        super(xmlPath);
    }
}

class RecentPaths{

    constructor(){
        var xmlPath = `file://${__dirname}/src/app/userData.xml`;
        this.userData = new XMLData();
    }

    add(path){

    }

    openXML(){

    }

    getAll(){
        
    }
}