const XMLFile = require('./xmlDatabase.js');

class UserData{

    constructor(userData){
        this.userData = userData;
    }

    static async build(){
        var xmlPath = `${__dirname}\\userData.xml`;
        var userData = await XMLFile.build(xmlPath);
        return new UserData(userData);
    }

    test(){
        this.userData.test();
    }

    addRecentDocument(path){

    }

    getDocuments(){

    }

    async save(){

    }
}

module.exports = UserData;