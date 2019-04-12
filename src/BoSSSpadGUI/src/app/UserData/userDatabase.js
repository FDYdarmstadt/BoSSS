const UserData = require("./userData.js");
const File = require('./file.js');

class UserDatabase{
    
    static async load(){
        var file = await File.initialize('./userData.xml');
        return new UserDatabase(file);
    }
    
    constructor(file){
        this.file = file;
        this.UserData = this.createUserData(file);
    }

    createUserData(file){
        var userData = file.getObject();
        if(userData == null){
            userData = new UserData();
        }
        return userData;
    }

    getUserData(){
        return this.UserData;
    }

    async save(){
        this.file.setObject(this.UserData);
        await this.file.save();
    };
}

module.exports = UserDatabase;