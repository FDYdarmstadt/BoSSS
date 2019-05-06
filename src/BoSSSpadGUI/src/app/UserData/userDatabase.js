const UserData = require("./userData.js");
const File = require('./file.js');
const Path = require('path');

class UserDatabase{
    
    static async load(path){
        var filePath = UserDatabase.getDatabasePath(path);
        var file = await File.initialize(filePath);
        return new UserDatabase(file);
    }
    
    constructor(file){
        this.file = file;
        this.UserData = this.createUserData(file);
    }

    static getDatabasePath(path){
        var userDataFolder = process.env.APPDATA || 
            (process.platform == 'darwin' ? process.env.HOME + 'Library/Preferences' : process.env.HOME + "/.local/share");
        var filePath = Path.join(userDataFolder, path);
        return filePath;
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