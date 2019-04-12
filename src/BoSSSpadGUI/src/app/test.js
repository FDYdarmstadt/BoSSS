const UserDatabase = require('./UserData/userDatabase.js');
function test(userDatabase){
    var userData = userDatabase.getUserData();
    console.log(userData);
    userDatabase.save();
}

UserDatabase.load().then(test);

