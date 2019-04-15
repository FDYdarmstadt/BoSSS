const UserDatabase = require('./UserData/userDatabase.js');
const RecentDocuments = require('./recentDocuments.js')

function test(userDatabase){
    var recentDocuments = new RecentDocuments(userDatabase.getUserData());
    recentDocuments.addRecentPath("SomeTestPath");
    userDatabase.save();
}

UserDatabase.load('./userData.txt').then(test);

