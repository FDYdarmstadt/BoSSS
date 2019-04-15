const RecentPaths = require('./UserData/recentPaths.js');

class RecentDocuments{

    constructor(userData, loadFunctionFactory){
        this.recentPaths = new RecentPaths(userData, 5);
        this.loadfunction;
    }

    getRecentDocuments(){
        var myRecentPaths = recentPaths.getRecentPaths(); 
        var item = 
        {
            label: path,
            click(){
                loadFunctionFactory(path);
            }
        }
        return item;
    }

    addRecentDocument(path){
        this.recentPaths.addRecentPath(path);
    };
}

module.exports = RecentDocuments;