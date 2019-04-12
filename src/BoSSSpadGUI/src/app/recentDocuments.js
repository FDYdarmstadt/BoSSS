class RecentDocuments{
    constructor(userData){
        this.path = userData.paths;
    }

    getRecentPaths(){
        return paths;
    }

    addRecentPath(path){
        this.paths.push(path);
    }
}

module.exports = RecentDocuments;