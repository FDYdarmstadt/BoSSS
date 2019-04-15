class RecentDocuments{

    constructor(userData, loadFunctionFactory){
        this.paths = userData.paths;
        this.maxPathNumber = 5;
        this.loadfunction;
    }

    getRecentPaths(){
        return paths;
    }

    addRecentPath(path){
        if(path.length >= this.maxPathNumber){
            this.paths.pop();
        }
        this.paths.unshift(path);

    }

    getMenuItem(path){
        var item = 
        {
            label: path,
            click(){
                loadFunctionFactory(path);
            }
        }
        return item;
    }
}

module.exports = RecentDocuments;