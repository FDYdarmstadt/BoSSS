class RecentPaths{

    constructor(userData, maxNumber){
        this.paths = userData.paths;
        this.maxPathNumber = maxNumber;
    }

    clearRecentPaths(){
        this.paths.length = 0;
    }

    getRecentPaths(){
        return this.paths;
    }

    addRecentPath(path){
        if(this.notInPaths(path)){
            if(this.paths.length >= this.maxPathNumber){
                this.paths.pop();
            }
            this.paths.unshift(path);
        }
        else{
            this.movePathToTop(path);
        }
    }

    notInPaths(path){
        var inPaths = this.paths.includes(path);
        return !inPaths;
    }

    movePathToTop(path){
        console.log("MovingPathToTop");
        var index = this.paths.indexOf(path);
        this.paths.splice(index, 1);
        this.paths.unshift(path);
    }
}

module.exports = RecentPaths;