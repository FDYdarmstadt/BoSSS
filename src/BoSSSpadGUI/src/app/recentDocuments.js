const RecentPaths = require('./UserData/recentPaths.js');

class RecentDocuments{

    constructor(userData){
        var numberOfRecentDocuments = 5;
        this.recentPaths = new RecentPaths(userData, numberOfRecentDocuments);
        this.onAddRecentDocument = null;
    }

    clearRecentDocuments(){
        this.recentPaths.clearRecentPaths();
    }

    getRecentDocuments( onClickFunction){
        var paths = this.recentPaths.getRecentPaths(); 
        var recentDocuments;
        if(paths.length > 0){
            recentDocuments = this.createRecentDocuments(paths, onClickFunction);
        }else{
            recentDocuments = [{
                label: "..."
            }];
        }
        this.attachClearRecentMenu(recentDocuments);
        return recentDocuments;
    }

    createRecentDocuments(paths, onClickFunction){
        var recentDocument = [];
        for(var i = 0; i < paths.length; ++i){
            var path = paths[i].repeat(1);
            recentDocument.push(
            {
                label: this.createLabel(path),
                click(){
                    onClickFunction(path);
                }
            });
        }
        return recentDocument;
    }

    createLabel(path){
        var start = 1;
        start = Math.max(start, path.length - 35);
        var label = path.slice(start,-1);
        if(start == 1){
            return label;
        }
        else{
            return "..." + label;
        }
    }

    attachClearRecentMenu(menu){
        var that = this;
        menu.push({
            type: 'separator'
        });
        menu.push(
        {
            label: 'Clear Recent',
            click() {
                that.clearRecentDocuments();
                if(that.onAddRecentDocument != null){
                    that.onAddRecentDocument();
                }
            }
        });
    }

    addRecentDocument(path){
        this.recentPaths.addRecentPath(path);
        if(this.onAddRecentDocument != null){
            this.onAddRecentDocument();
        }
    };
}

module.exports = RecentDocuments;