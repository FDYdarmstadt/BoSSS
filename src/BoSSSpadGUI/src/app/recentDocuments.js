const RecentPaths = require('./UserData/recentPaths.js');

class RecentDocuments{

    constructor(userData){
        this.recentPaths = new RecentPaths(userData, 5);
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
            var path = paths[i];
            recentDocument.push(
            {
                label: path.slice(1,-1),
                click(){
                    onClickFunction(path);
                }
            });
        }
        return recentDocument;
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