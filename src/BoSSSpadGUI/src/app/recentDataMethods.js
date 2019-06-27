const DataMethods = require('./dataMethods.js');

class RecentDataMethods extends DataMethods{

    constructor(mainWindow, recentDocuments){
        super(mainWindow);
        this.recentDocuments = recentDocuments;
    }

    open(filePath){
        this.recentDocuments.addRecentDocument(filePath);
        return super.open(filePath);      
    }

    save(filePath){
        this.recentDocuments.addRecentDocument(filePath);
        return super.save(filePath);
    }

    openFileFromPath(filePath){
        var that = this;
        that.AreYouSure_Save(() => 
        {
            try{
                that.open(filePath);
            }
            catch(err) {
                console.log(err);
            }
        });
    }
}

module.exports = RecentDataMethods;