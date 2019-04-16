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
}

module.exports = RecentDataMethods;