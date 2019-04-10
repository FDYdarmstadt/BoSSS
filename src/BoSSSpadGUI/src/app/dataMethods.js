const electron = require('electron');

class BoSSSDataMethods{
    
    constructor(mainWindow){
        this.mainWindow = mainWindow;
        this.savePath = null;
        this.dialog = electron.dialog;
    }

    openFile(){
        //bookmarks is mac stuff
        var that = this;
        this.AreYouSure_Save(() => 
        {
            function bosssPadOpenFile(filePaths, bookmarks ){
                try{
                    that.setPathAndOpen(filePaths[0]);
                }
                catch(err) {
                    console.log(err);
                }
            }
            
            that.dialog.showOpenDialog({
                properties: ['openFile'], 
                filters: [
                    {name: 'Supported Files', extensions: ['bws', 'tex']},
                    {name: 'BoSSSPad Worksheet', extensions: ['bws']},
                    {name: 'TeX Files', extensions: ['tex']},
                    {name: 'All Files', extensions: ['*']}
                ]
            }, bosssPadOpenFile);
        });
    }

    setPathAndOpen(filePath){
        var that = this;
        var filePath = '"' + filePath.replace(/\\/g, '/')+ '"';
        var command = 'BoSSSpad.openFile(' + filePath + ');';
        that.mainWindow.webContents.executeJavaScript( command );
        that.savePath = filePath;
        electron.app.addRecentDocument(filePath);
    }

	async saveFile(){
        var that = this;

        if (that.savePath === null){
            await that.saveFileAs();
        }else{
            await that.save( that.savePath );
            console.log("Saved!");
        }
	}

	saveFileAs(){
        //bookmarks is mac stuff
        var that = this;
        var saveAs = new Promise(function(resolve, reject)
        {
            that.dialog.showSaveDialog(
                {
                    filters: [
                        {name: 'BoSSSPad Worksheet', extensions: ['bws']},
                        {name: 'TeX Files', extensions: ['tex']},
                        {name: 'All Files', extensions: ['*']}
                    ]
                },
                bosssPadSaveFile
            );
            
            function bosssPadSaveFile(fileName, bookmarks ){		
                try{
                    that.setPathAndSave(fileName).then( a => {console.log("Saved!"); resolve()});
                }
                catch(err){
                    console.log(err);
                    reject();
                }
            }
        });
        return saveAs;
    }

    setPathAndSave(fileName){
        var that = this;
        var filePath = '"' + fileName.replace(/\\/g, '/')+ '"';
        that.savePath = filePath;
        return that.save(filePath);
    }
    
    save(filePath){
        var that = this;
        electron.app.addRecentDocument(filePath);
        var command = 'BoSSSpad.saveFile(' + filePath + ');';
        return that.mainWindow.webContents.executeJavaScript( command );
    }

	newFile(){
        this.AreYouSure_Save(() =>{
            this.savePath = null;
            this.mainWindow.webContents.executeJavaScript(
                'BoSSSpad.resetFile();'
            );
        });
    }

    AreYouSure_Save(func){
        this.mainWindow.webContents.executeJavaScript( 'BoSSSpad.hasChanged()').then(
            async(changed) =>
            {
                if(changed)
                {
                    var response = electron.dialog.showMessageBox({
                        type: "question", 
                        buttons: ["save", "discard"],
                        message: "Do you want to save the changes you have made?",
                        defaultId: 0,
                        cancelId: 2
                    });
                    //Save changes
                    if(response == 0)
                    {
                        await this.saveFile();
                        func();
                    }
                    //Discard changes
                    else if(response == 1){
                        func();
                    }
                    //cancel: keep open
                    else if(response == 2){
                        //Do nothing
                    }
                }else{
                    func();
                }
            }
        )
    }
}

module.exports = BoSSSDataMethods;