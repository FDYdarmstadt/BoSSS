class CommandActions{

    constructor(mainWindow){
        this.mainWindow = mainWindow;
    }

    executeFromHere(){
        var command = 'BoSSSpad.executeFromHere();';
        this.mainWindow.webContents.executeJavaScript( command );
    }

    executeEntireWorksheet(){
        var command = 'BoSSSpad.executeEntireWorksheet();';
        this.mainWindow.webContents.executeJavaScript( command );
    }

    executeUntilHere(){
        var command = 'BoSSSpad.executeUntilHere();';
        this.mainWindow.webContents.executeJavaScript( command );
    }

    interruptCurrentCommand(){
        var command = 'BoSSSpad.interruptCurrentCommand();';
        this.mainWindow.webContents.executeJavaScript( command );
    }

    deQueueAllPendingCommands(){
        var command = 'BoSSSpad.deQueueAllPendingCommands();';
        this.mainWindow.webContents.executeJavaScript( command );
    }

}

module.exports = CommandActions;