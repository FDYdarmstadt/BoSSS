const electron = require('electron');

class BoSSSMenu{
	constructor(commandActions, dataMethods, recentDocuments){
        this.dataMethods = dataMethods;
        this.commandActions = commandActions;
        this.recentDocuments = recentDocuments;
        this.recentDocuments.onAddRecentDocument = this.createMenu.bind(this);
        this.createMenu();
	}
			
	createMenu(){
		this.template = this.createTemplate();
        const menu = electron.Menu.buildFromTemplate(this.template);
		electron.Menu.setApplicationMenu(menu);
    }

    createTemplate(){
        var that = this;
        var file = this.buildFileMenuTemplate();
        var commands = {
            label: 'Commands',
            submenu:[
                {
                    label: 'Execute from here...',
                    accelerator: 'CmdOrCtrl+F5',
                    click(){
                        that.commandActions.executeFromHere();
                    }
                },
                {
                    label: 'Execute entire worksheet',
                    accelerator: 'F5',
                    click(){
                        that.commandActions.executeEntireWorksheet();
                    }
                },
                {
                    label: 'Execute until here',
                    accelerator: 'CmdOrCtrl+Shift+F5',
                    click(){
                        that.commandActions.executeUntilHere();
                    }
                },
                {
                    label: 'Interrupt current command',
                    click(){
                        that.commandActions.interruptCurrentCommand();
                    }
                },
                {
                    label: 'De-Queue all pending commands',
                    click(){
                        that.commandActions.deQueueAllPendingCommands();
                    }
                },
            ]
        };
        var dev = {	
            label: 'Dev',
            submenu:[
                {role: 'reload'},
                {role: 'toggledevtools'}
            ]
        };
        var help = {
            label: 'Help',
            submenu: [
                {
                    label: 'Documentation',
                    click () { 
                        electron.shell.openExternal('https://github.com/FDYdarmstadt/BoSSS');
                    } 
                },
                {
                    label: 'About',
                    click () { 
                        electron.shell.openExternal('http://www.fdy.tu-darmstadt.de/fdy/fdyresearch/bossscode/framework/framework.de.jsp');
                    } 
                }
            ]
        };
        var template = [
			file,
            commands,
            dev,
			help
		];
        return template;
    }

    buildFileMenuTemplate(){
        var that = this;
        var recentDocuments = this.recentDocuments.getRecentDocuments(
            this.dataMethods.openFileFromPath.bind(this.dataMethods)
            );
        var file = {
            label: 'File',
            submenu: [
                {
                    label: 'New File',
                    accelerator: 'CmdOrCtrl+N',
                    click() {
                        that.dataMethods.newFile();
                    }
                },
                {type: 'separator'},
                {
                    label: 'Open File',
                    accelerator: 'CmdOrCtrl+O',
                    click() {
                        that.dataMethods.openFile();
                    }
                },
                {
                    label: 'Open Recent', 
                    submenu: recentDocuments
                        
                },
                {type: 'separator'},
                {
                    label: 'Save File',
                    accelerator: 'CmdOrCtrl+S',
                    click() {
                        try{
                            that.dataMethods.saveFile();
                        }catch(e){
                            console.log(e);
                        }
                    }
                },
                {
                    label: 'Save File As...',
                    click() {
                        try{
                            that.dataMethods.saveFileAs();
                        }catch(e){
                            console.log(e);
                        }
                    }
                }
            ]
        };
        return file;
    }
}

module.exports = BoSSSMenu;
