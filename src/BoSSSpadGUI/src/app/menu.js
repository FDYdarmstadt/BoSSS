const electron = require('electron');
const BoSSSDataMethods = require('./dataMethods.js');

class BoSSSMenu{
	constructor(commandActions, dataMethods){
        this.Menu = electron.Menu;
        this.dataMethods = dataMethods;
        this.commandActions = commandActions;
        this.createMenu();
	}
			
	createMenu(){
		const template = this.createTemplate();
		const menu = this.Menu.buildFromTemplate(template);
		this.Menu.setApplicationMenu(menu);
    }

    createTemplate(){
        var that = this;
        var template = [
			{
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
						submenu: [
                            {
                                label: "Irgend/Ein/Path.bws"
                            },
                            {
                                type: 'separator'
                            },
                            {
                                label: 'Clear Recent',
                            }
						]
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
            },
            {
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
            },
            {	
				label: 'Dev',
				submenu:[
					{role: 'reload'},
					{role: 'toggledevtools'}
				]
			},
			{
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
			}
		];
        return template;
    }
}

module.exports = BoSSSMenu;
