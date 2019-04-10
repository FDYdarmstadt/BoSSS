const electron = require('electron')
const app = electron.app;
const BrowserWindow = electron.BrowserWindow;

let mainWindow;
let menu;

function createWindow () {
	mainWindow = new BrowserWindow({width: 800, height: 600, icon: __dirname + '/src/img/logo/logo_fdy.ico'})
	mainWindow.loadURL(`file://${__dirname}/dist/index.html`)
	mainWindow.on('closed', function(){ mainWindow = null}); 
    menu = new BoSSSMenu(mainWindow, electron);

    var closeBoSSSpad = false;
    mainWindow.on('close', (event) =>{
        if(!closeBoSSSpad)
        {
            event.preventDefault();
            AreYouSure_Save(() => {closeBoSSSpad = true; mainWindow.close(); });
        }
    });   
}

app.on('ready', createWindow);

app.on('window-all-closed', function () {
	if (process.platform !== 'darwin') {
		app.quit()
	}
});

app.on('activate', function () {
	if (mainWindow === null) {
		createWindow()
	}
});

/* Mac Stuff
app.on('activate', () => {
    // Unter macOS ist es üblich ein neues Fenster der App zu erstellen, wenn
    // das Dock Icon angeklickt wird und keine anderen Fenster offen sind.
    if (win === null) {
      createWindow()
    }
})
*/

function AreYouSure_Save(func){
    mainWindow.webContents.executeJavaScript( 'BoSSSpad.hasChanged()').then(
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
                    await menu.saveFile();
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

class BoSSSMenu{
	constructor(mainWindow, electron){
		this.mainWindow = mainWindow;
		this.Menu = electron.Menu;
		this.electron = electron;
		this.dialog = electron.dialog;
		this.savePath = null;
		this.createMenu();
	}
			
	createMenu(){
		var that = this;
		const template = [
			{
				label: 'File',
				submenu: [
					{
                        label: 'New File',
                        accelerator: 'CmdOrCtrl+N',
						click() {
							that.newFile();
						}
					},
					{type: 'separator'},
					{
                        label: 'Open File',
                        accelerator: 'CmdOrCtrl+O',
						click() {
							that.openFile();
						}
					},
					{
						label: 'Open Recent', 
						role: 'recentdocuments',
						submenu: [{
							label: 'Clear Recent',
							role: 'clearrecentdocuments'
						}]
					},
					{type: 'separator'},
					{
                        label: 'Save File',
                        accelerator: 'CmdOrCtrl+S',
						click() {
                            try{
                                that.saveFile();
                            }catch(e){
                                console.log(e);
                            }
						}
					},
					{
						label: 'Save File As...',
						click() {
                            try{
                                that.saveFileAs();
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
                            that.executeFromHere();
                        }
                    },
                    {
                        label: 'Execute entire worksheet',
                        accelerator: 'F5',
                        click(){
                            that.executeEntireWorksheet();
                        }
                    },
                    {
                        label: 'Execute until here',
                        accelerator: 'CmdOrCtrl+Shift+F5',
                        click(){
                            that.executeUntilHere();
                        }
                    },
                    {
                        label: 'Interrupt current command',
                        click(){
                            that.interruptCurrentCommand();
                        }
                    },
                    {
                        label: 'De-Queue all pending commands',
                        click(){
                            that.deQueueAllPendingCommands();
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
		const menu = this.Menu.buildFromTemplate(template);
		this.Menu.setApplicationMenu(menu);
    }
    
    //Command actions
    //===========================================================================
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

    //File actions
    //===========================================================================
	openFile(){
        //bookmarks is mac stuff
        var that = this;
        AreYouSure_Save(() => {
            function bosssPadOpenFile(filePaths, bookmarks ){
                try{
                    var filePath = '"' + filePaths[0].replace(/\\/g, '/')+ '"';
                    var command = 'BoSSSpad.openFile(' + filePath + ');';
                    that.mainWindow.webContents.executeJavaScript( command );
                    that.savePath = filePath;
                    that.electron.app.addRecentDocument(filePath);
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

	async saveFile(){
        var that = this;

        if (that.savePath === null){
            await that.saveFileAs();
        }else{
            var command = 'BoSSSpad.saveFile(' + that.savePath + ');';
            await that.mainWindow.webContents.executeJavaScript( command )
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
                    var filePath = '"' + fileName.replace(/\\/g, '/')+ '"';
                    that.savePath = filePath;
                    var command = 'BoSSSpad.saveFile(' + filePath + ');';
                    that.mainWindow.webContents.executeJavaScript( command ).then( a => {console.log("Saved!"); resolve()});
                }
                catch(err){
                    console.log(err);
                    reject();
                }
            }
        });
        return saveAs;
	}

	newFile(){
        AreYouSure_Save(() =>{
            this.savePath = null;
            this.mainWindow.webContents.executeJavaScript(
                'BoSSSpad.resetFile();'
            );
        });
	}
}