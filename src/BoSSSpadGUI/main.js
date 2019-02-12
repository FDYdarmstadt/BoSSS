const electron = require('electron')
const app = electron.app;
const BrowserWindow = electron.BrowserWindow;

let mainWindow;

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

app.on('ready', createWindow)

app.on('window-all-closed', function () {
	if (process.platform !== 'darwin') {
		app.quit()
	}
})

app.on('activate', function () {
	if (mainWindow === null) {
		createWindow()
	}
})

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
        changed =>
        {
            if(changed)
            {
                var response = electron.dialog.showMessageBox({
                    type: "question", 
                    buttons: ["save", "discard"],
                    message: "Do you want to save the changes you have made?",
                    defaultId: 0,
                    cancelId: 1
                });
                //Save changes
                if(response == 0)
                {
                    menu.saveFile().then(func());
                }
                //Discard changes
                else if(response == 1){
                    func();
                }
                //cancel: keep open
                else if(response == 2){
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
							that.saveFile();
						}
					},
					{
						label: 'Save File As...',
						click() {
							that.saveFileAs();
						}
					}
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
			},
			{	
				label: 'Dev',
				submenu:[
					{role: 'reload'},
					{role: 'toggledevtools'}
				]
			}
		];
		const menu = this.Menu.buildFromTemplate(template);
		this.Menu.setApplicationMenu(menu);
	}

	openFile(){
		//bookmarks is mac stuff
		function bosssPadOpenFile(filePaths, bookmarks ){
			try{
				var filePath = '"' + filePaths[0].replace(/\\/g, '/')+ '"';
				var command = 'BoSSSpad.openFile(' + filePath + ');';
				this.mainWindow.webContents.executeJavaScript( command );
				this.savePath = filePath;
				this.electron.app.addRecentDocument(filePath);
			}
			catch(err) {
				console.log(err);
			}
		}
		
		this.dialog.showOpenDialog({
			properties: ['openFile'], 
			filters: [
				{name: 'Supported Files', extensions: ['bws', 'tex']},
				{name: 'BoSSSPad Worksheet', extensions: ['bws']},
				{name: 'TeX Files', extensions: ['tex']},
				{name: 'All Files', extensions: ['*']}
			]
		}, bosssPadOpenFile.bind(this));
	}

	saveFile(){
        var that = this;
        var save = new Promise(function (resolve, reject)
        {
            if (that.savePath === null){
                that.saveFileAs().then(resolve());
            }else{
                var command = 'BoSSSpad.saveFile(' + that.savePath + ');';
                that.mainWindow.webContents.executeJavaScript( command ).then(resolve());
            }
        });
        return save;
	}

	saveFileAs(){
        //bookmarks is mac stuff
        var that = this;
        var saveAs = new Promise(function(resolve, reject)
        {
            function bosssPadSaveFile(fileName, bookmarks ){		
                try{
                    var filePath = '"' + fileName.replace(/\\/g, '/')+ '"';
                    that.savePath = filePath;
                    var command = 'BoSSSpad.saveFile(' + filePath + ');';
                    that.mainWindow.webContents.executeJavaScript( command ).then(resolve());
                }
                catch(err){
                    console.log(err);
                    reject();
                }
            }
        
            var path = that.dialog.showSaveDialog({
                filters: [
                    {name: 'BoSSSPad Worksheet', extensions: ['bws']},
                    {name: 'TeX Files', extensions: ['tex']},
                    {name: 'All Files', extensions: ['*']}
                ]
            });
            bosssPadSaveFile(path);
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