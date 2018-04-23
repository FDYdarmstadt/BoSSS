const electron = require('electron')
const app = electron.app
const BrowserWindow = electron.BrowserWindow;

let mainWindow;

function createWindow () {
	mainWindow = new BrowserWindow({width: 800, height: 600})
	mainWindow.loadURL(`file://${__dirname}/dist/index.html`)
	mainWindow.on('closed', function () {
		mainWindow = null
	})
	menu = new BoSSSMenu(mainWindow, electron);
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
		if (this.savePath === null){
			this.saveFileAs();
		}else{
			var command = 'BoSSSpad.saveFile(' + this.savePath + ');';
			this.mainWindow.webContents.executeJavaScript( command );
		}
	}

	saveFileAs(){
		//bookmarks is mac stuff
		function bosssPadSaveFile(fileName, bookmarks ){		
			try{
				var filePath = '"' + fileName.replace(/\\/g, '/')+ '"';
				this.savePath = filePath;
				var command = 'BoSSSpad.saveFile(' + filePath + ');';
				this.mainWindow.webContents.executeJavaScript( command );
		
			}
			catch(err){
				console.log(err);
			}
		}
	
		this.dialog.showSaveDialog({
			filters: [
				{name: 'BoSSSPad Worksheet', extensions: ['bws']},
				{name: 'TeX Files', extensions: ['tex']},
				{name: 'All Files', extensions: ['*']}
			]
		}, bosssPadSaveFile.bind(this));
	}

	newFile(){
		this.savePath = null;
		this.mainWindow.webContents.executeJavaScript(
			'BoSSSpad.resetFile();'
		);
	}
}