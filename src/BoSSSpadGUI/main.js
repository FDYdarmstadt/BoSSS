const electron = require('electron');
const app = electron.app;
const BrowserWindow = electron.BrowserWindow;
const addFunctionality = require('./src/app/functionality.js');

let mainWindow;

function createWindow () {
	mainWindow = new BrowserWindow({width: 800, height: 600, icon: __dirname + '/src/app/logo/logo_fdy.ico'});
	mainWindow.loadURL(`file://${__dirname}/boSSSpad/index.html`);
    mainWindow.on('closed', function(){ mainWindow = null}); 
    addFunctionality(mainWindow);
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

