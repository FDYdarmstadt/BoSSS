const BoSSSMenu = require('./boSSSMenu.js');
const BoSSSDataMethods = require('./dataMethods.js');
const CommandActions = require('./commandActions.js');
const UserDatabase = require('./UserData/userDatabase.js');
const RecentDocuments = require('./recentDocuments.js')

function addFunctionality(mainWindow){
    
    function createMenu(dataMethods, userDatabase){
        var commandActions = new CommandActions(mainWindow);
        var recentDocuments = new RecentDocuments(userDatabase);
        var menu = new BoSSSMenu(commandActions, dataMethods, recentDocuments);
        return menu;
    }

    function addCloseAction(dataMethods){
        mainWindow.on('close', (event) =>{
            if(!closeBoSSSpad)
            {
                event.preventDefault();
                dataMethods.AreYouSure_Save(() => {closeBoSSSpad = true; mainWindow.close(); });
            }
        });   
    }

    var dataMethods = new BoSSSDataMethods(mainWindow);
    var userDatabase = UserDatabase.load();
    createMenu(dataMethods, userDatabase);
    addCloseAction(dataMethods);
    var closeBoSSSpad = false;
}

module.exports = addFunctionality;