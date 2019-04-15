const BoSSSMenu = require('./boSSSMenu.js');
const BoSSSDataMethods = require('./recentDataMethods.js');
const CommandActions = require('./commandActions.js');
const UserDatabase = require('./UserData/userDatabase.js');
const RecentDocuments = require('./recentDocuments.js')

async function addFunctionality(mainWindow){
    
    function createMenu(){
        var commandActions = new CommandActions(mainWindow);
        var menu = new BoSSSMenu(commandActions, dataMethods, recentDocuments);
        return menu;
    }

    function addCloseAction(){
        mainWindow.on('close', (event) =>{
            if(!closeBoSSSpad)
            {
                event.preventDefault();
                dataMethods.AreYouSure_Save(() => {closeBoSSSpad = true; mainWindow.close(); });
            }
        });   
    }

    var userDatabase = await UserDatabase.load('./userData.xml');
    //var recentDocuments = new RecentDocuments(userDatabase.getUserData());
    //var dataMethods = new BoSSSDataMethods(mainWindow, recentDocuments);
    //createMenu();
    //addCloseAction();
    //var closeBoSSSpad = false;
}

module.exports = addFunctionality;