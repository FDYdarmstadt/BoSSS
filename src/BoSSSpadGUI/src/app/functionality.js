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
        var closeBoSSSpad = false;
        mainWindow.on('close', async (event) =>{
            if(!closeBoSSSpad)
            {
                await userDatabase.save();
                event.preventDefault();
                dataMethods.AreYouSure_Save(() => {
                    closeBoSSSpad = true; 
                    mainWindow.close();
                });
            }
        });   
    }

    var userDatabase = await UserDatabase.load('./userData.txt');
    var recentDocuments = new RecentDocuments(userDatabase.getUserData());
    var dataMethods = new BoSSSDataMethods(mainWindow, recentDocuments);
    createMenu();
    addCloseAction();
}

module.exports = addFunctionality;