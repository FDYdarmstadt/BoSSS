const BoSSSMenu = require('./boSSSMenu.js');
const DataMethods = require('./recentDataMethods.js');
const CommandMethods = require('./commandActions.js');
const UserDatabase = require('./UserData/userDatabase.js');
const RecentDocuments = require('./recentDocuments.js');

async function addFunctionality(mainWindow){
    
    var userDatabase = await UserDatabase.load('./userData.txt');
    var recentDocuments = new RecentDocuments(userDatabase.getUserData());
    var dataMethods = new DataMethods(mainWindow, recentDocuments);

    function createMenu(){
        var commandMethods = new CommandMethods(mainWindow);
        var menu = new BoSSSMenu(commandMethods, dataMethods, recentDocuments);
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
    createMenu();
    addCloseAction();
}

module.exports = addFunctionality;