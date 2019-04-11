const BoSSSMenu = require('./boSSSMenu.js');
const BoSSSDataMethods = require('./dataMethods.js');
const CommandActions = require('./commandActions.js');
const UserData = require('./UserData/userData.js');
const RecentDocuments = require('./recentDocuments.js')

function addFunctionality(mainWindow){
    
    function createMenu(dataMethods, userData){
        var commandActions = new CommandActions(mainWindow);
        var userData = new RecentDocuments(userData);
        var menu = new BoSSSMenu(commandActions, dataMethods, userData);
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
    var userData = UserData.build();
    createMenu(dataMethods, userData);
    addCloseAction(dataMethods);
    var closeBoSSSpad = false;
}

module.exports = addFunctionality;