const BoSSSMenu = require('./boSSSMenu.js');
const BoSSSDataMethods = require('./dataMethods.js');
const CommandActions = require('./commandActions.js');
const UserData = require('./commandActions.js');

function addFunctionality(mainWindow){
    
    function createMenu(dataMethods){
        var commandActions = new CommandActions(mainWindow);
        var userData = new UserData();
        var menu = new BoSSSMenu(commandActions, dataMethods);
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
    createMenu(dataMethods);
    addCloseAction(dataMethods);
    var closeBoSSSpad = false;
}

module.exports = addFunctionality;