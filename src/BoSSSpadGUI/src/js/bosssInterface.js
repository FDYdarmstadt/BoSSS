var edge = require('electron-edge-js');
var path = require('path');

//Singleton Approach
class BoSSS{
    constructor(){
        if(! this.constructor.prototype.instance ){
            var BoSSS_DLL_path = path.join(__dirname, '../src/cs/bin/Debug/ElectronWorksheet.dll');
            var requireBoSSS = edge.func(
            {
                assemblyFile: BoSSS_DLL_path,
                typeName: 'BoSSS.Application.BoSSSpad.ElectronInterface',
                methodName: 'Invoke', // This must be Func<object,Task<object>>
            });
            var BoSSS_Dir = path.join(__dirname, '../src/cs/bin/Debug/');
            this.BoSSSRuntime = requireBoSSS( BoSSS_Dir, true);
            this.constructor.prototype.instance = this;
        }
        return this.constructor.prototype.instance;
    }

    provideAutoComplete(valueString){
        var that = this;
        var runPromise = new Promise(
            function(resolve, reject){
                that.BoSSSRuntime.getAutoCompleteSuggestions(
                    valueString,
                    async function(error, result) {
                        if (error){
                            reject(error);
                        }
                        resolve(result);
                    }
                );
            }
        );
        return runPromise;
    }

    runCommand(commandString){
        var that = this;
        var runPromise = new Promise(
            function(resolve, reject){
                that.BoSSSRuntime.runCommand(
                    commandString,
                    async function(error, result) {
                        if (error){
                            reject(error);
                        }
                        resolve(result);
                    }
                );
            }
        );
        return runPromise;
    }
    save(data){
        //data {path: string, commands: string[], results: string[]}
        var that = this;
        var runPromise = new Promise(
            function(resolve, reject){
                that.BoSSSRuntime.save(
                    data,
                    async function(error, result) {
                        if (error){
                            reject(error);
                        }
                        resolve(result);
                    }
                );
            }
        );
        return runPromise;
    }

    load(path){
        var that = this;
        var runPromise = new Promise(
            function(resolve, reject){
                that.BoSSSRuntime.load(
                    path,
                    async function(error, result) {
                        if (error){
                            reject(error);
                        }
                        resolve(result);
                    }
                );
            }
        );
        return runPromise;
    }

}

const instance = new BoSSS();
Object.freeze(instance);

export default instance;






