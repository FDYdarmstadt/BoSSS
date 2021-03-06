﻿using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using System.Diagnostics;
using BoSSS.Foundation.IO;

namespace BoSSS.Application.BoSSSpad{

    /// <summary>
    /// Driver class for <see cref="ElectronWorksheet"/>, used by the Electron-GUI
    /// </summary>
    public class ElectronInterface{

        static ElectronWorksheet worksheet;
        static CancellationTokenSource runCommandManager = null;

        /// <summary>
        /// Entrypoint of Electron BoSSSpad, i.e. the electron-edge-js package
        /// </summary>
        /// <param name="input">
        /// Path to the ElectronWorksheet.dll, ElectronBoSSSpad.exe and affiliated DLLs
        /// </param>
        /// <returns></returns>
        public async Task<object> Invoke(object input){
            worksheet = new ElectronWorksheet(input.ToString());
            return new{
                runCommand = (Func<object, Task<object>>)(async (i) => 
                {
                    runCommandManager = new CancellationTokenSource();
                    return await Task.Run(() => ElectronInterface.RunCommand(i), runCommandManager.Token);
                }),
                save = (Func<object, Task<object>>)(async (i) => 
                {
                    return await Task.Run(() => ElectronInterface.Save(i));
                }),
                load = (Func<object, Task<object>>)(async (i) => 
                {
                    return await Task.Run(() => ElectronInterface.Load(i));
                }),
                getAutoCompleteSuggestions = (Func<object, Task<object>>)(async (i) => 
                {
                    return await Task.Run(() => ElectronInterface.GetAutoCompleteSuggestions(i));
                }),
                forceAbort = (Func<object, Task<object>>)(async (i) =>
                {
                    return await Task.Run(() => ElectronInterface.ForceAbort());
                })
            };
        }

        static object RunCommand(object input){
            Tuple<string, string> output = worksheet.RunCommand(input.ToString());
            return output;
        }

        static bool Save(dynamic input){
            string path = (string)input.path;
            object[] commands = (object[])input.commands;
            object[] results = (object[])input.results;

            string[] stringCommands = new string[commands.Length];
            for(int i = 0; i < commands.Length; ++i){
                stringCommands[i] = commands[i].ToString();
            }
            string[] stringResults = new string[results.Length];
            for (int i = 0; i < results.Length; ++i)
            {
                stringResults[i] = results[i].ToString();
            }
            worksheet.Save(path, stringCommands, stringResults);
            return true;
        }

        static object Load(object path){
            return worksheet.Load((string)path);
        }

        static object GetAutoCompleteSuggestions(object input){
            return worksheet.GetAutoCompleteSuggestions(input.ToString());
        }

        static object ForceAbort()
        {
            if (runCommandManager != null)
                runCommandManager.Cancel();
            return null;
        }
    }
}
