using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using System.Diagnostics;
using BoSSS.Foundation.IO;

namespace BoSSS.Application.BoSSSpad{

    public class ElectronInterface{

        static ElectronWorksheet worksheet;

        public async Task<object> Invoke(object input){
            worksheet = ElectronWorksheet.Instance;
            return new{
                runCommand = (Func<object, Task<object>>)(async (i) => {
                    return await Task.Run(() => ElectronInterface.RunCommand(i));
                }),

            };
        }

        static string RunCommand(object input){
            string output = worksheet.RunCommand(input.ToString());
            return output;
        }
    }
}
