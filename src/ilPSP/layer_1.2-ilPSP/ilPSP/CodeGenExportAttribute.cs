using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ilPSP.Connectors {
    /// <summary>
    /// Attribute to mark functions that should be mapped in external language proxies (<see cref="IForeignLanguageProxy"/>
    /// </summary>
    public class CodeGenExportAttribute : Attribute {
    }
}
