using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ilPSP.Connectors {

    /// <summary>
    /// Has to be implemented by classes which should be available through external 
    /// language bindings to BoSSS.
    /// </summary>
    public interface IForeignLanguageProxy {

        /// <summary>
        /// Stores a reference/pointer/handle to some proxy object in a foreign language
        /// </summary>
        /// <param name="ptr">
        /// pointer, reference, etc. to a proxy object in a foreign language
        /// </param>
        void _SetForeignPointer(IntPtr ptr);


        /// <summary>
        /// Getter
        /// </summary>
        IntPtr _GetForeignPointer();

    }
}
