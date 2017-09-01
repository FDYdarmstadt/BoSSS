/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using System;
using System.Windows.Forms;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Creates a WindowsForm to plot data. Should be started as a new Thread. 
    /// </summary>
    /// 
    public static class AutonomuousPlotter {

        /// <summary>
        /// The thread of the Windows-Form. Depending on the passed arguments,
        /// a live, or a static plotting form is started.  
        /// </summary>
        /// <param name="formFunction"></param>
        [STAThread]
        public static void DisplayWindow(Func<Form> formFunction) {
            System.Windows.Forms.Application.EnableVisualStyles();
            System.Windows.Forms.Application.SetCompatibleTextRenderingDefault(false);
            System.Windows.Forms.Application.Run(formFunction());
        }  

    }
}
