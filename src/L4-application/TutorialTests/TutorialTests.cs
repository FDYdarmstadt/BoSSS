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
using System.Collections.Generic;
using System.IO;
using System.Reflection;
using ilPSP;
using ilPSP.Connectors.Matlab;
using MPI.Wrappers;
using NUnit.Framework;
using NUnitLite;

namespace BoSSS.Application.TutorialTests {

    static class TutorialTestsMain {

        static int Main(string[] args) {
            AllUpTest.DirectoryOffset = Path.Combine("..", "..", "..", "..", "..", "doc", "handbook");
            BoSSS.Solution.Application.InitMPI(new string[0]);
            
          
            //var losScriptos = GetListOfScripts();
            //int r = 0;
            //AllUpTest.RunWorksheet(losScriptos[int.Parse(args[0])]);
            //foreach(var s in losScriptos) {
            //    AllUpTest.RunWorksheet(s);
            //}




            var tr = new TextRunner(typeof(TutorialTestsMain).Assembly);
            int r = tr.Execute(new[] { "--result=result-TutorialTests.xml" });


            csMPI.Raw.mpiFinalize();

            
            return r;
        }

        static string[] GetListOfScripts() {
            var r = new List<string>();
            

            var mmm = typeof(AllUpTest).GetMethods();

            foreach(var m in mmm) {

                if(m.GetCustomAttribute(typeof(TestAttribute)) != null) {
                    var dc = m.GetCustomAttribute(typeof(NUnitFileToCopyHackAttribute)) as NUnitFileToCopyHackAttribute;

                    if(dc != null) {
                        r.AddRange(dc.SomeFileNames);
                    }
                }
            }

            return r.ToArray();
        }
    }


}
    
