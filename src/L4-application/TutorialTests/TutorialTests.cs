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
using System.Drawing;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Threading.Tasks;
using ilPSP;
using ilPSP.Connectors.Matlab;
using MPI.Wrappers;
using NUnit.Framework;
using NUnitLite;

namespace BoSSS.Application.TutorialTests {

    static class TutorialTestsMain {

        static int Main(string[] args) {
            /*
            BoSSS.Solution.Application.InitMPI(new string[0]);

 ;
            */
            

            AllUpTest.DirectoryOffset = Path.Combine("..", "..", "..", ".." ,"..", ".." , "doc", "handbook");
            
            if(!Directory.Exists(Path.Combine(Directory.GetCurrentDirectory(), AllUpTest.DirectoryOffset)))
                throw new IOException("not running in nexpected directory");

            // if we enter Main, it seems we are executing the tutorial tests locally...
            // so delete any local tex files since we want to run the scripts from doc/handbook
            foreach(string ext in new[] { "*.tex", "*.ipynb" }) {
                var localTexFiles = (new DirectoryInfo(Directory.GetCurrentDirectory())).GetFiles(ext);
                foreach(var f in localTexFiles) {
                    f.Delete();
                }
            }


            BoSSS.Solution.Application.InitMPI(new string[0]);

            BoSSS.Application.TutorialTests.WorkflowMgmTests.Run__WorkFlowManagerFunctionalityCheck_1();
              
            // start the minibatchprocessor which is used internally
            //bool iStartedThisShit = AllUpTest.OneTimeSetUp();

            var losScriptos = GetListOfScripts();
            int r = 0;
/*
            int i = 1;
            foreach(var s in losScriptos) {
                AllUpTest.RunWorksheet(s);
                //Console.WriteLine($"{i} : {s}");
                i++;
            }
            

            //AllUpTest.OneTimeTearDown(iStartedThisShit);
            */
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
    
