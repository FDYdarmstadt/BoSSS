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

            // if we enter Main, it seems we are executing the tutorial tests locally...
            // so delete any local tex files since we want to run the scripts from doc/handbook
            var localTexFiles = (new DirectoryInfo(Directory.GetCurrentDirectory())).GetFiles("*.tex");
            foreach (var f in localTexFiles) {
                f.Delete();
            }


            BoSSS.Solution.Application.InitMPI(new string[0]);
            
            // start the minibatchprocessor which is used internally
            bool iStartedThisShit = AllUpTest.OneTimeSetUp();

            var losScriptos = GetListOfScripts();
            int r = 0;
            int i = 1;
            foreach(var s in losScriptos) {
                AllUpTest.RunWorksheet(s);
                //Console.WriteLine($"{i} : {s}");
                i++;
            }

            

            
            //var tr = new TextRunner(typeof(TutorialTestsMain).Assembly);
            
            //int r = tr.Execute(new[] { "--result=result-TutorialTests.xml"
            //    //, "--test=BoSSS.Application.TutorialTests.AllUpTest.Run__BoundaryAndInitialData" 
            //});
            
/*
            int r = 0;
            AllUpTest.Run__BoundaryAndInitialData();
            AllUpTest.Run__InitialValues();
            AllUpTest.Run__channel();
            AllUpTest.Run__GridGeneration();
            AllUpTest.Run__IsentropicVortex();
            AllUpTest.Run__MetaJobManager();
            AllUpTest.Run__tutorialMatlab();
            AllUpTest.Run__ue2Basics();
#if !DEBUG
            AllUpTest.Run__CsharpAndBoSSSpad();
            AllUpTest.Run__convStudy();
            AllUpTest.Run__Poisson();
            AllUpTest.Run__sip();
            AllUpTest.Run__StokesEq();
            AllUpTest.Run__SpatialOperatorNexpTimeInt();
            AllUpTest.Run__ue6ScalarConvStability();
            AllUpTest.Run__ue5NumFluxConv();
#endif
*/
            AllUpTest.OneTimeTearDown(iStartedThisShit);
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
    
