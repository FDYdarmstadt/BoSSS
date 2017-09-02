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
using ilPSP;
using ilPSP.Connectors.Matlab;

namespace BoSSS.Application.TutorialTests {

    static class TutorialTestsMain {
        
        static void Main(string[] args) {
            AllUpTest.TestFixtureSetUp();
            AllUpTest.RunWorksheets("quickStartCNS/IsentropicVortex.tex");
            AllUpTest.RunWorksheets("MetaJobManager/MetaJobManager.tex");
            //AllUpTest.RunWorksheets("quickStartIBM/channel.tex");
            //AllUpTest.RunWorksheets("tutorial2/uebung2tutorial.tex");
            //AllUpTest.RunWorksheets("tutorial4/tutorial4.tex");
            //AllUpTest.RunWorksheets("tutorial5/uebung5tutorial.tex");
            //AllUpTest.RunWorksheets("tutorial6/tutorial6.tex");
            //AllUpTest.RunWorksheets("tutorial9-SIP/sip.tex");
            //AllUpTest.RunWorksheets("tutorial10-PoissonSystem/Poisson.tex");
            //AllUpTest.RunWorksheets("tutorial11-Stokes/StokesEq.tex");
            AllUpTest.TestFixtureTearDown();
        }

    }
}
