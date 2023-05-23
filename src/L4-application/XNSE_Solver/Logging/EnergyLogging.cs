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

using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.Logging {

    /// <summary>
    /// in-situ post-processing for <see cref="CapillaryRise"/>
    /// </summary>
    [Serializable]
    public class EnergyLogging : EnergyLogging<XNSE_Control> {
    }

    /// <summary>
    /// For single droplet/bubble computations with origin in (0,0,0):
    /// extraction of the major axis lengths (theta = 0 and theta = 90, where theta describes the inclination angle) and volume  
    /// </summary>
    [Serializable]
    public class EnergyLogging<T> : XNSEinSituPostProcessingModule<T> where T : XNSE_Control, new() {

        /// <summary>
        /// filename
        /// </summary>
        protected override string LogFileName => "EnergyLogValues";


        /// <summary>
        /// CSV header
        /// </summary>
        protected override void WriteHeader(TextWriter textWriter) {

            string header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", "#timestep", "time", "kineticEnergy", "surfaceEnergy", "totalEnergy", "kineticDissipation");
            Log.WriteLine(header);

        }


        /// <summary>
        /// 
        /// </summary>
        protected override void PerformTimestepPostProcessing(int TimestepNo, double physTime) {
            using(new FuncTrace()) {

                double[] rhoS = new double[] {Control.PhysicalParameters.rho_A, Control.PhysicalParameters.rho_B };
                int quadOrder = CurrentVel[0].Basis.Degree * CurrentVel[0].Basis.Degree;
                double kineticEnergy = EnergyUtils.GetKineticEnergy(LsTrk, CurrentVel, rhoS, quadOrder);

                double sigma = Control.PhysicalParameters.Sigma;
                double surfaceEnergy = EnergyUtils.GetSurfaceEnergy(LsTrk, sigma, quadOrder);

                double totalEnergy = kineticEnergy + surfaceEnergy;

                double[] muS = new double[] { Control.PhysicalParameters.mu_A, Control.PhysicalParameters.mu_B };
                double kineticDissipation = EnergyUtils.GetKineticDissipation(LsTrk, CurrentVel, muS, quadOrder);

                string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}\t{5}", TimestepNo, physTime, kineticEnergy, surfaceEnergy, totalEnergy, kineticDissipation);
                Log.WriteLine(line);
                Log.Flush();
            }
        }
    }

}
