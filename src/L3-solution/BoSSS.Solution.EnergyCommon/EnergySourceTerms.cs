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
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.XNSECommon.Operator.Viscosity;
using ilPSP;
using ilPSP.Utils;

namespace BoSSS.Solution.EnergyCommon {


    public class PowerofGravity : IVolumeForm, ISpeciesFilter {


        /// <summary>
        /// Spatial dimension;
        /// </summary>
        int m_D;

        double rho;


        public PowerofGravity(int SpatDim, SpeciesId spcId, double _rho) {
            m_D = SpatDim;
            m_spcId = spcId;
            rho = _rho;
        }


        SpeciesId m_spcId;

        public SpeciesId validSpeciesId {
            get { return m_spcId; }
        }


        public IList<String> ArgumentOrdering {
            get {
                return new string[] { };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.GravityVector(m_D));
            }
        }


        public double VolumeForm(ref CommonParamsVol cpv, Double[] U, Double[,] GradU, Double V, Double[] GradV) {

            double[] Vel = cpv.Parameters.GetSubVector(0, m_D);
            double[] Grav = cpv.Parameters.GetSubVector(m_D, m_D);

            double ret = 0;

            for (int d = 0; d < m_D; d++) {
                ret += Grav[d] * Vel[d];
            }

            return -rho * ret * V;
        }


        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.V;
            }
        }


    }



    public class SurfaceEnergy : ILevelSetForm {

        LevelSetTracker m_LsTrk;

        public SurfaceEnergy(int _D, LevelSetTracker LsTrk, double _sigma) {
            m_LsTrk = LsTrk;
            this.m_D = _D;
            this.sigma = _sigma;
        }

        int m_D;

        double sigma;


        public double LevelSetForm(ref CommonParamsLs cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            //Debug.Assert(cp.ParamsPos[0] == cp.ParamsNeg[0], "interface velocityX must be continuous across interface");
            //Debug.Assert(cp.ParamsPos[1] == cp.ParamsNeg[1], "interface velocityY must be continuous across interface");
            //Debug.Assert(cp.ParamsPos[2] == cp.ParamsNeg[2], "curvature must be continuous across interface");

            double curvature = cp.ParamsPos[m_D];
            double[] Vel = cp.ParamsPos.GetSubVector(0, m_D);
            double[] Normal = cp.n;

            double surfE = 0;
            for (int d = 0; d < m_D; d++) {
                surfE -= curvature * sigma * (Vel[d] * Normal[d]);
            }

            double FlxNeg = -0.5 * surfE;
            double FlxPos = +0.5 * surfE;

            Debug.Assert(!(double.IsNaN(FlxNeg) || double.IsInfinity(FlxNeg)));
            Debug.Assert(!(double.IsNaN(FlxPos) || double.IsInfinity(FlxPos)));

            return FlxNeg * vA - FlxPos * vB;
        }

        public IList<string> ArgumentOrdering {
            get {
                return new string[] { };
            }
        }


        public IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0MeanVector(m_D), VariableNames.Curvature);
            }
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public SpeciesId PositiveSpecies {
            get { return this.m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return this.m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get { return TermActivationFlags.V; }
        }


    }



}
