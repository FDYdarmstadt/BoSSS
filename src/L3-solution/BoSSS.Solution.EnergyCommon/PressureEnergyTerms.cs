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
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.XNSECommon.Operator.Viscosity;
using ilPSP;
using ilPSP.Utils;

namespace BoSSS.Solution.EnergyCommon {

    public class DivergencePressureEnergyInSpeciesBulk : LinearFlux, ISpeciesFilter {

        /// <summary>
        /// Spatial dimension;
        /// </summary>
        int m_D;


        /// <summary>
        /// Mapping from edge tags to boundary values.<br/>
        /// 1st index: edge tag;<br/>
        /// 2nd index: spatial direction
        /// </summary>
        protected Func<double[], double, double>[,] VelocFunction;

        protected Func<double[], double, double>[] PressFunction;


        public DivergencePressureEnergyInSpeciesBulk(int SpatDim, IncompressibleMultiphaseBoundaryCondMap _bcmap, string spcName, SpeciesId spcId) {
            m_D = SpatDim;
            m_bcMap = _bcmap;
            m_spcId = spcId;


            VelocFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDim];
            for (int d = 0; d < m_D; d++)
                VelocFunction.SetColumn(m_bcMap.bndFunction[VariableNames.Velocity_d(d) + "#" + spcName], d);

            PressFunction = m_bcMap.bndFunction[VariableNames.Pressure + "#" + spcName];
        }


        IncompressibleMultiphaseBoundaryCondMap m_bcMap;


        SpeciesId m_spcId;

        public SpeciesId validSpeciesId {
            get { return m_spcId; }
        }


        public override IList<String> ArgumentOrdering {
            get {
                return new string[] { };
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.Pressure0);
            }
        }


        protected override void Flux(ref CommonParamsVol inp, Double[] U, Double[] output) {

            double[] Vel = inp.Parameters.GetSubVector(0, m_D);
            double Press = inp.Parameters[m_D];

            for (int d = 0; d < m_D; d++) {
                output[d] = Press * Vel[d];        // pressure term
            }
            //output.ScaleV(-1.0);
        }


        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, Double[] Uin) {

            double[] Vel_IN = inp.Parameters_IN.GetSubVector(0, m_D);
            double Press_IN = inp.Parameters_IN[m_D];

            double acc = 0;

            IncompressibleBcType edgType = m_bcMap.EdgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case IncompressibleBcType.Wall: {
                        //for (int d = 0; d < m_D; d++) {
                        //    acc -= Press_IN * Vel_IN[d] * inp.Normale[d];
                        //}
                        break;
                    }
                case IncompressibleBcType.Velocity_Inlet: {
                        for (int d = 0; d < m_D; d++) {
                            double VelD = VelocFunction[inp.EdgeTag,d](inp.X, inp.time);
                            acc -= Press_IN * VelD * inp.Normale[d];
                        }
                        break;
                    }
                case IncompressibleBcType.Pressure_Outlet: {
                        double pD = PressFunction[inp.EdgeTag](inp.X, inp.time);
                        for (int d = 0; d < m_D; d++) {
                            acc -= pD * Vel_IN[d] * inp.Normale[d];
                        }
                        break;
                    }
                default: {
                        throw new NotImplementedException("ToDo");
                    }
            }

            return -acc;

        }

        protected override double InnerEdgeFlux(ref CommonParams inp, Double[] Uin, Double[] Uout) {

            double[] Vel_IN = inp.Parameters_IN.GetSubVector(0, m_D);
            double[] Vel_OUT = inp.Parameters_OUT.GetSubVector(0, m_D);
            double Press_IN = inp.Parameters_IN[m_D];
            double Press_OUT = inp.Parameters_OUT[m_D];

            double acc = 0;

            for (int d = 0; d < m_D; d++) {
                acc -= 0.5 * (Press_IN * Vel_IN[d] + Press_OUT * Vel_OUT[d]) * inp.Normale[d];
            }

            return -acc;
        }


        override public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.GradV;
            }
        }

        override public TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.V;
            }
        }

        override public TermActivationFlags InnerEdgeTerms {
            get {
                return TermActivationFlags.V;
            }
        }
    }


    public class DivergencePressureEnergyAtLevelSet : ILevelSetForm {

        LevelSetTracker m_LsTrk;

        public DivergencePressureEnergyAtLevelSet(LevelSetTracker lstrk) {
            this.m_LsTrk = lstrk;
            this.m_D = lstrk.GridDat.SpatialDimension;
        }

        int m_D;


        public Double LevelSetForm(ref CommonParamsLs inp, Double[] uA, Double[] uB, Double[,] Grad_uA, Double[,] Grad_uB, Double vA, Double vB, Double[] Grad_vA, Double[] Grad_vB) {

            double[] Vel_A = inp.ParamsNeg.GetSubVector(0, m_D);
            double[] Vel_B = inp.ParamsPos.GetSubVector(0, m_D);
            double p_A = inp.ParamsNeg[m_D];
            double p_B = inp.ParamsPos[m_D];

            double ret = 0.0;

            for (int d = 0; d < m_D; d++) {
                ret += 0.5 * (p_A * Vel_A[d] + p_B * Vel_B[d]) * inp.n[d];  // pressure
            }

            ret *= (vA - vB);

            return ret;
        }


        public int LevelSetIndex {
            get { return 0; }
        }

        public IList<string> ArgumentOrdering {
            get { return new string[] { }; }
        }

        public SpeciesId PositiveSpecies {
            get { return m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.V;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.Pressure0);
            }
        }
    }


    public class PressureGradientConvection : IVolumeForm, ISpeciesFilter {


        /// <summary>
        /// Spatial dimension;
        /// </summary>
        int m_D;


        public PressureGradientConvection(int SpatDim, SpeciesId spcId) {
            m_D = SpatDim;
            m_spcId = spcId;
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
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.PressureGradient(m_D));
            }
        }


        public double VolumeForm(ref CommonParamsVol cpv, Double[] U, Double[,] GradU, Double V, Double[] GradV) {

            double[] Vel = cpv.Parameters.GetSubVector(0, m_D);
            double[] PressGrad = cpv.Parameters.GetSubVector(m_D, m_D);

            double ret = 0;

            for (int d = 0; d < m_D; d++) {
                ret += PressGrad[d] * Vel[d];
            }

            return ret * V;
        }


        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.V;
            }
        }


    }


}
