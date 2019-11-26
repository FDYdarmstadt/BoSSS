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

    // Laplace of kinetic energy
    // =========================

    public class KineticEnergyLaplaceInSpeciesBulk : ViscosityInSpeciesBulk_GradUTerm {

        public KineticEnergyLaplaceInSpeciesBulk(double penalty, double sw, IncompressibleMultiphaseBoundaryCondMap bcMap,
            string spcName, SpeciesId spcId, int _D, double _muA, double _muB, double _betaS = 0.0)
            : base(penalty, sw, bcMap, spcName, spcId, 0, _D, _muA, _muB, _betaS) {
        }

        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.KineticEnergy };
            }
        }

        protected override double g_Diri(double[] X, double time, int EdgeTag, int d) {

            double kinE_Diri = 0.0;
            for (int i = 0; i < m_D; i++) {
                Func<double[], double, double> boundVel = this.velFunction[i][EdgeTag];
                kinE_Diri += boundVel(X, time) * boundVel(X, time);
            }

            double ret = kinE_Diri / 2.0;

            return ret;

        }

    }

    public class KineticEnergyLaplaceAtInterface : ILevelSetForm, ILevelSetEquationComponentCoefficient {

        LevelSetTracker m_LsTrk;

        public KineticEnergyLaplaceAtInterface(LevelSetTracker lstrk, double _muA, double _muB, double _penalty) {
            this.m_LsTrk = lstrk;
            this.muA = _muA;
            this.muB = _muB;
            this.penalty = _penalty;

        }

        double muA;
        double muB;

        double penalty;


        /// <summary>
        /// default-implementation
        /// </summary>
        public double LevelSetForm(ref CommonParamsLs inp,
        //public override double EdgeForm(ref Linear2ndDerivativeCouplingFlux.CommonParams inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            double[] N = inp.n;
            double hCellMin = this.m_LsTrk.GridDat.Cells.h_min[inp.jCell];

            int D = N.Length;
            //Debug.Assert(this.ArgumentOrdering.Count == D);
            Debug.Assert(Grad_uA.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uB.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uA.GetLength(1) == D);
            Debug.Assert(Grad_uB.GetLength(1) == D);

            double Grad_uA_xN = 0, Grad_uB_xN = 0, Grad_vA_xN = 0, Grad_vB_xN = 0;
            for (int d = 0; d < D; d++) {
                Grad_uA_xN += Grad_uA[0, d] * N[d];
                Grad_uB_xN += Grad_uB[0, d] * N[d];
                Grad_vA_xN += Grad_vA[d] * N[d];
                Grad_vB_xN += Grad_vB[d] * N[d];
            }

            double PosCellLengthScale = PosLengthScaleS[inp.jCell];
            double NegCellLengthScale = NegLengthScaleS[inp.jCell];

            double hCutCellMin = Math.Min(NegCellLengthScale, PosCellLengthScale);
            Debug.Assert(!(double.IsInfinity(hCutCellMin) || double.IsNaN(hCutCellMin)));

            if (hCutCellMin <= 1.0e-10 * hCellMin)
                // very small cell -- clippling
                hCutCellMin = hCellMin;

            double Ret = 0.0;

            Ret -= 0.5 * (muA * Grad_uA_xN + muB * Grad_uB_xN) * (vA - vB);                           // consistency term
            Ret -= 0.5 * (muA * Grad_vA_xN + muB * Grad_vB_xN) * (uA[0] - uB[0]);     // symmetry term
            Ret += (penalty / hCutCellMin) * (uA[0] - uB[0]) * (vA - vB) * (Math.Abs(muA) > Math.Abs(muB) ? muA : muB); // penalty term


            Debug.Assert(!(double.IsInfinity(Ret) || double.IsNaN(Ret)));
            return Ret;
        }

        MultidimensionalArray PosLengthScaleS;
        MultidimensionalArray NegLengthScaleS;

        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            NegLengthScaleS = csA.CellLengthScales;
            PosLengthScaleS = csB.CellLengthScales;
        }



        public int LevelSetIndex {
            get { return 0; }
        }

        public IList<string> ArgumentOrdering {
            get { return new string[] { VariableNames.KineticEnergy }; }
        }

        public SpeciesId PositiveSpecies {
            get { return m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        } 
    }


    // Divergence of stress tensor
    // ===========================


    public class StressDivergenceInSpeciesBulk : LinearFlux, ISpeciesFilter {

        /// <summary>
        /// Spatial dimension
        /// </summary>
        int m_D;

        double mu;

        public StressDivergenceInSpeciesBulk(int SpatDim, IncompressibleMultiphaseBoundaryCondMap bcMap, SpeciesId spcId, double _mu) {
            m_D = SpatDim;
            this.m_bcMap = bcMap;
            m_spcId = spcId;
            mu = _mu;
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
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.VelocityX_GradientVector(), VariableNames.VelocityY_GradientVector()); //, VariableNames.Pressure);
            }
        }

        static double[,] VelociytGradient(double[] GradVelX, double[] GradVelY) {
            Debug.Assert(GradVelX.Length == 2);
            Debug.Assert(GradVelY.Length == 2);

            int D = GradVelX.Length;
            double[,] GradVel = new double[D, D];

            GradVel[0, 0] = GradVelX[0];
            GradVel[0, 1] = GradVelX[1];
            GradVel[1, 0] = GradVelY[0];
            GradVel[1, 1] = GradVelY[1];

            return GradVel;

        }


        protected override void Flux(ref CommonParamsVol inp, Double[] U, Double[] output) {

            double[] Vel = inp.Parameters.GetSubVector(0, m_D);
            double[,] GradVel = VelociytGradient(inp.Parameters.GetSubVector(m_D, m_D), inp.Parameters.GetSubVector(2 * m_D, m_D));
            //double Press = inp.Parameters[3 * m_D];

            for (int d = 0; d < m_D; d++) {
                output[d] = 0; // -Press * Vel[d];        // pressure term
                for (int dd = 0; dd < m_D; dd++) {
                    output[d] += mu * GradVel[d, dd] * Vel[dd];      // velocity grad
                    //output[d] += mu * GradVel[dd, d] * Vel[dd];      // velocity grad transposed term
                }
            }
            output.ScaleV(-1.0);
        }


        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, Double[] Uin) {

            double[] Vel_IN = inp.Parameters_IN.GetSubVector(0, m_D);
            double[,] GradVel_IN = VelociytGradient(inp.Parameters_IN.GetSubVector(m_D, m_D), inp.Parameters_IN.GetSubVector(2 * m_D, m_D));
            //double Press_IN = inp.Parameters_IN[3 * m_D];

            double acc = 0;

            IncompressibleBcType edgType = m_bcMap.EdgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case IncompressibleBcType.Wall: {
                        for (int d = 0; d < m_D; d++) {
                            //acc -= Press_IN * Vel_IN[d] * inp.Normale[d];
                            for (int dd = 0; dd < m_D; dd++) {
                                acc += mu * (GradVel_IN[d, dd] * Vel_IN[dd]) * inp.Normale[d];
                                //acc += mu * (GradVel_IN[dd, d] * Vel_IN[dd]) * inp.Normale[d];  // transposed term
                            }
                        }
                        break;
                    }
                case IncompressibleBcType.Velocity_Inlet: {
                        for (int d = 0; d < m_D; d++) {
                            //acc -= Press_IN * Vel_IN[d] * inp.Normale[d];
                            for (int dd = 0; dd < m_D; dd++) {
                                acc += mu * (GradVel_IN[d, dd] * Vel_IN[dd]) * inp.Normale[d];
                                //acc += mu * (GradVel_IN[dd, d] * Vel_IN[dd]) * inp.Normale[d];  // transposed term
                            }
                        }
                        break;
                    }
                case IncompressibleBcType.Pressure_Outlet:
                case IncompressibleBcType.Pressure_Dirichlet: {
                        for (int d = 0; d < m_D; d++) {
                            //acc -= Press_IN * Vel_IN[d] * inp.Normale[d];
                            for (int dd = 0; dd < m_D; dd++) {
                                acc += mu * (GradVel_IN[d, dd] * Vel_IN[dd]) * inp.Normale[d];
                                //acc += mu * (GradVel_IN[dd, d] * Vel_IN[dd]) * inp.Normale[d];  // transposed term
                            }
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
            double[,] GradVel_IN = VelociytGradient(inp.Parameters_IN.GetSubVector(m_D, m_D), inp.Parameters_IN.GetSubVector(2 * m_D, m_D));
            double[,] GradVel_OUT = VelociytGradient(inp.Parameters_OUT.GetSubVector(m_D, m_D), inp.Parameters_OUT.GetSubVector(2 * m_D, m_D));
            //double Press_IN = inp.Parameters_IN[3 * m_D];
            //double Press_OUT = inp.Parameters_OUT[3 * m_D];

            double acc = 0;

            for (int d = 0; d < m_D; d++) {
                //acc -= 0.5 * (Press_IN * Vel_IN[d] + Press_OUT * Vel_OUT[d]) * inp.Normale[d];
                for (int dd = 0; dd < m_D; dd++) {
                    acc += 0.5 * mu * (GradVel_IN[d, dd] * Vel_IN[dd] + GradVel_OUT[d, dd] * Vel_OUT[dd]) * inp.Normale[d];
                    //acc += 0.5 * mu * (GradVel_IN[dd, d] * Vel_IN[dd] + GradVel_OUT[dd, d] * Vel_OUT[dd]) * inp.Normale[d];  // transposed term
                }
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


    public class StressDivergenceAtLevelSet : ILevelSetForm {

        LevelSetTracker m_LsTrk;

        public StressDivergenceAtLevelSet(LevelSetTracker lstrk, double _muA, double _muB) {
            this.m_LsTrk = lstrk;
            this.muA = _muA;
            this.muB = _muB;
            this.m_D = lstrk.GridDat.SpatialDimension;
        }

        double muA;
        double muB;

        int m_D;



        static double[,] VelociytGradient(double[] GradVelX, double[] GradVelY) {
            Debug.Assert(GradVelX.Length == 2);
            Debug.Assert(GradVelY.Length == 2);

            int D = GradVelX.Length;
            double[,] GradVel = new double[D, D];

            GradVel[0, 0] = GradVelX[0];
            GradVel[0, 1] = GradVelX[1];
            GradVel[1, 0] = GradVelY[0];
            GradVel[1, 1] = GradVelY[1];

            return GradVel;

        }


        public Double LevelSetForm(ref CommonParamsLs inp, Double[] uA, Double[] uB, Double[,] Grad_uA, Double[,] Grad_uB, Double vA, Double vB, Double[] Grad_vA, Double[] Grad_vB) {

            double[] Vel_A = inp.ParamsNeg.GetSubVector(0, m_D);
            double[] Vel_B = inp.ParamsPos.GetSubVector(0, m_D);
            double[,] GradVel_A = VelociytGradient(inp.ParamsNeg.GetSubVector(m_D, m_D), inp.ParamsNeg.GetSubVector(2 * m_D, m_D));
            double[,] GradVel_B = VelociytGradient(inp.ParamsPos.GetSubVector(m_D, m_D), inp.ParamsPos.GetSubVector(2 * m_D, m_D));
            //double p_A = inp.ParamsNeg[3 * m_D];
            //double p_B = inp.ParamsPos[3 * m_D];

            double ret = 0.0;

            for (int d = 0; d < m_D; d++) {
                //ret += 0.5 * (p_A * Vel_A[d] + p_B * Vel_B[d]) * inp.n[d];  // pressure
                for (int dd = 0; dd < m_D; dd++) {
                    ret -= 0.5 * (muA * GradVel_A[d, dd] * Vel_A[dd] + muB * GradVel_B[d, dd] * Vel_B[dd]) * inp.n[d];  // gradU
                    //ret -= 0.5 * (muA * GradVel_A[dd, d] * Vel_A[dd] + muB * GradVel_B[dd, d] * Vel_B[dd]) * inp.n[d];  // gradU transposed
                }
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
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.VelocityX_GradientVector(), VariableNames.VelocityY_GradientVector()); //, VariableNames.Pressure);
            }
        }


    }


    // Dissipation
    // ===========


    public class Dissipation : IVolumeForm, ISpeciesFilter {


        /// <summary>
        /// Spatial dimension;
        /// </summary>
        int m_D;

        double mu;

        public Dissipation(int SpatDim, double _mu, SpeciesId spcId) {
            m_D = SpatDim;
            mu = _mu;
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
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.VelocityX_GradientVector(), VariableNames.VelocityY_GradientVector());
            }
        }


        static double[,] VelociytGradient(double[] GradVelX, double[] GradVelY) {
            Debug.Assert(GradVelX.Length == 2);
            Debug.Assert(GradVelY.Length == 2);

            int D = GradVelX.Length;
            double[,] GradVel = new double[D, D];

            GradVel[0, 0] = GradVelX[0];
            GradVel[0, 1] = GradVelX[1];
            GradVel[1, 0] = GradVelY[0];
            GradVel[1, 1] = GradVelY[1];

            return GradVel;

        }


        public double VolumeForm(ref CommonParamsVol cpv, Double[] U, Double[,] GradU, Double V, Double[] GradV) {

            double[] Vel = cpv.Parameters.GetSubVector(0, m_D);
            double[,] GradVel = VelociytGradient(cpv.Parameters.GetSubVector(m_D, m_D), cpv.Parameters.GetSubVector(2 * m_D, m_D));

            double ret = 0;

            for (int d = 0; d < m_D; d++) {
                for (int dd = 0; dd < m_D; dd++) {
                    ret += GradVel[d, dd] * GradVel[dd, d];
                    ret += GradVel[dd, d] * GradVel[dd, d];  // transposed term
                }
            }

            return mu * ret * V;
        }


        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.V;
            }
        }


    }



}
