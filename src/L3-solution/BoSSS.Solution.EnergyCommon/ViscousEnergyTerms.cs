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

    // Laplace of kinetic energy
    // =========================

    public class KineticEnergyLaplaceInSpeciesBulk : ViscosityInSpeciesBulk_GradUTerm {

        public KineticEnergyLaplaceInSpeciesBulk(double penalty, double sw, IncompressibleMultiphaseBoundaryCondMap bcMap,
            string spcName, SpeciesId spcId, int _D, double _muA, double _muB, double _betaS = 0.0)
            : base(penalty, sw, bcMap, spcName, 0, _D, _muA, _muB, _betaS) {
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
        public double InnerEdgeForm(ref CommonParams inp,
        //public override double EdgeForm(ref Linear2ndDerivativeCouplingFlux.CommonParams inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            double[] N = inp.Normal;
            double hCellMin = this.m_LsTrk.GridDat.Cells.h_min[inp.jCellIn];

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

            double PosCellLengthScale = PosLengthScaleS[inp.jCellIn];
            double NegCellLengthScale = NegLengthScaleS[inp.jCellIn];

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

        public string PositiveSpecies {
            get { return "B"; }
        }

        public string NegativeSpecies {
            get { return "A"; }
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

        bool transposedTerm; 

        /// <summary>
        /// Mapping from edge tags to boundary values.<br/>
        /// 1st index: edge tag;<br/>
        /// 2nd index: spatial direction
        /// </summary>
        protected Func<double[], double, double>[,] VelFunction;

        public StressDivergenceInSpeciesBulk(int SpatDim, IncompressibleMultiphaseBoundaryCondMap bcMap, 
            string spcName, SpeciesId spcId, double _mu, bool transposed = false) {

            m_D = SpatDim;
            this.m_bcMap = bcMap;
            m_spcId = spcId;
            mu = _mu;
            ValidSpecies = spcName;

            transposedTerm = transposed;

            VelFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDim];
            for (int d = 0; d < m_D; d++)
                VelFunction.SetColumn(m_bcMap.bndFunction[VariableNames.Velocity_d(d) + "#" + spcName], d);
        }


        IncompressibleMultiphaseBoundaryCondMap m_bcMap;

        SpeciesId m_spcId;

        public string ValidSpecies {
            get;
            private set;
        }


        public override IList<String> ArgumentOrdering {
            get {
                return new string[] { };
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.VelocityVector(m_D), VariableNames.VelocityX_GradientVector(), VariableNames.VelocityY_GradientVector()); 
            }
        }

        static double[,] VelocityGradient(double[] GradVelX, double[] GradVelY) {
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


        protected override void Flux(ref CommonParamsVol inp, double[] U, double[] output) {

            double[] Vel = inp.Parameters.GetSubVector(0, m_D);
            double[,] GradVel = VelocityGradient(inp.Parameters.GetSubVector(m_D, m_D), inp.Parameters.GetSubVector(2 * m_D, m_D));

            for (int d = 0; d < m_D; d++) {
                output[d] = 0; 
                for (int dd = 0; dd < m_D; dd++) {
                    output[d] -= mu * GradVel[d, dd] * Vel[dd];      // velocity grad
                    if(transposedTerm)
                        output[d] -= mu * GradVel[dd, d] * Vel[dd];      // velocity grad transposed term
                }
            }

        }


        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {

            double[] Vel_IN = inp.Parameters_IN.GetSubVector(0, m_D);
            double[,] GradVel_IN = VelocityGradient(inp.Parameters_IN.GetSubVector(m_D, m_D), inp.Parameters_IN.GetSubVector(2 * m_D, m_D));

            double acc = 0;

            IncompressibleBcType edgType = m_bcMap.EdgeTag2Type[inp.EdgeTag];

            switch (edgType) {
                case IncompressibleBcType.Wall: {
                        break;
                    }
                case IncompressibleBcType.Velocity_Inlet: {
                        for (int d = 0; d < m_D; d++) {
                            for (int dd = 0; dd < m_D; dd++) {
                                double VelD = VelFunction[inp.EdgeTag, dd](inp.X, inp.time);
                                acc += mu * (GradVel_IN[d, dd] * VelD) * inp.Normal[d];
                                if (transposedTerm) {
                                    acc += mu * (GradVel_IN[dd, d] * VelD) * inp.Normal[d];  // transposed term
                                }
                            }
                        }
                        break;
                    }
                case IncompressibleBcType.Pressure_Outlet:
                case IncompressibleBcType.Pressure_Dirichlet: {
                        for (int d = 0; d < m_D; d++) {
                            for (int dd = 0; dd < m_D; dd++) {
                                acc += mu * (GradVel_IN[d, dd] * Vel_IN[dd]) * inp.Normal[d];
                                if (transposedTerm) {
                                    acc += mu * (GradVel_IN[dd, d] * Vel_IN[dd]) * inp.Normal[d];  // transposed term
                                }
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

        protected override double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout) {

            double[] Vel_IN = inp.Parameters_IN.GetSubVector(0, m_D);
            double[] Vel_OUT = inp.Parameters_OUT.GetSubVector(0, m_D);
            double[,] GradVel_IN = VelocityGradient(inp.Parameters_IN.GetSubVector(m_D, m_D), inp.Parameters_IN.GetSubVector(2 * m_D, m_D));
            double[,] GradVel_OUT = VelocityGradient(inp.Parameters_OUT.GetSubVector(m_D, m_D), inp.Parameters_OUT.GetSubVector(2 * m_D, m_D));

            double acc = 0;

            for (int d = 0; d < m_D; d++) {
                for (int dd = 0; dd < m_D; dd++) {
                    acc += 0.5 * mu * (GradVel_IN[d, dd] * Vel_IN[dd] + GradVel_OUT[d, dd] * Vel_OUT[dd]) * inp.Normal[d];
                    if (transposedTerm) {
                        acc += 0.5 * mu * (GradVel_IN[dd, d] * Vel_IN[dd] + GradVel_OUT[dd, d] * Vel_OUT[dd]) * inp.Normal[d];  // transposed term
                    }
                }
            }

            return -acc;
        }


        public override TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.GradV;
            }
        }

        public override TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.V;
            }
        }

        public override TermActivationFlags InnerEdgeTerms {
            get {
                return TermActivationFlags.V;
            }
        }

    }


    public class StressDivergenceAtLevelSet : ILevelSetForm {

        //LevelSetTracker m_LsTrk;

        public StressDivergenceAtLevelSet(LevelSetTracker lstrk, double _muA, double _muB, bool transposed = false) {
            //this.m_LsTrk = lstrk;
            this.muA = _muA;
            this.muB = _muB;
            this.m_D = lstrk.GridDat.SpatialDimension;

            transposedTerm = transposed;
        }

        double muA;
        double muB;

        int m_D;

        bool transposedTerm;


        static double[,] VelocityGradient(double[] GradVelX, double[] GradVelY) {
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


        public double InnerEdgeForm(ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            double[] Vel_A = inp.Parameters_IN.GetSubVector(0, m_D);
            double[] Vel_B = inp.Parameters_OUT.GetSubVector(0, m_D);
            double[,] GradVel_A = VelocityGradient(inp.Parameters_IN.GetSubVector(m_D, m_D), inp.Parameters_IN.GetSubVector(2 * m_D, m_D));
            double[,] GradVel_B = VelocityGradient(inp.Parameters_OUT.GetSubVector(m_D, m_D), inp.Parameters_OUT.GetSubVector(2 * m_D, m_D));
            //double p_A = inp.ParamsNeg[3 * m_D];
            //double p_B = inp.ParamsPos[3 * m_D];

            double ret = 0.0;

            for (int d = 0; d < m_D; d++) {
                //ret += 0.5 * (p_A * Vel_A[d] + p_B * Vel_B[d]) * inp.n[d];  // pressure
                for (int dd = 0; dd < m_D; dd++) {
                    ret -= 0.5 * (muA * GradVel_A[d, dd] * Vel_A[dd] + muB * GradVel_B[d, dd] * Vel_B[dd]) * inp.Normal[d];  // gradU
                    if(transposedTerm)
                        ret -= 0.5 * (muA * GradVel_A[dd, d] * Vel_A[dd] + muB * GradVel_B[dd, d] * Vel_B[dd]) * inp.Normal[d];  // gradU transposed
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

        public string PositiveSpecies {
            get { return "B"; }
        }

        public string NegativeSpecies {
            get { return "A"; }
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



    public class StressDivergence_Local : IVolumeForm, ISpeciesFilter {


        /// <summary>
        /// Spatial dimension;
        /// </summary>
        int m_D;

        double mu;

        bool transposedTerm;


        public StressDivergence_Local(int SpatDim, double _mu, string spcNmn, SpeciesId spcId, bool transposed = false) {
            m_D = SpatDim;
            mu = _mu;
            m_spcId = spcId;
            transposedTerm = transposed;
            ValidSpecies = spcNmn;
        }


        SpeciesId m_spcId;

        public string ValidSpecies {
            get;
            private set;
        }

        public IList<String> ArgumentOrdering {
            get {
                return new string[] { };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.VelocityX_GradientVector(), VariableNames.VelocityY_GradientVector(), 
                    new string[] { "VelocityXGradX_GradientX", "VelocityXGradX_GradientY" },
                    new string[] { "VelocityXGradY_GradientX", "VelocityXGradY_GradientY" },
                    new string[] { "VelocityYGradX_GradientX", "VelocityYGradX_GradientY" },
                    new string[] { "VelocityYGradY_GradientX", "VelocityYGradY_GradientY" });
            }
        }


        static double[,] VelocityGradient(double[] GradVelX, double[] GradVelY) {
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


        //static double[] LaplaceU(double uxx, double uyy, double vxx, double vyy) {

        //    double[] LapU = new double[2];

        //    LapU[0] = uxx + uyy;
        //    LapU[1] = vxx + vyy;

        //    return LapU;
        //}
 

        public double VolumeForm(ref CommonParamsVol cpv, Double[] U, Double[,] GradU, Double V, Double[] GradV) {

            double[] Vel = cpv.Parameters.GetSubVector(0, m_D);
            double[,] GradVel = VelocityGradient(cpv.Parameters.GetSubVector(m_D, m_D), cpv.Parameters.GetSubVector(2 * m_D, m_D));

            double[] LapU = new double[2];
            LapU[0] = cpv.Parameters[3 * m_D] + cpv.Parameters[4 * m_D + 1];
            LapU[1] = cpv.Parameters[5 * m_D] + cpv.Parameters[6 * m_D + 1];

            double[] DivGradUT = new double[2];
            DivGradUT[0] = cpv.Parameters[3 * m_D] + cpv.Parameters[5 * m_D + 1];
            DivGradUT[1] = cpv.Parameters[4 * m_D] + cpv.Parameters[6 * m_D + 1];

            double ret = 0;

            for (int d = 0; d < m_D; d++) {
                ret += DivGradUT[d] * Vel[d];
                if (transposedTerm)
                    ret += LapU[d] * Vel[d];
                for (int dd = 0; dd < m_D; dd++) {
                    ret += GradVel[dd, d] * GradVel[d, dd];
                    if(transposedTerm)
                            ret += GradVel[d, dd] * GradVel[d, dd];
                    }
            }

            return -mu * ret * V;
        }


        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.V;
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

        bool withPressure;

        public Dissipation(int SpatDim, double _mu, string spcNmn, SpeciesId spcId, bool _withPressure) {
            m_D = SpatDim;
            mu = _mu;
            m_spcId = spcId;
            ValidSpecies = spcNmn;
            this.withPressure = _withPressure;
        }


        SpeciesId m_spcId;

        public string ValidSpecies {
            get;
            private set;
        }

        public IList<String> ArgumentOrdering {
            get {
                return new string[] { };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.VelocityVector(m_D), 
                    VariableNames.VelocityX_GradientVector(), VariableNames.VelocityY_GradientVector(), 
                    VariableNames.Pressure);
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


        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            double[] Vel = cpv.Parameters.GetSubVector(0, m_D);
            double[,] GradVel = VelociytGradient(cpv.Parameters.GetSubVector(m_D, m_D), cpv.Parameters.GetSubVector(2 * m_D, m_D));

            double ret = (withPressure) ? -cpv.Parameters[3 * m_D] * (GradVel[0, 0] + GradVel[1, 1]) : 0.0;

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
