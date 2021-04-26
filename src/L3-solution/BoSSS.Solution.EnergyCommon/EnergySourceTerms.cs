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


        public PowerofGravity(int SpatDim, string spcNmn, SpeciesId spcId, double _rho) {
            m_D = SpatDim;
            m_spcId = spcId;
            ValidSpecies = spcNmn;
            rho = _rho;
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
                return ArrayTools.Cat(VariableNames.VelocityVector(m_D), VariableNames.GravityVector(m_D));
            }
        }


        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {

            double[] Vel = cpv.Parameters.GetSubVector(0, m_D);
            double[] Grav = cpv.Parameters.GetSubVector(m_D, m_D);

            double ret = 0;

            for (int d = 0; d < m_D; d++) {
                ret -= Grav[d] * Vel[d];
            }

            return rho * ret * V;
        }


        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.V;
            }
        }


    }



    public class SurfaceEnergy : ILevelSetForm {

        //LevelSetTracker m_LsTrk;

        public SurfaceEnergy(int _D, double _sigma, double _rhoA, double _rhoB) {
            //m_LsTrk = LsTrk;
            this.m_D = _D;
            this.sigma = _sigma;
            this.rhoA = _rhoA;
            this.rhoB = _rhoB;
        }

        int m_D;

        double sigma;
        double rhoA;
        double rhoB;


        private double[] GetInterfaceValue(double[] ValA, double[] ValB) {

            double[] ValI = new double[m_D];
            for (int d = 0; d < m_D; d++) {
                ValI[d] = (rhoA * ValA[d] + rhoB * ValB[d]) / (rhoA + rhoB);
            }

            return ValI;
        }

        //protected static double[] SurfaceNormal(double[] param) {

        //    double[] N = new double[param.Length];

        //    for (int d = 0; d < param.Length; d++) {
        //        N[d] = param[d];
        //    }

        //    return N.Normalize();
        //}

        protected static double[,] SurfaceProjection(double[] Nsurf) {

            int D = Nsurf.Length;
            double[,] P = new double[D, D];

            for (int d = 0; d < D; d++) {
                for (int dd = 0; dd < D; dd++) {
                    if (dd == d)
                        P[d, dd] = (1.0 - Nsurf[d] * Nsurf[dd]);
                    else
                        P[d, dd] = (0.0 - Nsurf[d] * Nsurf[dd]);
                }
            }

            return P;
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


        public double InnerEdgeForm(ref CommonParams cp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            //Debug.Assert(cp.ParamsPos[0] == cp.ParamsNeg[0], "interface velocityX must be continuous across interface");
            //Debug.Assert(cp.ParamsPos[1] == cp.ParamsNeg[1], "interface velocityY must be continuous across interface");
            Debug.Assert(cp.Parameters_OUT[3 * m_D] == cp.Parameters_IN[3 * m_D], "curvature must be continuous across interface");

            double curvature = cp.Parameters_IN[3 * m_D];
            double[] VelI = GetInterfaceValue(cp.Parameters_IN.GetSubVector(0, m_D), cp.Parameters_OUT.GetSubVector(0, m_D));
            double[] Normal = cp.Normal;
            double[,] Psurf = SurfaceProjection(Normal);
            double[] GradVelXI = GetInterfaceValue(cp.Parameters_IN.GetSubVector(m_D, m_D), cp.Parameters_OUT.GetSubVector(m_D, m_D));
            double[] GradVelYI = GetInterfaceValue(cp.Parameters_IN.GetSubVector(2 * m_D, m_D), cp.Parameters_OUT.GetSubVector(2 * m_D, m_D));
            double[,] GradVelI = VelocityGradient(GradVelXI, GradVelYI);


            double surfE = 0;
            for (int d = 0; d < m_D; d++) {
                surfE -= curvature * (VelI[d] * Normal[d]);
                //for (int dd = 0; dd < m_D; dd++) {
                //    surfE -= Psurf[d, dd] * GradVelI[dd, d];
                //}
            }
            surfE *= sigma;

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
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.VelocityX_GradientVector(), VariableNames.VelocityY_GradientVector(), VariableNames.Curvature);
            }
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public string PositiveSpecies {
            get { return "B"; }
        }

        public string NegativeSpecies {
            get { return "A"; }
        }

        public TermActivationFlags LevelSetTerms {
            get { return TermActivationFlags.V; }
        }


    }



}
