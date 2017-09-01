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
using System.Linq;
using System.Text;
using BoSSS.Solution.Utils;
using ilPSP.Utils;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.NSECommon {

    /// <summary>
    /// Numerical flux for advection of scalar quantity
    /// using the coupled Lax-Friedrichs flux for velocity and scalar,
    /// cf. F. Pochet, K. Hillewaert, P. Geuzaine, J.-F. Remacle, and È. Marchandise,
    /// “A 3D strongly coupled implicit discontinuous Galerkin level set-based method for modeling two-phase flows,”
    /// Comput. Fluids, 2013.
    /// </summary>
    public class CoupledLaxFriedrichsScalar : LinearFlux {

        int SpatDimension;
        MaterialLaw EoS;

        BoundaryCondMap<IncompressibleBcType> bcmap;
        Func<double[], double, double>[,] velFunction;
        Func<double[], double, double>[] scalarFunction;

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="SpatDimension"></param>
        /// <param name="EoS"></param>
        /// <param name="bcmap"></param>
        public CoupledLaxFriedrichsScalar(int SpatDimension, MaterialLaw EoS, IncompressibleBoundaryCondMap bcmap) {
            this.SpatDimension = SpatDimension;
            this.EoS = EoS;
            this.bcmap = bcmap;

            velFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, SpatDimension];
            for (int d = 0; d < SpatDimension; d++)
                velFunction.SetColumn(bcmap.bndFunction[VariableNames.Velocity_d(d)], d);

            scalarFunction = bcmap.bndFunction[VariableNames.LevelSet];
        }

        protected override double BorderEdgeFlux(ref Foundation.CommonParamsBnd inp, double[] Uin) {

            IncompressibleBcType edgeType = bcmap.EdgeTag2Type[inp.EdgeTag];

            switch (edgeType) {
                case IncompressibleBcType.Wall:
                case IncompressibleBcType.Velocity_Inlet: {

                        Foundation.CommonParams inp2;
                        inp2.GridDat = inp.GridDat;
                        inp2.Normale = inp.Normale;
                        inp2.iEdge = inp.iEdge;
                        inp2.Parameters_IN = inp.Parameters_IN;
                        inp2.X = inp.X;
                        inp2.time = inp.time;

                        // Boundary conditions scalar
                        // ==========================

                        inp2.Parameters_OUT = new double[2 * SpatDimension + 2];

                        double ScalarOut = scalarFunction[inp.EdgeTag](inp.X, 0.0);
                        inp2.Parameters_OUT[2 * SpatDimension] = scalarFunction[inp.EdgeTag](inp.X, inp.time);
                        // ScalarMeanOut = ScalarMeanIn
                        inp2.Parameters_OUT[2 * SpatDimension + 1] = inp.Parameters_IN[2 * SpatDimension + 1];

                        // Boundary conditions velocity
                        // ============================

                        for (int j = 0; j < SpatDimension; j++) {
                            inp2.Parameters_OUT[j] = velFunction[inp.EdgeTag, j](inp.X, inp.time);
                            // VelocityMeanOut = VelocityMeanIn
                            inp2.Parameters_OUT[SpatDimension + j] = inp.Parameters_IN[SpatDimension + j];
                        }

                        return InnerEdgeFlux(ref inp2, Uin, new double[] { ScalarOut });
                    }
                case IncompressibleBcType.Pressure_Outlet: {
                        double r = 0.0;

                        double Scalar0 = inp.Parameters_IN[2 * SpatDimension];
                        double rho = EoS.GetDensity(Scalar0);
                        double Scalar = Uin[0];

                        for (int j = 0; j < SpatDimension; j++) {
                            double u_j = inp.Parameters_IN[j];
                            r += rho * Scalar * u_j * inp.Normale[j];
                        }

                        return r;
                    }
                default:
                    throw new NotImplementedException("Boundary condition not implemented!");
            }
        }

        protected override double InnerEdgeFlux(ref Foundation.CommonParams inp, double[] Uin, double[] Uout) {
            double res = 0.0;

            // Calculate Lambda
            // ================

            IList<double> _VelocityMeanIn = new List<double>();
            IList<double> _VelocityMeanOut = new List<double>();
            for (int d = 0; d < SpatDimension; d++) {
                _VelocityMeanIn.Add(inp.Parameters_IN[SpatDimension + d]);
                _VelocityMeanOut.Add(inp.Parameters_OUT[SpatDimension + d]);
            }
            double[] VelocityMeanIn = _VelocityMeanIn.ToArray();
            double[] VelocityMeanOut = _VelocityMeanOut.ToArray();

            double ScalarMeanIn = inp.Parameters_IN[2 * SpatDimension + 1];
            double ScalarMeanOut = inp.Parameters_OUT[2 * SpatDimension + 1];

            double LambdaIn = EoS.GetLambda(VelocityMeanIn, inp.Normale, ScalarMeanIn);
            double LambdaOut = EoS.GetLambda(VelocityMeanOut, inp.Normale, ScalarMeanOut);

            double Lambda = Math.Max(LambdaIn, LambdaOut);

            // Calculate central part
            // ======================

            double Scalar0In = inp.Parameters_IN[2 * SpatDimension];
            double Scalar0Out = inp.Parameters_OUT[2 * SpatDimension];

            double rhoIn = EoS.GetDensity(Scalar0In);
            double rhoOut = EoS.GetDensity(Scalar0Out);

            double ScalarIn = Uin[0];
            double ScalarOut = Uout[0];

            for (int j = 0; j < SpatDimension; j++) {
                double u_j_In = inp.Parameters_IN[j];
                double u_j_Out = inp.Parameters_OUT[j];
                res += 0.5 * rhoIn * ScalarIn * u_j_In * inp.Normale[j];
                res += 0.5 * rhoOut * ScalarOut * u_j_Out * inp.Normale[j];
            }

            // Calculate dissipative part
            // ==========================

            // Jump scalar
            double Jump_Scalar = ScalarIn - ScalarOut;
            double TransformationFactorScalar = GetTransformationFactorScalar(ScalarMeanIn, ScalarMeanOut);
            res += 0.5 * Lambda * TransformationFactorScalar * Jump_Scalar;

            return res;
        }

        private double GetTransformationFactorScalar(double ScalarMeanIn, double ScalarMeanOut) {
            double Scalar_Average = 0.5 * (ScalarMeanIn + ScalarMeanOut);
            double Rho_Average = EoS.GetDensity(Scalar_Average);
            double DrhoDT_Average = EoS.DiffRho_Temp(Scalar_Average);
            return (Rho_Average + Scalar_Average * DrhoDT_Average);
        }

        protected override void Flux(ref Foundation.CommonParamsVol inp, double[] U, double[] output) {

            double rho = EoS.GetDensity(inp.Parameters[2 * SpatDimension]);
            double Temp = U[0];
            double u_1 = inp.Parameters[0];
            double u_2 = inp.Parameters[1];

            output[0] = rho * Temp * u_1;
            output[1] = rho * Temp * u_2;

            if (SpatDimension == 3) {
                double u_3 = inp.Parameters[2];
                output[2] = rho * Temp * u_3;
            }
        }

        public override IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.LevelSet };
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(SpatDimension),
                    VariableNames.Velocity0MeanVector(SpatDimension),
                    VariableNames.Phi0,
                    VariableNames.Phi0Mean);
            }
        }
    }
}
