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
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Utils;
using ilPSP.Utils;
using BoSSS.Solution.LevelSetTools.Smoothing;
using BoSSS.Solution.Timestepping;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using System.Diagnostics;
using ilPSP.Tracing;
using ilPSP;
using BoSSS.Platform;
using BoSSS.Solution.NSECommon;
using ilPSP.LinSolvers;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.LevelSetTools.Advection {
    /// <summary>
    /// A abstract class for calculating the numerical flux included in the spatial discretization of the level set advection equation
    /// The edge Fluxes must be implemented
    /// </summary>
    abstract class LevelSetAdvectionFlux : LinearFlux {
        public LevelSetAdvectionFlux(GridData GridDat, IncompressibleBoundaryCondMap BcMap) {
            this.GridDat = GridDat;
            D = this.GridDat.SpatialDimension;
            this.BcMap = BcMap;

            VelFunction = new Func<double[], double, double>[GridCommons.FIRST_PERIODIC_BC_TAG, D];
            for (int d = 0; d < D; d++)
                VelFunction.SetColumn(this.BcMap.bndFunction[VariableNames.Velocity_d(d)], d);

            LevelSetFunction = this.BcMap.bndFunction[VariableNames.LevelSet];
        }

        /// <summary>
        /// Spatial dimension
        /// </summary>
        GridData GridDat;
        protected int D;
        protected IncompressibleBoundaryCondMap BcMap;
        Func<double[], double, double>[,] VelFunction;
        protected Func<double[], double, double>[] LevelSetFunction;

        abstract protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin);


        abstract protected override double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout);

        protected override void Flux(ref CommonParamsVol inp, double[] U, double[] output) {
            for (int d = D - 1; d >= 0; d--)
                output[d] = inp.Parameters[d] * U[0];
        }

        public override IList<string> ArgumentOrdering
        {
            get
            {
                return new string[] { "LevelSet" };
            }
        }

        public override IList<string> ParameterOrdering
        {
            get
            {
                //var ParamList = D.ForLoop(d => string.Format("U[{0}]", d));
                //ParamList = ParamList.Cat(D.ForLoop(d => string.Format("U_Mean[{0}]", d)));
                //return ParamList;
                return ArrayTools.Cat(VariableNames.Velocity0Vector(D), VariableNames.Velocity0MeanVector(D));
            }
        }
    }

    class LevelSetUpwindFlux : LevelSetAdvectionFlux {
        public LevelSetUpwindFlux(GridData GridDat, IncompressibleBoundaryCondMap BcMap) : base(GridDat, BcMap) {
        }

        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {
            //double flux;

            //double U_dot_Normal = 0;
            //for (int d = D - 1; d >= 0; d--) {
            //    U_dot_Normal += (inp.Parameters_IN[d] * inp.Normale[d]);
            //}
            //flux = Uin[0] * U_dot_Normal;

            //return flux;



            double flux;
            var normal = inp.Normale;
            var EdgeTag = inp.EdgeTag;
            var x = inp.X;

            IncompressibleBcType edgeType = BcMap.EdgeTag2Type[EdgeTag];

            double U_dot_Normal = 0;
            //for (int d = D - 1; d >= 0; d--) {
            //    U_dot_Normal += (inp.Parameters_IN[d] * normal[d]);
            //}

            //if (U_dot_Normal >= 0) {
            //    // flux is 'going out of this cell'
            //    flux = Uin[0] * U_dot_Normal;
            //}
            //else {
            //    // flux is 'coming into this cell'
            //    flux = Uin[0] * U_dot_Normal; // original
            //}

            switch (edgeType) {
                case IncompressibleBcType.Wall:
                    flux = 0.0;
                    break;
                case IncompressibleBcType.Velocity_Inlet:

                    double Umean_dot_Normal = 0;
                    for (int d = D - 1; d >= 0; d--) {
                        Umean_dot_Normal += inp.Parameters_IN[D + d] * normal[d];
                    }

                    flux = 0;
                    if (Umean_dot_Normal >= 0) {
                        for (int d = D - 1; d >= 0; d--) {
                            flux += inp.Parameters_IN[d] * normal[d];
                        }
                        flux *= Uin[0];
                    }
                    else {
                        for (int d = D - 1; d >= 0; d--) {
                            flux += inp.Parameters_IN[d] * normal[d];
                        }
                        flux *= LevelSetFunction[EdgeTag](x, 0);
                    }
                    //for (int d = D - 1; d >= 0; d--) {
                    //    U_dot_Normal += (VelFunction[EdgeTag, d](x) * normal[d]);
                    //}
                    //flux = LevelSetFunction[EdgeTag](x) * U_dot_Normal;
                    break;
                case IncompressibleBcType.Outflow:
                    for (int d = D - 1; d >= 0; d--) {
                        U_dot_Normal += (inp.Parameters_IN[d] * normal[d]);
                    }
                    flux = Uin[0] * U_dot_Normal;
                    break;
                case IncompressibleBcType.Pressure_Outlet:
                    // In some bad cases this might cause similar stability problems as the inflow without a BC
                    // Probably Dirichlet BC's are more suitable here
                    for (int d = D - 1; d >= 0; d--) {
                        U_dot_Normal += (inp.Parameters_IN[d] * normal[d]);
                    }
                    flux = Uin[0] * U_dot_Normal;
                    break;
                default:
                    throw new NotImplementedException("boundary condition not implemented!");
            }


            return flux;

        }

        protected override double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout) {

            double Umean_dot_Normal = 0;
            for (int d = D - 1; d >= 0; d--) {
                Umean_dot_Normal += (inp.Parameters_IN[d] + inp.Parameters_OUT[d]) * inp.Normale[d];
            }
            Umean_dot_Normal *= 0.5;

            double flux = 0;
            if (Umean_dot_Normal >= 0) {
                for (int d = D - 1; d >= 0; d--) {
                    flux += inp.Parameters_IN[d] * inp.Normale[d];
                }
                flux *= Uin[0];
            }
            else {
                for (int d = D - 1; d >= 0; d--) {
                    flux += inp.Parameters_OUT[d] * inp.Normale[d];
                }
                flux *= Uout[0];
            }
            return flux;
        }
    }

    class LevelSetLLFFlux : LevelSetAdvectionFlux {
        public LevelSetLLFFlux(GridData GridDat, IncompressibleBoundaryCondMap BcMap) : base(GridDat, BcMap) {
        }
        protected override double BorderEdgeFlux(ref CommonParamsBnd inp, double[] Uin) {
            //double flux;

            //double U_dot_Normal = 0;
            //for (int d = D - 1; d >= 0; d--) {
            //    U_dot_Normal += (inp.Parameters_IN[d] * inp.Normale[d]);
            //}
            //flux = Uin[0] * U_dot_Normal;

            //return flux;



            double flux;
            var normal = inp.Normale;
            var EdgeTag = inp.EdgeTag;
            var x = inp.X;

            //IncompressibleBcType edgeType = BcMap.EdgeTag2Type[EdgeTag];

            double U_dot_Normal = 0;
            //for (int d = D - 1; d >= 0; d--) {
            //    U_dot_Normal += (inp.Parameters_IN[d] * normal[d]);
            //}

            //if (U_dot_Normal >= 0) {
            //    // flux is 'going out of this cell'
            //    flux = Uin[0] * U_dot_Normal;
            //}
            //else {
            //    // flux is 'coming into this cell'
            //    flux = Uin[0] * U_dot_Normal; // original
            //}


            // In some bad cases this might cause similar stability problems as the inflow without a BC
            // Probably Dirichlet BC's are more suitable here
            for (int d = D - 1; d >= 0; d--) {
                U_dot_Normal += (inp.Parameters_IN[d] * normal[d]);
            }
            flux = Uin[0] * U_dot_Normal;

            //switch (edgeType) {
            //    case IncompressibleBcType.FreeSlip:
            //    case IncompressibleBcType.Wall:
            //    case IncompressibleBcType.Velocity_Inlet:
            //    case IncompressibleBcType.Outflow:
            //    case IncompressibleBcType.Pressure_Outlet:
            //        // In some bad cases this might cause similar stability problems as the inflow without a BC
            //        // Probably Dirichlet BC's are more suitable here
            //        for (int d = D - 1; d >= 0; d--) {
            //            U_dot_Normal += (inp.Parameters_IN[d] * normal[d]);
            //        }
            //        flux = Uin[0] * U_dot_Normal;
            //        break;

            //    default:
            //        throw new NotImplementedException("boundary condition not implemented!");
            //}


            return flux;

        }

        protected override double InnerEdgeFlux(ref CommonParams inp, double[] Uin, double[] Uout) {
            double flux = 0;

            //Central Part
            flux += Uin[0] * (inp.Parameters_IN[0] * inp.Normale[0] + inp.Parameters_IN[1] * inp.Normale[1]);
            flux += Uout[0] * (inp.Parameters_OUT[0] * inp.Normale[0] + inp.Parameters_OUT[1] * inp.Normale[1]);
            if (D == 3) {
                flux += Uin[0] * inp.Parameters_IN[2] * inp.Normale[2] + Uout[0] * inp.Parameters_OUT[2] * inp.Normale[2];
            }

            //Dissipative Part
            double[] VelocityMeanIn = new double[D];
            double[] VelocityMeanOut = new double[D];

            for (int d = 0; d < D; d++) {
                VelocityMeanIn[d] = inp.Parameters_IN[D + d];
                VelocityMeanOut[d] = inp.Parameters_OUT[D + d];
            }


            double LambdaIn;
            double LambdaOut;
            LambdaIn = LambdaConvection.GetLambda(VelocityMeanIn, inp.Normale, false);
            LambdaOut = LambdaConvection.GetLambda(VelocityMeanOut, inp.Normale, false);

            double Lambda = Math.Max(LambdaIn, LambdaOut);
            double Scalar_Jump = Uin[0] - Uout[0];

            flux += Lambda * Scalar_Jump;

            flux *= 0.5;
            return flux;
        }
    }

    class FextSource : LinearSource {

        protected override double Source(double[] x, double[] parameters, double[] U) {
            //double divU = U[1];
            //double Phi = U[0];
            //return -divU * Phi;
            return -parameters[0] * U[0];
        }

        public override IList<string> ArgumentOrdering
        {
            get
            {
                return new string[] { "LevelSet" };
            }
        }

        public override IList<string> ParameterOrdering
        {
            get
            {
                return new string[] { "div(U)" };
            }
        }
    }
}
