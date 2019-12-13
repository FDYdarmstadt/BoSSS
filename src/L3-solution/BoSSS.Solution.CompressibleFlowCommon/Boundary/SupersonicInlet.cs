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
using BoSSS.Platform.LinAlg;
using BoSSS.Solution.CompressibleFlowCommon;
using System.Diagnostics;
using ilPSP;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;

namespace BoSSS.Solution.CompressibleFlowCommon.Boundary {

    /// <summary>
    /// Implementation of the boundary condition for a supersonic inlet (i.e.
    /// Ma > 1.0). In this case, all values at the boundary have to be
    /// prescribed. For this matter, we choose the density, the momentum and
    /// the pressure as input.
    /// </summary>
    public class SupersonicInlet : BoundaryCondition {

        /// <summary>
        /// The prescribed density at the inlet
        /// </summary>
        public readonly Func<double[], double, double> DensityFunction;

        /// <summary>
        /// The prescribed velocity components at the inlet
        /// </summary>
        public readonly Func<double[], double, double>[] VelocityFunctions;

        /// <summary>
        /// The prescribed pressure at the inlet
        /// </summary>
        public Func<double[], double, double> PressureFunction;

        /// <summary>
        /// Constructs a new supersonic inlet using the values defined by
        /// <paramref name="densityFunction"/>,
        /// <paramref name="velocityFunctions"/> and
        /// <paramref name="pressureFunction"/> as values at the boundary.
        /// </summary>
        /// <param name="config"><see cref="BoundaryCondition"/></param>
        /// <param name="densityFunction">
        /// The prescribed density at the inlet
        /// </param>
        /// <param name="velocityFunctions">
        /// The prescribed velocity components at the inlet
        /// </param>
        /// <param name="pressureFunction">
        /// The prescribed pressure at the inlet
        /// </param>
        public SupersonicInlet(MaterialProperty.Material config, Func<double[], double, double> densityFunction, Func<double[], double, double>[] velocityFunctions, Func<double[], double, double> pressureFunction)
            : base(config) {
            this.DensityFunction = densityFunction;
            this.VelocityFunctions = velocityFunctions;
            this.PressureFunction = pressureFunction;
        }

        /// <summary>
        /// Calculates the boundary values for a supersonic inlet (Mach number
        /// greater than 1). Theoretically, we have to impose five conditions
        /// (in the inviscid as well as in the viscid case!) which is why
        /// calculate the complete state from the given density, momentum and
        /// pressure
        /// </summary>
        /// <returns>
        /// \f$ (\rho^*, u_0^*[, u_1^*[, u_2^*]], p*)^T\f$ 
        /// </returns>
        public override StateVector GetBoundaryState(double time, Vector x, Vector normal, StateVector stateIn) {
            //Convection.OptimizedHLLCFlux.SupersonicInlet.Start();

            Debug.Assert(x.Dim == stateIn.Dimension);
            int D = x.Dim;

            Vector uOut = new Vector(D);
            for (int i = 0; i < D; i++) {
                uOut[i] = VelocityFunctions[i](x, time);
            }

            StateVector retval = StateVector.FromPrimitiveQuantities(
                stateIn.Material,
                DensityFunction(x, time),
                uOut,
                PressureFunction(x, time));
            //Convection.OptimizedHLLCFlux.SupersonicInlet.Stop();
            return retval;
        }

        
      
        /// <summary>
        /// Vectorized implementation of <see cref="GetBoundaryState(double, Vector, Vector, StateVector)"/>
        /// </summary>
        public override void GetBoundaryState(MultidimensionalArray[] StateOut, double time, MultidimensionalArray X, MultidimensionalArray Normals, MultidimensionalArray[] StateIn, int Offset, int NoOfEdges, bool normalFlipped, Material material) {
            //Convection.OptimizedHLLCFlux.SupersonicInlet.Start();
            if (X.Dimension != 3)
                throw new ArgumentException();
            int D = X.GetLength(2);
            int NoOfNodes = X.GetLength(1);
            
            if (StateIn.Length != D + 2)
                throw new ArgumentException();
            if (StateOut.Length != D + 2)
                throw new ArgumentException();
            bool is3D = D >= 3;
            if (D < 2) {
                // base-implementation supports also 1D;
                // The implementation here is (a bit) tuned for performance, we only support 2D and 3D;
                // in 1D, performance is not so relevant anyway.
                base.GetBoundaryState(StateOut, time, X, Normals, StateIn, Offset, NoOfEdges, normalFlipped, material);
                return;
            }

            var Density = StateOut[0];
            var Energy = StateOut[D + 1];
            var MomentumX = StateOut[1];
            var MomentumY = StateOut[2];
            var MomentumZ = is3D ? StateOut[3] : null;

            var DensityIn = StateIn[0];
            //var EnergyIn = StateIn[D + 1];
            //var MomentumXin = StateIn[1];
            //var MomentumYin = StateIn[2];
            //var MomentumZin = is3D ? StateIn[3] : null;

            double MachScaling = material.EquationOfState.HeatCapacityRatio * material.MachNumber * material.MachNumber;

            double[] xLocal = new double[D];
            //Vector uOut = new Vector(D);
            for (int e = 0; e < NoOfEdges; e++) {
                int edge = e + Offset;

                // Loop over nodes
                for (int n = 0; n < NoOfNodes; n++) {
                    xLocal[0] = X[edge, n, 0];
                    xLocal[1] = X[edge, n, 1];
                    double uOut_x = VelocityFunctions[0](xLocal, time);
                    double uOut_y = VelocityFunctions[1](xLocal, time);
                    double uOut_z = 0.0;
                    double velocity_AbsSquare = uOut_x * uOut_x + uOut_y * uOut_y;
                    if(is3D) {
                        xLocal[2] = X[edge, n, 2];
                        uOut_z = VelocityFunctions[2](xLocal, time);
                        velocity_AbsSquare += uOut_z * uOut_z;
                    }

#if DEBUG
                    var uOutVec = new Vector(D);
                    uOutVec.x = uOut_x;
                    uOutVec.y = uOut_y;
                    uOutVec.z = uOut_z;
                    StateVector stateBoundary = StateVector.FromPrimitiveQuantities(
                        material,
                        DensityFunction(xLocal, time),
                        uOutVec,
                        PressureFunction(xLocal, time));
#endif
                    double density = DensityFunction(xLocal, time);
                    double pressure = PressureFunction(xLocal, time);

                    Density[edge, n] = density;
                    MomentumX[edge, n] = uOut_x*density;
                    MomentumY[edge, n] = uOut_y*density;
                    if(is3D)
                        MomentumZ[edge, n] = uOut_z*density;
                    Energy[edge, n] = material.EquationOfState.GetInnerEnergy(density, pressure) + 0.5 * MachScaling * density * velocity_AbsSquare;
#if DEBUG
                    Debug.Assert(Density[edge, n].RelErrorSmallerEps(stateBoundary.Density), "density error");
                    for(int d = 0; d < D; d++)
                        Debug.Assert(StateOut[d + 1][edge, n].RelErrorSmallerEps(stateBoundary.Momentum[d]), "momentum error");
                    Debug.Assert(Energy[edge, n].RelErrorSmallerEps(stateBoundary.Energy), "energy error");

#endif
                }
            }

            //Convection.OptimizedHLLCFlux.SupersonicInlet.Stop();


        }
    }
}
