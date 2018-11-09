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
using BoSSS.Platform;
using BoSSS.Foundation;
using BoSSS.Solution.NSECommon;
using ilPSP.LinSolvers;

namespace NSE_SIMPLE {

    /// <summary>
    /// Collection of matrix assemblies for variable density flows.
    /// </summary>
    public class MatrixFactoryVariableDensityFlowField {

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="SolverConf"></param>
        /// <param name="OperatorsFlowField"></param>
        /// <param name="Rho"></param>        
        /// <param name="BDF"></param>
        /// <param name="VelocityMapping"></param>
        /// <param name="VelocityVectorMapping"></param>
        public MatrixFactoryVariableDensityFlowField(SolverConfiguration SolverConf, OperatorFactoryFlowFieldVariableDensity OperatorsFlowField,
            BlockDiagonalMatrix Rho, BDFScheme BDF,
            UnsetteledCoordinateMapping VelocityMapping, UnsetteledCoordinateMapping VelocityVectorMapping) {

            // Initialize Predictor
            // ====================

            ViscSplit = new SIMPLEMatrixAssembly[SolverConf.SpatialDimension, SolverConf.SpatialDimension];
            for (int i = 0; i < SolverConf.SpatialDimension; i++) {
                for (int j = 0; j < SolverConf.SpatialDimension; j++) {
                    switch (SolverConf.Control.PhysicsMode) {
                        case PhysicsMode.LowMach:
                            ViscSplit[i, j] = new MatrixAssemblyViscSplit(OperatorsFlowField.Swip2[i], OperatorsFlowField.Swip3[i],
                                j, VelocityMapping, VelocityVectorMapping);
                            break;
                        case PhysicsMode.Multiphase:
                            ViscSplit[i, j] = new MatrixAssemblyViscSplit(OperatorsFlowField.Swip2[i], null,
                                j, VelocityMapping, VelocityVectorMapping);
                            break;
                        case PhysicsMode.Incompressible:
                            throw new ApplicationException("Using wrong matrix factory for incompressible flows");
                        default:
                            throw new NotImplementedException();
                    }
                }
            }

            //SaveMatricesToTextFile(OperatorsFlowField, VelocityMapping, VelocityVectorMapping);

            Predictor = new SIMPLEMatrixAssembly[SolverConf.SpatialDimension];
            PredictorApprox = new SIMPLEMatrixAssembly[SolverConf.SpatialDimension];
            PredictorApproxInv = new SIMPLEMatrixAssembly[SolverConf.SpatialDimension];
            for (int comp = 0; comp < SolverConf.SpatialDimension; comp++) {
                Predictor[comp] = new MatrixAssemblyVariableDensityPredictor(comp, OperatorsFlowField.Convection, OperatorsFlowField.Visc, ViscSplit);
                PredictorApprox[comp] = new MatrixAssemblyApprox(SolverConf, VelocityMapping.GridDat.iLogicalCells.NoOfLocalUpdatedCells, Predictor[comp], BDF, Rho, 1 + 2 * SolverConf.Control.PredictorApproximationUpdateCycle);
                PredictorApproxInv[comp] = new MatrixAssemblyApproxInv(SolverConf.Control, VelocityMapping.GridDat.iLogicalCells.NoOfLocalUpdatedCells, PredictorApprox[comp], 1 + SolverConf.Control.PredictorApproximationUpdateCycle);
            }

            // Initialize Corrector
            // ====================

            switch (SolverConf.Control.PhysicsMode) {
                case PhysicsMode.LowMach:
                    Corrector = new MatrixAssemblyVariableDensityCorrector(OperatorsFlowField.DivergenceConti, PredictorApproxInv, OperatorsFlowField.PressureGradient);
                    break;
                case PhysicsMode.Multiphase:
                    Corrector = new MatrixAssemblyVariableDensityCorrector(OperatorsFlowField.DivergenceConti, PredictorApproxInv, OperatorsFlowField.PressureGradient, SolverConf.Control.PredictorApproximationUpdateCycle);
                    break;
                case PhysicsMode.Incompressible:
                    throw new ApplicationException("Using wrong matrix factory for incompressible flows");
                default:
                    throw new NotImplementedException();
            }
        }

        /// <summary>
        /// Splitted matrices for viscosity term2 and term3.
        /// </summary>
        public SIMPLEMatrixAssembly[,] ViscSplit {
            get;
            private set;
        }

        /// <summary>
        /// Predictor.
        /// </summary>
        public SIMPLEMatrixAssembly[] Predictor {
            get;
            private set;
        }

        /// <summary>
        /// Approximation of Predictor.
        /// </summary>
        public SIMPLEMatrixAssembly[] PredictorApprox {
            get;
            private set;
        }

        /// <summary>
        /// Inverse of approximation of Predictor.
        /// </summary>
        public SIMPLEMatrixAssembly[] PredictorApproxInv {
            get;
            private set;
        }

        /// <summary>
        /// Corrector.
        /// </summary>
        public SIMPLEMatrixAssembly Corrector {
            get;
            private set;
        }

        /// <summary>
        /// Test code for debugging.
        /// </summary>
        /// <param name="OperatorsFlowField"></param>
        /// <param name="VelocityMapping"></param>
        /// <param name="VelocityVectorMapping"></param>
        private void SaveMatricesToTextFile(OperatorFactoryFlowFieldVariableDensity OperatorsFlowField,
            UnsetteledCoordinateMapping VelocityMapping, UnsetteledCoordinateMapping VelocityVectorMapping) {

            OperatorsFlowField.Swip2[0].OperatorMatrix.SaveToTextFileSparse("C:\\tmptest\\Swip20.txt");
            OperatorsFlowField.Swip2[1].OperatorMatrix.SaveToTextFileSparse("C:\\tmptest\\Swip21.txt");
            OperatorsFlowField.Swip3[0].OperatorMatrix.SaveToTextFileSparse("C:\\tmptest\\Swip30.txt");
            OperatorsFlowField.Swip3[1].OperatorMatrix.SaveToTextFileSparse("C:\\tmptest\\Swip31.txt");

            ViscSplit[0, 0].AssemblyMatrix.SaveToTextFileSparse("C:\\tmptest\\ViscSplit00.txt");
            ViscSplit[0, 1].AssemblyMatrix.SaveToTextFileSparse("C:\\tmptest\\ViscSplit01.txt");
            ViscSplit[1, 0].AssemblyMatrix.SaveToTextFileSparse("C:\\tmptest\\ViscSplit10.txt");
            ViscSplit[1, 1].AssemblyMatrix.SaveToTextFileSparse("C:\\tmptest\\ViscSplit11.txt");

            int[] IndicesVelocity = VelocityMapping.GetSubvectorIndices(true, 0);
            int[] IndicesVelocityVector0 = VelocityVectorMapping.GetSubvectorIndices(true, 0);
            int[] IndicesVelocityVector1 = VelocityVectorMapping.GetSubvectorIndices(true, 1);

            MsrMatrix Swip2Mtx = new MsrMatrix(VelocityVectorMapping);

            OperatorsFlowField.Swip2[0].OperatorMatrix.WriteSubMatrixTo<IList<int>, IList<int>, IList<int>, IList<int>>(Swip2Mtx,
                IndicesVelocity,
                IndicesVelocityVector0,
                IndicesVelocityVector0,
                IndicesVelocityVector0);
            OperatorsFlowField.Swip2[0].OperatorMatrix.AccSubMatrixTo<IList<int>, IList<int>, IList<int>, IList<int>>(1.0,
                Swip2Mtx,
                IndicesVelocity,
                IndicesVelocityVector0,
                IndicesVelocityVector1,
                IndicesVelocityVector1);
            OperatorsFlowField.Swip2[1].OperatorMatrix.AccSubMatrixTo<IList<int>, IList<int>, IList<int>, IList<int>>(1.0,
                Swip2Mtx,
                IndicesVelocity,
                IndicesVelocityVector1,
                IndicesVelocityVector0,
                IndicesVelocityVector0);
            OperatorsFlowField.Swip2[1].OperatorMatrix.AccSubMatrixTo<IList<int>, IList<int>, IList<int>, IList<int>>(1.0,
                Swip2Mtx,
                IndicesVelocity,
                IndicesVelocityVector1,
                IndicesVelocityVector1,
                IndicesVelocityVector1);

            Swip2Mtx.SaveToTextFileSparse("C:\\tmptest\\Swip2.txt");
        }
    }
}
