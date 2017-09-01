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

namespace NSE_SIMPLE {

    /// <summary>
    /// Collection of matrix assemblies for incompressible flows.
    /// </summary>
    public class MatrixFactoryIncompressibleFlows {

        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="SolverConf"></param>
        /// <param name="IncompressibleOperators"></param>        
        /// <param name="PressureMapping"></param>
        /// <param name="Pressure"></param>
        /// <param name="BDF"></param>
        public MatrixFactoryIncompressibleFlows(SolverConfiguration SolverConf, OperatorFactoryFlowFieldIncompressible IncompressibleOperators,
            UnsetteledCoordinateMapping PressureMapping, SinglePhaseField Pressure, BDFScheme BDF) {

            // Initialize Predictor
            Predictor = new SIMPLEMatrixAssembly[SolverConf.SpatialDimension];
            for (int d = 0; d < SolverConf.SpatialDimension; d++)
                Predictor[d] = new MatrixAssemblyIncompressiblePredictor(IncompressibleOperators.Convection[d], IncompressibleOperators.Visc[d], 2, 1);

            PredictorApprox = new MatrixAssemblyApprox(
                SolverConf, PressureMapping.GridDat.iLogicalCells.NoOfLocalUpdatedCells, Predictor[0], BDF, 1 + (1 + SolverConf.SpatialDimension) * SolverConf.Control.PredictorApproximationUpdateCycle);

            // Initialize Corrector
            PredictorApproxInv = new MatrixAssemblyApproxInv(
                SolverConf.Control, PressureMapping.GridDat.iLogicalCells.NoOfLocalUpdatedCells, PredictorApprox, SolverConf.SpatialDimension + SolverConf.SpatialDimension * SolverConf.Control.PredictorApproximationUpdateCycle);

            switch (SolverConf.Control.PredictorApproximation) {
                case PredictorApproximations.Identity:
                case PredictorApproximations.Diagonal:
                case PredictorApproximations.BlockDiagonal:
                    Corrector = new MatrixAssemblyIncompressibleCorrector(IncompressibleOperators.VelocityDivergence,
                            PredictorApproxInv,
                            IncompressibleOperators.PressureGradient,
                            SolverConf.Control.PredictorApproximationUpdateCycle,
                            IncompressibleOperators.PressureStabilization);                                          
                    break;
                case PredictorApproximations.Identity_IP1:
                    Corrector = new MatrixAssemblyCorrectorIP1(IncompressibleOperators.IP1PressureCorrection,
                        IncompressibleOperators.PressureStabilization,
                        SolverConf,
                        BDF);
                    break;                
                default:
                    throw new NotSupportedException("Unknown option for extended property 'Option_Approximation_Predictor'.");
            }
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
        public SIMPLEMatrixAssembly PredictorApprox {
            get;
            private set;
        }

        /// <summary>
        /// Inverse of approximation of Predictor.
        /// </summary>
        public SIMPLEMatrixAssembly PredictorApproxInv {
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
    }
}
