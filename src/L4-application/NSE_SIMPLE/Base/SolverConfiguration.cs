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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using NSE_SIMPLE.LowMach;
using System;

namespace NSE_SIMPLE {

    public class SolverConfiguration {

        public readonly SIMPLEControl Control;

        public SolverConfiguration(
            SIMPLEControl control, IGridData GridDat, VariableSet workingSet, int MPIRank) {

            this.SpatialDimension = GridDat.SpatialDimension;
            this.MPIRank = MPIRank;
            this.Control = control;
            this.BcMap = new IncompressibleBoundaryCondMap(GridDat, Control.BoundaryValues, Control.PhysicsMode);

            CheckSetupPressure(SpatialDimension);
            if (Control.PressureReferencePoint != null) {
                UnsetteledCoordinateMapping MappingPressure = new UnsetteledCoordinateMapping(
                    new Basis(GridDat, workingSet.Pressure.Basis.Degree));

                PressureReferencePointIndex = BoSSS.Solution.NSECommon.SolverUtils.GetIndexOfPressureReferencePoint(Control.PressureReferencePoint, MappingPressure, 0);
                if (PressureReferencePointIndex < 0) {
                    throw new ApplicationException("Unknown error finding index for pressure reference point.");
                }
            }

            // SIP penalty for viscous part of momentum equation
            this.PenaltyViscMomentum = Control.ViscousPenaltyScaling;

            // SIP penalty for SIP pressure correction
            if ((Control.PredictorApproximation == PredictorApproximations.Identity_IP1)) {
                this.PenaltyPressureCorrection = Control.PressureCorrectionPenaltyScaling * GetSIPPenaltyBase(workingSet.Pressure.Basis.Degree, GridDat);
            } else if (Control.PressureCorrectionPenaltyScaling != 1.0) {
                throw new NotSupportedException("Extended property 'penalty_PressureCorrection' is only available for Identity_IP option.");
            }

            if (control.PhysicsMode == PhysicsMode.LowMach) {
                LowMachSIMPLEControl lowMachControl = control as LowMachSIMPLEControl;
                this.PenaltyHeatConduction = lowMachControl.PenaltyHeatConductionScaling * GetSIPPenaltyBase(workingSet.Pressure.Basis.Degree, GridDat);
            }

            this.Velocity = workingSet.Velocity;

            dt = control.Endtime / control.NoOfTimesteps;
        }

        /// <summary>
        /// MPIRank of current process.
        /// </summary>
        public readonly int MPIRank;

        /// <summary>
        /// Spatial dimension of geometry.
        /// </summary>
        public readonly int SpatialDimension;

        /// <summary>
        /// Penalty factor SIP discretization of viscous terms.
        /// </summary>
        public readonly double PenaltyViscMomentum;

        /// <summary>
        /// Penalty factor SIP discretization of SIP pressure correction.
        /// </summary>
        public readonly double PenaltyPressureCorrection;

        /// <summary>
        /// Penalty for SIP discretization of heat conduction.
        /// </summary>
        public readonly double PenaltyHeatConduction;

        /// <summary>
        /// Boundary condition map.
        /// </summary>
        public readonly IncompressibleBoundaryCondMap BcMap;

        /// <summary>
        /// Global unique index for reference point of the pressure.
        /// </summary>
        public readonly int PressureReferencePointIndex = -1;

         
        /// <summary>
        /// Calculation of SIP penalty base, cf. Chapter 3 in 
        /// K. Hillewaert, “Development of the discontinuous Galerkin method for high-resolution, large scale CFD and acoustics in industrial geometries”,
        /// Université catholique de Louvain, 2013.
        /// </summary>
        /// <param name="DegreeBasis"></param>
        /// <param name="GridDat"></param>
        /// <returns>
        /// Base part of SIP penalty.
        /// </returns>
        private double GetSIPPenaltyBase(int DegreeBasis, IGridData GridDat) {
            if (!(GridDat is BoSSS.Foundation.Grid.Classic.GridData))
                throw new NotImplementedException("SIP penalty is only implemented for Classic grid.");

            if (GridDat.iGeomCells.RefElements.Length > 1)
                throw new NotImplementedException("SIP penalty implemented only for one type of RefElements.");

            Type GridType = GridDat.iGeomCells.RefElements[0].GetType();

            double PenaltyBase;

            if ((GridType == typeof(Triangle)) || (GridType == typeof(Tetra))) {
                PenaltyBase = (DegreeBasis + 1.0) * (DegreeBasis + GridDat.SpatialDimension) / GridDat.SpatialDimension;
            } else if ((GridType == typeof(Square)) || (GridType == typeof(Cube))) {
                PenaltyBase = (DegreeBasis + 1.0) * (DegreeBasis + 1.0);
            } else {
                throw new NotImplementedException("Unknown RefElement.");
            }

            return PenaltyBase;
        }
        

        /// <summary>
        /// Checks correct setup of pressure constraints / boundary conditions.
        /// </summary>
        /// <param name="D">Spatial Dimension</param>
        private void CheckSetupPressure(int D) {

            //check dimensions
            if ((Control.PressureReferencePoint != null) && (Control.PressureReferencePoint.Length != D))
                throw new ApplicationException("wrong dimension (number of elements) of extended property \"RefPtPressure\": grid is "
                            + D + "-dimensional, but vector length is " + Control.PressureReferencePoint.Length + ".");

            if ((Control.PressureGradientSource != null) && (Control.PressureGradientSource.Length != D))
                throw new ApplicationException("wrong dimension (number of elements) of extended property \"SourcePressureGradient\": grid is "
                            + D + "-dimensional, but vector length is " + Control.PressureGradientSource.Length + ".");

            //check pressure constraints
            if (this.BcMap.DirichletPressureBoundary) {
                if (Control.PressureReferencePoint != null)
                    throw new ApplicationException("There are too much constraints for the pressure. "
                        + "Only one 'pressure_outlet' or one reference point for the pressure ('RefPtPressure') "
                        + "can be given.");
            } else {
                if (Control.PressureReferencePoint == null)
                    throw new ApplicationException("Missing constraint for the pressure. "
                        + "Either one 'pressure_outlet' or one reference point for the pressure ('RefPtPressure') "
                        + "has to be specified.");
            }

            if ((Control.PressureGradientSource != null) && (!this.BcMap.ContainsPeriodicEdge)) {
                Console.WriteLine("WARNING: Control-File contains extended property 'SourcePressureGradient', but the grid has got no periodic edges.");
            }

            if (this.BcMap.ContainsPeriodicEdge && (Control.PressureGradientSource == null)) {
                Console.WriteLine("WARNING: Grid contains periodic edges, but the extended property 'SourcePressureGradient' has not been used.");
            }
        }

        // Used for calculating current BDF order
        VectorFieldHistory<SinglePhaseField> Velocity;

        public readonly double dt;

        /// <summary>
        /// Minimum of desired BDF order and populated length.
        /// </summary>
        public int BDFOrder {
            get {
                int DesiredOrder = Control.TimeOrder;
                int PopulatedLength = this.Velocity.GetPopulatedLength();
                return Math.Min(DesiredOrder, PopulatedLength);
            }
        }
    }
}
