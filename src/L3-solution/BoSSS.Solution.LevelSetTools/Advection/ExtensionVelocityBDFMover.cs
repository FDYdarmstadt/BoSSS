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
using ilPSP;
using ilPSP.Utils;
using BoSSS.Platform;
using BoSSS.Foundation;
using BoSSS.Solution.TimeStepping;
using BoSSS.Solution.NSECommon;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.LevelSetTools.EllipticExtension;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.Grid;

namespace BoSSS.Solution.LevelSetTools.Advection {
    /// <summary>
    /// A stripped-down class providing the motion algorithm for Extension Velocity
    /// The Extension Velocity is based on the density-averaged value at the interface
    /// </summary>
    public class ExtensionVelocityBDFMover:ILevelSetAdvection {

        Extender[] VelocityExtender;
        VectorField<DGField> Velocity;
        SinglePhaseField LevelSet;
        bool nearfield;
        VectorField<SinglePhaseField> VectorExtension;
        VectorField<SinglePhaseField> LevelSetGradient;

        LevelSetTracker LSTrk;


        SinglePhaseField OldRHS;
        GridData GridDat;
        int D;
        // Mean Velocity for choosing Flux Direction
        VectorField<SinglePhaseField> MeanVelocity;
        SpatialOperator AdvectionSpatialOperator;
        DGField divU;

        BDFTimestepper myBDFTimestepper;

        SubGrid subGrid;

        public ExtensionVelocityBDFMover(LevelSetTracker LSTrk,
                    SinglePhaseField LevelSet,
                    VectorField<SinglePhaseField> LevelSetGradient,
                    VectorField<DGField> Velocity,
                    EllipticExtVelAlgoControl Control,
                    IncompressibleBoundaryCondMap bcMap,
                    int BDForder,
                    VectorField<SinglePhaseField> VectorExtension,
                    double[] Density = null,
                    bool AssumeDivergenceFreeVelocity = false, SubGrid subGrid = null){
            this.GridDat = LSTrk.GridDat;
            D = GridDat.SpatialDimension;
            this.LevelSetGradient = LevelSetGradient;
            this.LSTrk = LSTrk;
            this.LevelSet = LevelSet;

            this.Velocity = Velocity;
            this.OldRHS = LevelSet.CloneAs();
            this.AdvectionSpatialOperator = CreateAdvectionSpatialOperator(bcMap);

            this.subGrid = subGrid;
            this.nearfield = subGrid != null;

            Basis NonXVelocityBasis;
            if (Velocity == null) {
                throw new ArgumentException("Velocity Field not initialized!");
            }

            // Initialize Extension Velocity Algorithm
            double PenaltyBase = Control.PenaltyMultiplierInterface * ((double)((LevelSet.Basis.Degree + 1) * (LevelSet.Basis.Degree + D))) / ((double)D);
            ILevelSetComponent InterfaceFlux;

            

            //VectorExtension = new VectorField<SinglePhaseField>(D, Velocity[0].Basis, "ExtVel", SinglePhaseField.Factory);
            if (Velocity[0].GetType() == typeof(SinglePhaseField)) {
                NonXVelocityBasis = ((SinglePhaseField)Velocity[0]).Basis;
                InterfaceFlux = new SingleComponentInterfaceForm(PenaltyBase, LSTrk);
            }
            else if (Velocity[0].GetType() == typeof(XDGField)) {
                NonXVelocityBasis = ((XDGField)Velocity[0]).Basis.NonX_Basis;
                InterfaceFlux = new DensityWeightedExtVel(PenaltyBase, LSTrk, Density);
            }
            else {
                throw new ArgumentException("VelocityField must be either a SinglePhaseField or a XDGField!");
            };
            //VectorExtension = new VectorField<SinglePhaseField>(D, NonXVelocityBasis, "ExtVel", SinglePhaseField.Factory);
            this.VectorExtension = VectorExtension;



            VelocityExtender = new Extender[D];
            for (int d = 0; d < D; d++) {
                VelocityExtender[d] = new Extender(VectorExtension[d], LSTrk, InterfaceFlux, new List<DGField> { Velocity[d] }, LevelSetGradient, Control);
                VelocityExtender[d].ConstructExtension(new List<DGField> { Velocity[d] }, Control.subGridRestriction);
            }
#if DEBUG
            VectorExtension.CheckForNanOrInf();
#endif

            // Initialize Advection Algorithm
            divU = new SinglePhaseField(NonXVelocityBasis);
            divU.Identification = "Divergence";
            divU.Clear();
            divU.Divergence(1.0, VectorExtension);
            MeanVelocity = new VectorField<SinglePhaseField>(D.ForLoop(d => new SinglePhaseField(new Basis(GridDat, 0), VariableNames.Velocity0MeanVector(D)[d])));
            MeanVelocity.Clear();
            MeanVelocity.AccLaidBack(1.0, VectorExtension);
            myBDFTimestepper = new BDFTimestepper(AdvectionSpatialOperator, new List<DGField>() { LevelSet }, ArrayTools.Cat(VectorExtension, this.MeanVelocity, this.divU), BDForder, Control.solverFactory, false, subGrid);
        }

        SpatialOperator CreateAdvectionSpatialOperator(IncompressibleBoundaryCondMap bcMap) {

            Func<int[], int[], int[], int> QuadOrderFunction = QuadOrderFunc.SumOfMaxDegrees();

            string[] parameterList;
            parameterList = ArrayTools.Cat(
                VariableNames.Velocity0Vector(D),
                VariableNames.Velocity0MeanVector(D),
                "div(U)");

            SpatialOperator SO = new SpatialOperator(new string[] { "LevelSet" },
                    parameterList,
                    new string[] { "Phi-Evo" },
                    QuadOrderFunc.NonLinear(2));
            //div(u*phi)
            SO.EquationComponents["Phi-Evo"].Add(new LevelSetLLFFlux(GridDat, bcMap));
            //-phi*div(u)
            SO.EquationComponents["Phi-Evo"].Add(new FextSource());
            SO.Commit();
            return SO;
        }

        /// <summary>
        /// Setup and solve the System, No evaluation of RHS, this is done in <see cref="FinishTimeStep"/>
        /// </summary>
        /// <param name="dt">TimestepSize</param>
        public void Advect(double dt) {
            // Construct the Extension Velocity
            VectorExtension.Clear();
            for (int d = 0; d < D; d++) {
                VelocityExtender[d].ConstructExtension(new List<DGField> { Velocity[d] }, this.nearfield);
            }
            VectorExtension.CheckForNanOrInf(true, true,true);
            // Advect with that Extension Velocity
            this.divU.Clear();

            // Compute Divergence locally,
            // the values at the singularities might spoil accuracy in the other cells,
            // if a flux formulation is used
            this.divU.Divergence(1.0, this.VectorExtension, subGrid?.VolumeMask);
                        
            MeanVelocity.Clear();
            MeanVelocity.AccLaidBack(1.0, VectorExtension, subGrid?.VolumeMask);

            myBDFTimestepper.Perform(dt);
            
        }

        /// <summary>
        /// Update the RHS of the Crank-Nicholson after the timestep is finished
        /// NOTE: only call this after the timestep,
        /// NOT after each time <see cref="Advect(double)"/> is called,
        /// since this may be included in a iteration loop
        /// </summary>
        public void FinishTimeStep() {
            // Update the RHS according to the new level-set
            myBDFTimestepper.FinishTimeStep();
        }
    }
}
