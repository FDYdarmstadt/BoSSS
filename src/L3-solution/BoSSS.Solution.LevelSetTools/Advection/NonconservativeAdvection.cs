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
using BoSSS.Solution.TimeStepping;

namespace BoSSS.Solution.LevelSetTools.Advection {

    /// <summary>
    /// Nonconservative Advection: 
    /// \f[
    /// \textrm{d}_t \varphi + \vec{u} \cdot \operatorname{grad} \varphi = 0, \qquad
    /// \operatorname{div} \vec{u} \neq 0
    /// \f]
    /// Now with Runge-Kutta Timestepping!
    /// </summary>
    public class ExplicitNonconservativeAdvection : VectorVelocityLevelSetAdvection  {
        
        double m_dt;
        RungeKutta TimeEvo;
                
        VectorFieldHistory<SinglePhaseField> Velocity;
        VectorField<SinglePhaseField> VelocityInterpolation;

        /// <summary>
        /// Create Operators for Nonconservative Advection
        /// </summary>
        public ExplicitNonconservativeAdvection(SinglePhaseField LevelSet, 
            VectorFieldHistory<SinglePhaseField> Velocity, 
            IncompressibleBoundaryCondMap BcMap, bool AssumeDivergenceFreeVelocity = false) 
            : base (LevelSet, BcMap, AssumeDivergenceFreeVelocity) {

            this.Velocity = Velocity;

            VelocityInterpolation = new VectorField<SinglePhaseField>(D.ForLoop(d => new SinglePhaseField(Velocity.Current[d].Basis, "VelocityInterpolation" + d)));

            MeanVelocity = new VectorField<SinglePhaseField>(D.ForLoop(d => new SinglePhaseField(new Basis(GridDat, 0), VariableNames.Velocity0MeanVector(D)[d])));

            SpatialOperator Lsevo = CreateAdvectionSpatialOperator(BcMap);

            DGField[] ParamFields;
            if (!divUzero) {
                ParamFields = ArrayTools.Cat(this.Velocity.Current, this.MeanVelocity, this.divU);
            } else {
                ParamFields = ArrayTools.Cat(this.Velocity.Current, this.MeanVelocity);
            }

            TimeEvo = new RungeKutta(RungeKuttaScheme.TVD3, Lsevo, this.LevelSet.Mapping, new CoordinateMapping(ParamFields));

            TimeEvo.OnBeforeComputeChangeRate += RKParamterUpdate;

        }

         /// <summary>
        /// Perform Advection
        /// </summary>
        /// <param name="dt"> Timestepsize</param>
        public override void Advect(double dt) {
            m_dt = dt;
            MeanVelocity.Clear();
            MeanVelocity.AccLaidBack(1.0, Velocity.Current);

            TimeEvo.Perform(dt);
        }

        public override void FinishTimeStep() {
            // Do nothing
        }



        void RKParamterUpdate(double AbsTime, double RelTime) {
            var gDat = this.LevelSet.GridDat;
            int D = gDat.SpatialDimension;

            // Compute level set gradient
            // ==========================

            var PhiGradient = new VectorField<SinglePhaseField>(D, LevelSet.Basis, SinglePhaseField.Factory);
            PhiGradient.GradientByFlux(1.0, this.LevelSet);

            // linear interpolation of filtered velocity
            // =========================================


            this.VelocityInterpolation.Clear();
            double a0 = 1.0 - RelTime / this.m_dt, a1 = RelTime / this.m_dt;
            if (a0 != 0.0)
                this.VelocityInterpolation.Acc(a0, this.Velocity[0]);
            if (a1 != 0.0)
                this.VelocityInterpolation.Acc(a1, this.Velocity[1]);


            // Projection in non-cut cells
            // ===========================
            if (!divUzero) { 
                this.divU.Clear();
                // Compute Divergence locally,
                // the values at the singularities might spoil accuracy in the other cells,
                // if a flux formulation is used
                this.divU.Divergence(1.0, this.VelocityInterpolation);
            }
        }
    }

    /// <summary>
    /// Nonconservative Advection:
    /// \f[
    /// \textrm{d}_t \varphi + \vec{u} \cdot \operatorname{grad} \varphi = 0, \qquad
    /// \operatorname{div} \vec{u} \neq 0
    /// \f]
    /// With Crank-Nicholson Timestepping:
    /// 
    /// u_new/dt + 1/2 \vec{u_new} \cdot \operatorname{grad} \varphi_new = u_old/dt- 1/2 \vec{u_old} \cdot \operatorname{grad} \varphi_old
    /// 
    /// 
    /// </summary>
    public class BDFNonconservativeAdvection : VectorVelocityLevelSetAdvection {

        internal VectorField<SinglePhaseField> Velocity;
        internal SinglePhaseField OldRHS;

        BDFTimestepper myBDFTimestepper;

        /// <summary>
        /// Setup and initial evaluation of RHS
        /// </summary>
        /// <param name="LevelSet"></param>
        /// <param name="Velocity"></param>
        /// <param name="bcMap">Boundary Conditions for LevelSet</param>
        /// <param name="AssumeDivergenceFreeVelocity">Switch for the source term on the rhs arising from a non divergence-free velocity</param>
        public BDFNonconservativeAdvection(SinglePhaseField LevelSet, VectorField<SinglePhaseField> Velocity, IncompressibleBoundaryCondMap bcMap, int BDForder, SubGrid subGrid = null,  bool AssumeDivergenceFreeVelocity = false)
            : base(LevelSet, bcMap, AssumeDivergenceFreeVelocity) {

            this.Velocity = Velocity;
            this.OldRHS = LevelSet.CloneAs();
            this.SO = CreateAdvectionSpatialOperator(bcMap);

            if (Velocity == null) {
                throw new ArgumentException("Velocity Field not initialized!");
            }

            MeanVelocity = new VectorField<SinglePhaseField>(D.ForLoop(d => new SinglePhaseField(new Basis(GridDat, 0), VariableNames.Velocity0MeanVector(D)[d])));
            MeanVelocity.Clear();
            MeanVelocity.AccLaidBack(1.0, Velocity);

            myBDFTimestepper = new BDFTimestepper(SO, new List<DGField>() { LevelSet }, ArrayTools.Cat(this.Velocity, this.MeanVelocity, this.divU), BDForder, () => new ilPSP.LinSolvers.MUMPS.MUMPSSolver(),false, subGrid);

        }

        /// <summary>
        /// Setup and solve the System, No evaluation of RHS, this is done in <see cref="FinishTimeStep"/>
        /// </summary>
        /// <param name="dt">TimestepSize</param>
        public override void Advect(double dt) {
            if (!divUzero) {
                this.divU.Clear();
                // Compute Divergence locally,
                // the values at the singularities might spoil accuracy in the other cells,
                // if a flux formulation is used
                this.divU.Divergence(1.0, this.Velocity);
            }

            MeanVelocity.Clear();
            MeanVelocity.AccLaidBack(1.0, Velocity);

            var OldLevelSet = LevelSet.CloneAs();
            myBDFTimestepper.Perform(dt);
            OldLevelSet.Acc(-1.0, LevelSet);
        }

        /// <summary>
        /// Update the RHS of the Crank-Nicholson after the timestep is finished
        /// NOTE: only call this after the timestep,
        /// NOT after each time <see cref="Advect(double)"/> is called,
        /// since this may be included in a iteration loop
        /// </summary>
        public override void FinishTimeStep() {
            // Update the RHS according to the new level-set
            myBDFTimestepper.FinishTimeStep();
        }
    }







    /// <summary>
    /// This class provides the spatial operators, but no timestepping scheme
    /// </summary>
    public abstract class VectorVelocityLevelSetAdvection: ILevelSetAdvection{
        internal GridData GridDat;
        internal int D;
        internal VectorField<SinglePhaseField> MeanVelocity;
        internal SinglePhaseField LevelSet;
        internal SpatialOperator SO;
        internal SinglePhaseField divU;
        internal bool divUzero;
        
        /// <summary>
        /// ctr
        /// </summary>
        /// <param name="LevelSet">the level-set which is to be moved</param>
        /// <param name="bcMap">boundary conditions</param>
        /// <param name="AssumeDivergenceFreeVelocity">switch for the nonlinear term on the rhs due to non-divergence free velocity fields</param>
        public VectorVelocityLevelSetAdvection(SinglePhaseField LevelSet, IncompressibleBoundaryCondMap bcMap, bool AssumeDivergenceFreeVelocity = false) {
            GridDat = (GridData)(LevelSet.Basis.GridDat);
            D = GridDat.SpatialDimension;
            this.LevelSet = LevelSet;
            divUzero = AssumeDivergenceFreeVelocity;
            divU = new SinglePhaseField(LevelSet.Basis);

        }

        public abstract void Advect(double dt);
        public abstract void FinishTimeStep();

        /// <summary>
        /// Spatial Operators and Matrices
        /// </summary>
        /// with
        /// <param name="bcMap">Boundary Conditions</param>
        public SpatialOperator CreateAdvectionSpatialOperator(IncompressibleBoundaryCondMap bcMap) {

            Func<int[], int[], int[], int> QuadOrderFunction = QuadOrderFunc.SumOfMaxDegrees();


            string[] parameterList;
            parameterList = ArrayTools.Cat(VariableNames.Velocity0Vector(D),
                    VariableNames.Velocity0MeanVector(D));
            if (!divUzero) {
                parameterList = ArrayTools.Cat(parameterList, "div(U)");
            }

            SpatialOperator SO = new SpatialOperator(new string[] { "LevelSet" },
                    parameterList,
                    new string[] { "Phi-Evo" },
                    QuadOrderFunc.NonLinear(2));
            //div(u.phi)
            //SO.EquationComponents["Phi-Evo"].Add(new LevelSetUpwindFlux(GridDat, bcMap));
            SO.EquationComponents["Phi-Evo"].Add(new LevelSetLLFFlux(GridDat, bcMap));
            //bcMap.PhysMode = PhysicsMode.Multiphase;
            //SO.EquationComponents["Phi-Evo"].Add(new LinearizedScalarConvection(D, bcMap,null));
            //SO.EquationComponents["Phi-Evo"].Add(new LevelSetAdvectionCentralFlux(D));
            //-phi*div(u)
            if (!divUzero) {
                SO.EquationComponents["Phi-Evo"].Add(new FextSource());
            }
            //penalization
            //Lsevo.EquationComponents["Phi-Evo"].Add(new JumpPenalization.GradientJumpForm2());

            SO.Commit();
            return SO;
        }

    }

}
