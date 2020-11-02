using BoSSS.Foundation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.XNSECommon;
using BoSSS.Solution.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using ilPSP.Utils;
using System.Diagnostics;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Control;

namespace BoSSS.Application.XNSE_Solver
{
    class XNSE : XdgApplicationWithSolver<XNSE_Control>
    {
        LevelSet levelSet;

        protected override LevelSetHandling LevelSetHandling => this.Control.Timestepper_LevelSetHandling;

        protected override LevelSetTracker InstantiateTracker() {
            levelSet = new LevelSet(new Basis(GridData, Control.FieldOptions["Phi"].Degree), "Phi");
            levelSet.ProjectField(Control.InitialValues_Evaluators["Phi"]);
            LevelSetTracker myLsTrk = new LevelSetTracker((GridData)GridData, Control.CutCellQuadratureType, 1, new string[] { "A", "B" }, levelSet);
            return myLsTrk;
        }

        VectorField<XDGField> gravity;
        
        protected override IEnumerable<DGField> CreateAdditionalFields() {
            int D = this.GridData.SpatialDimension;
            gravity = new VectorField<XDGField>(D.ForLoop(d => new XDGField(this.CurrentState.BasisS[d] as XDGBasis, VariableNames.Gravity_d(d))));
            return gravity;
        }

        protected override XSpatialOperatorMk2 GetOperatorInstance(int D) {

            //QuadOrder
            int degU = default;
            if (Control.FieldOptions.TryGetValue("Velocity*", out FieldOpts field))
            {
                degU = field.Degree;
            }
            else if (Control.FieldOptions.TryGetValue(BoSSS.Solution.NSECommon.VariableNames.VelocityX, out FieldOpts field1))
            {
                degU = field1.Degree;
            }
            else
            {
                throw new Exception("Velocity not found!");
            }
            int quadOrder = degU * (this.Control.PhysicalParameters.IncludeConvection ? 3 : 2);
            if (this.Control.solveKineticEnergyEquation)
                quadOrder *= 2;


            XNSFE_OperatorConfiguration config = new XNSFE_OperatorConfiguration(this.Control);
            IncompressibleMultiphaseBoundaryCondMap boundaryMap = new IncompressibleMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, this.LsTrk.SpeciesNames.ToArray());

            //Build Equations
            OperatorFactory opFactory = new OperatorFactory();
            for (int d = 0; d < D; ++d) {
                opFactory.AddEquation(new NavierStokes("A", d, LsTrk, D, boundaryMap, config));
                opFactory.AddParameter(Gravity.CreateFrom("A", d, D, Control, Control.PhysicalParameters.rho_A));
                opFactory.AddEquation(new NavierStokes("B", d, LsTrk, D, boundaryMap, config));
                opFactory.AddParameter(Gravity.CreateFrom("B", d, D, Control, Control.PhysicalParameters.rho_B));
                opFactory.AddEquation(new NSEInterface("A", "B", d, D, boundaryMap, LsTrk, config));
                opFactory.AddEquation(new NSESurfaceTensionForce("A", "B", d, D, boundaryMap, LsTrk, config));
            }
            opFactory.AddParameter(new Velocity0(D));
            opFactory.AddParameter(new Velocity0Mean(D, LsTrk));
            opFactory.AddParameter(new Normals(D, LsTrk));
            opFactory.AddParameter(Curvature.CreateFrom(Control, config, LsTrk, quadOrder));
            
            opFactory.AddEquation(new Continuity(config, D, "A", LsTrk.GetSpeciesId("A"), boundaryMap));
            opFactory.AddEquation(new Continuity(config, D, "B", LsTrk.GetSpeciesId("B"), boundaryMap));
            opFactory.AddEquation(new InterfaceContinuity(config, D, LsTrk));


            //Get Spatial Operator
            XSpatialOperatorMk2 XOP = opFactory.GetSpatialOperator(quadOrder);
            //final settings
            //XOP.LinearizationHint = LinearizationHint.AdHoc;
            XOP.Commit();

            // test the ordering
            Debug.Assert(XOP.DomainVar.IndexOf(VariableNames.VelocityX) < XOP.DomainVar.IndexOf(VariableNames.VelocityY));
            Debug.Assert(XOP.DomainVar.IndexOf(VariableNames.VelocityY) < XOP.DomainVar.IndexOf(VariableNames.Pressure));
            Debug.Assert(XOP.CodomainVar.IndexOf(EquationNames.MomentumEquationX) < XOP.CodomainVar.IndexOf(EquationNames.MomentumEquationY));
            Debug.Assert(XOP.CodomainVar.IndexOf(EquationNames.MomentumEquationY) < XOP.CodomainVar.IndexOf(EquationNames.ContinuityEquation));
            if(D > 2) {
                Debug.Assert(XOP.DomainVar.IndexOf(VariableNames.VelocityY) < XOP.DomainVar.IndexOf(VariableNames.VelocityZ));
                Debug.Assert(XOP.DomainVar.IndexOf(VariableNames.VelocityZ) < XOP.DomainVar.IndexOf(VariableNames.Pressure));
                Debug.Assert(XOP.CodomainVar.IndexOf(EquationNames.MomentumEquationY) < XOP.CodomainVar.IndexOf(EquationNames.MomentumEquationZ));
                Debug.Assert(XOP.CodomainVar.IndexOf(EquationNames.MomentumEquationZ) < XOP.CodomainVar.IndexOf(EquationNames.ContinuityEquation));
            }
            return XOP;
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt)
        {
            //Update Calls
            dt = GetFixedTimestep();
            Timestepping.Solve(phystime, dt, Control.SkipSolveAndEvaluateResidual);
            return dt;
        }
    }
}
