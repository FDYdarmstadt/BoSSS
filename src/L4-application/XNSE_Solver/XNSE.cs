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

namespace BoSSS.Application.XNSE_Solver
{
    class XNSE : XdgApplicationWithSolver<XNSE_Control>
    {
        VectorField<XDGField> velocity;

        XDGField pressure;

        LevelSet levelSet;

        protected override LevelSetHandling LevelSetHandling => this.Control.Timestepper_LevelSetHandling;

        protected override IEnumerable<DGField> InstantiateSolutionFields()
        {
            pressure = new XDGField(new XDGBasis(this.LsTrk, this.Control.FieldOptions[VariableNames.Pressure].Degree), VariableNames.Pressure);
            pressure.GetSpeciesShadowField("A").ProjectField(Control.InitialValues_Evaluators["Pressure#A"]);
            pressure.GetSpeciesShadowField("B").ProjectField(Control.InitialValues_Evaluators["Pressure#B"]);
            XDGField vx = new XDGField(new XDGBasis(this.LsTrk, 1), VariableNames.VelocityX);
            vx.GetSpeciesShadowField("A").ProjectField(Control.InitialValues_Evaluators["VelocityX#A"]);
            vx.GetSpeciesShadowField("B").ProjectField(Control.InitialValues_Evaluators["VelocityX#B"]);
            XDGField vy = new XDGField(new XDGBasis(this.LsTrk, 1), VariableNames.VelocityY);
            vy.GetSpeciesShadowField("A").ProjectField(Control.InitialValues_Evaluators["VelocityY#A"]);
            vy.GetSpeciesShadowField("B").ProjectField(Control.InitialValues_Evaluators["VelocityY#B"]);
            velocity = new VectorField<XDGField>(vx, vy);

            return new DGField[] { vx, vy, pressure };
        }

        protected override void SetInitial()
        {
            
        }

        protected override LevelSetTracker InstantiateTracker()
        {
            levelSet = new LevelSet(new Basis(GridData, Control.FieldOptions["Phi"].Degree), VariableNames.LevelSet);
            levelSet.ProjectField(Control.InitialValues_Evaluators["Phi"]);
            LevelSetTracker myLsTrk = new LevelSetTracker((GridData)GridData, XQuadFactoryHelper.MomentFittingVariants.Saye, 1, new string[] { "A", "B"}, levelSet);
            return myLsTrk;
        }

        protected override XSpatialOperatorMk2 GetOperatorInstance(int D)
        {
            int degU = velocity.Max(field => field.Basis.Degree);
            int order = degU * 2 + 1;
            IXNSE_Configuration config = new XNSFE_OperatorConfiguration(this.Control);
            IncompressibleMultiphaseBoundaryCondMap boundaryMap = new IncompressibleMultiphaseBoundaryCondMap(this.GridData, this.Control.BoundaryValues, this.LsTrk.SpeciesNames.ToArray());

            //Build Equations
            SystemOfEquations equationSystem = new SystemOfEquations();
            for (int d = 0; d < D; ++d)
            {
                equationSystem.AddEquation(new NavierStokes("A", d, LsTrk, D, boundaryMap, config));
                equationSystem.AddEquation(new NavierStokes("B", d, LsTrk, D, boundaryMap, config));
                equationSystem.AddEquation(new NSEInterface("A", "B", d, D, boundaryMap, LsTrk, config));
                equationSystem.AddEquation(new NSESurfaceTensionForce("A", "B", d, D, degU, boundaryMap, LsTrk, config));

            }
            equationSystem.AddEquation(new Continuity(config, D, "A", LsTrk.GetSpeciesId("A"), boundaryMap));
            equationSystem.AddEquation(new Continuity(config, D, "B", LsTrk.GetSpeciesId("B"), boundaryMap));
            equationSystem.AddEquation(new InterfaceContinuity(config, D, LsTrk));

            //Get Spatial Operator
            XSpatialOperatorMk2 xSpatialOperator = equationSystem.GetSpatialOperator(order);

            xSpatialOperator.Commit();
            return xSpatialOperator;
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt)
        {
            //Update Calls
            dt = base.GetFixedTimestep();
            Timestepping.Solve(phystime, 1E100, true);
            return dt;
        }
    }
}
