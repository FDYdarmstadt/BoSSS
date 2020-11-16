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
using System.Runtime.CompilerServices;
using NUnit.Framework.Constraints;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Foundation.IO;

namespace BoSSS.Application.XNSE_Solver
{
    class XNSE : XdgApplicationWithSolver<XNSE_Control>
    {
        LevelSetUpdater lsUpdater;

        protected override LevelSetHandling LevelSetHandling => this.Control.Timestepper_LevelSetHandling;

        int QuadOrder()
        {
            //QuadOrder
            int degU;
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
            return quadOrder;
        }

        protected override LevelSetTracker InstantiateTracker() 
        {
            LevelSet levelSet = new LevelSet(new Basis(GridData, Control.FieldOptions["Phi"].Degree), "Phi");
            levelSet.ProjectField(Control.InitialValues_Evaluators["Phi"]);
            lsUpdater = new LevelSetUpdater((GridData)GridData, Control.CutCellQuadratureType, 1, new string[] { "A", "B" }, levelSet);
            switch (Control.Option_LevelSetEvolution)
            {
                case LevelSetEvolution.Fourier:
                    lsUpdater.AddEvolver("Phi", new FourierEvolver(Control, QuadOrder()));
                    break;
                case LevelSetEvolution.FastMarching:
                case LevelSetEvolution.None:
                    lsUpdater.AddEvolver("Phi", new FastMarcher(Control, QuadOrder()));
                    break;
                default:
                    throw new NotImplementedException();
            }
            //RegisterField(levelSet);
            lsUpdater.Tracker.UpdateTracker(0.0);
            lsUpdater.Tracker.PushStacks();
            return lsUpdater.Tracker;
        }
        
        public override double UpdateLevelset(DGField[] domainFields, double time, double dt, double UnderRelax, bool incremental)
        {
            var DomainVarsDict = new Dictionary<string, DGField>(domainFields.Length);
            for (int iVar = 0; iVar < domainFields.Length; iVar++)
            {
                DomainVarsDict.Add(Operator.DomainVar[iVar], domainFields[iVar]);
            }

            var parameterFields = Timestepping.Parameters;
            var ParameterVarsDict = new Dictionary<string, DGField>(parameterFields.Count());
            for (int iVar = 0; iVar < parameterFields.Count(); iVar++)
            {
                ParameterVarsDict.Add(Operator.ParameterVar[iVar], parameterFields[iVar]);
            }
            double residual = lsUpdater.UpdateLevelSets(DomainVarsDict, ParameterVarsDict, time, dt, UnderRelax, incremental);
            return residual;
        }

        protected override XSpatialOperatorMk2 GetOperatorInstance(int D) {

            int quadOrder = QuadOrder();
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
            opFactory.AddParameter(new Velocity0Mean(D, LsTrk, quadOrder));
            opFactory.AddParameter(new Normals(D, LsTrk));
            opFactory.AddParameter(Curvature.CreateFrom(Control, config, LsTrk, quadOrder));
            
            if (config.isContinuity)
            {
                opFactory.AddEquation(new Continuity(config, D, "A", LsTrk.GetSpeciesId("A"), boundaryMap));
                opFactory.AddEquation(new Continuity(config, D, "B", LsTrk.GetSpeciesId("B"), boundaryMap));
                opFactory.AddEquation(new InterfaceContinuity(config, D, LsTrk));
            }

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
            var parameterFields = Timestepping.Parameters;
            var ParameterVarsDict = new Dictionary<string, DGField>(parameterFields.Count());
            for (int iVar = 0; iVar < parameterFields.Count(); iVar++)
            {
                ParameterVarsDict.Add(Operator.ParameterVar[iVar], parameterFields[iVar]);
            }
            //if (TimestepNo == 1)
            //    lsUpdater.Initialize(ParameterVarsDict, phystime);

            //Update Calls
            dt = GetFixedTimestep();
            Timestepping.Solve(phystime, dt, Control.SkipSolveAndEvaluateResidual);
            return dt;
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 1)
        {
            if(timestepNo == 0)
            {
                this.LsTrk.UpdateTracker(0.0);
                XDGField velX = (XDGField)this.CurrentStateVector.Fields[0];
                var vxA = velX.GetSpeciesShadowField("A");
                //vxA.ProjectField((double[] X) => 0.0);
                var cxB = velX.GetSpeciesShadowField("B");
                //cxB.ProjectField((double[] X) => 20.0);
                Tecplot.PlotFields(new DGField[] { velX, vxA, cxB }, "grid" + timestepNo, physTime, 4);
            }
            Tecplot.PlotFields(this.CurrentStateVector.Fields, "XNSE_Solver" + timestepNo, physTime, superSampling);
        }
    }
}
