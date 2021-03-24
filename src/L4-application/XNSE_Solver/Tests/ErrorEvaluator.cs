using BoSSS.Application.XNSE_Solver;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Control;
using BoSSS.Solution.XdgTimestepping;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ilPSP;
using BoSSS.Solution.NSECommon;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Utils;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Foundation.Quadrature;
using BoSSS.Solution;
using BoSSS.Solution.XNSECommon;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Application.XNSE_Solver.Tests
{

    /// <summary>
    /// base class for error evaluation within tests
    /// </summary>
    public abstract class ErrorEvaluator {

        protected IApplication solver;

        public ErrorEvaluator(IApplication solver) {
            this.solver = solver;
        }

        /// <summary>
        /// Computes the L2 Error of the actual solution against the exact solution in the control object 
        /// </summary>
        public abstract double[] ComputeL2Error(double time, XNSE_Control control);


        /// <summary>
        /// number of cells in the grid/mesh
        /// </summary>
        public double GetGridCells() {
            return solver.GridData.CellPartitioning.TotalLength;
        }

        /// <summary>
        /// characteristic grid/mesh length scale
        /// </summary>
        public double GetGrid_h() {
            return ((Foundation.Grid.Classic.GridData)(solver.GridData)).Cells.h_minGlobal;
        }
    }

    /// <summary>
    /// error evaluator for XNSE-based tests 
    /// computes error for: velocity, pressure 
    /// </summary>
    /// <remarks>
    /// Seems redundant with <see cref="L2ErrorLogger"/>
    /// </remarks>
    public class XNSEErrorEvaluator<T> : ErrorEvaluator where T: XNSE_Control, new() {

        public XNSEErrorEvaluator(IApplication solver) : base(solver){

        }

        new protected ApplicationWithSolver<T> solver {
            get {
                return (ApplicationWithSolver<T>)(base.solver);
            }
        }

        public double[] ComputeVelocityError(IDictionary<string, Func<double[], double, double>[]> exactVelocity, double time)
        {
            int D = solver.GridData.SpatialDimension;
            var FluidSpecies = exactVelocity.Keys.ToArray();
            double[] Ret = new double[D];

            int order = 0;
            if (solver.LsTrk.GetCachedOrders().Count > 0)
            {
                order = solver.LsTrk.GetCachedOrders().Max();
            }
            else
            {
                order = 1;
            }

            var SchemeHelper = solver.LsTrk.GetXDGSpaceMetrics(solver.LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;

            Dictionary<string, double[]> L2Error_Species = new Dictionary<string, double[]>();
            double[] L2Error = new double[D];

            foreach (var spc in FluidSpecies)
            {
                L2Error_Species.Add(spc, new double[D]);

                SpeciesId spId = solver.LsTrk.GetSpeciesId(spc);
                var scheme = SchemeHelper.GetVolumeQuadScheme(spId);


                for (int d = 0; d < D; d++)
                {
                    string velocity = VariableNames.VelocityVector(D)[d];
                    ConventionalDGField Vel_d = ((XDGField)solver.CurrentStateVector.Mapping.Single(Field => Field.Identification == velocity)).GetSpeciesShadowField(spc);

                    L2Error_Species[spc][d] = Vel_d.L2Error(exactVelocity[spc][d].Vectorize(time), order, scheme);
                    L2Error[d] += L2Error_Species[spc][d].Pow2();

                    solver.QueryHandler.ValueQuery("L2err_" + VariableNames.Velocity_d(d) + "#" + spc, L2Error_Species[spc][d], true);
                }
            }
            L2Error = L2Error.Select(x => x.Sqrt()).ToArray();

            for (int d = 0; d < D; d++)
            {
                solver.QueryHandler.ValueQuery("L2err_" + VariableNames.Velocity_d(d), L2Error[d], true);
                Ret[d] = L2Error[d];
            }
            return L2Error;
        }

        public double ComputePressureError(IDictionary<string, Func<double[], double, double>> exactPressure, double time)
        {
            int D = solver.GridData.SpatialDimension;
            var FluidSpecies = exactPressure.Keys.ToArray();

            string pressureName = VariableNames.Pressure;
            XDGField pressure = (XDGField)solver.CurrentStateVector.Mapping.Single(Field => Field.Identification == pressureName);

            int order = 0;
            if (solver.LsTrk.GetCachedOrders().Count > 0) {
                order = solver.LsTrk.GetCachedOrders().Max();
            } else {
                order = 1;
            }
            //order = Math.Max(pressure.Basis.Degree * 2, order);
            Console.WriteLine("pressure error with order " + order);

            var SchemeHelper = solver.LsTrk.GetXDGSpaceMetrics(solver.LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;

            // pass 1: mean value of pressure difference
            double DiffInt = 0;
            double Volume = 0;
            foreach (var spc in FluidSpecies)
            {

                SpeciesId spId = solver.LsTrk.GetSpeciesId(spc);
                var scheme = SchemeHelper.GetVolumeQuadScheme(spId);
                var rule = scheme.Compile(solver.GridData, order);
               
                DiffInt += pressure.GetSpeciesShadowField(spc).LxError(exactPressure[spc].Vectorize(time), (X, a, b) => (a - b), rule);
                Volume += pressure.GetSpeciesShadowField(spc).LxError(exactPressure[spc].Vectorize(time), (X, a, b) => 1.0, rule);
            }
            double PressureDiffMean = DiffInt / Volume;
            Console.WriteLine("Mean Pressure diff: " + PressureDiffMean);

            double L2Error = 0;
            Dictionary<string, double> L2Error_Species = new Dictionary<string, double>();

            foreach (var spc in FluidSpecies)
            {

                SpeciesId spId = solver.LsTrk.GetSpeciesId(spc);
                var scheme = SchemeHelper.GetVolumeQuadScheme(spId);
                var rule = scheme.Compile(solver.GridData, order);

                double IdV = pressure.GetSpeciesShadowField(spc).LxError(exactPressure[spc].Vectorize(time), (X, a, b) => (a - b - PressureDiffMean).Pow2(), rule);
                L2Error += IdV;
                L2Error_Species.Add(spc, IdV.Sqrt());

                solver.QueryHandler.ValueQuery("L2err_" + VariableNames.Pressure + "#" + spc, L2Error_Species[spc], true);
            }


            L2Error = L2Error.Sqrt();
            solver.QueryHandler.ValueQuery("L2err_" + VariableNames.Pressure, L2Error, true);
            return L2Error;
        }


        /// <summary>
        /// Computes the L2 Error of the actual solution against the exact solution in the control object 
        /// (<see cref="XNSE_Control.ExactSolutionVelocity"/> and <see cref="XNSE_Control.ExactSolutionPressure"/>).
        /// </summary>
        public override double[] ComputeL2Error(double time, XNSE_Control control) {
            int D = solver.GridData.SpatialDimension;
            double[] Ret = new double[D + 1];

            if(control.ExactSolutionVelocity != null) {
                double[] error = ComputeVelocityError(control.ExactSolutionVelocity, time);
                error.CopyTo(Ret, 0);
            }
            if(control.ExactSolutionPressure != null) {
                double error = ComputePressureError(control.ExactSolutionPressure, time);
                Ret[D] = error;
            }
            return Ret;
        }
    }

    /// <summary>
    /// error evaluator for LevelSet-based tests 
    /// computes error fields for: Phi, PhiDG, gradient of PhiDG 
    /// integral values: area, length
    /// </summary>
    public class LevelSetErrorEvaluator<T> : ErrorEvaluator where T: XNSE_Control, new() {

        public LevelSetErrorEvaluator(IApplication solver) : base(solver) {

        }

        new protected SolverWithLevelSetUpdater<T> solver {
            get {
                return (SolverWithLevelSetUpdater<T>)(base.solver);
            }
        }


        LevelSet exactPhi;

        LevelSetTracker exactLsTrk;

        void SetExactPhiAndLevelSetTracker(double time) {
            // exact level-set field
            Func<double[], double, double> phiExactFunc = solver.Control.Phi;
            exactPhi = new LevelSet(new Basis(solver.GridData, 8), "exactLevelSet");
            exactPhi.Clear();
            exactPhi.ProjectField(NonVectorizedScalarFunction.Vectorize(phiExactFunc, time));
            // exact level-set tracker
            exactLsTrk = new LevelSetTracker((GridData)solver.GridData,
            XQuadFactoryHelper.MomentFittingVariants.Saye, 1, solver.LsTrk.SpeciesNames.ToArray(), exactPhi);
            exactLsTrk.UpdateTracker(time);
        }


        /// <summary>
        /// computes the error against the continuous level set field "Phi"
        /// </summary>
        /// <param name="exactLevelSetFunc"></param>
        /// <param name="time"></param>
        /// <param name="cm"></param>
        /// <returns></returns>
        public double ComputeLevelSetError(CellMask cm) {

            SinglePhaseField PhiCG = solver.LsUpdater.LevelSets[VariableNames.LevelSetCG].CGLevelSet;

            double L2Error = PhiCG.L2Error(exactPhi, cm);

            solver.QueryHandler.ValueQuery("L2err_Phi", L2Error, true);

            return L2Error;
        }

        /// <summary>
        /// computes the error against the discontinuous level set field "PhiDG"
        /// </summary>
        /// <param name="exactLevelSetFunc"></param>
        /// <param name="time"></param>
        /// <param name="cm"></param>
        /// <returns></returns>
        public double ComputeDGLevelSetError(CellMask cm) {


            SinglePhaseField PhiDG = solver.LsUpdater.LevelSets[VariableNames.LevelSetCG].DGLevelSet;

            double L2Error = PhiDG.L2Error(exactPhi, cm);

            solver.QueryHandler.ValueQuery("L2err_PhiDG", L2Error, true);

            return L2Error;
        }

        /// <summary>
        /// computes the error of the signed distance property
        /// </summary>
        /// <param name="cm"></param>
        /// <returns></returns>
        public double ComputeDGLevelSetGradientError(CellMask cm) {

            SinglePhaseField PhiDG = solver.LsUpdater.LevelSets[VariableNames.LevelSetCG].DGLevelSet;

            int D = solver.GridData.SpatialDimension;
            var GradientPhiDG = new VectorField<SinglePhaseField>(D.ForLoop(d => new SinglePhaseField(PhiDG.Basis, "dPhiDG_dx[" + d + "]")));
            GradientPhiDG.Gradient(1.0, PhiDG, cm);

            double L2Norm = GradientPhiDG.L2Norm();
            double L2error = L2Norm - 1.0;

            return L2error;
        }

        /// <summary>
        /// computes the error of the interface points with respect to the given interface form
        /// </summary>
        /// <returns></returns>
        public double ComputeInterfacePointsError(double time) {

            SinglePhaseField PhiCG = solver.LsUpdater.LevelSets[VariableNames.LevelSetCG].CGLevelSet;
            SubGrid sbgrd = solver.LsTrk.Regions.GetCutCellSubGrid();
            MultidimensionalArray interfaceP = XNSEUtils.GetInterfacePoints(solver.LsTrk, PhiCG, sbgrd);

            Func<double[], double, double> phiExact = solver.Control.Phi;
            double error = 0.0;
            for (int i = 0; i < interfaceP.Lengths[0]; i++) {
                double dist = phiExact(interfaceP.ExtractSubArrayShallow(i, -1).To1DArray(), time);
                error += dist.Pow2();
            }

            return error.Sqrt();
        }


        /// <summary>
        /// computes the interface length in 2D and area in 3D
        /// </summary>
        /// <returns></returns>
        public double ComputeInterfaceSizeError() {

            double exactInterfaceSize = XNSEUtils.GetInterfaceLength(exactLsTrk);
            double interfaceSize = XNSEUtils.GetInterfaceLength(solver.LsTrk);

            return Math.Abs(interfaceSize- exactInterfaceSize);
        }

        /// <summary>
        /// computes the species area
        /// </summary>
        /// <returns></returns>
        public Dictionary<SpeciesId, double> ComputeSpeciesDomainSizeError() {

            Dictionary<SpeciesId, double> spcArea = new Dictionary<SpeciesId, double>();

            foreach (SpeciesId spcId in solver.LsTrk.SpeciesIdS) {
                double exactArea = XNSEUtils.GetSpeciesArea(exactLsTrk, spcId);
                double area = XNSEUtils.GetSpeciesArea(solver.LsTrk, spcId);
                double error = Math.Abs(area - exactArea);
                spcArea.Add(spcId, error);
            }

            return spcArea;
        }


        /// <summary>
        /// Computes the L2 Error of the actual solution against the exact solution in the control object 
        /// (<see cref="XNSE_Control.ExactSolutionVelocity"/> and <see cref="XNSE_Control.ExactSolutionPressure"/>).
        /// </summary>
        public override double[] ComputeL2Error(double time, XNSE_Control control) {

            SetExactPhiAndLevelSetTracker(time);

            CellMask cm = solver.LsTrk.Regions.GetCutCellMask();
            //CellMask cm = solver.LsTrk.Regions.GetNearFieldMask(1);

            var spcIds = solver.LsTrk.SpeciesIdS;

            double[] Ret = new double[4 + spcIds.Count];

            if (control.Phi != null) {
                //Ret[0] = ComputeInterfacePointsError(time);
                Ret[0] = ComputeLevelSetError(cm);
                Ret[1] = ComputeDGLevelSetError(cm);
            }
            Ret[2] = ComputeDGLevelSetGradientError(cm);
            Ret[3] = ComputeInterfaceSizeError();

            Dictionary<SpeciesId, double> spcArea = ComputeSpeciesDomainSizeError();
            int n = 4;
            foreach (SpeciesId spc in spcIds) {
                Ret[n] = spcArea[spc];
                n++;
            }

            return Ret;
        }

    }


    class XHeatErrorEvaluator<T> : ErrorEvaluator where T: XNSE_Control, new() {

        public XHeatErrorEvaluator(IApplication solver) : base(solver) {

        }
        new protected SolverWithLevelSetUpdater<T> solver {
            get {
                return (SolverWithLevelSetUpdater<T>)(base.solver);
            }
        }


        public double ComputeTemperatureError(IDictionary<string, Func<double[], double, double>> exactTemperature, double time) {
            int D = solver.GridData.SpatialDimension;

            int order = 0;
            if (solver.LsTrk.GetCachedOrders().Count > 0) {
                order = solver.LsTrk.GetCachedOrders().Max();
            } else {
                order = 1;
            }

            var SchemeHelper = solver.LsTrk.GetXDGSpaceMetrics(solver.LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;

            double L2Error = 0;
            Dictionary<string, double> L2Error_Species = new Dictionary<string, double>();

            foreach (var spc in solver.LsTrk.SpeciesNames) {

                SpeciesId spId = solver.LsTrk.GetSpeciesId(spc);
                var scheme = SchemeHelper.GetVolumeQuadScheme(spId);

                string temperatureName = VariableNames.Temperature;
                ConventionalDGField temperature = ((XDGField)solver.CurrentStateVector.Mapping.Single(Field => Field.Identification == temperatureName)).GetSpeciesShadowField(spc);
                var rule = scheme.Compile(solver.GridData, order);

                double IdV = temperature.LxError(exactTemperature[spc].Vectorize(time), (X, a, b) => (a - b).Pow2(), rule);
                L2Error += IdV;
                L2Error_Species.Add(spc, L2Error.Sqrt());

                solver.QueryHandler.ValueQuery("L2err_" + VariableNames.Temperature + "#" + spc, L2Error_Species[spc], true);
            }


            L2Error = L2Error.Sqrt();
            solver.QueryHandler.ValueQuery("L2err_" + VariableNames.Temperature, L2Error, true);
            return L2Error;            
        }

        public double ComputeEnergyError(Func<double, double> exactEnergy, double time) {
            int D = solver.GridData.SpatialDimension;

            int order = 0;
            if (solver.LsTrk.GetCachedOrders().Count > 0) {
                order = solver.LsTrk.GetCachedOrders().Max();
            } else {
                order = 1;
            }

            var SchemeHelper = solver.LsTrk.GetXDGSpaceMetrics(solver.LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;

            double TotalEnergy = 0;

            foreach (var spc in solver.LsTrk.SpeciesNames) {

                double c, rho;
                switch (spc) {
                    case "A": { c = solver.Control.ThermalParameters.c_A; rho = solver.Control.ThermalParameters.rho_A; break; }
                    case "B": { c = solver.Control.ThermalParameters.c_B; rho = solver.Control.ThermalParameters.rho_B; break; }
                    default: { throw new ArgumentException(); }
                }

                SpeciesId spId = solver.LsTrk.GetSpeciesId(spc);
                var scheme = SchemeHelper.GetVolumeQuadScheme(spId);

                string temperatureName = VariableNames.Temperature;
                ConventionalDGField temperature = ((XDGField)solver.CurrentStateVector.Mapping.Single(Field => Field.Identification == temperatureName)).GetSpeciesShadowField(spc);
                var rule = scheme.Compile(solver.GridData, order);

                double E = 0.0;
                CellQuadrature.GetQuadrature(new int[] { 1 }, solver.LsTrk.GridDat, rule,
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    temperature.Evaluate(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1,-1,0));
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        E += c * rho * ResultsOfIntegration[i, 0];
                }).Execute();
                TotalEnergy += E;
            }

            double EnergyError = (exactEnergy(time) - TotalEnergy).Abs() / exactEnergy(time).Abs(); // relative Error in Energy, scaled by some factor
            return EnergyError;
        }

        /// <summary>
        /// Computes the L2 Error of the actual solution against the exact solution in the control object 
        /// (<see cref="XNSE_Control.ExactSolutionTemperature"/>.
        /// </summary>
        public override double[] ComputeL2Error(double time, XNSE_Control control) {
            double[] Ret = new double[1];

            if (control.ExactSolutionTemperature != null) {
                double error = ComputeTemperatureError(control.ExactSolutionTemperature, time);
                Ret[0] = error;
            }
            
            return Ret;
        }
    }
}
