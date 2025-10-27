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
using ilPSP.Utils;

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


        /// <summary>
        /// Identifies for which species an exact solution, name starting with <paramref name="SolPrefix"/>, is provided by the control object
        /// </summary>
        /// <param name="SolPrefix">Beginning (prefix) of the variable name to look at, e.g., 'Velocity'</param>
        /// <param name="Control"></param>
        /// <returns>
        /// - null, if no field with name starting with <paramref name="SolPrefix"/> can be found in <see cref="AppControl.ExactSolutions"/>
        /// - otherwise, a list of species names for which this solution is specified
        /// </returns>
        protected IEnumerable<string> GetSpecies(string SolPrefix, XNSE_Control Control) {
            var species = new HashSet<string>();

            bool bfound = false;
            foreach(string exSolName in Control.ExactSolutions_Evaluators_TimeDep.Keys) {
                if(exSolName.StartsWith(SolPrefix)) {
                    bfound = true;

                    var Split_exSolName = exSolName.Split('#', StringSplitOptions.RemoveEmptyEntries);
                    if(Split_exSolName.Length > 1)
                        species.Add(Split_exSolName[1]);
                }
            }


            if(bfound == false)
                return null; // no exact solution specified -> indicate this with null

            if(species.Count == 0) {
                species.AddRange(solver.LsTrk.SpeciesNames);
            }

            return species;
        }


        protected IDictionary<string, Func<double[], double, double>> GetExactSolution(XNSE_Control control, string VariableName) {

            var FluidSpecies = GetSpecies(VariableName, control);
            if(FluidSpecies == null)
                return null;

            var ret = new Dictionary<string, Func<double[], double, double>>();

            foreach(var s in FluidSpecies) {


                Func<double[], double, double> exSolImpl = null;

                if(control.ExactSolutions_Evaluators_TimeDep.TryGetValue($"{VariableName}#{s}", out exSolImpl)) {
                    // found species-specific solution
                } else if(control.ExactSolutions_Evaluators_TimeDep.TryGetValue(VariableName, out exSolImpl)) {
                    // found common solution for both species
                } else {
                    // default value
                    exSolImpl = (X, t) => 0.0;
                }

                ret.Add(s, exSolImpl);
            }


            return ret;

        }

        /// <summary>
        /// Returns the longest common prefix of the given strings.
        /// Returns empty string if input is null/empty or contains only null/empty items.
        /// </summary>
        static string GetCommonPrefix(IEnumerable<string> items, bool ignoreCase = false) {
            if(items == null)
                return string.Empty;

            var arr = items.Where(s => !string.IsNullOrEmpty(s)).ToArray();
            if(arr.Length == 0)
                return string.Empty;
            if(arr.Length == 1)
                return arr[0];

            int minLen = arr.Min(s => s.Length);
            int i = 0;
            for(; i < minLen; i++) {
                char c = arr[0][i];
                for(int j = 1; j < arr.Length; j++) {
                    char cj = arr[j][i];
                    if(ignoreCase) {
                        if(char.ToUpperInvariant(c) != char.ToUpperInvariant(cj)) {
                            goto end;
                        }
                    } else {
                        if(c != cj) {
                            goto end;
                        }
                    }
                }
            }
        end:
            return arr[0].Substring(0, i);
        }

        protected IDictionary<string, Func<double[], double, double>[]> GetExactSolution(XNSE_Control control, string[] VectorVariableNames) {
            int D = VectorVariableNames.Length;

            var FluidSpecies = GetSpecies(GetCommonPrefix(VectorVariableNames), control);
            if(FluidSpecies == null)
                return null;

            var ret = new Dictionary<string, Func<double[], double, double>[]>();

            foreach(var s in FluidSpecies) {
                var VelVec = new Func<double[], double, double>[D];
                for(int d = 0; d < D; d++) {

                    Func<double[], double, double> exSolImpl = null;

                    if(control.ExactSolutions_Evaluators_TimeDep.TryGetValue($"{VectorVariableNames[d]}#{s}", out exSolImpl)) {
                        // found species-specific solution
                    } else if(control.ExactSolutions_Evaluators_TimeDep.TryGetValue(VectorVariableNames[d], out exSolImpl)) {
                        // found common solution for both species
                    } else {
                        // default value
                        exSolImpl = (X, t) => 0.0;
                    }

                    VelVec[d] = exSolImpl;
                }
                ret.Add(s, VelVec);
            }

            return ret;
        }


    }

    /// <summary>
    /// error evaluator for XNSE-based tests 
    /// computes error for: velocity, pressure 
    /// </summary>
    /// <remarks>
    /// Seems redundant with <see cref="L2ErrorLogger"/>
    /// </remarks>
    public class XNSEErrorEvaluator<T> : ErrorEvaluator where T : XNSE_Control, new() {

        public XNSEErrorEvaluator(IApplication solver) : base(solver) {

        }

        new protected ApplicationWithSolver<T> solver {
            get {
                return (ApplicationWithSolver<T>)(base.solver);
            }
        }

        public double[] ComputeVelocityError(IDictionary<string, Func<double[], double, double>[]> exactVelocity, double time) {
            int D = solver.GridData.SpatialDimension;
            var FluidSpecies = exactVelocity.Keys.ToArray();
            double[] Ret = new double[D];

            int order = 0;
            if(solver.LsTrk.GetCachedOrders().Count > 0) {
                order = solver.LsTrk.GetCachedOrders().Max();
            } else {
                order = 1;
            }

            var SchemeHelper = solver.LsTrk.GetXDGSpaceMetrics(solver.LsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;

            Dictionary<string, double[]> L2Error_Species = new Dictionary<string, double[]>();
            double[] L2Error = new double[D];

            foreach(var spc in FluidSpecies) {
                L2Error_Species.Add(spc, new double[D]);

                SpeciesId spId = solver.LsTrk.GetSpeciesId(spc);
                var scheme = SchemeHelper.GetVolumeQuadScheme(spId);


                for(int d = 0; d < D; d++) {
                    string velocity = VariableNames.VelocityVector(D)[d];
                    ConventionalDGField Vel_d = ((XDGField)solver.CurrentStateVector.Mapping.Single(Field => Field.Identification == velocity)).GetSpeciesShadowField(spc);

                    double LxError = Vel_d.LxError(exactVelocity[spc][d].Vectorize(time), null, scheme.SaveCompile(Vel_d.GridDat, order));
                    LxError = (LxError > -1e-12 && LxError < 0) ? 0.0 : LxError; // Avoid nans if solution is too close to the analytic one
                    L2Error_Species[spc][d] = LxError.Sqrt();

                    L2Error[d] += L2Error_Species[spc][d].Pow2();

                    solver.QueryHandler.ValueQuery("L2err_" + VariableNames.Velocity_d(d) + "#" + spc, L2Error_Species[spc][d], true);
                }
            }
            L2Error = L2Error.Select(x => x.Sqrt()).ToArray();

            for(int d = 0; d < D; d++) {
                solver.QueryHandler.ValueQuery("L2err_" + VariableNames.Velocity_d(d), L2Error[d], true);
                Ret[d] = L2Error[d];
            }
            return L2Error;
        }

        public double ComputePressureError(IDictionary<string, Func<double[], double, double>> exactPressure, double time) {
            int D = solver.GridData.SpatialDimension;
            var FluidSpecies = exactPressure.Keys.ToArray();

            string pressureName = VariableNames.Pressure;
            XDGField pressure = (XDGField)solver.CurrentStateVector.Mapping.Single(Field => Field.Identification == pressureName);

            int order = 0;
            if(solver.LsTrk.GetCachedOrders().Count > 0) {
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
            foreach(var spc in FluidSpecies) {

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

            foreach(var spc in FluidSpecies) {

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


        IDictionary<string, Func<double[], double, double>[]> GetExactSolutionVelocity(XNSE_Control control) {
            int D = solver.GridData.SpatialDimension;
            return base.GetExactSolution(control, VariableNames.VelocityVector(D));
        }

        IDictionary<string, Func<double[], double, double>> GetExactSolutionPressure(XNSE_Control control) {
            return base.GetExactSolution(control, VariableNames.Pressure);
        }



        /// <summary>
        /// Computes the L2 Error of the actual solution against the exact solution in the control object 
        /// (<see cref="AppControl.ExactSolutions"/>).
        /// </summary>
        public override double[] ComputeL2Error(double time, XNSE_Control control) {
            int D = solver.GridData.SpatialDimension;
            double[] Ret = new double[D + 1];

            IDictionary<string, Func<double[], double, double>[]> ExactSolutionVelocity = this.GetExactSolutionVelocity(control);
            IDictionary<string, Func<double[], double, double>> ExactSolutionPressure = this.GetExactSolutionPressure(control);

            if(ExactSolutionVelocity != null) {
                double[] error = ComputeVelocityError(ExactSolutionVelocity, time);
                error.CopyTo(Ret, 0);
            }
            if(ExactSolutionPressure != null) {
                double error = ComputePressureError(ExactSolutionPressure, time);
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
            Func<double[], double, double> phiExactFunc;// = solver.Control.Phi;
            solver.Control.ExactSolutions_Evaluators_TimeDep.TryGetValue("Phi", out phiExactFunc);
            exactPhi = new LevelSet(new Basis(solver.GridData, 8), "exactLevelSet");
            exactPhi.Clear();
            exactPhi.ProjectField(NonVectorizedScalarFunction.Vectorize(phiExactFunc, time));
            // exact level-set tracker
            exactLsTrk = new LevelSetTracker((GridData)solver.GridData,
            CutCellQuadratureMethod.Saye, 1, solver.LsTrk.SpeciesNames.ToArray(), exactPhi);
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
        /// computes the error of the interface points with respect to the given interface form
        /// </summary>
        /// <returns></returns>
        public double ComputeInterfacePointsError(double time) {

            SinglePhaseField PhiCG = solver.LsUpdater.LevelSets[VariableNames.LevelSetCG].CGLevelSet;
            SubGrid sbgrd = solver.LsTrk.Regions.GetCutCellSubGrid();
            MultidimensionalArray interfaceP = XNSEUtils.GetInterfacePoints(solver.LsTrk, PhiCG, sbgrd);

            Func<double[], double, double> phiExact;// = solver.Control.Phi;
            solver.Control.ExactSolutions_Evaluators_TimeDep.TryGetValue("Phi", out phiExact);
            double error = 0.0;
            for (int i = 0; i < interfaceP.Lengths[0]; i++) {
                double dist = phiExact(interfaceP.GetRow(i), time);
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

            double[] Ret = new double[2]; //4 + spcIds.Count];

            if (control.ExactSolutions.ContainsKey("Phi")) {
                //Ret[0] = ComputeInterfacePointsError(time);
                Ret[0] = ComputeLevelSetError(cm);
                Ret[1] = ComputeDGLevelSetError(cm);
            }
            //Ret[2] = ComputeDGLevelSetGradientError(cm);
            //Ret[3] = ComputeInterfaceSizeError();

            //Dictionary<SpeciesId, double> spcArea = ComputeSpeciesDomainSizeError();
            //int n = 4;
            //foreach (SpeciesId spc in spcIds) {
            //    Ret[n] = spcArea[spc];
            //    n++;
            //}

            return Ret;
        }

    }   
}
