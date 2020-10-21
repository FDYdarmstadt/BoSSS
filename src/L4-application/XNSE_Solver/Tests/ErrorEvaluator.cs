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

namespace BoSSS.Application.XNSE_Solver.Tests
{
    class ErrorEvaluator
    {

        XNSE solver;

        public ErrorEvaluator(XNSE solver) 
        {
            this.solver = solver;
        }

        /// <summary>
        /// Computes the L2 Error of the actual solution against the exact solution in the control object 
        /// (<see cref="XNSE_Control.ExactSolutionVelocity"/> and <see cref="XNSE_Control.ExactSolutionPressure"/>).
        /// </summary>
        public double[] ComputeL2Error(double time, XNSE_Control control)
        {
            int D = solver.GridData.SpatialDimension;
            double[] Ret = new double[D + 1];

            if (control.ExactSolutionVelocity == null && control.ExactSolutionPressure == null)
                // nothing to do
                return Ret;

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


            // Velocity error
            // ==============
            if (control.ExactSolutionVelocity != null)
            {
                Dictionary<string, double[]> L2Error_Species = new Dictionary<string, double[]>();
                double[] L2Error = new double[D];

                foreach (var spc in solver.LsTrk.SpeciesNames)
                {
                    L2Error_Species.Add(spc, new double[D]);

                    SpeciesId spId = solver.LsTrk.GetSpeciesId(spc);
                    var scheme = SchemeHelper.GetVolumeQuadScheme(spId);


                    for (int d = 0; d < D; d++)
                    {
                        string velocity = VariableNames.VelocityVector(D)[d];
                        ConventionalDGField Vel_d = ((XDGField)solver.CurrentStateVector.Mapping.Single(Field => Field.Identification == velocity)).GetSpeciesShadowField(spc);

                        L2Error_Species[spc][d] = Vel_d.L2Error(control.ExactSolutionVelocity[spc][d].Vectorize(time), order, scheme);
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

            }


            // pressure error
            // ==============
            if (control.ExactSolutionPressure != null)
            {

                // pass 1: mean value of pressure difference
                double DiffInt = 0;
                foreach (var spc in solver.LsTrk.SpeciesNames)
                {

                    SpeciesId spId = solver.LsTrk.GetSpeciesId(spc);
                    var scheme = SchemeHelper.GetVolumeQuadScheme(spId);
                    var rule = scheme.Compile(solver.GridData, order);

                    string pressureName = VariableNames.Pressure;
                    XDGField pressure = (XDGField)solver.CurrentStateVector.Mapping.Single(Field => Field.Identification == pressureName);
                    DiffInt += pressure.GetSpeciesShadowField(spc).LxError(control.ExactSolutionPressure[spc].Vectorize(time), (X, a, b) => (a - b), rule);
                }
                double Volume2 = (new SubGrid(CellMask.GetFullMask(solver.GridData))).Volume;
                double PressureDiffMean = DiffInt / Volume2;


                double L2Error = 0;
                Dictionary<string, double> L2Error_Species = new Dictionary<string, double>();

                foreach (var spc in solver.LsTrk.SpeciesNames)
                {

                    SpeciesId spId = solver.LsTrk.GetSpeciesId(spc);
                    var scheme = SchemeHelper.GetVolumeQuadScheme(spId);
                    var rule = scheme.Compile(solver.GridData, order);

                    string pressureName = VariableNames.Pressure;
                    XDGField pressure = (XDGField)solver.CurrentStateVector.Mapping.Single(Field => Field.Identification == pressureName);
                    double IdV = pressure.GetSpeciesShadowField(spc).LxError(control.ExactSolutionPressure[spc].Vectorize(time), (X, a, b) => (a - b - PressureDiffMean).Pow2(), rule);
                    L2Error += IdV;
                    L2Error_Species.Add(spc, IdV.Sqrt());

                    solver.QueryHandler.ValueQuery("L2err_" + VariableNames.Pressure + "#" + spc, L2Error_Species[spc], true);
                }


                L2Error = L2Error.Sqrt();
                solver.QueryHandler.ValueQuery("L2err_" + VariableNames.Pressure, L2Error, true);
                Ret[D] = L2Error;

            } //*/

            return Ret;
        }

    }
}
