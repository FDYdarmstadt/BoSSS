using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.Tracing;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Diagnostics.Contracts;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Solution.XdgTimestepping;

namespace BoSSS.Solution.LevelSetTools.PhasefieldLevelSet
{
    /// <summary>
    /// Movement of the Phasefield
    /// </summary>
    partial class Phasefield : DgApplicationWithSolver<PhasefieldControl>
    {
        /// <summary>
        /// Move Phasefield
        /// </summary>
        /// <param name="_TimestepNo"></param>
        /// <param name="_dt"></param>
        /// <param name="_phystime"></param>
        public void MovePhasefield(int _TimestepNo,double _dt, double _phystime)
        {
            RunSolverOneStep(_TimestepNo, _dt, _phystime);
        }

        /// <summary>
        /// Initial step to initialize level set near its equilibrium
        /// TODO make this to calculate initial steady state without convection
        /// This has been replaced by <see cref="ReInit(double, double)"/>
        /// </summary>
        /// <param name="_TimestepNo"></param>
        /// <param name="_dt"></param>
        /// <param name="_phystime"></param>
        public void RelaxationStep(int _TimestepNo = -1, double _dt = 10.0, double _phystime = 0.0)
        {
            foreach (var f in this.Velocity)
                f.Clear();
            RunSolverOneStep(_TimestepNo, _dt, _phystime);
        }

        protected override double RunSolverOneStep(int _TimestepNo, double _dt, double _phystime)
        {
            using (new FuncTrace())
            {
                if (_TimestepNo == -1.0)
                {
                    Console.WriteLine("Initializing Phasefield");
                }
                else
                {
                    Console.WriteLine("Moving Phasefield in timestep #{0}, t = {1}, dt = {2} ...", _TimestepNo, _phystime, _dt);
                }

                // Perform timestep
                // ================                

                //if (this.Control.CurvatureCorrectionType == PhasefieldControl.CurvatureCorrection.DirectCoupledOnce)
                //{
                //    phi0.Clear();
                //    phi0.Acc(1.0, phi);
                //    gradPhi0.Clear();
                //    gradPhi0.Gradient(1.0, phi0);

                //    VectorField<SinglePhaseField> filtgrad;
                //    CurvatureAlgorithmsForLevelSet.CurvatureDriver(
                //                    CurvatureAlgorithmsForLevelSet.SurfaceStressTensor_IsotropicMode.Curvature_Projected,
                //                    CurvatureAlgorithmsForLevelSet.FilterConfiguration.Phasefield,
                //                    this.DCurvature, out filtgrad, CorrectionLsTrk,
                //                    this.DCurvature.Basis.Degree * 2,
                //                    phi0);

                //}                

                //PlotCurrentState(_phystime, new Foundation.IO.TimestepNumber(new int[] { _TimestepNo , 0}), 2);
                base.Timestepping.Solve(_phystime, _dt);
                //PlotCurrentState(_phystime, new Foundation.IO.TimestepNumber(new int[] { _TimestepNo }), 2);

                // algebraic correction
                switch (this.Control.CorrectionType)
                {
                    case PhasefieldControl.Correction.Concentration:
                        ConservativityCorrection();
                        break;
                    case PhasefieldControl.Correction.Mass:
                        MassCorrection();
                        break;
                    case PhasefieldControl.Correction.None:
                    default:
                        break;
                }

                // update DG LevelSet
                DGLevSet.Clear();
                DGLevSet.Acc(1.0, phi);

                CorrectionLevSet.Clear();
                CorrectionLevSet.Acc(1.0, phi);
                this.CorrectionLsTrk.UpdateTracker(0.0);

                // return
                // ======
                WriteLogLine(_TimestepNo, _phystime + _dt);

                //PlotCurrentState(_phystime, new Foundation.IO.TimestepNumber(new int[] { _TimestepNo }), 2);

                Console.WriteLine("done moving Phasefield in timestep #{0}.", _TimestepNo);

                return _dt;
            }
        }



        /// <summary>
        /// initializes the format of the Log File
        /// careful causes an exception, as the filestream is not properly closed when the Phasefield is reinitialized
        /// </summary>
        /// <param name="sessionID"></param>
        public void InitLogFile(Guid sessionID)
        {
            if(this.MPIRank != 0)
                return;

            File.Create("Phasefield_Quantities.txt");
            return;
        }

        /// <summary>
        /// writes one line to the Log File
        /// </summary>
        public void WriteLogLine(TimestepNumber TimestepNo, double phystime)
        {


            double[] BQnts = ComputeBenchmarkQuantities();

            if (this.MPIRank != 0)
                return;

            using (var Log = new StreamWriter("Phasefield_Quantities.txt", true))
            {
                string line = String.Format($"{TimestepNo}\t{phystime}\t{BQnts[0]}\t{BQnts[1]}\t{BQnts[2]}\t{BQnts[3]}\t{BQnts[4]}\t{BQnts[5]}\t{BQnts[6]}");
                Log.WriteLine(line);
                Log.Flush();
            }

            return;
        }

        private void MassCorrection()
        {
            double[] Qnts_old = ComputeBenchmarkQuantities();            

            CorrectionLevSet.Clear();
            CorrectionLevSet.Acc(1.0, phi);
            this.CorrectionLsTrk.UpdateTracker(0.0);

            double[] Qnts = ComputeBenchmarkQuantities();

            double massDiff = Qnts_old[0] - Qnts[0];

            // we assume the current phasefield is close to the equilibrium tangenshyperbolicus form
            SinglePhaseField phiNew = new SinglePhaseField(phi.Basis);
            GridData GridDat = (GridData)(phi.GridDat);
            double mass_uc = Qnts[0];

            int i = 0;
            while (massDiff.Abs() > 1e-6)
            {
                // calculated for a cone, one could include the shape e.g. by using the circularity
                // correction guess
                double correction = Math.Sign(massDiff) * 1e-10;//-Qnts[1] / (4 * Qnts[0] * Math.PI) * massDiff * 1e-5;

                // take the correction guess and calculate a forward difference to approximate the derivative
                phiNew.ProjectField(
                    (ScalarFunctionEx)delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result)
                    { // ScalarFunction2
                        Debug.Assert(result.Dimension == 2);
                        Debug.Assert(Len == result.GetLength(0));
                        int K = result.GetLength(1); // number of nodes

                        // evaluate Phi
                        // -----------------------------
                        phi.Evaluate(j0, Len, NS, result);

                        // compute the pointwise values of the new level set
                        // -----------------------------

                        result.ApplyAll(x => 0.5 * Math.Log(Math.Max(1 + x, 1e-10) / Math.Max(1 - x, 1e-10)) * Math.Sqrt(2) * this.Control.cahn);
                        result.ApplyAll(x => Math.Tanh((x + correction) / (Math.Sqrt(2) * this.Control.cahn)));
                    }
                );

                // update LsTracker
                CorrectionLevSet.Clear();
                CorrectionLevSet.Acc(1.0, phiNew);
                this.CorrectionLsTrk.UpdateTracker(0.0);

                Qnts = ComputeBenchmarkQuantities();

                correction = -(massDiff) / ((Qnts_old[0] - Qnts[0] - massDiff) / (correction));

                double initial = massDiff;
                bool finished = false;
                int k = 0;
                //while (massDiff.Abs() - initial.Abs() >= 0.0 && step > 1e-12)
                while(!finished)
                {
                    double step = Math.Pow(0.5, k);
                    // compute and project 
                    // step one calculate distance field phiDist = 0.5 * log(Max(1+c, eps)/Max(1-c, eps)) * sqrt(2) * Cahn
                    // step two project the new phasefield phiNew = tanh((cDist + correction)/(sqrt(2) * Cahn))
                    // ===================
                    phiNew.ProjectField(
                        (ScalarFunctionEx)delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result)
                        { // ScalarFunction2
                        Debug.Assert(result.Dimension == 2);
                            Debug.Assert(Len == result.GetLength(0));
                            int K = result.GetLength(1); // number of nodes

                        // evaluate Phi
                        // -----------------------------
                        phi.Evaluate(j0, Len, NS, result);

                        // compute the pointwise values of the new level set
                        // -----------------------------

                        result.ApplyAll(x => 0.5 * Math.Log(Math.Max(1 + x, 1e-10) / Math.Max(1 - x, 1e-10)) * Math.Sqrt(2) * this.Control.cahn);
                        result.ApplyAll(x => Math.Tanh((x + correction * step) / (Math.Sqrt(2) * this.Control.cahn)));
                        }
                    );

                    // update LsTracker
                    CorrectionLevSet.Clear();
                    CorrectionLevSet.Acc(1.0, phiNew);
                    this.CorrectionLsTrk.UpdateTracker(0.0);

                    Qnts = ComputeBenchmarkQuantities();
                    massDiff = Qnts_old[0] - Qnts[0];

                    if (massDiff.Abs() < (1 - 1e-4 * step) * initial.Abs())
                    {
                        finished = true;

                        // update field
                        phi.Clear();
                        phi.Acc(1.0, phiNew);

                        // update LsTracker
                        CorrectionLevSet.Clear();
                        CorrectionLevSet.Acc(1.0, phi);
                        this.CorrectionLsTrk.UpdateTracker(0.0);
                        Console.WriteLine($"" +
                            $"converged with stepsize:  {step}, correction: {correction}\n" +
                            $"                     dM:  {massDiff}");
                        if (k > 0)
                            Console.WriteLine($"Finished Linesearch in {k} iterations");
                    } 
                    else if (Math.Abs(correction * step) < 1e-15)
                    {
                        // reset LsTracker
                        CorrectionLevSet.Clear();
                        CorrectionLevSet.Acc(1.0, phi);
                        this.CorrectionLsTrk.UpdateTracker(0.0);

                        Qnts = ComputeBenchmarkQuantities();
                        massDiff = Qnts_old[0] - Qnts[0];

                        Console.WriteLine($" Linesearch failed after {k} iterations");
                        goto Failed;
                    }
                    else
                    {
                        k++;
                    }                    
                }
                i++;                   
            }

        Failed:           

            Console.WriteLine($"Performed Mass Correction in {i} iteratins: \n" +
                $"\told mass:           {Qnts_old[0]:N4}\n" +
                $"\tuncorrected mass:   {mass_uc:N4}\n" +
                $"\tcorrected mass:     {Qnts[0]:N4}");

        }       

        private void ConservativityCorrection()
        {
            throw new NotImplementedException();
        }

        private double[] ComputeBenchmarkQuantities()
        {
            int order = 0;
            if (CorrectionLsTrk.GetCachedOrders().Count > 0)
            {
                order = CorrectionLsTrk.GetCachedOrders().Max();
            }
            else
            {
                order = 1;
            }
            var SchemeHelper = CorrectionLsTrk.GetXDGSpaceMetrics(CorrectionLsTrk.SpeciesIdS.ToArray(), order, 1).XQuadSchemeHelper;

            // area of bubble
            double area = 0.0;
            SpeciesId spcId = CorrectionLsTrk.SpeciesIdS[1];
            var vqs = SchemeHelper.GetVolumeQuadScheme(spcId);
            CellQuadrature.GetQuadrature(new int[] { 1 }, CorrectionLsTrk.GridDat,
                vqs.Compile(CorrectionLsTrk.GridDat, order),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        area += ResultsOfIntegration[i, 0];
                }
            ).Execute();
            area = area.MPISum();

            // surface
            double surface = 0.0;
            //CellQuadratureScheme cqs = SchemeHelper.GetLevelSetquadScheme(0, LsTrk.Regions.GetCutCellMask());
            var surfElemVol = SchemeHelper.Get_SurfaceElement_VolumeQuadScheme(spcId, 0);
            CellQuadrature.GetQuadrature(new int[] { 1 }, CorrectionLsTrk.GridDat,
                surfElemVol.Compile(CorrectionLsTrk.GridDat, this.m_HMForder),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        surface += ResultsOfIntegration[i, 0];
                }
            ).Execute();
            surface = surface.MPISum();

            // circularity
            double diamtr_c = Math.Sqrt(4 * area / Math.PI);
            double perimtr_b = surface;

            double circ = Math.PI * diamtr_c / perimtr_b;

            // total concentration, careful above values are "old" when CorrectionTracker is not updated, this value is always "new"
            double concentration = 0.0;
            var tqs = new CellQuadratureScheme();
            CellQuadrature.GetQuadrature(new int[] { 1 }, phi.GridDat,
                tqs.Compile(phi.GridDat, phi.Basis.Degree * 2 + 2),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    phi.Evaluate(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        concentration += ResultsOfIntegration[i, 0];
                }
            ).Execute();
            concentration = concentration.MPISum();

            // total mixing energy
            double energy = 0.0;
            var eqs = new CellQuadratureScheme();

            int D = phi.GridDat.SpatialDimension;
            SinglePhaseField[] PhiGrad = new SinglePhaseField[D];
            for (int d = 0; d < D; d++)
            {
                PhiGrad[d] = new SinglePhaseField(phi.Basis, string.Format("G_{0}", d));
                PhiGrad[d].Derivative(1.0, phi, d);
            }

            MultidimensionalArray _Phi = new MultidimensionalArray(2);
            MultidimensionalArray _GradPhi = new MultidimensionalArray(3);
            MultidimensionalArray _NormGrad = new MultidimensionalArray(2);

            CellQuadrature.GetQuadrature(new int[] { 1 }, phi.GridDat,
                eqs.Compile(phi.GridDat, phi.Basis.Degree * 2 + 2),
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    int K = EvalResult.GetLength(1);
                    // alloc buffers
                    // -------------

                    if (_Phi.GetLength(0) != Length || _Phi.GetLength(1) != K)
                    {
                        _Phi.Allocate(Length, K);
                        _GradPhi.Allocate(Length, K, D);
                        _NormGrad.Allocate(Length, K);
                    }
                    else
                    {
                        _Phi.Clear();
                        _GradPhi.Clear();
                        _NormGrad.Clear();
                    }

                    // chemical potential
                    phi.Evaluate(i0, Length, QR.Nodes, _Phi.ExtractSubArrayShallow(-1, -1));
                    _Phi.ApplyAll(x => 0.25 / (this.Control.cahn.Pow2()) * (x.Pow2() - 1.0).Pow2());

                    for (int d = 0; d < D; d++)
                    {
                        PhiGrad[d].Evaluate(i0, Length, QR.Nodes, _GradPhi.ExtractSubArrayShallow(-1, -1, d));
                    }

                    // free surface energy
                    for (int d = 0; d < D; d++)
                    {
                        var GradPhi_d = _GradPhi.ExtractSubArrayShallow(-1, -1, d);
                        _NormGrad.Multiply(1.0, GradPhi_d, GradPhi_d, 1.0, "ik", "ik", "ik");
                    }
                    _NormGrad.ApplyAll(x => 0.5 * x);

                    EvalResult.ExtractSubArrayShallow(-1, -1, 0).Acc(1.0, _Phi);
                    EvalResult.ExtractSubArrayShallow(-1, -1, 0).Acc(1.0, _NormGrad);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++)
                        energy += ResultsOfIntegration[i, 0];
                }
            ).Execute();
            energy = energy.MPISum();

            // replace 1.96 (RB_TC2) with actual surface tension
            // see Yue (2004)
            double lambda = this.Control.cahn.Pow2() * 1.96 * surface / energy;            

            return new double[] { area, surface, circ, concentration, energy, lambda, this.Control.diff };
        }
    }
}
