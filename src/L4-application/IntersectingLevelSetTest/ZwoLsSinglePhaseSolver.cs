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
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using NUnit.Framework;
using System;

namespace IntersectingLevelSetTest {

    /// <summary>
    /// This guy tests the basic functionality of the XDG framework
    /// if more than one level-set is involved.
    /// </summary>
    internal class ZwoLsSinglePhaseSolver<T> : BoSSS.Solution.Application<T> where T : BoSSS.Solution.Control.AppControl, new() {
        internal XQuadFactoryHelper.MomentFittingVariants MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;

        protected override IGrid CreateOrLoadGrid() {
            var t = Triangle.Instance;
            var grd = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-0.5, 0.5, 6), GenericBlas.Linspace(-0.5, 0.5, 6));
            return grd;
        }

        public override void Init(BoSSS.Solution.Control.AppControl control) {
            BoSSS.Solution.Application.DeleteOldPlotFiles();
            base.Init(control);
        }

        private LevelSet Phi0;

        private LevelSet Phi1;

        private SinglePhaseField u;

        private SinglePhaseField du_dx;

        private SinglePhaseField du_dx_Exact;

        private SinglePhaseField ERR;

        private SinglePhaseField Amarker;
        private SinglePhaseField Bmarker;
        private SinglePhaseField Cmarker;

        protected override void CreateFields() {
            Phi0 = new LevelSet(new Basis(this.GridData, 3), "Phi_0");
            Phi1 = new LevelSet(new Basis(this.GridData, 3), "Phi_1");

            string[,] speciesTable = new string[2, 2];
            speciesTable[0, 0] = "B"; // Liquid
            speciesTable[0, 1] = "C"; // Solid
            speciesTable[1, 0] = "A"; // Gas
            speciesTable[1, 1] = "C"; // Solid

            base.LsTrk = new LevelSetTracker((GridData)(this.GridData), MomentFittingVariant, 1, speciesTable, Phi0, Phi1);

            u = new SinglePhaseField(new Basis(GridData, DEGREE), "U");
            du_dx = new SinglePhaseField(new Basis(GridData, DEGREE), "du_dx");
            du_dx_Exact = new SinglePhaseField(new Basis(GridData, DEGREE), "du_dx_exact");
            ERR = new SinglePhaseField(new Basis(GridData, DEGREE), "ERROR");

            Amarker = new SinglePhaseField(new Basis(this.GridData, 0), "Amarker");
            Bmarker = new SinglePhaseField(new Basis(this.GridData, 0), "Bmarker");
            Cmarker = new SinglePhaseField(new Basis(this.GridData, 0), "Cmarker");

            base.m_RegisteredFields.Add(Phi0);
            base.m_RegisteredFields.Add(Phi1);
            base.m_RegisteredFields.Add(u);
            base.m_RegisteredFields.Add(du_dx);
            base.m_RegisteredFields.Add(du_dx_Exact);
            base.m_RegisteredFields.Add(ERR);
            base.m_RegisteredFields.Add(Amarker);
            base.m_RegisteredFields.Add(Bmarker);
            base.m_RegisteredFields.Add(Cmarker);
        }

        /// <summary>
        /// DG polynomial degree
        /// </summary>
        internal int DEGREE = 1;

        private void LsUpdate(double t) {
            Console.WriteLine("LSUpdate t = " + t);
            SetLs(t);
            
            LsTrk.UpdateTracker(t);
            LsTrk.PushStacks();

            Amarker.Clear();
            Bmarker.Clear();
            Cmarker.Clear();
            Amarker.AccConstant(1.0, LsTrk.Regions.GetSpeciesSubGrid("A").VolumeMask);
            Bmarker.AccConstant(1.0, LsTrk.Regions.GetSpeciesSubGrid("B").VolumeMask);
            Cmarker.AccConstant(1.0, LsTrk.Regions.GetSpeciesSubGrid("C").VolumeMask);
        }

        private void SetLs(double t) {
            t = t / 90 * Math.PI;
            double phi0(double x, double y) {
                return (x - (Math.Tan(t) * y));
            }
            double phi1(double x, double y) {
                return (Math.Tan(t) * x) - (y);
            }
            
            Phi0.ProjectField(phi0);
            Phi1.ProjectField(phi1);
        }

        private void SetLs2(double t) {

            double phi0(double x, double y) {
                return - Math.Sqrt((x- 1 - t * 0.1).Pow2() + y.Pow2()) + 1.2 + t * 0.1;
            }
            double phi1(double x, double y) {
                return (x - 0.1);
            }

            Phi0.ProjectField(phi0);
            Phi1.ProjectField(phi1);
        }

        private void SetLs3(double t) {

            double phi0(double x, double y) {
                return (x + 0.1 + t * 0.01);
            }
            double phi1(double x, double y) {
                return (x - 0.1);
            }

            Phi0.ProjectField(phi0);
            Phi1.ProjectField(phi1);
        }

        protected override void SetInitial(double t) {
            this.LsUpdate(t);

            u.ProjectField((x, y) => x );
            du_dx_Exact.ProjectField((x, y) => 1);
        }

        private XDifferentialOperatorMk2 Op;

        private int QuadOrder {
            get {
                return this.DEGREE * 2 + 2;
            }
        }

        protected override void CreateEquationsAndSolvers(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L) {
            Op = new XDifferentialOperatorMk2(1, 0, 1,
                QuadOrderFunc: (int[] DomDegs, int[] ParamDegs, int[] CoDomDegs) => QuadOrder,
                __Species: new[] { "A", "B"},
                __varnames: new[] { "u", "c1" });

            Op.EquationComponents["c1"].Add(new LevelSetJumpFlux()); // Flux in Bulk Phase;
            Op.EquationComponents["c1"].Add(new LevSetJump_AB()); // flux am lev-set 0
            Op.EquationComponents["c1"].Add(new LevSetFlx_CA()); // flux am lev-set 1
            Op.EquationComponents["c1"].Add(new LevSetFlx_CB()); // flux am lev-set 1

            Op.Commit();
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            Console.WriteLine("    Timestep # " + TimestepNo + ", phystime = " + phystime);

            //phystime = 1.8;
            LsUpdate(phystime);

            // operator-matrix assemblieren
            MsrMatrix OperatorMatrix = new MsrMatrix(u.Mapping, u.Mapping);
            double[] Affine = new double[OperatorMatrix.RowPartitioning.LocalLength];

            // operator matrix assembly
            XDifferentialOperatorMk2.XEvaluatorLinear mtxBuilder = Op.GetMatrixBuilder(base.LsTrk, u.Mapping, null, u.Mapping);
            mtxBuilder.time = 0.0;
            mtxBuilder.ComputeMatrix(OperatorMatrix, Affine);
            OperatorMatrix.SaveToTextFile("operatorMatrix.txt");

            // mass matrix factory
            SpeciesId[] species = new SpeciesId[] { LsTrk.GetSpeciesId("A") , LsTrk.GetSpeciesId("B") };
            var Mfact = LsTrk.GetXDGSpaceMetrics(species, QuadOrder, 1).MassMatrixFactory;// new MassMatrixFactory(u.Basis, Agg);


            // Mass matrix/Inverse Mass matrix
            var Mass = Mfact.GetMassMatrix(u.Mapping, new double[] { 1.0 }, false, species);
            var MassInv = Mass.InvertBlocks(OnlyDiagonal: true, Subblocks: true, ignoreEmptyBlocks: true, SymmetricalInversion: false);

            // test that operator depends only on B-species values
            //double DepTest = LsTrk.Regions.GetSpeciesSubGrid("B").TestMatrixDependency(OperatorMatrix, u.Mapping, u.Mapping);
            //Console.WriteLine("Matrix dependency test: " + DepTest);
            //Assert.LessOrEqual(DepTest, 0.0);

            // operator auswerten:
            double[] x = new double[Affine.Length];
            BLAS.daxpy(x.Length, 1.0, Affine, 1, x, 1);
            OperatorMatrix.SpMVpara(1.0, u.CoordinateVector, 1.0, x);
            MassInv.SpMV(1.0, x, 0.0, du_dx.CoordinateVector);

            

            // compute integrals 
            Integrals integrals = new Integrals();
            //integrals.Evaluate(LsTrk, 2, LsTrk.GetSpeciesId("B"), LsTrk.GetSpeciesId("A"));


            // compute error
            ERR.Clear();
            ERR.Acc(1.0, du_dx_Exact);
            ERR.Acc(-1.0, du_dx);
            //ERR.GetSpeciesShadowField("C").Clear();
            double L2Err = ERR.L2Norm();
            Console.WriteLine("L2 Error: " + L2Err);

            // check error
            double ErrorThreshold = 1.0e-1;
            if (this.MomentFittingVariant == XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes)
                ErrorThreshold = 1.0e-6; // HMF is designed for such integrands and should perform close to machine accuracy; on general integrands, the precision is different.

            bool IsPassed = (L2Err <= ErrorThreshold);
            if (IsPassed) {
                Console.WriteLine("Test PASSED");
            } else {
                Console.WriteLine("Test FAILED: check errors.");
                //PlotCurrentState(phystime, TimestepNo, 3);
            }

            if (TimestepNo > 1) {
                //Assert.LessOrEqual(xL2Err, ErrorThreshold, "XDG L2 error of computing du_dx");
            }

            // return/Ende
            base.NoOfTimesteps = 20;
            //base.NoOfTimesteps = 2;
            dt = 1;
            return dt;
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            string filename = "ZwoLsTest." + timestepNo;
            Tecplot.PlotFields(new DGField[] { u, du_dx, Phi0, Phi1, Amarker, Bmarker, Cmarker, du_dx_Exact, ERR }, filename, 0, superSampling);
        }
    }
}