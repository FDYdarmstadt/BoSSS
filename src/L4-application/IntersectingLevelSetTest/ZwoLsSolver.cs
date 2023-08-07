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
using System.Runtime.Serialization;

namespace IntersectingLevelSetTest {

    /// <summary>
    /// This guy tests the basic functionality of the XDG framework
    /// if more than one level-set is involved.
    /// </summary>
    internal class ZwoLsSolver<T> : BoSSS.Solution.Application<T> where T : BoSSS.Solution.Control.AppControl, new() {
        internal XQuadFactoryHelper.MomentFittingVariants MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes;

        int resolution;
        int dimension;
        double errorThreshold;
        public double errorBack;
        Func<double, double, double, double> levelSet0;
        Func<double, double, double, double> levelSet1;
        Func<double, double, double, double, double> levelSet3D_0;
        Func<double, double, double, double, double> levelSet3D_1;

        protected override IGrid CreateOrLoadGrid() {
            var t = Triangle.Instance;
            if (dimension == 2) {
                var grd = Grid2D.Cartesian2DGrid(GenericBlas.Linspace(-0.5, 0.5, resolution), GenericBlas.Linspace(-0.5, 0.5, resolution));
                return grd;
            }
            else if (dimension == 3) {
                var grd = Grid3D.Cartesian3DGrid(GenericBlas.Linspace(-0.5, 0.5, resolution), GenericBlas.Linspace(-0.5, 0.5, resolution), GenericBlas.Linspace(-0.5, 0.5, resolution));
                return grd;
            }
            else {
                throw new ArgumentException("please give spatial dimension value 2 or 3");
            }
        }

        public override void Init(BoSSS.Solution.Control.AppControl control) {
            //BoSSS.Solution.Application.DeleteOldPlotFiles();
            base.Init(control);
            if (control is TestControl testControl) {
                levelSet0 = testControl.LevelSet0;
                levelSet1 = testControl.LevelSet1;
                levelSet3D_0 = testControl.LevelSet3D_0;
                levelSet3D_1 = testControl.LevelSet3D_1;
                dimension = testControl.Dimension;
                resolution = testControl.Resolution;
                errorThreshold = testControl.ErrorThreshold;
            }

        }

        private LevelSet Phi0;

        private LevelSet Phi1;

        private SinglePhaseField u;

        private SinglePhaseField du_dx;

        private SinglePhaseField du_dx_Exact;

        private SinglePhaseField ERR;

        private XDGField XERR;

        private SinglePhaseField Amarker;
        private SinglePhaseField Bmarker;
        private SinglePhaseField Cmarker;

        protected override void CreateFields() {
            //Phi0 = new LevelSet(new Basis(this.GridData, 3), "Phi_0");
            //Phi1 = new LevelSet(new Basis(this.GridData, 3), "Phi_1");
            Phi0 = new LevelSet(new Basis(this.GridData, DEGREE), "Phi_0");
            Phi1 = new LevelSet(new Basis(this.GridData, DEGREE), "Phi_1");

            string[,] speciesTable = new string[2, 2];
            speciesTable[0, 0] = "B"; // Liquid
            speciesTable[0, 1] = "C"; // Solid
            speciesTable[1, 0] = "A"; // Gas
            speciesTable[1, 1] = "C"; // Solid

            base.LsTrk = new LevelSetTracker((GridData)(this.GridData), MomentFittingVariant, 1, speciesTable, Phi0, Phi1);

            u = new SinglePhaseField(new Basis(this.GridData, DEGREE), "U");
            du_dx = new SinglePhaseField(new Basis(this.GridData, DEGREE), "du_dx");
            du_dx_Exact = new SinglePhaseField(new Basis(this.GridData, DEGREE), "du_dx_exact");
            ERR = new SinglePhaseField(new Basis(this.GridData, DEGREE), "ERROR");
            XERR = new XDGField(new XDGBasis(LsTrk, DEGREE), "ERROR_xdg");

            Amarker = new SinglePhaseField(new Basis(this.GridData, 0), "Amarker");
            Bmarker = new SinglePhaseField(new Basis(this.GridData, 0), "Bmarker");
            Cmarker = new SinglePhaseField(new Basis(this.GridData, 0), "Cmarker");

            base.m_RegisteredFields.Add(Phi0);
            base.m_RegisteredFields.Add(Phi1);
            base.m_RegisteredFields.Add(u);
            base.m_RegisteredFields.Add(du_dx);
            base.m_RegisteredFields.Add(du_dx_Exact);
            base.m_RegisteredFields.Add(ERR);
            base.m_RegisteredFields.Add(XERR);
            base.m_RegisteredFields.Add(Amarker);
            base.m_RegisteredFields.Add(Bmarker);
            base.m_RegisteredFields.Add(Cmarker);
        }

        /// <summary>
        /// DG polynomial degree
        /// </summary>
        internal int DEGREE;

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

            if (dimension == 2) {

                double phi0(double x, double y) {
                    return levelSet0(x, y, t);
                }

                double phi1(double x, double y) {
                    return levelSet1(x, y, t);
                }

                Phi0.ProjectField(phi0);
                Phi1.ProjectField(phi1);
            }
            else if (dimension == 3) {

                double phi0(double x, double y, double z) {
                    return levelSet3D_0(x, y, z, t);
                }

                double phi1(double x, double y, double z) {
                    return levelSet3D_1(x, y, z, t);
                }

                Phi0.ProjectField(phi0);
                Phi1.ProjectField(phi1);
            }

        }


        protected override void SetInitial(double t) {
            this.LsUpdate(t);

            if (dimension == 2) {
                //u.ProjectField((x, y) => x * x);
                //du_dx_Exact.ProjectField((x, y) => 2 * x);
                //u.ProjectField((x, y) => y * y);
                //du_dx_Exact.ProjectField((x, y) => 0);
                //u.ProjectField((x, y) => Math.Sin(x));
                //du_dx_Exact.ProjectField((x, y) => Math.Cos(x));
                u.ProjectField((x, y) => Math.Sin(x) * Math.Cos(y));
                du_dx_Exact.ProjectField((x, y) => Math.Cos(x) * Math.Cos(y));
                //u.ProjectField((x, y) => x * x * y);
                //du_dx_Exact.ProjectField((x, y) => 2 * x * y);
                //u.ProjectField((x, y) => 1);
                //du_dx_Exact.ProjectField((x, y) => 0);
            }
            else if (dimension == 3) {
                u.ProjectField((x, y, z) => x * x);
                du_dx_Exact.ProjectField((x, y, z) => 2 * x);
            }
        }

        private XSpatialOperatorMk2 Op;

        private int QuadOrder {
            get {
                return this.DEGREE * 2 + Math.Max(Phi0.Basis.Degree,Phi1.Basis.Degree); 
                // degree of ansatz function + degree of test function + Max level set degree. 
            }
        }

        protected override void CreateEquationsAndSolvers(BoSSS.Solution.LoadBalancing.GridUpdateDataVaultBase L) {
            Op = new XSpatialOperatorMk2(1, 0, 1,
                QuadOrderFunc: (int[] DomDegs, int[] ParamDegs, int[] CoDomDegs) => QuadOrder,
                __Species: new[] { "B" },
                __varnames: new[] { "u", "c1" });

            Op.EquationComponents["c1"].Add(new DxFlux()); // Flux in Bulk Phase;
            Op.EquationComponents["c1"].Add(new LevSetFlx_AB()); // flux am lev-set 0
            Op.EquationComponents["c1"].Add(new LevSetFlx_CB()); // flux am lev-set 1
            Op.EquationComponents["c1"].Add(new LevSetFlx_CA()); // flux am lev-set 1

            //Op.EquationComponents["c1"].Add(new DxBroken());

            Op.Commit();
        }

        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            Console.WriteLine("    Timestep # " + TimestepNo + ", phystime = " + phystime);


            //reset current solution

            //phystime = 1.8;
            LsUpdate(phystime);

            //var map = du_dx.Mapping;
            var map = u.Mapping;

            // operator-matrix assemblieren
            MsrMatrix OperatorMatrix = new MsrMatrix(map, map);
            double[] Affine = new double[OperatorMatrix.RowPartitioning.LocalLength];

            // operator matrix assembly
            XSpatialOperatorMk2.XEvaluatorLinear mtxBuilder = Op.GetMatrixBuilder(base.LsTrk, map, null, map);
            mtxBuilder.time = 0.0;
            mtxBuilder.ComputeMatrix(OperatorMatrix, Affine);

            // mass matrix factory
            var Mfact = LsTrk.GetXDGSpaceMetrics(new SpeciesId[] { LsTrk.GetSpeciesId("B") }, QuadOrder, 1).MassMatrixFactory;// new MassMatrixFactory(u.Basis, Agg);

            // Mass matrix/Inverse Mass matrix
            var Mass = Mfact.GetMassMatrix(map, new double[] { 1.0 }, false, LsTrk.GetSpeciesId("B"));
            var MassInv = Mass.InvertBlocks(OnlyDiagonal: true, Subblocks: true, ignoreEmptyBlocks: true, SymmetricalInversion: false);

            // test that operator depends only on B-species values
            double DepTest = LsTrk.Regions.GetSpeciesSubGrid("B").TestMatrixDependency(OperatorMatrix, map, map);
            Console.WriteLine("Matrix dependency test: " + DepTest);
            Assert.LessOrEqual(DepTest, 0.0);

            // operator auswerten:
            double[] x = new double[Affine.Length];
            BLAS.daxpy(x.Length, 1.0, Affine, 1, x, 1);
            OperatorMatrix.SpMVpara(1.0, u.CoordinateVector, 1.0, x);
            du_dx.Clear();
            MassInv.SpMV(1.0, x, 0.0, du_dx.CoordinateVector);

            OperatorMatrix.SaveToTextFile("matrix.txt");

            // compute integrals 
            Integrals integrals = new Integrals();
            //integrals.Evaluate(LsTrk, 2, LsTrk.GetSpeciesId("B"), LsTrk.GetSpeciesId("A"));


            // compute error
            ERR.Clear();
            ERR.Acc(1.0, du_dx_Exact, LsTrk.Regions.GetSpeciesSubGrid("B").VolumeMask);
            ERR.Acc(-1.0, du_dx, LsTrk.Regions.GetSpeciesSubGrid("B").VolumeMask);
            double L2Err = ERR.L2Norm(LsTrk.Regions.GetSpeciesSubGrid("B").VolumeMask);
            Console.WriteLine("L2 Error: " + L2Err);

            XERR.Clear();
            XERR.GetSpeciesShadowField("B").Acc(1.0, ERR, LsTrk.Regions.GetSpeciesSubGrid("B").VolumeMask);
            double xL2Err = XERR.L2Norm();
            Console.WriteLine("L2 Error (in XDG space): " + xL2Err);

            // check error
            //double ErrorThreshold = 1.0e-1;
            double ErrorThreshold = errorThreshold;
            if (this.MomentFittingVariant == XQuadFactoryHelper.MomentFittingVariants.OneStepGaussAndStokes)
                ErrorThreshold = 1.0e-6; // HMF is designed for such integrands and should perform close to machine accuracy; on general integrands, the precision is different.

            bool IsPassed = (L2Err <= ErrorThreshold || xL2Err <= ErrorThreshold);
            if (IsPassed) {
                Console.WriteLine("Test PASSED");
                //PlotCurrentState(phystime, TimestepNo, 4);
            }
            else {
                Console.WriteLine("Test FAILED: check errors.");
                //PlotCurrentState(phystime, TimestepNo, 4);
            }

            if (TimestepNo > 0) {
                Assert.LessOrEqual(xL2Err, ErrorThreshold, "XDG L2 error of computing du_dx");
            }

            errorBack = xL2Err;

            // return/Ende
            //base.NoOfTimesteps = timeStep;
            //base.NoOfTimesteps = 181;
            dt = 1;
            return dt;
        }

        protected override void PlotCurrentState(double physTime, TimestepNumber timestepNo, int superSampling = 0) {
            string filename = "ZwoLsTest." + timestepNo;
            Tecplot.PlotFields(new DGField[] { u, du_dx, Phi0, Phi1, Amarker, Bmarker, Cmarker, du_dx_Exact, ERR, XERR }, filename, 0, superSampling);
        }
    }
}