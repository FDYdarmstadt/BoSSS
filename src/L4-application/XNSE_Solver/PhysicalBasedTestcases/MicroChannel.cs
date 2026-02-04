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
using System.Threading.Tasks;

using ilPSP;
using ilPSP.Utils;
using BoSSS.Solution.Control;
using BoSSS.Solution.AdvancedSolvers;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.XNSECommon;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.GridImport;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Application.XNSE_Solver.Logging;
using MathNet.Numerics.Statistics.Mcmc;
using Microsoft.CodeAnalysis.Operations;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Solution.Gnuplot;
using System.Configuration;
using ilPSP.LinSolvers.HYPRE;
using System.Drawing.Drawing2D;
using System.Runtime.Serialization;
using NUnit.Options;


namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {

    [DataContract]
    [Serializable]
    public class SharpCornerElement : IBMelement {

        [DataMember]
        MultidimensionalArray m_points;

        public SharpCornerElement(MultidimensionalArray points) {

            m_points = points;

        }

        public override double GetLevelSetFunction(double[] X) {

            double[] X2D = new double[] { X[0], X[1] };   // in case of 3D extrude in z-direction

            double dist = double.MaxValue;
            double sign = -1.0;
            for (int p = 0; p < m_points.Lengths[0] - 1; p++) {
                (double[] pA, double[] pB) = GetSegment(p);
                (double[] proj, bool isWithin) = ProjectionOnSegment(pA, pB, X2D);
                double distSegment = X2D.L2Dist(proj);
                double signSegment = (double)SignToLine(pA, pB, X2D);
                if (distSegment < dist) {
                    dist = distSegment;
                    sign = isWithin ? signSegment : -1.0;
                }
            }

            return sign * dist;
        }

        private (double[] pA, double[] pB) GetSegment(int pInd) {
            if (pInd < 0 || pInd > m_points.Lengths[0])
                throw new ArgumentOutOfRangeException();

            double[] pA = GetPoint(pInd);
            double[] pB = pInd == m_points.Lengths[0] - 1 ? GetPoint(0) : GetPoint(pInd + 1);

            return (pA, pB);
        }
        private double[] GetPoint(int pInd) {
            if (pInd < 0 || pInd > m_points.Lengths[0])
                throw new ArgumentOutOfRangeException();

            return m_points.ExtractSubArrayShallow(pInd, -1).To1DArray();
        }

        private static (double[] proj, bool isWithin) ProjectionOnSegment(double[] pA, double[] pB, double[] pP) {

            double l2 = pA.L2DistPow2(pB);

            //const float t = max(0, min(1, dot(p - v, w - v) / l2));
            //const vec2 projection = v + t * (w - v);  // Projection falls on the segment
            //return distance(p, projection);

            double[] AP = pP.CloneAs();
            AP.AccV(-1.0, pA);
            double[] AB = pB.CloneAs();
            AB.AccV(-1.0, pA);

            double dotP = AP.InnerProd(AB);
            double t = dotP / l2;
            double t01 = Math.Max(0.0, Math.Min(1.0, t));

            double[] proj = pA.CloneAs();
            proj.AccV(t01, AB);

            return (proj, t >= 0.0 && t <= 1.0);
        }

        private static int SignToLine(double[] pA, double[] pB, double[] pP) {

            //position = sign((Bx - Ax) * (Y - Ay) - (By - Ay) * (X - Ax))

            double det = (pB[0] - pA[0]) * (pP[1] - pA[1]) - (pB[1] - pA[1]) * (pP[0] - pA[0]);

            return Math.Sign(det) == 0 ? -1 : -Math.Sign(det);

        }

    }



    public static class MicroChannel {


        public static XNSE_Control SharpCornerMicroChannel3D_IBM(int k = 3) {

            XNSE_Control C = new XNSE_Control();

            C.SetDGdegree(k);

            // basic database options
            // ======================
            #region db

            string _DbPath = null; // @"D:\local\local_test_db";
            C.DbPath = _DbPath;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/MicroChannelIBM";
            C.ProjectDescription = "Micro Channel flow for Particle Separation";

            #endregion


            // Physical Parameters
            // ===================
            #region physics

            C.PhysicalParameters.rho_A = 1000;
            C.PhysicalParameters.mu_A = 1e-3;

            C.PhysicalParameters.IncludeConvection = false;

            #endregion


            // grid generation
            // ===============
            #region grid

            double SizeScale = 1.0e-6;

            double H = 2.0 * SizeScale;
            double H_center = 0.5 * SizeScale;  // defines resolution
            double H_wedge = 0.75 * SizeScale;
            double L_inflow = 2.0 * SizeScale;
            double L_wedge = 0.75 * SizeScale;
            double L_outflow = 3.0 * SizeScale;
            double W = H / 4.0;

            int resolution = 1;
            int kelemH = 8 * resolution;
            int kelemL = 23 * resolution;
            int kelemW = 2 * resolution;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(-L_inflow, L_wedge + L_outflow, kelemL + 1);
                double[] Ynodes = GenericBlas.Linspace(-H / 2.0, H / 2.0, kelemH + 1);
                //var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);
                double[] Znodes = GenericBlas.Linspace(-W / 2.0, W / 2.0, kelemW + 1);
                var grd = Grid3D.Cartesian3DGrid(Xnodes, Ynodes, Znodes);


                grd.EdgeTagNames.Add(1, "wall");
                grd.EdgeTagNames.Add(2, "velocity_inlet_left");
                grd.EdgeTagNames.Add(3, "pressure_outlet_right");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[1] + (H / 2.0)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - (H / 2.0)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] + L_inflow) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[0] - (L_wedge + L_outflow)) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[2] + (W / 2.0)) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[2] - (W / 2.0)) <= 1.0e-8)
                        et = 1;

                    return et;
                });

                return grd;
            };

            #endregion


            // Particle Properties
            // ===================
            #region particle

            C.UseImmersedBoundary = true;

            List<IBMelement> IBMelements = new List<IBMelement>();

            // lower wegde
            MultidimensionalArray pointsL = MultidimensionalArray.Create(3, 2);
            pointsL.SetSubVector(new double[] { 0.0 * SizeScale, -1.0 * SizeScale }, new int[] { 0, -1 });
            pointsL.SetSubVector(new double[] { 0.0 * SizeScale, -0.25 * SizeScale }, new int[] { 1, -1 });
            pointsL.SetSubVector(new double[] { 0.75 * SizeScale, -1.0 * SizeScale }, new int[] { 2, -1 });
            IBMelements.Add(new SharpCornerElement(pointsL));

            // upper wegde
            MultidimensionalArray pointsU = MultidimensionalArray.Create(3, 2);
            pointsU.SetSubVector(new double[] { 0.75 * SizeScale, 1.0 * SizeScale }, new int[] { 0, -1 });
            pointsU.SetSubVector(new double[] { 0.0 * SizeScale, 0.25 * SizeScale }, new int[] { 1, -1 });
            pointsU.SetSubVector(new double[] { 0.0 * SizeScale, 1.0 * SizeScale }, new int[] { 2, -1 });
            IBMelements.Add(new SharpCornerElement(pointsU));

            // particle 
            //particles.Add(new ParticleDisk(noMotion, particel_Radius, particle_InitPosition));

            //double levelSet(double[] X) {
            //    double levelSetFunction = int.MinValue;
            //    foreach (SharpCornerElement element in IBMelements) {
            //        if (levelSetFunction < element.GetLevelSetFunction(X))
            //            levelSetFunction = element.GetLevelSetFunction(X);
            //    }
            //    return levelSetFunction;
            //}
            //C.InitialValues_Evaluators.Add(VariableNames.LevelSetCGidx(1), levelSet);
            C.m_IBMelements = IBMelements;

            C.Option_LevelSetEvolution = Solution.LevelSetTools.LevelSetEvolution.None;
            //C.Option_LevelSetEvolution2 = Solution.LevelSetTools.LevelSetEvolution.None;  // via init

            //C.AdaptiveMeshRefinement = true;
            //C.activeAMRlevelIndicators.Add(new AMRaroundRigidObject(particles, H / kelemH ) { maxRefinementLevel = 1 });
            //C.AMR_startUpSweeps = 1;

            #endregion


            // Initial Values
            // ==============
            #region init


            #endregion


            // boundary conditions
            // ===================
            #region BC

            C.AddBoundaryValue("wall");

            double U = 0.6;
            //C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", (X, t) => ((-4.0 * U / H.Pow2()) * X[1].Pow2() + U)); // * Math.Sin(2.0*Math.PI*(t/T)));
            C.AddBoundaryValue("velocity_inlet_left", "VelocityX#A", (X, t) => U); // * Math.Sin(2.0*Math.PI*(t/T)));

            C.AddBoundaryValue("pressure_outlet_right");

            #endregion


            // misc. solver options
            // ====================
            #region solver


            #endregion


            // Timestepping
            // ============
            #region time

            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;

            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;

            #endregion


            return C;
        }



        public static XNSE_Control SharpCornerMicroChannel3D() {

            XNSE_Control C = new XNSE_Control();


            double SizeScale = 1.0e-6;
            double Ls = 80.0 * SizeScale;
            double Le = 200.0 * SizeScale;
            double Wi = 42.5 * SizeScale;
            double We = 200.0 * SizeScale;
            double H = 50.0 * SizeScale;

            int Resolution = 1;
            int Ls_Res = 1; // 4 * Resolution;
            int Le_Res = 8 * Resolution;
            int Wi_Res = 1; // 2 * Resolution;
            int halfWe_Res = 4 * Resolution;
            int H_Res = 1; // 2 * Resolution;

            double Linlet = 200.0 * SizeScale;
            int Linlet_Res = 8 * Resolution;

            int nSections = 1;
            double Lsection = Ls + Le;
            double Lsystem = nSections * Lsection;

            double Loutlet = 200.0 * SizeScale;
            int Loutlet_Res = 4 * Resolution;

            double Lend = Lsystem + Loutlet;

            //List<GridCommons> SystemParts = new List<GridCommons>();
            //SystemParts.Add(ExpansionGrid(-Linlet, 0.0, Wi, We, Linlet_Res, Wi_Res, halfWe_Res, false));

            //GridCommons[] Sections = new GridCommons[nSections];
            //for (int n = 0; n < nSections; n++) {
            //    Sections[n] = SectionGrid(n * Lsection, Ls, Le, Wi, We, Ls_Res, Le_Res, Wi_Res, halfWe_Res);
            //}
            //SystemParts.AddRange(Sections);

            //SystemParts.Add(ExpansionGrid(Lsystem, Lend, Wi, We, Loutlet_Res, Wi_Res, halfWe_Res));

            //var gridSystem = GridCommons.MergeLogically(SystemParts.ToArray());
            //var grd = GridCommons.Seal(gridSystem, 4);


            //grd.DefineEdgeTags(delegate (Vector X) {
            //    string ret = null;

            //    for (int n = 0; n < nSections; n++) {
            //        if (X.y <= (-Wi / 2.0) - X.x * (We - Wi) / 2.0 && (X.x - (n * Lsection)) >= -1e-8 && (X.x - (n * Lsection + Ls)) <= 1e-8)    // lower wedge
            //            ret = IncompressibleBcType.Wall.ToString();
            //        if (X.y >= (Wi / 2.0) + X.x * (We - Wi) / 2.0 && (X.x - (n * Lsection)) >= -1e-8 && (X.x - (n * Lsection + Ls)) <= 1e-8)     // upper wedge
            //            ret = IncompressibleBcType.Wall.ToString();
            //    }

            //    //if (X.y <= (-Wi / 2.0) - X.x * (We - Wi) / 2.0 && X.x >= Lsection && X.x <= Lsection + Ls)    // lower wedge
            //    //    ret = IncompressibleBcType.Wall.ToString();
            //    //if (X.y >= (Wi / 2.0) + X.x * (We - Wi) / 2.0 && X.x >= Lsection && X.x <= Lsection + Ls)     // upper wedge
            //    //    ret = IncompressibleBcType.Wall.ToString();

            //    if ((X.y + We / 2.0).Abs() <= 1e-8)
            //        ret = IncompressibleBcType.Wall.ToString();
            //    if ((X.y - We / 2.0).Abs() <= 1e-8)
            //        ret = IncompressibleBcType.Wall.ToString();

            //    if ((X.x + Linlet).Abs() <= 1e-8)
            //        ret = IncompressibleBcType.Velocity_Inlet.ToString();
            //    if ((X.x - Lend).Abs() <= 1e-8)
            //        ret = IncompressibleBcType.Pressure_Outlet.ToString();

            //    return ret;
            //});

            //GridCommons grd = LeadingUpperTriangleGrid3D(0.0, Ls, 0.0, Wi, 0.0, H, Ls_Res, Wi_Res, H_Res);
            //grd.InvalidateGridData();

            Gmsh gmshGrid = new Gmsh(@"D:\BoSSS-experimental\public\examples\MicroFluidChannel\Grids\OneSectionSystem3D.msh");
            GridCommons grd = gmshGrid.GenerateBoSSSGrid();
            grd.DefineEdgeTags(delegate (Vector X) {
                string ret = null;

                for (int n = 0; n < nSections; n++) {
                    if (X.y <= (-Wi / 2.0) - X.x * (We - Wi) / 2.0 && (X.x - (n * Lsection)) >= -1e-8 && (X.x - (n * Lsection + Ls)) <= 1e-8)    // lower wedge
                        ret = IncompressibleBcType.Wall.ToString();
                    if (X.y >= (Wi / 2.0) + X.x * (We - Wi) / 2.0 && (X.x - (n * Lsection)) >= -1e-8 && (X.x - (n * Lsection + Ls)) <= 1e-8)     // upper wedge
                        ret = IncompressibleBcType.Wall.ToString();
                }

                if ((X.y + We / 2.0).Abs() <= 1e-8)
                    ret = IncompressibleBcType.Wall.ToString();
                if ((X.y - We / 2.0).Abs() <= 1e-8)
                    ret = IncompressibleBcType.Wall.ToString();

                if ((X.x + Linlet).Abs() <= 1e-8)
                    ret = IncompressibleBcType.Velocity_Inlet.ToString();
                if ((X.x - Lend).Abs() <= 1e-8)
                    ret = IncompressibleBcType.Pressure_Outlet.ToString();

                if ((X.z + H / 2.0).Abs() <= 1e-8)
                    ret = IncompressibleBcType.Wall.ToString();
                if ((X.z - H / 2.0).Abs() <= 1e-8)
                    ret = IncompressibleBcType.Wall.ToString();

                return ret;
            });

            C.GridFunc = delegate () { return grd; };

            int k = 2;
            C.SetDGdegree(k);

            Formula InletVelocity = new Formula("X => 1.0");

            C.AddBoundaryValue(IncompressibleBcType.Velocity_Inlet.ToString(), "VelocityX#A", InletVelocity);

            C.PhysicalParameters.rho_A = 1000.0;
            C.PhysicalParameters.mu_A = 1.0e-3;
            C.PhysicalParameters.IncludeConvection = false;

            C.CutCellQuadratureType = Foundation.XDG.CutCellQuadratureMethod.OneStepGaussAndStokes;

            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;

            C.SessionName = "SCMC3D_OneSectionSystemTest";

            return C;
        }


        //public static Grid3D LeadingUpperTriangleGrid3D(double L0, double Lend, double W0, double Wend, double H0, double Hend, int L_Res, int W_Res, int H_Res) {

        //    Grid3D grid = new Grid3D(Tetra.Instance);

        //    int J = L_Res * W_Res * H_Res * 6;
        //    grid.Cells = new Cell[J];

        //    int num_xNodes = L_Res + 1;
        //    double[] xNodes = GenericBlas.Linspace(L0, Lend, num_xNodes);
        //    int num_yNodes = W_Res + 1;
        //    double[] yNodes = GenericBlas.Linspace(W0, Wend, num_yNodes);
        //    int num_zNodes = H_Res + 1;
        //    double[] zNodes = GenericBlas.Linspace(H0, Hend, num_zNodes);

        //    return Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes);

        //    int globalIndex = 0;

        //    for (int zId = 0; zId < num_zNodes - 1; zId++) {
        //        for (int yId = 0; yId < num_yNodes - 1; yId++) {
        //            for (int xId = 0; xId < num_xNodes - 1; xId++) {

        //                {
        //                    Cell Cj = new Cell();

        //                    Cj.GlobalID = globalIndex;
        //                    Cj.Type = CellType.Tetra_Linear;
        //                    Cj.TransformationParams = MultidimensionalArray.Create(4, 3);

        //                    Cj.TransformationParams[0, 0] = xNodes[xId];
        //                    Cj.TransformationParams[0, 1] = yNodes[yId + 1];
        //                    Cj.TransformationParams[0, 2] = zNodes[zId];

        //                    Cj.TransformationParams[1, 0] = xNodes[xId + 1];
        //                    Cj.TransformationParams[1, 1] = yNodes[yId + 1];
        //                    Cj.TransformationParams[1, 2] = zNodes[zId];

        //                    Cj.TransformationParams[2, 0] = xNodes[xId];
        //                    Cj.TransformationParams[2, 1] = yNodes[yId];
        //                    Cj.TransformationParams[2, 2] = zNodes[zId];

        //                    Cj.TransformationParams[3, 0] = xNodes[xId];
        //                    Cj.TransformationParams[3, 1] = yNodes[yId + 1];
        //                    Cj.TransformationParams[3, 2] = zNodes[zId + 1];

        //                    Cj.NodeIndices = new long[] {
        //                                zId * (num_xNodes + num_yNodes) + (yId + 1) * num_xNodes + xId,
        //                                zId * (num_xNodes + num_yNodes) + (yId + 1) * num_xNodes + (xId + 1),
        //                                zId * (num_xNodes + num_yNodes) + yId * num_xNodes + xId,
        //                                (zId + 1) * (num_xNodes + num_yNodes) + (yId + 1) * num_xNodes + xId
        //                    };

        //                    grid.Cells[globalIndex] = Cj;
        //                    globalIndex++;
        //                }

        //                {
        //                    Cell Cj = new Cell();

        //                    Cj.GlobalID = globalIndex;
        //                    Cj.Type = CellType.Tetra_Linear;
        //                    Cj.TransformationParams = MultidimensionalArray.Create(4, 3);

        //                    Cj.TransformationParams[0, 0] = xNodes[xId];
        //                    Cj.TransformationParams[0, 1] = yNodes[yId];
        //                    Cj.TransformationParams[0, 2] = zNodes[zId];

        //                    Cj.TransformationParams[1, 0] = xNodes[xId];
        //                    Cj.TransformationParams[1, 1] = yNodes[yId];
        //                    Cj.TransformationParams[1, 2] = zNodes[zId + 1];

        //                    Cj.TransformationParams[2, 0] = xNodes[xId];
        //                    Cj.TransformationParams[2, 1] = yNodes[yId + 1];
        //                    Cj.TransformationParams[2, 2] = zNodes[zId + 1];

        //                    Cj.TransformationParams[3, 0] = xNodes[xId + 1];
        //                    Cj.TransformationParams[3, 1] = yNodes[yId + 1];
        //                    Cj.TransformationParams[3, 2] = zNodes[zId + 1];

        //                    Cj.NodeIndices = new long[] {
        //                                (zId) * (num_xNodes + num_yNodes) + (yId) * num_xNodes + xId,
        //                                (zId + 1) * (num_xNodes + num_yNodes) + (yId) * num_xNodes + (xId),
        //                                (zId + 1) * (num_xNodes + num_yNodes) + (yId + 1) * num_xNodes + xId,
        //                                (zId + 1) * (num_xNodes + num_yNodes) + (yId + 1) * num_xNodes + (xId + 1)
        //                    };

        //                    grid.Cells[globalIndex] = Cj;
        //                    globalIndex++;
        //                }

        //                {
        //                    Cell Cj = new Cell();

        //                    Cj.GlobalID = globalIndex;
        //                    Cj.Type = CellType.Tetra_Linear;
        //                    Cj.TransformationParams = MultidimensionalArray.Create(4, 3);

        //                    Cj.TransformationParams[0, 0] = xNodes[xId];
        //                    Cj.TransformationParams[0, 1] = yNodes[yId + 1];
        //                    Cj.TransformationParams[0, 2] = zNodes[zId + 1];

        //                    Cj.TransformationParams[1, 0] = xNodes[xId + 1];
        //                    Cj.TransformationParams[1, 1] = yNodes[yId + 1];
        //                    Cj.TransformationParams[1, 2] = zNodes[zId + 1];

        //                    Cj.TransformationParams[2, 0] = xNodes[xId + 1];
        //                    Cj.TransformationParams[2, 1] = yNodes[yId + 1];
        //                    Cj.TransformationParams[2, 2] = zNodes[zId];

        //                    Cj.TransformationParams[3, 0] = xNodes[xId];
        //                    Cj.TransformationParams[3, 1] = yNodes[yId];
        //                    Cj.TransformationParams[3, 2] = zNodes[zId];

        //                    Cj.NodeIndices = new long[] {
        //                                (zId + 1) * (num_xNodes + num_yNodes) + (yId + 1) * num_xNodes + (xId),
        //                                (zId + 1) * (num_xNodes + num_yNodes) + (yId + 1) * num_xNodes + (xId + 1),
        //                                (zId) * (num_xNodes + num_yNodes) + (yId + 1) * num_xNodes + (xId + 1),
        //                                (zId) * (num_xNodes + num_yNodes) + (yId) * num_xNodes + (xId)
        //                    };

        //                    grid.Cells[globalIndex] = Cj;
        //                    globalIndex++;
        //                }

        //                {
        //                    Cell Cj = new Cell();

        //                    Cj.GlobalID = globalIndex;
        //                    Cj.Type = CellType.Tetra_Linear;
        //                    Cj.TransformationParams = MultidimensionalArray.Create(4, 3);

        //                    Cj.TransformationParams[0, 0] = xNodes[xId + 1];
        //                    Cj.TransformationParams[0, 1] = yNodes[yId + 1];
        //                    Cj.TransformationParams[0, 2] = zNodes[zId];

        //                    Cj.TransformationParams[1, 0] = xNodes[xId + 1];
        //                    Cj.TransformationParams[1, 1] = yNodes[yId];
        //                    Cj.TransformationParams[1, 2] = zNodes[zId];

        //                    Cj.TransformationParams[2, 0] = xNodes[xId];
        //                    Cj.TransformationParams[2, 1] = yNodes[yId];
        //                    Cj.TransformationParams[2, 2] = zNodes[zId];

        //                    Cj.TransformationParams[3, 0] = xNodes[xId + 1];
        //                    Cj.TransformationParams[3, 1] = yNodes[yId];
        //                    Cj.TransformationParams[3, 2] = zNodes[zId + 1];

        //                    Cj.NodeIndices = new long[] {
        //                                                        (zId) * (num_xNodes + num_yNodes) + (yId + 1) * num_xNodes + (xId + 1),
        //                                                        (zId) * (num_xNodes + num_yNodes) + (yId) * num_xNodes + (xId + 1),
        //                                                        (zId) * (num_xNodes + num_yNodes) + (yId) * num_xNodes + (xId),
        //                                                        (zId + 1) * (num_xNodes + num_yNodes) + (yId) * num_xNodes + (xId + 1)
        //                                            };

        //                    grid.Cells[globalIndex] = Cj;
        //                    globalIndex++;
        //                }

        //                {
        //                    Cell Cj = new Cell();

        //                    Cj.GlobalID = globalIndex;
        //                    Cj.Type = CellType.Tetra_Linear;
        //                    Cj.TransformationParams = MultidimensionalArray.Create(4, 3);

        //                    Cj.TransformationParams[0, 0] = xNodes[xId];
        //                    Cj.TransformationParams[0, 1] = yNodes[yId];
        //                    Cj.TransformationParams[0, 2] = zNodes[zId];

        //                    Cj.TransformationParams[1, 0] = xNodes[xId + 1];
        //                    Cj.TransformationParams[1, 1] = yNodes[yId];
        //                    Cj.TransformationParams[1, 2] = zNodes[zId + 1];

        //                    Cj.TransformationParams[2, 0] = xNodes[xId];
        //                    Cj.TransformationParams[2, 1] = yNodes[yId];
        //                    Cj.TransformationParams[2, 2] = zNodes[zId + 1];

        //                    Cj.TransformationParams[3, 0] = xNodes[xId + 1];
        //                    Cj.TransformationParams[3, 1] = yNodes[yId + 1];
        //                    Cj.TransformationParams[3, 2] = zNodes[zId + 1];

        //                    Cj.NodeIndices = new long[] {
        //                                                        (zId) * (num_xNodes + num_yNodes) + (yId) * num_xNodes + (xId),
        //                                                        (zId + 1) * (num_xNodes + num_yNodes) + (yId) * num_xNodes + (xId + 1),
        //                                                        (zId + 1) * (num_xNodes + num_yNodes) + (yId) * num_xNodes + (xId),
        //                                                        (zId + 1) * (num_xNodes + num_yNodes) + (yId + 1) * num_xNodes + (xId + 1)
        //                                            };

        //                    grid.Cells[globalIndex] = Cj;
        //                    globalIndex++;
        //                }

        //                {
        //                    Cell Cj = new Cell();

        //                    Cj.GlobalID = globalIndex;
        //                    Cj.Type = CellType.Tetra_Linear;
        //                    Cj.TransformationParams = MultidimensionalArray.Create(4, 3);


        //                    Cj.TransformationParams[0, 0] = xNodes[xId + 1];
        //                    Cj.TransformationParams[0, 1] = yNodes[yId + 1];
        //                    Cj.TransformationParams[0, 2] = zNodes[zId + 1];

        //                    Cj.TransformationParams[1, 0] = xNodes[xId + 1];
        //                    Cj.TransformationParams[1, 1] = yNodes[yId];
        //                    Cj.TransformationParams[1, 2] = zNodes[zId + 1];

        //                    Cj.TransformationParams[2, 0] = xNodes[xId + 1];
        //                    Cj.TransformationParams[2, 1] = yNodes[yId + 1];
        //                    Cj.TransformationParams[2, 2] = zNodes[zId];

        //                    Cj.TransformationParams[3, 0] = xNodes[xId];
        //                    Cj.TransformationParams[3, 1] = yNodes[yId];
        //                    Cj.TransformationParams[3, 2] = zNodes[zId];

        //                    Cj.NodeIndices = new long[] {
        //                                                        (zId + 1) * (num_xNodes + num_yNodes) + (yId + 1) * num_xNodes + (xId + 1),
        //                                                        (zId + 1) * (num_xNodes + num_yNodes) + (yId) * num_xNodes + (xId + 1),
        //                                                        (zId) * (num_xNodes + num_yNodes) + (yId + 1) * num_xNodes + (xId + 1),
        //                                                        (zId) * (num_xNodes + num_yNodes) + (yId) * num_xNodes + (xId)
        //                                            };

        //                    grid.Cells[globalIndex] = Cj;
        //                    globalIndex++;
        //                }

        //            }
        //        }
        //    }

        //    return grid;
        //}




        public static XNSE_Control SharpCornerMicroChannel() {

            XNSE_Control C = new XNSE_Control();


            double SizeScale = 1.0e-6;
            double Ls = 80.0 * SizeScale;
            double Le = 200.0 * SizeScale;
            double Wi = 42.5 * SizeScale;
            double We = 200.0 * SizeScale;

            int Resolution = 1;
            int Ls_Res = 4 * Resolution;
            int Le_Res = 8 * Resolution;
            int Wi_Res = 2 * Resolution;
            int halfWe_Res = 4 * Resolution;

            double Linlet = 200.0 * SizeScale;
            int Linlet_Res = 8 * Resolution;

            int nSections = 1;
            double Lsection = Ls + Le;
            double Lsystem = nSections * Lsection;

            double Loutlet = 200.0 * SizeScale;
            int Loutlet_Res = 4 * Resolution;

            double Lend = Lsystem + Loutlet;

            //List<GridCommons> SystemParts = new List<GridCommons>();
            //SystemParts.Add(ExpansionGrid(-Linlet, 0.0, Wi, We, Linlet_Res, Wi_Res, halfWe_Res, false));

            //GridCommons[] Sections = new GridCommons[nSections];
            //for (int n = 0; n < nSections; n++) {
            //    Sections[n] = SectionGrid(n * Lsection, Ls, Le, Wi, We, Ls_Res, Le_Res, Wi_Res, halfWe_Res);
            //}
            //SystemParts.AddRange(Sections);

            //SystemParts.Add(ExpansionGrid(Lsystem, Lend, Wi, We, Loutlet_Res, Wi_Res, halfWe_Res));

            //var gridSystem = GridCommons.MergeLogically(SystemParts.ToArray());
            //var grd = GridCommons.Seal(gridSystem, 4);

            Gmsh gmshGrid = new Gmsh(@"D:\BoSSS-experimental\public\examples\MicroFluidChannel\Grids\OneSectionSystem2D.msh");
            GridCommons grd = gmshGrid.GenerateBoSSSGrid();


            grd.DefineEdgeTags(delegate (Vector X) {
                string ret = null;

                for (int n = 0; n < nSections; n++) {
                    if (X.y <= (-Wi / 2.0) - X.x * (We - Wi) / 2.0 && (X.x - (n * Lsection)) >= -1e-8 && (X.x - (n * Lsection + Ls)) <= 1e-8)    // lower wedge
                        ret = IncompressibleBcType.Wall.ToString();
                    if (X.y >= (Wi / 2.0) + X.x * (We - Wi) / 2.0 && (X.x - (n * Lsection)) >= -1e-8 && (X.x - (n * Lsection + Ls)) <= 1e-8)     // upper wedge
                        ret = IncompressibleBcType.Wall.ToString();
                }

                //if (X.y <= (-Wi / 2.0) - X.x * (We - Wi) / 2.0 && X.x >= Lsection && X.x <= Lsection + Ls)    // lower wedge
                //    ret = IncompressibleBcType.Wall.ToString();
                //if (X.y >= (Wi / 2.0) + X.x * (We - Wi) / 2.0 && X.x >= Lsection && X.x <= Lsection + Ls)     // upper wedge
                //    ret = IncompressibleBcType.Wall.ToString();

                if ((X.y + We / 2.0).Abs() <= 1e-8)
                    ret = IncompressibleBcType.Wall.ToString();
                if ((X.y - We / 2.0).Abs() <= 1e-8)
                    ret = IncompressibleBcType.Wall.ToString();

                if ((X.x + Linlet).Abs() <= 1e-8)
                    ret = IncompressibleBcType.Velocity_Inlet.ToString();
                if ((X.x - Lend).Abs() <= 1e-8)
                    ret = IncompressibleBcType.Pressure_Outlet.ToString();

                return ret;
            });

            C.GridFunc = delegate () { return grd; };

            int k = 3;
            C.SetDGdegree(k);

            Formula InletVelocity_AtExpansion = new Formula(
                "VelX",
                false,
                "double VelX(double[] X) { " +
                "double WeHalf = 100.0e-6;" +
                "double U_max = (3.0/2.0)*0.25;" +
                "return (-U_max / WeHalf.Pow2()) * X[1].Pow2() + U_max; } "
            );

            Formula InletVelocity_AtContraction = new Formula(
                "VelX",
                false,
                "double VelX(double[] X) { " +
                "double WiHalf = 21.25e-6;" +
                "double U_max = (3.0/2.0)*1.0;" +
                "return (-U_max / WiHalf.Pow2()) * X[1].Pow2() + U_max; } "
            );

            C.AddBoundaryValue(IncompressibleBcType.Velocity_Inlet.ToString(), "VelocityX#A", InletVelocity_AtExpansion);

            C.PhysicalParameters.rho_A = 1000.0;
            C.PhysicalParameters.mu_A = 1.0e-3;
            C.PhysicalParameters.IncludeConvection = true;

            C.CutCellQuadratureType = Foundation.XDG.CutCellQuadratureMethod.OneStepGaussAndStokes;

            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;

            C.SessionName = "SCMC_SetupTest_TriangleGridGmsh";

            return C;
        }


        public static GridCommons SectionGrid(double L0, double Ls, double Le, double Wi, double We, int Ls_Res, int Le_Res, int Wi_Res, int halfWe_Res) {

            var gridCenter = CenterlineGrid(L0, L0 + Ls, Wi, Ls_Res, Wi_Res);

            var gridUpperExp = UpperExpansionGrid(L0, L0 + Ls, Wi / 2.0, We / 2.0, Ls_Res, halfWe_Res);

            var gridLowerExp = LowerExpansionGrid(L0, L0 + Ls, Wi / 2.0, We / 2.0, Ls_Res, halfWe_Res);

            var gridExpansion = ExpansionGridSymmetric(L0 + Ls, L0 + (Ls + Le), Wi, We, Le_Res, Wi_Res, halfWe_Res);


            var gridSection = GridCommons.MergeLogically(new GridCommons[] { gridCenter, gridUpperExp, gridLowerExp, gridExpansion });
            var grd = GridCommons.Seal(gridSection, 4);

            return grd;
        }


        public static GridCommons UpperExpansionGrid(double L0, double Lend, double halfWi, double halfWe, int L_Res, int halfW_Res) {

            Grid2D grid = new Grid2D(Triangle.Instance);

            int J = L_Res * halfW_Res;
            grid.Cells = new Cell[J];

            int num_xNodes = L_Res + 1;
            double[] xNodes = GenericBlas.Linspace(L0, Lend, num_xNodes);
            int num_yNodes = halfW_Res + 1;
            double[] yNodes = GenericBlas.Linspace(halfWi, halfWe, num_yNodes);  

            int globalIndex = 0;

            int xId_start = 0;
            for (int yId = 0; yId < num_yNodes - 1; yId++) {
                bool skipFirstFlag = true;
                for (int xId = xId_start; xId < num_xNodes - 1; xId++) {

                    if (!skipFirstFlag) {
                        Cell Cj = new Cell();

                        Cj.GlobalID = globalIndex;
                        Cj.Type = CellType.Triangle_3;
                        Cj.TransformationParams = MultidimensionalArray.Create(3, 2);

                        Cj.TransformationParams[0, 0] = xNodes[xId];
                        Cj.TransformationParams[0, 1] = yNodes[yId];
                        Cj.TransformationParams[1, 0] = xNodes[xId + 1];
                        Cj.TransformationParams[1, 1] = yNodes[yId + 1];
                        Cj.TransformationParams[2, 0] = xNodes[xId];
                        Cj.TransformationParams[2, 1] = yNodes[yId + 1];

                        Cj.NodeIndices = new long[] {
                                        yId * num_xNodes + xId,
                                        (yId + 1) * num_xNodes + (xId + 1),
                                        (yId + 1) * num_xNodes + xId };

                        grid.Cells[globalIndex] = Cj;
                        globalIndex++;
                    }
                    skipFirstFlag = false;

                    {
                        Cell Cj = new Cell();

                        Cj.GlobalID = globalIndex;
                        Cj.Type = CellType.Triangle_3;
                        Cj.TransformationParams = MultidimensionalArray.Create(3, 2);

                        Cj.TransformationParams[0, 0] = xNodes[xId];
                        Cj.TransformationParams[0, 1] = yNodes[yId];
                        Cj.TransformationParams[1, 0] = xNodes[xId + 1];
                        Cj.TransformationParams[1, 1] = yNodes[yId];
                        Cj.TransformationParams[2, 0] = xNodes[xId + 1];
                        Cj.TransformationParams[2, 1] = yNodes[yId + 1];

                        Cj.NodeIndices = new long[] {
                                        yId * num_xNodes + xId,
                                        yId * num_xNodes + (xId + 1),
                                        (yId + 1) * num_xNodes + (xId + 1) };

                        grid.Cells[globalIndex] = Cj;
                        globalIndex++;
                    }
                }
                xId_start++;
            }

            return grid;
        }

        public static GridCommons LowerExpansionGrid(double L0, double Lend, double halfWi, double halfWe, int L_Res, int halfW_Res) {

            Grid2D grid = new Grid2D(Triangle.Instance);

            int J = L_Res * halfW_Res;
            grid.Cells = new Cell[J];

            int num_xNodes = L_Res + 1;
            double[] xNodes = GenericBlas.Linspace(L0, Lend, num_xNodes);
            int num_yNodes = halfW_Res + 1;
            double[] yNodes = GenericBlas.Linspace(-halfWe, -halfWi, num_yNodes);

            int globalIndex = 0;

            int xId_start = num_xNodes - 2;
            for (int yId = 0; yId < num_yNodes - 1; yId++) {
                bool skipFirstFlag = true;
                for (int xId = xId_start; xId < num_xNodes - 1; xId++) {

                    if (!skipFirstFlag) {
                        Cell Cj = new Cell();

                        Cj.GlobalID = globalIndex;
                        Cj.Type = CellType.Triangle_3;
                        Cj.TransformationParams = MultidimensionalArray.Create(3, 2);

                        Cj.TransformationParams[0, 0] = xNodes[xId];
                        Cj.TransformationParams[0, 1] = yNodes[yId];
                        Cj.TransformationParams[1, 0] = xNodes[xId + 1];
                        Cj.TransformationParams[1, 1] = yNodes[yId];
                        Cj.TransformationParams[2, 0] = xNodes[xId];
                        Cj.TransformationParams[2, 1] = yNodes[yId + 1];

                        Cj.NodeIndices = new long[] {
                                        ((num_yNodes - 1) + yId) * num_xNodes + xId,
                                        ((num_yNodes - 1) + yId) * num_xNodes + (xId + 1),
                                        ((num_yNodes - 1) + (yId + 1)) * num_xNodes + xId };

                        grid.Cells[globalIndex] = Cj;
                        globalIndex++;
                    }
                    skipFirstFlag = false;

                    {
                        Cell Cj = new Cell();

                        Cj.GlobalID = globalIndex;
                        Cj.Type = CellType.Triangle_3;
                        Cj.TransformationParams = MultidimensionalArray.Create(3, 2);

                        Cj.TransformationParams[0, 0] = xNodes[xId + 1];
                        Cj.TransformationParams[0, 1] = yNodes[yId];
                        Cj.TransformationParams[1, 0] = xNodes[xId + 1];
                        Cj.TransformationParams[1, 1] = yNodes[yId + 1];
                        Cj.TransformationParams[2, 0] = xNodes[xId];
                        Cj.TransformationParams[2, 1] = yNodes[yId + 1];

                        Cj.NodeIndices = new long[] {
                                        ((num_yNodes - 1) + yId) * num_xNodes + (xId + 1),
                                        ((num_yNodes - 1) + (yId + 1)) * num_xNodes + (xId + 1),
                                        ((num_yNodes - 1) + (yId + 1)) * num_xNodes + xId };

                        grid.Cells[globalIndex] = Cj;
                        globalIndex++;
                    }
                }
                xId_start--;
            }

            return grid;
        }


        public static GridCommons ExpansionGrid(double L0, double Lend, double Wi, double We, int L_Res, int Wi_Res, int halfW_Res, bool forwardCenterlineOrientation = true) {

            var CenterGrid = CenterlineGrid(L0, Lend, Wi, L_Res, Wi_Res, forwardCenterlineOrientation);

            var LowerSlabGrid = LeadingUpperTriangleGrid(L0, Lend, -We / 2.0, -Wi / 2.0, L_Res, halfW_Res);

            var UpperSlabGrid = LeadingLowerTriangleGrid(L0, Lend, Wi / 2.0, We / 2.0, L_Res, halfW_Res);

            var gridExpansion = GridCommons.MergeLogically(new GridCommons[] { CenterGrid, LowerSlabGrid, UpperSlabGrid });
            var grd = GridCommons.Seal(gridExpansion, 4);

            return grd;
        }

        public static GridCommons ExpansionGridSymmetric(double L0, double Lend, double Wi, double We, int L_Res, int Wi_Res, int halfW_Res, bool forwardCenterlineOrientation = true) {

            double Lend1 = L0 + (Lend - L0) / 2.0;

            var CenterGrid1 = CenterlineGrid(L0, Lend1, Wi, L_Res / 2, Wi_Res, forwardCenterlineOrientation);

            var LowerSlabGrid1 = LeadingUpperTriangleGrid(L0, Lend1, -We / 2.0, -Wi / 2.0, L_Res / 2, halfW_Res);

            var UpperSlabGrid1 = LeadingLowerTriangleGrid(L0, Lend1, Wi / 2.0, We / 2.0, L_Res / 2, halfW_Res);

            var CenterGrid2 = CenterlineGrid(Lend1, Lend, Wi, L_Res / 2, Wi_Res, !forwardCenterlineOrientation);

            var LowerSlabGrid2 = LeadingLowerTriangleGrid(Lend1, Lend, -We / 2.0, -Wi / 2.0, L_Res / 2, halfW_Res);

            var UpperSlabGrid2 = LeadingUpperTriangleGrid(Lend1, Lend, Wi / 2.0, We / 2.0, L_Res / 2, halfW_Res); 

            var gridExpansion = GridCommons.MergeLogically(new GridCommons[] { CenterGrid1, LowerSlabGrid1, UpperSlabGrid1, CenterGrid2, LowerSlabGrid2, UpperSlabGrid2 });
            var grd = GridCommons.Seal(gridExpansion, 4);

            return grd;
        }

        public static GridCommons CenterlineGrid(double L0, double Lend, double Wi, int L_Res, int Wi_Res, bool forwardOrientation = true) {

            GridCommons LowerSlabGrid;
            GridCommons UpperSlabGrid;
            if (forwardOrientation) {
                LowerSlabGrid = LeadingUpperTriangleGrid(L0, Lend, -Wi / 2.0, 0.0, L_Res, Wi_Res / 2);
                UpperSlabGrid = LeadingLowerTriangleGrid(L0, Lend, 0.0, Wi / 2.0, L_Res, Wi_Res / 2);
            } else {
                LowerSlabGrid = LeadingLowerTriangleGrid(L0, Lend, -Wi / 2.0, 0.0, L_Res, Wi_Res / 2); 
                UpperSlabGrid = LeadingUpperTriangleGrid(L0, Lend, 0.0, Wi / 2.0, L_Res, Wi_Res / 2);
            }

            var gridCenter = GridCommons.MergeLogically(new GridCommons[] { LowerSlabGrid, UpperSlabGrid });
            var grd = GridCommons.Seal(gridCenter, 4);

            return grd;
        }


        public static Grid2D LeadingUpperTriangleGrid(double L0, double Lend, double W0, double Wend, int L_Res, int W_Res) {

            Grid2D grid = new Grid2D(Triangle.Instance);

            int J = L_Res * W_Res * 2;
            grid.Cells = new Cell[J];

            int num_xNodes = L_Res + 1;
            double[] xNodes = GenericBlas.Linspace(L0, Lend, num_xNodes);
            int num_yNodes = W_Res + 1;
            double[] yNodes = GenericBlas.Linspace(W0, Wend, num_yNodes);   

            int globalIndex = 0;

            for (int yId = 0; yId < num_yNodes - 1; yId++) {
                for (int xId = 0; xId < num_xNodes - 1; xId++) {

                    {
                        Cell Cj = new Cell();

                        Cj.GlobalID = globalIndex;
                        Cj.Type = CellType.Triangle_3;
                        Cj.TransformationParams = MultidimensionalArray.Create(3, 2);

                        Cj.TransformationParams[0, 0] = xNodes[xId];
                        Cj.TransformationParams[0, 1] = yNodes[yId];
                        Cj.TransformationParams[1, 0] = xNodes[xId + 1];
                        Cj.TransformationParams[1, 1] = yNodes[yId + 1];
                        Cj.TransformationParams[2, 0] = xNodes[xId];
                        Cj.TransformationParams[2, 1] = yNodes[yId + 1];

                        Cj.NodeIndices = new long[] {
                                        yId * num_xNodes + xId,
                                        (yId + 1) * num_xNodes + (xId + 1),
                                        (yId + 1) * num_xNodes + xId };

                        grid.Cells[globalIndex] = Cj;
                        globalIndex++;
                    }

                    {
                        Cell Cj = new Cell();

                        Cj.GlobalID = globalIndex;
                        Cj.Type = CellType.Triangle_3;
                        Cj.TransformationParams = MultidimensionalArray.Create(3, 2);

                        Cj.TransformationParams[0, 0] = xNodes[xId];
                        Cj.TransformationParams[0, 1] = yNodes[yId];
                        Cj.TransformationParams[1, 0] = xNodes[xId + 1];
                        Cj.TransformationParams[1, 1] = yNodes[yId];
                        Cj.TransformationParams[2, 0] = xNodes[xId + 1];
                        Cj.TransformationParams[2, 1] = yNodes[yId + 1];

                        Cj.NodeIndices = new long[] {
                                        yId * num_xNodes + xId,
                                        yId * num_xNodes + (xId + 1),
                                        (yId + 1) * num_xNodes + (xId + 1) };

                        grid.Cells[globalIndex] = Cj;
                        globalIndex++;
                    }
                }
            }

            return grid;
        }

        public static Grid2D LeadingLowerTriangleGrid(double L0, double Lend, double W0, double Wend, int L_Res, int W_Res) {

            Grid2D grid = new Grid2D(Triangle.Instance);

            int J = L_Res * W_Res * 2;
            grid.Cells = new Cell[J];

            int num_xNodes = L_Res + 1;
            double[] xNodes = GenericBlas.Linspace(L0, Lend, num_xNodes);
            int num_yNodes = W_Res + 1;
            double[] yNodes = GenericBlas.Linspace(W0, Wend, num_yNodes);

            int globalIndex = 0;

            for (int yId = 0; yId < num_yNodes - 1; yId++) {
                for (int xId = 0; xId < num_xNodes - 1; xId++) {

                    {
                        Cell Cj = new Cell();

                        Cj.GlobalID = globalIndex;
                        Cj.Type = CellType.Triangle_3;
                        Cj.TransformationParams = MultidimensionalArray.Create(3, 2);

                        Cj.TransformationParams[0, 0] = xNodes[xId];
                        Cj.TransformationParams[0, 1] = yNodes[yId];
                        Cj.TransformationParams[1, 0] = xNodes[xId + 1];
                        Cj.TransformationParams[1, 1] = yNodes[yId];
                        Cj.TransformationParams[2, 0] = xNodes[xId];
                        Cj.TransformationParams[2, 1] = yNodes[yId + 1];

                        Cj.NodeIndices = new long[] {
                                        ((num_yNodes - 1) + yId) * num_xNodes + xId,
                                        ((num_yNodes - 1) + yId) * num_xNodes + (xId + 1),
                                        ((num_yNodes - 1) + (yId + 1)) * num_xNodes + xId };

                        grid.Cells[globalIndex] = Cj;
                        globalIndex++;
                    }

                    {
                        Cell Cj = new Cell();

                        Cj.GlobalID = globalIndex;
                        Cj.Type = CellType.Triangle_3;
                        Cj.TransformationParams = MultidimensionalArray.Create(3, 2);

                        Cj.TransformationParams[0, 0] = xNodes[xId + 1];
                        Cj.TransformationParams[0, 1] = yNodes[yId];
                        Cj.TransformationParams[1, 0] = xNodes[xId + 1];
                        Cj.TransformationParams[1, 1] = yNodes[yId + 1];
                        Cj.TransformationParams[2, 0] = xNodes[xId];
                        Cj.TransformationParams[2, 1] = yNodes[yId + 1];

                        Cj.NodeIndices = new long[] {
                                        ((num_yNodes - 1) + yId) * num_xNodes + (xId + 1),
                                        ((num_yNodes - 1) + (yId + 1)) * num_xNodes + (xId + 1),
                                        ((num_yNodes - 1) + (yId + 1)) * num_xNodes + xId };

                        grid.Cells[globalIndex] = Cj;
                        globalIndex++;
                    }
                }
            }

            return grid;
        }



        public static XNSE_Control MOFFMicroChannel() {

            XNSE_Control C = new XNSE_Control();

            double density = 1000.0;
            double viscosity = 1.0e-3;
            double averBulkVelocity = 1.0;
            double SizeScale = 1.0e-6;
            double H = 50.0 * SizeScale;   // channel height
            double We = 200.0 * SizeScale;
            double Wi = 42.5 * SizeScale;


            double Le = 200.0 * SizeScale;
            double Ls = 100.0 * SizeScale;
            double L = Le + Ls;

            int FlowResolution = 1;

            double Linlet = 200.0 * SizeScale;

            int InletStreamWiseResolution = 2 * FlowResolution;
            int ContractionCrossWiseResolution = 1 * FlowResolution;
            int ContractionStreamWiseResolution = 2 * FlowResolution;
            int OuterExpansionCrossWiseResolution = 1 * FlowResolution;

            GridCommons gridInlet;
            {
                double[] xNodes1 = GenericBlas.Linspace(-Linlet, 0.0, InletStreamWiseResolution + 1);
                double[] yNodes1 = GenericBlas.Linspace(-Wi / 2.0, Wi / 2.0, ContractionCrossWiseResolution + 1);
                var grd1 = Grid2D.Cartesian2DGrid(xNodes1, yNodes1);

                double[] xNodes2 = xNodes1;
                double[] yNodes2 = GenericBlas.Linspace(-We / 2.0, -Wi / 2.0, OuterExpansionCrossWiseResolution + 1);
                var grd2 = Grid2D.Cartesian2DGrid(xNodes2, yNodes2);

                double[] xNodes3 = xNodes1;
                double[] yNodes3 = GenericBlas.Linspace(Wi / 2.0, We / 2.0, OuterExpansionCrossWiseResolution + 1);
                var grd3 = Grid2D.Cartesian2DGrid(xNodes3, yNodes3);

                gridInlet = GridCommons.MergeLogically(new GridCommons[] { grd1, grd2, grd3 });
            }
            var gridInletSealed = GridCommons.Seal(gridInlet, 4);

            gridInletSealed.DefineEdgeTags(delegate (Vector X) {
                string ret = null;

                if ((X.y + We / 2.0).Abs() <= 1e-8)
                    ret = IncompressibleBcType.Wall.ToString();
                if ((X.y - We / 2.0).Abs() <= 1e-8)
                    ret = IncompressibleBcType.Wall.ToString();

                // inlet
                if ((X.x + Linlet).Abs() <= 1e-8)
                    ret = IncompressibleBcType.Velocity_Inlet.ToString();
                if (X.x.Abs() <= 1e-8 && (X.y.Abs() - Wi / 2.0) < 1e-8)
                    ret = IncompressibleBcType.Pressure_Outlet.ToString(); //"inner edge";
                if (X.x.Abs() <= 1e-8 && (X.y.Abs() - Wi / 2.0) >= 1e-8)
                    ret = IncompressibleBcType.Wall.ToString();

                return ret;
            });


            double Loutlet = 200.0 * SizeScale;

            int OutletStreamWiseResolution = 2 * FlowResolution;

            GridCommons gridOutlet;
            {
                double[] xNodes1 = GenericBlas.Linspace(0.0, Loutlet, OutletStreamWiseResolution + 1);
                double[] yNodes1 = GenericBlas.Linspace(-Wi / 2.0, Wi / 2.0, ContractionCrossWiseResolution + 1);
                var grd1 = Grid2D.Cartesian2DGrid(xNodes1, yNodes1);

                double[] xNodes2 = xNodes1;
                double[] yNodes2 = GenericBlas.Linspace(-We / 2.0, -Wi / 2.0, OuterExpansionCrossWiseResolution + 1);
                var grd2 = Grid2D.Cartesian2DGrid(xNodes2, yNodes2);

                double[] xNodes3 = xNodes1;
                double[] yNodes3 = GenericBlas.Linspace(Wi / 2.0, We / 2.0, OuterExpansionCrossWiseResolution + 1);
                var grd3 = Grid2D.Cartesian2DGrid(xNodes3, yNodes3);

                gridOutlet = GridCommons.MergeLogically(new GridCommons[] { grd1, grd2, grd3 });
            }
            var gridOutletSealed = GridCommons.Seal(gridOutlet, 4);

            gridOutletSealed.DefineEdgeTags(delegate (Vector X) {
                string ret = null;

                if ((X.y + We / 2.0).Abs() <= 1e-8)
                    ret = IncompressibleBcType.Wall.ToString();
                if ((X.y - We / 2.0).Abs() <= 1e-8)
                    ret = IncompressibleBcType.Wall.ToString();

                // outlet
                if ((X.x - Loutlet).Abs() <= 1e-8)
                    ret = IncompressibleBcType.Pressure_Outlet.ToString();
                if (X.x.Abs() <= 1e-8 && (X.y.Abs() - Wi / 2.0) < 1e-8)
                    ret = "inner edge";
                if (X.x.Abs() <= 1e-8 && (X.y.Abs() - Wi / 2.0) >= 1e-8)
                    ret = IncompressibleBcType.Wall.ToString();

                return ret;
            });

            var gridSystem = GridCommons.MergeLogically(new GridCommons[] { gridInletSealed, gridOutletSealed });
            var grd = GridCommons.Seal(gridSystem, 4);

            C.GridFunc = delegate () { return grd; };
            int k = 2;
            C.SetDGdegree(k);

            C.PhysicalParameters.rho_A = density;
            C.PhysicalParameters.mu_A = viscosity;
            C.PhysicalParameters.IncludeConvection = true;

            Formula ExpansionInletVelocity = new Formula(
                "VelX",
                false,
                "double VelX(double[] X) { " +
                "double WeHalf = 100.0e-6;" +
                "double U_max = (3.0/2.0)*0.25;" +
                "return (-U_max / WeHalf.Pow2()) * X[1].Pow2() + U_max; } "
            );

            C.AddBoundaryValue(IncompressibleBcType.Velocity_Inlet.ToString(), "VelocityX#A", ExpansionInletVelocity);

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;

            C.SessionName = "SCMC_SetupTestInlet";

            return C;
        }

    }
}
