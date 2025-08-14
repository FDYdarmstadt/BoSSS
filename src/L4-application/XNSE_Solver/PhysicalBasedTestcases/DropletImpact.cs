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
using BoSSS.Solution.XNSECommon;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Solution.LevelSetTools.FourierLevelSet;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Application.XNSE_Solver.Logging;
using System.Configuration;
using static BoSSS.Solution.AMRLevelIndicatorLibrary;
using BoSSS.Solution;
using BoSSS.Application.XNSE_Solver.SpecificSolutions;
using System.IO;
using BoSSS.Solution.NSECommon;
using MathNet.Numerics.Statistics.Mcmc;
using BoSSS.Foundation;
using static BoSSS.Foundation.Grid.Classic.GridCommons;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {

    /// <summary>
    /// class providing Controls for droplet impact testcases
    /// </summary>
    public static class DropletImpact {


        public static XNSE_Control MovingWallBoundaryLayer(int k = 3, int numCells = 4) {

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            C.DbPath = null;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/MovingWallBundaryLayer";
            C.ProjectDescription = "test setup for resolving the boundary layer at a moving wall";

            #endregion


            // DG degrees
            // ==========
            C.SetDGdegree(k);


            // physical parameters
            // ===================
            #region physics

            //C.Tags.Add("Water at 1atm");
            //C.PhysicalParameters.rho_A = 1e3;
            //C.PhysicalParameters.mu_A = 1e-3;

            C.Tags.Add("Air at 1atm");
            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.mu_A = 1.48e-5;


            C.PhysicalParameters.IncludeConvection = true;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double domainSize = 3.0e-3;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, domainSize, numCells + 1);
                double[] Ynodes = GenericBlas.Linspace(0, 10 * domainSize, 10 * numCells + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: true);

                grd.EdgeTagNames.Add(1, "wall_lower");
                grd.EdgeTagNames.Add(2, "velocity_inlet_upper");
                //grd.EdgeTagNames.Add(2, "freestream_upper");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    //if (Math.Abs(X[0]) <= 1.0e-8)
                    //    et = 1;
                    //if (Math.Abs(X[0] - domainSize) <= 1.0e-8)
                    //    et = 2;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[1] - 10 * domainSize) <= 1.0e-8)
                        et = 2;

                    return et;
                });

                return grd;
            };

            #endregion


            // initial values
            // ==============
            #region init

            Func<double[], double> PhiFunc = X => -1.0;

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            //double g = 9.81;
            //C.InitialValues_Evaluators.Add("GravityY#A", X => -g);
            //C.InitialValues_Evaluators.Add("GravityY#B", X => -g);

            #endregion


            // boundary conditions
            // ===================
            #region BC

            double Uwall = -0.1;
            C.AddBoundaryValue("wall_lower", "VelocityX#A", X => Uwall);
            C.AddBoundaryValue("velocity_inlet_upper");
            //C.AddBoundaryValue("freestream_upper");

            #endregion


            // misc. solver options
            // ====================
            #region solver


            C.Option_LevelSetEvolution = LevelSetEvolution.None;

            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
            C.NonLinearSolver.ConvergenceCriterion = 1e-9;

            C.AdaptiveMeshRefinement = false;
            C.activeAMRlevelIndicators.Add(new AMRonBoundary(new byte[] { 1 }) { maxRefinementLevel = 1 });
            C.AMR_startUpSweeps = 1;

            #endregion


            // Timestepping
            // ============
            #region time


            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;

            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;

            #endregion


            return C;

        }


        static public Grid3D RotatingDiskSector_Linearized(double radiusOP, double l_radial, double l_azimuthal, double h_axial, int res_azimuthal, int res_radial, int res_axial) {

            double[] xNodes = GenericBlas.Linspace(-l_radial / 2.0, l_radial / 2.0, res_radial + 1);    // radial direction
            double[] yNodes = GenericBlas.Linspace(-l_azimuthal / 2.0, l_azimuthal / 2.0, res_azimuthal + 1);    // azimuthal direction
            double[] zNodes = GenericBlas.Linspace(0.0, h_axial, res_axial + 1);    // axial direction

            var grd = Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes, periodicY: false);
            grd.Name = $"RotatingDiskSector3D_Linearized_{res_radial}x{res_azimuthal}x{res_axial}";

            grd.EdgeTagNames.Add(1, "wall_rotatingDisk");
            grd.EdgeTagNames.Add(2, "velocity_inlet_top");
            grd.EdgeTagNames.Add(3, "velocity_inlet_front");
            grd.EdgeTagNames.Add(4, "velocity_inlet_back");
            grd.EdgeTagNames.Add(5, "velocity_inlet_upstream");
            grd.EdgeTagNames.Add(6, "velocity_inlet_downstream");

            grd.DefineEdgeTags(delegate (Vector X) {
                byte et = 0;
                if (X.z.Abs() <= 1e-8)
                    et = 1;
                if ((X.z - h_axial).Abs() <= 1e-8)
                    et = 2;
                if ((X.x + (l_radial / 2.0)).Abs() <= 1e-8)
                    et = 3;
                if ((X.x - (l_radial / 2.0)).Abs() <= 1e-8)
                    et = 4;
                if ((X.y + (l_azimuthal / 2.0)).Abs() <= 1e-8)
                    et = 5;
                if ((X.y - (l_azimuthal / 2.0)).Abs() <= 1e-8)
                    et = 6;

                return et;
            });

            return grd;
        }


        public static XNSE_Control RotatingDiskBoundaryLayer(int k = 3, int gridRes = 2) {

            var C = new XNSE_Control();

            C.SetDGdegree(k);


            double radiusOP = 100; // operating point -> Re = radiusOp / Lstar
            double viscosity = 0.01; // kinematic viscosity
            double omega = 0.01; // rotation rate
            double Lstar = Math.Sqrt(viscosity / omega);

            double z = 10;
            double zstar = z * Lstar;


            // physical parameters
            double density = 1.0;
            C.PhysicalParameters.rho_A = density;
            C.PhysicalParameters.mu_A = density * viscosity;

            C.PhysicalParameters.IncludeConvection = true;


            // set grid
            double l_radial = zstar / 10.0;
            double l_azimuthal = zstar / 5.0;
            double h_axial = zstar;

            int res_global = 2;
            int res_radial = 1 * res_global;
            int res_azimuthal = 2 * res_global;
            int res_axial = 10 * res_global;

            Grid3D grd = RotatingDiskSector_Linearized(radiusOP, l_radial, l_azimuthal, h_axial, res_azimuthal, res_radial, res_axial);
            C.GridFunc = delegate() { return grd; };


            //C.UseRotInertForceTerms = true;
            //C.AngularVelocity = new double[] { 0.0, 0.0, omega };
            //C.UseCylindricalCoords = true;


            // boundary condition
            Formula RotatingDisk = new Formula("VelY", false,
                "double VelY(double[] X) { double radiusOP = 100.0; double omega = 0.01; return (X[1] + radiusOP) * omega; } ");
            C.AddBoundaryValue("wall_rotatingDisk", "VelocityX#A", RotatingDisk);


            double velW_top = -0.884473;
            double velWstar_top = velW_top * Math.Sqrt(viscosity * omega);
            Formula VelocityZ_top = new Formula("VelZ", false, "double VelZ(double[] X) { return -0.00884473; } ");

            C.AddBoundaryValue("velocity_inlet_top", "VelocityZ#A", VelocityZ_top);

            var vonKarman_sol = new RotatingDiskBoundaryLayerConditions();

            C.AddBoundaryValue("velocity_inlet_front", "VelocityX#A", vonKarman_sol);
            C.AddBoundaryValue("velocity_inlet_front", "VelocityY#A", vonKarman_sol);
            C.AddBoundaryValue("velocity_inlet_front", "VelocityZ#A", vonKarman_sol);

            C.AddBoundaryValue("velocity_inlet_back", "VelocityX#A", vonKarman_sol);
            C.AddBoundaryValue("velocity_inlet_back", "VelocityY#A", vonKarman_sol);
            C.AddBoundaryValue("velocity_inlet_back", "VelocityZ#A", vonKarman_sol);

            C.AddBoundaryValue("velocity_inlet_upstream", "VelocityX#A", vonKarman_sol);
            C.AddBoundaryValue("velocity_inlet_upstream", "VelocityY#A", vonKarman_sol);
            C.AddBoundaryValue("velocity_inlet_upstream", "VelocityZ#A", vonKarman_sol);

            C.AddBoundaryValue("velocity_inlet_downstream", "VelocityX#A", vonKarman_sol);
            C.AddBoundaryValue("velocity_inlet_downstream", "VelocityY#A", vonKarman_sol);
            C.AddBoundaryValue("velocity_inlet_downstream", "VelocityZ#A", vonKarman_sol);


            C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            // C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            // C.TimeSteppingScheme = TimeSteppingScheme.BDF3;
            // C.dtFixed = Math.PI * 1.0e-3;
            // C.NoOfTimesteps = 2000;

            //{
            //    C.AdaptiveMeshRefinement = true;
            //    int AMRlevel = 2;
            //    C.activeAMRlevelIndicators.Add(new AMRLevelIndicatorLibrary.AMRonBoundary(new byte[] { 1 }) { maxRefinementLevel = AMRlevel });
            //    C.AMR_startUpSweeps = AMRlevel;
            //}

            //Func<double[], double> PhiFunc = X => -1.0;
            //C.InitialValues_Evaluators.Add("Phi", PhiFunc);

            C.Option_LevelSetEvolution = LevelSetEvolution.None;


            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
            C.NonLinearSolver.ConvergenceCriterion = 1e-9;

            return C;
        }


        static public Grid3D RotatingDiskSector_CartesianCutOut(double radiusOP, double l_radial, double l_azimuthal, double h_axial, int res_radial, int res_azimuthal, int res_axial) {

            double[] xNodes = GenericBlas.Linspace(radiusOP - (l_radial / 2.0), radiusOP + (l_radial / 2.0), res_radial + 1);    // radial direction
            double[] yNodes = GenericBlas.Linspace(-l_azimuthal / 2.0, l_azimuthal / 2.0, res_azimuthal + 1);    // azimuthal direction
            // double[] yNodes = GenericBlas.Linspace(-l_azimuthal / 3.0, (2.0 * l_azimuthal) / 3.0, res_azimuthal + 1);    // azimuthal direction
            double[] zNodes = GenericBlas.Linspace(0.0, h_axial, res_axial + 1);    // axial direction

            var grd = Grid3D.Cartesian3DGrid(xNodes, yNodes, zNodes, periodicY: false);
            grd.Name = $"RotatingDiskSector3D_CartesianCutOut_{res_radial}x{res_azimuthal}x{res_axial}";

            grd.EdgeTagNames.Add(1, "velocity_inlet_rotatingDisk");
            grd.EdgeTagNames.Add(2, "velocity_inlet_top");
            grd.EdgeTagNames.Add(3, "velocity_inlet_back");
            grd.EdgeTagNames.Add(4, "velocity_inlet_front");
            grd.EdgeTagNames.Add(5, "velocity_inlet_upstream");
            grd.EdgeTagNames.Add(6, "pressure_outlet_downstream");

            grd.DefineEdgeTags(delegate (Vector X) {
                byte et = 0;
                if (X.z.Abs() <= 1e-8)
                    et = 1;
                if ((X.z - h_axial).Abs() <= 1e-8)
                    et = 2;
                if (((X.x - radiusOP) + (l_radial / 2.0)).Abs() <= 1e-8)
                    et = 3;
                if (((X.x - radiusOP) - (l_radial / 2.0)).Abs() <= 1e-8)
                    et = 4;
                if ((X.y + (l_azimuthal / 2.0)).Abs() <= 1e-8)
                    et = 5;
                if ((X.y - (l_azimuthal / 2.0)).Abs() <= 1e-8)
                    et = 6;
                // if ((X.y + (l_azimuthal / 3.0)).Abs() <= 1e-8)
                //     et = 5;
                // if ((X.y - (2.0 * l_azimuthal / 3.0)).Abs() <= 1e-8)
                //     et = 6;

                return et;
            });

            return grd;
        }


        static public GridCommons RotatingDiskSector_CartesianCutOutExtended(double radiusOP, double l_radial, double l_upstream, double l_downstream, double h_axial, int res_radial, int res_azimuthal, int res_axial,
                                                                    int radialExtLayer = 0, bool bothRadialDir = false, int azimuthalExtLayer = 0, bool bothAzimuthalDir = false, int axialExtLayer = 0)
        {

            int maxLayer = (new int[] { radialExtLayer, azimuthalExtLayer, axialExtLayer }).Max();
            GridBox[] grdBoxes = new GridBox[maxLayer + 1];

            grdBoxes[maxLayer] = new GridBox( new double[] { radiusOP - (l_radial / 2.0), -l_upstream, 0.0 },
                new double[] { radiusOP + (l_radial / 2.0), l_downstream , h_axial},
                res_radial, res_azimuthal, res_axial );

            double dx_radial = l_radial / res_radial;
            double dx_azimuthal = (l_upstream + l_downstream) / res_azimuthal;
            double dx_axial = h_axial / res_axial;

            for (int l = 1; l <= maxLayer; l++) {
                var embBB = grdBoxes[(maxLayer + 1) - l].boundingBox;
                double xMin = embBB.Min[0];
                double xMax = embBB.Max[0];
                double yMin = embBB.Min[1];
                double yMax = embBB.Max[1];
                double zMin = 0.0;
                double zMax = embBB.Max[2];

                var cellsBB = grdBoxes[(maxLayer + 1) - l].numOfCells;

                int resExtLayer_radial = (cellsBB[0] / 2);
                if (l <= radialExtLayer) {
                    if (bothRadialDir) {
                        resExtLayer_radial += 2;
                        xMax += dx_radial * 2.0.Pow(l);
                        xMin -= dx_radial * 2.0.Pow(l);
                        // ensure 2:1
                        int rEL_radMod = (resExtLayer_radial / 2) % 2;
                        if (rEL_radMod == 1) {
                            resExtLayer_radial += 2;
                            xMax += dx_radial * 2.0.Pow(l);
                            xMin -= dx_radial * 2.0.Pow(l);
                        }
                    } else {
                        resExtLayer_radial += 1;
                        xMax += dx_radial * 2.0.Pow(l);
                        // ensure 2:1
                        int rEL_radMod = resExtLayer_radial % 2;
                        if (rEL_radMod == 1) {
                            resExtLayer_radial += 1;
                            xMax += dx_radial * 2.0.Pow(l);
                        }
                    }
                }

                int resExtLayer_azimuthal = (cellsBB[1] / 2);
                if (l <= azimuthalExtLayer) {
                    if (bothAzimuthalDir) {
                        resExtLayer_azimuthal += 2;
                        yMax += dx_azimuthal * 2.0.Pow(l);
                        yMin -= dx_azimuthal * 2.0.Pow(l);
                        // ensure 2:1
                        int rEL_aziMod = (resExtLayer_azimuthal / 2) % 2;
                        if (rEL_aziMod == 1) {
                            resExtLayer_azimuthal += 2;
                            yMax += dx_azimuthal * 2.0.Pow(l);
                            yMin -= dx_azimuthal * 2.0.Pow(l);
                        }
                    }
                    else {
                        resExtLayer_azimuthal += 1;
                        yMax += dx_azimuthal * 2.0.Pow(l);
                        // ensure 2:1
                        int rEL_aziMod = resExtLayer_azimuthal % 2;
                        if (rEL_aziMod == 1) {
                            resExtLayer_azimuthal += 1;
                            yMax += dx_azimuthal * 2.0.Pow(l);
                        }
                    }
                }

                int resExtLayer_axial = (cellsBB[2] / 2);
                if (l <= axialExtLayer) {
                    resExtLayer_axial += 1;
                    zMax += dx_axial * 2.0.Pow(l); 
                    // ensure 2:1
                    int rEL_axiMod = resExtLayer_axial % 2;
                    if (rEL_axiMod == 1) {
                        resExtLayer_axial += 1;
                        zMax += dx_axial * 2.0.Pow(l);  
                    }
                }

                grdBoxes[maxLayer - l] = new GridBox( new double[] { xMin, yMin, zMin },
                    new double[] { xMax, yMax, zMax },
                    resExtLayer_radial, resExtLayer_azimuthal, resExtLayer_axial);
            }


            var grd = Grid3D.HangingNodes3D(false, false, false, grdBoxes);
            grd.Name = $"RotatingDiskSector3D_CartesianCutOutExtended_{res_radial}x{res_azimuthal}x{res_axial}";

            grd.EdgeTagNames.Add(1, "velocity_inlet_rotatingDisk");
            grd.EdgeTagNames.Add(2, "dong_outflow_top");
            grd.EdgeTagNames.Add(3, "velocity_inlet_back");
            grd.EdgeTagNames.Add(4, "Dong_OutFlow_front");
            grd.EdgeTagNames.Add(5, "velocity_inlet_upstream");
            grd.EdgeTagNames.Add(6, "dong_outflow_downstream");

            grd.DefineEdgeTags(delegate (Vector X) {
                byte et = 0;
                if ((X.z - grdBoxes[0].boundingBox.Min[2]).Abs() <= 1e-8)
                    et = 1;
                if ((X.z - grdBoxes[0].boundingBox.Max[2]).Abs() <= 1e-8)
                    et = 2;
                if ((X.x - grdBoxes[0].boundingBox.Min[0]).Abs() <= 1e-8)
                    et = 3;
                if ((X.x - grdBoxes[0].boundingBox.Max[0]).Abs() <= 1e-8)
                    et = 4;
                if ((X.y - grdBoxes[0].boundingBox.Min[1]).Abs() <= 1e-8)
                    et = 5;
                if ((X.y - grdBoxes[0].boundingBox.Max[1]).Abs() <= 1e-8)
                    et = 6;

                return et;
            });

            return grd;
        }


        static public GridCommons RotatingDiskSector2D_CartesianCutOutExtended(double radiusOP, double l_upstream, double l_downstream, double h_axial, int res_azimuthal, int res_axial,
                                                                    int azimuthalExtLayer = 0, bool bothAzimuthalDir = false, int axialExtLayer = 0)
        {

            int maxLayer = (new int[] { azimuthalExtLayer, axialExtLayer }).Max();
            GridBox[] grdBoxes = new GridBox[maxLayer + 1];

            grdBoxes[maxLayer] = new GridBox(new double[] { -l_upstream, 0.0 },
                new double[] { l_downstream, h_axial },
                res_azimuthal, res_axial);

            double dx_azimuthal = (l_upstream + l_downstream) / res_azimuthal;
            double dx_axial = h_axial / res_axial;

            for (int l = 1; l <= maxLayer; l++)
            {
                var embBB = grdBoxes[(maxLayer + 1) - l].boundingBox;
                double xMin = embBB.Min[0];
                double xMax = embBB.Max[0];
                double yMin = 0.0;
                double yMax = embBB.Max[1];

                var cellsBB = grdBoxes[(maxLayer + 1) - l].numOfCells;


                int resExtLayer_azimuthal = (cellsBB[0] / 2);
                if (l <= azimuthalExtLayer)
                {
                    if (bothAzimuthalDir)
                    {
                        resExtLayer_azimuthal += 2;
                        xMax += dx_azimuthal * 2.0.Pow(l);
                        xMin -= dx_azimuthal * 2.0.Pow(l);
                        // ensure 2:1
                        int rEL_aziMod = (resExtLayer_azimuthal / 2) % 2;
                        if (rEL_aziMod == 1)
                        {
                            resExtLayer_azimuthal += 2;
                            xMax += dx_azimuthal * 2.0.Pow(l);
                            xMin -= dx_azimuthal * 2.0.Pow(l);
                        }
                    }
                    else
                    {
                        resExtLayer_azimuthal += 1;
                        xMax += dx_azimuthal * 2.0.Pow(l);
                        // ensure 2:1
                        int rEL_aziMod = resExtLayer_azimuthal % 2;
                        if (rEL_aziMod == 1)
                        {
                            resExtLayer_azimuthal += 1;
                            xMax += dx_azimuthal * 2.0.Pow(l);
                        }
                    }
                }

                int resExtLayer_axial = (cellsBB[1] / 2);
                if (l <= axialExtLayer)
                {
                    resExtLayer_axial += 1;
                    yMax += dx_axial * 2.0.Pow(l);
                    // ensure 2:1
                    int rEL_axiMod = resExtLayer_axial % 2;
                    if (rEL_axiMod == 1)
                    {
                        resExtLayer_axial += 1;
                        yMax += dx_axial * 2.0.Pow(l);
                    }
                }

                grdBoxes[maxLayer - l] = new GridBox(new double[] { xMin, yMin },
                    new double[] { xMax, yMax },
                    resExtLayer_azimuthal, resExtLayer_axial);
            }


            var grd = Grid2D.HangingNodes2D(false, false, grdBoxes);
            grd.Name = $"RotatingDiskSector2D_CartesianCutOutExtended_{res_azimuthal}x{res_axial}";

            grd.EdgeTagNames.Add(1, "velocity_inlet_rotatingDisk");
            grd.EdgeTagNames.Add(2, "dong_outflow_top");
            grd.EdgeTagNames.Add(3, "velocity_inlet_upstream");
            grd.EdgeTagNames.Add(4, "dong_outflow_downstream");

            grd.DefineEdgeTags(delegate (Vector X) {
                byte et = 0;
                if ((X.y - grdBoxes[0].boundingBox.Min[1]).Abs() <= 1e-8)
                    et = 1;
                if ((X.y - grdBoxes[0].boundingBox.Max[1]).Abs() <= 1e-8)
                    et = 2;
                if ((X.x - grdBoxes[0].boundingBox.Min[0]).Abs() <= 1e-8)
                    et = 3;
                if ((X.x - grdBoxes[0].boundingBox.Max[0]).Abs() <= 1e-8)
                    et = 4;

                return et;
            });

            return grd;
        }


        public static XNSE_Control DropletRebound(int k = 3, int gridRes = 8) {

            var C = new XNSE_Control();

            C.SetDGdegree(k);


            // physical parameters
            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.mu_A = 1.0;

            C.PhysicalParameters.rho_B = 1.0;
            C.PhysicalParameters.mu_B = 0.01;

            C.PhysicalParameters.Sigma = 72.9e-3;

            C.PhysicalParameters.IncludeConvection = true;


            // set grid
            double zTopStar = 10;

            double l_radial = zTopStar;
            double l_upstream = zTopStar / 2.0;
            double l_downstream = 3.0 * zTopStar / 2.0;
            double h_axial = zTopStar;

            int res_global = gridRes;
            int res_radial = 1 * res_global;
            int res_azimuthal = (4 * res_global) / 2;
            int res_axial = 1 * res_global;

            double radiusOP = 100;
            //Grid3D grd = RotatingDiskSector_CartesianCutOut(radiusOP, l_radial, l_azimuthal, h_axial, res_radial, res_azimuthal, res_axial);
            GridCommons grd = RotatingDiskSector_CartesianCutOutExtended(radiusOP, l_radial, l_upstream, l_downstream, h_axial, 
                res_radial, res_azimuthal, res_axial, 2, true, 3, true, 2);
            C.GridFunc = delegate () { return grd; };


            // boundary conditions
            //C.AddBoundaryValue("velocity_inlet_rotatingDisk", "VelocityX#B", RotatingDiskVelocityX);
            //C.AddBoundaryValue("velocity_inlet_rotatingDisk", "VelocityY#B", RotatingDiskVelocityY);

            //C.AddBoundaryValue("velocity_inlet_top", "VelocityX#B", vonKarmanHAM_velX);
            //C.AddBoundaryValue("velocity_inlet_top", "VelocityY#B", vonKarmanHAM_velY);
            //C.AddBoundaryValue("velocity_inlet_top", "VelocityZ#B", vonKarmanHAM_velZ);

            //C.AddBoundaryValue("velocity_inlet_front", "VelocityX#B", vonKarmanHAM_velX);
            //C.AddBoundaryValue("velocity_inlet_front", "VelocityY#B", vonKarmanHAM_velY);
            //C.AddBoundaryValue("velocity_inlet_front", "VelocityZ#B", vonKarmanHAM_velZ);

            //C.AddBoundaryValue("velocity_inlet_back", "VelocityX#B", vonKarmanHAM_velX);
            //C.AddBoundaryValue("velocity_inlet_back", "VelocityY#B", vonKarmanHAM_velY);
            //C.AddBoundaryValue("velocity_inlet_back", "VelocityZ#B", vonKarmanHAM_velZ);

            //C.AddBoundaryValue("velocity_inlet_upstream", "VelocityX#B", vonKarmanHAM_velX);
            //C.AddBoundaryValue("velocity_inlet_upstream", "VelocityY#B", vonKarmanHAM_velY);
            //C.AddBoundaryValue("velocity_inlet_upstream", "VelocityZ#B", vonKarmanHAM_velZ);

            // C.AddBoundaryValue("velocity_inlet_downstream", "VelocityX#B", vonKarman_velX);
            // C.AddBoundaryValue("velocity_inlet_downstream", "VelocityY#B", vonKarman_velY);
            // C.AddBoundaryValue("velocity_inlet_downstream", "VelocityZ#B", vonKarman_velZ);
            // C.AddBoundaryValue("pressure_dirichlet_downstream", "Pressure#B", vonKarman_P);


            // initial conditions
            //C.AddInitialValue("VelocityX#B", vonKarmanHAM_velX);
            //C.AddInitialValue("VelocityY#B", vonKarmanHAM_velY);
            //C.AddInitialValue("VelocityZ#B", vonKarmanHAM_velZ);
            // C.AddInitialValue("Pressure#B", vonKarman_P);

            Formula PhiFunc = new Formula(
                "Phi",
                false,
                "double Phi(double[] X) { " +
                "double radiusOP = 100;" +
                "double radiusDrop = 1;" +
                "double initHeight = 7.5;" +
                "return Math.Sqrt((X[0] - radiusOP).Pow2() + X[1].Pow2() + (X[2] - initHeight).Pow2()) - radiusDrop; } "
            );
            //C.AddInitialValue("Phi", PhiFunc);
            C.AddInitialValue("Phi", new Formula("X => -1.0", false));

            // C.AddInitialValue("VelocityZ#A", InitVelocity);

            //C.UseRampUpForInitialVelocityField = true;
            //C.rampUpValues_VelocityZ = new double[] { 0.0, -0.112, -0.56, -1.12 };
            //C.AddInitialRampUpValue("VelocityZ#A", new ConstantValue(0.0));
            //C.AddInitialRampUpValue("VelocityZ#A", new ConstantValue(-0.112));
            //C.AddInitialRampUpValue("VelocityZ#A", new ConstantValue(-0.56));
            //C.AddInitialRampUpValue("VelocityZ#A", new ConstantValue(-1.12));


            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            C.LSContiProjectionMethod = ContinuityProjectionOption.ConstrainedDG;


            C.SkipSolveAndEvaluateResidual = true;

            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
            C.NonLinearSolver.ConvergenceCriterion = 1e-9;
            C.NonLinearSolver.MaxSolverIterations = 20;

            // C.LinearSolver = LinearSolverCode.exp_Kcycle_schwarz.GetConfig();


            //C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            // C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;
            // C.Option_LevelSetEvolution = LevelSetEvolution.StokesExtension;
            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.dtFixed = 1.0e-6;
            C.NoOfTimesteps = 1;

            {
                C.AdaptiveMeshRefinement = true;
                int AMRlevel = 2;
                //C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = AMRlevel });
                C.activeAMRlevelIndicators.Add(new AMRLevelIndicatorLibrary.AMRonBoundaryByResolution(new byte[] { 1 }, (2.0 * h_axial / res_axial) * 1.5) { maxRefinementLevel = AMRlevel });
                C.AMR_startUpSweeps = AMRlevel + 1;
            }

            return C;
        }


        public static XNSE_Control DropletRebound2D(int k = 3, int gridRes = 8)
        {

            var C = new XNSE_Control();

            C.SetDGdegree(k);


            // physical parameters
            C.PhysicalParameters.rho_A = 1.0;
            C.PhysicalParameters.mu_A = 1.0;

            C.PhysicalParameters.rho_B = 1.0;
            C.PhysicalParameters.mu_B = 0.01;

            C.PhysicalParameters.Sigma = 72.9e-3;

            C.PhysicalParameters.IncludeConvection = true;


            // set grid
            double zTopStar = 10;

            double l_upstream = zTopStar / 2.0;
            double l_downstream = 3.0 * zTopStar / 2.0;
            double h_axial = zTopStar;

            int res_global = gridRes;
            int res_azimuthal = (4 * res_global) / 2;
            int res_axial = 1 * res_global;

            double radiusOP = 100;
            //Grid3D grd = RotatingDiskSector_CartesianCutOut(radiusOP, l_radial, l_azimuthal, h_axial, res_radial, res_azimuthal, res_axial);
            GridCommons grd = RotatingDiskSector2D_CartesianCutOutExtended(radiusOP, l_upstream, l_downstream, h_axial,
               res_azimuthal, res_axial, 3, true, 2);
            C.GridFunc = delegate () { return grd; };


            // boundary conditions
            //C.AddBoundaryValue("velocity_inlet_rotatingDisk", "VelocityX#B", RotatingDiskVelocityX);
            //C.AddBoundaryValue("velocity_inlet_rotatingDisk", "VelocityY#B", RotatingDiskVelocityY);

            //C.AddBoundaryValue("velocity_inlet_top", "VelocityX#B", vonKarmanHAM_velX);
            //C.AddBoundaryValue("velocity_inlet_top", "VelocityY#B", vonKarmanHAM_velY);
            //C.AddBoundaryValue("velocity_inlet_top", "VelocityZ#B", vonKarmanHAM_velZ);

            //C.AddBoundaryValue("velocity_inlet_front", "VelocityX#B", vonKarmanHAM_velX);
            //C.AddBoundaryValue("velocity_inlet_front", "VelocityY#B", vonKarmanHAM_velY);
            //C.AddBoundaryValue("velocity_inlet_front", "VelocityZ#B", vonKarmanHAM_velZ);

            //C.AddBoundaryValue("velocity_inlet_back", "VelocityX#B", vonKarmanHAM_velX);
            //C.AddBoundaryValue("velocity_inlet_back", "VelocityY#B", vonKarmanHAM_velY);
            //C.AddBoundaryValue("velocity_inlet_back", "VelocityZ#B", vonKarmanHAM_velZ);

            //C.AddBoundaryValue("velocity_inlet_upstream", "VelocityX#B", vonKarmanHAM_velX);
            //C.AddBoundaryValue("velocity_inlet_upstream", "VelocityY#B", vonKarmanHAM_velY);
            //C.AddBoundaryValue("velocity_inlet_upstream", "VelocityZ#B", vonKarmanHAM_velZ);

            // C.AddBoundaryValue("velocity_inlet_downstream", "VelocityX#B", vonKarman_velX);
            // C.AddBoundaryValue("velocity_inlet_downstream", "VelocityY#B", vonKarman_velY);
            // C.AddBoundaryValue("velocity_inlet_downstream", "VelocityZ#B", vonKarman_velZ);
            // C.AddBoundaryValue("pressure_dirichlet_downstream", "Pressure#B", vonKarman_P);


            // initial conditions
            //C.AddInitialValue("VelocityX#B", vonKarmanHAM_velX);
            //C.AddInitialValue("VelocityY#B", vonKarmanHAM_velY);
            //C.AddInitialValue("VelocityZ#B", vonKarmanHAM_velZ);
            // C.AddInitialValue("Pressure#B", vonKarman_P);

            Formula PhiFunc = new Formula(
                "Phi",
                false,
                "double Phi(double[] X) { " +
                "double radiusOP = 100;" +
                "double radiusDrop = 1;" +
                "double initHeight = 7.5;" +
                "return Math.Sqrt(X[0].Pow2() + (X[1] - initHeight).Pow2()) - radiusDrop; } "
            );
            //C.AddInitialValue("Phi", PhiFunc);
            C.AddInitialValue("Phi", new Formula("X => -1.0", false));

            // C.AddInitialValue("VelocityZ#A", InitVelocity);

            //C.UseRampUpForInitialVelocityField = true;
            //C.rampUpValues_VelocityZ = new double[] { 0.0, -0.112, -0.56, -1.12 };
            //C.AddInitialRampUpValue("VelocityZ#A", new ConstantValue(0.0));
            //C.AddInitialRampUpValue("VelocityZ#A", new ConstantValue(-0.112));
            //C.AddInitialRampUpValue("VelocityZ#A", new ConstantValue(-0.56));
            //C.AddInitialRampUpValue("VelocityZ#A", new ConstantValue(-1.12));


            C.AdvancedDiscretizationOptions.SST_isotropicMode = SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;
            C.LSContiProjectionMethod = ContinuityProjectionOption.ConstrainedDG;


            C.SkipSolveAndEvaluateResidual = true;

            C.NonLinearSolver.SolverCode = NonLinearSolverCode.Picard;
            C.NonLinearSolver.ConvergenceCriterion = 1e-9;
            C.NonLinearSolver.MaxSolverIterations = 20;

            // C.LinearSolver = LinearSolverCode.exp_Kcycle_schwarz.GetConfig();


            //C.TimesteppingMode = AppControl._TimesteppingMode.Steady;
            C.Timestepper_LevelSetHandling = LevelSetHandling.None;
            C.Option_LevelSetEvolution = LevelSetEvolution.None;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            // C.Timestepper_LevelSetHandling = LevelSetHandling.Coupled_Once;
            // C.Option_LevelSetEvolution = LevelSetEvolution.StokesExtension;
            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;
            C.dtFixed = 1.0e-6;
            C.NoOfTimesteps = 1;

            {
                C.AdaptiveMeshRefinement = true;
                int AMRlevel = 3;
                //C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = AMRlevel });
                C.activeAMRlevelIndicators.Add(new AMRLevelIndicatorLibrary.AMRonBoundaryByResolution(new byte[] { 1 }, (h_axial/res_axial)*1.5) { maxRefinementLevel = AMRlevel });
                C.AMR_startUpSweeps = AMRlevel + 1;
            }

            return C;
        }


        public static XNSE_Control DropletReboundGauthier_Restart() {

            string _DbPath = @"\\dc3\userspace\smuda\hpccluster\DropletRebound_Gauthier";
            var rstID = new Guid("636c0c52-ae2b-46dd-9547-3a18d8f6d1a0");
            int RestartTimestep = 1;

            string sessionDir = Path.Combine(_DbPath, StandardFsDriver.SessionsDir, rstID.ToString());
            string path_obj = Path.Combine(sessionDir, "Control-obj.txt");

            string ctrlfileContent = File.ReadAllText(path_obj);

            var ctrl = (XNSE_Control)AppControl.Deserialize(ctrlfileContent).CloneAs(); ;

            ctrl.InitialValues.Clear();
            ctrl.InitialValues_Evaluators.Clear();

            ctrl.RestartInfo = new Tuple<Guid, TimestepNumber>(rstID, RestartTimestep);
            ctrl.GridGuid = Guid.Empty;
            ctrl.DbPath = _DbPath;
            ctrl.savetodb = false;

            ctrl.SessionName = "DropletReboundGauthier_restart_DebugRun";

            return ctrl;

        }

        /// <summary>
        /// cases for testing on jenkins
        /// DropletReboundGauthier_8x12x8AMR1_k2_ReI4_restart3 - 4x4 multicore simulation
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control DropletReboundGauthier_8x12x8AMR1_k2_ReI4_Restart(int numTimeSteps = 1) {

            string _DbPath = @"\\fdygitrunner\ValidationTests\databases\DropletRebound_Gauthier";
            var rstID = new Guid("2cd1d33d-5099-414d-b532-95096879a027");
            int RestartTimestep = 77;

            string sessionDir = Path.Combine(_DbPath, StandardFsDriver.SessionsDir, rstID.ToString());
            string path_obj = Path.Combine(sessionDir, "Control-obj.txt");

            string ctrlfileContent = File.ReadAllText(path_obj);

            var ctrl = (XNSE_Control)AppControl.Deserialize(ctrlfileContent).CloneAs(); ;

            ctrl.InitialValues.Clear();
            ctrl.InitialValues_Evaluators.Clear();

            ctrl.RestartInfo = new Tuple<Guid, TimestepNumber>(rstID, RestartTimestep);
            ctrl.ReInitTimestepIndex = RestartTimestep;
            ctrl.NoOfTimesteps = numTimeSteps;

            ctrl.Option_LevelSetEvolution = LevelSetEvolution.StokesExtension;

            ctrl.DbPath = _DbPath;
            ctrl.savetodb = false;

            ctrl.SessionName = "DropletReboundGauthier_8x12x8AMR1_k2_ReI4_restart4";

            return ctrl;

        }


        /// <summary>
        /// cases for testing on jenkins
        /// DropletReboundGauthier_8x8x8AMR1_k3_ReI4_restart5 - 4x4 multicore simulation
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control DropletReboundGauthier_8x8x8AMR1_k3_ReI4_Restart(int numTimeSteps = 1) {

            string _DbPath = @"\\fdygitrunner\ValidationTests\databases\DropletRebound_Gauthier";
            var rstID = new Guid("4048fd03-ba72-488c-a363-66df1b140d15");
            int RestartTimestep = -1;

            string sessionDir = Path.Combine(_DbPath, StandardFsDriver.SessionsDir, rstID.ToString());
            string path_obj = Path.Combine(sessionDir, "Control-obj.txt");

            string ctrlfileContent = File.ReadAllText(path_obj);

            var ctrl = (XNSE_Control)AppControl.Deserialize(ctrlfileContent).CloneAs(); ;

            ctrl.InitialValues.Clear();
            ctrl.InitialValues_Evaluators.Clear();

            ctrl.RestartInfo = new Tuple<Guid, TimestepNumber>(rstID, RestartTimestep);
            ctrl.ReInitTimestepIndex = RestartTimestep;
            ctrl.NoOfTimesteps = numTimeSteps;

            ctrl.Option_LevelSetEvolution = LevelSetEvolution.StokesExtension;

            ctrl.DbPath = _DbPath;
            ctrl.savetodb = false;

            ctrl.SessionName = "DropletReboundGauthier_8x8x8AMR1_k3_ReI4_restart6";

            return ctrl;

        }

        /// <summary>
        /// cases for testing on jenkins
        /// case 1: DropletReboundGauthier_8x8x8AMR0_k3_ReI4_restart1 - single core simulation
        /// case 2: DropletReboundGauthier_8x8x8AMR0_k3_ReI4_2cores_restart1 - 2 core simulation
        /// </summary>
        /// <returns></returns>
        public static XNSE_Control DropletReboundGauthier_8x8x8AMR0_k3_ReI4_Restart(int caseNum = 1, int numTimeSteps = 1) {

            string _DbPath = @"\\fdygitrunner\ValidationTests\databases\DropletRebound_Gauthier";

            Guid rstID;
            int RestartTimestep;
            string RestartName;
            switch (caseNum) {
                case 1: {
                        rstID = new Guid("57e8f3ff-209a-4b66-bfcc-049279c44627");   // single core simulation
                        RestartTimestep = 66;
                        RestartName = "DropletReboundGauthier_8x8x8AMR0_k3_ReI4_restart2_DongBC";
                        break;
                    }
                case 2: {
                        rstID = new Guid("1a8cfde2-27a6-4c58-b4ed-0607a243731a");   // 2 core simulation
                        RestartTimestep = 77;
                        RestartName = "DropletReboundGauthier_8x8x8AMR0_k3_ReI4_2cores_restart2";
                        break;
                    }
                default: {
                        throw new ArgumentOutOfRangeException();
                    }
            }

            string sessionDir = Path.Combine(_DbPath, StandardFsDriver.SessionsDir, rstID.ToString());
            string path_obj = Path.Combine(sessionDir, "Control-obj.txt");

            string ctrlfileContent = File.ReadAllText(path_obj);

            var ctrl = (XNSE_Control)AppControl.Deserialize(ctrlfileContent).CloneAs();

            ctrl.InitialValues.Clear();
            ctrl.InitialValues_Evaluators.Clear();

            ctrl.RestartInfo = new Tuple<Guid, TimestepNumber>(rstID, RestartTimestep);
            ctrl.ReInitTimestepIndex = RestartTimestep;
            ctrl.NoOfTimesteps = numTimeSteps;

            ctrl.Option_LevelSetEvolution = LevelSetEvolution.StokesExtension;

            ctrl.NonLinearSolver.Globalization = Newton.GlobalizationOption.LineSearch;

            ctrl.ChangeBoundaryCondition("pressure_outlet_top", "Dong_OutFlow_top");
            ctrl.ChangeBoundaryCondition("pressure_outlet_front", "Dong_OutFlow_front");
            ctrl.ChangeBoundaryCondition("pressure_outlet_downstream", "Dong_OutFlow_downstream");

            ctrl.DbPath = _DbPath;
            ctrl.savetodb = false;

            ctrl.SessionName = RestartName;

            return ctrl;

        }


        public static XNSE_Control DropletImpactTest_hydrophilicSurface(int k = 3, int numCells = 8) {

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            C.DbPath = null;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/DropletImpact";
            C.ProjectDescription = "droplet impact on solid wall (hydrophilic)";

            #endregion


            // DG degrees
            // ==========
            C.SetDGdegree(k);


            // physical parameters
            // ===================
            #region physics

            C.Tags.Add("Water/Air at 1atm");
            C.PhysicalParameters.rho_A = 1e3;
            C.PhysicalParameters.rho_B = 1.0;
            double mu_scl = 1.0;
            C.PhysicalParameters.mu_A = 1e-3 * mu_scl;
            C.PhysicalParameters.mu_B = 1.48e-5 * mu_scl;
            double sigma = 0.07;
            C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.betaS_A = 0.0;
            C.PhysicalParameters.betaS_B = 0.0;

            C.PhysicalParameters.betaL = 0.0;
            C.PhysicalParameters.theta_e = (60.0 / 180.0) * Math.PI;

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double domainSize = 3.0e-3;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, domainSize, numCells + 1);
                double[] Ynodes = GenericBlas.Linspace(0, domainSize, numCells + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "slipsymmetry_left");
                grd.EdgeTagNames.Add(2, "pressure_outlet_right");
                grd.EdgeTagNames.Add(3, "navierslip_linear_surface");
                grd.EdgeTagNames.Add(4, "pressure_outlet_upper");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] - domainSize) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[1] - domainSize) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // initial values
            // ==============
            #region init

            double radius = 1.0e-3;
            double offset = 0.0; 

            Func<double[], double> PhiFunc = (X => ((X[0]).Pow2() + (X[1] - (radius + offset)).Pow2()).Sqrt() - radius);

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);


            C.InitialValues_Evaluators.Add("Pressure#A", X => sigma / radius);
            C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);

            double U0 = 0.2;
            C.InitialValues_Evaluators.Add("VelocityY#A", X => -U0);

            double g = 9.81;
            C.InitialValues_Evaluators.Add("GravityY#A", X => -g);
            C.InitialValues_Evaluators.Add("GravityY#B", X => -g);

            #endregion


            // boundary conditions
            // ===================
            #region BC


            C.AddBoundaryValue("slipsymmetry_left");
            C.AddBoundaryValue("pressure_outlet_right");
            C.AddBoundaryValue("navierslip_linear_surface");
            C.AddBoundaryValue("pressure_outlet_upper");


            #endregion


            // misc. solver options
            // ====================
            #region solver


            C.Option_LevelSetEvolution = LevelSetEvolution.FastMarching;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.LevelSetTools.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            C.AdaptiveMeshRefinement = true;
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 1 });
            C.activeAMRlevelIndicators.Add(new AMRonBoundary(new byte[] { 3 }) { maxRefinementLevel = 1 });
            C.AMR_startUpSweeps = 1;


            #endregion


            // Timestepping
            // ============
            #region time


            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;

            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            double dt = 1e-6;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 1000;
            C.saveperiod = 1;

            #endregion


            return C;

        }


        public static XNSE_Control DropletImpactTest_hydrophobicSurface(int k = 3, int numCells = 8) {

            XNSE_Control C = new XNSE_Control();

            // basic database options
            // ======================
            #region db

            C.DbPath = null;
            C.savetodb = C.DbPath != null;
            C.ProjectName = "XNSE/DropletImpact";
            C.ProjectDescription = "droplet impact on solid wall (hydrophobic)";

            //C.PostprocessingModules.Add(new MovingContactLineLogging() { LogPeriod = 1, printToConsole = true });

            #endregion


            // DG degrees
            // ==========
            C.SetDGdegree(k);


            // physical parameters
            // ===================
            #region physics

            C.Tags.Add("Water/Air at 1atm");
            C.PhysicalParameters.rho_A = 1e3;
            C.PhysicalParameters.rho_B = 1.0;
            double mu_scl = 1.0;
            C.PhysicalParameters.mu_A = 1e-3 * mu_scl;
            C.PhysicalParameters.mu_B = 1.48e-5 * mu_scl;
            double sigma = 0.07;
            C.PhysicalParameters.Sigma = sigma;

            C.PhysicalParameters.betaS_A = 5.0;
            C.PhysicalParameters.betaS_B = 5.0;

            C.PhysicalParameters.betaL = 0.0;
            C.PhysicalParameters.theta_e = (150.0 / 180.0) * Math.PI;

            C.PhysicalParameters.IncludeConvection = false;
            C.PhysicalParameters.Material = true;

            #endregion


            // grid generation
            // ===============
            #region grid

            double domainSize = 3.15e-3;

            C.GridFunc = delegate () {
                double[] Xnodes = GenericBlas.Linspace(0, domainSize, numCells + 1);
                double[] Ynodes = GenericBlas.Linspace(0, domainSize, numCells + 1);
                var grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes);

                grd.EdgeTagNames.Add(1, "slipsymmetry_left");
                grd.EdgeTagNames.Add(2, "pressure_outlet_right");
                grd.EdgeTagNames.Add(3, "navierslip_linear_surface");
                grd.EdgeTagNames.Add(4, "pressure_outlet_upper");

                grd.DefineEdgeTags(delegate (double[] X) {
                    byte et = 0;
                    if (Math.Abs(X[0]) <= 1.0e-8)
                        et = 1;
                    if (Math.Abs(X[0] - domainSize) <= 1.0e-8)
                        et = 2;
                    if (Math.Abs(X[1]) <= 1.0e-8)
                        et = 3;
                    if (Math.Abs(X[1] - domainSize) <= 1.0e-8)
                        et = 4;

                    return et;
                });

                return grd;
            };

            #endregion


            // initial values
            // ==============
            #region init

            double radius = 1.05e-3;
            double offset = -1.0e-6;

            Func<double[], double> PhiFunc = (X => ((X[0]).Pow2() + (X[1] - (radius + offset)).Pow2()).Sqrt() - radius);
            //Func<double[], double> PhiFunc = (X => ((X[0]).Pow2() + (X[1] - (radius + offset)).Pow2()) - radius.Pow2());

            C.InitialValues_Evaluators.Add("Phi", PhiFunc);


            C.InitialValues_Evaluators.Add("Pressure#A", X => sigma / radius);
            C.InitialValues_Evaluators.Add("Pressure#B", X => 0.0);

            double U0 = 0.0;
            C.InitialValues_Evaluators.Add("VelocityY#A", X => -U0);

            double g = 9.81;
            C.InitialValues_Evaluators.Add("GravityY#A", X => -g);
            C.InitialValues_Evaluators.Add("GravityY#B", X => -g);

            #endregion


            // boundary conditions
            // ===================
            #region BC


            C.AddBoundaryValue("slipsymmetry_left");
            C.AddBoundaryValue("pressure_outlet_right");
            C.AddBoundaryValue("navierslip_linear_surface");
            C.AddBoundaryValue("pressure_outlet_upper");


            #endregion


            // misc. solver options
            // ====================
            #region solver


            C.Option_LevelSetEvolution = LevelSetEvolution.StokesExtension;

            C.AdvancedDiscretizationOptions.SST_isotropicMode = Solution.LevelSetTools.SurfaceStressTensor_IsotropicMode.LaplaceBeltrami_ContactLine;

            C.AdaptiveMeshRefinement = true;
            C.activeAMRlevelIndicators.Add(new AMRonNarrowband() { maxRefinementLevel = 2 });
            //C.activeAMRlevelIndicators.Add(new AMRonBoundary(new byte[] { 3 }) { maxRefinementLevel = 3 });
            C.AMR_startUpSweeps = 3;

            //C.ReInitPeriod = 1;
            //C.InitSignedDistance = true;  // Deprecated


            #endregion


            // Timestepping
            // ============
            #region time


            C.TimeSteppingScheme = TimeSteppingScheme.ImplicitEuler;

            C.Timestepper_LevelSetHandling = LevelSetHandling.LieSplitting;

            C.TimesteppingMode = AppControl._TimesteppingMode.Transient;
            double dt = 5e-6;
            C.dtMax = dt;
            C.dtMin = dt;
            C.Endtime = 1000;
            C.NoOfTimesteps = 1000;
            C.saveperiod = 1;

            #endregion


            return C;

        }

    }
}
