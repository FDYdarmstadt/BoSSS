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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.XdgTimestepping;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Runtime.Serialization;



namespace BoSSS.Application.FSI_Solver {
    [DataContract]
    [Serializable]
    public class FSI_Control : IBM_Solver.IBM_Control {

        /// <summary>
        /// ctor to remain compatible with old control files
        /// </summary>
        public FSI_Control() : base() {
            this.Particles = new List<Particle>();
        }

        /// <summary>
        /// ctor
        /// </summary>
        public FSI_Control(int degree, string projectName, string projectDescription = null, List<string> tags = null) : base() {
            this.Particles = new List<Particle>();
            ProjectName = projectName;
            SessionName = projectName;
            ProjectDescription = projectDescription;
            if (tags != null) {
                for (int i = 0; i < tags.Count(); i++) {
                    Tags.Add(tags[i]);
                }
            }
            SetDGdegree(degree);
        }

        /// <summary>
        /// Set true during restart.
        /// </summary>
        public bool IsRestart = false;

        /// <summary>
        /// List of the boundary values at the domain boundary.
        /// </summary>
        readonly List<string> m_BoundaryValues = new List<string>();

        /// <summary>
        /// Overall domain volume
        /// </summary>
        public double DomainVolume;

        /// <summary>
        /// The position of all boundaries, independent of boundary condition.
        /// </summary>
        public double[][] BoundaryPositionPerDimension;

        /// <summary>
        /// The position of all boundaries, only for walls.
        /// </summary>
        public double[][] WallPositionPerDimension;

        /// <summary>
        /// True for periodic walls
        /// </summary>
        public bool[] BoundaryIsPeriodic;

        /// <summary>
        /// Max grid length
        /// </summary>
        public double MaxGridLength;

        /// <summary>
        /// Min grid length
        /// </summary>
        public double MinGridLength;

        /// <summary>
        /// Setting <see cref="Solution.Control.AppControl.FieldOptions"/>
        /// </summary>
        public override void SetDGdegree(int k) {
            if (k < 1)
                throw new ArgumentOutOfRangeException("DG polynomial degree must be at least 1.");
            int k_phiDG = Math.Max(4, k);
            int k_phi = k_phiDG * 2;
            base.FieldOptions.Clear();
            this.AddFieldOption("Velocity*", k);
            this.AddFieldOption("Pressure", k - 1);
            this.AddFieldOption("PhiDG", k_phiDG);
            this.AddFieldOption("Phi", k_phi);
        }

        public void SetSaveOptions(string dataBasePath = null, int savePeriod = 1) {
            if (dataBasePath != null) {
                savetodb = true;
                DbPath = dataBasePath;
                saveperiod = savePeriod;
            }
            else
                savetodb = false;
        }
        
        public void SetAddaptiveMeshRefinement(int amrLevel) {
            if (amrLevel == 0)
                return;
            AdaptiveMeshRefinement = true;
            RefinementLevel = amrLevel;
            AMR_startUpSweeps = amrLevel;
        }

        public void SetBoundaries(List<string> boundaryValues) {
            if (boundaryValues.Count() > 4)
                throw new NotImplementedException("max 4 boundary values");
            for (int i = 0; i < boundaryValues.Count(); i++) {
                AddBoundaryValue(boundaryValues[i]);
                m_BoundaryValues.Add(boundaryValues[i]);
            }
        }

        public void SetGrid(double lengthX, double lengthY, double cellsPerUnitLength, bool periodicX = false, bool periodicY = false) {
            MaxGridLength = 1 / cellsPerUnitLength;
            BoundaryPositionPerDimension = new double[2][];
            WallPositionPerDimension = new double[2][];
            WallPositionPerDimension[0] = new double[2];
            WallPositionPerDimension[1] = new double[2];
            BoundaryIsPeriodic = new bool[2];
            BoundaryPositionPerDimension[0] = new double[] { -lengthX / 2, lengthX / 2 };
            BoundaryPositionPerDimension[1] = new double[] { -lengthY / 2, lengthY / 2 };
            BoundaryIsPeriodic[0] = periodicX;
            BoundaryIsPeriodic[1] = periodicY;
            if (IsRestart)
                return;
            if (m_BoundaryValues.IsNullOrEmpty() && !BoundaryIsPeriodic[0] && !BoundaryIsPeriodic[1])
                SetBoundaries(new List<string> { "Wall" });
            GridFunc = delegate {
                DomainVolume = lengthX * lengthY;
                int q = new int(); // #Cells in x-dircetion + 1
                int r = new int(); // #Cells in y-dircetion + 1

                q = (int)(cellsPerUnitLength * lengthX);
                r = (int)(cellsPerUnitLength * lengthY);

                double[] Xnodes = GenericBlas.Linspace(-lengthX / 2, lengthX / 2, q + 1);
                double[] Ynodes = GenericBlas.Linspace(-lengthY / 2, lengthY / 2, r + 1);

                Grid2D grd = Grid2D.Cartesian2DGrid(Xnodes, Ynodes, periodicX: periodicX, periodicY: periodicY);

                for (int i = 0; i < m_BoundaryValues.Count(); i++) {
                    byte iB = (byte)(i + 1);
                    grd.EdgeTagNames.Add(iB, m_BoundaryValues[i]);
                }

                if (m_BoundaryValues.Count() == 0 && periodicX == false && periodicY == false)
                    throw new Exception("Please specify boundaries before creating the grid");

                if (m_BoundaryValues.Count() == 1) {
                    grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 1;
                        if (m_BoundaryValues[0].Contains("wall") || m_BoundaryValues[0].Contains("Wall")) {
                            WallPositionPerDimension[0][0] = -lengthX / 2;
                            WallPositionPerDimension[0][1] = lengthX / 2;
                            WallPositionPerDimension[1][0] = -lengthY / 2;
                            WallPositionPerDimension[1][1] = lengthY / 2;
                        }
                        return et;
                    });
                }
                else {
                    grd.DefineEdgeTags(delegate (double[] X) {
                        byte et = 0;
                        if (Math.Abs(X[0] - (-lengthX / 2)) <= 1.0e-8) {
                            for (int i = 0; i < m_BoundaryValues.Count(); i++) {
                                if (m_BoundaryValues[i].Contains("left") || m_BoundaryValues[i].Contains("Left")) {
                                    et = (byte)(i + 1);
                                    if (m_BoundaryValues[i].Contains("wall") || m_BoundaryValues[i].Contains("Wall")) {
                                        WallPositionPerDimension[0][0] = -lengthX / 2;
                                    }
                                }
                            }
                        }
                        if (Math.Abs(X[0] + (-lengthX / 2)) <= 1.0e-8) {
                            for (int i = 0; i < m_BoundaryValues.Count(); i++) {
                                if (m_BoundaryValues[i].Contains("right") || m_BoundaryValues[i].Contains("Right")) {
                                    et = (byte)(i + 1);
                                    if (m_BoundaryValues[i].Contains("wall") || m_BoundaryValues[i].Contains("Wall")) {
                                        WallPositionPerDimension[0][1] = lengthX / 2;
                                    }
                                }
                            }
                        }

                        if (Math.Abs(X[1] - (-lengthY / 2)) <= 1.0e-8) {
                            for (int i = 0; i < m_BoundaryValues.Count(); i++) {
                                if (m_BoundaryValues[i].Contains("lower") || m_BoundaryValues[i].Contains("Lower")) {
                                    et = (byte)(i + 1);
                                    if (m_BoundaryValues[i].Contains("wall") || m_BoundaryValues[i].Contains("Wall")) {
                                        WallPositionPerDimension[1][0] = -lengthY / 2;
                                    }
                                }
                            }
                        }
                        if (Math.Abs(X[1] + (-lengthY / 2)) <= 1.0e-8) {
                            for (int i = 0; i < m_BoundaryValues.Count(); i++) {
                                if (m_BoundaryValues[i].Contains("upper") || m_BoundaryValues[i].Contains("Upper")) {
                                    et = (byte)(i + 1);
                                    if (m_BoundaryValues[i].Contains("wall") || m_BoundaryValues[i].Contains("Wall")) {
                                        WallPositionPerDimension[1][1] = lengthY / 2;
                                    }
                                }
                            }
                        }
                        Debug.Assert(et != 0);
                        return et;
                    });
                }
                Console.WriteLine("Cells:" + grd.NumberOfCells);
                return grd;
            };
        }


        public void SetTimesteps(double dt, int noOfTimesteps, bool staticTimestep = true) {
            dtMax = dt;
            dtMin = dt;
            Endtime = noOfTimesteps * dt;
            NoOfTimesteps = noOfTimesteps;
            this.staticTimestep = staticTimestep;
        }


        /// <summary>
        /// Set true if the coupling between fluid and particle should be calculated iterative, while using Lie-Splitting.
        /// </summary>
        [DataMember]
        public int fullyCoupledSplittingMaxIterations = 100000;

        /// <summary>
        /// Set true if the complete output should be printed to the console.
        /// </summary>
        [DataMember]
        public bool FullOutputToConsole = false;

        /// <summary>
        /// Function describing the fixed level-set movement
        /// </summary>
        [DataMember]
        public Func<double[], double, double> MovementFunc;
        
        /// <summary>
        /// The termination criterion for fully coupled/implicit level-set evolution.
        /// </summary>
        [DataMember]
        public double hydrodynamicsConvergenceCriterion = 1.0e-6;

        /// <summary>
        /// under-relaxation of the level set movement in case of coupled iterative
        /// </summary>
        [DataMember]
        public double LSunderrelax = 1.0;

        /// <summary>
        /// coefficient of restitution
        /// </summary>
        [DataMember]
        public double CoefficientOfRestitution = 1.0;
        
        /// <summary>
        /// Gravity acting on the particles, zero by default.
        /// </summary>
       // [DataMember]
        public Vector gravity = new Vector(0, 0);

        /// <summary>
        /// See <see cref="LevelSetHandling"/>, Lie-Splitting with iterative coupling by default.
        /// </summary>
        [DataMember]
        public LevelSetHandling Timestepper_LevelSetHandling = LevelSetHandling.FSILieSplittingFullyCoupled;

        /// <summary>
        /// All particles in the FSI
        /// </summary>
        [DataMember]
        public List<Particle> Particles {
            get;
            set;
        }

        /// <summary>
        /// Fix all particles
        /// </summary>
        [DataMember]
        public bool fixPosition = false;

        /// <summary>
        /// For Collisions
        /// </summary>
        [DataMember]
        public double minDistanceThreshold = 0;
        
        /// <summary>
        /// if true the flow solver is turned off
        /// </summary>
        [DataMember]
        public bool pureDryCollisions = false;

        public override Type GetSolverType() {
            return typeof(FSI_SolverMain);
        }
    }
}
