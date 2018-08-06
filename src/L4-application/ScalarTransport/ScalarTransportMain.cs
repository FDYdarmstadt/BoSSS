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
using BoSSS.Solution;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.Tecplot;
using BoSSS.Foundation.Comm;
using BoSSS.Foundation.IO;
using System.Diagnostics;
using BoSSS.Foundation.Quadrature;
using MPI.Wrappers;
using BoSSS.Platform;
using ilPSP.Tracing;
using BoSSS.Solution.Utils;
using System.IO;
using BoSSS.Foundation.Grid.Classic;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.RefElements;
using ilPSP;
using BoSSS.Foundation.Caching;
using BoSSS.Foundation.Grid.Aggregation;

namespace BoSSS.Application.ScalarTransport {

    /// <summary>
    /// Application class of my own BoSSS Application
    /// </summary>
    class ScalarTransportMain : BoSSS.Solution.Application {

        /// <summary>
        /// Entry point of my program
        /// </summary>
        /// <param name="args"></param>
        static void Main(string[] args) {
            BoSSS.Solution.Application._Main(args, true, delegate() {
                return new ScalarTransportMain();
            });
        }

        /// <summary>
        /// creates a simple 2d/3d Cartesian grid within the domain \f$ (-7,7)^D \f$
        /// </summary>
        protected override GridCommons CreateOrLoadGrid() {
            /*
            double[] xnodes = GenericBlas.Linspace(-7, 7, 11);
            double[] ynodes = GenericBlas.Linspace(-7, 7, 11);
            GridCommons grd = Grid2D.Cartesian2DGrid(xnodes, ynodes, type: CellType.Square_Linear);
            this.Control.NoOfMultigridLevels = 3;

            return grd;
            */
            //double[] xnodes = GenericBlas.Linspace(-7, 7, 25);
            //double[] ynodes = GenericBlas.Linspace(-7, 7, 25);
            //double[] znodes = GenericBlas.Linspace(-7, 7, 25);
            //var grd = Grid3D.Cartesian3DGrid(xnodes, ynodes, znodes);

            
            double[] xNodes = GenericBlas.Linspace(-7, 7, 4);
            double[] yNodes = GenericBlas.Linspace(-7, 7, 5);

            var baseGrid = Grid2D.UnstructuredTriangleGrid(xNodes, yNodes);
            var baseGdat = new GridData(baseGrid);
            var aggGrid = CoarseningAlgorithms.Coarsen(baseGdat, 2);
            base.AggGrid = aggGrid;
                  

            return null;
            
        }


        /// <summary>
        /// the scalar property that is convected
        /// </summary>
        DGField u;

        /// <summary>
        /// transport velocity
        /// </summary>
        VectorField<SinglePhaseField> Velocity;

        /// <summary>
        /// stores the mpi processor rank - just for illustration
        /// </summary>
        DGField mpi_rank;

        /// <summary>
        /// creates the field <see cref="u"/>;
        /// </summary>
        protected override void CreateFields() {
            u = new SinglePhaseField(new Basis(this.GridData, 5), "u");
            mpi_rank = new SinglePhaseField(new Basis(this.GridData, 0), "MPI_rank");
            Velocity = new VectorField<SinglePhaseField>(this.GridData.SpatialDimension, new Basis(this.GridData, 2), "Velocity", SinglePhaseField.Factory);

            m_IOFields.Add(u);

            m_RegisteredFields.Add(u);
            m_RegisteredFields.AddRange(Velocity);
            m_RegisteredFields.AddRange(mpi_rank);
        }

        ///// <summary>
        ///// time-stepping algorithm for this method;
        ///// </summary>
        //ROCK4 Timestepper;

        ExplicitEuler Timestepper;

        SpatialOperator diffOp;
        
        /// <summary>
        /// 
        /// </summary>
        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L) {
            diffOp = new SpatialOperator(new string[] { "u" },
                Solution.NSECommon.VariableNames.VelocityVector(this.GridData.SpatialDimension),
                new string[] { "codom1" },
                QuadOrderFunc.Linear());


            switch (this.GridData.SpatialDimension) {
                case 2:
                    diffOp.EquationComponents["codom1"].Add(new ScalarTransportFlux());
                    break;
                case 3:
                    diffOp.EquationComponents["codom1"].Add(new ScalarTransportFlux3D());
                    break;
                default:
                    throw new NotImplementedException("spatial dim. not supported.");
            }
            diffOp.Commit();

            Timestepper = new RungeKutta(RungeKuttaScheme.TVD3,
                                   diffOp,
                                   u.Mapping, Velocity.Mapping);

            Solution.Multigrid.AggregationGridBasis.CreateSequence(this.MultigridSequence, u.Mapping.BasisS);


            //Timestepper = new ROCK4(diffOp, u.CoordinateVector, Velocity.Mapping);
        }

        /// <summary>
        /// Auto-tuning
        /// </summary>
        public void PerformanceVsCachesize() {
            double[] dummy = new double[this.u.CoordinateVector.Count];

            var eval = diffOp.GetEvaluatorEx(new DGField[] { this.u }, this.Velocity.ToArray(), this.u.Mapping,
                edgeQrCtx:new EdgeQuadratureScheme(false, EdgeMask.GetEmptyMask(this.GridData)),
                volQrCtx:new CellQuadratureScheme(true,null));

            Stopwatch stw = new Stopwatch();
            int NoOfRuns = 1;
            Console.WriteLine("BlkSz\tNoChks\tRunTime");
            //var testsizes = new int[] { 1, 2, 4, 8, 16, 32, 64, 256, 512, 1024, 2048, 8192, 16384, 16384 * 2, 16384 * 4, 100000000 };
            var testsizes = new int[] { 1, 2, 4, 1024, 16384 * 2, 16384 * 4, 100000000 };
            foreach(int bulkSize in testsizes) {
                Quadrature_Bulksize.CHUNK_DATA_LIMIT = bulkSize;

                stw.Reset();
                stw.Start();
                for(int i = 0; i < NoOfRuns; i++) {
                    
                    eval.Evaluate(1.0,0.0, dummy);
                   
                }
                 stw.Stop();
                 
                double runTime = stw.Elapsed.TotalSeconds;
                Console.WriteLine("{0}\t{2}\t{1}", bulkSize, runTime.ToStringDot("0.####E-00"),this.GridData.iGeomCells.Count/bulkSize);
            }
        }


        public void SimplifiedPerformance() {
            int J = this.GridData.iGeomCells.Count;
            var NodeSet = this.GridData.iGeomCells.RefElements[0].GetQuadratureRule(12).Nodes;

            int[] ChunkSize = new[] { J / 2, J };
            if (J % 2 != 0)
                throw new Exception();

            var stw = new Stopwatch();


            foreach(int sz in ChunkSize) {
                Console.WriteLine("No Of chunks: " + (J/sz));

                for (int sweep = 0; sweep < 2; sweep++) {
                    double dummyRes = 0;
                    Console.WriteLine("  Sweep # {0}...", sweep);
                   
                //Cache.ClearCache();

                    stw.Reset();
                    stw.Start();
                    for (int j = 0; j < J; j += sz) {
                        //var invJac = this.GridData.InverseJacobian.GetValue_Cell(NodeSet, j, sz);
                        var detJac = this.GridData.JacobianDeterminat.GetValue_Cell(NodeSet, j, sz);

                        //dummyRes += invJac.Storage.L2NormPow2();
                        dummyRes += detJac.Storage.L2NormPow2();
                    }
                    stw.Stop();

                    Console.WriteLine("   done. Runtime {0} msec.", stw.ElapsedMilliseconds);
                    Console.WriteLine("   Cache hits/misses: {0}/{1}", Cache.Hits, Cache.Misses);

                    Console.WriteLine("                dummy result (" + dummyRes + ")");
                }

            }


        }


        protected override int[] ComputeNewCellDistribution(int TimeStepNo, double physTime) {
            if(this.MPISize == 2) {
                int J = this.GridData.iLogicalCells.NoOfLocalUpdatedCells;
                int[] part = new int[J];

                int R0Cells = 0;
                int R1Cells = 0;

                for(int j = 0; j < J; j++) {
                    if (this.u.GetMeanValue(j) > 0.3) {
                        part[j] = 1;
                        R1Cells++;
                    } else {
                        part[j] = 0;
                        R0Cells++;
                    }
                }

                R0Cells = R0Cells.MPISum();
                R1Cells = R1Cells.MPISum();
                Console.WriteLine("Redist: rank 0 {0} cells, rank 1 {1} cells.", R0Cells, R1Cells);
                //Debugger.Launch();

                return part;

            } else {
                return null;
            }
        }


        /// <summary>
        /// performs some timesteps
        /// </summary>
        protected override double RunSolverOneStep(int TimestepNo, double phystime, double dt) {
            using (var ft = new FuncTrace()) {
                //PerformanceVsCachesize();
                //SimplifiedPerformance();

                this.u.ProjectField((x, y) => x * y);
                
                base.TerminationKey = true;
                return 0.0;


                //u.ProjectField((_2D)((x, y) => x+y));
                //var f = new SinglePhaseField(u.Basis, "f");
                //var sgrd = new SubGrid(new CellMask(this.GridDat, Chunk.GetSingleElementChunk(200)));
                //f.Clear();
                //f.DerivativeByFlux(1.0, u, 0, sgrd, true);
                //Tecplot plt1 = new Tecplot(GridDat, true, false, 0);
                //plt1.PlotFields("ggg" , "Scalar Transport", phystime, u, f, mpi_rank);
                //base.TerminationKey = true;


                double dtCFL;
                if(this.GridData is GridData) {
                    dtCFL = ((GridData)(this.GridData)).ComputeCFLTime(this.Velocity, 1.0e10);
                } else {
                    Console.WriteLine("Nix CFL");
                    dtCFL = 1e-3;
                }
                if (dt <= 0) {
                    // if dt <= 0, we're free to set the timestep on our own
                    //base.NoOfTimesteps = -1;
                    //dt = 1;

                    base.NoOfTimesteps = 10;
                    base.EndTime = 5.0;
                    dtCFL *= 1.0 / (((double)u.Basis.Degree).Pow2());

                }


                Console.Write("Timestp. #" + TimestepNo + " of " + base.NoOfTimesteps + " ... \t");
                dt = dtCFL * 0.5;
                Timestepper.Perform(dt);

                // set mpi_rank
                double rank = GridData.MpiRank;
                int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
                for (int j = 0; j < J; j++) {
                    mpi_rank.SetMeanValue(j, rank);
                }


                //dt = Timestepper.PerformAdaptive(1.0e2);

                Console.WriteLine("finished (dt = {0:0.###E-00}, CFL frac = {1:0.###E-00})!", dt, dt / dtCFL);


                return dt;
            }
        }

        /// <summary>
        /// performs Tecplot output of field <see cref="u"/>
        /// </summary>
        protected override void PlotCurrentState(double phystime, TimestepNumber timestepNo, int superSampling) {
            Tecplot plt1 = new Tecplot((GridData)GridData, true, false, (uint) superSampling);
            plt1.PlotFields("transport." + timestepNo, phystime, u, mpi_rank);
        }

        /// <summary>
        /// sets some initial value for field <see cref="u"/>;
        /// </summary>
        protected override void SetInitial() {

            switch (GridData.SpatialDimension) {
                case 2:

                
                
                u.ProjectField((_2D)(delegate(double x, double y) {
                    
                    double r = Math.Sqrt(x*x + y*y);


                    return Math.Exp(-r * r);
                }));
                               

                Velocity[0].ProjectField((_2D)((x, y) => 1.0));
                Velocity[1].ProjectField((_2D)((x, y) => 0.1));

                break;

                case 3:
                u.ProjectField((x, y, z) => ((Math.Abs(x) <= 3.0 && Math.Abs(y) <= 3.0 && Math.Abs(z) <= 3.0) ? 1 : 0));
                Velocity[0].ProjectField((_3D)((x, y, z) => 1));
                Velocity[1].ProjectField((_3D)((x, y, z) => 0));
                Velocity[2].ProjectField((_3D)((x, y, z) => 0));

                
                break;

                default:
                throw new NotImplementedException();
            }




            double rank = GridData.MpiRank;
            int J = GridData.iLogicalCells.NoOfLocalUpdatedCells;
            for (int j = 0; j < J; j++) {
                mpi_rank.SetMeanValue(j, rank);
            }
        }
    }
}

