using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases {


    /// <summary>
    /// Logging of evaporation, <see cref="HeatedWall"/>
    /// </summary>
    [Serializable]
    public class EvaporationLogging : XNSEinSituPostProcessingModule {
        
        /// <summary>
        /// Probably some specialization for the Fourier level set.
        /// </summary>
        public enum Mode {

            /// <summary>
            /// Line interface ?
            /// </summary>
            LineInterface = 1,

            /// <summary>
            /// circular interface ?
            /// </summary>
            CircleInterface = 2
        }


        /// <summary>
        /// 
        /// </summary>
        [DataMember]
        public Mode mode = Mode.LineInterface;

        /// <summary>
        /// 
        /// </summary>
        public const string LogfileName = "Evaporation";

        /// <summary>
        /// Evaporation
        /// </summary>
        protected override string LogFileName => LogfileName;

        /// <summary>
        /// CSV first line
        /// </summary>
        protected override void WriteHeader(TextWriter textWriter) {
            string header = String.Format("{0}\t{1}\t{2}\t{3}\t{4}", "'#timestep", "time", "interfacePosition", "meanInterfaceVelocity", "meanMassFlux"); //, "Temperatureprofile");
            Log.WriteLine(header);
            Log.Flush();
        }

        /// <summary>
        /// 
        /// </summary>
        protected override void PerformTimestepPostProcessing(int TimestepNo, double phystime) {
            using(new FuncTrace()) {

                MultidimensionalArray InterfacePoints = XNSEUtils.GetInterfacePoints(this.LsTrk, this.LevSet);
                double nNodes = InterfacePoints.GetLength(0);

                double posI = 0.0;
                if(this.mode == Mode.LineInterface)
                    posI = InterfacePoints.ExtractSubArrayShallow(-1, 1).To1DArray().Sum() / nNodes;

                double EvapVelocMean = this.ComputeEvapVelocityMean();
                double hVap = this.Control.ThermalParameters.hVap;
                double MassFlux = EvapVelocMean * hVap;

                string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}", TimestepNo, phystime, posI, EvapVelocMean, MassFlux);

                //// temperature profile
                //int N = 100;
                //double[] tempP = new double[N + 1];

                //double[] probe = new double[N + 1];
                ////if (this.Control.LogValues == XNSE_Control.LoggingValues.EvaporationL) {
                //double L = this.Control.AdditionalParameters[0];
                //double x_probe = this.Control.AdditionalParameters[1];
                //for (int i = 0; i <= N; i++) {
                //    probe = new double[] { x_probe, i * (L / (double)N) };
                //    try {
                //        tempP[i] = this.Temperature.ProbeAt(probe);
                //    } catch {
                //        tempP[i] = 0.0;
                //    }
                //}
                ////}

                //line = line + "\t" + String.Join("\t", tempP.Select(ip => ip.ToString()).ToArray());

                Log.WriteLine(line);
                Log.Flush();


            }
        }

        ConventionalDGField[] ConstructEvaporativeVelocity(DGField[] EvoVelocity) {
            
            if(base.SolverMain is XNSE_SolverMain oldSolver) {
                return oldSolver.ConstructEvaporativeVelocity(oldSolver.GetMeanVelocityFromXDGField(EvoVelocity));
            } else if(base.SolverMain is XNSE newSolver) {
                IList<string> velocityName = BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable("Phi", BoSSS.Solution.NSECommon.VariableNames.VelocityVector(3));
                IReadOnlyDictionary<string, DGField> parameters = newSolver.LsUpdater.Parameters;

                List<ConventionalDGField> velocity = new List<ConventionalDGField>(3);
                for(int i = 0; i < 3; ++i) {
                    
                    if (parameters.TryGetValue(velocityName[i], out DGField velocityField)) {
                        velocity.Add((ConventionalDGField)velocityField);
                    }
                }
                return velocity.ToArray();
                
            } else {
                throw new NotSupportedException();
            }
        }

        double ComputeEvapVelocityMean() {
            double EvapVelocMean = 0.0;

            var evapVelocity = ConstructEvaporativeVelocity(this.CurrentVel);

            int p = evapVelocity[0].Basis.Degree;
            SubGrid sgrd = LsTrk.Regions.GetCutCellSubgrid4LevSet(0);
            NodeSet[] Nodes = LsTrk.GridDat.Grid.RefElements.Select(Kref => Kref.GetQuadratureRule(p * 2).Nodes).ToArray();

            var cp = new ClosestPointFinder(LsTrk, 0, sgrd, Nodes);
            MultidimensionalArray[] VelocityEval = evapVelocity.Select(sf => cp.EvaluateAtCp(sf)).ToArray();
            double nNodes = VelocityEval[0].Length;

            if(this.mode == Mode.LineInterface) {
                double evapVelY = VelocityEval[1].Sum() / nNodes;
                EvapVelocMean = evapVelY;
            } else if(this.mode == Mode.CircleInterface) {
                EvapVelocMean = 0.0;
                for(int s = 0; s < sgrd.GlobalNoOfCells; s++) {
                    for(int n = 0; n < Nodes.Length; n++) {
                        double velX = VelocityEval[0].To2DArray()[s, n];
                        double velY = VelocityEval[1].To2DArray()[s, n];
                        EvapVelocMean += Math.Sqrt(velX.Pow2() + velY.Pow2());
                    }
                }
                EvapVelocMean /= nNodes;
            }

            return EvapVelocMean;
        }
    }
}
