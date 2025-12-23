using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using MPI.Wrappers;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Security.Cryptography.X509Certificates;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Solution.Statistic;
using ilPSP;


namespace ZwoLevelSetSolver.Tests {

    [TestFixture]
    public static class ThreePhaseTests {


        [Test]
        public static void ThreePhaseContactLine_TensionBalance([Values(true)] bool UseGravity) {


            var C = ZwoLevelSetSolver.ControlFiles.Droplet.Wiki(p: 4, AMRlvl: 2, SlipLength: 0, useGravity: UseGravity);
            C.NoOfTimesteps = 1;
            
            //if(UseGravity == false) {
            //    // turn off gravity
            //    var toRemove = C.InitialValues.Keys.Where(name => name.StartsWith("Gravity")).ToArray();
            //    foreach(var name in toRemove) 
            //        C.InitialValues.Remove(name);
            //}
           

            using(var S = new ZLS()) {
                S.Init(C);
                S.RunSolverMode();

                double[] L2_vel = new double[2];
                for(int d = 0; d < 2; d++) {
                    L2_vel[d] = S.Velocity[d].L2Norm();
                    Console.WriteLine($"Velocity {d} norm: " + L2_vel[d]);
                }
                // with gravity
                // Velocity 0 norm: 9.070818979018125E-06
                // Velocity 1 norm: 1.5000248940334127E-05

                if(UseGravity) {
                    Assert.Less(L2_vel[0], 5.0e-5, "VelocityX out of bound");
                    Assert.Less(L2_vel[1], 5.0e-5, "VelocityY out of bound");
                }
            }


        }


        /// <summary>
        /// Test whether parallel and serial runs produce roughly the same results for the Aland-testcase.
        ///
        /// Especially before the three-phase code was merged with the master branch again, large jumps seem to occur at MPI boundaries.
        /// </summary>
        [Test]
        public static void AlandTest() {
            // --test=ZwoLevelSetSolver.Tests.ThreePhaseTests.AlandTest
            bool ReferenceRun = csMPI.Size_World == 1;

            var dbPath = "bosss_db_debug";
            if(ReferenceRun) {
                if(Directory.Exists(dbPath))
                    throw new IOException("database seems to exist already; expecting a plain field for the reference run.");
                BoSSS.Foundation.IO.DatabaseUtils.CreateDatabase(dbPath);
            }

            var db = DatabaseInfo.Open(dbPath);

            var c = ZwoLevelSetSolver.ControlFiles.Droplet.AlandSL3D(2, 2, 0);
            c.DbPath = db.Path;
            c.savetodb = true;
            c.NoOfTimesteps = 3;

            using(var p = new ZLS()) {
                p.Init(c);
                p.RunSolverMode();

                foreach(var field in p.CurrentState.Fields) {
                    double totalJumpNorm = field.JumpNorm();
                    double inrprJumpNorm = field.JumpNorm(p.GridData.GetInterprocessEdges());

                    Console.WriteLine($"Jump norm of {field.Identification}: \t{totalJumpNorm:g6} \t(interprocess: {inrprJumpNorm:g7})");
                    Assert.Less(inrprJumpNorm, totalJumpNorm*0.1, "inter-process jump seems suspiciously high");
                }
            }

            if(!ReferenceRun) {
                // +++++++++++++++++++++++++++++++++++++++++
                // compare 8-Ranks run to Singe Rank run
                // +++++++++++++++++++++++++++++++++++++++++
                var si1 = db.Sessions.Single(si => si.ComputeNodeNames.Count() == 1);
                var si8 = db.Sessions.Single(si => si.ComputeNodeNames.Count() != 1);

                string[] FieldsToCompare = new string[] { "VelocityX", "VelocityY", "VelocityZ", "Pressure", "DisplacementX", "DisplacementY", "DisplacementZ" };

                int cnt = 1;
                foreach(var ts in si8.Timesteps) {
                    var ref_ts = si1.Timesteps.Single(_ts => _ts.TimeStepNumber.Equals(ts.TimeStepNumber));
                    Console.WriteLine();
                    Console.WriteLine();
                    Console.WriteLine($"Comparing {ts.TimeStepNumber}, t = {ts.PhysicalTime} (cnt = {cnt}): ");
                    Console.WriteLine("     " + ts.ID + $"\t({ts.Grid.NumberOfCells} cells)");
                    Console.WriteLine(" vs. " + ref_ts.ID + $"\t({ts.Grid.NumberOfCells} cells)");
                    Console.WriteLine();

                    var fields2Plot = new List<DGField>();
                    foreach(var id in FieldsToCompare) {
                        var f8 = ts.Fields.Single(f => f.Identification == id) as XDGField;
                        var f1 = ref_ts.Fields.Single(f => f.Identification == id) as XDGField;

                        var trk = f8.Basis.Tracker;
                        var f1_inj = f8.CloneAs();
                        f1_inj.Identification = f1_inj.Identification + "Ref";
                        f1_inj.Clear();
                        f1_inj.ProjectFromForeignGrid(1.0, f1);

                        
                        double L2f8 = f8.L2NormAllSpecies();
                        double L2f1_inj = f1_inj.L2NormAllSpecies();

                        var f1_err = f1_inj.CloneAs();
                        f1_err.Identification = f1_err.Identification + "ERR";
                        f1_err.Acc(-1.0, f8);

                        double L2err = f1_err.L2NormAllSpecies();


                        Assert.IsFalse(L2f8.IsNaNorInf(), "NAN/INF in f8, " + id);
                        Assert.IsFalse(L2f1_inj.IsNaNorInf(), "NAN/INF in f1_inj, " + id);
                        Assert.IsFalse(L2err.IsNaNorInf(), "NAN/INF in err, " + id);


                        Console.WriteLine($"L2 Error for {f8.Identification}: \t{L2err:g6} \t**{L2err / Math.Max(L2f1_inj, L2f8):g6}** \t ({L2f1_inj:g5} \t{L2f8:g5})");


                        NUnit.Framework.Assert.LessOrEqual(L2err, Math.Max(L2f1_inj, L2f8) * 0.03, $"Difference for {id} to large");


                        fields2Plot.Add(f8);
                        fields2Plot.Add(f1_inj);
                        fields2Plot.Add(f1_err);
                    }

                    fields2Plot.Add(ts.Fields.Single(f => f.Identification == "Phi"));
                    fields2Plot.Add(ts.Fields.Single(f => f.Identification == "Phi2"));
                    fields2Plot.Add(ts.Fields.Single(f => f.Identification == "MPIrank"));
                    //BoSSS.Solution.Tecplot.Tecplot.PlotFields(fields2Plot, "c:\\tmp\\err" + cnt, 0.0, 4);
                    
                    
                    
                    
                    
                    cnt++;
                }


            }


        }
    }
}
