using BoSSS.Solution.AdvancedSolvers.Testing;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ZwoLevelSetSolver.Tests;
using ilPSP.LinSolvers;
using ilPSP.Connectors.Matlab;
using BoSSS.Solution.Statistic;
using NUnit.Framework;
using System.Diagnostics;
using System.IO;
using BoSSS.Foundation.IO;
using ilPSP;
using BoSSS.Foundation.Grid;
using System.Numerics;

namespace ZwoLevelSetSolver {
    class ZLSmain {

        static void Main(string[] args) {
            //int L = Vector<double>.Count;



            BoSSS.Solution.Application.InitMPI(num_threads:1);
            //BoSSS.Solution.Application.DeleteOldPlotFiles();

            //RunSolver(args);
            ReferenceRun();
            //ComparisonRun();
            
            BoSSS.Solution.Application.FinalizeMPI();
        }

        static void RunSolver(string[] args) {
            //Debugger.Launch();
            ZLS._Main(args, false, delegate () {
                //Control file from runtime via args
                var p = new ZLS();
                return p;
            });
        }

        static void ReferenceRun() {
            var dbPath = "bosss_db_debug";
            if(!Directory.Exists(dbPath))
                BoSSS.Foundation.IO.DatabaseUtils.CreateDatabase(dbPath);
            var db = DatabaseInfo.Open(dbPath);
            
            var c = ZwoLevelSetSolver.ControlFiles.Droplet.AlandSL3D(2, 2, 0);
            c.DbPath = db.Path;
            c.savetodb = true;
            c.NoOfTimesteps = 1;

            using(var p = new ZLS()) {
                p.Init(c);
                p.RunSolverMode();

                foreach(var field in p.CurrentState.Fields) {
                    double totalJumpNorm = field.JumpNorm();
                    double inrprJumpNorm = field.JumpNorm(p.GridData.GetInterprocessEdges());

                    Console.WriteLine($"Jump norm of {field.Identification}: \t{totalJumpNorm:g6} \t(interprocess: {inrprJumpNorm:g7})");
                }
            }

           
        }

       

    }

}

