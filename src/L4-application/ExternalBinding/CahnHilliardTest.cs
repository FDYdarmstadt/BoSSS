using System;
using System.Collections;
using System.Collections.Generic;
using BoSSS.Foundation.Grid;

namespace BoSSS.Application.ExternalBinding {

    static public class CahnHilliardTest {

        public static void Main() {

            Init();
            GridImportTest.ConvertFOAMGrid();
            // Console.WriteLine("Running Cahn-Hilliard Test");
            // var chOp = new FixedOperators();
            // OpenFOAMGrid grd = GridImportTest.GenerateFOAMGrid();
            // OpenFoamDGField f = new(grd, 2, 2);
            // OpenFoamMatrix mtx = new(grd, f);
            // // var EdgeValues = new List<List<double>>();
            // // foreach (var val in new double[]{1, -1, 0}){
            // //     EdgeValues.Add(new List<double>{val});
            // // }
            // // OpenFoamPatchField cPtch = new(grd, 1, new int[]{1,2,3}, new string[]{"dirichlet","dirichlet","neumann"}, new double[]{1,-1,0});
            // OpenFoamPatchField cPtch;
            // int[] safeEts = new[] { 1, 2, 3 };
            // // int* ets = (int*)safeEts[0];
            // string[] safeEtyps = new[] { "dirichlet", "dirichlet", "neumann" };
            // // int* eTyps = (int*)safeEtyps[0];
            // double[] safeVals = new[] { 1.0, -1.0, 0 };
            // cPtch = new(grd, 1, safeEts, safeEtyps, safeVals);
            // // unsafe {
            // //     int[] safeEts = new[]{1, 2, 3};
            // //     // int* ets = (int*)safeEts[0];
            // //     int[] safeEtyps = new[]{1, 1, 0};
            // //     // int* eTyps = (int*)safeEtyps[0];
            // //     double[] safeVals = new[]{1.0, -1.0, 0};
            // //     fixed (double* eVals = safeVals){
            // //         fixed (int* ets = safeEts)
            // //         {
            // //             fixed (int* eTyps = safeEtyps)
            // //             {
            // //                 cPtch = new(grd, 1, 3, ets, eTyps, eVals);
            // //             }
            // //         }
            // //     }
            // // }
            // chOp.CahnHilliard(mtx, null, cPtch, null);
            // Cleanup();

        }

        static Initializer MyInit;

        /// <summary>
        /// MPI Init
        /// </summary>
        public static void Init() {
            MyInit = new Initializer();
            MyInit.BoSSSInitialize();
        }

        /// <summary>
        /// MPI shutdown
        /// </summary>
        public static void Cleanup() {
            MyInit.BoSSSFinalize();
        }
    }
}
