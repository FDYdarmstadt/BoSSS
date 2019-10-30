using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Connectors.Matlab;
using BoSSS.Foundation.Grid;
using MPI.Wrappers;
using ilPSP.Utils;
using System.Diagnostics;

namespace BoSSS.Solution {

    public class OpAnalysisBase {

        /// <summary>
        /// Delegate for computation of operator matrix with:
        /// </summary>
        /// <param name="OpMtx">
        /// operator matrix </param>
        /// <param name="OpAffine">
        /// affine vector</param>
        /// <param name="Mapping">
        /// mapping</param>
        /// <param name="CurrentState">
        /// current state</param>
        /// <param name="AgglomeratedCellLengthScales"></param>
        /// <param name="time"></param>
        public delegate void DelComputeOperatorMatrix(BlockMsrMatrix OpMtx, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double time);


        /// <summary>
        /// Ctor for operator analysis class
        /// </summary>
        /// <param name="delComputeOperatorMatrix"></param>
        /// <param name="Mapping"></param>
        /// <param name="CurrentState"></param>
        /// <param name="AgglomeratedCellLengthScales"></param>
        /// <param name="time"></param>
        public OpAnalysisBase(DelComputeOperatorMatrix delComputeOperatorMatrix, UnsetteledCoordinateMapping Mapping, DGField[] CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double time){

            //System.Threading.Thread.Sleep(10000);

            m_OpMtx = new BlockMsrMatrix(Mapping, Mapping); //operator matrix
            localRHS = new double[Mapping.LocalLength]; //right hand side
            RHSlen = Mapping.TotalLength;
            m_map = Mapping; // mapping

            VarGroup = m_map.BasisS.Count.ForLoop(i => i); //default: all dependent variables are included in operator matrix
            delComputeOperatorMatrix(m_OpMtx, localRHS, Mapping, CurrentState, AgglomeratedCellLengthScales, time); // delegate for computing the operator matrix
        }


        public OpAnalysisBase(BlockMsrMatrix Mtx, double[] RHS, UnsetteledCoordinateMapping Mapping) {

            //System.Threading.Thread.Sleep(10000);

            m_OpMtx = new BlockMsrMatrix(Mapping, Mapping); //operator matrix
            localRHS = new double[Mapping.LocalLength]; //right hand side
            RHSlen = Mapping.TotalLength;
            m_map = Mapping; // mapping

            VarGroup = m_map.BasisS.Count.ForLoop(i => i); //default: all dependent variables are included in operator matrix
            m_OpMtx = Mtx;
            localRHS = RHS;


        }


        BlockMsrMatrix m_OpMtx;
        int RHSlen;
        UnsetteledCoordinateMapping m_map;
        int[] _VarGroup;
        double[] localRHS;

        /// <summary>
        /// user-defined indices of dependend variables, if not the full matrix should be analyzed, e.g. 0 = u_x, 1=u_y, 2=u_z, 3=p ...
        /// </summary>
        public int[] VarGroup{
            get{
                return _VarGroup;
            }
            set{
                _VarGroup = value;
                if (_VarGroup.Distinct().Count() != _VarGroup.Count()){
                    throw new ArithmeticException("duplicate found!");
                }

                foreach (int k in _VarGroup){
                    if (k >= m_map.BasisS.Count){
                        throw new ArithmeticException("dependent variables out of range!!");
                    }
                    else if (k < 0){
                        throw new ArithmeticException("dependent variables out of range!!");
                    }
                }
            }
        }


        /// <summary>
        /// fully analyses the matrix with
        /// - condition number,
        /// - symmetric positive definiteness,
        /// - maximum and minimum eigenvalues and
        /// - the existence of a unique solution
        /// </summary>
        public void Analyse() {

            Console.WriteLine("Starting matrix Analysis");
            Console.WriteLine("Calculating condition number");
            double[] CondNum_Write = CondNum();
            Console.WriteLine("Doing symmetry test");
            bool[] Symmetry_Write = Symmetry();
            Console.WriteLine("Doing eigenvalues test");
            double[] Eigenval_Write = Eigenval();

            Console.WriteLine("");
            Console.WriteLine("==================================================================");
            Console.WriteLine("Log of the analysis");
            Console.WriteLine("==================================================================");
            Console.WriteLine("Condition number:");
            Console.WriteLine("full matrix: {0}", CondNum_Write[0]);
            Console.WriteLine("inner matrix: {0}", CondNum_Write[1]);
            Console.WriteLine("==================================================================");
            Console.WriteLine("Symmetry and positive definiteness:");
            Console.WriteLine("is symmetric: {0}", Symmetry_Write[0]);
            Console.WriteLine("is positive definite: {0}", Symmetry_Write[1]);
            Console.WriteLine("==================================================================");
            Console.WriteLine("Eigenvalues:");
            Console.WriteLine("maximal eigenvalue: {0}", Eigenval_Write[0]);
            Console.WriteLine("minimal eigenvalue: {0}", Eigenval_Write[1]);

            rankAnalysis(m_OpMtx, localRHS);

            Console.WriteLine("");

        }


        /// <summary>
        /// According to the Rouché-Capelli theorem, the system is inconsistent if rank(augMatrix) > rank(Matrix). 
        /// If rank(augMatrix) == rank(Matrix), the system has at least one solution.
        /// Additionally, if the rank is equal to the number of variables, the solution of the system is unique.
        /// Remark: This requires the whole RHS not only local!!!
        /// </summary>
        /// <param name="OpMatrix"></param>
        /// <param name="RHS"></param>
        public void rankAnalysis(BlockMsrMatrix OpMatrix, double[] RHS){

            MultidimensionalArray outputArray = MultidimensionalArray.Create(2, 1); // The two rank values

            //At this point OpMatrix and RHS are local, they are collected within bmc on proc rank==0
            using (var bmc = new BatchmodeConnector()){
                bmc.PutSparseMatrix(OpMatrix, "OpMatrix");
                bmc.PutVector(RHS, "RHS");
                bmc.Cmd("output = zeros(2,1)"); // First value is rank(OpMatrix), second the rank of the augmented matrix = rank([Matrix|RHS])
                bmc.Cmd("");
                bmc.Cmd("fullMtx = full(OpMatrix);");
                bmc.Cmd("augmentedMtx = [fullMtx RHS];");
                bmc.Cmd("output(1) = rank(fullMtx)");
                bmc.Cmd("output(2) = rank(augmentedMtx)");
                bmc.GetMatrix(outputArray, "output");
                bmc.Execute(false);
            }


            double[] output = new double[2];
            output[0] = outputArray[0, 0]; //Rank matrix 
            output[1] = outputArray[1, 0]; //Rank augmented matrix ( [Matrix|RHS])

            double rnkMtx = output[0];
            double rnkAugmentedMtx = output[1];

            // Some tests
            Console.WriteLine("==================================================================");

            //Output
            {
                Console.WriteLine("Results of rank analysis:");
                if (rnkAugmentedMtx > rnkMtx)
                {
                    //throw new Exception("The rank of the augmented matrix shouldn't be greater than the one of the original matrix!!"); 
                    Console.WriteLine("======================================================");
                    Console.WriteLine("WARNING!!!!!!! The rank of the augmented matrix shouldn't be greater than the one of the original matrix!!");
                    Console.WriteLine("This means that the system doesnt have a solution!");
                }

                if (rnkAugmentedMtx == rnkMtx)
                {
                    Console.WriteLine("The system has at least a solution");
                }

                //RHS and OpMatrix will be collected in Bmc, so total length has to be considered for RHS: RHS.length will lead to errors in parallel execution
                if (rnkMtx < RHSlen)
                {
                    Console.WriteLine("The rank of the matrix is smaller than the number of variables. There are {0} free parameters", (RHS.Length - rnkMtx));
                }

                else if (rnkMtx == RHSlen)
                {
                    Console.WriteLine("The system has a unique solution :) ");
                }

                else
                {
                    throw new Exception("what? should not happen");
                }

                Console.WriteLine("Rank of the matrix : {0} \n" + "Rank of the augmented matrix : {1} \n" + "Number of variables: {2}", output[0], output[1], RHS.Length);
                Console.WriteLine("==================================================================");
            }

            Debug.Assert(output[0].MPIEquals(), "value does not match on procs");
            Debug.Assert(output[1].MPIEquals(), "value does not match on procs");
        }


        /// <summary>
        /// returns the condition number of the full matrix and the inner matrix without boundary terms
        /// </summary>
        /// <returns>Array condestOut=[ConditionNumberFullOp, ConditionNumberInnerOp]</returns>
        public double[] CondNum(){

            int[] DepVars = this.VarGroup;
            var grd = m_map.GridDat;
            int NoOfCells = grd.Grid.NumberOfCells;
            int NoOfBdryCells = grd.GetBoundaryCells().NoOfItemsLocally_WithExternal;

            // only for full matrix, if there are no inner cells
            if (NoOfCells == NoOfBdryCells){

                Console.WriteLine("");
                Console.WriteLine("Since there are only boundary cells, the condest of the non-existing inner matrix is set to zero!");

                double[] Full_0Vars = this.m_map.GetSubvectorIndices(true, DepVars).Select(i => i + 1.0).ToArray();
                MultidimensionalArray output = MultidimensionalArray.Create(1, 1);

                using (BatchmodeConnector bmc = new BatchmodeConnector()){

                    // if Octave should be used instead of Matlab....
                    //BatchmodeConnector.Flav = BatchmodeConnector.Flavor.Octave;

                    bmc.PutSparseMatrix(m_OpMtx, "FullMatrix");
                    bmc.PutVector(Full_0Vars, "Full_0Vars");
                    bmc.Cmd("output = zeros(1,1)");
                    bmc.Cmd("output = condest(FullMatrix(Full_0Vars,Full_0Vars));");
                    bmc.GetMatrix(output, "output");
                    bmc.Execute(false);
                }

                double condestFull = output[0, 0];
                double condestInner = 0;

                double[] condestOut = new double[] { condestFull, condestInner };
                return condestOut;
            }
            // for full and inner matrix
            else{

                CellMask InnerCellsMask = grd.GetBoundaryCells().Complement();
                SubGrid InnerCells = new SubGrid(InnerCellsMask);

                double[] Full_0Vars = this.m_map.GetSubvectorIndices(true, DepVars).Select(i => i + 1.0).ToArray();
                double[] Inner_0Vars = this.m_map.GetSubvectorIndices(InnerCells, true, DepVars).Select(i => i + 1.0).ToArray();

                MultidimensionalArray output = MultidimensionalArray.Create(2, 1);
                string[] names = new string[] { "Full_0Vars", "Inner_0Vars" };

                using (BatchmodeConnector bmc = new BatchmodeConnector()){

                    // if Octave should be used instead of Matlab....
                    // BatchmodeConnector.Flav = BatchmodeConnector.Flavor.Octave;

                    bmc.PutSparseMatrix(m_OpMtx, "FullMatrix");
                    bmc.PutVector(Inner_0Vars, "Inner_0Vars");
                    bmc.PutVector(Full_0Vars, "Full_0Vars");
                    bmc.Cmd("output = zeros(2,1)");

                    int k = 1;
                    foreach (var s in names){
                        bmc.Cmd("output({1}) = condest(FullMatrix({0},{0}));", s,k);
                        k++;
                    }

                    bmc.GetMatrix(output, "output");
                    bmc.Execute(false);

                    double condestFull = output[0, 0];
                    double condestInner = output[1, 0];

                    double[] condestOut = new double[] { condestFull, condestInner };


                    Debug.Assert(condestOut[0].MPIEquals(),"value does not match on procs");
                    Debug.Assert(condestOut[1].MPIEquals(), "value does not match on procs");
                    return condestOut;
                }
            }
        }


        /// <summary>
        /// Test if the matrix is symmetric positive definite
        /// </summary>
        /// <returns>bool array res=[symmetry, positive definit]</returns>
        public bool[] Symmetry(){

            bool[] res = new bool[2];

            //extract submatrix for selected dependent variables
            int[] DepVars = this.VarGroup;
            double[] DepVars_subvec = this.m_map.GetSubvectorIndices(true, DepVars).Select(i => i + 1.0).ToArray();

            //MsrMatrix OpMtxMSR = m_OpMtx.ToMsrMatrix();

            int[] SubMatrixIdx_Row = m_map.GetSubvectorIndices(false, DepVars);
            int[] SubMatrixIdx_Cols = m_map.GetSubvectorIndices(false, DepVars);
            int L = SubMatrixIdx_Row.Length;

            MsrMatrix SubOpMtx = new MsrMatrix(L, L, 1, 1);

            m_OpMtx.WriteSubMatrixTo(SubOpMtx, SubMatrixIdx_Row, default(int[]), SubMatrixIdx_Cols, default(int[]));


            // symmetry test by calculation of symmetry deviation
            bool sym = false;
            double SymmDev = m_OpMtx.SymmetryDeviation();
            if (SymmDev < 1e-5)
                sym = true;

            res[0] = sym;


            // positive definite test by Cholesky decomposition
            var FullyPopulatedMatrix = m_OpMtx.ToFullMatrixOnProc0();

            bool posDef = true;
            // only proc 0 gets info so the following is executed exclusively on rank 0
            if (ilPSP.Environment.MPIEnv.MPI_Rank == 0)
            {
                try
                {
                    FullyPopulatedMatrix.Cholesky();
                }
                catch (ArithmeticException)
                {
                    posDef = false;
                }
            }
            res[1] = MPIEnviroment.Broadcast<bool>(posDef, 0, ilPSP.Environment.MPIEnv.Mpi_comm);

            Debug.Assert(res[0].MPIEquals(), "value does not match on procs");
            Debug.Assert(res[1].MPIEquals(), "value does not match on procs");
            return res;
        }


        /// <summary>
        /// returns the maximum and minimum eigenvalues of the matrix
        /// </summary>
        /// <returns>Array myeigs =[MaximumEig, MinimumEig]</returns>
        public double[] Eigenval(){

            double[] eigenvalues = new double[2];
            MultidimensionalArray eigs = MultidimensionalArray.Create(1, 2);
            MultidimensionalArray output = MultidimensionalArray.Create(2, 1);

            int[] DepVars = this.VarGroup;
            double[] DepVars_subvec = this.m_map.GetSubvectorIndices(true, DepVars).Select(i => i + 1.0).ToArray();

            using (BatchmodeConnector bmc = new BatchmodeConnector()){

                // if Octave should be used instead of Matlab....
                //BatchmodeConnector.Flav = BatchmodeConnector.Flavor.Octave;

                bmc.PutSparseMatrix(m_OpMtx, "FullMatrix");
                bmc.PutVector(DepVars_subvec, "DepVars_subvec");
                bmc.Cmd("output = zeros(2,1)");
                bmc.Cmd("output(1) = eigs(FullMatrix(DepVars_subvec,DepVars_subvec),1,'lm');");
                bmc.Cmd("output(2) = eigs(FullMatrix(DepVars_subvec,DepVars_subvec),1,'sm');");
                bmc.GetMatrix(output, "output");
                bmc.Execute(false);     
            }

            double[] myeigs = new double[] { output[0, 0], output[1, 0] };
            Debug.Assert(output[0, 0].MPIEquals(), "value does not match on procs");
            Debug.Assert(output[1, 0].MPIEquals(), "value does not match on procs");
            return myeigs;
        }

    }
}
