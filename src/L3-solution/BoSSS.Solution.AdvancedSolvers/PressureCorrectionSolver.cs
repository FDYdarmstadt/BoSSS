using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.LinSolvers.PARDISO;
using ilPSP.Tracing;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Text;

namespace BoSSS.Solution.AdvancedSolvers {
    class PressureCorrectionSolver : ISolverSmootherTemplate {
        public int IterationsInNested => throw new NotImplementedException();

        public int ThisLevelIterations => throw new NotImplementedException();

        public bool Converged => throw new NotImplementedException();

        public object Clone() {
            throw new NotImplementedException();
        }

        public void Dispose() {
            if(PressureSolver != null) {
                PressureSolver.Dispose();
                PressureSolver = null;
            }

            velCompsMask = null;
            prsCompsMask = null;
            VelDiv = null;
            RHSmatrix = null;
        }

        public void Init(MultigridOperator op) {
            using (new FuncTrace()) {
                ICoordinateMapping map = op.Mapping;
                var Mtx = op.OperatorMatrix;

                int D = op.Mapping.SpatialDimension;

                SubBlockSelector velComps = new SubBlockSelector(map);
                velComps.SetVariableSelector(D.ForLoop(iVar => iVar));
                velCompsMask = new BlockMask(velComps);

                SubBlockSelector prsComps = new SubBlockSelector(map);
                prsComps.SetVariableSelector(new int[] { D });
                prsCompsMask = new BlockMask(prsComps);


                // convection-diffusion part
                var CD = velCompsMask.GetSubBlockMatrix(Mtx, Mtx.MPI_Comm);

                // pressure gradient
                var PrGrad = velCompsMask.GetSubBlockMatrix(Mtx, prsCompsMask, Mtx.MPI_Comm);

                // pressure gradient
                VelDiv = prsCompsMask.GetSubBlockMatrix(Mtx, velCompsMask, Mtx.MPI_Comm);

                

                var Poisson = BlockMsrMatrix.Multiply(VelDiv, PrGrad);
                PressureSolver = new PARDISOSolver();
                PressureSolver.DefineMatrix(Poisson);
                RHSmatrix = BlockMsrMatrix.Multiply(VelDiv, CD);
            }
        }

        BlockMask velCompsMask;
        
        BlockMask prsCompsMask;
        
        PARDISOSolver PressureSolver;

        BlockMsrMatrix RHSmatrix;
        BlockMsrMatrix VelDiv;


        public void ResetStat() {
            throw new NotImplementedException();
        }

        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> //
        {
            using(new FuncTrace()) {

                var comm = VelDiv.MPI_Comm;

                int Lvel = velCompsMask.LocalLength;
                int Lprs = prsCompsMask.LocalLength;


                double[] B_vel = velCompsMask.GetSubVec(B);
                double[] RHSpoisson = new double[Lprs];
                VelDiv.SpMV(1.0, B_vel, 0.0, RHSpoisson);

                double[] X_vel = velCompsMask.GetSubVec(X);
                if(X_vel.MPI_L2Norm(comm) > 0) {
                    throw new NotImplementedException();
                }

                double[] Pres = new double[Lprs];
                PressureSolver.Solve(Pres, RHSpoisson);

                prsCompsMask.AccSubVec(Pres, X);

            }
        }

        public long UsedMemory() {
            throw new NotImplementedException();
        }
    }



    class VelocityPredictionSolver : ISolverSmootherTemplate {
        public int IterationsInNested => throw new NotImplementedException();

        public int ThisLevelIterations => throw new NotImplementedException();

        public bool Converged => throw new NotImplementedException();

        public object Clone() {
            throw new NotImplementedException();
        }

        public void Dispose() {
            if (VelocitySolver != null) {
                VelocitySolver.Dispose();
                VelocitySolver = null;
            }

            velCompsMask = null;
            prsCompsMask = null;
        }

        public void Init(MultigridOperator op) {
            using (new FuncTrace()) {
                ICoordinateMapping map = op.Mapping;
                var Mtx = op.OperatorMatrix;

                int D = op.Mapping.SpatialDimension;

                SubBlockSelector velComps = new SubBlockSelector(map);
                velComps.SetVariableSelector(D.ForLoop(iVar => iVar));
                velCompsMask = new BlockMask(velComps);

                SubBlockSelector prsComps = new SubBlockSelector(map);
                prsComps.SetVariableSelector(new int[] { D });
                prsCompsMask = new BlockMask(prsComps);


                // convection-diffusion part
                var CD = velCompsMask.GetSubBlockMatrix(Mtx, Mtx.MPI_Comm);

                VelocitySolver = new PARDISOSolver();
                VelocitySolver.DefineMatrix(CD);
            }
        }

        BlockMask velCompsMask;

        BlockMask prsCompsMask;

        //PARDISOSolver PressureSolver;

        //BlockMsrMatrix RHSmatrix;
        //BlockMsrMatrix VelDiv;

        PARDISOSolver VelocitySolver;


        public void ResetStat() {
            throw new NotImplementedException();
        }

        public void Solve<U, V>(U X, V B)
            where U : IList<double>
            where V : IList<double> //
        {
            using (new FuncTrace()) {

                var comm = VelocitySolver.MpiComm;

                int Lvel = velCompsMask.LocalLength;
                int Lprs = prsCompsMask.LocalLength;


                double[] B_vel = velCompsMask.GetSubVec(B);
                double[] RHSpoisson = B_vel;

                double[] X_prs = prsCompsMask.GetSubVec(X);
                if (X_prs.MPI_L2Norm(comm) > 0) {
                    throw new NotImplementedException();
                }

                double[] Vel = new double[Lvel];
                VelocitySolver.Solve(Vel, RHSpoisson);

                velCompsMask.AccSubVec(Vel, X);

            }
        }

        public long UsedMemory() {
            throw new NotImplementedException();
        }
    }
}
