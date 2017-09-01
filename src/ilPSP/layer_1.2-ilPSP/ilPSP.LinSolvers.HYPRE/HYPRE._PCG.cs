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
//using System.Reflection;
using System.Runtime.InteropServices;
using System.Collections.Generic;
using MPI.Wrappers;


namespace ilPSP.LinSolvers.HYPRE {
    /// <summary>
    /// PCG Solver
    /// </summary>
    public class PCG : Solver, IImplicitPrecondSupport {        
        /// <summary>
        /// creates the PCG solver object
        /// </summary>
        protected override void CreateSolver() {
            HypreException.Check(Wrappers.ParCSRPCG.Create(csMPI.Raw._COMM.WORLD, out m_Solver));
        }

        /// <summary>
        /// calls <see cref="Dispose"/>
        /// </summary>
        ~PCG() {
            Dispose();
        }


        /// <summary>
        /// destroys the hypre PCG solver object, and the nested preconditioner,
        /// if not allready done;
        /// </summary>
        public override void Dispose() {
            if (m_Solver.p != IntPtr.Zero) {
                IDisposable pc = m_NestedPrecond as IDisposable;
                if (pc != null)
                    pc.Dispose();

                HypreException.Check(Wrappers.ParCSRPCG.Destroy(m_Solver));
                m_Solver.p = IntPtr.Zero;
            }
            base.Dispose();
        }



        /// <summary>
        /// turns two-norm residual on/off;
        /// </summary>
        public bool TwoNorm {
            set {
                int b = value ? 1:0;
                HypreException.Check(Wrappers.PCG.SetTwoNorm(m_Solver, b));
            }
            get {
                int r;
                HypreException.Check(Wrappers.PCG.GetTwoNorm(m_Solver, out r));
                return (r!=0);
            }
        }  

        ///// <summary>
        ///// sets a <see cref="BoomerAMG"/> solver as preconditioner for this method
        ///// </summary>
        ///// <param name="precond"></param>
        //private void SetPrecond(BoomerAMG precond) {

            

        //    //Type t = typeof(Wrappers.BoomerAMG);
        //    //MethodInfo mi = t.GetMethod("HYPRE_BoomerAMGSetup");
        //    //object[] attribs = mi.GetCustomAttributes(true);


        //    //solve_setup solve = Wrappers.BoomerAMG.HYPRE_BoomerAMGSolve;
        //    //solve_setup setup = Wrappers.BoomerAMG.HYPRE_BoomerAMGSetup;
        //    //solve_setup2 setup = boomersetup;
            



        //    //Wrappers.PCG.SetPrecond(m_PCGSolver,
        //    //                        solve,
        //    //                        Marshal.GetFunctionPointerForDelegate(setup),
        //    //                        precond.HypreSolver);
        //    IntPtr name = Marshal.StringToHGlobalAnsi("boomeramg");
        //    int r = Wrappers.PCG.SetPrecond(m_PCGSolver,
        //        name,
        //        precond.HypreSolver);
        //    Marshal.FreeHGlobal(name);
        //    if (r == -1)
        //        throw new ApplicationException("internal error");

        //}


        /// <summary>
        /// get/set residual treshold for convergence criterion.
        /// wether this threshold means an absolute or a relative residual
        /// is determined by <see cref="Solver.ConvergenceType"/>;
        /// </summary>
        public double Tolerance {
            get {
                double Tol;
                if (_ConvType == ConvergenceTypes.Relative)
                    HypreException.Check(Wrappers.PCG.GetTol(m_Solver, out Tol));
                else if (_ConvType == ConvergenceTypes.Absolute)
                    HypreException.Check(Wrappers.PCG.GetAbsoluteTol(m_Solver, out Tol));
                else throw new ApplicationException("internal error; unknown convergence type.");

                return Tol;
            }
            set {
                double Tol = value;
                if (_ConvType == ConvergenceTypes.Relative) {
                    HypreException.Check(Wrappers.PCG.SetTol(m_Solver, Tol));
                    HypreException.Check(Wrappers.PCG.SetAbsoluteTol(m_Solver, 0.0));
                } else if (_ConvType == ConvergenceTypes.Absolute) {
                    HypreException.Check(Wrappers.PCG.SetTol(m_Solver, 0.0));
                    HypreException.Check(Wrappers.PCG.SetAbsoluteTol(m_Solver, Tol));
                } else throw new ApplicationException("internal error; unknown convergence type.");

            }
        }

        /// <summary>
        /// get/set maximum number of solver iterations
        /// </summary>
        public int MaxIterations {
            get {
                int MaxIter;
                if (m_Solver.p == IntPtr.Zero)
                    throw new ApplicationException("solver not initialized");

                HypreException.Check(Wrappers.PCG.GetMaxIter(m_Solver, out MaxIter));
                return MaxIter;
            }
            set {
                if (m_Solver.p == IntPtr.Zero)
                    throw new ApplicationException("solver not initialized");


                int MaxIter = value;
                HypreException.Check(Wrappers.PCG.SetMaxIter(m_Solver, MaxIter));
            }
        }

        /// <summary>
        /// turns hypre-internal printing of PCG informative messages
        /// on/off;
        /// </summary>
        public int PrintLevel {
            get {
                int d;
                HypreException.Check(Wrappers.PCG.GetPrintLevel(m_Solver,out d));
                return d;
            }
            set {
                int d = value;
                if (d < 0 || d > 2)
                    throw new ArgumentOutOfRangeException();
                HypreException.Check(Wrappers.PCG.SetPrintLevel(m_Solver, d));
            }
        }

        /// <summary>
        /// calls the PCG solver
        /// </summary>
        /// <param name="Converged">true if converged</param>
        /// <param name="NoOfIter">no of iterations done by solver</param>
        /// <param name="Unknowns"></param>
        /// <param name="Rhs"></param>
        protected override void CallSolver(out int NoOfIter, out bool Converged, IJVector Unknowns, IJVector Rhs) {
            if (m_Solver.p == IntPtr.Zero)
                throw new ApplicationException("solver not initialized");

            HypreException.Check(Wrappers.PCG.Setup(m_Solver, m_Matrix.m_ParCSR_matrix, Rhs.ParCRS_vector, Unknowns.ParCRS_vector));
            Wrappers.PCG.Solve(m_Solver, m_Matrix.m_ParCSR_matrix, Rhs.ParCRS_vector, Unknowns.ParCRS_vector);
            
            // We don't want to raise an exception for
            // a 'method did not converge'- or 'numerical breakdown' - Error
            // that may occur in HYPRE
            Wrappers.Utilities.HYPRE_ClearAllErrors();                                            

            int Con;
            HypreException.Check(Wrappers.PCG.GetConverged(m_Solver, out Con));
            if (Con != 0) Converged =  true;
            else Converged = false;

            HypreException.Check(Wrappers.PCG.GetNumIterations(m_Solver, out NoOfIter));
        }


        /// <summary>
        /// <see cref="IImplicitPrecondSupport.SupportedPrecond"/>;
        /// </summary>
        public ICollection<Type> SupportedPrecond {
            get {
                return new Type[] { typeof(BoomerAMG), typeof(Euclid), typeof(ParaSails) };
            }
        }

        IImplicitPrecond m_NestedPrecond;
        
        /// <summary>
        /// <see cref="IImplicitPrecondSupport.SupportedPrecond"/>;
        /// </summary>
        public IImplicitPrecond NestedPrecond {
            get {
                return (IImplicitPrecond)m_NestedPrecond;
            }
            set {
                ICollection<Type> Supp = this.SupportedPrecond;
                if (Supp.Contains(value.GetType())) {
                    
                    Solver slv = (Solver) value;

                    if (slv.m_NativeSolverFuncPtr == IntPtr.Zero)
                        throw new ApplicationException("internal error: 'Solver'-function delegate on preconditioner not set.");
                    if (slv.m_NativeSetupFuncPtr == IntPtr.Zero)
                        throw new ApplicationException("internal error: 'Setup'-function delegate on preconditioner not set.");
                    
                    HypreException.Check(
                        Wrappers.PCG.SetPrecond(m_Solver,
                        slv.m_NativeSolverFuncPtr, slv.m_NativeSetupFuncPtr,
                        slv.m_Solver));

                } else {
                    throw new ArgumentException("type '" + value.GetType().Name + "' is not supported as preconditioner for Hyper PCG solver.");
                }
                m_NestedPrecond = (IImplicitPrecond)value;
            }
        }

    }
}
