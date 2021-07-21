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
using System.Reflection;
using System.Runtime.InteropServices;
using MPI.Wrappers;
using System.Collections.Generic;

namespace ilPSP.LinSolvers.HYPRE
{
    /// <summary>
    ///  GMRES Solver
    /// </summary>
    public class GMRES : Solver, IImplicitPrecondSupport
    {
        /// <summary>
        /// create GMRES solver object
        /// </summary>
        protected override void CreateSolver()
        {
            HypreException.Check(Wrappers.ParCSRGMRES.CreateGMRES(csMPI.Raw._COMM.WORLD, out m_Solver));
            PrintLevel = 0;
        }

        /// <summary>
        /// calls <see cref="Dispose"/>
        /// </summary>
        ~GMRES() {
            Dispose();
        }

        /// <summary>
        /// destroys the HYPRE object, if not already done;
        /// </summary>
        public override void Dispose() {
            if (m_Solver.p != IntPtr.Zero) {
                HypreException.Check(Wrappers.ParCSRGMRES.Destroy(m_Solver));
                m_Solver.p = IntPtr.Zero;
            }
            base.Dispose();
        }


        /// <summary>
        /// Setup absolute or relative convergence tolerance
        /// </summary>
        public double Tolerance
        {           
            get
            {
                double _Tolerance;
                switch (_ConvType)
                {
                    case ConvergenceTypes.Relative:
                        HypreException.Check(Wrappers.GMRES.GetRelativeTolerance(m_Solver, out _Tolerance));
                        break;
                    case ConvergenceTypes.Absolute:
                        HypreException.Check(Wrappers.GMRES.GetAbsoluteTolerance(m_Solver, out _Tolerance));
                        break;
                    default:
                        throw new ApplicationException("internal error; unknown convergence type.");
                }          
                return _Tolerance;
            }
            set
            {
                switch (_ConvType)
                {
                    case ConvergenceTypes.Relative:
                        HypreException.Check(Wrappers.GMRES.SetRelativeTolerance(m_Solver, value));
                        break;
                    case ConvergenceTypes.Absolute:
                        HypreException.Check(Wrappers.GMRES.SetAbsoluteTolerance(m_Solver, value));
                        HypreException.Check(Wrappers.GMRES.SetRelativeTolerance(m_Solver, 0.0));
                        break;
                    default:
                        throw new ApplicationException("internal error; unknown convergence type.");
                }
            }
        }

        /// <summary>
        /// (Optional)Set and get maximal number of iterations for the solver;
        /// </summary>
        public int MaxIterations
        {
            get
            {
                int _MaxIterations;
                HypreException.Check(Wrappers.GMRES.GetMaxIterations(m_Solver, out _MaxIterations));
                return _MaxIterations;
            }
            set
            {
                HypreException.Check(Wrappers.GMRES.SetMaxIterations(m_Solver, value));
            }
        }
        /// <summary>
        ///  Set the minimal number of iterations
        /// </summary>
        public int MinIterations
        {
            get
            {
                int _MinIterations;
                HypreException.Check(Wrappers.GMRES.GetMinIterations(m_Solver, out _MinIterations));
                return _MinIterations;
            }
            set {
                HypreException.Check(Wrappers.GMRES.SetMinIterations(m_Solver, value));
            }
        }
        /// <summary>
        /// (Optional) Set and get maximum size of Krylov  space
        /// </summary>
        public int KrylovSpaceMaxDimension
        {
            get {
                int _KrylovSpaceMaxDimension;
                HypreException.Check(Wrappers.GMRES.GetKrylovSpaceDim(m_Solver, out _KrylovSpaceMaxDimension));
                return _KrylovSpaceMaxDimension;
            }
            set {
                HypreException.Check(Wrappers.GMRES.SetKrylovSpaceDim(m_Solver, value));
            }
        }

        /// <summary>
        ///  (Optional) Additionally require that the relative difference in successive iterates be small
        /// </summary>
        public int RelChange
        {
            get {
                int _RelChange;
                HypreException.Check(Wrappers.GMRES.GetRelChange(m_Solver, out _RelChange));
                return _RelChange;
            }
            set {
                HypreException.Check(Wrappers.GMRES.SetRelChange(m_Solver, value));
            }
        }

        int m_PrintLevel = 0;

        /// <summary>
        /// turns hypre-internal printing of GMRES informative messages
        /// on/off;
        /// </summary>
        public int PrintLevel {
            get {
                //int d;
                //HypreException.Check(Wrappers.GMRES.GetPrintLevel(m_Solver, out d));
                return m_PrintLevel;
            }
            set {
                int d = value;
                if (d < 0 || d > 2)
                    throw new ArgumentOutOfRangeException();
                HypreException.Check(Wrappers.GMRES.SetPrintLevel(m_Solver, d));
                m_PrintLevel = d;
            }
        }

        /// <summary>
        /// Solve equation
        /// </summary>
        /// <param name="NoOfIter">placeholder for number of iteration taken</param>
        /// <param name="Converged">placeholder for convergence state</param>
        /// <param name="Unknowns"></param>
        /// <param name="Rhs"></param>
        protected override void CallSolver(out int NoOfIter, out bool Converged, IJVector Unknowns, IJVector Rhs)
        {
            if (m_Solver.p == IntPtr.Zero)
                throw new ApplicationException("solver not initialized");

            

            HypreException.Check(Wrappers.GMRES.Setup(m_Solver, m_Matrix.m_ParCSR_matrix, Rhs.ParCRS_vector, Unknowns.ParCRS_vector));
            Wrappers.GMRES.Solve(m_Solver, m_Matrix.m_ParCSR_matrix, Rhs.ParCRS_vector, Unknowns.ParCRS_vector);
            // We don't want to raise an exception for
            // a 'method did not converge'- or 'numerical breakdown' - Error
            // that may occur in HYPRE
            Wrappers.Utilities.HYPRE_ClearAllErrors();                                            


            int _Converged;
            HypreException.Check(Wrappers.GMRES.GetConverged(m_Solver, out _Converged));
            Converged = _Converged == 0 ? false : true;

            HypreException.Check(Wrappers.GMRES.GetNumOfIterations(m_Solver, out NoOfIter)); 
         }

        #region INestedSolver Members

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
                return m_NestedPrecond;
            }
            set {
                ICollection<Type> Supp = this.SupportedPrecond;
                if (Supp.Contains(value.GetType())) {

                    Solver slv = (Solver)value;

                    if (slv.m_NativeSolverFuncPtr == IntPtr.Zero)
                        throw new ApplicationException("internal error: 'Solver'-function delegate on preconditioner not set.");
                    if (slv.m_NativeSetupFuncPtr == IntPtr.Zero)
                        throw new ApplicationException("internal error: 'Setup'-function delegate on preconditioner not set.");

                    int r = Wrappers.GMRES.SetPrecond(m_Solver,
                        slv.m_NativeSolverFuncPtr, slv.m_NativeSetupFuncPtr,
                        slv.m_Solver);

                    if (r == -1)
                        throw new ApplicationException("HYPRE error: return code " + r);
                } else {
                    throw new ArgumentException("type '" + value.GetType().Name + "' is not supported as preconditioner for Hyper PCG solver.");
                }
                m_NestedPrecond = value;
            }
        }


        #endregion
    }
}
