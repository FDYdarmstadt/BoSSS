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
using System.Collections.Generic;
using System.IO;
using System.Globalization;


namespace ilPSP.LinSolvers {

    /// <summary>
    /// returned status info for calls to linear sparse solvers
    /// </summary>
    public class SolverResult {

        /// <summary>
        /// runtime of the solver (wall-time);
        /// This is the only performance criterion to compare different solvers;
        /// </summary>
        public TimeSpan RunTime;

        /// <summary>
        /// number of iterations the solver made
        /// </summary>
        public int NoOfIterations;

        /// <summary>
        /// if this is true, the solver has finished successfully
        /// </summary>
        public bool Converged;

        /// <summary>
        /// prints the statistics to a string
        /// </summary>
        /// <returns></returns>
        public override string ToString() {
            return ToString(NumberFormatInfo.CurrentInfo);
        }

        /// <summary>
        /// prints the statistics to a string
        /// </summary>
        /// <param name="fi">
        /// format provider for floating point numbers in the output
        /// </param>
        public string ToString(NumberFormatInfo fi) {
            StringWriter stw = new StringWriter();
            stw.Write("Converged: ");
            stw.Write(Converged.ToString(fi));
            stw.Write(", Number of Iter.: ");
            stw.Write(NoOfIterations);
            stw.Write(", Runtime: ");
            stw.Write(RunTime.TotalSeconds.ToString(fi));
            stw.Write(";");

            return stw.ToString();

        }
    }

    /// <summary>
    /// Termination options for <see cref="ISparseSolver"/>;
    /// </summary>
    public enum ConvergenceTypes {
        
        /// <summary>
        /// terminate iterations if an absolute residual norm criterion has
        /// been reached;
        /// </summary>
        Absolute,

        /// <summary>
        /// terminate iterations if an relative residual norm criterion has
        /// been reached;
        /// </summary>
        Relative,

        /// <summary>
        /// solver - specific;
        /// </summary>
        Other
    }

    /// <summary>
    /// This kind of preconditioners alter the matrix (of the linear system)
    /// and pass an altered system to a so-called 'nested solver' (see <see cref="NestedSolver"/>).
    /// The explicit preconditioner looks like a sparse solver itself 
    /// (not that it also implements the <see cref="ISparseSolver"/> - interface).
    /// </summary>
    /// <remarks>
    /// Mathematically, 
    /// for an explicit preconditioner, the preconditioning
    /// matrix is known explicity.
    /// Technically, such an precond. coud be eiter implemented as
    /// <see cref="IImplicitPrecond"/> or as <see cref="IExplicitPrecond"/>
    /// with both versions offering situation-dependent advantaged and disadvantages.
    /// </remarks>
    /// <seealso cref="IImplicitPrecond"/>
    public interface IExplicitPrecond : ISparseSolver {

        /// <summary>
        /// Defines the solver that is used for solving the altered system;
        /// </summary>
        ISparseSolver NestedSolver { get; set; }
    }

    /// <summary>
    /// This kind of preconditioners are nested in sparse solvers,
    /// i.e. their operation is applied in every iteration
    /// of the outer solver;<br/>
    /// A solver that wants to use such an object must implement the 
    /// <see cref="IImplicitPrecondSupport"/>-interface;
    /// </summary>
    public interface IImplicitPrecond {

    }


    /// <summary>
    /// Impelented by solvers that support implicit preconditioners 
    /// (see <see cref="IImplicitPrecond"/>,nested preconditioners), i.e.
    /// solvers which can use other solvers as preconditioners.
    /// </summary>
    public interface IImplicitPrecondSupport : ISparseSolver {

        /// <summary>
        /// Types of all implicit prconditioners (classes that implement <see cref="IImplicitPrecond"/>)
        /// which can be used with this solver
        /// </summary>
        ICollection<Type> SupportedPrecond { get; }

        /// <summary>
        /// Sets the preconditioner; The type of the value (for setting)
        /// must be contained in the <see cref="SupportedPrecond"/>-collection;
        /// </summary>
        IImplicitPrecond NestedPrecond { get; set; }
    }

    
    /// <summary>
    /// A common minimal interface for linear (sparse) solvers
    /// (e.g. Aztec, PETSc, hypre), higher code levels build on;
    /// (The abstraction of solvers with this interface makes them exchangeable 
    /// in higher layers)
    /// </summary>
    public interface ISparseSolver : IDisposable {

        /// <summary>
        /// defines the matrix of the equation system that should be solved;
        /// (usually, implementers of that interface need to convert the matrix to some
        /// solver-dependent format.);
        /// </summary>
        /// <param name="M"></param>
        /// <remarks>
        /// can be called only once in the lifetime of an <see cref="ISparseSolver"/>-object
        /// </remarks>
        void DefineMatrix(IMutableMatrixEx M);

        /// <summary>
        /// solves the equation system 
        /// M*<paramref name="x"/>=<paramref name="rhs"/>,
        /// where M denotes the <see cref="MsrMatrix"/>
        /// that was provided by <see cref="DefineMatrix"/>;
        /// </summary>
        /// <typeparam name="Tunknowns"></typeparam>
        /// <typeparam name="Trhs"></typeparam>
        /// <param name="x">vector for the unknowns</param>
        /// <param name="rhs">right hand side</param>
        SolverResult Solve<Tunknowns, Trhs>(Tunknowns x, Trhs rhs)
            where Tunknowns : IList<double>
            where Trhs : IList<double>;
    }

    /// <summary>
    /// Extended interface for sparse solvers, where the matrix values (not the occupation)
    /// can be altered after matrix assembly (i.e. after calling <see cref="ISparseSolver.DefineMatrix"/>);<br/>
    /// This is necessary e.g. for solving the Heat Equation with implicit Euler: 
    /// \f[  (id - \Delta t \ M) \cdot u^{n+1} = u^n. \f]
    /// </summary>
    public interface ISparseSolverExt : ISparseSolver {

        /// <summary>
        /// returns an interface to the internally used matrix,
        /// exposed as some <see cref="ISparseMatrix"/>-implementation.
        /// The entries of this matrix are equal to the entries of those matrix which was provided
        /// in <see cref="ISparseSolver.DefineMatrix"/>.
        /// </summary>
        /// <returns></returns>
        ISparseMatrix GetMatrix();


        /// <summary>
        /// solves the equation system 
        /// (diag(<paramref name="d"/>) + <paramref name="Scale"/>*M)*<paramref name="x"/>=<paramref name="rhs"/>,
        /// where diag(<paramref name="d"/>) denotes a diagonal matrix with the diagonal vector <paramref name="d"/>
        /// and  M denotes the <see cref="MsrMatrix"/>
        /// that was provided by <see cref="ISparseSolver.DefineMatrix"/>;
        /// </summary>
        /// <typeparam name="Tunknowns"></typeparam>
        /// <typeparam name="Trhs"></typeparam>
        /// <typeparam name="Tdiag"></typeparam>
        /// <param name="Scale">
        /// scaling of the matrix (does not apply to <paramref name="d"/>);
        /// must be unequal to 0;
        /// this scaling does not alter the matrix-it's the same after the return of this method.
        /// </param>
        /// <param name="d">
        /// diagonal vector; can be null;
        /// length must be equal to the Nupdate-length of the matrix mapping; 
        /// may be altered during execution;
        /// </param>
        /// <param name="x">on exit, the approximate solution to the equation system;
        /// </param>
        /// <param name="rhs">
        /// right hand side of the equation system; may be overwritten;
        /// </param>
        SolverResult Solve<Tdiag, Tunknowns, Trhs>(double Scale, Tdiag d, Tunknowns x, Trhs rhs)
            where Tunknowns : IList<double>
            where Tdiag : IList<double>
            where Trhs : IList<double>;
    }
}
