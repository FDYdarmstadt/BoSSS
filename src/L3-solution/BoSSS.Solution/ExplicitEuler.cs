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
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;

namespace BoSSS.Solution.Timestepping {

    /// <summary>
    /// Explicit Euler timestepping with auxiliary condition (that is a
    /// non-time dependent equation).
    /// </summary>
    /// <seealso cref="RungeKutta">
    /// This class also serves as a base-class for Runge-Kutta (and maybe
    /// other) timesteppers
    /// </seealso>
    public class ExplicitEuler : ITimeStepper {

        /// <summary>
        /// The spatial operator this time stepper is based on
        /// </summary>
        public SpatialOperator Operator {
            get;
            private set;
        }

        /// <summary>
        /// mapping of DG fields that are affected by this timestepper
        /// </summary>
        public CoordinateMapping Mapping {
            get;
            private set;
        }

        /// <summary>
        /// mapping of parameter fields
        /// </summary>
        public CoordinateMapping ParameterMapping {
            get;
            private set;
        }

        /// <summary>
        /// The DG coordinates of the <see cref="BoSSS.Foundation.DGField"/>'s
        /// in <see cref="Mapping"/> (see
        /// <see cref="BoSSS.Foundation.CoordinateMapping.Fields"/>).
        /// </summary>
        public CoordinateVector DGCoordinates {
            get;
            private set;
        }

        /// <summary>
        /// List of various time step constraints, e.g CFL condition
        /// </summary>
        public IList<TimeStepConstraint> timeStepConstraints {
            get;
            protected set;
        }

        protected SubGrid Subgrid;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="spatialOp"></param>
        /// <param name="Fields"></param>
        /// <remarks>
        /// This constructor does not support parameter fields; 
        /// </remarks>
        public ExplicitEuler(SpatialOperator spatialOp, params DGField[] Fields)
            : this(spatialOp, new CoordinateMapping(Fields), null) {
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="spatialOp"></param>
        /// <param name="Fieldsmap"></param>
        /// <param name="Parameters"></param>
        /// <param name="timeStepConstraints">
        /// optional list of time step constraints <see cref="TimeStepConstraint"/>
        /// </param>
        /// <param name="sgrd"></param>
        /// <remarks>
        /// This constructor sets
        /// <see cref="SpatialOperator.SubGridBoundaryModes"/> to default value
        /// BoundaryEdge 
        /// </remarks>
        public ExplicitEuler(SpatialOperator spatialOp, CoordinateMapping Fieldsmap, CoordinateMapping Parameters, IList<TimeStepConstraint> timeStepConstraints = null, SubGrid sgrd = null)
            : this(spatialOp, Fieldsmap, Parameters, SpatialOperator.SubGridBoundaryModes.BoundaryEdge, timeStepConstraints, sgrd) {
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="spatialOp"></param>
        /// <param name="Fieldsmap"></param>
        /// <param name="Parameters">
        /// optional parameter fields, can be null if
        /// <paramref name="spatialOp"/> contains no parameters; must match
        /// the parameter field list of <paramref name="spatialOp"/>, see
        /// <see cref="BoSSS.Foundation.SpatialOperator.ParameterVar"/>
        /// </param>
        /// <param name="sgrdBnd">
        /// Options for the treatment of edges at the boundary of a SubGrid,
        /// <see cref="SpatialOperator.SubGridBoundaryModes"/></param>
        /// <param name="timeStepConstraints">
        /// optional list of time step constraints <see cref="TimeStepConstraint"/>
        /// </param>
        /// <param name="sgrd">
        /// optional restriction to computational domain
        /// </param>
        public ExplicitEuler(SpatialOperator spatialOp, CoordinateMapping Fieldsmap, CoordinateMapping Parameters, SpatialOperator.SubGridBoundaryModes sgrdBnd, IList<TimeStepConstraint> timeStepConstraints = null, SubGrid sgrd = null) {
            using (new ilPSP.Tracing.FuncTrace()) {

                // verify input
                // ============

                if (!spatialOp.ContainsNonlinear && !(spatialOp.ContainsLinear()))
                    throw new ArgumentException("spatial differential operator seems to contain no components.", "spatialOp");
                if (spatialOp.DomainVar.Count != spatialOp.CodomainVar.Count)
                    throw new ArgumentException("spatial differential operator must have the same number of domain and codomain variables.", "spatialOp");
                if (Fieldsmap.Fields.Count != spatialOp.CodomainVar.Count)
                    throw new ArgumentException("the number of fields in the coordinate mapping must be equal to the number of domain/codomain variables of the spatial differential operator", "fields");
                if (Parameters == null) {
                    if (spatialOp.ParameterVar.Count != 0)
                        throw new ArgumentException("the number of fields in the parameter mapping must be equal to the number of parameter variables of the spatial differential operator", "Parameters");
                } else {
                    if (Parameters.Fields.Count != spatialOp.ParameterVar.Count)
                        throw new ArgumentException("the number of fields in the parameter mapping must be equal to the number of parameter variables of the spatial differential operator", "Parameters");
                }

                Mapping = Fieldsmap;
                DGCoordinates = new CoordinateVector(Mapping);
                ParameterMapping = Parameters;
                IList<DGField> ParameterFields =
                    (ParameterMapping == null) ? (new DGField[0]) : ParameterMapping.Fields;

                this.timeStepConstraints = timeStepConstraints;
                Subgrid = sgrd;

                // generate Evaluator
                // ==================

                CellMask cm = (sgrd == null) ? null : sgrd.VolumeMask;
                EdgeMask em = (sgrd == null) ? null : sgrd.AllEdgesMask;

                Operator = spatialOp;
                m_Evaluator = new Lazy<SpatialOperator.Evaluator>(() => spatialOp.GetEvaluatorEx(
                    Fieldsmap, ParameterFields, Fieldsmap,
                    new EdgeQuadratureScheme(true, em),
                    new CellQuadratureScheme(true, cm),
                    sgrd,
                    sgrdBnd));
            }
        }

        /// <summary>
        /// see <see cref="Evaluator"/>;
        /// </summary>
        protected Lazy<SpatialOperator.Evaluator> m_Evaluator;

        /// <summary>
        /// Evaluator for the spatial operator;
        /// </summary>
        public SpatialOperator.Evaluator Evaluator {
            get {
                return m_Evaluator.Value;
            }
        }

        /// <summary>
        /// computes the change rate for the time evolution problem;
        /// </summary>
        /// <param name="k">
        /// the "change rate" is <b>accumulated</b> here;<br/>
        /// Indices into this array can be computed by <see cref="Mapping"/>;
        /// </param>
        /// <param name="AbsTime">
        /// absolute time at which the change rate will be evaluated
        /// </param>
        /// <param name="RelTime">
        /// time, relative to the point of time at which the initial value holds.
        /// </param>
        /// <remarks>
        /// Assume an ODE
        /// <br/>
        /// d/dt u = f(u)
        /// <br/>
        /// which is discretized by an explicit Euler scheme
        /// <br/>
        /// (u_1 - u_0)/delta_t = f(u_0).
        /// <br/>
        /// Here, f(u) denotes the discretization of the spatial operator by the DG method.
        /// The purpose of this method is to evaluate f(u_0).<br/>
        /// It does so by calling <see cref="SpatialOperator.Evaluator.Evaluate{Tout}"/>
        /// of <see cref="Evaluator"/>.<br/>
        /// Override this method e.g. for the implementation of (some types of) limiters.
        /// </remarks>
        virtual internal protected void ComputeChangeRate(double[] k, double AbsTime, double RelTime, double[] edgeFluxes = null) {
            using (var tr = new ilPSP.Tracing.FuncTrace()) {
                RaiseOnBeforeComputechangeRate(AbsTime, RelTime);

                // k += F(u0)
                Evaluator.Evaluate<double[]>(1.0, 1.0, k, AbsTime, outputBndEdge: edgeFluxes);
            }
        }

        protected void RaiseOnBeforeComputechangeRate(double AbsTime, double RelTime) {
            if (OnBeforeComputeChangeRate != null) {
                OnBeforeComputeChangeRate(AbsTime, RelTime);
            }
        }

        /// <summary>
        /// The matrix representing the affine linear part of <see cref="Operator"/>.
        /// </summary>
        public MsrMatrix Matrix {
            get;
            private set;
        }

        /// <summary>
        /// The vector representing the offset of the affine linear part of
        /// <see cref="Operator"/>.
        /// </summary>
        public double[] AffineOffset {
            get;
            private set;
        }

        /// <summary>
        /// fired immediately after entering <see cref="ComputeChangeRate"/>; can be used e.g. to update some parameters.
        /// </summary>
        public event ChangeRateCallback OnBeforeComputeChangeRate;

        /// <summary>
        /// fired after the DG fields have been updated; can be used e.g. to apply some filter/limiter.
        /// </summary>
        public event AfterUpdateCallback OnAfterFieldUpdate;

        /// <summary>
        /// see <see cref="OnBeforeComputeChangeRate"/>;
        /// </summary>
        /// <param name="AbsTime">
        /// absolute time, i.e. all timesteps that were performed so far on the current object.
        /// </param>
        /// <param name="RelTime">
        /// relative time since the start of the current timestep.
        /// </param>
        public delegate void ChangeRateCallback(double AbsTime, double RelTime);

        /// <summary>
        /// see <see cref="OnAfterFieldUpdate"/>;
        /// </summary>
        /// <param name="AbsTime">
        /// absolute time, i.e. all timesteps that were performed so far on the current object.
        /// </param>
        /// <param name="Fields">
        /// the DG fields which were updated.
        /// </param>
        public delegate void AfterUpdateCallback(double AbsTime, CoordinateMapping Fields);


        /// <summary>
        /// performs one Explicit Euler timestep
        /// </summary>
        /// <param name="dt">size of timestep</param>
        public virtual double Perform(double dt) {
            using (new ilPSP.Tracing.FuncTrace()) {
                if (timeStepConstraints != null) {
                    dt = CalculateTimeStep();
                }
                double[] k = new double[Mapping.LocalLength];
                ComputeChangeRate(k, m_Time, 0);
                DGCoordinates.axpy<double[]>(k, -dt);

                ApplyFilter(dt);

                m_Time += dt;
            }
            return dt;
        }

        /// <summary>
        /// Calculates a stable time step according to the time step constraints
        /// The individual constraints are connected by a harmonic sum 
        /// </summary>
        /// <param name="dt"></param>
        /// <returns></returns>
        virtual protected double CalculateTimeStep() {
            double dt;
            if (timeStepConstraints.First().dtMin != timeStepConstraints.First().dtMax) {
                // Use "harmonic sum" of step - sizes, see
                // WatkinsAsthanaJameson2016 for the reasoning
                dt = 1.0 / timeStepConstraints.Sum(
                        c => 1.0 / c.GetGloballyAdmissibleStepSize(Subgrid));
                if (dt == 0.0) {
                    throw new ArgumentException(
                        "Time-step size is exactly zero.");
                } else if (double.IsNaN(dt)) {
                    throw new ArgumentException(
                        "Could not determine stable time-step size. This indicates illegal values in some cells.");
                }

                dt = Math.Min(dt, timeStepConstraints.First().Endtime - Time);
                dt = Math.Min(Math.Max(dt, timeStepConstraints.First().dtMin), timeStepConstraints.First().dtMax);
            } else {
                dt = timeStepConstraints.First().dtMin;
                dt = Math.Min(dt, timeStepConstraints.First().Endtime - Time);
            }

            return dt;
        }

        /// <summary>
        /// invokes the <see cref="OnAfterFieldUpdate"/>- event
        /// </summary>
        protected void ApplyFilter(double dt) {
            if (OnAfterFieldUpdate != null) {
                OnAfterFieldUpdate(m_Time + dt, this.DGCoordinates.Mapping);
            }
        }

        #region ITimeStepper Members

        /// <summary>
        /// <see cref="ITimeStepper.Time"/>
        /// </summary>
        public double Time {
            get {
                return m_Time;
            }
        }

        /// <summary>
        /// physical time 
        /// </summary>
        protected double m_Time = 0.0;

        /// <summary>
        /// <see cref="ITimeStepper.ResetTime"/>
        /// </summary>
        /// <param name="NewTime"></param>
        public virtual void ResetTime(double NewTime) {
            m_Time = NewTime;
        }

        #endregion
    }
}
