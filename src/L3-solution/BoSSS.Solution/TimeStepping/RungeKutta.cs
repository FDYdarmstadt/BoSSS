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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using ilPSP.Utils;
using System;
using System.Linq;
using System.Collections.Generic;

namespace BoSSS.Solution.Timestepping {

    /// <summary>
    /// Runge-Kutta time stepping
    /// </summary>
    public class RungeKutta : ExplicitEuler {

        /// <summary>
        /// Supported Runge-Kutta schemes. Using this enum allows for a simple
        /// instantiation of a selected scheme from a configuration string by
        /// means of
        /// <see cref="RungeKutta.SchemeFactory"/>.
        /// </summary>
        public enum RungeKuttaSchemes {

            /// <summary>
            /// <see cref="RungeKutta.RungeKuttaScheme.ExplicitEuler"/>
            /// </summary>
            ExplicitEuler,

            /// <summary>
            /// <see cref="RungeKutta.RungeKuttaScheme.RungeKutta1901"/>
            /// </summary>
            RungeKutta1901,

            /// <summary>
            /// <see cref="RungeKutta.RungeKuttaScheme.Heun"/>
            /// </summary>
            Heun,

            /// <summary>
            /// <see cref="RungeKutta.RungeKuttaScheme.Middlepoint"/>
            /// </summary>
            Middlepoint,

            /// <summary>
            /// <see cref="RungeKutta.RungeKuttaScheme.TVD3"/>
            /// </summary>
            TVD3,

            /// <summary>
            /// <see cref="RungeKutta.RungeKuttaScheme.ThreeOverEight"/>
            /// </summary>
            ThreeOverEight,

            /// <summary>
            /// <see cref="RungeKutta.RungeKuttaScheme.SSP54"/>
            /// </summary>
            SSP54
        }

        /// <summary>
        /// one more ctor.
        /// </summary>
        public RungeKutta(RungeKuttaSchemes schemeEnum, SpatialOperator spatialOp, CoordinateMapping Fieldsmap, CoordinateMapping Parameters, IList<TimeStepConstraint> timeStepConstraints = null, SubGrid sgrd = null)
            : this(SchemeFactory(schemeEnum), spatialOp, Fieldsmap, Parameters, timeStepConstraints, sgrd) {
        }

        /// <summary>
        /// yet another ctor.
        /// </summary>
        public RungeKutta(RungeKuttaScheme scheme, SpatialOperator spatialOp, CoordinateMapping Fieldsmap, CoordinateMapping Parameters, SubGrid sgrd=null)
            : this(scheme, spatialOp, Fieldsmap, Parameters, null, sgrd) {
        }

        /// <summary>
        /// another ctor.
        /// </summary>
        public RungeKutta(RungeKuttaScheme scheme, SpatialOperator spatialOp, params DGField[] Fields)
            : this(scheme, spatialOp, new CoordinateMapping(Fields),null, null) {
        }

        /// <summary>
        /// another ctor.
        /// </summary>
        public RungeKutta(RungeKuttaSchemes scheme, SpatialOperator spatialOp, params DGField[] Fields)
            : this(scheme, spatialOp, new CoordinateMapping(Fields),null, null) {
        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="spatialOp">Spatial operator</param>
        /// <param name="Fieldsmap"></param>
        /// <param name="scheme">Runge-Kutta scheme</param>
        /// <param name="Parameters">
        /// Optional parameter fields, can be null if
        /// <paramref name="spatialOp"/> contains no parameters;
        /// must match the parameter field list of
        /// <paramref name="spatialOp"/>, see
        /// <see cref="BoSSS.Foundation.SpatialOperator.ParameterVar"/>
        /// </param>
        /// <param name="timeStepConstraints">
        /// optional list of time step constraints <see cref="TimeStepConstraint"/>
        /// </param>
        /// <param name="sgrd">
        /// optional restriction to computational domain
        /// </param>
        public RungeKutta(RungeKuttaScheme scheme, SpatialOperator spatialOp, CoordinateMapping Fieldsmap, CoordinateMapping Parameters, IList<TimeStepConstraint> timeStepConstraints, SubGrid sgrd = null)
            : base(spatialOp, Fieldsmap, Parameters, timeStepConstraints, sgrd) {
            using (new ilPSP.Tracing.FuncTrace()) {
                m_Scheme = scheme;

                m_Scheme.Verify();
                if (!m_Scheme.IsExplicit)
                    throw new ArgumentException("Implicit Runge-Kutta -- schemes are not supported.");
            }
        }

        /// <summary>
        /// Returns the default Runge-Kutta scheme for a given order
        /// </summary>
        /// <param name="order">
        /// The desired order of the scheme.
        /// </param>
        /// <returns>
        /// Depending on <paramref name="order"/>, the following scheme will be
        /// chosen:
        /// <list type="bullet">
        ///     <item>
        ///         <term>Order 1</term>
        ///         <description>
        ///             <see cref="RungeKuttaScheme.ExplicitEuler"/>
        ///         </description>
        ///     </item>
        ///     <item>
        ///         <term>Order 2</term>
        ///         <description>
        ///             <see cref="RungeKuttaScheme.Middlepoint"/>
        ///         </description>
        ///     </item>
        ///     <item>
        ///         <term>Order 3</term>
        ///         <description>
        ///             <see cref="RungeKuttaScheme.Heun"/>
        ///         </description>
        ///     </item>
        ///     <item>
        ///         <term>Order 4</term>
        ///         <description>
        ///             <see cref="RungeKuttaScheme.RungeKutta1901"/>
        ///         </description>
        ///     </item>
        /// </list>
        /// </returns>
        public static RungeKuttaSchemes GetDefaultScheme(int order) {
            switch (order) {
                case 1:
                    return RungeKuttaSchemes.ExplicitEuler;

                case 2:
                    return RungeKuttaSchemes.Middlepoint;

                case 3:
                    return RungeKuttaSchemes.Heun;

                case 4:
                    return RungeKuttaSchemes.RungeKutta1901;

                default:
                    throw new ArgumentException(
                        "Order not supported. Order must be between 1 and 4, but was " + order, "order");
            }
        }

        /// <summary>
        /// Utility method that helps creating a Runge-Kutta scheme from the
        /// enum <see cref="RungeKuttaSchemes"/>. Allows for a simple
        /// instantiation of the correct scheme from a configuration string.
        /// </summary>
        /// <param name="scheme"></param>
        /// <returns>
        /// An instance of a Runge-Kutta scheme the selected
        /// <paramref name="scheme"/>.
        /// </returns>
        public static RungeKuttaScheme SchemeFactory(RungeKuttaSchemes scheme) {
            switch (scheme) {
                case RungeKuttaSchemes.ExplicitEuler:
                    return RungeKuttaScheme.ExplicitEuler;

                case RungeKuttaSchemes.Heun:
                    return RungeKuttaScheme.Heun;

                case RungeKuttaSchemes.Middlepoint:
                    return RungeKuttaScheme.Middlepoint;

                case RungeKuttaSchemes.RungeKutta1901:
                    return RungeKuttaScheme.RungeKutta1901;

                case RungeKuttaSchemes.ThreeOverEight:
                    return RungeKuttaScheme.ThreeOverEight;

                case RungeKuttaSchemes.TVD3:
                    return RungeKuttaScheme.TVD3;

                default:
                    throw new NotImplementedException();
            }
        }

        /// <summary>
        /// the scheme
        /// </summary>
        RungeKuttaScheme m_Scheme;

        /// <summary>
        /// the scheme
        /// </summary>
        public RungeKuttaScheme Scheme {
            get {
                return m_Scheme;
            }
        }

        /// <summary>
        /// performs one timestep
        /// </summary>
        /// <param name="dt">size of timestep</param>
        public override double Perform(double dt) {

            using (var tr = new ilPSP.Tracing.FuncTrace()) {

                if (TimeStepConstraints != null) {
                    dt = CalculateTimeStep();
                }

                double[][] k = new double[m_Scheme.Stages][];
                for (int i = 0; i < m_Scheme.Stages; i++)
                    k[i] = new double[Mapping.LocalLength];

                double[] y0 = new double[Mapping.LocalLength];
                CurrentState.CopyTo(y0, 0);

                // logging
                tr.Info("Runge-Kutta Scheme with " + m_Scheme.Stages + " stages.");
                double time0 = m_Time;
                tr.Info("time = " + time0 + ", dt = " + dt);

                // berechne k[0]
                ComputeChangeRate(k[0], m_Time, 0);// invokes MPI communication

                for (int s = 1; s < m_Scheme.Stages; s++) {
                    PerformStage(y0, s, k, dt);

                    m_Time = time0 + m_Scheme.c[s] * dt;
                    ComputeChangeRate(k[s], m_Time, m_Scheme.c[s] * dt);// invokes MPI communication
                }

                // next timestep
                CurrentState.Clear();
                CurrentState.CopyFrom(y0, 0);
                Array.Clear(y0, 0, y0.Length);
                for (int s = 0; s < m_Scheme.Stages; s++) {
                    if (m_Scheme.b[s] != 0.0) {
                        BLAS.daxpy(y0.Length, -m_Scheme.b[s] * dt, k[s], 1, y0, 1);
                    }
                }
                CurrentState.axpy<double[]>(y0, 1.0);
                base.ApplyFilter(dt);

                m_Time = time0 + dt;

            }
            return dt;
        }

        /// <summary>
        /// Performs a single stage of a timestep.
        /// </summary>
        /// <param name="y0">Initial DG coordinates</param>
        /// <param name="s">The current stage</param>
        /// <param name="k">The current change rate</param>
        /// <param name="dt">The size of the timestep</param>
        protected virtual void PerformStage(double[] y0, int s, double[][] k, double dt) {
            CurrentState.Clear();
            CurrentState.axpy<double[]>(y0, 1.0);
            double relTime = 0;
            for (int r = 0; r < s; r++) {
                if (m_Scheme.a[s, r] != 0.0) {
                    CurrentState.axpy<double[]>(k[r], -dt * m_Scheme.a[s, r]);
                    relTime += dt * m_Scheme.a[s, r];
                }
            }

            base.ApplyFilter(dt * relTime);
        }


    }
}
