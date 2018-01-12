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
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using ilPSP.Utils;

namespace BoSSS.Solution.Timestepping {

    /// <summary>
    /// Adams-Bashforth time-stepping;
    /// </summary>
    public class AdamsBashforth : ExplicitEuler {

        /// <summary>
        /// Stores History of time steps
        /// </summary>
        protected internal Queue<double> HistoryTime {
            get;
            private set;
        }

        /// <summary>
        /// stores ChangeRate values of the previous time steps
        /// </summary>
        protected internal Queue<double[]> HistoryChangeRate {
            get;
            private set;
        }

        /// <summary>
        /// stores computed change rate (Euler step)
        /// </summary>
        protected internal double[] CurrentChangeRate {
            get;
            set;
        }

        /// <summary>
        /// stores change rate after AdamsBashforth step (combination of
        /// currentChangRate with older steps)
        /// </summary>
        protected internal double[] CompleteChangeRate {
            get;
            set;
        }

        /// <summary>
        /// stores AdamsBashforth coefficients for variable time step size
        /// </summary>
        protected internal double[] ABCoefficients;

        /// <summary>
        /// The order of the scheme
        /// </summary>
        protected int order;

        /// <summary>
        /// Needed for the start-up, i.e. if no old time-steps are available
        /// </summary>
        public RungeKutta RungeKuttaScheme;

        /// <summary>
        /// Adams-Bashforth Constructor
        /// </summary>
        /// <param name="spatialOp">Spatial operator</param>
        /// <param name="Fieldsmap"></param>
        /// <param name="Parameters">
        /// optional parameter fields, can be null if
        /// <paramref name="spatialOp"/> contains no parameters; must match the
        /// parameter field list of <paramref name="spatialOp"/>, see
        /// <see cref="BoSSS.Foundation.SpatialOperator.ParameterVar"/>
        /// </param>
        /// <param name="order">Adams-Bashforth order</param>
        /// <param name="timeStepConstraints">
        /// optional list of time step constraints <see cref="TimeStepConstraint"/>
        /// </param>
        /// <param name="sgrd">
        /// optional restriction to computational domain
        /// </param>
        public AdamsBashforth(SpatialOperator spatialOp, CoordinateMapping Fieldsmap, CoordinateMapping Parameters, int order, IList<TimeStepConstraint> timeStepConstraints = null, SubGrid sgrd = null)
            : base(spatialOp, Fieldsmap, Parameters, SpatialOperator.SubGridBoundaryModes.InnerEdgeLTS, timeStepConstraints, sgrd) {
            if (order > 3 || order == 0) {
                throw new ArgumentException("Order not supported. Order must be between 1 and 3, but was " + order, "order");
            }
            this.order = order;

            HistoryChangeRate = new Queue<double[]>(order);

            RungeKuttaScheme = new RungeKutta(
                RungeKutta.SchemeFactory(RungeKutta.GetDefaultScheme(order)),
                spatialOp,
                Fieldsmap,
                Parameters,
                timeStepConstraints,
                sgrd);

            CurrentChangeRate = new double[Mapping.LocalLength];
            CompleteChangeRate = new double[Mapping.LocalLength];
            HistoryTime = new Queue<double>(order);
            ABCoefficients = new double[order];
        }

        /// <summary>
        /// another ctor.
        /// </summary>
        /// <param name="order"></param>
        /// <param name="spatialOp"></param>
        /// <param name="Fields"></param>
        public AdamsBashforth(int order, SpatialOperator spatialOp, params DGField[] Fields)
            : this(spatialOp, new CoordinateMapping(Fields), null, order) {
        }

        /// <summary>
        /// performs one time step
        /// </summary>
        /// <param name="dt">size of time step</param>
        public override double Perform(double dt) {
            using (new ilPSP.Tracing.FuncTrace()) {

                if (TimeStepConstraints != null) {
                    dt = CalculateTimeStep();
                }

                // new array object needed, because currentChangeRate from last time steps are stored in historyCR queue (in queues are only references of the array)
                CurrentChangeRate = new double[Mapping.LocalLength];
                CompleteChangeRate = new double[Mapping.LocalLength];

                // AB only possible if history with "order-1" entries exists
                // e.g: order=3 needs t_n (currentChangeRate), t_n-1 (history(2)), t_n-2 (history(1)) --> history has 2 entries 
                if (HistoryChangeRate.Count >= order - 1) {

                    ABCoefficients = ComputeCoefficients(dt, HistoryTime.ToArray());

                    ComputeChangeRate(CurrentChangeRate, m_Time, 0);
                    //y  <--  alpha*x + y
                    BLAS.daxpy(CompleteChangeRate.Length, ABCoefficients[0], CurrentChangeRate, 1, CompleteChangeRate, 1);
                    // calculate completeChangeRate
                    int i = 1;
                    foreach (double[] oldRate in HistoryChangeRate) {
                        BLAS.daxpy(CompleteChangeRate.Length, ABCoefficients[i], oldRate, 1, CompleteChangeRate, 1);
                        i++;
                    }

                    // perform time step
                    CurrentState.axpy<double[]>(CompleteChangeRate, -1);
                    m_Time += dt;


                    HistoryTime.Enqueue(m_Time);
                    HistoryTime.Dequeue();


                    HistoryChangeRate.Enqueue(CurrentChangeRate);
                    HistoryChangeRate.Dequeue();

                } else {

                    if (HistoryTime.Count == 0)
                        HistoryTime.Enqueue(RungeKuttaScheme.Time);

                    // Needed for the history
                    ComputeChangeRate(CurrentChangeRate, m_Time, 0);
                    // RungeKutta with order+1 starts the time stepping
                    RungeKuttaScheme.Perform(dt);

                    // Disposal RK after "order" of time steps
                    if (HistoryChangeRate.Count == order - 1) {
                        RungeKuttaScheme = null;
                    }

                    m_Time = RungeKuttaScheme.Time;
                    HistoryTime.Enqueue(RungeKuttaScheme.Time);
                    HistoryChangeRate.Enqueue(CurrentChangeRate);
                }
            }

            return dt;
        }


        /// <summary>
        /// Computes the Adams-Bashforth coefficients
        /// </summary>
        /// <param name="dt">Time step size</param>
        /// <param name="histTime">Time history</param>
        /// <remarks>
        /// Coefficients depend on the previous time steps, saved in historyTime
        /// </remarks>
        protected double[] ComputeCoefficients(double dt, double[] histTime) {
            double[] NewCoeff = new double[order];

            switch (order) {
                case 1:
                    NewCoeff[0] = dt;
                    return NewCoeff;
                case 2:
                    // DtHist[0] = t_n-1
                    // DtHist[1] = t_n
                    double deltaT_n = histTime[1] - histTime[0];
                    NewCoeff[0] = dt * (dt * 0.5 / deltaT_n + 1); //F_n
                    NewCoeff[1] = dt / deltaT_n * (-0.5 * dt); //F_n-1
                    return NewCoeff;
                case 3:
                    // DtHist[0] = t_n-2
                    // DtHist[1] = t_n-1
                    // DtHist[2] = t_n
                    double r1 = dt / (histTime[2] - histTime[1]);
                    double r0 = (histTime[2] - histTime[1]) / (histTime[1] - histTime[0]);
                    double rho = r0 / (r0 + 1);

                    NewCoeff[0] = dt / 6.0 * (2.0 * rho * r1 * r1 + 3.0 * r1 * (rho + 1.0) + 6.0); // F_n
                    NewCoeff[1] = dt / 6.0 * (r1 * r0 * rho * (2.0 * r1 + 3.0)); // F_n-2
                    NewCoeff[2] = dt / 6.0 * (-r1 * r0 * (3.0 / rho + 2 * r1)); // F_n-1
                    return NewCoeff;
                default:
                    throw new ArgumentException("Order not supported. Order must be between 2 or 3, but was " + order);
            }

        }
    }
}
