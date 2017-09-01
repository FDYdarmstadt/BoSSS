/*
 *
 * Copyright (c) 2010, Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsmechanik
 *
 * This file is part of the BoSSS software. 
 * The software (source code or binaries compiled from the source code) may not
 * be copied, compiled ore executed, partly or as a whole, without an explicit 
 * written permission from the Fachgebiet fuer Stroemungsmechanik, TU Darmstadt.
 *
 */
using System;
using System.Collections.Generic;
using BoSSS.Foundation;
using BoSSS.Platform;
using BoSSS.Foundation.Grid;
using ilPSP.LinSolvers;
using ilPSP.Utils;

//tst

namespace BoSSS.Solution.Timestepping{

    /// <summary>
    /// Runge-Kutta time-stepping for narrow band applications;
    /// </summary>
    public class RungeKuttaSubgrid : ExplicitEulerSubgrid{

        private RungeKutta.RungeKuttaScheme m_RKscheme;

        /// <summary>
        /// the Runge-Kutta scheme
        /// </summary>
        public RungeKutta.RungeKuttaScheme Scheme {
            get {
                return m_RKscheme;
            }
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="ctx">Context</param>
        /// <param name="subgridMapping">Coordinate Mapping that may be restricted to the subgrid</param>
        /// <param name="operatorMatrix">Matrix associated with the differential operator</param>
        /// <param name="affine">Affine part of the operator matrix</param>
        /// <param name="scheme">Runge Kutta scheme, available in <see cref="RungeKutta"/></param>
        /// <remarks>
        /// This constructor does not support parameter fields; 
        /// </remarks>
        public RungeKuttaSubgrid(Context ctx,  SubgridCoordinateMapping subgridMapping, MsrMatrix operatorMatrix, double[] affine, RungeKutta.RungeKuttaScheme scheme)
            : base(ctx,  subgridMapping, null, operatorMatrix, affine) {
                m_RKscheme = scheme;
        }

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="ctx">Context</param>
        /// <param name="subgridMapping">Coordinate Mapping that may be restricted to the subgrid</param>
        /// <param name="operatorMatrix">Matrix associated with the differential operator</param>
        /// <param name="affine">Affine part of the operator matrix</param>
        /// <param name="scheme">Runge-Kutta scheme, available in <see cref="RungeKutta"/></param>
        /// <param name="Parameters">
        /// optional parameter fields, can be null if the spatial operator contains no parameters;
        /// must match the parameter field list of the spatial operator, see <see cref="BoSSS.Foundation.SpatialOperator.ParameterVar"/>
        /// </param>
        public RungeKuttaSubgrid(Context ctx, SubgridCoordinateMapping subgridMapping, CoordinateMapping Parameters, 
            MsrMatrix operatorMatrix, double[] affine, RungeKutta.RungeKuttaScheme scheme):base(ctx,subgridMapping,Parameters,operatorMatrix,affine){
                m_RKscheme = scheme; 
        }

        /// <summary>
        /// performs one timestep
        /// </summary>
        /// <param name="dt">size of timestep</param>
        public override void Perform(double dt) {

            using (var tr = new ilPSP.Tracing.FuncTrace()) {


                double[][] k = new double[m_RKscheme.Stages][];
                for (int i = 0; i < m_RKscheme.Stages; i++)
                    k[i] = new double[m_SubgridMapping.subgridCoordinates.Length];

                double[] y0 = new double[m_SubgridMapping.subgridCoordinates.Length];
                m_DGCoordinates.CopyTo(y0, 0);

                // logging
                tr.Info("Runge-Kutta Scheme with " + m_RKscheme.Stages + " stages.");
                double time0 = m_Time;
                tr.Info("time = " + time0 + ", dt = " + dt);

                // berechne k[0]
                Evaluate(k[0], m_Time);// invokes MPI communication
               
                for (int s = 1; s < m_RKscheme.Stages; s++) {
                    PerformStage(y0, s, k, dt);

                    m_Time = time0 + m_RKscheme.c[s] * dt;
                    //ComputeChangeRate(k[s], m_Time, m_RKscheme.c[s] * dt);// invokes MPI communication
                    Evaluate(k[s], m_Time);
                }

                // next timestep
               
                y0.CopyTo(m_DGCoordinates, 0);
                Array.Clear(y0, 0, y0.Length);
                for (int s = 0; s < m_RKscheme.Stages; s++)
                    BLAS.daxpy(y0.Length, -m_RKscheme.b[s] * dt, k[s], 1, y0, 1);
                //m_DGCoordinates.axpy<double[]>(y0, 1.0);
                BLAS.daxpy(m_DGCoordinates.Length, 1.0, y0, 1, m_DGCoordinates, 1);
                m_Time = time0 + dt;

            }
        }
        /// <summary>
        /// Performs a single stage of a timestep.
        /// </summary>
        /// <param name="y0">Initial DG coordinates</param>
        /// <param name="s">The current stage</param>
        /// <param name="k">The current change rate</param>
        /// <param name="dt">The size of the timestep</param>
        protected virtual void PerformStage(double[] y0, int s, double[][] k, double dt) {
          
            Array.Clear(m_DGCoordinates, 0, DGCoordinates.Length);
            //m_DGCoordinates.axpy<double[]>(y0, 1.0);
            BLAS.daxpy(m_DGCoordinates.Length, 1.0, y0, 1, m_DGCoordinates, 1);
            for (int r = 0; r < s; r++) {
                if (m_RKscheme.a[s, r] != 0.0) {
                    //m_DGCoordinates.axpy<double[]>(k[r], -dt * m_Scheme.a[s, r]);
                    BLAS.daxpy(m_DGCoordinates.Length, -dt*m_RKscheme.a[s,r], k[r], 1, m_DGCoordinates, 1);
                }
            }
        }

    }
}
