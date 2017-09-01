/*
 *
 * Copyright (c) 2010, Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)
 *
 * This file is part of the BoSSS software. 
 * The software (source code or binaries compiled from the source code) may not
 * be copied, compiled ore executed, partly or as a whole, without an explicit 
 * written permission from the Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics), TU Darmstadt.
 *
 */
using System;
using System.Collections.Generic;
using System.IO;
using System.Globalization;

using BoSSS.Foundation;
using BoSSS.Foundation.Comm;
using BoSSS.Foundation.Grid;
using System.Collections;
using ilPSP.LinSolvers;
using BoSSS.Platform;
using ilPSP.Utils;

namespace BoSSS.Solution.Timestepping {

#pragma warning disable 1572
#pragma warning disable 1573
#pragma warning disable 1574
#pragma warning disable 1591
#pragma warning disable 1734

    /// <summary>
    /// explicit Euler timestepping with auxiliary condition (that is a non-time dependent equation);
    /// this class also serves as a baseclass for Runge-Kutta (and maybe other) timesteppers (<see cref="RungeKutta"/>).
    /// </summary>
    public class ExplicitEulerSubgrid {

        /// <summary>
        /// the BoSSS kernel
        /// </summary>
        protected Context m_Context;

        /// <summary>
        /// Mapping of DG fields that are affected by this timestepper restricted to a subgrid.
        /// All entries outside of the subgrid domain are set to 0.
        /// </summary>
        public SubgridCoordinateMapping SubgridMapping {
            get {
                m_SubgridMapping.subgridCoordinates = DGCoordinates;
                m_SubgridMapping.Decompress();
                return m_SubgridMapping; }
        }

        /// <summary>
        /// Mapping of fields for which evolution-equations are defined, with restriction to a narrow band.
        /// </summary>
        protected SubgridCoordinateMapping m_SubgridMapping;

        /// <summary>
        /// mapping of parameter fields
        /// </summary>
        protected CoordinateMapping m_Parameters;

        /// <summary>
        /// <see cref="DGCoordinates"/>
        /// </summary>
        internal double[] m_DGCoordinates;
        /// <summary>
        /// Restriction of the operator matrix to a subgrid.
        /// </summary>
        MsrMatrix subgridMatrix;
        double[] subgridAffine;

        /// <summary>
        /// The DG coordinates of the <see cref="BoSSS.Foundation.Field"/>'s in the original coordinate mapping
        /// (see <see cref="BoSSS.Foundation.CoordinateMapping.Fields"/>);
        /// </summary>
        public double[] DGCoordinates {
            get { return m_DGCoordinates; }
        }

        /// <summary>
        /// Constructor for an explicit Euler scheme operating on subgrids(without parameters).
        /// </summary>
        /// <param name="ctx"></param>
        /// <param name="subgridMapping">Coordinate Mapping on the subgrid</param>
        /// <param name="operatorMatrix">Matrix of the differential operator</param>
        /// <param name="affine">Affine part of the operator matrix</param>
        public ExplicitEulerSubgrid(Context ctx,  SubgridCoordinateMapping subgridMapping,MsrMatrix operatorMatrix, double[] affine)
            : this(ctx,  subgridMapping, null, operatorMatrix, affine) { }

         /// <summary>
         ///  Constructor for an explicit Euler scheme operating on subgrids supporting parameters.
         /// </summary>
         /// <param name="ctx"></param>
         /// <param name="subgridMapping">Coordinate Mapping on the subgrid</param>
         /// <param name="Parameters">Optional parameters which have to match the matrix dimensions</param>
         /// <param name="operatorMatrix">Matrix of the differential operator</param>
         /// <param name="affine">Affine part of the operator matrix</param>
        public ExplicitEulerSubgrid(Context ctx, SubgridCoordinateMapping subgridMapping, CoordinateMapping Parameters, MsrMatrix operatorMatrix, double[] affine) {

            using (new ilPSP.Tracing.FuncTrace()) {
                m_Context = ctx;
                m_SubgridMapping = subgridMapping;
                m_DGCoordinates = subgridMapping.subgridCoordinates;
                m_Parameters = Parameters;
                IList<Field> ParameterFields = (m_Parameters == null) ? (new Field[0]) : m_Parameters.Fields;
              
                m_SubgridMapping.SetupSubmatrix(affine, operatorMatrix, out subgridAffine, out subgridMatrix);
            }
        }
        
        public delegate void ChangeRateCallback(double AbsTime, double RelTime);


        virtual protected void ComputeChangeRate(double[] k, double time) {
           Evaluate(k, time);
        }

   
        //Include in Euler Time Stepping......
        /// <summary>
        /// Evaluation of the operator on the subgrid by matrix-vector product. There might be a more efficient method....
        /// </summary>
        /// <param name="k">results of Ay+b</param>
        /// <param name="dt">optional scaling by time step size</param>
        protected void Evaluate(double[] k,double dt) {

            subgridMatrix.SpMVpara<double[], double[]>(dt, m_SubgridMapping.subgridCoordinates, 1.0, k);
            //shift OperatorMatrix*subgridCoordinates+subgridAffine into mappingCopy
            BLAS.daxpy(subgridAffine.Length, 1.0, subgridAffine, 1, k, 1);
        
        }
    
        /// <summary>
        /// performs one Explicit Euler timestep
        /// </summary>
        /// <param name="dt">size of timestep</param>
   

        public virtual void Perform(double dt) {

            using (new ilPSP.Tracing.FuncTrace()) {

                double[] k = new double[m_SubgridMapping.subgridCoordinates.Length];
                Evaluate(k,dt);
                BLAS.daxpy(DGCoordinates.Length, -1.0, k, 1, DGCoordinates, 1);
            }
        }


        #region ITimeStepper Members

        /// <summary>
        /// <see cref="ITimeStepper.Time"/>
        /// </summary>
        public double Time {
            get { return m_Time; }
        }

        /// <summary>
        /// physical time 
        /// </summary>
        protected double m_Time = 0.0;

        /// <summary>
        /// <see cref="ITimeStepper.ResetTime"/>
        /// </summary>
        /// <param name="NewTime"></param>
        public void ResetTime(double NewTime) {
            m_Time = NewTime;
        }

        #endregion
    }

#pragma warning restore 1572
#pragma warning restore 1573
#pragma warning restore 1591
#pragma warning restore 1734
#pragma warning restore 1574


}
