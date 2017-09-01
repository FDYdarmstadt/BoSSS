using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Platform;
using ilPSP.LinSolvers;
using BoSSS.Foundation.Grid;
using ilPSP.Utils;

namespace BoSSS.Solution {
    public abstract class OperatorSplitting {
        /// <summary>
        /// The matrix for the operator on  the narrow band
        /// </summary>
        public MsrMatrix SubgridOperatorMatr;
        /// <summary>
        /// The affine part of the operator matrix on the narrow band
        /// </summary>
        public double[] SubgridAffine;
        /// <summary>
        /// The variable which the operator splitting is applied to
        /// </summary>
        public SubgridCoordinateMapping sourcevariable;
        /// <summary>
        /// The context of the application
        /// </summary>
        internal Context m_Context;

        public Context context {
            get { return m_Context; }
        }
      

        public SubGrid m_Subgrid;
        public  double[] affine;
        public MsrMatrix opmatr;
        public OperatorSplitting(Context context, SubGrid subgrid
            , SubgridCoordinateMapping v){
            m_Context=context;

            sourcevariable=v;
            sourcevariable.Compress();
            Partition part = new Partition(v.NUpdate, 1);
            affine=new double[v.NUpdate];
            //opmatr= new MsrMatrix(part,v.MaxTotalNoOfCoordinatesPerCell * (int)m_Context.GridDat.GlobalNoOfCells);
            opmatr = new MsrMatrix(part);
            m_Subgrid = subgrid;
            
        }

       
        /// <summary>
        /// Returns the source which contains Au_{n-1}+b for the Operator and the field provided in the constructor
        /// </summary>
        /// <param name="k">results of Ay+b</param>
        /// <param name="dt">optional scaling by time step size</param>
        public void ExplicitEulerSource(SubgridCoordinateMapping u, out double[] convectiveSource, double dt) {
            
            convectiveSource = new double[u.subgridCoordinates.Length];
            u.Compress();
            SubgridOperatorMatr.SpMVpara<double[], double[]>(dt, u.subgridCoordinates, 1.0, convectiveSource);
            BLAS.daxpy(SubgridAffine.Length, dt, SubgridAffine, 1, convectiveSource, 1);

        }
        /// <summary>
        /// Returns the source which contains Au_{n-1}+b for the Operator and the field provided in the constructor
        /// </summary>
        /// <param name="k">results of Ay+b</param>
        /// <param name="dt">optional scaling by time step size</param>
        public void ExplicitEulerSource(double[] subgridCoordinates, out double[] convectiveSource, double dt) {

            convectiveSource = new double[subgridCoordinates.Length];
           
            SubgridOperatorMatr.SpMVpara<double[], double[]>(dt, subgridCoordinates, 1.0, convectiveSource);
            BLAS.daxpy(SubgridAffine.Length, dt, SubgridAffine, 1, convectiveSource, 1);

        }
        public abstract void CreateMatrix(MsrMatrix matr, double[] affine);
    }


}
