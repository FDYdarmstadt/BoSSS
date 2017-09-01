using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;

using ilPSP.LinSolvers;
using BoSSS.Platform;

namespace BoSSS.Solution{
    public class ExplicitConvection {
        SpatialOperator convectionop;
        string name;
        Context m_Context;
        SinglePhaseField[] creepingFlow;
        MsrMatrix SubgridOperatorMatr;
        double[] SubgridAffine;
        public ExplicitConvection(Context context, string variablename,SubGrid subgrid,int basisdegree
            , SubgridCoordinateMapping v){
            m_Context=context;
            name=variablename;
            convectionop=new SpatialOperator(1,3,1,name,"u","v","w");
            convectionop.EquationComponents["Density"].Add(new SurfaceConvectionUpwinding(new string[]{"u","v","w"}, subgrid, new string[]{name},0));
            convectionop.Commit();
            CreepingFlowFactory fact = new CreepingFlowFactory(m_Context, basisdegree);
            fact.CreateFlowFields(out creepingFlow);
            Partition part = new Partition(v.NUpdate);
            double[] affine=new double[v.NUpdate];
            MsrMatrix opmatr= new MsrMatrix(part,v.MaxTotalNoOfCoordinatesPerCell * (int)m_Context.GridDat.GlobalNoOfCells);

            convectionop.ComputeMatrixEx(m_Context, v, new Field[] { creepingFlow[0], creepingFlow[1], creepingFlow[2] }, v, opmatr, affine, false, subgrid);

            v.SetupSubmatrix(affine, opmatr, out SubgridAffine, out SubgridOperatorMatr);
        }
        public void ExplicitEulerSource(SubgridCoordinateMapping u, out double[] convectiveSource, double dt){
            convectiveSource= new double[u.subgridCoordinates.Length];
            Evaluate(u, convectiveSource, dt);
        }

        /// <summary>
        /// Evaluation of the operator on the subgrid by matrix-vector product. There might be a more efficient method....
        /// </summary>
        /// <param name="k">results of Ay+b</param>
        /// <param name="dt">optional scaling by time step size</param>
        protected void Evaluate(SubgridCoordinateMapping u,double[] k, double dt) {

            SubgridOperatorMatr.SpMVpara<double[], double[]>(dt, u.subgridCoordinates, 1.0, k);
            BLAS.daxpy(SubgridAffine.Length, dt, SubgridAffine, 1, k, 1);

        }
    
    }
}
