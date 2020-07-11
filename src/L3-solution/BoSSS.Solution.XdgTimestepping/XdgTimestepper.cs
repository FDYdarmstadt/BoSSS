using BoSSS.Foundation;
using BoSSS.Solution.Timestepping;
using ilPSP.LinSolvers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.XdgTimestepping {


    public enum TimeSteppingScheme {
        ExplicitEuler = 0,

        /// <summary>
        /// equivalent of BDF1, using the BDF-implementation
        /// </summary>
        ImplicitEuler = 1,

        CrankNicolson = 2000,

        BDF2 = 2,

        BDF3 = 3,

        BDF4 = 4,

        BDF5 = 5,

        BDF6 = 6,

        RK4 = 104,

        RK3 = 103,

        RK2 = 102,

        RK1 = 101,

        RK1u1 = 110,

        /// <summary>
        /// Implicit Euler using the implicit Runge-Kutta implementation
        /// </summary>
        RK_ImplicitEuler = 201,

        /// <summary>
        /// Crank Nicolson using the implicit Runge-Kutta implementation
        /// </summary>
        RK_CrankNic = 202,

        RK_IMEX3 = 203
    }


    /// <summary>
    /// Driver class which provides a simplified interface to <see cref="XdgBDFTimestepping"/> and <see cref="XdgRKTimestepping"/>
    /// </summary>
    public class XdgTimestepper {

        public TimeSteppingScheme Scheme {
            get;
            private set;
        }

        public XdgTimestepper(
            IEnumerable<DGField> Fields,
            IEnumerable<DGField> IterationResiduals,
            TimeSteppingScheme __Scheme,
            LevelSetTracker LsTrk,
            bool DelayInit,
            DelComputeOperatorMatrix _ComputeOperatorMatrix,
            DelComputeMassMatrix _ComputeMassMatrix,
            DelUpdateLevelset _UpdateLevelset,
            int BDForder,
            LevelSetHandling _LevelSetHandling,
            MassMatrixShapeandDependence _MassMatrixShapeandDependence,
            SpatialOperatorType _SpatialOperatorType,
            IDictionary<SpeciesId, IEnumerable<double>> _MassScale,
            MultigridOperator.ChangeOfBasisConfig[][] _MultigridOperatorConfig,
            AggregationGridData[] _MultigridSequence,
            SpeciesId[] _SpId,
            int _CutCellQuadOrder,
            double _AgglomerationThreshold, bool _useX, Control.NonLinearSolverConfig nonlinconfig,
            Control.LinearSolverConfig linearconfig) 
        {
            this.Scheme = __Scheme;


            RungeKuttaScheme rksch = null;
            int bdfOrder= -1000;
            if (this.Scheme == TimeSteppingScheme.CrankNicolson)
                bdfOrder = -1;
            else if (this.Scheme == TimeSteppingScheme.ExplicitEuler)
                bdfOrder = 0;
            else if (this.Scheme == TimeSteppingScheme.ImplicitEuler)
                bdfOrder = 1;
            else if (this.Scheme.ToString().StartsWith("BDF"))
                bdfOrder = Convert.ToInt32(this.Scheme.ToString().Substring(3));
            else if (this.Scheme == TimeSteppingScheme.RK1)
                rksch = RungeKuttaScheme.ExplicitEuler;
            else if (this.Scheme == TimeSteppingScheme.RK1u1)
                rksch = RungeKuttaScheme.ExplicitEuler2;
            else if (this.Scheme == TimeSteppingScheme.RK2)
                rksch = RungeKuttaScheme.Heun2;
            else if (this.Scheme == TimeSteppingScheme.RK3)
                rksch = RungeKuttaScheme.TVD3;
            else if (this.Scheme == TimeSteppingScheme.RK4)
                rksch = RungeKuttaScheme.RungeKutta1901;
            else if (this.Scheme == TimeSteppingScheme.RK_ImplicitEuler)
                rksch = RungeKuttaScheme.ImplicitEuler;
            else if (this.Scheme == TimeSteppingScheme.RK_CrankNic)
                rksch = RungeKuttaScheme.CrankNicolson;
            else if(this.Scheme == TimeSteppingScheme.RK_IMEX3)
                rksch = RungeKuttaScheme.IMEX3;
            else
                throw new NotImplementedException();
            

            if (bdfOrder > -1000) {
                m_BDF_Timestepper = new XdgBDFTimestepping(Fields, IterationResiduals, 
                    LsTrk, true,
                    DelComputeOperatorMatrix, null, DelUpdateLevelset,
                    bdfOrder,
                    lsh,
                    MassMatrixShapeandDependence.IsTimeDependent,
                    SpatialOperatorType.LinearTimeDependent,
                    MassScale,
                    null, base.MultigridSequence,
                    this.LsTrk.SpeciesIdS.ToArray(), quadOrder,
                    this.Control.AgglomerationThreshold, false,
                    this.Control.NonLinearSolver,
                    this.Control.LinearSolver);
            } else {
                m_RK_Timestepper = new XdgRKTimestepping(Fields, IterationResiduals,
                    LsTrk,
                    DelComputeOperatorMatrix, DelUpdateLevelset,
                    rksch,
                    lsh,
                    MassMatrixShapeandDependence.IsTimeDependent,
                    SpatialOperatorType.LinearTimeDependent,
                    MassScale,
                    null, base.MultigridSequence,
                    this.LsTrk.SpeciesIdS.ToArray(), quadOrder,
                    this.Control.AgglomerationThreshold, false,
                    this.Control.NonLinearSolver,
                    this.Control.LinearSolver);
            }

        }


        void DelComputeOperatorMatrix(BlockMsrMatrix OpMtx, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double time) {

        }




        XdgBDFTimestepping m_BDF_Timestepper;



        XdgRKTimestepping m_RK_Timestepper;


        public void Solve(double phystime, double dt) {
            if((m_BDF_Timestepper == null) == (m_RK_Timestepper == null))
                throw new ApplicationException();

            if(m_BDF_Timestepper != null) {
                m_BDF_Timestepper.Solve(phystime, dt);
            } else {
                m_RK_Timestepper.Solve(phystime, dt);
            }

        }
    }
}
