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
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.Quadrature;
using BoSSS.Solution;
using BoSSS.Solution.Utils;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using ilPSP.Utils;
using System.Diagnostics;
using MPI.Wrappers;
using BoSSS.Platform;
using ilPSP;
using System.Linq;
using BoSSS.Foundation.SpecFEM;
using BoSSS.Solution.Queries;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Solution.AdvancedSolvers;
using ilPSP.Connectors.Matlab;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Solution.XdgTimestepping;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution.NSECommon;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Control;

namespace BoSSS.Solution.LevelSetTools.PhasefieldLevelSet
{
    /// <summary>
    /// Base class for a level set based on the Cahn-Hilliard Phasefield equation
    /// Insert Latex style equations here
    /// </summary>
    public partial class Phasefield : Application
    {
#pragma warning disable 649
        /// <summary>
        /// concentration aka the level set
        /// </summary>
        protected SinglePhaseField phi;

        /// <summary>
        /// concentration (linearization point)
        /// </summary>
        protected SinglePhaseField phi0;

        /// <summary>
        /// chemical potential
        /// </summary>
        protected SinglePhaseField mu;

        /// <summary>
        /// residual of 'c(phi)'-equation
        /// </summary>
        protected SinglePhaseField phi_Resi;

        /// <summary>
        /// residual of 'mu'-equation
        /// </summary>
        protected SinglePhaseField mu_Resi;


        /// <summary>
        /// Transport velocity
        /// </summary>
        protected VectorField<SinglePhaseField> Velocity;
#pragma warning restore 649

        /// <summary>
        /// LevelSet
        /// </summary>
        LevelSet LevSet;

        /// <summary>
        /// DG LevelSet
        /// </summary>
        protected SinglePhaseField DGLevSet;

        /// <summary>
        /// LevelSetTracker
        /// </summary>
        LevelSetTracker LsTrk;

        /// <summary>
        /// DummyLevelSet, as the Phasefield is in DG not XDG
        /// </summary>
        LevelSet DummyLevSet;

        /// <summary>
        /// DummyLevelSetTracker
        /// </summary>
        LevelSetTracker DummyLsTrk;

        /// <summary>
        /// Cahn number determines interface thickness
        /// </summary>
        double Cahn;

        /// <summary>
        /// (Bulk) Diffusion coefficient
        /// </summary>
        double Diff;

        /// <summary>
        /// Peclet number, ratio of convective to diffusive timescale has to be set so that diffusive timescale is smaller
        /// </summary>
        double Peclet;

        /// <summary>
        /// Control for Bulk vs Surface Diffusion
        /// </summary>
        double Lambda;

        /// <summary>
        /// Cahn-Hilliard spatial Operator
        /// </summary>
        SpatialOperator CHOp;

        /// <summary>
        /// Boundary Condition map for Cahn Hilliard
        /// </summary>
        BoundaryCondMap<BoundaryType> m_bcMap;

        /// <summary>
        /// Level Set Timestepper
        /// </summary>
        Solution.XdgTimestepping.XdgBDFTimestepping m_Timestepper;

        /// <summary>
        /// GridData
        /// </summary>
        IGridData GridData;

        /// <summary>
        /// MultigridSequence
        /// </summary>
        AggregationGridData[] mgSeq;

        /// <summary>
        /// Control of base solver that calls this level set method
        /// </summary>
        AppControl ParentControl;

        /// <summary>
        /// Phasefield instantiation, constructor
        /// </summary>
        public Phasefield(LevelSet _LevSet, SinglePhaseField _DGLevSet,  LevelSetTracker _LsTrk, VectorField<SinglePhaseField> _Velocity, IGridData _GridData, AppControl _control, AggregationGridData[] _mgSeq, double _peclet = 1e3)
        {
            LevSet = _LevSet;
            LsTrk = _LsTrk;
            DGLevSet = _DGLevSet;
            Velocity = _Velocity;
            GridData = _GridData;
            ParentControl = _control;
            mgSeq = _mgSeq;
            Peclet = _peclet;
        }

        /// <summary>
        /// Updating Phasefield after changes in Grid and/or Basis
        /// </summary>
        public void UpdateFields(LevelSet _LevSet, SinglePhaseField _DGLevSet, LevelSetTracker _LsTrk, VectorField<SinglePhaseField> _Velocity, IGridData _GridData, AppControl _control, AggregationGridData[] _mgSeq)
        {
            // check if signature of external and local Grid changed
            if (this.GridData != _GridData)
            {
                Console.WriteLine("Grid changed, adapting Phasefield");
                LevSet = _LevSet;
                LsTrk = _LsTrk;
                DGLevSet = _DGLevSet;
                Velocity = _Velocity;
                GridData = _GridData;
                ParentControl = _control;
                mgSeq = _mgSeq;

                CreateFields();
                CreateEquationsAndSolvers(null);
            }
        }

        /// <summary>
        /// Create Fields with same basis as DG Level Set
        /// </summary>
        protected override void CreateFields()
        {
            phi0 = new SinglePhaseField(DGLevSet.Basis, "phi0");
            phi = new SinglePhaseField(DGLevSet.Basis, "phi");
            phi.Acc(1.0, DGLevSet);
            mu = new SinglePhaseField(DGLevSet.Basis, "mu");
            phi_Resi = new SinglePhaseField(DGLevSet.Basis, "phi_Resi");
            mu_Resi = new SinglePhaseField(DGLevSet.Basis, "mu_Resi");

            DummyLevSet = new LevelSet(new Basis(this.GridData, 1), "Levset");
            DummyLevSet.AccConstant(-1);
            this.DummyLsTrk = new LevelSetTracker((GridData)(this.GridData), XQuadFactoryHelper.MomentFittingVariants.Saye, 1, new string[] { "A", "B" }, DummyLevSet);
            this.DummyLsTrk.UpdateTracker();
        }

        /// <summary>
        /// Includes assembly of the matrix.
        /// </summary>
        /// <param name="L"></param>
        protected override void CreateEquationsAndSolvers(GridUpdateDataVaultBase L)
        {
            using (FuncTrace tr = new FuncTrace())
            {

                // create operator
                // ===============
                {

                    int D = this.GridData.SpatialDimension;
                    double _D = D;
                    double penalty_base = (phi.Basis.Degree + 1) * (phi.Basis.Degree + _D) / _D;

                    // Get this from where?
                    double penalty_factor = 2.6 * penalty_base;

                    //BoundaryCondMap<BoundaryType> PoissonBcMap = new BoundaryCondMap<BoundaryType>(this.GridData, this.Control.BoundaryValues, "T");

                    CHOp = new SpatialOperator(
                        new string[] { "phi", "mu" },
                        VariableNames.VelocityVector(D).Cat("phi0"),
                        new string[] { "Res_phi", "Res_mu" },
                        QuadOrderFunc.NonLinear(3)
                        );

                    MultidimensionalArray LengthScales;
                    if (this.GridData is GridData)
                    {
                        LengthScales = ((GridData)this.GridData).Cells.cj;
                    }
                    else if (this.GridData is AggregationGridData)
                    {
                        LengthScales = ((AggregationGridData)this.GridData).AncestorGrid.Cells.cj;
                    }
                    else
                    {
                        throw new NotImplementedException();
                    }

                    m_bcMap = new BoundaryCondMap<BoundaryType>(this.GridData, BoundaryTranslator(ParentControl.BoundaryValues), "phi");

                    SetCHCoefficents(out Cahn, out Diff, out Lambda);

                    CHOp.EquationComponents["Res_phi"].Add(
                    new phi_Diffusion(D, penalty_factor, LengthScales, Diff, Lambda, m_bcMap)
                    );


                    CHOp.EquationComponents["Res_phi"].Add(
                    new phi_Flux(D, m_bcMap)
                    );


                    CHOp.EquationComponents["Res_mu"].Add(
                        new mu_Diffusion(D, penalty_factor, LengthScales, Cahn, m_bcMap)
                        );

                    CHOp.EquationComponents["Res_mu"].Add(
                        new mu_Source(true, Cahn)
                        );

                    CHOp.Commit();
                }


                // create solver
                // =============

                {              

                    m_Timestepper = new XdgBDFTimestepping(
                        new DGField[] { phi, mu },
                        new DGField[] { phi_Resi, mu_Resi },
                        this.DummyLsTrk,
                        false,
                        DelComputeOperatorMatrix, null, null,
                        1, // BDF order
                        LevelSetHandling.None, MassMatrixShapeandDependence.IsTimeDependent, SpatialOperatorType.Nonlinear,
                        this.MassScale, this.MgConfig, mgSeq, new[] { this.DummyLsTrk.GetSpeciesId("A") }, phi.Basis.Degree * 2, 0.0, false,
                        GetNonLinearSolver(), GetLinearSolver()
                        );


                }
            }
        }

        /// <summary>
        /// Generate linear solver config for Cahn-Hilliard Level Set
        /// </summary>
        /// <returns></returns>
        private LinearSolverConfig GetLinearSolver()
        {
            LinearSolverConfig LinConfig = new LinearSolverConfig();
            return LinConfig;
        }

        /// <summary>
        /// Generate nonlinear solver config for Cahn-Hilliard Level Set
        /// </summary>
        /// <returns></returns>
        private NonLinearSolverConfig GetNonLinearSolver()
        {
            NonLinearSolverConfig NonLinConfig = new NonLinearSolverConfig();
            return NonLinConfig;
        }

        MultigridOperator.ChangeOfBasisConfig[][] MgConfig
        {
            get
            {
                int p = this.phi.Basis.Degree;
                int NoOfLevels = mgSeq.Length;
                var config = new MultigridOperator.ChangeOfBasisConfig[NoOfLevels][];

                for (int iLevel = 0; iLevel < NoOfLevels; iLevel++)
                {
                    config[iLevel] = new MultigridOperator.ChangeOfBasisConfig[] {
                        new MultigridOperator.ChangeOfBasisConfig() {
                            VarIndex = new int[] {0, 1},
                            mode = MultigridOperator.Mode.SymPart_DiagBlockEquilib,
                            DegreeS = new int[] { Math.Max(1, p - iLevel), Math.Max(1, p - iLevel) }
                        }
                    };

                }

                return config;
            }

        }

        /// <summary>
        /// Computation of operator matrix
        /// </summary>
        void DelComputeOperatorMatrix(BlockMsrMatrix OpMtx, double[] OpAffine, UnsetteledCoordinateMapping Mapping, DGField[] CurrentState, Dictionary<SpeciesId, MultidimensionalArray> AgglomeratedCellLengthScales, double time)
        {
            SinglePhaseField Current_phi = (SinglePhaseField)(CurrentState[0]);
            SinglePhaseField Current_mu = (SinglePhaseField)(CurrentState[1]);

            phi0.Clear();
            phi0.Acc(1.0, Current_phi);
            
            OpMtx.Clear();
            OpAffine.ClearEntries();
            var mb = CHOp.GetMatrixBuilder(Mapping, this.Velocity.ToArray().Cat(phi0), Mapping);
            mb.ComputeMatrix(OpMtx, OpAffine);
        }


        /// <summary>
        /// Block scaling of the mass matrix: for each species $\frakS$, a vector $(\rho_\frakS, \ldots, \rho_frakS, 0 )$.
        /// </summary>
        protected IDictionary<SpeciesId, IEnumerable<double>> MassScale
        {
            get
            {

                Dictionary<SpeciesId, IEnumerable<double>> R = new Dictionary<SpeciesId, IEnumerable<double>>();
                R.Add(this.DummyLsTrk.GetSpeciesId("A"), new double[] { 1.0, 0.0 });

                return R;
            }
        }

        /// <summary>
        /// default plotting
        /// </summary>
        protected override void PlotCurrentState(double phystime, TimestepNumber timestepNo, int superSampling = 0)
        {
            string caseStr = "";
            
            DGField[] Fields = new DGField[0];
            Fields = Fields.Cat(this.phi, this.mu, this.Velocity);
            BoSSS.Solution.Tecplot.Tecplot.PlotFields(Fields, "Phasefield-" + timestepNo + caseStr, phystime, superSampling);
        }

    }

    /// <summary>
    /// Interior Penalty Flux, with Dirichlet boundary conditions for variable 'phi'
    /// </summary>
    class phi_Diffusion : BoSSS.Solution.NSECommon.SIPLaplace
    {

        public phi_Diffusion(int D, double penalty_const, MultidimensionalArray cj, double __diff, double __lambda, BoundaryCondMap<BoundaryType> __boundaryCondMap)
            : base(penalty_const, cj, "mu") // note: in the equation for 'phi', we have the Laplacian of 'mu'
        {
            m_D = D;
            m_diff = __diff;
            m_boundaryCondMap = __boundaryCondMap;
            m_lambda = __lambda;
        }

        double m_diff;
        double m_lambda;
        BoundaryCondMap<BoundaryType> m_boundaryCondMap;

        int m_D;
        public override IList<string> ParameterOrdering => new[] { "phi0" }.Cat(VariableNames.VelocityVector(m_D));

        protected override double g_Diri(ref CommonParamsBnd inp)
        {
            double UxN = 0;
            for (int d = 0; d < m_D; d++)
            {
                UxN += (inp.Parameters_IN[d + 1]) * inp.Normal[d];
            }

            double v;
            if (UxN >= 0)
            {
                v = 1.0;
            }
            else
            {
                v = 0.0;
            }

            return v;
        }

        protected override double g_Neum(ref CommonParamsBnd inp)
        {
            return 0.0;
        }

        protected override bool IsDirichlet(ref CommonParamsBnd inp)
        {
            BoundaryType edgeType = m_boundaryCondMap.EdgeTag2Type[inp.EdgeTag];
            switch (edgeType)
            {
                case BoundaryType.Wall:
                case BoundaryType.Slip:
                case BoundaryType.SlipSymmetry:
                    // a Dicirchlet b.c. for 'c' mean a Neumann b.c. for 'phi'
                    return false;

                case BoundaryType.Inlet:
                case BoundaryType.Outlet:
                case BoundaryType.Outflow:
                case BoundaryType.Pressure_Dirichlet:
                default:
                    return true;
            }
        }

        public override double Nu(double[] x, double[] p, int jCell)
        {
            if (m_lambda == 0.0)
            {
                return -m_diff;
            }
            else if (m_lambda > 0.0 && m_lambda <= 1.0)
            {
                return -m_diff * Math.Max(1 - m_lambda * Math.Pow(p[0], 2.0),0.0);
            }
            else
            {
                throw new ArgumentOutOfRangeException();
            }
        }

        /// <summary>
        /// Integrand on boundary mesh edges of the SIP
        /// </summary>
        public override double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA)
        {
            double Acc = 0.0;

            double pnlty = 2 * this.GetPenalty(inp.jCellIn, -1);//, inp.GridDat.Cells.cj);
            double nuA = this.Nu(inp.X, inp.Parameters_IN, inp.jCellIn);

            if (this.IsDirichlet(ref inp))
            {
                // inhom. Dirichlet b.c.
                // +++++++++++++++++++++

                double g_D = this.g_Diri(ref inp);

                // inflow / outflow
                if (g_D == 0)
                {
                    for (int d = 0; d < inp.D; d++)
                    {
                        double nd = inp.Normal[d];
                        Acc += (nuA * _Grad_uA[0, d]) * (_vA) * nd;        // consistency
                        Acc += (nuA * _Grad_vA[d]) * (_uA[0] - g_D) * nd;  // symmetry
                    }
                    Acc *= this.m_alpha;

                    Acc -= nuA * (_uA[0] - g_D) * (_vA - 0) * pnlty; // penalty
                }
                else
                {
                    for (int d = 0; d < inp.D; d++)
                    {
                        double nd = inp.Normal[d];
                        //Acc += (nuA * _Grad_uA[0, d]) * (_vA) * nd;        // consistency 
                        //Acc += (nuA * 0.0) * (_vA) * nd;        // consistency, homogenous Neumann
                    }
                    Acc *= this.m_alpha;

                }
            }
            else
            {

                double g_N = this.g_Neum(ref inp);

                Acc += nuA * g_N * _vA * this.m_alpha;
            }
            return Acc;
        }
    }

    /// <summary>
    /// Transport flux for Cahn-Hilliard
    /// </summary>
    class phi_Flux : IVolumeForm, IEdgeForm
    {
        public phi_Flux(int D, BoundaryCondMap<BoundaryType> __boundaryCondMap)
        {
            m_D = D;
            m_boundaryCondMap = __boundaryCondMap;
            m_bndFunc = m_boundaryCondMap.bndFunction["phi"];
        }

        int m_D;
        BoundaryCondMap<BoundaryType> m_boundaryCondMap;
        Func<double[], double, double>[] m_bndFunc;

        public TermActivationFlags VolTerms => TermActivationFlags.UxGradV;

        public IList<string> ArgumentOrdering => new string[] { "phi" };

        public IList<string> ParameterOrdering => VariableNames.VelocityVector(m_D);

        public TermActivationFlags BoundaryEdgeTerms => InnerEdgeTerms | TermActivationFlags.V;

        public TermActivationFlags InnerEdgeTerms => TermActivationFlags.UxV;

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uIN, double[,] _Grad_uIN, double _vIN, double[] _Grad_vIN)
        {
            // expand for treatment of input functions, for now hardcode to -1.0
            double UxN = 0;
            for (int d = 0; d < m_D; d++)
            {
                UxN += (inp.Parameters_IN[d]) * inp.Normal[d];
            }

            double phi;
            if (UxN >= 0)
            {
                phi = _uIN[0];
            }
            else
            {
                phi = -1.0;//m_bndFunc[inp.EdgeTag](inp.X, inp.time);
            }

            return phi * UxN * _vIN;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT)
        {
            double UxN = 0;
            for (int d = 0; d < m_D; d++)
            {
                UxN += 0.5 * (inp.Parameters_IN[d] + inp.Parameters_OUT[d]) * inp.Normal[d];
            }

            double phi;
            if (UxN >= 0)
            {
                phi = _uIN[0];
            }
            else
            {
                phi = _uOUT[0];
            }

            return phi * UxN * (_vIN - _vOUT);
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV)
        {
            double acc = 0;
            double phi = U[0];
            for (int d = 0; d < m_D; d++)
            {
                acc += phi * cpv.Parameters[d] * GradV[d];
            }

            return -acc;
        }
    }

    /// <summary>
    /// Interior Penalty Flux, with Dirichlet boundary conditions for variable 'mu'
    /// </summary>
    class mu_Diffusion : BoSSS.Solution.NSECommon.SIPLaplace
    {

        public mu_Diffusion(int D, double penalty_const, MultidimensionalArray cj, double __cahn, BoundaryCondMap<BoundaryType> __boundaryCondMap)
            : base(penalty_const, cj, "phi") // note: in the equation for 'mu', we have the Laplacian of 'phi'
        {
            m_D = D;
            m_cahn = __cahn * __cahn;
            m_boundaryCondMap = __boundaryCondMap;
            m_bndFunc = m_boundaryCondMap.bndFunction["phi"];
        }

        double m_cahn;
        int m_D;
        public override IList<string> ParameterOrdering => VariableNames.VelocityVector(m_D).Cat("phi0");
        BoundaryCondMap<BoundaryType> m_boundaryCondMap;
        Func<double[], double, double>[] m_bndFunc;

        protected override double g_Diri(ref CommonParamsBnd inp)
        {
            double UxN = 0;
            for (int d = 0; d < m_D; d++)
            {
                UxN += (inp.Parameters_IN[d]) * inp.Normal[d];
            }

            double v;
            if (UxN >= 0)
            {
                v = 0.0;
            }
            else
            {
                // for now set inflow value to -1.0
                v = -1.0; // m_bndFunc[inp.EdgeTag](inp.X, inp.time);
            }
            return v;
        }

        protected override double g_Neum(ref CommonParamsBnd inp)
        {
            return 0.0;
        }

        protected override bool IsDirichlet(ref CommonParamsBnd inp)
        {
            BoundaryType edgeType = m_boundaryCondMap.EdgeTag2Type[inp.EdgeTag];
            switch (edgeType)
            {
                case BoundaryType.Wall:
                case BoundaryType.Slip:
                case BoundaryType.SlipSymmetry:
                    // a Dicirchlet b.c. for 'c' mean a Neumann b.c. for 'phi'
                    return false;

                case BoundaryType.Inlet:
                case BoundaryType.Outlet:
                case BoundaryType.Outflow:
                case BoundaryType.Pressure_Dirichlet:
                default:
                    return true;
            }
        }

        public override double Nu(double[] x, double[] p, int jCell)
        {
            return -m_cahn;
        }

        /// <summary>
        /// Integrand on boundary mesh edges of the SIP
        /// </summary>
        public override double BoundaryEdgeForm(ref Foundation.CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA)
        {
            double Acc = 0.0;

            double pnlty = 2 * this.GetPenalty(inp.jCellIn, -1);//, inp.GridDat.Cells.cj);
            double nuA = this.Nu(inp.X, inp.Parameters_IN, inp.jCellIn);

            if (this.IsDirichlet(ref inp))
            {
                // inhom. Dirichlet b.c.
                // +++++++++++++++++++++

                double g_D = this.g_Diri(ref inp);                

                // inflow / outflow
                if (g_D != 0)
                {
                    for (int d = 0; d < inp.D; d++)
                    {
                        double nd = inp.Normal[d];
                        Acc += (nuA * _Grad_uA[0, d]) * (_vA) * nd;        // consistency
                        Acc += (nuA * _Grad_vA[d]) * (_uA[0] - g_D) * nd;  // symmetry
                    }
                    Acc *= this.m_alpha;

                    Acc -= nuA * (_uA[0] - g_D) * (_vA - 0) * pnlty; // penalty
                }
                else
                {
                    for (int d = 0; d < inp.D; d++)
                    {
                        double nd = inp.Normal[d];
                        //Acc += (nuA * _Grad_uA[0, d]) * (_vA) * nd;        // consistency 
                        //Acc += (nuA * 0.0) * (_vA) * nd;        // consistency Neumann
                    }

                    //double D0 = 1.0;
                    //Acc -= (nuA * _vA) * D0 * (_uA[0]-inp.Parameters_IN[inp.D])/ 0.01;        // Outflow boundary condition see S. Dong

                    Acc *= this.m_alpha;

                }
            }
            else
            {

                double g_N = this.g_Neum(ref inp);

                Acc += nuA * g_N * _vA * this.m_alpha;
            }
            return Acc;
        }
    }

    /// <summary>
    /// nonlinear source term in the 'phi'-equation
    /// </summary>
    class mu_Source : IVolumeForm
    {

        public mu_Source(bool __inc, double _cahn = 0.0)
        {
            m_inc = __inc;
            m_offset = 0.0;//_cahn * 0.488;
        }

        bool m_inc;
        double m_offset;

        public TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public IList<string> ArgumentOrdering => new[] { "mu", "phi" };

        public IList<string> ParameterOrdering => new[] { "phi0" };

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV)
        {

            double mu = U[0];
            // values seem shifted without this offset hack
            double phi = U[1] + m_offset;
            double phi0 = cpv.Parameters[0];

            double Acc = 0;
            if (m_inc == false)
            {
                Acc += -mu;
            }
            else
            {
                Acc += -mu;

                Acc += (3 * phi0 * phi0 - 1.0) * phi - 2 * Math.Pow(phi0, 3.0); // linearized around phi0 (Taylor expansion)
            }


            return Acc * V;
        }
    }


}

