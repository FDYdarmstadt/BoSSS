using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Aggregation;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.Voronoi;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.SipPoisson
{
    class SmoothIpFlux : IEdgeForm, IVolumeForm, ISupportsJacobianComponent
    {
        string argumentName;

        MultidimensionalArray inverseLengthScales;

        double penalty;

        public SmoothIpFlux(
            string argumentName,
            double penalty)
        {
            this.argumentName = argumentName;
            this.penalty = penalty;
        }

        public void SetGrid(AggregationGrid grid)
        {
            inverseLengthScales = ((AggregationGridData)grid.iGridData).AncestorGrid.Cells.cj;
        }

        public void SetGrid(GridCommons grid)
        {
            inverseLengthScales = grid.GridData.Cells.cj;
        }

        IList<string> IEquationComponent.ArgumentOrdering {
            get {
                return new string[] { argumentName };
            }
        }

        IList<string> IEquationComponent.ParameterOrdering {
            get {
                return new string[] { "Tex" };
            }
        }

        TermActivationFlags IEdgeForm.BoundaryEdgeTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV;
            }
        }

        TermActivationFlags IEdgeForm.InnerEdgeTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV;
            }
        }

        TermActivationFlags IVolumeForm.VolTerms {
            get {
                return (TermActivationFlags.GradUxGradV);
            }
        }

        protected virtual double GetPenalty(int jCellIn, int jCellOut)
        {
            double cj_in = inverseLengthScales[jCellIn];
            double mu = penalty * cj_in;
            if (jCellOut >= 0)
            {
                double cj_out = inverseLengthScales[jCellOut];
                mu = Math.Max(mu, penalty * cj_out);
            }
            return mu;
        }

        public double Nu(double[] x, double[] p, int jCell)
        {
            return 1.0;
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA)
        {
            return 0.0;
        }

        virtual public IEquationComponent[] GetJacobianComponents(int SpatialDimension)
        {
            return new[] { this };
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB)
        {
            double Acc = 0.0;

            double pnlty = this.GetPenalty(inp.jCellIn, inp.jCellOut);
            double nuA = this.Nu(inp.X, inp.Parameters_IN, inp.jCellIn);
            double nuB = this.Nu(inp.X, inp.Parameters_OUT, inp.jCellOut);
            double edgeValue = (inp.Parameters_IN[0] + inp.Parameters_OUT[0]) / 2;

            for (int d = 0; d < inp.D; d++)
            {
                Acc += 0.5 * (nuA * _Grad_uA[0, d] + nuB * _Grad_uB[0, d]) * (_vA - _vB) * inp.Normal[d];  // consistency term
                Acc += 0.5 * (nuA * _Grad_vA[d] + nuB * _Grad_vB[d]) * (_uA[0] - _uB[0]) * inp.Normal[d];  // symmetry term
            }

            double nuMax = (Math.Abs(nuA) > Math.Abs(nuB)) ? nuA : nuB;


            Acc -= (_uA[0] - _uB[0]) * (_vA - _vB) * pnlty * nuMax; // penalty term


            return Acc;
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV)
        {
            double acc = 0;
            for (int d = 0; d < cpv.D; d++)
                acc += GradU[0, d] * GradV[d] * this.Nu(cpv.Xglobal, cpv.Parameters, cpv.jCell);
            return acc;
        }
    }
}
