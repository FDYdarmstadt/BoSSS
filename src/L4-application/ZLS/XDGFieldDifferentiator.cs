using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using ilPSP;
using ilPSP.Tracing;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using static BoSSS.Foundation.XDG.XDGField;

namespace ZwoLevelSetSolver {
    static class XDGFieldDifferentiator {

        public static void LaplacianByFlux(double alpha, XDGField target, XDGField source, string species, XDGField tmp = null) {
            using(new FuncTrace()) {
                if(tmp == null)
                    tmp = new XDGField(target.Basis, "tmp");

                int D = target.GridDat.SpatialDimension;
                for(int d = 0; d < D; d++) {
                    tmp.Clear();
                    DerivativeByFlux(alpha, tmp, source, source.Basis.Tracker, d, species);
                    DerivativeByFlux(alpha, target, tmp, source.Basis.Tracker, d, species);
                }
            }
        }

        public static void DivergenceByFlux(double alpha, XDGField target, XDGField[] source, LevelSetTracker lstrkr, string species) {
            for(int d = 0; d < source.Length; d++)
                DerivativeByFlux(alpha, target, source[d], lstrkr, d, species);
        }

        public static void DerivativeByFlux(double alpha, XDGField target, XDGField source, LevelSetTracker lstrkr,  int d, string species) {
            int D = target.Basis.GridDat.SpatialDimension;
            if(d < 0 || d >= D)
                throw new ArgumentException("spatial dimension out of range.", "d");
            MPICollectiveWatchDog.Watch(csMPI.Raw._COMM.WORLD);

            XSpatialOperatorMk2 d_dx = new XSpatialOperatorMk2(1, 1, QuadOrderFunc.Linear(), new string[] { species }, "in", "out");

            var flux = new DerivativeFlux(d);
            d_dx.EquationComponents["out"].Add(flux);
            foreach(int i in lstrkr.GetLevelSetIndicesOfSpecies(species)) {
                var pairs = lstrkr.GetSpeciesPairsSeparatedByLevSet(i);
                foreach(var pair in pairs) {
                    string otherSpecies = pair.Item1 != species ? pair.Item1 : pair.Item2;
                    var lsFlux = new DerivativeFluxBoundary(d, i, species, otherSpecies);
                    d_dx.EquationComponents["out"].Add(lsFlux);
                }
            }
            
            d_dx.Commit();

            var ev = d_dx.GetEvaluatorEx(
                lstrkr,
                new CoordinateMapping(source), null, target.Mapping);

            ev.Evaluate<CoordinateVector>(alpha, 1.0, target.CoordinateVector);
        }
    }

    class DerivativeFlux : IVolumeForm, IEdgeForm {

        int d;
        public DerivativeFlux(int d) {
            this.d = d;
        }
        
        public TermActivationFlags VolTerms => TermActivationFlags.UxGradV;

        public IList<string> ArgumentOrdering => new string[] { "in"};

        public IList<string> ParameterOrdering => new string[] {};

        public TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.UxV;

        public TermActivationFlags InnerEdgeTerms => TermActivationFlags.UxV;

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return _uA[0] * inp.Normal[d] * _vA;
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            return 0.5 *(_uIN[0]+_uOUT[0]) * inp.Normal[d] * (_vIN - _vOUT);
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            return -U[0] * GradV[d];
        }
    }

    class DerivativeFluxBoundary : ILevelSetForm {

        int d;
        string species;
        string otherSpecies;
        int levelSetNumber;
        public DerivativeFluxBoundary(int d, int levelSetNumber, string species, string otherSpecies) {
            this.d = d;
            this.species = species;
            this.otherSpecies = otherSpecies;
            this.levelSetNumber = levelSetNumber;
        }

        public int LevelSetIndex => levelSetNumber;

        public string PositiveSpecies => species;

        public string NegativeSpecies => otherSpecies;

        public TermActivationFlags LevelSetTerms => TermActivationFlags.UxV;

        public IList<string> ArgumentOrdering => new string[] { "in" };

        public IList<string> ParameterOrdering => new string[] {};

        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            return _uOUT[0] * inp.Normal[d] * -_vOUT;
        }
    }
}
