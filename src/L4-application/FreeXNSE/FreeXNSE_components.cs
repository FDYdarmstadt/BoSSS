using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.LevelSetTools.PhasefieldLevelSet;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.Utils;
using MathNet.Numerics.Distributions;
using Microsoft.CodeAnalysis.Operations;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FreeXNSE {

    public abstract class BulkComponent_FreeXNSE : IEdgeForm, IVolumeForm, ISpeciesFilter, ISupportsJacobianComponent {

        protected int D;
        protected string[] spc;
        protected FreeXNSE_BoundaryCondMap BcMap;

        protected BulkComponent_FreeXNSE(int D, string[] spc, FreeXNSE_BoundaryCondMap BcMap) {
            this.D = D;
            this.spc = spc;
            this.BcMap = BcMap;
        }

        public abstract IList<string> ArgumentOrdering { get; }

        public abstract IList<string> ParameterOrdering { get; }

        public abstract TermActivationFlags BoundaryEdgeTerms { get; }

        public abstract TermActivationFlags InnerEdgeTerms { get; }

        public abstract TermActivationFlags VolTerms { get; }

        public string ValidSpecies => spc.ElementAt(0);

        public abstract double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV);
        public abstract double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA);
        public abstract double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT);
        public abstract IEquationComponent[] GetJacobianComponents(int SpatialDimension);
    }

    public abstract class InterfaceComponent_FreeXNSE : ILevelSetForm, ISupportsJacobianComponent {

        public static bool FixedInterface = false;

        protected int D;
        protected string[] spc;
        protected FreeXNSE_BoundaryCondMap BcMap;

        protected InterfaceComponent_FreeXNSE(int D, string[] spc, FreeXNSE_BoundaryCondMap BcMap) {
            this.D = D;
            this.spc = spc;
            this.BcMap = BcMap;
        }

        public abstract IList<string> ArgumentOrdering { get; }

        public abstract IList<string> ParameterOrdering { get; }

        public abstract TermActivationFlags LevelSetTerms { get; }

        public string PositiveSpecies => spc.ElementAt(1);

        public string NegativeSpecies => spc.ElementAt(0);

        public int LevelSetIndex => 0;

        public abstract double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT);

        public abstract IEquationComponent[] GetJacobianComponents(int SpatialDimension);
    }

    public class BulkComponent_VelocityDivergence_Central : BulkComponent_FreeXNSE {

        public BulkComponent_VelocityDivergence_Central(int D, string[] spc, FreeXNSE_BoundaryCondMap BcMap) : base(D, spc, BcMap) {
        }

        public override IList<string> ArgumentOrdering => VariableNames.VelocityVector(D);

        public override IList<string> ParameterOrdering => new string[0];

        public override TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public override TermActivationFlags InnerEdgeTerms => TermActivationFlags.UxV;

        public override TermActivationFlags VolTerms => TermActivationFlags.GradUxV;

        public override double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            for(int d = 0; d < D; d++) {
                acc += GradU[d,d];
            }
            return acc * V;
        }
        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            FreeXNSE_BcType edgType = BcMap.EdgeTag2Type[inp.EdgeTag];
            switch(edgType) {
                default:
                case FreeXNSE_BcType.Dirichlet: {
                    double acc = 0;
                    for(int d = 0; d < D; d++) {
                        double uD = BcMap.bndFunction[ArgumentOrdering[d] + "#" + ValidSpecies][inp.EdgeTag](inp.X, inp.time);
                        acc -= (_uA[d] - uD) * inp.Normal[d];
                    }
                    return acc * _vA;
                }
                case FreeXNSE_BcType.Neumann: {
                    return 0.0;
                }
                case FreeXNSE_BcType.Robin: {
                    goto case FreeXNSE_BcType.Dirichlet;
                }
            }
        }
        public override double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double acc = 0;
            for(int d = 0; d < D; d++) {
                acc -= (_uIN[d] - _uOUT[d]) * inp.Normal[d];
            }
            return acc * 0.5 * (_vIN + _vOUT);
        }
        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }

    public class BulkComponent_PressurePenalty : BulkComponent_FreeXNSE, IEquationComponentCoefficient {

        double Re;
        public BulkComponent_PressurePenalty(int D, string[] spc, FreeXNSE_BoundaryCondMap BcMap, double Re) : base(D, spc, BcMap) {
            this.Re = Re;
        }

        public override IList<string> ArgumentOrdering => new string[] { VariableNames.Pressure };

        public override IList<string> ParameterOrdering => new string[0];

        public override TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.None;

        public override TermActivationFlags InnerEdgeTerms => TermActivationFlags.UxV;

        public override TermActivationFlags VolTerms => TermActivationFlags.None;

        public override double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            return 0.0;
        }
        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return 0.0;
        }
        public override double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double acc = 0.0;
            double h = Math.Max(cj[inp.jCellIn], cj[inp.jCellOut]);
            acc += h * (_uIN[0] - _uOUT[0]) * (_vIN - _vOUT);
            return Re / Mu(inp.jCellIn, inp.X) * acc;
        }

        /// <summary>
        /// Cell-wise length scales for the penalty computation.
        /// </summary>
        MultidimensionalArray cj;

        Func<int, double[], double> Mu;

        /// <summary>
        /// Update of penalty length scales.
        /// </summary>
        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            cj = cs.CellLengthScales;

            if(cs.UserDefinedValues.Keys.Contains(Coefficientnames.Bulkviscosityfield)) {
                Mu = (Func<int, double[], double>)cs.UserDefinedValues[Coefficientnames.Bulkviscosityfield];
            }
        }

        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }

    public class InterfaceComponent_VelocityDivergence_Central : InterfaceComponent_FreeXNSE {

        public InterfaceComponent_VelocityDivergence_Central(int D, string[] spc, FreeXNSE_BoundaryCondMap BcMap) : base(D, spc, BcMap) {
        }

        public override IList<string> ArgumentOrdering => VariableNames.VelocityVector(D);

        public override IList<string> ParameterOrdering => new string[0];

        public override TermActivationFlags LevelSetTerms => FixedInterface ? TermActivationFlags.UxV | TermActivationFlags.V : TermActivationFlags.None;

        public override double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {

            if(FixedInterface) {
                double acc = 0;
                for(int d = 0; d < inp.D; d++) {
                    acc += (_uIN[d] - 0.0) * inp.Normal[d] * _vIN;
                }
                return acc;
            }

            return 0.0;
        }
        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }

    public class BulkComponent_Convective_LaxFriedrich : BulkComponent_FreeXNSE {

        int d;

        public BulkComponent_Convective_LaxFriedrich(int d, int D, string[] spc, FreeXNSE_BoundaryCondMap BcMap) : base(D, spc, BcMap) {
            this.d = d;
        }

        public override IList<string> ArgumentOrdering => VariableNames.VelocityVector(D);

        public override IList<string> ParameterOrdering => VariableNames.Velocity0MeanVector(D);

        public override TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public override TermActivationFlags InnerEdgeTerms => TermActivationFlags.UxV;

        public override TermActivationFlags VolTerms => TermActivationFlags.UxGradV;

        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            FreeXNSE_BcType edgType = BcMap.EdgeTag2Type[inp.EdgeTag];
            switch(edgType) {
                default:
                case FreeXNSE_BcType.Dirichlet: {
                    double acc = 0;
                    double[] uD = new double[D];
                    for(int dd = 0; dd < D; dd++) {
                        uD[dd] = BcMap.bndFunction[ArgumentOrdering[dd] + "#" + ValidSpecies][inp.EdgeTag](inp.X, inp.time);
                    }
                    for(int dd = 0; dd < D; dd++) {
                        acc += 0.5 * (_uA[d] * _uA[dd] + uD[d] * uD[dd]) * inp.Normal[dd];
                    }

                    double lambda = LambdaConvection.GetLambda(inp.Parameters_IN, inp.Normal, false);
                    acc += lambda * (_uA[d] - uD[d]); // no additional tuning for LF scheme width (1.0)

                    return acc * (_vA);
                }
                case FreeXNSE_BcType.Neumann: {
                    double acc = 0;
                   
                    for(int dd = 0; dd < D; dd++) {
                        acc += (_uA[d] * _uA[dd]) * inp.Normal[dd];
                    }

                    return acc * (_vA);
                }
                case FreeXNSE_BcType.Robin: {
                    goto case FreeXNSE_BcType.Dirichlet;
                }
            }

        }
       
        public override double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double acc = 0;
            for(int dd = 0; dd < D; dd++) {
                acc += 0.5 * (_uIN[d] * _uIN[dd] + _uOUT[d] * _uOUT[dd]) * inp.Normal[dd];
            }

            double lambda = Math.Max(LambdaConvection.GetLambda(inp.Parameters_IN, inp.Normal, false), LambdaConvection.GetLambda(inp.Parameters_OUT, inp.Normal, false));
            acc += lambda * (_uIN[d] - _uOUT[d]); // no additional tuning for LF scheme width (1.0)

            return acc * (_vIN - _vOUT);
        }

        public override double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            for(int dd = 0; dd < D; dd++) {
                acc += -U[d] * U[dd] * GradV[dd];
            }
            return acc;
        }

        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] {new VolumeFormDifferentiator(this, SpatialDimension), new EdgeFormDifferentiator(this, SpatialDimension)};
        }

    }

    public class InterfaceComponent_Convective_LaxFriedrich : InterfaceComponent_FreeXNSE {

        int d;

        public InterfaceComponent_Convective_LaxFriedrich(int d, int D, string[] spc, FreeXNSE_BoundaryCondMap BcMap) : base(D, spc, BcMap) {
            this.d = d;
        }

        public override IList<string> ArgumentOrdering => VariableNames.VelocityVector(D);

        public override IList<string> ParameterOrdering => VariableNames.Velocity0MeanVector(D);

        public override TermActivationFlags LevelSetTerms => FixedInterface ? TermActivationFlags.UxV | TermActivationFlags.V : TermActivationFlags.UxV;

        public override double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double acc = 0;

            if(FixedInterface) {
                double[] uD = new double[D];

                for(int dd = 0; dd < D; dd++) {
                    acc += 0.5 * (_uIN[d] * _uIN[dd] + uD[d] * uD[dd]) * inp.Normal[dd];
                }

                double lambda = LambdaConvection.GetLambda(inp.Parameters_IN, inp.Normal, false);
                acc += lambda * (_uIN[d] - uD[d]); // no additional tuning for LF scheme width (1.0)

                return acc * (_vIN);
            }

            for(int dd = 0; dd < D; dd++) {
                acc += (_uIN[d] * _uIN[dd]) * inp.Normal[dd];
            }
            return acc * _vIN;
        }

        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { new LevelSetFormDifferentiator(this, SpatialDimension)};
        }

    }

    public class BulkComponent_Convective_Temam : BulkComponent_FreeXNSE {

        int d;

        public BulkComponent_Convective_Temam(int d, int D, string[] spc, FreeXNSE_BoundaryCondMap BcMap) : base(D, spc, BcMap) {
            this.d = d;
        }

        public override IList<string> ArgumentOrdering => VariableNames.VelocityVector(D);

        public override IList<string> ParameterOrdering => new string[0];

        public override TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public override TermActivationFlags InnerEdgeTerms => TermActivationFlags.UxV;

        public override TermActivationFlags VolTerms => TermActivationFlags.UxV | TermActivationFlags.GradUxV;

        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            FreeXNSE_BcType edgType = BcMap.EdgeTag2Type[inp.EdgeTag];
            switch(edgType) {
                default:
                case FreeXNSE_BcType.Dirichlet: {
                    double acc = 0;
                    double[] uD = new double[D];
                    for(int dd = 0; dd < D; dd++) {
                        uD[dd] = BcMap.bndFunction[ArgumentOrdering[dd] + "#" + ValidSpecies][inp.EdgeTag](inp.X, inp.time);
                    }
                    for(int dd = 0; dd < D; dd++) {
                        acc -= 0.5 * ((_uA[dd] - uD[dd]) * inp.Normal[dd] * 0.5 * (_uA[d] * _vA));
                    }

                    return acc;
                }
                case FreeXNSE_BcType.Neumann: {
                    double acc = 0;

                    return acc;
                }
                case FreeXNSE_BcType.Robin: {
                    goto case FreeXNSE_BcType.Dirichlet;
                }
            }

        }

        public override double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double acc = 0;
            for(int dd = 0; dd < D; dd++) {
                acc -= 0.5 * (_uIN[dd] + _uOUT[dd]) * inp.Normal[dd] * (_uIN[d] - _uOUT[d]) * 0.5 * (_vIN + _vOUT);
            }
            for(int dd = 0; dd < D; dd++) {
                acc -= 0.5 * ((_uIN[dd] - _uOUT[dd]) * inp.Normal[dd] * 0.5 * (_uIN[d] * _vIN + _uOUT[d] * _vOUT));
            }
            return acc;
        }

        public override double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            for(int dd = 0; dd < D; dd++) {
                acc += U[dd] * GradU[d, dd] * V;
            }
            for(int dd = 0; dd < D; dd++) {
                acc += 0.5 * GradU[dd,dd] * U[d] * V;
            }
            return acc;
        }

        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { new VolumeFormDifferentiator(this, SpatialDimension), new EdgeFormDifferentiator(this, SpatialDimension) };
        }

    }

    public class BulkComponent_Convective_ConservativeTemam : BulkComponent_FreeXNSE {

        int d;

        public BulkComponent_Convective_ConservativeTemam(int d, int D, string[] spc, FreeXNSE_BoundaryCondMap BcMap) : base(D, spc, BcMap) {
            this.d = d;
        }

        public override IList<string> ArgumentOrdering => VariableNames.VelocityVector(D);

        public override IList<string> ParameterOrdering => new string[0];

        public override TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public override TermActivationFlags InnerEdgeTerms => TermActivationFlags.UxV;

        public override TermActivationFlags VolTerms => TermActivationFlags.UxGradV;

        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            FreeXNSE_BcType edgType = BcMap.EdgeTag2Type[inp.EdgeTag];
            switch(edgType) {
                default:
                case FreeXNSE_BcType.Dirichlet: {
                    double acc = 0;
                    double[] uD = new double[D];
                    for(int dd = 0; dd < D; dd++) {
                        uD[dd] = BcMap.bndFunction[ArgumentOrdering[dd] + "#" + ValidSpecies][inp.EdgeTag](inp.X, inp.time);
                    }
                    for(int dd = 0; dd < D; dd++) {
                        acc += 0.5 * (uD[dd] * uD[dd]) * inp.Normal[d] * (_vA);
                    }

                    return acc;
                }
                case FreeXNSE_BcType.Neumann: {
                    double acc = 0;
                    for(int dd = 0; dd < D; dd++) {
                        acc += _uA[dd] * _uA[d] * inp.Normal[dd] * _vA;
                    }
                    for(int dd = 0; dd < D; dd++) {
                        acc += 0.5 * _uA[dd] * _uA[dd] * inp.Normal[d] * _vA;
                    }
                    return acc;
                }
                case FreeXNSE_BcType.Robin: {
                    goto case FreeXNSE_BcType.Dirichlet;
                }
            }

        }

        public override double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double acc = 0;
            for(int dd = 0; dd < D; dd++) {
                acc += 0.5 * (_uIN[dd] + _uOUT[dd]) * inp.Normal[dd] * 0.5 * (_uIN[d] + _uOUT[d]) * (_vIN - _vOUT);
            }
            for(int dd = 0; dd < D; dd++) {
                acc += 0.5 * 0.5 * (_uIN[dd] * _uIN[dd] + _uOUT[dd]*_uOUT[dd]) * inp.Normal[d] * (_vIN - _vOUT);
            }
            return acc;
        }

        public override double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            for(int dd = 0; dd < D; dd++) {
                acc -= U[d] * U[dd] * GradV[dd];
            }
            for(int dd = 0; dd < D; dd++) {
                acc -= 0.5 * U[dd] * U[dd] * GradV[d];
            }
            return acc;
        }

        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { new VolumeFormDifferentiator(this, SpatialDimension), new EdgeFormDifferentiator(this, SpatialDimension) };
        }

    }


    public class InterfaceComponent_Convective_Temam : InterfaceComponent_FreeXNSE {

        int d;

        public InterfaceComponent_Convective_Temam(int d, int D, string[] spc, FreeXNSE_BoundaryCondMap BcMap) : base(D, spc, BcMap) {
            this.d = d;
        }

        public override IList<string> ArgumentOrdering => VariableNames.VelocityVector(D);

        public override IList<string> ParameterOrdering => new string[0];

        public override TermActivationFlags LevelSetTerms => FixedInterface ? TermActivationFlags.UxV | TermActivationFlags.V : TermActivationFlags.None;

        public override double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {

            if(FixedInterface) {
                double acc = 0;
                double[] uD = new double[D];

                for(int dd = 0; dd < D; dd++) {
                    acc -= 0.5 * ((_uIN[dd] - uD[dd]) * inp.Normal[dd] * 0.5 * (_uIN[d] * _vIN));
                }

                return acc;
            }

            return 0.0;
        }

        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

    }

    public class InterfaceComponent_Convective_ConservativeTemam : InterfaceComponent_FreeXNSE {

        int d;

        public InterfaceComponent_Convective_ConservativeTemam(int d, int D, string[] spc, FreeXNSE_BoundaryCondMap BcMap) : base(D, spc, BcMap) {
            this.d = d;
        }

        public override IList<string> ArgumentOrdering => VariableNames.VelocityVector(D);

        public override IList<string> ParameterOrdering => new string[0];

        public override TermActivationFlags LevelSetTerms => FixedInterface ? TermActivationFlags.V : TermActivationFlags.UxV;

        public override double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double acc = 0;

            if(FixedInterface) {
                double[] uD = new double[D];

                for(int dd = 0; dd < D; dd++) {
                    acc += 0.5 * (uD[dd] * uD[dd]) * inp.Normal[d] * (_vIN);
                }

                return acc;
            }

            for(int dd = 0; dd < D; dd++) {
                acc += 0.5 * _uIN[dd] * _uIN[dd] * inp.Normal[d] * _vIN;
            }
            return acc;
        }

        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { new LevelSetFormDifferentiator(this, SpatialDimension) };
        }

    }

    public class BulkComponent_PressureGradient_Central : BulkComponent_FreeXNSE {

        int d;

        public BulkComponent_PressureGradient_Central(int d, int D, string[] spc, FreeXNSE_BoundaryCondMap BcMap) :  base(D, spc, BcMap) {
            this.d = d;
        }

        public override IList<string> ArgumentOrdering => new string[] { VariableNames.Pressure };

        public override IList<string> ParameterOrdering => new string[0];

        public override TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public override TermActivationFlags InnerEdgeTerms => TermActivationFlags.UxV;

        public override TermActivationFlags VolTerms => TermActivationFlags.UxGradV;

        public override double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;            
            acc += GradV[d];            
            return -acc * U[0];
        }
        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            FreeXNSE_BcType edgType = BcMap.EdgeTag2Type[inp.EdgeTag];
            switch(edgType) {
                default:
                case FreeXNSE_BcType.Dirichlet: {
                    double acc = 0;
                    acc += _vA * inp.Normal[d];
                    return acc * _uA[0];
                }
                case FreeXNSE_BcType.Neumann: {
                    return 0.0; // inhomogenous component implemented in viscous term
                }
                case FreeXNSE_BcType.Robin: {
                    goto case FreeXNSE_BcType.Dirichlet;
                }
            }
        }
        public override double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double acc = 0;
            acc += (_vIN - _vOUT) * inp.Normal[d];            
            return acc * 0.5 * (_uIN[0] + _uOUT[0]);
        }

        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }

    public class InterfaceComponent_PressureGradient_Central : InterfaceComponent_FreeXNSE {

        int d;

        public InterfaceComponent_PressureGradient_Central(int d, int D, string[] spc, FreeXNSE_BoundaryCondMap BcMap) : base(D, spc, BcMap) {
            this.d = d;
        }

        public override IList<string> ArgumentOrdering => new string[] { VariableNames.Pressure };

        public override IList<string> ParameterOrdering => new string[0];

        public override TermActivationFlags LevelSetTerms => FixedInterface ? TermActivationFlags.UxV : TermActivationFlags.None;

        public override double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {

            if(FixedInterface) {
                double acc = 0;
                acc += _vIN * inp.Normal[d];
                return acc * _uIN[0];
            }

            return 0.0;
        }
        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }

    public class BulkComponent_Viscous_IP : BulkComponent_FreeXNSE, IEquationComponentCoefficient {

        int d;
        double symmetry; // 1 SIP, -1 NIP, 0 IIP
        double penaltyBase; 
        double Re;

        public BulkComponent_Viscous_IP(int d, int D, string[] spc, FreeXNSE_BoundaryCondMap BcMap, double Re, double symmetry = 1.0, double penalty = 4.0) : base(D, spc, BcMap) {
            this.d = d;
            this.symmetry = symmetry;
            this.penaltyBase = penalty;
            this.Re = Re;
        }

        public override IList<string> ArgumentOrdering => VariableNames.VelocityVector(D);

        public override IList<string> ParameterOrdering => new string[0] {};

        public override TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.V | TermActivationFlags.GradV;

        public override TermActivationFlags InnerEdgeTerms => TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.V;

        public override TermActivationFlags VolTerms => TermActivationFlags.GradUxGradV;

        public override double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;
            for(int dd = 0; dd < D; dd++) {
                acc += GradU[d, dd] * GradV[dd];
            }            
            return 1/Re * Mu(cpv.jCell, cpv.Xglobal) * acc;
        }
        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            FreeXNSE_BcType edgType = BcMap.EdgeTag2Type[inp.EdgeTag];
            switch(edgType) {
                default:
                case FreeXNSE_BcType.Dirichlet: {
                    double acc = 0.0;
                    double uD = BcMap.bndFunction[ArgumentOrdering[d] + "#" + ValidSpecies][inp.EdgeTag](inp.X, inp.time);

                    for(int dd = 0; dd < inp.D; dd++) {
                        acc -= (_Grad_uA[d, dd]) * _vA * inp.Normal[dd];  // consistency term  
                        acc -= symmetry * (_Grad_vA[dd]) * (_uA[d] - uD) * inp.Normal[dd];  // symmetry term
                    }

                    if(penaltyBase > 0) {
                        double pnlty = this.penalty(inp.GridDat, inp.jCellIn, -1, inp.iEdge);
                        acc += pnlty * (_uA[d] - uD) * _vA; // penalty term 
                    }
                    return 1.0 / Re * Mu(inp.jCellIn, inp.X) * acc;
                }
                case FreeXNSE_BcType.Neumann: {
                    double acc = 0;
                    double uD = BcMap.bndFunction[ArgumentOrdering[d] + "#" + ValidSpecies][inp.EdgeTag](inp.X, inp.time);
                    acc += -uD * inp.Normal[d] * _vA;                    
                    return acc;
                }
                case FreeXNSE_BcType.Robin: {
                    double acc = 0.0;
                    double m_beta = Beta(inp.jCellIn, inp.X);
                    if(m_beta < 0 || m_beta == double.PositiveInfinity)
                        goto case FreeXNSE_BcType.Dirichlet;

                    for(int dN = 0; dN < inp.D; dN++) {
                        double uD = BcMap.bndFunction[ArgumentOrdering[dN] + "#" + ValidSpecies][inp.EdgeTag](inp.X, inp.time);

                        for(int dD = 0; dD < inp.D; dD++) {
                            // consistency
                            acc -= 1.0 / Re * Mu(inp.jCellIn, inp.X) * (inp.Normal[dN] * _Grad_uA[dN, dD] * inp.Normal[dD]) * (_vA * inp.Normal[d]);
                            // symmetry
                            acc -= symmetry * 1.0 / Re * Mu(inp.jCellIn, inp.X) * (inp.Normal[d] * _Grad_vA[dD] * inp.Normal[dD]) * (_uA[dN] - uD) * inp.Normal[dN];
                        }

                        // penalty
                        if(penaltyBase > 0) {
                            double pnlty = this.penalty(inp.GridDat, inp.jCellIn, -1, inp.iEdge);
                            acc += 1.0 / Re * Mu(inp.jCellIn, inp.X) * ((_uA[dN] - uD) * inp.Normal[dN]) * ((_vA - 0) * inp.Normal[d]) * pnlty;
                        }
                    }


                    double[,] P = new double[D, D];
                    for(int d1 = 0; d1 < D; d1++) {
                        for(int d2 = 0; d2 < D; d2++) {
                            double nn = inp.Normal[d1] * inp.Normal[d2];
                            if(d1 == d2) {
                                P[d1, d2] = 1 - nn;
                            } else {
                                P[d1, d2] = -nn;
                            }
                        }
                    }

                    // tangential dissipation force term
                    for(int d1 = 0; d1 < D; d1++) {
                        for(int d2 = 0; d2 < D; d2++) {
                            double uD = BcMap.bndFunction[ArgumentOrdering[d2] + "#" + ValidSpecies][inp.EdgeTag](inp.X, inp.time);
                            acc += (m_beta * P[d1, d2] * (_uA[d2] - uD)) * (P[d1, d] * _vA);
                        }
                    }

                    return acc;
                }
            }
        }
        public override double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double acc = 0.0;

            for(int dd = 0; dd < inp.D; dd++) {
                acc -= 0.5 * (_Grad_uIN[d, dd] + _Grad_uOUT[d, dd]) * (_vIN - _vOUT) * inp.Normal[dd];  // consistency term  
                acc -= symmetry * 0.5 * (_Grad_vIN[dd] + _Grad_vOUT[dd]) * (_uIN[d] - _uOUT[d]) * inp.Normal[dd];  // symmetry term
            }
            if(penaltyBase > 0) {
                double pnlty = this.penalty(inp.GridDat, inp.jCellIn, inp.jCellOut, inp.iEdge);
                acc += pnlty * (_uIN[d] - _uOUT[d]) * (_vIN - _vOUT); // penalty term 
            }
            return 1/Re * Mu(inp.jCellIn, inp.X) * acc;
        }

        /// <summary>
        /// penalty adapted for spatial dimension and DG-degree
        /// </summary>
        double m_penalty;

        /// <summary>
        /// computation of penalty parameter according to:
        /// An explicit expression for the penalty parameter of the
        /// interior penalty method, K. Shahbazi, J. of Comp. Phys. 205 (2004) 401-407,
        /// look at formula (7) in cited paper
        /// </summary>
        protected double penalty(IGridData g, int jCellIn, int jCellOut, int iEdge) {
            double penaltySizeFactor_A = 1.0 / cj[jCellIn];
            double penaltySizeFactor_B = jCellOut >= 0 ? 1.0 / cj[jCellOut] : 0;

            double penaltySizeFactor = Math.Max(penaltySizeFactor_A, penaltySizeFactor_B);

            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penaltyBase));
            Debug.Assert(!double.IsInfinity(m_penalty));

            double µ = penaltySizeFactor * m_penalty * penaltyBase;
            if(µ <= 0 || µ.IsNaNorInf()) {
                string errStr = ($"Inf/NaN/Non-Positive in penalty comp: {µ}; (m_penalty = {m_penalty}, m_penalty_base = {penaltyBase}, jCellIn = {jCellIn}, jCellOut = {jCellOut}, cj_in = {cj[jCellIn]}, cj_out = {(jCellOut >= 0 ? 1.0 / cj[jCellOut] : 0)}, penaltySizeFactor_A = {penaltySizeFactor_A}, penaltySizeFactor_B = {penaltySizeFactor_B})");
                throw new ArithmeticException(errStr);
            }
            return µ;
        }

        /// <summary>
        /// Cell-wise length scales for the penalty computation.
        /// </summary>
        MultidimensionalArray cj;

        /// <summary>
        /// Evaluation of friction at a point in a cell
        /// </summary>
        Func<int, double[], double> Beta;

        /// <summary>
        /// Scaling of the viscosity
        /// </summary>
        Func<int, double[], double> Mu;

        /// <summary>
        /// Update of penalty length scales.
        /// </summary>
        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            double _p = DomainDGdeg.Max();

            double penalty_deg_tri = (_p + 1) * (_p + D) / D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes

            m_penalty = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            cj = cs.CellLengthScales;

            if(cs.UserDefinedValues.Keys.Contains(Coefficientnames.Bulkfriction)) {
                Beta = (Func<int, double[], double>)cs.UserDefinedValues[Coefficientnames.Bulkfriction];
            }

            if(cs.UserDefinedValues.Keys.Contains(Coefficientnames.Bulkviscosityfield)) {
                Mu = (Func<int, double[], double>)cs.UserDefinedValues[Coefficientnames.Bulkviscosityfield];
            }
        }


        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }

    public class InterfaceComponent_Viscous_IP : InterfaceComponent_FreeXNSE, IEquationComponentCoefficient {

        int d;
        double symmetry; // 1 SIP, -1 NIP, 0 IIP
        double penaltyBase;
        double Re;

        public InterfaceComponent_Viscous_IP(int d, int D, string[] spc, FreeXNSE_BoundaryCondMap BcMap, double Re, double symmetry = 1.0, double penalty = 4.0) : base(D, spc, BcMap) {
            this.d = d;
            this.symmetry = symmetry;
            this.penaltyBase = penalty;
            this.Re = Re;
        }

        public override IList<string> ArgumentOrdering => VariableNames.VelocityVector(D);

        public override IList<string> ParameterOrdering => new string[0];

        public override TermActivationFlags LevelSetTerms => FixedInterface ? TermActivationFlags.GradUxV | TermActivationFlags.UxGradV | TermActivationFlags.UxV | TermActivationFlags.V | TermActivationFlags.GradV : TermActivationFlags.None;

        public override double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {

            if(FixedInterface) {
                // like freeslip, but with fixed normal velocity!
                double acc = 0.0;
                for(int dN = 0; dN < inp.D; dN++) {
                    double uD = 0.0;

                    for(int dD = 0; dD < inp.D; dD++) {
                        // consistency
                        acc -= 1.0 / Re * Mu(inp.jCellIn, inp.X) * (inp.Normal[dN] * _Grad_uIN[dN, dD] * inp.Normal[dD]) * (_vIN * inp.Normal[d]); // in first approximation we do not care about the pressure!
                        // symmetry
                        acc -= symmetry * 1.0 / Re * Mu(inp.jCellIn, inp.X) * (inp.Normal[d] * _Grad_vIN[dD] * inp.Normal[dD]) * (_uIN[dN] - uD) * inp.Normal[dN];
                    }

                    // penalty
                    if(penaltyBase > 0) {
                        double pnlty = this.penalty(inp.GridDat, inp.jCellIn, -1, inp.iEdge);
                        acc += 1.0 / Re * Mu(inp.jCellIn, inp.X) * ((_uIN[dN] - uD) * inp.Normal[dN]) * ((_vIN - 0) * inp.Normal[d]) * pnlty;
                    }
                }
                return acc;
            }
            return 0.0;
        }

        /// <summary>
        /// penalty adapted for spatial dimension and DG-degree
        /// </summary>
        double m_penalty;

        /// <summary>
        /// computation of penalty parameter according to:
        /// An explicit expression for the penalty parameter of the
        /// interior penalty method, K. Shahbazi, J. of Comp. Phys. 205 (2004) 401-407,
        /// look at formula (7) in cited paper
        /// </summary>
        protected double penalty(IGridData g, int jCellIn, int jCellOut, int iEdge) {
            double penaltySizeFactor_A = 1.0 / cj[jCellIn];
            double penaltySizeFactor_B = jCellOut >= 0 ? 1.0 / cj[jCellOut] : 0;

            double penaltySizeFactor = Math.Max(penaltySizeFactor_A, penaltySizeFactor_B);

            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penaltyBase));
            Debug.Assert(!double.IsInfinity(m_penalty));

            double µ = penaltySizeFactor * m_penalty * penaltyBase;
            if(µ <= 0 || µ.IsNaNorInf()) {
                string errStr = ($"Inf/NaN/Non-Positive in penalty comp: {µ}; (m_penalty = {m_penalty}, m_penalty_base = {penaltyBase}, jCellIn = {jCellIn}, jCellOut = {jCellOut}, cj_in = {cj[jCellIn]}, cj_out = {(jCellOut >= 0 ? 1.0 / cj[jCellOut] : 0)}, penaltySizeFactor_A = {penaltySizeFactor_A}, penaltySizeFactor_B = {penaltySizeFactor_B})");
                throw new ArithmeticException(errStr);
            }
            return µ;
        }

        /// <summary>
        /// Cell-wise length scales for the penalty computation.
        /// </summary>
        MultidimensionalArray cj;

        /// <summary>
        /// Scaling of the viscosity
        /// </summary>
        Func<int, double[], double> Mu;

        /// <summary>
        /// Update of penalty length scales.
        /// </summary>
        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            double _p = DomainDGdeg.Max();

            double penalty_deg_tri = (_p + 1) * (_p + D) / D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes

            m_penalty = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            cj = cs.CellLengthScales;

            if(cs.UserDefinedValues.Keys.Contains(Coefficientnames.Bulkviscosityfield)) {
                Mu = (Func<int, double[], double>)cs.UserDefinedValues[Coefficientnames.Bulkviscosityfield];
            }
        }


        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }
    /// <summary>
    /// Bulk component with respect to the interface
    /// </summary>
    public class InterfaceComponent_SurfaceTension_LaplaceBeltrami : BulkComponent_FreeXNSE, IEquationComponentCoefficient {

        int d;
        double We;

        public InterfaceComponent_SurfaceTension_LaplaceBeltrami(int d, int D, string[] spc, FreeXNSE_BoundaryCondMap bcmap, double We) : base(D, spc, bcmap) {
            this.d = d;
            this.We = We;
        }

        public override IList<string> ArgumentOrdering => VariableNames.VelocityVector(D);

        public override IList<string> ParameterOrdering => VariableNames.NormalVector(D);

        public override TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.UxV | TermActivationFlags.V;

        public override TermActivationFlags InnerEdgeTerms => TermActivationFlags.V;

        public override TermActivationFlags VolTerms => TermActivationFlags.GradV;

        public override double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;

            double[] Nsurf = SurfaceNormal(cpv.Parameters);
            double[,] Psurf = SurfaceProjection(Nsurf);

            for(int dd = 0; dd < D; dd++) {
                acc += 1 / (We) * Sigma(cpv.Xglobal) * Psurf[d, dd] * GradV[dd];
            }

            return acc;
        }


        public override double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {

            double[] EdgeNormal = inp.Normal;
            double[] SurfaceNormal_IN = SurfaceNormal(inp.Parameters_IN);
            double[] SurfaceNormal_OUT = SurfaceNormal(inp.Parameters_OUT);

            double[] Tangente_IN = Tangent(SurfaceNormal_IN, EdgeNormal);
            double[] Tangente_OUT = Tangent(SurfaceNormal_OUT, EdgeNormal);

            double acc = 0.5 * 1/ (We) * Sigma(inp.X) * (Tangente_IN[d] + Tangente_OUT[d]) * (_vA - _vB);

            return -acc;
        }


        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {

            double Flx_InCell = 0;

            FreeXNSE_BcType edgType = BcMap.EdgeTag2Type[inp.EdgeTag];

            switch(edgType) {
                case FreeXNSE_BcType.Dirichlet:
                case FreeXNSE_BcType.Neumann: {

                    double[] EdgeNormal = inp.Normal;
                    double[] SurfaceNormal_IN = SurfaceNormal(inp.Parameters_IN);
                    double[] Tangente_IN = Tangent(SurfaceNormal_IN, EdgeNormal);

                    Flx_InCell = -1/ (We) * Sigma(inp.X) * Tangente_IN[d];

                    break;
                }
                case FreeXNSE_BcType.Robin: {

                    double[] EdgeNormal = inp.Normal;
                    double[] SurfaceNormal_IN = SurfaceNormal(inp.Parameters_IN);
                    double[] Tangente_IN = Tangent(SurfaceNormal_IN, EdgeNormal);

                    double[] PSnI = new double[D]; // projection of surface/level-set normal onto domain boundary tangent
                    for(int d1 = 0; d1 < D; d1++) {
                        for(int d2 = 0; d2 < D; d2++) {
                            double nn = EdgeNormal[d1] * EdgeNormal[d2];
                            if(d1 == d2) {
                                PSnI[d1] += (1 - nn) * SurfaceNormal_IN[d2];
                            } else {
                                PSnI[d1] += -nn * SurfaceNormal_IN[d2];
                            }
                        }
                    }
                    double PSnINorm = PSnI.L2Norm();
                    double[] PSnINormal_IN = PSnI.Normalize(); // line normal: tangential to domain boundary & normal on contact line


                    // isotropic surface tension terms
                    for(int dd = 0; dd < D; dd++) {
                        Flx_InCell -= 1/(We) * Sigma(inp.X) * (EdgeNormal[dd] * Tangente_IN[dd]) * EdgeNormal[d];
                    }


                    if(edgType == FreeXNSE_BcType.Robin) {

                        double alpha = Alpha[0](inp.X);
                        double gamma = Alpha[1](inp.X);
                        double theta0 = Theta(inp.X);
                        double thetaadv = ThetaAdv(inp.X);
                        double thetarec = ThetaRec(inp.X);

                        double CosTheta = -SurfaceNormal_IN.InnerProd(EdgeNormal);
                        bool pinned = false;
                        if(CosTheta < Math.Cos(Math.Max(thetarec, 0)) && CosTheta > Math.Cos(Math.Min(thetaadv, Math.PI))) {
                            pinned = true;
                        }

                        if(pinned) {
                            // Free contactangle condition
                            for(int dd = 0; dd < D; dd++) {
                                Flx_InCell -= 1 / (We) * Sigma(inp.X) * (PSnINormal_IN[dd] * Tangente_IN[dd]) * PSnINormal_IN[d];
                            }

                            // Contactlinevelocity has to be zero
                            // dissipative contact line force
                            // beta*(u*nL)

                            for(int dd = 0; dd < D; dd++) {
                                double g_D = BcMap.bndFunction[ArgumentOrdering[dd] + "#" + ValidSpecies][inp.EdgeTag](inp.X, inp.time);
                                Flx_InCell += 100 * ((_uA[dd] - g_D) * PSnINormal_IN[dd]) * PSnINormal_IN[d]; // suitable penalty??
                            }
                        } else {
                            if(alpha < 0 || alpha == double.PositiveInfinity) {
                                // Free contactangle condition
                                for(int dd = 0; dd < D; dd++) {
                                    Flx_InCell -= 1 / (We) * Sigma(inp.X) * (PSnINormal_IN[dd] * Tangente_IN[dd]) * PSnINormal_IN[d];
                                }
                            } else {
                                // Slipcondition
                                // Young's relation (static contact angle)
                                Flx_InCell -= 1 / (We) * Sigma(inp.X) * Math.Cos(theta0) * PSnINormal_IN[d];

                                // dissipative contact line force
                                // beta*(u*nL)

                                double Ucl = 0;
                                for(int dd = 0; dd < D; dd++) {
                                    double g_D = BcMap.bndFunction[ArgumentOrdering[dd] + "#" + ValidSpecies][inp.EdgeTag](inp.X, inp.time);
                                    Ucl += ((_uA[dd] - g_D) * PSnINormal_IN[dd]);
                                }

                                if(gamma == 1.0) {
                                    Flx_InCell += alpha * Ucl * PSnINormal_IN[d];
                                } else {
                                    Flx_InCell += alpha * Math.Sign(Ucl) * Math.Pow(Math.Abs(Ucl), gamma) * PSnINormal_IN[d];
                                }


                                // idea match contactline force with viscous force, seems not to work... maybe add integration of that stress at the contactline with a small artificial radius?
                                //for(int d1 = 0; d1 < D; d1++) {
                                //    for(int d2 = 0; d2 < D; d2++) {
                                //        Flx_InCell += alpha * (-_uA[D] * (d1 == d2 ? 1.0 : 0.0) + /* 1/Re * */ _Grad_uA[d1, d2]) * EdgeNormal[d2] * PSnINormal_IN[d1] * PSnINormal_IN[d];
                                //    }
                                //}
                            }
                        }
                    }

                    break;
                }
                default:
                break;
            }

            return Flx_InCell * _vA;
        }


        protected static double[] SurfaceNormal(double[] param) {

            double[] N = new double[param.Length];

            for(int d = 0; d < param.Length; d++) {
                N[d] = param[d];
            }

            return N.Normalize();
        }

        protected static double[,] SurfaceProjection(double[] Nsurf) {

            int D = Nsurf.Length;
            double[,] P = new double[D, D];

            for(int d = 0; d < D; d++) {
                for(int dd = 0; dd < D; dd++) {
                    if(dd == d)
                        P[d, dd] = (1 - Nsurf[d] * Nsurf[dd]);
                    else
                        P[d, dd] = (0 - Nsurf[d] * Nsurf[dd]);
                }
            }

            return P;
        }

        protected static double[] Tangent(double[] Nsurf, double[] Nedge) {
            Debug.Assert(Nsurf.Length == Nedge.Length);

            int D = Nsurf.Length;

            double[] tau = new double[D];
            for(int d1 = 0; d1 < D; d1++) {
                for(int d2 = 0; d2 < D; d2++) {
                    double nn = Nsurf[d1] * Nsurf[d2];
                    if(d1 == d2) {
                        tau[d1] += (1 - nn) * Nedge[d2];
                    } else {
                        tau[d1] += -nn * Nedge[d2];
                    }
                }
            }

            return tau.Normalize();
        }

        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            //if(Alpha[1](SpatialDimension.ForLoop(i => 0.0)) == 1 || Alpha[1](SpatialDimension.ForLoop(i => 0.0)) == 0) { // null reference as CoefficientUpdate is not yet invoked.
            //    return new IEquationComponent[] { this }; // only parameter dependent, not present in jacobian, boundary term if slip is on is dependent on velocity
            //} else {
                return new IEquationComponent[] { new VolumeFormDifferentiator(this, SpatialDimension), new EdgeFormDifferentiator(this, SpatialDimension) }; // actually only called if solver level is already marked nonlinear
            //}
        }

        Func<double[], double> Theta;
        Func<double[], double> ThetaAdv;
        Func<double[], double> ThetaRec;
        Func<double[], double>[] Alpha;
        Func<double[], double> Sigma;

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            if(cs.UserDefinedValues.Keys.Contains(Coefficientnames.Contactangle)) {
                Theta = (Func<double[], double>)cs.UserDefinedValues[Coefficientnames.Contactangle];
            }

            if(cs.UserDefinedValues.Keys.Contains(Coefficientnames.Contactangle + "Adv")) {
                ThetaAdv = (Func<double[], double>)cs.UserDefinedValues[Coefficientnames.Contactangle + "Adv"];
            }

            if(cs.UserDefinedValues.Keys.Contains(Coefficientnames.Contactangle + "Rec")) {
                ThetaRec = (Func<double[], double>)cs.UserDefinedValues[Coefficientnames.Contactangle + "Rec"];
            }

            if(cs.UserDefinedValues.Keys.Contains(Coefficientnames.Contactlinefriction)) {
                Alpha = (Func<double[], double>[])cs.UserDefinedValues[Coefficientnames.Contactlinefriction];
            }

            if(cs.UserDefinedValues.Keys.Contains(Coefficientnames.Surfacetensionfield)) {
                Sigma = (Func<double[], double>)cs.UserDefinedValues[Coefficientnames.Surfacetensionfield];
            }
        }
    }

    /// <summary>
    /// Bulk component with respect to the interface
    /// </summary>
    public class InterfaceComponent_SurfaceTension_LaplaceBeltrami_BoussinesqScriven : BulkComponent_FreeXNSE, IEquationComponentCoefficient {

        int d;
        double We;

        public InterfaceComponent_SurfaceTension_LaplaceBeltrami_BoussinesqScriven(int d, int D, string[] spc, FreeXNSE_BoundaryCondMap bcmap, double We) : base(D, spc, bcmap) {
            this.d = d;
            this.We = We;
        }

        public override IList<string> ArgumentOrdering => VariableNames.VelocityVector(D);

        public override IList<string> ParameterOrdering => VariableNames.NormalVector(D);

        public override TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.GradUxV | TermActivationFlags.UxV | TermActivationFlags.V;

        public override TermActivationFlags InnerEdgeTerms => TermActivationFlags.GradUxV | TermActivationFlags.V;

        public override TermActivationFlags VolTerms => TermActivationFlags.GradUxGradV | TermActivationFlags.GradV;

        public override double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double acc = 0;

            double[] Nsurf = SurfaceNormal(cpv.Parameters);
            double[,] Psurf = SurfaceProjection(Nsurf);

            // isotropic part
            for(int dd = 0; dd < D; dd++) {
                acc += 1 / (We) * Sigma(cpv.Xglobal) * Psurf[d, dd] * GradV[dd];
            }

            // surface dilatational part
            double surfDivU = 0.0;
            for(int d1 = 0; d1 < D; d1++) {
                for(int d2 = 0; d2 < D; d2++) {
                    surfDivU += Psurf[d1, d2] * GradU[d2, d1];
                }
            }
            for(int dd = 0; dd < D; dd++) {
                acc += 1 / (We) * (Lambda(cpv.Xglobal) - Mu(cpv.Xglobal)) * surfDivU * Psurf[d, dd] * GradV[dd];
            }

            // surface shear part
            double[,] surfGradU  = new double[D,D];
            for(int d1 = 0; d1 < D; d1++) {
                for(int d2 = 0; d2 < D; d2++) {
                    for(int dd = 0; dd < D; dd++) {
                        surfGradU[d1, d2] += Psurf[d1, dd] * GradU[dd, d2];
                    }
                }
            }
            for(int d1 = 0; d1 < D; d1++) {
                for(int d2 = 0; d2 < D; d2++) {
                    for(int dd = 0; dd < D; dd++) {
                        acc += 1 / (We) * Mu(cpv.Xglobal) * surfGradU[d,d2] * Psurf[d2, d1] * GradV[d1];
                    }
                }
            }

            return acc;
        }


        public override double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {

            double[] EdgeNormal = inp.Normal;
            double[] SurfaceNormal_IN = SurfaceNormal(inp.Parameters_IN);
            double[] SurfaceNormal_OUT = SurfaceNormal(inp.Parameters_OUT);

            double[] Tangente_IN = Tangent(SurfaceNormal_IN, EdgeNormal);
            double[] Tangente_OUT = Tangent(SurfaceNormal_OUT, EdgeNormal);

            double[,] Psurf_IN = SurfaceProjection(SurfaceNormal_IN);
            double[,] Psurf_OUT = SurfaceProjection(SurfaceNormal_OUT);

            // isotropic
            double acc = -0.5 * 1 / (We) * Sigma(inp.X) * (Tangente_IN[d] + Tangente_OUT[d]) * (_vA - _vB);

            // dilatational
            double surfDivU_IN = 0.0;
            double surfDivU_OUT = 0.0;
            for(int d1 = 0; d1 < inp.D; d1++) {
                for(int dd = 0; dd < inp.D; dd++) {
                    surfDivU_IN += Psurf_IN[d1, dd] * _Grad_uA[dd, d1];
                    surfDivU_OUT += Psurf_OUT[d1, dd] * _Grad_uB[dd, d1];
                }
            }
            acc += -0.5 * 1 / (We) * (Lambda(inp.X) - Mu(inp.X)) * (surfDivU_IN * Tangente_IN[d] + surfDivU_OUT * Tangente_OUT[d]) * (_vA - _vB);

            // shear
            for(int d1 = 0; d1 < D; d1++) {
                for(int d2 = 0; d2 < D; d2++) {
                    acc += -0.5 * 1 / (We) * Mu(inp.X) * (Psurf_IN[d, d2] * _Grad_uA[d2, d1] * Tangente_IN[d1] + Psurf_OUT[d, d2] * _Grad_uB[d2, d1] * Tangente_OUT[d1]) * (_vA - _vB);
                }
            }

            return acc;
        }


        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {

            double Flx_InCell = 0;

            FreeXNSE_BcType edgType = BcMap.EdgeTag2Type[inp.EdgeTag];

            switch(edgType) {
                case FreeXNSE_BcType.Dirichlet:
                case FreeXNSE_BcType.Neumann: {

                    double[] EdgeNormal = inp.Normal;
                    double[] SurfaceNormal_IN = SurfaceNormal(inp.Parameters_IN);
                    double[] Tangente_IN = Tangent(SurfaceNormal_IN, EdgeNormal);
                    double[,] Psurf_IN = SurfaceProjection(SurfaceNormal_IN);

                    // isotropic
                    Flx_InCell += -1 / (We) * Sigma(inp.X) * Tangente_IN[d];

                    // dilatational
                    double surfDivU_IN = 0.0;
                    for(int d1 = 0; d1 < inp.D; d1++) {
                        for(int dd = 0; dd < inp.D; dd++) {
                            surfDivU_IN += Psurf_IN[d1, dd] * _Grad_uA[dd, d1];
                        }
                    }
                    Flx_InCell += - 1 / (We) * (Lambda(inp.X) - Mu(inp.X)) * surfDivU_IN * Tangente_IN[d] * (_vA - 0);

                    // shear
                    for(int d1 = 0; d1 < D; d1++) {
                        for(int d2 = 0; d2 < D; d2++) {
                            Flx_InCell += - 1 / (We) * Mu(inp.X) * Psurf_IN[d, d2] * _Grad_uA[d2, d1] * Tangente_IN[d1] * (_vA - 0);
                        }
                    }

                    break;
                }
                case FreeXNSE_BcType.Robin: {

                    double[] EdgeNormal = inp.Normal;
                    double[] SurfaceNormal_IN = SurfaceNormal(inp.Parameters_IN);
                    double[] Tangente_IN = Tangent(SurfaceNormal_IN, EdgeNormal);
                    double[,] Psurf_IN = SurfaceProjection(SurfaceNormal_IN);

                    double[] PSnI = new double[D]; // projection of surface/level-set normal onto domain boundary tangent
                    for(int d1 = 0; d1 < D; d1++) {
                        for(int d2 = 0; d2 < D; d2++) {
                            double nn = EdgeNormal[d1] * EdgeNormal[d2];
                            if(d1 == d2) {
                                PSnI[d1] += (1 - nn) * SurfaceNormal_IN[d2];
                            } else {
                                PSnI[d1] += -nn * SurfaceNormal_IN[d2];
                            }
                        }
                    }
                    double PSnINorm = PSnI.L2Norm();
                    double[] PSnINormal_IN = PSnI.Normalize(); // line normal: tangential to domain boundary & normal on contact line


                    // isotropic surface tension terms
                    for(int dd = 0; dd < D; dd++) {
                        Flx_InCell -= 1 / (We) * Sigma(inp.X) * (EdgeNormal[dd] * Tangente_IN[dd]) * EdgeNormal[d];
                    }

                    // for now only the isotropic part will take on the dynamic contact angle condition
                    // dilatational
                    double surfDivU_IN = 0.0;
                    for(int d1 = 0; d1 < inp.D; d1++) {
                        for(int dd = 0; dd < inp.D; dd++) {
                            surfDivU_IN += Psurf_IN[d1, dd] * _Grad_uA[dd, d1];
                        }
                    }
                    Flx_InCell += -1 / (We) * (Lambda(inp.X) - Mu(inp.X)) * surfDivU_IN * Tangente_IN[d] * (_vA - 0);

                    // shear
                    for(int d1 = 0; d1 < D; d1++) {
                        for(int d2 = 0; d2 < D; d2++) {
                            Flx_InCell += -1 / (We) * Mu(inp.X) * Psurf_IN[d, d2] * _Grad_uA[d2, d1] * Tangente_IN[d1] * (_vA - 0);
                        }
                    }


                    if(edgType == FreeXNSE_BcType.Robin) {

                        double alpha = Alpha(inp.X);
                        double theta0 = Theta(inp.X);

                        if(alpha < 0 || alpha == double.PositiveInfinity) {
                            // Free contactangle condition
                            for(int dd = 0; dd < D; dd++) {
                                Flx_InCell -= 1 / (We) * Sigma(inp.X) * (PSnINormal_IN[dd] * Tangente_IN[dd]) * PSnINormal_IN[d];
                            }
                        } else {
                            // Slipcondition
                            // Young's relation (static contact angle)
                            Flx_InCell -= 1 / (We) * Sigma(inp.X) * Math.Cos(theta0) * PSnINormal_IN[d];

                            // dissipative contact line force
                            // beta*(u*nL)
                            double g_D = BcMap.bndFunction[ArgumentOrdering[d] + "#" + ValidSpecies][inp.EdgeTag](inp.X, inp.time);

                            for(int d = 0; d < D; d++) {
                                Flx_InCell += alpha * ((_uA[d] - g_D) * PSnINormal_IN[d]) * PSnINormal_IN[d];
                            }
                        }
                    }

                    break;
                }
                default:
                break;
            }

            return Flx_InCell * _vA;
        }


        protected static double[] SurfaceNormal(double[] param) {

            double[] N = new double[param.Length];

            for(int d = 0; d < param.Length; d++) {
                N[d] = param[d];
            }

            return N.Normalize();
        }

        protected static double[,] SurfaceProjection(double[] Nsurf) {

            int D = Nsurf.Length;
            double[,] P = new double[D, D];

            for(int d = 0; d < D; d++) {
                for(int dd = 0; dd < D; dd++) {
                    if(dd == d)
                        P[d, dd] = (1 - Nsurf[d] * Nsurf[dd]);
                    else
                        P[d, dd] = (0 - Nsurf[d] * Nsurf[dd]);
                }
            }

            return P;
        }

        protected static double[] Tangent(double[] Nsurf, double[] Nedge) {
            Debug.Assert(Nsurf.Length == Nedge.Length);

            int D = Nsurf.Length;

            double[] tau = new double[D];
            for(int d1 = 0; d1 < D; d1++) {
                for(int d2 = 0; d2 < D; d2++) {
                    double nn = Nsurf[d1] * Nsurf[d2];
                    if(d1 == d2) {
                        tau[d1] += (1 - nn) * Nedge[d2];
                    } else {
                        tau[d1] += -nn * Nedge[d2];
                    }
                }
            }

            return tau.Normalize();
        }

        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this }; // only parameter dependent, not present in jacobian, boundary term if slip is on is dependent on velocity
        }

        Func<double[], double> Theta;
        Func<double[], double> Alpha;
        Func<double[], double> Sigma;
        Func<double[], double> Lambda;
        Func<double[], double> Mu;

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            if(cs.UserDefinedValues.Keys.Contains(Coefficientnames.Contactangle)) {
                Theta = (Func<double[], double>)cs.UserDefinedValues[Coefficientnames.Contactangle];
            }

            if(cs.UserDefinedValues.Keys.Contains(Coefficientnames.Contactlinefriction)) {
                Alpha = (Func<double[], double>)cs.UserDefinedValues[Coefficientnames.Contactlinefriction];
            }

            if(cs.UserDefinedValues.Keys.Contains(Coefficientnames.Surfacetensionfield)) {
                Sigma = (Func<double[], double>)cs.UserDefinedValues[Coefficientnames.Surfacetensionfield];
            }

            if(cs.UserDefinedValues.Keys.Contains(Coefficientnames.Surfaceshearviscosityfield)) {
                Mu = (Func<double[], double>)cs.UserDefinedValues[Coefficientnames.Surfaceshearviscosityfield];
            }

            if(cs.UserDefinedValues.Keys.Contains(Coefficientnames.Surfacedilatationalviscosityfield)) {
                Lambda = (Func<double[], double>)cs.UserDefinedValues[Coefficientnames.Surfacedilatationalviscosityfield];
            }
        }
    }

    public class BulkComponent_VolumeForce : BulkComponent_FreeXNSE {

        public override IList<string> ArgumentOrdering => new string[0];

        public override IList<string> ParameterOrdering => new string[0];

        public override TermActivationFlags VolTerms => TermActivationFlags.V;

        public override TermActivationFlags BoundaryEdgeTerms => TermActivationFlags.None;

        public override TermActivationFlags InnerEdgeTerms => TermActivationFlags.None;

        int d;
        Forcing VolForce;
        double Fr;

        public BulkComponent_VolumeForce(int d, double Fr, Forcing VolForce, int D, string[] spc, FreeXNSE_BoundaryCondMap BcMap) : base(D, spc, BcMap) {
            this.d = d;
            this.VolForce = VolForce;
            this.Fr = Fr;
        }

        public override double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            return -1.0/Fr * VolForce.Evaluate(cpv.Xglobal, cpv.time)[d] * V;
        }
        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return 0.0;
        }

        public override double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            return 0.0;
        }

        public override IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }

    public class BulkComponent_Dummy : IVolumeForm, ISpeciesFilter, ISupportsJacobianComponent {
        public TermActivationFlags VolTerms => TermActivationFlags.UxV;

        public IList<string> ArgumentOrdering => new string[] { dom };

        public IList<string> ParameterOrdering => new string[0];

        public string ValidSpecies => spc;

        string spc;
        string dom;
        public BulkComponent_Dummy(string spc, string dom) {
            this.spc = spc;
            this.dom = dom;
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            return U[0] * V;
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }
}
