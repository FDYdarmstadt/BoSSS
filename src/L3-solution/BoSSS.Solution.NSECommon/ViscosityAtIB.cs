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
using System.Linq;
using System.Text;
using BoSSS.Foundation.XDG;
using System.Diagnostics;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;
using BoSSS.Platform;
using ilPSP;
using BoSSS.Foundation;

namespace BoSSS.Solution.NSECommon.Operator.Viscosity {

    /// <summary>
    /// Viscosity at an immersed boundary;
    /// </summary>
    public class ViscosityAtIB : BoSSS.Foundation.XDG.ILevelSetForm, ISupportsJacobianComponent, ILevelSetEquationComponentCoefficient {

        //LevelSetTracker m_LsTrk;

        public ViscosityAtIB(int _d, int _D, double penalty_base, double _muA, int iLevSet, string FluidSpc, string SolidSpecies, bool UseLevelSetVelocityParameter) {

            this.m_penalty_base = penalty_base;
            //this.m_LsTrk = t;
            this.muA = _muA;
            this.component = _d;
            this.m_D = _D;
            this.m_iLevSet = iLevSet;
            this.m_SolidSpecies = SolidSpecies;
            this.m_FluidSpc = FluidSpc;
            this.m_UseLevelSetVelocityParameter = UseLevelSetVelocityParameter;
        }
        int m_iLevSet;
        string m_FluidSpc;
        string m_SolidSpecies;
        int component;
        int m_D;
        bool m_UseLevelSetVelocityParameter;

        /// <summary>
        /// Viskosity in species A
        /// </summary>
        double muA;

        /// <summary>
        /// safety factor
        /// </summary>
        double m_penalty_base;

        /// <summary>
        /// degree and spatial dimension
        /// </summary>
        double m_penalty_degree;


        /// <summary>
        /// length scale for negative species
        /// </summary>
        MultidimensionalArray NegLengthScaleS;

        /// <summary>
        /// <see cref="ILevelSetEquationComponentCoefficient.CoefficientUpdate"/>
        /// </summary>
        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {

            double _D = csA.GrdDat.SpatialDimension;
            double _p = DomainDGdeg.Max();

            double penalty_deg_tri = (_p + 1) * (_p + _D) / _D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes

            m_penalty_degree = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            NegLengthScaleS = csA.CellLengthScales;
        }


        /// <summary>
        /// computation of penalty parameter according to: $` \mathrm{SafetyFactor} \cdot k^2 / h `$
        /// </summary>
        protected double Penalty(int jCellIn) {

            double penaltySizeFactor_A = 1.0 / NegLengthScaleS[jCellIn];

            double penaltySizeFactor = penaltySizeFactor_A;

            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(m_penalty_degree));

            double scaledPenalty = penaltySizeFactor * m_penalty_degree * m_penalty_base;
            if(scaledPenalty.IsNaNorInf())
                throw new ArithmeticException("NaN/Inf detected for penalty parameter.");
            return scaledPenalty;

        }

        /// <summary>
        /// default-implementation
        /// </summary>
        public double InnerEdgeForm(ref CommonParams inp,
        //public override double EdgeForm(ref Linear2ndDerivativeCouplingFlux.CommonParams inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double[] N = inp.Normal;


            double _penalty = Penalty(inp.jCellIn);
            int D = N.Length;

            //var parameters_P = m_getParticleParams(inp.X, inp.time);
            //double[] uLevSet = new double[] { parameters_P[0], parameters_P[1] };
            //double wLevSet = parameters_P[2];
            //pRadius = parameters_P[3];


            Debug.Assert(this.ArgumentOrdering.Count == D);
            Debug.Assert(Grad_uA.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uB.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uA.GetLength(1) == D);
            Debug.Assert(Grad_uB.GetLength(1) == D);

            double Grad_uA_xN = 0, Grad_vA_xN = 0;
            for (int d = 0; d < D; d++) {
                Grad_uA_xN += Grad_uA[component, d] * N[d];
                Grad_vA_xN += Grad_vA[d] * N[d];
            }


            // Evaluate the complete velocity as a sum of translation and angular velocity
            double Ret = 0.0;
            
            double uAFict = 0;
            if(m_UseLevelSetVelocityParameter) {
                uAFict = inp.Parameters_IN[0];
            }

            Ret -= Grad_uA_xN * (vA);                           // consistency term
            Ret -= Grad_vA_xN * (uA[component] - uAFict);     // symmetry term

            Ret += _penalty * (uA[component] - uAFict) * (vA); // penalty term

            Debug.Assert(!(double.IsInfinity(Ret) || double.IsNaN(Ret)));

            return Ret * muA;
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            throw new NotSupportedException();
        }


        public int LevelSetIndex {
            get { return m_iLevSet; }
        }

        public IList<string> ArgumentOrdering {
            get { return VariableNames.VelocityVector(this.m_D); }
        }
    
        /// <summary>
        /// Species ID of the solid
        /// </summary>
        public string PositiveSpecies {
            get { return m_SolidSpecies; }
        }

        /// <summary>
        /// Species ID of the fluid; 
        /// </summary>
        public string NegativeSpecies {
            get { return m_FluidSpc; }
        }

        /// <summary>
        /// %
        /// </summary>
        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.V | TermActivationFlags.GradV;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                if(m_UseLevelSetVelocityParameter)
                    return new[] { VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(m_iLevSet), VariableNames.Velocity_d(component)) };
                else
                    return null;
            }
        }

        /// <summary>
        /// Linear component - returns this object itself.
        /// </summary>
        virtual public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

    }

    /// <summary>
    /// Viscosity at an immersed boundary;
    /// </summary>
    public class ViscosityAtIB_FullySymmetric : BoSSS.Foundation.XDG.ILevelSetForm, ISupportsJacobianComponent, ILevelSetEquationComponentCoefficient {

        public ViscosityAtIB_FullySymmetric(int _d, int _D, double penalty_base, double _muA, int iLevSet, string FluidSpc, string SolidSpecies, bool UseLevelSetVelocityParameter) {

            this.m_penalty_base = penalty_base;
            this.muA = _muA;
            this.component = _d;
            this.m_D = _D;
            this.m_iLevSet = iLevSet;
            this.m_SolidSpecies = SolidSpecies;
            this.m_FluidSpc = FluidSpc;
            this.m_UseLevelSetVelocityParameter = UseLevelSetVelocityParameter;

            this.NegativeSpecies = FluidSpc;
            this.PositiveSpecies = SolidSpecies;
        }
        int m_iLevSet;
        string m_FluidSpc;
        string m_SolidSpecies;
        int component;
        int m_D;
        bool m_UseLevelSetVelocityParameter;

        /// <summary>
        /// Viskosity in species A
        /// </summary>
        double muA;

        /// <summary>
        /// safety factor
        /// </summary>
        double m_penalty_base;

        /// <summary>
        /// degree and spatial dimension
        /// </summary>
        double m_penalty_degree;


        /// <summary>
        /// length scale for negative species
        /// </summary>
        MultidimensionalArray NegLengthScaleS;

        /// <summary>
        /// <see cref="ILevelSetEquationComponentCoefficient.CoefficientUpdate"/>
        /// </summary>
        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {

            double _D = csA.GrdDat.SpatialDimension;
            double _p = DomainDGdeg.Max();

            double penalty_deg_tri = (_p + 1) * (_p + _D) / _D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes

            m_penalty_degree = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            NegLengthScaleS = csA.CellLengthScales;
        }


        /// <summary>
        /// computation of penalty parameter according to: $` \mathrm{SafetyFactor} \cdot k^2 / h `$
        /// </summary>
        protected double Penalty(int jCellIn) {

            double penaltySizeFactor_A = 1.0 / NegLengthScaleS[jCellIn];

            double penaltySizeFactor = penaltySizeFactor_A;

            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(m_penalty_degree));

            double scaledPenalty = penaltySizeFactor * m_penalty_degree * m_penalty_base;
            if (scaledPenalty.IsNaNorInf())
                throw new ArithmeticException("NaN/Inf detected for penalty parameter.");
            return scaledPenalty;

        }

        /// <summary>
        /// default-implementation
        /// </summary>
        public double InnerEdgeForm(ref CommonParams inp,
        //public override double EdgeForm(ref Linear2ndDerivativeCouplingFlux.CommonParams inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double[] N = inp.Normal;


            double _penalty = Penalty(inp.jCellIn);

            int D = N.Length;

            //var parameters_P = m_getParticleParams(inp.X, inp.time);
            //double[] uLevSet = new double[] { parameters_P[0], parameters_P[1] };
            //double wLevSet = parameters_P[2];
            //pRadius = parameters_P[3];


            Debug.Assert(this.ArgumentOrdering.Count == D);
            Debug.Assert(Grad_uA.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uB.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uA.GetLength(1) == D);
            Debug.Assert(Grad_uB.GetLength(1) == D);

            double Grad_uA_xN = 0, Grad_vA_xN = 0;
            for (int d = 0; d < D; d++) {
                Grad_uA_xN += Grad_uA[component, d] * N[d];
                Grad_vA_xN += Grad_vA[d] * N[d];
            }


            // Evaluate the complete velocity as a sum of translation and angular velocity
            double Ret = 0.0;

            Vector uAFict = m_UseLevelSetVelocityParameter ? new Vector(inp.Parameters_IN, 0, D) : new Vector(D);

            Ret -= Grad_uA_xN * (vA);                           // consistency term
            Ret -= Grad_vA_xN * (uA[component] - uAFict[component]);     // symmetry term
            Ret += _penalty * (uA[component] - uAFict[component]) * (vA); // penalty term

            // transpose term
            for (int d = 0; d < D; d++) {
                Ret -= Grad_uA[d, component] * N[d] * (vA);                           // consistency term
                Ret -= Grad_vA[d] * (uA[d] - uAFict[d]) * N[component];     // symmetry term
            }

            Debug.Assert(!(double.IsInfinity(Ret) || double.IsNaN(Ret)));
            return Ret * muA;
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            throw new NotSupportedException();
        }


        public int LevelSetIndex {
            get { return m_iLevSet; }
        }

        public IList<string> ArgumentOrdering {
            get { return VariableNames.VelocityVector(this.m_D); }
        }

        /// <summary>
        /// Species of the solid
        /// </summary>
        public string PositiveSpecies {
            get;
            private set;
        }

        /// <summary>
        /// Species of the fluid; 
        /// </summary>
        public string NegativeSpecies {
            get;
            private set;
        }

        /// <summary>
        /// %
        /// </summary>
        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.V | TermActivationFlags.GradV;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                if (m_UseLevelSetVelocityParameter)
                    return VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(m_iLevSet), VariableNames.VelocityVector(m_D));
                else
                    return null;
            }
        }

        /// <summary>
        /// Linear component - returns this object itself.
        /// </summary>
        virtual public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

    }

    /// <summary>
    /// Slip Viscosity at an immersed boundary;
    /// </summary>
    public class ViscosityAtSlipIB : BoSSS.Foundation.XDG.ILevelSetForm, ISupportsJacobianComponent, ILevelSetEquationComponentCoefficient {

        //LevelSetTracker m_LsTrk;

        public ViscosityAtSlipIB(int _d, int _D, double penalty_base, double _muA, int iLevSet, string FluidSpc, string SolidSpecies, bool UseLevelSetVelocityParameter) {

            this.m_penalty_base = penalty_base;
            //this.m_LsTrk = t;
            this.muA = _muA;
            this.component = _d;
            this.m_D = _D;
            this.m_iLevSet = iLevSet;
            this.m_SolidSpecies = SolidSpecies;
            this.m_FluidSpc = FluidSpc;
            this.m_UseLevelSetVelocityParameter = UseLevelSetVelocityParameter;
        }
        int m_iLevSet;
        string m_FluidSpc;
        string m_SolidSpecies;
        int component;
        int m_D;
        bool m_UseLevelSetVelocityParameter;

        /// <summary>
        /// Viskosity in species A
        /// </summary>
        double muA;

        /// <summary>
        /// safety factor
        /// </summary>
        double m_penalty_base;

        /// <summary>
        /// degree and spatial dimension
        /// </summary>
        double m_penalty_degree;


        /// <summary>
        /// length scale for negative species
        /// </summary>
        MultidimensionalArray NegLengthScaleS;

        /// <summary>
        /// <see cref="ILevelSetEquationComponentCoefficient.CoefficientUpdate"/>
        /// </summary>
        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {

            double _D = csA.GrdDat.SpatialDimension;
            double _p = DomainDGdeg.Max();

            double penalty_deg_tri = (_p + 1) * (_p + _D) / _D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes

            m_penalty_degree = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            NegLengthScaleS = csA.CellLengthScales;
        }


        /// <summary>
        /// computation of penalty parameter according to: $` \mathrm{SafetyFactor} \cdot k^2 / h `$
        /// </summary>
        protected double Penalty(int jCellIn) {

            double penaltySizeFactor_A = 1.0 / NegLengthScaleS[jCellIn];

            double penaltySizeFactor = penaltySizeFactor_A;

            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(m_penalty_degree));

            double scaledPenalty = penaltySizeFactor * m_penalty_degree * m_penalty_base;
            if (scaledPenalty.IsNaNorInf())
                throw new ArithmeticException("NaN/Inf detected for penalty parameter.");
            return scaledPenalty;

        }

        /// <summary>
        /// default-implementation
        /// </summary>
        public double InnerEdgeForm(ref CommonParams inp,
        //public override double EdgeForm(ref Linear2ndDerivativeCouplingFlux.CommonParams inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double[] N = inp.Normal;


            double penalty = Penalty(inp.jCellIn);

            int D = N.Length;

            Debug.Assert(this.ArgumentOrdering.Count == D);
            Debug.Assert(Grad_uA.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uB.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uA.GetLength(1) == D);
            Debug.Assert(Grad_uB.GetLength(1) == D);

            double g_D = 0.0;
            double Acc = 0.0;

            for (int dN = 0; dN < D; dN++) {

                for (int dD = 0; dD < D; dD++) {
                    // consistency
                    Acc += muA * (inp.Normal[dN] * Grad_uA[dN, dD] * inp.Normal[dD]) * (vA * inp.Normal[component]);
                    // symmetry
                    Acc += muA * (inp.Normal[component] * Grad_vA[dD] * inp.Normal[dD]) * (uA[dN] - g_D) * inp.Normal[dN];
                }

                // penalty
                Acc -= muA * ((uA[dN] - g_D) * inp.Normal[dN]) * ((vA - 0) * inp.Normal[component]) * penalty;
                //Acc = 0;
            }

            return -Acc;
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            throw new NotSupportedException();
        }


        public int LevelSetIndex {
            get { return m_iLevSet; }
        }

        public IList<string> ArgumentOrdering {
            get { return VariableNames.VelocityVector(this.m_D); }
        }

        /// <summary>
        /// Species ID of the solid
        /// </summary>
        public string PositiveSpecies {
            get { return m_SolidSpecies; }
        }

        /// <summary>
        /// Species ID of the fluid; 
        /// </summary>
        public string NegativeSpecies {
            get { return m_FluidSpc; }
        }

        /// <summary>
        /// %
        /// </summary>
        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.V | TermActivationFlags.GradV;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                if (m_UseLevelSetVelocityParameter)
                    return new[] { VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(m_iLevSet), VariableNames.Velocity_d(component)) };
                else
                    return null;
            }
        }

        /// <summary>
        /// Linear component - returns this object itself.
        /// </summary>
        virtual public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

    }


    /// <summary>
    /// Viscosity at an immersed boundary;
    /// </summary>
    public class ViscosityAtIB_FullySymmetric_LowMach : BoSSS.Foundation.XDG.ILevelSetForm, ISupportsJacobianComponent, ILevelSetEquationComponentCoefficient {

        public ViscosityAtIB_FullySymmetric_LowMach(int _d, int _D, double penalty_base, ViscosityOption _ViscosityMode, PhysicsMode physmode , double constantViscosityValue, double reynolds, MaterialLaw EoS, int iLevSet, string FluidSpc, string SolidSpecies, bool UseLevelSetVelocityParameter) {

            this.m_penalty_base = penalty_base;
            this.muA = constantViscosityValue;
            this.component = _d;
            this.m_D = _D;
            this.m_iLevSet = iLevSet;
            this.m_SolidSpecies = SolidSpecies;
            this.m_FluidSpc = FluidSpc;
            this.m_UseLevelSetVelocityParameter = UseLevelSetVelocityParameter;
            this.m_ViscosityMode = _ViscosityMode;

            switch (_ViscosityMode) {
                case ViscosityOption.ConstantViscosity:
                    if (double.IsNaN(constantViscosityValue))
                        throw new ArgumentException("constantViscosityValue is missing!");
                    this.m_constantViscosityValue = constantViscosityValue;
                    break;
                case ViscosityOption.ConstantViscosityDimensionless:
                    if (double.IsNaN(reynolds))
                        throw new ArgumentException("Reynolds number is missing!");
                    this.m_reynolds = reynolds;
                    break;
                case ViscosityOption.VariableViscosity:
                    this.m_EoS = EoS;
                    break;
                case ViscosityOption.VariableViscosityDimensionless:
                    if (double.IsNaN(reynolds))
                        throw new ArgumentException("Reynolds number is missing!");
                    this.m_reynolds = reynolds;
                    this.m_EoS = EoS;
                    break;
                default:
                    throw new NotImplementedException();

            }

            switch (physmode) {
        
                case PhysicsMode.MixtureFraction:
                    m_ArgumentOrdering = ArrayTools.Cat(VariableNames.VelocityVector(m_D), VariableNames.MixtureFraction);
                    break;
                case PhysicsMode.Combustion:
                    m_ArgumentOrdering = ArrayTools.Cat(VariableNames.VelocityVector(m_D), VariableNames.Temperature); 
                    break;
                default:
                    throw new NotImplementedException();

            }

            this.NegativeSpecies = FluidSpc;
            this.PositiveSpecies = SolidSpecies;
        }


        int m_iLevSet;
        string m_FluidSpc;
        string m_SolidSpecies;
        int component;
        int m_D;
        bool m_UseLevelSetVelocityParameter;

        /// <summary>
        /// Viskosity in species A
        /// </summary>
        double muA;

        /// <summary>
        /// safety factor
        /// </summary>
        double m_penalty_base;

        /// <summary>
        /// degree and spatial dimension
        /// </summary>
        double m_penalty_degree;


        /// <summary>
        /// length scale for negative species
        /// </summary>
        MultidimensionalArray NegLengthScaleS;

        /// <summary>
        /// <see cref="ILevelSetEquationComponentCoefficient.CoefficientUpdate"/>
        /// </summary>
        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {

            double _D = csA.GrdDat.SpatialDimension;
            double _p = DomainDGdeg.Max();

            double penalty_deg_tri = (_p + 1) * (_p + _D) / _D; // formula for triangles/tetras
            double penalty_deg_sqr = (_p + 1.0) * (_p + 1.0); // formula for squares/cubes

            m_penalty_degree = Math.Max(penalty_deg_tri, penalty_deg_sqr); // the conservative choice

            NegLengthScaleS = csA.CellLengthScales;
        }


        /// <summary>
        /// computation of penalty parameter according to: $` \mathrm{SafetyFactor} \cdot k^2 / h `$
        /// </summary>
        protected double Penalty(int jCellIn) {

            double penaltySizeFactor_A = 1.0 / NegLengthScaleS[jCellIn];

            double penaltySizeFactor = penaltySizeFactor_A;

            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(m_penalty_degree));

            double scaledPenalty = penaltySizeFactor * m_penalty_degree * m_penalty_base;
            if (scaledPenalty.IsNaNorInf())
                throw new ArithmeticException("NaN/Inf detected for penalty parameter.");
            return scaledPenalty;

        }
        /// <summary>
        /// see <see cref="ViscosityOption"/>
        /// </summary>
        protected ViscosityOption m_ViscosityMode;

        /// <summary>
        /// in the case of constant viscosity, the value of the viscosity
        /// </summary>
        double m_constantViscosityValue = double.NaN;


        /// <summary>
        /// Reynolds number for dimensionless formulation.
        /// </summary>
        protected double m_reynolds = double.NaN;

        /// <summary>
        /// Optional material law for calculating the viscosity
        /// as a function of temperature / level-set.
        /// </summary>
        MaterialLaw m_EoS = null;

        /// <summary>
        /// the molecular viscosity
        /// </summary>
        virtual protected double Viscosity(double[] Parameters) {
            switch (m_ViscosityMode) {
                case ViscosityOption.ConstantViscosity:
                    return m_constantViscosityValue;
                case ViscosityOption.ConstantViscosityDimensionless:
                    return (1.0 / m_reynolds);
                case ViscosityOption.VariableViscosity:
                    if (m_EoS == null) {
                        return Parameters[0];
                    } else {
                        if (Parameters.Length > 1) {
                            return m_EoS.GetViscosity(Parameters);
                        } else {
                            return m_EoS.GetViscosity(Parameters[0]);
                        }
                    }
                case ViscosityOption.VariableViscosityDimensionless:
                    if (m_EoS == null) {
                        return (Parameters[0] / m_reynolds);
                    } else {
                        if (Parameters.Length > 1) {
                            return (m_EoS.GetViscosity(Parameters) / m_reynolds);
                        } else {
                            return (m_EoS.GetViscosity(Parameters[0]) / m_reynolds);
                        }

                    }
                default:
                    throw new NotImplementedException();
            }
        }


        /// <summary>
        /// default-implementation
        /// </summary>
        public double InnerEdgeForm(ref CommonParams inp,
        //public override double EdgeForm(ref Linear2ndDerivativeCouplingFlux.CommonParams inp,
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double[] N = inp.Normal;


            double _penalty = Penalty(inp.jCellIn);

            int D = N.Length;
            int NOargs = uA.Length;
            double[] viscArgsIN =  uA.GetSubVector(inp.D, NOargs - inp.D) ;
            double muA = this.Viscosity(viscArgsIN);
            //var parameters_P = m_getParticleParams(inp.X, inp.time);
            //double[] uLevSet = new double[] { parameters_P[0], parameters_P[1] };
            //double wLevSet = parameters_P[2];
            //pRadius = parameters_P[3];


            Debug.Assert(Grad_uA.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uB.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uA.GetLength(1) == D);
            Debug.Assert(Grad_uB.GetLength(1) == D);

            double Grad_uA_xN = 0, Grad_vA_xN = 0;
            for (int d = 0; d < D; d++) {
                Grad_uA_xN += Grad_uA[component, d] * N[d];
                Grad_vA_xN += Grad_vA[d] * N[d];
            }


            // Evaluate the complete velocity as a sum of translation and angular velocity
            double Ret = 0.0;

            Vector uAFict = m_UseLevelSetVelocityParameter ? new Vector(inp.Parameters_IN, 0, D) : new Vector(D);

            Ret -= Grad_uA_xN * (vA);                           // consistency term
            Ret -= Grad_vA_xN * (uA[component] - uAFict[component]);     // symmetry term
            Ret += _penalty * (uA[component] - uAFict[component]) * (vA); // penalty term

            // transpose term
            for (int d = 0; d < D; d++) {
                Ret -= Grad_uA[d, component] * N[d] * (vA);                           // consistency term
                Ret -= Grad_vA[d] * (uA[d] - uAFict[d]) * N[component];     // symmetry term
            }

            Debug.Assert(!(double.IsInfinity(Ret) || double.IsNaN(Ret)));
            return Ret * muA;
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            throw new NotSupportedException();
        }


        public int LevelSetIndex {
            get { return m_iLevSet; }
        }
        /// <summary>
        /// A list with all the arguments of this flux
        /// </summary>
        string[] m_ArgumentOrdering;
        public IList<string> ArgumentOrdering {
            get {
                return m_ArgumentOrdering;
            }
        }

        /// <summary>
        /// Species of the solid
        /// </summary>
        public string PositiveSpecies {
            get;
            private set;
        }

        /// <summary>
        /// Species of the fluid; 
        /// </summary>
        public string NegativeSpecies {
            get;
            private set;
        }

        /// <summary>
        /// %
        /// </summary>
        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.UxGradV | TermActivationFlags.GradUxV | TermActivationFlags.V | TermActivationFlags.GradV;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                if (m_UseLevelSetVelocityParameter)
                    return VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(m_iLevSet), VariableNames.VelocityVector(m_D));
                else
                    return null;
            }
        }

        /// <summary>
        /// Linear component - returns this object itself.
        /// </summary>
        virtual public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

    }
}