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
using BoSSS.Foundation.XDG;
using System.Diagnostics;
using BoSSS.Foundation;
using ilPSP;
using System.Linq;
using BoSSS.Solution.NSECommon;

namespace BoSSS.Application.XNSERO_Solver {
    public class ViscosityAtIB : BoSSS.Foundation.XDG.ILevelSetForm, ISupportsJacobianComponent, ILevelSetEquationComponentCoefficient {

        LevelSetTracker m_LsTrk;

        public ViscosityAtIB(int _d, int _D, LevelSetTracker t, double penalty_base, double _muA, int iLevSet, string FluidSpc, string SolidSpecies, bool UseLevelSetVelocityParameter, bool UsePhoretic) {

            this.m_penalty_base = penalty_base;
            this.m_LsTrk = t;
            this.FluidViscosity = _muA;
            Component = _d;
            this.m_D = _D;
            this.m_iLevSet = iLevSet;
            this.m_SolidSpecies = SolidSpecies;
            this.m_FluidSpc = FluidSpc;
            this.m_UseLevelSetVelocityParameter = UseLevelSetVelocityParameter;
            this.m_UsePhoretic = UsePhoretic;
        }
        readonly int m_iLevSet;
        readonly string m_FluidSpc;
        readonly string m_SolidSpecies;
        readonly int Component;
        readonly int m_D;
        readonly bool m_UseLevelSetVelocityParameter;
        readonly bool m_UsePhoretic;

        /// <summary>
        /// Viskosity in species A
        /// </summary>
        double FluidViscosity;

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

        private enum BoundaryConditionType {
            passive = 0,
            active = 1
        }

        static public bool write = true;

        /// <summary>
        /// default-implementation
        /// </summary>
        public double InnerEdgeForm(ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {

            Vector normalVector = inp.Normal;
            double _penalty = Penalty(inp.jCellIn);
            int dim = normalVector.Dim;

            

            Debug.Assert(uA.Length == ArgumentOrdering.Count);
            Debug.Assert(uB.Length == ArgumentOrdering.Count);
            Debug.Assert(Grad_uA.GetLength(0) == ArgumentOrdering.Count);
            Debug.Assert(Grad_uB.GetLength(0) == ArgumentOrdering.Count);
            Debug.Assert(Grad_uA.GetLength(1) == dim);
            Debug.Assert(Grad_uB.GetLength(1) == dim);
            if(inp.X.IsNullOrEmpty())
                throw new Exception("X is null or empty");
            if(inp.X.Abs() < 0)
                throw new ArithmeticException("invalid length of position vector");


            Vector uAFict = new Vector(inp.Parameters_IN[0], inp.Parameters_IN[1]);
            Vector activeStressVector = new Vector(inp.Parameters_IN[2], inp.Parameters_IN[3]);
            BoundaryConditionType bcType = activeStressVector.Abs() <= 1e-8 ? BoundaryConditionType.passive : BoundaryConditionType.active;

            
            if(m_UsePhoretic) {
                Debug.Assert(ArgumentOrdering.Count == dim + 1);
                double phoreticVal = uA[dim];
                Vector tangential = normalVector.Rotate2D(Math.PI * 0.5);

                // todo: add computation of slip velocity.
                // Note: if the relation is non-linear, special treatment is required!
                if(write) {
                    Console.WriteLine("todo: add computation of slip velocity; phoretic value is " + phoreticVal);
                    write = false; // prevent end-less output; 
                }
            } else {
                Debug.Assert(ArgumentOrdering.Count == dim);
            }


            // Gradient of u and v 
            // =====================
            double Grad_uA_xN = 0, Grad_vA_xN = 0;
            for(int d = 0; d < dim; d++) {
                Grad_uA_xN += Grad_uA[Component, d] * inp.Normal[d];
                Grad_vA_xN += Grad_vA[d] * inp.Normal[d];
            }

            double returnValue = 0.0;

            // 3D for IBM_Solver
            // =====================
            if(dim == 3) {
                returnValue -= Grad_uA_xN * (vA);                                                    // consistency term
                returnValue -= Grad_vA_xN * (uA[Component] - 0);                                     // symmetry term
                returnValue += _penalty * (uA[Component] - 0) * (vA);                           // penalty term
                Debug.Assert(!(double.IsInfinity(returnValue) || double.IsNaN(returnValue)));
                return returnValue * FluidViscosity;
            }

            // 2D
            // =====================
            switch(bcType) {
                case BoundaryConditionType.passive: {
                    for(int d = 0; d < dim; d++) {
                        returnValue -= FluidViscosity * Grad_uA[Component, d] * vA * normalVector[d];
                        returnValue -= FluidViscosity * Grad_vA[d] * (uA[Component] - uAFict[Component]) * normalVector[d];
                    }
                    returnValue += FluidViscosity * (uA[Component] - uAFict[Component]) * vA * _penalty;
                    break;
                }

                case BoundaryConditionType.active: {
                    // normal direction, solid wall
                    for(int dN = 0; dN < dim; dN++) {
                        for(int dD = 0; dD < dim; dD++) {
                            // consistency term
                            returnValue -= FluidViscosity * (normalVector[dN] * Grad_uA[dN, dD] * normalVector[dD]) * vA * normalVector[Component];
                            // symmetry term
                            returnValue -= FluidViscosity * (Grad_vA[dD] * normalVector[dD]) * (normalVector[dN] * uA[dN] - normalVector[dN] * uAFict[dN]) * normalVector[Component];
                        }
                        // penalty term
                        returnValue += FluidViscosity * inp.Normal[dN] * (uA[dN] - uAFict[dN]) * normalVector[Component] * vA * _penalty;
                    }
                    // tangential direction, active part
                    double[,] P = new double[dim, dim];
                    for(int d1 = 0; d1 < dim; d1++) {
                        for(int d2 = 0; d2 < dim; d2++) {
                            if(d1 == d2) {
                                P[d1, d2] = 1 - normalVector[d1] * normalVector[d2];
                            } else {
                                P[d1, d2] = normalVector[d1] * normalVector[d2];
                            }
                        }
                    }
                    for(int d1 = 0; d1 < dim; d1++) {
                        for(int d2 = 0; d2 < dim; d2++) {
                            returnValue -= P[d1, d2] * activeStressVector[d2] * (P[d1, Component] * vA);
                        }
                    }
                    break;
                }
            }
            Debug.Assert(!(double.IsInfinity(returnValue) || double.IsNaN(returnValue)));
            return returnValue;
        }

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            throw new NotSupportedException();
        }


        public int LevelSetIndex {
            get { return m_iLevSet; }
        }

        public IList<string> ArgumentOrdering {
            get { 
                var ret = VariableNames.VelocityVector(this.m_D); 
                if(this.m_UsePhoretic) {
                    ret = ret.Append(VariableNames.Phoretic).ToArray();
                }
                return ret;
            }
        }

        /// <summary>
        /// Species ID of the solid
        /// </summary>
        public SpeciesId PositiveSpecies {
            get { return m_LsTrk.GetSpeciesId(m_SolidSpecies); }
        }

        /// <summary>
        /// Species ID of the fluid; 
        /// </summary>
        public SpeciesId NegativeSpecies {
            get { return m_LsTrk.GetSpeciesId(m_FluidSpc); }
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
                return new[] { VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(m_iLevSet), VariableNames.Velocity_d(0)),
                               VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(m_iLevSet), VariableNames.Velocity_d(1)),
                               VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(m_iLevSet), VariableNames.SurfaceForceComponent(0)),
                               VariableNames.AsLevelSetVariable(VariableNames.LevelSetCGidx(m_iLevSet), VariableNames.SurfaceForceComponent(1))};
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