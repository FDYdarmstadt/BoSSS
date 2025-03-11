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

using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;

namespace BoSSS.Solution.XNSECommon.Operator.Continuity {

    /// <summary>
    /// velocity jump penalty for the divergence operator, on edges
    /// </summary>
    public class DivergenceInSpeciesBulk_Edge : Divergence_DerivativeSource_Flux, ISpeciesFilter {

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="_component">
        /// component of the divergence
        /// </param>
        /// <param name="_bcmap"></param>
        public DivergenceInSpeciesBulk_Edge(int _component, IncompressibleMultiphaseBoundaryCondMap _bcmap, string spcName,  
            double _rho, double _vorZeichen, bool _RescaleConti)
            : base(_component, _bcmap) {

            rho = _rho;
            //m_spcId = spcId;
            ValidSpecies = spcName;

            this.RescaleConti = _RescaleConti;
            scale = _vorZeichen / ((RescaleConti) ? rho : 1.0);

            int d = base.component;
            base.bndFunction = _bcmap.bndFunction[VariableNames.Velocity_d(d) + "#" + spcName];
            
            this.m_bcmap = _bcmap;
        }

        IncompressibleMultiphaseBoundaryCondMap m_bcmap;

        double rho;

        bool RescaleConti;
        double scale;

        //SpeciesId m_spcId;

        public string ValidSpecies {
            get;
            private set;
        }



        protected override void InnerEdgeFlux(ref Foundation.CommonParams inp, double[] Uin, double[] Uout, out double FluxInCell, out double FluxOutCell) {
            base.InnerEdgeFlux(ref inp, Uin, Uout, out FluxInCell, out FluxOutCell);
            FluxInCell *= scale;
            FluxOutCell *= scale;
        }

        protected override void BorderEdgeFlux_(ref BoSSS.Foundation.CommonParamsBnd inp, double[] Uin, out double FluxInCell) {
            Debug.Assert(Uin.Length == 1);
            Debug.Assert(base.ArgumentOrdering.Count == 1);

            base.BorderEdgeFlux_(ref inp, Uin, out FluxInCell);

            FluxInCell *= scale;
        }



    }

    /// <summary>
    /// volume term for the Divergence / Continuity equation
    /// </summary>
    public class DivergenceInSpeciesBulk_Volume : Divergence_DerivativeSource, ISpeciesFilter {

        public DivergenceInSpeciesBulk_Volume(int _component, int _D, string spcName, double _rho, double _vorZeichen, bool _RescaleConti)
            : base(_component, _D) {

            ValidSpecies = spcName;
            scale = _vorZeichen / ((_RescaleConti) ? _rho : 1.0);
        }

        double scale;

        
        public string ValidSpecies {
            get;
            private set;
        }


        public override double _DerivativeSource(ilPSP.Vector x, double[] Parameters, double[,] GradientU) {
            return base._DerivativeSource(x, Parameters, GradientU) * scale;
        }


    }

    /// <summary>
    /// Variable density bulk term. Formulation suitable for use of the Newton solver
    /// </summary>
    public class DivergenceInSpeciesBulk_CentralDifferenceNewton : Divergence_CentralDifferenceJacobian, ISpeciesFilter {   

        public DivergenceInSpeciesBulk_CentralDifferenceNewton(string spcName, int Component, IncompressibleBoundaryCondMap Bcmap, int SpatDim, MaterialLaw EoS, int NumberOfChemicalSpecies) 
            : base(Component, Bcmap, SpatDim, EoS, NumberOfChemicalSpecies) {
            ValidSpecies = spcName;
        }
        /// <summary>
        /// constructor for incompressible
        /// </summary>
        /// <param name="spcName"></param>
        /// <param name="Component"></param>
        /// <param name="Bcmap"></param>
        /// <param name="SpatDim"></param>
        /// <param name="EoS"></param>
        /// <param name="NumberOfChemicalSpecies"></param>
        public DivergenceInSpeciesBulk_CentralDifferenceNewton(string spcName, int Component, IncompressibleBoundaryCondMap Bcmap, int SpatDim) : base(Component, Bcmap, SpatDim) {
            ValidSpecies = spcName;
        }

        public string ValidSpecies {
            get;
            private set;
        }
    }


    /// <summary>
    /// Variable density bulk term. Formulation suitable for use of the Newton solver
    /// </summary>
    public class DivergenceInSpeciesBulk_CentralDifference : Divergence_CentralDifference, ISpeciesFilter {

        public DivergenceInSpeciesBulk_CentralDifference(string spcName, int Component, IncompressibleBoundaryCondMap Bcmap, int SpatDim, MaterialLaw EoS, int NumberOfChemicalSpecies)
            : base(Component, Bcmap, EoS, NumberOfChemicalSpecies) {
            ValidSpecies = spcName;
        }
        /// <summary>
        /// constructor for incompressible
        /// </summary>
        /// <param name="spcName"></param>
        /// <param name="Component"></param>
        /// <param name="Bcmap"></param>
        /// <param name="SpatDim"></param>
        /// <param name="EoS"></param>
        /// <param name="NumberOfChemicalSpecies"></param>
        public DivergenceInSpeciesBulk_CentralDifference(string spcName, int Component, IncompressibleBoundaryCondMap Bcmap, int SpatDim) : base(Component, Bcmap) {
            ValidSpecies = spcName;
        }

        public string ValidSpecies {
            get;
            private set;
        }
    }



    /// <summary>
    /// velocity jump penalty for the divergence operator, on the level set
    /// </summary>
    public class DivergenceAtLevelSetLowMach : ILevelSetForm, ISupportsJacobianComponent {

        LevelSetTracker m_lsTrk;

        public DivergenceAtLevelSetLowMach(int _D, LevelSetTracker lsTrk, double _rhoA, double _rhoB,
            bool _MaterialInterface, double vorZeichen, bool RescaleConti,
            double _wA = 1.0, double _wB = 1.0) {
            this.D = _D;
            this.rhoA = _rhoA;
            this.rhoB = _rhoB;
            this.m_lsTrk = lsTrk;
            MaterialInterface = _MaterialInterface;

            scaleA = vorZeichen;
            scaleB = vorZeichen;
            //RescaleConti = false;
            if (RescaleConti) {
                scaleA *= rhoA;
                scaleB *= rhoB;
            }

            this.wA = _wA;
            this.wB = _wB;
        }

        bool MaterialInterface;
        int D;
        double rhoA;
        double rhoB;

        double scaleA;
        double scaleB;

        double wA;
        double wB;

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }

        public double InnerEdgeForm(ref CommonParams cp,
    double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB,
    double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            
            double uAxN = GenericBlas.InnerProd(U_Neg, cp.Normal);
            double uBxN = GenericBlas.InnerProd(U_Pos, cp.Normal);

            //// transform from species B to A: we call this the "A-fictitious" value
            //double uAxN_fict = uBxN;
            //// transform from species A to B: we call this the "B-fictitious" value
            //double uBxN_fict = uAxN;

            //// compute the fluxes: note that for the continuity equation, we use not a real flux,
            //// but some kind of penalization, therefore the fluxes have opposite signs!

            //double FlxNeg = -Flux(uAxN, uAxN_fict, 1.0, 1.0); // flux on A-side
            //double FlxPos = -Flux(uBxN_fict , uBxN , 1.0, 1.0);  // flux on B-side

            //FlxNeg *= scaleA;
            //FlxPos *= scaleB;





            double res = +0.5 * (uAxN + uBxN) * (this.rhoA * vA - this.rhoB * vB); // This is correct WITHOUT mass evaporation flux
            //double res = 0.5 * (this.rhoA * uAxN + this.rhoB * uBxN) * (vA - vB); // This is correct WITH mass evaporation flux

            //0.5 * (densityIn * Uin[Component] + densityOut * Uout[Component]) * inp.Normal[Component]
            return res;
            //return FlxNeg * vA - FlxPos * vB;
        }

        /// <summary>
        /// the penalty flux
        /// </summary>
        static double Flux(double UxN_in, double UxN_out, double w_in, double w_out) {
            return (UxN_in + UxN_out) * w_in / (w_in + w_out);
        }

        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

        public IList<string> ArgumentOrdering {
            get {
                return VariableNames.VelocityVector(this.D);
            }
        }

        public IList<string> ParameterOrdering {
            get { return null; }
        }


        public int LevelSetIndex {
            get { return 0; }
        }

        //public SpeciesId PositiveSpecies {
        //    get { return this.m_lsTrk.GetSpeciesId("B"); }
        //}

        //public SpeciesId NegativeSpecies {
        //    get { return this.m_lsTrk.GetSpeciesId("A"); }
        //}
        public string PositiveSpecies
        {
            get { return "B"; }
        }

        public string NegativeSpecies
        {
            get { return "A"; }
        }


    }


}
