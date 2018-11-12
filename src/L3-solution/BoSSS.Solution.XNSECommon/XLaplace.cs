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
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using System.Diagnostics;
using ilPSP;
using BoSSS.Solution.Control;

namespace BoSSS.Solution.XNSECommon {

    public class XLaplaceBCs {

        /// <summary>
        /// Function which determines which part of the domain boundary is of Dirichlet type
        /// </summary>
        public Func<CommonParamsBnd, bool> IsDirichlet;

        /// <summary>
        /// Dirichlet boundary value
        /// </summary>
        public Func<CommonParamsBnd, double> g_Diri;

        /// <summary>
        /// Neumann boundary value
        /// </summary>
        public Func<CommonParamsBnd, double> g_Neum;
    }

    /// <summary>
    /// Laplace operator in the bulk.
    /// </summary>
    public class XLaplace_Bulk : BoSSS.Solution.NSECommon.SIPLaplace, IEquationComponentSpeciesNotification, IEquationComponentCoefficient {

        public XLaplace_Bulk(LevelSetTracker __LsTrk, double __penatly_baseFactor, string n, XLaplaceBCs boundaries, double sw, double _muA, double _muB, MultidimensionalArray PenaltyLengthScales, XLaplace_Interface.Mode _m)
            : base(__penatly_baseFactor, PenaltyLengthScales, n) {
            muA = _muA;
            muB = _muB;
            base.m_alpha = sw;
            this.boundaries = boundaries;
            this.LsTrk = __LsTrk;
            this.penatly_baseFactor = __penatly_baseFactor;
            this.m_Mode = _m;
        }

        double muA;
        double muB;
        SpeciesId SpcId;
        LevelSetTracker LsTrk;
        double penatly_baseFactor;
        XLaplace_Interface.Mode m_Mode;
        MultidimensionalArray m_LenScales;

        string current_species;

        public void SetParameter(string speciesName, SpeciesId __SpcId) {
            switch(speciesName) {
                case "A": species_Mu = muA; otherSpecies_Mu = muB; SpcId = __SpcId; break;
                case "B": species_Mu = muB; otherSpecies_Mu = muA; SpcId = __SpcId; break;
                default: throw new ArgumentException("Unknown species.");
            }
            current_species = speciesName;
        }

        double species_Mu;
        double otherSpecies_Mu;

        public override double Nu(double[] x, double[] parameters, int jCell) {
            return this.species_Mu;
        }

        XLaplaceBCs boundaries;

        protected override double g_Diri(ref CommonParamsBnd inp) {
            return boundaries.g_Diri(inp);
        }

        protected override bool IsDirichlet(ref CommonParamsBnd inp) {
            return boundaries.IsDirichlet(inp);
        }
        protected override double g_Neum(ref CommonParamsBnd inp) {
            return boundaries.g_Neum(inp);
        }

        double GetPenalty(ref CommonParams inp) {
            //double penaltySizeFactor_A = 1.0 / this.ccBB.Get_hminBB(this.SpcId, inp.jCellIn);
            //double penaltySizeFactor_B = 1.0 / this.ccBB.Get_hminBB(this.SpcId, inp.jCellOut);
            double penaltySizeFactor_A = 1.0 / this.m_LenScales[inp.jCellIn];
            double penaltySizeFactor_B = 1.0 / this.m_LenScales[inp.jCellOut];
            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
            double penaltySizeFactor = Math.Max(penaltySizeFactor_A, penaltySizeFactor_B);


            double penalty_muFactor;
            switch(this.m_Mode) {
                case XLaplace_Interface.Mode.SIP:
                //case XLaplace_Interface.Mode.VolumeScaled:
                case XLaplace_Interface.Mode.SWIP:
                penalty_muFactor = species_Mu;
                break;

                default:
                throw new NotImplementedException();
            }

            return this.penatly_baseFactor * penaltySizeFactor * penalty_muFactor;
        }
        

        override public double InnerEdgeForm(ref CommonParams inp, double[] _uA, double[] _uB, double[,] _Grad_uA, double[,] _Grad_uB, double _vA, double _vB, double[] _Grad_vA, double[] _Grad_vB) {
            double Acc = 0.0;


            double muA = this.Nu(inp.X, inp.Parameters_IN, inp.jCellIn);
            double muB = this.Nu(inp.X, inp.Parameters_OUT, inp.jCellOut);

            double scaleIN, scaleOT;
            ComputeScaling(ref inp, out scaleIN, out scaleOT);

            
            for (int d = 0; d < inp.D; d++) {
                Acc += (scaleIN * muA * _Grad_uA[0, d] + scaleOT * muB * _Grad_uB[0, d]) * (_vA - _vB) * inp.Normale[d];  // consistency term
                Acc += (scaleIN * muA * _Grad_vA[d] + scaleOT * muB * _Grad_vB[d]) * (_uA[0] - _uB[0]) * inp.Normale[d];  // symmetry term
            }
            Acc *= this.m_alpha;

            
            Acc -= (_uA[0] - _uB[0]) * (_vA - _vB) * this.GetPenalty(ref inp); // penalty term

            return Acc;
        }

        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
           
            return base.BoundaryEdgeForm(ref inp, _uA, _Grad_uA, _vA, _Grad_vA);
        }


        void ComputeScaling(ref CommonParams inp, out double scaleIN, out double scaleOT) {
            switch(this.m_Mode) {
                case XLaplace_Interface.Mode.SIP:
                case XLaplace_Interface.Mode.SWIP: {
                    scaleIN = 0.5;
                    scaleOT = 0.5;
                    return;
                }

                /*
                case XLaplace_Interface.Mode.VolumeScaled: {
                    double volIN = LsTrk._Regions.GetSpeciesVolume(inp.jCellIn, this.SpcId);
                    double volOT = LsTrk._Regions.GetSpeciesVolume(inp.jCellOut, this.SpcId);

                    scaleIN = volIN / (volIN + volOT);
                    scaleOT = volOT / (volIN + volOT);

                    Debug.Assert(Math.Abs(scaleIN + scaleOT - 1.0) <= 1.0e-8);
                    return;
                }
                */
                default:
                throw new NotImplementedException();
            }
        }

        public void CoefficientUpdate(CoefficientSet cs, int[] DomainDGdeg, int TestDGdeg) {
            this.m_LenScales = cs.CellLengthScales;
        }
    }

    /// <summary>
    /// Laplace operator at the interface
    /// </summary>
    public class XLaplace_Interface : ILevelSetForm, ILevelSetEquationComponentCoefficient {

        public enum Mode {
            SWIP,
            SIP
        }

        protected LevelSetTracker m_LsTrk; 

        public XLaplace_Interface(LevelSetTracker lstrk, double _muA, double _muB, double __penatly_baseFactor, Mode _m) {
            this.m_LsTrk = lstrk;
            this.muA = _muA;
            this.muB = _muB;
            this.penatly_baseFactor = __penatly_baseFactor;
            this.m_mode = _m;
        }


        protected double muA;
        protected double muB;
        protected double penatly_baseFactor;
        protected Mode m_mode;

        
        public virtual double LevelSetForm(ref CommonParamsLs inp, 
            double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB,
            double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double[] N = inp.n;
            double Grad_uA_xN = 0, Grad_uB_xN = 0, Grad_vA_xN = 0, Grad_vB_xN = 0;
            int D = N.Length;
            Debug.Assert(Grad_uA.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uB.GetLength(0) == this.ArgumentOrdering.Count);
            Debug.Assert(Grad_uA.GetLength(1) == D);
            Debug.Assert(Grad_uB.GetLength(1) == D);

            for(int d = 0; d < D; d++) {
                Grad_uA_xN += Grad_uA[0, d] * N[d];
                Grad_uB_xN += Grad_uB[0, d] * N[d];
                Grad_vA_xN += Grad_vA[d] * N[d];
                Grad_vB_xN += Grad_vB[d] * N[d];
            }

            double omega_A, omega_B;
            ComputeScaling(ref inp, out omega_A, out omega_B);
            
            double Ret = 0.0;
            Ret += (muA * omega_A * Grad_uA_xN + muB * omega_B * Grad_uB_xN) * (vA - vB);
            Ret += (muA * omega_A * Grad_vA_xN + muB * omega_B * Grad_vB_xN) * (uA[0] - uB[0]);

            //
            Ret -= GetPenalty(ref inp) * (uA[0] - uB[0]) * (vA - vB);

            return Ret;
        }

        protected double GetPenalty(ref CommonParamsLs inp) {
            //double penaltySizeFactor_A = 1.0 / this.ccBB.Get_hminBB(this.NegativeSpecies, inp.jCell);
            //double penaltySizeFactor_B = 1.0 / this.ccBB.Get_hminBB(this.PositiveSpecies, inp.jCell);

            double PosCellLengthScale = PosLengthScaleS[inp.jCell];
            double NegCellLengthScale = NegLengthScaleS[inp.jCell];

            double penaltySizeFactor_A = 1.0 / NegCellLengthScale;
            double penaltySizeFactor_B = 1.0 / PosCellLengthScale;
            Debug.Assert(!double.IsNaN(penaltySizeFactor_A));
            Debug.Assert(!double.IsNaN(penaltySizeFactor_B));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_A));
            Debug.Assert(!double.IsInfinity(penaltySizeFactor_B));
            double penaltySizeFactor = Math.Max(penaltySizeFactor_A, penaltySizeFactor_B);


            double penalty_muFactor;
            switch(this.m_mode) {
                case Mode.SWIP:
                penalty_muFactor = (2.0 * muA * muB) / (muA + muB);
                break;


                case Mode.SIP:
                //case Mode.VolumeScaled:
                penalty_muFactor = Math.Max(Math.Abs(muA), Math.Abs(muB)) * Math.Sign(muA);
                break;

                default:
                throw new NotImplementedException();
            }

            return this.penatly_baseFactor * penaltySizeFactor * penalty_muFactor;
        }

        //List<int> cellElo = new List<int>();


        protected void ComputeScaling(ref CommonParamsLs inp, out double scaleIN, out double scaleOT) {
            Debug.Assert(Math.Sign(muA) == Math.Sign(muB));
                    
            switch(this.m_mode) {
                case Mode.SWIP: {
                    scaleIN = muB / (muA + muB);
                    scaleOT = muA / (muA + muB);
                    return;
                }

                case Mode.SIP: {
                    // Konventionell:
                    scaleIN = 0.5;
                    scaleOT = 0.5;
                    return;
                }
                          /*
                case Mode.VolumeScaled: {
                    double volIN = this.m_LsTrk._Regions.GetSpeciesVolume(inp.jCell, this.NegativeSpecies);
                    double volOT = this.m_LsTrk._Regions.GetSpeciesVolume(inp.jCell, this.PositiveSpecies);

                    scaleIN = volIN / (volIN + volOT);
                    scaleOT = volOT / (volIN + volOT);
                    Debug.Assert(Math.Abs(scaleIN + scaleOT - 1.0) <= 1.0e-8);
                    return;
                }
                */

                default:
                throw new NotImplementedException();
            }
        }

        MultidimensionalArray PosLengthScaleS;
        MultidimensionalArray NegLengthScaleS;

        public void CoefficientUpdate(CoefficientSet csA, CoefficientSet csB, int[] DomainDGdeg, int TestDGdeg) {
            NegLengthScaleS = csA.CellLengthScales;
            PosLengthScaleS = csB.CellLengthScales;
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public IList<string> ArgumentOrdering {
            get { return new string[] { "u" }; }
        }

        public SpeciesId PositiveSpecies {
            get { return m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.GradUxV | TermActivationFlags.UxGradV;
            }
        }

        public IList<string> ParameterOrdering {
            get { return null; }
        }
    }


}
