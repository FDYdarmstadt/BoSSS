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
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;
using BoSSS.Platform;
using System.Diagnostics;
using BoSSS.Foundation;

namespace BoSSS.Solution.XNSECommon.Operator.Convection {

    class ConvectionAtLevelSet_LLF : ILevelSetForm {

        LevelSetTracker m_LsTrk;

        bool movingmesh;

        public ConvectionAtLevelSet_LLF(int _d, int _D, LevelSetTracker LsTrk, double _rhoA, double _rhoB, double _LFFA, double _LFFB, bool _MaterialInterface, IncompressibleMultiphaseBoundaryCondMap _bcmap, bool _movingmesh) {
            m_D = _D;
            m_d = _d;
            rhoA = _rhoA;
            rhoB = _rhoB;
            m_LsTrk = LsTrk;
            //varMode = _varMode;
            MaterialInterface = _MaterialInterface;
            movingmesh = _movingmesh;

            NegFlux = new ConvectionInBulk_LLF(_D, _bcmap, _d, _rhoA, _rhoB, _LFFA, double.NaN, LsTrk);
            NegFlux.SetParameter("A", LsTrk.GetSpeciesId("A"), null);
            PosFlux = new ConvectionInBulk_LLF(_D, _bcmap, _d, _rhoA, _rhoB, double.NaN, _LFFB, LsTrk);
            PosFlux.SetParameter("B", LsTrk.GetSpeciesId("B"), null);

        }

        bool MaterialInterface;
        double rhoA;
        double rhoB;
        int m_D;
        int m_d;
        //EquationAndVarMode varMode;

        // Use Fluxes as in Bulk Convection
        ConvectionInBulk_LLF NegFlux;
        ConvectionInBulk_LLF PosFlux;
        
        /*
        // Flux over interface
        public override void DerivativVar_LevelSetFlux(out double FlxNeg, out double FlxPos, 
            ref CommonParamsLs  cp,
            double[] U_Neg, double[] U_Pos, double[,] GradU_Neg, double[,] GradU_Pos) {

            double[] U_NegFict, U_PosFict;
            this.TransformU(ref U_Neg, ref U_Pos, out U_NegFict, out U_PosFict);

            double[] ParamsNeg = cp.ParamsNeg;
            double[] ParamsPos = cp.ParamsPos;
            double[] ParamsPosFict, ParamsNegFict;
            this.TransformU(ref ParamsNeg, ref ParamsPos, out ParamsNegFict, out ParamsPosFict);
            //Flux for negativ side
            {
                //double flx = 0.0;
                //for (int d = m_D - 1; d >= 0; d--)
                //    flx += cp.ParamsNeg[d] * cp.n[d];
                //flx *= U_Neg[0];
                //FlxNeg = flx;
                                
                BoSSS.Foundation.CommonParams inp; // = default(BoSSS.Foundation.InParams);
                inp.Parameters_IN = ParamsNeg;
                inp.Parameters_OUT = ParamsNegFict;
                inp.Normale = cp.n;
                inp.iEdge = int.MinValue;
                inp.GridDat = base.m_LsTrk.GridDat;
                inp.X = cp.x;
                inp.time = cp.time;
                //inp.jCellIn = cp.jCell;
                //inp.jCellOut = cp.jCell;


                FlxNeg = this.NegFlux.IEF(ref inp, U_Neg, U_NegFict);
            }
            // Flux for positive side
            {
                //double flx = 0.0;
                //for (int d = m_D - 1; d >= 0; d--)
                //    flx += cp.ParamsPos[d] * cp.n[d];
                //flx *= U_Pos[0];
                //FlxPos = flx;

                BoSSS.Foundation.CommonParams inp; // = default(BoSSS.Foundation.InParams);
                inp.Parameters_IN = ParamsPosFict;
                inp.Parameters_OUT = ParamsPos;
                inp.Normale = cp.n;
                inp.iEdge = int.MinValue;
                inp.GridDat = base.m_LsTrk.GridDat;
                inp.X = cp.x;
                inp.time = cp.time;
                //inp.jCellIn = cp.jCell;
                //inp.jCellOut = cp.jCell;

                FlxPos = this.PosFlux.IEF(ref inp, U_PosFict, U_Pos);
            }
        }
        */

        void TransformU(ref double[] U_Neg, ref double[] U_Pos, out double[] U_NegFict, out double[] U_PosFict) {
            if (this.MaterialInterface) {
                //if (varMode == EquationAndVarMode.mom_p)
                //    throw new NotImplementedException();

                U_NegFict = U_Pos;
                U_PosFict = U_Neg;

            } else {
                throw new NotImplementedException();
            }
        }


        /*
        public override void PrimalVar_LevelSetFlux(out double FlxNeg, out double FlxPos,
            ref CommonParamsLs cp, 
            double[] U_Neg, double[] U_Pos) {
            FlxNeg = 0;
            FlxPos = 0;
        }

        public override void FluxPotential(out double G, double[] U) {
            G = 0;
        }

        public override void Nu(out double NuNeg, out double NuPos,
            ref CommonParamsLs cp) {
            NuPos = 1.0;
            NuNeg = 1.0;
        }
        */

        public double LevelSetForm(ref CommonParamsLs cp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {
            double[] U_NegFict, U_PosFict;

            this.TransformU(ref U_Neg, ref U_Pos, out U_NegFict, out U_PosFict);

            double[] ParamsNeg = cp.ParamsNeg;
            double[] ParamsPos = cp.ParamsPos;
            double[] ParamsPosFict, ParamsNegFict;
            this.TransformU(ref ParamsNeg, ref ParamsPos, out ParamsNegFict, out ParamsPosFict);
            //Flux for negativ side
            double FlxNeg;
            {
                //double flx = 0.0;
                //for (int d = m_D - 1; d >= 0; d--)
                //    flx += cp.ParamsNeg[d] * cp.n[d];
                //flx *= U_Neg[0];
                //FlxNeg = flx;

                BoSSS.Foundation.CommonParams inp; // = default(BoSSS.Foundation.InParams);
                inp.Parameters_IN = ParamsNeg;
                inp.Parameters_OUT = ParamsNegFict;
                inp.Normale = cp.n;
                inp.iEdge = int.MinValue;
                inp.GridDat = this.m_LsTrk.GridDat;
                inp.X = cp.x;
                inp.time = cp.time;
                //inp.jCellIn = cp.jCell;
                //inp.jCellOut = cp.jCell;


                FlxNeg = this.NegFlux.IEF(ref inp, U_Neg, U_NegFict);
            }
            // Flux for positive side
            double FlxPos;
            {
                //double flx = 0.0;
                //for (int d = m_D - 1; d >= 0; d--)
                //    flx += cp.ParamsPos[d] * cp.n[d];
                //flx *= U_Pos[0];
                //FlxPos = flx;

                BoSSS.Foundation.CommonParams inp; // = default(BoSSS.Foundation.InParams);
                inp.Parameters_IN = ParamsPosFict;
                inp.Parameters_OUT = ParamsPos;
                inp.Normale = cp.n;
                inp.iEdge = int.MinValue;
                inp.GridDat = this.m_LsTrk.GridDat;
                inp.X = cp.x;
                inp.time = cp.time;
                //inp.jCellIn = cp.jCell;
                //inp.jCellOut = cp.jCell;

                FlxPos = this.PosFlux.IEF(ref inp, U_PosFict, U_Pos);
            }

            if (movingmesh)
                return 0.0;
            else 
                return FlxNeg * v_Neg - FlxPos * v_Pos;
        }

        public IList<string> ArgumentOrdering {
            get {
                return new string[] { VariableNames.Velocity_d(m_d) };
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return ArrayTools.Cat(VariableNames.Velocity0Vector(m_D), VariableNames.Velocity0MeanVector(m_D));  
            }
        }

        public int LevelSetIndex {
            get { return 0; }
        }

        public SpeciesId PositiveSpecies {
            get { return this.m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return this.m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }
    }
}
