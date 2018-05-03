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
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Utils;

namespace BoSSS.Solution.LevelSetTools.EllipticExtension {


    /// <summary>
    /// The Laplace Operator for the Elliptic Extension Velocity as per Paper
    /// Utz,Kummer - A high-order discontinuous Galerkin method for extension problems, IJNME 2017
    /// 
    ///with the addition of an isotropic viscosity term \epsilon
    ///\div ((\epsilon I + \nabla \phi \otimes \nabla \phi ) \cdot \nabla u ) = 0
    /// </summary>
    public abstract class EllipticExtVelForm : BoSSS.Foundation.IEdgeForm, BoSSS.Foundation.IVolumeForm, IObserver<LevelSetTracker.LevelSetRegions>{
        /// <summary>
        /// CTOR
        /// </summary>
        /// <param name="PenaltyBase">Penalty Mutliplier based on Polynomial Degree</param>
        /// <param name="IsotropicViscosity"></param>
        /// <param name="LSTrck">Level-Set Tracker</param>
        public EllipticExtVelForm(double PenaltyBase, double IsotropicViscosity , LevelSetTracker LSTrck) {
            this.PenaltyBase = PenaltyBase;
            this.IsotropicViscosity = IsotropicViscosity;
            this.LSTrck = LSTrck;
            this.D = LSTrck.GridDat.SpatialDimension;
            this.cj = LSTrck.GridDat.Cells.cj;
            LSTrck.Subscribe(this);
            CutCellMask = LSTrck.Regions.GetCutCellMask().GetBitMaskWithExternal();
        }

#region Observer for Level-Set Tracker
        protected BitArray CutCellMask;
        public void OnNext(LevelSetTracker.LevelSetRegions value)
        {
            CutCellMask = LSTrck.Regions.GetCutCellMask().GetBitMaskWithExternal();
        }

        public void OnError(Exception error)
        {
            // Do Nothing
        }

        public void OnCompleted()
        {
            // Do nothing
        }
#endregion

        /// <summary>
        /// Length Scales for Penalty Parameter
        /// </summary>
        protected MultidimensionalArray cj;

        /// <summary>
        /// A Small number to avoid division by zero
        /// </summary>
        protected static double myEps = Math.Sqrt(double.Epsilon);


        protected double IsotropicViscosity = 1e-2;


        /// <summary>
        /// Spatial Degree
        /// </summary>
        protected int D;
        /// <summary>
        /// Penalty Multiplier based on polynomial degrees
        /// </summary>
        protected double PenaltyBase;

        /// <summary>
        /// Terms, which are to be evaluated in Volume Integral
        /// </summary>
        public TermActivationFlags VolTerms {
            get {
                return (TermActivationFlags.GradUxGradV);
            }
        }

        /// <summary>
        /// The Volume Form, see <see cref="BoSSS.Foundation.IVolumeForm"/>
        /// </summary>
        /// <param name="cpv"></param>
        /// <param name="U"></param>
        /// <param name="GradU"></param>
        /// <param name="V"></param>
        /// <param name="GradV"></param>
        /// <returns></returns>
        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            Debug.Assert(GradU.GetLength(0) == 1);
            Debug.Assert(GradU.GetLength(1) == D);
            double Acc = 0;

            double GradPhiGradU = 0;
            double GradPhiGradV = 0;
            double AbsGradPhiSquared = 0;
            double StandardForm = 0;

            /// Weak Form
            for (int d = 0; d < D; d++) {
                GradPhiGradU += cpv.Parameters[d] * GradU[0, d];
                GradPhiGradV += cpv.Parameters[d] * GradV[d];
                AbsGradPhiSquared += cpv.Parameters[d] * cpv.Parameters[d];
                StandardForm += GradU[0, d] * GradV[d];
            }
            Acc += GradPhiGradU * GradPhiGradV / (AbsGradPhiSquared + myEps) ;
            Acc += IsotropicViscosity * StandardForm ;

            /// Strong Form, added due to Chessa et al. 2002
            //Acc -= GradPhiGradU / (Math.Sqrt(AbsGradPhiSquared) + myEps) * Math.Sign(cpv.Parameters[D]) * V;

            return Acc;
        }

        /// <summary>
        /// <see cref="IEquationComponent.ArgumentOrdering"/>
        /// </summary>
        public IList<string> ArgumentOrdering {
            get {
                return new string[] { "Extension" };
            }
        }

        /// <summary>
        /// <see cref="IEquationComponent.ParameterOrdering"/>
        /// </summary>
        public virtual IList<string> ParameterOrdering {
            get {
                return VariableNames.LevelSetGradient(D);
            }
        }

        /// <summary>
        /// <see cref="IEdgeForm.BoundaryEdgeTerms"/>
        /// </summary>
        public TermActivationFlags BoundaryEdgeTerms {
            get {
                return (TermActivationFlags.GradUxV | TermActivationFlags.UxV | TermActivationFlags.V);
            }
        }

        /// <summary>
        /// <see cref="IEdgeForm.InnerEdgeTerms"/>
        /// </summary>
        public TermActivationFlags InnerEdgeTerms {
            get {
                return (TermActivationFlags.UxV | TermActivationFlags.GradUxV | TermActivationFlags.UxGradV);
            }
        }

        /// <summary>
        /// The ubiquitous <see cref="LevelSetTracker"/>
        /// </summary>
        protected LevelSetTracker LSTrck;

        /// <summary>
        /// <see cref="swipViscosity_Term1"/>
        /// </summary>
        /// <param name="inp">
        /// </param>
        /// <param name="uA">
        /// ExtensionVelocity
        /// </param>
        /// <param name="uB"></param>
        /// <param name="Grad_uA"></param>
        /// <param name="Grad_uB"></param>
        /// <param name="vA"></param>
        /// <param name="vB"></param>
        /// <param name="Grad_vA"></param>
        /// <param name="Grad_vB"></param>
        /// <returns></returns>
        public abstract double InnerEdgeForm(ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB);


        /// <summary>
        /// The Idea here is, that the artificial viscosity term needs to have a free boundary, i.e.
        /// the flux at the boundary edge must be controlled by
        /// </summary>
        /// <param name="inp"></param>
        /// <param name="uA"></param>
        /// <param name="Grad_uA"></param>
        /// <param name="vA"></param>
        /// <param name="Grad_vA"></param>
        /// <returns></returns>
        public virtual double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] uA, double[,] Grad_uA, double vA, double[] Grad_vA) {
            double Acc = 0; // => return
            // Initialize all values, which are calculoated in the flux itself
            double GradUGradPhi = 0;
            double GradVGradPhi = 0;
            double nGradPhi = 0;
            double AbsGradPhi = 0;  // => divide by this to normalize GradPhi

            for (int d = 0; d < inp.D; d++)
            {
                GradUGradPhi += Grad_uA[0, d] * inp.Parameters_IN[d];
                GradVGradPhi += Grad_vA[d] * inp.Parameters_IN[d];
                nGradPhi += inp.Normale[d] * inp.Parameters_IN[d];
                AbsGradPhi += inp.Parameters_IN[d] * inp.Parameters_IN[d];
            }
            AbsGradPhi = Math.Sqrt(AbsGradPhi) + myEps;
            GradUGradPhi /= AbsGradPhi;
            GradVGradPhi /= AbsGradPhi;
            nGradPhi /= AbsGradPhi;


            // Central Differences for Isotropic Viscosity
            //============================================
            for (int d = 0; d < inp.D; d++){
                Acc += IsotropicViscosity *  (Grad_uA[0, d]) * (vA) * inp.Normale[d] ;
             //   Acc += IsotropicViscosity *  (Grad_vA[d]   ) * (uA - uOut[0]) * inp.Normale[d];
            }
            Acc -= IsotropicViscosity * GradUGradPhi * vA * nGradPhi;
            return -Acc;
        }


        protected double GetPenalty(int jCellIn, int jCellOut, MultidimensionalArray cj) {
            double cj_in = cj[jCellIn];
            double mu = PenaltyBase * cj_in;
            if (jCellOut >= 0) {
                double cj_out = cj[jCellOut];
                mu = Math.Max(mu, PenaltyBase * cj_out);
            }
            return mu;
        }

    }


    /// <summary>
    /// An SIP-Based flux for the Elliptic Extension Velocity,
    /// This one is based on a cenmtral difference idea and thus fully coupled
    /// </summary>
    public class EllipticExtVelFormSWIP: EllipticExtVelForm {

        /// <summary>
        /// CTOR
        /// </summary>
        /// <param name="PenaltyBase">penalty multiplier based on polynomials</param>
        /// <param name="LSTrck"><see cref="LevelSetTracker"/></param>
        public EllipticExtVelFormSWIP(double PenaltyBase, double IsotropicViscosity , LevelSetTracker LSTrck):base(PenaltyBase, IsotropicViscosity, LSTrck) {
        }

        /// <summary>
        /// <see cref="swipViscosity_Term1"/>
        /// </summary>
        /// <param name="inp">
        /// </param>
        /// <param name="uA">
        /// ExtensionVelocity
        /// </param>
        /// <param name="uB"></param>
        /// <param name="Grad_uA"></param>
        /// <param name="Grad_uB"></param>
        /// <param name="vA"></param>
        /// <param name="vB"></param>
        /// <param name="Grad_vA"></param>
        /// <param name="Grad_vB"></param>
        /// <returns></returns>
        public override double InnerEdgeForm(ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            double Acc = 0.0;
            double pnlty = GetPenalty(inp.jCellIn, inp.jCellOut, this.cj);

            double GradUAGradPhi = 0;
            double GradUBGradPhi = 0;
            double GradVAGradPhi = 0;
            double GradVBGradPhi = 0;
            double nGradPhi = 0;
            for (int d = 0; d < inp.D; d++) {
                GradUAGradPhi += Grad_uA[0, d] * inp.Parameters_IN[d];
                GradVAGradPhi += Grad_vA[d] * inp.Parameters_IN[d];
                GradVBGradPhi += Grad_vB[d] * inp.Parameters_OUT[d];
                GradUBGradPhi += Grad_uB[0, d] * inp.Parameters_OUT[d];
                nGradPhi += inp.Normale[d] * 0.5 * (inp.Parameters_IN[d] + inp.Parameters_OUT[d]);
            }
            Acc += 0.5 * (GradUAGradPhi + GradUBGradPhi) * (vA - vB) * nGradPhi;
            Acc += 0.5 * (GradVAGradPhi + GradVBGradPhi) * (uA[0] - uB[0]) * nGradPhi;


            //penaltyFactor in normal direction, see Ern,Stephansen,Zunino 2008 - eqns 2.12 - 2.14
            double deltaKA = 0;
            double deltaKB = 0;
            for (int d = 0; d < D; d++) {
                deltaKA = inp.Parameters_IN[d] * inp.Normale[d];
                deltaKB = inp.Parameters_OUT[d] * inp.Normale[d];
            }
            deltaKA = deltaKA * deltaKA;
            deltaKB = deltaKB * deltaKB;

            //harmonic average
            double Viscosity = (deltaKA == 0 && deltaKB == 0) ? 0 : deltaKA * deltaKB / (deltaKA + deltaKB);

            Acc -= (uA[0] - uB[0]) * (vA - vB) * pnlty * Viscosity; // penalty term


            return -Acc;
            //}
        }



    }


    /// <summary>
    /// Extension-Velocity Flux, the direction is based on the gradient of the level-set
    /// </summary>
    public class EllipticExtVelFormDirected : EllipticExtVelForm, IObserver<LevelSetTracker.LevelSetRegions>    {

        
        /// <summary>
        /// CTOR
        /// </summary>
        /// <param name="PenaltyBase">penalty multiplier based on polynomials</param>
        /// <param name="IsotropicViscosity">Value of the global isotropic viscosity</param>
        /// <param name="LSTrck"><see cref="LevelSetTracker"/></param>
        public EllipticExtVelFormDirected(double PenaltyBase, double IsotropicViscosity, Foundation.XDG.ILevelSetForm InterfaceFlux,  LevelSetTracker LSTrck) : base(PenaltyBase,IsotropicViscosity, LSTrck){

            _ParameterOrdering = ArrayTools.Cat(VariableNames.LevelSetGradient(D), VariableNames.LevelSet, VariableNames.MeanLevelSetGradient(D), InterfaceFlux.ParameterOrdering);
            this.InterfaceFlux = InterfaceFlux;

            // this is really ugly:
            // In the boundary-edge flux we evaluate the flux at the interface for inflow boundaries.
            // A clean solution would be to implement separate fluxes,
            //for each boundary condition or even an own IEquationcomponent just for the boundary.
            // However, this seems to work for now.
            if (InterfaceFlux.ParameterOrdering.First() == "InterfaceValue") {
                BoundaryFunc = inp => inp.Parameters_IN[2 * D + 1];
            }
            else if (InterfaceFlux.ParameterOrdering.First() == VariableNames.VelocityX) {
                // I am not sure, if this is correct...
                BoundaryFunc = ScalarVelocity;
            }
            else throw new NotImplementedException("Up to now: only SingleComponent and Scalar Velocity is supported");
        }

        double ScalarVelocity(CommonParamsBnd inp) {
            double S0 = 0; // velocity in normal direction, at the interface
            for (int d = 0; d < D; d++) {
                S0 += inp.Normale[d] * inp.Parameters_IN[d];
            }
            return S0;
        }

        Func<CommonParamsBnd, double> BoundaryFunc;

        Foundation.XDG.ILevelSetForm InterfaceFlux;

        IList<string> _ParameterOrdering;
        /// <summary>
        /// <see cref="Foundation.ILevelSetComponent.ParameterOrdering"/>
        /// </summary>
        public override IList<string> ParameterOrdering {
            get => _ParameterOrdering;
        }

        /// <summary>
        /// <see cref="swipViscosity_Term1"/>
        /// </summary>
        /// <param name="inp">
        /// </param>
        /// <param name="uIn">
        /// ExtensionVelocity
        /// </param>
        /// <param name="uOut"></param>
        /// <param name="Grad_uIn"></param>
        /// <param name="Grad_uOut"></param>
        /// <param name="vIn"></param>
        /// <param name="vOut"></param>
        /// <param name="Grad_vIn"></param>
        /// <param name="Grad_vOut"></param>
        /// <returns></returns>
        public override double InnerEdgeForm(ref CommonParams inp, double[] uIn, double[] uOut, double[,] Grad_uIn, double[,] Grad_uOut, double vIn, double vOut, double[] Grad_vIn, double[] Grad_vOut) {

            // Initialize all values, which are calculoated in the flux itself
            double Acc = 0.0; // => return
            double GradUInGradPhi = 0;
            double GradUOutGradPhi = 0;
            double GradVInGradPhi = 0;
            double GradVOutGradPhi = 0;
            double nGradPhiIn = 0;
            double nGradPhiOut = 0;
            double AbsGradPhiIn = 0;  // => divide by this to normalize GradPhi
            double AbsGradPhiOut = 0; // => divide by this to normalize GradPhi
            double pnlty = GetPenalty(inp.jCellIn, inp.jCellOut, this.cj);

            for (int d = 0; d < inp.D; d++) {
                GradUInGradPhi += Grad_uIn[0, d] * inp.Parameters_IN[d];
                GradVInGradPhi += Grad_vIn[d] * inp.Parameters_IN[d];
                GradVOutGradPhi += Grad_vOut[d] * inp.Parameters_OUT[d];
                GradUOutGradPhi += Grad_uOut[0, d] * inp.Parameters_OUT[d];
                nGradPhiIn += inp.Normale[d] * inp.Parameters_IN[d];
                nGradPhiOut += inp.Normale[d] * inp.Parameters_OUT[d];
                AbsGradPhiIn += inp.Parameters_IN[d] * inp.Parameters_IN[d];
                AbsGradPhiOut += inp.Parameters_OUT[d] * inp.Parameters_OUT[d];
            }
            AbsGradPhiIn = Math.Sqrt(AbsGradPhiIn) + myEps;
            AbsGradPhiOut = Math.Sqrt(AbsGradPhiOut) + myEps;
            GradUInGradPhi /= AbsGradPhiIn;
            GradUOutGradPhi /= AbsGradPhiOut;
            GradVInGradPhi /= AbsGradPhiIn;
            GradVOutGradPhi /= AbsGradPhiOut;
            nGradPhiIn /= AbsGradPhiIn;
            nGradPhiOut /= AbsGradPhiOut;

            // Calculate the Flow-Direction by the mean Value of GradPhi
            // DirectionSelector > 0 => In -> Out
            // DirectionSelector < 0 => Out -> In
            // ==========================================================
            double DirectionSelector = 0;
            double AbsGradPhiMeanIn = 0;
            double AbsGradPhiMeanOut = 0;
            double Direction_IN = 0;
            double Direction_OUT = 0;
            for (int d = 0; d < inp.D; d++){
                Direction_IN += inp.Parameters_IN[d + D + 1] * inp.Normale[d] * Math.Sign(inp.Parameters_IN[D]);
                Direction_OUT += inp.Parameters_OUT[d + D + 1] * inp.Normale[d] * Math.Sign(inp.Parameters_OUT[D]);
                AbsGradPhiMeanIn += inp.Parameters_IN[d + D + 1] * inp.Parameters_IN[d + D + 1];
                AbsGradPhiMeanOut += inp.Parameters_OUT[d + D + 1] * inp.Parameters_OUT[d + D + 1];
            }
            AbsGradPhiMeanIn = Math.Sqrt(AbsGradPhiMeanIn);
            AbsGradPhiMeanOut = Math.Sqrt(AbsGradPhiMeanOut);
            Direction_IN /= (AbsGradPhiMeanIn + myEps);
            Direction_OUT /= (AbsGradPhiMeanOut + myEps);



            // The information is always advancing from the cut-cell into the outer cells
            bool CellInCut = CutCellMask[inp.jCellIn];
            bool CellOutCut = CutCellMask[inp.jCellOut];

            // Bidirectional Cuppling in Cut-Cells => Central Differences Flux
            // ===============================================================
            if (CellInCut && CellOutCut)
            {
                Acc += 0.5 * (GradUInGradPhi + GradUOutGradPhi) * (vIn - vOut) * 0.5 * (nGradPhiIn + nGradPhiOut);
                Acc += 0.5 * (GradVInGradPhi + GradVOutGradPhi) * (uIn[0] - uOut[0]) * 0.5 * (nGradPhiIn + nGradPhiOut);
                Acc -= (uIn[0] - uOut[0]) * (vIn - vOut) * 1.5 * pnlty;

                for (int d = 0; d < inp.D; d++) {
                    Acc += IsotropicViscosity * 0.5 * (Grad_uIn[0, d] + Grad_uOut[0, d]) * (vIn - vOut) * inp.Normale[d];
                    Acc += IsotropicViscosity * 0.5 * (Grad_vIn[d] + Grad_vOut[d]) * (uIn[0] - uOut[0]) * inp.Normale[d];
                }
                Acc -= IsotropicViscosity * (uIn[0] - uOut[0]) * (vIn - vOut) * 1.5 * pnlty;

                return -Acc;
            }

            // No Backflow of information from Uncut to Cut Cells
            // =================================================
            if (CellInCut && !CellOutCut) {
                DirectionSelector = 1;
            }
            else if (!CellInCut && CellOutCut) {
                DirectionSelector = -1;
            }
            // Else: Choice of Flow-Direction based on Direction Selector
            // ==========================================================
            else {
                DirectionSelector = Direction_IN + Direction_OUT;
            }
            // Uwind-Coupling for Extension
            //=============================
            if (DirectionSelector > 0) {
                Acc += 0.5 * (GradUInGradPhi + GradUOutGradPhi) * (-vOut) * 0.5 * (nGradPhiIn + nGradPhiOut);
                Acc += (GradVOutGradPhi) * (uIn[0] - uOut[0]) * 0.5 * (nGradPhiIn + nGradPhiOut);
                Acc -= (uIn[0] - uOut[0]) * (-vOut) * pnlty;
            }
            else {
                Acc += 0.5 * (GradUInGradPhi + GradUOutGradPhi) * (vIn) * 0.5 * (nGradPhiIn + nGradPhiOut);
                Acc += (GradVInGradPhi) * (uIn[0] - uOut[0]) * 0.5 * (nGradPhiIn + nGradPhiOut);
                Acc -= (uIn[0] - uOut[0]) * (vIn) * pnlty ;
            }

            // Central Differences for Isotropic Viscosity
            //============================================
            for (int d = 0; d < inp.D; d++) {
                Acc += IsotropicViscosity * 0.5 * (Grad_uIn[0, d] + Grad_uOut[0, d]) * (vIn - vOut) * inp.Normale[d];
                Acc += IsotropicViscosity * 0.5 * (Grad_vIn[d] + Grad_vOut[d]) * (uIn[0] - uOut[0]) * inp.Normale[d];
            }
            Acc -= IsotropicViscosity * (uIn[0] - uOut[0]) * (vIn - vOut) * 1.5 * pnlty;

            return -Acc;
            
        }


        /// <summary>
        /// The Idea here is, that the artificial viscosity term needs to have a free boundary, i.e.
        /// the flux at the boundary edge must be controlled
        /// In Addition, we add an inflow boundary condition.
        /// </summary>
        /// <param name="inp"></param>
        /// <param name="uA"></param>
        /// <param name="Grad_uA"></param>
        /// <param name="vA"></param>
        /// <param name="Grad_vA"></param>
        /// <returns></returns>
        public override double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] uA, double[,] Grad_uA, double vA, double[] Grad_vA)
        {
            double Acc = 0; // => return
            // Initialize all values, which are calculoated in the flux itself
            double GradUGradPhi = 0;
            double GradVGradPhi = 0;
            double nGradPhi = 0;
            double AbsGradPhi = 0;  // => divide by this to normalize GradPhi

            for (int d = 0; d < inp.D; d++)
            {
                GradUGradPhi += Grad_uA[0, d] * inp.Parameters_IN[d];                    
                GradVGradPhi += Grad_vA[d] * inp.Parameters_IN[d];
                nGradPhi += inp.Normale[d] * inp.Parameters_IN[d];
                AbsGradPhi += inp.Parameters_IN[d] * inp.Parameters_IN[d];
            }
            AbsGradPhi = Math.Sqrt(AbsGradPhi) + myEps;
            GradUGradPhi /= AbsGradPhi;
            GradVGradPhi /= AbsGradPhi;
            nGradPhi /= AbsGradPhi;


            // Central Differences for Isotropic Viscosity
            //============================================
            for (int d = 0; d < inp.D; d++)
            {
                Acc += IsotropicViscosity * (Grad_uA[0, d]) * (vA) * inp.Normale[d];
                //   Acc += IsotropicViscosity *  (Grad_vA[d]   ) * (uA - uOut[0]) * inp.Normale[d];
            }
            Acc -= IsotropicViscosity * GradUGradPhi * vA * nGradPhi;



            // Calculate the Flow-Direction by the mean Value of GradPhi
            // DirectionSelector > 0 => In -> Out
            // DirectionSelector < 0 => Out -> In
            // ==========================================================

            double Direction = 0; 
            for (int d = 0; d < inp.D; d++){
                // No normalization required, since we are only using one component.
                //Thus the sign is indepoendent of the absolute value.
                Direction += inp.Parameters_IN[d + D + 1] * inp.Normale[d] * Math.Sign(inp.Parameters_IN[D]);
            }
            
            


            // Uwind-Coupling for Extension
            //=============================
            if (Direction < 0){
                double pnlty = 3 * GetPenalty(inp.jCellIn, -1, this.cj);
                Acc +=  GradUGradPhi * vA * nGradPhi;
                Acc += (GradVGradPhi) * (uA[0] - inp.Parameters_IN[2*D+1] ) * nGradPhi;
                Acc -= (uA[0] - BoundaryFunc(inp)) * (vA) * pnlty;
            }

            return -Acc;
        }
    }

    /// <summary>
    /// Extension-Velocity Flux, the direction is based on the cell-average of the level-set, which is equivalent to a finite-difference gradient
    /// </summary>
    public class EllipticExtVelFormLevelSetBased : EllipticExtVelForm {

        /// <summary>
        /// CTOR
        /// </summary>
        /// <param name="PenaltyBase">penalty multiplier based on polynomials</param>
        /// <param name="LSTrck"><see cref="LevelSetTracker"/></param>
        public EllipticExtVelFormLevelSetBased(double PenaltyBase, double IsotropicViscosity, LevelSetTracker LSTrck) : base(PenaltyBase, IsotropicViscosity, LSTrck) {
            CutCellBitMask = LSTrck.Regions.GetCutCellMask().GetBitMask();
        }

        BitArray CutCellBitMask;

        /// <summary>
        /// <see cref="IEquationComponent.ParameterOrdering"/>
        /// </summary>
        public override IList<string> ParameterOrdering
        {
            get
            {
                return ArrayTools.Cat(VariableNames.LevelSetGradient(D), VariableNames.LevelSet, "MeanLevelSet");
            }
        }

        /// <summary>
        /// <see cref="swipViscosity_Term1"/>
        ///Direction Selector:
        /// in Cut-Cells: Central Differences
        /// from Cut-Cells into Non-Cut-Cells: Upwinding away from Cut-cell
        /// Between Non-Cut-Cells: Upwinding according to Level-Set Direction
        /// </summary>
        /// <param name="inp">
        /// </param>
        /// <param name="uIn">
        /// ExtensionVelocity
        /// </param>
        /// <param name="uOut"></param>
        /// <param name="Grad_uIn"></param>
        /// <param name="Grad_uOut"></param>
        /// <param name="vIn"></param>
        /// <param name="vOut"></param>
        /// <param name="Grad_vIn"></param>
        /// <param name="Grad_vOut"></param>
        /// <returns></returns>
        public override double InnerEdgeForm(ref CommonParams inp, double[] uIn, double[] uOut, double[,] Grad_uIn, double[,] Grad_uOut, double vIn, double vOut, double[] Grad_vIn, double[] Grad_vOut) {
            double Acc = 0.0;
            double pnlty = 2* GetPenalty(inp.jCellIn, inp.jCellOut, this.cj);

            // Calculate Gradient \cdot XYZ
            double GradUInGradPhi = 0;
            double GradUOutGradPhi = 0;
            double GradVInGradPhi = 0;
            double GradVOutGradPhi = 0;
            double nGradPhiIn = 0;
            double nGradPhiOut = 0;
            double AbsGradPhiIn = 0;
            double AbsGradPhiOut = 0;
            for (int d = 0; d < inp.D; d++)
            {
                GradUInGradPhi += Grad_uIn[0, d] * inp.Parameters_IN[d];
                GradVInGradPhi += Grad_vIn[d] * inp.Parameters_IN[d];
                GradVOutGradPhi += Grad_vOut[d] * inp.Parameters_OUT[d];
                GradUOutGradPhi += Grad_uOut[0, d] * inp.Parameters_OUT[d];
                nGradPhiIn += inp.Normale[d] * inp.Parameters_IN[d];
                nGradPhiOut += inp.Normale[d] * inp.Parameters_OUT[d];
                AbsGradPhiIn += inp.Parameters_IN[d] * inp.Parameters_IN[d];
                AbsGradPhiOut += inp.Parameters_OUT[d] * inp.Parameters_OUT[d];
            }
            AbsGradPhiIn = Math.Sqrt(AbsGradPhiIn) + myEps;
            AbsGradPhiOut = Math.Sqrt(AbsGradPhiOut) + myEps;
            GradUInGradPhi /= AbsGradPhiIn;
            GradUOutGradPhi /= AbsGradPhiOut;
            GradVInGradPhi /= AbsGradPhiIn;
            GradVOutGradPhi /= AbsGradPhiOut;
            nGradPhiIn /= AbsGradPhiIn;
            nGradPhiOut /= AbsGradPhiOut;


            //penaltyFactor in normal direction, see Ern,Stephansen,Zunino 2008 - eqns 2.12 - 2.14
            double deltaKA = nGradPhiIn * nGradPhiIn;
            double deltaKB = nGradPhiOut * nGradPhiOut;

            double Viscosity = (deltaKA == 0 && deltaKB == 0) ? 0 : deltaKA * deltaKB / (deltaKA + deltaKB);


            // Direction Selector:
            // in Cut-Cells: Central Differences
            // from Cut-Cells into Non-Cut-Cells: Upwinding away from Cut-cell
            // Between Non-Cut-Cells: Upwinding according to Level-Set Direction


            // The information is always advancing from the cut-cell into the outer cells
            // Central Differences on Cut-Cells
            bool CellInCut = CutCellBitMask[inp.jCellIn];
            bool CellOutCut = CutCellBitMask[inp.jCellOut];

            if (CellInCut && CellOutCut) {
                Acc += 0.5 * (GradUInGradPhi + GradUOutGradPhi) * (vIn-vOut) * 0.5*(nGradPhiIn+nGradPhiOut);
                Acc += 0.5*(GradVInGradPhi+GradVOutGradPhi) * (vIn - vOut) * 0.5 * (nGradPhiIn + nGradPhiOut);
                Acc -= (uIn[0] - uOut[0]) * (vIn-vOut) * 1.5 * pnlty * Viscosity;
                return -Acc;
            }
            // levelSet*gradPhi*n
            double DirectionSelector = 0;
            //Upwinding

            if (CellInCut && !CellOutCut) {
                DirectionSelector = 1;
            }
            else if (!CellInCut && CellOutCut) {
                DirectionSelector = -1;
            }
            // Compute Flow-Direction for
            // Far-Cells and between Cut-Cells
            else {
                //for (int d = 0; d < inp.D; d++) {
                    // = (PhiMeanOut - PhiMeanIn)*(PhiMeanOut+PhiMeanIn) -> Idea: (FiniteDifferenceGradient)*(Mean(MeanLevelSetValue))
                    DirectionSelector += (inp.Parameters_OUT[D + 1] * inp.Parameters_OUT[D + 1] - inp.Parameters_IN[D + 1] * inp.Parameters_IN[D + 1]);
            }



            if (DirectionSelector > 0) {
                Acc += 0.5*(GradUInGradPhi+GradUOutGradPhi) * (-vOut) * nGradPhiIn;
                Acc +=      GradVOutGradPhi * (vIn-vOut) * nGradPhiIn;
                Acc -= (uIn[0] - uOut[0]) * (-vOut) * 1.5 * pnlty * Viscosity;
            }
            else {
                Acc += 0.5 * (GradUInGradPhi + GradUOutGradPhi) * (vIn) * nGradPhiIn;
                Acc += GradVInGradPhi * (vIn - vOut) * nGradPhiIn;
                Acc -= (uIn[0] - uOut[0]) * (vIn) *  pnlty * Viscosity;
            }
            return -Acc;

        }

    }

    /// <summary>
    /// <see cref="EllipticExtVelFormDirected"/>, but for a limited cellmask,
    /// all cells outside this mask are treated as boundary edges
    /// </summary>
    public class ExtVelForm_bulk : EllipticExtVelFormDirected {

        /// <summary>
        /// CTOR, <see cref="EllipticExtVelFormDirected"/>, but for a limited cellmask
        /// </summary>
        public ExtVelForm_bulk(double PenaltyBase, double IsotropicViscosity, Foundation.XDG.ILevelSetForm InterfaceFlux, LevelSetTracker LSTrck, BitArray mask) :base(PenaltyBase, IsotropicViscosity, InterfaceFlux, LSTrck) {
            this.mask = mask;
        }
        BitArray mask;

        /// <summary>
        /// <see cref="EllipticExtVelFormDirected.InnerEdgeForm(ref CommonParams, double[], double[], double[,], double[,], double, double, double[], double[])"/>
        /// </summary>
        public override double InnerEdgeForm(ref CommonParams inp, double[] uA, double[] uB, double[,] Grad_uA, double[,] Grad_uB, double vA, double vB, double[] Grad_vA, double[] Grad_vB) {
            if (mask[inp.jCellIn] && mask[inp.jCellOut]) {
                return base.InnerEdgeForm(ref inp, uA, uB, Grad_uA, Grad_uB, vA, vB, Grad_vA, Grad_vB);
            }
            else {
                // boundary of sub-domain
                return 0;
            }
        }
    }


}
