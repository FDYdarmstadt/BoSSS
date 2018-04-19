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

using ilPSP;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using BoSSS.Solution.Utils;
using ilPSP.Utils;

namespace BoSSS.Solution.LevelSetTools.Reinit.FastMarch {
    
    /// <summary>
    /// Computation of level-set gradients
    /// </summary>
    public class GradientModule {

        /// <summary>
        /// Derivative computation in a single cell, see <see cref="GradientUpdate(int, BitArray, SinglePhaseField, VectorField{SinglePhaseField})"/>,
        ///  downwind/upwind based on 'accepted'-bitmask.
        /// </summary>
        class Gradient : BoSSS.Foundation.IEdgeForm, BoSSS.Foundation.IVolumeForm {
            public Gradient(int Direction, int jCell, BitArray AcceptedMask) {
                m_Direction = Direction;
                m_jCell = jCell;
                m_AcceptedMask = AcceptedMask;
                Debug.Assert(AcceptedMask[jCell] == false);
            }

            int m_Direction;
            int m_jCell;
            BitArray m_AcceptedMask;


            public TermActivationFlags BoundaryEdgeTerms {
                get {
                    return TermActivationFlags.UxV | TermActivationFlags.V;
                }
            }

            public TermActivationFlags InnerEdgeTerms {
                get {
                    return TermActivationFlags.UxV | TermActivationFlags.V;
                }
            }

            public double InnerEdgeForm(ref CommonParams inp, double[] PhiIn, double[] PhiOt, double[,] GradPhiIn, double[,] GradPhiOt, double vIn, double vOt, double[] Grad_vIn, double[] Grad_vOt) {

                if(this.m_AcceptedMask[inp.jCellIn] == true && this.m_AcceptedMask[inp.jCellOut] == false) {
                    // IN-cell is accepted // OUT-cell should be computed
                    // => the boundary value is given as IN-parameter
                    // => flux penalizes the OUT-Cell

                    Debug.Assert(inp.jCellOut == m_jCell);

                    return PhiIn[0] * (-inp.Normale[this.m_Direction]) * vOt;

                } else if(this.m_AcceptedMask[inp.jCellIn] == false && this.m_AcceptedMask[inp.jCellOut] == true) {
                    // ... vice-versa

                    Debug.Assert(inp.jCellIn == m_jCell);

                    return PhiOt[0] * (inp.Normale[this.m_Direction]) * vIn;

                } else {
                    
                    if(this.m_jCell == inp.jCellOut) {

                        return PhiOt[0] * (-inp.Normale[this.m_Direction]) * vOt;

                    } else if(this.m_jCell == inp.jCellIn) {

                        return PhiIn[0] * (+inp.Normale[this.m_Direction]) * vIn;

                    } else {
                        throw new ApplicationException();
                    }
                }
            }

            public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] PhiIn, double[,] GradPhiIn, double vIn, double[] Grad_vIn) {
                return PhiIn[0] * (+inp.Normale[this.m_Direction]) * vIn;
            }

            public IList<string> ArgumentOrdering {
                get {
                    return new string[] { "Phi" };
                }
            }

            public IList<string> ParameterOrdering {
                get {
                    return null;
                }
            }

            public TermActivationFlags VolTerms {
                get {
                    return TermActivationFlags.UxGradV | TermActivationFlags.GradV;
                }
            }

            public double VolumeForm(ref CommonParamsVol cpv, double[] Phi, double[,] GradPhi, double V, double[] GradV) {
                return (-Phi[0] * GradV[this.m_Direction]);
            }
        }


        /// <summary>
        /// Derivative computation, see <see cref="GradientUpdate(SubGrid, double[], SinglePhaseField, VectorField{SinglePhaseField})"/>,
        /// downwind/upwind based on mean value of level-set-field.
        /// </summary>
        class Gradient2 : BoSSS.Foundation.IEdgeForm, BoSSS.Foundation.IVolumeForm {
            public Gradient2(int Direction, double[] PhiMean) {
                m_PhiMean = PhiMean;
                m_Direction = Direction;
            }

            double[] m_PhiMean;
            int m_Direction;


            public TermActivationFlags BoundaryEdgeTerms {
                get {
                    return TermActivationFlags.UxV | TermActivationFlags.V;
                }
            }

            public TermActivationFlags InnerEdgeTerms {
                get {
                    return TermActivationFlags.UxV | TermActivationFlags.V;
                }
            }

            public double InnerEdgeForm(ref CommonParams inp, double[] PhiIn, double[] PhiOt, double[,] GradPhiIn, double[,] GradPhiOt, double vIn, double vOt, double[] Grad_vIn, double[] Grad_vOt) {
                
                double signIn = Math.Sign(m_PhiMean[inp.jCellIn]);
                double signOt = Math.Sign(m_PhiMean[inp.jCellOut]);
                double Flux;

                if(signIn != signOt || (signIn * signOt == 0)) {
                    // central difference

                    Flux = 0.5 * (PhiIn[0] + PhiOt[0]) * (inp.Normale[this.m_Direction]);
                } else {
                    Debug.Assert(Math.Abs(signIn) == 1);
                    Debug.Assert(Math.Abs(signOt) == 1);
                    Debug.Assert(signIn == signOt);
                    double sign = signIn;

                    double delta = (m_PhiMean[inp.jCellIn] - m_PhiMean[inp.jCellOut]) * sign;


                    if(delta < 0)
                        Flux = PhiIn[0] * (inp.Normale[this.m_Direction]);
                    else
                        Flux = PhiOt[0] * (inp.Normale[this.m_Direction]);
                }
                
                return Flux * (vIn - vOt);

            }

            public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] PhiIn, double[,] GradPhiIn, double vIn, double[] Grad_vIn) {
                return PhiIn[0] * (+inp.Normale[this.m_Direction]) * vIn;
            }

            public IList<string> ArgumentOrdering {
                get {
                    return new string[] { "Phi" };
                }
            }

            public IList<string> ParameterOrdering {
                get {
                    return null;
                }
            }

            public TermActivationFlags VolTerms {
                get {
                    return TermActivationFlags.UxGradV | TermActivationFlags.GradV;
                }
            }

            public double VolumeForm(ref CommonParamsVol cpv, double[] Phi, double[,] GradPhi, double V, double[] GradV) {
                return (-Phi[0] * GradV[this.m_Direction]);
            }
        }

        /// <summary>
        /// Update of level-set gradient in 'upwind'-direction, i.e. at boundaries towards accepted cells,
        /// the outer value is taken.
        /// </summary>
        /// <param name="jCell">Cell index to update.</param>
        /// <param name="AcceptedMask"></param>
        /// <param name="Phi">Input: the level-set</param>
        /// <param name="gradPhi">Output: gradient of <paramref name="Phi"/> in cell <paramref name="jCell"/>.</param>
        public void GradientUpdate(int jCell, BitArray AcceptedMask, SinglePhaseField Phi, VectorField<SinglePhaseField> gradPhi) {
            var GridDat = Phi.GridDat;
            
            if(m_gradEvo == null || jCell != m_gradEvo_jCell) {
                var Sgrd = new SubGrid(new CellMask(GridDat, Chunk.GetSingleElementChunk(jCell)));

                

                SpatialOperator op = new SpatialOperator(1, 2, QuadOrderFunc.Linear(), "Phi", "g0", "g1");
                op.EquationComponents["g0"].Add(new Gradient(0, jCell, AcceptedMask));
                op.EquationComponents["g1"].Add(new Gradient(1, jCell, AcceptedMask));
                op.Commit();

                m_gradEvo = op.GetEvaluatorEx(
                    Phi.Mapping, null, gradPhi.Mapping,
                    edgeQrCtx: (new EdgeQuadratureScheme(domain: Sgrd.AllEdgesMask)),
                    volQrCtx: (new CellQuadratureScheme(domain: Sgrd.VolumeMask)),
                    subGridBoundaryTreatment: SpatialOperator.SubGridBoundaryModes.BoundaryEdge, sgrd:Sgrd);

                m_gradEvo_jCell = jCell;
            }

            foreach(var f in gradPhi) {
                f.Coordinates.ClearRow(jCell);
            }
            m_gradEvo.Evaluate(1.0, 1.0, gradPhi.CoordinateVector, 0.0, MPIexchange: false);
        }


        /// <summary>
        /// 
        /// </summary>
        public void GradientUpdate(SubGrid Sgrd, double[] PhiMean, SinglePhaseField Phi, VectorField<SinglePhaseField> gradPhi) {
            var GridDat = Phi.GridDat;
            
            
            gradPhi.Clear(Sgrd.VolumeMask);

            SpatialOperator op = new SpatialOperator(1, 2, QuadOrderFunc.Linear(), "Phi", "g0", "g1");
            op.EquationComponents["g0"].Add(new Gradient2(0, PhiMean));
            op.EquationComponents["g1"].Add(new Gradient2(1, PhiMean));
            op.Commit();

            var gradEvo = op.GetEvaluatorEx(
                Phi.Mapping, null, gradPhi.Mapping,
                (new EdgeQuadratureScheme(domain: Sgrd.AllEdgesMask)),
                (new CellQuadratureScheme(domain: Sgrd.VolumeMask)));

            gradEvo.ActivateSubgridBoundary(Sgrd.VolumeMask, SpatialOperator.SubGridBoundaryModes.BoundaryEdge);

            //Sgrd.VolumeMask.ToTxtFile("nar.csv", false);

            
            gradPhi.Clear(Sgrd.VolumeMask);
            gradEvo.time = 0.0;
            gradEvo.MPITtransceive = false;
            gradEvo.Evaluate(1.0, 0.0, gradPhi.CoordinateVector);
            //gradPhi.GradientByFlux(1.0, Phi, optionalSubGrid:Sgrd , bndMode: SpatialOperator.SubGridBoundaryModes.BoundaryEdge);

            
        }






        IEvaluatorNonLin m_gradEvo;
        int m_gradEvo_jCell;
    }
}
