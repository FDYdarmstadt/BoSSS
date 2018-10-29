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
using BoSSS.Platform.LinAlg;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using BoSSS.Solution.Utils;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using CNS.Convection;
using CNS.Diffusion;
using CNS.MaterialProperty;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;

namespace CNS.IBM {

    /// <summary>
    /// Abstract base class for source terms that represent boundary conditions
    /// at immersed interfaces. 
    /// </summary>
    public abstract class BoundaryConditionSource : IEquationComponent {

        /// <summary>
        /// Configuration options
        /// </summary>
        protected CNSControl config;

        /// <summary>
        /// Information about the location of the immersed boundary
        /// </summary>
        protected ISpeciesMap speciesMap;

        /// <summary>
        /// The boundary condition to be enforced
        /// </summary>
        protected BoundaryCondition boundaryCondition;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="config"></param>
        /// <param name="speciesMap"></param>
        /// <param name="boundaryCondition"></param>
        public BoundaryConditionSource(CNSControl config, ISpeciesMap speciesMap, BoundaryCondition boundaryCondition) {
            this.config = config;
            this.speciesMap = speciesMap;
            this.boundaryCondition = boundaryCondition;
        }

        #region IEquationComponent Members

        /// <summary>
        /// <see cref="CNSEnvironment.PrimalArgumentOrdering"/>
        /// </summary>
        public IList<string> ArgumentOrdering {
            get {
                return CNSEnvironment.PrimalArgumentOrdering;
            }
        }

        /// <summary>
        /// The level set gradient
        /// </summary>
        public IList<string> ParameterOrdering {
            get {
                return Enumerable.
                    Range(0, CNSEnvironment.NumberOfDimensions).
                    Select(i => "levelSetGradient" + i).
                    ToList();
            }
        }

        #endregion
    }

    /// <summary>
    /// Helper methods for <see cref="BoundaryConditionSource"/>
    /// </summary>
    public static class IEquationComponentExtensions {

        /// <summary>
        /// Dispatcher to determine the correct sub-class of
        /// <see cref="BoundaryConditionSource"/> when creating a source term
        /// for a specific equation component.
        /// </summary>
        public static BoundaryConditionSource CreateBoundaryConditionSource(
            this IEquationComponent equationComponent, CNSControl control, ISpeciesMap speciesMap, BoundaryCondition boundaryCondition) {
            BoundaryConditionSource result;

            if (equationComponent is INonlinearFlux) {
                result = new BoundaryConditionSourceFromINonlinearFlux(
                    control, speciesMap, boundaryCondition, (INonlinearFlux)equationComponent);
            } else if (equationComponent is SIPGFlux) {
                result = new BoundaryConditionSourceFromSIPGFlux(
                    control, speciesMap, boundaryCondition, (SIPGFlux)equationComponent);
            } else if (equationComponent is INonlinear2ndOrderForm) {
                result = new BoundaryConditionSourceFromINonlinear2ndOrderForm(
                    control, speciesMap, boundaryCondition, (INonlinear2ndOrderForm)equationComponent);
            } else {
                throw new NotImplementedException("To do");
            }

            return result;
        }
    }

    /// <summary>
    /// Implements a source term that enforces a boundary condition at the zero
    /// iso-contour of a level set. Note that this source term does not include
    /// the multiplication with the delta distribution (i.e., this has to be
    /// accounted for during quadrature)
    /// </summary>
    public class BoundaryConditionSourceFromINonlinearFlux : BoundaryConditionSource, INonlinearSource {

        /// <summary>
        /// The flux function to be used to evaluate the flux across the zero
        /// level set (by making use of
        /// <see cref="EulerFlux.InnerEdgeFlux(double[], double, StateVector, StateVector, ref Vector, int)"/>)
        /// </summary>
        private INonlinearFlux fluxFunction;

        /// <summary>
        /// Constructs a new source term for a given
        /// <paramref name="boundaryCondition"/>
        /// </summary>
        /// <param name="config">
        /// Configuration options
        /// </param>
        /// <param name="speciesMap">
        /// Information about the location of the immersed boundary
        /// </param>
        /// <param name="boundaryCondition">
        /// The boundary condition to be enforced
        /// </param>
        /// <param name="fluxFunction">
        /// The flux function to be used to evaluate the flux across the zero
        /// level set (by making use of
        /// <see cref="NonlinearFlux.InnerEdgeFlux(double, int, MultidimensionalArray, MultidimensionalArray, MultidimensionalArray[], MultidimensionalArray[], int, int, MultidimensionalArray)"/>)
        /// </param>
        public BoundaryConditionSourceFromINonlinearFlux(
            CNSControl config, ISpeciesMap speciesMap, BoundaryCondition boundaryCondition, INonlinearFlux fluxFunction)
            : base(config, speciesMap, boundaryCondition) {
            this.fluxFunction = fluxFunction;
        }

        #region INonlinearSource Members

        /// <summary>
        /// Evaluates the configured flux function (see
        /// <see cref="BoundaryConditionSourceFromINonlinearFlux.BoundaryConditionSourceFromINonlinearFlux(CNSControl, ISpeciesMap, BoundaryCondition, INonlinearFlux)"/>
        /// using the present flow state <paramref name="U"/> and the boundary
        /// value provided by the configured boundary condition.
        /// </summary>
        /// <param name="time"></param>
        /// <param name="x"></param>
        /// <param name="U"></param>
        /// <param name="IndexOffset"></param>
        /// <param name="FirstCellInd"></param>
        /// <param name="Lenght"></param>
        /// <param name="Output"></param>
        public void Source(double time, MultidimensionalArray x, MultidimensionalArray[] U, int IndexOffset, int FirstCellInd, int Lenght, MultidimensionalArray Output) {
            int D = CNSEnvironment.NumberOfDimensions;
            int noOfNodes = x.GetLength(1);
            int noOfVariables = D + 2;

            MultidimensionalArray normal = MultidimensionalArray.Create(Lenght, noOfNodes, D);
            MultidimensionalArray[] Uout = new MultidimensionalArray[noOfVariables];
            for (int i = 0; i < noOfVariables; i++) {
                Uout[i] = MultidimensionalArray.Create(Lenght, noOfNodes);
            }

            double[] xLocal = new double[D];
            Material material = speciesMap.GetMaterial(double.NaN);
            for (int i = 0; i < Lenght; i++) {
                for (int j = 0; j < noOfNodes; j++) {
                    StateVector stateIn = new StateVector(material, U, i, j, CNSEnvironment.NumberOfDimensions);

                    Vector levelSetNormal = new Vector(CNSEnvironment.NumberOfDimensions);
                    int offset = CNSEnvironment.NumberOfDimensions + 2;
                    for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                        levelSetNormal[d] = U[offset + d][i + IndexOffset, j];
                    }
                    levelSetNormal.Normalize();
                    Debug.Assert(Math.Abs(levelSetNormal.Abs() - 1.0) < 1e-13, "Abnormal normal vector");

                    for (int d = 0; d < D; d++) {
                        xLocal[d] = x[i + IndexOffset, j, d];
                        normal[i, j, d] = levelSetNormal[d];
                    }

                    StateVector stateOut = boundaryCondition.GetBoundaryState(time, xLocal, levelSetNormal, stateIn);
                    Debug.Assert(stateOut.IsValid, "Invalid boundary state");

                    Uout[0][i, j] = stateOut.Density;
                    for (int d = 0; d < D; d++) {
                        Uout[d + 1][i, j] = stateOut.Momentum[d];
                    }
                    Uout[D + 1][i, j] = stateOut.Energy;
                }
            }

            fluxFunction.InnerEdgeFlux(time, -1, x, normal, U, Uout, IndexOffset, Lenght, Output);
        }

        #endregion
    }

    /// <summary>
    /// Variant of <see cref="BoundaryConditionSourceFromINonlinearFlux"/> for
    /// sub-classes of <see cref="SIPGFlux"/>
    /// </summary>
    public class BoundaryConditionSourceFromSIPGFlux : BoundaryConditionSource, IVolumeForm {

        private SIPGFlux fluxFunction;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="config"></param>
        /// <param name="speciesMap"></param>
        /// <param name="boundaryCondition"></param>
        /// <param name="fluxFunction"></param>
        public BoundaryConditionSourceFromSIPGFlux(
            CNSControl config, ISpeciesMap speciesMap, BoundaryCondition boundaryCondition, SIPGFlux fluxFunction)
            : base(config, speciesMap, boundaryCondition) {
            this.fluxFunction = fluxFunction;
        }

        /// <summary>
        /// See <see cref="SIPGFlux.BoundaryEdgeTerms"/>
        /// </summary>
        public TermActivationFlags VolTerms {
            get {
                return fluxFunction.BoundaryEdgeTerms;
            }
        }

        /// <summary>
        /// Evaluates the configured flux function (see
        /// <see cref="BoundaryConditionSourceFromSIPGFlux.BoundaryConditionSourceFromSIPGFlux(CNSControl, ISpeciesMap, BoundaryCondition, SIPGFlux)"/>
        /// using the present flow state <paramref name="U"/> and the boundary
        /// value provided by the configured boundary condition.
        /// </summary>
        /// <param name="cpv"></param>
        /// <param name="U"></param>
        /// <param name="GradU"></param>
        /// <param name="V"></param>
        /// <param name="GradV"></param>
        /// <returns></returns>
        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            int D = CNSEnvironment.NumberOfDimensions;

            double[] normal = new double[D];
            double abs = 0.0;
            for (int d = 0; d < D; d++) {
                normal[d] = cpv.Parameters[d];
                abs += normal[d] * normal[d];
            }
            abs = Math.Sqrt(abs);

            Debug.Assert(abs > 1e-10, "Extremely flat level set gradient");

            for (int d = 0; d < D; d++) {
                normal[d] /= abs;
            }

            Material material = speciesMap.GetMaterial(double.NaN);
            StateVector stateIn = new StateVector(U, material);
            StateVector stateOut = boundaryCondition.GetBoundaryState(cpv.time, cpv.Xglobal, normal, stateIn);
            Debug.Assert(stateOut.IsValid, "Invalid boundary state");

            CommonParams InParams = new CommonParams() {
                GridDat = cpv.GridDat,
                iEdge = Math.Abs(cpv.GridDat.iLogicalCells.Cells2Edges[cpv.jCell][0]) - 1, // TO BE CHANGED
                Normale = normal,
                Parameters_IN = cpv.Parameters,
                Parameters_OUT = cpv.Parameters,
                time = cpv.time,
                X = cpv.Xglobal
            };

            // cf. SIPGFlux, line 206
            double[,] GradUBoundary = GradU;

            // cf. SIPGFlux, line 207 & 210
            double VBoundary = 0.0;

            // cf. SIPGFlux, line 203
            // BEWARE: _Not_ zero since we want {gradV} = gradVIn
            double[] GradVBoundary = GradV;



            //return fluxFunction.InnerEdgeForm(
            //    ref InParams, U, stateOut.ToArray(), GradU, GradUBoundary, V, VBoundary, GradV, GradVBoundary);


            SIPGFlux.EVIL_HACK_CELL_INDEX = cpv.jCell;
            double flux = fluxFunction.InnerEdgeForm(
                ref InParams, U, stateOut.ToArray(), GradU, GradUBoundary, V, VBoundary, GradV, GradVBoundary);
            SIPGFlux.EVIL_HACK_CELL_INDEX = -1;
            return flux;
        }
    }

    /// <summary>
    /// Variant of <see cref="BoundaryConditionSourceFromINonlinearFlux"/> for
    /// sub-classes of <see cref="INonlinear2ndOrderForm"/>
    /// </summary>
    public class BoundaryConditionSourceFromINonlinear2ndOrderForm : BoundaryConditionSource, INonlinVolumeForm_GradV, INonlinVolumeForm_V {

        private INonlinear2ndOrderForm fluxFunction;

        private bool adiaWall;

        /// <summary>
        /// Constructor
        /// </summary>
        /// <param name="config"></param>
        /// <param name="speciesMap"></param>
        /// <param name="boundaryCondition"></param>
        /// <param name="fluxFunction"></param>
        public BoundaryConditionSourceFromINonlinear2ndOrderForm(
            CNSControl config, ISpeciesMap speciesMap, BoundaryCondition boundaryCondition, INonlinear2ndOrderForm fluxFunction)
            : base(config, speciesMap, boundaryCondition) {
            this.fluxFunction = fluxFunction;

            if (boundaryCondition is AdiabaticWall) {
                adiaWall = true;
            }
        }

        /// <summary>
        /// See <see cref="IVolumeForm.VolTerms"/>
        /// </summary>
        TermActivationFlags IVolumeForm.VolTerms {
            get {
                return fluxFunction.InnerEdgeTerms;
            }
        }

        /// <summary>
        /// Not used
        /// </summary>
        /// <param name="cpv"></param>
        /// <param name="U"></param>
        /// <param name="GradU"></param>
        /// <param name="V"></param>
        /// <param name="GradV"></param>
        /// <returns></returns>
        double IVolumeForm.VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            // Not used
            throw new NotImplementedException();
        }

        /// <summary>
        /// Passes the given parameters to <see cref="INonlinEdgeForm_V.InternalEdge"/>
        /// </summary>
        /// <param name="prm"></param>
        /// <param name="U"></param>
        /// <param name="GradU"></param>
        /// <param name="f"></param>
        void INonlinVolumeForm_V.Form(ref VolumFormParams prm, MultidimensionalArray[] U, MultidimensionalArray[] GradU, MultidimensionalArray f) {
            INonlinEdgeForm_V flux = fluxFunction;

            MultidimensionalArray[] UBoundary;
            MultidimensionalArray normals;
            EdgeFormParams efp;
            AdaptParameters(ref prm, U, GradU, out efp, out UBoundary, out normals);

            MultidimensionalArray[] GradUBoundary = GradU; // cf. SIPGFlux, line 206

            // Set fBoundary to zero
            MultidimensionalArray fBoundary = MultidimensionalArray.Create(
                U[0].GetLength(0), prm.Xglobal.GetLength(1));


            OptimizedSIPGEnergyFlux.EVIL_HACK_CELL_INDEX = prm.j0;
            OptimizedSIPGMomentumFlux.EVIL_HACK_CELL_INDEX = prm.j0;
            fluxFunction.AdiabaticWall = this.adiaWall;
            flux.InternalEdge(ref efp, U, UBoundary, GradU, GradUBoundary, f, fBoundary);
            OptimizedSIPGEnergyFlux.EVIL_HACK_CELL_INDEX = -1;
            OptimizedSIPGMomentumFlux.EVIL_HACK_CELL_INDEX = -1;
        }

        /// <summary>
        /// Passes the given parameters to <see cref="INonlinEdgeForm_GradV.InternalEdge"/>
        /// </summary>
        /// <param name="prm"></param>
        /// <param name="U"></param>
        /// <param name="GradU"></param>
        /// <param name="f"></param>
        void INonlinVolumeForm_GradV.Form(ref VolumFormParams prm, MultidimensionalArray[] U, MultidimensionalArray[] GradU, MultidimensionalArray f) {
            INonlinEdgeForm_GradV flux = fluxFunction;

            MultidimensionalArray[] UBoundary;
            MultidimensionalArray normals;
            EdgeFormParams efp;
            AdaptParameters(ref prm, U, GradU, out efp, out UBoundary, out normals);

            MultidimensionalArray[] GradUBoundary = GradU; // cf. SIPGFlux, line 206

            // Set fBoundary to zero
            MultidimensionalArray fBoundary = MultidimensionalArray.Create(
                U[0].GetLength(0), prm.Xglobal.GetLength(1), CNSEnvironment.NumberOfDimensions);
            fluxFunction.AdiabaticWall = this.adiaWall;
            flux.InternalEdge(ref efp, U, UBoundary, GradU, GradUBoundary, f, fBoundary);
        }

        /// <summary>
        /// Reformulates the given parameters into <paramref name="efp"/>,
        /// <paramref name="UBoundary"/> and <paramref name="normals"/>, which
        /// are in the form required by
        /// <see cref="INonlinEdgeForm_GradV.InternalEdge"/>
        /// and <see cref="INonlinEdgeForm_V.InternalEdge"/>
        /// </summary>
        /// <param name="prm"></param>
        /// <param name="U"></param>
        /// <param name="GradU"></param>
        /// <param name="efp"></param>
        /// <param name="UBoundary"></param>
        /// <param name="normals"></param>
        private void AdaptParameters(ref VolumFormParams prm, MultidimensionalArray[] U, MultidimensionalArray[] GradU, out EdgeFormParams efp, out MultidimensionalArray[] UBoundary, out MultidimensionalArray normals) {

            Debug.Assert(U[0].GetLength(0) == 1, "Number of cells must be 1");
            Debug.Assert(prm.Len == 1, "Number of cells must be 1");

            INonlinEdgeForm_GradV flux = fluxFunction;
            int noOfCells = 1;
            int noOfNodesPerCell = prm.Xglobal.GetLength(1);

            UBoundary = new MultidimensionalArray[U.Length];
            for (int k = 0; k < U.Length; k++) {
                UBoundary[k] = MultidimensionalArray.Create(noOfCells, noOfNodesPerCell);
            }

            normals = MultidimensionalArray.Create(
                noOfCells, noOfNodesPerCell, CNSEnvironment.NumberOfDimensions);
            Material material = speciesMap.GetMaterial(double.NaN);
            for (int j = 0; j < noOfNodesPerCell; j++) {
                double[] x = new double[CNSEnvironment.NumberOfDimensions];
                double[] normal = new double[CNSEnvironment.NumberOfDimensions];

                double abs = 0.0;
                for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                    x[d] = prm.Xglobal[0, j, d];
                    normal[d] = prm.ParameterVars[d][0, j];
                    abs += normal[d] * normal[d];
                }
                abs = Math.Sqrt(abs);

                Debug.Assert(abs > 1e-10, "Extremely flat level set gradient");

                for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                    normal[d] /= abs;
                }

                StateVector stateIn = new StateVector(material, U, 0, j, CNSEnvironment.NumberOfDimensions);
                StateVector stateBoundary = boundaryCondition.GetBoundaryState(
                    prm.time, x, normal, stateIn);
                Debug.Assert(stateBoundary.IsValid, "Invalid boundary state");

                double[] UBoundaryLocal = stateBoundary.ToArray();
                for (int k = 0; k < U.Length; k++) {
                    UBoundary[k][0, j] = UBoundaryLocal[k];
                }

                for (int d = 0; d < CNSEnvironment.NumberOfDimensions; d++) {
                    normals[0, j, d] = normal[d];
                }
            }

            efp = new EdgeFormParams() {
                e0 = Math.Abs(prm.GridDat.iLogicalCells.Cells2Edges[prm.j0][0]) - 1, // THIS IS AN EVIL HACK; NEEDS TO BE CHANGED
                GridDat = prm.GridDat,
                Len = prm.Len,
                NodesGlobal = prm.Xglobal,
                Normals = normals,
                ParameterVars_IN = prm.ParameterVars,
                ParameterVars_OUT = prm.ParameterVars,
                time = prm.time
            };
        }
    }
}
