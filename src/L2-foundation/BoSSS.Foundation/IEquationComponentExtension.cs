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
using ilPSP.Utils;
using BoSSS.Foundation.Quadrature;
using System.Collections;
using ilPSP.LinSolvers;
using BoSSS.Foundation.Grid;

namespace BoSSS.Foundation {
    
    /// <summary>
    /// extension methods for the <see cref="IEquationComponent"/>-interface
    /// </summary>
    public static class IEquationComponentExtension {

        /// <summary>
        /// creates the spatial operator that consists only of component <paramref name="c"/>
        /// </summary>
        public static SpatialOperator Operator(this IEquationComponent c, int DegreeOfNonlinearity = 1)
        {
            return Operator(c, QuadOrderFunc.NonLinear(DegreeOfNonlinearity));
        }


        /// <summary>
        /// creates the spatial operator that consists only of component <paramref name="c"/>
        /// </summary>
        public static SpatialOperator Operator(this IEquationComponent c, Func<int[], int[], int[], int> quadOrderFunc) {

            string[] Codomain = new string[] { "v1" };
            string[] Domain = c.ArgumentOrdering.ToArray();
            string[] Param = (c.ParameterOrdering != null) ? c.ParameterOrdering.ToArray() : new string[0];

            SpatialOperator ret = new SpatialOperator(Domain, Param, Codomain, quadOrderFunc);
            ret.EquationComponents[Codomain[0]].Add(c);
            ret.Commit();

            return ret;
        }

       

        /// <summary>
        /// computes the affine offset and/or matrix of the operator, expert
        /// version;
        /// </summary>
        /// <param name="DomainMap">
        /// the mapping which is used to compute column indices into
        /// <paramref name="Matrix"/>;
        /// </param>
        /// <param name="Parameters">
        /// The parameter variables (of this differential operator);
        /// The number of elements in the list must match the parameter count
        /// of the differential operator (see
        /// <see cref="SpatialOperator.ParameterVar"/>);  It is allowed to set
        /// an entry to 'null', in this case the values of the parameter field
        /// are assumed to be 0.0; If the differential operator contains no
        /// parameters, this argument can be null;
        /// </param>
        /// <param name="CodomainMap">
        /// the mapping which is used to compute row indices into
        /// <paramref name="Matrix"/> and <paramref name="AffineOffset"/>.
        /// </param>
        /// <param name="Matrix">
        /// Acc output: the matrix which represents the linear part of this
        /// operator, according to the mapping given by
        /// <paramref name="DomainMap"/> and <paramref name="CodomainMap"/>,
        /// is <b>ACCUMULATED</b> here; <br/>
        /// Setting all matrix entries to 0.0 is left to the user;
        /// </param>
        /// <param name="AffineOffset">
        /// Acc output: the vector which represents the affine part of this
        /// operator, according to the mapping given by
        /// <paramref name="DomainMap"/> and <paramref name="CodomainMap"/>,
        /// is <b>ACCUMULATED</b> here; <br/>
        /// Setting all vector entries to 0.0 is left to the user;
        /// </param>
        /// <param name="OnlyAffine">
        /// If true, only the <paramref name="AffineOffset"/> is computed, and
        /// the <paramref name="Matrix"/> is not touched (can be null);
        /// </param>
        /// <remarks>
        /// The operator assembly must be finalized before by calling
        /// <see cref="Commit"/> before this method can be called.
        /// </remarks>
        /// <param name="edgeRule">
        /// Quadrature rule and domain for edge integration; specifying this is exclusive with <paramref name="edgeQuadScheme"/>, i.e. both cannot be unequal null at the same time.
        /// </param>
        /// <param name="edgeQuadScheme">
        /// Quadrature scheme for edge integration; specifying this is exclusive with <paramref name="edgeRule"/>, i.e. both cannot be unequal null at the same time.
        /// </param>
        /// <param name="volRule">
        /// Quadrature rule and domain for volume integration; specifying this is exclusive with <paramref name="volQuadScheme"/>, i.e. both cannot be unequal null at the same time.
        /// </param>
        /// <param name="volQuadScheme">
        /// Quadrature scheme for volume integration; specifying this is exclusive with <paramref name="volRule"/>, i.e. both cannot be unequal null at the same time.
        /// </param>
        /// <param name="SubGridBoundaryMask">
        /// </param>
        /// <param name="ParameterMPIExchange">
        /// Determines whether parameter fields have to exchange ghost cell
        /// data before the assembly of the operator.
        /// </param>
        /// <param name="time"></param>
        /// <param name="op"></param>
        static public void ComputeMatrixEx<M, V>(this SpatialOperator op,
            UnsetteledCoordinateMapping DomainMap, IList<DGField> Parameters, UnsetteledCoordinateMapping CodomainMap,
            M Matrix, V AffineOffset,bool OnlyAffine = false,
            double time = 0.0, 
            EdgeQuadratureScheme edgeQuadScheme = null, CellQuadratureScheme volQuadScheme = null,
            //ICompositeQuadRule<QuadRule> edgeRule = null, ICompositeQuadRule<QuadRule> volRule = null,
            BitArray SubGridBoundaryMask = null,
            bool ParameterMPIExchange = true)
            where M : IMutableMatrixEx
            where V : IList<double> //
        {
           
            var ev = op.GetMatrixBuilder(DomainMap, Parameters, CodomainMap, edgeQuadScheme, volQuadScheme);
            ev.time = time;
            ev.MPITtransceive = ParameterMPIExchange;
            if(SubGridBoundaryMask != null) {
                throw new NotSupportedException();
                //ev.ActivateSubgridBoundary(new Grid.CellMask(ev.GridData, SubGridBoundaryMask));
            }

            if(OnlyAffine)
                ev.ComputeAffine(AffineOffset);
            else
                ev.ComputeMatrix(Matrix, AffineOffset);
        }

        /// <summary>
        /// another legacy interface
        /// </summary>
        static public void ComputeAffine<V>(
            this SpatialOperator op,
            UnsetteledCoordinateMapping DomainMap,
            IList<DGField> Parameters,
            UnsetteledCoordinateMapping CodomainMap,
            V AffineOffset,
            bool OnlyBoundaryEdges = true, double time = 0.0,
            EdgeQuadratureScheme edgeQr = null, CellQuadratureScheme volQr = null)
            where V : IList<double> {

            var GridDat = CodomainMap.GridDat;
            if (Parameters != null)
                foreach (var prm in Parameters)
                    if (!object.ReferenceEquals(prm.GridDat, GridDat))
                        throw new ArgumentException(string.Format("parameter field {0} is assigned to a different grid.", prm.Identification));

            //Using order zero for DomainMap will lead to inconsistent (and possibly insufficient) quadrature order!!! 
            //UnsetteledCoordinateMapping DomainMap;
            //Basis b = new Basis(GridDat, 0);
            //Basis[] B = new Basis[this.DomainVar.Count];
            //B.SetAll(b);
            //DomainMap = new UnsetteledCoordinateMapping(B);


            if (OnlyBoundaryEdges) {
                if (edgeQr != null)
                    throw new ArgumentException("If 'OnlyBoundaryEdges == true', 'edgeQr' must be null!", "edgeQr");
                if (volQr != null)
                    throw new ArgumentException("If 'OnlyBoundaryEdges == true', 'volQr' must be null!", "volQr");

                volQr = new CellQuadratureScheme(true, CellMask.GetEmptyMask(GridDat));
                edgeQr = new EdgeQuadratureScheme(true, GridDat.GetBoundaryEdgeMask());
            }

            op.ComputeMatrixEx(
                DomainMap, Parameters, CodomainMap,
                default(MsrMatrix), AffineOffset,
                OnlyAffine: true, time:time,
                volQuadScheme: volQr, edgeQuadScheme: edgeQr);
        }


        ///// <summary>
        ///// Legacy interface
        ///// </summary>
        //static public void ComputeMatrixEx<M, V>(
        //    this SpatialOperator op,
        //    UnsetteledCoordinateMapping DomainMap, IList<DGField> Parameters, UnsetteledCoordinateMapping CodomainMap,
        //    M Matrix, V AffineOffset, bool OnlyAffine = false,
        //    double time = 0.0, 
        //    EdgeQuadratureScheme edgeQuadScheme = null, CellQuadratureScheme volQuadScheme = null,
        //    ICompositeQuadRule<QuadRule> edgeRule = null, ICompositeQuadRule<QuadRule> volRule = null,
        //    BitArray SubGridBoundaryMask = null,
        //    bool ParameterMPIExchange = true)
        //    where M : IMutableMatrix
        //    where V : IList<double> //
        //{
        //    //

        //    if(edgeQuadScheme != null && edgeRule != null)
        //        throw new ArgumentException();
        //    if(volQuadScheme != null && volRule != null)
        //        throw new ArgumentException();

        //    //if (edgeQrCtx == null)
        //    //    edgeQrCtx = new EdgeQuadratureScheme(true);

        //    var GridDat = CheckArguments(DomainMap, Parameters, CodomainMap);

        //    int order = 0;
        //    if(edgeRule == null || volRule == null) 
        //        order = this.GetOrderFromQuadOrderFunction(DomainMap, Parameters, CodomainMap);
            
        //    if(edgeRule == null)
        //        edgeRule = edgeQuadScheme.SaveCompile(GridDat, order);
        //    if(volRule == null)
        //        volRule = volQuadScheme.SaveCompile(GridDat, order);

        //    Internal_ComputeMatrixEx(GridDat, DomainMap, Parameters, CodomainMap, Matrix, AffineOffset, OnlyAffine, time,
        //        edgeRule,
        //        volRule,
        //        SubGridBoundaryMask, ParameterMPIExchange);
        //}
    }
}
