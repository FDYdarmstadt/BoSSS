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
using static BoSSS.Foundation.DifferentialOperator;
using ilPSP.Tracing;

namespace BoSSS.Foundation {
    
    /// <summary>
    /// extension methods for the <see cref="IEquationComponent"/>-interface
    /// </summary>
    public static class IEquationComponentExtension {

        /// <summary>
        /// creates the spatial operator that consists only of component <paramref name="c"/>
        /// </summary>
        public static DifferentialOperator Operator(this IEquationComponent c, int DegreeOfNonlinearity = 1, Func<IGridData,EdgeQuadratureScheme> edgeSchemeProvider = null,  Func<IGridData,CellQuadratureScheme> cellSchemeProvider = null, bool doCommit = true)
        {
            return Operator(c, QuadOrderFunc.NonLinear(DegreeOfNonlinearity), edgeSchemeProvider, cellSchemeProvider, doCommit : doCommit);
        }


        /// <summary>
        /// creates the spatial operator that consists only of component <paramref name="c"/>
        /// </summary>
        public static DifferentialOperator Operator(this IEquationComponent c, Func<int[], int[], int[], int> quadOrderFunc, Func<IGridData,EdgeQuadratureScheme> edgeSchemeProvider = null,  Func<IGridData,CellQuadratureScheme> cellSchemeProvider = null, bool doCommit = true) {

            string[] Codomain = new string[] { "v1" };
            string[] Domain = c.ArgumentOrdering.ToArray();
            string[] Param = (c.ParameterOrdering != null) ? c.ParameterOrdering.ToArray() : new string[0];

            DifferentialOperator ret = new DifferentialOperator(Domain, Param, Codomain, quadOrderFunc);
            if(edgeSchemeProvider != null)
                ret.EdgeQuadraturSchemeProvider = edgeSchemeProvider;
            if(cellSchemeProvider != null)
                ret.VolumeQuadraturSchemeProvider = cellSchemeProvider;

            ret.EquationComponents[Codomain[0]].Add(c);
            if(doCommit)
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
        /// <see cref="DifferentialOperator.ParameterVar"/>);  It is allowed to set
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
        /// <see cref="IDifferentialOperator.Commit"/> before this method can be called.
        /// </remarks>
        /// <param name="edgeQuadScheme">
        /// Quadrature scheme for edge integration; 
        /// </param>
        /// <param name="volQuadScheme">
        /// Quadrature scheme for volume integration;
        /// </param>
        /// <param name="SubGridBoundaryMask">
        /// </param>
        /// <param name="ParameterMPIExchange">
        /// Determines whether parameter fields have to exchange ghost cell
        /// data before the assembly of the operator.
        /// </param>
        /// <param name="time"></param>
        /// <param name="op"></param>
        static public void ComputeMatrixEx<M, V>(this DifferentialOperator op,
            UnsetteledCoordinateMapping DomainMap, IList<DGField> Parameters, UnsetteledCoordinateMapping CodomainMap,
            M Matrix, V AffineOffset, bool OnlyAffine = false,
            double time = 0.0, 
            EdgeQuadratureScheme edgeQuadScheme = null, CellQuadratureScheme volQuadScheme = null,
            //ICompositeQuadRule<QuadRule> edgeRule = null, ICompositeQuadRule<QuadRule> volRule = null,
            BitArray SubGridBoundaryMask = null,
            bool ParameterMPIExchange = true)
            where M : IMutableMatrixEx
            where V : IList<double> //
        {

            var bkup = op.LegacySupport_ModifyQuadSchemProvider(edgeQuadScheme, volQuadScheme);

            var ev = op.GetMatrixBuilder(DomainMap, Parameters, CodomainMap);
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

            op.LegacySupport_RestoreQuadSchemeProvider(bkup);
        }

        /// <summary>
        /// another legacy interface
        /// </summary>
        static public void ComputeAffine<V>(
            this DifferentialOperator op,
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


        /// <summary>
        /// Evaluates this operator for given DG fields;
        /// </summary>
        /// <param name="op"></param>
        /// <param name="time"></param>
        /// <param name="DomainMapping">
        /// the domain variables, or "input data" for the operator; the number
        /// of elements must be equal to the number of elements in the <see cref="IDifferentialOperator.DomainVar"/>-list;
        /// </param>
        /// <param name="Params">
        /// List of parameter fields; May be null
        /// </param>
        /// <param name="CodomainMapping">
        /// the co-domain variables, or "output" for the evaluation of the
        /// operator; the number of elements must be equal to the number of elements in the <see cref="IDifferentialOperator.CodomainVar"/>-list;
        /// </param>
        /// <param name="alpha">
        /// scaling of the operator 
        /// </param>
        /// <param name="beta">
        /// scaling of the accumulator (<paramref name="CodomainMapping"/>);
        /// </param>
        /// <param name="sgrd">
        /// subgrid, for restricted evaluation; null indicates evaluation on
        /// the full grid.
        /// </param>
        /// <param name="bndMode">
        /// Treatment of subgrid boundaries, if <paramref name="sgrd"/> is not null. See <see cref="IEvaluator.SubGridBoundaryTreatment"/>
        /// </param>
        /// <param name="qInsEdge">
        /// Optional definition of the edge quadrature scheme. Since this
        /// already implies a domain of integration, must be null if
        /// <paramref name="sgrd"/> is not null.
        /// </param>
        /// <param name="qInsVol">
        /// Optional definition of the volume quadrature scheme. Since this
        /// already implies a domain of integration, must be null if
        /// <paramref name="sgrd"/> is not null.
        /// </param>
        /// <remarks>
        /// If some of the input data, <paramref name="DomainMapping"/>, is
        /// contained in the  output data, <paramref name="CodomainMapping"/>,
        /// these DG fields will be cloned to ensure correct operation of the
        /// operator evaluation.
        /// </remarks>
        static public void Evaluate(
                             this DifferentialOperator op,
                             double alpha, double beta,
                             CoordinateMapping DomainMapping, IList<DGField> Params, CoordinateMapping CodomainMapping,
                             SubGrid sgrd = null,
                             EdgeQuadratureScheme qInsEdge = null,
                             CellQuadratureScheme qInsVol = null,
                             SubGridBoundaryModes bndMode = SubGridBoundaryModes.OpenBoundary,
                             double time = double.NaN) //
        {
            using(new FuncTrace()) {
                if(sgrd != null && (qInsEdge != null || qInsVol != null))
                    throw new ArgumentException("Specification of Subgrid and quadrature schemes is exclusive: not allowed to specify both at the same time.", "sgrd");
#if DEBUG
                op.Verify(false); 
#endif
                

                IList<DGField> _DomainFields = DomainMapping.Fields;
                IList<DGField> _CodomainFields = CodomainMapping.Fields;
                DGField[] _DomainFieldsRevisited = new DGField[_DomainFields.Count];

                bool a = false;
                for(int i = 0; i < _DomainFields.Count; i++) {
                    DGField f = _DomainFields[i];

                    if(_CodomainFields.Contains(f)) {
                        // some of the domain variables (input data)
                        // is also a member of the codomain variables (output);
                        // the data need to be cloned to provide correct results
                        a = true;
                        _DomainFieldsRevisited[i] = (DGField)f.Clone();
                    } else
                        _DomainFieldsRevisited[i] = f;
                }

                CoordinateMapping domainMappingRevisited;
                if(a)
                    domainMappingRevisited = new CoordinateMapping(_DomainFieldsRevisited);
                else
                    domainMappingRevisited = DomainMapping;

                if(sgrd != null) {
                    CellMask cm = (sgrd == null) ? null : sgrd.VolumeMask;
                    EdgeMask em = (sgrd == null) ? null : sgrd.AllEdgesMask;

                    qInsEdge = new EdgeQuadratureScheme(true, em);
                    qInsVol = new CellQuadratureScheme(true, cm);
                }

                var bkup = op.LegacySupport_ModifyQuadSchemProvider(qInsEdge, qInsVol);


                var ev = op.GetEvaluatorEx(
                    domainMappingRevisited, Params, CodomainMapping);

                if(sgrd != null)
                    ev.ActivateSubgridBoundary(sgrd.VolumeMask, bndMode);
                CoordinateVector outp = new CoordinateVector(CodomainMapping);
                ev.time = time;
                ev.Evaluate<CoordinateVector>(alpha, beta, outp);

                op.LegacySupport_RestoreQuadSchemeProvider(bkup);
            }
        }

        /// <summary>
        /// Another wrapper for operator evaluation
        /// </summary>
        /// <param name="DomainFields">
        /// the domain variables; the number of elements
        /// must be equal to the number of elements in the <see cref="IDifferentialOperator.DomainVar"/>-list;
        /// </param>
        /// <param name="CodomainFields">
        /// the codomain variables; the number of elements
        /// must be equal to the number of elements in the <see cref="IDifferentialOperator.CodomainVar"/>-list;
        /// </param>
        /// <remarks>
        /// If some of the input field, i.e. some element of <paramref name="DomainFields"/>, is contained in the 
        /// output fields, i.e. in the list <paramref name="CodomainFields"/>, these DG fields will be cloned to ensure 
        /// correct operation of the operator evaluation.<br/>
        /// It is not a good choice to use this 
        /// function if this operator should be evaluated multiple times
        /// and 
        /// contains linear components (i.e. <see cref="IDifferentialOperator.IsLinear"/> returns true);
        /// If the later one is the case, the matrix which represents
        /// the linear components of the operator must be computed first, which is computational- and
        /// memory - intensive;
        /// After execution of this method, the matrix will be lost;
        /// If multiple evaluation is desired, the <see cref="IDifferentialOperator.GetEvaluatorEx"/>-method should be used.
        /// </remarks>
        /// <param name="op"></param>
        /// <param name="time"></param>
        static public void Evaluate(this DifferentialOperator op, IList<DGField> DomainFields, IList<DGField> CodomainFields, double time = double.NaN) {
            if(DomainFields.Count != op.DomainVar.Count)
                throw new ArgumentException("wrong number of domain fields", "DomainFields");
            if(CodomainFields.Count != op.CodomainVar.Count)
                throw new ArgumentException("wrong number of domain fields", "CodomainFields");

            CoordinateMapping inp = new CoordinateMapping(DomainFields);
            CoordinateMapping outp = new CoordinateMapping(CodomainFields);

            Evaluate(op, 1.0, 0.0, inp, null, outp, null, null, null, SubGridBoundaryModes.OpenBoundary, time);
        }

        /// <summary>
        /// And another wrapper.
        /// </summary>
        static public void Evaluate(this DifferentialOperator op, params DGField[] f) {
            Evaluate(op, double.NaN, f);
        }


        /// <summary>
        /// And another wrapper.
        /// </summary>
        static public void Evaluate(this DifferentialOperator op, double time, params DGField[] f) {
            Evaluate(op, time, null, SubGridBoundaryModes.OpenBoundary, f);
        }

        /// <summary>
        /// And another wrapper.
        /// </summary>
        static public void Evaluate(this DifferentialOperator op, double time, SubGrid subGrid = null, SubGridBoundaryModes subGridBoundaryMode = SubGridBoundaryModes.OpenBoundary, params DGField[] f) {
            if(op.DomainVar.Count + op.ParameterVar.Count + op.CodomainVar.Count != f.Length)
                throw new ArgumentException("wrong number of domain/parameter/codomain fields", "f");

            CoordinateMapping inp = new CoordinateMapping(f.GetSubVector(0, op.DomainVar.Count));
            DGField[] Parameters = f.GetSubVector(op.DomainVar.Count, op.ParameterVar.Count);
            CoordinateMapping outp = new CoordinateMapping(f.GetSubVector(op.DomainVar.Count + op.ParameterVar.Count, op.CodomainVar.Count));

            Evaluate(op, 1.0, 0.0, inp, Parameters, outp, subGrid, null, null, subGridBoundaryMode, time);
        }


        /// <summary>
        /// And another wrapper.
        /// </summary>
        static public void Evaluate(this DifferentialOperator op, uint NoOfDomainVar, uint NoOfCodomainVar, double time, params DGField[] fields) {
            if(fields.Length != (NoOfDomainVar + NoOfCodomainVar))
                throw new ArgumentException("wrong number of fields.", "fields");

            DGField[] DomF = new DGField[NoOfDomainVar];
            Array.Copy(fields, 0, DomF, 0, NoOfDomainVar);
            DGField[] CoDomF = new DGField[NoOfCodomainVar];
            Array.Copy(fields, NoOfDomainVar, CoDomF, 0, NoOfCodomainVar);

            Evaluate(op, DomF, CoDomF, time);
        }


        /// <summary>
        /// simplified version of 
        /// <see cref="ComputeMatrixEx{M,V}"/>;
        /// </summary>
        static public void ComputeMatrix<M, V>(this DifferentialOperator op,
            UnsetteledCoordinateMapping DomainMap, IList<DGField> Parameters, UnsetteledCoordinateMapping CodomainMap,
            M Matrix, V AffineOffset, bool OnlyAffine, double time = 0.0,
            SubGrid sgrd = null)
            where M : IMutableMatrixEx
            where V : IList<double> {

            var bkup = op.LegacySupport_ModifyQuadSchemProvider(
                new EdgeQuadratureScheme(true, sgrd == null ? null : sgrd.AllEdgesMask),
                new CellQuadratureScheme(true, sgrd == null ? null : sgrd.VolumeMask));



            var ev = op.GetMatrixBuilder(DomainMap, Parameters, CodomainMap);

            ev.time = time;

            if(OnlyAffine) {
                ev.ComputeAffine(AffineOffset);
            } else {
                ev.ComputeMatrix(Matrix, AffineOffset);
            }

            op.LegacySupport_RestoreQuadSchemeProvider(bkup);
        }

        private static IGridData CheckArguments(UnsetteledCoordinateMapping DomainMap, IList<DGField> Parameters, UnsetteledCoordinateMapping CodomainMap) {
            var GridDat = DomainMap.GridDat;
            if(!object.ReferenceEquals(GridDat, CodomainMap.GridDat))
                throw new ArgumentException("Domain and codomain map must be assigend to the same grid.");
            if(Parameters != null)
                foreach(var prm in Parameters)
                    if(prm != null && (!object.ReferenceEquals(prm.GridDat, GridDat)))
                        throw new ArgumentException(string.Format("parameter field {0} is assigned to a different grid.", prm.Identification));
            return GridDat;
        }


        /// <summary>
        /// simplified version of 
        /// <see cref="ComputeMatrixEx{M,V}"/>;
        /// </summary>
        static public BlockMsrMatrix ComputeMatrix(this DifferentialOperator op,
            UnsetteledCoordinateMapping DomainMap, IList<DGField> Parameters, UnsetteledCoordinateMapping CodomainMap, double time = 0.0, SubGrid SubGrid = null) {

            //int RowBlkSize = (CodomainMap.MaxTotalNoOfCoordinatesPerCell == CodomainMap.MinTotalNoOfCoordinatesPerCell) ? CodomainMap.MaxTotalNoOfCoordinatesPerCell : 1;
            //int ColBlkSize = (DomainMap.MaxTotalNoOfCoordinatesPerCell == DomainMap.MinTotalNoOfCoordinatesPerCell) ? DomainMap.MaxTotalNoOfCoordinatesPerCell : 1;

            BlockMsrMatrix Matrix = new BlockMsrMatrix(CodomainMap, DomainMap);

            double[] dummyRHS = new double[Matrix.RowPartitioning.LocalLength];
            ComputeMatrix(op, 
                DomainMap, Parameters, CodomainMap,
                Matrix, dummyRHS,
                false, time, SubGrid);

            return Matrix;
        }

        /// <summary>
        /// legacy interface
        /// </summary>
        static public double[] ComputeAffine(this DifferentialOperator op,
            UnsetteledCoordinateMapping DomainMap, IList<DGField> Parameters, UnsetteledCoordinateMapping CodomainMap, double time = 0.0) {
            
            int RowBlkSize = (CodomainMap.MaxTotalNoOfCoordinatesPerCell == CodomainMap.MinTotalNoOfCoordinatesPerCell) ? CodomainMap.MaxTotalNoOfCoordinatesPerCell : 1;
            int ColBlkSize = (DomainMap.MaxTotalNoOfCoordinatesPerCell == DomainMap.MinTotalNoOfCoordinatesPerCell) ? DomainMap.MaxTotalNoOfCoordinatesPerCell : 1;

            MsrMatrix Matrix = new MsrMatrix(CodomainMap.LocalLength, (int)DomainMap.GlobalCount, RowBlkSize, ColBlkSize);

            double[] RHS = new double[Matrix.RowPartitioning.LocalLength];
            //ComputeMatrix(op, 
            //    DomainMap, Parameters, CodomainMap,
            //    Matrix, dummyRHS,
            //    false, time, SubGrid);
            op.ComputeAffine(DomainMap, Parameters, CodomainMap, RHS, false, time);

            return RHS;
        }

       
    }
}
