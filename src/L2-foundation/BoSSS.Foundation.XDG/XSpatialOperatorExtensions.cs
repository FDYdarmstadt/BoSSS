using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// Extensions/legacy interfaced for the spatial operator
    /// </summary>
    public static class XSpatialOperatorExtensions {

        /// <summary>
        /// Legacy interface, preserved as static function.
        /// </summary>
        public static void ComputeMatrixEx<M, V>(
            this XSpatialOperator xOp,
            LevelSetTracker lsTrk,
            UnsetteledCoordinateMapping DomainMap, IList<DGField> Parameters, UnsetteledCoordinateMapping CodomainMap,
            M Matrix, V AffineOffset, bool OnlyAffine, double time, bool ParameterMPIExchange,
            IDictionary<SpeciesId, MultidimensionalArray> CellLengthScales,
            SubGrid SubGrid, params SpeciesId[] whichSpc)
            where M : IMutableMatrixEx
            where V : IList<double> //
        {
            var SpeciesDictionary = new Dictionary<SpeciesId, XSpatialOperator.QrSchemPair>();
            foreach(var sp in whichSpc)
                SpeciesDictionary.Add(sp, new XSpatialOperator.QrSchemPair());

            ComputeMatrixEx<M, V>(
                xOp,
                lsTrk,
                DomainMap, Parameters, CodomainMap,
                Matrix, AffineOffset, OnlyAffine, time,
                ParameterMPIExchange, SpeciesDictionary, CellLengthScales,
                //agg, out mass,
                SubGrid);

        }

        /// <summary>
        /// Legacy interface, preserved as static function.
        /// </summary>
        static public void ComputeMatrixEx<M, V>(
            this XSpatialOperator xOp,
            LevelSetTracker lsTrk,
            UnsetteledCoordinateMapping DomainMap, IList<DGField> Parameters, UnsetteledCoordinateMapping CodomainMap,
            M Matrix, V AffineOffset, bool OnlyAffine,
            double time, bool ParameterMPIExchange,
            IDictionary<SpeciesId, MultidimensionalArray> CellLengthScales,
            params SpeciesId[] whichSpc)
            where M : IMutableMatrixEx
            where V : IList<double> //
        {
            ComputeMatrixEx<M, V>(
                xOp,
                lsTrk,
                DomainMap, Parameters, CodomainMap,
                Matrix, AffineOffset, OnlyAffine, time,
                ParameterMPIExchange,
                CellLengthScales,
                null, whichSpc);
        }

        /// <summary>
        /// Legacy interface, preserved as static function.
        /// </summary>
        static public void ComputeMatrixEx<M, V>(
            this XSpatialOperator xOp,
            LevelSetTracker lsTrk,
            UnsetteledCoordinateMapping DomainMap, IList<DGField> Parameters, UnsetteledCoordinateMapping CodomainMap,
            M Matrix, V AffineOffset, bool OnlyAffine, double time, bool MPIParameterExchange,
            params SpeciesId[] whichSpc)
            where M : IMutableMatrixEx
            where V : IList<double> //
        {
            ComputeMatrixEx<M, V>(
               xOp,
               lsTrk,
               DomainMap, Parameters, CodomainMap,
               Matrix, AffineOffset,
               OnlyAffine, time, MPIParameterExchange, subGrid: null, whichSpc: whichSpc
               );
        }


        /// <summary>
        /// Legacy interface, preserved as static function.
        /// </summary>
        static public void ComputeMatrixEx<M, V>(
            this XSpatialOperator xOp,
            LevelSetTracker lsTrk,
            UnsetteledCoordinateMapping DomainMap, IList<DGField> Parameters, UnsetteledCoordinateMapping CodomainMap,
            M Matrix, V AffineOffset, bool OnlyAffine,
            double time,
            bool MPIParameterExchange,
            SubGrid subGrid, params SpeciesId[] whichSpc)
            where M : IMutableMatrixEx
            where V : IList<double>  //
        {
            Console.WriteLine("Warning: Legacy interface!");
            Console.WriteLine("         Usage of XDG without an agglomerated penalty!");


            int order = xOp.GetOrderFromQuadOrderFunction(DomainMap, Parameters, CodomainMap);
            MultiphaseCellAgglomerator dummy = lsTrk.GetAgglomerator(lsTrk.SpeciesIdS.ToArray(), order, 0.0);

            var bla = new Dictionary<SpeciesId, XSpatialOperator.QrSchemPair>();
            foreach(var sp in whichSpc)
                bla.Add(sp, new XSpatialOperator.QrSchemPair());

            ComputeMatrixEx<M, V>(
                xOp,
                lsTrk,
                DomainMap, Parameters, CodomainMap,
                Matrix, AffineOffset,
                OnlyAffine, time, MPIParameterExchange, bla,
                dummy.CellLengthScales, subGrid);

            Debug.Assert(dummy.TotalNumberOfAgglomerations <= 0, "internal error");

        }


        /// <summary>
        /// Legacy interface, preserved as static function.
        /// </summary>
        static public void ComputeMatrixEx<M, V>(
            this XSpatialOperator xOp,
            LevelSetTracker lsTrk,
            UnsetteledCoordinateMapping DomainMap, IList<DGField> Parameters, UnsetteledCoordinateMapping CodomainMap,
            M Matrix, V AffineOffset, bool OnlyAffine, double time, bool MPIParameterExchange,
            IDictionary<SpeciesId, XSpatialOperator.QrSchemPair> SpeciesSchemes, IDictionary<SpeciesId, MultidimensionalArray> CellLengthScales,
            SubGrid SubGrid = null)
            where M : IMutableMatrixEx
            where V : IList<double> //
        {
            CellMask SubGridCellMask = null;
            EdgeMask SubGridEdgeMask = null;
            if(SubGrid != null) {
                SubGridCellMask = SubGrid.VolumeMask;
                /// I don't know why, but this seems to work:
                SubGridEdgeMask = SubGrid.AllEdgesMask;
                /// And this does not:
                //SubGridEdgeMask = SubGrid.InnerEdgesMask;
            }

            SpeciesId[] ReqSpecies = SpeciesSchemes.Keys.ToArray();

            //int order = xOp.GetOrderFromQuadOrderFunction(DomainMap, Parameters, CodomainMap);
            //var SchemeHelper = lsTrk.GetXDGSpaceMetrics(ReqSpecies, order, 1).XQuadSchemeHelper;

            //var SpeciesSchemes_out = new Dictionary<SpeciesId, XSpatialOperator.QrSchemPair>();

            //foreach(var SpeciesId in SpeciesSchemes.Keys) {
            //    EdgeQuadratureScheme edgeScheme;
            //    CellQuadratureScheme cellScheme;

            //    var qrSchemes = SpeciesSchemes[SpeciesId];

            //    bool AssembleOnFullGrid = (SubGrid == null);
            //    if(qrSchemes.EdgeScheme == null) {
            //        edgeScheme = SchemeHelper.GetEdgeQuadScheme(SpeciesId, AssembleOnFullGrid, SubGridEdgeMask);
            //    } else {
            //        //edgeScheme = qrSchemes.EdgeScheme;
            //        //throw new ArgumentException();
            //    }

            //    if(qrSchemes.CellScheme == null) {
            //        cellScheme = SchemeHelper.GetVolumeQuadScheme(SpeciesId, AssembleOnFullGrid, SubGridCellMask);
            //    } else {
            //        //cellScheme = qrSchemes.CellScheme;
            //        throw new ArgumentException();
            //    }

            //    SpeciesSchemes_out.Add(SpeciesId,
            //        new XSpatialOperator.QrSchemPair() {
            //            CellScheme = cellScheme,
            //            EdgeScheme = edgeScheme
            //        });
      
            //}


            var ev = xOp.GetMatrixBuilder(lsTrk, DomainMap, Parameters, CodomainMap, SpeciesSchemes);


            ev.time = time;
            ev.MPITtransceive = MPIParameterExchange;

            foreach(var s in CellLengthScales.Keys) {
                if (ev.SpeciesOperatorCoefficients.ContainsKey(s)) {
                    ev.SpeciesOperatorCoefficients[s].CellLengthScales = CellLengthScales[s];
                } else {
                    ev.SpeciesOperatorCoefficients.Add(s,
                        new CoefficientSet() {
                            CellLengthScales = CellLengthScales[s]
                        });
                }
            }

            if(OnlyAffine)
                ev.ComputeAffine(AffineOffset);
            else
                ev.ComputeMatrix(Matrix, AffineOffset);

        }
    }
}
