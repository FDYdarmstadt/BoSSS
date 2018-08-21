using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using System;
using System.Collections;
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
            IDictionary<SpeciesId, MultidimensionalArray> InterfaceLengthScales,
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
                InterfaceLengthScales,
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
            IDictionary<SpeciesId, MultidimensionalArray> InterfaceLengthScales,
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
                InterfaceLengthScales,
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
            Dictionary<SpeciesId, MultidimensionalArray> Idummy = dummy.XDGSpaceMetrics.CutCellMetrics.InterfaceArea;


            var bla = new Dictionary<SpeciesId, XSpatialOperator.QrSchemPair>();
            foreach(var sp in whichSpc)
                bla.Add(sp, new XSpatialOperator.QrSchemPair());

            ComputeMatrixEx<M, V>(
                xOp,
                lsTrk,
                DomainMap, Parameters, CodomainMap,
                Matrix, AffineOffset,
                OnlyAffine, time, MPIParameterExchange, bla,
                dummy.CellLengthScales, Idummy, subGrid);

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
            IDictionary<SpeciesId, MultidimensionalArray> InterfaceLengthScales,
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
                            CellLengthScales = CellLengthScales[s],
                            GrdDat = lsTrk.GridDat
                        });
                }
                ev.SpeciesOperatorCoefficients[s].UserDefinedValues["InterfaceLengths"] = InterfaceLengthScales[s];

                MultidimensionalArray slipLengths = lsTrk.GridDat.Cells.h_min.CloneAs();
                slipLengths.Clear();
                double hmin = lsTrk.GridDat.Cells.h_min.To1DArray().Min();
                int degU = DomainMap.BasisS.ToArray()[0].Degree;
                CellMask ContactArea = lsTrk.Regions.GetNearFieldMask(1).Intersect(lsTrk.GridDat.BoundaryCells.VolumeMask);
                foreach(Chunk cnk in ContactArea) {
                    for(int i = cnk.i0; i < cnk.JE; i++) {
                        slipLengths[i] = hmin / (degU + 1);
                    }
                }
                ev.SpeciesOperatorCoefficients[s].UserDefinedValues["SlipLengths"] = slipLengths;
            }


            if(OnlyAffine)
                ev.ComputeAffine(AffineOffset);
            else
                ev.ComputeMatrix(Matrix, AffineOffset);

        }
    }
}
