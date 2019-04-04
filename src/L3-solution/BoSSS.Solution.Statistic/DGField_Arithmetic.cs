using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Statistic;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.Statistic {
    /// <summary>
    /// Arithmetic operation for DG fields on different meshes
    /// </summary>
    public static class DGField_Arithmetic {


        /// <summary>
        /// Addition/Subtraction of DG fields on different grids
        /// </summary>
        /// <param name="A"></param>
        /// <param name="scaleA"></param>
        /// <param name="B"></param>
        /// <param name="scaleB"></param>
        /// <returns>
        /// <paramref name="scaleA"/>*<paramref name="A"/> + <paramref name="scaleB"/>*<paramref name="B"/>
        /// </returns>
        static public DGField ScaledSummation(this DGField A, double scaleA, DGField B, double scaleB) {
            if(object.ReferenceEquals(A.GridDat, B.GridDat)) {
                // ++++++++++++++++++++++++++++
                // both fields on the same grid
                // ++++++++++++++++++++++++++++

                DGField sum;

                if(A.Basis.IsSubBasis(B.Basis)) {
                    sum = (DGField)B.Clone();
                    sum.Scale(scaleB);
                    sum.AccLaidBack(1.0, A);
                } else if(B.Basis.IsSubBasis(A.Basis)) {
                    sum = (DGField)A.Clone();
                    sum.Scale(scaleA);
                    sum.AccLaidBack(1.0, B);
                } else {
                    throw new ApplicationException("can't add the two fields, because their basis are incompatible");
                }


                sum.Identification = "(" + scaleA + "*" + A.Identification + "+" + scaleB + "*" + B.Identification + ")";

                return sum;
            } else {
                // ++++++++++++++++++++++++++
                // fields on different grids
                // ++++++++++++++++++++++++++

                DGField fine, coarse;
                double aF, aC;
                if(A.GridDat.CellPartitioning.TotalLength > B.GridDat.CellPartitioning.TotalLength) {
                    fine = A;
                    coarse = B;
                    aF = scaleA;
                    aC = scaleB;
                } else {
                    coarse = A;
                    fine = B;
                    aC = scaleA;
                    aF = scaleB;
                }

                Foundation.Grid.Classic.GridData fineGridData = ExtractGridData(fine.GridDat);
                Foundation.Grid.Classic.GridData coarseGridData = ExtractGridData(coarse.GridDat);

                DGFieldComparison.ComputeFine2CoarseMap(
                        fineGridData,
                        coarseGridData,
                        out var Fine2CoarseMapS);

                DGField injected;
                if(coarse is ConventionalDGField) {
                    ConventionalDGField _injected = new SinglePhaseField(
                        new Basis(fine.GridDat, Math.Max(coarse.Basis.Degree, fine.Basis.Degree)), 
                        coarse.Identification);

                    DGFieldComparison.InjectDGField(Fine2CoarseMapS, _injected, coarse as ConventionalDGField);
                    injected = _injected;
                } else if(coarse is XDGField) {
                    XDGField _injected = new XDGField(
                        new XDGBasis((fine as XDGField).Basis.Tracker, Math.Max(coarse.Basis.Degree, fine.Basis.Degree)), 
                        coarse.Identification);

                    DGFieldComparison.InjectXDGField(Fine2CoarseMapS, _injected, coarse as XDGField);
                    injected = _injected;
                } else {
                    throw new NotSupportedException();
                }

                return ScaledSummation(injected, aC, fine, aF);
            }
        }

        static Foundation.Grid.Classic.GridData ExtractGridData(Foundation.Grid.IGridData gridData)
        {
            if(gridData is Foundation.Grid.Classic.GridData classicGridData)
            {
                return classicGridData;
            }
            else if (gridData is Foundation.Grid.Aggregation.AggregationGridData aggGridData)
            {
                Foundation.Grid.Classic.GridData data = aggGridData.AncestorGrid;
                return aggGridData.AncestorGrid;
            }
            else
            {
                throw new NotImplementedException("ToDo: Grid Type not yet implemented");
            }
        }

    }
}
