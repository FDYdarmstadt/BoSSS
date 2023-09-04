using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using ilPSP;
using System;
using System.Collections.Generic;
using BoSSS.Solution.Utils;
using System.Collections;
using BoSSS.Solution.Statistic;
using System.Runtime.CompilerServices;
using BoSSS.Solution.LevelSetTools;
using static BoSSS.Solution.GridImport.NASTRAN.NastranFile;
using BoSSS.Platform;
using System.Linq;
using ilPSP.LinSolvers;

namespace ApplicationWithIDT.OptiLevelSets {
    /// <summary>
    /// Standard DG Field as Level Set
    /// </summary>
    public class SinglePhaseFieldOptiLevelSet : LevelSet, IOptiLevelSet {
        /// <summary>
        /// the Transformation matrix between LsTBO and this one
        /// </summary>
        public MultidimensionalArray m_TransMat;
        /// <summary>
        /// LevelSet used by DGCode
        /// </summary>
        public LevelSet m_LevelSet;
        /// <summary>
        /// LSTracker for the LevelSet used by DGCode
        /// </summary>
        public LevelSetTracker m_Tracker;
        /// <summary>
        /// LSTracker for the OptiLevelSet
        /// </summary>
        public LevelSetTracker m_thisTracker;
        public int m_LevelSetNumber;

        public IGridData GetGrid() {
            return GridDat;
        }

        public SinglePhaseFieldOptiLevelSet(Basis __Basis, string __Identification, LevelSet LsTBO, LevelSetTracker trk, int num) : base(__Basis, __Identification) {
            AssembleTransMat(LsTBO);
            m_LevelSet = LsTBO;
            m_Tracker = trk;
            m_LevelSetNumber = num;
        }
        /// <summary>
        /// cannot be done in the Constructor as this gives an infinite Loop ^^
        /// </summary>
        public void AssembleTracker() {
            var pair = m_Tracker.GetSpeciesSeparatedByLevSet(m_LevelSetNumber);
            m_thisTracker = new LevelSetTracker((GridData)Basis.GridDat, m_Tracker.CutCellQuadratureType, 1, pair.ToArray(), this);
        }
        public override LevelSet CloneAs() {
            return (SinglePhaseFieldOptiLevelSet)base.CloneAs();
        }
        public override SinglePhaseFieldOptiLevelSet Clone() {
            var r = new SinglePhaseFieldOptiLevelSet(Basis, Identification, m_LevelSet, m_Tracker, m_LevelSetNumber);
            r.CoordinateVector.Acc(1.0, CoordinateVector);
            return r;
        }
        public void AccToParam(int index, double acc) {
            CoordinateVector[index] += acc;
        }

        public void AssembleTransMat(LevelSet targetLS) {
            //m_TransMat = MultidimensionalArray.Create(targetLS.CoordinateVector.Length, this.GetLength());
            //SinglePhaseField mp = this.CloneAs();
            //LevelSet LSbackup = targetLS.CloneAs();
            ////Set all params of mp to zero
            //tmp.Clear();
            //for(int i = 0; i < this.GetLength(); i++) {
            //    tmp.CoordinateVector[i] = 1;
            //    LSbackup.Clear();
            //    LSbackup.ProjectFromForeignGrid(1.0,mp);
            //    for(int j = 0; j < LSbackup.CoordinateVector.Length; j++) {
            //        m_TransMat[j, i] = LSbackup.CoordinateVector[j];
            //    }
            //    tmp.CoordinateVector[i] = 0;
            //}
        }
        /// <summary>
        /// Checks if the param index corresponds to a cell in the Near-field of the MemberGrid (not the grid of the internal BOSSS LevelSet)
        /// </summary>
        /// <param name="index">i param index</param>
        /// <returns></returns>
        public bool IsInNearBand(int index) {
            //this gives us the Cell the index corresponds to
            Mapping.LocalFieldCoordinateIndex(index, out int iField, out int jCell, out int nMode);
            //then we get the Cell indices from the NearField
            var NFmask = m_thisTracker.Regions.GetNearFieldMask(jCell);
            //we check if the cell is in the near-field
            return NFmask.Contains(jCell);
        }
        public void CopyParamsFrom(IOptiLevelSet source) {
            if(GetLength() == source.GetLength()) {
                for(int i = 0; i < GetLength(); i++) {
                    CoordinateVector[i] = source.GetParam(i);
                }
            } else {
                throw new ArgumentException("OptiLevelSet.CopyParamsFrom: source length, " + source.GetLength() + "  differs from target length: " + GetLength());
            }
        }

        public int GetLength() {
            return CoordinateVector.Length;
        }

        public double GetParam(int index) {
            return CoordinateVector[index];
        }

        public string GetParamName(int index) {
            int j;
            int i;
            int n;
            Mapping.LocalFieldCoordinateIndex(index, out i, out j, out n);

            return "Field " + i + ", Cell " + j + ", Mode " + n;
        }

        public double[] GetParamsAsArray() {
            return CoordinateVector.ToArray();
        }

        public void Print() {
            for(int i = 0; i < GetLength(); i++) {
                Console.Write(CoordinateVector[i] + " ");
            }
        }

        public void ProjectFromForeignLevelSet(SinglePhaseField sourceLS) {
            Clear();
            this.ProjectFromForeignGrid(1.0, sourceLS);
        }

        public void ProjectFromFunction(Func<double[], double> initialShockPostion) {
            this.ProjectField(initialShockPostion);
        }

        public void ProjectFromLevelSet(ConventionalDGField sourceLS) {
            Clear();
            this.ProjectFromForeignGrid(1.0, sourceLS);
        }

        public void ProjectOntoLevelSet(LevelSet targetLS) {
            targetLS.Clear();
            targetLS.ProjectFromForeignGrid(1.0, CloneAs());
        }

        public void SetParam(int index, double val) {
            CoordinateVector[index] = val;
        }

        public bool TestOrthonormality() {
            return true;
        }

        public SinglePhaseField ToSinglePhaseField(int degree) {
            SinglePhaseField new_field;
            if(degree == GetDegree()) {
                new_field = new SinglePhaseField(Basis);
                new_field.CoordinateVector.Acc(1.0, CoordinateVector);
            } else {
                var new_basis = new Basis(GridDat, degree);
                new_field = new SinglePhaseField(new_basis);
                new_field.ProjectFromForeignGrid(1.0, this);
            }
            return new_field;
        }

        public LevelSet ToLevelSet(int degree) {
            LevelSet new_field;
            if(degree == GetDegree()) {
                new_field = new LevelSet(Basis, "OptiLevelSet");
                new_field.CoordinateVector.Acc(1.0, CoordinateVector);
            } else {
                var new_basis = new Basis(GridDat, degree);
                new_field = new LevelSet(new_basis, "OptiLevelSet");
                new_field.ProjectFromForeignGrid(1.0, this);
            }
            return new_field;
        }

        public int GetDegree() {
            return Basis.Degree;
        }

        double IOptiLevelSet.GetMeanValue(int i) {
            return GetMeanValue(i);
        }

        /// <summary>
        /// here one could do a continuity projection ??, potentially also interface straightening by projection onto p1? but how to identify an oscillating interface?
        /// </summary>
        public void Reinitialize(double L, double kappa_S) {

        }
        public double Norm(double[] levelSetStepCoordinates) {
            //Old version of calculation of norm (works better for ECCOMAS Test Case Scalar Advection)
            //double norm = 0;
            //for(int i = 0; i < this.GetLength(); i++) {
            //    norm += levelSetStepCoordinates[i];
            //}
            return levelSetStepCoordinates.L2Norm();
            
        }

        public MsrMatrix GetRegMatrix() {
            var length_l = this.GetLength();
            var reg = new MsrMatrix(length_l, length_l, 1, 1);
            IMutuableMatrixEx_Extensions.AccEyeSp(reg, 1.0);
            return reg;
        }
    }
}
