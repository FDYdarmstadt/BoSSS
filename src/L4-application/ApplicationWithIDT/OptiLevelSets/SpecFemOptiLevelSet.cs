using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using ilPSP;
using System;
using BoSSS.Solution.Utils;
using BoSSS.Solution.Statistic;
using BoSSS.Foundation.SpecFEM;
using System.Linq;
using ilPSP.LinSolvers;

namespace ApplicationWithIDT.OptiLevelSets {
    /// <summary>
    /// Spectral Finite Element Level Set:
    /// - defined on a cartesian grid -> must have at least degree 2
    /// - is continuous
    /// </summary>
    public class SpecFemOptiLevelSet : SpecFemField, IOptiLevelSet {
        /// <summary>
        /// the Transformation matrix between LsTBO and this one
        /// </summary>
        public MultidimensionalArray m_TransMat;

        /// <summary>
        /// This ist the LevelSet which is used internally by the method to compute Residuals etc.
        /// </summary>
        public LevelSet m_LevelSet;

        // this is a DG_Field with the same Basis as base.ContainingDGBasis
        public SinglePhaseField m_DGField;
        /// <summary>
        /// LSTracker for the LevelSet used by DGCode
        /// </summary>
        public LevelSetTracker m_Tracker;
        /// <summary>
        /// LSTracker for the OptiLevelSet
        /// </summary>
        public LevelSetTracker m_thisTracker;
        /// <summary>
        /// Number of the LevelSet in LSTRacker logic(usually 0 or 1)
        /// </summary>
        public int m_LevelSetNumber;
        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="__Basis"></param>
        /// <param name="LsTBO"></param>
        /// <param name="withAssembly"></param>
        /// <param name="trk"></param>
        /// <param name="num"></param>
        public SpecFemOptiLevelSet(SpecFemBasis __Basis, LevelSet LsTBO, bool withAssembly, LevelSetTracker trk, int num) : base(__Basis) {
            m_DGField = new SinglePhaseField(Basis.ContainingDGBasis);

            m_LevelSet = LsTBO;
            if(withAssembly) {
                AssembleTransMat(LsTBO);
            }
            m_Tracker = trk;
            m_LevelSetNumber = num;
        }
        /// <summary>
        /// Return this Grid
        /// </summary>
        /// <returns></returns>
        public IGridData GetGrid() {
            return Basis.GridDat;
        }

        /// <summary>
        /// Clones this Field
        /// </summary>
        /// <returns></returns>
        public SpecFemOptiLevelSet CloneAs() {
            return (SpecFemOptiLevelSet)Clone();
        }
        /// <summary>
        /// Accumalates onto a coordinate
        /// </summary>
        /// <param name="index"> index of that coordinate</param>
        /// <param name="acc">value to be accumlated</param>
        public void AccToParam(int index, double acc) {
            Coordinates[index] += acc;
        }

        /// <summary>
        /// Clears all the Coordinates for this Field
        /// </summary>
        public void Clear() {
            Coordinates.Clear();
        }

        /// <summary>
        /// method that assembles a matrix that Transforms this basis onto Basis of targetLS
        /// </summary>
        /// <param name="targetLS"></param>
        public void AssembleTransMat(LevelSet targetLS) {
            m_TransMat = MultidimensionalArray.Create(targetLS.CoordinateVector.Length, GetLength());
            SpecFemOptiLevelSet tmp = CloneAs();
            LevelSet LSbackup = targetLS.CloneAs();
            //Set all params of tmp to zero
            tmp.Clear();
            for(int i = 0; i < GetLength(); i++) {
                // set one Coordinate to 1
                tmp.Coordinates[i] = 1;

                // project onto m_DGField
                m_DGField.Clear();
                AccToDGField(1.0, m_DGField);

                // project onto LSbackup
                LSbackup.Clear();
                LSbackup.ProjectFromForeignGrid(1.0, m_DGField);

                //then the coordinates of LSBackup are equlling the entries of the general transformation matrix (bc of orthonormality).
                for(int j = 0; j < LSbackup.CoordinateVector.Length; j++) {
                    m_TransMat[j, i] = LSbackup.CoordinateVector[j];
                }
                // reset tmp to zero
                tmp.Coordinates[i] = 0;
            }
        }

        public void CopyParamsFrom(IOptiLevelSet source) {
            if(GetLength() == source.GetLength()) {
                for(int i = 0; i < GetLength(); i++) {
                    Coordinates[i] = source.GetParam(i);
                }
            } else {
                throw new ArgumentException("OptiLevelSet.CopyParamsFrom: source length, " + source.GetLength() + "  differs from target length: " + GetLength());
            }
        }

        public int GetLength() {
            return Coordinates.Length;
        }

        public double GetParam(int index) {
            return Coordinates[index];
        }

        public string GetParamName(int index) {
            return index.ToString();
        }

        public double[] GetParamsAsArray() {
            return Coordinates.To1DArray();
        }

        public void Print() {
            for(int i = 0; i < GetLength(); i++) {
                Console.Write(Coordinates[i] + " ");
            }
        }

        public void ProjectFromForeignLevelSet(SinglePhaseField sourceLS) {
            //// unneeded?
            //m_DGField.Clear();
            //m_DGField.ProjectFromForeignGrid(1.0, sourceLS);
            //base.ProjectDGField(1.0, m_DGField);
            Clear();
            ProjectDGField(1.0, sourceLS);
        }

        public void ProjectFromFunction(Func<double[], double> initialShockPostion) {
            Clear();
            m_DGField.ProjectField(initialShockPostion);
            ProjectDGField(1.0, m_DGField);
        }

        public void ProjectFromLevelSet(ConventionalDGField sourceLS) {
            Clear();
            ProjectDGField(1.0, sourceLS);
        }

        public void ProjectOntoLevelSet(LevelSet targetLS) {
            targetLS.Clear();
            m_DGField.Clear();
            AccToDGField(1.0, m_DGField);
            targetLS.ProjectFromForeignGrid(1.0, m_DGField);
        }
        /// <summary>
        /// Set a param
        /// </summary>
        /// <param name="index"></param>
        /// <param name="val"></param>
        public void SetParam(int index, double val) {
            Coordinates[index] = val;
        }
        /// <summary>
        /// indicates if this basis is orthonormal
        /// </summary>
        /// <returns></returns>
        public bool TestOrthonormality() {
            return false;
        }
        /// <summary>
        /// Converts to a SinglePhaseField of varaible Degree
        /// </summary>
        /// <param name="degree"></param>
        /// <returns></returns>
        public SinglePhaseField ToSinglePhaseField(int degree) {
            SinglePhaseField new_field;
            var new_basis = new Basis(Basis.GridDat, degree);
            new_field = new SinglePhaseField(new_basis);
            m_DGField.Clear();
            AccToDGField(1.0, m_DGField);
            new_field.ProjectFromForeignGrid(1.0, m_DGField);

            return new_field;
        }

        public LevelSet ToLevelSet(int degree) {
            LevelSet new_field;
            var new_basis = new Basis(Basis.GridDat, degree);
            new_field = new LevelSet(new_basis, "SpecFemLevelSetCopy");
            m_DGField.Clear();
            AccToDGField(1.0, m_DGField);
            new_field.ProjectFromForeignGrid(1.0, m_DGField);
            return new_field;
        }

        public int GetDegree() {
            return Basis.MaxDegree;
        }
        /// <summary>
        /// Get the mean Value
        /// </summary>
        /// <param name="i"></param>
        /// <returns></returns>
        double IOptiLevelSet.GetMeanValue(int i) {
            m_DGField.Clear();
            AccToDGField(1.0, m_DGField);
            return m_DGField.GetMeanValue(i);
        }

        /// <summary>
        /// Right now I dont have any Idea how to check if the coordinate influences the residual. 
        /// </summary>
        /// <param name="index"></param>
        /// <returns></returns>
        public bool IsInNearBand(int index) {
            return true;
        }
        /// <summary>
        /// Assembles a Tracker for this field
        /// </summary>
        public void AssembleTracker() {
            var dummy_LS = new LevelSet(Basis.ContainingDGBasis, "dummy_LS");
            m_thisTracker = new LevelSetTracker(Basis.GridDat, m_Tracker.CutCellQuadratureType, GetDegree(), m_Tracker.GetSpeciesSeparatedByLevSet(m_LevelSetNumber).ToArray(), dummy_LS);
        }

        public object Clone() {
            var r = new SpecFemOptiLevelSet(this.Basis, m_LevelSet, false, m_Tracker, m_LevelSetNumber);
            r.Coordinates.Acc(1.0, Coordinates);
            r.m_DGField = new SinglePhaseField(r.Basis.ContainingDGBasis);
            //m_TransMat can be shared as it is the same

            r.m_TransMat = m_TransMat;
            return r;
        }

        public double Norm(double[] levelSetStepCoordinates) {
            //Old version of calculation of norm (works better for ECCOMAS Test Case Scalar Advection)
            //double norm = 0;
            //for(int i = 0; i < this.GetLength(); i++) {
            //    norm += levelSetStepCoordinates[i];
            //}
            return levelSetStepCoordinates.L2Norm();
        }

        // <summary>
        /// here one could do interface straightening by projection onto p=1? but how to identify an oscilating interface?
        /// </summary>
        public void Reinitialize(double L, double kappa_S) {

        }

        public MsrMatrix GetRegMatrix() {
            var length_l = this.GetLength();
            var reg = new MsrMatrix(length_l, length_l, 1, 1);
            IMutuableMatrixEx_Extensions.AccEyeSp(reg, 1.0);
            return reg;
        }
    }
}

