using BoSSS.Foundation;
using ilPSP;
using ilPSP.Connectors;
using ilPSP.LinSolvers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Foundation.Grid;

namespace BoSSS.Application.ExternalBinding {
    public class OpenFoamMatrix : BlockMsrMatrix, IForeignLanguageProxy {

        public OpenFoamDGField[] Fields;

        /// <summary>
        /// Constructor for quadratic matrices
        /// </summary>
        [CodeGenExport]
        public OpenFoamMatrix(OpenFOAMGrid grd, OpenFoamDGField f) :
            base(f.Mapping, f.Mapping) //
        {
            Fields = new[]{f};
            RowMap = f.Mapping;
            ColMap = f.Mapping;
            m_SolBuffer = f;
        }

        /// <summary>
        /// 
        /// </summary>
        public UnsetteledCoordinateMapping RowMap {
            get;
            private set;
        }

        /// <summary>
        /// 
        /// </summary>
        public UnsetteledCoordinateMapping ColMap {
            get;
            private set;
        }


        double[] m_RHSbuffer;

        /// <summary>
        /// additional buffer to store the RHS of the system
        /// </summary>
        public double[] RHSbuffer {
            get {
                if(m_RHSbuffer == null)
                    m_RHSbuffer = new double[_RowPartitioning.LocalLength];
                return m_RHSbuffer;
            }
        }

        /// <summary>
        /// Returns an entry of the RHS
        /// </summary>
        /// <param name="f">field/basis index</param>
        /// <param name="j">cell index</param>
        /// <param name="n">mode index</param>
        /// <returns></returns>
        [CodeGenExport]
        public double GetRHScoordinate(int f, int j, int n) {
            return RHSbuffer[this.RowMap.LocalUniqueCoordinateIndex(f, j, n)];
        }

        OpenFoamDGField m_SolBuffer;

        /// <summary>
        /// additional buffer to store the Solution of the system
        /// </summary>
        public OpenFoamDGField SolBuffer {
            get {
                return m_SolBuffer;
            }
        }

        /// <summary>
        /// Returns an entry of the Solution; <see cref="Solve"/> must be called first.
        /// </summary>
        /// <param name="f">field/basis index</param>
        /// <param name="j">cell index</param>
        /// <param name="n">mode index</param>
        /// <returns></returns>
        [CodeGenExport]
        public double GetSolCoordinate(int f, int j, int n) {
            // return SolBuffer[f][this.RowMap.LocalUniqueCoordinateIndex(f, j, n)];
            return SolBuffer[this.RowMap.LocalUniqueCoordinateIndex(0, j, n)];
        }


        /// <summary>
        /// Solve!
        /// </summary>
        /// <remarks>
        /// At this point, only the direct solver is used; later, this should be upgraded to use the advanced solver framework.
        /// </remarks>
        [CodeGenExport]
        public void Solve() {
            try {
                this.Solve_Direct(this.SolBuffer, RHSbuffer);
            } catch(Exception e) {
                Console.WriteLine(e.GetType().FullName + ": " + e.Message);
            }
        }




        /// <summary>
        /// Number of occupied rows in a block
        /// </summary>
        /// <param name="iLocal">local (within current MPI process) cell index/Block index</param>
        /// <returns></returns>
        [CodeGenExport]
        public int GetNoOfRowsInBlock(int iLocal) {
            long iGlob = _RowPartitioning.i0 + iLocal;
            _RowPartitioning.TestIfInLocalRange(iGlob);
            //_RowPartitioning.GetBlockLen(); // may also include un-used rows, e.g. for XDG

            // determine the index of the **last occupied** row 
            int Bt = _RowPartitioning.GetBlockType(iGlob);
            int[] i0S = _RowPartitioning.GetSubblk_i0(Bt);
            int[] Lns = _RowPartitioning.GetSubblkLen(Bt);
            Debug.Assert(i0S.Length == Lns.Length);
            int ke = i0S.Length - 1;

            int sz = i0S[ke] + Lns[ke];
            return sz;
        }

        /// <summary>
        /// Number of occupied columns in a block
        /// </summary>
        /// <param name="iLocal">local (within current MPI process) cell index/Block index</param>
        /// <returns></returns>
        [CodeGenExport]
        public int GetNoOfColsInBlock(int iLocal) {
            long iGlob = _ColPartitioning.i0 + iLocal;
            _ColPartitioning.TestIfInLocalRange(iGlob);
            //_RowPartitioning.GetBlockLen(); // may also include un-used rows, e.g. for XDG

            // determine the index of the **last occupied** row 
            int Bt = _ColPartitioning.GetBlockType(iGlob);
            int[] i0S = _ColPartitioning.GetSubblk_i0(Bt);
            int[] Lns = _ColPartitioning.GetSubblkLen(Bt);
            Debug.Assert(i0S.Length == Lns.Length);
            int ke = i0S.Length - 1;

            int sz = i0S[ke] + Lns[ke];
            return sz;
        }

        /// <summary>
        /// Read block from matrix;
        /// size is determined by <see cref="GetNoOfRowsInBlock(int)"/> <see cref="GetNoOfColsInBlock(int)"/>
        /// </summary>
        /// <param name="iRowLocal"></param>
        /// <param name="jColLocal"></param>
        /// <param name="InputReadBuffer">
        /// A buffer of <see cref="GetNoOfRowsInBlock(int)"/>*<see cref="GetNoOfColsInBlock(int)"/> double values,
        /// arranged in row major order (row index cycling fastest).
        /// </param>
        [CodeGenExport]
        unsafe public void GetBlock(int iRowLocal, int jColLocal, double* InputReadBuffer) {
            long iRowGlob = _RowPartitioning.i0 + iRowLocal;
            long jColGlob = _ColPartitioning.i0 + jColLocal;

            int N = GetNoOfRowsInBlock(iRowLocal);
            int M = GetNoOfColsInBlock(jColLocal);
            if(M <= 0 || N <= 0)
                return; // nothing to do;

            var Tmp = MultidimensionalArray.Create(N, M);

            this.ReadBlock(iRowGlob, jColGlob, Tmp);

            unsafe {
                fixed(double* pS = Tmp.Storage) {
                    // check the memory layout
                    Debug.Assert(Tmp.Index(0, 0) == 0);
                    Debug.Assert(Tmp.Index(0, 1) == 1);

                    double* pDst = InputReadBuffer;
                    double* pSrc = pS;
                    for(int cnt = N*M; cnt > 0; cnt--) {
                        *pDst = *pSrc;
                        pSrc++;
                        pDst++;
                    }
                }
            }

        }

        /// <summary>
        /// Accumulates a block to this matrix;
        /// size is determined by <see cref="GetNoOfRowsInBlock(int)"/> <see cref="GetNoOfColsInBlock(int)"/>
        /// </summary>
        /// <param name="iRowLocal"></param>
        /// <param name="jColLocal"></param>
        /// <param name="InputReadBuffer">
        /// A buffer of <see cref="GetNoOfRowsInBlock(int)"/>*<see cref="GetNoOfColsInBlock(int)"/> double values,
        /// arranged in row major order (row index cycling fastest).
        /// </param>
        [CodeGenExport]
        unsafe public void AccBlock(int iRowLocal, int jColLocal, double alpha, double* InputReadBuffer) {
            long iRowGlob = _RowPartitioning.i0 + iRowLocal;
            long jColGlob = _ColPartitioning.i0 + jColLocal;

            int N = GetNoOfRowsInBlock(iRowLocal);
            int M = GetNoOfColsInBlock(jColLocal);
            if(M <= 0 || N <= 0)
                return; // nothing to do;

            var Tmp = MultidimensionalArray.Create(N, M);

            unsafe {
                fixed(double* pS = Tmp.Storage) {
                    // check the memory layout
                    Debug.Assert(Tmp.Index(0, 0) == 0);
                    Debug.Assert(Tmp.Index(0, 1) == 1);

                    double* pSrc = InputReadBuffer;
                    double* pDst = pS;
                    for(int cnt = N*M; cnt > 0; cnt--) {
                        *pDst = *pSrc;
                        pSrc++;
                        pDst++;
                    }
                }
            }

            this.AccBlock(iRowGlob, jColGlob, 1.0, Tmp);
        }

        /// <summary>
        /// Sets an entire block to zero;
        /// size is determined by <see cref="GetNoOfRowsInBlock(int)"/> <see cref="GetNoOfColsInBlock(int)"/>
        /// </summary>
        /// <param name="iRowLocal"></param>
        /// <param name="jColLocal"></param>
        [CodeGenExport]
        unsafe public void ClearBlock(int iRowLocal, int jColLocal) {
            long iRowGlob = _RowPartitioning.i0 + iRowLocal;
            long jColGlob = _ColPartitioning.i0 + jColLocal;

            int N = GetNoOfRowsInBlock(iRowLocal);
            int M = GetNoOfColsInBlock(jColLocal);
            if(M <= 0 || N <= 0)
                return; // nothing to do;

            

            base.ClearBlock(iRowGlob, jColGlob, N, M);
        }




        IntPtr m_ForeignPtr;

        /// <summary>
        /// %
        /// </summary>
        public void _SetForeignPointer(IntPtr ptr) {
            if(ptr == IntPtr.Zero) {
                m_ForeignPtr = IntPtr.Zero;
            } else {

                if(m_ForeignPtr != IntPtr.Zero) {
                    throw new ApplicationException("already registered");
                }
                m_ForeignPtr = ptr;
            }
        }

        /// <summary>
        /// %
        /// </summary>
        public IntPtr _GetForeignPointer() {
            return m_ForeignPtr;
        }

        // public override string ToString(){
        //     var ret = new System.String("");
        //     ret += "OpenFoamMatrix with " + this.NoOfRows + " rows and " + this.NoOfCols + " columns\n";
        //     for (int i = 0; i < this.NoOfRows; i++){
        //         for (int j = 0; j < this.NoOfCols; j++){
        //             ret += string.Format("{0:#.#e+0} ", this[i, j]);
        //         }
        //         ret += "\n";
        //     }
        //     return ret;
        // }

    }
}
