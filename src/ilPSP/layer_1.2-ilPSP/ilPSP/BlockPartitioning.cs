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

using ilPSP.Utils;
using MPI.Wrappers;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ilPSP {

    /// <summary>
    /// An implementation of <see cref="IBlockPartitioning"/>.
    /// </summary>
    public class BlockPartitioning : Partitioning, IBlockPartitioning {

        public BlockPartitioning(int LocalLength, IEnumerable<int> BlockI0, IEnumerable<int> BlockLen, MPI_Comm MpiComm, bool i0isLocal = false) : base(LocalLength, MpiComm) {

            //var enu_I0 = BlockI0.GetEnumerator();
            //var enuLen = BlockLen.GetEnumerator();

            Debug.Assert(base.LocalLength == LocalLength);

            List<int> _SubblkLen1 = new List<int>();
            List<int> _SubblkLen2 = new List<int>();
            List<int> _BlockType = new List<int>();

            //int NoOfBlocks = 0;
            //while(enu_I0.MoveNext()) {
            //    if (!enuLen.MoveNext())
            //        throw new ArgumentException("Both enumerations must contain the same number of elements.");

            //    int N = enuLen.Current;
            //    int i0Block = enu_I0.Current;

            //    if (i0isLocal)
            //        i0Block += base.i0;
            //    base.IsInLocalRange(i0Block);
            //    base.IsInLocalRange(Math.Max(i0Block, i0Block + N - 1));




            //    NoOfBlocks++;
            //}
            //if (enuLen.MoveNext())
            //    throw new ArgumentException("Both enumerations must contain the same number of elements.");


            int[] _BlockI0 = BlockI0.ToArray();
            int[] _BlockLen = BlockLen.ToArray();

            if(_BlockI0.Length != _BlockLen.Length)
                throw new ArgumentException("Both enumerations must contain the same number of elements.");

            int NoOfBlocks = _BlockI0.Length;
            int[] BlockType = new int[NoOfBlocks];

            int FrameBlockSize = -1;

            for(int iBlock = 0; iBlock < NoOfBlocks; iBlock++) {
                int i0Block = _BlockI0[iBlock];
                if (i0isLocal)
                    i0Block += base.i0;

                int N = _BlockLen[iBlock];

                base.IsInLocalRange(i0Block);
                base.IsInLocalRange(Math.Max(i0Block, i0Block + N - 1));

                int iEBlock;
                if( iBlock < NoOfBlocks - 1) {
                    iEBlock = _BlockI0[iBlock + 1];
                    if (i0isLocal)
                        iEBlock += base.i0;
                } else {
                    iEBlock = LocalLength + base.i0;
                }

                if (i0Block + N > iEBlock)
                    throw new ArgumentException("Block Length exceeds i0 of next block.");

                int NE = iEBlock - i0Block;

                if (iBlock == 0) {
                    FrameBlockSize = NE;
                } else {
                    if (NE != FrameBlockSize)
                        FrameBlockSize = -1;
                }

                int iBlockType;
                bool bFound = false;
                for(iBlockType = 0; iBlockType < _SubblkLen1.Count; iBlockType++) {
                    if(N == _SubblkLen1[iBlockType] && NE == _SubblkLen2[iBlockType]) {
                        bFound = true;
                        break;
                    }
                }
                Debug.Assert(bFound == (iBlockType < _SubblkLen1.Count));
                if(!bFound) {
                    _SubblkLen1.Add(N);
                    _SubblkLen2.Add(NE);
                    iBlockType = _SubblkLen1.Count - 1;
                }

                Debug.Assert(N == _SubblkLen1[iBlockType]);
                Debug.Assert(NE == _SubblkLen2[iBlockType]);
                BlockType[iBlock] = iBlockType;
            }

            //

            int NoOfBlockTypes = _SubblkLen1.Count;
            int[][] i0_Sblk = NoOfBlocks.ForLoop(iBlkType => new int[] { 0 });
            int[][] LenSblk = NoOfBlocks.ForLoop(iBlkType => new int[] { _SubblkLen1[iBlkType] });
            ConstructorCommon(FrameBlockSize, i0_Sblk, LenSblk, BlockType, MpiComm);

        }



        /// <summary>
        /// Shallow constructor - basically everything has to be provided
        /// </summary>
        /// <param name="LocalLength"></param>
        /// <param name="FrameBlockSize">
        /// The size of the frame block:
        /// - if negative, not used; in this case, the length of each block can be different and is defined by the sub-blocking, i.e. by <paramref name="_Subblk_i0"/>, <paramref name="_SubblkLen"/> and <paramref name="_BlockType"/>. In this case <see cref="AllBlockSizesEqual"/> is false.
        /// - if positive, a frame of constant size is assumed around each block, i.e. the index range of the i-th block starts at index i*<paramref name="FrameBlockSize"/>; then <see cref="AllBlockSizesEqual"/> is true. However, the sub-blocking can still differ from block to block.
        /// </param>
        /// <param name="_Subblk_i0">
        /// see <see cref="GetSubblk_i0(int)"/>.
        /// </param>
        /// <param name="_SubblkLen">
        /// see <see cref="GetSubblkLen(int)"/>.
        /// </param>
        /// <param name="_BlockType">
        /// see <see cref="GetBlockType(int)"/>
        /// </param>
        /// <param name="MpiComm"></param>
        public BlockPartitioning(
            int LocalLength, 
            int FrameBlockSize,
            int[][] _Subblk_i0, int[][] _SubblkLen,
            int[] _BlockType,
            MPI_Comm MpiComm) : base(LocalLength, MpiComm) {
            ConstructorCommon(FrameBlockSize, _Subblk_i0, _SubblkLen, _BlockType, MpiComm);
        }

        private void ConstructorCommon(int FrameBlockSize, int[][] _Subblk_i0, int[][] _SubblkLen, int[] _BlockType, MPI_Comm MpiComm) {
            MPICollectiveWatchDog.Watch(MpiComm);
            int LocalLength = base.LocalLength;

            // ===============
            // check arguments
            // ===============
            if (_Subblk_i0.Length != _SubblkLen.Length)
                throw new ArgumentException();

            for (int i = _Subblk_i0.Length - 1; i >= 0; i--) {
                if (_Subblk_i0[i].Length != _SubblkLen[i].Length)
                    throw new ArgumentException();
                int i0min = 0;
                for (int iSblk = 0; iSblk < _Subblk_i0[i].Length; iSblk++) {
                    if (_Subblk_i0[i][iSblk] < i0min)
                        throw new ArgumentException("Potential sub-block overlapping.");
                    i0min += _SubblkLen[i][iSblk];
                }
            }

            this.m_Subblk_i0 = _Subblk_i0;
            this.m_SubblkLen = _SubblkLen;

            if (FrameBlockSize >= 1) {
                for (int BlockType = 0; BlockType < _Subblk_i0.Length; BlockType++) {
                    int NoOfSubblocks = _Subblk_i0[BlockType].Length;
                    int BlockLen = _Subblk_i0[BlockType][NoOfSubblocks - 1] + _SubblkLen[BlockType][NoOfSubblocks - 1];

                    if (BlockLen > FrameBlockSize)
                        throw new ArgumentException();
                }

                if (LocalLength % FrameBlockSize != 0) {
                    throw new ArgumentException("'FrameBlockSize', if specified, must be a divider of 'LocalLength'.");
                }

                if (_BlockType != null) {
                    if (_BlockType.Length != LocalLength / FrameBlockSize)
                        throw new ArgumentException("Mismatch between number of blocks specified by '_BlockType' and 'LocalLength/FrameBlockSize'.");
                }
            }

#if DEBUG
            if (_BlockType != null) {
                for (int j = 0; j < _BlockType.Length; j++) {
                    if (_BlockType[j] < 0 || _BlockType[j] > _Subblk_i0.Length)
                        throw new ArgumentException();
                }
            }
#endif

            // ====================
            // set internal members
            // ====================

            int LocalNoOfBlocks;
            if (_BlockType != null) {
                LocalNoOfBlocks = _BlockType.Length;

                if (FrameBlockSize > 0) {
                    //throw new ArgumentException("If no BlockTape is given, the FrameBlockSize must be specified.");
                    if (LocalLength / FrameBlockSize != _BlockType.Length)
                        throw new ArgumentException("FrameBlockSize does not match _BlockType.");
                }
            } else {
                if (FrameBlockSize <= 0)
                    throw new ArgumentException("If no BlockTape is given, the FrameBlockSize must be specified.");
                LocalNoOfBlocks = LocalLength / FrameBlockSize;
            }

            m_BlocksPartition = new Partitioning(LocalNoOfBlocks, MpiComm);

            if (_Subblk_i0.Length > 1) {
                this.m_BlockType = _BlockType;
            } else {
                this.m_BlockType = null;
            }
            int J = _BlockType.Length;

            if (FrameBlockSize == 0) {
                throw new ArgumentException();
            } else if (FrameBlockSize < 0) {
                int[] BlockI0 = new int[J + 1];
                BlockI0[0] = base.i0;

                int FrameBlockSizeMaybe = -1;
                for (int j = 0; j < J; j++) {
                    int bT = _BlockType[j];
                    int B = _Subblk_i0[bT].Length;
                    int BlockLen = _Subblk_i0[bT][B - 1] + _SubblkLen[bT][B - 1];
                    BlockI0[j + 1] = BlockI0[j] + BlockLen;
                    if (j == 0) {
                        FrameBlockSizeMaybe = BlockLen;
                    } else {
                        if (FrameBlockSizeMaybe != BlockLen)
                            FrameBlockSizeMaybe = -1;
                    }
                }

                Debug.Assert(BlockI0[J] == base.LocalLength + base.i0);

                int FrameBlockSizeMaybeMin = FrameBlockSizeMaybe.MPIMin(MpiComm);
                int FrameBlockSizeMaybeMax = FrameBlockSizeMaybe.MPIMax(MpiComm);
                if (FrameBlockSizeMaybeMin > 0 && FrameBlockSizeMaybe == FrameBlockSizeMaybeMin && FrameBlockSizeMaybe == FrameBlockSizeMaybeMax) {
                    // constant blocking, although the stupid caller told us different
                    m_FrameBlockSize = FrameBlockSizeMaybe;
                    m_Block_i0 = null;
                } else {
                    m_Block_i0 = BlockI0;
                    m_FrameBlockSize = -1;
                }
            } else {
                Debug.Assert(FrameBlockSize >= 1);
                m_FrameBlockSize = FrameBlockSize;

            }
        }

        Partitioning m_BlocksPartition;


        int[][] m_SubblkLen;
#if DEBUG
        int[][] m_SubblkLen_Backup;
#endif

        /// <summary>
        /// Lengths of sub-blocks.
        /// </summary>
        public int[] GetSubblkLen(int blockType) {

#if DEBUG
            // make sure that no one outside modifies this member
            if (m_SubblkLen_Backup == null) {
                m_SubblkLen_Backup = m_SubblkLen.Select(bl => (int[])(bl.Clone())).ToArray();
            } else {
                Debug.Assert(m_SubblkLen.Length == m_SubblkLen_Backup.Length);
                for (int i = 0; i < m_SubblkLen.Length; i++) {
                    Debug.Assert(ArrayTools.ListEquals(m_SubblkLen[i], m_SubblkLen_Backup[i]), "Some moron messed with blocking indices.");
                }
            }
#endif
            return m_SubblkLen[blockType];

        }

        int[][] m_Subblk_i0;
#if DEBUG
        int[][] m_Subblki0_Backup;
#endif

        /// <summary>
        /// Start-indices of sub-blocks.
        /// </summary>
        public int[] GetSubblk_i0(int blockType) {

#if DEBUG
            // make sure that no one outside modifies this member
            if (m_Subblki0_Backup == null) {
                m_Subblki0_Backup = m_Subblk_i0.Select(bl => (int[])(bl.Clone())).ToArray();
            } else {
                Debug.Assert(m_Subblk_i0.Length == m_Subblki0_Backup.Length);
                for (int i = 0; i < m_SubblkLen.Length; i++) {
                    Debug.Assert(ArrayTools.ListEquals(m_Subblk_i0[i], m_Subblki0_Backup[i]), "Some moron messed with blocking indices.");
                }
            }
#endif
            return m_Subblk_i0[blockType];

        }

 

        int[] m_BlockType;
        int[] m_Block_i0;

        /// <summary>
        /// Block type.
        /// </summary>
        /// <param name="iBlock">
        /// Block index.
        /// </param>
        /// <returns>
        /// The block type, 
        /// an index into <see cref="BlockLen"/>, <see cref="Subblk_i0"/>, <see cref="SubblkLen"/>.
        /// </returns>
        public int GetBlockType(int iBlock) {
            m_BlocksPartition.TestIfInLocalRange(iBlock);
            if(m_BlockType == null) {
                return 0;
            } else {
                return m_BlockType[iBlock - m_BlocksPartition.i0];
            }
        }


        int m_FrameBlockSize = 0;

        /// <summary>
        /// Index offset for some block.
        /// </summary>
        /// <param name="iBlock">
        /// Block index.
        /// </param>
        /// <returns>
        /// A global index of the first element in the block.
        /// </returns>
        public int GetBlockI0(int iBlock) {
            m_BlocksPartition.TestIfInLocalRange(iBlock);
            Debug.Assert(m_FrameBlockSize != 0);
            if (m_FrameBlockSize > 0) {
                Debug.Assert(m_Block_i0 == null);
                Debug.Assert(m_FrameBlockSize * m_BlocksPartition.i0 == this.i0);
                return iBlock * m_FrameBlockSize;
            } else {
                return m_Block_i0[iBlock - m_BlocksPartition.i0];
            }
            //return base.BlockSize * iBlock;
        }

        /// <summary>
        /// Terminal index offset for some block.
        /// </summary>
        /// <param name="iBlock">
        /// Block index.
        /// </param>
        /// <returns>
        /// A global index of the last element in the block, *plus one*.
        /// </returns>
        public int GetBlockLen(int iBlock) {
            m_BlocksPartition.TestIfInLocalRange(iBlock);
            Debug.Assert(m_FrameBlockSize != 0);
            if (m_FrameBlockSize > 0) {
                Debug.Assert(m_Block_i0 == null);
                return m_FrameBlockSize;
            } else {
                int iBlockLoc = iBlock - m_BlocksPartition.i0;
                return m_Block_i0[iBlockLoc + 1] - m_Block_i0[iBlockLoc];
            }
        }

        /// <summary>
        /// Converts an index into a global index.
        /// </summary>
        /// <param name="i">
        /// Global index.
        /// </param>
        /// <returns>
        /// Global block index.
        /// </returns>
        public int GetBlockIndex(int i) {
            this.TestIfInLocalRange(i);
            Debug.Assert(m_FrameBlockSize != 0);
            if (m_FrameBlockSize > 0) {
                Debug.Assert(m_Block_i0 == null);
                Debug.Assert(m_FrameBlockSize * m_BlocksPartition.i0 == this.i0);

                return i / m_FrameBlockSize;
            } else {

                int iBlockLoc = Array.BinarySearch<int>(m_Block_i0, i);
                if (iBlockLoc < 0) {
                    iBlockLoc = (~iBlockLoc) - 1;
                }

                // Play it safe in case of potentially empty blocks
                Debug.Assert(this.m_Block_i0.Length - 1 == m_BlocksPartition.LocalLength);
                while (m_Block_i0[iBlockLoc] == m_Block_i0[iBlockLoc + 1] && iBlockLoc < this.m_Block_i0.Length - 2) {
                    iBlockLoc++;
                }

                Debug.Assert(iBlockLoc < m_BlocksPartition.LocalLength);
                Debug.Assert(i >= m_Block_i0[iBlockLoc]);
                Debug.Assert(i < m_Block_i0[iBlockLoc + 1]);
                return iBlockLoc + m_BlocksPartition.i0;
            }
        }

        virtual public IBlockPartitioning GetImmutableBlockPartitioning() {
            if (this.IsMutable) {
                return new BlockPartitioning(this.LocalLength, m_FrameBlockSize,
                    this.m_Subblk_i0.Select(a => (int[])a.Clone()).ToArray(), this.m_SubblkLen.Select(a => (int[])a.Clone()).ToArray(),
                    (int[])this.m_Block_i0.Clone(),
                    this.MPI_Comm);
            } else {
                return this;
            }
        }

        public int GetFirstBlock(int proc) {
            return m_BlocksPartition.GetI0Offest(proc);
        }

        public int GetLocalNoOfBlocks(int proc) {
            return m_BlocksPartition.GetLocalLength(proc);
        }

        public int FindProcessForBlock(int iBlk) {
            return m_BlocksPartition.FindProcess(iBlk);
        }


        /// <summary>
        /// Total number of blocks (over all MPI processors).
        /// </summary>
        public int TotalNoOfBlocks {
            get {
                return m_BlocksPartition.TotalLength;
            }
        }

        /// <summary>
        /// Local number of blocks (on this MPI process).
        /// </summary>
        public int LocalNoOfBlocks {
            get {
                Debug.Assert(this.m_BlockType == null || this.m_BlockType.Length == m_BlocksPartition.LocalLength);
                return m_BlocksPartition.LocalLength;
            }
        }

        /// <summary>
        /// Global block index of first block on this MPI process.
        /// </summary>
        public int FirstBlock  {
            get {
                return this.m_BlocksPartition.i0;
            }
        }

        /// <summary>
        /// whether the blocking can change over time, or not.
        /// </summary>
        public override bool IsMutable {
            get {
                return false;
            }
        }

        public bool AllBlockSizesEqual {
            get {
                Debug.Assert(m_FrameBlockSize != 0);
                return m_FrameBlockSize > 0;
            }
        }
    }
}
