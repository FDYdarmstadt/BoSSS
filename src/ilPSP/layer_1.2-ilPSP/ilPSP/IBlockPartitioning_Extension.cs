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
    /// Extension methods.
    /// </summary>
    public static class IBlockPartitioning_Extension {

        /// <summary>
        /// Gets a list of all occupied indices in the blocking.
        /// </summary>
        public static int[] GetOccupiedIndicesList(this IBlockPartitioning B) {
            List<int> R = new List<int>();

            for(int iBlk = B.FirstBlock; iBlk < B.LocalNoOfBlocks + B.FirstBlock; iBlk++) {
                int i0 = B.GetBlockI0(iBlk);

                int bt = B.GetBlockType(iBlk);
                int[] S_i0 = B.GetSubblk_i0(bt);
                int[] SLen = B.GetSubblkLen(bt);
                Debug.Assert(S_i0.Length == SLen.Length);
                int NoSbl = S_i0.Length;

                for(int iSbl = 0; iSbl < NoSbl; iSbl++) {
                    int si0 = S_i0[iSbl];
                    int siE = si0 + SLen[iSbl];
                    for (int i = si0; i < siE; i++) {
                        R.Add(i + i0);
                    }
                }
            }

            return R.ToArray();
        }



        /// <summary>
        /// Constructs a block partitioning from a sub-vector index list.
        /// </summary>
        /// <param name="B"></param>
        /// <param name="SubvectorIdx">
        /// Indices into <paramref name="B"/> that should be in the result.
        /// </param>
        /// <param name="comm"></param>
        /// <param name="FrameBlockSize">
        /// <see cref="IBlockPartitioning.AllBlockSizesEqual"/>;
        /// The size of the frame block:
        /// - if negative, not used; in this case, the length of each block can be different and is defined by the sub-blocking, i.e. by <paramref name="_Subblk_i0"/>, <paramref name="_SubblkLen"/> and <paramref name="_BlockType"/>. In this case <see cref="AllBlockSizesEqual"/> is false.
        /// - if positive, a frame of constant size is assumed around each block, i.e. the index range of the i-th block starts at index i*<paramref name="FrameBlockSize"/>; then <see cref="AllBlockSizesEqual"/> is true. However, the sub-blocking can still differ from block to block.
        /// </param>
        /// <returns></returns>
        public static BlockPartitioning GetSubBlocking<T>(this IBlockPartitioning B, T SubvectorIdx, MPI_Comm comm, int FrameBlockSize = -1) 
            where T : IEnumerable<int>
        {
            List<int[]> Subblk_i0 = new List<int[]>();
            List<int[]> SubblkLen = new List<int[]>();

            List<int> BlockType = new List<int>();

            int iCheck = -1;

            int iCurrentBlock = -1,  CurrentBlockType = -12, B_i0 = -1, BLen = -1;
            int[] B_sblk_i0 = null;
            int[] B_sblkLen = null;
            bool[] SblkMarker = null;
            int LocalLength = 0;
            List<int> tmp_i0 = new List<int>();
            List<int> tmpLen = new List<int>();

            int FrameBlkAuto = 0;
            foreach (int i in SubvectorIdx) {
                B.TestIfInLocalRange(i);
                if (i < iCheck)
                    throw new ArgumentException("Sub-vector indices must be sorted.");
                iCheck = i;

                int iBlock = B.GetBlockIndex(i);
                if(iBlock != iCurrentBlock) {
                    // block change - need to do something

                    // record actual block
                    if (iCurrentBlock >= 0) {
                        RecordBlock(Subblk_i0, SubblkLen, BlockType, B_sblk_i0, B_sblkLen, SblkMarker, tmp_i0, tmpLen, ref FrameBlkAuto);
                    }

                    // move on to next block
                    iCurrentBlock = iBlock;
                    CurrentBlockType = B.GetBlockType(iBlock);
                    B_sblk_i0 = B.GetSubblk_i0(CurrentBlockType);
                    B_sblkLen = B.GetSubblkLen(CurrentBlockType);
                    B_i0 = B.GetBlockI0(iBlock);
                    BLen = B.GetBlockLen(BLen);
                    if (SblkMarker == null || SblkMarker.Length < BLen)
                        SblkMarker = new bool[BLen];
                    else
                        Array.Clear(SblkMarker, 0, SblkMarker.Length);
                }

                int iWithin = i - B_i0;
                Debug.Assert(iWithin >= 0);
                Debug.Assert(iWithin < BLen);

                SblkMarker[iWithin] = true;
                LocalLength++;
            }

            // record actual block
            if (iCurrentBlock >= 0) {
                RecordBlock(Subblk_i0, SubblkLen, BlockType, B_sblk_i0, B_sblkLen, SblkMarker, tmp_i0, tmpLen, ref FrameBlkAuto);
            }

            Debug.Assert(LocalLength == SubvectorIdx.Count());
            Debug.Assert(BlockType.Count == B.LocalNoOfBlocks);



            // finalize data structure
            // -----------------------
            if (FrameBlockSize == 0)
                FrameBlockSize = FrameBlkAuto;
            if (FrameBlockSize > 0)
                LocalLength = FrameBlockSize * BlockType.Count;
            return new BlockPartitioning(LocalLength, FrameBlockSize, 
                Subblk_i0.ToArray(), SubblkLen.ToArray(), BlockType.ToArray(), comm);
        }

        private static void RecordBlock(List<int[]> Subblk_i0, List<int[]> SubblkLen, List<int> BlockType, int[] B_sblk_i0, int[] B_sblkLen, bool[] SblkMarker, List<int> tmp_i0, List<int> tmpLen, ref int MaxBlockSize) {
            int NoOfSblk = 0;
            Debug.Assert(B_sblk_i0.Length == B_sblkLen.Length);
            int BNoOfSblk = B_sblk_i0.Length;

            // build sub-blocking
            tmp_i0.Clear();
            tmpLen.Clear();

            int i0New = 0;
            for (int bBlk = 0; bBlk < BNoOfSblk; bBlk++) { // loop over the original sub-blocks
                int iESblk = B_sblk_i0[bBlk] + B_sblkLen[bBlk];

                int _i0New = i0New;
                int LenNew = 0;
                for (int iSblk = B_sblk_i0[bBlk]; iSblk < iESblk; iSblk++) { // loop over indices in original sub-block...
                    if (SblkMarker[iSblk]) {
                        // index is also in the new block partition
                        LenNew++;
                    }
                }
                if (LenNew > 0) {
                    tmp_i0.Add(_i0New);
                    tmpLen.Add(LenNew);
                    NoOfSblk++;
                }

                i0New += LenNew;
            }

            MaxBlockSize = Math.Max(MaxBlockSize, tmp_i0[NoOfSblk - 1] + tmpLen[NoOfSblk - 1]);

            // record block type
            Debug.Assert(Subblk_i0.Count == SubblkLen.Count);
            int iBlockType = -1, L = Subblk_i0.Count;
            for (int l = 0; l < L; l++) {
                int[] _i0 = Subblk_i0[l];
                int[] Len = SubblkLen[l];
                Debug.Assert(_i0.Length == Len.Length);
                if (NoOfSblk != _i0.Length)
                    continue;
                bool bFound = true;
                for (int k = 0; k < NoOfSblk; k++) {
                    if (tmp_i0[k] != _i0[k]) {
                        bFound = false;
                        break;
                    }
                    if (tmpLen[k] != Len[k]) {
                        bFound = false;
                        break;
                    }
                }

                if (bFound) {
                    iBlockType = l;
                    break;
                }
            }
            if (iBlockType < 0) {
                Subblk_i0.Add(tmp_i0.ToArray());
                SubblkLen.Add(tmpLen.ToArray());
                iBlockType = Subblk_i0.Count - 1;
            }
            BlockType.Add(iBlockType);
        }
    }
}
