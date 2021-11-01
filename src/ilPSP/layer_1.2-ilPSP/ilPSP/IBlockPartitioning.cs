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
using System.Threading.Tasks;

namespace ilPSP {

    /// <summary>
    /// In addition to the partitioning of an index over a range of MPI processors, which is defined by <see cref="IPartitioning"/>,
    /// this interface defines *two* additional levels of partitioning:
    /// - into blocks, which can have different sizes.
    /// - into sub-blocks which can also vary from block to block.
    /// The general assumption behind this design is that the number of different block types
    /// (block length, sub-block partitioning) is rather small in comparison to the total number of blocks.
    /// This is usually the case for a DG method - a typical scenario would be local adaption of the polynomial degree:
    /// In this case, one would have maybe three or four different polynomial degrees. Each degree corresponds to a 
    /// different block type. Then each cell -- usually much more than hundred, up to millions -- corresponds to a block in this partition,
    /// and one of the three resp. four block types.
    /// </summary>
    public interface IBlockPartitioning : IPartitioning {

        /// <summary>
        /// Length of sub-blocks.
        /// </summary>
        /// <param name="blockType">
        /// The block type index, e.g. the return value of <see cref="GetBlockType(long)"/>.
        /// </param>
        /// <returns>
        /// Lengths of sub-blocks, for the respective block type.
        /// - index: sub-block index within the respective block type
        /// - content: length of sub-block
        /// </returns>
        int[] GetSubblkLen(int blockType);

        /// <summary>
        /// Start-indices of sub-blocks.
        /// </summary>
        /// <param name="blockType">
        /// The block type index, e.g. the return value of <see cref="GetBlockType"/>.
        /// </param>
        /// <returns>
        /// Start-indices of sub-blocks, for the respective block type.
        /// - index: sub-block index within the respective block type
        /// - content: start index of sub-block within block
        /// </returns>
        int[] GetSubblk_i0(int blockType);

        /// <summary>
        /// Block type.
        /// </summary>
        /// <param name="iBlock">
        /// Block index.
        /// </param>
        /// <returns>
        /// The block type, which can e.g. be used as an argument for 
        /// an index into <see cref="GetSubblkLen"/> and <see cref="GetSubblk_i0(int)"/>.
        /// </returns>
        int GetBlockType(long iBlock);

        /// <summary>
        /// Index offset for some block.
        /// </summary>
        /// <param name="iBlock">
        /// Block index.
        /// </param>
        /// <returns>
        /// A global index of the first element in the block.
        /// </returns>
        long GetBlockI0(long iBlock);


        /// <summary>
        /// The Length of some block.
        /// </summary>
        /// <param name="iBlock">
        /// Block index.
        /// </param>
        /// <returns>
        /// The length of the respective block.
        /// </returns>
        int GetBlockLen(long iBlock);

        /// <summary>
        /// Converts an index into a global index.
        /// </summary>
        /// <param name="i">
        /// Global index.
        /// </param>
        /// <returns>
        /// Global block index.
        /// </returns>
        long GetBlockIndex(long i);


        /// <summary>
        /// Total number of blocks (over all MPI processors).
        /// </summary>
        long TotalNoOfBlocks {
            get;
        }

        /// <summary>
        /// Local number of blocks (on this MPI process).
        /// </summary>
        int LocalNoOfBlocks {
            get;
        }

        /// <summary>
        /// Global block index of first block on this MPI process.
        /// </summary>
        long FirstBlock  {
            get;
        }

        /// <summary>
        /// Global block index of first block for some processor
        /// </summary>
        /// <param name="proc">
        /// process index; 
        /// </param>
        /// <returns>index of the first permutation entry stored by processor <paramref name="proc"/></returns>
        long GetFirstBlock(int proc);

        /// <summary>
        /// Returns the number of blocks which are stored by
        /// processor <paramref name="proc"/>;
        /// </summary>
        /// <param name="proc"></param>
        /// <returns></returns>
        int GetLocalNoOfBlocks(int proc);

        /// <summary>
        /// MPI rank for a block.
        /// </summary>
        int FindProcessForBlock(long iBlk);
        
        /// <summary>
        /// Well, true if all blocks are the same; otherwise, maybe something different. 
        /// </summary>
        bool AllBlockSizesEqual {
            get;
        }


        /// <summary>
        /// If this partition is mutable (see <see cref="IPartitioning.IsMutable"/>), this method should 
        /// return a frozen version of this partitioning.
        /// </summary>
        IBlockPartitioning GetImmutableBlockPartitioning(); 

    }



    
}
