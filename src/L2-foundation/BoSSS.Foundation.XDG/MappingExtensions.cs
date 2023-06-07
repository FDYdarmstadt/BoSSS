using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.XDG {

    /// <summary>
    /// XDG-related extensions for <see cref="UnsetteledCoordinateMapping"/> and <see cref="CoordinateMapping"/>
    /// </summary>
    static public class MappingExtensions {


        /// <summary>
        /// computes a local unique coordinate index ("local" means local on this processor);
        /// this index is unique over all fields (in this mapping), over all cells, over all basis functions, 
        /// but it's only locally (on this processor) valid.
        /// </summary>
        /// <param name="find">
        /// the field or basis index (see <see cref="UnsetteledCoordinateMapping.BasisS"/>);
        /// </param>
        /// <param name="j">local cell index</param>
        /// <param name="n">DG mode index</param>
        /// <param name="lsTrk"></param>
        /// <param name="map"></param>
        /// <param name="spc">species</param>
        /// <returns>
        /// A (MPI-) local index in the update range 
        /// (between 0 and smaller than <see cref="ilPSP.IPartitioning.LocalLength"/>).
        /// It can be converted into 
        /// a (MPI-) global index by adding <see cref="ilPSP.IPartitioning.i0"/>.
        /// </returns>
        public static int LocalUniqueCoordinateIndex(this UnsetteledCoordinateMapping map, LevelSetTracker lsTrk, int find, int j, SpeciesId spc, int n) {
            int spcIdx = lsTrk.Regions.GetSpeciesIndex(spc, j);
            return LocalUniqueCoordinateIndex(map, lsTrk, find, j, spcIdx, n);
        }


        /// <summary>
        /// Number of DG modes for a specific variable, cell and species.
        /// </summary>
        /// <param name="find">
        /// the field or basis index (see <see cref="UnsetteledCoordinateMapping.BasisS"/>);
        /// </param>
        /// <param name="j">local cell index</param>
        /// <param name="lsTrk"></param>
        /// <param name="map"></param>
        /// <param name="spc">species</param>
        public static int GetNumberOfModes(this UnsetteledCoordinateMapping map, LevelSetTracker lsTrk, int find, int j, SpeciesId spc) {
            int spcIdx = lsTrk.Regions.GetSpeciesIndex(spc, j);
            return GetNumberOfModes(map, lsTrk, find, j, spcIdx);
        }

        /// <summary>
        /// Number of DG modes for a specific variable, cell and species.
        /// </summary>
        /// <param name="find">
        /// the field or basis index (see <see cref="UnsetteledCoordinateMapping.BasisS"/>);
        /// </param>
        /// <param name="j">local cell index</param>
        /// <param name="lsTrk"></param>
        /// <param name="map"></param>
        /// <param name="spcIdx">index</param>
        public static int GetNumberOfModes(this UnsetteledCoordinateMapping map, LevelSetTracker lsTrk, int find, int j, int spcIdx) {
            int NoOfSpc = lsTrk.Regions.GetNoOfSpecies(j, out var rrc);
            if(spcIdx < 0 || spcIdx >= NoOfSpc)
                throw new IndexOutOfRangeException($"Species index out of range (species index is {spcIdx}, Number of species is {NoOfSpc}, cell {j})");

            Basis b = map.BasisS[find];
            XDGBasis xb = b as XDGBasis;

            if(xb == null)
                return b.Length;
            else
                return xb.DOFperSpeciesPerCell;
        }


        /// <summary>
        /// computes a local unique coordinate index ("local" means local on this processor);
        /// this index is unique over all fields (in this mapping), over all cells, over all basis functions, 
        /// but it's only locally (on this processor) valid.
        /// </summary>
        /// <param name="find">
        /// the field or basis index (see <see cref="UnsetteledCoordinateMapping.BasisS"/>);
        /// </param>
        /// <param name="j">local cell index</param>
        /// <param name="n">DG mode index</param>
        /// <param name="lsTrk"></param>
        /// <param name="map"></param>
        /// <param name="spcIdx">species index</param>
        /// <returns>
        /// A (MPI-) local index in the update range 
        /// (between 0 and smaller than <see cref="ilPSP.IPartitioning.LocalLength"/>).
        /// It can be converted into 
        /// a (MPI-) global index by adding <see cref="ilPSP.IPartitioning.i0"/>.
        /// </returns>
        public static int LocalUniqueCoordinateIndex(this UnsetteledCoordinateMapping map, LevelSetTracker lsTrk, int find, int j, int spcIdx, int n) {
            //;
            
            int NoOfSpc = lsTrk.Regions.GetNoOfSpecies(j, out var rrc);
            if(spcIdx < 0 || spcIdx >= NoOfSpc)
                throw new IndexOutOfRangeException($"Species index out of range (species index is {spcIdx}, Number of species is {NoOfSpc}, cell {j})");

            Basis b = map.BasisS[find];
            XDGBasis xb = b as XDGBasis;

            if(xb == null) {
                return map.LocalUniqueCoordinateIndex(find, j, n);


            } else {
                int i0 = map.LocalUniqueCoordinateIndex(find, j, 0);
                int Ns = xb.NonX_Basis.GetLength(j);
                if(n < 0 || n >= Ns) {
                    throw new IndexOutOfRangeException($"DG mode index (within species) of range (DG mode index is {n}, but non-XDG basis length is {Ns}, cell {j})");
                }

                return i0 + Ns * spcIdx + n;
            }
        }


        /// <summary>
        /// Returns the length (number of DG modes for all species) 
        /// </summary>
        public static int[] GetBasisLengths(this UnsetteledCoordinateMapping map, int j) {
            var basisS = map.BasisS;
            int[] N = new int[basisS.Count];
            for(int f = 0; f< N.Length; f++) {
                N[f] = basisS[f].GetLength(j);
            }
            return N;
        }
        
        /// <summary>
        /// Returns the length (number of DG modes for all species) 
        /// </summary>
        public static int[] GetNonXBasisLengths(this UnsetteledCoordinateMapping map, int j) {
            var basisS = map.BasisS;
            int[] N = new int[basisS.Count];
            for(int f = 0; f< N.Length; f++) {
                var b = basisS[f];
                if(b is XDGBasis xb) {
                    var _b = xb.NonX_Basis;
                    N[f] = _b.GetLength(j);
                } else {
                    N[f] = b.GetLength(j);
                }
            }
            return N;
        }

        /// <summary>
        /// Returns a bool-array to indicate whether the individual basis components <see cref="UnsetteledCoordinateMapping.BasisS"/> are XDG or normal DG
        /// </summary>
        public static bool[] XorNonXbasis(this UnsetteledCoordinateMapping map) {
            var basisS = map.BasisS;
            bool[] N = new bool[basisS.Count];
            for(int f = 0; f< N.Length; f++) {
                var b = basisS[f];
                N[f] = b is XDGBasis xb;
            }
            return N;
        }
    }
}
