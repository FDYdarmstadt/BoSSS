using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ilPSP.HilbertCurve {

    public static class HilbertCurve {




        // suchen/ersetzen
        // bitmask_t ---> ulong
        // unsigned ---> int
        // halfmask_t ---> uint
        // int must be used for shift operation

        static ulong bitTranspose(int nDims, int nBits, ulong inCoords) {
            int nDims1 = nDims - 1;
            int inB = nBits;
            int utB;
            ulong inFieldEnds = 1;
            ulong inMask = ((((ulong)2) << (inB - 1)) - 1);
            ulong coords = 0;

            while ((utB = inB / 2) != 0) {
                int shiftAmt = nDims1 * utB;
                ulong utFieldEnds =
           inFieldEnds | (inFieldEnds << (shiftAmt + utB));
                ulong utMask =
           (utFieldEnds << utB) - utFieldEnds;
                ulong utCoords = 0;
                int d;
                if ((inB & 1) != 0) {
                    ulong inFieldStarts = inFieldEnds << (inB - 1);
                    int oddShift = 2 * shiftAmt;
                    for (d = 0; d < nDims; ++d) {
                        ulong _in = inCoords & inMask;
                        inCoords >>= inB;
                        coords |= (_in & inFieldStarts) << oddShift++;
                        _in &= ~inFieldStarts;
                        _in = (_in | (_in << shiftAmt)) & utMask;
                        utCoords |= _in << (d * utB);
                    }
                } else {
                    for (d = 0; d < nDims; ++d) {
                        ulong _in = inCoords & inMask;
                        inCoords >>= inB;
                        _in = (_in | (_in << shiftAmt)) & utMask;
                        utCoords |= _in << (d * utB);
                    }
                }
                inCoords = utCoords;
                inB = utB;
                inFieldEnds = utFieldEnds;
                inMask = utMask;
            }
            coords |= inCoords;
            return coords;
        }

        static void CheckCoords(int nBits, ulong[] coord) {
            int nDims=coord.Length;
            if (nDims < 1)
                throw new ArgumentNullException();
            if (nBits < 1)
                throw new ArgumentOutOfRangeException();
            if (nBits*nDims > 64)   //64=8*sizeof(ulong)
                throw new ArgumentOutOfRangeException();
            ulong CoordMax = (ulong)1 << (nBits);
            for (int d = nDims - 1; d >= 0; d--) {
                if (coord[d] < 0)
                    throw new ArgumentOutOfRangeException();
                if (coord[d] >= CoordMax)
                    throw new ArgumentOutOfRangeException();
            }
        }


        public static ulong hilbert_c2i(int nBits, ulong[] coord) {
            CheckCoords(nBits, coord);
            int nDims = coord.Length;
            if (nDims > 1) {
                int nDimsBits = nDims * nBits;
                ulong index;
                int d;
                ulong coords = 0;
                for (d = nDims - 1; d >= 0; d--) {
                    coords <<= nBits;
                    coords |= coord[d];
                }

                if (nBits > 1) {
                    uint ndOnes = ((((uint)2) << (nDims - 1)) - 1);
                    uint nd1Ones = ndOnes >> 1;
                    int b = nDimsBits;
                    int rotation = 0;
                    uint flipBit = 0;
                    ulong nthbits = ((((ulong)2) << (nDimsBits - 1)) - 1) / ndOnes;
                    coords = bitTranspose(nDims, nBits, coords);
                    coords ^= coords >> nDims;
                    index = 0;
                    do {
                        uint bits = (uint)((coords >> (b -= nDims)) & ndOnes);
                        bits = (uint)((((flipBit ^ bits) >> (rotation)) | ((flipBit ^ bits) << ((nDims) - (rotation)))) & ((((ulong)2) << (nDims - 1)) - 1));
                        index <<= nDims;
                        index |= bits;
                        flipBit = (uint)1 << rotation;
                        do {
                            bits &= (uint)(-bits & nd1Ones);
                            while (bits != 0) {
                                bits >>= 1;
                                ++rotation;
                            }
                            if (++rotation >= nDims)
                                rotation -= nDims;
                        } while (false);
                    } while (b != 0);
                    index ^= nthbits >> 1;
                } else
                    index = coords;
                for (d = 1; d < nDimsBits; d *= 2)
                    index ^= index >> d;
                return index;
            } else
                return coord[0];
        }

        static public void hilbert_i2c(int nBits, ulong index, ulong[] coord) {
            int nDims = coord.Length;
            if (nDims > 1) {
                ulong coords;
                uint nbOnes = ((((uint)2) << (nBits - 1)) - 1);
                int d;

                if (nBits > 1) {
                    int nDimsBits = nDims * nBits;
                    uint ndOnes = ((((uint)2) << (nDims - 1)) - 1);
                    uint nd1Ones = ndOnes >> 1;
                    int b = nDimsBits;
                    int rotation = 0;
                    uint flipBit = 0;
                    ulong nthbits = ((((ulong)2) << (nDimsBits - 1)) - 1) / ndOnes;
                    index ^= (index ^ nthbits) >> 1;
                    coords = 0;
                    do {
                        uint bits = (uint)((index >> (b -= nDims)) & ndOnes);
                        coords <<= nDims;
                        coords |= ((((bits) << (rotation)) | ((bits) >> ((nDims) - (rotation)))) & ((((ulong)2) << (nDims - 1)) - 1)) ^ flipBit;
                        flipBit = (uint)1 << rotation;
                        do {
                            bits &= (uint)(-bits & nd1Ones);
                            while (bits != 0) {
                                bits >>= 1;
                                ++rotation;
                            }
                            if (++rotation >= nDims)
                                rotation -= nDims;
                        } while (false);
                    } while (b != 0);
                    for (b = nDims; b < nDimsBits; b *= 2)
                        coords ^= coords >> b;
                    coords = bitTranspose(nBits, nDims, coords);
                } else
                    coords = index ^ (index >> 1);

                for (d = 0; d < nDims; ++d) {
                    coord[d] = coords & nbOnes;
                    coords >>= nBits;
                }
            } else {
                coord[0] = index;
            }
#if DEBUG
            CheckCoords(nBits, coord);
#endif
        }

    }

}
