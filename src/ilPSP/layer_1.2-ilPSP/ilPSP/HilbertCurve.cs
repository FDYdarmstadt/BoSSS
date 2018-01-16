using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ilPSP.HilbertCurve {
    public static class HilbertCurve {

        // suchen/ersetzen
        // bitmask_t ---> long
        // unsigned ---> int
        // halfmask_t ---> int
        // (shift operations seem to work only with signed types in c#)

        static long bitTranspose(int nDims, int nBits, long inCoords) {
            int nDims1 = nDims - 1;
            int inB = nBits;
            int utB;
            long inFieldEnds = 1;
            long inMask = ((((long)2) << (inB - 1)) - 1);
            long coords = 0;

            while ((utB = inB / 2) != 0) {
                int shiftAmt = nDims1 * utB;
                long utFieldEnds =
           inFieldEnds | (inFieldEnds << (shiftAmt + utB));
                long utMask =
           (utFieldEnds << utB) - utFieldEnds;
                long utCoords = 0;
                int d;
                if ((inB & 1) != 0) {
                    long inFieldStarts = inFieldEnds << (inB - 1);
                    int oddShift = 2 * shiftAmt;
                    for (d = 0; d < nDims; ++d) {
                        long _in = inCoords & inMask;
                        inCoords >>= inB;
                        coords |= (_in & inFieldStarts) << oddShift++;
                        _in &= ~inFieldStarts;
                        _in = (_in | (_in << shiftAmt)) & utMask;
                        utCoords |= _in << (d * utB);
                    }
                } else {
                    for (d = 0; d < nDims; ++d) {
                        long _in = inCoords & inMask;
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

        static void CheckCoords(int nBits, long[] coord) {
            if (nBits < 1)
                throw new ArgumentOutOfRangeException();
            if (nBits >= 64)
                throw new ArgumentOutOfRangeException();
            long CoordMax;
            if (nBits == 63) {
                CoordMax = long.MaxValue; // "(long)1 << (nBits)" would evaluate to -long.MaxValue
            } else {
                CoordMax = (long)1 << (nBits);
            }
            int nDims = coord.Length;
            for (int d = nDims - 1; d >= 0; d--) {
                if (coord[d] < 0)
                    throw new ArgumentOutOfRangeException();
                if (coord[d] >= CoordMax)
                    throw new ArgumentOutOfRangeException();
            }
        }


        public static long hilbert_c2i(int nBits, long[] coord) {
            CheckCoords(nBits, coord);
            int nDims = coord.Length;
            if (nDims > 1) {
                int nDimsBits = nDims * nBits;
                long index;
                int d;
                long coords = 0;
                for (d = nDims - 1; d >= 0; d--) {
                    coords <<= nBits;
                    coords |= coord[d];
                }

                if (nBits > 1) {
                    int ndOnes = ((((int)2) << (nDims - 1)) - 1);
                    int nd1Ones = ndOnes >> 1;
                    int b = nDimsBits;
                    int rotation = 0;
                    int flipBit = 0;
                    long nthbits = ((((long)2) << (nDimsBits - 1)) - 1) / ndOnes;
                    coords = bitTranspose(nDims, nBits, coords);
                    coords ^= coords >> nDims;
                    index = 0;
                    do {
                        int bits = (int)((coords >> (b -= nDims)) & ndOnes);
                        bits = (int)((((flipBit ^ bits) >> (rotation)) | ((flipBit ^ bits) << ((nDims) - (rotation)))) & ((((long)2) << (nDims - 1)) - 1));
                        index <<= nDims;
                        index |= bits;
                        flipBit = (int)1 << rotation;
                        do {
                            bits &= -bits & nd1Ones;
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

        static public void hilbert_i2c(int nBits, long index, long[] coord) {
            int nDims = coord.Length;
            if (nDims > 1) {
                long coords;
                int nbOnes = ((((int)2) << (nBits - 1)) - 1);
                int d;

                if (nBits > 1) {
                    int nDimsBits = nDims * nBits;
                    int ndOnes = ((((int)2) << (nDims - 1)) - 1);
                    int nd1Ones = ndOnes >> 1;
                    int b = nDimsBits;
                    int rotation = 0;
                    int flipBit = 0;
                    long nthbits = ((((long)2) << (nDimsBits - 1)) - 1) / ndOnes;
                    index ^= (index ^ nthbits) >> 1;
                    coords = 0;
                    do {
                        int bits = (int)((index >> (b -= nDims)) & ndOnes);
                        coords <<= nDims;
                        coords |= ((((bits) << (rotation)) | ((bits) >> ((nDims) - (rotation)))) & ((((long)2) << (nDims - 1)) - 1)) ^ flipBit;
                        flipBit = (int)1 << rotation;
                        do {
                            bits &= -bits & nd1Ones;
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