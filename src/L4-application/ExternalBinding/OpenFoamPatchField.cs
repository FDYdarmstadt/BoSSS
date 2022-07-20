using ilPSP.Connectors;
using ilPSP.Utils;
using ilPSP;
using System;
using System.Collections.Generic;
using BoSSS.Solution.Control;
using BoSSS.Foundation.Grid;


namespace BoSSS.Application.ExternalBinding {

    public class OpenFoamPatchField : IForeignLanguageProxy {

        /// <summary>
        /// Ctor
        /// </summary>
        [CodeGenExport]
        unsafe public OpenFoamPatchField(OpenFOAMGrid grdDat, int nBoundaries, int* edgeTags, int* edgeTypes, double* edgeValues)
            {

                this.Values = new double[nBoundaries];
                this.EdgeTags = new int[nBoundaries];
                this.EdgeTypes = new string[nBoundaries];
                for (int i = 0; i < nBoundaries; i++){
                    this.Values[i] = edgeValues[i];
                    this.EdgeTags[i] = edgeTags[i];
                    this.EdgeTypes[i] = IntToBCType(edgeTypes[i]);
                }
        }

        // [CodeGenExport]
        // unsafe public int BCTypeToInt(char* BCType)
        //     {
        //         string BCTypeS = new string(BCType);
        //         return BCTypes[BCTypeS];
        // }

        public string IntToBCType(int BCTypeInt)
            {
                return BCTypes[BCTypeInt];
        }

        Dictionary<int, string> BCTypes = new Dictionary<int, string>() {
            {0, "dirichlet"},
            {1, "neumann"},
            {-1, "empty"}
        };

        public bool IsDirichlet(int EdgeTag){
            int FaceIndex = -1;
            int i = 0;
            bool found = false;
            foreach (var et in this.EdgeTags){
                if (et == EdgeTag){
                    if (found){
                        throw new ApplicationException("nonunique EdgeTag encountered");
                    }
                    FaceIndex = i;
                    found = true;
                }
                i++;
            }
            if (!found){
                throw new ApplicationException("EdgeTag " + EdgeTag + " not found");
            }
            return this.EdgeTypes[FaceIndex] == "dirichlet";
        }

        public double[] Values;
        public int[] EdgeTags;
        public string[] EdgeTypes;

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
    }
}
