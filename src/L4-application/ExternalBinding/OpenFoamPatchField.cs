using ilPSP.Connectors;
using ilPSP.Utils;
using ilPSP;
using System;
using System.Collections.Generic;
using BoSSS.Solution.Control;
using BoSSS.Foundation.Grid;
using System.Linq;


namespace BoSSS.Application.ExternalBinding {

    public class OpenFoamPatchField : IForeignLanguageProxy {

        // /// <summary>
        // /// Ctor
        // /// </summary>
        // [CodeGenExport]
        // unsafe public OpenFoamPatchField(OpenFOAMGrid grdDat, int nBoundaries, int* edgeTags, int* edgeTypes, double* edgeValues)
        //     {

        //         this.Values = new List<List<double>>();
        //         this.EdgeTags = new int[nBoundaries];
        //         this.EdgeTypes = new string[nBoundaries];
        //         for (int i = 0; i < nBoundaries; i++){
        //             this.Values.Add(new List<double>{edgeValues[i]});
        //             this.EdgeTags[i] = edgeTags[i];
        //             this.EdgeTypes[i] = IntToBCType(edgeTypes[i]);
        //         }
        // }
        public OpenFOAMGrid Grid;

        /// <summary>
        /// Ctor
        /// </summary>
        [CodeGenExport]
        unsafe public OpenFoamPatchField(OpenFOAMGrid grdDat, int dim, int nBoundaries, int* edgeTags, int* edgeTypes, double* edgeValues)
            {
                this.Grid = grdDat;
                this.Values = new List<List<double>>();
                this.EdgeTags = new int[nBoundaries];
                this.EdgeTypes = new string[nBoundaries];
                for (int i = 0; i < nBoundaries; i++){
                    this.Values.Add(new List<double>());
                    for (int d = 0; d < dim; d++) {
                    this.Values[i].Add(edgeValues[dim * i + d]);
                }
                    this.EdgeTags[i] = edgeTags[i];
                    this.EdgeTypes[i] = IntToBCType(edgeTypes[i]);
                }
        }

        public OpenFoamPatchField(OpenFOAMGrid grdDat, int dim, int[] edgeTags, string[] edgeTypes, double[] edgeValues)
            {
                this.Grid = grdDat;
                this.EdgeTags = edgeTags;
                this.Values = new List<List<double>>();
                for (int i = 0; i < edgeTags.Length; i++){
                    this.Values.Add(new List<double>());
                    for (int d = 0; d < dim; d++)
                    {
                        try { this.Values[i].Add(edgeValues[dim * i + d]);
                        } catch (Exception e) {
                            Console.WriteLine("Index i : " + i);
                            Console.WriteLine("Index d : " + d);
                            Console.WriteLine("dim : " + dim);
                            Console.WriteLine("Index dim * i + d : " + (int)(dim * i + d));
                            Console.WriteLine("Length of edgeValues: " + edgeValues.Length);
                            throw e;
                        }
                    }
                }
                this.EdgeTypes = edgeTypes;
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
            // foreach (var et in this.EdgeTags){
            //     Console.WriteLine(et);
            // }
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

        public List<List<double>> Values;
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
