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
using System.Runtime.Serialization;

namespace BoSSS.Foundation.IO {

    /// <summary>
    /// Items of the data-vector that is stored for one time-step; an instance
    /// of this class contains basically the DG coordinates of multiple DG
    /// fields, for one cell. 
    /// </summary>
    [Serializable]
    [DataContract]
    public class CellFieldDataSet {

        /// <summary>
        /// Global ID of the cell
        /// </summary>
        [DataMember]
        public long GlobalID;

        /// <summary>
        ///  a buffer of data -- must be altered with unsafe methods -- Har, har!
        ///  - 0th int: No of stored DG fields, Q;
        ///  - 1st int to Q-th int: number of DG coordinates for the respective field
        ///  after that: doubles that store the DG data
        /// </summary>
        [DataMember]
        private byte[] DGData;

        /// <summary>
        /// Appends DG coordinate data to <see cref="DGData"/>
        /// </summary>
        /// <param name="Coords">DG coordinates in respective cell</param>
        public void AppendDGCoordinates(double[] Coords) {
            if (DGData == null)
                DGData = new byte[sizeof(int)];
            int Nofldold = this.NoOfFields;
            int NofldNew = Nofldold + 1;

            int N = Coords.Length;
            byte[] NewDGData = new byte[this.DGData.Length + sizeof(int) + sizeof(double) * N];

            unsafe {
                fixed (byte* pOldDGData = DGData, pNewDGData = NewDGData) {
                    *((int*)pNewDGData) = NofldNew;

                    Array.Copy(this.DGData, sizeof(int), NewDGData, sizeof(int), sizeof(int) * Nofldold);
                    ((int*)pNewDGData)[NofldNew] = N;

                    int oOld = sizeof(int) * (Nofldold + 1);
                    int oNew = sizeof(int) * (NofldNew + 1);
                    int Len = DGData.Length - oOld;
                    Array.Copy(this.DGData, oOld, NewDGData, oNew, Len);

                    fixed (double* pCoordsOrg = Coords) {
                        double* pCoordsTrg = (double*)(pNewDGData + sizeof(int) * (NofldNew + 1) + Len);

                        for (int n = 0; n < N; n++) {
                            pCoordsTrg[n] = pCoordsOrg[n];
                        }
                    }
                }
            }

            this.DGData = NewDGData;
        }

        /// <summary>
        /// number of fields stored in this data set.
        /// </summary>
        public int NoOfFields {
            get {
                if (DGData == null || DGData.Length <= 0) {
                    return 0;
                }

                unsafe {
                    fixed (void* pDGData = DGData) {
                        int ret;
                        ret = *((int*)pDGData);
                        return ret;
                    }
                }
            }
        }

        /// <summary>
        /// returns the DG coordinates of the <paramref name="ifld"/>-th field.
        /// </summary>
        public double[] GetDGCoordinates(int ifld) {
            int Nofld = this.NoOfFields;
            if (ifld >= Nofld || ifld < 0) {
                throw new IndexOutOfRangeException();
            }

            unsafe {
                fixed (byte* pDGData = DGData) {
                    int offset = 0;
                    for (int i = 0; i < ifld; i++)
                        offset += ((int*)pDGData)[i + 1];
                    int N = ((int*)pDGData)[ifld + 1];
                    double[] Coords = new double[N];
                    fixed (double* pCoordsTarg = Coords) {
                        byte* pCoordSrc = ((byte*)pDGData + sizeof(int) * (Nofld + 1) + sizeof(double) * offset);
                        byte* pCoordTrg = (byte*)pCoordsTarg;

                        for (int i = 0; i < sizeof(double) * N; i++) {
                            pCoordTrg[i] = pCoordSrc[i];
                        }
                    }
                    return Coords;
                }
            }
        }
    }
}
