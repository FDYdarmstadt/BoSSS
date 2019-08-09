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
using System.Diagnostics;
using System.Runtime.InteropServices;
using System.Runtime.CompilerServices;

namespace ilPSP {

    
    /// <summary>
    /// Manages two buffers for temporary use;
    /// use with care -- never return the buffer objects!
    /// </summary>
    public static class TempBuffer {

        const int MAX_BUFFERS = 8;

        static WeakReference[] BigBuf = new WeakReference[MAX_BUFFERS];
        static bool[] BufLok = new bool[MAX_BUFFERS];

        /// <summary>
        /// Only use this if you know what you are doing: you're about to theave tha save .NET-world.
        /// Some <c>stackalloc</c> would be ideal for this, but it may create stack overflow...
        /// </summary>
        public static double[] GetTempBuffer(out int _iBuf, int size) {
            lock(BufLok) {
                double[] BB = null;
                int iAlloc = -1;
                int iBuf = -1;

                // try to find a buffer which already fits...
                for(int iB = 0; iB < MAX_BUFFERS; iB++) {
                    if(BufLok[iB] == false) {
                        if(BigBuf[iB] != null && BigBuf[iB].IsAlive) {
                            double[] C = (double[])BigBuf[iB].Target;
                            if(C != null && C.Length >= size) {
                                if(BB == null) {
                                    // passenden Buffer gefunden
                                    BB = C;
                                    iBuf = iB;
                                } else {
                                    if(C.Length < BB.Length) {
                                        // noch besser passenden Buffer gefunden
                                        BB = C;
                                        iBuf = iB;
                                    }
                                }
                            } else {
                                // Buffer ist tot, oder zu klein
                                // => könnte (re-)allokiert werden, falls nötig.
                                iAlloc = iB;
                            }
                        } else {
                            // Buffer ist tot
                            // => könnte reallokiert werden, falls nötig.
                            iAlloc = iB;
                        }
                    }
                }

                if(BB != null) {
                    _iBuf = iBuf;
                    BufLok[iBuf] = true;
					return BB;
                }

                // (re-) allocation necessary
                if(iAlloc < 0)
                    throw new ApplicationException("running out of temp buffers.");

                _iBuf = iAlloc;
                BB = new double[size];
                BigBuf[iAlloc] = new WeakReference(BB);
                BufLok[iAlloc] = true;
				return BB;
            }
        }

        /// <summary>
        /// Removes the lock on a buffer (see <see cref="BufLok"/>)
        /// </summary>
        /// <param name="ibuf"></param>
        public static void FreeTempBuffer(int ibuf) {
            lock(BufLok) {
                // "Sperre" aufheben
                Debug.Assert(BufLok[ibuf] == true);
                BufLok[ibuf] = false;
				Debug.Assert(BufLok[ibuf] == false);
            }
        }

        /// <summary>
        /// Only use this if you know what you are doing: you're about to theave tha save .NET-world.
        /// </summary>
        public static MultidimensionalArray GetTempMultidimensionalarray(out int iBuf, params int[] Lengths) {
            int size = 1;
            for(int i = Lengths.Length - 1; i >= 0; i--)
                size *= Lengths[i];

            double[] buffer = GetTempBuffer(out iBuf, size);

            var R = MultidimensionalArray.CreateWrapper(buffer, Lengths);
            R.Clear();
            return R;
        }

    }
}
