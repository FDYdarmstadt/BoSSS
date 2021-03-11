using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Foundation.Quadrature {
    
    /// <summary>
    /// (Mainly internal) Utility routines for the integration on periodic boundaries
    /// </summary>
    public static class PeriodicBoundaryUtils {


        /// <summary>
        /// <see cref="ISpatialOperator.VectorFieldIndices"/>
        /// </summary>
        static public IEnumerable<int[]> GetVectorFieldIndices(IEnumerable<string> DomVar) {

            var ret = new List<int[]>();

            var difference = new char[][] {
                    new char[] { 'X', 'Y', 'Z' },
                    new char[] { 'x', 'y', 'z' },
                    new char[] { '0', '1', '2' }
                };

            bool IsComponentChar( int d, char c) {
                foreach(char[] V in difference) {
                    if(V[d] == c)
                        return true;
                }

                return false;
            }


            string[] _DomVar = DomVar.ToArray();

            for(int iVar = 0; iVar < _DomVar.Length; iVar++) {
                var DVi = _DomVar[iVar];
                if(DVi == null)
                    continue;
                var curList = new List<int>();
                curList.Add(iVar);






                for(int jVar = iVar + 1; jVar < _DomVar.Length; jVar++) {
                    var DVj = _DomVar[jVar];
                    if(DVj == null)
                        continue;

                    if(DVi.Length == DVj.Length) {
                        int NoOfCharsDiff = 0;
                        int DiffPosition = -1;
                        for(int k = 0; k < DVi.Length; k++) {
                            if(DVi[k] != DVj[k]) {
                                NoOfCharsDiff++;
                                DiffPosition = k;
                            }
                        }

                        if(NoOfCharsDiff != 1)
                            continue;

                        bool ok = true;
                        for(int d = 0; d < curList.Count; d++) {
                            if(!IsComponentChar(d, _DomVar[curList[d]][DiffPosition])) {
                                ok = false;
                            }
                        }
                        if(!IsComponentChar(curList.Count, DVj[DiffPosition])) {
                            ok = false;
                        }

                        if(!ok)
                            continue;

                        curList.Add(jVar);
                    }
                }

                if(curList.Count > 1) {
                    ret.Add(curList.ToArray());
                    string[] fucker = curList.Select(i => _DomVar[i]).ToArray();
                    foreach(int jVar in curList) {
                        _DomVar[jVar] = null;
                    }
                }
            }

            return ret.ToArray();

        }

    }
}
