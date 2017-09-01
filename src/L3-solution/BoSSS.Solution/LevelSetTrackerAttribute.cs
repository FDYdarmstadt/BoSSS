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

namespace BoSSS.Solution {
    /// <summary>
    /// Attribute for decorating the <see cref="BoSSS.Foundation.XDG.LevelSetTracker"/>;
    /// </summary>
    [AttributeUsage(AttributeTargets.Field, AllowMultiple = false, Inherited = false)]
    public class LevelSetTrackerAttribute : Attribute {

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="speciesTable">
        /// specified which level set signs map to which species;
        /// e.g. for three species, one may set<c>---:A --+:B -+-:A -++:C +**:D</c>.
        /// </param>
        /// <param name="NearCellWidth"></param>
        public LevelSetTrackerAttribute(string speciesTable, int NearCellWidth) {
            m_SpeciesTable = speciesTable;
            m_NearCellWidth = NearCellWidth;
        }

        string m_SpeciesTable;

        /// <summary>
        /// desired width of near region for the level-set tracker
        /// </summary>
        internal int m_NearCellWidth;


        /// <summary>
        /// computes the species table array from the string that was provided in the constructor;
        /// </summary>
        internal Array GetSpeciesTable(int dim) {

            int[] len = new int[dim];
            for(int d = 0; d < dim; d++)
                len[d] = 2;

            Array speciestable = Array.CreateInstance(typeof(string), len);
            Array touch = Array.CreateInstance(typeof(bool), len);

            try {
                string[] parts = m_SpeciesTable.Split(new char[] { '\n', '\t', ' ', '\r' }, StringSplitOptions.RemoveEmptyEntries);

                int touchcnt = 0;
                foreach(var p in parts) {
                    string[] code_N_name = p.Split(new char[] { ':' }, StringSplitOptions.RemoveEmptyEntries);
                    if(code_N_name.Length != 2)
                        throw new Exception();

                    if(code_N_name[0].Length != dim)
                        throw new Exception();

                    char[] signs = code_N_name[0].ToCharArray();


                    SetRec(dim, signs, code_N_name[1], new int[dim], 0, speciestable, touch, ref touchcnt);
                }


                if(touchcnt != Math.Round(Math.Pow(2, dim)))
                    throw new Exception();
            } catch(Exception) {
                throw new ApplicationException("error in species table");
            }


            return speciestable;

        }

        static void SetRec(int D, char[] signs, string SpecName, int[] idx, int d, Array speciesTable, Array touch, ref int touchcnt) {
            if(d == D) {
                if((bool)(touch.GetValue(idx)))
                    throw new Exception();

                touch.SetValue(true, idx);
                speciesTable.SetValue(SpecName, idx);
                touchcnt++;

            } else {
                switch(signs[d]) {
                    case '-':
                    idx[d] = 0;
                    SetRec(D, signs, SpecName, idx, d + 1, speciesTable, touch, ref touchcnt);
                    return;

                    case '+':
                    idx[d] = 1;
                    SetRec(D, signs, SpecName, idx, d + 1, speciesTable, touch, ref touchcnt);
                    return;

                    case '*':
                    idx[d] = 0;
                    SetRec(D, signs, SpecName, idx, d + 1, speciesTable, touch, ref touchcnt);
                    idx[d] = 1;
                    SetRec(D, signs, SpecName, idx, d + 1, speciesTable, touch, ref touchcnt);
                    return;


                    default: throw new Exception();
                }
            }
        }
    }

}
