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
using System.IO;
using System.Globalization;

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using ilPSP.Utils;

namespace BoSSS.Solution.GridImport.NASTRAN {
    
    
    /// <summary>
    /// 
    /// </summary>
    public class NastranFile {

        /// <summary>
        /// grid points
        /// </summary>
        public List<GRID> m_GRID_List = new List<GRID>();

        /// <summary>
        /// Triangle list
        /// </summary>
        public List<CTRIA3> m_CTRIA3_List = new List<CTRIA3>();


        /// <summary>
        /// renames all ID's in <see cref="m_GRID_List"/> (see <see cref="GRID.ID"/>)
        /// to be equal to the index in <see cref="m_GRID_List"/>
        /// and applies that cahnege also to the Grid id's in <see cref="m_CTRIA3_List"/> (
        /// see <see cref="CTRIA3.Grid1"/>, <see cref="CTRIA3.Grid3"/>, <see cref="CTRIA3.Grid3"/>);
        /// </summary>
        private void RenameIds() {

            SortedDictionary<int, int> id2index = new SortedDictionary<int, int>();

            for (int i = 0; i < m_GRID_List.Count; i++) {
                id2index.Add(m_GRID_List[i].ID, i);
                GRID g = m_GRID_List[i];
                g.ID = i;
                m_GRID_List[i] = g;
            }



            for (int j = 0; j < m_CTRIA3_List.Count; j++) {
                CTRIA3 t = m_CTRIA3_List[j];
                t.Grid1 = id2index[t.Grid1];
                t.Grid2 = id2index[t.Grid2];
                t.Grid3 = id2index[t.Grid3];
                t.ID = j;
                m_CTRIA3_List[j] = t;
            }
        }

        /// <summary>
        /// creates some test grid
        /// </summary>
        /// <param name="N">resolution of testgrid</param>
        public NastranFile(int N) {
            double[] nodes = GenericBlas.Linspace(-3, 3, N);
            double dh = nodes[1] - nodes[0];
            dh *= 0.5;

            int cnt = 0;
            for (int i = 0; i < N; i++) {
                double even = 0;
                if (i % 2 == 0) even = dh;

                for (int j = 0; j < N; j++) {
                    GRID g = new GRID();

                    g.X1 = (float)(nodes[j] + even);
                    g.X2 = (float)nodes[i];
                    g.ID = cnt;

                    m_GRID_List.Add(g);
                    cnt++;
                }
            }


            cnt = 0;
            for (int k = 0; k < (N - 1); k++) {
                int even = 0;
                int odd = 0;
                if (k % 2 == 0) even = 1;
                else odd = 1;

                for (int i = k * N; i < (k * N + N - 1); i++) {

                    CTRIA3 tup = new CTRIA3();
                    tup.Grid1 = i + N + even;
                    tup.Grid2 = i;
                    tup.Grid3 = i + 1;
                    tup.ID = cnt;
                    cnt = cnt + 1;
                    m_CTRIA3_List.Add(tup);


                    CTRIA3 tdw = new CTRIA3();
                    tdw.Grid1 = i + odd;
                    tdw.Grid2 = i + N + 1;
                    tdw.Grid3 = i + N;
                    tup.ID = cnt;
                    cnt = cnt + 1;
                    m_CTRIA3_List.Add(tdw);



                }
            }



        }




        /// <summary>
        /// imports a nastran file 
        /// </summary>
        /// <param name="FilePath"></param>
        public NastranFile(string FilePath) {

            StreamReader rd = null;

            rd = new StreamReader(FilePath);

            string line = rd.ReadLine();
            while (line != null) {

                if (line.StartsWith("GRID")) {
                    GRID g;
                    GRID.Parse(line, out g);
                    m_GRID_List.Add(g);

                }

                if (line.StartsWith("CTRIA3")) {
                    CTRIA3 t;
                    CTRIA3.Parse(line, out t);
                    m_CTRIA3_List.Add(t);
                }



                line = rd.ReadLine();
            }




            rd.Close();


            RenameIds();
        }



        /// <summary>
        /// One line of a NASTRAN file is organized in columns, each one exactly
        /// 8 characters wide. This routine splits the <paramref name="Line"/>
        /// into <paramref name="NoOfColumns"/> substrings which represent the seperate columns.
        /// </summary>
        /// <param name="Line"></param>
        /// <param name="NoOfColumns"></param>
        /// <returns></returns>
        public static string[] SplitNastranLine(string Line, int NoOfColumns) {
            string[] r = new String[NoOfColumns];

            string[] whitespaces = new string[] { " ", "\n", "\r" };

            for (int i = 0; i < NoOfColumns; i++) {
                int ist = i * 8;
                int ien = ist + 8;
                if (ien > Line.Length) ien = Line.Length;
                int l = ien - ist;

                if (l <= 0)
                    r[i] = "";
                else {
                    string s = Line.Substring(ist, ien - ist);

                    string[] pp = s.Split(whitespaces, StringSplitOptions.RemoveEmptyEntries);
                    if (pp.Length > 1)
                        throw new ApplicationException("unknown whitespaces.");
                    if (pp.Length > 0)
                        r[i] = pp[0];
                    else
                        r[i] = "";

                }
            }

            return r;
        }



        /// <summary>
        /// parses a floating point value in a nastran file
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
        public static float ParseFloat(string s) {


            int i = 0;
            if (s[0] == '-') i++;
            if (s[i] == '.') {
                s = s.Insert(i, "0");
                i++;
                i++;
            }

            while (i < s.Length) {
                if (s[i] == '-') {
                    s = s.Insert(i, "E");
                    break;
                }
                i++;
            }
            return float.Parse(s, NumberFormatInfo.InvariantInfo);
        }




        /// <summary>
        /// a NASTRAN GRID - line
        /// </summary>
        public struct GRID {

            /// <summary>
            /// Grid point identification number
            /// </summary>
            public int ID;

            /// <summary>
            /// Identification number of the coordinate system in which the location of the grid points
            /// is defined.
            /// </summary>
            public int CP;

            /// <summary>
            /// First coordinate of the grid point in coordinate system <see cref="CP"/>;
            /// </summary>
            public float X1;

            /// <summary>
            /// Second coordinate of the grid point in coordinate system <see cref="CP"/>;
            /// </summary>
            public float X2;

            /// <summary>
            /// Third coordinate of the grid point in coordinate system <see cref="CP"/>;
            /// </summary>
            public float X3;

            /// <summary>
            /// identification number of coordinate system in which the displacements, degrees of freedom,
            /// constraints and solution vectors are defined at the grid point.
            /// </summary>
            public int CD;

            /// <summary>
            /// Permanent single-point constraits associated with the grid point
            /// </summary>
            public int PS;


            /// <summary>
            /// superelement identification number
            /// </summary>
            public int SEID;


            /// <summary>
            /// 
            /// </summary>
            /// <param name="Line"></param>
            /// <param name="g"></param>
            static public void Parse(string Line, out GRID g) {
                string[] Parts = SplitNastranLine(Line, 9);
                if (!Parts[0].Equals("GRID"))
                    throw new ArgumentException("line must start with GRID", "Line");
                g.ID = int.Parse(Parts[1]);
                g.CP = int.Parse(Parts[2]);
                g.X1 = ParseFloat(Parts[3]);
                g.X2 = ParseFloat(Parts[4]);
                g.X3 = ParseFloat(Parts[5]);
                if (Parts[6].Length > 0) g.CD = int.Parse(Parts[6]); else g.CD = 0;
                if (Parts[7].Length > 0) g.PS = int.Parse(Parts[7]); else g.PS = 0;
                if (Parts[8].Length > 0) g.SEID = int.Parse(Parts[8]); else g.SEID = 0;
            }

        }


        /// <summary>
        /// struct representing the NASTRAN CTRIA3 - line
        /// </summary>
        public struct CTRIA3 {

            /// <summary>
            /// ID of this triangle
            /// </summary>
            public int ID;





            /// <summary>
            /// ID of the frist point of this triangle; (<see cref="GRID.ID"/>);
            /// </summary>
            public int Grid1;


            /// <summary>
            /// ID of the second point of this triangle; (<see cref="GRID.ID"/>);
            /// </summary>
            public int Grid2;


            /// <summary>
            /// ID of the third point of this triangle; (<see cref="GRID.ID"/>);
            /// </summary>
            public int Grid3;


            /// <summary>
            /// 
            /// </summary>
            /// <param name="Line"></param>
            /// <param name="t"></param>
            public static void Parse(string Line, out CTRIA3 t) {
                string[] Parts = SplitNastranLine(Line, 9);
                if (!Parts[0].Equals("CTRIA3"))
                    throw new ArgumentException("line must start with CTRIA3", "Line");
                t.ID = int.Parse(Parts[1]);

                t.Grid1 = int.Parse(Parts[3]);
                t.Grid2 = int.Parse(Parts[4]);
                t.Grid3 = int.Parse(Parts[5]);
            }


        }
    }
}