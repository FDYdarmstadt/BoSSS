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
using BoSSS.Foundation.Grid;
using System.IO;
using BoSSS.Platform;
using ilPSP.Utils;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.GridImport {
    
     public class Cgns : IGridImporter {

        public CGNSBase_t[] base_t;

        /// <summary>
        /// reads the CGNS tree from a file
        /// </summary>
        public Cgns(string filePath) {
            if(!File.Exists(filePath))
                throw new FileNotFoundException("Given grid file does not exist -- no import started.", filePath);

            bool hdf5 = filePath.ToLowerInvariant().EndsWith(".hdf") || filePath.ToLowerInvariant().EndsWith(".hdf5");

            using (CgnsDriver cg = new CgnsDriver(hdf5)) {
                
                // open file
                int index_file = 1;
                if (cg.open(filePath.ToNullTermCharAry(), (int)Mode.MODE_READ, out index_file) != (int)Error.CG_OK) {
                    Console.WriteLine("Critical error opening cgns file -- application may be unrecoverable.");
                    Console.Out.Flush();
                    ThrowError(index_file, cg);
                }

                //try {

                    // load all base-nodes
                    int nbases;
                    if (cg.nbases(index_file, out nbases) != (int)Error.CG_OK) {
                        ThrowError(index_file, cg);
                    }

                    this.base_t = new CGNSBase_t[nbases];
                    for (int i = 1; i <= nbases; i++) {
                        this.base_t[i - 1] = new CGNSBase_t(cg, index_file, i);
                    }

                //} catch (Exception e) {
                //    throw e;
                //} finally {
                    // close file
                    if (cg.close(index_file) != (int)Error.CG_OK) {
                        ThrowError(cg);
                    }
                //}
            }
        }

        


        public class CGNSBase_t {
            public int CellDimension;

            public int PhysicalDimension;

            public string BaseName;

            public Zone_t[] zones;
            
            internal CGNSBase_t(CgnsDriver cg, int index_file, int index_base) {

                char[] basename = new char[1000];
                if (cg.base_read(index_file, index_base, basename, out CellDimension, out PhysicalDimension) != (int)Error.CG_OK) {
                    ThrowError(index_file, cg);
                }
                this.BaseName = basename.FromNullTermCharAry();

                int nzones;
                if (cg.nzones(index_file, index_base, out nzones) != (int)Error.CG_OK) {
                    ThrowError(index_file, cg);
                }

                zones = new Zone_t[nzones];
                for (int index_zone = 1; index_zone <= nzones; index_zone++) {
                    ZoneType_t zone_type;
                    if (cg.zone_type(index_file, index_base, index_zone, out zone_type) != (int)Error.CG_OK) {
                        ThrowError(index_file, cg);
                    }
                    
                    switch (zone_type) {
                        case ZoneType_t.Unstructured: zones[index_zone - 1] = new UnstructuredZone_t(this, cg, index_file, index_base, index_zone); break;

                        default:
                        throw new NotImplementedException("Zone type " + zone_type + " not implemented for CGNS.");
                    }
                }

                
            }
        }

        /// <summary>
        /// baseclass for both, structured and unstructured zones
        /// </summary>
        public abstract class Zone_t {

            public ZoneType_t zone_type;

            public string ZoneName;
        }


        public class UnstructuredZone_t : Zone_t {
                       

            /// <summary>
            /// http://www.grc.nasa.gov/WWW/cgns/CGNS_docs_rel25/midlevel/structural.html
            /// </summary>
            int[,] isize = new int[3, 3];


            public int ncoords;

            public CGNSBase_t m_owner;

            public string[] coordname;
            
            /// <summary>
            /// data-type in CGNS file; in the C#-part, we only store double precision, and convert if necessary.
            /// </summary>
            public DataType_t[] data_type;

            /// <summary>
            /// 1st index: spatial direction
            /// 2nd index: node
            /// </summary>
            public double[][] coord;

            internal UnstructuredZone_t(CGNSBase_t owner, CgnsDriver cg, int index_file, int index_base, int index_zone) {
                m_owner = owner;

                // general info about the zone
                // ===========================
                {
                    char[] zonename = new char[1000];
                    if (cg.zone_read(index_file, index_base, index_zone, zonename, isize) != (int)Error.CG_OK) {
                        ThrowError(index_file, cg);
                    }
                    ZoneName = zonename.FromNullTermCharAry();

                    if (cg.zone_type(index_file, index_base, index_zone, out zone_type) != (int)Error.CG_OK) {
                        ThrowError(index_file, cg);
                    }
                    if (this.zone_type != ZoneType_t.Unstructured) {
                        throw new NotSupportedException("supporting only unstructured grids");
                    }

                    if (cg.ncoords(index_file, index_base, index_zone, out ncoords) != (int)Error.CG_OK) {
                        ThrowError(index_file, cg);
                    }
                }

                // x,y, and z - Nodes
                // ==================
                {
                    int D = m_owner.CellDimension;

                    coordname = new string[3];
                    data_type = new DataType_t[3];

                    for (int d = 0; d < 3; d++) {
                        char[] _coordname = new char[1000];
                        if (cg.coord_info(index_file, index_base, index_zone, d + 1, out data_type[d], _coordname) != (int)Error.CG_OK) {
                            ThrowError(index_file, cg);
                        }
                        coordname[d] = _coordname.FromNullTermCharAry();
                    }

                    coord = new double[3][];
                    for (int d = 0; d < 3; d++) {
                        if (data_type[d] != DataType_t.RealDouble)
                            throw new NotSupportedException();

                        coord[d] = new double[isize[0, 0]];

                        int[] irmin = new int[] { 1 };
                        int[] irmax = new int[] { coord[d].Length };

                        if (cg.coord_read1(index_file, index_base, index_zone, coordname[d].ToNullTermCharAry(), data_type[d], irmin, irmax, coord[d]) != (int)Error.CG_OK) {
                            ThrowError(index_file, cg);
                        }
                    }
                    if (D < 3) {
                        double[] coordNorm = coord.Select(coord_d => coord_d.L2Norm()).ToArray();

                        int SkipDimFound = 0;
                        for (int d = 0; d < 3; d++) {
                            if (coordNorm[d] == 0) {
                                SkipDimFound++;

                                ArrayTools.RemoveAt(ref coordname, d);
                                ArrayTools.RemoveAt(ref data_type, d);
                                ArrayTools.RemoveAt(ref coord, d);
                            }
                        }

                        if ((SkipDimFound + D) != 3) {
                            throw new ArgumentException("CGNS node index out of range.");
                        }

                    }
                }

                // sections
                // ========

                {
                    int nsections;
                    if (cg.nsections(index_file, index_base, index_zone, out nsections) != (int)Error.CG_OK) {
                        ThrowError(index_file, cg);
                    }

                    elements = new Elements_t[nsections];
                    for (int index_section = 1; index_section <= nsections; index_section++) {
                        elements[index_section - 1] = new Elements_t(this, cg, index_file, index_base, index_zone, index_section);
                    }
                }

                // boundary conditions
                // ===================
                {
                    int nbocos;
                    if (cg.nbocos(index_file, index_base, index_zone, out nbocos) != (int)Error.CG_OK) {
                        ThrowError(index_file, cg);
                    }

                    bcs = new BC_t[nbocos];
                    for (int index_bc = 1; index_bc <= nbocos; index_bc++) {
                        bcs[index_bc - 1] = new BC_t(this, cg, index_file, index_base, index_zone, index_bc);
                    }

                }

            }

            public Elements_t[] elements;

            public BC_t[] bcs;
        }


        public class BC_t {
            Zone_t m_Owner;

            public string Name;

            public int[] normal_index;

            public BCType_t bc_type;

            public PointSetType_t pointset_type;

            public int[] bcelement;

            internal BC_t(Zone_t Owner, CgnsDriver cg, int index_file, int index_base, int index_zone, int index_bc) {
                m_Owner = Owner;

                char[] boconame = new char[1000];
                normal_index = new int[3];
                int normal_list_flag, nbcelement, ndataset; // bc definition
                DataType_t data_type;
                if (cg.boco_info(index_file, index_base, index_zone, index_bc, boconame, out bc_type, out pointset_type, out nbcelement, normal_index, out normal_list_flag, out data_type, out ndataset) != (int)Error.CG_OK) {
                    ThrowError(index_file, cg);
                }
                Name = boconame.FromNullTermCharAry();

                bcelement = new int[nbcelement];
                if (cg.boco_read1(index_file, index_base, index_zone, index_bc, bcelement, null) != (int)Error.CG_OK) {
                    ThrowError(index_file, cg);
                }

                
            }

        }


        public class Elements_t {
            Zone_t m_Owner;

            public ElementType_t element_type;

            public string Name;

            public int start;

            public int end;
            
            public int nbndry;
            
            public int parent_flag;

            public int[,] ielem;

            internal Elements_t(Zone_t Owner, CgnsDriver cg, int index_file, int index_base, int index_zone, int index_section) {
                m_Owner = Owner;

                char[] secname = new char[1000];
                if (cg.section_read(index_file, index_base, index_zone, index_section, secname, out element_type, out start, out end, out nbndry, out parent_flag) != (int)Error.CG_OK) {
                    ThrowError(index_file, cg);
                }
                this.Name = secname.FromNullTermCharAry();

                ielem = new int[end - start + 1, element_type.NoOfNodes()];

                if (cg.elements_read2(index_file, index_base, index_zone, index_section, ielem, null) != (int)Error.CG_OK) {
                    ThrowError(index_file, cg);
                }
            }
        }




        static private void ThrowError(CgnsDriver cg) {
            throw new ApplicationException(cg.get_error());
        }

        static private void ThrowError(int index_file, CgnsDriver cg) {
            String s = cg.get_error();
            throw new ApplicationException(s);
        }

        static private void ThrowError(String s, int index_file, CgnsDriver cg) {
            throw new ApplicationException(s);
        }

        #region IImporter Members

        public GridCommons GenerateBoSSSGrid() {
            return CGNS2BoSSSGrid.ToBoSSSGrid(this);
        }

        #endregion
     }
}
