using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using BoSSS.Foundation.Grid;

namespace BoSSS.Solution.GridImport {
#pragma warning disable 1591
    public class CgnsFile {

        public int celldim, physdim;

        /// <summary>
        /// loads or saves an adf or hdf5 file from or to a specified path
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="hdf5"></param>
        public CgnsFile(string filePath, bool hdf5) {
            using (CgnsDriver cg = new CgnsDriver(hdf5)) {
                int index_file = 0;
                int index_base = 0;
                int index_zone = 0;
                //int index_section = 0;
                //int index_bc = 0;
                int nbases;
                int nzones;
                //int ncoords;
                //int nsections;
                //int nbocos;
                char[] basename = new char[100];
                //char[] zonename = new char[100];
                //char[] coordnameX = new char[100];
                //char[] coordnameY = new char[100];
                //char[] coordnameZ = new char[100];
                //char[] flow = new char[100];
                //char[] secname = new char[100];
                //char[] boconame;
                //char[] solname = new char[100];
                //int[,] isize = new int[3, 3];
                //int[,] ielem;
                //int[] irmin = new int[3];
                //int[] irmax = new int[3];
                // zone definition
                //int start, end, nbndry, parent_flag; // element definition
                //int normal_list_flag, nbcelement, ndataset; // bc definition
                //DataType_t data_type;
                //ElementType_t element_type;
                //ZoneType_t zone_type;
                //BCType_t bc_type;
                //GridLocation_t gridlocation_type;
                //PointSetType_t pointset_type;

                index_file = index_base = 1;

                if (cg.open(filePath.ToNullTermCharAry(), (int)Mode.MODE_READ, out index_file) != (int)Error.CG_OK) {
                    ThrowError(index_file, cg);
                }
                if (cg.nbases(index_file, out nbases) != (int)Error.CG_OK) {
                    ThrowError(index_file, cg);
                }
                if (nbases != 1) {
                    ThrowError("Only one base allowed", index_file, cg);
                }
                if (cg.base_read(index_file, index_base, basename, out celldim, out physdim) != (int)Error.CG_OK) {
                    ThrowError(index_file, cg);
                }
                if (cg.nzones(index_file, index_base, out nzones) != (int)Error.CG_OK) {
                    ThrowError(index_file, cg);
                }

                //x_of_all_zones = new double[nzones][, ,];
                //y_of_all_zones = new double[nzones][, ,];
                //z_of_all_zones = new double[nzones][, ,];
                //bcelements_of_all_zones = new int[nzones][][][];
                //normal_indices_of_all_zones = new int[nzones][][][];
                //starts_of_all_zones = new int[nzones][];
                //ends_of_all_zones = new int[nzones][];
                //nbndrys_of_all_zones = new int[nzones][];
                //pointset_types_of_all_zones = new PointSetType_t[nzones][][];
                //gridlocation_types_of_all_zones = new GridLocation_t[nzones][][];
                //bc_types_of_all_zones = new BCType_t[nzones][][];
                //element_types_of_all_zones = new ElementType_t[nzones][];
                //elements_of_all_zones = new int[nzones][][,];
                //boconames_of_all_zones = new char[nzones][][][];
                //zone_types_of_all_zones = new ZoneType_t[nzones];


                for (index_zone = 1; index_zone <= nzones; index_zone++) {
                    this.Zones.Add(new CgnsZone(this, index_file, index_base, index_zone, cg));
                }


                if (cg.close(index_file) != (int)Error.CG_OK) {
                    ThrowError(cg);
                }
            }
        }


        public IList<CgnsZone> Zones = new List<CgnsZone>();



        public class CgnsZone {

            public double[, ,] x = null;
            public double[, ,] y = null;
            public double[, ,] z = null;

            CgnsFile owner;


            public int[][][] bcelements = null;
            public int[][][] normal_indices = null;
            public int[] starts = null;
            public int[] ends = null;
            public int[] nbndrys = null;
            public PointSetType_t[][] pointset_types = null;
            public BCType_t[][] bc_types = null;
            public ElementType_t[] element_types = null;
            public char[][][] boconames = null;
            public GridLocation_t[][] gridlocation_types = null;
            public int[][,] elements = null;
            public ZoneType_t zone_type;


            internal CgnsZone(CgnsFile owna, int index_file, int index_base, int index_zone, CgnsDriver cg) {
                owner = owna;
                int index_section = 0;
                int index_bc = 0;
                int ncoords;
                int nsections;
                int nbocos;
                char[] basename = new char[100];
                char[] zonename = new char[100];
                char[] coordnameX = new char[100];
                char[] coordnameY = new char[100];
                char[] coordnameZ = new char[100];
                char[] flow = new char[100];
                char[] secname = new char[100];
                char[] boconame;
                char[] solname = new char[100];
                int[,] isize = new int[3, 3];
                int[,] ielem;
                int[] irmin = new int[3];
                int[] irmax = new int[3];
                int start, end, nbndry, parent_flag; // element definition
                int normal_list_flag, nbcelement, ndataset; // bc definition
                DataType_t data_type;
                ElementType_t element_type;

                BCType_t bc_type;
                GridLocation_t gridlocation_type;
                PointSetType_t pointset_type;

                //index_zone = 4; int zone_offset = index_zone - 1;
                //{

                if (cg.zone_type(index_file, index_base, index_zone, out zone_type) != (int)Error.CG_OK) {
                    owner.ThrowError(index_file, cg);
                }
                //zone_types_of_all_zones[index_zone - 1 - zone_offset] = zone_type;
                if (cg.zone_read(index_file, index_base, index_zone, zonename, isize) != (int)Error.CG_OK) {
                    owner.ThrowError(index_file, cg);
                }
                if (cg.ncoords(index_file, index_base, index_zone, out ncoords) != (int)Error.CG_OK) {
                    owner.ThrowError(index_file, cg);
                }
                if (cg.coord_info(index_file, index_base, index_zone, 1, out data_type, coordnameX) != (int)Error.CG_OK) {
                    owner.ThrowError(index_file, cg);
                }
                if (cg.coord_info(index_file, index_base, index_zone, 2, out data_type, coordnameY) != (int)Error.CG_OK) {
                    owner.ThrowError(index_file, cg);
                }
                if (owner.celldim > 2) {
                    if (cg.coord_info(index_file, index_base, index_zone, 3, out data_type, coordnameZ) != (int)Error.CG_OK) {
                        owner.ThrowError(index_file, cg);
                    }
                }

                irmin[0] = 1;
                irmin[1] = 1;
                irmin[2] = 1;
                irmax[0] = isize[0, 0];
                irmax[1] = isize[0, 1];
                irmax[2] = isize[0, 2];



                if ((int)zone_type == (int)ZoneType_t.Unstructured) {
                    x = new double[isize[0, 0], 1, 1];
                    y = new double[isize[0, 0], 1, 1];
                    z = new double[isize[0, 0], 1, 1];
                } else {
                    x = new double[isize[0, 2], isize[0, 1], isize[0, 0]];
                    y = new double[isize[0, 2], isize[0, 1], isize[0, 0]];
                    z = new double[isize[0, 2], isize[0, 1], isize[0, 0]];
                }

                char[][] boconametemp = null;
                int[][] bcelementtemp = null;
                int[] bcelement = null;
                int[][] normal_indextemp = null;
                int[] normal_index = null;
                PointSetType_t[] pointset_typtemp = null;
                GridLocation_t[] gridlocation_typtemp = null;
                BCType_t[] bc_typtemp = null;


                // read coordinates of this zone
                // =============================

                // x - coordinates
                if (cg.coord_read3(index_file, index_base, index_zone, coordnameX, data_type, irmin, irmax, x) != (int)Error.CG_OK) {
                    owner.ThrowError(index_file, cg);
                }

                // y - coordinates
                if (cg.coord_read3(index_file, index_base, index_zone, coordnameY, data_type, irmin, irmax, y) != (int)Error.CG_OK) {
                    owner.ThrowError(index_file, cg);
                }

                // optionally, if running in 3D, z - coordinates
                if (owner.celldim > 2) {
                    if (cg.coord_read3(index_file, index_base, index_zone, coordnameZ, data_type, irmin, irmax, z) != (int)Error.CG_OK) {
                        owner.ThrowError(index_file, cg);
                    }
                }
                if (cg.nsections(index_file, index_base, index_zone, out nsections) != (int)Error.CG_OK) {
                    owner.ThrowError(index_file, cg);
                }

                if ((int)zone_type == (int)ZoneType_t.Unstructured) {


                    //{
                    //    // fk test anfang
                    //    int nBoco = -1;
                    //    cg.nbocos(index_file, index_base, index_zone, out nBoco);

                    //    for (int iBoco = 1; iBoco <= nBoco; iBoco++) {
                    //        char[] _Boconame = new char[1000];
                    //        BCType_t BocoType;
                    //        PointSetType_t BocoPointsetType;
                    //        int nPnts;
                    //        int[] NormalIndex = new int[666]; ArrayTools.Set(NormalIndex, -666);
                    //        int NormalListFlag;
                    //        DataType_t NormalDataType;
                    //        int _ndataset;
                    //        cg.boco_info(index_file, index_base, index_zone, iBoco, _Boconame, out BocoType, out BocoPointsetType, out nPnts, NormalIndex,
                    //            out NormalListFlag, out NormalDataType, out _ndataset);

                    //        //string sBocoName = new string(_Boconame,);


                    //        Console.WriteLine();
                    //    }

                    //    // fk test ende
                    //}


                    boconames = new char[nsections][][];
                    bcelements = new int[nsections][][];
                    normal_indices = new int[nsections][][];
                    pointset_types = new PointSetType_t[nsections][];
                    gridlocation_types = new GridLocation_t[nsections][];
                    bc_types = new BCType_t[nsections][];
                    starts = new int[nsections];
                    ends = new int[nsections];
                    nbndrys = new int[nsections];
                    element_types = new ElementType_t[nsections];
                    elements = new int[nsections][,];


                    //



                    // Original Scheerer for reading bc
                    for (index_section = 1; index_section <= nsections; index_section++) {
                        if (cg.section_read(index_file, index_base, index_zone, index_section, secname, out element_type, out start, out end, out nbndry, out parent_flag) != (int)Error.CG_OK) {
                            owner.ThrowError(index_file, cg);
                        }

                        starts[index_section - 1] = start;
                        ends[index_section - 1] = end;
                        nbndrys[index_section - 1] = nbndry;
                        element_types[index_section - 1] = element_type;

                        char[] c = element_type.ToString().ToNullTermCharAry();

                        ielem = new int[end - start + 1, c[c.Length - 2] - 48];

                        if (cg.elements_read2(index_file, index_base, index_zone, index_section, ielem, null) != (int)Error.CG_OK) {
                            owner.ThrowError(index_file, cg);
                        }

                        elements[index_section - 1] = ielem;

                        if (cg.nbocos(index_file, index_base, index_zone, out nbocos) != (int)Error.CG_OK) {
                            owner.ThrowError(index_file, cg);
                        }
                        boconametemp = new char[nbocos][];
                        bcelementtemp = new int[nbocos][];
                        normal_indextemp = new int[nbocos][];
                        pointset_typtemp = new PointSetType_t[nbocos];
                        gridlocation_typtemp = new GridLocation_t[nbocos];
                        bc_typtemp = new BCType_t[nbocos];

                        int depht = 3;
                        char[] label0 = "Zone_t".ToNullTermCharAry(), label1 = "ZoneBC_t".ToNullTermCharAry(), label2 = "BC_t".ToNullTermCharAry();
                        string[] labels = new string[depht];
                        labels[0] = new String(label0);
                        labels[1] = new String(label1);
                        labels[2] = new String(label2);

                        int[] num = new int[depht];

                        num[0] = index_zone;
                        num[1] = 1;
                        for (index_bc = 1; index_bc <= nbocos; index_bc++) {
                            boconame = new char[100];
                            normal_index = new int[3];
                            if (cg.boco_info(index_file, index_base, index_zone, index_bc, boconame, out bc_type, out pointset_type, out nbcelement, normal_index, out normal_list_flag, out data_type, out ndataset) != (int)Error.CG_OK) {
                                owner.ThrowError(index_file, cg);
                            }
                            bcelement = new int[nbcelement];
                            if (cg.boco_read1(index_file, index_base, index_zone, index_bc, bcelement, null) != (int)Error.CG_OK) {
                                owner.ThrowError(index_file, cg);
                            }
                            num[2] = index_bc;
                            if (pointset_type.Equals(PointSetType_t.ElementList) || pointset_type.Equals(PointSetType_t.ElementRange)) {
                                gridlocation_type = GridLocation_t.FaceCenter;
                            } else {
                                if (cg.golist(index_file, index_base, depht, labels, num) != (int)Error.CG_OK) {
                                    owner.ThrowError(index_file, cg);
                                }
                                if (cg.gridlocation_read(out gridlocation_type) != (int)Error.CG_OK) {
                                    owner.ThrowError(index_file, cg);
                                }
                            }
                            boconametemp[index_bc - 1] = boconame;
                            bcelementtemp[index_bc - 1] = bcelement;
                            normal_indextemp[index_bc - 1] = normal_index;
                            pointset_typtemp[index_bc - 1] = pointset_type;
                            bc_typtemp[index_bc - 1] = bc_type;
                            gridlocation_typtemp[index_bc - 1] = gridlocation_type;
                        }
                        boconames[index_section - 1] = boconametemp;
                        bcelements[index_section - 1] = bcelementtemp;
                        normal_indices[index_section - 1] = normal_indextemp;
                        pointset_types[index_section - 1] = pointset_typtemp;
                        bc_types[index_section - 1] = bc_typtemp;
                        gridlocation_types[index_section - 1] = gridlocation_typtemp;
                    }
                } else {
                    //need_struct_conversion = true;
                    // insert here ...

                    throw new NotImplementedException("no structured CGNS grids supported until now");
                }



                //bcelements_of_all_zones[index_zone - 1 - zone_offset] = bcelements;
                //normal_indices_of_all_zones[index_zone - 1 - zone_offset] = normal_indices;
                //starts_of_all_zones[index_zone - 1 - zone_offset] = starts;
                //ends_of_all_zones[index_zone - 1 - zone_offset] = ends;
                //nbndrys_of_all_zones[index_zone - 1 - zone_offset] = nbndrys;
                //pointset_types_of_all_zones[index_zone - 1 - zone_offset] = pointset_types;
                //bc_types_of_all_zones[index_zone - 1 - zone_offset] = bc_types;
                //element_types_of_all_zones[index_zone - 1 - zone_offset] = element_types;
                //elements_of_all_zones[index_zone - 1 - zone_offset] = elements;
                //boconames_of_all_zones[index_zone - 1 - zone_offset] = boconames;
                //gridlocation_types_of_all_zones[index_zone - 1 - zone_offset] = gridlocation_types;

            }


            public GridCommons GenerateBoSSSGrid() {
                CgnsZoneToBoSSSGrid converter = new CgnsZoneToBoSSSGrid(this);
                return converter.GenerateBoSSSGrid();
            }


        }

        private void ThrowError(CgnsDriver cg) {
            throw new ApplicationException(cg.get_error());
        }

        private static void ThrowError(String s) {
            throw new ApplicationException(s);
        }

        private void ThrowError(int index_file, CgnsDriver cg) {
            String s = cg.get_error();

            if (cg.close(index_file) != (int)Error.CG_OK) {
                s = s + " " + cg.get_error();
            }

            throw new ApplicationException(s);
        }

        private void ThrowError(String s, int index_file, CgnsDriver cg) {
            if (cg.close(index_file) != (int)Error.CG_OK) {
                s = s + " " + cg.get_error();
            }

            throw new ApplicationException(s);
        }


        public GridCommons GenerateBoSSSGrid() {
            GridCommons ret = null;
            foreach (var zone in this.Zones) {
                GridCommons g = zone.GenerateBoSSSGrid();
                if (ret == null) {
                    ret = g;
                } else {
                    ret = GridCommons.MergeGrids(ret, g);
                    ret.MergeUnspecifiedEdges();
                }
            }
            return ret;
        }
    }

#pragma warning restore  1591
}
