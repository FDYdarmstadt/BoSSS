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
using MPI.Wrappers.Utils;
using System.Runtime.InteropServices;

namespace BoSSS.Solution.GridImport {
    /// <summary>
    /// some extension methods for ...
    /// </summary>
    public static class _666 {

        /// <summary>
        /// converts a string to a 0-terminated character (16-bit) array.
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
        public static char[] ToNullTermCharAry(this string s) {
            var c = s.ToCharArray();
            var ret = new char[c.Length + 1];
            Array.Copy(c, ret, c.Length);
            return ret;
        }

        /// <summary>
        /// converts a NULL-terminated character array int a string
        /// </summary>
        /// <param name="s"></param>
        /// <returns></returns>
        public static string FromNullTermCharAry(this char[] s) {
            int l;
            for (l = 0; l < s.Length; l++) {
                if (s[l] == '\0')
                    break;
            }
            if (l == s.Length)
                throw new OverflowException("found no NULL-termination");

            return new string(s, 0, l);
        }

        /// <summary>
        /// number of nodes for element of type <paramref name="t"/>.
        /// </summary>
        public static int NoOfNodes(this ElementType_t t) {
            switch (t) {
                case ElementType_t.ElementTypeNull: throw new NotSupportedException();
                case ElementType_t.ElementTypeUserDefined: throw new NotSupportedException();
                case ElementType_t.NODE: return 1;
                case ElementType_t.BAR_2: return 2;
                case ElementType_t.BAR_3: return 3; 			
                case ElementType_t.TRI_3: return 3;
                case ElementType_t.TRI_6: return 6;				
                case ElementType_t.QUAD_4: return 4;
                case ElementType_t.QUAD_8: return 8;
                case ElementType_t.QUAD_9: return 9;			
                case ElementType_t.TETRA_4: return 4;
                case ElementType_t.TETRA_10: return 10; 		
                case ElementType_t.PYRA_5: return 5;
                case ElementType_t.PYRA_14: return 14; 			
                case ElementType_t.PENTA_6: return 6;
                case ElementType_t.PENTA_15: return 15;
                case ElementType_t.PENTA_18: return 18;	
                case ElementType_t.HEXA_8: return 8;
                case ElementType_t.HEXA_20: return 20;
                case ElementType_t.HEXA_27: return 27;	
                case ElementType_t.MIXED: throw new NotSupportedException();
                case ElementType_t.NGON_n: throw new NotSupportedException();			
                default: throw new NotImplementedException();
            }
        }

        /// <summary>
        /// Spatial dimension of CGNS element <paramref name="t"/>
        /// </summary>
        public static int Dimension(this ElementType_t t) {
            switch (t) {
                case ElementType_t.ElementTypeNull: throw new NotSupportedException();
                case ElementType_t.ElementTypeUserDefined: throw new NotSupportedException();
                case ElementType_t.NODE: return 0;
                case ElementType_t.BAR_2: return 1;
                case ElementType_t.BAR_3: return 1;
                case ElementType_t.TRI_3: return 2;
                case ElementType_t.TRI_6: return 2;
                case ElementType_t.QUAD_4: return 2;
                case ElementType_t.QUAD_8: return 2;
                case ElementType_t.QUAD_9: return 2;
                case ElementType_t.TETRA_4: return 3;
                case ElementType_t.TETRA_10: return 3;
                case ElementType_t.PYRA_5: return 3;
                case ElementType_t.PYRA_14: return 3;
                case ElementType_t.PENTA_6: return 3;
                case ElementType_t.PENTA_15: return 3;
                case ElementType_t.PENTA_18: return 3;
                case ElementType_t.HEXA_8: return 3;
                case ElementType_t.HEXA_20: return 3;
                case ElementType_t.HEXA_27: return 3;
                case ElementType_t.MIXED: throw new NotSupportedException();
                case ElementType_t.NGON_n: throw new NotSupportedException();
                default: throw new NotImplementedException();
            }
        }

        
    }

#pragma warning disable        1591

    /// <summary>
    /// errors for cgns file;
    /// </summary>
    public enum Error {
        CG_OK = 0,
        CG_ERROR = 1,
        CG_NODE_NOT_FOUND = 2,
        CG_INCORRECT_PATH = 3,
        CG_NO_INDEX_DIM = 4
    };

    /// <summary>
    /// modes for cgns file;
    /// </summary>
    public enum Mode {
        MODE_READ = 0,
        MODE_WRITE = 1,
        MODE_CLOSED = 2,
        MODE_MODIFY = 3
    };

    /// <summary>
    /// Simulation types;
    /// </summary>
    public enum SimulationType_t {
        SimulationTypeNull,
        SimulationTypeUserDefined,
        TimeAccurate,
        NonTimeAccurate
    };

    /// <summary>
    /// Data types:  Can not add data types and stay forward compatible;
    /// </summary>
    public enum DataType_t {
        DataTypeNull,
        DataTypeUserDefined,
        Integer,
        RealSingle,
        RealDouble,
        Character
    };

    /// <summary>
    /// Zone types;
    /// </summary>
    public enum ZoneType_t {
        ZoneTypeNull,
        ZoneTypeUserDefined,
        Structured,
        Unstructured
    };

    /// <summary>
    /// Element types;
    /// </summary>
    public enum ElementType_t {
        ElementTypeNull,
        ElementTypeUserDefined,	/* 0, 1,	*/
        NODE,
        BAR_2,
        BAR_3, 				/* 2, 3, 4, 	*/
        TRI_3,
        TRI_6,					/* 5, 6,	*/
        QUAD_4,
        QUAD_8,
        QUAD_9,				/* 7, 8, 9,	*/
        TETRA_4,
        TETRA_10, 				/* 10, 11,	*/
        PYRA_5,
        PYRA_14, 				/* 12, 13,	*/
        PENTA_6,
        PENTA_15,
        PENTA_18,			/* 14, 15, 16,	*/
        HEXA_8,
        HEXA_20,
        HEXA_27, 			/* 17, 18, 19,	*/
        MIXED,
        NGON_n					/* 20, 21+	*/
    };

    /// <summary>
    /// Boundary Condition Types;
    /// </summary>
    public enum BCType_t {
        BCTypeNull,
        BCTypeUserDefined,
        BCAxisymmetricWedge,
        BCDegenerateLine,
        BCDegeneratePoint,
        BCDirichlet,
        BCExtrapolate,
        BCFarfield,
        BCGeneral,
        BCInflow,
        BCInflowSubsonic,
        BCInflowSupersonic,
        BCNeumann,
        BCOutflow,
        BCOutflowSubsonic,
        BCOutflowSupersonic,
        BCSymmetryPlane,
        BCSymmetryPolar,
        BCTunnelInflow,
        BCTunnelOutflow,
        BCWall,
        BCWallInviscid,
        BCWallViscous,
        BCWallViscousHeatFlux,
        BCWallViscousIsothermal,
        FamilySpecified
    };

    /// <summary>
    /// Point Set Types: Can't add types and stay forward compatible;
    /// </summary>
    public enum PointSetType_t {
        PointSetTypeNull,
        PointSetTypeUserDefined,
        PointList,
        PointListDonor,
        PointRange,
        PointRangeDonor,
        ElementRange,
        ElementList,
        CellListDonor
    };

    /// <summary>
    /// Grid Location;
    /// </summary>
    public enum GridLocation_t {
        GridLocationNull,
        GridLocationUserDefined,
        Vertex,
        CellCenter,
        FaceCenter,
        IFaceCenter,
        JFaceCenter,
        KFaceCenter,
        EdgeCenter
    };


    public class CgnsDriver : DynLibLoader {

        static string WinMangling(string _name) {
            string full = "BoSSS_cg_" + _name;

            return full.TrimEnd('1', '2', '3');
        }

        public CgnsDriver(bool hdf5) :
            base(hdf5 ? new string[] { "cgnsHdf5.dll" } : new string[] { "cgnsAdf.dll" },
                 new string[1][][],
                 new GetNameMangling[] { WinMangling },
                 new PlatformID[] { PlatformID.Win32NT },
                 new int[] { -1 }) {
        }

        public delegate int decl_cg_open(char[] filename, int mode, out int fn);
        public delegate int decl_cg_close(int fn);
        public delegate int decl_cg_golist(int fn, int B, int depth, string[] label, int[] num);
        public delegate int decl_cg_gopath(int fn, char[] path);
        public delegate int decl_cg_where(ref int fn, ref int B, out int depth, out IntPtr label, [Out] int[] num);
        public delegate int decl_cg_nbases(int fn, out int nbases);
        public delegate int decl_cg_base_write(int file_number, char[] basename, int cell_dim, int phys_dim, out int B);
        public delegate int decl_cg_base_read(int file_number, int B, [Out] char[] basename, out int cell_dim, out int phys_dim);
        public delegate int decl_cg_nzones(int fn, int B, out int nzones);
        public delegate int decl_cg_zone_type(int fn, int B, int Z, out ZoneType_t zonetype);
        public delegate int decl_cg_zone_write(int fn, int B, char[] zonename, int[,] size, ZoneType_t type, out int Z);
        public delegate int decl_cg_zone_read(int fn, int B, int Z, [Out] char[] zonename, [Out] int[,] size);
        public delegate int decl_cg_nsections(int fn, int B, int Z, out int nsections);
        public delegate int decl_cg_section_write3(int file_number, int B, int Z, char[] SectionName, ElementType_t type, int start, int end, int nbndry, int[, ,] elements, out int S);
        public delegate int decl_cg_section_write2(int file_number, int B, int Z, char[] SectionName, ElementType_t type, int start, int end, int nbndry, int[,] elements, out int S);
        public delegate int decl_cg_section_write1(int file_number, int B, int Z, char[] SectionName, ElementType_t type, int start, int end, int nbndry, int[] elements, out int S);
        public delegate int decl_cg_section_read(int file_number, int B, int Z, int S, [Out] char[] SectionName, out ElementType_t type, out int start, out int end, out int nbndry, out int parent_flag);
        public delegate int decl_cg_elements_read2(int fn, int B, int Z, int S, [Out] int[,] elements, [Out] int[] parent_data);
        public delegate int decl_cg_elements_read1(int fn, int B, int Z, int S, [Out] int[] elements, [Out] int[] parent_data);
        public delegate int decl_cg_gridlocation_write(GridLocation_t GridLocation);
        public delegate int decl_cg_gridlocation_read(out GridLocation_t GridLocation);
        public delegate int decl_cg_ngrids(int fn, int B, int Z, out int ngrids);
        public delegate int decl_cg_ncoords(int fn, int B, int Z, out int ncoords);
        public delegate int decl_cg_coord_write3(int fn, int B, int Z, DataType_t type, char[] coordname, double[, ,] coord_ptr, out int C);
        public delegate int decl_cg_coord_write2(int fn, int B, int Z, DataType_t type, char[] coordname, double[,] coord_ptr, out int C);
        public delegate int decl_cg_coord_write1(int fn, int B, int Z, DataType_t type, char[] coordname, double[] coord_ptr, out int C);
        public delegate int decl_cg_coord_read3(int fn, int B, int Z, char[] coordname, DataType_t type, int[] rmin, int[] rmax, [Out] double[, ,] coord);
        public delegate int decl_cg_coord_read2(int fn, int B, int Z, char[] coordname, DataType_t type, int[] rmin, int[] rmax, [Out] double[,] coord);
        public delegate int decl_cg_coord_read1(int fn, int B, int Z, char[] coordname, DataType_t type, int[] rmin, int[] rmax, [Out] double[] coord);
        public delegate int decl_cg_coord_info(int fn, int B, int Z, int C, out DataType_t type, [Out] char[] coordname);
        public delegate int decl_cg_nsols(int fn, int B, int Z, out int nsols);
        public delegate int decl_cg_sol_info(int fn, int B, int Z, int S, [Out] char[] solname, out GridLocation_t location);
        public delegate int decl_cg_sol_write(int fn, int B, int Z, char[] solname, GridLocation_t location, out int S);
        public delegate int decl_cg_nfields(int fn, int B, int Z, int S, out int nfields);
        public delegate int decl_cg_field_info(int fn, int B, int Z, int S, int F, out DataType_t type, [Out] char[] fieldname);
        public delegate int decl_cg_field_write3(int fn, int B, int Z, int S, DataType_t type, char[] fieldname, double[, ,] field_ptr, out int F);
        public delegate int decl_cg_field_write2(int fn, int B, int Z, int S, DataType_t type, char[] fieldname, double[,] field_ptr, out int F);
        public delegate int decl_cg_field_write1(int fn, int B, int Z, int S, DataType_t type, char[] fieldname, double[] field_ptr, out int F);
        public delegate int decl_cg_field_read3(int fn, int B, int Z, int S, char[] fieldname, DataType_t type, int[] rmin, int[] rmax, [Out] double[, ,] field_ptr);
        public delegate int decl_cg_field_read2(int fn, int B, int Z, int S, char[] fieldname, DataType_t type, int[] rmin, int[] rmax, [Out] double[,] field_ptr);
        public delegate int decl_cg_field_read1(int fn, int B, int Z, int S, char[] fieldname, DataType_t type, int[] rmin, int[] rmax, [Out] double[] field_ptr);
        public delegate int decl_cg_nbocos(int fn, int B, int Z, out int nbocos);
        public delegate int decl_cg_boco_write3(int file_number, int B, int Z, char[] boconame, BCType_t bocotype, PointSetType_t ptset_type, int npnts, int[, ,] pnts, out int BC);
        public delegate int decl_cg_boco_write2(int file_number, int B, int Z, char[] boconame, BCType_t bocotype, PointSetType_t ptset_type, int npnts, int[,] pnts, out int BC);
        public delegate int decl_cg_boco_write1(int file_number, int B, int Z, char[] boconame, BCType_t bocotype, PointSetType_t ptset_type, int npnts, int[] pnts, out int BC);
        public delegate int decl_cg_boco_read3(int fn, int B, int Z, int BC, [Out] int[, ,] pnts, [Out] double[, ,] NormalList);
        public delegate int decl_cg_boco_read2(int fn, int B, int Z, int BC, [Out] int[,] pnts, [Out] double[, ,] NormalList);
        public delegate int decl_cg_boco_read1(int fn, int B, int Z, int BC, [Out] int[] pnts, [Out] double[, ,] NormalList);
        public delegate int decl_cg_boco_info(int fn, int B, int Z, int BC, [Out] char[] boconame, out BCType_t bocotype, out PointSetType_t ptset_type, out int npnts, [Out] int[] NormalIndex, out int NormalListFlag, out DataType_t NormalDataType, out int ndataset);
        public delegate int decl_cg_ndescriptors(out int ndescriptors);
        public delegate int decl_cg_descriptor_write(char[] descr_name, char[] descr_text);
        public delegate int decl_cg_descriptor_read(int D, [Out] char[] descr_name, [Out] char[,] descr_text);
        public delegate int decl_cg_simulation_type_write(int file_number, int B, SimulationType_t type);
        public delegate int decl_cg_simulation_type_read(int file_number, int B, out SimulationType_t type);
        public delegate int decl_cg_biter_write(int file_number, int B, char[] bitername, int nsteps);
        public delegate int decl_cg_biter_read(int file_number, int B, [Out] char[] bitername, out int nsteps);
        public delegate int decl_cg_ziter_write(int file_number, int B, int Z, char[] zitername);
        public delegate int decl_cg_ziter_read(int file_number, int B, int Z, [Out] char[] zitername);
        public delegate int decl_cg_rind_write(int[] rind_data);
        public delegate int decl_cg_rind_read([Out] int[] rind_data);
        public delegate String decl_cg_get_error();

        public decl_cg_open open;
        public decl_cg_close close;
        public decl_cg_golist golist;
        public decl_cg_gopath gopath;
        public decl_cg_where where;
        public decl_cg_nbases nbases;
        public decl_cg_base_write base_write;
        public decl_cg_base_read base_read;
        public decl_cg_nzones nzones;
        public decl_cg_zone_type zone_type;
        public decl_cg_zone_write zone_write;
        public decl_cg_zone_read zone_read;
        public decl_cg_nsections nsections;
        public decl_cg_section_write3 section_write3;
        public decl_cg_section_write2 section_write2;
        public decl_cg_section_write1 section_write1;
        public decl_cg_section_read section_read;
        public decl_cg_elements_read2 elements_read2;
        public decl_cg_elements_read1 elements_read1;
        public decl_cg_gridlocation_write gridlocation_write;
        public decl_cg_gridlocation_read gridlocation_read;
        public decl_cg_ngrids ngrids;
        public decl_cg_ncoords ncoords;
        public decl_cg_coord_write3 coord_write3;
        public decl_cg_coord_write2 coord_write2;
        public decl_cg_coord_write1 coord_write1;
        public decl_cg_coord_read3 coord_read3;
        public decl_cg_coord_read2 coord_read2;
        public decl_cg_coord_read1 coord_read1;
        public decl_cg_coord_info coord_info;
        public decl_cg_nsols nsols;
        public decl_cg_sol_info sol_info;
        public decl_cg_sol_write sol_write;
        public decl_cg_nfields nfields;
        public decl_cg_field_info field_info;
        public decl_cg_field_write3 field_write3;
        public decl_cg_field_write2 field_write2;
        public decl_cg_field_write1 field_write1;
        public decl_cg_field_read3 field_read3;
        public decl_cg_field_read2 field_read2;
        public decl_cg_field_read1 field_read1;
        public decl_cg_nbocos nbocos;
        public decl_cg_boco_write3 boco_write3;
        public decl_cg_boco_write2 boco_write2;
        public decl_cg_boco_write1 boco_write1;
        public decl_cg_boco_read3 boco_read3;
        public decl_cg_boco_read2 boco_read2;
        public decl_cg_boco_read1 boco_read1;
        public decl_cg_boco_info boco_info;
        public decl_cg_ndescriptors ndescriptors;
        public decl_cg_descriptor_write descriptor_write;
        public decl_cg_descriptor_read descriptor_read;
        public decl_cg_simulation_type_write simulation_type_write;
        public decl_cg_simulation_type_read simulation_type_read;
        public decl_cg_biter_write biter_write;
        public decl_cg_biter_read biter_read;
        public decl_cg_ziter_write ziter_write;
        public decl_cg_ziter_read ziter_read;
        public decl_cg_rind_write rind_write;
        public decl_cg_rind_read rind_read;
        public decl_cg_get_error get_error;

    }
#pragma warning restore  1591

}
