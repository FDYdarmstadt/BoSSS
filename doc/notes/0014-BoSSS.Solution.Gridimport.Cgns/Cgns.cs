using System;
using System.Diagnostics;
using System.Collections.Generic;
using System.Collections;
using System.Runtime.InteropServices;
using System.Text.RegularExpressions;
using System.Text;

using BoSSS.Platform;
using BoSSS.Foundation;
using BoSSS.Solution;
using BoSSS.Foundation.Grid;

namespace BoSSS.Solution.GridImport {

    /// <summary>
    /// errors for cgns file;
    /// </summary>
    enum Error {
        CG_OK = 0,
        CG_ERROR = 1,
        CG_NODE_NOT_FOUND = 2,
        CG_INCORRECT_PATH = 3,
        CG_NO_INDEX_DIM = 4
    };

    /// <summary>
    /// modes for cgns file;
    /// </summary>
    enum Mode {
        MODE_READ = 0,
        MODE_WRITE = 1,
        MODE_CLOSED = 2,
        MODE_MODIFY = 3
    };

    /// <summary>
    /// Simulation types;
    /// </summary>
    enum SimulationType_t {
        SimulationTypeNull,
        SimulationTypeUserDefined,
        TimeAccurate,
        NonTimeAccurate
    };

    /// <summary>
    /// Data types:  Can not add data types and stay forward compatible;
    /// </summary>
    enum DataType_t {
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
    enum ZoneType_t {
        ZoneTypeNull,
        ZoneTypeUserDefined,
        Structured,
        Unstructured
    };

    /// <summary>
    /// Element types;
    /// </summary>
    enum ElementType_t {
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
    enum BCType_t {
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
    enum PointSetType_t {
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
    enum GridLocation_t {
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

    /// <summary>
    /// a cgns file (at least the sections that are important for BoSSS)
    /// </summary>
    public class Cgns : PlotDriver {

        [DllImport("kernel32.dll", SetLastError = true, CallingConvention = CallingConvention.Winapi)]
        [return: MarshalAs(UnmanagedType.Bool)]
        public static extern bool IsWow64Process([In] IntPtr hProcess, [Out] out bool lpSystemInfo);

        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_open")]
        private static extern int Adf_cg_open(char[] filename, int mode, out int fn);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_close")]
        private static extern int Adf_cg_close(int fn);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_gopath")]
        private static extern int Adf_cg_gopath(int fn, char[] path);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_golist")]
        private static extern int Adf_cg_golist(int fn, int B, int depth, string[] label, int[] num);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_where")]
        private static extern int Adf_cg_where(ref int fn, ref int B, out int depth, out IntPtr label, [Out] int[] num);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_nbases")]
        private static extern int Adf_cg_nbases(int fn, out int nbases);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_base_write")]
        private static extern int Adf_cg_base_write(int file_number, char[] basename, int cell_dim, int phys_dim, out int B);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_base_read")]
        private static extern int Adf_cg_base_read(int file_number, int B, [Out] char[] basename, out int cell_dim, out int phys_dim);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_nzones")]
        private static extern int Adf_cg_nzones(int fn, int B, out int nzones);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_zone_type")]
        private static extern int Adf_cg_zone_type(int fn, int B, int Z, out ZoneType_t zonetype);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_zone_write")]
        private static extern int Adf_cg_zone_write(int fn, int B, char[] zonename, int[,] size, ZoneType_t type, out int Z);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_zone_read")]
        private static extern int Adf_cg_zone_read(int fn, int B, int Z, [Out] char[] zonename, [Out] int[,] size);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_nsections")]
        private static extern int Adf_cg_nsections(int fn, int B, int Z, out int nsections);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_section_write")]
        private static extern int Adf_cg_section_write2(int file_number, int B, int Z, char[] SectionName, ElementType_t type, int start, int end, int nbndry, int[,] elements, out int S);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_section_write")]
        private static extern int Adf_cg_section_write1(int file_number, int B, int Z, char[] SectionName, ElementType_t type, int start, int end, int nbndry, int[] elements, out int S);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_section_read")]
        private static extern int Adf_cg_section_read(int file_number, int B, int Z, int S, [Out] char[] SectionName, out ElementType_t type, out int start, out int end, out int nbndry, out int parent_flag);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_elements_read")]
        private static extern int Adf_cg_elements_read2(int fn, int B, int Z, int S, [Out] int[,] elements, [Out] int[] parent_data);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_elements_read")]
        private static extern int Adf_cg_elements_read1(int fn, int B, int Z, int S, [Out] int[] elements, [Out] int[] parent_data);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_ngrids")]
        private static extern int Adf_cg_ngrids(int fn, int B, int Z, out int ngrids);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_gridlocation_write")]
        private static extern int Adf_cg_gridlocation_write(GridLocation_t GridLocation);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_gridlocation_read")]
        private static extern int Adf_cg_gridlocation_read(out GridLocation_t GridLocation);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_ncoords")]
        private static extern int Adf_cg_ncoords(int fn, int B, int Z, out int ncoords);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_coord_write")]
        private static extern int Adf_cg_coord_write3(int fn, int B, int Z, DataType_t type, char[] coordname, double[, ,] coord_ptr, out int C);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_coord_read")]
        private static extern int Adf_cg_coord_write2(int fn, int B, int Z, DataType_t type, char[] coordname, double[,] coord_ptr, out int C);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_coord_read")]
        private static extern int Adf_cg_coord_write1(int fn, int B, int Z, DataType_t type, char[] coordname, double[] coord_ptr, out int C);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_coord_read")]
        private static extern int Adf_cg_coord_read3(int fn, int B, int Z, char[] coordname, DataType_t type, int[] rmin, int[] rmax, [Out] double[, ,] coord);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_coord_write")]
        private static extern int Adf_cg_coord_read2(int fn, int B, int Z, char[] coordname, DataType_t type, int[] rmin, int[] rmax, [Out] double[,] coord);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_coord_write")]
        private static extern int Adf_cg_coord_read1(int fn, int B, int Z, char[] coordname, DataType_t type, int[] rmin, int[] rmax, [Out] double[] coord);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_coord_info")]
        private static extern int Adf_cg_coord_info(int fn, int B, int Z, int C, out DataType_t type, [Out] char[] coordname);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_nsols")]
        private static extern int Adf_cg_nsols(int fn, int B, int Z, out int nsols);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_sol_info")]
        private static extern int Adf_cg_sol_info(int fn, int B, int Z, int S, [Out] char[] solname, out GridLocation_t location);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_sol_write")]
        private static extern int Adf_cg_sol_write(int fn, int B, int Z, char[] solname, GridLocation_t location, out int S);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_sol_nfields")]
        private static extern int Adf_cg_nfields(int fn, int B, int Z, int S, out int nfields);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_field_info")]
        private static extern int Adf_cg_field_info(int fn, int B, int Z, int S, int F, out DataType_t type, [Out] char[] fieldname);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_field_write")]
        private static extern int Adf_cg_field_write3(int fn, int B, int Z, int S, DataType_t type, char[] fieldname, double[, ,] field_ptr, out int F);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_field_write")]
        private static extern int Adf_cg_field_write2(int fn, int B, int Z, int S, DataType_t type, char[] fieldname, double[,] field_ptr, out int F);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_field_write")]
        private static extern int Adf_cg_field_write1(int fn, int B, int Z, int S, DataType_t type, char[] fieldname, double[] field_ptr, out int F);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_field_read")]
        private static extern int Adf_cg_field_read3(int fn, int B, int Z, int S, char[] fieldname, DataType_t type, int[] rmin, int[] rmax, [Out] double[, ,] field_ptr);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_field_read")]
        private static extern int Adf_cg_field_read2(int fn, int B, int Z, int S, char[] fieldname, DataType_t type, int[] rmin, int[] rmax, [Out] double[,] field_ptr);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_field_read")]
        private static extern int Adf_cg_field_read1(int fn, int B, int Z, int S, char[] fieldname, DataType_t type, int[] rmin, int[] rmax, [Out] double[] field_ptr);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_nbocos")]
        private static extern int Adf_cg_nbocos(int fn, int B, int Z, out int nbocos);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_boco_write")]
        private static extern int Adf_cg_boco_write3(int file_number, int B, int Z, char[] boconame, BCType_t bocotype, PointSetType_t ptset_type, int npnts, int[, ,] pnts, out int BC);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_boco_write")]
        private static extern int Adf_cg_boco_write2(int file_number, int B, int Z, char[] boconame, BCType_t bocotype, PointSetType_t ptset_type, int npnts, int[,] pnts, out int BC);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_boco_write")]
        private static extern int Adf_cg_boco_write1(int file_number, int B, int Z, char[] boconame, BCType_t bocotype, PointSetType_t ptset_type, int npnts, int[] pnts, out int BC);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_boco_read")]
        private static extern int Adf_cg_boco_read3(int fn, int B, int Z, int BC, [Out] int[, ,] pnts, [Out] double[, ,] NormalList);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_boco_read")]
        private static extern int Adf_cg_boco_read2(int fn, int B, int Z, int BC, [Out] int[,] pnts, [Out] double[, ,] NormalList);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_boco_read")]
        private static extern int Adf_cg_boco_read1(int fn, int B, int Z, int BC, [Out] int[] pnts, [Out] double[, ,] NormalList);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_boco_info")]
        private static extern int Adf_cg_boco_info(int fn, int B, int Z, int BC, [Out] char[] boconame, out BCType_t bocotype, out PointSetType_t ptset_type, out int npnts, [Out] int[] NormalIndex, out int NormalListFlag, out DataType_t NormalDataType, out int ndataset);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_ndescriptors")]
        private static extern int Adf_cg_ndescriptors(out int ndescriptors);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_descriptor_write")]
        private static extern int Adf_cg_descriptor_write(char[] descr_name, char[] descr_text);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_descriptor_read")]
        private static extern int Adf_cg_descriptor_read(int D, [Out] char[] descr_name, [Out] char[,] descr_text);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_simulation_type_write")]
        private static extern int Adf_cg_simulation_type_write(int file_number, int B, SimulationType_t type);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_simulation_type_read")]
        private static extern int Adf_cg_simulation_type_read(int file_number, int B, out SimulationType_t type);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_biter_write")]
        private static extern int Adf_cg_biter_write(int file_number, int B, char[] bitername, int nsteps);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_biter_read")]
        private static extern int Adf_cg_biter_read(int file_number, int B, [Out] char[] bitername, out int nsteps);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_ziter_write")]
        private static extern int Adf_cg_ziter_write(int file_number, int B, int Z, char[] zitername);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_ziter_read")]
        private static extern int Adf_cg_ziter_read(int file_number, int B, int Z, [Out] char[] zitername);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_rind_write")]
        private static extern int Adf_cg_rind_write(int[] rind_data);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_rind_read")]
        private static extern int Adf_cg_rind_read([Out] int[] rind_data);
        [DllImport("cgnsAdf.dll", EntryPoint = "BoSSS_cg_get_error")]
        private static extern String Adf_cg_get_error();

        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_open")]
        private static extern int Hdf5_cg_open(char[] filename, int mode, out int fn);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_close")]
        private static extern int Hdf5_cg_close(int fn);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_gopath")]
        private static extern int Hdf5_cg_gopath(int fn, char[] path);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_golist")]
        private static extern int Hdf5_cg_golist(int fn, int B, int depth, string[] label, int[] num);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_where")]
        private static extern int Hdf5_cg_where(ref int fn, ref int B, out int depth, out IntPtr label, [Out] int[] num);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_nbases")]
        private static extern int Hdf5_cg_nbases(int fn, out int nbases);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_base_write")]
        private static extern int Hdf5_cg_base_write(int file_number, char[] basename, int cell_dim, int phys_dim, out int B);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_base_read")]
        private static extern int Hdf5_cg_base_read(int file_number, int B, [Out] char[] basename, out int cell_dim, out int phys_dim);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_nzones")]
        private static extern int Hdf5_cg_nzones(int fn, int B, out int nzones);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_zone_type")]
        private static extern int Hdf5_cg_zone_type(int fn, int B, int Z, out ZoneType_t zonetype);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_zone_write")]
        private static extern int Hdf5_cg_zone_write(int fn, int B, char[] zonename, int[,] size, ZoneType_t type, out int Z);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_zone_read")]
        private static extern int Hdf5_cg_zone_read(int fn, int B, int Z, [Out] char[] zonename, [Out] int[,] size);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_nsections")]
        private static extern int Hdf5_cg_nsections(int fn, int B, int Z, out int nsections);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_section_write")]
        private static extern int Hdf5_cg_section_write2(int file_number, int B, int Z, char[] SectionName, ElementType_t type, int start, int end, int nbndry, int[,] elements, out int S);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_section_write")]
        private static extern int Hdf5_cg_section_write1(int file_number, int B, int Z, char[] SectionName, ElementType_t type, int start, int end, int nbndry, int[] elements, out int S);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_section_read")]
        private static extern int Hdf5_cg_section_read(int file_number, int B, int Z, int S, [Out] char[] SectionName, out ElementType_t type, out int start, out int end, out int nbndry, out int parent_flag);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_elements_read")]
        private static extern int Hdf5_cg_elements_read2(int fn, int B, int Z, int S, [Out] int[,] elements, [Out] int[] parent_data);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_elements_read")]
        private static extern int Hdf5_cg_elements_read1(int fn, int B, int Z, int S, [Out] int[] elements, [Out] int[] parent_data);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_gridlocation_write")]
        private static extern int Hdf5_cg_gridlocation_write(GridLocation_t GridLocation);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_gridlocation_read")]
        private static extern int Hdf5_cg_gridlocation_read(out GridLocation_t GridLocation);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_ncoords")]
        private static extern int Hdf5_cg_ngrids(int fn, int B, int Z, out int ngrids);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_ncoords")]
        private static extern int Hdf5_cg_ncoords(int fn, int B, int Z, out int ncoords);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_coord_write")]
        private static extern int Hdf5_cg_coord_write3(int fn, int B, int Z, DataType_t type, char[] coordname, double[, ,] coord_ptr, out int C);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_coord_read")]
        private static extern int Hdf5_cg_coord_write2(int fn, int B, int Z, DataType_t type, char[] coordname, double[,] coord_ptr, out int C);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_coord_read")]
        private static extern int Hdf5_cg_coord_write1(int fn, int B, int Z, DataType_t type, char[] coordname, double[] coord_ptr, out int C);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_coord_read")]
        private static extern int Hdf5_cg_coord_read3(int fn, int B, int Z, char[] coordname, DataType_t type, int[] rmin, int[] rmax, [Out] double[, ,] coord);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_coord_write")]
        private static extern int Hdf5_cg_coord_read2(int fn, int B, int Z, char[] coordname, DataType_t type, int[] rmin, int[] rmax, [Out] double[,] coord);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_coord_write")]
        private static extern int Hdf5_cg_coord_read1(int fn, int B, int Z, char[] coordname, DataType_t type, int[] rmin, int[] rmax, [Out] double[] coord);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_coord_info")]
        private static extern int Hdf5_cg_coord_info(int fn, int B, int Z, int C, out DataType_t type, [Out] char[] coordname);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_nsols")]
        private static extern int Hdf5_cg_nsols(int fn, int B, int Z, out int nsols);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_sol_info")]
        private static extern int Hdf5_cg_sol_info(int fn, int B, int Z, int S, [Out] char[] solname, out GridLocation_t location);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_sol_write")]
        private static extern int Hdf5_cg_sol_write(int fn, int B, int Z, char[] solname, GridLocation_t location, out int S);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_sol_nfields")]
        private static extern int Hdf5_cg_nfields(int fn, int B, int Z, int S, out int nfields);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_field_info")]
        private static extern int Hdf5_cg_field_info(int fn, int B, int Z, int S, int F, out DataType_t type, [Out] char[] fieldname);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_field_write")]
        private static extern int Hdf5_cg_field_write3(int fn, int B, int Z, int S, DataType_t type, char[] fieldname, double[, ,] field_ptr, out int F);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_field_write")]
        private static extern int Hdf5_cg_field_write2(int fn, int B, int Z, int S, DataType_t type, char[] fieldname, double[,] field_ptr, out int F);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_field_write")]
        private static extern int Hdf5_cg_field_write1(int fn, int B, int Z, int S, DataType_t type, char[] fieldname, double[] field_ptr, out int F);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_field_read")]
        private static extern int Hdf5_cg_field_read3(int fn, int B, int Z, int S, char[] fieldname, DataType_t type, int[] rmin, int[] rmax, [Out] double[, ,] field_ptr);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_field_read")]
        private static extern int Hdf5_cg_field_read2(int fn, int B, int Z, int S, char[] fieldname, DataType_t type, int[] rmin, int[] rmax, [Out] double[,] field_ptr);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_field_read")]
        private static extern int Hdf5_cg_field_read1(int fn, int B, int Z, int S, char[] fieldname, DataType_t type, int[] rmin, int[] rmax, [Out] double[] field_ptr);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_nbocos")]
        private static extern int Hdf5_cg_nbocos(int fn, int B, int Z, out int nbocos);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_boco_write")]
        private static extern int Hdf5_cg_boco_write3(int file_number, int B, int Z, char[] boconame, BCType_t bocotype, PointSetType_t ptset_type, int npnts, int[, ,] pnts, out int BC);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_boco_write")]
        private static extern int Hdf5_cg_boco_write2(int file_number, int B, int Z, char[] boconame, BCType_t bocotype, PointSetType_t ptset_type, int npnts, int[,] pnts, out int BC);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_boco_write")]
        private static extern int Hdf5_cg_boco_write1(int file_number, int B, int Z, char[] boconame, BCType_t bocotype, PointSetType_t ptset_type, int npnts, int[] pnts, out int BC);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_boco_read")]
        private static extern int Hdf5_cg_boco_read3(int fn, int B, int Z, int BC, [Out] int[, ,] pnts, [Out] double[, ,] NormalList);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_boco_read")]
        private static extern int Hdf5_cg_boco_read2(int fn, int B, int Z, int BC, [Out] int[,] pnts, [Out] double[, ,] NormalList);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_boco_read")]
        private static extern int Hdf5_cg_boco_read1(int fn, int B, int Z, int BC, [Out] int[] pnts, [Out] double[, ,] NormalList);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_boco_info")]
        private static extern int Hdf5_cg_boco_info(int fn, int B, int Z, int BC, [Out] char[] boconame, out BCType_t bocotype, out PointSetType_t ptset_type, out int npnts, [Out] int[] NormalIndex, out int NormalListFlag, out DataType_t NormalDataType, out int ndataset);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_ndescriptors")]
        private static extern int Hdf5_cg_ndescriptors(out int ndescriptors);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_descriptor_write")]
        private static extern int Hdf5_cg_descriptor_write(char[] descr_name, char[] descr_text);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_descriptor_read")]
        private static extern int Hdf5_cg_descriptor_read(int D, [Out] char[] descr_name, [Out] char[,] descr_text);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_simulation_type_write")]
        private static extern int Hdf5_cg_simulation_type_write(int file_number, int B, SimulationType_t type);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_simulation_type_read")]
        private static extern int Hdf5_cg_simulation_type_read(int file_number, int B, out SimulationType_t type);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_biter_write")]
        private static extern int Hdf5_cg_biter_write(int file_number, int B, char[] bitername, int nsteps);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_biter_read")]
        private static extern int Hdf5_cg_biter_read(int file_number, int B, [Out] char[] bitername, out int nsteps);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_ziter_write")]
        private static extern int Hdf5_cg_ziter_write(int file_number, int B, int Z, char[] zitername);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_ziter_read")]
        private static extern int Hdf5_cg_ziter_read(int file_number, int B, int Z, [Out] char[] zitername);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_rind_write")]
        private static extern int Hdf5_cg_rind_write(int[] rind_data);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_rind_read")]
        private static extern int Hdf5_cg_rind_read([Out] int[] rind_data);
        [DllImport("cgnsHdf5.dll", EntryPoint = "BoSSS_cg_get_error")]
        private static extern String Hdf5_cg_get_error();

        private delegate int decl_cg_open(char[] filename, int mode, out int fn);
        private delegate int decl_cg_close(int fn);
        private delegate int decl_cg_golist(int fn, int B, int depth, string[] label, int[] num);
        private delegate int decl_cg_gopath(int fn, char[] path);
        private delegate int decl_cg_where(ref int fn, ref int B, out int depth, out IntPtr label, [Out] int[] num);
        private delegate int decl_cg_nbases(int fn, out int nbases);
        private delegate int decl_cg_base_write(int file_number, char[] basename, int cell_dim, int phys_dim, out int B);
        private delegate int decl_cg_base_read(int file_number, int B, [Out] char[] basename, out int cell_dim, out int phys_dim);
        private delegate int decl_cg_nzones(int fn, int B, out int nzones);
        private delegate int decl_cg_zone_type(int fn, int B, int Z, out ZoneType_t zonetype);
        private delegate int decl_cg_zone_write(int fn, int B, char[] zonename, int[,] size, ZoneType_t type, out int Z);
        private delegate int decl_cg_zone_read(int fn, int B, int Z, [Out] char[] zonename, [Out] int[,] size);
        private delegate int decl_cg_nsections(int fn, int B, int Z, out int nsections);
        private delegate int decl_cg_section_write2(int file_number, int B, int Z, char[] SectionName, ElementType_t type, int start, int end, int nbndry, int[,] elements, out int S);
        private delegate int decl_cg_section_write1(int file_number, int B, int Z, char[] SectionName, ElementType_t type, int start, int end, int nbndry, int[] elements, out int S);
        private delegate int decl_cg_section_read(int file_number, int B, int Z, int S, [Out] char[] SectionName, out ElementType_t type, out int start, out int end, out int nbndry, out int parent_flag);
        private delegate int decl_cg_elements_read2(int fn, int B, int Z, int S, [Out] int[,] elements, [Out] int[] parent_data);
        private delegate int decl_cg_elements_read1(int fn, int B, int Z, int S, [Out] int[] elements, [Out] int[] parent_data);
        private delegate int decl_cg_gridlocation_write(GridLocation_t GridLocation);
        private delegate int decl_cg_gridlocation_read(out GridLocation_t GridLocation);
        private delegate int decl_cg_ngrids(int fn, int B, int Z, out int ngrids);
        private delegate int decl_cg_ncoords(int fn, int B, int Z, out int ncoords);
        private delegate int decl_cg_coord_write3(int fn, int B, int Z, DataType_t type, char[] coordname, double[, ,] coord_ptr, out int C);
        private delegate int decl_cg_coord_write2(int fn, int B, int Z, DataType_t type, char[] coordname, double[,] coord_ptr, out int C);
        private delegate int decl_cg_coord_write1(int fn, int B, int Z, DataType_t type, char[] coordname, double[] coord_ptr, out int C);
        private delegate int decl_cg_coord_read3(int fn, int B, int Z, char[] coordname, DataType_t type, int[] rmin, int[] rmax, [Out] double[, ,] coord);
        private delegate int decl_cg_coord_read2(int fn, int B, int Z, char[] coordname, DataType_t type, int[] rmin, int[] rmax, [Out] double[,] coord);
        private delegate int decl_cg_coord_read1(int fn, int B, int Z, char[] coordname, DataType_t type, int[] rmin, int[] rmax, [Out] double[] coord);
        private delegate int decl_cg_coord_info(int fn, int B, int Z, int C, out DataType_t type, [Out] char[] coordname);
        private delegate int decl_cg_nsols(int fn, int B, int Z, out int nsols);
        private delegate int decl_cg_sol_info(int fn, int B, int Z, int S, [Out] char[] solname, out GridLocation_t location);
        private delegate int decl_cg_sol_write(int fn, int B, int Z, char[] solname, GridLocation_t location, out int S);
        private delegate int decl_cg_nfields(int fn, int B, int Z, int S, out int nfields);
        private delegate int decl_cg_field_info(int fn, int B, int Z, int S, int F, out DataType_t type, [Out] char[] fieldname);
        private delegate int decl_cg_field_write3(int fn, int B, int Z, int S, DataType_t type, char[] fieldname, double[, ,] field_ptr, out int F);
        private delegate int decl_cg_field_write2(int fn, int B, int Z, int S, DataType_t type, char[] fieldname, double[,] field_ptr, out int F);
        private delegate int decl_cg_field_write1(int fn, int B, int Z, int S, DataType_t type, char[] fieldname, double[] field_ptr, out int F);
        private delegate int decl_cg_field_read3(int fn, int B, int Z, int S, char[] fieldname, DataType_t type, int[] rmin, int[] rmax, [Out] double[, ,] field_ptr);
        private delegate int decl_cg_field_read2(int fn, int B, int Z, int S, char[] fieldname, DataType_t type, int[] rmin, int[] rmax, [Out] double[,] field_ptr);
        private delegate int decl_cg_field_read1(int fn, int B, int Z, int S, char[] fieldname, DataType_t type, int[] rmin, int[] rmax, [Out] double[] field_ptr);
        private delegate int decl_cg_nbocos(int fn, int B, int Z, out int nbocos);
        private delegate int decl_cg_boco_write3(int file_number, int B, int Z, char[] boconame, BCType_t bocotype, PointSetType_t ptset_type, int npnts, int[, ,] pnts, out int BC);
        private delegate int decl_cg_boco_write2(int file_number, int B, int Z, char[] boconame, BCType_t bocotype, PointSetType_t ptset_type, int npnts, int[,] pnts, out int BC);
        private delegate int decl_cg_boco_write1(int file_number, int B, int Z, char[] boconame, BCType_t bocotype, PointSetType_t ptset_type, int npnts, int[] pnts, out int BC);
        private delegate int decl_cg_boco_read3(int fn, int B, int Z, int BC, [Out] int[, ,] pnts, [Out] double[, ,] NormalList);
        private delegate int decl_cg_boco_read2(int fn, int B, int Z, int BC, [Out] int[,] pnts, [Out] double[, ,] NormalList);
        private delegate int decl_cg_boco_read1(int fn, int B, int Z, int BC, [Out] int[] pnts, [Out] double[, ,] NormalList);
        private delegate int decl_cg_boco_info(int fn, int B, int Z, int BC, [Out] char[] boconame, out BCType_t bocotype, out PointSetType_t ptset_type, out int npnts, [Out] int[] NormalIndex, out int NormalListFlag, out DataType_t NormalDataType, out int ndataset);
        private delegate int decl_cg_ndescriptors(out int ndescriptors);
        private delegate int decl_cg_descriptor_write(char[] descr_name, char[] descr_text);
        private delegate int decl_cg_descriptor_read(int D, [Out] char[] descr_name, [Out] char[,] descr_text);
        private delegate int decl_cg_simulation_type_write(int file_number, int B, SimulationType_t type);
        private delegate int decl_cg_simulation_type_read(int file_number, int B, out SimulationType_t type);
        private delegate int decl_cg_biter_write(int file_number, int B, char[] bitername, int nsteps);
        private delegate int decl_cg_biter_read(int file_number, int B, [Out] char[] bitername, out int nsteps);
        private delegate int decl_cg_ziter_write(int file_number, int B, int Z, char[] zitername);
        private delegate int decl_cg_ziter_read(int file_number, int B, int Z, [Out] char[] zitername);
        private delegate int decl_cg_rind_write(int[] rind_data);
        private delegate int decl_cg_rind_read([Out] int[] rind_data);
        private delegate String decl_cg_get_error();

        decl_cg_open cg_open;
        decl_cg_close cg_close;
        decl_cg_golist cg_golist;
        decl_cg_gopath cg_gopath;
        decl_cg_where cg_where;
        decl_cg_nbases cg_nbases;
        decl_cg_base_write cg_base_write;
        decl_cg_base_read cg_base_read;
        decl_cg_nzones cg_nzones;
        decl_cg_zone_type cg_zone_type;
        decl_cg_zone_write cg_zone_write;
        decl_cg_zone_read cg_zone_read;
        decl_cg_nsections cg_nsections;
        decl_cg_section_write2 cg_section_write2;
        decl_cg_section_write1 cg_section_write1;
        decl_cg_section_read cg_section_read;
        decl_cg_elements_read2 cg_elements_read2;
        decl_cg_elements_read1 cg_elements_read1;
        decl_cg_gridlocation_write cg_gridlocation_write;
        decl_cg_gridlocation_read cg_gridlocation_read;
        decl_cg_ngrids cg_ngrids;
        decl_cg_ncoords cg_ncoords;
        decl_cg_coord_write3 cg_coord_write3;
        decl_cg_coord_write2 cg_coord_write2;
        decl_cg_coord_write1 cg_coord_write1;
        decl_cg_coord_read3 cg_coord_read3;
        decl_cg_coord_read2 cg_coord_read2;
        decl_cg_coord_read1 cg_coord_read1;
        decl_cg_coord_info cg_coord_info;
        decl_cg_nsols cg_nsols;
        decl_cg_sol_info cg_sol_info;
        decl_cg_sol_write cg_sol_write;
        decl_cg_nfields cg_nfields;
        decl_cg_field_info cg_field_info;
        decl_cg_field_write3 cg_field_write3;
        decl_cg_field_write2 cg_field_write2;
        decl_cg_field_write1 cg_field_write1;
        decl_cg_field_read3 cg_field_read3;
        decl_cg_field_read2 cg_field_read2;
        decl_cg_field_read1 cg_field_read1;
        decl_cg_nbocos cg_nbocos;
        decl_cg_boco_write3 cg_boco_write3;
        decl_cg_boco_write2 cg_boco_write2;
        decl_cg_boco_write1 cg_boco_write1;
        decl_cg_boco_read3 cg_boco_read3;
        decl_cg_boco_read2 cg_boco_read2;
        decl_cg_boco_read1 cg_boco_read1;
        decl_cg_boco_info cg_boco_info;
        decl_cg_ndescriptors cg_ndescriptors;
        decl_cg_descriptor_write cg_descriptor_write;
        decl_cg_descriptor_read cg_descriptor_read;
        decl_cg_simulation_type_write cg_simulation_type_write;
        decl_cg_simulation_type_read cg_simulation_type_read;
        decl_cg_biter_write cg_biter_write;
        decl_cg_biter_read cg_biter_read;
        decl_cg_ziter_write cg_ziter_write;
        decl_cg_ziter_read cg_ziter_read;
        decl_cg_rind_write cg_rind_write;
        decl_cg_rind_read cg_rind_read;
        decl_cg_get_error cg_get_error;

        private double[][, ,] x_of_all_zones;
        private double[][, ,] y_of_all_zones;
        private double[][, ,] z_of_all_zones;
        private int[][][][] bcelements_of_all_zones;
        private int[][][][] normal_indices_of_all_zones;
        private int[][] starts_of_all_zones;
        private int[][] ends_of_all_zones;
        private int[][] nbndrys_of_all_zones;
        private PointSetType_t[][][] pointset_types_of_all_zones;
        private GridLocation_t[][][] gridlocation_types_of_all_zones;
        private BCType_t[][][] bc_types_of_all_zones;
        private ElementType_t[][] element_types_of_all_zones;
        private int[][][,] elements_of_all_zones;
        private char[][][][] boconames_of_all_zones;
        private ZoneType_t[] zone_types_of_all_zones;

        private double[] x;
        private double[] y;
        private double[] z;
        private int[][] bcelements; // Either elements or points!
        private int[][] normal_indices;
        private PointSetType_t[] pointset_types;
        private GridLocation_t[] gridlocation_types;
        private BCType_t[] bc_types;
        private ElementType_t[] element_types;
        private int[][] elements;
        private int[][][] faces;
        private int[,] bcfaces;
        private char[][] boconames;
        private ZoneType_t[] zone_types;
        private int[][] to_node_index;
        private int[][][] to_element_index;

        private Boolean need_struct_conversion;
        private int number_of_coord_dimensions = 2;
        private List<int>[] node_to_element_indices;
        private List<int>[] bc_to_element_indices;
        private int bosss_elements;
        private int number_of_face_edges;
        private int number_of_faces_per_element;
        private int number_of_faces_per_element_mask;
        private int[] cgns_element_index_to_bosss_element_index;
        private int[] bosss_element_index_to_cgns_element_index;
        private ElementType_t bosss_element_type;
        private long[,] cell_neighbours;

        private bool hdf5;
        private GridCommons grd;
        private IEnumerable<BoSSS.Foundation.Field> fields;
        private double time;
        private string title;

        private bool Is64Bit() {
            if (IntPtr.Size == 8 || (IntPtr.Size == 4 && Is32BitProcessOn64BitProcessor())) {
                return true;
            } else {
                return false;
            }
        }

        private bool Is32BitProcessOn64BitProcessor() {
            bool retVal;

            IsWow64Process(Process.GetCurrentProcess().Handle, out retVal);

            return retVal;
        }

        /// <summary>
        /// Implement this methods by exporting/plotting all fields contained
        /// in <paramref name="fieldsToPlot"/> while regarding the settings
        /// <see cref="PlotDriver.superSampling"/>, <see cref="PlotDriver.ghostZone"/> and <see cref="PlotDriver.showJumps"/>.
        /// </summary>
        /// <param name="fileNameBase">
        /// The base name of the file (if any) that should be created. If
        /// called in parallel, rank information will be appended. If a file
        /// with the given name already exists, the resulting file name will
        /// contain an additional counter.
        /// </param>
        /// <param name="title">The title of the plot</param>
        /// <param name="time">
        /// The solution time represented by <paramref name="fieldsToPlot"/>
        /// </param>
        /// <param name="fieldsToPlot">
        /// A list of fields which should be plotted
        /// </param>
        public override void PlotFields(string fileNameBase, string title, double time, IEnumerable<BoSSS.Foundation.Field> fieldsToPlot) {
            string filename = GenerateFileName(fileNameBase);
            this.time = time;
            this.fields = fieldsToPlot;
            this.title = title;
            Build(filename, hdf5, true, context);
        }
        /// <summary>
        /// Implement this method by returning the file ending that should be
        /// appended at the end of a file name.
        /// </summary>
        protected override string FileEnding {
            get {
                if (hdf5) {
                    return "hdf5";
                } else {
                    return "cgns";
                }
            }
        }

        /// <summary>
        /// Permutation tables for the conversion between Tecplot ordering and
        /// BoSSS ordering of the nodes in an element.
        /// </summary>
        private Dictionary<ElementType_t, int[]> elementTypeToPermutationTableMap = new Dictionary<ElementType_t, int[]> {
            {ElementType_t.BAR_2, new int[] {0, 1}},
            {ElementType_t.TRI_3, new int[] {0, 1, 2}},
            {ElementType_t.QUAD_4, new int[] {0, 1, 3, 2}},
            {ElementType_t.TETRA_4, new int[] {0, 1, 2, 3}},
            {ElementType_t.HEXA_8, new int[] {0, 1, 4, 2, 3, 6, 5, 7}}
        };

        private void SwitchNode(ref int source, ref int target) {
            int buffer = source;

            source = target;
            target = buffer;
        }

        private void SwitchXZ() {
            double buffer;

            for (int i = 0; i < x.Length; i++) {
                buffer = x[i];
                x[i] = z[i];
                z[i] = buffer;
            }
        }

        private void SwitchYZ() {
            double buffer;

            for (int i = 0; i < x.Length; i++) {
                buffer = y[i];
                y[i] = z[i];
                z[i] = buffer;
            }
        }
        
        private void SwitchXY() {
            double buffer;

            for (int i = 0; i < x.Length; i++) {
                buffer = x[i];
                x[i] = y[i];
                y[i] = buffer;
            }
        }

        private void ThrowError() {
            throw new ApplicationException(cg_get_error());
        }

        private static void ThrowError(String s) {
            throw new ApplicationException(s);
        }

        private void ThrowError(int index_file) {
            String s = cg_get_error();

            if (cg_close(index_file) != (int)Error.CG_OK) {
                s = s + " " + cg_get_error();
            }

            throw new ApplicationException(s);
        }

        private void ThrowError(String s, int index_file) {
            if (cg_close(index_file) != (int)Error.CG_OK) {
                s = s + " " + cg_get_error();
            }

            throw new ApplicationException(s);
        }

        private int NumberOfFaceEdges(ElementType_t type) {
            if (Is3D(type)) {
                if (type.Equals(ElementType_t.HEXA_8)) {
                    return 4;
                }
                return 3;
            } else {
                if (number_of_coord_dimensions == 3) {
                    if (type.Equals(ElementType_t.TRI_3)) {
                        return 3;
                    }
                    if (type.Equals(ElementType_t.QUAD_4)) {
                        return 4;
                    }
                }
                return 2;
            }
        }

        private int NumberOfFacesPerElement(ElementType_t type) {
            if (number_of_coord_dimensions == 3) {
                if (type.Equals(ElementType_t.TRI_3) || type.Equals(ElementType_t.QUAD_4)) {
                    return 1;
                }
            }
            if (type.Equals(ElementType_t.HEXA_8)) {
                return 6;
            }
            char[] number = type.ToString().ToCharArray();
            return number[number.Length - 1] - 48;
        }

        private static int NumberOfNodesPerElement(ElementType_t type) {
            char[] number = type.ToString().ToCharArray();
            return number[number.Length - 1] - 48;
        }

        private static bool Is3D(ElementType_t type) {
            return type.Equals(ElementType_t.TETRA_4) || type.Equals(ElementType_t.HEXA_8);
        }

        private static bool BoSSSCompliantElementType(ElementType_t type) {
            return type.Equals(ElementType_t.TRI_3) || type.Equals(ElementType_t.TETRA_4) || type.Equals(ElementType_t.QUAD_4) || type.Equals(ElementType_t.HEXA_8);
        }

        private bool ChoosenBoSSSCompliantElementType(ElementType_t type) {
            return type.Equals(bosss_element_type);
        }

        private static String CleanString(string strIn) {
            return Regex.Replace(strIn, @"[^a-zA-Z0-9]", string.Empty);
        }

        /// <summary>
        /// Convert a structured grid to an unstructured grid
        /// </summary>
        private void ConvertStructure2Unstructure() {
            // insert here
            need_struct_conversion = false;
        }

        /// <summary>
        /// Generates zoneless arrays
        /// </summary>
        private void GenerateZonelessAndSectionlessArrays() {
            if (need_struct_conversion) {
                ThrowError("Struct conversion needed");
            }

            int i, j, k, l, m, n, a, b, c;

            zone_types = zone_types_of_all_zones;

            int array_length = 0;

            for (i = 0; i < zone_types.Length; i++) {
                array_length += x_of_all_zones[i].GetLength(0);
            }

            x = new double[array_length];
            y = new double[array_length];
            z = new double[array_length];

            node_to_element_indices = new List<int>[array_length];
            to_node_index = new int[zone_types.Length][];
            to_element_index = new int[zone_types.Length][][];
            for (i = 0, j = 0; i < zone_types.Length; i++) // over zone indices
            {
                a = x_of_all_zones[i].GetLength(0);
                to_node_index[i] = new int[a];
                for (k = 0; k < a; j++, k++) // over coordinate indices per zone
                {
                    x[j] = x_of_all_zones[i][k, 0, 0];
                    y[j] = y_of_all_zones[i][k, 0, 0];
                    z[j] = z_of_all_zones[i][k, 0, 0];
                    to_node_index[i][k] = j; // i|k index -> j index
                    node_to_element_indices[j] = new List<int>();
                }
            }

            array_length = 0;
            for (i = 0; i < elements_of_all_zones.Length; i++) {
                a = elements_of_all_zones[i].Length;
                for (k = 0; k < a; k++) {
                    array_length += elements_of_all_zones[i][k].GetLength(0);
                }
            }

            element_types = new ElementType_t[array_length];
            elements = new int[array_length][];
            for (i = 0, j = 0, m = 0; i < elements_of_all_zones.Length; i++) // over zone indices
            {
                a = elements_of_all_zones[i].Length;
                to_element_index[i] = new int[a][];
                for (k = 0; k < a; j++, k++) // over section indices per zone
                {
                    b = elements_of_all_zones[i][k].GetLength(0);
                    c = NumberOfNodesPerElement(element_types_of_all_zones[i][k]);
                    to_element_index[i][k] = new int[b];
                    for (l = 0; l < b; l++, m++) // over element indices (per section) per zone
                    {
                        to_element_index[i][k][l] = m; // i|k|l index -> m index
                        elements[m] = new int[c];
                        element_types[m] = element_types_of_all_zones[i][k];
                        if (Is3D(element_types[m])) {
                            number_of_coord_dimensions = 3;
                        }
                        if (BoSSSCompliantElementType(element_types[m])) {
                            if (verticesPerCell < c) {
                                verticesPerCell = c;
                                bosss_element_type = element_types[m];
                            }
                            if (number_of_face_edges < NumberOfFaceEdges(element_types[m])) {
                                number_of_face_edges = NumberOfFaceEdges(element_types[m]);
                                bosss_element_type = element_types[m];
                            }
                            if (number_of_faces_per_element < NumberOfFacesPerElement(element_types[m])) {
                                number_of_faces_per_element = NumberOfFacesPerElement(element_types[m]);
                                number_of_faces_per_element_mask = (1 << number_of_faces_per_element) - 1;
                                bosss_element_type = element_types[m];
                            }
                        }
                        for (n = 0; n < c; n++) // over coordinate indices per elements (per section) per zone
                        {
                            node_to_element_indices[elements[m][n] = to_node_index[i][elements_of_all_zones[i][k][l, n] - 1]].Add(m);
                        }
                    }
                }
            }
            for (m = 0; m < elements.Length; m++) {
                if (BoSSSCompliantElementType(element_types[m]) && (number_of_coord_dimensions == 2 || Is3D(element_types[m]))) {
                    c = NumberOfNodesPerElement(element_types[m]);
                    if (verticesPerCell > c) {
                        verticesPerCell = c;
                        bosss_element_type = element_types[m];
                    }
                    if (number_of_face_edges > NumberOfFaceEdges(element_types[m])) {
                        number_of_face_edges = NumberOfFaceEdges(element_types[m]);
                        bosss_element_type = element_types[m];
                    }
                    if (number_of_faces_per_element > NumberOfFacesPerElement(element_types[m])) {
                        number_of_faces_per_element = NumberOfFacesPerElement(element_types[m]);
                        number_of_faces_per_element_mask = (1 << number_of_faces_per_element) - 1;
                        bosss_element_type = element_types[m];
                    }
                }
            }
            cgns_element_index_to_bosss_element_index = new int[elements.Length];
            bosss_element_index_to_cgns_element_index = new int[elements.Length];
            for (m = 0; m < elements.Length; m++) {
                if (ChoosenBoSSSCompliantElementType(element_types[m])) {
                    bosss_element_index_to_cgns_element_index[cgns_element_index_to_bosss_element_index[m] = bosss_elements++] = m;
                }
            }
            array_length = 0;
            for (i = 0; i < bcelements_of_all_zones.Length; i++) // over zone indices
            {
                a = bcelements_of_all_zones[i].Length;
                for (k = 0; k < a; k++) // over section indices per zone
                {
                    array_length += bcelements_of_all_zones[i][k].Length;
                }
            }

            gridlocation_types = new GridLocation_t[array_length];
            pointset_types = new PointSetType_t[array_length];
            bc_types = new BCType_t[array_length];
            boconames = new char[array_length][];
            bcelements = new int[array_length][];
            normal_indices = new int[array_length][];
            for (i = 0, j = 0, m = 0; i < bcelements_of_all_zones.Length; i++) // over zone indices
            {
                a = bcelements_of_all_zones[i].Length;
                for (k = 0; k < a; j++, k++) // over section indices per zone
                {
                    b = bcelements_of_all_zones[i][k].Length;
                    for (l = 0; l < b; l++, m++) // over boundary condition indices (per section) per zone
                    {
                        gridlocation_types[m] = gridlocation_types_of_all_zones[i][k][l];

                        if (!gridlocation_types[m].Equals(GridLocation_t.FaceCenter) && !gridlocation_types[m].Equals(GridLocation_t.Vertex)) {
                            ThrowError("Unsupported gridlocation");
                        }

                        pointset_types[m] = pointset_types_of_all_zones[i][k][l];
                        bc_types[m] = bc_types_of_all_zones[i][k][l];
                        boconames[m] = boconames_of_all_zones[i][k][l];
                        normal_indices[m] = normal_indices_of_all_zones[i][k][l];
                        c = bcelements_of_all_zones[i][k][l].Length;
                        bcelements[m] = new int[c];
                        for (n = 0; n < c; n++) // over coordinate indices per boundary condition (per section) per zone
                        {
                            if (gridlocation_types[m].Equals(GridLocation_t.Vertex)) {
                                bcelements[m][n] = to_node_index[i][bcelements_of_all_zones[i][k][l][n] - 1];
                            } else {
                                bcelements[m][n] = bcelements_of_all_zones[i][k][l][n] - 1;
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Reorder the coordinates inside an element
        /// </summary>
        private void VertexIndexNumberingTransformation() {
            for (int i = 0; i < elements.Length; i++) {
                if (element_types[i].Equals(ElementType_t.QUAD_4)) {
                    SwitchNode(ref elements[i][2], ref elements[i][3]);
                }
                if (!element_types[i].Equals(ElementType_t.HEXA_8)) {
                    continue;
                }
                SwitchNode(ref elements[i][2], ref elements[i][3]);
                SwitchNode(ref elements[i][3], ref elements[i][4]);
                SwitchNode(ref elements[i][5], ref elements[i][6]);
            }
        }

        /// <summary>
        /// Reorder the coordinates inside an element
        /// </summary>
        private void CoordinateSystemTransformation() {
            // insert here
        }

        /// <summary>
        /// Builds faces from elements
        /// </summary>
        private void BuildFaces() {
            int i, j;

            faces = new int[elements.Length][][];

            bc_to_element_indices = new List<int>[bc_types.Length];

            for (i = 0; i < bc_types.Length; i++) // over elements
            {
                bc_to_element_indices[i] = new List<int>();
            }
            for (i = 0; i < faces.Length; i++) // over elements
            {
                if (element_types[i].Equals(ElementType_t.BAR_2)) {
                    faces[i] = new int[1][];
                    faces[i][0] = new int[2];
                    faces[i][0][0] = elements[i][0];
                    faces[i][0][1] = elements[i][1];
                } else {
                    faces[i] = new int[NumberOfFacesPerElement(element_types[i])][];

                    for (j = 0; j < faces[i].Length; j++) {
                        faces[i][j] = new int[NumberOfFaceEdges(element_types[i])];
                    }
                    if (element_types[i].Equals(ElementType_t.TRI_3)) {
                        if (number_of_coord_dimensions == 2) {
                            faces[i][0][0] = elements[i][0];
                            faces[i][0][1] = elements[i][1];
                            faces[i][1][0] = elements[i][1];
                            faces[i][1][1] = elements[i][2];
                            faces[i][2][0] = elements[i][2];
                            faces[i][2][1] = elements[i][0];
                        } else {
                            faces[i][0][0] = elements[i][0];
                            faces[i][0][1] = elements[i][1];
                            faces[i][0][2] = elements[i][2];
                        }
                    } else if (element_types[i].Equals(ElementType_t.TETRA_4)) {
                        faces[i][0][0] = elements[i][1];
                        faces[i][0][1] = elements[i][2];
                        faces[i][0][2] = elements[i][3];
                        faces[i][1][0] = elements[i][0];
                        faces[i][1][1] = elements[i][2];
                        faces[i][1][2] = elements[i][3];
                        faces[i][2][0] = elements[i][0];
                        faces[i][2][1] = elements[i][1];
                        faces[i][2][2] = elements[i][2];
                        faces[i][3][0] = elements[i][0];
                        faces[i][3][1] = elements[i][1];
                        faces[i][3][2] = elements[i][3];
                    } else if (element_types[i].Equals(ElementType_t.QUAD_4)) {
                        if (number_of_coord_dimensions == 2) {
                            faces[i][0][0] = elements[i][0];
                            faces[i][0][1] = elements[i][2];
                            faces[i][1][0] = elements[i][1];
                            faces[i][1][1] = elements[i][3];
                            faces[i][2][0] = elements[i][2];
                            faces[i][2][1] = elements[i][3];
                            faces[i][3][0] = elements[i][0];
                            faces[i][3][1] = elements[i][1];
                        } else {
                            faces[i][0][0] = elements[i][0];
                            faces[i][0][1] = elements[i][1];
                            faces[i][0][2] = elements[i][2];
                            faces[i][0][3] = elements[i][3];
                        }
                    } else if (element_types[i].Equals(ElementType_t.HEXA_8)) {
                        faces[i][0][0] = elements[i][0];
                        faces[i][0][1] = elements[i][3];
                        faces[i][0][2] = elements[i][7];
                        faces[i][0][3] = elements[i][2];
                        faces[i][1][0] = elements[i][1];
                        faces[i][1][1] = elements[i][4];
                        faces[i][1][2] = elements[i][5];
                        faces[i][1][3] = elements[i][6];
                        faces[i][2][0] = elements[i][2];
                        faces[i][2][1] = elements[i][4];
                        faces[i][2][2] = elements[i][5];
                        faces[i][2][3] = elements[i][7];
                        faces[i][3][0] = elements[i][0];
                        faces[i][3][1] = elements[i][1];
                        faces[i][3][2] = elements[i][6];
                        faces[i][3][3] = elements[i][3];
                        faces[i][4][0] = elements[i][3];
                        faces[i][4][1] = elements[i][7];
                        faces[i][4][2] = elements[i][5];
                        faces[i][4][3] = elements[i][6];
                        faces[i][5][0] = elements[i][0];
                        faces[i][5][1] = elements[i][2];
                        faces[i][5][2] = elements[i][4];
                        faces[i][5][3] = elements[i][1];
                    }
                }
            }
        }

        /// <summary>
        /// Builds and shares the boundary conditions of one point with all other points positioned at the same location
        /// </summary>
        private void BuildAndShareBoundaryConditions() {
            int i, j, k, l, m, n;
            int range_from;
            int range_to;
            int index;
            int index2;
            int edge_counter;
            int verify;

            bcfaces = new int[elements.Length, number_of_faces_per_element];
            for (i = 0; i < bcfaces.GetLength(0); i++) {
                for (j = 0; j < bcfaces.GetLength(1); j++) {
                    bcfaces[i, j] = -1;
                }
            }
            for (i = 0; i < bcelements.Length; i++) // over boundary condition indices
            {
                range_from = 0;
                range_to = bcelements[i].Length;

                if (pointset_types[i].Equals(PointSetType_t.ElementRange) || pointset_types[i].Equals(PointSetType_t.PointRange)) {
                    range_from = bcelements[i][0];
                    range_to = bcelements[i][1] + 1;
                }
                for (j = range_from; j < range_to; j++) // over coordinate indices per boundary condition (if branch) or element indices per boundary condition (else branch)
                {
                    if (pointset_types[i].Equals(PointSetType_t.ElementRange) || pointset_types[i].Equals(PointSetType_t.PointRange)) {
                        index = j;
                    } else {
                        index = bcelements[i][j];
                    }
                    if (gridlocation_types[i].Equals(GridLocation_t.Vertex)) {
                        foreach (int z in node_to_element_indices[index]) // over elements
                        {
                            for (k = 0; k < faces[z].Length; k++) // over face indices per element
                            {
                                if (k >= bcfaces.GetLength(1)) {
                                    continue;
                                }
                                for (m = edge_counter = 0; m < faces[z][k].Length; m++) // over coordinate indices per faces
                                {
                                    for (n = range_from; n < range_to; n++) // over coordinate indices per boundary condition
                                    {
                                        if (pointset_types[i].Equals(PointSetType_t.ElementRange) || pointset_types[i].Equals(PointSetType_t.PointRange)) {
                                            index2 = n;
                                        } else {
                                            index2 = bcelements[i][n];
                                        }
                                        if (index2 == faces[z][k][m]) {
                                            edge_counter++;
                                        }
                                    }
                                }
                                if (edge_counter == number_of_face_edges) {
                                    bcfaces[z, k] = i;
                                }
                            }
                        }
                    } else {
                        for (k = 0; k < faces[j].Length; k++) // over face indices per element
                        {
                            if (k >= bcfaces.GetLength(1)) {
                                continue;
                            }
                            bcfaces[j, k] = i;
                        }
                    }
                }
                for (j = 0; j < faces.Length; j++) // over element indices
                {
                    for (k = 0; k < faces[j].Length; k++) // over face indices per element
                    {
                        foreach (int z in node_to_element_indices[faces[j][k][0]]) // over elements
                        {
                            if (z == j || !ChoosenBoSSSCompliantElementType(element_types[z])) {
                                continue;
                            }
                            for (l = 0; l < faces[z].Length; l++) // over face indices per element
                            {
                                if (l >= bcfaces.GetLength(1)) {
                                    continue;
                                }
                                for (m = edge_counter = 0; m < faces[z][l].Length; m++) // over coordinate indices per faces
                                {
                                    for (n = 0; n < faces[j][k].Length; n++) // over coordinate indices per faces
                                    {
                                        if (faces[j][k][n] == faces[z][l][m]) {
                                            edge_counter++;
                                        }
                                    }
                                }
                                if (edge_counter == number_of_face_edges) {
                                    if (bcfaces[z, l] == -1) {
                                        bcfaces[z, l] = bcfaces[j, k];
                                    }
                                }
                            }
                        }
                    }
                }
            }
            if (bc_to_element_indices.Length == 0) {
                return;
            }
            for (i = 0; i < faces.Length; i++) // over elements
            {
                if (!ChoosenBoSSSCompliantElementType(element_types[i])) {
                    continue;
                }
                for (j = 0; j < faces[i].Length; j++) // over faces
                {
                    if ((verify = bcfaces[i, j]) >= 0) {
                        bc_to_element_indices[verify].Add(i << number_of_faces_per_element | j);
                    }
                }
            }
        }

        /// <summary>
        /// Build neighbour relations
        /// </summary>
        private void BuildNeighbourship() {
            int i, j, k, l, m, n;
            int edge_counter;

            cell_neighbours = new long[bosss_elements, number_of_faces_per_element];

            for (i = 0; i < bosss_elements; i++) // over elements
            {
                for (j = 0; j < number_of_faces_per_element; j++) // over face indices per element
                {
                    cell_neighbours[i, j] = -1;
                }
            }
            for (n = 0; n < bosss_elements; n++) // over elements
            {
                i = bosss_element_index_to_cgns_element_index[n];
                if (!ChoosenBoSSSCompliantElementType(element_types[i])) {
                    continue;
                }
                for (j = 0; j < faces[i].Length; j++) // over face indices per element
                {
                    if (cell_neighbours[cgns_element_index_to_bosss_element_index[i], j] > -1) {
                        continue;
                    }
                    foreach (int z in node_to_element_indices[faces[i][j][0]]) // over elements
                    {
                        if (z == i || !ChoosenBoSSSCompliantElementType(element_types[z])) {
                            continue;
                        }
                        for (l = 0; l < faces[z].Length; l++) // over face indices per element
                        {
                            for (m = edge_counter = 0; m < faces[z][l].Length; m++) // over coordinate indices per faces
                            {
                                for (k = 0; k < faces[i][j].Length; k++) // over coordinate indices per faces
                                {
                                    if (faces[i][j][k] == faces[z][l][m]) {
                                        edge_counter++;
                                    }
                                }
                            }
                            if (edge_counter == number_of_face_edges) {
                                cell_neighbours[cgns_element_index_to_bosss_element_index[i], j] = cgns_element_index_to_bosss_element_index[z];
                            }
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Generates a grid commons object
        /// </summary>
        /// <returns>the grid commons object</returns>
        public GridCommons GenerateBoSSSGrid() {
            int i, j, k, l;
            int max_bc_index = -1;

            GenerateZonelessAndSectionlessArrays();
            CoordinateSystemTransformation();
            VertexIndexNumberingTransformation();
            BuildFaces();
            BuildAndShareBoundaryConditions();
            BuildNeighbourship();

            double[, ,] vertices = new double[bosss_elements, verticesPerCell, number_of_coord_dimensions];
            long[] globalID = new long[bosss_elements];
            byte[,] edgeTags = new byte[bosss_elements, number_of_faces_per_element];

            for (l = 0, k = 0; l < bosss_elements; l++) {
                j = bosss_element_index_to_cgns_element_index[l];
                for (int n = 0; n < verticesPerCell; n++) {
                    k = cgns_element_index_to_bosss_element_index[j];

                    int node_no = elements[j][n];

                    vertices[k, n, 0] = x[node_no];
                    vertices[k, n, 1] = y[node_no];
                    if (number_of_coord_dimensions == 3) {
                        vertices[k, n, 2] = z[node_no];
                    }
                }
            }

            for (j = 0; j < bosss_elements; j++) {
                globalID[j] = j;
            }

            GridCommons grd = null;
            if (bosss_element_type.Equals(ElementType_t.QUAD_4)) {
                grd = new Cartesian2DGrid();
            } else if (bosss_element_type.Equals(ElementType_t.TRI_3)) {
                grd = new UnstructuredTriangleGrid();
            } else if (bosss_element_type.Equals(ElementType_t.HEXA_8)) {
                grd = new Cartesian3DGrid();
            } else if (bosss_element_type.Equals(ElementType_t.TETRA_4)) {
                grd = new UnstructuredTetraGrid();
            } else {
                ThrowError("Unsupported BoSSS element type");
            }

            int tag_number = 1;
            for (i = 0; i < bcelements.Length; i++) // over boundary condition indices
            {
                foreach (int z in bc_to_element_indices[i]) // over elements
                {
                    j = z >> number_of_faces_per_element;
                    j = cgns_element_index_to_bosss_element_index[j];
                    k = z & number_of_faces_per_element_mask;
                    edgeTags[j, k] = (byte)tag_number;
                    if (tag_number > max_bc_index) {
                        max_bc_index = tag_number;
                    }
                }
                tag_number++;
            }

            tag_number = 1;
            for (i = 0; i < max_bc_index; i++) // over boundary condition indices
            {
                grd.EdgeTagsNames.Add((byte)tag_number, CleanString(new string(boconames[i])));
                tag_number++;
            }

            grd.Description = "Created from BoSSS";
            grd.GlobalID = globalID;
            grd.CellNeighbours = cell_neighbours;
            grd.Vertices = vertices;
            grd.EdgeTags = edgeTags;

            return grd;
        }

        /// <summary>
        /// Constructs a new PlotDriver
        /// </summary>
        /// <param name="context">The omnipresent context</param>
        /// <param name="hdf5"></param>
        /// <param name="showJumps">
        /// <param name="ghostZone"><see cref="PlotDriver"/></param>
        /// If true, the real DG data (including discontinuities) should be
        /// exported. Otherwise, the mean value of all values in one node
        /// should be calculated
        /// </param>
        /// <param name="superSampling"><see cref="PlotDriver.superSampling"/></param>
        public Cgns(Context context, bool hdf5, bool showJumps, bool ghostZone, uint superSampling)
            : base(context, showJumps, ghostZone, superSampling) {
            this.hdf5 = hdf5;
        }

        /// <summary>
        /// loads or saves an adf or hdf5 file from or to a specified path
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="hdf5"></param>
        /// <param name="save"></param>
        /// <param name="context"></param>
        public Cgns(string filePath, bool hdf5, bool save, Context context)
            : base(null, false, false, 0) {
            Build(filePath, hdf5, save, context);
        }

        /// <summary>
        /// loads or saves an adf or hdf5 file from or to a specified path
        /// </summary>
        /// <param name="filePath"></param>
        /// <param name="hdf5"></param>
        /// <param name="save"></param>
        /// <param name="context"></param>
        public void Build(string filePath, bool hdf5, bool save, Context context) {
            int index_file = 0;
            int index_base = 0;
            int index_zone = 0;
            int index_section = 0;
            int index_bc = 0;
            int nbases;
            int nzones;
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
            int celldim, physdim; // zone definition
            int start, end, nbndry, parent_flag; // element definition
            int normal_list_flag, nbcelement, ndataset; // bc definition
            DataType_t data_type;
            ElementType_t element_type;
            ZoneType_t zone_type;
            BCType_t bc_type;
            GridLocation_t gridlocation_type;
            PointSetType_t pointset_type;

            if (Is32BitProcessOn64BitProcessor()) {
                Console.WriteLine("Cgns runs as 32bit process ...");
            } else {
                Console.WriteLine("Cgns runs as 64bit process ...");
            }

            if (hdf5) {
                cg_open = Hdf5_cg_open;
                cg_close = Hdf5_cg_close;
                cg_gopath = Hdf5_cg_gopath;
                cg_golist = Hdf5_cg_golist;
                cg_where = Hdf5_cg_where;
                cg_nbases = Hdf5_cg_nbases;
                cg_base_write = Hdf5_cg_base_write;
                cg_base_read = Hdf5_cg_base_read;
                cg_nzones = Hdf5_cg_nzones;
                cg_zone_type = Hdf5_cg_zone_type;
                cg_zone_write = Hdf5_cg_zone_write;
                cg_zone_read = Hdf5_cg_zone_read;
                cg_nsections = Hdf5_cg_nsections;
                cg_section_write2 = Hdf5_cg_section_write2;
                cg_section_write1 = Hdf5_cg_section_write1;
                cg_section_read = Hdf5_cg_section_read;
                cg_elements_read2 = Hdf5_cg_elements_read2;
                cg_elements_read1 = Hdf5_cg_elements_read1;
                cg_gridlocation_write = Hdf5_cg_gridlocation_write;
                cg_gridlocation_read = Hdf5_cg_gridlocation_read;
                cg_ngrids = Hdf5_cg_ngrids;
                cg_ncoords = Hdf5_cg_ncoords;
                cg_coord_write3 = Hdf5_cg_coord_write3;
                cg_coord_write2 = Hdf5_cg_coord_write2;
                cg_coord_write1 = Hdf5_cg_coord_write1;
                cg_coord_read3 = Hdf5_cg_coord_read3;
                cg_coord_read2 = Hdf5_cg_coord_read2;
                cg_coord_read1 = Hdf5_cg_coord_read1;
                cg_coord_info = Hdf5_cg_coord_info;
                cg_nsols = Hdf5_cg_nsols;
                cg_sol_info = Hdf5_cg_sol_info;
                cg_sol_write = Hdf5_cg_sol_write;
                cg_nfields = Hdf5_cg_nfields;
                cg_field_info = Hdf5_cg_field_info;
                cg_field_write3 = Hdf5_cg_field_write3;
                cg_field_write2 = Hdf5_cg_field_write2;
                cg_field_write1 = Hdf5_cg_field_write1;
                cg_field_read3 = Hdf5_cg_field_read3;
                cg_field_read2 = Hdf5_cg_field_read2;
                cg_field_read1 = Hdf5_cg_field_read1;
                cg_nbocos = Hdf5_cg_nbocos;
                cg_boco_write3 = Hdf5_cg_boco_write3;
                cg_boco_write2 = Hdf5_cg_boco_write2;
                cg_boco_write1 = Hdf5_cg_boco_write1;
                cg_boco_read3 = Hdf5_cg_boco_read3;
                cg_boco_read2 = Hdf5_cg_boco_read2;
                cg_boco_read1 = Hdf5_cg_boco_read1;
                cg_boco_info = Hdf5_cg_boco_info;
                cg_ndescriptors = Hdf5_cg_ndescriptors;
                cg_descriptor_write = Hdf5_cg_descriptor_write;
                cg_descriptor_read = Hdf5_cg_descriptor_read;
                cg_simulation_type_write = Hdf5_cg_simulation_type_write;
                cg_simulation_type_read = Hdf5_cg_simulation_type_read;
                cg_biter_write = Hdf5_cg_biter_write;
                cg_biter_read = Hdf5_cg_biter_read;
                cg_ziter_write = Hdf5_cg_ziter_write;
                cg_ziter_read = Hdf5_cg_ziter_read;
                cg_rind_write = Hdf5_cg_rind_write;
                cg_rind_read = Hdf5_cg_rind_read;
                cg_get_error = Hdf5_cg_get_error;
            } else {
                cg_open = Adf_cg_open;
                cg_close = Adf_cg_close;
                cg_gopath = Adf_cg_gopath;
                cg_golist = Adf_cg_golist;
                cg_where = Adf_cg_where;
                cg_nbases = Adf_cg_nbases;
                cg_base_write = Adf_cg_base_write;
                cg_base_read = Adf_cg_base_read;
                cg_nzones = Adf_cg_nzones;
                cg_zone_type = Adf_cg_zone_type;
                cg_zone_write = Adf_cg_zone_write;
                cg_zone_read = Adf_cg_zone_read;
                cg_nsections = Adf_cg_nsections;
                cg_section_write2 = Adf_cg_section_write2;
                cg_section_write1 = Adf_cg_section_write1;
                cg_section_read = Adf_cg_section_read;
                cg_elements_read2 = Adf_cg_elements_read2;
                cg_elements_read1 = Adf_cg_elements_read1;
                cg_gridlocation_write = Adf_cg_gridlocation_write;
                cg_gridlocation_read = Adf_cg_gridlocation_read;
                cg_ngrids = Adf_cg_ngrids;
                cg_ncoords = Adf_cg_ncoords;
                cg_coord_write3 = Adf_cg_coord_write3;
                cg_coord_write2 = Adf_cg_coord_write2;
                cg_coord_write1 = Adf_cg_coord_write1;
                cg_coord_read3 = Adf_cg_coord_read3;
                cg_coord_read2 = Adf_cg_coord_read2;
                cg_coord_read1 = Adf_cg_coord_read1;
                cg_coord_info = Adf_cg_coord_info;
                cg_nsols = Adf_cg_nsols;
                cg_sol_info = Adf_cg_sol_info;
                cg_sol_write = Adf_cg_sol_write;
                cg_nfields = Adf_cg_nfields;
                cg_field_info = Adf_cg_field_info;
                cg_field_write3 = Adf_cg_field_write3;
                cg_field_write2 = Adf_cg_field_write2;
                cg_field_write1 = Adf_cg_field_write1;
                cg_field_read3 = Adf_cg_field_read3;
                cg_field_read2 = Adf_cg_field_read2;
                cg_field_read1 = Adf_cg_field_read1;
                cg_nbocos = Adf_cg_nbocos;
                cg_boco_write3 = Adf_cg_boco_write3;
                cg_boco_write2 = Adf_cg_boco_write2;
                cg_boco_write1 = Adf_cg_boco_write1;
                cg_boco_read3 = Adf_cg_boco_read3;
                cg_boco_read2 = Adf_cg_boco_read2;
                cg_boco_read1 = Adf_cg_boco_read1;
                cg_boco_info = Adf_cg_boco_info;
                cg_ndescriptors = Adf_cg_ndescriptors;
                cg_descriptor_write = Adf_cg_descriptor_write;
                cg_descriptor_read = Adf_cg_descriptor_read;
                cg_simulation_type_write = Adf_cg_simulation_type_write;
                cg_simulation_type_read = Adf_cg_simulation_type_read;
                cg_biter_write = Adf_cg_biter_write;
                cg_biter_read = Adf_cg_biter_read;
                cg_ziter_write = Adf_cg_ziter_write;
                cg_ziter_read = Adf_cg_ziter_read;
                cg_rind_write = Adf_cg_rind_write;
                cg_rind_read = Adf_cg_rind_read;
                cg_get_error = Adf_cg_get_error;
            }
            if (save) {
                basename = "Base".ToCharArray();
                zonename = "Zone".ToCharArray();
                coordnameX = "CoordinateX".ToCharArray();
                coordnameY = "CoordinateY".ToCharArray();
                coordnameZ = "CoordinateZ".ToCharArray();
                flow = "VectorMagnitude".ToCharArray();
                secname = "Section".ToCharArray();
                boconame = "BoundaryCondition".ToCharArray();
                solname = "FlowSolution".ToCharArray();
                data_type = DataType_t.RealDouble;

                int index_coord;
                int index_flow;
                int index_field;

                if (cg_open(filePath.ToCharArray(), (int)Mode.MODE_WRITE, out index_file) != (int)Error.CG_OK) {
                    ThrowError(index_file);
                }

                celldim = physdim = context.Grid.SpatialDimension;

                if (cg_base_write(index_file, basename, celldim, physdim, out index_base) != (int)Error.CG_OK) {
                    ThrowError(index_file);
                }

                cg_gopath(index_file, "/Base".ToCharArray());

                cg_descriptor_write("Information".ToCharArray(), title.ToCharArray());

                grd = context.Grid;

                if (grd.GridSimplex.GetType() == typeof(Line)) {
                    bosss_element_type = ElementType_t.BAR_2;
                } else if (grd.GridSimplex.GetType() == typeof(Square)) {
                    bosss_element_type = ElementType_t.QUAD_4;
                } else if (grd.GridSimplex.GetType() == typeof(Triangle)) {
                    bosss_element_type = ElementType_t.TRI_3;
                } else if (grd.GridSimplex.GetType() == typeof(Cube)) {
                    bosss_element_type = ElementType_t.HEXA_8;
                } else if (grd.GridSimplex.GetType() == typeof(Tetra)) {
                    bosss_element_type = ElementType_t.TETRA_4;
                } else {
                    ThrowError("Unsupported BoSSS element type", index_file);
                }
                
                isize[0, 0] = totalVertices;
                isize[0, 1] = totalCells;
                isize[0, 2] = 0;

                double[, ,] x = null;
                double[, ,] y = null;
                double[, ,] z = null;

                if (cg_zone_write(index_file, index_base, zonename, isize, ZoneType_t.Unstructured, out index_zone) != (int)Error.CG_OK) {
                    ThrowError(index_file);
                }

                x = new double[isize[0, 0], 1, 1];
                for (int vertex = 0; vertex < isize[0, 0]; vertex++) {
                    x[vertex, 0, 0] = vertices[vertex, 0];
                }
                if (cg_coord_write3(index_file, index_base, index_zone, data_type, coordnameX, x, out index_coord) != (int)Error.CG_OK) {
                    ThrowError(index_file);
                }
                y = new double[isize[0, 0], 1, 1];
                for (int vertex = 0; vertex < isize[0, 0]; vertex++) {
                    y[vertex, 0, 0] = vertices[vertex, 1];
                }
                if (cg_coord_write3(index_file, index_base, index_zone, data_type, coordnameY, y, out index_coord) != (int)Error.CG_OK) {
                    ThrowError(index_file);
                }
                if (celldim > 2) {
                    z = new double[isize[0, 0], 1, 1];
                    for (int vertex = 0; vertex < isize[0, 0]; vertex++) {
                        z[vertex, 0, 0] = vertices[vertex, 2];
                    }
                    if (cg_coord_write3(index_file, index_base, index_zone, data_type, coordnameZ, z, out index_coord) != (int)Error.CG_OK) {
                        ThrowError(index_file);
                    }
                }

                if (cg_sol_write(index_file, index_base, index_zone, solname, GridLocation_t.Vertex, out index_flow) != (int)Error.CG_OK) {
                    ThrowError(index_file);
                }

                foreach (Field field in fields) {
                    calculateField3(field);

                    if (!showJumps) {
                        cg_field_write3(index_file, index_base, index_zone, index_flow, data_type, field.ToString().ToCharArray(), smoothedResult3, out index_field);
                    } else {
                        cg_field_write3(index_file, index_base, index_zone, index_flow, data_type, field.ToString().ToCharArray(), notSmoothedResult3, out index_field);
                    }
                }

                int[,] permutedConnectivity = new int[isize[0, 1], connectivity.GetLength(1)];
                int[] permutationTable = elementTypeToPermutationTableMap[bosss_element_type];
                for (int i = 0; i < isize[0, 1]; i++) {
                    for (int j = 0; j < permutationTable.Length; j++) {
                        permutedConnectivity[i, j] = 1 + connectivity[i, permutationTable[j]];
                    }
                }

                start = 1;
                end = isize[0, 1];
                nbndry = 0;

                if (cg_section_write2(index_file, index_base, index_zone, secname, bosss_element_type, start, end, nbndry, permutedConnectivity, out index_section) != (int)Error.CG_OK) {
                    ThrowError(index_file);
                }
            } else {
                index_file = index_base = 1;

                if (cg_open(filePath.ToCharArray(), (int)Mode.MODE_READ, out index_file) != (int)Error.CG_OK) {
                    ThrowError(index_file);
                }
                if (cg_nbases(index_file, out nbases) != (int)Error.CG_OK) {
                    ThrowError(index_file);
                }
                if (nbases != 1) {
                    ThrowError("Only one base allowed", index_file);
                }
                if (cg_base_read(index_file, index_base, basename, out celldim, out physdim) != (int)Error.CG_OK) {
                    ThrowError(index_file);
                }
                if (cg_nzones(index_file, index_base, out nzones) != (int)Error.CG_OK) {
                    ThrowError(index_file);
                }

                x_of_all_zones = new double[nzones][, ,];
                y_of_all_zones = new double[nzones][, ,];
                z_of_all_zones = new double[nzones][, ,];
                bcelements_of_all_zones = new int[nzones][][][];
                normal_indices_of_all_zones = new int[nzones][][][];
                starts_of_all_zones = new int[nzones][];
                ends_of_all_zones = new int[nzones][];
                nbndrys_of_all_zones = new int[nzones][];
                pointset_types_of_all_zones = new PointSetType_t[nzones][][];
                gridlocation_types_of_all_zones = new GridLocation_t[nzones][][];
                bc_types_of_all_zones = new BCType_t[nzones][][];
                element_types_of_all_zones = new ElementType_t[nzones][];
                elements_of_all_zones = new int[nzones][][,];
                boconames_of_all_zones = new char[nzones][][][];
                zone_types_of_all_zones = new ZoneType_t[nzones];

                for (index_zone = 1; index_zone <= nzones; index_zone++) {
                    if (cg_zone_type(index_file, index_base, index_zone, out zone_type) != (int)Error.CG_OK) {
                        ThrowError(index_file);
                    }
                    zone_types_of_all_zones[index_zone - 1] = zone_type;
                    if (cg_zone_read(index_file, index_base, index_zone, zonename, isize) != (int)Error.CG_OK) {
                        ThrowError(index_file);
                    }
                    if (cg_ncoords(index_file, index_base, index_zone, out ncoords) != (int)Error.CG_OK) {
                        ThrowError(index_file);
                    }
                    if (cg_coord_info(index_file, index_base, index_zone, 1, out data_type, coordnameX) != (int)Error.CG_OK) {
                        ThrowError(index_file);
                    }
                    if (cg_coord_info(index_file, index_base, index_zone, 2, out data_type, coordnameY) != (int)Error.CG_OK) {
                        ThrowError(index_file);
                    }
                    if (cg_coord_info(index_file, index_base, index_zone, 3, out data_type, coordnameZ) != (int)Error.CG_OK) {
                        ThrowError(index_file);
                    }

                    irmin[0] = 1;
                    irmin[1] = 1;
                    irmin[2] = 1;
                    irmax[0] = isize[0, 0];
                    irmax[1] = isize[0, 1];
                    irmax[2] = isize[0, 2];

                    double[, ,] x = null;
                    double[, ,] y = null;
                    double[, ,] z = null;

                    if ((int)zone_type == (int)ZoneType_t.Unstructured) {
                        x = new double[isize[0, 0], 1, 1];
                        y = new double[isize[0, 0], 1, 1];
                        z = new double[isize[0, 0], 1, 1];
                    } else {
                        x = new double[isize[0, 2], isize[0, 1], isize[0, 0]];
                        y = new double[isize[0, 2], isize[0, 1], isize[0, 0]];
                        z = new double[isize[0, 2], isize[0, 1], isize[0, 0]];
                    }

                    char[][][] boconames = null;
                    char[][] boconametemp = null;
                    int[][][] bcelements = null;
                    int[][] bcelementtemp = null;
                    int[] bcelement = null;
                    int[][][] normal_indices = null;
                    int[][] normal_indextemp = null;
                    int[] normal_index = null;
                    int[] starts = null;
                    int[] ends = null;
                    int[] nbndrys = null;
                    PointSetType_t[][] pointset_types = null;
                    PointSetType_t[] pointset_typtemp = null;
                    GridLocation_t[][] gridlocation_types = null;
                    GridLocation_t[] gridlocation_typtemp = null;
                    BCType_t[][] bc_types = null;
                    BCType_t[] bc_typtemp = null;
                    ElementType_t[] element_types = null;
                    int[][,] elements = null;

                    if (cg_coord_read3(index_file, index_base, index_zone, coordnameX, data_type, irmin, irmax, x) != (int)Error.CG_OK) {
                        ThrowError(index_file);
                    }
                    if (cg_coord_read3(index_file, index_base, index_zone, coordnameY, data_type, irmin, irmax, y) != (int)Error.CG_OK) {
                        ThrowError(index_file);
                    }
                    if (cg_coord_read3(index_file, index_base, index_zone, coordnameZ, data_type, irmin, irmax, z) != (int)Error.CG_OK) {
                        ThrowError(index_file);
                    }
                    if (cg_nsections(index_file, index_base, index_zone, out nsections) != (int)Error.CG_OK) {
                        ThrowError(index_file);
                    }

                    if ((int)zone_type == (int)ZoneType_t.Unstructured) {
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
                        for (index_section = 1; index_section <= nsections; index_section++) {
                            if (cg_section_read(index_file, index_base, index_zone, index_section, secname, out element_type, out start, out end, out nbndry, out parent_flag) != (int)Error.CG_OK) {
                                ThrowError(index_file);
                            }

                            starts[index_section - 1] = start;
                            ends[index_section - 1] = end;
                            nbndrys[index_section - 1] = nbndry;
                            element_types[index_section - 1] = element_type;

                            char[] c = element_type.ToString().ToCharArray();

                            ielem = new int[end - start + 1, c[c.Length - 1] - 48];

                            if (cg_elements_read2(index_file, index_base, index_zone, index_section, ielem, null) != (int)Error.CG_OK) {
                                ThrowError(index_file);
                            }

                            elements[index_section - 1] = ielem;

                            if (cg_nbocos(index_file, index_base, index_zone, out nbocos) != (int)Error.CG_OK) {
                                ThrowError(index_file);
                            }
                            boconametemp = new char[nbocos][];
                            bcelementtemp = new int[nbocos][];
                            normal_indextemp = new int[nbocos][];
                            pointset_typtemp = new PointSetType_t[nbocos];
                            gridlocation_typtemp = new GridLocation_t[nbocos];
                            bc_typtemp = new BCType_t[nbocos];

                            int depht = 3;
                            char[] label0 = "Zone_t".ToCharArray(), label1 = "ZoneBC_t".ToCharArray(), label2 = "BC_t".ToCharArray();
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
                                if (cg_boco_info(index_file, index_base, index_zone, index_bc, boconame, out bc_type, out pointset_type, out nbcelement, normal_index, out normal_list_flag, out data_type, out ndataset) != (int)Error.CG_OK) {
                                    ThrowError(index_file);
                                }
                                bcelement = new int[nbcelement];
                                if (cg_boco_read1(index_file, index_base, index_zone, index_bc, bcelement, null) != (int)Error.CG_OK) {
                                    ThrowError(index_file);
                                }
                                num[2] = index_bc;
                                if (pointset_type.Equals(PointSetType_t.ElementList) || pointset_type.Equals(PointSetType_t.ElementRange)) {
                                    gridlocation_type = GridLocation_t.FaceCenter;
                                } else {
                                    if (cg_golist(index_file, index_base, depht, labels, num) != (int)Error.CG_OK) {
                                        ThrowError(index_file);
                                    }
                                    if (cg_gridlocation_read(out gridlocation_type) != (int)Error.CG_OK) {
                                        ThrowError(index_file);
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
                        need_struct_conversion = true;
                        // insert here ...
                    }
                    x_of_all_zones[index_zone - 1] = x;
                    y_of_all_zones[index_zone - 1] = y;
                    z_of_all_zones[index_zone - 1] = z;
                    bcelements_of_all_zones[index_zone - 1] = bcelements;
                    normal_indices_of_all_zones[index_zone - 1] = normal_indices;
                    starts_of_all_zones[index_zone - 1] = starts;
                    ends_of_all_zones[index_zone - 1] = ends;
                    nbndrys_of_all_zones[index_zone - 1] = nbndrys;
                    pointset_types_of_all_zones[index_zone - 1] = pointset_types;
                    bc_types_of_all_zones[index_zone - 1] = bc_types;
                    element_types_of_all_zones[index_zone - 1] = element_types;
                    elements_of_all_zones[index_zone - 1] = elements;
                    boconames_of_all_zones[index_zone - 1] = boconames;
                    gridlocation_types_of_all_zones[index_zone - 1] = gridlocation_types;
                }
            }
            if (cg_close(index_file) != (int)Error.CG_OK) {
                ThrowError();
            }
        }
    }
}
