using System;
using System.Collections.Generic;
using ilPSP.LinSolvers;
using System.Runtime.InteropServices;
using System.Xml;
using System.IO;
using System.Runtime.Serialization.Formatters.Binary;

namespace ilPSP.ExternalBinding {
	public static class Common_ {
				
		/// <summary>
		/// releases all C# - references on an object, making it ready for garbage-collection.
		/// </summary>
		public static void ReleaseObject(ref int _ref, out int ierr) {
			ierr = 0;
			try {
				object o = Infrastructure.GetObject(_ref);	
				if( o is IDisposable) 
					((IDisposable)o).Dispose();
				Infrastructure.ForgetObject(_ref);
			} catch (Exception e) {
				ierr = Infrastructure.ErrorHandler(e);
			}
		}
		
		public static void Test(ref int  dummy, out int ierr)  {
			ierr = 0;
			Console.WriteLine("Hello from Test: " + (dummy));	
		}	
		
		static bool mustFinalizeMPI;
		
		public static void ilPSPInitialize() {
            ilPSP.Environment.Bootstrap(new string[0], out mustFinalizeMPI);
			//Console.WriteLine("elo from " +  Enviroment.MPIEnv.MPI_Rank + " of " + Enviroment.MPIEnv.MPI_Size);
		}
		
		public static void ilPSPFinalize() {
			if(mustFinalizeMPI)
				MPI.Wrappers.csMPI.Raw.mpiFinalize();
		}
	}
	
	public static class MsrMatrix_ {
		
		#region MsrMatrix_Bindings
		
		public static void New(out int _ref, ref int LocalNumberOfRows, ref int NoOfColumns, ref int _RowsPerBlock, ref int _ColPerBlock, out int ierr) {
			ierr = 0;
			try {
				_ref = Infrastructure.RegisterObject(new MsrMatrix(LocalNumberOfRows,NoOfColumns,_RowsPerBlock,_ColPerBlock));
			} catch (Exception e) {
				ierr = Infrastructure.ErrorHandler(e);
				_ref = -1;
			}
		}
		
		#endregion
		
	}
	
	public static class ISparseMatrix_  {
		public static void GetI0(ref int _ref, out int i0, out int ierr)  {
			ierr = 0;
			i0 = -1;
			try {
				ISparseMatrix M = (ISparseMatrix) Infrastructure.GetObject(_ref);
				i0 = (int) M.RowPartitioning.i0;
			} catch (Exception e) {
				ierr = Infrastructure.ErrorHandler(e);
			}
		}

		public static void GetLocLen(ref int _ref, out int LocLen, out int ierr)  {
			ierr = 0;
			LocLen = -1;
			try {
				ISparseMatrix M = (ISparseMatrix) Infrastructure.GetObject(_ref);
				LocLen = M.RowPartitioning.LocalLength;
			} catch (Exception e) {
				ierr = Infrastructure.ErrorHandler(e);
			}
		}
		
		public static void GetRowPart(ref int MtxRef, out int PartRef, out int ierr) {
			ierr = 0;
			PartRef = -1;
			try {
				ISparseMatrix M = (ISparseMatrix) Infrastructure.GetObject(MtxRef);
				PartRef = Infrastructure.RegisterObject(M.RowPartitioning);
			} catch (Exception e) {
				ierr = Infrastructure.ErrorHandler(e);
			}
		}
	}	
	
	public static class Partition_ {
		public static void GetI0(ref int _ref, out int i0, out int ierr)  {
			ierr = 0;
			i0 = -1;
			try {
				IPartitioning p = (IPartitioning) Infrastructure.GetObject(_ref);
				i0 = (int) p.i0;
			} catch (Exception e) {
				ierr = Infrastructure.ErrorHandler(e);
			}
		}

		public static void GetLocLen(ref int _ref, out int LocLen, out int ierr)  {
			ierr = 0;
			LocLen = -1;
			try {
				IPartitioning p = (IPartitioning) Infrastructure.GetObject(_ref);
				LocLen = (int) p.LocalLength;
			} catch (Exception e) {
				ierr = Infrastructure.ErrorHandler(e);
			}
		}		

		public static void GetMPIrank(ref int _ref, out int MPIrank, out int ierr)  {
			ierr = 0;
			MPIrank = -1;
			try {
				IPartitioning p = (IPartitioning) Infrastructure.GetObject(_ref);
				MPIrank = (int) p.MpiRank;
			} catch (Exception e) {
				ierr = Infrastructure.ErrorHandler(e);
			}
		}		

		public static void GetMPIsize(ref int _ref, out int MPIsize, out int ierr)  {
			ierr = 0;
			MPIsize = -1;
			try {
				IPartitioning p = (IPartitioning) Infrastructure.GetObject(_ref);
				MPIsize = (int) p.Size;
			} catch (Exception e) {
				ierr = Infrastructure.ErrorHandler(e);
			}
		}		
		
	}
	
	
	public static class IMutableMatrix_ {
		
		#region IMutableMatrix_Bindings
		
		public static void SetEntry(ref int _ref, ref int i, ref int j, ref double val, out int ierr) {
			ierr = 0;
			try {
				IMutableMatrix M = (IMutableMatrix) Infrastructure.GetObject(_ref);
				M[i,j] = val;
			} catch (Exception e) {
				ierr = Infrastructure.ErrorHandler(e);
			}
		}
		#endregion
	}
	
	public static class IMutableMatrixEx_ {
		
		unsafe public static void SaveToTextFileSparse(ref int _ref, byte* path, out int ierr) {
			ierr = 0;
			try {
				string _paht = Marshal.PtrToStringAnsi((IntPtr) path);
				IMutableMatrixEx M = (IMutableMatrixEx) Infrastructure.GetObject(_ref);
				//Console.WriteLine("saving: >" + _paht + "<");
				M.SaveToTextFileSparse(_paht);
			} catch (Exception e) {
				ierr = Infrastructure.ErrorHandler(e);
			}
		}
		
		unsafe public static void SaveToTextFileSparseF(ref int _ref, byte* termchar, byte* path, out int ierr) {
			ierr = 0;
			try {
				string _path  = Infrastructure.FortranToDotNET(path,*termchar);
				
				IMutableMatrixEx M = (IMutableMatrixEx) Infrastructure.GetObject(_ref);
				M.SaveToTextFileSparse(_path);
			} catch (Exception e) {
				ierr = Infrastructure.ErrorHandler(e);
			}
		}
	}
	
	public static class ISparseSolver_ {
		
		static void FromXMLInternal(out int ref_, string xmlString) {
			xmlString = "<?xml version=\"1.0\" encoding=\"utf-8\"?> \n" +
					"<document>\n" +
					xmlString +
				    "</document>";
			//Console.WriteLine(xmlString);
     		XmlDocument docu = new XmlDocument();
			docu.LoadXml(xmlString);
				
			XmlNode node = docu.SelectSingleNode("/document/sparsesolver");
				             
			var config = ilPSP.LinSolvers.SolverConfigurationParser.parseXMLConfigurationElements(node);
				
			var solver = SolverFactory.CreateSolver<ISparseSolver>(config);
				
			ref_ = Infrastructure.RegisterObject(solver);
		}
	
		unsafe public static void FromXML(out int ref_, byte* xmlcode, out int ierr) {
			ref_ = -1;
			ierr = 0;
			try {
				string xmlString = Marshal.PtrToStringAnsi((IntPtr) xmlcode);
				FromXMLInternal(out ref_, xmlString);
			} catch (Exception e) {
				ierr = Infrastructure.ErrorHandler(e);
			}
		}
		
		/// <summary>
		/// 
		/// </summary>
		/// <param name="lineTerm">
		/// character that terminates the line in the subsequent calls to <see cref="FromXMLSubmit"/>
		/// </param>
		/// <param name="ierr">
		/// 0 in case of success. otherwise an error code.
		/// </param>
		unsafe public static void FromXMLBegin(byte* lineTerm, out int ierr) { 
			ierr = 0;
			try {
				LineTerm = *lineTerm;
				XMLCode = "";
			} catch (Exception e) {
				ierr = Infrastructure.ErrorHandler(e);	
			}
		}
		
		static string XMLCode = null;
		static byte LineTerm = 0;
		
		unsafe public static void FromXMLSubmit(byte* line, out int ierr) {
			ierr = 0;
			try {
				if(XMLCode == null)
					throw new ApplicationException("FromXMLBegin(...) must be called before.");
				
				
				string mLine = Infrastructure.FortranToDotNET(line,LineTerm);
				mLine += "\n";
				XMLCode += mLine;
								
			} catch (Exception e) {
				ierr = Infrastructure.ErrorHandler(e);	
			}
		}

		public static void FromXMLEnd(out int SolverRef, out int ierr) {
			ierr = 0;
			SolverRef = -1;
			try {
				if(XMLCode == null)
					throw new ApplicationException("FromXMLBegin(...) must be called before.");
				
				FromXMLInternal(out SolverRef, XMLCode);
				
			} catch (Exception e) {
				ierr = Infrastructure.ErrorHandler(e);	
			} finally {
				XMLCode = null;	
			}
		}
		
		static int FromXMLFileInternal(string file, string name) {
			XmlDocument docu = new XmlDocument();
			docu.Load(file);
				
			XmlNode node = docu.SelectSingleNode("/SolverConfig/sparsesolver[@name='" + name + "']");
				             
			var config = ilPSP.LinSolvers.SolverConfigurationParser.parseXMLConfigurationElements(node);
			var solver = SolverFactory.CreateSolver<ISparseSolver>(config);
			return Infrastructure.RegisterObject(solver);
		}
		
		public unsafe static void FromXMLFileF(out int SolverRef, byte* termChar, byte* File, byte* solvername, out int ierr) {
			ierr = 0;
			SolverRef = -1;
			try {
				string _File = Infrastructure.FortranToDotNET(File,*termChar);
				string _Name = Infrastructure.FortranToDotNET(solvername,*termChar);
				
				SolverRef = FromXMLFileInternal(_File,_Name);
				
			} catch (Exception e) {
				ierr = Infrastructure.ErrorHandler(e);	
			} 
		
		}
		
		public unsafe  static void FromXMLFile(out int SolverRef, byte* File, byte* solvername, out int ierr) {
			ierr = 0;
			SolverRef = -1;
			try {
				string _File = Marshal.PtrToStringAnsi((IntPtr)File);
				string _Name = Marshal.PtrToStringAnsi((IntPtr)solvername);
				
				SolverRef = FromXMLFileInternal(_File,_Name);
				
			} catch (Exception e) {
				ierr = Infrastructure.ErrorHandler(e);	
			} 
		}
				
		
		unsafe public static void DefineMatrix(ref int SolverRef, ref int MatrixRef, out int ierr) {
			ierr = 0;
			try {
				ISparseSolver solver = (ISparseSolver) Infrastructure.GetObject(SolverRef);
				IMutableMatrixEx matrix = (IMutableMatrixEx) Infrastructure.GetObject(MatrixRef);
				
				solver.DefineMatrix(matrix);
			
			} catch (Exception e) {
				ierr = Infrastructure.ErrorHandler(e);
			}
		}
		
		
		
		unsafe public static void Solve(ref int SolverRef, ref int N, double* x, double* rhs, out int ierr) {
			ierr = 0;
			try {
				ISparseSolver solver = (ISparseSolver) Infrastructure.GetObject(SolverRef);
				// create managed buffers from unmanaged ones
				// (ginge auch effizienter, aber ...)
				double[] _x = new double[N];
				double[] _rhs = new double[N];
				unsafe {
					Marshal.Copy((IntPtr)x,_x,0,N);	
					Marshal.Copy((IntPtr)rhs,_rhs,0,N);
						
					// call solver	
					solver.Solve(_x, _rhs);
				
					Marshal.Copy(_x,0,(IntPtr)x,N);
					Marshal.Copy(_rhs,0,(IntPtr)rhs,N);
				}
				
			} catch (Exception e) {
				ierr = Infrastructure.ErrorHandler(e);
			}
		}
	
	}
}

