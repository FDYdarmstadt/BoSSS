using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Runtime.InteropServices;
using System.Linq.Expressions;
using MPI.Wrappers;

namespace BoSSS.Application.ExternalBinding
{
	public static  class Infrastructure {


		/*
		static SortedDictionary<string,Delegate> delegates = new SortedDictionary<string, Delegate>();
		
		public static IntPtr Addr(string name) {
			
	        Delegate del = null;
			
			delegates.TryGetValue(name, out del);
			
			if( del == null) {
    	        var methInf = (typeof(ExternalBinding)).GetMethod(name);
            
				Console.WriteLine("elo: " + name);
				Console.WriteLine("method found: " + (methInf != null));
				
	            var paramInfoS = methInf.GetParameters();
            	Type[] typeArgs = new Type[paramInfoS.Length];
        	    for (int i = 0; i < paramInfoS.Length; i++) {
    	            typeArgs[i] = paramInfoS[i].ParameterType;
	            }
			
        	    // builds a delegate type
    	        Type delegateType = Expression.GetActionType(typeArgs);
	            del = Delegate.CreateDelegate(delegateType, methInf);
				delegates.Add(name,del);
			}
			
            IntPtr ret = Marshal.GetFunctionPointerForDelegate(del);
            Console.WriteLine("return value = " + ret);
            return ret;
        }
		*/
		
		/// <summary>
		/// all objects; in future, maybe a smarter thing than an infinity growing list should be used
		/// </summary>
		private static List<object> objects = new List<object>();
		
		
		internal static int RegisterObject(object o) {
			if(o==null)
				throw new ArgumentNullException();
			objects.Add(o);
			return objects.Count -1;
		}
		
		internal static void ForgetObject(int oIndex) {
			objects[oIndex] = null;	
		}
		
		internal static object GetObject(int oIndex) {
			return objects[oIndex];	
		}
		
		internal static int ErrorHandler(Exception e) {
            int rank, size;
			csMPI.Raw.Comm_Rank(csMPI.Raw._COMM.WORLD, out rank);
            csMPI.Raw.Comm_Size(csMPI.Raw._COMM.WORLD, out size);
			Console.Error.WriteLine("=========================================================");
			Console.Error.WriteLine(".NET- Exception at MPI rank " + rank + " of " + size);
			Console.Error.WriteLine(e.GetType().FullName + ": " + e.Message);
			Console.Error.WriteLine("Stacktrace:");
			Console.Error.WriteLine(e.StackTrace.ToString());
			Console.Error.WriteLine("=========================================================");
			return -1;	
		}

        unsafe internal static string FortranToDotNET(byte* line, byte termChar) {
            byte[] mini = new byte[2];
            string mLine = "";
            fixed (byte* pMini = mini) {
                while (*line != termChar) {
                    *pMini = *line;
                    string s = Marshal.PtrToStringAnsi((IntPtr)pMini);
                    mLine += s;
                    line++;
                }
            }
            return mLine;
        }

	}
}

