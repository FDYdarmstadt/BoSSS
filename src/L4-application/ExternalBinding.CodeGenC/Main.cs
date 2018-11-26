using System;
using ilPSP.ExternalBinding;
using System.IO;
using System.Reflection;
using System.Collections.Generic;

namespace ilPSP.ExternalBinding.CodeGenC
{
	class MainClass {
		public static void Main (string[] args) {
		
			var TypesToExprt = new Type[] { typeof(MsrMatrix_), typeof(IMutableMatrix_), typeof(Common_), typeof(IMutableMatrixEx_), typeof(ISparseMatrix_), typeof(ISparseSolver_), typeof(Partition_) };
			
			BindingDotC_1.WriteLine("/* Autogenerated code - do not modify */");
			BindingDotC_1.WriteLine("#include <assert.h>");
			//BindingDotC_1.WriteLine("#include <glib-2.0/glib.h>");
			BindingDotC_1.WriteLine("#include <mono/jit/jit.h>");
			BindingDotC_1.WriteLine("#include <mono/metadata/assembly.h>");
			BindingDotC_1.WriteLine("#include <mono/metadata/debug-helpers.h>");
			BindingDotC_1.WriteLine("#define DEFINE_MONKEY_INTERNALS");
			BindingDotC_1.WriteLine("#define DECLARE_BINDING");
			BindingDotC_1.WriteLine("#include \"monkey.h\"");
			
			BindingDotH.WriteLine("/* Autogenerated code - do not modify */");
			BindingDotH.WriteLine("#ifndef DECLARE_BINDING");
			BindingDotH.WriteLine("#define BINDING");
			BindingDotH.WriteLine("#else");
			BindingDotH.WriteLine("#define BINDING extern");
			BindingDotH.WriteLine("#endif");
			
			BindingDotC_2.WriteLine("void InitBinding() {");
            BindingDotC_2.WriteLine("MonoClass *klass;");
            BindingDotC_2.WriteLine("MonoMethodDesc* desc_X;");
			
			
			foreach(var t in TypesToExprt)
				BindingsForType(t);
			
			BindingDotH.WriteLine();
			BindingDotC_1.WriteLine();
			BindingDotC_2.WriteLine("}");
			BindingDotC_2.WriteLine();
			
			StreamWriter bindC = new StreamWriter("binding.c");
			bindC.Write(BindingDotC_1.ToString());
			bindC.Write(BindingDotC_2.ToString());
			bindC.WriteLine();
			bindC.Close();
			
			StreamWriter bindh = new StreamWriter("binding.h");
			bindh.Write(BindingDotH.ToString());
			bindh.Close();
		}
		
		static void BindingsForType(Type t) {
			// load class code:
			BindingDotC_2.WriteLine("klass = mono_class_from_name (img, \"" + t.Namespace + "\", \"" + t.Name + "\");");
			BindingDotC_2.WriteLine("assert(klass);");

			// get pointers for all methods:
			var methods = t.GetMethods(BindingFlags.Static | BindingFlags.Public);
			List<string> methodNames = new List<string>();
			foreach(var mi in methods) {
				
				EncodeMethod(mi, methodNames);
			}
		}
		
		
		
		static void EncodeMethod(MethodInfo mi, List<string> allNames) {
			// test 
			
			if(allNames.Contains(mi.Name))
				throw new ApplicationException("overloading (i mean: equal name -- different signature) is not allowed.");
			
			if(!mi.IsStatic)
				throw new ApplicationException("only static methods are supported.");

			// arguments
			StringWriter searchString = new StringWriter();
			searchString.Write(mi.DeclaringType.FullName + ":" + mi.Name + "(");
			
			var parameters = mi.GetParameters();
			string[] CArgs = new string[parameters.Length];
			string[] sStrArgs = new string[parameters.Length];
			for(int i = 0; i < parameters.Length; i++) {
				ParameterInfo pi = parameters[i];
				
				if(!pi.ParameterType.IsPointer && !pi.ParameterType.IsByRef) {
					throw new ApplicationException("only pointers or reference (ref, out) are allowed: '" + mi.DeclaringType.FullName + ":" + mi.Name + "'");	
				}
				
				switch(pi.ParameterType.Name) {
				case "Int32&": CArgs[i] = "int*"; sStrArgs[i] = "int&"; break;
				case "Double&": CArgs[i] = "double*"; sStrArgs[i] = "double&"; break;
				case "Byte*": CArgs[i] = "char*"; sStrArgs[i] = "byte*"; break;
				case "Double*": CArgs[i] = "double*"; sStrArgs[i] = "double*"; break;
				default: throw new ApplicationException("unsupported parameter type '" + pi.ParameterType.Name + "', method '" + mi.DeclaringType.FullName + ":" + mi.Name + "'");
				}
				
				CArgs[i] += " " + pi.Name;
				searchString.Write(sStrArgs[i]);
				if((i+1) < parameters.Length)
					searchString.Write(",");
			}
			searchString.Write(")");
			
			// declare mono variable
			string methname = mi.DeclaringType.Name + mi.Name;
			string varName = "p" + methname;
			BindingDotC_1.WriteLine("MonoMethod *" + varName + ";");
			
			// code for mono lookup
			BindingDotC_2.WriteLine( "desc_X = mono_method_desc_new (\"" + searchString + "\", 1);");
            BindingDotC_2.WriteLine( "assert(desc_X);");
		    BindingDotC_2.WriteLine( varName + " = mono_method_desc_search_in_class (desc_X, klass);");
            BindingDotC_2.WriteLine( "assert(" + varName + ");");
			
			// function signature
			string sig;
			{
			 	StringWriter sigW = new StringWriter();	
				sigW.Write("void " + methname + "(");
				for( int i = 0; i < CArgs.Length; i++) {
					sigW.Write(CArgs[i]);
					if((i+1) < parameters.Length)
						sigW.Write(",");
				}
				sigW.Write(")");
				sig = sigW.ToString();
			}
			
			// c-wrapper
			BindingDotC_1.WriteLine(sig + "{");
			BindingDotC_1.WriteLine("void* args["+ CArgs.Length + "];");
			for( int i = 0; i < CArgs.Length; i++)
				BindingDotC_1.WriteLine("args [" + i + "] = "+ parameters[i].Name +";");
    		BindingDotC_1.WriteLine("mono_runtime_invoke ("+varName+", NULL, args, NULL);");
			BindingDotC_1.WriteLine("}");
			
			// c-protoype
			BindingDotH.WriteLine(sig + ";");
			
			// all kinds of name mangling (to support fortran)
			NameManglingImpl(methname.ToLowerInvariant(), methname,CArgs,parameters);
			NameManglingImpl(methname.ToUpperInvariant(), methname,CArgs,parameters);
			NameManglingImpl(methname.ToLowerInvariant() + "_", methname,CArgs,parameters);
			NameManglingImpl(methname.ToUpperInvariant() + "_", methname,CArgs,parameters);
			NameManglingImpl("_" + methname.ToLowerInvariant(), methname,CArgs,parameters);
			NameManglingImpl("_" + methname.ToUpperInvariant(), methname,CArgs,parameters);
			NameManglingImpl("_" + methname, methname,CArgs,parameters);
			NameManglingImpl(methname + "_", methname,CArgs,parameters);
		}
		
		static void NameManglingImpl(string newName, string orgName, string[] Cargs, ParameterInfo[] Para) {
			BindingDotC_1.Write("void " + newName + "(");
			for( int i = 0; i < Cargs.Length; i++) {
				BindingDotC_1.Write(Cargs[i]);
				if((i+1) < Para.Length)
					BindingDotC_1.Write(",");
			}
			BindingDotC_1.Write(")");
			BindingDotC_1.Write("{");
			BindingDotC_1.Write(orgName);
			BindingDotC_1.Write("(");
			for( int i = 0; i < Cargs.Length; i++) {
				BindingDotC_1.Write( Para[i].Name );
				if( i < Cargs.Length - 1)
					BindingDotC_1.Write(",");
			}
			BindingDotC_1.WriteLine(");}");
		}
		
		
		
		static StringWriter BindingDotH = new StringWriter();

		static StringWriter BindingDotC_1 = new StringWriter();
		static StringWriter BindingDotC_2 = new StringWriter();
		
		
		
	}
}
