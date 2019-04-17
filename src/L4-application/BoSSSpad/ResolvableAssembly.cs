using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.BoSSSpad
{
    /// <summary>
    /// Enables user to resolve Dlls from arbitrary path
    /// </summary>
    public class ResolvableAssembly
    {
        /// <summary>
        /// Enables user to resolve Dlls from arbitrary path
        /// </summary>
        /// <param name="path">to Dll folder</param>
        public ResolvableAssembly(string path)
        {
            this.path = path;
            //Find dlls in own folder if called from ElectronBoSSSpad
            AppDomain.CurrentDomain.AssemblyResolve += CurrentDomain_AssemblyResolve;
        }

        string path = null;

        /// <summary>
        /// Resolves assembly not found exceptions. 
        /// This happens when the exe must look for dlls in additional folders
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="args"></param>
        /// <returns></returns>
        private System.Reflection.Assembly CurrentDomain_AssemblyResolve(object sender, ResolveEventArgs args)
        {
            // Ignore missing resources
            if (args.Name.Contains(".resources"))
                return null;

            // check for assemblies already loaded

            System.Reflection.Assembly assembly = AppDomain.CurrentDomain.GetAssemblies().
                FirstOrDefault(a => a.FullName == args.Name);
            if (assembly != null)
                return assembly;

            // Try to load by filename - split out the filename of the full assembly name
            // and append the base path of the original assembly (ie. look in the same dir)
            string filename = args.Name.Split(',')[0] + ".dll".ToLower();

            string asmFile = System.IO.Path.Combine(path, filename);

            try
            {
                return System.Reflection.Assembly.LoadFrom(asmFile);
            }
            catch (Exception)
            {
                return null;
            }
        }
    }
}
