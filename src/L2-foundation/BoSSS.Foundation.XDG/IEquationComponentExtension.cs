using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace BoSSS.Foundation.XDG {
    
    /// <summary>
    /// extension methods regarding XDG operators
    /// </summary>
    static public class IEquationComponentExtension_XDG {

        /// <summary>
        /// creates the spatial operator that consists only of component <paramref name="c"/>
        /// </summary>
        public static XSpatialOperator XOperator(this IEquationComponent c) {

            string[] Codomain = new string[] { "v1" };
            string[] Domain = c.ArgumentOrdering.ToArray();
            string[] Param = (c.ParameterOrdering != null) ? c.ParameterOrdering.ToArray() : new string[0];

            XSpatialOperator ret = new XSpatialOperator(Domain, Param, Codomain);
            ret.EquationComponents[Codomain[0]].Add(c);
            ret.Commit();

            return ret;
        }

    }
}
