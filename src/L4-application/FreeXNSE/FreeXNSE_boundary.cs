using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Control;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace FreeXNSE {

    /// <summary>
    /// types of boundary conditions for incompressible Navier-Stokes solvers, Free surface edition;
    /// </summary>
    public enum FreeXNSE_BcType {

        /// <summary>
        /// homogenous or inhomogenous Dirichlet
        /// </summary>
        Dirichlet = 0,

        /// <summary>
        /// homogenous or inhomogenous Neumann
        /// </summary>
        Neumann = 2,

        /// <summary>
        /// Robin
        /// </summary>
        Robin = 3   

    }

    public class FreeXNSE_BoundaryCondMap : BoundaryCondMap<FreeXNSE_BcType> {
        static string[] BndFunctions(IGridData g, string[] SpeciesNames) {
            int D = g.SpatialDimension;
            List<string> scalarFields = new List<string>();

            foreach(var S in SpeciesNames) {
                for(int d = 0; d < D; d++) {
                    scalarFields.Add(VariableNames.Velocity_d(d) + "#" + S);                    
                }
                scalarFields.Add(VariableNames.Pressure + "#" + S);
            }

            return scalarFields.ToArray();
        }

        /// <summary>
        /// Loops over all boundary conditions:
        /// If e.g. `VelocityX` is defined, but not `VelocityX#A`, the value for `VelocityX#A` is inferred from `VelocityX`.
        /// </summary>
        static IDictionary<string, AppControl.BoundaryValueCollection> BndyModify(IDictionary<string, AppControl.BoundaryValueCollection> b, string[] SpeciesNames) {

            bool isNonXname(string s) {
                foreach(var sn in SpeciesNames) {
                    string end = "#" + sn;
                    if(s.Length > 2 && s.EndsWith(end))
                        return false;
                }
                return true;
            }

            var ret = new Dictionary<string, AppControl.BoundaryValueCollection>();

            foreach(var kv in b) {
                var coll = kv.Value.CloneAs();
                ret.Add(kv.Key, coll);

                string[] definedKeys = coll.Evaluators.Keys.ToArray();
                foreach(var varName in definedKeys) {
                    if(isNonXname(varName)) {
                        foreach(var spc in SpeciesNames) {
                            string XvarName = varName + "#" + spc;
                            if(!kv.Value.Evaluators.ContainsKey(XvarName)) {
                                coll.Evaluators.Add(XvarName, coll.Evaluators[varName]);
                            }
                        }


                    }
                }
            }


            return ret;
        }

        /// <summary>
        /// ctor
        /// </summary>
        public FreeXNSE_BoundaryCondMap(IGridData f, IDictionary<string, AppControl.BoundaryValueCollection> b, string[] species)
            : base(f, BndyModify(b, species), BndFunctions(f, species)) {
            foreach(var pair in f.Grid.EdgeTagNames) {
                string legacyKey = pair.Value.ToLowerInvariant().Replace("dirichlet", "dirichlet_velocity_inlet").Replace("neumann", "neumann_pressure_outlet").Replace("robin", "robin_velocity_inlet");
                legacyEdgeTagNames.Add(pair.Key, legacyKey);
            }
            f.Grid.EdgeTagNames.Clear();

            // pretty dirty, expand edgenames to contain both old and new labeling
            foreach(var pair in legacyEdgeTagNames) {
                f.Grid.EdgeTagNames.Add(pair);
            }

            foreach(var entry in b) {
                string legacyKey = entry.Key.ToLowerInvariant().Replace("dirichlet", "dirichlet_velocity_inlet").Replace("neumann", "neumann_pressure_outlet").Replace("robin", "robin_velocity_inlet");
                //if(legacyKey.Contains("Robin")) throw new NotImplementedException();
                legacyCollection.Add(legacyKey, entry.Value);
            }
            
        }

        /// <summary>
        /// true if there is at least one pressure-outlet edge:
        /// <see cref="FreeXNSE_BcType.Neumann"/>;
        /// </summary>
        public bool DirichletPressureBoundary {
            get {
                return ((BCTypeUseCount[FreeXNSE_BcType.Neumann] > 0));
            }
        }

        IDictionary<byte, string> legacyEdgeTagNames = new Dictionary<byte, string> ();

        IDictionary<string, AppControl.BoundaryValueCollection> legacyCollection = new Dictionary<string, AppControl.BoundaryValueCollection>();
        public IncompressibleBoundaryCondMap TranslateToLegacy(IGridData grd) {
            IncompressibleBoundaryCondMap legacyBndMap = new IncompressibleMultiphaseBoundaryCondMap(grd, this.legacyCollection, new string[] { "A", "B" });
            return legacyBndMap;
        }
    }
}
