using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Text;

namespace BoSSS.Solution.Utils {
    
    /// <summary>
    /// Driver routines for integral evaluation
    /// </summary>
    public static class IntegrationDriverRoutine {


        /// <summary>
        /// Driver interface for numerical integration;
        /// 
        /// For each cell $ K_j $ in <paramref name="cm"/>, this computes
        /// \[
        ///   I_j = \int_{K_j \cap S} f(\vec{x}, u_1(\vec{x}), ... , u_n(\vec{x}) ) dV 
        /// \]
        /// where:
        /// - $ S $ is the domain of the species <paramref name="Spc"/>
        /// - the position $ \vec{x} $
        /// - $ u_i $ are the DG fields provided by <paramref name="U"/>
        /// - $ f(\vec{x}, u_1(\vec{x}), ... , u_n(\vec{x}) ) $ is the function <paramref name="func"/>
        /// <param name="lsTrk">
        /// required XDG tracker <see cref="LevelSetTracker.GetXDGSpaceMetrics"/>
        /// </param>
        /// </summary>
        /// <param name="Spc">
        /// species index
        /// </param>
        /// <param name="cm">
        /// optional cell-mask; if null, the entire domain is taken
        /// </param>
        /// <param name="func">
        /// the integrand
        /// </param>
        /// <param name="U">
        /// DG fields on which <paramref name="func"/> might depend on
        /// </param>
        /// <returns>
        /// the integrand value for each cell
        /// </returns>
        public double[] VolumeIntegral(this LevelSetTracker lsTrk, SpeciesId Spc, CellMask cm, Func<Vector, double[], double> func, int quadOrder, params DGField[] U) {

        }

        /// <summary>
        /// Driver interface for numerical integration;
        /// 
        /// For each cell $ K_j $ in <paramref name="cm"/>, this computes
        /// \[
        ///   I_j = \oint_{\partial K_j \cap S} f(\vec{x}, \vec{n}, u_1(\vec{x}), ... , u_n(\vec{x}) ) dV 
        /// \]
        /// where:
        /// - $ S $ is the domain of the species <paramref name="Spc"/>
        /// - the position $ \vec{x} $
        /// - the outer cell normal $ \vec{n} $
        /// - $ u_i $ are the DG fields provided by <paramref name="U"/>
        /// - $ f(\vec{x}, \vec{n}, u_1(\vec{x}), ... , u_n(\vec{x}) ) $ is the function <paramref name="func"/>, which might depend on the position $ \vec{x} $ and the values of the DG fields
        /// </summary>
        /// <param name="lsTrk">
        /// required XDG tracker <see cref="LevelSetTracker.GetXDGSpaceMetrics"/>
        /// </param>
        /// <param name="Spc">
        /// species index
        /// </param>
        /// <param name="cm">
        /// optional cell-mask; if null, the entire domain is taken
        /// </param>
        /// <param name="func">
        /// the integrand
        /// </param>
        /// <param name="U">
        /// DG fields on which <paramref name="func"/> might depend on
        /// </param>
        /// <returns>
        /// the integrand value for each cell
        /// </returns>
        public double[] CellSurfaceIntegral(this LevelSetTracker lsTrk, SpeciesId Spc, CellMask cm, Func<Vector, Vector, double[], double> func, int quadOrder, params DGField[] U) {

        }


        /// <summary>
        /// Driver interface for numerical integration;
        /// 
        /// For each cell $ K_j $ in <paramref name="cm"/>, this computes
        /// \[
        ///   I_j = \oint_{\mathfrak{I} \cap \overline{S_1} \cap \overline{S_2}} f(\vec{x}, \vec{n}, u_1(\vec{x}), ... , u_n(\vec{x}) ) dV 
        /// \]
        /// where:
        /// - $ S $ is the domain of the species <paramref name="Spc"/>
        /// - the position $ \vec{x} $
        /// - the level-set/interface normal $ \vec{n} $, pointing from species <paramref name="SpcA"/> to <paramref name="SpcB"/>
        /// - $ u_i $ are the DG fields provided by <paramref name="U"/>
        /// - $ f(\vec{x}, \vec{n}, u_1(\vec{x}), ... , u_n(\vec{x}) ) $ is the function <paramref name="func"/>, which might depend on the position $ \vec{x} $ and the values of the DG fields
        /// </summary>
        /// <param name="lsTrk">
        /// required XDG tracker <see cref="LevelSetTracker.GetXDGSpaceMetrics"/>
        /// </param>
        /// <param name="iLevSet">
        /// level-set/interface index
        /// </param>
        /// <param name="SpcA">
        /// first species index
        /// </param>
        /// <param name="SpcA">
        /// second species index
        /// </param>
        /// <param name="cm">
        /// optional cell-mask; if null, the entire domain is taken
        /// </param>
        /// <param name="func">
        /// the integrand
        /// </param>
        /// <param name="U">
        /// DG fields on which <paramref name="func"/> might depend on
        /// </param>
        /// <returns>
        /// the integrand value for each cell
        /// </returns>
        public double[] InterfaceIntegral(this LevelSetTracker lsTrk, int iLevSet, SpeciesId SpcA, SpeciesId SpcB, CellMask cm, Func<Vector, Vector, double[], double> func, params DGField[] U) {

        }






    }
}
