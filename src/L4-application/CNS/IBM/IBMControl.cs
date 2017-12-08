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

using BoSSS.Foundation.XDG;
using BoSSS.Platform.LinAlg;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using CNS.EquationSystem;
using System;

namespace CNS.IBM {

    /// <summary>
    /// Specialized control file for immersed boundary runs
    /// </summary>
    public class IBMControl : CNSControl {

        public IBMControl() : base() {
            Console.WriteLine("Warning: Auto-adding default estimator factory for IBM. Is this what you want?");
            this.DynamicLoadBalancing_CellCostEstimatorFactories.Add(delegate (IApplication<AppControl> app, int performanceClassCount) {
                if (performanceClassCount != 3) {
                    throw new ArgumentException(
                        "Only valid for exactly three performance classes",
                        nameof(performanceClassCount));
                }

                int[] map = new int[] { 1, 15, 0 };
                return new StaticCellCostEstimator(map);
            });
        }

        /// <summary>
        /// The name of the 'fluid' in the void domain
        /// </summary>
        public string VoidSpeciesName = "void";

        /// <summary>
        /// The name of the fluid in the physical domain
        /// </summary>
        public string FluidSpeciesName = "fluid";

        /// <summary>
        /// The level set function at a given point in time
        /// </summary>
        public Func<double[], double, double> LevelSetFunction;

        /// <summary>
        /// Determines which part of the domain is considered as physical
        /// domain. For example, <see cref="LevelsetSign.Positive"/> indicates
        /// that the positive level set region is occupied by fluid, while the
        /// negative level set region is considered void.
        /// </summary>
        public LevelsetSign FluidSpeciesSign = LevelsetSign.Positive;

        /// <summary>
        /// The edge tag associated with the zero level set
        /// </summary>
        public string LevelSetBoundaryTag;

        /// <summary>
        /// The (HMF) quadrature order to be used for the evaluation of volume
        /// and edge operators
        /// </summary>
        public int LevelSetQuadratureOrder;

        /// <summary>
        /// Volume fraction threshold for cell agglomeration
        /// </summary>
        public double AgglomerationThreshold = 0.0;

        /// <summary>
        /// If true, information about the number of agglomerated cells is
        /// displayed on the console
        /// </summary>
        public bool PrintAgglomerationInfo = false;

        /// <summary>
        /// If true, a file containing information about agglomeration pairs
        /// will be saved to the database
        /// </summary>
        public bool SaveAgglomerationPairs = false;

        /// <summary>
        /// The quadrature variant to be used
        /// </summary>
        public XQuadFactoryHelper.MomentFittingVariants MomentFittingVariant =
            XQuadFactoryHelper.MomentFittingVariants.Classic;

        public bool SurfaceHMF_ProjectNodesToLevelSet = false;

        public bool SurfaceHMF_RestrictNodes = false;

        public bool SurfaceHMF_UseGaussNodes = true;

        public bool VolumeHMF_RestrictNodes = false;

        public double VolumeHMF_NodeCountSafetyFactor = 1.0;

        public bool VolumeHMF_UseGaussNodes = true;

        public TimesteppingStrategies TimesteppingStrategy = TimesteppingStrategies.LieSplitting;

        public Func<double[], double, Vector3D> LevelSetVelocity;

        /// <summary>
        /// Verifies the control file
        /// </summary>
        public override void Verify() {
            base.Verify();

            if (DomainType == DomainTypes.Standard) {
                throw new Exception(
                    "Wrong domain type for immersed boundary runs");
            }

            if (LevelSetFunction == null) {
                if (RestartInfo != null && DomainType != DomainTypes.MovingImmersedBoundary) {
                    // Everything is fine (level set data is read from database
                    // and is not moving)
                } else {
                    throw new Exception(
                        "A level set function is required when running in IBM mode");
                }
            }

            if (!FieldOptions.ContainsKey(IBMVariables.LevelSet)) {
                throw new Exception(
                    "A level set is required when running in IBM mode");
            }
        }
    }
}
