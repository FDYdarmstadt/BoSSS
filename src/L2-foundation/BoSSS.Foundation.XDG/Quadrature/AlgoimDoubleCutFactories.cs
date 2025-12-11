// Ignore Spelling: Algoim

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG.Quadrature.HMF;
using ilPSP;
using ilPSP.Utils;
using Newtonsoft.Json.Linq;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Linq;
using System.Security.Policy;

namespace BoSSS.Foundation.XDG.Quadrature.Algoim {
	/// <summary>
	/// Provides a factory configuration for processing double cut cells based on specified parameters.
	/// There are two types of integration: surface and volume
	/// There are two types of frames: cell inner domain and cell boundary (will be then transformed into edge rules)
	/// </summary>
	public class AlgoimDoubleCutFactories {


        /// <summary>
        /// This would return a factory object with the configuration 
        /// </summary>
        /// <param name="lvlsets">level set data</param>
        /// <param name="e">reference element</param>
        public AlgoimDoubleCutFactories(LevelSetTracker.LevelSetData[] lvlsets, RefElement e) {
            refElement = e;
            lsData = lvlsets;
			if (lsData.Length != 2)
				throw new ArgumentOutOfRangeException("Double cut cells need 2 level sets, not more not less!");

			// Get the double cut cells (the calculation is done for all double cut cells and species at once and cached)
			// The later queries are done returned from cache
			var cutDom1 = lsData[0].Region.GetCutCellMask4LevSet(0);
			var cutDom2 = lsData[1].Region.GetCutCellMask4LevSet(1);
			var cutDom = cutDom1.Intersect(cutDom2).ToGeometicalMask();
			var _cutDom = cutDom.Intersect(lsData[0].GridDat.Cells.GetCells4Refelement(e));
			allDoubleCutCells = _cutDom;
		}

		ExecutionMask allDoubleCutCells;

		LevelSetTracker.LevelSetData[] lsData;

        RefElement refElement;

        /// <summary>
        /// key: quadrature order <br/>
        /// value: quadrature rule
        /// </summary>
        Dictionary<int, ChunkRulePair<QuadRule>[][]> m_SurfaceRules = new Dictionary<int, ChunkRulePair<QuadRule>[][]>();

        /// <summary>
        /// key: quadrature order <br/>
        /// value: quadrature rule
        /// </summary>
        Dictionary<int, ChunkRulePair<QuadRule>[][]> m_VolumeRules = new Dictionary<int, ChunkRulePair<QuadRule>[][]>();

		/// <summary>
		/// key: quadrature order <br/>
		/// value: quadrature rule
		/// </summary>
		Dictionary<int, ChunkRulePair<QuadRule>[]> m_EdgeRules = new Dictionary<int, ChunkRulePair<QuadRule>[]>();

		/// <summary>
		/// key: quadrature order <br/>
		/// value: quadrature rule
		/// </summary>
		Dictionary<int, ChunkRulePair<CellBoundaryQuadRule>[][]> m_CellBoundarySurfaceRules = new Dictionary<int, ChunkRulePair<CellBoundaryQuadRule>[][]>();

		/// <summary>
		/// key: quadrature order <br/>
		/// value: quadrature rule
		/// </summary>
		Dictionary<int, ChunkRulePair<CellBoundaryQuadRule>[][]> m_CellBoundaryVolumeRules = new Dictionary<int, ChunkRulePair<CellBoundaryQuadRule>[][]>();

		public IQuadRuleFactory<QuadRule> GetSurfaceFactory(JumpTypes[] jumps) {
            var factory = new AlgoimDoubleCutSurfaceFactory() {
                m_Owner = this,
                m_Rules = this.m_SurfaceRules,
				m_Jumps = jumps
			};
            //factory.useMetrics = true;
            factory.m_CalculateQuadRule = ilPSP.Utils.Algoim.GetSurfaceQuadratureRules;
            return factory;
        }

        public IQuadRuleFactory<QuadRule> GetVolumeFactory(JumpTypes[] jumps) {
            var factory = new AlgoimDoubleCutVolumeFactory() {
                m_Owner = this,
                m_Rules = this.m_VolumeRules,
                m_Jumps = jumps
            };
            factory.m_CalculateQuadRule = ilPSP.Utils.Algoim.GetVolumeQuadratureRules;
            return factory;
        }

		public IQuadRuleFactory<CellBoundaryQuadRule> GetCellBoundarySurfaceFactory(JumpTypes[] jumps) {
			var factory = new AlgoimDoubleCutCellBoundarySurfFactory() {
				m_Owner = this,
				m_Rules = this.m_CellBoundarySurfaceRules,
				m_Jumps = jumps
			};
			//factory.useMetrics = true;
			factory.m_CalculateQuadRule = ilPSP.Utils.Algoim.GetSurfaceQuadratureRules;
            factory.m_CalculateQuadRuleSingle = ilPSP.Utils.Algoim.GetSurfaceQuadratureRules;

            return factory;
		}

		public IQuadRuleFactory<CellBoundaryQuadRule> GetCellBoundaryVolumeFactory(JumpTypes[] jumps) {
			var factory = new AlgoimDoubleCutCellBoundaryVolFactory() {
				m_Owner = this,
				m_Rules = this.m_CellBoundaryVolumeRules,
				m_Jumps = jumps
			};
			factory.m_CalculateQuadRule = ilPSP.Utils.Algoim.GetVolumeQuadratureRules;
			return factory;
		}

		abstract class AlgoimDoubleCutBaseFactory<TQuadRule> : IQuadRuleFactory<TQuadRule> where TQuadRule : QuadRule {
            internal AlgoimDoubleCutFactories m_Owner;

            /// <summary>
            /// cache for computed rules
            /// - key: quadrature order
            /// - values: quadrature rules
            ///   - 1st index: species index defined by <see cref="GetSpecies"/>
            ///   - 2nd index: enumeration over the items of <see cref="allDoubleCutCells"/>
            /// </summary>
            internal Dictionary<int, ChunkRulePair<TQuadRule>[][]> m_Rules;

            internal JumpTypes[] m_Jumps;

			public RefElement RefElement => m_Owner.refElement;

            internal ExecutionMask allDoubleCutCells => m_Owner.allDoubleCutCells;

			internal int spaceDim => RefElement.SpatialDimension;

			virtual internal int ruleDim => RefElement.SpatialDimension;

			internal int noOfSpeciesPermutations = 4;   //volume rules: number of species (all possible)
                                                        //surface rules: number of combinations for level sets (e.g., ls0=0; ls1>1)


            //public bool useMetrics = false;

            // only valid for single quadrature rule calculations (either surface or volume)
            internal GetQuadratureRule m_CalculateQuadRule;
            internal delegate UnsafeAlgoim.QuadScheme[] GetQuadratureRule(int dim, int p1, int p2, int q, int[] sizes1, int[] sizes2, double[] coordinates1, double[] coordinates2, double[] LSvalues1, double[] LSvalues2);

            


            public LevelSetTracker.LevelSetData[] lsData => m_Owner.lsData;

            public int[] GetCachedRuleOrders() {
                throw new NotImplementedException();
            }

            /// <summary>
            /// Returns the quadrature rule of given order for each cell in mask.
            /// </summary>
            /// <param name="mask">Cell mask for which quadrature nodes are requested</param>
            /// <param name="RequestedOrder">Order of quadrature</param>
            /// <returns>An array of quadrature rules for cells in mask</returns>
            /// <exception cref="ArgumentException"></exception>
            public IEnumerable<IChunkRulePair<TQuadRule>> GetQuadRuleSet(ExecutionMask mask, int RequestedOrder) {
                if (!(mask is CellMask))
                    throw new ArgumentException("Expecting a cell mask.");
                if(mask.MaskType != MaskType.Geometrical)
                    throw new ArgumentException("expecting a geometrical mask.");

                if(!m_Rules.ContainsKey(RequestedOrder)) {
                    CalculateQuadRuleSetSingle(RequestedOrder);
                }

				//which part of the quadRule (a quadRule for a cell is returned for all species/permutation. This is to indicate which one we need)
				int speciesIndex = GetSpecies(); 

				// check if all the mask or a submask is requested
				if (mask.NoOfItemsLocally == allDoubleCutCells.NoOfItemsLocally) {
                    if(mask.Except(allDoubleCutCells).NoOfItemsLocally > 0)
                        throw new ArgumentException("requested quadrature rules outside of supported range;");

					return m_Rules[RequestedOrder][speciesIndex];
				} else {
					var Rule = m_Rules[RequestedOrder][speciesIndex];
					int localLength = mask.NoOfItemsLocally, totalLength = Rule.Length;
					var Ret = new ChunkRulePair<TQuadRule>[localLength];
					int t = 0;
					int jsub = 0;
					foreach (int jCell in mask.ItemEnum) {
						Debug.Assert(Rule[t].Chunk.Len == 1);
						while (jCell > Rule[t].Chunk.i0) {
							t++;
						}

						Debug.Assert(jCell == Rule[t].Chunk.i0);
						Ret[jsub] = Rule[t];
#if DEBUG
						Ret[jsub].Rule.Weights.CheckForNanOrInf();
						Ret[jsub].Rule.Nodes.CheckForNanOrInf();
#endif
						jsub++;
					}
					Debug.Assert(jsub == localLength);

					return Ret;
				}
            }

            /// <summary>
            /// Get array index of the species requested for Algoim data. 
            /// <see cref="ilPSP.Utils.Algoim.GetSurfaceQuadratureRules(int, int, int, int, int[], int[], double[], double[], double[], double[])"/>
            /// <see cref="ilPSP.Utils.Algoim.GetVolumeQuadratureRules(int, int, int, int, int[], int[], double[], double[], double[], double[])"/>
            /// </summary>
            /// <returns></returns>
            abstract internal int GetSpecies();

			internal virtual (int, int[], double[], double[]) CreatePhiData(LevelSetTracker.LevelSetData lsData, int jCell) {
				//number of nodes in 1d for level set interpolation = degree + 1
				int n = (lsData.LevelSet as LevelSet).Basis.Degree + 1;

				//create Chebyshev nodes (must be identical with Algoim, otherwise leads to interpolation errors for high orders)
				double[] points = GenericBlas.ChebyshevNodesSecondKind(-1.0, 1.0, n);

				//Cartesian pair product for the points
				int numberOfCombinations = (int)Math.Pow(points.Length, spaceDim);
				MultidimensionalArray combinations = MultidimensionalArray.CreateCartesianPairProduct(points, spaceDim);

				//Get level set values on combinations
				NodeSet NS = new NodeSet(RefElement, combinations, false);
				var ret = lsData.GetLevSetValues(NS, jCell, 1);

				//Convert to 1D Array for the Wrapper
				double[] y = new double[ret.Length];
				for (int i = 0; i < ret.Length; i++)
					y[i] = ret[0, i];

				// create the double array for the coordinates that level set are queried (1D version, for details see AlgoimWrapper)
				double[] x = Enumerable.Repeat(points, spaceDim).SelectMany(i => i).ToArray();

				// create the double array for the number of points in each coordinate
				int[] sizes = Enumerable.Repeat(points.Length, spaceDim).ToArray();

				return (n, sizes, x, y);
			}

            protected void GetMonomialExpression(LevelSetTracker.LevelSetData lsData, int j) {
                int N = 10; // 1,   x, y, z,   x*x, x*y, y*y, x*z, y*z, z*z 
                int D = this.RefElement.SpatialDimension;
                var Nodes = new NodeSet(this.RefElement, N, D, false);
                var rnd = new Random(0);
                for(int n = 0; n < N; n++) {
                    for(int d = 0; d < D; d++) {
                        Nodes[n, d] = rnd.NextDouble();
                    }
                }
                Nodes.LockForever();

                var lsVal = lsData.GetLevSetValues(Nodes, j, 1);

                var LHS = MultidimensionalArray.Create(N, 10);
                var RHS = new double[N];
                for(int n = 0; n < N; n++) {
                    var x = Nodes[n, 0];
                    var y = Nodes[n, 1];
                    var z = Nodes[n, 2];
                    LHS.SetRow(n, new double[] { 1, x, y, z, x * x, x * y, y * y, x * z, y * z, z * z });
                    RHS[n] = lsVal[0, n];
                }

                var MomomialCoeffs = LHS.Solve(RHS);

                Console.WriteLine($"  1*({MomomialCoeffs[0]}) + x*({MomomialCoeffs[1]}) + y*({MomomialCoeffs[2]}) + z*({MomomialCoeffs[3]}) + x*x*({MomomialCoeffs[4]}) + x*y*({MomomialCoeffs[5]}) + y*y*({MomomialCoeffs[6]}) + x*z*({MomomialCoeffs[7]}) + y*z*({MomomialCoeffs[8]}) + z*z*({MomomialCoeffs[9]})");

            }

            /// <summary>
            /// slight modified version to query points on the face
            /// </summary>
            /// <param name="lsData"></param>
            /// <param name="jCell"></param>
            /// <param name="faceIndex"></param>
            /// <returns></returns>
            internal (int, int[], double[], double[], int sign) CreatePhiDataForFace(LevelSetTracker.LevelSetData lsData, int jCell, int faceIndex) {
				//number of nodes in 1d for level set interpolation = degree + 1
				int n = (lsData.LevelSet as LevelSet).Basis.Degree + 1;

				//create Chebyshev nodes (must be identical with Algoim, otherwise leads to interpolation errors for high orders)
				double[] points = GenericBlas.ChebyshevNodesSecondKind(-1.0, 1.0, n);

				int numberOfCombinations = points.Length.Pow(ruleDim); //Cartesian pair product for the points
				MultidimensionalArray combinationsOnEdge = MultidimensionalArray.CreateCartesianPairProduct(points, ruleDim);

				// Transform edge based coordinates to cell based coordinates (to query level sets)
				MultidimensionalArray combinations = MultidimensionalArray.Create(numberOfCombinations, spaceDim);  // to transform edge coordinates to cell, create an array in the original space dim
				RefElement.TransformFaceCoordinates(faceIndex, combinationsOnEdge, combinations);

				NodeSet NS = new NodeSet(RefElement, combinations, false);
				var ret = lsData.GetLevSetValues(NS, jCell, 1);

				//Convert to 1D Array for the Wrapper
				double[] y = new double[ret.Length];
                bool isPos = false;
                bool isNeg = false;
                for(int i = 0; i < ret.Length; i++) {
                    double yi = ret[0, i];
                    if(yi < 0)
                        isNeg = true;
                    if(yi > 0)
                        isPos = true;
                    y[i] = yi;
                }
                int sign;
                if(isPos == true && isNeg == true) {
                    sign = 0;
                } else if(isPos == true && isNeg == false) {
                    sign = +1;
                } else if(isPos == false && isNeg == true) {
                    sign = -1;
                } else {
                    // (isPos == false && isNeg == false) 
                    throw new ArithmeticException($"level-set is exactly zero - should have been caught by the special case handling for level-sets which coincide with mesh edges; level-set-values are between {y.Min()} and {y.Max()}");
                }

                // create the double array for the coordinates that level set are queried (1D version, for details see AlgoimWrapper)
                double[] x = Enumerable.Repeat(points, spaceDim).SelectMany(i => i).ToArray();

				// create the double array for the number of points in each coordinate
				int[] sizes = Enumerable.Repeat(points.Length, spaceDim).ToArray();

				return (n, sizes, x, y, sign);
			}

            internal void CalculateQuadRuleSetSingle(int RequestedOrder) {
                List<ChunkRulePair<TQuadRule>>[] ret = new List<ChunkRulePair<TQuadRule>>[4];

                foreach(int cell in allDoubleCutCells.ItemEnum) {
                    var quadRules = GetNodesAndWeights(cell, RequestedOrder);
                    for(int q = 0; q < quadRules.Length; q++) {
                        if(ret[q] == null)
                            ret[q] = new List<ChunkRulePair<TQuadRule>>();
                        quadRules[q].OrderOfPrecision = RequestedOrder;
                        ret[q].Add(new ChunkRulePair<TQuadRule>(Chunk.GetSingleElementChunk(cell), quadRules[q]));
                    }
                }
                var retArray = ret.Select(list => list?.ToArray() ?? new ChunkRulePair<TQuadRule>[0]).ToArray();
                m_Rules.Add(RequestedOrder, retArray);
            }

			/// <summary>
			/// Get the quadrature nodes and weights for a cell. (either surface or volume, declared at the instantiation of the factory)
			/// </summary>
			/// <param name="jCell">local index of the cell</param>
			/// <param name="RequestedOrder">requested order for quadrature</param>
			/// <returns></returns>
			abstract internal TQuadRule[] GetNodesAndWeights(int jCell, int RequestedOrder);

			internal void SetNodesAndWeight(UnsafeAlgoim.QuadScheme qs, QuadRule quadRule) {
				for (int row = 0; row < qs.length; row++) {
					quadRule.Weights[row] = qs.weights[row];

					for (int d = 0; d < qs.dimension; d++) { // map 1d array back to 2d
						int ind = row * qs.dimension + d;
						quadRule.Nodes[row, d] = qs.nodes[ind];
					}
				}
				quadRule.Nodes.LockForever();
			}
		}

        class AlgoimDoubleCutSurfaceFactory : AlgoimDoubleCutBaseFactory<QuadRule> {

            override internal int GetSpecies() {
                if (m_Jumps[0] == JumpTypes.Implicit && m_Jumps[1] == JumpTypes.OneMinusHeaviside)
                    return 0;
                else if (m_Jumps[0] == JumpTypes.Implicit && m_Jumps[1] == JumpTypes.Heaviside)
                    return 1;
                else if (m_Jumps[0] == JumpTypes.OneMinusHeaviside && m_Jumps[1] == JumpTypes.Implicit)
                    return 2;
                else if (m_Jumps[0] == JumpTypes.Heaviside && m_Jumps[1] == JumpTypes.Implicit)
                    return 3;
                else
                    throw new ArgumentOutOfRangeException();
            }

            




            /// <summary>
            /// Get the quadrature nodes and weights for a cell. (either surface or volume, declared at the instantiation of the factory)
            /// </summary>
            /// <param name="jCell">local index of the cell</param>
            /// <param name="RequestedOrder">requested order for quadrature</param>
            /// <returns></returns>
            internal override QuadRule[] GetNodesAndWeights(int jCell, int RequestedOrder) {
                (int n1, int[] sizes1, double[] x1, double[] y1) = CreatePhiData(lsData[0], jCell);
                (int n2, int[] sizes2, double[] x2, double[] y2) = CreatePhiData(lsData[1], jCell);

                



                // get the quadrature rule from the wrapper
                UnsafeAlgoim.QuadScheme[] qsArray = m_CalculateQuadRule(spaceDim, n1, n2, RequestedOrder, sizes1, sizes2, x1, x2, y1, y2);
                Debug.Assert(qsArray.Length == 2);


                /*
                UnsafeAlgoim.QuadScheme[] qsArray_alt = null;
                if(this.RefElement.SpatialDimension == 3)
                    qsArray_alt = m_CalculateQuadRule(spaceDim, n2, n1, RequestedOrder, sizes2, sizes1, x2, x1, y2, y1);
                */


                //Returned rules need to be categorized as ls0=0 ls1<0; ls0=0 ls1>0; ls0<0 ls1=0; ls0>0 ls1=0 
                QuadRule[] quadRuleArray = new QuadRule[4]; //for all possible combinations of level sets with each other
                for (int i = 0; i < qsArray.Length; i++) {
                    var qs = qsArray[i];
                    // If quadrature rule is empty, return.
                    if (qs.length < 1) {
                        quadRuleArray[2 * i] = CreateEmptyQuadRule();
                        quadRuleArray[2 * i + 1] = CreateEmptyQuadRule();
                        continue;
                    }

                    // Create quadrature rule and copy from the scheme
                    QuadRule quadRule = QuadRule.CreateBlank(RefElement, qs.length, qs.dimension);

                    for(int row = 0; row < qs.length; row++) {
                        quadRule.Weights[row] = qs.weights[row];
                        for (int d = 0; d < qs.dimension; d++) { // map 1d array back to 2d
                            int ind = row * qs.dimension + d;
                            quadRule.Nodes[row, d] = qs.nodes[ind];
                        }
                    }
                    quadRule.Nodes.LockForever();

                    /*
                    if(this.RefElement.SpatialDimension == 3) {
                        // Create quadrature rule and copy from the scheme
                        //var qs_alt = qsArray_alt[i];
                        //QuadRule quadRule_alt = QuadRule.CreateBlank(RefElement, qs_alt.length, qs_alt.dimension);

                        //for(int row = 0; row < qs_alt.length; row++) {
                        //    quadRule_alt.Weights[row] = qs_alt.weights[row];
                        //    for(int d = 0; d < qs_alt.dimension; d++) { // map 1d array back to 2d
                        //        int ind = row * qs_alt.dimension + d;
                        //        quadRule_alt.Nodes[row, d] = qs_alt.nodes[ind];
                        //    }
                        //}
                        //quadRule_alt.Nodes.LockForever();

                        GetMonomialExpression(lsData[0], jCell);
                        GetMonomialExpression(lsData[1], jCell);



                        var gd = this.allDoubleCutCells.GridData;
                        var test = new CompositeQuadRule<QuadRule>(null, gd);
                        test.chunkRulePairs.Add(new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(jCell), quadRule));
                        ICompositeQuadRule_Ext.SaveToTextFileCell(test, lsData[0].Region.GetCutCellMask().GridData, "j" + jCell + "_" + i + ".csv", false);

                        //var test_alt = new CompositeQuadRule<QuadRule>(null, gd);
                        //test_alt.chunkRulePairs.Add(new ChunkRulePair<QuadRule>(Chunk.GetSingleElementChunk(jCell), quadRule_alt));
                        //ICompositeQuadRule_Ext.SaveToTextFileCell(test_alt, lsData[0].Region.GetCutCellMask().GridData, "j" + jCell + "_" + i + "_alt.csv", false);

                    }*/
                    

                    //ApplyMetrics(quadRule, jCell);

                    var (negativeRule, positiveRule) = DivideQuadRules(lsData[1 - i], jCell, quadRule); //1-i: the other level set since we have two level sets
                    quadRuleArray[2 * i] = negativeRule;
                    quadRuleArray[2 * i + 1] = positiveRule;
                }

                return quadRuleArray;
            }

            /// <summary>
            /// Divides the quadrature rule depending on the value of level set into two
            /// </summary>
            /// <param name="lsData"></param>
            /// <param name="jCell"></param>
            /// <param name="quadRule"></param>
            /// <returns></returns>
            private (QuadRule negativeRule, QuadRule positiveRule) DivideQuadRules(LevelSetTracker.LevelSetData lsData, int jCell, QuadRule quadRule) {
                //weight and nodes at negative level set values
                List<double> negativeWeight = new List<double>();
                List<double[]> negativeNodes = new List<double[]>();
                //weight and nodes at positive level set values
                List<double> positiveWeight = new List<double>();
                List<double[]> positiveNodes = new List<double[]>();

                var ret = lsData.GetLevSetValues(quadRule.Nodes, jCell, 1);
                for (int i = 0; i < ret.Length; i++) {
                    var ls = ret[0, i];
                    if (ls <= 0) {
                        negativeWeight.Add(quadRule.Weights[i]);
                        negativeNodes.Add(quadRule.Nodes.GetRow(i));
                    } else if (ls > 0) {
                        positiveWeight.Add(quadRule.Weights[i]);
                        positiveNodes.Add(quadRule.Nodes.GetRow(i));
                    }
                }
                QuadRule negRule = QuadRule.CreateBlank(RefElement, negativeWeight.Count(), spaceDim);
                for (int n = 0; n < negativeWeight.Count(); n++) {
                    negRule.Weights[n] = negativeWeight[n];
                    negRule.Nodes.SetRow(n, negativeNodes[n]);
                }
                negRule.Nodes.LockForever();
                QuadRule posRule = QuadRule.CreateBlank(RefElement, positiveWeight.Count(), spaceDim);
                for (int n = 0; n < positiveWeight.Count(); n++) {
                    posRule.Weights[n] = positiveWeight[n];
                    posRule.Nodes.SetRow(n, positiveNodes[n]);
                }
                posRule.Nodes.LockForever();
                Debug.Assert(positiveWeight.Count() + negativeWeight.Count() == ret.Length);
                return (negRule, posRule);
            }

            private QuadRule CreateEmptyQuadRule() {
                QuadRule quadRuleEmpty = QuadRule.CreateBlank(RefElement, 1, spaceDim);
                quadRuleEmpty.Nodes.LockForever();
                return quadRuleEmpty;
            }
		}

		class AlgoimDoubleCutVolumeFactory : AlgoimDoubleCutBaseFactory<QuadRule> {

            override internal int GetSpecies() {
                if (m_Jumps[0] == JumpTypes.OneMinusHeaviside && m_Jumps[1] == JumpTypes.OneMinusHeaviside)
                    return 0;
                else if (m_Jumps[0] == JumpTypes.OneMinusHeaviside && m_Jumps[1] == JumpTypes.Heaviside)
					return 1;
				else if (m_Jumps[0] == JumpTypes.Heaviside && m_Jumps[1] == JumpTypes.OneMinusHeaviside)
					return 2;
				else if (m_Jumps[0] == JumpTypes.Heaviside && m_Jumps[1] == JumpTypes.Heaviside)
                    return 3;
                else
                    throw new ArgumentOutOfRangeException();
            }

			/// <summary>
			/// Get the quadrature nodes and weights for a cell. (either surface or volume, declared at the instantiation of the factory)
			/// </summary>
			/// <param name="jCell">local index of the cell</param>
			/// <param name="RequestedOrder">requested order for quadrature</param>
			/// <returns></returns>
			internal override QuadRule[] GetNodesAndWeights(int jCell, int RequestedOrder) {
				(int n1, int[] sizes1, double[] x1, double[] y1) = CreatePhiData(lsData[0], jCell);
				(int n2, int[] sizes2, double[] x2, double[] y2) = CreatePhiData(lsData[1], jCell);

                /*
                if(this.RefElement.SpatialDimension == 3) {
                    base.GetMonomialExpression(lsData[0], jCell);
                    base.GetMonomialExpression(lsData[1], jCell);
                }*/

                // get the quadrature rule from the wrapper
                UnsafeAlgoim.QuadScheme[] qsArray = m_CalculateQuadRule(spaceDim, n1, n2, RequestedOrder, sizes1, sizes2, x1, x2, y1, y2);
				QuadRule[] quadRuleArray = new QuadRule[qsArray.Length];
                for(int i=0; i < qsArray.Length; i++) { //for each permutation/species
                    var qs = qsArray[i];
                    // If quadrature rule is empty, return.
                    if (qs.length < 1) {
						quadRuleArray[i] = CreateEmptyQuadRule();
                        continue;
                    }

					// Create quadrature rule and copy from the scheme
					QuadRule quadRule = QuadRule.CreateBlank(RefElement, qs.length, qs.dimension);
                    SetNodesAndWeight(qs, quadRule);
					quadRuleArray[i] = quadRule;
				}
				return quadRuleArray;
			}

			private QuadRule CreateEmptyQuadRule() {
				QuadRule quadRuleEmpty = QuadRule.CreateBlank(RefElement, 1, spaceDim);
				quadRuleEmpty.Nodes.LockForever();
				return quadRuleEmpty;
			}
		}

		class AlgoimDoubleCutCellBoundaryVolFactory : AlgoimDoubleCutBaseFactory<CellBoundaryQuadRule> {

			override internal int GetSpecies() {
				if (m_Jumps[0] == JumpTypes.OneMinusHeaviside && m_Jumps[1] == JumpTypes.OneMinusHeaviside)
					return 0;
				else if (m_Jumps[0] == JumpTypes.OneMinusHeaviside && m_Jumps[1] == JumpTypes.Heaviside)
					return 1;
				else if (m_Jumps[0] == JumpTypes.Heaviside && m_Jumps[1] == JumpTypes.OneMinusHeaviside)
					return 2;
				else if (m_Jumps[0] == JumpTypes.Heaviside && m_Jumps[1] == JumpTypes.Heaviside)
					return 3;
				else
					throw new ArgumentOutOfRangeException();
			}

			// Edge rule is D-1 for the original space D
			internal override int ruleDim => spaceDim - 1;

			/// <summary>
			/// Combines different edge rules into a single CellBoundaryQuadRule rule.
			/// </summary>
			/// <param name="rules">An array of edge rules to be combined.</param>
			/// <returns>A combined <c>CellBoundaryQuadRule</c> with arranged <c>NumbersOfNodesPerFace</c>.</returns>
			CellBoundaryQuadRule CombineEdgeRules(QuadRule[] rules) {
				int numberOfNodes = 0;
				foreach (QuadRule rule in rules) {
					numberOfNodes += rule.NoOfNodes;
				}
				CellBoundaryQuadRule combinedRule = CellBoundaryQuadRule.CreateEmpty(RefElement, numberOfNodes, spaceDim, RefElement.NoOfFaces);
				int subarrayPointer = 0;
				for (int i = 0; i < rules.Length; ++i) {
					int subNumberOfNodes = rules[i].NoOfNodes;
					combinedRule.NumbersOfNodesPerFace[i] = subNumberOfNodes;
					for (int j = 0; j < subNumberOfNodes; ++j) {
						for (int d = 0; d < spaceDim; ++d) {
							combinedRule.Nodes[subarrayPointer + j, d] = rules[i].Nodes[j, d];
						}
						combinedRule.Weights[subarrayPointer + j] = rules[i].Weights[j];
					}
					subarrayPointer += subNumberOfNodes;
				}
				combinedRule.Nodes.LockForever();
				return combinedRule;
			}

			internal override CellBoundaryQuadRule[] GetNodesAndWeights(int jCell, int RequestedOrder) {
				int noOfEdges = RefElement.NoOfFaces;

				//we have two indices: first one is for the species (e.g., A, B or between A-B and A-C) and second one for edge index (e.g., top, bottom)
				QuadRule[][] allEdgesAndSpecies = new QuadRule[noOfSpeciesPermutations][];
				for (int specI = 0; specI < noOfSpeciesPermutations; specI++) {
					allEdgesAndSpecies[specI] = new QuadRule[noOfEdges];
				}

				for (int e = 0; e < noOfEdges; e++) {
					var retAllSpecies = GetNodesAndWeightsOnEdge(jCell, RequestedOrder, e); //edge specific query (notice that this is not the global edge index but only the type of edges as in RefElement)
					for (int specI = 0; specI < retAllSpecies.Length; specI++) {
						allEdgesAndSpecies[specI][e] = retAllSpecies[specI];
					}
				}
				//combine edges in for each permutation/species
				CellBoundaryQuadRule[] ret = new CellBoundaryQuadRule[noOfSpeciesPermutations];
				for (int specI = 0; specI < noOfSpeciesPermutations; specI++) {
					var combinedEdgeRules = CombineEdgeRules(allEdgesAndSpecies[specI]);
					ret[specI] = combinedEdgeRules;
				}
				return ret;
			}

			CellBoundaryQuadRule[] GetNodesAndWeightsOnEdge(int jCell, int RequestedOrder, int edgeIndex) {
				(int n1, int[] sizes1, double[] x1, double[] y1, int signPhi1) = CreatePhiDataForFace(lsData[0], jCell, edgeIndex);
				(int n2, int[] sizes2, double[] x2, double[] y2, int signPhi2) = CreatePhiDataForFace(lsData[1], jCell, edgeIndex);

				// get the quadrature rule from the wrapper
				UnsafeAlgoim.QuadScheme[] qsArray = m_CalculateQuadRule(ruleDim, n1, n2, RequestedOrder, sizes1, sizes2, x1, x2, y1, y2);
				CellBoundaryQuadRule[] quadRuleArray = new CellBoundaryQuadRule[qsArray.Length];

				for (int i = 0; i < qsArray.Length; i++) { //for each permutation/species
					var qs = qsArray[i];

					// If quadrature rule is empty, return.
					if (qs.length < 1) {
						quadRuleArray[i] = CreateEmptyQuadRule();
						continue;
					}

                    // Create rule (edge based)
					QuadRule quadRuleOnEdge = QuadRule.CreateBlank(RefElement.FaceRefElement, qs.length, qs.dimension);
                    SetNodesAndWeight(qs, quadRuleOnEdge);

					// Algoim returns an edge based rule, it must be converted to a CellBoundaryQuadRule (cell based coordinates)
					CellBoundaryQuadRule quadRule = CellBoundaryQuadRule.CreateEmpty(RefElement, qs.length, spaceDim, RefElement.NoOfFaces);
					quadRule.Weights = quadRuleOnEdge.Weights;
					RefElement.TransformFaceCoordinates(edgeIndex, quadRuleOnEdge.Nodes, quadRule.Nodes); // to transform edge rule back to cell coordinates (since we are creating CellBoundaryQuadRule, cell based rule)
					quadRule.Nodes.LockForever();

					//if (useMetrics) {
					//	var metrics = lsData[0].GetLevelSetNormalReferenceToPhysicalMetrics(quadRule.Nodes, jCell, 1);
					//	for (int k = 0; k < qs.length; k++)
					//		quadRule.Weights[k] /= metrics[0, k];
					//}

					quadRuleArray[i] = quadRule;
				}

				return quadRuleArray;
			}


			CellBoundaryQuadRule CreateEmptyQuadRule() {
				CellBoundaryQuadRule quadRuleEmpty = CellBoundaryQuadRule.CreateEmpty(RefElement, 1, spaceDim, RefElement.NoOfFaces);
				quadRuleEmpty.Nodes.LockForever();
				return quadRuleEmpty;
			}
		}

		class AlgoimDoubleCutCellBoundarySurfFactory : AlgoimDoubleCutBaseFactory<CellBoundaryQuadRule> {

            internal GetQuadratureRuleSingle m_CalculateQuadRuleSingle;
            internal delegate UnsafeAlgoim.QuadScheme GetQuadratureRuleSingle(int dim, int p, int q, int[] lengths, double[] x, double[] y);


            override internal int GetSpecies() {
				if (m_Jumps[0] == JumpTypes.Implicit && m_Jumps[1] == JumpTypes.OneMinusHeaviside)
					return 0;
				else if (m_Jumps[0] == JumpTypes.Implicit && m_Jumps[1] == JumpTypes.Heaviside)
					return 1;
				else if (m_Jumps[0] == JumpTypes.OneMinusHeaviside && m_Jumps[1] == JumpTypes.Implicit)
					return 2;
				else if (m_Jumps[0] == JumpTypes.Heaviside && m_Jumps[1] == JumpTypes.Implicit)
					return 3;
				else
					throw new ArgumentOutOfRangeException();
			}

			// Edge rule is D-1 for the original space D
			internal override int ruleDim => spaceDim - 1;

			/// <summary>
			/// Combines different edge rules into a single CellBoundaryQuadRule rule.
			/// </summary>
			/// <param name="rules">An array of edge rules to be combined.</param>
			/// <returns>A combined <c>CellBoundaryQuadRule</c> with arranged <c>NumbersOfNodesPerFace</c>.</returns>
			CellBoundaryQuadRule CombineEdgeRules(QuadRule[] rules) {
				int numberOfNodes = 0;
				foreach (QuadRule rule in rules) {
					numberOfNodes += rule?.NoOfNodes ?? 0;
				}
				CellBoundaryQuadRule combinedRule = CellBoundaryQuadRule.CreateEmpty(RefElement, numberOfNodes, spaceDim, RefElement.NoOfFaces);
				int subarrayPointer = 0;
				for (int i = 0; i < rules.Length; ++i) {
                    if(rules[i] == null)
                        continue;
					int subNumberOfNodes = rules[i].NoOfNodes;
					combinedRule.NumbersOfNodesPerFace[i] = subNumberOfNodes;
					for (int j = 0; j < subNumberOfNodes; ++j) {
						for (int d = 0; d < spaceDim; ++d) {
							combinedRule.Nodes[subarrayPointer + j, d] = rules[i].Nodes[j, d];
						}
						combinedRule.Weights[subarrayPointer + j] = rules[i].Weights[j];
					}
					subarrayPointer += subNumberOfNodes;
				}
				combinedRule.Nodes.LockForever();
				return combinedRule;
			}

			internal override CellBoundaryQuadRule[] GetNodesAndWeights(int jCell, int RequestedOrder) {
				//we have two indices: first one is for the species (e.g., A, B or between A-B and A-C) and second one for face index (e.g., top, bottom)
				QuadRule[][] allEdgesAndSpecies = new QuadRule[noOfSpeciesPermutations][];
				int noOfFaces = RefElement.NoOfFaces;
                for(int i = 0; i < allEdgesAndSpecies.Length; i++) {
                    allEdgesAndSpecies[i] = new QuadRule[noOfFaces];
                }


				for (int e = 0; e < noOfFaces; e++) {
					var retAllSpecies = GetNodesAndWeightsOnFace(jCell, RequestedOrder, e); //edge specific query (notice that this is not the global edge index but only the type of edges as in RefElement)
					for (int specI = 0; specI < retAllSpecies.Length; specI++) {
						allEdgesAndSpecies[specI][e] = retAllSpecies[specI];
					}
				}
				//combine edges in for each permutation/species
				CellBoundaryQuadRule[] ret = new CellBoundaryQuadRule[noOfSpeciesPermutations];
				for (int specI = 0; specI < noOfSpeciesPermutations; specI++) {
					var combinedEdgeRules = CombineEdgeRules(allEdgesAndSpecies[specI]);
					ret[specI] = combinedEdgeRules;
				}
                return ret;
			}

			CellBoundaryQuadRule[] GetNodesAndWeightsOnFace(int jCell, int RequestedOrder, int faceIndex) {
				(int n1, int[] sizes1, double[] x1, double[] y1, int signPhi1) = CreatePhiDataForFace(lsData[0], jCell, faceIndex);
				(int n2, int[] sizes2, double[] x2, double[] y2, int signPhi2) = CreatePhiDataForFace(lsData[1], jCell, faceIndex);

                //Returned rules need to be categorized as ls0=0 ls1<0; ls0=0 ls1>0; ls0<0 ls1=0; ls0>0 ls1=0 
                CellBoundaryQuadRule[] quadRuleArray = new CellBoundaryQuadRule[4]; //for all possible combinations of level sets with each other

                if(signPhi1 == 0 && signPhi2 == 0) {

                    // get the quadrature rule from the wrapper
                    UnsafeAlgoim.QuadScheme[] qsArray = m_CalculateQuadRule(ruleDim, n1, n2, RequestedOrder, sizes1, sizes2, x1, x2, y1, y2);
                    
                    for(int i = 0; i < qsArray.Length; i++) { //for each permutation/species
                        var qs = qsArray[i];

                        // If quadrature rule is empty, continue.
                        if(qs.length < 1) {
                            quadRuleArray[2 * i] = CreateEmptyQuadRule();
                            quadRuleArray[2 * i + 1] = CreateEmptyQuadRule();
                            continue;
                        }

                        // Create rule (edge based)
                        QuadRule quadRuleOnEdge = QuadRule.CreateBlank(RefElement.FaceRefElement, qs.length, qs.dimension);
                        SetNodesAndWeight(qs, quadRuleOnEdge);

                        // Algoim returns an edge based rule, it must be converted to a CellBoundaryQuadRule (cell based coordinates)
                        CellBoundaryQuadRule quadRule = TransfromFromFace2Cell(faceIndex, qs, quadRuleOnEdge);

                        var (negativeRule, positiveRule) = DivideQuadRules(lsData[1 - i], jCell, quadRule); //1-i: the other level set since we have two level sets
                        quadRuleArray[2 * i] = negativeRule;
                        quadRuleArray[2 * i + 1] = positiveRule;
                    }

                    

                } else if(signPhi1 == 0 && signPhi2 != 0) {
                    // ++++++++++++++++++++++++++++++++++++++
                    // second level-set is empty on this face
                    // ++++++++++++++++++++++++++++++++++++++

                    UnsafeAlgoim.QuadScheme qs = m_CalculateQuadRuleSingle(ruleDim, n1, RequestedOrder, sizes1, x1, y1);

                    CellBoundaryQuadRule quadRule;
                    if(qs.length < 1) {
                        quadRule = CreateEmptyQuadRule();
                    } else {
                        // Create rule (edge based)
                        QuadRule quadRuleOnEdge = QuadRule.CreateBlank(RefElement.FaceRefElement, qs.length, qs.dimension);
                        SetNodesAndWeight(qs, quadRuleOnEdge);

                        // Algoim returns an edge based rule, it must be converted to a CellBoundaryQuadRule (cell based coordinates)
                        quadRule = TransfromFromFace2Cell(faceIndex, qs, quadRuleOnEdge);
                    }

                    quadRuleArray[2] = CreateEmptyQuadRule();
                    quadRuleArray[3] = CreateEmptyQuadRule();

                    if(signPhi2 > 0) {
                        // ls1>0
                        quadRuleArray[0] = CreateEmptyQuadRule();
                        quadRuleArray[1] = quadRule;
                    } else {
                        // ls1<0
                        quadRuleArray[0] = quadRule;
                        quadRuleArray[1] = CreateEmptyQuadRule(); 
                    }

                } else if(signPhi1 != 0 && signPhi2 == 0) {
                    // ++++++++++++++++++++++++++++++++++++++
                    // first level-set is empty on this face
                    // ++++++++++++++++++++++++++++++++++++++

                    UnsafeAlgoim.QuadScheme qs = m_CalculateQuadRuleSingle(ruleDim, n2, RequestedOrder, sizes2, x2, y2);

                    CellBoundaryQuadRule quadRule;
                    if(qs.length < 1) {
                        quadRule = CreateEmptyQuadRule();
                    } else {
                        // Create rule (edge based)
                        QuadRule quadRuleOnEdge = QuadRule.CreateBlank(RefElement.FaceRefElement, qs.length, qs.dimension);
                        SetNodesAndWeight(qs, quadRuleOnEdge);

                        // Algoim returns an edge based rule, it must be converted to a CellBoundaryQuadRule (cell based coordinates)
                        quadRule = TransfromFromFace2Cell(faceIndex, qs, quadRuleOnEdge);
                    }

                    quadRuleArray[0] = CreateEmptyQuadRule();
                    quadRuleArray[1] = CreateEmptyQuadRule();

                    if(signPhi1 > 0) {
                        // ls0>0
                        quadRuleArray[2] = CreateEmptyQuadRule();
                        quadRuleArray[3] = quadRule;
                    } else {
                        // ls0<0
                        quadRuleArray[2] = quadRule;
                        quadRuleArray[3] = CreateEmptyQuadRule();
                    }


                } else if(signPhi1 != 0 && signPhi2 != 0) {
                    // ++++++++++++++++++++++++++++++++++++++
                    // both level-set are empty on this face
                    // ++++++++++++++++++++++++++++++++++++++


                    quadRuleArray[0] = CreateEmptyQuadRule();
                    quadRuleArray[1] = CreateEmptyQuadRule();
                    quadRuleArray[2] = CreateEmptyQuadRule();
                    quadRuleArray[3] = CreateEmptyQuadRule();
                }

                return quadRuleArray;

                CellBoundaryQuadRule TransfromFromFace2Cell(int faceIndex, UnsafeAlgoim.QuadScheme qs, QuadRule quadRuleOnEdge) {
                    CellBoundaryQuadRule quadRule = CellBoundaryQuadRule.CreateEmpty(RefElement, qs.length, spaceDim, RefElement.NoOfFaces);
                    quadRule.Weights = quadRuleOnEdge.Weights;
                    RefElement.TransformFaceCoordinates(faceIndex, quadRuleOnEdge.Nodes, quadRule.Nodes); // to transform edge rule back to cell coordinates (since we are creating CellBoundaryQuadRule, cell based rule)
                    quadRule.NumbersOfNodesPerFace[faceIndex] = quadRule.Weights.Length;
                    quadRule.Nodes.LockForever();
                    return quadRule;
                }
            }

            /// <summary>
            /// Divides the quadrature rule depending on the value of level set into two
            /// </summary>
            /// <param name="lsData"></param>
            /// <param name="jCell"></param>
            /// <param name="quadRule"></param>
            /// <returns></returns>
            private (CellBoundaryQuadRule negativeRule, CellBoundaryQuadRule positiveRule) DivideQuadRules(LevelSetTracker.LevelSetData lsData, int jCell, CellBoundaryQuadRule quadRule) {
                //weight and nodes at negative level set values
                List<double> negativeWeight = new List<double>();
                List<double[]> negativeNodes = new List<double[]>();
                //weight and nodes at positive level set values
                List<double> positiveWeight = new List<double>();
                List<double[]> positiveNodes = new List<double[]>();
                int NoofFaces = quadRule.NumbersOfNodesPerFace.Length;
                int[] negative_NumbersOfNodesPerFace = new int[NoofFaces];
                int[] positive_NumbersOfNodesPerFace = new int[NoofFaces];

                var ret = lsData.GetLevSetValues(quadRule.Nodes, jCell, 1);
                int k = 0;
                for(int iFace = 0; iFace < quadRule.NumbersOfNodesPerFace.Length; iFace++) {
                    for(int i = 0; i < quadRule.NumbersOfNodesPerFace[iFace]; i++) {
                        var ls = ret[0, k];
                        if(ls < 0) {
                            negativeWeight.Add(quadRule.Weights[i]);
                            negativeNodes.Add(quadRule.Nodes.GetRow(i));
                            negative_NumbersOfNodesPerFace[iFace]++;
                        } else if(ls > 0) {
                            positiveWeight.Add(quadRule.Weights[i]);
                            positiveNodes.Add(quadRule.Nodes.GetRow(i));
                            positive_NumbersOfNodesPerFace[iFace]++;
                        }
                        k++;
                    }
                }
                CellBoundaryQuadRule negRule = CellBoundaryQuadRule.CreateEmpty(RefElement, negativeWeight.Count(), spaceDim, NoofFaces);
                for(int n = 0; n < negativeWeight.Count(); n++) {
                    negRule.Weights[n] = negativeWeight[n];
                    negRule.Nodes.SetRow(n, negativeNodes[n]);
                }
                negRule.NumbersOfNodesPerFace = negative_NumbersOfNodesPerFace;
                negRule.Nodes.LockForever();
                CellBoundaryQuadRule posRule = CellBoundaryQuadRule.CreateEmpty(RefElement, positiveWeight.Count(), spaceDim, NoofFaces);
                for(int n = 0; n < positiveWeight.Count(); n++) {
                    posRule.Weights[n] = positiveWeight[n];
                    posRule.Nodes.SetRow(n, positiveNodes[n]);
                }
                negRule.NumbersOfNodesPerFace = negative_NumbersOfNodesPerFace;
                posRule.Nodes.LockForever();
                Debug.Assert(positiveWeight.Count() + negativeWeight.Count() == ret.Length);
                return (negRule, posRule);
            }



            CellBoundaryQuadRule CreateEmptyQuadRule() {
				CellBoundaryQuadRule quadRuleEmpty = CellBoundaryQuadRule.CreateEmpty(RefElement, 0, spaceDim, RefElement.NoOfFaces);
				quadRuleEmpty.Nodes.LockForever();
				return quadRuleEmpty;
			}
		}

	}





}


