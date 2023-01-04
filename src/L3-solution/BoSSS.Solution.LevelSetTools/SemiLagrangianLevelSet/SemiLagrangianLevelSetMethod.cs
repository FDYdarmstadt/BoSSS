using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.Statistic;
using ilPSP;
using MPI.Wrappers;

namespace BoSSS.Solution.LevelSetTool.SemiLagrangianLevelSet
{
    public abstract class SemiLagrangianLevelSetMethod<PointsofCellType, PointType, CorrectionMaskType>
        where PointsofCellType : SemiLagrangianLevelSetMethod<PointsofCellType, PointType, CorrectionMaskType>.ItemsofCell
        where PointType : SemiLagrangianLevelSetMethod<PointsofCellType, PointType, CorrectionMaskType>.SinglePoint
        where CorrectionMaskType : SemiLagrangianLevelSetMethod<PointsofCellType, PointType, CorrectionMaskType>.CorrectionMask
    {
        public readonly int SpatialDimension;
        protected readonly VectorField<SinglePhaseField> Velocity_Current;
        protected readonly VectorField<SinglePhaseField> Velocity_Next;
        protected readonly SinglePhaseField Interface;
        protected readonly LevelSetTracker InterfaceTracker;
        protected readonly VectorField<SinglePhaseField> InterfaceLevelSetGradient;
        protected readonly FieldEvaluation Evaluator;

        protected readonly int ReseedingInterval;
        public enum MinimalDistanceSearchMode { FullSearch, GlobalTreeFilterSearch, CellTreeFilterSearch, TREEFILTERTEST};
        public enum TopologyMergingMode { On, Off};
        public enum NormalVectorDampingMode { Off, MaxAngle90, MaxAngle45, MaxAngle45_avg}
        protected readonly MinimalDistanceSearchMode MinimalDistanceSearch;
        protected readonly TopologyMergingMode TopologyMerging;
        protected static NormalVectorDampingMode NormalVectorDamping;
        protected readonly int NarrowBandWidth;
        public IGridData Grid;
        protected Random CoordGen = new Random();
        private readonly int MaxCellIndexSearchDepth = 100;
        protected CellMask NarrowBandCells_Old;
        protected CellMask NarrowBandCells;
        protected CellMask ReseedingBandCells;
        protected int ReseedingWdith;
        protected K_D_TreeFilter SearchTree;

        protected readonly int TargetNbrofPointsPerCellPerDimension;
        public readonly int UpperLimitPointsPerCellPerDimension;
        public readonly int LowerLimitPointsPerCellPerDimension;

        protected Stopwatch sw = new Stopwatch ();

        public SemiLagrangianLevelSetMethod (VectorField<SinglePhaseField> Velocity_Next, VectorField<SinglePhaseField> Velocity_Current, SinglePhaseField Interface, LevelSetTracker InterfaceTracker,
            VectorField<SinglePhaseField> InterfaceLevelSetGradient, int NarrowBandWidth, IGridData Grid, int ReseedingInterval,
            MinimalDistanceSearchMode MinimalDistanceSearch, TopologyMergingMode TopologyMerging, NormalVectorDampingMode NormalVectorDamping,
            int TargetNbrofPointsPerCellPerDimension = 0,int UpperLimitPointsPerCellPerDimension = 0, int LowerLimitPointsPerCellPerDimension = 0)
        {

            if (ReseedingInterval <= -2) throw new Exception("The Reseeding Interval must be higher than zero");
            if (NarrowBandWidth > InterfaceTracker.NearRegionWidth) throw new Exception("FieldWidth must not be larger than the LevelSetTrackers Region Width");

            SpatialDimension = Grid.SpatialDimension;
            this.Velocity_Next = Velocity_Next;
            this.Velocity_Current = Velocity_Current;
            this.Interface = Interface;
            this.InterfaceTracker = InterfaceTracker;
            this.InterfaceLevelSetGradient = InterfaceLevelSetGradient;
            Evaluator = new FieldEvaluation((GridData)Grid);

            this.ReseedingInterval = ReseedingInterval;
            this.NarrowBandWidth = NarrowBandWidth;
            this.Grid = Grid;
            this.MinimalDistanceSearch = MinimalDistanceSearch;
            this.TopologyMerging = TopologyMerging;
            SemiLagrangianLevelSetMethod<PointsofCellType, PointType, CorrectionMaskType>.NormalVectorDamping = NormalVectorDamping;

            this.TargetNbrofPointsPerCellPerDimension = TargetNbrofPointsPerCellPerDimension;
            this.UpperLimitPointsPerCellPerDimension = UpperLimitPointsPerCellPerDimension;
            this.LowerLimitPointsPerCellPerDimension = LowerLimitPointsPerCellPerDimension;
        }

        /// <summary>
    	/// ItemsofCell
	    /// Basic class for storing the various types of points of one cell.
    	/// Each type of points is represented by a dictionary entry with the key as the type identifier and the value as the point storage.
	    /// The MaxAdvectionStepLength stores the maximum step length of all points after each advection step.
	    /// AllPointsInactive is used to mark a cell in which all points are inactive.
    	/// If a k-d-Tree is constructed based on the points in one cell it is stored in TreeFilter. 
	    /// </summary>
        public abstract class ItemsofCell
        {
            public int CellIndex; 
            public BoundingBox Box;
            public double MaxAdvectionStepLength;
            public bool AllPointsInactive = false;
            public K_D_TreeFilter TreeFilter;
            public Dictionary<int, List<PointType>> ItemList;

	        /// <summary>
	        /// Initializes a new instance of the
	        /// <see cref="T:SemiLagrangianLevelSet.SemiLagrangianLevelSetMethod`3.ItemsofCell"/> class.
	        /// The types of points can be defined by putting the identifiers the PointSignsArray.
	        /// </summary>
	        /// <param name="CellIndex">Cell index.</param>
	        /// <param name="PointSigns">Point signs.</param>
            protected ItemsofCell(int CellIndex, int[] PointSigns = null)
            {
                this.CellIndex = CellIndex;
                if (PointSigns != null)
                {
                    ItemList = new Dictionary<int, List<PointType>>();
                    foreach (int sign in PointSigns)
                    {
			            if (ItemList.ContainsKey (sign))
			                throw new Exception ("Each type of point needs to have a different identifier");
                        ItemList.Add(sign, new List<PointType>());
                    }
                }
            }

	        /// <summary>
	        /// This function is used together with PointsCoordinateArray to evaluate all points of a cell in one step	
	        /// It takes the resulting values of a evaluation at each point and
	        /// adjusts the parameters of each points based on the resulting values
	        /// Has to be implemented seperately because different methods have different point parameters. 
	        /// </summary>
	        /// <param name="ResultValues">Result values.</param>
	        /// <param name="NbrofPointsPerType">Nbrof points per type.</param>
            public abstract void AdjustParametersofCellPoints(MultidimensionalArray ResultValues, Dictionary<int, int> NbrofPointsPerType);

	        /// <summary>
	        /// This function is used together with AdjustParametersofCellPoints to evaluate all points of a cell in one step
	        /// It returns the Coordinates of all Points in an Array.
	        /// </summary>
	        /// <param name="PointsCoordinates">Points coordinates.</param>
	        /// <param name="NbrofPointsPerType">Nbrof points per type.</param>
            public void PointsCoordinateArray(out MultidimensionalArray PointsCoordinates,out Dictionary<int,int> NbrofPointsPerType)
            {
                int SpatialDimension=0;
                int TotalNbrOfPoints = 0;
                NbrofPointsPerType = new Dictionary<int, int>();
                foreach (KeyValuePair<int, List<PointType>> DictEntry in ItemList)
                {
                    TotalNbrOfPoints += DictEntry.Value.Count;
                    NbrofPointsPerType.Add(DictEntry.Key, DictEntry.Value.Count);
                    if(DictEntry.Value.Count != 0)
                        SpatialDimension = DictEntry.Value[0].Coordinates.GetLength(1);
                }

                PointsCoordinates = MultidimensionalArray.Create(TotalNbrOfPoints,SpatialDimension);
                int indx = 0;
                foreach (KeyValuePair<int, List<PointType>> DictEntry in ItemList)
                {
                    for (int k = 0; k < NbrofPointsPerType[DictEntry.Key]; k++)
                    {
                        for (int dim = 0; dim < SpatialDimension; dim++)
                        {
                            PointsCoordinates[indx, dim] = DictEntry.Value[k].Coordinates[0, dim];
                        }
                        indx++;
                    }
                }
            }

	        /// <summary>
	        /// This function is used in combination with AdjustActive to check for possible merging of two level sets
	        /// For each Point two extension points are created along an extension vector. The extension vector is computed 
            /// by using the velocity vector of the last adevection step but limiting it if it differ for more than an angle of
            /// 45 degrees from the normal vector of the point
	        /// This function returns an array of all the extension points
	        /// </summary>
	        /// <param name="PointsCoordinates">Points coordinates.</param>
	        /// <param name="NbrofPointsPerType">Nbrof points per type.</param>
	        public void MergingTestCoordinatesArray(out MultidimensionalArray PointsCoordinates, out Dictionary<int, int> NbrofPointsPerType)
            {
                int SpatialDimension = 0;
                int TotalNbrOfPoints = 0;
                NbrofPointsPerType = new Dictionary<int, int>();
                foreach (KeyValuePair<int, List<PointType>> DictEntry in ItemList)
                {
                    TotalNbrOfPoints += DictEntry.Value.Count;
                    NbrofPointsPerType.Add(DictEntry.Key, DictEntry.Value.Count);
                    if (DictEntry.Value.Count != 0)
                        SpatialDimension = DictEntry.Value[0].Coordinates.GetLength(1);
                }
                PointsCoordinates = MultidimensionalArray.Create(2*TotalNbrOfPoints, SpatialDimension);
                int indx = 0;
                MultidimensionalArray VelocityVector;
                MultidimensionalArray TestDirectionVector;
                foreach (KeyValuePair<int, List<PointType>> DictEntry in ItemList)
                {
                    for (int k = 0; k < NbrofPointsPerType[DictEntry.Key]; k++)
                    {
                        TestDirectionVector = DictEntry.Value[k].Normal;
                        
                        if (DictEntry.Value[k].Coordinates == DictEntry.Value[k].Coordinates_prev)
                        {
                            TestDirectionVector = DictEntry.Value[k].Normal;
                        }
                        else
                        {
                            VelocityVector = ComputeVectorbetweenPoints(DictEntry.Value[k].Coordinates, DictEntry.Value[k].Coordinates_prev);

                            double Angle = VectorAngle(VelocityVector, DictEntry.Value[k].Normal);
                            if (Angle > (Math.PI / 4) && Angle < (3 * (Math.PI / 4)))
                            {
                                MultidimensionalArray CoordVecX = MultidimensionalArray.Create(1, 3);
                                MultidimensionalArray CoordVecY = MultidimensionalArray.Create(1, 3);
                                for (int dim = 0; dim < SpatialDimension; dim++)
                                {
                                    CoordVecX[0, dim] = DictEntry.Value[k].Normal[0, dim];
                                    CoordVecY[0, dim] = VelocityVector[0, dim];
                                }
                                CoordVecY = NormalizeVector(CoordVecY);
                                MultidimensionalArray CoordVecZ = NormalizeVector(CrossProduct(CoordVecX, CoordVecY));
                                CoordVecY = NormalizeVector(CrossProduct(CoordVecX, CoordVecZ));

                                MultidimensionalArray[] TestDirection = { MultidimensionalArray.Create(1, 3), MultidimensionalArray.Create(1, 3) };
                                TestDirection[0][0, 0] = CoordVecX[0, 0] * Math.Cos(0.25 * Math.PI) + CoordVecY[0, 0] * Math.Sin(0.25 * Math.PI);
                                TestDirection[0][0, 1] = CoordVecX[0, 1] * Math.Cos(0.25 * Math.PI) + CoordVecY[0, 1] * Math.Sin(0.25 * Math.PI);
                                TestDirection[0][0, 2] = CoordVecX[0, 2] * Math.Cos(0.25 * Math.PI) + CoordVecY[0, 2] * Math.Sin(0.25 * Math.PI);

                                TestDirection[1][0, 0] = CoordVecX[0, 0] * Math.Cos(-0.25 * Math.PI) + CoordVecY[0, 0] * Math.Sin(-0.25 * Math.PI);
                                TestDirection[1][0, 1] = CoordVecX[0, 1] * Math.Cos(-0.25 * Math.PI) + CoordVecY[0, 1] * Math.Sin(-0.25 * Math.PI);
                                TestDirection[1][0, 2] = CoordVecX[0, 2] * Math.Cos(-0.25 * Math.PI) + CoordVecY[0, 2] * Math.Sin(-0.25 * Math.PI);

                                MultidimensionalArray[] _TestDirection = { MultidimensionalArray.Create(1, SpatialDimension), MultidimensionalArray.Create(1, SpatialDimension) };
                                for (int dim = 0; dim < SpatialDimension; dim++)
                                {
                                    _TestDirection[0][0, dim] = TestDirection[0][0, dim];
                                    _TestDirection[1][0, dim] = TestDirection[1][0, dim];
                                }

                                double[] TestVectorAngle = new double[2];
                                TestVectorAngle[0] = VectorAngle(VelocityVector, _TestDirection[0]);
                                TestVectorAngle[1] = VectorAngle(VelocityVector, _TestDirection[1]);
                                if (TestVectorAngle[0] < TestVectorAngle[1])
                                {
                                    TestDirectionVector = _TestDirection[0];
                                }
                                else if (TestVectorAngle[0] > TestVectorAngle[1])
                                {
                                    TestDirectionVector = _TestDirection[1];
                                }
                                else
                                    throw new Exception("Failure in 45°Vector Computation");
                            }
                            else
                            {
                                TestDirectionVector = NormalizeVector (VelocityVector);
                            }
                        }
                        
                        for (int dim = 0; dim < SpatialDimension; dim++)
                        {
                            PointsCoordinates[indx, dim] = DictEntry.Value[k].Coordinates[0, dim] -
                                TestDirectionVector[0,dim] * 2 * DictEntry.Value[k].AdvectionStepLength;
                        }
                        indx++;
                        for (int dim = 0; dim < SpatialDimension; dim++)
                        {
                            PointsCoordinates[indx, dim] = DictEntry.Value[k].Coordinates[0, dim] +
                                TestDirectionVector [0, dim] * 2 * DictEntry.Value[k].AdvectionStepLength;
                        }
                        indx++;
                    }
                }
            }

	        /// <summary>
	        /// This function is used in combination with MergingTestCoordinatesArray to check for possible merging of two level sets
	        /// Uses the Results of the extension points and deactivates a point if the sign of the level set at the two extension points is equal
	        /// </summary>
	        /// <param name="ResultValues">Result values.</param>
	        /// <param name="NbrofPointsPerType">Nbrof points per type.</param>
            public void AdjustActive(MultidimensionalArray ResultValues, Dictionary<int, int> NbrofPointsPerType)
            {
                int SpatialDimension = ResultValues.GetLength(0) - 1;
                MultidimensionalArray LevelSet = MultidimensionalArray.Create(1, 1);
                MultidimensionalArray LevelSetGradient = MultidimensionalArray.Create(1, SpatialDimension);
                int index = 0;

                foreach (KeyValuePair<int, int> DictEntry in NbrofPointsPerType)
                {
                    for (int k = 0; k < DictEntry.Value; k++)
                    {
                        if (ItemList[DictEntry.Key][k].Active)
                        {
                            if (Math.Sign(ResultValues[0, index]) == Math.Sign(ResultValues[0, index + 1]))
                            {
                                Console.WriteLine("Deactivated");
                                ItemList[DictEntry.Key][k].Active = false;
                                ItemList[DictEntry.Key][k].DeactivationSign = Math.Sign(ResultValues[0, index]);
                            }
                            index += 2;
                        }
                    }
                }
            }

	        /// <summary>
	        /// Sets AllPointsInactive = true if all the points in a cell are inactive
	        /// </summary>
            public void CheckforInactivePoints()
            {
                bool AllPointsInactive = true;
                foreach(KeyValuePair<int,List<PointType>> DictEntry in ItemList)
                {
                    foreach(PointType Point in DictEntry.Value)
                    {
                        if (AllPointsInactive)
                        {
                            AllPointsInactive = (Point.Active) ? false : true;
                        }
                        else
                            break;
                    }
                    if (!AllPointsInactive)
                        break;
                }
                if (AllPointsInactive)
                    this.AllPointsInactive = AllPointsInactive;
            }

	        /// <summary>
	        /// Initializes a new instance of the
	        /// <see cref="T:SemiLagrangianLevelSet.SemiLagrangianLevelSetMethod`3.ItemsofCell"/> class.
	        /// </summary>
	        /// <param name="CellIndex">Cell index.</param>
            protected ItemsofCell (int CellIndex)
            {
                this.CellIndex = CellIndex;
            }

	        /// <summary>
	        /// Gets the number of points in a cell.
	        /// </summary>
	        /// <returns>The number of pointsin cell.</returns>
            public int GetNumberOfPointsinCell()
            {
                int PointCount = 0;
                foreach (KeyValuePair<int, List<PointType>> DictEntry in ItemList)
                    PointCount += DictEntry.Value.Count;
                return PointCount;
            }

	        /// <summary>
	        /// Prints the point data of a cell.
	        /// </summary>
	        /// <returns>The to string.</returns>
            public abstract string PrintToString();

	        /// <summary>
	        /// Prints the point data to a csv file.
	        /// </summary>
	        /// <returns>The to string.</returns>
	        public abstract string PrintToCSV(out string header);

        }
        protected abstract PointsofCellType NewPointsofCellObject(int CellIndex,IGridData Grid,MinimalDistanceSearchMode LevelSetCorrection);

	    /// <summary>
	    /// Basic class for a point.
	    /// It stores the coordinates, the previous coordinates, the normal vector based on the level set,
    	/// the level set value and the length of the last advection step.
	    /// Additionally the activation boolean and the sign of a deactivation is stored.
	    /// </summary>
        [Serializable]
        public abstract class SinglePoint
        {
            public MultidimensionalArray Coordinates;
            public MultidimensionalArray Coordinates_prev;
            public MultidimensionalArray Normal;
            public MultidimensionalArray Phi;
            public double AdvectionStepLength;
            public bool Active;
            public int DeactivationSign = 0;

	        /// <summary>
	        /// Initializes a new instance of the
	        /// <see cref="T:SemiLagrangianLevelSet.SemiLagrangianLevelSetMethod`3.SinglePoint"/> class.
	        /// with initial coordinates and a initial advectionstepLength
	        /// </summary>
	        /// <param name="Coordinates">Coordinates.</param>
	        /// <param name="AdvectionStepLengthINI">Advection step length ini.</param>
            protected SinglePoint(MultidimensionalArray Coordinates, double AdvectionStepLengthINI)
            {
                this.Coordinates = Coordinates;
                Coordinates_prev = Coordinates;
                Normal = null;
                Phi = null;
                Active = true;
                AdvectionStepLength = AdvectionStepLengthINI;
            }

	        /// <summary>
	        /// Takes the LevelSet Value and the LevelSet Gradient Values and 
	        /// adjusts the points parameters.
	        /// Depending on the settings the normal vector adjustment is dampened.
	        /// </summary>
	        /// <param name="LevelSet">Level set.</param>
	        /// <param name="LevelSetGradient">Level set gradient.</param>
            public virtual void AdjustParameters(MultidimensionalArray LevelSet, MultidimensionalArray LevelSetGradient)
            {
                Phi = LevelSet;
                if (Normal != null)
                {
                    if (NormalVectorDamping == NormalVectorDampingMode.Off)
                    {
                        Normal = NormalizeVector(LevelSetGradient);
                    }
                    else
                    {
                        MultidimensionalArray NewNormal = NormalizeVector(LevelSetGradient);
                        if (NormalVectorDamping == NormalVectorDampingMode.MaxAngle90)
                        {
                            if (DotProduct(NewNormal, Normal) >= 0)
                            {
                                Normal = NewNormal;
                            }
                        }
                        else if (NormalVectorDamping == NormalVectorDampingMode.MaxAngle45)
                        {
                            if (VectorAngle(NewNormal, Normal) <= 0.25 * Math.PI)
                            {
                                Normal = NewNormal;
                            }
                        }
                        else if (NormalVectorDamping == NormalVectorDampingMode.MaxAngle45_avg)
                        {
                            if (DotProduct(NewNormal, Normal) >= 0 && VectorAngle(NewNormal, Normal) >= 0.25 * Math.PI)
                            {
                                Normal = NormalizeVector(VectorAddition(1.0, NewNormal, 1.0, Normal));
                            }
                        }
                    }
                }
                else
                {
                    Normal = NormalizeVector(LevelSetGradient);
                }
            }

            /// <summary>
            /// Prints Point properties to a string
            /// </summary>
            /// <returns>The to string.</returns>
            public string PrintToString()
            {
                string line = null;
                string[] coord = new string[3] { "x", "y", "z" };
                for(int dim=0;dim<Coordinates.GetLength(1);dim++)
                {
                    line += "|" + coord[dim] + string.Format(":{0,8:0.00000}",Coordinates[0, dim]);
                }
                line += "| ";
                string[] normal = new string[3] { "n_x", "n_y", "n_z" };
                for (int dim = 0; dim < Coordinates.GetLength(1); dim++)
                {
                    line += "|" + normal[dim] + string.Format(":{0,8:0.00000}", Normal[0, dim]);
                }
                line += "| ";
                line += "|phi:"+string.Format("{0,8:0.00000}", Phi[0, 0])+"| ";
                return line;
            }

            /// <summary>
            /// Prints Point properties to a csv file.
            /// </summary>
            /// <returns>The to csv.</returns>
            /// <param name="header">Header.</param>
            public virtual string PrintToCSV(out string header)
            {
                string line = null;
                for (int dim = 0; dim < Coordinates.GetLength(1); dim++)
                {
                    line += Coordinates[0, dim] + " , ";
                }
                if (Phi != null) line += Phi[0, 0] + " , ";
                else line += " NaN, ";
                string[] coord = new string[3] { "x", "y", "z" };
                header = null;
                for (int dim = 0; dim < Coordinates.GetLength(1); dim++)
                {
                    header += coord[dim] + ",";
                }
                header += "phi,";
                line += " ";
                return line;
            }
        }

        /// <summary>
        /// Basic class for a k-d-Tree to find the nearest neighbour of a point
        /// </summary>
        public abstract class K_D_TreeFilter
        {
            protected readonly int SpatialDimension;
            private readonly int NbrOfPointsinNode;
            private enum RelativeNodePosition { InsideStrictNodeBound, InsideNodeBound, OutsideNodeBound };
            private readonly Node Root;

            /// <summary>
            /// Initializes a new instance of the
            /// <see cref="T:SemiLagrangianLevelSet.SemiLagrangianLevelSetMethod`3.K_D_TreeFilter"/> class.
            /// Takes the list of points and the maximum number of points at final node.
            /// </summary>
            /// <param name="Points">Points.</param>
            /// <param name="NbrOfPointsinNode">Nbr of pointsin node.</param>
            public K_D_TreeFilter(List<PointType> Points,int NbrOfPointsinNode)
            {
                SpatialDimension = Points[0].Coordinates.GetLength(1);
                this.NbrOfPointsinNode = NbrOfPointsinNode;

                Root = new Node (null,Points);
                if (Points.Count > NbrOfPointsinNode)
                {
                    BuildTree(Root);
                }
            }

            /// <summary>
            /// Class for a Node in the tree
            /// It stores a pointer to the parent node and the level in the tree.
            /// The DiscriminationKey is the dimension based on which the points are seperated with the DiscriminationValue as the decision value
            /// BoxBoundary stores the min and max Limit of the Points in the node and all its child nodes for all dimensions. The zero index represents the min
            /// value while the one index stores the max value
            /// BoxBoundaryFree works like BoxBoundary while storing bools determining if there is another node in the tree on the other side of a boundary. If thats the case 
            /// the boolean is false.
            /// Child_Lower and Child_Upper store the child nodes.
            /// </summary>
            private class Node
            {
                public Node Parent = null;
                public int TreeLevel;
                public int DiscriminationKey;
                public List<PointType> Points = null;

                public double DiscriminationValue;
                public double[][] BoxBoundary = new double[2][];
                public bool[][] BoxBoundaryFree = new bool[2][];

                public Node Child_Lower = null;
                public Node Child_Upper = null;

                /// <summary>
                /// Initializes a new instance of the
                /// <see cref="T:SemiLagrangianLevelSet.SemiLagrangianLevelSetMethod`3.K_D_TreeFilter.Node"/> class.
                /// Builds a new node taking a list of points as parameters
                /// </summary>
                /// <param name="Parent">Parent.</param>
                /// <param name="Points">Points.</param>
                public Node(Node Parent,List<PointType> Points)
                {
                    int SpatialDimension = Points[0].Coordinates.GetLength(1);
                    this.Parent = Parent;
                    this.Points = Points;
                    if (Parent == null)
                    {
                        for (int k = 0; k < 2; k++)
                        {
                            BoxBoundary[k] = new double[SpatialDimension];
                            BoxBoundaryFree[k] = new bool[SpatialDimension];
                        }
                        for (int dim = 0; dim < SpatialDimension; dim++)
                        {
                            BoxBoundaryFree[0][dim] = true;
                            BoxBoundaryFree[1][dim] = true;
                            Points.Sort(new PointComparer(dim));
                            BoxBoundary[0][dim] = Points.First().Coordinates[0, dim];
                            BoxBoundary[1][dim] = Points.Last().Coordinates[0, dim];
                        }
                        TreeLevel = 0;
                    }
                    else
                    {
                        for(int k=0;k<2;k++)
                        {
                            BoxBoundary[k] = new double[SpatialDimension];
                            BoxBoundaryFree[k] = new bool[SpatialDimension];
                            for(int dim=0;dim<SpatialDimension;dim++)
                            {
                                BoxBoundary[k][dim] = Parent.BoxBoundary[k][dim];
                                BoxBoundaryFree[k][dim] = Parent.BoxBoundaryFree[k][dim];
                            }
                        }
                        TreeLevel = Parent.TreeLevel + 1;
                        BoxBoundary[0][Parent.DiscriminationKey] = this.Points.First().Coordinates[0, Parent.DiscriminationKey];
                        BoxBoundary[1][Parent.DiscriminationKey] = this.Points.Last().Coordinates[0, Parent.DiscriminationKey];

                        if (BoxBoundary[0][Parent.DiscriminationKey] < Parent.DiscriminationValue &&
                            BoxBoundary[1][Parent.DiscriminationKey] < Parent.DiscriminationValue)
                            BoxBoundaryFree[1][Parent.DiscriminationKey] = false;
                        else 
                        if (BoxBoundary[0][Parent.DiscriminationKey] > Parent.DiscriminationValue &&
                            BoxBoundary[1][Parent.DiscriminationKey] > Parent.DiscriminationValue)
                            BoxBoundaryFree[0][Parent.DiscriminationKey] = false;
                        else
                            throw new Exception("Node does overlap with Parents DiscriminationValue");
                    }

                    DiscriminationKey = TreeLevel % SpatialDimension;
                    this.Points.Sort(new PointComparer(DiscriminationKey));
                    int TotalNbr = this.Points.Count;
                    int MedianNbr = this.Points.Count / 2;
                    if (MedianNbr != 0)
                    {
                        DiscriminationValue = (this.Points[MedianNbr - 1].Coordinates[0, DiscriminationKey] +
                                               this.Points[MedianNbr].Coordinates[0, DiscriminationKey]) / 2;
                    }
                    else
                    {
                        DiscriminationValue = this.Points[MedianNbr].Coordinates[0, DiscriminationKey];
                    }
                }
            }

            /// <summary>
            /// Takes provisionally final node and builds new nodes if the number of points in the node is higher than the limit
            /// </summary>
            /// <param name="RootNode">Root node.</param>
            private void BuildTree(Node RootNode)
            {
                RootNode.Points.Sort (new PointComparer (RootNode.DiscriminationKey));
                int TotalNbr = RootNode.Points.Count;
                int MedianNbr = RootNode.Points.Count / 2;

                int ChildLevel = RootNode.TreeLevel + 1;
                List<PointType> LowerList = RootNode.Points.GetRange(0, MedianNbr).ToList();
                List<PointType> UpperList = RootNode.Points.GetRange(MedianNbr, TotalNbr - MedianNbr).ToList();
                if (LowerList.Count <= 1 || UpperList.Count <= 1)
                    Console.WriteLine("Ende");
                RootNode.Child_Lower = new Node (RootNode, LowerList);
                RootNode.Child_Upper = new Node (RootNode, UpperList);

                if (LowerList.Count > NbrOfPointsinNode && UpperList.Count > NbrOfPointsinNode)
                {
                    BuildTree (RootNode.Child_Upper);
                    BuildTree (RootNode.Child_Lower);
                }
                RootNode.Points = null;
            }

            /// <summary>
            /// Builds a method to compare two points based on a given DiscriminationKey
            /// </summary>
            private class PointComparer : IComparer<PointType>
            {
                readonly int DiscriminationKey;
                public PointComparer(int DiscriminationKey)
                {
                    this.DiscriminationKey = DiscriminationKey;
                }
                public int Compare (PointType x, PointType y)
                {
                    double KV_x, KV_y;
                    KV_x = x.Coordinates [0, DiscriminationKey];
                    KV_y = y.Coordinates [0, DiscriminationKey];
                    if (KV_x < KV_y)
                        return -1;
                    else if (KV_x > KV_y)
                        return 1;
                    else
                        return 0;
                }
            }

            /// <summary>
            /// Takes a node, a point and a radius
            /// Checks the relation of the radius around a node with the boxboundary of a node
            /// </summary>
            /// <param name="Node"></param>
            /// <param name="GlobalPoint"></param>
            /// <param name="Radius"></param>
            /// <returns>
            ///  2: Circle with Radius around Point is fully confined by Box
            ///  1: Circle with Radius around Point is fully confined by Box with openBorder
            ///  0: Circle with Radius around Point does overlap with Box
            /// -1: Circle with Radius around Point does not overlap with Box 
            /// </returns>
            private int IsInsideNode(Node Node, MultidimensionalArray GlobalPoint, double Radius)
            {
                RelativeNodePosition GlobalNodePos;
                GlobalNodePos = RelativeNodePosition.InsideStrictNodeBound;
                for(int dim=0;dim<SpatialDimension;dim++)
                {
                    if (!(Node.BoxBoundary[0][dim] < GlobalPoint[dim] && GlobalPoint[dim] < Node.BoxBoundary[1][dim]))
                    {
                        if(Node.BoxBoundary[0][dim] > GlobalPoint[dim] && Node.BoxBoundaryFree[0][dim] == true)
                        {
                            GlobalNodePos = (GlobalNodePos != RelativeNodePosition.OutsideNodeBound) ?
                                RelativeNodePosition.InsideNodeBound : RelativeNodePosition.OutsideNodeBound;
                        }
                        else if (Node.BoxBoundary[1][dim] < GlobalPoint[dim] && Node.BoxBoundaryFree[1][dim] == true)
                        {
                            GlobalNodePos = (GlobalNodePos != RelativeNodePosition.OutsideNodeBound) ?
                                RelativeNodePosition.InsideNodeBound : RelativeNodePosition.OutsideNodeBound;
                        }
                        else
                        {
                            GlobalNodePos = RelativeNodePosition.OutsideNodeBound;
                        }
                    }
                }

                if (GlobalNodePos == RelativeNodePosition.InsideStrictNodeBound)
                {
                    int CirclePos = 2;
                    for (int dim = 0; dim < SpatialDimension; dim++)
                    {
                        for (int i = 0; i < 2; i++)
                        {
                            if (Math.Abs(Node.BoxBoundary[i][dim] - GlobalPoint[dim]) < Radius)
                            {
                                if(Node.BoxBoundaryFree[i][dim])
                                {
                                    CirclePos = (CirclePos != 0) ? 1 : CirclePos;
                                }
                                else
                                {
                                    CirclePos = 0;
                                }
                            }
                        }
                    }
                    return CirclePos;
                }
                else if (GlobalNodePos == RelativeNodePosition.InsideNodeBound)
                {
                    for (int dim = 0; dim < SpatialDimension; dim++)
                    {
                        for (int i = 0; i < 2; i++)
                        {

                            if ((Math.Abs(Node.BoxBoundary[i][dim] - GlobalPoint[dim]) < Radius) && !Node.BoxBoundaryFree[i][dim])
                                    return 0;
                        }
                    }
                    return 1;
                }
                else if (GlobalNodePos == RelativeNodePosition.OutsideNodeBound)
                {
                    bool OutsideBoxPlusRadius = false;
                    for(int dim=0;dim<SpatialDimension;dim++)
                    {
                        if ((GlobalPoint[dim] <= Node.BoxBoundary[0][dim] - Radius) || (Node.BoxBoundary[1][dim] + Radius <= GlobalPoint[dim]))
                            OutsideBoxPlusRadius = true;
                    }
                    if (OutsideBoxPlusRadius)
                        return -1;

                    int[] CornerIdent = new int[SpatialDimension];
                    for (int dim = 0; dim < SpatialDimension; dim++)
                    {
                        if (GlobalPoint[dim] <= Node.BoxBoundary[0][dim])
                            CornerIdent[dim] = 0;
                        else if (Node.BoxBoundary[1][dim] <= GlobalPoint[dim])
                            CornerIdent[dim] = 1;
                        else
                            return 0;
                    }

                    double[] Corner = new double[SpatialDimension];
                    for (int dim = 0; dim < SpatialDimension; dim++)
                    {
                        Corner[dim] = Node.BoxBoundary[CornerIdent[dim]][dim];
                    }

                    double dist = 0;
                    for (int dim = 0; dim < SpatialDimension; dim++)
                    {
                        dist += (Corner[dim] - GlobalPoint[dim]).Pow2();
                    }
                    dist = Math.Sqrt(dist);

                    if (dist >= Radius)
                        return -1;

                    return 0;
                }
                else
                    throw new Exception("GlobalNodePos was set to value not implemented");
            }

            /// <summary>
            /// Takes CoordinatePoint as Parameter and returns the nearest Neigbour in the k-d-Tree
            /// </summary>
            /// <returns>The nearest neighbour.</returns>
            /// <param name="GlobalNode">Global node.</param>
            public PointType FindNearestNeighbour(MultidimensionalArray GlobalNode)
            {
                MultidimensionalArray GlobalNode_cpy;
                GlobalNode_cpy = GlobalNode.CloneAs();

                NearestPointInTreeArea(GlobalNode, out PointType NearestPoint, out double MinRadius,out Node TerminalNode);

                GlobalNode = GlobalNode_cpy;
                Node CurrentNode = TerminalNode;
                bool RadiusIsInside = false;
                
                while (CurrentNode.Parent != null && RadiusIsInside == false)
                {
                    Node OriginNode = CurrentNode;
                    CurrentNode = CurrentNode.Parent;
                    Node OppositNode;
                    if (CurrentNode.Child_Lower == OriginNode)
                        OppositNode = CurrentNode.Child_Upper;
                    else
                        OppositNode = CurrentNode.Child_Lower;

                    int InsideNode = IsInsideNode(CurrentNode, GlobalNode, MinRadius);
                    if (InsideNode >= 1)
                        RadiusIsInside = true;

                    InsideNode = IsInsideNode(OppositNode, GlobalNode, MinRadius);
                    if (InsideNode >= 0)
                        MinRadiusDownwardSearch(GlobalNode, ref NearestPoint, ref MinRadius, OppositNode);
                }
                MultidimensionalArray NormaltoGlobalNode = MultidimensionalArray.Create(1, SpatialDimension);
                for (int dim = 0; dim < SpatialDimension; dim++)
                    NormaltoGlobalNode[0, dim] = GlobalNode[dim];
                MultidimensionalArray NormaltoNearestPoint = NearestPoint.Coordinates.CloneAs();

                double FoundDistance = VectorNorm(ComputeVectorbetweenPoints(NormaltoGlobalNode, NormaltoNearestPoint));
                NormaltoGlobalNode = NormalizeVector(NormaltoGlobalNode);
                for (int dim = 0; dim < SpatialDimension; dim++)
                    NormaltoGlobalNode[0, dim]*=0.7;

                double LevelSetDistance = Math.Sqrt(GlobalNode[0].Pow2() + GlobalNode[1].Pow2()) - 0.7;

                return NearestPoint;
            }

            /// <summary>
            /// Based on a Radius around a CoordinatePoint this Function moves down the Tree from a specific node searching for a point inside the radius
            /// </summary>
            /// <param name="GlobalNode">Global node.</param>
            /// <param name="NearestPoint">Nearest point.</param>
            /// <param name="MinRadius">Minimum radius.</param>
            /// <param name="CurrentNode">Current node.</param>
            private void MinRadiusDownwardSearch(MultidimensionalArray GlobalNode, ref PointType NearestPoint,ref double MinRadius,Node CurrentNode)
            {
                if(CurrentNode.Points == null)
                {
                    if (IsInsideNode(CurrentNode.Child_Upper, GlobalNode, MinRadius) >= 0)
                        MinRadiusDownwardSearch(GlobalNode, ref NearestPoint, ref MinRadius, CurrentNode.Child_Upper);
                    if (IsInsideNode(CurrentNode.Child_Lower, GlobalNode, MinRadius) >= 0)
                        MinRadiusDownwardSearch(GlobalNode, ref NearestPoint, ref MinRadius, CurrentNode.Child_Lower);
                }
                else
                    NearestPointAndMinRadiusFromList(CurrentNode.Points, GlobalNode, ref NearestPoint, ref MinRadius);                
            }

            /// <summary>
            /// Finds the nearest neighbour in a k-d Tree.
            /// </summary>
            /// <param name="GlobalNode">Global node.</param>
            /// <param name="NearestPoint">Nearest point.</param>
            /// <param name="MinRadius">Minimum radius.</param>
            /// <param name="TerminalNode">Terminal node.</param>
            private void NearestPointInTreeArea(MultidimensionalArray GlobalNode,out PointType NearestPoint,out double MinRadius,out Node TerminalNode)
            {
                Node CurrentNode = Root;
                while (CurrentNode.Points == null)
                {
                    if (GlobalNode[CurrentNode.DiscriminationKey] <= CurrentNode.DiscriminationValue)
                        CurrentNode = CurrentNode.Child_Lower;
                    else
                        CurrentNode = CurrentNode.Child_Upper;
                }

                NearestPoint = null;
                MinRadius = double.MaxValue;
                NearestPointAndMinRadiusFromList(CurrentNode.Points, GlobalNode, ref NearestPoint, ref MinRadius);
                TerminalNode = CurrentNode;
            }

            /// <summary>
            /// Abstract fucntion to compute the nearest neighbour from a list of points
            /// </summary>
            /// <param name="Points">Points.</param>
            /// <param name="GlobalNode">Global node.</param>
            /// <param name="NearestPoint">Nearest point.</param>
            /// <param name="MinRadius">Minimum radius.</param>
            protected abstract void NearestPointAndMinRadiusFromList(List<PointType> Points, MultidimensionalArray GlobalNode, ref PointType NearestPoint, ref double MinRadius);
        }

        public class OldNearestPointMemory
        {
            List<PointType> PointerReference;
            List<MultidimensionalArray> CoordinateMemory;
            Dictionary<int, Dictionary<int, int>> OldNearestNeighbours;

            public OldNearestPointMemory()
            {
                PointerReference = new List<PointType>();
                CoordinateMemory = new List<MultidimensionalArray>();
                OldNearestNeighbours = new Dictionary<int, Dictionary<int, int>>();
            }
            public void AddNode(int cell,int NodeNbr,PointType Point)
            {
                int index = PointerReference.IndexOf(Point);
                if (index == -1)
                {
                    PointerReference.Add(Point);
                    CoordinateMemory.Add(Point.Coordinates.CloneAs());
                    index = CoordinateMemory.Count;
                }
                if (OldNearestNeighbours.ContainsKey(cell))
                {
                    OldNearestNeighbours[cell].Add(NodeNbr, index);
                }
                else
                {
                    OldNearestNeighbours.Add(cell, new Dictionary<int, int>());
                    OldNearestNeighbours[cell].Add(NodeNbr, index);
                }
            }
            public MultidimensionalArray FindOldNearestNeighbour(int cell,int NodeNbr)
            {
                return CoordinateMemory[OldNearestNeighbours[cell][NodeNbr]];
            }       
            public PointType FindOldNearestNeighbourPoint(int cell,int NodeNbr)
            {
                return PointerReference[OldNearestNeighbours[cell][NodeNbr]];
            }
        }

        public abstract class CorrectionMask
        {
            public int CellIndex;
            public Dictionary<int, int> MovedPoints;

            public CorrectionMask(int CellIndex, int[] PointsSign)
            {
                this.CellIndex = CellIndex;
                MovedPoints = new Dictionary<int, int>();
                foreach (int i in PointsSign)
                    MovedPoints.Add(i, 0);
            }

            public bool Compare(CorrectionMask Other)
            {
                if (Other.CellIndex == this.CellIndex) return true;
                else return false;
            }

        }
        protected abstract CorrectionMaskType NewFullCorrectionMask(int CellIndex);
        protected abstract CorrectionMaskType NewCorrectionMask(int CellIndex);

        /// <summary>
        /// Stores all the Cell Objects that store the points
        /// </summary>
        protected List<PointsofCellType> Points;

        /// <summary>
        /// Stores the cell index as keys and the index in the Points List as the value
        /// </summary>
        protected Dictionary<int, int> CelltoArrayindex;

        public abstract void Initialize ();

        public void PerformTimestep(double dt, int substeps, int Timestep)
        {
            NarrowBandCells = InterfaceTracker.Regions.GetNearFieldMask(NarrowBandWidth);
            PointCount();
            InterfaceLevelSetGradient.Clear();
            InterfaceLevelSetGradient.Gradient(1, Interface, NarrowBandCells);
            //InterfaceTracker.UpdateTracker();
        }

        /// <summary>
        /// Setup function for cell objects and Dictionary
        /// </summary>
        /// <param name="PointMask">Point mask.</param>
        /// <param name="CellEval">If set to <c>true</c> cell eval.</param>
        protected void InitializePointstoLevelSet(out List<CorrectionMaskType> PointMask, ref bool CellEval)
        {
            NarrowBandCells = NarrowBandCells = InterfaceTracker.Regions.GetNearFieldMask(NarrowBandWidth);
            ReseedingBandCells = InterfaceTracker.Regions.GetCutCellMask();
            Points = new List<PointsofCellType>(0);
            CelltoArrayindex = new Dictionary<int, int>(0);
            Console.WriteLine(ReseedingBandCells.ItemEnum.Count() + " Cells in NarrowBandWidth:"+ ReseedingWdith + " of Interface");
            foreach (int CellIndex in ReseedingBandCells.ItemEnum)
            {
                Points.Add(NewPointsofCellObject(CellIndex, Grid, MinimalDistanceSearch));
                CelltoArrayindex[CellIndex] = Points.Count - 1;
            }
            ReseedingPoints(out PointMask);
            CellEval = false;
        }

        /// <summary>
        /// Reseeding function on all cut-cells. Function reseeds or removes points if the point number if higer or lower than the upper and lower limit.
        /// Function returns a Mask that contains information which cells were treated.
        /// </summary>
        /// <param name="PointMask">Point mask.</param>
        protected void ReseedingPoints(out List<CorrectionMaskType> PointMask)
        {
            ReseedingBandCells = InterfaceTracker.Regions.GetCutCellMask();
            foreach (int CellIndex in ReseedingBandCells.ItemEnum)
            {
                if (!CelltoArrayindex.ContainsKey(CellIndex))
                {
                    Points.Add(NewPointsofCellObject(CellIndex, Grid, MinimalDistanceSearch));
                    CelltoArrayindex[CellIndex] = Points.Count - 1;
                }
            }
            PointMask = new List<CorrectionMaskType>();
            foreach (PointsofCellType CellObj in Points)
            {
                int PointCount = CellObj.GetNumberOfPointsinCell();
                int LowerLimitPoints = (LowerLimitPointsPerCellPerDimension == 1) ? LowerLimitPointsPerCellPerDimension * SpatialDimension :
                    (int)(((double)(LowerLimitPointsPerCellPerDimension)).Pow(SpatialDimension));
                int UpperLimitPoints = (UpperLimitPointsPerCellPerDimension == 1) ? UpperLimitPointsPerCellPerDimension * SpatialDimension :
                    (int)(((double)(UpperLimitPointsPerCellPerDimension)).Pow(SpatialDimension));

                if(PointCount < LowerLimitPoints || PointCount > UpperLimitPoints)
                {
                    RemoveAllPoints(CellObj);
                    AddandAttractPoints(CellObj, PointMask, MovePointtoTarget);
                }
                /*
                PointCount = CellObj.GetNumberOfPointsinCell();

                if (PointCount < LowerLimitPoints)
                {
                    Console.WriteLine(PointCount + "<" + LowerLimitPoints);
                    AddandAttractPoints(CellObj, PointMask, MovePointtoTarget);
                }
                if (PointCount > UpperLimitPoints)
                {
                    Console.WriteLine(PointCount + ">" + UpperLimitPoints);
                    RemoveSurplusPoints(CellObj, PointMask);
                    CorrectionMaskType other = NewCorrectionMask(CellObj.CellIndex);
                    int index;
                    if ((index = PointMask.FindIndex(other.Compare)) != -1)
                    {
                        foreach(KeyValuePair<int,int> DictEntry in PointMask[index].MovedPoints)
                        PointMask[index].MovedPoints[1] = CellObj.ItemList[1].Count;
                        PointMask[index].MovedPoints[-1] = CellObj.ItemList[-1].Count;
                    }
                }
                */
            }
        }

        /// <summary>
        /// Function to reduce the number of points to a certain number.
        /// </summary>
        /// <param name="TreatedCell">Treated cell.</param>
        /// <param name="PointMask">Particle mask.</param>
        protected abstract void RemoveSurplusPoints(PointsofCellType TreatedCell, List<CorrectionMaskType> PointMask);

        /// <summary>
        /// Function to remove all points from an cell object
        /// </summary>
        /// <param name="TreatedCell"></param>
        protected void RemoveAllPoints(PointsofCellType TreatedCell)
        {
            foreach(KeyValuePair<int,List<PointType>> DictEntry in TreatedCell.ItemList)
            {
                DictEntry.Value.Clear();
            }
        }

        /// <summary>
        /// Function to add points to every cut-cell and attract them to their target position.
        /// </summary>
        /// <param name="TreatedCell">Treated cell.</param>
        /// <param name="PointMask">Point mask.</param>
        /// <param name="Attract">Attract.</param>
        public void AddandAttractPoints(PointsofCellType TreatedCell, List<CorrectionMaskType> PointMask, Func<PointsofCellType, int, MultidimensionalArray, List<CorrectionMaskType>, bool> Attract)
        {
            Dictionary<int, int> PointstoAdd = new Dictionary<int, int>();
            int PointNum = 0;
            foreach (KeyValuePair<int, List<PointType>> DictEntry in TreatedCell.ItemList.ToList())
            {
                PointstoAdd[DictEntry.Key] = 0;
                PointNum += DictEntry.Value.Count;
            }
            int TotalPointstoAdd = ((TargetNbrofPointsPerCellPerDimension == 1) ? TargetNbrofPointsPerCellPerDimension * SpatialDimension :
                    (int)(((double)(TargetNbrofPointsPerCellPerDimension)).Pow(SpatialDimension))) - PointNum;
            foreach(KeyValuePair<int,int> DictEntry in PointstoAdd.ToList())
            {
                PointstoAdd[DictEntry.Key] = TotalPointstoAdd / PointstoAdd.Count;
            }

            foreach (KeyValuePair<int, int> AddDict in PointstoAdd)
            {                
                MultidimensionalArray LocalPointIn = MultidimensionalArray.Create(AddDict.Value, SpatialDimension);
                if (AddDict.Value == 0)
                    continue;
                else if (AddDict.Value == SpatialDimension)
                {
                    double[] Coord1d = new double[AddDict.Value];
                    double interval = 2.0 / AddDict.Value;
                    for (int i = 0; i < AddDict.Value; i++)
                        Coord1d[i] = -1 + 0.5 * interval + interval * i;
                    for (int k = 0; k < LocalPointIn.GetLength(0); k++)
                    {
                        for (int dim = 0; dim < SpatialDimension; dim++)
                        {
                            LocalPointIn[k, dim] = Coord1d[k] + (CoordGen.NextDouble() * interval - 0.5 * interval);
                        }
                    }
                }
                else
                {                    
                    int AddNumberPerDim = (int)Math.Round(Math.Pow(AddDict.Value, 1.0/SpatialDimension));
                    double[] Coord1d = new double[AddNumberPerDim];
                    double interval = 2.0 / AddNumberPerDim;
                    for (int i = 0; i < AddNumberPerDim; i++)
                        Coord1d[i] = -1 + 0.5*interval + interval * i;

                    int index;
                    int stepsize;
                    for (int dim = 0; dim < SpatialDimension; dim++)
                    {
                        index = Coord1d.Length - 1;
                        stepsize = (dim == 0) ? 1 : Coord1d.Length * dim;
                        for (int k = 0; k < LocalPointIn.GetLength(0); k++)
                        {
                            if (k % stepsize == 0)
                                index = (index == Coord1d.Length - 1) ? 0 : index + 1;
                            LocalPointIn[k, dim] = Coord1d[index] + (CoordGen.NextDouble() * interval - 0.5 * interval);
                        }
                    }
                }
                int failCount = 0;
                MultidimensionalArray GlobalPointsOut = MultidimensionalArray.Create(LocalPointIn.GetLength(0), SpatialDimension);
                Grid.TransformLocal2Global(LocalPointIn, GlobalPointsOut, TreatedCell.CellIndex);
                for (int nbr = 0; nbr < GlobalPointsOut.GetLength(0); nbr++)
                {
                    MultidimensionalArray GlobalPointOut = MultidimensionalArray.Create(1, SpatialDimension);
                    for (int dim = 0; dim < SpatialDimension; dim++)
                        GlobalPointOut[0, dim] = GlobalPointsOut[nbr, dim];
                    if (!Attract(TreatedCell, AddDict.Key, GlobalPointOut, PointMask))
                    {
                        failCount++;
                    }
                    if (failCount > 2 * AddDict.Value)
                    {
                        bool CellPosSure = false;
                        MultidimensionalArray Phi = Evaluate(new DGField[] { Interface }, GlobalPointOut, TreatedCell.CellIndex, CellPosSure);
                        throw new Exception("Multiple failures in point attraction to level set! Something is wrong in Cell: " + TreatedCell.CellIndex +
                            "at Coord:[" + GlobalPointOut[0, 0] + "|" + GlobalPointOut[0, 1] + "] and Phi:" + Phi[0,0]);
                    }
                }
                //Console.WriteLine("Added " + GlobalPointsOut.GetLength(0) + " Points");
            }

        }

        /// <summary>
        /// Abstract function to define the attraction method of the points depending on the desired target position.
        /// </summary>
        /// <returns><c>true</c> if point was succesfully attracted and <c>false</c> otherwise.</returns>
        /// <param name="TreatedCell">Treated cell.</param>
        /// <param name="signKey">Sign key.</param>
        /// <param name="Coordinate">Coordinate.</param>
        /// <param name="ParticleMask">Particle mask.</param>
        protected abstract bool MovePointtoTarget(PointsofCellType TreatedCell, int signKey, MultidimensionalArray Coordinate, List<CorrectionMaskType> ParticleMask);

        /// <summary>
        /// Function to update the point parameters of every point and checking if points might be deactivated because of topology merging.
        /// </summary>
        /// <param name="timestepNo">Timestep no.</param>
        /// <param name="CellPosSure">If set to <c>true</c> cell position sure.</param>
        public void AdjustPointParamters(int timestepNo,bool CellPosSure)
        {
            DGField[] Fields = new DGField[1 + SpatialDimension];
            Fields[0] = Interface;
            for (int dim = 0; dim < SpatialDimension; dim++)
            {
                Fields[1 + dim] = InterfaceLevelSetGradient[dim];
            }
            MultidimensionalArray Results;

            for (int i = 0; i < Points.Count; i++)
            {
                if (!NarrowBandCells.Contains(Points[i].CellIndex) && !Points[i].AllPointsInactive)
                    throw new Exception("Advection step has moved Markers outside of the Narrow Band." +
                        " Consider increasing the NarrowBand Width or reducing the TimestepSize." +
                        " This Limitation is necessary to enable Topology change." +
                        " Either reduce the timestep size or enlarge the Narrrow Band");

                Points[i].PointsCoordinateArray(out MultidimensionalArray PointsCoordinates, out Dictionary<int, int> NbrofPointsPerType);
                Results = MultiEvaluatePointsInCell(Fields, PointsCoordinates, Points[i].CellIndex, CellPosSure);
                Points[i].AdjustParametersofCellPoints(Results, NbrofPointsPerType);
            }

            Dictionary<int, MultidimensionalArray []> NormalField = null;
            Dictionary<int, MultidimensionalArray> NormalInCell = null;

            if (TopologyMerging == TopologyMergingMode.On && timestepNo>0)
            {
                // Evaluates the normal vector of the level set function in a regular pattern in certain cells.
                // Averages the normal vector in a cell to one representative.
                // begin
                double PointsPerDimension = 1;
                int nbrofTestPoints = (int)(PointsPerDimension.Pow (SpatialDimension));
                double [,] nds = new double [nbrofTestPoints,SpatialDimension];

                double [] Coord1d = new double [(int)PointsPerDimension];
                double interval = 2.0 / (int)PointsPerDimension;
                for (int i = 0; i < PointsPerDimension; i++)
                    Coord1d [i] = -1 + 0.5 * interval + interval * i;

                int index;
                int stepsize;
                for (int dim = 0; dim < SpatialDimension; dim++)
                {
                    index = Coord1d.Length - 1;
                    stepsize = (dim == 0) ? 1 : Coord1d.Length * dim;
                    for (int k = 0; k < nbrofTestPoints; k++)
                    {
                        if (k % stepsize == 0)
                            index = (index == Coord1d.Length - 1) ? 0 : index + 1;
                        nds [k, dim] = Coord1d [index];
                    }
                }

                MultidimensionalArray[] result = new MultidimensionalArray [SpatialDimension];
                for (int dim = 0; dim < SpatialDimension; dim++)
                    result [dim] = MultidimensionalArray.Create (1,nbrofTestPoints);

                NormalField = new Dictionary<int, MultidimensionalArray []> ();

                foreach (int cell in NarrowBandCells.ItemEnum) {
                    NodeSet NormalEvalPoints = new NodeSet(Grid.iGeomCells.GetRefElement(cell), nds, false);
                    NormalField.Add(cell, new MultidimensionalArray[nbrofTestPoints]);
                    for (int dim = 0; dim < SpatialDimension; dim++) {
                        InterfaceLevelSetGradient[dim].Evaluate(cell, 1, NormalEvalPoints, result[dim]);
                    }
                    for (int k = 0; k < nbrofTestPoints; k++) {
                        NormalField[cell][k] = MultidimensionalArray.Create(1, SpatialDimension);
                        for (int dim = 0; dim < SpatialDimension; dim++) {
                            NormalField[cell][k][0, dim] = result[dim][0, k];
                        }
                        NormalField[cell][k] = NormalizeVector(NormalField[cell][k]);
                    }
                }

                NormalInCell = new Dictionary<int, MultidimensionalArray>();

                foreach(KeyValuePair<int,MultidimensionalArray[]> DictEntry in NormalField)
                {
                    MultidimensionalArray Normals = MultidimensionalArray.Create(1, SpatialDimension);
                    for(int nbr=0;nbr<DictEntry.Value.Count();nbr++)
                    {
                        for(int dim=0;dim<SpatialDimension;dim++)
                        {
                            Normals[0, dim] += (DictEntry.Value[nbr][0, dim]/DictEntry.Value.Count());
                        }
                    }
                    NormalInCell[DictEntry.Key] = Normals;
                }
                // end


                // Removes cell items with inactive points if their neighbour points are inactive too
                // begin
                foreach (PointsofCellType CellObj in Points)
                {
                    CellObj.CheckforInactivePoints();
                }                
                for(int i=0;i<Points.Count;i++)
                {
                    if (Points[i].AllPointsInactive)
                    {
                        bool AllNeighboursCellPointsInactive = true;
                        Grid.GetCellNeighbours(Points[i].CellIndex, GetCellNeighbours_Mode.ViaVertices, out int[] Neighbours, out int[] Empty);
                        foreach(int cell in Neighbours)
                        {
                            if (CelltoArrayindex.ContainsKey(cell) && !Points[CelltoArrayindex[cell]].AllPointsInactive)
                                AllNeighboursCellPointsInactive = false;
                        }
                        if(AllNeighboursCellPointsInactive)
                        {
                            List<int> Keys = CelltoArrayindex.Keys.ToList();
                            for (int k = 0; k < Keys.Count; k++)
                            {
                                if (CelltoArrayindex[Keys[k]] > i)
                                {
                                    CelltoArrayindex[Keys[k]] = CelltoArrayindex[Keys[k]] - 1;
                                }
                            }
                            CelltoArrayindex = new Dictionary<int, int>(CelltoArrayindex);
                            CelltoArrayindex.Remove(Points[i].CellIndex);
                            Points.RemoveAt(i--);
                        }
                    }
                }
                // end

                for (int i = 0; i < Points.Count; i++)
                {
                    // Removes cellitems from list if they are empty
                    //begin
                    bool EmptyCell = true;
                    foreach (KeyValuePair<int, List<PointType>> DictEntry in Points[i].ItemList)
                    {
                        if (!DictEntry.Value.IsNullOrEmpty())
                            EmptyCell = false;
                    }

                    if (EmptyCell)
                    {
                        List<int> Keys = CelltoArrayindex.Keys.ToList();
                        for (int k = 0; k < Keys.Count; k++)
                        {
                            if (CelltoArrayindex[Keys[k]] > i)
                            {
                                CelltoArrayindex[Keys[k]] = CelltoArrayindex[Keys[k]] - 1;
                            }
                        }
                        CelltoArrayindex = new Dictionary<int, int>(CelltoArrayindex);
                        CelltoArrayindex.Remove(Points[i].CellIndex);
                        Points.RemoveAt(i--);
                        continue;
                    }
                    //end

                    // Compare normal vectors in one cell if dotproduct < 0 -> TestForMerging True
                    bool TestForMerging = false;
                    int cell = Points[i].CellIndex;
                    for (int k0 = 0; k0 < NormalField[cell].Length - 1; k0++)
                    {
                        for (int k1 = k0 + 1; k1 < NormalField[cell].Length; k1++)
                        {
                            if (DotProduct(NormalField[cell][k0], NormalField[cell][k1]) < 0)
                            {
                                TestForMerging = true;
                                break;
                            }
                        }
                        if (TestForMerging)
                            break;
                    }

                    // Compare normal vectors neighboring cells if dotproduct < 0 -> TestForMerging True
                    if (!TestForMerging)
                    {
                        Grid.GetCellNeighbours(Points[i].CellIndex, GetCellNeighbours_Mode.ViaVertices, out int[] Neighbours, out int[] Connecting);
                        foreach (int j in Neighbours)
                        {
                            if (NarrowBandCells.Contains(j))
                            {
                                for (int k0 = 0; k0 < NormalField[cell].Length; k0++)
                                {
                                    for (int k1 = 0; k1 < NormalField[j].Length; k1++)
                                    {
                                        if (DotProduct(NormalField[cell][k0], NormalField[j][k1]) < 0)
                                        {
                                            MultidimensionalArray NormalCell, NormalNeigh;
                                            NormalCell = NormalField[cell][k0];
                                            NormalNeigh = NormalField[j][k1];
                                            BoundingBox BoxCell = new BoundingBox(SpatialDimension), BoxNeigh = new BoundingBox(SpatialDimension);
                                            Grid.iGeomCells.GetCellBoundingBox(cell, BoxCell);
                                            Grid.iGeomCells.GetCellBoundingBox(j, BoxNeigh);
                                            TestForMerging = true;
                                            Console.WriteLine("Merging at cell: " + cell);
                                            break;
                                        }
                                    }
                                    if (TestForMerging)
                                        break;
                                }
                                if (TestForMerging)
                                    break;
                            }
                        }
                    }
                    if (TestForMerging)
                    {
                        //Points[i].MergingTestCoordinatesArray(out MultidimensionalArray MergingTestCoordinates, out Dictionary<int, int> _NbrofPointsPerType);
                        //bool _CellPosSure = false;
                        //Results = MultiEvaluatePointsInCell(new DGField[] { Interface }, MergingTestCoordinates, Points[i].CellIndex, _CellPosSure);
                        //Points[i].AdjustActive(Results, _NbrofPointsPerType);
                    }
                }
            }
        }

        /// <summary>
        /// Function to evaluate an array of Fields at a coordinate.
        /// The function takes a cellindex guess and a boolean that determines if the cellindex is the correct one.
        /// If the cell index is surely the correct one the function is much faster
        /// </summary>
        /// <returns>A multidimensional array with the dimension[1,Fields.Length]</returns>
        /// <param name="Fields">Fields.</param>
        /// <param name="GlobalCoordinates">Global coordinates.</param>
        /// <param name="CellIndexGuess">Cell index guess.</param>
        /// <param name="CellPosSure">If set to <c>true</c> cell position sure.</param>
        protected MultidimensionalArray Evaluate(DGField[] Fields, MultidimensionalArray GlobalCoordinates, int CellIndexGuess, bool CellPosSure)
        {
            MultidimensionalArray FieldsResults = MultidimensionalArray.Create(1, Fields.Length);

            if (CellPosSure || ((GridData)Grid).Cells.IsInCell(GlobalCoordinates.ExtractSubArrayShallow(0, -1).To1DArray(), CellIndexGuess))
            {
                MultidimensionalArray logicCoords = MultidimensionalArray.Create(1, 1, SpatialDimension);
                Grid.TransformGlobal2Local(GlobalCoordinates, logicCoords, CellIndexGuess, 1, 0);
                MultidimensionalArray buffer = MultidimensionalArray.Create(1, 1);
                for (int f = 0; f < Fields.Length; f++)
                {
                    buffer[0, 0] = 0;
                    Fields[f].Evaluate(CellIndexGuess, 1, new NodeSet(((GridData)Grid).iGeomCells.GetRefElement(CellIndexGuess), logicCoords.ExtractSubArrayShallow(0, -1, -1), false), buffer, 1.0);
                    if (double.IsNaN(buffer[0, 0]) || double.IsInfinity(buffer[0, 0]))
                        throw new Exception("NaN of Infinity Values evaluated");
                    FieldsResults[0, f] = buffer[0, 0];
                }
                return FieldsResults;
            }

            Grid.GetCellNeighbours(CellIndexGuess, GetCellNeighbours_Mode.ViaVertices, out int[] NeighbourIndex, out int[] ConnectingEntities);
            foreach (int j in NeighbourIndex)
            {
                if (((GridData)Grid).Cells.IsInCell(GlobalCoordinates.ExtractSubArrayShallow(0, -1).To1DArray(), j))
                {
                    MultidimensionalArray logicCoords = MultidimensionalArray.Create(1, 1, SpatialDimension);
                    Grid.TransformGlobal2Local(GlobalCoordinates, logicCoords, j, 1, 0);
                    MultidimensionalArray buffer = MultidimensionalArray.Create(1, 1);
                    for (int f = 0; f < Fields.Length; f++)
                    {
                        buffer[0, 0] = 0;
                        Fields[f].Evaluate(j, 1, new NodeSet(((GridData)Grid).iGeomCells.GetRefElement(j), logicCoords.ExtractSubArrayShallow(0, -1, -1), false), buffer, 1.0);
                        if (double.IsNaN(buffer[0, 0]) || double.IsNaN(buffer[0, 0]))
                            throw new Exception("NaN of Infinity Values evaluated");
                        FieldsResults[0, f] = buffer[0, 0];
                    }
                    return FieldsResults;
                }
            }

            int NoOfLocallyUnlocatedPts = 1;
            int NoOfGloballyUnlocatedPts = NoOfLocallyUnlocatedPts.MPISum();
            while (NoOfGloballyUnlocatedPts > 0)
            {
                // do Comm


                NoOfGloballyUnlocatedPts = NoOfLocallyUnlocatedPts.MPISum();
            }
        

            Evaluator.Evaluate(1.0, Fields, GlobalCoordinates, 1.0, FieldsResults);
            return FieldsResults;
        }

        /// <summary>
        /// Function to evaluate an array of Fields at multiple points.
        /// The function takes a cellindex guess and a boolean that determines if the cellindex is the correct one for all points.
        /// </summary>
        /// <returns> A Multidimensional Array with the dimension[Fields.Length,Number of Points]</returns>
        /// <param name="Fields">Fields.</param>
        /// <param name="GlobalPointsCoordinates">Global points coordinates.</param>
        /// <param name="CellIndexGuess">Cell index guess.</param>
        /// <param name="CellPosSure">If set to <c>true</c> cell position sure.</param>
        protected MultidimensionalArray MultiEvaluatePointsInCell(DGField[] Fields, MultidimensionalArray GlobalPointsCoordinates, int CellIndexGuess, bool CellPosSure)
        {
            if(GlobalPointsCoordinates.GetLength(0) == 0)
            {
                return MultidimensionalArray.Create(Fields.Length, 0);
            }
            int NbrofPoints = GlobalPointsCoordinates.GetLength(0);

            MultidimensionalArray Result = MultidimensionalArray.Create(Fields.Length,1,NbrofPoints);

            if (CellPosSure)
            {
                MultidimensionalArray LogicPointsCoordinates = MultidimensionalArray.Create(1, NbrofPoints, SpatialDimension);
                Grid.TransformGlobal2Local(GlobalPointsCoordinates, LogicPointsCoordinates, CellIndexGuess, 1, 0);
                for (int f = 0; f < Fields.Length; f++)
                {
                    Fields[f].Evaluate(CellIndexGuess, 1, new NodeSet(((GridData)Grid).iGeomCells.GetRefElement(CellIndexGuess),
                        LogicPointsCoordinates.ExtractSubArrayShallow(0, -1, -1), false), Result.ExtractSubArrayShallow(f,-1, -1));
                }
            }
            else
            {
                BitArray InCellPointMask = new BitArray(NbrofPoints);
                int NbrInCell = 0;
                for (int index = 0; index < NbrofPoints; index++)
                {
                    if ((((GridData)Grid).Cells.IsInCell(GlobalPointsCoordinates.ExtractSubArrayShallow(index, -1).To1DArray(), CellIndexGuess)))
                    {
                        InCellPointMask[index] = true;
                        NbrInCell++;
                    }
                    else
                        InCellPointMask[index] = false;
                }
                MultidimensionalArray GlobalPointsInCellCoordinates = MultidimensionalArray.Create(NbrInCell, SpatialDimension);
                MultidimensionalArray ResultInCell = MultidimensionalArray.Create(Fields.Length,1, NbrInCell);
                int indexInCell = 0;
                for (int index = 0; index < NbrofPoints; index++)
                {
                    if (InCellPointMask[index] == true)
                    {
                        for (int dim = 0; dim < SpatialDimension; dim++)
                        {
                            GlobalPointsInCellCoordinates[indexInCell, dim] = GlobalPointsCoordinates[index, dim];
                        }
                        indexInCell++;
                    }
                    else if (InCellPointMask[index] == false)
                    {
                        bool PointInNeighbourhood = false;
                        Grid.GetCellNeighbours(CellIndexGuess, GetCellNeighbours_Mode.ViaVertices, out int[] NeighbourIndex, out int[] ConnectingEntities);
                        foreach (int j in NeighbourIndex)
                        {
                            if (((GridData)Grid).Cells.IsInCell(GlobalPointsCoordinates.ExtractSubArrayShallow(index, -1).To1DArray(), j))
                            {
                                PointInNeighbourhood = true;
                                MultidimensionalArray logicCoords = MultidimensionalArray.Create(1, 1, SpatialDimension);
                                MultidimensionalArray Coord = MultidimensionalArray.Create(1, SpatialDimension);
                                for (int dim = 0; dim < SpatialDimension; dim++) Coord[0, dim] = GlobalPointsCoordinates[index, dim];
                                Grid.TransformGlobal2Local(Coord, logicCoords, j, 1, 0);
                                MultidimensionalArray Res = MultidimensionalArray.Create(1, 1);
                                for (int f = 0; f < Fields.Length; f++)
                                {
                                    Res[0, 0] = 0;
                                    Fields[f].Evaluate(j, 1, new NodeSet(((GridData)Grid).iGeomCells.GetRefElement(j),
                                        logicCoords.ExtractSubArrayShallow(0, -1, -1), false), Res, 1.0);
                                    if (double.IsNaN(Res[0,0]) || double.IsInfinity(Res[0,0]))
                                        throw new Exception("NaN or Infinity not allowed");
                                    Result[f, 0, index] = Res[0, 0];
                                }
                            }
                        }
                        if (!PointInNeighbourhood)
                        {
                            MultidimensionalArray Res = MultidimensionalArray.Create(1, Fields.Length);
                            MultidimensionalArray Coord = MultidimensionalArray.Create(1, SpatialDimension);
                            for (int dim = 0; dim < SpatialDimension; dim++)
                                Coord[0, dim] = GlobalPointsCoordinates[index, dim];
                            Evaluator.Evaluate(1.0, Fields, Coord, 1.0,Res);
                            for (int f = 0; f < Fields.Length; f++)
                            {
                                if (double.IsNaN(Res[0,f]) || double.IsInfinity(Res[0,f]))
                                    throw new Exception("NaN or Infinity not allowed");
                                Result[f, 0, index] = Res[0, f];
                            }
                        }
                    }
                }
                if (NbrInCell > 0)
                {
                    MultidimensionalArray LogicPointsInCellCoordinates = MultidimensionalArray.Create(1, NbrInCell, SpatialDimension);
                    Grid.TransformGlobal2Local(GlobalPointsInCellCoordinates, LogicPointsInCellCoordinates, CellIndexGuess, 1, 0);
                    for (int f = 0; f < Fields.Count(); f++)
                    {
                        Fields[f].Evaluate(CellIndexGuess, 1, new NodeSet(((GridData)Grid).iGeomCells.GetRefElement(CellIndexGuess),
                            LogicPointsInCellCoordinates.ExtractSubArrayShallow(0, -1, -1), false), ResultInCell.ExtractSubArrayShallow(f, -1, -1));
                    }

                    indexInCell = 0;
                    for (int index = 0; index < NbrofPoints; index++)
                    {
                        if (InCellPointMask[index])
                        {
                            for (int f = 0; f < Fields.Count(); f++)
                            {
                                if (double.IsNaN(ResultInCell[f, 0, indexInCell]) || double.IsInfinity(ResultInCell[f, 0, indexInCell]))
                                    throw new Exception("NaN or Infinity not allowed");
                                Result[f, 0, index] = ResultInCell[f, 0, indexInCell];
                            }
                            indexInCell++;
                        }
                    }
                }
            }

            return Result.ExtractSubArrayShallow(-1,0,-1).CloneAs();
        }

        /// <summary>
        /// Computes the normal at a given point in the LevelSet Field
        /// </summary>
        /// <returns>The normal.</returns>
        /// <param name="Point">Point.</param>
        /// <param name="CellIndex">Cell index.</param>
        /// <param name="CellPosSure">If set to <c>true</c> cell position sure.</param>
        protected MultidimensionalArray EvaluateNormal(MultidimensionalArray Point, int CellIndex, bool CellPosSure)
        {
            MultidimensionalArray Buffer;
            MultidimensionalArray Normal = MultidimensionalArray.Create(1, SpatialDimension);
            for (int d = 0; d < SpatialDimension; d++)
            {
                Buffer = Evaluate(new DGField[] { InterfaceLevelSetGradient[d] }, Point, CellIndex, CellPosSure);
                Normal[0, d] = Buffer[0, 0];
            }
            // Calculation of eulerian gradient norm
            double Gradient_Norm = 0;
            for (int d = 0; d < SpatialDimension; d++)
            {
                Gradient_Norm += Normal[0, d] * Normal[0, d];
            }
            Gradient_Norm = Math.Sqrt(Gradient_Norm);
            for (int d = 0; d < SpatialDimension; d++)
            {
                Normal[0, d] /= Gradient_Norm;
            }
            return Normal;
        }

        /// <summary>
        /// Advection of all points with simple explicit euler
        /// </summary>
        /// <param name="timestep">Timestep.</param>
        /// <param name="PointMask">Point mask.</param>
        /// <param name="CellPosSure">If set to <c>true</c> cell position sure.</param>
        /// <param name="substeps">Substeps.</param>
        protected void Advect(double timestep, out List<CorrectionMaskType> PointMask, ref bool CellPosSure, int substeps = 1)
        {
            PointMask = new List<CorrectionMaskType>();

            DGField[] VelocityField = new DGField[2 * SpatialDimension];
            int index = 0;
            for (int dim = 0; dim < SpatialDimension; dim++, index += 2)
            {                
                VelocityField[index] = Velocity_Next[dim];
                VelocityField[index + 1] = Velocity_Current[dim];
            }

            MultidimensionalArray V;

            foreach (PointsofCellType CellObj in Points)
            {
                CellObj.MaxAdvectionStepLength = double.MinValue;
                int TotalNbrOfPoints = 0;
                Dictionary<int, int> NbrofPointsperType = new Dictionary<int, int>();
                foreach(KeyValuePair<int,List<PointType>> DictEntry in CellObj.ItemList)
                {
                    TotalNbrOfPoints += DictEntry.Value.Count;
                    NbrofPointsperType.Add(DictEntry.Key, DictEntry.Value.Count);
                }
                MultidimensionalArray PointsCoordinates = MultidimensionalArray.Create(TotalNbrOfPoints, SpatialDimension);
                int indx = 0;
                foreach (KeyValuePair<int, List<PointType>> DictEntry in CellObj.ItemList)
                {
                    for (int k=0;k<NbrofPointsperType[DictEntry.Key];k++)
                    {
                        for (int dim = 0; dim < SpatialDimension; dim++)
                        {
                            PointsCoordinates[indx, dim] = DictEntry.Value[k].Coordinates[0, dim];
                        }
                        indx++;
                    }
                }

                bool _CellPosSure = CellPosSure;
                _CellPosSure = false;
                for (int n = 0; n < substeps; n++, _CellPosSure = false)
                {
                    V = MultiEvaluatePointsInCell(VelocityField, PointsCoordinates, CellObj.CellIndex, _CellPosSure);

                    indx = 0;
                    foreach (KeyValuePair<int, List<PointType>> DictEntry in CellObj.ItemList)
                    {
                        foreach (PointType SinglePoint in DictEntry.Value)
                        {
                            if (SinglePoint.Active)
                            {
                                if (n == 0)
                                    SinglePoint.Coordinates_prev = SinglePoint.Coordinates.CloneAs();

                                // Loop over all substeps

                                for (int dim = 0; dim < SpatialDimension; dim++)
                                {
                                    // Numerical Timestep...weighted averaging over both velocitys and simple explicit Euler
                                    if (n == substeps - 1)
                                    {
                                        SinglePoint.Coordinates[0, dim] = PointsCoordinates[indx, dim] + (0.5 / (substeps * substeps)) * timestep * ((2 * (double)n + 1) * V[2 * dim, indx]
                                            + (2 * (double)(substeps - n) - 1) * V[2 * dim + 1, indx]);
                                    }
                                    else
                                    {
                                        PointsCoordinates[indx, dim] += (0.5 / (substeps * substeps)) * timestep * ((2 * (double)n + 1) * V[2 * dim, indx]
                                            + (2 * (double)(substeps - n) - 1) * V[2 * dim + 1, indx]);
                                    }
                                    CellPosSure = false;
                                }
                                if (n == substeps - 1)
                                {
                                    double AdvectionStepLength = 0;
                                    for (int dim = 0; dim < SpatialDimension; dim++)
                                        AdvectionStepLength += (SinglePoint.Coordinates[0, dim] - SinglePoint.Coordinates_prev[0, dim]).Pow2();
                                    SinglePoint.AdvectionStepLength = Math.Sqrt(AdvectionStepLength);
                                    if (CellObj.MaxAdvectionStepLength < AdvectionStepLength)
                                        CellObj.MaxAdvectionStepLength = AdvectionStepLength;
                                }
                            }
                            indx++;
                        }
                    }
                }
                PointMask.Add(NewFullCorrectionMask(CellObj.CellIndex));
            }
        }

        /// <summary>
        /// Checks the cell location of all Points and moves them to another cell if necessary
        /// </summary>
        /// <param name="CellPointItemMask">Particle mask.</param>
        /// <param name="CellPosSure">If set to <c>true</c> cell position sure.</param>
        protected void CorrectCellLocationOfPoints(List<CorrectionMaskType> CellPointItemMask, ref bool CellPosSure)
        {
            MPICollectiveWatchDog.Watch();

            //MPI start
            int PointListIndex;
            foreach(CorrectionMaskType SingleMask in CellPointItemMask)
            {
                PointListIndex = CelltoArrayindex[SingleMask.CellIndex];
                foreach(KeyValuePair<int,int> DictEntry in SingleMask.MovedPoints)
                {
                    if (CorrectCellLocation(PointListIndex,DictEntry.Key,DictEntry.Value)) break;
                }
            }
            CellPosSure = true;
        }

        /// <summary>
        /// Corrects the cell location of all the last "NbrofPointstoCorrect" points in a cell
        /// </summary>
        /// <returns><c>true</c>, if cell is empty after the function <c>false</c> otherwise.</returns>
        /// <param name="PointListIndex">Point list index.</param>
        /// <param name="signKey">Sign key.</param>
        /// <param name="NbrofPointstoCorrect">Nbrof particlesto correct.</param>
        private bool CorrectCellLocation(int PointListIndex,int signKey,int NbrofPointstoCorrect)
        {
            int listStart, listEnd;
            listStart = Points[PointListIndex].ItemList[signKey].Count - NbrofPointstoCorrect;
            listEnd = Points[PointListIndex].ItemList[signKey].Count;

            for (int singlePointIndex = listStart; singlePointIndex < listEnd; singlePointIndex++)
            {
                if (!((GridData)(Grid)).Cells.IsInCell(Points[PointListIndex].ItemList[signKey][singlePointIndex].Coordinates.ExtractSubArrayShallow(0, -1).To1DArray(),
                    Points[PointListIndex].CellIndex))
                {
                    int newCellIndex = FindCellIndexofPoint(Points[PointListIndex].ItemList[signKey][singlePointIndex].Coordinates, Points[PointListIndex].CellIndex, MaxCellIndexSearchDepth);
                    if (newCellIndex != -1)
                    {
                        if (MoveParticleMemoryToOtherCell(PointListIndex, ref singlePointIndex, ref listEnd, newCellIndex, signKey))
                            return true;
                    }
                    else
                    {
                        Points[PointListIndex].ItemList[signKey].RemoveAt(singlePointIndex--);
                        listEnd--;

                        bool EmptyCell = false;
                        foreach (KeyValuePair<int, List<PointType>> DictEntry in Points[PointListIndex].ItemList)
                            EmptyCell = DictEntry.Value.IsNullOrEmpty() ? true : false;

                        if (EmptyCell)
                        {
                            List<int> Keys = CelltoArrayindex.Keys.ToList();
                            for (int k = 0; k < Keys.Count; k++)
                            {
                                if (CelltoArrayindex[Keys[k]] > PointListIndex)
                                {
                                    CelltoArrayindex[Keys[k]] = CelltoArrayindex[Keys[k]] - 1;
                                }
                            }
                            CelltoArrayindex = new Dictionary<int, int>(CelltoArrayindex);
                            CelltoArrayindex.Remove(Points[PointListIndex].CellIndex);
                            Points.RemoveAt(PointListIndex--);
                            return true;
                        }
                    }
                }
            }
            return false;
        }

        /// <summary>
        /// Moves the point object from one cell object to another
        /// </summary>
        /// <returns><c>true</c> if cell object is empty afterwards and <c>false</c> otherwise.</returns>
        /// <param name="PointListIndex">Point list index.</param>
        /// <param name="OriginPositionofParticle">Origin positionof particle.</param>
        /// <param name="listEnd">List end.</param>
        /// <param name="TargetCellIndex">Target cell index.</param>
        /// <param name="signKey">Sign key.</param>
        private bool MoveParticleMemoryToOtherCell(int PointListIndex, ref int OriginPositionofParticle, ref int listEnd, int TargetCellIndex, int signKey)
        {
            if (!CelltoArrayindex.TryGetValue(TargetCellIndex, out int TargetCellListIndex))
            {
                Points.Add(NewPointsofCellObject(TargetCellIndex, Grid, MinimalDistanceSearch));
                CelltoArrayindex[TargetCellIndex] = Points.Count - 1;
                TargetCellListIndex = CelltoArrayindex[TargetCellIndex];
            }
            Points[TargetCellListIndex].ItemList[signKey].Insert(0, Points[PointListIndex].ItemList[signKey][OriginPositionofParticle]);
            Points[PointListIndex].ItemList[signKey].RemoveAt(OriginPositionofParticle--);
            listEnd--;

            bool EmptyCell = true;
            foreach (KeyValuePair<int, List<PointType>> DictEntry in Points[PointListIndex].ItemList)
                EmptyCell = DictEntry.Value.IsNullOrEmpty() ? EmptyCell : false;

            if (EmptyCell)
            {
                List<int> Keys = CelltoArrayindex.Keys.ToList();
                for (int k = 0; k < Keys.Count; k++)
                {
                    if (CelltoArrayindex[Keys[k]] > PointListIndex)
                    {
                        CelltoArrayindex[Keys[k]] = CelltoArrayindex[Keys[k]] - 1;
                    }
                }
                CelltoArrayindex = new Dictionary<int, int>(CelltoArrayindex);

                CelltoArrayindex.Remove(Points[PointListIndex].CellIndex);
                Points.RemoveAt(PointListIndex--);
                return true;
            }
            else return false;
        }

        /// <summary>
        /// Function to remove all Points that are inactive
        /// </summary>
        /// <param name="RemovedPointsCell">Removed points cell.</param>
        protected void RemoveNonActivePoints(out List<int> RemovedPointsCell)
        {
            RemovedPointsCell = new List<int>();
            for (int index = 0; index < Points.Count; index++)
            {
                bool PointsRemoved = false;
                bool CellIsEmpty = true;
                foreach (KeyValuePair<int, List<PointType>> DictEntry in Points[index].ItemList)
                {
                    for (int k = 0; k < DictEntry.Value.Count; k++)
                    {
                        if(DictEntry.Value[k].Active == false)
                        {
                            PointsRemoved = true;
                            DictEntry.Value.RemoveAt(k--);
                        }
                    }
                    if (!DictEntry.Value.IsNullOrEmpty())
                        CellIsEmpty = false;
                }
                if (CellIsEmpty)
                {
                    List<int> Keys = CelltoArrayindex.Keys.ToList();
                    for (int k = 0; k < Keys.Count; k++)
                    {
                        if (CelltoArrayindex[Keys[k]] > index)
                        {
                            CelltoArrayindex[Keys[k]] = CelltoArrayindex[Keys[k]] - 1;
                        }
                    }
                    CelltoArrayindex = new Dictionary<int, int>(CelltoArrayindex);
                    CelltoArrayindex.Remove(Points[index].CellIndex);
                    Points.RemoveAt(index--);
                }
                if (!CellIsEmpty && PointsRemoved)
                    RemovedPointsCell.Add(Points[index].CellIndex);
            }
        }

        /// <summary>
        /// Correction of the LevelSet
        /// </summary>
        /// <param name="timestepNo">Timestep no.</param>
        protected abstract void CorrectLevelSet(int timestepNo,CellMask CellsToCorrect = null,bool ReIteration = false);

        /// <summary>
        /// Projects a Point an the given LevelSet Isocontour in the normal direction of the LevelSet Field
        /// </summary>
        /// <returns><c>true</c>, if on level set value was projected, <c>false</c> otherwise.</returns>
        /// <param name="LevelSetValue">Level set value.</param>
        /// <param name="Coordinates">Coordinates.</param>
        /// <param name="CellIndex">Cell index.</param>
        /// <param name="CellPosSure">If set to <c>true</c> cell position sure.</param>
        /// <param name="LevelSet">Level set.</param>
        protected bool ProjectOnLevelSetValue(double LevelSetValue, MultidimensionalArray Coordinates, int CellIndex, ref bool CellPosSure, SinglePhaseField LevelSet = null)
        {
            SinglePhaseField Interface = this.Interface;
            SinglePhaseField[] InterfaceGrad = this.InterfaceLevelSetGradient.ToArray();
            if (LevelSet != null)
                Interface = LevelSet;

            BoundingBox CellBox = new BoundingBox(SpatialDimension);
            Grid.iGeomCells.GetCellBoundingBox(CellIndex, CellBox);
            double ERROR = CellBox.Diameter * 1e-3;

            MultidimensionalArray Coordinates_cpy = Coordinates;

            MultidimensionalArray Phi = Evaluate(new DGField[] { Interface }, Coordinates, CellIndex, CellPosSure);
            MultidimensionalArray PhiGrad = Evaluate(InterfaceGrad, Coordinates, CellIndex, CellPosSure);
            double stepWidth = Math.Max(PhiGrad.L2Norm(), 1.0);
            int signPhi = Math.Sign(Phi[0, 0]);
            if (Math.Abs(Phi[0, 0] - LevelSetValue) < 0) return true;
            MultidimensionalArray VectorToLevelSet = CalculateNormaltoLevelSetValue(LevelSetValue, Coordinates, CellIndex, CellPosSure, Interface);
            for (int dim = 0; dim < SpatialDimension; dim++) Coordinates[0, dim] = Coordinates[0, dim] + 1 / stepWidth * VectorToLevelSet[0, dim];
            CellIndex = FindCellIndexofPoint(Coordinates, CellIndex);
            if (CellIndex == -1) return false;
            CellPosSure = true;

            Grid.iGeomCells.GetCellBoundingBox(CellIndex, CellBox);
            ERROR = CellBox.Diameter * 1e-3;

            Phi = Evaluate(new DGField[] { Interface }, Coordinates, CellIndex, CellPosSure);

            int count = 0;
            MultidimensionalArray Coordinates_New = Coordinates.CloneAs();
            while (Math.Abs(Phi[0, 0] - LevelSetValue) > ERROR && count++ < 100)
            {
                VectorToLevelSet = CalculateNormaltoLevelSetValue(LevelSetValue, Coordinates, CellIndex, CellPosSure);
                PhiGrad = Evaluate(InterfaceGrad, Coordinates, CellIndex, CellPosSure);
                stepWidth = Math.Max(PhiGrad.L2Norm(), 1.0);
                for (int dim = 0; dim < SpatialDimension; dim++)
                    Coordinates_New[0, dim] = Coordinates[0, dim] + 1 / stepWidth * VectorToLevelSet[0, dim];
                CellPosSure = false;
                CellIndex = FindCellIndexofPoint(Coordinates_New, CellIndex);
                if (CellIndex == -1) return false;
                CellPosSure = true;
                MultidimensionalArray Phi_New;
                Phi_New = Evaluate(new DGField[] { Interface }, Coordinates_New, CellIndex, CellPosSure);

                if (Math.Abs(Phi_New[0, 0] - LevelSetValue) > Math.Abs(Phi[0, 0] - LevelSetValue) * 0.5 && Math.Abs(Phi_New[0, 0] - LevelSetValue) < Math.Abs(Phi[0, 0] - LevelSetValue))
                {
                    //double iterationStepSize = Math.Abs(Phi_New[0, 0] - LevelSetValue) / Math.Abs(Phi[0, 0] - LevelSetValue);
                    //for (int dim = 0; dim < SpatialDimension; dim++)
                    //    Coordinates_New[0, dim] = Coordinates[0, dim] + (1 / iterationStepSize) * VectorToLevelSet[0, dim];
                    if (Math.Sign(Phi_New[0, 0] - LevelSetValue) != Math.Sign(Phi[0, 0] - LevelSetValue))
                    {
                        double PhiDist = Math.Abs(Phi[0, 0] - LevelSetValue) + Math.Abs(Phi_New[0, 0] - LevelSetValue);
                        for (int dim = 0; dim < SpatialDimension; dim++)
                            Coordinates_New[0, dim] = Coordinates[0, dim] + (Coordinates_New[0, dim] - Coordinates[0, dim]) * (Math.Abs(Phi[0, 0] - LevelSetValue)/ PhiDist);
                    }
                    else
                    {
                        double PhiDist = Math.Abs(Phi[0, 0] - LevelSetValue) - Math.Abs(Phi_New[0, 0] - LevelSetValue);
                        for (int dim = 0; dim < SpatialDimension; dim++)
                            Coordinates_New[0, dim] = Coordinates[0, dim] + (Coordinates_New[0, dim] - Coordinates[0, dim]) * (Math.Abs(Phi[0, 0] - LevelSetValue) / PhiDist);
                    }
                    CellIndex = FindCellIndexofPoint(Coordinates_New, CellIndex);
                    if (CellIndex == -1) return false;
                    CellPosSure = true;
                    Phi_New = Evaluate(new DGField[] { Interface }, Coordinates_New, CellIndex, CellPosSure);
                }
                else if (Math.Abs(Phi_New[0, 0] - LevelSetValue) > Math.Abs(Phi[0, 0] - LevelSetValue))
                {
                    if (Math.Sign(Phi_New[0, 0] - LevelSetValue) != Math.Sign(Phi[0, 0] - LevelSetValue))
                    {
                        double PhiDist = Math.Abs(Phi[0, 0] - LevelSetValue) + Math.Abs(Phi_New[0, 0] - LevelSetValue);
                        //for (int dim = 0; dim < SpatialDimension; dim++)
                        //    Coordinates_New[0, dim] = Coordinates[0, dim] + (Math.Abs((Phi[0, 0] - LevelSetValue) / PhiDist)) * VectorToLevelSet[0, dim];
                        for (int dim = 0; dim < SpatialDimension; dim++)
                            Coordinates_New[0, dim] = Coordinates[0, dim] + (Coordinates_New[0, dim] - Coordinates[0, dim]) * (Math.Abs((Phi[0, 0] - LevelSetValue) / PhiDist));
                        CellIndex = FindCellIndexofPoint(Coordinates_New, CellIndex);
                        if (CellIndex == -1) return false;
                        CellPosSure = true;
                        Phi_New = Evaluate(new DGField[] { Interface }, Coordinates_New, CellIndex, CellPosSure);
                    }

                }

                Phi = Phi_New;
                for (int dim = 0; dim < SpatialDimension; dim++)
                    Coordinates[0, dim] = Coordinates_New[0, dim];
                CellPosSure = false;
                CellIndex = FindCellIndexofPoint(Coordinates_New, CellIndex);
                if (CellIndex == -1) return false;
                CellPosSure = true;
            }

            if (count == 101)
            {
                if (Math.Abs(Phi[0, 0] - LevelSetValue) < ERROR * 1e+2)
                {
                    return true;
                }
                else
                {
                    for (int dim = 0; dim < SpatialDimension; dim++)
                        Coordinates[0, dim] = Coordinates_cpy[0, dim];
                    return false;
                }
            }
            return true;
        }

        /// <summary>
        /// Calculates a vector from a point to the given LevelSetValue in the LevelSet Field.
        /// The vector direction is the direction of the LevelSet normal at the point.
        /// If the sign distance field is not completly correct the vector will not be not be completly correct too.
        /// </summary>
        /// <returns>The normalto level set value.</returns>
        /// <param name="LevelSetValue">Level set value.</param>
        /// <param name="Point">Coordinates.</param>
        /// <param name="CellIndex">Cell index.</param>
        /// <param name="CellPosSure">If set to <c>true</c> cell position sure.</param>
        /// <param name="LevelSet">Level set.</param>
        private MultidimensionalArray CalculateNormaltoLevelSetValue(double LevelSetValue, MultidimensionalArray Point, int CellIndex, bool CellPosSure, SinglePhaseField LevelSet = null)
        {
            SinglePhaseField Interface = this.Interface;
            if (LevelSet != null)
                Interface = LevelSet;

            //MultidimensionalArray Phi = Evaluate(new DGField[] { Interface }, Point, CellIndex, CellPosSure);
            //MultidimensionalArray Gradient = MultidimensionalArray.Create(1, 3);
            //MultidimensionalArray Buffer;
            //double DirectionLen = LevelSetValue - Phi[0, 0];
            //for (int d = 0; d < SpatialDimension; d++)
            //{
            //    Buffer = Evaluate(new DGField[] { InterfaceLevelSetGradient[d] }, Point, CellIndex, CellPosSure);
            //    Gradient[0, d] = Buffer[0, 0];
            //    //Console.WriteLine("Gradient " + d + " " + Buffer[0, 0]);
            //}
            //// Calculation of eulerian gradient norm
            //double Gradient_Norm = 0;
            //for (int d = 0; d < SpatialDimension; d++)
            //{
            //    Gradient_Norm += Gradient[0, d].Pow2();
            //}
            //Gradient_Norm = Math.Sqrt(Gradient_Norm);
            //MultidimensionalArray Normal = MultidimensionalArray.Create(1, SpatialDimension);
            //for (int d = 0; d < SpatialDimension; d++)
            //{
            //    Normal[0, d] = (Gradient[0, d] * DirectionLen) / Gradient_Norm;
            //}
            //return Normal;
            MultidimensionalArray Phi = Evaluate(new DGField[] { Interface }, Point, CellIndex, CellPosSure);
            MultidimensionalArray Gradient = Evaluate(InterfaceLevelSetGradient.ToArray(), Point, CellIndex, CellPosSure);
            double DirectionLen = LevelSetValue - Phi[0, 0];
            MultidimensionalArray Normal = Gradient.CloneAs();
            Normal.Scale(DirectionLen / Gradient.L2Norm());
            return Normal;
        }

        /// <summary>
        /// Searches for the cell of a point and returns its index starting by a index guess
        /// If the cell is not found -1 is returned.
        /// If no initial guess is given the whole grid is searched.
        /// The search range can be limited by setting maxSearchRange
        /// </summary>
        /// <returns>The cell index of the point.</returns>
        /// <param name="Point">Coordinate.</param>
        /// <param name="CellIndexGuess">Cell index guess.</param>
        /// <param name="maxSearchRange">Max search range.</param>
        protected int FindCellIndexofPoint (MultidimensionalArray Point, int CellIndexGuess = 0, int maxSearchRange = -1)
        {
            if (!((GridData)(Grid)).Cells.IsInCell(Point.ExtractSubArrayShallow(0, -1).To1DArray(), CellIndexGuess))
            {
                List<int> CellIndexFrontList = new List<int> { CellIndexGuess };
                List<int> CellIndexDoneList = new List<int> ();
                int SearchDepth = 0;
                while (SearchDepth < maxSearchRange || maxSearchRange == -1)
                {
                    SearchDepth++;
                    foreach(int j in CellIndexFrontList.ToList())
                    {
                        Grid.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaVertices, out int[] ResOut, out int[] ConnectingEntities);
                        CellIndexFrontList.Remove(j);
                        CellIndexDoneList.Add(j);
                        foreach(int i in ResOut)
                        {
                            if (!CellIndexDoneList.Contains(i) && !CellIndexFrontList.Contains(i))
                                CellIndexFrontList.Add(i);
                        }
                    }
                    if (CellIndexFrontList.IsNullOrEmpty())
                        break;
                    foreach(int cell in CellIndexFrontList)
                    {
                        if (((GridData)(Grid)).Cells.IsInCell(Point.ExtractSubArrayShallow(0, -1).To1DArray(), cell))
                        {
                            return cell;
                        }
                    }
                    if (SearchDepth > 100) Console.WriteLine("CellIndex of Point Search Algorithm is looping to often! Attention may loop until infinity if Point is not inside Mesh Boundaries");
                }
                return -1;
            }
            else
                return CellIndexGuess;
        }

        /// <summary>
        /// Returns all cells in a Radius around the BasisCell.
        /// </summary>
        /// <returns>The in radius around cell.</returns>
        /// <param name="BasicCellIndex">Cell index 0.</param>
        /// <param name="Radius">Radius.</param>
        protected List<int> CellsInRadiusAroundCell(int BasicCellIndex,double Radius)
        {
            BoundingBox Box = new BoundingBox(SpatialDimension);
            BoundingBox Box_Cell = new BoundingBox(SpatialDimension);
            Grid.iGeomCells.GetCellBoundingBox (BasicCellIndex, Box);
            double[] Max_0 = Box.Max;
            double[] Min_0 = Box.Min;
            double [] Max_Cell, Min_Cell;
            List<int> CellIndexDoneList = new List<int> ();
            List<int> CellIndexFrontList = new List<int> { BasicCellIndex };

            bool RadiusReached = false;
            while(!RadiusReached)
            {
                bool AllOutsideRadius = true;
                foreach(int FrontCell in CellIndexFrontList)
                {
                    Grid.iGeomCells.GetCellBoundingBox (FrontCell, Box_Cell);
                    Max_Cell = Box_Cell.Max;
                    Min_Cell = Box_Cell.Min;
                    bool OneDimOutsideRadius = false;
                    for(int dim=0;dim<SpatialDimension;dim++)
                    {
                        if (Math.Abs(Min_Cell[dim] - Max_0[dim]) > Radius && Math.Abs(Min_0[dim] - Max_Cell[dim]) > Radius)
                        {
                            OneDimOutsideRadius = true;
                        }
                    }
                    if (FrontCell == BasicCellIndex)
                        OneDimOutsideRadius = false;
                    if (!OneDimOutsideRadius)
                        //Console.WriteLine(FrontCell + " is inside Radius");
                    AllOutsideRadius = (OneDimOutsideRadius) ? AllOutsideRadius : false;
                }

                RadiusReached = (AllOutsideRadius) ? true : false;

                if (!RadiusReached)
                {
                    foreach (int FrontCell in CellIndexFrontList.ToList())
                    {
                        Grid.GetCellNeighbours(FrontCell, GetCellNeighbours_Mode.ViaVertices, out int[] Neighbours, out int[] ConnectingEntities);
                        CellIndexDoneList.Add(FrontCell);
                        CellIndexFrontList.Remove(FrontCell);
                        foreach (int neighbourcell in Neighbours)
                        {
                            if (!CellIndexDoneList.Contains(neighbourcell) && !CellIndexFrontList.Contains(neighbourcell))
                                CellIndexFrontList.Add(neighbourcell);
                        }
                    }
                }
            }
            return CellIndexDoneList;
        }

        /// <summary>
        /// Evaluates the LevelSet at a Point
        /// </summary>
        /// <returns>The phi.</returns>
        /// <param name="Point">Particle.</param>
        /// <param name="CellIndex">Cell index.</param>
        private MultidimensionalArray AdjustPhi(MultidimensionalArray Point, int CellIndex)
        {
            MultidimensionalArray Phi = Evaluate(new DGField[] { Interface }, Point, CellIndex, false);
            return Phi;
        }

        /// <summary>
        /// Evaluates the Normal at a Point in the level set field.
        /// </summary>
        /// <returns>The normal.</returns>
        /// <param name="Point">Particle.</param>
        /// <param name="CellIndex">Cell index.</param>
        /// <param name="CellPosSure">If set to <c>true</c> cell position sure.</param>
        public MultidimensionalArray AdjustNormal(MultidimensionalArray Point, int CellIndex, bool CellPosSure)
        {
            MultidimensionalArray Buffer;
            MultidimensionalArray Normal = MultidimensionalArray.Create(1, SpatialDimension);
            for (int d = 0; d < SpatialDimension; d++)
            {
                Buffer = Evaluate(new DGField[] { InterfaceLevelSetGradient[d] }, Point, CellIndex, CellPosSure);
                if(double.IsNaN(Buffer[0,0]) || double.IsInfinity(Buffer[0,0]))
                    throw new Exception("NaN or Infinity value not allowed");
                Normal[0, d] = Buffer[0, 0];
            }
            // Calculation of eulerian gradient norm
            return NormalizeVector (Normal);
        }

        /// <summary>
        /// Calculates the radius of the level set zero level contour by comparing normal vectors
        /// Only used for FrontTracking Approach
        /// </summary>
        /// <returns>The radius.</returns>
        /// <param name="Coordinate">Coordinate.</param>
        /// <param name="CellIndexGuess">Cell index guess.</param>
        protected double CalculateRadius (MultidimensionalArray Coordinate, int CellIndexGuess)
        {
            BoundingBox Box = new BoundingBox(SpatialDimension);
            Grid.iGeomCells.GetCellBoundingBox (CellIndexGuess, Box);
            double delta_dist = 0.01*Box.Diameter;
            MultidimensionalArray Nz = AdjustNormal (Coordinate, CellIndexGuess, false);
            MultidimensionalArray N_p;
            double theta;
            double R = Double.MaxValue;

            if (SpatialDimension == 1) throw new Exception ("Interface Curvature does not exist in 0-dimensional Surfaces");
            else if (SpatialDimension == 2) {
                MultidimensionalArray Coord_p = MultidimensionalArray.Create (1, SpatialDimension);
                int sign = 1;
                for (int i = 0; i < 2; i++, sign -= 2) {
                    Coord_p [0, 0] = Coordinate [0, 0] + sign * Nz [0, 1]*delta_dist;
                    Coord_p [0, 1] = Coordinate [0, 1] + sign * -Nz [0, 0]*delta_dist;
                    bool CellPosSure = false;
                    if (ProjectOnLevelSetValue (0.0,Coord_p, CellIndexGuess, ref CellPosSure)) {
                        N_p = AdjustNormal (Coord_p, CellIndexGuess, false);
                        theta = VectorAngle (Nz, N_p);
                        double distance = VectorNorm(ComputeVectorbetweenPoints(Coordinate, Coord_p));
                        R = Math.Min ((delta_dist / (2 * Math.Sin (theta / 2))), R);
                    }
                }
                return R;
            } else if (SpatialDimension == 3) {
                MultidimensionalArray Ny, Nx;
                MultidimensionalArray V1 = MultidimensionalArray.CreateWrapper (new double [3] { 0, 0, 1 }, new int [2] { 1, 3 });
                MultidimensionalArray V2 = MultidimensionalArray.CreateWrapper (new double [3] { 0, 1, 0 }, new int [2] { 1, 3 });
                if (VectorAngle (Nz, V1) > Math.PI / 6) {
                    Ny = NormalizeVector (CrossProduct (Nz, V1));
                } else if (VectorAngle (Nz, V2) > Math.PI / 6) {
                    Ny = NormalizeVector (CrossProduct (Nz, V2));
                } else throw new Exception (" Failure in Compuation of Coordinate System");
                Nx = NormalizeVector (CrossProduct (Nz, Ny));

                MultidimensionalArray Coord_p = MultidimensionalArray.Create (1, 3);
                double sigma = 0.0;
                for (int i = 0; i < 8; i++, sigma += (Math.PI / 4)) {
                    Coord_p [0, 0] = Coordinate [0, 0] + (Nx [0, 0] * Math.Cos (sigma) + Ny [0, 0] * Math.Sin (sigma)) * delta_dist;
                    Coord_p [0, 1] = Coordinate [0, 1] + (Nx [0, 1] * Math.Cos (sigma) + Ny [0, 1] * Math.Sin (sigma)) * delta_dist;
                    Coord_p [0, 2] = Coordinate [0, 2] + (Nx [0, 2] * Math.Cos (sigma) + Ny [0, 2] * Math.Sin (sigma)) * delta_dist;
                    bool CellPosSure = false;
                    if (ProjectOnLevelSetValue (0.0,Coord_p, CellIndexGuess, ref CellPosSure)) {
                        N_p = AdjustNormal (Coord_p, CellIndexGuess, false);
                        theta = VectorAngle (Nz, N_p);
                        R = Math.Min ((delta_dist / (2 * Math.Sin (theta / 2))), R);
                    }
                }
                return R;
            } else throw new Exception ("Curvature Calculation for this Dimension not implemented");
        }

        /// <summary>
        /// Returns the vector from V2 to V1.
        /// </summary>
        /// <returns>The vector between points.</returns>
        /// <param name="V1">V1.</param>
        /// <param name="V2">V2.</param>
        public static MultidimensionalArray ComputeVectorbetweenPoints (MultidimensionalArray V1, MultidimensionalArray V2)
        {
            if (V1.GetLength (0) != 1 || V2.GetLength (0) != 1) throw new Exception ("Only one Vector per Parameter allowed!");
            if (V1.GetLength (1) != V2.GetLength (1)) throw new Exception ("Inconsitent Dimension in Connection!");

            MultidimensionalArray Connection = MultidimensionalArray.Create (1, V1.GetLength (1));
            for (int dim = 0; dim < V1.GetLength (1); dim++)
            {
                if (double.IsNaN(V1[0, dim]) || double.IsInfinity(V1[0, dim])) throw new Exception("Double NaN or Infinity not allowed");
                if (double.IsNaN(V2[0, dim]) || double.IsInfinity(V2[0, dim])) throw new Exception("Double NaN or Infinity not allowed");
                Connection[0, dim] = V1 [0, dim] - V2 [0, dim];
            }
            return Connection;
        }

        /// <summary>
        /// Returns the dotproduct of two vectors.
        /// </summary>
        /// <returns>The product.</returns>
        /// <param name="V1">V1.</param>
        /// <param name="V2">V2.</param>
        public static double DotProduct (MultidimensionalArray V1, MultidimensionalArray V2)
        {
            if (V1.GetLength (0) != 1 || V2.GetLength (0) != 1) throw new Exception ("Only one Vector allowed!");
            if (V1.GetLength (1) != V2.GetLength (1)) throw new Exception ("Inconsitent Dimension in DotProduct!");
            double result = 0;
            for (int dim = 0; dim < V1.GetLength(1); dim++)
            {
                if (double.IsNaN(V1[0, dim]) || double.IsInfinity(V1[0, dim])) throw new Exception("Double NaN or Infinity not allowed");
                if (double.IsNaN(V2[0, dim]) || double.IsInfinity(V2[0, dim])) throw new Exception("Double NaN or Infinity not allowed");
                result += V1 [0, dim] * V2 [0, dim];
            }
            return result;
        }

        /// <summary>
        /// Multiplys two vectors with a scalar each and adds them
        /// </summary>
        /// <returns>The addition.</returns>
        /// <param name="alpha1">Alpha1.</param>
        /// <param name="V1">V1.</param>
        /// <param name="alpha2">Alpha2.</param>
        /// <param name="V2">V2.</param>
        public static MultidimensionalArray VectorAddition (double alpha1, MultidimensionalArray V1, double alpha2, MultidimensionalArray V2)
        {
            if (V1.GetLength (0) != 1) throw new Exception("First Vector Paramter has more than one Vector");
            if (V2.GetLength (0) != 1) throw new Exception ("Second Vector Paramter has more than one Vector");
            if (V1.GetLength (1) != V2.GetLength (1)) throw new Exception ("Inconsitent Dimension in Addition!");
            MultidimensionalArray result = MultidimensionalArray.Create (1, V1.GetLength (1));
            for (int dim = 0; dim < V1.GetLength (1); dim++)
            {
                if (double.IsNaN(V1[0, dim]) || double.IsInfinity(V1[0, dim])) throw new Exception("Double NaN or Infinity not allowed");
                if (double.IsNaN(V2[0, dim]) || double.IsInfinity(V2[0, dim])) throw new Exception("Double NaN or Infinity not allowed");
                result[0,dim] = alpha1 * V1 [0, dim] + alpha2 * V2 [0, dim];
                //Console.WriteLine(alpha1 + " * " + V1[0, dim] + " + " + alpha2 + " * " + V2[0, dim] + " = " + result[0, dim]);
            }
            return result;
        }

        /// <summary>
        /// Inverts the direction of a vector
        /// </summary>
        /// <returns>The vector.</returns>
        /// <param name="V1">V1.</param>
        public static MultidimensionalArray InvertVector (MultidimensionalArray V1)
        {
            for (int dim = 0; dim < V1.GetLength (1); dim++)
            {
                if (double.IsNaN(V1[0, dim]) || double.IsInfinity(V1[0, dim])) throw new Exception("Double NaN or Infinity not allowed");
                V1[0, dim] *= -1;
            }
            return V1;
        }

        /// <summary>
        /// Computes cross product of two vectors
        /// </summary>
        /// <returns>The product.</returns>
        /// <param name="V1">V1.</param>
        /// <param name="V2">V2.</param>
        public static MultidimensionalArray CrossProduct(MultidimensionalArray V1, MultidimensionalArray V2)
        {
            if (V1.GetLength (0) != 1 || V2.GetLength (0) != 1) throw new Exception ("Only one Vector allowed!");
            if (V1.GetLength (1) != V2.GetLength (1)) throw new Exception ("Inconsitent Dimension in DotProduct!");
            if (V1.GetLength (1) != 3 || V2.GetLength (1) != 3) throw new Exception ("Vector for Cross Product must be 3 dimensional");
            MultidimensionalArray result = MultidimensionalArray.Create(1, 3);
            result [0,0] = V1 [0,1] * V2 [0,2] - V1 [0,2] * V2 [0,1];
            result [0,1] = V1 [0,2] * V2 [0,0] - V2 [0,2] * V1 [0,0];
            result [0,2] = V1 [0,0] * V2 [0,1] - V1 [0,1] * V2 [0,0];
            return result;
        }

        /// <summary>
        /// Computes eulerian norm of vector.
        /// </summary>
        /// <returns>The norm.</returns>
        /// <param name="V1">V1.</param>
        public static double VectorNorm(MultidimensionalArray V1)
        {
            if (V1.Dimension == 1)
            {
                double result = 0;
                for (int dim = 0; dim < V1.GetLength(0); dim++)
                {
                    if (double.IsNaN(V1[dim]) || double.IsInfinity(V1[dim])) throw new Exception("Double NaN or Infinity not allowed");
                    result += V1[dim].Pow2();
                }
                return Math.Sqrt(result);
            }
            else if (V1.Dimension == 2)
            {
                double result = 0;
                for (int dim = 0; dim < V1.GetLength(1); dim++)
                {
                    if (double.IsNaN(V1[0, dim]) || double.IsInfinity(V1[0, dim])) throw new Exception("Double NaN or Infinity not allowed");
                    result += V1[0, dim].Pow2();
                }
                return Math.Sqrt(result);
            }
            return 0;
        }

        /// <summary>
        /// Normalizes vector to length one.
        /// </summary>
        /// <returns>The vector.</returns>
        /// <param name="V1">V1.</param>
        public static MultidimensionalArray NormalizeVector(MultidimensionalArray V1)
        {
            if (V1.Dimension == 1)
            {
                double norm = VectorNorm(V1);
                if (norm == 0.0)
                    throw new Exception("Vector has Length 0; Division with zero");
                for (int dim = 0; dim < V1.GetLength(0); dim++)
                {
                    if (double.IsNaN(V1[dim]) || double.IsInfinity(V1[dim])) throw new Exception("Double NaN or Infinity not allowed");
                    V1[dim] /= norm;
                }
                if (Math.Abs(VectorNorm(V1) - 1) > 1e-2) Console.WriteLine(VectorNorm(V1));
                return V1;
            }
            else if (V1.Dimension == 2)
            {
                double norm = VectorNorm(V1);
                if (norm == 0.0)
                    throw new Exception("Vector has Length 0; Division with zero");
                for (int dim = 0; dim < V1.GetLength(1); dim++)
                {
                    if (double.IsNaN(V1[0, dim]) || double.IsInfinity(V1[0, dim])) throw new Exception("Double NaN or Infinity not allowed");
                    V1[0, dim] /= norm;
                }
                if (Math.Abs(VectorNorm(V1) - 1) > 1e-2) Console.WriteLine(VectorNorm(V1));
                return V1;
            }
            return null;
        }

        /// <summary>
        /// Computes the angle between two vectors.
        /// </summary>
        /// <returns>the angle.</returns>
        /// <param name="V1">V1.</param>
        /// <param name="V2">V2.</param>
        public static double VectorAngle(MultidimensionalArray V1, MultidimensionalArray V2)
        {
            return Math.Acos (DotProduct (V1, V2) / (VectorNorm (V1) * VectorNorm (V2)));
        }

        /// <summary>
        /// Prints the Point Mask.
        /// </summary>
        /// <param name="PointMask">Point mask.</param>
        public void PrintPointMask(List<CorrectionMaskType> PointMask)
        {
            Console.WriteLine("------------Point Mask------------");
            foreach (CorrectionMaskType SingleMask in PointMask)
            {

                Console.Write("CellIndex:"+SingleMask.CellIndex + "|");
                foreach(KeyValuePair<int,int> DictEntry in SingleMask.MovedPoints)
                {
                    Console.Write(" (" + DictEntry.Key + ") " + DictEntry.Value);
                }
                Console.WriteLine();
            }
            Console.WriteLine("---------------------------------");
        }

        /// <summary>
        /// Prints the total number of Points.
        /// </summary>
        public void PointCount()
        {
            int add = 0;
            foreach (PointsofCellType CellObj in Points)
            {
                foreach(KeyValuePair<int,List<PointType>> DictEntry in CellObj.ItemList)
                {
                    add += DictEntry.Value.Count;
                }
            }
            Console.WriteLine("Number of Particles: " + add + " in " + CelltoArrayindex.Count() + " Cells");
        }

        /// <summary>
        /// Prints the Points Dictionary and List.
        /// </summary>
        public void PrintDictionaryandList()
        {
            Console.WriteLine("-----------------------------------------------------Point Dictionary and List------------------------------------------------------");
            foreach (KeyValuePair<int, int> entry in CelltoArrayindex)
                Console.WriteLine("CellIndex: " + entry.Key + " Listindex: " + entry.Value);

            foreach(PointsofCellType CellObj in Points)
            {
                Console.Write(CellObj.PrintToString());
            }
            Console.WriteLine("------------------------------------------------------------------------------------------------------------------------------------");
        }

        /// <summary>
        /// Prints the number of points in a cell for all cells with points inside.
        /// </summary>
        public void PrintPointinCell()
        {
            Console.WriteLine("-----------Point in Cell-----------");
            foreach (PointsofCellType CellObj in Points)
            {
                Console.Write("CellIndex:" + CellObj.CellIndex + "|");
                foreach (KeyValuePair<int,List<PointType>> DictEntry in CellObj.ItemList)
                {
                    Console.Write(" (" + DictEntry.Key + ") " + DictEntry.Value.Count());
                }
                Console.WriteLine();
            }
            Console.WriteLine("-----------------------------------");
        }

        /// <summary>
        /// Prints a .csv file of all Points with properties.
        /// </summary>
        /// <param name="Name">Name.</param>
        public void PrintCSV(string Name)
        {
            StringBuilder file = new StringBuilder();
            StringBuilder text = new StringBuilder();
            string header=null;
            foreach (PointsofCellType CellObj in Points)
            {
                text.Append(CellObj.PrintToCSV(out header));
            }
            file.Append(text);
            File.WriteAllText(Name + ".csv", file.ToString());
        }        
    }
}