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
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.Utils.Geom;
using BoSSS.Solution.Statistic;
using ilPSP;

namespace BoSSS.Solution.LevelSetTool.SemiLagrangianLevelSet
{
    /// <summary>
    /// Implementation of cell class for marker cells.
    /// Point Type Key is 1
    /// </summary>
    public class MarkersofCell : SemiLagrangianLevelSetMethod
        <MarkersofCell, SingleLvlSetMarker, MarkerCorrectionMask>.ItemsofCell
    {
        /// <summary>
        /// Initializes a new instance of the <see cref="T:SemiLagrangianLevelSet.MarkersofCell"/> class.
        /// </summary>
        /// <param name="CellIndex">Cell index.</param>
        /// <param name="Grid">Grid.</param>
        /// <param name="Correction">Correction.</param>
        public MarkersofCell(int CellIndex, IGridData Grid, MarkerLevelSet.MinimalDistanceSearchMode Correction)
            : base(CellIndex, new int[1] { 1 })
        {
            Box = new BoundingBox(Grid.SpatialDimension);
            Grid.iGeomCells.GetCellBoundingBox(CellIndex, Box);
        }

        /// <summary>
        /// Adjusting LevelSet value and LevelSet normal.
        /// </summary>
        /// <param name="ResultValues">Result values.</param>
        /// <param name="NbrofPointsPerType">Nbrof points per type.</param>
        public override void AdjustParametersofCellPoints(MultidimensionalArray ResultValues, Dictionary<int, int> NbrofPointsPerType)
        {
            int SpatialDimension = ResultValues.GetLength(0) - 1;
            MultidimensionalArray LevelSet = MultidimensionalArray.Create(1, 1);
            MultidimensionalArray LevelSetGradient = MultidimensionalArray.Create(1, SpatialDimension);
            int index = 0;

            foreach (KeyValuePair<int, int> DictEntry in NbrofPointsPerType)
            {
                if (NbrofPointsPerType[DictEntry.Key] != ItemList[DictEntry.Key].Count)
                    throw new Exception("Number of Points in Parameter adjustment does not fit");
                for (int k = 0; k < ItemList[DictEntry.Key].Count; k++)
                {
                    LevelSet[0, 0] = ResultValues[0, index];
                    for (int dim = 0; dim < SpatialDimension; dim++)
                    {
                        LevelSetGradient[0, dim] = ResultValues[dim + 1, index];
                    }
                    ItemList[DictEntry.Key][k].AdjustParameters(LevelSet, LevelSetGradient);
                    index++;
                }
            }
        }

        /// <summary>
        /// Printing marker properties to csv.
        /// </summary>
        /// <returns>The to csv.</returns>
        /// <param name="header">Header.</param>
        public override string PrintToCSV(out string header)
        {
            string lines = null;
            string header_part0 = "CellIndex,";
            string header_part1 = null;
            foreach (KeyValuePair<int, List<SingleLvlSetMarker>> DictEntry in ItemList)
            {
                foreach (SingleLvlSetMarker Point in DictEntry.Value)
                {
                    lines += CellIndex.ToString() + " , ";
                    lines += Point.PrintToCSV(out header_part1);
                    lines += "\n";
                }
            }
            header = header_part0 + header_part1 + "\n";

            return lines;
        }

        /// <summary>
        /// Printing marker properties to string.
        /// </summary>
        /// <returns>The to string.</returns>
        public override string PrintToString()
        {
            string lines = null;
            foreach (KeyValuePair<int, List<SingleLvlSetMarker>> DictEntry in ItemList)
            {
                foreach (SingleLvlSetMarker Point in DictEntry.Value)
                {
                    lines += "CellIndex:" + string.Format("{0,4}", CellIndex);
                    lines += Point.PrintToString();
                    lines += "\n";
                }
            }
            return lines;
        }
    }

    /// <summary>
    /// Implementation of Point class for markers.
    /// </summary>
    public class SingleLvlSetMarker : SemiLagrangianLevelSetMethod
        <MarkersofCell, SingleLvlSetMarker, MarkerCorrectionMask>.SinglePoint
    {

        /// <summary>
        /// Initializes a new instance of the <see cref="T:SemiLagrangianLevelSet.SingleLvlSetMarker"/> class.
        /// </summary>
        /// <param name="Coordinates">Coordinates.</param>
        /// <param name="AdvectionStepLengthINI">Advection step length ini.</param>
        public SingleLvlSetMarker(MultidimensionalArray Coordinates,double AdvectionStepLengthINI)
            :base(Coordinates,AdvectionStepLengthINI)
        { }

        /// <summary>
        /// Adjusts the parameters of the marker
        /// </summary>
        /// <param name="LevelSet">Level set.</param>
        /// <param name="LevelSetGradient">Level set gradient.</param>
        public new void AdjustParameters(MultidimensionalArray LevelSet,MultidimensionalArray LevelSetGradient)
        {
            base.AdjustParameters(LevelSet, LevelSetGradient);
            MultidimensionalArray Normal_old=null;
            if (Normal != null)
            {
                Normal_old = Normal.CloneAs();
                int sign = Math.Sign(SemiLagrangianLevelSetMethod<MarkersofCell, SingleLvlSetMarker, MarkerCorrectionMask>.
                    DotProduct(Normal, Normal_old));
                if (sign <= 0 && Active && Normal != null)
                {
                    Active = false;
                    Console.WriteLine("Normal Inversion");
                }
            }
        }

        /// <summary>
        /// Prints marker properties 
        /// </summary>
        /// <returns>String of marker properties</returns>
        public new string PrintToString()
        {
            string line = base.PrintToString();
            line += "|Active:" + Active + " ";
            return line;
        }

        /// <summary>
        /// Prints marker properties to csv.
        /// </summary>
        /// <returns>String of marker properties</returns>
        /// <param name="header">Header.</param>
        public new string PrintToCSV(out string header)
        {
            string line = base.PrintToCSV(out header);
            if (Active)
                line += 1;
            else
                line += 0;
            line += ",  " + DeactivationSign;
            header += ",active,DeactivationSign";
            return line;
        }
    }

    /// <summary>
    /// Marker correction mask implementation
    /// </summary>
    public class MarkerCorrectionMask : SemiLagrangianLevelSetMethod
        <MarkersofCell, SingleLvlSetMarker, MarkerCorrectionMask>.CorrectionMask
    {

        /// <summary>
        /// Initializes a new instance of the <see cref="T:SemiLagrangianLevelSet.MarkerCorrectionMask"/> class.
        /// Empty initialization.
        /// </summary>
        /// <param name="CellIndex">Cell index.</param>
        public MarkerCorrectionMask(int CellIndex)
            : base(CellIndex, new int[1] { 1 })
        { }

        /// <summary>
        /// Initializes a new instance of the <see cref="T:SemiLagrangianLevelSet.MarkerCorrectionMask"/> class.
        /// Initilization with values.
        /// </summary>
        /// <param name="CellIndex">Cell index.</param>
        /// <param name="NbrMovedPoints">Nbr moved points.</param>
        public MarkerCorrectionMask(int CellIndex, int[] NbrMovedPoints)
            : base(CellIndex, new int[1]{1})
        {
            if (NbrMovedPoints.Length != 1)
                throw new Exception("MarkerCorrectionMask constructor needs exactly one initial values");
            MovedPoints[1] = NbrMovedPoints[0];
        }
    }

    /// <summary>
    /// Marker k-d-tree filter.
    /// </summary>
    public class Marker_K_D_TreeFilter : SemiLagrangianLevelSetMethod
        <MarkersofCell, SingleLvlSetMarker, MarkerCorrectionMask>.K_D_TreeFilter
    {

        /// <summary>
        /// Initializes a new instance of the <see cref="T:SemiLagrangianLevelSet.Marker_K_D_TreeFilter"/> class.
        /// </summary>
        /// <param name="Points">Points.</param>
        /// <param name="NbrOfPointsinNode">Nbr of pointsin node.</param>
        public Marker_K_D_TreeFilter(List<SingleLvlSetMarker> Points, int NbrOfPointsinNode)
            : base(Points, NbrOfPointsinNode)
        { }

        /// <summary>
        /// Nearest point calculation for markers
        /// </summary>
        /// <param name="Points">Points.</param>
        /// <param name="GlobalNode">Global node.</param>
        /// <param name="NearestPoint">Nearest point.</param>
        /// <param name="MinRadius">Minimum radius.</param>
        protected override void NearestPointAndMinRadiusFromList(List<SingleLvlSetMarker> Points, MultidimensionalArray GlobalNode, ref SingleLvlSetMarker NearestPoint, ref double MinRadius)
        {
            foreach (SingleLvlSetMarker Point in Points)
            {
                double dist = 0;
                for (int dim = 0; dim < SpatialDimension; dim++)
                {
                    dist += (Point.Coordinates[0, dim] - GlobalNode[dim]).Pow2();
                }
                dist = Math.Sqrt(dist);
                if (MinRadius > dist)
                {
                    NearestPoint = Point;
                    MinRadius = dist;
                }
            }
        }
    }

    public class MarkerLevelSet : SemiLagrangianLevelSetMethod
        <MarkersofCell, SingleLvlSetMarker, MarkerCorrectionMask>
    {
        public enum MarkerCorrectionMode { Gaussian, MinDistanceWithNormal, MinDistanceWithSignUpdate };
        protected readonly MarkerCorrectionMode MarkerCorrection;

        private readonly double DELTA_PHI;
        private readonly double KERNEL_RADIUS;
        private OldNearestPointMemory PreviousNearestMarker;
        
        public MarkerLevelSet(VectorField<SinglePhaseField> Velocity_Next, VectorField<SinglePhaseField> Velocity_Current, SinglePhaseField Interface, LevelSetTracker InterfaceTrck,
             VectorField<SinglePhaseField> LevelSetGradient, int TargetNbrMarkersPerCell, int UpperLimitMarkersPerCell, int LowerLimitMarkersPerCell,
             int max_AddandAttract_Iteration, int NarrowBandWidth, IGridData Grid, int ReseedingInterval, MinimalDistanceSearchMode LevelSetCorrection,
             MarkerCorrectionMode MarkerCorrection, TopologyMergingMode TopologyMerging, NormalVectorDampingMode NormalVectorDamping,
             double DELTA_PHI=0, double KERNEL_RADIUS=0)

            : base(Velocity_Next, Velocity_Current, Interface, InterfaceTrck, LevelSetGradient, NarrowBandWidth, Grid,
                 ReseedingInterval, LevelSetCorrection, TopologyMerging, NormalVectorDamping, TargetNbrMarkersPerCell, UpperLimitMarkersPerCell, LowerLimitMarkersPerCell)
        {
            if (TargetNbrMarkersPerCell <= 0) throw new Exception("Target number of markers per cell must not be lower than 1.");
            if (TargetNbrMarkersPerCell > UpperLimitMarkersPerCell) throw new Exception("The upper limit of markers per cell must be equal or higher than the target number of particles per cell.");
            if (TargetNbrMarkersPerCell < LowerLimitMarkersPerCell) throw new Exception("The lower limit of markers per cell must be equal or lower than the target number of particles per cell.");
            this.MarkerCorrection = MarkerCorrection;
            this.DELTA_PHI = DELTA_PHI;
            this.KERNEL_RADIUS = KERNEL_RADIUS;
            ReseedingWdith = 0;
            if (NarrowBandWidth < 1)
                throw new Exception("Narrow Band Width must be greater than 0");
        }

        /// <summary>
        /// Initialize this instance.
        /// </summary>
        public override void Initialize()
        {
            bool CellPosSure = false;
            Console.WriteLine("Particles initialized");
            sw.Start();
            InitializePointstoLevelSet(out List<MarkerCorrectionMask> CorrectCellPosMask, ref CellPosSure);
            sw.Stop();
            Console.WriteLine("InitializeParticles and Reseed Elapsed:  {0}", sw.Elapsed);
            sw.Reset();
            //PrintCSV("Marker_Init");
            sw.Start();
            CorrectCellLocationOfPoints(CorrectCellPosMask, ref CellPosSure);
            sw.Stop();
            Console.WriteLine("Correct Cell Location Elapsed:           {0}", sw.Elapsed);
            //PrintCSV("Marker_Corrected");
            sw.Reset();
            //PrintPointinCell();
            sw.Start();
            AdjustPointParamters(0, CellPosSure);
            sw.Stop();
            Console.WriteLine("Adjust Parameters Elapsed:               {0}", sw.Elapsed);
            //PrintCSV("Marker_" + 0);
            //Console.WriteLine("Wrote Points to file:                    Marker_{0}.csv", 0);
            RemoveNonActivePoints(out _);
        }

        public new void PerformTimestep(double dt, int substeps, int Timestep)
        {
            bool Reseed = true;
            base.PerformTimestep(dt, substeps, Timestep);
            bool CellPosSure = true;
            sw.Start();
            Advect(dt, out List<MarkerCorrectionMask> CorrectCellPosMask, ref CellPosSure, substeps);
            sw.Stop();
            Console.WriteLine("Advected Markers:                        {0}", sw.Elapsed);
            sw.Reset();
            sw.Start();
            CorrectCellLocationOfPoints(CorrectCellPosMask, ref CellPosSure);
            sw.Stop();
            Console.WriteLine("Correct Cell Location Elapsed:           {0}", sw.Elapsed);
            sw.Reset();
            PrintCSV("Marker_" + Timestep);
            //Console.WriteLine("Wrote Points to file:                    Marker_{0}.csv", Timestep);
            sw.Start();
            CorrectLevelSet(Timestep);
            sw.Stop();
            Console.WriteLine("Correct Level Set Elapsed:               {0}", sw.Elapsed);
            sw.Reset();
            sw.Start();
            AdjustPointParamters(Timestep, CellPosSure);
            sw.Stop();
            Console.WriteLine("Adjust Parameters Elapsed:               {0}", sw.Elapsed);
            sw.Reset();
            if (Points.Count == 0)
                Reseed = false;

            if ((Timestep % ReseedingInterval) == 0 && ReseedingInterval != -1 && Reseed)
            {
                Console.WriteLine("--------------------Reseeding--------------------");
                sw.Start();
                ReseedingPoints(out List<MarkerCorrectionMask> ReseedingMarkerMask);
                sw.Stop();
                Console.WriteLine("Reseed Markers Elapsed:                  {0}", sw.Elapsed);
                sw.Reset();
                sw.Start();
                CorrectCellLocationOfPoints(ReseedingMarkerMask, ref CellPosSure);
                sw.Stop();
                Console.WriteLine("Correct Cell Location Elapsed:           {0}", sw.Elapsed);
                sw.Reset();
                sw.Start();
                AdjustPointParamters(Timestep, CellPosSure);
                sw.Stop();
                Console.WriteLine("Adjust Parameters Elapsed:               {0}", sw.Elapsed);
                sw.Reset();
            }
            
            if (ReseedingInterval != -1 && (((Timestep + ReseedingInterval / 2) % ReseedingInterval) == 0) && Reseed)
            {
                Console.WriteLine("--------------------Reattraction--------------------");
                sw.Start();
                ReattractMarkers(null);
                sw.Stop();
                Console.WriteLine("Reattract Markers Elapsed:               {0}", sw.Elapsed);
                sw.Reset();                
            }            
        }

        /// <summary>
        /// Removes the surplus points in a random pattern
        /// </summary>
        /// <param name="TreatedCell">Treated cell.</param>
        /// <param name="ParticleMask">Particle mask.</param>
        protected override void RemoveSurplusPoints(MarkersofCell TreatedCell, List<MarkerCorrectionMask> ParticleMask)
        {
            int NbrofSurplusMarkers = TreatedCell.ItemList[1].Count - TargetNbrofPointsPerCellPerDimension * SpatialDimension;
            int count = 0;
            for (int k = 0; k < NbrofSurplusMarkers; k++)
            {
                int indx = (int)(CoordGen.NextDouble() * TreatedCell.ItemList[1].Count());
                TreatedCell.ItemList[1].RemoveAt(indx);
                count++;
            }
            //Console.WriteLine("Removed " + count + " Points");
        }

        /// <summary>
        /// Projects a Point on the LevelSet zero Isocontour.
        /// </summary>
        /// <returns><c>true</c> was succesfully projected and <c>false</c> otherwise.</returns>
        /// <param name="TreatedCell">Treated cell.</param>
        /// <param name="signKey">Sign key.</param>
        /// <param name="Coordinate">Coordinate.</param>
        /// <param name="ParticleMask">Particle mask.</param>
        protected override bool MovePointtoTarget(MarkersofCell TreatedCell, int signKey, MultidimensionalArray Coordinate, List<MarkerCorrectionMask> ParticleMask)
        {
            BoundingBox Box = new BoundingBox(SpatialDimension);
            Grid.iGeomCells.GetCellBoundingBox(TreatedCell.CellIndex, Box);
            double AdvectionStepLengthINI = Box.Diameter * 1e-2;

            bool CellPosSure = true;
            if (ProjectOnLevelSetValue(0.0, Coordinate, TreatedCell.CellIndex, ref CellPosSure))
            {
                TreatedCell.ItemList[1].Add(new SingleLvlSetMarker(Coordinate, AdvectionStepLengthINI));
                MarkerCorrectionMask other = new MarkerCorrectionMask(TreatedCell.CellIndex, new int[] { 1 });
                int index;
                if ((index = ParticleMask.FindIndex(other.Compare)) == -1)
                {
                    ParticleMask.Add(new MarkerCorrectionMask(TreatedCell.CellIndex));
                    ParticleMask.Last().MovedPoints[signKey]++;
                }
                else ParticleMask[index].MovedPoints[signKey]++;
                return true;
            }
            else return false;
        }

        /// <summary>
        /// Reattracts all markers to the LevelSet zero Isocountour.
        /// </summary>
        /// <param name="Cells">Removed points cell.</param>
        private void ReattractMarkers(List<int> Cells)
        {
            if (Cells == null)
            {
                Cells = new List<int>();
                foreach (KeyValuePair<int,int> DictEntry in CelltoArrayindex)
                {
                    Cells.Add(DictEntry.Key);
                }
            }
            foreach (int CellIndex in Cells)
            {
                MarkersofCell CellObj = Points[CelltoArrayindex[CellIndex]];
                /*BoundingBox CellBox = new BoundingBox(SpatialDimension);
                Grid.iGeomCells.GetCellBoundingBox(CellObj.CellIndex, CellBox);
                double ERROR = CellBox.Diameter * 0.5e-1;*/
                for(int i=0;i<CellObj.ItemList[1].Count();i++)
                {
                    if (Math.Abs(CellObj.ItemList[1][i].Phi[0, 0]) > CellObj.MaxAdvectionStepLength)
                    {
                        bool CellPosSure = true;
                        if (!ProjectOnLevelSetValue(0.0, CellObj.ItemList[1][i].Coordinates, CellObj.CellIndex, ref CellPosSure))
                        {
                            MultidimensionalArray Coord = CellObj.ItemList[1][i].Coordinates;
                            Console.WriteLine("Could not reproject Marker: Coord[" + Coord[0, 0] + "|" + Coord[0, 1] + "] mit Phi:" + CellObj.ItemList[1][i].Phi[0, 0]);
                            CellObj.ItemList[1].RemoveAt(i--);
                        }
                    }
                }
                if(CellObj.ItemList[1].IsNullOrEmpty())
                {
                    List<int> Keys = CelltoArrayindex.Keys.ToList();
                    for (int k = 0; k < Keys.Count; k++)
                    {
                        if (CelltoArrayindex[Keys[k]] > CelltoArrayindex[CellIndex])
                        {
                            CelltoArrayindex[Keys[k]] = CelltoArrayindex[Keys[k]] - 1;
                        }
                    }
                    CelltoArrayindex = new Dictionary<int, int>(CelltoArrayindex);
                    Points.RemoveAt(CelltoArrayindex[CellIndex]--);
                    CelltoArrayindex.Remove(Points[CelltoArrayindex[CellIndex]].CellIndex);
                }
            }
        }

        protected override MarkersofCell NewPointsofCellObject(int CellIndex, IGridData Grid, MinimalDistanceSearchMode LevelSetCorrection)
        {
            return new MarkersofCell(CellIndex, Grid, LevelSetCorrection);
        }

        protected override MarkerCorrectionMask NewFullCorrectionMask(int CellIndex)
        {
            int ArrayIndex = CelltoArrayindex[CellIndex];
            int PointNbr = Points[ArrayIndex].ItemList[1].Count;
            return new MarkerCorrectionMask(CellIndex, new int[1] { Points[CelltoArrayindex[CellIndex]].ItemList[1].Count });
        }

        protected override MarkerCorrectionMask NewCorrectionMask(int CellIndex)
        {
            return new MarkerCorrectionMask(CellIndex);
        }

        protected override void CorrectLevelSet(int timestepNo, CellMask CellsToCorrect = null,bool ReIteration = false)
        {
            if (CelltoArrayindex.IsNullOrEmpty())
                return;
            DGField Interface_old = Interface.CloneAs();

            if (!ReIteration)
            {
                int NumberOfMarkersPerNode = ((TargetNbrofPointsPerCellPerDimension * SpatialDimension / 5) > 5) ?
                    (TargetNbrofPointsPerCellPerDimension * SpatialDimension / 5) : 5;
                if (MinimalDistanceSearch == MinimalDistanceSearchMode.GlobalTreeFilterSearch || MinimalDistanceSearch == MinimalDistanceSearchMode.TREEFILTERTEST)
                {
                    List<SingleLvlSetMarker> AllMarker = new List<SingleLvlSetMarker>();
                    foreach (MarkersofCell CellObj in Points)
                    {
                        foreach (KeyValuePair<int, List<SingleLvlSetMarker>> DictEntry in CellObj.ItemList)
                        {
                            foreach (SingleLvlSetMarker SingleMarker in DictEntry.Value)
                            {
                                if (SingleMarker.Active == true)
                                    AllMarker.Add(SingleMarker);
                            }
                        }
                    }
                    SearchTree = new Marker_K_D_TreeFilter(AllMarker, NumberOfMarkersPerNode);
                }
                else if (MinimalDistanceSearch == MinimalDistanceSearchMode.CellTreeFilterSearch)
                {
                    foreach (MarkersofCell CellObj in Points)
                    {
                        foreach (KeyValuePair<int, List<SingleLvlSetMarker>> DictEntry in CellObj.ItemList)
                        {
                            List<SingleLvlSetMarker> CellMarkers = new List<SingleLvlSetMarker>();
                            foreach (SingleLvlSetMarker SingleMarker in DictEntry.Value)
                            {
                                if (SingleMarker.Active == true)
                                    CellMarkers.Add(SingleMarker);
                            }
                            CellObj.TreeFilter = new Marker_K_D_TreeFilter(CellMarkers, NumberOfMarkersPerNode);
                        }
                    }
                }
            }

            if (MarkerCorrection == MarkerCorrectionMode.MinDistanceWithNormal)
            {
                ScalarFunctionEx Marker_LevelSetDistancewithNormalCorrection = delegate (int cell0, int Len, NodeSet Ns, MultidimensionalArray result)
                {
                    Interface_old.Evaluate(cell0, Len, Ns, result);
                    MultidimensionalArray GlobalNodes = this.Grid.GlobalNodes.GetValue_Cell(Ns, cell0, Len);
                    int c = 0;
                    BoundingBox Box = new BoundingBox(SpatialDimension);
                    for (int cell = cell0; cell < cell0 + Len; cell++, c++)
                    {

                        if (CellsToCorrect != null && !CellsToCorrect.Contains(cell))
                        {
                            continue;
                        }

                        if (CellsToCorrect == null && !NarrowBandCells.Contains(cell) && !CelltoArrayindex.ContainsKey(cell))
                        {
                            for (int i = 0; i < Ns.NoOfNodes; i++)
                            {
                                if (result[c, i] > 0.0)
                                    result[c, i] = 1;
                                else
                                    result[c, i] = -1;
                            }
                            continue;
                        }

                        if (MinimalDistanceSearch == MinimalDistanceSearchMode.FullSearch)
                        {
                            List<int> MarkerIndexList = new List<int>();
                            if (CelltoArrayindex.ContainsKey(cell))
                            {
                                Grid.iGeomCells.GetCellBoundingBox(cell, Box);
                                double Radius = Box.Diameter;
                                List<int> CellsInRadius = CellsInRadiusAroundCell(cell, Radius);
                                foreach (int j in CellsInRadius)
                                    if (CelltoArrayindex.ContainsKey(j))
                                        MarkerIndexList.Add(j);
                            }
                            else
                            {
                                bool Markerfound = false;
                                List<int> CellIndexFrontList = new List<int> { cell };
                                List<int> CellIndexDoneList = new List<int>();
                                while (Markerfound == false)
                                {
                                    foreach (int j in CellIndexFrontList.ToList())
                                    {
                                        Grid.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaVertices, out int[] ResOut, out int[] ConnectingEntities);
                                        CellIndexFrontList.Remove(j);
                                        CellIndexDoneList.Add(j);
                                        foreach (int i in ResOut)
                                        {
                                            if (!CellIndexDoneList.Contains(i) && !CellIndexFrontList.Contains(i))
                                                CellIndexFrontList.Add(i);
                                        }
                                    }
                                    if (CellIndexFrontList.IsNullOrEmpty())
                                        break;

                                    foreach (int i in CellIndexFrontList)
                                    {
                                        if (CelltoArrayindex.ContainsKey(i))
                                        {
                                            MarkerIndexList.Add(i);
                                            Markerfound = true;
                                            break;
                                        }                                    
                                    }                                   
                                    if (Markerfound == true)
                                        break;
                                }
                                MultidimensionalArray MaxDist = MultidimensionalArray.Create(1, SpatialDimension);
                                BoundingBox CellBox = new BoundingBox(SpatialDimension);
                                Grid.iGeomCells.GetCellBoundingBox(cell, CellBox);
                                for (int dim = 0; dim < SpatialDimension; dim++)
                                {
                                    double dist1 = Math.Abs(CellBox.Min[dim] - Points[CelltoArrayindex[MarkerIndexList[0]]].Box.Max[dim]);
                                    double dist2 = Math.Abs(CellBox.Max[dim] - Points[CelltoArrayindex[MarkerIndexList[0]]].Box.Min[dim]);
                                    MaxDist[0, dim] = Math.Max(dist1, dist2);
                                }
                                double Radius = VectorNorm(MaxDist);
                                List<int> CellsInRadius = CellsInRadiusAroundCell(cell, Radius);

                                MarkerIndexList.Clear();
                                foreach (int Index in CellsInRadius)
                                    if (CelltoArrayindex.ContainsKey(Index))
                                        MarkerIndexList.Add(Index);

                            }

                            for (int k = 0; k < MarkerIndexList.Count(); k++)
                                MarkerIndexList[k] = CelltoArrayindex[MarkerIndexList[k]];

                            for (int i = 0; i < Ns.NoOfNodes; i++)
                            {
                                double Phi_abs = double.MaxValue;
                                MultidimensionalArray NearestPassivePoint = null;
                                int DeactivationSign = 0;
                                MultidimensionalArray NearestPoint = null;
                                MultidimensionalArray NearestPointNormal = null;

                                foreach (int index in MarkerIndexList)
                                {
                                    if (MinimalDistanceSearch == MinimalDistanceSearchMode.FullSearch)
                                    {
                                        foreach (KeyValuePair<int, List<SingleLvlSetMarker>> DictEntry in Points[index].ItemList)
                                        {
                                            if (DictEntry.Value.IsNullOrEmpty())
                                                throw new Exception("No Markers available for GlobalNode in Full Search");
                                            foreach (SingleLvlSetMarker Marker in DictEntry.Value)
                                            {
                                                if (Marker.Active)
                                                {
                                                    double dist = 0;
                                                    for (int dim = 0; dim < SpatialDimension; dim++)
                                                        dist += (GlobalNodes[c, i, dim] - Marker.Coordinates[0, dim]).Pow2();

                                                    dist = Math.Sqrt(dist);
                                                    if (Phi_abs > dist)
                                                    {
                                                        Phi_abs = dist;
                                                        NearestPoint = Marker.Coordinates;
                                                        NearestPointNormal = Marker.Normal;
                                                    }
                                                }
                                                else
                                                {
                                                    double dist = 0;
                                                    for (int dim = 0; dim < SpatialDimension; dim++)
                                                        dist += (GlobalNodes[c, i, dim] - Marker.Coordinates[0, dim]).Pow2();

                                                    dist = Math.Sqrt(dist);
                                                    if (Phi_abs > dist)
                                                    {
                                                        NearestPassivePoint = Marker.Coordinates;
                                                        DeactivationSign = Marker.DeactivationSign;
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }

                                if (NearestPoint != null && NearestPassivePoint == null)
                                {
                                    MultidimensionalArray VectorNodetoPoint = MultidimensionalArray.Create(1, SpatialDimension);
                                    for (int dim = 0; dim < SpatialDimension; dim++)
                                        VectorNodetoPoint[0, dim] = GlobalNodes[c, i, dim] - NearestPoint[0, dim];
                                    int sign = Math.Sign(DotProduct(VectorNodetoPoint, NearestPointNormal));
                                    result[c, i] = Phi_abs * sign;
                                }
                                else if (NearestPoint != null && NearestPassivePoint != null)
                                {
                                    int sign;
                                    MultidimensionalArray VectorNodetoActivePoint = MultidimensionalArray.Create(1, SpatialDimension);
                                    MultidimensionalArray VectorNodetoPassivePoint = MultidimensionalArray.Create(1, SpatialDimension);
                                    for (int dim = 0; dim < SpatialDimension; dim++)
                                    {
                                        VectorNodetoActivePoint[0, dim] = GlobalNodes[c, i, dim] - NearestPoint[0, dim];
                                        VectorNodetoPassivePoint[0, dim] = GlobalNodes[c, i, dim] - NearestPassivePoint[0, dim];
                                    }
                                    sign = Math.Sign(DotProduct(VectorNodetoActivePoint, NearestPointNormal));

                                    double distPassive = VectorNorm(VectorNodetoPassivePoint);
                                    if (Phi_abs > distPassive)
                                        sign = DeactivationSign;
                                    else if (DotProduct(VectorNodetoActivePoint, VectorNodetoPassivePoint) < 0)
                                        sign = DeactivationSign;

                                    result[c, i] = Phi_abs * sign;
                                }
                                else if (NearestPoint == null && NearestPassivePoint != null)
                                {
                                    result[c, i] = DeactivationSign;
                                }
                            }
                            continue;
                        }
                        else if (MinimalDistanceSearch == MinimalDistanceSearchMode.CellTreeFilterSearch)
                        {
                            List<int> MarkerIndexList = new List<int>();
                            if (CelltoArrayindex.ContainsKey(cell))
                            {
                                Grid.iGeomCells.GetCellBoundingBox(cell, Box);
                                double Radius = Box.Diameter;
                                List<int> CellsInRadius = CellsInRadiusAroundCell(cell, Radius);
                                foreach (int j in CellsInRadius)
                                    if (CelltoArrayindex.ContainsKey(j))
                                        MarkerIndexList.Add(j);
                            }
                            else
                            {
                                bool Markerfound = false;
                                List<int> indexFront = new List<int> { cell };
                                List<int> indexDone = new List<int>();
                                while (Markerfound == false)
                                {
                                    foreach (int j in indexFront.ToList())
                                    {
                                        Grid.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaVertices, out int[] NeighbourIndex, out int[] ConnectingEntities);
                                        indexDone.Add(j);
                                        indexFront.Remove(j);
                                        foreach (int i in NeighbourIndex)
                                        {
                                            if (CelltoArrayindex.ContainsKey(i))
                                            {
                                                MarkerIndexList.Add(i);
                                                Markerfound = true;
                                                break;
                                            }
                                            if (!indexFront.Contains(i) && !indexDone.Contains(i))
                                                indexFront.Add(i);
                                        }
                                        if (Markerfound == true) break;
                                    }
                                }
                                //Console.WriteLine("Markers found in cell: " + MarkerIndexList[0]);
                                MultidimensionalArray MaxDist = MultidimensionalArray.Create(1, SpatialDimension);
                                for (int dim = 0; dim < SpatialDimension; dim++)
                                {
                                    double dist1 = Points[CelltoArrayindex[MarkerIndexList[0]]].Box.Min[dim] - Points[CelltoArrayindex[MarkerIndexList[0]]].Box.Max[dim];
                                    double dist2 = Points[CelltoArrayindex[MarkerIndexList[0]]].Box.Max[dim] - Points[CelltoArrayindex[MarkerIndexList[0]]].Box.Min[dim];
                                    MaxDist[0, dim] = Math.Max(dist1, dist2);
                                }
                                double Radius = VectorNorm(MaxDist);
                                //Console.WriteLine("Radius to MarkerIndex: " + Radius);
                                List<int> CellsInRadius = CellsInRadiusAroundCell(cell, Radius);

                                MarkerIndexList.Clear();
                                foreach (int Index in CellsInRadius)
                                    if (CelltoArrayindex.ContainsKey(Index))
                                        MarkerIndexList.Add(Index);

                            }

                            for (int k = 0; k < MarkerIndexList.Count(); k++)
                                MarkerIndexList[k] = CelltoArrayindex[MarkerIndexList[k]];

                            for (int i = 0; i < Ns.NoOfNodes; i++)
                            {
                                double Phi_abs = double.MaxValue;
                                MultidimensionalArray NearestPoint = null;
                                MultidimensionalArray NearestPointNormal = null;
                                List<SingleLvlSetMarker> NearestMarkers = new List<SingleLvlSetMarker>();
                                foreach (int index in MarkerIndexList)
                                {
                                    NearestMarkers.Add(Points[index].TreeFilter.FindNearestNeighbour(GlobalNodes.ExtractSubArrayShallow(c, i, -1)));
                                }
                                foreach (SingleLvlSetMarker Marker in NearestMarkers)
                                {
                                    double dist = 0;
                                    for (int dim = 0; dim < SpatialDimension; dim++)
                                        dist += (GlobalNodes[c, i, dim] - Marker.Coordinates[0, dim]).Pow2();

                                    dist = Math.Sqrt(dist);
                                    if (Phi_abs > dist)
                                    {
                                        Phi_abs = dist;
                                        NearestPoint = Marker.Coordinates;
                                        NearestPointNormal = Marker.Normal;
                                    }
                                }

                            if (NearestPoint != null)
                                {
                                    MultidimensionalArray VectorNodetoPoint = MultidimensionalArray.Create(1, SpatialDimension);
                                    for (int dim = 0; dim < SpatialDimension; dim++)
                                        VectorNodetoPoint[0, dim] = GlobalNodes[c, i, dim] - NearestPoint[0, dim];
                                    int sign = Math.Sign(DotProduct(VectorNodetoPoint, NearestPointNormal));
                                    result[c, i] = Phi_abs * sign;
                                    if (Phi_abs > 4)
                                        Console.WriteLine("CellIndex with Phi_abs > 4: " + cell);
                                    //Console.WriteLine("New Node Value: " + result[c, i]);
                                }
                            }
                            continue;
                        }
                        else if (MinimalDistanceSearch == MinimalDistanceSearchMode.GlobalTreeFilterSearch || MinimalDistanceSearch == MinimalDistanceSearchMode.TREEFILTERTEST)
                        {
                            SingleLvlSetMarker NearestMarker;

                            for (int i = 0; i < Ns.NoOfNodes; i++)
                            {
                                NearestMarker = SearchTree.FindNearestNeighbour(GlobalNodes.ExtractSubArrayShallow(c, i, -1));

                                if (MinimalDistanceSearch == MinimalDistanceSearchMode.TREEFILTERTEST)
                                {
                                    double Phi_abs = double.MaxValue;
                                    double distance;
                                    SingleLvlSetMarker NearestPointFull = null;
                                    foreach (MarkersofCell CellObj in Points)
                                    {
                                        foreach (KeyValuePair<int, List<SingleLvlSetMarker>> DictEntry in CellObj.ItemList)
                                        {
                                            foreach (SingleLvlSetMarker Marker in DictEntry.Value)
                                            {
                                                distance = 0;
                                                for (int dim = 0; dim < SpatialDimension; dim++)
                                                {
                                                    distance += (GlobalNodes[c, i, dim] - Marker.Coordinates[0, dim]).Pow2();
                                                }
                                                distance = Math.Sqrt(distance);
                                                if (Phi_abs > distance)
                                                {
                                                    Phi_abs = distance;
                                                    NearestPointFull = Marker;
                                                }
                                            }
                                        }
                                    }

                                    double distTree = 0, distFull = 0;
                                    for (int dim = 0; dim < SpatialDimension; dim++)
                                    {
                                        distTree += (GlobalNodes[c, i, dim] - NearestMarker.Coordinates[0, dim]).Pow2();
                                        distFull += (GlobalNodes[c, i, dim] - NearestPointFull.Coordinates[0, dim]).Pow2();
                                    }
                                    distTree = Math.Sqrt(distTree);
                                    distFull = Math.Sqrt(distFull);

                                    if (NearestMarker != NearestPointFull)
                                        throw new Exception("Test failed");
                                    NearestMarker = SearchTree.FindNearestNeighbour(GlobalNodes.ExtractSubArrayShallow(c, i, -1));
                                }


                                MultidimensionalArray VectorNodetoPoint = MultidimensionalArray.Create(1, SpatialDimension);
                                for (int dim = 0; dim < SpatialDimension; dim++)
                                    VectorNodetoPoint[0, dim] = GlobalNodes[c, i, dim] - NearestMarker.Coordinates[0, dim];

                                int sign = Math.Sign(DotProduct(VectorNodetoPoint, NearestMarker.Normal));
                                double dist = 0;
                                for (int dim = 0; dim < SpatialDimension; dim++)
                                    dist += (GlobalNodes[c, i, dim] - NearestMarker.Coordinates[0, dim]).Pow2();

                                dist = Math.Sqrt(dist);
                                result[c, i] = dist * sign;
                            }
                            continue;
                        }
                        Console.WriteLine("Cell still there: " + cell);
                    }
                    
                };
                Interface.ProjectField(Marker_LevelSetDistancewithNormalCorrection);
            }
            else if (MarkerCorrection == MarkerCorrectionMode.MinDistanceWithSignUpdate)
            {
                throw new NotImplementedException();
                SinglePhaseField Interface_Unchanged = Interface.CloneAs();
                ScalarFunctionEx Marker_LevelSetAbsDistanceCorrection = delegate (int cell0, int Len, NodeSet Ns, MultidimensionalArray result)
                {
                    Interface_Unchanged.Evaluate(cell0, Len, Ns, result);
                    MultidimensionalArray GlobalNodes = this.Grid.GlobalNodes.GetValue_Cell(Ns, cell0, Len);
                    int c = 0;
                    List<int> MarkersPerNode = new List<int>();
                    int MarkerCounter = 0;
                    for (int cell = cell0; cell < cell0 + Len; cell++, c++)
                    {
                        if (MinimalDistanceSearch == MinimalDistanceSearchMode.FullSearch)
                        {
                            if (NarrowBandCells.Contains(cell))
                            {
                                List<int> MarkerIndexList = new List<int>();
                                if (CelltoArrayindex.ContainsKey(cell))
                                {
                                    MarkerIndexList.Add(cell);
                                    Grid.GetCellNeighbours(cell, GetCellNeighbours_Mode.ViaVertices, out int[] NeighbourIndex, out int[] ConnectingEntities);
                                    foreach (int j in NeighbourIndex) if (CelltoArrayindex.ContainsKey(j)) MarkerIndexList.Add(j);
                                }
                                else
                                {
                                    bool Markerfound = false;
                                    List<int> indexFront = new List<int> { cell };
                                    List<int> indexDone = new List<int>();
                                    while (Markerfound == false)
                                    {
                                        foreach (int j in indexFront.ToList())
                                        {
                                            Grid.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaVertices, out int[] NeighbourIndex, out int[] ConnectingEntities);
                                            indexDone.Add(j);
                                            indexFront.Remove(j);
                                            foreach (int i in NeighbourIndex)
                                            {
                                                if (CelltoArrayindex.ContainsKey(i))
                                                {
                                                    MarkerIndexList.Add(i);
                                                    Markerfound = true;
                                                    break;
                                                }
                                                if (!indexFront.Contains(i) && !indexDone.Contains(i))
                                                    indexFront.Add(i);
                                            }
                                            if (Markerfound == true) break;
                                        }
                                    }
                                    MultidimensionalArray MaxDist = MultidimensionalArray.Create(1, SpatialDimension);
                                    for (int dim = 0; dim < SpatialDimension; dim++)
                                    {
                                        double dist1 = Points[CelltoArrayindex[MarkerIndexList[0]]].Box.Min[dim] - Points[CelltoArrayindex[MarkerIndexList[0]]].Box.Max[dim];
                                        double dist2 = Points[CelltoArrayindex[MarkerIndexList[0]]].Box.Max[dim] - Points[CelltoArrayindex[MarkerIndexList[0]]].Box.Min[dim];
                                        MaxDist[0, dim] = Math.Max(dist1, dist2);
                                    }
                                    double Radius = VectorNorm(MaxDist);
                                    List<int> CellsInRadius = CellsInRadiusAroundCell(cell, Radius);

                                    MarkerIndexList.Clear();
                                    foreach (int Index in CellsInRadius)
                                        if (CelltoArrayindex.ContainsKey(Index))
                                            MarkerIndexList.Add(Index);

                                }

                                for (int k = 0; k < MarkerIndexList.Count(); k++)
                                    MarkerIndexList[k] = CelltoArrayindex[MarkerIndexList[k]];

                                for (int i = 0; i < Ns.NoOfNodes; i++)
                                {
                                    double Phi_abs = double.MaxValue;
                                    MultidimensionalArray NearestPoint = null;
                                    MultidimensionalArray NearestPointNormal = null;
                                    MarkerCounter = 0;
                                    foreach (int index in MarkerIndexList)
                                    {
                                        if (MinimalDistanceSearch == MinimalDistanceSearchMode.FullSearch)
                                        {
                                            foreach (KeyValuePair<int, List<SingleLvlSetMarker>> DictEntry in Points[index].ItemList)
                                            {
                                                if (DictEntry.Value.IsNullOrEmpty())
                                                    throw new Exception("No Markers available for GlobalNode in Full Search");
                                                foreach (SingleLvlSetMarker Marker in DictEntry.Value)
                                                {
                                                    if (Marker.Active)
                                                    {
                                                        double dist = 0;
                                                        for (int dim = 0; dim < SpatialDimension; dim++)
                                                            dist += (GlobalNodes[c, i, dim] - Marker.Coordinates[0, dim]).Pow2();

                                                        MarkerCounter++;
                                                        dist = Math.Sqrt(dist);
                                                        if (Phi_abs > dist)
                                                        {
                                                            Phi_abs = dist;
                                                        }
                                                    }
                                                }
                                            }
                                        }

                                        if (Phi_abs != double.MaxValue)
                                        {
                                            double sign = Math.Sign(result[c, i]);
                                            result[c, i] = Phi_abs * sign;
                                        }
                                    }
                                }
                            }
                            else
                            {
                                for (int i = 0; i < Ns.NoOfNodes; i++)
                                {
                                    if (result[c, i] > 0.0)
                                        result[c, i] = 1;
                                    else
                                        result[c, i] = -1;
                                }
                            }
                        }
                        else if (MinimalDistanceSearch == MinimalDistanceSearchMode.CellTreeFilterSearch)
                        {
                            if (NarrowBandCells.Contains(cell))
                            {
                                List<int> MarkerIndexList = new List<int>();
                                if (CelltoArrayindex.ContainsKey(cell))
                                {
                                    MarkerIndexList.Add(cell);
                                    Grid.GetCellNeighbours(cell, GetCellNeighbours_Mode.ViaVertices, out int[] NeighbourIndex, out int[] ConnectingEntities);
                                    foreach (int j in NeighbourIndex) if (CelltoArrayindex.ContainsKey(j)) MarkerIndexList.Add(j);
                                }
                                else
                                {
                                    bool Markerfound = false;
                                    List<int> indexFront = new List<int> { cell };
                                    List<int> indexDone = new List<int>();
                                    while (Markerfound == false)
                                    {
                                        foreach (int j in indexFront.ToList())
                                        {
                                            Grid.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaVertices, out int[] NeighbourIndex, out int[] ConnectingEntities);
                                            indexDone.Add(j);
                                            indexFront.Remove(j);
                                            foreach (int i in NeighbourIndex)
                                            {
                                                if (CelltoArrayindex.ContainsKey(i))
                                                {
                                                    MarkerIndexList.Add(i);
                                                    Markerfound = true;
                                                    break;
                                                }
                                                if (!indexFront.Contains(i) && !indexDone.Contains(i))
                                                    indexFront.Add(i);
                                            }
                                            if (Markerfound == true) break;
                                        }
                                    }
                                    MultidimensionalArray MaxDist = MultidimensionalArray.Create(1, SpatialDimension);
                                    for (int dim = 0; dim < SpatialDimension; dim++)
                                    {
                                        double dist1 = Points[CelltoArrayindex[MarkerIndexList[0]]].Box.Min[dim] - Points[CelltoArrayindex[MarkerIndexList[0]]].Box.Max[dim];
                                        double dist2 = Points[CelltoArrayindex[MarkerIndexList[0]]].Box.Max[dim] - Points[CelltoArrayindex[MarkerIndexList[0]]].Box.Min[dim];
                                        MaxDist[0, dim] = Math.Max(dist1, dist2);
                                    }
                                    double Radius = VectorNorm(MaxDist);
                                    List<int> CellsInRadius = CellsInRadiusAroundCell(cell, Radius);

                                    MarkerIndexList.Clear();
                                    foreach (int Index in CellsInRadius)
                                        if (CelltoArrayindex.ContainsKey(Index))
                                            MarkerIndexList.Add(Index);

                                }

                                for (int k = 0; k < MarkerIndexList.Count(); k++)
                                    MarkerIndexList[k] = CelltoArrayindex[MarkerIndexList[k]];

                                for (int i = 0; i < Ns.NoOfNodes; i++)
                                {
                                    double Phi_abs = double.MaxValue;
                                    MultidimensionalArray NearestPoint = null;
                                    MultidimensionalArray NearestPointNormal = null;
                                    MarkerCounter = 0;
                                    List<SingleLvlSetMarker> NearestMarkers = new List<SingleLvlSetMarker>();
                                    foreach (int index in MarkerIndexList)
                                    {
                                        NearestMarkers.Add(Points[index].TreeFilter.FindNearestNeighbour(GlobalNodes.ExtractSubArrayShallow(c, i, -1)));
                                    }
                                    foreach (SingleLvlSetMarker Marker in NearestMarkers)
                                    {
                                        double dist = 0;
                                        for (int dim = 0; dim < SpatialDimension; dim++)
                                            dist += (GlobalNodes[c, i, dim] - Marker.Coordinates[0, dim]).Pow2();

                                        MarkerCounter++;
                                        dist = Math.Sqrt(dist);
                                        if (Phi_abs > dist)
                                        {
                                            Phi_abs = dist;
                                        }
                                    }

                                    if (Phi_abs != double.MaxValue)
                                    {
                                        if (Phi_abs != double.MaxValue)
                                        {
                                            double sign = Math.Sign(result[c, i]);
                                            result[c, i] = Phi_abs * sign;
                                        }
                                    }
                                }
                            }
                            else
                            {
                                for (int i = 0; i < Ns.NoOfNodes; i++)
                                {
                                    if (result[c, i] > 0.0)
                                        result[c, i] = 1;
                                    else
                                        result[c, i] = -1;
                                }
                            }
                        }
                        else if (MinimalDistanceSearch == MinimalDistanceSearchMode.GlobalTreeFilterSearch || MinimalDistanceSearch == MinimalDistanceSearchMode.TREEFILTERTEST)
                        {
                            if (NarrowBandCells.Contains(cell))
                            {
                                SingleLvlSetMarker NearestMarker;

                                for (int i = 0; i < Ns.NoOfNodes; i++)
                                {
                                    NearestMarker = SearchTree.FindNearestNeighbour(GlobalNodes.ExtractSubArrayShallow(c, i, -1));

                                    if (MinimalDistanceSearch == MinimalDistanceSearchMode.TREEFILTERTEST)
                                    {
                                        double Phi_abs = double.MaxValue;
                                        double distance;
                                        SingleLvlSetMarker NearestPointFull = null;
                                        foreach (MarkersofCell CellObj in Points)
                                        {
                                            foreach (KeyValuePair<int, List<SingleLvlSetMarker>> DictEntry in CellObj.ItemList)
                                            {
                                                foreach (SingleLvlSetMarker Marker in DictEntry.Value)
                                                {
                                                    distance = 0;
                                                    for (int dim = 0; dim < SpatialDimension; dim++)
                                                    {
                                                        distance += (GlobalNodes[c, i, dim] - Marker.Coordinates[0, dim]).Pow2();
                                                    }
                                                    distance = Math.Sqrt(distance);
                                                    if (Phi_abs > distance)
                                                    {
                                                        Phi_abs = distance;
                                                        NearestPointFull = Marker;
                                                    }
                                                }
                                            }
                                        }

                                        double distTree = 0, distFull = 0;
                                        for (int dim = 0; dim < SpatialDimension; dim++)
                                        {
                                            distTree += (GlobalNodes[c, i, dim] - NearestMarker.Coordinates[0, dim]).Pow2();
                                            distFull += (GlobalNodes[c, i, dim] - NearestPointFull.Coordinates[0, dim]).Pow2();
                                        }
                                        distTree = Math.Sqrt(distTree);
                                        distFull = Math.Sqrt(distFull);

                                        if (NearestMarker != NearestPointFull)
                                            throw new Exception("Test failed");
                                        NearestMarker = SearchTree.FindNearestNeighbour(GlobalNodes.ExtractSubArrayShallow(c, i, -1));
                                    }

                                    int sign = Math.Sign(result[c, i]);
                                    double dist = 0;
                                    for (int dim = 0; dim < SpatialDimension; dim++)
                                        dist += (GlobalNodes[c, i, dim] - NearestMarker.Coordinates[0, dim]).Pow2();
                                    dist = Math.Sqrt(dist);

                                    MarkerCounter++;
                                    result[c, i] = dist * sign;
                                }
                            }
                            else
                            {
                                for (int i = 0; i < Ns.NoOfNodes; i++)
                                {
                                    if (result[c, i] > 0.0)
                                        result[c, i] = 1;
                                    else
                                        result[c, i] = -1;
                                }
                            }
                        }
                    }
                };
                Interface.ProjectField(Marker_LevelSetAbsDistanceCorrection);
                Console.WriteLine("Distance corrected");

                InterfaceLevelSetGradient.Gradient(1, Interface, NarrowBandCells);
                SinglePhaseField Interface_UpdtAbsDist = Interface.CloneAs();
                ScalarFunctionEx Marker_LevelSetSignCorrection = delegate (int cell0, int Len, NodeSet Ns, MultidimensionalArray result)
                {
                    Interface_UpdtAbsDist.Evaluate(cell0, Len, Ns, result);
                    MultidimensionalArray result_unchanged = result.CloneAs();
                    Interface_Unchanged.Evaluate(cell0, Len, Ns, result_unchanged);
                    MultidimensionalArray GlobalNodes = this.Grid.GlobalNodes.GetValue_Cell(Ns, cell0, Len);
                    int c = 0;
                    BoundingBox Box = new BoundingBox(SpatialDimension);
                    for (int cell = cell0; cell < cell0 + Len; cell++, c++)
                    {
                        double MaxAdvectionLength = double.MinValue;
                        double delta_X;
                        if (NarrowBandCells.Contains(cell))
                        {
                            Console.Write(cell);
                            Grid.iGeomCells.GetCellBoundingBox(cell, Box);
                            delta_X = Box.h_min;
                            if (CelltoArrayindex.ContainsKey(cell))
                                MaxAdvectionLength = Points[CelltoArrayindex[cell]].MaxAdvectionStepLength;
                            Grid.GetCellNeighbours(cell, GetCellNeighbours_Mode.ViaVertices, out int[] NeighbourIndex, out int[] ConnectingEntities);
                            foreach (int j in NeighbourIndex)
                            {
                                if (CelltoArrayindex.ContainsKey(j))
                                    if (MaxAdvectionLength < Points[CelltoArrayindex[j]].MaxAdvectionStepLength)
                                        MaxAdvectionLength = Points[CelltoArrayindex[j]].MaxAdvectionStepLength;
                                Grid.iGeomCells.GetCellBoundingBox(j, Box);
                                if (Box.h_min < delta_X)
                                    delta_X = Box.h_min;
                            }
                            if (MaxAdvectionLength == double.MinValue)
                                continue;
                            double distChange = double.MinValue;
                            for (int i = 0; i < Ns.NoOfNodes; i++)
                            {
                                distChange = Math.Abs(result[c, i] - result_unchanged[c, i]);
                                MaxAdvectionLength = (MaxAdvectionLength < distChange) ? distChange : MaxAdvectionLength;
                            }
                            int CountNodes = 0;
                            List<int> NodesNearInterface = new List<int>();
                            for (int i = 0; i < Ns.NoOfNodes; i++)
                            {
                                if (Math.Abs(result[c, i]) <= 1.1 * MaxAdvectionLength)
                                {
                                    CountNodes++;
                                    NodesNearInterface.Add(i);
                                }
                            }
                            if (CountNodes == 0)
                                continue;

                            Console.WriteLine("---" + CountNodes);
                            MultidimensionalArray Nodes = MultidimensionalArray.Create(CountNodes, SpatialDimension);
                            int index = 0;
                            foreach (int i in NodesNearInterface)
                            {
                                for (int dim = 0; dim < SpatialDimension; dim++)
                                {
                                    Nodes[index, dim] = GlobalNodes[c, i, dim];
                                }
                                index++;
                            }
                            bool CellPosSure = true;
                            DGField[] Fields = new DGField[SpatialDimension];
                            for (int dim = 0; dim < SpatialDimension; dim++)
                                Fields[dim] = InterfaceLevelSetGradient[dim];
                            MultidimensionalArray NodeNormals = MultiEvaluatePointsInCell(Fields, Nodes, cell, CellPosSure);
                            for (int inde = 0; inde < NodeNormals.GetLength(1); inde++)
                            {
                                MultidimensionalArray Normal = NormalizeVector(NodeNormals.ExtractSubArrayShallow(-1, inde));
                                for (int dim = 0; dim < SpatialDimension; dim++)
                                    NodeNormals[dim, inde] = Normal[dim];
                            }
                            Dictionary<int, MultidimensionalArray> DistLevelSetPoints = new Dictionary<int, MultidimensionalArray>();
                            DistLevelSetPoints.Add(-1, MultidimensionalArray.Create(CountNodes, 1, SpatialDimension));
                            DistLevelSetPoints.Add(1, MultidimensionalArray.Create(CountNodes, 1, SpatialDimension));
                            foreach (KeyValuePair<int, MultidimensionalArray> DictEntry in DistLevelSetPoints)
                            {
                                for (int ze = 0; ze < CountNodes; ze++)
                                {
                                    for (int dim = 0; dim < SpatialDimension; dim++)
                                    {
                                        DictEntry.Value[ze, 0, dim] = Nodes[ze, dim] + DictEntry.Key * 1.5 * MaxAdvectionLength * NodeNormals[dim, ze];
                                    }
                                }
                            }
                            BitArray bitArray = new BitArray(CountNodes);
                            bitArray.SetAll(true);
                            foreach (KeyValuePair<int, MultidimensionalArray> DictEntry in DistLevelSetPoints)
                            {
                                for (int zi = 0; zi < CountNodes; zi++)
                                {
                                    bool CellSure = false;
                                    bool ProjectSuccesful = ProjectOnLevelSetValue(DictEntry.Key * delta_X, DictEntry.Value.ExtractSubArrayShallow(zi, -1, -1),
                                        cell, ref CellSure, Interface_UpdtAbsDist);
                                    bitArray[zi] = (bitArray[zi]) ? ProjectSuccesful : false;
                                }
                            }

                            double dist_minus = 0, dist_plus = 0;
                            int ind = 0, z = 0;
                            foreach (int i in NodesNearInterface)
                            {
                                for (int dim = 0; dim < SpatialDimension; dim++)
                                {
                                    dist_minus += (Nodes[ind, dim] - DistLevelSetPoints[-1][z, 0, dim]).Pow2();
                                    dist_plus += (Nodes[ind, dim] - DistLevelSetPoints[1][z, 0, dim]).Pow2();
                                }
                                dist_minus = Math.Sqrt(dist_minus);
                                dist_plus = Math.Sqrt(dist_plus);

                                int sign;
                                if (dist_minus < dist_plus)
                                    sign = Math.Sign(dist_minus - delta_X);
                                else
                                    sign = Math.Sign(delta_X - dist_plus);

                                result[c, i] = sign * Math.Abs(result[c, i]);

                                ind++;
                                z++;
                            }
                        }
                    }
                };
                Interface.ProjectField(Marker_LevelSetSignCorrection);
                Console.WriteLine("Sign corrected");
            }
            else
            {
                throw new NotImplementedException();
                OldNearestPointMemory NewNearestPointMemory = new OldNearestPointMemory();

                SinglePhaseField Interface_Unchanged = Interface.CloneAs();
                ScalarFunctionEx Marker_LevelSetAbsDistanceCorrection = delegate (int cell0, int Len, NodeSet Ns, MultidimensionalArray result)
                {
                    Interface_Unchanged.Evaluate(cell0, Len, Ns, result);
                    MultidimensionalArray GlobalNodes = this.Grid.GlobalNodes.GetValue_Cell(Ns, cell0, Len);
                    int c = 0;
                    List<int> MarkersPerNode = new List<int>();
                    int MarkerCounter = 0;
                    for (int cell = cell0; cell < cell0 + Len; cell++, c++)
                    {
                        if (MinimalDistanceSearch == MinimalDistanceSearchMode.FullSearch)
                        {
                            if (NarrowBandCells.Contains(cell))
                            {
                                List<int> MarkerIndexList = new List<int>();
                                if (CelltoArrayindex.ContainsKey(cell))
                                {
                                    MarkerIndexList.Add(cell);
                                    Grid.GetCellNeighbours(cell, GetCellNeighbours_Mode.ViaVertices, out int[] NeighbourIndex, out int[] ConnectingEntities);
                                    foreach (int j in NeighbourIndex) if (CelltoArrayindex.ContainsKey(j)) MarkerIndexList.Add(j);
                                }
                                else
                                {
                                    bool Markerfound = false;
                                    List<int> indexFront = new List<int> { cell };
                                    List<int> indexDone = new List<int>();
                                    while (Markerfound == false)
                                    {
                                        foreach (int j in indexFront.ToList())
                                        {
                                            Grid.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaVertices, out int[] NeighbourIndex, out int[] ConnectingEntities);
                                            indexDone.Add(j);
                                            indexFront.Remove(j);
                                            foreach (int i in NeighbourIndex)
                                            {
                                                if (CelltoArrayindex.ContainsKey(i))
                                                {
                                                    MarkerIndexList.Add(i);
                                                    Markerfound = true;
                                                    break;
                                                }
                                                if (!indexFront.Contains(i) && !indexDone.Contains(i))
                                                    indexFront.Add(i);
                                            }
                                            if (Markerfound == true) break;
                                        }
                                    }
                                    MultidimensionalArray MaxDist = MultidimensionalArray.Create(1, SpatialDimension);
                                    for (int dim = 0; dim < SpatialDimension; dim++)
                                    {
                                        double dist1 = Points[CelltoArrayindex[MarkerIndexList[0]]].Box.Min[dim] - Points[CelltoArrayindex[MarkerIndexList[0]]].Box.Max[dim];
                                        double dist2 = Points[CelltoArrayindex[MarkerIndexList[0]]].Box.Max[dim] - Points[CelltoArrayindex[MarkerIndexList[0]]].Box.Min[dim];
                                        MaxDist[0, dim] = Math.Max(dist1, dist2);
                                    }
                                    double Radius = VectorNorm(MaxDist);
                                    List<int> CellsInRadius = CellsInRadiusAroundCell(cell, Radius);

                                    MarkerIndexList.Clear();
                                    foreach (int Index in CellsInRadius)
                                        if (CelltoArrayindex.ContainsKey(Index))
                                            MarkerIndexList.Add(Index);

                                }

                                for (int k = 0; k < MarkerIndexList.Count(); k++)
                                    MarkerIndexList[k] = CelltoArrayindex[MarkerIndexList[k]];

                                for (int i = 0; i < Ns.NoOfNodes; i++)
                                {
                                    double Phi_abs = double.MaxValue;
                                    SingleLvlSetMarker NearestMarker = null;
                                    MarkerCounter = 0;
                                    foreach (int index in MarkerIndexList)
                                    {
                                        if (MinimalDistanceSearch == MinimalDistanceSearchMode.FullSearch)
                                        {
                                            foreach (KeyValuePair<int, List<SingleLvlSetMarker>> DictEntry in Points[index].ItemList)
                                            {
                                                if (DictEntry.Value.IsNullOrEmpty())
                                                    throw new Exception("No Markers available for GlobalNode in Full Search");
                                                foreach (SingleLvlSetMarker Marker in DictEntry.Value)
                                                {
                                                    if (Marker.Active)
                                                    {
                                                        double dist = 0;
                                                        for (int dim = 0; dim < SpatialDimension; dim++)
                                                            dist += (GlobalNodes[c, i, dim] - Marker.Coordinates[0, dim]).Pow2();

                                                        MarkerCounter++;
                                                        dist = Math.Sqrt(dist);
                                                        if (Phi_abs > dist)
                                                        {
                                                            Phi_abs = dist;
                                                            NearestMarker = Marker;
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                        if (Phi_abs != double.MaxValue)
                                        {
                                            double sign = Math.Sign(result[c, i]);
                                            result[c, i] = Phi_abs * sign;
                                            NewNearestPointMemory.AddNode(cell, i, NearestMarker);
                                        }
                                    }
                                }
                            }
                            else
                            {
                                for (int i = 0; i < Ns.NoOfNodes; i++)
                                {
                                    if (result[c, i] > 0.0)
                                        result[c, i] = 1;
                                    else
                                        result[c, i] = -1;
                                }
                            }
                        }
                        else if (MinimalDistanceSearch == MinimalDistanceSearchMode.CellTreeFilterSearch)
                        {
                            if (NarrowBandCells.Contains(cell))
                            {
                                //FindCells.Start();
                                List<int> MarkerIndexList = new List<int>();
                                if (CelltoArrayindex.ContainsKey(cell))
                                {
                                    MarkerIndexList.Add(cell);
                                    Grid.GetCellNeighbours(cell, GetCellNeighbours_Mode.ViaVertices, out int[] NeighbourIndex, out int[] ConnectingEntities);
                                    foreach (int j in NeighbourIndex) if (CelltoArrayindex.ContainsKey(j)) MarkerIndexList.Add(j);
                                }
                                else
                                {
                                    bool Markerfound = false;
                                    List<int> indexFront = new List<int> { cell };
                                    List<int> indexDone = new List<int>();
                                    while (Markerfound == false)
                                    {
                                        foreach (int j in indexFront.ToList())
                                        {
                                            Grid.GetCellNeighbours(j, GetCellNeighbours_Mode.ViaVertices, out int[] NeighbourIndex, out int[] ConnectingEntities);
                                            indexDone.Add(j);
                                            indexFront.Remove(j);
                                            foreach (int i in NeighbourIndex)
                                            {
                                                if (CelltoArrayindex.ContainsKey(i))
                                                {
                                                    MarkerIndexList.Add(i);
                                                    Markerfound = true;
                                                    break;
                                                }
                                                if (!indexFront.Contains(i) && !indexDone.Contains(i))
                                                    indexFront.Add(i);
                                            }
                                            if (Markerfound == true) break;
                                        }
                                    }
                                    //Console.WriteLine("Markers found in cell: " + MarkerIndexList[0]);
                                    MultidimensionalArray MaxDist = MultidimensionalArray.Create(1, SpatialDimension);
                                    for (int dim = 0; dim < SpatialDimension; dim++)
                                    {
                                        double dist1 = Points[CelltoArrayindex[MarkerIndexList[0]]].Box.Min[dim] - Points[CelltoArrayindex[MarkerIndexList[0]]].Box.Max[dim];
                                        double dist2 = Points[CelltoArrayindex[MarkerIndexList[0]]].Box.Max[dim] - Points[CelltoArrayindex[MarkerIndexList[0]]].Box.Min[dim];
                                        MaxDist[0, dim] = Math.Max(dist1, dist2);
                                    }
                                    double Radius = VectorNorm(MaxDist);
                                    //Console.WriteLine("Radius to MarkerIndex: " + Radius);
                                    List<int> CellsInRadius = CellsInRadiusAroundCell(cell, Radius);

                                    MarkerIndexList.Clear();
                                    foreach (int Index in CellsInRadius)
                                        if (CelltoArrayindex.ContainsKey(Index))
                                            MarkerIndexList.Add(Index);

                                }

                                for (int k = 0; k < MarkerIndexList.Count(); k++)
                                    MarkerIndexList[k] = CelltoArrayindex[MarkerIndexList[k]];

                                for (int i = 0; i < Ns.NoOfNodes; i++)
                                {
                                    double Phi_abs = double.MaxValue;
                                    SingleLvlSetMarker NearestMarker = null;

                                    List<SingleLvlSetMarker> NearestMarkers = new List<SingleLvlSetMarker>();
                                    foreach (int index in MarkerIndexList)
                                    {
                                        NearestMarkers.Add(Points[index].TreeFilter.FindNearestNeighbour(GlobalNodes.ExtractSubArrayShallow(c, i, -1)));
                                    }
                                    foreach (SingleLvlSetMarker Marker in NearestMarkers)
                                    {
                                        double dist = 0;
                                        for (int dim = 0; dim < SpatialDimension; dim++)
                                            dist += (GlobalNodes[c, i, dim] - Marker.Coordinates[0, dim]).Pow2();

                                        dist = Math.Sqrt(dist);
                                        if (Phi_abs > dist)
                                        {
                                            Phi_abs = dist;
                                            NearestMarker = Marker;
                                        }
                                    }

                                    if (Phi_abs != double.MaxValue)
                                    {
                                        double sign = Math.Sign(result[c, i]);
                                        result[c, i] = Phi_abs * sign;
                                        NewNearestPointMemory.AddNode(cell, i, NearestMarker);
                                    }
                                }
                            }
                            else
                            {
                                for (int i = 0; i < Ns.NoOfNodes; i++)
                                {
                                    if (result[c, i] > 0.0)
                                        result[c, i] = 1;
                                    else
                                        result[c, i] = -1;
                                }
                            }
                        }
                        else if (MinimalDistanceSearch == MinimalDistanceSearchMode.GlobalTreeFilterSearch || MinimalDistanceSearch == MinimalDistanceSearchMode.TREEFILTERTEST)
                        {
                            if (NarrowBandCells.Contains(cell))
                            {
                                SingleLvlSetMarker NearestMarker;

                                for (int i = 0; i < Ns.NoOfNodes; i++)
                                {
                                    NearestMarker = SearchTree.FindNearestNeighbour(GlobalNodes.ExtractSubArrayShallow(c, i, -1));

                                    if (MinimalDistanceSearch == MinimalDistanceSearchMode.TREEFILTERTEST)
                                    {
                                        double Phi_abs = double.MaxValue;
                                        double distance;
                                        SingleLvlSetMarker NearestPointFull = null;
                                        foreach (MarkersofCell CellObj in Points)
                                        {
                                            foreach (KeyValuePair<int, List<SingleLvlSetMarker>> DictEntry in CellObj.ItemList)
                                            {
                                                foreach (SingleLvlSetMarker Marker in DictEntry.Value)
                                                {
                                                    distance = 0;
                                                    for (int dim = 0; dim < SpatialDimension; dim++)
                                                    {
                                                        distance += (GlobalNodes[c, i, dim] - Marker.Coordinates[0, dim]).Pow2();
                                                    }
                                                    distance = Math.Sqrt(distance);
                                                    if (Phi_abs > distance)
                                                    {
                                                        Phi_abs = distance;
                                                        NearestPointFull = Marker;
                                                    }
                                                }
                                            }
                                        }

                                        double distTree = 0, distFull = 0;
                                        for (int dim = 0; dim < SpatialDimension; dim++)
                                        {
                                            distTree += (GlobalNodes[c, i, dim] - NearestMarker.Coordinates[0, dim]).Pow2();
                                            distFull += (GlobalNodes[c, i, dim] - NearestPointFull.Coordinates[0, dim]).Pow2();
                                        }
                                        distTree = Math.Sqrt(distTree);
                                        distFull = Math.Sqrt(distFull);

                                        if (NearestMarker != NearestPointFull)
                                            throw new Exception("Test failed");
                                        NearestMarker = SearchTree.FindNearestNeighbour(GlobalNodes.ExtractSubArrayShallow(c, i, -1));
                                    }

                                    int sign = Math.Sign(result[c, i]);
                                    double dist = 0;
                                    for (int dim = 0; dim < SpatialDimension; dim++)
                                        dist += (GlobalNodes[c, i, dim] - NearestMarker.Coordinates[0, dim]).Pow2();
                                    dist = Math.Sqrt(dist);

                                    result[c, i] = dist * sign;
                                    NewNearestPointMemory.AddNode(cell, i, NearestMarker);
                                }
                            }
                            else
                            {
                                for (int i = 0; i < Ns.NoOfNodes; i++)
                                {
                                    if (result[c, i] > 0.0)
                                        result[c, i] = 1;
                                    else
                                        result[c, i] = -1;
                                }
                            }
                        }
                    }
                };
                Interface.ProjectField(Marker_LevelSetAbsDistanceCorrection);
                Console.WriteLine("Distance corrected");

                SinglePhaseField Interface_UpdtAbsDist = Interface.CloneAs();
                ScalarFunctionEx Marker_LevelSetSignCorrection = delegate (int cell0, int Len, NodeSet Ns, MultidimensionalArray result)
                {
                    Interface_UpdtAbsDist.Evaluate(cell0, Len, Ns, result);
                    MultidimensionalArray result_unchanged = result.CloneAs();
                    Interface_Unchanged.Evaluate(cell0, Len, Ns, result_unchanged);
                    MultidimensionalArray GlobalNodes = this.Grid.GlobalNodes.GetValue_Cell(Ns, cell0, Len);
                    int c = 0;
                    for (int cell = cell0; cell < cell0 + Len; cell++, c++)
                    {
                        double MaxAdvectionLength = double.MinValue;
                        if (NarrowBandCells.Contains(cell))
                        {
                            if (CelltoArrayindex.ContainsKey(cell))
                                MaxAdvectionLength = Points[CelltoArrayindex[cell]].MaxAdvectionStepLength;
                            Grid.GetCellNeighbours(cell, GetCellNeighbours_Mode.ViaVertices, out int[] NeighbourIndex, out int[] ConnectingEntities);
                            foreach (int j in NeighbourIndex)
                            {
                                if (CelltoArrayindex.ContainsKey(j))
                                    if (MaxAdvectionLength < Points[CelltoArrayindex[j]].MaxAdvectionStepLength)
                                        MaxAdvectionLength = Points[CelltoArrayindex[j]].MaxAdvectionStepLength;
                            }
                            if (MaxAdvectionLength == double.MinValue)
                                continue;

                            double distChange = double.MinValue;
                            for (int i = 0; i < Ns.NoOfNodes; i++)
                            {
                                distChange = Math.Abs(result[c, i] - result_unchanged[c, i]);
                                MaxAdvectionLength = (MaxAdvectionLength < distChange) ? distChange : MaxAdvectionLength;
                            }
                            int CountNodes = 0;
                            List<int> NodesNearInterface = new List<int>();
                            for (int i = 0; i < Ns.NoOfNodes; i++)
                            {
                                if (Math.Abs(result[c, i]) <= 1.5 * MaxAdvectionLength)
                                {
                                    CountNodes++;
                                    NodesNearInterface.Add(i);
                                }
                            }
                            if (CountNodes == 0)
                                continue;

                            MultidimensionalArray[] GlobalNode = new MultidimensionalArray[CountNodes];
                            int index = 0;
                            foreach (int i in NodesNearInterface)
                            {
                                GlobalNode[index] = MultidimensionalArray.Create(1, SpatialDimension);
                                for (int dim = 0; dim < SpatialDimension; dim++)
                                {
                                    GlobalNode[index][0, dim] = GlobalNodes[c, i, dim];
                                }
                                index++;
                            }

                            index = 0;
                            foreach (int i in NodesNearInterface)
                            {
                                MultidimensionalArray PreviousNearestMarkerCoord = PreviousNearestMarker.FindOldNearestNeighbour(cell, i);
                                MultidimensionalArray NewNearestMarkerCoord = NewNearestPointMemory.FindOldNearestNeighbour(cell, i);

                                MultidimensionalArray VectorToPrevious = ComputeVectorbetweenPoints(PreviousNearestMarkerCoord, GlobalNode[index]);
                                MultidimensionalArray VectorToNew = ComputeVectorbetweenPoints(NewNearestMarkerCoord, GlobalNode[index]);

                                bool signChange = false;
                                double signDir = DotProduct(VectorToPrevious, VectorToNew);
                                if (signDir < 0)
                                    signChange = true;
                                else if (signDir > 0)
                                    signChange = false;
                                else
                                {
                                    double PrevLen = VectorNorm(VectorToPrevious);
                                    double NewLen = VectorNorm(VectorToNew);
                                    if (NewLen == 0)
                                        signChange = false;
                                    else if (PrevLen == 0)
                                        throw new Exception("New thinking needed");
                                }

                                if (signChange)
                                    result[c, i] *= -1;
                            }
                        }
                    }
                };
                Interface.ProjectField(Marker_LevelSetSignCorrection);
                Console.WriteLine("Sign corrected");
            }

            ScalarFunctionEx Marker_LevelSetGaussianCorrection = delegate (int cell0, int Len, NodeSet Ns, MultidimensionalArray result)
            {
                Interface_old.Evaluate(cell0, Len, Ns, result);
                MultidimensionalArray GlobalNodes = this.Grid.GlobalNodes.GetValue_Cell(Ns, cell0, Len);

                int c = 0;
                for (int cell = cell0; cell < cell0 + Len; cell++, c++)
                {
                    if (NarrowBandCells.Contains(cell))
                    {
                        //Console.WriteLine("Cell Nbr: " + cell + " is corrected");

                        List<int> MarkerIndexList = new List<int>();
                        List<int> CellIndexListChecked = new List<int>(cell);
                        List<int> CellIndexListFront = new List<int>(cell);
                        if (CelltoArrayindex.ContainsKey(cell))
                        {
                            MarkerIndexList.Add(cell);
                        }
                        bool Continue = true;
                        while (Continue)
                        {
                            foreach (int cellIndex in CellIndexListFront)
                            {
                                Grid.GetCellNeighbours(cellIndex, GetCellNeighbours_Mode.ViaVertices, out int[] NeighbourIndex, out int[] ConnectingEntities);
                                CellIndexListFront.Remove(cellIndex);
                                foreach (int j in NeighbourIndex)
                                {
                                    if (!CellIndexListChecked.Contains(j))
                                    {
                                        double[] min = new double[3];
                                        double dist = 0;
                                        for (int dim = 0; dim < SpatialDimension; dim++)
                                        {
                                            min[dim] = Math.Abs(Points[CelltoArrayindex[j]].Box.Max[dim] - Points[CelltoArrayindex[cell]].Box.Min[dim]);
                                            double d = Math.Abs(Points[CelltoArrayindex[j]].Box.Min[dim] - Points[CelltoArrayindex[cell]].Box.Max[dim]);
                                            min[dim] = (d < min[dim]) ? d : min[dim];
                                            dist = min[dim].Pow2();
                                        }
                                        dist = Math.Sqrt(dist);
                                        if (dist < KERNEL_RADIUS)
                                        {
                                            CellIndexListFront.Add(j);
                                            CellIndexListChecked.Add(j);
                                            if (CelltoArrayindex.ContainsKey(j))
                                            {
                                                MarkerIndexList.Add(j);
                                            }
                                        }
                                    }
                                }
                            }
                            if (CellIndexListFront.IsNullOrEmpty()) Continue = false;
                        }

                        for (int k = 0; k < MarkerIndexList.Count(); k++)
                            MarkerIndexList[k] = CelltoArrayindex[MarkerIndexList[k]];

                        for (int i = 0; i < Ns.NoOfNodes; i++)
                        {
                            int len = 0;
                            foreach (int k in MarkerIndexList)
                            {
                                len += Points[k].ItemList[1].Count;
                            }
                            double[] w = new double[len];
                            double[] Phi_Markers = new double[len];
                            int f = 0;
                            foreach (int k in MarkerIndexList)
                            {
                                for (int j = 0; j < Points[k].ItemList[1].Count; j++)
                                {
                                    double q = 0;
                                    double dist = 0;
                                    for (int dim = 0; dim < SpatialDimension; dim++)
                                    {
                                        dist += (GlobalNodes[c, i, dim] - Points[k].ItemList[1][j].Coordinates[0, dim]).Pow2();
                                    }
                                    dist = Math.Sqrt(dist);
                                    q = dist / KERNEL_RADIUS;
                                    double c_ = 0.5;
                                    w[f] = (Math.Pow(Math.E, -(q / c_).Pow2()) - Math.Pow(Math.E, -(1 / c_.Pow2()))) / (-Math.Pow(Math.E, -(1 / c_.Pow2())));
                                    Phi_Markers[f++] = Evaluate(new DGField[] { Interface }, Points[k].ItemList[1][j].Coordinates, k, true).To1DArray()[0];
                                }
                            }
                            double sum_w = 0;
                            double sum_w_phi = 0;
                            for (int z = 0; z < Phi_Markers.Length; z++)
                            {
                                sum_w += w[z];
                                sum_w_phi += w[z] * Phi_Markers[z];
                            }
                            result[c, i] -= sum_w_phi / sum_w;
                        }
                    }
                    //else Console.WriteLine("Cell Nbr: " + cell + " stays the same");
                }
            };

            if (!ReIteration)
            {
                //InterfaceTracker.UpdateTracker();
                NarrowBandCells_Old = NarrowBandCells;
                NarrowBandCells = InterfaceTracker.Regions.GetNearFieldMask(NarrowBandWidth);
                CellMask ReIterationMask = NarrowBandCells.Except(NarrowBandCells_Old);
                CorrectLevelSet(timestepNo, ReIterationMask, true);
                InterfaceLevelSetGradient.Clear();
                InterfaceLevelSetGradient.Gradient(1, Interface, NarrowBandCells);
            }
        }
    }
}
