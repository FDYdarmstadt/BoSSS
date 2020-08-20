using System;
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
    public class VerticesofCell : SemiLagrangianLevelSetMethod
    <VerticesofCell, SingleVertice, VerticeCorrectionMask>.ItemsofCell
    {
        public VerticesofCell (int CellIndex, IGridData Grid, SemiLagrangianLevelSetMethod
            <VerticesofCell, SingleVertice, VerticeCorrectionMask>.MinimalDistanceSearchMode Correction, int FilterWidth = 10)
            : base (CellIndex, new int [1] { 1 })
        {
            Box = new BoundingBox (Grid.SpatialDimension);
            Grid.iGeomCells.GetCellBoundingBox (CellIndex, Box);
        }

        public override void AdjustParametersofCellPoints(MultidimensionalArray ResultValues, Dictionary<int, int> NbrofPointsPerType)
        {
            throw new NotImplementedException();
        }

        public override string PrintToCSV (out string header)
        {
            string lines = null;
            string header_part0 = "CellIndex,";
            string header_part1 = null;
            foreach (KeyValuePair<int, List<SingleVertice>> DictEntry in ItemList)
            {
                foreach (SingleVertice Point in DictEntry.Value)
                {
                    lines += CellIndex.ToString () + " , ";
                    lines += Point.PrintToCSV (out header_part1);
                    lines += "\n";
                }
            }
            header = header_part0 + header_part1 + "\n";

            return lines;
        }

        public override string PrintToString ()
        {
            string lines = null;
            foreach (KeyValuePair<int, List<SingleVertice>> DictEntry in ItemList)
            {
                foreach (SingleVertice Point in DictEntry.Value)
                {
                    lines += "CellIndex:" + string.Format ("{0,4}", CellIndex);
                    lines += Point.PrintToString ();
                    lines += "\n";
                }
            }
            return lines;
        }
    }

    public class EdgesofCell : SemiLagrangianLevelSetMethod
        <VerticesofCell,SingleVertice,VerticeCorrectionMask>.ItemsofCell
    {
        public new Dictionary<int, List<Edge>> ItemList;

        public EdgesofCell(int CellIndex, IGridData Grid, int FilterWidth = 10)
            : base(CellIndex, new int[1] { 1 })
        {
            ItemList = new Dictionary<int, List<Edge>>();
            ItemList[1] = new List<Edge>();
        }

        public override void AdjustParametersofCellPoints(MultidimensionalArray ResultValues, Dictionary<int, int> NbrofPointsPerType)
        {
            throw new NotImplementedException();
        }

        public override string PrintToCSV (out string header)
        {
            throw new NotImplementedException ();
        }

        public override string PrintToString ()
        {
            throw new NotImplementedException ();
        }
    }

    public class SingleVertice : SemiLagrangianLevelSetMethod
        <VerticesofCell, SingleVertice, VerticeCorrectionMask>.SinglePoint
    {
        public SingleVertice(MultidimensionalArray Coordinates, double AdvectionStepLengthINI)
            : base(Coordinates, AdvectionStepLengthINI)
        {
        }

        public new string PrintToString ()
        {
            string line = base.PrintToString ();
            line += "|Active:" + Active + " ";
            return line;
        }

        public new string PrintToCSV (out string header)
        {
            string line = base.PrintToCSV (out header);
            line += (Active ? 1 : 0);
            header += ",active";
            return line;
        }
    }

    public class Edge
    {
        public SingleVertice[] Vertice = new SingleVertice[2];
        public List<Edge> [] NeighbourEdge = new List<Edge> [2];
        public SingleVertice Center;
        public int CenterCellIndex;
        public MultidimensionalArray Normal;

        public Edge(SingleVertice Vertice_0,SingleVertice Vertice_1)
        {
            Vertice [0] = Vertice_0;
            NeighbourEdge [0] = new List<Edge> ();
            Vertice [1] = Vertice_1;
            NeighbourEdge [1] = new List<Edge> ();
            Center = new SingleVertice (GetCenter (),0);
            Normal = GetNormal();
        }

        public Edge(Edge OriginEdge,int VerticeIndex,SingleVertice NewVertice)
        {
            int NewEdgeVerticeConnectingIndex = (VerticeIndex == 0) ? 1 : 0;

            NeighbourEdge[NewEdgeVerticeConnectingIndex] = new List<Edge> { OriginEdge };
            OriginEdge.NeighbourEdge[VerticeIndex] = new List<Edge> { this };

            Vertice[NewEdgeVerticeConnectingIndex] = OriginEdge.Vertice[VerticeIndex];
            Vertice[VerticeIndex] = NewVertice;
            Center = new SingleVertice(GetCenter(),0);
            Normal = GetNormal();
        }

        public Edge(Edge Edge_1,int VerticeIndex_1,Edge Edge_2,int VerticeIndex_2)
        {
            Edge_1.NeighbourEdge[VerticeIndex_1] = new List<Edge> { this };
            Edge_2.NeighbourEdge[VerticeIndex_2] = new List<Edge> { this };
            NeighbourEdge[VerticeIndex_2] = new List<Edge> { Edge_1 };
            NeighbourEdge[VerticeIndex_1] = new List<Edge> { Edge_2 };
            Vertice[VerticeIndex_2] = Edge_1.Vertice[VerticeIndex_1];
            Vertice[VerticeIndex_1] = Edge_2.Vertice[VerticeIndex_2];
            Center = new SingleVertice(GetCenter(),0);
            Normal = GetNormal();
        }

        public MultidimensionalArray EdgeExtensionVector(int FrontIndex)
        {
            MultidimensionalArray OppositePoint = null; ;
            for(int i=0;i<Vertice.Length;i++)
            {
                if (i == FrontIndex) continue;
                OppositePoint = Vertice[i].Coordinates;
            }
            if (OppositePoint == null)
                throw new Exception("Failure Edge must not have more than two Vertices");
            MultidimensionalArray EdgeExtension;
            EdgeExtension = SemiLagrangianLevelSetMethod<VerticesofCell, SingleVertice, VerticeCorrectionMask>.
                ComputeVectorbetweenPoints(Vertice[FrontIndex].Coordinates, OppositePoint);
            return SemiLagrangianLevelSetMethod<VerticesofCell, SingleVertice, VerticeCorrectionMask>.NormalizeVector(EdgeExtension);
        }

        public MultidimensionalArray GetCenter ()
        {
            MultidimensionalArray Center = MultidimensionalArray.Create (1, Vertice [0].Coordinates.GetLength (1));
            for(int i=0;i<Vertice.Length;i++)
            {
                for(int dim=0;dim<Vertice[0].Coordinates.GetLength(1);dim++)
                {
                    Center [0, dim] += Vertice [i].Coordinates [0, dim] / Vertice.Length;
                }
            }
            return Center;
        }

        public double GetEdgeLength()
        {
            return SemiLagrangianLevelSetMethod
                <VerticesofCell, SingleVertice, VerticeCorrectionMask>.VectorNorm(ToVector());
        }

        public MultidimensionalArray GetNormal ()
        {
            MultidimensionalArray Vector = ToVector ();
            MultidimensionalArray Normal = MultidimensionalArray.Create (1, Vector.GetLength (1));
            Normal [0, 0] = -Vector [0, 1];
            Normal [0, 1] = Vector [0, 0];
            Normal = SemiLagrangianLevelSetMethod
            <VerticesofCell, SingleVertice, VerticeCorrectionMask>.NormalizeVector (Normal);
            if (SemiLagrangianLevelSetMethod<VerticesofCell, SingleVertice, VerticeCorrectionMask>.DotProduct (Vertice [0].Normal, Vertice [1].Normal) < 0)
                throw new Exception ("Failure Normal Vector of both Points in Edge must not be ");
            if (SemiLagrangianLevelSetMethod<VerticesofCell, SingleVertice, VerticeCorrectionMask>.DotProduct (Vertice [0].Normal, Normal) < 0)
                Normal = SemiLagrangianLevelSetMethod<VerticesofCell, SingleVertice, VerticeCorrectionMask>.InvertVector (Normal);

            return SemiLagrangianLevelSetMethod<VerticesofCell, SingleVertice, VerticeCorrectionMask>.NormalizeVector (Normal);
        }

        public void SetNormal()
        {
            MultidimensionalArray Normal = GetNormal ();
        
            if (SemiLagrangianLevelSetMethod<VerticesofCell, SingleVertice, VerticeCorrectionMask>.DotProduct (Normal, this.Normal) < 0)
                Console.WriteLine ("Something is wrong here: Normal Vector inverted");
            this.Normal = Normal;
        }

        public void SetCenter()
        {
            Center = new SingleVertice(GetCenter(),0);
        }

        public MultidimensionalArray ToVector()
        {
            MultidimensionalArray Vector = MultidimensionalArray.Create (1, Vertice [0].Coordinates.GetLength (1));
            Vector = SemiLagrangianLevelSetMethod<VerticesofCell, SingleVertice, VerticeCorrectionMask>
            .ComputeVectorbetweenPoints (Vertice [0].Coordinates, Vertice [1].Coordinates);
            return Vector;
        }

        public double DistanceToPoint(MultidimensionalArray Point)
        {
            double dist;
            double t;
            t = - (
                  SemiLagrangianLevelSetMethod<VerticesofCell, SingleVertice, VerticeCorrectionMask>.DotProduct (
                  SemiLagrangianLevelSetMethod<VerticesofCell, SingleVertice, VerticeCorrectionMask>.ComputeVectorbetweenPoints (Vertice [0].Coordinates, Point),
                  SemiLagrangianLevelSetMethod<VerticesofCell, SingleVertice, VerticeCorrectionMask>.ComputeVectorbetweenPoints (Vertice [1].Coordinates, Vertice [0].Coordinates))
                  )/(
                  SemiLagrangianLevelSetMethod<VerticesofCell, SingleVertice, VerticeCorrectionMask>.VectorNorm(
                  SemiLagrangianLevelSetMethod<VerticesofCell, SingleVertice, VerticeCorrectionMask>.ComputeVectorbetweenPoints (Vertice [1].Coordinates, Vertice [0].Coordinates))
                  ).Pow2 ();
            if(t > 1)
            {
                dist =
                SemiLagrangianLevelSetMethod<VerticesofCell, SingleVertice, VerticeCorrectionMask>.VectorNorm (
                SemiLagrangianLevelSetMethod<VerticesofCell, SingleVertice, VerticeCorrectionMask>.ComputeVectorbetweenPoints (Vertice [1].Coordinates, Point));
            }
            else if(t<0)
            {
                dist =
                SemiLagrangianLevelSetMethod<VerticesofCell, SingleVertice, VerticeCorrectionMask>.VectorNorm (
                SemiLagrangianLevelSetMethod<VerticesofCell, SingleVertice, VerticeCorrectionMask>.ComputeVectorbetweenPoints (Vertice [0].Coordinates, Point));
            }
            else
            {
                dist =
                SemiLagrangianLevelSetMethod<VerticesofCell, SingleVertice, VerticeCorrectionMask>.VectorNorm (
                SemiLagrangianLevelSetMethod<VerticesofCell, SingleVertice, VerticeCorrectionMask>.ComputeVectorbetweenPoints (
                SemiLagrangianLevelSetMethod<VerticesofCell, SingleVertice, VerticeCorrectionMask>.VectorAddition (t,
                SemiLagrangianLevelSetMethod<VerticesofCell, SingleVertice, VerticeCorrectionMask>.ComputeVectorbetweenPoints (
                Vertice [1].Coordinates, Vertice [0].Coordinates), 1, Vertice [0].Coordinates),
                Point));
            }
            return dist;
        }

        public double SignedDistanceToPoint(MultidimensionalArray Point)
        {
            int sign;
            MultidimensionalArray CenterToPoint = SemiLagrangianLevelSetMethod<VerticesofCell, SingleVertice, VerticeCorrectionMask>.ComputeVectorbetweenPoints (Point, Center.Coordinates);
            double DotProd = SemiLagrangianLevelSetMethod<VerticesofCell, SingleVertice, VerticeCorrectionMask>.DotProduct (Normal, CenterToPoint);
            sign = (DotProd > 0) ? 1 : -1;
            return sign * DistanceToPoint (Point);
        }
    }

    public class VerticeCorrectionMask : SemiLagrangianLevelSetMethod
        <VerticesofCell, SingleVertice, VerticeCorrectionMask>.CorrectionMask
    {
        public VerticeCorrectionMask (int CellIndex)
            : base (CellIndex, new int [1] { 1 })
        { }

        public VerticeCorrectionMask (int CellIndex, int [] NbrMovedPoints)
            : base (CellIndex, new int [1] { 1 })
        {
            if (NbrMovedPoints.Length != 1)
                throw new Exception ("VerticeCorrectionMask constructor needs exactly one initial values");
            MovedPoints [1] = NbrMovedPoints [0];
        }
    }

    public class FrontTrackingLevelSet2D : SemiLagrangianLevelSetMethod
        <VerticesofCell, SingleVertice, VerticeCorrectionMask>
    {
        protected  readonly double EdgeLengthToCurvatureFraction;
        protected double LongestEdge;

        protected readonly double MaxEdgeLengthToCurvatureFraction;
        protected readonly double MinEdgeLengthToCurvatureFraction;

        protected readonly double MaxEdgeLengthToCellFraction;
        protected readonly double MinEdgeLengthToCellFraction;

        protected Dictionary<int, int> CelltoArrayindexEdge;
        protected List<EdgesofCell> Edges;

        protected List<Edge> AllEdges;

        public struct FrontVertice
        {
            public int CellIndexGuess;
            public int VerticeIndex;
            public Edge Edge;
            public FrontVertice(Edge Edge, int VerticeIndex,int CellIndexGuess)
            {
                this.CellIndexGuess = CellIndexGuess;
                this.Edge = Edge;
                this.VerticeIndex = VerticeIndex;
            }
        }
        protected List<FrontVertice> ActiveFront;
        protected List<FrontVertice> PassiveFront;

        public FrontTrackingLevelSet2D (VectorField<SinglePhaseField> Velocity_New, VectorField<SinglePhaseField> Velocity_Old, SinglePhaseField Interface, LevelSetTracker InterfaceTrck,
            VectorField<SinglePhaseField> LevelSetGradient,int max_AddandAttract_Iteration, int NarrowBandWidth, IGridData Grid, bool LevelSetCorrectionWithNeighbourCells,
            int ReseedingInterval, MinimalDistanceSearchMode LevelSetCorrection, TopologyMergingMode TopologyMerging,
             NormalVectorDampingMode NormalVectorDamping, double EdgeLengthToCurvatureFraction)

            : base (Velocity_New, Velocity_Old, Interface, InterfaceTrck, LevelSetGradient, NarrowBandWidth, Grid, ReseedingInterval,
                 LevelSetCorrection, TopologyMerging, NormalVectorDamping,1)
        {
            this.EdgeLengthToCurvatureFraction = EdgeLengthToCurvatureFraction;
        }

        public override void Initialize ()
        {
            bool CellPosSure = false;
            sw.Start();
            TriangulateLevelSet ();
            sw.Stop();
            Console.WriteLine("Triangulate:                              {0}", sw.Elapsed);
            sw.Reset();
            sw.Start();
            AdjustVerticesParamters_Fast(CellPosSure);
            sw.Stop();
            Console.WriteLine("Adjusted Vertices Parameters:             {0}", sw.Elapsed);
            sw.Reset();
            sw.Start();
            ActualizeEdges();
            sw.Stop();
            Console.WriteLine("Actualize Edges:                          {0}", sw.Elapsed);
            sw.Reset();

        }

        public new void PerformTimestep (double dt, int substeps, int Timestep)
        {
            base.PerformTimestep(dt, substeps, Timestep);
            bool CellPosSure = false;
            sw.Start();
            Advect(dt, out List<VerticeCorrectionMask> CorrectCellPosMask, ref CellPosSure, substeps);
            sw.Stop();
            Console.WriteLine("Advected Vertices:                        {0}", sw.Elapsed);
            sw.Reset();
            sw.Start();
            ActualizeEdges();
            sw.Stop();
            Console.WriteLine("Actualize Edges:                          {0}", sw.Elapsed);
            sw.Reset();
            CorrectLevelSet(Timestep);
            sw.Start();
            AdjustVerticesParamters_Fast(CellPosSure);
            sw.Stop();
            Console.WriteLine("Adjusted Vertices Parameters:             {0}", sw.Elapsed);
            sw.Reset();
        }

        public void TriangulateLevelSet()
        {
            int nbrofVertices = 0;
            CelltoArrayindex = new Dictionary<int, int>();
            CelltoArrayindexEdge = new Dictionary<int, int>();

            Points = new List<VerticesofCell>();
            Edges = new List<EdgesofCell>();
            AllEdges = new List<Edge>();

            ActiveFront = new List<FrontVertice>();
            PassiveFront = new List<FrontVertice>();

            // Seed initial Edge
            CellMask CutCells = InterfaceTracker.Regions.GetCutCellMask ();
            int InitialCell = CutCells.ItemEnum.First();
            //Console.WriteLine ("InitialCell: " + InitialCell);

            MultidimensionalArray LocalPointIn = MultidimensionalArray.Create (1, SpatialDimension);
            MultidimensionalArray GlobalPointOut = MultidimensionalArray.Create (1, SpatialDimension);
            bool CellPosSure = true;

            MultidimensionalArray NewPoint;
            double CurvatureRadius;
            double CurvatureRadiusNewPoint;
            SingleVertice InitialVertice;
            MultidimensionalArray Tangential;
            do
            {
                do
                {
                    for (int d = 0; d < SpatialDimension; d++)
                    {
                        LocalPointIn [0, d] = (CoordGen.NextDouble () - 0.5) * 2;
                    }
                    Grid.TransformLocal2Global (LocalPointIn, GlobalPointOut, InitialCell);
                }
                while (!ProjectOnLevelSetValue (0.0,GlobalPointOut, InitialCell, ref CellPosSure));
                //Console.WriteLine("Initial Point: " + GlobalPointOut[0, 0] + " | " + GlobalPointOut[0, 1]);
                CellPosSure = false;
                InitialVertice = new SingleVertice (GlobalPointOut,0);
                //Console.WriteLine("Initial Point: "+InitialVertice.Coordinates[0,0]+" | "+InitialVertice.Coordinates[0,1]);
                InitialVertice.Normal = AdjustNormal (InitialVertice.Coordinates, InitialCell, CellPosSure);
                CurvatureRadius = CalculateRadius (InitialVertice.Coordinates, InitialCell);
                //Console.WriteLine("Curvature Radius: " + CurvatureRadius);
                Tangential = MultidimensionalArray.Create (1, SpatialDimension);
                Tangential [0, 0] = -InitialVertice.Normal [0, 1];
                Tangential [0, 1] = InitialVertice.Normal [0, 0];
                NewPoint = VectorAddition (1, InitialVertice.Coordinates, CurvatureRadius * EdgeLengthToCurvatureFraction, Tangential);
            }
            while (!ProjectOnLevelSetValue (0.0,NewPoint, InitialCell, ref CellPosSure));

            CurvatureRadiusNewPoint = CalculateRadius(InitialVertice.Coordinates, InitialCell);
            NewPoint = VectorAddition (1, InitialVertice.Coordinates, (CurvatureRadius + CurvatureRadiusNewPoint) / 2 * EdgeLengthToCurvatureFraction, Tangential);
            if (!ProjectOnLevelSetValue (0.0,NewPoint, InitialCell, ref CellPosSure))
                throw new Exception ("Failure in Initial Point projection");
            //Console.WriteLine("Second Point: " + NewPoint[0, 0] + " | " + NewPoint[0, 1]);

            SingleVertice SecondVertice = new SingleVertice (NewPoint,0);
            SecondVertice.Normal = AdjustNormal (SecondVertice.Coordinates, InitialCell, false);

            if (!AddVerticeToMemory (InitialVertice, InitialCell, MinimalDistanceSearch))
                throw new Exception ("Initial Vertice must be outside the Grid");
            nbrofVertices++;

            int CellIndex = FindCellIndexofPoint(SecondVertice.Coordinates);
            if (!AddVerticeToMemory (SecondVertice, CellIndex, MinimalDistanceSearch))
                throw new Exception ("Second Vertice must be outsinde the Grid");
            nbrofVertices++;

            Edge SeedEdge = new Edge (InitialVertice, SecondVertice);
            LongestEdge = SeedEdge.GetEdgeLength();

            CellIndex = FindCellIndexofPoint(SeedEdge.Center.Coordinates);
            //Console.WriteLine("Edge Center Index: " + CellIndex);
            if (!AddEdgeToMemory (SeedEdge, CellIndex))
                throw new Exception ("Seed Edge must be outside the Grid");

            ActiveFront.Add(new FrontVertice(SeedEdge, 0,InitialCell));
            ActiveFront.Add(new FrontVertice(SeedEdge, 1,InitialCell));

            // Continue Triangulation
            Edge NextEdge;
            int EdgeCellIndex;
            //Console.WriteLine("Start Triangulation:" + ActiveFront.Count);
            while (!ActiveFront.IsNullOrEmpty())
            {
                for(int k=0;k<ActiveFront.Count;k++)
                {
                    //Console.WriteLine();
                    MultidimensionalArray FrontVerticeCoord = ActiveFront[k].Edge.Vertice[ActiveFront[k].VerticeIndex].Coordinates;
                    Edge FrontEdge = ActiveFront[k].Edge;
                    //Console.WriteLine("Front Vertice:" + FrontVerticeCoord[0, 0]+ " | " + FrontVerticeCoord[0, 1]);
                    MultidimensionalArray NextPoint;
                    CellPosSure = false;
                    double multi = 1;
                    MultidimensionalArray ExtensionVector;
                    double CurvatureRadius_0 = CalculateRadius(ActiveFront[k].Edge.Vertice[ActiveFront[k].VerticeIndex].Coordinates, ActiveFront[k].CellIndexGuess);
                    do
                    {
                        ExtensionVector = FrontEdge.EdgeExtensionVector(ActiveFront[k].VerticeIndex);
                        //Console.WriteLine("Extension: " + ExtensionVector[0, 0] + " | " + ExtensionVector[0, 1]);
                        //Console.WriteLine("Origin Edge Length: " + FrontEdge.GetEdgeLength());
                        double ExtensionVecLen = FrontEdge.GetEdgeLength() * multi;
                        NextPoint = VectorAddition(1, FrontVerticeCoord,ExtensionVecLen,ExtensionVector);

                        //Console.WriteLine("NewPoint: " + NextPoint[0, 0] + " | " + NextPoint[0, 1]);
                        multi /= 2;
                    }
                    while (!ProjectOnLevelSetValue(0.0,NextPoint, ActiveFront[k].CellIndexGuess, ref CellPosSure));
                    double CurvatureRadius_1 = CalculateRadius(NextPoint, ActiveFront[k].CellIndexGuess);
                    NextPoint = VectorAddition(1, ActiveFront[k].Edge.Vertice[ActiveFront[k].VerticeIndex].Coordinates,
                            (CurvatureRadius_0+CurvatureRadius_1)*0.5*EdgeLengthToCurvatureFraction, ExtensionVector);
                    if (!ProjectOnLevelSetValue(0.0,NextPoint, ActiveFront[k].CellIndexGuess, ref CellPosSure))
                        throw new Exception("Failure in NextPoint Attraction while Triangulation");

                    //Console.WriteLine("Next Vertice:" + NextPoint[0, 0] + " | " + NextPoint[0, 1]);
                    bool EdgeDiscarded = false;
                    for(int i=0;i<ActiveFront.Count;i++)
                    {
                        if (i == k) continue;
                        double dist = VectorNorm(ComputeVectorbetweenPoints(NextPoint, ActiveFront[i].Edge.Vertice[ActiveFront[i].VerticeIndex].Coordinates));
                        if(dist < LongestEdge)
                        {
                            //Console.WriteLine(i);
                            double OriginEdgeLen = ActiveFront[i].Edge.GetEdgeLength();
                            double NeighbourEdgeLen = ActiveFront[k].Edge.GetEdgeLength();
                            if(dist<OriginEdgeLen*0.5 || dist<NeighbourEdgeLen*0.5)
                            {
                                PassiveFront.Add(ActiveFront[k]);
                                ActiveFront.RemoveAt(k--);
                                EdgeDiscarded = true;
                            }
                        }
                    }

                    for (int i = 0; i < PassiveFront.Count && EdgeDiscarded == false ; k++)
                    {
                        double dist = VectorNorm(ComputeVectorbetweenPoints(NextPoint, PassiveFront[i].Edge.Vertice[PassiveFront[i].VerticeIndex].Coordinates));
                        if (dist < LongestEdge)
                        {
                            double OriginEdgeLen = PassiveFront[i].Edge.GetEdgeLength();
                            double NeighbourEdgeLen = ActiveFront[k].Edge.GetEdgeLength();
                            if (dist < OriginEdgeLen*0.5 || dist < NeighbourEdgeLen*0.5)
                            {
                                PassiveFront.Add(ActiveFront[k]);
                                ActiveFront.RemoveAt(k--);
                                EdgeDiscarded = true;
                            }
                        }
                    }
                    //Console.WriteLine("New Edge accpeted: " + !EdgeDiscarded);

                    if(!EdgeDiscarded)
                    {
                        CellIndex = FindCellIndexofPoint(NextPoint, ActiveFront[k].CellIndexGuess);
                        SingleVertice NextVertice = new SingleVertice(NextPoint,0);
                        NextVertice.Normal = AdjustNormal(NextVertice.Coordinates, CellIndex, false);

                        if (!AddVerticeToMemory(NextVertice, CellIndex, MinimalDistanceSearch))
                            throw new Exception("Next Vertice must be outsinde the Grid");
                        nbrofVertices++;

                        NextEdge = new Edge(ActiveFront[k].Edge, ActiveFront[k].VerticeIndex, NextVertice);
                        double EdgeLength = NextEdge.GetEdgeLength();
                        LongestEdge = (LongestEdge < EdgeLength) ? EdgeLength : LongestEdge;

                        EdgeCellIndex = FindCellIndexofPoint(NextEdge.Center.Coordinates, CellIndex);
                        if (!AddEdgeToMemory(NextEdge, EdgeCellIndex))
                            throw new Exception("Next Edge must be outside the Grid");

                        FrontVertice NextFrontVertice = new FrontVertice(NextEdge, ActiveFront[k].VerticeIndex, CellIndex);
                        ActiveFront.RemoveAt(k--);
                        ActiveFront.Add(NextFrontVertice);
                    }
                }
            }

            //Close Remaining Gap

            if (PassiveFront.Count != 2)
                throw new Exception("More than two Ends to connect");
            NextEdge = new Edge(PassiveFront[0].Edge,PassiveFront[0].VerticeIndex,PassiveFront[1].Edge,PassiveFront[1].VerticeIndex);
            EdgeCellIndex = FindCellIndexofPoint(NextEdge.Center.Coordinates, PassiveFront[0].CellIndexGuess);
            //Console.WriteLine("Center of final Edge: " + NextEdge.Center.Coordinates[0, 0] + " | " + NextEdge.Center.Coordinates[0, 1]);
            if (!AddEdgeToMemory(NextEdge, EdgeCellIndex))
                throw new Exception("Final Edge must be outside the Grid");

            /*
            Console.WriteLine("Last Edge");
            Console.WriteLine(NextEdge.Vertice[0].Coordinates[0, 0] + " | " + NextEdge.Vertice[0].Coordinates[0, 1]);
            Console.WriteLine(NextEdge.Vertice[1].Coordinates[0, 0] + " | " + NextEdge.Vertice[1].Coordinates[0, 1]);

            Console.WriteLine("Nbr of Vertices added: " + nbrofVertices);
            */
        }

        public bool AddVerticeToMemory(SingleVertice Vertice,int CellIndexGuess, MinimalDistanceSearchMode LevelSetCorrection)
        {
            int CellIndex;
            if ((CellIndex = FindCellIndexofPoint (Vertice.Coordinates, CellIndexGuess)) == -1)
                return false;
            if (CelltoArrayindex.ContainsKey(CellIndex))
            {
                Points[CelltoArrayindex[CellIndex]].ItemList[1].Add(Vertice);
                //Console.Write("Added Vertice to Cell: " + CellIndex);
            }
            else
            {
                Points.Add(new VerticesofCell(CellIndex, Grid, LevelSetCorrection));
                CelltoArrayindex[CellIndex] = Points.Count - 1;
                Points[CelltoArrayindex[CellIndex]].ItemList[1].Add(Vertice);
            }
                return true;
        }

        public bool AddEdgeToMemory(Edge NewEdge,int CellIndexGuess)
        {
            NewEdge.Center.Coordinates = NewEdge.GetCenter ();
            int CellIndex;
            if ((CellIndex = FindCellIndexofPoint (NewEdge.Center.Coordinates, CellIndexGuess)) == -1)
                return false;

            AllEdges.Add (NewEdge);

            if (CelltoArrayindexEdge.ContainsKey (CellIndex))
                Edges [CelltoArrayindexEdge [CellIndex]].ItemList [1].Add (NewEdge);
            else
            {
                Edges.Add (new EdgesofCell (CellIndex, Grid));
                CelltoArrayindexEdge [CellIndex] = Edges.Count - 1;
                Edges[CelltoArrayindexEdge[CellIndex]].ItemList[1].Add(NewEdge);
            }
            return true;
        }

        protected void AdaptTriangulation ()
        {
            foreach (Edge Edge in AllEdges)
            {
                double CurvatureRadius = double.MaxValue;
                foreach (List<Edge> NeighbourEdgeList in Edge.NeighbourEdge)
                {
                    foreach (Edge Neighbour in NeighbourEdgeList)
                    {
                        double dist = VectorNorm (ComputeVectorbetweenPoints (Neighbour.Center.Coordinates, Edge.Center.Coordinates));
                        double normalAngle = VectorAngle (Edge.Normal, Neighbour.Normal);
                        double _CurvatureRadius = dist / (2 * Math.Sin (normalAngle / 2));
                        CurvatureRadius = (CurvatureRadius < _CurvatureRadius) ? CurvatureRadius : _CurvatureRadius;
                    }
                }
                double EdgeLength = Edge.GetEdgeLength ();
                double TargetEdgeLength = CurvatureRadius * EdgeLengthToCurvatureFraction;
                BoundingBox Box = new BoundingBox (SpatialDimension);
                Grid.iGeomCells.GetCellBoundingBox (Edge.CenterCellIndex, Box);
                double CellDiameter = Box.Diameter;
                if((EdgeLength/TargetEdgeLength) < MinEdgeLengthToCurvatureFraction)
                { }
            }
        }

        protected void ActualizeEdges()
        {
            foreach(Edge Edge in AllEdges)
            {
                Edge.SetNormal();
                Edge.GetCenter();
            }
        }

        public void AdjustVerticesParamters_Fast(bool CellPosSure)
        {
            DGField[] Fields = new DGField[1 + SpatialDimension];
            Fields[0] = Interface;
            for (int dim = 0; dim < SpatialDimension; dim++)
            {
                Fields[1 + dim] = InterfaceLevelSetGradient[dim];
            }
            MultidimensionalArray Results;

            foreach (VerticesofCell CellObj in Points)
            {
                int TotalNbrOfPoints = 0;
                Dictionary<int, int> NbrofPointsperType = new Dictionary<int, int>();
                foreach (KeyValuePair<int, List<SingleVertice>> DictEntry in CellObj.ItemList)
                {
                    TotalNbrOfPoints += DictEntry.Value.Count;
                    NbrofPointsperType.Add(DictEntry.Key, DictEntry.Value.Count);
                }
                MultidimensionalArray PointsCoordinates = MultidimensionalArray.Create(TotalNbrOfPoints, SpatialDimension);
                int indx = 0;
                foreach (KeyValuePair<int, List<SingleVertice>> DictEntry in CellObj.ItemList)
                {
                    for (int k = 0; k < NbrofPointsperType[DictEntry.Key]; k++)
                    {
                        for (int dim = 0; dim < SpatialDimension; dim++)
                        {
                            PointsCoordinates[indx, dim] = DictEntry.Value[k].Coordinates[0, dim];
                        }
                        indx++;
                    }
                }
                Results = MultiEvaluatePointsInCell(Fields, PointsCoordinates, CellObj.CellIndex, CellPosSure);
                indx = 0;
                foreach (KeyValuePair<int, List<SingleVertice>> DictEntry in CellObj.ItemList)
                {
                    foreach (SingleVertice SingleVertice in DictEntry.Value)
                    {
                        // Set Phi
                        if (SingleVertice.Phi == null) SingleVertice.Phi = MultidimensionalArray.Create(1, 1);
                        SingleVertice.Phi[0, 0] = Results[0, indx];

                        //Set Normal
                        MultidimensionalArray Normal_new = MultidimensionalArray.Create(1, SpatialDimension);
                        for (int dim = 0; dim < SpatialDimension; dim++) Normal_new[0, dim] = Results[dim + 1, indx];
                        Normal_new = NormalizeVector(Normal_new);
                        if (SingleVertice.Normal == null)
                            SingleVertice.Normal = Normal_new;
                        else
                        {
                            if (DotProduct(SingleVertice.Normal, Normal_new) < 0.0) SingleVertice.Active = false;
                            SingleVertice.Normal = Normal_new;
                        }
                        indx++;
                    }
                }
            }
        }

        protected override VerticesofCell NewPointsofCellObject (int CellIndex, IGridData Grid, MinimalDistanceSearchMode LevelSetCorrection)
        {
            return new VerticesofCell (CellIndex, Grid, LevelSetCorrection);
        }

        protected override VerticeCorrectionMask NewFullCorrectionMask (int CellIndex)
        {
            return new VerticeCorrectionMask(CellIndex, new int[1] { Points[CelltoArrayindex[CellIndex]].ItemList[1].Count });
        }

        protected override VerticeCorrectionMask NewCorrectionMask (int CellIndex)
        {
            return new VerticeCorrectionMask(CellIndex);
        }

        protected override void RemoveSurplusPoints (VerticesofCell TreatedCell, List<VerticeCorrectionMask> ParticleMask)
        {
            throw new NotImplementedException ();
        }

        protected override bool MovePointtoTarget (VerticesofCell TreatedCell, int signKey, MultidimensionalArray Coordinate, List<VerticeCorrectionMask> ParticleMask)
        {
            throw new NotImplementedException ();
        }

        protected override void CorrectLevelSet (int timestepNo, CellMask CellsToCorrect = null, bool ReIteration = false)
        {
            DGField Interface_old = Interface.CloneAs ();

            ScalarFunctionEx FrontTracking_LevelSetDistanceCorrectionFull = delegate (int cell0, int Len, NodeSet Ns, MultidimensionalArray result)
            {
                Interface_old.Evaluate (cell0, Len, Ns, result);
                MultidimensionalArray GlobalNodes = this.Grid.GlobalNodes.GetValue_Cell (Ns, cell0, Len);

                int c = 0;
                for (int cell = cell0; cell < cell0 + Len; cell++, c++)
                {
                    if (NarrowBandCells.Contains (cell))
                    //if(true)
                    {
                        for (int i = 0; i < Ns.NoOfNodes; i++)
                        {
                            MultidimensionalArray GlobalNode = MultidimensionalArray.Create (1, SpatialDimension);
                            for (int dim = 0; dim < SpatialDimension; dim++)
                                GlobalNode [0, dim] = GlobalNodes [c, i, dim];

                            double Phi_abs = AllEdges[0].SignedDistanceToPoint (GlobalNode);

                            foreach (Edge Edge in AllEdges)
                            {
                                double dist = Edge.SignedDistanceToPoint (GlobalNode);
                                if (Math.Abs(Phi_abs) > Math.Abs(dist))
                                {
                                    Phi_abs = dist;
                                }
                            }
                            result[c, i] = Phi_abs;
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
            };

            Interface.ProjectField(FrontTracking_LevelSetDistanceCorrectionFull);
        }
    }
}