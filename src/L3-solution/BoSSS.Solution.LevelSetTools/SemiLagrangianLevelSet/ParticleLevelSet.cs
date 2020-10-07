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
using BoSSS.Solution.LevelSetTools.Reinit.FastMarch;
using BoSSS.Solution.Statistic;
using ilPSP;

namespace BoSSS.Solution.LevelSetTool.SemiLagrangianLevelSet
{
    public class ParticlesofCell : SemiLagrangianLevelSetMethod
        <ParticlesofCell, SingleLvlSetParticle, ParticleCorrectionMask>.ItemsofCell
    {
        private readonly double Radius_min;
        private readonly double Radius_max;

        public ParticlesofCell(int CellIndex, IGridData Grid, ParticleLevelSet.MinimalDistanceSearchMode Correction, int FilterDepth = 3)
            : base(CellIndex, new int[2] { -1, 1 })
        {}

        public override void AdjustParametersofCellPoints(MultidimensionalArray ResultValues, Dictionary<int, int> NbrofPointsPerType)
        {
            int SpatialDimension = ResultValues.GetLength(0)-1;
            MultidimensionalArray LevelSet = MultidimensionalArray.Create(1, 1);
            MultidimensionalArray LevelSetGradient = MultidimensionalArray.Create(1, SpatialDimension);
            int index = 0;

            foreach(KeyValuePair<int,int> DictEntry in NbrofPointsPerType)
            {
                for(int k=0;k<DictEntry.Value;k++)
                {
                    LevelSet[0, 0] = ResultValues[0, index];
                    for(int dim=0;dim<SpatialDimension;dim++)
                    {
                        LevelSetGradient[0,dim] = ResultValues[dim + 1, index];
                    }
                    ItemList[DictEntry.Key][k].AdjustParameters(DictEntry.Key, Radius_max, Radius_min, LevelSet, LevelSetGradient);
                    index++;
                }
            }
        }

        public override string PrintToString()
        {
            string lines = null;
            foreach (KeyValuePair<int, List<SingleLvlSetParticle>> DictEntry in ItemList)
            {
                foreach (SingleLvlSetParticle Point in DictEntry.Value)
                {
                    lines += "CellIndex:"+string.Format("{0,4}", CellIndex);
                    lines += "|Sign:"+string.Format("{0,2}", DictEntry.Key)+"| ";
                    lines += Point.PrintToString();
                    lines += "\n";
                }
            }
            return lines;
        }

        public override string PrintToCSV(out string header)
        {
            string lines = null;
            string header_part0 = "CellIndex,sign,";
            string header_part1 = null;
            foreach (KeyValuePair<int, List<SingleLvlSetParticle>> DictEntry in ItemList)
            {
                foreach (SingleLvlSetParticle Point in DictEntry.Value)
                {
                    lines += CellIndex.ToString() + " , " + DictEntry.Key.ToString() + " , ";
                    lines += Point.PrintToCSV(out header_part1);
                    lines += "\n";
                }
            }
            header = header_part0 + header_part1 + "\n";

            return lines;
        }
    }

    public class SingleLvlSetParticle : SemiLagrangianLevelSetMethod
        <ParticlesofCell, SingleLvlSetParticle,ParticleCorrectionMask>.SinglePoint
    {
        public double Radius;
        public bool Escaped;

        public SingleLvlSetParticle(MultidimensionalArray Coordinates,double AdvectionStepLengthINI)
            : base(Coordinates,AdvectionStepLengthINI)
        {
            Escaped = false;
        }

        public void AdjustParameters(int sign,double Radius_max,double Radius_min,MultidimensionalArray LevelSet, MultidimensionalArray LevelSetGradient)
        {
            // Set Radius
            if (sign * LevelSet[0,0] < Radius_min) Radius = Radius_min;
            else if (sign * LevelSet[0,0] > Radius_max) Radius = Radius_max;
            else Radius = sign * LevelSet[0,0];

            MultidimensionalArray Normal_old = Normal.CloneAs();
            base.AdjustParameters(LevelSet, LevelSetGradient);
            int Direction = Math.Sign(SemiLagrangianLevelSetMethod
                <MarkersofCell, SingleLvlSetMarker, MarkerCorrectionMask>.
                DotProduct(Normal, Normal_old));
            if (Direction <= 0)
                Active = false;

            // Manage Escaped
            if (Escaped == false && -1 * sign * LevelSet[0, 0] >= Radius)
            {
                Escaped = true;
            }
            if (Escaped == true && -1 * sign * LevelSet[0, 0] < Radius)
            {
                Escaped = false;
            }
        }

        public new string PrintToString()
        {
            string line = base.PrintToString();
            line += "|Radius:"+string.Format("{0,7:0.00000}",Radius);
            line += "|Escaped:" + Escaped;
            line += "|Active:" + Active;
            return line;
        }

        public new string PrintToCSV(out string header)
        {
            string line = base.PrintToCSV(out header);
            line += Radius + " , ";
            line += (Escaped ? 1 : 0) + " , ";
            line += (Active ? 1 : 0);
                header += "radius,escaped,active";
            return line;
        }
    }

    public class ParticleCorrectionMask : SemiLagrangianLevelSetMethod
        <ParticlesofCell,SingleLvlSetParticle,ParticleCorrectionMask>.CorrectionMask
    { 
        public ParticleCorrectionMask(int CellIndex)
            : base(CellIndex,new int[2]{-1,1})
        { }

        public ParticleCorrectionMask(int CellIndex,int[] NbrMovedPoints)
            :base(CellIndex,new int[2]{-1,1})
        {
            if (NbrMovedPoints.Length != 2)
                throw new Exception("ParticleCorrectionMask constructor needs exactly two initial values");
            MovedPoints[-1] = NbrMovedPoints[0];
            MovedPoints[1] = NbrMovedPoints[1];
        }
    }

    public class Particle_K_D_TreeFilter : SemiLagrangianLevelSetMethod
        <ParticlesofCell, SingleLvlSetParticle, ParticleCorrectionMask>.K_D_TreeFilter
    {
        private readonly double Radius_max;
        private readonly double Radius_min;
        private readonly int sign;
         
        public Particle_K_D_TreeFilter(List<SingleLvlSetParticle> Points, int NbrOfPointsinNode,double Radius_max,double Radius_min,int sign)
            : base(Points, NbrOfPointsinNode)
        {
            this.Radius_max = Radius_max;
            this.Radius_min = Radius_min;
            this.sign = sign;
        }

        protected override void NearestPointAndMinRadiusFromList(List<SingleLvlSetParticle> Points, MultidimensionalArray GlobalNode, ref SingleLvlSetParticle NearestPoint, ref double MinRadius)
        {
            double Phi_correction = double.MinValue;
            foreach (SingleLvlSetParticle Particle in Points)
            {
                double dist = 0;
                for (int dim = 0; dim < SpatialDimension; dim++)
                {
                    dist += (GlobalNode[dim] - Particle.Coordinates[0, dim]).Pow2();
                }
                dist = Math.Sqrt(dist);

                double Phi_p = sign * (Particle.Radius - dist);
                Phi_correction = (sign * Phi_correction > sign * Phi_p) ? Phi_correction : Phi_p;
            }
        }
    }

    public class ParticleLevelSet : SemiLagrangianLevelSetMethod
        <ParticlesofCell,SingleLvlSetParticle,ParticleCorrectionMask>
    {
        public readonly double Band_max;
        public readonly double Band_min;
        public readonly double Radius_max;
        public readonly double Radius_min;
        public readonly int max_AddandAttract_Iteration;

        private readonly int MaxSearchDepth = 3;
        public readonly int FilterDepth;
        private CellMask CorrectionCells;
        private FastMarchReinit Reinitialization;

        public ParticleLevelSet(VectorField<SinglePhaseField> Velocity_New, VectorField<SinglePhaseField> Velocity_Old,
            SinglePhaseField Interface, LevelSetTracker InterfaceTrck, VectorField<SinglePhaseField> LevelSetGradient,
            int TargetNbrParticlesPerCell, int UpperLimitParticlesPerCell, int LowerLimitParticlesPerCell, double Band_max, double Band_min,
            double Radius_max, double Radius_min, int max_AddandAttract_Iteration, int NarrowBandWidth, IGridData Grid,
            bool LevelSetCorrectionWithNeighbourCells, int ReseedingInterval,MinimalDistanceSearchMode LevelSetCorrection,
            TopologyMergingMode TopologyMerging,NormalVectorDampingMode NormalVectorDamping)

            : base(Velocity_New,Velocity_Old,Interface,InterfaceTrck,LevelSetGradient,NarrowBandWidth,Grid,ReseedingInterval,LevelSetCorrection,
                  TopologyMerging, NormalVectorDamping,TargetNbrParticlesPerCell,UpperLimitParticlesPerCell,LowerLimitParticlesPerCell)
        {
            if (TargetNbrParticlesPerCell <= 0) throw new Exception ("Target number of particles per cell must not be lower than 1.");
            if (TargetNbrParticlesPerCell > UpperLimitParticlesPerCell) throw new Exception ("The upper limit of particles per cell must be equal or higher than the target number of particles per cell.");
            if (TargetNbrParticlesPerCell < LowerLimitParticlesPerCell) throw new Exception ("The lower limit of particles per cell must be equal or lower than the target number of particles per cell.");
            this.Band_max = Band_max;
            this.Band_min = Band_min;
            this.Radius_max = Radius_max;
            this.Radius_min = Radius_min;
            this.max_AddandAttract_Iteration = max_AddandAttract_Iteration;
            CorrectionCells = InterfaceTracker.Regions.GetNearFieldMask(1);
            ReseedingWdith = NarrowBandWidth;
            Reinitialization = new FastMarchReinit(Interface.Basis);
        }

        public override void Initialize()
        {
            bool CellPosSure = true;
            sw.Start();
            InitializePointstoLevelSet(out List<ParticleCorrectionMask> CorrectCellPosMask, ref CellPosSure);
            sw.Stop();
            Console.WriteLine("InitializeParticles and Reseed Elapsed:  {0}", sw.Elapsed);
            sw.Reset();
            sw.Start();
            CorrectCellLocationOfPoints(CorrectCellPosMask, ref CellPosSure);
            sw.Stop();
            Console.WriteLine("Correct Cell Location Elapsed:           {0}", sw.Elapsed);
            sw.Reset();
            sw.Start();
            AdjustParticleParamters(CellPosSure);
            sw.Stop();
            Console.WriteLine("Adjust Parameters Elapsed:               {0}", sw.Elapsed);
            sw.Reset();
            PointCount();
            RemoveNonActivePoints(out _);
        }

        public new void PerformTimestep(double dt, int substeps, int Timestep)
        {
            base.PerformTimestep(dt, substeps, Timestep);
            bool CellPosSure = false;
            sw.Start();
            Advect(dt, out List<ParticleCorrectionMask> CorrectCellPosMask, ref CellPosSure, substeps);
            sw.Stop();
            Console.WriteLine("Advected Particles:                      {0}", sw.Elapsed);
            sw.Reset();
            sw.Start();
            CorrectCellLocationOfPoints(CorrectCellPosMask, ref CellPosSure);
            sw.Stop();
            Console.WriteLine("Correct Cell Location Elapsed:           {0}", sw.Elapsed);
            sw.Reset();
            sw.Start();
            CheckForEscapedParticles(MinimalDistanceSearch, CellPosSure);
            sw.Stop();
            Console.WriteLine("Check for Escaped Particles Elapsed:     {0}", sw.Elapsed);
            sw.Reset();
            sw.Start();
            CorrectLevelSet(Timestep);
            sw.Stop();
            Console.WriteLine("Correct Level Set Elapsed:               {0}", sw.Elapsed);
            sw.Reset();
            sw.Start();
            Reinitialization.FirstOrderReinit(InterfaceTracker, "A", NarrowBandWidth);
            sw.Stop();
            Console.WriteLine("Reinitialize Narrrow Band:               {0}", sw.Elapsed);
            sw.Reset();
            sw.Start();
            CorrectLevelSet(Timestep);
            sw.Stop();
            Console.WriteLine("Correct Level Set Elapsed:               {0}", sw.Elapsed);
            sw.Reset();
            sw.Start();
            AdjustPointParamters(Timestep,CellPosSure);
            sw.Stop();
            Console.WriteLine("Adjust Parameters Elapsed:               {0}", sw.Elapsed);
            sw.Reset();
            PointCount();
            RemoveNonActivePoints(out _);

            if ((Timestep % ReseedingInterval) == 0)
            {
                Console.WriteLine("--------------------Reseeding--------------------");
                sw.Start();
                ReseedingPoints(out List<ParticleCorrectionMask> ReseedingParticleMask);
                sw.Stop();
                Console.WriteLine("Reseed Particles Elapsed:                {0}", sw.Elapsed);
                sw.Reset();
                sw.Start();
                CorrectCellLocationOfPoints(ReseedingParticleMask, ref CellPosSure);
                sw.Stop();
                Console.WriteLine("Correct Cell Location Elapsed:           {0}", sw.Elapsed);
                sw.Reset();
                sw.Start();
                AdjustPointParamters(Timestep,CellPosSure);
                sw.Stop();
                Console.WriteLine("Adjust Parameters Elapsed:               {0}", sw.Elapsed);
                sw.Reset();
                sw.Start();
                RemoveFarParticles(ref CellPosSure);
                sw.Stop();
                Console.WriteLine("Remove Far Elapsed:                      {0}", sw.Elapsed);
                sw.Reset();
                Console.WriteLine("------------------------------------------------");
            }
        }

        protected override bool MovePointtoTarget(ParticlesofCell TreatedCell, int signKey, MultidimensionalArray Coordinate, List<ParticleCorrectionMask> ParticleMask)
        {
            BoundingBox Box = new BoundingBox(SpatialDimension);
            Grid.iGeomCells.GetCellBoundingBox(TreatedCell.CellIndex, Box);
            double AdvectionStepLengthINI = Box.Diameter * 1e-2;

            MultidimensionalArray Phi;
            MultidimensionalArray Normal;
            double lambda = 1;
            bool cellPosSure = true;
            for (int it = 0; it < max_AddandAttract_Iteration - 1; it++, cellPosSure = false)
            {
                Phi = Evaluate(new DGField[]{Interface}, Coordinate,TreatedCell.CellIndex,cellPosSure);
                //Console.WriteLine("Initial Phi: " + Phi[0, 0]);
                if (signKey * Phi[0, 0] < Band_max && signKey * Phi[0, 0] > Band_min)
                {
                    TreatedCell.ItemList[signKey].Add(new SingleLvlSetParticle(Coordinate,AdvectionStepLengthINI));

                    ParticleCorrectionMask other = new ParticleCorrectionMask(TreatedCell.CellIndex, new int[] { -1, 1 });
                    int index;
                    if ((index = ParticleMask.FindIndex(other.Compare)) == -1)
                    {
                        ParticleMask.Add(new ParticleCorrectionMask(TreatedCell.CellIndex));
                        ParticleMask.Last().MovedPoints[signKey]++;
                    }
                    else ParticleMask[index].MovedPoints[signKey]++;
                    return true;
                }
                double Phi_goal = CoordGen.NextDouble() * signKey * (Band_max - Band_min) + signKey * Band_min;
                Normal = EvaluateNormal(Coordinate, TreatedCell.CellIndex, cellPosSure);
                for (int d = 0; d < SpatialDimension; d++)
                {
                    Coordinate[0, d] += lambda * (Phi_goal - Phi[0, 0]) * Normal[0, d];
                }
                cellPosSure = false;
                //PrintParticleSpecs(GlobalPointOut, Normal, TreatedCell.CellIndex, "It_2 " + it + " Plus");
                lambda /= 1 + 1 / max_AddandAttract_Iteration;
            }
            return false;
        }

        /// <summary>
        /// Removes all surplus particles of a cell.
        /// Sorts all particles for both plus and minus using the ParticleComparer.
        /// Removes any non-escaped particle starting with the end of the List.
        /// The Phi values at any point have to be calculated at the current coordinates
        /// </summary>
        protected override void RemoveSurplusPoints(ParticlesofCell TreatedCell,List<ParticleCorrectionMask> ParticleMask)
        {
            Dictionary<int, int> PointstoRemove = new Dictionary<int, int>();
            int ParticleCount = 0;
            foreach(KeyValuePair<int, List<SingleLvlSetParticle>> DictEntry in TreatedCell.ItemList)
            {
                PointstoRemove.Add(DictEntry.Key,0);
                ParticleCount += DictEntry.Value.Count;
            }
            int Faktor = 0;
            Faktor = TreatedCell.ItemList[1].Count() / (ParticleCount);

            PointstoRemove[1] = (1 - Faktor) * TreatedCell.ItemList[1].Count();
            PointstoRemove[-1] = ParticleCount - PointstoRemove[1] - TargetNbrofPointsPerCellPerDimension*SpatialDimension;

            ParticleComparer Pc = new ParticleComparer();
            foreach(KeyValuePair<int,int>DictEntry in PointstoRemove)
            {
                if (DictEntry.Value != 0)
                {
                    TreatedCell.ItemList[DictEntry.Key].Sort(Pc);
                    int indx = TreatedCell.ItemList[DictEntry.Key].Count - 1;
                    for (int i = 0; i < DictEntry.Value; i++)
                    {
                        if (TreatedCell.ItemList[DictEntry.Key][indx].Escaped == false) TreatedCell.ItemList[DictEntry.Key].RemoveAt(indx);
                        indx--;
                    }
                }
            }
        }

        /// <summary>
        /// Compares Particles
        /// Uses the function: Sign(p) Phi(x) - Radius
        /// If the function of a particle is greater than that of another the particle is greater
        /// </summary>
        public class ParticleComparer : IComparer<SingleLvlSetParticle>
        {
            public int Compare(SingleLvlSetParticle A, SingleLvlSetParticle B)
            {
                double a, b;
                if (A.Phi[0, 0] > 0 || B.Phi[0, 0] > 0)
                {
                    a = A.Phi[0, 0] - A.Radius;
                    b = B.Phi[0, 0] - B.Radius;
                }
                else if (A.Phi[0, 0] < 0 || B.Phi[0, 0] < 0)
                {
                    a = -A.Phi[0, 0] - A.Radius;
                    b = -B.Phi[0, 0] - B.Radius;
                }
                else
                    return 0;

                if (a < b)
                    return -1;
                else if (b < a)
                    return 1;
                else
                    return 0;
            }
        }

        /// <summary>
        /// Removes Particles that are further away from the Interface than Band max and not escaped.
        /// If both particle lists are empty after the process the particle cell is deleted.
        /// Method must not be used if the cell index of any particle is doubtful!
        /// </summary>
        private void RemoveFarParticles(ref bool CellPosSure)
        {
            foreach(ParticlesofCell CellObj in Points)
            {
                foreach(KeyValuePair<int,List<SingleLvlSetParticle>> PointList in CellObj.ItemList)
                {
                    for(int nbr=0;nbr<PointList.Value.Count;nbr++)
                    {
                        if ((PointList.Value[nbr].Phi[0, 0] > PointList.Key*Band_max && !PointList.Value[nbr].Escaped) || (!PointList.Value[nbr].Active && PointList.Value[nbr].Escaped))
                        {
                            PointList.Value.RemoveAt(nbr--);
                        }
                    }
                }
            }
        }

        /// <summary>
        /// Adjusts the radius of every particle in the Particle class.
        /// This method must not be used if the cell index of any particle is doubtful.
        /// </summary>
        private void CorrectRadiusofAllParticles(bool CellPosSure = false)
        {
            foreach (ParticlesofCell CellObj in Points)
            {
                foreach (KeyValuePair<int, List<SingleLvlSetParticle>> DictEntry in CellObj.ItemList)
                {
                    foreach (SingleLvlSetParticle SingleParticle in DictEntry.Value)
                    {
                        // Berechnung der lokalen Level Set Werte eines jeden Partikels in einer Zelle
                        SingleParticle.Phi = Evaluate(new DGField[] { Interface }, SingleParticle.Coordinates, CellObj.CellIndex, CellPosSure);
                        if (DictEntry.Key * SingleParticle.Phi[0, 0] < Radius_min) SingleParticle.Radius = Radius_min;
                        else if (DictEntry.Key * SingleParticle.Phi[0, 0] > Radius_max) SingleParticle.Radius = Radius_max;
                        else SingleParticle.Radius = DictEntry.Key * SingleParticle.Phi[0, 0];
                        //Console.WriteLine("CellIndex: " + CellObj.CellIndex + " Plus Radius: " + SingleParticle.Radius);
                    }
                }
            }
        }

        public void CheckForEscapedParticles(MinimalDistanceSearchMode Correction, bool CellPosSure = false)
        {
            foreach (ParticlesofCell CellObj in Points)
            {
                foreach (KeyValuePair<int, List<SingleLvlSetParticle>> DictEntry in CellObj.ItemList)
                {
                    foreach (SingleLvlSetParticle SingleParticle in DictEntry.Value)
                    {
                        SingleParticle.Phi = Evaluate(new DGField[] { Interface }, SingleParticle.Coordinates, CellObj.CellIndex, CellPosSure);
                        MultidimensionalArray Normal_new = AdjustNormal(SingleParticle.Coordinates, CellObj.CellIndex, CellPosSure);
                        if (DotProduct(SingleParticle.Normal, Normal_new) < 0.0) SingleParticle.Active = false;
                        SingleParticle.Normal = Normal_new;

                        if (SingleParticle.Escaped == false && -1 * DictEntry.Key * SingleParticle.Phi[0, 0] >= SingleParticle.Radius)
                        {
                            SingleParticle.Escaped = true;
                        }
                        if (SingleParticle.Escaped == true && -1 * DictEntry.Key * SingleParticle.Phi[0, 0] < SingleParticle.Radius)
                        {
                            SingleParticle.Escaped = false;
                        }

                    }
                }
            }
        }

        public void CheckForEscapedParticles_Fast(MinimalDistanceSearchMode Correction, bool CellPosSure = false)
        {
            DGField[] Fields = new DGField[1+SpatialDimension];
            Fields[0] = Interface;
            for(int dim=0;dim<SpatialDimension;dim++)
            {
                Fields[1 + dim] = InterfaceLevelSetGradient[dim];
            }
            MultidimensionalArray Results;

            foreach (ParticlesofCell CellObj in Points)
            {
                int TotalNbrOfPoints = 0;
                Dictionary<int, int> NbrofPointsperType = new Dictionary<int, int>();
                foreach (KeyValuePair<int, List<SingleLvlSetParticle>> DictEntry in CellObj.ItemList)
                {
                    TotalNbrOfPoints += DictEntry.Value.Count;
                    NbrofPointsperType.Add(DictEntry.Key, DictEntry.Value.Count);
                }
                MultidimensionalArray PointsCoordinates = MultidimensionalArray.Create(TotalNbrOfPoints, SpatialDimension);
                int indx = 0;
                foreach (KeyValuePair<int, List<SingleLvlSetParticle>> DictEntry in CellObj.ItemList)
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
                foreach (KeyValuePair<int,List<SingleLvlSetParticle>> DictEntry in CellObj.ItemList)
                {
                    foreach (SingleLvlSetParticle SingleParticle in DictEntry.Value)
                    {
                        SingleParticle.Phi[0, 0] = Results[0, indx];
                        MultidimensionalArray Normal_new = MultidimensionalArray.Create(1, SpatialDimension);
                        for (int dim = 0; dim < SpatialDimension; dim++) Normal_new[0, dim] = Results[dim + 1, indx];
                        Normal_new = NormalizeVector(Normal_new);
                        if (DotProduct(SingleParticle.Normal, Normal_new) < 0.0) SingleParticle.Active = false;
                        SingleParticle.Normal = Normal_new;

                        if (SingleParticle.Escaped == false && -1 * DictEntry.Key * SingleParticle.Phi[0, 0] >= SingleParticle.Radius)
                        {
                            SingleParticle.Escaped = true;
                        }
                        if (SingleParticle.Escaped == true && -1 * DictEntry.Key * SingleParticle.Phi[0, 0] < SingleParticle.Radius)
                        {
                            SingleParticle.Escaped = false;
                        }
                        indx++;
                    }
                }
            }
        }


        /// <summary>
        /// Correction of LevelSet with escaped particles.
        /// First all particle are assigned with a recalculated level set phi
        /// Then they are marked as escaped depending on their position in the level set DGField.
        /// Finally a spherical local level set function is constructed around every escaped particle.
        /// These functions are overlapped using the following rules:
        /// For the plus(minus) particles the maximum(minimum) of the function is used.
        /// Positiv and negative particles are overlapped by prioritizing the level set corrections that are nearer to the interface.
        /// Method has to be called after the advection step but before the radius adjustment
        /// </summary>
        protected override void CorrectLevelSet(int timestepNo, CellMask CellsToCorrect = null, bool ReIteration = false)
        {
            DGField Interface_old = Interface.CloneAs();

            ScalarFunctionEx Particle_LevelSetCorrection = delegate (int cell0, int Len, NodeSet Ns, MultidimensionalArray result)
            {
                Interface_old.Evaluate(cell0, Len, Ns, result);

                MultidimensionalArray GlobalNodes = this.Grid.GlobalNodes.GetValue_Cell(Ns, cell0, Len);

                int c = 0;
                for (int cell = cell0; cell < cell0 + Len; cell++, c++)
                {
                    if (NarrowBandCells.Contains(cell))
                    {
                        List<int> listindex = new List<int>();
                        if (CelltoArrayindex.ContainsKey(cell)) listindex.Add(cell);

                        Grid.GetCellNeighbours(cell, GetCellNeighbours_Mode.ViaVertices, out int[] NeighbourIndex, out int[] ConnectingEntities);
                        foreach (int j in NeighbourIndex)
                        {
                            if (CelltoArrayindex.ContainsKey(j))
                            {
                                listindex.Add(j);
                            }
                        }
                        if (listindex.IsNullOrEmpty()) continue;

                        for (int l = 0; l < listindex.Count; l++)
                        {
                            listindex[l] = CelltoArrayindex[listindex[l]];
                        }


                        for (int i = 0; i < Ns.NoOfNodes; i++)
                        {
                            Dictionary<int, double> Phi_correction = new Dictionary<int, double>
                            {
                                [-1] = result[c, i],
                                [1] = result[c, i]
                            };
                            if (MinimalDistanceSearch == MinimalDistanceSearchMode.FullSearch)
                            {
                                foreach (int index in listindex)
                                {
                                    foreach (KeyValuePair<int, List<SingleLvlSetParticle>> DictEntry in Points[index].ItemList)
                                    {
                                        foreach (SingleLvlSetParticle SingleParticle in DictEntry.Value)
                                        {
                                            if (SingleParticle.Escaped == false || SingleParticle.Active == false) continue;
                                            double dist = 0;
                                            for (int dim = 0; dim < SpatialDimension; dim++)
                                            {
                                                dist += (GlobalNodes[c, i, dim] - SingleParticle.Coordinates[0, dim]).Pow2();
                                            }
                                            dist = Math.Sqrt(dist);
                                            double Phi_p = DictEntry.Key * (SingleParticle.Radius - dist);
                                            Phi_correction[DictEntry.Key] = (DictEntry.Key * Phi_correction[DictEntry.Key] > DictEntry.Key * Phi_p)
                                            ? Phi_correction[DictEntry.Key] : Phi_p;
                                        }
                                    }
                                }
                            }
                            if (Math.Abs(Phi_correction[1]) <= Math.Abs(Phi_correction[-1]))
                                result[c, i] = Phi_correction[1];
                            else
                                result[c, i] = Phi_correction[-1];
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

            Interface.ProjectField(Particle_LevelSetCorrection);
        }

        protected override ParticlesofCell NewPointsofCellObject(int CellIndex, IGridData Grid, MinimalDistanceSearchMode LevelSetCorrection)
        {
            return new ParticlesofCell(CellIndex, Grid, LevelSetCorrection, FilterDepth);
        }

        protected override ParticleCorrectionMask NewFullCorrectionMask(int CellIndex)
        {
            int[] NbrofPointsMoved = new int[2];
            int indexinPointList = CelltoArrayindex[CellIndex];
            NbrofPointsMoved[0] = Points[indexinPointList].ItemList[-1].Count;
            NbrofPointsMoved[1] = Points[indexinPointList].ItemList[1].Count;
            return new ParticleCorrectionMask(CellIndex, NbrofPointsMoved);
        }

        protected override ParticleCorrectionMask NewCorrectionMask(int CellIndex)
        {
            int[] NbrofPointsMoved = new int[2];
            int indexinPointList = CelltoArrayindex[CellIndex];
            NbrofPointsMoved[0] = 0;
            NbrofPointsMoved[1] = 0;
            return new ParticleCorrectionMask(CellIndex, NbrofPointsMoved);
        }

        public void AdjustParticleParamters(bool CellPosSure)
        {
            foreach (ParticlesofCell CellObj in Points)
            {
                foreach (KeyValuePair<int,List<SingleLvlSetParticle>> DictEntry in CellObj.ItemList)
                {
                    foreach (SingleLvlSetParticle SingleParticle in DictEntry.Value)
                    {
                        // Evaluate particle properties on field
                        MultidimensionalArray FieldsResults = MultidimensionalArray.Create(1, 1 + SpatialDimension);
                        DGField[] FieldstoEval = new DGField[1 + SpatialDimension];
                        FieldstoEval[0] = Interface;
                        for (int dim = 0; dim < SpatialDimension; dim++) FieldstoEval[1 + dim] = InterfaceLevelSetGradient[dim];
                        FieldsResults = Evaluate(FieldstoEval, SingleParticle.Coordinates, CellObj.CellIndex, CellPosSure);

                        // Set Phi
                        if (SingleParticle.Phi == null) SingleParticle.Phi = MultidimensionalArray.Create(1, 1);
                        SingleParticle.Phi[0, 0] = FieldsResults[0, 0];

                        // Set Radius
                        if (DictEntry.Key * SingleParticle.Phi[0, 0] < Radius_min) SingleParticle.Radius = Radius_min;
                        else if (DictEntry.Key * SingleParticle.Phi[0, 0] > Radius_max) SingleParticle.Radius = Radius_max;
                        else SingleParticle.Radius = DictEntry.Key * SingleParticle.Phi[0, 0];

                        // Manage Escaped
                        if (SingleParticle.Escaped == false && -1 * DictEntry.Key * SingleParticle.Phi[0, 0] >= SingleParticle.Radius)
                        {
                            SingleParticle.Escaped = true;
                            //PrintParticleSpecs(SingleParticle.Coordinates, SingleParticle.Normal, CellObj.CellIndex, "Escaped " + SIGN);
                        }
                        if (SingleParticle.Escaped == true && -1 * DictEntry.Key * SingleParticle.Phi[0, 0] < SingleParticle.Radius)
                        {
                            SingleParticle.Escaped = false;
                            //PrintParticleSpecs(SingleParticle.Coordinates, SingleParticle.Normal, CellObj.CellIndex, "NonEscaped " + SIGN);
                        }

                        //Set Normal
                        MultidimensionalArray Normal = MultidimensionalArray.Create(1, SpatialDimension);
                        // Calculation of eulerian gradient norm
                        for (int dim = 0; dim < SpatialDimension; dim++)
                        {
                            Normal[0, dim] = FieldsResults[0, dim+1];
                        }
                        SingleParticle.Normal = NormalizeVector(Normal); 
                    }
                }
            }
        }
    }
}
