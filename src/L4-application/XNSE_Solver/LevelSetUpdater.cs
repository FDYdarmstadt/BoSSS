using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.SpecFEM;
using BoSSS.Foundation.XDG;
using BoSSS.Solution;
using BoSSS.Solution.Control;
using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XdgTimestepping;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xml.Serialization;

namespace BoSSS.Application.XNSE_Solver
{
    struct DualLevelSet
    {
        public int LevelSetIndex;

        public LevelSet CGLevelSet;

        public LevelSet DGLevelSet;
    }

    class LevelSetUpdater
    {
        public LevelSetTracker Tracker;

        public XQuadFactoryHelper.MomentFittingVariants CutCellQuadratureType = XQuadFactoryHelper.MomentFittingVariants.Saye;
        
        Dictionary<string, ILevelSetEvolver> evolvers;

        DualLevelSet[] interfaces;

        ContinuityProjection[] enforcers;

        public LevelSetUpdater(GridData BackgroundGrid, XQuadFactoryHelper.MomentFittingVariants cutCellquadType, int __NearRegionWidth, string[] _SpeciesTable, LevelSet dgLevelSet)
        {
            LevelSet levelSet = ContinuityProjection.CreateField(
                    dgLevelSet, BackgroundGrid, ContinuityProjectionOption.SpecFEM);
            levelSet.AccLaidBack(1.0, dgLevelSet);
            interfaces = new DualLevelSet[]
            {
                new DualLevelSet
                {
                        LevelSetIndex = 0,
                        CGLevelSet = levelSet,
                        DGLevelSet = dgLevelSet,
                }
            };
            ContinuityProjection enforcer = new ContinuityProjection(
                levelSet.Basis,
                dgLevelSet.Basis,
                BackgroundGrid,
                ContinuityProjectionOption.SpecFEM);
            enforcers = new ContinuityProjection[]
            {
                enforcer
            };
            Tracker = new LevelSetTracker(BackgroundGrid, cutCellquadType, __NearRegionWidth, _SpeciesTable, levelSet);
            evolvers = new Dictionary<string, ILevelSetEvolver>(4);
        }

        public void AddEvolver(string name, ILevelSetEvolver evolver)
        {
            evolvers.Add(name, evolver);
        }

        public double UpdateLevelSets(
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields, 
            double time, 
            double dt, 
            double underRelax,
            bool incremental) 
        {
            UpdateInterfaces(interfaces, Tracker);
            double residual = 0;
            for(int i = 0; i < interfaces.Length;  ++i)
            {
                var singleInterface = interfaces[i];
                var enforcer = enforcers[i];
                if (evolvers.TryGetValue(singleInterface.CGLevelSet.Identification, out ILevelSetEvolver evolver))
                {
                    LevelSet ls = singleInterface.CGLevelSet;
                    LevelSet lsBkUp = ls.CloneAs();

                    //Move LevelSet and update Params
                    UpdateLevelSet(
                        evolver,
                        enforcer,
                        singleInterface,
                        DomainVarFields,
                        ParameterVarFields,
                        time,
                        dt,
                        underRelax,
                        incremental);

                    //Calculate Residual
                    CellMask oldCC = Tracker.Regions.GetCutCellMask4LevSet(singleInterface.LevelSetIndex);
                    var newCC = Tracker.Regions.GetCutCellMask();
                    lsBkUp.Acc(-1.0, ls);
                    double levSetResidual = lsBkUp.L2Norm(newCC.Union(oldCC));
                    residual += levSetResidual;
                }
            }
            Tracker.UpdateTracker(time + dt, incremental: true);
            return residual;
        }

        static void UpdateInterfaces(DualLevelSet[] interfaces, LevelSetTracker tracker)
        {
            for (int i = 0; i < interfaces.Length; ++i)
            {
                DualLevelSet combo = interfaces[i];
                combo.CGLevelSet = (LevelSet)tracker.LevelSets[i];
                combo.DGLevelSet.Clear();
                combo.DGLevelSet.AccLaidBack(1.0, combo.CGLevelSet);
            }
        }

        void UpdateLevelSet(
            ILevelSetEvolver evolver,
            ContinuityProjection enforcer,
            DualLevelSet phaseInterface,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields,
            double time,
            double dt,
            double underRelax,
            bool incremental)
        {
            LevelSet dglsBkUp = default;
            if (underRelax < 1.0)
            {
                
                dglsBkUp = phaseInterface.DGLevelSet.CloneAs();
            }
            
            evolver.UpdatePhaseInterface(
                phaseInterface,
                Tracker,
                time, 
                dt,
                incremental,
                DomainVarFields, 
                ParameterVarFields);
            
            //UnderRelax
            if (underRelax < 1.0)
            {
                LevelSet dgLs = phaseInterface.DGLevelSet;
                dgLs.Scale(underRelax);
                dgLs.Acc((1.0 - underRelax), dglsBkUp);
            }
            CellMask Near1 = Tracker.Regions.GetNearMask4LevSet(phaseInterface.LevelSetIndex, 1);
            CellMask PosFF = Tracker.Regions.GetLevelSetWing(phaseInterface.LevelSetIndex, +1).VolumeMask;
            enforcer.MakeContinuous(phaseInterface.DGLevelSet, phaseInterface.CGLevelSet, Near1, PosFF);
        }
    }
}
