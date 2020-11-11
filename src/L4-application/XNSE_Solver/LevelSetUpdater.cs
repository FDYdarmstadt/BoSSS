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
        
        Dictionary<string, ILevelSetHandler> lsHandlers;

        DualLevelSet[] currentInterfaces;

        ContinuityProjection[] enforcers;

        public LevelSetUpdater(GridData BackgroundGrid, XQuadFactoryHelper.MomentFittingVariants cutCellquadType, int __NearRegionWidth, string[] _SpeciesTable, LevelSet dgLevelSet)
        {
            LevelSet levelSet = ContinuityProjection.CreateField(
                    dgLevelSet, BackgroundGrid, ContinuityProjectionOption.SpecFEM);
            levelSet.AccLaidBack(1.0, dgLevelSet);
            currentInterfaces = new DualLevelSet[]
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
            lsHandlers = new Dictionary<string, ILevelSetHandler>(4);
        }

        public void AddEvolver(string name, ILevelSetHandler evolver)
        {
            lsHandlers.Add(name, evolver);
        }

        public double UpdateLevelSets(
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields, 
            double time, 
            double dt, 
            double underRelax,
            bool incremental) 
        {
            SetCurrentInterfaces(currentInterfaces, Tracker);
            double residual = 0;
            for(int i = 0; i < currentInterfaces.Length;  ++i)
            {
                var singleInterface = currentInterfaces[i];
                var enforcer = enforcers[i];
                if (lsHandlers.TryGetValue(singleInterface.CGLevelSet.Identification, out ILevelSetHandler handler))
                {
                    LevelSet ls = singleInterface.CGLevelSet;
                    LevelSet lsBkUp = ls.CloneAs();

                    //Move LevelSet and update Params
                    MoveLevelSet(
                        handler,
                        enforcer,
                        singleInterface,
                        DomainVarFields,
                        ParameterVarFields,
                        time,
                        dt,
                        underRelax,
                        incremental);

                    handler.UpdateParameters(
                        singleInterface,
                        Tracker,
                        time,
                        ParameterVarFields);

                    //Calculate Residual
                    CellMask oldCC = Tracker.Regions.GetCutCellMask4LevSet(singleInterface.LevelSetIndex);
                    var newCC = Tracker.Regions.GetCutCellMask();
                    lsBkUp.Acc(-1.0, ls);
                    double levSetResidual = lsBkUp.L2Norm(newCC.Union(oldCC));
                    residual += levSetResidual;
                }
                else
                {
                    throw new Exception($"LevelSet #{i} does not have a registered handler");
                }
            }
            Tracker.UpdateTracker(time + dt, incremental: true);
            return residual;
        }

        static void SetCurrentInterfaces(DualLevelSet[] interfaces, LevelSetTracker tracker)
        {
            for (int i = 0; i < interfaces.Length; ++i)
            {
                DualLevelSet combo = interfaces[i];
                combo.CGLevelSet = (LevelSet)tracker.LevelSets[i];
                combo.DGLevelSet.Clear();
                combo.DGLevelSet.AccLaidBack(1.0, combo.CGLevelSet);
            }
        }

        void MoveLevelSet(
            ILevelSetHandler lsHandler,
            ContinuityProjection enforcer,
            DualLevelSet phaseInterface,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields,
            double time,
            double dt,
            double underRelax,
            bool incremental)
        {
            LevelSet dglsBkUp = null;
            if (underRelax < 1.0)
            {
                dglsBkUp = phaseInterface.DGLevelSet.CloneAs();
            }

            lsHandler.MovePhaseInterface(
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
            //Make Continuous
            CellMask Near1 = Tracker.Regions.GetNearMask4LevSet(phaseInterface.LevelSetIndex, 1);
            CellMask PosFF = Tracker.Regions.GetLevelSetWing(phaseInterface.LevelSetIndex, +1).VolumeMask;
            enforcer.MakeContinuous(phaseInterface.DGLevelSet, phaseInterface.CGLevelSet, Near1, PosFF);
        }

        public void Initialize(
            IReadOnlyDictionary<string, DGField> ParameterVarFields,
            double time)
        {
            for (int i = 0; i < currentInterfaces.Length; ++i)
            {
                var singleInterface = currentInterfaces[i];
                if (lsHandlers.TryGetValue(singleInterface.CGLevelSet.Identification, out ILevelSetHandler handler))
                {
                    handler.UpdateParameters(
                        singleInterface,
                        Tracker,
                        time,
                        ParameterVarFields);
                }
                else
                {
                    throw new Exception($"LevelSet #{i} does not have a registered handler");
                }
            }
        }
    }
}
