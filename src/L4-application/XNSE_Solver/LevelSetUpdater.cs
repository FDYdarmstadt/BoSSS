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
        
        Dictionary<string, ILevelSetMover> lsMovers;

        Dictionary<string, ICollection<ILevelSetParameter>> lsParameterUpdaters;

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
            lsMovers = new Dictionary<string, ILevelSetMover>(4);
            lsParameterUpdaters = new Dictionary<string, ICollection<ILevelSetParameter>>(4);
            Tracker.UpdateTracker(0.0);
        }

        public void AddLevelSetParameter(string levelSetName, ILevelSetParameter updater)
        {
            if(lsParameterUpdaters.TryGetValue(levelSetName, out ICollection<ILevelSetParameter> parameters))
            {
                parameters.Add(updater);
            }
            else
            {
                ICollection<ILevelSetParameter> parameterCollection = new LinkedList<ILevelSetParameter>();
                parameterCollection.Add(updater);
                lsParameterUpdaters.Add(levelSetName, parameterCollection);
            }
        }

        public void AddEvolver(string levelSetName, ILevelSetMover evolver)
        {
            lsMovers.Add(levelSetName, evolver);
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
            for (int i = 0; i < currentInterfaces.Length; ++i)
            {
                var singleInterface = currentInterfaces[i];
                var enforcer = enforcers[i];
                if (lsMovers.TryGetValue(singleInterface.CGLevelSet.Identification, out ILevelSetMover handler))
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
                    //Calculate Residual
                    CellMask oldCC = Tracker.Regions.GetCutCellMask4LevSet(singleInterface.LevelSetIndex);
                    var newCC = Tracker.Regions.GetCutCellMask();
                    lsBkUp.Acc(-1.0, ls);
                    double levSetResidual = lsBkUp.L2Norm(newCC.Union(oldCC));
                    residual += levSetResidual;
                }
            }
            UpdateParameters(ParameterVarFields, time + dt);
            Tracker.UpdateTracker(time + dt, -1, incremental: true);
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
            ILevelSetMover lsHandler,
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
            
            enforcer.SetFarField(phaseInterface.DGLevelSet, Near1, PosFF);
            enforcer.MakeContinuous(phaseInterface.DGLevelSet, phaseInterface.CGLevelSet, Near1, PosFF);
        }

        public void UpdateParameters(
            IReadOnlyDictionary<string, DGField> ParameterVarFields,
            double time)
        {
            for (int i = 0; i < currentInterfaces.Length; ++i)
            {
                var singleInterface = currentInterfaces[i];
                if (lsParameterUpdaters.TryGetValue(singleInterface.CGLevelSet.Identification, out ICollection<ILevelSetParameter> parameters))
                {
                    UpdateLevelSetParameter(parameters, singleInterface, ParameterVarFields, time);
                }
                else
                {
                    Console.WriteLine($"Warning: LevelSet #{i} does not have a registered updater");
                }
            }
        }

        void UpdateLevelSetParameter(ICollection<ILevelSetParameter> parameters, DualLevelSet levelSet, IReadOnlyDictionary<string, DGField> ParameterVarFields, double time)
        {
            foreach (ILevelSetParameter parameter in parameters)
            {
                bool parametermatch = true;
                foreach (string pName in parameter.ParameterNames)
                {
                    if (!ParameterVarFields.ContainsKey(pName))
                    {
                        parametermatch = false;
                        break;
                    }
                }
                if (parametermatch)
                {
                    parameter.UpdateParameters(
                    levelSet,
                    Tracker,
                    time,
                    ParameterVarFields);
                }
                else
                {
                    string names = "";
                    foreach (string name in parameter.ParameterNames)
                    {
                        names += name + ",";
                    }
                    Console.WriteLine($"Warning: LevelSet Parameters {names}  not found in global Parameters.");
                }
            }
        }
    }
}
