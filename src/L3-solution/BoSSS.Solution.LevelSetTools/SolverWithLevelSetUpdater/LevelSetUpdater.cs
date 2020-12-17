using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using System;
using System.Collections.Generic;

namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {
    public struct DualLevelSet
    {
        public int LevelSetIndex;

        public LevelSet CGLevelSet;

        public LevelSet DGLevelSet;

        public LevelSetTracker Tracker;
    }

    public interface ILevelSetParameter {
        IList<string> ParameterNames { get; }

        (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields);

        void LevelSetParameterUpdate(
            DualLevelSet levelSet,
            double time,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields);
    }

    public interface ILevelSetEvolver {
        IList<string> ParameterNames { get; }

        IList<string> VariableNames { get; }

        void MovePhaseInterface(
            DualLevelSet levelSet,
            double time,
            double dt,
            bool incremental,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields);
    }

    class SingleLevelSetUpdater
    {
        ILevelSetEvolver lsMover;

        ContinuityProjection enforcer;

        DualLevelSet phaseInterface;

        ICollection<ILevelSetParameter> lsParameters;

        public ICollection<ILevelSetParameter> LevelSetParameters
        { 
            get{ return lsParameters; }
        }

        public SingleLevelSetUpdater(DualLevelSet phaseInterface, ContinuityProjection enforcer)
        {
            this.enforcer = enforcer;
            this.phaseInterface = phaseInterface;
            lsParameters = new List<ILevelSetParameter>(10);
        }

        public void SetLevelSetEvolver(ILevelSetEvolver evolver)
        {
            if(lsMover != null)
            {
                throw new Exception("Only one evolver allowed for each levelSet.");
            }
            lsMover = evolver;
        }

        public void AddLevelSetParameter(ILevelSetParameter parameter)
        {
            //Check if already registered
            foreach(ILevelSetParameter registeredParameter in lsParameters)
            {
                foreach(string registeredName in registeredParameter.ParameterNames)
                {
                    foreach(string name in parameter.ParameterNames)
                    {
                        if(registeredName == name)
                        {
                            throw new Exception("Parameter can only be registered once.");
                        }
                    }
                }
            }
            lsParameters.Add(parameter);
        }

        static void UpdateCurrentInterfaces(DualLevelSet phaseInterface)
        {
            DualLevelSet combo = phaseInterface;
            combo.CGLevelSet = (LevelSet)phaseInterface.Tracker.LevelSets[phaseInterface.LevelSetIndex];
            combo.DGLevelSet.Clear();
            combo.DGLevelSet.AccLaidBack(1.0, combo.CGLevelSet);
        }

        public double UpdateLevelSet(
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields,
            double time,
            double dt,
            double underRelax,
            bool incremental)
        {
            UpdateCurrentInterfaces(phaseInterface);
            LevelSet ls = phaseInterface.CGLevelSet;
            LevelSet lsBkUp = ls.CloneAs();

            //Move LevelSet and update Params
            MoveLevelSet(
                phaseInterface,
                DomainVarFields,
                ParameterVarFields,
                time,
                dt,
                underRelax,
                incremental);
            //Calculate Residual
            CellMask oldCC = phaseInterface.Tracker.Regions.GetCutCellMask4LevSet(phaseInterface.LevelSetIndex);
            var newCC = phaseInterface.Tracker.Regions.GetCutCellMask();
            lsBkUp.Acc(-1.0, ls);
            double levSetResidual = lsBkUp.L2Norm(newCC.Union(oldCC));
            return levSetResidual;
        }

        void MoveLevelSet(
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

            lsMover.MovePhaseInterface(
                phaseInterface,
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
            LevelSetTracker Tracker = phaseInterface.Tracker;
            CellMask Near1 = Tracker.Regions.GetNearMask4LevSet(phaseInterface.LevelSetIndex, 1);
            CellMask PosFF = Tracker.Regions.GetLevelSetWing(phaseInterface.LevelSetIndex, +1).VolumeMask;

            enforcer.SetFarField(phaseInterface.DGLevelSet, Near1, PosFF);
            enforcer.MakeContinuous(phaseInterface.DGLevelSet, phaseInterface.CGLevelSet, Near1, PosFF);
        }

        public void UpdateParameters(
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields,
            double time)
        {
            foreach (ILevelSetParameter parameter in lsParameters)
            {
                parameter.LevelSetParameterUpdate(
                    phaseInterface,
                    time,
                    DomainVarFields,
                    ParameterVarFields);
            }
        }
    }

    public class LevelSetUpdater
    {
        public LevelSetTracker Tracker;

        Dictionary<string, SingleLevelSetUpdater> lsUpdaters;

        Dictionary<string, DGField> lsParameterFields;

        public LevelSetUpdater(GridData BackgroundGrid, XQuadFactoryHelper.MomentFittingVariants cutCellquadType, 
            int __NearRegionWidth, string[] _SpeciesTable, LevelSet dgLevelSet, ContinuityProjectionOption continuityMode = ContinuityProjectionOption.ConstrainedDG)
        {
            LevelSet levelSet = ContinuityProjection.CreateField(
                    dgLevelSet, BackgroundGrid, continuityMode);
            levelSet.AccLaidBack(1.0, dgLevelSet);
            
            ContinuityProjection enforcer = new ContinuityProjection(
                levelSet.Basis,
                dgLevelSet.Basis,
                BackgroundGrid,
                continuityMode);
            Tracker = new LevelSetTracker(BackgroundGrid, cutCellquadType, __NearRegionWidth, _SpeciesTable, levelSet);
            DualLevelSet currentInterface = new DualLevelSet
            {
                LevelSetIndex = 0,
                CGLevelSet = levelSet,
                DGLevelSet = dgLevelSet,
                Tracker = Tracker,
            };
            lsUpdaters = new Dictionary<string, SingleLevelSetUpdater>(4);
            lsUpdaters.Add(dgLevelSet.Identification, new SingleLevelSetUpdater( currentInterface, enforcer));
            Tracker.UpdateTracker(0.0);
        }

        public void AddLevelSetParameter(string levelSetName, ILevelSetParameter levelSetParameter)
        {
            if (lsUpdaters.TryGetValue(levelSetName, out SingleLevelSetUpdater mover))
            {
                mover.AddLevelSetParameter(levelSetParameter);
            }
            else
            {
                throw new Exception("LevelSet not registered");
            }
        }

        public void AddEvolver(string levelSetName, ILevelSetEvolver evolver)
        {
            if(lsUpdaters.TryGetValue(levelSetName, out SingleLevelSetUpdater mover))
            {
                mover.SetLevelSetEvolver(evolver);
            }  else
            {
                throw new Exception("LevelSet not registered");
            }
        }

        public double UpdateLevelSets(
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields, 
            double time, 
            double dt, 
            double underRelax,
            bool incremental) 
        {
            var InnerParameterFields = Combine(ParameterVarFields, this.lsParameterFields);
            double residual = 0;
            UpdateParameters(DomainVarFields, InnerParameterFields, time);
            foreach (SingleLevelSetUpdater updater in lsUpdaters.Values) 
            {
                residual += updater.UpdateLevelSet(
                    DomainVarFields,
                    InnerParameterFields,
                    time,
                    dt,
                    underRelax,
                    incremental);
            }
            Tracker.UpdateTracker(time + dt, -1, incremental: true);
            UpdateParameters(DomainVarFields, InnerParameterFields, time + dt);            
            return residual;
        }

        public void InitializeParameters(
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields, 
            double time = 0.0)
        {
            foreach (SingleLevelSetUpdater updater in lsUpdaters.Values)
            {
                InitializeParameters(updater.LevelSetParameters, DomainVarFields, ParameterVarFields);
            }
            var InnerParameterFields = Combine(ParameterVarFields, this.lsParameterFields);
            UpdateParameters(DomainVarFields, InnerParameterFields, 0.0);
        }

        void InitializeParameters(
            ICollection<ILevelSetParameter> parameters,
            IReadOnlyDictionary<string, DGField> DomainVarFields, 
            IReadOnlyDictionary<string, DGField> ParameterVarFields)
            {
            lsParameterFields = new Dictionary<string, DGField>(10);
            foreach (ILevelSetParameter parameter in parameters)
            {
                LinkedList<string> notFound = new LinkedList<string>();
                foreach (string pName in parameter.ParameterNames)
                {
                    if (!ParameterVarFields.ContainsKey(pName))
                    {
                        notFound.AddLast(pName);
                    }
                }
                if(notFound.Count > 0)
                {
                    (string name, DGField field)[] parameterFields = parameter.ParameterFactory(DomainVarFields);
                    while (notFound.Count > 0)
                    {
                        string current = notFound.First.Value;
                        notFound.RemoveFirst();

                        for (int i = 0; i < parameterFields.Length; ++i)
                        {
                            if (parameterFields[i].name == current)
                            {
                                lsParameterFields.Add(parameterFields[i].name, parameterFields[i].field);
                            }
                        }
                    }
                }
            }
        }

        IReadOnlyDictionary<string, DGField> Combine(IReadOnlyDictionary<string, DGField> a, IReadOnlyDictionary<string, DGField> b)
        {
            Dictionary<string, DGField> combination = new Dictionary<string, DGField>(a.Count + b.Count + 5);
            foreach( var entry in a)
            {
                combination.Add(entry.Key, entry.Value);
            }
            foreach (var entry in b)
            {
                combination.Add(entry.Key, entry.Value);
            }
            return combination;
        }

        void UpdateParameters(
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields,
            double time)
        {
            //IReadOnlyDictionary<string, DGField> combinedParameters = Combine(ParameterVarFields, lsParameterFields);
            foreach (SingleLevelSetUpdater updater in lsUpdaters.Values)
            {
                updater.UpdateParameters(DomainVarFields, ParameterVarFields, time);
            }
        }
    }
}
