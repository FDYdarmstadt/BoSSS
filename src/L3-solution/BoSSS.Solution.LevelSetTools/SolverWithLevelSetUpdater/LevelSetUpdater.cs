using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using ilPSP.Utils;
using System;
using System.Collections.Generic;

namespace BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater {
    public struct DualLevelSet
    {
        public string Identification;

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

    class SingleLevelSetUpdater
    {
        ILevelSetEvolver lsMover;

        ContinuityProjection enforcer;

        public DualLevelSet phaseInterface;

        ICollection<ILevelSetParameter> lsParameters;

        public ICollection<ILevelSetParameter> LevelSetParameters {
            get { return lsParameters; }
        }

        public SingleLevelSetUpdater(DualLevelSet phaseInterface, ContinuityProjection enforcer) {
            this.enforcer = enforcer;
            this.phaseInterface = phaseInterface;
            lsParameters = new List<ILevelSetParameter>(10);
        }

        public void SetLevelSetEvolver(ILevelSetEvolver evolver) {
            if (lsMover != null) {
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

        static void UpdateCurrentInterfaces(DualLevelSet phaseInterface) {
            // what for? 
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
            bool incremental){
            
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
            AssertAllRequiredFieldsArePresent(lsMover.ParameterNames, DomainVarFields, ParameterVarFields);
            AssertAllRequiredFieldsArePresent(lsMover.VariableNames, DomainVarFields, ParameterVarFields);
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

            //enforcer.SetFarField(phaseInterface.DGLevelSet, Near1, PosFF);
            enforcer.MakeContinuous(phaseInterface.DGLevelSet, phaseInterface.CGLevelSet, Near1, PosFF);
        }

        public void UpdateParameters(
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields,
            double time)
        {
            foreach (ILevelSetParameter parameter in lsParameters) {
                parameter.LevelSetParameterUpdate(
                    phaseInterface,
                    time,
                    DomainVarFields,
                    ParameterVarFields);
            }
        }

        void AssertAllRequiredFieldsArePresent(IList<string> names, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            if(names!= null) {
                foreach (string name in names) {
                    if (!(DomainVarFields.ContainsKey(name) || ParameterVarFields.ContainsKey(name))) {
                        throw new Exception($"LevelSetUpdater is missing field with name {name}");
                    }
                }
            }
        }
    }

    public class LevelSetUpdater {

        public LevelSetTracker Tracker;

        Dictionary<string, SingleLevelSetUpdater> lsUpdaters;

        Dictionary<string, DGField> lsParameterFields;

        public LevelSetUpdater(GridData backgroundGrid, XQuadFactoryHelper.MomentFittingVariants cutCellquadType,
            int __NearRegionWidth, string[] _SpeciesTable, LevelSet dgLevelSet, string interfaceName) {
            ContinuityProjectionOption continuityMode = ContinuityProjectionOption.ConstrainedDG;
            LevelSet cgLevelSet = ContinuityProjection.CreateField(
                    dgLevelSet, backgroundGrid, continuityMode);
            cgLevelSet.AccLaidBack(1.0, dgLevelSet);
            
            Tracker = new LevelSetTracker(backgroundGrid, cutCellquadType, __NearRegionWidth, _SpeciesTable, cgLevelSet);

            lsUpdaters = new Dictionary<string, SingleLevelSetUpdater>(1);

            DualLevelSet levelSet0 = new DualLevelSet {
                Identification = interfaceName,
                LevelSetIndex = 0,
                CGLevelSet = cgLevelSet,
                DGLevelSet = dgLevelSet,
                Tracker = Tracker,
            };
            SingleLevelSetUpdater singleUpdater = CreateSingleLevelSetUpdater(levelSet0, backgroundGrid, continuityMode);
            lsUpdaters.Add(levelSet0.Identification, singleUpdater);
            
            Tracker.UpdateTracker(0.0);
        }

        public LevelSetUpdater(GridData backgroundGrid, XQuadFactoryHelper.MomentFittingVariants cutCellquadType,
            int __NearRegionWidth, string[,] _SpeciesTable, LevelSet dgLevelSet0, string interfaceName0, LevelSet dgLevelSet1, string interfaceName1) {
            
            ContinuityProjectionOption continuityMode = ContinuityProjectionOption.ConstrainedDG;
            LevelSet[] dgLevelSets = new LevelSet[] { dgLevelSet0, dgLevelSet1 };
            string[] interfaceNames = new string[] { interfaceName0, interfaceName1 };

            LevelSet[] cgLevelSets = new LevelSet[2];

            for (int i = 0; i < 2; ++i) {
                cgLevelSets[i] = ContinuityProjection.CreateField(dgLevelSets[i], backgroundGrid, continuityMode);
                cgLevelSets[i].AccLaidBack(1.0, dgLevelSets[i]);
            }
            Tracker = new LevelSetTracker(backgroundGrid, cutCellquadType, __NearRegionWidth, _SpeciesTable, cgLevelSets[0], cgLevelSets[1]);
            
            lsUpdaters = new Dictionary<string, SingleLevelSetUpdater>(2);
            for(int i = 0; i < 2; ++i) {
                DualLevelSet dualLevelSet = new DualLevelSet {
                    Identification = interfaceNames[i],
                    LevelSetIndex = i,
                    CGLevelSet = cgLevelSets[i],
                    DGLevelSet = dgLevelSets[i],
                    Tracker = Tracker,
                };
                SingleLevelSetUpdater singleUpdater = CreateSingleLevelSetUpdater(dualLevelSet, backgroundGrid, continuityMode);
                lsUpdaters.Add(dualLevelSet.Identification, singleUpdater);
            }

            Tracker.UpdateTracker(0.0);
        }

        public LevelSetUpdater(GridData backgroundGrid, XQuadFactoryHelper.MomentFittingVariants cutCellquadType,
            int __NearRegionWidth, string[,,] _SpeciesTable, LevelSet dgLevelSet0, string interfaceName0, LevelSet dgLevelSet1, string interfaceName1, LevelSet dgLevelSet2, string interfaceName2) {
            ContinuityProjectionOption continuityMode = ContinuityProjectionOption.ConstrainedDG;

            LevelSet[] dgLevelSets = new LevelSet[] { dgLevelSet0, dgLevelSet1, dgLevelSet2 };
            LevelSet[] cgLevelSets = new LevelSet[2];
            string[] interfaceNames = new string[] { interfaceName0, interfaceName1, interfaceName2 };

            for (int i = 0; i < 2; ++i) {
                cgLevelSets[i] = ContinuityProjection.CreateField(dgLevelSets[i], backgroundGrid, continuityMode);
                cgLevelSets[i].AccLaidBack(1.0, dgLevelSets[i]);
            }
            Tracker = new LevelSetTracker(backgroundGrid, cutCellquadType, __NearRegionWidth, _SpeciesTable, cgLevelSets[0], cgLevelSets[1], cgLevelSets[2]);

            lsUpdaters = new Dictionary<string, SingleLevelSetUpdater>(2);
            for (int i = 0; i < 2; ++i) {
                DualLevelSet dualLevelSet = new DualLevelSet {
                    Identification = interfaceNames[i],
                    LevelSetIndex = i,
                    CGLevelSet = cgLevelSets[i],
                    DGLevelSet = dgLevelSets[i],
                    Tracker = Tracker,
                };
                SingleLevelSetUpdater singleUpdater = CreateSingleLevelSetUpdater(dualLevelSet, backgroundGrid, continuityMode);
                lsUpdaters.Add(dualLevelSet.Identification, singleUpdater);
            }

            Tracker.UpdateTracker(0.0);
        }

        public LevelSetUpdater(GridData backgroundGrid, XQuadFactoryHelper.MomentFittingVariants cutCellquadType,
            int __NearRegionWidth, string[,,,] _SpeciesTable, 
            LevelSet dgLevelSet0, string interfaceName0, LevelSet dgLevelSet1, string interfaceName1, LevelSet dgLevelSet2, string interfaceName2, LevelSet dgLevelSet3, string interfaceName3) {
            ContinuityProjectionOption continuityMode = ContinuityProjectionOption.ConstrainedDG;

            LevelSet[] dgLevelSets = new LevelSet[] { dgLevelSet0, dgLevelSet1, dgLevelSet2, dgLevelSet3 };
            LevelSet[] cgLevelSets = new LevelSet[3];
            string[] interfaceNames = new string[] { interfaceName0, interfaceName1, interfaceName2, interfaceName3};

            for (int i = 0; i < 2; ++i) {
                cgLevelSets[i] = ContinuityProjection.CreateField(dgLevelSets[i], backgroundGrid, continuityMode);
                cgLevelSets[i].AccLaidBack(1.0, dgLevelSets[i]);
            }
            Tracker = new LevelSetTracker(backgroundGrid, cutCellquadType, __NearRegionWidth, _SpeciesTable, cgLevelSets[0], cgLevelSets[1], cgLevelSets[2], cgLevelSets[3]);

            lsUpdaters = new Dictionary<string, SingleLevelSetUpdater>(2);
            for (int i = 0; i < 2; ++i) {
                DualLevelSet dualLevelSet = new DualLevelSet {
                    Identification = interfaceNames[i],
                    LevelSetIndex = i,
                    CGLevelSet = cgLevelSets[i],
                    DGLevelSet = dgLevelSets[i],
                    Tracker = Tracker,
                };
                SingleLevelSetUpdater singleUpdater = CreateSingleLevelSetUpdater(dualLevelSet, backgroundGrid, continuityMode);
                lsUpdaters.Add(dualLevelSet.Identification, singleUpdater);
            }

            Tracker.UpdateTracker(0.0);
        }

        static SingleLevelSetUpdater CreateSingleLevelSetUpdater( DualLevelSet levelSet, GridData grid, ContinuityProjectionOption continuityMode) {
            ContinuityProjection enforcer1 = new ContinuityProjection(
                levelSet.CGLevelSet.Basis,
                levelSet.DGLevelSet.Basis,
                grid,
                continuityMode);
            return new SingleLevelSetUpdater(levelSet, enforcer1);
        }

        public IReadOnlyDictionary<string, DGField> Parameters {
            get { return lsParameterFields; }
        }

        public IReadOnlyDictionary<string, DualLevelSet> LevelSets {
            get { 
                var levelSets = new Dictionary<string, DualLevelSet>(lsUpdaters.Count);
                foreach(var lsUpater in lsUpdaters) {
                    string lsName = lsUpater.Key;
                    DualLevelSet ls = lsUpater.Value.phaseInterface;
                    levelSets.Add(lsName, ls);
                }
                return levelSets;
            }
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
            //Tecplot.Tecplot.PlotFields(ArrayTools.Cat(DomainVarFields.Values, InnerParameterFields.Values), "beforeUpdateParameters", time, 2);
            UpdateParameters(DomainVarFields, InnerParameterFields, time);
            //Tecplot.Tecplot.PlotFields(ArrayTools.Cat(DomainVarFields.Values, InnerParameterFields.Values), "afterUpdateParameters", time, 2);
            foreach (SingleLevelSetUpdater updater in lsUpdaters.Values) {
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
            foreach (SingleLevelSetUpdater updater in lsUpdaters.Values){
                InitializeParameters(updater.LevelSetParameters, DomainVarFields, ParameterVarFields);
            }
            var InnerParameterFields = Combine(ParameterVarFields, this.lsParameterFields);
            UpdateParameters(DomainVarFields, InnerParameterFields, 0.0);
            //Tecplot.Tecplot.PlotFields(ArrayTools.Cat(DomainVarFields.Values, InnerParameterFields.Values), "afterUpdateParametersInitial", 0.0, 2);

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
            foreach (SingleLevelSetUpdater updater in lsUpdaters.Values)
            {
                updater.UpdateParameters(DomainVarFields, ParameterVarFields, time);
            }
        }
    }
}
