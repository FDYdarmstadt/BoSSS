using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.Control;
using BoSSS.Solution.Utils;
using BoSSS.Solution.XheatCommon;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using ilPSP.Tracing;
using BoSSS.Solution.NSECommon;
using BoSSS.Solution.LevelSetTools.Advection;
using BoSSS.Solution.LevelSetTools.Reinit.FastMarch;
using BoSSS.Foundation.IO;

namespace BoSSS.Application.XNSE_Solver {
    public static class FromControl {
        public static BeltramiGradient BeltramiGradient(XNSE_Control control, string levelSetName, int D) {
            string levelSet = levelSetName;
            int levelSetDegree;
            if (control.FieldOptions.TryGetValue(levelSet, out FieldOpts lsOpts)) {
                var levelSetSource = control.AdvancedDiscretizationOptions.FilterConfiguration.LevelSetSource;
                levelSetDegree = (levelSetSource == CurvatureAlgorithms.LevelSetSource.fromDG) ? lsOpts.Degree : lsOpts.Degree + 1;
            } else {
                throw new Exception("Level set options not found in FieldOptions");
            }

            DoNotTouchParameters AdvancedDiscretizationOptions = control.AdvancedDiscretizationOptions;
            return new BeltramiGradient(D, AdvancedDiscretizationOptions, levelSetDegree);
        }

        public static BeltramiGradientAndCurvature BeltramiGradientAndCurvature(XNSE_Control control, string levelSetName, int m_HMForder, int D) {
            string curvature = BoSSS.Solution.NSECommon.VariableNames.Curvature;
            int curvatureDegree;
            if (control.FieldOptions.TryGetValue(curvature, out FieldOpts opts)) {
                curvatureDegree = opts.Degree;
            } else {
                throw new Exception("Curvature options not found in FieldOptions");
            }
            
            string levelSet = levelSetName;
            int levelSetDegree;
            if (control.FieldOptions.TryGetValue(levelSet, out FieldOpts lsOpts)) {
                var levelSetSource = control.AdvancedDiscretizationOptions.FilterConfiguration.LevelSetSource;
                levelSetDegree = (levelSetSource == CurvatureAlgorithms.LevelSetSource.fromDG) ? lsOpts.Degree : lsOpts.Degree + 1;
            } else {
                levelSetDegree = 1;
            }
            DoNotTouchParameters AdvancedDiscretizationOptions = control.AdvancedDiscretizationOptions;
            return new BeltramiGradientAndCurvature(curvatureDegree, levelSetDegree, m_HMForder, AdvancedDiscretizationOptions, D);
        }
    }

    /// <summary>
    /// Computation of the fluid interface velocity for material interfaces:
    /// Due to the discontinuous approximation at the interface, 
    /// and the weak enforcement of the velocity jump condition `$ [[\vec{u}]] = 0 `$
    /// the velocities of both phases do not match exactly in the discrete setting.
    /// (The are equal in the continuous setting, however.)
    /// 
    /// Therefore, the phase velocities are averaged according to <see cref="XNSE_Control.InterfaceVelocityAveraging"/>
    /// to obtain a single interface velocity.
    /// </summary>
    public class LevelSetVelocity : ILevelSetParameter {
        protected int D;

        protected IList<string> parameters;

        protected int degree;

        public IList<string> ParameterNames => parameters;

        protected XNSE_Control.InterfaceVelocityAveraging averagingMode;

        protected PhysicalParameters physicalParameters;

        /// <summary>
        /// ctor.
        /// </summary>
        public LevelSetVelocity(string levelSetName, int D, int degree, XNSE_Control.InterfaceVelocityAveraging averagingMode, PhysicalParameters physicalParameters) {
            this.averagingMode = averagingMode;
            this.physicalParameters = physicalParameters;
            this.D = D;
            this.degree = degree;
            parameters = BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D));
        }

        /// <summary>
        /// averaging velocity at interface
        /// </summary>
        public virtual void LevelSetParameterUpdate(
            DualLevelSet levelSet, double time,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            using(new FuncTrace()) {
                int D = levelSet.Tracker.GridDat.SpatialDimension;

                //Mean Velocity
                XDGField[] EvoVelocity; // = new XDGField[]
                try {
                    EvoVelocity = D.ForLoop(
                        d => (XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.Velocity_d(d)]
                        );
                } catch {
                    Console.Error.WriteLine("Velocity not registered as Domainvar, using Velocity from Parametervars");
                    EvoVelocity = D.ForLoop(
                        d => (XDGField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.Velocity0_d(d)]
                        );
                }

                DGField[] meanVelocity = new ConventionalDGField[D];

                double rho_A = physicalParameters.rho_A, rho_B = physicalParameters.rho_B;
                double mu_A = physicalParameters.mu_A, mu_B = physicalParameters.mu_B;
                LevelSetTracker lsTrkr = levelSet.Tracker;
                CellMask CC = lsTrkr.Regions.GetCutCellMask4LevSet(0);
                CellMask Neg = lsTrkr.Regions.GetLevelSetWing(0, -1).VolumeMask;
                CellMask Pos = lsTrkr.Regions.GetLevelSetWing(0, +1).VolumeMask;
                CellMask posNear = lsTrkr.Regions.GetNearMask4LevSet(0, 1).Except(Neg);
                CellMask negNear = lsTrkr.Regions.GetNearMask4LevSet(0, 1).Except(Pos);

                for(int d = 0; d < D; d++) {
                    //Basis b = EvoVelocity[d].Basis.NonX_Basis;
                    meanVelocity[d] = ParameterVarFields[ParameterNames[d]];
                    meanVelocity[d].Clear();


                    foreach(string spc in lsTrkr.SpeciesNames) {
                        double rhoSpc;
                        double muSpc;
                        switch(spc) {
                            case "A": rhoSpc = rho_A; muSpc = mu_A; break;
                            case "B": rhoSpc = rho_B; muSpc = mu_B; break;
                            case "C": continue; // solid phase; ignore
                            default: throw new NotSupportedException("Unknown species name '" + spc + "'");
                        }

                        double scale = 1.0;
                        switch(averagingMode) {
                            case XNSE_Control.InterfaceVelocityAveraging.mean: {
                                scale = 0.5;
                                break;
                            }
                            case XNSE_Control.InterfaceVelocityAveraging.density: {
                                scale = rhoSpc / (rho_A + rho_B);
                                break;
                            }
                            case XNSE_Control.InterfaceVelocityAveraging.viscosity: {
                                scale = muSpc / (mu_A + mu_B);
                                break;
                            }
                            case XNSE_Control.InterfaceVelocityAveraging.phaseA: {
                                scale = (spc == "A") ? 1.0 : 0.0;
                                break;
                            }
                            case XNSE_Control.InterfaceVelocityAveraging.phaseB: {
                                scale = (spc == "B") ? 1.0 : 0.0;
                                break;
                            }
                        }

                        meanVelocity[d].Acc(scale, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(spc), CC);
                        switch(spc) {
                            //case "A": meanVelocity[d].Acc(1.0, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(spc), Neg.Except(CC)); break;
                            case "A": {
                                if(averagingMode != XNSE_Control.InterfaceVelocityAveraging.phaseB)
                                    meanVelocity[d].Acc(1.0, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(spc), negNear);
                                break;
                            }
                            case "B": {
                                if(averagingMode != XNSE_Control.InterfaceVelocityAveraging.phaseA)
                                    meanVelocity[d].Acc(1.0, ((XDGField)EvoVelocity[d]).GetSpeciesShadowField(spc), posNear);
                                break;
                            }
                            default: throw new NotSupportedException("Unknown species name '" + spc + "'");
                        }
                    }
                }
            }
        }

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var velocties = new (string, DGField)[D];
            for(int d = 0; d < D; ++d) {
                Basis basis = new Basis(DomainVarFields.First().Value.GridDat, degree);
                string paramName = ParameterNames[d];
                DGField lsVelocity = new SinglePhaseField(basis, paramName);
                velocties[d] = (paramName, lsVelocity);
            }
            return velocties;
        }
    }


    /// <summary>
    /// Level set velocity, i.e. parameters with name <see cref="BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(string, IList{string})"/>
    /// </summary>
    class ExplicitLevelSetVelocity :  ParameterS, ILevelSetParameter {

        string[] m_ParameterNames;

        public override IList<string> ParameterNames {
            get {
                return m_ParameterNames;
            }
        }

        // delegate (string ParameterName, DGField ParamField)[] DelParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields)
        public override DelParameterFactory Factory => ParameterFactory;

        public ExplicitLevelSetVelocity(string levelSetName, ScalarFunctionTimeDep[] components) : base() {
            this.D = components.Length;
            m_ParameterNames = BoSSS.Solution.NSECommon.VariableNames.AsLevelSetVariable(levelSetName, BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)).ToArray();
            m_components = components.CloneAs();
            
        }

        ScalarFunctionTimeDep[] m_components;
        int D;

        public override DelPartialParameterUpdate Update {
            get {
                return InternalParameterUpdate;
            }
        }

        void InternalParameterUpdate(double t, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            using(new FuncTrace()) {

                DGField[] meanVelocity = new ConventionalDGField[D];
                for(int d = 0; d < D; d++) {
                    //Basis b = EvoVelocity[d].Basis.NonX_Basis;
                    meanVelocity[d] = ParameterVarFields[ParameterNames[d]];
                    meanVelocity[d].Clear();
                    if(m_components[d] != null)
                        meanVelocity[d].ProjectField(m_components[d].SetTime(t));
                }
            }
        }

        public void LevelSetParameterUpdate(DualLevelSet levelSet, double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            InternalParameterUpdate(time, DomainVarFields, ParameterVarFields);
        }

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var velocties = new (string, DGField)[D];
            var bv = DomainVarFields[VariableNames.VelocityX].Basis;
            var b = new Basis(bv.GridDat, bv.Degree);

            for(int d = 0; d < D; ++d) {
                string paramName = ParameterNames[d];
                DGField lsVelocity = new SinglePhaseField(b, paramName);
                velocties[d] = (paramName, lsVelocity);
            }
            return velocties;
        }
    }



    class LevelSetVelocityEvaporative : LevelSetVelocity {

        ThermalParameters thermalParameters;
        XNSFE_OperatorConfiguration config;

        public LevelSetVelocityEvaporative(string levelSetName, int D, int degree, XNSE_Control.InterfaceVelocityAveraging averagingMode, PhysicalParameters physicalParameters, XNSFE_OperatorConfiguration config) : base(levelSetName, D, degree, averagingMode, physicalParameters) {
            this.thermalParameters = config.getThermParams;
            this.config = config;
        }

        public override void LevelSetParameterUpdate(
            DualLevelSet levelSet, double time,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields) {

            // Evaporative part
            BitArray EvapMicroRegion = new BitArray(levelSet.Tracker.GridDat.Cells.Count);  //this.LsTrk.GridDat.GetBoundaryCells().GetBitMask();

            double kA = 0, kB = 0, rhoA = 0, rhoB = 0;

            foreach(var spc in levelSet.Tracker.SpeciesNames) {
                switch (spc) {
                    case "A": { kA = thermalParameters.k_A; rhoA = thermalParameters.rho_A; break; }
                    case "B": { kB = thermalParameters.k_B; rhoB = thermalParameters.rho_B; break; }
                    default: { throw new ArgumentException("unknown species"); }
                }
            }

            XDGField temperature = (XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.Temperature];
            for (int d = 0; d < D; d++) {
                DGField evapVelocity = ParameterVarFields[ParameterNames[d]];

                int order = evapVelocity.Basis.Degree * evapVelocity.Basis.Degree + 2;
                evapVelocity.Clear();
                evapVelocity.ProjectField(1.0,
                   delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                       int K = result.GetLength(1); // No nof Nodes

                       MultidimensionalArray VelA = MultidimensionalArray.Create(Len, K, D);
                       MultidimensionalArray VelB = MultidimensionalArray.Create(Len, K, D);

                       for (int dd = 0; dd < D; dd++) {
                           ((XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[dd]]).GetSpeciesShadowField("A").Evaluate(j0, Len, NS, VelA.ExtractSubArrayShallow(new int[] { -1, -1, dd }));
                           ((XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[dd]]).GetSpeciesShadowField("B").Evaluate(j0, Len, NS, VelB.ExtractSubArrayShallow(new int[] { -1, -1, dd }));
                       }

                       MultidimensionalArray GradTempA_Res = MultidimensionalArray.Create(Len, K, D);
                       MultidimensionalArray GradTempB_Res = MultidimensionalArray.Create(Len, K, D);

                       temperature.GetSpeciesShadowField("A").EvaluateGradient(j0, Len, NS, GradTempA_Res);
                       temperature.GetSpeciesShadowField("B").EvaluateGradient(j0, Len, NS, GradTempB_Res);

                       MultidimensionalArray HeatFluxA_Res = MultidimensionalArray.Create(Len, K, D);
                       MultidimensionalArray HeatFluxB_Res = MultidimensionalArray.Create(Len, K, D);

                       //for (int dd = 0; dd < D; dd++) {
                       //    ((XDGField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.HeatFlux0VectorComponent(dd)]).GetSpeciesShadowField("A").Evaluate(j0, Len, NS, HeatFluxA_Res.ExtractSubArrayShallow(new int[] { -1, -1, dd }));
                       //    ((XDGField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.HeatFlux0VectorComponent(dd)]).GetSpeciesShadowField("B").Evaluate(j0, Len, NS, HeatFluxB_Res.ExtractSubArrayShallow(new int[] { -1, -1, dd }));
                       //}


                       MultidimensionalArray TempA_Res = MultidimensionalArray.Create(Len, K);
                       MultidimensionalArray TempB_Res = MultidimensionalArray.Create(Len, K);
                       //MultidimensionalArray Curv_Res = MultidimensionalArray.Create(Len, K);
                       //MultidimensionalArray Pdisp_Res = MultidimensionalArray.Create(Len, K);

                       temperature.GetSpeciesShadowField("A").Evaluate(j0, Len, NS, TempA_Res);
                       temperature.GetSpeciesShadowField("B").Evaluate(j0, Len, NS, TempB_Res);
                       //ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.Curvature].Evaluate(j0, Len, NS, Curv_Res);
                       //this.DisjoiningPressure.Evaluate(j0, Len, NS, Pdisp_Res);

                       var Normals = levelSet.Tracker.DataHistories[0].Current.GetLevelSetNormals(NS, j0, Len);

                       for (int j = 0; j < Len; j++) {

                           MultidimensionalArray globCoord = MultidimensionalArray.Create(K, D);
                           levelSet.Tracker.GridDat.TransformLocal2Global(NS, globCoord, j);

                           for (int k = 0; k < K; k++) {

                               double qEvap = 0.0;
                               if (EvapMicroRegion[j]) {
                                   throw new NotImplementedException("Check consistency for micro regions");
                                   // micro region
                                   //double Tsat = this.Control.ThermalParameters.T_sat;
                                   //double pc = this.Control.ThermalParameters.pc;
                                   //double pc0 = (pc < 0.0) ? this.Control.PhysicalParameters.Sigma * Curv_Res[j, k] + Pdisp_Res[j, k] : pc;
                                   //double f = this.Control.ThermalParameters.fc;
                                   //double R = this.Control.ThermalParameters.Rc;
                                   //if (this.Control.ThermalParameters.hVap_A > 0) {
                                   //    hVap = this.Control.ThermalParameters.hVap_A;
                                   //    rho_l = this.Control.PhysicalParameters.rho_A;
                                   //    rho_v = this.Control.PhysicalParameters.rho_B;
                                   //    double TintMin = Tsat * (1 + (pc0 / (hVap * rho_l)));
                                   //    double Rint = ((2.0 - f) / (2 * f)) * Tsat * Math.Sqrt(2 * Math.PI * R * Tsat) / (rho_v * hVap.Pow2());
                                   //    if (TempA_Res[j, k] > TintMin)
                                   //        qEvap = -(TempA_Res[j, k] - TintMin) / Rint;
                                   //} else {
                                   //    hVap = -this.Control.ThermalParameters.hVap_A;
                                   //    rho_l = this.Control.PhysicalParameters.rho_B;
                                   //    rho_v = this.Control.PhysicalParameters.rho_A;
                                   //    double TintMin = Tsat * (1 + (pc0 / (hVap * rho_l)));
                                   //    double Rint = ((2.0 - f) / (2 * f)) * Tsat * Math.Sqrt(2 * Math.PI * R * Tsat) / (rho_v * hVap.Pow2());
                                   //    if (TempB_Res[j, k] > TintMin)
                                   //        qEvap = (TempB_Res[j, k] - TintMin) / Rint;
                                   //}

                               } else {
                                   //macro region
                                   for (int dd = 0; dd < D; dd++) {
                                       //qEvap += (HeatFluxB_Res[j, k, dd] - HeatFluxA_Res[j, k, dd]) * Normals[j, k, dd];
                                       qEvap += ((-kA) * GradTempA_Res[j, k, dd] - (-kB) * GradTempB_Res[j, k, dd]) * Normals[j, k, dd];
                                   }
                               }
                               //Console.WriteLine("qEvap delUpdateLevelSet = {0}", qEvap);
                               double[] globX = new double[] { globCoord[k, 0], globCoord[k, 1] };
                               double mEvap = (config.prescribedMassflux != null) ? config.prescribedMassflux(globX, time) : qEvap / thermalParameters.hVap; // mass flux
                                                                                                                                                             //double mEvap = qEvap / this.Control.ThermalParameters.hVap;
                                                                                                                                                             //Console.WriteLine("mEvap - delUpdateLevelSet = {0}", mEvap);
                                                                                                                                                             //Console.WriteLine("prescribedMassFlux = {0}", this.XOpConfig.prescribedMassflux(globX, hack_Phystime));

                               double sNeg = VelA[j, k, d] - mEvap * (1 / rhoA) * Normals[j, k, d];
                               //Console.WriteLine("sNeg = {0}", sNeg);
                               double sPos = VelB[j, k, d] - mEvap * (1 / rhoB) * Normals[j, k, d];
                               //Console.WriteLine("sPos = {0}", sPos);

                               result[j, k] = (rhoA * sNeg + rhoB * sPos) / (rhoA + rhoB);   // density averaged evap velocity 
                           }
                       }
                   }, (new CellQuadratureScheme(false, levelSet.Tracker.Regions.GetCutCellMask())).AddFixedOrderRules(levelSet.Tracker.GridDat, order));

            } 

        }
        
    }

    public class LevelSetVelocityGeneralNonMaterial : LevelSetVelocity {
        ThermalParameters thermalParameters;
        XNSFE_OperatorConfiguration config;
        public LevelSetVelocityGeneralNonMaterial(string levelSetName, int D, int degree, XNSE_Control.InterfaceVelocityAveraging averagingMode, PhysicalParameters physicalParameters, XNSFE_OperatorConfiguration config) : base(levelSetName, D, degree, averagingMode, physicalParameters) {
            this.thermalParameters = config.getThermParams;
            this.config = config;
        }

        public override void LevelSetParameterUpdate(
            DualLevelSet levelSet, double time,
            IReadOnlyDictionary<string, DGField> DomainVarFields,
            IReadOnlyDictionary<string, DGField> ParameterVarFields) {

            double rhoA = 0, rhoB = 0;

            foreach (var spc in levelSet.Tracker.SpeciesNames) {
                switch (spc) {
                    case "A": { rhoA = physicalParameters.rho_A; break; }
                    case "B": { rhoB = physicalParameters.rho_B; break; }
                    case "C": { break; }
                    default: { throw new ArgumentException("unknown species"); }
                }
            }

            for (int d = 0; d < D; d++) {
                DGField interfaceVelocity = ParameterVarFields[ParameterNames[d]];

                int order = interfaceVelocity.Basis.Degree * interfaceVelocity.Basis.Degree + 2;
                interfaceVelocity.Clear();
                interfaceVelocity.ProjectField(1.0,

                delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                    int K = result.GetLength(1); // No nof Nodes

                    MultidimensionalArray VelA = MultidimensionalArray.Create(Len, K, D);
                    MultidimensionalArray VelB = MultidimensionalArray.Create(Len, K, D);
                    MultidimensionalArray MassFlux = MultidimensionalArray.Create(Len, K);
                    MultidimensionalArray MassFlux_d = MultidimensionalArray.Create(Len, K);

                    for (int dd = 0; dd < D; dd++) {
                        ((XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[dd]]).GetSpeciesShadowField("A").Evaluate(j0, Len, NS, VelA.ExtractSubArrayShallow(new int[] { -1, -1, dd }));
                        ((XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[dd]]).GetSpeciesShadowField("B").Evaluate(j0, Len, NS, VelB.ExtractSubArrayShallow(new int[] { -1, -1, dd }));
                    }
                    
                    ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.MassFluxExtension].Evaluate(j0, Len, NS, MassFlux.ExtractSubArrayShallow(new int[] { -1, -1}));

                    if (config.isEvaporation) {
                        XDGField temperature = (XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.Temperature];

                        MultidimensionalArray GradTempA_Res = MultidimensionalArray.Create(Len, K, D);
                        MultidimensionalArray GradTempB_Res = MultidimensionalArray.Create(Len, K, D);

                        temperature.GetSpeciesShadowField("A").EvaluateGradient(j0, Len, NS, GradTempA_Res);
                        temperature.GetSpeciesShadowField("B").EvaluateGradient(j0, Len, NS, GradTempB_Res);

                        MassFlux_d.Acc(-thermalParameters.k_A / thermalParameters.hVap, GradTempA_Res.ExtractSubArrayShallow(-1, -1, d));
                        MassFlux_d.Acc(thermalParameters.k_B / thermalParameters.hVap, GradTempB_Res.ExtractSubArrayShallow(-1, -1, d));
                    }

                    var Normals = levelSet.Tracker.DataHistories[0].Current.GetLevelSetNormals(NS, j0, Len);

                    for (int j = 0; j < Len; j++) {

                        for (int k = 0; k < K; k++) {
                            double uNeg = VelA[j, k, d];
                            double uPos = VelB[j, k, d];
                            bool useNormal = false;
                            double sNeg, sPos;

                            if (config.prescribedMassflux != null) {
                                MultidimensionalArray globCoord = MultidimensionalArray.Create(K, D);
                                levelSet.Tracker.GridDat.TransformLocal2Global(NS, globCoord, j);
                                double[] globX = new double[] { globCoord[k, 0], globCoord[k, 1] };

                                sNeg = uNeg - (1 / rhoA) * config.prescribedMassflux(globX, time) * Normals[j, k, d];
                                sPos = uPos - (1 / rhoB) * config.prescribedMassflux(globX, time) * Normals[j, k, d];
                            } else {
                                if (!useNormal) {
                                    sNeg = uNeg - (1 / rhoA) * MassFlux_d[j, k];
                                    sPos = uPos - (1 / rhoB) * MassFlux_d[j, k];
                                } else {
                                    sNeg = uNeg - (1 / rhoA) * MassFlux[j, k] * Normals[j, k, d];
                                    sPos = uPos - (1 / rhoB) * MassFlux[j, k] * Normals[j, k, d];
                                }
                            }

                            switch (this.averagingMode) {
                                case (XNSE_Control.InterfaceVelocityAveraging.mean):
                                {
                                    result[j, k] = 0.5 * (sNeg + sPos);
                                    break;
                                }
                                case (XNSE_Control.InterfaceVelocityAveraging.phaseA): {
                                    result[j, k] = sNeg;
                                    break;
                                }
                                case (XNSE_Control.InterfaceVelocityAveraging.phaseB): {
                                    result[j, k] = sPos;
                                    break;
                                }
                                case (XNSE_Control.InterfaceVelocityAveraging.density):
                                default: {
                                    result[j, k] = (rhoA * sNeg + rhoB * sPos)/(rhoA + rhoB);
                                    break;
                                }
                            }

                            //if (rhoA != rhoB) {
                            //    result[j, k] = (rhoA * VelA[j, k, d] - rhoB * VelB[j, k, d]) / (rhoA - rhoB);   // interface velocity for arbitrary mass flux
                            //} else if (Math.Abs(VelA[j, k, d] - VelB[j, k, d]) <= 1e-6) {
                            //    double rho = 0.5 * (rhoA + rhoB);                             
                            //    result[j, k] = 0.5 * (VelA[j, k, d] + VelB[j, k, d]) - MassFlux[j, k] / rho; // (relative) interface velocity depends solely on the mass flux
                            //} else {
                            //    throw new ApplicationException("unable to calculate Interface velocity");
                            //}
                        }
                    }
                }, (new CellQuadratureScheme(false, levelSet.Tracker.Regions.GetCutCellMask())).AddFixedOrderRules(levelSet.Tracker.GridDat, order));
            }
        }
    }

    /// <summary>
    /// Update of Massflux Parameter for XNSFE;
    /// Massflux is defined to be positive when pointing in LS normal direction, i.e.
    /// mass flows from - phase to + phase.
    /// Keep in mind maybe it is more stable to only update massflux once per timestep, not in every nonlinear iteration.
    /// </summary>
    public class MassFluxExtension_Evaporation : ParameterS, ILevelSetParameter {

        XNSFE_OperatorConfiguration config;
        DualLevelSet levelSet;
        double time;

        public override IList<string> ParameterNames => new string[] { BoSSS.Solution.NSECommon.VariableNames.MassFluxExtension };

        public override DelParameterFactory Factory => ParameterFactory;

        public override DelPartialParameterUpdate Update {
            get {                
                //return MassFluxExtension_Evaporation_Update;  // seems more stable to update once per timestep
                return null;
            }
        }

        public MassFluxExtension_Evaporation(XNSFE_OperatorConfiguration config) {
            this.config = config;
        }

        public (string, DGField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var massfluxext = new (string, DGField)[1];
            string paramName = BoSSS.Solution.NSECommon.VariableNames.MassFluxExtension;
            Basis b = new Basis(DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityX].Basis.GridDat, DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityX].Basis.Degree);
            DGField MassFluxExtension = new SinglePhaseField(b, paramName);
            massfluxext[0] = (paramName, MassFluxExtension);
            return massfluxext;
        }
        public void MassFluxExtension_Evaporation_Update(double phystime, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            MassFluxExtension_Evaporation_Update(DomainVarFields, ParameterVarFields);
        }
        public void MassFluxExtension_Evaporation_Update(IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            var thermalParams = config.getThermParams;
            double kA = 0, kB = 0;
            foreach (var spc in levelSet.Tracker.SpeciesNames) {
                switch (spc) {
                    case "A": { kA = thermalParams.k_A; break; }
                    case "B": { kB = thermalParams.k_B; break; }
                    case "C": { break; }
                    default: { throw new ArgumentException("unknown species"); }
                }
            }
            var paramName = BoSSS.Solution.NSECommon.VariableNames.MassFluxExtension;

            SinglePhaseField MassFluxField = new SinglePhaseField(ParameterVarFields[paramName].Basis);
            int order = MassFluxField.Basis.Degree * MassFluxField.Basis.Degree + 2;

            XDGField temperature = (XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.Temperature];

            int D = levelSet.Tracker.GridDat.SpatialDimension;
            MassFluxField.Clear();
            MassFluxField.ProjectField(1.0,
                delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                    int K = result.GetLength(1); // No nof Nodes

                    MultidimensionalArray GradTempA_Res = MultidimensionalArray.Create(Len, K, D);
                    MultidimensionalArray GradTempB_Res = MultidimensionalArray.Create(Len, K, D);

                    temperature.GetSpeciesShadowField("A").EvaluateGradient(j0, Len, NS, GradTempA_Res);
                    temperature.GetSpeciesShadowField("B").EvaluateGradient(j0, Len, NS, GradTempB_Res);

                    //MultidimensionalArray HeatFluxA_Res = MultidimensionalArray.Create(Len, K, D);
                    //MultidimensionalArray HeatFluxB_Res = MultidimensionalArray.Create(Len, K, D);

                    //for (int dd = 0; dd < D; dd++) {
                    //    ((XDGField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.HeatFlux0VectorComponent(dd)]).GetSpeciesShadowField("A").Evaluate(j0, Len, NS, HeatFluxA_Res.ExtractSubArrayShallow(new int[] { -1, -1, dd }));
                    //    ((XDGField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.HeatFlux0VectorComponent(dd)]).GetSpeciesShadowField("B").Evaluate(j0, Len, NS, HeatFluxB_Res.ExtractSubArrayShallow(new int[] { -1, -1, dd }));
                    //}

                    var Normals = levelSet.Tracker.DataHistories[0].Current.GetLevelSetNormals(NS, j0, Len);

                    for (int j = 0; j < Len; j++) {

                        MultidimensionalArray globCoord = MultidimensionalArray.Create(K, D);
                        levelSet.Tracker.GridDat.TransformLocal2Global(NS, globCoord, j);

                        for (int k = 0; k < K; k++) {

                            double qEvap = 0.0;
                            //macro region                            
                            for (int dd = 0; dd < D; dd++) {
                                qEvap += ((-kA) * GradTempA_Res[j, k, dd] - (-kB) * GradTempB_Res[j, k, dd]) * Normals[j, k, dd];
                                //qEvap += (HeatFluxB_Res[j, k, dd] - HeatFluxA_Res[j, k, dd]) * Normals[j, k, dd];
                            }

                            //Console.WriteLine("qEvap delUpdateLevelSet = {0}", qEvap);
                            double[] globX = new double[] { globCoord[k, 0], globCoord[k, 1] };
                            double mEvap = (config.prescribedMassflux != null) ? config.prescribedMassflux(globX, time) : qEvap / thermalParams.hVap; // mass flux

                            //Console.WriteLine("mEvap - delUpdateLevelSet = {0}", mEvap);
                            result[j, k] = mEvap;

                        }
                    }
                }, (new CellQuadratureScheme(false, levelSet.Tracker.Regions.GetCutCellMask())).AddFixedOrderRules(levelSet.Tracker.GridDat, order));


            // no extension
            ParameterVarFields[paramName].Clear();
            ParameterVarFields[paramName].Acc(1.0, MassFluxField);

            // extension
            //SubGrid CCgrid = levelSet.Tracker.Regions.GetCutCellSubGrid();
            //CellMask CC = levelSet.Tracker.Regions.GetCutCellMask();
            //CellMask NEAR = levelSet.Tracker.Regions.GetNearFieldMask(1);
            //int J = this.levelSet.Tracker.GridDat.Cells.NoOfLocalUpdatedCells;
            //double[][] MassFluxMin = new double[1][];
            //double[][] MassFluxMax = new double[1][];
            //MassFluxMin[0] = new double[J];
            //MassFluxMax[0] = new double[J];

            //VectorField<SinglePhaseField> DGLevSetGradient = new VectorField<SinglePhaseField>(D.ForLoop(d => new SinglePhaseField(levelSet.DGLevelSet.Basis)));
            //DGLevSetGradient.Gradient(1.0, levelSet.DGLevelSet);

            //NarrowMarchingBand.ConstructExtVel_PDE(levelSet.Tracker, CCgrid, new SinglePhaseField[] { (SinglePhaseField)ParameterVarFields[paramName] }, new SinglePhaseField[] { MassFluxField },
            //    levelSet.DGLevelSet, DGLevSetGradient, MassFluxMin, MassFluxMax, order);

            //var marcher = new FastMarchReinit(levelSet.DGLevelSet.Basis);
            //marcher.ConstructExtension(levelSet.DGLevelSet, NEAR.Except(CC), CC, new SinglePhaseField[] { (SinglePhaseField)ParameterVarFields[paramName] },
            //    MassFluxMin, MassFluxMax, DGLevSetGradient, 0);

            ParameterVarFields[paramName].CheckForNanOrInf(true, true, true);
        }

        public void LevelSetParameterUpdate(DualLevelSet levelSet, double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            this.levelSet = levelSet;
            this.time = time;

            MassFluxExtension_Evaporation_Update(DomainVarFields, ParameterVarFields);
        }
    }

    class Temperature0 : ParameterS {
        public override IList<string> ParameterNames => new string[] { BoSSS.Solution.NSECommon.VariableNames.Temperature0 };

        public override DelParameterFactory Factory => Temperature0Factory;

        public Temperature0() {
        }

        (string, DGField)[] Temperature0Factory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var temperature0 = new (string, DGField)[1];
            
            string temperaturename = BoSSS.Solution.NSECommon.VariableNames.Temperature;
            DGField temperature = DomainVarFields[temperaturename];
            string paramName = BoSSS.Solution.NSECommon.VariableNames.Temperature0;
            temperature0[0] = (paramName, temperature);
            
            return temperature0;
        }
    }

    class HeatFlux0 : ParameterS {

        int D;
        string[] parameters;
        LevelSetTracker lstrk;
        Dictionary<string, double> k_spc;
        public override IList<string> ParameterNames => parameters;

        public override DelParameterFactory Factory => HeatFlux0Factory;

        public HeatFlux0(int D, LevelSetTracker lstrk, XNSFE_OperatorConfiguration config) {
            this.D = D;
            this.lstrk = lstrk;
            parameters = BoSSS.Solution.NSECommon.VariableNames.HeatFlux0Vector(D);
            k_spc = new Dictionary<string, double>();
            var thermParams = config.getThermParams;
            foreach(var spc in lstrk.SpeciesNames) {
                switch (spc) {
                    case "A": { k_spc.Add(spc, thermParams.k_A); break; }
                    case "B": { k_spc.Add(spc, thermParams.k_B); break; }
                    default: { throw new ArgumentException("Unknown species.");}                    
                }
            }
            config_conductMode = config.conductMode;
            //if(config.conductMode == ConductivityInSpeciesBulk.ConductivityMode.SIP)
            //    Update = HeatFlux0Update;
        }

        ConductivityInSpeciesBulk.ConductivityMode config_conductMode;

        public override DelPartialParameterUpdate Update {
            get {
                if(config_conductMode == ConductivityInSpeciesBulk.ConductivityMode.SIP)
                    return HeatFlux0Update;
                else
                    return null;

            }
        }
        (string, DGField)[] HeatFlux0Factory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var heatflux0 = new (string, DGField)[D];

            for (int d = 0; d < D; d++) {

                string heatfluxname = BoSSS.Solution.NSECommon.VariableNames.HeatFluxVectorComponent(d);
                string paramName = BoSSS.Solution.NSECommon.VariableNames.HeatFlux0VectorComponent(d);

                XDGField heatflux;
                if (DomainVarFields.ContainsKey(heatfluxname)) {
                    // if the heatflux is a domainvariable, we can directly use him, no special update necessary
                    heatflux = (XDGField)DomainVarFields[heatfluxname];
                } else {
                    heatflux = new XDGField((XDGBasis)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.Temperature].Basis, paramName);
                }
                heatflux0[d] = (paramName, heatflux);
            }          

            return heatflux0;
        }
        private void HeatFlux0Update(double t, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            var varname = BoSSS.Solution.NSECommon.VariableNames.Temperature;
            DGField temperaure = DomainVarFields[varname];

            SubGrid Sgrd = lstrk.Regions.GetCutCellSubGrid();

            for (int d = 0; d < D; d++) {
                foreach(var spc in lstrk.SpeciesNames) {
                    
                    DGField temperature_Spc = ((temperaure as XDGField).GetSpeciesShadowField(spc));

                    double aSpc = (k_spc != null && k_spc.ContainsKey(spc)) ? k_spc[spc] : 1.0;

                    var paramname = BoSSS.Solution.NSECommon.VariableNames.HeatFlux0VectorComponent(d);
                    DGField heatflux = ParameterVarFields[paramname];
                    heatflux.Clear();
                    // q=-k*grad(T)
                    (heatflux as XDGField).GetSpeciesShadowField(spc).DerivativeByFlux(-aSpc, temperature_Spc, d, Sgrd);
                }
            }            
        }

    }

    /// <summary>
    /// Heat source Parameter, note a few specials:
    /// 1.  The sign of this heat source set in the controlfile is that of resultant source. 
    ///     However here we need to project that field by factor -1. I.e. when looking at the field it is exactly opposite in sign to that speciified in the controlfile.
    ///     This is due to the implementation, where the source term is brought to the left-hand-side of the  Heat Equation
    ///     <see cref="Solution.XNSECommon.Operator.MultiPhaseSource"/>
    /// </summary>
    public class HeatSource : ParameterS {
        int degree;

        ScalarFunctionTimeDep initial;

        string[] names;

        public override DelParameterFactory Factory => ParameterFactory;

        public override DelPartialParameterUpdate Update => ParameterUpdate;

        public HeatSource(string species, ScalarFunctionTimeDep initial, int degree) {
            this.degree = degree;

            names = new string[1];
            string source = BoSSS.Solution.NSECommon.VariableNames.HeatSource;
            names[0] = source + "#" + species;
            this.initial = initial;
        }

        public static HeatSource CreateFrom(string species, AppControl control, ScalarFunctionTimeDep sourceFunc) {
            string source = BoSSS.Solution.NSECommon.VariableNames.HeatSource;
            string sourceOfSpecies = source + "#" + species;            

            int sourceDegree;
            if (control.FieldOptions.TryGetValue(source, out FieldOpts opts)) {
                sourceDegree = Math.Max(0, opts.Degree);
            } else if (control.FieldOptions.TryGetValue("Velocity*", out FieldOpts velOpts)) {
                sourceDegree = velOpts.Degree;
            } else {
                sourceDegree = 0;
            }

            return new HeatSource(species, sourceFunc, sourceDegree);
        }

        public override IList<string> ParameterNames => names;

        public (string, DGField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            Basis basis = new Basis(DomainVarFields.First().Value.GridDat, degree);
            DGField source = new SinglePhaseField(basis, names[0]);
            source.Clear();
            if (initial != null)
                source.ProjectField(-1.0, initial.SetTime(0.0));

            return new (string, DGField)[] { (names[0], source) };
        }

        public void ParameterUpdate(double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            if (initial != null) {
                DGField source = ParameterVarFields[names[0]];
                source.Clear();
                source.ProjectField(-1.0, initial.SetTime(time));
            }
        }
    }
}

