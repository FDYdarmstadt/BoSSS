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
using BoSSS.Application.XNSE_Solver;

namespace BoSSS.Application.XNSFE_Solver {

    public class LevelSetVelocityEvaporative : LevelSetVelocity {

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
            using (new FuncTrace()) {
                // Evaporative part
                BitArray EvapMicroRegion = new BitArray(levelSet.Tracker.GridDat.Cells.Count);  //this.LsTrk.GridDat.GetBoundaryCells().GetBitMask();

                double kA = 0, kB = 0, rhoA = 0, rhoB = 0;

                foreach (var spc in levelSet.Tracker.SpeciesNames) {
                    switch (spc) {
                        case "A": { kA = thermalParameters.k_A; rhoA = thermalParameters.rho_A; break; }
                        case "B": { kB = thermalParameters.k_B; rhoB = thermalParameters.rho_B; break; }
                        default: { throw new ArgumentException("unknown species"); }
                    }
                }

                XDGField temperature = (XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.Temperature];

                int dgDeg = ParameterVarFields[ParameterNames[0]].Basis.Degree;
                int order = dgDeg*dgDeg + 2;
                var quadrule = (new CellQuadratureScheme(false, levelSet.Tracker.Regions.GetCutCellMask())).AddFixedOrderRules(levelSet.Tracker.GridDat, order).Compile(levelSet.Tracker.GridDat, order);

                VectorField<DGField> evapVecVel = new VectorField<DGField>(D.ForLoop(d => ParameterVarFields[ParameterNames[d]]));

                // `GetSpeciesShadowField(...)` is an MPI-Collective function (it calls `MPICollectiveWatchDog.Watch()`) , calling it from inside the integrand might cause some deadlock
                var VelVecA = D.ForLoop(dd => ((XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[dd]]).GetSpeciesShadowField("A"));
                var VelVecB = D.ForLoop(dd => ((XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.VelocityVector(D)[dd]]).GetSpeciesShadowField("B"));
                var temperatureA = temperature.GetSpeciesShadowField("A");
                var temperatureB = temperature.GetSpeciesShadowField("B");

                //for (int d = 0; d < D; d++) {
                //    DGField evapVelocity = ParameterVarFields[ParameterNames[d]];

                evapVecVel.Clear();
                evapVecVel.ProjectField(1.0,
                    delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                        int K = result.GetLength(1); // No nof Nodes

                        MultidimensionalArray VelA = MultidimensionalArray.Create(Len, K, D);
                        MultidimensionalArray VelB = MultidimensionalArray.Create(Len, K, D);


                        for (int dd = 0; dd < D; dd++) {
                            VelVecA[dd].Evaluate(j0, Len, NS, VelA.ExtractSubArrayShallow(new int[] { -1, -1, dd }));
                            VelVecB[dd].Evaluate(j0, Len, NS, VelB.ExtractSubArrayShallow(new int[] { -1, -1, dd }));
                        }

                        MultidimensionalArray GradTempA_Res = MultidimensionalArray.Create(Len, K, D);
                        MultidimensionalArray GradTempB_Res = MultidimensionalArray.Create(Len, K, D);

                        temperatureA.EvaluateGradient(j0, Len, NS, GradTempA_Res);
                        temperatureB.EvaluateGradient(j0, Len, NS, GradTempB_Res);

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

                        temperatureA.Evaluate(j0, Len, NS, TempA_Res);
                        temperatureB.Evaluate(j0, Len, NS, TempB_Res);
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
                                for (int d = 0; d < D; d++) {
                                    double sNeg = VelA[j, k, d] - mEvap * (1 / rhoA) * Normals[j, k, d];
                                    //Console.WriteLine("sNeg = {0}", sNeg);
                                    double sPos = VelB[j, k, d] - mEvap * (1 / rhoB) * Normals[j, k, d];
                                    //Console.WriteLine("sPos = {0}", sPos);

                                    result[j, k, d] = (rhoA * sNeg + rhoB * sPos) / (rhoA + rhoB);   // density averaged evap velocity 
                                }
                            }
                        }
                    }, quadrule);

                //}

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
                    MultidimensionalArray MassFlux_D = MultidimensionalArray.Create(Len, K, D);

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
                        for (int dd = 0; dd < D; dd++) {
                            MassFlux_D.ExtractSubArrayShallow(-1, -1, dd).Acc(-thermalParameters.k_A / thermalParameters.hVap, GradTempA_Res.ExtractSubArrayShallow(-1, -1, dd));
                            MassFlux_D.ExtractSubArrayShallow(-1, -1, dd).Acc(thermalParameters.k_B / thermalParameters.hVap, GradTempB_Res.ExtractSubArrayShallow(-1, -1, dd));
                        }
                    }

                    var Normals = levelSet.Tracker.DataHistories[0].Current.GetLevelSetNormals(NS, j0, Len);

                    for (int j = 0; j < Len; j++) {

                        for (int k = 0; k < K; k++) {
                            
                            bool useNormal = false;
                            double[] sNeg = new double[D];
                            double[] sPos = new double[D];
                            double[] s = new double[D];

                            for (int dd = 0; dd < D; dd++) {
                                double uNeg = VelA[j, k, dd];
                                double uPos = VelB[j, k, dd];
                                if (config.prescribedMassflux != null) {
                                    MultidimensionalArray globCoord = MultidimensionalArray.Create(K, D);
                                    levelSet.Tracker.GridDat.TransformLocal2Global(NS, globCoord, j);
                                    double[] globX = new double[] { globCoord[k, 0], globCoord[k, 1] };

                                    sNeg[dd] = uNeg - (1 / rhoA) * config.prescribedMassflux(globX, time) * Normals[j, k, dd];
                                    sPos[dd] = uPos - (1 / rhoB) * config.prescribedMassflux(globX, time) * Normals[j, k, dd];
                                } else {
                                    if (!useNormal) {
                                        sNeg[dd] = uNeg - (1 / rhoA) * MassFlux_D[j, k, dd];
                                        sPos[dd] = uPos - (1 / rhoB) * MassFlux_D[j, k, dd];
                                    } else {
                                        sNeg[dd] = uNeg - (1 / rhoA) * MassFlux[j, k] * Normals[j, k, dd];
                                        sPos[dd] = uPos - (1 / rhoB) * MassFlux[j, k] * Normals[j, k, dd];
                                    }
                                }

                                switch (this.averagingMode) {
                                    case (XNSE_Control.InterfaceVelocityAveraging.mean): {
                                        s[dd] = 0.5 * (sNeg[dd] + sPos[dd]);
                                        break;
                                    }
                                    case (XNSE_Control.InterfaceVelocityAveraging.phaseA): {
                                        s[dd] = sNeg[dd];
                                        break;
                                    }
                                    case (XNSE_Control.InterfaceVelocityAveraging.phaseB): {
                                        s[dd] = sPos[dd];
                                        break;
                                    }
                                    case (XNSE_Control.InterfaceVelocityAveraging.density):
                                    default: {
                                        s[dd] = (rhoA * sNeg[dd] + rhoB * sPos[dd]) / (rhoA + rhoB);
                                        break;
                                    }
                                }
                            }

                            

                            // project in normal direction
                            if (false) {
                                result[j, k] = Normals[j, k, d];
                                double SxN = 0.0;
                                for(int dd = 0; dd < D; dd++) {
                                    SxN += s[dd] * Normals[j, k, dd];
                                }
                                result[j, k] *= SxN;
                            } else {
                                result[j, k] = s[d];
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
                }, (new CellQuadratureScheme(false, levelSet.Tracker.Regions.GetCutCellMask4LevSet(levelSet.LevelSetIndex))).AddFixedOrderRules(levelSet.Tracker.GridDat, order));
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

