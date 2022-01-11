using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.LevelSetTools.SolverWithLevelSetUpdater;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;

namespace BoSSS.Application.XNSEC {

    internal class ThermodynamicPressure : ParameterS {
        private ThermodynamicPressureMode m_mode;
        private double m_initialMass = 1.0;
        private MaterialLaw_MultipleSpecies m_EoS;

        public ThermodynamicPressure(double initialMass, ThermodynamicPressureMode mode, MaterialLaw EoS) {
            this.m_initialMass = initialMass;
            this.m_mode = mode;
            this.m_EoS = (MaterialLaw_MultipleSpecies)EoS;
        }

        public override IList<string> ParameterNames => new string[] { BoSSS.Solution.NSECommon.VariableNames.ThermodynamicPressure };

        public override DelParameterFactory Factory => ThermodynamicPressureFactory;

        public

        (string, DGField)[] ThermodynamicPressureFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var ThermodynamicPressure0 = new (string, DGField)[1];
            string ThermodynamicPressureName = BoSSS.Solution.NSECommon.VariableNames.ThermodynamicPressure;
            //var grd = DomainVarFields[VariableNames.Temperature].GridDat;
            //     var b = new Basis(grd, 0); // polynomial degree zero
            //var ThermodynamicPressure = new SinglePhaseField(b, ThermodynamicPressureName); // Should i define it as singlephasefield or XDGfield?
            var ThermodynamicPressure = new XDGField((XDGBasis)DomainVarFields[VariableNames.Temperature].Basis, ThermodynamicPressureName);
            ThermodynamicPressure0[0] = (ThermodynamicPressureName, ThermodynamicPressure);

            if(!m_EoS.IsInitialized) {
                XDGField Temperature = (XDGField)DomainVarFields[VariableNames.Temperature];
                double init_p0 = m_EoS.GetMassDeterminedThermodynamicPressure(m_initialMass, Temperature);

                m_EoS.Initialize(init_p0);
            }

            return ThermodynamicPressure0;
        }

        public override DelPartialParameterUpdate Update => ThermodynamicPressureUpdate;

        private void ThermodynamicPressureUpdate(double t, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            XDGField ThermodynamicPressure = (XDGField)ParameterVarFields[VariableNames.ThermodynamicPressure];
            XDGField Temperature = (XDGField)DomainVarFields[VariableNames.Temperature];
            ThermodynamicPressure.Clear();

            double p0;
            switch(m_mode) {
                case ThermodynamicPressureMode.Constant:
                p0 = 1.0;
                m_EoS.ThermodynamicPressure = p0;
                break;

                case ThermodynamicPressureMode.MassDetermined:
                p0 = m_EoS.GetMassDeterminedThermodynamicPressure(m_initialMass, Temperature);
                //Console.WriteLine("Calculated Thermodynamic Pressure: " + p0);
                //m_EoS.ThermodynamicPressure.Current.Clear();
                //m_EoS.ThermodynamicPressure.Current.AccConstant(p0);
                m_EoS.ThermodynamicPressure = p0;
                break;

                default:
                throw new NotImplementedException();
            }
            //Console.WriteLine("Updated Thermodynamic Pressure value: " + p0);
            ThermodynamicPressure.Clear();
            ThermodynamicPressure.AccConstant(p0); //////// How should p0 be defined in a multiphase setting?

            //var p0actual = ParameterVarFields[VariableNames.ThermodynamicPressure];
            //var p0Old= ParameterVarFields[VariableNames.ThermodynamicPressure + "_t0"];

            //p0actual.Acc(-1.0, p0Old);

            //Console.WriteLine( "p0-p0_0: " + p0actual.GetMeanValueTotal(null) );
        }
    }

    //internal class ThermodynamicPressure_Oldtimestep : ParameterS, ILevelSetParameter {
    //    private ThermodynamicPressureMode m_mode;
    //    private double m_initialMass = 1.0;
    //    private MaterialLawMultiSpecies m_EoS;

    //    public ThermodynamicPressure_Oldtimestep(double initialMass, ThermodynamicPressureMode mode, MaterialLawMultiSpecies EoS) {
    //        this.m_initialMass = initialMass;
    //        this.m_mode = mode;
    //        this.m_EoS = EoS;
    //        paramName = VariableNames.ThermodynamicPressure + "_t0";
    //    }

    //    private string paramName;
    //    public override IList<string> ParameterNames => new string[] { paramName };

    //    public override DelParameterFactory Factory => ParameterFactory;

    //    public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
    //        var ThermodynamicPressure0_Old = new (string, DGField)[1];
    //        string ThermodynamicPressureName = paramName;
    //        //var grd = DomainVarFields[VariableNames.Temperature].GridDat;
    //        //     var b = new Basis(grd, 0); // polynomial degree zero
    //        //var ThermodynamicPressure = new SinglePhaseField(b, ThermodynamicPressureName); // Should i define it as singlephasefield or XDGfield?
    //        var ThermodynamicPressureOld = new XDGField((XDGBasis)DomainVarFields[VariableNames.Temperature].Basis, ThermodynamicPressureName);
    //        ThermodynamicPressure0_Old[0] = (ThermodynamicPressureName, ThermodynamicPressureOld);
    //        return ThermodynamicPressure0_Old;
    //    }

    //    public void LevelSetParameterUpdate(DualLevelSet levelSet, double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
    //        XDGField ThermodynamicPressureOld = (XDGField)ParameterVarFields[paramName];

    //        double p0;
    //        switch(m_mode) {
    //            case ThermodynamicPressureMode.Constant:
    //            p0 = 1.0;
    //            m_EoS.ThermodynamicPressure = p0;
    //            break;

    //            case ThermodynamicPressureMode.MassDetermined:
    //            //XDGField Temperature = (XDGField)DomainVarFields[VariableNames.Temperature];
    //            //p0 = m_EoS.GetMassDeterminedThermodynamicPressure(m_initialMass, Temperature);
    //            p0 = m_EoS.ThermodynamicPressure; // Should give the same result as calculating it again...
    //            break;

    //            default:
    //            throw new NotImplementedException();
    //        }

    //        ThermodynamicPressureOld.Clear();
    //        ThermodynamicPressureOld.AccConstant(p0); //////// How should p0 be defined in a multiphase setting?
    //    }
    //}

    internal class ThermodynamicPressure_Oldtimestep : ParameterS {
        private ThermodynamicPressureMode m_mode;
        private double m_initialMass = 1.0;
        private MaterialLaw_MultipleSpecies m_EoS;
        private double t_previous = -1;

        public ThermodynamicPressure_Oldtimestep(double initialMass, ThermodynamicPressureMode mode, MaterialLaw_MultipleSpecies EoS) {
            this.m_initialMass = initialMass;
            this.m_mode = mode;
            paramName = VariableNames.ThermodynamicPressure + "_t0";
            this.m_EoS = (MaterialLaw_MultipleSpecies)EoS;
        }

        private string paramName;
        public override IList<string> ParameterNames => new string[] { paramName };

        public override DelParameterFactory Factory => ThermodynamicPressure_Oldtimestep_Factory;

        public (string ParameterName, DGField ParamField)[] ThermodynamicPressure_Oldtimestep_Factory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var ThermodynamicPressure0_Old = new (string, DGField)[1];
            string ThermodynamicPressureName = paramName;
            //var grd = DomainVarFields[VariableNames.Temperature].GridDat;
            //     var b = new Basis(grd, 0); // polynomial degree zero
            //var ThermodynamicPressure = new SinglePhaseField(b, ThermodynamicPressureName); // Should i define it as singlephasefield or XDGfield?
            var ThermodynamicPressureOld = new XDGField((XDGBasis)DomainVarFields[VariableNames.Temperature].Basis, ThermodynamicPressureName);
            ThermodynamicPressure0_Old[0] = (ThermodynamicPressureName, ThermodynamicPressureOld);
            return ThermodynamicPressure0_Old;
        }

        public override DelPartialParameterUpdate Update => OldThermodynamicPressureUpdate;

        private void OldThermodynamicPressureUpdate(double t, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            //if(this.t_previous != t && t != 0.0) { // Update only once at the beggining from each timestep.
            //    XDGField ThermodynamicPressure = (XDGField)ParameterVarFields[VariableNames.ThermodynamicPressure];
            //    XDGField ThermodynamicPressure_old = (XDGField)ParameterVarFields[VariableNames.ThermodynamicPressure + "_t0"];
            //    ThermodynamicPressure_old.Clear();
            //    ThermodynamicPressure_old.Acc(1.0, ThermodynamicPressure);
            //    Console.WriteLine("p0 : {0}", ThermodynamicPressure_old.GetMeanValueTotal(null));

            //}
            //this.t_previous = t;
        }
    }

    /// <summary>
    /// Time derivative of the thermodynamic pressure used in the energy equation for closed systems.
    /// See A comparative study of high-order variable-property segregated
    /// algorithms for unsteady low Mach number flows,    R.Knikker∗
    /// </summary>
    internal class dp0dt : ParameterS {
        private string paramName;
        private MaterialLaw_MultipleSpecies m_EoS;
        private double m_Reynolds;
        private double m_Prandtl;

        public dp0dt(MaterialLaw_MultipleSpecies EoS, double Reynolds, double Prandtl) {
            paramName = "dp0dt";
            this.m_EoS = (MaterialLaw_MultipleSpecies)EoS;
            m_Reynolds = Reynolds;
            m_Prandtl = Prandtl;
        }

        public override IList<string> ParameterNames => new string[] { paramName };

        public override DelParameterFactory Factory => dp0dt_Factory;

        public (string ParameterName, DGField ParamField)[] dp0dt_Factory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var dp0dt = new (string, DGField)[1];
            Basis b = new Basis(DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.Temperature].Basis.GridDat, DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.Temperature].Basis.Degree);
            DGField dp0dt_field = new SinglePhaseField(b, paramName);
            //DGField dp0dt_field =                 new XDGField((XDGBasis)DomainVarFields[VariableNames.Temperature].Basis, paramName);
            dp0dt[0] = (paramName, dp0dt_field);
            return dp0dt;
        }

        public override DelPartialParameterUpdate Update => dp0dt_update;

        private void dp0dt_update(double t, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            //dp0dt can be calculated as
            //dp0dt = \frac{\gamma - 1 }{V} \int{\frac{\partial}{\partial x_j}\left(k * \frac{\partial T}{\partial x_j} \right)    }_V

            var temperature = (XDGField)DomainVarFields[VariableNames.Temperature];

            //var aux_field_0 = (XDGField)ParameterVarFields["dp0dt"].CloneAs();
            //aux_field_0.Clear();
            //var aux_field_1 = (XDGField)ParameterVarFields["dp0dt"].CloneAs();
            //aux_field_1.Clear();
            //var auxVector = new XDGField[] { aux_field_0, aux_field_1 };
            //int D = 2;

            SinglePhaseField aux_field_0 = new SinglePhaseField(ParameterVarFields["dp0dt"].Basis);
            SinglePhaseField aux_field_1 = new SinglePhaseField(ParameterVarFields["dp0dt"].Basis);

            var auxVector = new SinglePhaseField[] { aux_field_0, aux_field_1 };
            int D = 2;

            //================================
            //Calculation of k*dT/dx and k*dT/dy
            //================================
            for(int i = 0; i < D; i++) {
                auxVector[i].ProjectField(1.0,
                  delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                      int K = result.GetLength(1); // No of Nodes
                      MultidimensionalArray tempA = MultidimensionalArray.Create(Len, K);
                      MultidimensionalArray tempB = MultidimensionalArray.Create(Len, K);
                      temperature.GetSpeciesShadowField("A").Evaluate(j0, Len, NS, tempA);
                      temperature.GetSpeciesShadowField("B").Evaluate(j0, Len, NS, tempB);
                      MultidimensionalArray GradTempA_Res = MultidimensionalArray.Create(Len, K, D);
                      MultidimensionalArray GradTempB_Res = MultidimensionalArray.Create(Len, K, D);
                      temperature.GetSpeciesShadowField("A").EvaluateGradient(j0, Len, NS, GradTempA_Res);
                      temperature.GetSpeciesShadowField("B").EvaluateGradient(j0, Len, NS, GradTempB_Res);
                      for(int j = 0; j < Len; j++) {
                          for(int k = 0; k < K; k++) {
                              double lambda = m_EoS.GetHeatConductivity(tempA[j, k]);
                              result[j, k] = lambda * GradTempA_Res[j, k, i];
                          }
                      }
                  }, new Foundation.Quadrature.CellQuadratureScheme(true, null));
            }

            var auxfield_0_x = aux_field_0.CloneAs();
            auxfield_0_x.Clear();
            var auxfield_1_y = aux_field_1.CloneAs();
            auxfield_1_y.Clear();

            auxfield_0_x.Derivative(1.0, aux_field_0, 0);
            auxfield_1_y.Derivative(1.0, aux_field_1, 1);

            var sum_auxFields = ParameterVarFields["dp0dt"].CloneAs();
            sum_auxFields.Clear();
            sum_auxFields.Acc(1.0, auxfield_0_x);
            sum_auxFields.Acc(1.0, auxfield_1_y);
            double integralTotal = sum_auxFields.IntegralOver(null);

            //================================
            //Calculation of the total volume
            //================================

            CellMask cellMask = CellMask.GetFullMask(temperature.GridDat);
            double volumeTotal = 0;
            foreach(int cell in cellMask.ItemEnum) {
                double cellvolume = temperature.GridDat.iGeomCells.GetCellVolume(cell);
                volumeTotal += cellvolume;
            }

            double gamma = 1.4;
            double dp0dtValue = (gamma - 1) / volumeTotal * integralTotal;

            dp0dtValue *= 1.0 / (m_Reynolds * m_Prandtl); // should be scaled?? todo revisar

            var dp0dtfield = ParameterVarFields["dp0dt"];
            dp0dtfield.Clear();
            dp0dtfield.AccConstant(dp0dtValue);

            //Console.WriteLine(integralTotal);
        }
    }

    internal class Density_t0 : ParameterS/*, ILevelSetParameter*/ {
        private string paramName;
        private int m_NumberOfChemicalSpecies;
        private MaterialLaw_MultipleSpecies m_EoS;

        public Density_t0(int NumberOfChemicalSpecies, MaterialLaw_MultipleSpecies EoS) {
            paramName = "Density_t0";
            m_NumberOfChemicalSpecies = NumberOfChemicalSpecies;
            m_EoS = EoS;
        }

        public override IList<string> ParameterNames => new string[] { paramName };

        public void LevelSetParameterUpdate(DualLevelSet levelSet, double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            string[] Fields = ArrayTools.Cat(new string[] { VariableNames.Temperature }, VariableNames.MassFractions(m_NumberOfChemicalSpecies));
            XDGField[] allFields = new XDGField[Fields.Length];
            for(int i = 0; i < Fields.Length; i++) {
                allFields[i] = (XDGField)DomainVarFields[Fields[i]].CloneAs();
            }

            XDGField rhoOldTimeStep = (XDGField)ParameterVarFields[paramName];

            var p0Field = (XDGField)ParameterVarFields[VariableNames.ThermodynamicPressure];
            double p0 = p0Field.GetSpeciesShadowField("A").GetMeanValueTotal(null); // should already contain the "correct" value from the last timestep...

            //if(p0 <= 0)
            //    Console.WriteLine("??");

            rhoOldTimeStep.Clear();
            //rhoOldTimeStep.GetSpeciesShadowField("A").ProjectField(1.0,
            //                   delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
            //                       int K = result.GetLength(1);
            //                       MultidimensionalArray tempT = MultidimensionalArray.Create(Len, K);
            //                       var tempMFs = new MultidimensionalArray[m_NumberOfChemicalSpecies];

            //                       allFields[0].GetSpeciesShadowField("A").Evaluate(j0, Len, NS, tempT);
            //                       for(int i = 0; i < m_NumberOfChemicalSpecies; i++) {
            //                           tempMFs[i] = MultidimensionalArray.Create(Len, K);
            //                           allFields[i + 1].GetSpeciesShadowField("A").Evaluate(j0, Len, NS, tempMFs[i]);
            //                       }

            //                       for(int j = 0; j < Len; j++) {
            //                           for(int k = 0; k < K; k++) {
            //                               double[] densityArguments = new double[1 + m_NumberOfChemicalSpecies];
            //                               densityArguments[0] = tempT[j, k];
            //                               for(int i = 0; i < m_NumberOfChemicalSpecies; i++) {
            //                                   densityArguments[i + 1] = tempMFs[i][j, k];
            //                               }
            //                               result[j, k] = m_EoS.GetDensity(p0, densityArguments);
            //                           }
            //                       }
            //                   }, new Foundation.Quadrature.CellQuadratureScheme(true, null));

            double min; double max;
            rhoOldTimeStep.GetSpeciesShadowField("A").GetExtremalValues(out min, out max);

            //XDGField rhoOld = (XDGField)ParameterVarFields[paramName];
            //XDGField rho = (XDGField)ParameterVarFields[VariableNames.Rho];
            //rhoOld = rho.CloneAs();
        }

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var density = new (string, DGField)[1];
            var temperature = (XDGField)DomainVarFields[VariableNames.Temperature]; // use same basis as Temperature field =================> is this correct?
            //ConventionalDGField densityField = new SinglePhaseField(temperature.Basis.NonX_Basis, paramName);
            var densityField = new XDGField(temperature.Basis, paramName);

            density[0] = (paramName, densityField);
            return density;
        }

        public override DelParameterFactory Factory => ParameterFactory;
    }

    /// <summary>
    ///
    /// </summary>
    internal class Density : ParameterS {
        public override IList<string> ParameterNames => new string[] { paramName };

        public override DelParameterFactory Factory => DensityFactory;
        private string paramName;

        private MaterialLaw m_EoS_A;
        private MaterialLaw m_EoS_B;
        private int m_NumberOfChemicalComponents;

        public Density(MaterialLaw EoS_A, MaterialLaw EoS_B, int NumberOfChemicalComponents) {
            paramName = VariableNames.Rho;
            m_EoS_A = EoS_A;
            m_EoS_B = EoS_B;
            m_NumberOfChemicalComponents = NumberOfChemicalComponents;
        }

        private (string, DGField)[] DensityFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var density = new (string, DGField)[1];
            var temperature = (XDGField)DomainVarFields[VariableNames.Temperature]; // use same basis as Temperature field
            var densityField = new XDGField(temperature.Basis, VariableNames.Rho);
            density[0] = (paramName, densityField);
            return density;
        }

        public override DelPartialParameterUpdate Update => DensityUpdate;

        private void DensityUpdate(double t, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            var rho = (XDGField)ParameterVarFields[VariableNames.Rho];
            // obtain density arguments
            var Temperature = (XDGField)DomainVarFields[VariableNames.Temperature];
            var massfractions = new XDGField[m_NumberOfChemicalComponents];
            for(int i = 0; i < m_NumberOfChemicalComponents; i++) {
                massfractions[i] = (XDGField)DomainVarFields[VariableNames.MassFraction_n(i)];
            }
            rho.Clear();

            string[] species = new string[] { "A", "B" };

            foreach(var sp in species) {
                MaterialLaw EoS = sp == "A" ? m_EoS_A : m_EoS_B;               
                   

                rho.GetSpeciesShadowField(sp).ProjectField(1.0,
                   delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                       int K = result.GetLength(1);
                       MultidimensionalArray tempT = MultidimensionalArray.Create(Len, K);
                       var tempMFs = new MultidimensionalArray[m_NumberOfChemicalComponents];

                       Temperature.GetSpeciesShadowField(sp).Evaluate(j0, Len, NS, tempT);
                       for(int i = 0; i < m_NumberOfChemicalComponents; i++) {
                           tempMFs[i] = MultidimensionalArray.Create(Len, K);
                           massfractions[i].GetSpeciesShadowField(sp).Evaluate(j0, Len, NS, tempMFs[i]);
                       }
                       for(int j = 0; j < Len; j++) {
                           for(int k = 0; k < K; k++) {
                               double[] densityArguments = new double[1 + m_NumberOfChemicalComponents];
                               densityArguments[0] = tempT[j, k];

                               for(int i = 0; i < m_NumberOfChemicalComponents; i++) {
                                   densityArguments[1 + i] = tempMFs[i][j, k];
                               }

                               result[j, k] = EoS.GetDensity(densityArguments);
                           }
                       }
                   }, new Foundation.Quadrature.CellQuadratureScheme(true, null));
            }
        }
    }

    /// <summary>
    ///
    /// </summary>
    internal class ReactionRate : ParameterS {
        public override IList<string> ParameterNames => new string[] { paramName };

        public override DelParameterFactory Factory => ReactionRateFactory;
        private string paramName;

        private MaterialLaw_MultipleSpecies m_EoS;

        public ReactionRate(MaterialLaw_MultipleSpecies EoS) {
            paramName = "kReact";
            m_EoS = EoS;
        }

        private (string, DGField)[] ReactionRateFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var kReact = new (string, DGField)[1];
            var temperature = (XDGField)DomainVarFields[VariableNames.Temperature]; // use same basis as Temperature field

            var kReactField = new XDGField(temperature.Basis, paramName);
            kReact[0] = (paramName, kReactField);
            return kReact;
        }

        public override DelPartialParameterUpdate Update => kreactUpdate;

        private void kreactUpdate(double t, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            var kReact = (XDGField)ParameterVarFields["kReact"];
            var Temperature = (XDGField)DomainVarFields[VariableNames.Temperature];
            var MassFraction0 = (XDGField)DomainVarFields[VariableNames.MassFraction0];
            var MassFraction1 = (XDGField)DomainVarFields[VariableNames.MassFraction1];
            var MassFraction2 = (XDGField)DomainVarFields[VariableNames.MassFraction2];
            var MassFraction3 = (XDGField)DomainVarFields[VariableNames.MassFraction3];
            string sp = "A";//....
            double Ta = 15900 / 300;
            // obtain density arguments
            kReact.GetSpeciesShadowField(sp).ProjectField(1.0,
                                delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                                    int K = result.GetLength(1);
                                    MultidimensionalArray Y0 = MultidimensionalArray.Create(Len, K);
                                    MultidimensionalArray Y1 = MultidimensionalArray.Create(Len, K);
                                    MultidimensionalArray Y2 = MultidimensionalArray.Create(Len, K);
                                    MultidimensionalArray Y3 = MultidimensionalArray.Create(Len, K);

                                    MultidimensionalArray T = MultidimensionalArray.Create(Len, K);
                                    MassFraction0.GetSpeciesShadowField(sp).Evaluate(j0, Len, NS, Y0);
                                    MassFraction1.GetSpeciesShadowField(sp).Evaluate(j0, Len, NS, Y1);
                                    MassFraction2.GetSpeciesShadowField(sp).Evaluate(j0, Len, NS, Y2);
                                    MassFraction3.GetSpeciesShadowField(sp).Evaluate(j0, Len, NS, Y3);
                                    Temperature.GetSpeciesShadowField(sp).Evaluate(j0, Len, NS, T);

                                    for(int j = 0; j < Len; j++) {
                                        for(int k = 0; k < K; k++) {
                                            double Temp = T[j, k];
                                            double _Y0 = Y0[j, k];
                                            double _Y1 = Y1[j, k];
                                            double _Y2 = Y2[j, k];
                                            double _Y3 = Y3[j, k];
                                            double rho = m_EoS.GetDensity(new double[] { Temp, _Y0, _Y1, _Y2, _Y3 });
                                            double TA = m_EoS.m_ChemModel.getTa(_Y0, _Y1);
                                            result[j, k] = Math.Exp(-TA / Temp * 300) * (rho * _Y0 / m_EoS.MolarMasses[0]) * (rho * _Y1 / m_EoS.MolarMasses[1]);
                                        }
                                    }
                                }, new Foundation.Quadrature.CellQuadratureScheme(true, null));

            double min; double max;
            kReact.GetExtremalValues(out min, out max);
        }
    }

    /// <summary>
    ///
    /// </summary>
    internal class Viscosity : ParameterS {
        public override IList<string> ParameterNames => new string[] { paramName };

        public override DelParameterFactory Factory => ViscosityFactory;
        private string paramName;

        private MaterialLaw m_EoS_A;
        private MaterialLaw m_EoS_B;
        private int m_NumberOfChemicalComponents;

        public Viscosity(MaterialLaw EoS_A, MaterialLaw EoS_B) {
            paramName = VariableNames.Mu;
            m_EoS_A = EoS_A;
            m_EoS_B = EoS_B;
        }

        private (string, DGField)[] ViscosityFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var viscosity = new (string, DGField)[1];
            var temperature = (XDGField)DomainVarFields[VariableNames.Temperature]; // use same basis as Temperature field
            var ViscosityField = new XDGField(temperature.Basis, VariableNames.Mu);
            viscosity[0] = (paramName, ViscosityField);
            return viscosity;
        }

        public override DelPartialParameterUpdate Update => ViscosityUpdate;

        private void ViscosityUpdate(double t, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            var mu = (XDGField)ParameterVarFields[VariableNames.Mu];

            var Temperature = (XDGField)DomainVarFields[VariableNames.Temperature];
            mu.Clear();

            string[] species = new string[] { "A", "B" };

            foreach(var sp in species) {
                MaterialLaw EoS = sp == "A" ? m_EoS_A : m_EoS_B;

                mu.GetSpeciesShadowField(sp).ProjectField(1.0,
                   delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                       int K = result.GetLength(1);
                       MultidimensionalArray tempT = MultidimensionalArray.Create(Len, K);
                       Temperature.GetSpeciesShadowField(sp).Evaluate(j0, Len, NS, tempT);

                       for(int j = 0; j < Len; j++) {
                           for(int k = 0; k < K; k++) {
                               double Temperature = tempT[j, k];
                               result[j, k] = EoS.GetViscosity(Temperature);
                           }
                       }
                   }, new Foundation.Quadrature.CellQuadratureScheme(true, null));
            }
        }
    }

    /// <summary>
    ///
    /// </summary>
    internal class HeatCapacity : ParameterS {
        public override IList<string> ParameterNames => new string[] { paramName };

        public override DelParameterFactory Factory => HeatCapacityFactory;
        private string paramName;

        private MaterialLaw m_EoS_A;
        private MaterialLaw m_EoS_B;
        public HeatCapacity(MaterialLaw EoS_A, MaterialLaw EoS_B ) {
            paramName = VariableNames.cp;
            m_EoS_A = EoS_A;
            m_EoS_B = EoS_B;
        }

        private (string, DGField)[] HeatCapacityFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var cp = new (string, DGField)[1];
            var temperature = (XDGField)DomainVarFields[VariableNames.Temperature]; // use same basis as Temperature field
            var cpField = new XDGField(temperature.Basis, VariableNames.cp);
            cp[0] = (paramName, cpField);
            return cp;
        }

        public override DelPartialParameterUpdate Update => HeatCapacityUpdate;

        private void HeatCapacityUpdate(double t, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            var cp = (XDGField)ParameterVarFields[VariableNames.cp];

            var Temperature = (XDGField)DomainVarFields[VariableNames.Temperature];
            cp.Clear();

            string[] species = new string[] { "A", "B" };

            foreach(var sp in species) {
                MaterialLaw EoS = sp == "A" ? m_EoS_A : m_EoS_B;

                cp.GetSpeciesShadowField(sp).ProjectField(1.0,
                   delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                       int K = result.GetLength(1);
                       MultidimensionalArray tempT = MultidimensionalArray.Create(Len, K);
                       Temperature.GetSpeciesShadowField(sp).Evaluate(j0, Len, NS, tempT);

                       for(int j = 0; j < Len; j++) {
                           for(int k = 0; k < K; k++) {
                               double Temperature = tempT[j, k];
                               result[j, k] = EoS.GetMixtureHeatCapacity(new double[] { Temperature, 0.0, 0.0, 0.0, 0.0/*, 1.0 */});
                           }
                       }
                   }, new Foundation.Quadrature.CellQuadratureScheme(true, null));
            }
        }
    }
}