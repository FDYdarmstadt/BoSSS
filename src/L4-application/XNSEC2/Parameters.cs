using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Foundation.XDG.OperatorFactory;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.Linq;

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

        public (string, DGField)[] ThermodynamicPressureFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var ThermodynamicPressure0 = new (string, DGField)[1];
            string ThermodynamicPressureName = BoSSS.Solution.NSECommon.VariableNames.ThermodynamicPressure;
            //var grd = DomainVarFields[VariableNames.Temperature].GridDat;
            //     var b = new Basis(grd, 0); // polynomial degree zero
            //var ThermodynamicPressure = new SinglePhaseField(b, ThermodynamicPressureName); // Should i define it as singlephasefield or XDGfield?
            var ThermodynamicPressure = new XDGField((XDGBasis)DomainVarFields[VariableNames.VelocityX].Basis, ThermodynamicPressureName);
            ThermodynamicPressure0[0] = (ThermodynamicPressureName, ThermodynamicPressure);

            if (!m_EoS.IsInitialized) {
                XDGField Temperature = (XDGField)DomainVarFields[VariableNames.Temperature];
                double init_p0 = m_EoS.GetMassDeterminedThermodynamicPressure(m_initialMass, Temperature);

                m_EoS.Initialize(init_p0);
            }

            return ThermodynamicPressure0;
        }

        public override DelPartialParameterUpdate Update => ThermodynamicPressureUpdate;

        private void ThermodynamicPressureUpdate(double t, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            XDGField ThermodynamicPressure = (XDGField)ParameterVarFields[VariableNames.ThermodynamicPressure];
            ThermodynamicPressure.Clear();

            double p0;
            switch (m_mode) {
                case ThermodynamicPressureMode.Constant:
                    p0 = 1.0;
                    m_EoS.ThermodynamicPressure = p0;
                    break;

                case ThermodynamicPressureMode.MassDetermined:
                    XDGField Temperature = (XDGField)DomainVarFields[VariableNames.Temperature];
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
            for (int i = 0; i < D; i++) {
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
                      for (int j = 0; j < Len; j++) {
                          for (int k = 0; k < K; k++) {
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
            foreach (int cell in cellMask.ItemEnum) {
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

    /// <summary>
    /// rho_t-2
    /// </summary>
    internal class Density_t00 : ParameterS/*, ILevelSetParameter*/ {
        private string paramName;
        private MaterialLaw m_EoS;

        public Density_t00(int NumberOfChemicalSpecies, MaterialLaw EoS) {
            paramName = "Density_t00";
            m_EoS = EoS;
        }

        public override IList<string> ParameterNames => new string[] { paramName };

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var density = new (string, DGField)[1];
            XDGField field;
            if (m_EoS is MaterialLawMixtureFractionNew) {
                field = (XDGField)DomainVarFields[VariableNames.MixtureFraction]; // use same basis as Temperature field =================> is this correct?
            } else {
                field = (XDGField)DomainVarFields[VariableNames.Temperature]; // use same basis as Temperature field =================> is this correct?
            }

            var densityField = new XDGField(field.Basis, paramName);

            density[0] = (paramName, densityField);
            return density;
        }

        public override DelParameterFactory Factory => ParameterFactory;
    }

    /// <summary>
    /// rho_t-1
    /// </summary>
    internal class Density_t0 : ParameterS/*, ILevelSetParameter*/ {
        private string paramName;
        private MaterialLaw m_EoS;

        public Density_t0(int NumberOfChemicalSpecies, MaterialLaw EoS) {
            paramName = "Density_t0";
            m_EoS = EoS;
        }

        public override IList<string> ParameterNames => new string[] { paramName };

        public (string ParameterName, DGField ParamField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var density = new (string, DGField)[1];
            XDGField field;
            if (m_EoS is MaterialLawMixtureFractionNew) {
                field = (XDGField)DomainVarFields[VariableNames.MixtureFraction]; // use same basis as Temperature field =================> is this correct?
            } else {
                field = (XDGField)DomainVarFields[VariableNames.Temperature]; // use same basis as Temperature field =================> is this correct?
            }

            var densityField = new XDGField(field.Basis, paramName);

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
            for (int i = 0; i < m_NumberOfChemicalComponents; i++) {
                massfractions[i] = (XDGField)DomainVarFields[VariableNames.MassFraction_n(i)];
            }
            rho.Clear();

            string[] species = new string[] { "A", "B" };

            foreach (var sp in species) {
                MaterialLaw EoS = sp == "A" ? m_EoS_A : m_EoS_B;

                rho.GetSpeciesShadowField(sp).ProjectField(1.0,
                   delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                       int K = result.GetLength(1);
                       MultidimensionalArray tempT = MultidimensionalArray.Create(Len, K);
                       var tempMFs = new MultidimensionalArray[m_NumberOfChemicalComponents];

                       Temperature.GetSpeciesShadowField(sp).Evaluate(j0, Len, NS, tempT);
                       for (int i = 0; i < m_NumberOfChemicalComponents; i++) {
                           tempMFs[i] = MultidimensionalArray.Create(Len, K);
                           massfractions[i].GetSpeciesShadowField(sp).Evaluate(j0, Len, NS, tempMFs[i]);
                       }
                       for (int j = 0; j < Len; j++) {
                           for (int k = 0; k < K; k++) {
                               double[] densityArguments = new double[1 + m_NumberOfChemicalComponents];
                               densityArguments[0] = tempT[j, k] == 0 ? 1.0 : tempT[j, k];


                               for (int i = 0; i < m_NumberOfChemicalComponents; i++) {
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
    internal class DensityMF : ParameterS {
        public override IList<string> ParameterNames => new string[] { paramName };

        public override DelParameterFactory Factory => DensityFactory;
        private string paramName;

        private MaterialLaw m_EoS_A;
        private MaterialLaw m_EoS_B;
        private int m_NumberOfChemicalComponents;

        public DensityMF(MaterialLaw EoS_A, MaterialLaw EoS_B, int NumberOfChemicalComponents) {
            paramName = VariableNames.Rho;
            m_EoS_A = EoS_A;
            m_EoS_B = EoS_B;
            m_NumberOfChemicalComponents = NumberOfChemicalComponents;
        }

        private (string, DGField)[] DensityFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var density = new (string, DGField)[1];
            var MixtureFraction = (XDGField)DomainVarFields[VariableNames.MixtureFraction];
            var densityField = new XDGField(MixtureFraction.Basis, VariableNames.Rho);
            density[0] = (paramName, densityField);
            return density;
        }

        public override DelPartialParameterUpdate Update => DensityUpdate;

        private void DensityUpdate(double t, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            var rho = (XDGField)ParameterVarFields[VariableNames.Rho];
            // obtain density arguments
            var MixtureFraction = (XDGField)DomainVarFields[VariableNames.MixtureFraction];

            rho.Clear();

            string[] species = new string[] { "A", "B" };

            foreach (var sp in species) {
                MaterialLaw EoS = sp == "A" ? m_EoS_A : m_EoS_B;

                rho.GetSpeciesShadowField(sp).ProjectField(1.0,
                   delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                       int K = result.GetLength(1);
                       MultidimensionalArray tempZ = MultidimensionalArray.Create(Len, K);

                       MixtureFraction.GetSpeciesShadowField(sp).Evaluate(j0, Len, NS, tempZ);

                       for (int j = 0; j < Len; j++) {
                           for (int k = 0; k < K; k++) {
                               double[] densityArguments = new double[1];
                               densityArguments[0] = tempZ[j, k];

                               result[j, k] = EoS.getDensityFromZ(densityArguments[0]);
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
        private double Ta;
        private double Da;
        private int noOfChemComponents;

        public ReactionRate(MaterialLaw_MultipleSpecies EoS, double[] ReactionRateConstants, int _noOfChemComponents) {
            paramName = "kReact";
            m_EoS = EoS;
            Da = ReactionRateConstants[0];
            Ta = ReactionRateConstants[1];
            noOfChemComponents = _noOfChemComponents;
        }

        private (string, DGField)[] ReactionRateFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var kReact = new (string, DGField)[1];
            var temperature = (XDGField)DomainVarFields[VariableNames.Temperature]; // use same basis as Temperature field
            var b = new Basis(temperature.GridDat, temperature.Basis.Degree * 4); // polynomial degree zero
            var kReactField = new XDGField(temperature.Basis, paramName);
            kReact[0] = (paramName, kReactField);
            return kReact;
        }

        public override DelPartialParameterUpdate Update => kreactUpdate;

        private void kreactUpdate(double t, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            var kReact = (XDGField)ParameterVarFields["kReact"];
            XDGField Temperature = (XDGField)DomainVarFields[VariableNames.Temperature];

            XDGField[] MassFractions = new XDGField[noOfChemComponents];

            for (int i = 0; i < noOfChemComponents; i++) {
                MassFractions[i] = (XDGField)DomainVarFields[VariableNames.MassFraction_n(i)];
            }

            string[] species = new string[] { "A", "B" };

            foreach (string sp in species) {
                kReact.GetSpeciesShadowField(sp).ProjectField(1.0,
                                    delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                                        int K = result.GetLength(1);
                                        MultidimensionalArray T = MultidimensionalArray.Create(Len, K);
                                        MultidimensionalArray[] Ys = new MultidimensionalArray[noOfChemComponents];

                                        for (int i = 0; i < noOfChemComponents; i++) {
                                            Ys[i] = MultidimensionalArray.Create(Len, K);
                                            MassFractions[i].GetSpeciesShadowField(sp).Evaluate(j0, Len, NS, Ys[i]);
                                        }

                                        Temperature.GetSpeciesShadowField(sp).Evaluate(j0, Len, NS, T);

                                        for (int j = 0; j < Len; j++) {
                                            for (int k = 0; k < K; k++) {
                                                double[] _Ys = new double[noOfChemComponents];
                                                for (int i = 0; i < noOfChemComponents; i++) {
                                                    _Ys[i] = (Ys[i])[j, k];
                                                }
                                                double Temp = T[j, k];
                                                double[] phi = (new double[] { Temp }).Concat(_Ys).ToArray();
                                                double rho = m_EoS.GetDensity(phi);
                                            //double TA = m_EoS.m_ChemModel.getTa(_Y0, _Y1);
                                            result[j, k] = Da * Math.Exp(-Ta / Temp) * (rho * _Ys[0] / m_EoS.MolarMasses[0]) * (rho * _Ys[1] / m_EoS.MolarMasses[1]);
                                            //result[j, k] = (rho * _Y0 / m_EoS.MolarMasses[0]) * (rho * _Y1 / m_EoS.MolarMasses[1]);
                                        }
                                        }
                                    }, new Foundation.Quadrature.CellQuadratureScheme(true, null));
            }
            double min; double max;
            kReact.GetExtremalValues(out min, out max);
        }
    }

    internal class Temperature0 : ParameterS {
        public override IList<string> ParameterNames => new string[] { BoSSS.Solution.NSECommon.VariableNames.Temperature0 };

        public override DelParameterFactory Factory => Temperature0Factory;

        public Temperature0() {
        }

        private (string, DGField)[] Temperature0Factory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var temperature0 = new (string, DGField)[1];

            string temperaturename = BoSSS.Solution.NSECommon.VariableNames.Temperature;
            DGField temperature = DomainVarFields[temperaturename];
            string paramName = BoSSS.Solution.NSECommon.VariableNames.Temperature0;
            temperature0[0] = (paramName, temperature);

            return temperature0;
        }
    }

    internal class MassFraction0_0 : ParameterS {
        public override IList<string> ParameterNames => new string[] { BoSSS.Solution.NSECommon.VariableNames.MassFraction0_0 };

        public override DelParameterFactory Factory => MassFraction0_0Factory;

        public MassFraction0_0() {
        }

        private (string, DGField)[] MassFraction0_0Factory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            var massFraction0_0 = new (string, DGField)[1];

            string mfname = BoSSS.Solution.NSECommon.VariableNames.MassFraction0;
            DGField mf0 = DomainVarFields[mfname];
            string paramName = BoSSS.Solution.NSECommon.VariableNames.MassFraction0_0;
            massFraction0_0[0] = (paramName, mf0);

            return massFraction0_0;
        }
    }

    /// <summary>
    /// Cell-wise mean value, required for the for the
    /// localized Lax-Friedrichs flux <see cref="XNSECommon.Operator.Convection.ConvectionInBulk_LLF"/>,
    /// to have a constant Eigenvalue (aka. flow direction) along an edge.
    /// </summary>
    public class Temperature0Mean : ParameterS {
        protected int cutCellQuadOrder;

        protected LevelSetTracker LsTrk;

        public Temperature0Mean(LevelSetTracker LsTrk, int cutCellQuadOrder) {
            this.cutCellQuadOrder = cutCellQuadOrder;
            this.LsTrk = LsTrk;
        }

        public override DelParameterFactory Factory => ParameterFactory;

        public override IList<string> ParameterNames => new string[] { BoSSS.Solution.NSECommon.VariableNames.Temperature0Mean };

        public (string, DGField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            this.LsTrk = ((XDGBasis)DomainVarFields.First().Value.Basis).Tracker;
            var Temperature0Mean = new (string, DGField)[0];

            XDGBasis T0meanBasis = new XDGBasis(LsTrk, 0);
            string paramName = BoSSS.Solution.NSECommon.VariableNames.Temperature0Mean;
            XDGField U0mean = new XDGField(T0meanBasis, paramName);
            Temperature0Mean[0] = (paramName, U0mean);

            return Temperature0Mean;
        }

        protected IList<string> SpeciesNames;

        protected LevelSetTracker.LevelSetRegions regions;

        protected IDictionary<string, SpeciesId> speciesMap;

        protected XQuadSchemeHelper schemeHelper;

        protected double minvol;

        public override DelPartialParameterUpdate Update {
            get {
                return T0MeanUpdate;
            }
        }

        protected virtual void T0MeanUpdate(double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            using (new FuncTrace()) {
                foreach (string speciesName in SpeciesNames) {
                    XDGField paramMeanTemperature = (XDGField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.Temperature0Mean];
                    DGField speciesParam = paramMeanTemperature.GetSpeciesShadowField(speciesName);

                    XDGField temperature = (XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.Temperature];
                    DGField speciesTemperature = temperature.GetSpeciesShadowField(speciesName);

                    //Uncut
                    speciesParam.SetMeanValueTo(speciesTemperature);

                    //Cut
                    CellMask cutCells = regions.GetSpeciesMask(speciesName);
                    SpeciesId speciesId = speciesMap[speciesName];
                    CellQuadratureScheme scheme = schemeHelper.GetVolumeQuadScheme(speciesId, IntegrationDomain: cutCells);
                    SetMeanValueToMeanOf(speciesParam, speciesTemperature, minvol, cutCellQuadOrder, scheme);
                }
            }
        }

        protected static void SetMeanValueToMeanOf(DGField target, DGField source, double minvol, int order, CellQuadratureScheme scheme) {
            //Cut
            int D = source.GridDat.SpatialDimension;
            var rule = scheme.Compile(source.GridDat, order);
            CellQuadrature.GetQuadrature(new int[] { 2 },
                source.GridDat,
                rule,
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.Clear();
                    source.Evaluate(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    var Vol = EvalResult.ExtractSubArrayShallow(-1, -1, 1);
                    Vol.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++) {
                        int jCell = i + i0;

                        double Volume = ResultsOfIntegration[i, 1];
                        if (Math.Abs(Volume) < minvol * 1.0e-12) {
                            // keep current value
                            // since the volume of species 'Spc' in cell 'jCell' is 0.0, the value in this cell should have no effect
                        } else {
                            double IntVal = ResultsOfIntegration[i, 0];
                            target.SetMeanValue(jCell, IntVal / Volume);
                        }
                    }
                }).Execute();
        }
    }

    /// <summary>
    /// Cell-wise mean value, required for the for the
    /// localized Lax-Friedrichs flux <see cref="XNSECommon.Operator.Convection.ConvectionInBulk_LLF"/>,
    /// to have a constant Eigenvalue (aka. flow direction) along an edge.
    /// </summary>
    public class MassFraction0Mean : ParameterS {
        protected int cutCellQuadOrder;

        protected LevelSetTracker LsTrk;

        public MassFraction0Mean(LevelSetTracker LsTrk, int cutCellQuadOrder) {
            this.cutCellQuadOrder = cutCellQuadOrder;
            this.LsTrk = LsTrk;
        }

        public override DelParameterFactory Factory => ParameterFactory;

        public override IList<string> ParameterNames => new string[] { BoSSS.Solution.NSECommon.VariableNames.MassFraction0Mean };

        public (string, DGField)[] ParameterFactory(IReadOnlyDictionary<string, DGField> DomainVarFields) {
            this.LsTrk = ((XDGBasis)DomainVarFields.First().Value.Basis).Tracker;
            var MassFraction0Mean = new (string, DGField)[0];

            XDGBasis T0meanBasis = new XDGBasis(LsTrk, 0);
            string paramName = BoSSS.Solution.NSECommon.VariableNames.MassFraction0Mean;
            XDGField U0mean = new XDGField(T0meanBasis, paramName);
            MassFraction0Mean[0] = (paramName, U0mean);

            return MassFraction0Mean;
        }

        protected IList<string> SpeciesNames;

        protected LevelSetTracker.LevelSetRegions regions;

        protected IDictionary<string, SpeciesId> speciesMap;

        protected XQuadSchemeHelper schemeHelper;

        protected double minvol;

        public override DelPartialParameterUpdate Update {
            get {
                return T0MeanUpdate;
            }
        }

        protected virtual void T0MeanUpdate(double time, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            using (new FuncTrace()) {
                foreach (string speciesName in SpeciesNames) {
                    XDGField paramMeanTemperature = (XDGField)ParameterVarFields[BoSSS.Solution.NSECommon.VariableNames.MassFraction0Mean];
                    DGField speciesParam = paramMeanTemperature.GetSpeciesShadowField(speciesName);

                    XDGField temperature = (XDGField)DomainVarFields[BoSSS.Solution.NSECommon.VariableNames.MassFraction0];
                    DGField speciesTemperature = temperature.GetSpeciesShadowField(speciesName);

                    //Uncut
                    speciesParam.SetMeanValueTo(speciesTemperature);

                    //Cut
                    CellMask cutCells = regions.GetSpeciesMask(speciesName);
                    SpeciesId speciesId = speciesMap[speciesName];
                    CellQuadratureScheme scheme = schemeHelper.GetVolumeQuadScheme(speciesId, IntegrationDomain: cutCells);
                    SetMeanValueToMeanOf(speciesParam, speciesTemperature, minvol, cutCellQuadOrder, scheme);
                }
            }
        }

        protected static void SetMeanValueToMeanOf(DGField target, DGField source, double minvol, int order, CellQuadratureScheme scheme) {
            //Cut
            int D = source.GridDat.SpatialDimension;
            var rule = scheme.Compile(source.GridDat, order);
            CellQuadrature.GetQuadrature(new int[] { 2 },
                source.GridDat,
                rule,
                delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                    EvalResult.Clear();
                    source.Evaluate(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    var Vol = EvalResult.ExtractSubArrayShallow(-1, -1, 1);
                    Vol.SetAll(1.0);
                },
                delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    for (int i = 0; i < Length; i++) {
                        int jCell = i + i0;

                        double Volume = ResultsOfIntegration[i, 1];
                        if (Math.Abs(Volume) < minvol * 1.0e-12) {
                            // keep current value
                            // since the volume of species 'Spc' in cell 'jCell' is 0.0, the value in this cell should have no effect
                        } else {
                            double IntVal = ResultsOfIntegration[i, 0];
                            target.SetMeanValue(jCell, IntVal / Volume);
                        }
                    }
                }).Execute();
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

            foreach (var sp in species) {
                MaterialLaw EoS = sp == "A" ? m_EoS_A : m_EoS_B;

                mu.GetSpeciesShadowField(sp).ProjectField(1.0,
                   delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                       int K = result.GetLength(1);
                       MultidimensionalArray tempT = MultidimensionalArray.Create(Len, K);
                       Temperature.GetSpeciesShadowField(sp).Evaluate(j0, Len, NS, tempT);

                       for (int j = 0; j < Len; j++) {
                           for (int k = 0; k < K; k++) {
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

        public HeatCapacity(MaterialLaw EoS_A, MaterialLaw EoS_B) {
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

            foreach (var sp in species) {
                MaterialLaw EoS = sp == "A" ? m_EoS_A : m_EoS_B;

                cp.GetSpeciesShadowField(sp).ProjectField(1.0,
                   delegate (int j0, int Len, NodeSet NS, MultidimensionalArray result) {
                       int K = result.GetLength(1);
                       MultidimensionalArray tempT = MultidimensionalArray.Create(Len, K);
                       Temperature.GetSpeciesShadowField(sp).Evaluate(j0, Len, NS, tempT);

                       for (int j = 0; j < Len; j++) {
                           for (int k = 0; k < K; k++) {
                               double Temperature = tempT[j, k];
                               result[j, k] = EoS.GetMixtureHeatCapacity(new double[] { Temperature, 0.0, 0.0, 0.0, 0.0/*, 1.0 */});
                           }
                       }
                   }, new Foundation.Quadrature.CellQuadratureScheme(true, null));
            }
        }
    }
}