using BoSSS.Foundation;
using BoSSS.Solution.XNSECommon;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Web;

namespace BoSSS.Application.XNSE_Solver.Logging {

    public abstract class XNSELogValue {

        abstract public string[] LogValues { get; }

        /// <summary>
        /// Main implementation point for the post-processing routine after each time-step
        /// </summary>
        abstract public double[] ComputeLogValues<T>(XNSEinSituPostProcessingModule<T> PPmodule) where T : XNSE_Control, new();
    }


    #region energy values


    public enum EnergyLogValues { 
    
        Kinetic = 0,

        KineticPhaseA = 1,

        KineticPhaseB = 2,

        Surface = 3,

        SurfaceDivergence = 4,

        KineticDissipation = 5,

        KineticDissipationPhaseA = 6,

        KineticDissipationPhaseB = 7,

    }



    public class XNSELogValue_Energy : XNSELogValue {

        public XNSELogValue_Energy() { }

        List<EnergyLogValues> EnergyLogValues;

        public XNSELogValue_Energy(List<EnergyLogValues> LogValues) {
            EnergyLogValues = LogValues;
        }


        public void AddEnergyLogValue(EnergyLogValues newLogValue) {
            this.EnergyLogValues.Add(newLogValue);
        }


        public override string[] LogValues => GetLogValueNames();


        protected string[] GetLogValueNames() {
            string[] ret = new string[EnergyLogValues.Count()];
            for(int i = 0; i < EnergyLogValues.Count(); i++) {
                ret[i] = LogValues[i].ToString();
            }
            return ret;
        }


        public override double[] ComputeLogValues<T>(XNSEinSituPostProcessingModule<T> PPmodule) {
            double[] val = new double[EnergyLogValues.Count()];
            for(int i = 0; i < EnergyLogValues.Count(); i++) {
                //val[i] = ComputeEnergyLogValue(PPmodule, EnergyLogValues[i]);
            }

            return val;
        }

        //protected double ComputeEnergyLogValue<T>(XNSEinSituPostProcessingModule<T> PPmodule, EnergyLogValues logValue) where T : XNSE_Control, new() {

        //    var PhysParam = PPmodule.Control.PhysicalParameters;
        //    var CurrVel = PPmodule.CurrentVel;

        //    int quadOrder = CurrVel[0].Basis.Degree * CurrVel[0].Basis.Degree;

        //    double val = 0;

        //    switch(logValue) {
        //        case Logging.EnergyLogValues.Kinetic: {
        //            double[] rhoS = new double[] { PhysParam.rho_A, PhysParam.rho_B };
        //            val = EnergyUtils.GetKineticEnergy(PPmodule.LsTrk, CurrVel, rhoS, quadOrder);
        //            break;
        //        }
        //        case Logging.EnergyLogValues.KineticPhaseA: {
        //            double[] rhoS = new double[] { PhysParam.rho_A, 0.0 };
        //            val = EnergyUtils.GetKineticEnergy(PPmodule.LsTrk, CurrVel, rhoS, quadOrder);
        //            break;
        //        }
        //        case Logging.EnergyLogValues.KineticPhaseB: {
        //            double[] rhoS = new double[] { 0.0, PhysParam.rho_B };
        //            val = EnergyUtils.GetKineticEnergy(PPmodule.LsTrk, CurrVel, rhoS, quadOrder);
        //            break;
        //        }
        //        case Logging.EnergyLogValues.Surface: {
        //            double[] rhoS = new double[] { 0.0, PhysParam.rho_B };
        //            val = EnergyUtils.GetSurfaceEnergy(PPmodule.LsTrk, PhysParam.Sigma, quadOrder);
        //            break;
        //        }
        //        case Logging.EnergyLogValues.SurfaceDivergence: {
        //            ConventionalDGField[] meanVelocity = XNSEUtils.GetMeanVelocity(CurrVel, PPmodule.LsTrk, PhysParam.rho_A, PhysParam.rho_B);
        //            val = EnergyUtils.GetSurfaceChangerate(PPmodule.LsTrk, meanVelocity, quadOrder);
        //            break;
        //        }
        //        case Logging.EnergyLogValues.KineticDissipation: {
        //            double[] muS = new double[] { PhysParam.mu_A, PhysParam.mu_B };
        //            val = EnergyUtils.GetKineticDissipation(PPmodule.LsTrk, CurrVel, muS, quadOrder);
        //            break;
        //        }
        //        case Logging.EnergyLogValues.KineticDissipationPhaseA: {
        //            double[] muS = new double[] { PhysParam.mu_A, 0.0 };
        //            val = EnergyUtils.GetKineticDissipation(PPmodule.LsTrk, CurrVel, muS, quadOrder);
        //            break;
        //        }
        //        case Logging.EnergyLogValues.KineticDissipationPhaseB: {
        //            double[] muS = new double[] { 0.0, PhysParam.mu_B };
        //            val = EnergyUtils.GetKineticDissipation(PPmodule.LsTrk, CurrVel, muS, quadOrder);
        //            break;
        //        }
        //    }

        //    return val;
        //}

    }

    #endregion
}
