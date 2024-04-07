using BoSSS.Foundation;
using BoSSS.Foundation.IO;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Gnuplot;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using System.Security.Cryptography.X509Certificates;

namespace ApplicationWithIDT {

    /// <summary>
    /// Time-step data which contains OptimizationHistory
    /// </summary>
    [Serializable]
    [DataContract]
    public class IDTTimeStepInfo : TimestepInfo {


        /// <summary>
        /// empty constructor for serialization
        /// </summary>
        protected IDTTimeStepInfo() : base() { }

        /// <summary>
        /// 
        /// </summary>
        public IDTTimeStepInfo(double physTime, ISessionInfo session, TimestepNumber TimestepNo,
            IEnumerable<DGField> fields, double gamma,double kappa, double alpha,
            double mu, double residual, double enrichedresidual,
            List<double> AlphaHistory, List<double> GammaHistory, List<double> KappaHistory, List<double> MuHistory,
            List<double> ResHistory, List<double> EnResHistory, List<double> TimeStepNumbers, double[] LevelSetParams, double[] splineLevelSetYPoints =null)
            : base(physTime, session, TimestepNo, fields) //
        {
            Gamma = gamma;
            Kappa=kappa;
            Alpha = alpha;
            Mu = mu;
            Residual = residual;
            EnrichedResidual = enrichedresidual;
            this.AlphaHistory = AlphaHistory.ToArray();
            this.GammaHistory = GammaHistory;
            this.KappaHistory = KappaHistory;
            this.MuHistory = MuHistory;
            this.ResHistory = ResHistory;
            this.EnResHistory = EnResHistory;
            this.TimeStepNumbers = TimeStepNumbers;
            this.LevelSetParams = LevelSetParams;
            SplineLevelSetYPoints = splineLevelSetYPoints;
        }

        [DataMember]
        public double Gamma;
        [DataMember]
        public double Kappa;
        [DataMember]
        public double Alpha;
        [DataMember]
        public double Mu;
        [DataMember]
        public double Residual;
        [DataMember]
        public double EnrichedResidual;
        [DataMember]
        public double[] AlphaHistory;
        [DataMember] 
        public List<double> GammaHistory;
        [DataMember]
        public List<double> KappaHistory;
        [DataMember] 
        public List<double> MuHistory;
        [DataMember] 
        public List<double> ResHistory;
        [DataMember] 
        public List<double> EnResHistory;
        [DataMember]
        public List<double> TimeStepNumbers;
        [DataMember]
        public double[] LevelSetParams;
        [DataMember]
        public double[] SplineLevelSetYPoints;

        /// <summary>
        /// this utility function returns the ConservativeVars dependent on the equation
        /// </summary>
        /// <param name="equation"></param>
        /// <returns></returns>
        public XDGField[] GetConsFields(string equation) {
            if(equation == "iEuler2D") {
                XDGField[] ConsFields = new XDGField[4];
                ConsFields[0] = (XDGField) this.Fields.Single(f => f.Identification.Equals("rho"));
                ConsFields[1] = (XDGField)this.Fields.Single(f => f.Identification.Equals("m0"));
                ConsFields[2] = (XDGField)this.Fields.Single(f => f.Identification.Equals("m1"));
                ConsFields[3] = (XDGField)this.Fields.Single(f => f.Identification.Equals("rhoE"));
                return ConsFields;
            //} else if(equation == "iEuler3D") {
            //    XDGField[] ConsFields = new XDGField[5];
            //    return ConsFields;
            //} else if(equation == "iEuler1D") {
            //    XDGField[] ConsFields = new XDGField[3];
            //    return ConsFields;
            } else {
                XDGField[] ConsFields = new XDGField[1];
                ConsFields[0] = (XDGField)this.Fields.Single(f => f.Identification.Equals("c"));
                return ConsFields; 
            }
            
        }
        /// <summary>
        /// Gives a Plot where the two Residuals are plotted against the iteration numbers in a Y-log plot
        /// </summary>
        /// <returns></returns>
        public Plot2Ddata GetResPlot() {
            var resplot = new Plot2Ddata();
            var Fmt = new PlotFormat();
            Fmt.PointType = PointTypes.OpenCircle;
            //Fmt.PointSize = 1.2;
            //Fmt.LineWidth = 1;
            Fmt.Style = Styles.LinesPoints;
            Fmt.LineColor = LineColors.Blue;
            Fmt.PointType = PointTypes.Diamond;

            var Fmt2 = Fmt.CloneAs();
            Fmt2.PointType = PointTypes.Asterisk;
            Fmt2.LineColor = LineColors.Red;

            resplot.AddDataGroup("||R(z)||_2", this.TimeStepNumbers, this.EnResHistory, Fmt);
            resplot.AddDataGroup("||r(z)||_2", this.TimeStepNumbers, this.ResHistory, Fmt2);

            resplot.Xlabel = "Iteration";
            resplot.LogY = true;
            resplot.ShowXtics = true;

            return resplot;
        }

        public List<DGField> GetShadowFields() {
            var fields = this.Fields.ToList();
            var sFields= new List<DGField>();   
            foreach(DGField field in fields) {
                 if(field is XDGField xdgfield) {
                        var shadowL = xdgfield.GetSpeciesShadowField("L");
                        var shadowR = xdgfield.GetSpeciesShadowField("R");
                        sFields.Add(shadowL);
                        sFields.Add(shadowR);
                 } else {
                        sFields.Add(field);
                 }
             }
            return sFields;
        }
        
    }
    
}

