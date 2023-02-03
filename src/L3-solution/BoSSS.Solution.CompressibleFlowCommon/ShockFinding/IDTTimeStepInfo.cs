using BoSSS.Foundation;
using BoSSS.Foundation.IO;
using System;
using System.Collections.Generic;
using System.Runtime.Serialization;

namespace BoSSS.Solution.CompressibleFlowCommon.IDTUtils {


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
            IEnumerable<DGField> fields, double gamma, double alpha,
            double mu, double residual, double enrichedresidual,
            List<double> AlphaHistory, List<double> GammaHistory, List<double> MuHistory,
            List<double> ResHistory, List<double> EnResHistory, List<double> TimeStepNumbers, List<double> levelSetParams)
            : base(physTime, session, TimestepNo, fields) //
        {
            Gamma = gamma;
            Alpha = alpha;
            Mu = mu;
            Residual = residual;
            EnrichedResidual = enrichedresidual;
            this.AlphaHistory = AlphaHistory.ToArray();
            this.GammaHistory = GammaHistory;
            this.MuHistory = MuHistory;
            this.ResHistory = ResHistory;
            this.EnResHistory = EnResHistory;
            this.TimeStepNumbers = TimeStepNumbers;
            LevelSetParams = levelSetParams;
        }

        [DataMember]
        public double Gamma;
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
        public List<double> MuHistory;
        [DataMember]
        public List<double> ResHistory;
        [DataMember]
        public List<double> EnResHistory;
        [DataMember]
        public List<double> TimeStepNumbers;
        [DataMember]
        public List<double> LevelSetParams;
    }
}

