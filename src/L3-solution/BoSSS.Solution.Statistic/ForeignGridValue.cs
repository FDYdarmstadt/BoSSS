using BoSSS.Foundation;
using BoSSS.Foundation.IO;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.Statistic {


    /// <summary>
    /// 
    /// </summary>
    [DataContract]
    [Serializable]
    public class ForeignGridValue : BoSSS.Solution.Control.IBoundaryAndInitialData {

        public ForeignGridValue(ITimestepInfo tst, string VarName) : this(tst.ID, VarName) {
            throw new NotImplementedException("todo");
        }


        public ForeignGridValue(Guid TimestepID, string VarName) {
            this.TimestepID = TimestepID;
            throw new NotImplementedException("todo");
        }

        [DataMember]
        public Guid TimestepID {
            get;
            private set;
        }

        [DataMember]
        public string VariName {
            get;
            private set;
        }

        IDatabaseInfo m_Database;

        [NonSerialized]
        ITimestepInfo m_Timestep;

        [NonSerialized]
        DGField dGField;

        void Init() {

            //m_Database.Controller.Get


            dGField = m_Timestep.Fields.Single(f => f.Identification.Equals(VariName));
        }



        public double Evaluate(double[] X, double t) {
            throw new NotImplementedException();
        }

        public void Evaluate(MultidimensionalArray input, double time, MultidimensionalArray output) {
            throw new NotImplementedException();
        }
    }
}
