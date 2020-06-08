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
    /// This class allows to use DG fields defined on 
    ///  different mesh to be used 
    /// for boundary or initial values.
    /// </summary>
    [DataContract]
    [Serializable]
    public class ForeignGridValue : BoSSS.Solution.Control.IBoundaryAndInitialData {

        /// <summary>
        /// Private ctor for 
        /// </summary>
        private ForeignGridValue() {

        }


        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="tst"></param>
        /// <param name="__FieldName">
        /// <see cref="FieldName"/>
        /// </param>
        public ForeignGridValue(ITimestepInfo tst, string __FieldName) {
            FieldName = __FieldName;
            TimestepID = tst.ID;
            SessionID = tst.Session.ID;

            var bla = tst.Database.AlternateDbPaths.ToList();
            bla.Insert(0, (tst.Database.Path, null));

            this.DbPaths = bla.ToArray();
        }

        /// <summary>
        /// Paths to search for the Database
        /// </summary>
        [DataMember]
        (string DbPath, string MachineFilter)[] DbPaths;



        /// <summary>
        /// <see cref="IDatabaseEntityInfo{T}.ID"/>
        /// </summary>
        [DataMember]
        public Guid TimestepID {
            get;
            private set;
        }


        /// <summary>
        /// <see cref="IDatabaseEntityInfo{T}.ID"/>
        /// </summary>
        [DataMember]
        public Guid SessionID {
            get;
            private set;
        }

        /// <summary>
        /// Field name within the provided timestep 
        /// (<see cref="ITimestepInfo.Fields"/>, <see cref="DGField.Identification"/>).
        /// </summary>
        [DataMember]
        public string FieldName {
            get;
            private set;
        }

        [NonSerialized]
        IDatabaseInfo m_Database;

        [NonSerialized]
        ISessionInfo m_Session;

        [NonSerialized]
        ITimestepInfo m_Timestep;

        [NonSerialized]
        DGField m_dGField;

        [NonSerialized]
        FieldEvaluation m_eval;

        void Init() {
            if(m_Database == null) {
                m_Database = DatabaseInfo.Open(this.DbPaths);
            }

            if(m_Session == null) {
                m_Session = m_Database.Sessions.Single(si => si.ID.Equals(SessionID));
            }

            if(m_Timestep == null) {
                m_Timestep = m_Session.Timesteps.Single(ti => ti.ID.Equals(TimestepID));
            }

            if(m_dGField == null) {
                m_dGField = m_Timestep.Fields.Single(f => f.Identification.Equals(FieldName));
            }

            if(m_eval == null) {
                m_eval = new FieldEvaluation(GridHelper.ExtractGridData(m_dGField.GridDat));
            }
        }


        /// <summary>
        /// scalar evaluation
        /// </summary>
        public double Evaluate(double[] X, double t) {
            Init();
            return m_eval.Evaluate(X, m_dGField);
        }

        /// <summary>
        /// vectorized evaluation
        /// </summary>
        public void Evaluate(MultidimensionalArray input, double time, MultidimensionalArray output) {
            Init();

            int L = input.GetLength(0);
            m_eval.Evaluate(1.0, new DGField[] { m_dGField }, input, 0.0, output.ResizeShallow(L, 1));
        }


        /// <summary>
        /// %
        /// </summary>
        public override bool Equals(object obj) {
            var odha = obj as ForeignGridValue;
            if(odha == null)
                return false;

            if(!odha.SessionID.Equals(this.SessionID))
                return false;

            if(!odha.TimestepID.Equals(this.TimestepID))
                return false;

            return true;
        }

        /// <summary>
        /// %
        /// </summary>
        public override int GetHashCode() {
            return SessionID.GetHashCode();
        }
    }
}

