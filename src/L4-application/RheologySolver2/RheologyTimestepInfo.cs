using BoSSS.Foundation;
using BoSSS.Foundation.IO;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Runtime.Serialization;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.Rheology {

    /// <summary>
    /// Time-step data which contains additional particle information
    /// </summary>
    [Serializable]
    [DataContract]
    public class RheologyTimestepInfo : TimestepInfo {

        /// <summary>
        /// empty constructor for serialization
        /// </summary>
        protected RheologyTimestepInfo() : base() { }


        /// <summary>
        /// 
        /// </summary>
        public RheologyTimestepInfo(double physTime, ISessionInfo session, TimestepNumber TimestepNo, IEnumerable<DGField> fields, double __currentWeissenbergNumber)
            : base(physTime, session, TimestepNo, fields) //
        {
            currentWeissenbergNumber = __currentWeissenbergNumber;
        }

        /// <summary>
        /// Weissenberg number at which the timestep was computed (<see cref="Rheology.currentWeissenberg"/>)
        /// </summary>
        [DataMember]
        public double currentWeissenbergNumber;
    }
}
