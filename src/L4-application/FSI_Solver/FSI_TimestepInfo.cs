using BoSSS.Foundation;
using BoSSS.Foundation.IO;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;

namespace BoSSS.Application.FSI_Solver {

    /// <summary>
    /// Time-step data which contains additional particle information
    /// </summary>
    public class FSI_TimestepInfo : TimestepInfo {

        /// <summary>
        /// empty constructor for serialization
        /// </summary>
        protected FSI_TimestepInfo() : base() { }


        /// <summary>
        /// 
        /// </summary>
        public FSI_TimestepInfo(double physTime, ISessionInfo session, TimestepNumber TimestepNo, IEnumerable<DGField> fields, IEnumerable<Particle> _particles)
            : base(physTime, session, TimestepNo, fields) //
        {
            Particles = _particles.ToArray().CloneAs();
        }

        /// <summary>
        /// particle state (<see cref="FSI_SolverMain.m_Particles"/>)
        /// </summary>
        public Particle[] Particles;
    }
}
