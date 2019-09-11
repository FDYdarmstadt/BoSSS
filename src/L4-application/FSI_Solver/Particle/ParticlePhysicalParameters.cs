using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.FSI_Solver {
    public class ParticlePhysicalParameters {
        /// <summary>
        /// Constructor for physical parameters of the particle
        /// </summary>
        /// <param name="particleDensity"></param>
        /// <param name="particleActiveStress"></param>
        public ParticlePhysicalParameters(double particleDensity, double particleActiveStress) {
            density = particleDensity;
            activeStress = particleActiveStress;
        }

        /// <summary>
        /// Density of the particle
        /// </summary>
        public double density = 1;

        /// <summary>
        /// Active stress on the particle surface. Positive value indicates a pusher, negative value a puller.
        /// </summary>
        public double activeStress = 0;
    }
}
