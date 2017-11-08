using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Control;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;

namespace BoSSS.Application.RefineAndLoadBal {
    public class RefineAndLoadBalControl : AppControl {
        

       

        /// <summary>
        /// Equation coefficient, species A.
        /// </summary>
        public double alpha_A = 0.1;

        /// <summary>
        /// Equation coefficient, species B.
        /// </summary>
        public double alpha_B = 3.5;

        /// <summary>
        /// Time-dependent level set.
        /// </summary>
        public Func<double[], double, double> LevelSet;


        /// <summary>
        /// Solution in domain A.
        /// </summary>
        public double uEx_A(double[] X, double t) {
            return X[0] + alpha_A * t;
        }

        /// <summary>
        /// Solution in domain B.
        /// </summary>
        public double uEx_B(double[] X, double t) {
            return X[0] + alpha_B * t;
        }
    }
}
