using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {
    public class Solid {
        public double PoissonsRatio = 0;

        //In GPA
        public double ModulusOfElasticity = 0;

        public double Viscosity = 0;

        //In g/cc
        public double Density = 0;

        public double Lame1 = 0;

        //Schubmodul
        public double Lame2 = 0;

    }

    public class IncompressibleViscoElastic : Solid {
        public IncompressibleViscoElastic(double modulusOfElasticity, double viscosity, double density) {
            PoissonsRatio = 0.5;
            ModulusOfElasticity = modulusOfElasticity;
            Lame1 = 0;
            Lame2 = 0.5 / (1 + PoissonsRatio) * ModulusOfElasticity;
            Viscosity = viscosity;
            Density = density;
        }
    }

    class SoftSiliconeRubber : IncompressibleViscoElastic {
        public SoftSiliconeRubber() : base(0.000005, 1, 300) {
        }
    }

    class MediumSiliconeRubber : IncompressibleViscoElastic {
        public MediumSiliconeRubber() : base(0.9, 1, 0.7 ) {
        }
    }

    class HardSiliconeRubber : IncompressibleViscoElastic {
        public HardSiliconeRubber() : base(19, 1, 3.8) {
        }
    }

    class Beam : IncompressibleViscoElastic {
        public Beam() : base(30, 0.1, 1) {
        }
    }

    class ConvergenceTest : IncompressibleViscoElastic {

        public ConvergenceTest():base(3, 0.1, 1){
        }
    }

    public class AllOne : Solid {
        public AllOne() {
            PoissonsRatio = 1;
            ModulusOfElasticity = 1;
            Density = 1;
            Lame1 = 1;
            Lame2 = 1;
            Viscosity = 1;
        }
    }
}
