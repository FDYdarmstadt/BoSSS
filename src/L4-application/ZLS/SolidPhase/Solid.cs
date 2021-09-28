using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.SolidPhase {
    public class Solid {
        public double PoissonsRatio { get; protected set; }
        
        //In GPA
        public double ModulusOfElasticity { get; protected set; }

        public double Viscosity { get; protected set; }

        //In g/cc
        public double Density { get; protected set; }

        public double Lame1 { get; protected set; }

        //Schubmodul
        public double Lame2 { get; protected set; }

    }

    class SoftSiliconeRubber : Solid {
        public SoftSiliconeRubber() {
            PoissonsRatio = 0.5;
            ModulusOfElasticity = 0.000005;
            Density = 300;
            Lame1 = ModulusOfElasticity * PoissonsRatio / ((1 + PoissonsRatio) * (1 - PoissonsRatio));
            Lame2 = 0.5 / (1 + PoissonsRatio) * ModulusOfElasticity;
        }
    }

    class MediumSiliconeRubber : Solid {
        public MediumSiliconeRubber() {
            PoissonsRatio = 0.5;
            ModulusOfElasticity = 0.9;
            Density = 0.7;
            Lame1 = ModulusOfElasticity * PoissonsRatio / ((1 + PoissonsRatio) * (1 - PoissonsRatio));
            Lame2 = 0.5 / (1 + PoissonsRatio) * ModulusOfElasticity;
        }
    }

    class HardSiliconeRubber : Solid {
        public HardSiliconeRubber() {
            PoissonsRatio = 0.5;
            ModulusOfElasticity = 19;
            Density = 3.8;
            Lame1 = ModulusOfElasticity * PoissonsRatio / ((1 + PoissonsRatio) * (1 - PoissonsRatio));
            Lame2 = 0.5 / (1 + PoissonsRatio) * ModulusOfElasticity;
            Viscosity = 1;
        }
    }

    class Beam : Solid {

        public Beam() {
            PoissonsRatio = 0.5;
            ModulusOfElasticity = 300;
            Density = 10;
            Lame1 = ModulusOfElasticity * PoissonsRatio / ((1 + PoissonsRatio) * (1 - PoissonsRatio));
            Lame2 = 0.5 / (1 + PoissonsRatio) * ModulusOfElasticity;
        }
    }

    class ConvergenceTest : Solid {
        public ConvergenceTest() {
            PoissonsRatio = 0.5;
            ModulusOfElasticity = 3;
            Density = 1;
            Viscosity = 0;
            Lame1 = ModulusOfElasticity * PoissonsRatio / ((1 + PoissonsRatio) * (1 - PoissonsRatio));
            Lame2 = 0.5 / (1 + PoissonsRatio) * ModulusOfElasticity;
        }
    }

    public class AllOne : Solid {
        public AllOne() {
            PoissonsRatio = 1;
            ModulusOfElasticity = 1;
            Density = 1;
            Lame1 = 1;
            Lame2 = 1;
        }
    }
}
