namespace BoSSS.Application.ExternalBinding {

    /// <summary>
    /// In this struct, the parameters controlling the Cahn-Hilliard calculation using FixedOperators are set.
    /// One may think of this as a poor man's control object for this solver.
    /// </summary>
    public readonly struct CahnHilliardParameters {
        public CahnHilliardParameters(double _cahn = 0.1, double _diffusion = 0.1, double _dt = 1e5, double _endT = 1.5e5, bool _stationary = false){
            Cahn = _cahn;
            Diffusion = _diffusion;
            dt = _dt;
            endT = _endT;
            if (_stationary){
                dt = 1e5;
                endT = 1.5*dt;
            }
            // Convection = _convection;
        }
        public double Cahn {get;}
        public double Diffusion {get;}
        public double dt {get;}
        public double endT {get;}
        // public bool Convection {get;}
    }
}
