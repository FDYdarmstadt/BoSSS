namespace BoSSS.Application.ExternalBinding {

    /// <summary>
    /// In this struct, the parameters controlling the Cahn-Hilliard calculation using FixedOperators are set.
    /// One may think of this as a poor man's control object for this solver.
    /// </summary>
    public readonly struct CahnHilliardParameters {
        public CahnHilliardParameters(double _cahn = 0.1, double _diffusion = 0.1){
            Cahn = _cahn;
            Diffusion = _diffusion;
        }
        public double Cahn {get;}
        public double Diffusion {get;}
    }
}
