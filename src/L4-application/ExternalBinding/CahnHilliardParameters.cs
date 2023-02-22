namespace BoSSS.Application.ExternalBinding {
    public readonly struct CahnHilliardParameters {
        public CahnHilliardParameters(double _cahn = 0.1, double _diffusion = 0.1){
            Cahn = _cahn;
            Diffusion = _diffusion;
        }
        public double Cahn {get;}
        public double Diffusion {get;}
    }
}
