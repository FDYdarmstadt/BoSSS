using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Security.Cryptography.X509Certificates;
using BoSSS.Application.BoSSSpad;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using ilPSP;
using ilPSP.Utils;
using static BoSSS.Solution.GridImport.NASTRAN.NastranFile;

namespace XESF.Fluxes {
    public abstract class RoeBaseFlux : IEdgeForm, IVolumeForm, ISupportsJacobianComponent {

        /// <summary>
        /// <see cref="HLLCDensityFlux.HLLCDensityFlux"/>
        /// </summary>
        protected IBoundaryConditionMap boundaryMap;

        protected Material material;
        protected int D;

        protected string speciesName;

        protected double s_alpha;
        /// <summary>
        /// Constructs a new flux
        /// </summary>
        /// <param name="boundaryMap">
        /// Mapping for boundary conditions
        /// </param>
        public RoeBaseFlux(IBoundaryConditionMap boundaryMap, Material material, double s_alpha, int D=2) {
            this.boundaryMap = boundaryMap;
            this.material = material;
            this.s_alpha = s_alpha;
            this.D = D;
        }
        public TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.UxV | TermActivationFlags.V;
            }
        }
        public TermActivationFlags InnerEdgeTerms {
            get {
                return TermActivationFlags.UxV;
            }
        }
        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.UxGradV;
            }
        }
        public IList<string> ArgumentOrdering {
            get {
                string[] variables;
                if(D == 1) {
                    variables = new string[] { CompressibleVariables.Density, CompressibleVariables.Momentum.xComponent, CompressibleVariables.Energy };
                } else if(D == 2) {
                    variables = new string[] { CompressibleVariables.Density, CompressibleVariables.Momentum.xComponent, CompressibleVariables.Momentum.yComponent, CompressibleVariables.Energy };
                } else {
                    variables = new string[] { CompressibleVariables.Density, CompressibleVariables.Momentum.xComponent, CompressibleVariables.Momentum.yComponent, CompressibleVariables.Momentum.zComponent, CompressibleVariables.Energy };
                }
                return variables;
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

        public double[] GetExternalState(ref CommonParamsBnd inp, double[] Uin) {
            double[] Uout = new double[Uin.Length];

            var xdgBoudaryMap = this.boundaryMap as XDGCompressibleBoundaryCondMap;
            var boundaryMap = this.boundaryMap as CompressibleBoundaryCondMap;
            if(xdgBoudaryMap == null && boundaryMap == null)
                throw new NotSupportedException("This type of boundary condition map is not supported.");

            //get EdgeTag
            byte EdgeTag = inp.EdgeTag;

            // Get boundary condition on this edge
            BoundaryCondition boundaryCondition;
            if(xdgBoudaryMap != null)
                boundaryCondition = xdgBoudaryMap.GetBoundaryConditionForSpecies(EdgeTag, this.speciesName);
            else
                boundaryCondition = boundaryMap.GetBoundaryCondition(EdgeTag);

            //create state vector from Uin
            Vector momentum = new Vector(inp.D);
            for(int i = 0; i < inp.D; i++) {
                momentum[i] = Uin[i + 1];
            }
            var stateIn = new StateVector(material, Uin[0], momentum, Uin[inp.D + 1]);

            //Call state evaluation
            var stateOut = boundaryCondition.GetBoundaryState(inp.time, inp.X, inp.Normal, stateIn);

            //Creat Uout from stateOut
            Uout[0] = stateOut.Density;
            for(int i = 0; i < inp.D; i++) {
                Uout[1 + i] = stateOut.Momentum[i];
            }
            Uout[inp.D + 1] = stateOut.Energy;

            return Uout;
        }
        //public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] Uin, double[,] _Grad_Uin, double Vin, double[] _Grad_Vin) {


        //    double[] Uout = GetExternalState(ref inp, Uin);


        //    var inpIn = new CommonParams();
        //    inpIn.Normal = inp.Normal;
        //    inpIn.X = inp.X;
        //    inpIn.time = inp.time;
        //    inpIn.iEdge = inp.iEdge;
        //    inpIn.GridDat = inp.GridDat;
        //    inpIn.EdgeTag = 0;
        //    //inpIn.D = inp.D;  <- doesn't work read only....



        //    //return InnerEdgeForm(ref inp, Uin, Uout, _Grad_Uin, _Grad_Uout, Vin,Vout, _Grad_Vin, _Grad_Vout); <- is this even possible somehow?
        //    return 0;
        //    }

        public (Vector vAverage, double speedOfSoundAverage, double enthalpyAverage, double[] eigenVals, MultidimensionalArray V0_inv) Setup(double[] normal, double[] Uin, double[] Uout, int D) {
        //public (MultidimensionalArray V0_inv, MultidimensionalArray eigenVals, MultidimensionalArray V0) Setup(double[] normal, double[] Uin, double[] Uout, int D)
            
            MultidimensionalArray V0_inv = MultidimensionalArray.Create(D + 2, D + 2);
            //MultidimensionalArray eigenVals = MultidimensionalArray.Create(D + 2, D + 2);
            //MultidimensionalArray V0 = MultidimensionalArray.Create(D + 2, D + 2);

            double densityIn = Uin[0];
            double densityOut = Uout[0];
            double[] vIn = new double[D];
            double[] vOut = new double[D];

            double momentumSquaredIn = 0;
            double momentumSquaredOut = 0;
            for(int d = 0; d < D; d++) {
                vIn[d] = Uin[d + 1] / densityIn;
                vOut[d] = Uout[d + 1] / densityOut;
                momentumSquaredIn += Math.Pow(Uin[d + 1], 2);
                momentumSquaredOut += Math.Pow(Uout[d + 1], 2);
            }
            double energyIn = Uin[D + 1];
            double energyOut = Uout[D + 1];

            // copied from HLLCEnergyFlux.cs, to calculate the pressure
            double gamma = this.material.EquationOfState.HeatCapacityRatio;
            double Mach = this.material.MachNumber;
            double gammaMachSquared = gamma * Mach * Mach;

            double pressureIn = (gamma - 1.0) * (energyIn - gammaMachSquared * 0.5 * momentumSquaredIn / densityIn);
            double pressureOut = (gamma - 1.0) * (energyOut - gammaMachSquared * 0.5 * momentumSquaredOut / densityOut);

            // Enthalpy: H = (E + pressure) / roh ; Toro p.347
            double enthalpyIn = (energyIn + pressureIn) / densityIn;
            double enthalpyOut = (energyOut + pressureOut) / densityOut;

            Vector vAverage = new Vector(D);

            for(int i = 0; i < D; i++) {
                // 0:uAverage, 1:vAverage, 2:wAverage
                vAverage[i] = (Math.Sqrt(densityIn) * vIn[i] + Math.Sqrt(densityOut) * vOut[i])
                    / (Math.Sqrt(densityIn) + Math.Sqrt(densityOut));
            }
            double velocityAverage = vAverage.L2Norm();

            double enthalpyAverage = (Math.Sqrt(densityIn) * enthalpyIn + Math.Sqrt(densityOut) * enthalpyOut)
                / (Math.Sqrt(densityIn) + Math.Sqrt(densityOut));

            double speedOfSoundAverage = Math.Sqrt((gamma - 1) * (enthalpyAverage - 0.5 * Math.Pow(velocityAverage, 2)));

            // Toro p.374, eq. (11.58), EigenValues[0...4]
            double vN = vAverage * (new Vector(normal));

            //for Eigenvals
            //s_alpha = 10;
            double SmoothedAbs(double x) {
                return x*Math.Tanh(s_alpha * x);
            }

            /*//Eigenvals
            eigenVals[0,0]= SmoothedAbs(vN - speedOfSoundAverage);
            for(int i = 1; i < D+1; i++) {
                // if D=2 or D=3 set the values; else it's = 0
                eigenVals[i,i] = SmoothedAbs(vN);
            }
            eigenVals[D+1, D+1] = SmoothedAbs(vN + speedOfSoundAverage);
            */
            //Eigenvals 
            double[] eigenVals = new double[D+2];
            eigenVals[0] = SmoothedAbs(vN - speedOfSoundAverage);
            for(int i = 1; i < D + 1; i++) {
                // if D=2 or D=3 set the values; else it's = 0
                eigenVals[i] = SmoothedAbs(vN);
            }
            eigenVals[D + 1] = SmoothedAbs(vN + speedOfSoundAverage);
            

            /*//+++++++++++++++++++ V0 +++++++++++++++++++++++
            //first row
            V0[0, 0] = 1;
            for(int j = 1; j < D + 1; j++) {
                V0[0, j] = normal[j - 1];
            }
            V0[0, D + 1] = 1;

            //midle rows
            for(int i = 1; i < D+1; i++) {
                //first col
                V0[i, 0] = vAverage[i-1] - speedOfSoundAverage * normal[i-1];

                //middle part
                for(int j = 1; j < D+1; j++) {
                    if(i == j) {
                        V0[i, j] = (vAverage[i - 1] - normal[i - 1]) * normal[j - 1] + 1;
                    } else {
                        V0[i, j] = (vAverage[i - 1] - normal[i - 1]) * normal[j - 1];
                    }
                }
                //last col
                V0[i, D+1] = vAverage[i - 1] + speedOfSoundAverage * normal[i - 1];
            }

            //last row
            V0[D + 1, 0] = enthalpyAverage - vN * speedOfSoundAverage;
            for(int j = 1; j < D + 1; j++) {
                V0[D + 1, j] = vAverage[j-1] + (vAverage.L2NormPow2()/2-vN)*normal[j - 1];
            }
            V0[D + 1, D + 1] = enthalpyAverage + vN * speedOfSoundAverage;
            
            //+++++++++++++++++++ END +++++++++++++++++++++++
            */

            //+++++++++++++++++++ V0_inv +++++++++++++++++++++++
            //first row
            double cdivgam = speedOfSoundAverage / (gamma - 1);
            V0_inv[0, 0] = vAverage.L2NormPow2()/2+vN* cdivgam;
            for(int j = 1; j < D + 1; j++) {
                V0_inv[0, j] = -vAverage[j-1] - cdivgam * normal[j - 1];
            }
            V0_inv[0, D + 1] = 1;

            //midle rows
            for(int i = 1; i < D + 1; i++) {
                //first col
                V0_inv[i, 0] = (2*speedOfSoundAverage* cdivgam) * ((vN+1) * normal[i - 1] - vAverage[i - 1]) - vAverage.L2NormPow2() * normal[i-1];

                //middle part
                for(int j = 1; j < D + 1; j++) {
                    if(i == j) {
                        V0_inv[i, j] = 2 * normal[i - 1] * vAverage[j - 1] +                //xxx i oder j jetzt I'm confused
                            (2 * speedOfSoundAverage * cdivgam) * (1 - normal[i - 1] * normal[j- 1]);
                    } else {
                        V0_inv[i, j] = 2 * normal[i - 1] * vAverage[j - 1] +
                            (2 * speedOfSoundAverage * cdivgam) * (-normal[i - 1] * normal[j - 1]);
                    }

                }
                //last col
                V0_inv[i, D+1] = -2 * normal[i - 1];
            }

            //last row
            V0_inv[D + 1, 0] = vAverage.L2NormPow2()/2 - vN * cdivgam;
            for(int j = 1; j < D + 1; j++) {
                V0_inv[D + 1, j] = -vAverage[j - 1] + cdivgam * normal[j - 1];
            }
            V0_inv[D + 1, D + 1] = 1;

            V0_inv.Scale((gamma - 1) / (2 * speedOfSoundAverage * speedOfSoundAverage));
            //+++++++++++++++++++ END +++++++++++++++++++++++

            
            return (vAverage, speedOfSoundAverage, enthalpyAverage, eigenVals, V0_inv);
        }
        public (MultidimensionalArray V0_inv, MultidimensionalArray eigenVals, MultidimensionalArray V0) OldSetup(double[] normal, double[] Uin, double[] Uout, int D) {

            MultidimensionalArray V0_inv = MultidimensionalArray.Create(D + 2, D + 2);
            MultidimensionalArray eigenVals = MultidimensionalArray.Create(D + 2, D + 2);
            MultidimensionalArray V0 = MultidimensionalArray.Create(D + 2, D + 2);

            double densityIn = Uin[0];
            double densityOut = Uout[0];
            double[] vIn = new double[D];
            double[] vOut = new double[D];

            double momentumSquaredIn = 0;
            double momentumSquaredOut = 0;
            for(int d = 0; d < D; d++) {
                vIn[d] = Uin[d + 1] / densityIn;
                vOut[d] = Uout[d + 1] / densityOut;
                momentumSquaredIn += Math.Pow(Uin[d + 1], 2);
                momentumSquaredOut += Math.Pow(Uout[d + 1], 2);
            }
            double energyIn = Uin[D + 1];
            double energyOut = Uout[D + 1];

            // copied from HLLCEnergyFlux.cs, to calculate the pressure
            double gamma = this.material.EquationOfState.HeatCapacityRatio;
            double Mach = this.material.MachNumber;
            double gammaMachSquared = gamma * Mach * Mach;

            double pressureIn = (gamma - 1.0) * (energyIn - gammaMachSquared * 0.5 * momentumSquaredIn / densityIn);
            double pressureOut = (gamma - 1.0) * (energyOut - gammaMachSquared * 0.5 * momentumSquaredOut / densityOut);

            // Enthalpy: H = (E + pressure) / roh ; Toro p.347
            double enthalpyIn = (energyIn + pressureIn) / densityIn;
            double enthalpyOut = (energyOut + pressureOut) / densityOut;

            Vector vAverage = new Vector(D);

            for(int i = 0; i < D; i++) {
                // 0:uAverage, 1:vAverage, 2:wAverage
                vAverage[i] = (Math.Sqrt(densityIn) * vIn[i] + Math.Sqrt(densityOut) * vOut[i])
                    / (Math.Sqrt(densityIn) + Math.Sqrt(densityOut));
            }
            double velocityAverage = vAverage.L2Norm();

            double enthalpyAverage = (Math.Sqrt(densityIn) * enthalpyIn + Math.Sqrt(densityOut) * enthalpyOut)
                / (Math.Sqrt(densityIn) + Math.Sqrt(densityOut));

            double speedOfSoundAverage = Math.Sqrt((gamma - 1) * (enthalpyAverage - 0.5 * Math.Pow(velocityAverage, 2)));

            // Toro p.374, eq. (11.58), EigenValues[0...4]
            double vN = vAverage * (new Vector(normal));

            double s_alpha = 10;
            double SmoothedAbs(double x) {
                return x * Math.Tanh(s_alpha * x);
            }
            //Eigenvals
            eigenVals[0, 0] = SmoothedAbs(vN - speedOfSoundAverage);
            for(int i = 1; i < D + 1; i++) {
                // if D=2 or D=3 set the values; else it's = 0
                eigenVals[i, i] = SmoothedAbs(vN);
            }
            eigenVals[D + 1, D + 1] = SmoothedAbs(vN + speedOfSoundAverage);

            //+++++++++++++++++++ V0 +++++++++++++++++++++++
            //first row
            V0[0, 0] = 1;
            for(int j = 1; j < D + 1; j++) {
                V0[0, j] = normal[j - 1];
            }
            V0[0, D + 1] = 1;

            //midle rows
            for(int i = 1; i < D + 1; i++) {
                //first col
                V0[i, 0] = vAverage[i - 1] - speedOfSoundAverage * normal[i - 1];

                //middle part
                for(int j = 1; j < D + 1; j++) {
                    if(i == j) {                                        
                        V0[i, j] = (vAverage[i - 1] - normal[i - 1]) * normal[j - 1] + 1;
                    } else {
                        V0[i, j] = (vAverage[i - 1] - normal[i - 1]) * normal[j - 1];
                    }
                }
                //last col
                V0[i, D + 1] = vAverage[i - 1] + speedOfSoundAverage * normal[i - 1];
            }

            //last row
            V0[D + 1, 0] = enthalpyAverage - vN * speedOfSoundAverage;
            for(int j = 1; j < D + 1; j++) {
                V0[D + 1, j] = vAverage[j - 1] + (vAverage.L2NormPow2() / 2 - vN) * normal[j - 1];
            }
            V0[D + 1, D + 1] = enthalpyAverage + vN * speedOfSoundAverage;
            //+++++++++++++++++++ END +++++++++++++++++++++++


            //+++++++++++++++++++ V0_inv +++++++++++++++++++++++
            //first row
            double cdivgam = speedOfSoundAverage / (gamma - 1);
            V0_inv[0, 0] = vAverage.L2NormPow2() / 2 + vN * cdivgam;
            for(int j = 1; j < D + 1; j++) {
                V0_inv[0, j] = -vAverage[j - 1] - cdivgam * normal[j - 1];
            }
            V0_inv[0, D + 1] = 1;

            //midle rows
            for(int i = 1; i < D + 1; i++) {
                //first col
                V0_inv[i, 0] = (2 * speedOfSoundAverage * cdivgam) * ((vN + 1) * normal[i - 1] - vAverage[i - 1]) - vAverage.L2NormPow2() * normal[i - 1];

                //middle part
                for(int j = 1; j < D + 1; j++) {
                    if(i == j) {
                        V0_inv[i, j] = 2 * normal[i - 1] * vAverage[j - 1] +
                            (2 * speedOfSoundAverage * cdivgam) * (1 - normal[i - 1] * normal[j - 1]);
                    } else {
                        V0_inv[i, j] = 2 * normal[i - 1] * vAverage[j - 1] +
                            (2 * speedOfSoundAverage * cdivgam) * (-normal[i - 1] * normal[j - 1]);
                    }

                }
                //last col
                V0_inv[i, D + 1] = -2 * normal[i - 1];
            }

            //last row
            V0_inv[D + 1, 0] = vAverage.L2NormPow2() / 2 - vN * cdivgam;
            for(int j = 1; j < D + 1; j++) {
                V0_inv[D + 1, j] = -vAverage[j - 1] + cdivgam * normal[j - 1];
            }
            V0_inv[D + 1, D + 1] = 1;

            V0_inv.Scale((gamma - 1) / (2 * speedOfSoundAverage * speedOfSoundAverage));
            //+++++++++++++++++++ END +++++++++++++++++++++++
            return (V0_inv, eigenVals, V0);
        }

        // public double InnerEdgeFlux(double[] normal, double[] Uin, double[] Uout, int D) ==>
        public (double[] waveStrength,double[] eigenVals, double[,] eigenVec) SetupXSplit(double[] normal, double[] Uin, double[] Uout, int D) {

            // xxxx mal probieren die Variablen alle hier rein zu bringen -- vllt über static !?
            //---------------------------------------------------------------------------------- 
            // Uin[i] and Uout[i] and Averages
            double densityIn = Uin[0];
            double densityOut = Uout[0];
            /*
            vIn/OUT => velocity In/Out => 0:u, 1:v, 2:w
            momentum => 0:x, 1:y, 2:z
            Note: if D = 2 or 3, the value for the respective Index remains 0 (e.x. v=0,w=0)
                and the overall dimension for later inputs are same for all dimensions 
            */
            double[] vIn = new double[3];
            double[] vOut = new double[3];
            Vector momentumIn = new Vector(3);
            Vector momentumOut = new Vector(3);
            double momentumSquaredIn = 0;
            double momentumSquaredOut = 0;
            for(int d = 0; d < D; d++) {
                momentumIn[d] = Uin[d + 1];
                momentumOut[d] = Uout[d + 1];
                //vIn[d] = momentumIn[d] / densityIn;
                //vOut[d] = momentumOut[d] / densityOut;
                vIn[d] = momentumIn[d] * normal[d] / densityIn; // xxx maybe wrong
                vOut[d] = momentumOut[d] * normal[d] / densityOut;
                momentumSquaredIn += Math.Pow(momentumIn[d], 2);
                momentumSquaredOut += Math.Pow(momentumOut[d], 2);
            }
            double energyIn = Uin[D + 1];
            double energyOut = Uout[D + 1];

            double[] waveStrength = new double[5];

            // vAverage => 0:uAverage, 1:vAverage, 2:wAverage 
            double[] vAverage = new double[3];

            // 5x Eigenvalues/Vectors (U = [mass,x,y,z-mom.,energy] )
            double[] eigenVal = new double[5];

            // copied from HLLCEnergyFlux.cs, to calculate the pressure
            double gamma = this.material.EquationOfState.HeatCapacityRatio;
            double Mach = this.material.MachNumber;
            double gammaMachSquared = gamma * Mach * Mach;

            double pressureIn = (gamma - 1.0) * (energyIn - gammaMachSquared * 0.5 * momentumSquaredIn / densityIn);
            double pressureOut = (gamma - 1.0) * (energyOut - gammaMachSquared * 0.5 * momentumSquaredOut / densityOut);

            // Enthalpy: H = (E + pressure) / roh ; Toro p.347
            double enthalpyIn = (energyIn + pressureIn) / densityIn;
            double enthalpyOut = (energyOut + pressureOut) / densityOut;

            /* Toro p.375, eq. (11.60), Averages (u,v,w,H,a[speed of sound])
            Note: if D = 2 or 3, the value for the respective Index remains 0 (e.x. v=0,w=0)
                and the overall dimension for later inputs are same for all dimensions 
            */

            for(int i = 0; i < D; i++) {
                // 0:uAverage, 1:vAverage, 2:wAverage
                vAverage[i] = (Math.Sqrt(densityIn) * vIn[i] + Math.Sqrt(densityOut) * vOut[i])
                    / (Math.Sqrt(densityIn) + Math.Sqrt(densityOut));
            }
            double enthalpyAverage = (Math.Sqrt(densityIn) * enthalpyIn + Math.Sqrt(densityOut) * enthalpyOut)
                / (Math.Sqrt(densityIn) + Math.Sqrt(densityOut));
            double velocityAverage = Math.Sqrt(Math.Pow(vAverage[0], 2) + Math.Pow(vAverage[1], 2) + Math.Pow(vAverage[2], 2));
            double speedOfSoundAverage = Math.Pow((gamma - 1) * (enthalpyAverage - 0.5 * Math.Pow(velocityAverage, 2)), 0.5);

            // Toro p.374, eq. (11.58), EigenValues[0...4]
            eigenVal[0] = vAverage[0] - speedOfSoundAverage;
            eigenVal[1] = vAverage[0];
            for(int i = 2; i <= D; i++) {
                // if D=2 or D=3 set the values; else it's = 0
                eigenVal[i] = eigenVal[1];
            }
            eigenVal[4] = vAverage[0] + speedOfSoundAverage;


            // ------------------------------------------------------------------------
            // Toro p.375, eq. (11.59), averaged right Eigenvector
            double[,] eigenVec = new double[5, 5] {
                { 1, vAverage[0] - speedOfSoundAverage, vAverage[1], vAverage[2],
                enthalpyAverage - vAverage[0] * speedOfSoundAverage },
                { 1, vAverage[0], vAverage[1], vAverage[2], 0.5 * Math.Pow(vAverage[1], 2) },
                { 0, 0, 1, 0, vAverage[1] },
                { 0, 0, 0, 1, vAverage[2] },
                { 1, vAverage[0] + speedOfSoundAverage, vAverage[1],
                vAverage[2], enthalpyAverage + vAverage[0] * speedOfSoundAverage }
            };
            /*xxx TEST TEST TEST TEST TEST TEST TEST TEST 
            Vector eigenVec_0 = new Vector { 1, vAverage[0] - speedOfSoundAverage, vAverage[1], vAverage[2],
                enthalpyAverage - vAverage[0] * speedOfSoundAverage };
            Vector eigenVec_1 = new Vector { 1, vAverage[0], vAverage[1], vAverage[2], 0.5 * Math.Pow(vAverage[1], 2) };
            Vector eigenVec_2 = new Vector { 0, 0, 1, 0, vAverage[1] };
            Vector eigenVec_3 = new Vector { 0, 0, 0, 1, vAverage[2] };
            Vector eigenVec_4 = new Vector  { 1, vAverage[0] + speedOfSoundAverage, vAverage[1],
                vAverage[2], enthalpyAverage + vAverage[0] * speedOfSoundAverage }; 
            */ // --------------------------------------------------------------------------

            // Toro p.376 eq.(11.68-11.70) "wave strengths"
            double u5_medium = (energyOut - energyIn) - ((momentumOut[1] - momentumIn[1]) - vAverage[1] * (densityOut - densityIn)) * vAverage[1] -
                ((momentumOut[2] - momentumIn[2]) - vAverage[2] * (densityOut - densityIn)) * vAverage[2];
            // xxx für u5_medium => kp wie der Strich heißt, wenn es nicht average heißt
            // waveStrength [0...4] =>Averaged !!
            waveStrength[2] = (momentumOut[1] - momentumIn[1]) - vAverage[1] * (densityOut - densityIn);
            waveStrength[3] = (momentumOut[2] - momentumIn[2]) - vAverage[2] * (densityOut - densityIn);
            waveStrength[1] = (gamma - 1) / Math.Pow(speedOfSoundAverage, 2) *
                ((momentumOut[0] - momentumIn[0]) * (enthalpyAverage - Math.Pow(vAverage[0], 2) +
                vAverage[0] * (momentumOut[0] - momentumIn[0]) - u5_medium));
            waveStrength[0] = (1 / (2 * speedOfSoundAverage)) * ((densityOut - densityIn) * (vAverage[0] + speedOfSoundAverage) -
                (momentumOut[0] - momentumIn[0]) - speedOfSoundAverage * waveStrength[1]);
            waveStrength[4] = (densityOut - densityIn) - (waveStrength[0] + waveStrength[1]);
            // -------------------------------------------------------------------------------------
            return (waveStrength, eigenVal, eigenVec);
        }



        public abstract double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] Uin, double[,] _Grad_Uin, double Vin, double[] _Grad_Vin);


        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            if(SpatialDimension != 2)
                throw new NotImplementedException("Only supporting 2D.");

            return new IEquationComponent[] {
                new EdgeFormDifferentiator(this, SpatialDimension),
                new VolumeFormDifferentiator(this, SpatialDimension)
            };
        }

        public abstract double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT);

        public abstract double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV);
    }
}
