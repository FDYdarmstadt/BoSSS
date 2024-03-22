using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Security.Cryptography.X509Certificates;
using BoSSS.Application.BoSSSpad;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.CompressibleFlowCommon;
using BoSSS.Solution.CompressibleFlowCommon.Boundary;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;
using ilPSP;
using ilPSP.Utils;
using XESF.Fluxes;
using static BoSSS.Solution.GridImport.NASTRAN.NastranFile;

namespace XESTSF.Fluxes {
    public abstract class RoeSTBaseFlux : IEdgeForm, IVolumeForm, ISupportsJacobianComponent {

        /// <summary>
        /// <see cref="OptimizedHLLCDensityFlux.OptimizedHLLCDensityFlux"/>
        /// </summary>
        protected IBoundaryConditionMap boundaryMap;
        /// <summary>
        /// <see cref="Material"/>
        /// </summary>
        protected Material material;
        /// <summary>
        /// Dimension
        /// </summary>
        protected int D;
        /// <summary>
        /// XDGField which is used to prescribe boundary values, if space-time boundary condition is chosen
        /// </summary>
        protected DGField[] previous_u;
        /// <summary>
        /// Name of species this flux is used
        /// </summary>
        protected string speciesName;
        
        /// <summary>
        /// smoothing parameter used to smooth the absolute value in Roe Flux ( |EigenVal| )
        /// </summary>
        protected double s_alpha=10;
        /// <summary>
        /// Constructs a new flux
        /// </summary>
        /// <param name="boundaryMap">
        /// Mapping for boundary conditions
        /// </param>
        public RoeSTBaseFlux(IBoundaryConditionMap boundaryMap, Material material, double s_alpha, DGField[] _previous_u, int D=2) {
            this.boundaryMap = boundaryMap;
            this.material = material;
            this.s_alpha = s_alpha;
            this.D = D;
            this.previous_u = _previous_u;

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

        /// <summary>
        /// Obtains the externatl state Uout for boudary edges depending on boundary-condition-type chosen and using prescribed functions/fields
        /// </summary>
        /// <param name="inp"></param>
        /// <param name="Uin"></param>
        /// <returns></returns>
        /// <exception cref="NotSupportedException"></exception>
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

            if(boundaryCondition is SpaceTimeBoundaryCondition) {
                Uout = EvalPrevState(inp);
                return Uout;
            } else {
                //create state vector from Uin
                Vector momentum = new Vector(inp.D - 1);
                for(int i = 0; i < inp.D - 1; i++) {
                    momentum[i] = Uin[i + 1];
                }
                var stateIn = new StateVector(material, Uin[0], momentum, Uin[inp.D]);

                //Call state evaluation

                //reduce dimension
                double[] xRed = new double[inp.D - 1];
                double[] nRed = new double[inp.D - 1];
                for(int i = 0; i < inp.D - 1; i++) {
                    xRed[i] = inp.X[i];
                    nRed[i] = inp.Normal[i];
                }

                var stateOut = boundaryCondition.GetBoundaryState(inp.X[inp.D - 1], xRed, nRed, stateIn);
                //Creat Uout from stateOut
                Uout[0] = stateOut.Density;
                for(int i = 0; i < inp.D - 1; i++) {
                    Uout[1 + i] = stateOut.Momentum[i];
                }
                Uout[inp.D] = stateOut.Energy;

                return Uout;
            }
            
        }
        /// <summary>
        /// Evaluates the previous state u_previous at the point inp.X
        /// </summary>
        /// <param name="inp"></param>
        /// <returns></returns>
        private double[] EvalPrevState(CommonParamsBnd inp) {
            
            double[] Uout = new double[inp.D+1];

            //space time point to evaluate
            var x = inp.X;
            double[] xarr;
            if(x.Dim == 1) {
                xarr = new double[] { x.x };
            } else if(x.Dim == 2) {
                xarr= new double[] { x.x, x.y };
            } else {
                xarr=new double[] { x.x, x.y,x.z };
            }
            var grddat = (GridData)previous_u[0].GridDat;

            bool CheckIfPointInRefSquare(MultidimensionalArray _pt_local){
                double tol = 1e-14;
                if(-1.0 - tol < _pt_local[0, 0] && _pt_local[0, 0] < 1.0 + tol) {
                    if(-1.0 - tol < _pt_local[0, 1] && _pt_local[0, 1] < 1.0 + tol) {
                        return true;
                    }
                }
                return false;
            }
            ////first searching the cell and then evaluating
            foreach(var cell in grddat.GetBoundaryCells().ItemEnum) {
                    MultidimensionalArray _pt = MultidimensionalArray.CreateWrapper(xarr, 1, D+1);          // point to search for
                    MultidimensionalArray _pt_local = MultidimensionalArray.CreateWrapper(xarr.CloneAs(), 1, D+1);
                    bool[] bb = new bool[1];
                    grddat.TransformGlobal2Local(_pt, _pt_local, cell, bb);
                    
                if(CheckIfPointInRefSquare(_pt_local)) {
                    for(int iField = 0; iField < previous_u.Length; iField++) {
                        NodeSet NS = new NodeSet(previous_u[iField].GridDat.iGeomCells.RefElements[0], _pt_local, true);
                        var result = MultidimensionalArray.Create(1, 1);
                        previous_u[iField].Evaluate(cell, 1, NS, result);
                        Uout[iField] = result[0, 0];
                    }
                    return Uout;
                }
            }

            //using global evaluate
            //for(int iField=0;iField < previous_u.Length; iField++) {
            //    NodeSet NSglob = new NodeSet(previous_u[iField].GridDat.iGeomCells.RefElements[0], x, true);
            //    var ret = previous_u[iField].Evaluate(NSglob, grddat.GetBoundaryCells());
            //    Uout[iField] = ret[0, 0];
            //}
            Console.WriteLine("SpaceTimeBoundary - no cell found for point x=(" + x.x+","+x.y + " Point not found - Uout is zero");
            return Uout;

        }

        /// <summary>
        /// Computes some quantities relevant for Roe flux:
        /// - eigenVals - eigen values of projected Jacobian
        /// ....
        /// </summary>
        /// <param name="normal"></param>
        /// <param name="Uin"></param>
        /// <param name="Uout"></param>
        /// <param name="D"></param>
        /// <param name="n_t"></param>
        /// <returns></returns>
        public (Vector vAverage, double speedOfSoundAverage, double enthalpyAverage, double[] eigenVals, MultidimensionalArray V0_inv) Setup(double[] normal, double[] Uin, double[] Uout, int D,double n_t) {
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
                //Debug.Assert(!vAverage[i].IsNaN());
            }
            double velocityAverage = vAverage.L2Norm();
            
            double enthalpyAverage = (Math.Sqrt(densityIn) * enthalpyIn + Math.Sqrt(densityOut) * enthalpyOut)
                / (Math.Sqrt(densityIn) + Math.Sqrt(densityOut));
            
            double speedOfSoundAverage = Math.Sqrt((gamma - 1) * (enthalpyAverage - 0.5 * Math.Pow(velocityAverage, 2)));

            //Debug.Assert(!velocityAverage.IsNaN());
            //Debug.Assert(!enthalpyAverage.IsNaN());
            //Debug.Assert(!speedOfSoundAverage.IsNaN());
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
            eigenVals[0] = SmoothedAbs((vN - speedOfSoundAverage)*normal.L2Norm() + n_t);
            for(int i = 1; i < D + 1; i++) {
                // if D=2 or D=3 set the values; else it's = 0
                eigenVals[i] = SmoothedAbs(vN * normal.L2Norm() + n_t);
            }
            eigenVals[D + 1] = SmoothedAbs((vN + speedOfSoundAverage) * normal.L2Norm() + n_t);
            

            /*//+++++++++++++++++++ Full V0 +++++++++++++++++++++++
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

        public abstract double InnerEdgeFlux(double[] normal, double[] Uin, double[] Uout, int D);

        public abstract Vector Flux(double[] U, int D);

        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] Uin, double[,] _Grad_Uin, double Vin, double[] _Grad_Vin) {
            // Copied from HLLCEnergyFlux
            double[] Uout = GetExternalState(ref inp, Uin);
            double ret = InnerEdgeFlux(inp.Normal, Uin, Uout, inp.D);
            return ret * (Vin);
        }


        public double InnerEdgeForm(ref CommonParams inp, double[] Uin, double[] Uout, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            double ret = InnerEdgeFlux(inp.Normal, Uin, Uout, inp.D);
            return ret * (_vIN - _vOUT);
        }

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            Vector fluxvector = Flux(U, cpv.D);
            return -1.0 * fluxvector * GradV;
        }


        public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            if(SpatialDimension != 2)
                throw new NotImplementedException("Only supporting 2D.");

            return new IEquationComponent[] {
                new EdgeFormDifferentiator(this, SpatialDimension),
                new VolumeFormDifferentiator(this, SpatialDimension)
            };
        }
    }
}
