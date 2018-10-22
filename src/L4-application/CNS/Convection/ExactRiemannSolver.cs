/* =======================================================================
Copyright 2017 Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

using BoSSS.Platform.LinAlg;
using BoSSS.Solution.CompressibleFlowCommon;
using CNS.MaterialProperty;
using System;
using BoSSS.Solution.CompressibleFlowCommon.MaterialProperty;

namespace CNS.Convection {

    /// <summary>
    /// C# port of E.F. Toro's exact Riemann solver for the Euler equations
    /// originally written in Fortran taken from
    /// Toro, E.F - "Riemann Solvers and Numerical Methods for Fluid Dynamics"
    /// </summary>
    public class ExactRiemannSolver {

        /// <summary>
        /// Termination condition for the solution of the pressure in the
        /// star region in <see cref="GetStarRegionValues"/>.
        /// </summary>
        private const double PRESSURE_TOLERANCE = 1e-6;

        /// <summary>
        /// Maximum number of Newton-Raphson iterations in
        /// <see cref="GetStarRegionValues"/>
        /// </summary>
        private const double MAX_PRESSURE_ITERATIONS = 20;

        /// <summary>
        /// Heat capacity ratios for gas on the left and on the right
        /// </summary>
        private double kappaLeft, kappaRight;

        /// <summary>
        /// Various intermediate values for the gas on the left
        /// </summary>
        private double G1L, G2L, G3L, G4L, G5L, G6L, G7L, G8L;

        /// <summary>
        /// Various intermediate values for the gas on the right
        /// </summary>
        private double G1R, G2R, G3R, G4R, G5R, G6R, G7R, G8R;

        /// <summary>
        /// Density, velocity, pressure and speed of sound on the left
        /// </summary>
        private double densityLeft, normalVelocityLeft, pressureLeft, speedOfSoundLeft;

        /// <summary>
        /// Density, velocity and pressure on the right
        /// </summary>
        private double densityRight, normalVelocityRight, pressureRight, speedOfSoundRight;

        /// <summary>
        /// Density, velocity and pressure in the star region
        /// </summary>
        private double densityStar, normalVelocityStar, pressureStar;

        /// <summary>
        /// Some mean values
        /// </summary>

        /// <summary>
        /// Full states left and right of the discontinuity
        /// </summary>
        private StateVector stateLeft, stateRight;

        /// <summary>
        /// Velocities in the edge coordinate system defined by the given edge
        /// normal.
        /// </summary>
        private Vector velocityLeft, velocityRight, edgeNormal;

        /// <summary>
        /// Constructor for the multiphase case, initializes values
        /// </summary>
        /// <param name="stateLeft">
        /// State left of the discontinuity
        /// </param>
        /// <param name="stateRight">
        /// State right of the discontinuity
        /// </param>
        /// <param name="edgeNormal">
        /// The normal vector of the edge between the 'left' and the 'right'
        /// state (pointing from 'left' to 'right'). 
        /// </param>
        public ExactRiemannSolver(StateVector stateLeft, StateVector stateRight, Vector edgeNormal) {
            IdealGas gasLeft = stateLeft.Material.EquationOfState as IdealGas;
            IdealGas gasRight = stateRight.Material.EquationOfState as IdealGas;
            if (gasLeft == null || gasRight == null) {
                throw new ArgumentException(
                    "Exact Riemann solver only implemented for ideal gases");
            }

            this.edgeNormal = edgeNormal;
            this.stateLeft = stateLeft;
            this.stateRight = stateRight;
            this.kappaLeft = gasLeft.HeatCapacityRatio;
            this.kappaRight = gasRight.HeatCapacityRatio;
            this.densityLeft = stateLeft.Density;
            this.pressureLeft = stateLeft.Pressure;
            this.densityRight = stateRight.Density;
            this.pressureRight = stateRight.Pressure;

            velocityLeft = stateLeft.Velocity.ToEdgeCoordinates(edgeNormal);
            velocityRight = stateRight.Velocity.ToEdgeCoordinates(edgeNormal);
            this.normalVelocityLeft = velocityLeft[0];
            this.normalVelocityRight = velocityRight[0];

            speedOfSoundLeft = Math.Sqrt(kappaLeft * pressureLeft / densityLeft);
            speedOfSoundRight = Math.Sqrt(kappaRight * pressureRight / densityRight);

            // left side
            G1L = (kappaLeft - 1.0) / (2.0 * kappaLeft);
            G2L = (kappaLeft + 1.0) / (2.0 * kappaLeft);
            G3L = 2.0 * kappaLeft / (kappaLeft - 1.0);
            G4L = 2.0 / (kappaLeft - 1.0);
            G5L = 2.0 / (kappaLeft + 1.0);
            G6L = (kappaLeft - 1.0) / (kappaLeft + 1.0);
            G7L = (kappaLeft - 1.0) / 2.0;
            G8L = kappaLeft - 1.0;

            // right side
            G1R = (kappaRight - 1.0) / (2.0 * kappaRight);
            G2R = (kappaRight + 1.0) / (2.0 * kappaRight);
            G3R = 2.0 * kappaRight / (kappaRight - 1.0);
            G4R = 2.0 / (kappaRight - 1.0);
            G5R = 2.0 / (kappaRight + 1.0);
            G6R = (kappaRight - 1.0) / (kappaRight + 1.0);
            G7R = (kappaRight - 1.0) / 2.0;
            G8R = kappaRight - 1.0;

            // the pressure positivity condition is tested for
            if (G4L * (speedOfSoundLeft + speedOfSoundRight) <= (normalVelocityRight - normalVelocityLeft)) {
                throw new ArgumentOutOfRangeException("the initial data is such that vacuum is generated");
            }
        }

        /// <summary>
        /// Computes the exact solution of the Riemann problem, samples it at
        /// <paramref name="x"/> (coordinate normal to discontinuity position)
        /// at time <paramref name="t"/>
        /// </summary>
        /// <param name="meanPressure"></param>
        /// <param name="meanVelocity"></param>
        /// <param name="x"></param>
        /// <param name="t"></param>
        /// <returns></returns>
        public StateVector GetState(double meanPressure, double meanVelocity, double x, double t) {
            double S = x / t;
            Sample(meanPressure, meanVelocity, S, out densityStar, out normalVelocityStar, out pressureStar);

            // Return exact solution in conservative variables
            Vector u;
            Material material;
            if (S <= meanVelocity) {
                // left of the interface
                u = velocityLeft;
                material = stateLeft.Material;
            } else {
                // right of the interface
                u = velocityRight;
                material = stateRight.Material;
            }
            u[0] = normalVelocityStar;

            return StateVector.FromPrimitiveQuantities(
                material,
                densityStar,
                u.FromEdgeCoordinates(edgeNormal),
                pressureStar);
        }

        /// <summary>
        /// Computes the exact solution of the Riemann problem, samples it at
        /// x=0.
        /// </summary>
        /// <returns>
        /// The density, the momentum and the energy at x=0.
        /// </returns>
        public StateVector GetCentralState() {
            double meanPressure, meanVelocity;
            GetStarRegionValues(out meanPressure, out meanVelocity);
            return GetState(meanPressure, meanVelocity, 0.0, 1.0);
        }

        /// <summary>
        /// compute the solution for pressure and velocity in the Star Region
        /// </summary>
        /// <param name="P"></param>
        /// <param name="U"></param>
        /// <returns></returns>
        public void GetStarRegionValues(out double P, out double U) {
            double CHANGE, FL, FLD, FR, FRD, POLD, UDIFF, PCOMP;
            CHANGE = double.NaN;
            FL = double.NaN;
            FR = double.NaN;
            PCOMP = double.NaN;

            // guessed value PSTART is computed
            double PSTART = EstimatePressure();

            POLD = PSTART;
            UDIFF = normalVelocityRight - normalVelocityLeft;

            for (int i = 1; i <= MAX_PRESSURE_ITERATIONS; i++) {
                EvaluatePressureFunctions(out FL, out FLD, POLD, densityLeft, pressureLeft, speedOfSoundLeft, true);
                EvaluatePressureFunctions(out FR, out FRD, POLD, densityRight, pressureRight, speedOfSoundRight, false);

                PCOMP = POLD - (FL + FR + UDIFF) / (FLD + FRD);

                CHANGE = 2.0 * Math.Abs((PCOMP - POLD) / (PCOMP + POLD));
                if (CHANGE <= PRESSURE_TOLERANCE) {
                    break;
                }
                if (PCOMP < 0.0) {
                    PCOMP = PRESSURE_TOLERANCE;
                }
                POLD = PCOMP;
            }

            if (CHANGE > PRESSURE_TOLERANCE) {
                throw new Exception("No convergence in Newton-Raphson iteration");
            }

            if (double.IsNaN(CHANGE) | double.IsNaN(FL) | double.IsNaN(FR) | double.IsNaN(PCOMP)) {
                throw new Exception("Error in Iteration");
            }

            // Compute velocity in Star Region
            U = 0.5 * (normalVelocityLeft + normalVelocityRight + FR - FL);
            P = PCOMP;
        }

        /// <summary>
        /// Provides a guess value for the pressure PM in the Star Region. The
        /// choice is made according to adaptive Riemann solver using the PVRS,
        /// TRRS, and TSRS approximate Riemann solvers.
        /// See Sect. 9.5 of Chapt. 9 of Ref. 1
        /// </summary>
        /// <returns>
        /// An estimate for the pressure in the star region
        /// </returns>
        private double EstimatePressure() {
            double CUP, GEL, GER, PMAX, PMIN, PPV, PQ, PTL, PTR, QMAX, QUSER, UM, P;

            QUSER = 2.0;

            //  Compute guess pressure from PVRS Riemann solver
            CUP = 0.25 * (densityLeft + densityRight) * (speedOfSoundLeft + speedOfSoundRight);
            PPV = 0.5 * (pressureLeft + pressureRight) + 0.5 * (normalVelocityLeft - normalVelocityRight) * CUP;
            PPV = Math.Max(0.0, PPV);
            PMIN = Math.Min(pressureLeft, pressureRight);
            PMAX = Math.Max(pressureLeft, pressureRight);
            QMAX = PMAX / PMIN;

            if (QMAX <= QUSER && (PMIN <= PPV && PPV <= PMAX)) {
                //slect PVRS Riemann solver
                P = PPV;
            } else {
                if (PPV < PMIN) {
                    // Select Two-Rarefaction Riemann solver

                    PQ = Math.Pow((pressureLeft / pressureRight), G1L);
                    UM = (PQ * normalVelocityLeft / speedOfSoundLeft + normalVelocityRight / speedOfSoundRight + G4L * (PQ - 1.0)) / (PQ / speedOfSoundLeft + 1.0 / speedOfSoundRight);
                    PTL = 1.0 + G7L * (normalVelocityLeft - UM) / speedOfSoundLeft;
                    PTR = 1.0 + G7R * (UM - normalVelocityRight) / speedOfSoundRight;
                    P = 0.5 * (pressureLeft * Math.Pow(PTL, G3L) + pressureRight * Math.Pow(PTR, G3R));
                } else {
                    // Select two-shock riemann solver with PVRS as estimate
                    GEL = Math.Sqrt((G5L / densityLeft) / (G6L * pressureLeft + PPV));
                    GER = Math.Sqrt((G5R / densityRight) / (G6R * pressureRight + PPV));
                    P = (GEL * pressureLeft + GER * pressureRight - (normalVelocityRight - normalVelocityLeft)) / (GEL + GER);
                }
            }

            return P;
        }

        /// <summary>
        /// evaluate the pressure functions FL and FR in exact Riemann solver
        /// </summary>
        private void EvaluatePressureFunctions(out double F, out double FD, double P, double DK, double PK, double CK, bool leftSide) {
            double AK, BK, PRAT, QRT;
            double G1, G2, G4, G5, G6;
            if (leftSide) {
                G1 = G1L;
                G2 = G2L;
                G4 = G4L;
                G5 = G5L;
                G6 = G6L;
            } else {
                G1 = G1R;
                G2 = G2R;
                G4 = G4R;
                G5 = G5R;
                G6 = G6R;
            }

            if (P <= PK) {
                // rarefaction wave
                PRAT = P / PK;
                F = G4 * CK * (Math.Pow(PRAT, G1) - 1.0);
                FD = (1.0 / (DK * CK)) * Math.Pow(PRAT, -G2);
            } else {
                // shock wave
                AK = G5 / DK;
                BK = G6 * PK;
                QRT = Math.Sqrt(AK / (BK + P));
                F = (P - PK) * QRT;
                FD = (1.0 - 0.5 * (P - PK) / (BK + P)) * QRT;
            }
        }

        /// <summary>
        /// Samples the solution throughout the wave pattern for given pressure
        /// PM and velocity UM in the Star Region. Sampling is performed in
        /// terms of the 'speed' S=X/T. Sampled values are D,U,P
        /// </summary>
        private void Sample(double PM, double UM, double S, out double D, out double U, out double P) {
            double C, CML, STL, PML, SL, SHL, SHR, CMR, STR, PMR, SR;
            if (S <= UM) {
                // Sampling point lies to the left of the contact discontinuity

                if (PM <= pressureLeft) {
                    // left rarefaction
                    SHL = normalVelocityLeft - speedOfSoundLeft;

                    if (S <= SHL) {

                        //sampled point is left data state
                        D = densityLeft;
                        U = normalVelocityLeft;
                        P = pressureLeft;
                    } else {
                        CML = speedOfSoundLeft * Math.Pow(PM / pressureLeft, G1L);
                        STL = UM - CML;

                        if (S > STL) {
                            // sampled point is Star Left state
                            D = densityLeft * Math.Pow(PM / pressureLeft, 1.0 / kappaLeft);
                            U = UM;
                            P = PM;
                        } else {
                            // sampled point is inside left fan
                            U = G5L * (speedOfSoundLeft + G7L * normalVelocityLeft + S); // (UL + S), mistake by Toro?
                            C = G5L * (speedOfSoundLeft + G7L * (normalVelocityLeft - S));
                            D = densityLeft * Math.Pow(C / speedOfSoundLeft, G4L);
                            P = pressureLeft * Math.Pow(C / speedOfSoundLeft, G3L);
                        }
                    }
                } else {
                    //left shock
                    PML = PM / pressureLeft;
                    SL = normalVelocityLeft - speedOfSoundLeft * Math.Sqrt(G2L * PML + G1L);

                    if (S <= SL) {
                        // sampled point is left data state
                        D = densityLeft;
                        U = normalVelocityLeft;
                        P = pressureLeft;

                    } else {
                        // sampled point is LStar Left state
                        D = densityLeft * (PML + G6L) / (PML * G6L + 1.0);
                        U = UM;
                        P = PM;
                    }
                }
            } else {
                // sampling point lies to the right of the contact discontinuity
                if (PM > pressureRight) {
                    // right shock
                    PMR = PM / pressureRight;
                    SR = normalVelocityRight + speedOfSoundRight * Math.Sqrt(G2R * PMR + G1R);

                    if (S >= SR) {
                        // sampled point is right data state
                        D = densityRight;
                        U = normalVelocityRight;
                        P = pressureRight;
                    } else {
                        // sampled point is Star Right state
                        D = densityRight * (PMR + G6R) / (PMR * G6R + 1.0);
                        U = UM;
                        P = PM;
                    }
                } else {
                    // right rarefaction
                    SHR = normalVelocityRight + speedOfSoundRight;

                    if (S >= SHR) {
                        // sampled point is right data state
                        D = densityRight;
                        U = normalVelocityRight;
                        P = pressureRight;
                    } else {
                        CMR = speedOfSoundRight * Math.Pow(PM / pressureRight, G1R);
                        STR = UM + CMR;

                        if (S <= STR) {
                            // sampled point is Star Right state
                            D = densityRight * Math.Pow(PM / pressureRight, 1.0 / kappaRight);
                            U = UM;
                            P = PM;
                        } else {
                            // sampled point is inside left fan
                            U = G5R * (-speedOfSoundRight + G7R * normalVelocityRight + S);
                            C = G5R * (speedOfSoundRight - G7R * (normalVelocityRight - S));
                            D = densityRight * Math.Pow(C / speedOfSoundRight, G4R);
                            P = pressureRight * Math.Pow((C / speedOfSoundRight), G3R);
                        }
                    }
                }
            }
        }
    }
}
