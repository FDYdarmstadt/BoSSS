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

using System;
using System.Collections.Generic;
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Timestepping;
using BoSSS.Solution.Utils;
using BoSSS.Foundation.Grid;
using ilPSP.Tracing;
using System.Globalization;
using BoSSS.Foundation.Quadrature;
using ilPSP.Utils;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Solution.LevelSetTools.Reinit.Iterative {
    /// <summary>
    /// A class for obtaining the time integrated spatial discretization of the hyperbolic level set reinitialization equation
    /// See Dissertation of Roozbeh Mousavi 20124
    /// </summary>
    public class LevelSetReInitialization : IReinitialisationAlgorithm {
        /// <summary>
        /// Initialize Internal Fields
        /// </summary>
        /// <param name="LS"></param>
        public LevelSetReInitialization(LevelSet LS) {
            m_ctx = (GridData)(LS.Basis.GridDat);
            int D = m_ctx.SpatialDimension;

            LSUG = new VectorField<SinglePhaseField>(D, LS.Basis, "LSUG", SinglePhaseField.Factory);
            LSDG = new VectorField<SinglePhaseField>(D, LS.Basis, "LSDG", SinglePhaseField.Factory);
            LSCG = new VectorField<SinglePhaseField>(D, LS.Basis, "LSCG", SinglePhaseField.Factory);
            LSCGV = new SinglePhaseField(LS.Basis, "LSCGV");
        }

        /// <summary>
        /// The internal context of the class which is connected to the outside via the constructor
        /// </summary>
        GridData m_ctx;

        /// <summary>
        /// Level set gradient using upwind flux based on the direction of the edge or face normal vector
        /// </summary>
        VectorField<SinglePhaseField> LSUG;

        /// <summary>
        /// Level set gradient using downwind flux based on the direction of the edge or face normal vector
        /// </summary>
        VectorField<SinglePhaseField> LSDG;

        /// <summary>
        /// Level set gradient using central flux
        /// </summary>
        VectorField<SinglePhaseField> LSCG;

        /// <summary>
        /// Value (magnitude) of the level set gradient using central flux
        /// </summary>
        SinglePhaseField LSCGV;

        RungeKuttaScheme m_Scheme = RungeKuttaScheme.TVD3;

        /// <summary>
        /// Runge-Kutta method that is used for Pseudo-timestepping
        /// </summary>
        public RungeKuttaScheme Scheme {
            get {
                return m_Scheme;
            }
            set {
                m_Scheme = value;
            }
        }


        /// <summary>
        /// Instantiating the class corresponding to the specified flux calculation method
        /// </summary>
        /// <param name="ctx"> The context </param>
        /// <param name="f"> Specifying the method of flux calculation </param>
        /// <param name="i"> Component of the flux vector </param>
        /// <returns></returns>
        GradientFlux CreateFlux(GridData ctx, string f, int i) {
            switch (f) {
                case "Upwind":
                    return new UpwindGradientFlux(ctx, i);
                case "Downwind":
                    return new DownwindGradientFlux(ctx, i);
                case "Central":
                    return new CentralGradientFlux(ctx, i);
                default:
                    throw new ApplicationException();
            }
        }

        /// <summary>
        /// Calculating the level set gradient using the specified scheme in a narrow band around the zero level set, therein the calculions are performed
        /// </summary>
        /// <param name="LS"> The level set function </param>
        /// <param name="LSG"> Gradient of the level set function </param>
        /// <param name="f"> Specifying the method of flux calculation </param>
        /// <param name="Restriction"> The narrow band around the zero level set wherein the calculations are performed </param>
        void CalculateLevelSetGradient(LevelSet LS, VectorField<SinglePhaseField> LSG, string f, SubGrid Restriction) {
            SpatialOperator SO;
            CoordinateMapping CoDom;

            if (m_ctx.SpatialDimension == 2) {
                SO = new SpatialOperator(1, 2,QuadOrderFunc.Linear(), new string[] { "LS", "LSG[0]", "LSG[1]" });
                SO.EquationComponents["LSG[0]"].Add(CreateFlux(m_ctx, f, 0));
                SO.EquationComponents["LSG[1]"].Add(CreateFlux(m_ctx, f, 1));
                SO.Commit();
                CoDom = new CoordinateMapping(LSG[0], LSG[1]);
            } else if (m_ctx.SpatialDimension == 3) {
                SO = new SpatialOperator(1, 3,QuadOrderFunc.Linear(), new string[] { "LS", "LSG[0]", "LSG[1]", "LSG[2]" });
                SO.EquationComponents["LSG[0]"].Add(CreateFlux(m_ctx, f, 0));
                SO.EquationComponents["LSG[1]"].Add(CreateFlux(m_ctx, f, 1));
                SO.EquationComponents["LSG[2]"].Add(CreateFlux(m_ctx, f, 2));
                SO.Commit();
                CoDom = new CoordinateMapping(LSG[0], LSG[1], LSG[2]);
            } else {
                throw new NotSupportedException();
            }

            SO.Evaluate(1.0, 0.0, LS.Mapping, null, CoDom, sgrd: Restriction, bndMode: SpatialOperator.SubGridBoundaryModes.OpenBoundary);
        }


        /// <summary>
        /// This only exists to implement the standard Interface
        /// Do ReInit using 5000 (sic!) Timesteps
        /// </summary>
        /// <param name="LS"></param>
        /// <param name="Restriction"></param>
        public void ReInitialize(LevelSet LS, SubGrid Restriction = null) {
            Console.WriteLine("You are running the hyperbolic ReInit with 5000 Timesteps! Are you sure?");
            GridData grd = (GridData)(LS.GridDat);
            double k = LS.Basis.Degree;
            double h = grd.Cells.h_minGlobal;
            double TimestepSize = 0.5 * h / (k * k);  // CFL criterion, 0.5 is the CFL constant
            double thickness = h;

            this.ReInitialize(LS, Restriction, thickness, TimestepSize, 5000);

        }

        /// <summary>
        /// Do ReInit
        /// Configurable
        /// </summary>
        /// <param name="LS"></param>
        /// <param name="NumberOfTimeSteps"></param>
        /// <param name="Restriction"></param>
        public void ReInitialize(LevelSet LS,int NumberOfTimeSteps = 5000, SubGrid Restriction = null) {
            var grd = LS.GridDat;
            
            double k = LS.Basis.Degree;
            double h = ((GridData)grd).Cells.h_minGlobal;
            double TimestepSize = 0.5*h/(k*k);  // CFL criterion, 0.5 is the CFL constant
            double thickness = h;

            this.ReInitialize(LS, Restriction, thickness, TimestepSize, NumberOfTimeSteps);

        }



        /// <summary>
        /// Obtaining the time integrated spatial discretization of the reinitialization equation in a narrow band around the zero level set, based on a Godunov's numerical Hamiltonian calculation
        /// </summary>
        /// <param name="LS"> The level set function </param>
        /// <param name="Restriction"> The narrow band around the zero level set </param>
        /// <param name="NumberOfTimesteps">
        /// maximum number of pseudo-timesteps 
        /// </param>
        /// <param name="thickness">
        /// The smoothing width of the signum function.
        /// This is the main stabilization parameter for re-initialization. 
        /// It should be set to approximately 3 cells.
        /// </param>
        /// <param name="TimestepSize">
        /// size of the pseudo-timestep
        /// </param>
        public void ReInitialize(LevelSet LS, SubGrid Restriction, double thickness, double TimestepSize, int NumberOfTimesteps) {
            using (var tr = new FuncTrace()) {
                
                // log parameters:
                tr.Info("thickness: " + thickness.ToString(NumberFormatInfo.InvariantInfo));
                tr.Info("TimestepSize: " + TimestepSize.ToString(NumberFormatInfo.InvariantInfo));
                tr.Info("NumberOfTimesteps: " + NumberOfTimesteps);

                ExplicitEuler TimeIntegrator;

                SpatialOperator SO;
                Func<int[], int[], int[], int> QuadratureOrder = QuadOrderFunc.NonLinear(3);
                if (m_ctx.SpatialDimension == 2) {
                    SO = new SpatialOperator(1, 5, 1, QuadratureOrder, new string[] { "LS", "LSCGV", "LSDG[0]", "LSUG[0]", "LSDG[1]", "LSUG[1]", "Result" });
                    SO.EquationComponents["Result"].Add(new GodunovHamiltonian(m_ctx, thickness));
                    SO.Commit();
                    TimeIntegrator = new RungeKutta(m_Scheme, SO, new CoordinateMapping(LS), new CoordinateMapping(LSCGV, LSDG[0], LSUG[0], LSDG[1], LSUG[1]), sgrd:Restriction);
                } else {
                    SO = new SpatialOperator(1, 7, 1, QuadratureOrder, new string[] { "LS", "LSCGV", "LSDG[0]", "LSUG[0]", "LSDG[1]", "LSUG[1]", "LSDG[2]", "LSUG[2]", "Result" });
                    SO.EquationComponents["Result"].Add(new GodunovHamiltonian(m_ctx, thickness));
                    SO.Commit();
                    TimeIntegrator = new RungeKutta(m_Scheme, SO, new CoordinateMapping(LS), new CoordinateMapping(LSCGV, LSDG[0], LSUG[0], LSDG[1], LSUG[1], LSDG[2], LSUG[2]), sgrd:Restriction);
                }



                // Calculating the gradients in each sub-stage of a Runge-Kutta integration procedure
                ExplicitEuler.ChangeRateCallback EvalGradients = delegate(double t1, double t2) {
                    LSUG.Clear();
                    CalculateLevelSetGradient(LS, LSUG, "Upwind", Restriction);

                    LSDG.Clear();
                    CalculateLevelSetGradient(LS, LSDG, "Downwind", Restriction);

                    LSCG.Clear();
                    CalculateLevelSetGradient(LS, LSCG, "Central", Restriction);

                    LSCGV.Clear();
                    var VolMask = (Restriction != null) ? Restriction.VolumeMask : null;
                    LSCGV.ProjectAbs(1.0, VolMask, LSCG.ToArray());
                };
                TimeIntegrator.OnBeforeComputeChangeRate += EvalGradients;

                
                {
                    EvalGradients(0, 0);
                    var GodunovResi = new SinglePhaseField(LS.Basis, "Residual");
                    SO.Evaluate(1.0, 0.0, LS.Mapping, TimeIntegrator.ParameterMapping.Fields, GodunovResi.Mapping, Restriction);
                    
                    //Tecplot.Tecplot.PlotFields(ArrayTools.Cat<DGField>( LSUG, LSDG, LS, GodunovResi), "Residual", 0, 3); 
                }
                

                

                // pseudo-timestepping
                // ===================
                double factor = 1.0;
                double time = 0;
                LevelSet prevLevSet = new LevelSet(LS.Basis, "prevLevSet");

                CellMask RestrictionMask = (Restriction == null) ? null : Restriction.VolumeMask;

                for (int i = 0; (i < NumberOfTimesteps); i++) {
                    tr.Info("Level set reinitialization pseudo-timestepping, timestep " + i);

                    // backup old Levelset
                    // -------------------
                    prevLevSet.Clear();
                    prevLevSet.Acc(1.0, LS, RestrictionMask);

                    // time integration
                    // ----------------
                    double dt = TimestepSize * factor;
                    tr.Info("dt = " + dt.ToString(NumberFormatInfo.InvariantInfo) + " (factor = " + factor.ToString(NumberFormatInfo.InvariantInfo) + ")");
                    TimeIntegrator.Perform(dt);
                    time += dt;

                    // change norm
                    // ------

                    prevLevSet.Acc(-1.0, LS, RestrictionMask);
                    double ChangeNorm = prevLevSet.L2Norm(RestrictionMask);
                    Console.WriteLine("Reinit: PseudoTime: {0}  - Changenorm: {1}", i, ChangeNorm);

                    //Tecplot.Tecplot.PlotFields(new SinglePhaseField[] { LS }, m_ctx, "Reinit-" + i, "Reinit-" + i, i, 3); 
                }

                //*/
            }
        }

        /// <summary>
        /// some 
        /// </summary>
        public void One(MultidimensionalArray I, MultidimensionalArray O) {
            double x;
            double y;

            for (int i = 0; i < I.GetLength(0); i++) {
                x = I[i, 0];
                y = I[i, 1];

                O[i] = 1.0;
            }

        }
    }

    /// <summary>
    /// Calculation of the level set gradient using upwind (taking the value on the left side, resp., closer to negative infinity) 
    /// flux based on the direction of the face normal vector
    /// </summary>
    class UpwindGradientFlux : GradientFlux {
        public UpwindGradientFlux(GridData ctx, int i)
            : base(ctx, i) {
        }

        protected override double InnerEdgeFlux(double time, double[] x, double[] normal, double[] Uin, double[] Uout, int jEdge) {
            double n_i = normal[m_i];
            return Math.Max(n_i, 0.0) * Uin[0] + Math.Min(n_i, 0.0) * Uout[0];
        }
    }

    /// <summary>
    /// Calculation of the level set gradient using downwind flux (taking the value on the right side, resp., closer to positive infinity) 
    /// based on the direction of the face normal vector
    /// </summary>
    class DownwindGradientFlux : GradientFlux {
        public DownwindGradientFlux(GridData ctx, int i)
            : base(ctx, i) {
        }

        protected override double InnerEdgeFlux(double time, double[] x, double[] normal, double[] Uin, double[] Uout, int jEdge) {
            double n_i = normal[m_i];
            return Math.Max(n_i, 0.0) * Uout[0] + Math.Min(n_i, 0.0) * Uin[0];
        }
    }

    /// <summary>
    /// Calculation of the level set gradient using central flux based on the direction of the face normal vector
    /// </summary>
    class CentralGradientFlux : GradientFlux {
        public CentralGradientFlux(GridData ctx, int i)
            : base(ctx, i) {
        }

        protected override double InnerEdgeFlux(double time, double[] x, double[] normal, double[] Uin, double[] Uout, int jEdge) {
            double n_i = normal[m_i];
            return n_i * 0.5 * (Uout[0] + Uin[0]);
        }
    }

    /// <summary>
    /// Calculation of the numerical fluxes included in the spatial discretization of the level set gradients
    /// </summary>
    abstract class GradientFlux : NonlinearFlux {
        public GradientFlux(GridData ctx, int i) {
            m_ctx = ctx;
            D = m_ctx.SpatialDimension;
            m_i = i;
        }

        /// <summary>
        /// Internal context of the class which is connected to the outside via the constructor
        /// </summary>
        GridData m_ctx;

        /// <summary>
        /// Spatial dimension
        /// </summary>
        protected int D;

        /// <summary>
        /// Spatial component
        /// </summary>
        protected int m_i;

        string[] m_ArgumentOrdering = new string[] { "LS" };

        protected override double BorderEdgeFlux(double time, double[] x, double[] normal, byte EdgeTag, double[] Uin, int jEdge) {
            return normal[m_i] * Uin[0];
        }

        protected override void Flux(double time, double[] x, double[] U, double[] output) {
            for (int d = D - 1; d >= 0; d--)
                output[d] = 0;
            output[m_i] = U[0];
        }

        public override IList<string> ArgumentOrdering {
            get {
                return m_ArgumentOrdering;
            }
        }
    }

    /// <summary>
    /// Calculation of the Godunov's numerical Hamiltonian as a source term made by a combination of the level set gradients
    /// </summary>
    class GodunovHamiltonian : NonlinearSource {
        public GodunovHamiltonian(GridData ctx, double _thickness_RI) {
            m_ctx = ctx;
            thickness_RI = _thickness_RI;
        }

        GridData m_ctx;
        double thickness_RI;

        string[] m_AO = new string[] { "LS" };

        string[] m_PO2D = new string[] { "LSCGV", "LSDG[0]", "LSUG[0]", "LSDG[1]", "LSUG[1]" };

        string[] m_PO3D = new string[] { "LSCGV", "LSDG[0]", "LSUG[0]", "LSDG[1]", "LSUG[1]", "LSDG[2]", "LSUG[2]" };

        // A smoothed sign function
        double Sign(double LS, double LSCGV, double dx, double[] X) {
            // Narrow band width
            double w = 2.0 * dx;
            return LS / Math.Sqrt(Square(LS) + Square(Math.Abs(LSCGV) * (w / 2.0)));
            //double x = X[0], y = X[1];
            //double phi = 2.2 - (x*x/0.9 + y*y/1.1);
            //return phi / Math.Sqrt(Square(phi) + (w / 2.0));
        }

        // Limiting the level set slop to be one inside a band around the interface, and zero outside the band
        double Limit(double LS, double dx) {
            #region Limiting the level set function outside a narrow band: The commented lines bellow should be uncommented
            /*
            // Narrow band width
            double w = 2.0 * dx;

            // Normalized level set function
            double xi = Math.Abs(LS) / (w / 2.0);

            if (xi <= 1.0) {
                return 1.0;
            } else if (xi >= 1.25) {
                return 0.0;
            } else {
                // A smoothed transition from inside to outside the band
                return -4 * (xi - 1.25) - Math.Sin(4 * Math.PI * (xi - 1.25)) / Math.PI;
            }
            */
            #endregion

            return 1.0;
        }

        protected override double Source(double time, int j, double[] X, double[] U) {
            // Maximum local grid resolution
            //double dx = m_ctx.GridDat.h_max[j];
            double dx = thickness_RI;



            if (m_ctx.SpatialDimension == 2) {
                if (Sign(U[0], U[1], dx, X) >= 0.0) {
                    return Sign(U[0], U[1], dx, X)*(Math.Sqrt(Math.Max(Square(Math.Min(U[2], 0.0))
                                                                     , Square(Math.Max(U[3], 0.0)))
                                                            + Math.Max(Square(Math.Min(U[4], 0.0))
                                                                     , Square(Math.Max(U[5], 0.0)))) - Limit(U[0], dx));
                } else {
                    return Sign(U[0], U[1], dx, X)*(Math.Sqrt(Math.Max(Square(Math.Max(U[2], 0.0))
                                                                     , Square(Math.Min(U[3], 0.0)))
                                                            + Math.Max(Square(Math.Max(U[4], 0.0))
                                                                     , Square(Math.Min(U[5], 0.0)))) - Limit(U[0], dx));
                }
            } else {
                if (Sign(U[0], U[1], dx, X) >= 0.0) {
                    return Sign(U[0], U[1], dx, X)*(Math.Sqrt(Math.Max(Square(Math.Min(U[2], 0.0))
                                                                     , Square(Math.Max(U[3], 0.0)))
                                                            + Math.Max(Square(Math.Min(U[4], 0.0))
                                                                     , Square(Math.Max(U[5], 0.0)))
                                                            + Math.Max(Square(Math.Min(U[6], 0.0))
                                                                     , Square(Math.Max(U[7], 0.0)))) - Limit(U[0], dx));
                } else {
                    return Sign(U[0], U[1], dx, X)*(Math.Sqrt(Math.Max(Square(Math.Max(U[2], 0.0))
                                                                     , Square(Math.Min(U[3], 0.0)))
                                                            + Math.Max(Square(Math.Max(U[4], 0.0))
                                                                     , Square(Math.Min(U[5], 0.0)))
                                                            + Math.Max(Square(Math.Max(U[6], 0.0))
                                                                     , Square(Math.Min(U[7], 0.0)))) - Limit(U[0], dx));
                }
            }
        }

        static double Square(double x) {
            return x * x;
        }

        public override IList<string> ArgumentOrdering {
            get {
                return m_AO;
            }
        }

        public override IList<string> ParameterOrdering {
            get {
                if (m_ctx.SpatialDimension == 2) {
                    return m_PO2D;
                } else {
                    return m_PO3D;
                }
            }
        }
    }
}