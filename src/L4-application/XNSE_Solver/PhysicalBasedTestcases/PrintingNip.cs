using System;
using System.Collections.Generic;
using System.Linq;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.Grid.RefElements;
using BoSSS.Foundation.IO;
using BoSSS.Solution.Control;
using BoSSS.Foundation;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.EnergyCommon;
using BoSSS.Solution.XNSECommon;
using MathNet.Numerics;
using MathNet.Numerics.Interpolation;
using MathNet.Numerics.RootFinding;




/// <summary>
/// This namespace contains various helper routines for the printing nip simulations for SFB1194 K65
/// </summary>
namespace BoSSS.Application.XNSE_Solver.PhysicalBasedTestcases.PrintingNip {

    #region Boundary Values
    /// <summary>
    /// Setting the boundary values
    /// </summary>
    public static class BoundaryValueFactory {
        public static string GetPrefixCode(double _delta, double RotVel, double Radius = 0.1) {
            using (var stw = new System.IO.StringWriter()) {

                stw.WriteLine("static class BoundaryValues {");
                stw.WriteLine("  static public Vector VelVec(Vector X) {");
                stw.WriteLine("    var Org = new Vector(0, Math.Sign(X[1]) * " + (Radius + _delta / 2.0) + ");");
                stw.WriteLine("    var Dir = Org - X;");
                stw.WriteLine("    var R   = Dir.Rotate2D(-Math.Sign(X[1]) * Math.PI/2);");
                stw.WriteLine("    R.NormalizeInPlace();");
                stw.WriteLine("    R = R*" + RotVel + ";");
                stw.WriteLine("    return R;");
                stw.WriteLine("  }");
                stw.WriteLine(" ");
                stw.WriteLine("  static public double VelX(double[] X) {");
                stw.WriteLine("    return VelVec(X).x;");
                stw.WriteLine("  }");
                stw.WriteLine("  static public double VelY(double[] X) {");
                stw.WriteLine("    return VelVec(X).y;");
                stw.WriteLine("  }");
                stw.WriteLine("}");

                return stw.ToString();
            }
        }

        static public Formula Get_VelX(double Radius, double _delta, double RotVel) {
            return new Formula("BoundaryValues.VelX", AdditionalPrefixCode: GetPrefixCode(Radius, _delta, RotVel));
        }

        static public Formula Get_VelY(double Radius, double _delta, double RotVel) {
            return new Formula("BoundaryValues.VelY", AdditionalPrefixCode: GetPrefixCode(Radius, _delta, RotVel));
        }
    }
    #endregion

    #region Grid
    public static class GridFactory {

        #region StaticVariables

        // Set myDb from Notebook, if grid saving is intended
        public static IDatabaseInfo myDb = null;
        public static CellType cellType = CellType.Square_49;
        #endregion


        #region GridCreation

        // create XNodes and shift them for an appropriate grid resolution
        public static double[] GetXNodes(int Res, double xi0, double xi1, double c0, double c1) {
            double[] xNodes = GenericBlas.Linspace(xi0, xi1, Res * 20 + 1);
            xNodes = xNodes.Select(x => 2 * MathNet.Numerics.Trig.Acot(c0 / c1 * MathNet.Numerics.Trig.Cot(x / 2))).ToArray();
            xNodes = xNodes.Select(x => x < 0.0 ? x + 2 * Math.PI : x).ToArray();
            return xNodes;
        }

        static Dictionary<Tuple<int, double>, Grid2D> GridCache = new Dictionary<Tuple<int, double>, Grid2D>();

        public static Grid2D GenerateGrid(int Res, double _delta, double R = 0.1) {
            var id = Tuple.Create(Res, _delta);
            if (GridCache.ContainsKey(id))
                return GridCache[id];

            var Trf = new GridTrafo(_delta, R);

            double c = Trf.c;
            //  var xNodes = GenericBlas.Linspace(-Trf.R, Trf.R, 20*Res + 1);
            var xNodes = GenericBlas.Linspace(-1, 1, 20 * Res + 1);
            xNodes = xNodes.Select(x => Math.Sign(x) * x * x * Trf.R).ToArray();
            var xiNodes = xNodes.Select(x => Math.Atan(x / Trf.c) - Math.Atan(x / -Trf.c) + Math.PI).ToArray();
            var etaNodes = GenericBlas.Linspace(Trf.eta0, Trf.eta1, Res + 1);

            var grd = Grid2D.Cartesian2DGrid(xiNodes, etaNodes, periodicX: false, NonlinearGridTrafo: Trf.T, type: cellType);

            grd.DefineEdgeTags(delegate (Vector X) {
                if (Math.Abs(Math.Sqrt(X[0] * X[0] + Math.Pow(-Trf.epsilon - Trf.R - X[1], 2)) - Trf.R) < 1e-12)
                    return "wall_substrat";
                else if (Math.Abs(Math.Sqrt(X[0] * X[0] + Math.Pow(+Trf.epsilon + Trf.R - X[1], 2)) - Trf.R) < 1e-12)
                    return "wall_walze";
                else if (X[0] < 0.0)
                    return "pressure_outlet_in";
                else //(X[0] > 0.0)
                    return "pressure_outlet_out";
            });

            if (myDb != null) {
                IGrid grid = grd;
                myDb.Controller.DBDriver.SaveGridIfUnique(ref grid, out _, myDb);// this takes forever
            }
            GridCache[id] = grd;
            return grd;
        }

        #endregion

        #region GridTransformation 
        // Grid transformation 
        class GridTrafo {
            private double m_Radius = 0.1; // assume unit Meters
            private double m_epsilon; // nip width between cylinder and plate
            private double m_c; // position of bipolar focal points
            private double m_eta0, m_eta1;
            private double m_xi0, m_xi1;

            public GridTrafo(double epsilon, double R = 0.1) {
                this.m_epsilon = epsilon / 2.0; // half nip width on both sides of symmetry plane
                this.m_Radius = R;
                this.m_c = m_Radius * Math.Abs(MathNet.Numerics.Trig.Sinh(MathNet.Numerics.Trig.Acosh(m_epsilon / m_Radius + 1)));
                this.m_eta0 = -0.5 * Math.Log((Math.Pow(m_epsilon + m_c, 2)) / (Math.Pow(m_epsilon - m_c, 2))); ;
                this.m_eta1 = 0.5 * Math.Log((Math.Pow(m_epsilon + m_c, 2)) / (Math.Pow(m_epsilon - m_c, 2)));
                this.m_xi0 = MathNet.Numerics.Trig.Atan(-m_Radius / m_c) - MathNet.Numerics.Trig.Atan(-m_Radius / -m_c);
                this.m_xi1 = -this.m_xi0;
            }

            public double R {
                get {
                    double val = m_Radius;
                    return val;
                }
            }
            public double c {
                get {
                    double val = m_c;
                    return val;
                }
            }
            public double epsilon {
                get {
                    double val = m_epsilon;
                    return val;
                }
            }
            public double eta0 {
                get {
                    double val = m_eta0;
                    return val;
                }
            }
            public double eta1 {
                get {
                    double val = m_eta1;
                    return val;
                }
            }

            public double xi0 {
                get {
                    double val = m_xi0;
                    return val + Math.PI; // need to shift this due, to calculation of Atan
                }
            }
            public double xi1 {
                get {
                    double val = m_xi1;
                    return val + Math.PI; // need to shift this due, to calculation of Atan
                }
            }

            public Vector T(Vector Q) {

                double xi = Q.x;
                double nu = Q.y;

                var R = new Vector(2);
                R.x = -c * Math.Sin(xi) / (MathNet.Numerics.Trig.Cosh(nu) - Math.Cos(xi));
                R.y = c * MathNet.Numerics.Trig.Sinh(nu) / (MathNet.Numerics.Trig.Cosh(nu) - Math.Cos(xi));
                return R;
            }
        }
        #endregion

    }
    #endregion

    #region Postprocessing

        public static class Postprocessing {

            // Helper Function to accelerate postprocessing
            static Dictionary<Tuple<double, double, double, double>, IEnumerable<DGField>> FieldCache = new Dictionary<Tuple<double, double, double, double>, IEnumerable<DGField>>();

            static IEnumerable<DGField> LoadFields(ISessionInfo si) {
                Tuple<double, double, double, double> key = Tuple.Create(Convert.ToDouble(si.KeysAndQueries["id:Radius"]), Convert.ToDouble(si.KeysAndQueries["id:delta"]), Convert.ToDouble(si.KeysAndQueries["id:V_Wall"]), Convert.ToDouble(si.KeysAndQueries["id:P_Diff"]));

                if (FieldCache.TryGetValue(key, out var fields)) {
                    return fields;
                } else {
                    FieldCache[key] = si.Timesteps.Last().Fields;
                    Console.WriteLine("Loaded Fields for {0}", si.Name);
                    return FieldCache[key];
                }
            }

            static int InterpolationDegree(IGridData grd) {
                return ((GridData)grd).Cells.GetInterpolationDegree(0);
            }


            // Helper function to create splines on edges
            static public CubicSpline SplineOnEdge(EdgeMask em, DGField field, int d, out double lower_Bound, out double upper_Bound, double offset = 0.0) {
                var grd = field.GridDat;

                List<double> nodes = new List<double>();
                List<double> values = new List<double>();

                int TransformDegree = grd.SpatialDimension * (InterpolationDegree(grd) - 1);

                EdgeQuadrature.GetQuadrature(new int[] { 1 }, grd,
                    (new EdgeQuadratureScheme(true, em)).Compile(grd, field.Basis.Degree * 2 + TransformDegree),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {

                        MultidimensionalArray DummyIN = MultidimensionalArray.Create(Length, QR.NoOfNodes);
                        MultidimensionalArray DummyOT = MultidimensionalArray.Create(Length, QR.NoOfNodes);

                        MultidimensionalArray GlobalNodes = MultidimensionalArray.Create(Length, QR.NoOfNodes, 2);

                        if (field is XDGField xField) {
                            xField.GetSpeciesShadowField("A").EvaluateEdge(i0, Length, QR.Nodes, DummyIN, DummyOT, null, null, null, null, 0, 0.0);
                        } else {
                            field.EvaluateEdge(i0, Length, QR.Nodes, DummyIN, DummyOT, null, null, null, null, 0, 0.0);
                        }

                        for (int i = 0; i < Length; i++) {
                            int iTrafo = ((GridData)grd).Edges.Edge2CellTrafoIndex[i0 + i, 0];
                            NodeSet volNodeSet = QR.Nodes.GetVolumeNodeSet(grd, iTrafo, false);
                            int jCell = ((GridData)grd).Edges.CellIndices[i0 + i, 0];
                            grd.TransformLocal2Global(volNodeSet, jCell, 1, GlobalNodes, i);
                            int K = QR.NoOfNodes;
                            for (int k = 0; k < K; k++) {
                                nodes.Add(GlobalNodes[i, k, d]);
                                values.Add(DummyIN[i, k] - offset);
                            }
                        }

                    }, delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                    }).Execute();

                lower_Bound = nodes.Min();
                upper_Bound = nodes.Max();

                return CubicSpline.InterpolateAkima(nodes.ToArray(), values.ToArray());
            }

            // Finds the position where the specified field assumes the given value
            static public double PosOfValueOnEdge(double value, EdgeMask em, DGField field, int d = 0) {

                var Spline = SplineOnEdge(em, field, 0, out double lB, out double uB, value);
                double root = Bisection.FindRoot(t => Spline.Interpolate(t), lB, uB);

                return root;
            }

            // Finds the position where the specified field assumes an extremum, in the given direction
            static public double PosOfExtremumOnEdge(EdgeMask em, DGField field, int d) {

                var dField = field.CloneAs();
                dField.Clear();
                dField.Derivative(1.0, field, d);

                double root = PosOfValueOnEdge(0.0, em, dField, d);

                return root;
            }

            // only look at 3 cm before and past the nip center
            // to search for a pressure maximum
            static public bool NipSelectionFunction(ilPSP.Vector X) {
                double x = X[0];
                if (double.IsNaN(x))
                    throw new ArithmeticException();
                if (Math.Abs(x) < 0.03)
                    return true;
                return false;
            }

            static public object PressureRange(ISessionInfo si) {
                var fields = LoadFields(si);
                DGField Pressure = fields.Single(f => f.Identification == "Pressure");
                CellMask NipRegionMask = CellMask.GetCellMask(Pressure.GridDat as GridData, NipSelectionFunction);
                //Console.WriteLine("n = " + NipRegionMask.NoOfItemsLocally);

                double mini, maxi;
                long jMini, jMaxi;
                Pressure.GetExtremalValues(out mini, out maxi, out jMini, out jMaxi, NipRegionMask);
                double range = maxi - mini;
                //Console.WriteLine("mini = " + mini);
                //Console.WriteLine("maxi = " + maxi);
                //Console.WriteLine("PressureRange for " + si.Name + " done!");

                object ret = range;
                return ret;
            }

            // Maximum Velocity in nip vicinity
            static public object VelocityXMax(ISessionInfo si) {
                var fields = LoadFields(si);
                DGField VelocityX = fields.Single(f => f.Identification == "VelocityX");
                CellMask NipRegionMask = CellMask.GetCellMask(VelocityX.GridDat as GridData, NipSelectionFunction);
                //Console.WriteLine("n = " + NipRegionMask.NoOfItemsLocally);

                double mini, maxi;
                long jMini, jMaxi;
                VelocityX.GetExtremalValues(out mini, out maxi, out jMini, out jMaxi, NipRegionMask);
                double value = maxi;
                //Console.WriteLine("mini = " + mini);
                //Console.WriteLine("maxi = " + maxi);
                //Console.WriteLine("VelocityMax for " + si.Name + " done!");

                object ret = value;
                return ret;
            }

            // massflux through nip
            static public object Massflux(ISessionInfo si) {
                var fields = LoadFields(si);
                DGField VelocityX = fields.Single(f => f.Identification == "VelocityX");
                var grd = VelocityX.GridDat;

                double rho = Convert.ToDouble(si.KeysAndQueries["PhysicalParameters.rho_A"]);
                double delta = Convert.ToDouble(si.KeysAndQueries["id:delta"]);

                EdgeMask em = new EdgeMask(grd, X => Math.Abs(X[0]) < 1e-10); // edges need to be aligned in the nip!
                int TransformDegree = grd.SpatialDimension * (InterpolationDegree(grd) - 1);

                double massflux = 0.0;
                EdgeQuadrature.GetQuadrature(new int[] { 1 }, grd,
                    (new EdgeQuadratureScheme(true, em)).Compile(grd, VelocityX.Basis.Degree * 2 + TransformDegree),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {

                        MultidimensionalArray VelocityXIN = MultidimensionalArray.Create(Length, QR.NoOfNodes, 1);
                        MultidimensionalArray VelocityXOT = MultidimensionalArray.Create(Length, QR.NoOfNodes, 1);

                        ((XDGField)VelocityX).GetSpeciesShadowField("A").EvaluateEdge(i0, Length, QR.Nodes, VelocityXIN.ExtractSubArrayShallow(-1, -1, 0), VelocityXOT.ExtractSubArrayShallow(-1, -1, 0), null, null, null, null, 0, 0.0);

                        EvalResult.Clear();
                        EvalResult.Acc(0.5, VelocityXIN);
                        EvalResult.Acc(0.5, VelocityXOT);
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for (int _i = 0; _i < Length; _i++) {
                            massflux += rho * ResultsOfIntegration[_i, 0];
                        }
                    }
                ).Execute();
                //Console.WriteLine("Massflux for " + si.Name + " done!");

                object ret = massflux;
                return ret;
            }

            #region Shear Rates on cylinder, substrate and symmetry plane in the nip

            static public object NipShearRateSubstrate(ISessionInfo si) {
                var fields = LoadFields(si);
                DGField VelocityX = fields.Single(f => f.Identification == "VelocityX");
                DGField VelocityY = fields.Single(f => f.Identification == "VelocityY");
                var grd = (GridData)VelocityX.GridDat;

                var ShearRate = (DGField)VelocityX.Clone();
                ShearRate.Clear();
                ShearRate.Derivative(1.0, VelocityX, 1);
                ShearRate.Derivative(1.0, VelocityY, 0);

                double R = Convert.ToDouble(si.KeysAndQueries["id:Radius"]);
                double delta = Convert.ToDouble(si.KeysAndQueries["id:delta"]);

                EdgeMask em = new EdgeMask(grd, X => Math.Abs(Math.Sqrt(X[0] * X[0] + Math.Pow(-delta / 2.0 - R + X[1], 2)) - R) < 1e-12); // lower cylinder
                var spline = SplineOnEdge(em, ShearRate, 0, out _, out _);
                double shearrate = spline.Interpolate(0.0);

                //Console.WriteLine("NipShearRateSubstrate for " + si.Name + " done!");

                object ret = Math.Abs(shearrate);
                return ret;
            }

            static public object NipShearRateCylinder(ISessionInfo si) {
                var fields = LoadFields(si);
                DGField VelocityX = fields.Single(f => f.Identification == "VelocityX");
                DGField VelocityY = fields.Single(f => f.Identification == "VelocityY");
                var grd = (GridData)VelocityX.GridDat;

                var ShearRate = (DGField)VelocityX.Clone();
                ShearRate.Clear();
                ShearRate.Derivative(1.0, VelocityX, 1);
                ShearRate.Derivative(1.0, VelocityY, 0);

                double R = Convert.ToDouble(si.KeysAndQueries["id:Radius"]);
                double delta = Convert.ToDouble(si.KeysAndQueries["id:delta"]);

                EdgeMask em = new EdgeMask(grd, X => Math.Abs(Math.Sqrt(X[0] * X[0] + Math.Pow(+delta / 2.0 + R - X[1], 2)) - R) < 1e-12); // lower cylinder
                var spline = SplineOnEdge(em, ShearRate, 0, out _, out _);
                double shearrate = spline.Interpolate(0.0);

                //Console.WriteLine("NipShearRateCylinder for " + si.Name + " done!");

                object ret = Math.Abs(shearrate);
                return ret;
            }

            static public object NipShearRate(ISessionInfo si) {
                var fields = LoadFields(si);
                DGField VelocityX = fields.Single(f => f.Identification == "VelocityX");
                DGField VelocityY = fields.Single(f => f.Identification == "VelocityY");
                var grd = (GridData)VelocityX.GridDat;

                var ShearRate = (DGField)VelocityX.Clone();
                ShearRate.Clear();
                ShearRate.Derivative(1.0, VelocityX, 1);
                ShearRate.Derivative(1.0, VelocityY, 0);

                double R = Convert.ToDouble(si.KeysAndQueries["id:Radius"]);
                double delta = Convert.ToDouble(si.KeysAndQueries["id:delta"]);

                EdgeMask em = new EdgeMask(grd, X => Math.Abs(X[1]) < 1e-12); // symmetry plane
                var spline = SplineOnEdge(em, ShearRate, 0, out _, out _);
                double shearrate = spline.Interpolate(0.0);

                //Console.WriteLine("NipShearRate for " + si.Name + " done!");

                object ret = Math.Abs(shearrate);
                return ret;
            }

            #endregion

            #region Integrated Shear stresses on both surfaces and Viscous Dissipation

            static public object ShearStressSubstrate(ISessionInfo si) {
                var fields = LoadFields(si);
                DGField VelocityX = fields.Single(f => f.Identification == "VelocityX");
                DGField VelocityY = fields.Single(f => f.Identification == "VelocityY");
                var grd = (GridData)VelocityX.GridDat;

                double R = Convert.ToDouble(si.KeysAndQueries["id:Radius"]);
                double delta = Convert.ToDouble(si.KeysAndQueries["id:delta"]);
                double Mu = Convert.ToDouble(si.KeysAndQueries["PhysicalParameters.mu_A"]);


                var ShearRateXY = (DGField)VelocityX.Clone();
                ShearRateXY.Clear();
                ShearRateXY.Derivative(1.0, VelocityX, 1);
                ShearRateXY.Derivative(1.0, VelocityY, 0);
                var ShearRateXX = (DGField)VelocityX.Clone();
                ShearRateXX.Clear();
                ShearRateXX.Derivative(2.0, VelocityX, 0);
                var ShearRateYY = (DGField)VelocityX.Clone();
                ShearRateYY.Clear();
                ShearRateYY.Derivative(2.0, VelocityY, 1);

                Func<double[], double[], double, double, double, double> ShearRate = (normal, tangent, t_xx, t_xy, t_yy) => {
                    double[] n = normal;
                    double[] t = tangent;

                    // stress in wall normal direction
                    double[] s_n = new double[] {
                    t_xx * n[0] + t_xy * n[1],
                    t_xy * n[0] + t_yy * n[1]
                };
                    double[] s_t = new double[] {
                    s_n[0] - (s_n[0] * n[0] + s_n[1] * n[1]) * n[0],
                    s_n[1] - (s_n[0] * n[0] + s_n[1] * n[1]) * n[1]
                };

                    return s_t[0] * t[0] + s_t[1] * t[1];
                };

                double eta0 = -MathNet.Numerics.Trig.Acosh(0.5 * delta / R + 1.0);
                double c = Math.Abs(MathNet.Numerics.Trig.Sinh(eta0)) * R;
                Func<double[], double> Eta = X => {
                    var x = X[0];
                    var y = X[1];
                    return 0.5 * Math.Log((Math.Pow(y + c, 2) + Math.Pow(x, 2)) / (Math.Pow(y - c, 2) + Math.Pow(x, 2)));
                };

                //EdgeMask em = new EdgeMask(grd, X => Math.Abs(Eta(X) - eta0) < 1e-8); // cylinder
                EdgeMask em = new EdgeMask(grd, X => Math.Abs(Math.Sqrt(X[0] * X[0] + Math.Pow(-delta / 2.0 - R + X[1], 2)) - R) < 1e-12); // lower cylinder
                int TransformDegree = grd.SpatialDimension * (InterpolationDegree(grd) - 1);

                double shearstress = 0.0;
                EdgeQuadrature.GetQuadrature(new int[] { 1 }, grd,
                    (new EdgeQuadratureScheme(true, em)).Compile(grd, VelocityX.Basis.Degree * 2 + TransformDegree),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        var EdgeNormals = grd.Edges.NormalsCache.GetNormals_Edge(QR.Nodes, i0, Length);

                        MultidimensionalArray ShearRateIN = MultidimensionalArray.Create(Length, QR.NoOfNodes, 3);
                        MultidimensionalArray ShearRateOT = MultidimensionalArray.Create(Length, QR.NoOfNodes, 3);

                        ((XDGField)ShearRateXX).GetSpeciesShadowField("A").EvaluateEdge(i0, Length, QR.Nodes, ShearRateIN.ExtractSubArrayShallow(-1, -1, 0), ShearRateOT.ExtractSubArrayShallow(-1, -1, 0), null, null, null, null, 0, 0.0);
                        ((XDGField)ShearRateXY).GetSpeciesShadowField("A").EvaluateEdge(i0, Length, QR.Nodes, ShearRateIN.ExtractSubArrayShallow(-1, -1, 1), ShearRateOT.ExtractSubArrayShallow(-1, -1, 1), null, null, null, null, 0, 0.0);
                        ((XDGField)ShearRateYY).GetSpeciesShadowField("A").EvaluateEdge(i0, Length, QR.Nodes, ShearRateIN.ExtractSubArrayShallow(-1, -1, 2), ShearRateOT.ExtractSubArrayShallow(-1, -1, 2), null, null, null, null, 0, 0.0);

                        int K = QR.NoOfNodes;
                        for (int i = 0; i < Length; i++) {
                            for (int k = 0; k < K; k++) {
                                double[] Tangent = new double[2];
                                double[] Normal = EdgeNormals.ExtractSubArrayShallow(i, k, -1).To1DArray();
                                Tangent[0] = Normal[1];
                                Tangent[1] = -Normal[0];

                                EvalResult[i, k] = ShearRate(Normal, Tangent, ShearRateIN[i, k, 0], ShearRateIN[i, k, 1], ShearRateIN[i, k, 2]);
                            }
                        }
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for (int _i = 0; _i < Length; _i++) {
                            shearstress += Mu * ResultsOfIntegration[_i, 0];
                        }
                    }
                ).Execute();
                //Console.WriteLine("ShearStressCylinder for " + si.Name + " done!");

                object ret = Math.Abs(shearstress);
                return ret;
            }

            static public object ShearStressCylinder(ISessionInfo si) {
                var fields = LoadFields(si);
                DGField VelocityX = fields.Single(f => f.Identification == "VelocityX");
                DGField VelocityY = fields.Single(f => f.Identification == "VelocityY");
                var grd = (GridData)VelocityX.GridDat;

                double R = Convert.ToDouble(si.KeysAndQueries["id:Radius"]);
                double delta = Convert.ToDouble(si.KeysAndQueries["id:delta"]);
                double Mu = Convert.ToDouble(si.KeysAndQueries["PhysicalParameters.mu_A"]);


                var ShearRateXY = (DGField)VelocityX.Clone();
                ShearRateXY.Clear();
                ShearRateXY.Derivative(1.0, VelocityX, 1);
                ShearRateXY.Derivative(1.0, VelocityY, 0);
                var ShearRateXX = (DGField)VelocityX.Clone();
                ShearRateXX.Clear();
                ShearRateXX.Derivative(2.0, VelocityX, 0);
                var ShearRateYY = (DGField)VelocityX.Clone();
                ShearRateYY.Clear();
                ShearRateYY.Derivative(2.0, VelocityY, 1);

                Func<double[], double[], double, double, double, double> ShearRate = (normal, tangent, t_xx, t_xy, t_yy) => {
                    double[] n = normal;
                    double[] t = tangent;

                    // stress in wall normal direction
                    double[] s_n = new double[] {
                    t_xx * n[0] + t_xy * n[1],
                    t_xy * n[0] + t_yy * n[1]
                };
                    double[] s_t = new double[] {
                    s_n[0] - (s_n[0] * n[0] + s_n[1] * n[1]) * n[0],
                    s_n[1] - (s_n[0] * n[0] + s_n[1] * n[1]) * n[1]
                };

                    return s_t[0] * t[0] + s_t[1] * t[1];
                };

                double eta0 = MathNet.Numerics.Trig.Acosh(0.5 * delta / R + 1.0);
                double c = Math.Abs(MathNet.Numerics.Trig.Sinh(eta0)) * R;
                Func<double[], double> Eta = X => {
                    var x = X[0];
                    var y = X[1];
                    return 0.5 * Math.Log((Math.Pow(y + c, 2) + Math.Pow(x, 2)) / (Math.Pow(y - c, 2) + Math.Pow(x, 2)));
                };

                //EdgeMask em = new EdgeMask(grd, X => Math.Abs(Eta(X) - eta0) < 1e-8); // cylinder
                EdgeMask em = new EdgeMask(grd, X => Math.Abs(Math.Sqrt(X[0] * X[0] + Math.Pow(delta / 2.0 + R - X[1], 2)) - R) < 1e-12); // upper cylinder
                int TransformDegree = grd.SpatialDimension * (InterpolationDegree(grd) - 1);

                double shearstress = 0.0;
                EdgeQuadrature.GetQuadrature(new int[] { 1 }, grd,
                    (new EdgeQuadratureScheme(true, em)).Compile(grd, VelocityX.Basis.Degree * 2 + TransformDegree),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        var EdgeNormals = grd.Edges.NormalsCache.GetNormals_Edge(QR.Nodes, i0, Length);

                        MultidimensionalArray ShearRateIN = MultidimensionalArray.Create(Length, QR.NoOfNodes, 3);
                        MultidimensionalArray ShearRateOT = MultidimensionalArray.Create(Length, QR.NoOfNodes, 3);

                        ((XDGField)ShearRateXX).GetSpeciesShadowField("A").EvaluateEdge(i0, Length, QR.Nodes, ShearRateIN.ExtractSubArrayShallow(-1, -1, 0), ShearRateOT.ExtractSubArrayShallow(-1, -1, 0), null, null, null, null, 0, 0.0);
                        ((XDGField)ShearRateXY).GetSpeciesShadowField("A").EvaluateEdge(i0, Length, QR.Nodes, ShearRateIN.ExtractSubArrayShallow(-1, -1, 1), ShearRateOT.ExtractSubArrayShallow(-1, -1, 1), null, null, null, null, 0, 0.0);
                        ((XDGField)ShearRateYY).GetSpeciesShadowField("A").EvaluateEdge(i0, Length, QR.Nodes, ShearRateIN.ExtractSubArrayShallow(-1, -1, 2), ShearRateOT.ExtractSubArrayShallow(-1, -1, 2), null, null, null, null, 0, 0.0);

                        int K = QR.NoOfNodes;
                        for (int i = 0; i < Length; i++) {
                            for (int k = 0; k < K; k++) {
                                double[] Tangent = new double[2];
                                double[] Normal = EdgeNormals.ExtractSubArrayShallow(i, k, -1).To1DArray();
                                Tangent[0] = Normal[1];
                                Tangent[1] = -Normal[0];

                                EvalResult[i, k] = ShearRate(Normal, Tangent, ShearRateIN[i, k, 0], ShearRateIN[i, k, 1], ShearRateIN[i, k, 2]);
                            }
                        }
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for (int _i = 0; _i < Length; _i++) {
                            shearstress += Mu * ResultsOfIntegration[_i, 0];
                        }
                    }
                ).Execute();
                //Console.WriteLine("ShearStressCylinder for " + si.Name + " done!");

                object ret = Math.Abs(shearstress);
                return ret;
            }

            // Frobenius product of Shear-stress tensor and velocity gradient
            static public object ViscousDissipation(ISessionInfo si) {
                var fields = LoadFields(si);
                DGField VelocityX = fields.Single(f => f.Identification == "VelocityX");
                DGField VelocityY = fields.Single(f => f.Identification == "VelocityY");
                var grd = (GridData)VelocityX.GridDat;

                double R = Convert.ToDouble(si.KeysAndQueries["id:Radius"]);
                double delta = Convert.ToDouble(si.KeysAndQueries["id:delta"]);
                double Mu = Convert.ToDouble(si.KeysAndQueries["PhysicalParameters.mu_A"]);


                ScalarFunctionEx changerate_dEspc = EnergyUtils.GetSpeciesKineticDissipationFunc(new DGField[] { VelocityX, VelocityY }, Mu);
                int TransformDegree = grd.SpatialDimension * (InterpolationDegree(grd) - 1);

                double diss = 0.0;
                CellQuadrature.GetQuadrature(new int[] { 1 }, grd,
                    (new CellQuadratureScheme(true, CellMask.GetFullMask(grd))).Compile(grd, VelocityX.Basis.Degree * 2 + TransformDegree),
                    delegate (int i0, int Length, QuadRule QR, MultidimensionalArray EvalResult) {
                        changerate_dEspc(i0, Length, QR.Nodes, EvalResult.ExtractSubArrayShallow(-1, -1, 0));
                    },
                    delegate (int i0, int Length, MultidimensionalArray ResultsOfIntegration) {
                        for (int i = 0; i < Length; i++)
                            diss += ResultsOfIntegration[i, 0];
                    }
                ).Execute();

                //Console.WriteLine("ViscousDissipation for " + si.Name + " done!");

                object ret = Math.Abs(diss);
                return ret;
            }

            #endregion

            #region Characteristic Flow Positions

            static public object PositionOfSynchronousFlow(ISessionInfo si) {
                var fields = LoadFields(si);
                DGField VelocityX = fields.Single(f => f.Identification == "VelocityX");
                var grd = (GridData)VelocityX.GridDat;

                double V_Soll = Convert.ToDouble(si.KeysAndQueries["id:V_Wall"]);

                EdgeMask em = new EdgeMask(grd, X => Math.Abs(X[1]) < 1e-12 && X[0] > 0); // symmetry axis in positive halfplane

                double x = PosOfValueOnEdge(V_Soll, em, VelocityX);

                //Console.WriteLine("PositionOfSynchronousFlow for " + si.Name + " done!");

                object ret = x;
                return ret;
            }

            static public object PositionOfStagnatingFlow(ISessionInfo si) {
                var fields = LoadFields(si);
                DGField VelocityX = fields.Single(f => f.Identification == "VelocityX");
                var grd = (GridData)VelocityX.GridDat;

                double V_Soll = 0.0;

                EdgeMask em = new EdgeMask(grd, X => Math.Abs(X[1]) < 1e-12 && X[0] > 0); // symmetry axis in positive halfplane

                double x = PosOfValueOnEdge(V_Soll, em, VelocityX);


                //Console.WriteLine("PositionOfStagnatingFlow for " + si.Name + " done!");

                object ret = x;
                return ret;
            }

            #endregion

            #region Pressure Gradient at Characteristic Flow Positions

            static public object dPdXatSynchronousPoint(ISessionInfo si) {
                var fields = LoadFields(si);
                DGField VelocityX = fields.Single(f => f.Identification == "VelocityX");
                DGField Pressure = fields.Single(f => f.Identification == "Pressure");

                var grd = (GridData)VelocityX.GridDat;

                var dP = Pressure.CloneAs();
                dP.Clear();
                dP.Derivative(1.0, Pressure, 0);

                double V_Soll = Convert.ToDouble(si.KeysAndQueries["id:V_Wall"]);

                EdgeMask em = new EdgeMask(grd, X => Math.Abs(X[1]) < 1e-12 && X[0] > 0); // symmetry axis in positive halfplane

                double x = PosOfValueOnEdge(V_Soll, em, VelocityX);
                var spline = SplineOnEdge(em, dP, 0, out _, out _);
                double p = spline.Interpolate(x);

                //Console.WriteLine("dPdXatSynchronousPoint for " + si.Name + " done!");

                object ret = Math.Abs(p);
                return ret;
            }

            static public object dPdXatStagnationPoint(ISessionInfo si) {
                var fields = LoadFields(si);
                DGField VelocityX = fields.Single(f => f.Identification == "VelocityX");
                DGField Pressure = fields.Single(f => f.Identification == "Pressure");

                var grd = (GridData)VelocityX.GridDat;

                var dP = Pressure.CloneAs();
                dP.Clear();
                dP.Derivative(1.0, Pressure, 0);

                double V_Soll = 0.0;

                EdgeMask em = new EdgeMask(grd, X => Math.Abs(X[1]) < 1e-12 && X[0] > 0); // symmetry axis in positive halfplane

                double x = PosOfValueOnEdge(V_Soll, em, VelocityX);
                var spline = SplineOnEdge(em, dP, 0, out _, out _);
                double p = spline.Interpolate(x);

                //Console.WriteLine("dPdXatStagnationPoint for " + si.Name + " done!");

                object ret = Math.Abs(p);
                return ret;
            }

            static public object dPdXatNip(ISessionInfo si) {
                var fields = LoadFields(si);
                DGField VelocityX = fields.Single(f => f.Identification == "VelocityX");
                DGField Pressure = fields.Single(f => f.Identification == "Pressure");

                var grd = (GridData)VelocityX.GridDat;

                var dP = Pressure.CloneAs();
                dP.Clear();
                dP.Derivative(1.0, Pressure, 0);

                EdgeMask em = new EdgeMask(grd, X => Math.Abs(X[1]) < 1e-12); // symmetry axis in positive halfplane

                var spline = SplineOnEdge(em, dP, 0, out _, out _);
                double p = spline.Interpolate(0.0);

                //Console.WriteLine("dPdXatNip for " + si.Name + " done!");

                object ret = Math.Abs(p);
                return ret;
            }

            static public object dPdXatConstantX(ISessionInfo si) {
                var fields = LoadFields(si);
                DGField VelocityX = fields.Single(f => f.Identification == "VelocityX");
                DGField Pressure = fields.Single(f => f.Identification == "Pressure");

                var grd = (GridData)VelocityX.GridDat;

                double x = 0.009;

                var dP = Pressure.CloneAs();
                dP.Clear();
                dP.Derivative(1.0, Pressure, 0);

                EdgeMask em = new EdgeMask(grd, X => Math.Abs(X[1]) < 1e-12); // symmetry axis in positive halfplane

                var spline = SplineOnEdge(em, dP, 0, out _, out _);
                double p = spline.Interpolate(x);

                //Console.WriteLine("dPdXat " + x + " for " + si.Name + " done!");

                object ret = Math.Abs(p);
                return ret;
            }
            #endregion
        }
    


    #endregion


}
