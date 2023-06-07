using BoSSS.Foundation;
using BoSSS.Solution.Control;
using BoSSS.Solution.Utils;
using ilPSP;
using ilPSP.Utils;
using Newtonsoft.Json;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.XNSE_Solver.SpecificSolutions {

    /// <summary>
    /// Rotationally symmetric initial values for velocity in polar coordinates
    /// </summary>
    /// <remarks>
    /// Implemented for the DACH-project on droplet oscillation (cooperation w. Prof Brenn / TU Graz),
    /// in order to compare analytical solutions to numerical simulations.
    /// The initial values from the analytical solution are provided in polar coordinates.
    /// </remarks>
    [Serializable]
    public class PolarAxiallySymmetricInitialValues : IBoundaryAndInitialData {


        /// <summary>
        /// All angles at which velocity values are provided
        /// </summary>
        [JsonProperty]
        double[] Theta;

        /// <summary>
        /// For each angle, all radii at which velocity values are provided
        /// - 1st index: correlates with <see cref="Theta"/>
        /// - 2nd index: enumeration of radii samples
        /// </summary>
        [JsonProperty]
        double[][] Radius;

        /// <summary>
        /// polar velocity component, indexing correlates with <see cref="Radius"/>
        /// </summary>
        [JsonProperty]
        double[][] PolarVel;

        /// <summary>
        /// radial velocity component, indexing correlates with <see cref="Radius"/>
        /// </summary>
        [JsonProperty]
        double[][] RadialVel;

        /// <summary>
        /// spatial direction
        /// </summary>
        [JsonIgnore]
        public int VelocityComponent {
            get {
                return m_VelocityComponent;
            }
            set {
                if(value < 0 || value >= 3)
                    throw new ArgumentOutOfRangeException("illegal spatial dimension index (expecting either 0, 1 or 2).");
                m_VelocityComponent = value;
            }
        }

        [JsonProperty]
        int m_VelocityComponent = 0;

        /// <summary>
        /// 
        /// </summary>
        public void SetData(double[] thetas, double[] radii, double[] polarVel, double[] radialVel) {
            int L = thetas.Length;
            if(radii.Length != L)
                throw new ArgumentException("Expecting all input arrays to be of equal length (radii)");
            if(polarVel.Length != L)
                throw new ArgumentException("Expecting all input arrays to be of equal length (polarVel)");
            if(radialVel.Length != L)
                throw new ArgumentException("Expecting all input arrays to be of equal length (radialVel)");


            var tempData = new SortedDictionary<double, SortedDictionary<double, (double PolarVel, double radialVel)>>();

            for(int l = 0; l < L; l++) {
                double th = thetas[l];
                double ra = radii[l];
                double vp = polarVel[l];
                double vr = radialVel[l];

                if(!tempData.TryGetValue(th, out var tempData_th)) {
                    tempData_th = new SortedDictionary<double, (double PolarVel, double radialVel)>();
                    tempData.Add(th, tempData_th);
                }

                tempData_th.Add(ra, (vp, vr));
            }

            this.Theta = tempData.Keys.ToArray();
            this.Radius = new double[this.Theta.Length][];
            this.PolarVel = new double[this.Theta.Length][];
            this.RadialVel = new double[this.Theta.Length][];
            for(int i = 0; i < Theta.Length; i++) {
                if(i > 0) {
                    if(this.Theta[i - 1] >= this.Theta[i])
                        throw new ApplicationException();
                }

                var tempData_th = tempData[this.Theta[i]];

                this.Radius[i] = tempData_th.Keys.ToArray();
                this.PolarVel[i] = this.Radius[i].Select(ra => tempData_th[ra].PolarVel).ToArray();
                this.RadialVel[i] = this.Radius[i].Select(ra => tempData_th[ra].radialVel).ToArray();
            }


        }


        /// <summary>
        /// evaluation
        /// </summary>
        public double Evaluate(double[] _X, double t) {
            Vector X = _X;
            if(X.Dim != 3)
                throw new NotSupportedException();

            double xy = Math.Sqrt(X.x.Pow2() + X.y.Pow2());
            double z = X.z;
            double Dx = X.x / xy;
            double Dy = X.y / xy;

            double theta = Math.Atan2(xy, z); // yes this is verified.
            double radius = X.L2Norm();


            (int i1, int i2) = binarySearch(theta, Theta);
            (int j_11, int j_12) = binarySearch(radius, Radius[i1]);
            (int j_21, int j_22) = binarySearch(radius, Radius[i2]);

            double polar1 = Interpolate(radius, Radius[i1][j_11], Radius[i1][j_12], PolarVel[i1][j_11], PolarVel[i1][j_12]);
            double polar2 = Interpolate(radius, Radius[i2][j_21], Radius[i2][j_22], PolarVel[i2][j_21], PolarVel[i2][j_22]);
            double v_theta = Interpolate(theta, Theta[i1], Theta[i2], polar1, polar2);

            double radial1 = Interpolate(radius, Radius[i1][j_11], Radius[i1][j_12], RadialVel[i1][j_11], RadialVel[i1][j_12]);
            double radial2 = Interpolate(radius, Radius[i2][j_21], Radius[i2][j_22], RadialVel[i2][j_21], RadialVel[i2][j_22]);
            double v_radial = Interpolate(theta, Theta[i1], Theta[i2], radial1, radial2);


            //double Vxy = Math.Cos(theta)*v_radial
            //    + Math.Sin(theta)*v_theta;
            //double Vz = Math.Sin(theta)*v_radial
            //    - Math.Cos(theta)*v_theta;
            double Vxy = v_radial * Math.Sin(theta)
                + v_theta * Math.Cos(theta);
            double Vz = v_radial * Math.Cos(theta)
                - v_theta * Math.Sin(theta);
     

            var Vel = new Vector(Dx * Vxy, Dy * Vxy, Vz);

            var ret = Vel[VelocityComponent];

            if(ret.IsNaNorInf()) {
                throw new ArithmeticException();
            }

            return ret;
        }

        /// <summary>
        /// vectorized evaluation
        /// </summary>
        public void EvaluateV(MultidimensionalArray input, double time, MultidimensionalArray output) {
            ScalarFunction sF = NonVectorizedScalarFunction.Vectorize(this.Evaluate, time);
            sF(input, output);
        }


        double Interpolate(double x, double x1, double x2, double y1, double y2) {
            if(Math.Abs(x2 - x1) <= Math.Abs(y2 - y1) * 1e-10)
                return 0.5 * (y2 + y1);

            double k = (y2 - y1) / (x2 - x1);
            double y = k * (x - x1) + y1;

            return y;
        }



        /// <summary>
        /// this function finds the interval which x belongs to
        /// </summary>
        /// <param name="x">
        /// The node that is to be placed between the interpolation points
        /// </param>
        /// <returns>
        /// To which interval this node belongs taking into account the interpolation points
        /// </returns>
        static private (int i1, int i2) binarySearch(double x, double[] nodes) {
            if(x <= nodes[0]) {
                return (0, 0);
            } else if(x >= nodes[nodes.Length - 1]) {
                return (nodes.Length - 1, nodes.Length - 1);
            } else {

                int idx = binarySearch_impl(x, nodes);
                if(x > nodes[idx + 1])
                    throw new ApplicationException("overshoot");
                if(x < nodes[idx])
                    throw new ApplicationException("undershoot");
                return (idx, idx + 1);
            }
        }

        static private int binarySearch_impl(double x, double[] nodes) {
            //if(binS_mem >= 0 && this.nodes[binS_mem + 1] >= x && this.nodes[binS_mem] < x)
            //    return binS_mem;

            int min = 0;
            int max = nodes.Length - 1;
            int mid = 0;
            do {
                mid = (min + max) / 2;
                if (x > nodes[mid]) {
                    min = mid + 1;
                } else {
                    max = mid - 1;
                }
            }
            while (mid > 0 && min <= max && (!(nodes[mid] >= x && nodes[mid - 1] < x)) && (!(nodes[mid] < x && nodes[mid + 1] >= x)));
            if (mid > 0 && nodes[mid] >= x && nodes[mid - 1] < x) {
                //binS_mem = mid - 1;
                return (mid - 1);
            }
            //binS_mem = mid;
            return mid;
        }

        /// <summary>
        /// equality
        /// </summary>
        public override bool Equals(object obj) {
            var other = obj as PolarAxiallySymmetricInitialValues;
            if(other == null)
                return false;

            if(!Theta.ListEquals(other.Theta))
                return false;

            if(m_VelocityComponent != other.m_VelocityComponent)
                return false;

            if(!PolarVel.ListEquals(other.PolarVel, (l1, l2) => l1.ListEquals(l2)))
                return false;
            if(!RadialVel.ListEquals(other.RadialVel, (l1, l2) => l1.ListEquals(l2)))
                return false;
            if(!Radius.ListEquals(other.Radius, (l1, l2) => l1.ListEquals(l2)))
                return false;


            return true;
        }

        /// <summary>
        /// 
        /// </summary>
        public override int GetHashCode() {
            if(Theta == null)
                return 0;
            else
                return (int) Math.Round(Theta.L2Norm() * 1e5);
        }
    }
}
