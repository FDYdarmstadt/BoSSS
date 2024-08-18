using BoSSS.Solution.LevelSetTools;
using BoSSS.Solution.NSECommon;
using ilPSP;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Runtime.Serialization;
using Newtonsoft.Json;
using BoSSS.Platform.LinAlg;
using System.Drawing;
using BoSSS.Solution.Control;

namespace BoSSS.Application.XNSE_Solver {

    public enum Shape {
        None = 0,
        Sphere = 1,
        Cube = 2,
        Torus = 3,
        CollidingSpheres = 4, //a temporary test case (highly experimental) TODO: Move this to SetRigidLevelSet and EvolveRigidLevelSet
        Popcorn = 5,
        MovingSphere = 6
    }

    [DataContract]
    [Serializable]
    public class XRigid {

        [DataMember]
        private double[] m_pos; //center of mass
        [DataMember]
        private double m_angleVelocity = -1.0;
        [DataMember]
        private double m_partRadius = -1.0;
        [DataMember]
        private double m_rateOfRadius = 0.0;
        [DataMember]
        private int m_SpaceDim = 0;
        [DataMember]
        private double m_ringRadius = -1.0;
        [DataMember]
        private double m_tiltDegree = 0.0;
        [DataMember]
        private double[] m_tiltVector = new double[] {0, 1, 0};
        [DataMember]
        private double[] m_RotationCenter;
        [DataMember]
        private Shape theShape = Shape.None;
        [DataMember]
        private bool m_staticShape = false;
        [DataMember]
        private bool m_isThereExactSolution = false;
        [DataMember]
        private string m_RotationAxis = "z";

        [NonSerialized]
        private XNSE_Control m_ctrl;


        public bool IsInitialized() {
            return m_SpaceDim > 0;
        }

        /// <summary>
        /// TODO: Move this to SetRigidLevelSet and EvolveRigidLevelSet
        /// </summary>
        public XRigid() {

        }


        /// <summary>
        /// TODO: Move this to SetRigidLevelSet and EvolveRigidLevelSet
        /// </summary>
        public void SetParameters(double[] pos, double angleVelocity, double majorRadius, int SpaceDim, double minorRadius = 0.0, double rateOfRadius = 0.0, bool staticShape = false) {
            m_pos = pos;
            m_angleVelocity = angleVelocity;
            m_partRadius = majorRadius;
            m_SpaceDim = SpaceDim;
            m_ringRadius = minorRadius;
            m_rateOfRadius = rateOfRadius;
            m_staticShape = staticShape;
        }

        public void SpecifyShape(Shape shape) {
            theShape = shape;
            if(shape == Shape.None)
                m_SpaceDim = 0;
        }

        public void SetRotationAxis(string Axis) {
            double[] rotCenter = m_pos;
            SetRotation(Axis, rotCenter);
        }

        public void SetRotation(string Axis, double[] RotationCenter) {
            if (m_SpaceDim != RotationCenter.Length)
                throw new ArgumentException($"RotationCenter should be in {m_SpaceDim}-dimension");

            switch (m_SpaceDim) {
                case 3:
                var verify = new string[] { "x", "y", "z" }.Take(m_SpaceDim);
                bool isSupported = verify.Any(s => s == Axis);
                if(!isSupported)
                    throw new NotSupportedException("Axis not supported: " + Axis);
                break;

                case 2:
                if(Axis != "z")
                    throw new NotSupportedException("Axis not supported: " + Axis + " (In 2D one can only rotate around z)");
                break;
            }
            m_RotationCenter = RotationCenter;
            m_RotationAxis = Axis;
        }

        /// <summary>
        /// Tilt the object around an axis
        /// </summary>
        /// <param name="TiltAxis">Vector around which the body is tilted</param>
        /// <param name="TiltDegree">Degree in radians within [0,2pi]</param>
        /// <exception cref="NotImplementedException"></exception>
        /// <exception cref="ArgumentOutOfRangeException"></exception>
        public void SetTilt(double[] TiltAxis, double TiltDegree) {
            if (TiltAxis.Length != 3)
                throw new NotImplementedException("Tilting the rigid body is only supported in 3d");

            if (TiltDegree > 2 * Math.PI || TiltDegree < 0)
                throw new ArgumentOutOfRangeException("Tilt degree should be in radians with in the range [0,2pi]");

            m_tiltVector = TiltAxis;
            m_tiltDegree = TiltDegree;
        }

        public void ArrangeAll(XNSE_Control ctrl) {
            m_ctrl = ctrl;
            switch (theShape) {
                case Shape.Sphere:
                    DefineSphere();
                    break;
                case Shape.Cube:
                    DefineCube();
                    break;
                case Shape.Torus:
                    DefineTorus();
                    break;
                case Shape.CollidingSpheres:
                    DefineCollidingSpheres();
                    break;
                case Shape.Popcorn:
                    DefinePopcorn();
                    break;
                case Shape.MovingSphere:
                    DefineMovingSphere();
                    break;
                default:
                    throw new NotSupportedException();
            }
            SetVelocityAtIB();
            if (m_isThereExactSolution)
                SetExactSolutionForSteadyRotatingSphere();
        }

        public void AssignExactSolutionForSteadyRotatingSphere()
        {
            m_isThereExactSolution = true;
        }

        public void SetExactSolutionForSteadyRotatingSphere() {
            m_ctrl.ExactSolutionVelocity = new Dictionary<string, Func<double[], double, double>[]>();
            //m_ctrl.ExactSolutionPressure = new Dictionary<string, Func<double[], double, double>>();

            m_ctrl.ExactSolutionVelocity.Add("A", new Func<double[], double, double>[] { (X, t) => -m_angleVelocity * m_partRadius * m_partRadius * X[1] / (X[0] * X[0] + X[1] * X[1]), (X, t) => m_angleVelocity * m_partRadius * m_partRadius * X[0] / (X[0] * X[0] + X[1] * X[1]) });
            m_ctrl.ExactSolutionVelocity.Add("C", new Func<double[], double, double>[] { (X, t) => 0, (X, t) => 0 });
        }

        private void DefineSphere() {
            var pos = m_pos;
            var anglevelocity = m_angleVelocity;
            var SpaceDim = m_SpaceDim;
            var particleRad = m_partRadius;
            var rateOfRadius = m_rateOfRadius;
            var RotationCenter = m_RotationCenter;
            var RotationAxis = m_RotationAxis;
            m_ctrl.Tags.Add("Sphere");
            m_ctrl.LSContiProjectionMethod = ContinuityProjectionOption.ConstrainedDG;

            Func<double[], double, double> PhiFunc = delegate (double[] X, double t) {
                if (m_staticShape) //static shape -> stay at the initial position
                    t = 0;

                double[] RotationArm = new double[SpaceDim];
                double angle = m_staticShape ? 0 : -(anglevelocity * t) % (2 * Math.PI);
                double dynamicRadius = rateOfRadius == 0.0 ? particleRad : Math.Max((1 + rateOfRadius * t) * particleRad, 0.0);
                Vector rotAxis;
                switch (RotationAxis) {
                    case "x":
                        rotAxis = new Vector(1, 0, 0);
                        break;
                    case "y":
                        rotAxis = new Vector(0, 1, 0);
                        break;
                    case "z":
                        rotAxis = new Vector(0, 0, 1);
                        break;
                    default:
                        throw new NotSupportedException();
                }

                if (!Enumerable.SequenceEqual(pos, RotationCenter))
                    RotationArm = pos.Zip(RotationCenter, (p, RC) => p - RC).ToArray(); // subtraction: pos - RotationCenter 
                AffineTrafo affineTrafoFinal;
                if (SpaceDim == 2) {
                    affineTrafoFinal = AffineTrafo.Some2DRotation(angle);
                } else {
                    affineTrafoFinal = AffineTrafo.Rotation3D(rotAxis, angle);
                }

                double[] rotated_arm = affineTrafoFinal.Transform(RotationArm);
                double[] rotated_pos = RotationCenter.Zip(rotated_arm, (RC, rA) => RC + rA).ToArray(); // sum: RotationCenter + affineTrafoFinal.Transform(RotationArm) 

                switch (SpaceDim) {
                    case 2:
                        // circle
                        return -(X[0] - rotated_pos[0]) * (X[0] - rotated_pos[0]) - (X[1] - rotated_pos[1]) * (X[1] - rotated_pos[1]) + dynamicRadius * dynamicRadius;

                    case 3:
                        // sphere
                        return -(X[0] - rotated_pos[0]) * (X[0] - rotated_pos[0]) - (X[1] - rotated_pos[1]) * (X[1] - rotated_pos[1]) - (X[2] - pos[2]) * (X[2] - rotated_pos[2]) + dynamicRadius * dynamicRadius;

                    default:
                        throw new NotImplementedException();
                }
            };
            SetPhi(PhiFunc);
        }

        private void DefineCollidingSpheres() {
            //var pos = m_pos;
            //var anglevelocity = m_angleVelocity;
            var SpaceDim = m_SpaceDim;
            var particleRad = m_partRadius;
            var rateOfRadius = m_rateOfRadius;
            //var RotationCenter = m_RotationCenter;
            //var RotationAxis = m_RotationAxis;
            m_ctrl.Tags.Add("CollidingSphere");
            m_ctrl.LSContiProjectionMethod = ContinuityProjectionOption.ConstrainedDG;

            Func<double[], double, double> PhiFunc = delegate (double[] X, double t) {
                if (m_staticShape) //static shape -> stay at the initial position
                    t = 0;

                double[] posL = new double[SpaceDim];
                posL[0] = - 1.5 * particleRad  + (rateOfRadius * t) * particleRad; // (initial pos + change)

                double[] posR = new double[SpaceDim];
                posR[0] = 1.5 * particleRad - (rateOfRadius * t) * particleRad;
                double L, R;

                switch (SpaceDim) {
                    case 2:
                        // circle
                        L = -(X[0] - posL[0]) * (X[0] - posL[0]) - (X[1] - posL[1]) * (X[1] - posL[1]) + particleRad * particleRad;
                        R = -(X[0] - posR[0]) * (X[0] - posR[0]) - (X[1] - posR[1]) * (X[1] - posR[1]) + particleRad * particleRad;

                        return Math.Max(L, R);
                    case 3:
                        // sphere
                        L = -(X[0] - posL[0]) * (X[0] - posL[0]) - (X[1] - posL[1]) * (X[1] - posL[1]) - (X[2] - posL[2]) * (X[2] - posL[2]) + particleRad * particleRad;
                        R = -(X[0] - posR[0]) * (X[0] - posR[0]) - (X[1] - posR[1]) * (X[1] - posR[1]) - (X[2] - posR[2]) * (X[2] - posR[2]) + particleRad * particleRad;

                        return Math.Max(L, R);
                    default:
                    throw new NotImplementedException();
                }
            };
            SetPhi(PhiFunc);
        }

        private void DefineMovingSphere()
        {
            //var pos = m_pos;
            //var anglevelocity = m_angleVelocity;
            var SpaceDim = m_SpaceDim;
            var particleRad = m_partRadius;
            var rateOfRadius = m_rateOfRadius;
            //var RotationCenter = m_RotationCenter;
            //var RotationAxis = m_RotationAxis;
            m_ctrl.Tags.Add("MovingSphere");
            m_ctrl.LSContiProjectionMethod = ContinuityProjectionOption.ConstrainedDG;

            Func<double[], double, double> PhiFunc = delegate (double[] X, double t) {
                if (m_staticShape) //static shape -> stay at the initial position
                    t = 0;

                double[] posL = new double[SpaceDim];
                posL[0] = -1.5 * particleRad + (rateOfRadius * t) * particleRad; // (initial pos + change)

                double L;

                switch (SpaceDim)
                {
                    case 2:
                        // circle
                        L = -(X[0] - posL[0]) * (X[0] - posL[0]) - (X[1] - posL[1]) * (X[1] - posL[1]) + particleRad * particleRad;

                        return L;
                    case 3:
                        // sphere
                        L = -(X[0] - posL[0]) * (X[0] - posL[0]) - (X[1] - posL[1]) * (X[1] - posL[1]) - (X[2] - posL[2]) * (X[2] - posL[2]) + particleRad * particleRad;

                        return L;
                    default:
                        throw new NotImplementedException();
                }
            };
            SetPhi(PhiFunc);
        }

        private void DefineCube() {
            var pos = m_pos;
            var anglevelocity = m_angleVelocity;
            var SpaceDim = m_SpaceDim;
            var particleRad = m_partRadius;
            var tiltVector = m_tiltVector;
            var tiltDegree = m_tiltDegree;
            var RotationAxis = m_RotationAxis;

            m_ctrl.Tags.Add("Cube");
            m_ctrl.LSContiProjectionMethod = ContinuityProjectionOption.ConstrainedDG;

            Func<double[], double, double> PhiFunc = delegate (double[] x, double t) {
                if (m_staticShape) //static shape -> stay at the initial position
                    t = 0;

                Vector TiltVector = new Vector(tiltVector);

                double angle = -(anglevelocity * t) % (2 * Math.PI);
                Vector rotAxis;

                switch (RotationAxis) {
                    case "x":
                        rotAxis = new Vector(1, 0, 0);
                        break;
                    case "y":
                        rotAxis = new Vector(0, 1, 0);
                        break;
                    case "z":
                        rotAxis = new Vector(0, 0, 1);
                        break;
                    default:
                        throw new NotSupportedException();
                }

                AffineTrafo affineTrafoTilt;
                if (tiltDegree == 0.0) {
                    affineTrafoTilt = AffineTrafo.Identity(TiltVector.Dim);
                } else {
                    affineTrafoTilt = AffineTrafo.Rotation3D(TiltVector, tiltDegree);
                }

                double[] X = new double[] { 0, 0, 0 };

                var affineTrafoRot = AffineTrafo.Rotation3D(rotAxis, angle);
                var affineTrafoFinal = affineTrafoRot * affineTrafoTilt;

                X = affineTrafoFinal.Transform(x);

                var CubeObject = new BoSSS.Solution.LevelSetTools.TestCases.Cube(particleRad);
                switch (SpaceDim) {
                    case 2:
                        return CubeObject.SignedDistance2D(X);
                    case 3:
                        return CubeObject.SignedDistance(X);
                    //return Math.Pow(particleRad - Math.Sqrt(X[0].Pow2() + X[1].Pow2()), 2) + X[2] * X[2] - ringRad.Pow2();
                    default:
                        throw new NotImplementedException();
                }
            };
            //Func<double[], double, double> PhiFunc = delegate (double[] X, double t) {
            //    double angle = -(anglevelocity * t) % (2 * Math.PI);
            //    switch (SpaceDim) {
            //        case 2:
            //        // Inf-Norm square
            //        return -Math.Max(Math.Abs((X[0] - pos[0]) * Math.Cos(angle) - (X[1] - pos[1]) * Math.Sin(angle)),
            //            Math.Abs((X[0] - pos[0]) * Math.Sin(angle) + (X[1] - pos[1]) * Math.Cos(angle)))
            //            + particleRad;

            //        case 3:
            //            switch (m_RotationAxis) {
            //            case "x":
            //            return -Math.Max(Math.Abs(X[0] - pos[0]),
            //                        Math.Max(Math.Abs((X[1] - pos[1]) * Math.Cos(angle) - (X[2] - pos[2]) * Math.Sin(angle)),
            //                        Math.Abs((X[1] - pos[1]) * Math.Sin(angle) + (X[2] - pos[2]) * Math.Cos(angle))))
            //                        + particleRad;
            //            case "y":
            //            return -Math.Max(Math.Abs((X[0] - pos[0]) * Math.Cos(angle) + (X[2] - pos[2]) * Math.Sin(angle)),
            //                        Math.Max(Math.Abs(X[1] - pos[1]),
            //                        Math.Abs(-(X[0] - pos[0]) * Math.Sin(angle) + (X[2] - pos[2]) * Math.Cos(angle))))
            //                        + particleRad;
            //            case "z":
            //            return -Math.Max(Math.Abs((X[0] - pos[0]) * Math.Cos(angle) - (X[1] - pos[1]) * Math.Sin(angle)),
            //                        Math.Max(Math.Abs((X[0] - pos[0]) * Math.Sin(angle) + (X[1] - pos[1]) * Math.Cos(angle)),
            //                        Math.Abs(X[2] - pos[2])))
            //                        + particleRad;
            //            default:
            //                throw new NotSupportedException();
            //            }

            //        default:
            //        throw new NotImplementedException();
            //    }
            //};
            SetPhi(PhiFunc);
        }

        private void DefineTorus() {
            var pos = m_pos;
            var anglevelocity = m_angleVelocity;
            var SpaceDim = m_SpaceDim;
            var particleRad = m_partRadius;
            var rateOfRadius = m_rateOfRadius;
            var ringRad = m_ringRadius;
            var tiltVector = m_tiltVector;
            var tiltDegree = m_tiltDegree;
            var RotationAxis = m_RotationAxis;

            m_ctrl.Tags.Add("Torus");
            m_ctrl.LSContiProjectionMethod = ContinuityProjectionOption.ConstrainedDG;

            Func<double[], double, double> PhiFunc = delegate (double[] x, double t) {
                if (m_staticShape) //static shape -> stay at the initial position
                    t = 0;

                Vector TiltVector = new Vector(tiltVector);

                double angle = -(anglevelocity * t) % (2 * Math.PI);
                double dynamicRadius = rateOfRadius == 0.0 ? particleRad : Math.Max((1 + rateOfRadius * t) * particleRad, 0.0);

                Vector rotAxis;

                switch (RotationAxis) {
                    case "x":
                        rotAxis = new Vector(1, 0, 0);
                        break;
                    case "y":
                        rotAxis = new Vector(0, 1, 0);
                        break;
                    case "z":
                        rotAxis = new Vector(0, 0, 1);
                        break;
                    default:
                        throw new NotSupportedException();
                }

                AffineTrafo affineTrafoTilt;
                if (tiltDegree == 0.0) {
                    affineTrafoTilt = AffineTrafo.Identity(TiltVector.Dim);
                } else {
                    affineTrafoTilt = AffineTrafo.Rotation3D(TiltVector, tiltDegree);
                }

                var affineTrafoRot = AffineTrafo.Rotation3D(rotAxis, angle);
                var affineTrafoFinal = affineTrafoRot * affineTrafoTilt;

                double[] X = new double[] { 0, 0, 0 };

                // Update intermediate variable X w.r.t. input x 
                for (int d=0; d<x.Length; d++)
                    X[d] = x[d];

                X = affineTrafoFinal.Transform(X);

                var TorusObject = new BoSSS.Solution.LevelSetTools.TestCases.Torus(dynamicRadius, ringRad);
                switch (SpaceDim) {
                    case 2:
                        return -TorusObject.SignedDistance2D(X);
                    case 3:
                        return -TorusObject.SignedDistance(X);
                        //return Math.Pow(particleRad - Math.Sqrt(X[0].Pow2() + X[1].Pow2()), 2) + X[2] * X[2] - ringRad.Pow2();
                    default:
                        throw new NotImplementedException();
                }
            };
            SetPhi(PhiFunc);
        }

        private void DefinePopcorn() {
            var pos = m_pos;
            var anglevelocity = m_angleVelocity;
            var SpaceDim = m_SpaceDim;
            var particleRad = m_partRadius;
            var rateOfRadius = m_rateOfRadius;
            var ringRad = m_ringRadius;
            var tiltVector = m_tiltVector;
            var tiltDegree = m_tiltDegree;
            var RotationAxis = m_RotationAxis;

            m_ctrl.Tags.Add("Popcorn");
            m_ctrl.LSContiProjectionMethod = ContinuityProjectionOption.ConstrainedDG;

            Func<double[], double, double> PhiFunc = delegate (double[] x, double t) {
                if (m_staticShape) //static shape -> stay at the initial position
                    t = 0;

                Vector TiltVector = new Vector(tiltVector);

                double angle = -(anglevelocity * t) % (2 * Math.PI);
                double dynamicRadius = rateOfRadius == 0.0 ? particleRad : Math.Max((1 + rateOfRadius * t) * particleRad, 0.0);

                Vector rotAxis;

                switch (RotationAxis) {
                    case "x":
                        rotAxis = new Vector(1, 0, 0);
                        break;
                    case "y":
                        rotAxis = new Vector(0, 1, 0);
                        break;
                    case "z":
                        rotAxis = new Vector(0, 0, 1);
                        break;
                    default:
                        throw new NotSupportedException();
                }

                AffineTrafo affineTrafoTilt;
                if (tiltDegree == 0.0) {
                    affineTrafoTilt = AffineTrafo.Identity(TiltVector.Dim);
                } else {
                    affineTrafoTilt = AffineTrafo.Rotation3D(TiltVector, tiltDegree);
                }

                var affineTrafoRot = AffineTrafo.Rotation3D(rotAxis, angle);
                var affineTrafoFinal = affineTrafoRot * affineTrafoTilt;

                double[] X = new double[] { 0, 0, 0 };

                // Update intermediate variable X w.r.t. input x 
                for (int d = 0; d < x.Length; d++)
                    X[d] = x[d];

                X = affineTrafoFinal.Transform(X);

                var PopcornObject = new BoSSS.Solution.LevelSetTools.TestCases.Popcorn(dynamicRadius);
                switch (SpaceDim) {
                    case 2:
                        return -PopcornObject.LevelSetFunction2D(X); // (A<0, I=0, C>0)
                    case 3:
                        return -PopcornObject.LevelSetFunction(X); // (A<0, I=0, C>0)
                    //return Math.Pow(particleRad - Math.Sqrt(X[0].Pow2() + X[1].Pow2()), 2) + X[2] * X[2] - ringRad.Pow2();
                    default:
                        throw new NotImplementedException();
                }
            };
            SetPhi(PhiFunc);
        }

        private void SetVelocityAtIB() {
            var pos = m_pos;
            var rotCenter = m_RotationCenter;
            var anglevelocity = m_angleVelocity;
            var SpaceDim = m_SpaceDim;
            var tiltVector = m_tiltVector;
            var tiltDegree = m_tiltDegree;
            var particleRad = m_partRadius;
            var rateOfRadius = m_rateOfRadius;
            Func<double[], double, double[]> VelocityAtIB;

            if (theShape == Shape.CollidingSpheres)
            {

                VelocityAtIB = delegate (double[] X, double t) {
                    if (pos.Length != X.Length)
                        throw new ArgumentException("check dimension of center of mass");


                    double[] posL = new double[SpaceDim];
                    posL[0] = -1.5 * particleRad + (rateOfRadius * t) * particleRad; // (initial pos + change)

                    double[] posR = new double[SpaceDim];
                    posR[0] = 1.5 * particleRad - (rateOfRadius * t) * particleRad;
                    double L, R;

                    switch (SpaceDim) {
                        case 2:
                            // circle
                            L = -(X[0] - posL[0]) * (X[0] - posL[0]) - (X[1] - posL[1]) * (X[1] - posL[1]) + particleRad * particleRad;
                            R = -(X[0] - posR[0]) * (X[0] - posR[0]) - (X[1] - posR[1]) * (X[1] - posR[1]) + particleRad * particleRad;
                            break;
                        case 3:
                            // sphere
                            L = -(X[0] - posL[0]) * (X[0] - posL[0]) - (X[1] - posL[1]) * (X[1] - posL[1]) - (X[2] - posL[2]) * (X[2] - posL[2]) + particleRad * particleRad;
                            R = -(X[0] - posR[0]) * (X[0] - posR[0]) - (X[1] - posR[1]) * (X[1] - posR[1]) - (X[2] - posR[2]) * (X[2] - posR[2]) + particleRad * particleRad;
                            break;
                        default:
                            throw new NotImplementedException();
                    }

                    double[] pointVel = new double[SpaceDim];
                    double tol = -1E-8;
                    if (L > tol)
                    {
                        pointVel[0] += rateOfRadius;
                    }

                    if (R > tol) {
                        pointVel[0] -= rateOfRadius;
                    }

                    Vector pointVelocity = new Vector(pointVel);
                    return pointVelocity;
                };

            } else if(theShape == Shape.MovingSphere) {

                VelocityAtIB = delegate (double[] X, double t) {
                    if (pos.Length != X.Length)
                        throw new ArgumentException("check dimension of center of mass");


                    double[] posL = new double[SpaceDim];
                    posL[0] = -1.5 * particleRad + (rateOfRadius * t) * particleRad; // (initial pos + change)

                    double L;

                    switch (SpaceDim)
                    {
                        case 2:
                            // circle
                            L = -(X[0] - posL[0]) * (X[0] - posL[0]) - (X[1] - posL[1]) * (X[1] - posL[1]) + particleRad * particleRad;
                            break;
                        case 3:
                            // sphere
                            L = -(X[0] - posL[0]) * (X[0] - posL[0]) - (X[1] - posL[1]) * (X[1] - posL[1]) - (X[2] - posL[2]) * (X[2] - posL[2]) + particleRad * particleRad;
                            break;
                        default:
                            throw new NotImplementedException();
                    }

                    double[] pointVel = new double[SpaceDim];
                    double tol = -1E-8;
                    if (L > tol)
                    {
                        pointVel[0] += rateOfRadius;
                    }

                    Vector pointVelocity = new Vector(pointVel);
                    return pointVelocity;
                };

            } else {

            VelocityAtIB = delegate (double[] x, double time) {
                if (pos.Length != x.Length)
                    throw new ArgumentException("check dimension of center of mass");
                double[] X = x;

                //if (tiltDegree != 0.0 && SpaceDim == 3) {
                //    Vector TiltVector = new Vector(tiltVector);

                //    AffineTrafo affineTrafoTilt;
                //    if (tiltDegree == 0.0) {
                //        affineTrafoTilt = AffineTrafo.Identity(TiltVector.Dim);
                //    } else {
                //        affineTrafoTilt = AffineTrafo.Rotation3D(TiltVector, tiltDegree); //to define the rotation around the original axis of the rigid body, reverse the tilt
                //    }
                //    X = affineTrafoTilt.Transform(x);
                //}

                Vector angVelo = new Vector(new double[] { 0, 0, 0 });
                switch (m_RotationAxis) {
                    case "x":
                        angVelo.x = anglevelocity;
                        break;
                    case "y":
                        angVelo.y = anglevelocity;
                        break;
                    case "z":
                        angVelo.z = anglevelocity;
                        break;
                    default:
                        throw new NotSupportedException("Axis not supported");
                }

                Vector CenterofMass = new Vector(rotCenter); //rotation around the given point
                Vector radialVector = new Vector(X) - CenterofMass;
                Vector transVelocity = new Vector(new double[SpaceDim]);
                Vector pointVelocity;

                switch (SpaceDim) {
                    case 2:
                    pointVelocity = new Vector(transVelocity[0] - angVelo[2] * radialVector[1], transVelocity[1] + angVelo[2] * radialVector[0]);
                    break;
                    case 3:
                    pointVelocity = transVelocity + angVelo.CrossProduct(radialVector);
                    break;
                    default:
                    throw new NotImplementedException("this number of dimensions is not supported");
                }

                return pointVelocity;
            };
            }

            Func<double[], double, double> VelocityX = delegate (double[] X, double time) { return VelocityAtIB(X, time)[0]; };
            Func<double[], double, double> VelocityY = delegate (double[] X, double time) { return VelocityAtIB(X, time)[1]; };
            Func<double[], double, double> VelocityZ = delegate (double[] X, double time) { return VelocityAtIB(X, time)[2]; };
            m_ctrl.InitialValues_Evaluators_TimeDep.Add("VelocityX@Phi2", VelocityX);
            m_ctrl.InitialValues_Evaluators_TimeDep.Add("VelocityY@Phi2", VelocityY);
            if (SpaceDim == 3)
                m_ctrl.InitialValues_Evaluators_TimeDep.Add("VelocityZ@Phi2", VelocityZ);
        }

        private void SetPhi(Func<double[], double, double> PhiFunc) {
            m_ctrl.UseImmersedBoundary = true;
            m_ctrl.InitialValues_Evaluators_TimeDep.Add(VariableNames.LevelSetCGidx(1), PhiFunc);
        }
    }
}
