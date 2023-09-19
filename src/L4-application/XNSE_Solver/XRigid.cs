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

namespace BoSSS.Application.XNSE_Solver {

    public enum Shape {
        None = 0,
        Sphere = 1,
        Cube = 2,
        Torus = 3
    }

    [DataContract]
    [Serializable]
    public class XRigid {

        [DataMember]
        private double[] m_pos; //center of mass
        [DataMember]
        private double m_anglevelocity = -1.0;
        [DataMember]
        private double m_partRadius = -1.0;
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

        [NonSerialized]
        private XNSE_Control m_ctrl;


        [DataMember]
        private string m_RotationAxis = "z";

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
        public void SetParameters(double[] pos, double anglevelocity, double partRadius, int SpaceDim, double ringRadius = 0) {
            m_pos = pos;
            m_anglevelocity = anglevelocity;
            m_partRadius = partRadius;
            m_SpaceDim = SpaceDim;
            m_ringRadius = ringRadius;
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
                default:
                    throw new NotSupportedException();
            }
            SetVelocityAtIB();
        }

        private void DefineSphere() {
            var pos = m_pos;
            var anglevelocity = m_anglevelocity;
            var SpaceDim = m_SpaceDim;
            var particleRad = m_partRadius;
            var RotationCenter = m_RotationCenter;
            var RotationAxis = m_RotationAxis;
            m_ctrl.Tags.Add("Sphere");
            m_ctrl.LSContiProjectionMethod = ContinuityProjectionOption.None;

            Func<double[], double, double> PhiFunc = delegate (double[] X, double t) {
                double[] RotationArm = new double[SpaceDim];
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

                if (!Enumerable.SequenceEqual(pos, RotationCenter))
                    RotationArm = pos.Zip(RotationCenter, (p, RC) => p - RC).ToArray(); // subtraction: pos - RotationCenter 
                AffineTrafo affineTrafoFinal;
                if (SpaceDim == 2) {
                    affineTrafoFinal = AffineTrafo.Some2DRotation(angle);
                } else {
                    affineTrafoFinal = AffineTrafo.Rotation3D(rotAxis, angle);
                }

                double[] rotated_arm = affineTrafoFinal.Transform(RotationArm);
                double[] rotated_pos = RotationCenter.Zip(rotated_arm, (RC,rA) => RC + rA).ToArray(); // sum: RotationCenter + affineTrafoFinal.Transform(RotationArm) 

                switch (SpaceDim) {
                    case 2:
                    // circle
                    return -(X[0] - rotated_pos[0]) * (X[0] - rotated_pos[0]) - (X[1] - rotated_pos[1]) * (X[1] - rotated_pos[1]) + particleRad * particleRad;

                    case 3:
                    // sphere
                    return -(X[0] - rotated_pos[0]) * (X[0] - rotated_pos[0]) - (X[1] - rotated_pos[1]) * (X[1] - rotated_pos[1]) - (X[2] - pos[2]) * (X[2] - rotated_pos[2]) + particleRad * particleRad;

                    default:
                    throw new NotImplementedException();
                }
            };
            SetPhi(PhiFunc);
        }

        private void DefineCube() {
            var pos = m_pos;
            var anglevelocity = m_anglevelocity;
            var SpaceDim = m_SpaceDim;
            var particleRad = m_partRadius;
            var tiltVector = m_tiltVector;
            var tiltDegree = m_tiltDegree;
            var RotationAxis = m_RotationAxis;

            m_ctrl.Tags.Add("Cube");
            m_ctrl.LSContiProjectionMethod = ContinuityProjectionOption.ConstrainedDG;

            Func<double[], double, double> PhiFunc = delegate (double[] x, double t) {
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
            var anglevelocity = m_anglevelocity;
            var SpaceDim = m_SpaceDim;
            var particleRad = m_partRadius;
            var ringRad = m_ringRadius;
            var tiltVector = m_tiltVector;
            var tiltDegree = m_tiltDegree;
            var RotationAxis = m_RotationAxis;

            m_ctrl.Tags.Add("Torus");
            m_ctrl.LSContiProjectionMethod = ContinuityProjectionOption.ConstrainedDG;

            Func<double[], double, double> PhiFunc = delegate (double[] x, double t) {
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

                var affineTrafoRot = AffineTrafo.Rotation3D(rotAxis, angle);
                var affineTrafoFinal = affineTrafoRot * affineTrafoTilt;

                double[] X = new double[] { 0, 0, 0 };

                // Update intermediate variable X w.r.t. input x 
                for (int d=0; d<x.Length; d++)
                    X[d] = x[d];

                X = affineTrafoFinal.Transform(X);

                var TorusObject = new BoSSS.Solution.LevelSetTools.TestCases.Torus(particleRad, ringRad);
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

        private void SetVelocityAtIB() {
            var pos = m_pos;
            var rotCenter = m_RotationCenter;
            var anglevelocity = m_anglevelocity;
            var SpaceDim = m_SpaceDim;
            var tiltVector = m_tiltVector;
            var tiltDegree = m_tiltDegree;

            Func<double[], double, double[]> VelocityAtIB = delegate (double[] x, double time) {
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
