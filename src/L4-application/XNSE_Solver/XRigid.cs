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

namespace BoSSS.Application.XNSE_Solver {

    public enum Shape {
        None = 0,
        Sphere = 1,
        Cube = 2
    }

    [DataContract]
    [Serializable]
    public class XRigid {
        
        [DataMember]
        private double[] m_pos;
        [DataMember]
        private double m_anglevelocity = -1.0;
        [DataMember]
        private double m_partRadius = -1.0;
        [DataMember]
        private int m_SpaceDim = 0;
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

        public void SetParameters(double[] pos, double anglevelocity, double partRadius, int SpaceDim) {
            m_pos = pos;
            m_anglevelocity = anglevelocity;
            m_partRadius = partRadius;
            m_SpaceDim = SpaceDim;
        }

        public void SpecifyShape(Shape shape) {
            theShape = shape;
        }

        public void SetRotationAxis(string Axis) {
            switch (m_SpaceDim) {
                case 3:
                    string[] verify = new string[] { "x", "y", "z" };
                    bool isSupported = verify.Any(s => s == Axis);
                    if (!isSupported)
                        throw new NotSupportedException("Axis not supported: " + Axis);
                break;
                case 2:
                    if (Axis!="z")
                        throw new NotSupportedException("Axis not supported: " + Axis);
                break;
            }
            m_RotationAxis = Axis;
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
            m_ctrl.Tags.Add("Sphere");
            m_ctrl.LSContiProjectionMethod = ContinuityProjectionOption.None;
            Func<double[], double, double> PhiFunc = delegate (double[] X, double t) {
                switch (SpaceDim) {
                    case 2:
                    // circle
                    return -(X[0] - pos[0]) * (X[0] - pos[0]) - (X[1] - pos[1]) * (X[1] - pos[1]) + particleRad * particleRad;

                    case 3:
                    // sphere
                    return -(X[0] - pos[0]) * (X[0] - pos[0]) - (X[1] - pos[1]) * (X[1] - pos[1]) - (X[2] - pos[2]) * (X[2] - pos[2]) + particleRad * particleRad;

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
            m_ctrl.Tags.Add("Cube");
            m_ctrl.LSContiProjectionMethod = ContinuityProjectionOption.ConstrainedDG;
            Func<double[], double, double> PhiFunc = delegate (double[] X, double t) {
                double angle = -(anglevelocity * t) % (2 * Math.PI);
                switch (SpaceDim) {
                    case 2:
                    // Inf-Norm square
                    return -Math.Max(Math.Abs((X[0] - pos[0]) * Math.Cos(angle) - (X[1] - pos[1]) * Math.Sin(angle)),
                        Math.Abs((X[0] - pos[0]) * Math.Sin(angle) + (X[1] - pos[1]) * Math.Cos(angle)))
                        + particleRad;

                    case 3:
                        switch (m_RotationAxis) {
                        case "x":
                        return -Math.Max(Math.Abs(X[0] - pos[0]),
                                    Math.Max(Math.Abs((X[1] - pos[1]) * Math.Cos(angle) - (X[2] - pos[2]) * Math.Sin(angle)),
                                    Math.Abs((X[1] - pos[1]) * Math.Sin(angle) + (X[2] - pos[2]) * Math.Cos(angle))))
                                    + particleRad;
                        case "y":
                        return -Math.Max(Math.Abs((X[0] - pos[0]) * Math.Cos(angle) + (X[2] - pos[2]) * Math.Sin(angle)),
                                    Math.Max(Math.Abs(X[1] - pos[1]),
                                    Math.Abs(-(X[0] - pos[0]) * Math.Sin(angle) + (X[2] - pos[2]) * Math.Cos(angle))))
                                    + particleRad;
                        case "z":
                        return -Math.Max(Math.Abs((X[0] - pos[0]) * Math.Cos(angle) - (X[1] - pos[1]) * Math.Sin(angle)),
                                    Math.Max(Math.Abs((X[0] - pos[0]) * Math.Sin(angle) + (X[1] - pos[1]) * Math.Cos(angle)),
                                    Math.Abs(X[2] - pos[2])))
                                    + particleRad;
                        default:
                            throw new NotSupportedException();
                        }
                        
                    default:
                    throw new NotImplementedException();
                }
            };
            SetPhi(PhiFunc);
        }

        private void SetVelocityAtIB() {
            var pos = m_pos;
            var anglevelocity = m_anglevelocity;
            var SpaceDim = m_SpaceDim;
            Func<double[], double, double[]> VelocityAtIB = delegate (double[] X, double time) {

                if (pos.Length != X.Length)
                    throw new ArgumentException("check dimension of center of mass");

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
                        throw new NotSupportedException("Axis not suppored");
                }

                Vector CenterofMass = new Vector(pos);
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
            m_ctrl.InitialValues_Evaluators_TimeDep.Add(VariableNames.LevelSetCGidx(1), PhiFunc);
        }
    }
}
