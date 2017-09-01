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
using System.Linq;
using System.Collections.Generic;
using BoSSS.Platform.LinAlg;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Solution.Utils;
using System.Collections;

namespace NSE_SIMPLE {

    /// <summary>
    /// Class for calculating drag and lift forces.
    /// </summary>
    class Force {

        NSE_SIMPLEMain m_app;

        SinglePhaseField DuDx;
        SinglePhaseField DuDy;
        SinglePhaseField DuDz;
        SinglePhaseField DvDx;
        SinglePhaseField DvDy;
        SinglePhaseField DvDz;
        SinglePhaseField DwDx;
        SinglePhaseField DwDy;
        SinglePhaseField DwDz;

        int NoOfEdges;
        string[] edgeTagNames;

        EdgeIntegral[] XForceIntegral;
        EdgeIntegral[] YForceIntegral;
        EdgeIntegral[] ZForceIntegral;

        EdgeIntegral[] XMomentIntegral;
        EdgeIntegral[] YMomentIntegral;
        EdgeIntegral[] ZMomentIntegral;
        
        /// <summary>
        /// Ctor.
        /// </summary>
        /// <param name="_app"></param>
        /// <param name="_edgeTagNames">
        /// Edgetags, which should be considered for calculating the forces.
        /// </param>
        public Force(NSE_SIMPLEMain _app, string[] _edgeTagNames) {
            m_app = _app;
            CreateFields();

            NoOfEdges = _edgeTagNames.Length;
            edgeTagNames = _edgeTagNames;

            XForceIntegral = new EdgeIntegral[NoOfEdges];
            YForceIntegral = new EdgeIntegral[NoOfEdges];
            ZMomentIntegral = new EdgeIntegral[NoOfEdges];


            if (m_app.GridData.SpatialDimension == 3) {
                ZForceIntegral = new EdgeIntegral[NoOfEdges];
                XMomentIntegral = new EdgeIntegral[NoOfEdges];
                YMomentIntegral = new EdgeIntegral[NoOfEdges];
            }

            InitEdgeIntegrals();
        }

        void InitEdgeIntegrals() {
            if (m_app.GridData.SpatialDimension == 2) {
                for (int edge = 0; edge < NoOfEdges; edge++) {
                    CoordinateMapping m_CoordinateMapping = new CoordinateMapping(m_app.WorkingSet.Pressure, DuDx, DuDy, DvDx, DvDy);
                    var QuadratureOrder = m_CoordinateMapping.BasisS.Max(basis => basis.Degree) * 2 + 1;
                    XForceIntegral[edge] = new EdgeIntegral(m_app.GridData,
                        edgeTagNames[edge],
                        new ForceFlux2D(0, m_app),
                        m_CoordinateMapping,
                        QuadratureOrder);

                    YForceIntegral[edge] = new EdgeIntegral(m_app.GridData,
                        edgeTagNames[edge],
                        new ForceFlux2D(1, m_app),
                        m_CoordinateMapping,
                        QuadratureOrder);

                    ZMomentIntegral[edge] = new EdgeIntegral(m_app.GridData,
                        edgeTagNames[edge],
                        new MomentFlux2D(m_app),
                        m_CoordinateMapping,
                        QuadratureOrder);
                }
            } else {
                for (int edge = 0; edge < NoOfEdges; edge++) {
                    var m_CoordMap = new CoordinateMapping(m_app.WorkingSet.Pressure, DuDx, DuDy, DuDz, DvDx, DvDy, DvDz, DwDx, DwDy, DwDz);
                    var QuadOrder = m_CoordMap.BasisS.Max(basis => basis.Degree) * 2 + 1;
                    XForceIntegral[edge] = new EdgeIntegral(m_app.GridData,
                        edgeTagNames[edge],
                        new ForceFlux3D(0, m_app),
                        m_CoordMap,
                        QuadOrder);

                    YForceIntegral[edge] = new EdgeIntegral(m_app.GridData,
                        edgeTagNames[edge],
                        new ForceFlux3D(1, m_app),
                        m_CoordMap,
                        QuadOrder);

                    ZForceIntegral[edge] = new EdgeIntegral(m_app.GridData,
                        edgeTagNames[edge],
                        new ForceFlux3D(2, m_app),
                        m_CoordMap,
                        QuadOrder);

                    XMomentIntegral[edge] = new EdgeIntegral(m_app.GridData,
                        edgeTagNames[edge],
                        new MomentFlux3D(0, m_app),
                        m_CoordMap,
                        QuadOrder);

                    YMomentIntegral[edge] = new EdgeIntegral(m_app.GridData,
                        edgeTagNames[edge],
                        new MomentFlux3D(1, m_app),
                        m_CoordMap,
                        QuadOrder);

                    ZMomentIntegral[edge] = new EdgeIntegral(m_app.GridData,
                        edgeTagNames[edge],
                        new MomentFlux3D(2, m_app),
                        m_CoordMap,
                        QuadOrder);

                }
            }
        }
     
        void CreateFields() {
            int DGOrder = m_app.WorkingSet.VelBasis.Degree;
            Basis DGBasis = new Basis(m_app.GridData, DGOrder);

            DuDx = new SinglePhaseField(DGBasis);
            DuDy = new SinglePhaseField(DGBasis);
            DvDx = new SinglePhaseField(DGBasis);
            DvDy = new SinglePhaseField(DGBasis);

            if (m_app.GridData.SpatialDimension == 3) {
                DuDz = new SinglePhaseField(DGBasis);
                DvDz = new SinglePhaseField(DGBasis);
                DwDx = new SinglePhaseField(DGBasis);
                DwDy = new SinglePhaseField(DGBasis);
                DwDz = new SinglePhaseField(DGBasis);
            }
        }

        void CalculateDerivative() {
            DuDx.Clear();
            DuDy.Clear();
            DvDx.Clear();
            DvDy.Clear();

            DuDx.Derivative(1.0, m_app.WorkingSet.Velocity.Current[0], 0);
            DuDy.Derivative(1.0, m_app.WorkingSet.Velocity.Current[0], 1);
            DvDx.Derivative(1.0, m_app.WorkingSet.Velocity.Current[1], 0);
            DvDy.Derivative(1.0, m_app.WorkingSet.Velocity.Current[1], 1);

            if (m_app.GridData.SpatialDimension == 3) {
                DuDz.Clear();
                DvDz.Clear();
                DwDx.Clear();
                DwDy.Clear();
                DwDz.Clear();

                DuDz.Derivative(1.0, m_app.WorkingSet.Velocity.Current[0], 2);
                DvDz.Derivative(1.0, m_app.WorkingSet.Velocity.Current[1], 2);
                DwDx.Derivative(1.0, m_app.WorkingSet.Velocity.Current[2], 0);
                DwDy.Derivative(1.0, m_app.WorkingSet.Velocity.Current[2], 1);
                DwDz.Derivative(1.0, m_app.WorkingSet.Velocity.Current[2], 2);
            }
        }

        void CalculateDerivativeByFlux() {
            DuDx.Clear();
            DuDy.Clear();
            DvDx.Clear();
            DvDy.Clear();

            DuDx.DerivativeByFlux(1.0, m_app.WorkingSet.Velocity.Current[0], 0);
            DuDy.DerivativeByFlux(1.0, m_app.WorkingSet.Velocity.Current[0], 1);
            DvDx.DerivativeByFlux(1.0, m_app.WorkingSet.Velocity.Current[1], 0);
            DvDy.DerivativeByFlux(1.0, m_app.WorkingSet.Velocity.Current[1], 1);

            if (m_app.GridData.SpatialDimension == 3) {
                DuDz.Clear();
                DvDz.Clear();
                DwDx.Clear();
                DwDy.Clear();
                DwDz.Clear();

                DuDz.DerivativeByFlux(1.0, m_app.WorkingSet.Velocity.Current[0], 2);
                DvDz.DerivativeByFlux(1.0, m_app.WorkingSet.Velocity.Current[1], 2);
                DwDx.DerivativeByFlux(1.0, m_app.WorkingSet.Velocity.Current[2], 0);
                DwDy.DerivativeByFlux(1.0, m_app.WorkingSet.Velocity.Current[2], 1);
                DwDz.DerivativeByFlux(1.0, m_app.WorkingSet.Velocity.Current[2], 2);
            }
        }

        #region Property Encapsulation
        double m_XForce;

        /// <summary>
        /// Drag force.
        /// </summary>
        public double XForce {
            get { 
                return m_XForce;
            }
        }

        double m_YForce;

        /// <summary>
        /// Lift force.
        /// </summary>
        public double YForce {
            get {
                return m_YForce;
            }
        }

        double m_ZForce;

        /// <summary>
        /// Lateral force.
        /// </summary>
        public double ZForce {
            get {
                return m_ZForce;
            }
        }

        double m_XMoment;

        /// <summary>
        /// Drag Moment.
        /// </summary>
        public double XMoment {
            get {
                return m_XMoment;
            }
        }

        double m_YMoment;

        /// <summary>
        /// Lift Moment.
        /// </summary>
        public double YMoment {
            get {
                return m_YMoment;
            }
        }

        double m_ZMoment;

        /// <summary>
        /// Lateral Moment.
        /// </summary>
        public double ZMoment {
            get {
                return m_ZMoment;
            }
        }
        #endregion

        public void CalculateForce( ) {
            m_XForce = 0.0;
            m_YForce = 0.0;
            m_ZForce = 0.0;
            m_XMoment = 0.0;
            m_YMoment = 0.0;
            m_ZMoment = 0.0;
            
            //CalculateDerivative();
            //Console.WriteLine("Warning: Derivatives Calculated by Flux instead of Symbolic - Symbolic not implemented for Curved Elements");
            CalculateDerivativeByFlux();


            for (int edge = 0; edge < NoOfEdges; edge++) {

                double LocalXForce = 0.0;
                double GlobalXForce = double.NaN;

                double LocalYForce = 0.0;
                double GlobalYForce = double.NaN;

                double LocalZForce = 0.0;
                double GlobalZForce = double.NaN;

                double LocalXMoment = 0.0;
                double GlobalXMoment = double.NaN;

                double LocalYMoment = 0.0;
                double GlobalYMoment = double.NaN;

                double LocalZMoment = 0.0;
                double GlobalZMoment = double.NaN;


                LocalXForce = XForceIntegral[edge].Evaluate();                
                LocalYForce = YForceIntegral[edge].Evaluate();
                LocalZMoment = ZMomentIntegral[edge].Evaluate();

                unsafe {
                    MPI.Wrappers.csMPI.Raw.Allreduce((IntPtr)(&LocalXForce), (IntPtr)(&GlobalXForce), 1, MPI.Wrappers.csMPI.Raw._DATATYPE.DOUBLE, MPI.Wrappers.csMPI.Raw._OP.SUM, MPI.Wrappers.csMPI.Raw._COMM.WORLD);
                    MPI.Wrappers.csMPI.Raw.Allreduce((IntPtr)(&LocalYForce), (IntPtr)(&GlobalYForce), 1, MPI.Wrappers.csMPI.Raw._DATATYPE.DOUBLE, MPI.Wrappers.csMPI.Raw._OP.SUM, MPI.Wrappers.csMPI.Raw._COMM.WORLD);
                    MPI.Wrappers.csMPI.Raw.Allreduce((IntPtr)(&LocalZMoment), (IntPtr)(&GlobalZMoment), 1, MPI.Wrappers.csMPI.Raw._DATATYPE.DOUBLE, MPI.Wrappers.csMPI.Raw._OP.SUM, MPI.Wrappers.csMPI.Raw._COMM.WORLD);
                }

                m_XForce += GlobalXForce;
                m_YForce += GlobalYForce;
                m_ZMoment += GlobalZMoment;

                if (m_app.GridData.SpatialDimension == 3) {
                    LocalZForce = ZForceIntegral[edge].Evaluate();
                    LocalXMoment = XMomentIntegral[edge].Evaluate();
                    LocalYMoment = YMomentIntegral[edge].Evaluate();


                    unsafe {
                        MPI.Wrappers.csMPI.Raw.Allreduce((IntPtr)(&LocalZForce), (IntPtr)(&GlobalZForce), 1, MPI.Wrappers.csMPI.Raw._DATATYPE.DOUBLE, MPI.Wrappers.csMPI.Raw._OP.SUM, MPI.Wrappers.csMPI.Raw._COMM.WORLD);
                        MPI.Wrappers.csMPI.Raw.Allreduce((IntPtr)(&LocalXMoment), (IntPtr)(&GlobalXMoment), 1, MPI.Wrappers.csMPI.Raw._DATATYPE.DOUBLE, MPI.Wrappers.csMPI.Raw._OP.SUM, MPI.Wrappers.csMPI.Raw._COMM.WORLD);
                        MPI.Wrappers.csMPI.Raw.Allreduce((IntPtr)(&LocalYMoment), (IntPtr)(&GlobalYMoment), 1, MPI.Wrappers.csMPI.Raw._DATATYPE.DOUBLE, MPI.Wrappers.csMPI.Raw._OP.SUM, MPI.Wrappers.csMPI.Raw._COMM.WORLD);
                    }

                    m_ZForce  += GlobalZForce;
                    m_XMoment += GlobalXMoment;
                    m_YMoment += GlobalYMoment;
                }
            }
        }
    }

    /// <summary>
    /// Forces in 2D.
    /// </summary>
    class ForceFlux2D : EdgeIntegral.EdgeFlux {        

        NSE_SIMPLEMain m_app;
        int m_Component;
        double m_Reynolds;

        public ForceFlux2D(int _Component, NSE_SIMPLEMain _app) {
            m_app = _app;
            m_Component = _Component;            
            m_Reynolds = m_app.SolverConf.Control.Reynolds;
        }

        protected override double BorderEdgeFlux(double time, double[] x, double[] normal, byte EdgeTag, double[] Uin, int jEdge) {
            switch (m_Component) {
                case 0:
                    return -Uin[0] * normal[0] + (2.0 / m_Reynolds) * Uin[1] * normal[0]
                                               + (1.0 / m_Reynolds) * (Uin[3] + Uin[2]) * normal[1];
                case 1:
                    return -Uin[0] * normal[1] + (2.0 / m_Reynolds) * Uin[4] * normal[1]
                                               + (1.0 / m_Reynolds) * (Uin[3] + Uin[2]) * normal[0];
                default:
                    throw new ArgumentException();
            }            
        }

        string[] m_ArgumentOrdering = new string[] { "Pressure", "DuDx", "DuDy", "DvDx", "DvDy" };

        public override IList<string> ArgumentOrdering {
            get { return m_ArgumentOrdering; }
        }
    }






    /// <summary>
    /// Moment about zAxis in 2D.
    /// </summary>
    class MomentFlux2D : EdgeIntegral.EdgeFlux {

        NSE_SIMPLEMain m_app;
        double m_Reynolds;

        public MomentFlux2D( NSE_SIMPLEMain _app) {
            m_app = _app;
            m_Reynolds = m_app.SolverConf.Control.Reynolds;
        }

        protected override double BorderEdgeFlux(double time, double[] x, double[] normal, byte EdgeTag, double[] Uin, int jEdge) {
            
                double xForce =  -Uin[0] * normal[0] + (2.0 / m_Reynolds) * Uin[1] * normal[0]
                                               + (1.0 / m_Reynolds) * (Uin[3] + Uin[2]) * normal[1];
                double yForce = -Uin[0] * normal[1] + (2.0 / m_Reynolds) * Uin[4] * normal[1]
                                               + (1.0 / m_Reynolds) * (Uin[3] + Uin[2]) * normal[0];
                return x[0]*yForce - x[1]*xForce;
        }

        string[] m_ArgumentOrdering = new string[] { "Pressure", "DuDx", "DuDy", "DvDx", "DvDy" };

        public override IList<string> ArgumentOrdering {
            get { return m_ArgumentOrdering; }
        }
    }


    /// <summary>
    /// Forces in 3D.
    /// </summary>
    class ForceFlux3D : EdgeIntegral.EdgeFlux {        

        NSE_SIMPLEMain m_app;
        int m_Component;
        double m_Reynolds;

        public ForceFlux3D(int _Component, NSE_SIMPLEMain _app) {
            m_app = _app;
            m_Component = _Component;            
            m_Reynolds = m_app.SolverConf.Control.Reynolds;
        }

        protected override double BorderEdgeFlux(double time, double[] x, double[] normal, byte EdgeTag, double[] Uin, int jEdge) {
            switch (m_Component) {
                case 0:
                    return -Uin[0] * normal[0] + (2.0 / m_Reynolds) * Uin[1] * normal[0]
                                               + (1.0 / m_Reynolds) * (Uin[2] + Uin[4]) * normal[1]
                                               + (1.0 / m_Reynolds) * (Uin[3] + Uin[7]) * normal[2];
                case 1:
                    return -Uin[0] * normal[1] + (2.0 / m_Reynolds) * Uin[5] * normal[1]
                                               + (1.0 / m_Reynolds) * (Uin[6] + Uin[8]) * normal[2]
                                               + (1.0 / m_Reynolds) * (Uin[4] + Uin[2]) * normal[0];
                case 2:
                    return -Uin[0] * normal[2] + (2.0 / m_Reynolds) * Uin[9] * normal[2]
                                               + (1.0 / m_Reynolds) * (Uin[7] + Uin[3]) * normal[0]
                                               + (1.0 / m_Reynolds) * (Uin[8] + Uin[6]) * normal[1];
                default:
                    throw new ArgumentException();
            }            
        }

        string[] m_ArgumentOrdering = new string[] { "Pressure", "DuDx", "DuDy", "DuDz", "DvDx", "DvDy", "DvDz", "DwDx", "DwDy", "DwDz" };

        public override IList<string> ArgumentOrdering {
            get { return m_ArgumentOrdering; }
        }
    }


    /// <summary>
    /// Moments in 3D.
    /// </summary>
    class MomentFlux3D : EdgeIntegral.EdgeFlux {

        NSE_SIMPLEMain m_app;
        int m_Component;
        double m_Reynolds;

        public MomentFlux3D(int _Component, NSE_SIMPLEMain _app) {
            m_app = _app;
            m_Component = _Component;
            m_Reynolds = m_app.SolverConf.Control.Reynolds;
        }

        protected override double BorderEdgeFlux(double time, double[] x, double[] normal, byte EdgeTag, double[] Uin, int jEdge) {
            double xForce = -Uin[0] * normal[0] + (2.0 / m_Reynolds) * Uin[1] * normal[0]
                                            + (1.0 / m_Reynolds) * (Uin[2] + Uin[4]) * normal[1]
                                            + (1.0 / m_Reynolds) * (Uin[3] + Uin[7]) * normal[2];

            double yForce = -Uin[0] * normal[1] + (2.0 / m_Reynolds) * Uin[5] * normal[1]
                                           + (1.0 / m_Reynolds) * (Uin[6] + Uin[8]) * normal[2]
                                           + (1.0 / m_Reynolds) * (Uin[4] + Uin[2]) * normal[0];

            double zForce = -Uin[0] * normal[2] + (2.0 / m_Reynolds) * Uin[9] * normal[2]
                                           + (1.0 / m_Reynolds) * (Uin[7] + Uin[3]) * normal[0]
                                           + (1.0 / m_Reynolds) * (Uin[8] + Uin[6]) * normal[1];

            switch (m_Component) {
                case 0:
                    return x[1] * zForce - x[2] * yForce;

                case 1:
                    return x[2] * xForce - x[0] * zForce;
                
                case 2:               
                        return x[0]*yForce - x[1]*xForce;
                default:
                    throw new ArgumentException();
            }
        }

        string[] m_ArgumentOrdering = new string[] { "Pressure", "DuDx", "DuDy", "DuDz", "DvDx", "DvDy", "DvDz", "DwDx", "DwDy", "DwDz" };

        public override IList<string> ArgumentOrdering {
            get { return m_ArgumentOrdering; }
        }
    }






















}