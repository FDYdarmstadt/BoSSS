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
using System.Linq;
using System.Diagnostics;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using ilPSP.Tracing;

namespace BoSSS.Foundation.Grid.RefElements {


    /// <summary>
    /// The line reference element, i.e. \f$ K^{\textrm{line}} = ( -1,1 ) \f$.
    /// </summary>
    public class Line : RefElement {


        /// <summary>
        /// Indicates the "faces" -1 and 1 .
        /// </summary>
        public enum Edge {
            
            /// <summary>
            /// -1
            /// </summary>
            Left = 0,

            /// <summary>
            /// 1
            /// </summary>
            Right = 1

        }

        ///// <summary>
        ///// Nodes for volume (internal) quadrature in local coordinates;
        ///// 1st index: quadrature order; 
        ///// </summary>
        ///// <see cref="m_QuadratureWeights"/>
        //public double[][,] m_QuadratureNodes;

        ///// <summary>
        ///// Weights for  quadrature in local coordinates;
        ///// 1st index: quadrature order; 2nd index: node index; 
        ///// </summary>
        ///// <see cref="m_QuadratureNodes"/>
        //public double[][] m_QuadratureWeights;



        private static Line instance = null;
        private static readonly object padlock = new object();

        /// <summary>
        /// Access to the single, global instance.
        /// </summary>
        public static Line Instance {
            get {
                lock(padlock) {
                    if(instance == null) {
                        instance = new Line();
                    }
                    return instance;
                }
            }
        }



        /// <summary>
        /// default constructor
        /// </summary>
        private Line() {
            using (new FuncTrace()) { 
            // ===============
            // define vertices
            // ===============

            var _Vertices = new double[2, 1] { { -1 }, { 1 } };
            this.m_Vertices = new NodeSet(this, 2, 1);
            this.m_Vertices.InitializeFrom(_Vertices);
            this.m_Vertices.LockForever();

            m_NoOfFaces = 2;

            // ============
            // edge simplex
            // ============

            m_FaceRefElement = Point.Instance;


            // ===================================
            // define Quadrature nodes and weights
            // ===================================

            //m_QuadratureNodes   = new double[10][,];
            //m_QuadratureWeights = new double[10][];

            {
                var qrTemp = QuadRuleResource.DecodeFromBase64(Resource.LineQuadRules_bin);
                foreach (var q in qrTemp) {
                    var realQr = QuadRule.CreateEmpty(this, q.Item2.GetLength(0), this.SpatialDimension);
                    realQr.Nodes.Set2DArray(q.Item2);
                    realQr.Weights.SetVector(q.Item3);
                    realQr.OrderOfPrecision = q.Item1;
                    realQr.Nodes.LockForever();
                    realQr.Weights.LockForever();
                    base.m_QuadRules.Add(realQr);
                }
            }

            /*
            foreach (int NoOfNodes in QuadratureRules1D.AvailableNoOfNodes) {

                // quadrature rule
                // ---------------


                QuadratureRules1D gq = new QuadratureRules1D(NoOfNodes);


                // setup integral nodes
                // --------------------

                QuadRule qr = QuadRule.CreateEmpty(this, NoOfNodes, 1);
                var VolumeNodes = qr.Nodes;
                var VolumeWeights = qr.Weights;

                for (int i = 0; i < NoOfNodes; i++) {

                    int ino = i - NoOfNodes / 2;
                    double isng = 1.0;
                    if (ino < 0)
                        isng = -1.0;
                    if (ino < 0) {
                        ino *= -1;
                        if (NoOfNodes % 2 == 0)
                            ino--;
                    }

                    VolumeNodes[i, 0] = gq.m_Nodes[ino] * isng;
                    VolumeWeights[i] = gq.m_Weights[ino];


                }

                qr.OrderOfPrecision = gq.precision;
                qr.Nodes.LockForever();
                m_QuadRules.Add(qr);

            }
            */


            // ==============================
            // define orthonormal polynomials
            // ==============================
#pragma warning disable 612
            #region POLYDEF
            OrthonormalPolynomials = new Polynomial[21];
            Polynomial p;

            p = new Polynomial(new Guid("{D74D00DA-D680-43bd-827D-7A404505CDBA}"));
            OrthonormalPolynomials[0] = p;
            p.AddCoeff(7.07106781186547524401e-01, new int[] { 0 });

            p = new Polynomial(new Guid("{CDCCFAC6-6AEF-42fb-BDF0-9915FFC6D982}"));
            OrthonormalPolynomials[1] = p;
            p.AddCoeff(1.22474487139158904910e+00, new int[] { 1 });

            p = new Polynomial(new Guid("{1578163B-6EF8-45ce-B373-841C03B77B51}"));
            OrthonormalPolynomials[2] = p;
            p.AddCoeff(-7.90569415042094833000e-01, new int[] { 0 });
            p.AddCoeff(2.37170824512628449900e+00, new int[] { 2 });

            p = new Polynomial(new Guid("{A87ACEDB-BD12-49f2-B493-89AA9B868B9E}"));
            OrthonormalPolynomials[3] = p;
            p.AddCoeff(-2.80624304008045603919e+00, new int[] { 1 });
            p.AddCoeff(4.67707173346742673198e+00, new int[] { 3 });

            p = new Polynomial(new Guid("{0F6B1BB8-8E0D-4610-80D1-72CC27A95CD5}"));
            OrthonormalPolynomials[4] = p;
            p.AddCoeff(7.95495128834865964951e-01, new int[] { 0 });
            p.AddCoeff(-7.95495128834865964951e+00, new int[] { 2 });
            p.AddCoeff(9.28077650307343625776e+00, new int[] { 4 });

            p = new Polynomial(new Guid("{E18976B3-677A-4b53-95BD-063AA8E11CE5}"));
            OrthonormalPolynomials[5] = p;
            p.AddCoeff(4.39726477483446520741e+00, new int[] { 1 });
            p.AddCoeff(-2.05205689492275043012e+01, new int[] { 3 });
            p.AddCoeff(1.84685120543047538711e+01, new int[] { 5 });

            p = new Polynomial(new Guid("{94C8333B-68CF-4936-A3D4-F2AB27E5E814}"));
            OrthonormalPolynomials[6] = p;
            p.AddCoeff(-7.96721798998872629692e-01, new int[] { 0 });
            p.AddCoeff(1.67311577789763252235e+01, new int[] { 2 });
            p.AddCoeff(-5.01934733369289756706e+01, new int[] { 4 });
            p.AddCoeff(3.68085471137479154918e+01, new int[] { 6 });

            p = new Polynomial(new Guid("{9508DB63-E7A9-485b-B05A-89D43FF23BFE}"));
            OrthonormalPolynomials[7] = p;
            p.AddCoeff(-5.99071547271275436594e+00, new int[] { 1 });
            p.AddCoeff(5.39164392544147892934e+01, new int[] { 3 });
            p.AddCoeff(-1.18616166359712536446e+02, new int[] { 5 });
            p.AddCoeff(7.34290553655363320853e+01, new int[] { 7 });

            p = new Polynomial(new Guid("{D838BD06-CF79-4f2c-8ECC-1F5F0EA339D5}"));
            OrthonormalPolynomials[8] = p;
            p.AddCoeff(7.97200454373380923752e-01, new int[] { 0 });
            p.AddCoeff(-2.86992163574417132551e+01, new int[] { 2 });
            p.AddCoeff(1.57845689965929422903e+02, new int[] { 4 });
            p.AddCoeff(-2.73599195940944333032e+02, new int[] { 6 });
            p.AddCoeff(1.46570997825505892696e+02, new int[] { 8 });

            p = new Polynomial(new Guid("{141C0FAA-C58B-4516-9935-3C56FA62D4D9}"));
            OrthonormalPolynomials[9] = p;
            p.AddCoeff(7.58511879271573274152e+00, new int[] { 1 });
            p.AddCoeff(-1.11248408959830746876e+02, new int[] { 3 });
            p.AddCoeff(4.33868794943339912815e+02, new int[] { 5 });
            p.AddCoeff(-6.19812564204771304021e+02, new int[] { 7 });
            p.AddCoeff(2.92689266430030893566e+02, new int[] { 9 });

            p = new Polynomial(new Guid("{6C16E97A-6A32-456d-841D-CA98681F6F5C}"));
            OrthonormalPolynomials[10] = p;
            p.AddCoeff(-7.97434890624404676857e-01, new int[] { 0 });
            p.AddCoeff(4.38589189843422572271e+01, new int[] { 2 });
            p.AddCoeff(-3.80110631197632895969e+02, new int[] { 4 });
            p.AddCoeff(1.14033189359289868791e+03, new int[] { 6 });
            p.AddCoeff(-1.38468872793423412103e+03, new int[] { 8 });
            p.AddCoeff(5.84646351794454406656e+02, new int[] { 10 });




            p = new Polynomial(new Guid("{3C81368B-C77A-4F38-9AFA-50BB43DA4E2F}"));
            OrthonormalPolynomials[11] = p;
            p.AddCoeff(-9.17998960606603675854e+00, new int[] { 1 });
            p.AddCoeff(1.98899774798097463102e+02, new int[] { 3 });
            p.AddCoeff(-1.19339864878858477861e+03, new int[] { 5 });
            p.AddCoeff(2.89825386134370589091e+03, new int[] { 7 });
            p.AddCoeff(-3.05926796475168955152e+03, new int[] { 9 });
            p.AddCoeff(1.16808413199609964694e+03, new int[] { 11 });

            p = new Polynomial(new Guid("{38257CD0-8D4E-41B6-8BA5-585BA6EE5A6A}"));
            OrthonormalPolynomials[12] = p;
            p.AddCoeff(7.97566730732873428401e-01, new int[] { 0 });
            p.AddCoeff(-6.22102049971641274153e+01, new int[] { 2 });
            p.AddCoeff(7.77627562464551592691e+02, new int[] { 4 });
            p.AddCoeff(-3.52524494983930055353e+03, new int[] { 6 });
            p.AddCoeff(7.17639150503000469827e+03, new int[] { 8 });
            p.AddCoeff(-6.69796540469467105171e+03, new int[] { 10 });
            p.AddCoeff(2.33413945921177930590e+03, new int[] { 12 });

            p = new Polynomial(new Guid("{D5FB0641-6C99-4196-8EE3-9057B5A3549D}"));
            OrthonormalPolynomials[13] = p;
            p.AddCoeff(1.07751235804363532650e+01, new int[] { 1 });
            p.AddCoeff(-3.23253707413090597949e+02, new int[] { 3 });
            p.AddCoeff(2.74765651301127008257e+03, new int[] { 5 });
            p.AddCoeff(-9.94389976137412029882e+03, new int[] { 7 });
            p.AddCoeff(1.74018245824047105229e+04, new int[] { 9 });
            p.AddCoeff(-1.45542532871021215283e+04, new int[] { 11 });
            p.AddCoeff(4.66482477150709023342e+03, new int[] { 13 });

            p = new Polynomial(new Guid("{805C0075-690E-44B9-8C7C-5CBE4D3A6F9E}"));
            OrthonormalPolynomials[14] = p;
            p.AddCoeff(-7.97648110941312659802e-01, new int[] { 0 });
            p.AddCoeff(8.37530516488378292792e+01, new int[] { 2 });
            p.AddCoeff(-1.42380187803024309775e+03, new int[] { 4 });
            p.AddCoeff(9.01741189419153961906e+03, new int[] { 6 });
            p.AddCoeff(-2.70522356825746188572e+04, new int[] { 8 });
            p.AddCoeff(4.14800947132810822477e+04, new int[] { 10 });
            p.AddCoeff(-3.14243141767280926119e+04, new int[] { 12 });
            p.AddCoeff(9.32369761287536813759e+03, new int[] { 14 });

            p = new Polynomial(new Guid("{29D06AEC-8F47-4487-96DC-343DF735A531}"));
            OrthonormalPolynomials[15] = p;
            p.AddCoeff(-1.23704200852700204862e+01, new int[] { 1 });
            p.AddCoeff(4.90693330049044145955e+02, new int[] { 3 });
            p.AddCoeff(-5.59390396255910326388e+03, new int[] { 5 });
            p.AddCoeff(2.79695198127955163194e+04, new int[] { 7 });
            p.AddCoeff(-7.14776617438107639274e+04, new int[] { 9 });
            p.AddCoeff(9.74695387415601326282e+04, new int[] { 11 });
            p.AddCoeff(-6.74789114364647072042e+04, new int[] { 13 });
            p.AddCoeff(1.86370326824521572278e+04, new int[] { 15 });

            p = new Polynomial(new Guid("{8DB3325C-CAFA-4A43-A43C-BB60FC025B9E}"));
            OrthonormalPolynomials[16] = p;
            p.AddCoeff(7.97701830045050123895e-01, new int[] { 0 });
            p.AddCoeff(-1.08487448886126816850e+02, new int[] { 2 });
            p.AddCoeff(2.40480511697581110683e+03, new int[] { 4 });
            p.AddCoeff(-2.02003629825968132974e+04, new int[] { 6 });
            p.AddCoeff(8.29657765356654831858e+04, new int[] { 8 });
            p.AddCoeff(-1.84368392301478851524e+05, new int[] { 10 });
            p.AddCoeff(2.26270299642724045052e+05, new int[] { 12 });
            p.AddCoeff(-1.44216234937120819923e+05, new int[] { 14 });
            p.AddCoeff(3.72558606920895451469e+04, new int[] { 16 });

            p = new Polynomial(new Guid("{A3333639-DF1A-4217-AE29-F988AC4A5C46}"));
            OrthonormalPolynomials[17] = p;
            p.AddCoeff(1.39658239139854728044e+01, new int[] { 1 });
            p.AddCoeff(-7.07601744975263955424e+02, new int[] { 3 });
            p.AddCoeff(1.04017456511363801447e+04, new int[] { 5 });
            p.AddCoeff(-6.83543285646104980940e+04, new int[] { 7 });
            p.AddCoeff(2.37341418627119785049e+05, new int[] { 9 });
            p.AddCoeff(-4.66052240213253396095e+05, new int[] { 11 });
            p.AddCoeff(5.19827498699398018722e+05, new int[] { 13 });
            p.AddCoeff(-3.06945761136787401531e+05, new int[] { 15 });
            p.AddCoeff(7.44794861581910606656e+04, new int[] { 17 });

            p = new Polynomial(new Guid("{F83B3BFF-6FB2-48B1-A131-B09331C8CBF3}"));
            OrthonormalPolynomials[18] = p;
            p.AddCoeff(-7.97739132849907901080e-01, new int[] { 0 });
            p.AddCoeff(1.36413391717334251085e+02, new int[] { 2 });
            p.AddCoeff(-3.81957496808535903037e+03, new int[] { 4 });
            p.AddCoeff(4.09967713241161869260e+04, new int[] { 6 });
            p.AddCoeff(-2.19625560664908144246e+05, new int[] { 8 });
            p.AddCoeff(6.58876681994724432739e+05, new int[] { 10 });
            p.AddCoeff(-1.15802568350587930603e+06, new int[] { 12 });
            p.AddCoeff(1.18347679742908544462e+06, new int[] { 14 });
            p.AddCoeff(-6.50912238585996994541e+05, new int[] { 16 });
            p.AddCoeff(1.48901492486992776529e+05, new int[] { 18 });

            p = new Polynomial(new Guid("{F79EE3E8-13E6-413B-9F46-D81AC135BE58}"));
            OrthonormalPolynomials[19] = p;
            p.AddCoeff(-1.55613022863318221426e+01, new int[] { 1 });
            p.AddCoeff(9.80362044038904794983e+02, new int[] { 3 });
            p.AddCoeff(-1.80386616103158482277e+04, new int[] { 5 });
            p.AddCoeff(1.50322180085965401897e+05, new int[] { 7 });
            p.AddCoeff(-6.76449810386844308538e+05, new int[] { 9 });
            p.AddCoeff(1.78336768192895317706e+06, new int[] { 11 });
            p.AddCoeff(-2.83509734050243838404e+06, new int[] { 13 });
            p.AddCoeff(2.67309177818801333352e+06, new int[] { 15 });
            p.AddCoeff(-1.37585606230265392167e+06, new int[] { 17 });
            p.AddCoeff(2.97699849738001140945e+05, new int[] { 19 });

            p = new Polynomial(new Guid("{5F643FB3-EF93-4F4F-8F12-ACA609DE66ED}"));
            OrthonormalPolynomials[20] = p;
            p.AddCoeff(7.97766083041055939798e-01, new int[] { 0 });
            p.AddCoeff(-1.67530877438621747358e+02, new int[] { 2 });
            p.AddCoeff(5.77981527163245028384e+03, new int[] { 4 });
            p.AddCoeff(-7.70642036217660037845e+04, new int[] { 6 });
            p.AddCoeff(5.20183374446920525545e+05, new int[] { 8 });
            p.AddCoeff(-2.01137571452809269878e+06, new int[] { 10 });
            p.AddCoeff(4.72368539017961164106e+06, new int[] { 12 });
            p.AddCoeff(-6.85193924729350260022e+06, new int[] { 14 });
            p.AddCoeff(5.99544684138181477520e+06, new int[] { 16 });
            p.AddCoeff(-2.89975860302126989127e+06, new int[] { 18 });
            p.AddCoeff(5.95213607988576451366e+05, new int[] { 20 });
#pragma warning restore 612
            #endregion POLYDEF
        }
        }

        /*
        #region QUADRULEDEF

        /// <summary>
        /// A collection of
        /// gauss quadrature rules for the intervall [-1,1]
        /// from different sources
        /// </summary>
        internal class QuadratureRules1D {

            /// <summary>
            /// nodes in one half of the domain (nodes are symmetric), in the intervall [0,1];
            /// nodes are sorted ascending;
            /// </summary>
            internal double[] m_Nodes;

            /// <summary>
            /// Weights for quadrature nodes
            /// </summary>
            internal double[] m_Weights;


            /// <summary>
            /// Number of nodes: either (Length of <see cref="m_Nodes"/>)*2, or (Length of <see cref="m_Nodes"/>)*2 - 1,
            /// depending wether <see cref="m_Nodes"/> contains 0 or not.
            /// </summary>
            internal int NoOfNodes;

            /// <summary>
            /// precision of quadrature rule: polynomials up to this degree are integrated exact, up to
            /// rounding errors.
            /// </summary>
            internal int precision;


            /// <summary>
            /// A collection of all available number of nodes
            /// </summary>
            public static int[] AvailableNoOfNodes = Enumerable.Range(1, 60).ToArray();

            /// <summary>
            /// A collection of 1D gaussian quadrature rules created using
            /// Florian's GaussRule.mw
            /// </summary>
            /// <param name="__NoOfNodes">
            /// The number of evaluation points to be evaluated
            /// </param>
            internal QuadratureRules1D(int __NoOfNodes) {

                m_Nodes = new double[__NoOfNodes / 2 + __NoOfNodes % 2];
                m_Weights = (double[])m_Nodes.Clone();

                NoOfNodes = __NoOfNodes;

                precision = 2 * __NoOfNodes - 1;

                //Gauss rules created using Florian's GaussRule.mw
                switch (__NoOfNodes) {
                    case 1:
                        m_Nodes[0] = 0.0;
                        m_Weights[0] = 2.00000000000000000000e+00;
                        break;

                    case 2:
                        m_Nodes[0] = 5.77350269189625764509e-01;
                        m_Weights[0] = 1.00000000000000000000e+00;
                        break;

                    case 3:
                        m_Nodes[0] = 0.0;
                        m_Nodes[1] = 7.74596669241483377036e-01;
                        m_Weights[0] = 8.88888888888888888889e-01;
                        m_Weights[1] = 5.55555555555555555556e-01;
                        break;

                    case 4:
                        m_Nodes[0] = 3.39981043584856264803e-01;
                        m_Nodes[1] = 8.61136311594052575224e-01;
                        m_Weights[0] = 6.52145154862546142627e-01;
                        m_Weights[1] = 3.47854845137453857373e-01;
                        break;

                    case 5:
                        m_Nodes[0] = 0.0;
                        m_Nodes[1] = 5.38469310105683091036e-01;
                        m_Nodes[2] = 9.06179845938663992798e-01;
                        m_Weights[0] = 5.68888888888888888889e-01;
                        m_Weights[1] = 4.78628670499366468041e-01;
                        m_Weights[2] = 2.36926885056189087514e-01;
                        break;

                    case 6:
                        m_Nodes[0] = 2.38619186083196908631e-01;
                        m_Nodes[1] = 6.61209386466264513661e-01;
                        m_Nodes[2] = 9.32469514203152027812e-01;
                        m_Weights[0] = 4.67913934572691047390e-01;
                        m_Weights[1] = 3.60761573048138607570e-01;
                        m_Weights[2] = 1.71324492379170345040e-01;
                        break;

                    case 7:
                        m_Nodes[0] = 0.0;
                        m_Nodes[1] = 4.05845151377397166907e-01;
                        m_Nodes[2] = 7.41531185599394439864e-01;
                        m_Nodes[3] = 9.49107912342758524526e-01;
                        m_Weights[0] = 4.17959183673469387755e-01;
                        m_Weights[1] = 3.81830050505118944950e-01;
                        m_Weights[2] = 2.79705391489276667901e-01;
                        m_Weights[3] = 1.29484966168869693271e-01;
                        break;

                    case 8:
                        m_Nodes[0] = 1.83434642495649804939e-01;
                        m_Nodes[1] = 5.25532409916328985818e-01;
                        m_Nodes[2] = 7.96666477413626739592e-01;
                        m_Nodes[3] = 9.60289856497536231684e-01;
                        m_Weights[0] = 3.62683783378361982965e-01;
                        m_Weights[1] = 3.13706645877887287338e-01;
                        m_Weights[2] = 2.22381034453374470544e-01;
                        m_Weights[3] = 1.01228536290376259153e-01;
                        break;

                    case 9:
                        m_Nodes[0] = 0.0;
                        m_Nodes[1] = 3.24253423403808929039e-01;
                        m_Nodes[2] = 6.13371432700590397309e-01;
                        m_Nodes[3] = 8.36031107326635794299e-01;
                        m_Nodes[4] = 9.68160239507626089836e-01;
                        m_Weights[0] = 3.30239355001259763165e-01;
                        m_Weights[1] = 3.12347077040002840069e-01;
                        m_Weights[2] = 2.60610696402935462319e-01;
                        m_Weights[3] = 1.80648160694857404058e-01;
                        m_Weights[4] = 8.12743883615744119719e-02;
                        break;

                    case 10:
                        m_Nodes[0] = 1.48874338981631210885e-01;
                        m_Nodes[1] = 4.33395394129247190799e-01;
                        m_Nodes[2] = 6.79409568299024406234e-01;
                        m_Nodes[3] = 8.65063366688984510732e-01;
                        m_Nodes[4] = 9.73906528517171720078e-01;
                        m_Weights[0] = 2.95524224714752870174e-01;
                        m_Weights[1] = 2.69266719309996355091e-01;
                        m_Weights[2] = 2.19086362515982043996e-01;
                        m_Weights[3] = 1.49451349150580593146e-01;
                        m_Weights[4] = 6.66713443086881375936e-02;
                        break;

                    case 11:
                        m_Nodes[0] = 0.0;
                        m_Nodes[1] = 2.69543155952344972332e-01;
                        m_Nodes[2] = 5.19096129206811815926e-01;
                        m_Nodes[3] = 7.30152005574049324093e-01;
                        m_Nodes[4] = 8.87062599768095299075e-01;
                        m_Nodes[5] = 9.78228658146056992804e-01;
                        m_Weights[0] = 2.72925086777900630714e-01;
                        m_Weights[1] = 2.62804544510246662181e-01;
                        m_Weights[2] = 2.33193764591990479919e-01;
                        m_Weights[3] = 1.86290210927734251426e-01;
                        m_Weights[4] = 1.25580369464904624635e-01;
                        m_Weights[5] = 5.56685671161736664828e-02;
                        break;

                    case 12:
                        m_Nodes[0] = 1.25233408511468915472e-01;
                        m_Nodes[1] = 3.67831498998180193753e-01;
                        m_Nodes[2] = 5.87317954286617447297e-01;
                        m_Nodes[3] = 7.69902674194304687037e-01;
                        m_Nodes[4] = 9.04117256370474856678e-01;
                        m_Nodes[5] = 9.81560634246719250691e-01;
                        m_Weights[0] = 2.49147045813402785001e-01;
                        m_Weights[1] = 2.33492536538354808761e-01;
                        m_Weights[2] = 2.03167426723065921749e-01;
                        m_Weights[3] = 1.60078328543346226335e-01;
                        m_Weights[4] = 1.06939325995318430960e-01;
                        m_Weights[5] = 4.71753363865118271946e-02;
                        break;

                    case 13:
                        m_Nodes[0] = 0.0;
                        m_Nodes[1] = 2.30458315955134794066e-01;
                        m_Nodes[2] = 4.48492751036446852878e-01;
                        m_Nodes[3] = 6.42349339440340220644e-01;
                        m_Nodes[4] = 8.01578090733309912794e-01;
                        m_Nodes[5] = 9.17598399222977965207e-01;
                        m_Nodes[6] = 9.84183054718588149473e-01;
                        m_Weights[0] = 2.32551553230873910195e-01;
                        m_Weights[1] = 2.26283180262897238412e-01;
                        m_Weights[2] = 2.07816047536888502313e-01;
                        m_Weights[3] = 1.78145980761945738280e-01;
                        m_Weights[4] = 1.38873510219787238464e-01;
                        m_Weights[5] = 9.21214998377284479144e-02;
                        m_Weights[6] = 4.04840047653158795200e-02;
                        break;

                    case 14:
                        m_Nodes[0] = 1.08054948707343662066e-01;
                        m_Nodes[1] = 3.19112368927889760436e-01;
                        m_Nodes[2] = 5.15248636358154091965e-01;
                        m_Nodes[3] = 6.87292904811685470148e-01;
                        m_Nodes[4] = 8.27201315069764993190e-01;
                        m_Nodes[5] = 9.28434883663573517336e-01;
                        m_Nodes[6] = 9.86283808696812338842e-01;
                        m_Weights[0] = 2.15263853463157790196e-01;
                        m_Weights[1] = 2.05198463721295603966e-01;
                        m_Weights[2] = 1.85538397477937813742e-01;
                        m_Weights[3] = 1.57203167158193534570e-01;
                        m_Weights[4] = 1.21518570687903184689e-01;
                        m_Weights[5] = 8.01580871597602098056e-02;
                        m_Weights[6] = 3.51194603317518630318e-02;
                        break;

                    case 15:
                        m_Nodes[0] = 0.0;
                        m_Nodes[1] = 2.01194093997434522301e-01;
                        m_Nodes[2] = 3.94151347077563369897e-01;
                        m_Nodes[3] = 5.70972172608538847537e-01;
                        m_Nodes[4] = 7.24417731360170047416e-01;
                        m_Nodes[5] = 8.48206583410427216201e-01;
                        m_Nodes[6] = 9.37273392400705904308e-01;
                        m_Nodes[7] = 9.87992518020485428490e-01;
                        m_Weights[0] = 2.02578241925561272881e-01;
                        m_Weights[1] = 1.98431485327111576456e-01;
                        m_Weights[2] = 1.86161000015562211027e-01;
                        m_Weights[3] = 1.66269205816993933553e-01;
                        m_Weights[4] = 1.39570677926154314448e-01;
                        m_Weights[5] = 1.07159220467171935012e-01;
                        m_Weights[6] = 7.03660474881081247093e-02;
                        m_Weights[7] = 3.07532419961172683546e-02;
                        break;

                    case 16:
                        m_Nodes[0] = 9.50125098376374401853e-02;
                        m_Nodes[1] = 2.81603550779258913230e-01;
                        m_Nodes[2] = 4.58016777657227386342e-01;
                        m_Nodes[3] = 6.17876244402643748447e-01;
                        m_Nodes[4] = 7.55404408355003033895e-01;
                        m_Nodes[5] = 8.65631202387831743880e-01;
                        m_Nodes[6] = 9.44575023073232576078e-01;
                        m_Nodes[7] = 9.89400934991649932596e-01;
                        m_Weights[0] = 1.89450610455068496285e-01;
                        m_Weights[1] = 1.82603415044923588867e-01;
                        m_Weights[2] = 1.69156519395002538189e-01;
                        m_Weights[3] = 1.49595988816576732082e-01;
                        m_Weights[4] = 1.24628971255533872052e-01;
                        m_Weights[5] = 9.51585116824927848099e-02;
                        m_Weights[6] = 6.22535239386478928628e-02;
                        m_Weights[7] = 2.71524594117540948518e-02;
                        break;

                    case 17:
                        m_Nodes[0] = 0.0;
                        m_Nodes[1] = 1.78484181495847855851e-01;
                        m_Nodes[2] = 3.51231763453876315297e-01;
                        m_Nodes[3] = 5.12690537086476967886e-01;
                        m_Nodes[4] = 6.57671159216690765850e-01;
                        m_Nodes[5] = 7.81514003896801406925e-01;
                        m_Nodes[6] = 8.80239153726985902123e-01;
                        m_Nodes[7] = 9.50675521768767761223e-01;
                        m_Nodes[8] = 9.90575475314417335675e-01;
                        m_Weights[0] = 1.79446470356206525458e-01;
                        m_Weights[1] = 1.76562705366992646325e-01;
                        m_Weights[2] = 1.68004102156450044510e-01;
                        m_Weights[3] = 1.54045761076810288081e-01;
                        m_Weights[4] = 1.35136368468525473286e-01;
                        m_Weights[5] = 1.11883847193403971095e-01;
                        m_Weights[6] = 8.50361483171791808835e-02;
                        m_Weights[7] = 5.54595293739872011294e-02;
                        m_Weights[8] = 2.41483028685479319601e-02;
                        break;

                    case 18:
                        m_Nodes[0] = 8.47750130417353012423e-02;
                        m_Nodes[1] = 2.51886225691505509589e-01;
                        m_Nodes[2] = 4.11751161462842646036e-01;
                        m_Nodes[3] = 5.59770831073947534608e-01;
                        m_Nodes[4] = 6.91687043060353207875e-01;
                        m_Nodes[5] = 8.03704958972523115682e-01;
                        m_Nodes[6] = 8.92602466497555739206e-01;
                        m_Nodes[7] = 9.55823949571397755181e-01;
                        m_Nodes[8] = 9.91565168420930946730e-01;
                        m_Weights[0] = 1.69142382963143591841e-01;
                        m_Weights[1] = 1.64276483745832722986e-01;
                        m_Weights[2] = 1.54684675126265244925e-01;
                        m_Weights[3] = 1.40642914670650651205e-01;
                        m_Weights[4] = 1.22555206711478460185e-01;
                        m_Weights[5] = 1.00942044106287165563e-01;
                        m_Weights[6] = 7.64257302548890565291e-02;
                        m_Weights[7] = 4.97145488949697964533e-02;
                        m_Weights[8] = 2.16160135264833103133e-02;
                        break;

                    case 19:
                        m_Nodes[0] = 0.0;
                        m_Nodes[1] = 1.60358645640225375868e-01;
                        m_Nodes[2] = 3.16564099963629831990e-01;
                        m_Nodes[3] = 4.64570741375960945717e-01;
                        m_Nodes[4] = 6.00545304661681023470e-01;
                        m_Nodes[5] = 7.20966177335229378617e-01;
                        m_Nodes[6] = 8.22714656537142824979e-01;
                        m_Nodes[7] = 9.03155903614817901643e-01;
                        m_Nodes[8] = 9.60208152134830030853e-01;
                        m_Nodes[9] = 9.92406843843584403189e-01;
                        m_Weights[0] = 1.61054449848783695979e-01;
                        m_Weights[1] = 1.58968843393954347650e-01;
                        m_Weights[2] = 1.52766042065859666779e-01;
                        m_Weights[3] = 1.42606702173606611776e-01;
                        m_Weights[4] = 1.28753962539336227676e-01;
                        m_Weights[5] = 1.11566645547333994716e-01;
                        m_Weights[6] = 9.14900216224499994645e-02;
                        m_Weights[7] = 6.90445427376412265807e-02;
                        m_Weights[8] = 4.48142267656996003328e-02;
                        m_Weights[9] = 1.94617882297264770363e-02;
                        break;

                    case 20:
                        m_Nodes[0] = 7.65265211334973337546e-02;
                        m_Nodes[1] = 2.27785851141645078080e-01;
                        m_Nodes[2] = 3.73706088715419560673e-01;
                        m_Nodes[3] = 5.10867001950827098004e-01;
                        m_Nodes[4] = 6.36053680726515025453e-01;
                        m_Nodes[5] = 7.46331906460150792614e-01;
                        m_Nodes[6] = 8.39116971822218823395e-01;
                        m_Nodes[7] = 9.12234428251325905868e-01;
                        m_Nodes[8] = 9.63971927277913791268e-01;
                        m_Nodes[9] = 9.93128599185094924786e-01;
                        m_Weights[0] = 1.52753387130725850698e-01;
                        m_Weights[1] = 1.49172986472603746788e-01;
                        m_Weights[2] = 1.42096109318382051329e-01;
                        m_Weights[3] = 1.31688638449176626898e-01;
                        m_Weights[4] = 1.18194531961518417312e-01;
                        m_Weights[5] = 1.01930119817240435037e-01;
                        m_Weights[6] = 8.32767415767047487248e-02;
                        m_Weights[7] = 6.26720483341090635695e-02;
                        m_Weights[8] = 4.06014298003869413310e-02;
                        m_Weights[9] = 1.76140071391521183119e-02;
                        break;

                    case 21:
                        m_Nodes[0] = 0.0;
                        m_Nodes[1] = 1.45561854160895090937e-01;
                        m_Nodes[2] = 2.88021316802401096601e-01;
                        m_Nodes[3] = 4.24342120207438783574e-01;
                        m_Nodes[4] = 5.51618835887219807059e-01;
                        m_Nodes[5] = 6.67138804197412319306e-01;
                        m_Nodes[6] = 7.68439963475677908616e-01;
                        m_Nodes[7] = 8.53363364583317283647e-01;
                        m_Nodes[8] = 9.20099334150400828790e-01;
                        m_Nodes[9] = 9.67226838566306294317e-01;
                        m_Nodes[10] = 9.93752170620389500260e-01;
                        m_Weights[0] = 1.46081133649690427192e-01;
                        m_Weights[1] = 1.44524403989970059064e-01;
                        m_Weights[2] = 1.39887394791073154722e-01;
                        m_Weights[3] = 1.32268938633337461781e-01;
                        m_Weights[4] = 1.21831416053728534195e-01;
                        m_Weights[5] = 1.08797299167148377663e-01;
                        m_Weights[6] = 9.34444234560338615533e-02;
                        m_Weights[7] = 7.61001136283793020171e-02;
                        m_Weights[8] = 5.71344254268572082836e-02;
                        m_Weights[9] = 3.69537897708524938000e-02;
                        m_Weights[10] = 1.60172282577743333242e-02;
                        break;

                    case 22:
                        m_Nodes[0] = 6.97392733197222212138e-02;
                        m_Nodes[1] = 2.07860426688221285479e-01;
                        m_Nodes[2] = 3.41935820892084225158e-01;
                        m_Nodes[3] = 4.69355837986757026406e-01;
                        m_Nodes[4] = 5.87640403506911592959e-01;
                        m_Nodes[5] = 6.94487263186682780051e-01;
                        m_Nodes[6] = 7.87816805979208162004e-01;
                        m_Nodes[7] = 8.65812577720300136536e-01;
                        m_Nodes[8] = 9.26956772187174000521e-01;
                        m_Nodes[9] = 9.70060497835428727124e-01;
                        m_Nodes[10] = 9.94294585482399292073e-01;
                        m_Weights[0] = 1.39251872855631993375e-01;
                        m_Weights[1] = 1.36541498346015171353e-01;
                        m_Weights[2] = 1.31173504787062370733e-01;
                        m_Weights[3] = 1.23252376810512424286e-01;
                        m_Weights[4] = 1.12932296080539218393e-01;
                        m_Weights[5] = 1.00414144442880964932e-01;
                        m_Weights[6] = 8.59416062170677274144e-02;
                        m_Weights[7] = 6.97964684245204880950e-02;
                        m_Weights[8] = 5.22933351526832859403e-02;
                        m_Weights[9] = 3.37749015848141547933e-02;
                        m_Weights[10] = 1.46279952982722006850e-02;
                        break;

                    case 23:
                        m_Nodes[0] = 0.00000000000000000000e+00;
                        m_Nodes[1] = 1.33256824298466110932e-01;
                        m_Nodes[2] = 2.64135680970344930534e-01;
                        m_Nodes[3] = 3.90301038030290831421e-01;
                        m_Nodes[4] = 5.09501477846007549690e-01;
                        m_Nodes[5] = 6.19609875763646156385e-01;
                        m_Nodes[6] = 7.18661363131950194462e-01;
                        m_Nodes[7] = 8.04888401618839892151e-01;
                        m_Nodes[8] = 8.76752358270441667378e-01;
                        m_Nodes[9] = 9.32971086826016102349e-01;
                        m_Nodes[10] = 9.72542471218115231956e-01;
                        m_Nodes[11] = 9.94769334997552123524e-01;
                        m_Weights[0] = 1.33654572186106175351e-01;
                        m_Weights[1] = 1.32462039404696617372e-01;
                        m_Weights[2] = 1.28905722188082149979e-01;
                        m_Weights[3] = 1.23049084306729530468e-01;
                        m_Weights[4] = 1.14996640222411364942e-01;
                        m_Weights[5] = 1.04892091464541410074e-01;
                        m_Weights[6] = 9.29157660600351474770e-02;
                        m_Weights[7] = 7.92814117767189549229e-02;
                        m_Weights[8] = 6.42324214085258521272e-02;
                        m_Weights[9] = 4.80376717310846685716e-02;
                        m_Weights[10] = 3.09880058569794443107e-02;
                        m_Weights[11] = 1.34118594871417720813e-02;
                        break;

                    case 24:
                        m_Nodes[0] = 6.40568928626056260850e-02;
                        m_Nodes[1] = 1.91118867473616309159e-01;
                        m_Nodes[2] = 3.15042679696163374387e-01;
                        m_Nodes[3] = 4.33793507626045138487e-01;
                        m_Nodes[4] = 5.45421471388839535658e-01;
                        m_Nodes[5] = 6.48093651936975569252e-01;
                        m_Nodes[6] = 7.40124191578554364244e-01;
                        m_Nodes[7] = 8.20001985973902921954e-01;
                        m_Nodes[8] = 8.86415527004401034213e-01;
                        m_Nodes[9] = 9.38274552002732758524e-01;
                        m_Nodes[10] = 9.74728555971309498198e-01;
                        m_Nodes[11] = 9.95187219997021360180e-01;
                        m_Weights[0] = 1.27938195346752156974e-01;
                        m_Weights[1] = 1.25837456346828296121e-01;
                        m_Weights[2] = 1.21670472927803391204e-01;
                        m_Weights[3] = 1.15505668053725601353e-01;
                        m_Weights[4] = 1.07444270115965634783e-01;
                        m_Weights[5] = 9.76186521041138882699e-02;
                        m_Weights[6] = 8.61901615319532759172e-02;
                        m_Weights[7] = 7.33464814110803057340e-02;
                        m_Weights[8] = 5.92985849154367807464e-02;
                        m_Weights[9] = 4.42774388174198061686e-02;
                        m_Weights[10] = 2.85313886289336631813e-02;
                        m_Weights[11] = 1.23412297999871995468e-02;
                        break;

                    case 25:
                        m_Nodes[0] = 0.00000000000000000000e+00;
                        m_Nodes[1] = 1.22864692610710396387e-01;
                        m_Nodes[2] = 2.43866883720988432045e-01;
                        m_Nodes[3] = 3.61172305809387837736e-01;
                        m_Nodes[4] = 4.73002731445714960522e-01;
                        m_Nodes[5] = 5.77662930241222967724e-01;
                        m_Nodes[6] = 6.73566368473468364485e-01;
                        m_Nodes[7] = 7.59259263037357630577e-01;
                        m_Nodes[8] = 8.33442628760834001421e-01;
                        m_Nodes[9] = 8.94991997878275368851e-01;
                        m_Nodes[10] = 9.42974571228974339414e-01;
                        m_Nodes[11] = 9.76663921459517511498e-01;
                        m_Nodes[12] = 9.95556969790498097909e-01;
                        m_Weights[0] = 1.23176053726715451204e-01;
                        m_Weights[1] = 1.22242442990310041689e-01;
                        m_Weights[2] = 1.19455763535784772228e-01;
                        m_Weights[3] = 1.14858259145711648339e-01;
                        m_Weights[4] = 1.08519624474263653116e-01;
                        m_Weights[5] = 1.00535949067050644202e-01;
                        m_Weights[6] = 9.10282619829636498115e-02;
                        m_Weights[7] = 8.01407003350010180132e-02;
                        m_Weights[8] = 6.80383338123569172072e-02;
                        m_Weights[9] = 5.49046959758351919259e-02;
                        m_Weights[10] = 4.09391567013063126556e-02;
                        m_Weights[11] = 2.63549866150321372619e-02;
                        m_Weights[12] = 1.13937985010262879479e-02;
                        break;

                    case 26:
                        m_Nodes[0] = 5.92300934293132070937e-02;
                        m_Nodes[1] = 1.76858820356890183969e-01;
                        m_Nodes[2] = 2.92004839485956895143e-01;
                        m_Nodes[3] = 4.03051755123486306481e-01;
                        m_Nodes[4] = 5.08440714824505717696e-01;
                        m_Nodes[5] = 6.06692293017618063232e-01;
                        m_Nodes[6] = 6.96427260419957264864e-01;
                        m_Nodes[7] = 7.76385948820678856193e-01;
                        m_Nodes[8] = 8.45445942788498018798e-01;
                        m_Nodes[9] = 9.02637861984307074218e-01;
                        m_Nodes[10] = 9.47159066661714250136e-01;
                        m_Nodes[11] = 9.78385445956470991101e-01;
                        m_Nodes[12] = 9.95885701145616929003e-01;
                        m_Weights[0] = 1.18321415279262276516e-01;
                        m_Weights[1] = 1.16660443485296582045e-01;
                        m_Weights[2] = 1.13361816546319666549e-01;
                        m_Weights[3] = 1.08471840528576590657e-01;
                        m_Weights[4] = 1.02059161094425423238e-01;
                        m_Weights[5] = 9.42138003559141484637e-02;
                        m_Weights[6] = 8.50458943134852392104e-02;
                        m_Weights[7] = 7.46841497656597458871e-02;
                        m_Weights[8] = 6.32740463295748355395e-02;
                        m_Weights[9] = 5.09758252971478119983e-02;
                        m_Weights[10] = 3.79623832943627639503e-02;
                        m_Weights[11] = 2.44178510926319087896e-02;
                        m_Weights[12] = 1.05513726173430071557e-02;
                        break;

                    case 27:
                        m_Nodes[0] = 0.00000000000000000000e+00;
                        m_Nodes[1] = 1.13972585609529966933e-01;
                        m_Nodes[2] = 2.26459365439536858857e-01;
                        m_Nodes[3] = 3.35993903638508899730e-01;
                        m_Nodes[4] = 4.41148251750026880586e-01;
                        m_Nodes[5] = 5.40551564579456894900e-01;
                        m_Nodes[6] = 6.32907971946495140928e-01;
                        m_Nodes[7] = 7.17013473739423699295e-01;
                        m_Nodes[8] = 7.91771639070508227144e-01;
                        m_Nodes[9] = 8.56207908018294490303e-01;
                        m_Nodes[10] = 9.09482320677491104301e-01;
                        m_Nodes[11] = 9.50900557814705006852e-01;
                        m_Nodes[12] = 9.79923475961501222856e-01;
                        m_Nodes[13] = 9.96179262888988566939e-01;
                        m_Weights[0] = 1.14220867378956989045e-01;
                        m_Weights[1] = 1.13476346108965148620e-01;
                        m_Weights[2] = 1.11252488356845192672e-01;
                        m_Weights[3] = 1.07578285788533187212e-01;
                        m_Weights[4] = 1.02501637817745798671e-01;
                        m_Weights[5] = 9.60887273700285075657e-02;
                        m_Weights[6] = 8.84231585437569501943e-02;
                        m_Weights[7] = 7.96048677730577712631e-02;
                        m_Weights[8] = 6.97488237662455929843e-02;
                        m_Weights[9] = 5.89835368598335991103e-02;
                        m_Weights[10] = 4.74494125206150627040e-02;
                        m_Weights[11] = 3.52970537574197110226e-02;
                        m_Weights[12] = 2.26862315961806231960e-02;
                        m_Weights[13] = 9.79899605129436026115e-03;
                        break;

                    case 28:
                        m_Nodes[0] = 5.50792898840342704265e-02;
                        m_Nodes[1] = 1.64569282133380771281e-01;
                        m_Nodes[2] = 2.72061627635178077677e-01;
                        m_Nodes[3] = 3.76251516089078710221e-01;
                        m_Nodes[4] = 4.75874224955118261034e-01;
                        m_Nodes[5] = 5.69720471811401719308e-01;
                        m_Nodes[6] = 6.56651094038864961220e-01;
                        m_Nodes[7] = 7.35610878013631772028e-01;
                        m_Nodes[8] = 8.05641370917179171448e-01;
                        m_Nodes[9] = 8.65892522574395048942e-01;
                        m_Nodes[10] = 9.15633026392132073870e-01;
                        m_Nodes[11] = 9.54259280628938197254e-01;
                        m_Nodes[12] = 9.81303165370872753695e-01;
                        m_Nodes[13] = 9.96442497573954449950e-01;
                        m_Weights[0] = 1.10047013016475196282e-01;
                        m_Weights[1] = 1.08711192258294135254e-01;
                        m_Weights[2] = 1.06055765922846417910e-01;
                        m_Weights[3] = 1.02112967578060769814e-01;
                        m_Weights[4] = 9.69306579979299158505e-02;
                        m_Weights[5] = 9.05717443930328409422e-02;
                        m_Weights[6] = 8.31134172289012183904e-02;
                        m_Weights[7] = 7.46462142345687790239e-02;
                        m_Weights[8] = 6.52729239669995957934e-02;
                        m_Weights[9] = 5.51073456757167454315e-02;
                        m_Weights[10] = 4.42729347590042278396e-02;
                        m_Weights[11] = 3.29014277823043799776e-02;
                        m_Weights[12] = 2.11321125927712597515e-02;
                        m_Weights[13] = 9.12428259309451773881e-03;
                        break;

                    case 29:
                        m_Nodes[0] = 0.00000000000000000000e+00;
                        m_Nodes[1] = 1.06278230132679230171e-01;
                        m_Nodes[2] = 2.11352286166001074506e-01;
                        m_Nodes[3] = 3.14031637867639934948e-01;
                        m_Nodes[4] = 4.13152888174008663891e-01;
                        m_Nodes[5] = 5.07592955124227642103e-01;
                        m_Nodes[6] = 5.96281797138227820380e-01;
                        m_Nodes[7] = 6.78214537602686515156e-01;
                        m_Nodes[8] = 7.52462851734477133913e-01;
                        m_Nodes[9] = 8.18185487615252444990e-01;
                        m_Nodes[10] = 8.74637804920102790418e-01;
                        m_Nodes[11] = 9.21180232953058785094e-01;
                        m_Nodes[12] = 9.57285595778087725798e-01;
                        m_Nodes[13] = 9.82545505261413174871e-01;
                        m_Nodes[14] = 9.96679442260596586163e-01;
                        m_Weights[0] = 1.06479381718314244247e-01;
                        m_Weights[1] = 1.05876155097320941407e-01;
                        m_Weights[2] = 1.04073310077729373913e-01;
                        m_Weights[3] = 1.01091273759914966122e-01;
                        m_Weights[4] = 9.69638340944086063019e-02;
                        m_Weights[5] = 9.17377571392587633480e-02;
                        m_Weights[6] = 8.54722573661725275453e-02;
                        m_Weights[7] = 7.82383271357637838281e-02;
                        m_Weights[8] = 7.01179332550512785696e-02;
                        m_Weights[9] = 6.12030906570791385422e-02;
                        m_Weights[10] = 5.15948269024979239128e-02;
                        m_Weights[11] = 4.14020625186828361047e-02;
                        m_Weights[12] = 3.07404922020936226447e-02;
                        m_Weights[13] = 1.97320850561227059838e-02;
                        m_Weights[14] = 8.51690387874640965428e-03;
                        break;

                    case 30:
                        m_Nodes[0] = 5.14718425553176958330e-02;
                        m_Nodes[1] = 1.53869913608583546964e-01;
                        m_Nodes[2] = 2.54636926167889846440e-01;
                        m_Nodes[3] = 3.52704725530878113471e-01;
                        m_Nodes[4] = 4.47033769538089176781e-01;
                        m_Nodes[5] = 5.36624148142019899264e-01;
                        m_Nodes[6] = 6.20526182989242861140e-01;
                        m_Nodes[7] = 6.97850494793315796932e-01;
                        m_Nodes[8] = 7.67777432104826194918e-01;
                        m_Nodes[9] = 8.29565762382768397443e-01;
                        m_Nodes[10] = 8.82560535792052681543e-01;
                        m_Nodes[11] = 9.26200047429274325879e-01;
                        m_Nodes[12] = 9.60021864968307512217e-01;
                        m_Nodes[13] = 9.83668123279747209970e-01;
                        m_Nodes[14] = 9.96893484074649540272e-01;
                        m_Weights[0] = 1.02852652893558840341e-01;
                        m_Weights[1] = 1.01762389748405504596e-01;
                        m_Weights[2] = 9.95934205867952670628e-02;
                        m_Weights[3] = 9.63687371746442596395e-02;
                        m_Weights[4] = 9.21225222377861287176e-02;
                        m_Weights[5] = 8.68997872010829798024e-02;
                        m_Weights[6] = 8.07558952294202153547e-02;
                        m_Weights[7] = 7.37559747377052062683e-02;
                        m_Weights[8] = 6.59742298821804951281e-02;
                        m_Weights[9] = 5.74931562176190664817e-02;
                        m_Weights[10] = 4.84026728305940529029e-02;
                        m_Weights[11] = 3.87991925696270495970e-02;
                        m_Weights[12] = 2.87847078833233693499e-02;
                        m_Weights[13] = 1.84664683110909591420e-02;
                        m_Weights[14] = 7.96819249616660561542e-03;
                        break;

                    case 31:
                        m_Nodes[0] = 0.00000000000000000000e+00;
                        m_Nodes[1] = 9.95553121523415203252e-02;
                        m_Nodes[2] = 1.98121199335570628772e-01;
                        m_Nodes[3] = 2.94718069981701616618e-01;
                        m_Nodes[4] = 3.88385901608232943061e-01;
                        m_Nodes[5] = 4.78193782044902480441e-01;
                        m_Nodes[6] = 5.63249161407149262721e-01;
                        m_Nodes[7] = 6.42706722924260346184e-01;
                        m_Nodes[8] = 7.15776784586853283906e-01;
                        m_Nodes[9] = 7.81733148416624940406e-01;
                        m_Nodes[10] = 8.39920320146267340087e-01;
                        m_Nodes[11] = 8.89760029948271043374e-01;
                        m_Nodes[12] = 9.30756997896648164957e-01;
                        m_Nodes[13] = 9.62503925092949661789e-01;
                        m_Nodes[14] = 9.84685909665152484002e-01;
                        m_Nodes[15] = 9.97087481819477074056e-01;
                        m_Weights[0] = 9.97205447934264514275e-02;
                        m_Weights[1] = 9.92250112266723078749e-02;
                        m_Weights[2] = 9.77433353863287250935e-02;
                        m_Weights[3] = 9.52902429123195128072e-02;
                        m_Weights[4] = 9.18901138936414782154e-02;
                        m_Weights[5] = 8.75767406084778761262e-02;
                        m_Weights[6] = 8.23929917615892639038e-02;
                        m_Weights[7] = 7.63903865987766164264e-02;
                        m_Weights[8] = 6.96285832354103661677e-02;
                        m_Weights[9] = 6.21747865610284269104e-02;
                        m_Weights[10] = 5.41030824249168537116e-02;
                        m_Weights[11] = 4.54937075272011029027e-02;
                        m_Weights[12] = 3.64322739123854640251e-02;
                        m_Weights[13] = 2.70090191849794218006e-02;
                        m_Weights[14] = 1.73186207903105824628e-02;
                        m_Weights[15] = 7.47083157924877585870e-03;
                        break;

                    case 32:
                        m_Nodes[0] = 4.83076656877383162348e-02;
                        m_Nodes[1] = 1.44471961582796493485e-01;
                        m_Nodes[2] = 2.39287362252137074545e-01;
                        m_Nodes[3] = 3.31868602282127649780e-01;
                        m_Nodes[4] = 4.21351276130635345364e-01;
                        m_Nodes[5] = 5.06899908932229390024e-01;
                        m_Nodes[6] = 5.87715757240762329041e-01;
                        m_Nodes[7] = 6.63044266930215200975e-01;
                        m_Nodes[8] = 7.32182118740289680387e-01;
                        m_Nodes[9] = 7.94483795967942406963e-01;
                        m_Nodes[10] = 8.49367613732569970134e-01;
                        m_Nodes[11] = 8.96321155766052123965e-01;
                        m_Nodes[12] = 9.34906075937739689171e-01;
                        m_Nodes[13] = 9.64762255587506430774e-01;
                        m_Nodes[14] = 9.85611511545268335400e-01;
                        m_Nodes[15] = 9.97263861849481563545e-01;
                        m_Weights[0] = 9.65400885147278005668e-02;
                        m_Weights[1] = 9.56387200792748594191e-02;
                        m_Weights[2] = 9.38443990808045656392e-02;
                        m_Weights[3] = 9.11738786957638847129e-02;
                        m_Weights[4] = 8.76520930044038111428e-02;
                        m_Weights[5] = 8.33119242269467552222e-02;
                        m_Weights[6] = 7.81938957870703064717e-02;
                        m_Weights[7] = 7.23457941088485062254e-02;
                        m_Weights[8] = 6.58222227763618468376e-02;
                        m_Weights[9] = 5.86840934785355471452e-02;
                        m_Weights[10] = 5.09980592623761761971e-02;
                        m_Weights[11] = 4.28358980222266806542e-02;
                        m_Weights[12] = 3.42738629130214330967e-02;
                        m_Weights[13] = 2.53920653092620594547e-02;
                        m_Weights[14] = 1.62743947309056706054e-02;
                        m_Weights[15] = 7.01861000947009660046e-03;
                        break;

                    case 33:
                        m_Nodes[0] = 0.00000000000000000000e+00;
                        m_Nodes[1] = 9.36310658547333856707e-02;
                        m_Nodes[2] = 1.86439298827991572336e-01;
                        m_Nodes[3] = 2.77609097152497029403e-01;
                        m_Nodes[4] = 3.66339257748073341070e-01;
                        m_Nodes[5] = 4.51850017272450695726e-01;
                        m_Nodes[6] = 5.33389904786347643549e-01;
                        m_Nodes[7] = 6.10242345836379027307e-01;
                        m_Nodes[8] = 6.81731959969742786268e-01;
                        m_Nodes[9] = 7.47230496449562157859e-01;
                        m_Nodes[10] = 8.06162356274166589796e-01;
                        m_Nodes[11] = 8.58009652676504064643e-01;
                        m_Nodes[12] = 9.02316767743433583041e-01;
                        m_Nodes[13] = 9.38694372611168350356e-01;
                        m_Nodes[14] = 9.66822909689992768928e-01;
                        m_Nodes[15] = 9.86455726230642488110e-01;
                        m_Nodes[16] = 9.97424694246455217266e-01;
                        m_Weights[0] = 9.37684461602099965673e-02;
                        m_Weights[1] = 9.33564260655961161610e-02;
                        m_Weights[2] = 9.21239866433168462132e-02;
                        m_Weights[3] = 9.00819586606385772397e-02;
                        m_Weights[4] = 8.72482876188443376073e-02;
                        m_Weights[5] = 8.36478760670387076139e-02;
                        m_Weights[6] = 7.93123647948867383639e-02;
                        m_Weights[7] = 7.42798548439541493425e-02;
                        m_Weights[8] = 6.85945728186567128060e-02;
                        m_Weights[9] = 6.23064825303174800316e-02;
                        m_Weights[10] = 5.54708466316635612852e-02;
                        m_Weights[11] = 4.81477428187116956680e-02;
                        m_Weights[12] = 4.04015413316695915654e-02;
                        m_Weights[13] = 3.23003586323289532835e-02;
                        m_Weights[14] = 2.39155481017494803495e-02;
                        m_Weights[15] = 1.53217015129346761276e-02;
                        m_Weights[16] = 6.60622784758737805902e-03;
                        break;

                    case 34:
                        m_Nodes[0] = 4.55098219531025427491e-02;
                        m_Nodes[1] = 1.36152357259182975894e-01;
                        m_Nodes[2] = 2.25666691616449483869e-01;
                        m_Nodes[3] = 3.13311081339463247458e-01;
                        m_Nodes[4] = 3.98359277758645940631e-01;
                        m_Nodes[5] = 4.80106545190327034194e-01;
                        m_Nodes[6] = 5.57875500669746642736e-01;
                        m_Nodes[7] = 6.31021727080528545318e-01;
                        m_Nodes[8] = 6.98939113216262907933e-01;
                        m_Nodes[9] = 7.61064876629873014187e-01;
                        m_Nodes[10] = 8.16884227900933664592e-01;
                        m_Nodes[11] = 8.65934638334564469264e-01;
                        m_Nodes[12] = 9.07809677718324468801e-01;
                        m_Nodes[13] = 9.42162397405107091632e-01;
                        m_Nodes[14] = 9.68708262533344281765e-01;
                        m_Nodes[15] = 9.87227816406309485050e-01;
                        m_Nodes[16] = 9.97571753790841919243e-01;
                        m_Weights[0] = 9.09567403302598736153e-02;
                        m_Weights[1] = 9.02030443706407295739e-02;
                        m_Weights[2] = 8.87018978356938692871e-02;
                        m_Weights[3] = 8.64657397470357497842e-02;
                        m_Weights[4] = 8.35130996998456551870e-02;
                        m_Weights[5] = 7.98684443397718447388e-02;
                        m_Weights[6] = 7.55619746600319312708e-02;
                        m_Weights[7] = 7.06293758142557249991e-02;
                        m_Weights[8] = 6.51115215540764113783e-02;
                        m_Weights[9] = 5.90541358275244931938e-02;
                        m_Weights[10] = 5.25074145726781061670e-02;
                        m_Weights[11] = 4.55256115233532724534e-02;
                        m_Weights[12] = 3.81665937963875163188e-02;
                        m_Weights[13] = 3.04913806384461318145e-02;
                        m_Weights[14] = 2.25637219854949700802e-02;
                        m_Weights[15] = 1.44501627485950354180e-02;
                        m_Weights[16] = 6.22914055590868471921e-03;
                        break;

                    case 35:
                        m_Nodes[0] = 0.00000000000000000000e+00;
                        m_Nodes[1] = 8.83713432756592636009e-02;
                        m_Nodes[2] = 1.76051061165989569974e-01;
                        m_Nodes[3] = 2.62352941209296057971e-01;
                        m_Nodes[4] = 3.46601554430813945877e-01;
                        m_Nodes[5] = 4.28137541517814254188e-01;
                        m_Nodes[6] = 5.06322773241488615024e-01;
                        m_Nodes[7] = 5.80545344749764509935e-01;
                        m_Nodes[8] = 6.50224364665890388676e-01;
                        m_Nodes[9] = 7.14814501556628783264e-01;
                        m_Nodes[10] = 7.73810252286912555267e-01;
                        m_Nodes[11] = 8.26749899092225406834e-01;
                        m_Nodes[12] = 8.73219125025222331523e-01;
                        m_Nodes[13] = 9.12854261359317614465e-01;
                        m_Nodes[14] = 9.45345148207827329539e-01;
                        m_Nodes[15] = 9.70437616039229833215e-01;
                        m_Nodes[16] = 9.87935764443851498035e-01;
                        m_Nodes[17] = 9.97706569099600297260e-01;
                        m_Weights[0] = 8.84867949071042906382e-02;
                        m_Weights[1] = 8.81405304302754629707e-02;
                        m_Weights[2] = 8.71044469971835342433e-02;
                        m_Weights[3] = 8.53866533920991252259e-02;
                        m_Weights[4] = 8.30005937288565883799e-02;
                        m_Weights[5] = 7.99649422423242629327e-02;
                        m_Weights[6] = 7.63034571554420535387e-02;
                        m_Weights[7] = 7.20447947725600646655e-02;
                        m_Weights[8] = 6.72222852690869039644e-02;
                        m_Weights[9] = 6.18736719660801888873e-02;
                        m_Weights[10] = 5.60408162123701285785e-02;
                        m_Weights[11] = 4.97693704013535298006e-02;
                        m_Weights[12] = 4.31084223261702187938e-02;
                        m_Weights[13] = 3.61101158634633805399e-02;
                        m_Weights[14] = 2.88292601088942540886e-02;
                        m_Weights[15] = 2.13229799114835808525e-02;
                        m_Weights[16] = 1.36508283483614922706e-02;
                        m_Weights[17] = 5.88343342044308497574e-03;
                        break;

                    case 36:
                        m_Nodes[0] = 4.30181984737086072270e-02;
                        m_Nodes[1] = 1.28736103809384788652e-01;
                        m_Nodes[2] = 2.13500892316865578943e-01;
                        m_Nodes[3] = 2.96684995344028270503e-01;
                        m_Nodes[4] = 3.77672547119689216323e-01;
                        m_Nodes[5] = 4.55863944433420267207e-01;
                        m_Nodes[6] = 5.30680285926245161641e-01;
                        m_Nodes[7] = 6.01567658135980535080e-01;
                        m_Nodes[8] = 6.68001236585521062097e-01;
                        m_Nodes[9] = 7.29489171593556582090e-01;
                        m_Nodes[10] = 7.85576230132206512828e-01;
                        m_Nodes[11] = 8.35847166992475306419e-01;
                        m_Nodes[12] = 8.79929800890397131982e-01;
                        m_Nodes[13] = 9.17497774515659066076e-01;
                        m_Nodes[14] = 9.48272984399507545202e-01;
                        m_Nodes[15] = 9.72027691049697949336e-01;
                        m_Nodes[16] = 9.88586478902212238073e-01;
                        m_Nodes[17] = 9.97830462484085836199e-01;
                        m_Weights[0] = 8.59832756703947474901e-02;
                        m_Weights[1] = 8.53466857393386274919e-02;
                        m_Weights[2] = 8.40782189796619349335e-02;
                        m_Weights[3] = 8.21872667043397095172e-02;
                        m_Weights[4] = 7.96878289120716019087e-02;
                        m_Weights[5] = 7.65984106458706745288e-02;
                        m_Weights[6] = 7.29418850056530613539e-02;
                        m_Weights[7] = 6.87453238357364426136e-02;
                        m_Weights[8] = 6.40397973550154895568e-02;
                        m_Weights[9] = 5.88601442453248173083e-02;
                        m_Weights[10] = 5.32447139777599190929e-02;
                        m_Weights[11] = 4.72350834902659784121e-02;
                        m_Weights[12] = 4.08757509236448954638e-02;
                        m_Weights[13] = 3.42138107703072299983e-02;
                        m_Weights[14] = 2.72986214985687790906e-02;
                        m_Weights[15] = 2.01815152977354714793e-02;
                        m_Weights[16] = 1.29159472840655744133e-02;
                        m_Weights[17] = 5.56571966424504535855e-03;
                        break;

                    case 37:
                        m_Nodes[0] = 0.00000000000000000000e+00;
                        m_Nodes[1] = 8.36704089547699019430e-02;
                        m_Nodes[2] = 1.66753930239851976969e-01;
                        m_Nodes[3] = 2.48667792791365758806e-01;
                        m_Nodes[4] = 3.28837429883706999498e-01;
                        m_Nodes[5] = 4.06700509318326110101e-01;
                        m_Nodes[6] = 4.81710877803205554147e-01;
                        m_Nodes[7] = 5.53342391861581781235e-01;
                        m_Nodes[8] = 6.21092608408924483148e-01;
                        m_Nodes[9] = 6.84486309130959357446e-01;
                        m_Nodes[10] = 7.43078833981965262547e-01;
                        m_Nodes[11] = 7.96459200509902293393e-01;
                        m_Nodes[12] = 8.44252987340555967987e-01;
                        m_Nodes[13] = 8.86124962155486078946e-01;
                        m_Nodes[14] = 9.21781437412463742668e-01;
                        m_Nodes[15] = 9.50972343262094821329e-01;
                        m_Nodes[16] = 9.73493030056485744329e-01;
                        m_Nodes[17] = 9.89185963214319186684e-01;
                        m_Nodes[18] = 9.97944582477913648941e-01;
                        m_Weights[0] = 8.37683609931389047970e-02;
                        m_Weights[1] = 8.34745736258627872523e-02;
                        m_Weights[2] = 8.25952722364372508912e-02;
                        m_Weights[3] = 8.11366245084650305099e-02;
                        m_Weights[4] = 7.91088618375293807672e-02;
                        m_Weights[5] = 7.65262075705292378859e-02;
                        m_Weights[6] = 7.34067772484881727246e-02;
                        m_Weights[7] = 6.97724515557003448851e-02;
                        m_Weights[8] = 6.56487228727512494840e-02;
                        m_Weights[9] = 6.10645165232259861301e-02;
                        m_Weights[10] = 5.60519879982749178069e-02;
                        m_Weights[11] = 5.06462976548246015938e-02;
                        m_Weights[12] = 4.48853646624371665978e-02;
                        m_Weights[13] = 3.88096025019345445047e-02;
                        m_Weights[14] = 3.24616398475214811211e-02;
                        m_Weights[15] = 2.58860369905589337262e-02;
                        m_Weights[16] = 1.91290444890839659320e-02;
                        m_Weights[17] = 1.22387801003075565798e-02;
                        m_Weights[18] = 5.27305727949793934979e-03;
                        break;

                    case 38:
                        m_Nodes[0] = 4.07851479045782399133e-02;
                        m_Nodes[1] = 1.22084025337867419870e-01;
                        m_Nodes[2] = 2.02570453892116703204e-01;
                        m_Nodes[3] = 2.81708809790165261360e-01;
                        m_Nodes[4] = 3.58972440479435013257e-01;
                        m_Nodes[5] = 4.33847169432376484373e-01;
                        m_Nodes[6] = 5.05834717927931103241e-01;
                        m_Nodes[7] = 5.74456021047807081133e-01;
                        m_Nodes[8] = 6.39254415829681707180e-01;
                        m_Nodes[9] = 6.99798680379184355913e-01;
                        m_Nodes[10] = 7.55685903753970680738e-01;
                        m_Nodes[11] = 8.06544167605316815552e-01;
                        m_Nodes[12] = 8.52035021932362188860e-01;
                        m_Nodes[13] = 8.91855739004632216795e-01;
                        m_Nodes[14] = 9.25741332048584396825e-01;
                        m_Nodes[15] = 9.53466330933529595671e-01;
                        m_Nodes[16] = 9.74846328590153507641e-01;
                        m_Nodes[17] = 9.89739454266385571944e-01;
                        m_Nodes[18] = 9.98049930535687619813e-01;
                        m_Weights[0] = 8.15250292803857866992e-02;
                        m_Weights[1] = 8.09824937705971006233e-02;
                        m_Weights[2] = 7.99010332435278215860e-02;
                        m_Weights[3] = 7.82878446582109480754e-02;
                        m_Weights[4] = 7.61536635484463960660e-02;
                        m_Weights[5] = 7.35126925847434571452e-02;
                        m_Weights[6] = 7.03825070668989547393e-02;
                        m_Weights[7] = 6.67839379791404119350e-02;
                        m_Weights[8] = 6.27409333921330540525e-02;
                        m_Weights[9] = 5.82803991469972060225e-02;
                        m_Weights[10] = 5.34320199103323199724e-02;
                        m_Weights[11] = 4.82280618607586833957e-02;
                        m_Weights[12] = 4.27031585046744343133e-02;
                        m_Weights[13] = 3.68940815940247381371e-02;
                        m_Weights[14] = 3.08395005451750546765e-02;
                        m_Weights[15] = 2.45797397382323759012e-02;
                        m_Weights[16] = 1.81565777096132369760e-02;
                        m_Weights[17] = 1.16134447164686742626e-02;
                        m_Weights[18] = 5.00288074963934566306e-03;
                        break;

                    case 39:
                        m_Nodes[0] = 0.00000000000000000000e+00;
                        m_Nodes[1] = 7.94438046087554775819e-02;
                        m_Nodes[2] = 1.58385339997837799923e-01;
                        m_Nodes[3] = 2.36325512461835767336e-01;
                        m_Nodes[4] = 3.12771559248185922536e-01;
                        m_Nodes[5] = 3.87240163971561455854e-01;
                        m_Nodes[6] = 4.59260512309136048663e-01;
                        m_Nodes[7] = 5.28377268660437473896e-01;
                        m_Nodes[8] = 5.94153454957277988693e-01;
                        m_Nodes[9] = 6.56173213432010910734e-01;
                        m_Nodes[10] = 7.14044435894534679134e-01;
                        m_Nodes[11] = 7.67401242931063499832e-01;
                        m_Nodes[12] = 8.15906297430143104353e-01;
                        m_Nodes[13] = 8.59252937999906153914e-01;
                        m_Nodes[14] = 8.97167119292992887848e-01;
                        m_Nodes[15] = 9.29409148486738229698e-01;
                        m_Nodes[16] = 9.55775212324652277111e-01;
                        m_Nodes[17] = 9.76098709333471053845e-01;
                        m_Nodes[18] = 9.90251536854685983640e-01;
                        m_Nodes[19] = 9.98147383066432906005e-01;
                        m_Weights[0] = 7.95276221394428524174e-02;
                        m_Weights[1] = 7.92762225683684710102e-02;
                        m_Weights[2] = 7.85236132873711767251e-02;
                        m_Weights[3] = 7.72745525446820167285e-02;
                        m_Weights[4] = 7.55369373228360577048e-02;
                        m_Weights[5] = 7.33217534142686173812e-02;
                        m_Weights[6] = 7.06430059706087607701e-02;
                        m_Weights[7] = 6.75176309662312653630e-02;
                        m_Weights[8] = 6.39653881386823889867e-02;
                        m_Weights[9] = 6.00087360885961495760e-02;
                        m_Weights[10] = 5.56726903409162999049e-02;
                        m_Weights[11] = 5.09846652921294051904e-02;
                        m_Weights[12] = 4.59743011089166319041e-02;
                        m_Weights[13] = 4.06732768479338442186e-02;
                        m_Weights[14] = 3.51151114981313311066e-02;
                        m_Weights[15] = 2.93349559839033771339e-02;
                        m_Weights[16] = 2.33693848321781631398e-02;
                        m_Weights[17] = 1.72562290937249194889e-02;
                        m_Weights[18] = 1.10347889391645942057e-02;
                        m_Weights[19] = 4.75294469163510133601e-03;
                        break;

                    case 40:
                        m_Nodes[0] = 3.87724175060508219332e-02;
                        m_Nodes[1] = 1.16084070675255208483e-01;
                        m_Nodes[2] = 1.92697580701371099716e-01;
                        m_Nodes[3] = 2.68152185007253681141e-01;
                        m_Nodes[4] = 3.41994090825758473007e-01;
                        m_Nodes[5] = 4.13779204371605001525e-01;
                        m_Nodes[6] = 4.83075801686178712909e-01;
                        m_Nodes[7] = 5.49467125095128202076e-01;
                        m_Nodes[8] = 6.12553889667980237953e-01;
                        m_Nodes[9] = 6.71956684614179548379e-01;
                        m_Nodes[10] = 7.27318255189927103281e-01;
                        m_Nodes[11] = 7.78305651426519387695e-01;
                        m_Nodes[12] = 8.24612230833311663196e-01;
                        m_Nodes[13] = 8.65959503212259503821e-01;
                        m_Nodes[14] = 9.02098806968874296728e-01;
                        m_Nodes[15] = 9.32812808278676533361e-01;
                        m_Nodes[16] = 9.57916819213791655805e-01;
                        m_Nodes[17] = 9.77259949983774262663e-01;
                        m_Nodes[18] = 9.90726238699457006453e-01;
                        m_Nodes[19] = 9.98237709710559200350e-01;
                        m_Weights[0] = 7.75059479784248112637e-02;
                        m_Weights[1] = 7.70398181642479655883e-02;
                        m_Weights[2] = 7.61103619006262423716e-02;
                        m_Weights[3] = 7.47231690579682642002e-02;
                        m_Weights[4] = 7.28865823958040590605e-02;
                        m_Weights[5] = 7.06116473912867796955e-02;
                        m_Weights[6] = 6.79120458152339038257e-02;
                        m_Weights[7] = 6.48040134566010380745e-02;
                        m_Weights[8] = 6.13062424929289391669e-02;
                        m_Weights[9] = 5.74397690993915513708e-02;
                        m_Weights[10] = 5.32278469839368243467e-02;
                        m_Weights[11] = 4.86958076350722320979e-02;
                        m_Weights[12] = 4.38709081856732718647e-02;
                        m_Weights[13] = 3.87821679744720174194e-02;
                        m_Weights[14] = 3.34601952825478475071e-02;
                        m_Weights[15] = 2.79370069800234013281e-02;
                        m_Weights[16] = 2.22458491941669581622e-02;
                        m_Weights[17] = 1.64210583819078897106e-02;
                        m_Weights[18] = 1.04982845311528138790e-02;
                        m_Weights[19] = 4.52127709853319125116e-03;
                        break;

                    case 41:
                        m_Nodes[0] = 0.00000000000000000000e+00;
                        m_Nodes[1] = 7.56232589891629969238e-02;
                        m_Nodes[2] = 1.50813354863992163574e-01;
                        m_Nodes[3] = 2.25139605633422775606e-01;
                        m_Nodes[4] = 2.98176277341824865923e-01;
                        m_Nodes[5] = 3.69505022640481441428e-01;
                        m_Nodes[6] = 4.38717277051407088517e-01;
                        m_Nodes[7] = 5.05416599199406032708e-01;
                        m_Nodes[8] = 5.69220941610215869655e-01;
                        m_Nodes[9] = 6.29764839072196320489e-01;
                        m_Nodes[10] = 6.86701502034951289585e-01;
                        m_Nodes[11] = 7.39704803069926181060e-01;
                        m_Nodes[12] = 7.88471145047409372736e-01;
                        m_Nodes[13] = 8.32721200401361331244e-01;
                        m_Nodes[14] = 8.72201511692441408834e-01;
                        m_Nodes[15] = 9.06685944758101172958e-01;
                        m_Nodes[16] = 9.35976987497853825682e-01;
                        m_Nodes[17] = 9.59906891730346226099e-01;
                        m_Nodes[18] = 9.78338673561083384469e-01;
                        m_Nodes[19] = 9.91167109699016308250e-01;
                        m_Nodes[20] = 9.98321588574771441519e-01;
                        m_Weights[0] = 7.56955356472983723188e-02;
                        m_Weights[1] = 7.54787470927158240272e-02;
                        m_Weights[2] = 7.48296231762215518913e-02;
                        m_Weights[3] = 7.37518820272234699393e-02;
                        m_Weights[4] = 7.22516968610230733963e-02;
                        m_Weights[5] = 7.03376606208174974817e-02;
                        m_Weights[6] = 6.80207367608767667355e-02;
                        m_Weights[7] = 6.53141964535274104361e-02;
                        m_Weights[8] = 6.22335425809663164717e-02;
                        m_Weights[9] = 5.87964209498719449964e-02;
                        m_Weights[10] = 5.50225192425787418718e-02;
                        m_Weights[11] = 5.09334542946174948026e-02;
                        m_Weights[12] = 4.65526483690143421792e-02;
                        m_Weights[13] = 4.19051951959096898003e-02;
                        m_Weights[14] = 3.70177167035079892704e-02;
                        m_Weights[15] = 3.19182117316992838204e-02;
                        m_Weights[16] = 2.66358992071104385887e-02;
                        m_Weights[17] = 2.12010633687795562621e-02;
                        m_Weights[18] = 1.56449384078185855775e-02;
                        m_Weights[19] = 9.99993877390594649668e-03;
                        m_Weights[20] = 4.30614035816488777136e-03;
                        break;

                    case 42:
                        m_Nodes[0] = 3.69489431653517758131e-02;
                        m_Nodes[1] = 1.10645027208519868349e-01;
                        m_Nodes[2] = 1.83736806564854550853e-01;
                        m_Nodes[3] = 2.55825079342879083966e-01;
                        m_Nodes[4] = 3.26516124465411512197e-01;
                        m_Nodes[5] = 3.95423852042975057677e-01;
                        m_Nodes[6] = 4.62171912070421929759e-01;
                        m_Nodes[7] = 5.26395749931192287593e-01;
                        m_Nodes[8] = 5.87744597485109322841e-01;
                        m_Nodes[9] = 6.45883388869247833957e-01;
                        m_Nodes[10] = 7.00494590556171213742e-01;
                        m_Nodes[11] = 7.51279935689480489568e-01;
                        m_Nodes[12] = 7.97962053255487413233e-01;
                        m_Nodes[13] = 8.40285983261816900925e-01;
                        m_Nodes[14] = 8.78020569812172742712e-01;
                        m_Nodes[15] = 9.10959724904127452584e-01;
                        m_Nodes[16] = 9.38923557354988178533e-01;
                        m_Nodes[17] = 9.61759365338204488747e-01;
                        m_Nodes[18] = 9.79342508063748193709e-01;
                        m_Nodes[19] = 9.91577288340860919792e-01;
                        m_Nodes[20] = 9.98399618990062415023e-01;
                        m_Weights[0] = 7.38642342321728799964e-02;
                        m_Weights[1] = 7.34608134534675282640e-02;
                        m_Weights[2] = 7.26561752438041048879e-02;
                        m_Weights[3] = 7.14547142651709829218e-02;
                        m_Weights[4] = 6.98629924925941597662e-02;
                        m_Weights[5] = 6.78897033765219448554e-02;
                        m_Weights[6] = 6.55456243649089789270e-02;
                        m_Weights[7] = 6.28435580450025764094e-02;
                        m_Weights[8] = 5.97982622275866543126e-02;
                        m_Weights[9] = 5.64263693580183816421e-02;
                        m_Weights[10] = 5.27462956991740703803e-02;
                        m_Weights[11] = 4.87781407928032449864e-02;
                        m_Weights[12] = 4.45435777719658779666e-02;
                        m_Weights[13] = 4.00657351806922605131e-02;
                        m_Weights[14] = 3.53690710975921099161e-02;
                        m_Weights[15] = 3.04792406996034699827e-02;
                        m_Weights[16] = 2.54229595261130526284e-02;
                        m_Weights[17] = 2.02278695690526458464e-02;
                        m_Weights[18] = 1.49224436973574973475e-02;
                        m_Weights[19] = 9.53622030174850286204e-03;
                        m_Weights[20] = 4.10599860464908460589e-03;
                        break;

                    case 43:
                        m_Nodes[0] = 0.00000000000000000000e+00;
                        m_Nodes[1] = 7.21529908745862354223e-02;
                        m_Nodes[2] = 1.43929809510713310770e-01;
                        m_Nodes[3] = 2.14956244860518209015e-01;
                        m_Nodes[4] = 2.84861998032913627106e-01;
                        m_Nodes[5] = 3.53282612864303806645e-01;
                        m_Nodes[6] = 4.19861376029269252487e-01;
                        m_Nodes[7] = 4.84251176785734724070e-01;
                        m_Nodes[8] = 5.46116316660084719140e-01;
                        m_Nodes[9] = 6.05134259639600935725e-01;
                        m_Nodes[10] = 6.60997313751498133165e-01;
                        m_Nodes[11] = 7.13414235268957054852e-01;
                        m_Nodes[12] = 7.62111747194955121460e-01;
                        m_Nodes[13] = 8.06835964136938635279e-01;
                        m_Nodes[14] = 8.47353716209315048999e-01;
                        m_Nodes[15] = 8.83453765218616863338e-01;
                        m_Nodes[16] = 9.14947907206138729456e-01;
                        m_Nodes[17] = 9.41671956847637861818e-01;
                        m_Nodes[18] = 9.63486613014079993410e-01;
                        m_Nodes[19] = 9.80278220980255331506e-01;
                        m_Nodes[20] = 9.91959557593244146421e-01;
                        m_Nodes[21] = 9.98472332242507713518e-01;
                        m_Weights[0] = 7.22157516937989879775e-02;
                        m_Weights[1] = 7.20275019714219743453e-02;
                        m_Weights[2] = 7.14637342525141412976e-02;
                        m_Weights[3] = 7.05273877650850281263e-02;
                        m_Weights[4] = 6.92233441936566842823e-02;
                        m_Weights[5] = 6.75584022293651691924e-02;
                        m_Weights[6] = 6.55412421263227974912e-02;
                        m_Weights[7] = 6.31823804493961123256e-02;
                        m_Weights[8] = 6.04941152499912945198e-02;
                        m_Weights[9] = 5.74904619569105194291e-02;
                        m_Weights[10] = 5.41870803188817868642e-02;
                        m_Weights[11] = 5.06011927843901564769e-02;
                        m_Weights[12] = 4.67514947543465801070e-02;
                        m_Weights[13] = 4.26580571979820822240e-02;
                        m_Weights[14] = 3.83422221941326565420e-02;
                        m_Weights[15] = 3.38264920868602934092e-02;
                        m_Weights[16] = 2.91344132614984979310e-02;
                        m_Weights[17] = 2.42904566138388411664e-02;
                        m_Weights[18] = 1.93199014236839051659e-02;
                        m_Weights[19] = 1.42487564315764852451e-02;
                        m_Weights[20] = 9.10399663740139727940e-03;
                        m_Weights[21] = 3.91949025384412795798e-03;
                        break;

                    case 44:
                        m_Nodes[0] = 3.52892369641353590582e-02;
                        m_Nodes[1] = 1.05691901708653247117e-01;
                        m_Nodes[2] = 1.75568014775516785747e-01;
                        m_Nodes[3] = 2.44569456928201251507e-01;
                        m_Nodes[4] = 3.12352466502785812237e-01;
                        m_Nodes[5] = 3.78579352014707132512e-01;
                        m_Nodes[6] = 4.42920174525411483835e-01;
                        m_Nodes[7] = 5.05054391388202317983e-01;
                        m_Nodes[8] = 5.64672453185470768425e-01;
                        m_Nodes[9] = 6.21477345903575847802e-01;
                        m_Nodes[10] = 6.75186070666122365334e-01;
                        m_Nodes[11] = 7.25531053660717002607e-01;
                        m_Nodes[12] = 7.72261479248755899018e-01;
                        m_Nodes[13] = 8.15144539645135010487e-01;
                        m_Nodes[14] = 8.53966595004710378728e-01;
                        m_Nodes[15] = 8.88534238286043202338e-01;
                        m_Nodes[16] = 9.18675259984175774323e-01;
                        m_Nodes[17] = 9.44239509118194099203e-01;
                        m_Nodes[18] = 9.65099650422493139394e-01;
                        m_Nodes[19] = 9.81151833077913966663e-01;
                        m_Nodes[20] = 9.92316392138515808483e-01;
                        m_Nodes[21] = 9.98540200636774224936e-01;
                        m_Weights[0] = 7.05491577893540688113e-02;
                        m_Weights[1] = 7.01976854735582125871e-02;
                        m_Weights[2] = 6.94964918615725780371e-02;
                        m_Weights[3] = 6.84490702693666609855e-02;
                        m_Weights[4] = 6.70606389062936523957e-02;
                        m_Weights[5] = 6.53381148791814349842e-02;
                        m_Weights[6] = 6.32900797332038549501e-02;
                        m_Weights[7] = 6.09267367015619680385e-02;
                        m_Weights[8] = 5.82598598775954953347e-02;
                        m_Weights[9] = 5.53027355637280525508e-02;
                        m_Weights[10] = 5.20700960917044618268e-02;
                        m_Weights[11] = 4.85780464483520376037e-02;
                        m_Weights[12] = 4.48439840819700310701e-02;
                        m_Weights[13] = 4.08865123103462166629e-02;
                        m_Weights[14] = 3.67253478138088881245e-02;
                        m_Weights[15] = 3.23812228120698355718e-02;
                        m_Weights[16] = 2.78757828212809970438e-02;
                        m_Weights[17] = 2.32314819020192005843e-02;
                        m_Weights[18] = 1.84714817368147974507e-02;
                        m_Weights[19] = 1.36195867555799710634e-02;
                        m_Weights[20] = 8.70048136752485912194e-03;
                        m_Weights[21] = 3.74540480311277660786e-03;
                        break;

                    case 45:
                        m_Nodes[0] = 0.00000000000000000000e+00;
                        m_Nodes[1] = 6.89869801631441724904e-02;
                        m_Nodes[2] = 1.37645205983253028757e-01;
                        m_Nodes[3] = 2.05647489783263745720e-01;
                        m_Nodes[4] = 2.72669769752377560609e-01;
                        m_Nodes[5] = 3.38392654250602161643e-01;
                        m_Nodes[6] = 4.02502943858541914078e-01;
                        m_Nodes[7] = 4.64695123919635098580e-01;
                        m_Nodes[8] = 5.24672820462916067091e-01;
                        m_Nodes[9] = 5.82150212569353186681e-01;
                        m_Nodes[10] = 6.36853394453223359271e-01;
                        m_Nodes[11] = 6.88521680771200525232e-01;
                        m_Nodes[12] = 7.36908848945490352624e-01;
                        m_Nodes[13] = 7.81784312593906291312e-01;
                        m_Nodes[14] = 8.22934220502086337036e-01;
                        m_Nodes[15] = 8.60162475960664225339e-01;
                        m_Nodes[16] = 8.93291671753241738465e-01;
                        m_Nodes[17] = 9.22163936719000388097e-01;
                        m_Nodes[18] = 9.46641690995629061785e-01;
                        m_Nodes[19] = 9.66608310396894604736e-01;
                        m_Nodes[20] = 9.81968715034540568239e-01;
                        m_Nodes[21] = 9.92649998447203741749e-01;
                        m_Nodes[22] = 9.98603645181936638157e-01;
                        m_Weights[0] = 6.90418248292320201108e-02;
                        m_Weights[1] = 6.88773169776613228820e-02;
                        m_Weights[2] = 6.83845773786696745317e-02;
                        m_Weights[3] = 6.75659541636075362709e-02;
                        m_Weights[4] = 6.64253484498425280829e-02;
                        m_Weights[5] = 6.49681957507234308538e-02;
                        m_Weights[6] = 6.32014400738199377500e-02;
                        m_Weights[7] = 6.11335008310665225018e-02;
                        m_Weights[8] = 5.87742327188417385745e-02;
                        m_Weights[9] = 5.61348787597864766443e-02;
                        m_Weights[10] = 5.32280167312689519474e-02;
                        m_Weights[11] = 5.00674992379520299838e-02;
                        m_Weights[12] = 4.66683877183733651211e-02;
                        m_Weights[13] = 4.30468807091649695019e-02;
                        m_Weights[14] = 3.92202367293024536267e-02;
                        m_Weights[15] = 3.52066922016090273831e-02;
                        m_Weights[16] = 3.10253749345154645249e-02;
                        m_Weights[17] = 2.66962139675777466475e-02;
                        m_Weights[18] = 2.22398475505789082868e-02;
                        m_Weights[19] = 1.76775352579376024562e-02;
                        m_Weights[20] = 1.30311049915827421078e-02;
                        m_Weights[21] = 8.32318929621823236028e-03;
                        m_Weights[22] = 3.58266315528355871619e-03;
                        break;

                    case 46:
                        m_Nodes[0] = 3.37721900160520415196e-02;
                        m_Nodes[1] = 1.01162475305584239516e-01;
                        m_Nodes[2] = 1.68091179467103528607e-01;
                        m_Nodes[3] = 2.34252922206269768626e-01;
                        m_Nodes[4] = 2.99345822701870015483e-01;
                        m_Nodes[5] = 3.63072877020995710124e-01;
                        m_Nodes[6] = 4.25143313282828397322e-01;
                        m_Nodes[7] = 4.85273918388164662772e-01;
                        m_Nodes[8] = 5.43190330261802635271e-01;
                        m_Nodes[9] = 5.98628289712715153177e-01;
                        m_Nodes[10] = 6.51334846201997715106e-01;
                        m_Nodes[11] = 7.01069512020405697512e-01;
                        m_Nodes[12] = 7.47605359615666054000e-01;
                        m_Nodes[13] = 7.90730057075274255189e-01;
                        m_Nodes[14] = 8.30246837066066053032e-01;
                        m_Nodes[15] = 8.65975394866858062916e-01;
                        m_Nodes[16] = 8.97752711533941965701e-01;
                        m_Nodes[17] = 9.25433798806753950977e-01;
                        m_Nodes[18] = 9.48892363446089795622e-01;
                        m_Nodes[19] = 9.68021391853991942738e-01;
                        m_Nodes[20] = 9.82733669804166863478e-01;
                        m_Nodes[21] = 9.92962348906174364073e-01;
                        m_Nodes[22] = 9.98663042133817981128e-01;
                        m_Weights[0] = 6.75186858490364588202e-02;
                        m_Weights[1] = 6.72106136006781758624e-02;
                        m_Weights[2] = 6.65958747684548873758e-02;
                        m_Weights[3] = 6.56772742677812073788e-02;
                        m_Weights[4] = 6.44590034671390695883e-02;
                        m_Weights[5] = 6.29466210643945081790e-02;
                        m_Weights[6] = 6.11470277246504810153e-02;
                        m_Weights[7] = 5.90684345955463148074e-02;
                        m_Weights[8] = 5.67203258439912358172e-02;
                        m_Weights[9] = 5.41134153858567544918e-02;
                        m_Weights[10] = 5.12595980071430213740e-02;
                        m_Weights[11] = 4.81718951017122007472e-02;
                        m_Weights[12] = 4.48643952773181274887e-02;
                        m_Weights[13] = 4.13521901096787299843e-02;
                        m_Weights[14] = 3.76513053573860705015e-02;
                        m_Weights[15] = 3.37786279991069746453e-02;
                        m_Weights[16] = 2.97518295522027144815e-02;
                        m_Weights[17] = 2.55892863971300274571e-02;
                        m_Weights[18] = 2.13099987541363110043e-02;
                        m_Weights[19] = 1.69335140078360420565e-02;
                        m_Weights[20] = 1.24798837709885947975e-02;
                        m_Weights[21] = 7.96989822972459412567e-03;
                        m_Weights[22] = 3.43030086810705304333e-03;
                        break;

                    case 47:
                        m_Nodes[0] = 0.00000000000000000000e+00;
                        m_Nodes[1] = 6.60869239163556751605e-02;
                        m_Nodes[2] = 1.31884866554514897054e-01;
                        m_Nodes[3] = 1.97106110279111807961e-01;
                        m_Nodes[4] = 2.61465459214974570307e-01;
                        m_Nodes[5] = 3.24681486337735902211e-01;
                        m_Nodes[6] = 3.86477764084667139583e-01;
                        m_Nodes[7] = 4.46584073104855702725e-01;
                        m_Nodes[8] = 5.04737583863577919774e-01;
                        m_Nodes[9] = 5.60684005934664194483e-01;
                        m_Nodes[10] = 6.14178699956373608595e-01;
                        m_Nodes[11] = 6.64987747390332729137e-01;
                        m_Nodes[12] = 7.12888973409064301662e-01;
                        m_Nodes[13] = 7.57672918445438633574e-01;
                        m_Nodes[14] = 7.99143754167741942916e-01;
                        m_Nodes[15] = 8.37120139899902121278e-01;
                        m_Nodes[16] = 8.71436015796896316941e-01;
                        m_Nodes[17] = 9.01941329438525356867e-01;
                        m_Nodes[18] = 9.28502693012360648197e-01;
                        m_Nodes[19] = 9.51003969257708442590e-01;
                        m_Nodes[20] = 9.69346787326564497146e-01;
                        m_Nodes[21] = 9.83451003071623708765e-01;
                        m_Nodes[22] = 9.93255210987768634692e-01;
                        m_Nodes[23] = 9.98718728584212109184e-01;
                        m_Weights[0] = 6.61351296236554796534e-02;
                        m_Weights[1] = 6.59905335888104745336e-02;
                        m_Weights[2] = 6.55573777665497402511e-02;
                        m_Weights[3] = 6.48375562389457267026e-02;
                        m_Weights[4] = 6.38342166057170306313e-02;
                        m_Weights[5] = 6.25517462209216626406e-02;
                        m_Weights[6] = 6.09957530087396453307e-02;
                        m_Weights[7] = 5.91730409423388759762e-02;
                        m_Weights[8] = 5.70915802932315402218e-02;
                        m_Weights[9] = 5.47604727815302259581e-02;
                        m_Weights[10] = 5.21899117800571448771e-02;
                        m_Weights[11] = 4.93911377473611694569e-02;
                        m_Weights[12] = 4.63763890865059104698e-02;
                        m_Weights[13] = 4.31588486484795406705e-02;
                        m_Weights[14] = 3.97525861225310060912e-02;
                        m_Weights[15] = 3.61724965841750289423e-02;
                        m_Weights[16] = 3.24342355151846521706e-02;
                        m_Weights[17] = 2.85541507006435212618e-02;
                        m_Weights[18] = 2.45492116596588758323e-02;
                        m_Weights[19] = 2.04369381476685138064e-02;
                        m_Weights[20] = 1.62353331464331174243e-02;
                        m_Weights[21] = 1.19628484643122047109e-02;
                        m_Weights[22] = 7.63861629584886345853e-03;
                        m_Weights[23] = 3.28745384252803037776e-03;
                        break;

                    case 48:
                        m_Nodes[0] = 3.23801709628693620333e-02;
                        m_Nodes[1] = 9.70046992094626989301e-02;
                        m_Nodes[2] = 1.61222356068891718056e-01;
                        m_Nodes[3] = 2.24763790394689061225e-01;
                        m_Nodes[4] = 2.87362487355455576736e-01;
                        m_Nodes[5] = 3.48755886292160738160e-01;
                        m_Nodes[6] = 4.08686481990716729916e-01;
                        m_Nodes[7] = 4.66902904750958404545e-01;
                        m_Nodes[8] = 5.23160974722233033678e-01;
                        m_Nodes[9] = 5.77224726083972703818e-01;
                        m_Nodes[10] = 6.28867396776513623995e-01;
                        m_Nodes[11] = 6.77872379632663905212e-01;
                        m_Nodes[12] = 7.24034130923814654674e-01;
                        m_Nodes[13] = 7.67159032515740339254e-01;
                        m_Nodes[14] = 8.07066204029442627083e-01;
                        m_Nodes[15] = 8.43588261624393530711e-01;
                        m_Nodes[16] = 8.76572020274247885906e-01;
                        m_Nodes[17] = 9.05879136715569672822e-01;
                        m_Nodes[18] = 9.31386690706554333114e-01;
                        m_Nodes[19] = 9.52987703160430860723e-01;
                        m_Nodes[20] = 9.70591592546247250461e-01;
                        m_Nodes[21] = 9.84124583722826857745e-01;
                        m_Nodes[22] = 9.93530172266350757548e-01;
                        m_Nodes[23] = 9.98771007252426118601e-01;
                        m_Weights[0] = 6.47376968126839225030e-02;
                        m_Weights[1] = 6.44661644359500822065e-02;
                        m_Weights[2] = 6.39242385846481866239e-02;
                        m_Weights[3] = 6.31141922862540256571e-02;
                        m_Weights[4] = 6.20394231598926639042e-02;
                        m_Weights[5] = 6.07044391658938800530e-02;
                        m_Weights[6] = 5.91148396983956357465e-02;
                        m_Weights[7] = 5.72772921004032157050e-02;
                        m_Weights[8] = 5.51995036999841628678e-02;
                        m_Weights[9] = 5.28901894851936670763e-02;
                        m_Weights[10] = 5.03590355538544749253e-02;
                        m_Weights[11] = 4.76166584924904748974e-02;
                        m_Weights[12] = 4.46745608566942823776e-02;
                        m_Weights[13] = 4.15450829434647598352e-02;
                        m_Weights[14] = 3.82413510658307487254e-02;
                        m_Weights[15] = 3.47772225647702792702e-02;
                        m_Weights[16] = 3.11672278327979467870e-02;
                        m_Weights[17] = 2.74265097083572354059e-02;
                        m_Weights[18] = 2.35707608393248588637e-02;
                        m_Weights[19] = 1.96161604573556160746e-02;
                        m_Weights[20] = 1.55793157229442039210e-02;
                        m_Weights[21] = 1.14772345792343483453e-02;
                        m_Weights[22] = 7.32755390127637745050e-03;
                        m_Weights[23] = 3.15334605230574971921e-03;
                        break;

                    case 49:
                        m_Nodes[0] = 0.00000000000000000000e+00;
                        m_Nodes[1] = 6.34206849826867860288e-02;
                        m_Nodes[2] = 1.26585997269672051068e-01;
                        m_Nodes[3] = 1.89241592461813586485e-01;
                        m_Nodes[4] = 2.51135178612577273507e-01;
                        m_Nodes[5] = 3.12017532119748762208e-01;
                        m_Nodes[6] = 3.71643501262284888864e-01;
                        m_Nodes[7] = 4.29772993341576524659e-01;
                        m_Nodes[8] = 4.86171941452492042177e-01;
                        m_Nodes[9] = 5.40613246991726066558e-01;
                        m_Nodes[10] = 5.92877694108900712456e-01;
                        m_Nodes[11] = 6.42754832419237664057e-01;
                        m_Nodes[12] = 6.90043824425132113505e-01;
                        m_Nodes[13] = 7.34554254237402696214e-01;
                        m_Nodes[14] = 7.76106894345446635018e-01;
                        m_Nodes[15] = 8.14534427359855431540e-01;
                        m_Nodes[16] = 8.49682119844165701035e-01;
                        m_Nodes[17] = 8.81408445573008910037e-01;
                        m_Nodes[18] = 9.09585655828073285213e-01;
                        m_Nodes[19] = 9.34100294755810149059e-01;
                        m_Nodes[20] = 9.54853658674137233555e-01;
                        m_Nodes[21] = 9.71762200901555380140e-01;
                        m_Nodes[22] = 9.84757895914213004359e-01;
                        m_Nodes[23] = 9.93788661944167790760e-01;
                        m_Nodes[24] = 9.98820150606635379362e-01;
                        m_Weights[0] = 6.34632814047905977183e-02;
                        m_Weights[1] = 6.33355092964917485908e-02;
                        m_Weights[2] = 6.29527074651956994744e-02;
                        m_Weights[3] = 6.23164173200572674011e-02;
                        m_Weights[4] = 6.14292009791929362968e-02;
                        m_Weights[5] = 6.02946309531520173031e-02;
                        m_Weights[6] = 5.89172757600272660245e-02;
                        m_Weights[7] = 5.73026815301874754852e-02;
                        m_Weights[8] = 5.54573496748035886902e-02;
                        m_Weights[9] = 5.33887107082589685293e-02;
                        m_Weights[10] = 5.11050943301445906822e-02;
                        m_Weights[11] = 4.86156958878282402586e-02;
                        m_Weights[12] = 4.59305393555958535590e-02;
                        m_Weights[13] = 4.30604369812596001318e-02;
                        m_Weights[14] = 4.00169457663730359297e-02;
                        m_Weights[15] = 3.68123209630006975084e-02;
                        m_Weights[16] = 3.34594667916225282493e-02;
                        m_Weights[17] = 2.99718846205841918967e-02;
                        m_Weights[18] = 2.63636189270667827533e-02;
                        m_Weights[19] = 2.26492015874478649166e-02;
                        m_Weights[20] = 1.88435958530910878406e-02;
                        m_Weights[21] = 1.49621449356250475306e-02;
                        m_Weights[22] = 1.10205510315889101049e-02;
                        m_Weights[23] = 7.03509959008655034763e-03;
                        m_Weights[24] = 3.02727898892302399310e-03;
                        break;

                    case 50:
                        m_Nodes[0] = 3.10983383271888761123e-02;
                        m_Nodes[1] = 9.31747015600861408545e-02;
                        m_Nodes[2] = 1.54890589998145902072e-01;
                        m_Nodes[3] = 2.16007236876041756847e-01;
                        m_Nodes[4] = 2.76288193779531990328e-01;
                        m_Nodes[5] = 3.35500245419437356837e-01;
                        m_Nodes[6] = 3.93414311897565127394e-01;
                        m_Nodes[7] = 4.49806334974038789147e-01;
                        m_Nodes[8] = 5.04458144907464201651e-01;
                        m_Nodes[9] = 5.57158304514650054316e-01;
                        m_Nodes[10] = 6.07702927184950239180e-01;
                        m_Nodes[11] = 6.55896465685439360782e-01;
                        m_Nodes[12] = 7.01552468706822251090e-01;
                        m_Nodes[13] = 7.44494302226068538261e-01;
                        m_Nodes[14] = 7.84555832900399263905e-01;
                        m_Nodes[15] = 8.21582070859335948356e-01;
                        m_Nodes[16] = 8.55429769429946084611e-01;
                        m_Nodes[17] = 8.85967979523613048638e-01;
                        m_Nodes[18] = 9.13078556655791893090e-01;
                        m_Nodes[19] = 9.36656618944877933781e-01;
                        m_Nodes[20] = 9.56610955242807942998e-01;
                        m_Nodes[21] = 9.72864385106692073713e-01;
                        m_Nodes[22] = 9.85354084048005882309e-01;
                        m_Nodes[23] = 9.94031969432090712585e-01;
                        m_Nodes[24] = 9.98866404420071050185e-01;
                        m_Weights[0] = 6.21766166553472623210e-02;
                        m_Weights[1] = 6.19360674206832433841e-02;
                        m_Weights[2] = 6.14558995903166637564e-02;
                        m_Weights[3] = 6.07379708417702160318e-02;
                        m_Weights[4] = 5.97850587042654575096e-02;
                        m_Weights[5] = 5.86008498132224458351e-02;
                        m_Weights[6] = 5.71899256477283837230e-02;
                        m_Weights[7] = 5.55577448062125176235e-02;
                        m_Weights[8] = 5.37106218889962465229e-02;
                        m_Weights[9] = 5.16557030695811384807e-02;
                        m_Weights[10] = 4.94009384494663148878e-02;
                        m_Weights[11] = 4.69550513039484321600e-02;
                        m_Weights[12] = 4.43275043388032746178e-02;
                        m_Weights[13] = 4.15284630901477122165e-02;
                        m_Weights[14] = 3.85687566125876245635e-02;
                        m_Weights[15] = 3.54598356151459885886e-02;
                        m_Weights[16] = 3.22137282235775828392e-02;
                        m_Weights[17] = 2.88429935805359982379e-02;
                        m_Weights[18] = 2.53606735700112212159e-02;
                        m_Weights[19] = 2.17802431701224471060e-02;
                        m_Weights[20] = 1.81155607134922973158e-02;
                        m_Weights[21] = 1.43808227614832769517e-02;
                        m_Weights[22] = 1.05905483836497511811e-02;
                        m_Weights[23] = 6.75979919574541553150e-03;
                        m_Weights[24] = 2.90862255315546895921e-03;
                        break;

                    case 51:
                        m_Nodes[0] = 0.00000000000000000000e+00;
                        m_Nodes[1] = 6.09611001505787247342e-02;
                        m_Nodes[2] = 1.21695421018888766964e-01;
                        m_Nodes[3] = 1.81977026957077545324e-01;
                        m_Nodes[4] = 2.41581666447798703847e-01;
                        m_Nodes[5] = 3.00287606335331939530e-01;
                        m_Nodes[6] = 3.57876456688409509775e-01;
                        m_Nodes[7] = 4.14133983226303877937e-01;
                        m_Nodes[8] = 4.68850904286041063610e-01;
                        m_Nodes[9] = 5.21823669366185842514e-01;
                        m_Nodes[10] = 5.72855216351303836522e-01;
                        m_Nodes[11] = 6.21755704600723273755e-01;
                        m_Nodes[12] = 6.68343221175370086864e-01;
                        m_Nodes[13] = 7.12444457577036644581e-01;
                        m_Nodes[14] = 7.53895354485375525764e-01;
                        m_Nodes[15] = 7.92541712099381205234e-01;
                        m_Nodes[16] = 8.28239763823064832855e-01;
                        m_Nodes[17] = 8.60856711182292371473e-01;
                        m_Nodes[18] = 8.90271218029527303278e-01;
                        m_Nodes[19] = 9.16373862309780230824e-01;
                        m_Nodes[20] = 9.39067544002962383435e-01;
                        m_Nodes[21] = 9.58267848613908194558e-01;
                        m_Nodes[22] = 9.73903368019323867232e-01;
                        m_Nodes[23] = 9.85915991735902996584e-01;
                        m_Nodes[24] = 9.94261260436752574621e-01;
                        m_Nodes[25] = 9.98909990848903495169e-01;
                        m_Weights[0] = 6.09989248412058801598e-02;
                        m_Weights[1] = 6.08854648448563438812e-02;
                        m_Weights[2] = 6.05455069347377951381e-02;
                        m_Weights[3] = 5.99803157775032520901e-02;
                        m_Weights[4] = 5.91919939229615437835e-02;
                        m_Weights[5] = 5.81834739825921405984e-02;
                        m_Weights[6] = 5.69585077202586621001e-02;
                        m_Weights[7] = 5.55216520957386930167e-02;
                        m_Weights[8] = 5.38782523130455614345e-02;
                        m_Weights[9] = 5.20344219366970875689e-02;
                        m_Weights[10] = 4.99970201500574098088e-02;
                        m_Weights[11] = 4.77736262406231021432e-02;
                        m_Weights[12] = 4.53725114076500680896e-02;
                        m_Weights[13] = 4.28026079978800929731e-02;
                        m_Weights[14] = 4.00734762854964951904e-02;
                        m_Weights[15] = 3.71952689232601606778e-02;
                        m_Weights[16] = 3.41786932041877154352e-02;
                        m_Weights[17] = 3.10349712901594219890e-02;
                        m_Weights[18] = 2.77757985941611661089e-02;
                        m_Weights[19] = 2.44133005737776275450e-02;
                        m_Weights[20] = 2.09599884017101241569e-02;
                        m_Weights[21] = 1.74287147234191760313e-02;
                        m_Weights[22] = 1.38326340064700654811e-02;
                        m_Weights[23] = 1.01851912978255691539e-02;
                        m_Weights[24] = 6.50033778325918863879e-03;
                        m_Weights[25] = 2.79680717108977765296e-03;
                        break;

                    case 52:
                        m_Nodes[0] = 2.99141097973387660437e-02;
                        m_Nodes[1] = 8.96352446489005654889e-02;
                        m_Nodes[2] = 1.49035508606949180489e-01;
                        m_Nodes[3] = 2.07902264156366059686e-01;
                        m_Nodes[4] = 2.66024783605001827473e-01;
                        m_Nodes[5] = 3.23195003434807825501e-01;
                        m_Nodes[6] = 3.79208269116093669247e-01;
                        m_Nodes[7] = 4.33864067718761670309e-01;
                        m_Nodes[8] = 4.86966745698096077782e-01;
                        m_Nodes[9] = 5.38326209285827438376e-01;
                        m_Nodes[10] = 5.87758604979579069902e-01;
                        m_Nodes[11] = 6.35086977695245924298e-01;
                        m_Nodes[12] = 6.80141904227167702092e-01;
                        m_Nodes[13] = 7.22762099749983193677e-01;
                        m_Nodes[14] = 7.62794995193744960279e-01;
                        m_Nodes[15] = 8.00097283430468324335e-01;
                        m_Nodes[16] = 8.34535432326734534962e-01;
                        m_Nodes[17] = 8.65986162846067585244e-01;
                        m_Nodes[18] = 8.94336890534495322521e-01;
                        m_Nodes[19] = 9.19486128916424539894e-01;
                        m_Nodes[20] = 9.41343853641359056844e-01;
                        m_Nodes[21] = 9.59831826933086552532e-01;
                        m_Nodes[22] = 9.74883884221744503141e-01;
                        m_Nodes[23] = 9.86446195651549840645e-01;
                        m_Nodes[24] = 9.94477590929216029245e-01;
                        m_Nodes[25] = 9.98951111103950278091e-01;
                        m_Weights[0] = 5.98103657452918602478e-02;
                        m_Weights[1] = 5.95962601712481582583e-02;
                        m_Weights[2] = 5.91688154660429703693e-02;
                        m_Weights[3] = 5.85295617718138685503e-02;
                        m_Weights[4] = 5.76807874525268276539e-02;
                        m_Weights[5] = 5.66255309023685971908e-02;
                        m_Weights[6] = 5.53675696693026525491e-02;
                        m_Weights[7] = 5.39114069327572647508e-02;
                        m_Weights[8] = 5.22622553839069930353e-02;
                        m_Weights[9] = 5.04260185663423772385e-02;
                        m_Weights[10] = 4.84092697440748969829e-02;
                        m_Weights[11] = 4.62192283727847946823e-02;
                        m_Weights[12] = 4.38637342590004137384e-02;
                        m_Weights[13] = 4.13512195005602685293e-02;
                        m_Weights[14] = 3.86906783104240049477e-02;
                        m_Weights[15] = 3.58916348350973408570e-02;
                        m_Weights[16] = 3.29641090897185638392e-02;
                        m_Weights[17] = 2.99185811471436309955e-02;
                        m_Weights[18] = 2.67659537465134452626e-02;
                        m_Weights[19] = 2.35175135539915041807e-02;
                        m_Weights[20] = 2.01848915079238250403e-02;
                        m_Weights[21] = 1.67800233962651848517e-02;
                        m_Weights[22] = 1.33151149823097482618e-02;
                        m_Weights[23] = 9.80263457948108670461e-03;
                        m_Weights[24] = 6.25552396297798482625e-03;
                        m_Weights[25] = 2.69131695004765195184e-03;
                        break;

                    case 53:
                        m_Nodes[0] = 0.00000000000000000000e+00;
                        m_Nodes[1] = 5.86850543002594650227e-02;
                        m_Nodes[2] = 1.17167809071955150140e-01;
                        m_Nodes[3] = 1.75246662155325750730e-01;
                        m_Nodes[4] = 2.32721403724272593643e-01;
                        m_Nodes[5] = 2.89393906451626206427e-01;
                        m_Nodes[6] = 3.45068808495722356694e-01;
                        m_Nodes[7] = 3.99554186953952977393e-01;
                        m_Nodes[8] = 4.52662219461845791383e-01;
                        m_Nodes[9] = 5.04209831657133437039e-01;
                        m_Nodes[10] = 5.54019328277067881015e-01;
                        m_Nodes[11] = 6.01919005713769327464e-01;
                        m_Nodes[12] = 6.47743743916510068751e-01;
                        m_Nodes[13] = 6.91335575601366723541e-01;
                        m_Nodes[14] = 7.32544230807510253782e-01;
                        m_Nodes[15] = 7.71227654925532307866e-01;
                        m_Nodes[16] = 8.07252498416895478220e-01;
                        m_Nodes[17] = 8.40494576545801375430e-01;
                        m_Nodes[18] = 8.70839297558241351602e-01;
                        m_Nodes[19] = 8.98182057875426625926e-01;
                        m_Nodes[20] = 9.22428603042812128268e-01;
                        m_Nodes[21] = 9.43495353464441879021e-01;
                        m_Nodes[22] = 9.61309694623136332370e-01;
                        m_Nodes[23] = 9.75810233714984581633e-01;
                        m_Nodes[24] = 9.86947035023371521720e-01;
                        m_Nodes[25] = 9.94681919308007078636e-01;
                        m_Nodes[26] = 9.98989947776328227121e-01;
                        m_Weights[0] = 5.87187941511643645255e-02;
                        m_Weights[1] = 5.86175862327202633181e-02;
                        m_Weights[2] = 5.83143113622560075563e-02;
                        m_Weights[3] = 5.78100149917131963197e-02;
                        m_Weights[4] = 5.71064355362671917734e-02;
                        m_Weights[5] = 5.62059983817397098087e-02;
                        m_Weights[6] = 5.51118075239335990024e-02;
                        m_Weights[7] = 5.38276348687310290421e-02;
                        m_Weights[8] = 5.23579072298727181995e-02;
                        m_Weights[9] = 5.07076910692927152935e-02;
                        m_Weights[10] = 4.88826750326991403937e-02;
                        m_Weights[11] = 4.68891503407503139402e-02;
                        m_Weights[12] = 4.47339891036728137902e-02;
                        m_Weights[13] = 4.24246206345200091487e-02;
                        m_Weights[14] = 3.99690058435403972934e-02;
                        m_Weights[15] = 3.73756098034828501516e-02;
                        m_Weights[16] = 3.46533725835338223121e-02;
                        m_Weights[17] = 3.18116784590141299897e-02;
                        m_Weights[18] = 2.88603236178228598770e-02;
                        m_Weights[19] = 2.58094825107842735192e-02;
                        m_Weights[20] = 2.26696730570647297634e-02;
                        m_Weights[21] = 1.94517211077173879909e-02;
                        m_Weights[22] = 1.61667252567194161574e-02;
                        m_Weights[23] = 1.28260261443115041283e-02;
                        m_Weights[24] = 9.44120228492797865626e-03;
                        m_Weights[25] = 6.02427622695195442016e-03;
                        m_Weights[26] = 2.59168372056699694426e-03;
                        break;

                    case 54:
                        m_Nodes[0] = 2.88167481993417776562e-02;
                        m_Nodes[1] = 8.63545182632482152854e-02;
                        m_Nodes[2] = 1.43605427316256153947e-01;
                        m_Nodes[3] = 2.00379293606213569779e-01;
                        m_Nodes[4] = 2.56487520069997300077e-01;
                        m_Nodes[5] = 3.11743720834468228883e-01;
                        m_Nodes[6] = 3.65964340372191181984e-01;
                        m_Nodes[7] = 4.18969263255204528036e-01;
                        m_Nodes[8] = 4.70582412481382283683e-01;
                        m_Nodes[9] = 5.20632334385933073327e-01;
                        m_Nodes[10] = 5.68952768195209429732e-01;
                        m_Nodes[11] = 6.15383198331127370730e-01;
                        m_Nodes[12] = 6.59769387631983124692e-01;
                        m_Nodes[13] = 7.01963889719172919386e-01;
                        m_Nodes[14] = 7.41826538809184316286e-01;
                        m_Nodes[15] = 7.79224915346254021536e-01;
                        m_Nodes[16] = 8.14034785913567835470e-01;
                        m_Nodes[17] = 8.46140515970772949426e-01;
                        m_Nodes[18] = 8.75435454065568939418e-01;
                        m_Nodes[19] = 9.01822286284701580757e-01;
                        m_Nodes[20] = 9.25213359866651486256e-01;
                        m_Nodes[21] = 9.45530975164995853764e-01;
                        m_Nodes[22] = 9.62707645785923583257e-01;
                        m_Nodes[23] = 9.76686328857903237200e-01;
                        m_Nodes[24] = 9.87420637397343558552e-01;
                        m_Nodes[25] = 9.94875117018338884919e-01;
                        m_Nodes[26] = 9.99026666867340983851e-01;
                        m_Weights[0] = 5.76175367071470246724e-02;
                        m_Weights[1] = 5.74261370541121148593e-02;
                        m_Weights[2] = 5.70439735587945985678e-02;
                        m_Weights[3] = 5.64723157306259650310e-02;
                        m_Weights[4] = 5.57130625605899876834e-02;
                        m_Weights[5] = 5.47687362130579863062e-02;
                        m_Weights[6] = 5.36424736475536112721e-02;
                        m_Weights[7] = 5.23380161982987446656e-02;
                        m_Weights[8] = 5.08596971461881443198e-02;
                        m_Weights[9] = 4.92124273245288860762e-02;
                        m_Weights[10] = 4.74016788064449912244e-02;
                        m_Weights[11] = 4.54334667282767146538e-02;
                        m_Weights[12] = 4.33143293095970006374e-02;
                        m_Weights[13] = 4.10513061366450513072e-02;
                        m_Weights[14] = 3.86519147821024006323e-02;
                        m_Weights[15] = 3.61241258403831234439e-02;
                        m_Weights[16] = 3.34763364643717031746e-02;
                        m_Weights[17] = 3.07173424978658405181e-02;
                        m_Weights[18] = 2.78563093106382757926e-02;
                        m_Weights[19] = 2.49027414672489451000e-02;
                        m_Weights[20] = 2.18664514227868652117e-02;
                        m_Weights[21] = 1.87575276217524920516e-02;
                        m_Weights[22] = 1.55863030359834590434e-02;
                        m_Weights[23] = 1.23633281288816883635e-02;
                        m_Weights[24] = 9.09936945531942763902e-03;
                        m_Weights[25] = 5.80561101517794882020e-03;
                        m_Weights[26] = 2.49748183576593864607e-03;
                        break;

                    case 55:
                        m_Nodes[0] = 0.00000000000000000000e+00;
                        m_Nodes[1] = 5.65727538183367763273e-02;
                        m_Nodes[2] = 1.12964288059329266588e-01;
                        m_Nodes[3] = 1.68993963646873208283e-01;
                        m_Nodes[4] = 2.24482300647845483400e-01;
                        m_Nodes[5] = 2.79251553200806538550e-01;
                        m_Nodes[6] = 3.33126278890023885189e-01;
                        m_Nodes[7] = 3.85933900740979429756e-01;
                        m_Nodes[8] = 4.37505260037174591808e-01;
                        m_Nodes[9] = 4.87675158187474097208e-01;
                        m_Nodes[10] = 5.36282885908343296721e-01;
                        m_Nodes[11] = 5.83172738026032102974e-01;
                        m_Nodes[12] = 6.28194512249928140091e-01;
                        m_Nodes[13] = 6.71203990319826395796e-01;
                        m_Nodes[14] = 7.12063399986637838909e-01;
                        m_Nodes[15] = 7.50641856348021908675e-01;
                        m_Nodes[16] = 7.86815781127622365898e-01;
                        m_Nodes[17] = 8.20469298559320912454e-01;
                        m_Nodes[18] = 8.51494606617154471460e-01;
                        m_Nodes[19] = 8.79792322419895506068e-01;
                        m_Nodes[20] = 9.05271800744000025782e-01;
                        m_Nodes[21] = 9.27851424720791696816e-01;
                        m_Nodes[22] = 9.47458868041210741860e-01;
                        m_Nodes[23] = 9.64031328593135198779e-01;
                        m_Nodes[24] = 9.77515735503989208859e-01;
                        m_Nodes[25] = 9.87868941198889198522e-01;
                        m_Nodes[26] = 9.95057977847411875043e-01;
                        m_Nodes[27] = 9.99061419564818541479e-01;
                        m_Weights[0] = 5.66029764445604254401e-02;
                        m_Weights[1] = 5.65123182497720014007e-02;
                        m_Weights[2] = 5.62406340710843680283e-02;
                        m_Weights[3] = 5.57887941952840871029e-02;
                        m_Weights[4] = 5.51582460025086875967e-02;
                        m_Weights[5] = 5.43510093299111020703e-02;
                        m_Weights[6] = 5.33696700016054727236e-02;
                        m_Weights[7] = 5.22173715456320845644e-02;
                        m_Weights[8] = 5.08978051244939792244e-02;
                        m_Weights[9] = 4.94151977115517394791e-02;
                        m_Weights[10] = 4.77742985512006954609e-02;
                        m_Weights[11] = 4.59803639462838379433e-02;
                        m_Weights[12] = 4.40391404216065877684e-02;
                        m_Weights[13] = 4.19568463177187917625e-02;
                        m_Weights[14] = 3.97401518743370627695e-02;
                        m_Weights[15] = 3.73961578679665561827e-02;
                        m_Weights[16] = 3.49323728735906396512e-02;
                        m_Weights[17] = 3.23566892261737468587e-02;
                        m_Weights[18] = 2.96773577651920758008e-02;
                        m_Weights[19] = 2.69029614563915779383e-02;
                        m_Weights[20] = 2.40423880097862987501e-02;
                        m_Weights[21] = 2.11048016678193363542e-02;
                        m_Weights[22] = 1.80996145205638159777e-02;
                        m_Weights[23] = 1.50364583336639681984e-02;
                        m_Weights[24] = 1.19251607197528994055e-02;
                        m_Weights[25] = 8.77574610698129293288e-03;
                        m_Weights[26] = 5.59863226650010313854e-03;
                        m_Weights[27] = 2.40832361997383766105e-03;
                        break;

                    case 56:
                        m_Nodes[0] = 2.77970352872754370941e-02;
                        m_Nodes[1] = 8.33051868224353744403e-02;
                        m_Nodes[2] = 1.38555846810376242013e-01;
                        m_Nodes[3] = 1.93378238635275258240e-01;
                        m_Nodes[4] = 2.47602909434337203973e-01;
                        m_Nodes[5] = 3.01062253867220669053e-01;
                        m_Nodes[6] = 3.53591032174954520970e-01;
                        m_Nodes[7] = 4.05026880927091278119e-01;
                        m_Nodes[8] = 4.55210814878459578949e-01;
                        m_Nodes[9] = 5.03987718384381714195e-01;
                        m_Nodes[10] = 5.51206824855534618754e-01;
                        m_Nodes[11] = 5.96722182770663320104e-01;
                        m_Nodes[12] = 6.40393106807006894268e-01;
                        m_Nodes[13] = 6.82084612694470455502e-01;
                        m_Nodes[14] = 7.21667834450188083523e-01;
                        m_Nodes[15] = 7.59020422705128902202e-01;
                        m_Nodes[16] = 7.94026922893866498030e-01;
                        m_Nodes[17] = 8.26579132142881651672e-01;
                        m_Nodes[18] = 8.56576433762748635403e-01;
                        m_Nodes[19] = 8.83926108327827540789e-01;
                        m_Nodes[20] = 9.08543620420655490846e-01;
                        m_Nodes[21] = 9.30352880247496300547e-01;
                        m_Nodes[22] = 9.49286479561962635647e-01;
                        m_Nodes[23] = 9.65285901905490183626e-01;
                        m_Nodes[24] = 9.78301709140256383377e-01;
                        m_Nodes[25] = 9.88293715540161511090e-01;
                        m_Nodes[26] = 9.95231226081069747216e-01;
                        m_Nodes[27] = 9.99094343801465584353e-01;
                        m_Weights[0] = 5.55797463065143958463e-02;
                        m_Weights[1] = 5.54079525032451232178e-02;
                        m_Weights[2] = 5.50648959017624257963e-02;
                        m_Weights[3] = 5.45516368708894210618e-02;
                        m_Weights[4] = 5.38697618657144857090e-02;
                        m_Weights[5] = 5.30213785240107639680e-02;
                        m_Weights[6] = 5.20091091517413998431e-02;
                        m_Weights[7] = 5.08360826177984805599e-02;
                        m_Weights[8] = 4.95059246830475789198e-02;
                        m_Weights[9] = 4.80227467936002581209e-02;
                        m_Weights[10] = 4.63911333730018969092e-02;
                        m_Weights[11] = 4.46161276526922834544e-02;
                        m_Weights[12] = 4.27032160846670917147e-02;
                        m_Weights[13] = 4.06583113847445581723e-02;
                        m_Weights[14] = 3.84877342592475693483e-02;
                        m_Weights[15] = 3.61981938723151070304e-02;
                        m_Weights[16] = 3.37967671156132636088e-02;
                        m_Weights[17] = 3.12908767472932681686e-02;
                        m_Weights[18] = 2.86882684738345040396e-02;
                        m_Weights[19] = 2.59969870581615490006e-02;
                        m_Weights[20] = 2.32253515624930098133e-02;
                        m_Weights[21] = 2.03819298823885634570e-02;
                        m_Weights[22] = 1.74755129108389190902e-02;
                        m_Weights[23] = 1.45150892783214887700e-02;
                        m_Weights[24] = 1.15098243399773246657e-02;
                        m_Weights[25] = 8.46906316343052210625e-03;
                        m_Weights[26] = 5.40252224618215717641e-03;
                        m_Weights[27] = 2.32385537577709808707e-03;
                        break;

                    case 57:
                        m_Nodes[0] = 0.00000000000000000000e+00;
                        m_Nodes[1] = 5.46071510016468242198e-02;
                        m_Nodes[2] = 1.09051332808787800979e-01;
                        m_Nodes[3] = 1.63170062591264251043e-01;
                        m_Nodes[4] = 2.16801828796124036414e-01;
                        m_Nodes[5] = 2.69786573161838765763e-01;
                        m_Nodes[6] = 3.21966168395378640590e-01;
                        m_Nodes[7] = 3.73184890086594458552e-01;
                        m_Nodes[8] = 4.23289881451563950960e-01;
                        m_Nodes[9] = 4.72131609517975709588e-01;
                        m_Nodes[10] = 5.19564311391187606315e-01;
                        m_Nodes[11] = 5.65446429269236759019e-01;
                        m_Nodes[12] = 6.09641032908715365424e-01;
                        m_Nodes[13] = 6.52016228280976891249e-01;
                        m_Nodes[14] = 6.92445551199517739040e-01;
                        m_Nodes[15] = 7.30808344744523322827e-01;
                        m_Nodes[16] = 7.66990119359450195492e-01;
                        m_Nodes[17] = 8.00882894547218242076e-01;
                        m_Nodes[18] = 8.32385521150439120829e-01;
                        m_Nodes[19] = 8.61403983262046944722e-01;
                        m_Nodes[20] = 8.87851678882221329513e-01;
                        m_Nodes[21] = 9.11649678521391212728e-01;
                        m_Nodes[22] = 9.32726961067101696101e-01;
                        m_Nodes[23] = 9.51020626447876741912e-01;
                        m_Nodes[24] = 9.66476085171886679115e-01;
                        m_Nodes[25] = 9.79047226709468713799e-01;
                        m_Nodes[26] = 9.88696577650222048850e-01;
                        m_Nodes[27] = 9.95395523678430311135e-01;
                        m_Nodes[28] = 9.99125565625262850572e-01;
                        m_Weights[0] = 5.46343287565840240628e-02;
                        m_Weights[1] = 5.45528036047618864801e-02;
                        m_Weights[2] = 5.43084714524986431387e-02;
                        m_Weights[3] = 5.39020614832985746428e-02;
                        m_Weights[4] = 5.33347865848191584266e-02;
                        m_Weights[5] = 5.26083397291774324402e-02;
                        m_Weights[6] = 5.17248889205178247206e-02;
                        m_Weights[7] = 5.06870707249274086566e-02;
                        m_Weights[8] = 4.94979824020196789947e-02;
                        m_Weights[9] = 4.81611726616877512627e-02;
                        m_Weights[10] = 4.66806310736415038748e-02;
                        m_Weights[11] = 4.50607761613811580396e-02;
                        m_Weights[12] = 4.33064422162151951437e-02;
                        m_Weights[13] = 4.14228648708011258959e-02;
                        m_Weights[14] = 3.94156654754800342547e-02;
                        m_Weights[15] = 3.72908343244183916953e-02;
                        m_Weights[16] = 3.50547127823078376442e-02;
                        m_Weights[17] = 3.27139743663816265624e-02;
                        m_Weights[18] = 3.02756048426975466749e-02;
                        m_Weights[19] = 2.77468814021770372706e-02;
                        m_Weights[20] = 2.51353509906072527822e-02;
                        m_Weights[21] = 2.24488078901399589033e-02;
                        m_Weights[22] = 1.96952707013096095336e-02;
                        m_Weights[23] = 1.68829590238662969854e-02;
                        m_Weights[24] = 1.40202707916951086270e-02;
                        m_Weights[25] = 1.11157637285989974894e-02;
                        m_Weights[26] = 8.17816006624464887120e-03;
                        m_Weights[27] = 5.21653347442161669185e-03;
                        m_Weights[28] = 2.24375387225370085542e-03;
                        break;

                    case 58:
                        m_Nodes[0] = 2.68470123659423558033e-02;
                        m_Nodes[1] = 8.04636302141427293098e-02;
                        m_Nodes[2] = 1.33848250595466857022e-01;
                        m_Nodes[3] = 1.86846951835761321374e-01;
                        m_Nodes[4] = 2.39306924966153454429e-01;
                        m_Nodes[5] = 2.91076914311109189533e-01;
                        m_Nodes[6] = 3.42007653597995261248e-01;
                        m_Nodes[7] = 3.91952296330753150371e-01;
                        m_Nodes[8] = 4.40766839186839565194e-01;
                        m_Nodes[9] = 4.88310537216718463616e-01;
                        m_Nodes[10] = 5.34446309648847586400e-01;
                        m_Nodes[11] = 5.79041135130225030490e-01;
                        m_Nodes[12] = 6.21966435263079111034e-01;
                        m_Nodes[13] = 6.63098445332125266433e-01;
                        m_Nodes[14] = 7.02318571153908113480e-01;
                        m_Nodes[15] = 7.39513731020042267847e-01;
                        m_Nodes[16] = 7.74576681749652745266e-01;
                        m_Nodes[17] = 8.07406327913088141049e-01;
                        m_Nodes[18] = 8.37908013339373316352e-01;
                        m_Nodes[19] = 8.65993794074807479275e-01;
                        m_Nodes[20] = 8.91582692022030176400e-01;
                        m_Nodes[21] = 9.14600928564352540687e-01;
                        m_Nodes[22] = 9.34982137588259348481e-01;
                        m_Nodes[23] = 9.52667557518869091443e-01;
                        m_Nodes[24] = 9.67606202502924090153e-01;
                        m_Nodes[25] = 9.79755014694350309108e-01;
                        m_Nodes[26] = 9.89079008248442636500e-01;
                        m_Nodes[27] = 9.95551476597290902603e-01;
                        m_Nodes[28] = 9.99155200407386606443e-01;
                        m_Weights[0] = 5.36811198633348488639e-02;
                        m_Weights[1] = 5.35263433040582521006e-02;
                        m_Weights[2] = 5.32172364465790141035e-02;
                        m_Weights[3] = 5.27546905263708334296e-02;
                        m_Weights[4] = 5.21400391836698189713e-02;
                        m_Weights[5] = 5.13750546182857254745e-02;
                        m_Weights[6] = 5.04619424799531252977e-02;
                        m_Weights[7] = 4.94033355089623928665e-02;
                        m_Weights[8] = 4.82022859454177484068e-02;
                        m_Weights[9] = 4.68622567290263469178e-02;
                        m_Weights[10] = 4.53871115148198024782e-02;
                        m_Weights[11] = 4.37811035336402506103e-02;
                        m_Weights[12] = 4.20488633295821279094e-02;
                        m_Weights[13] = 4.01953854098678140929e-02;
                        m_Weights[14] = 3.82260138458589344018e-02;
                        m_Weights[15] = 3.61464268670850578684e-02;
                        m_Weights[16] = 3.39626204934170449282e-02;
                        m_Weights[17] = 3.16808912538765053905e-02;
                        m_Weights[18] = 2.93078180443092722278e-02;
                        m_Weights[19] = 2.68502431820246812138e-02;
                        m_Weights[20] = 2.43152527238323207342e-02;
                        m_Weights[21] = 2.17101561388448252463e-02;
                        m_Weights[22] = 1.90424654626858901973e-02;
                        m_Weights[23] = 1.63198742329197681814e-02;
                        m_Weights[24] = 1.35502371008215554008e-02;
                        m_Weights[25] = 1.07415535364291349797e-02;
                        m_Weights[26] = 7.90197385054208242183e-03;
                        m_Weights[27] = 5.03998161307001535966e-03;
                        m_Weights[28] = 2.16772324960401910983e-03;
                        break;

                    case 59:
                        m_Nodes[0] = 0.00000000000000000000e+00;
                        m_Nodes[1] = 5.27734840883100039517e-02;
                        m_Nodes[2] = 1.05399879016344143837e-01;
                        m_Nodes[3] = 1.57732505587857968115e-01;
                        m_Nodes[4] = 2.09625503392036544923e-01;
                        m_Nodes[5] = 2.60934237342811711611e-01;
                        m_Nodes[6] = 3.11515700803013700318e-01;
                        m_Nodes[7] = 3.61228914169794809992e-01;
                        m_Nodes[8] = 4.09935317810418966723e-01;
                        m_Nodes[9] = 4.57499158253266690226e-01;
                        m_Nodes[10] = 5.03787866557717978768e-01;
                        m_Nodes[11] = 5.48672427808396384372e-01;
                        m_Nodes[12] = 5.92027740704030144464e-01;
                        m_Nodes[13] = 6.33732966238850097513e-01;
                        m_Nodes[14] = 6.73671864504937227024e-01;
                        m_Nodes[15] = 7.11733118677197731596e-01;
                        m_Nodes[16] = 7.47810645278640231889e-01;
                        m_Nodes[17] = 7.81803889862360905633e-01;
                        m_Nodes[18] = 8.13618107288211571435e-01;
                        m_Nodes[19] = 8.43164625816872201470e-01;
                        m_Nodes[20] = 8.70361094292882260963e-01;
                        m_Nodes[21] = 8.95131711743472085364e-01;
                        m_Nodes[22] = 9.17407438788155281353e-01;
                        m_Nodes[23] = 9.37126190353453859405e-01;
                        m_Nodes[24] = 9.54233009376951055863e-01;
                        m_Nodes[25] = 9.68680221681781531354e-01;
                        m_Nodes[26] = 9.80427573956715688450e-01;
                        m_Nodes[27] = 9.89442365133730931782e-01;
                        m_Nodes[28] = 9.95699640383245964687e-01;
                        m_Nodes[29] = 9.99183353909294683756e-01;
                        m_Weights[0] = 5.27980126219904214155e-02;
                        m_Weights[1] = 5.27244338591279319613e-02;
                        m_Weights[2] = 5.25039026478287390509e-02;
                        m_Weights[3] = 5.21370336483753913840e-02;
                        m_Weights[4] = 5.16248493908914821464e-02;
                        m_Weights[5] = 5.09687774253939168502e-02;
                        m_Weights[6] = 5.01706463429969028107e-02;
                        m_Weights[7] = 4.92326806793619857796e-02;
                        m_Weights[8] = 4.81574947146064403877e-02;
                        m_Weights[9] = 4.69480851869620191877e-02;
                        m_Weights[10] = 4.56078229405097699987e-02;
                        m_Weights[11] = 4.41404435302973826226e-02;
                        m_Weights[12] = 4.25500368110676350160e-02;
                        m_Weights[13] = 4.08410355386866648591e-02;
                        m_Weights[14] = 3.90182030161598562521e-02;
                        m_Weights[15] = 3.70866198188695483151e-02;
                        m_Weights[16] = 3.50516696363947953811e-02;
                        m_Weights[17] = 3.29190242710190865853e-02;
                        m_Weights[18] = 3.06946278360602616231e-02;
                        m_Weights[19] = 2.83846802007990724754e-02;
                        m_Weights[20] = 2.59956197312110924992e-02;
                        m_Weights[21] = 2.35341053975598488312e-02;
                        m_Weights[22] = 2.10069982849595579577e-02;
                        m_Weights[23] = 1.84213427515433621229e-02;
                        m_Weights[24] = 1.57843472958292560893e-02;
                        m_Weights[25] = 1.31033662983893471438e-02;
                        m_Weights[26] = 1.03858855044647770569e-02;
                        m_Weights[27] = 7.63952946559766230717e-03;
                        m_Weights[28] = 4.87223916674961755492e-03;
                        m_Weights[29] = 2.09549228518365580567e-03;
                        break;

                    case 60:
                        m_Nodes[0] = 2.59597723012477985892e-02;
                        m_Nodes[1] = 7.78093339495365694193e-02;
                        m_Nodes[2] = 1.29449135396945003146e-01;
                        m_Nodes[3] = 1.80739964873425417241e-01;
                        m_Nodes[4] = 2.31543551376029338010e-01;
                        m_Nodes[5] = 2.81722937423261691691e-01;
                        m_Nodes[6] = 3.31142848268448194252e-01;
                        m_Nodes[7] = 3.79670056576797977155e-01;
                        m_Nodes[8] = 4.27173741583078389307e-01;
                        m_Nodes[9] = 4.73525841761707111108e-01;
                        m_Nodes[10] = 5.18601400058569747418e-01;
                        m_Nodes[11] = 5.62278900753944539178e-01;
                        m_Nodes[12] = 6.04440597048510363444e-01;
                        m_Nodes[13] = 6.44972828489477067813e-01;
                        m_Nodes[14] = 6.83766327381355437223e-01;
                        m_Nodes[15] = 7.20716513355730399436e-01;
                        m_Nodes[16] = 7.55723775306585686869e-01;
                        m_Nodes[17] = 7.88693739932264054570e-01;
                        m_Nodes[18] = 8.19537526162145759369e-01;
                        m_Nodes[19] = 8.48171984785929632491e-01;
                        m_Nodes[20] = 8.74519922646898315129e-01;
                        m_Nodes[21] = 8.98510310810045941938e-01;
                        m_Nodes[22] = 9.20078476177627552857e-01;
                        m_Nodes[23] = 9.39166276116423249495e-01;
                        m_Nodes[24] = 9.55722255839996107397e-01;
                        m_Nodes[25] = 9.69701788765052733722e-01;
                        m_Nodes[26] = 9.81067201752598185619e-01;
                        m_Nodes[27] = 9.89787895222221717367e-01;
                        m_Nodes[28] = 9.95840525118838173877e-01;
                        m_Nodes[29] = 9.99210123227436022034e-01;
                        m_Weights[0] = 5.19078776312206397329e-02;
                        m_Weights[1] = 5.17679431749101875438e-02;
                        m_Weights[2] = 5.14884515009809339950e-02;
                        m_Weights[3] = 5.10701560698556274045e-02;
                        m_Weights[4] = 5.05141845325093745982e-02;
                        m_Weights[5] = 4.98220356905501810112e-02;
                        m_Weights[6] = 4.89955754557568353895e-02;
                        m_Weights[7] = 4.80370318199711809639e-02;
                        m_Weights[8] = 4.69489888489122048453e-02;
                        m_Weights[9] = 4.57343797161144866611e-02;
                        m_Weights[10] = 4.43964787957871135521e-02;
                        m_Weights[11] = 4.29388928359356407911e-02;
                        m_Weights[12] = 4.13655512355847656599e-02;
                        m_Weights[13] = 3.96806954523807840600e-02;
                        m_Weights[14] = 3.78888675692432494529e-02;
                        m_Weights[15] = 3.59948980510897090503e-02;
                        m_Weights[16] = 3.40038927249354348002e-02;
                        m_Weights[17] = 3.19212190192277788480e-02;
                        m_Weights[18] = 2.97524915007394421873e-02;
                        m_Weights[19] = 2.75035567505178887851e-02;
                        m_Weights[20] = 2.51804776256335791938e-02;
                        m_Weights[21] = 2.27895169511396608670e-02;
                        m_Weights[22] = 2.03371207111078837563e-02;
                        m_Weights[23] = 1.78299010325763888436e-02;
                        m_Weights[24] = 1.52746186135325220373e-02;
                        m_Weights[25] = 1.26781664598714390463e-02;
                        m_Weights[26] = 1.00475571494532367681e-02;
                        m_Weights[27] = 7.38993116228854420668e-03;
                        m_Weights[28] = 4.71272992143352500090e-03;
                        m_Weights[29] = 2.02681196863778363879e-03;
                        break;

                    default:
                        throw new NotSupportedException("No Quadrature with " + NoOfNodes + " Nodes available.");
                }
            }
        }

        #endregion
        */


        /// <summary>
        /// transforms some vertices from local edge coordinates to the local coorinate system of this
        /// simplex
        /// </summary>
        /// <param name="EdgeIndex">specifyes the edge; <see cref="Edge"/></param>
        /// <param name="EdgeVertices">
        /// Input; Vertices in the local coordinate system of the edge;
        /// 1st index: vertex index; 2nd index: spatial coordinate index, the only valid index is 0;
        /// </param>
        /// <param name="VolumeVertices">
        /// On exit, the <paramref name="EdgeVertices"/> transformed to the local coordinate system of this simplex;
        /// </param>
        public override void TransformFaceCoordinates(int EdgeIndex, MultidimensionalArray EdgeVertices, MultidimensionalArray VolumeVertices) {
            if (EdgeVertices.Dimension != 2)
                throw new ArgumentOutOfRangeException("dimension of EdgeVertices must be 2.");
            if (VolumeVertices.Dimension != 2)
                throw new ArgumentOutOfRangeException("dimension of VolumeVertices must be 2.");
            if (EdgeVertices.GetLength(0) != VolumeVertices.GetLength(0))
                throw new ArgumentException("mismatch between number of untransfomed edges and number of transfomred edges.");
            if (EdgeVertices.GetLength(1) != 1)
                throw new ArgumentException("wrong spatial dimension for EdgeVertices-argument.");
            if (VolumeVertices.GetLength(1) != 1)
                throw new ArgumentException("wrong spatial dimension for VolumeVertices-argument.");

            int N = VolumeVertices.GetLength(0);



            if (EdgeIndex == (int)Edge.Left) {
                for (int n = 0; n < N; n++)
                    VolumeVertices[n, 0] = EdgeVertices[n, 0] - 1.0;
            } else if (EdgeIndex == (int)Edge.Right) {
                for (int n = 0; n < N; n++)
                    VolumeVertices[n, 0] = EdgeVertices[n, 0] + 1.0;
            } else {
                throw new ArgumentOutOfRangeException("illegal edge index.");
            }


        }

        /// <summary>
        /// splits the line element int two sub-lines
        /// </summary>
        /// <returns></returns>
        public override AffineTrafo[] GetSubdivision() {
            AffineTrafo[] ret = new AffineTrafo[] { AffineTrafo.Identity(1), AffineTrafo.Identity(1) };
            ret[0].Matrix.Scale(0.5);
            ret[0].Affine[0] = -0.5;

            ret[1].Matrix.Scale(0.5);
            ret[1].Affine[0] = 0.5;

            return ret;
        }

        /// <summary>
        /// Tests whether <paramref name="pt"/> is within the convex hull of
        /// vertices or not;
        /// </summary>
        /// <param name="pt"></param>
        /// <param name="tolerance"></param>
        /// <returns></returns>
        public override bool IsWithin(double[] pt, double tolerance) {
            if (pt.Length != 1)
                throw new ArgumentException("wrong spatial dimension.", "pt");

            if (pt[0] < -1.0 - tolerance || pt[0] > 1.0 + tolerance)
                return false;
            else
                return true;
        }

        /// <summary>
        /// see <see cref="RefElement.GetNodeSet(int,out MultidimensionalArray,out int[])"/>
        /// </summary>
        protected override void GetNodeSet(int px, out MultidimensionalArray Nodes, out int[] Type) {
            if (px < 2)
                throw new ArgumentOutOfRangeException("at least two nodes in each direction are required.");

            Nodes = MultidimensionalArray.Create(px, 1);
            Type = new int[Nodes.GetLength(0)];

            var Nodes1D = ilPSP.Utils.GenericBlas.Linspace(-1, 1, px);
            int cnt = 0;
            for (int i = 0; i < px; i++) {
                int i_edge = (i == 0 || i == px - 1) ? 1 : 0;
                Nodes[cnt, 0] = Nodes1D[i];
                Type[cnt] = i_edge;
                cnt++;
            }
        }

        /// <summary>
        /// <see cref="RefElement.GetInterpolationNodes_NonLin"/>
        /// </summary>
        protected override void GetInterpolationNodes_NonLin(CellType Type, out NodeSet InterpolationNodes, out PolynomialList InterpolationPolynomials, out int[] NodeType, out int[] EntityIndex) {
            switch (Type) {
                case CellType.Line_2: {
                        base.SelectNodalPolynomials(2, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                        return;
                    }
                case CellType.Line_3: {
                        base.SelectNodalPolynomials(3, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                        return;
                    }
                case CellType.Line_4: {
                        base.SelectNodalPolynomials(4, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                        return;
                    }
                case CellType.Line_5: {
                        base.SelectNodalPolynomials(5, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                        return;
                    }
                case CellType.Line_6: {
                        base.SelectNodalPolynomials(6, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                        return;
                    }
                default:
                    throw new NotImplementedException();
            }
        }

        /// <summary>
        /// <see cref="RefElement.GetForeignElementType"/>
        /// </summary>
        /// <param name="Type"></param>
        /// <param name="conv"></param>
        /// <param name="ForeignName"></param>
        /// <param name="ForeignTypeConstant"></param>
        public override void GetForeignElementType(CellType Type, RefElement.ExchangeFormats conv, out string ForeignName, out int ForeignTypeConstant) {
            ForeignName = "Bar";
            ForeignTypeConstant = 0;
            if (conv == ExchangeFormats.Gmsh) {
                if (Type == CellType.Line_2) {
                    ForeignTypeConstant = 1;
                } else if (Type == CellType.Line_3) {
                    ForeignTypeConstant = 8;
                } else if (Type == CellType.Line_4) {
                    ForeignTypeConstant = 26;
                } else if (Type == CellType.Line_5) {
                    ForeignTypeConstant = 27;
                } else if (Type == CellType.Line_6) {
                    ForeignTypeConstant = 28;
                } else {
                    throw new NotSupportedException("Wrong minor cell type");
                }
            } else if (conv == ExchangeFormats.CGNS) {
                if (Type == CellType.Line_2) {
                    ForeignTypeConstant = 3;
                } else if (Type == CellType.Line_3) {
                    ForeignTypeConstant = 4;
                } else {
                    throw new NotSupportedException("Wrong minor cell type");
                }
            } else if (conv == ExchangeFormats.GambitNeutral) {
                throw new NotSupportedException("Wrong minor cell type");
            }
        }


        double[] m_FaceTrafoGramianSqrt = new double[] { 1, 1 };

        /// <summary>
        /// Always 1.0 for the line-element.
        /// </summary>
        public override double[] FaceTrafoGramianSqrt {
            get {
                return m_FaceTrafoGramianSqrt;
            }
        }
    }
}
