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
using System.Diagnostics;
using BoSSS.Foundation.Quadrature;
using BoSSS.Platform;
using BoSSS.Platform.LinAlg;
using ilPSP.Utils;
using System.Collections.Generic;
using ilPSP.Tracing;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Foundation.Grid.RefElements {

    /// <summary>
    /// a simplex [-1,1]x[-1,1] for 2D Cartesian grids; 
    /// </summary>
    public class Square : RefElement {

        /// <summary>
        /// Indicates all four edges of the square
        /// </summary>
        public enum Edge {
            /// <summary>
            /// edge between this cell the neighbor cell with x - coordinates closer to negative infinity
            /// </summary>
            Left = 0,

            /// <summary>
            /// edge between this cell the neighbor cell with x - coordinates closer to positive infinity
            /// </summary>
            Right = 1,

            /// <summary>
            /// edge between this cell the neighbor cell with y - coordinates closer to positive infinity
            /// </summary>
            Top = 2,

            /// <summary>
            /// edge between this cell the neighbor cell with y - coordinates closer to negative infinity
            /// </summary>
            Bottom = 3

        }

        private static Square instance = null;
        private static readonly object padlock = new object();

        /// <summary>
        /// Access to the single, global instance.
        /// </summary>
        public static Square Instance {
            get {
                lock(padlock) {
                    if(instance == null) {
                        instance = new Square();
                    }
                    return instance;
                }
            }
        }


        /// <summary>
        /// constructs a new square simplex
        /// </summary>
        private Square() {
            using (new FuncTrace()) {
                // ===============
                // define vertices
                // ===============

                var _Vertices = new double[4, 2] { { -1, -1 }, { 1, -1 }, { -1, 1 }, { 1, 1 } };
                this.m_Vertices = new NodeSet(this, 4, 2);
                this.m_Vertices.InitializeFrom(_Vertices);
                this.m_Vertices.LockForever();

                m_NoOfFaces = 4;

                // ============
                // edge simplex
                // ============

                m_FaceRefElement = Line.Instance;


                // ===================================
                // define Quadrature nodes and weights
                // ===================================

                //m_QuadratureNodes   = new double[9][,];
                //m_QuadratureWeights = new double[9][];

                {
                    var qrTemp1D = QuadRuleResource.DecodeFromBase64(Resource.LineQuadRules_bin);
                    foreach(var q in qrTemp1D) {
                        int NN = q.Item2.GetLength(0);
                        int D = this.SpatialDimension;
                        var realQr = QuadRule.CreateEmpty(this, NN*NN, D);

                        for(int i = 0; i < NN; i++) {
                            for(int j = 0; j < NN; j++) {
                                realQr.Nodes[i * NN + j, 0] = q.Item2[j, 0];
                                realQr.Nodes[i * NN + j, 1] = q.Item2[i, 0];
                                realQr.Weights[i * NN + j] = q.Item3[i] * q.Item3[j];
                            }
                        }
                        
                        realQr.OrderOfPrecision = q.Item1;
                        realQr.Nodes.LockForever();
                        realQr.Weights.LockForever();
                        base.m_QuadRules.Add(realQr);
                    }
                }


                /*
                foreach (int NoOfNodes in Line.QuadratureRules1D.AvailableNoOfNodes) {

                    // quadrature rule
                    // ---------------


                    Line.QuadratureRules1D gq = new Line.QuadratureRules1D(NoOfNodes);


                    // setup volume integral nodes
                    // ---------------------------

                    QuadRule quadRule = QuadRule.CreateEmpty(this, NoOfNodes * NoOfNodes, 2);
                    var VolumeNodes = quadRule.Nodes;
                    var VolumeWeights = quadRule.Weights;

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


                        for (int j = 0; j < NoOfNodes; j++) {
                            int jno = j - NoOfNodes / 2;
                            double jsng = 1.0;
                            if (jno < 0)
                                jsng = -1.0;
                            if (jno < 0) {
                                jno *= -1;
                                if (NoOfNodes % 2 == 0)
                                    jno--;
                            }

                            VolumeNodes[i * NoOfNodes + j, 0] = gq.m_Nodes[ino] * isng;
                            VolumeNodes[i * NoOfNodes + j, 1] = gq.m_Nodes[jno] * jsng;
                            VolumeWeights[i * NoOfNodes + j] = gq.m_Weights[jno] * gq.m_Weights[ino];
                        }

                        ino++;
                    }

                    quadRule.Nodes.LockForever();
                    quadRule.OrderOfPrecision = gq.precision;
                    m_QuadRules.Add(quadRule);

                }
                 */

                // ==============================
                // define orthonormal polynomials
                // ==============================
#pragma warning disable 612
                #region POLY_DEF

                OrthonormalPolynomials = new Polynomial[300];
                Polynomial p;

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{B3587384-C1AF-49e4-A041-BBA0D8F3DA25}"));
                OrthonormalPolynomials[0] = p;
                p.AddCoeff(5.0000000000000000e-01, new int[] { 0, 0 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{66425A82-3BF1-4210-8CBD-53160BD60DEB}"));
                OrthonormalPolynomials[1] = p;
                p.AddCoeff(8.6602540378443865e-01, new int[] { 0, 1 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{63AA2458-2C4D-4608-A94D-4826F3377BB2}"));
                OrthonormalPolynomials[2] = p;
                p.AddCoeff(8.6602540378443865e-01, new int[] { 1, 0 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{A4FAA87C-70FF-4288-8275-4D6DEB0C3246}"));
                OrthonormalPolynomials[3] = p;
                p.AddCoeff(-5.5901699437494742e-01, new int[] { 0, 0 });
                p.AddCoeff(1.6770509831248423e+00, new int[] { 0, 2 });


                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{9C23E677-7E9E-4575-89B8-00245DD6C00B}"));
                OrthonormalPolynomials[4] = p;
                p.AddCoeff(1.5000000000000000e+00, new int[] { 1, 1 });


                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{A06D89AE-E12E-4ccb-8F8A-40125DE33D28}"));
                OrthonormalPolynomials[5] = p;
                p.AddCoeff(-5.5901699437494742e-01, new int[] { 0, 0 });
                p.AddCoeff(1.6770509831248423e+00, new int[] { 2, 0 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{929D6D36-9F76-400f-9600-FB5CDB279CE1}"));
                OrthonormalPolynomials[6] = p;
                p.AddCoeff(-1.9843134832984429e+00, new int[] { 0, 1 });
                p.AddCoeff(3.3071891388307382e+00, new int[] { 0, 3 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{DA3A05CE-648C-4fa7-A5D1-97353283C64F}"));
                OrthonormalPolynomials[7] = p;
                p.AddCoeff(-9.6824583655185422e-01, new int[] { 1, 0 });
                p.AddCoeff(2.9047375096555627e+00, new int[] { 1, 2 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{7B1626B4-AD2A-4eaa-9278-4994864783CA}"));
                OrthonormalPolynomials[8] = p;
                p.AddCoeff(-9.6824583655185422e-01, new int[] { 0, 1 });
                p.AddCoeff(2.9047375096555627e+00, new int[] { 2, 1 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{C0549E56-FC4F-4449-87EF-A75F7934FFF2}"));
                OrthonormalPolynomials[9] = p;
                p.AddCoeff(-1.9843134832984429e+00, new int[] { 1, 0 });
                p.AddCoeff(3.3071891388307382e+00, new int[] { 3, 0 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{8DA5F391-6E93-4fac-AFDC-CA6FC2D2B095}"));
                OrthonormalPolynomials[10] = p;
                p.AddCoeff(5.6250000000000000e-01, new int[] { 0, 0 });
                p.AddCoeff(-5.6250000000000000e+00, new int[] { 0, 2 });
                p.AddCoeff(6.5625000000000000e+00, new int[] { 0, 4 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{258A2C8C-6FE2-4b16-9FDD-07E078E912AF}"));
                OrthonormalPolynomials[11] = p;
                p.AddCoeff(-3.4369317712168800e+00, new int[] { 1, 1 });
                p.AddCoeff(5.7282196186948000e+00, new int[] { 1, 3 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{28D11649-C455-4e06-B661-99F99EB422FB}"));
                OrthonormalPolynomials[12] = p;
                p.AddCoeff(6.2500000000000000e-01, new int[] { 0, 0 });
                p.AddCoeff(-1.8750000000000000e+00, new int[] { 0, 2 });
                p.AddCoeff(-1.8750000000000000e+00, new int[] { 2, 0 });
                p.AddCoeff(5.6250000000000000e+00, new int[] { 2, 2 });
                p = new Polynomial(new Guid());

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{9B579400-2A8B-429b-8CFC-7E3D238C6411}"));
                OrthonormalPolynomials[13] = p;
                p.AddCoeff(-3.4369317712168800e+00, new int[] { 1, 1 });
                p.AddCoeff(5.7282196186948000e+00, new int[] { 3, 1 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{537AF7A0-5DA7-418a-8819-CA53D0720457}"));
                OrthonormalPolynomials[14] = p;
                p.AddCoeff(5.6250000000000000e-01, new int[] { 0, 0 });
                p.AddCoeff(-5.6250000000000000e+00, new int[] { 2, 0 });
                p.AddCoeff(6.5625000000000000e+00, new int[] { 4, 0 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{1794963D-F560-4e60-B1CF-4DD825736D21}"));
                OrthonormalPolynomials[15] = p;
                p.AddCoeff(3.1093357409581874e+00, new int[] { 0, 1 });
                p.AddCoeff(-1.4510233457804874e+01, new int[] { 0, 3 });
                p.AddCoeff(1.3059210112024387e+01, new int[] { 0, 5 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{125D392B-B374-4665-BAA4-606F50898396}"));
                OrthonormalPolynomials[16] = p;
                p.AddCoeff(9.7427857925749348e-01, new int[] { 1, 0 });
                p.AddCoeff(-9.7427857925749348e+00, new int[] { 1, 2 });
                p.AddCoeff(1.1366583424670757e+01, new int[] { 1, 4 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{3AB63474-7C87-41ec-B664-614DDD8B0E7B}"));
                OrthonormalPolynomials[17] = p;
                p.AddCoeff(2.2185299186623560e+00, new int[] { 0, 1 });
                p.AddCoeff(-3.6975498644372600e+00, new int[] { 0, 3 });
                p.AddCoeff(-6.6555897559870680e+00, new int[] { 2, 1 });
                p.AddCoeff(1.1092649593311780e+01, new int[] { 2, 3 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{1181A815-32B0-40d7-8279-AC86E468B9D6}"));
                OrthonormalPolynomials[18] = p;
                p.AddCoeff(2.2185299186623560e+00, new int[] { 1, 0 });
                p.AddCoeff(-6.6555897559870680e+00, new int[] { 1, 2 });
                p.AddCoeff(-3.6975498644372600e+00, new int[] { 3, 0 });
                p.AddCoeff(1.1092649593311780e+01, new int[] { 3, 2 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{9C7ABB84-D6FC-49c1-A97E-0785C1325C93}"));
                OrthonormalPolynomials[19] = p;
                p.AddCoeff(9.7427857925749348e-01, new int[] { 0, 1 });
                p.AddCoeff(-9.7427857925749348e+00, new int[] { 2, 1 });
                p.AddCoeff(1.1366583424670757e+01, new int[] { 4, 1 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{FAB608DB-44D3-43e2-9ECE-0C259BF2795A}"));
                OrthonormalPolynomials[20] = p;
                p.AddCoeff(3.1093357409581874e+00, new int[] { 1, 0 });
                p.AddCoeff(-1.4510233457804874e+01, new int[] { 3, 0 });
                p.AddCoeff(1.3059210112024387e+01, new int[] { 5, 0 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{B0396907-AF6E-43e5-AC3C-DA343F6C4E21}"));
                OrthonormalPolynomials[21] = p;
                p.AddCoeff(-5.6336738679124833e-01, new int[] { 0, 0 });
                p.AddCoeff(1.1830715122616215e+01, new int[] { 0, 2 });
                p.AddCoeff(-3.5492145367848645e+01, new int[] { 0, 4 });
                p.AddCoeff(2.6027573269755673e+01, new int[] { 0, 6 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{110E997F-A9FF-4404-A4A0-BAADE6629626}"));
                OrthonormalPolynomials[22] = p;
                p.AddCoeff(5.3855274811294019e+00, new int[] { 1, 1 });
                p.AddCoeff(-2.5132461578603875e+01, new int[] { 1, 3 });
                p.AddCoeff(2.2619215420743488e+01, new int[] { 1, 5 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{745B2E6A-5A08-45d3-8967-4B4C5719FB44}"));
                OrthonormalPolynomials[23] = p;
                p.AddCoeff(-6.2889411867181585e-01, new int[] { 0, 0 });
                p.AddCoeff(6.2889411867181585e+00, new int[] { 0, 2 });
                p.AddCoeff(1.8866823560154476e+00, new int[] { 2, 0 });
                p.AddCoeff(-7.3370980511711849e+00, new int[] { 0, 4 });
                p.AddCoeff(-1.8866823560154476e+01, new int[] { 2, 2 });
                p.AddCoeff(2.2011294153513555e+01, new int[] { 2, 4 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{D9C4256F-DAEF-4eb2-A9FC-03C7F2C42BD1}"));
                OrthonormalPolynomials[24] = p;
                p.AddCoeff(7.8750000000000000e+00, new int[] { 1, 1 });
                p.AddCoeff(-1.3125000000000000e+01, new int[] { 1, 3 });
                p.AddCoeff(-1.3125000000000000e+01, new int[] { 3, 1 });
                p.AddCoeff(2.1875000000000000e+01, new int[] { 3, 3 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{D9AC893A-24C5-4bc8-9311-4F21E8FBB72E}"));
                OrthonormalPolynomials[25] = p;
                p.AddCoeff(-6.2889411867181585e-01, new int[] { 0, 0 });
                p.AddCoeff(1.8866823560154476e+00, new int[] { 0, 2 });
                p.AddCoeff(6.2889411867181585e+00, new int[] { 2, 0 });
                p.AddCoeff(-1.8866823560154476e+01, new int[] { 2, 2 });
                p.AddCoeff(-7.3370980511711849e+00, new int[] { 4, 0 });
                p.AddCoeff(2.2011294153513555e+01, new int[] { 4, 2 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{EF2FE15B-2C8D-485c-90C9-A78B4570FE34}"));
                OrthonormalPolynomials[26] = p;
                p.AddCoeff(5.3855274811294019e+00, new int[] { 1, 1 });
                p.AddCoeff(-2.5132461578603875e+01, new int[] { 3, 1 });
                p.AddCoeff(2.2619215420743488e+01, new int[] { 5, 1 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{C81743DD-2482-45fc-A590-93905077964E}"));
                OrthonormalPolynomials[27] = p;
                p.AddCoeff(-5.6336738679124833e-01, new int[] { 0, 0 });
                p.AddCoeff(1.1830715122616215e+01, new int[] { 2, 0 });
                p.AddCoeff(-3.5492145367848645e+01, new int[] { 4, 0 });
                p.AddCoeff(2.6027573269755673e+01, new int[] { 6, 0 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{23E82537-5A7D-41a7-8D6F-1D081F7A5810}"));
                OrthonormalPolynomials[28] = p;
                p.AddCoeff(-4.2360755349143622e+00, new int[] { 0, 1 });
                p.AddCoeff(3.8124679814229260e+01, new int[] { 0, 3 });
                p.AddCoeff(-8.3874295591304372e+01, new int[] { 0, 5 });
                p.AddCoeff(5.1922182985093183e+01, new int[] { 0, 7 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{D5DC0D50-365B-456c-8F71-1097F7C337E5}"));
                OrthonormalPolynomials[29] = p;
                p.AddCoeff(-9.7578093724974972e-01, new int[] { 1, 0 });
                p.AddCoeff(2.0491399682244744e+01, new int[] { 1, 2 });
                p.AddCoeff(-6.1474199046734232e+01, new int[] { 1, 4 });
                p.AddCoeff(4.5081079300938437e+01, new int[] { 1, 6 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{246F1B57-32BD-4a0c-B1FE-541443CE9169}"));
                OrthonormalPolynomials[30] = p;
                p.AddCoeff(-3.4763430408260920e+00, new int[] { 0, 1 });
                p.AddCoeff(1.6222934190521763e+01, new int[] { 0, 3 });
                p.AddCoeff(1.0429029122478276e+01, new int[] { 2, 1 });
                p.AddCoeff(-1.4600640771469586e+01, new int[] { 0, 5 });
                p.AddCoeff(-4.8668802571565288e+01, new int[] { 2, 3 });
                p.AddCoeff(4.3801922314408759e+01, new int[] { 2, 5 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{E9D3AE0F-5975-44f0-94F8-8C1660052897}"));
                OrthonormalPolynomials[31] = p;
                p.AddCoeff(-2.2323526687107483e+00, new int[] { 1, 0 });
                p.AddCoeff(2.2323526687107483e+01, new int[] { 1, 2 });
                p.AddCoeff(3.7205877811845805e+00, new int[] { 3, 0 });
                p.AddCoeff(-2.6044114468292064e+01, new int[] { 1, 4 });
                p.AddCoeff(-3.7205877811845805e+01, new int[] { 3, 2 });
                p.AddCoeff(4.3406857447153439e+01, new int[] { 3, 4 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{0E4B27A6-EE5C-499a-BF87-67E969425B40}"));
                OrthonormalPolynomials[32] = p;
                p.AddCoeff(-2.2323526687107483e+00, new int[] { 0, 1 });
                p.AddCoeff(3.7205877811845805e+00, new int[] { 0, 3 });
                p.AddCoeff(2.2323526687107483e+01, new int[] { 2, 1 });
                p.AddCoeff(-3.7205877811845805e+01, new int[] { 2, 3 });
                p.AddCoeff(-2.6044114468292064e+01, new int[] { 4, 1 });
                p.AddCoeff(4.3406857447153439e+01, new int[] { 4, 3 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{BB3D9ADB-B5BE-42ad-A907-2F6395E1D295}"));
                OrthonormalPolynomials[33] = p;
                p.AddCoeff(-3.4763430408260920e+00, new int[] { 1, 0 });
                p.AddCoeff(1.0429029122478276e+01, new int[] { 1, 2 });
                p.AddCoeff(1.6222934190521763e+01, new int[] { 3, 0 });
                p.AddCoeff(-4.8668802571565288e+01, new int[] { 3, 2 });
                p.AddCoeff(-1.4600640771469586e+01, new int[] { 5, 0 });
                p.AddCoeff(4.3801922314408759e+01, new int[] { 5, 2 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{BF2C55B3-8660-4af6-A500-9796E789D557}"));
                OrthonormalPolynomials[34] = p;
                p.AddCoeff(-9.7578093724974972e-01, new int[] { 0, 1 });
                p.AddCoeff(2.0491399682244744e+01, new int[] { 2, 1 });
                p.AddCoeff(-6.1474199046734232e+01, new int[] { 4, 1 });
                p.AddCoeff(4.5081079300938437e+01, new int[] { 6, 1 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{8CA68E6C-CA64-4148-B321-F0A79A851DDE}"));
                OrthonormalPolynomials[35] = p;
                p.AddCoeff(-4.2360755349143622e+00, new int[] { 1, 0 });
                p.AddCoeff(3.8124679814229260e+01, new int[] { 3, 0 });
                p.AddCoeff(-8.3874295591304372e+01, new int[] { 5, 0 });
                p.AddCoeff(5.1922182985093183e+01, new int[] { 7, 0 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{C42959E4-85B4-4c51-ADF6-B815BCA95BF7}"));
                OrthonormalPolynomials[36] = p;
                p.AddCoeff(5.6370584725241453e-01, new int[] { 0, 0 });
                p.AddCoeff(-2.0293410501086923e+01, new int[] { 0, 2 });
                p.AddCoeff(1.1161375775597808e+02, new int[] { 0, 4 });
                p.AddCoeff(-1.9346384677702867e+02, new int[] { 0, 6 });
                p.AddCoeff(1.0364134648769393e+02, new int[] { 0, 8 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{D8A8543D-C267-49a8-BDC2-8102F9A75A6A}"));
                OrthonormalPolynomials[37] = p;
                p.AddCoeff(-7.3370980511711849e+00, new int[] { 1, 1 });
                p.AddCoeff(6.6033882460540664e+01, new int[] { 1, 3 });
                p.AddCoeff(-1.4527454141318946e+02, new int[] { 1, 5 });
                p.AddCoeff(8.9931858970069667e+01, new int[] { 1, 7 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{C7AC8D44-5BBF-4bd6-8DCF-03A2E47C2305}"));
                OrthonormalPolynomials[38] = p;
                p.AddCoeff(6.2986388658582419e-01, new int[] { 0, 0 });
                p.AddCoeff(-1.3227141618302308e+01, new int[] { 0, 2 });
                p.AddCoeff(-1.8895916597574726e+00, new int[] { 2, 0 });
                p.AddCoeff(3.9681424854906924e+01, new int[] { 0, 4 });
                p.AddCoeff(3.9681424854906924e+01, new int[] { 2, 2 });
                p.AddCoeff(-2.9099711560265078e+01, new int[] { 0, 6 });
                p.AddCoeff(-1.1904427456472077e+02, new int[] { 2, 4 });
                p.AddCoeff(8.7299134680795233e+01, new int[] { 2, 6 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{205B1A0C-7FB0-4e5b-8504-557DE06CA20E}"));
                OrthonormalPolynomials[39] = p;
                p.AddCoeff(-1.2339793669770172e+01, new int[] { 1, 1 });
                p.AddCoeff(5.7585703792260801e+01, new int[] { 1, 3 });
                p.AddCoeff(2.0566322782950286e+01, new int[] { 3, 1 });
                p.AddCoeff(-5.1827133413034721e+01, new int[] { 1, 5 });
                p.AddCoeff(-9.5976172987101335e+01, new int[] { 3, 3 });
                p.AddCoeff(8.6378555688391202e+01, new int[] { 3, 5 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{52C8CDA0-AC55-4dd9-9258-5DE2696E8ED2}"));
                OrthonormalPolynomials[40] = p;
                p.AddCoeff(6.3281250000000000e-01, new int[] { 0, 0 });
                p.AddCoeff(-6.3281250000000000e+00, new int[] { 0, 2 });
                p.AddCoeff(-6.3281250000000000e+00, new int[] { 2, 0 });
                p.AddCoeff(7.3828125000000000e+00, new int[] { 0, 4 });
                p.AddCoeff(6.3281250000000000e+01, new int[] { 2, 2 });
                p.AddCoeff(7.3828125000000000e+00, new int[] { 4, 0 });
                p.AddCoeff(-7.3828125000000000e+01, new int[] { 2, 4 });
                p.AddCoeff(-7.3828125000000000e+01, new int[] { 4, 2 });
                p.AddCoeff(8.6132812500000000e+01, new int[] { 4, 4 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{4AD216E8-BFEB-4cdb-8B1F-2F8C6A6C9CFF}"));
                OrthonormalPolynomials[41] = p;
                p.AddCoeff(-1.2339793669770172e+01, new int[] { 1, 1 });
                p.AddCoeff(2.0566322782950286e+01, new int[] { 1, 3 });
                p.AddCoeff(5.7585703792260801e+01, new int[] { 3, 1 });
                p.AddCoeff(-9.5976172987101335e+01, new int[] { 3, 3 });
                p.AddCoeff(-5.1827133413034721e+01, new int[] { 5, 1 });
                p.AddCoeff(8.6378555688391202e+01, new int[] { 5, 3 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{0B4280A3-6D75-4fba-B744-B5A24E4B8606}"));
                OrthonormalPolynomials[42] = p;
                p.AddCoeff(6.2986388658582419e-01, new int[] { 0, 0 });
                p.AddCoeff(-1.8895916597574726e+00, new int[] { 0, 2 });
                p.AddCoeff(-1.3227141618302308e+01, new int[] { 2, 0 });
                p.AddCoeff(3.9681424854906924e+01, new int[] { 2, 2 });
                p.AddCoeff(3.9681424854906924e+01, new int[] { 4, 0 });
                p.AddCoeff(-1.1904427456472077e+02, new int[] { 4, 2 });
                p.AddCoeff(-2.9099711560265078e+01, new int[] { 6, 0 });
                p.AddCoeff(8.7299134680795233e+01, new int[] { 6, 2 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{A67729BC-6837-40e6-AE2A-D3BA76491EFE}"));
                OrthonormalPolynomials[43] = p;
                p.AddCoeff(-7.3370980511711849e+00, new int[] { 1, 1 });
                p.AddCoeff(6.6033882460540664e+01, new int[] { 3, 1 });
                p.AddCoeff(-1.4527454141318946e+02, new int[] { 5, 1 });
                p.AddCoeff(8.9931858970069667e+01, new int[] { 7, 1 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{BEB61B32-DABD-426b-A4D9-6414F4ECD954}"));
                OrthonormalPolynomials[44] = p;
                p.AddCoeff(5.6370584725241453e-01, new int[] { 0, 0 });
                p.AddCoeff(-2.0293410501086923e+01, new int[] { 2, 0 });
                p.AddCoeff(1.1161375775597808e+02, new int[] { 4, 0 });
                p.AddCoeff(-1.9346384677702867e+02, new int[] { 6, 0 });
                p.AddCoeff(1.0364134648769393e+02, new int[] { 8, 0 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{827F5652-506B-4539-8EFA-48A5E6833BFA}"));
                OrthonormalPolynomials[45] = p;
                p.AddCoeff(5.3634889344348132e+00, new int[] { 0, 1 });
                p.AddCoeff(-7.8664504371710593e+01, new int[] { 0, 3 });
                p.AddCoeff(3.0679156704967131e+02, new int[] { 0, 5 });
                p.AddCoeff(-4.3827366721381616e+02, new int[] { 0, 7 });
                p.AddCoeff(2.0696256507319096e+02, new int[] { 0, 9 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{1004C233-B625-4f0c-91BC-C906400559B4}"));
                OrthonormalPolynomials[46] = p;
                p.AddCoeff(9.7636716796484277e-01, new int[] { 1, 0 });
                p.AddCoeff(-3.5149218046734340e+01, new int[] { 1, 2 });
                p.AddCoeff(1.9332069925703887e+02, new int[] { 1, 4 });
                p.AddCoeff(-3.3508921204553404e+02, new int[] { 1, 6 });
                p.AddCoeff(1.7951207788153609e+02, new int[] { 1, 8 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{6677CC67-50A5-4bab-A38A-0725D253BD63}"));
                OrthonormalPolynomials[47] = p;
                p.AddCoeff(4.7360764269461488e+00, new int[] { 0, 1 });
                p.AddCoeff(-4.2624687842515340e+01, new int[] { 0, 3 });
                p.AddCoeff(-1.4208229280838447e+01, new int[] { 2, 1 });
                p.AddCoeff(9.3774313253533747e+01, new int[] { 0, 5 });
                p.AddCoeff(1.2787406352754602e+02, new int[] { 2, 3 });
                p.AddCoeff(-5.8050765347425653e+01, new int[] { 0, 7 });
                p.AddCoeff(-2.8132293976060124e+02, new int[] { 2, 5 });
                p.AddCoeff(1.7415229604227696e+02, new int[] { 2, 7 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{7B298595-CC46-4564-B931-F40B30BC0B6C}"));
                OrthonormalPolynomials[48] = p;
                p.AddCoeff(2.2357950033209664e+00, new int[] { 1, 0 });
                p.AddCoeff(-4.6951695069740294e+01, new int[] { 1, 2 });
                p.AddCoeff(-3.7263250055349439e+00, new int[] { 3, 0 });
                p.AddCoeff(1.4085508520922088e+02, new int[] { 1, 4 });
                p.AddCoeff(7.8252825116233823e+01, new int[] { 3, 2 });
                p.AddCoeff(-1.0329372915342865e+02, new int[] { 1, 6 });
                p.AddCoeff(-2.3475847534870147e+02, new int[] { 3, 4 });
                p.AddCoeff(1.7215621525571441e+02, new int[] { 3, 6 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{CCCFB548-366C-46bb-A5C1-B4321FCD4BDF}"));
                OrthonormalPolynomials[49] = p;
                p.AddCoeff(3.4980027085779608e+00, new int[] { 0, 1 });
                p.AddCoeff(-1.6324012640030484e+01, new int[] { 0, 3 });
                p.AddCoeff(-3.4980027085779608e+01, new int[] { 2, 1 });
                p.AddCoeff(1.4691611376027435e+01, new int[] { 0, 5 });
                p.AddCoeff(1.6324012640030484e+02, new int[] { 2, 3 });
                p.AddCoeff(4.0810031600076209e+01, new int[] { 4, 1 });
                p.AddCoeff(-1.4691611376027435e+02, new int[] { 2, 5 });
                p.AddCoeff(-1.9044681413368898e+02, new int[] { 4, 3 });
                p.AddCoeff(1.7140213272032008e+02, new int[] { 4, 5 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{DCD6926D-3D6C-412f-8C75-28EF1A7436FC}"));
                OrthonormalPolynomials[50] = p;
                p.AddCoeff(3.4980027085779608e+00, new int[] { 1, 0 });
                p.AddCoeff(-3.4980027085779608e+01, new int[] { 1, 2 });
                p.AddCoeff(-1.6324012640030484e+01, new int[] { 3, 0 });
                p.AddCoeff(4.0810031600076209e+01, new int[] { 1, 4 });
                p.AddCoeff(1.6324012640030484e+02, new int[] { 3, 2 });
                p.AddCoeff(1.4691611376027435e+01, new int[] { 5, 0 });
                p.AddCoeff(-1.9044681413368898e+02, new int[] { 3, 4 });
                p.AddCoeff(-1.4691611376027435e+02, new int[] { 5, 2 });
                p.AddCoeff(1.7140213272032008e+02, new int[] { 5, 4 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{9C43BB79-2809-4848-B584-62B12FE0C3F4}"));
                OrthonormalPolynomials[51] = p;
                p.AddCoeff(2.2357950033209664e+00, new int[] { 0, 1 });
                p.AddCoeff(-3.7263250055349439e+00, new int[] { 0, 3 });
                p.AddCoeff(-4.6951695069740294e+01, new int[] { 2, 1 });
                p.AddCoeff(7.8252825116233823e+01, new int[] { 2, 3 });
                p.AddCoeff(1.4085508520922088e+02, new int[] { 4, 1 });
                p.AddCoeff(-2.3475847534870147e+02, new int[] { 4, 3 });
                p.AddCoeff(-1.0329372915342865e+02, new int[] { 6, 1 });
                p.AddCoeff(1.7215621525571441e+02, new int[] { 6, 3 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{C15713F7-BC51-4e35-B39D-091268C2498E}"));
                OrthonormalPolynomials[52] = p;
                p.AddCoeff(4.7360764269461488e+00, new int[] { 1, 0 });
                p.AddCoeff(-1.4208229280838447e+01, new int[] { 1, 2 });
                p.AddCoeff(-4.2624687842515340e+01, new int[] { 3, 0 });
                p.AddCoeff(1.2787406352754602e+02, new int[] { 3, 2 });
                p.AddCoeff(9.3774313253533747e+01, new int[] { 5, 0 });
                p.AddCoeff(-2.8132293976060124e+02, new int[] { 5, 2 });
                p.AddCoeff(-5.8050765347425653e+01, new int[] { 7, 0 });
                p.AddCoeff(1.7415229604227696e+02, new int[] { 7, 2 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{0A8CD397-6C9F-4298-809C-86D8E5E16D6E}"));
                OrthonormalPolynomials[53] = p;
                p.AddCoeff(9.7636716796484277e-01, new int[] { 0, 1 });
                p.AddCoeff(-3.5149218046734340e+01, new int[] { 2, 1 });
                p.AddCoeff(1.9332069925703887e+02, new int[] { 4, 1 });
                p.AddCoeff(-3.3508921204553404e+02, new int[] { 6, 1 });
                p.AddCoeff(1.7951207788153609e+02, new int[] { 8, 1 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{9FAE154C-16C5-4f92-85B4-B5E1DB79A91E}"));
                OrthonormalPolynomials[54] = p;
                p.AddCoeff(5.3634889344348132e+00, new int[] { 1, 0 });
                p.AddCoeff(-7.8664504371710593e+01, new int[] { 3, 0 });
                p.AddCoeff(3.0679156704967131e+02, new int[] { 5, 0 });
                p.AddCoeff(-4.3827366721381616e+02, new int[] { 7, 0 });
                p.AddCoeff(2.0696256507319096e+02, new int[] { 9, 0 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{499D63EF-82CA-4e66-B96C-8EBDB4AFAEE3}"));
                OrthonormalPolynomials[55] = p;
                p.AddCoeff(-5.6387161871526938e-01, new int[] { 0, 0 });
                p.AddCoeff(3.1012939029339816e+01, new int[] { 0, 2 });
                p.AddCoeff(-2.6877880492094507e+02, new int[] { 0, 4 });
                p.AddCoeff(8.0633641476283521e+02, new int[] { 0, 6 });
                p.AddCoeff(-9.7912278935487132e+02, new int[] { 0, 8 });
                p.AddCoeff(4.1340739994983456e+02, new int[] { 0, 10 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{18A10FFA-87EC-4c34-B3C4-A9CB0234A9D2}"));
                OrthonormalPolynomials[56] = p;
                p.AddCoeff(9.2898353402745553e+00, new int[] { 1, 1 });
                p.AddCoeff(-1.3625091832402681e+02, new int[] { 1, 3 });
                p.AddCoeff(5.3137858146370456e+02, new int[] { 1, 5 });
                p.AddCoeff(-7.5911225923386366e+02, new int[] { 1, 7 });
                p.AddCoeff(3.5846967797154673e+02, new int[] { 1, 9 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{ED8C6208-922E-4faf-8794-D1DDAA7AF226}"));
                OrthonormalPolynomials[57] = p;
                p.AddCoeff(-6.3024229688525597e-01, new int[] { 0, 0 });
                p.AddCoeff(2.2688722687869215e+01, new int[] { 0, 2 });
                p.AddCoeff(1.8907268906557679e+00, new int[] { 2, 0 });
                p.AddCoeff(-1.2478797478328068e+02, new int[] { 0, 4 });
                p.AddCoeff(-6.8066168063607645e+01, new int[] { 2, 2 });
                p.AddCoeff(2.1629915629101985e+02, new int[] { 0, 6 });
                p.AddCoeff(3.7436392434984205e+02, new int[] { 2, 4 });
                p.AddCoeff(-1.1587454801304635e+02, new int[] { 0, 8 });
                p.AddCoeff(-6.4889746887305955e+02, new int[] { 2, 6 });
                p.AddCoeff(3.4762364403913904e+02, new int[] { 2, 8 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{C9F82908-D159-4602-9A9D-A92869AEC04C}"));
                OrthonormalPolynomials[58] = p;
                p.AddCoeff(1.6811403600402466e+01, new int[] { 1, 1 });
                p.AddCoeff(-1.5130263240362219e+02, new int[] { 1, 3 });
                p.AddCoeff(-2.8019006000670777e+01, new int[] { 3, 1 });
                p.AddCoeff(3.3286579128796883e+02, new int[] { 1, 5 });
                p.AddCoeff(2.5217105400603699e+02, new int[] { 3, 3 });
                p.AddCoeff(-2.0605977555921880e+02, new int[] { 1, 7 });
                p.AddCoeff(-5.5477631881328138e+02, new int[] { 3, 5 });
                p.AddCoeff(3.4343295926536466e+02, new int[] { 3, 7 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{2F9A9504-268A-4701-B2B7-1520647DF397}"));
                OrthonormalPolynomials[59] = p;
                p.AddCoeff(-6.3378831014015437e-01, new int[] { 0, 0 });
                p.AddCoeff(1.3309554512943242e+01, new int[] { 0, 2 });
                p.AddCoeff(6.3378831014015437e+00, new int[] { 2, 0 });
                p.AddCoeff(-3.9928663538829725e+01, new int[] { 0, 4 });
                p.AddCoeff(-1.3309554512943242e+02, new int[] { 2, 2 });
                p.AddCoeff(-7.3941969516351343e+00, new int[] { 4, 0 });
                p.AddCoeff(2.9281019928475132e+01, new int[] { 0, 6 });
                p.AddCoeff(3.9928663538829725e+02, new int[] { 2, 4 });
                p.AddCoeff(1.5527813598433782e+02, new int[] { 4, 2 });
                p.AddCoeff(-2.9281019928475132e+02, new int[] { 2, 6 });
                p.AddCoeff(-4.6583440795301346e+02, new int[] { 4, 4 });
                p.AddCoeff(3.4161189916554320e+02, new int[] { 4, 6 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{9B62F049-7003-4f1a-9DCB-35C4B39320E1}"));
                OrthonormalPolynomials[60] = p;
                p.AddCoeff(1.9335937500000000e+01, new int[] { 1, 1 });
                p.AddCoeff(-9.0234375000000000e+01, new int[] { 1, 3 });
                p.AddCoeff(-9.0234375000000000e+01, new int[] { 3, 1 });
                p.AddCoeff(8.1210937500000000e+01, new int[] { 1, 5 });
                p.AddCoeff(4.2109375000000000e+02, new int[] { 3, 3 });
                p.AddCoeff(8.1210937500000000e+01, new int[] { 5, 1 });
                p.AddCoeff(-3.7898437500000000e+02, new int[] { 3, 5 });
                p.AddCoeff(-3.7898437500000000e+02, new int[] { 5, 3 });
                p.AddCoeff(3.4108593750000000e+02, new int[] { 5, 5 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{2967DE6E-8D46-476d-96C5-77ED932DC6C4}"));
                OrthonormalPolynomials[61] = p;
                p.AddCoeff(-6.3378831014015437e-01, new int[] { 0, 0 });
                p.AddCoeff(6.3378831014015437e+00, new int[] { 0, 2 });
                p.AddCoeff(1.3309554512943242e+01, new int[] { 2, 0 });
                p.AddCoeff(-7.3941969516351343e+00, new int[] { 0, 4 });
                p.AddCoeff(-1.3309554512943242e+02, new int[] { 2, 2 });
                p.AddCoeff(-3.9928663538829725e+01, new int[] { 4, 0 });
                p.AddCoeff(1.5527813598433782e+02, new int[] { 2, 4 });
                p.AddCoeff(3.9928663538829725e+02, new int[] { 4, 2 });
                p.AddCoeff(2.9281019928475132e+01, new int[] { 6, 0 });
                p.AddCoeff(-4.6583440795301346e+02, new int[] { 4, 4 });
                p.AddCoeff(-2.9281019928475132e+02, new int[] { 6, 2 });
                p.AddCoeff(3.4161189916554320e+02, new int[] { 6, 4 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{CD04265A-5A67-4ea7-989C-C234F68D78BF}"));
                OrthonormalPolynomials[62] = p;
                p.AddCoeff(1.6811403600402466e+01, new int[] { 1, 1 });
                p.AddCoeff(-2.8019006000670777e+01, new int[] { 1, 3 });
                p.AddCoeff(-1.5130263240362219e+02, new int[] { 3, 1 });
                p.AddCoeff(2.5217105400603699e+02, new int[] { 3, 3 });
                p.AddCoeff(3.3286579128796883e+02, new int[] { 5, 1 });
                p.AddCoeff(-5.5477631881328138e+02, new int[] { 5, 3 });
                p.AddCoeff(-2.0605977555921880e+02, new int[] { 7, 1 });
                p.AddCoeff(3.4343295926536466e+02, new int[] { 7, 3 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{C75BB612-8B0D-4c9e-B6BA-AF7551745A57}"));
                OrthonormalPolynomials[63] = p;
                p.AddCoeff(-6.3024229688525597e-01, new int[] { 0, 0 });
                p.AddCoeff(1.8907268906557679e+00, new int[] { 0, 2 });
                p.AddCoeff(2.2688722687869215e+01, new int[] { 2, 0 });
                p.AddCoeff(-6.8066168063607645e+01, new int[] { 2, 2 });
                p.AddCoeff(-1.2478797478328068e+02, new int[] { 4, 0 });
                p.AddCoeff(3.7436392434984205e+02, new int[] { 4, 2 });
                p.AddCoeff(2.1629915629101985e+02, new int[] { 6, 0 });
                p.AddCoeff(-6.4889746887305955e+02, new int[] { 6, 2 });
                p.AddCoeff(-1.1587454801304635e+02, new int[] { 8, 0 });
                p.AddCoeff(3.4762364403913904e+02, new int[] { 8, 2 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{5AFB47B6-A4BE-4b65-89C9-1BBAD0B47710}"));
                OrthonormalPolynomials[64] = p;
                p.AddCoeff(9.2898353402745553e+00, new int[] { 1, 1 });
                p.AddCoeff(-1.3625091832402681e+02, new int[] { 3, 1 });
                p.AddCoeff(5.3137858146370456e+02, new int[] { 5, 1 });
                p.AddCoeff(-7.5911225923386366e+02, new int[] { 7, 1 });
                p.AddCoeff(3.5846967797154673e+02, new int[] { 9, 1 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{B327BC8F-5CB1-4943-A876-6B8D19127AF1}"));
                OrthonormalPolynomials[65] = p;
                p.AddCoeff(-5.6387161871526938e-01, new int[] { 0, 0 });
                p.AddCoeff(3.1012939029339816e+01, new int[] { 2, 0 });
                p.AddCoeff(-2.6877880492094507e+02, new int[] { 4, 0 });
                p.AddCoeff(8.0633641476283521e+02, new int[] { 6, 0 });
                p.AddCoeff(-9.7912278935487132e+02, new int[] { 8, 0 });
                p.AddCoeff(4.1340739994983456e+02, new int[] { 10, 0 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{F8B81CB7-8BE4-498C-A2D6-43BE682A7E5D}"));
                OrthonormalPolynomials[66] = p;
                p.AddCoeff(-6.4912329016713177e+00, new int[] { 0, 1 });
                p.AddCoeff(1.4064337953621188e+02, new int[] { 0, 3 });
                p.AddCoeff(-8.4386027721727130e+02, new int[] { 0, 5 });
                p.AddCoeff(2.0493749589562303e+03, new int[] { 0, 7 });
                p.AddCoeff(-2.1632291233426875e+03, new int[] { 0, 9 });
                p.AddCoeff(8.2596021073084433e+02, new int[] { 0, 11 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{4611FBAF-0B0F-4D59-B30E-30FB8DA30B83}"));
                OrthonormalPolynomials[67] = p;
                p.AddCoeff(-9.7665429256095239e-01, new int[] { 1, 0 });
                p.AddCoeff(5.3715986090852381e+01, new int[] { 1, 2 });
                p.AddCoeff(-4.6553854612072064e+02, new int[] { 1, 4 });
                p.AddCoeff(1.3966156383621619e+03, new int[] { 1, 6 });
                p.AddCoeff(-1.6958904180111966e+03, new int[] { 1, 8 });
                p.AddCoeff(7.1604262093806079e+02, new int[] { 1, 10 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{1B4F2501-898F-4C50-B9D5-B520150456D9}"));
                OrthonormalPolynomials[68] = p;
                p.AddCoeff(-5.9965629269820774e+00, new int[] { 0, 1 });
                p.AddCoeff(8.7949589595737135e+01, new int[] { 0, 3 });
                p.AddCoeff(1.7989688780946232e+01, new int[] { 2, 1 });
                p.AddCoeff(-3.4300339942337483e+02, new int[] { 0, 5 });
                p.AddCoeff(-2.6384876878721141e+02, new int[] { 2, 3 });
                p.AddCoeff(4.9000485631910690e+02, new int[] { 0, 7 });
                p.AddCoeff(1.0290101982701245e+03, new int[] { 2, 5 });
                p.AddCoeff(-2.3139118215068937e+02, new int[] { 0, 9 });
                p.AddCoeff(-1.4700145689573207e+03, new int[] { 2, 7 });
                p.AddCoeff(6.9417354645206810e+02, new int[] { 2, 9 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{8105F7BB-48D6-49D8-B1B9-DFD1564B77DA}"));
                OrthonormalPolynomials[69] = p;
                p.AddCoeff(-2.2371382266342774e+00, new int[] { 1, 0 });
                p.AddCoeff(8.0536976158833985e+01, new int[] { 1, 2 });
                p.AddCoeff(3.7285637110571289e+00, new int[] { 3, 0 });
                p.AddCoeff(-4.4295336887358692e+02, new int[] { 1, 4 });
                p.AddCoeff(-1.3422829359805664e+02, new int[] { 3, 2 });
                p.AddCoeff(7.6778583938088399e+02, new int[] { 1, 6 });
                p.AddCoeff(7.3825561478931153e+02, new int[] { 3, 4 });
                p.AddCoeff(-4.1131384252547357e+02, new int[] { 1, 8 });
                p.AddCoeff(-1.2796430656348067e+03, new int[] { 3, 6 });
                p.AddCoeff(6.8552307087578928e+02, new int[] { 3, 8 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{4070443C-1815-48D2-B3A4-DE52A7F3C548}"));
                OrthonormalPolynomials[70] = p;
                p.AddCoeff(-4.7655849767786575e+00, new int[] { 0, 1 });
                p.AddCoeff(4.2890264791007917e+01, new int[] { 0, 3 });
                p.AddCoeff(4.7655849767786575e+01, new int[] { 2, 1 });
                p.AddCoeff(-9.4358582540217418e+01, new int[] { 0, 5 });
                p.AddCoeff(-4.2890264791007917e+02, new int[] { 2, 3 });
                p.AddCoeff(-5.5598491395751004e+01, new int[] { 4, 1 });
                p.AddCoeff(5.8412455858229830e+01, new int[] { 0, 7 });
                p.AddCoeff(9.4358582540217418e+02, new int[] { 2, 5 });
                p.AddCoeff(5.0038642256175904e+02, new int[] { 4, 3 });
                p.AddCoeff(-5.8412455858229830e+02, new int[] { 2, 7 });
                p.AddCoeff(-1.1008501296358699e+03, new int[] { 4, 5 });
                p.AddCoeff(6.8147865167934802e+02, new int[] { 4, 7 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{27292CDF-D254-49EC-9804-90898A4966A5}"));
                OrthonormalPolynomials[71] = p;
                p.AddCoeff(-3.5033967020804877e+00, new int[] { 1, 0 });
                p.AddCoeff(7.3571330743690242e+01, new int[] { 1, 2 });
                p.AddCoeff(1.6349184609708943e+01, new int[] { 3, 0 });
                p.AddCoeff(-2.2071399223107073e+02, new int[] { 1, 4 });
                p.AddCoeff(-3.4333287680388779e+02, new int[] { 3, 2 });
                p.AddCoeff(-1.4714266148738048e+01, new int[] { 5, 0 });
                p.AddCoeff(1.6185692763611853e+02, new int[] { 1, 6 });
                p.AddCoeff(1.0299986304116634e+03, new int[] { 3, 4 });
                p.AddCoeff(3.0899958912349902e+02, new int[] { 5, 2 });
                p.AddCoeff(-7.5533232896855315e+02, new int[] { 3, 6 });
                p.AddCoeff(-9.2699876737049705e+02, new int[] { 5, 4 });
                p.AddCoeff(6.7979909607169783e+02, new int[] { 5, 6 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{2501D2C1-AF74-4EEB-9DA6-76865AC1CA72}"));
                OrthonormalPolynomials[72] = p;
                p.AddCoeff(-3.5033967020804877e+00, new int[] { 0, 1 });
                p.AddCoeff(1.6349184609708943e+01, new int[] { 0, 3 });
                p.AddCoeff(7.3571330743690242e+01, new int[] { 2, 1 });
                p.AddCoeff(-1.4714266148738048e+01, new int[] { 0, 5 });
                p.AddCoeff(-3.4333287680388779e+02, new int[] { 2, 3 });
                p.AddCoeff(-2.2071399223107073e+02, new int[] { 4, 1 });
                p.AddCoeff(3.0899958912349902e+02, new int[] { 2, 5 });
                p.AddCoeff(1.0299986304116634e+03, new int[] { 4, 3 });
                p.AddCoeff(1.6185692763611853e+02, new int[] { 6, 1 });
                p.AddCoeff(-9.2699876737049705e+02, new int[] { 4, 5 });
                p.AddCoeff(-7.5533232896855315e+02, new int[] { 6, 3 });
                p.AddCoeff(6.7979909607169783e+02, new int[] { 6, 5 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{92726165-A4CF-4E7A-9BC0-E7C49FFAFD21}"));
                OrthonormalPolynomials[73] = p;
                p.AddCoeff(-4.7655849767786575e+00, new int[] { 1, 0 });
                p.AddCoeff(4.7655849767786575e+01, new int[] { 1, 2 });
                p.AddCoeff(4.2890264791007917e+01, new int[] { 3, 0 });
                p.AddCoeff(-5.5598491395751004e+01, new int[] { 1, 4 });
                p.AddCoeff(-4.2890264791007917e+02, new int[] { 3, 2 });
                p.AddCoeff(-9.4358582540217418e+01, new int[] { 5, 0 });
                p.AddCoeff(5.0038642256175904e+02, new int[] { 3, 4 });
                p.AddCoeff(9.4358582540217418e+02, new int[] { 5, 2 });
                p.AddCoeff(5.8412455858229830e+01, new int[] { 7, 0 });
                p.AddCoeff(-1.1008501296358699e+03, new int[] { 5, 4 });
                p.AddCoeff(-5.8412455858229830e+02, new int[] { 7, 2 });
                p.AddCoeff(6.8147865167934802e+02, new int[] { 7, 4 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{5DD5FCEA-50EA-4677-BA10-058EDDE836CD}"));
                OrthonormalPolynomials[74] = p;
                p.AddCoeff(-2.2371382266342774e+00, new int[] { 0, 1 });
                p.AddCoeff(3.7285637110571289e+00, new int[] { 0, 3 });
                p.AddCoeff(8.0536976158833985e+01, new int[] { 2, 1 });
                p.AddCoeff(-1.3422829359805664e+02, new int[] { 2, 3 });
                p.AddCoeff(-4.4295336887358692e+02, new int[] { 4, 1 });
                p.AddCoeff(7.3825561478931153e+02, new int[] { 4, 3 });
                p.AddCoeff(7.6778583938088399e+02, new int[] { 6, 1 });
                p.AddCoeff(-1.2796430656348067e+03, new int[] { 6, 3 });
                p.AddCoeff(-4.1131384252547357e+02, new int[] { 8, 1 });
                p.AddCoeff(6.8552307087578928e+02, new int[] { 8, 3 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{9FDA31F9-F928-41F6-AD48-207523903B0F}"));
                OrthonormalPolynomials[75] = p;
                p.AddCoeff(-5.9965629269820774e+00, new int[] { 1, 0 });
                p.AddCoeff(1.7989688780946232e+01, new int[] { 1, 2 });
                p.AddCoeff(8.7949589595737135e+01, new int[] { 3, 0 });
                p.AddCoeff(-2.6384876878721141e+02, new int[] { 3, 2 });
                p.AddCoeff(-3.4300339942337483e+02, new int[] { 5, 0 });
                p.AddCoeff(1.0290101982701245e+03, new int[] { 5, 2 });
                p.AddCoeff(4.9000485631910690e+02, new int[] { 7, 0 });
                p.AddCoeff(-1.4700145689573207e+03, new int[] { 7, 2 });
                p.AddCoeff(-2.3139118215068937e+02, new int[] { 9, 0 });
                p.AddCoeff(6.9417354645206810e+02, new int[] { 9, 2 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{86AF4D2D-5FBD-4DD2-8B2B-598EC00DF57E}"));
                OrthonormalPolynomials[76] = p;
                p.AddCoeff(-9.7665429256095239e-01, new int[] { 0, 1 });
                p.AddCoeff(5.3715986090852381e+01, new int[] { 2, 1 });
                p.AddCoeff(-4.6553854612072064e+02, new int[] { 4, 1 });
                p.AddCoeff(1.3966156383621619e+03, new int[] { 6, 1 });
                p.AddCoeff(-1.6958904180111966e+03, new int[] { 8, 1 });
                p.AddCoeff(7.1604262093806079e+02, new int[] { 10, 1 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{5962B38B-F8FC-4FCC-9A0B-F133FC46DCEB}"));
                OrthonormalPolynomials[77] = p;
                p.AddCoeff(-6.4912329016713177e+00, new int[] { 1, 0 });
                p.AddCoeff(1.4064337953621188e+02, new int[] { 3, 0 });
                p.AddCoeff(-8.4386027721727130e+02, new int[] { 5, 0 });
                p.AddCoeff(2.0493749589562303e+03, new int[] { 7, 0 });
                p.AddCoeff(-2.1632291233426875e+03, new int[] { 9, 0 });
                p.AddCoeff(8.2596021073084433e+02, new int[] { 11, 0 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{BD4773CB-9AF7-4D38-B436-7EA2B8EF2974}"));
                OrthonormalPolynomials[78] = p;
                p.AddCoeff(5.6396484375000000e-01, new int[] { 0, 0 });
                p.AddCoeff(-4.3989257812500000e+01, new int[] { 0, 2 });
                p.AddCoeff(5.4986572265625000e+02, new int[] { 0, 4 });
                p.AddCoeff(-2.4927246093750000e+03, new int[] { 0, 6 });
                p.AddCoeff(5.0744750976562500e+03, new int[] { 0, 8 });
                p.AddCoeff(-4.7361767578125000e+03, new int[] { 0, 10 });
                p.AddCoeff(1.6504858398437500e+03, new int[] { 0, 12 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{42947D6B-F3CE-474F-9684-F07F36136A5B}"));
                OrthonormalPolynomials[79] = p;
                p.AddCoeff(-1.1243145189457472e+01, new int[] { 1, 1 });
                p.AddCoeff(2.4360147910491190e+02, new int[] { 1, 3 });
                p.AddCoeff(-1.4616088746294714e+03, new int[] { 1, 5 });
                p.AddCoeff(3.5496215526715734e+03, new int[] { 1, 7 });
                p.AddCoeff(-3.7468227500422164e+03, new int[] { 1, 9 });
                p.AddCoeff(1.4306050500161190e+03, new int[] { 1, 11 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{8E7CEF1B-315B-4899-82CA-B46CCF64F7CF}"));
                OrthonormalPolynomials[80] = p;
                p.AddCoeff(6.3042763501509248e-01, new int[] { 0, 0 });
                p.AddCoeff(-3.4673519925830086e+01, new int[] { 0, 2 });
                p.AddCoeff(-1.8912829050452774e+00, new int[] { 2, 0 });
                p.AddCoeff(3.0050383935719408e+02, new int[] { 0, 4 });
                p.AddCoeff(1.0402055977749026e+02, new int[] { 2, 2 });
                p.AddCoeff(-9.0151151807158224e+02, new int[] { 0, 6 });
                p.AddCoeff(-9.0151151807158224e+02, new int[] { 2, 4 });
                p.AddCoeff(1.0946925576583499e+03, new int[] { 0, 8 });
                p.AddCoeff(2.7045345542147467e+03, new int[] { 2, 6 });
                p.AddCoeff(-4.6220352434463661e+02, new int[] { 0, 10 });
                p.AddCoeff(-3.2840776729750496e+03, new int[] { 2, 8 });
                p.AddCoeff(1.3866105730339098e+03, new int[] { 2, 10 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{24A83BB8-3365-429B-81DF-6639F2F6E473}"));
                OrthonormalPolynomials[81] = p;
                p.AddCoeff(-2.1285686820241996e+01, new int[] { 1, 1 });
                p.AddCoeff(3.1219007336354928e+02, new int[] { 1, 3 });
                p.AddCoeff(3.5476144700403327e+01, new int[] { 3, 1 });
                p.AddCoeff(-1.2175412861178422e+03, new int[] { 1, 5 });
                p.AddCoeff(-5.2031678893924880e+02, new int[] { 3, 3 });
                p.AddCoeff(1.7393446944540603e+03, new int[] { 1, 7 });
                p.AddCoeff(2.0292354768630703e+03, new int[] { 3, 5 });
                p.AddCoeff(-8.2135721682552846e+02, new int[] { 1, 9 });
                p.AddCoeff(-2.8989078240901005e+03, new int[] { 3, 7 });
                p.AddCoeff(1.3689286947092141e+03, new int[] { 3, 9 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{807B3E16-68F5-4B41-BB5C-8CAC384BFFFB}"));
                OrthonormalPolynomials[82] = p;
                p.AddCoeff(6.3416907815896634e-01, new int[] { 0, 0 });
                p.AddCoeff(-2.2830086813722788e+01, new int[] { 0, 2 });
                p.AddCoeff(-6.3416907815896634e+00, new int[] { 2, 0 });
                p.AddCoeff(1.2556547747547534e+02, new int[] { 0, 4 });
                p.AddCoeff(2.2830086813722788e+02, new int[] { 2, 2 });
                p.AddCoeff(7.3986392451879407e+00, new int[] { 4, 0 });
                p.AddCoeff(-2.1764682762415725e+02, new int[] { 0, 6 });
                p.AddCoeff(-1.2556547747547534e+03, new int[] { 2, 4 });
                p.AddCoeff(-2.6635101282676586e+02, new int[] { 4, 2 });
                p.AddCoeff(1.1659651479865567e+02, new int[] { 0, 8 });
                p.AddCoeff(2.1764682762415725e+03, new int[] { 2, 6 });
                p.AddCoeff(1.4649305705472123e+03, new int[] { 4, 4 });
                p.AddCoeff(-1.1659651479865567e+03, new int[] { 2, 8 });
                p.AddCoeff(-2.5392129889485012e+03, new int[] { 4, 6 });
                p.AddCoeff(1.3602926726509828e+03, new int[] { 4, 8 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{ACEDD6E7-A434-4771-A0F4-679322372703}"));
                OrthonormalPolynomials[83] = p;
                p.AddCoeff(-2.6342762124215597e+01, new int[] { 1, 1 });
                p.AddCoeff(2.3708485911794037e+02, new int[] { 1, 3 });
                p.AddCoeff(1.2293288991300612e+02, new int[] { 3, 1 });
                p.AddCoeff(-5.2158669005946881e+02, new int[] { 1, 5 });
                p.AddCoeff(-1.1063960092170551e+03, new int[] { 3, 3 });
                p.AddCoeff(-1.1063960092170551e+02, new int[] { 5, 1 });
                p.AddCoeff(3.2288699860824260e+02, new int[] { 1, 7 });
                p.AddCoeff(2.4340712202775211e+03, new int[] { 3, 5 });
                p.AddCoeff(9.9575640829534955e+02, new int[] { 5, 3 });
                p.AddCoeff(-1.5068059935051321e+03, new int[] { 3, 7 });
                p.AddCoeff(-2.1906640982497690e+03, new int[] { 5, 5 });
                p.AddCoeff(1.3561253941546189e+03, new int[] { 5, 7 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{53EC9FA6-5C91-46FE-88DD-764BE34C0712}"));
                OrthonormalPolynomials[84] = p;
                p.AddCoeff(6.3476562500000000e-01, new int[] { 0, 0 });
                p.AddCoeff(-1.3330078125000000e+01, new int[] { 0, 2 });
                p.AddCoeff(-1.3330078125000000e+01, new int[] { 2, 0 });
                p.AddCoeff(3.9990234375000000e+01, new int[] { 0, 4 });
                p.AddCoeff(2.7993164062500000e+02, new int[] { 2, 2 });
                p.AddCoeff(3.9990234375000000e+01, new int[] { 4, 0 });
                p.AddCoeff(-2.9326171875000000e+01, new int[] { 0, 6 });
                p.AddCoeff(-8.3979492187500000e+02, new int[] { 2, 4 });
                p.AddCoeff(-8.3979492187500000e+02, new int[] { 4, 2 });
                p.AddCoeff(-2.9326171875000000e+01, new int[] { 6, 0 });
                p.AddCoeff(6.1584960937500000e+02, new int[] { 2, 6 });
                p.AddCoeff(2.5193847656250000e+03, new int[] { 4, 4 });
                p.AddCoeff(6.1584960937500000e+02, new int[] { 6, 2 });
                p.AddCoeff(-1.8475488281250000e+03, new int[] { 4, 6 });
                p.AddCoeff(-1.8475488281250000e+03, new int[] { 6, 4 });
                p.AddCoeff(1.3548691406250000e+03, new int[] { 6, 6 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{AEADA46B-B631-4846-BBB4-78ABD90B43CA}"));
                OrthonormalPolynomials[85] = p;
                p.AddCoeff(-2.6342762124215597e+01, new int[] { 1, 1 });
                p.AddCoeff(1.2293288991300612e+02, new int[] { 1, 3 });
                p.AddCoeff(2.3708485911794037e+02, new int[] { 3, 1 });
                p.AddCoeff(-1.1063960092170551e+02, new int[] { 1, 5 });
                p.AddCoeff(-1.1063960092170551e+03, new int[] { 3, 3 });
                p.AddCoeff(-5.2158669005946881e+02, new int[] { 5, 1 });
                p.AddCoeff(9.9575640829534955e+02, new int[] { 3, 5 });
                p.AddCoeff(2.4340712202775211e+03, new int[] { 5, 3 });
                p.AddCoeff(3.2288699860824260e+02, new int[] { 7, 1 });
                p.AddCoeff(-2.1906640982497690e+03, new int[] { 5, 5 });
                p.AddCoeff(-1.5068059935051321e+03, new int[] { 7, 3 });
                p.AddCoeff(1.3561253941546189e+03, new int[] { 7, 5 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{07C90A8E-30C2-428D-A0CB-167B6FDD5702}"));
                OrthonormalPolynomials[86] = p;
                p.AddCoeff(6.3416907815896634e-01, new int[] { 0, 0 });
                p.AddCoeff(-6.3416907815896634e+00, new int[] { 0, 2 });
                p.AddCoeff(-2.2830086813722788e+01, new int[] { 2, 0 });
                p.AddCoeff(7.3986392451879407e+00, new int[] { 0, 4 });
                p.AddCoeff(2.2830086813722788e+02, new int[] { 2, 2 });
                p.AddCoeff(1.2556547747547534e+02, new int[] { 4, 0 });
                p.AddCoeff(-2.6635101282676586e+02, new int[] { 2, 4 });
                p.AddCoeff(-1.2556547747547534e+03, new int[] { 4, 2 });
                p.AddCoeff(-2.1764682762415725e+02, new int[] { 6, 0 });
                p.AddCoeff(1.4649305705472123e+03, new int[] { 4, 4 });
                p.AddCoeff(2.1764682762415725e+03, new int[] { 6, 2 });
                p.AddCoeff(1.1659651479865567e+02, new int[] { 8, 0 });
                p.AddCoeff(-2.5392129889485012e+03, new int[] { 6, 4 });
                p.AddCoeff(-1.1659651479865567e+03, new int[] { 8, 2 });
                p.AddCoeff(1.3602926726509828e+03, new int[] { 8, 4 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{4FD83456-373A-4FD3-9196-A8907BA663B8}"));
                OrthonormalPolynomials[87] = p;
                p.AddCoeff(-2.1285686820241996e+01, new int[] { 1, 1 });
                p.AddCoeff(3.5476144700403327e+01, new int[] { 1, 3 });
                p.AddCoeff(3.1219007336354928e+02, new int[] { 3, 1 });
                p.AddCoeff(-5.2031678893924880e+02, new int[] { 3, 3 });
                p.AddCoeff(-1.2175412861178422e+03, new int[] { 5, 1 });
                p.AddCoeff(2.0292354768630703e+03, new int[] { 5, 3 });
                p.AddCoeff(1.7393446944540603e+03, new int[] { 7, 1 });
                p.AddCoeff(-2.8989078240901005e+03, new int[] { 7, 3 });
                p.AddCoeff(-8.2135721682552846e+02, new int[] { 9, 1 });
                p.AddCoeff(1.3689286947092141e+03, new int[] { 9, 3 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{A120F6E8-682F-4ADD-8821-E709718E44EE}"));
                OrthonormalPolynomials[88] = p;
                p.AddCoeff(6.3042763501509248e-01, new int[] { 0, 0 });
                p.AddCoeff(-1.8912829050452774e+00, new int[] { 0, 2 });
                p.AddCoeff(-3.4673519925830086e+01, new int[] { 2, 0 });
                p.AddCoeff(1.0402055977749026e+02, new int[] { 2, 2 });
                p.AddCoeff(3.0050383935719408e+02, new int[] { 4, 0 });
                p.AddCoeff(-9.0151151807158224e+02, new int[] { 4, 2 });
                p.AddCoeff(-9.0151151807158224e+02, new int[] { 6, 0 });
                p.AddCoeff(2.7045345542147467e+03, new int[] { 6, 2 });
                p.AddCoeff(1.0946925576583499e+03, new int[] { 8, 0 });
                p.AddCoeff(-3.2840776729750496e+03, new int[] { 8, 2 });
                p.AddCoeff(-4.6220352434463661e+02, new int[] { 10, 0 });
                p.AddCoeff(1.3866105730339098e+03, new int[] { 10, 2 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{2D3311A9-A092-4953-9FD3-D74BCACDAAA2}"));
                OrthonormalPolynomials[89] = p;
                p.AddCoeff(-1.1243145189457472e+01, new int[] { 1, 1 });
                p.AddCoeff(2.4360147910491190e+02, new int[] { 3, 1 });
                p.AddCoeff(-1.4616088746294714e+03, new int[] { 5, 1 });
                p.AddCoeff(3.5496215526715734e+03, new int[] { 7, 1 });
                p.AddCoeff(-3.7468227500422164e+03, new int[] { 9, 1 });
                p.AddCoeff(1.4306050500161190e+03, new int[] { 11, 1 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{01B4E4B9-60BC-44A4-A2CE-B90A09A08166}"));
                OrthonormalPolynomials[90] = p;
                p.AddCoeff(5.6396484375000000e-01, new int[] { 0, 0 });
                p.AddCoeff(-4.3989257812500000e+01, new int[] { 2, 0 });
                p.AddCoeff(5.4986572265625000e+02, new int[] { 4, 0 });
                p.AddCoeff(-2.4927246093750000e+03, new int[] { 6, 0 });
                p.AddCoeff(5.0744750976562500e+03, new int[] { 8, 0 });
                p.AddCoeff(-4.7361767578125000e+03, new int[] { 10, 0 });
                p.AddCoeff(1.6504858398437500e+03, new int[] { 12, 0 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{94D14835-507F-43F3-BF0E-16334ABA2701}"));
                OrthonormalPolynomials[91] = p;
                p.AddCoeff(7.6191629518496170e+00, new int[] { 0, 1 });
                p.AddCoeff(-2.2857488855548851e+02, new int[] { 0, 3 });
                p.AddCoeff(1.9428865527216523e+03, new int[] { 0, 5 });
                p.AddCoeff(-7.0313989527069322e+03, new int[] { 0, 7 });
                p.AddCoeff(1.2304948167237131e+04, new int[] { 0, 9 });
                p.AddCoeff(-1.0291411194416510e+04, new int[] { 0, 11 });
                p.AddCoeff(3.2985292289796506e+03, new int[] { 0, 13 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{04785C4B-25BF-441A-84CA-6B4BE6C4F6E2}"));
                OrthonormalPolynomials[92] = p;
                p.AddCoeff(9.7681576305764320e-01, new int[] { 1, 0 });
                p.AddCoeff(-7.6191629518496170e+01, new int[] { 1, 2 });
                p.AddCoeff(9.5239536898120212e+02, new int[] { 1, 4 });
                p.AddCoeff(-4.3175256727147829e+03, new int[] { 1, 6 });
                p.AddCoeff(8.7892486908836653e+03, new int[] { 1, 8 });
                p.AddCoeff(-8.2032987781580876e+03, new int[] { 1, 10 });
                p.AddCoeff(2.8587253317823639e+03, new int[] { 1, 12 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{CC0C27F5-4595-44B3-B5A0-7E8C6337E257}"));
                OrthonormalPolynomials[93] = p;
                p.AddCoeff(7.2574190129601373e+00, new int[] { 0, 1 });
                p.AddCoeff(-1.5724407861413631e+02, new int[] { 0, 3 });
                p.AddCoeff(-2.1772257038880412e+01, new int[] { 2, 1 });
                p.AddCoeff(9.4346447168481784e+02, new int[] { 0, 5 });
                p.AddCoeff(4.7173223584240892e+02, new int[] { 2, 3 });
                p.AddCoeff(-2.2912708598059862e+03, new int[] { 0, 7 });
                p.AddCoeff(-2.8303934150544535e+03, new int[] { 2, 5 });
                p.AddCoeff(2.4185636853507632e+03, new int[] { 0, 9 });
                p.AddCoeff(6.8738125794179586e+03, new int[] { 2, 7 });
                p.AddCoeff(-9.2345158895210959e+02, new int[] { 0, 11 });
                p.AddCoeff(-7.2556910560522896e+03, new int[] { 2, 9 });
                p.AddCoeff(2.7703547668563288e+03, new int[] { 2, 11 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{5B453737-5E91-4B10-9AFC-4691CA7D460E}"));
                OrthonormalPolynomials[94] = p;
                p.AddCoeff(2.2377961117320553e+00, new int[] { 1, 0 });
                p.AddCoeff(-1.2307878614526304e+02, new int[] { 1, 2 });
                p.AddCoeff(-3.7296601862200922e+00, new int[] { 3, 0 });
                p.AddCoeff(1.0666828132589464e+03, new int[] { 1, 4 });
                p.AddCoeff(2.0513131024210507e+02, new int[] { 3, 2 });
                p.AddCoeff(-3.2000484397768391e+03, new int[] { 1, 6 });
                p.AddCoeff(-1.7778046887649106e+03, new int[] { 3, 4 });
                p.AddCoeff(3.8857731054433047e+03, new int[] { 1, 8 });
                p.AddCoeff(5.3334140662947319e+03, new int[] { 3, 6 });
                p.AddCoeff(-1.6406597556316175e+03, new int[] { 1, 10 });
                p.AddCoeff(-6.4762885090721744e+03, new int[] { 3, 8 });
                p.AddCoeff(2.7344329260526959e+03, new int[] { 3, 10 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{589ACD96-55A1-4C43-862D-476A0A4E56F3}"));
                OrthonormalPolynomials[95] = p;
                p.AddCoeff(6.0339250512391648e+00, new int[] { 0, 1 });
                p.AddCoeff(-8.8497567418174417e+01, new int[] { 0, 3 });
                p.AddCoeff(-6.0339250512391648e+01, new int[] { 2, 1 });
                p.AddCoeff(3.4514051293088023e+02, new int[] { 0, 5 });
                p.AddCoeff(8.8497567418174417e+02, new int[] { 2, 3 });
                p.AddCoeff(7.0395792264456923e+01, new int[] { 4, 1 });
                p.AddCoeff(-4.9305787561554318e+02, new int[] { 0, 7 });
                p.AddCoeff(-3.4514051293088023e+03, new int[] { 2, 5 });
                p.AddCoeff(-1.0324716198787015e+03, new int[] { 4, 3 });
                p.AddCoeff(2.3283288570733984e+02, new int[] { 0, 9 });
                p.AddCoeff(4.9305787561554318e+03, new int[] { 2, 7 });
                p.AddCoeff(4.0266393175269360e+03, new int[] { 4, 5 });
                p.AddCoeff(-2.3283288570733984e+03, new int[] { 2, 9 });
                p.AddCoeff(-5.7523418821813371e+03, new int[] { 4, 7 });
                p.AddCoeff(2.7163836665856314e+03, new int[] { 4, 9 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{F56C18F1-77E8-4048-90DA-8E21D0E75D4F}"));
                OrthonormalPolynomials[96] = p;
                p.AddCoeff(3.5055014764980982e+00, new int[] { 1, 0 });
                p.AddCoeff(-1.2619805315393154e+02, new int[] { 1, 2 });
                p.AddCoeff(-1.6359006890324458e+01, new int[] { 3, 0 });
                p.AddCoeff(6.9408929234662345e+02, new int[] { 1, 4 });
                p.AddCoeff(5.8892424805168050e+02, new int[] { 3, 2 });
                p.AddCoeff(1.4723106201292013e+01, new int[] { 5, 0 });
                p.AddCoeff(-1.2030881067341473e+03, new int[] { 1, 6 });
                p.AddCoeff(-3.2390833642842428e+03, new int[] { 3, 4 });
                p.AddCoeff(-5.3003182324651245e+02, new int[] { 5, 2 });
                p.AddCoeff(6.4451148575043606e+02, new int[] { 1, 8 });
                p.AddCoeff(5.6144111647593541e+03, new int[] { 3, 6 });
                p.AddCoeff(2.9151750278558185e+03, new int[] { 5, 4 });
                p.AddCoeff(-3.0077202668353683e+03, new int[] { 3, 8 });
                p.AddCoeff(-5.0529700482834187e+03, new int[] { 5, 6 });
                p.AddCoeff(2.7069482401518314e+03, new int[] { 5, 8 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{F2E9E993-F572-4A84-8C3D-F4700C22CC3A}"));
                OrthonormalPolynomials[97] = p;
                p.AddCoeff(4.7729336087100873e+00, new int[] { 0, 1 });
                p.AddCoeff(-4.2956402478390786e+01, new int[] { 0, 3 });
                p.AddCoeff(-1.0023160578291183e+02, new int[] { 2, 1 });
                p.AddCoeff(9.4504085452459729e+01, new int[] { 0, 5 });
                p.AddCoeff(9.0208445204620650e+02, new int[] { 2, 3 });
                p.AddCoeff(3.0069481734873550e+02, new int[] { 4, 1 });
                p.AddCoeff(-5.8502529089617927e+01, new int[] { 0, 7 });
                p.AddCoeff(-1.9845857945016543e+03, new int[] { 2, 5 });
                p.AddCoeff(-2.7062533561386195e+03, new int[] { 4, 3 });
                p.AddCoeff(-2.2050953272240603e+02, new int[] { 6, 1 });
                p.AddCoeff(1.2285531108819765e+03, new int[] { 2, 7 });
                p.AddCoeff(5.9537573835049629e+03, new int[] { 4, 5 });
                p.AddCoeff(1.9845857945016543e+03, new int[] { 6, 3 });
                p.AddCoeff(-3.6856593326459294e+03, new int[] { 4, 7 });
                p.AddCoeff(-4.3660887479036395e+03, new int[] { 6, 5 });
                p.AddCoeff(2.7028168439403482e+03, new int[] { 6, 7 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{B9FB3542-22AE-4B70-B726-B2BDCE6C98F8}"));
                OrthonormalPolynomials[98] = p;
                p.AddCoeff(4.7729336087100873e+00, new int[] { 1, 0 });
                p.AddCoeff(-1.0023160578291183e+02, new int[] { 1, 2 });
                p.AddCoeff(-4.2956402478390786e+01, new int[] { 3, 0 });
                p.AddCoeff(3.0069481734873550e+02, new int[] { 1, 4 });
                p.AddCoeff(9.0208445204620650e+02, new int[] { 3, 2 });
                p.AddCoeff(9.4504085452459729e+01, new int[] { 5, 0 });
                p.AddCoeff(-2.2050953272240603e+02, new int[] { 1, 6 });
                p.AddCoeff(-2.7062533561386195e+03, new int[] { 3, 4 });
                p.AddCoeff(-1.9845857945016543e+03, new int[] { 5, 2 });
                p.AddCoeff(-5.8502529089617927e+01, new int[] { 7, 0 });
                p.AddCoeff(1.9845857945016543e+03, new int[] { 3, 6 });
                p.AddCoeff(5.9537573835049629e+03, new int[] { 5, 4 });
                p.AddCoeff(1.2285531108819765e+03, new int[] { 7, 2 });
                p.AddCoeff(-4.3660887479036395e+03, new int[] { 5, 6 });
                p.AddCoeff(-3.6856593326459294e+03, new int[] { 7, 4 });
                p.AddCoeff(2.7028168439403482e+03, new int[] { 7, 6 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{CBD61C58-4100-4071-B5F8-26597BE865CB}"));
                OrthonormalPolynomials[99] = p;
                p.AddCoeff(3.5055014764980982e+00, new int[] { 0, 1 });
                p.AddCoeff(-1.6359006890324458e+01, new int[] { 0, 3 });
                p.AddCoeff(-1.2619805315393154e+02, new int[] { 2, 1 });
                p.AddCoeff(1.4723106201292013e+01, new int[] { 0, 5 });
                p.AddCoeff(5.8892424805168050e+02, new int[] { 2, 3 });
                p.AddCoeff(6.9408929234662345e+02, new int[] { 4, 1 });
                p.AddCoeff(-5.3003182324651245e+02, new int[] { 2, 5 });
                p.AddCoeff(-3.2390833642842428e+03, new int[] { 4, 3 });
                p.AddCoeff(-1.2030881067341473e+03, new int[] { 6, 1 });
                p.AddCoeff(2.9151750278558185e+03, new int[] { 4, 5 });
                p.AddCoeff(5.6144111647593541e+03, new int[] { 6, 3 });
                p.AddCoeff(6.4451148575043606e+02, new int[] { 8, 1 });
                p.AddCoeff(-5.0529700482834187e+03, new int[] { 6, 5 });
                p.AddCoeff(-3.0077202668353683e+03, new int[] { 8, 3 });
                p.AddCoeff(2.7069482401518314e+03, new int[] { 8, 5 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{9A38141A-0328-49D4-B10B-F1A44FBA59C8}"));
                OrthonormalPolynomials[100] = p;
                p.AddCoeff(6.0339250512391648e+00, new int[] { 1, 0 });
                p.AddCoeff(-6.0339250512391648e+01, new int[] { 1, 2 });
                p.AddCoeff(-8.8497567418174417e+01, new int[] { 3, 0 });
                p.AddCoeff(7.0395792264456923e+01, new int[] { 1, 4 });
                p.AddCoeff(8.8497567418174417e+02, new int[] { 3, 2 });
                p.AddCoeff(3.4514051293088023e+02, new int[] { 5, 0 });
                p.AddCoeff(-1.0324716198787015e+03, new int[] { 3, 4 });
                p.AddCoeff(-3.4514051293088023e+03, new int[] { 5, 2 });
                p.AddCoeff(-4.9305787561554318e+02, new int[] { 7, 0 });
                p.AddCoeff(4.0266393175269360e+03, new int[] { 5, 4 });
                p.AddCoeff(4.9305787561554318e+03, new int[] { 7, 2 });
                p.AddCoeff(2.3283288570733984e+02, new int[] { 9, 0 });
                p.AddCoeff(-5.7523418821813371e+03, new int[] { 7, 4 });
                p.AddCoeff(-2.3283288570733984e+03, new int[] { 9, 2 });
                p.AddCoeff(2.7163836665856314e+03, new int[] { 9, 4 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{A55AABB5-A09F-476A-8DB7-68CC1A695DEB}"));
                OrthonormalPolynomials[101] = p;
                p.AddCoeff(2.2377961117320553e+00, new int[] { 0, 1 });
                p.AddCoeff(-3.7296601862200922e+00, new int[] { 0, 3 });
                p.AddCoeff(-1.2307878614526304e+02, new int[] { 2, 1 });
                p.AddCoeff(2.0513131024210507e+02, new int[] { 2, 3 });
                p.AddCoeff(1.0666828132589464e+03, new int[] { 4, 1 });
                p.AddCoeff(-1.7778046887649106e+03, new int[] { 4, 3 });
                p.AddCoeff(-3.2000484397768391e+03, new int[] { 6, 1 });
                p.AddCoeff(5.3334140662947319e+03, new int[] { 6, 3 });
                p.AddCoeff(3.8857731054433047e+03, new int[] { 8, 1 });
                p.AddCoeff(-6.4762885090721744e+03, new int[] { 8, 3 });
                p.AddCoeff(-1.6406597556316175e+03, new int[] { 10, 1 });
                p.AddCoeff(2.7344329260526959e+03, new int[] { 10, 3 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{1E891FDC-2A3A-4D4D-851F-651E423F094C}"));
                OrthonormalPolynomials[102] = p;
                p.AddCoeff(7.2574190129601373e+00, new int[] { 1, 0 });
                p.AddCoeff(-2.1772257038880412e+01, new int[] { 1, 2 });
                p.AddCoeff(-1.5724407861413631e+02, new int[] { 3, 0 });
                p.AddCoeff(4.7173223584240892e+02, new int[] { 3, 2 });
                p.AddCoeff(9.4346447168481784e+02, new int[] { 5, 0 });
                p.AddCoeff(-2.8303934150544535e+03, new int[] { 5, 2 });
                p.AddCoeff(-2.2912708598059862e+03, new int[] { 7, 0 });
                p.AddCoeff(6.8738125794179586e+03, new int[] { 7, 2 });
                p.AddCoeff(2.4185636853507632e+03, new int[] { 9, 0 });
                p.AddCoeff(-7.2556910560522896e+03, new int[] { 9, 2 });
                p.AddCoeff(-9.2345158895210959e+02, new int[] { 11, 0 });
                p.AddCoeff(2.7703547668563288e+03, new int[] { 11, 2 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{B96F4081-6113-485C-B7E6-9EE7C02AE1E0}"));
                OrthonormalPolynomials[103] = p;
                p.AddCoeff(9.7681576305764320e-01, new int[] { 0, 1 });
                p.AddCoeff(-7.6191629518496170e+01, new int[] { 2, 1 });
                p.AddCoeff(9.5239536898120212e+02, new int[] { 4, 1 });
                p.AddCoeff(-4.3175256727147829e+03, new int[] { 6, 1 });
                p.AddCoeff(8.7892486908836653e+03, new int[] { 8, 1 });
                p.AddCoeff(-8.2032987781580876e+03, new int[] { 10, 1 });
                p.AddCoeff(2.8587253317823639e+03, new int[] { 12, 1 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{8A312D13-8AEC-41EE-858C-F9727F5593EC}"));
                OrthonormalPolynomials[104] = p;
                p.AddCoeff(7.6191629518496170e+00, new int[] { 1, 0 });
                p.AddCoeff(-2.2857488855548851e+02, new int[] { 3, 0 });
                p.AddCoeff(1.9428865527216523e+03, new int[] { 5, 0 });
                p.AddCoeff(-7.0313989527069322e+03, new int[] { 7, 0 });
                p.AddCoeff(1.2304948167237131e+04, new int[] { 9, 0 });
                p.AddCoeff(-1.0291411194416510e+04, new int[] { 11, 0 });
                p.AddCoeff(3.2985292289796506e+03, new int[] { 13, 0 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{E3E9D69F-C922-4DE3-BF19-20A8419961E9}"));
                OrthonormalPolynomials[105] = p;
                p.AddCoeff(-5.6402238824724176e-01, new int[] { 0, 0 });
                p.AddCoeff(5.9222350765960384e+01, new int[] { 0, 2 });
                p.AddCoeff(-1.0067799630213265e+03, new int[] { 0, 4 });
                p.AddCoeff(6.3762730991350680e+03, new int[] { 0, 6 });
                p.AddCoeff(-1.9128819297405204e+04, new int[] { 0, 8 });
                p.AddCoeff(2.9330856256021313e+04, new int[] { 0, 10 });
                p.AddCoeff(-2.2220345648500995e+04, new int[] { 0, 12 });
                p.AddCoeff(6.5928498077969984e+03, new int[] { 0, 14 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{D1758F63-69C3-48BF-BDF9-F84476E97CAA}"));
                OrthonormalPolynomials[106] = p;
                p.AddCoeff(1.3196777343750000e+01, new int[] { 1, 1 });
                p.AddCoeff(-3.9590332031250000e+02, new int[] { 1, 3 });
                p.AddCoeff(3.3651782226562500e+03, new int[] { 1, 5 });
                p.AddCoeff(-1.2178740234375000e+04, new int[] { 1, 7 });
                p.AddCoeff(2.1312795410156250e+04, new int[] { 1, 9 });
                p.AddCoeff(-1.7825247070312500e+04, new int[] { 1, 11 });
                p.AddCoeff(5.7132202148437500e+03, new int[] { 1, 13 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{42929A62-19D7-4D27-BA61-6778B56E7A0A}"));
                OrthonormalPolynomials[107] = p;
                p.AddCoeff(-6.3053186377252371e-01, new int[] { 0, 0 });
                p.AddCoeff(4.9181485374256849e+01, new int[] { 0, 2 });
                p.AddCoeff(1.8915955913175711e+00, new int[] { 2, 0 });
                p.AddCoeff(-6.1476856717821061e+02, new int[] { 0, 4 });
                p.AddCoeff(-1.4754445612277055e+02, new int[] { 2, 2 });
                p.AddCoeff(2.7869508378745548e+03, new int[] { 0, 6 });
                p.AddCoeff(1.8443057015346318e+03, new int[] { 2, 4 });
                p.AddCoeff(-5.6734356342446294e+03, new int[] { 0, 8 });
                p.AddCoeff(-8.3608525136236643e+03, new int[] { 2, 6 });
                p.AddCoeff(5.2952065919616541e+03, new int[] { 0, 10 });
                p.AddCoeff(1.7020306902733888e+04, new int[] { 2, 8 });
                p.AddCoeff(-1.8452992668957279e+03, new int[] { 0, 12 });
                p.AddCoeff(-1.5885619775884962e+04, new int[] { 2, 10 });
                p.AddCoeff(5.5358978006871838e+03, new int[] { 2, 12 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{5211634F-B43F-43F4-A3A6-0E55542134CE}"));
                OrthonormalPolynomials[108] = p;
                p.AddCoeff(2.5761281940033743e+01, new int[] { 1, 1 });
                p.AddCoeff(-5.5816110870073110e+02, new int[] { 1, 3 });
                p.AddCoeff(-4.2935469900056238e+01, new int[] { 3, 1 });
                p.AddCoeff(3.3489666522043866e+03, new int[] { 1, 5 });
                p.AddCoeff(9.3026851450121850e+02, new int[] { 3, 3 });
                p.AddCoeff(-8.1332047267820817e+03, new int[] { 1, 7 });
                p.AddCoeff(-5.5816110870073110e+03, new int[] { 3, 5 });
                p.AddCoeff(8.5850494338255307e+03, new int[] { 1, 9 });
                p.AddCoeff(1.3555341211303470e+04, new int[] { 3, 7 });
                p.AddCoeff(-3.2779279656424754e+03, new int[] { 1, 11 });
                p.AddCoeff(-1.4308415723042551e+04, new int[] { 3, 9 });
                p.AddCoeff(5.4632132760707923e+03, new int[] { 3, 11 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{304962D1-02BC-4223-87B8-44A6CE7A7101}"));
                OrthonormalPolynomials[109] = p;
                p.AddCoeff(-6.3435557105467805e-01, new int[] { 0, 0 });
                p.AddCoeff(3.4889556408007293e+01, new int[] { 0, 2 });
                p.AddCoeff(6.3435557105467805e+00, new int[] { 2, 0 });
                p.AddCoeff(-3.0237615553606320e+02, new int[] { 0, 4 });
                p.AddCoeff(-3.4889556408007293e+02, new int[] { 2, 2 });
                p.AddCoeff(-7.4008149956379106e+00, new int[] { 4, 0 });
                p.AddCoeff(9.0712846660818961e+02, new int[] { 0, 6 });
                p.AddCoeff(3.0237615553606320e+03, new int[] { 2, 4 });
                p.AddCoeff(4.0704482476008508e+02, new int[] { 4, 2 });
                p.AddCoeff(-1.1015131380242302e+03, new int[] { 0, 8 });
                p.AddCoeff(-9.0712846660818961e+03, new int[] { 2, 6 });
                p.AddCoeff(-3.5277218145874040e+03, new int[] { 4, 4 });
                p.AddCoeff(4.6508332494356388e+02, new int[] { 0, 10 });
                p.AddCoeff(1.1015131380242302e+04, new int[] { 2, 8 });
                p.AddCoeff(1.0583165443762212e+04, new int[] { 4, 6 });
                p.AddCoeff(-4.6508332494356388e+03, new int[] { 2, 10 });
                p.AddCoeff(-1.2850986610282686e+04, new int[] { 4, 8 });
                p.AddCoeff(5.4259721243415786e+03, new int[] { 4, 10 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{5375DA86-487D-474C-A123-9BB3BA741297}"));
                OrthonormalPolynomials[110] = p;
                p.AddCoeff(3.3353775680143817e+01, new int[] { 1, 1 });
                p.AddCoeff(-4.8918870997544265e+02, new int[] { 1, 3 });
                p.AddCoeff(-1.5565095317400448e+02, new int[] { 3, 1 });
                p.AddCoeff(1.9078359689042263e+03, new int[] { 1, 5 });
                p.AddCoeff(2.2828806465520657e+03, new int[] { 3, 3 });
                p.AddCoeff(1.4008585785660403e+02, new int[] { 5, 1 });
                p.AddCoeff(-2.7254799555774662e+03, new int[] { 1, 7 });
                p.AddCoeff(-8.9032345215530562e+03, new int[] { 3, 5 });
                p.AddCoeff(-2.0545925818968591e+03, new int[] { 5, 3 });
                p.AddCoeff(1.2870322012449146e+03, new int[] { 1, 9 });
                p.AddCoeff(1.2718906459361509e+04, new int[] { 3, 7 });
                p.AddCoeff(8.0129110693977506e+03, new int[] { 5, 5 });
                p.AddCoeff(-6.0061502724762681e+03, new int[] { 3, 9 });
                p.AddCoeff(-1.1447015813425358e+04, new int[] { 5, 7 });
                p.AddCoeff(5.4055352452286413e+03, new int[] { 5, 9 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{0BB37EB9-450C-4C55-BF16-EA1C14D0CB83}"));
                OrthonormalPolynomials[111] = p;
                p.AddCoeff(-6.3514698017107873e-01, new int[] { 0, 0 });
                p.AddCoeff(2.2865291286158834e+01, new int[] { 0, 2 });
                p.AddCoeff(1.3338086583592653e+01, new int[] { 2, 0 });
                p.AddCoeff(-1.2575910207387359e+02, new int[] { 0, 4 });
                p.AddCoeff(-4.8017111700933552e+02, new int[] { 2, 2 });
                p.AddCoeff(-4.0014259750777960e+01, new int[] { 4, 0 });
                p.AddCoeff(2.1798244359471422e+02, new int[] { 0, 6 });
                p.AddCoeff(2.6409411435513453e+03, new int[] { 2, 4 });
                p.AddCoeff(1.4405133510280066e+03, new int[] { 4, 2 });
                p.AddCoeff(2.9343790483903837e+01, new int[] { 6, 0 });
                p.AddCoeff(-1.1677630906859690e+02, new int[] { 0, 8 });
                p.AddCoeff(-4.5776313154889986e+03, new int[] { 2, 6 });
                p.AddCoeff(-7.9228234306540360e+03, new int[] { 4, 4 });
                p.AddCoeff(-1.0563764574205381e+03, new int[] { 6, 2 });
                p.AddCoeff(2.4523024904405350e+03, new int[] { 2, 8 });
                p.AddCoeff(1.3732893946466996e+04, new int[] { 4, 6 });
                p.AddCoeff(5.8100705158129598e+03, new int[] { 6, 4 });
                p.AddCoeff(-7.3569074713216049e+03, new int[] { 4, 8 });
                p.AddCoeff(-1.0070788894075797e+04, new int[] { 6, 6 });
                p.AddCoeff(5.3950654789691769e+03, new int[] { 6, 8 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{639D619F-5D29-4A18-8EAC-0F17F93B2A0D}"));
                OrthonormalPolynomials[112] = p;
                p.AddCoeff(3.5888671875000000e+01, new int[] { 1, 1 });
                p.AddCoeff(-3.2299804687500000e+02, new int[] { 1, 3 });
                p.AddCoeff(-3.2299804687500000e+02, new int[] { 3, 1 });
                p.AddCoeff(7.1059570312500000e+02, new int[] { 1, 5 });
                p.AddCoeff(2.9069824218750000e+03, new int[] { 3, 3 });
                p.AddCoeff(7.1059570312500000e+02, new int[] { 5, 1 });
                p.AddCoeff(-4.3989257812500000e+02, new int[] { 1, 7 });
                p.AddCoeff(-6.3953613281250000e+03, new int[] { 3, 5 });
                p.AddCoeff(-6.3953613281250000e+03, new int[] { 5, 3 });
                p.AddCoeff(-4.3989257812500000e+02, new int[] { 7, 1 });
                p.AddCoeff(3.9590332031250000e+03, new int[] { 3, 7 });
                p.AddCoeff(1.4069794921875000e+04, new int[] { 5, 5 });
                p.AddCoeff(3.9590332031250000e+03, new int[] { 7, 3 });
                p.AddCoeff(-8.7098730468750000e+03, new int[] { 5, 7 });
                p.AddCoeff(-8.7098730468750000e+03, new int[] { 7, 5 });
                p.AddCoeff(5.3918261718750000e+03, new int[] { 7, 7 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{35014D41-442E-4CB6-8B13-1BFA729764C9}"));
                OrthonormalPolynomials[113] = p;
                p.AddCoeff(-6.3514698017107873e-01, new int[] { 0, 0 });
                p.AddCoeff(1.3338086583592653e+01, new int[] { 0, 2 });
                p.AddCoeff(2.2865291286158834e+01, new int[] { 2, 0 });
                p.AddCoeff(-4.0014259750777960e+01, new int[] { 0, 4 });
                p.AddCoeff(-4.8017111700933552e+02, new int[] { 2, 2 });
                p.AddCoeff(-1.2575910207387359e+02, new int[] { 4, 0 });
                p.AddCoeff(2.9343790483903837e+01, new int[] { 0, 6 });
                p.AddCoeff(1.4405133510280066e+03, new int[] { 2, 4 });
                p.AddCoeff(2.6409411435513453e+03, new int[] { 4, 2 });
                p.AddCoeff(2.1798244359471422e+02, new int[] { 6, 0 });
                p.AddCoeff(-1.0563764574205381e+03, new int[] { 2, 6 });
                p.AddCoeff(-7.9228234306540360e+03, new int[] { 4, 4 });
                p.AddCoeff(-4.5776313154889986e+03, new int[] { 6, 2 });
                p.AddCoeff(-1.1677630906859690e+02, new int[] { 8, 0 });
                p.AddCoeff(5.8100705158129598e+03, new int[] { 4, 6 });
                p.AddCoeff(1.3732893946466996e+04, new int[] { 6, 4 });
                p.AddCoeff(2.4523024904405350e+03, new int[] { 8, 2 });
                p.AddCoeff(-1.0070788894075797e+04, new int[] { 6, 6 });
                p.AddCoeff(-7.3569074713216049e+03, new int[] { 8, 4 });
                p.AddCoeff(5.3950654789691769e+03, new int[] { 8, 6 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{BF877274-97BC-40D7-B0E4-19D00ECC2DBA}"));
                OrthonormalPolynomials[114] = p;
                p.AddCoeff(3.3353775680143817e+01, new int[] { 1, 1 });
                p.AddCoeff(-1.5565095317400448e+02, new int[] { 1, 3 });
                p.AddCoeff(-4.8918870997544265e+02, new int[] { 3, 1 });
                p.AddCoeff(1.4008585785660403e+02, new int[] { 1, 5 });
                p.AddCoeff(2.2828806465520657e+03, new int[] { 3, 3 });
                p.AddCoeff(1.9078359689042263e+03, new int[] { 5, 1 });
                p.AddCoeff(-2.0545925818968591e+03, new int[] { 3, 5 });
                p.AddCoeff(-8.9032345215530562e+03, new int[] { 5, 3 });
                p.AddCoeff(-2.7254799555774662e+03, new int[] { 7, 1 });
                p.AddCoeff(8.0129110693977506e+03, new int[] { 5, 5 });
                p.AddCoeff(1.2718906459361509e+04, new int[] { 7, 3 });
                p.AddCoeff(1.2870322012449146e+03, new int[] { 9, 1 });
                p.AddCoeff(-1.1447015813425358e+04, new int[] { 7, 5 });
                p.AddCoeff(-6.0061502724762681e+03, new int[] { 9, 3 });
                p.AddCoeff(5.4055352452286413e+03, new int[] { 9, 5 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{0B09C9CB-76C7-4E05-86D6-20E7E42CEA56}"));
                OrthonormalPolynomials[115] = p;
                p.AddCoeff(-6.3435557105467805e-01, new int[] { 0, 0 });
                p.AddCoeff(6.3435557105467805e+00, new int[] { 0, 2 });
                p.AddCoeff(3.4889556408007293e+01, new int[] { 2, 0 });
                p.AddCoeff(-7.4008149956379106e+00, new int[] { 0, 4 });
                p.AddCoeff(-3.4889556408007293e+02, new int[] { 2, 2 });
                p.AddCoeff(-3.0237615553606320e+02, new int[] { 4, 0 });
                p.AddCoeff(4.0704482476008508e+02, new int[] { 2, 4 });
                p.AddCoeff(3.0237615553606320e+03, new int[] { 4, 2 });
                p.AddCoeff(9.0712846660818961e+02, new int[] { 6, 0 });
                p.AddCoeff(-3.5277218145874040e+03, new int[] { 4, 4 });
                p.AddCoeff(-9.0712846660818961e+03, new int[] { 6, 2 });
                p.AddCoeff(-1.1015131380242302e+03, new int[] { 8, 0 });
                p.AddCoeff(1.0583165443762212e+04, new int[] { 6, 4 });
                p.AddCoeff(1.1015131380242302e+04, new int[] { 8, 2 });
                p.AddCoeff(4.6508332494356388e+02, new int[] { 10, 0 });
                p.AddCoeff(-1.2850986610282686e+04, new int[] { 8, 4 });
                p.AddCoeff(-4.6508332494356388e+03, new int[] { 10, 2 });
                p.AddCoeff(5.4259721243415786e+03, new int[] { 10, 4 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{68FDA7BE-CA14-4F93-B57C-98A432921EDD}"));
                OrthonormalPolynomials[116] = p;
                p.AddCoeff(2.5761281940033743e+01, new int[] { 1, 1 });
                p.AddCoeff(-4.2935469900056238e+01, new int[] { 1, 3 });
                p.AddCoeff(-5.5816110870073110e+02, new int[] { 3, 1 });
                p.AddCoeff(9.3026851450121850e+02, new int[] { 3, 3 });
                p.AddCoeff(3.3489666522043866e+03, new int[] { 5, 1 });
                p.AddCoeff(-5.5816110870073110e+03, new int[] { 5, 3 });
                p.AddCoeff(-8.1332047267820817e+03, new int[] { 7, 1 });
                p.AddCoeff(1.3555341211303470e+04, new int[] { 7, 3 });
                p.AddCoeff(8.5850494338255307e+03, new int[] { 9, 1 });
                p.AddCoeff(-1.4308415723042551e+04, new int[] { 9, 3 });
                p.AddCoeff(-3.2779279656424754e+03, new int[] { 11, 1 });
                p.AddCoeff(5.4632132760707923e+03, new int[] { 11, 3 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{531C6823-E198-404D-AF92-948343D6C041}"));
                OrthonormalPolynomials[117] = p;
                p.AddCoeff(-6.3053186377252371e-01, new int[] { 0, 0 });
                p.AddCoeff(1.8915955913175711e+00, new int[] { 0, 2 });
                p.AddCoeff(4.9181485374256849e+01, new int[] { 2, 0 });
                p.AddCoeff(-1.4754445612277055e+02, new int[] { 2, 2 });
                p.AddCoeff(-6.1476856717821061e+02, new int[] { 4, 0 });
                p.AddCoeff(1.8443057015346318e+03, new int[] { 4, 2 });
                p.AddCoeff(2.7869508378745548e+03, new int[] { 6, 0 });
                p.AddCoeff(-8.3608525136236643e+03, new int[] { 6, 2 });
                p.AddCoeff(-5.6734356342446294e+03, new int[] { 8, 0 });
                p.AddCoeff(1.7020306902733888e+04, new int[] { 8, 2 });
                p.AddCoeff(5.2952065919616541e+03, new int[] { 10, 0 });
                p.AddCoeff(-1.5885619775884962e+04, new int[] { 10, 2 });
                p.AddCoeff(-1.8452992668957279e+03, new int[] { 12, 0 });
                p.AddCoeff(5.5358978006871838e+03, new int[] { 12, 2 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{93448F52-F0AE-438B-869E-864E6BC0D02A}"));
                OrthonormalPolynomials[118] = p;
                p.AddCoeff(1.3196777343750000e+01, new int[] { 1, 1 });
                p.AddCoeff(-3.9590332031250000e+02, new int[] { 3, 1 });
                p.AddCoeff(3.3651782226562500e+03, new int[] { 5, 1 });
                p.AddCoeff(-1.2178740234375000e+04, new int[] { 7, 1 });
                p.AddCoeff(2.1312795410156250e+04, new int[] { 9, 1 });
                p.AddCoeff(-1.7825247070312500e+04, new int[] { 11, 1 });
                p.AddCoeff(5.7132202148437500e+03, new int[] { 13, 1 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{1068F3B4-EA06-4960-AAF4-F6D72F9F17C3}"));
                OrthonormalPolynomials[119] = p;
                p.AddCoeff(-5.6402238824724176e-01, new int[] { 0, 0 });
                p.AddCoeff(5.9222350765960384e+01, new int[] { 2, 0 });
                p.AddCoeff(-1.0067799630213265e+03, new int[] { 4, 0 });
                p.AddCoeff(6.3762730991350680e+03, new int[] { 6, 0 });
                p.AddCoeff(-1.9128819297405204e+04, new int[] { 8, 0 });
                p.AddCoeff(2.9330856256021313e+04, new int[] { 10, 0 });
                p.AddCoeff(-2.2220345648500995e+04, new int[] { 12, 0 });
                p.AddCoeff(6.5928498077969984e+03, new int[] { 14, 0 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{C3C15257-DC1D-4FE1-A48A-88FF2A13E3DF}"));
                OrthonormalPolynomials[120] = p;
                p.AddCoeff(-8.7472079284207009e+00, new int[] { 0, 1 });
                p.AddCoeff(3.4697258116068780e+02, new int[] { 0, 3 });
                p.AddCoeff(-3.9554874252318410e+03, new int[] { 0, 5 });
                p.AddCoeff(1.9777437126159205e+04, new int[] { 0, 7 });
                p.AddCoeff(-5.0542339322406857e+04, new int[] { 0, 9 });
                p.AddCoeff(6.8921371803282077e+04, new int[] { 0, 11 });
                p.AddCoeff(-4.7714795863810669e+04, new int[] { 0, 13 });
                p.AddCoeff(1.3178372190957232e+04, new int[] { 0, 15 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{86924C8E-B987-4D24-B103-2FC8CB07ECF7}"));
                OrthonormalPolynomials[121] = p;
                p.AddCoeff(-9.7691543305056193e-01, new int[] { 1, 0 });
                p.AddCoeff(1.0257612047030900e+02, new int[] { 1, 2 });
                p.AddCoeff(-1.7437940479952530e+03, new int[] { 1, 4 });
                p.AddCoeff(1.1044028970636603e+04, new int[] { 1, 6 });
                p.AddCoeff(-3.3132086911909808e+04, new int[] { 1, 8 });
                p.AddCoeff(5.0802533264928372e+04, new int[] { 1, 10 });
                p.AddCoeff(-3.8486767624945736e+04, new int[] { 1, 12 });
                p.AddCoeff(1.1419150833775109e+04, new int[] { 1, 14 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{7EA4BF2E-0552-43BC-A208-0A2195341812}"));
                OrthonormalPolynomials[122] = p;
                p.AddCoeff(-8.5184831459918503e+00, new int[] { 0, 1 });
                p.AddCoeff(2.5555449437975551e+02, new int[] { 0, 3 });
                p.AddCoeff(2.5555449437975551e+01, new int[] { 2, 1 });
                p.AddCoeff(-2.1722132022279218e+03, new int[] { 0, 5 });
                p.AddCoeff(-7.6666348313926652e+02, new int[] { 2, 3 });
                p.AddCoeff(7.8613430175867647e+03, new int[] { 0, 7 });
                p.AddCoeff(6.5166396066837655e+03, new int[] { 2, 5 });
                p.AddCoeff(-1.3757350280776838e+04, new int[] { 0, 9 });
                p.AddCoeff(-2.3584029052760294e+04, new int[] { 2, 7 });
                p.AddCoeff(1.1506147507558810e+04, new int[] { 0, 11 });
                p.AddCoeff(4.1272050842330515e+04, new int[] { 2, 9 });
                p.AddCoeff(-3.6878677908842340e+03, new int[] { 0, 13 });
                p.AddCoeff(-3.4518442522676430e+04, new int[] { 2, 11 });
                p.AddCoeff(1.1063603372652702e+04, new int[] { 2, 13 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{C2C3CB5A-A70D-49FD-BE30-56098902C091}"));
                OrthonormalPolynomials[123] = p;
                p.AddCoeff(-2.2381660871188492e+00, new int[] { 1, 0 });
                p.AddCoeff(1.7457695479527024e+02, new int[] { 1, 2 });
                p.AddCoeff(3.7302768118647487e+00, new int[] { 3, 0 });
                p.AddCoeff(-2.1822119349408780e+03, new int[] { 1, 4 });
                p.AddCoeff(-2.9096159132545040e+02, new int[] { 3, 2 });
                p.AddCoeff(9.8926941050653135e+03, new int[] { 1, 6 });
                p.AddCoeff(3.6370198915681300e+03, new int[] { 3, 4 });
                p.AddCoeff(-2.0138698713882960e+04, new int[] { 1, 8 });
                p.AddCoeff(-1.6487823508442189e+04, new int[] { 3, 6 });
                p.AddCoeff(1.8796118799624096e+04, new int[] { 1, 10 });
                p.AddCoeff(3.3564497856471600e+04, new int[] { 3, 8 });
                p.AddCoeff(-6.5501626119902152e+03, new int[] { 1, 12 });
                p.AddCoeff(-3.1326864666040160e+04, new int[] { 3, 10 });
                p.AddCoeff(1.0916937686650359e+04, new int[] { 3, 12 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{04579983-BC12-4033-BA17-0055E1FCD9CB}"));
                OrthonormalPolynomials[124] = p;
                p.AddCoeff(-7.3026370143802324e+00, new int[] { 0, 1 });
                p.AddCoeff(1.5822380197823837e+02, new int[] { 0, 3 });
                p.AddCoeff(7.3026370143802324e+01, new int[] { 2, 1 });
                p.AddCoeff(-9.4934281186943021e+02, new int[] { 0, 5 });
                p.AddCoeff(-1.5822380197823837e+03, new int[] { 2, 3 });
                p.AddCoeff(-8.5197431834436044e+01, new int[] { 4, 1 });
                p.AddCoeff(2.3055468288257591e+03, new int[] { 0, 7 });
                p.AddCoeff(9.4934281186943021e+03, new int[] { 2, 5 });
                p.AddCoeff(1.8459443564127810e+03, new int[] { 4, 3 });
                p.AddCoeff(-2.4336327637605235e+03, new int[] { 0, 9 });
                p.AddCoeff(-2.3055468288257591e+04, new int[] { 2, 7 });
                p.AddCoeff(-1.1075666138476686e+04, new int[] { 4, 5 });
                p.AddCoeff(9.2920523707219987e+02, new int[] { 0, 11 });
                p.AddCoeff(2.4336327637605235e+04, new int[] { 2, 9 });
                p.AddCoeff(2.6898046336300523e+04, new int[] { 4, 7 });
                p.AddCoeff(-9.2920523707219987e+03, new int[] { 2, 11 });
                p.AddCoeff(-2.8392382243872774e+04, new int[] { 4, 9 });
                p.AddCoeff(1.0840727765842332e+04, new int[] { 4, 11 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{332E60FF-ED2F-4B81-B585-1204931FEB6F}"));
                OrthonormalPolynomials[125] = p;
                p.AddCoeff(-3.5065323547666692e+00, new int[] { 1, 0 });
                p.AddCoeff(1.9285927951216681e+02, new int[] { 1, 2 });
                p.AddCoeff(1.6363817655577790e+01, new int[] { 3, 0 });
                p.AddCoeff(-1.6714470891054457e+03, new int[] { 1, 4 });
                p.AddCoeff(-9.0000997105677843e+02, new int[] { 3, 2 });
                p.AddCoeff(-1.4727435890020011e+01, new int[] { 5, 0 });
                p.AddCoeff(5.0143412673163370e+03, new int[] { 1, 6 });
                p.AddCoeff(7.8000864158254131e+03, new int[] { 3, 4 });
                p.AddCoeff(8.1000897395110059e+02, new int[] { 5, 2 });
                p.AddCoeff(-6.0888429674555521e+03, new int[] { 1, 8 });
                p.AddCoeff(-2.3400259247476239e+04, new int[] { 3, 6 });
                p.AddCoeff(-7.0200777742428718e+03, new int[] { 5, 4 });
                p.AddCoeff(2.5708448084812331e+03, new int[] { 1, 10 });
                p.AddCoeff(2.8414600514792576e+04, new int[] { 3, 8 });
                p.AddCoeff(2.1060233322728615e+04, new int[] { 5, 6 });
                p.AddCoeff(-1.1997275772912421e+04, new int[] { 3, 10 });
                p.AddCoeff(-2.5573140463313319e+04, new int[] { 5, 8 });
                p.AddCoeff(1.0797548195621179e+04, new int[] { 5, 10 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{894DB9C0-6561-4DA1-8549-9F6870AAB0C6}"));
                OrthonormalPolynomials[126] = p;
                p.AddCoeff(-6.0432294901526354e+00, new int[] { 0, 1 });
                p.AddCoeff(8.8634032522238653e+01, new int[] { 0, 3 });
                p.AddCoeff(1.2690781929320534e+02, new int[] { 2, 1 });
                p.AddCoeff(-3.4567272683673075e+02, new int[] { 0, 5 });
                p.AddCoeff(-1.8613146829670117e+03, new int[] { 2, 3 });
                p.AddCoeff(-3.8072345787961603e+02, new int[] { 4, 1 });
                p.AddCoeff(4.9381818119532964e+02, new int[] { 0, 7 });
                p.AddCoeff(7.2591272635713457e+03, new int[] { 2, 5 });
                p.AddCoeff(5.5839440489010352e+03, new int[] { 4, 3 });
                p.AddCoeff(2.7919720244505176e+02, new int[] { 6, 1 });
                p.AddCoeff(-2.3319191889779455e+02, new int[] { 0, 9 });
                p.AddCoeff(-1.0370181805101922e+04, new int[] { 2, 7 });
                p.AddCoeff(-2.1777381790714037e+04, new int[] { 4, 5 });
                p.AddCoeff(-4.0948923025274258e+03, new int[] { 6, 3 });
                p.AddCoeff(4.8970302968536856e+03, new int[] { 2, 9 });
                p.AddCoeff(3.1110545415305767e+04, new int[] { 4, 7 });
                p.AddCoeff(1.5970079979856961e+04, new int[] { 6, 5 });
                p.AddCoeff(-1.4691090890561057e+04, new int[] { 4, 9 });
                p.AddCoeff(-2.2814399971224229e+04, new int[] { 6, 7 });
                p.AddCoeff(1.0773466653078108e+04, new int[] { 6, 9 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{FB9DA9C9-3BF0-425A-806B-83B88248CD24}"));
                OrthonormalPolynomials[127] = p;
                p.AddCoeff(-4.7758010968682513e+00, new int[] { 1, 0 });
                p.AddCoeff(1.7192883948725705e+02, new int[] { 1, 2 });
                p.AddCoeff(4.2982209871814261e+01, new int[] { 3, 0 });
                p.AddCoeff(-9.4560861717991375e+02, new int[] { 1, 4 });
                p.AddCoeff(-1.5473595553853134e+03, new int[] { 3, 2 });
                p.AddCoeff(-9.4560861717991375e+01, new int[] { 5, 0 });
                p.AddCoeff(1.6390549364451838e+03, new int[] { 1, 6 });
                p.AddCoeff(8.5104775546192238e+03, new int[] { 3, 4 });
                p.AddCoeff(3.4041910218476895e+03, new int[] { 5, 2 });
                p.AddCoeff(5.8537676301613708e+01, new int[] { 7, 0 });
                p.AddCoeff(-8.7806514452420563e+02, new int[] { 1, 8 });
                p.AddCoeff(-1.4751494428006655e+04, new int[] { 3, 6 });
                p.AddCoeff(-1.8723050620162292e+04, new int[] { 5, 4 });
                p.AddCoeff(-2.1073563468580935e+03, new int[] { 7, 2 });
                p.AddCoeff(7.9025863007178506e+03, new int[] { 3, 8 });
                p.AddCoeff(3.2453287741614640e+04, new int[] { 5, 6 });
                p.AddCoeff(1.1590459907719514e+04, new int[] { 7, 4 });
                p.AddCoeff(-1.7385689861579271e+04, new int[] { 5, 8 });
                p.AddCoeff(-2.0090130506713825e+04, new int[] { 7, 6 });
                p.AddCoeff(1.0762569914310978e+04, new int[] { 7, 8 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{0DDFBF71-FC04-4BD5-924F-92DB214FA974}"));
                OrthonormalPolynomials[128] = p;
                p.AddCoeff(-4.7758010968682513e+00, new int[] { 0, 1 });
                p.AddCoeff(4.2982209871814261e+01, new int[] { 0, 3 });
                p.AddCoeff(1.7192883948725705e+02, new int[] { 2, 1 });
                p.AddCoeff(-9.4560861717991375e+01, new int[] { 0, 5 });
                p.AddCoeff(-1.5473595553853134e+03, new int[] { 2, 3 });
                p.AddCoeff(-9.4560861717991375e+02, new int[] { 4, 1 });
                p.AddCoeff(5.8537676301613708e+01, new int[] { 0, 7 });
                p.AddCoeff(3.4041910218476895e+03, new int[] { 2, 5 });
                p.AddCoeff(8.5104775546192238e+03, new int[] { 4, 3 });
                p.AddCoeff(1.6390549364451838e+03, new int[] { 6, 1 });
                p.AddCoeff(-2.1073563468580935e+03, new int[] { 2, 7 });
                p.AddCoeff(-1.8723050620162292e+04, new int[] { 4, 5 });
                p.AddCoeff(-1.4751494428006655e+04, new int[] { 6, 3 });
                p.AddCoeff(-8.7806514452420563e+02, new int[] { 8, 1 });
                p.AddCoeff(1.1590459907719514e+04, new int[] { 4, 7 });
                p.AddCoeff(3.2453287741614640e+04, new int[] { 6, 5 });
                p.AddCoeff(7.9025863007178506e+03, new int[] { 8, 3 });
                p.AddCoeff(-2.0090130506713825e+04, new int[] { 6, 7 });
                p.AddCoeff(-1.7385689861579271e+04, new int[] { 8, 5 });
                p.AddCoeff(1.0762569914310978e+04, new int[] { 8, 7 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{0315CF74-CC25-4660-BCC2-281CA2E28B73}"));
                OrthonormalPolynomials[129] = p;
                p.AddCoeff(-6.0432294901526354e+00, new int[] { 1, 0 });
                p.AddCoeff(1.2690781929320534e+02, new int[] { 1, 2 });
                p.AddCoeff(8.8634032522238653e+01, new int[] { 3, 0 });
                p.AddCoeff(-3.8072345787961603e+02, new int[] { 1, 4 });
                p.AddCoeff(-1.8613146829670117e+03, new int[] { 3, 2 });
                p.AddCoeff(-3.4567272683673075e+02, new int[] { 5, 0 });
                p.AddCoeff(2.7919720244505176e+02, new int[] { 1, 6 });
                p.AddCoeff(5.5839440489010352e+03, new int[] { 3, 4 });
                p.AddCoeff(7.2591272635713457e+03, new int[] { 5, 2 });
                p.AddCoeff(4.9381818119532964e+02, new int[] { 7, 0 });
                p.AddCoeff(-4.0948923025274258e+03, new int[] { 3, 6 });
                p.AddCoeff(-2.1777381790714037e+04, new int[] { 5, 4 });
                p.AddCoeff(-1.0370181805101922e+04, new int[] { 7, 2 });
                p.AddCoeff(-2.3319191889779455e+02, new int[] { 9, 0 });
                p.AddCoeff(1.5970079979856961e+04, new int[] { 5, 6 });
                p.AddCoeff(3.1110545415305767e+04, new int[] { 7, 4 });
                p.AddCoeff(4.8970302968536856e+03, new int[] { 9, 2 });
                p.AddCoeff(-2.2814399971224229e+04, new int[] { 7, 6 });
                p.AddCoeff(-1.4691090890561057e+04, new int[] { 9, 4 });
                p.AddCoeff(1.0773466653078108e+04, new int[] { 9, 6 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{3761744B-508A-4AB1-BF08-7A403F52FE2A}"));
                OrthonormalPolynomials[130] = p;
                p.AddCoeff(-3.5065323547666692e+00, new int[] { 0, 1 });
                p.AddCoeff(1.6363817655577790e+01, new int[] { 0, 3 });
                p.AddCoeff(1.9285927951216681e+02, new int[] { 2, 1 });
                p.AddCoeff(-1.4727435890020011e+01, new int[] { 0, 5 });
                p.AddCoeff(-9.0000997105677843e+02, new int[] { 2, 3 });
                p.AddCoeff(-1.6714470891054457e+03, new int[] { 4, 1 });
                p.AddCoeff(8.1000897395110059e+02, new int[] { 2, 5 });
                p.AddCoeff(7.8000864158254131e+03, new int[] { 4, 3 });
                p.AddCoeff(5.0143412673163370e+03, new int[] { 6, 1 });
                p.AddCoeff(-7.0200777742428718e+03, new int[] { 4, 5 });
                p.AddCoeff(-2.3400259247476239e+04, new int[] { 6, 3 });
                p.AddCoeff(-6.0888429674555521e+03, new int[] { 8, 1 });
                p.AddCoeff(2.1060233322728615e+04, new int[] { 6, 5 });
                p.AddCoeff(2.8414600514792576e+04, new int[] { 8, 3 });
                p.AddCoeff(2.5708448084812331e+03, new int[] { 10, 1 });
                p.AddCoeff(-2.5573140463313319e+04, new int[] { 8, 5 });
                p.AddCoeff(-1.1997275772912421e+04, new int[] { 10, 3 });
                p.AddCoeff(1.0797548195621179e+04, new int[] { 10, 5 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{FBB3CE8C-B7BF-4688-BBC7-C787C4CEA1B3}"));
                OrthonormalPolynomials[131] = p;
                p.AddCoeff(-7.3026370143802324e+00, new int[] { 1, 0 });
                p.AddCoeff(7.3026370143802324e+01, new int[] { 1, 2 });
                p.AddCoeff(1.5822380197823837e+02, new int[] { 3, 0 });
                p.AddCoeff(-8.5197431834436044e+01, new int[] { 1, 4 });
                p.AddCoeff(-1.5822380197823837e+03, new int[] { 3, 2 });
                p.AddCoeff(-9.4934281186943021e+02, new int[] { 5, 0 });
                p.AddCoeff(1.8459443564127810e+03, new int[] { 3, 4 });
                p.AddCoeff(9.4934281186943021e+03, new int[] { 5, 2 });
                p.AddCoeff(2.3055468288257591e+03, new int[] { 7, 0 });
                p.AddCoeff(-1.1075666138476686e+04, new int[] { 5, 4 });
                p.AddCoeff(-2.3055468288257591e+04, new int[] { 7, 2 });
                p.AddCoeff(-2.4336327637605235e+03, new int[] { 9, 0 });
                p.AddCoeff(2.6898046336300523e+04, new int[] { 7, 4 });
                p.AddCoeff(2.4336327637605235e+04, new int[] { 9, 2 });
                p.AddCoeff(9.2920523707219987e+02, new int[] { 11, 0 });
                p.AddCoeff(-2.8392382243872774e+04, new int[] { 9, 4 });
                p.AddCoeff(-9.2920523707219987e+03, new int[] { 11, 2 });
                p.AddCoeff(1.0840727765842332e+04, new int[] { 11, 4 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{0D3D3DA6-7291-49B4-AA2C-0EE6AFE0FAD5}"));
                OrthonormalPolynomials[132] = p;
                p.AddCoeff(-2.2381660871188492e+00, new int[] { 0, 1 });
                p.AddCoeff(3.7302768118647487e+00, new int[] { 0, 3 });
                p.AddCoeff(1.7457695479527024e+02, new int[] { 2, 1 });
                p.AddCoeff(-2.9096159132545040e+02, new int[] { 2, 3 });
                p.AddCoeff(-2.1822119349408780e+03, new int[] { 4, 1 });
                p.AddCoeff(3.6370198915681300e+03, new int[] { 4, 3 });
                p.AddCoeff(9.8926941050653135e+03, new int[] { 6, 1 });
                p.AddCoeff(-1.6487823508442189e+04, new int[] { 6, 3 });
                p.AddCoeff(-2.0138698713882960e+04, new int[] { 8, 1 });
                p.AddCoeff(3.3564497856471600e+04, new int[] { 8, 3 });
                p.AddCoeff(1.8796118799624096e+04, new int[] { 10, 1 });
                p.AddCoeff(-3.1326864666040160e+04, new int[] { 10, 3 });
                p.AddCoeff(-6.5501626119902152e+03, new int[] { 12, 1 });
                p.AddCoeff(1.0916937686650359e+04, new int[] { 12, 3 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{3ADB7D66-24C3-4C74-9713-CBE07B97B846}"));
                OrthonormalPolynomials[133] = p;
                p.AddCoeff(-8.5184831459918503e+00, new int[] { 1, 0 });
                p.AddCoeff(2.5555449437975551e+01, new int[] { 1, 2 });
                p.AddCoeff(2.5555449437975551e+02, new int[] { 3, 0 });
                p.AddCoeff(-7.6666348313926652e+02, new int[] { 3, 2 });
                p.AddCoeff(-2.1722132022279218e+03, new int[] { 5, 0 });
                p.AddCoeff(6.5166396066837655e+03, new int[] { 5, 2 });
                p.AddCoeff(7.8613430175867647e+03, new int[] { 7, 0 });
                p.AddCoeff(-2.3584029052760294e+04, new int[] { 7, 2 });
                p.AddCoeff(-1.3757350280776838e+04, new int[] { 9, 0 });
                p.AddCoeff(4.1272050842330515e+04, new int[] { 9, 2 });
                p.AddCoeff(1.1506147507558810e+04, new int[] { 11, 0 });
                p.AddCoeff(-3.4518442522676430e+04, new int[] { 11, 2 });
                p.AddCoeff(-3.6878677908842340e+03, new int[] { 13, 0 });
                p.AddCoeff(1.1063603372652702e+04, new int[] { 13, 2 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{C8094705-8C86-427F-974F-D8862D38DAB4}"));
                OrthonormalPolynomials[134] = p;
                p.AddCoeff(-9.7691543305056193e-01, new int[] { 0, 1 });
                p.AddCoeff(1.0257612047030900e+02, new int[] { 2, 1 });
                p.AddCoeff(-1.7437940479952530e+03, new int[] { 4, 1 });
                p.AddCoeff(1.1044028970636603e+04, new int[] { 6, 1 });
                p.AddCoeff(-3.3132086911909808e+04, new int[] { 8, 1 });
                p.AddCoeff(5.0802533264928372e+04, new int[] { 10, 1 });
                p.AddCoeff(-3.8486767624945736e+04, new int[] { 12, 1 });
                p.AddCoeff(1.1419150833775109e+04, new int[] { 14, 1 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{2DB7F033-D8B8-4749-AED2-CA4B729799ED}"));
                OrthonormalPolynomials[135] = p;
                p.AddCoeff(-8.7472079284207009e+00, new int[] { 1, 0 });
                p.AddCoeff(3.4697258116068780e+02, new int[] { 3, 0 });
                p.AddCoeff(-3.9554874252318410e+03, new int[] { 5, 0 });
                p.AddCoeff(1.9777437126159205e+04, new int[] { 7, 0 });
                p.AddCoeff(-5.0542339322406857e+04, new int[] { 9, 0 });
                p.AddCoeff(6.8921371803282077e+04, new int[] { 11, 0 });
                p.AddCoeff(-4.7714795863810669e+04, new int[] { 13, 0 });
                p.AddCoeff(1.3178372190957232e+04, new int[] { 15, 0 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{DB093DB6-D827-4487-9F27-92EF1FEB68B8}"));
                OrthonormalPolynomials[136] = p;
                p.AddCoeff(5.6406037338977378e-01, new int[] { 0, 0 });
                p.AddCoeff(-7.6712210781009234e+01, new int[] { 0, 2 });
                p.AddCoeff(1.7004540056457047e+03, new int[] { 0, 4 });
                p.AddCoeff(-1.4283813647423919e+04, new int[] { 0, 6 });
                p.AddCoeff(5.8665663194776812e+04, new int[] { 0, 8 });
                p.AddCoeff(-1.3036814043283736e+05, new int[] { 0, 10 });
                p.AddCoeff(1.5999726325848221e+05, new int[] { 0, 12 });
                p.AddCoeff(-1.0197627768123042e+05, new int[] { 0, 14 });
                p.AddCoeff(2.6343871734317859e+04, new int[] { 0, 16 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{7F4D5BDC-2D81-4361-8F29-505E861CE37B}"));
                OrthonormalPolynomials[137] = p;
                p.AddCoeff(-1.5150608556393961e+01, new int[] { 1, 1 });
                p.AddCoeff(6.0097413940362713e+02, new int[] { 1, 3 });
                p.AddCoeff(-6.8511051892013493e+03, new int[] { 1, 5 });
                p.AddCoeff(3.4255525946006746e+04, new int[] { 1, 7 });
                p.AddCoeff(-8.7541899639795019e+04, new int[] { 1, 9 });
                p.AddCoeff(1.1937531769062957e+05, new int[] { 1, 11 });
                p.AddCoeff(-8.2644450708897395e+04, new int[] { 1, 13 });
                p.AddCoeff(2.2825610195790709e+04, new int[] { 1, 15 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{FC849285-60FA-4D64-9D2B-EC03A487CB46}"));
                OrthonormalPolynomials[138] = p;
                p.AddCoeff(6.3059620047630551e-01, new int[] { 0, 0 });
                p.AddCoeff(-6.6212601050012079e+01, new int[] { 0, 2 });
                p.AddCoeff(-1.8917886014289165e+00, new int[] { 2, 0 });
                p.AddCoeff(1.1256142178502053e+03, new int[] { 0, 4 });
                p.AddCoeff(1.9863780315003624e+02, new int[] { 2, 2 });
                p.AddCoeff(-7.1288900463846338e+03, new int[] { 0, 6 });
                p.AddCoeff(-3.3768426535506160e+03, new int[] { 2, 4 });
                p.AddCoeff(2.1386670139153901e+04, new int[] { 0, 8 });
                p.AddCoeff(2.1386670139153901e+04, new int[] { 2, 6 });
                p.AddCoeff(-3.2792894213369316e+04, new int[] { 0, 10 });
                p.AddCoeff(-6.4160010417461704e+04, new int[] { 2, 8 });
                p.AddCoeff(2.4843101676794936e+04, new int[] { 0, 12 });
                p.AddCoeff(9.8378682640107947e+04, new int[] { 2, 10 });
                p.AddCoeff(-7.3710301678402558e+03, new int[] { 0, 14 });
                p.AddCoeff(-7.4529305030384808e+04, new int[] { 2, 12 });
                p.AddCoeff(2.2113090503520767e+04, new int[] { 2, 14 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{F06DF6E1-0339-4BDB-B8AE-4C4B8FA98EC2}"));
                OrthonormalPolynomials[139] = p;
                p.AddCoeff(-3.0237615553606320e+01, new int[] { 1, 1 });
                p.AddCoeff(9.0712846660818961e+02, new int[] { 1, 3 });
                p.AddCoeff(5.0396025922677200e+01, new int[] { 3, 1 });
                p.AddCoeff(-7.7105919661696117e+03, new int[] { 1, 5 });
                p.AddCoeff(-1.5118807776803160e+03, new int[] { 3, 3 });
                p.AddCoeff(2.7904999496613833e+04, new int[] { 1, 7 });
                p.AddCoeff(1.2850986610282686e+04, new int[] { 3, 5 });
                p.AddCoeff(-4.8833749119074207e+04, new int[] { 1, 9 });
                p.AddCoeff(-4.6508332494356388e+04, new int[] { 3, 7 });
                p.AddCoeff(4.0842771990498428e+04, new int[] { 1, 11 });
                p.AddCoeff(8.1389581865123679e+04, new int[] { 3, 9 });
                p.AddCoeff(-1.3090632048236676e+04, new int[] { 1, 13 });
                p.AddCoeff(-6.8071286650830713e+04, new int[] { 3, 11 });
                p.AddCoeff(2.1817720080394459e+04, new int[] { 3, 13 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{2F6BD435-8171-4514-A13C-E1356D1E1647}"));
                OrthonormalPolynomials[140] = p;
                p.AddCoeff(6.3446044921875000e-01, new int[] { 0, 0 });
                p.AddCoeff(-4.9487915039062500e+01, new int[] { 0, 2 });
                p.AddCoeff(-6.3446044921875000e+00, new int[] { 2, 0 });
                p.AddCoeff(6.1859893798828125e+02, new int[] { 0, 4 });
                p.AddCoeff(4.9487915039062500e+02, new int[] { 2, 2 });
                p.AddCoeff(7.4020385742187500e+00, new int[] { 4, 0 });
                p.AddCoeff(-2.8043151855468750e+03, new int[] { 0, 6 });
                p.AddCoeff(-6.1859893798828125e+03, new int[] { 2, 4 });
                p.AddCoeff(-5.7735900878906250e+02, new int[] { 4, 2 });
                p.AddCoeff(5.7087844848632812e+03, new int[] { 0, 8 });
                p.AddCoeff(2.8043151855468750e+04, new int[] { 2, 6 });
                p.AddCoeff(7.2169876098632812e+03, new int[] { 4, 4 });
                p.AddCoeff(-5.3281988525390625e+03, new int[] { 0, 10 });
                p.AddCoeff(-5.7087844848632812e+04, new int[] { 2, 8 });
                p.AddCoeff(-3.2717010498046875e+04, new int[] { 4, 6 });
                p.AddCoeff(1.8567965698242188e+03, new int[] { 0, 12 });
                p.AddCoeff(5.3281988525390625e+04, new int[] { 2, 10 });
                p.AddCoeff(6.6602485656738281e+04, new int[] { 4, 8 });
                p.AddCoeff(-1.8567965698242188e+04, new int[] { 2, 12 });
                p.AddCoeff(-6.2162319946289062e+04, new int[] { 4, 10 });
                p.AddCoeff(2.1662626647949219e+04, new int[] { 4, 12 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{62569AA1-33A0-4C07-8BC9-0825BEFDED04}"));
                OrthonormalPolynomials[141] = p;
                p.AddCoeff(-4.0366844928100702e+01, new int[] { 1, 1 });
                p.AddCoeff(8.7461497344218188e+02, new int[] { 1, 3 });
                p.AddCoeff(1.8837860966446994e+02, new int[] { 3, 1 });
                p.AddCoeff(-5.2476898406530913e+03, new int[] { 1, 5 });
                p.AddCoeff(-4.0815365427301821e+03, new int[] { 3, 3 });
                p.AddCoeff(-1.6954074869802295e+02, new int[] { 5, 1 });
                p.AddCoeff(1.2744389613014650e+04, new int[] { 1, 7 });
                p.AddCoeff(2.4489219256381093e+04, new int[] { 3, 5 });
                p.AddCoeff(3.6733828884571639e+03, new int[] { 5, 3 });
                p.AddCoeff(-1.3452411258182131e+04, new int[] { 1, 9 });
                p.AddCoeff(-5.9473818194068368e+04, new int[] { 3, 7 });
                p.AddCoeff(-2.2040297330742983e+04, new int[] { 5, 5 });
                p.AddCoeff(5.1363752076695409e+03, new int[] { 1, 11 });
                p.AddCoeff(6.2777919204849944e+04, new int[] { 3, 9 });
                p.AddCoeff(5.3526436374661531e+04, new int[] { 5, 7 });
                p.AddCoeff(-2.3969750969124524e+04, new int[] { 3, 11 });
                p.AddCoeff(-5.6500127284364949e+04, new int[] { 5, 9 });
                p.AddCoeff(2.1572775872212072e+04, new int[] { 5, 11 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{B69D391D-B0C7-4B2F-8895-4B25E70F54B1}"));
                OrthonormalPolynomials[142] = p;
                p.AddCoeff(6.3533376064274492e-01, new int[] { 0, 0 });
                p.AddCoeff(-3.4943356835350971e+01, new int[] { 0, 2 });
                p.AddCoeff(-1.3342008973497643e+01, new int[] { 2, 0 });
                p.AddCoeff(3.0284242590637508e+02, new int[] { 0, 4 });
                p.AddCoeff(7.3381049354237039e+02, new int[] { 2, 2 });
                p.AddCoeff(4.0026026920492930e+01, new int[] { 4, 0 });
                p.AddCoeff(-9.0852727771912524e+02, new int[] { 0, 6 });
                p.AddCoeff(-6.3596909440338767e+03, new int[] { 2, 4 });
                p.AddCoeff(-2.2014314806271112e+03, new int[] { 4, 2 });
                p.AddCoeff(-2.9352419741694815e+01, new int[] { 6, 0 });
                p.AddCoeff(1.1032116943732235e+03, new int[] { 0, 8 });
                p.AddCoeff(1.9079072832101630e+04, new int[] { 2, 6 });
                p.AddCoeff(1.9079072832101630e+04, new int[] { 4, 4 });
                p.AddCoeff(1.6143830857932148e+03, new int[] { 6, 2 });
                p.AddCoeff(-4.6580049317980548e+02, new int[] { 0, 10 });
                p.AddCoeff(-2.3167445581837694e+04, new int[] { 2, 8 });
                p.AddCoeff(-5.7237218496304890e+04, new int[] { 4, 6 });
                p.AddCoeff(-1.3991320076874529e+04, new int[] { 6, 4 });
                p.AddCoeff(9.7818103567759151e+03, new int[] { 2, 10 });
                p.AddCoeff(6.9502336745513081e+04, new int[] { 4, 8 });
                p.AddCoeff(4.1973960230623586e+04, new int[] { 6, 6 });
                p.AddCoeff(-2.9345431070327745e+04, new int[] { 4, 10 });
                p.AddCoeff(-5.0968380280042926e+04, new int[] { 6, 8 });
                p.AddCoeff(2.1519982784907013e+04, new int[] { 6, 10 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{7ABBDE06-D123-41BE-9D57-3619ACBA6D7B}"));
                OrthonormalPolynomials[143] = p;
                p.AddCoeff(-4.5440288513886428e+01, new int[] { 1, 1 });
                p.AddCoeff(6.6645756487033427e+02, new int[] { 1, 3 });
                p.AddCoeff(4.0896259662497785e+02, new int[] { 3, 1 });
                p.AddCoeff(-2.5991845029943037e+03, new int[] { 1, 5 });
                p.AddCoeff(-5.9981180838330084e+03, new int[] { 3, 3 });
                p.AddCoeff(-8.9971771257495127e+02, new int[] { 5, 1 });
                p.AddCoeff(3.7131207185632909e+03, new int[] { 1, 7 });
                p.AddCoeff(2.3392660526948733e+04, new int[] { 3, 5 });
                p.AddCoeff(1.3195859784432619e+04, new int[] { 5, 3 });
                p.AddCoeff(5.5696810778449364e+02, new int[] { 7, 1 });
                p.AddCoeff(-1.7534181170993318e+03, new int[] { 1, 9 });
                p.AddCoeff(-3.3418086467069618e+04, new int[] { 3, 7 });
                p.AddCoeff(-5.1463853159287212e+04, new int[] { 5, 5 });
                p.AddCoeff(-8.1688655808392401e+03, new int[] { 7, 3 });
                p.AddCoeff(1.5780763053893986e+04, new int[] { 3, 9 });
                p.AddCoeff(7.3519790227553161e+04, new int[] { 5, 7 });
                p.AddCoeff(3.1858575765273036e+04, new int[] { 7, 5 });
                p.AddCoeff(-3.4717678718566770e+04, new int[] { 5, 9 });
                p.AddCoeff(-4.5512251093247195e+04, new int[] { 7, 7 });
                p.AddCoeff(2.1491896349588953e+04, new int[] { 7, 9 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{30D7E352-8E6F-464D-81F2-ECA273319F6F}"));
                OrthonormalPolynomials[144] = p;
                p.AddCoeff(6.3552856445312500e-01, new int[] { 0, 0 });
                p.AddCoeff(-2.2879028320312500e+01, new int[] { 0, 2 });
                p.AddCoeff(-2.2879028320312500e+01, new int[] { 2, 0 });
                p.AddCoeff(1.2583465576171875e+02, new int[] { 0, 4 });
                p.AddCoeff(8.2364501953125000e+02, new int[] { 2, 2 });
                p.AddCoeff(1.2583465576171875e+02, new int[] { 4, 0 });
                p.AddCoeff(-2.1811340332031250e+02, new int[] { 0, 6 });
                p.AddCoeff(-4.5300476074218750e+03, new int[] { 2, 4 });
                p.AddCoeff(-4.5300476074218750e+03, new int[] { 4, 2 });
                p.AddCoeff(-2.1811340332031250e+02, new int[] { 6, 0 });
                p.AddCoeff(1.1684646606445312e+02, new int[] { 0, 8 });
                p.AddCoeff(7.8520825195312500e+03, new int[] { 2, 6 });
                p.AddCoeff(2.4915261840820312e+04, new int[] { 4, 4 });
                p.AddCoeff(7.8520825195312500e+03, new int[] { 6, 2 });
                p.AddCoeff(1.1684646606445312e+02, new int[] { 8, 0 });
                p.AddCoeff(-4.2064727783203125e+03, new int[] { 2, 8 });
                p.AddCoeff(-4.3186453857421875e+04, new int[] { 4, 6 });
                p.AddCoeff(-4.3186453857421875e+04, new int[] { 6, 4 });
                p.AddCoeff(-4.2064727783203125e+03, new int[] { 8, 2 });
                p.AddCoeff(2.3135600280761719e+04, new int[] { 4, 8 });
                p.AddCoeff(7.4856520019531250e+04, new int[] { 6, 6 });
                p.AddCoeff(2.3135600280761719e+04, new int[] { 8, 4 });
                p.AddCoeff(-4.0101707153320312e+04, new int[] { 6, 8 });
                p.AddCoeff(-4.0101707153320312e+04, new int[] { 8, 6 });
                p.AddCoeff(2.1483057403564453e+04, new int[] { 8, 8 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{23DF02AC-7487-49B0-A45B-07135F884FDE}"));
                OrthonormalPolynomials[145] = p;
                p.AddCoeff(-4.5440288513886428e+01, new int[] { 1, 1 });
                p.AddCoeff(4.0896259662497785e+02, new int[] { 1, 3 });
                p.AddCoeff(6.6645756487033427e+02, new int[] { 3, 1 });
                p.AddCoeff(-8.9971771257495127e+02, new int[] { 1, 5 });
                p.AddCoeff(-5.9981180838330084e+03, new int[] { 3, 3 });
                p.AddCoeff(-2.5991845029943037e+03, new int[] { 5, 1 });
                p.AddCoeff(5.5696810778449364e+02, new int[] { 1, 7 });
                p.AddCoeff(1.3195859784432619e+04, new int[] { 3, 5 });
                p.AddCoeff(2.3392660526948733e+04, new int[] { 5, 3 });
                p.AddCoeff(3.7131207185632909e+03, new int[] { 7, 1 });
                p.AddCoeff(-8.1688655808392401e+03, new int[] { 3, 7 });
                p.AddCoeff(-5.1463853159287212e+04, new int[] { 5, 5 });
                p.AddCoeff(-3.3418086467069618e+04, new int[] { 7, 3 });
                p.AddCoeff(-1.7534181170993318e+03, new int[] { 9, 1 });
                p.AddCoeff(3.1858575765273036e+04, new int[] { 5, 7 });
                p.AddCoeff(7.3519790227553161e+04, new int[] { 7, 5 });
                p.AddCoeff(1.5780763053893986e+04, new int[] { 9, 3 });
                p.AddCoeff(-4.5512251093247195e+04, new int[] { 7, 7 });
                p.AddCoeff(-3.4717678718566770e+04, new int[] { 9, 5 });
                p.AddCoeff(2.1491896349588953e+04, new int[] { 9, 7 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{729D0768-341A-4E24-BC31-8F340C25E9F4}"));
                OrthonormalPolynomials[146] = p;
                p.AddCoeff(6.3533376064274492e-01, new int[] { 0, 0 });
                p.AddCoeff(-1.3342008973497643e+01, new int[] { 0, 2 });
                p.AddCoeff(-3.4943356835350971e+01, new int[] { 2, 0 });
                p.AddCoeff(4.0026026920492930e+01, new int[] { 0, 4 });
                p.AddCoeff(7.3381049354237039e+02, new int[] { 2, 2 });
                p.AddCoeff(3.0284242590637508e+02, new int[] { 4, 0 });
                p.AddCoeff(-2.9352419741694815e+01, new int[] { 0, 6 });
                p.AddCoeff(-2.2014314806271112e+03, new int[] { 2, 4 });
                p.AddCoeff(-6.3596909440338767e+03, new int[] { 4, 2 });
                p.AddCoeff(-9.0852727771912524e+02, new int[] { 6, 0 });
                p.AddCoeff(1.6143830857932148e+03, new int[] { 2, 6 });
                p.AddCoeff(1.9079072832101630e+04, new int[] { 4, 4 });
                p.AddCoeff(1.9079072832101630e+04, new int[] { 6, 2 });
                p.AddCoeff(1.1032116943732235e+03, new int[] { 8, 0 });
                p.AddCoeff(-1.3991320076874529e+04, new int[] { 4, 6 });
                p.AddCoeff(-5.7237218496304890e+04, new int[] { 6, 4 });
                p.AddCoeff(-2.3167445581837694e+04, new int[] { 8, 2 });
                p.AddCoeff(-4.6580049317980548e+02, new int[] { 10, 0 });
                p.AddCoeff(4.1973960230623586e+04, new int[] { 6, 6 });
                p.AddCoeff(6.9502336745513081e+04, new int[] { 8, 4 });
                p.AddCoeff(9.7818103567759151e+03, new int[] { 10, 2 });
                p.AddCoeff(-5.0968380280042926e+04, new int[] { 8, 6 });
                p.AddCoeff(-2.9345431070327745e+04, new int[] { 10, 4 });
                p.AddCoeff(2.1519982784907013e+04, new int[] { 10, 6 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{C195A9E6-BB9F-451F-948B-E08B3F100681}"));
                OrthonormalPolynomials[147] = p;
                p.AddCoeff(-4.0366844928100702e+01, new int[] { 1, 1 });
                p.AddCoeff(1.8837860966446994e+02, new int[] { 1, 3 });
                p.AddCoeff(8.7461497344218188e+02, new int[] { 3, 1 });
                p.AddCoeff(-1.6954074869802295e+02, new int[] { 1, 5 });
                p.AddCoeff(-4.0815365427301821e+03, new int[] { 3, 3 });
                p.AddCoeff(-5.2476898406530913e+03, new int[] { 5, 1 });
                p.AddCoeff(3.6733828884571639e+03, new int[] { 3, 5 });
                p.AddCoeff(2.4489219256381093e+04, new int[] { 5, 3 });
                p.AddCoeff(1.2744389613014650e+04, new int[] { 7, 1 });
                p.AddCoeff(-2.2040297330742983e+04, new int[] { 5, 5 });
                p.AddCoeff(-5.9473818194068368e+04, new int[] { 7, 3 });
                p.AddCoeff(-1.3452411258182131e+04, new int[] { 9, 1 });
                p.AddCoeff(5.3526436374661531e+04, new int[] { 7, 5 });
                p.AddCoeff(6.2777919204849944e+04, new int[] { 9, 3 });
                p.AddCoeff(5.1363752076695409e+03, new int[] { 11, 1 });
                p.AddCoeff(-5.6500127284364949e+04, new int[] { 9, 5 });
                p.AddCoeff(-2.3969750969124524e+04, new int[] { 11, 3 });
                p.AddCoeff(2.1572775872212072e+04, new int[] { 11, 5 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{3C191B18-05E5-4839-A855-9D060F156703}"));
                OrthonormalPolynomials[148] = p;
                p.AddCoeff(6.3446044921875000e-01, new int[] { 0, 0 });
                p.AddCoeff(-6.3446044921875000e+00, new int[] { 0, 2 });
                p.AddCoeff(-4.9487915039062500e+01, new int[] { 2, 0 });
                p.AddCoeff(7.4020385742187500e+00, new int[] { 0, 4 });
                p.AddCoeff(4.9487915039062500e+02, new int[] { 2, 2 });
                p.AddCoeff(6.1859893798828125e+02, new int[] { 4, 0 });
                p.AddCoeff(-5.7735900878906250e+02, new int[] { 2, 4 });
                p.AddCoeff(-6.1859893798828125e+03, new int[] { 4, 2 });
                p.AddCoeff(-2.8043151855468750e+03, new int[] { 6, 0 });
                p.AddCoeff(7.2169876098632812e+03, new int[] { 4, 4 });
                p.AddCoeff(2.8043151855468750e+04, new int[] { 6, 2 });
                p.AddCoeff(5.7087844848632812e+03, new int[] { 8, 0 });
                p.AddCoeff(-3.2717010498046875e+04, new int[] { 6, 4 });
                p.AddCoeff(-5.7087844848632812e+04, new int[] { 8, 2 });
                p.AddCoeff(-5.3281988525390625e+03, new int[] { 10, 0 });
                p.AddCoeff(6.6602485656738281e+04, new int[] { 8, 4 });
                p.AddCoeff(5.3281988525390625e+04, new int[] { 10, 2 });
                p.AddCoeff(1.8567965698242188e+03, new int[] { 12, 0 });
                p.AddCoeff(-6.2162319946289062e+04, new int[] { 10, 4 });
                p.AddCoeff(-1.8567965698242188e+04, new int[] { 12, 2 });
                p.AddCoeff(2.1662626647949219e+04, new int[] { 12, 4 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{F320E4FF-8101-4FA4-9F19-7FFDE690A30F}"));
                OrthonormalPolynomials[149] = p;
                p.AddCoeff(-3.0237615553606320e+01, new int[] { 1, 1 });
                p.AddCoeff(5.0396025922677200e+01, new int[] { 1, 3 });
                p.AddCoeff(9.0712846660818961e+02, new int[] { 3, 1 });
                p.AddCoeff(-1.5118807776803160e+03, new int[] { 3, 3 });
                p.AddCoeff(-7.7105919661696117e+03, new int[] { 5, 1 });
                p.AddCoeff(1.2850986610282686e+04, new int[] { 5, 3 });
                p.AddCoeff(2.7904999496613833e+04, new int[] { 7, 1 });
                p.AddCoeff(-4.6508332494356388e+04, new int[] { 7, 3 });
                p.AddCoeff(-4.8833749119074207e+04, new int[] { 9, 1 });
                p.AddCoeff(8.1389581865123679e+04, new int[] { 9, 3 });
                p.AddCoeff(4.0842771990498428e+04, new int[] { 11, 1 });
                p.AddCoeff(-6.8071286650830713e+04, new int[] { 11, 3 });
                p.AddCoeff(-1.3090632048236676e+04, new int[] { 13, 1 });
                p.AddCoeff(2.1817720080394459e+04, new int[] { 13, 3 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{1574C3D4-970C-4B45-94DB-3E680FF12C25}"));
                OrthonormalPolynomials[150] = p;
                p.AddCoeff(6.3059620047630551e-01, new int[] { 0, 0 });
                p.AddCoeff(-1.8917886014289165e+00, new int[] { 0, 2 });
                p.AddCoeff(-6.6212601050012079e+01, new int[] { 2, 0 });
                p.AddCoeff(1.9863780315003624e+02, new int[] { 2, 2 });
                p.AddCoeff(1.1256142178502053e+03, new int[] { 4, 0 });
                p.AddCoeff(-3.3768426535506160e+03, new int[] { 4, 2 });
                p.AddCoeff(-7.1288900463846338e+03, new int[] { 6, 0 });
                p.AddCoeff(2.1386670139153901e+04, new int[] { 6, 2 });
                p.AddCoeff(2.1386670139153901e+04, new int[] { 8, 0 });
                p.AddCoeff(-6.4160010417461704e+04, new int[] { 8, 2 });
                p.AddCoeff(-3.2792894213369316e+04, new int[] { 10, 0 });
                p.AddCoeff(9.8378682640107947e+04, new int[] { 10, 2 });
                p.AddCoeff(2.4843101676794936e+04, new int[] { 12, 0 });
                p.AddCoeff(-7.4529305030384808e+04, new int[] { 12, 2 });
                p.AddCoeff(-7.3710301678402558e+03, new int[] { 14, 0 });
                p.AddCoeff(2.2113090503520767e+04, new int[] { 14, 2 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{1E8C8D8E-4016-4AC1-AFA3-C7DCCE4D4155}"));
                OrthonormalPolynomials[151] = p;
                p.AddCoeff(-1.5150608556393961e+01, new int[] { 1, 1 });
                p.AddCoeff(6.0097413940362713e+02, new int[] { 3, 1 });
                p.AddCoeff(-6.8511051892013493e+03, new int[] { 5, 1 });
                p.AddCoeff(3.4255525946006746e+04, new int[] { 7, 1 });
                p.AddCoeff(-8.7541899639795019e+04, new int[] { 9, 1 });
                p.AddCoeff(1.1937531769062957e+05, new int[] { 11, 1 });
                p.AddCoeff(-8.2644450708897395e+04, new int[] { 13, 1 });
                p.AddCoeff(2.2825610195790709e+04, new int[] { 15, 1 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{34AF1AB9-CDC5-42A0-90FC-1E4919CA8542}"));
                OrthonormalPolynomials[152] = p;
                p.AddCoeff(5.6406037338977378e-01, new int[] { 0, 0 });
                p.AddCoeff(-7.6712210781009234e+01, new int[] { 2, 0 });
                p.AddCoeff(1.7004540056457047e+03, new int[] { 4, 0 });
                p.AddCoeff(-1.4283813647423919e+04, new int[] { 6, 0 });
                p.AddCoeff(5.8665663194776812e+04, new int[] { 8, 0 });
                p.AddCoeff(-1.3036814043283736e+05, new int[] { 10, 0 });
                p.AddCoeff(1.5999726325848221e+05, new int[] { 12, 0 });
                p.AddCoeff(-1.0197627768123042e+05, new int[] { 14, 0 });
                p.AddCoeff(2.6343871734317859e+04, new int[] { 16, 0 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{A58BB949-04E0-4097-B400-5436EB967306}"));
                OrthonormalPolynomials[153] = p;
                p.AddCoeff(9.8753287944363784e+00, new int[] { 0, 1 });
                p.AddCoeff(-5.0034999225144317e+02, new int[] { 0, 3 });
                p.AddCoeff(7.3551448860962147e+03, new int[] { 0, 5 });
                p.AddCoeff(-4.8333809251489411e+04, new int[] { 0, 7 });
                p.AddCoeff(1.6782572656767156e+05, new int[] { 0, 9 });
                p.AddCoeff(-3.2954869944197325e+05, new int[] { 0, 11 });
                p.AddCoeff(3.6757354937758555e+05, new int[] { 0, 13 });
                p.AddCoeff(-2.1704342915628861e+05, new int[] { 0, 15 });
                p.AddCoeff(5.2664949721746501e+04, new int[] { 0, 17 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{E04DCE39-FE0D-48E3-9CB1-6A393D0D2A70}"));
                OrthonormalPolynomials[154] = p;
                p.AddCoeff(9.7698122524736014e-01, new int[] { 1, 0 });
                p.AddCoeff(-1.3286944663364098e+02, new int[] { 1, 2 });
                p.AddCoeff(2.9452727337123750e+03, new int[] { 1, 4 });
                p.AddCoeff(-2.4740290963183950e+04, new int[] { 1, 6 });
                p.AddCoeff(1.0161190931307694e+05, new int[] { 1, 8 });
                p.AddCoeff(-2.2580424291794875e+05, new int[] { 1, 10 });
                p.AddCoeff(2.7712338903566438e+05, new int[] { 1, 12 });
                p.AddCoeff(-1.7662809411064323e+05, new int[] { 1, 14 });
                p.AddCoeff(4.5628924311916168e+04, new int[] { 1, 16 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{2F2BC13B-791E-44B7-9D01-3A28E7658C9D}"));
                OrthonormalPolynomials[155] = p;
                p.AddCoeff(9.7796757706369010e+00, new int[] { 0, 1 });
                p.AddCoeff(-3.8792713890193041e+02, new int[] { 0, 3 });
                p.AddCoeff(-2.9339027311910703e+01, new int[] { 2, 1 });
                p.AddCoeff(4.4223693834820066e+03, new int[] { 0, 5 });
                p.AddCoeff(1.1637814167057912e+03, new int[] { 2, 3 });
                p.AddCoeff(-2.2111846917410033e+04, new int[] { 0, 7 });
                p.AddCoeff(-1.3267108150446020e+04, new int[] { 2, 5 });
                p.AddCoeff(5.6508053233381196e+04, new int[] { 0, 9 });
                p.AddCoeff(6.6335540752230099e+04, new int[] { 2, 7 });
                p.AddCoeff(-7.7056436227337994e+04, new int[] { 0, 11 });
                p.AddCoeff(-1.6952415970014359e+05, new int[] { 2, 9 });
                p.AddCoeff(5.3346763542003227e+04, new int[] { 0, 13 });
                p.AddCoeff(2.3116930868201398e+05, new int[] { 2, 11 });
                p.AddCoeff(-1.4733868025886605e+04, new int[] { 0, 15 });
                p.AddCoeff(-1.6004029062600968e+05, new int[] { 2, 13 });
                p.AddCoeff(4.4201604077659816e+04, new int[] { 2, 15 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{29BED423-3EE8-4253-9F23-C1F25E5DAEED}"));
                OrthonormalPolynomials[156] = p;
                p.AddCoeff(2.2383944597623821e+00, new int[] { 1, 0 });
                p.AddCoeff(-2.3503141827505012e+02, new int[] { 1, 2 });
                p.AddCoeff(-3.7306574329373035e+00, new int[] { 3, 0 });
                p.AddCoeff(3.9955341106758521e+03, new int[] { 1, 4 });
                p.AddCoeff(3.9171903045841687e+02, new int[] { 3, 2 });
                p.AddCoeff(-2.5305049367613730e+04, new int[] { 1, 6 });
                p.AddCoeff(-6.6592235177930868e+03, new int[] { 3, 4 });
                p.AddCoeff(7.5915148102841189e+04, new int[] { 1, 8 });
                p.AddCoeff(4.2175082279356216e+04, new int[] { 3, 6 });
                p.AddCoeff(-1.1640322709102316e+05, new int[] { 1, 10 });
                p.AddCoeff(-1.2652524683806865e+05, new int[] { 3, 8 });
                p.AddCoeff(8.8184262947744816e+04, new int[] { 1, 12 });
                p.AddCoeff(1.9400537848503859e+05, new int[] { 3, 10 });
                p.AddCoeff(-2.6164561533946264e+04, new int[] { 1, 14 });
                p.AddCoeff(-1.4697377157957469e+05, new int[] { 3, 12 });
                p.AddCoeff(4.3607602556577107e+04, new int[] { 3, 14 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{3D852180-C38B-4600-AC82-0C0097720DB5}"));
                OrthonormalPolynomials[157] = p;
                p.AddCoeff(8.5715583208308191e+00, new int[] { 0, 1 });
                p.AddCoeff(-2.5714674962492457e+02, new int[] { 0, 3 });
                p.AddCoeff(-8.5715583208308191e+01, new int[] { 2, 1 });
                p.AddCoeff(2.1857473718118589e+03, new int[] { 0, 5 });
                p.AddCoeff(2.5714674962492457e+03, new int[] { 2, 3 });
                p.AddCoeff(1.0000151374302622e+02, new int[] { 4, 1 });
                p.AddCoeff(-7.9103238217952988e+03, new int[] { 0, 7 });
                p.AddCoeff(-2.1857473718118589e+04, new int[] { 2, 5 });
                p.AddCoeff(-3.0000454122907867e+03, new int[] { 4, 3 });
                p.AddCoeff(1.3843066688141773e+04, new int[] { 0, 9 });
                p.AddCoeff(7.9103238217952988e+04, new int[] { 2, 7 });
                p.AddCoeff(2.5500386004471687e+04, new int[] { 4, 5 });
                p.AddCoeff(-1.1577837593718574e+04, new int[] { 0, 11 });
                p.AddCoeff(-1.3843066688141773e+05, new int[] { 2, 9 });
                p.AddCoeff(-9.2287111254278485e+04, new int[] { 4, 7 });
                p.AddCoeff(3.7108453826021069e+03, new int[] { 0, 13 });
                p.AddCoeff(1.1577837593718574e+05, new int[] { 2, 11 });
                p.AddCoeff(1.6150244469498735e+05, new int[] { 4, 9 });
                p.AddCoeff(-3.7108453826021069e+04, new int[] { 2, 13 });
                p.AddCoeff(-1.3507477192671669e+05, new int[] { 4, 11 });
                p.AddCoeff(4.3293196130357914e+04, new int[] { 4, 13 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{47B93A7B-13F0-496E-A0EC-EC9DAF40BEEB}"));
                OrthonormalPolynomials[158] = p;
                p.AddCoeff(3.5071120906315492e+00, new int[] { 1, 0 });
                p.AddCoeff(-2.7355474306926084e+02, new int[] { 1, 2 });
                p.AddCoeff(-1.6366523089613896e+01, new int[] { 3, 0 });
                p.AddCoeff(3.4194342883657605e+03, new int[] { 1, 4 });
                p.AddCoeff(1.2765888009898839e+03, new int[] { 3, 2 });
                p.AddCoeff(1.4729870780652507e+01, new int[] { 5, 0 });
                p.AddCoeff(-1.5501435440591448e+04, new int[] { 1, 6 });
                p.AddCoeff(-1.5957360012373549e+04, new int[] { 3, 4 });
                p.AddCoeff(-1.1489299208908955e+03, new int[] { 5, 2 });
                p.AddCoeff(3.1556493575489732e+04, new int[] { 1, 8 });
                p.AddCoeff(7.2340032056093422e+04, new int[] { 3, 6 });
                p.AddCoeff(1.4361624011136194e+04, new int[] { 5, 4 });
                p.AddCoeff(-2.9452727337123750e+04, new int[] { 1, 10 });
                p.AddCoeff(-1.4726363668561875e+05, new int[] { 3, 8 });
                p.AddCoeff(-6.5106028850484080e+04, new int[] { 5, 6 });
                p.AddCoeff(1.0263829223543125e+04, new int[] { 1, 12 });
                p.AddCoeff(1.3744606090657750e+05, new int[] { 3, 10 });
                p.AddCoeff(1.3253727301705688e+05, new int[] { 5, 8 });
                p.AddCoeff(-4.7897869709867917e+04, new int[] { 3, 12 });
                p.AddCoeff(-1.2370145481591975e+05, new int[] { 5, 10 });
                p.AddCoeff(4.3108082738881125e+04, new int[] { 5, 12 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{0A8BAE86-E881-419E-9DA1-8F8BD021E4F8}"));
                OrthonormalPolynomials[159] = p;
                p.AddCoeff(7.3138978337358849e+00, new int[] { 0, 1 });
                p.AddCoeff(-1.5846778639761084e+02, new int[] { 0, 3 });
                p.AddCoeff(-1.5359185450845358e+02, new int[] { 2, 1 });
                p.AddCoeff(9.5080671838566503e+02, new int[] { 0, 5 });
                p.AddCoeff(3.3278235143498276e+03, new int[] { 2, 3 });
                p.AddCoeff(4.6077556352536075e+02, new int[] { 4, 1 });
                p.AddCoeff(-2.3091020303651865e+03, new int[] { 0, 7 });
                p.AddCoeff(-1.9966941086098966e+04, new int[] { 2, 5 });
                p.AddCoeff(-9.9834705430494828e+03, new int[] { 4, 3 });
                p.AddCoeff(-3.3790207991859788e+02, new int[] { 6, 1 });
                p.AddCoeff(2.4373854764965858e+03, new int[] { 0, 9 });
                p.AddCoeff(4.8491142637668917e+04, new int[] { 2, 7 });
                p.AddCoeff(5.9900823258296897e+04, new int[] { 4, 5 });
                p.AddCoeff(7.3212117315696208e+03, new int[] { 6, 3 });
                p.AddCoeff(-9.3063809102596911e+02, new int[] { 0, 11 });
                p.AddCoeff(-5.1185095006428301e+04, new int[] { 2, 9 });
                p.AddCoeff(-1.4547342791300675e+05, new int[] { 4, 7 });
                p.AddCoeff(-4.3927270389417725e+04, new int[] { 6, 5 });
                p.AddCoeff(1.9543399911545351e+04, new int[] { 2, 11 });
                p.AddCoeff(1.5355528501928490e+05, new int[] { 4, 9 });
                p.AddCoeff(1.0668051380287162e+05, new int[] { 6, 7 });
                p.AddCoeff(-5.8630199734636054e+04, new int[] { 4, 11 });
                p.AddCoeff(-1.1260720901414226e+05, new int[] { 6, 9 });
                p.AddCoeff(4.2995479805399773e+04, new int[] { 6, 11 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{0CDEEC43-528D-450C-8928-336AA67E1E70}"));
                OrthonormalPolynomials[160] = p;
                p.AddCoeff(4.7772055377446240e+00, new int[] { 1, 0 });
                p.AddCoeff(-2.6274630457595432e+02, new int[] { 1, 2 });
                p.AddCoeff(-4.2994849839701616e+01, new int[] { 3, 0 });
                p.AddCoeff(2.2771346396582708e+03, new int[] { 1, 4 });
                p.AddCoeff(2.3647167411835889e+03, new int[] { 3, 2 });
                p.AddCoeff(9.4588669647343556e+01, new int[] { 5, 0 });
                p.AddCoeff(-6.8314039189748124e+03, new int[] { 1, 6 });
                p.AddCoeff(-2.0494211756924437e+04, new int[] { 3, 4 });
                p.AddCoeff(-5.2023768306038956e+03, new int[] { 5, 2 });
                p.AddCoeff(-5.8554890734069820e+01, new int[] { 7, 0 });
                p.AddCoeff(8.2952761873265579e+03, new int[] { 1, 8 });
                p.AddCoeff(6.1482635270773311e+04, new int[] { 3, 6 });
                p.AddCoeff(4.5087265865233762e+04, new int[] { 5, 4 });
                p.AddCoeff(3.2205189903738401e+03, new int[] { 7, 2 });
                p.AddCoeff(-3.5024499457601022e+03, new int[] { 1, 10 });
                p.AddCoeff(-7.4657485685939021e+04, new int[] { 3, 8 });
                p.AddCoeff(-1.3526179759570129e+05, new int[] { 5, 6 });
                p.AddCoeff(-2.7911164583239948e+04, new int[] { 7, 4 });
                p.AddCoeff(3.1522049511840920e+04, new int[] { 3, 10 });
                p.AddCoeff(1.6424646850906585e+05, new int[] { 5, 8 });
                p.AddCoeff(8.3733493749719843e+04, new int[] { 7, 6 });
                p.AddCoeff(-6.9348508926050024e+04, new int[] { 5, 10 });
                p.AddCoeff(-1.0167638526751695e+05, new int[] { 7, 8 });
                p.AddCoeff(4.2930029335173824e+04, new int[] { 7, 10 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{E5C06B9C-BAFC-4208-A1C1-AC66A79262B4}"));
                OrthonormalPolynomials[161] = p;
                p.AddCoeff(6.0468601480290527e+00, new int[] { 0, 1 });
                p.AddCoeff(-8.8687282171092773e+01, new int[] { 0, 3 });
                p.AddCoeff(-2.1768696532904590e+02, new int[] { 2, 1 });
                p.AddCoeff(3.4588040046726181e+02, new int[] { 0, 5 });
                p.AddCoeff(3.1927421581593398e+03, new int[] { 2, 3 });
                p.AddCoeff(1.1972783093097524e+03, new int[] { 4, 1 });
                p.AddCoeff(-4.9411485781037402e+02, new int[] { 0, 7 });
                p.AddCoeff(-1.2451694416821425e+04, new int[] { 2, 5 });
                p.AddCoeff(-1.7560081869876369e+04, new int[] { 4, 3 });
                p.AddCoeff(-2.0752824028035709e+03, new int[] { 6, 1 });
                p.AddCoeff(2.3333201618823218e+02, new int[] { 0, 9 });
                p.AddCoeff(1.7788134881173465e+04, new int[] { 2, 7 });
                p.AddCoeff(6.8484319292517839e+04, new int[] { 4, 5 });
                p.AddCoeff(3.0437475241119040e+04, new int[] { 6, 3 });
                p.AddCoeff(1.1117584300733415e+03, new int[] { 8, 1 });
                p.AddCoeff(-8.3999525827763583e+03, new int[] { 2, 9 });
                p.AddCoeff(-9.7834741846454056e+04, new int[] { 4, 7 });
                p.AddCoeff(-1.1870615344036425e+05, new int[] { 6, 5 });
                p.AddCoeff(-1.6305790307742343e+04, new int[] { 8, 3 });
                p.AddCoeff(4.6199739205269971e+04, new int[] { 4, 9 });
                p.AddCoeff(1.6958021920052036e+05, new int[] { 6, 7 });
                p.AddCoeff(6.3592582200195136e+04, new int[] { 8, 5 });
                p.AddCoeff(-8.0079547955801283e+04, new int[] { 6, 9 });
                p.AddCoeff(-9.0846546000278766e+04, new int[] { 8, 7 });
                p.AddCoeff(4.2899757833464973e+04, new int[] { 8, 9 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{D6942D2C-0D77-4F11-8229-80BE9F0FF75D}"));
                OrthonormalPolynomials[162] = p;
                p.AddCoeff(6.0468601480290527e+00, new int[] { 1, 0 });
                p.AddCoeff(-2.1768696532904590e+02, new int[] { 1, 2 });
                p.AddCoeff(-8.8687282171092773e+01, new int[] { 3, 0 });
                p.AddCoeff(1.1972783093097524e+03, new int[] { 1, 4 });
                p.AddCoeff(3.1927421581593398e+03, new int[] { 3, 2 });
                p.AddCoeff(3.4588040046726181e+02, new int[] { 5, 0 });
                p.AddCoeff(-2.0752824028035709e+03, new int[] { 1, 6 });
                p.AddCoeff(-1.7560081869876369e+04, new int[] { 3, 4 });
                p.AddCoeff(-1.2451694416821425e+04, new int[] { 5, 2 });
                p.AddCoeff(-4.9411485781037402e+02, new int[] { 7, 0 });
                p.AddCoeff(1.1117584300733415e+03, new int[] { 1, 8 });
                p.AddCoeff(3.0437475241119040e+04, new int[] { 3, 6 });
                p.AddCoeff(6.8484319292517839e+04, new int[] { 5, 4 });
                p.AddCoeff(1.7788134881173465e+04, new int[] { 7, 2 });
                p.AddCoeff(2.3333201618823218e+02, new int[] { 9, 0 });
                p.AddCoeff(-1.6305790307742343e+04, new int[] { 3, 8 });
                p.AddCoeff(-1.1870615344036425e+05, new int[] { 5, 6 });
                p.AddCoeff(-9.7834741846454056e+04, new int[] { 7, 4 });
                p.AddCoeff(-8.3999525827763583e+03, new int[] { 9, 2 });
                p.AddCoeff(6.3592582200195136e+04, new int[] { 5, 8 });
                p.AddCoeff(1.6958021920052036e+05, new int[] { 7, 6 });
                p.AddCoeff(4.6199739205269971e+04, new int[] { 9, 4 });
                p.AddCoeff(-9.0846546000278766e+04, new int[] { 7, 8 });
                p.AddCoeff(-8.0079547955801283e+04, new int[] { 9, 6 });
                p.AddCoeff(4.2899757833464973e+04, new int[] { 9, 8 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{75DCE2F8-1A84-4D5E-93EB-D215DF844621}"));
                OrthonormalPolynomials[163] = p;
                p.AddCoeff(4.7772055377446240e+00, new int[] { 0, 1 });
                p.AddCoeff(-4.2994849839701616e+01, new int[] { 0, 3 });
                p.AddCoeff(-2.6274630457595432e+02, new int[] { 2, 1 });
                p.AddCoeff(9.4588669647343556e+01, new int[] { 0, 5 });
                p.AddCoeff(2.3647167411835889e+03, new int[] { 2, 3 });
                p.AddCoeff(2.2771346396582708e+03, new int[] { 4, 1 });
                p.AddCoeff(-5.8554890734069820e+01, new int[] { 0, 7 });
                p.AddCoeff(-5.2023768306038956e+03, new int[] { 2, 5 });
                p.AddCoeff(-2.0494211756924437e+04, new int[] { 4, 3 });
                p.AddCoeff(-6.8314039189748124e+03, new int[] { 6, 1 });
                p.AddCoeff(3.2205189903738401e+03, new int[] { 2, 7 });
                p.AddCoeff(4.5087265865233762e+04, new int[] { 4, 5 });
                p.AddCoeff(6.1482635270773311e+04, new int[] { 6, 3 });
                p.AddCoeff(8.2952761873265579e+03, new int[] { 8, 1 });
                p.AddCoeff(-2.7911164583239948e+04, new int[] { 4, 7 });
                p.AddCoeff(-1.3526179759570129e+05, new int[] { 6, 5 });
                p.AddCoeff(-7.4657485685939021e+04, new int[] { 8, 3 });
                p.AddCoeff(-3.5024499457601022e+03, new int[] { 10, 1 });
                p.AddCoeff(8.3733493749719843e+04, new int[] { 6, 7 });
                p.AddCoeff(1.6424646850906585e+05, new int[] { 8, 5 });
                p.AddCoeff(3.1522049511840920e+04, new int[] { 10, 3 });
                p.AddCoeff(-1.0167638526751695e+05, new int[] { 8, 7 });
                p.AddCoeff(-6.9348508926050024e+04, new int[] { 10, 5 });
                p.AddCoeff(4.2930029335173824e+04, new int[] { 10, 7 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{45198E35-E108-41B4-88FA-788C66E9EEC1}"));
                OrthonormalPolynomials[164] = p;
                p.AddCoeff(7.3138978337358849e+00, new int[] { 1, 0 });
                p.AddCoeff(-1.5359185450845358e+02, new int[] { 1, 2 });
                p.AddCoeff(-1.5846778639761084e+02, new int[] { 3, 0 });
                p.AddCoeff(4.6077556352536075e+02, new int[] { 1, 4 });
                p.AddCoeff(3.3278235143498276e+03, new int[] { 3, 2 });
                p.AddCoeff(9.5080671838566503e+02, new int[] { 5, 0 });
                p.AddCoeff(-3.3790207991859788e+02, new int[] { 1, 6 });
                p.AddCoeff(-9.9834705430494828e+03, new int[] { 3, 4 });
                p.AddCoeff(-1.9966941086098966e+04, new int[] { 5, 2 });
                p.AddCoeff(-2.3091020303651865e+03, new int[] { 7, 0 });
                p.AddCoeff(7.3212117315696208e+03, new int[] { 3, 6 });
                p.AddCoeff(5.9900823258296897e+04, new int[] { 5, 4 });
                p.AddCoeff(4.8491142637668917e+04, new int[] { 7, 2 });
                p.AddCoeff(2.4373854764965858e+03, new int[] { 9, 0 });
                p.AddCoeff(-4.3927270389417725e+04, new int[] { 5, 6 });
                p.AddCoeff(-1.4547342791300675e+05, new int[] { 7, 4 });
                p.AddCoeff(-5.1185095006428301e+04, new int[] { 9, 2 });
                p.AddCoeff(-9.3063809102596911e+02, new int[] { 11, 0 });
                p.AddCoeff(1.0668051380287162e+05, new int[] { 7, 6 });
                p.AddCoeff(1.5355528501928490e+05, new int[] { 9, 4 });
                p.AddCoeff(1.9543399911545351e+04, new int[] { 11, 2 });
                p.AddCoeff(-1.1260720901414226e+05, new int[] { 9, 6 });
                p.AddCoeff(-5.8630199734636054e+04, new int[] { 11, 4 });
                p.AddCoeff(4.2995479805399773e+04, new int[] { 11, 6 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{0AA0A244-C8AD-4E2B-85C1-09FA6B34BBEB}"));
                OrthonormalPolynomials[165] = p;
                p.AddCoeff(3.5071120906315492e+00, new int[] { 0, 1 });
                p.AddCoeff(-1.6366523089613896e+01, new int[] { 0, 3 });
                p.AddCoeff(-2.7355474306926084e+02, new int[] { 2, 1 });
                p.AddCoeff(1.4729870780652507e+01, new int[] { 0, 5 });
                p.AddCoeff(1.2765888009898839e+03, new int[] { 2, 3 });
                p.AddCoeff(3.4194342883657605e+03, new int[] { 4, 1 });
                p.AddCoeff(-1.1489299208908955e+03, new int[] { 2, 5 });
                p.AddCoeff(-1.5957360012373549e+04, new int[] { 4, 3 });
                p.AddCoeff(-1.5501435440591448e+04, new int[] { 6, 1 });
                p.AddCoeff(1.4361624011136194e+04, new int[] { 4, 5 });
                p.AddCoeff(7.2340032056093422e+04, new int[] { 6, 3 });
                p.AddCoeff(3.1556493575489732e+04, new int[] { 8, 1 });
                p.AddCoeff(-6.5106028850484080e+04, new int[] { 6, 5 });
                p.AddCoeff(-1.4726363668561875e+05, new int[] { 8, 3 });
                p.AddCoeff(-2.9452727337123750e+04, new int[] { 10, 1 });
                p.AddCoeff(1.3253727301705688e+05, new int[] { 8, 5 });
                p.AddCoeff(1.3744606090657750e+05, new int[] { 10, 3 });
                p.AddCoeff(1.0263829223543125e+04, new int[] { 12, 1 });
                p.AddCoeff(-1.2370145481591975e+05, new int[] { 10, 5 });
                p.AddCoeff(-4.7897869709867917e+04, new int[] { 12, 3 });
                p.AddCoeff(4.3108082738881125e+04, new int[] { 12, 5 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{61D7540D-8E1E-4B7A-9489-D6BE61126637}"));
                OrthonormalPolynomials[166] = p;
                p.AddCoeff(8.5715583208308191e+00, new int[] { 1, 0 });
                p.AddCoeff(-8.5715583208308191e+01, new int[] { 1, 2 });
                p.AddCoeff(-2.5714674962492457e+02, new int[] { 3, 0 });
                p.AddCoeff(1.0000151374302622e+02, new int[] { 1, 4 });
                p.AddCoeff(2.5714674962492457e+03, new int[] { 3, 2 });
                p.AddCoeff(2.1857473718118589e+03, new int[] { 5, 0 });
                p.AddCoeff(-3.0000454122907867e+03, new int[] { 3, 4 });
                p.AddCoeff(-2.1857473718118589e+04, new int[] { 5, 2 });
                p.AddCoeff(-7.9103238217952988e+03, new int[] { 7, 0 });
                p.AddCoeff(2.5500386004471687e+04, new int[] { 5, 4 });
                p.AddCoeff(7.9103238217952988e+04, new int[] { 7, 2 });
                p.AddCoeff(1.3843066688141773e+04, new int[] { 9, 0 });
                p.AddCoeff(-9.2287111254278485e+04, new int[] { 7, 4 });
                p.AddCoeff(-1.3843066688141773e+05, new int[] { 9, 2 });
                p.AddCoeff(-1.1577837593718574e+04, new int[] { 11, 0 });
                p.AddCoeff(1.6150244469498735e+05, new int[] { 9, 4 });
                p.AddCoeff(1.1577837593718574e+05, new int[] { 11, 2 });
                p.AddCoeff(3.7108453826021069e+03, new int[] { 13, 0 });
                p.AddCoeff(-1.3507477192671669e+05, new int[] { 11, 4 });
                p.AddCoeff(-3.7108453826021069e+04, new int[] { 13, 2 });
                p.AddCoeff(4.3293196130357914e+04, new int[] { 13, 4 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{1E879529-DAEB-4FFA-AB6D-3FF8BDE398DD}"));
                OrthonormalPolynomials[167] = p;
                p.AddCoeff(2.2383944597623821e+00, new int[] { 0, 1 });
                p.AddCoeff(-3.7306574329373035e+00, new int[] { 0, 3 });
                p.AddCoeff(-2.3503141827505012e+02, new int[] { 2, 1 });
                p.AddCoeff(3.9171903045841687e+02, new int[] { 2, 3 });
                p.AddCoeff(3.9955341106758521e+03, new int[] { 4, 1 });
                p.AddCoeff(-6.6592235177930868e+03, new int[] { 4, 3 });
                p.AddCoeff(-2.5305049367613730e+04, new int[] { 6, 1 });
                p.AddCoeff(4.2175082279356216e+04, new int[] { 6, 3 });
                p.AddCoeff(7.5915148102841189e+04, new int[] { 8, 1 });
                p.AddCoeff(-1.2652524683806865e+05, new int[] { 8, 3 });
                p.AddCoeff(-1.1640322709102316e+05, new int[] { 10, 1 });
                p.AddCoeff(1.9400537848503859e+05, new int[] { 10, 3 });
                p.AddCoeff(8.8184262947744816e+04, new int[] { 12, 1 });
                p.AddCoeff(-1.4697377157957469e+05, new int[] { 12, 3 });
                p.AddCoeff(-2.6164561533946264e+04, new int[] { 14, 1 });
                p.AddCoeff(4.3607602556577107e+04, new int[] { 14, 3 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{59E81E2E-35AD-4049-931E-3F17589155B1}"));
                OrthonormalPolynomials[168] = p;
                p.AddCoeff(9.7796757706369010e+00, new int[] { 1, 0 });
                p.AddCoeff(-2.9339027311910703e+01, new int[] { 1, 2 });
                p.AddCoeff(-3.8792713890193041e+02, new int[] { 3, 0 });
                p.AddCoeff(1.1637814167057912e+03, new int[] { 3, 2 });
                p.AddCoeff(4.4223693834820066e+03, new int[] { 5, 0 });
                p.AddCoeff(-1.3267108150446020e+04, new int[] { 5, 2 });
                p.AddCoeff(-2.2111846917410033e+04, new int[] { 7, 0 });
                p.AddCoeff(6.6335540752230099e+04, new int[] { 7, 2 });
                p.AddCoeff(5.6508053233381196e+04, new int[] { 9, 0 });
                p.AddCoeff(-1.6952415970014359e+05, new int[] { 9, 2 });
                p.AddCoeff(-7.7056436227337994e+04, new int[] { 11, 0 });
                p.AddCoeff(2.3116930868201398e+05, new int[] { 11, 2 });
                p.AddCoeff(5.3346763542003227e+04, new int[] { 13, 0 });
                p.AddCoeff(-1.6004029062600968e+05, new int[] { 13, 2 });
                p.AddCoeff(-1.4733868025886605e+04, new int[] { 15, 0 });
                p.AddCoeff(4.4201604077659816e+04, new int[] { 15, 2 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{A9090946-BB14-4381-9A7C-3B5341D3F453}"));
                OrthonormalPolynomials[169] = p;
                p.AddCoeff(9.7698122524736014e-01, new int[] { 0, 1 });
                p.AddCoeff(-1.3286944663364098e+02, new int[] { 2, 1 });
                p.AddCoeff(2.9452727337123750e+03, new int[] { 4, 1 });
                p.AddCoeff(-2.4740290963183950e+04, new int[] { 6, 1 });
                p.AddCoeff(1.0161190931307694e+05, new int[] { 8, 1 });
                p.AddCoeff(-2.2580424291794875e+05, new int[] { 10, 1 });
                p.AddCoeff(2.7712338903566438e+05, new int[] { 12, 1 });
                p.AddCoeff(-1.7662809411064323e+05, new int[] { 14, 1 });
                p.AddCoeff(4.5628924311916168e+04, new int[] { 16, 1 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{31B1496C-1742-4253-BB7B-47E1589C646E}"));
                OrthonormalPolynomials[170] = p;
                p.AddCoeff(9.8753287944363784e+00, new int[] { 1, 0 });
                p.AddCoeff(-5.0034999225144317e+02, new int[] { 3, 0 });
                p.AddCoeff(7.3551448860962147e+03, new int[] { 5, 0 });
                p.AddCoeff(-4.8333809251489411e+04, new int[] { 7, 0 });
                p.AddCoeff(1.6782572656767156e+05, new int[] { 9, 0 });
                p.AddCoeff(-3.2954869944197325e+05, new int[] { 11, 0 });
                p.AddCoeff(3.6757354937758555e+05, new int[] { 13, 0 });
                p.AddCoeff(-2.1704342915628861e+05, new int[] { 15, 0 });
                p.AddCoeff(5.2664949721746501e+04, new int[] { 17, 0 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{F768B235-15E8-4616-A3EB-2F577ECEC55C}"));
                OrthonormalPolynomials[171] = p;
                p.AddCoeff(-5.6408675045604599e-01, new int[] { 0, 0 });
                p.AddCoeff(9.6458834327983865e+01, new int[] { 0, 2 });
                p.AddCoeff(-2.7008473611835482e+03, new int[] { 0, 4 });
                p.AddCoeff(2.8989095010036751e+04, new int[] { 0, 6 });
                p.AddCoeff(-1.5529872326805402e+05, new int[] { 0, 8 });
                p.AddCoeff(4.6589616980416207e+05, new int[] { 0, 10 });
                p.AddCoeff(-8.1884781359519394e+05, new int[] { 0, 12 });
                p.AddCoeff(8.3684446883904435e+05, new int[] { 0, 14 });
                p.AddCoeff(-4.6026445786147439e+05, new int[] { 0, 16 });
                p.AddCoeff(1.0528925506635035e+05, new int[] { 0, 18 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{0E0C43CF-467D-4511-91EF-F4EA84357580}"));
                OrthonormalPolynomials[172] = p;
                p.AddCoeff(1.7104571213411717e+01, new int[] { 1, 1 });
                p.AddCoeff(-8.6663160814619365e+02, new int[] { 1, 3 });
                p.AddCoeff(1.2739484639749047e+04, new int[] { 1, 5 });
                p.AddCoeff(-8.3716613346922306e+04, new int[] { 1, 7 });
                p.AddCoeff(2.9068268523236912e+05, new int[] { 1, 9 });
                p.AddCoeff(-5.7079509100174300e+05, new int[] { 1, 11 });
                p.AddCoeff(6.3665606304040565e+05, new int[] { 1, 13 });
                p.AddCoeff(-3.7593024674766810e+05, new int[] { 1, 15 });
                p.AddCoeff(9.1218368696125347e+04, new int[] { 1, 17 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{1F129701-9C4B-440E-8AF2-2DAD54F380E9}"));
                OrthonormalPolynomials[173] = p;
                p.AddCoeff(-6.3063866915672383e-01, new int[] { 0, 0 });
                p.AddCoeff(8.5766859005314440e+01, new int[] { 0, 2 });
                p.AddCoeff(1.8919160074701715e+00, new int[] { 2, 0 });
                p.AddCoeff(-1.9011653746178034e+03, new int[] { 0, 4 });
                p.AddCoeff(-2.5730057701594332e+02, new int[] { 2, 2 });
                p.AddCoeff(1.5969789146789549e+04, new int[] { 0, 6 });
                p.AddCoeff(5.7034961238534103e+03, new int[] { 2, 4 });
                p.AddCoeff(-6.5590205424314218e+04, new int[] { 0, 8 });
                p.AddCoeff(-4.7909367440368646e+04, new int[] { 2, 6 });
                p.AddCoeff(1.4575601205403160e+05, new int[] { 0, 10 });
                p.AddCoeff(1.9677061627294265e+05, new int[] { 2, 8 });
                p.AddCoeff(-1.7888237842994787e+05, new int[] { 0, 12 });
                p.AddCoeff(-4.3726803616209479e+05, new int[] { 2, 10 });
                p.AddCoeff(1.1401294449381293e+05, new int[] { 0, 14 });
                p.AddCoeff(5.3664713528984360e+05, new int[] { 2, 12 });
                p.AddCoeff(-2.9453343994235006e+04, new int[] { 0, 16 });
                p.AddCoeff(-3.4203883348143878e+05, new int[] { 2, 14 });
                p.AddCoeff(8.8360031982705018e+04, new int[] { 2, 16 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{2B0765F9-938A-48DE-AB3C-1C28C751DE1E}"));
                OrthonormalPolynomials[174] = p;
                p.AddCoeff(3.4714405267160477e+01, new int[] { 1, 1 });
                p.AddCoeff(-1.3770047422640322e+03, new int[] { 1, 3 });
                p.AddCoeff(-5.7857342111934128e+01, new int[] { 3, 1 });
                p.AddCoeff(1.5697854061809967e+04, new int[] { 1, 5 });
                p.AddCoeff(2.2950079037733871e+03, new int[] { 3, 3 });
                p.AddCoeff(-7.8489270309049837e+04, new int[] { 1, 7 });
                p.AddCoeff(-2.6163090103016612e+04, new int[] { 3, 5 });
                p.AddCoeff(2.0058369078979403e+05, new int[] { 1, 9 });
                p.AddCoeff(1.3081545051508306e+05, new int[] { 3, 7 });
                p.AddCoeff(-2.7352321471335549e+05, new int[] { 1, 11 });
                p.AddCoeff(-3.3430615131632338e+05, new int[] { 3, 9 });
                p.AddCoeff(1.8936222557078457e+05, new int[] { 1, 13 });
                p.AddCoeff(4.5587202452225916e+05, new int[] { 3, 11 });
                p.AddCoeff(-5.2300043252883358e+04, new int[] { 1, 15 });
                p.AddCoeff(-3.1560370928464095e+05, new int[] { 3, 13 });
                p.AddCoeff(8.7166738754805597e+04, new int[] { 3, 15 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{97954A97-C507-4424-A69C-A1A0229B79F1}"));
                OrthonormalPolynomials[175] = p;
                p.AddCoeff(-6.3452518677814697e-01, new int[] { 0, 0 });
                p.AddCoeff(6.6625144611705432e+01, new int[] { 0, 2 });
                p.AddCoeff(6.3452518677814697e+00, new int[] { 2, 0 });
                p.AddCoeff(-1.1326274583989923e+03, new int[] { 0, 4 });
                p.AddCoeff(-6.6625144611705432e+02, new int[] { 2, 2 });
                p.AddCoeff(-7.4027938457450480e+00, new int[] { 4, 0 });
                p.AddCoeff(7.1733072365269515e+03, new int[] { 0, 6 });
                p.AddCoeff(1.1326274583989923e+04, new int[] { 2, 4 });
                p.AddCoeff(7.7729335380323004e+02, new int[] { 4, 2 });
                p.AddCoeff(-2.1519921709580855e+04, new int[] { 0, 8 });
                p.AddCoeff(-7.1733072365269515e+04, new int[] { 2, 6 });
                p.AddCoeff(-1.3213987014654911e+04, new int[] { 4, 4 });
                p.AddCoeff(3.2997213288023977e+04, new int[] { 0, 10 });
                p.AddCoeff(2.1519921709580855e+05, new int[] { 2, 8 });
                p.AddCoeff(8.3688584426147768e+04, new int[] { 4, 6 });
                p.AddCoeff(-2.4997888854563619e+04, new int[] { 0, 12 });
                p.AddCoeff(-3.2997213288023977e+05, new int[] { 2, 10 });
                p.AddCoeff(-2.5106575327844330e+05, new int[] { 4, 8 });
                p.AddCoeff(7.4169560337716232e+03, new int[] { 0, 14 });
                p.AddCoeff(2.4997888854563619e+05, new int[] { 2, 12 });
                p.AddCoeff(3.8496748836027973e+05, new int[] { 4, 10 });
                p.AddCoeff(-7.4169560337716232e+04, new int[] { 2, 14 });
                p.AddCoeff(-2.9164203663657556e+05, new int[] { 4, 12 });
                p.AddCoeff(8.6531153727335604e+04, new int[] { 4, 14 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{9F550CC1-89E1-459D-8B6D-F2AA3AAD818B}"));
                OrthonormalPolynomials[176] = p;
                p.AddCoeff(4.7381071364740997e+01, new int[] { 1, 1 });
                p.AddCoeff(-1.4214321409422299e+03, new int[] { 1, 3 });
                p.AddCoeff(-2.2111166636879132e+02, new int[] { 3, 1 });
                p.AddCoeff(1.2082173198008954e+04, new int[] { 1, 5 });
                p.AddCoeff(6.6333499910637396e+03, new int[] { 3, 3 });
                p.AddCoeff(1.9900049973191219e+02, new int[] { 5, 1 });
                p.AddCoeff(-4.3725960145175263e+04, new int[] { 1, 7 });
                p.AddCoeff(-5.6383474924041787e+04, new int[] { 3, 5 });
                p.AddCoeff(-5.9700149919573657e+03, new int[] { 5, 3 });
                p.AddCoeff(7.6520430254056711e+04, new int[] { 1, 9 });
                p.AddCoeff(2.0405448067748456e+05, new int[] { 3, 7 });
                p.AddCoeff(5.0745127431637608e+04, new int[] { 5, 5 });
                p.AddCoeff(-6.3998905303392886e+04, new int[] { 1, 11 });
                p.AddCoeff(-3.5709534118559798e+05, new int[] { 3, 9 });
                p.AddCoeff(-1.8364903260973611e+05, new int[] { 5, 7 });
                p.AddCoeff(2.0512469648523361e+04, new int[] { 1, 13 });
                p.AddCoeff(2.9866155808250013e+05, new int[] { 3, 11 });
                p.AddCoeff(3.2138580706703819e+05, new int[] { 5, 9 });
                p.AddCoeff(-9.5724858359775683e+04, new int[] { 3, 13 });
                p.AddCoeff(-2.6879540227425012e+05, new int[] { 5, 11 });
                p.AddCoeff(8.6152372523798115e+04, new int[] { 5, 13 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{3EA5CC88-2893-4ABA-820D-0C17EC883289}"));
                OrthonormalPolynomials[177] = p;
                p.AddCoeff(-6.3543880053114435e-01, new int[] { 0, 0 });
                p.AddCoeff(4.9564226441429260e+01, new int[] { 0, 2 });
                p.AddCoeff(1.3344214811154031e+01, new int[] { 2, 0 });
                p.AddCoeff(-6.1955283051786574e+02, new int[] { 0, 4 });
                p.AddCoeff(-1.0408487552700145e+03, new int[] { 2, 2 });
                p.AddCoeff(-4.0032644433462094e+01, new int[] { 4, 0 });
                p.AddCoeff(2.8086394983476580e+03, new int[] { 0, 6 });
                p.AddCoeff(1.3010609440875181e+04, new int[] { 2, 4 });
                p.AddCoeff(3.1225462658100434e+03, new int[] { 4, 2 });
                p.AddCoeff(2.9357272584538869e+01, new int[] { 6, 0 });
                p.AddCoeff(-5.7175875502077324e+03, new int[] { 0, 8 });
                p.AddCoeff(-5.8981429465300819e+04, new int[] { 2, 6 });
                p.AddCoeff(-3.9031828322625542e+04, new int[] { 4, 4 });
                p.AddCoeff(-2.2898672615940318e+03, new int[] { 6, 2 });
                p.AddCoeff(5.3364150468605503e+03, new int[] { 0, 10 });
                p.AddCoeff(1.2006933855436238e+05, new int[] { 2, 8 });
                p.AddCoeff(1.7694428839590246e+05, new int[] { 4, 6 });
                p.AddCoeff(2.8623340769925397e+04, new int[] { 6, 4 });
                p.AddCoeff(-1.8596597890574645e+03, new int[] { 0, 12 });
                p.AddCoeff(-1.1206471598407156e+05, new int[] { 2, 10 });
                p.AddCoeff(-3.6020801566308714e+05, new int[] { 4, 8 });
                p.AddCoeff(-1.2975914482366180e+05, new int[] { 6, 6 });
                p.AddCoeff(3.9052855570206754e+04, new int[] { 2, 12 });
                p.AddCoeff(3.3619414795221467e+05, new int[] { 4, 10 });
                p.AddCoeff(2.6415254481959724e+05, new int[] { 6, 8 });
                p.AddCoeff(-1.1715856671062026e+05, new int[] { 4, 12 });
                p.AddCoeff(-2.4654237516495742e+05, new int[] { 6, 10 });
                p.AddCoeff(8.5916282254454859e+04, new int[] { 6, 12 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{14EB5BBB-7FA2-42AA-BD1F-4ED07926CEDC}"));
                OrthonormalPolynomials[178] = p;
                p.AddCoeff(5.4994705772402069e+01, new int[] { 1, 1 });
                p.AddCoeff(-1.1915519584020448e+03, new int[] { 1, 3 });
                p.AddCoeff(-4.9495235195161862e+02, new int[] { 3, 1 });
                p.AddCoeff(7.1493117504122690e+03, new int[] { 1, 5 });
                p.AddCoeff(1.0723967625618403e+04, new int[] { 3, 3 });
                p.AddCoeff(1.0888951742935610e+03, new int[] { 5, 1 });
                p.AddCoeff(-1.7362614251001225e+04, new int[] { 1, 7 });
                p.AddCoeff(-6.4343805753710421e+04, new int[] { 3, 5 });
                p.AddCoeff(-2.3592728776360488e+04, new int[] { 5, 3 });
                p.AddCoeff(-6.7407796503887108e+02, new int[] { 7, 1 });
                p.AddCoeff(1.8327203931612404e+04, new int[] { 1, 9 });
                p.AddCoeff(1.5626352825901102e+05, new int[] { 3, 7 });
                p.AddCoeff(1.4155637265816293e+05, new int[] { 5, 5 });
                p.AddCoeff(1.4605022575842207e+04, new int[] { 7, 3 });
                p.AddCoeff(-6.9976596829792815e+03, new int[] { 1, 11 });
                p.AddCoeff(-1.6494483538451163e+05, new int[] { 3, 9 });
                p.AddCoeff(-3.4377976216982425e+05, new int[] { 5, 7 });
                p.AddCoeff(-8.7630135455053240e+04, new int[] { 7, 5 });
                p.AddCoeff(6.2978937146813533e+04, new int[] { 3, 11 });
                p.AddCoeff(3.6287863784592560e+05, new int[] { 5, 9 });
                p.AddCoeff(2.1281604324798644e+05, new int[] { 7, 7 });
                p.AddCoeff(-1.3855366172298977e+05, new int[] { 5, 11 });
                p.AddCoeff(-2.2463915676176346e+05, new int[] { 7, 9 });
                p.AddCoeff(8.5771314399946050e+04, new int[] { 7, 11 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{41DC59D4-438E-400B-B4F9-A091677D9351}"));
                OrthonormalPolynomials[179] = p;
                p.AddCoeff(-6.3571545713896273e-01, new int[] { 0, 0 });
                p.AddCoeff(3.4964350142642950e+01, new int[] { 0, 2 });
                p.AddCoeff(2.2885756457002658e+01, new int[] { 2, 0 });
                p.AddCoeff(-3.0302436790290557e+02, new int[] { 0, 4 });
                p.AddCoeff(-1.2587166051351462e+03, new int[] { 2, 2 });
                p.AddCoeff(-1.2587166051351462e+02, new int[] { 4, 0 });
                p.AddCoeff(9.0907310370871670e+02, new int[] { 0, 6 });
                p.AddCoeff(1.0908877244504600e+04, new int[] { 2, 4 });
                p.AddCoeff(6.9229413282433041e+03, new int[] { 4, 2 });
                p.AddCoeff(2.1817754489009201e+02, new int[] { 6, 0 });
                p.AddCoeff(-1.1038744830748703e+03, new int[] { 0, 8 });
                p.AddCoeff(-3.2726631733513801e+04, new int[] { 2, 6 });
                p.AddCoeff(-5.9998824844775302e+04, new int[] { 4, 4 });
                p.AddCoeff(-1.1999764968955060e+04, new int[] { 6, 2 });
                p.AddCoeff(-1.1688082761969215e+02, new int[] { 8, 0 });
                p.AddCoeff(4.6608033729827856e+02, new int[] { 0, 10 });
                p.AddCoeff(3.9739481390695330e+04, new int[] { 2, 8 });
                p.AddCoeff(1.7999647453432591e+05, new int[] { 4, 6 });
                p.AddCoeff(1.0399796306427719e+05, new int[] { 6, 4 });
                p.AddCoeff(6.4284455190830681e+03, new int[] { 8, 2 });
                p.AddCoeff(-1.6778892142738028e+04, new int[] { 2, 10 });
                p.AddCoeff(-2.1856714764882432e+05, new int[] { 4, 8 });
                p.AddCoeff(-3.1199388919283157e+05, new int[] { 6, 6 });
                p.AddCoeff(-5.5713194498719924e+04, new int[] { 8, 4 });
                p.AddCoeff(9.2283906785059155e+04, new int[] { 4, 10 });
                p.AddCoeff(3.7884972259129548e+05, new int[] { 6, 8 });
                p.AddCoeff(1.6713958349615977e+05, new int[] { 8, 6 });
                p.AddCoeff(-1.5995877176076920e+05, new int[] { 6, 10 });
                p.AddCoeff(-2.0295520853105115e+05, new int[] { 8, 8 });
                p.AddCoeff(8.5692199157554930e+04, new int[] { 8, 10 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{3B7881F2-5090-4B68-9E3A-07B0495ADF3D}"));
                OrthonormalPolynomials[180] = p;
                p.AddCoeff(5.7534027099609375e+01, new int[] { 1, 1 });
                p.AddCoeff(-8.4383239746093750e+02, new int[] { 1, 3 });
                p.AddCoeff(-8.4383239746093750e+02, new int[] { 3, 1 });
                p.AddCoeff(3.2909463500976562e+03, new int[] { 1, 5 });
                p.AddCoeff(1.2376208496093750e+04, new int[] { 3, 3 });
                p.AddCoeff(3.2909463500976562e+03, new int[] { 5, 1 });
                p.AddCoeff(-4.7013519287109375e+03, new int[] { 1, 7 });
                p.AddCoeff(-4.8267213134765625e+04, new int[] { 3, 5 });
                p.AddCoeff(-4.8267213134765625e+04, new int[] { 5, 3 });
                p.AddCoeff(-4.7013519287109375e+03, new int[] { 7, 1 });
                p.AddCoeff(2.2200828552246094e+03, new int[] { 1, 9 });
                p.AddCoeff(6.8953161621093750e+04, new int[] { 3, 7 });
                p.AddCoeff(1.8824213122558594e+05, new int[] { 5, 5 });
                p.AddCoeff(6.8953161621093750e+04, new int[] { 7, 3 });
                p.AddCoeff(2.2200828552246094e+03, new int[] { 9, 1 });
                p.AddCoeff(-3.2561215209960938e+04, new int[] { 3, 9 });
                p.AddCoeff(-2.6891733032226562e+05, new int[] { 5, 7 });
                p.AddCoeff(-2.6891733032226562e+05, new int[] { 7, 5 });
                p.AddCoeff(-3.2561215209960938e+04, new int[] { 9, 3 });
                p.AddCoeff(1.2698873931884766e+05, new int[] { 5, 9 });
                p.AddCoeff(3.8416761474609375e+05, new int[] { 7, 7 });
                p.AddCoeff(1.2698873931884766e+05, new int[] { 9, 5 });
                p.AddCoeff(-1.8141248474121094e+05, new int[] { 7, 9 });
                p.AddCoeff(-1.8141248474121094e+05, new int[] { 9, 7 });
                p.AddCoeff(8.5667006683349609e+04, new int[] { 9, 9 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{81D93784-7166-4549-9923-4A8A005382CB}"));
                OrthonormalPolynomials[181] = p;
                p.AddCoeff(-6.3571545713896273e-01, new int[] { 0, 0 });
                p.AddCoeff(2.2885756457002658e+01, new int[] { 0, 2 });
                p.AddCoeff(3.4964350142642950e+01, new int[] { 2, 0 });
                p.AddCoeff(-1.2587166051351462e+02, new int[] { 0, 4 });
                p.AddCoeff(-1.2587166051351462e+03, new int[] { 2, 2 });
                p.AddCoeff(-3.0302436790290557e+02, new int[] { 4, 0 });
                p.AddCoeff(2.1817754489009201e+02, new int[] { 0, 6 });
                p.AddCoeff(6.9229413282433041e+03, new int[] { 2, 4 });
                p.AddCoeff(1.0908877244504600e+04, new int[] { 4, 2 });
                p.AddCoeff(9.0907310370871670e+02, new int[] { 6, 0 });
                p.AddCoeff(-1.1688082761969215e+02, new int[] { 0, 8 });
                p.AddCoeff(-1.1999764968955060e+04, new int[] { 2, 6 });
                p.AddCoeff(-5.9998824844775302e+04, new int[] { 4, 4 });
                p.AddCoeff(-3.2726631733513801e+04, new int[] { 6, 2 });
                p.AddCoeff(-1.1038744830748703e+03, new int[] { 8, 0 });
                p.AddCoeff(6.4284455190830681e+03, new int[] { 2, 8 });
                p.AddCoeff(1.0399796306427719e+05, new int[] { 4, 6 });
                p.AddCoeff(1.7999647453432591e+05, new int[] { 6, 4 });
                p.AddCoeff(3.9739481390695330e+04, new int[] { 8, 2 });
                p.AddCoeff(4.6608033729827856e+02, new int[] { 10, 0 });
                p.AddCoeff(-5.5713194498719924e+04, new int[] { 4, 8 });
                p.AddCoeff(-3.1199388919283157e+05, new int[] { 6, 6 });
                p.AddCoeff(-2.1856714764882432e+05, new int[] { 8, 4 });
                p.AddCoeff(-1.6778892142738028e+04, new int[] { 10, 2 });
                p.AddCoeff(1.6713958349615977e+05, new int[] { 6, 8 });
                p.AddCoeff(3.7884972259129548e+05, new int[] { 8, 6 });
                p.AddCoeff(9.2283906785059155e+04, new int[] { 10, 4 });
                p.AddCoeff(-2.0295520853105115e+05, new int[] { 8, 8 });
                p.AddCoeff(-1.5995877176076920e+05, new int[] { 10, 6 });
                p.AddCoeff(8.5692199157554930e+04, new int[] { 10, 8 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{92AF5E39-D6B3-44D4-B487-EF0123331A32}"));
                OrthonormalPolynomials[182] = p;
                p.AddCoeff(5.4994705772402069e+01, new int[] { 1, 1 });
                p.AddCoeff(-4.9495235195161862e+02, new int[] { 1, 3 });
                p.AddCoeff(-1.1915519584020448e+03, new int[] { 3, 1 });
                p.AddCoeff(1.0888951742935610e+03, new int[] { 1, 5 });
                p.AddCoeff(1.0723967625618403e+04, new int[] { 3, 3 });
                p.AddCoeff(7.1493117504122690e+03, new int[] { 5, 1 });
                p.AddCoeff(-6.7407796503887108e+02, new int[] { 1, 7 });
                p.AddCoeff(-2.3592728776360488e+04, new int[] { 3, 5 });
                p.AddCoeff(-6.4343805753710421e+04, new int[] { 5, 3 });
                p.AddCoeff(-1.7362614251001225e+04, new int[] { 7, 1 });
                p.AddCoeff(1.4605022575842207e+04, new int[] { 3, 7 });
                p.AddCoeff(1.4155637265816293e+05, new int[] { 5, 5 });
                p.AddCoeff(1.5626352825901102e+05, new int[] { 7, 3 });
                p.AddCoeff(1.8327203931612404e+04, new int[] { 9, 1 });
                p.AddCoeff(-8.7630135455053240e+04, new int[] { 5, 7 });
                p.AddCoeff(-3.4377976216982425e+05, new int[] { 7, 5 });
                p.AddCoeff(-1.6494483538451163e+05, new int[] { 9, 3 });
                p.AddCoeff(-6.9976596829792815e+03, new int[] { 11, 1 });
                p.AddCoeff(2.1281604324798644e+05, new int[] { 7, 7 });
                p.AddCoeff(3.6287863784592560e+05, new int[] { 9, 5 });
                p.AddCoeff(6.2978937146813533e+04, new int[] { 11, 3 });
                p.AddCoeff(-2.2463915676176346e+05, new int[] { 9, 7 });
                p.AddCoeff(-1.3855366172298977e+05, new int[] { 11, 5 });
                p.AddCoeff(8.5771314399946050e+04, new int[] { 11, 7 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{94776049-9AD2-4436-ADF9-397C1FFCC7A9}"));
                OrthonormalPolynomials[183] = p;
                p.AddCoeff(-6.3543880053114435e-01, new int[] { 0, 0 });
                p.AddCoeff(1.3344214811154031e+01, new int[] { 0, 2 });
                p.AddCoeff(4.9564226441429260e+01, new int[] { 2, 0 });
                p.AddCoeff(-4.0032644433462094e+01, new int[] { 0, 4 });
                p.AddCoeff(-1.0408487552700145e+03, new int[] { 2, 2 });
                p.AddCoeff(-6.1955283051786574e+02, new int[] { 4, 0 });
                p.AddCoeff(2.9357272584538869e+01, new int[] { 0, 6 });
                p.AddCoeff(3.1225462658100434e+03, new int[] { 2, 4 });
                p.AddCoeff(1.3010609440875181e+04, new int[] { 4, 2 });
                p.AddCoeff(2.8086394983476580e+03, new int[] { 6, 0 });
                p.AddCoeff(-2.2898672615940318e+03, new int[] { 2, 6 });
                p.AddCoeff(-3.9031828322625542e+04, new int[] { 4, 4 });
                p.AddCoeff(-5.8981429465300819e+04, new int[] { 6, 2 });
                p.AddCoeff(-5.7175875502077324e+03, new int[] { 8, 0 });
                p.AddCoeff(2.8623340769925397e+04, new int[] { 4, 6 });
                p.AddCoeff(1.7694428839590246e+05, new int[] { 6, 4 });
                p.AddCoeff(1.2006933855436238e+05, new int[] { 8, 2 });
                p.AddCoeff(5.3364150468605503e+03, new int[] { 10, 0 });
                p.AddCoeff(-1.2975914482366180e+05, new int[] { 6, 6 });
                p.AddCoeff(-3.6020801566308714e+05, new int[] { 8, 4 });
                p.AddCoeff(-1.1206471598407156e+05, new int[] { 10, 2 });
                p.AddCoeff(-1.8596597890574645e+03, new int[] { 12, 0 });
                p.AddCoeff(2.6415254481959724e+05, new int[] { 8, 6 });
                p.AddCoeff(3.3619414795221467e+05, new int[] { 10, 4 });
                p.AddCoeff(3.9052855570206754e+04, new int[] { 12, 2 });
                p.AddCoeff(-2.4654237516495742e+05, new int[] { 10, 6 });
                p.AddCoeff(-1.1715856671062026e+05, new int[] { 12, 4 });
                p.AddCoeff(8.5916282254454859e+04, new int[] { 12, 6 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{E92F8AA7-0AFE-4AF8-9A3E-56B2EF838161}"));
                OrthonormalPolynomials[184] = p;
                p.AddCoeff(4.7381071364740997e+01, new int[] { 1, 1 });
                p.AddCoeff(-2.2111166636879132e+02, new int[] { 1, 3 });
                p.AddCoeff(-1.4214321409422299e+03, new int[] { 3, 1 });
                p.AddCoeff(1.9900049973191219e+02, new int[] { 1, 5 });
                p.AddCoeff(6.6333499910637396e+03, new int[] { 3, 3 });
                p.AddCoeff(1.2082173198008954e+04, new int[] { 5, 1 });
                p.AddCoeff(-5.9700149919573657e+03, new int[] { 3, 5 });
                p.AddCoeff(-5.6383474924041787e+04, new int[] { 5, 3 });
                p.AddCoeff(-4.3725960145175263e+04, new int[] { 7, 1 });
                p.AddCoeff(5.0745127431637608e+04, new int[] { 5, 5 });
                p.AddCoeff(2.0405448067748456e+05, new int[] { 7, 3 });
                p.AddCoeff(7.6520430254056711e+04, new int[] { 9, 1 });
                p.AddCoeff(-1.8364903260973611e+05, new int[] { 7, 5 });
                p.AddCoeff(-3.5709534118559798e+05, new int[] { 9, 3 });
                p.AddCoeff(-6.3998905303392886e+04, new int[] { 11, 1 });
                p.AddCoeff(3.2138580706703819e+05, new int[] { 9, 5 });
                p.AddCoeff(2.9866155808250013e+05, new int[] { 11, 3 });
                p.AddCoeff(2.0512469648523361e+04, new int[] { 13, 1 });
                p.AddCoeff(-2.6879540227425012e+05, new int[] { 11, 5 });
                p.AddCoeff(-9.5724858359775683e+04, new int[] { 13, 3 });
                p.AddCoeff(8.6152372523798115e+04, new int[] { 13, 5 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{E27F0DA5-0BF6-4D45-8BD1-C462279F0EA3}"));
                OrthonormalPolynomials[185] = p;
                p.AddCoeff(-6.3452518677814697e-01, new int[] { 0, 0 });
                p.AddCoeff(6.3452518677814697e+00, new int[] { 0, 2 });
                p.AddCoeff(6.6625144611705432e+01, new int[] { 2, 0 });
                p.AddCoeff(-7.4027938457450480e+00, new int[] { 0, 4 });
                p.AddCoeff(-6.6625144611705432e+02, new int[] { 2, 2 });
                p.AddCoeff(-1.1326274583989923e+03, new int[] { 4, 0 });
                p.AddCoeff(7.7729335380323004e+02, new int[] { 2, 4 });
                p.AddCoeff(1.1326274583989923e+04, new int[] { 4, 2 });
                p.AddCoeff(7.1733072365269515e+03, new int[] { 6, 0 });
                p.AddCoeff(-1.3213987014654911e+04, new int[] { 4, 4 });
                p.AddCoeff(-7.1733072365269515e+04, new int[] { 6, 2 });
                p.AddCoeff(-2.1519921709580855e+04, new int[] { 8, 0 });
                p.AddCoeff(8.3688584426147768e+04, new int[] { 6, 4 });
                p.AddCoeff(2.1519921709580855e+05, new int[] { 8, 2 });
                p.AddCoeff(3.2997213288023977e+04, new int[] { 10, 0 });
                p.AddCoeff(-2.5106575327844330e+05, new int[] { 8, 4 });
                p.AddCoeff(-3.2997213288023977e+05, new int[] { 10, 2 });
                p.AddCoeff(-2.4997888854563619e+04, new int[] { 12, 0 });
                p.AddCoeff(3.8496748836027973e+05, new int[] { 10, 4 });
                p.AddCoeff(2.4997888854563619e+05, new int[] { 12, 2 });
                p.AddCoeff(7.4169560337716232e+03, new int[] { 14, 0 });
                p.AddCoeff(-2.9164203663657556e+05, new int[] { 12, 4 });
                p.AddCoeff(-7.4169560337716232e+04, new int[] { 14, 2 });
                p.AddCoeff(8.6531153727335604e+04, new int[] { 14, 4 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{DB857498-2D5E-4577-891A-F5214E299ECF}"));
                OrthonormalPolynomials[186] = p;
                p.AddCoeff(3.4714405267160477e+01, new int[] { 1, 1 });
                p.AddCoeff(-5.7857342111934128e+01, new int[] { 1, 3 });
                p.AddCoeff(-1.3770047422640322e+03, new int[] { 3, 1 });
                p.AddCoeff(2.2950079037733871e+03, new int[] { 3, 3 });
                p.AddCoeff(1.5697854061809967e+04, new int[] { 5, 1 });
                p.AddCoeff(-2.6163090103016612e+04, new int[] { 5, 3 });
                p.AddCoeff(-7.8489270309049837e+04, new int[] { 7, 1 });
                p.AddCoeff(1.3081545051508306e+05, new int[] { 7, 3 });
                p.AddCoeff(2.0058369078979403e+05, new int[] { 9, 1 });
                p.AddCoeff(-3.3430615131632338e+05, new int[] { 9, 3 });
                p.AddCoeff(-2.7352321471335549e+05, new int[] { 11, 1 });
                p.AddCoeff(4.5587202452225916e+05, new int[] { 11, 3 });
                p.AddCoeff(1.8936222557078457e+05, new int[] { 13, 1 });
                p.AddCoeff(-3.1560370928464095e+05, new int[] { 13, 3 });
                p.AddCoeff(-5.2300043252883358e+04, new int[] { 15, 1 });
                p.AddCoeff(8.7166738754805597e+04, new int[] { 15, 3 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{2E7FA6DA-2106-4D0A-A27B-EA8084B6D786}"));
                OrthonormalPolynomials[187] = p;
                p.AddCoeff(-6.3063866915672383e-01, new int[] { 0, 0 });
                p.AddCoeff(1.8919160074701715e+00, new int[] { 0, 2 });
                p.AddCoeff(8.5766859005314440e+01, new int[] { 2, 0 });
                p.AddCoeff(-2.5730057701594332e+02, new int[] { 2, 2 });
                p.AddCoeff(-1.9011653746178034e+03, new int[] { 4, 0 });
                p.AddCoeff(5.7034961238534103e+03, new int[] { 4, 2 });
                p.AddCoeff(1.5969789146789549e+04, new int[] { 6, 0 });
                p.AddCoeff(-4.7909367440368646e+04, new int[] { 6, 2 });
                p.AddCoeff(-6.5590205424314218e+04, new int[] { 8, 0 });
                p.AddCoeff(1.9677061627294265e+05, new int[] { 8, 2 });
                p.AddCoeff(1.4575601205403160e+05, new int[] { 10, 0 });
                p.AddCoeff(-4.3726803616209479e+05, new int[] { 10, 2 });
                p.AddCoeff(-1.7888237842994787e+05, new int[] { 12, 0 });
                p.AddCoeff(5.3664713528984360e+05, new int[] { 12, 2 });
                p.AddCoeff(1.1401294449381293e+05, new int[] { 14, 0 });
                p.AddCoeff(-3.4203883348143878e+05, new int[] { 14, 2 });
                p.AddCoeff(-2.9453343994235006e+04, new int[] { 16, 0 });
                p.AddCoeff(8.8360031982705018e+04, new int[] { 16, 2 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{52A125DB-7699-4791-9C28-6C093AC3B95D}"));
                OrthonormalPolynomials[188] = p;
                p.AddCoeff(1.7104571213411717e+01, new int[] { 1, 1 });
                p.AddCoeff(-8.6663160814619365e+02, new int[] { 3, 1 });
                p.AddCoeff(1.2739484639749047e+04, new int[] { 5, 1 });
                p.AddCoeff(-8.3716613346922306e+04, new int[] { 7, 1 });
                p.AddCoeff(2.9068268523236912e+05, new int[] { 9, 1 });
                p.AddCoeff(-5.7079509100174300e+05, new int[] { 11, 1 });
                p.AddCoeff(6.3665606304040565e+05, new int[] { 13, 1 });
                p.AddCoeff(-3.7593024674766810e+05, new int[] { 15, 1 });
                p.AddCoeff(9.1218368696125347e+04, new int[] { 17, 1 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{5F19381D-6841-4B7A-8CD7-28E1249DB13E}"));
                OrthonormalPolynomials[189] = p;
                p.AddCoeff(-5.6408675045604599e-01, new int[] { 0, 0 });
                p.AddCoeff(9.6458834327983865e+01, new int[] { 2, 0 });
                p.AddCoeff(-2.7008473611835482e+03, new int[] { 4, 0 });
                p.AddCoeff(2.8989095010036751e+04, new int[] { 6, 0 });
                p.AddCoeff(-1.5529872326805402e+05, new int[] { 8, 0 });
                p.AddCoeff(4.6589616980416207e+05, new int[] { 10, 0 });
                p.AddCoeff(-8.1884781359519394e+05, new int[] { 12, 0 });
                p.AddCoeff(8.3684446883904435e+05, new int[] { 14, 0 });
                p.AddCoeff(-4.6026445786147439e+05, new int[] { 16, 0 });
                p.AddCoeff(1.0528925506635035e+05, new int[] { 18, 0 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{370846CB-0146-4134-AE1C-8DB479E38F82}"));
                OrthonormalPolynomials[190] = p;
                p.AddCoeff(-1.1003502370758957e+01, new int[] { 0, 1 });
                p.AddCoeff(6.9322064935781432e+02, new int[] { 0, 3 });
                p.AddCoeff(-1.2755259948183784e+04, new int[] { 0, 5 });
                p.AddCoeff(1.0629383290153153e+05, new int[] { 0, 7 });
                p.AddCoeff(-4.7832224805689188e+05, new int[] { 0, 9 });
                p.AddCoeff(1.2610313812408968e+06, new int[] { 0, 11 });
                p.AddCoeff(-2.0047165547932205e+06, new int[] { 0, 13 });
                p.AddCoeff(1.8901613230907508e+06, new int[] { 0, 15 });
                p.AddCoeff(-9.7287715159082760e+05, new int[] { 0, 17 });
                p.AddCoeff(2.1050558250795685e+05, new int[] { 0, 19 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{F0B2ACA6-15FD-499D-81B0-E23465B26982}"));
                OrthonormalPolynomials[191] = p;
                p.AddCoeff(-9.7702691166629822e-01, new int[] { 1, 0 });
                p.AddCoeff(1.6707160189493700e+02, new int[] { 1, 2 });
                p.AddCoeff(-4.6780048530582359e+03, new int[] { 1, 4 });
                p.AddCoeff(5.0210585422825065e+04, new int[] { 1, 6 });
                p.AddCoeff(-2.6898527905084856e+05, new int[] { 1, 8 });
                p.AddCoeff(8.0695583715254569e+05, new int[] { 1, 10 });
                p.AddCoeff(-1.4182860168135652e+06, new int[] { 1, 12 });
                p.AddCoeff(1.4494571380622149e+06, new int[] { 1, 14 });
                p.AddCoeff(-7.9720142593421822e+05, new int[] { 1, 16 });
                p.AddCoeff(1.8236633926599763e+05, new int[] { 1, 18 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{1078CA75-7844-45CA-A93F-09DC454DBF33}"));
                OrthonormalPolynomials[192] = p;
                p.AddCoeff(-1.1040953242260395e+01, new int[] { 0, 1 });
                p.AddCoeff(5.5940829760785999e+02, new int[] { 0, 3 });
                p.AddCoeff(3.3122859726781184e+01, new int[] { 2, 1 });
                p.AddCoeff(-8.2233019748355419e+03, new int[] { 0, 5 });
                p.AddCoeff(-1.6782248928235800e+03, new int[] { 2, 3 });
                p.AddCoeff(5.4038841548919275e+04, new int[] { 0, 7 });
                p.AddCoeff(2.4669905924506626e+04, new int[] { 2, 5 });
                p.AddCoeff(-1.8763486648930304e+05, new int[] { 0, 9 });
                p.AddCoeff(-1.6211652464675783e+05, new int[] { 2, 7 });
                p.AddCoeff(3.6844664692444960e+05, new int[] { 0, 11 });
                p.AddCoeff(5.6290459946790912e+05, new int[] { 2, 9 });
                p.AddCoeff(-4.1095972156957840e+05, new int[] { 0, 13 });
                p.AddCoeff(-1.1053399407733488e+06, new int[] { 2, 11 });
                p.AddCoeff(2.4266193083156058e+05, new int[] { 0, 15 });
                p.AddCoeff(1.2328791647087352e+06, new int[] { 2, 13 });
                p.AddCoeff(-5.8881203804716906e+04, new int[] { 0, 17 });
                p.AddCoeff(-7.2798579249468174e+05, new int[] { 2, 15 });
                p.AddCoeff(1.7664361141415072e+05, new int[] { 2, 17 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{0F0279BA-2948-4232-9E49-78659EBEB683}"));
                OrthonormalPolynomials[193] = p;
                p.AddCoeff(-2.2385452086233647e+00, new int[] { 1, 0 });
                p.AddCoeff(3.0444214837277760e+02, new int[] { 1, 2 });
                p.AddCoeff(3.7309086810389412e+00, new int[] { 3, 0 });
                p.AddCoeff(-6.7484676222632369e+03, new int[] { 1, 4 });
                p.AddCoeff(-5.0740358062129600e+02, new int[] { 3, 2 });
                p.AddCoeff(5.6687128027011190e+04, new int[] { 1, 6 });
                p.AddCoeff(1.1247446037105395e+04, new int[] { 3, 4 });
                p.AddCoeff(-2.3282213296808167e+05, new int[] { 1, 8 });
                p.AddCoeff(-9.4478546711685316e+04, new int[] { 3, 6 });
                p.AddCoeff(5.1738251770684816e+05, new int[] { 1, 10 });
                p.AddCoeff(3.8803688828013612e+05, new int[] { 3, 8 });
                p.AddCoeff(-6.3496945354931365e+05, new int[] { 1, 12 });
                p.AddCoeff(-8.6230419617808027e+05, new int[] { 3, 10 });
                p.AddCoeff(4.0470580555890321e+05, new int[] { 1, 14 });
                p.AddCoeff(1.0582824225821894e+06, new int[] { 3, 12 });
                p.AddCoeff(-1.0454899976938333e+05, new int[] { 1, 16 });
                p.AddCoeff(-6.7450967593150534e+05, new int[] { 3, 14 });
                p.AddCoeff(1.7424833294897221e+05, new int[] { 3, 16 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{63E7641A-CCF6-4E29-8EC9-C0FE9073ECFB}"));
                OrthonormalPolynomials[194] = p;
                p.AddCoeff(-9.8406089194732886e+00, new int[] { 0, 1 });
                p.AddCoeff(3.9034415380577378e+02, new int[] { 0, 3 });
                p.AddCoeff(9.8406089194732886e+01, new int[] { 2, 1 });
                p.AddCoeff(-4.4499233533858211e+03, new int[] { 0, 5 });
                p.AddCoeff(-3.9034415380577378e+03, new int[] { 2, 3 });
                p.AddCoeff(-1.1480710406052170e+02, new int[] { 4, 1 });
                p.AddCoeff(2.2249616766929105e+04, new int[] { 0, 7 });
                p.AddCoeff(4.4499233533858211e+04, new int[] { 2, 5 });
                p.AddCoeff(4.5540151277340274e+03, new int[] { 4, 3 });
                p.AddCoeff(-5.6860131737707714e+04, new int[] { 0, 9 });
                p.AddCoeff(-2.2249616766929105e+05, new int[] { 2, 7 });
                p.AddCoeff(-5.1915772456167913e+04, new int[] { 4, 5 });
                p.AddCoeff(7.7536543278692337e+04, new int[] { 0, 11 });
                p.AddCoeff(5.6860131737707714e+05, new int[] { 2, 9 });
                p.AddCoeff(2.5957886228083956e+05, new int[] { 4, 7 });
                p.AddCoeff(-5.3679145346787003e+04, new int[] { 0, 13 });
                p.AddCoeff(-7.7536543278692337e+05, new int[] { 2, 11 });
                p.AddCoeff(-6.6336820360659000e+05, new int[] { 4, 9 });
                p.AddCoeff(1.4825668714826886e+04, new int[] { 0, 15 });
                p.AddCoeff(5.3679145346787003e+05, new int[] { 2, 13 });
                p.AddCoeff(9.0459300491807727e+05, new int[] { 4, 11 });
                p.AddCoeff(-1.4825668714826886e+05, new int[] { 2, 15 });
                p.AddCoeff(-6.2625669571251503e+05, new int[] { 4, 13 });
                p.AddCoeff(1.7296613500631368e+05, new int[] { 4, 15 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{CD3B74DE-00DB-4DE0-9365-C62867339F77}"));
                OrthonormalPolynomials[195] = p;
                p.AddCoeff(-3.5074699409554877e+00, new int[] { 1, 0 });
                p.AddCoeff(3.6828434380032621e+02, new int[] { 1, 2 });
                p.AddCoeff(1.6368193057792276e+01, new int[] { 3, 0 });
                p.AddCoeff(-6.2608338446055456e+03, new int[] { 1, 4 });
                p.AddCoeff(-1.7186602710681890e+03, new int[] { 3, 2 });
                p.AddCoeff(-1.4731373752013048e+01, new int[] { 5, 0 });
                p.AddCoeff(3.9651947682501789e+04, new int[] { 1, 6 });
                p.AddCoeff(2.9217224608159213e+04, new int[] { 3, 4 });
                p.AddCoeff(1.5467942439613701e+03, new int[] { 5, 2 });
                p.AddCoeff(-1.1895584304750537e+05, new int[] { 1, 8 });
                p.AddCoeff(-1.8504242251834168e+05, new int[] { 3, 6 });
                p.AddCoeff(-2.6295502147343292e+04, new int[] { 5, 4 });
                p.AddCoeff(1.8239895933950823e+05, new int[] { 1, 10 });
                p.AddCoeff(5.5512726755502504e+05, new int[] { 3, 8 });
                p.AddCoeff(1.6653818026650751e+05, new int[] { 5, 6 });
                p.AddCoeff(-1.3818102980265775e+05, new int[] { 1, 12 });
                p.AddCoeff(-8.5119514358437173e+05, new int[] { 3, 10 });
                p.AddCoeff(-4.9961454079952254e+05, new int[] { 5, 8 });
                p.AddCoeff(4.0998767084305046e+04, new int[] { 1, 14 });
                p.AddCoeff(6.4484480574573616e+05, new int[] { 3, 12 });
                p.AddCoeff(7.6607562922593456e+05, new int[] { 5, 10 });
                p.AddCoeff(-1.9132757972675688e+05, new int[] { 3, 14 });
                p.AddCoeff(-5.8036032517116255e+05, new int[] { 5, 12 });
                p.AddCoeff(1.7219482175408119e+05, new int[] { 5, 14 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{FA4DCC48-F734-4729-B0A3-6B2AA143C1FF}"));
                OrthonormalPolynomials[196] = p;
                p.AddCoeff(-8.5847758434404250e+00, new int[] { 0, 1 });
                p.AddCoeff(2.5754327530321275e+02, new int[] { 0, 3 });
                p.AddCoeff(1.8028029271224893e+02, new int[] { 2, 1 });
                p.AddCoeff(-2.1891178400773084e+03, new int[] { 0, 5 });
                p.AddCoeff(-5.4084087813674678e+03, new int[] { 2, 3 });
                p.AddCoeff(-5.4084087813674678e+02, new int[] { 4, 1 });
                p.AddCoeff(7.9225217069464494e+03, new int[] { 0, 7 });
                p.AddCoeff(4.5971474641623476e+04, new int[] { 2, 5 });
                p.AddCoeff(1.6225226344102403e+04, new int[] { 4, 3 });
                p.AddCoeff(3.9661664396694764e+02, new int[] { 6, 1 });
                p.AddCoeff(-1.3864412987156286e+04, new int[] { 0, 9 });
                p.AddCoeff(-1.6637295584587544e+05, new int[] { 2, 7 });
                p.AddCoeff(-1.3791442392487043e+05, new int[] { 4, 5 });
                p.AddCoeff(-1.1898499319008429e+04, new int[] { 6, 3 });
                p.AddCoeff(1.1595690861985258e+04, new int[] { 0, 11 });
                p.AddCoeff(2.9115267273028201e+05, new int[] { 2, 9 });
                p.AddCoeff(4.9911886753762631e+05, new int[] { 4, 7 });
                p.AddCoeff(1.0113724421157165e+05, new int[] { 6, 5 });
                p.AddCoeff(-3.7165675839696339e+03, new int[] { 0, 13 });
                p.AddCoeff(-2.4350950810169041e+05, new int[] { 2, 11 });
                p.AddCoeff(-8.7345801819084604e+05, new int[] { 4, 9 });
                p.AddCoeff(-3.6602050286092596e+05, new int[] { 6, 7 });
                p.AddCoeff(7.8047919263362312e+04, new int[] { 2, 13 });
                p.AddCoeff(7.3052852430507124e+05, new int[] { 4, 11 });
                p.AddCoeff(6.4053588000662043e+05, new int[] { 6, 9 });
                p.AddCoeff(-2.3414375779008693e+05, new int[] { 4, 13 });
                p.AddCoeff(-5.3572091782371891e+05, new int[] { 6, 11 });
                p.AddCoeff(1.7170542237939709e+05, new int[] { 6, 13 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{0DC0FBF8-0073-4B9C-A17F-05CFF8BEAC11}"));
                OrthonormalPolynomials[197] = p;
                p.AddCoeff(-4.7779953543223519e+00, new int[] { 1, 0 });
                p.AddCoeff(3.7268363763714345e+02, new int[] { 1, 2 });
                p.AddCoeff(4.3001958188901167e+01, new int[] { 3, 0 });
                p.AddCoeff(-4.6585454704642931e+03, new int[] { 1, 4 });
                p.AddCoeff(-3.3541527387342910e+03, new int[] { 3, 2 });
                p.AddCoeff(-9.4604308015582568e+01, new int[] { 5, 0 });
                p.AddCoeff(2.1118739466104795e+04, new int[] { 1, 6 });
                p.AddCoeff(4.1926909234178638e+04, new int[] { 3, 4 });
                p.AddCoeff(7.3791360252154403e+03, new int[] { 5, 2 });
                p.AddCoeff(5.8564571628693971e+01, new int[] { 7, 0 });
                p.AddCoeff(-4.2991719627427619e+04, new int[] { 1, 8 });
                p.AddCoeff(-1.9006865519494316e+05, new int[] { 3, 6 });
                p.AddCoeff(-9.2239200315193004e+04, new int[] { 5, 4 });
                p.AddCoeff(-4.5680365870381297e+03, new int[] { 7, 2 });
                p.AddCoeff(4.0125604985599111e+04, new int[] { 1, 10 });
                p.AddCoeff(3.8692547664684857e+05, new int[] { 3, 8 });
                p.AddCoeff(4.1815104142887495e+05, new int[] { 5, 6 });
                p.AddCoeff(5.7100457337976621e+04, new int[] { 7, 4 });
                p.AddCoeff(-1.3983165373769387e+04, new int[] { 1, 12 });
                p.AddCoeff(-3.6113044487039200e+05, new int[] { 3, 10 });
                p.AddCoeff(-8.5123604862306686e+05, new int[] { 5, 8 });
                p.AddCoeff(-2.5885540659882735e+05, new int[] { 7, 6 });
                p.AddCoeff(1.2584848836392449e+05, new int[] { 3, 12 });
                p.AddCoeff(7.9448697871486241e+05, new int[] { 5, 10 });
                p.AddCoeff(5.2695564914761282e+05, new int[] { 7, 8 });
                p.AddCoeff(-2.7686667440063387e+05, new int[] { 5, 12 });
                p.AddCoeff(-4.9182527253777197e+05, new int[] { 7, 10 });
                p.AddCoeff(1.7139365558134478e+05, new int[] { 7, 12 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{5CB7B231-280C-4113-A42D-765970CFA233}"));
                OrthonormalPolynomials[198] = p;
                p.AddCoeff(-7.3182918850987587e+00, new int[] { 0, 1 });
                p.AddCoeff(1.5856299084380644e+02, new int[] { 0, 3 });
                p.AddCoeff(2.6345850786355531e+02, new int[] { 2, 1 });
                p.AddCoeff(-9.5137794506283863e+02, new int[] { 0, 5 });
                p.AddCoeff(-5.7082676703770318e+03, new int[] { 2, 3 });
                p.AddCoeff(-1.4490217932495542e+03, new int[] { 4, 1 });
                p.AddCoeff(2.3104892951526081e+03, new int[] { 0, 7 });
                p.AddCoeff(3.4249606022262191e+04, new int[] { 2, 5 });
                p.AddCoeff(3.1395472187073675e+04, new int[] { 4, 3 });
                p.AddCoeff(2.5116377749658940e+03, new int[] { 6, 1 });
                p.AddCoeff(-2.4388498115499752e+03, new int[] { 0, 9 });
                p.AddCoeff(-8.3177614625493891e+04, new int[] { 2, 7 });
                p.AddCoeff(-1.8837283312244205e+05, new int[] { 4, 5 });
                p.AddCoeff(-5.4418818457594369e+04, new int[] { 6, 3 });
                p.AddCoeff(-1.3455202365888718e+03, new int[] { 8, 1 });
                p.AddCoeff(9.3119720077362690e+02, new int[] { 0, 11 });
                p.AddCoeff(8.7798593215799107e+04, new int[] { 2, 9 });
                p.AddCoeff(4.5747688044021640e+05, new int[] { 4, 7 });
                p.AddCoeff(3.2651291074556622e+05, new int[] { 6, 5 });
                p.AddCoeff(2.9152938459425555e+04, new int[] { 8, 3 });
                p.AddCoeff(-3.3523099227850568e+04, new int[] { 2, 11 });
                p.AddCoeff(-4.8289226268689509e+05, new int[] { 4, 9 });
                p.AddCoeff(-7.9295992609637510e+05, new int[] { 6, 7 });
                p.AddCoeff(-1.7491763075655333e+05, new int[] { 8, 5 });
                p.AddCoeff(1.8437704575317813e+05, new int[] { 4, 11 });
                p.AddCoeff(8.3701325532395149e+05, new int[] { 6, 9 });
                p.AddCoeff(4.2479996040877237e+05, new int[] { 8, 7 });
                p.AddCoeff(-3.1958687930550875e+05, new int[] { 6, 11 });
                p.AddCoeff(-4.4839995820925973e+05, new int[] { 8, 9 });
                p.AddCoeff(1.7120725677080826e+05, new int[] { 8, 11 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{CD13A4C3-C4FE-4C2B-9105-024C77F6B564}"));
                OrthonormalPolynomials[199] = p;
                p.AddCoeff(-6.0486383748423868e+00, new int[] { 1, 0 });
                p.AddCoeff(3.3267511061633127e+02, new int[] { 1, 2 });
                p.AddCoeff(8.8713362831021673e+01, new int[] { 3, 0 });
                p.AddCoeff(-2.8831842920082044e+03, new int[] { 1, 4 });
                p.AddCoeff(-4.8792349557061920e+03, new int[] { 3, 2 });
                p.AddCoeff(-3.4598211504098452e+02, new int[] { 5, 0 });
                p.AddCoeff(8.6495528760246131e+03, new int[] { 1, 6 });
                p.AddCoeff(4.2286702949453664e+04, new int[] { 3, 4 });
                p.AddCoeff(1.9029016327254149e+04, new int[] { 5, 2 });
                p.AddCoeff(4.9426016434426361e+02, new int[] { 7, 0 });
                p.AddCoeff(-1.0503028492315602e+04, new int[] { 1, 8 });
                p.AddCoeff(-1.2686010884836099e+05, new int[] { 3, 6 });
                p.AddCoeff(-1.6491814150286929e+05, new int[] { 5, 4 });
                p.AddCoeff(-2.7184309038934498e+04, new int[] { 7, 2 });
                p.AddCoeff(-2.3340063316256893e+02, new int[] { 9, 0 });
                p.AddCoeff(4.4346120300888096e+03, new int[] { 1, 10 });
                p.AddCoeff(1.5404441788729549e+05, new int[] { 3, 8 });
                p.AddCoeff(4.9475442450860787e+05, new int[] { 5, 6 });
                p.AddCoeff(2.3559734500409899e+05, new int[] { 7, 4 });
                p.AddCoeff(1.2837034823941291e+04, new int[] { 9, 2 });
                p.AddCoeff(-6.5040976441302540e+04, new int[] { 3, 10 });
                p.AddCoeff(-6.0077322976045241e+05, new int[] { 5, 8 });
                p.AddCoeff(-7.0679203501229696e+05, new int[] { 7, 6 });
                p.AddCoeff(-1.1125430180749119e+05, new int[] { 9, 4 });
                p.AddCoeff(2.5365980812107991e+05, new int[] { 5, 10 });
                p.AddCoeff(8.5824747108636059e+05, new int[] { 7, 8 });
                p.AddCoeff(3.3376290542247356e+05, new int[] { 9, 6 });
                p.AddCoeff(-3.6237115445868558e+05, new int[] { 7, 10 });
                p.AddCoeff(-4.0528352801300361e+05, new int[] { 9, 8 });
                p.AddCoeff(1.7111971182771264e+05, new int[] { 9, 10 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{14933C74-5896-4CCF-BDC0-79B28CFB2D7E}"));
                OrthonormalPolynomials[200] = p;
                p.AddCoeff(-6.0486383748423868e+00, new int[] { 0, 1 });
                p.AddCoeff(8.8713362831021673e+01, new int[] { 0, 3 });
                p.AddCoeff(3.3267511061633127e+02, new int[] { 2, 1 });
                p.AddCoeff(-3.4598211504098452e+02, new int[] { 0, 5 });
                p.AddCoeff(-4.8792349557061920e+03, new int[] { 2, 3 });
                p.AddCoeff(-2.8831842920082044e+03, new int[] { 4, 1 });
                p.AddCoeff(4.9426016434426361e+02, new int[] { 0, 7 });
                p.AddCoeff(1.9029016327254149e+04, new int[] { 2, 5 });
                p.AddCoeff(4.2286702949453664e+04, new int[] { 4, 3 });
                p.AddCoeff(8.6495528760246131e+03, new int[] { 6, 1 });
                p.AddCoeff(-2.3340063316256893e+02, new int[] { 0, 9 });
                p.AddCoeff(-2.7184309038934498e+04, new int[] { 2, 7 });
                p.AddCoeff(-1.6491814150286929e+05, new int[] { 4, 5 });
                p.AddCoeff(-1.2686010884836099e+05, new int[] { 6, 3 });
                p.AddCoeff(-1.0503028492315602e+04, new int[] { 8, 1 });
                p.AddCoeff(1.2837034823941291e+04, new int[] { 2, 9 });
                p.AddCoeff(2.3559734500409899e+05, new int[] { 4, 7 });
                p.AddCoeff(4.9475442450860787e+05, new int[] { 6, 5 });
                p.AddCoeff(1.5404441788729549e+05, new int[] { 8, 3 });
                p.AddCoeff(4.4346120300888096e+03, new int[] { 10, 1 });
                p.AddCoeff(-1.1125430180749119e+05, new int[] { 4, 9 });
                p.AddCoeff(-7.0679203501229696e+05, new int[] { 6, 7 });
                p.AddCoeff(-6.0077322976045241e+05, new int[] { 8, 5 });
                p.AddCoeff(-6.5040976441302540e+04, new int[] { 10, 3 });
                p.AddCoeff(3.3376290542247356e+05, new int[] { 6, 9 });
                p.AddCoeff(8.5824747108636059e+05, new int[] { 8, 7 });
                p.AddCoeff(2.5365980812107991e+05, new int[] { 10, 5 });
                p.AddCoeff(-4.0528352801300361e+05, new int[] { 8, 9 });
                p.AddCoeff(-3.6237115445868558e+05, new int[] { 10, 7 });
                p.AddCoeff(1.7111971182771264e+05, new int[] { 10, 9 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{50468800-053F-4E73-9F48-266D847788DD}"));
                OrthonormalPolynomials[201] = p;
                p.AddCoeff(-7.3182918850987587e+00, new int[] { 1, 0 });
                p.AddCoeff(2.6345850786355531e+02, new int[] { 1, 2 });
                p.AddCoeff(1.5856299084380644e+02, new int[] { 3, 0 });
                p.AddCoeff(-1.4490217932495542e+03, new int[] { 1, 4 });
                p.AddCoeff(-5.7082676703770318e+03, new int[] { 3, 2 });
                p.AddCoeff(-9.5137794506283863e+02, new int[] { 5, 0 });
                p.AddCoeff(2.5116377749658940e+03, new int[] { 1, 6 });
                p.AddCoeff(3.1395472187073675e+04, new int[] { 3, 4 });
                p.AddCoeff(3.4249606022262191e+04, new int[] { 5, 2 });
                p.AddCoeff(2.3104892951526081e+03, new int[] { 7, 0 });
                p.AddCoeff(-1.3455202365888718e+03, new int[] { 1, 8 });
                p.AddCoeff(-5.4418818457594369e+04, new int[] { 3, 6 });
                p.AddCoeff(-1.8837283312244205e+05, new int[] { 5, 4 });
                p.AddCoeff(-8.3177614625493891e+04, new int[] { 7, 2 });
                p.AddCoeff(-2.4388498115499752e+03, new int[] { 9, 0 });
                p.AddCoeff(2.9152938459425555e+04, new int[] { 3, 8 });
                p.AddCoeff(3.2651291074556622e+05, new int[] { 5, 6 });
                p.AddCoeff(4.5747688044021640e+05, new int[] { 7, 4 });
                p.AddCoeff(8.7798593215799107e+04, new int[] { 9, 2 });
                p.AddCoeff(9.3119720077362690e+02, new int[] { 11, 0 });
                p.AddCoeff(-1.7491763075655333e+05, new int[] { 5, 8 });
                p.AddCoeff(-7.9295992609637510e+05, new int[] { 7, 6 });
                p.AddCoeff(-4.8289226268689509e+05, new int[] { 9, 4 });
                p.AddCoeff(-3.3523099227850568e+04, new int[] { 11, 2 });
                p.AddCoeff(4.2479996040877237e+05, new int[] { 7, 8 });
                p.AddCoeff(8.3701325532395149e+05, new int[] { 9, 6 });
                p.AddCoeff(1.8437704575317813e+05, new int[] { 11, 4 });
                p.AddCoeff(-4.4839995820925973e+05, new int[] { 9, 8 });
                p.AddCoeff(-3.1958687930550875e+05, new int[] { 11, 6 });
                p.AddCoeff(1.7120725677080826e+05, new int[] { 11, 8 });
                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{2F66FEBE-8061-40E9-9721-CA395F2CAE33}"));
                OrthonormalPolynomials[202] = p;
                p.AddCoeff(-4.7779953543223519e+00, new int[] { 0, 1 });
                p.AddCoeff(4.3001958188901167e+01, new int[] { 0, 3 });
                p.AddCoeff(3.7268363763714345e+02, new int[] { 2, 1 });
                p.AddCoeff(-9.4604308015582568e+01, new int[] { 0, 5 });
                p.AddCoeff(-3.3541527387342910e+03, new int[] { 2, 3 });
                p.AddCoeff(-4.6585454704642931e+03, new int[] { 4, 1 });
                p.AddCoeff(5.8564571628693971e+01, new int[] { 0, 7 });
                p.AddCoeff(7.3791360252154403e+03, new int[] { 2, 5 });
                p.AddCoeff(4.1926909234178638e+04, new int[] { 4, 3 });
                p.AddCoeff(2.1118739466104795e+04, new int[] { 6, 1 });
                p.AddCoeff(-4.5680365870381297e+03, new int[] { 2, 7 });
                p.AddCoeff(-9.2239200315193004e+04, new int[] { 4, 5 });
                p.AddCoeff(-1.9006865519494316e+05, new int[] { 6, 3 });
                p.AddCoeff(-4.2991719627427619e+04, new int[] { 8, 1 });
                p.AddCoeff(5.7100457337976621e+04, new int[] { 4, 7 });
                p.AddCoeff(4.1815104142887495e+05, new int[] { 6, 5 });
                p.AddCoeff(3.8692547664684857e+05, new int[] { 8, 3 });
                p.AddCoeff(4.0125604985599111e+04, new int[] { 10, 1 });
                p.AddCoeff(-2.5885540659882735e+05, new int[] { 6, 7 });
                p.AddCoeff(-8.5123604862306686e+05, new int[] { 8, 5 });
                p.AddCoeff(-3.6113044487039200e+05, new int[] { 10, 3 });
                p.AddCoeff(-1.3983165373769387e+04, new int[] { 12, 1 });
                p.AddCoeff(5.2695564914761282e+05, new int[] { 8, 7 });
                p.AddCoeff(7.9448697871486241e+05, new int[] { 10, 5 });
                p.AddCoeff(1.2584848836392449e+05, new int[] { 12, 3 });
                p.AddCoeff(-4.9182527253777197e+05, new int[] { 10, 7 });
                p.AddCoeff(-2.7686667440063387e+05, new int[] { 12, 5 });
                p.AddCoeff(1.7139365558134478e+05, new int[] { 12, 7 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{29A32116-EFFC-46BE-BB6C-CB193628648F}"));
                OrthonormalPolynomials[203] = p;
                p.AddCoeff(-8.5847758434404250e+00, new int[] { 1, 0 });
                p.AddCoeff(1.8028029271224893e+02, new int[] { 1, 2 });
                p.AddCoeff(2.5754327530321275e+02, new int[] { 3, 0 });
                p.AddCoeff(-5.4084087813674678e+02, new int[] { 1, 4 });
                p.AddCoeff(-5.4084087813674678e+03, new int[] { 3, 2 });
                p.AddCoeff(-2.1891178400773084e+03, new int[] { 5, 0 });
                p.AddCoeff(3.9661664396694764e+02, new int[] { 1, 6 });
                p.AddCoeff(1.6225226344102403e+04, new int[] { 3, 4 });
                p.AddCoeff(4.5971474641623476e+04, new int[] { 5, 2 });
                p.AddCoeff(7.9225217069464494e+03, new int[] { 7, 0 });
                p.AddCoeff(-1.1898499319008429e+04, new int[] { 3, 6 });
                p.AddCoeff(-1.3791442392487043e+05, new int[] { 5, 4 });
                p.AddCoeff(-1.6637295584587544e+05, new int[] { 7, 2 });
                p.AddCoeff(-1.3864412987156286e+04, new int[] { 9, 0 });
                p.AddCoeff(1.0113724421157165e+05, new int[] { 5, 6 });
                p.AddCoeff(4.9911886753762631e+05, new int[] { 7, 4 });
                p.AddCoeff(2.9115267273028201e+05, new int[] { 9, 2 });
                p.AddCoeff(1.1595690861985258e+04, new int[] { 11, 0 });
                p.AddCoeff(-3.6602050286092596e+05, new int[] { 7, 6 });
                p.AddCoeff(-8.7345801819084604e+05, new int[] { 9, 4 });
                p.AddCoeff(-2.4350950810169041e+05, new int[] { 11, 2 });
                p.AddCoeff(-3.7165675839696339e+03, new int[] { 13, 0 });
                p.AddCoeff(6.4053588000662043e+05, new int[] { 9, 6 });
                p.AddCoeff(7.3052852430507124e+05, new int[] { 11, 4 });
                p.AddCoeff(7.8047919263362312e+04, new int[] { 13, 2 });
                p.AddCoeff(-5.3572091782371891e+05, new int[] { 11, 6 });
                p.AddCoeff(-2.3414375779008693e+05, new int[] { 13, 4 });
                p.AddCoeff(1.7170542237939709e+05, new int[] { 13, 6 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{C6561BF4-9023-40FE-B237-231D82145EED}"));
                OrthonormalPolynomials[204] = p;
                p.AddCoeff(-3.5074699409554877e+00, new int[] { 0, 1 });
                p.AddCoeff(1.6368193057792276e+01, new int[] { 0, 3 });
                p.AddCoeff(3.6828434380032621e+02, new int[] { 2, 1 });
                p.AddCoeff(-1.4731373752013048e+01, new int[] { 0, 5 });
                p.AddCoeff(-1.7186602710681890e+03, new int[] { 2, 3 });
                p.AddCoeff(-6.2608338446055456e+03, new int[] { 4, 1 });
                p.AddCoeff(1.5467942439613701e+03, new int[] { 2, 5 });
                p.AddCoeff(2.9217224608159213e+04, new int[] { 4, 3 });
                p.AddCoeff(3.9651947682501789e+04, new int[] { 6, 1 });
                p.AddCoeff(-2.6295502147343292e+04, new int[] { 4, 5 });
                p.AddCoeff(-1.8504242251834168e+05, new int[] { 6, 3 });
                p.AddCoeff(-1.1895584304750537e+05, new int[] { 8, 1 });
                p.AddCoeff(1.6653818026650751e+05, new int[] { 6, 5 });
                p.AddCoeff(5.5512726755502504e+05, new int[] { 8, 3 });
                p.AddCoeff(1.8239895933950823e+05, new int[] { 10, 1 });
                p.AddCoeff(-4.9961454079952254e+05, new int[] { 8, 5 });
                p.AddCoeff(-8.5119514358437173e+05, new int[] { 10, 3 });
                p.AddCoeff(-1.3818102980265775e+05, new int[] { 12, 1 });
                p.AddCoeff(7.6607562922593456e+05, new int[] { 10, 5 });
                p.AddCoeff(6.4484480574573616e+05, new int[] { 12, 3 });
                p.AddCoeff(4.0998767084305046e+04, new int[] { 14, 1 });
                p.AddCoeff(-5.8036032517116255e+05, new int[] { 12, 5 });
                p.AddCoeff(-1.9132757972675688e+05, new int[] { 14, 3 });
                p.AddCoeff(1.7219482175408119e+05, new int[] { 14, 5 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{79A2578B-B717-4120-9773-BFF97AAFBCC7}"));
                OrthonormalPolynomials[205] = p;
                p.AddCoeff(-9.8406089194732886e+00, new int[] { 1, 0 });
                p.AddCoeff(9.8406089194732886e+01, new int[] { 1, 2 });
                p.AddCoeff(3.9034415380577378e+02, new int[] { 3, 0 });
                p.AddCoeff(-1.1480710406052170e+02, new int[] { 1, 4 });
                p.AddCoeff(-3.9034415380577378e+03, new int[] { 3, 2 });
                p.AddCoeff(-4.4499233533858211e+03, new int[] { 5, 0 });
                p.AddCoeff(4.5540151277340274e+03, new int[] { 3, 4 });
                p.AddCoeff(4.4499233533858211e+04, new int[] { 5, 2 });
                p.AddCoeff(2.2249616766929105e+04, new int[] { 7, 0 });
                p.AddCoeff(-5.1915772456167913e+04, new int[] { 5, 4 });
                p.AddCoeff(-2.2249616766929105e+05, new int[] { 7, 2 });
                p.AddCoeff(-5.6860131737707714e+04, new int[] { 9, 0 });
                p.AddCoeff(2.5957886228083956e+05, new int[] { 7, 4 });
                p.AddCoeff(5.6860131737707714e+05, new int[] { 9, 2 });
                p.AddCoeff(7.7536543278692337e+04, new int[] { 11, 0 });
                p.AddCoeff(-6.6336820360659000e+05, new int[] { 9, 4 });
                p.AddCoeff(-7.7536543278692337e+05, new int[] { 11, 2 });
                p.AddCoeff(-5.3679145346787003e+04, new int[] { 13, 0 });
                p.AddCoeff(9.0459300491807727e+05, new int[] { 11, 4 });
                p.AddCoeff(5.3679145346787003e+05, new int[] { 13, 2 });
                p.AddCoeff(1.4825668714826886e+04, new int[] { 15, 0 });
                p.AddCoeff(-6.2625669571251503e+05, new int[] { 13, 4 });
                p.AddCoeff(-1.4825668714826886e+05, new int[] { 15, 2 });
                p.AddCoeff(1.7296613500631368e+05, new int[] { 15, 4 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{948D0809-A1A6-4E31-80FA-991F2F4DFB0B}"));
                OrthonormalPolynomials[206] = p;
                p.AddCoeff(-2.2385452086233647e+00, new int[] { 0, 1 });
                p.AddCoeff(3.7309086810389412e+00, new int[] { 0, 3 });
                p.AddCoeff(3.0444214837277760e+02, new int[] { 2, 1 });
                p.AddCoeff(-5.0740358062129600e+02, new int[] { 2, 3 });
                p.AddCoeff(-6.7484676222632369e+03, new int[] { 4, 1 });
                p.AddCoeff(1.1247446037105395e+04, new int[] { 4, 3 });
                p.AddCoeff(5.6687128027011190e+04, new int[] { 6, 1 });
                p.AddCoeff(-9.4478546711685316e+04, new int[] { 6, 3 });
                p.AddCoeff(-2.3282213296808167e+05, new int[] { 8, 1 });
                p.AddCoeff(3.8803688828013612e+05, new int[] { 8, 3 });
                p.AddCoeff(5.1738251770684816e+05, new int[] { 10, 1 });
                p.AddCoeff(-8.6230419617808027e+05, new int[] { 10, 3 });
                p.AddCoeff(-6.3496945354931365e+05, new int[] { 12, 1 });
                p.AddCoeff(1.0582824225821894e+06, new int[] { 12, 3 });
                p.AddCoeff(4.0470580555890321e+05, new int[] { 14, 1 });
                p.AddCoeff(-6.7450967593150534e+05, new int[] { 14, 3 });
                p.AddCoeff(-1.0454899976938333e+05, new int[] { 16, 1 });
                p.AddCoeff(1.7424833294897221e+05, new int[] { 16, 3 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{A2D162C9-6D04-418F-869B-430385DAF98B}"));
                OrthonormalPolynomials[207] = p;
                p.AddCoeff(-1.1040953242260395e+01, new int[] { 1, 0 });
                p.AddCoeff(3.3122859726781184e+01, new int[] { 1, 2 });
                p.AddCoeff(5.5940829760785999e+02, new int[] { 3, 0 });
                p.AddCoeff(-1.6782248928235800e+03, new int[] { 3, 2 });
                p.AddCoeff(-8.2233019748355419e+03, new int[] { 5, 0 });
                p.AddCoeff(2.4669905924506626e+04, new int[] { 5, 2 });
                p.AddCoeff(5.4038841548919275e+04, new int[] { 7, 0 });
                p.AddCoeff(-1.6211652464675783e+05, new int[] { 7, 2 });
                p.AddCoeff(-1.8763486648930304e+05, new int[] { 9, 0 });
                p.AddCoeff(5.6290459946790912e+05, new int[] { 9, 2 });
                p.AddCoeff(3.6844664692444960e+05, new int[] { 11, 0 });
                p.AddCoeff(-1.1053399407733488e+06, new int[] { 11, 2 });
                p.AddCoeff(-4.1095972156957840e+05, new int[] { 13, 0 });
                p.AddCoeff(1.2328791647087352e+06, new int[] { 13, 2 });
                p.AddCoeff(2.4266193083156058e+05, new int[] { 15, 0 });
                p.AddCoeff(-7.2798579249468174e+05, new int[] { 15, 2 });
                p.AddCoeff(-5.8881203804716906e+04, new int[] { 17, 0 });
                p.AddCoeff(1.7664361141415072e+05, new int[] { 17, 2 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{58F1E8B9-C9D6-4203-9794-8E99400F7B93}"));
                OrthonormalPolynomials[208] = p;
                p.AddCoeff(-9.7702691166629822e-01, new int[] { 0, 1 });
                p.AddCoeff(1.6707160189493700e+02, new int[] { 2, 1 });
                p.AddCoeff(-4.6780048530582359e+03, new int[] { 4, 1 });
                p.AddCoeff(5.0210585422825065e+04, new int[] { 6, 1 });
                p.AddCoeff(-2.6898527905084856e+05, new int[] { 8, 1 });
                p.AddCoeff(8.0695583715254569e+05, new int[] { 10, 1 });
                p.AddCoeff(-1.4182860168135652e+06, new int[] { 12, 1 });
                p.AddCoeff(1.4494571380622149e+06, new int[] { 14, 1 });
                p.AddCoeff(-7.9720142593421822e+05, new int[] { 16, 1 });
                p.AddCoeff(1.8236633926599763e+05, new int[] { 18, 1 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{746345E4-49D7-4A24-9B1F-C9A12AD50928}"));
                OrthonormalPolynomials[209] = p;
                p.AddCoeff(-1.1003502370758957e+01, new int[] { 1, 0 });
                p.AddCoeff(6.9322064935781432e+02, new int[] { 3, 0 });
                p.AddCoeff(-1.2755259948183784e+04, new int[] { 5, 0 });
                p.AddCoeff(1.0629383290153153e+05, new int[] { 7, 0 });
                p.AddCoeff(-4.7832224805689188e+05, new int[] { 9, 0 });
                p.AddCoeff(1.2610313812408968e+06, new int[] { 11, 0 });
                p.AddCoeff(-2.0047165547932205e+06, new int[] { 13, 0 });
                p.AddCoeff(1.8901613230907508e+06, new int[] { 15, 0 });
                p.AddCoeff(-9.7287715159082760e+05, new int[] { 17, 0 });
                p.AddCoeff(2.1050558250795685e+05, new int[] { 19, 0 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{CEB03273-9549-4BFC-B39A-0221DDD0224A}"));
                OrthonormalPolynomials[210] = p;
                p.AddCoeff(5.6410580711896104e-01, new int[] { 0, 0 });
                p.AddCoeff(-1.1846221949498182e+02, new int[] { 0, 2 });
                p.AddCoeff(4.0869465725768728e+03, new int[] { 0, 4 });
                p.AddCoeff(-5.4492620967691637e+04, new int[] { 0, 6 });
                p.AddCoeff(3.6782519153191855e+05, new int[] { 0, 8 });
                p.AddCoeff(-1.4222574072567517e+06, new int[] { 0, 10 });
                p.AddCoeff(3.3401499715878260e+06, new int[] { 0, 12 });
                p.AddCoeff(-4.8450527060394839e+06, new int[] { 0, 14 });
                p.AddCoeff(4.2394211177845484e+06, new int[] { 0, 16 });
                p.AddCoeff(-2.0504389720003698e+06, new int[] { 0, 18 });
                p.AddCoeff(4.2087957846323380e+05, new int[] { 0, 20 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{7344050B-F8FB-40EE-81A0-8575B29692A6}"));
                OrthonormalPolynomials[211] = p;
                p.AddCoeff(-1.9058625167359108e+01, new int[] { 1, 1 });
                p.AddCoeff(1.2006933855436238e+03, new int[] { 1, 3 });
                p.AddCoeff(-2.2092758294002678e+04, new int[] { 1, 5 });
                p.AddCoeff(1.8410631911668898e+05, new int[] { 1, 7 });
                p.AddCoeff(-8.2847843602510043e+05, new int[] { 1, 9 });
                p.AddCoeff(2.1841704222479920e+06, new int[] { 1, 11 });
                p.AddCoeff(-3.4722709276762950e+06, new int[] { 1, 13 });
                p.AddCoeff(3.2738554460947925e+06, new int[] { 1, 15 });
                p.AddCoeff(-1.6850726560782020e+06, new int[] { 1, 17 });
                p.AddCoeff(3.6460636418066359e+05, new int[] { 1, 19 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{5EE2C832-BC10-47EF-84CA-469ECF0E5473}"));
                OrthonormalPolynomials[212] = p;
                p.AddCoeff(6.3066815961333967e-01, new int[] { 0, 0 });
                p.AddCoeff(-1.0784425529388108e+02, new int[] { 0, 2 });
                p.AddCoeff(-1.8920044788400190e+00, new int[] { 2, 0 });
                p.AddCoeff(3.0196391482286703e+03, new int[] { 0, 4 });
                p.AddCoeff(3.2353276588164325e+02, new int[] { 2, 2 });
                p.AddCoeff(-3.2410793524321062e+04, new int[] { 0, 6 });
                p.AddCoeff(-9.0589174446860110e+03, new int[] { 2, 4 });
                p.AddCoeff(1.7362925102314854e+05, new int[] { 0, 8 });
                p.AddCoeff(9.7232380572963185e+04, new int[] { 2, 6 });
                p.AddCoeff(-5.2088775306944563e+05, new int[] { 0, 10 });
                p.AddCoeff(-5.2088775306944563e+05, new int[] { 2, 8 });
                p.AddCoeff(9.1549968721296505e+05, new int[] { 0, 12 });
                p.AddCoeff(1.5626632592083369e+06, new int[] { 2, 10 });
                p.AddCoeff(-9.3562055945940384e+05, new int[] { 0, 14 });
                p.AddCoeff(-2.7464990616388951e+06, new int[] { 2, 12 });
                p.AddCoeff(5.1459130770267211e+05, new int[] { 0, 16 });
                p.AddCoeff(2.8068616783782115e+06, new int[] { 2, 14 });
                p.AddCoeff(-1.1771696581433676e+05, new int[] { 0, 18 });
                p.AddCoeff(-1.5437739231080163e+06, new int[] { 2, 16 });
                p.AddCoeff(3.5315089744301027e+05, new int[] { 2, 18 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{4F6A9D7C-8021-4538-BB82-AE434143555A}"));
                OrthonormalPolynomials[213] = p;
                p.AddCoeff(-3.9191496157610927e+01, new int[] { 1, 1 });
                p.AddCoeff(1.9857024719856203e+03, new int[] { 1, 3 });
                p.AddCoeff(6.5319160262684878e+01, new int[] { 3, 1 });
                p.AddCoeff(-2.9189826338188618e+04, new int[] { 1, 5 });
                p.AddCoeff(-3.3095041199760338e+03, new int[] { 3, 3 });
                p.AddCoeff(1.9181885879381092e+05, new int[] { 1, 7 });
                p.AddCoeff(4.8649710563647697e+04, new int[] { 3, 5 });
                p.AddCoeff(-6.6603770414517680e+05, new int[] { 1, 9 });
                p.AddCoeff(-3.1969809798968487e+05, new int[] { 3, 7 });
                p.AddCoeff(1.3078558554123472e+06, new int[] { 1, 11 });
                p.AddCoeff(1.1100628402419613e+06, new int[] { 3, 9 });
                p.AddCoeff(-1.4587623002676180e+06, new int[] { 1, 13 });
                p.AddCoeff(-2.1797597590205786e+06, new int[] { 3, 11 });
                p.AddCoeff(8.6136440587230777e+05, new int[] { 1, 15 });
                p.AddCoeff(2.4312705004460300e+06, new int[] { 3, 13 });
                p.AddCoeff(-2.0900753966019233e+05, new int[] { 1, 17 });
                p.AddCoeff(-1.4356073431205130e+06, new int[] { 3, 15 });
                p.AddCoeff(3.4834589943365388e+05, new int[] { 3, 17 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{982C2ED3-BBC2-41C4-BBE9-E079AE91BFAA}"));
                OrthonormalPolynomials[214] = p;
                p.AddCoeff(6.3456792006349550e-01, new int[] { 0, 0 });
                p.AddCoeff(-8.6301237128635388e+01, new int[] { 0, 2 });
                p.AddCoeff(-6.3456792006349550e+00, new int[] { 2, 0 });
                p.AddCoeff(1.9130107563514178e+03, new int[] { 0, 4 });
                p.AddCoeff(8.6301237128635388e+02, new int[] { 2, 2 });
                p.AddCoeff(7.4032924007407809e+00, new int[] { 4, 0 });
                p.AddCoeff(-1.6069290353351909e+04, new int[] { 0, 6 });
                p.AddCoeff(-1.9130107563514178e+04, new int[] { 2, 4 });
                p.AddCoeff(-1.0068477665007462e+03, new int[] { 4, 2 });
                p.AddCoeff(6.5998871094123913e+04, new int[] { 0, 8 });
                p.AddCoeff(1.6069290353351909e+05, new int[] { 2, 6 });
                p.AddCoeff(2.2318458824099874e+04, new int[] { 4, 4 });
                p.AddCoeff(-1.4666415798694203e+05, new int[] { 0, 10 });
                p.AddCoeff(-6.5998871094123913e+05, new int[] { 2, 8 });
                p.AddCoeff(-1.8747505412243894e+05, new int[] { 4, 6 });
                p.AddCoeff(1.7999692116579249e+05, new int[] { 0, 12 });
                p.AddCoeff(1.4666415798694203e+06, new int[] { 2, 10 });
                p.AddCoeff(7.6998682943144565e+05, new int[] { 4, 8 });
                p.AddCoeff(-1.1472331239138422e+05, new int[] { 0, 14 });
                p.AddCoeff(-1.7999692116579249e+06, new int[] { 2, 12 });
                p.AddCoeff(-1.7110818431809903e+06, new int[] { 4, 10 });
                p.AddCoeff(2.9636855701107591e+04, new int[] { 0, 16 });
                p.AddCoeff(1.1472331239138422e+06, new int[] { 2, 14 });
                p.AddCoeff(2.0999640802675791e+06, new int[] { 4, 12 });
                p.AddCoeff(-2.9636855701107591e+05, new int[] { 2, 16 });
                p.AddCoeff(-1.3384386445661493e+06, new int[] { 4, 14 });
                p.AddCoeff(3.4576331651292190e+05, new int[] { 4, 16 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{D6316433-319B-4A7A-87FA-4837DCCC83EA}"));
                OrthonormalPolynomials[215] = p;
                p.AddCoeff(-5.4396012490862623e+01, new int[] { 1, 1 });
                p.AddCoeff(2.1577084954708840e+03, new int[] { 1, 3 });
                p.AddCoeff(2.5384805829069224e+02, new int[] { 3, 1 });
                p.AddCoeff(-2.4597876848368078e+04, new int[] { 1, 5 });
                p.AddCoeff(-1.0069306312197459e+04, new int[] { 3, 3 });
                p.AddCoeff(-2.2846325246162301e+02, new int[] { 5, 1 });
                p.AddCoeff(1.2298938424184039e+05, new int[] { 1, 7 });
                p.AddCoeff(1.1479009195905103e+05, new int[] { 3, 5 });
                p.AddCoeff(9.0623756809777129e+03, new int[] { 5, 3 });
                p.AddCoeff(-3.1430620417359211e+05, new int[] { 1, 9 });
                p.AddCoeff(-5.7395045979525515e+05, new int[] { 3, 7 });
                p.AddCoeff(-1.0331108276314593e+05, new int[] { 5, 5 });
                p.AddCoeff(4.2859936932762560e+05, new int[] { 1, 11 });
                p.AddCoeff(1.4667622861434298e+06, new int[] { 3, 9 });
                p.AddCoeff(5.1655541381572964e+05, new int[] { 5, 7 });
                p.AddCoeff(-2.9672264030374080e+05, new int[] { 1, 13 });
                p.AddCoeff(-2.0001303901955861e+06, new int[] { 3, 11 });
                p.AddCoeff(-1.3200860575290868e+06, new int[] { 5, 9 });
                p.AddCoeff(8.1951967321985554e+04, new int[] { 1, 15 });
                p.AddCoeff(1.3847056547507904e+06, new int[] { 3, 13 });
                p.AddCoeff(1.8001173511760275e+06, new int[] { 5, 11 });
                p.AddCoeff(-3.8244251416926592e+05, new int[] { 3, 15 });
                p.AddCoeff(-1.2462350892757114e+06, new int[] { 5, 13 });
                p.AddCoeff(3.4419826275233933e+05, new int[] { 5, 15 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{7034D419-F82D-41E8-B6F3-9897B81ABDD4}"));
                OrthonormalPolynomials[216] = p;
                p.AddCoeff(6.3550363791721496e-01, new int[] { 0, 0 });
                p.AddCoeff(-6.6727881981307571e+01, new int[] { 0, 2 });
                p.AddCoeff(-1.3345576396261514e+01, new int[] { 2, 0 });
                p.AddCoeff(1.1343739936822287e+03, new int[] { 0, 4 });
                p.AddCoeff(1.4012855216074590e+03, new int[] { 2, 2 });
                p.AddCoeff(4.0036729188784543e+01, new int[] { 4, 0 });
                p.AddCoeff(-7.1843686266541151e+03, new int[] { 0, 6 });
                p.AddCoeff(-2.3821853867326803e+04, new int[] { 2, 4 });
                p.AddCoeff(-4.2038565648223770e+03, new int[] { 4, 2 });
                p.AddCoeff(-2.9360268071775331e+01, new int[] { 6, 0 });
                p.AddCoeff(2.1553105879962345e+04, new int[] { 0, 8 });
                p.AddCoeff(1.5087174115973642e+05, new int[] { 2, 6 });
                p.AddCoeff(7.1465561601980408e+04, new int[] { 4, 4 });
                p.AddCoeff(3.0828281475364098e+03, new int[] { 6, 2 });
                p.AddCoeff(-3.3048095682608930e+04, new int[] { 0, 10 });
                p.AddCoeff(-4.5261522347920925e+05, new int[] { 2, 8 });
                p.AddCoeff(-4.5261522347920925e+05, new int[] { 4, 6 });
                p.AddCoeff(-5.2408078508118966e+04, new int[] { 6, 4 });
                p.AddCoeff(2.5036436123188583e+04, new int[] { 0, 12 });
                p.AddCoeff(6.9401000933478752e+05, new int[] { 2, 10 });
                p.AddCoeff(1.3578456704376278e+06, new int[] { 4, 8 });
                p.AddCoeff(3.3191783055142012e+05, new int[] { 6, 6 });
                p.AddCoeff(-7.4283931354515576e+03, new int[] { 0, 14 });
                p.AddCoeff(-5.2576515858696024e+05, new int[] { 2, 12 });
                p.AddCoeff(-2.0820300280043626e+06, new int[] { 4, 10 });
                p.AddCoeff(-9.9575349165426036e+05, new int[] { 6, 8 });
                p.AddCoeff(1.5599625584448271e+05, new int[] { 2, 14 });
                p.AddCoeff(1.5772954757608807e+06, new int[] { 4, 12 });
                p.AddCoeff(1.5268220205365325e+06, new int[] { 6, 10 });
                p.AddCoeff(-4.6798876753344813e+05, new int[] { 4, 14 });
                p.AddCoeff(-1.1566833488913125e+06, new int[] { 6, 12 });
                p.AddCoeff(3.4319176285786196e+05, new int[] { 6, 14 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{EA99F77F-E0F4-493E-8168-9957C12B447C}"));
                OrthonormalPolynomials[217] = p;
                p.AddCoeff(-6.4550699553712114e+01, new int[] { 1, 1 });
                p.AddCoeff(1.9365209866113634e+03, new int[] { 1, 3 });
                p.AddCoeff(5.8095629598340903e+02, new int[] { 3, 1 });
                p.AddCoeff(-1.6460428386196589e+04, new int[] { 1, 5 });
                p.AddCoeff(-1.7428688879502271e+04, new int[] { 3, 3 });
                p.AddCoeff(-1.2781038511634999e+03, new int[] { 5, 1 });
                p.AddCoeff(5.9571074159568608e+04, new int[] { 1, 7 });
                p.AddCoeff(1.4814385547576930e+05, new int[] { 3, 5 });
                p.AddCoeff(3.8343115534904996e+04, new int[] { 5, 3 });
                p.AddCoeff(7.9120714595835706e+02, new int[] { 7, 1 });
                p.AddCoeff(-1.0424937977924506e+05, new int[] { 1, 9 });
                p.AddCoeff(-5.3613966743611748e+05, new int[] { 3, 7 });
                p.AddCoeff(-3.2591648204669247e+05, new int[] { 5, 5 });
                p.AddCoeff(-2.3736214378750712e+04, new int[] { 7, 3 });
                p.AddCoeff(8.7190390360823145e+04, new int[] { 1, 11 });
                p.AddCoeff(9.3824441801320558e+05, new int[] { 3, 9 });
                p.AddCoeff(1.1795072683594584e+06, new int[] { 5, 7 });
                p.AddCoeff(2.0175782221938105e+05, new int[] { 7, 5 });
                p.AddCoeff(-2.7945637936161264e+04, new int[] { 1, 13 });
                p.AddCoeff(-7.8471351324740831e+05, new int[] { 3, 11 });
                p.AddCoeff(-2.0641377196290523e+06, new int[] { 5, 9 });
                p.AddCoeff(-7.3017116612728380e+05, new int[] { 7, 7 });
                p.AddCoeff(2.5151074142545138e+05, new int[] { 3, 13 });
                p.AddCoeff(1.7263697291442983e+06, new int[] { 5, 11 });
                p.AddCoeff(1.2777995407227467e+06, new int[] { 7, 9 });
                p.AddCoeff(-5.5332363113599304e+05, new int[] { 5, 13 });
                p.AddCoeff(-1.0687050704226608e+06, new int[] { 7, 11 });
                p.AddCoeff(3.4253367641751950e+05, new int[] { 7, 13 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{7E4BEC70-F40B-4627-8863-7AB26426D4D0}"));
                OrthonormalPolynomials[218] = p;
                p.AddCoeff(6.3582056013333865e-01, new int[] { 0, 0 });
                p.AddCoeff(-4.9594003690400415e+01, new int[] { 0, 2 });
                p.AddCoeff(-2.2889540164800191e+01, new int[] { 2, 0 });
                p.AddCoeff(6.1992504613000519e+02, new int[] { 0, 4 });
                p.AddCoeff(1.7853841328544149e+03, new int[] { 2, 2 });
                p.AddCoeff(1.2589247090640105e+02, new int[] { 4, 0 });
                p.AddCoeff(-2.8103268757893568e+03, new int[] { 0, 6 });
                p.AddCoeff(-2.2317301660680187e+04, new int[] { 2, 4 });
                p.AddCoeff(-9.8196127306992821e+03, new int[] { 4, 2 });
                p.AddCoeff(-2.1821361623776183e+02, new int[] { 6, 0 });
                p.AddCoeff(5.7210225685711907e+03, new int[] { 0, 8 });
                p.AddCoeff(1.0117176752841685e+05, new int[] { 2, 6 });
                p.AddCoeff(1.2274515913374103e+05, new int[] { 4, 4 });
                p.AddCoeff(1.7020662066545422e+04, new int[] { 6, 2 });
                p.AddCoeff(1.1690015155594384e+02, new int[] { 8, 0 });
                p.AddCoeff(-5.3396210639997780e+03, new int[] { 0, 10 });
                p.AddCoeff(-2.0595681246856287e+05, new int[] { 2, 8 });
                p.AddCoeff(-5.5644472140629266e+05, new int[] { 4, 6 });
                p.AddCoeff(-2.1275827583181778e+05, new int[] { 6, 4 });
                p.AddCoeff(-9.1182118213636191e+03, new int[] { 8, 2 });
                p.AddCoeff(1.8607770374544681e+03, new int[] { 0, 12 });
                p.AddCoeff(1.9222635830399201e+05, new int[] { 2, 10 });
                p.AddCoeff(1.1327624685770958e+06, new int[] { 4, 8 });
                p.AddCoeff(9.6450418377090727e+05, new int[] { 6, 6 });
                p.AddCoeff(1.1397764776704524e+05, new int[] { 8, 4 });
                p.AddCoeff(-6.6987973348360851e+04, new int[] { 2, 12 });
                p.AddCoeff(-1.0572449706719560e+06, new int[] { 4, 10 });
                p.AddCoeff(-1.9634549455336327e+06, new int[] { 6, 8 });
                p.AddCoeff(-5.1669866987727175e+05, new int[] { 8, 6 });
                p.AddCoeff(3.6843385341598468e+05, new int[] { 4, 12 });
                p.AddCoeff(1.8325579491647238e+06, new int[] { 6, 10 });
                p.AddCoeff(1.0518508636787318e+06, new int[] { 8, 8 });
                p.AddCoeff(-6.3861867925437345e+05, new int[] { 6, 12 });
                p.AddCoeff(-9.8172747276681633e+05, new int[] { 8, 10 });
                p.AddCoeff(3.4211714960055720e+05, new int[] { 8, 12 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{CF957A02-0ACC-41D6-B6DB-059B056C26AC}"));
                OrthonormalPolynomials[219] = p;
                p.AddCoeff(-6.9631311677906592e+01, new int[] { 1, 1 });
                p.AddCoeff(1.5086784196879762e+03, new int[] { 1, 3 });
                p.AddCoeff(1.0212592379426300e+03, new int[] { 3, 1 });
                p.AddCoeff(-9.0520705181278569e+03, new int[] { 1, 5 });
                p.AddCoeff(-2.2127283488756984e+04, new int[] { 3, 3 });
                p.AddCoeff(-3.9829110279762570e+03, new int[] { 5, 1 });
                p.AddCoeff(2.1983599829739081e+04, new int[] { 1, 7 });
                p.AddCoeff(1.3276370093254190e+05, new int[] { 3, 5 });
                p.AddCoeff(8.6296405606152236e+04, new int[] { 5, 3 });
                p.AddCoeff(5.6898728971089386e+03, new int[] { 7, 1 });
                p.AddCoeff(-2.3204910931391252e+04, new int[] { 1, 9 });
                p.AddCoeff(-3.2242613083617319e+05, new int[] { 3, 7 });
                p.AddCoeff(-5.1777843363691342e+05, new int[] { 5, 5 });
                p.AddCoeff(-1.2328057943736034e+05, new int[] { 7, 3 });
                p.AddCoeff(-2.6868844236347766e+03, new int[] { 9, 1 });
                p.AddCoeff(8.8600569010766600e+03, new int[] { 1, 11 });
                p.AddCoeff(3.4033869366040503e+05, new int[] { 3, 9 });
                p.AddCoeff(1.2574619102610754e+06, new int[] { 5, 7 });
                p.AddCoeff(7.3968347662416202e+05, new int[] { 7, 5 });
                p.AddCoeff(5.8215829178753493e+04, new int[] { 9, 3 });
                p.AddCoeff(-1.2994750121579101e+05, new int[] { 3, 11 });
                p.AddCoeff(-1.3273209052755796e+06, new int[] { 5, 9 });
                p.AddCoeff(-1.7963741575158221e+06, new int[] { 7, 7 });
                p.AddCoeff(-3.4929497507252096e+05, new int[] { 9, 5 });
                p.AddCoeff(5.0679525474158495e+05, new int[] { 5, 11 });
                p.AddCoeff(1.8961727218222566e+06, new int[] { 7, 9 });
                p.AddCoeff(8.4828779660469375e+05, new int[] { 9, 7 });
                p.AddCoeff(-7.2399322105940707e+05, new int[] { 7, 11 });
                p.AddCoeff(-8.9541489641606562e+05, new int[] { 9, 9 });
                p.AddCoeff(3.4188568772249778e+05, new int[] { 9, 11 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{70CB0917-6A77-47D2-995F-5F8A4E579DA0}"));
                OrthonormalPolynomials[220] = p;
                p.AddCoeff(6.3590240478515625e-01, new int[] { 0, 0 });
                p.AddCoeff(-3.4974632263183594e+01, new int[] { 0, 2 });
                p.AddCoeff(-3.4974632263183594e+01, new int[] { 2, 0 });
                p.AddCoeff(3.0311347961425781e+02, new int[] { 0, 4 });
                p.AddCoeff(1.9236047744750977e+03, new int[] { 2, 2 });
                p.AddCoeff(3.0311347961425781e+02, new int[] { 4, 0 });
                p.AddCoeff(-9.0934043884277344e+02, new int[] { 0, 6 });
                p.AddCoeff(-1.6671241378784180e+04, new int[] { 2, 4 });
                p.AddCoeff(-1.6671241378784180e+04, new int[] { 4, 2 });
                p.AddCoeff(-9.0934043884277344e+02, new int[] { 6, 0 });
                p.AddCoeff(1.1041991043090820e+03, new int[] { 0, 8 });
                p.AddCoeff(5.0013724136352539e+04, new int[] { 2, 6 });
                p.AddCoeff(1.4448409194946289e+05, new int[] { 4, 4 });
                p.AddCoeff(5.0013724136352539e+04, new int[] { 6, 2 });
                p.AddCoeff(1.1041991043090820e+03, new int[] { 8, 0 });
                p.AddCoeff(-4.6621739959716797e+02, new int[] { 0, 10 });
                p.AddCoeff(-6.0730950736999512e+04, new int[] { 2, 8 });
                p.AddCoeff(-4.3345227584838867e+05, new int[] { 4, 6 });
                p.AddCoeff(-4.3345227584838867e+05, new int[] { 6, 4 });
                p.AddCoeff(-6.0730950736999512e+04, new int[] { 8, 2 });
                p.AddCoeff(-4.6621739959716797e+02, new int[] { 10, 0 });
                p.AddCoeff(2.5641956977844238e+04, new int[] { 2, 10 });
                p.AddCoeff(5.2633490638732910e+05, new int[] { 4, 8 });
                p.AddCoeff(1.3003568275451660e+06, new int[] { 6, 6 });
                p.AddCoeff(5.2633490638732910e+05, new int[] { 8, 4 });
                p.AddCoeff(2.5641956977844238e+04, new int[] { 10, 2 });
                p.AddCoeff(-2.2223029380798340e+05, new int[] { 4, 10 });
                p.AddCoeff(-1.5790047191619873e+06, new int[] { 6, 8 });
                p.AddCoeff(-1.5790047191619873e+06, new int[] { 8, 6 });
                p.AddCoeff(-2.2223029380798340e+05, new int[] { 10, 4 });
                p.AddCoeff(6.6669088142395020e+05, new int[] { 6, 10 });
                p.AddCoeff(1.9173628732681274e+06, new int[] { 8, 8 });
                p.AddCoeff(6.6669088142395020e+05, new int[] { 10, 6 });
                p.AddCoeff(-8.0955321315765381e+05, new int[] { 8, 10 });
                p.AddCoeff(-8.0955321315765381e+05, new int[] { 10, 8 });
                p.AddCoeff(3.4181135666656494e+05, new int[] { 10, 10 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{7759EB66-2688-4218-9DFB-1F4ADC2ACD07}"));
                OrthonormalPolynomials[221] = p;
                p.AddCoeff(-6.9631311677906592e+01, new int[] { 1, 1 });
                p.AddCoeff(1.0212592379426300e+03, new int[] { 1, 3 });
                p.AddCoeff(1.5086784196879762e+03, new int[] { 3, 1 });
                p.AddCoeff(-3.9829110279762570e+03, new int[] { 1, 5 });
                p.AddCoeff(-2.2127283488756984e+04, new int[] { 3, 3 });
                p.AddCoeff(-9.0520705181278569e+03, new int[] { 5, 1 });
                p.AddCoeff(5.6898728971089386e+03, new int[] { 1, 7 });
                p.AddCoeff(8.6296405606152236e+04, new int[] { 3, 5 });
                p.AddCoeff(1.3276370093254190e+05, new int[] { 5, 3 });
                p.AddCoeff(2.1983599829739081e+04, new int[] { 7, 1 });
                p.AddCoeff(-2.6868844236347766e+03, new int[] { 1, 9 });
                p.AddCoeff(-1.2328057943736034e+05, new int[] { 3, 7 });
                p.AddCoeff(-5.1777843363691342e+05, new int[] { 5, 5 });
                p.AddCoeff(-3.2242613083617319e+05, new int[] { 7, 3 });
                p.AddCoeff(-2.3204910931391252e+04, new int[] { 9, 1 });
                p.AddCoeff(5.8215829178753493e+04, new int[] { 3, 9 });
                p.AddCoeff(7.3968347662416202e+05, new int[] { 5, 7 });
                p.AddCoeff(1.2574619102610754e+06, new int[] { 7, 5 });
                p.AddCoeff(3.4033869366040503e+05, new int[] { 9, 3 });
                p.AddCoeff(8.8600569010766600e+03, new int[] { 11, 1 });
                p.AddCoeff(-3.4929497507252096e+05, new int[] { 5, 9 });
                p.AddCoeff(-1.7963741575158221e+06, new int[] { 7, 7 });
                p.AddCoeff(-1.3273209052755796e+06, new int[] { 9, 5 });
                p.AddCoeff(-1.2994750121579101e+05, new int[] { 11, 3 });
                p.AddCoeff(8.4828779660469375e+05, new int[] { 7, 9 });
                p.AddCoeff(1.8961727218222566e+06, new int[] { 9, 7 });
                p.AddCoeff(5.0679525474158495e+05, new int[] { 11, 5 });
                p.AddCoeff(-8.9541489641606562e+05, new int[] { 9, 9 });
                p.AddCoeff(-7.2399322105940707e+05, new int[] { 11, 7 });
                p.AddCoeff(3.4188568772249778e+05, new int[] { 11, 9 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{437A514E-0FE1-4D95-905C-419EA09D70EF}"));
                OrthonormalPolynomials[222] = p;
                p.AddCoeff(6.3582056013333865e-01, new int[] { 0, 0 });
                p.AddCoeff(-2.2889540164800191e+01, new int[] { 0, 2 });
                p.AddCoeff(-4.9594003690400415e+01, new int[] { 2, 0 });
                p.AddCoeff(1.2589247090640105e+02, new int[] { 0, 4 });
                p.AddCoeff(1.7853841328544149e+03, new int[] { 2, 2 });
                p.AddCoeff(6.1992504613000519e+02, new int[] { 4, 0 });
                p.AddCoeff(-2.1821361623776183e+02, new int[] { 0, 6 });
                p.AddCoeff(-9.8196127306992821e+03, new int[] { 2, 4 });
                p.AddCoeff(-2.2317301660680187e+04, new int[] { 4, 2 });
                p.AddCoeff(-2.8103268757893568e+03, new int[] { 6, 0 });
                p.AddCoeff(1.1690015155594384e+02, new int[] { 0, 8 });
                p.AddCoeff(1.7020662066545422e+04, new int[] { 2, 6 });
                p.AddCoeff(1.2274515913374103e+05, new int[] { 4, 4 });
                p.AddCoeff(1.0117176752841685e+05, new int[] { 6, 2 });
                p.AddCoeff(5.7210225685711907e+03, new int[] { 8, 0 });
                p.AddCoeff(-9.1182118213636191e+03, new int[] { 2, 8 });
                p.AddCoeff(-2.1275827583181778e+05, new int[] { 4, 6 });
                p.AddCoeff(-5.5644472140629266e+05, new int[] { 6, 4 });
                p.AddCoeff(-2.0595681246856287e+05, new int[] { 8, 2 });
                p.AddCoeff(-5.3396210639997780e+03, new int[] { 10, 0 });
                p.AddCoeff(1.1397764776704524e+05, new int[] { 4, 8 });
                p.AddCoeff(9.6450418377090727e+05, new int[] { 6, 6 });
                p.AddCoeff(1.1327624685770958e+06, new int[] { 8, 4 });
                p.AddCoeff(1.9222635830399201e+05, new int[] { 10, 2 });
                p.AddCoeff(1.8607770374544681e+03, new int[] { 12, 0 });
                p.AddCoeff(-5.1669866987727175e+05, new int[] { 6, 8 });
                p.AddCoeff(-1.9634549455336327e+06, new int[] { 8, 6 });
                p.AddCoeff(-1.0572449706719560e+06, new int[] { 10, 4 });
                p.AddCoeff(-6.6987973348360851e+04, new int[] { 12, 2 });
                p.AddCoeff(1.0518508636787318e+06, new int[] { 8, 8 });
                p.AddCoeff(1.8325579491647238e+06, new int[] { 10, 6 });
                p.AddCoeff(3.6843385341598468e+05, new int[] { 12, 4 });
                p.AddCoeff(-9.8172747276681633e+05, new int[] { 10, 8 });
                p.AddCoeff(-6.3861867925437345e+05, new int[] { 12, 6 });
                p.AddCoeff(3.4211714960055720e+05, new int[] { 12, 8 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{F3227A85-C52A-44E0-A984-7217DE1257BC}"));
                OrthonormalPolynomials[223] = p;
                p.AddCoeff(-6.4550699553712114e+01, new int[] { 1, 1 });
                p.AddCoeff(5.8095629598340903e+02, new int[] { 1, 3 });
                p.AddCoeff(1.9365209866113634e+03, new int[] { 3, 1 });
                p.AddCoeff(-1.2781038511634999e+03, new int[] { 1, 5 });
                p.AddCoeff(-1.7428688879502271e+04, new int[] { 3, 3 });
                p.AddCoeff(-1.6460428386196589e+04, new int[] { 5, 1 });
                p.AddCoeff(7.9120714595835706e+02, new int[] { 1, 7 });
                p.AddCoeff(3.8343115534904996e+04, new int[] { 3, 5 });
                p.AddCoeff(1.4814385547576930e+05, new int[] { 5, 3 });
                p.AddCoeff(5.9571074159568608e+04, new int[] { 7, 1 });
                p.AddCoeff(-2.3736214378750712e+04, new int[] { 3, 7 });
                p.AddCoeff(-3.2591648204669247e+05, new int[] { 5, 5 });
                p.AddCoeff(-5.3613966743611748e+05, new int[] { 7, 3 });
                p.AddCoeff(-1.0424937977924506e+05, new int[] { 9, 1 });
                p.AddCoeff(2.0175782221938105e+05, new int[] { 5, 7 });
                p.AddCoeff(1.1795072683594584e+06, new int[] { 7, 5 });
                p.AddCoeff(9.3824441801320558e+05, new int[] { 9, 3 });
                p.AddCoeff(8.7190390360823145e+04, new int[] { 11, 1 });
                p.AddCoeff(-7.3017116612728380e+05, new int[] { 7, 7 });
                p.AddCoeff(-2.0641377196290523e+06, new int[] { 9, 5 });
                p.AddCoeff(-7.8471351324740831e+05, new int[] { 11, 3 });
                p.AddCoeff(-2.7945637936161264e+04, new int[] { 13, 1 });
                p.AddCoeff(1.2777995407227467e+06, new int[] { 9, 7 });
                p.AddCoeff(1.7263697291442983e+06, new int[] { 11, 5 });
                p.AddCoeff(2.5151074142545138e+05, new int[] { 13, 3 });
                p.AddCoeff(-1.0687050704226608e+06, new int[] { 11, 7 });
                p.AddCoeff(-5.5332363113599304e+05, new int[] { 13, 5 });
                p.AddCoeff(3.4253367641751950e+05, new int[] { 13, 7 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{D5F1FE9E-4913-4948-ACB5-4735801CAB8C}"));
                OrthonormalPolynomials[224] = p;
                p.AddCoeff(6.3550363791721496e-01, new int[] { 0, 0 });
                p.AddCoeff(-1.3345576396261514e+01, new int[] { 0, 2 });
                p.AddCoeff(-6.6727881981307571e+01, new int[] { 2, 0 });
                p.AddCoeff(4.0036729188784543e+01, new int[] { 0, 4 });
                p.AddCoeff(1.4012855216074590e+03, new int[] { 2, 2 });
                p.AddCoeff(1.1343739936822287e+03, new int[] { 4, 0 });
                p.AddCoeff(-2.9360268071775331e+01, new int[] { 0, 6 });
                p.AddCoeff(-4.2038565648223770e+03, new int[] { 2, 4 });
                p.AddCoeff(-2.3821853867326803e+04, new int[] { 4, 2 });
                p.AddCoeff(-7.1843686266541151e+03, new int[] { 6, 0 });
                p.AddCoeff(3.0828281475364098e+03, new int[] { 2, 6 });
                p.AddCoeff(7.1465561601980408e+04, new int[] { 4, 4 });
                p.AddCoeff(1.5087174115973642e+05, new int[] { 6, 2 });
                p.AddCoeff(2.1553105879962345e+04, new int[] { 8, 0 });
                p.AddCoeff(-5.2408078508118966e+04, new int[] { 4, 6 });
                p.AddCoeff(-4.5261522347920925e+05, new int[] { 6, 4 });
                p.AddCoeff(-4.5261522347920925e+05, new int[] { 8, 2 });
                p.AddCoeff(-3.3048095682608930e+04, new int[] { 10, 0 });
                p.AddCoeff(3.3191783055142012e+05, new int[] { 6, 6 });
                p.AddCoeff(1.3578456704376278e+06, new int[] { 8, 4 });
                p.AddCoeff(6.9401000933478752e+05, new int[] { 10, 2 });
                p.AddCoeff(2.5036436123188583e+04, new int[] { 12, 0 });
                p.AddCoeff(-9.9575349165426036e+05, new int[] { 8, 6 });
                p.AddCoeff(-2.0820300280043626e+06, new int[] { 10, 4 });
                p.AddCoeff(-5.2576515858696024e+05, new int[] { 12, 2 });
                p.AddCoeff(-7.4283931354515576e+03, new int[] { 14, 0 });
                p.AddCoeff(1.5268220205365325e+06, new int[] { 10, 6 });
                p.AddCoeff(1.5772954757608807e+06, new int[] { 12, 4 });
                p.AddCoeff(1.5599625584448271e+05, new int[] { 14, 2 });
                p.AddCoeff(-1.1566833488913125e+06, new int[] { 12, 6 });
                p.AddCoeff(-4.6798876753344813e+05, new int[] { 14, 4 });
                p.AddCoeff(3.4319176285786196e+05, new int[] { 14, 6 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{2CFE9C79-9354-456F-BC24-78C94B0769DE}"));
                OrthonormalPolynomials[225] = p;
                p.AddCoeff(-5.4396012490862623e+01, new int[] { 1, 1 });
                p.AddCoeff(2.5384805829069224e+02, new int[] { 1, 3 });
                p.AddCoeff(2.1577084954708840e+03, new int[] { 3, 1 });
                p.AddCoeff(-2.2846325246162301e+02, new int[] { 1, 5 });
                p.AddCoeff(-1.0069306312197459e+04, new int[] { 3, 3 });
                p.AddCoeff(-2.4597876848368078e+04, new int[] { 5, 1 });
                p.AddCoeff(9.0623756809777129e+03, new int[] { 3, 5 });
                p.AddCoeff(1.1479009195905103e+05, new int[] { 5, 3 });
                p.AddCoeff(1.2298938424184039e+05, new int[] { 7, 1 });
                p.AddCoeff(-1.0331108276314593e+05, new int[] { 5, 5 });
                p.AddCoeff(-5.7395045979525515e+05, new int[] { 7, 3 });
                p.AddCoeff(-3.1430620417359211e+05, new int[] { 9, 1 });
                p.AddCoeff(5.1655541381572964e+05, new int[] { 7, 5 });
                p.AddCoeff(1.4667622861434298e+06, new int[] { 9, 3 });
                p.AddCoeff(4.2859936932762560e+05, new int[] { 11, 1 });
                p.AddCoeff(-1.3200860575290868e+06, new int[] { 9, 5 });
                p.AddCoeff(-2.0001303901955861e+06, new int[] { 11, 3 });
                p.AddCoeff(-2.9672264030374080e+05, new int[] { 13, 1 });
                p.AddCoeff(1.8001173511760275e+06, new int[] { 11, 5 });
                p.AddCoeff(1.3847056547507904e+06, new int[] { 13, 3 });
                p.AddCoeff(8.1951967321985554e+04, new int[] { 15, 1 });
                p.AddCoeff(-1.2462350892757114e+06, new int[] { 13, 5 });
                p.AddCoeff(-3.8244251416926592e+05, new int[] { 15, 3 });
                p.AddCoeff(3.4419826275233933e+05, new int[] { 15, 5 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{31BA5D8B-C7DA-4D74-A015-CEE99C6667FF}"));
                OrthonormalPolynomials[226] = p;
                p.AddCoeff(6.3456792006349550e-01, new int[] { 0, 0 });
                p.AddCoeff(-6.3456792006349550e+00, new int[] { 0, 2 });
                p.AddCoeff(-8.6301237128635388e+01, new int[] { 2, 0 });
                p.AddCoeff(7.4032924007407809e+00, new int[] { 0, 4 });
                p.AddCoeff(8.6301237128635388e+02, new int[] { 2, 2 });
                p.AddCoeff(1.9130107563514178e+03, new int[] { 4, 0 });
                p.AddCoeff(-1.0068477665007462e+03, new int[] { 2, 4 });
                p.AddCoeff(-1.9130107563514178e+04, new int[] { 4, 2 });
                p.AddCoeff(-1.6069290353351909e+04, new int[] { 6, 0 });
                p.AddCoeff(2.2318458824099874e+04, new int[] { 4, 4 });
                p.AddCoeff(1.6069290353351909e+05, new int[] { 6, 2 });
                p.AddCoeff(6.5998871094123913e+04, new int[] { 8, 0 });
                p.AddCoeff(-1.8747505412243894e+05, new int[] { 6, 4 });
                p.AddCoeff(-6.5998871094123913e+05, new int[] { 8, 2 });
                p.AddCoeff(-1.4666415798694203e+05, new int[] { 10, 0 });
                p.AddCoeff(7.6998682943144565e+05, new int[] { 8, 4 });
                p.AddCoeff(1.4666415798694203e+06, new int[] { 10, 2 });
                p.AddCoeff(1.7999692116579249e+05, new int[] { 12, 0 });
                p.AddCoeff(-1.7110818431809903e+06, new int[] { 10, 4 });
                p.AddCoeff(-1.7999692116579249e+06, new int[] { 12, 2 });
                p.AddCoeff(-1.1472331239138422e+05, new int[] { 14, 0 });
                p.AddCoeff(2.0999640802675791e+06, new int[] { 12, 4 });
                p.AddCoeff(1.1472331239138422e+06, new int[] { 14, 2 });
                p.AddCoeff(2.9636855701107591e+04, new int[] { 16, 0 });
                p.AddCoeff(-1.3384386445661493e+06, new int[] { 14, 4 });
                p.AddCoeff(-2.9636855701107591e+05, new int[] { 16, 2 });
                p.AddCoeff(3.4576331651292190e+05, new int[] { 16, 4 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{53ECEE92-3851-4FC2-BC22-5ABA35D739CE}"));
                OrthonormalPolynomials[227] = p;
                p.AddCoeff(-3.9191496157610927e+01, new int[] { 1, 1 });
                p.AddCoeff(6.5319160262684878e+01, new int[] { 1, 3 });
                p.AddCoeff(1.9857024719856203e+03, new int[] { 3, 1 });
                p.AddCoeff(-3.3095041199760338e+03, new int[] { 3, 3 });
                p.AddCoeff(-2.9189826338188618e+04, new int[] { 5, 1 });
                p.AddCoeff(4.8649710563647697e+04, new int[] { 5, 3 });
                p.AddCoeff(1.9181885879381092e+05, new int[] { 7, 1 });
                p.AddCoeff(-3.1969809798968487e+05, new int[] { 7, 3 });
                p.AddCoeff(-6.6603770414517680e+05, new int[] { 9, 1 });
                p.AddCoeff(1.1100628402419613e+06, new int[] { 9, 3 });
                p.AddCoeff(1.3078558554123472e+06, new int[] { 11, 1 });
                p.AddCoeff(-2.1797597590205786e+06, new int[] { 11, 3 });
                p.AddCoeff(-1.4587623002676180e+06, new int[] { 13, 1 });
                p.AddCoeff(2.4312705004460300e+06, new int[] { 13, 3 });
                p.AddCoeff(8.6136440587230777e+05, new int[] { 15, 1 });
                p.AddCoeff(-1.4356073431205130e+06, new int[] { 15, 3 });
                p.AddCoeff(-2.0900753966019233e+05, new int[] { 17, 1 });
                p.AddCoeff(3.4834589943365388e+05, new int[] { 17, 3 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{9BD4504A-E952-4643-BEE1-63C86E899F63}"));
                OrthonormalPolynomials[228] = p;
                p.AddCoeff(6.3066815961333967e-01, new int[] { 0, 0 });
                p.AddCoeff(-1.8920044788400190e+00, new int[] { 0, 2 });
                p.AddCoeff(-1.0784425529388108e+02, new int[] { 2, 0 });
                p.AddCoeff(3.2353276588164325e+02, new int[] { 2, 2 });
                p.AddCoeff(3.0196391482286703e+03, new int[] { 4, 0 });
                p.AddCoeff(-9.0589174446860110e+03, new int[] { 4, 2 });
                p.AddCoeff(-3.2410793524321062e+04, new int[] { 6, 0 });
                p.AddCoeff(9.7232380572963185e+04, new int[] { 6, 2 });
                p.AddCoeff(1.7362925102314854e+05, new int[] { 8, 0 });
                p.AddCoeff(-5.2088775306944563e+05, new int[] { 8, 2 });
                p.AddCoeff(-5.2088775306944563e+05, new int[] { 10, 0 });
                p.AddCoeff(1.5626632592083369e+06, new int[] { 10, 2 });
                p.AddCoeff(9.1549968721296505e+05, new int[] { 12, 0 });
                p.AddCoeff(-2.7464990616388951e+06, new int[] { 12, 2 });
                p.AddCoeff(-9.3562055945940384e+05, new int[] { 14, 0 });
                p.AddCoeff(2.8068616783782115e+06, new int[] { 14, 2 });
                p.AddCoeff(5.1459130770267211e+05, new int[] { 16, 0 });
                p.AddCoeff(-1.5437739231080163e+06, new int[] { 16, 2 });
                p.AddCoeff(-1.1771696581433676e+05, new int[] { 18, 0 });
                p.AddCoeff(3.5315089744301027e+05, new int[] { 18, 2 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{3685001D-173B-4C59-9B45-56260EC4282B}"));
                OrthonormalPolynomials[229] = p;
                p.AddCoeff(-1.9058625167359108e+01, new int[] { 1, 1 });
                p.AddCoeff(1.2006933855436238e+03, new int[] { 3, 1 });
                p.AddCoeff(-2.2092758294002678e+04, new int[] { 5, 1 });
                p.AddCoeff(1.8410631911668898e+05, new int[] { 7, 1 });
                p.AddCoeff(-8.2847843602510043e+05, new int[] { 9, 1 });
                p.AddCoeff(2.1841704222479920e+06, new int[] { 11, 1 });
                p.AddCoeff(-3.4722709276762950e+06, new int[] { 13, 1 });
                p.AddCoeff(3.2738554460947925e+06, new int[] { 15, 1 });
                p.AddCoeff(-1.6850726560782020e+06, new int[] { 17, 1 });
                p.AddCoeff(3.6460636418066359e+05, new int[] { 19, 1 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{731E1926-3D6F-49D1-9EEA-007048856014}"));
                OrthonormalPolynomials[230] = p;
                p.AddCoeff(5.6410580711896104e-01, new int[] { 0, 0 });
                p.AddCoeff(-1.1846221949498182e+02, new int[] { 2, 0 });
                p.AddCoeff(4.0869465725768728e+03, new int[] { 4, 0 });
                p.AddCoeff(-5.4492620967691637e+04, new int[] { 6, 0 });
                p.AddCoeff(3.6782519153191855e+05, new int[] { 8, 0 });
                p.AddCoeff(-1.4222574072567517e+06, new int[] { 10, 0 });
                p.AddCoeff(3.3401499715878260e+06, new int[] { 12, 0 });
                p.AddCoeff(-4.8450527060394839e+06, new int[] { 14, 0 });
                p.AddCoeff(4.2394211177845484e+06, new int[] { 16, 0 });
                p.AddCoeff(-2.0504389720003698e+06, new int[] { 18, 0 });
                p.AddCoeff(4.2087957846323380e+05, new int[] { 20, 0 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{BC07BFC7-25F3-4E48-9EA8-9091AE57BA45}"));
                OrthonormalPolynomials[231] = p;
                p.AddCoeff(1.2131714034993529e+01, new int[] { 0, 1 });
                p.AddCoeff(-9.3009807601617055e+02, new int[] { 0, 3 });
                p.AddCoeff(2.0927206710363837e+04, new int[] { 0, 5 });
                p.AddCoeff(-2.1525126902088518e+05, new int[] { 0, 7 });
                p.AddCoeff(1.2137779892011026e+06, new int[] { 0, 9 });
                p.AddCoeff(-4.1047764725710014e+06, new int[] { 0, 11 });
                p.AddCoeff(8.6831809996694260e+06, new int[] { 0, 13 });
                p.AddCoeff(-1.1577574666225901e+07, new int[] { 0, 15 });
                p.AddCoeff(9.4493440290520225e+06, new int[] { 0, 17 });
                p.AddCoeff(-4.3102271009710980e+06, new int[] { 0, 19 });
                p.AddCoeff(8.4152052923721436e+05, new int[] { 0, 21 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{05F9F35E-3E21-466D-BEF3-8B884AD1B1E8}"));
                OrthonormalPolynomials[232] = p;
                p.AddCoeff(9.7705991877468981e-01, new int[] { 1, 0 });
                p.AddCoeff(-2.0518258294268486e+02, new int[] { 1, 2 });
                p.AddCoeff(7.0787991115226276e+03, new int[] { 1, 4 });
                p.AddCoeff(-9.4383988153635035e+04, new int[] { 1, 6 });
                p.AddCoeff(6.3709192003703649e+05, new int[] { 1, 8 });
                p.AddCoeff(-2.4634220908098744e+06, new int[] { 1, 10 });
                p.AddCoeff(5.7853094556898566e+06, new int[] { 1, 12 });
                p.AddCoeff(-8.3918774522094623e+06, new int[] { 1, 14 });
                p.AddCoeff(7.3428927706832795e+06, new int[] { 1, 16 });
                p.AddCoeff(-3.5514644773239391e+06, new int[] { 1, 18 });
                p.AddCoeff(7.2898481376649277e+05, new int[] { 1, 20 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{C0BD32B8-12C3-4EFB-BF12-C07191EC8E44}"));
                OrthonormalPolynomials[233] = p;
                p.AddCoeff(1.2302289645798562e+01, new int[] { 0, 1 });
                p.AddCoeff(-7.7504424768530938e+02, new int[] { 0, 3 });
                p.AddCoeff(-3.6906868937395685e+01, new int[] { 2, 1 });
                p.AddCoeff(1.4260814157409693e+04, new int[] { 0, 5 });
                p.AddCoeff(2.3251327430559281e+03, new int[] { 2, 3 });
                p.AddCoeff(-1.1884011797841410e+05, new int[] { 0, 7 });
                p.AddCoeff(-4.2782442472229078e+04, new int[] { 2, 5 });
                p.AddCoeff(5.3478053090286347e+05, new int[] { 0, 9 });
                p.AddCoeff(3.5652035393524231e+05, new int[] { 2, 7 });
                p.AddCoeff(-1.4098759451075491e+06, new int[] { 0, 11 });
                p.AddCoeff(-1.6043415927085904e+06, new int[] { 2, 9 });
                p.AddCoeff(2.2413412460684115e+06, new int[] { 0, 13 });
                p.AddCoeff(4.2296278353226474e+06, new int[] { 2, 11 });
                p.AddCoeff(-2.1132646034359308e+06, new int[] { 0, 15 });
                p.AddCoeff(-6.7240237382052344e+06, new int[] { 2, 13 });
                p.AddCoeff(1.0877097223567291e+06, new int[] { 0, 17 });
                p.AddCoeff(6.3397938103077924e+06, new int[] { 2, 15 });
                p.AddCoeff(-2.3535239606549109e+05, new int[] { 0, 19 });
                p.AddCoeff(-3.2631291670701873e+06, new int[] { 2, 17 });
                p.AddCoeff(7.0605718819647327e+05, new int[] { 2, 19 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{47E104C7-46C6-4976-BBF1-4CF76A7DC610}"));
                OrthonormalPolynomials[234] = p;
                p.AddCoeff(2.2386498893598723e+00, new int[] { 1, 0 });
                p.AddCoeff(-3.8280913108053817e+02, new int[] { 1, 2 });
                p.AddCoeff(-3.7310831489331206e+00, new int[] { 3, 0 });
                p.AddCoeff(1.0718655670255069e+04, new int[] { 1, 4 });
                p.AddCoeff(6.3801521846756362e+02, new int[] { 3, 2 });
                p.AddCoeff(-1.1504690419407107e+05, new int[] { 1, 6 });
                p.AddCoeff(-1.7864426117091781e+04, new int[] { 3, 4 });
                p.AddCoeff(6.1632270103966645e+05, new int[] { 1, 8 });
                p.AddCoeff(1.9174484032345179e+05, new int[] { 3, 6 });
                p.AddCoeff(-1.8489681031189994e+06, new int[] { 1, 10 });
                p.AddCoeff(-1.0272045017327774e+06, new int[] { 3, 8 });
                p.AddCoeff(3.2497015145727868e+06, new int[] { 1, 12 });
                p.AddCoeff(3.0816135051983323e+06, new int[] { 3, 10 });
                p.AddCoeff(-3.3211235258820788e+06, new int[] { 1, 14 });
                p.AddCoeff(-5.4161691909546446e+06, new int[] { 3, 12 });
                p.AddCoeff(1.8266179392351433e+06, new int[] { 1, 16 });
                p.AddCoeff(5.5352058764701313e+06, new int[] { 3, 14 });
                p.AddCoeff(-4.1785377694921579e+05, new int[] { 1, 18 });
                p.AddCoeff(-3.0443632320585722e+06, new int[] { 3, 16 });
                p.AddCoeff(6.9642296158202632e+05, new int[] { 3, 18 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{C213617E-9881-4D32-B909-5DB83AB39B81}"));
                OrthonormalPolynomials[235] = p;
                p.AddCoeff(1.1109744893740926e+01, new int[] { 0, 1 });
                p.AddCoeff(-5.6289374128287357e+02, new int[] { 0, 3 });
                p.AddCoeff(-1.1109744893740926e+02, new int[] { 2, 1 });
                p.AddCoeff(8.2745379968582415e+03, new int[] { 0, 5 });
                p.AddCoeff(5.6289374128287357e+03, new int[] { 2, 3 });
                p.AddCoeff(1.2961369042697747e+02, new int[] { 4, 1 });
                p.AddCoeff(-5.4375535407925587e+04, new int[] { 0, 7 });
                p.AddCoeff(-8.2745379968582415e+04, new int[] { 2, 5 });
                p.AddCoeff(-6.5670936483001917e+03, new int[] { 4, 3 });
                p.AddCoeff(1.8880394238863051e+05, new int[] { 0, 9 });
                p.AddCoeff(5.4375535407925587e+05, new int[] { 2, 7 });
                p.AddCoeff(9.6536276630012817e+04, new int[] { 4, 5 });
                p.AddCoeff(-3.7074228687221991e+05, new int[] { 0, 11 });
                p.AddCoeff(-1.8880394238863051e+06, new int[] { 2, 9 });
                p.AddCoeff(-6.3438124642579851e+05, new int[] { 4, 7 });
                p.AddCoeff(4.1352024304978375e+05, new int[] { 0, 13 });
                p.AddCoeff(3.7074228687221991e+06, new int[] { 2, 11 });
                p.AddCoeff(2.2027126612006893e+06, new int[] { 4, 9 });
                p.AddCoeff(-2.4417385780082469e+05, new int[] { 0, 15 });
                p.AddCoeff(-4.1352024304978375e+06, new int[] { 2, 13 });
                p.AddCoeff(-4.3253266801758990e+06, new int[] { 4, 11 });
                p.AddCoeff(5.9248068436964814e+04, new int[] { 0, 17 });
                p.AddCoeff(2.4417385780082469e+06, new int[] { 2, 15 });
                p.AddCoeff(4.8244028355808104e+06, new int[] { 4, 13 });
                p.AddCoeff(-5.9248068436964814e+05, new int[] { 2, 17 });
                p.AddCoeff(-2.8486950076762880e+06, new int[] { 4, 15 });
                p.AddCoeff(6.9122746509792283e+05, new int[] { 4, 17 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{48D22519-6570-41F6-BED7-D486EF512253}"));
                OrthonormalPolynomials[236] = p;
                p.AddCoeff(3.5077061580780882e+00, new int[] { 1, 0 });
                p.AddCoeff(-4.7704803749861999e+02, new int[] { 1, 2 });
                p.AddCoeff(-1.6369295404364411e+01, new int[] { 3, 0 });
                p.AddCoeff(1.0574564831219410e+04, new int[] { 1, 4 });
                p.AddCoeff(2.2262241749935600e+03, new int[] { 3, 2 });
                p.AddCoeff(1.4732365863927970e+01, new int[] { 5, 0 });
                p.AddCoeff(-8.8826344582243042e+04, new int[] { 1, 6 });
                p.AddCoeff(-4.9347969212357246e+04, new int[] { 3, 4 });
                p.AddCoeff(-2.0036017574942040e+03, new int[] { 5, 2 });
                p.AddCoeff(3.6482248667706964e+05, new int[] { 1, 8 });
                p.AddCoeff(4.1452294138380086e+05, new int[] { 3, 6 });
                p.AddCoeff(4.4413172291121521e+04, new int[] { 5, 4 });
                p.AddCoeff(-8.1071663706015475e+05, new int[] { 1, 10 });
                p.AddCoeff(-1.7025049378263250e+06, new int[] { 3, 8 });
                p.AddCoeff(-3.7307064724542078e+05, new int[] { 5, 6 });
                p.AddCoeff(9.9497041821018992e+05, new int[] { 1, 12 });
                p.AddCoeff(3.7833443062807222e+06, new int[] { 3, 10 });
                p.AddCoeff(1.5322544440436925e+06, new int[] { 5, 8 });
                p.AddCoeff(-6.3415696984825292e+05, new int[] { 1, 14 });
                p.AddCoeff(-4.6431952849808863e+06, new int[] { 3, 12 });
                p.AddCoeff(-3.4050098756526500e+06, new int[] { 5, 10 });
                p.AddCoeff(1.6382388387746534e+05, new int[] { 1, 16 });
                p.AddCoeff(2.9593991926251803e+06, new int[] { 3, 14 });
                p.AddCoeff(4.1788757564827977e+06, new int[] { 5, 12 });
                p.AddCoeff(-7.6451145809483824e+05, new int[] { 3, 16 });
                p.AddCoeff(-2.6634592733626623e+06, new int[] { 5, 14 });
                p.AddCoeff(6.8806031228535441e+05, new int[] { 5, 16 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{4B077613-E8F5-40D1-8FE1-EAE68F3E5B43}"));
                OrthonormalPolynomials[237] = p;
                p.AddCoeff(9.8557833447081181e+00, new int[] { 0, 1 });
                p.AddCoeff(-3.9094607267342202e+02, new int[] { 0, 3 });
                p.AddCoeff(-2.0697145023887048e+02, new int[] { 2, 1 });
                p.AddCoeff(4.4567852284770110e+03, new int[] { 0, 5 });
                p.AddCoeff(8.2098675261418624e+03, new int[] { 2, 3 });
                p.AddCoeff(6.2091435071661144e+02, new int[] { 4, 1 });
                p.AddCoeff(-2.2283926142385055e+04, new int[] { 0, 7 });
                p.AddCoeff(-9.3592489798017231e+04, new int[] { 2, 5 });
                p.AddCoeff(-2.4629602578425587e+04, new int[] { 4, 3 });
                p.AddCoeff(-4.5533719052551506e+02, new int[] { 6, 1 });
                p.AddCoeff(5.6947811252761807e+04, new int[] { 0, 9 });
                p.AddCoeff(4.6796244899008615e+05, new int[] { 2, 7 });
                p.AddCoeff(2.8077746939405169e+05, new int[] { 4, 5 });
                p.AddCoeff(1.8061708557512097e+04, new int[] { 6, 3 });
                p.AddCoeff(-7.7656106253766101e+04, new int[] { 0, 11 });
                p.AddCoeff(-1.1959040363079980e+06, new int[] { 2, 9 });
                p.AddCoeff(-1.4038873469702585e+06, new int[] { 4, 7 });
                p.AddCoeff(-2.0590347755563791e+05, new int[] { 6, 5 });
                p.AddCoeff(5.3761919714145762e+04, new int[] { 0, 13 });
                p.AddCoeff(1.6307782313290881e+06, new int[] { 2, 11 });
                p.AddCoeff(3.5877121089239939e+06, new int[] { 4, 9 });
                p.AddCoeff(1.0295173877781895e+06, new int[] { 6, 7 });
                p.AddCoeff(-1.4848530206764068e+04, new int[] { 0, 15 });
                p.AddCoeff(-1.1290003139970610e+06, new int[] { 2, 13 });
                p.AddCoeff(-4.8923346939872643e+06, new int[] { 4, 11 });
                p.AddCoeff(-2.6309888798775955e+06, new int[] { 6, 9 });
                p.AddCoeff(3.1181913434204542e+05, new int[] { 2, 15 });
                p.AddCoeff(3.3870009419911830e+06, new int[] { 4, 13 });
                p.AddCoeff(3.5877121089239939e+06, new int[] { 6, 11 });
                p.AddCoeff(-9.3545740302613626e+05, new int[] { 4, 15 });
                p.AddCoeff(-2.4838006907935342e+06, new int[] { 6, 13 });
                p.AddCoeff(6.8600209555249992e+05, new int[] { 6, 15 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{A303DB49-1CDE-4E2F-BB92-13CF942348A0}"));
                OrthonormalPolynomials[238] = p;
                p.AddCoeff(4.7784828799962214e+00, new int[] { 1, 0 });
                p.AddCoeff(-5.0174070239960325e+02, new int[] { 1, 2 });
                p.AddCoeff(-4.3006345919965993e+01, new int[] { 3, 0 });
                p.AddCoeff(8.5295919407932552e+03, new int[] { 1, 4 });
                p.AddCoeff(4.5156663215964292e+03, new int[] { 3, 2 });
                p.AddCoeff(9.4613961023925184e+01, new int[] { 5, 0 });
                p.AddCoeff(-5.4020748958357283e+04, new int[] { 1, 6 });
                p.AddCoeff(-7.6766327467139297e+04, new int[] { 3, 4 });
                p.AddCoeff(-9.9344659075121443e+03, new int[] { 5, 2 });
                p.AddCoeff(-5.8570547300525114e+01, new int[] { 7, 0 });
                p.AddCoeff(1.6206224687507185e+05, new int[] { 1, 8 });
                p.AddCoeff(4.8618674062521555e+05, new int[] { 3, 6 });
                p.AddCoeff(1.6888592042770645e+05, new int[] { 5, 4 });
                p.AddCoeff(6.1499074665551370e+03, new int[] { 7, 2 });
                p.AddCoeff(-2.4849544520844350e+05, new int[] { 1, 10 });
                p.AddCoeff(-1.4585602218756466e+06, new int[] { 3, 8 });
                p.AddCoeff(-1.0696108293754742e+06, new int[] { 5, 6 });
                p.AddCoeff(-1.0454842693143733e+05, new int[] { 7, 4 });
                p.AddCoeff(1.8825412515791174e+05, new int[] { 1, 12 });
                p.AddCoeff(2.2364590068759915e+06, new int[] { 3, 10 });
                p.AddCoeff(3.2088324881264226e+06, new int[] { 5, 8 });
                p.AddCoeff(6.6214003723243641e+05, new int[] { 7, 6 });
                p.AddCoeff(-5.5855619552347440e+04, new int[] { 1, 14 });
                p.AddCoeff(-1.6942871264212057e+06, new int[] { 3, 12 });
                p.AddCoeff(-4.9202098151271813e+06, new int[] { 5, 10 });
                p.AddCoeff(-1.9864201116973092e+06, new int[] { 7, 8 });
                p.AddCoeff(5.0270057597112696e+05, new int[] { 3, 14 });
                p.AddCoeff(3.7274316781266525e+06, new int[] { 5, 12 });
                p.AddCoeff(3.0458441712692075e+06, new int[] { 7, 10 });
                p.AddCoeff(-1.1059412671364793e+06, new int[] { 5, 14 });
                p.AddCoeff(-2.3074577055069754e+06, new int[] { 7, 12 });
                p.AddCoeff(6.8463030822734434e+05, new int[] { 7, 14 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{85ADC914-7063-432E-AD91-4B14992B7396}"));
                OrthonormalPolynomials[239] = p;
                p.AddCoeff(8.5899334142531919e+00, new int[] { 0, 1 });
                p.AddCoeff(-2.5769800242759576e+02, new int[] { 0, 3 });
                p.AddCoeff(-3.0923760291311491e+02, new int[] { 2, 1 });
                p.AddCoeff(2.1904330206345639e+03, new int[] { 0, 5 });
                p.AddCoeff(9.2771280873934473e+03, new int[] { 2, 3 });
                p.AddCoeff(1.7008068160221320e+03, new int[] { 4, 1 });
                p.AddCoeff(-7.9272814080108028e+03, new int[] { 0, 7 });
                p.AddCoeff(-7.8855588742844302e+04, new int[] { 2, 5 });
                p.AddCoeff(-5.1024204480663960e+04, new int[] { 4, 3 });
                p.AddCoeff(-2.9480651477716955e+03, new int[] { 6, 1 });
                p.AddCoeff(1.3872742464018905e+04, new int[] { 0, 9 });
                p.AddCoeff(2.8538213068838890e+05, new int[] { 2, 7 });
                p.AddCoeff(4.3370573808564366e+05, new int[] { 4, 5 });
                p.AddCoeff(8.8441954433150864e+04, new int[] { 6, 3 });
                p.AddCoeff(1.5793206148776940e+03, new int[] { 8, 1 });
                p.AddCoeff(-1.1602657333543084e+04, new int[] { 0, 11 });
                p.AddCoeff(-4.9941872870468058e+05, new int[] { 2, 9 });
                p.AddCoeff(-1.5696017187861390e+06, new int[] { 4, 7 });
                p.AddCoeff(-7.5175661268178235e+05, new int[] { 6, 5 });
                p.AddCoeff(-4.7379618446330820e+04, new int[] { 8, 3 });
                p.AddCoeff(3.7188004274176552e+03, new int[] { 0, 13 });
                p.AddCoeff(4.1769566400755103e+05, new int[] { 2, 11 });
                p.AddCoeff(2.7468030078757432e+06, new int[] { 4, 9 });
                p.AddCoeff(2.7206429792293075e+06, new int[] { 6, 7 });
                p.AddCoeff(4.0272675679381197e+05, new int[] { 8, 5 });
                p.AddCoeff(-1.3387681538703559e+05, new int[] { 2, 13 });
                p.AddCoeff(-2.2973261520415307e+06, new int[] { 4, 11 });
                p.AddCoeff(-4.7611252136512882e+06, new int[] { 6, 9 });
                p.AddCoeff(-1.4574873103014148e+06, new int[] { 8, 7 });
                p.AddCoeff(7.3632248462869573e+05, new int[] { 4, 13 });
                p.AddCoeff(3.9820319968719865e+06, new int[] { 6, 11 });
                p.AddCoeff(2.5506027930274758e+06, new int[] { 8, 9 });
                p.AddCoeff(-1.2762923066897393e+06, new int[] { 6, 13 });
                p.AddCoeff(-2.1332314268957070e+06, new int[] { 8, 11 });
                p.AddCoeff(6.8372802144093175e+05, new int[] { 8, 13 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{DD560F72-67DC-423A-90A6-5A348F2FA9D4}"));
                OrthonormalPolynomials[240] = p;
                p.AddCoeff(6.0496383977267668e+00, new int[] { 1, 0 });
                p.AddCoeff(-4.7187179502268781e+02, new int[] { 1, 2 });
                p.AddCoeff(-8.8728029833325913e+01, new int[] { 3, 0 });
                p.AddCoeff(5.8983974377835976e+03, new int[] { 1, 4 });
                p.AddCoeff(6.9207863269994212e+03, new int[] { 3, 2 });
                p.AddCoeff(3.4603931634997106e+02, new int[] { 5, 0 });
                p.AddCoeff(-2.6739401717952309e+04, new int[] { 1, 6 });
                p.AddCoeff(-8.6509829087492765e+04, new int[] { 3, 4 });
                p.AddCoeff(-2.6991066675297743e+04, new int[] { 5, 2 });
                p.AddCoeff(-4.9434188049995866e+02, new int[] { 7, 0 });
                p.AddCoeff(5.4433782068688630e+04, new int[] { 1, 8 });
                p.AddCoeff(3.9217789186330054e+05, new int[] { 3, 6 });
                p.AddCoeff(3.3738833344122178e+05, new int[] { 5, 4 });
                p.AddCoeff(3.8558666678996775e+04, new int[] { 7, 2 });
                p.AddCoeff(2.3343922134720270e+02, new int[] { 9, 0 });
                p.AddCoeff(-5.0804863264109388e+04, new int[] { 1, 10 });
                p.AddCoeff(-7.9836213700743323e+05, new int[] { 3, 8 });
                p.AddCoeff(-1.5294937782668721e+06, new int[] { 5, 6 });
                p.AddCoeff(-4.8198333348745969e+05, new int[] { 7, 4 });
                p.AddCoeff(-1.8208259265081811e+04, new int[] { 9, 2 });
                p.AddCoeff(1.7704725076886605e+04, new int[] { 1, 12 });
                p.AddCoeff(7.4513799454027102e+05, new int[] { 3, 10 });
                p.AddCoeff(3.1136123343289896e+06, new int[] { 5, 8 });
                p.AddCoeff(2.1849911118098173e+06, new int[] { 7, 6 });
                p.AddCoeff(2.2760324081352263e+05, new int[] { 9, 4 });
                p.AddCoeff(-2.5966930112767020e+05, new int[] { 3, 12 });
                p.AddCoeff(-2.9060381787070570e+06, new int[] { 5, 10 });
                p.AddCoeff(-4.4480176204699852e+06, new int[] { 7, 8 });
                p.AddCoeff(-1.0318013583546359e+06, new int[] { 9, 6 });
                p.AddCoeff(1.0127102743979138e+06, new int[] { 5, 12 });
                p.AddCoeff(4.1514831124386528e+06, new int[] { 7, 10 });
                p.AddCoeff(2.1004527652219374e+06, new int[] { 9, 8 });
                p.AddCoeff(-1.4467289634255911e+06, new int[] { 7, 12 });
                p.AddCoeff(-1.9604225808738083e+06, new int[] { 9, 10 });
                p.AddCoeff(6.8317756606208470e+05, new int[] { 9, 12 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{842A2429-2985-4221-B84C-019D050428D3}"));
                OrthonormalPolynomials[241] = p;
                p.AddCoeff(7.3204440074464418e+00, new int[] { 0, 1 });
                p.AddCoeff(-1.5860962016133957e+02, new int[] { 0, 3 });
                p.AddCoeff(-4.0262442040955430e+02, new int[] { 2, 1 });
                p.AddCoeff(9.5165772096803743e+02, new int[] { 0, 5 });
                p.AddCoeff(8.7235291088736765e+03, new int[] { 2, 3 });
                p.AddCoeff(3.4894116435494706e+03, new int[] { 4, 1 });
                p.AddCoeff(-2.3111687509223766e+03, new int[] { 0, 7 });
                p.AddCoeff(-5.2341174653242059e+04, new int[] { 2, 5 });
                p.AddCoeff(-7.5603918943571863e+04, new int[] { 4, 3 });
                p.AddCoeff(-1.0468234930648412e+04, new int[] { 6, 1 });
                p.AddCoeff(2.4395670148625087e+03, new int[] { 0, 9 });
                p.AddCoeff(1.2711428130073071e+05, new int[] { 2, 7 });
                p.AddCoeff(4.5362351366143118e+05, new int[] { 4, 5 });
                p.AddCoeff(2.2681175683071559e+05, new int[] { 6, 3 });
                p.AddCoeff(1.2711428130073071e+04, new int[] { 8, 1 });
                p.AddCoeff(-9.3147104203841240e+02, new int[] { 0, 11 });
                p.AddCoeff(-1.3417618581743798e+05, new int[] { 2, 9 });
                p.AddCoeff(-1.1016571046063329e+06, new int[] { 4, 7 });
                p.AddCoeff(-1.3608705409842935e+06, new int[] { 6, 5 });
                p.AddCoeff(-2.7541427615158321e+05, new int[] { 8, 3 });
                p.AddCoeff(-5.3670474326975191e+03, new int[] { 10, 1 });
                p.AddCoeff(5.1230907312112682e+04, new int[] { 2, 11 });
                p.AddCoeff(1.1628602770844625e+06, new int[] { 4, 9 });
                p.AddCoeff(3.3049713138189986e+06, new int[] { 6, 7 });
                p.AddCoeff(1.6524856569094993e+06, new int[] { 8, 5 });
                p.AddCoeff(1.1628602770844625e+05, new int[] { 10, 3 });
                p.AddCoeff(-4.4400119670497658e+05, new int[] { 4, 11 });
                p.AddCoeff(-3.4885808312533874e+06, new int[] { 6, 9 });
                p.AddCoeff(-4.0131794524944983e+06, new int[] { 8, 7 });
                p.AddCoeff(-6.9771616625067748e+05, new int[] { 10, 5 });
                p.AddCoeff(1.3320035901149297e+06, new int[] { 6, 11 });
                p.AddCoeff(4.2361338665219704e+06, new int[] { 8, 9 });
                p.AddCoeff(1.6944535466087882e+06, new int[] { 10, 7 });
                p.AddCoeff(-1.6174329308538432e+06, new int[] { 8, 11 });
                p.AddCoeff(-1.7885898547537208e+06, new int[] { 10, 9 });
                p.AddCoeff(6.8291612636051159e+05, new int[] { 10, 11 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{36FD7C00-4980-4C07-A723-810750FA2E1F}"));
                OrthonormalPolynomials[242] = p;
                p.AddCoeff(7.3204440074464418e+00, new int[] { 1, 0 });
                p.AddCoeff(-4.0262442040955430e+02, new int[] { 1, 2 });
                p.AddCoeff(-1.5860962016133957e+02, new int[] { 3, 0 });
                p.AddCoeff(3.4894116435494706e+03, new int[] { 1, 4 });
                p.AddCoeff(8.7235291088736765e+03, new int[] { 3, 2 });
                p.AddCoeff(9.5165772096803743e+02, new int[] { 5, 0 });
                p.AddCoeff(-1.0468234930648412e+04, new int[] { 1, 6 });
                p.AddCoeff(-7.5603918943571863e+04, new int[] { 3, 4 });
                p.AddCoeff(-5.2341174653242059e+04, new int[] { 5, 2 });
                p.AddCoeff(-2.3111687509223766e+03, new int[] { 7, 0 });
                p.AddCoeff(1.2711428130073071e+04, new int[] { 1, 8 });
                p.AddCoeff(2.2681175683071559e+05, new int[] { 3, 6 });
                p.AddCoeff(4.5362351366143118e+05, new int[] { 5, 4 });
                p.AddCoeff(1.2711428130073071e+05, new int[] { 7, 2 });
                p.AddCoeff(2.4395670148625087e+03, new int[] { 9, 0 });
                p.AddCoeff(-5.3670474326975191e+03, new int[] { 1, 10 });
                p.AddCoeff(-2.7541427615158321e+05, new int[] { 3, 8 });
                p.AddCoeff(-1.3608705409842935e+06, new int[] { 5, 6 });
                p.AddCoeff(-1.1016571046063329e+06, new int[] { 7, 4 });
                p.AddCoeff(-1.3417618581743798e+05, new int[] { 9, 2 });
                p.AddCoeff(-9.3147104203841240e+02, new int[] { 11, 0 });
                p.AddCoeff(1.1628602770844625e+05, new int[] { 3, 10 });
                p.AddCoeff(1.6524856569094993e+06, new int[] { 5, 8 });
                p.AddCoeff(3.3049713138189986e+06, new int[] { 7, 6 });
                p.AddCoeff(1.1628602770844625e+06, new int[] { 9, 4 });
                p.AddCoeff(5.1230907312112682e+04, new int[] { 11, 2 });
                p.AddCoeff(-6.9771616625067748e+05, new int[] { 5, 10 });
                p.AddCoeff(-4.0131794524944983e+06, new int[] { 7, 8 });
                p.AddCoeff(-3.4885808312533874e+06, new int[] { 9, 6 });
                p.AddCoeff(-4.4400119670497658e+05, new int[] { 11, 4 });
                p.AddCoeff(1.6944535466087882e+06, new int[] { 7, 10 });
                p.AddCoeff(4.2361338665219704e+06, new int[] { 9, 8 });
                p.AddCoeff(1.3320035901149297e+06, new int[] { 11, 6 });
                p.AddCoeff(-1.7885898547537208e+06, new int[] { 9, 10 });
                p.AddCoeff(-1.6174329308538432e+06, new int[] { 11, 8 });
                p.AddCoeff(6.8291612636051159e+05, new int[] { 11, 10 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{17D8A986-A5BD-450E-A449-4793F83B753B}"));
                OrthonormalPolynomials[243] = p;
                p.AddCoeff(6.0496383977267668e+00, new int[] { 0, 1 });
                p.AddCoeff(-8.8728029833325913e+01, new int[] { 0, 3 });
                p.AddCoeff(-4.7187179502268781e+02, new int[] { 2, 1 });
                p.AddCoeff(3.4603931634997106e+02, new int[] { 0, 5 });
                p.AddCoeff(6.9207863269994212e+03, new int[] { 2, 3 });
                p.AddCoeff(5.8983974377835976e+03, new int[] { 4, 1 });
                p.AddCoeff(-4.9434188049995866e+02, new int[] { 0, 7 });
                p.AddCoeff(-2.6991066675297743e+04, new int[] { 2, 5 });
                p.AddCoeff(-8.6509829087492765e+04, new int[] { 4, 3 });
                p.AddCoeff(-2.6739401717952309e+04, new int[] { 6, 1 });
                p.AddCoeff(2.3343922134720270e+02, new int[] { 0, 9 });
                p.AddCoeff(3.8558666678996775e+04, new int[] { 2, 7 });
                p.AddCoeff(3.3738833344122178e+05, new int[] { 4, 5 });
                p.AddCoeff(3.9217789186330054e+05, new int[] { 6, 3 });
                p.AddCoeff(5.4433782068688630e+04, new int[] { 8, 1 });
                p.AddCoeff(-1.8208259265081811e+04, new int[] { 2, 9 });
                p.AddCoeff(-4.8198333348745969e+05, new int[] { 4, 7 });
                p.AddCoeff(-1.5294937782668721e+06, new int[] { 6, 5 });
                p.AddCoeff(-7.9836213700743323e+05, new int[] { 8, 3 });
                p.AddCoeff(-5.0804863264109388e+04, new int[] { 10, 1 });
                p.AddCoeff(2.2760324081352263e+05, new int[] { 4, 9 });
                p.AddCoeff(2.1849911118098173e+06, new int[] { 6, 7 });
                p.AddCoeff(3.1136123343289896e+06, new int[] { 8, 5 });
                p.AddCoeff(7.4513799454027102e+05, new int[] { 10, 3 });
                p.AddCoeff(1.7704725076886605e+04, new int[] { 12, 1 });
                p.AddCoeff(-1.0318013583546359e+06, new int[] { 6, 9 });
                p.AddCoeff(-4.4480176204699852e+06, new int[] { 8, 7 });
                p.AddCoeff(-2.9060381787070570e+06, new int[] { 10, 5 });
                p.AddCoeff(-2.5966930112767020e+05, new int[] { 12, 3 });
                p.AddCoeff(2.1004527652219374e+06, new int[] { 8, 9 });
                p.AddCoeff(4.1514831124386528e+06, new int[] { 10, 7 });
                p.AddCoeff(1.0127102743979138e+06, new int[] { 12, 5 });
                p.AddCoeff(-1.9604225808738083e+06, new int[] { 10, 9 });
                p.AddCoeff(-1.4467289634255911e+06, new int[] { 12, 7 });
                p.AddCoeff(6.8317756606208470e+05, new int[] { 12, 9 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{34D5AE1F-FEC3-4443-8750-5DD00F85D3A4}"));
                OrthonormalPolynomials[244] = p;
                p.AddCoeff(8.5899334142531919e+00, new int[] { 1, 0 });
                p.AddCoeff(-3.0923760291311491e+02, new int[] { 1, 2 });
                p.AddCoeff(-2.5769800242759576e+02, new int[] { 3, 0 });
                p.AddCoeff(1.7008068160221320e+03, new int[] { 1, 4 });
                p.AddCoeff(9.2771280873934473e+03, new int[] { 3, 2 });
                p.AddCoeff(2.1904330206345639e+03, new int[] { 5, 0 });
                p.AddCoeff(-2.9480651477716955e+03, new int[] { 1, 6 });
                p.AddCoeff(-5.1024204480663960e+04, new int[] { 3, 4 });
                p.AddCoeff(-7.8855588742844302e+04, new int[] { 5, 2 });
                p.AddCoeff(-7.9272814080108028e+03, new int[] { 7, 0 });
                p.AddCoeff(1.5793206148776940e+03, new int[] { 1, 8 });
                p.AddCoeff(8.8441954433150864e+04, new int[] { 3, 6 });
                p.AddCoeff(4.3370573808564366e+05, new int[] { 5, 4 });
                p.AddCoeff(2.8538213068838890e+05, new int[] { 7, 2 });
                p.AddCoeff(1.3872742464018905e+04, new int[] { 9, 0 });
                p.AddCoeff(-4.7379618446330820e+04, new int[] { 3, 8 });
                p.AddCoeff(-7.5175661268178235e+05, new int[] { 5, 6 });
                p.AddCoeff(-1.5696017187861390e+06, new int[] { 7, 4 });
                p.AddCoeff(-4.9941872870468058e+05, new int[] { 9, 2 });
                p.AddCoeff(-1.1602657333543084e+04, new int[] { 11, 0 });
                p.AddCoeff(4.0272675679381197e+05, new int[] { 5, 8 });
                p.AddCoeff(2.7206429792293075e+06, new int[] { 7, 6 });
                p.AddCoeff(2.7468030078757432e+06, new int[] { 9, 4 });
                p.AddCoeff(4.1769566400755103e+05, new int[] { 11, 2 });
                p.AddCoeff(3.7188004274176552e+03, new int[] { 13, 0 });
                p.AddCoeff(-1.4574873103014148e+06, new int[] { 7, 8 });
                p.AddCoeff(-4.7611252136512882e+06, new int[] { 9, 6 });
                p.AddCoeff(-2.2973261520415307e+06, new int[] { 11, 4 });
                p.AddCoeff(-1.3387681538703559e+05, new int[] { 13, 2 });
                p.AddCoeff(2.5506027930274758e+06, new int[] { 9, 8 });
                p.AddCoeff(3.9820319968719865e+06, new int[] { 11, 6 });
                p.AddCoeff(7.3632248462869573e+05, new int[] { 13, 4 });
                p.AddCoeff(-2.1332314268957070e+06, new int[] { 11, 8 });
                p.AddCoeff(-1.2762923066897393e+06, new int[] { 13, 6 });
                p.AddCoeff(6.8372802144093175e+05, new int[] { 13, 8 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{594BC0E3-82EA-4784-A243-F3A9A9A9497F}"));
                OrthonormalPolynomials[245] = p;
                p.AddCoeff(4.7784828799962214e+00, new int[] { 0, 1 });
                p.AddCoeff(-4.3006345919965993e+01, new int[] { 0, 3 });
                p.AddCoeff(-5.0174070239960325e+02, new int[] { 2, 1 });
                p.AddCoeff(9.4613961023925184e+01, new int[] { 0, 5 });
                p.AddCoeff(4.5156663215964292e+03, new int[] { 2, 3 });
                p.AddCoeff(8.5295919407932552e+03, new int[] { 4, 1 });
                p.AddCoeff(-5.8570547300525114e+01, new int[] { 0, 7 });
                p.AddCoeff(-9.9344659075121443e+03, new int[] { 2, 5 });
                p.AddCoeff(-7.6766327467139297e+04, new int[] { 4, 3 });
                p.AddCoeff(-5.4020748958357283e+04, new int[] { 6, 1 });
                p.AddCoeff(6.1499074665551370e+03, new int[] { 2, 7 });
                p.AddCoeff(1.6888592042770645e+05, new int[] { 4, 5 });
                p.AddCoeff(4.8618674062521555e+05, new int[] { 6, 3 });
                p.AddCoeff(1.6206224687507185e+05, new int[] { 8, 1 });
                p.AddCoeff(-1.0454842693143733e+05, new int[] { 4, 7 });
                p.AddCoeff(-1.0696108293754742e+06, new int[] { 6, 5 });
                p.AddCoeff(-1.4585602218756466e+06, new int[] { 8, 3 });
                p.AddCoeff(-2.4849544520844350e+05, new int[] { 10, 1 });
                p.AddCoeff(6.6214003723243641e+05, new int[] { 6, 7 });
                p.AddCoeff(3.2088324881264226e+06, new int[] { 8, 5 });
                p.AddCoeff(2.2364590068759915e+06, new int[] { 10, 3 });
                p.AddCoeff(1.8825412515791174e+05, new int[] { 12, 1 });
                p.AddCoeff(-1.9864201116973092e+06, new int[] { 8, 7 });
                p.AddCoeff(-4.9202098151271813e+06, new int[] { 10, 5 });
                p.AddCoeff(-1.6942871264212057e+06, new int[] { 12, 3 });
                p.AddCoeff(-5.5855619552347440e+04, new int[] { 14, 1 });
                p.AddCoeff(3.0458441712692075e+06, new int[] { 10, 7 });
                p.AddCoeff(3.7274316781266525e+06, new int[] { 12, 5 });
                p.AddCoeff(5.0270057597112696e+05, new int[] { 14, 3 });
                p.AddCoeff(-2.3074577055069754e+06, new int[] { 12, 7 });
                p.AddCoeff(-1.1059412671364793e+06, new int[] { 14, 5 });
                p.AddCoeff(6.8463030822734434e+05, new int[] { 14, 7 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{68DF0C34-0A7B-4CCC-9FDA-00C784FB3075}"));
                OrthonormalPolynomials[246] = p;
                p.AddCoeff(9.8557833447081181e+00, new int[] { 1, 0 });
                p.AddCoeff(-2.0697145023887048e+02, new int[] { 1, 2 });
                p.AddCoeff(-3.9094607267342202e+02, new int[] { 3, 0 });
                p.AddCoeff(6.2091435071661144e+02, new int[] { 1, 4 });
                p.AddCoeff(8.2098675261418624e+03, new int[] { 3, 2 });
                p.AddCoeff(4.4567852284770110e+03, new int[] { 5, 0 });
                p.AddCoeff(-4.5533719052551506e+02, new int[] { 1, 6 });
                p.AddCoeff(-2.4629602578425587e+04, new int[] { 3, 4 });
                p.AddCoeff(-9.3592489798017231e+04, new int[] { 5, 2 });
                p.AddCoeff(-2.2283926142385055e+04, new int[] { 7, 0 });
                p.AddCoeff(1.8061708557512097e+04, new int[] { 3, 6 });
                p.AddCoeff(2.8077746939405169e+05, new int[] { 5, 4 });
                p.AddCoeff(4.6796244899008615e+05, new int[] { 7, 2 });
                p.AddCoeff(5.6947811252761807e+04, new int[] { 9, 0 });
                p.AddCoeff(-2.0590347755563791e+05, new int[] { 5, 6 });
                p.AddCoeff(-1.4038873469702585e+06, new int[] { 7, 4 });
                p.AddCoeff(-1.1959040363079980e+06, new int[] { 9, 2 });
                p.AddCoeff(-7.7656106253766101e+04, new int[] { 11, 0 });
                p.AddCoeff(1.0295173877781895e+06, new int[] { 7, 6 });
                p.AddCoeff(3.5877121089239939e+06, new int[] { 9, 4 });
                p.AddCoeff(1.6307782313290881e+06, new int[] { 11, 2 });
                p.AddCoeff(5.3761919714145762e+04, new int[] { 13, 0 });
                p.AddCoeff(-2.6309888798775955e+06, new int[] { 9, 6 });
                p.AddCoeff(-4.8923346939872643e+06, new int[] { 11, 4 });
                p.AddCoeff(-1.1290003139970610e+06, new int[] { 13, 2 });
                p.AddCoeff(-1.4848530206764068e+04, new int[] { 15, 0 });
                p.AddCoeff(3.5877121089239939e+06, new int[] { 11, 6 });
                p.AddCoeff(3.3870009419911830e+06, new int[] { 13, 4 });
                p.AddCoeff(3.1181913434204542e+05, new int[] { 15, 2 });
                p.AddCoeff(-2.4838006907935342e+06, new int[] { 13, 6 });
                p.AddCoeff(-9.3545740302613626e+05, new int[] { 15, 4 });
                p.AddCoeff(6.8600209555249992e+05, new int[] { 15, 6 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{507C63C4-CD7C-41A8-BEB6-287A92C37E72}"));
                OrthonormalPolynomials[247] = p;
                p.AddCoeff(3.5077061580780882e+00, new int[] { 0, 1 });
                p.AddCoeff(-1.6369295404364411e+01, new int[] { 0, 3 });
                p.AddCoeff(-4.7704803749861999e+02, new int[] { 2, 1 });
                p.AddCoeff(1.4732365863927970e+01, new int[] { 0, 5 });
                p.AddCoeff(2.2262241749935600e+03, new int[] { 2, 3 });
                p.AddCoeff(1.0574564831219410e+04, new int[] { 4, 1 });
                p.AddCoeff(-2.0036017574942040e+03, new int[] { 2, 5 });
                p.AddCoeff(-4.9347969212357246e+04, new int[] { 4, 3 });
                p.AddCoeff(-8.8826344582243042e+04, new int[] { 6, 1 });
                p.AddCoeff(4.4413172291121521e+04, new int[] { 4, 5 });
                p.AddCoeff(4.1452294138380086e+05, new int[] { 6, 3 });
                p.AddCoeff(3.6482248667706964e+05, new int[] { 8, 1 });
                p.AddCoeff(-3.7307064724542078e+05, new int[] { 6, 5 });
                p.AddCoeff(-1.7025049378263250e+06, new int[] { 8, 3 });
                p.AddCoeff(-8.1071663706015475e+05, new int[] { 10, 1 });
                p.AddCoeff(1.5322544440436925e+06, new int[] { 8, 5 });
                p.AddCoeff(3.7833443062807222e+06, new int[] { 10, 3 });
                p.AddCoeff(9.9497041821018992e+05, new int[] { 12, 1 });
                p.AddCoeff(-3.4050098756526500e+06, new int[] { 10, 5 });
                p.AddCoeff(-4.6431952849808863e+06, new int[] { 12, 3 });
                p.AddCoeff(-6.3415696984825292e+05, new int[] { 14, 1 });
                p.AddCoeff(4.1788757564827977e+06, new int[] { 12, 5 });
                p.AddCoeff(2.9593991926251803e+06, new int[] { 14, 3 });
                p.AddCoeff(1.6382388387746534e+05, new int[] { 16, 1 });
                p.AddCoeff(-2.6634592733626623e+06, new int[] { 14, 5 });
                p.AddCoeff(-7.6451145809483824e+05, new int[] { 16, 3 });
                p.AddCoeff(6.8806031228535441e+05, new int[] { 16, 5 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{0FAC08FA-3F88-43FD-A005-24EDDADA531E}"));
                OrthonormalPolynomials[248] = p;
                p.AddCoeff(1.1109744893740926e+01, new int[] { 1, 0 });
                p.AddCoeff(-1.1109744893740926e+02, new int[] { 1, 2 });
                p.AddCoeff(-5.6289374128287357e+02, new int[] { 3, 0 });
                p.AddCoeff(1.2961369042697747e+02, new int[] { 1, 4 });
                p.AddCoeff(5.6289374128287357e+03, new int[] { 3, 2 });
                p.AddCoeff(8.2745379968582415e+03, new int[] { 5, 0 });
                p.AddCoeff(-6.5670936483001917e+03, new int[] { 3, 4 });
                p.AddCoeff(-8.2745379968582415e+04, new int[] { 5, 2 });
                p.AddCoeff(-5.4375535407925587e+04, new int[] { 7, 0 });
                p.AddCoeff(9.6536276630012817e+04, new int[] { 5, 4 });
                p.AddCoeff(5.4375535407925587e+05, new int[] { 7, 2 });
                p.AddCoeff(1.8880394238863051e+05, new int[] { 9, 0 });
                p.AddCoeff(-6.3438124642579851e+05, new int[] { 7, 4 });
                p.AddCoeff(-1.8880394238863051e+06, new int[] { 9, 2 });
                p.AddCoeff(-3.7074228687221991e+05, new int[] { 11, 0 });
                p.AddCoeff(2.2027126612006893e+06, new int[] { 9, 4 });
                p.AddCoeff(3.7074228687221991e+06, new int[] { 11, 2 });
                p.AddCoeff(4.1352024304978375e+05, new int[] { 13, 0 });
                p.AddCoeff(-4.3253266801758990e+06, new int[] { 11, 4 });
                p.AddCoeff(-4.1352024304978375e+06, new int[] { 13, 2 });
                p.AddCoeff(-2.4417385780082469e+05, new int[] { 15, 0 });
                p.AddCoeff(4.8244028355808104e+06, new int[] { 13, 4 });
                p.AddCoeff(2.4417385780082469e+06, new int[] { 15, 2 });
                p.AddCoeff(5.9248068436964814e+04, new int[] { 17, 0 });
                p.AddCoeff(-2.8486950076762880e+06, new int[] { 15, 4 });
                p.AddCoeff(-5.9248068436964814e+05, new int[] { 17, 2 });
                p.AddCoeff(6.9122746509792283e+05, new int[] { 17, 4 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{BFE789E7-5773-4678-87D3-87242349BCAE}"));
                OrthonormalPolynomials[249] = p;
                p.AddCoeff(2.2386498893598723e+00, new int[] { 0, 1 });
                p.AddCoeff(-3.7310831489331206e+00, new int[] { 0, 3 });
                p.AddCoeff(-3.8280913108053817e+02, new int[] { 2, 1 });
                p.AddCoeff(6.3801521846756362e+02, new int[] { 2, 3 });
                p.AddCoeff(1.0718655670255069e+04, new int[] { 4, 1 });
                p.AddCoeff(-1.7864426117091781e+04, new int[] { 4, 3 });
                p.AddCoeff(-1.1504690419407107e+05, new int[] { 6, 1 });
                p.AddCoeff(1.9174484032345179e+05, new int[] { 6, 3 });
                p.AddCoeff(6.1632270103966645e+05, new int[] { 8, 1 });
                p.AddCoeff(-1.0272045017327774e+06, new int[] { 8, 3 });
                p.AddCoeff(-1.8489681031189994e+06, new int[] { 10, 1 });
                p.AddCoeff(3.0816135051983323e+06, new int[] { 10, 3 });
                p.AddCoeff(3.2497015145727868e+06, new int[] { 12, 1 });
                p.AddCoeff(-5.4161691909546446e+06, new int[] { 12, 3 });
                p.AddCoeff(-3.3211235258820788e+06, new int[] { 14, 1 });
                p.AddCoeff(5.5352058764701313e+06, new int[] { 14, 3 });
                p.AddCoeff(1.8266179392351433e+06, new int[] { 16, 1 });
                p.AddCoeff(-3.0443632320585722e+06, new int[] { 16, 3 });
                p.AddCoeff(-4.1785377694921579e+05, new int[] { 18, 1 });
                p.AddCoeff(6.9642296158202632e+05, new int[] { 18, 3 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{8F32062E-AFAC-480F-8F3F-9FA76B419C82}"));
                OrthonormalPolynomials[250] = p;
                p.AddCoeff(1.2302289645798562e+01, new int[] { 1, 0 });
                p.AddCoeff(-3.6906868937395685e+01, new int[] { 1, 2 });
                p.AddCoeff(-7.7504424768530938e+02, new int[] { 3, 0 });
                p.AddCoeff(2.3251327430559281e+03, new int[] { 3, 2 });
                p.AddCoeff(1.4260814157409693e+04, new int[] { 5, 0 });
                p.AddCoeff(-4.2782442472229078e+04, new int[] { 5, 2 });
                p.AddCoeff(-1.1884011797841410e+05, new int[] { 7, 0 });
                p.AddCoeff(3.5652035393524231e+05, new int[] { 7, 2 });
                p.AddCoeff(5.3478053090286347e+05, new int[] { 9, 0 });
                p.AddCoeff(-1.6043415927085904e+06, new int[] { 9, 2 });
                p.AddCoeff(-1.4098759451075491e+06, new int[] { 11, 0 });
                p.AddCoeff(4.2296278353226474e+06, new int[] { 11, 2 });
                p.AddCoeff(2.2413412460684115e+06, new int[] { 13, 0 });
                p.AddCoeff(-6.7240237382052344e+06, new int[] { 13, 2 });
                p.AddCoeff(-2.1132646034359308e+06, new int[] { 15, 0 });
                p.AddCoeff(6.3397938103077924e+06, new int[] { 15, 2 });
                p.AddCoeff(1.0877097223567291e+06, new int[] { 17, 0 });
                p.AddCoeff(-3.2631291670701873e+06, new int[] { 17, 2 });
                p.AddCoeff(-2.3535239606549109e+05, new int[] { 19, 0 });
                p.AddCoeff(7.0605718819647327e+05, new int[] { 19, 2 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{30E0202A-6576-4A70-925A-CB8AA9C40C9A}"));
                OrthonormalPolynomials[251] = p;
                p.AddCoeff(9.7705991877468981e-01, new int[] { 0, 1 });
                p.AddCoeff(-2.0518258294268486e+02, new int[] { 2, 1 });
                p.AddCoeff(7.0787991115226276e+03, new int[] { 4, 1 });
                p.AddCoeff(-9.4383988153635035e+04, new int[] { 6, 1 });
                p.AddCoeff(6.3709192003703649e+05, new int[] { 8, 1 });
                p.AddCoeff(-2.4634220908098744e+06, new int[] { 10, 1 });
                p.AddCoeff(5.7853094556898566e+06, new int[] { 12, 1 });
                p.AddCoeff(-8.3918774522094623e+06, new int[] { 14, 1 });
                p.AddCoeff(7.3428927706832795e+06, new int[] { 16, 1 });
                p.AddCoeff(-3.5514644773239391e+06, new int[] { 18, 1 });
                p.AddCoeff(7.2898481376649277e+05, new int[] { 20, 1 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{1A7FD837-435D-4721-BC38-290D4AD5898A}"));
                OrthonormalPolynomials[252] = p;
                p.AddCoeff(1.2131714034993529e+01, new int[] { 1, 0 });
                p.AddCoeff(-9.3009807601617055e+02, new int[] { 3, 0 });
                p.AddCoeff(2.0927206710363837e+04, new int[] { 5, 0 });
                p.AddCoeff(-2.1525126902088518e+05, new int[] { 7, 0 });
                p.AddCoeff(1.2137779892011026e+06, new int[] { 9, 0 });
                p.AddCoeff(-4.1047764725710014e+06, new int[] { 11, 0 });
                p.AddCoeff(8.6831809996694260e+06, new int[] { 13, 0 });
                p.AddCoeff(-1.1577574666225901e+07, new int[] { 15, 0 });
                p.AddCoeff(9.4493440290520225e+06, new int[] { 17, 0 });
                p.AddCoeff(-4.3102271009710980e+06, new int[] { 19, 0 });
                p.AddCoeff(8.4152052923721436e+05, new int[] { 21, 0 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{A6A35057-7B87-46DF-8567-A471BCA0ED7E}"));
                OrthonormalPolynomials[253] = p;
                p.AddCoeff(-5.6412002045046031e-01, new int[] { 0, 0 });
                p.AddCoeff(1.4272236517396646e+02, new int[] { 0, 2 });
                p.AddCoeff(-5.9467652155819357e+03, new int[] { 0, 4 });
                p.AddCoeff(9.6337596492427359e+04, new int[] { 0, 6 });
                p.AddCoeff(-7.9822579950868383e+05, new int[] { 0, 8 });
                p.AddCoeff(3.8492221887418754e+06, new int[] { 0, 10 });
                p.AddCoeff(-1.1547666566225626e+07, new int[] { 0, 12 });
                p.AddCoeff(2.2207051088895435e+07, new int[] { 0, 14 });
                p.AddCoeff(-2.7388696342971036e+07, new int[] { 0, 16 });
                p.AddCoeff(2.0944297203448439e+07, new int[] { 0, 18 });
                p.AddCoeff(-9.0391177404356423e+06, new int[] { 0, 20 });
                p.AddCoeff(1.6826063326352061e+06, new int[] { 0, 22 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{42079783-3AAE-496C-B8C8-4F90E2D7B8B3}"));
                OrthonormalPolynomials[254] = p;
                p.AddCoeff(2.1012745091505225e+01, new int[] { 1, 1 });
                p.AddCoeff(-1.6109771236820672e+03, new int[] { 1, 3 });
                p.AddCoeff(3.6246985282846512e+04, new int[] { 1, 5 });
                p.AddCoeff(-3.7282613433784984e+05, new int[] { 1, 7 });
                p.AddCoeff(2.1023251464050977e+06, new int[] { 1, 9 });
                p.AddCoeff(-7.1096814042063305e+06, new int[] { 1, 11 });
                p.AddCoeff(1.5039710662744161e+07, new int[] { 1, 13 });
                p.AddCoeff(-2.0052947550325547e+07, new int[] { 1, 15 });
                p.AddCoeff(1.6366743956515704e+07, new int[] { 1, 17 });
                p.AddCoeff(-7.4655323310422510e+06, new int[] { 1, 19 });
                p.AddCoeff(1.4575563122511062e+06, new int[] { 1, 21 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{C4C49B2F-72B4-452C-8E86-117FC941D807}"));
                OrthonormalPolynomials[255] = p;
                p.AddCoeff(-6.3068946561019085e-01, new int[] { 0, 0 });
                p.AddCoeff(1.3244478777814008e+02, new int[] { 0, 2 });
                p.AddCoeff(1.8920683968305725e+00, new int[] { 2, 0 });
                p.AddCoeff(-4.5693451783458327e+03, new int[] { 0, 4 });
                p.AddCoeff(-3.9733436333442023e+02, new int[] { 2, 2 });
                p.AddCoeff(6.0924602377944436e+04, new int[] { 0, 6 });
                p.AddCoeff(1.3708035535037498e+04, new int[] { 2, 4 });
                p.AddCoeff(-4.1124106605112494e+05, new int[] { 0, 8 });
                p.AddCoeff(-1.8277380713383331e+05, new int[] { 2, 6 });
                p.AddCoeff(1.5901321220643498e+06, new int[] { 0, 10 });
                p.AddCoeff(1.2337231981533748e+06, new int[] { 2, 8 });
                p.AddCoeff(-3.7344011957571851e+06, new int[] { 0, 12 });
                p.AddCoeff(-4.7703963661930493e+06, new int[] { 2, 10 });
                p.AddCoeff(5.4169336026367959e+06, new int[] { 0, 14 });
                p.AddCoeff(1.1203203587271555e+07, new int[] { 2, 12 });
                p.AddCoeff(-4.7398169023071964e+06, new int[] { 0, 16 });
                p.AddCoeff(-1.6250800807910388e+07, new int[] { 2, 14 });
                p.AddCoeff(2.2924604625538074e+06, new int[] { 0, 18 });
                p.AddCoeff(1.4219450706921589e+07, new int[] { 2, 16 });
                p.AddCoeff(-4.7055767389262363e+05, new int[] { 0, 20 });
                p.AddCoeff(-6.8773813876614223e+06, new int[] { 2, 18 });
                p.AddCoeff(1.4116730216778709e+06, new int[] { 2, 20 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{C4FB3213-BF00-4E17-9BDC-4A47DBCF5388}"));
                OrthonormalPolynomials[256] = p;
                p.AddCoeff(4.3668796235606764e+01, new int[] { 1, 1 });
                p.AddCoeff(-2.7511341628432261e+03, new int[] { 1, 3 });
                p.AddCoeff(-7.2781327059344606e+01, new int[] { 3, 1 });
                p.AddCoeff(5.0620868596315361e+04, new int[] { 1, 5 });
                p.AddCoeff(4.5852236047387102e+03, new int[] { 3, 3 });
                p.AddCoeff(-4.2184057163596134e+05, new int[] { 1, 7 });
                p.AddCoeff(-8.4368114327192268e+04, new int[] { 3, 5 });
                p.AddCoeff(1.8982825723618260e+06, new int[] { 1, 9 });
                p.AddCoeff(7.0306761939326890e+05, new int[] { 3, 7 });
                p.AddCoeff(-5.0045631453175413e+06, new int[] { 1, 11 });
                p.AddCoeff(-3.1638042872697100e+06, new int[] { 3, 9 });
                p.AddCoeff(7.9559721797355785e+06, new int[] { 1, 13 });
                p.AddCoeff(8.3409385755292355e+06, new int[] { 3, 11 });
                p.AddCoeff(-7.5013451980364026e+06, new int[] { 1, 15 });
                p.AddCoeff(-1.3259953632892631e+07, new int[] { 3, 13 });
                p.AddCoeff(3.8609864989893249e+06, new int[] { 1, 17 });
                p.AddCoeff(1.2502241996727338e+07, new int[] { 3, 15 });
                p.AddCoeff(-8.3541813136026328e+05, new int[] { 1, 19 });
                p.AddCoeff(-6.4349774983155415e+06, new int[] { 3, 17 });
                p.AddCoeff(1.3923635522671055e+06, new int[] { 3, 19 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{5DE3594E-9F9F-4CCB-B8CB-9B2B0ADBB5AF}"));
                OrthonormalPolynomials[257] = p;
                p.AddCoeff(-6.3459759426305174e-01, new int[] { 0, 0 });
                p.AddCoeff(1.0851618861898185e+02, new int[] { 0, 2 });
                p.AddCoeff(6.3459759426305174e+00, new int[] { 2, 0 });
                p.AddCoeff(-3.0384532813314917e+03, new int[] { 0, 4 });
                p.AddCoeff(-1.0851618861898185e+03, new int[] { 2, 2 });
                p.AddCoeff(-7.4036385997356037e+00, new int[] { 4, 0 });
                p.AddCoeff(3.2612731886291345e+04, new int[] { 0, 6 });
                p.AddCoeff(3.0384532813314917e+04, new int[] { 2, 4 });
                p.AddCoeff(1.2660222005547882e+03, new int[] { 4, 2 });
                p.AddCoeff(-1.7471106367656077e+05, new int[] { 0, 8 });
                p.AddCoeff(-3.2612731886291345e+05, new int[] { 2, 6 });
                p.AddCoeff(-3.5448621615534070e+04, new int[] { 4, 4 });
                p.AddCoeff(5.2413319102968232e+05, new int[] { 0, 10 });
                p.AddCoeff(1.7471106367656077e+06, new int[] { 2, 8 });
                p.AddCoeff(3.8048187200673235e+05, new int[] { 4, 6 });
                p.AddCoeff(-9.2120379029459318e+05, new int[] { 0, 12 });
                p.AddCoeff(-5.2413319102968232e+06, new int[] { 2, 10 });
                p.AddCoeff(-2.0382957428932090e+06, new int[] { 4, 8 });
                p.AddCoeff(9.4145002744392490e+05, new int[] { 0, 14 });
                p.AddCoeff(9.2120379029459318e+06, new int[] { 2, 12 });
                p.AddCoeff(6.1148872286796271e+06, new int[] { 4, 10 });
                p.AddCoeff(-5.1779751509415869e+05, new int[] { 0, 16 });
                p.AddCoeff(-9.4145002744392490e+06, new int[] { 2, 14 });
                p.AddCoeff(-1.0747377553436920e+07, new int[] { 4, 12 });
                p.AddCoeff(1.1845041194964415e+05, new int[] { 0, 18 });
                p.AddCoeff(5.1779751509415869e+06, new int[] { 2, 16 });
                p.AddCoeff(1.0983583653512457e+07, new int[] { 4, 14 });
                p.AddCoeff(-1.1845041194964415e+06, new int[] { 2, 18 });
                p.AddCoeff(-6.0409710094318514e+06, new int[] { 4, 16 });
                p.AddCoeff(1.3819214727458484e+06, new int[] { 4, 18 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{00F9E7B4-B025-43AA-B992-FB695FD55972}"));
                OrthonormalPolynomials[258] = p;
                p.AddCoeff(6.1411425548509120e+01, new int[] { 1, 1 });
                p.AddCoeff(-3.1115122277911287e+03, new int[] { 1, 3 });
                p.AddCoeff(-2.8658665255970923e+02, new int[] { 3, 1 });
                p.AddCoeff(4.5739229748529592e+04, new int[] { 1, 5 });
                p.AddCoeff(1.4520390396358601e+04, new int[] { 3, 3 });
                p.AddCoeff(2.5792798730373830e+02, new int[] { 5, 1 });
                p.AddCoeff(-3.0057208120462304e+05, new int[] { 1, 7 });
                p.AddCoeff(-2.1344973882647143e+05, new int[] { 3, 5 });
                p.AddCoeff(-1.3068351356722741e+04, new int[] { 5, 3 });
                p.AddCoeff(1.0436530597382744e+06, new int[] { 1, 9 });
                p.AddCoeff(1.4026697122882408e+06, new int[] { 3, 7 });
                p.AddCoeff(1.9210476494382429e+05, new int[] { 5, 5 });
                p.AddCoeff(-2.0493550991224298e+06, new int[] { 1, 11 });
                p.AddCoeff(-4.8703809454452807e+06, new int[] { 3, 9 });
                p.AddCoeff(-1.2624027410594167e+06, new int[] { 5, 7 });
                p.AddCoeff(2.2858191490211717e+06, new int[] { 1, 13 });
                p.AddCoeff(9.5636571292380057e+06, new int[] { 3, 11 });
                p.AddCoeff(4.3833428509007526e+06, new int[] { 5, 9 });
                p.AddCoeff(-1.3497217832315490e+06, new int[] { 1, 15 });
                p.AddCoeff(-1.0667156028765468e+07, new int[] { 3, 13 });
                p.AddCoeff(-8.6072914163142051e+06, new int[] { 5, 11 });
                p.AddCoeff(3.2750602093118468e+05, new int[] { 1, 17 });
                p.AddCoeff(6.2987016550805620e+06, new int[] { 3, 15 });
                p.AddCoeff(9.6004404258889211e+06, new int[] { 5, 13 });
                p.AddCoeff(-1.5283614310121952e+06, new int[] { 3, 17 });
                p.AddCoeff(-5.6688314895725058e+06, new int[] { 5, 15 });
                p.AddCoeff(1.3755252879109757e+06, new int[] { 5, 17 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{548483C4-BC87-4DEF-BE48-7E10F012E3A3}"));
                OrthonormalPolynomials[259] = p;
                p.AddCoeff(-6.3554643709818528e-01, new int[] { 0, 0 });
                p.AddCoeff(8.6434315445353198e+01, new int[] { 0, 2 });
                p.AddCoeff(1.3346475179061891e+01, new int[] { 2, 0 });
                p.AddCoeff(-1.9159606590386626e+03, new int[] { 0, 4 });
                p.AddCoeff(-1.8151206243524172e+03, new int[] { 2, 2 });
                p.AddCoeff(-4.0039425537185673e+01, new int[] { 4, 0 });
                p.AddCoeff(1.6094069535924765e+04, new int[] { 0, 6 });
                p.AddCoeff(4.0235173839811914e+04, new int[] { 2, 4 });
                p.AddCoeff(5.4453618730572515e+03, new int[] { 4, 2 });
                p.AddCoeff(2.9362245393936160e+01, new int[] { 6, 0 });
                p.AddCoeff(-6.6100642736833858e+04, new int[] { 0, 8 });
                p.AddCoeff(-3.3797546025442008e+05, new int[] { 2, 6 });
                p.AddCoeff(-1.2070552151943574e+05, new int[] { 4, 4 });
                p.AddCoeff(-3.9932653735753178e+03, new int[] { 6, 2 });
                p.AddCoeff(1.4689031719296413e+05, new int[] { 0, 10 });
                p.AddCoeff(1.3881134974735110e+06, new int[] { 2, 8 });
                p.AddCoeff(1.0139263807632602e+06, new int[] { 4, 6 });
                p.AddCoeff(8.8517382447586210e+04, new int[] { 6, 4 });
                p.AddCoeff(-1.8027448019136507e+05, new int[] { 0, 12 });
                p.AddCoeff(-3.0846966610522467e+06, new int[] { 2, 10 });
                p.AddCoeff(-4.1643404924205331e+06, new int[] { 4, 8 });
                p.AddCoeff(-7.4354601255972417e+05, new int[] { 6, 6 });
                p.AddCoeff(1.1490021814394697e+05, new int[] { 0, 14 });
                p.AddCoeff(3.7857640840186664e+06, new int[] { 2, 12 });
                p.AddCoeff(9.2540899831567402e+06, new int[] { 4, 10 });
                p.AddCoeff(3.0538496944417243e+06, new int[] { 6, 8 });
                p.AddCoeff(-2.9682556353852966e+04, new int[] { 0, 16 });
                p.AddCoeff(-2.4129045810228863e+06, new int[] { 2, 14 });
                p.AddCoeff(-1.1357292252055999e+07, new int[] { 4, 12 });
                p.AddCoeff(-6.7863326543149428e+06, new int[] { 6, 10 });
                p.AddCoeff(6.2333368343091229e+05, new int[] { 2, 16 });
                p.AddCoeff(7.2387137430686589e+06, new int[] { 4, 14 });
                p.AddCoeff(8.3286809848410661e+06, new int[] { 6, 12 });
                p.AddCoeff(-1.8700010502927369e+06, new int[] { 4, 16 });
                p.AddCoeff(-5.3083900782503498e+06, new int[] { 6, 14 });
                p.AddCoeff(1.3713341035480070e+06, new int[] { 6, 16 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{F32B2A4C-F2EA-4715-844A-5B3917EFE577}"));
                OrthonormalPolynomials[260] = p;
                p.AddCoeff(7.4107667008783742e+01, new int[] { 1, 1 });
                p.AddCoeff(-2.9396041246817551e+03, new int[] { 1, 3 });
                p.AddCoeff(-6.6696900307905368e+02, new int[] { 3, 1 });
                p.AddCoeff(3.3511487021372008e+04, new int[] { 1, 5 });
                p.AddCoeff(2.6456437122135796e+04, new int[] { 3, 3 });
                p.AddCoeff(1.4673318067739181e+03, new int[] { 5, 1 });
                p.AddCoeff(-1.6755743510686004e+05, new int[] { 1, 7 });
                p.AddCoeff(-3.0160338319234807e+05, new int[] { 3, 5 });
                p.AddCoeff(-5.8204161668698751e+04, new int[] { 5, 3 });
                p.AddCoeff(-9.0834826133623501e+02, new int[] { 7, 1 });
                p.AddCoeff(4.2820233416197566e+05, new int[] { 1, 9 });
                p.AddCoeff(1.5080169159617404e+06, new int[] { 3, 7 });
                p.AddCoeff(6.6352744302316576e+05, new int[] { 5, 5 });
                p.AddCoeff(3.6031147699670655e+04, new int[] { 7, 3 });
                p.AddCoeff(-5.8391227385723954e+05, new int[] { 1, 11 });
                p.AddCoeff(-3.8538210074577809e+06, new int[] { 3, 9 });
                p.AddCoeff(-3.3176372151158288e+06, new int[] { 5, 7 });
                p.AddCoeff(-4.1075508377624547e+05, new int[] { 7, 5 });
                p.AddCoeff(4.0424695882424276e+05, new int[] { 1, 13 });
                p.AddCoeff(5.2552104647151558e+06, new int[] { 3, 11 });
                p.AddCoeff(8.4784062164071181e+06, new int[] { 5, 9 });
                p.AddCoeff(2.0537754188812274e+06, new int[] { 7, 7 });
                p.AddCoeff(-1.1164916005621943e+05, new int[] { 1, 15 });
                p.AddCoeff(-3.6382226294181848e+06, new int[] { 3, 13 });
                p.AddCoeff(-1.1561463022373343e+07, new int[] { 5, 11 });
                p.AddCoeff(-5.2485371815853588e+06, new int[] { 7, 9 });
                p.AddCoeff(1.0048424405059748e+06, new int[] { 3, 15 });
                p.AddCoeff(8.0040897847200066e+06, new int[] { 5, 13 });
                p.AddCoeff(7.1570961567073074e+06, new int[] { 7, 11 });
                p.AddCoeff(-2.2106533691131447e+06, new int[] { 5, 15 });
                p.AddCoeff(-4.9549127238742898e+06, new int[] { 7, 13 });
                p.AddCoeff(1.3684997046890896e+06, new int[] { 7, 15 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{EE7BD099-6B69-4652-9002-54F27E3FFB4C}"));
                OrthonormalPolynomials[261] = p;
                p.AddCoeff(-6.3588543647248341e-01, new int[] { 0, 0 });
                p.AddCoeff(6.6767970829610758e+01, new int[] { 0, 2 });
                p.AddCoeff(2.2891875713009403e+01, new int[] { 2, 0 });
                p.AddCoeff(-1.1350555041033829e+03, new int[] { 0, 4 });
                p.AddCoeff(-2.4036469498659873e+03, new int[] { 2, 2 });
                p.AddCoeff(-1.2590531642155171e+02, new int[] { 4, 0 });
                p.AddCoeff(7.1886848593214249e+03, new int[] { 0, 6 });
                p.AddCoeff(4.0861998147721784e+04, new int[] { 2, 4 });
                p.AddCoeff(1.3220058224262930e+04, new int[] { 4, 2 });
                p.AddCoeff(2.1823588179735631e+02, new int[] { 6, 0 });
                p.AddCoeff(-2.1566054577964275e+04, new int[] { 0, 8 });
                p.AddCoeff(-2.5879265493557130e+05, new int[] { 2, 6 });
                p.AddCoeff(-2.2474098981246981e+05, new int[] { 4, 4 });
                p.AddCoeff(-2.2914767588722412e+04, new int[] { 6, 2 });
                p.AddCoeff(-1.1691207953429802e+02, new int[] { 8, 0 });
                p.AddCoeff(3.3067950352878555e+04, new int[] { 0, 10 });
                p.AddCoeff(7.7637796480671389e+05, new int[] { 2, 8 });
                p.AddCoeff(1.4233596021456421e+06, new int[] { 4, 6 });
                p.AddCoeff(3.8955104900828101e+05, new int[] { 6, 4 });
                p.AddCoeff(1.2275768351101292e+04, new int[] { 8, 2 });
                p.AddCoeff(-2.5051477540059511e+04, new int[] { 0, 12 });
                p.AddCoeff(-1.1904462127036280e+06, new int[] { 2, 10 });
                p.AddCoeff(-4.2700788064369264e+06, new int[] { 4, 8 });
                p.AddCoeff(-2.4671566437191130e+06, new int[] { 6, 6 });
                p.AddCoeff(-2.0868806196872197e+05, new int[] { 8, 4 });
                p.AddCoeff(7.4328559734242506e+03, new int[] { 0, 14 });
                p.AddCoeff(9.0185319144214240e+05, new int[] { 2, 12 });
                p.AddCoeff(6.5474541698699538e+06, new int[] { 4, 10 });
                p.AddCoeff(7.4014699311573391e+06, new int[] { 6, 8 });
                p.AddCoeff(1.3216910591352391e+06, new int[] { 8, 6 });
                p.AddCoeff(-2.6758281504327302e+05, new int[] { 2, 14 });
                p.AddCoeff(-4.9601925529317832e+06, new int[] { 4, 12 });
                p.AddCoeff(-1.1348920561107920e+07, new int[] { 6, 10 });
                p.AddCoeff(-3.9650731774057174e+06, new int[] { 8, 8 });
                p.AddCoeff(1.4717054827380016e+06, new int[] { 4, 14 });
                p.AddCoeff(8.5976670917484242e+06, new int[] { 6, 12 });
                p.AddCoeff(6.0797788720221000e+06, new int[] { 8, 10 });
                p.AddCoeff(-2.5509561700792028e+06, new int[] { 6, 14 });
                p.AddCoeff(-4.6058930848652273e+06, new int[] { 8, 12 });
                p.AddCoeff(1.3665836625424301e+06, new int[] { 8, 14 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{E65E3076-BF25-44CE-95C5-AC4D28690811}"));
                OrthonormalPolynomials[262] = p;
                p.AddCoeff(8.1730592363802215e+01, new int[] { 1, 1 });
                p.AddCoeff(-2.4519177709140665e+03, new int[] { 1, 3 });
                p.AddCoeff(-1.1987153546690992e+03, new int[] { 3, 1 });
                p.AddCoeff(2.0841301052769565e+04, new int[] { 1, 5 });
                p.AddCoeff(3.5961460640072975e+04, new int[] { 3, 3 });
                p.AddCoeff(4.6749898832094867e+03, new int[] { 5, 1 });
                p.AddCoeff(-7.5425660952880330e+04, new int[] { 1, 7 });
                p.AddCoeff(-3.0567241544062029e+05, new int[] { 3, 5 });
                p.AddCoeff(-1.4024969649628460e+05, new int[] { 5, 3 });
                p.AddCoeff(-6.6785569760135525e+03, new int[] { 7, 1 });
                p.AddCoeff(1.3199490666754058e+05, new int[] { 1, 9 });
                p.AddCoeff(1.1062430273089115e+06, new int[] { 3, 7 });
                p.AddCoeff(1.1921224202184191e+06, new int[] { 5, 5 });
                p.AddCoeff(2.0035670928040657e+05, new int[] { 7, 3 });
                p.AddCoeff(3.1537630164508442e+03, new int[] { 9, 1 });
                p.AddCoeff(-1.1039574012194303e+05, new int[] { 1, 11 });
                p.AddCoeff(-1.9359252977905951e+06, new int[] { 3, 9 });
                p.AddCoeff(-4.3143478065047549e+06, new int[] { 5, 7 });
                p.AddCoeff(-1.7030320288834559e+06, new int[] { 7, 5 });
                p.AddCoeff(-9.4612890493525327e+04, new int[] { 9, 3 });
                p.AddCoeff(3.5383250039084304e+04, new int[] { 1, 13 });
                p.AddCoeff(1.6191375217884978e+06, new int[] { 3, 11 });
                p.AddCoeff(7.5501086613833211e+06, new int[] { 5, 9 });
                p.AddCoeff(6.1633540092925070e+06, new int[] { 7, 7 });
                p.AddCoeff(8.0420956919496528e+05, new int[] { 9, 5 });
                p.AddCoeff(-5.1895433390656979e+05, new int[] { 3, 13 });
                p.AddCoeff(-6.3146363349751412e+06, new int[] { 5, 11 });
                p.AddCoeff(-1.0785869516261887e+07, new int[] { 7, 9 });
                p.AddCoeff(-2.9104727266103505e+06, new int[] { 9, 7 });
                p.AddCoeff(2.0239219022356222e+06, new int[] { 5, 13 });
                p.AddCoeff(9.0209090499644875e+06, new int[] { 7, 11 });
                p.AddCoeff(5.0933272715681134e+06, new int[] { 9, 9 });
                p.AddCoeff(-2.8913170031937460e+06, new int[] { 7, 13 });
                p.AddCoeff(-4.2598737180387858e+06, new int[] { 9, 11 });
                p.AddCoeff(1.3653441403970467e+06, new int[] { 9, 13 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{C8AA8484-0DF4-48BE-A668-F2EB483644B2}"));
                OrthonormalPolynomials[263] = p;
                p.AddCoeff(-6.3600753868763294e-01, new int[] { 0, 0 });
                p.AddCoeff(4.9608588017635369e+01, new int[] { 0, 2 });
                p.AddCoeff(3.4980414627819812e+01, new int[] { 2, 0 });
                p.AddCoeff(-6.2010735022044212e+02, new int[] { 0, 4 });
                p.AddCoeff(-2.7284723409699453e+03, new int[] { 2, 2 });
                p.AddCoeff(-3.0316359344110503e+02, new int[] { 4, 0 });
                p.AddCoeff(2.8111533209993376e+03, new int[] { 0, 6 });
                p.AddCoeff(3.4105904262124316e+04, new int[] { 2, 4 });
                p.AddCoeff(2.3646760288406193e+04, new int[] { 4, 2 });
                p.AddCoeff(9.0949078032331510e+02, new int[] { 6, 0 });
                p.AddCoeff(-5.7227049748915087e+03, new int[] { 0, 8 });
                p.AddCoeff(-1.5461343265496357e+05, new int[] { 2, 6 });
                p.AddCoeff(-2.9558450360507741e+05, new int[] { 4, 4 });
                p.AddCoeff(-7.0940280865218578e+04, new int[] { 6, 2 });
                p.AddCoeff(-1.1043816618211683e+03, new int[] { 8, 0 });
                p.AddCoeff(5.3411913098987414e+03, new int[] { 0, 10 });
                p.AddCoeff(3.1474877361903298e+05, new int[] { 2, 8 });
                p.AddCoeff(1.3399830830096843e+06, new int[] { 4, 6 });
                p.AddCoeff(8.8675351081523222e+05, new int[] { 6, 4 });
                p.AddCoeff(8.6141769622051130e+04, new int[] { 8, 2 });
                p.AddCoeff(4.6629447943560441e+02, new int[] { 10, 0 });
                p.AddCoeff(-1.8613242443586523e+03, new int[] { 0, 12 });
                p.AddCoeff(-2.9376552204443078e+05, new int[] { 2, 10 });
                p.AddCoeff(-2.7278227046982858e+06, new int[] { 4, 8 });
                p.AddCoeff(-4.0199492490290528e+06, new int[] { 6, 6 });
                p.AddCoeff(-1.0767721202756391e+06, new int[] { 8, 4 });
                p.AddCoeff(-3.6370969395977144e+04, new int[] { 10, 2 });
                p.AddCoeff(1.0237283343972588e+05, new int[] { 2, 12 });
                p.AddCoeff(2.5459678577184001e+06, new int[] { 4, 10 });
                p.AddCoeff(8.1834681140948574e+06, new int[] { 6, 8 });
                p.AddCoeff(4.8813669452495641e+06, new int[] { 8, 6 });
                p.AddCoeff(4.5463711744971430e+05, new int[] { 10, 4 });
                p.AddCoeff(-8.8723122314429094e+05, new int[] { 4, 12 });
                p.AddCoeff(-7.6379035731552002e+06, new int[] { 6, 10 });
                p.AddCoeff(-9.9370684242580411e+06, new int[] { 8, 8 });
                p.AddCoeff(-2.0610215991053715e+06, new int[] { 10, 6 });
                p.AddCoeff(2.6616936694328728e+06, new int[] { 6, 12 });
                p.AddCoeff(9.2745971959741717e+06, new int[] { 8, 10 });
                p.AddCoeff(4.1956511124645062e+06, new int[] { 10, 8 });
                p.AddCoeff(-3.2320565985970598e+06, new int[] { 8, 12 });
                p.AddCoeff(-3.9159410383002058e+06, new int[] { 10, 10 });
                p.AddCoeff(1.3646461194076475e+06, new int[] { 10, 12 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{3FC6B4DA-25DD-4E63-B410-2BF2DCE44757}"));
                OrthonormalPolynomials[264] = p;
                p.AddCoeff(8.4272209167480469e+01, new int[] { 1, 1 });
                p.AddCoeff(-1.8258978652954102e+03, new int[] { 1, 3 });
                p.AddCoeff(-1.8258978652954102e+03, new int[] { 3, 1 });
                p.AddCoeff(1.0955387191772461e+04, new int[] { 1, 5 });
                p.AddCoeff(3.9561120414733887e+04, new int[] { 3, 3 });
                p.AddCoeff(1.0955387191772461e+04, new int[] { 5, 1 });
                p.AddCoeff(-2.6605940322875977e+04, new int[] { 1, 7 });
                p.AddCoeff(-2.3736672248840332e+05, new int[] { 3, 5 });
                p.AddCoeff(-2.3736672248840332e+05, new int[] { 5, 3 });
                p.AddCoeff(-2.6605940322875977e+04, new int[] { 7, 1 });
                p.AddCoeff(2.8084048118591309e+04, new int[] { 1, 9 });
                p.AddCoeff(5.7646204032897949e+05, new int[] { 3, 7 });
                p.AddCoeff(1.4242003349304199e+06, new int[] { 5, 5 });
                p.AddCoeff(5.7646204032897949e+05, new int[] { 7, 3 });
                p.AddCoeff(2.8084048118591309e+04, new int[] { 9, 1 });
                p.AddCoeff(-1.0723000190734863e+04, new int[] { 1, 11 });
                p.AddCoeff(-6.0848770923614502e+05, new int[] { 3, 9 });
                p.AddCoeff(-3.4587722419738770e+06, new int[] { 5, 7 });
                p.AddCoeff(-3.4587722419738770e+06, new int[] { 7, 5 });
                p.AddCoeff(-6.0848770923614502e+05, new int[] { 9, 3 });
                p.AddCoeff(-1.0723000190734863e+04, new int[] { 11, 1 });
                p.AddCoeff(2.3233167079925537e+05, new int[] { 3, 11 });
                p.AddCoeff(3.6509262554168701e+06, new int[] { 5, 9 });
                p.AddCoeff(8.3998754447937012e+06, new int[] { 7, 7 });
                p.AddCoeff(3.6509262554168701e+06, new int[] { 9, 5 });
                p.AddCoeff(2.3233167079925537e+05, new int[] { 11, 3 });
                p.AddCoeff(-1.3939900247955322e+06, new int[] { 5, 11 });
                p.AddCoeff(-8.8665351917266846e+06, new int[] { 7, 9 });
                p.AddCoeff(-8.8665351917266846e+06, new int[] { 9, 7 });
                p.AddCoeff(-1.3939900247955322e+06, new int[] { 11, 5 });
                p.AddCoeff(3.3854043459320068e+06, new int[] { 7, 11 });
                p.AddCoeff(9.3591204801559448e+06, new int[] { 9, 9 });
                p.AddCoeff(3.3854043459320068e+06, new int[] { 11, 7 });
                p.AddCoeff(-3.5734823651504517e+06, new int[] { 9, 11 });
                p.AddCoeff(-3.5734823651504517e+06, new int[] { 11, 9 });
                p.AddCoeff(1.3644205394210815e+06, new int[] { 11, 11 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{AEB306C6-2281-4D72-93A8-B137A5BC6134}"));
                OrthonormalPolynomials[265] = p;
                p.AddCoeff(-6.3600753868763294e-01, new int[] { 0, 0 });
                p.AddCoeff(3.4980414627819812e+01, new int[] { 0, 2 });
                p.AddCoeff(4.9608588017635369e+01, new int[] { 2, 0 });
                p.AddCoeff(-3.0316359344110503e+02, new int[] { 0, 4 });
                p.AddCoeff(-2.7284723409699453e+03, new int[] { 2, 2 });
                p.AddCoeff(-6.2010735022044212e+02, new int[] { 4, 0 });
                p.AddCoeff(9.0949078032331510e+02, new int[] { 0, 6 });
                p.AddCoeff(2.3646760288406193e+04, new int[] { 2, 4 });
                p.AddCoeff(3.4105904262124316e+04, new int[] { 4, 2 });
                p.AddCoeff(2.8111533209993376e+03, new int[] { 6, 0 });
                p.AddCoeff(-1.1043816618211683e+03, new int[] { 0, 8 });
                p.AddCoeff(-7.0940280865218578e+04, new int[] { 2, 6 });
                p.AddCoeff(-2.9558450360507741e+05, new int[] { 4, 4 });
                p.AddCoeff(-1.5461343265496357e+05, new int[] { 6, 2 });
                p.AddCoeff(-5.7227049748915087e+03, new int[] { 8, 0 });
                p.AddCoeff(4.6629447943560441e+02, new int[] { 0, 10 });
                p.AddCoeff(8.6141769622051130e+04, new int[] { 2, 8 });
                p.AddCoeff(8.8675351081523222e+05, new int[] { 4, 6 });
                p.AddCoeff(1.3399830830096843e+06, new int[] { 6, 4 });
                p.AddCoeff(3.1474877361903298e+05, new int[] { 8, 2 });
                p.AddCoeff(5.3411913098987414e+03, new int[] { 10, 0 });
                p.AddCoeff(-3.6370969395977144e+04, new int[] { 2, 10 });
                p.AddCoeff(-1.0767721202756391e+06, new int[] { 4, 8 });
                p.AddCoeff(-4.0199492490290528e+06, new int[] { 6, 6 });
                p.AddCoeff(-2.7278227046982858e+06, new int[] { 8, 4 });
                p.AddCoeff(-2.9376552204443078e+05, new int[] { 10, 2 });
                p.AddCoeff(-1.8613242443586523e+03, new int[] { 12, 0 });
                p.AddCoeff(4.5463711744971430e+05, new int[] { 4, 10 });
                p.AddCoeff(4.8813669452495641e+06, new int[] { 6, 8 });
                p.AddCoeff(8.1834681140948574e+06, new int[] { 8, 6 });
                p.AddCoeff(2.5459678577184001e+06, new int[] { 10, 4 });
                p.AddCoeff(1.0237283343972588e+05, new int[] { 12, 2 });
                p.AddCoeff(-2.0610215991053715e+06, new int[] { 6, 10 });
                p.AddCoeff(-9.9370684242580411e+06, new int[] { 8, 8 });
                p.AddCoeff(-7.6379035731552002e+06, new int[] { 10, 6 });
                p.AddCoeff(-8.8723122314429094e+05, new int[] { 12, 4 });
                p.AddCoeff(4.1956511124645062e+06, new int[] { 8, 10 });
                p.AddCoeff(9.2745971959741717e+06, new int[] { 10, 8 });
                p.AddCoeff(2.6616936694328728e+06, new int[] { 12, 6 });
                p.AddCoeff(-3.9159410383002058e+06, new int[] { 10, 10 });
                p.AddCoeff(-3.2320565985970598e+06, new int[] { 12, 8 });
                p.AddCoeff(1.3646461194076475e+06, new int[] { 12, 10 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{71AE3546-F9CC-447D-8F38-0766D4B8887A}"));
                OrthonormalPolynomials[266] = p;
                p.AddCoeff(8.1730592363802215e+01, new int[] { 1, 1 });
                p.AddCoeff(-1.1987153546690992e+03, new int[] { 1, 3 });
                p.AddCoeff(-2.4519177709140665e+03, new int[] { 3, 1 });
                p.AddCoeff(4.6749898832094867e+03, new int[] { 1, 5 });
                p.AddCoeff(3.5961460640072975e+04, new int[] { 3, 3 });
                p.AddCoeff(2.0841301052769565e+04, new int[] { 5, 1 });
                p.AddCoeff(-6.6785569760135525e+03, new int[] { 1, 7 });
                p.AddCoeff(-1.4024969649628460e+05, new int[] { 3, 5 });
                p.AddCoeff(-3.0567241544062029e+05, new int[] { 5, 3 });
                p.AddCoeff(-7.5425660952880330e+04, new int[] { 7, 1 });
                p.AddCoeff(3.1537630164508442e+03, new int[] { 1, 9 });
                p.AddCoeff(2.0035670928040657e+05, new int[] { 3, 7 });
                p.AddCoeff(1.1921224202184191e+06, new int[] { 5, 5 });
                p.AddCoeff(1.1062430273089115e+06, new int[] { 7, 3 });
                p.AddCoeff(1.3199490666754058e+05, new int[] { 9, 1 });
                p.AddCoeff(-9.4612890493525327e+04, new int[] { 3, 9 });
                p.AddCoeff(-1.7030320288834559e+06, new int[] { 5, 7 });
                p.AddCoeff(-4.3143478065047549e+06, new int[] { 7, 5 });
                p.AddCoeff(-1.9359252977905951e+06, new int[] { 9, 3 });
                p.AddCoeff(-1.1039574012194303e+05, new int[] { 11, 1 });
                p.AddCoeff(8.0420956919496528e+05, new int[] { 5, 9 });
                p.AddCoeff(6.1633540092925070e+06, new int[] { 7, 7 });
                p.AddCoeff(7.5501086613833211e+06, new int[] { 9, 5 });
                p.AddCoeff(1.6191375217884978e+06, new int[] { 11, 3 });
                p.AddCoeff(3.5383250039084304e+04, new int[] { 13, 1 });
                p.AddCoeff(-2.9104727266103505e+06, new int[] { 7, 9 });
                p.AddCoeff(-1.0785869516261887e+07, new int[] { 9, 7 });
                p.AddCoeff(-6.3146363349751412e+06, new int[] { 11, 5 });
                p.AddCoeff(-5.1895433390656979e+05, new int[] { 13, 3 });
                p.AddCoeff(5.0933272715681134e+06, new int[] { 9, 9 });
                p.AddCoeff(9.0209090499644875e+06, new int[] { 11, 7 });
                p.AddCoeff(2.0239219022356222e+06, new int[] { 13, 5 });
                p.AddCoeff(-4.2598737180387858e+06, new int[] { 11, 9 });
                p.AddCoeff(-2.8913170031937460e+06, new int[] { 13, 7 });
                p.AddCoeff(1.3653441403970467e+06, new int[] { 13, 9 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{BE818372-AF58-4858-8A7D-AB18962E2CA9}"));
                OrthonormalPolynomials[267] = p;
                p.AddCoeff(-6.3588543647248341e-01, new int[] { 0, 0 });
                p.AddCoeff(2.2891875713009403e+01, new int[] { 0, 2 });
                p.AddCoeff(6.6767970829610758e+01, new int[] { 2, 0 });
                p.AddCoeff(-1.2590531642155171e+02, new int[] { 0, 4 });
                p.AddCoeff(-2.4036469498659873e+03, new int[] { 2, 2 });
                p.AddCoeff(-1.1350555041033829e+03, new int[] { 4, 0 });
                p.AddCoeff(2.1823588179735631e+02, new int[] { 0, 6 });
                p.AddCoeff(1.3220058224262930e+04, new int[] { 2, 4 });
                p.AddCoeff(4.0861998147721784e+04, new int[] { 4, 2 });
                p.AddCoeff(7.1886848593214249e+03, new int[] { 6, 0 });
                p.AddCoeff(-1.1691207953429802e+02, new int[] { 0, 8 });
                p.AddCoeff(-2.2914767588722412e+04, new int[] { 2, 6 });
                p.AddCoeff(-2.2474098981246981e+05, new int[] { 4, 4 });
                p.AddCoeff(-2.5879265493557130e+05, new int[] { 6, 2 });
                p.AddCoeff(-2.1566054577964275e+04, new int[] { 8, 0 });
                p.AddCoeff(1.2275768351101292e+04, new int[] { 2, 8 });
                p.AddCoeff(3.8955104900828101e+05, new int[] { 4, 6 });
                p.AddCoeff(1.4233596021456421e+06, new int[] { 6, 4 });
                p.AddCoeff(7.7637796480671389e+05, new int[] { 8, 2 });
                p.AddCoeff(3.3067950352878555e+04, new int[] { 10, 0 });
                p.AddCoeff(-2.0868806196872197e+05, new int[] { 4, 8 });
                p.AddCoeff(-2.4671566437191130e+06, new int[] { 6, 6 });
                p.AddCoeff(-4.2700788064369264e+06, new int[] { 8, 4 });
                p.AddCoeff(-1.1904462127036280e+06, new int[] { 10, 2 });
                p.AddCoeff(-2.5051477540059511e+04, new int[] { 12, 0 });
                p.AddCoeff(1.3216910591352391e+06, new int[] { 6, 8 });
                p.AddCoeff(7.4014699311573391e+06, new int[] { 8, 6 });
                p.AddCoeff(6.5474541698699538e+06, new int[] { 10, 4 });
                p.AddCoeff(9.0185319144214240e+05, new int[] { 12, 2 });
                p.AddCoeff(7.4328559734242506e+03, new int[] { 14, 0 });
                p.AddCoeff(-3.9650731774057174e+06, new int[] { 8, 8 });
                p.AddCoeff(-1.1348920561107920e+07, new int[] { 10, 6 });
                p.AddCoeff(-4.9601925529317832e+06, new int[] { 12, 4 });
                p.AddCoeff(-2.6758281504327302e+05, new int[] { 14, 2 });
                p.AddCoeff(6.0797788720221000e+06, new int[] { 10, 8 });
                p.AddCoeff(8.5976670917484242e+06, new int[] { 12, 6 });
                p.AddCoeff(1.4717054827380016e+06, new int[] { 14, 4 });
                p.AddCoeff(-4.6058930848652273e+06, new int[] { 12, 8 });
                p.AddCoeff(-2.5509561700792028e+06, new int[] { 14, 6 });
                p.AddCoeff(1.3665836625424301e+06, new int[] { 14, 8 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{5CA0D9D7-08F3-4D11-B5AC-0274DE1A14EE}"));
                OrthonormalPolynomials[268] = p;
                p.AddCoeff(7.4107667008783742e+01, new int[] { 1, 1 });
                p.AddCoeff(-6.6696900307905368e+02, new int[] { 1, 3 });
                p.AddCoeff(-2.9396041246817551e+03, new int[] { 3, 1 });
                p.AddCoeff(1.4673318067739181e+03, new int[] { 1, 5 });
                p.AddCoeff(2.6456437122135796e+04, new int[] { 3, 3 });
                p.AddCoeff(3.3511487021372008e+04, new int[] { 5, 1 });
                p.AddCoeff(-9.0834826133623501e+02, new int[] { 1, 7 });
                p.AddCoeff(-5.8204161668698751e+04, new int[] { 3, 5 });
                p.AddCoeff(-3.0160338319234807e+05, new int[] { 5, 3 });
                p.AddCoeff(-1.6755743510686004e+05, new int[] { 7, 1 });
                p.AddCoeff(3.6031147699670655e+04, new int[] { 3, 7 });
                p.AddCoeff(6.6352744302316576e+05, new int[] { 5, 5 });
                p.AddCoeff(1.5080169159617404e+06, new int[] { 7, 3 });
                p.AddCoeff(4.2820233416197566e+05, new int[] { 9, 1 });
                p.AddCoeff(-4.1075508377624547e+05, new int[] { 5, 7 });
                p.AddCoeff(-3.3176372151158288e+06, new int[] { 7, 5 });
                p.AddCoeff(-3.8538210074577809e+06, new int[] { 9, 3 });
                p.AddCoeff(-5.8391227385723954e+05, new int[] { 11, 1 });
                p.AddCoeff(2.0537754188812274e+06, new int[] { 7, 7 });
                p.AddCoeff(8.4784062164071181e+06, new int[] { 9, 5 });
                p.AddCoeff(5.2552104647151558e+06, new int[] { 11, 3 });
                p.AddCoeff(4.0424695882424276e+05, new int[] { 13, 1 });
                p.AddCoeff(-5.2485371815853588e+06, new int[] { 9, 7 });
                p.AddCoeff(-1.1561463022373343e+07, new int[] { 11, 5 });
                p.AddCoeff(-3.6382226294181848e+06, new int[] { 13, 3 });
                p.AddCoeff(-1.1164916005621943e+05, new int[] { 15, 1 });
                p.AddCoeff(7.1570961567073074e+06, new int[] { 11, 7 });
                p.AddCoeff(8.0040897847200066e+06, new int[] { 13, 5 });
                p.AddCoeff(1.0048424405059748e+06, new int[] { 15, 3 });
                p.AddCoeff(-4.9549127238742898e+06, new int[] { 13, 7 });
                p.AddCoeff(-2.2106533691131447e+06, new int[] { 15, 5 });
                p.AddCoeff(1.3684997046890896e+06, new int[] { 15, 7 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{2358692F-F048-4D8C-8C44-54E279830965}"));
                OrthonormalPolynomials[269] = p;
                p.AddCoeff(-6.3554643709818528e-01, new int[] { 0, 0 });
                p.AddCoeff(1.3346475179061891e+01, new int[] { 0, 2 });
                p.AddCoeff(8.6434315445353198e+01, new int[] { 2, 0 });
                p.AddCoeff(-4.0039425537185673e+01, new int[] { 0, 4 });
                p.AddCoeff(-1.8151206243524172e+03, new int[] { 2, 2 });
                p.AddCoeff(-1.9159606590386626e+03, new int[] { 4, 0 });
                p.AddCoeff(2.9362245393936160e+01, new int[] { 0, 6 });
                p.AddCoeff(5.4453618730572515e+03, new int[] { 2, 4 });
                p.AddCoeff(4.0235173839811914e+04, new int[] { 4, 2 });
                p.AddCoeff(1.6094069535924765e+04, new int[] { 6, 0 });
                p.AddCoeff(-3.9932653735753178e+03, new int[] { 2, 6 });
                p.AddCoeff(-1.2070552151943574e+05, new int[] { 4, 4 });
                p.AddCoeff(-3.3797546025442008e+05, new int[] { 6, 2 });
                p.AddCoeff(-6.6100642736833858e+04, new int[] { 8, 0 });
                p.AddCoeff(8.8517382447586210e+04, new int[] { 4, 6 });
                p.AddCoeff(1.0139263807632602e+06, new int[] { 6, 4 });
                p.AddCoeff(1.3881134974735110e+06, new int[] { 8, 2 });
                p.AddCoeff(1.4689031719296413e+05, new int[] { 10, 0 });
                p.AddCoeff(-7.4354601255972417e+05, new int[] { 6, 6 });
                p.AddCoeff(-4.1643404924205331e+06, new int[] { 8, 4 });
                p.AddCoeff(-3.0846966610522467e+06, new int[] { 10, 2 });
                p.AddCoeff(-1.8027448019136507e+05, new int[] { 12, 0 });
                p.AddCoeff(3.0538496944417243e+06, new int[] { 8, 6 });
                p.AddCoeff(9.2540899831567402e+06, new int[] { 10, 4 });
                p.AddCoeff(3.7857640840186664e+06, new int[] { 12, 2 });
                p.AddCoeff(1.1490021814394697e+05, new int[] { 14, 0 });
                p.AddCoeff(-6.7863326543149428e+06, new int[] { 10, 6 });
                p.AddCoeff(-1.1357292252055999e+07, new int[] { 12, 4 });
                p.AddCoeff(-2.4129045810228863e+06, new int[] { 14, 2 });
                p.AddCoeff(-2.9682556353852966e+04, new int[] { 16, 0 });
                p.AddCoeff(8.3286809848410661e+06, new int[] { 12, 6 });
                p.AddCoeff(7.2387137430686589e+06, new int[] { 14, 4 });
                p.AddCoeff(6.2333368343091229e+05, new int[] { 16, 2 });
                p.AddCoeff(-5.3083900782503498e+06, new int[] { 14, 6 });
                p.AddCoeff(-1.8700010502927369e+06, new int[] { 16, 4 });
                p.AddCoeff(1.3713341035480070e+06, new int[] { 16, 6 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{70FA5DA0-98E3-488D-B90E-35D48753673E}"));
                OrthonormalPolynomials[270] = p;
                p.AddCoeff(6.1411425548509120e+01, new int[] { 1, 1 });
                p.AddCoeff(-2.8658665255970923e+02, new int[] { 1, 3 });
                p.AddCoeff(-3.1115122277911287e+03, new int[] { 3, 1 });
                p.AddCoeff(2.5792798730373830e+02, new int[] { 1, 5 });
                p.AddCoeff(1.4520390396358601e+04, new int[] { 3, 3 });
                p.AddCoeff(4.5739229748529592e+04, new int[] { 5, 1 });
                p.AddCoeff(-1.3068351356722741e+04, new int[] { 3, 5 });
                p.AddCoeff(-2.1344973882647143e+05, new int[] { 5, 3 });
                p.AddCoeff(-3.0057208120462304e+05, new int[] { 7, 1 });
                p.AddCoeff(1.9210476494382429e+05, new int[] { 5, 5 });
                p.AddCoeff(1.4026697122882408e+06, new int[] { 7, 3 });
                p.AddCoeff(1.0436530597382744e+06, new int[] { 9, 1 });
                p.AddCoeff(-1.2624027410594167e+06, new int[] { 7, 5 });
                p.AddCoeff(-4.8703809454452807e+06, new int[] { 9, 3 });
                p.AddCoeff(-2.0493550991224298e+06, new int[] { 11, 1 });
                p.AddCoeff(4.3833428509007526e+06, new int[] { 9, 5 });
                p.AddCoeff(9.5636571292380057e+06, new int[] { 11, 3 });
                p.AddCoeff(2.2858191490211717e+06, new int[] { 13, 1 });
                p.AddCoeff(-8.6072914163142051e+06, new int[] { 11, 5 });
                p.AddCoeff(-1.0667156028765468e+07, new int[] { 13, 3 });
                p.AddCoeff(-1.3497217832315490e+06, new int[] { 15, 1 });
                p.AddCoeff(9.6004404258889211e+06, new int[] { 13, 5 });
                p.AddCoeff(6.2987016550805620e+06, new int[] { 15, 3 });
                p.AddCoeff(3.2750602093118468e+05, new int[] { 17, 1 });
                p.AddCoeff(-5.6688314895725058e+06, new int[] { 15, 5 });
                p.AddCoeff(-1.5283614310121952e+06, new int[] { 17, 3 });
                p.AddCoeff(1.3755252879109757e+06, new int[] { 17, 5 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{2384CE3F-92A1-4C89-AF6E-AEAA091B0B72}"));
                OrthonormalPolynomials[271] = p;
                p.AddCoeff(-6.3459759426305174e-01, new int[] { 0, 0 });
                p.AddCoeff(6.3459759426305174e+00, new int[] { 0, 2 });
                p.AddCoeff(1.0851618861898185e+02, new int[] { 2, 0 });
                p.AddCoeff(-7.4036385997356037e+00, new int[] { 0, 4 });
                p.AddCoeff(-1.0851618861898185e+03, new int[] { 2, 2 });
                p.AddCoeff(-3.0384532813314917e+03, new int[] { 4, 0 });
                p.AddCoeff(1.2660222005547882e+03, new int[] { 2, 4 });
                p.AddCoeff(3.0384532813314917e+04, new int[] { 4, 2 });
                p.AddCoeff(3.2612731886291345e+04, new int[] { 6, 0 });
                p.AddCoeff(-3.5448621615534070e+04, new int[] { 4, 4 });
                p.AddCoeff(-3.2612731886291345e+05, new int[] { 6, 2 });
                p.AddCoeff(-1.7471106367656077e+05, new int[] { 8, 0 });
                p.AddCoeff(3.8048187200673235e+05, new int[] { 6, 4 });
                p.AddCoeff(1.7471106367656077e+06, new int[] { 8, 2 });
                p.AddCoeff(5.2413319102968232e+05, new int[] { 10, 0 });
                p.AddCoeff(-2.0382957428932090e+06, new int[] { 8, 4 });
                p.AddCoeff(-5.2413319102968232e+06, new int[] { 10, 2 });
                p.AddCoeff(-9.2120379029459318e+05, new int[] { 12, 0 });
                p.AddCoeff(6.1148872286796271e+06, new int[] { 10, 4 });
                p.AddCoeff(9.2120379029459318e+06, new int[] { 12, 2 });
                p.AddCoeff(9.4145002744392490e+05, new int[] { 14, 0 });
                p.AddCoeff(-1.0747377553436920e+07, new int[] { 12, 4 });
                p.AddCoeff(-9.4145002744392490e+06, new int[] { 14, 2 });
                p.AddCoeff(-5.1779751509415869e+05, new int[] { 16, 0 });
                p.AddCoeff(1.0983583653512457e+07, new int[] { 14, 4 });
                p.AddCoeff(5.1779751509415869e+06, new int[] { 16, 2 });
                p.AddCoeff(1.1845041194964415e+05, new int[] { 18, 0 });
                p.AddCoeff(-6.0409710094318514e+06, new int[] { 16, 4 });
                p.AddCoeff(-1.1845041194964415e+06, new int[] { 18, 2 });
                p.AddCoeff(1.3819214727458484e+06, new int[] { 18, 4 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{79A38D1F-4140-48EC-A6CF-4A05FCD7FE2E}"));
                OrthonormalPolynomials[272] = p;
                p.AddCoeff(4.3668796235606764e+01, new int[] { 1, 1 });
                p.AddCoeff(-7.2781327059344606e+01, new int[] { 1, 3 });
                p.AddCoeff(-2.7511341628432261e+03, new int[] { 3, 1 });
                p.AddCoeff(4.5852236047387102e+03, new int[] { 3, 3 });
                p.AddCoeff(5.0620868596315361e+04, new int[] { 5, 1 });
                p.AddCoeff(-8.4368114327192268e+04, new int[] { 5, 3 });
                p.AddCoeff(-4.2184057163596134e+05, new int[] { 7, 1 });
                p.AddCoeff(7.0306761939326890e+05, new int[] { 7, 3 });
                p.AddCoeff(1.8982825723618260e+06, new int[] { 9, 1 });
                p.AddCoeff(-3.1638042872697100e+06, new int[] { 9, 3 });
                p.AddCoeff(-5.0045631453175413e+06, new int[] { 11, 1 });
                p.AddCoeff(8.3409385755292355e+06, new int[] { 11, 3 });
                p.AddCoeff(7.9559721797355785e+06, new int[] { 13, 1 });
                p.AddCoeff(-1.3259953632892631e+07, new int[] { 13, 3 });
                p.AddCoeff(-7.5013451980364026e+06, new int[] { 15, 1 });
                p.AddCoeff(1.2502241996727338e+07, new int[] { 15, 3 });
                p.AddCoeff(3.8609864989893249e+06, new int[] { 17, 1 });
                p.AddCoeff(-6.4349774983155415e+06, new int[] { 17, 3 });
                p.AddCoeff(-8.3541813136026328e+05, new int[] { 19, 1 });
                p.AddCoeff(1.3923635522671055e+06, new int[] { 19, 3 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{079D70BC-6329-4911-8268-032D8FCC05A5}"));
                OrthonormalPolynomials[273] = p;
                p.AddCoeff(-6.3068946561019085e-01, new int[] { 0, 0 });
                p.AddCoeff(1.8920683968305725e+00, new int[] { 0, 2 });
                p.AddCoeff(1.3244478777814008e+02, new int[] { 2, 0 });
                p.AddCoeff(-3.9733436333442023e+02, new int[] { 2, 2 });
                p.AddCoeff(-4.5693451783458327e+03, new int[] { 4, 0 });
                p.AddCoeff(1.3708035535037498e+04, new int[] { 4, 2 });
                p.AddCoeff(6.0924602377944436e+04, new int[] { 6, 0 });
                p.AddCoeff(-1.8277380713383331e+05, new int[] { 6, 2 });
                p.AddCoeff(-4.1124106605112494e+05, new int[] { 8, 0 });
                p.AddCoeff(1.2337231981533748e+06, new int[] { 8, 2 });
                p.AddCoeff(1.5901321220643498e+06, new int[] { 10, 0 });
                p.AddCoeff(-4.7703963661930493e+06, new int[] { 10, 2 });
                p.AddCoeff(-3.7344011957571851e+06, new int[] { 12, 0 });
                p.AddCoeff(1.1203203587271555e+07, new int[] { 12, 2 });
                p.AddCoeff(5.4169336026367959e+06, new int[] { 14, 0 });
                p.AddCoeff(-1.6250800807910388e+07, new int[] { 14, 2 });
                p.AddCoeff(-4.7398169023071964e+06, new int[] { 16, 0 });
                p.AddCoeff(1.4219450706921589e+07, new int[] { 16, 2 });
                p.AddCoeff(2.2924604625538074e+06, new int[] { 18, 0 });
                p.AddCoeff(-6.8773813876614223e+06, new int[] { 18, 2 });
                p.AddCoeff(-4.7055767389262363e+05, new int[] { 20, 0 });
                p.AddCoeff(1.4116730216778709e+06, new int[] { 20, 2 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{FABFB73C-41FA-491A-A091-22DB8D2AB66B}"));
                OrthonormalPolynomials[274] = p;
                p.AddCoeff(2.1012745091505225e+01, new int[] { 1, 1 });
                p.AddCoeff(-1.6109771236820672e+03, new int[] { 3, 1 });
                p.AddCoeff(3.6246985282846512e+04, new int[] { 5, 1 });
                p.AddCoeff(-3.7282613433784984e+05, new int[] { 7, 1 });
                p.AddCoeff(2.1023251464050977e+06, new int[] { 9, 1 });
                p.AddCoeff(-7.1096814042063305e+06, new int[] { 11, 1 });
                p.AddCoeff(1.5039710662744161e+07, new int[] { 13, 1 });
                p.AddCoeff(-2.0052947550325547e+07, new int[] { 15, 1 });
                p.AddCoeff(1.6366743956515704e+07, new int[] { 17, 1 });
                p.AddCoeff(-7.4655323310422510e+06, new int[] { 19, 1 });
                p.AddCoeff(1.4575563122511062e+06, new int[] { 21, 1 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{55BDD4BA-A476-41FC-948C-22606130853A}"));
                OrthonormalPolynomials[275] = p;
                p.AddCoeff(-5.6412002045046031e-01, new int[] { 0, 0 });
                p.AddCoeff(1.4272236517396646e+02, new int[] { 2, 0 });
                p.AddCoeff(-5.9467652155819357e+03, new int[] { 4, 0 });
                p.AddCoeff(9.6337596492427359e+04, new int[] { 6, 0 });
                p.AddCoeff(-7.9822579950868383e+05, new int[] { 8, 0 });
                p.AddCoeff(3.8492221887418754e+06, new int[] { 10, 0 });
                p.AddCoeff(-1.1547666566225626e+07, new int[] { 12, 0 });
                p.AddCoeff(2.2207051088895435e+07, new int[] { 14, 0 });
                p.AddCoeff(-2.7388696342971036e+07, new int[] { 16, 0 });
                p.AddCoeff(2.0944297203448439e+07, new int[] { 18, 0 });
                p.AddCoeff(-9.0391177404356423e+06, new int[] { 20, 0 });
                p.AddCoeff(1.6826063326352061e+06, new int[] { 22, 0 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{BA94949D-D7CA-4036-8118-45AB87902CCD}"));
                OrthonormalPolynomials[276] = p;
                p.AddCoeff(-1.3259954110337796e+01, new int[] { 0, 1 });
                p.AddCoeff(1.2154957934476313e+03, new int[] { 0, 3 });
                p.AddCoeff(-3.2818386423086044e+04, new int[] { 0, 5 });
                p.AddCoeff(4.0788565982978369e+05, new int[] { 0, 7 });
                p.AddCoeff(-2.8098789899385099e+06, new int[] { 0, 9 });
                p.AddCoeff(1.1801491757741742e+07, new int[] { 0, 11 });
                p.AddCoeff(-3.1773247040073919e+07, new int[] { 0, 13 });
                p.AddCoeff(5.5981435261082620e+07, new int[] { 0, 15 });
                p.AddCoeff(-6.4213999270065358e+07, new int[] { 0, 17 });
                p.AddCoeff(4.6189017018818942e+07, new int[] { 0, 19 });
                p.AddCoeff(-1.8915502207706805e+07, new int[] { 0, 21 });
                p.AddCoeff(3.3644173887225542e+06, new int[] { 0, 23 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{316E8993-C34E-4A74-AF26-2476D7C05EAA}"));
                OrthonormalPolynomials[277] = p;
                p.AddCoeff(-9.7708453698699135e-01, new int[] { 1, 0 });
                p.AddCoeff(2.4720238785770881e+02, new int[] { 1, 2 });
                p.AddCoeff(-1.0300099494071200e+04, new int[] { 1, 4 });
                p.AddCoeff(1.6686161180395345e+05, new int[] { 1, 6 });
                p.AddCoeff(-1.3825676406613286e+06, new int[] { 1, 8 });
                p.AddCoeff(6.6670484005224066e+06, new int[] { 1, 10 });
                p.AddCoeff(-2.0001145201567220e+07, new int[] { 1, 12 });
                p.AddCoeff(3.8463740772244654e+07, new int[] { 1, 14 });
                p.AddCoeff(-4.7438613619101740e+07, new int[] { 1, 16 });
                p.AddCoeff(3.6276586885195448e+07, new int[] { 1, 18 });
                p.AddCoeff(-1.5656211182031720e+07, new int[] { 1, 20 });
                p.AddCoeff(2.9143596572613158e+06, new int[] { 1, 22 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{FB35DE56-6D02-482B-949D-EFAE7BA5211B}"));
                OrthonormalPolynomials[278] = p;
                p.AddCoeff(-1.3563668632916897e+01, new int[] { 0, 1 });
                p.AddCoeff(1.0398812618569621e+03, new int[] { 0, 3 });
                p.AddCoeff(4.0691005898750690e+01, new int[] { 2, 1 });
                p.AddCoeff(-2.3397328391781646e+04, new int[] { 0, 5 });
                p.AddCoeff(-3.1196437855708862e+03, new int[] { 2, 3 });
                p.AddCoeff(2.4065823488689694e+05, new int[] { 0, 7 });
                p.AddCoeff(7.0191985175344939e+04, new int[] { 2, 5 });
                p.AddCoeff(-1.3570450467233355e+06, new int[] { 0, 9 });
                p.AddCoeff(-7.2197470466069081e+05, new int[] { 2, 7 });
                p.AddCoeff(4.5892796125552800e+06, new int[] { 0, 11 });
                p.AddCoeff(4.0711351401700065e+06, new int[] { 2, 9 });
                p.AddCoeff(-9.7080914880977078e+06, new int[] { 0, 13 });
                p.AddCoeff(-1.3767838837665840e+07, new int[] { 2, 11 });
                p.AddCoeff(1.2944121984130277e+07, new int[] { 0, 15 });
                p.AddCoeff(2.9124274464293123e+07, new int[] { 2, 13 });
                p.AddCoeff(-1.0564687795871035e+07, new int[] { 0, 17 });
                p.AddCoeff(-3.8832365952390831e+07, new int[] { 2, 15 });
                p.AddCoeff(4.8189803981166124e+06, new int[] { 0, 19 });
                p.AddCoeff(3.1694063387613105e+07, new int[] { 2, 17 });
                p.AddCoeff(-9.4084855391800528e+05, new int[] { 0, 21 });
                p.AddCoeff(-1.4456941194349837e+07, new int[] { 2, 19 });
                p.AddCoeff(2.8225456617540158e+06, new int[] { 2, 21 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{CDEB1411-7F26-43B5-87FB-AD1E05DF9263}"));
                OrthonormalPolynomials[279] = p;
                p.AddCoeff(-2.2387255181462104e+00, new int[] { 1, 0 });
                p.AddCoeff(4.7013235881070418e+02, new int[] { 1, 2 });
                p.AddCoeff(3.7312091969103506e+00, new int[] { 3, 0 });
                p.AddCoeff(-1.6219566378969294e+04, new int[] { 1, 4 });
                p.AddCoeff(-7.8355393135117363e+02, new int[] { 3, 2 });
                p.AddCoeff(2.1626088505292392e+05, new int[] { 1, 6 });
                p.AddCoeff(2.7032610631615490e+04, new int[] { 3, 4 });
                p.AddCoeff(-1.4597609741072365e+06, new int[] { 1, 8 });
                p.AddCoeff(-3.6043480842153987e+05, new int[] { 3, 6 });
                p.AddCoeff(5.6444090998813143e+06, new int[] { 1, 10 });
                p.AddCoeff(2.4329349568453941e+06, new int[] { 3, 8 });
                p.AddCoeff(-1.3255809249721269e+07, new int[] { 1, 12 });
                p.AddCoeff(-9.4073484998021906e+06, new int[] { 3, 10 });
                p.AddCoeff(1.9228206823771510e+07, new int[] { 1, 14 });
                p.AddCoeff(2.2093015416202114e+07, new int[] { 3, 12 });
                p.AddCoeff(-1.6824680970800072e+07, new int[] { 1, 16 });
                p.AddCoeff(-3.2047011372952517e+07, new int[] { 3, 14 });
                p.AddCoeff(8.1374273976418647e+06, new int[] { 1, 18 });
                p.AddCoeff(2.8041134951333453e+07, new int[] { 3, 16 });
                p.AddCoeff(-1.6703140447791196e+06, new int[] { 1, 20 });
                p.AddCoeff(-1.3562378996069774e+07, new int[] { 3, 18 });
                p.AddCoeff(2.7838567412985327e+06, new int[] { 3, 20 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{E7640AF9-668D-46E8-85FC-FFEF0E34678D}"));
                OrthonormalPolynomials[280] = p;
                p.AddCoeff(-1.2378940167103827e+01, new int[] { 0, 1 });
                p.AddCoeff(7.7987323052754111e+02, new int[] { 0, 3 });
                p.AddCoeff(1.2378940167103827e+02, new int[] { 2, 1 });
                p.AddCoeff(-1.4349667441706756e+04, new int[] { 0, 5 });
                p.AddCoeff(-7.7987323052754111e+03, new int[] { 2, 3 });
                p.AddCoeff(-1.4442096861621132e+02, new int[] { 4, 1 });
                p.AddCoeff(1.1958056201422297e+05, new int[] { 0, 7 });
                p.AddCoeff(1.4349667441706756e+05, new int[] { 2, 5 });
                p.AddCoeff(9.0985210228213130e+03, new int[] { 4, 3 });
                p.AddCoeff(-5.3811252906400337e+05, new int[] { 0, 9 });
                p.AddCoeff(-1.1958056201422297e+06, new int[] { 2, 7 });
                p.AddCoeff(-1.6741278681991216e+05, new int[] { 4, 5 });
                p.AddCoeff(1.4186603038960089e+06, new int[] { 0, 11 });
                p.AddCoeff(5.3811252906400337e+06, new int[] { 2, 9 });
                p.AddCoeff(1.3951065568326013e+06, new int[] { 4, 7 });
                p.AddCoeff(-2.2553061241423731e+06, new int[] { 0, 13 });
                p.AddCoeff(-1.4186603038960089e+07, new int[] { 2, 11 });
                p.AddCoeff(-6.2779795057467059e+06, new int[] { 4, 9 });
                p.AddCoeff(2.1264314884770946e+06, new int[] { 0, 15 });
                p.AddCoeff(2.2553061241423731e+07, new int[] { 2, 13 });
                p.AddCoeff(1.6551036878786770e+07, new int[] { 4, 11 });
                p.AddCoeff(-1.0944867955396811e+06, new int[] { 0, 17 });
                p.AddCoeff(-2.1264314884770946e+07, new int[] { 2, 15 });
                p.AddCoeff(-2.6311904781661019e+07, new int[] { 4, 13 });
                p.AddCoeff(2.3681878032145146e+05, new int[] { 0, 19 });
                p.AddCoeff(1.0944867955396811e+07, new int[] { 2, 17 });
                p.AddCoeff(2.4808367365566104e+07, new int[] { 4, 15 });
                p.AddCoeff(-2.3681878032145146e+06, new int[] { 2, 19 });
                p.AddCoeff(-1.2769012614629612e+07, new int[] { 4, 17 });
                p.AddCoeff(2.7628857704169337e+06, new int[] { 4, 19 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{0DB19280-6188-4723-920C-3DFF15ADE8B1}"));
                OrthonormalPolynomials[281] = p;
                p.AddCoeff(-3.5078701883878918e+00, new int[] { 1, 0 });
                p.AddCoeff(5.9984580221432950e+02, new int[] { 1, 2 });
                p.AddCoeff(1.6370060879143495e+01, new int[] { 3, 0 });
                p.AddCoeff(-1.6795682462001226e+04, new int[] { 1, 4 });
                p.AddCoeff(-2.7992804103335377e+03, new int[] { 3, 2 });
                p.AddCoeff(-1.4733054791229146e+01, new int[] { 5, 0 });
                p.AddCoeff(1.8027365842547982e+05, new int[] { 1, 6 });
                p.AddCoeff(7.8379851489339054e+04, new int[] { 3, 4 });
                p.AddCoeff(2.5193523693001839e+03, new int[] { 5, 2 });
                p.AddCoeff(-9.6575174156507049e+05, new int[] { 1, 8 });
                p.AddCoeff(-8.4127707265223918e+05, new int[] { 3, 6 });
                p.AddCoeff(-7.0541866340405149e+04, new int[] { 5, 4 });
                p.AddCoeff(2.8972552246952115e+06, new int[] { 1, 10 });
                p.AddCoeff(4.5068414606369956e+06, new int[] { 3, 8 });
                p.AddCoeff(7.5714936538701526e+05, new int[] { 5, 6 });
                p.AddCoeff(-5.0921455464340080e+06, new int[] { 1, 12 });
                p.AddCoeff(-1.3520524381910987e+07, new int[] { 3, 10 });
                p.AddCoeff(-4.0561573145732961e+06, new int[] { 5, 8 });
                p.AddCoeff(5.2040608331688214e+06, new int[] { 1, 14 });
                p.AddCoeff(2.3763345883358704e+07, new int[] { 3, 12 });
                p.AddCoeff(1.2168471943719888e+07, new int[] { 5, 10 });
                p.AddCoeff(-2.8622334582428518e+06, new int[] { 1, 16 });
                p.AddCoeff(-2.4285617221454500e+07, new int[] { 3, 14 });
                p.AddCoeff(-2.1387011295022834e+07, new int[] { 5, 12 });
                p.AddCoeff(6.5475928783333210e+05, new int[] { 1, 18 });
                p.AddCoeff(1.3357089471799975e+07, new int[] { 3, 16 });
                p.AddCoeff(2.1857055499309050e+07, new int[] { 5, 14 });
                p.AddCoeff(-3.0555433432222165e+06, new int[] { 3, 18 });
                p.AddCoeff(-1.2021380524619977e+07, new int[] { 5, 16 });
                p.AddCoeff(2.7499890088999948e+06, new int[] { 5, 18 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{E4FCC893-25E7-47F0-A0CD-6375B77D86A2}"));
                OrthonormalPolynomials[282] = p;
                p.AddCoeff(-1.1126876353251982e+01, new int[] { 0, 1 });
                p.AddCoeff(5.6376173523143378e+02, new int[] { 0, 3 });
                p.AddCoeff(2.3366440341829163e+02, new int[] { 2, 1 });
                p.AddCoeff(-8.2872975079020766e+03, new int[] { 0, 5 });
                p.AddCoeff(-1.1838996439860109e+04, new int[] { 2, 3 });
                p.AddCoeff(-7.0099321025487490e+02, new int[] { 4, 1 });
                p.AddCoeff(5.4459383623356503e+04, new int[] { 0, 7 });
                p.AddCoeff(1.7403324766594361e+05, new int[] { 2, 5 });
                p.AddCoeff(3.5516989319580328e+04, new int[] { 4, 3 });
                p.AddCoeff(5.1406168752024159e+02, new int[] { 6, 1 });
                p.AddCoeff(-1.8909508202554341e+05, new int[] { 0, 9 });
                p.AddCoeff(-1.1436470560904866e+06, new int[] { 2, 7 });
                p.AddCoeff(-5.2209974299783082e+05, new int[] { 4, 5 });
                p.AddCoeff(-2.6045792167692241e+04, new int[] { 6, 3 });
                p.AddCoeff(3.7131397925015798e+05, new int[] { 0, 11 });
                p.AddCoeff(3.9709967225364117e+06, new int[] { 2, 9 });
                p.AddCoeff(3.4309411682714597e+06, new int[] { 4, 7 });
                p.AddCoeff(3.8287314486507594e+05, new int[] { 6, 5 });
                p.AddCoeff(-4.1415789993286851e+05, new int[] { 0, 13 });
                p.AddCoeff(-7.7975935642533175e+06, new int[] { 2, 11 });
                p.AddCoeff(-1.1912990167609235e+07, new int[] { 4, 9 });
                p.AddCoeff(-2.5160235233990704e+06, new int[] { 6, 7 });
                p.AddCoeff(2.4455037900797950e+05, new int[] { 0, 15 });
                p.AddCoeff(8.6973158985902387e+06, new int[] { 2, 13 });
                p.AddCoeff(2.3392780692759952e+07, new int[] { 4, 11 });
                p.AddCoeff(8.7361927895801057e+06, new int[] { 6, 9 });
                p.AddCoeff(-5.9339430200465614e+04, new int[] { 0, 17 });
                p.AddCoeff(-5.1355579591675695e+06, new int[] { 2, 15 });
                p.AddCoeff(-2.6091947695770716e+07, new int[] { 4, 13 });
                p.AddCoeff(-1.7154705841357298e+07, new int[] { 6, 11 });
                p.AddCoeff(1.2461280342097779e+06, new int[] { 2, 17 });
                p.AddCoeff(1.5406673877502709e+07, new int[] { 4, 15 });
                p.AddCoeff(1.9134094976898525e+07, new int[] { 6, 13 });
                p.AddCoeff(-3.7383841026293337e+06, new int[] { 4, 17 });
                p.AddCoeff(-1.1298227510168653e+07, new int[] { 6, 15 });
                p.AddCoeff(2.7414816752615114e+06, new int[] { 6, 17 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{22676291-22B0-4085-9DC4-C5A0B106DC97}"));
                OrthonormalPolynomials[283] = p;
                p.AddCoeff(-4.7788046958621617e+00, new int[] { 1, 0 });
                p.AddCoeff(6.4991743863725399e+02, new int[] { 1, 2 });
                p.AddCoeff(4.3009242262759455e+01, new int[] { 3, 0 });
                p.AddCoeff(-1.4406503223125797e+04, new int[] { 1, 4 });
                p.AddCoeff(-5.8492569477352859e+03, new int[] { 3, 2 });
                p.AddCoeff(-9.4620332978070802e+01, new int[] { 5, 0 });
                p.AddCoeff(1.2101462707425669e+05, new int[] { 1, 6 });
                p.AddCoeff(1.2965852900813217e+05, new int[] { 3, 4 });
                p.AddCoeff(1.2868365285017629e+04, new int[] { 5, 2 });
                p.AddCoeff(5.8574491843567639e+01, new int[] { 7, 0 });
                p.AddCoeff(-4.9702436119783999e+05, new int[] { 1, 8 });
                p.AddCoeff(-1.0891316436683102e+06, new int[] { 3, 6 });
                p.AddCoeff(-2.8524876381789078e+05, new int[] { 5, 4 });
                p.AddCoeff(-7.9661308907251989e+03, new int[] { 7, 2 });
                p.AddCoeff(1.1044985804396444e+06, new int[] { 1, 10 });
                p.AddCoeff(4.4732192507805599e+06, new int[] { 3, 8 });
                p.AddCoeff(2.3960896160702825e+06, new int[] { 5, 6 });
                p.AddCoeff(1.7658256807774191e+05, new int[] { 7, 4 });
                p.AddCoeff(-1.3555209850850182e+06, new int[] { 1, 12 });
                p.AddCoeff(-9.9404872239567998e+06, new int[] { 3, 10 });
                p.AddCoeff(-9.8410823517172318e+06, new int[] { 5, 8 });
                p.AddCoeff(-1.4832935718530320e+06, new int[] { 7, 6 });
                p.AddCoeff(8.6395843005418739e+05, new int[] { 1, 14 });
                p.AddCoeff(1.2199688865765163e+07, new int[] { 3, 12 });
                p.AddCoeff(2.1869071892704960e+07, new int[] { 5, 10 });
                p.AddCoeff(6.0920985986820959e+06, new int[] { 7, 8 });
                p.AddCoeff(-2.2318926109733174e+05, new int[] { 1, 16 });
                p.AddCoeff(-7.7756258704876865e+06, new int[] { 3, 14 });
                p.AddCoeff(-2.6839315504683359e+07, new int[] { 5, 12 });
                p.AddCoeff(-1.3537996885960213e+07, new int[] { 7, 10 });
                p.AddCoeff(2.0087033498759857e+06, new int[] { 3, 16 });
                p.AddCoeff(1.7106376915072910e+07, new int[] { 5, 14 });
                p.AddCoeff(1.6614814360042080e+07, new int[] { 7, 12 });
                p.AddCoeff(-4.4191473697271685e+06, new int[] { 5, 16 });
                p.AddCoeff(-1.0589661899807040e+07, new int[] { 7, 14 });
                p.AddCoeff(2.7356626574501519e+06, new int[] { 7, 16 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{DAAF42CF-4AA3-44FA-87FF-73C583C47DE8}"));
                OrthonormalPolynomials[284] = p;
                p.AddCoeff(-9.8617045127668579e+00, new int[] { 0, 1 });
                p.AddCoeff(3.9118094567308536e+02, new int[] { 0, 3 });
                p.AddCoeff(3.5502136245960689e+02, new int[] { 2, 1 });
                p.AddCoeff(-4.4594627806731732e+03, new int[] { 0, 5 });
                p.AddCoeff(-1.4082514044231073e+04, new int[] { 2, 3 });
                p.AddCoeff(-1.9526174935278379e+03, new int[] { 4, 1 });
                p.AddCoeff(2.2297313903365866e+04, new int[] { 0, 7 });
                p.AddCoeff(1.6054066010423423e+05, new int[] { 2, 5 });
                p.AddCoeff(7.7453827243270902e+04, new int[] { 4, 3 });
                p.AddCoeff(3.3845369887815856e+03, new int[] { 6, 1 });
                p.AddCoeff(-5.6982024419712768e+04, new int[] { 0, 9 });
                p.AddCoeff(-8.0270330052117117e+05, new int[] { 2, 7 });
                p.AddCoeff(-8.8297363057328828e+05, new int[] { 4, 5 });
                p.AddCoeff(-1.3425330055500290e+05, new int[] { 6, 3 });
                p.AddCoeff(-1.8131448154187066e+03, new int[] { 8, 1 });
                p.AddCoeff(7.7702760572335593e+04, new int[] { 0, 11 });
                p.AddCoeff(2.0513528791096597e+06, new int[] { 2, 9 });
                p.AddCoeff(4.4148681528664414e+06, new int[] { 4, 7 });
                p.AddCoeff(1.5304876263270330e+06, new int[] { 6, 5 });
                p.AddCoeff(7.1921411011608695e+04, new int[] { 8, 3 });
                p.AddCoeff(-5.3794218857770795e+04, new int[] { 0, 13 });
                p.AddCoeff(-2.7972993806040813e+06, new int[] { 2, 11 });
                p.AddCoeff(-1.1282440835103128e+07, new int[] { 4, 9 });
                p.AddCoeff(-7.6524381316351651e+06, new int[] { 6, 7 });
                p.AddCoeff(-8.1990408553233912e+05, new int[] { 8, 5 });
                p.AddCoeff(1.4857450922622410e+04, new int[] { 0, 15 });
                p.AddCoeff(1.9365918788797486e+06, new int[] { 2, 13 });
                p.AddCoeff(1.5385146593322447e+07, new int[] { 4, 11 });
                p.AddCoeff(1.9556230780845422e+07, new int[] { 6, 9 });
                p.AddCoeff(4.0995204276616956e+06, new int[] { 8, 7 });
                p.AddCoeff(-5.3486823321440676e+05, new int[] { 2, 15 });
                p.AddCoeff(-1.0651255333838617e+07, new int[] { 4, 13 });
                p.AddCoeff(-2.6667587428425575e+07, new int[] { 6, 11 });
                p.AddCoeff(-1.0476552204024333e+07, new int[] { 8, 9 });
                p.AddCoeff(2.9417752826792372e+06, new int[] { 4, 15 });
                p.AddCoeff(1.8462175911986937e+07, new int[] { 6, 13 });
                p.AddCoeff(1.4286207550942273e+07, new int[] { 8, 11 });
                p.AddCoeff(-5.0990771566440111e+06, new int[] { 6, 15 });
                p.AddCoeff(-9.8904513814215733e+06, new int[] { 8, 13 });
                p.AddCoeff(2.7316484767735774e+06, new int[] { 8, 15 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{0BDB2C20-2950-4091-B9C1-DAB20199F872}"));
                OrthonormalPolynomials[285] = p;
                p.AddCoeff(-6.0502556762751543e+00, new int[] { 1, 0 });
                p.AddCoeff(6.3527684600889121e+02, new int[] { 1, 2 });
                p.AddCoeff(8.8737083252035597e+01, new int[] { 3, 0 });
                p.AddCoeff(-1.0799706382151150e+04, new int[] { 1, 4 });
                p.AddCoeff(-9.3173937414637377e+03, new int[] { 3, 2 });
                p.AddCoeff(-3.4607462468293883e+02, new int[] { 5, 0 });
                p.AddCoeff(6.8398140420290620e+04, new int[] { 1, 6 });
                p.AddCoeff(1.5839569360488354e+05, new int[] { 3, 4 });
                p.AddCoeff(3.6337835591708577e+04, new int[] { 5, 2 });
                p.AddCoeff(4.9439232097562690e+02, new int[] { 7, 0 });
                p.AddCoeff(-2.0519442126087186e+05, new int[] { 1, 8 });
                p.AddCoeff(-1.0031727261642624e+06, new int[] { 3, 6 });
                p.AddCoeff(-6.1774320505904581e+05, new int[] { 5, 4 });
                p.AddCoeff(-5.1911193702440824e+04, new int[] { 7, 2 });
                p.AddCoeff(-2.3346304046071270e+02, new int[] { 9, 0 });
                p.AddCoeff(3.1463144593333685e+05, new int[] { 1, 10 });
                p.AddCoeff(3.0095181784927873e+06, new int[] { 3, 8 });
                p.AddCoeff(3.9123736320406234e+06, new int[] { 5, 6 });
                p.AddCoeff(8.8249029294149401e+05, new int[] { 7, 4 });
                p.AddCoeff(2.4513619248374834e+04, new int[] { 9, 2 });
                p.AddCoeff(-2.3835715601010367e+05, new int[] { 1, 12 });
                p.AddCoeff(-4.6145945403556071e+06, new int[] { 3, 10 });
                p.AddCoeff(-1.1737120896121870e+07, new int[] { 5, 8 });
                p.AddCoeff(-5.5891051886294621e+06, new int[] { 7, 6 });
                p.AddCoeff(-4.1673152722237217e+05, new int[] { 9, 4 });
                p.AddCoeff(7.0721353981019772e+04, new int[] { 1, 14 });
                p.AddCoeff(3.4959049548148539e+06, new int[] { 3, 12 });
                p.AddCoeff(1.7996918707386868e+07, new int[] { 5, 10 });
                p.AddCoeff(1.6767315565888386e+07, new int[] { 7, 8 });
                p.AddCoeff(2.6392996724083571e+06, new int[] { 9, 6 });
                p.AddCoeff(-1.0372465250549567e+06, new int[] { 3, 14 });
                p.AddCoeff(-1.3634029323777930e+07, new int[] { 5, 12 });
                p.AddCoeff(-2.5709883867695526e+07, new int[] { 7, 10 });
                p.AddCoeff(-7.9178990172250713e+06, new int[] { 9, 8 });
                p.AddCoeff(4.0452614477143309e+06, new int[] { 5, 14 });
                p.AddCoeff(1.9477184748254186e+07, new int[] { 7, 12 });
                p.AddCoeff(1.2140778493078443e+07, new int[] { 9, 10 });
                p.AddCoeff(-5.7789449253061871e+06, new int[] { 7, 14 });
                p.AddCoeff(-9.1975594644533656e+06, new int[] { 9, 12 });
                p.AddCoeff(2.7289462147279217e+06, new int[] { 9, 14 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{E3719162-6004-40E5-9073-BB4758F02724}"));
                OrthonormalPolynomials[286] = p;
                p.AddCoeff(-8.5924594938297071e+00, new int[] { 0, 1 });
                p.AddCoeff(2.5777378481489121e+02, new int[] { 0, 3 });
                p.AddCoeff(4.7258527216063389e+02, new int[] { 2, 1 });
                p.AddCoeff(-2.1910771709265753e+03, new int[] { 0, 5 });
                p.AddCoeff(-1.4177558164819017e+04, new int[] { 2, 3 });
                p.AddCoeff(-4.0957390253921604e+03, new int[] { 4, 1 });
                p.AddCoeff(7.9296126185914154e+03, new int[] { 0, 7 });
                p.AddCoeff(1.2050924440096164e+05, new int[] { 2, 5 });
                p.AddCoeff(1.2287217076176481e+05, new int[] { 4, 3 });
                p.AddCoeff(1.2287217076176481e+04, new int[] { 6, 1 });
                p.AddCoeff(-1.3876822082534977e+04, new int[] { 0, 9 });
                p.AddCoeff(-4.3612869402252785e+05, new int[] { 2, 7 });
                p.AddCoeff(-1.0444134514750009e+06, new int[] { 4, 5 });
                p.AddCoeff(-3.6861651228529443e+05, new int[] { 6, 3 });
                p.AddCoeff(-1.4920192163928584e+04, new int[] { 8, 1 });
                p.AddCoeff(1.1606069378120163e+04, new int[] { 0, 11 });
                p.AddCoeff(7.6322521453942373e+05, new int[] { 2, 9 });
                p.AddCoeff(3.7797820148619080e+06, new int[] { 4, 7 });
                p.AddCoeff(3.1332403544250027e+06, new int[] { 6, 5 });
                p.AddCoeff(4.4760576491785753e+05, new int[] { 8, 3 });
                p.AddCoeff(6.2996366914365133e+03, new int[] { 10, 1 });
                p.AddCoeff(-3.7198940314487700e+03, new int[] { 0, 13 });
                p.AddCoeff(-6.3833381579660894e+05, new int[] { 2, 11 });
                p.AddCoeff(-6.6146185260083390e+06, new int[] { 4, 9 });
                p.AddCoeff(-1.1339346044585724e+07, new int[] { 6, 7 });
                p.AddCoeff(-3.8046490018017890e+06, new int[] { 8, 5 });
                p.AddCoeff(-1.8898910074309540e+05, new int[] { 10, 3 });
                p.AddCoeff(2.0459417172968235e+05, new int[] { 2, 13 });
                p.AddCoeff(5.5322264035706108e+06, new int[] { 4, 11 });
                p.AddCoeff(1.9843855578025017e+07, new int[] { 6, 9 });
                p.AddCoeff(1.3769205911282665e+07, new int[] { 8, 7 });
                p.AddCoeff(1.6064073563163109e+06, new int[] { 10, 5 });
                p.AddCoeff(-1.7731494883239137e+06, new int[] { 4, 13 });
                p.AddCoeff(-1.6596679210711832e+07, new int[] { 6, 11 });
                p.AddCoeff(-2.4096110344744664e+07, new int[] { 8, 9 });
                p.AddCoeff(-5.8136647180971252e+06, new int[] { 10, 7 });
                p.AddCoeff(5.3194484649717412e+06, new int[] { 6, 13 });
                p.AddCoeff(2.0153110470150082e+07, new int[] { 8, 11 });
                p.AddCoeff(1.0173913256669969e+07, new int[] { 10, 9 });
                p.AddCoeff(-6.4593302788942571e+06, new int[] { 8, 13 });
                p.AddCoeff(-8.5090910873967014e+06, new int[] { 10, 11 });
                p.AddCoeff(2.7272727844220197e+06, new int[] { 10, 13 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{478A5530-DA80-4B2A-AE27-3F25167F1CF5}"));
                OrthonormalPolynomials[287] = p;
                p.AddCoeff(-7.3216542982718476e+00, new int[] { 1, 0 });
                p.AddCoeff(5.7108903526520411e+02, new int[] { 1, 2 });
                p.AddCoeff(1.5863584312922336e+02, new int[] { 3, 0 });
                p.AddCoeff(-7.1386129408150514e+03, new int[] { 1, 4 });
                p.AddCoeff(-1.2373595764079422e+04, new int[] { 3, 2 });
                p.AddCoeff(-9.5181505877534018e+02, new int[] { 5, 0 });
                p.AddCoeff(3.2361711998361566e+04, new int[] { 1, 6 });
                p.AddCoeff(1.5466994705099278e+05, new int[] { 3, 4 });
                p.AddCoeff(7.4241574584476534e+04, new int[] { 5, 2 });
                p.AddCoeff(2.3115508570258262e+03, new int[] { 7, 0 });
                p.AddCoeff(-6.5879199425236045e+04, new int[] { 1, 8 });
                p.AddCoeff(-7.0117042663116727e+05, new int[] { 3, 6 });
                p.AddCoeff(-9.2801968230595668e+05, new int[] { 5, 4 });
                p.AddCoeff(-1.8030096684801444e+05, new int[] { 7, 2 });
                p.AddCoeff(-2.4399703490828165e+03, new int[] { 9, 0 });
                p.AddCoeff(6.1487252796886976e+04, new int[] { 1, 10 });
                p.AddCoeff(1.4273826542134477e+06, new int[] { 3, 8 });
                p.AddCoeff(4.2070225597870036e+06, new int[] { 5, 6 });
                p.AddCoeff(2.2537620856001805e+06, new int[] { 7, 4 });
                p.AddCoeff(1.9031768722845969e+05, new int[] { 9, 2 });
                p.AddCoeff(9.3162504237707539e+02, new int[] { 11, 0 });
                p.AddCoeff(-2.1427375974672734e+04, new int[] { 1, 12 });
                p.AddCoeff(-1.3322238105992178e+06, new int[] { 3, 10 });
                p.AddCoeff(-8.5642959252806859e+06, new int[] { 5, 8 });
                p.AddCoeff(-1.0217054788054152e+07, new int[] { 7, 6 });
                p.AddCoeff(-2.3789710903557461e+06, new int[] { 9, 4 });
                p.AddCoeff(-7.2666753305411880e+04, new int[] { 11, 2 });
                p.AddCoeff(4.6425981278457590e+05, new int[] { 3, 12 });
                p.AddCoeff(7.9933428635953069e+06, new int[] { 5, 10 });
                p.AddCoeff(2.0799004389967380e+07, new int[] { 7, 8 });
                p.AddCoeff(1.0784668942946049e+07, new int[] { 9, 6 });
                p.AddCoeff(9.0833441631764851e+05, new int[] { 11, 4 });
                p.AddCoeff(-2.7855588767074554e+06, new int[] { 5, 12 });
                p.AddCoeff(-1.9412404097302888e+07, new int[] { 7, 10 });
                p.AddCoeff(-2.1954504633854457e+07, new int[] { 9, 8 });
                p.AddCoeff(-4.1177826873066732e+06, new int[] { 11, 6 });
                p.AddCoeff(6.7649287005752489e+06, new int[] { 7, 12 });
                p.AddCoeff(2.0490870991597493e+07, new int[] { 9, 10 });
                p.AddCoeff(8.3826290420171562e+06, new int[] { 11, 8 });
                p.AddCoeff(-7.1407580728294294e+06, new int[] { 9, 12 });
                p.AddCoeff(-7.8237871058826791e+06, new int[] { 11, 10 });
                p.AddCoeff(2.7264712641712367e+06, new int[] { 11, 12 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{516B78AB-0FF6-4BD2-A3C4-8CF022C5D848}"));
                OrthonormalPolynomials[288] = p;
                p.AddCoeff(-7.3216542982718476e+00, new int[] { 0, 1 });
                p.AddCoeff(1.5863584312922336e+02, new int[] { 0, 3 });
                p.AddCoeff(5.7108903526520411e+02, new int[] { 2, 1 });
                p.AddCoeff(-9.5181505877534018e+02, new int[] { 0, 5 });
                p.AddCoeff(-1.2373595764079422e+04, new int[] { 2, 3 });
                p.AddCoeff(-7.1386129408150514e+03, new int[] { 4, 1 });
                p.AddCoeff(2.3115508570258262e+03, new int[] { 0, 7 });
                p.AddCoeff(7.4241574584476534e+04, new int[] { 2, 5 });
                p.AddCoeff(1.5466994705099278e+05, new int[] { 4, 3 });
                p.AddCoeff(3.2361711998361566e+04, new int[] { 6, 1 });
                p.AddCoeff(-2.4399703490828165e+03, new int[] { 0, 9 });
                p.AddCoeff(-1.8030096684801444e+05, new int[] { 2, 7 });
                p.AddCoeff(-9.2801968230595668e+05, new int[] { 4, 5 });
                p.AddCoeff(-7.0117042663116727e+05, new int[] { 6, 3 });
                p.AddCoeff(-6.5879199425236045e+04, new int[] { 8, 1 });
                p.AddCoeff(9.3162504237707539e+02, new int[] { 0, 11 });
                p.AddCoeff(1.9031768722845969e+05, new int[] { 2, 9 });
                p.AddCoeff(2.2537620856001805e+06, new int[] { 4, 7 });
                p.AddCoeff(4.2070225597870036e+06, new int[] { 6, 5 });
                p.AddCoeff(1.4273826542134477e+06, new int[] { 8, 3 });
                p.AddCoeff(6.1487252796886976e+04, new int[] { 10, 1 });
                p.AddCoeff(-7.2666753305411880e+04, new int[] { 2, 11 });
                p.AddCoeff(-2.3789710903557461e+06, new int[] { 4, 9 });
                p.AddCoeff(-1.0217054788054152e+07, new int[] { 6, 7 });
                p.AddCoeff(-8.5642959252806859e+06, new int[] { 8, 5 });
                p.AddCoeff(-1.3322238105992178e+06, new int[] { 10, 3 });
                p.AddCoeff(-2.1427375974672734e+04, new int[] { 12, 1 });
                p.AddCoeff(9.0833441631764851e+05, new int[] { 4, 11 });
                p.AddCoeff(1.0784668942946049e+07, new int[] { 6, 9 });
                p.AddCoeff(2.0799004389967380e+07, new int[] { 8, 7 });
                p.AddCoeff(7.9933428635953069e+06, new int[] { 10, 5 });
                p.AddCoeff(4.6425981278457590e+05, new int[] { 12, 3 });
                p.AddCoeff(-4.1177826873066732e+06, new int[] { 6, 11 });
                p.AddCoeff(-2.1954504633854457e+07, new int[] { 8, 9 });
                p.AddCoeff(-1.9412404097302888e+07, new int[] { 10, 7 });
                p.AddCoeff(-2.7855588767074554e+06, new int[] { 12, 5 });
                p.AddCoeff(8.3826290420171562e+06, new int[] { 8, 11 });
                p.AddCoeff(2.0490870991597493e+07, new int[] { 10, 9 });
                p.AddCoeff(6.7649287005752489e+06, new int[] { 12, 7 });
                p.AddCoeff(-7.8237871058826791e+06, new int[] { 10, 11 });
                p.AddCoeff(-7.1407580728294294e+06, new int[] { 12, 9 });
                p.AddCoeff(2.7264712641712367e+06, new int[] { 12, 11 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{A96C4D1B-5486-4A5D-A463-C9A61F2AC6F1}"));
                OrthonormalPolynomials[289] = p;
                p.AddCoeff(-8.5924594938297071e+00, new int[] { 1, 0 });
                p.AddCoeff(4.7258527216063389e+02, new int[] { 1, 2 });
                p.AddCoeff(2.5777378481489121e+02, new int[] { 3, 0 });
                p.AddCoeff(-4.0957390253921604e+03, new int[] { 1, 4 });
                p.AddCoeff(-1.4177558164819017e+04, new int[] { 3, 2 });
                p.AddCoeff(-2.1910771709265753e+03, new int[] { 5, 0 });
                p.AddCoeff(1.2287217076176481e+04, new int[] { 1, 6 });
                p.AddCoeff(1.2287217076176481e+05, new int[] { 3, 4 });
                p.AddCoeff(1.2050924440096164e+05, new int[] { 5, 2 });
                p.AddCoeff(7.9296126185914154e+03, new int[] { 7, 0 });
                p.AddCoeff(-1.4920192163928584e+04, new int[] { 1, 8 });
                p.AddCoeff(-3.6861651228529443e+05, new int[] { 3, 6 });
                p.AddCoeff(-1.0444134514750009e+06, new int[] { 5, 4 });
                p.AddCoeff(-4.3612869402252785e+05, new int[] { 7, 2 });
                p.AddCoeff(-1.3876822082534977e+04, new int[] { 9, 0 });
                p.AddCoeff(6.2996366914365133e+03, new int[] { 1, 10 });
                p.AddCoeff(4.4760576491785753e+05, new int[] { 3, 8 });
                p.AddCoeff(3.1332403544250027e+06, new int[] { 5, 6 });
                p.AddCoeff(3.7797820148619080e+06, new int[] { 7, 4 });
                p.AddCoeff(7.6322521453942373e+05, new int[] { 9, 2 });
                p.AddCoeff(1.1606069378120163e+04, new int[] { 11, 0 });
                p.AddCoeff(-1.8898910074309540e+05, new int[] { 3, 10 });
                p.AddCoeff(-3.8046490018017890e+06, new int[] { 5, 8 });
                p.AddCoeff(-1.1339346044585724e+07, new int[] { 7, 6 });
                p.AddCoeff(-6.6146185260083390e+06, new int[] { 9, 4 });
                p.AddCoeff(-6.3833381579660894e+05, new int[] { 11, 2 });
                p.AddCoeff(-3.7198940314487700e+03, new int[] { 13, 0 });
                p.AddCoeff(1.6064073563163109e+06, new int[] { 5, 10 });
                p.AddCoeff(1.3769205911282665e+07, new int[] { 7, 8 });
                p.AddCoeff(1.9843855578025017e+07, new int[] { 9, 6 });
                p.AddCoeff(5.5322264035706108e+06, new int[] { 11, 4 });
                p.AddCoeff(2.0459417172968235e+05, new int[] { 13, 2 });
                p.AddCoeff(-5.8136647180971252e+06, new int[] { 7, 10 });
                p.AddCoeff(-2.4096110344744664e+07, new int[] { 9, 8 });
                p.AddCoeff(-1.6596679210711832e+07, new int[] { 11, 6 });
                p.AddCoeff(-1.7731494883239137e+06, new int[] { 13, 4 });
                p.AddCoeff(1.0173913256669969e+07, new int[] { 9, 10 });
                p.AddCoeff(2.0153110470150082e+07, new int[] { 11, 8 });
                p.AddCoeff(5.3194484649717412e+06, new int[] { 13, 6 });
                p.AddCoeff(-8.5090910873967014e+06, new int[] { 11, 10 });
                p.AddCoeff(-6.4593302788942571e+06, new int[] { 13, 8 });
                p.AddCoeff(2.7272727844220197e+06, new int[] { 13, 10 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{8DEBF36F-7F53-4C5A-A46F-3223D8D00577}"));
                OrthonormalPolynomials[290] = p;
                p.AddCoeff(-6.0502556762751543e+00, new int[] { 0, 1 });
                p.AddCoeff(8.8737083252035597e+01, new int[] { 0, 3 });
                p.AddCoeff(6.3527684600889121e+02, new int[] { 2, 1 });
                p.AddCoeff(-3.4607462468293883e+02, new int[] { 0, 5 });
                p.AddCoeff(-9.3173937414637377e+03, new int[] { 2, 3 });
                p.AddCoeff(-1.0799706382151150e+04, new int[] { 4, 1 });
                p.AddCoeff(4.9439232097562690e+02, new int[] { 0, 7 });
                p.AddCoeff(3.6337835591708577e+04, new int[] { 2, 5 });
                p.AddCoeff(1.5839569360488354e+05, new int[] { 4, 3 });
                p.AddCoeff(6.8398140420290620e+04, new int[] { 6, 1 });
                p.AddCoeff(-2.3346304046071270e+02, new int[] { 0, 9 });
                p.AddCoeff(-5.1911193702440824e+04, new int[] { 2, 7 });
                p.AddCoeff(-6.1774320505904581e+05, new int[] { 4, 5 });
                p.AddCoeff(-1.0031727261642624e+06, new int[] { 6, 3 });
                p.AddCoeff(-2.0519442126087186e+05, new int[] { 8, 1 });
                p.AddCoeff(2.4513619248374834e+04, new int[] { 2, 9 });
                p.AddCoeff(8.8249029294149401e+05, new int[] { 4, 7 });
                p.AddCoeff(3.9123736320406234e+06, new int[] { 6, 5 });
                p.AddCoeff(3.0095181784927873e+06, new int[] { 8, 3 });
                p.AddCoeff(3.1463144593333685e+05, new int[] { 10, 1 });
                p.AddCoeff(-4.1673152722237217e+05, new int[] { 4, 9 });
                p.AddCoeff(-5.5891051886294621e+06, new int[] { 6, 7 });
                p.AddCoeff(-1.1737120896121870e+07, new int[] { 8, 5 });
                p.AddCoeff(-4.6145945403556071e+06, new int[] { 10, 3 });
                p.AddCoeff(-2.3835715601010367e+05, new int[] { 12, 1 });
                p.AddCoeff(2.6392996724083571e+06, new int[] { 6, 9 });
                p.AddCoeff(1.6767315565888386e+07, new int[] { 8, 7 });
                p.AddCoeff(1.7996918707386868e+07, new int[] { 10, 5 });
                p.AddCoeff(3.4959049548148539e+06, new int[] { 12, 3 });
                p.AddCoeff(7.0721353981019772e+04, new int[] { 14, 1 });
                p.AddCoeff(-7.9178990172250713e+06, new int[] { 8, 9 });
                p.AddCoeff(-2.5709883867695526e+07, new int[] { 10, 7 });
                p.AddCoeff(-1.3634029323777930e+07, new int[] { 12, 5 });
                p.AddCoeff(-1.0372465250549567e+06, new int[] { 14, 3 });
                p.AddCoeff(1.2140778493078443e+07, new int[] { 10, 9 });
                p.AddCoeff(1.9477184748254186e+07, new int[] { 12, 7 });
                p.AddCoeff(4.0452614477143309e+06, new int[] { 14, 5 });
                p.AddCoeff(-9.1975594644533656e+06, new int[] { 12, 9 });
                p.AddCoeff(-5.7789449253061871e+06, new int[] { 14, 7 });
                p.AddCoeff(2.7289462147279217e+06, new int[] { 14, 9 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{1FDE930D-AC0C-4F61-AAC6-7EB59843DF87}"));
                OrthonormalPolynomials[291] = p;
                p.AddCoeff(-9.8617045127668579e+00, new int[] { 1, 0 });
                p.AddCoeff(3.5502136245960689e+02, new int[] { 1, 2 });
                p.AddCoeff(3.9118094567308536e+02, new int[] { 3, 0 });
                p.AddCoeff(-1.9526174935278379e+03, new int[] { 1, 4 });
                p.AddCoeff(-1.4082514044231073e+04, new int[] { 3, 2 });
                p.AddCoeff(-4.4594627806731732e+03, new int[] { 5, 0 });
                p.AddCoeff(3.3845369887815856e+03, new int[] { 1, 6 });
                p.AddCoeff(7.7453827243270902e+04, new int[] { 3, 4 });
                p.AddCoeff(1.6054066010423423e+05, new int[] { 5, 2 });
                p.AddCoeff(2.2297313903365866e+04, new int[] { 7, 0 });
                p.AddCoeff(-1.8131448154187066e+03, new int[] { 1, 8 });
                p.AddCoeff(-1.3425330055500290e+05, new int[] { 3, 6 });
                p.AddCoeff(-8.8297363057328828e+05, new int[] { 5, 4 });
                p.AddCoeff(-8.0270330052117117e+05, new int[] { 7, 2 });
                p.AddCoeff(-5.6982024419712768e+04, new int[] { 9, 0 });
                p.AddCoeff(7.1921411011608695e+04, new int[] { 3, 8 });
                p.AddCoeff(1.5304876263270330e+06, new int[] { 5, 6 });
                p.AddCoeff(4.4148681528664414e+06, new int[] { 7, 4 });
                p.AddCoeff(2.0513528791096597e+06, new int[] { 9, 2 });
                p.AddCoeff(7.7702760572335593e+04, new int[] { 11, 0 });
                p.AddCoeff(-8.1990408553233912e+05, new int[] { 5, 8 });
                p.AddCoeff(-7.6524381316351651e+06, new int[] { 7, 6 });
                p.AddCoeff(-1.1282440835103128e+07, new int[] { 9, 4 });
                p.AddCoeff(-2.7972993806040813e+06, new int[] { 11, 2 });
                p.AddCoeff(-5.3794218857770795e+04, new int[] { 13, 0 });
                p.AddCoeff(4.0995204276616956e+06, new int[] { 7, 8 });
                p.AddCoeff(1.9556230780845422e+07, new int[] { 9, 6 });
                p.AddCoeff(1.5385146593322447e+07, new int[] { 11, 4 });
                p.AddCoeff(1.9365918788797486e+06, new int[] { 13, 2 });
                p.AddCoeff(1.4857450922622410e+04, new int[] { 15, 0 });
                p.AddCoeff(-1.0476552204024333e+07, new int[] { 9, 8 });
                p.AddCoeff(-2.6667587428425575e+07, new int[] { 11, 6 });
                p.AddCoeff(-1.0651255333838617e+07, new int[] { 13, 4 });
                p.AddCoeff(-5.3486823321440676e+05, new int[] { 15, 2 });
                p.AddCoeff(1.4286207550942273e+07, new int[] { 11, 8 });
                p.AddCoeff(1.8462175911986937e+07, new int[] { 13, 6 });
                p.AddCoeff(2.9417752826792372e+06, new int[] { 15, 4 });
                p.AddCoeff(-9.8904513814215733e+06, new int[] { 13, 8 });
                p.AddCoeff(-5.0990771566440111e+06, new int[] { 15, 6 });
                p.AddCoeff(2.7316484767735774e+06, new int[] { 15, 8 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{CA5DD791-E462-4C9A-A6F8-A9C4AD2C91D8}"));
                OrthonormalPolynomials[292] = p;
                p.AddCoeff(-4.7788046958621617e+00, new int[] { 0, 1 });
                p.AddCoeff(4.3009242262759455e+01, new int[] { 0, 3 });
                p.AddCoeff(6.4991743863725399e+02, new int[] { 2, 1 });
                p.AddCoeff(-9.4620332978070802e+01, new int[] { 0, 5 });
                p.AddCoeff(-5.8492569477352859e+03, new int[] { 2, 3 });
                p.AddCoeff(-1.4406503223125797e+04, new int[] { 4, 1 });
                p.AddCoeff(5.8574491843567639e+01, new int[] { 0, 7 });
                p.AddCoeff(1.2868365285017629e+04, new int[] { 2, 5 });
                p.AddCoeff(1.2965852900813217e+05, new int[] { 4, 3 });
                p.AddCoeff(1.2101462707425669e+05, new int[] { 6, 1 });
                p.AddCoeff(-7.9661308907251989e+03, new int[] { 2, 7 });
                p.AddCoeff(-2.8524876381789078e+05, new int[] { 4, 5 });
                p.AddCoeff(-1.0891316436683102e+06, new int[] { 6, 3 });
                p.AddCoeff(-4.9702436119783999e+05, new int[] { 8, 1 });
                p.AddCoeff(1.7658256807774191e+05, new int[] { 4, 7 });
                p.AddCoeff(2.3960896160702825e+06, new int[] { 6, 5 });
                p.AddCoeff(4.4732192507805599e+06, new int[] { 8, 3 });
                p.AddCoeff(1.1044985804396444e+06, new int[] { 10, 1 });
                p.AddCoeff(-1.4832935718530320e+06, new int[] { 6, 7 });
                p.AddCoeff(-9.8410823517172318e+06, new int[] { 8, 5 });
                p.AddCoeff(-9.9404872239567998e+06, new int[] { 10, 3 });
                p.AddCoeff(-1.3555209850850182e+06, new int[] { 12, 1 });
                p.AddCoeff(6.0920985986820959e+06, new int[] { 8, 7 });
                p.AddCoeff(2.1869071892704960e+07, new int[] { 10, 5 });
                p.AddCoeff(1.2199688865765163e+07, new int[] { 12, 3 });
                p.AddCoeff(8.6395843005418739e+05, new int[] { 14, 1 });
                p.AddCoeff(-1.3537996885960213e+07, new int[] { 10, 7 });
                p.AddCoeff(-2.6839315504683359e+07, new int[] { 12, 5 });
                p.AddCoeff(-7.7756258704876865e+06, new int[] { 14, 3 });
                p.AddCoeff(-2.2318926109733174e+05, new int[] { 16, 1 });
                p.AddCoeff(1.6614814360042080e+07, new int[] { 12, 7 });
                p.AddCoeff(1.7106376915072910e+07, new int[] { 14, 5 });
                p.AddCoeff(2.0087033498759857e+06, new int[] { 16, 3 });
                p.AddCoeff(-1.0589661899807040e+07, new int[] { 14, 7 });
                p.AddCoeff(-4.4191473697271685e+06, new int[] { 16, 5 });
                p.AddCoeff(2.7356626574501519e+06, new int[] { 16, 7 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{BF947DD8-0C16-47B9-9E33-B8349C322EB8}"));
                OrthonormalPolynomials[293] = p;
                p.AddCoeff(-1.1126876353251982e+01, new int[] { 1, 0 });
                p.AddCoeff(2.3366440341829163e+02, new int[] { 1, 2 });
                p.AddCoeff(5.6376173523143378e+02, new int[] { 3, 0 });
                p.AddCoeff(-7.0099321025487490e+02, new int[] { 1, 4 });
                p.AddCoeff(-1.1838996439860109e+04, new int[] { 3, 2 });
                p.AddCoeff(-8.2872975079020766e+03, new int[] { 5, 0 });
                p.AddCoeff(5.1406168752024159e+02, new int[] { 1, 6 });
                p.AddCoeff(3.5516989319580328e+04, new int[] { 3, 4 });
                p.AddCoeff(1.7403324766594361e+05, new int[] { 5, 2 });
                p.AddCoeff(5.4459383623356503e+04, new int[] { 7, 0 });
                p.AddCoeff(-2.6045792167692241e+04, new int[] { 3, 6 });
                p.AddCoeff(-5.2209974299783082e+05, new int[] { 5, 4 });
                p.AddCoeff(-1.1436470560904866e+06, new int[] { 7, 2 });
                p.AddCoeff(-1.8909508202554341e+05, new int[] { 9, 0 });
                p.AddCoeff(3.8287314486507594e+05, new int[] { 5, 6 });
                p.AddCoeff(3.4309411682714597e+06, new int[] { 7, 4 });
                p.AddCoeff(3.9709967225364117e+06, new int[] { 9, 2 });
                p.AddCoeff(3.7131397925015798e+05, new int[] { 11, 0 });
                p.AddCoeff(-2.5160235233990704e+06, new int[] { 7, 6 });
                p.AddCoeff(-1.1912990167609235e+07, new int[] { 9, 4 });
                p.AddCoeff(-7.7975935642533175e+06, new int[] { 11, 2 });
                p.AddCoeff(-4.1415789993286851e+05, new int[] { 13, 0 });
                p.AddCoeff(8.7361927895801057e+06, new int[] { 9, 6 });
                p.AddCoeff(2.3392780692759952e+07, new int[] { 11, 4 });
                p.AddCoeff(8.6973158985902387e+06, new int[] { 13, 2 });
                p.AddCoeff(2.4455037900797950e+05, new int[] { 15, 0 });
                p.AddCoeff(-1.7154705841357298e+07, new int[] { 11, 6 });
                p.AddCoeff(-2.6091947695770716e+07, new int[] { 13, 4 });
                p.AddCoeff(-5.1355579591675695e+06, new int[] { 15, 2 });
                p.AddCoeff(-5.9339430200465614e+04, new int[] { 17, 0 });
                p.AddCoeff(1.9134094976898525e+07, new int[] { 13, 6 });
                p.AddCoeff(1.5406673877502709e+07, new int[] { 15, 4 });
                p.AddCoeff(1.2461280342097779e+06, new int[] { 17, 2 });
                p.AddCoeff(-1.1298227510168653e+07, new int[] { 15, 6 });
                p.AddCoeff(-3.7383841026293337e+06, new int[] { 17, 4 });
                p.AddCoeff(2.7414816752615114e+06, new int[] { 17, 6 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{38F944EA-7567-4CEA-B5EB-B0EFEC6304F8}"));
                OrthonormalPolynomials[294] = p;
                p.AddCoeff(-3.5078701883878918e+00, new int[] { 0, 1 });
                p.AddCoeff(1.6370060879143495e+01, new int[] { 0, 3 });
                p.AddCoeff(5.9984580221432950e+02, new int[] { 2, 1 });
                p.AddCoeff(-1.4733054791229146e+01, new int[] { 0, 5 });
                p.AddCoeff(-2.7992804103335377e+03, new int[] { 2, 3 });
                p.AddCoeff(-1.6795682462001226e+04, new int[] { 4, 1 });
                p.AddCoeff(2.5193523693001839e+03, new int[] { 2, 5 });
                p.AddCoeff(7.8379851489339054e+04, new int[] { 4, 3 });
                p.AddCoeff(1.8027365842547982e+05, new int[] { 6, 1 });
                p.AddCoeff(-7.0541866340405149e+04, new int[] { 4, 5 });
                p.AddCoeff(-8.4127707265223918e+05, new int[] { 6, 3 });
                p.AddCoeff(-9.6575174156507049e+05, new int[] { 8, 1 });
                p.AddCoeff(7.5714936538701526e+05, new int[] { 6, 5 });
                p.AddCoeff(4.5068414606369956e+06, new int[] { 8, 3 });
                p.AddCoeff(2.8972552246952115e+06, new int[] { 10, 1 });
                p.AddCoeff(-4.0561573145732961e+06, new int[] { 8, 5 });
                p.AddCoeff(-1.3520524381910987e+07, new int[] { 10, 3 });
                p.AddCoeff(-5.0921455464340080e+06, new int[] { 12, 1 });
                p.AddCoeff(1.2168471943719888e+07, new int[] { 10, 5 });
                p.AddCoeff(2.3763345883358704e+07, new int[] { 12, 3 });
                p.AddCoeff(5.2040608331688214e+06, new int[] { 14, 1 });
                p.AddCoeff(-2.1387011295022834e+07, new int[] { 12, 5 });
                p.AddCoeff(-2.4285617221454500e+07, new int[] { 14, 3 });
                p.AddCoeff(-2.8622334582428518e+06, new int[] { 16, 1 });
                p.AddCoeff(2.1857055499309050e+07, new int[] { 14, 5 });
                p.AddCoeff(1.3357089471799975e+07, new int[] { 16, 3 });
                p.AddCoeff(6.5475928783333210e+05, new int[] { 18, 1 });
                p.AddCoeff(-1.2021380524619977e+07, new int[] { 16, 5 });
                p.AddCoeff(-3.0555433432222165e+06, new int[] { 18, 3 });
                p.AddCoeff(2.7499890088999948e+06, new int[] { 18, 5 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{544852D0-EA38-4A93-8DD5-8E00FA39608E}"));
                OrthonormalPolynomials[295] = p;
                p.AddCoeff(-1.2378940167103827e+01, new int[] { 1, 0 });
                p.AddCoeff(1.2378940167103827e+02, new int[] { 1, 2 });
                p.AddCoeff(7.7987323052754111e+02, new int[] { 3, 0 });
                p.AddCoeff(-1.4442096861621132e+02, new int[] { 1, 4 });
                p.AddCoeff(-7.7987323052754111e+03, new int[] { 3, 2 });
                p.AddCoeff(-1.4349667441706756e+04, new int[] { 5, 0 });
                p.AddCoeff(9.0985210228213130e+03, new int[] { 3, 4 });
                p.AddCoeff(1.4349667441706756e+05, new int[] { 5, 2 });
                p.AddCoeff(1.1958056201422297e+05, new int[] { 7, 0 });
                p.AddCoeff(-1.6741278681991216e+05, new int[] { 5, 4 });
                p.AddCoeff(-1.1958056201422297e+06, new int[] { 7, 2 });
                p.AddCoeff(-5.3811252906400337e+05, new int[] { 9, 0 });
                p.AddCoeff(1.3951065568326013e+06, new int[] { 7, 4 });
                p.AddCoeff(5.3811252906400337e+06, new int[] { 9, 2 });
                p.AddCoeff(1.4186603038960089e+06, new int[] { 11, 0 });
                p.AddCoeff(-6.2779795057467059e+06, new int[] { 9, 4 });
                p.AddCoeff(-1.4186603038960089e+07, new int[] { 11, 2 });
                p.AddCoeff(-2.2553061241423731e+06, new int[] { 13, 0 });
                p.AddCoeff(1.6551036878786770e+07, new int[] { 11, 4 });
                p.AddCoeff(2.2553061241423731e+07, new int[] { 13, 2 });
                p.AddCoeff(2.1264314884770946e+06, new int[] { 15, 0 });
                p.AddCoeff(-2.6311904781661019e+07, new int[] { 13, 4 });
                p.AddCoeff(-2.1264314884770946e+07, new int[] { 15, 2 });
                p.AddCoeff(-1.0944867955396811e+06, new int[] { 17, 0 });
                p.AddCoeff(2.4808367365566104e+07, new int[] { 15, 4 });
                p.AddCoeff(1.0944867955396811e+07, new int[] { 17, 2 });
                p.AddCoeff(2.3681878032145146e+05, new int[] { 19, 0 });
                p.AddCoeff(-1.2769012614629612e+07, new int[] { 17, 4 });
                p.AddCoeff(-2.3681878032145146e+06, new int[] { 19, 2 });
                p.AddCoeff(2.7628857704169337e+06, new int[] { 19, 4 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{6B379F8C-34B7-46CF-B88F-39119C3BC44C}"));
                OrthonormalPolynomials[296] = p;
                p.AddCoeff(-2.2387255181462104e+00, new int[] { 0, 1 });
                p.AddCoeff(3.7312091969103506e+00, new int[] { 0, 3 });
                p.AddCoeff(4.7013235881070418e+02, new int[] { 2, 1 });
                p.AddCoeff(-7.8355393135117363e+02, new int[] { 2, 3 });
                p.AddCoeff(-1.6219566378969294e+04, new int[] { 4, 1 });
                p.AddCoeff(2.7032610631615490e+04, new int[] { 4, 3 });
                p.AddCoeff(2.1626088505292392e+05, new int[] { 6, 1 });
                p.AddCoeff(-3.6043480842153987e+05, new int[] { 6, 3 });
                p.AddCoeff(-1.4597609741072365e+06, new int[] { 8, 1 });
                p.AddCoeff(2.4329349568453941e+06, new int[] { 8, 3 });
                p.AddCoeff(5.6444090998813143e+06, new int[] { 10, 1 });
                p.AddCoeff(-9.4073484998021906e+06, new int[] { 10, 3 });
                p.AddCoeff(-1.3255809249721269e+07, new int[] { 12, 1 });
                p.AddCoeff(2.2093015416202114e+07, new int[] { 12, 3 });
                p.AddCoeff(1.9228206823771510e+07, new int[] { 14, 1 });
                p.AddCoeff(-3.2047011372952517e+07, new int[] { 14, 3 });
                p.AddCoeff(-1.6824680970800072e+07, new int[] { 16, 1 });
                p.AddCoeff(2.8041134951333453e+07, new int[] { 16, 3 });
                p.AddCoeff(8.1374273976418647e+06, new int[] { 18, 1 });
                p.AddCoeff(-1.3562378996069774e+07, new int[] { 18, 3 });
                p.AddCoeff(-1.6703140447791196e+06, new int[] { 20, 1 });
                p.AddCoeff(2.7838567412985327e+06, new int[] { 20, 3 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{5AF069FB-4EE1-492E-B5FA-30A8C657C859}"));
                OrthonormalPolynomials[297] = p;
                p.AddCoeff(-1.3563668632916897e+01, new int[] { 1, 0 });
                p.AddCoeff(4.0691005898750690e+01, new int[] { 1, 2 });
                p.AddCoeff(1.0398812618569621e+03, new int[] { 3, 0 });
                p.AddCoeff(-3.1196437855708862e+03, new int[] { 3, 2 });
                p.AddCoeff(-2.3397328391781646e+04, new int[] { 5, 0 });
                p.AddCoeff(7.0191985175344939e+04, new int[] { 5, 2 });
                p.AddCoeff(2.4065823488689694e+05, new int[] { 7, 0 });
                p.AddCoeff(-7.2197470466069081e+05, new int[] { 7, 2 });
                p.AddCoeff(-1.3570450467233355e+06, new int[] { 9, 0 });
                p.AddCoeff(4.0711351401700065e+06, new int[] { 9, 2 });
                p.AddCoeff(4.5892796125552800e+06, new int[] { 11, 0 });
                p.AddCoeff(-1.3767838837665840e+07, new int[] { 11, 2 });
                p.AddCoeff(-9.7080914880977078e+06, new int[] { 13, 0 });
                p.AddCoeff(2.9124274464293123e+07, new int[] { 13, 2 });
                p.AddCoeff(1.2944121984130277e+07, new int[] { 15, 0 });
                p.AddCoeff(-3.8832365952390831e+07, new int[] { 15, 2 });
                p.AddCoeff(-1.0564687795871035e+07, new int[] { 17, 0 });
                p.AddCoeff(3.1694063387613105e+07, new int[] { 17, 2 });
                p.AddCoeff(4.8189803981166124e+06, new int[] { 19, 0 });
                p.AddCoeff(-1.4456941194349837e+07, new int[] { 19, 2 });
                p.AddCoeff(-9.4084855391800528e+05, new int[] { 21, 0 });
                p.AddCoeff(2.8225456617540158e+06, new int[] { 21, 2 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{40880373-2DF2-4DF3-9BBC-FD55BEE7D9E6}"));
                OrthonormalPolynomials[298] = p;
                p.AddCoeff(-9.7708453698699135e-01, new int[] { 0, 1 });
                p.AddCoeff(2.4720238785770881e+02, new int[] { 2, 1 });
                p.AddCoeff(-1.0300099494071200e+04, new int[] { 4, 1 });
                p.AddCoeff(1.6686161180395345e+05, new int[] { 6, 1 });
                p.AddCoeff(-1.3825676406613286e+06, new int[] { 8, 1 });
                p.AddCoeff(6.6670484005224066e+06, new int[] { 10, 1 });
                p.AddCoeff(-2.0001145201567220e+07, new int[] { 12, 1 });
                p.AddCoeff(3.8463740772244654e+07, new int[] { 14, 1 });
                p.AddCoeff(-4.7438613619101740e+07, new int[] { 16, 1 });
                p.AddCoeff(3.6276586885195448e+07, new int[] { 18, 1 });
                p.AddCoeff(-1.5656211182031720e+07, new int[] { 20, 1 });
                p.AddCoeff(2.9143596572613158e+06, new int[] { 22, 1 });

                //----------------------------------------------------------------------------------------
                p = new Polynomial(new Guid("{27193D53-8C3B-4990-B912-855BB7127257}"));
                OrthonormalPolynomials[299] = p;
                p.AddCoeff(-1.3259954110337796e+01, new int[] { 1, 0 });
                p.AddCoeff(1.2154957934476313e+03, new int[] { 3, 0 });
                p.AddCoeff(-3.2818386423086044e+04, new int[] { 5, 0 });
                p.AddCoeff(4.0788565982978369e+05, new int[] { 7, 0 });
                p.AddCoeff(-2.8098789899385099e+06, new int[] { 9, 0 });
                p.AddCoeff(1.1801491757741742e+07, new int[] { 11, 0 });
                p.AddCoeff(-3.1773247040073919e+07, new int[] { 13, 0 });
                p.AddCoeff(5.5981435261082620e+07, new int[] { 15, 0 });
                p.AddCoeff(-6.4213999270065358e+07, new int[] { 17, 0 });
                p.AddCoeff(4.6189017018818942e+07, new int[] { 19, 0 });
                p.AddCoeff(-1.8915502207706805e+07, new int[] { 21, 0 });
                p.AddCoeff(3.3644173887225542e+06, new int[] { 23, 0 });
            #endregion POLY_DEF
#pragma warning restore 612
            }
        }
        /// <summary>
        /// transforms some vertices (<paramref name="EdgeVertices"/>) from the local 1D-coordinate system of either
        /// the top, bottom, left or right edge to the local coordinate system of the square;
        /// </summary>
        /// <param name="EdgeIndex">0, 1, 2, or 3; <see cref="Edge"/></param>
        /// <param name="EdgeVertices">input;</param>
        /// <param name="VolumeVertices">output;</param>
        public override void TransformFaceCoordinates(int EdgeIndex, MultidimensionalArray EdgeVertices, MultidimensionalArray VolumeVertices) {
            if (EdgeVertices.Dimension != 2)
                throw new ArgumentOutOfRangeException("dimension of EdgeVertices must be 2.");
            if (VolumeVertices.Dimension != 2)
                throw new ArgumentOutOfRangeException("dimension of VolumeVertices must be 2.");
            if (VolumeVertices.GetLength(1) != SpatialDimension)
                throw new ArgumentOutOfRangeException("wrong spatial dimension of output");
            if (EdgeVertices.GetLength(1) != SpatialDimension - 1)
                throw new ArgumentOutOfRangeException("wrong spatial dimension of input");
            if (EdgeVertices.GetLength(0) != VolumeVertices.GetLength(0))
                throw new ArgumentOutOfRangeException("mismatch in number of vertices between input and output.");

            int L = EdgeVertices.GetLength(0);

            switch (EdgeIndex) {
                case (int)Edge.Left:
                    for (int i = 0; i < L; i++) {
                        VolumeVertices[i, 0] = -1.0;
                        VolumeVertices[i, 1] = EdgeVertices[i, 0];
                    }
                    break;

                case (int)Edge.Right:
                    for (int i = 0; i < L; i++) {
                        VolumeVertices[i, 0] = 1.0;
                        VolumeVertices[i, 1] = EdgeVertices[i, 0];
                    }
                    break;

                case (int)Edge.Top:
                    for (int i = 0; i < L; i++) {
                        VolumeVertices[i, 0] = EdgeVertices[i, 0];
                        VolumeVertices[i, 1] = 1.0;
                    }
                    break;

                case (int)Edge.Bottom:
                    for (int i = 0; i < L; i++) {
                        VolumeVertices[i, 0] = EdgeVertices[i, 0];
                        VolumeVertices[i, 1] = -1.0;
                    }
                    break;

                default:
                    throw new ArgumentException("EdgeIndex out of range");
            }
        }

        /// <summary>
        /// divides the [-1,1]x[-1,1] - square into 4 squares of equal size;
        /// </summary>
        /// <returns></returns>
        public override AffineTrafo[] GetSubdivision() {
            AffineTrafo[] ret = new AffineTrafo[4];

            for (int i = 0; i < 4; i++) {
                ret[i] = new AffineTrafo(2);
                ret[i].Matrix.Clear();
                ret[i].Matrix.AccEye(0.5);
            }

            ret[0].Affine = new double[] { -0.5, -0.5 };
            ret[1].Affine = new double[] { 0.5, -0.5 };
            ret[2].Affine = new double[] { -0.5, 0.5 };
            ret[3].Affine = new double[] { 0.5, 0.5 };

            return ret;
        }


        /// <summary>
        /// tests whether <paramref name="pt"/> is within the convex hull of
        /// vertices or not;
        /// </summary>
        /// <param name="pt"></param>
        /// <param name="tolerance"></param>
        /// <returns></returns>
        public override bool IsWithin(double[] pt, double tolerance) {
            if (pt.Length != 2)
                throw new ArgumentException("wrong spatial dimension.", "pt");

            if ((pt[0] < -1.0 - tolerance) || (pt[0] > 1.0 + tolerance)
                || (pt[1] < -1.0 - tolerance) || (pt[1] > 1.0 + tolerance))
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

            Nodes = MultidimensionalArray.Create(px * px, 2);
            Type = new int[Nodes.GetLength(0)];

            var Nodes1D = GenericBlas.Linspace(-1, 1, px);
            int cnt = 0;
            for (int i = 0; i < px; i++) {
                int i_edge = (i == 0 || i == px - 1) ? 1 : 0;

                for (int j = 0; j < px; j++) {
                    int j_edge = (j == 0 || j == px - 1) ? 1 : 0;
                    Nodes[cnt, 0] = Nodes1D[i];
                    Nodes[cnt, 1] = Nodes1D[j];

                    Type[cnt] = i_edge + j_edge;

                    cnt++;
                }
            }
        }

        /// <summary>
        /// <see cref="RefElement.GetInterpolationNodes_NonLin"/>
        /// </summary>
        override protected void GetInterpolationNodes_NonLin(CellType Type, out NodeSet InterpolationNodes, out PolynomialList InterpolationPolynomials, out int[] NodeType, out int[] EntityIndex) {
            switch (Type) {
                case CellType.Square_4: {
                        base.SelectNodalPolynomials(2, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                        return;
                    }
                case CellType.Square_8: {
                        base.SelectNodalPolynomials(3, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex, ModalBasisSelector: null, NodeTypeFilter: new int[] { 0 });
                        Debug.Assert(NodeType.Length == 8);
                        return;
                    }
                case CellType.Square_9: {
                        base.SelectNodalPolynomials(3, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                        return;
                    }
                case CellType.Square_12: {
                        base.SelectNodalPolynomials(4, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex, NodeTypeFilter: new int[] { 0 },

                            // ill cond.:
                            // ModalBasisSelector:delegate(p) {
                            //    if (p.Degree <= 3)
                            //        return true;
                            //    if (p.Degree == 4 || p.Degree == 6)
                            //        return false;

                            //    if (p.Degree == 5) {
                            //        for (int i = 0; i < p.Coeff.Length; i++) {
                            //            int a1 = p.Exponents[i, 0];
                            //            int a2 = p.Exponents[i, 1];

                            //            if (a1 > 3 || a2 > 3)
                            //                return false;
                            //        }
                            //        return true;
                            //    }

                            //    return false;
                            //}

                             ModalBasisSelector: delegate(Polynomial p) {
                            if (p.AbsoluteDegree <= 3)
                                return true;
                            for (int i = 0; i < p.Coeff.Length; i++) {
                                int a1 = p.Exponents[i, 0];
                                int a2 = p.Exponents[i, 1];

                                if (a1 > 3 || a2 > 3)
                                    return false;


                                //if (a1 == 2 && a2 == 2 && p.Degree == 4)
                                //    return true;
                                //if (a1 == 3 && a2 == 3 && p.Degree == 6)
                                //    return true;

                                if (p.AbsoluteDegree == 4 && ((a1 == 3 && a2 == 1) || (a1 == 1 && a2 == 3)))
                                    return true;
                            }



                            return false;
                        }

                            );
                        Debug.Assert(NodeType.Length == 12);
                        return;
                    }
                case CellType.Square_16: {
                        base.SelectNodalPolynomials(4, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                        return;
                    }
                case CellType.Square_25: {
                        base.SelectNodalPolynomials(5, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                        return;
                    }
                case CellType.Square_36: {
                        base.SelectNodalPolynomials(6, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                        return;
                    }
                case CellType.Square_49: {
                        base.SelectNodalPolynomials(7, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                        return;
                    }
                case CellType.Square_64: {
                        base.SelectNodalPolynomials(8, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                        return;
                    }
                case CellType.Square_81: {
                        base.SelectNodalPolynomials(9, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                        return;
                    }
                case CellType.Square_100: {
                        base.SelectNodalPolynomials(10, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                        return;
                    }
                default:
                    throw new NotImplementedException();
            }
        }

        /// <summary>
        /// <see cref="RefElement.GetForeignElementMapping"/>
        /// </summary>
        /// <param name="Type"></param>
        /// <param name="conv"></param>
        /// <returns></returns>
        public override int[] GetForeignElementMapping(CellType Type, RefElement.ExchangeFormats conv) {
            int[] permutationArray = new int[Vertices.GetLength(0)];
            for (int i = 0; i < permutationArray.Length; i++) {
                permutationArray[i] = i;
            }
            if (conv == ExchangeFormats.Gmsh || conv == ExchangeFormats.CGNS) {
                SwitchNode(ref permutationArray[2], ref permutationArray[3]);
                return permutationArray;
            }
            return permutationArray;
        }

        /// <summary>
        /// <see cref="RefElement.GetForeignElementType"/>
        /// </summary>
        /// <param name="Type"></param>
        /// <param name="conv"></param>
        /// <param name="ForeignName"></param>
        /// <param name="ForeignTypeConstant"></param>
        public override void GetForeignElementType(CellType Type, RefElement.ExchangeFormats conv, out string ForeignName, out int ForeignTypeConstant) {
            ForeignName = "Quadrangle";
            ForeignTypeConstant = 0;
            if (conv == ExchangeFormats.Gmsh) {
                if (Type == CellType.Square_Linear || Type == CellType.Square_Linear) {
                    ForeignTypeConstant = 3;
                } else if (Type == CellType.Square_8) {
                    ForeignTypeConstant = 16;
                } else if (Type == CellType.Square_9) {
                    ForeignTypeConstant = 10;
                } else {
                    throw new NotSupportedException("Wrong minor cell type");
                }
            } else if (conv == ExchangeFormats.CGNS) {
                ForeignName = "Quadrilateral";
                if (Type == CellType.Square_4 || Type == 0) {
                    ForeignTypeConstant = 7;
                } else if (Type == CellType.Square_8) {
                    ForeignTypeConstant = 8;
                } else if (Type == CellType.Square_9) {
                    ForeignTypeConstant = 9;
                } else {
                    throw new NotSupportedException("Wrong minor cell type");
                }
            } else if (conv == ExchangeFormats.GambitNeutral) {
                ForeignName = "Quadrilateral";
                if (Type == CellType.Square_Linear || Type == CellType.Square_Linear) {
                    ForeignTypeConstant = 2;
                } else {
                    throw new NotSupportedException("Wrong minor cell type");
                }
            } else {
                throw new NotSupportedException("Wrong foreign convention type");
            }
        }
    }
}
