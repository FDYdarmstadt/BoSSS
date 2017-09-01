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

using BoSSS.Platform.LinAlg;
using ilPSP.Utils;
using BoSSS.Platform;
using System.Diagnostics;
using BoSSS.Foundation.Quadrature;
using System.Linq;
using System.Collections.Generic;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;

namespace BoSSS.Foundation.Grid.RefElements {

    /// <summary>
    /// The cubic reference element, i.e. \f$ K^{\textrm{cube} } = ( -1,1 )^3 \f$.
    /// Dies ist ein test 
    /// \f[
    /// \testBoSSSmacro
    /// \f]
    /// </summary>
    public class Cube : RefElement {
        

        /// <summary>
        /// The encoding to identify all six faces of the cube.
        /// </summary>
        public enum Edge {
            /// <summary>
            /// edge between this cell the neighbour cell with \f$ x \f$--coordinates closer to negative infinity.
            /// </summary>
            Left = 0,

            /// <summary>
            /// edge between this cell the neighbour cell with \f$ x \f$--coordinates closer to positive infinity.
            /// </summary>
            Right = 1,

            /// <summary>
            /// edge between this cell the neighbour cell with \f$ y \f$--coordinates closer to positive infinity.
            /// </summary>
            Top = 2,

            /// <summary>
            /// edge between this cell the neighbour cell with \f$ y \f$--coordinates closer to negative infinity.
            /// </summary>
            Bottom = 3,


            /// <summary>
            /// edge between this cell the neighbour cell with \f$ z \f$--coordinates closer to positive infinity
            /// </summary>
            Front = 4,

            /// <summary>
            /// edge between this cell the neighbour cell with \f$ z \f$--coordinates closer to negative infinity
            /// </summary>
            Back = 5

        }



        private static Cube instance = null;
        private static readonly object padlock = new object();
        
        /// <summary>
        /// Access to the single, global instance.
        /// </summary>
        public static Cube Instance {
            get {
                lock(padlock) {
                    if(instance == null) {
                        instance = new Cube();
                    }
                    return instance;
                }
            }
        }


        /// <summary>
        /// standard constructor
        /// </summary>
        private Cube() {


            // ===============
            // define vertices
            // ===============

            var _Vertices = new double[8, 3] { { -1, -1, -1 }, { 1, -1, -1 }, { -1, 1, -1 }, { -1, -1, 1 },
                                               { 1, 1, -1 }, { 1, 1, 1 }, { 1, -1, 1 }, { -1, 1, 1 }};
            this.m_Vertices = new NodeSet(this, 8, 3);
            this.m_Vertices.InitializeFrom(_Vertices);
            this.m_Vertices.LockForever();

            m_NoOfFaces = 6;

            // ============
            // edge simplex
            // ============

            m_FaceRefElement = Square.Instance;
            
            // ===================================
            // define Quadrature nodes and weights
            // ===================================

            {
                var qrTemp1D = QuadRuleResource.DecodeFromBase64( Resource.LineQuadRules_bin);
                foreach(var q in qrTemp1D) {
                    int NN = q.Item2.GetLength(0);
                    int D = this.SpatialDimension;
                    var realQr = QuadRule.CreateEmpty(this, NN*NN*NN, D);

                    for(int i = 0; i < NN; i++) {
                        for(int j = 0; j < NN; j++) {
                            for(int k = 0; k < NN; k++) {
                                realQr.Nodes[(i * NN + j) * NN + k, 0] = q.Item2[k, 0];
                                realQr.Nodes[(i * NN + j) * NN + k, 1] = q.Item2[j, 0];
                                realQr.Nodes[(i * NN + j) * NN + k, 2] = q.Item2[i, 0];
                                realQr.Weights[(i * NN + j) * NN + k] = q.Item3[i] * q.Item3[j] * q.Item3[k];
                            }
                        }
                    }

                    realQr.OrderOfPrecision = q.Item1;
                    realQr.Nodes.LockForever();
                    realQr.Weights.LockForever();
                    base.m_QuadRules.Add(realQr);
                }
            }

            // ==============================
            // define orthonormal polynomials
            // ==============================
#pragma warning disable 612
            #region POLY_DEF

            OrthonormalPolynomials = new Polynomial[816];
            Polynomial p;

            p = new Polynomial(new Guid("{FC462EAE-0F6A-49e9-BD17-64E0985F5598}"));
            OrthonormalPolynomials[0] = p;
            p.AddCoeff(3.5355339059327500e-01, new int[] { 0, 0, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{AAAA0C47-F67D-4bed-A0CE-3821B683C485}"));
            OrthonormalPolynomials[1] = p;
            p.AddCoeff(6.1237243569579500e-01, new int[] { 0, 0, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{AC095490-F1B1-43bd-92D0-2C9D3B0EDEDD}"));
            OrthonormalPolynomials[2] = p;
            p.AddCoeff(6.1237243569579500e-01, new int[] { 0, 1, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{6DD46C5B-9160-48b7-8211-F30D5FC9AD9E}"));
            OrthonormalPolynomials[3] = p;
            p.AddCoeff(6.1237243569579500e-01, new int[] { 1, 0, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{C974678B-5B58-406e-B1C1-26EFC968DBA3}"));
            OrthonormalPolynomials[4] = p;
            p.AddCoeff(-3.9528470752104800e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.1858541225631400e+00, new int[] { 0, 0, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{431467E6-BE94-420d-BC85-256817709055}"));
            OrthonormalPolynomials[5] = p;
            p.AddCoeff(1.0606601717798200e+00, new int[] { 0, 1, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{137BFE76-0DDB-46e1-8485-A078F45FD465}"));
            OrthonormalPolynomials[6] = p;
            p.AddCoeff(-3.9528470752104800e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.1858541225631400e+00, new int[] { 0, 2, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{51B9457F-6591-4f8b-85D3-AF72FA1285B2}"));
            OrthonormalPolynomials[7] = p;
            p.AddCoeff(1.0606601717798200e+00, new int[] { 1, 0, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{83F549A3-17E2-43dd-9046-1608089CEC96}"));
            OrthonormalPolynomials[8] = p;
            p.AddCoeff(1.0606601717798200e+00, new int[] { 1, 1, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{B86C4376-3A17-4e9e-8446-62C5F290850B}"));
            OrthonormalPolynomials[9] = p;
            p.AddCoeff(-3.9528470752104800e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.1858541225631400e+00, new int[] { 2, 0, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{3A9FC9AB-36AF-4561-8E86-DB8A28225542}"));
            OrthonormalPolynomials[10] = p;
            p.AddCoeff(-1.4031215200402300e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(2.3385358667337100e+00, new int[] { 0, 0, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{56F6D1AF-17BA-4a8e-94AE-B2B0BFA15DC1}"));
            OrthonormalPolynomials[11] = p;
            p.AddCoeff(-6.8465319688145800e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(2.0539595906443700e+00, new int[] { 0, 1, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{BFB8BF60-CAC2-4839-ADFA-DDF1DF4086DA}"));
            OrthonormalPolynomials[12] = p;
            p.AddCoeff(-6.8465319688145800e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(2.0539595906443700e+00, new int[] { 0, 2, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{38A18AE8-4ED8-4577-B057-43963E1EA372}"));
            OrthonormalPolynomials[13] = p;
            p.AddCoeff(-1.4031215200402300e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(2.3385358667337100e+00, new int[] { 0, 3, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{5CE4F162-0F7D-402e-9C2A-95FAF3ECB5FA}"));
            OrthonormalPolynomials[14] = p;
            p.AddCoeff(-6.8465319688145800e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(2.0539595906443700e+00, new int[] { 1, 0, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{0AB19904-430B-479d-B777-1B810CBC5F8E}"));
            OrthonormalPolynomials[15] = p;
            p.AddCoeff(1.8371173070873800e+00, new int[] { 1, 1, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{21C21D84-9B9D-42ad-B64F-609FEEFA6E52}"));
            OrthonormalPolynomials[16] = p;
            p.AddCoeff(-6.8465319688145800e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(2.0539595906443700e+00, new int[] { 1, 2, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{73080600-ABDA-43b9-AF86-376F4100409B}"));
            OrthonormalPolynomials[17] = p;
            p.AddCoeff(-6.8465319688145800e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(2.0539595906443700e+00, new int[] { 2, 0, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{E23A6B97-BE0A-428e-8691-8A555DEF9D82}"));
            OrthonormalPolynomials[18] = p;
            p.AddCoeff(-6.8465319688145800e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(2.0539595906443700e+00, new int[] { 2, 1, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{D89A37FC-826C-4150-8B02-2747F2E4E1BE}"));
            OrthonormalPolynomials[19] = p;
            p.AddCoeff(-1.4031215200402300e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(2.3385358667337100e+00, new int[] { 3, 0, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{76FB8F1C-D85F-4b71-A783-A28B308C19E2}"));
            OrthonormalPolynomials[20] = p;
            p.AddCoeff(3.9774756441743400e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-3.9774756441743400e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(4.6403882515367300e+00, new int[] { 0, 0, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{9250A297-688F-450f-BCE7-BB036081085D}"));
            OrthonormalPolynomials[21] = p;
            p.AddCoeff(-2.4302777619029500e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(4.0504629365049100e+00, new int[] { 0, 1, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{89CCDCEA-3643-4e6b-8F8D-62DEFC74ED71}"));
            OrthonormalPolynomials[22] = p;
            p.AddCoeff(4.4194173824159400e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-1.3258252147247800e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-1.3258252147247800e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(3.9774756441743400e+00, new int[] { 0, 2, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{DA7D30BD-230C-40cc-9455-28610320E8D5}"));
            OrthonormalPolynomials[23] = p;
            p.AddCoeff(-2.4302777619029500e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(4.0504629365049100e+00, new int[] { 0, 3, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{9B3F6BE2-15D8-4b52-9DD7-E3926288D2ED}"));
            OrthonormalPolynomials[24] = p;
            p.AddCoeff(3.9774756441743400e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-3.9774756441743400e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(4.6403882515367300e+00, new int[] { 0, 4, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{B3E841BC-9722-4f02-9AEC-B7C13F2A2E7D}"));
            OrthonormalPolynomials[25] = p;
            p.AddCoeff(-2.4302777619029500e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(4.0504629365049100e+00, new int[] { 1, 0, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{EC62F440-F476-4019-BD70-6F18A5ABCE35}"));
            OrthonormalPolynomials[26] = p;
            p.AddCoeff(-1.1858541225631400e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(3.5575623676894300e+00, new int[] { 1, 1, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{E6C35E8A-6732-45b4-84D2-BAB5273230D8}"));
            OrthonormalPolynomials[27] = p;
            p.AddCoeff(-1.1858541225631400e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(3.5575623676894300e+00, new int[] { 1, 2, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{7079C244-A823-451c-B95B-425D9879FBA0}"));
            OrthonormalPolynomials[28] = p;
            p.AddCoeff(-2.4302777619029500e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(4.0504629365049100e+00, new int[] { 1, 3, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{0B785F1A-61CA-4ffe-B438-AAA83E02DE5E}"));
            OrthonormalPolynomials[29] = p;
            p.AddCoeff(4.4194173824159400e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-1.3258252147247800e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-1.3258252147247800e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(3.9774756441743400e+00, new int[] { 2, 0, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{0E122787-86EB-4941-B157-6150B19F95A3}"));
            OrthonormalPolynomials[30] = p;
            p.AddCoeff(-1.1858541225631400e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(3.5575623676894300e+00, new int[] { 2, 1, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{F6D58569-33E2-4f40-AF6F-AFA4C6F7F6B2}"));
            OrthonormalPolynomials[31] = p;
            p.AddCoeff(4.4194173824159400e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-1.3258252147247800e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-1.3258252147247800e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(3.9774756441743400e+00, new int[] { 2, 2, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{262E414E-767E-4aed-9684-CEBCB1CA49BD}"));
            OrthonormalPolynomials[32] = p;
            p.AddCoeff(-2.4302777619029500e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(4.0504629365049100e+00, new int[] { 3, 0, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{3535CB47-1B12-4a07-94D2-70572E5E83CA}"));
            OrthonormalPolynomials[33] = p;
            p.AddCoeff(-2.4302777619029500e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(4.0504629365049100e+00, new int[] { 3, 1, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{5608C858-C4B1-4eac-B56A-28B5A65841BA}"));
            OrthonormalPolynomials[34] = p;
            p.AddCoeff(3.9774756441743400e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-3.9774756441743400e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(4.6403882515367300e+00, new int[] { 4, 0, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{FC26F6FB-BD79-4838-B792-37CCA0ABAF15}"));
            OrthonormalPolynomials[35] = p;
            p.AddCoeff(2.1986323874172300e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-1.0260284474613800e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(9.2342560271523800e+00, new int[] { 0, 0, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{678CC9B8-0504-4869-84E1-C3D62B14568B}"));
            OrthonormalPolynomials[36] = p;
            p.AddCoeff(6.8891899015776900e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(-6.8891899015776900e+00, new int[] { 0, 1, 2 });
            p.AddCoeff(8.0373882185073100e+00, new int[] { 0, 1, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{FF65B334-7756-4359-917C-3BB5CB48511E}"));
            OrthonormalPolynomials[37] = p;
            p.AddCoeff(1.5687375497513900e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-2.6145625829189900e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(-4.7062126492541800e+00, new int[] { 0, 2, 1 });
            p.AddCoeff(7.8436877487569600e+00, new int[] { 0, 2, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{E341CCC1-7DD0-46f9-9CC3-11C760E93913}"));
            OrthonormalPolynomials[38] = p;
            p.AddCoeff(1.5687375497513900e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-4.7062126492541800e+00, new int[] { 0, 1, 2 });
            p.AddCoeff(-2.6145625829189900e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(7.8436877487569600e+00, new int[] { 0, 3, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{2E7DE0FF-8F4D-4805-88D1-74BE7A9C4BC0}"));
            OrthonormalPolynomials[39] = p;
            p.AddCoeff(6.8891899015776900e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(-6.8891899015776900e+00, new int[] { 0, 2, 1 });
            p.AddCoeff(8.0373882185073100e+00, new int[] { 0, 4, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{EBB8DFE5-0A5E-43e7-8AAA-C48E68AAB27B}"));
            OrthonormalPolynomials[40] = p;
            p.AddCoeff(2.1986323874172300e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-1.0260284474613800e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(9.2342560271523800e+00, new int[] { 0, 5, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{2A260376-2B3E-431c-9B08-1E6C7589B0A4}"));
            OrthonormalPolynomials[41] = p;
            p.AddCoeff(6.8891899015776900e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(-6.8891899015776900e+00, new int[] { 1, 0, 2 });
            p.AddCoeff(8.0373882185073100e+00, new int[] { 1, 0, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{67B2AA4A-7825-491b-B124-A3AFD3A9CC75}"));
            OrthonormalPolynomials[42] = p;
            p.AddCoeff(-4.2093645601206800e+00, new int[] { 1, 1, 1 });
            p.AddCoeff(7.0156076002011400e+00, new int[] { 1, 1, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{973B0D40-8D9B-4e09-836F-ACC58D58D05D}"));
            OrthonormalPolynomials[43] = p;
            p.AddCoeff(7.6546554461974400e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(-2.2963966338592300e+00, new int[] { 1, 0, 2 });
            p.AddCoeff(-2.2963966338592300e+00, new int[] { 1, 2, 0 });
            p.AddCoeff(6.8891899015776900e+00, new int[] { 1, 2, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{E7B2DC42-B1FC-4885-854C-4D9643B90833}"));
            OrthonormalPolynomials[44] = p;
            p.AddCoeff(-4.2093645601206800e+00, new int[] { 1, 1, 1 });
            p.AddCoeff(7.0156076002011400e+00, new int[] { 1, 3, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{67F276FA-8A25-41c4-A420-6E1A3F774390}"));
            OrthonormalPolynomials[45] = p;
            p.AddCoeff(6.8891899015776900e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(-6.8891899015776900e+00, new int[] { 1, 2, 0 });
            p.AddCoeff(8.0373882185073100e+00, new int[] { 1, 4, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{1E17EA6C-8927-4c9d-98C0-3AC3C178B77E}"));
            OrthonormalPolynomials[46] = p;
            p.AddCoeff(1.5687375497513900e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-2.6145625829189900e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(-4.7062126492541800e+00, new int[] { 2, 0, 1 });
            p.AddCoeff(7.8436877487569600e+00, new int[] { 2, 0, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{58CCA6F6-794B-4fa5-A1C2-4A96C0AAD4C6}"));
            OrthonormalPolynomials[47] = p;
            p.AddCoeff(7.6546554461974400e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(-2.2963966338592300e+00, new int[] { 0, 1, 2 });
            p.AddCoeff(-2.2963966338592300e+00, new int[] { 2, 1, 0 });
            p.AddCoeff(6.8891899015776900e+00, new int[] { 2, 1, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{4EE8B930-BB89-427b-B0A7-2D7D8F0F283E}"));
            OrthonormalPolynomials[48] = p;
            p.AddCoeff(7.6546554461974400e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(-2.2963966338592300e+00, new int[] { 0, 2, 1 });
            p.AddCoeff(-2.2963966338592300e+00, new int[] { 2, 0, 1 });
            p.AddCoeff(6.8891899015776900e+00, new int[] { 2, 2, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{2FE217E4-1668-4592-97CB-71A4D9399672}"));
            OrthonormalPolynomials[49] = p;
            p.AddCoeff(1.5687375497513900e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-2.6145625829189900e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(-4.7062126492541800e+00, new int[] { 2, 1, 0 });
            p.AddCoeff(7.8436877487569600e+00, new int[] { 2, 3, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{36CCFB81-3A03-459a-B12B-97413E6C9A5A}"));
            OrthonormalPolynomials[50] = p;
            p.AddCoeff(1.5687375497513900e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-4.7062126492541800e+00, new int[] { 1, 0, 2 });
            p.AddCoeff(-2.6145625829189900e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(7.8436877487569600e+00, new int[] { 3, 0, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{A33648DD-6ABE-4328-BCD1-15AFF34B3629}"));
            OrthonormalPolynomials[51] = p;
            p.AddCoeff(-4.2093645601206800e+00, new int[] { 1, 1, 1 });
            p.AddCoeff(7.0156076002011400e+00, new int[] { 3, 1, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{7B846DC8-A908-4a30-8011-B1FCF0E294FC}"));
            OrthonormalPolynomials[52] = p;
            p.AddCoeff(1.5687375497513900e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-4.7062126492541800e+00, new int[] { 1, 2, 0 });
            p.AddCoeff(-2.6145625829189900e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(7.8436877487569600e+00, new int[] { 3, 2, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{79CE45FB-7CCD-462d-8C1D-3A0F1F472DE0}"));
            OrthonormalPolynomials[53] = p;
            p.AddCoeff(6.8891899015776900e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(-6.8891899015776900e+00, new int[] { 2, 0, 1 });
            p.AddCoeff(8.0373882185073100e+00, new int[] { 4, 0, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{E75BEE85-2941-4742-9884-BB6C6659538C}"));
            OrthonormalPolynomials[54] = p;
            p.AddCoeff(6.8891899015776900e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(-6.8891899015776900e+00, new int[] { 2, 1, 0 });
            p.AddCoeff(8.0373882185073100e+00, new int[] { 4, 1, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{945C5E09-8304-4f3c-98A6-81BD45A267CB}"));
            OrthonormalPolynomials[55] = p;
            p.AddCoeff(2.1986323874172300e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-1.0260284474613800e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(9.2342560271523800e+00, new int[] { 5, 0, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{86A2EB63-02AE-4e14-B7B3-FA1DC8E7D46F}"));
            OrthonormalPolynomials[56] = p;
            p.AddCoeff(-3.9836089949943600e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(8.3655788894881500e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-2.5096736668464500e+01, new int[] { 0, 0, 4 });
            p.AddCoeff(1.8404273556873900e+01, new int[] { 0, 0, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{DD601938-ADCA-4474-91BE-16B970A8778B}"));
            OrthonormalPolynomials[57] = p;
            p.AddCoeff(3.8081430021731100e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-1.7771334010141200e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(1.5994200609127000e+01, new int[] { 0, 1, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{380B17C3-077D-44cf-BA9C-E2790CDD6317}"));
            OrthonormalPolynomials[58] = p;
            p.AddCoeff(-4.4469529596117800e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(4.4469529596117800e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-5.1881117862137500e+00, new int[] { 0, 0, 4 });
            p.AddCoeff(1.3340858878835400e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-1.3340858878835400e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(1.5564335358641200e+01, new int[] { 0, 2, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{B8A5ABCC-3B77-429b-BCEF-9401F6DD4007}"));
            OrthonormalPolynomials[59] = p;
            p.AddCoeff(5.5684659018440800e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-9.2807765030734700e+00, new int[] { 0, 1, 3 });
            p.AddCoeff(-9.2807765030734700e+00, new int[] { 0, 3, 1 });
            p.AddCoeff(1.5467960838455800e+01, new int[] { 0, 3, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{C0C1CF9A-B0A6-4a6e-9998-64EFBEC4FCFF}"));
            OrthonormalPolynomials[60] = p;
            p.AddCoeff(-4.4469529596117800e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.3340858878835400e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(4.4469529596117800e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-1.3340858878835400e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(-5.1881117862137500e+00, new int[] { 0, 4, 0 });
            p.AddCoeff(1.5564335358641200e+01, new int[] { 0, 4, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{93B5D562-EDC4-4f21-86CF-02781AB61283}"));
            OrthonormalPolynomials[61] = p;
            p.AddCoeff(3.8081430021731100e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-1.7771334010141200e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(1.5994200609127000e+01, new int[] { 0, 5, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{72B82D32-1E0D-4141-9B2D-ABB536A7BBA0}"));
            OrthonormalPolynomials[62] = p;
            p.AddCoeff(-3.9836089949943600e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(8.3655788894881500e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-2.5096736668464500e+01, new int[] { 0, 4, 0 });
            p.AddCoeff(1.8404273556873900e+01, new int[] { 0, 6, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{3C2196C7-655F-4113-9D34-0FD0B9A9171B}"));
            OrthonormalPolynomials[63] = p;
            p.AddCoeff(3.8081430021731100e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-1.7771334010141200e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(1.5994200609127000e+01, new int[] { 1, 0, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{F5226CD9-62A7-46c3-A511-165681414D7E}"));
            OrthonormalPolynomials[64] = p;
            p.AddCoeff(1.1932426932523000e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-1.1932426932523000e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(1.3921164754610200e+01, new int[] { 1, 1, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{AE7946C5-112B-4d6c-97F1-00F896E85A27}"));
            OrthonormalPolynomials[65] = p;
            p.AddCoeff(2.7171331399105100e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-4.5285552331841900e+00, new int[] { 1, 0, 3 });
            p.AddCoeff(-8.1513994197315400e+00, new int[] { 1, 2, 1 });
            p.AddCoeff(1.3585665699552600e+01, new int[] { 1, 2, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{A3E0DE35-2D88-4c0b-BF72-4C7C01799EE2}"));
            OrthonormalPolynomials[66] = p;
            p.AddCoeff(2.7171331399105100e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-8.1513994197315400e+00, new int[] { 1, 1, 2 });
            p.AddCoeff(-4.5285552331841900e+00, new int[] { 1, 3, 0 });
            p.AddCoeff(1.3585665699552600e+01, new int[] { 1, 3, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{E007FA2F-9C4D-4486-8726-5FE30B776B0D}"));
            OrthonormalPolynomials[67] = p;
            p.AddCoeff(1.1932426932523000e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-1.1932426932523000e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(1.3921164754610200e+01, new int[] { 1, 4, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{A1C948CC-44E2-4187-8949-AB2E2D1F4A46}"));
            OrthonormalPolynomials[68] = p;
            p.AddCoeff(3.8081430021731100e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-1.7771334010141200e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(1.5994200609127000e+01, new int[] { 1, 5, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{CA1BDFB1-A111-4d4c-ADD0-4FBDA58ECDA5}"));
            OrthonormalPolynomials[69] = p;
            p.AddCoeff(-4.4469529596117800e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(4.4469529596117800e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-5.1881117862137500e+00, new int[] { 0, 0, 4 });
            p.AddCoeff(1.3340858878835400e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-1.3340858878835400e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(1.5564335358641200e+01, new int[] { 2, 0, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{172482C7-806F-4458-813B-7D03140212C6}"));
            OrthonormalPolynomials[70] = p;
            p.AddCoeff(2.7171331399105100e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-4.5285552331841900e+00, new int[] { 0, 1, 3 });
            p.AddCoeff(-8.1513994197315400e+00, new int[] { 2, 1, 1 });
            p.AddCoeff(1.3585665699552600e+01, new int[] { 2, 1, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{2D688D0E-AF35-4b3e-9B4A-028EA0E4D463}"));
            OrthonormalPolynomials[71] = p;
            p.AddCoeff(-4.9410588440130900e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.4823176532039300e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(1.4823176532039300e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-4.4469529596117800e+00, new int[] { 0, 2, 2 });
            p.AddCoeff(1.4823176532039300e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-4.4469529596117800e+00, new int[] { 2, 0, 2 });
            p.AddCoeff(-4.4469529596117800e+00, new int[] { 2, 2, 0 });
            p.AddCoeff(1.3340858878835400e+01, new int[] { 2, 2, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{24BE4BAE-9BF1-4b86-95F9-C64A71997C81}"));
            OrthonormalPolynomials[72] = p;
            p.AddCoeff(2.7171331399105100e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-4.5285552331841900e+00, new int[] { 0, 3, 1 });
            p.AddCoeff(-8.1513994197315400e+00, new int[] { 2, 1, 1 });
            p.AddCoeff(1.3585665699552600e+01, new int[] { 2, 3, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{DF93E2ED-905D-4fb8-B0DD-DC713D7FBF11}"));
            OrthonormalPolynomials[73] = p;
            p.AddCoeff(-4.4469529596117800e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(4.4469529596117800e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-5.1881117862137500e+00, new int[] { 0, 4, 0 });
            p.AddCoeff(1.3340858878835400e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-1.3340858878835400e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(1.5564335358641200e+01, new int[] { 2, 4, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{29921C04-A942-4cff-8991-FA100ABAD753}"));
            OrthonormalPolynomials[74] = p;
            p.AddCoeff(5.5684659018440800e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-9.2807765030734700e+00, new int[] { 1, 0, 3 });
            p.AddCoeff(-9.2807765030734700e+00, new int[] { 3, 0, 1 });
            p.AddCoeff(1.5467960838455800e+01, new int[] { 3, 0, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{C4FB9ACA-E9EE-4b9b-82ED-9B8241D7DDA1}"));
            OrthonormalPolynomials[75] = p;
            p.AddCoeff(2.7171331399105100e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-8.1513994197315400e+00, new int[] { 1, 1, 2 });
            p.AddCoeff(-4.5285552331841900e+00, new int[] { 3, 1, 0 });
            p.AddCoeff(1.3585665699552600e+01, new int[] { 3, 1, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{F47D934E-C022-46d0-BFF0-B160B86E698C}"));
            OrthonormalPolynomials[76] = p;
            p.AddCoeff(2.7171331399105100e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-8.1513994197315400e+00, new int[] { 1, 2, 1 });
            p.AddCoeff(-4.5285552331841900e+00, new int[] { 3, 0, 1 });
            p.AddCoeff(1.3585665699552600e+01, new int[] { 3, 2, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{F9338F7A-604E-499e-A1A5-BD17321A1D86}"));
            OrthonormalPolynomials[77] = p;
            p.AddCoeff(5.5684659018440800e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-9.2807765030734700e+00, new int[] { 1, 3, 0 });
            p.AddCoeff(-9.2807765030734700e+00, new int[] { 3, 1, 0 });
            p.AddCoeff(1.5467960838455800e+01, new int[] { 3, 3, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{D33CEE11-0518-4314-AC94-56E310E50674}"));
            OrthonormalPolynomials[78] = p;
            p.AddCoeff(-4.4469529596117800e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.3340858878835400e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(4.4469529596117800e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-1.3340858878835400e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(-5.1881117862137500e+00, new int[] { 4, 0, 0 });
            p.AddCoeff(1.5564335358641200e+01, new int[] { 4, 0, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{B2F21F09-D250-48b6-A138-393D5AB202C0}"));
            OrthonormalPolynomials[79] = p;
            p.AddCoeff(1.1932426932523000e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-1.1932426932523000e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(1.3921164754610200e+01, new int[] { 4, 1, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{0E33515F-F3DF-4154-A958-E2D585A31DF4}"));
            OrthonormalPolynomials[80] = p;
            p.AddCoeff(-4.4469529596117800e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.3340858878835400e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(4.4469529596117800e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-1.3340858878835400e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(-5.1881117862137500e+00, new int[] { 4, 0, 0 });
            p.AddCoeff(1.5564335358641200e+01, new int[] { 4, 2, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{90A6A705-93FE-48d3-9686-FDE5A5CBA9B2}"));
            OrthonormalPolynomials[81] = p;
            p.AddCoeff(3.8081430021731100e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-1.7771334010141200e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(1.5994200609127000e+01, new int[] { 5, 0, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{84CA8DC8-D06D-4ec9-B649-80076811824F}"));
            OrthonormalPolynomials[82] = p;
            p.AddCoeff(3.8081430021731100e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-1.7771334010141200e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(1.5994200609127000e+01, new int[] { 5, 1, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("{FCEEBED9-831E-4a91-A959-B78E48703CF8}"));
            OrthonormalPolynomials[83] = p;
            p.AddCoeff(-3.9836089949943600e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(8.3655788894881500e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-2.5096736668464500e+01, new int[] { 4, 0, 0 });
            p.AddCoeff(1.8404273556873900e+01, new int[] { 6, 0, 0 });




            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7dd0c036-0a85-4b7a-b254-0781af1206dc"));
            OrthonormalPolynomials[84] = p;
            p.AddCoeff(-2.9953577363563800e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(2.6958219627207400e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-5.9308083179856300e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(3.6714527682768200e+01, new int[] { 0, 0, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7d57ae11-1477-4d45-bf37-e335eba5b404"));
            OrthonormalPolynomials[85] = p;
            p.AddCoeff(-6.8998131768186300e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(1.4489607671319100e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-4.3468823013957400e+01, new int[] { 0, 1, 4 });
            p.AddCoeff(3.1877136876902100e+01, new int[] { 0, 1, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("db8748fb-fe35-41b4-8485-0b817791b798"));
            OrthonormalPolynomials[86] = p;
            p.AddCoeff(-2.4581457378987900e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(1.1471346776861000e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-1.0324212099174900e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(7.3744372136963700e+00, new int[] { 0, 2, 1 });
            p.AddCoeff(-3.4414040330583000e+01, new int[] { 0, 2, 3 });
            p.AddCoeff(3.0972636297524700e+01, new int[] { 0, 2, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("242c7985-021b-4b31-8b10-935b7f072823"));
            OrthonormalPolynomials[87] = p;
            p.AddCoeff(-1.5785117100452600e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(1.5785117100452600e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-1.8415969950528000e+01, new int[] { 0, 1, 4 });
            p.AddCoeff(2.6308528500754300e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(-2.6308528500754300e+01, new int[] { 0, 3, 2 });
            p.AddCoeff(3.0693283250880000e+01, new int[] { 0, 3, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c3e55d8e-bd78-4e3c-b48d-d960f6285a7a"));
            OrthonormalPolynomials[88] = p;
            p.AddCoeff(-1.5785117100452600e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(2.6308528500754300e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(1.5785117100452600e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-2.6308528500754300e+01, new int[] { 0, 2, 3 });
            p.AddCoeff(-1.8415969950528000e+01, new int[] { 0, 4, 1 });
            p.AddCoeff(3.0693283250880000e+01, new int[] { 0, 4, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3dd41d0e-ea6b-4544-999c-80fdfc1efe16"));
            OrthonormalPolynomials[89] = p;
            p.AddCoeff(-2.4581457378987900e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(7.3744372136963700e+00, new int[] { 0, 1, 2 });
            p.AddCoeff(1.1471346776861000e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-3.4414040330583000e+01, new int[] { 0, 3, 2 });
            p.AddCoeff(-1.0324212099174900e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(3.0972636297524700e+01, new int[] { 0, 5, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("ab55318d-84c9-49a7-8eff-b2dd30276fe7"));
            OrthonormalPolynomials[90] = p;
            p.AddCoeff(-6.8998131768186300e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(1.4489607671319100e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-4.3468823013957400e+01, new int[] { 0, 4, 1 });
            p.AddCoeff(3.1877136876902100e+01, new int[] { 0, 6, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f0f94146-67f0-4a6a-861b-b77363188762"));
            OrthonormalPolynomials[91] = p;
            p.AddCoeff(-2.9953577363563800e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(2.6958219627207400e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-5.9308083179856300e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(3.6714527682768200e+01, new int[] { 0, 7, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("9762a315-a9a8-458d-a02a-daa9a2ca4e10"));
            OrthonormalPolynomials[92] = p;
            p.AddCoeff(-6.8998131768186300e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(1.4489607671319100e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-4.3468823013957400e+01, new int[] { 1, 0, 4 });
            p.AddCoeff(3.1877136876902100e+01, new int[] { 1, 0, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3118a107-b68e-4eeb-95a8-96d5f6fed3e3"));
            OrthonormalPolynomials[93] = p;
            p.AddCoeff(6.5958971622517000e+00, new int[] { 1, 1, 1 });
            p.AddCoeff(-3.0780853423841300e+01, new int[] { 1, 1, 3 });
            p.AddCoeff(2.7702768081457100e+01, new int[] { 1, 1, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("0499ed9a-ecb5-4cab-a9db-fe554f025193"));
            OrthonormalPolynomials[94] = p;
            p.AddCoeff(-7.7023484649164000e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(7.7023484649164000e+00, new int[] { 1, 0, 2 });
            p.AddCoeff(-8.9860732090691300e+00, new int[] { 1, 0, 4 });
            p.AddCoeff(2.3107045394749200e+00, new int[] { 1, 2, 0 });
            p.AddCoeff(-2.3107045394749200e+01, new int[] { 1, 2, 2 });
            p.AddCoeff(2.6958219627207400e+01, new int[] { 1, 2, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("6881b09c-3206-49f5-a416-59fbe35b9c0f"));
            OrthonormalPolynomials[95] = p;
            p.AddCoeff(9.6448658622087700e+00, new int[] { 1, 1, 1 });
            p.AddCoeff(-1.6074776437014600e+01, new int[] { 1, 1, 3 });
            p.AddCoeff(-1.6074776437014600e+01, new int[] { 1, 3, 1 });
            p.AddCoeff(2.6791294061691000e+01, new int[] { 1, 3, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("55057bee-bb5e-48b6-8411-ea32230b761f"));
            OrthonormalPolynomials[96] = p;
            p.AddCoeff(-7.7023484649164000e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(2.3107045394749200e+00, new int[] { 1, 0, 2 });
            p.AddCoeff(7.7023484649164000e+00, new int[] { 1, 2, 0 });
            p.AddCoeff(-2.3107045394749200e+01, new int[] { 1, 2, 2 });
            p.AddCoeff(-8.9860732090691300e+00, new int[] { 1, 4, 0 });
            p.AddCoeff(2.6958219627207400e+01, new int[] { 1, 4, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("61ed9100-f80b-4fb8-8bd3-3223ea1507e3"));
            OrthonormalPolynomials[97] = p;
            p.AddCoeff(6.5958971622517000e+00, new int[] { 1, 1, 1 });
            p.AddCoeff(-3.0780853423841300e+01, new int[] { 1, 3, 1 });
            p.AddCoeff(2.7702768081457100e+01, new int[] { 1, 5, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f911f18f-3ec0-41bb-acb1-63a044297e2c"));
            OrthonormalPolynomials[98] = p;
            p.AddCoeff(-6.8998131768186300e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(1.4489607671319100e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-4.3468823013957400e+01, new int[] { 1, 4, 0 });
            p.AddCoeff(3.1877136876902100e+01, new int[] { 1, 6, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("08b0cf5f-54e2-42a0-a31a-b850cd35c144"));
            OrthonormalPolynomials[99] = p;
            p.AddCoeff(-2.4581457378987900e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(1.1471346776861000e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-1.0324212099174900e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(7.3744372136963700e+00, new int[] { 2, 0, 1 });
            p.AddCoeff(-3.4414040330583000e+01, new int[] { 2, 0, 3 });
            p.AddCoeff(3.0972636297524700e+01, new int[] { 2, 0, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("37fc3a98-ba3f-466b-a61f-4156e9208dcf"));
            OrthonormalPolynomials[100] = p;
            p.AddCoeff(-7.7023484649164000e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(7.7023484649164000e+00, new int[] { 0, 1, 2 });
            p.AddCoeff(-8.9860732090691300e+00, new int[] { 0, 1, 4 });
            p.AddCoeff(2.3107045394749200e+00, new int[] { 2, 1, 0 });
            p.AddCoeff(-2.3107045394749200e+01, new int[] { 2, 1, 2 });
            p.AddCoeff(2.6958219627207400e+01, new int[] { 2, 1, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("ad1386a2-5004-4b35-b07b-bc2e67808c81"));
            OrthonormalPolynomials[101] = p;
            p.AddCoeff(-1.7539019000502800e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(2.9231698334171400e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(5.2617057001508500e+00, new int[] { 0, 2, 1 });
            p.AddCoeff(-8.7695095002514200e+00, new int[] { 0, 2, 3 });
            p.AddCoeff(5.2617057001508500e+00, new int[] { 2, 0, 1 });
            p.AddCoeff(-8.7695095002514200e+00, new int[] { 2, 0, 3 });
            p.AddCoeff(-1.5785117100452600e+01, new int[] { 2, 2, 1 });
            p.AddCoeff(2.6308528500754300e+01, new int[] { 2, 2, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e79173c8-cea5-43d3-ad65-10f315cd20c9"));
            OrthonormalPolynomials[102] = p;
            p.AddCoeff(-1.7539019000502800e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(5.2617057001508500e+00, new int[] { 0, 1, 2 });
            p.AddCoeff(2.9231698334171400e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(-8.7695095002514200e+00, new int[] { 0, 3, 2 });
            p.AddCoeff(5.2617057001508500e+00, new int[] { 2, 1, 0 });
            p.AddCoeff(-1.5785117100452600e+01, new int[] { 2, 1, 2 });
            p.AddCoeff(-8.7695095002514200e+00, new int[] { 2, 3, 0 });
            p.AddCoeff(2.6308528500754300e+01, new int[] { 2, 3, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a9143cea-c951-48e4-8e1c-7c29f36030ff"));
            OrthonormalPolynomials[103] = p;
            p.AddCoeff(-7.7023484649164000e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(7.7023484649164000e+00, new int[] { 0, 2, 1 });
            p.AddCoeff(-8.9860732090691300e+00, new int[] { 0, 4, 1 });
            p.AddCoeff(2.3107045394749200e+00, new int[] { 2, 0, 1 });
            p.AddCoeff(-2.3107045394749200e+01, new int[] { 2, 2, 1 });
            p.AddCoeff(2.6958219627207400e+01, new int[] { 2, 4, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5f83b71c-0ca0-42bb-8c7a-fe794d680787"));
            OrthonormalPolynomials[104] = p;
            p.AddCoeff(-2.4581457378987900e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(1.1471346776861000e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-1.0324212099174900e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(7.3744372136963700e+00, new int[] { 2, 1, 0 });
            p.AddCoeff(-3.4414040330583000e+01, new int[] { 2, 3, 0 });
            p.AddCoeff(3.0972636297524700e+01, new int[] { 2, 5, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("1b94a84e-23c9-4427-a009-e585925e9db5"));
            OrthonormalPolynomials[105] = p;
            p.AddCoeff(-1.5785117100452600e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(1.5785117100452600e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-1.8415969950528000e+01, new int[] { 1, 0, 4 });
            p.AddCoeff(2.6308528500754300e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(-2.6308528500754300e+01, new int[] { 3, 0, 2 });
            p.AddCoeff(3.0693283250880000e+01, new int[] { 3, 0, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("02566e80-16a6-4a41-a8d0-f28618dbb98a"));
            OrthonormalPolynomials[106] = p;
            p.AddCoeff(9.6448658622087700e+00, new int[] { 1, 1, 1 });
            p.AddCoeff(-1.6074776437014600e+01, new int[] { 1, 1, 3 });
            p.AddCoeff(-1.6074776437014600e+01, new int[] { 3, 1, 1 });
            p.AddCoeff(2.6791294061691000e+01, new int[] { 3, 1, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f5f67512-5fae-4735-8c89-4ece63946673"));
            OrthonormalPolynomials[107] = p;
            p.AddCoeff(-1.7539019000502800e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(5.2617057001508500e+00, new int[] { 1, 0, 2 });
            p.AddCoeff(5.2617057001508500e+00, new int[] { 1, 2, 0 });
            p.AddCoeff(-1.5785117100452600e+01, new int[] { 1, 2, 2 });
            p.AddCoeff(2.9231698334171400e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(-8.7695095002514200e+00, new int[] { 3, 0, 2 });
            p.AddCoeff(-8.7695095002514200e+00, new int[] { 3, 2, 0 });
            p.AddCoeff(2.6308528500754300e+01, new int[] { 3, 2, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("cdce2d4a-8d4c-4f0c-ae05-2f75eaf59de3"));
            OrthonormalPolynomials[108] = p;
            p.AddCoeff(9.6448658622087700e+00, new int[] { 1, 1, 1 });
            p.AddCoeff(-1.6074776437014600e+01, new int[] { 1, 3, 1 });
            p.AddCoeff(-1.6074776437014600e+01, new int[] { 3, 1, 1 });
            p.AddCoeff(2.6791294061691000e+01, new int[] { 3, 3, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("da1ba734-ec83-40cf-9a7d-eaaa8c18f09d"));
            OrthonormalPolynomials[109] = p;
            p.AddCoeff(-1.5785117100452600e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(1.5785117100452600e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-1.8415969950528000e+01, new int[] { 1, 4, 0 });
            p.AddCoeff(2.6308528500754300e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(-2.6308528500754300e+01, new int[] { 3, 2, 0 });
            p.AddCoeff(3.0693283250880000e+01, new int[] { 3, 4, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("2e3be0f6-2c25-4595-b4c3-7239a86c5607"));
            OrthonormalPolynomials[110] = p;
            p.AddCoeff(-1.5785117100452600e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(2.6308528500754300e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(1.5785117100452600e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-2.6308528500754300e+01, new int[] { 2, 0, 3 });
            p.AddCoeff(-1.8415969950528000e+01, new int[] { 4, 0, 1 });
            p.AddCoeff(3.0693283250880000e+01, new int[] { 4, 0, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("775d6717-d53c-41a0-8bf8-ed8ecd8f9721"));
            OrthonormalPolynomials[111] = p;
            p.AddCoeff(-7.7023484649164000e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(2.3107045394749200e+00, new int[] { 0, 1, 2 });
            p.AddCoeff(7.7023484649164000e+00, new int[] { 2, 1, 0 });
            p.AddCoeff(-2.3107045394749200e+01, new int[] { 2, 1, 2 });
            p.AddCoeff(-8.9860732090691300e+00, new int[] { 4, 1, 0 });
            p.AddCoeff(2.6958219627207400e+01, new int[] { 4, 1, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5a3b94a6-d491-4973-9839-f7b7c252c2aa"));
            OrthonormalPolynomials[112] = p;
            p.AddCoeff(-7.7023484649164000e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(2.3107045394749200e+00, new int[] { 0, 2, 1 });
            p.AddCoeff(7.7023484649164000e+00, new int[] { 2, 0, 1 });
            p.AddCoeff(-2.3107045394749200e+01, new int[] { 2, 2, 1 });
            p.AddCoeff(-8.9860732090691300e+00, new int[] { 4, 0, 1 });
            p.AddCoeff(2.6958219627207400e+01, new int[] { 4, 2, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("94031ced-3f18-465e-be05-f5c1ad555618"));
            OrthonormalPolynomials[113] = p;
            p.AddCoeff(-1.5785117100452600e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(2.6308528500754300e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(1.5785117100452600e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-2.6308528500754300e+01, new int[] { 2, 3, 0 });
            p.AddCoeff(-1.8415969950528000e+01, new int[] { 4, 1, 0 });
            p.AddCoeff(3.0693283250880000e+01, new int[] { 4, 3, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3e7006e4-c207-4108-b6b9-15e4187b498f"));
            OrthonormalPolynomials[114] = p;
            p.AddCoeff(-2.4581457378987900e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(7.3744372136963700e+00, new int[] { 1, 0, 2 });
            p.AddCoeff(1.1471346776861000e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-3.4414040330583000e+01, new int[] { 3, 0, 2 });
            p.AddCoeff(-1.0324212099174900e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(3.0972636297524700e+01, new int[] { 5, 0, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("6f0d8d5f-d5ac-42b4-866c-1c41b496cb97"));
            OrthonormalPolynomials[115] = p;
            p.AddCoeff(6.5958971622517000e+00, new int[] { 1, 1, 1 });
            p.AddCoeff(-3.0780853423841300e+01, new int[] { 3, 1, 1 });
            p.AddCoeff(2.7702768081457100e+01, new int[] { 5, 1, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("1dfb234e-8a40-46f5-ac1d-5df052a4ba09"));
            OrthonormalPolynomials[116] = p;
            p.AddCoeff(-2.4581457378987900e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(7.3744372136963700e+00, new int[] { 1, 2, 0 });
            p.AddCoeff(1.1471346776861000e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-3.4414040330583000e+01, new int[] { 3, 2, 0 });
            p.AddCoeff(-1.0324212099174900e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(3.0972636297524700e+01, new int[] { 5, 2, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("138c1ac4-4bc1-4273-b169-a1db870bd588"));
            OrthonormalPolynomials[117] = p;
            p.AddCoeff(-6.8998131768186300e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(1.4489607671319100e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-4.3468823013957400e+01, new int[] { 4, 0, 1 });
            p.AddCoeff(3.1877136876902100e+01, new int[] { 6, 0, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b538b1e5-6c4c-48eb-a425-08dd61424710"));
            OrthonormalPolynomials[118] = p;
            p.AddCoeff(-6.8998131768186300e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(1.4489607671319100e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-4.3468823013957400e+01, new int[] { 4, 1, 0 });
            p.AddCoeff(3.1877136876902100e+01, new int[] { 6, 1, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b6a2c91a-fb89-4086-9075-41c0d18fac1c"));
            OrthonormalPolynomials[119] = p;
            p.AddCoeff(-2.9953577363563800e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(2.6958219627207400e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-5.9308083179856300e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(3.6714527682768200e+01, new int[] { 7, 0, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e9a02821-6558-499b-9cfa-97c5d6efec6b"));
            OrthonormalPolynomials[120] = p;
            p.AddCoeff(3.9860022718669000e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-1.4349608178720900e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(7.8922844982964700e+01, new int[] { 0, 0, 4 });
            p.AddCoeff(-1.3679959797047200e+02, new int[] { 0, 0, 6 });
            p.AddCoeff(7.3285498912752900e+01, new int[] { 0, 0, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("16bf8ed9-0a33-4074-81d0-27ab9c9b789a"));
            OrthonormalPolynomials[121] = p;
            p.AddCoeff(-5.1881117862137500e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(4.6693006075923700e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(-1.0272461336703200e+02, new int[] { 0, 1, 5 });
            p.AddCoeff(6.3591427322448500e+01, new int[] { 0, 1, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a14cc51b-2541-463d-9f43-6123d5c64020"));
            OrthonormalPolynomials[122] = p;
            p.AddCoeff(4.4538102542935200e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-9.3530015340163800e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(2.8059004602049100e+01, new int[] { 0, 0, 4 });
            p.AddCoeff(-2.0576603374836000e+01, new int[] { 0, 0, 6 });
            p.AddCoeff(-1.3361430762880500e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(2.8059004602049100e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(-8.4177013806147400e+01, new int[] { 0, 2, 4 });
            p.AddCoeff(6.1729810124508100e+01, new int[] { 0, 2, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("9f32ec9a-8126-49b2-ba46-5fdb78d14143"));
            OrthonormalPolynomials[123] = p;
            p.AddCoeff(-8.7255517823373500e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(4.0719241650907600e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(-3.6647317485816900e+01, new int[] { 0, 1, 5 });
            p.AddCoeff(1.4542586303895600e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(-6.7865402751512700e+01, new int[] { 0, 3, 3 });
            p.AddCoeff(6.1078862476361500e+01, new int[] { 0, 3, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f85af4e7-ad13-4200-89ee-706ce64ea3d9"));
            OrthonormalPolynomials[124] = p;
            p.AddCoeff(4.4746600996961400e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-4.4746600996961400e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(5.2204367829788300e+00, new int[] { 0, 0, 4 });
            p.AddCoeff(-4.4746600996961400e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(4.4746600996961400e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(-5.2204367829788300e+01, new int[] { 0, 2, 4 });
            p.AddCoeff(5.2204367829788300e+00, new int[] { 0, 4, 0 });
            p.AddCoeff(-5.2204367829788300e+01, new int[] { 0, 4, 2 });
            p.AddCoeff(6.0905095801419600e+01, new int[] { 0, 4, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c7ac3db3-d802-4a7f-9409-2717524f18c0"));
            OrthonormalPolynomials[125] = p;
            p.AddCoeff(-8.7255517823373500e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(1.4542586303895600e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(4.0719241650907600e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(-6.7865402751512700e+01, new int[] { 0, 3, 3 });
            p.AddCoeff(-3.6647317485816900e+01, new int[] { 0, 5, 1 });
            p.AddCoeff(6.1078862476361500e+01, new int[] { 0, 5, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("74c747b2-0b63-4d25-92d0-46b4be8cafb8"));
            OrthonormalPolynomials[126] = p;
            p.AddCoeff(4.4538102542935200e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-1.3361430762880500e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-9.3530015340163800e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(2.8059004602049100e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(2.8059004602049100e+01, new int[] { 0, 4, 0 });
            p.AddCoeff(-8.4177013806147400e+01, new int[] { 0, 4, 2 });
            p.AddCoeff(-2.0576603374836000e+01, new int[] { 0, 6, 0 });
            p.AddCoeff(6.1729810124508100e+01, new int[] { 0, 6, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c05fc5be-ea5d-4e37-a3d8-e88ecdbdab56"));
            OrthonormalPolynomials[127] = p;
            p.AddCoeff(-5.1881117862137500e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(4.6693006075923700e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(-1.0272461336703200e+02, new int[] { 0, 5, 1 });
            p.AddCoeff(6.3591427322448500e+01, new int[] { 0, 7, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("94b271bd-e9f3-4f4f-9366-f8b41b838d9f"));
            OrthonormalPolynomials[128] = p;
            p.AddCoeff(3.9860022718669000e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-1.4349608178720900e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(7.8922844982964700e+01, new int[] { 0, 4, 0 });
            p.AddCoeff(-1.3679959797047200e+02, new int[] { 0, 6, 0 });
            p.AddCoeff(7.3285498912752900e+01, new int[] { 0, 8, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("98719eeb-8d79-46a0-aac8-4fbad9239186"));
            OrthonormalPolynomials[129] = p;
            p.AddCoeff(-5.1881117862137500e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(4.6693006075923700e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(-1.0272461336703200e+02, new int[] { 1, 0, 5 });
            p.AddCoeff(6.3591427322448500e+01, new int[] { 1, 0, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c30464c6-6f00-464e-a5e9-ae805dd048ac"));
            OrthonormalPolynomials[130] = p;
            p.AddCoeff(-1.1950826984983100e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(2.5096736668464500e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(-7.5290210005393400e+01, new int[] { 1, 1, 4 });
            p.AddCoeff(5.5212820670621800e+01, new int[] { 1, 1, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("2275e567-8481-43ee-8d69-8d6b95cd1603"));
            OrthonormalPolynomials[131] = p;
            p.AddCoeff(-4.2576333104495900e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(1.9868955448764700e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(-1.7882059903888300e+01, new int[] { 1, 0, 5 });
            p.AddCoeff(1.2772899931348800e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(-5.9606866346294200e+01, new int[] { 1, 2, 3 });
            p.AddCoeff(5.3646179711664800e+01, new int[] { 1, 2, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a1fe14cd-6463-4e45-a27f-225c1ac90024"));
            OrthonormalPolynomials[132] = p;
            p.AddCoeff(-2.7340624821408200e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(2.7340624821408200e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(-3.1897395624976200e+01, new int[] { 1, 1, 4 });
            p.AddCoeff(4.5567708035680300e+00, new int[] { 1, 3, 0 });
            p.AddCoeff(-4.5567708035680300e+01, new int[] { 1, 3, 2 });
            p.AddCoeff(5.3162326041627000e+01, new int[] { 1, 3, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("95edd4a1-4cb3-4ff5-93c5-d3e71595e874"));
            OrthonormalPolynomials[133] = p;
            p.AddCoeff(-2.7340624821408200e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(4.5567708035680300e+00, new int[] { 1, 0, 3 });
            p.AddCoeff(2.7340624821408200e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(-4.5567708035680300e+01, new int[] { 1, 2, 3 });
            p.AddCoeff(-3.1897395624976200e+01, new int[] { 1, 4, 1 });
            p.AddCoeff(5.3162326041627000e+01, new int[] { 1, 4, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("0599bcb0-dc38-446e-aa6f-a498e89ec336"));
            OrthonormalPolynomials[134] = p;
            p.AddCoeff(-4.2576333104495900e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(1.2772899931348800e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(1.9868955448764700e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(-5.9606866346294200e+01, new int[] { 1, 3, 2 });
            p.AddCoeff(-1.7882059903888300e+01, new int[] { 1, 5, 0 });
            p.AddCoeff(5.3646179711664800e+01, new int[] { 1, 5, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5951b721-ce34-4038-a8ba-0849323315c2"));
            OrthonormalPolynomials[135] = p;
            p.AddCoeff(-1.1950826984983100e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(2.5096736668464500e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(-7.5290210005393400e+01, new int[] { 1, 4, 1 });
            p.AddCoeff(5.5212820670621800e+01, new int[] { 1, 6, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("86e5e20a-090b-4d74-b5ca-4232fe8782f1"));
            OrthonormalPolynomials[136] = p;
            p.AddCoeff(-5.1881117862137500e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(4.6693006075923700e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(-1.0272461336703200e+02, new int[] { 1, 5, 0 });
            p.AddCoeff(6.3591427322448500e+01, new int[] { 1, 7, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("4e56d166-0c01-4f40-ab44-97cf3bf0215f"));
            OrthonormalPolynomials[137] = p;
            p.AddCoeff(4.4538102542935200e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-9.3530015340163800e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(2.8059004602049100e+01, new int[] { 0, 0, 4 });
            p.AddCoeff(-2.0576603374836000e+01, new int[] { 0, 0, 6 });
            p.AddCoeff(-1.3361430762880500e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(2.8059004602049100e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(-8.4177013806147400e+01, new int[] { 2, 0, 4 });
            p.AddCoeff(6.1729810124508100e+01, new int[] { 2, 0, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("ff3ddf81-e66b-47e2-8068-752ddcd91346"));
            OrthonormalPolynomials[138] = p;
            p.AddCoeff(-4.2576333104495900e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(1.9868955448764700e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(-1.7882059903888300e+01, new int[] { 0, 1, 5 });
            p.AddCoeff(1.2772899931348800e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(-5.9606866346294200e+01, new int[] { 2, 1, 3 });
            p.AddCoeff(5.3646179711664800e+01, new int[] { 2, 1, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("274fb628-7f50-45be-bc96-1b5f444671e7"));
            OrthonormalPolynomials[139] = p;
            p.AddCoeff(4.9718445552179300e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-4.9718445552179300e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(5.8004853144209200e+00, new int[] { 0, 0, 4 });
            p.AddCoeff(-1.4915533665653800e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(1.4915533665653800e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(-1.7401455943262800e+01, new int[] { 0, 2, 4 });
            p.AddCoeff(-1.4915533665653800e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(1.4915533665653800e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(-1.7401455943262800e+01, new int[] { 2, 0, 4 });
            p.AddCoeff(4.4746600996961400e+00, new int[] { 2, 2, 0 });
            p.AddCoeff(-4.4746600996961400e+01, new int[] { 2, 2, 2 });
            p.AddCoeff(5.2204367829788300e+01, new int[] { 2, 2, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("00d910f8-2c40-43c6-9b9b-739f578cb41a"));
            OrthonormalPolynomials[140] = p;
            p.AddCoeff(-6.2257341434565000e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(1.0376223572427500e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(1.0376223572427500e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(-1.7293705954045800e+01, new int[] { 0, 3, 3 });
            p.AddCoeff(1.8677202430369500e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(-3.1128670717282500e+01, new int[] { 2, 1, 3 });
            p.AddCoeff(-3.1128670717282500e+01, new int[] { 2, 3, 1 });
            p.AddCoeff(5.1881117862137500e+01, new int[] { 2, 3, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("0d0761ee-e566-4ea7-9a3e-ca060632f154"));
            OrthonormalPolynomials[141] = p;
            p.AddCoeff(4.9718445552179300e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-1.4915533665653800e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-4.9718445552179300e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(1.4915533665653800e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(5.8004853144209200e+00, new int[] { 0, 4, 0 });
            p.AddCoeff(-1.7401455943262800e+01, new int[] { 0, 4, 2 });
            p.AddCoeff(-1.4915533665653800e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(4.4746600996961400e+00, new int[] { 2, 0, 2 });
            p.AddCoeff(1.4915533665653800e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(-4.4746600996961400e+01, new int[] { 2, 2, 2 });
            p.AddCoeff(-1.7401455943262800e+01, new int[] { 2, 4, 0 });
            p.AddCoeff(5.2204367829788300e+01, new int[] { 2, 4, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a4c34ebd-a7a2-4dfd-930f-48d9f588ba51"));
            OrthonormalPolynomials[142] = p;
            p.AddCoeff(-4.2576333104495900e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(1.9868955448764700e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(-1.7882059903888300e+01, new int[] { 0, 5, 1 });
            p.AddCoeff(1.2772899931348800e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(-5.9606866346294200e+01, new int[] { 2, 3, 1 });
            p.AddCoeff(5.3646179711664800e+01, new int[] { 2, 5, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("bd4ed659-09c7-4016-9bbb-ad76e9c3818b"));
            OrthonormalPolynomials[143] = p;
            p.AddCoeff(4.4538102542935200e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-9.3530015340163800e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(2.8059004602049100e+01, new int[] { 0, 4, 0 });
            p.AddCoeff(-2.0576603374836000e+01, new int[] { 0, 6, 0 });
            p.AddCoeff(-1.3361430762880500e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(2.8059004602049100e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(-8.4177013806147400e+01, new int[] { 2, 4, 0 });
            p.AddCoeff(6.1729810124508100e+01, new int[] { 2, 6, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("2d304046-b7b4-475a-b42a-ae3212efa284"));
            OrthonormalPolynomials[144] = p;
            p.AddCoeff(-8.7255517823373500e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(4.0719241650907600e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(-3.6647317485816900e+01, new int[] { 1, 0, 5 });
            p.AddCoeff(1.4542586303895600e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(-6.7865402751512700e+01, new int[] { 3, 0, 3 });
            p.AddCoeff(6.1078862476361500e+01, new int[] { 3, 0, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("72d816a5-3af6-4424-94fc-00fdb6202267"));
            OrthonormalPolynomials[145] = p;
            p.AddCoeff(-2.7340624821408200e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(2.7340624821408200e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(-3.1897395624976200e+01, new int[] { 1, 1, 4 });
            p.AddCoeff(4.5567708035680300e+00, new int[] { 3, 1, 0 });
            p.AddCoeff(-4.5567708035680300e+01, new int[] { 3, 1, 2 });
            p.AddCoeff(5.3162326041627000e+01, new int[] { 3, 1, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("8bcc338d-093a-4c7f-8aa9-c5a3590e18c8"));
            OrthonormalPolynomials[146] = p;
            p.AddCoeff(-6.2257341434565000e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(1.0376223572427500e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(1.8677202430369500e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(-3.1128670717282500e+01, new int[] { 1, 2, 3 });
            p.AddCoeff(1.0376223572427500e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(-1.7293705954045800e+01, new int[] { 3, 0, 3 });
            p.AddCoeff(-3.1128670717282500e+01, new int[] { 3, 2, 1 });
            p.AddCoeff(5.1881117862137500e+01, new int[] { 3, 2, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("8a610160-3094-45a4-8b76-4492583e289a"));
            OrthonormalPolynomials[147] = p;
            p.AddCoeff(-6.2257341434565000e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(1.8677202430369500e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(1.0376223572427500e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(-3.1128670717282500e+01, new int[] { 1, 3, 2 });
            p.AddCoeff(1.0376223572427500e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(-3.1128670717282500e+01, new int[] { 3, 1, 2 });
            p.AddCoeff(-1.7293705954045800e+01, new int[] { 3, 3, 0 });
            p.AddCoeff(5.1881117862137500e+01, new int[] { 3, 3, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("131c4a91-b622-4b70-9a81-5be1d6c9bdbb"));
            OrthonormalPolynomials[148] = p;
            p.AddCoeff(-2.7340624821408200e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(2.7340624821408200e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(-3.1897395624976200e+01, new int[] { 1, 4, 1 });
            p.AddCoeff(4.5567708035680300e+00, new int[] { 3, 0, 1 });
            p.AddCoeff(-4.5567708035680300e+01, new int[] { 3, 2, 1 });
            p.AddCoeff(5.3162326041627000e+01, new int[] { 3, 4, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("978748e8-9621-4976-a30d-df6bf8349328"));
            OrthonormalPolynomials[149] = p;
            p.AddCoeff(-8.7255517823373500e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(4.0719241650907600e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(-3.6647317485816900e+01, new int[] { 1, 5, 0 });
            p.AddCoeff(1.4542586303895600e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(-6.7865402751512700e+01, new int[] { 3, 3, 0 });
            p.AddCoeff(6.1078862476361500e+01, new int[] { 3, 5, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d4a37d1e-c6b0-4c63-90f1-3c964db5ab26"));
            OrthonormalPolynomials[150] = p;
            p.AddCoeff(4.4746600996961400e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-4.4746600996961400e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(5.2204367829788300e+00, new int[] { 0, 0, 4 });
            p.AddCoeff(-4.4746600996961400e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(4.4746600996961400e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(-5.2204367829788300e+01, new int[] { 2, 0, 4 });
            p.AddCoeff(5.2204367829788300e+00, new int[] { 4, 0, 0 });
            p.AddCoeff(-5.2204367829788300e+01, new int[] { 4, 0, 2 });
            p.AddCoeff(6.0905095801419600e+01, new int[] { 4, 0, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7981b960-5752-4afc-9b92-085d7360b36e"));
            OrthonormalPolynomials[151] = p;
            p.AddCoeff(-2.7340624821408200e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(4.5567708035680300e+00, new int[] { 0, 1, 3 });
            p.AddCoeff(2.7340624821408200e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(-4.5567708035680300e+01, new int[] { 2, 1, 3 });
            p.AddCoeff(-3.1897395624976200e+01, new int[] { 4, 1, 1 });
            p.AddCoeff(5.3162326041627000e+01, new int[] { 4, 1, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e14eac30-b8f9-4816-8e8a-590a419abb3d"));
            OrthonormalPolynomials[152] = p;
            p.AddCoeff(4.9718445552179300e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-1.4915533665653800e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-1.4915533665653800e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(4.4746600996961400e+00, new int[] { 0, 2, 2 });
            p.AddCoeff(-4.9718445552179300e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(1.4915533665653800e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(1.4915533665653800e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(-4.4746600996961400e+01, new int[] { 2, 2, 2 });
            p.AddCoeff(5.8004853144209200e+00, new int[] { 4, 0, 0 });
            p.AddCoeff(-1.7401455943262800e+01, new int[] { 4, 0, 2 });
            p.AddCoeff(-1.7401455943262800e+01, new int[] { 4, 2, 0 });
            p.AddCoeff(5.2204367829788300e+01, new int[] { 4, 2, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f025d439-b469-48ca-83a9-16772b8dd6d3"));
            OrthonormalPolynomials[153] = p;
            p.AddCoeff(-2.7340624821408200e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(4.5567708035680300e+00, new int[] { 0, 3, 1 });
            p.AddCoeff(2.7340624821408200e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(-4.5567708035680300e+01, new int[] { 2, 3, 1 });
            p.AddCoeff(-3.1897395624976200e+01, new int[] { 4, 1, 1 });
            p.AddCoeff(5.3162326041627000e+01, new int[] { 4, 3, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("4aba7006-9df9-4943-a78a-78805a48223b"));
            OrthonormalPolynomials[154] = p;
            p.AddCoeff(4.4746600996961400e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-4.4746600996961400e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(5.2204367829788300e+00, new int[] { 0, 4, 0 });
            p.AddCoeff(-4.4746600996961400e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(4.4746600996961400e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(-5.2204367829788300e+01, new int[] { 2, 4, 0 });
            p.AddCoeff(5.2204367829788300e+00, new int[] { 4, 0, 0 });
            p.AddCoeff(-5.2204367829788300e+01, new int[] { 4, 2, 0 });
            p.AddCoeff(6.0905095801419600e+01, new int[] { 4, 4, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("ee863ab9-9d05-4e6e-bf50-fb5939f22e07"));
            OrthonormalPolynomials[155] = p;
            p.AddCoeff(-8.7255517823373500e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(1.4542586303895600e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(4.0719241650907600e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(-6.7865402751512700e+01, new int[] { 3, 0, 3 });
            p.AddCoeff(-3.6647317485816900e+01, new int[] { 5, 0, 1 });
            p.AddCoeff(6.1078862476361500e+01, new int[] { 5, 0, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7d8d8635-f1da-4bd0-b52e-49e104a8a74c"));
            OrthonormalPolynomials[156] = p;
            p.AddCoeff(-4.2576333104495900e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(1.2772899931348800e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(1.9868955448764700e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(-5.9606866346294200e+01, new int[] { 3, 1, 2 });
            p.AddCoeff(-1.7882059903888300e+01, new int[] { 5, 1, 0 });
            p.AddCoeff(5.3646179711664800e+01, new int[] { 5, 1, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e50cd454-f6b2-4fcc-9847-11491648d18b"));
            OrthonormalPolynomials[157] = p;
            p.AddCoeff(-4.2576333104495900e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(1.2772899931348800e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(1.9868955448764700e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(-5.9606866346294200e+01, new int[] { 3, 2, 1 });
            p.AddCoeff(-1.7882059903888300e+01, new int[] { 5, 0, 1 });
            p.AddCoeff(5.3646179711664800e+01, new int[] { 5, 2, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("546bdc18-e473-4b6d-8b33-88c309baffa9"));
            OrthonormalPolynomials[158] = p;
            p.AddCoeff(-8.7255517823373500e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(1.4542586303895600e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(4.0719241650907600e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(-6.7865402751512700e+01, new int[] { 3, 3, 0 });
            p.AddCoeff(-3.6647317485816900e+01, new int[] { 5, 1, 0 });
            p.AddCoeff(6.1078862476361500e+01, new int[] { 5, 3, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("0ea377b6-6e1a-42ad-9a78-3c1e968b32db"));
            OrthonormalPolynomials[159] = p;
            p.AddCoeff(4.4538102542935200e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-1.3361430762880500e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-9.3530015340163800e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(2.8059004602049100e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(2.8059004602049100e+01, new int[] { 4, 0, 0 });
            p.AddCoeff(-8.4177013806147400e+01, new int[] { 4, 0, 2 });
            p.AddCoeff(-2.0576603374836000e+01, new int[] { 6, 0, 0 });
            p.AddCoeff(6.1729810124508100e+01, new int[] { 6, 0, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("8aeeac54-3891-482d-be68-f00412d57d8f"));
            OrthonormalPolynomials[160] = p;
            p.AddCoeff(-1.1950826984983100e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(2.5096736668464500e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(-7.5290210005393400e+01, new int[] { 4, 1, 1 });
            p.AddCoeff(5.5212820670621800e+01, new int[] { 6, 1, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("9d8090ea-781d-431d-955b-8fb6c58bb8c5"));
            OrthonormalPolynomials[161] = p;
            p.AddCoeff(4.4538102542935200e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-1.3361430762880500e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-9.3530015340163800e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(2.8059004602049100e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(2.8059004602049100e+01, new int[] { 4, 0, 0 });
            p.AddCoeff(-8.4177013806147400e+01, new int[] { 4, 2, 0 });
            p.AddCoeff(-2.0576603374836000e+01, new int[] { 6, 0, 0 });
            p.AddCoeff(6.1729810124508100e+01, new int[] { 6, 2, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3bfad4b6-0b8d-4c80-a07e-433b0b35593b"));
            OrthonormalPolynomials[162] = p;
            p.AddCoeff(-5.1881117862137500e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(4.6693006075923700e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(-1.0272461336703200e+02, new int[] { 5, 0, 1 });
            p.AddCoeff(6.3591427322448500e+01, new int[] { 7, 0, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5088572c-a2c3-40df-9934-196db332a5e4"));
            OrthonormalPolynomials[163] = p;
            p.AddCoeff(-5.1881117862137500e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(4.6693006075923700e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(-1.0272461336703200e+02, new int[] { 5, 1, 0 });
            p.AddCoeff(6.3591427322448500e+01, new int[] { 7, 1, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("13508120-cf28-48fe-bb52-01a2128f39ae"));
            OrthonormalPolynomials[164] = p;
            p.AddCoeff(3.9860022718669000e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-1.4349608178720900e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(7.8922844982964700e+01, new int[] { 4, 0, 0 });
            p.AddCoeff(-1.3679959797047200e+02, new int[] { 6, 0, 0 });
            p.AddCoeff(7.3285498912752900e+01, new int[] { 8, 0, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3c4ade16-d78b-4bb4-ac51-e92bbcdbea00"));
            OrthonormalPolynomials[165] = p;
            p.AddCoeff(3.7925593963578700e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-5.5624204479915400e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(2.1693439747167000e+02, new int[] { 0, 0, 5 });
            p.AddCoeff(-3.0990628210238600e+02, new int[] { 0, 0, 7 });
            p.AddCoeff(1.4634463321501600e+02, new int[] { 0, 0, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("bb65a625-f82b-40c7-9a1c-01296ef1196a"));
            OrthonormalPolynomials[166] = p;
            p.AddCoeff(6.9039584539584700e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(-2.4854250434250500e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(1.3669837738837800e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(-2.3694385413985500e+02, new int[] { 0, 1, 6 });
            p.AddCoeff(1.2693420757492200e+02, new int[] { 0, 1, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("6159b2dc-d36e-409f-bb87-d90b4515039b"));
            OrthonormalPolynomials[167] = p;
            p.AddCoeff(3.3489117577113800e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-3.0140205819402400e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(6.6308452802685300e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(-4.1048089830233800e+01, new int[] { 0, 0, 7 });
            p.AddCoeff(-1.0046735273134100e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(9.0420617458207200e+01, new int[] { 0, 2, 3 });
            p.AddCoeff(-1.9892535840805600e+02, new int[] { 0, 2, 5 });
            p.AddCoeff(1.2314426949070100e+02, new int[] { 0, 2, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e657ce88-a058-499d-a884-e3ca2397d9ce"));
            OrthonormalPolynomials[168] = p;
            p.AddCoeff(1.5809458081912500e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-3.3199861972016200e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(9.9599585916048800e+01, new int[] { 0, 1, 4 });
            p.AddCoeff(-7.3039696338435800e+01, new int[] { 0, 1, 6 });
            p.AddCoeff(-2.6349096803187500e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(5.5333103286693800e+01, new int[] { 0, 3, 2 });
            p.AddCoeff(-1.6599930986008100e+02, new int[] { 0, 3, 4 });
            p.AddCoeff(1.2173282723072600e+02, new int[] { 0, 3, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d9b12d18-0eb3-46a5-86fc-16e0bbecdc62"));
            OrthonormalPolynomials[169] = p;
            p.AddCoeff(2.4734614358443900e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-1.1542820033940500e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(1.0388538030546400e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(-2.4734614358443900e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(1.1542820033940500e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(-1.0388538030546400e+02, new int[] { 0, 2, 5 });
            p.AddCoeff(2.8857050084851200e+01, new int[] { 0, 4, 1 });
            p.AddCoeff(-1.3466623372930600e+02, new int[] { 0, 4, 3 });
            p.AddCoeff(1.2119961035637500e+02, new int[] { 0, 4, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("efa7f422-ab26-487d-99ba-87cee1cacf33"));
            OrthonormalPolynomials[170] = p;
            p.AddCoeff(2.4734614358443900e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-2.4734614358443900e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(2.8857050084851200e+01, new int[] { 0, 1, 4 });
            p.AddCoeff(-1.1542820033940500e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(1.1542820033940500e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(-1.3466623372930600e+02, new int[] { 0, 3, 4 });
            p.AddCoeff(1.0388538030546400e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(-1.0388538030546400e+02, new int[] { 0, 5, 2 });
            p.AddCoeff(1.2119961035637500e+02, new int[] { 0, 5, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("4472b9b9-a6ed-4b66-8b8f-202906aee318"));
            OrthonormalPolynomials[171] = p;
            p.AddCoeff(1.5809458081912500e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-2.6349096803187500e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(-3.3199861972016200e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(5.5333103286693800e+01, new int[] { 0, 2, 3 });
            p.AddCoeff(9.9599585916048800e+01, new int[] { 0, 4, 1 });
            p.AddCoeff(-1.6599930986008100e+02, new int[] { 0, 4, 3 });
            p.AddCoeff(-7.3039696338435800e+01, new int[] { 0, 6, 1 });
            p.AddCoeff(1.2173282723072600e+02, new int[] { 0, 6, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e593d6e2-d158-401b-ac17-03978795585b"));
            OrthonormalPolynomials[172] = p;
            p.AddCoeff(3.3489117577113800e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-1.0046735273134100e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-3.0140205819402400e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(9.0420617458207200e+01, new int[] { 0, 3, 2 });
            p.AddCoeff(6.6308452802685300e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(-1.9892535840805600e+02, new int[] { 0, 5, 2 });
            p.AddCoeff(-4.1048089830233800e+01, new int[] { 0, 7, 0 });
            p.AddCoeff(1.2314426949070100e+02, new int[] { 0, 7, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("522dec44-e355-426f-9297-35e4fddca9f3"));
            OrthonormalPolynomials[173] = p;
            p.AddCoeff(6.9039584539584700e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(-2.4854250434250500e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(1.3669837738837800e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(-2.3694385413985500e+02, new int[] { 0, 6, 1 });
            p.AddCoeff(1.2693420757492200e+02, new int[] { 0, 8, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("ea0c8c75-5506-4c1c-af14-23c405fe08c1"));
            OrthonormalPolynomials[174] = p;
            p.AddCoeff(3.7925593963578700e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-5.5624204479915400e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(2.1693439747167000e+02, new int[] { 0, 5, 0 });
            p.AddCoeff(-3.0990628210238600e+02, new int[] { 0, 7, 0 });
            p.AddCoeff(1.4634463321501600e+02, new int[] { 0, 9, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("26410743-ccb2-470b-b638-8d57cbc51be8"));
            OrthonormalPolynomials[175] = p;
            p.AddCoeff(6.9039584539584700e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(-2.4854250434250500e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(1.3669837738837800e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(-2.3694385413985500e+02, new int[] { 1, 0, 6 });
            p.AddCoeff(1.2693420757492200e+02, new int[] { 1, 0, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("99ea33d6-7553-4aaf-bd79-bbdfb69dc266"));
            OrthonormalPolynomials[176] = p;
            p.AddCoeff(-8.9860732090691300e+00, new int[] { 1, 1, 1 });
            p.AddCoeff(8.0874658881622200e+01, new int[] { 1, 1, 3 });
            p.AddCoeff(-1.7792424953956900e+02, new int[] { 1, 1, 5 });
            p.AddCoeff(1.1014358304830400e+02, new int[] { 1, 1, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("957f8743-0cce-4a56-901b-fff18479d889"));
            OrthonormalPolynomials[177] = p;
            p.AddCoeff(7.7142256477076200e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(-1.6199873860186000e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(4.8599621580558000e+01, new int[] { 1, 0, 4 });
            p.AddCoeff(-3.5639722492409200e+01, new int[] { 1, 0, 6 });
            p.AddCoeff(-2.3142676943122900e+00, new int[] { 1, 2, 0 });
            p.AddCoeff(4.8599621580558000e+01, new int[] { 1, 2, 2 });
            p.AddCoeff(-1.4579886474167400e+02, new int[] { 1, 2, 4 });
            p.AddCoeff(1.0691916747722800e+02, new int[] { 1, 2, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("757fe792-bf5c-426f-b37a-07c1adf37986"));
            OrthonormalPolynomials[178] = p;
            p.AddCoeff(-1.5113099011081400e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(7.0527795385046700e+01, new int[] { 1, 1, 3 });
            p.AddCoeff(-6.3475015846542000e+01, new int[] { 1, 1, 5 });
            p.AddCoeff(2.5188498351802400e+01, new int[] { 1, 3, 1 });
            p.AddCoeff(-1.1754632564174400e+02, new int[] { 1, 3, 3 });
            p.AddCoeff(1.0579169307757000e+02, new int[] { 1, 3, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e7f15a01-ac51-46d5-99f7-1841aec78b68"));
            OrthonormalPolynomials[179] = p;
            p.AddCoeff(7.7503386392749100e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(-7.7503386392749100e+00, new int[] { 1, 0, 2 });
            p.AddCoeff(9.0420617458207200e+00, new int[] { 1, 0, 4 });
            p.AddCoeff(-7.7503386392749100e+00, new int[] { 1, 2, 0 });
            p.AddCoeff(7.7503386392749100e+01, new int[] { 1, 2, 2 });
            p.AddCoeff(-9.0420617458207200e+01, new int[] { 1, 2, 4 });
            p.AddCoeff(9.0420617458207200e+00, new int[] { 1, 4, 0 });
            p.AddCoeff(-9.0420617458207200e+01, new int[] { 1, 4, 2 });
            p.AddCoeff(1.0549072036790800e+02, new int[] { 1, 4, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("4d2178b9-4ed9-4366-a818-f36aa0034e1a"));
            OrthonormalPolynomials[180] = p;
            p.AddCoeff(-1.5113099011081400e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(2.5188498351802400e+01, new int[] { 1, 1, 3 });
            p.AddCoeff(7.0527795385046700e+01, new int[] { 1, 3, 1 });
            p.AddCoeff(-1.1754632564174400e+02, new int[] { 1, 3, 3 });
            p.AddCoeff(-6.3475015846542000e+01, new int[] { 1, 5, 1 });
            p.AddCoeff(1.0579169307757000e+02, new int[] { 1, 5, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("215a7899-1b4f-43aa-8ea2-e734f0d9954e"));
            OrthonormalPolynomials[181] = p;
            p.AddCoeff(7.7142256477076200e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(-2.3142676943122900e+00, new int[] { 1, 0, 2 });
            p.AddCoeff(-1.6199873860186000e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(4.8599621580558000e+01, new int[] { 1, 2, 2 });
            p.AddCoeff(4.8599621580558000e+01, new int[] { 1, 4, 0 });
            p.AddCoeff(-1.4579886474167400e+02, new int[] { 1, 4, 2 });
            p.AddCoeff(-3.5639722492409200e+01, new int[] { 1, 6, 0 });
            p.AddCoeff(1.0691916747722800e+02, new int[] { 1, 6, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("17132c13-c0bf-4ff3-a7c6-7042cdc33cbe"));
            OrthonormalPolynomials[182] = p;
            p.AddCoeff(-8.9860732090691300e+00, new int[] { 1, 1, 1 });
            p.AddCoeff(8.0874658881622200e+01, new int[] { 1, 3, 1 });
            p.AddCoeff(-1.7792424953956900e+02, new int[] { 1, 5, 1 });
            p.AddCoeff(1.1014358304830400e+02, new int[] { 1, 7, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d9f5c6c1-bf31-4cc3-826e-86fd1fdf09d7"));
            OrthonormalPolynomials[183] = p;
            p.AddCoeff(6.9039584539584700e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(-2.4854250434250500e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(1.3669837738837800e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(-2.3694385413985500e+02, new int[] { 1, 6, 0 });
            p.AddCoeff(1.2693420757492200e+02, new int[] { 1, 8, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("361aee67-fbe7-4ebb-b33b-5645754d8264"));
            OrthonormalPolynomials[184] = p;
            p.AddCoeff(3.3489117577113800e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-3.0140205819402400e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(6.6308452802685300e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(-4.1048089830233800e+01, new int[] { 0, 0, 7 });
            p.AddCoeff(-1.0046735273134100e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(9.0420617458207200e+01, new int[] { 2, 0, 3 });
            p.AddCoeff(-1.9892535840805600e+02, new int[] { 2, 0, 5 });
            p.AddCoeff(1.2314426949070100e+02, new int[] { 2, 0, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("252a4fde-4d30-4fcc-ba97-372e4a9e44a8"));
            OrthonormalPolynomials[185] = p;
            p.AddCoeff(7.7142256477076200e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(-1.6199873860186000e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(4.8599621580558000e+01, new int[] { 0, 1, 4 });
            p.AddCoeff(-3.5639722492409200e+01, new int[] { 0, 1, 6 });
            p.AddCoeff(-2.3142676943122900e+00, new int[] { 2, 1, 0 });
            p.AddCoeff(4.8599621580558000e+01, new int[] { 2, 1, 2 });
            p.AddCoeff(-1.4579886474167400e+02, new int[] { 2, 1, 4 });
            p.AddCoeff(1.0691916747722800e+02, new int[] { 2, 1, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("07408081-258f-4778-acc9-ce5721c0499d"));
            OrthonormalPolynomials[186] = p;
            p.AddCoeff(2.7482904842715400e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-1.2825355593267200e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(1.1542820033940500e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(-8.2448714528146200e+00, new int[] { 0, 2, 1 });
            p.AddCoeff(3.8476066779801600e+01, new int[] { 0, 2, 3 });
            p.AddCoeff(-3.4628460101821400e+01, new int[] { 0, 2, 5 });
            p.AddCoeff(-8.2448714528146200e+00, new int[] { 2, 0, 1 });
            p.AddCoeff(3.8476066779801600e+01, new int[] { 2, 0, 3 });
            p.AddCoeff(-3.4628460101821400e+01, new int[] { 2, 0, 5 });
            p.AddCoeff(2.4734614358443900e+01, new int[] { 2, 2, 1 });
            p.AddCoeff(-1.1542820033940500e+02, new int[] { 2, 2, 3 });
            p.AddCoeff(1.0388538030546400e+02, new int[] { 2, 2, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("baac566c-6801-4a3d-91a7-95f86059c8ac"));
            OrthonormalPolynomials[187] = p;
            p.AddCoeff(1.7648297434703200e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-1.7648297434703200e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(2.0589680340487000e+01, new int[] { 0, 1, 4 });
            p.AddCoeff(-2.9413829057838600e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(2.9413829057838600e+01, new int[] { 0, 3, 2 });
            p.AddCoeff(-3.4316133900811700e+01, new int[] { 0, 3, 4 });
            p.AddCoeff(-5.2944892304109500e+00, new int[] { 2, 1, 0 });
            p.AddCoeff(5.2944892304109500e+01, new int[] { 2, 1, 2 });
            p.AddCoeff(-6.1769041021461100e+01, new int[] { 2, 1, 4 });
            p.AddCoeff(8.8241487173515800e+00, new int[] { 2, 3, 0 });
            p.AddCoeff(-8.8241487173515800e+01, new int[] { 2, 3, 2 });
            p.AddCoeff(1.0294840170243500e+02, new int[] { 2, 3, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("9849bc8c-7d33-41eb-a955-526bca875f30"));
            OrthonormalPolynomials[188] = p;
            p.AddCoeff(1.7648297434703200e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-2.9413829057838600e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(-1.7648297434703200e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(2.9413829057838600e+01, new int[] { 0, 2, 3 });
            p.AddCoeff(2.0589680340487000e+01, new int[] { 0, 4, 1 });
            p.AddCoeff(-3.4316133900811700e+01, new int[] { 0, 4, 3 });
            p.AddCoeff(-5.2944892304109500e+00, new int[] { 2, 0, 1 });
            p.AddCoeff(8.8241487173515800e+00, new int[] { 2, 0, 3 });
            p.AddCoeff(5.2944892304109500e+01, new int[] { 2, 2, 1 });
            p.AddCoeff(-8.8241487173515800e+01, new int[] { 2, 2, 3 });
            p.AddCoeff(-6.1769041021461100e+01, new int[] { 2, 4, 1 });
            p.AddCoeff(1.0294840170243500e+02, new int[] { 2, 4, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f992c4dd-e2b5-4cfb-8836-38787d88bd80"));
            OrthonormalPolynomials[189] = p;
            p.AddCoeff(2.7482904842715400e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-8.2448714528146200e+00, new int[] { 0, 1, 2 });
            p.AddCoeff(-1.2825355593267200e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(3.8476066779801600e+01, new int[] { 0, 3, 2 });
            p.AddCoeff(1.1542820033940500e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(-3.4628460101821400e+01, new int[] { 0, 5, 2 });
            p.AddCoeff(-8.2448714528146200e+00, new int[] { 2, 1, 0 });
            p.AddCoeff(2.4734614358443900e+01, new int[] { 2, 1, 2 });
            p.AddCoeff(3.8476066779801600e+01, new int[] { 2, 3, 0 });
            p.AddCoeff(-1.1542820033940500e+02, new int[] { 2, 3, 2 });
            p.AddCoeff(-3.4628460101821400e+01, new int[] { 2, 5, 0 });
            p.AddCoeff(1.0388538030546400e+02, new int[] { 2, 5, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d1d1de04-5650-42d0-afba-f735bf3030a6"));
            OrthonormalPolynomials[190] = p;
            p.AddCoeff(7.7142256477076200e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(-1.6199873860186000e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(4.8599621580558000e+01, new int[] { 0, 4, 1 });
            p.AddCoeff(-3.5639722492409200e+01, new int[] { 0, 6, 1 });
            p.AddCoeff(-2.3142676943122900e+00, new int[] { 2, 0, 1 });
            p.AddCoeff(4.8599621580558000e+01, new int[] { 2, 2, 1 });
            p.AddCoeff(-1.4579886474167400e+02, new int[] { 2, 4, 1 });
            p.AddCoeff(1.0691916747722800e+02, new int[] { 2, 6, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3edf9e9c-db1c-4311-a171-8eb8608c9e70"));
            OrthonormalPolynomials[191] = p;
            p.AddCoeff(3.3489117577113800e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-3.0140205819402400e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(6.6308452802685300e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(-4.1048089830233800e+01, new int[] { 0, 7, 0 });
            p.AddCoeff(-1.0046735273134100e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(9.0420617458207200e+01, new int[] { 2, 3, 0 });
            p.AddCoeff(-1.9892535840805600e+02, new int[] { 2, 5, 0 });
            p.AddCoeff(1.2314426949070100e+02, new int[] { 2, 7, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3de36b16-0d59-4bca-886f-97dfede7baf4"));
            OrthonormalPolynomials[192] = p;
            p.AddCoeff(1.5809458081912500e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-3.3199861972016200e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(9.9599585916048800e+01, new int[] { 1, 0, 4 });
            p.AddCoeff(-7.3039696338435800e+01, new int[] { 1, 0, 6 });
            p.AddCoeff(-2.6349096803187500e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(5.5333103286693800e+01, new int[] { 3, 0, 2 });
            p.AddCoeff(-1.6599930986008100e+02, new int[] { 3, 0, 4 });
            p.AddCoeff(1.2173282723072600e+02, new int[] { 3, 0, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("8eb13b35-2980-443c-97eb-cb9978b45bd5"));
            OrthonormalPolynomials[193] = p;
            p.AddCoeff(-1.5113099011081400e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(7.0527795385046700e+01, new int[] { 1, 1, 3 });
            p.AddCoeff(-6.3475015846542000e+01, new int[] { 1, 1, 5 });
            p.AddCoeff(2.5188498351802400e+01, new int[] { 3, 1, 1 });
            p.AddCoeff(-1.1754632564174400e+02, new int[] { 3, 1, 3 });
            p.AddCoeff(1.0579169307757000e+02, new int[] { 3, 1, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a53ebae5-64e1-422f-b17e-1c6179197f18"));
            OrthonormalPolynomials[194] = p;
            p.AddCoeff(1.7648297434703200e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-1.7648297434703200e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(2.0589680340487000e+01, new int[] { 1, 0, 4 });
            p.AddCoeff(-5.2944892304109500e+00, new int[] { 1, 2, 0 });
            p.AddCoeff(5.2944892304109500e+01, new int[] { 1, 2, 2 });
            p.AddCoeff(-6.1769041021461100e+01, new int[] { 1, 2, 4 });
            p.AddCoeff(-2.9413829057838600e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(2.9413829057838600e+01, new int[] { 3, 0, 2 });
            p.AddCoeff(-3.4316133900811700e+01, new int[] { 3, 0, 4 });
            p.AddCoeff(8.8241487173515800e+00, new int[] { 3, 2, 0 });
            p.AddCoeff(-8.8241487173515800e+01, new int[] { 3, 2, 2 });
            p.AddCoeff(1.0294840170243500e+02, new int[] { 3, 2, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("db612381-dec8-4e14-af7a-297cb0680b40"));
            OrthonormalPolynomials[195] = p;
            p.AddCoeff(-2.2099163940633600e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(3.6831939901056000e+01, new int[] { 1, 1, 3 });
            p.AddCoeff(3.6831939901056000e+01, new int[] { 1, 3, 1 });
            p.AddCoeff(-6.1386566501760000e+01, new int[] { 1, 3, 3 });
            p.AddCoeff(3.6831939901056000e+01, new int[] { 3, 1, 1 });
            p.AddCoeff(-6.1386566501760000e+01, new int[] { 3, 1, 3 });
            p.AddCoeff(-6.1386566501760000e+01, new int[] { 3, 3, 1 });
            p.AddCoeff(1.0231094416960000e+02, new int[] { 3, 3, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("37723cf6-712c-4ee1-970c-713694df426e"));
            OrthonormalPolynomials[196] = p;
            p.AddCoeff(1.7648297434703200e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-5.2944892304109500e+00, new int[] { 1, 0, 2 });
            p.AddCoeff(-1.7648297434703200e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(5.2944892304109500e+01, new int[] { 1, 2, 2 });
            p.AddCoeff(2.0589680340487000e+01, new int[] { 1, 4, 0 });
            p.AddCoeff(-6.1769041021461100e+01, new int[] { 1, 4, 2 });
            p.AddCoeff(-2.9413829057838600e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(8.8241487173515800e+00, new int[] { 3, 0, 2 });
            p.AddCoeff(2.9413829057838600e+01, new int[] { 3, 2, 0 });
            p.AddCoeff(-8.8241487173515800e+01, new int[] { 3, 2, 2 });
            p.AddCoeff(-3.4316133900811700e+01, new int[] { 3, 4, 0 });
            p.AddCoeff(1.0294840170243500e+02, new int[] { 3, 4, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e55ce590-e924-482e-94d4-9913c9e3673c"));
            OrthonormalPolynomials[197] = p;
            p.AddCoeff(-1.5113099011081400e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(7.0527795385046700e+01, new int[] { 1, 3, 1 });
            p.AddCoeff(-6.3475015846542000e+01, new int[] { 1, 5, 1 });
            p.AddCoeff(2.5188498351802400e+01, new int[] { 3, 1, 1 });
            p.AddCoeff(-1.1754632564174400e+02, new int[] { 3, 3, 1 });
            p.AddCoeff(1.0579169307757000e+02, new int[] { 3, 5, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f4d397ec-9c2b-49b6-9fde-bc2e5fbfec3a"));
            OrthonormalPolynomials[198] = p;
            p.AddCoeff(1.5809458081912500e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-3.3199861972016200e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(9.9599585916048800e+01, new int[] { 1, 4, 0 });
            p.AddCoeff(-7.3039696338435800e+01, new int[] { 1, 6, 0 });
            p.AddCoeff(-2.6349096803187500e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(5.5333103286693800e+01, new int[] { 3, 2, 0 });
            p.AddCoeff(-1.6599930986008100e+02, new int[] { 3, 4, 0 });
            p.AddCoeff(1.2173282723072600e+02, new int[] { 3, 6, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5f2bc6a2-d83a-4429-9f11-dc80629e3509"));
            OrthonormalPolynomials[199] = p;
            p.AddCoeff(2.4734614358443900e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-1.1542820033940500e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(1.0388538030546400e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(-2.4734614358443900e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(1.1542820033940500e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(-1.0388538030546400e+02, new int[] { 2, 0, 5 });
            p.AddCoeff(2.8857050084851200e+01, new int[] { 4, 0, 1 });
            p.AddCoeff(-1.3466623372930600e+02, new int[] { 4, 0, 3 });
            p.AddCoeff(1.2119961035637500e+02, new int[] { 4, 0, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e86da69d-fe84-4f6b-a2e3-8b456fdf56ec"));
            OrthonormalPolynomials[200] = p;
            p.AddCoeff(7.7503386392749100e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(-7.7503386392749100e+00, new int[] { 0, 1, 2 });
            p.AddCoeff(9.0420617458207200e+00, new int[] { 0, 1, 4 });
            p.AddCoeff(-7.7503386392749100e+00, new int[] { 2, 1, 0 });
            p.AddCoeff(7.7503386392749100e+01, new int[] { 2, 1, 2 });
            p.AddCoeff(-9.0420617458207200e+01, new int[] { 2, 1, 4 });
            p.AddCoeff(9.0420617458207200e+00, new int[] { 4, 1, 0 });
            p.AddCoeff(-9.0420617458207200e+01, new int[] { 4, 1, 2 });
            p.AddCoeff(1.0549072036790800e+02, new int[] { 4, 1, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d4480211-fb0b-48c5-9273-74cb75edc20e"));
            OrthonormalPolynomials[201] = p;
            p.AddCoeff(1.7648297434703200e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-2.9413829057838600e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(-5.2944892304109500e+00, new int[] { 0, 2, 1 });
            p.AddCoeff(8.8241487173515800e+00, new int[] { 0, 2, 3 });
            p.AddCoeff(-1.7648297434703200e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(2.9413829057838600e+01, new int[] { 2, 0, 3 });
            p.AddCoeff(5.2944892304109500e+01, new int[] { 2, 2, 1 });
            p.AddCoeff(-8.8241487173515800e+01, new int[] { 2, 2, 3 });
            p.AddCoeff(2.0589680340487000e+01, new int[] { 4, 0, 1 });
            p.AddCoeff(-3.4316133900811700e+01, new int[] { 4, 0, 3 });
            p.AddCoeff(-6.1769041021461100e+01, new int[] { 4, 2, 1 });
            p.AddCoeff(1.0294840170243500e+02, new int[] { 4, 2, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3accb26a-f05b-431f-8e18-c4877f6e563a"));
            OrthonormalPolynomials[202] = p;
            p.AddCoeff(1.7648297434703200e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-5.2944892304109500e+00, new int[] { 0, 1, 2 });
            p.AddCoeff(-2.9413829057838600e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(8.8241487173515800e+00, new int[] { 0, 3, 2 });
            p.AddCoeff(-1.7648297434703200e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(5.2944892304109500e+01, new int[] { 2, 1, 2 });
            p.AddCoeff(2.9413829057838600e+01, new int[] { 2, 3, 0 });
            p.AddCoeff(-8.8241487173515800e+01, new int[] { 2, 3, 2 });
            p.AddCoeff(2.0589680340487000e+01, new int[] { 4, 1, 0 });
            p.AddCoeff(-6.1769041021461100e+01, new int[] { 4, 1, 2 });
            p.AddCoeff(-3.4316133900811700e+01, new int[] { 4, 3, 0 });
            p.AddCoeff(1.0294840170243500e+02, new int[] { 4, 3, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d4d033f0-40ce-4df0-9b69-117a69f28181"));
            OrthonormalPolynomials[203] = p;
            p.AddCoeff(7.7503386392749100e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(-7.7503386392749100e+00, new int[] { 0, 2, 1 });
            p.AddCoeff(9.0420617458207200e+00, new int[] { 0, 4, 1 });
            p.AddCoeff(-7.7503386392749100e+00, new int[] { 2, 0, 1 });
            p.AddCoeff(7.7503386392749100e+01, new int[] { 2, 2, 1 });
            p.AddCoeff(-9.0420617458207200e+01, new int[] { 2, 4, 1 });
            p.AddCoeff(9.0420617458207200e+00, new int[] { 4, 0, 1 });
            p.AddCoeff(-9.0420617458207200e+01, new int[] { 4, 2, 1 });
            p.AddCoeff(1.0549072036790800e+02, new int[] { 4, 4, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("8528383f-fec0-4daa-a4be-1a587c9495b2"));
            OrthonormalPolynomials[204] = p;
            p.AddCoeff(2.4734614358443900e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-1.1542820033940500e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(1.0388538030546400e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(-2.4734614358443900e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(1.1542820033940500e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(-1.0388538030546400e+02, new int[] { 2, 5, 0 });
            p.AddCoeff(2.8857050084851200e+01, new int[] { 4, 1, 0 });
            p.AddCoeff(-1.3466623372930600e+02, new int[] { 4, 3, 0 });
            p.AddCoeff(1.2119961035637500e+02, new int[] { 4, 5, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("0986e831-e14a-41cc-9fa2-d246a77f06b8"));
            OrthonormalPolynomials[205] = p;
            p.AddCoeff(2.4734614358443900e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-2.4734614358443900e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(2.8857050084851200e+01, new int[] { 1, 0, 4 });
            p.AddCoeff(-1.1542820033940500e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(1.1542820033940500e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(-1.3466623372930600e+02, new int[] { 3, 0, 4 });
            p.AddCoeff(1.0388538030546400e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(-1.0388538030546400e+02, new int[] { 5, 0, 2 });
            p.AddCoeff(1.2119961035637500e+02, new int[] { 5, 0, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5a5e913c-fd69-4631-8351-105854279456"));
            OrthonormalPolynomials[206] = p;
            p.AddCoeff(-1.5113099011081400e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(2.5188498351802400e+01, new int[] { 1, 1, 3 });
            p.AddCoeff(7.0527795385046700e+01, new int[] { 3, 1, 1 });
            p.AddCoeff(-1.1754632564174400e+02, new int[] { 3, 1, 3 });
            p.AddCoeff(-6.3475015846542000e+01, new int[] { 5, 1, 1 });
            p.AddCoeff(1.0579169307757000e+02, new int[] { 5, 1, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d2ea4f1b-9e9a-4fe7-9688-1173bc569a83"));
            OrthonormalPolynomials[207] = p;
            p.AddCoeff(2.7482904842715400e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-8.2448714528146200e+00, new int[] { 1, 0, 2 });
            p.AddCoeff(-8.2448714528146200e+00, new int[] { 1, 2, 0 });
            p.AddCoeff(2.4734614358443900e+01, new int[] { 1, 2, 2 });
            p.AddCoeff(-1.2825355593267200e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(3.8476066779801600e+01, new int[] { 3, 0, 2 });
            p.AddCoeff(3.8476066779801600e+01, new int[] { 3, 2, 0 });
            p.AddCoeff(-1.1542820033940500e+02, new int[] { 3, 2, 2 });
            p.AddCoeff(1.1542820033940500e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(-3.4628460101821400e+01, new int[] { 5, 0, 2 });
            p.AddCoeff(-3.4628460101821400e+01, new int[] { 5, 2, 0 });
            p.AddCoeff(1.0388538030546400e+02, new int[] { 5, 2, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("2ac9875d-068e-4a96-b296-fde6362704fa"));
            OrthonormalPolynomials[208] = p;
            p.AddCoeff(-1.5113099011081400e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(2.5188498351802400e+01, new int[] { 1, 3, 1 });
            p.AddCoeff(7.0527795385046700e+01, new int[] { 3, 1, 1 });
            p.AddCoeff(-1.1754632564174400e+02, new int[] { 3, 3, 1 });
            p.AddCoeff(-6.3475015846542000e+01, new int[] { 5, 1, 1 });
            p.AddCoeff(1.0579169307757000e+02, new int[] { 5, 3, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7ec67f23-c81a-4889-81ff-0874f17816d6"));
            OrthonormalPolynomials[209] = p;
            p.AddCoeff(2.4734614358443900e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-2.4734614358443900e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(2.8857050084851200e+01, new int[] { 1, 4, 0 });
            p.AddCoeff(-1.1542820033940500e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(1.1542820033940500e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(-1.3466623372930600e+02, new int[] { 3, 4, 0 });
            p.AddCoeff(1.0388538030546400e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(-1.0388538030546400e+02, new int[] { 5, 2, 0 });
            p.AddCoeff(1.2119961035637500e+02, new int[] { 5, 4, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b44cefd3-eb00-4ac9-b602-0d3a61e126e2"));
            OrthonormalPolynomials[210] = p;
            p.AddCoeff(1.5809458081912500e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-2.6349096803187500e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(-3.3199861972016200e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(5.5333103286693800e+01, new int[] { 2, 0, 3 });
            p.AddCoeff(9.9599585916048800e+01, new int[] { 4, 0, 1 });
            p.AddCoeff(-1.6599930986008100e+02, new int[] { 4, 0, 3 });
            p.AddCoeff(-7.3039696338435800e+01, new int[] { 6, 0, 1 });
            p.AddCoeff(1.2173282723072600e+02, new int[] { 6, 0, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f6f596c4-45c0-4367-bfdf-5aa623b396ac"));
            OrthonormalPolynomials[211] = p;
            p.AddCoeff(7.7142256477076200e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(-2.3142676943122900e+00, new int[] { 0, 1, 2 });
            p.AddCoeff(-1.6199873860186000e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(4.8599621580558000e+01, new int[] { 2, 1, 2 });
            p.AddCoeff(4.8599621580558000e+01, new int[] { 4, 1, 0 });
            p.AddCoeff(-1.4579886474167400e+02, new int[] { 4, 1, 2 });
            p.AddCoeff(-3.5639722492409200e+01, new int[] { 6, 1, 0 });
            p.AddCoeff(1.0691916747722800e+02, new int[] { 6, 1, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("4766358e-b4b8-4ef3-8b45-f11f38f2f1dd"));
            OrthonormalPolynomials[212] = p;
            p.AddCoeff(7.7142256477076200e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(-2.3142676943122900e+00, new int[] { 0, 2, 1 });
            p.AddCoeff(-1.6199873860186000e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(4.8599621580558000e+01, new int[] { 2, 2, 1 });
            p.AddCoeff(4.8599621580558000e+01, new int[] { 4, 0, 1 });
            p.AddCoeff(-1.4579886474167400e+02, new int[] { 4, 2, 1 });
            p.AddCoeff(-3.5639722492409200e+01, new int[] { 6, 0, 1 });
            p.AddCoeff(1.0691916747722800e+02, new int[] { 6, 2, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("0218c2d2-9820-483f-a4dd-a3cd3511ef57"));
            OrthonormalPolynomials[213] = p;
            p.AddCoeff(1.5809458081912500e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-2.6349096803187500e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(-3.3199861972016200e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(5.5333103286693800e+01, new int[] { 2, 3, 0 });
            p.AddCoeff(9.9599585916048800e+01, new int[] { 4, 1, 0 });
            p.AddCoeff(-1.6599930986008100e+02, new int[] { 4, 3, 0 });
            p.AddCoeff(-7.3039696338435800e+01, new int[] { 6, 1, 0 });
            p.AddCoeff(1.2173282723072600e+02, new int[] { 6, 3, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c227fb8b-8eb2-47a1-b208-531878797517"));
            OrthonormalPolynomials[214] = p;
            p.AddCoeff(3.3489117577113800e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-1.0046735273134100e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-3.0140205819402400e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(9.0420617458207200e+01, new int[] { 3, 0, 2 });
            p.AddCoeff(6.6308452802685300e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(-1.9892535840805600e+02, new int[] { 5, 0, 2 });
            p.AddCoeff(-4.1048089830233800e+01, new int[] { 7, 0, 0 });
            p.AddCoeff(1.2314426949070100e+02, new int[] { 7, 0, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("62a4d411-ec56-44fe-b479-6919d09d5997"));
            OrthonormalPolynomials[215] = p;
            p.AddCoeff(-8.9860732090691300e+00, new int[] { 1, 1, 1 });
            p.AddCoeff(8.0874658881622200e+01, new int[] { 3, 1, 1 });
            p.AddCoeff(-1.7792424953956900e+02, new int[] { 5, 1, 1 });
            p.AddCoeff(1.1014358304830400e+02, new int[] { 7, 1, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d37a05f0-6347-4435-a66d-cd49caeb987e"));
            OrthonormalPolynomials[216] = p;
            p.AddCoeff(3.3489117577113800e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-1.0046735273134100e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-3.0140205819402400e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(9.0420617458207200e+01, new int[] { 3, 2, 0 });
            p.AddCoeff(6.6308452802685300e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(-1.9892535840805600e+02, new int[] { 5, 2, 0 });
            p.AddCoeff(-4.1048089830233800e+01, new int[] { 7, 0, 0 });
            p.AddCoeff(1.2314426949070100e+02, new int[] { 7, 2, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("8caab61d-f96c-468c-8c5d-906a1667dfb0"));
            OrthonormalPolynomials[217] = p;
            p.AddCoeff(6.9039584539584700e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(-2.4854250434250500e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(1.3669837738837800e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(-2.3694385413985500e+02, new int[] { 6, 0, 1 });
            p.AddCoeff(1.2693420757492200e+02, new int[] { 8, 0, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("60772c1c-f28c-4611-8601-1bab2cb66056"));
            OrthonormalPolynomials[218] = p;
            p.AddCoeff(6.9039584539584700e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(-2.4854250434250500e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(1.3669837738837800e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(-2.3694385413985500e+02, new int[] { 6, 1, 0 });
            p.AddCoeff(1.2693420757492200e+02, new int[] { 8, 1, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d8af0a09-b028-43c8-86bd-4f660d900fd7"));
            OrthonormalPolynomials[219] = p;
            p.AddCoeff(3.7925593963578700e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-5.5624204479915400e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(2.1693439747167000e+02, new int[] { 5, 0, 0 });
            p.AddCoeff(-3.0990628210238600e+02, new int[] { 7, 0, 0 });
            p.AddCoeff(1.4634463321501600e+02, new int[] { 9, 0, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7c301a3d-c3f9-4f1a-8666-11fc2c88a3c6"));
            OrthonormalPolynomials[220] = p;
            p.AddCoeff(-3.9871744531220200e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(2.1929459492171100e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(-1.9005531559881600e+02, new int[] { 0, 0, 4 });
            p.AddCoeff(5.7016594679644900e+02, new int[] { 0, 0, 6 });
            p.AddCoeff(-6.9234436396711700e+02, new int[] { 0, 0, 8 });
            p.AddCoeff(2.9232317589722700e+02, new int[] { 0, 0, 10 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b12623a0-91aa-4b03-8d04-7eee9fa2acb2"));
            OrthonormalPolynomials[221] = p;
            p.AddCoeff(6.5689055652145700e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-9.6343948289813700e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(3.7574139833027300e+02, new int[] { 0, 1, 5 });
            p.AddCoeff(-5.3677342618610500e+02, new int[] { 0, 1, 7 });
            p.AddCoeff(2.5347634014343800e+02, new int[] { 0, 1, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("81cc7c50-2809-427e-a656-f4cf50f6d6d7"));
            OrthonormalPolynomials[222] = p;
            p.AddCoeff(-4.4564860191815000e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.6043349669053400e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(-8.8238423179793700e+01, new int[] { 0, 0, 4 });
            p.AddCoeff(1.5294660017830900e+02, new int[] { 0, 0, 6 });
            p.AddCoeff(-8.1935678666951300e+01, new int[] { 0, 0, 8 });
            p.AddCoeff(1.3369458057544500e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-4.8130049007160200e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(2.6471526953938100e+02, new int[] { 0, 2, 4 });
            p.AddCoeff(-4.5883980053492700e+02, new int[] { 0, 2, 6 });
            p.AddCoeff(2.4580703600085400e+02, new int[] { 0, 2, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("822d2b0a-a840-498c-a5ec-17f5a4a94d00"));
            OrthonormalPolynomials[223] = p;
            p.AddCoeff(1.1887457487108500e+01, new int[] { 0, 1, 1 });
            p.AddCoeff(-1.0698711738397600e+02, new int[] { 0, 1, 3 });
            p.AddCoeff(2.3537165824474800e+02, new int[] { 0, 1, 5 });
            p.AddCoeff(-1.4570626462770100e+02, new int[] { 0, 1, 7 });
            p.AddCoeff(-1.9812429145180800e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(1.7831186230662700e+02, new int[] { 0, 3, 3 });
            p.AddCoeff(-3.9228609707458000e+02, new int[] { 0, 3, 5 });
            p.AddCoeff(2.4284377437950200e+02, new int[] { 0, 3, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("8e49c50e-365d-4fcb-81a6-9e8ccc48fda7"));
            OrthonormalPolynomials[224] = p;
            p.AddCoeff(-4.4815601193686500e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(9.4112762506741700e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-2.8233828752022500e+01, new int[] { 0, 0, 4 });
            p.AddCoeff(2.0704807751483200e+01, new int[] { 0, 0, 6 });
            p.AddCoeff(4.4815601193686500e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-9.4112762506741700e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(2.8233828752022500e+02, new int[] { 0, 2, 4 });
            p.AddCoeff(-2.0704807751483200e+02, new int[] { 0, 2, 6 });
            p.AddCoeff(-5.2284868059301000e+00, new int[] { 0, 4, 0 });
            p.AddCoeff(1.0979822292453200e+02, new int[] { 0, 4, 2 });
            p.AddCoeff(-3.2939466877359600e+02, new int[] { 0, 4, 4 });
            p.AddCoeff(2.4155609043397000e+02, new int[] { 0, 4, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b3167ea2-43b5-4ed9-859a-5ee4db90af5f"));
            OrthonormalPolynomials[225] = p;
            p.AddCoeff(1.3672572526849300e+01, new int[] { 0, 1, 1 });
            p.AddCoeff(-6.3805338458630100e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(5.7424804612767100e+01, new int[] { 0, 1, 5 });
            p.AddCoeff(-6.3805338458630100e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(2.9775824614027400e+02, new int[] { 0, 3, 3 });
            p.AddCoeff(-2.6798242152624600e+02, new int[] { 0, 3, 5 });
            p.AddCoeff(5.7424804612767100e+01, new int[] { 0, 5, 1 });
            p.AddCoeff(-2.6798242152624600e+02, new int[] { 0, 5, 3 });
            p.AddCoeff(2.4118417937362200e+02, new int[] { 0, 5, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b5ba58a7-19b3-4242-8360-eda9b6d2cd02"));
            OrthonormalPolynomials[226] = p;
            p.AddCoeff(-4.4815601193686500e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(4.4815601193686500e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-5.2284868059301000e+00, new int[] { 0, 0, 4 });
            p.AddCoeff(9.4112762506741700e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-9.4112762506741700e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(1.0979822292453200e+02, new int[] { 0, 2, 4 });
            p.AddCoeff(-2.8233828752022500e+01, new int[] { 0, 4, 0 });
            p.AddCoeff(2.8233828752022500e+02, new int[] { 0, 4, 2 });
            p.AddCoeff(-3.2939466877359600e+02, new int[] { 0, 4, 4 });
            p.AddCoeff(2.0704807751483200e+01, new int[] { 0, 6, 0 });
            p.AddCoeff(-2.0704807751483200e+02, new int[] { 0, 6, 2 });
            p.AddCoeff(2.4155609043397000e+02, new int[] { 0, 6, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d13b9ff8-6562-448e-8344-2b53e37a3fdf"));
            OrthonormalPolynomials[227] = p;
            p.AddCoeff(1.1887457487108500e+01, new int[] { 0, 1, 1 });
            p.AddCoeff(-1.9812429145180800e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(-1.0698711738397600e+02, new int[] { 0, 3, 1 });
            p.AddCoeff(1.7831186230662700e+02, new int[] { 0, 3, 3 });
            p.AddCoeff(2.3537165824474800e+02, new int[] { 0, 5, 1 });
            p.AddCoeff(-3.9228609707458000e+02, new int[] { 0, 5, 3 });
            p.AddCoeff(-1.4570626462770100e+02, new int[] { 0, 7, 1 });
            p.AddCoeff(2.4284377437950200e+02, new int[] { 0, 7, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("69b92f87-a3db-43c7-a6c6-2ab873443552"));
            OrthonormalPolynomials[228] = p;
            p.AddCoeff(-4.4564860191815000e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.3369458057544500e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(1.6043349669053400e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(-4.8130049007160200e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(-8.8238423179793700e+01, new int[] { 0, 4, 0 });
            p.AddCoeff(2.6471526953938100e+02, new int[] { 0, 4, 2 });
            p.AddCoeff(1.5294660017830900e+02, new int[] { 0, 6, 0 });
            p.AddCoeff(-4.5883980053492700e+02, new int[] { 0, 6, 2 });
            p.AddCoeff(-8.1935678666951300e+01, new int[] { 0, 8, 0 });
            p.AddCoeff(2.4580703600085400e+02, new int[] { 0, 8, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("400d0ce8-802d-42d8-81ce-bbdaf2e3587b"));
            OrthonormalPolynomials[229] = p;
            p.AddCoeff(6.5689055652145700e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-9.6343948289813700e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(3.7574139833027300e+02, new int[] { 0, 5, 1 });
            p.AddCoeff(-5.3677342618610500e+02, new int[] { 0, 7, 1 });
            p.AddCoeff(2.5347634014343800e+02, new int[] { 0, 9, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5037c02c-447e-4b0c-ab41-201b377e2e1a"));
            OrthonormalPolynomials[230] = p;
            p.AddCoeff(-3.9871744531220200e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(2.1929459492171100e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(-1.9005531559881600e+02, new int[] { 0, 4, 0 });
            p.AddCoeff(5.7016594679644900e+02, new int[] { 0, 6, 0 });
            p.AddCoeff(-6.9234436396711700e+02, new int[] { 0, 8, 0 });
            p.AddCoeff(2.9232317589722700e+02, new int[] { 0, 10, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("ef9f9c32-e538-4dfd-a4dc-b98e787ad651"));
            OrthonormalPolynomials[231] = p;
            p.AddCoeff(6.5689055652145700e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-9.6343948289813700e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(3.7574139833027300e+02, new int[] { 1, 0, 5 });
            p.AddCoeff(-5.3677342618610500e+02, new int[] { 1, 0, 7 });
            p.AddCoeff(2.5347634014343800e+02, new int[] { 1, 0, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c8f9e120-4841-4e7a-95b7-c17e4f7bde53"));
            OrthonormalPolynomials[232] = p;
            p.AddCoeff(1.1958006815600700e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-4.3048824536162600e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(2.3676853494889400e+02, new int[] { 1, 1, 4 });
            p.AddCoeff(-4.1039879391141600e+02, new int[] { 1, 1, 6 });
            p.AddCoeff(2.1985649673825900e+02, new int[] { 1, 1, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a7604f8e-cbfb-4ed5-be80-d9c6afe37736"));
            OrthonormalPolynomials[233] = p;
            p.AddCoeff(5.8004853144209200e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-5.2204367829788300e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(1.1484960922553400e+02, new int[] { 1, 0, 5 });
            p.AddCoeff(-7.1097377139616400e+01, new int[] { 1, 0, 7 });
            p.AddCoeff(-1.7401455943262800e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(1.5661310348936500e+02, new int[] { 1, 2, 3 });
            p.AddCoeff(-3.4454882767660300e+02, new int[] { 1, 2, 5 });
            p.AddCoeff(2.1329213141884900e+02, new int[] { 1, 2, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("bddbf04e-4d3f-4de8-b63f-d04e4230b681"));
            OrthonormalPolynomials[234] = p;
            p.AddCoeff(2.7382784638002900e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-5.7503847739806100e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(1.7251154321941800e+02, new int[] { 1, 1, 4 });
            p.AddCoeff(-1.2650846502757300e+02, new int[] { 1, 1, 6 });
            p.AddCoeff(-4.5637974396671500e+00, new int[] { 1, 3, 0 });
            p.AddCoeff(9.5839746233010100e+01, new int[] { 1, 3, 2 });
            p.AddCoeff(-2.8751923869903000e+02, new int[] { 1, 3, 4 });
            p.AddCoeff(2.1084744171262200e+02, new int[] { 1, 3, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a45df22c-fb41-4abd-a705-32ba6513b183"));
            OrthonormalPolynomials[235] = p;
            p.AddCoeff(4.2841608774447400e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-1.9992750761408800e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(1.7993475685267900e+01, new int[] { 1, 0, 5 });
            p.AddCoeff(-4.2841608774447400e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(1.9992750761408800e+02, new int[] { 1, 2, 3 });
            p.AddCoeff(-1.7993475685267900e+02, new int[] { 1, 2, 5 });
            p.AddCoeff(4.9981876903522000e+01, new int[] { 1, 4, 1 });
            p.AddCoeff(-2.3324875888310300e+02, new int[] { 1, 4, 3 });
            p.AddCoeff(2.0992388299479200e+02, new int[] { 1, 4, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e179661e-419d-4db0-b72f-f79236a59095"));
            OrthonormalPolynomials[236] = p;
            p.AddCoeff(4.2841608774447400e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-4.2841608774447400e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(4.9981876903522000e+01, new int[] { 1, 1, 4 });
            p.AddCoeff(-1.9992750761408800e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(1.9992750761408800e+02, new int[] { 1, 3, 2 });
            p.AddCoeff(-2.3324875888310300e+02, new int[] { 1, 3, 4 });
            p.AddCoeff(1.7993475685267900e+01, new int[] { 1, 5, 0 });
            p.AddCoeff(-1.7993475685267900e+02, new int[] { 1, 5, 2 });
            p.AddCoeff(2.0992388299479200e+02, new int[] { 1, 5, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("26f38bf1-0727-4c57-9d71-6e1fde4e1693"));
            OrthonormalPolynomials[237] = p;
            p.AddCoeff(2.7382784638002900e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-4.5637974396671500e+00, new int[] { 1, 0, 3 });
            p.AddCoeff(-5.7503847739806100e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(9.5839746233010100e+01, new int[] { 1, 2, 3 });
            p.AddCoeff(1.7251154321941800e+02, new int[] { 1, 4, 1 });
            p.AddCoeff(-2.8751923869903000e+02, new int[] { 1, 4, 3 });
            p.AddCoeff(-1.2650846502757300e+02, new int[] { 1, 6, 1 });
            p.AddCoeff(2.1084744171262200e+02, new int[] { 1, 6, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("64f4ca72-28e0-4486-8768-806c5316f68c"));
            OrthonormalPolynomials[238] = p;
            p.AddCoeff(5.8004853144209200e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-1.7401455943262800e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(-5.2204367829788300e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(1.5661310348936500e+02, new int[] { 1, 3, 2 });
            p.AddCoeff(1.1484960922553400e+02, new int[] { 1, 5, 0 });
            p.AddCoeff(-3.4454882767660300e+02, new int[] { 1, 5, 2 });
            p.AddCoeff(-7.1097377139616400e+01, new int[] { 1, 7, 0 });
            p.AddCoeff(2.1329213141884900e+02, new int[] { 1, 7, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d3200f2f-5436-4cc5-9b01-846e370dc08c"));
            OrthonormalPolynomials[239] = p;
            p.AddCoeff(1.1958006815600700e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-4.3048824536162600e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(2.3676853494889400e+02, new int[] { 1, 4, 1 });
            p.AddCoeff(-4.1039879391141600e+02, new int[] { 1, 6, 1 });
            p.AddCoeff(2.1985649673825900e+02, new int[] { 1, 8, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("000ca062-a1d0-43d7-9374-a9244d20df25"));
            OrthonormalPolynomials[240] = p;
            p.AddCoeff(6.5689055652145700e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-9.6343948289813700e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(3.7574139833027300e+02, new int[] { 1, 5, 0 });
            p.AddCoeff(-5.3677342618610500e+02, new int[] { 1, 7, 0 });
            p.AddCoeff(2.5347634014343800e+02, new int[] { 1, 9, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("be2cef55-c6ed-4682-ae20-dd4d3fbbf0ab"));
            OrthonormalPolynomials[241] = p;
            p.AddCoeff(-4.4564860191815000e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.6043349669053400e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(-8.8238423179793700e+01, new int[] { 0, 0, 4 });
            p.AddCoeff(1.5294660017830900e+02, new int[] { 0, 0, 6 });
            p.AddCoeff(-8.1935678666951300e+01, new int[] { 0, 0, 8 });
            p.AddCoeff(1.3369458057544500e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-4.8130049007160200e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(2.6471526953938100e+02, new int[] { 2, 0, 4 });
            p.AddCoeff(-4.5883980053492700e+02, new int[] { 2, 0, 6 });
            p.AddCoeff(2.4580703600085400e+02, new int[] { 2, 0, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a4a1d062-f760-4689-9f6a-b80e5766239e"));
            OrthonormalPolynomials[242] = p;
            p.AddCoeff(5.8004853144209200e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-5.2204367829788300e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(1.1484960922553400e+02, new int[] { 0, 1, 5 });
            p.AddCoeff(-7.1097377139616400e+01, new int[] { 0, 1, 7 });
            p.AddCoeff(-1.7401455943262800e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(1.5661310348936500e+02, new int[] { 2, 1, 3 });
            p.AddCoeff(-3.4454882767660300e+02, new int[] { 2, 1, 5 });
            p.AddCoeff(2.1329213141884900e+02, new int[] { 2, 1, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("10841864-2d21-4664-bb89-3e623435589d"));
            OrthonormalPolynomials[243] = p;
            p.AddCoeff(-4.9795112437429500e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.0456973611860200e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(-3.1370920835580600e+01, new int[] { 0, 0, 4 });
            p.AddCoeff(2.3005341946092400e+01, new int[] { 0, 0, 6 });
            p.AddCoeff(1.4938533731228800e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-3.1370920835580600e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(9.4112762506741700e+01, new int[] { 0, 2, 4 });
            p.AddCoeff(-6.9016025838277300e+01, new int[] { 0, 2, 6 });
            p.AddCoeff(1.4938533731228800e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-3.1370920835580600e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(9.4112762506741700e+01, new int[] { 2, 0, 4 });
            p.AddCoeff(-6.9016025838277300e+01, new int[] { 2, 0, 6 });
            p.AddCoeff(-4.4815601193686500e+00, new int[] { 2, 2, 0 });
            p.AddCoeff(9.4112762506741700e+01, new int[] { 2, 2, 2 });
            p.AddCoeff(-2.8233828752022500e+02, new int[] { 2, 2, 4 });
            p.AddCoeff(2.0704807751483200e+02, new int[] { 2, 2, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5821a0b4-c7cb-45e2-85c7-628dc8009b7d"));
            OrthonormalPolynomials[244] = p;
            p.AddCoeff(9.7554634632503400e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-4.5525496161834900e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(4.0972946545651400e+01, new int[] { 0, 1, 5 });
            p.AddCoeff(-1.6259105772083900e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(7.5875826936391600e+01, new int[] { 0, 3, 3 });
            p.AddCoeff(-6.8288244242752400e+01, new int[] { 0, 3, 5 });
            p.AddCoeff(-2.9266390389751000e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(1.3657648848550500e+02, new int[] { 2, 1, 3 });
            p.AddCoeff(-1.2291883963695400e+02, new int[] { 2, 1, 5 });
            p.AddCoeff(4.8777317316251700e+01, new int[] { 2, 3, 1 });
            p.AddCoeff(-2.2762748080917500e+02, new int[] { 2, 3, 3 });
            p.AddCoeff(2.0486473272825700e+02, new int[] { 2, 3, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("0b566338-86a3-4d17-864f-ac9b26de02e5"));
            OrthonormalPolynomials[245] = p;
            p.AddCoeff(-5.0028220795632600e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(5.0028220795632600e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-5.8366257594904700e+00, new int[] { 0, 0, 4 });
            p.AddCoeff(5.0028220795632600e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-5.0028220795632600e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(5.8366257594904700e+01, new int[] { 0, 2, 4 });
            p.AddCoeff(-5.8366257594904700e+00, new int[] { 0, 4, 0 });
            p.AddCoeff(5.8366257594904700e+01, new int[] { 0, 4, 2 });
            p.AddCoeff(-6.8093967194055400e+01, new int[] { 0, 4, 4 });
            p.AddCoeff(1.5008466238689800e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-1.5008466238689800e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(1.7509877278471400e+01, new int[] { 2, 0, 4 });
            p.AddCoeff(-1.5008466238689800e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(1.5008466238689800e+02, new int[] { 2, 2, 2 });
            p.AddCoeff(-1.7509877278471400e+02, new int[] { 2, 2, 4 });
            p.AddCoeff(1.7509877278471400e+01, new int[] { 2, 4, 0 });
            p.AddCoeff(-1.7509877278471400e+02, new int[] { 2, 4, 2 });
            p.AddCoeff(2.0428190158216600e+02, new int[] { 2, 4, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e86bcff8-1658-4e53-b899-f0469abad28b"));
            OrthonormalPolynomials[246] = p;
            p.AddCoeff(9.7554634632503400e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-1.6259105772083900e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(-4.5525496161834900e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(7.5875826936391600e+01, new int[] { 0, 3, 3 });
            p.AddCoeff(4.0972946545651400e+01, new int[] { 0, 5, 1 });
            p.AddCoeff(-6.8288244242752400e+01, new int[] { 0, 5, 3 });
            p.AddCoeff(-2.9266390389751000e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(4.8777317316251700e+01, new int[] { 2, 1, 3 });
            p.AddCoeff(1.3657648848550500e+02, new int[] { 2, 3, 1 });
            p.AddCoeff(-2.2762748080917500e+02, new int[] { 2, 3, 3 });
            p.AddCoeff(-1.2291883963695400e+02, new int[] { 2, 5, 1 });
            p.AddCoeff(2.0486473272825700e+02, new int[] { 2, 5, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("40e210f4-71c4-4eef-9dc3-7adb358919a8"));
            OrthonormalPolynomials[247] = p;
            p.AddCoeff(-4.9795112437429500e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.4938533731228800e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(1.0456973611860200e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(-3.1370920835580600e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(-3.1370920835580600e+01, new int[] { 0, 4, 0 });
            p.AddCoeff(9.4112762506741700e+01, new int[] { 0, 4, 2 });
            p.AddCoeff(2.3005341946092400e+01, new int[] { 0, 6, 0 });
            p.AddCoeff(-6.9016025838277300e+01, new int[] { 0, 6, 2 });
            p.AddCoeff(1.4938533731228800e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-4.4815601193686500e+00, new int[] { 2, 0, 2 });
            p.AddCoeff(-3.1370920835580600e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(9.4112762506741700e+01, new int[] { 2, 2, 2 });
            p.AddCoeff(9.4112762506741700e+01, new int[] { 2, 4, 0 });
            p.AddCoeff(-2.8233828752022500e+02, new int[] { 2, 4, 2 });
            p.AddCoeff(-6.9016025838277300e+01, new int[] { 2, 6, 0 });
            p.AddCoeff(2.0704807751483200e+02, new int[] { 2, 6, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("2fa5233d-0536-4cba-8d9f-3aab7d9a0360"));
            OrthonormalPolynomials[248] = p;
            p.AddCoeff(5.8004853144209200e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-5.2204367829788300e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(1.1484960922553400e+02, new int[] { 0, 5, 1 });
            p.AddCoeff(-7.1097377139616400e+01, new int[] { 0, 7, 1 });
            p.AddCoeff(-1.7401455943262800e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(1.5661310348936500e+02, new int[] { 2, 3, 1 });
            p.AddCoeff(-3.4454882767660300e+02, new int[] { 2, 5, 1 });
            p.AddCoeff(2.1329213141884900e+02, new int[] { 2, 7, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("0d315b13-f2a8-424e-b39f-462b08f6f443"));
            OrthonormalPolynomials[249] = p;
            p.AddCoeff(-4.4564860191815000e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.6043349669053400e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(-8.8238423179793700e+01, new int[] { 0, 4, 0 });
            p.AddCoeff(1.5294660017830900e+02, new int[] { 0, 6, 0 });
            p.AddCoeff(-8.1935678666951300e+01, new int[] { 0, 8, 0 });
            p.AddCoeff(1.3369458057544500e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-4.8130049007160200e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(2.6471526953938100e+02, new int[] { 2, 4, 0 });
            p.AddCoeff(-4.5883980053492700e+02, new int[] { 2, 6, 0 });
            p.AddCoeff(2.4580703600085400e+02, new int[] { 2, 8, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7617cdf2-8c09-48c1-9c0a-e1205319172e"));
            OrthonormalPolynomials[250] = p;
            p.AddCoeff(1.1887457487108500e+01, new int[] { 1, 0, 1 });
            p.AddCoeff(-1.0698711738397600e+02, new int[] { 1, 0, 3 });
            p.AddCoeff(2.3537165824474800e+02, new int[] { 1, 0, 5 });
            p.AddCoeff(-1.4570626462770100e+02, new int[] { 1, 0, 7 });
            p.AddCoeff(-1.9812429145180800e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(1.7831186230662700e+02, new int[] { 3, 0, 3 });
            p.AddCoeff(-3.9228609707458000e+02, new int[] { 3, 0, 5 });
            p.AddCoeff(2.4284377437950200e+02, new int[] { 3, 0, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("22f5e4c9-20c0-4aae-923b-f606895fbecc"));
            OrthonormalPolynomials[251] = p;
            p.AddCoeff(2.7382784638002900e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-5.7503847739806100e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(1.7251154321941800e+02, new int[] { 1, 1, 4 });
            p.AddCoeff(-1.2650846502757300e+02, new int[] { 1, 1, 6 });
            p.AddCoeff(-4.5637974396671500e+00, new int[] { 3, 1, 0 });
            p.AddCoeff(9.5839746233010100e+01, new int[] { 3, 1, 2 });
            p.AddCoeff(-2.8751923869903000e+02, new int[] { 3, 1, 4 });
            p.AddCoeff(2.1084744171262200e+02, new int[] { 3, 1, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d61af773-36bf-4c4e-83ce-f3e55b4981de"));
            OrthonormalPolynomials[252] = p;
            p.AddCoeff(9.7554634632503400e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-4.5525496161834900e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(4.0972946545651400e+01, new int[] { 1, 0, 5 });
            p.AddCoeff(-2.9266390389751000e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(1.3657648848550500e+02, new int[] { 1, 2, 3 });
            p.AddCoeff(-1.2291883963695400e+02, new int[] { 1, 2, 5 });
            p.AddCoeff(-1.6259105772083900e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(7.5875826936391600e+01, new int[] { 3, 0, 3 });
            p.AddCoeff(-6.8288244242752400e+01, new int[] { 3, 0, 5 });
            p.AddCoeff(4.8777317316251700e+01, new int[] { 3, 2, 1 });
            p.AddCoeff(-2.2762748080917500e+02, new int[] { 3, 2, 3 });
            p.AddCoeff(2.0486473272825700e+02, new int[] { 3, 2, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e6dfecc1-8ffe-4723-b80d-cf4afdcc0ae5"));
            OrthonormalPolynomials[253] = p;
            p.AddCoeff(6.2645241395745900e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-6.2645241395745900e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(7.3086114961703600e+01, new int[] { 1, 1, 4 });
            p.AddCoeff(-1.0440873565957700e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(1.0440873565957700e+02, new int[] { 1, 3, 2 });
            p.AddCoeff(-1.2181019160283900e+02, new int[] { 1, 3, 4 });
            p.AddCoeff(-1.0440873565957700e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(1.0440873565957700e+02, new int[] { 3, 1, 2 });
            p.AddCoeff(-1.2181019160283900e+02, new int[] { 3, 1, 4 });
            p.AddCoeff(1.7401455943262800e+01, new int[] { 3, 3, 0 });
            p.AddCoeff(-1.7401455943262800e+02, new int[] { 3, 3, 2 });
            p.AddCoeff(2.0301698600473200e+02, new int[] { 3, 3, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("29d1eef9-4698-4f20-8c03-2c2b610f459f"));
            OrthonormalPolynomials[254] = p;
            p.AddCoeff(6.2645241395745900e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-1.0440873565957700e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(-6.2645241395745900e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(1.0440873565957700e+02, new int[] { 1, 2, 3 });
            p.AddCoeff(7.3086114961703600e+01, new int[] { 1, 4, 1 });
            p.AddCoeff(-1.2181019160283900e+02, new int[] { 1, 4, 3 });
            p.AddCoeff(-1.0440873565957700e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(1.7401455943262800e+01, new int[] { 3, 0, 3 });
            p.AddCoeff(1.0440873565957700e+02, new int[] { 3, 2, 1 });
            p.AddCoeff(-1.7401455943262800e+02, new int[] { 3, 2, 3 });
            p.AddCoeff(-1.2181019160283900e+02, new int[] { 3, 4, 1 });
            p.AddCoeff(2.0301698600473200e+02, new int[] { 3, 4, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("4d9cc007-28fa-4f9b-b979-5f1d14a734df"));
            OrthonormalPolynomials[255] = p;
            p.AddCoeff(9.7554634632503400e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-2.9266390389751000e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(-4.5525496161834900e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(1.3657648848550500e+02, new int[] { 1, 3, 2 });
            p.AddCoeff(4.0972946545651400e+01, new int[] { 1, 5, 0 });
            p.AddCoeff(-1.2291883963695400e+02, new int[] { 1, 5, 2 });
            p.AddCoeff(-1.6259105772083900e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(4.8777317316251700e+01, new int[] { 3, 1, 2 });
            p.AddCoeff(7.5875826936391600e+01, new int[] { 3, 3, 0 });
            p.AddCoeff(-2.2762748080917500e+02, new int[] { 3, 3, 2 });
            p.AddCoeff(-6.8288244242752400e+01, new int[] { 3, 5, 0 });
            p.AddCoeff(2.0486473272825700e+02, new int[] { 3, 5, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f61966b0-e708-4101-bc35-970a55d820b6"));
            OrthonormalPolynomials[256] = p;
            p.AddCoeff(2.7382784638002900e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-5.7503847739806100e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(1.7251154321941800e+02, new int[] { 1, 4, 1 });
            p.AddCoeff(-1.2650846502757300e+02, new int[] { 1, 6, 1 });
            p.AddCoeff(-4.5637974396671500e+00, new int[] { 3, 0, 1 });
            p.AddCoeff(9.5839746233010100e+01, new int[] { 3, 2, 1 });
            p.AddCoeff(-2.8751923869903000e+02, new int[] { 3, 4, 1 });
            p.AddCoeff(2.1084744171262200e+02, new int[] { 3, 6, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("ea0a324d-4377-4727-a12c-2dafe1b710d2"));
            OrthonormalPolynomials[257] = p;
            p.AddCoeff(1.1887457487108500e+01, new int[] { 1, 1, 0 });
            p.AddCoeff(-1.0698711738397600e+02, new int[] { 1, 3, 0 });
            p.AddCoeff(2.3537165824474800e+02, new int[] { 1, 5, 0 });
            p.AddCoeff(-1.4570626462770100e+02, new int[] { 1, 7, 0 });
            p.AddCoeff(-1.9812429145180800e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(1.7831186230662700e+02, new int[] { 3, 3, 0 });
            p.AddCoeff(-3.9228609707458000e+02, new int[] { 3, 5, 0 });
            p.AddCoeff(2.4284377437950200e+02, new int[] { 3, 7, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3eefa7a3-6591-4725-a191-283302b99585"));
            OrthonormalPolynomials[258] = p;
            p.AddCoeff(-4.4815601193686500e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(9.4112762506741700e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-2.8233828752022500e+01, new int[] { 0, 0, 4 });
            p.AddCoeff(2.0704807751483200e+01, new int[] { 0, 0, 6 });
            p.AddCoeff(4.4815601193686500e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-9.4112762506741700e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(2.8233828752022500e+02, new int[] { 2, 0, 4 });
            p.AddCoeff(-2.0704807751483200e+02, new int[] { 2, 0, 6 });
            p.AddCoeff(-5.2284868059301000e+00, new int[] { 4, 0, 0 });
            p.AddCoeff(1.0979822292453200e+02, new int[] { 4, 0, 2 });
            p.AddCoeff(-3.2939466877359600e+02, new int[] { 4, 0, 4 });
            p.AddCoeff(2.4155609043397000e+02, new int[] { 4, 0, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("49f3c8a3-5d57-45b4-805e-6440ea6a1556"));
            OrthonormalPolynomials[259] = p;
            p.AddCoeff(4.2841608774447400e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-1.9992750761408800e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(1.7993475685267900e+01, new int[] { 0, 1, 5 });
            p.AddCoeff(-4.2841608774447400e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(1.9992750761408800e+02, new int[] { 2, 1, 3 });
            p.AddCoeff(-1.7993475685267900e+02, new int[] { 2, 1, 5 });
            p.AddCoeff(4.9981876903522000e+01, new int[] { 4, 1, 1 });
            p.AddCoeff(-2.3324875888310300e+02, new int[] { 4, 1, 3 });
            p.AddCoeff(2.0992388299479200e+02, new int[] { 4, 1, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("6dd03eb6-e58c-4d0b-a3e8-bffd4636abd7"));
            OrthonormalPolynomials[260] = p;
            p.AddCoeff(-5.0028220795632600e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(5.0028220795632600e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-5.8366257594904700e+00, new int[] { 0, 0, 4 });
            p.AddCoeff(1.5008466238689800e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-1.5008466238689800e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(1.7509877278471400e+01, new int[] { 0, 2, 4 });
            p.AddCoeff(5.0028220795632600e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-5.0028220795632600e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(5.8366257594904700e+01, new int[] { 2, 0, 4 });
            p.AddCoeff(-1.5008466238689800e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(1.5008466238689800e+02, new int[] { 2, 2, 2 });
            p.AddCoeff(-1.7509877278471400e+02, new int[] { 2, 2, 4 });
            p.AddCoeff(-5.8366257594904700e+00, new int[] { 4, 0, 0 });
            p.AddCoeff(5.8366257594904700e+01, new int[] { 4, 0, 2 });
            p.AddCoeff(-6.8093967194055400e+01, new int[] { 4, 0, 4 });
            p.AddCoeff(1.7509877278471400e+01, new int[] { 4, 2, 0 });
            p.AddCoeff(-1.7509877278471400e+02, new int[] { 4, 2, 2 });
            p.AddCoeff(2.0428190158216600e+02, new int[] { 4, 2, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d7844f23-eb13-470f-bea2-3b80570f7087"));
            OrthonormalPolynomials[261] = p;
            p.AddCoeff(6.2645241395745900e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-1.0440873565957700e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(-1.0440873565957700e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(1.7401455943262800e+01, new int[] { 0, 3, 3 });
            p.AddCoeff(-6.2645241395745900e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(1.0440873565957700e+02, new int[] { 2, 1, 3 });
            p.AddCoeff(1.0440873565957700e+02, new int[] { 2, 3, 1 });
            p.AddCoeff(-1.7401455943262800e+02, new int[] { 2, 3, 3 });
            p.AddCoeff(7.3086114961703600e+01, new int[] { 4, 1, 1 });
            p.AddCoeff(-1.2181019160283900e+02, new int[] { 4, 1, 3 });
            p.AddCoeff(-1.2181019160283900e+02, new int[] { 4, 3, 1 });
            p.AddCoeff(2.0301698600473200e+02, new int[] { 4, 3, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("2da0a267-2d5f-4ebf-b738-1f10bbf28c6f"));
            OrthonormalPolynomials[262] = p;
            p.AddCoeff(-5.0028220795632600e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.5008466238689800e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(5.0028220795632600e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-1.5008466238689800e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(-5.8366257594904700e+00, new int[] { 0, 4, 0 });
            p.AddCoeff(1.7509877278471400e+01, new int[] { 0, 4, 2 });
            p.AddCoeff(5.0028220795632600e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-1.5008466238689800e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(-5.0028220795632600e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(1.5008466238689800e+02, new int[] { 2, 2, 2 });
            p.AddCoeff(5.8366257594904700e+01, new int[] { 2, 4, 0 });
            p.AddCoeff(-1.7509877278471400e+02, new int[] { 2, 4, 2 });
            p.AddCoeff(-5.8366257594904700e+00, new int[] { 4, 0, 0 });
            p.AddCoeff(1.7509877278471400e+01, new int[] { 4, 0, 2 });
            p.AddCoeff(5.8366257594904700e+01, new int[] { 4, 2, 0 });
            p.AddCoeff(-1.7509877278471400e+02, new int[] { 4, 2, 2 });
            p.AddCoeff(-6.8093967194055400e+01, new int[] { 4, 4, 0 });
            p.AddCoeff(2.0428190158216600e+02, new int[] { 4, 4, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5d6d27a5-3f01-4f66-8b3c-899f13895aa5"));
            OrthonormalPolynomials[263] = p;
            p.AddCoeff(4.2841608774447400e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-1.9992750761408800e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(1.7993475685267900e+01, new int[] { 0, 5, 1 });
            p.AddCoeff(-4.2841608774447400e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(1.9992750761408800e+02, new int[] { 2, 3, 1 });
            p.AddCoeff(-1.7993475685267900e+02, new int[] { 2, 5, 1 });
            p.AddCoeff(4.9981876903522000e+01, new int[] { 4, 1, 1 });
            p.AddCoeff(-2.3324875888310300e+02, new int[] { 4, 3, 1 });
            p.AddCoeff(2.0992388299479200e+02, new int[] { 4, 5, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3cb74f2b-d6a0-47c2-a18b-dbee5d85e373"));
            OrthonormalPolynomials[264] = p;
            p.AddCoeff(-4.4815601193686500e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(9.4112762506741700e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-2.8233828752022500e+01, new int[] { 0, 4, 0 });
            p.AddCoeff(2.0704807751483200e+01, new int[] { 0, 6, 0 });
            p.AddCoeff(4.4815601193686500e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-9.4112762506741700e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(2.8233828752022500e+02, new int[] { 2, 4, 0 });
            p.AddCoeff(-2.0704807751483200e+02, new int[] { 2, 6, 0 });
            p.AddCoeff(-5.2284868059301000e+00, new int[] { 4, 0, 0 });
            p.AddCoeff(1.0979822292453200e+02, new int[] { 4, 2, 0 });
            p.AddCoeff(-3.2939466877359600e+02, new int[] { 4, 4, 0 });
            p.AddCoeff(2.4155609043397000e+02, new int[] { 4, 6, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("4dfeeede-5141-4391-b511-0b8f134803c9"));
            OrthonormalPolynomials[265] = p;
            p.AddCoeff(1.3672572526849300e+01, new int[] { 1, 0, 1 });
            p.AddCoeff(-6.3805338458630100e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(5.7424804612767100e+01, new int[] { 1, 0, 5 });
            p.AddCoeff(-6.3805338458630100e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(2.9775824614027400e+02, new int[] { 3, 0, 3 });
            p.AddCoeff(-2.6798242152624600e+02, new int[] { 3, 0, 5 });
            p.AddCoeff(5.7424804612767100e+01, new int[] { 5, 0, 1 });
            p.AddCoeff(-2.6798242152624600e+02, new int[] { 5, 0, 3 });
            p.AddCoeff(2.4118417937362200e+02, new int[] { 5, 0, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("8d3e8d9f-2af3-4b2e-9cdc-b3f15e83362c"));
            OrthonormalPolynomials[266] = p;
            p.AddCoeff(4.2841608774447400e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-4.2841608774447400e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(4.9981876903522000e+01, new int[] { 1, 1, 4 });
            p.AddCoeff(-1.9992750761408800e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(1.9992750761408800e+02, new int[] { 3, 1, 2 });
            p.AddCoeff(-2.3324875888310300e+02, new int[] { 3, 1, 4 });
            p.AddCoeff(1.7993475685267900e+01, new int[] { 5, 1, 0 });
            p.AddCoeff(-1.7993475685267900e+02, new int[] { 5, 1, 2 });
            p.AddCoeff(2.0992388299479200e+02, new int[] { 5, 1, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("1f686dca-52b2-4b69-95c6-998ac47842cc"));
            OrthonormalPolynomials[267] = p;
            p.AddCoeff(9.7554634632503400e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-1.6259105772083900e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(-2.9266390389751000e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(4.8777317316251700e+01, new int[] { 1, 2, 3 });
            p.AddCoeff(-4.5525496161834900e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(7.5875826936391600e+01, new int[] { 3, 0, 3 });
            p.AddCoeff(1.3657648848550500e+02, new int[] { 3, 2, 1 });
            p.AddCoeff(-2.2762748080917500e+02, new int[] { 3, 2, 3 });
            p.AddCoeff(4.0972946545651400e+01, new int[] { 5, 0, 1 });
            p.AddCoeff(-6.8288244242752400e+01, new int[] { 5, 0, 3 });
            p.AddCoeff(-1.2291883963695400e+02, new int[] { 5, 2, 1 });
            p.AddCoeff(2.0486473272825700e+02, new int[] { 5, 2, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c7b14287-9387-471c-a1f4-13ed03e8acce"));
            OrthonormalPolynomials[268] = p;
            p.AddCoeff(9.7554634632503400e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-2.9266390389751000e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(-1.6259105772083900e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(4.8777317316251700e+01, new int[] { 1, 3, 2 });
            p.AddCoeff(-4.5525496161834900e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(1.3657648848550500e+02, new int[] { 3, 1, 2 });
            p.AddCoeff(7.5875826936391600e+01, new int[] { 3, 3, 0 });
            p.AddCoeff(-2.2762748080917500e+02, new int[] { 3, 3, 2 });
            p.AddCoeff(4.0972946545651400e+01, new int[] { 5, 1, 0 });
            p.AddCoeff(-1.2291883963695400e+02, new int[] { 5, 1, 2 });
            p.AddCoeff(-6.8288244242752400e+01, new int[] { 5, 3, 0 });
            p.AddCoeff(2.0486473272825700e+02, new int[] { 5, 3, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("2998cc95-8ae4-4b5f-9574-a979a7ec0489"));
            OrthonormalPolynomials[269] = p;
            p.AddCoeff(4.2841608774447400e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-4.2841608774447400e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(4.9981876903522000e+01, new int[] { 1, 4, 1 });
            p.AddCoeff(-1.9992750761408800e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(1.9992750761408800e+02, new int[] { 3, 2, 1 });
            p.AddCoeff(-2.3324875888310300e+02, new int[] { 3, 4, 1 });
            p.AddCoeff(1.7993475685267900e+01, new int[] { 5, 0, 1 });
            p.AddCoeff(-1.7993475685267900e+02, new int[] { 5, 2, 1 });
            p.AddCoeff(2.0992388299479200e+02, new int[] { 5, 4, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("49242672-6da4-42b3-bed7-727ed25acfb1"));
            OrthonormalPolynomials[270] = p;
            p.AddCoeff(1.3672572526849300e+01, new int[] { 1, 1, 0 });
            p.AddCoeff(-6.3805338458630100e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(5.7424804612767100e+01, new int[] { 1, 5, 0 });
            p.AddCoeff(-6.3805338458630100e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(2.9775824614027400e+02, new int[] { 3, 3, 0 });
            p.AddCoeff(-2.6798242152624600e+02, new int[] { 3, 5, 0 });
            p.AddCoeff(5.7424804612767100e+01, new int[] { 5, 1, 0 });
            p.AddCoeff(-2.6798242152624600e+02, new int[] { 5, 3, 0 });
            p.AddCoeff(2.4118417937362200e+02, new int[] { 5, 5, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("259a7226-73ed-487f-b20d-c7938568d7cb"));
            OrthonormalPolynomials[271] = p;
            p.AddCoeff(-4.4815601193686500e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(4.4815601193686500e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-5.2284868059301000e+00, new int[] { 0, 0, 4 });
            p.AddCoeff(9.4112762506741700e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-9.4112762506741700e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(1.0979822292453200e+02, new int[] { 2, 0, 4 });
            p.AddCoeff(-2.8233828752022500e+01, new int[] { 4, 0, 0 });
            p.AddCoeff(2.8233828752022500e+02, new int[] { 4, 0, 2 });
            p.AddCoeff(-3.2939466877359600e+02, new int[] { 4, 0, 4 });
            p.AddCoeff(2.0704807751483200e+01, new int[] { 6, 0, 0 });
            p.AddCoeff(-2.0704807751483200e+02, new int[] { 6, 0, 2 });
            p.AddCoeff(2.4155609043397000e+02, new int[] { 6, 0, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("1975862b-f221-47eb-9ca9-f37c8cc568fa"));
            OrthonormalPolynomials[272] = p;
            p.AddCoeff(2.7382784638002900e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-4.5637974396671500e+00, new int[] { 0, 1, 3 });
            p.AddCoeff(-5.7503847739806100e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(9.5839746233010100e+01, new int[] { 2, 1, 3 });
            p.AddCoeff(1.7251154321941800e+02, new int[] { 4, 1, 1 });
            p.AddCoeff(-2.8751923869903000e+02, new int[] { 4, 1, 3 });
            p.AddCoeff(-1.2650846502757300e+02, new int[] { 6, 1, 1 });
            p.AddCoeff(2.1084744171262200e+02, new int[] { 6, 1, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("4d67bcc9-b7a1-474b-9f71-1f7b40274ca1"));
            OrthonormalPolynomials[273] = p;
            p.AddCoeff(-4.9795112437429500e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.4938533731228800e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(1.4938533731228800e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-4.4815601193686500e+00, new int[] { 0, 2, 2 });
            p.AddCoeff(1.0456973611860200e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(-3.1370920835580600e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(-3.1370920835580600e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(9.4112762506741700e+01, new int[] { 2, 2, 2 });
            p.AddCoeff(-3.1370920835580600e+01, new int[] { 4, 0, 0 });
            p.AddCoeff(9.4112762506741700e+01, new int[] { 4, 0, 2 });
            p.AddCoeff(9.4112762506741700e+01, new int[] { 4, 2, 0 });
            p.AddCoeff(-2.8233828752022500e+02, new int[] { 4, 2, 2 });
            p.AddCoeff(2.3005341946092400e+01, new int[] { 6, 0, 0 });
            p.AddCoeff(-6.9016025838277300e+01, new int[] { 6, 0, 2 });
            p.AddCoeff(-6.9016025838277300e+01, new int[] { 6, 2, 0 });
            p.AddCoeff(2.0704807751483200e+02, new int[] { 6, 2, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("8580e037-e066-4c8d-bcd8-d7f6ea6ecadc"));
            OrthonormalPolynomials[274] = p;
            p.AddCoeff(2.7382784638002900e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-4.5637974396671500e+00, new int[] { 0, 3, 1 });
            p.AddCoeff(-5.7503847739806100e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(9.5839746233010100e+01, new int[] { 2, 3, 1 });
            p.AddCoeff(1.7251154321941800e+02, new int[] { 4, 1, 1 });
            p.AddCoeff(-2.8751923869903000e+02, new int[] { 4, 3, 1 });
            p.AddCoeff(-1.2650846502757300e+02, new int[] { 6, 1, 1 });
            p.AddCoeff(2.1084744171262200e+02, new int[] { 6, 3, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b4c4a604-e5fb-4dad-b5b4-f75c0eeb487c"));
            OrthonormalPolynomials[275] = p;
            p.AddCoeff(-4.4815601193686500e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(4.4815601193686500e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-5.2284868059301000e+00, new int[] { 0, 4, 0 });
            p.AddCoeff(9.4112762506741700e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-9.4112762506741700e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(1.0979822292453200e+02, new int[] { 2, 4, 0 });
            p.AddCoeff(-2.8233828752022500e+01, new int[] { 4, 0, 0 });
            p.AddCoeff(2.8233828752022500e+02, new int[] { 4, 2, 0 });
            p.AddCoeff(-3.2939466877359600e+02, new int[] { 4, 4, 0 });
            p.AddCoeff(2.0704807751483200e+01, new int[] { 6, 0, 0 });
            p.AddCoeff(-2.0704807751483200e+02, new int[] { 6, 2, 0 });
            p.AddCoeff(2.4155609043397000e+02, new int[] { 6, 4, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b185075a-3dc0-4aff-a72e-ab84ef12e4f7"));
            OrthonormalPolynomials[276] = p;
            p.AddCoeff(1.1887457487108500e+01, new int[] { 1, 0, 1 });
            p.AddCoeff(-1.9812429145180800e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(-1.0698711738397600e+02, new int[] { 3, 0, 1 });
            p.AddCoeff(1.7831186230662700e+02, new int[] { 3, 0, 3 });
            p.AddCoeff(2.3537165824474800e+02, new int[] { 5, 0, 1 });
            p.AddCoeff(-3.9228609707458000e+02, new int[] { 5, 0, 3 });
            p.AddCoeff(-1.4570626462770100e+02, new int[] { 7, 0, 1 });
            p.AddCoeff(2.4284377437950200e+02, new int[] { 7, 0, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d896ea3b-8ca1-4cf1-b68a-4e2c1986ef14"));
            OrthonormalPolynomials[277] = p;
            p.AddCoeff(5.8004853144209200e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-1.7401455943262800e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(-5.2204367829788300e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(1.5661310348936500e+02, new int[] { 3, 1, 2 });
            p.AddCoeff(1.1484960922553400e+02, new int[] { 5, 1, 0 });
            p.AddCoeff(-3.4454882767660300e+02, new int[] { 5, 1, 2 });
            p.AddCoeff(-7.1097377139616400e+01, new int[] { 7, 1, 0 });
            p.AddCoeff(2.1329213141884900e+02, new int[] { 7, 1, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("9e6f2914-fa94-4a8b-bdf7-dbfc80c95fe3"));
            OrthonormalPolynomials[278] = p;
            p.AddCoeff(5.8004853144209200e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-1.7401455943262800e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(-5.2204367829788300e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(1.5661310348936500e+02, new int[] { 3, 2, 1 });
            p.AddCoeff(1.1484960922553400e+02, new int[] { 5, 0, 1 });
            p.AddCoeff(-3.4454882767660300e+02, new int[] { 5, 2, 1 });
            p.AddCoeff(-7.1097377139616400e+01, new int[] { 7, 0, 1 });
            p.AddCoeff(2.1329213141884900e+02, new int[] { 7, 2, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7c5666ce-973d-4d32-8384-728985f018c1"));
            OrthonormalPolynomials[279] = p;
            p.AddCoeff(1.1887457487108500e+01, new int[] { 1, 1, 0 });
            p.AddCoeff(-1.9812429145180800e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(-1.0698711738397600e+02, new int[] { 3, 1, 0 });
            p.AddCoeff(1.7831186230662700e+02, new int[] { 3, 3, 0 });
            p.AddCoeff(2.3537165824474800e+02, new int[] { 5, 1, 0 });
            p.AddCoeff(-3.9228609707458000e+02, new int[] { 5, 3, 0 });
            p.AddCoeff(-1.4570626462770100e+02, new int[] { 7, 1, 0 });
            p.AddCoeff(2.4284377437950200e+02, new int[] { 7, 3, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("1558c50b-119e-4618-af59-d577329118bb"));
            OrthonormalPolynomials[280] = p;
            p.AddCoeff(-4.4564860191815000e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.3369458057544500e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(1.6043349669053400e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(-4.8130049007160200e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(-8.8238423179793700e+01, new int[] { 4, 0, 0 });
            p.AddCoeff(2.6471526953938100e+02, new int[] { 4, 0, 2 });
            p.AddCoeff(1.5294660017830900e+02, new int[] { 6, 0, 0 });
            p.AddCoeff(-4.5883980053492700e+02, new int[] { 6, 0, 2 });
            p.AddCoeff(-8.1935678666951300e+01, new int[] { 8, 0, 0 });
            p.AddCoeff(2.4580703600085400e+02, new int[] { 8, 0, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("62518d2e-020c-4a5f-90fb-399c662fd96f"));
            OrthonormalPolynomials[281] = p;
            p.AddCoeff(1.1958006815600700e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-4.3048824536162600e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(2.3676853494889400e+02, new int[] { 4, 1, 1 });
            p.AddCoeff(-4.1039879391141600e+02, new int[] { 6, 1, 1 });
            p.AddCoeff(2.1985649673825900e+02, new int[] { 8, 1, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("84d14657-a18a-4689-b5ff-39d8e8e63da1"));
            OrthonormalPolynomials[282] = p;
            p.AddCoeff(-4.4564860191815000e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.3369458057544500e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(1.6043349669053400e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(-4.8130049007160200e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(-8.8238423179793700e+01, new int[] { 4, 0, 0 });
            p.AddCoeff(2.6471526953938100e+02, new int[] { 4, 2, 0 });
            p.AddCoeff(1.5294660017830900e+02, new int[] { 6, 0, 0 });
            p.AddCoeff(-4.5883980053492700e+02, new int[] { 6, 2, 0 });
            p.AddCoeff(-8.1935678666951300e+01, new int[] { 8, 0, 0 });
            p.AddCoeff(2.4580703600085400e+02, new int[] { 8, 2, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("ad1b0202-4ccd-485e-97c2-8b171aaf11e2"));
            OrthonormalPolynomials[283] = p;
            p.AddCoeff(6.5689055652145700e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-9.6343948289813700e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(3.7574139833027300e+02, new int[] { 5, 0, 1 });
            p.AddCoeff(-5.3677342618610500e+02, new int[] { 7, 0, 1 });
            p.AddCoeff(2.5347634014343800e+02, new int[] { 9, 0, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("cac19797-153f-4c7b-a9a1-280af8e7eaf6"));
            OrthonormalPolynomials[284] = p;
            p.AddCoeff(6.5689055652145700e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-9.6343948289813700e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(3.7574139833027300e+02, new int[] { 5, 1, 0 });
            p.AddCoeff(-5.3677342618610500e+02, new int[] { 7, 1, 0 });
            p.AddCoeff(2.5347634014343800e+02, new int[] { 9, 1, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("718202d3-d101-48db-b841-1bdbaf95c062"));
            OrthonormalPolynomials[285] = p;
            p.AddCoeff(-3.9871744531220200e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(2.1929459492171100e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(-1.9005531559881600e+02, new int[] { 4, 0, 0 });
            p.AddCoeff(5.7016594679644900e+02, new int[] { 6, 0, 0 });
            p.AddCoeff(-6.9234436396711700e+02, new int[] { 8, 0, 0 });
            p.AddCoeff(2.9232317589722700e+02, new int[] { 10, 0, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("49bc864d-0a83-48a0-98bf-172c897a181d"));
            OrthonormalPolynomials[286] = p;
            p.AddCoeff(-4.5899948030330200e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(9.9449887399048800e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-5.9669932439429300e+02, new int[] { 0, 0, 5 });
            p.AddCoeff(1.4491269306718500e+03, new int[] { 0, 0, 7 });
            p.AddCoeff(-1.5296339823758500e+03, new int[] { 0, 0, 9 });
            p.AddCoeff(5.8404206599805000e+02, new int[] { 0, 0, 11 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e97e3072-6382-4f5c-903b-444767e130c2"));
            OrthonormalPolynomials[287] = p;
            p.AddCoeff(-6.9059887314479900e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(3.7982938022964000e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-3.2918546286568800e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(9.8755638859706300e+02, new int[] { 0, 1, 6 });
            p.AddCoeff(-1.1991756147250100e+03, new int[] { 0, 1, 8 });
            p.AddCoeff(5.0631859288389100e+02, new int[] { 0, 1, 10 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("1858254c-4beb-4bdd-b61e-6ce4c366f160"));
            OrthonormalPolynomials[288] = p;
            p.AddCoeff(-4.2402103094808700e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(6.2189751205719500e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-2.4254002970230600e+02, new int[] { 0, 0, 5 });
            p.AddCoeff(3.4648575671758000e+02, new int[] { 0, 0, 7 });
            p.AddCoeff(-1.6361827400552400e+02, new int[] { 0, 0, 9 });
            p.AddCoeff(1.2720630928442600e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-1.8656925361715800e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(7.2762008910691800e+02, new int[] { 0, 2, 5 });
            p.AddCoeff(-1.0394572701527400e+03, new int[] { 0, 2, 7 });
            p.AddCoeff(4.9085482201657100e+02, new int[] { 0, 2, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("85669bed-69c4-4557-9791-3cd015be4bd7"));
            OrthonormalPolynomials[289] = p;
            p.AddCoeff(-1.5818956105047400e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(5.6948241978170800e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-3.1321533087993900e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(5.4290657352522800e+02, new int[] { 0, 1, 6 });
            p.AddCoeff(-2.9084280724565800e+02, new int[] { 0, 1, 8 });
            p.AddCoeff(2.6364926841745700e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(-9.4913736630284600e+01, new int[] { 0, 3, 2 });
            p.AddCoeff(5.2202555146656500e+02, new int[] { 0, 3, 4 });
            p.AddCoeff(-9.0484428920871300e+02, new int[] { 0, 3, 6 });
            p.AddCoeff(4.8473801207609600e+02, new int[] { 0, 3, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7cc81647-9f5c-423b-af5c-6e3f5af330a5"));
            OrthonormalPolynomials[290] = p;
            p.AddCoeff(-3.3697774534009200e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(3.0327997080608300e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-6.6721593577338300e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(4.1303843643114200e+01, new int[] { 0, 0, 7 });
            p.AddCoeff(3.3697774534009200e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-3.0327997080608300e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(6.6721593577338300e+02, new int[] { 0, 2, 5 });
            p.AddCoeff(-4.1303843643114200e+02, new int[] { 0, 2, 7 });
            p.AddCoeff(-3.9314070289677400e+01, new int[] { 0, 4, 1 });
            p.AddCoeff(3.5382663260709700e+02, new int[] { 0, 4, 3 });
            p.AddCoeff(-7.7841859173561300e+02, new int[] { 0, 4, 5 });
            p.AddCoeff(4.8187817583633200e+02, new int[] { 0, 4, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a6b03dea-6c2c-4f35-81ae-dfce15cd6681"));
            OrthonormalPolynomials[291] = p;
            p.AddCoeff(-2.4772755652277100e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(5.2022786869781800e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-1.5606836060934500e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(1.1445013111352000e+02, new int[] { 0, 1, 6 });
            p.AddCoeff(1.1560619304396000e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-2.4277300539231500e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(7.2831901617694500e+02, new int[] { 0, 3, 4 });
            p.AddCoeff(-5.3410061186309300e+02, new int[] { 0, 3, 6 });
            p.AddCoeff(-1.0404557373956400e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(2.1849570485308400e+02, new int[] { 0, 5, 2 });
            p.AddCoeff(-6.5548711455925100e+02, new int[] { 0, 5, 4 });
            p.AddCoeff(4.8069055067678400e+02, new int[] { 0, 5, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("4ebb7e75-bba3-4695-bea9-d9c7e3486fa6"));
            OrthonormalPolynomials[292] = p;
            p.AddCoeff(-2.4772755652277100e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(1.1560619304396000e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-1.0404557373956400e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(5.2022786869781800e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-2.4277300539231500e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(2.1849570485308400e+02, new int[] { 0, 2, 5 });
            p.AddCoeff(-1.5606836060934500e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(7.2831901617694500e+02, new int[] { 0, 4, 3 });
            p.AddCoeff(-6.5548711455925100e+02, new int[] { 0, 4, 5 });
            p.AddCoeff(1.1445013111352000e+02, new int[] { 0, 6, 1 });
            p.AddCoeff(-5.3410061186309300e+02, new int[] { 0, 6, 3 });
            p.AddCoeff(4.8069055067678400e+02, new int[] { 0, 6, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("fc23948a-4bca-4a3c-abab-e261c0cf6df3"));
            OrthonormalPolynomials[293] = p;
            p.AddCoeff(-3.3697774534009200e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(3.3697774534009200e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-3.9314070289677400e+01, new int[] { 0, 1, 4 });
            p.AddCoeff(3.0327997080608300e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-3.0327997080608300e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(3.5382663260709700e+02, new int[] { 0, 3, 4 });
            p.AddCoeff(-6.6721593577338300e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(6.6721593577338300e+02, new int[] { 0, 5, 2 });
            p.AddCoeff(-7.7841859173561300e+02, new int[] { 0, 5, 4 });
            p.AddCoeff(4.1303843643114200e+01, new int[] { 0, 7, 0 });
            p.AddCoeff(-4.1303843643114200e+02, new int[] { 0, 7, 2 });
            p.AddCoeff(4.8187817583633200e+02, new int[] { 0, 7, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5f0bcadd-ca08-4deb-99ce-a6ab6a91bf80"));
            OrthonormalPolynomials[294] = p;
            p.AddCoeff(-1.5818956105047400e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(2.6364926841745700e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(5.6948241978170800e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-9.4913736630284600e+01, new int[] { 0, 2, 3 });
            p.AddCoeff(-3.1321533087993900e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(5.2202555146656500e+02, new int[] { 0, 4, 3 });
            p.AddCoeff(5.4290657352522800e+02, new int[] { 0, 6, 1 });
            p.AddCoeff(-9.0484428920871300e+02, new int[] { 0, 6, 3 });
            p.AddCoeff(-2.9084280724565800e+02, new int[] { 0, 8, 1 });
            p.AddCoeff(4.8473801207609600e+02, new int[] { 0, 8, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("69f57e27-7e3f-47d4-bc16-b161669bb034"));
            OrthonormalPolynomials[295] = p;
            p.AddCoeff(-4.2402103094808700e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(1.2720630928442600e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(6.2189751205719500e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-1.8656925361715800e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(-2.4254002970230600e+02, new int[] { 0, 5, 0 });
            p.AddCoeff(7.2762008910691800e+02, new int[] { 0, 5, 2 });
            p.AddCoeff(3.4648575671758000e+02, new int[] { 0, 7, 0 });
            p.AddCoeff(-1.0394572701527400e+03, new int[] { 0, 7, 2 });
            p.AddCoeff(-1.6361827400552400e+02, new int[] { 0, 9, 0 });
            p.AddCoeff(4.9085482201657100e+02, new int[] { 0, 9, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f0805f91-594e-4197-8dd6-8a6cb5a7d135"));
            OrthonormalPolynomials[296] = p;
            p.AddCoeff(-6.9059887314479900e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(3.7982938022964000e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-3.2918546286568800e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(9.8755638859706300e+02, new int[] { 0, 6, 1 });
            p.AddCoeff(-1.1991756147250100e+03, new int[] { 0, 8, 1 });
            p.AddCoeff(5.0631859288389100e+02, new int[] { 0, 10, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("8229af25-6bdd-4c13-b4b2-8ba50ee86fd5"));
            OrthonormalPolynomials[297] = p;
            p.AddCoeff(-4.5899948030330200e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(9.9449887399048800e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-5.9669932439429300e+02, new int[] { 0, 5, 0 });
            p.AddCoeff(1.4491269306718500e+03, new int[] { 0, 7, 0 });
            p.AddCoeff(-1.5296339823758500e+03, new int[] { 0, 9, 0 });
            p.AddCoeff(5.8404206599805000e+02, new int[] { 0, 11, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("ef14e701-f3df-4d04-a998-3e83dee08e60"));
            OrthonormalPolynomials[298] = p;
            p.AddCoeff(-6.9059887314479900e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(3.7982938022964000e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-3.2918546286568800e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(9.8755638859706300e+02, new int[] { 1, 0, 6 });
            p.AddCoeff(-1.1991756147250100e+03, new int[] { 1, 0, 8 });
            p.AddCoeff(5.0631859288389100e+02, new int[] { 1, 0, 10 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3e32503d-4727-4af1-8659-c94f8451565c"));
            OrthonormalPolynomials[299] = p;
            p.AddCoeff(1.1377678189073600e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-1.6687261343974600e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(6.5080319241501000e+02, new int[] { 1, 1, 5 });
            p.AddCoeff(-9.2971884630715700e+02, new int[] { 1, 1, 7 });
            p.AddCoeff(4.3903389964504700e+02, new int[] { 1, 1, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("49961184-cda6-4454-a01b-c1fa5423d93b"));
            OrthonormalPolynomials[300] = p;
            p.AddCoeff(-7.7188602084427100e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(2.7787896750393800e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-1.5283343212716600e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(2.6491128235375400e+02, new int[] { 1, 0, 6 });
            p.AddCoeff(-1.4191675840379700e+02, new int[] { 1, 0, 8 });
            p.AddCoeff(2.3156580625328100e+00, new int[] { 1, 2, 0 });
            p.AddCoeff(-8.3363690251181300e+01, new int[] { 1, 2, 2 });
            p.AddCoeff(4.5850029638149700e+02, new int[] { 1, 2, 4 });
            p.AddCoeff(-7.9473384706126200e+02, new int[] { 1, 2, 6 });
            p.AddCoeff(4.2575027521139000e+02, new int[] { 1, 2, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("9628a878-e7ba-4bff-926b-013462fe10ec"));
            OrthonormalPolynomials[301] = p;
            p.AddCoeff(2.0589680340487000e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-1.8530712306438300e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(4.0767567074164300e+02, new int[] { 1, 1, 5 });
            p.AddCoeff(-2.5237065331625500e+02, new int[] { 1, 1, 7 });
            p.AddCoeff(-3.4316133900811700e+01, new int[] { 1, 3, 1 });
            p.AddCoeff(3.0884520510730500e+02, new int[] { 1, 3, 3 });
            p.AddCoeff(-6.7945945123607200e+02, new int[] { 1, 3, 5 });
            p.AddCoeff(4.2061775552709200e+02, new int[] { 1, 3, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("0a06ece6-b798-446a-b2a3-5376860b8b55"));
            OrthonormalPolynomials[302] = p;
            p.AddCoeff(-7.7622898239209600e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(1.6300808630234000e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-4.8902425890702100e+01, new int[] { 1, 0, 4 });
            p.AddCoeff(3.5861778986514800e+01, new int[] { 1, 0, 6 });
            p.AddCoeff(7.7622898239209600e+00, new int[] { 1, 2, 0 });
            p.AddCoeff(-1.6300808630234000e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(4.8902425890702100e+02, new int[] { 1, 2, 4 });
            p.AddCoeff(-3.5861778986514800e+02, new int[] { 1, 2, 6 });
            p.AddCoeff(-9.0560047945744600e+00, new int[] { 1, 4, 0 });
            p.AddCoeff(1.9017610068606400e+02, new int[] { 1, 4, 2 });
            p.AddCoeff(-5.7052830205819100e+02, new int[] { 1, 4, 4 });
            p.AddCoeff(4.1838742150934000e+02, new int[] { 1, 4, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("57be6699-a7ee-4f2e-b01f-5d7bbb30884f"));
            OrthonormalPolynomials[303] = p;
            p.AddCoeff(2.3681590286673300e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-1.1051408800447600e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(9.9462679204028000e+01, new int[] { 1, 1, 5 });
            p.AddCoeff(-1.1051408800447600e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(5.1573241068755200e+02, new int[] { 1, 3, 3 });
            p.AddCoeff(-4.6415916961879700e+02, new int[] { 1, 3, 5 });
            p.AddCoeff(9.9462679204028000e+01, new int[] { 1, 5, 1 });
            p.AddCoeff(-4.6415916961879700e+02, new int[] { 1, 5, 3 });
            p.AddCoeff(4.1774325265691700e+02, new int[] { 1, 5, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("65fa6cfb-748e-4abb-acf0-f2dadd61735c"));
            OrthonormalPolynomials[304] = p;
            p.AddCoeff(-7.7622898239209600e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(7.7622898239209600e+00, new int[] { 1, 0, 2 });
            p.AddCoeff(-9.0560047945744600e+00, new int[] { 1, 0, 4 });
            p.AddCoeff(1.6300808630234000e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-1.6300808630234000e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(1.9017610068606400e+02, new int[] { 1, 2, 4 });
            p.AddCoeff(-4.8902425890702100e+01, new int[] { 1, 4, 0 });
            p.AddCoeff(4.8902425890702100e+02, new int[] { 1, 4, 2 });
            p.AddCoeff(-5.7052830205819100e+02, new int[] { 1, 4, 4 });
            p.AddCoeff(3.5861778986514800e+01, new int[] { 1, 6, 0 });
            p.AddCoeff(-3.5861778986514800e+02, new int[] { 1, 6, 2 });
            p.AddCoeff(4.1838742150934000e+02, new int[] { 1, 6, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("2683f5a5-7f68-45af-b540-7211d8a3550a"));
            OrthonormalPolynomials[305] = p;
            p.AddCoeff(2.0589680340487000e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-3.4316133900811700e+01, new int[] { 1, 1, 3 });
            p.AddCoeff(-1.8530712306438300e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(3.0884520510730500e+02, new int[] { 1, 3, 3 });
            p.AddCoeff(4.0767567074164300e+02, new int[] { 1, 5, 1 });
            p.AddCoeff(-6.7945945123607200e+02, new int[] { 1, 5, 3 });
            p.AddCoeff(-2.5237065331625500e+02, new int[] { 1, 7, 1 });
            p.AddCoeff(4.2061775552709200e+02, new int[] { 1, 7, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("23448438-2f74-4eb4-8f39-0b2480cc264a"));
            OrthonormalPolynomials[306] = p;
            p.AddCoeff(-7.7188602084427100e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(2.3156580625328100e+00, new int[] { 1, 0, 2 });
            p.AddCoeff(2.7787896750393800e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-8.3363690251181300e+01, new int[] { 1, 2, 2 });
            p.AddCoeff(-1.5283343212716600e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(4.5850029638149700e+02, new int[] { 1, 4, 2 });
            p.AddCoeff(2.6491128235375400e+02, new int[] { 1, 6, 0 });
            p.AddCoeff(-7.9473384706126200e+02, new int[] { 1, 6, 2 });
            p.AddCoeff(-1.4191675840379700e+02, new int[] { 1, 8, 0 });
            p.AddCoeff(4.2575027521139000e+02, new int[] { 1, 8, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("96322322-673e-4f9e-bab7-e596b28459be"));
            OrthonormalPolynomials[307] = p;
            p.AddCoeff(1.1377678189073600e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-1.6687261343974600e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(6.5080319241501000e+02, new int[] { 1, 5, 1 });
            p.AddCoeff(-9.2971884630715700e+02, new int[] { 1, 7, 1 });
            p.AddCoeff(4.3903389964504700e+02, new int[] { 1, 9, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5cfa45f5-15ce-49ca-97a3-4410d5f00c3e"));
            OrthonormalPolynomials[308] = p;
            p.AddCoeff(-6.9059887314479900e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(3.7982938022964000e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-3.2918546286568800e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(9.8755638859706300e+02, new int[] { 1, 6, 0 });
            p.AddCoeff(-1.1991756147250100e+03, new int[] { 1, 8, 0 });
            p.AddCoeff(5.0631859288389100e+02, new int[] { 1, 10, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("81147005-6c3f-4281-beb1-0968c693dbbf"));
            OrthonormalPolynomials[309] = p;
            p.AddCoeff(-4.2402103094808700e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(6.2189751205719500e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-2.4254002970230600e+02, new int[] { 0, 0, 5 });
            p.AddCoeff(3.4648575671758000e+02, new int[] { 0, 0, 7 });
            p.AddCoeff(-1.6361827400552400e+02, new int[] { 0, 0, 9 });
            p.AddCoeff(1.2720630928442600e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-1.8656925361715800e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(7.2762008910691800e+02, new int[] { 2, 0, 5 });
            p.AddCoeff(-1.0394572701527400e+03, new int[] { 2, 0, 7 });
            p.AddCoeff(4.9085482201657100e+02, new int[] { 2, 0, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("582d16f6-bfa0-403c-bdb9-809a6cedc54a"));
            OrthonormalPolynomials[310] = p;
            p.AddCoeff(-7.7188602084427100e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(2.7787896750393800e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-1.5283343212716600e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(2.6491128235375400e+02, new int[] { 0, 1, 6 });
            p.AddCoeff(-1.4191675840379700e+02, new int[] { 0, 1, 8 });
            p.AddCoeff(2.3156580625328100e+00, new int[] { 2, 1, 0 });
            p.AddCoeff(-8.3363690251181300e+01, new int[] { 2, 1, 2 });
            p.AddCoeff(4.5850029638149700e+02, new int[] { 2, 1, 4 });
            p.AddCoeff(-7.9473384706126200e+02, new int[] { 2, 1, 6 });
            p.AddCoeff(4.2575027521139000e+02, new int[] { 2, 1, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e2659833-2ea9-43cd-91a3-e83e9558dab8"));
            OrthonormalPolynomials[311] = p;
            p.AddCoeff(-3.7441971704454700e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(3.3697774534009200e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-7.4135103974820300e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(4.5893159603460200e+01, new int[] { 0, 0, 7 });
            p.AddCoeff(1.1232591511336400e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-1.0109332360202800e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(2.2240531192446100e+02, new int[] { 0, 2, 5 });
            p.AddCoeff(-1.3767947881038100e+02, new int[] { 0, 2, 7 });
            p.AddCoeff(1.1232591511336400e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-1.0109332360202800e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(2.2240531192446100e+02, new int[] { 2, 0, 5 });
            p.AddCoeff(-1.3767947881038100e+02, new int[] { 2, 0, 7 });
            p.AddCoeff(-3.3697774534009200e+01, new int[] { 2, 2, 1 });
            p.AddCoeff(3.0327997080608300e+02, new int[] { 2, 2, 3 });
            p.AddCoeff(-6.6721593577338300e+02, new int[] { 2, 2, 5 });
            p.AddCoeff(4.1303843643114200e+02, new int[] { 2, 2, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("243852e5-a620-42d1-b31e-e10ed14110c3"));
            OrthonormalPolynomials[312] = p;
            p.AddCoeff(-1.7675511479294900e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(3.7118574106519400e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-1.1135572231955800e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(8.1660863034342600e+01, new int[] { 0, 1, 6 });
            p.AddCoeff(2.9459185798824900e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(-6.1864290177532300e+01, new int[] { 0, 3, 2 });
            p.AddCoeff(1.8559287053259700e+02, new int[] { 0, 3, 4 });
            p.AddCoeff(-1.3610143839057100e+02, new int[] { 0, 3, 6 });
            p.AddCoeff(5.3026534437884800e+00, new int[] { 2, 1, 0 });
            p.AddCoeff(-1.1135572231955800e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(3.3406716695867400e+02, new int[] { 2, 1, 4 });
            p.AddCoeff(-2.4498258910302800e+02, new int[] { 2, 1, 6 });
            p.AddCoeff(-8.8377557396474700e+00, new int[] { 2, 3, 0 });
            p.AddCoeff(1.8559287053259700e+02, new int[] { 2, 3, 2 });
            p.AddCoeff(-5.5677861159779100e+02, new int[] { 2, 3, 4 });
            p.AddCoeff(4.0830431517171300e+02, new int[] { 2, 3, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("cb27880d-4aaf-499f-895f-9d74259869a7"));
            OrthonormalPolynomials[313] = p;
            p.AddCoeff(-2.7654139551361400e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(1.2905265123968600e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-1.1614738611571800e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(2.7654139551361400e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-1.2905265123968600e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(1.1614738611571800e+02, new int[] { 0, 2, 5 });
            p.AddCoeff(-3.2263162809921600e+01, new int[] { 0, 4, 1 });
            p.AddCoeff(1.5056142644630100e+02, new int[] { 0, 4, 3 });
            p.AddCoeff(-1.3550528380167100e+02, new int[] { 0, 4, 5 });
            p.AddCoeff(8.2962418654084100e+00, new int[] { 2, 0, 1 });
            p.AddCoeff(-3.8715795371905900e+01, new int[] { 2, 0, 3 });
            p.AddCoeff(3.4844215834715300e+01, new int[] { 2, 0, 5 });
            p.AddCoeff(-8.2962418654084100e+01, new int[] { 2, 2, 1 });
            p.AddCoeff(3.8715795371905900e+02, new int[] { 2, 2, 3 });
            p.AddCoeff(-3.4844215834715300e+02, new int[] { 2, 2, 5 });
            p.AddCoeff(9.6789488429764800e+01, new int[] { 2, 4, 1 });
            p.AddCoeff(-4.5168427933890200e+02, new int[] { 2, 4, 3 });
            p.AddCoeff(4.0651585140501200e+02, new int[] { 2, 4, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f35f80b7-a1d2-4e41-86db-3f201b195570"));
            OrthonormalPolynomials[314] = p;
            p.AddCoeff(-2.7654139551361400e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(2.7654139551361400e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-3.2263162809921600e+01, new int[] { 0, 1, 4 });
            p.AddCoeff(1.2905265123968600e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-1.2905265123968600e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(1.5056142644630100e+02, new int[] { 0, 3, 4 });
            p.AddCoeff(-1.1614738611571800e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(1.1614738611571800e+02, new int[] { 0, 5, 2 });
            p.AddCoeff(-1.3550528380167100e+02, new int[] { 0, 5, 4 });
            p.AddCoeff(8.2962418654084100e+00, new int[] { 2, 1, 0 });
            p.AddCoeff(-8.2962418654084100e+01, new int[] { 2, 1, 2 });
            p.AddCoeff(9.6789488429764800e+01, new int[] { 2, 1, 4 });
            p.AddCoeff(-3.8715795371905900e+01, new int[] { 2, 3, 0 });
            p.AddCoeff(3.8715795371905900e+02, new int[] { 2, 3, 2 });
            p.AddCoeff(-4.5168427933890200e+02, new int[] { 2, 3, 4 });
            p.AddCoeff(3.4844215834715300e+01, new int[] { 2, 5, 0 });
            p.AddCoeff(-3.4844215834715300e+02, new int[] { 2, 5, 2 });
            p.AddCoeff(4.0651585140501200e+02, new int[] { 2, 5, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("50f71174-38b7-4719-ac05-64ffd2459ca8"));
            OrthonormalPolynomials[315] = p;
            p.AddCoeff(-1.7675511479294900e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(2.9459185798824900e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(3.7118574106519400e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-6.1864290177532300e+01, new int[] { 0, 2, 3 });
            p.AddCoeff(-1.1135572231955800e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(1.8559287053259700e+02, new int[] { 0, 4, 3 });
            p.AddCoeff(8.1660863034342600e+01, new int[] { 0, 6, 1 });
            p.AddCoeff(-1.3610143839057100e+02, new int[] { 0, 6, 3 });
            p.AddCoeff(5.3026534437884800e+00, new int[] { 2, 0, 1 });
            p.AddCoeff(-8.8377557396474700e+00, new int[] { 2, 0, 3 });
            p.AddCoeff(-1.1135572231955800e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(1.8559287053259700e+02, new int[] { 2, 2, 3 });
            p.AddCoeff(3.3406716695867400e+02, new int[] { 2, 4, 1 });
            p.AddCoeff(-5.5677861159779100e+02, new int[] { 2, 4, 3 });
            p.AddCoeff(-2.4498258910302800e+02, new int[] { 2, 6, 1 });
            p.AddCoeff(4.0830431517171300e+02, new int[] { 2, 6, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f1bf2eb1-9ba6-45c5-afd0-a471e4ce6a0a"));
            OrthonormalPolynomials[316] = p;
            p.AddCoeff(-3.7441971704454700e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(1.1232591511336400e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(3.3697774534009200e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-1.0109332360202800e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(-7.4135103974820300e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(2.2240531192446100e+02, new int[] { 0, 5, 2 });
            p.AddCoeff(4.5893159603460200e+01, new int[] { 0, 7, 0 });
            p.AddCoeff(-1.3767947881038100e+02, new int[] { 0, 7, 2 });
            p.AddCoeff(1.1232591511336400e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-3.3697774534009200e+01, new int[] { 2, 1, 2 });
            p.AddCoeff(-1.0109332360202800e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(3.0327997080608300e+02, new int[] { 2, 3, 2 });
            p.AddCoeff(2.2240531192446100e+02, new int[] { 2, 5, 0 });
            p.AddCoeff(-6.6721593577338300e+02, new int[] { 2, 5, 2 });
            p.AddCoeff(-1.3767947881038100e+02, new int[] { 2, 7, 0 });
            p.AddCoeff(4.1303843643114200e+02, new int[] { 2, 7, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d91c2e88-9ae6-4872-8694-608693d430c2"));
            OrthonormalPolynomials[317] = p;
            p.AddCoeff(-7.7188602084427100e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(2.7787896750393800e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-1.5283343212716600e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(2.6491128235375400e+02, new int[] { 0, 6, 1 });
            p.AddCoeff(-1.4191675840379700e+02, new int[] { 0, 8, 1 });
            p.AddCoeff(2.3156580625328100e+00, new int[] { 2, 0, 1 });
            p.AddCoeff(-8.3363690251181300e+01, new int[] { 2, 2, 1 });
            p.AddCoeff(4.5850029638149700e+02, new int[] { 2, 4, 1 });
            p.AddCoeff(-7.9473384706126200e+02, new int[] { 2, 6, 1 });
            p.AddCoeff(4.2575027521139000e+02, new int[] { 2, 8, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("4acb2a31-8269-4f03-b958-52c3c601cb4f"));
            OrthonormalPolynomials[318] = p;
            p.AddCoeff(-4.2402103094808700e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(6.2189751205719500e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-2.4254002970230600e+02, new int[] { 0, 5, 0 });
            p.AddCoeff(3.4648575671758000e+02, new int[] { 0, 7, 0 });
            p.AddCoeff(-1.6361827400552400e+02, new int[] { 0, 9, 0 });
            p.AddCoeff(1.2720630928442600e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-1.8656925361715800e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(7.2762008910691800e+02, new int[] { 2, 5, 0 });
            p.AddCoeff(-1.0394572701527400e+03, new int[] { 2, 7, 0 });
            p.AddCoeff(4.9085482201657100e+02, new int[] { 2, 9, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5c1ebbef-2ba6-42be-8afa-5a452a70dcf2"));
            OrthonormalPolynomials[319] = p;
            p.AddCoeff(-1.5818956105047400e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(5.6948241978170800e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-3.1321533087993900e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(5.4290657352522800e+02, new int[] { 1, 0, 6 });
            p.AddCoeff(-2.9084280724565800e+02, new int[] { 1, 0, 8 });
            p.AddCoeff(2.6364926841745700e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(-9.4913736630284600e+01, new int[] { 3, 0, 2 });
            p.AddCoeff(5.2202555146656500e+02, new int[] { 3, 0, 4 });
            p.AddCoeff(-9.0484428920871300e+02, new int[] { 3, 0, 6 });
            p.AddCoeff(4.8473801207609600e+02, new int[] { 3, 0, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3c0d1fc2-d5ee-40c2-b1e5-aff29cd8c98c"));
            OrthonormalPolynomials[320] = p;
            p.AddCoeff(2.0589680340487000e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-1.8530712306438300e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(4.0767567074164300e+02, new int[] { 1, 1, 5 });
            p.AddCoeff(-2.5237065331625500e+02, new int[] { 1, 1, 7 });
            p.AddCoeff(-3.4316133900811700e+01, new int[] { 3, 1, 1 });
            p.AddCoeff(3.0884520510730500e+02, new int[] { 3, 1, 3 });
            p.AddCoeff(-6.7945945123607200e+02, new int[] { 3, 1, 5 });
            p.AddCoeff(4.2061775552709200e+02, new int[] { 3, 1, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("135dc339-d95b-4cd0-b2e3-e4e65a508024"));
            OrthonormalPolynomials[321] = p;
            p.AddCoeff(-1.7675511479294900e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(3.7118574106519400e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-1.1135572231955800e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(8.1660863034342600e+01, new int[] { 1, 0, 6 });
            p.AddCoeff(5.3026534437884800e+00, new int[] { 1, 2, 0 });
            p.AddCoeff(-1.1135572231955800e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(3.3406716695867400e+02, new int[] { 1, 2, 4 });
            p.AddCoeff(-2.4498258910302800e+02, new int[] { 1, 2, 6 });
            p.AddCoeff(2.9459185798824900e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(-6.1864290177532300e+01, new int[] { 3, 0, 2 });
            p.AddCoeff(1.8559287053259700e+02, new int[] { 3, 0, 4 });
            p.AddCoeff(-1.3610143839057100e+02, new int[] { 3, 0, 6 });
            p.AddCoeff(-8.8377557396474700e+00, new int[] { 3, 2, 0 });
            p.AddCoeff(1.8559287053259700e+02, new int[] { 3, 2, 2 });
            p.AddCoeff(-5.5677861159779100e+02, new int[] { 3, 2, 4 });
            p.AddCoeff(4.0830431517171300e+02, new int[] { 3, 2, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a07033a8-1b04-4af9-9560-afe7e8ea3f44"));
            OrthonormalPolynomials[322] = p;
            p.AddCoeff(3.4628460101821400e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-1.6159948047516700e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(1.4543953242765000e+02, new int[] { 1, 1, 5 });
            p.AddCoeff(-5.7714100169702400e+01, new int[] { 1, 3, 1 });
            p.AddCoeff(2.6933246745861100e+02, new int[] { 1, 3, 3 });
            p.AddCoeff(-2.4239922071275000e+02, new int[] { 1, 3, 5 });
            p.AddCoeff(-5.7714100169702400e+01, new int[] { 3, 1, 1 });
            p.AddCoeff(2.6933246745861100e+02, new int[] { 3, 1, 3 });
            p.AddCoeff(-2.4239922071275000e+02, new int[] { 3, 1, 5 });
            p.AddCoeff(9.6190166949503900e+01, new int[] { 3, 3, 1 });
            p.AddCoeff(-4.4888744576435200e+02, new int[] { 3, 3, 3 });
            p.AddCoeff(4.0399870118791700e+02, new int[] { 3, 3, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("17b2a033-7b5e-44a2-8ce3-60c40c0ee324"));
            OrthonormalPolynomials[323] = p;
            p.AddCoeff(-1.7758256738009100e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(1.7758256738009100e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-2.0717966194344000e+01, new int[] { 1, 0, 4 });
            p.AddCoeff(1.7758256738009100e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-1.7758256738009100e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(2.0717966194344000e+02, new int[] { 1, 2, 4 });
            p.AddCoeff(-2.0717966194344000e+01, new int[] { 1, 4, 0 });
            p.AddCoeff(2.0717966194344000e+02, new int[] { 1, 4, 2 });
            p.AddCoeff(-2.4170960560068000e+02, new int[] { 1, 4, 4 });
            p.AddCoeff(2.9597094563348500e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(-2.9597094563348500e+01, new int[] { 3, 0, 2 });
            p.AddCoeff(3.4529943657240000e+01, new int[] { 3, 0, 4 });
            p.AddCoeff(-2.9597094563348500e+01, new int[] { 3, 2, 0 });
            p.AddCoeff(2.9597094563348500e+02, new int[] { 3, 2, 2 });
            p.AddCoeff(-3.4529943657240000e+02, new int[] { 3, 2, 4 });
            p.AddCoeff(3.4529943657240000e+01, new int[] { 3, 4, 0 });
            p.AddCoeff(-3.4529943657240000e+02, new int[] { 3, 4, 2 });
            p.AddCoeff(4.0284934266780000e+02, new int[] { 3, 4, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("074bd8d7-6222-45c1-96fc-8eceef79fbf1"));
            OrthonormalPolynomials[324] = p;
            p.AddCoeff(3.4628460101821400e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-5.7714100169702400e+01, new int[] { 1, 1, 3 });
            p.AddCoeff(-1.6159948047516700e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(2.6933246745861100e+02, new int[] { 1, 3, 3 });
            p.AddCoeff(1.4543953242765000e+02, new int[] { 1, 5, 1 });
            p.AddCoeff(-2.4239922071275000e+02, new int[] { 1, 5, 3 });
            p.AddCoeff(-5.7714100169702400e+01, new int[] { 3, 1, 1 });
            p.AddCoeff(9.6190166949503900e+01, new int[] { 3, 1, 3 });
            p.AddCoeff(2.6933246745861100e+02, new int[] { 3, 3, 1 });
            p.AddCoeff(-4.4888744576435200e+02, new int[] { 3, 3, 3 });
            p.AddCoeff(-2.4239922071275000e+02, new int[] { 3, 5, 1 });
            p.AddCoeff(4.0399870118791700e+02, new int[] { 3, 5, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("84c22f52-c419-4d6e-b12b-c5135930fcd0"));
            OrthonormalPolynomials[325] = p;
            p.AddCoeff(-1.7675511479294900e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(5.3026534437884800e+00, new int[] { 1, 0, 2 });
            p.AddCoeff(3.7118574106519400e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-1.1135572231955800e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(-1.1135572231955800e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(3.3406716695867400e+02, new int[] { 1, 4, 2 });
            p.AddCoeff(8.1660863034342600e+01, new int[] { 1, 6, 0 });
            p.AddCoeff(-2.4498258910302800e+02, new int[] { 1, 6, 2 });
            p.AddCoeff(2.9459185798824900e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(-8.8377557396474700e+00, new int[] { 3, 0, 2 });
            p.AddCoeff(-6.1864290177532300e+01, new int[] { 3, 2, 0 });
            p.AddCoeff(1.8559287053259700e+02, new int[] { 3, 2, 2 });
            p.AddCoeff(1.8559287053259700e+02, new int[] { 3, 4, 0 });
            p.AddCoeff(-5.5677861159779100e+02, new int[] { 3, 4, 2 });
            p.AddCoeff(-1.3610143839057100e+02, new int[] { 3, 6, 0 });
            p.AddCoeff(4.0830431517171300e+02, new int[] { 3, 6, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("ae1f8fe4-5c80-45bb-93fa-673202b80ae9"));
            OrthonormalPolynomials[326] = p;
            p.AddCoeff(2.0589680340487000e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-1.8530712306438300e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(4.0767567074164300e+02, new int[] { 1, 5, 1 });
            p.AddCoeff(-2.5237065331625500e+02, new int[] { 1, 7, 1 });
            p.AddCoeff(-3.4316133900811700e+01, new int[] { 3, 1, 1 });
            p.AddCoeff(3.0884520510730500e+02, new int[] { 3, 3, 1 });
            p.AddCoeff(-6.7945945123607200e+02, new int[] { 3, 5, 1 });
            p.AddCoeff(4.2061775552709200e+02, new int[] { 3, 7, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b283c2e9-ddfe-4201-a66f-852667177720"));
            OrthonormalPolynomials[327] = p;
            p.AddCoeff(-1.5818956105047400e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(5.6948241978170800e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-3.1321533087993900e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(5.4290657352522800e+02, new int[] { 1, 6, 0 });
            p.AddCoeff(-2.9084280724565800e+02, new int[] { 1, 8, 0 });
            p.AddCoeff(2.6364926841745700e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(-9.4913736630284600e+01, new int[] { 3, 2, 0 });
            p.AddCoeff(5.2202555146656500e+02, new int[] { 3, 4, 0 });
            p.AddCoeff(-9.0484428920871300e+02, new int[] { 3, 6, 0 });
            p.AddCoeff(4.8473801207609600e+02, new int[] { 3, 8, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("0b0dff78-7333-42de-ac25-dee050445f39"));
            OrthonormalPolynomials[328] = p;
            p.AddCoeff(-3.3697774534009200e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(3.0327997080608300e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-6.6721593577338300e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(4.1303843643114200e+01, new int[] { 0, 0, 7 });
            p.AddCoeff(3.3697774534009200e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-3.0327997080608300e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(6.6721593577338300e+02, new int[] { 2, 0, 5 });
            p.AddCoeff(-4.1303843643114200e+02, new int[] { 2, 0, 7 });
            p.AddCoeff(-3.9314070289677400e+01, new int[] { 4, 0, 1 });
            p.AddCoeff(3.5382663260709700e+02, new int[] { 4, 0, 3 });
            p.AddCoeff(-7.7841859173561300e+02, new int[] { 4, 0, 5 });
            p.AddCoeff(4.8187817583633200e+02, new int[] { 4, 0, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("20a268a6-75cf-41ce-923e-76fc1619bcca"));
            OrthonormalPolynomials[329] = p;
            p.AddCoeff(-7.7622898239209600e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(1.6300808630234000e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-4.8902425890702100e+01, new int[] { 0, 1, 4 });
            p.AddCoeff(3.5861778986514800e+01, new int[] { 0, 1, 6 });
            p.AddCoeff(7.7622898239209600e+00, new int[] { 2, 1, 0 });
            p.AddCoeff(-1.6300808630234000e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(4.8902425890702100e+02, new int[] { 2, 1, 4 });
            p.AddCoeff(-3.5861778986514800e+02, new int[] { 2, 1, 6 });
            p.AddCoeff(-9.0560047945744600e+00, new int[] { 4, 1, 0 });
            p.AddCoeff(1.9017610068606400e+02, new int[] { 4, 1, 2 });
            p.AddCoeff(-5.7052830205819100e+02, new int[] { 4, 1, 4 });
            p.AddCoeff(4.1838742150934000e+02, new int[] { 4, 1, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b40e1eb9-653f-467c-b794-c8e0ca7d46d6"));
            OrthonormalPolynomials[330] = p;
            p.AddCoeff(-2.7654139551361400e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(1.2905265123968600e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-1.1614738611571800e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(8.2962418654084100e+00, new int[] { 0, 2, 1 });
            p.AddCoeff(-3.8715795371905900e+01, new int[] { 0, 2, 3 });
            p.AddCoeff(3.4844215834715300e+01, new int[] { 0, 2, 5 });
            p.AddCoeff(2.7654139551361400e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-1.2905265123968600e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(1.1614738611571800e+02, new int[] { 2, 0, 5 });
            p.AddCoeff(-8.2962418654084100e+01, new int[] { 2, 2, 1 });
            p.AddCoeff(3.8715795371905900e+02, new int[] { 2, 2, 3 });
            p.AddCoeff(-3.4844215834715300e+02, new int[] { 2, 2, 5 });
            p.AddCoeff(-3.2263162809921600e+01, new int[] { 4, 0, 1 });
            p.AddCoeff(1.5056142644630100e+02, new int[] { 4, 0, 3 });
            p.AddCoeff(-1.3550528380167100e+02, new int[] { 4, 0, 5 });
            p.AddCoeff(9.6789488429764800e+01, new int[] { 4, 2, 1 });
            p.AddCoeff(-4.5168427933890200e+02, new int[] { 4, 2, 3 });
            p.AddCoeff(4.0651585140501200e+02, new int[] { 4, 2, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("16d3f789-694e-42b5-99e1-61a3714b838a"));
            OrthonormalPolynomials[331] = p;
            p.AddCoeff(-1.7758256738009100e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(1.7758256738009100e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-2.0717966194344000e+01, new int[] { 0, 1, 4 });
            p.AddCoeff(2.9597094563348500e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(-2.9597094563348500e+01, new int[] { 0, 3, 2 });
            p.AddCoeff(3.4529943657240000e+01, new int[] { 0, 3, 4 });
            p.AddCoeff(1.7758256738009100e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-1.7758256738009100e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(2.0717966194344000e+02, new int[] { 2, 1, 4 });
            p.AddCoeff(-2.9597094563348500e+01, new int[] { 2, 3, 0 });
            p.AddCoeff(2.9597094563348500e+02, new int[] { 2, 3, 2 });
            p.AddCoeff(-3.4529943657240000e+02, new int[] { 2, 3, 4 });
            p.AddCoeff(-2.0717966194344000e+01, new int[] { 4, 1, 0 });
            p.AddCoeff(2.0717966194344000e+02, new int[] { 4, 1, 2 });
            p.AddCoeff(-2.4170960560068000e+02, new int[] { 4, 1, 4 });
            p.AddCoeff(3.4529943657240000e+01, new int[] { 4, 3, 0 });
            p.AddCoeff(-3.4529943657240000e+02, new int[] { 4, 3, 2 });
            p.AddCoeff(4.0284934266780000e+02, new int[] { 4, 3, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d7490d4a-46ea-4c6a-be56-7943a08e5dc8"));
            OrthonormalPolynomials[332] = p;
            p.AddCoeff(-1.7758256738009100e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(2.9597094563348500e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(1.7758256738009100e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-2.9597094563348500e+01, new int[] { 0, 2, 3 });
            p.AddCoeff(-2.0717966194344000e+01, new int[] { 0, 4, 1 });
            p.AddCoeff(3.4529943657240000e+01, new int[] { 0, 4, 3 });
            p.AddCoeff(1.7758256738009100e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-2.9597094563348500e+01, new int[] { 2, 0, 3 });
            p.AddCoeff(-1.7758256738009100e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(2.9597094563348500e+02, new int[] { 2, 2, 3 });
            p.AddCoeff(2.0717966194344000e+02, new int[] { 2, 4, 1 });
            p.AddCoeff(-3.4529943657240000e+02, new int[] { 2, 4, 3 });
            p.AddCoeff(-2.0717966194344000e+01, new int[] { 4, 0, 1 });
            p.AddCoeff(3.4529943657240000e+01, new int[] { 4, 0, 3 });
            p.AddCoeff(2.0717966194344000e+02, new int[] { 4, 2, 1 });
            p.AddCoeff(-3.4529943657240000e+02, new int[] { 4, 2, 3 });
            p.AddCoeff(-2.4170960560068000e+02, new int[] { 4, 4, 1 });
            p.AddCoeff(4.0284934266780000e+02, new int[] { 4, 4, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("02103ec1-89af-4886-b6f4-60f6b3a07d5b"));
            OrthonormalPolynomials[333] = p;
            p.AddCoeff(-2.7654139551361400e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(8.2962418654084100e+00, new int[] { 0, 1, 2 });
            p.AddCoeff(1.2905265123968600e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-3.8715795371905900e+01, new int[] { 0, 3, 2 });
            p.AddCoeff(-1.1614738611571800e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(3.4844215834715300e+01, new int[] { 0, 5, 2 });
            p.AddCoeff(2.7654139551361400e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-8.2962418654084100e+01, new int[] { 2, 1, 2 });
            p.AddCoeff(-1.2905265123968600e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(3.8715795371905900e+02, new int[] { 2, 3, 2 });
            p.AddCoeff(1.1614738611571800e+02, new int[] { 2, 5, 0 });
            p.AddCoeff(-3.4844215834715300e+02, new int[] { 2, 5, 2 });
            p.AddCoeff(-3.2263162809921600e+01, new int[] { 4, 1, 0 });
            p.AddCoeff(9.6789488429764800e+01, new int[] { 4, 1, 2 });
            p.AddCoeff(1.5056142644630100e+02, new int[] { 4, 3, 0 });
            p.AddCoeff(-4.5168427933890200e+02, new int[] { 4, 3, 2 });
            p.AddCoeff(-1.3550528380167100e+02, new int[] { 4, 5, 0 });
            p.AddCoeff(4.0651585140501200e+02, new int[] { 4, 5, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7b52c93b-2967-4ed7-a8da-c9ca9b5a2025"));
            OrthonormalPolynomials[334] = p;
            p.AddCoeff(-7.7622898239209600e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(1.6300808630234000e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-4.8902425890702100e+01, new int[] { 0, 4, 1 });
            p.AddCoeff(3.5861778986514800e+01, new int[] { 0, 6, 1 });
            p.AddCoeff(7.7622898239209600e+00, new int[] { 2, 0, 1 });
            p.AddCoeff(-1.6300808630234000e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(4.8902425890702100e+02, new int[] { 2, 4, 1 });
            p.AddCoeff(-3.5861778986514800e+02, new int[] { 2, 6, 1 });
            p.AddCoeff(-9.0560047945744600e+00, new int[] { 4, 0, 1 });
            p.AddCoeff(1.9017610068606400e+02, new int[] { 4, 2, 1 });
            p.AddCoeff(-5.7052830205819100e+02, new int[] { 4, 4, 1 });
            p.AddCoeff(4.1838742150934000e+02, new int[] { 4, 6, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("34e4fe80-c576-4613-a9a5-434f3811feec"));
            OrthonormalPolynomials[335] = p;
            p.AddCoeff(-3.3697774534009200e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(3.0327997080608300e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-6.6721593577338300e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(4.1303843643114200e+01, new int[] { 0, 7, 0 });
            p.AddCoeff(3.3697774534009200e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-3.0327997080608300e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(6.6721593577338300e+02, new int[] { 2, 5, 0 });
            p.AddCoeff(-4.1303843643114200e+02, new int[] { 2, 7, 0 });
            p.AddCoeff(-3.9314070289677400e+01, new int[] { 4, 1, 0 });
            p.AddCoeff(3.5382663260709700e+02, new int[] { 4, 3, 0 });
            p.AddCoeff(-7.7841859173561300e+02, new int[] { 4, 5, 0 });
            p.AddCoeff(4.8187817583633200e+02, new int[] { 4, 7, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3df740ae-84c4-462f-92ef-deec4356d1e0"));
            OrthonormalPolynomials[336] = p;
            p.AddCoeff(-2.4772755652277100e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(5.2022786869781800e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-1.5606836060934500e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(1.1445013111352000e+02, new int[] { 1, 0, 6 });
            p.AddCoeff(1.1560619304396000e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-2.4277300539231500e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(7.2831901617694500e+02, new int[] { 3, 0, 4 });
            p.AddCoeff(-5.3410061186309300e+02, new int[] { 3, 0, 6 });
            p.AddCoeff(-1.0404557373956400e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(2.1849570485308400e+02, new int[] { 5, 0, 2 });
            p.AddCoeff(-6.5548711455925100e+02, new int[] { 5, 0, 4 });
            p.AddCoeff(4.8069055067678400e+02, new int[] { 5, 0, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("4be71f33-b2d4-4ff2-9eb9-6327e080d2f2"));
            OrthonormalPolynomials[337] = p;
            p.AddCoeff(2.3681590286673300e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-1.1051408800447600e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(9.9462679204028000e+01, new int[] { 1, 1, 5 });
            p.AddCoeff(-1.1051408800447600e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(5.1573241068755200e+02, new int[] { 3, 1, 3 });
            p.AddCoeff(-4.6415916961879700e+02, new int[] { 3, 1, 5 });
            p.AddCoeff(9.9462679204028000e+01, new int[] { 5, 1, 1 });
            p.AddCoeff(-4.6415916961879700e+02, new int[] { 5, 1, 3 });
            p.AddCoeff(4.1774325265691700e+02, new int[] { 5, 1, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("fdfc79c5-54cd-4ed9-b236-4e29bfa3002b"));
            OrthonormalPolynomials[338] = p;
            p.AddCoeff(-2.7654139551361400e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(2.7654139551361400e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-3.2263162809921600e+01, new int[] { 1, 0, 4 });
            p.AddCoeff(8.2962418654084100e+00, new int[] { 1, 2, 0 });
            p.AddCoeff(-8.2962418654084100e+01, new int[] { 1, 2, 2 });
            p.AddCoeff(9.6789488429764800e+01, new int[] { 1, 2, 4 });
            p.AddCoeff(1.2905265123968600e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-1.2905265123968600e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(1.5056142644630100e+02, new int[] { 3, 0, 4 });
            p.AddCoeff(-3.8715795371905900e+01, new int[] { 3, 2, 0 });
            p.AddCoeff(3.8715795371905900e+02, new int[] { 3, 2, 2 });
            p.AddCoeff(-4.5168427933890200e+02, new int[] { 3, 2, 4 });
            p.AddCoeff(-1.1614738611571800e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(1.1614738611571800e+02, new int[] { 5, 0, 2 });
            p.AddCoeff(-1.3550528380167100e+02, new int[] { 5, 0, 4 });
            p.AddCoeff(3.4844215834715300e+01, new int[] { 5, 2, 0 });
            p.AddCoeff(-3.4844215834715300e+02, new int[] { 5, 2, 2 });
            p.AddCoeff(4.0651585140501200e+02, new int[] { 5, 2, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("de2c4732-5c75-4c3d-aeeb-c3df2ea7946b"));
            OrthonormalPolynomials[339] = p;
            p.AddCoeff(3.4628460101821400e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-5.7714100169702400e+01, new int[] { 1, 1, 3 });
            p.AddCoeff(-5.7714100169702400e+01, new int[] { 1, 3, 1 });
            p.AddCoeff(9.6190166949503900e+01, new int[] { 1, 3, 3 });
            p.AddCoeff(-1.6159948047516700e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(2.6933246745861100e+02, new int[] { 3, 1, 3 });
            p.AddCoeff(2.6933246745861100e+02, new int[] { 3, 3, 1 });
            p.AddCoeff(-4.4888744576435200e+02, new int[] { 3, 3, 3 });
            p.AddCoeff(1.4543953242765000e+02, new int[] { 5, 1, 1 });
            p.AddCoeff(-2.4239922071275000e+02, new int[] { 5, 1, 3 });
            p.AddCoeff(-2.4239922071275000e+02, new int[] { 5, 3, 1 });
            p.AddCoeff(4.0399870118791700e+02, new int[] { 5, 3, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("88b8ace9-5ddc-40e4-84a4-97158e68e6d0"));
            OrthonormalPolynomials[340] = p;
            p.AddCoeff(-2.7654139551361400e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(8.2962418654084100e+00, new int[] { 1, 0, 2 });
            p.AddCoeff(2.7654139551361400e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-8.2962418654084100e+01, new int[] { 1, 2, 2 });
            p.AddCoeff(-3.2263162809921600e+01, new int[] { 1, 4, 0 });
            p.AddCoeff(9.6789488429764800e+01, new int[] { 1, 4, 2 });
            p.AddCoeff(1.2905265123968600e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-3.8715795371905900e+01, new int[] { 3, 0, 2 });
            p.AddCoeff(-1.2905265123968600e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(3.8715795371905900e+02, new int[] { 3, 2, 2 });
            p.AddCoeff(1.5056142644630100e+02, new int[] { 3, 4, 0 });
            p.AddCoeff(-4.5168427933890200e+02, new int[] { 3, 4, 2 });
            p.AddCoeff(-1.1614738611571800e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(3.4844215834715300e+01, new int[] { 5, 0, 2 });
            p.AddCoeff(1.1614738611571800e+02, new int[] { 5, 2, 0 });
            p.AddCoeff(-3.4844215834715300e+02, new int[] { 5, 2, 2 });
            p.AddCoeff(-1.3550528380167100e+02, new int[] { 5, 4, 0 });
            p.AddCoeff(4.0651585140501200e+02, new int[] { 5, 4, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("75e0711d-c44a-4d4d-8a3d-93a212eaf88a"));
            OrthonormalPolynomials[341] = p;
            p.AddCoeff(2.3681590286673300e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-1.1051408800447600e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(9.9462679204028000e+01, new int[] { 1, 5, 1 });
            p.AddCoeff(-1.1051408800447600e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(5.1573241068755200e+02, new int[] { 3, 3, 1 });
            p.AddCoeff(-4.6415916961879700e+02, new int[] { 3, 5, 1 });
            p.AddCoeff(9.9462679204028000e+01, new int[] { 5, 1, 1 });
            p.AddCoeff(-4.6415916961879700e+02, new int[] { 5, 3, 1 });
            p.AddCoeff(4.1774325265691700e+02, new int[] { 5, 5, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3995f636-4f7b-45ba-8cc1-a1723f9ea6ff"));
            OrthonormalPolynomials[342] = p;
            p.AddCoeff(-2.4772755652277100e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(5.2022786869781800e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-1.5606836060934500e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(1.1445013111352000e+02, new int[] { 1, 6, 0 });
            p.AddCoeff(1.1560619304396000e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-2.4277300539231500e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(7.2831901617694500e+02, new int[] { 3, 4, 0 });
            p.AddCoeff(-5.3410061186309300e+02, new int[] { 3, 6, 0 });
            p.AddCoeff(-1.0404557373956400e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(2.1849570485308400e+02, new int[] { 5, 2, 0 });
            p.AddCoeff(-6.5548711455925100e+02, new int[] { 5, 4, 0 });
            p.AddCoeff(4.8069055067678400e+02, new int[] { 5, 6, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("375ff7cb-ae25-4ac7-b7ed-ddd6ec735e6f"));
            OrthonormalPolynomials[343] = p;
            p.AddCoeff(-2.4772755652277100e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(1.1560619304396000e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-1.0404557373956400e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(5.2022786869781800e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-2.4277300539231500e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(2.1849570485308400e+02, new int[] { 2, 0, 5 });
            p.AddCoeff(-1.5606836060934500e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(7.2831901617694500e+02, new int[] { 4, 0, 3 });
            p.AddCoeff(-6.5548711455925100e+02, new int[] { 4, 0, 5 });
            p.AddCoeff(1.1445013111352000e+02, new int[] { 6, 0, 1 });
            p.AddCoeff(-5.3410061186309300e+02, new int[] { 6, 0, 3 });
            p.AddCoeff(4.8069055067678400e+02, new int[] { 6, 0, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("63a9de70-a437-4cbf-bfb1-819a3dabae16"));
            OrthonormalPolynomials[344] = p;
            p.AddCoeff(-7.7622898239209600e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(7.7622898239209600e+00, new int[] { 0, 1, 2 });
            p.AddCoeff(-9.0560047945744600e+00, new int[] { 0, 1, 4 });
            p.AddCoeff(1.6300808630234000e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-1.6300808630234000e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(1.9017610068606400e+02, new int[] { 2, 1, 4 });
            p.AddCoeff(-4.8902425890702100e+01, new int[] { 4, 1, 0 });
            p.AddCoeff(4.8902425890702100e+02, new int[] { 4, 1, 2 });
            p.AddCoeff(-5.7052830205819100e+02, new int[] { 4, 1, 4 });
            p.AddCoeff(3.5861778986514800e+01, new int[] { 6, 1, 0 });
            p.AddCoeff(-3.5861778986514800e+02, new int[] { 6, 1, 2 });
            p.AddCoeff(4.1838742150934000e+02, new int[] { 6, 1, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7dcc5e27-80ab-4223-a346-0d53abe6ba9e"));
            OrthonormalPolynomials[345] = p;
            p.AddCoeff(-1.7675511479294900e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(2.9459185798824900e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(5.3026534437884800e+00, new int[] { 0, 2, 1 });
            p.AddCoeff(-8.8377557396474700e+00, new int[] { 0, 2, 3 });
            p.AddCoeff(3.7118574106519400e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-6.1864290177532300e+01, new int[] { 2, 0, 3 });
            p.AddCoeff(-1.1135572231955800e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(1.8559287053259700e+02, new int[] { 2, 2, 3 });
            p.AddCoeff(-1.1135572231955800e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(1.8559287053259700e+02, new int[] { 4, 0, 3 });
            p.AddCoeff(3.3406716695867400e+02, new int[] { 4, 2, 1 });
            p.AddCoeff(-5.5677861159779100e+02, new int[] { 4, 2, 3 });
            p.AddCoeff(8.1660863034342600e+01, new int[] { 6, 0, 1 });
            p.AddCoeff(-1.3610143839057100e+02, new int[] { 6, 0, 3 });
            p.AddCoeff(-2.4498258910302800e+02, new int[] { 6, 2, 1 });
            p.AddCoeff(4.0830431517171300e+02, new int[] { 6, 2, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("419878e1-d95b-49cc-acca-370187605be4"));
            OrthonormalPolynomials[346] = p;
            p.AddCoeff(-1.7675511479294900e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(5.3026534437884800e+00, new int[] { 0, 1, 2 });
            p.AddCoeff(2.9459185798824900e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(-8.8377557396474700e+00, new int[] { 0, 3, 2 });
            p.AddCoeff(3.7118574106519400e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-1.1135572231955800e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(-6.1864290177532300e+01, new int[] { 2, 3, 0 });
            p.AddCoeff(1.8559287053259700e+02, new int[] { 2, 3, 2 });
            p.AddCoeff(-1.1135572231955800e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(3.3406716695867400e+02, new int[] { 4, 1, 2 });
            p.AddCoeff(1.8559287053259700e+02, new int[] { 4, 3, 0 });
            p.AddCoeff(-5.5677861159779100e+02, new int[] { 4, 3, 2 });
            p.AddCoeff(8.1660863034342600e+01, new int[] { 6, 1, 0 });
            p.AddCoeff(-2.4498258910302800e+02, new int[] { 6, 1, 2 });
            p.AddCoeff(-1.3610143839057100e+02, new int[] { 6, 3, 0 });
            p.AddCoeff(4.0830431517171300e+02, new int[] { 6, 3, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("44f71188-eadf-452b-be92-c5048af1fb4f"));
            OrthonormalPolynomials[347] = p;
            p.AddCoeff(-7.7622898239209600e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(7.7622898239209600e+00, new int[] { 0, 2, 1 });
            p.AddCoeff(-9.0560047945744600e+00, new int[] { 0, 4, 1 });
            p.AddCoeff(1.6300808630234000e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-1.6300808630234000e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(1.9017610068606400e+02, new int[] { 2, 4, 1 });
            p.AddCoeff(-4.8902425890702100e+01, new int[] { 4, 0, 1 });
            p.AddCoeff(4.8902425890702100e+02, new int[] { 4, 2, 1 });
            p.AddCoeff(-5.7052830205819100e+02, new int[] { 4, 4, 1 });
            p.AddCoeff(3.5861778986514800e+01, new int[] { 6, 0, 1 });
            p.AddCoeff(-3.5861778986514800e+02, new int[] { 6, 2, 1 });
            p.AddCoeff(4.1838742150934000e+02, new int[] { 6, 4, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a958ce57-3eb7-4144-bd17-8680f7f4a46d"));
            OrthonormalPolynomials[348] = p;
            p.AddCoeff(-2.4772755652277100e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(1.1560619304396000e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-1.0404557373956400e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(5.2022786869781800e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-2.4277300539231500e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(2.1849570485308400e+02, new int[] { 2, 5, 0 });
            p.AddCoeff(-1.5606836060934500e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(7.2831901617694500e+02, new int[] { 4, 3, 0 });
            p.AddCoeff(-6.5548711455925100e+02, new int[] { 4, 5, 0 });
            p.AddCoeff(1.1445013111352000e+02, new int[] { 6, 1, 0 });
            p.AddCoeff(-5.3410061186309300e+02, new int[] { 6, 3, 0 });
            p.AddCoeff(4.8069055067678400e+02, new int[] { 6, 5, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7a499ebd-6ca2-45db-9997-693070e8aeb5"));
            OrthonormalPolynomials[349] = p;
            p.AddCoeff(-3.3697774534009200e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(3.3697774534009200e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-3.9314070289677400e+01, new int[] { 1, 0, 4 });
            p.AddCoeff(3.0327997080608300e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-3.0327997080608300e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(3.5382663260709700e+02, new int[] { 3, 0, 4 });
            p.AddCoeff(-6.6721593577338300e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(6.6721593577338300e+02, new int[] { 5, 0, 2 });
            p.AddCoeff(-7.7841859173561300e+02, new int[] { 5, 0, 4 });
            p.AddCoeff(4.1303843643114200e+01, new int[] { 7, 0, 0 });
            p.AddCoeff(-4.1303843643114200e+02, new int[] { 7, 0, 2 });
            p.AddCoeff(4.8187817583633200e+02, new int[] { 7, 0, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("21bcb173-89c8-42bd-b00a-461edf0429df"));
            OrthonormalPolynomials[350] = p;
            p.AddCoeff(2.0589680340487000e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-3.4316133900811700e+01, new int[] { 1, 1, 3 });
            p.AddCoeff(-1.8530712306438300e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(3.0884520510730500e+02, new int[] { 3, 1, 3 });
            p.AddCoeff(4.0767567074164300e+02, new int[] { 5, 1, 1 });
            p.AddCoeff(-6.7945945123607200e+02, new int[] { 5, 1, 3 });
            p.AddCoeff(-2.5237065331625500e+02, new int[] { 7, 1, 1 });
            p.AddCoeff(4.2061775552709200e+02, new int[] { 7, 1, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("8fd40693-63ad-4c96-b7b7-669e8985c5ef"));
            OrthonormalPolynomials[351] = p;
            p.AddCoeff(-3.7441971704454700e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(1.1232591511336400e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(1.1232591511336400e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-3.3697774534009200e+01, new int[] { 1, 2, 2 });
            p.AddCoeff(3.3697774534009200e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-1.0109332360202800e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(-1.0109332360202800e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(3.0327997080608300e+02, new int[] { 3, 2, 2 });
            p.AddCoeff(-7.4135103974820300e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(2.2240531192446100e+02, new int[] { 5, 0, 2 });
            p.AddCoeff(2.2240531192446100e+02, new int[] { 5, 2, 0 });
            p.AddCoeff(-6.6721593577338300e+02, new int[] { 5, 2, 2 });
            p.AddCoeff(4.5893159603460200e+01, new int[] { 7, 0, 0 });
            p.AddCoeff(-1.3767947881038100e+02, new int[] { 7, 0, 2 });
            p.AddCoeff(-1.3767947881038100e+02, new int[] { 7, 2, 0 });
            p.AddCoeff(4.1303843643114200e+02, new int[] { 7, 2, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("cff5bd57-ffb7-45db-a1e4-83550aa06a1d"));
            OrthonormalPolynomials[352] = p;
            p.AddCoeff(2.0589680340487000e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-3.4316133900811700e+01, new int[] { 1, 3, 1 });
            p.AddCoeff(-1.8530712306438300e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(3.0884520510730500e+02, new int[] { 3, 3, 1 });
            p.AddCoeff(4.0767567074164300e+02, new int[] { 5, 1, 1 });
            p.AddCoeff(-6.7945945123607200e+02, new int[] { 5, 3, 1 });
            p.AddCoeff(-2.5237065331625500e+02, new int[] { 7, 1, 1 });
            p.AddCoeff(4.2061775552709200e+02, new int[] { 7, 3, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("9a163077-5308-44df-823c-3c333165adf9"));
            OrthonormalPolynomials[353] = p;
            p.AddCoeff(-3.3697774534009200e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(3.3697774534009200e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-3.9314070289677400e+01, new int[] { 1, 4, 0 });
            p.AddCoeff(3.0327997080608300e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-3.0327997080608300e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(3.5382663260709700e+02, new int[] { 3, 4, 0 });
            p.AddCoeff(-6.6721593577338300e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(6.6721593577338300e+02, new int[] { 5, 2, 0 });
            p.AddCoeff(-7.7841859173561300e+02, new int[] { 5, 4, 0 });
            p.AddCoeff(4.1303843643114200e+01, new int[] { 7, 0, 0 });
            p.AddCoeff(-4.1303843643114200e+02, new int[] { 7, 2, 0 });
            p.AddCoeff(4.8187817583633200e+02, new int[] { 7, 4, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("356ab5f4-8db1-4b87-914d-866e425f4248"));
            OrthonormalPolynomials[354] = p;
            p.AddCoeff(-1.5818956105047400e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(2.6364926841745700e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(5.6948241978170800e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-9.4913736630284600e+01, new int[] { 2, 0, 3 });
            p.AddCoeff(-3.1321533087993900e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(5.2202555146656500e+02, new int[] { 4, 0, 3 });
            p.AddCoeff(5.4290657352522800e+02, new int[] { 6, 0, 1 });
            p.AddCoeff(-9.0484428920871300e+02, new int[] { 6, 0, 3 });
            p.AddCoeff(-2.9084280724565800e+02, new int[] { 8, 0, 1 });
            p.AddCoeff(4.8473801207609600e+02, new int[] { 8, 0, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3b07a19f-aa22-4077-aa81-897186e8b1d4"));
            OrthonormalPolynomials[355] = p;
            p.AddCoeff(-7.7188602084427100e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(2.3156580625328100e+00, new int[] { 0, 1, 2 });
            p.AddCoeff(2.7787896750393800e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-8.3363690251181300e+01, new int[] { 2, 1, 2 });
            p.AddCoeff(-1.5283343212716600e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(4.5850029638149700e+02, new int[] { 4, 1, 2 });
            p.AddCoeff(2.6491128235375400e+02, new int[] { 6, 1, 0 });
            p.AddCoeff(-7.9473384706126200e+02, new int[] { 6, 1, 2 });
            p.AddCoeff(-1.4191675840379700e+02, new int[] { 8, 1, 0 });
            p.AddCoeff(4.2575027521139000e+02, new int[] { 8, 1, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e04c40cc-f782-4a9a-8c4d-84cfa548af75"));
            OrthonormalPolynomials[356] = p;
            p.AddCoeff(-7.7188602084427100e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(2.3156580625328100e+00, new int[] { 0, 2, 1 });
            p.AddCoeff(2.7787896750393800e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-8.3363690251181300e+01, new int[] { 2, 2, 1 });
            p.AddCoeff(-1.5283343212716600e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(4.5850029638149700e+02, new int[] { 4, 2, 1 });
            p.AddCoeff(2.6491128235375400e+02, new int[] { 6, 0, 1 });
            p.AddCoeff(-7.9473384706126200e+02, new int[] { 6, 2, 1 });
            p.AddCoeff(-1.4191675840379700e+02, new int[] { 8, 0, 1 });
            p.AddCoeff(4.2575027521139000e+02, new int[] { 8, 2, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("768972a2-200d-45e9-9dd2-f08755df3e16"));
            OrthonormalPolynomials[357] = p;
            p.AddCoeff(-1.5818956105047400e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(2.6364926841745700e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(5.6948241978170800e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-9.4913736630284600e+01, new int[] { 2, 3, 0 });
            p.AddCoeff(-3.1321533087993900e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(5.2202555146656500e+02, new int[] { 4, 3, 0 });
            p.AddCoeff(5.4290657352522800e+02, new int[] { 6, 1, 0 });
            p.AddCoeff(-9.0484428920871300e+02, new int[] { 6, 3, 0 });
            p.AddCoeff(-2.9084280724565800e+02, new int[] { 8, 1, 0 });
            p.AddCoeff(4.8473801207609600e+02, new int[] { 8, 3, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("eeb5a045-e2eb-4783-8a78-7e51d40d2b81"));
            OrthonormalPolynomials[358] = p;
            p.AddCoeff(-4.2402103094808700e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(1.2720630928442600e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(6.2189751205719500e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-1.8656925361715800e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(-2.4254002970230600e+02, new int[] { 5, 0, 0 });
            p.AddCoeff(7.2762008910691800e+02, new int[] { 5, 0, 2 });
            p.AddCoeff(3.4648575671758000e+02, new int[] { 7, 0, 0 });
            p.AddCoeff(-1.0394572701527400e+03, new int[] { 7, 0, 2 });
            p.AddCoeff(-1.6361827400552400e+02, new int[] { 9, 0, 0 });
            p.AddCoeff(4.9085482201657100e+02, new int[] { 9, 0, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e68ccac5-db3c-4be5-a217-2f50e6ff99e7"));
            OrthonormalPolynomials[359] = p;
            p.AddCoeff(1.1377678189073600e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-1.6687261343974600e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(6.5080319241501000e+02, new int[] { 5, 1, 1 });
            p.AddCoeff(-9.2971884630715700e+02, new int[] { 7, 1, 1 });
            p.AddCoeff(4.3903389964504700e+02, new int[] { 9, 1, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("6e9d0747-98ac-4539-a686-bb0e51ba57c8"));
            OrthonormalPolynomials[360] = p;
            p.AddCoeff(-4.2402103094808700e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(1.2720630928442600e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(6.2189751205719500e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-1.8656925361715800e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(-2.4254002970230600e+02, new int[] { 5, 0, 0 });
            p.AddCoeff(7.2762008910691800e+02, new int[] { 5, 2, 0 });
            p.AddCoeff(3.4648575671758000e+02, new int[] { 7, 0, 0 });
            p.AddCoeff(-1.0394572701527400e+03, new int[] { 7, 2, 0 });
            p.AddCoeff(-1.6361827400552400e+02, new int[] { 9, 0, 0 });
            p.AddCoeff(4.9085482201657100e+02, new int[] { 9, 2, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("19c22d53-74aa-428a-b534-d9b5d91bd515"));
            OrthonormalPolynomials[361] = p;
            p.AddCoeff(-6.9059887314479900e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(3.7982938022964000e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-3.2918546286568800e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(9.8755638859706300e+02, new int[] { 6, 0, 1 });
            p.AddCoeff(-1.1991756147250100e+03, new int[] { 8, 0, 1 });
            p.AddCoeff(5.0631859288389100e+02, new int[] { 10, 0, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("81b1c791-2d84-4895-beb8-9baaff39f169"));
            OrthonormalPolynomials[362] = p;
            p.AddCoeff(-6.9059887314479900e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(3.7982938022964000e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-3.2918546286568800e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(9.8755638859706300e+02, new int[] { 6, 1, 0 });
            p.AddCoeff(-1.1991756147250100e+03, new int[] { 8, 1, 0 });
            p.AddCoeff(5.0631859288389100e+02, new int[] { 10, 1, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("da42481f-0a9a-445e-882b-ce2fc09aa972"));
            OrthonormalPolynomials[363] = p;
            p.AddCoeff(-4.5899948030330200e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(9.9449887399048800e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-5.9669932439429300e+02, new int[] { 5, 0, 0 });
            p.AddCoeff(1.4491269306718500e+03, new int[] { 7, 0, 0 });
            p.AddCoeff(-1.5296339823758500e+03, new int[] { 9, 0, 0 });
            p.AddCoeff(5.8404206599805000e+02, new int[] { 11, 0, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("4dde8965-c72b-495e-a11a-4b62d896da17"));
            OrthonormalPolynomials[364] = p;
            p.AddCoeff(3.9878336536643800e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-3.1105102498582200e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(3.8881378123227700e+02, new int[] { 0, 0, 4 });
            p.AddCoeff(-1.7626224749196600e+03, new int[] { 0, 0, 6 });
            p.AddCoeff(3.5881957525150100e+03, new int[] { 0, 0, 8 });
            p.AddCoeff(-3.3489827023473500e+03, new int[] { 0, 0, 10 });
            p.AddCoeff(1.1670697296058900e+03, new int[] { 0, 0, 12 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("33e53dcc-1e8b-4b98-be80-f822f10a0d5a"));
            OrthonormalPolynomials[365] = p;
            p.AddCoeff(-7.9501042053302700e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(1.7225225778215600e+02, new int[] { 0, 1, 3 });
            p.AddCoeff(-1.0335135466929300e+03, new int[] { 0, 1, 5 });
            p.AddCoeff(2.5099614705399800e+03, new int[] { 0, 1, 7 });
            p.AddCoeff(-2.6494037744588700e+03, new int[] { 0, 1, 9 });
            p.AddCoeff(1.0115905320661200e+03, new int[] { 0, 1, 11 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("fa90b4e6-f44f-4cb8-949f-c6c38810d089"));
            OrthonormalPolynomials[366] = p;
            p.AddCoeff(4.4577965576656800e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-2.4517881067161300e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(2.1248830258206400e+02, new int[] { 0, 0, 4 });
            p.AddCoeff(-6.3746490774619300e+02, new int[] { 0, 0, 6 });
            p.AddCoeff(7.7406453083466300e+02, new int[] { 0, 0, 8 });
            p.AddCoeff(-3.2682724635241300e+02, new int[] { 0, 0, 10 });
            p.AddCoeff(-1.3373389672997100e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(7.3553643201483800e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(-6.3746490774619300e+02, new int[] { 0, 2, 4 });
            p.AddCoeff(1.9123947232385800e+03, new int[] { 0, 2, 6 });
            p.AddCoeff(-2.3221935925039900e+03, new int[] { 0, 2, 8 });
            p.AddCoeff(9.8048173905724000e+02, new int[] { 0, 2, 10 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c52472fa-6c3a-4be6-8c8c-a8f147544428"));
            OrthonormalPolynomials[367] = p;
            p.AddCoeff(-1.5051253492806200e+01, new int[] { 0, 1, 1 });
            p.AddCoeff(2.2075171789449200e+02, new int[] { 0, 1, 3 });
            p.AddCoeff(-8.6093169978851700e+02, new int[] { 0, 1, 5 });
            p.AddCoeff(1.2299024282693100e+03, new int[] { 0, 1, 7 });
            p.AddCoeff(-5.8078725779384100e+02, new int[] { 0, 1, 9 });
            p.AddCoeff(2.5085422488010400e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(-3.6791952982415300e+02, new int[] { 0, 3, 3 });
            p.AddCoeff(1.4348861663142000e+03, new int[] { 0, 3, 5 });
            p.AddCoeff(-2.0498373804488500e+03, new int[] { 0, 3, 7 });
            p.AddCoeff(9.6797876298973500e+02, new int[] { 0, 3, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e0d88e04-49e1-4b76-a6e0-7b1e024b0ef0"));
            OrthonormalPolynomials[368] = p;
            p.AddCoeff(4.4842525558502700e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-1.6143309201061000e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(8.8788200605835300e+01, new int[] { 0, 0, 4 });
            p.AddCoeff(-1.5389954771678100e+02, new int[] { 0, 0, 6 });
            p.AddCoeff(8.2446186276847100e+01, new int[] { 0, 0, 8 });
            p.AddCoeff(-4.4842525558502700e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(1.6143309201061000e+02, new int[] { 0, 2, 2 });
            p.AddCoeff(-8.8788200605835300e+02, new int[] { 0, 2, 4 });
            p.AddCoeff(1.5389954771678100e+03, new int[] { 0, 2, 6 });
            p.AddCoeff(-8.2446186276847100e+02, new int[] { 0, 2, 8 });
            p.AddCoeff(5.2316279818253100e+00, new int[] { 0, 4, 0 });
            p.AddCoeff(-1.8833860734571100e+02, new int[] { 0, 4, 2 });
            p.AddCoeff(1.0358623404014100e+03, new int[] { 0, 4, 4 });
            p.AddCoeff(-1.7954947233624500e+03, new int[] { 0, 4, 6 });
            p.AddCoeff(9.6187217322988200e+02, new int[] { 0, 4, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d5ee9599-559b-4f59-9392-86906192f77e"));
            OrthonormalPolynomials[369] = p;
            p.AddCoeff(-1.8627145733216900e+01, new int[] { 0, 1, 1 });
            p.AddCoeff(1.6764431159895200e+02, new int[] { 0, 1, 3 });
            p.AddCoeff(-3.6881748551769500e+02, new int[] { 0, 1, 5 });
            p.AddCoeff(2.2831558627285900e+02, new int[] { 0, 1, 7 });
            p.AddCoeff(8.6926680088345700e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(-7.8234012079511100e+02, new int[] { 0, 3, 3 });
            p.AddCoeff(1.7211482657492500e+03, new int[] { 0, 3, 5 });
            p.AddCoeff(-1.0654727359400100e+03, new int[] { 0, 3, 7 });
            p.AddCoeff(-7.8234012079511100e+01, new int[] { 0, 5, 1 });
            p.AddCoeff(7.0410610871560000e+02, new int[] { 0, 5, 3 });
            p.AddCoeff(-1.5490334391743200e+03, new int[] { 0, 5, 5 });
            p.AddCoeff(9.5892546234600800e+02, new int[] { 0, 5, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5d9cce32-6b90-4356-ab4d-c59130aaa3a4"));
            OrthonormalPolynomials[370] = p;
            p.AddCoeff(4.4884707790161900e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-9.4257886359339900e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(2.8277365907802000e+01, new int[] { 0, 0, 4 });
            p.AddCoeff(-2.0736734999054800e+01, new int[] { 0, 0, 6 });
            p.AddCoeff(-9.4257886359339900e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(1.9794156135461400e+02, new int[] { 0, 2, 2 });
            p.AddCoeff(-5.9382468406384100e+02, new int[] { 0, 2, 4 });
            p.AddCoeff(4.3547143498015000e+02, new int[] { 0, 2, 6 });
            p.AddCoeff(2.8277365907802000e+01, new int[] { 0, 4, 0 });
            p.AddCoeff(-5.9382468406384100e+02, new int[] { 0, 4, 2 });
            p.AddCoeff(1.7814740521915200e+03, new int[] { 0, 4, 4 });
            p.AddCoeff(-1.3064143049404500e+03, new int[] { 0, 4, 6 });
            p.AddCoeff(-2.0736734999054800e+01, new int[] { 0, 6, 0 });
            p.AddCoeff(4.3547143498015000e+02, new int[] { 0, 6, 2 });
            p.AddCoeff(-1.3064143049404500e+03, new int[] { 0, 6, 4 });
            p.AddCoeff(9.5803715695633100e+02, new int[] { 0, 6, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5b1cbdd5-40ef-49ea-aedd-37dc67566444"));
            OrthonormalPolynomials[371] = p;
            p.AddCoeff(-1.8627145733216900e+01, new int[] { 0, 1, 1 });
            p.AddCoeff(8.6926680088345700e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(-7.8234012079511100e+01, new int[] { 0, 1, 5 });
            p.AddCoeff(1.6764431159895200e+02, new int[] { 0, 3, 1 });
            p.AddCoeff(-7.8234012079511100e+02, new int[] { 0, 3, 3 });
            p.AddCoeff(7.0410610871560000e+02, new int[] { 0, 3, 5 });
            p.AddCoeff(-3.6881748551769500e+02, new int[] { 0, 5, 1 });
            p.AddCoeff(1.7211482657492500e+03, new int[] { 0, 5, 3 });
            p.AddCoeff(-1.5490334391743200e+03, new int[] { 0, 5, 5 });
            p.AddCoeff(2.2831558627285900e+02, new int[] { 0, 7, 1 });
            p.AddCoeff(-1.0654727359400100e+03, new int[] { 0, 7, 3 });
            p.AddCoeff(9.5892546234600800e+02, new int[] { 0, 7, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("34771612-8d02-429b-8910-782c1dbefcf5"));
            OrthonormalPolynomials[372] = p;
            p.AddCoeff(4.4842525558502700e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-4.4842525558502700e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(5.2316279818253100e+00, new int[] { 0, 0, 4 });
            p.AddCoeff(-1.6143309201061000e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(1.6143309201061000e+02, new int[] { 0, 2, 2 });
            p.AddCoeff(-1.8833860734571100e+02, new int[] { 0, 2, 4 });
            p.AddCoeff(8.8788200605835300e+01, new int[] { 0, 4, 0 });
            p.AddCoeff(-8.8788200605835300e+02, new int[] { 0, 4, 2 });
            p.AddCoeff(1.0358623404014100e+03, new int[] { 0, 4, 4 });
            p.AddCoeff(-1.5389954771678100e+02, new int[] { 0, 6, 0 });
            p.AddCoeff(1.5389954771678100e+03, new int[] { 0, 6, 2 });
            p.AddCoeff(-1.7954947233624500e+03, new int[] { 0, 6, 4 });
            p.AddCoeff(8.2446186276847100e+01, new int[] { 0, 8, 0 });
            p.AddCoeff(-8.2446186276847100e+02, new int[] { 0, 8, 2 });
            p.AddCoeff(9.6187217322988200e+02, new int[] { 0, 8, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("4b176015-7e96-41ef-bb21-1306513f96fa"));
            OrthonormalPolynomials[373] = p;
            p.AddCoeff(-1.5051253492806200e+01, new int[] { 0, 1, 1 });
            p.AddCoeff(2.5085422488010400e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(2.2075171789449200e+02, new int[] { 0, 3, 1 });
            p.AddCoeff(-3.6791952982415300e+02, new int[] { 0, 3, 3 });
            p.AddCoeff(-8.6093169978851700e+02, new int[] { 0, 5, 1 });
            p.AddCoeff(1.4348861663142000e+03, new int[] { 0, 5, 3 });
            p.AddCoeff(1.2299024282693100e+03, new int[] { 0, 7, 1 });
            p.AddCoeff(-2.0498373804488500e+03, new int[] { 0, 7, 3 });
            p.AddCoeff(-5.8078725779384100e+02, new int[] { 0, 9, 1 });
            p.AddCoeff(9.6797876298973500e+02, new int[] { 0, 9, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c0330d0e-7592-4717-9954-e7c88bbe0f7c"));
            OrthonormalPolynomials[374] = p;
            p.AddCoeff(4.4577965576656800e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-1.3373389672997100e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-2.4517881067161300e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(7.3553643201483800e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(2.1248830258206400e+02, new int[] { 0, 4, 0 });
            p.AddCoeff(-6.3746490774619300e+02, new int[] { 0, 4, 2 });
            p.AddCoeff(-6.3746490774619300e+02, new int[] { 0, 6, 0 });
            p.AddCoeff(1.9123947232385800e+03, new int[] { 0, 6, 2 });
            p.AddCoeff(7.7406453083466300e+02, new int[] { 0, 8, 0 });
            p.AddCoeff(-2.3221935925039900e+03, new int[] { 0, 8, 2 });
            p.AddCoeff(-3.2682724635241300e+02, new int[] { 0, 10, 0 });
            p.AddCoeff(9.8048173905724000e+02, new int[] { 0, 10, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f3012a3a-968e-41cd-b5ac-2c6603eafb94"));
            OrthonormalPolynomials[375] = p;
            p.AddCoeff(-7.9501042053302700e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(1.7225225778215600e+02, new int[] { 0, 3, 1 });
            p.AddCoeff(-1.0335135466929300e+03, new int[] { 0, 5, 1 });
            p.AddCoeff(2.5099614705399800e+03, new int[] { 0, 7, 1 });
            p.AddCoeff(-2.6494037744588700e+03, new int[] { 0, 9, 1 });
            p.AddCoeff(1.0115905320661200e+03, new int[] { 0, 11, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d7d2d622-cdc0-42dc-b041-e75dbac768d6"));
            OrthonormalPolynomials[376] = p;
            p.AddCoeff(3.9878336536643800e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-3.1105102498582200e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(3.8881378123227700e+02, new int[] { 0, 4, 0 });
            p.AddCoeff(-1.7626224749196600e+03, new int[] { 0, 6, 0 });
            p.AddCoeff(3.5881957525150100e+03, new int[] { 0, 8, 0 });
            p.AddCoeff(-3.3489827023473500e+03, new int[] { 0, 10, 0 });
            p.AddCoeff(1.1670697296058900e+03, new int[] { 0, 12, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("2765f4fd-47c1-41a0-9d6d-d465e980ffb3"));
            OrthonormalPolynomials[377] = p;
            p.AddCoeff(-7.9501042053302700e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(1.7225225778215600e+02, new int[] { 1, 0, 3 });
            p.AddCoeff(-1.0335135466929300e+03, new int[] { 1, 0, 5 });
            p.AddCoeff(2.5099614705399800e+03, new int[] { 1, 0, 7 });
            p.AddCoeff(-2.6494037744588700e+03, new int[] { 1, 0, 9 });
            p.AddCoeff(1.0115905320661200e+03, new int[] { 1, 0, 11 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d7605c53-0ed3-47ab-bf51-e4980e751dbd"));
            OrthonormalPolynomials[378] = p;
            p.AddCoeff(-1.1961523359366100e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(6.5788378476513400e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(-5.7016594679644900e+02, new int[] { 1, 1, 4 });
            p.AddCoeff(1.7104978403893500e+03, new int[] { 1, 1, 6 });
            p.AddCoeff(-2.0770330919013500e+03, new int[] { 1, 1, 8 });
            p.AddCoeff(8.7696952769168200e+02, new int[] { 1, 1, 10 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("625f4e78-fef2-444a-9f7a-6494b7dbabf6"));
            OrthonormalPolynomials[379] = p;
            p.AddCoeff(-7.3442596907982200e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(1.0771580879837400e+02, new int[] { 1, 0, 3 });
            p.AddCoeff(-4.2009165431365800e+02, new int[] { 1, 0, 5 });
            p.AddCoeff(6.0013093473379800e+02, new int[] { 1, 0, 7 });
            p.AddCoeff(-2.8339516362429300e+02, new int[] { 1, 0, 9 });
            p.AddCoeff(2.2032779072394700e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(-3.2314742639512200e+02, new int[] { 1, 2, 3 });
            p.AddCoeff(1.2602749629409700e+03, new int[] { 1, 2, 5 });
            p.AddCoeff(-1.8003928042013900e+03, new int[] { 1, 2, 7 });
            p.AddCoeff(8.5018549087288000e+02, new int[] { 1, 2, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("ca2edc49-5da6-4a79-972b-9528da0fc3f7"));
            OrthonormalPolynomials[380] = p;
            p.AddCoeff(-2.7399235696644100e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(9.8637248507918700e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(-5.4250486679355300e+02, new int[] { 1, 1, 4 });
            p.AddCoeff(9.4034176910882500e+02, new int[] { 1, 1, 6 });
            p.AddCoeff(-5.0375451916544200e+02, new int[] { 1, 1, 8 });
            p.AddCoeff(4.5665392827740100e+00, new int[] { 1, 3, 0 });
            p.AddCoeff(-1.6439541417986400e+02, new int[] { 1, 3, 2 });
            p.AddCoeff(9.0417477798925500e+02, new int[] { 1, 3, 4 });
            p.AddCoeff(-1.5672362818480400e+03, new int[] { 1, 3, 6 });
            p.AddCoeff(8.3959086527573700e+02, new int[] { 1, 3, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d74835d6-b5ad-4933-bdd7-90058edb3624"));
            OrthonormalPolynomials[381] = p;
            p.AddCoeff(-5.8366257594904700e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(5.2529631835414200e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(-1.1556519003791100e+02, new int[] { 1, 0, 5 });
            p.AddCoeff(7.1540355737754600e+01, new int[] { 1, 0, 7 });
            p.AddCoeff(5.8366257594904700e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(-5.2529631835414200e+02, new int[] { 1, 2, 3 });
            p.AddCoeff(1.1556519003791100e+03, new int[] { 1, 2, 5 });
            p.AddCoeff(-7.1540355737754600e+02, new int[] { 1, 2, 7 });
            p.AddCoeff(-6.8093967194055400e+01, new int[] { 1, 4, 1 });
            p.AddCoeff(6.1284570474649900e+02, new int[] { 1, 4, 3 });
            p.AddCoeff(-1.3482605504423000e+03, new int[] { 1, 4, 5 });
            p.AddCoeff(8.3463748360713700e+02, new int[] { 1, 4, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("13084201-862d-4bb3-a925-4f5a7635dd67"));
            OrthonormalPolynomials[382] = p;
            p.AddCoeff(-4.2907671433232800e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(9.0106110009788900e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(-2.7031833002936700e+02, new int[] { 1, 1, 4 });
            p.AddCoeff(1.9823344202153600e+02, new int[] { 1, 1, 6 });
            p.AddCoeff(2.0023580002175300e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(-4.2049518004568200e+02, new int[] { 1, 3, 2 });
            p.AddCoeff(1.2614855401370400e+03, new int[] { 1, 3, 4 });
            p.AddCoeff(-9.2508939610049900e+02, new int[] { 1, 3, 6 });
            p.AddCoeff(-1.8021222001957800e+01, new int[] { 1, 5, 0 });
            p.AddCoeff(3.7844566204111300e+02, new int[] { 1, 5, 2 });
            p.AddCoeff(-1.1353369861233400e+03, new int[] { 1, 5, 4 });
            p.AddCoeff(8.3258045649044900e+02, new int[] { 1, 5, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("84f71ff6-159c-48d9-85fe-400a3238c1f7"));
            OrthonormalPolynomials[383] = p;
            p.AddCoeff(-4.2907671433232800e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(2.0023580002175300e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(-1.8021222001957800e+01, new int[] { 1, 0, 5 });
            p.AddCoeff(9.0106110009788900e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(-4.2049518004568200e+02, new int[] { 1, 2, 3 });
            p.AddCoeff(3.7844566204111300e+02, new int[] { 1, 2, 5 });
            p.AddCoeff(-2.7031833002936700e+02, new int[] { 1, 4, 1 });
            p.AddCoeff(1.2614855401370400e+03, new int[] { 1, 4, 3 });
            p.AddCoeff(-1.1353369861233400e+03, new int[] { 1, 4, 5 });
            p.AddCoeff(1.9823344202153600e+02, new int[] { 1, 6, 1 });
            p.AddCoeff(-9.2508939610049900e+02, new int[] { 1, 6, 3 });
            p.AddCoeff(8.3258045649044900e+02, new int[] { 1, 6, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c09cc2ca-e170-4223-b811-f5ad1ce13125"));
            OrthonormalPolynomials[384] = p;
            p.AddCoeff(-5.8366257594904700e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(5.8366257594904700e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(-6.8093967194055400e+01, new int[] { 1, 1, 4 });
            p.AddCoeff(5.2529631835414200e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(-5.2529631835414200e+02, new int[] { 1, 3, 2 });
            p.AddCoeff(6.1284570474649900e+02, new int[] { 1, 3, 4 });
            p.AddCoeff(-1.1556519003791100e+02, new int[] { 1, 5, 0 });
            p.AddCoeff(1.1556519003791100e+03, new int[] { 1, 5, 2 });
            p.AddCoeff(-1.3482605504423000e+03, new int[] { 1, 5, 4 });
            p.AddCoeff(7.1540355737754600e+01, new int[] { 1, 7, 0 });
            p.AddCoeff(-7.1540355737754600e+02, new int[] { 1, 7, 2 });
            p.AddCoeff(8.3463748360713700e+02, new int[] { 1, 7, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a950cf33-1b52-4af7-bedb-f454988b7e7e"));
            OrthonormalPolynomials[385] = p;
            p.AddCoeff(-2.7399235696644100e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(4.5665392827740100e+00, new int[] { 1, 0, 3 });
            p.AddCoeff(9.8637248507918700e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(-1.6439541417986400e+02, new int[] { 1, 2, 3 });
            p.AddCoeff(-5.4250486679355300e+02, new int[] { 1, 4, 1 });
            p.AddCoeff(9.0417477798925500e+02, new int[] { 1, 4, 3 });
            p.AddCoeff(9.4034176910882500e+02, new int[] { 1, 6, 1 });
            p.AddCoeff(-1.5672362818480400e+03, new int[] { 1, 6, 3 });
            p.AddCoeff(-5.0375451916544200e+02, new int[] { 1, 8, 1 });
            p.AddCoeff(8.3959086527573700e+02, new int[] { 1, 8, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b0ea8003-0a1e-4d64-a8c0-fde90440d016"));
            OrthonormalPolynomials[386] = p;
            p.AddCoeff(-7.3442596907982200e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(2.2032779072394700e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(1.0771580879837400e+02, new int[] { 1, 3, 0 });
            p.AddCoeff(-3.2314742639512200e+02, new int[] { 1, 3, 2 });
            p.AddCoeff(-4.2009165431365800e+02, new int[] { 1, 5, 0 });
            p.AddCoeff(1.2602749629409700e+03, new int[] { 1, 5, 2 });
            p.AddCoeff(6.0013093473379800e+02, new int[] { 1, 7, 0 });
            p.AddCoeff(-1.8003928042013900e+03, new int[] { 1, 7, 2 });
            p.AddCoeff(-2.8339516362429300e+02, new int[] { 1, 9, 0 });
            p.AddCoeff(8.5018549087288000e+02, new int[] { 1, 9, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("9e9deb13-51cc-4d63-9528-541cd6f2b890"));
            OrthonormalPolynomials[387] = p;
            p.AddCoeff(-1.1961523359366100e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(6.5788378476513400e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(-5.7016594679644900e+02, new int[] { 1, 4, 1 });
            p.AddCoeff(1.7104978403893500e+03, new int[] { 1, 6, 1 });
            p.AddCoeff(-2.0770330919013500e+03, new int[] { 1, 8, 1 });
            p.AddCoeff(8.7696952769168200e+02, new int[] { 1, 10, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("76afbed2-7537-4244-959b-9d4b5999edec"));
            OrthonormalPolynomials[388] = p;
            p.AddCoeff(-7.9501042053302700e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(1.7225225778215600e+02, new int[] { 1, 3, 0 });
            p.AddCoeff(-1.0335135466929300e+03, new int[] { 1, 5, 0 });
            p.AddCoeff(2.5099614705399800e+03, new int[] { 1, 7, 0 });
            p.AddCoeff(-2.6494037744588700e+03, new int[] { 1, 9, 0 });
            p.AddCoeff(1.0115905320661200e+03, new int[] { 1, 11, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("014a26dd-2129-4270-9e78-0f01c69c4684"));
            OrthonormalPolynomials[389] = p;
            p.AddCoeff(4.4577965576656800e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-2.4517881067161300e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(2.1248830258206400e+02, new int[] { 0, 0, 4 });
            p.AddCoeff(-6.3746490774619300e+02, new int[] { 0, 0, 6 });
            p.AddCoeff(7.7406453083466300e+02, new int[] { 0, 0, 8 });
            p.AddCoeff(-3.2682724635241300e+02, new int[] { 0, 0, 10 });
            p.AddCoeff(-1.3373389672997100e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(7.3553643201483800e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(-6.3746490774619300e+02, new int[] { 2, 0, 4 });
            p.AddCoeff(1.9123947232385800e+03, new int[] { 2, 0, 6 });
            p.AddCoeff(-2.3221935925039900e+03, new int[] { 2, 0, 8 });
            p.AddCoeff(9.8048173905724000e+02, new int[] { 2, 0, 10 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7ed95b50-911a-4f23-8a12-755cc8f4301c"));
            OrthonormalPolynomials[390] = p;
            p.AddCoeff(-7.3442596907982200e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(1.0771580879837400e+02, new int[] { 0, 1, 3 });
            p.AddCoeff(-4.2009165431365800e+02, new int[] { 0, 1, 5 });
            p.AddCoeff(6.0013093473379800e+02, new int[] { 0, 1, 7 });
            p.AddCoeff(-2.8339516362429300e+02, new int[] { 0, 1, 9 });
            p.AddCoeff(2.2032779072394700e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(-3.2314742639512200e+02, new int[] { 2, 1, 3 });
            p.AddCoeff(1.2602749629409700e+03, new int[] { 2, 1, 5 });
            p.AddCoeff(-1.8003928042013900e+03, new int[] { 2, 1, 7 });
            p.AddCoeff(8.5018549087288000e+02, new int[] { 2, 1, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("804a04e9-e940-4ce7-8d8e-bc0d8136500d"));
            OrthonormalPolynomials[391] = p;
            p.AddCoeff(4.9825028398336300e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-1.7937010223401100e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(9.8653556228705900e+01, new int[] { 0, 0, 4 });
            p.AddCoeff(-1.7099949746309000e+02, new int[] { 0, 0, 6 });
            p.AddCoeff(9.1606873640941200e+01, new int[] { 0, 0, 8 });
            p.AddCoeff(-1.4947508519500900e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(5.3811030670203200e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(-2.9596066868611800e+02, new int[] { 0, 2, 4 });
            p.AddCoeff(5.1299849238927100e+02, new int[] { 0, 2, 6 });
            p.AddCoeff(-2.7482062092282400e+02, new int[] { 0, 2, 8 });
            p.AddCoeff(-1.4947508519500900e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(5.3811030670203200e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(-2.9596066868611800e+02, new int[] { 2, 0, 4 });
            p.AddCoeff(5.1299849238927100e+02, new int[] { 2, 0, 6 });
            p.AddCoeff(-2.7482062092282400e+02, new int[] { 2, 0, 8 });
            p.AddCoeff(4.4842525558502700e+00, new int[] { 2, 2, 0 });
            p.AddCoeff(-1.6143309201061000e+02, new int[] { 2, 2, 2 });
            p.AddCoeff(8.8788200605835300e+02, new int[] { 2, 2, 4 });
            p.AddCoeff(-1.5389954771678100e+03, new int[] { 2, 2, 6 });
            p.AddCoeff(8.2446186276847100e+02, new int[] { 2, 2, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("9480dee4-44f4-41f7-a157-dec5c276f33b"));
            OrthonormalPolynomials[392] = p;
            p.AddCoeff(-1.3290581510406700e+01, new int[] { 0, 1, 1 });
            p.AddCoeff(1.1961523359366100e+02, new int[] { 0, 1, 3 });
            p.AddCoeff(-2.6315351390605400e+02, new int[] { 0, 1, 5 });
            p.AddCoeff(1.6290455622755700e+02, new int[] { 0, 1, 7 });
            p.AddCoeff(2.2150969184011200e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(-1.9935872265610100e+02, new int[] { 0, 3, 3 });
            p.AddCoeff(4.3858918984342300e+02, new int[] { 0, 3, 5 });
            p.AddCoeff(-2.7150759371259500e+02, new int[] { 0, 3, 7 });
            p.AddCoeff(3.9871744531220200e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(-3.5884570078098200e+02, new int[] { 2, 1, 3 });
            p.AddCoeff(7.8946054171816100e+02, new int[] { 2, 1, 5 });
            p.AddCoeff(-4.8871366868267100e+02, new int[] { 2, 1, 7 });
            p.AddCoeff(-6.6452907552033700e+01, new int[] { 2, 3, 1 });
            p.AddCoeff(5.9807616796830300e+02, new int[] { 2, 3, 3 });
            p.AddCoeff(-1.3157675695302700e+03, new int[] { 2, 3, 5 });
            p.AddCoeff(8.1452278113778500e+02, new int[] { 2, 3, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b6fb43b7-ba68-4693-b5e9-143871423385"));
            OrthonormalPolynomials[393] = p;
            p.AddCoeff(5.0105365360802100e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-1.0522126725768400e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(3.1566380177305300e+01, new int[] { 0, 0, 4 });
            p.AddCoeff(-2.3148678796690500e+01, new int[] { 0, 0, 6 });
            p.AddCoeff(-5.0105365360802100e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(1.0522126725768400e+02, new int[] { 0, 2, 2 });
            p.AddCoeff(-3.1566380177305300e+02, new int[] { 0, 2, 4 });
            p.AddCoeff(2.3148678796690500e+02, new int[] { 0, 2, 6 });
            p.AddCoeff(5.8456259587602400e+00, new int[] { 0, 4, 0 });
            p.AddCoeff(-1.2275814513396500e+02, new int[] { 0, 4, 2 });
            p.AddCoeff(3.6827443540189500e+02, new int[] { 0, 4, 4 });
            p.AddCoeff(-2.7006791929472300e+02, new int[] { 0, 4, 6 });
            p.AddCoeff(-1.5031609608240600e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(3.1566380177305300e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(-9.4699140531915900e+01, new int[] { 2, 0, 4 });
            p.AddCoeff(6.9446036390071600e+01, new int[] { 2, 0, 6 });
            p.AddCoeff(1.5031609608240600e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(-3.1566380177305300e+02, new int[] { 2, 2, 2 });
            p.AddCoeff(9.4699140531915900e+02, new int[] { 2, 2, 4 });
            p.AddCoeff(-6.9446036390071600e+02, new int[] { 2, 2, 6 });
            p.AddCoeff(-1.7536877876280700e+01, new int[] { 2, 4, 0 });
            p.AddCoeff(3.6827443540189500e+02, new int[] { 2, 4, 2 });
            p.AddCoeff(-1.1048233062056900e+03, new int[] { 2, 4, 4 });
            p.AddCoeff(8.1020375788416900e+02, new int[] { 2, 4, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("13e9f1ea-8ee6-4e40-a54f-a38ff45d37a3"));
            OrthonormalPolynomials[394] = p;
            p.AddCoeff(-1.5286400798665500e+01, new int[] { 0, 1, 1 });
            p.AddCoeff(7.1336537060439000e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(-6.4202883354395100e+01, new int[] { 0, 1, 5 });
            p.AddCoeff(7.1336537060439000e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(-3.3290383961538200e+02, new int[] { 0, 3, 3 });
            p.AddCoeff(2.9961345565384400e+02, new int[] { 0, 3, 5 });
            p.AddCoeff(-6.4202883354395100e+01, new int[] { 0, 5, 1 });
            p.AddCoeff(2.9961345565384400e+02, new int[] { 0, 5, 3 });
            p.AddCoeff(-2.6965211008846000e+02, new int[] { 0, 5, 5 });
            p.AddCoeff(4.5859202395996500e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(-2.1400961118131700e+02, new int[] { 2, 1, 3 });
            p.AddCoeff(1.9260865006318500e+02, new int[] { 2, 1, 5 });
            p.AddCoeff(-2.1400961118131700e+02, new int[] { 2, 3, 1 });
            p.AddCoeff(9.9871151884614700e+02, new int[] { 2, 3, 3 });
            p.AddCoeff(-8.9884036696153200e+02, new int[] { 2, 3, 5 });
            p.AddCoeff(1.9260865006318500e+02, new int[] { 2, 5, 1 });
            p.AddCoeff(-8.9884036696153200e+02, new int[] { 2, 5, 3 });
            p.AddCoeff(8.0895633026537900e+02, new int[] { 2, 5, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("014e860a-6d1c-4908-a663-005ba7e15426"));
            OrthonormalPolynomials[395] = p;
            p.AddCoeff(5.0105365360802100e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-5.0105365360802100e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(5.8456259587602400e+00, new int[] { 0, 0, 4 });
            p.AddCoeff(-1.0522126725768400e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(1.0522126725768400e+02, new int[] { 0, 2, 2 });
            p.AddCoeff(-1.2275814513396500e+02, new int[] { 0, 2, 4 });
            p.AddCoeff(3.1566380177305300e+01, new int[] { 0, 4, 0 });
            p.AddCoeff(-3.1566380177305300e+02, new int[] { 0, 4, 2 });
            p.AddCoeff(3.6827443540189500e+02, new int[] { 0, 4, 4 });
            p.AddCoeff(-2.3148678796690500e+01, new int[] { 0, 6, 0 });
            p.AddCoeff(2.3148678796690500e+02, new int[] { 0, 6, 2 });
            p.AddCoeff(-2.7006791929472300e+02, new int[] { 0, 6, 4 });
            p.AddCoeff(-1.5031609608240600e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(1.5031609608240600e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(-1.7536877876280700e+01, new int[] { 2, 0, 4 });
            p.AddCoeff(3.1566380177305300e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(-3.1566380177305300e+02, new int[] { 2, 2, 2 });
            p.AddCoeff(3.6827443540189500e+02, new int[] { 2, 2, 4 });
            p.AddCoeff(-9.4699140531915900e+01, new int[] { 2, 4, 0 });
            p.AddCoeff(9.4699140531915900e+02, new int[] { 2, 4, 2 });
            p.AddCoeff(-1.1048233062056900e+03, new int[] { 2, 4, 4 });
            p.AddCoeff(6.9446036390071600e+01, new int[] { 2, 6, 0 });
            p.AddCoeff(-6.9446036390071600e+02, new int[] { 2, 6, 2 });
            p.AddCoeff(8.1020375788416900e+02, new int[] { 2, 6, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("20a94497-988e-497b-b196-9e2a7d947982"));
            OrthonormalPolynomials[396] = p;
            p.AddCoeff(-1.3290581510406700e+01, new int[] { 0, 1, 1 });
            p.AddCoeff(2.2150969184011200e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(1.1961523359366100e+02, new int[] { 0, 3, 1 });
            p.AddCoeff(-1.9935872265610100e+02, new int[] { 0, 3, 3 });
            p.AddCoeff(-2.6315351390605400e+02, new int[] { 0, 5, 1 });
            p.AddCoeff(4.3858918984342300e+02, new int[] { 0, 5, 3 });
            p.AddCoeff(1.6290455622755700e+02, new int[] { 0, 7, 1 });
            p.AddCoeff(-2.7150759371259500e+02, new int[] { 0, 7, 3 });
            p.AddCoeff(3.9871744531220200e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(-6.6452907552033700e+01, new int[] { 2, 1, 3 });
            p.AddCoeff(-3.5884570078098200e+02, new int[] { 2, 3, 1 });
            p.AddCoeff(5.9807616796830300e+02, new int[] { 2, 3, 3 });
            p.AddCoeff(7.8946054171816100e+02, new int[] { 2, 5, 1 });
            p.AddCoeff(-1.3157675695302700e+03, new int[] { 2, 5, 3 });
            p.AddCoeff(-4.8871366868267100e+02, new int[] { 2, 7, 1 });
            p.AddCoeff(8.1452278113778500e+02, new int[] { 2, 7, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("885ec7c4-9d6d-4d3d-88fd-eed92ea461a2"));
            OrthonormalPolynomials[397] = p;
            p.AddCoeff(4.9825028398336300e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-1.4947508519500900e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-1.7937010223401100e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(5.3811030670203200e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(9.8653556228705900e+01, new int[] { 0, 4, 0 });
            p.AddCoeff(-2.9596066868611800e+02, new int[] { 0, 4, 2 });
            p.AddCoeff(-1.7099949746309000e+02, new int[] { 0, 6, 0 });
            p.AddCoeff(5.1299849238927100e+02, new int[] { 0, 6, 2 });
            p.AddCoeff(9.1606873640941200e+01, new int[] { 0, 8, 0 });
            p.AddCoeff(-2.7482062092282400e+02, new int[] { 0, 8, 2 });
            p.AddCoeff(-1.4947508519500900e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(4.4842525558502700e+00, new int[] { 2, 0, 2 });
            p.AddCoeff(5.3811030670203200e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(-1.6143309201061000e+02, new int[] { 2, 2, 2 });
            p.AddCoeff(-2.9596066868611800e+02, new int[] { 2, 4, 0 });
            p.AddCoeff(8.8788200605835300e+02, new int[] { 2, 4, 2 });
            p.AddCoeff(5.1299849238927100e+02, new int[] { 2, 6, 0 });
            p.AddCoeff(-1.5389954771678100e+03, new int[] { 2, 6, 2 });
            p.AddCoeff(-2.7482062092282400e+02, new int[] { 2, 8, 0 });
            p.AddCoeff(8.2446186276847100e+02, new int[] { 2, 8, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3292e2aa-9f19-43dd-b207-9814497eb09f"));
            OrthonormalPolynomials[398] = p;
            p.AddCoeff(-7.3442596907982200e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(1.0771580879837400e+02, new int[] { 0, 3, 1 });
            p.AddCoeff(-4.2009165431365800e+02, new int[] { 0, 5, 1 });
            p.AddCoeff(6.0013093473379800e+02, new int[] { 0, 7, 1 });
            p.AddCoeff(-2.8339516362429300e+02, new int[] { 0, 9, 1 });
            p.AddCoeff(2.2032779072394700e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(-3.2314742639512200e+02, new int[] { 2, 3, 1 });
            p.AddCoeff(1.2602749629409700e+03, new int[] { 2, 5, 1 });
            p.AddCoeff(-1.8003928042013900e+03, new int[] { 2, 7, 1 });
            p.AddCoeff(8.5018549087288000e+02, new int[] { 2, 9, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("263b7020-ed54-4859-b150-73e007b36c6d"));
            OrthonormalPolynomials[399] = p;
            p.AddCoeff(4.4577965576656800e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-2.4517881067161300e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(2.1248830258206400e+02, new int[] { 0, 4, 0 });
            p.AddCoeff(-6.3746490774619300e+02, new int[] { 0, 6, 0 });
            p.AddCoeff(7.7406453083466300e+02, new int[] { 0, 8, 0 });
            p.AddCoeff(-3.2682724635241300e+02, new int[] { 0, 10, 0 });
            p.AddCoeff(-1.3373389672997100e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(7.3553643201483800e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(-6.3746490774619300e+02, new int[] { 2, 4, 0 });
            p.AddCoeff(1.9123947232385800e+03, new int[] { 2, 6, 0 });
            p.AddCoeff(-2.3221935925039900e+03, new int[] { 2, 8, 0 });
            p.AddCoeff(9.8048173905724000e+02, new int[] { 2, 10, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("6f68c2fa-a0ef-4a45-9c2c-fe40b4278011"));
            OrthonormalPolynomials[400] = p;
            p.AddCoeff(-1.5051253492806200e+01, new int[] { 1, 0, 1 });
            p.AddCoeff(2.2075171789449200e+02, new int[] { 1, 0, 3 });
            p.AddCoeff(-8.6093169978851700e+02, new int[] { 1, 0, 5 });
            p.AddCoeff(1.2299024282693100e+03, new int[] { 1, 0, 7 });
            p.AddCoeff(-5.8078725779384100e+02, new int[] { 1, 0, 9 });
            p.AddCoeff(2.5085422488010400e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(-3.6791952982415300e+02, new int[] { 3, 0, 3 });
            p.AddCoeff(1.4348861663142000e+03, new int[] { 3, 0, 5 });
            p.AddCoeff(-2.0498373804488500e+03, new int[] { 3, 0, 7 });
            p.AddCoeff(9.6797876298973500e+02, new int[] { 3, 0, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e4e4d3f3-276d-4795-9fb4-2b729cb62889"));
            OrthonormalPolynomials[401] = p;
            p.AddCoeff(-2.7399235696644100e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(9.8637248507918700e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(-5.4250486679355300e+02, new int[] { 1, 1, 4 });
            p.AddCoeff(9.4034176910882500e+02, new int[] { 1, 1, 6 });
            p.AddCoeff(-5.0375451916544200e+02, new int[] { 1, 1, 8 });
            p.AddCoeff(4.5665392827740100e+00, new int[] { 3, 1, 0 });
            p.AddCoeff(-1.6439541417986400e+02, new int[] { 3, 1, 2 });
            p.AddCoeff(9.0417477798925500e+02, new int[] { 3, 1, 4 });
            p.AddCoeff(-1.5672362818480400e+03, new int[] { 3, 1, 6 });
            p.AddCoeff(8.3959086527573700e+02, new int[] { 3, 1, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("8656bd05-1f07-486e-a181-937a23630fcc"));
            OrthonormalPolynomials[402] = p;
            p.AddCoeff(-1.3290581510406700e+01, new int[] { 1, 0, 1 });
            p.AddCoeff(1.1961523359366100e+02, new int[] { 1, 0, 3 });
            p.AddCoeff(-2.6315351390605400e+02, new int[] { 1, 0, 5 });
            p.AddCoeff(1.6290455622755700e+02, new int[] { 1, 0, 7 });
            p.AddCoeff(3.9871744531220200e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(-3.5884570078098200e+02, new int[] { 1, 2, 3 });
            p.AddCoeff(7.8946054171816100e+02, new int[] { 1, 2, 5 });
            p.AddCoeff(-4.8871366868267100e+02, new int[] { 1, 2, 7 });
            p.AddCoeff(2.2150969184011200e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(-1.9935872265610100e+02, new int[] { 3, 0, 3 });
            p.AddCoeff(4.3858918984342300e+02, new int[] { 3, 0, 5 });
            p.AddCoeff(-2.7150759371259500e+02, new int[] { 3, 0, 7 });
            p.AddCoeff(-6.6452907552033700e+01, new int[] { 3, 2, 1 });
            p.AddCoeff(5.9807616796830300e+02, new int[] { 3, 2, 3 });
            p.AddCoeff(-1.3157675695302700e+03, new int[] { 3, 2, 5 });
            p.AddCoeff(8.1452278113778500e+02, new int[] { 3, 2, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("37921e8d-b6ca-4fe4-8dff-3d1a1443c325"));
            OrthonormalPolynomials[403] = p;
            p.AddCoeff(-6.2741841671161200e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(1.3175786750943800e+02, new int[] { 1, 1, 2 });
            p.AddCoeff(-3.9527360252831500e+02, new int[] { 1, 1, 4 });
            p.AddCoeff(2.8986730852076500e+02, new int[] { 1, 1, 6 });
            p.AddCoeff(1.0456973611860200e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(-2.1959644584906400e+02, new int[] { 1, 3, 2 });
            p.AddCoeff(6.5878933754719200e+02, new int[] { 1, 3, 4 });
            p.AddCoeff(-4.8311218086794100e+02, new int[] { 1, 3, 6 });
            p.AddCoeff(1.0456973611860200e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(-2.1959644584906400e+02, new int[] { 3, 1, 2 });
            p.AddCoeff(6.5878933754719200e+02, new int[] { 3, 1, 4 });
            p.AddCoeff(-4.8311218086794100e+02, new int[] { 3, 1, 6 });
            p.AddCoeff(-1.7428289353100300e+01, new int[] { 3, 3, 0 });
            p.AddCoeff(3.6599407641510700e+02, new int[] { 3, 3, 2 });
            p.AddCoeff(-1.0979822292453200e+03, new int[] { 3, 3, 4 });
            p.AddCoeff(8.0518696811323500e+02, new int[] { 3, 3, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f395b5cb-2be9-4168-8bf0-7d5656f29a51"));
            OrthonormalPolynomials[404] = p;
            p.AddCoeff(-9.8162457551295200e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(4.5809146857271100e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(-4.1228232171544000e+01, new int[] { 1, 0, 5 });
            p.AddCoeff(9.8162457551295200e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(-4.5809146857271100e+02, new int[] { 1, 2, 3 });
            p.AddCoeff(4.1228232171544000e+02, new int[] { 1, 2, 5 });
            p.AddCoeff(-1.1452286714317800e+02, new int[] { 1, 4, 1 });
            p.AddCoeff(5.3444004666816300e+02, new int[] { 1, 4, 3 });
            p.AddCoeff(-4.8099604200134700e+02, new int[] { 1, 4, 5 });
            p.AddCoeff(1.6360409591882500e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(-7.6348578095451800e+01, new int[] { 3, 0, 3 });
            p.AddCoeff(6.8713720285906600e+01, new int[] { 3, 0, 5 });
            p.AddCoeff(-1.6360409591882500e+02, new int[] { 3, 2, 1 });
            p.AddCoeff(7.6348578095451800e+02, new int[] { 3, 2, 3 });
            p.AddCoeff(-6.8713720285906600e+02, new int[] { 3, 2, 5 });
            p.AddCoeff(1.9087144523863000e+02, new int[] { 3, 4, 1 });
            p.AddCoeff(-8.9073341111360500e+02, new int[] { 3, 4, 3 });
            p.AddCoeff(8.0166007000224400e+02, new int[] { 3, 4, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("2ce921e8-8d07-4821-ad3a-ea2dc77c0a56"));
            OrthonormalPolynomials[405] = p;
            p.AddCoeff(-9.8162457551295200e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(9.8162457551295200e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(-1.1452286714317800e+02, new int[] { 1, 1, 4 });
            p.AddCoeff(4.5809146857271100e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(-4.5809146857271100e+02, new int[] { 1, 3, 2 });
            p.AddCoeff(5.3444004666816300e+02, new int[] { 1, 3, 4 });
            p.AddCoeff(-4.1228232171544000e+01, new int[] { 1, 5, 0 });
            p.AddCoeff(4.1228232171544000e+02, new int[] { 1, 5, 2 });
            p.AddCoeff(-4.8099604200134700e+02, new int[] { 1, 5, 4 });
            p.AddCoeff(1.6360409591882500e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(-1.6360409591882500e+02, new int[] { 3, 1, 2 });
            p.AddCoeff(1.9087144523863000e+02, new int[] { 3, 1, 4 });
            p.AddCoeff(-7.6348578095451800e+01, new int[] { 3, 3, 0 });
            p.AddCoeff(7.6348578095451800e+02, new int[] { 3, 3, 2 });
            p.AddCoeff(-8.9073341111360500e+02, new int[] { 3, 3, 4 });
            p.AddCoeff(6.8713720285906600e+01, new int[] { 3, 5, 0 });
            p.AddCoeff(-6.8713720285906600e+02, new int[] { 3, 5, 2 });
            p.AddCoeff(8.0166007000224400e+02, new int[] { 3, 5, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("42def108-37e8-422d-aa4a-978a4b857b83"));
            OrthonormalPolynomials[406] = p;
            p.AddCoeff(-6.2741841671161200e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(1.0456973611860200e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(1.3175786750943800e+02, new int[] { 1, 2, 1 });
            p.AddCoeff(-2.1959644584906400e+02, new int[] { 1, 2, 3 });
            p.AddCoeff(-3.9527360252831500e+02, new int[] { 1, 4, 1 });
            p.AddCoeff(6.5878933754719200e+02, new int[] { 1, 4, 3 });
            p.AddCoeff(2.8986730852076500e+02, new int[] { 1, 6, 1 });
            p.AddCoeff(-4.8311218086794100e+02, new int[] { 1, 6, 3 });
            p.AddCoeff(1.0456973611860200e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(-1.7428289353100300e+01, new int[] { 3, 0, 3 });
            p.AddCoeff(-2.1959644584906400e+02, new int[] { 3, 2, 1 });
            p.AddCoeff(3.6599407641510700e+02, new int[] { 3, 2, 3 });
            p.AddCoeff(6.5878933754719200e+02, new int[] { 3, 4, 1 });
            p.AddCoeff(-1.0979822292453200e+03, new int[] { 3, 4, 3 });
            p.AddCoeff(-4.8311218086794100e+02, new int[] { 3, 6, 1 });
            p.AddCoeff(8.0518696811323500e+02, new int[] { 3, 6, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e4a36c22-0620-4f70-9550-c49ee07a8fbe"));
            OrthonormalPolynomials[407] = p;
            p.AddCoeff(-1.3290581510406700e+01, new int[] { 1, 1, 0 });
            p.AddCoeff(3.9871744531220200e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(1.1961523359366100e+02, new int[] { 1, 3, 0 });
            p.AddCoeff(-3.5884570078098200e+02, new int[] { 1, 3, 2 });
            p.AddCoeff(-2.6315351390605400e+02, new int[] { 1, 5, 0 });
            p.AddCoeff(7.8946054171816100e+02, new int[] { 1, 5, 2 });
            p.AddCoeff(1.6290455622755700e+02, new int[] { 1, 7, 0 });
            p.AddCoeff(-4.8871366868267100e+02, new int[] { 1, 7, 2 });
            p.AddCoeff(2.2150969184011200e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(-6.6452907552033700e+01, new int[] { 3, 1, 2 });
            p.AddCoeff(-1.9935872265610100e+02, new int[] { 3, 3, 0 });
            p.AddCoeff(5.9807616796830300e+02, new int[] { 3, 3, 2 });
            p.AddCoeff(4.3858918984342300e+02, new int[] { 3, 5, 0 });
            p.AddCoeff(-1.3157675695302700e+03, new int[] { 3, 5, 2 });
            p.AddCoeff(-2.7150759371259500e+02, new int[] { 3, 7, 0 });
            p.AddCoeff(8.1452278113778500e+02, new int[] { 3, 7, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("87d27a2c-8a2c-47a8-be7b-2d76cac31e53"));
            OrthonormalPolynomials[408] = p;
            p.AddCoeff(-2.7399235696644100e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(9.8637248507918700e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(-5.4250486679355300e+02, new int[] { 1, 4, 1 });
            p.AddCoeff(9.4034176910882500e+02, new int[] { 1, 6, 1 });
            p.AddCoeff(-5.0375451916544200e+02, new int[] { 1, 8, 1 });
            p.AddCoeff(4.5665392827740100e+00, new int[] { 3, 0, 1 });
            p.AddCoeff(-1.6439541417986400e+02, new int[] { 3, 2, 1 });
            p.AddCoeff(9.0417477798925500e+02, new int[] { 3, 4, 1 });
            p.AddCoeff(-1.5672362818480400e+03, new int[] { 3, 6, 1 });
            p.AddCoeff(8.3959086527573700e+02, new int[] { 3, 8, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7f4ac321-3e1a-43f9-a285-e38d3e8f7c18"));
            OrthonormalPolynomials[409] = p;
            p.AddCoeff(-1.5051253492806200e+01, new int[] { 1, 1, 0 });
            p.AddCoeff(2.2075171789449200e+02, new int[] { 1, 3, 0 });
            p.AddCoeff(-8.6093169978851700e+02, new int[] { 1, 5, 0 });
            p.AddCoeff(1.2299024282693100e+03, new int[] { 1, 7, 0 });
            p.AddCoeff(-5.8078725779384100e+02, new int[] { 1, 9, 0 });
            p.AddCoeff(2.5085422488010400e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(-3.6791952982415300e+02, new int[] { 3, 3, 0 });
            p.AddCoeff(1.4348861663142000e+03, new int[] { 3, 5, 0 });
            p.AddCoeff(-2.0498373804488500e+03, new int[] { 3, 7, 0 });
            p.AddCoeff(9.6797876298973500e+02, new int[] { 3, 9, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("4384de1c-b627-4945-982b-803bc5e05795"));
            OrthonormalPolynomials[410] = p;
            p.AddCoeff(4.4842525558502700e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-1.6143309201061000e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(8.8788200605835300e+01, new int[] { 0, 0, 4 });
            p.AddCoeff(-1.5389954771678100e+02, new int[] { 0, 0, 6 });
            p.AddCoeff(8.2446186276847100e+01, new int[] { 0, 0, 8 });
            p.AddCoeff(-4.4842525558502700e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(1.6143309201061000e+02, new int[] { 2, 0, 2 });
            p.AddCoeff(-8.8788200605835300e+02, new int[] { 2, 0, 4 });
            p.AddCoeff(1.5389954771678100e+03, new int[] { 2, 0, 6 });
            p.AddCoeff(-8.2446186276847100e+02, new int[] { 2, 0, 8 });
            p.AddCoeff(5.2316279818253100e+00, new int[] { 4, 0, 0 });
            p.AddCoeff(-1.8833860734571100e+02, new int[] { 4, 0, 2 });
            p.AddCoeff(1.0358623404014100e+03, new int[] { 4, 0, 4 });
            p.AddCoeff(-1.7954947233624500e+03, new int[] { 4, 0, 6 });
            p.AddCoeff(9.6187217322988200e+02, new int[] { 4, 0, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("2925d18a-2bde-4e52-b6cb-5733727ad135"));
            OrthonormalPolynomials[411] = p;
            p.AddCoeff(-5.8366257594904700e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(5.2529631835414200e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(-1.1556519003791100e+02, new int[] { 0, 1, 5 });
            p.AddCoeff(7.1540355737754600e+01, new int[] { 0, 1, 7 });
            p.AddCoeff(5.8366257594904700e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(-5.2529631835414200e+02, new int[] { 2, 1, 3 });
            p.AddCoeff(1.1556519003791100e+03, new int[] { 2, 1, 5 });
            p.AddCoeff(-7.1540355737754600e+02, new int[] { 2, 1, 7 });
            p.AddCoeff(-6.8093967194055400e+01, new int[] { 4, 1, 1 });
            p.AddCoeff(6.1284570474649900e+02, new int[] { 4, 1, 3 });
            p.AddCoeff(-1.3482605504423000e+03, new int[] { 4, 1, 5 });
            p.AddCoeff(8.3463748360713700e+02, new int[] { 4, 1, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3bdfdac3-7040-4965-84f5-372e91c232e4"));
            OrthonormalPolynomials[412] = p;
            p.AddCoeff(5.0105365360802100e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-1.0522126725768400e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(3.1566380177305300e+01, new int[] { 0, 0, 4 });
            p.AddCoeff(-2.3148678796690500e+01, new int[] { 0, 0, 6 });
            p.AddCoeff(-1.5031609608240600e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(3.1566380177305300e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(-9.4699140531915900e+01, new int[] { 0, 2, 4 });
            p.AddCoeff(6.9446036390071600e+01, new int[] { 0, 2, 6 });
            p.AddCoeff(-5.0105365360802100e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(1.0522126725768400e+02, new int[] { 2, 0, 2 });
            p.AddCoeff(-3.1566380177305300e+02, new int[] { 2, 0, 4 });
            p.AddCoeff(2.3148678796690500e+02, new int[] { 2, 0, 6 });
            p.AddCoeff(1.5031609608240600e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(-3.1566380177305300e+02, new int[] { 2, 2, 2 });
            p.AddCoeff(9.4699140531915900e+02, new int[] { 2, 2, 4 });
            p.AddCoeff(-6.9446036390071600e+02, new int[] { 2, 2, 6 });
            p.AddCoeff(5.8456259587602400e+00, new int[] { 4, 0, 0 });
            p.AddCoeff(-1.2275814513396500e+02, new int[] { 4, 0, 2 });
            p.AddCoeff(3.6827443540189500e+02, new int[] { 4, 0, 4 });
            p.AddCoeff(-2.7006791929472300e+02, new int[] { 4, 0, 6 });
            p.AddCoeff(-1.7536877876280700e+01, new int[] { 4, 2, 0 });
            p.AddCoeff(3.6827443540189500e+02, new int[] { 4, 2, 2 });
            p.AddCoeff(-1.1048233062056900e+03, new int[] { 4, 2, 4 });
            p.AddCoeff(8.1020375788416900e+02, new int[] { 4, 2, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f99ca43b-951f-4515-8afc-fc8e4becad0e"));
            OrthonormalPolynomials[413] = p;
            p.AddCoeff(-9.8162457551295200e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(4.5809146857271100e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(-4.1228232171544000e+01, new int[] { 0, 1, 5 });
            p.AddCoeff(1.6360409591882500e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(-7.6348578095451800e+01, new int[] { 0, 3, 3 });
            p.AddCoeff(6.8713720285906600e+01, new int[] { 0, 3, 5 });
            p.AddCoeff(9.8162457551295200e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(-4.5809146857271100e+02, new int[] { 2, 1, 3 });
            p.AddCoeff(4.1228232171544000e+02, new int[] { 2, 1, 5 });
            p.AddCoeff(-1.6360409591882500e+02, new int[] { 2, 3, 1 });
            p.AddCoeff(7.6348578095451800e+02, new int[] { 2, 3, 3 });
            p.AddCoeff(-6.8713720285906600e+02, new int[] { 2, 3, 5 });
            p.AddCoeff(-1.1452286714317800e+02, new int[] { 4, 1, 1 });
            p.AddCoeff(5.3444004666816300e+02, new int[] { 4, 1, 3 });
            p.AddCoeff(-4.8099604200134700e+02, new int[] { 4, 1, 5 });
            p.AddCoeff(1.9087144523863000e+02, new int[] { 4, 3, 1 });
            p.AddCoeff(-8.9073341111360500e+02, new int[] { 4, 3, 3 });
            p.AddCoeff(8.0166007000224400e+02, new int[] { 4, 3, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("880a7217-d8b3-42b3-b81a-7508a49f30ee"));
            OrthonormalPolynomials[414] = p;
            p.AddCoeff(5.0339926121581500e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-5.0339926121581500e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(5.8729913808511800e+00, new int[] { 0, 0, 4 });
            p.AddCoeff(-5.0339926121581500e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(5.0339926121581500e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(-5.8729913808511800e+01, new int[] { 0, 2, 4 });
            p.AddCoeff(5.8729913808511800e+00, new int[] { 0, 4, 0 });
            p.AddCoeff(-5.8729913808511800e+01, new int[] { 0, 4, 2 });
            p.AddCoeff(6.8518232776597100e+01, new int[] { 0, 4, 4 });
            p.AddCoeff(-5.0339926121581500e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(5.0339926121581500e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(-5.8729913808511800e+01, new int[] { 2, 0, 4 });
            p.AddCoeff(5.0339926121581500e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(-5.0339926121581500e+02, new int[] { 2, 2, 2 });
            p.AddCoeff(5.8729913808511800e+02, new int[] { 2, 2, 4 });
            p.AddCoeff(-5.8729913808511800e+01, new int[] { 2, 4, 0 });
            p.AddCoeff(5.8729913808511800e+02, new int[] { 2, 4, 2 });
            p.AddCoeff(-6.8518232776597100e+02, new int[] { 2, 4, 4 });
            p.AddCoeff(5.8729913808511800e+00, new int[] { 4, 0, 0 });
            p.AddCoeff(-5.8729913808511800e+01, new int[] { 4, 0, 2 });
            p.AddCoeff(6.8518232776597100e+01, new int[] { 4, 0, 4 });
            p.AddCoeff(-5.8729913808511800e+01, new int[] { 4, 2, 0 });
            p.AddCoeff(5.8729913808511800e+02, new int[] { 4, 2, 2 });
            p.AddCoeff(-6.8518232776597100e+02, new int[] { 4, 2, 4 });
            p.AddCoeff(6.8518232776597100e+01, new int[] { 4, 4, 0 });
            p.AddCoeff(-6.8518232776597100e+02, new int[] { 4, 4, 2 });
            p.AddCoeff(7.9937938239363300e+02, new int[] { 4, 4, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("338e7af9-4b1f-45d8-a131-09663064d6a0"));
            OrthonormalPolynomials[415] = p;
            p.AddCoeff(-9.8162457551295200e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(1.6360409591882500e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(4.5809146857271100e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(-7.6348578095451800e+01, new int[] { 0, 3, 3 });
            p.AddCoeff(-4.1228232171544000e+01, new int[] { 0, 5, 1 });
            p.AddCoeff(6.8713720285906600e+01, new int[] { 0, 5, 3 });
            p.AddCoeff(9.8162457551295200e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(-1.6360409591882500e+02, new int[] { 2, 1, 3 });
            p.AddCoeff(-4.5809146857271100e+02, new int[] { 2, 3, 1 });
            p.AddCoeff(7.6348578095451800e+02, new int[] { 2, 3, 3 });
            p.AddCoeff(4.1228232171544000e+02, new int[] { 2, 5, 1 });
            p.AddCoeff(-6.8713720285906600e+02, new int[] { 2, 5, 3 });
            p.AddCoeff(-1.1452286714317800e+02, new int[] { 4, 1, 1 });
            p.AddCoeff(1.9087144523863000e+02, new int[] { 4, 1, 3 });
            p.AddCoeff(5.3444004666816300e+02, new int[] { 4, 3, 1 });
            p.AddCoeff(-8.9073341111360500e+02, new int[] { 4, 3, 3 });
            p.AddCoeff(-4.8099604200134700e+02, new int[] { 4, 5, 1 });
            p.AddCoeff(8.0166007000224400e+02, new int[] { 4, 5, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("0e1eb7b4-35cd-4aba-aeed-04fec0ee8d9f"));
            OrthonormalPolynomials[416] = p;
            p.AddCoeff(5.0105365360802100e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-1.5031609608240600e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-1.0522126725768400e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(3.1566380177305300e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(3.1566380177305300e+01, new int[] { 0, 4, 0 });
            p.AddCoeff(-9.4699140531915900e+01, new int[] { 0, 4, 2 });
            p.AddCoeff(-2.3148678796690500e+01, new int[] { 0, 6, 0 });
            p.AddCoeff(6.9446036390071600e+01, new int[] { 0, 6, 2 });
            p.AddCoeff(-5.0105365360802100e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(1.5031609608240600e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(1.0522126725768400e+02, new int[] { 2, 2, 0 });
            p.AddCoeff(-3.1566380177305300e+02, new int[] { 2, 2, 2 });
            p.AddCoeff(-3.1566380177305300e+02, new int[] { 2, 4, 0 });
            p.AddCoeff(9.4699140531915900e+02, new int[] { 2, 4, 2 });
            p.AddCoeff(2.3148678796690500e+02, new int[] { 2, 6, 0 });
            p.AddCoeff(-6.9446036390071600e+02, new int[] { 2, 6, 2 });
            p.AddCoeff(5.8456259587602400e+00, new int[] { 4, 0, 0 });
            p.AddCoeff(-1.7536877876280700e+01, new int[] { 4, 0, 2 });
            p.AddCoeff(-1.2275814513396500e+02, new int[] { 4, 2, 0 });
            p.AddCoeff(3.6827443540189500e+02, new int[] { 4, 2, 2 });
            p.AddCoeff(3.6827443540189500e+02, new int[] { 4, 4, 0 });
            p.AddCoeff(-1.1048233062056900e+03, new int[] { 4, 4, 2 });
            p.AddCoeff(-2.7006791929472300e+02, new int[] { 4, 6, 0 });
            p.AddCoeff(8.1020375788416900e+02, new int[] { 4, 6, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("ac03749c-8005-4bc0-a775-6c2fd4cedf0a"));
            OrthonormalPolynomials[417] = p;
            p.AddCoeff(-5.8366257594904700e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(5.2529631835414200e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(-1.1556519003791100e+02, new int[] { 0, 5, 1 });
            p.AddCoeff(7.1540355737754600e+01, new int[] { 0, 7, 1 });
            p.AddCoeff(5.8366257594904700e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(-5.2529631835414200e+02, new int[] { 2, 3, 1 });
            p.AddCoeff(1.1556519003791100e+03, new int[] { 2, 5, 1 });
            p.AddCoeff(-7.1540355737754600e+02, new int[] { 2, 7, 1 });
            p.AddCoeff(-6.8093967194055400e+01, new int[] { 4, 1, 1 });
            p.AddCoeff(6.1284570474649900e+02, new int[] { 4, 3, 1 });
            p.AddCoeff(-1.3482605504423000e+03, new int[] { 4, 5, 1 });
            p.AddCoeff(8.3463748360713700e+02, new int[] { 4, 7, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7bfa9253-3565-44a1-9af3-b248244fc211"));
            OrthonormalPolynomials[418] = p;
            p.AddCoeff(4.4842525558502700e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-1.6143309201061000e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(8.8788200605835300e+01, new int[] { 0, 4, 0 });
            p.AddCoeff(-1.5389954771678100e+02, new int[] { 0, 6, 0 });
            p.AddCoeff(8.2446186276847100e+01, new int[] { 0, 8, 0 });
            p.AddCoeff(-4.4842525558502700e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(1.6143309201061000e+02, new int[] { 2, 2, 0 });
            p.AddCoeff(-8.8788200605835300e+02, new int[] { 2, 4, 0 });
            p.AddCoeff(1.5389954771678100e+03, new int[] { 2, 6, 0 });
            p.AddCoeff(-8.2446186276847100e+02, new int[] { 2, 8, 0 });
            p.AddCoeff(5.2316279818253100e+00, new int[] { 4, 0, 0 });
            p.AddCoeff(-1.8833860734571100e+02, new int[] { 4, 2, 0 });
            p.AddCoeff(1.0358623404014100e+03, new int[] { 4, 4, 0 });
            p.AddCoeff(-1.7954947233624500e+03, new int[] { 4, 6, 0 });
            p.AddCoeff(9.6187217322988200e+02, new int[] { 4, 8, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("0fa5296c-684f-42c3-8f8e-6392eed0d146"));
            OrthonormalPolynomials[419] = p;
            p.AddCoeff(-1.8627145733216900e+01, new int[] { 1, 0, 1 });
            p.AddCoeff(1.6764431159895200e+02, new int[] { 1, 0, 3 });
            p.AddCoeff(-3.6881748551769500e+02, new int[] { 1, 0, 5 });
            p.AddCoeff(2.2831558627285900e+02, new int[] { 1, 0, 7 });
            p.AddCoeff(8.6926680088345700e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(-7.8234012079511100e+02, new int[] { 3, 0, 3 });
            p.AddCoeff(1.7211482657492500e+03, new int[] { 3, 0, 5 });
            p.AddCoeff(-1.0654727359400100e+03, new int[] { 3, 0, 7 });
            p.AddCoeff(-7.8234012079511100e+01, new int[] { 5, 0, 1 });
            p.AddCoeff(7.0410610871560000e+02, new int[] { 5, 0, 3 });
            p.AddCoeff(-1.5490334391743200e+03, new int[] { 5, 0, 5 });
            p.AddCoeff(9.5892546234600800e+02, new int[] { 5, 0, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("54a19f14-55b9-44df-8d90-eddc551e3285"));
            OrthonormalPolynomials[420] = p;
            p.AddCoeff(-4.2907671433232800e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(9.0106110009788900e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(-2.7031833002936700e+02, new int[] { 1, 1, 4 });
            p.AddCoeff(1.9823344202153600e+02, new int[] { 1, 1, 6 });
            p.AddCoeff(2.0023580002175300e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(-4.2049518004568200e+02, new int[] { 3, 1, 2 });
            p.AddCoeff(1.2614855401370400e+03, new int[] { 3, 1, 4 });
            p.AddCoeff(-9.2508939610049900e+02, new int[] { 3, 1, 6 });
            p.AddCoeff(-1.8021222001957800e+01, new int[] { 5, 1, 0 });
            p.AddCoeff(3.7844566204111300e+02, new int[] { 5, 1, 2 });
            p.AddCoeff(-1.1353369861233400e+03, new int[] { 5, 1, 4 });
            p.AddCoeff(8.3258045649044900e+02, new int[] { 5, 1, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("0b583e03-d62b-43a7-b0d6-6da9f884516c"));
            OrthonormalPolynomials[421] = p;
            p.AddCoeff(-1.5286400798665500e+01, new int[] { 1, 0, 1 });
            p.AddCoeff(7.1336537060439000e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(-6.4202883354395100e+01, new int[] { 1, 0, 5 });
            p.AddCoeff(4.5859202395996500e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(-2.1400961118131700e+02, new int[] { 1, 2, 3 });
            p.AddCoeff(1.9260865006318500e+02, new int[] { 1, 2, 5 });
            p.AddCoeff(7.1336537060439000e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(-3.3290383961538200e+02, new int[] { 3, 0, 3 });
            p.AddCoeff(2.9961345565384400e+02, new int[] { 3, 0, 5 });
            p.AddCoeff(-2.1400961118131700e+02, new int[] { 3, 2, 1 });
            p.AddCoeff(9.9871151884614700e+02, new int[] { 3, 2, 3 });
            p.AddCoeff(-8.9884036696153200e+02, new int[] { 3, 2, 5 });
            p.AddCoeff(-6.4202883354395100e+01, new int[] { 5, 0, 1 });
            p.AddCoeff(2.9961345565384400e+02, new int[] { 5, 0, 3 });
            p.AddCoeff(-2.6965211008846000e+02, new int[] { 5, 0, 5 });
            p.AddCoeff(1.9260865006318500e+02, new int[] { 5, 2, 1 });
            p.AddCoeff(-8.9884036696153200e+02, new int[] { 5, 2, 3 });
            p.AddCoeff(8.0895633026537900e+02, new int[] { 5, 2, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("fe219bbd-fd8d-42f9-87fb-20c47da2ead1"));
            OrthonormalPolynomials[422] = p;
            p.AddCoeff(-9.8162457551295200e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(9.8162457551295200e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(-1.1452286714317800e+02, new int[] { 1, 1, 4 });
            p.AddCoeff(1.6360409591882500e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(-1.6360409591882500e+02, new int[] { 1, 3, 2 });
            p.AddCoeff(1.9087144523863000e+02, new int[] { 1, 3, 4 });
            p.AddCoeff(4.5809146857271100e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(-4.5809146857271100e+02, new int[] { 3, 1, 2 });
            p.AddCoeff(5.3444004666816300e+02, new int[] { 3, 1, 4 });
            p.AddCoeff(-7.6348578095451800e+01, new int[] { 3, 3, 0 });
            p.AddCoeff(7.6348578095451800e+02, new int[] { 3, 3, 2 });
            p.AddCoeff(-8.9073341111360500e+02, new int[] { 3, 3, 4 });
            p.AddCoeff(-4.1228232171544000e+01, new int[] { 5, 1, 0 });
            p.AddCoeff(4.1228232171544000e+02, new int[] { 5, 1, 2 });
            p.AddCoeff(-4.8099604200134700e+02, new int[] { 5, 1, 4 });
            p.AddCoeff(6.8713720285906600e+01, new int[] { 5, 3, 0 });
            p.AddCoeff(-6.8713720285906600e+02, new int[] { 5, 3, 2 });
            p.AddCoeff(8.0166007000224400e+02, new int[] { 5, 3, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e16d8061-c202-40cd-a304-438310d0c086"));
            OrthonormalPolynomials[423] = p;
            p.AddCoeff(-9.8162457551295200e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(1.6360409591882500e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(9.8162457551295200e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(-1.6360409591882500e+02, new int[] { 1, 2, 3 });
            p.AddCoeff(-1.1452286714317800e+02, new int[] { 1, 4, 1 });
            p.AddCoeff(1.9087144523863000e+02, new int[] { 1, 4, 3 });
            p.AddCoeff(4.5809146857271100e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(-7.6348578095451800e+01, new int[] { 3, 0, 3 });
            p.AddCoeff(-4.5809146857271100e+02, new int[] { 3, 2, 1 });
            p.AddCoeff(7.6348578095451800e+02, new int[] { 3, 2, 3 });
            p.AddCoeff(5.3444004666816300e+02, new int[] { 3, 4, 1 });
            p.AddCoeff(-8.9073341111360500e+02, new int[] { 3, 4, 3 });
            p.AddCoeff(-4.1228232171544000e+01, new int[] { 5, 0, 1 });
            p.AddCoeff(6.8713720285906600e+01, new int[] { 5, 0, 3 });
            p.AddCoeff(4.1228232171544000e+02, new int[] { 5, 2, 1 });
            p.AddCoeff(-6.8713720285906600e+02, new int[] { 5, 2, 3 });
            p.AddCoeff(-4.8099604200134700e+02, new int[] { 5, 4, 1 });
            p.AddCoeff(8.0166007000224400e+02, new int[] { 5, 4, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("72710c00-c4c0-4255-ba28-54069e357a76"));
            OrthonormalPolynomials[424] = p;
            p.AddCoeff(-1.5286400798665500e+01, new int[] { 1, 1, 0 });
            p.AddCoeff(4.5859202395996500e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(7.1336537060439000e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(-2.1400961118131700e+02, new int[] { 1, 3, 2 });
            p.AddCoeff(-6.4202883354395100e+01, new int[] { 1, 5, 0 });
            p.AddCoeff(1.9260865006318500e+02, new int[] { 1, 5, 2 });
            p.AddCoeff(7.1336537060439000e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(-2.1400961118131700e+02, new int[] { 3, 1, 2 });
            p.AddCoeff(-3.3290383961538200e+02, new int[] { 3, 3, 0 });
            p.AddCoeff(9.9871151884614700e+02, new int[] { 3, 3, 2 });
            p.AddCoeff(2.9961345565384400e+02, new int[] { 3, 5, 0 });
            p.AddCoeff(-8.9884036696153200e+02, new int[] { 3, 5, 2 });
            p.AddCoeff(-6.4202883354395100e+01, new int[] { 5, 1, 0 });
            p.AddCoeff(1.9260865006318500e+02, new int[] { 5, 1, 2 });
            p.AddCoeff(2.9961345565384400e+02, new int[] { 5, 3, 0 });
            p.AddCoeff(-8.9884036696153200e+02, new int[] { 5, 3, 2 });
            p.AddCoeff(-2.6965211008846000e+02, new int[] { 5, 5, 0 });
            p.AddCoeff(8.0895633026537900e+02, new int[] { 5, 5, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("343e0969-6255-4cd1-8cea-4f88cb34799c"));
            OrthonormalPolynomials[425] = p;
            p.AddCoeff(-4.2907671433232800e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(9.0106110009788900e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(-2.7031833002936700e+02, new int[] { 1, 4, 1 });
            p.AddCoeff(1.9823344202153600e+02, new int[] { 1, 6, 1 });
            p.AddCoeff(2.0023580002175300e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(-4.2049518004568200e+02, new int[] { 3, 2, 1 });
            p.AddCoeff(1.2614855401370400e+03, new int[] { 3, 4, 1 });
            p.AddCoeff(-9.2508939610049900e+02, new int[] { 3, 6, 1 });
            p.AddCoeff(-1.8021222001957800e+01, new int[] { 5, 0, 1 });
            p.AddCoeff(3.7844566204111300e+02, new int[] { 5, 2, 1 });
            p.AddCoeff(-1.1353369861233400e+03, new int[] { 5, 4, 1 });
            p.AddCoeff(8.3258045649044900e+02, new int[] { 5, 6, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7e785497-7392-430c-854f-f320762b23e1"));
            OrthonormalPolynomials[426] = p;
            p.AddCoeff(-1.8627145733216900e+01, new int[] { 1, 1, 0 });
            p.AddCoeff(1.6764431159895200e+02, new int[] { 1, 3, 0 });
            p.AddCoeff(-3.6881748551769500e+02, new int[] { 1, 5, 0 });
            p.AddCoeff(2.2831558627285900e+02, new int[] { 1, 7, 0 });
            p.AddCoeff(8.6926680088345700e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(-7.8234012079511100e+02, new int[] { 3, 3, 0 });
            p.AddCoeff(1.7211482657492500e+03, new int[] { 3, 5, 0 });
            p.AddCoeff(-1.0654727359400100e+03, new int[] { 3, 7, 0 });
            p.AddCoeff(-7.8234012079511100e+01, new int[] { 5, 1, 0 });
            p.AddCoeff(7.0410610871560000e+02, new int[] { 5, 3, 0 });
            p.AddCoeff(-1.5490334391743200e+03, new int[] { 5, 5, 0 });
            p.AddCoeff(9.5892546234600800e+02, new int[] { 5, 7, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("341805cb-d9b7-414a-9404-30e7e6083858"));
            OrthonormalPolynomials[427] = p;
            p.AddCoeff(4.4884707790161900e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-9.4257886359339900e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(2.8277365907802000e+01, new int[] { 0, 0, 4 });
            p.AddCoeff(-2.0736734999054800e+01, new int[] { 0, 0, 6 });
            p.AddCoeff(-9.4257886359339900e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(1.9794156135461400e+02, new int[] { 2, 0, 2 });
            p.AddCoeff(-5.9382468406384100e+02, new int[] { 2, 0, 4 });
            p.AddCoeff(4.3547143498015000e+02, new int[] { 2, 0, 6 });
            p.AddCoeff(2.8277365907802000e+01, new int[] { 4, 0, 0 });
            p.AddCoeff(-5.9382468406384100e+02, new int[] { 4, 0, 2 });
            p.AddCoeff(1.7814740521915200e+03, new int[] { 4, 0, 4 });
            p.AddCoeff(-1.3064143049404500e+03, new int[] { 4, 0, 6 });
            p.AddCoeff(-2.0736734999054800e+01, new int[] { 6, 0, 0 });
            p.AddCoeff(4.3547143498015000e+02, new int[] { 6, 0, 2 });
            p.AddCoeff(-1.3064143049404500e+03, new int[] { 6, 0, 4 });
            p.AddCoeff(9.5803715695633100e+02, new int[] { 6, 0, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("be7a3a30-94cf-4c07-a6c0-8f1d2ae16d90"));
            OrthonormalPolynomials[428] = p;
            p.AddCoeff(-4.2907671433232800e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(2.0023580002175300e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(-1.8021222001957800e+01, new int[] { 0, 1, 5 });
            p.AddCoeff(9.0106110009788900e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(-4.2049518004568200e+02, new int[] { 2, 1, 3 });
            p.AddCoeff(3.7844566204111300e+02, new int[] { 2, 1, 5 });
            p.AddCoeff(-2.7031833002936700e+02, new int[] { 4, 1, 1 });
            p.AddCoeff(1.2614855401370400e+03, new int[] { 4, 1, 3 });
            p.AddCoeff(-1.1353369861233400e+03, new int[] { 4, 1, 5 });
            p.AddCoeff(1.9823344202153600e+02, new int[] { 6, 1, 1 });
            p.AddCoeff(-9.2508939610049900e+02, new int[] { 6, 1, 3 });
            p.AddCoeff(8.3258045649044900e+02, new int[] { 6, 1, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("18594afb-d5a7-4f98-a89d-00091ca6b719"));
            OrthonormalPolynomials[429] = p;
            p.AddCoeff(5.0105365360802100e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-5.0105365360802100e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(5.8456259587602400e+00, new int[] { 0, 0, 4 });
            p.AddCoeff(-1.5031609608240600e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(1.5031609608240600e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(-1.7536877876280700e+01, new int[] { 0, 2, 4 });
            p.AddCoeff(-1.0522126725768400e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(1.0522126725768400e+02, new int[] { 2, 0, 2 });
            p.AddCoeff(-1.2275814513396500e+02, new int[] { 2, 0, 4 });
            p.AddCoeff(3.1566380177305300e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(-3.1566380177305300e+02, new int[] { 2, 2, 2 });
            p.AddCoeff(3.6827443540189500e+02, new int[] { 2, 2, 4 });
            p.AddCoeff(3.1566380177305300e+01, new int[] { 4, 0, 0 });
            p.AddCoeff(-3.1566380177305300e+02, new int[] { 4, 0, 2 });
            p.AddCoeff(3.6827443540189500e+02, new int[] { 4, 0, 4 });
            p.AddCoeff(-9.4699140531915900e+01, new int[] { 4, 2, 0 });
            p.AddCoeff(9.4699140531915900e+02, new int[] { 4, 2, 2 });
            p.AddCoeff(-1.1048233062056900e+03, new int[] { 4, 2, 4 });
            p.AddCoeff(-2.3148678796690500e+01, new int[] { 6, 0, 0 });
            p.AddCoeff(2.3148678796690500e+02, new int[] { 6, 0, 2 });
            p.AddCoeff(-2.7006791929472300e+02, new int[] { 6, 0, 4 });
            p.AddCoeff(6.9446036390071600e+01, new int[] { 6, 2, 0 });
            p.AddCoeff(-6.9446036390071600e+02, new int[] { 6, 2, 2 });
            p.AddCoeff(8.1020375788416900e+02, new int[] { 6, 2, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("98ee6e9b-eea0-486e-9c82-2979100e9e26"));
            OrthonormalPolynomials[430] = p;
            p.AddCoeff(-6.2741841671161200e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(1.0456973611860200e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(1.0456973611860200e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(-1.7428289353100300e+01, new int[] { 0, 3, 3 });
            p.AddCoeff(1.3175786750943800e+02, new int[] { 2, 1, 1 });
            p.AddCoeff(-2.1959644584906400e+02, new int[] { 2, 1, 3 });
            p.AddCoeff(-2.1959644584906400e+02, new int[] { 2, 3, 1 });
            p.AddCoeff(3.6599407641510700e+02, new int[] { 2, 3, 3 });
            p.AddCoeff(-3.9527360252831500e+02, new int[] { 4, 1, 1 });
            p.AddCoeff(6.5878933754719200e+02, new int[] { 4, 1, 3 });
            p.AddCoeff(6.5878933754719200e+02, new int[] { 4, 3, 1 });
            p.AddCoeff(-1.0979822292453200e+03, new int[] { 4, 3, 3 });
            p.AddCoeff(2.8986730852076500e+02, new int[] { 6, 1, 1 });
            p.AddCoeff(-4.8311218086794100e+02, new int[] { 6, 1, 3 });
            p.AddCoeff(-4.8311218086794100e+02, new int[] { 6, 3, 1 });
            p.AddCoeff(8.0518696811323500e+02, new int[] { 6, 3, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("957f2c15-48b6-4c7c-8d1e-ad6ea30b0d45"));
            OrthonormalPolynomials[431] = p;
            p.AddCoeff(5.0105365360802100e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-1.5031609608240600e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-5.0105365360802100e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(1.5031609608240600e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(5.8456259587602400e+00, new int[] { 0, 4, 0 });
            p.AddCoeff(-1.7536877876280700e+01, new int[] { 0, 4, 2 });
            p.AddCoeff(-1.0522126725768400e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(3.1566380177305300e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(1.0522126725768400e+02, new int[] { 2, 2, 0 });
            p.AddCoeff(-3.1566380177305300e+02, new int[] { 2, 2, 2 });
            p.AddCoeff(-1.2275814513396500e+02, new int[] { 2, 4, 0 });
            p.AddCoeff(3.6827443540189500e+02, new int[] { 2, 4, 2 });
            p.AddCoeff(3.1566380177305300e+01, new int[] { 4, 0, 0 });
            p.AddCoeff(-9.4699140531915900e+01, new int[] { 4, 0, 2 });
            p.AddCoeff(-3.1566380177305300e+02, new int[] { 4, 2, 0 });
            p.AddCoeff(9.4699140531915900e+02, new int[] { 4, 2, 2 });
            p.AddCoeff(3.6827443540189500e+02, new int[] { 4, 4, 0 });
            p.AddCoeff(-1.1048233062056900e+03, new int[] { 4, 4, 2 });
            p.AddCoeff(-2.3148678796690500e+01, new int[] { 6, 0, 0 });
            p.AddCoeff(6.9446036390071600e+01, new int[] { 6, 0, 2 });
            p.AddCoeff(2.3148678796690500e+02, new int[] { 6, 2, 0 });
            p.AddCoeff(-6.9446036390071600e+02, new int[] { 6, 2, 2 });
            p.AddCoeff(-2.7006791929472300e+02, new int[] { 6, 4, 0 });
            p.AddCoeff(8.1020375788416900e+02, new int[] { 6, 4, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5ba3402f-871a-4ee2-acfb-bc6e9dfafbb9"));
            OrthonormalPolynomials[432] = p;
            p.AddCoeff(-4.2907671433232800e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(2.0023580002175300e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(-1.8021222001957800e+01, new int[] { 0, 5, 1 });
            p.AddCoeff(9.0106110009788900e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(-4.2049518004568200e+02, new int[] { 2, 3, 1 });
            p.AddCoeff(3.7844566204111300e+02, new int[] { 2, 5, 1 });
            p.AddCoeff(-2.7031833002936700e+02, new int[] { 4, 1, 1 });
            p.AddCoeff(1.2614855401370400e+03, new int[] { 4, 3, 1 });
            p.AddCoeff(-1.1353369861233400e+03, new int[] { 4, 5, 1 });
            p.AddCoeff(1.9823344202153600e+02, new int[] { 6, 1, 1 });
            p.AddCoeff(-9.2508939610049900e+02, new int[] { 6, 3, 1 });
            p.AddCoeff(8.3258045649044900e+02, new int[] { 6, 5, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e17df5f0-df3e-4899-8b4e-b7c0d9852913"));
            OrthonormalPolynomials[433] = p;
            p.AddCoeff(4.4884707790161900e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-9.4257886359339900e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(2.8277365907802000e+01, new int[] { 0, 4, 0 });
            p.AddCoeff(-2.0736734999054800e+01, new int[] { 0, 6, 0 });
            p.AddCoeff(-9.4257886359339900e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(1.9794156135461400e+02, new int[] { 2, 2, 0 });
            p.AddCoeff(-5.9382468406384100e+02, new int[] { 2, 4, 0 });
            p.AddCoeff(4.3547143498015000e+02, new int[] { 2, 6, 0 });
            p.AddCoeff(2.8277365907802000e+01, new int[] { 4, 0, 0 });
            p.AddCoeff(-5.9382468406384100e+02, new int[] { 4, 2, 0 });
            p.AddCoeff(1.7814740521915200e+03, new int[] { 4, 4, 0 });
            p.AddCoeff(-1.3064143049404500e+03, new int[] { 4, 6, 0 });
            p.AddCoeff(-2.0736734999054800e+01, new int[] { 6, 0, 0 });
            p.AddCoeff(4.3547143498015000e+02, new int[] { 6, 2, 0 });
            p.AddCoeff(-1.3064143049404500e+03, new int[] { 6, 4, 0 });
            p.AddCoeff(9.5803715695633100e+02, new int[] { 6, 6, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5b3a6bc0-2013-4e57-9e7a-2eb4871bc6b0"));
            OrthonormalPolynomials[434] = p;
            p.AddCoeff(-1.8627145733216900e+01, new int[] { 1, 0, 1 });
            p.AddCoeff(8.6926680088345700e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(-7.8234012079511100e+01, new int[] { 1, 0, 5 });
            p.AddCoeff(1.6764431159895200e+02, new int[] { 3, 0, 1 });
            p.AddCoeff(-7.8234012079511100e+02, new int[] { 3, 0, 3 });
            p.AddCoeff(7.0410610871560000e+02, new int[] { 3, 0, 5 });
            p.AddCoeff(-3.6881748551769500e+02, new int[] { 5, 0, 1 });
            p.AddCoeff(1.7211482657492500e+03, new int[] { 5, 0, 3 });
            p.AddCoeff(-1.5490334391743200e+03, new int[] { 5, 0, 5 });
            p.AddCoeff(2.2831558627285900e+02, new int[] { 7, 0, 1 });
            p.AddCoeff(-1.0654727359400100e+03, new int[] { 7, 0, 3 });
            p.AddCoeff(9.5892546234600800e+02, new int[] { 7, 0, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b5413e87-cd96-4ebd-ab00-c52cbefae550"));
            OrthonormalPolynomials[435] = p;
            p.AddCoeff(-5.8366257594904700e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(5.8366257594904700e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(-6.8093967194055400e+01, new int[] { 1, 1, 4 });
            p.AddCoeff(5.2529631835414200e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(-5.2529631835414200e+02, new int[] { 3, 1, 2 });
            p.AddCoeff(6.1284570474649900e+02, new int[] { 3, 1, 4 });
            p.AddCoeff(-1.1556519003791100e+02, new int[] { 5, 1, 0 });
            p.AddCoeff(1.1556519003791100e+03, new int[] { 5, 1, 2 });
            p.AddCoeff(-1.3482605504423000e+03, new int[] { 5, 1, 4 });
            p.AddCoeff(7.1540355737754600e+01, new int[] { 7, 1, 0 });
            p.AddCoeff(-7.1540355737754600e+02, new int[] { 7, 1, 2 });
            p.AddCoeff(8.3463748360713700e+02, new int[] { 7, 1, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b39c1dd2-a6b1-4faf-b7a9-f133acfbec44"));
            OrthonormalPolynomials[436] = p;
            p.AddCoeff(-1.3290581510406700e+01, new int[] { 1, 0, 1 });
            p.AddCoeff(2.2150969184011200e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(3.9871744531220200e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(-6.6452907552033700e+01, new int[] { 1, 2, 3 });
            p.AddCoeff(1.1961523359366100e+02, new int[] { 3, 0, 1 });
            p.AddCoeff(-1.9935872265610100e+02, new int[] { 3, 0, 3 });
            p.AddCoeff(-3.5884570078098200e+02, new int[] { 3, 2, 1 });
            p.AddCoeff(5.9807616796830300e+02, new int[] { 3, 2, 3 });
            p.AddCoeff(-2.6315351390605400e+02, new int[] { 5, 0, 1 });
            p.AddCoeff(4.3858918984342300e+02, new int[] { 5, 0, 3 });
            p.AddCoeff(7.8946054171816100e+02, new int[] { 5, 2, 1 });
            p.AddCoeff(-1.3157675695302700e+03, new int[] { 5, 2, 3 });
            p.AddCoeff(1.6290455622755700e+02, new int[] { 7, 0, 1 });
            p.AddCoeff(-2.7150759371259500e+02, new int[] { 7, 0, 3 });
            p.AddCoeff(-4.8871366868267100e+02, new int[] { 7, 2, 1 });
            p.AddCoeff(8.1452278113778500e+02, new int[] { 7, 2, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a2931403-6413-4565-82a3-afe77660b3a7"));
            OrthonormalPolynomials[437] = p;
            p.AddCoeff(-1.3290581510406700e+01, new int[] { 1, 1, 0 });
            p.AddCoeff(3.9871744531220200e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(2.2150969184011200e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(-6.6452907552033700e+01, new int[] { 1, 3, 2 });
            p.AddCoeff(1.1961523359366100e+02, new int[] { 3, 1, 0 });
            p.AddCoeff(-3.5884570078098200e+02, new int[] { 3, 1, 2 });
            p.AddCoeff(-1.9935872265610100e+02, new int[] { 3, 3, 0 });
            p.AddCoeff(5.9807616796830300e+02, new int[] { 3, 3, 2 });
            p.AddCoeff(-2.6315351390605400e+02, new int[] { 5, 1, 0 });
            p.AddCoeff(7.8946054171816100e+02, new int[] { 5, 1, 2 });
            p.AddCoeff(4.3858918984342300e+02, new int[] { 5, 3, 0 });
            p.AddCoeff(-1.3157675695302700e+03, new int[] { 5, 3, 2 });
            p.AddCoeff(1.6290455622755700e+02, new int[] { 7, 1, 0 });
            p.AddCoeff(-4.8871366868267100e+02, new int[] { 7, 1, 2 });
            p.AddCoeff(-2.7150759371259500e+02, new int[] { 7, 3, 0 });
            p.AddCoeff(8.1452278113778500e+02, new int[] { 7, 3, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7656fd84-87b6-4210-a8d6-1da8f278bb29"));
            OrthonormalPolynomials[438] = p;
            p.AddCoeff(-5.8366257594904700e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(5.8366257594904700e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(-6.8093967194055400e+01, new int[] { 1, 4, 1 });
            p.AddCoeff(5.2529631835414200e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(-5.2529631835414200e+02, new int[] { 3, 2, 1 });
            p.AddCoeff(6.1284570474649900e+02, new int[] { 3, 4, 1 });
            p.AddCoeff(-1.1556519003791100e+02, new int[] { 5, 0, 1 });
            p.AddCoeff(1.1556519003791100e+03, new int[] { 5, 2, 1 });
            p.AddCoeff(-1.3482605504423000e+03, new int[] { 5, 4, 1 });
            p.AddCoeff(7.1540355737754600e+01, new int[] { 7, 0, 1 });
            p.AddCoeff(-7.1540355737754600e+02, new int[] { 7, 2, 1 });
            p.AddCoeff(8.3463748360713700e+02, new int[] { 7, 4, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("57f6c223-57c9-4a8f-b7a8-ec390d0d115a"));
            OrthonormalPolynomials[439] = p;
            p.AddCoeff(-1.8627145733216900e+01, new int[] { 1, 1, 0 });
            p.AddCoeff(8.6926680088345700e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(-7.8234012079511100e+01, new int[] { 1, 5, 0 });
            p.AddCoeff(1.6764431159895200e+02, new int[] { 3, 1, 0 });
            p.AddCoeff(-7.8234012079511100e+02, new int[] { 3, 3, 0 });
            p.AddCoeff(7.0410610871560000e+02, new int[] { 3, 5, 0 });
            p.AddCoeff(-3.6881748551769500e+02, new int[] { 5, 1, 0 });
            p.AddCoeff(1.7211482657492500e+03, new int[] { 5, 3, 0 });
            p.AddCoeff(-1.5490334391743200e+03, new int[] { 5, 5, 0 });
            p.AddCoeff(2.2831558627285900e+02, new int[] { 7, 1, 0 });
            p.AddCoeff(-1.0654727359400100e+03, new int[] { 7, 3, 0 });
            p.AddCoeff(9.5892546234600800e+02, new int[] { 7, 5, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("afcf8c37-af69-4a0d-a2d0-0c4194e8bec7"));
            OrthonormalPolynomials[440] = p;
            p.AddCoeff(4.4842525558502700e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-4.4842525558502700e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(5.2316279818253100e+00, new int[] { 0, 0, 4 });
            p.AddCoeff(-1.6143309201061000e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(1.6143309201061000e+02, new int[] { 2, 0, 2 });
            p.AddCoeff(-1.8833860734571100e+02, new int[] { 2, 0, 4 });
            p.AddCoeff(8.8788200605835300e+01, new int[] { 4, 0, 0 });
            p.AddCoeff(-8.8788200605835300e+02, new int[] { 4, 0, 2 });
            p.AddCoeff(1.0358623404014100e+03, new int[] { 4, 0, 4 });
            p.AddCoeff(-1.5389954771678100e+02, new int[] { 6, 0, 0 });
            p.AddCoeff(1.5389954771678100e+03, new int[] { 6, 0, 2 });
            p.AddCoeff(-1.7954947233624500e+03, new int[] { 6, 0, 4 });
            p.AddCoeff(8.2446186276847100e+01, new int[] { 8, 0, 0 });
            p.AddCoeff(-8.2446186276847100e+02, new int[] { 8, 0, 2 });
            p.AddCoeff(9.6187217322988200e+02, new int[] { 8, 0, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("93673c80-433c-4756-8c7b-13ac8e94ce50"));
            OrthonormalPolynomials[441] = p;
            p.AddCoeff(-2.7399235696644100e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(4.5665392827740100e+00, new int[] { 0, 1, 3 });
            p.AddCoeff(9.8637248507918700e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(-1.6439541417986400e+02, new int[] { 2, 1, 3 });
            p.AddCoeff(-5.4250486679355300e+02, new int[] { 4, 1, 1 });
            p.AddCoeff(9.0417477798925500e+02, new int[] { 4, 1, 3 });
            p.AddCoeff(9.4034176910882500e+02, new int[] { 6, 1, 1 });
            p.AddCoeff(-1.5672362818480400e+03, new int[] { 6, 1, 3 });
            p.AddCoeff(-5.0375451916544200e+02, new int[] { 8, 1, 1 });
            p.AddCoeff(8.3959086527573700e+02, new int[] { 8, 1, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("4b2a8a41-4b47-4acc-a810-b3af974fbcdc"));
            OrthonormalPolynomials[442] = p;
            p.AddCoeff(4.9825028398336300e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-1.4947508519500900e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-1.4947508519500900e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(4.4842525558502700e+00, new int[] { 0, 2, 2 });
            p.AddCoeff(-1.7937010223401100e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(5.3811030670203200e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(5.3811030670203200e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(-1.6143309201061000e+02, new int[] { 2, 2, 2 });
            p.AddCoeff(9.8653556228705900e+01, new int[] { 4, 0, 0 });
            p.AddCoeff(-2.9596066868611800e+02, new int[] { 4, 0, 2 });
            p.AddCoeff(-2.9596066868611800e+02, new int[] { 4, 2, 0 });
            p.AddCoeff(8.8788200605835300e+02, new int[] { 4, 2, 2 });
            p.AddCoeff(-1.7099949746309000e+02, new int[] { 6, 0, 0 });
            p.AddCoeff(5.1299849238927100e+02, new int[] { 6, 0, 2 });
            p.AddCoeff(5.1299849238927100e+02, new int[] { 6, 2, 0 });
            p.AddCoeff(-1.5389954771678100e+03, new int[] { 6, 2, 2 });
            p.AddCoeff(9.1606873640941200e+01, new int[] { 8, 0, 0 });
            p.AddCoeff(-2.7482062092282400e+02, new int[] { 8, 0, 2 });
            p.AddCoeff(-2.7482062092282400e+02, new int[] { 8, 2, 0 });
            p.AddCoeff(8.2446186276847100e+02, new int[] { 8, 2, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("9f78f951-9d8f-404d-bf7e-74290a6af662"));
            OrthonormalPolynomials[443] = p;
            p.AddCoeff(-2.7399235696644100e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(4.5665392827740100e+00, new int[] { 0, 3, 1 });
            p.AddCoeff(9.8637248507918700e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(-1.6439541417986400e+02, new int[] { 2, 3, 1 });
            p.AddCoeff(-5.4250486679355300e+02, new int[] { 4, 1, 1 });
            p.AddCoeff(9.0417477798925500e+02, new int[] { 4, 3, 1 });
            p.AddCoeff(9.4034176910882500e+02, new int[] { 6, 1, 1 });
            p.AddCoeff(-1.5672362818480400e+03, new int[] { 6, 3, 1 });
            p.AddCoeff(-5.0375451916544200e+02, new int[] { 8, 1, 1 });
            p.AddCoeff(8.3959086527573700e+02, new int[] { 8, 3, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("57225f59-2d00-46bd-98c4-432658f6ab7d"));
            OrthonormalPolynomials[444] = p;
            p.AddCoeff(4.4842525558502700e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-4.4842525558502700e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(5.2316279818253100e+00, new int[] { 0, 4, 0 });
            p.AddCoeff(-1.6143309201061000e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(1.6143309201061000e+02, new int[] { 2, 2, 0 });
            p.AddCoeff(-1.8833860734571100e+02, new int[] { 2, 4, 0 });
            p.AddCoeff(8.8788200605835300e+01, new int[] { 4, 0, 0 });
            p.AddCoeff(-8.8788200605835300e+02, new int[] { 4, 2, 0 });
            p.AddCoeff(1.0358623404014100e+03, new int[] { 4, 4, 0 });
            p.AddCoeff(-1.5389954771678100e+02, new int[] { 6, 0, 0 });
            p.AddCoeff(1.5389954771678100e+03, new int[] { 6, 2, 0 });
            p.AddCoeff(-1.7954947233624500e+03, new int[] { 6, 4, 0 });
            p.AddCoeff(8.2446186276847100e+01, new int[] { 8, 0, 0 });
            p.AddCoeff(-8.2446186276847100e+02, new int[] { 8, 2, 0 });
            p.AddCoeff(9.6187217322988200e+02, new int[] { 8, 4, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e7b406ad-4789-4f7c-9874-5d49aa568934"));
            OrthonormalPolynomials[445] = p;
            p.AddCoeff(-1.5051253492806200e+01, new int[] { 1, 0, 1 });
            p.AddCoeff(2.5085422488010400e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(2.2075171789449200e+02, new int[] { 3, 0, 1 });
            p.AddCoeff(-3.6791952982415300e+02, new int[] { 3, 0, 3 });
            p.AddCoeff(-8.6093169978851700e+02, new int[] { 5, 0, 1 });
            p.AddCoeff(1.4348861663142000e+03, new int[] { 5, 0, 3 });
            p.AddCoeff(1.2299024282693100e+03, new int[] { 7, 0, 1 });
            p.AddCoeff(-2.0498373804488500e+03, new int[] { 7, 0, 3 });
            p.AddCoeff(-5.8078725779384100e+02, new int[] { 9, 0, 1 });
            p.AddCoeff(9.6797876298973500e+02, new int[] { 9, 0, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d4ded5d4-534a-4912-a2e5-9d57cc93202b"));
            OrthonormalPolynomials[446] = p;
            p.AddCoeff(-7.3442596907982200e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(2.2032779072394700e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(1.0771580879837400e+02, new int[] { 3, 1, 0 });
            p.AddCoeff(-3.2314742639512200e+02, new int[] { 3, 1, 2 });
            p.AddCoeff(-4.2009165431365800e+02, new int[] { 5, 1, 0 });
            p.AddCoeff(1.2602749629409700e+03, new int[] { 5, 1, 2 });
            p.AddCoeff(6.0013093473379800e+02, new int[] { 7, 1, 0 });
            p.AddCoeff(-1.8003928042013900e+03, new int[] { 7, 1, 2 });
            p.AddCoeff(-2.8339516362429300e+02, new int[] { 9, 1, 0 });
            p.AddCoeff(8.5018549087288000e+02, new int[] { 9, 1, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("ad09d6cc-3caf-4914-aab5-a1bdc2fb845b"));
            OrthonormalPolynomials[447] = p;
            p.AddCoeff(-7.3442596907982200e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(2.2032779072394700e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(1.0771580879837400e+02, new int[] { 3, 0, 1 });
            p.AddCoeff(-3.2314742639512200e+02, new int[] { 3, 2, 1 });
            p.AddCoeff(-4.2009165431365800e+02, new int[] { 5, 0, 1 });
            p.AddCoeff(1.2602749629409700e+03, new int[] { 5, 2, 1 });
            p.AddCoeff(6.0013093473379800e+02, new int[] { 7, 0, 1 });
            p.AddCoeff(-1.8003928042013900e+03, new int[] { 7, 2, 1 });
            p.AddCoeff(-2.8339516362429300e+02, new int[] { 9, 0, 1 });
            p.AddCoeff(8.5018549087288000e+02, new int[] { 9, 2, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("824cdd31-2e34-4302-b92a-fc8b0cb64bd6"));
            OrthonormalPolynomials[448] = p;
            p.AddCoeff(-1.5051253492806200e+01, new int[] { 1, 1, 0 });
            p.AddCoeff(2.5085422488010400e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(2.2075171789449200e+02, new int[] { 3, 1, 0 });
            p.AddCoeff(-3.6791952982415300e+02, new int[] { 3, 3, 0 });
            p.AddCoeff(-8.6093169978851700e+02, new int[] { 5, 1, 0 });
            p.AddCoeff(1.4348861663142000e+03, new int[] { 5, 3, 0 });
            p.AddCoeff(1.2299024282693100e+03, new int[] { 7, 1, 0 });
            p.AddCoeff(-2.0498373804488500e+03, new int[] { 7, 3, 0 });
            p.AddCoeff(-5.8078725779384100e+02, new int[] { 9, 1, 0 });
            p.AddCoeff(9.6797876298973500e+02, new int[] { 9, 3, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3781a894-cdd7-4032-918a-4797d61da601"));
            OrthonormalPolynomials[449] = p;
            p.AddCoeff(4.4577965576656800e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-1.3373389672997100e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-2.4517881067161300e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(7.3553643201483800e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(2.1248830258206400e+02, new int[] { 4, 0, 0 });
            p.AddCoeff(-6.3746490774619300e+02, new int[] { 4, 0, 2 });
            p.AddCoeff(-6.3746490774619300e+02, new int[] { 6, 0, 0 });
            p.AddCoeff(1.9123947232385800e+03, new int[] { 6, 0, 2 });
            p.AddCoeff(7.7406453083466300e+02, new int[] { 8, 0, 0 });
            p.AddCoeff(-2.3221935925039900e+03, new int[] { 8, 0, 2 });
            p.AddCoeff(-3.2682724635241300e+02, new int[] { 10, 0, 0 });
            p.AddCoeff(9.8048173905724000e+02, new int[] { 10, 0, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5f20c6d8-b40f-4640-a89d-0cebc84f6c2a"));
            OrthonormalPolynomials[450] = p;
            p.AddCoeff(-1.1961523359366100e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(6.5788378476513400e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(-5.7016594679644900e+02, new int[] { 4, 1, 1 });
            p.AddCoeff(1.7104978403893500e+03, new int[] { 6, 1, 1 });
            p.AddCoeff(-2.0770330919013500e+03, new int[] { 8, 1, 1 });
            p.AddCoeff(8.7696952769168200e+02, new int[] { 10, 1, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("ddfc939b-4bd6-44a3-97ca-1942b71e6d3d"));
            OrthonormalPolynomials[451] = p;
            p.AddCoeff(4.4577965576656800e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-1.3373389672997100e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-2.4517881067161300e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(7.3553643201483800e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(2.1248830258206400e+02, new int[] { 4, 0, 0 });
            p.AddCoeff(-6.3746490774619300e+02, new int[] { 4, 2, 0 });
            p.AddCoeff(-6.3746490774619300e+02, new int[] { 6, 0, 0 });
            p.AddCoeff(1.9123947232385800e+03, new int[] { 6, 2, 0 });
            p.AddCoeff(7.7406453083466300e+02, new int[] { 8, 0, 0 });
            p.AddCoeff(-2.3221935925039900e+03, new int[] { 8, 2, 0 });
            p.AddCoeff(-3.2682724635241300e+02, new int[] { 10, 0, 0 });
            p.AddCoeff(9.8048173905724000e+02, new int[] { 10, 2, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d22ec12d-2d1d-446b-8e2c-d075513ad3ab"));
            OrthonormalPolynomials[452] = p;
            p.AddCoeff(-7.9501042053302700e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(1.7225225778215600e+02, new int[] { 3, 0, 1 });
            p.AddCoeff(-1.0335135466929300e+03, new int[] { 5, 0, 1 });
            p.AddCoeff(2.5099614705399800e+03, new int[] { 7, 0, 1 });
            p.AddCoeff(-2.6494037744588700e+03, new int[] { 9, 0, 1 });
            p.AddCoeff(1.0115905320661200e+03, new int[] { 11, 0, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("61432642-b95f-4fad-8b40-122639abbdf4"));
            OrthonormalPolynomials[453] = p;
            p.AddCoeff(-7.9501042053302700e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(1.7225225778215600e+02, new int[] { 3, 1, 0 });
            p.AddCoeff(-1.0335135466929300e+03, new int[] { 5, 1, 0 });
            p.AddCoeff(2.5099614705399800e+03, new int[] { 7, 1, 0 });
            p.AddCoeff(-2.6494037744588700e+03, new int[] { 9, 1, 0 });
            p.AddCoeff(1.0115905320661200e+03, new int[] { 11, 1, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c5db8a84-775d-4063-a370-af194c816d93"));
            OrthonormalPolynomials[454] = p;
            p.AddCoeff(3.9878336536643800e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(-3.1105102498582200e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(3.8881378123227700e+02, new int[] { 4, 0, 0 });
            p.AddCoeff(-1.7626224749196600e+03, new int[] { 6, 0, 0 });
            p.AddCoeff(3.5881957525150100e+03, new int[] { 8, 0, 0 });
            p.AddCoeff(-3.3489827023473500e+03, new int[] { 10, 0, 0 });
            p.AddCoeff(1.1670697296058900e+03, new int[] { 12, 0, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("8614a656-ecee-4d58-bf59-1eb4e5a57402"));
            OrthonormalPolynomials[455] = p;
            p.AddCoeff(5.3875617902181800e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-1.6162685370654500e+02, new int[] { 0, 0, 3 });
            p.AddCoeff(1.3738282565056400e+03, new int[] { 0, 0, 5 });
            p.AddCoeff(-4.9719498806870600e+03, new int[] { 0, 0, 7 });
            p.AddCoeff(8.7009122912023700e+03, new int[] { 0, 0, 9 });
            p.AddCoeff(-7.2771266435510700e+03, new int[] { 0, 0, 11 });
            p.AddCoeff(2.3324123857535500e+03, new int[] { 0, 0, 13 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a01e93be-d48e-4de7-9161-de81267c6d76"));
            OrthonormalPolynomials[456] = p;
            p.AddCoeff(6.9071305002797200e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(-5.3875617902181800e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(6.7344522377727300e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(-3.0529516811236400e+03, new int[] { 0, 1, 6 });
            p.AddCoeff(6.2149373508588200e+03, new int[] { 0, 1, 8 });
            p.AddCoeff(-5.8006081941349100e+03, new int[] { 0, 1, 10 });
            p.AddCoeff(2.0214240676530700e+03, new int[] { 0, 1, 12 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("29e58f80-e510-40aa-b284-e44ba2805181"));
            OrthonormalPolynomials[457] = p;
            p.AddCoeff(5.1317701979762900e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-1.1118835428948600e+02, new int[] { 0, 0, 3 });
            p.AddCoeff(6.6713012573691800e+02, new int[] { 0, 0, 5 });
            p.AddCoeff(-1.6201731625039400e+03, new int[] { 0, 0, 7 });
            p.AddCoeff(1.7101827826430500e+03, new int[] { 0, 0, 9 });
            p.AddCoeff(-6.5297888064552900e+02, new int[] { 0, 0, 11 });
            p.AddCoeff(-1.5395310593928900e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(3.3356506286845900e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(-2.0013903772107500e+03, new int[] { 0, 2, 5 });
            p.AddCoeff(4.8605194875118300e+03, new int[] { 0, 2, 7 });
            p.AddCoeff(-5.1305483479291600e+03, new int[] { 0, 2, 9 });
            p.AddCoeff(1.9589366419365900e+03, new int[] { 0, 2, 11 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5690406d-5067-41ee-bfcd-e2f380aae333"));
            OrthonormalPolynomials[458] = p;
            p.AddCoeff(1.5823608055186300e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-8.7029844303524500e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(7.5425865063054500e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(-2.2627759518916400e+03, new int[] { 0, 1, 6 });
            p.AddCoeff(2.7476565130112700e+03, new int[] { 0, 1, 8 });
            p.AddCoeff(-1.1601216388269800e+03, new int[] { 0, 1, 10 });
            p.AddCoeff(-2.6372680091977100e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(1.4504974050587400e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(-1.2570977510509100e+03, new int[] { 0, 3, 4 });
            p.AddCoeff(3.7712932531527300e+03, new int[] { 0, 3, 6 });
            p.AddCoeff(-4.5794275216854500e+03, new int[] { 0, 3, 8 });
            p.AddCoeff(1.9335360647116400e+03, new int[] { 0, 3, 10 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("19f489a1-099d-4475-a1df-f13b06133f6c"));
            OrthonormalPolynomials[459] = p;
            p.AddCoeff(4.2666293209026000e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-6.2577230039904800e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(2.4405119715562900e+02, new int[] { 0, 0, 5 });
            p.AddCoeff(-3.4864456736518400e+02, new int[] { 0, 0, 7 });
            p.AddCoeff(1.6463771236689200e+02, new int[] { 0, 0, 9 });
            p.AddCoeff(-4.2666293209026000e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(6.2577230039904800e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(-2.4405119715562900e+03, new int[] { 0, 2, 5 });
            p.AddCoeff(3.4864456736518400e+03, new int[] { 0, 2, 7 });
            p.AddCoeff(-1.6463771236689200e+03, new int[] { 0, 2, 9 });
            p.AddCoeff(4.9777342077197000e+01, new int[] { 0, 4, 1 });
            p.AddCoeff(-7.3006768379889000e+02, new int[] { 0, 4, 3 });
            p.AddCoeff(2.8472639668156700e+03, new int[] { 0, 4, 5 });
            p.AddCoeff(-4.0675199525938100e+03, new int[] { 0, 4, 7 });
            p.AddCoeff(1.9207733109470800e+03, new int[] { 0, 4, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("11d31350-15c1-4597-b893-f3f12d1872fd"));
            OrthonormalPolynomials[460] = p;
            p.AddCoeff(2.4787638654912600e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-8.9235499157685300e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(4.9079524536726900e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(-8.5071175863660000e+02, new int[] { 0, 1, 6 });
            p.AddCoeff(4.5573844212675000e+02, new int[] { 0, 1, 8 });
            p.AddCoeff(-1.1567564705625900e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(4.1643232940253100e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(-2.2903778117139200e+03, new int[] { 0, 3, 4 });
            p.AddCoeff(3.9699882069708000e+03, new int[] { 0, 3, 6 });
            p.AddCoeff(-2.1267793965915000e+03, new int[] { 0, 3, 8 });
            p.AddCoeff(1.0410808235063300e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(-3.7478909646227800e+02, new int[] { 0, 5, 2 });
            p.AddCoeff(2.0613400305425300e+03, new int[] { 0, 5, 4 });
            p.AddCoeff(-3.5729893862737200e+03, new int[] { 0, 5, 6 });
            p.AddCoeff(1.9141014569323500e+03, new int[] { 0, 5, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e92b4d6f-549c-43f2-ba9d-4881163012d9"));
            OrthonormalPolynomials[461] = p;
            p.AddCoeff(3.3749737208720800e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-3.0374763487848700e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(6.6824479673267200e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(-4.1367535035832100e+01, new int[] { 0, 0, 7 });
            p.AddCoeff(-7.0874448138313700e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(6.3787003324482400e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(-1.4033140731386100e+03, new int[] { 0, 2, 5 });
            p.AddCoeff(8.6871823575247400e+02, new int[] { 0, 2, 7 });
            p.AddCoeff(2.1262334441494100e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(-1.9136100997344700e+03, new int[] { 0, 4, 3 });
            p.AddCoeff(4.2099422194158400e+03, new int[] { 0, 4, 5 });
            p.AddCoeff(-2.6061547072574200e+03, new int[] { 0, 4, 7 });
            p.AddCoeff(-1.5592378590429000e+02, new int[] { 0, 6, 1 });
            p.AddCoeff(1.4033140731386100e+03, new int[] { 0, 6, 3 });
            p.AddCoeff(-3.0872909609049500e+03, new int[] { 0, 6, 5 });
            p.AddCoeff(1.9111801186554400e+03, new int[] { 0, 6, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a29b8c89-4316-408e-890c-c4b0fb153a90"));
            OrthonormalPolynomials[462] = p;
            p.AddCoeff(3.3749737208720800e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-7.0874448138313700e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(2.1262334441494100e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(-1.5592378590429000e+02, new int[] { 0, 1, 6 });
            p.AddCoeff(-3.0374763487848700e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(6.3787003324482400e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(-1.9136100997344700e+03, new int[] { 0, 3, 4 });
            p.AddCoeff(1.4033140731386100e+03, new int[] { 0, 3, 6 });
            p.AddCoeff(6.6824479673267200e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(-1.4033140731386100e+03, new int[] { 0, 5, 2 });
            p.AddCoeff(4.2099422194158400e+03, new int[] { 0, 5, 4 });
            p.AddCoeff(-3.0872909609049500e+03, new int[] { 0, 5, 6 });
            p.AddCoeff(-4.1367535035832100e+01, new int[] { 0, 7, 0 });
            p.AddCoeff(8.6871823575247400e+02, new int[] { 0, 7, 2 });
            p.AddCoeff(-2.6061547072574200e+03, new int[] { 0, 7, 4 });
            p.AddCoeff(1.9111801186554400e+03, new int[] { 0, 7, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("31544b2d-e670-4ff4-b492-8a74c3898cdb"));
            OrthonormalPolynomials[463] = p;
            p.AddCoeff(2.4787638654912600e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-1.1567564705625900e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(1.0410808235063300e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(-8.9235499157685300e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(4.1643232940253100e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(-3.7478909646227800e+02, new int[] { 0, 2, 5 });
            p.AddCoeff(4.9079524536726900e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(-2.2903778117139200e+03, new int[] { 0, 4, 3 });
            p.AddCoeff(2.0613400305425300e+03, new int[] { 0, 4, 5 });
            p.AddCoeff(-8.5071175863660000e+02, new int[] { 0, 6, 1 });
            p.AddCoeff(3.9699882069708000e+03, new int[] { 0, 6, 3 });
            p.AddCoeff(-3.5729893862737200e+03, new int[] { 0, 6, 5 });
            p.AddCoeff(4.5573844212675000e+02, new int[] { 0, 8, 1 });
            p.AddCoeff(-2.1267793965915000e+03, new int[] { 0, 8, 3 });
            p.AddCoeff(1.9141014569323500e+03, new int[] { 0, 8, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("ef428e8f-801c-4d3b-8adc-bb6b94508e3d"));
            OrthonormalPolynomials[464] = p;
            p.AddCoeff(4.2666293209026000e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-4.2666293209026000e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(4.9777342077197000e+01, new int[] { 0, 1, 4 });
            p.AddCoeff(-6.2577230039904800e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(6.2577230039904800e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(-7.3006768379889000e+02, new int[] { 0, 3, 4 });
            p.AddCoeff(2.4405119715562900e+02, new int[] { 0, 5, 0 });
            p.AddCoeff(-2.4405119715562900e+03, new int[] { 0, 5, 2 });
            p.AddCoeff(2.8472639668156700e+03, new int[] { 0, 5, 4 });
            p.AddCoeff(-3.4864456736518400e+02, new int[] { 0, 7, 0 });
            p.AddCoeff(3.4864456736518400e+03, new int[] { 0, 7, 2 });
            p.AddCoeff(-4.0675199525938100e+03, new int[] { 0, 7, 4 });
            p.AddCoeff(1.6463771236689200e+02, new int[] { 0, 9, 0 });
            p.AddCoeff(-1.6463771236689200e+03, new int[] { 0, 9, 2 });
            p.AddCoeff(1.9207733109470800e+03, new int[] { 0, 9, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("efcd185c-e1e6-4ebb-8605-bc093a27d708"));
            OrthonormalPolynomials[465] = p;
            p.AddCoeff(1.5823608055186300e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-2.6372680091977100e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(-8.7029844303524500e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(1.4504974050587400e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(7.5425865063054500e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(-1.2570977510509100e+03, new int[] { 0, 4, 3 });
            p.AddCoeff(-2.2627759518916400e+03, new int[] { 0, 6, 1 });
            p.AddCoeff(3.7712932531527300e+03, new int[] { 0, 6, 3 });
            p.AddCoeff(2.7476565130112700e+03, new int[] { 0, 8, 1 });
            p.AddCoeff(-4.5794275216854500e+03, new int[] { 0, 8, 3 });
            p.AddCoeff(-1.1601216388269800e+03, new int[] { 0, 10, 1 });
            p.AddCoeff(1.9335360647116400e+03, new int[] { 0, 10, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3be70ada-5fe1-4228-ba7c-07b00caba7e9"));
            OrthonormalPolynomials[466] = p;
            p.AddCoeff(5.1317701979762900e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-1.5395310593928900e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-1.1118835428948600e+02, new int[] { 0, 3, 0 });
            p.AddCoeff(3.3356506286845900e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(6.6713012573691800e+02, new int[] { 0, 5, 0 });
            p.AddCoeff(-2.0013903772107500e+03, new int[] { 0, 5, 2 });
            p.AddCoeff(-1.6201731625039400e+03, new int[] { 0, 7, 0 });
            p.AddCoeff(4.8605194875118300e+03, new int[] { 0, 7, 2 });
            p.AddCoeff(1.7101827826430500e+03, new int[] { 0, 9, 0 });
            p.AddCoeff(-5.1305483479291600e+03, new int[] { 0, 9, 2 });
            p.AddCoeff(-6.5297888064552900e+02, new int[] { 0, 11, 0 });
            p.AddCoeff(1.9589366419365900e+03, new int[] { 0, 11, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("8d83a219-bfaa-4e3a-9582-e92d31c6ac0c"));
            OrthonormalPolynomials[467] = p;
            p.AddCoeff(6.9071305002797200e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(-5.3875617902181800e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(6.7344522377727300e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(-3.0529516811236400e+03, new int[] { 0, 6, 1 });
            p.AddCoeff(6.2149373508588200e+03, new int[] { 0, 8, 1 });
            p.AddCoeff(-5.8006081941349100e+03, new int[] { 0, 10, 1 });
            p.AddCoeff(2.0214240676530700e+03, new int[] { 0, 12, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b544b0dd-36fd-4fdc-aab2-4b3b1265352e"));
            OrthonormalPolynomials[468] = p;
            p.AddCoeff(5.3875617902181800e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-1.6162685370654500e+02, new int[] { 0, 3, 0 });
            p.AddCoeff(1.3738282565056400e+03, new int[] { 0, 5, 0 });
            p.AddCoeff(-4.9719498806870600e+03, new int[] { 0, 7, 0 });
            p.AddCoeff(8.7009122912023700e+03, new int[] { 0, 9, 0 });
            p.AddCoeff(-7.2771266435510700e+03, new int[] { 0, 11, 0 });
            p.AddCoeff(2.3324123857535500e+03, new int[] { 0, 13, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("612dda05-e9d4-4ee1-8d60-5fe8aeacfe2e"));
            OrthonormalPolynomials[469] = p;
            p.AddCoeff(6.9071305002797200e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(-5.3875617902181800e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(6.7344522377727300e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(-3.0529516811236400e+03, new int[] { 1, 0, 6 });
            p.AddCoeff(6.2149373508588200e+03, new int[] { 1, 0, 8 });
            p.AddCoeff(-5.8006081941349100e+03, new int[] { 1, 0, 10 });
            p.AddCoeff(2.0214240676530700e+03, new int[] { 1, 0, 12 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("72371156-948b-42fe-87bc-c58e0a35e5b3"));
            OrthonormalPolynomials[470] = p;
            p.AddCoeff(-1.3769984409099100e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(2.9834966219714600e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(-1.7900979731828800e+03, new int[] { 1, 1, 5 });
            p.AddCoeff(4.3473807920155600e+03, new int[] { 1, 1, 7 });
            p.AddCoeff(-4.5889019471275400e+03, new int[] { 1, 1, 9 });
            p.AddCoeff(1.7521261979941500e+03, new int[] { 1, 1, 11 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("bd4327fb-a4d8-4c3c-9b92-ef2f7db71701"));
            OrthonormalPolynomials[471] = p;
            p.AddCoeff(7.7211301276826300e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(-4.2466215702254500e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(3.6804053608620600e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(-1.1041216082586200e+03, new int[] { 1, 0, 6 });
            p.AddCoeff(1.3407190957426100e+03, new int[] { 1, 0, 8 });
            p.AddCoeff(-5.6608139598021100e+02, new int[] { 1, 0, 10 });
            p.AddCoeff(-2.3163390383047900e+00, new int[] { 1, 2, 0 });
            p.AddCoeff(1.2739864710676300e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(-1.1041216082586200e+03, new int[] { 1, 2, 4 });
            p.AddCoeff(3.3123648247758500e+03, new int[] { 1, 2, 6 });
            p.AddCoeff(-4.0221572872278200e+03, new int[] { 1, 2, 8 });
            p.AddCoeff(1.6982441879406300e+03, new int[] { 1, 2, 10 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("95a88e9d-529a-41b9-9491-704a06887a35"));
            OrthonormalPolynomials[472] = p;
            p.AddCoeff(-2.6069535767139000e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(3.8235319125137100e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(-1.4911774458803500e+03, new int[] { 1, 1, 5 });
            p.AddCoeff(2.1302534941147800e+03, new int[] { 1, 1, 7 });
            p.AddCoeff(-1.0059530388875400e+03, new int[] { 1, 1, 9 });
            p.AddCoeff(4.3449226278564900e+01, new int[] { 1, 3, 1 });
            p.AddCoeff(-6.3725531875228600e+02, new int[] { 1, 3, 3 });
            p.AddCoeff(2.4852957431339100e+03, new int[] { 1, 3, 5 });
            p.AddCoeff(-3.5504224901913100e+03, new int[] { 1, 3, 7 });
            p.AddCoeff(1.6765883981458900e+03, new int[] { 1, 3, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d6188273-0ee4-4ccb-8267-f407339cc753"));
            OrthonormalPolynomials[473] = p;
            p.AddCoeff(7.7669532607032800e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(-2.7961031738531800e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(1.5378567456192500e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(-2.6656183590733600e+02, new int[] { 1, 0, 6 });
            p.AddCoeff(1.4280098352178700e+02, new int[] { 1, 0, 8 });
            p.AddCoeff(-7.7669532607032800e+00, new int[] { 1, 2, 0 });
            p.AddCoeff(2.7961031738531800e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(-1.5378567456192500e+03, new int[] { 1, 2, 4 });
            p.AddCoeff(2.6656183590733600e+03, new int[] { 1, 2, 6 });
            p.AddCoeff(-1.4280098352178700e+03, new int[] { 1, 2, 8 });
            p.AddCoeff(9.0614454708204900e+00, new int[] { 1, 4, 0 });
            p.AddCoeff(-3.2621203694953800e+02, new int[] { 1, 4, 2 });
            p.AddCoeff(1.7941662032224600e+03, new int[] { 1, 4, 4 });
            p.AddCoeff(-3.1098880855855900e+03, new int[] { 1, 4, 6 });
            p.AddCoeff(1.6660114744208500e+03, new int[] { 1, 4, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c5967740-e675-4819-87d5-f22da4331dad"));
            OrthonormalPolynomials[474] = p;
            p.AddCoeff(-3.2263162809921600e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(2.9036846528929400e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(-6.3881062363644800e+02, new int[] { 1, 1, 5 });
            p.AddCoeff(3.9545419558446800e+02, new int[] { 1, 1, 7 });
            p.AddCoeff(1.5056142644630100e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(-1.3550528380167100e+03, new int[] { 1, 3, 3 });
            p.AddCoeff(2.9811162436367600e+03, new int[] { 1, 3, 5 });
            p.AddCoeff(-1.8454529127275200e+03, new int[] { 1, 3, 7 });
            p.AddCoeff(-1.3550528380167100e+02, new int[] { 1, 5, 1 });
            p.AddCoeff(1.2195475542150400e+03, new int[] { 1, 5, 3 });
            p.AddCoeff(-2.6830046192730800e+03, new int[] { 1, 5, 5 });
            p.AddCoeff(1.6609076214547600e+03, new int[] { 1, 5, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7821eb09-4728-4284-9f24-65ccee888a32"));
            OrthonormalPolynomials[475] = p;
            p.AddCoeff(7.7742594375442700e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(-1.6325944818843000e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(4.8977834456528900e+01, new int[] { 1, 0, 4 });
            p.AddCoeff(-3.5917078601454500e+01, new int[] { 1, 0, 6 });
            p.AddCoeff(-1.6325944818843000e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(3.4284484119570200e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(-1.0285345235871100e+03, new int[] { 1, 2, 4 });
            p.AddCoeff(7.5425865063054500e+02, new int[] { 1, 2, 6 });
            p.AddCoeff(4.8977834456528900e+01, new int[] { 1, 4, 0 });
            p.AddCoeff(-1.0285345235871100e+03, new int[] { 1, 4, 2 });
            p.AddCoeff(3.0856035707613200e+03, new int[] { 1, 4, 4 });
            p.AddCoeff(-2.2627759518916400e+03, new int[] { 1, 4, 6 });
            p.AddCoeff(-3.5917078601454500e+01, new int[] { 1, 6, 0 });
            p.AddCoeff(7.5425865063054500e+02, new int[] { 1, 6, 2 });
            p.AddCoeff(-2.2627759518916400e+03, new int[] { 1, 6, 4 });
            p.AddCoeff(1.6593690313872000e+03, new int[] { 1, 6, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b48c5e70-d394-4552-9534-80c172c59c9b"));
            OrthonormalPolynomials[476] = p;
            p.AddCoeff(-3.2263162809921600e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(1.5056142644630100e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(-1.3550528380167100e+02, new int[] { 1, 1, 5 });
            p.AddCoeff(2.9036846528929400e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(-1.3550528380167100e+03, new int[] { 1, 3, 3 });
            p.AddCoeff(1.2195475542150400e+03, new int[] { 1, 3, 5 });
            p.AddCoeff(-6.3881062363644800e+02, new int[] { 1, 5, 1 });
            p.AddCoeff(2.9811162436367600e+03, new int[] { 1, 5, 3 });
            p.AddCoeff(-2.6830046192730800e+03, new int[] { 1, 5, 5 });
            p.AddCoeff(3.9545419558446800e+02, new int[] { 1, 7, 1 });
            p.AddCoeff(-1.8454529127275200e+03, new int[] { 1, 7, 3 });
            p.AddCoeff(1.6609076214547600e+03, new int[] { 1, 7, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e1fe03be-fac4-4787-8d72-76695c7b5111"));
            OrthonormalPolynomials[477] = p;
            p.AddCoeff(7.7669532607032800e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(-7.7669532607032800e+00, new int[] { 1, 0, 2 });
            p.AddCoeff(9.0614454708204900e+00, new int[] { 1, 0, 4 });
            p.AddCoeff(-2.7961031738531800e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(2.7961031738531800e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(-3.2621203694953800e+02, new int[] { 1, 2, 4 });
            p.AddCoeff(1.5378567456192500e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(-1.5378567456192500e+03, new int[] { 1, 4, 2 });
            p.AddCoeff(1.7941662032224600e+03, new int[] { 1, 4, 4 });
            p.AddCoeff(-2.6656183590733600e+02, new int[] { 1, 6, 0 });
            p.AddCoeff(2.6656183590733600e+03, new int[] { 1, 6, 2 });
            p.AddCoeff(-3.1098880855855900e+03, new int[] { 1, 6, 4 });
            p.AddCoeff(1.4280098352178700e+02, new int[] { 1, 8, 0 });
            p.AddCoeff(-1.4280098352178700e+03, new int[] { 1, 8, 2 });
            p.AddCoeff(1.6660114744208500e+03, new int[] { 1, 8, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("03bc8239-65ee-4948-896b-2c70341af4b5"));
            OrthonormalPolynomials[478] = p;
            p.AddCoeff(-2.6069535767139000e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(4.3449226278564900e+01, new int[] { 1, 1, 3 });
            p.AddCoeff(3.8235319125137100e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(-6.3725531875228600e+02, new int[] { 1, 3, 3 });
            p.AddCoeff(-1.4911774458803500e+03, new int[] { 1, 5, 1 });
            p.AddCoeff(2.4852957431339100e+03, new int[] { 1, 5, 3 });
            p.AddCoeff(2.1302534941147800e+03, new int[] { 1, 7, 1 });
            p.AddCoeff(-3.5504224901913100e+03, new int[] { 1, 7, 3 });
            p.AddCoeff(-1.0059530388875400e+03, new int[] { 1, 9, 1 });
            p.AddCoeff(1.6765883981458900e+03, new int[] { 1, 9, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("8c2e84df-9312-4c2d-8c09-47ca753548d1"));
            OrthonormalPolynomials[479] = p;
            p.AddCoeff(7.7211301276826300e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(-2.3163390383047900e+00, new int[] { 1, 0, 2 });
            p.AddCoeff(-4.2466215702254500e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(1.2739864710676300e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(3.6804053608620600e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(-1.1041216082586200e+03, new int[] { 1, 4, 2 });
            p.AddCoeff(-1.1041216082586200e+03, new int[] { 1, 6, 0 });
            p.AddCoeff(3.3123648247758500e+03, new int[] { 1, 6, 2 });
            p.AddCoeff(1.3407190957426100e+03, new int[] { 1, 8, 0 });
            p.AddCoeff(-4.0221572872278200e+03, new int[] { 1, 8, 2 });
            p.AddCoeff(-5.6608139598021100e+02, new int[] { 1, 10, 0 });
            p.AddCoeff(1.6982441879406300e+03, new int[] { 1, 10, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3967099f-b274-494b-bed6-8318ab399510"));
            OrthonormalPolynomials[480] = p;
            p.AddCoeff(-1.3769984409099100e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(2.9834966219714600e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(-1.7900979731828800e+03, new int[] { 1, 5, 1 });
            p.AddCoeff(4.3473807920155600e+03, new int[] { 1, 7, 1 });
            p.AddCoeff(-4.5889019471275400e+03, new int[] { 1, 9, 1 });
            p.AddCoeff(1.7521261979941500e+03, new int[] { 1, 11, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("0e27986f-81e3-4046-a0ac-f8b8b8bb72ff"));
            OrthonormalPolynomials[481] = p;
            p.AddCoeff(6.9071305002797200e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(-5.3875617902181800e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(6.7344522377727300e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(-3.0529516811236400e+03, new int[] { 1, 6, 0 });
            p.AddCoeff(6.2149373508588200e+03, new int[] { 1, 8, 0 });
            p.AddCoeff(-5.8006081941349100e+03, new int[] { 1, 10, 0 });
            p.AddCoeff(2.0214240676530700e+03, new int[] { 1, 12, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b59ff6fe-059f-4120-81f2-f0d473df7859"));
            OrthonormalPolynomials[482] = p;
            p.AddCoeff(5.1317701979762900e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-1.1118835428948600e+02, new int[] { 0, 0, 3 });
            p.AddCoeff(6.6713012573691800e+02, new int[] { 0, 0, 5 });
            p.AddCoeff(-1.6201731625039400e+03, new int[] { 0, 0, 7 });
            p.AddCoeff(1.7101827826430500e+03, new int[] { 0, 0, 9 });
            p.AddCoeff(-6.5297888064552900e+02, new int[] { 0, 0, 11 });
            p.AddCoeff(-1.5395310593928900e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(3.3356506286845900e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(-2.0013903772107500e+03, new int[] { 2, 0, 5 });
            p.AddCoeff(4.8605194875118300e+03, new int[] { 2, 0, 7 });
            p.AddCoeff(-5.1305483479291600e+03, new int[] { 2, 0, 9 });
            p.AddCoeff(1.9589366419365900e+03, new int[] { 2, 0, 11 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e07bbf1d-54c7-4984-b8ed-c9eab9fd7441"));
            OrthonormalPolynomials[483] = p;
            p.AddCoeff(7.7211301276826300e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(-4.2466215702254500e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(3.6804053608620600e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(-1.1041216082586200e+03, new int[] { 0, 1, 6 });
            p.AddCoeff(1.3407190957426100e+03, new int[] { 0, 1, 8 });
            p.AddCoeff(-5.6608139598021100e+02, new int[] { 0, 1, 10 });
            p.AddCoeff(-2.3163390383047900e+00, new int[] { 2, 1, 0 });
            p.AddCoeff(1.2739864710676300e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(-1.1041216082586200e+03, new int[] { 2, 1, 4 });
            p.AddCoeff(3.3123648247758500e+03, new int[] { 2, 1, 6 });
            p.AddCoeff(-4.0221572872278200e+03, new int[] { 2, 1, 8 });
            p.AddCoeff(1.6982441879406300e+03, new int[] { 2, 1, 10 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d8dc0546-4c0c-47d0-bf76-ee82a6d77128"));
            OrthonormalPolynomials[484] = p;
            p.AddCoeff(4.7406992454473400e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-6.9530255599894300e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(2.7116799683958800e+02, new int[] { 0, 0, 5 });
            p.AddCoeff(-3.8738285262798200e+02, new int[] { 0, 0, 7 });
            p.AddCoeff(1.8293079151876900e+02, new int[] { 0, 0, 9 });
            p.AddCoeff(-1.4222097736342000e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(2.0859076679968300e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(-8.1350399051876300e+02, new int[] { 0, 2, 5 });
            p.AddCoeff(1.1621485578839500e+03, new int[] { 0, 2, 7 });
            p.AddCoeff(-5.4879237455630800e+02, new int[] { 0, 2, 9 });
            p.AddCoeff(-1.4222097736342000e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(2.0859076679968300e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(-8.1350399051876300e+02, new int[] { 2, 0, 5 });
            p.AddCoeff(1.1621485578839500e+03, new int[] { 2, 0, 7 });
            p.AddCoeff(-5.4879237455630800e+02, new int[] { 2, 0, 9 });
            p.AddCoeff(4.2666293209026000e+01, new int[] { 2, 2, 1 });
            p.AddCoeff(-6.2577230039904800e+02, new int[] { 2, 2, 3 });
            p.AddCoeff(2.4405119715562900e+03, new int[] { 2, 2, 5 });
            p.AddCoeff(-3.4864456736518400e+03, new int[] { 2, 2, 7 });
            p.AddCoeff(1.6463771236689200e+03, new int[] { 2, 2, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5b550fbe-9328-4fd6-98fe-d0e0732a1a5f"));
            OrthonormalPolynomials[485] = p;
            p.AddCoeff(1.7686130591985700e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-6.3670070131148600e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(3.5018538572131700e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(-6.0698800191695000e+02, new int[] { 0, 1, 6 });
            p.AddCoeff(3.2517214388408000e+02, new int[] { 0, 1, 8 });
            p.AddCoeff(-2.9476884319976200e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(1.0611678355191400e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(-5.8364230953552800e+02, new int[] { 0, 3, 4 });
            p.AddCoeff(1.0116466698615800e+03, new int[] { 0, 3, 6 });
            p.AddCoeff(-5.4195357314013400e+02, new int[] { 0, 3, 8 });
            p.AddCoeff(-5.3058391775957100e+00, new int[] { 2, 1, 0 });
            p.AddCoeff(1.9101021039344600e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(-1.0505561571639500e+03, new int[] { 2, 1, 4 });
            p.AddCoeff(1.8209640057508500e+03, new int[] { 2, 1, 6 });
            p.AddCoeff(-9.7551643165224000e+02, new int[] { 2, 1, 8 });
            p.AddCoeff(8.8430652959928600e+00, new int[] { 2, 3, 0 });
            p.AddCoeff(-3.1835035065574300e+02, new int[] { 2, 3, 2 });
            p.AddCoeff(1.7509269286065900e+03, new int[] { 2, 3, 4 });
            p.AddCoeff(-3.0349400095847500e+03, new int[] { 2, 3, 6 });
            p.AddCoeff(1.6258607194204000e+03, new int[] { 2, 3, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("ea1d3d42-f722-49be-9c4c-bb0d04dc0eae"));
            OrthonormalPolynomials[486] = p;
            p.AddCoeff(3.7675257274253000e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-3.3907731546827700e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(7.4597009403021000e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(-4.6179101059013000e+01, new int[] { 0, 0, 7 });
            p.AddCoeff(-3.7675257274253000e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(3.3907731546827700e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(-7.4597009403021000e+02, new int[] { 0, 2, 5 });
            p.AddCoeff(4.6179101059013000e+02, new int[] { 0, 2, 7 });
            p.AddCoeff(4.3954466819961800e+01, new int[] { 0, 4, 1 });
            p.AddCoeff(-3.9559020137965700e+02, new int[] { 0, 4, 3 });
            p.AddCoeff(8.7029844303524500e+02, new int[] { 0, 4, 5 });
            p.AddCoeff(-5.3875617902181800e+02, new int[] { 0, 4, 7 });
            p.AddCoeff(-1.1302577182275900e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(1.0172319464048300e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(-2.2379102820906300e+02, new int[] { 2, 0, 5 });
            p.AddCoeff(1.3853730317703900e+02, new int[] { 2, 0, 7 });
            p.AddCoeff(1.1302577182275900e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(-1.0172319464048300e+03, new int[] { 2, 2, 3 });
            p.AddCoeff(2.2379102820906300e+03, new int[] { 2, 2, 5 });
            p.AddCoeff(-1.3853730317703900e+03, new int[] { 2, 2, 7 });
            p.AddCoeff(-1.3186340045988600e+02, new int[] { 2, 4, 1 });
            p.AddCoeff(1.1867706041389700e+03, new int[] { 2, 4, 3 });
            p.AddCoeff(-2.6108953291057300e+03, new int[] { 2, 4, 5 });
            p.AddCoeff(1.6162685370654500e+03, new int[] { 2, 4, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a455fec9-8753-4bc8-8ef9-c5886f0de79b"));
            OrthonormalPolynomials[487] = p;
            p.AddCoeff(2.7696782814241800e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-5.8163243909907700e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(1.7448973172972300e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(-1.2795913660179700e+02, new int[] { 0, 1, 6 });
            p.AddCoeff(-1.2925165313312800e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(2.7142847157956900e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(-8.1428541473870800e+02, new int[] { 0, 3, 4 });
            p.AddCoeff(5.9714263747505300e+02, new int[] { 0, 3, 6 });
            p.AddCoeff(1.1632648781981500e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(-2.4428562442161200e+02, new int[] { 0, 5, 2 });
            p.AddCoeff(7.3285687326483700e+02, new int[] { 0, 5, 4 });
            p.AddCoeff(-5.3742837372754700e+02, new int[] { 0, 5, 6 });
            p.AddCoeff(-8.3090348442725300e+00, new int[] { 2, 1, 0 });
            p.AddCoeff(1.7448973172972300e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(-5.2346919518916900e+02, new int[] { 2, 1, 4 });
            p.AddCoeff(3.8387740980539100e+02, new int[] { 2, 1, 6 });
            p.AddCoeff(3.8775495939938500e+01, new int[] { 2, 3, 0 });
            p.AddCoeff(-8.1428541473870800e+02, new int[] { 2, 3, 2 });
            p.AddCoeff(2.4428562442161200e+03, new int[] { 2, 3, 4 });
            p.AddCoeff(-1.7914279124251600e+03, new int[] { 2, 3, 6 });
            p.AddCoeff(-3.4897946345944600e+01, new int[] { 2, 5, 0 });
            p.AddCoeff(7.3285687326483700e+02, new int[] { 2, 5, 2 });
            p.AddCoeff(-2.1985706197945100e+03, new int[] { 2, 5, 4 });
            p.AddCoeff(1.6122851211826400e+03, new int[] { 2, 5, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("8c0b2de2-fdf8-41d1-9027-2b90d637ea40"));
            OrthonormalPolynomials[488] = p;
            p.AddCoeff(2.7696782814241800e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-1.2925165313312800e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(1.1632648781981500e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(-5.8163243909907700e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(2.7142847157956900e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(-2.4428562442161200e+02, new int[] { 0, 2, 5 });
            p.AddCoeff(1.7448973172972300e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(-8.1428541473870800e+02, new int[] { 0, 4, 3 });
            p.AddCoeff(7.3285687326483700e+02, new int[] { 0, 4, 5 });
            p.AddCoeff(-1.2795913660179700e+02, new int[] { 0, 6, 1 });
            p.AddCoeff(5.9714263747505300e+02, new int[] { 0, 6, 3 });
            p.AddCoeff(-5.3742837372754700e+02, new int[] { 0, 6, 5 });
            p.AddCoeff(-8.3090348442725300e+00, new int[] { 2, 0, 1 });
            p.AddCoeff(3.8775495939938500e+01, new int[] { 2, 0, 3 });
            p.AddCoeff(-3.4897946345944600e+01, new int[] { 2, 0, 5 });
            p.AddCoeff(1.7448973172972300e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(-8.1428541473870800e+02, new int[] { 2, 2, 3 });
            p.AddCoeff(7.3285687326483700e+02, new int[] { 2, 2, 5 });
            p.AddCoeff(-5.2346919518916900e+02, new int[] { 2, 4, 1 });
            p.AddCoeff(2.4428562442161200e+03, new int[] { 2, 4, 3 });
            p.AddCoeff(-2.1985706197945100e+03, new int[] { 2, 4, 5 });
            p.AddCoeff(3.8387740980539100e+02, new int[] { 2, 6, 1 });
            p.AddCoeff(-1.7914279124251600e+03, new int[] { 2, 6, 3 });
            p.AddCoeff(1.6122851211826400e+03, new int[] { 2, 6, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("cddc6c7a-72e4-4b37-bde2-0a9692a2b986"));
            OrthonormalPolynomials[489] = p;
            p.AddCoeff(3.7675257274253000e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-3.7675257274253000e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(4.3954466819961800e+01, new int[] { 0, 1, 4 });
            p.AddCoeff(-3.3907731546827700e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(3.3907731546827700e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(-3.9559020137965700e+02, new int[] { 0, 3, 4 });
            p.AddCoeff(7.4597009403021000e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(-7.4597009403021000e+02, new int[] { 0, 5, 2 });
            p.AddCoeff(8.7029844303524500e+02, new int[] { 0, 5, 4 });
            p.AddCoeff(-4.6179101059013000e+01, new int[] { 0, 7, 0 });
            p.AddCoeff(4.6179101059013000e+02, new int[] { 0, 7, 2 });
            p.AddCoeff(-5.3875617902181800e+02, new int[] { 0, 7, 4 });
            p.AddCoeff(-1.1302577182275900e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(1.1302577182275900e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(-1.3186340045988600e+02, new int[] { 2, 1, 4 });
            p.AddCoeff(1.0172319464048300e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(-1.0172319464048300e+03, new int[] { 2, 3, 2 });
            p.AddCoeff(1.1867706041389700e+03, new int[] { 2, 3, 4 });
            p.AddCoeff(-2.2379102820906300e+02, new int[] { 2, 5, 0 });
            p.AddCoeff(2.2379102820906300e+03, new int[] { 2, 5, 2 });
            p.AddCoeff(-2.6108953291057300e+03, new int[] { 2, 5, 4 });
            p.AddCoeff(1.3853730317703900e+02, new int[] { 2, 7, 0 });
            p.AddCoeff(-1.3853730317703900e+03, new int[] { 2, 7, 2 });
            p.AddCoeff(1.6162685370654500e+03, new int[] { 2, 7, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5c2f61b1-5314-4823-a49e-3f341eda5a57"));
            OrthonormalPolynomials[490] = p;
            p.AddCoeff(1.7686130591985700e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-2.9476884319976200e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(-6.3670070131148600e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(1.0611678355191400e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(3.5018538572131700e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(-5.8364230953552800e+02, new int[] { 0, 4, 3 });
            p.AddCoeff(-6.0698800191695000e+02, new int[] { 0, 6, 1 });
            p.AddCoeff(1.0116466698615800e+03, new int[] { 0, 6, 3 });
            p.AddCoeff(3.2517214388408000e+02, new int[] { 0, 8, 1 });
            p.AddCoeff(-5.4195357314013400e+02, new int[] { 0, 8, 3 });
            p.AddCoeff(-5.3058391775957100e+00, new int[] { 2, 0, 1 });
            p.AddCoeff(8.8430652959928600e+00, new int[] { 2, 0, 3 });
            p.AddCoeff(1.9101021039344600e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(-3.1835035065574300e+02, new int[] { 2, 2, 3 });
            p.AddCoeff(-1.0505561571639500e+03, new int[] { 2, 4, 1 });
            p.AddCoeff(1.7509269286065900e+03, new int[] { 2, 4, 3 });
            p.AddCoeff(1.8209640057508500e+03, new int[] { 2, 6, 1 });
            p.AddCoeff(-3.0349400095847500e+03, new int[] { 2, 6, 3 });
            p.AddCoeff(-9.7551643165224000e+02, new int[] { 2, 8, 1 });
            p.AddCoeff(1.6258607194204000e+03, new int[] { 2, 8, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("6d3050a5-2418-468a-992b-7f44fb729296"));
            OrthonormalPolynomials[491] = p;
            p.AddCoeff(4.7406992454473400e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-1.4222097736342000e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-6.9530255599894300e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(2.0859076679968300e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(2.7116799683958800e+02, new int[] { 0, 5, 0 });
            p.AddCoeff(-8.1350399051876300e+02, new int[] { 0, 5, 2 });
            p.AddCoeff(-3.8738285262798200e+02, new int[] { 0, 7, 0 });
            p.AddCoeff(1.1621485578839500e+03, new int[] { 0, 7, 2 });
            p.AddCoeff(1.8293079151876900e+02, new int[] { 0, 9, 0 });
            p.AddCoeff(-5.4879237455630800e+02, new int[] { 0, 9, 2 });
            p.AddCoeff(-1.4222097736342000e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(4.2666293209026000e+01, new int[] { 2, 1, 2 });
            p.AddCoeff(2.0859076679968300e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(-6.2577230039904800e+02, new int[] { 2, 3, 2 });
            p.AddCoeff(-8.1350399051876300e+02, new int[] { 2, 5, 0 });
            p.AddCoeff(2.4405119715562900e+03, new int[] { 2, 5, 2 });
            p.AddCoeff(1.1621485578839500e+03, new int[] { 2, 7, 0 });
            p.AddCoeff(-3.4864456736518400e+03, new int[] { 2, 7, 2 });
            p.AddCoeff(-5.4879237455630800e+02, new int[] { 2, 9, 0 });
            p.AddCoeff(1.6463771236689200e+03, new int[] { 2, 9, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b60833ba-2ecb-4968-82c7-98d55f4a9e06"));
            OrthonormalPolynomials[492] = p;
            p.AddCoeff(7.7211301276826300e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(-4.2466215702254500e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(3.6804053608620600e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(-1.1041216082586200e+03, new int[] { 0, 6, 1 });
            p.AddCoeff(1.3407190957426100e+03, new int[] { 0, 8, 1 });
            p.AddCoeff(-5.6608139598021100e+02, new int[] { 0, 10, 1 });
            p.AddCoeff(-2.3163390383047900e+00, new int[] { 2, 0, 1 });
            p.AddCoeff(1.2739864710676300e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(-1.1041216082586200e+03, new int[] { 2, 4, 1 });
            p.AddCoeff(3.3123648247758500e+03, new int[] { 2, 6, 1 });
            p.AddCoeff(-4.0221572872278200e+03, new int[] { 2, 8, 1 });
            p.AddCoeff(1.6982441879406300e+03, new int[] { 2, 10, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d1af2279-6b92-4d13-8039-dd5e1b099547"));
            OrthonormalPolynomials[493] = p;
            p.AddCoeff(5.1317701979762900e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-1.1118835428948600e+02, new int[] { 0, 3, 0 });
            p.AddCoeff(6.6713012573691800e+02, new int[] { 0, 5, 0 });
            p.AddCoeff(-1.6201731625039400e+03, new int[] { 0, 7, 0 });
            p.AddCoeff(1.7101827826430500e+03, new int[] { 0, 9, 0 });
            p.AddCoeff(-6.5297888064552900e+02, new int[] { 0, 11, 0 });
            p.AddCoeff(-1.5395310593928900e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(3.3356506286845900e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(-2.0013903772107500e+03, new int[] { 2, 5, 0 });
            p.AddCoeff(4.8605194875118300e+03, new int[] { 2, 7, 0 });
            p.AddCoeff(-5.1305483479291600e+03, new int[] { 2, 9, 0 });
            p.AddCoeff(1.9589366419365900e+03, new int[] { 2, 11, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5904694b-309c-4abb-ac5f-70d16da1112c"));
            OrthonormalPolynomials[494] = p;
            p.AddCoeff(1.5823608055186300e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-8.7029844303524500e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(7.5425865063054500e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(-2.2627759518916400e+03, new int[] { 1, 0, 6 });
            p.AddCoeff(2.7476565130112700e+03, new int[] { 1, 0, 8 });
            p.AddCoeff(-1.1601216388269800e+03, new int[] { 1, 0, 10 });
            p.AddCoeff(-2.6372680091977100e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(1.4504974050587400e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(-1.2570977510509100e+03, new int[] { 3, 0, 4 });
            p.AddCoeff(3.7712932531527300e+03, new int[] { 3, 0, 6 });
            p.AddCoeff(-4.5794275216854500e+03, new int[] { 3, 0, 8 });
            p.AddCoeff(1.9335360647116400e+03, new int[] { 3, 0, 10 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("576ebd82-01d7-4ebd-8513-29bcfe1607d0"));
            OrthonormalPolynomials[495] = p;
            p.AddCoeff(-2.6069535767139000e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(3.8235319125137100e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(-1.4911774458803500e+03, new int[] { 1, 1, 5 });
            p.AddCoeff(2.1302534941147800e+03, new int[] { 1, 1, 7 });
            p.AddCoeff(-1.0059530388875400e+03, new int[] { 1, 1, 9 });
            p.AddCoeff(4.3449226278564900e+01, new int[] { 3, 1, 1 });
            p.AddCoeff(-6.3725531875228600e+02, new int[] { 3, 1, 3 });
            p.AddCoeff(2.4852957431339100e+03, new int[] { 3, 1, 5 });
            p.AddCoeff(-3.5504224901913100e+03, new int[] { 3, 1, 7 });
            p.AddCoeff(1.6765883981458900e+03, new int[] { 3, 1, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("8b877bf6-c1a5-4532-bc4b-38627a854c3e"));
            OrthonormalPolynomials[496] = p;
            p.AddCoeff(1.7686130591985700e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-6.3670070131148600e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(3.5018538572131700e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(-6.0698800191695000e+02, new int[] { 1, 0, 6 });
            p.AddCoeff(3.2517214388408000e+02, new int[] { 1, 0, 8 });
            p.AddCoeff(-5.3058391775957100e+00, new int[] { 1, 2, 0 });
            p.AddCoeff(1.9101021039344600e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(-1.0505561571639500e+03, new int[] { 1, 2, 4 });
            p.AddCoeff(1.8209640057508500e+03, new int[] { 1, 2, 6 });
            p.AddCoeff(-9.7551643165224000e+02, new int[] { 1, 2, 8 });
            p.AddCoeff(-2.9476884319976200e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(1.0611678355191400e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(-5.8364230953552800e+02, new int[] { 3, 0, 4 });
            p.AddCoeff(1.0116466698615800e+03, new int[] { 3, 0, 6 });
            p.AddCoeff(-5.4195357314013400e+02, new int[] { 3, 0, 8 });
            p.AddCoeff(8.8430652959928600e+00, new int[] { 3, 2, 0 });
            p.AddCoeff(-3.1835035065574300e+02, new int[] { 3, 2, 2 });
            p.AddCoeff(1.7509269286065900e+03, new int[] { 3, 2, 4 });
            p.AddCoeff(-3.0349400095847500e+03, new int[] { 3, 2, 6 });
            p.AddCoeff(1.6258607194204000e+03, new int[] { 3, 2, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d5a78221-72bd-420f-b812-73309aa6bdf6"));
            OrthonormalPolynomials[497] = p;
            p.AddCoeff(-4.7176884347612900e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(4.2459195912851600e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(-9.3410231008273600e+02, new int[] { 1, 1, 5 });
            p.AddCoeff(5.7825381100359800e+02, new int[] { 1, 1, 7 });
            p.AddCoeff(7.8628140579354900e+01, new int[] { 1, 3, 1 });
            p.AddCoeff(-7.0765326521419400e+02, new int[] { 1, 3, 3 });
            p.AddCoeff(1.5568371834712300e+03, new int[] { 1, 3, 5 });
            p.AddCoeff(-9.6375635167266400e+02, new int[] { 1, 3, 7 });
            p.AddCoeff(7.8628140579354900e+01, new int[] { 3, 1, 1 });
            p.AddCoeff(-7.0765326521419400e+02, new int[] { 3, 1, 3 });
            p.AddCoeff(1.5568371834712300e+03, new int[] { 3, 1, 5 });
            p.AddCoeff(-9.6375635167266400e+02, new int[] { 3, 1, 7 });
            p.AddCoeff(-1.3104690096559100e+02, new int[] { 3, 3, 1 });
            p.AddCoeff(1.1794221086903200e+03, new int[] { 3, 3, 3 });
            p.AddCoeff(-2.5947286391187100e+03, new int[] { 3, 3, 5 });
            p.AddCoeff(1.6062605861211100e+03, new int[] { 3, 3, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("15e55495-fea7-4db8-bb07-2d94d20e4df5"));
            OrthonormalPolynomials[498] = p;
            p.AddCoeff(1.7785640342151600e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-3.7349844718518300e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(1.1204953415555500e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(-8.2169658380740200e+01, new int[] { 1, 0, 6 });
            p.AddCoeff(-1.7785640342151600e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(3.7349844718518300e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(-1.1204953415555500e+03, new int[] { 1, 2, 4 });
            p.AddCoeff(8.2169658380740200e+02, new int[] { 1, 2, 6 });
            p.AddCoeff(2.0749913732510200e+01, new int[] { 1, 4, 0 });
            p.AddCoeff(-4.3574818838271300e+02, new int[] { 1, 4, 2 });
            p.AddCoeff(1.3072445651481400e+03, new int[] { 1, 4, 4 });
            p.AddCoeff(-9.5864601444196900e+02, new int[] { 1, 4, 6 });
            p.AddCoeff(-2.9642733903585900e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(6.2249741197530500e+01, new int[] { 3, 0, 2 });
            p.AddCoeff(-1.8674922359259100e+02, new int[] { 3, 0, 4 });
            p.AddCoeff(1.3694943063456700e+02, new int[] { 3, 0, 6 });
            p.AddCoeff(2.9642733903585900e+01, new int[] { 3, 2, 0 });
            p.AddCoeff(-6.2249741197530500e+02, new int[] { 3, 2, 2 });
            p.AddCoeff(1.8674922359259100e+03, new int[] { 3, 2, 4 });
            p.AddCoeff(-1.3694943063456700e+03, new int[] { 3, 2, 6 });
            p.AddCoeff(-3.4583189554183600e+01, new int[] { 3, 4, 0 });
            p.AddCoeff(7.2624698063785500e+02, new int[] { 3, 4, 2 });
            p.AddCoeff(-2.1787409419135700e+03, new int[] { 3, 4, 4 });
            p.AddCoeff(1.5977433574032800e+03, new int[] { 3, 4, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b4c4fb53-ce2a-40e1-8971-7f9351c844b0"));
            OrthonormalPolynomials[499] = p;
            p.AddCoeff(-5.4261340032805700e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(2.5321958681976000e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(-2.2789762813778400e+02, new int[] { 1, 1, 5 });
            p.AddCoeff(2.5321958681976000e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(-1.1816914051588800e+03, new int[] { 1, 3, 3 });
            p.AddCoeff(1.0635222646429900e+03, new int[] { 1, 3, 5 });
            p.AddCoeff(-2.2789762813778400e+02, new int[] { 1, 5, 1 });
            p.AddCoeff(1.0635222646429900e+03, new int[] { 1, 5, 3 });
            p.AddCoeff(-9.5717003817869200e+02, new int[] { 1, 5, 5 });
            p.AddCoeff(9.0435566721342800e+01, new int[] { 3, 1, 1 });
            p.AddCoeff(-4.2203264469960000e+02, new int[] { 3, 1, 3 });
            p.AddCoeff(3.7982938022964000e+02, new int[] { 3, 1, 5 });
            p.AddCoeff(-4.2203264469960000e+02, new int[] { 3, 3, 1 });
            p.AddCoeff(1.9694856752648000e+03, new int[] { 3, 3, 3 });
            p.AddCoeff(-1.7725371077383200e+03, new int[] { 3, 3, 5 });
            p.AddCoeff(3.7982938022964000e+02, new int[] { 3, 5, 1 });
            p.AddCoeff(-1.7725371077383200e+03, new int[] { 3, 5, 3 });
            p.AddCoeff(1.5952833969644900e+03, new int[] { 3, 5, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("6996a6fd-bbc2-41c5-835a-6b26b69a9386"));
            OrthonormalPolynomials[500] = p;
            p.AddCoeff(1.7785640342151600e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-1.7785640342151600e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(2.0749913732510200e+01, new int[] { 1, 0, 4 });
            p.AddCoeff(-3.7349844718518300e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(3.7349844718518300e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(-4.3574818838271300e+02, new int[] { 1, 2, 4 });
            p.AddCoeff(1.1204953415555500e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(-1.1204953415555500e+03, new int[] { 1, 4, 2 });
            p.AddCoeff(1.3072445651481400e+03, new int[] { 1, 4, 4 });
            p.AddCoeff(-8.2169658380740200e+01, new int[] { 1, 6, 0 });
            p.AddCoeff(8.2169658380740200e+02, new int[] { 1, 6, 2 });
            p.AddCoeff(-9.5864601444196900e+02, new int[] { 1, 6, 4 });
            p.AddCoeff(-2.9642733903585900e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(2.9642733903585900e+01, new int[] { 3, 0, 2 });
            p.AddCoeff(-3.4583189554183600e+01, new int[] { 3, 0, 4 });
            p.AddCoeff(6.2249741197530500e+01, new int[] { 3, 2, 0 });
            p.AddCoeff(-6.2249741197530500e+02, new int[] { 3, 2, 2 });
            p.AddCoeff(7.2624698063785500e+02, new int[] { 3, 2, 4 });
            p.AddCoeff(-1.8674922359259100e+02, new int[] { 3, 4, 0 });
            p.AddCoeff(1.8674922359259100e+03, new int[] { 3, 4, 2 });
            p.AddCoeff(-2.1787409419135700e+03, new int[] { 3, 4, 4 });
            p.AddCoeff(1.3694943063456700e+02, new int[] { 3, 6, 0 });
            p.AddCoeff(-1.3694943063456700e+03, new int[] { 3, 6, 2 });
            p.AddCoeff(1.5977433574032800e+03, new int[] { 3, 6, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("aec9f384-9abe-4e3c-9e4c-4062b001139b"));
            OrthonormalPolynomials[501] = p;
            p.AddCoeff(-4.7176884347612900e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(7.8628140579354900e+01, new int[] { 1, 1, 3 });
            p.AddCoeff(4.2459195912851600e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(-7.0765326521419400e+02, new int[] { 1, 3, 3 });
            p.AddCoeff(-9.3410231008273600e+02, new int[] { 1, 5, 1 });
            p.AddCoeff(1.5568371834712300e+03, new int[] { 1, 5, 3 });
            p.AddCoeff(5.7825381100359800e+02, new int[] { 1, 7, 1 });
            p.AddCoeff(-9.6375635167266400e+02, new int[] { 1, 7, 3 });
            p.AddCoeff(7.8628140579354900e+01, new int[] { 3, 1, 1 });
            p.AddCoeff(-1.3104690096559100e+02, new int[] { 3, 1, 3 });
            p.AddCoeff(-7.0765326521419400e+02, new int[] { 3, 3, 1 });
            p.AddCoeff(1.1794221086903200e+03, new int[] { 3, 3, 3 });
            p.AddCoeff(1.5568371834712300e+03, new int[] { 3, 5, 1 });
            p.AddCoeff(-2.5947286391187100e+03, new int[] { 3, 5, 3 });
            p.AddCoeff(-9.6375635167266400e+02, new int[] { 3, 7, 1 });
            p.AddCoeff(1.6062605861211100e+03, new int[] { 3, 7, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("2b260b8d-af7c-425e-99bb-96ae05fa0bda"));
            OrthonormalPolynomials[502] = p;
            p.AddCoeff(1.7686130591985700e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-5.3058391775957100e+00, new int[] { 1, 0, 2 });
            p.AddCoeff(-6.3670070131148600e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(1.9101021039344600e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(3.5018538572131700e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(-1.0505561571639500e+03, new int[] { 1, 4, 2 });
            p.AddCoeff(-6.0698800191695000e+02, new int[] { 1, 6, 0 });
            p.AddCoeff(1.8209640057508500e+03, new int[] { 1, 6, 2 });
            p.AddCoeff(3.2517214388408000e+02, new int[] { 1, 8, 0 });
            p.AddCoeff(-9.7551643165224000e+02, new int[] { 1, 8, 2 });
            p.AddCoeff(-2.9476884319976200e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(8.8430652959928600e+00, new int[] { 3, 0, 2 });
            p.AddCoeff(1.0611678355191400e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(-3.1835035065574300e+02, new int[] { 3, 2, 2 });
            p.AddCoeff(-5.8364230953552800e+02, new int[] { 3, 4, 0 });
            p.AddCoeff(1.7509269286065900e+03, new int[] { 3, 4, 2 });
            p.AddCoeff(1.0116466698615800e+03, new int[] { 3, 6, 0 });
            p.AddCoeff(-3.0349400095847500e+03, new int[] { 3, 6, 2 });
            p.AddCoeff(-5.4195357314013400e+02, new int[] { 3, 8, 0 });
            p.AddCoeff(1.6258607194204000e+03, new int[] { 3, 8, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("953e8db4-b094-4128-95cb-ea003b4ec83e"));
            OrthonormalPolynomials[503] = p;
            p.AddCoeff(-2.6069535767139000e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(3.8235319125137100e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(-1.4911774458803500e+03, new int[] { 1, 5, 1 });
            p.AddCoeff(2.1302534941147800e+03, new int[] { 1, 7, 1 });
            p.AddCoeff(-1.0059530388875400e+03, new int[] { 1, 9, 1 });
            p.AddCoeff(4.3449226278564900e+01, new int[] { 3, 1, 1 });
            p.AddCoeff(-6.3725531875228600e+02, new int[] { 3, 3, 1 });
            p.AddCoeff(2.4852957431339100e+03, new int[] { 3, 5, 1 });
            p.AddCoeff(-3.5504224901913100e+03, new int[] { 3, 7, 1 });
            p.AddCoeff(1.6765883981458900e+03, new int[] { 3, 9, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3eeb16c0-352a-4d77-ba98-9ead5cc72fe5"));
            OrthonormalPolynomials[504] = p;
            p.AddCoeff(1.5823608055186300e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-8.7029844303524500e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(7.5425865063054500e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(-2.2627759518916400e+03, new int[] { 1, 6, 0 });
            p.AddCoeff(2.7476565130112700e+03, new int[] { 1, 8, 0 });
            p.AddCoeff(-1.1601216388269800e+03, new int[] { 1, 10, 0 });
            p.AddCoeff(-2.6372680091977100e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(1.4504974050587400e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(-1.2570977510509100e+03, new int[] { 3, 4, 0 });
            p.AddCoeff(3.7712932531527300e+03, new int[] { 3, 6, 0 });
            p.AddCoeff(-4.5794275216854500e+03, new int[] { 3, 8, 0 });
            p.AddCoeff(1.9335360647116400e+03, new int[] { 3, 10, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d21b8356-b540-4e53-b1aa-baed9debd77e"));
            OrthonormalPolynomials[505] = p;
            p.AddCoeff(4.2666293209026000e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-6.2577230039904800e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(2.4405119715562900e+02, new int[] { 0, 0, 5 });
            p.AddCoeff(-3.4864456736518400e+02, new int[] { 0, 0, 7 });
            p.AddCoeff(1.6463771236689200e+02, new int[] { 0, 0, 9 });
            p.AddCoeff(-4.2666293209026000e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(6.2577230039904800e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(-2.4405119715562900e+03, new int[] { 2, 0, 5 });
            p.AddCoeff(3.4864456736518400e+03, new int[] { 2, 0, 7 });
            p.AddCoeff(-1.6463771236689200e+03, new int[] { 2, 0, 9 });
            p.AddCoeff(4.9777342077197000e+01, new int[] { 4, 0, 1 });
            p.AddCoeff(-7.3006768379889000e+02, new int[] { 4, 0, 3 });
            p.AddCoeff(2.8472639668156700e+03, new int[] { 4, 0, 5 });
            p.AddCoeff(-4.0675199525938100e+03, new int[] { 4, 0, 7 });
            p.AddCoeff(1.9207733109470800e+03, new int[] { 4, 0, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("4b1d5285-94c9-4d22-94d4-6fddc0a6beab"));
            OrthonormalPolynomials[506] = p;
            p.AddCoeff(7.7669532607032800e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(-2.7961031738531800e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(1.5378567456192500e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(-2.6656183590733600e+02, new int[] { 0, 1, 6 });
            p.AddCoeff(1.4280098352178700e+02, new int[] { 0, 1, 8 });
            p.AddCoeff(-7.7669532607032800e+00, new int[] { 2, 1, 0 });
            p.AddCoeff(2.7961031738531800e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(-1.5378567456192500e+03, new int[] { 2, 1, 4 });
            p.AddCoeff(2.6656183590733600e+03, new int[] { 2, 1, 6 });
            p.AddCoeff(-1.4280098352178700e+03, new int[] { 2, 1, 8 });
            p.AddCoeff(9.0614454708204900e+00, new int[] { 4, 1, 0 });
            p.AddCoeff(-3.2621203694953800e+02, new int[] { 4, 1, 2 });
            p.AddCoeff(1.7941662032224600e+03, new int[] { 4, 1, 4 });
            p.AddCoeff(-3.1098880855855900e+03, new int[] { 4, 1, 6 });
            p.AddCoeff(1.6660114744208500e+03, new int[] { 4, 1, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("9d2c713b-28da-4ae7-a148-4e89b6e5362d"));
            OrthonormalPolynomials[507] = p;
            p.AddCoeff(3.7675257274253000e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-3.3907731546827700e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(7.4597009403021000e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(-4.6179101059013000e+01, new int[] { 0, 0, 7 });
            p.AddCoeff(-1.1302577182275900e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(1.0172319464048300e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(-2.2379102820906300e+02, new int[] { 0, 2, 5 });
            p.AddCoeff(1.3853730317703900e+02, new int[] { 0, 2, 7 });
            p.AddCoeff(-3.7675257274253000e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(3.3907731546827700e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(-7.4597009403021000e+02, new int[] { 2, 0, 5 });
            p.AddCoeff(4.6179101059013000e+02, new int[] { 2, 0, 7 });
            p.AddCoeff(1.1302577182275900e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(-1.0172319464048300e+03, new int[] { 2, 2, 3 });
            p.AddCoeff(2.2379102820906300e+03, new int[] { 2, 2, 5 });
            p.AddCoeff(-1.3853730317703900e+03, new int[] { 2, 2, 7 });
            p.AddCoeff(4.3954466819961800e+01, new int[] { 4, 0, 1 });
            p.AddCoeff(-3.9559020137965700e+02, new int[] { 4, 0, 3 });
            p.AddCoeff(8.7029844303524500e+02, new int[] { 4, 0, 5 });
            p.AddCoeff(-5.3875617902181800e+02, new int[] { 4, 0, 7 });
            p.AddCoeff(-1.3186340045988600e+02, new int[] { 4, 2, 1 });
            p.AddCoeff(1.1867706041389700e+03, new int[] { 4, 2, 3 });
            p.AddCoeff(-2.6108953291057300e+03, new int[] { 4, 2, 5 });
            p.AddCoeff(1.6162685370654500e+03, new int[] { 4, 2, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("4e962c31-9ea7-4701-85cb-8131646128bb"));
            OrthonormalPolynomials[508] = p;
            p.AddCoeff(1.7785640342151600e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-3.7349844718518300e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(1.1204953415555500e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(-8.2169658380740200e+01, new int[] { 0, 1, 6 });
            p.AddCoeff(-2.9642733903585900e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(6.2249741197530500e+01, new int[] { 0, 3, 2 });
            p.AddCoeff(-1.8674922359259100e+02, new int[] { 0, 3, 4 });
            p.AddCoeff(1.3694943063456700e+02, new int[] { 0, 3, 6 });
            p.AddCoeff(-1.7785640342151600e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(3.7349844718518300e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(-1.1204953415555500e+03, new int[] { 2, 1, 4 });
            p.AddCoeff(8.2169658380740200e+02, new int[] { 2, 1, 6 });
            p.AddCoeff(2.9642733903585900e+01, new int[] { 2, 3, 0 });
            p.AddCoeff(-6.2249741197530500e+02, new int[] { 2, 3, 2 });
            p.AddCoeff(1.8674922359259100e+03, new int[] { 2, 3, 4 });
            p.AddCoeff(-1.3694943063456700e+03, new int[] { 2, 3, 6 });
            p.AddCoeff(2.0749913732510200e+01, new int[] { 4, 1, 0 });
            p.AddCoeff(-4.3574818838271300e+02, new int[] { 4, 1, 2 });
            p.AddCoeff(1.3072445651481400e+03, new int[] { 4, 1, 4 });
            p.AddCoeff(-9.5864601444196900e+02, new int[] { 4, 1, 6 });
            p.AddCoeff(-3.4583189554183600e+01, new int[] { 4, 3, 0 });
            p.AddCoeff(7.2624698063785500e+02, new int[] { 4, 3, 2 });
            p.AddCoeff(-2.1787409419135700e+03, new int[] { 4, 3, 4 });
            p.AddCoeff(1.5977433574032800e+03, new int[] { 4, 3, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("eb088234-fe41-4d75-a56a-411d92f0028c"));
            OrthonormalPolynomials[509] = p;
            p.AddCoeff(2.7826441153249400e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-1.2985672538183000e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(1.1687105284364700e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(-2.7826441153249400e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(1.2985672538183000e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(-1.1687105284364700e+02, new int[] { 0, 2, 5 });
            p.AddCoeff(3.2464181345457600e+01, new int[] { 0, 4, 1 });
            p.AddCoeff(-1.5149951294546900e+02, new int[] { 0, 4, 3 });
            p.AddCoeff(1.3634956165092200e+02, new int[] { 0, 4, 5 });
            p.AddCoeff(-2.7826441153249400e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(1.2985672538183000e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(-1.1687105284364700e+02, new int[] { 2, 0, 5 });
            p.AddCoeff(2.7826441153249400e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(-1.2985672538183000e+03, new int[] { 2, 2, 3 });
            p.AddCoeff(1.1687105284364700e+03, new int[] { 2, 2, 5 });
            p.AddCoeff(-3.2464181345457600e+02, new int[] { 2, 4, 1 });
            p.AddCoeff(1.5149951294546900e+03, new int[] { 2, 4, 3 });
            p.AddCoeff(-1.3634956165092200e+03, new int[] { 2, 4, 5 });
            p.AddCoeff(3.2464181345457600e+01, new int[] { 4, 0, 1 });
            p.AddCoeff(-1.5149951294546900e+02, new int[] { 4, 0, 3 });
            p.AddCoeff(1.3634956165092200e+02, new int[] { 4, 0, 5 });
            p.AddCoeff(-3.2464181345457600e+02, new int[] { 4, 2, 1 });
            p.AddCoeff(1.5149951294546900e+03, new int[] { 4, 2, 3 });
            p.AddCoeff(-1.3634956165092200e+03, new int[] { 4, 2, 5 });
            p.AddCoeff(3.7874878236367200e+02, new int[] { 4, 4, 1 });
            p.AddCoeff(-1.7674943176971300e+03, new int[] { 4, 4, 3 });
            p.AddCoeff(1.5907448859274200e+03, new int[] { 4, 4, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d0533952-bc80-4cd3-b54b-831f70f2fc7e"));
            OrthonormalPolynomials[510] = p;
            p.AddCoeff(2.7826441153249400e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-2.7826441153249400e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(3.2464181345457600e+01, new int[] { 0, 1, 4 });
            p.AddCoeff(-1.2985672538183000e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(1.2985672538183000e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(-1.5149951294546900e+02, new int[] { 0, 3, 4 });
            p.AddCoeff(1.1687105284364700e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(-1.1687105284364700e+02, new int[] { 0, 5, 2 });
            p.AddCoeff(1.3634956165092200e+02, new int[] { 0, 5, 4 });
            p.AddCoeff(-2.7826441153249400e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(2.7826441153249400e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(-3.2464181345457600e+02, new int[] { 2, 1, 4 });
            p.AddCoeff(1.2985672538183000e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(-1.2985672538183000e+03, new int[] { 2, 3, 2 });
            p.AddCoeff(1.5149951294546900e+03, new int[] { 2, 3, 4 });
            p.AddCoeff(-1.1687105284364700e+02, new int[] { 2, 5, 0 });
            p.AddCoeff(1.1687105284364700e+03, new int[] { 2, 5, 2 });
            p.AddCoeff(-1.3634956165092200e+03, new int[] { 2, 5, 4 });
            p.AddCoeff(3.2464181345457600e+01, new int[] { 4, 1, 0 });
            p.AddCoeff(-3.2464181345457600e+02, new int[] { 4, 1, 2 });
            p.AddCoeff(3.7874878236367200e+02, new int[] { 4, 1, 4 });
            p.AddCoeff(-1.5149951294546900e+02, new int[] { 4, 3, 0 });
            p.AddCoeff(1.5149951294546900e+03, new int[] { 4, 3, 2 });
            p.AddCoeff(-1.7674943176971300e+03, new int[] { 4, 3, 4 });
            p.AddCoeff(1.3634956165092200e+02, new int[] { 4, 5, 0 });
            p.AddCoeff(-1.3634956165092200e+03, new int[] { 4, 5, 2 });
            p.AddCoeff(1.5907448859274200e+03, new int[] { 4, 5, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("93e17e6c-d373-43d0-a966-9fe84439091f"));
            OrthonormalPolynomials[511] = p;
            p.AddCoeff(1.7785640342151600e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-2.9642733903585900e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(-3.7349844718518300e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(6.2249741197530500e+01, new int[] { 0, 2, 3 });
            p.AddCoeff(1.1204953415555500e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(-1.8674922359259100e+02, new int[] { 0, 4, 3 });
            p.AddCoeff(-8.2169658380740200e+01, new int[] { 0, 6, 1 });
            p.AddCoeff(1.3694943063456700e+02, new int[] { 0, 6, 3 });
            p.AddCoeff(-1.7785640342151600e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(2.9642733903585900e+01, new int[] { 2, 0, 3 });
            p.AddCoeff(3.7349844718518300e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(-6.2249741197530500e+02, new int[] { 2, 2, 3 });
            p.AddCoeff(-1.1204953415555500e+03, new int[] { 2, 4, 1 });
            p.AddCoeff(1.8674922359259100e+03, new int[] { 2, 4, 3 });
            p.AddCoeff(8.2169658380740200e+02, new int[] { 2, 6, 1 });
            p.AddCoeff(-1.3694943063456700e+03, new int[] { 2, 6, 3 });
            p.AddCoeff(2.0749913732510200e+01, new int[] { 4, 0, 1 });
            p.AddCoeff(-3.4583189554183600e+01, new int[] { 4, 0, 3 });
            p.AddCoeff(-4.3574818838271300e+02, new int[] { 4, 2, 1 });
            p.AddCoeff(7.2624698063785500e+02, new int[] { 4, 2, 3 });
            p.AddCoeff(1.3072445651481400e+03, new int[] { 4, 4, 1 });
            p.AddCoeff(-2.1787409419135700e+03, new int[] { 4, 4, 3 });
            p.AddCoeff(-9.5864601444196900e+02, new int[] { 4, 6, 1 });
            p.AddCoeff(1.5977433574032800e+03, new int[] { 4, 6, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("35297403-f7a6-4d42-9bc3-fcd524b6bc0f"));
            OrthonormalPolynomials[512] = p;
            p.AddCoeff(3.7675257274253000e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-1.1302577182275900e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-3.3907731546827700e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(1.0172319464048300e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(7.4597009403021000e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(-2.2379102820906300e+02, new int[] { 0, 5, 2 });
            p.AddCoeff(-4.6179101059013000e+01, new int[] { 0, 7, 0 });
            p.AddCoeff(1.3853730317703900e+02, new int[] { 0, 7, 2 });
            p.AddCoeff(-3.7675257274253000e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(1.1302577182275900e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(3.3907731546827700e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(-1.0172319464048300e+03, new int[] { 2, 3, 2 });
            p.AddCoeff(-7.4597009403021000e+02, new int[] { 2, 5, 0 });
            p.AddCoeff(2.2379102820906300e+03, new int[] { 2, 5, 2 });
            p.AddCoeff(4.6179101059013000e+02, new int[] { 2, 7, 0 });
            p.AddCoeff(-1.3853730317703900e+03, new int[] { 2, 7, 2 });
            p.AddCoeff(4.3954466819961800e+01, new int[] { 4, 1, 0 });
            p.AddCoeff(-1.3186340045988600e+02, new int[] { 4, 1, 2 });
            p.AddCoeff(-3.9559020137965700e+02, new int[] { 4, 3, 0 });
            p.AddCoeff(1.1867706041389700e+03, new int[] { 4, 3, 2 });
            p.AddCoeff(8.7029844303524500e+02, new int[] { 4, 5, 0 });
            p.AddCoeff(-2.6108953291057300e+03, new int[] { 4, 5, 2 });
            p.AddCoeff(-5.3875617902181800e+02, new int[] { 4, 7, 0 });
            p.AddCoeff(1.6162685370654500e+03, new int[] { 4, 7, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("9fa4735d-1043-4891-8fa4-970e64c6d2a3"));
            OrthonormalPolynomials[513] = p;
            p.AddCoeff(7.7669532607032800e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(-2.7961031738531800e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(1.5378567456192500e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(-2.6656183590733600e+02, new int[] { 0, 6, 1 });
            p.AddCoeff(1.4280098352178700e+02, new int[] { 0, 8, 1 });
            p.AddCoeff(-7.7669532607032800e+00, new int[] { 2, 0, 1 });
            p.AddCoeff(2.7961031738531800e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(-1.5378567456192500e+03, new int[] { 2, 4, 1 });
            p.AddCoeff(2.6656183590733600e+03, new int[] { 2, 6, 1 });
            p.AddCoeff(-1.4280098352178700e+03, new int[] { 2, 8, 1 });
            p.AddCoeff(9.0614454708204900e+00, new int[] { 4, 0, 1 });
            p.AddCoeff(-3.2621203694953800e+02, new int[] { 4, 2, 1 });
            p.AddCoeff(1.7941662032224600e+03, new int[] { 4, 4, 1 });
            p.AddCoeff(-3.1098880855855900e+03, new int[] { 4, 6, 1 });
            p.AddCoeff(1.6660114744208500e+03, new int[] { 4, 8, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("0601c234-39da-42f6-a7f6-161b4a70a8e3"));
            OrthonormalPolynomials[514] = p;
            p.AddCoeff(4.2666293209026000e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-6.2577230039904800e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(2.4405119715562900e+02, new int[] { 0, 5, 0 });
            p.AddCoeff(-3.4864456736518400e+02, new int[] { 0, 7, 0 });
            p.AddCoeff(1.6463771236689200e+02, new int[] { 0, 9, 0 });
            p.AddCoeff(-4.2666293209026000e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(6.2577230039904800e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(-2.4405119715562900e+03, new int[] { 2, 5, 0 });
            p.AddCoeff(3.4864456736518400e+03, new int[] { 2, 7, 0 });
            p.AddCoeff(-1.6463771236689200e+03, new int[] { 2, 9, 0 });
            p.AddCoeff(4.9777342077197000e+01, new int[] { 4, 1, 0 });
            p.AddCoeff(-7.3006768379889000e+02, new int[] { 4, 3, 0 });
            p.AddCoeff(2.8472639668156700e+03, new int[] { 4, 5, 0 });
            p.AddCoeff(-4.0675199525938100e+03, new int[] { 4, 7, 0 });
            p.AddCoeff(1.9207733109470800e+03, new int[] { 4, 9, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("57932b38-4a33-44e9-bdd3-072475c9be78"));
            OrthonormalPolynomials[515] = p;
            p.AddCoeff(2.4787638654912600e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-8.9235499157685300e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(4.9079524536726900e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(-8.5071175863660000e+02, new int[] { 1, 0, 6 });
            p.AddCoeff(4.5573844212675000e+02, new int[] { 1, 0, 8 });
            p.AddCoeff(-1.1567564705625900e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(4.1643232940253100e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(-2.2903778117139200e+03, new int[] { 3, 0, 4 });
            p.AddCoeff(3.9699882069708000e+03, new int[] { 3, 0, 6 });
            p.AddCoeff(-2.1267793965915000e+03, new int[] { 3, 0, 8 });
            p.AddCoeff(1.0410808235063300e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(-3.7478909646227800e+02, new int[] { 5, 0, 2 });
            p.AddCoeff(2.0613400305425300e+03, new int[] { 5, 0, 4 });
            p.AddCoeff(-3.5729893862737200e+03, new int[] { 5, 0, 6 });
            p.AddCoeff(1.9141014569323500e+03, new int[] { 5, 0, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5cde000e-777c-4c34-a32d-37c7aae2903e"));
            OrthonormalPolynomials[516] = p;
            p.AddCoeff(-3.2263162809921600e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(2.9036846528929400e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(-6.3881062363644800e+02, new int[] { 1, 1, 5 });
            p.AddCoeff(3.9545419558446800e+02, new int[] { 1, 1, 7 });
            p.AddCoeff(1.5056142644630100e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(-1.3550528380167100e+03, new int[] { 3, 1, 3 });
            p.AddCoeff(2.9811162436367600e+03, new int[] { 3, 1, 5 });
            p.AddCoeff(-1.8454529127275200e+03, new int[] { 3, 1, 7 });
            p.AddCoeff(-1.3550528380167100e+02, new int[] { 5, 1, 1 });
            p.AddCoeff(1.2195475542150400e+03, new int[] { 5, 1, 3 });
            p.AddCoeff(-2.6830046192730800e+03, new int[] { 5, 1, 5 });
            p.AddCoeff(1.6609076214547600e+03, new int[] { 5, 1, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("0e5b1411-b7da-45a5-a755-0d39c52870fd"));
            OrthonormalPolynomials[517] = p;
            p.AddCoeff(2.7696782814241800e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-5.8163243909907700e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(1.7448973172972300e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(-1.2795913660179700e+02, new int[] { 1, 0, 6 });
            p.AddCoeff(-8.3090348442725300e+00, new int[] { 1, 2, 0 });
            p.AddCoeff(1.7448973172972300e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(-5.2346919518916900e+02, new int[] { 1, 2, 4 });
            p.AddCoeff(3.8387740980539100e+02, new int[] { 1, 2, 6 });
            p.AddCoeff(-1.2925165313312800e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(2.7142847157956900e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(-8.1428541473870800e+02, new int[] { 3, 0, 4 });
            p.AddCoeff(5.9714263747505300e+02, new int[] { 3, 0, 6 });
            p.AddCoeff(3.8775495939938500e+01, new int[] { 3, 2, 0 });
            p.AddCoeff(-8.1428541473870800e+02, new int[] { 3, 2, 2 });
            p.AddCoeff(2.4428562442161200e+03, new int[] { 3, 2, 4 });
            p.AddCoeff(-1.7914279124251600e+03, new int[] { 3, 2, 6 });
            p.AddCoeff(1.1632648781981500e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(-2.4428562442161200e+02, new int[] { 5, 0, 2 });
            p.AddCoeff(7.3285687326483700e+02, new int[] { 5, 0, 4 });
            p.AddCoeff(-5.3742837372754700e+02, new int[] { 5, 0, 6 });
            p.AddCoeff(-3.4897946345944600e+01, new int[] { 5, 2, 0 });
            p.AddCoeff(7.3285687326483700e+02, new int[] { 5, 2, 2 });
            p.AddCoeff(-2.1985706197945100e+03, new int[] { 5, 2, 4 });
            p.AddCoeff(1.6122851211826400e+03, new int[] { 5, 2, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a36b4e77-18bd-405e-855f-38b6f5d2d322"));
            OrthonormalPolynomials[518] = p;
            p.AddCoeff(-5.4261340032805700e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(2.5321958681976000e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(-2.2789762813778400e+02, new int[] { 1, 1, 5 });
            p.AddCoeff(9.0435566721342800e+01, new int[] { 1, 3, 1 });
            p.AddCoeff(-4.2203264469960000e+02, new int[] { 1, 3, 3 });
            p.AddCoeff(3.7982938022964000e+02, new int[] { 1, 3, 5 });
            p.AddCoeff(2.5321958681976000e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(-1.1816914051588800e+03, new int[] { 3, 1, 3 });
            p.AddCoeff(1.0635222646429900e+03, new int[] { 3, 1, 5 });
            p.AddCoeff(-4.2203264469960000e+02, new int[] { 3, 3, 1 });
            p.AddCoeff(1.9694856752648000e+03, new int[] { 3, 3, 3 });
            p.AddCoeff(-1.7725371077383200e+03, new int[] { 3, 3, 5 });
            p.AddCoeff(-2.2789762813778400e+02, new int[] { 5, 1, 1 });
            p.AddCoeff(1.0635222646429900e+03, new int[] { 5, 1, 3 });
            p.AddCoeff(-9.5717003817869200e+02, new int[] { 5, 1, 5 });
            p.AddCoeff(3.7982938022964000e+02, new int[] { 5, 3, 1 });
            p.AddCoeff(-1.7725371077383200e+03, new int[] { 5, 3, 3 });
            p.AddCoeff(1.5952833969644900e+03, new int[] { 5, 3, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a83ea9eb-f58f-4b76-85b9-5951d3fc6972"));
            OrthonormalPolynomials[519] = p;
            p.AddCoeff(2.7826441153249400e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-2.7826441153249400e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(3.2464181345457600e+01, new int[] { 1, 0, 4 });
            p.AddCoeff(-2.7826441153249400e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(2.7826441153249400e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(-3.2464181345457600e+02, new int[] { 1, 2, 4 });
            p.AddCoeff(3.2464181345457600e+01, new int[] { 1, 4, 0 });
            p.AddCoeff(-3.2464181345457600e+02, new int[] { 1, 4, 2 });
            p.AddCoeff(3.7874878236367200e+02, new int[] { 1, 4, 4 });
            p.AddCoeff(-1.2985672538183000e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(1.2985672538183000e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(-1.5149951294546900e+02, new int[] { 3, 0, 4 });
            p.AddCoeff(1.2985672538183000e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(-1.2985672538183000e+03, new int[] { 3, 2, 2 });
            p.AddCoeff(1.5149951294546900e+03, new int[] { 3, 2, 4 });
            p.AddCoeff(-1.5149951294546900e+02, new int[] { 3, 4, 0 });
            p.AddCoeff(1.5149951294546900e+03, new int[] { 3, 4, 2 });
            p.AddCoeff(-1.7674943176971300e+03, new int[] { 3, 4, 4 });
            p.AddCoeff(1.1687105284364700e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(-1.1687105284364700e+02, new int[] { 5, 0, 2 });
            p.AddCoeff(1.3634956165092200e+02, new int[] { 5, 0, 4 });
            p.AddCoeff(-1.1687105284364700e+02, new int[] { 5, 2, 0 });
            p.AddCoeff(1.1687105284364700e+03, new int[] { 5, 2, 2 });
            p.AddCoeff(-1.3634956165092200e+03, new int[] { 5, 2, 4 });
            p.AddCoeff(1.3634956165092200e+02, new int[] { 5, 4, 0 });
            p.AddCoeff(-1.3634956165092200e+03, new int[] { 5, 4, 2 });
            p.AddCoeff(1.5907448859274200e+03, new int[] { 5, 4, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("1a63dd20-4dae-4a0c-8c71-9ab0a021860d"));
            OrthonormalPolynomials[520] = p;
            p.AddCoeff(-5.4261340032805700e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(9.0435566721342800e+01, new int[] { 1, 1, 3 });
            p.AddCoeff(2.5321958681976000e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(-4.2203264469960000e+02, new int[] { 1, 3, 3 });
            p.AddCoeff(-2.2789762813778400e+02, new int[] { 1, 5, 1 });
            p.AddCoeff(3.7982938022964000e+02, new int[] { 1, 5, 3 });
            p.AddCoeff(2.5321958681976000e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(-4.2203264469960000e+02, new int[] { 3, 1, 3 });
            p.AddCoeff(-1.1816914051588800e+03, new int[] { 3, 3, 1 });
            p.AddCoeff(1.9694856752648000e+03, new int[] { 3, 3, 3 });
            p.AddCoeff(1.0635222646429900e+03, new int[] { 3, 5, 1 });
            p.AddCoeff(-1.7725371077383200e+03, new int[] { 3, 5, 3 });
            p.AddCoeff(-2.2789762813778400e+02, new int[] { 5, 1, 1 });
            p.AddCoeff(3.7982938022964000e+02, new int[] { 5, 1, 3 });
            p.AddCoeff(1.0635222646429900e+03, new int[] { 5, 3, 1 });
            p.AddCoeff(-1.7725371077383200e+03, new int[] { 5, 3, 3 });
            p.AddCoeff(-9.5717003817869200e+02, new int[] { 5, 5, 1 });
            p.AddCoeff(1.5952833969644900e+03, new int[] { 5, 5, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("6eb5fc0f-00c3-4e4f-9301-62a1ed1c8e7c"));
            OrthonormalPolynomials[521] = p;
            p.AddCoeff(2.7696782814241800e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-8.3090348442725300e+00, new int[] { 1, 0, 2 });
            p.AddCoeff(-5.8163243909907700e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(1.7448973172972300e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(1.7448973172972300e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(-5.2346919518916900e+02, new int[] { 1, 4, 2 });
            p.AddCoeff(-1.2795913660179700e+02, new int[] { 1, 6, 0 });
            p.AddCoeff(3.8387740980539100e+02, new int[] { 1, 6, 2 });
            p.AddCoeff(-1.2925165313312800e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(3.8775495939938500e+01, new int[] { 3, 0, 2 });
            p.AddCoeff(2.7142847157956900e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(-8.1428541473870800e+02, new int[] { 3, 2, 2 });
            p.AddCoeff(-8.1428541473870800e+02, new int[] { 3, 4, 0 });
            p.AddCoeff(2.4428562442161200e+03, new int[] { 3, 4, 2 });
            p.AddCoeff(5.9714263747505300e+02, new int[] { 3, 6, 0 });
            p.AddCoeff(-1.7914279124251600e+03, new int[] { 3, 6, 2 });
            p.AddCoeff(1.1632648781981500e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(-3.4897946345944600e+01, new int[] { 5, 0, 2 });
            p.AddCoeff(-2.4428562442161200e+02, new int[] { 5, 2, 0 });
            p.AddCoeff(7.3285687326483700e+02, new int[] { 5, 2, 2 });
            p.AddCoeff(7.3285687326483700e+02, new int[] { 5, 4, 0 });
            p.AddCoeff(-2.1985706197945100e+03, new int[] { 5, 4, 2 });
            p.AddCoeff(-5.3742837372754700e+02, new int[] { 5, 6, 0 });
            p.AddCoeff(1.6122851211826400e+03, new int[] { 5, 6, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("2f087a60-8cb8-4fc6-a320-44d92a9039dc"));
            OrthonormalPolynomials[522] = p;
            p.AddCoeff(-3.2263162809921600e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(2.9036846528929400e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(-6.3881062363644800e+02, new int[] { 1, 5, 1 });
            p.AddCoeff(3.9545419558446800e+02, new int[] { 1, 7, 1 });
            p.AddCoeff(1.5056142644630100e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(-1.3550528380167100e+03, new int[] { 3, 3, 1 });
            p.AddCoeff(2.9811162436367600e+03, new int[] { 3, 5, 1 });
            p.AddCoeff(-1.8454529127275200e+03, new int[] { 3, 7, 1 });
            p.AddCoeff(-1.3550528380167100e+02, new int[] { 5, 1, 1 });
            p.AddCoeff(1.2195475542150400e+03, new int[] { 5, 3, 1 });
            p.AddCoeff(-2.6830046192730800e+03, new int[] { 5, 5, 1 });
            p.AddCoeff(1.6609076214547600e+03, new int[] { 5, 7, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("ac9a9869-077d-4810-b3ec-4bb49e615690"));
            OrthonormalPolynomials[523] = p;
            p.AddCoeff(2.4787638654912600e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-8.9235499157685300e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(4.9079524536726900e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(-8.5071175863660000e+02, new int[] { 1, 6, 0 });
            p.AddCoeff(4.5573844212675000e+02, new int[] { 1, 8, 0 });
            p.AddCoeff(-1.1567564705625900e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(4.1643232940253100e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(-2.2903778117139200e+03, new int[] { 3, 4, 0 });
            p.AddCoeff(3.9699882069708000e+03, new int[] { 3, 6, 0 });
            p.AddCoeff(-2.1267793965915000e+03, new int[] { 3, 8, 0 });
            p.AddCoeff(1.0410808235063300e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(-3.7478909646227800e+02, new int[] { 5, 2, 0 });
            p.AddCoeff(2.0613400305425300e+03, new int[] { 5, 4, 0 });
            p.AddCoeff(-3.5729893862737200e+03, new int[] { 5, 6, 0 });
            p.AddCoeff(1.9141014569323500e+03, new int[] { 5, 8, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("ce2cde6e-a5b3-4e90-afdd-2e6669f15fd6"));
            OrthonormalPolynomials[524] = p;
            p.AddCoeff(3.3749737208720800e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-3.0374763487848700e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(6.6824479673267200e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(-4.1367535035832100e+01, new int[] { 0, 0, 7 });
            p.AddCoeff(-7.0874448138313700e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(6.3787003324482400e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(-1.4033140731386100e+03, new int[] { 2, 0, 5 });
            p.AddCoeff(8.6871823575247400e+02, new int[] { 2, 0, 7 });
            p.AddCoeff(2.1262334441494100e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(-1.9136100997344700e+03, new int[] { 4, 0, 3 });
            p.AddCoeff(4.2099422194158400e+03, new int[] { 4, 0, 5 });
            p.AddCoeff(-2.6061547072574200e+03, new int[] { 4, 0, 7 });
            p.AddCoeff(-1.5592378590429000e+02, new int[] { 6, 0, 1 });
            p.AddCoeff(1.4033140731386100e+03, new int[] { 6, 0, 3 });
            p.AddCoeff(-3.0872909609049500e+03, new int[] { 6, 0, 5 });
            p.AddCoeff(1.9111801186554400e+03, new int[] { 6, 0, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("07c94a88-51df-4f5f-9688-20ca44cbbb49"));
            OrthonormalPolynomials[525] = p;
            p.AddCoeff(7.7742594375442700e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(-1.6325944818843000e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(4.8977834456528900e+01, new int[] { 0, 1, 4 });
            p.AddCoeff(-3.5917078601454500e+01, new int[] { 0, 1, 6 });
            p.AddCoeff(-1.6325944818843000e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(3.4284484119570200e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(-1.0285345235871100e+03, new int[] { 2, 1, 4 });
            p.AddCoeff(7.5425865063054500e+02, new int[] { 2, 1, 6 });
            p.AddCoeff(4.8977834456528900e+01, new int[] { 4, 1, 0 });
            p.AddCoeff(-1.0285345235871100e+03, new int[] { 4, 1, 2 });
            p.AddCoeff(3.0856035707613200e+03, new int[] { 4, 1, 4 });
            p.AddCoeff(-2.2627759518916400e+03, new int[] { 4, 1, 6 });
            p.AddCoeff(-3.5917078601454500e+01, new int[] { 6, 1, 0 });
            p.AddCoeff(7.5425865063054500e+02, new int[] { 6, 1, 2 });
            p.AddCoeff(-2.2627759518916400e+03, new int[] { 6, 1, 4 });
            p.AddCoeff(1.6593690313872000e+03, new int[] { 6, 1, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("240683a3-4828-4574-a9d1-c2d3c8b11303"));
            OrthonormalPolynomials[526] = p;
            p.AddCoeff(2.7696782814241800e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-1.2925165313312800e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(1.1632648781981500e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(-8.3090348442725300e+00, new int[] { 0, 2, 1 });
            p.AddCoeff(3.8775495939938500e+01, new int[] { 0, 2, 3 });
            p.AddCoeff(-3.4897946345944600e+01, new int[] { 0, 2, 5 });
            p.AddCoeff(-5.8163243909907700e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(2.7142847157956900e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(-2.4428562442161200e+02, new int[] { 2, 0, 5 });
            p.AddCoeff(1.7448973172972300e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(-8.1428541473870800e+02, new int[] { 2, 2, 3 });
            p.AddCoeff(7.3285687326483700e+02, new int[] { 2, 2, 5 });
            p.AddCoeff(1.7448973172972300e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(-8.1428541473870800e+02, new int[] { 4, 0, 3 });
            p.AddCoeff(7.3285687326483700e+02, new int[] { 4, 0, 5 });
            p.AddCoeff(-5.2346919518916900e+02, new int[] { 4, 2, 1 });
            p.AddCoeff(2.4428562442161200e+03, new int[] { 4, 2, 3 });
            p.AddCoeff(-2.1985706197945100e+03, new int[] { 4, 2, 5 });
            p.AddCoeff(-1.2795913660179700e+02, new int[] { 6, 0, 1 });
            p.AddCoeff(5.9714263747505300e+02, new int[] { 6, 0, 3 });
            p.AddCoeff(-5.3742837372754700e+02, new int[] { 6, 0, 5 });
            p.AddCoeff(3.8387740980539100e+02, new int[] { 6, 2, 1 });
            p.AddCoeff(-1.7914279124251600e+03, new int[] { 6, 2, 3 });
            p.AddCoeff(1.6122851211826400e+03, new int[] { 6, 2, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5faa57e1-1d12-46ce-9980-3c4f79222819"));
            OrthonormalPolynomials[527] = p;
            p.AddCoeff(1.7785640342151600e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-1.7785640342151600e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(2.0749913732510200e+01, new int[] { 0, 1, 4 });
            p.AddCoeff(-2.9642733903585900e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(2.9642733903585900e+01, new int[] { 0, 3, 2 });
            p.AddCoeff(-3.4583189554183600e+01, new int[] { 0, 3, 4 });
            p.AddCoeff(-3.7349844718518300e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(3.7349844718518300e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(-4.3574818838271300e+02, new int[] { 2, 1, 4 });
            p.AddCoeff(6.2249741197530500e+01, new int[] { 2, 3, 0 });
            p.AddCoeff(-6.2249741197530500e+02, new int[] { 2, 3, 2 });
            p.AddCoeff(7.2624698063785500e+02, new int[] { 2, 3, 4 });
            p.AddCoeff(1.1204953415555500e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(-1.1204953415555500e+03, new int[] { 4, 1, 2 });
            p.AddCoeff(1.3072445651481400e+03, new int[] { 4, 1, 4 });
            p.AddCoeff(-1.8674922359259100e+02, new int[] { 4, 3, 0 });
            p.AddCoeff(1.8674922359259100e+03, new int[] { 4, 3, 2 });
            p.AddCoeff(-2.1787409419135700e+03, new int[] { 4, 3, 4 });
            p.AddCoeff(-8.2169658380740200e+01, new int[] { 6, 1, 0 });
            p.AddCoeff(8.2169658380740200e+02, new int[] { 6, 1, 2 });
            p.AddCoeff(-9.5864601444196900e+02, new int[] { 6, 1, 4 });
            p.AddCoeff(1.3694943063456700e+02, new int[] { 6, 3, 0 });
            p.AddCoeff(-1.3694943063456700e+03, new int[] { 6, 3, 2 });
            p.AddCoeff(1.5977433574032800e+03, new int[] { 6, 3, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("6c768ead-1e6d-4c64-9d80-b8b3fc7063d1"));
            OrthonormalPolynomials[528] = p;
            p.AddCoeff(1.7785640342151600e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-2.9642733903585900e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(-1.7785640342151600e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(2.9642733903585900e+01, new int[] { 0, 2, 3 });
            p.AddCoeff(2.0749913732510200e+01, new int[] { 0, 4, 1 });
            p.AddCoeff(-3.4583189554183600e+01, new int[] { 0, 4, 3 });
            p.AddCoeff(-3.7349844718518300e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(6.2249741197530500e+01, new int[] { 2, 0, 3 });
            p.AddCoeff(3.7349844718518300e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(-6.2249741197530500e+02, new int[] { 2, 2, 3 });
            p.AddCoeff(-4.3574818838271300e+02, new int[] { 2, 4, 1 });
            p.AddCoeff(7.2624698063785500e+02, new int[] { 2, 4, 3 });
            p.AddCoeff(1.1204953415555500e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(-1.8674922359259100e+02, new int[] { 4, 0, 3 });
            p.AddCoeff(-1.1204953415555500e+03, new int[] { 4, 2, 1 });
            p.AddCoeff(1.8674922359259100e+03, new int[] { 4, 2, 3 });
            p.AddCoeff(1.3072445651481400e+03, new int[] { 4, 4, 1 });
            p.AddCoeff(-2.1787409419135700e+03, new int[] { 4, 4, 3 });
            p.AddCoeff(-8.2169658380740200e+01, new int[] { 6, 0, 1 });
            p.AddCoeff(1.3694943063456700e+02, new int[] { 6, 0, 3 });
            p.AddCoeff(8.2169658380740200e+02, new int[] { 6, 2, 1 });
            p.AddCoeff(-1.3694943063456700e+03, new int[] { 6, 2, 3 });
            p.AddCoeff(-9.5864601444196900e+02, new int[] { 6, 4, 1 });
            p.AddCoeff(1.5977433574032800e+03, new int[] { 6, 4, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("398b9740-ef49-41d1-9564-8f45f6e7cf4d"));
            OrthonormalPolynomials[529] = p;
            p.AddCoeff(2.7696782814241800e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-8.3090348442725300e+00, new int[] { 0, 1, 2 });
            p.AddCoeff(-1.2925165313312800e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(3.8775495939938500e+01, new int[] { 0, 3, 2 });
            p.AddCoeff(1.1632648781981500e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(-3.4897946345944600e+01, new int[] { 0, 5, 2 });
            p.AddCoeff(-5.8163243909907700e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(1.7448973172972300e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(2.7142847157956900e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(-8.1428541473870800e+02, new int[] { 2, 3, 2 });
            p.AddCoeff(-2.4428562442161200e+02, new int[] { 2, 5, 0 });
            p.AddCoeff(7.3285687326483700e+02, new int[] { 2, 5, 2 });
            p.AddCoeff(1.7448973172972300e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(-5.2346919518916900e+02, new int[] { 4, 1, 2 });
            p.AddCoeff(-8.1428541473870800e+02, new int[] { 4, 3, 0 });
            p.AddCoeff(2.4428562442161200e+03, new int[] { 4, 3, 2 });
            p.AddCoeff(7.3285687326483700e+02, new int[] { 4, 5, 0 });
            p.AddCoeff(-2.1985706197945100e+03, new int[] { 4, 5, 2 });
            p.AddCoeff(-1.2795913660179700e+02, new int[] { 6, 1, 0 });
            p.AddCoeff(3.8387740980539100e+02, new int[] { 6, 1, 2 });
            p.AddCoeff(5.9714263747505300e+02, new int[] { 6, 3, 0 });
            p.AddCoeff(-1.7914279124251600e+03, new int[] { 6, 3, 2 });
            p.AddCoeff(-5.3742837372754700e+02, new int[] { 6, 5, 0 });
            p.AddCoeff(1.6122851211826400e+03, new int[] { 6, 5, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("671d0c7e-05bd-4027-9917-643c25f16ba9"));
            OrthonormalPolynomials[530] = p;
            p.AddCoeff(7.7742594375442700e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(-1.6325944818843000e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(4.8977834456528900e+01, new int[] { 0, 4, 1 });
            p.AddCoeff(-3.5917078601454500e+01, new int[] { 0, 6, 1 });
            p.AddCoeff(-1.6325944818843000e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(3.4284484119570200e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(-1.0285345235871100e+03, new int[] { 2, 4, 1 });
            p.AddCoeff(7.5425865063054500e+02, new int[] { 2, 6, 1 });
            p.AddCoeff(4.8977834456528900e+01, new int[] { 4, 0, 1 });
            p.AddCoeff(-1.0285345235871100e+03, new int[] { 4, 2, 1 });
            p.AddCoeff(3.0856035707613200e+03, new int[] { 4, 4, 1 });
            p.AddCoeff(-2.2627759518916400e+03, new int[] { 4, 6, 1 });
            p.AddCoeff(-3.5917078601454500e+01, new int[] { 6, 0, 1 });
            p.AddCoeff(7.5425865063054500e+02, new int[] { 6, 2, 1 });
            p.AddCoeff(-2.2627759518916400e+03, new int[] { 6, 4, 1 });
            p.AddCoeff(1.6593690313872000e+03, new int[] { 6, 6, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("084ab2c4-dfc4-4e41-9bf2-1e1d1bcdeccb"));
            OrthonormalPolynomials[531] = p;
            p.AddCoeff(3.3749737208720800e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-3.0374763487848700e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(6.6824479673267200e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(-4.1367535035832100e+01, new int[] { 0, 7, 0 });
            p.AddCoeff(-7.0874448138313700e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(6.3787003324482400e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(-1.4033140731386100e+03, new int[] { 2, 5, 0 });
            p.AddCoeff(8.6871823575247400e+02, new int[] { 2, 7, 0 });
            p.AddCoeff(2.1262334441494100e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(-1.9136100997344700e+03, new int[] { 4, 3, 0 });
            p.AddCoeff(4.2099422194158400e+03, new int[] { 4, 5, 0 });
            p.AddCoeff(-2.6061547072574200e+03, new int[] { 4, 7, 0 });
            p.AddCoeff(-1.5592378590429000e+02, new int[] { 6, 1, 0 });
            p.AddCoeff(1.4033140731386100e+03, new int[] { 6, 3, 0 });
            p.AddCoeff(-3.0872909609049500e+03, new int[] { 6, 5, 0 });
            p.AddCoeff(1.9111801186554400e+03, new int[] { 6, 7, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("09c53106-2dab-440e-b9f4-f1762c8a6eba"));
            OrthonormalPolynomials[532] = p;
            p.AddCoeff(3.3749737208720800e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-7.0874448138313700e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(2.1262334441494100e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(-1.5592378590429000e+02, new int[] { 1, 0, 6 });
            p.AddCoeff(-3.0374763487848700e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(6.3787003324482400e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(-1.9136100997344700e+03, new int[] { 3, 0, 4 });
            p.AddCoeff(1.4033140731386100e+03, new int[] { 3, 0, 6 });
            p.AddCoeff(6.6824479673267200e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(-1.4033140731386100e+03, new int[] { 5, 0, 2 });
            p.AddCoeff(4.2099422194158400e+03, new int[] { 5, 0, 4 });
            p.AddCoeff(-3.0872909609049500e+03, new int[] { 5, 0, 6 });
            p.AddCoeff(-4.1367535035832100e+01, new int[] { 7, 0, 0 });
            p.AddCoeff(8.6871823575247400e+02, new int[] { 7, 0, 2 });
            p.AddCoeff(-2.6061547072574200e+03, new int[] { 7, 0, 4 });
            p.AddCoeff(1.9111801186554400e+03, new int[] { 7, 0, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3adb1bbb-58ac-42cd-b673-5a746a2d1e14"));
            OrthonormalPolynomials[533] = p;
            p.AddCoeff(-3.2263162809921600e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(1.5056142644630100e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(-1.3550528380167100e+02, new int[] { 1, 1, 5 });
            p.AddCoeff(2.9036846528929400e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(-1.3550528380167100e+03, new int[] { 3, 1, 3 });
            p.AddCoeff(1.2195475542150400e+03, new int[] { 3, 1, 5 });
            p.AddCoeff(-6.3881062363644800e+02, new int[] { 5, 1, 1 });
            p.AddCoeff(2.9811162436367600e+03, new int[] { 5, 1, 3 });
            p.AddCoeff(-2.6830046192730800e+03, new int[] { 5, 1, 5 });
            p.AddCoeff(3.9545419558446800e+02, new int[] { 7, 1, 1 });
            p.AddCoeff(-1.8454529127275200e+03, new int[] { 7, 1, 3 });
            p.AddCoeff(1.6609076214547600e+03, new int[] { 7, 1, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5e99e1cc-47a7-44ce-85d6-118011404e7e"));
            OrthonormalPolynomials[534] = p;
            p.AddCoeff(3.7675257274253000e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-3.7675257274253000e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(4.3954466819961800e+01, new int[] { 1, 0, 4 });
            p.AddCoeff(-1.1302577182275900e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(1.1302577182275900e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(-1.3186340045988600e+02, new int[] { 1, 2, 4 });
            p.AddCoeff(-3.3907731546827700e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(3.3907731546827700e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(-3.9559020137965700e+02, new int[] { 3, 0, 4 });
            p.AddCoeff(1.0172319464048300e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(-1.0172319464048300e+03, new int[] { 3, 2, 2 });
            p.AddCoeff(1.1867706041389700e+03, new int[] { 3, 2, 4 });
            p.AddCoeff(7.4597009403021000e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(-7.4597009403021000e+02, new int[] { 5, 0, 2 });
            p.AddCoeff(8.7029844303524500e+02, new int[] { 5, 0, 4 });
            p.AddCoeff(-2.2379102820906300e+02, new int[] { 5, 2, 0 });
            p.AddCoeff(2.2379102820906300e+03, new int[] { 5, 2, 2 });
            p.AddCoeff(-2.6108953291057300e+03, new int[] { 5, 2, 4 });
            p.AddCoeff(-4.6179101059013000e+01, new int[] { 7, 0, 0 });
            p.AddCoeff(4.6179101059013000e+02, new int[] { 7, 0, 2 });
            p.AddCoeff(-5.3875617902181800e+02, new int[] { 7, 0, 4 });
            p.AddCoeff(1.3853730317703900e+02, new int[] { 7, 2, 0 });
            p.AddCoeff(-1.3853730317703900e+03, new int[] { 7, 2, 2 });
            p.AddCoeff(1.6162685370654500e+03, new int[] { 7, 2, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("9eb06d13-da7a-41fc-8637-c710bdb980af"));
            OrthonormalPolynomials[535] = p;
            p.AddCoeff(-4.7176884347612900e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(7.8628140579354900e+01, new int[] { 1, 1, 3 });
            p.AddCoeff(7.8628140579354900e+01, new int[] { 1, 3, 1 });
            p.AddCoeff(-1.3104690096559100e+02, new int[] { 1, 3, 3 });
            p.AddCoeff(4.2459195912851600e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(-7.0765326521419400e+02, new int[] { 3, 1, 3 });
            p.AddCoeff(-7.0765326521419400e+02, new int[] { 3, 3, 1 });
            p.AddCoeff(1.1794221086903200e+03, new int[] { 3, 3, 3 });
            p.AddCoeff(-9.3410231008273600e+02, new int[] { 5, 1, 1 });
            p.AddCoeff(1.5568371834712300e+03, new int[] { 5, 1, 3 });
            p.AddCoeff(1.5568371834712300e+03, new int[] { 5, 3, 1 });
            p.AddCoeff(-2.5947286391187100e+03, new int[] { 5, 3, 3 });
            p.AddCoeff(5.7825381100359800e+02, new int[] { 7, 1, 1 });
            p.AddCoeff(-9.6375635167266400e+02, new int[] { 7, 1, 3 });
            p.AddCoeff(-9.6375635167266400e+02, new int[] { 7, 3, 1 });
            p.AddCoeff(1.6062605861211100e+03, new int[] { 7, 3, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e94b8f2a-9bcf-4904-9635-780516a885f7"));
            OrthonormalPolynomials[536] = p;
            p.AddCoeff(3.7675257274253000e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-1.1302577182275900e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-3.7675257274253000e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(1.1302577182275900e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(4.3954466819961800e+01, new int[] { 1, 4, 0 });
            p.AddCoeff(-1.3186340045988600e+02, new int[] { 1, 4, 2 });
            p.AddCoeff(-3.3907731546827700e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(1.0172319464048300e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(3.3907731546827700e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(-1.0172319464048300e+03, new int[] { 3, 2, 2 });
            p.AddCoeff(-3.9559020137965700e+02, new int[] { 3, 4, 0 });
            p.AddCoeff(1.1867706041389700e+03, new int[] { 3, 4, 2 });
            p.AddCoeff(7.4597009403021000e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(-2.2379102820906300e+02, new int[] { 5, 0, 2 });
            p.AddCoeff(-7.4597009403021000e+02, new int[] { 5, 2, 0 });
            p.AddCoeff(2.2379102820906300e+03, new int[] { 5, 2, 2 });
            p.AddCoeff(8.7029844303524500e+02, new int[] { 5, 4, 0 });
            p.AddCoeff(-2.6108953291057300e+03, new int[] { 5, 4, 2 });
            p.AddCoeff(-4.6179101059013000e+01, new int[] { 7, 0, 0 });
            p.AddCoeff(1.3853730317703900e+02, new int[] { 7, 0, 2 });
            p.AddCoeff(4.6179101059013000e+02, new int[] { 7, 2, 0 });
            p.AddCoeff(-1.3853730317703900e+03, new int[] { 7, 2, 2 });
            p.AddCoeff(-5.3875617902181800e+02, new int[] { 7, 4, 0 });
            p.AddCoeff(1.6162685370654500e+03, new int[] { 7, 4, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a3d3e715-1978-4f1d-a143-a225bc9e7e52"));
            OrthonormalPolynomials[537] = p;
            p.AddCoeff(-3.2263162809921600e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(1.5056142644630100e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(-1.3550528380167100e+02, new int[] { 1, 5, 1 });
            p.AddCoeff(2.9036846528929400e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(-1.3550528380167100e+03, new int[] { 3, 3, 1 });
            p.AddCoeff(1.2195475542150400e+03, new int[] { 3, 5, 1 });
            p.AddCoeff(-6.3881062363644800e+02, new int[] { 5, 1, 1 });
            p.AddCoeff(2.9811162436367600e+03, new int[] { 5, 3, 1 });
            p.AddCoeff(-2.6830046192730800e+03, new int[] { 5, 5, 1 });
            p.AddCoeff(3.9545419558446800e+02, new int[] { 7, 1, 1 });
            p.AddCoeff(-1.8454529127275200e+03, new int[] { 7, 3, 1 });
            p.AddCoeff(1.6609076214547600e+03, new int[] { 7, 5, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("455dee67-72f5-433c-a8ce-ecd05004c5db"));
            OrthonormalPolynomials[538] = p;
            p.AddCoeff(3.3749737208720800e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-7.0874448138313700e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(2.1262334441494100e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(-1.5592378590429000e+02, new int[] { 1, 6, 0 });
            p.AddCoeff(-3.0374763487848700e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(6.3787003324482400e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(-1.9136100997344700e+03, new int[] { 3, 4, 0 });
            p.AddCoeff(1.4033140731386100e+03, new int[] { 3, 6, 0 });
            p.AddCoeff(6.6824479673267200e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(-1.4033140731386100e+03, new int[] { 5, 2, 0 });
            p.AddCoeff(4.2099422194158400e+03, new int[] { 5, 4, 0 });
            p.AddCoeff(-3.0872909609049500e+03, new int[] { 5, 6, 0 });
            p.AddCoeff(-4.1367535035832100e+01, new int[] { 7, 0, 0 });
            p.AddCoeff(8.6871823575247400e+02, new int[] { 7, 2, 0 });
            p.AddCoeff(-2.6061547072574200e+03, new int[] { 7, 4, 0 });
            p.AddCoeff(1.9111801186554400e+03, new int[] { 7, 6, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("93510d92-7e6e-4bc0-8dc1-1b35a217ece5"));
            OrthonormalPolynomials[539] = p;
            p.AddCoeff(2.4787638654912600e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-1.1567564705625900e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(1.0410808235063300e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(-8.9235499157685300e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(4.1643232940253100e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(-3.7478909646227800e+02, new int[] { 2, 0, 5 });
            p.AddCoeff(4.9079524536726900e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(-2.2903778117139200e+03, new int[] { 4, 0, 3 });
            p.AddCoeff(2.0613400305425300e+03, new int[] { 4, 0, 5 });
            p.AddCoeff(-8.5071175863660000e+02, new int[] { 6, 0, 1 });
            p.AddCoeff(3.9699882069708000e+03, new int[] { 6, 0, 3 });
            p.AddCoeff(-3.5729893862737200e+03, new int[] { 6, 0, 5 });
            p.AddCoeff(4.5573844212675000e+02, new int[] { 8, 0, 1 });
            p.AddCoeff(-2.1267793965915000e+03, new int[] { 8, 0, 3 });
            p.AddCoeff(1.9141014569323500e+03, new int[] { 8, 0, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("43b251fd-f4d7-41cc-bcdf-04e986ef0638"));
            OrthonormalPolynomials[540] = p;
            p.AddCoeff(7.7669532607032800e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(-7.7669532607032800e+00, new int[] { 0, 1, 2 });
            p.AddCoeff(9.0614454708204900e+00, new int[] { 0, 1, 4 });
            p.AddCoeff(-2.7961031738531800e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(2.7961031738531800e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(-3.2621203694953800e+02, new int[] { 2, 1, 4 });
            p.AddCoeff(1.5378567456192500e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(-1.5378567456192500e+03, new int[] { 4, 1, 2 });
            p.AddCoeff(1.7941662032224600e+03, new int[] { 4, 1, 4 });
            p.AddCoeff(-2.6656183590733600e+02, new int[] { 6, 1, 0 });
            p.AddCoeff(2.6656183590733600e+03, new int[] { 6, 1, 2 });
            p.AddCoeff(-3.1098880855855900e+03, new int[] { 6, 1, 4 });
            p.AddCoeff(1.4280098352178700e+02, new int[] { 8, 1, 0 });
            p.AddCoeff(-1.4280098352178700e+03, new int[] { 8, 1, 2 });
            p.AddCoeff(1.6660114744208500e+03, new int[] { 8, 1, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("8fa6d65c-62c1-4604-a307-0e369275322b"));
            OrthonormalPolynomials[541] = p;
            p.AddCoeff(1.7686130591985700e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-2.9476884319976200e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(-5.3058391775957100e+00, new int[] { 0, 2, 1 });
            p.AddCoeff(8.8430652959928600e+00, new int[] { 0, 2, 3 });
            p.AddCoeff(-6.3670070131148600e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(1.0611678355191400e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(1.9101021039344600e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(-3.1835035065574300e+02, new int[] { 2, 2, 3 });
            p.AddCoeff(3.5018538572131700e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(-5.8364230953552800e+02, new int[] { 4, 0, 3 });
            p.AddCoeff(-1.0505561571639500e+03, new int[] { 4, 2, 1 });
            p.AddCoeff(1.7509269286065900e+03, new int[] { 4, 2, 3 });
            p.AddCoeff(-6.0698800191695000e+02, new int[] { 6, 0, 1 });
            p.AddCoeff(1.0116466698615800e+03, new int[] { 6, 0, 3 });
            p.AddCoeff(1.8209640057508500e+03, new int[] { 6, 2, 1 });
            p.AddCoeff(-3.0349400095847500e+03, new int[] { 6, 2, 3 });
            p.AddCoeff(3.2517214388408000e+02, new int[] { 8, 0, 1 });
            p.AddCoeff(-5.4195357314013400e+02, new int[] { 8, 0, 3 });
            p.AddCoeff(-9.7551643165224000e+02, new int[] { 8, 2, 1 });
            p.AddCoeff(1.6258607194204000e+03, new int[] { 8, 2, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("9a105932-1958-48cd-8707-c1d88ce4014d"));
            OrthonormalPolynomials[542] = p;
            p.AddCoeff(1.7686130591985700e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-5.3058391775957100e+00, new int[] { 0, 1, 2 });
            p.AddCoeff(-2.9476884319976200e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(8.8430652959928600e+00, new int[] { 0, 3, 2 });
            p.AddCoeff(-6.3670070131148600e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(1.9101021039344600e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(1.0611678355191400e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(-3.1835035065574300e+02, new int[] { 2, 3, 2 });
            p.AddCoeff(3.5018538572131700e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(-1.0505561571639500e+03, new int[] { 4, 1, 2 });
            p.AddCoeff(-5.8364230953552800e+02, new int[] { 4, 3, 0 });
            p.AddCoeff(1.7509269286065900e+03, new int[] { 4, 3, 2 });
            p.AddCoeff(-6.0698800191695000e+02, new int[] { 6, 1, 0 });
            p.AddCoeff(1.8209640057508500e+03, new int[] { 6, 1, 2 });
            p.AddCoeff(1.0116466698615800e+03, new int[] { 6, 3, 0 });
            p.AddCoeff(-3.0349400095847500e+03, new int[] { 6, 3, 2 });
            p.AddCoeff(3.2517214388408000e+02, new int[] { 8, 1, 0 });
            p.AddCoeff(-9.7551643165224000e+02, new int[] { 8, 1, 2 });
            p.AddCoeff(-5.4195357314013400e+02, new int[] { 8, 3, 0 });
            p.AddCoeff(1.6258607194204000e+03, new int[] { 8, 3, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5c27e109-fb15-4dab-b1e3-2e765a57e0d0"));
            OrthonormalPolynomials[543] = p;
            p.AddCoeff(7.7669532607032800e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(-7.7669532607032800e+00, new int[] { 0, 2, 1 });
            p.AddCoeff(9.0614454708204900e+00, new int[] { 0, 4, 1 });
            p.AddCoeff(-2.7961031738531800e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(2.7961031738531800e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(-3.2621203694953800e+02, new int[] { 2, 4, 1 });
            p.AddCoeff(1.5378567456192500e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(-1.5378567456192500e+03, new int[] { 4, 2, 1 });
            p.AddCoeff(1.7941662032224600e+03, new int[] { 4, 4, 1 });
            p.AddCoeff(-2.6656183590733600e+02, new int[] { 6, 0, 1 });
            p.AddCoeff(2.6656183590733600e+03, new int[] { 6, 2, 1 });
            p.AddCoeff(-3.1098880855855900e+03, new int[] { 6, 4, 1 });
            p.AddCoeff(1.4280098352178700e+02, new int[] { 8, 0, 1 });
            p.AddCoeff(-1.4280098352178700e+03, new int[] { 8, 2, 1 });
            p.AddCoeff(1.6660114744208500e+03, new int[] { 8, 4, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d4d7d488-5c99-4b17-8f29-d40a10db87af"));
            OrthonormalPolynomials[544] = p;
            p.AddCoeff(2.4787638654912600e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-1.1567564705625900e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(1.0410808235063300e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(-8.9235499157685300e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(4.1643232940253100e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(-3.7478909646227800e+02, new int[] { 2, 5, 0 });
            p.AddCoeff(4.9079524536726900e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(-2.2903778117139200e+03, new int[] { 4, 3, 0 });
            p.AddCoeff(2.0613400305425300e+03, new int[] { 4, 5, 0 });
            p.AddCoeff(-8.5071175863660000e+02, new int[] { 6, 1, 0 });
            p.AddCoeff(3.9699882069708000e+03, new int[] { 6, 3, 0 });
            p.AddCoeff(-3.5729893862737200e+03, new int[] { 6, 5, 0 });
            p.AddCoeff(4.5573844212675000e+02, new int[] { 8, 1, 0 });
            p.AddCoeff(-2.1267793965915000e+03, new int[] { 8, 3, 0 });
            p.AddCoeff(1.9141014569323500e+03, new int[] { 8, 5, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("229eb65e-045c-4aaa-b8ad-bc40fe74a9e0"));
            OrthonormalPolynomials[545] = p;
            p.AddCoeff(4.2666293209026000e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-4.2666293209026000e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(4.9777342077197000e+01, new int[] { 1, 0, 4 });
            p.AddCoeff(-6.2577230039904800e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(6.2577230039904800e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(-7.3006768379889000e+02, new int[] { 3, 0, 4 });
            p.AddCoeff(2.4405119715562900e+02, new int[] { 5, 0, 0 });
            p.AddCoeff(-2.4405119715562900e+03, new int[] { 5, 0, 2 });
            p.AddCoeff(2.8472639668156700e+03, new int[] { 5, 0, 4 });
            p.AddCoeff(-3.4864456736518400e+02, new int[] { 7, 0, 0 });
            p.AddCoeff(3.4864456736518400e+03, new int[] { 7, 0, 2 });
            p.AddCoeff(-4.0675199525938100e+03, new int[] { 7, 0, 4 });
            p.AddCoeff(1.6463771236689200e+02, new int[] { 9, 0, 0 });
            p.AddCoeff(-1.6463771236689200e+03, new int[] { 9, 0, 2 });
            p.AddCoeff(1.9207733109470800e+03, new int[] { 9, 0, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a32aeaa6-6a48-4a2f-b311-704d118d5863"));
            OrthonormalPolynomials[546] = p;
            p.AddCoeff(-2.6069535767139000e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(4.3449226278564900e+01, new int[] { 1, 1, 3 });
            p.AddCoeff(3.8235319125137100e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(-6.3725531875228600e+02, new int[] { 3, 1, 3 });
            p.AddCoeff(-1.4911774458803500e+03, new int[] { 5, 1, 1 });
            p.AddCoeff(2.4852957431339100e+03, new int[] { 5, 1, 3 });
            p.AddCoeff(2.1302534941147800e+03, new int[] { 7, 1, 1 });
            p.AddCoeff(-3.5504224901913100e+03, new int[] { 7, 1, 3 });
            p.AddCoeff(-1.0059530388875400e+03, new int[] { 9, 1, 1 });
            p.AddCoeff(1.6765883981458900e+03, new int[] { 9, 1, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f79a1b23-7db4-4841-8dd4-558d5e94f1cb"));
            OrthonormalPolynomials[547] = p;
            p.AddCoeff(4.7406992454473400e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-1.4222097736342000e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-1.4222097736342000e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(4.2666293209026000e+01, new int[] { 1, 2, 2 });
            p.AddCoeff(-6.9530255599894300e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(2.0859076679968300e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(2.0859076679968300e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(-6.2577230039904800e+02, new int[] { 3, 2, 2 });
            p.AddCoeff(2.7116799683958800e+02, new int[] { 5, 0, 0 });
            p.AddCoeff(-8.1350399051876300e+02, new int[] { 5, 0, 2 });
            p.AddCoeff(-8.1350399051876300e+02, new int[] { 5, 2, 0 });
            p.AddCoeff(2.4405119715562900e+03, new int[] { 5, 2, 2 });
            p.AddCoeff(-3.8738285262798200e+02, new int[] { 7, 0, 0 });
            p.AddCoeff(1.1621485578839500e+03, new int[] { 7, 0, 2 });
            p.AddCoeff(1.1621485578839500e+03, new int[] { 7, 2, 0 });
            p.AddCoeff(-3.4864456736518400e+03, new int[] { 7, 2, 2 });
            p.AddCoeff(1.8293079151876900e+02, new int[] { 9, 0, 0 });
            p.AddCoeff(-5.4879237455630800e+02, new int[] { 9, 0, 2 });
            p.AddCoeff(-5.4879237455630800e+02, new int[] { 9, 2, 0 });
            p.AddCoeff(1.6463771236689200e+03, new int[] { 9, 2, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("24ab709f-ce0a-4f3d-aba1-361f9c899675"));
            OrthonormalPolynomials[548] = p;
            p.AddCoeff(-2.6069535767139000e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(4.3449226278564900e+01, new int[] { 1, 3, 1 });
            p.AddCoeff(3.8235319125137100e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(-6.3725531875228600e+02, new int[] { 3, 3, 1 });
            p.AddCoeff(-1.4911774458803500e+03, new int[] { 5, 1, 1 });
            p.AddCoeff(2.4852957431339100e+03, new int[] { 5, 3, 1 });
            p.AddCoeff(2.1302534941147800e+03, new int[] { 7, 1, 1 });
            p.AddCoeff(-3.5504224901913100e+03, new int[] { 7, 3, 1 });
            p.AddCoeff(-1.0059530388875400e+03, new int[] { 9, 1, 1 });
            p.AddCoeff(1.6765883981458900e+03, new int[] { 9, 3, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b737c35d-5ee6-476c-ade1-524894fbb50e"));
            OrthonormalPolynomials[549] = p;
            p.AddCoeff(4.2666293209026000e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-4.2666293209026000e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(4.9777342077197000e+01, new int[] { 1, 4, 0 });
            p.AddCoeff(-6.2577230039904800e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(6.2577230039904800e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(-7.3006768379889000e+02, new int[] { 3, 4, 0 });
            p.AddCoeff(2.4405119715562900e+02, new int[] { 5, 0, 0 });
            p.AddCoeff(-2.4405119715562900e+03, new int[] { 5, 2, 0 });
            p.AddCoeff(2.8472639668156700e+03, new int[] { 5, 4, 0 });
            p.AddCoeff(-3.4864456736518400e+02, new int[] { 7, 0, 0 });
            p.AddCoeff(3.4864456736518400e+03, new int[] { 7, 2, 0 });
            p.AddCoeff(-4.0675199525938100e+03, new int[] { 7, 4, 0 });
            p.AddCoeff(1.6463771236689200e+02, new int[] { 9, 0, 0 });
            p.AddCoeff(-1.6463771236689200e+03, new int[] { 9, 2, 0 });
            p.AddCoeff(1.9207733109470800e+03, new int[] { 9, 4, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a517671a-398a-42f0-be92-c0a856deb305"));
            OrthonormalPolynomials[550] = p;
            p.AddCoeff(1.5823608055186300e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(-2.6372680091977100e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(-8.7029844303524500e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(1.4504974050587400e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(7.5425865063054500e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(-1.2570977510509100e+03, new int[] { 4, 0, 3 });
            p.AddCoeff(-2.2627759518916400e+03, new int[] { 6, 0, 1 });
            p.AddCoeff(3.7712932531527300e+03, new int[] { 6, 0, 3 });
            p.AddCoeff(2.7476565130112700e+03, new int[] { 8, 0, 1 });
            p.AddCoeff(-4.5794275216854500e+03, new int[] { 8, 0, 3 });
            p.AddCoeff(-1.1601216388269800e+03, new int[] { 10, 0, 1 });
            p.AddCoeff(1.9335360647116400e+03, new int[] { 10, 0, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("8abbea5e-69b4-4545-8c31-c3c5410e12a1"));
            OrthonormalPolynomials[551] = p;
            p.AddCoeff(7.7211301276826300e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(-2.3163390383047900e+00, new int[] { 0, 1, 2 });
            p.AddCoeff(-4.2466215702254500e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(1.2739864710676300e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(3.6804053608620600e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(-1.1041216082586200e+03, new int[] { 4, 1, 2 });
            p.AddCoeff(-1.1041216082586200e+03, new int[] { 6, 1, 0 });
            p.AddCoeff(3.3123648247758500e+03, new int[] { 6, 1, 2 });
            p.AddCoeff(1.3407190957426100e+03, new int[] { 8, 1, 0 });
            p.AddCoeff(-4.0221572872278200e+03, new int[] { 8, 1, 2 });
            p.AddCoeff(-5.6608139598021100e+02, new int[] { 10, 1, 0 });
            p.AddCoeff(1.6982441879406300e+03, new int[] { 10, 1, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("82dd4086-0d61-4b2b-81ab-74dc3c905525"));
            OrthonormalPolynomials[552] = p;
            p.AddCoeff(7.7211301276826300e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(-2.3163390383047900e+00, new int[] { 0, 2, 1 });
            p.AddCoeff(-4.2466215702254500e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(1.2739864710676300e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(3.6804053608620600e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(-1.1041216082586200e+03, new int[] { 4, 2, 1 });
            p.AddCoeff(-1.1041216082586200e+03, new int[] { 6, 0, 1 });
            p.AddCoeff(3.3123648247758500e+03, new int[] { 6, 2, 1 });
            p.AddCoeff(1.3407190957426100e+03, new int[] { 8, 0, 1 });
            p.AddCoeff(-4.0221572872278200e+03, new int[] { 8, 2, 1 });
            p.AddCoeff(-5.6608139598021100e+02, new int[] { 10, 0, 1 });
            p.AddCoeff(1.6982441879406300e+03, new int[] { 10, 2, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("48d54c31-225d-47ec-b500-7964592dfb15"));
            OrthonormalPolynomials[553] = p;
            p.AddCoeff(1.5823608055186300e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(-2.6372680091977100e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(-8.7029844303524500e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(1.4504974050587400e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(7.5425865063054500e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(-1.2570977510509100e+03, new int[] { 4, 3, 0 });
            p.AddCoeff(-2.2627759518916400e+03, new int[] { 6, 1, 0 });
            p.AddCoeff(3.7712932531527300e+03, new int[] { 6, 3, 0 });
            p.AddCoeff(2.7476565130112700e+03, new int[] { 8, 1, 0 });
            p.AddCoeff(-4.5794275216854500e+03, new int[] { 8, 3, 0 });
            p.AddCoeff(-1.1601216388269800e+03, new int[] { 10, 1, 0 });
            p.AddCoeff(1.9335360647116400e+03, new int[] { 10, 3, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("49c1d740-05a1-4a0e-909b-34372f20fa1e"));
            OrthonormalPolynomials[554] = p;
            p.AddCoeff(5.1317701979762900e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-1.5395310593928900e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-1.1118835428948600e+02, new int[] { 3, 0, 0 });
            p.AddCoeff(3.3356506286845900e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(6.6713012573691800e+02, new int[] { 5, 0, 0 });
            p.AddCoeff(-2.0013903772107500e+03, new int[] { 5, 0, 2 });
            p.AddCoeff(-1.6201731625039400e+03, new int[] { 7, 0, 0 });
            p.AddCoeff(4.8605194875118300e+03, new int[] { 7, 0, 2 });
            p.AddCoeff(1.7101827826430500e+03, new int[] { 9, 0, 0 });
            p.AddCoeff(-5.1305483479291600e+03, new int[] { 9, 0, 2 });
            p.AddCoeff(-6.5297888064552900e+02, new int[] { 11, 0, 0 });
            p.AddCoeff(1.9589366419365900e+03, new int[] { 11, 0, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("ab0dc30c-8547-4c9a-9d8b-c8b1b36bdab3"));
            OrthonormalPolynomials[555] = p;
            p.AddCoeff(-1.3769984409099100e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(2.9834966219714600e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(-1.7900979731828800e+03, new int[] { 5, 1, 1 });
            p.AddCoeff(4.3473807920155600e+03, new int[] { 7, 1, 1 });
            p.AddCoeff(-4.5889019471275400e+03, new int[] { 9, 1, 1 });
            p.AddCoeff(1.7521261979941500e+03, new int[] { 11, 1, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("30d6e7e3-b234-4813-b5bf-66959aba4002"));
            OrthonormalPolynomials[556] = p;
            p.AddCoeff(5.1317701979762900e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-1.5395310593928900e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-1.1118835428948600e+02, new int[] { 3, 0, 0 });
            p.AddCoeff(3.3356506286845900e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(6.6713012573691800e+02, new int[] { 5, 0, 0 });
            p.AddCoeff(-2.0013903772107500e+03, new int[] { 5, 2, 0 });
            p.AddCoeff(-1.6201731625039400e+03, new int[] { 7, 0, 0 });
            p.AddCoeff(4.8605194875118300e+03, new int[] { 7, 2, 0 });
            p.AddCoeff(1.7101827826430500e+03, new int[] { 9, 0, 0 });
            p.AddCoeff(-5.1305483479291600e+03, new int[] { 9, 2, 0 });
            p.AddCoeff(-6.5297888064552900e+02, new int[] { 11, 0, 0 });
            p.AddCoeff(1.9589366419365900e+03, new int[] { 11, 2, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c434b708-4653-4b1a-8727-fd7f66b222c5"));
            OrthonormalPolynomials[557] = p;
            p.AddCoeff(6.9071305002797200e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(-5.3875617902181800e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(6.7344522377727300e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(-3.0529516811236400e+03, new int[] { 6, 0, 1 });
            p.AddCoeff(6.2149373508588200e+03, new int[] { 8, 0, 1 });
            p.AddCoeff(-5.8006081941349100e+03, new int[] { 10, 0, 1 });
            p.AddCoeff(2.0214240676530700e+03, new int[] { 12, 0, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("2ac3d810-94a0-4da1-81e2-e817e4d8bc1a"));
            OrthonormalPolynomials[558] = p;
            p.AddCoeff(6.9071305002797200e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(-5.3875617902181800e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(6.7344522377727300e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(-3.0529516811236400e+03, new int[] { 6, 1, 0 });
            p.AddCoeff(6.2149373508588200e+03, new int[] { 8, 1, 0 });
            p.AddCoeff(-5.8006081941349100e+03, new int[] { 10, 1, 0 });
            p.AddCoeff(2.0214240676530700e+03, new int[] { 12, 1, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d45aa9ea-270b-4d03-a386-2a273df54bc1"));
            OrthonormalPolynomials[559] = p;
            p.AddCoeff(5.3875617902181800e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(-1.6162685370654500e+02, new int[] { 3, 0, 0 });
            p.AddCoeff(1.3738282565056400e+03, new int[] { 5, 0, 0 });
            p.AddCoeff(-4.9719498806870600e+03, new int[] { 7, 0, 0 });
            p.AddCoeff(8.7009122912023700e+03, new int[] { 9, 0, 0 });
            p.AddCoeff(-7.2771266435510700e+03, new int[] { 11, 0, 0 });
            p.AddCoeff(2.3324123857535500e+03, new int[] { 13, 0, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("49d7c23f-e571-4e90-9780-df5437b961f1"));
            OrthonormalPolynomials[560] = p;
            p.AddCoeff(-3.9882405547065600e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(4.1876525824418900e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(-7.1190093901512200e+02, new int[] { 0, 0, 4 });
            p.AddCoeff(4.5087059470957700e+03, new int[] { 0, 0, 6 });
            p.AddCoeff(-1.3526117841287300e+04, new int[] { 0, 0, 8 });
            p.AddCoeff(2.0740047356640600e+04, new int[] { 0, 0, 10 });
            p.AddCoeff(-1.5712157088364100e+04, new int[] { 0, 0, 12 });
            p.AddCoeff(4.6618488064376900e+03, new int[] { 0, 0, 14 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("fa58f4fc-bba2-4e3e-bcc6-7a2f2da0b411"));
            OrthonormalPolynomials[561] = p;
            p.AddCoeff(9.3315307495746500e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-2.7994592248724000e+02, new int[] { 0, 1, 3 });
            p.AddCoeff(2.3795403411415300e+03, new int[] { 0, 1, 5 });
            p.AddCoeff(-8.6116698060360400e+03, new int[] { 0, 1, 7 });
            p.AddCoeff(1.5070422160563000e+04, new int[] { 0, 1, 9 });
            p.AddCoeff(-1.2604353079743700e+04, new int[] { 0, 1, 11 });
            p.AddCoeff(4.0398567563281000e+03, new int[] { 0, 1, 13 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5fa176d1-e18c-4211-8727-95ba81b03124"));
            OrthonormalPolynomials[562] = p;
            p.AddCoeff(-4.4585335662774400e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(3.4776561816964000e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(-4.3470702271204900e+02, new int[] { 0, 0, 4 });
            p.AddCoeff(1.9706718362946300e+03, new int[] { 0, 0, 6 });
            p.AddCoeff(-4.0117248095997700e+03, new int[] { 0, 0, 8 });
            p.AddCoeff(3.7442764889597800e+03, new int[] { 0, 0, 10 });
            p.AddCoeff(-1.3048236249405400e+03, new int[] { 0, 0, 12 });
            p.AddCoeff(1.3375600698832300e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-1.0432968545089200e+02, new int[] { 0, 2, 2 });
            p.AddCoeff(1.3041210681361500e+03, new int[] { 0, 2, 4 });
            p.AddCoeff(-5.9120155088838900e+03, new int[] { 0, 2, 6 });
            p.AddCoeff(1.2035174428799300e+04, new int[] { 0, 2, 8 });
            p.AddCoeff(-1.1232829466879400e+04, new int[] { 0, 2, 10 });
            p.AddCoeff(3.9144708748216000e+03, new int[] { 0, 2, 12 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("29cdc3c8-8c01-46d0-a1d3-84eed0612cc4"));
            OrthonormalPolynomials[563] = p;
            p.AddCoeff(1.8215977151856400e+01, new int[] { 0, 1, 1 });
            p.AddCoeff(-3.9467950495688900e+02, new int[] { 0, 1, 3 });
            p.AddCoeff(2.3680770297413400e+03, new int[] { 0, 1, 5 });
            p.AddCoeff(-5.7510442150861000e+03, new int[] { 0, 1, 7 });
            p.AddCoeff(6.0705466714797800e+03, new int[] { 0, 1, 9 });
            p.AddCoeff(-2.3178450927468200e+03, new int[] { 0, 1, 11 });
            p.AddCoeff(-3.0359961919760700e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(6.5779917492814900e+02, new int[] { 0, 3, 3 });
            p.AddCoeff(-3.9467950495688900e+03, new int[] { 0, 3, 5 });
            p.AddCoeff(9.5850736918101700e+03, new int[] { 0, 3, 7 });
            p.AddCoeff(-1.0117577785799600e+04, new int[] { 0, 3, 9 });
            p.AddCoeff(3.8630751545780400e+03, new int[] { 0, 3, 11 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b2f3596e-1cb5-4d5a-a40b-a3c2d5ab151b"));
            OrthonormalPolynomials[564] = p;
            p.AddCoeff(-4.4855712597622800e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(2.4670641928692500e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(-2.1381223004866800e+02, new int[] { 0, 0, 4 });
            p.AddCoeff(6.4143669014600500e+02, new int[] { 0, 0, 6 });
            p.AddCoeff(-7.7888740946301000e+02, new int[] { 0, 0, 8 });
            p.AddCoeff(3.2886357288438100e+02, new int[] { 0, 0, 10 });
            p.AddCoeff(4.4855712597622800e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-2.4670641928692500e+02, new int[] { 0, 2, 2 });
            p.AddCoeff(2.1381223004866800e+03, new int[] { 0, 2, 4 });
            p.AddCoeff(-6.4143669014600500e+03, new int[] { 0, 2, 6 });
            p.AddCoeff(7.7888740946301000e+03, new int[] { 0, 2, 8 });
            p.AddCoeff(-3.2886357288438100e+03, new int[] { 0, 2, 10 });
            p.AddCoeff(-5.2331664697226600e+00, new int[] { 0, 4, 0 });
            p.AddCoeff(2.8782415583474600e+02, new int[] { 0, 4, 2 });
            p.AddCoeff(-2.4944760172344700e+03, new int[] { 0, 4, 4 });
            p.AddCoeff(7.4834280517033600e+03, new int[] { 0, 4, 6 });
            p.AddCoeff(-9.0870197770684300e+03, new int[] { 0, 4, 8 });
            p.AddCoeff(3.8367416836511000e+03, new int[] { 0, 4, 10 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("4145cea7-8799-4f75-8037-15510c592b69"));
            OrthonormalPolynomials[565] = p;
            p.AddCoeff(2.3584680961604700e+01, new int[] { 0, 1, 1 });
            p.AddCoeff(-3.4590865410353500e+02, new int[] { 0, 1, 3 });
            p.AddCoeff(1.3490437510037900e+03, new int[] { 0, 1, 5 });
            p.AddCoeff(-1.9272053585768400e+03, new int[] { 0, 1, 7 });
            p.AddCoeff(9.1006919710573000e+02, new int[] { 0, 1, 9 });
            p.AddCoeff(-1.1006184448748800e+02, new int[] { 0, 3, 1 });
            p.AddCoeff(1.6142403858165000e+03, new int[] { 0, 3, 3 });
            p.AddCoeff(-6.2955375046843400e+03, new int[] { 0, 3, 5 });
            p.AddCoeff(8.9936250066919200e+03, new int[] { 0, 3, 7 });
            p.AddCoeff(-4.2469895864934000e+03, new int[] { 0, 3, 9 });
            p.AddCoeff(9.9055660038739600e+01, new int[] { 0, 5, 1 });
            p.AddCoeff(-1.4528163472348500e+03, new int[] { 0, 5, 3 });
            p.AddCoeff(5.6659837542159100e+03, new int[] { 0, 5, 5 });
            p.AddCoeff(-8.0942625060227200e+03, new int[] { 0, 5, 7 });
            p.AddCoeff(3.8222906278440600e+03, new int[] { 0, 5, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("4d98d42e-583b-4bdc-8c5d-73a869efff04"));
            OrthonormalPolynomials[566] = p;
            p.AddCoeff(-4.4911673672912700e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.6168202522248600e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(-8.8925113872367100e+01, new int[] { 0, 0, 4 });
            p.AddCoeff(1.5413686404543600e+02, new int[] { 0, 0, 6 });
            p.AddCoeff(-8.2573320024340800e+01, new int[] { 0, 0, 8 });
            p.AddCoeff(9.4314514713116600e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-3.3953225296722000e+02, new int[] { 0, 2, 2 });
            p.AddCoeff(1.8674273913197100e+03, new int[] { 0, 2, 4 });
            p.AddCoeff(-3.2368741449541600e+03, new int[] { 0, 2, 6 });
            p.AddCoeff(1.7340397205111600e+03, new int[] { 0, 2, 8 });
            p.AddCoeff(-2.8294354413935000e+01, new int[] { 0, 4, 0 });
            p.AddCoeff(1.0185967589016600e+03, new int[] { 0, 4, 2 });
            p.AddCoeff(-5.6022821739591300e+03, new int[] { 0, 4, 4 });
            p.AddCoeff(9.7106224348624800e+03, new int[] { 0, 4, 6 });
            p.AddCoeff(-5.2021191615334600e+03, new int[] { 0, 4, 8 });
            p.AddCoeff(2.0749193236885600e+01, new int[] { 0, 6, 0 });
            p.AddCoeff(-7.4697095652788300e+02, new int[] { 0, 6, 2 });
            p.AddCoeff(4.1083402609033600e+03, new int[] { 0, 6, 4 });
            p.AddCoeff(-7.1211231188991500e+03, new int[] { 0, 6, 6 });
            p.AddCoeff(3.8148873851245400e+03, new int[] { 0, 6, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("dd02b279-f01e-4bee-9a78-57957b4f3d81"));
            OrthonormalPolynomials[567] = p;
            p.AddCoeff(2.5377123250591500e+01, new int[] { 0, 1, 1 });
            p.AddCoeff(-2.2839410925532400e+02, new int[] { 0, 1, 3 });
            p.AddCoeff(5.0246704036171200e+02, new int[] { 0, 1, 5 });
            p.AddCoeff(-3.1105102498582200e+02, new int[] { 0, 1, 7 });
            p.AddCoeff(-2.2839410925532400e+02, new int[] { 0, 3, 1 });
            p.AddCoeff(2.0555469832979100e+03, new int[] { 0, 3, 3 });
            p.AddCoeff(-4.5222033632554100e+03, new int[] { 0, 3, 5 });
            p.AddCoeff(2.7994592248724000e+03, new int[] { 0, 3, 7 });
            p.AddCoeff(5.0246704036171200e+02, new int[] { 0, 5, 1 });
            p.AddCoeff(-4.5222033632554100e+03, new int[] { 0, 5, 3 });
            p.AddCoeff(9.9488473991619000e+03, new int[] { 0, 5, 5 });
            p.AddCoeff(-6.1588102947192700e+03, new int[] { 0, 5, 7 });
            p.AddCoeff(-3.1105102498582200e+02, new int[] { 0, 7, 1 });
            p.AddCoeff(2.7994592248724000e+03, new int[] { 0, 7, 3 });
            p.AddCoeff(-6.1588102947192700e+03, new int[] { 0, 7, 5 });
            p.AddCoeff(3.8125968491119300e+03, new int[] { 0, 7, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("39a9f762-4b57-41b6-b44f-9994ed0aa0b6"));
            OrthonormalPolynomials[568] = p;
            p.AddCoeff(-4.4911673672912700e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(9.4314514713116600e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-2.8294354413935000e+01, new int[] { 0, 0, 4 });
            p.AddCoeff(2.0749193236885600e+01, new int[] { 0, 0, 6 });
            p.AddCoeff(1.6168202522248600e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(-3.3953225296722000e+02, new int[] { 0, 2, 2 });
            p.AddCoeff(1.0185967589016600e+03, new int[] { 0, 2, 4 });
            p.AddCoeff(-7.4697095652788300e+02, new int[] { 0, 2, 6 });
            p.AddCoeff(-8.8925113872367100e+01, new int[] { 0, 4, 0 });
            p.AddCoeff(1.8674273913197100e+03, new int[] { 0, 4, 2 });
            p.AddCoeff(-5.6022821739591300e+03, new int[] { 0, 4, 4 });
            p.AddCoeff(4.1083402609033600e+03, new int[] { 0, 4, 6 });
            p.AddCoeff(1.5413686404543600e+02, new int[] { 0, 6, 0 });
            p.AddCoeff(-3.2368741449541600e+03, new int[] { 0, 6, 2 });
            p.AddCoeff(9.7106224348624800e+03, new int[] { 0, 6, 4 });
            p.AddCoeff(-7.1211231188991500e+03, new int[] { 0, 6, 6 });
            p.AddCoeff(-8.2573320024340800e+01, new int[] { 0, 8, 0 });
            p.AddCoeff(1.7340397205111600e+03, new int[] { 0, 8, 2 });
            p.AddCoeff(-5.2021191615334600e+03, new int[] { 0, 8, 4 });
            p.AddCoeff(3.8148873851245400e+03, new int[] { 0, 8, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("76299edb-5155-46a6-a06f-84da9c4eb1dc"));
            OrthonormalPolynomials[569] = p;
            p.AddCoeff(2.3584680961604700e+01, new int[] { 0, 1, 1 });
            p.AddCoeff(-1.1006184448748800e+02, new int[] { 0, 1, 3 });
            p.AddCoeff(9.9055660038739600e+01, new int[] { 0, 1, 5 });
            p.AddCoeff(-3.4590865410353500e+02, new int[] { 0, 3, 1 });
            p.AddCoeff(1.6142403858165000e+03, new int[] { 0, 3, 3 });
            p.AddCoeff(-1.4528163472348500e+03, new int[] { 0, 3, 5 });
            p.AddCoeff(1.3490437510037900e+03, new int[] { 0, 5, 1 });
            p.AddCoeff(-6.2955375046843400e+03, new int[] { 0, 5, 3 });
            p.AddCoeff(5.6659837542159100e+03, new int[] { 0, 5, 5 });
            p.AddCoeff(-1.9272053585768400e+03, new int[] { 0, 7, 1 });
            p.AddCoeff(8.9936250066919200e+03, new int[] { 0, 7, 3 });
            p.AddCoeff(-8.0942625060227200e+03, new int[] { 0, 7, 5 });
            p.AddCoeff(9.1006919710573000e+02, new int[] { 0, 9, 1 });
            p.AddCoeff(-4.2469895864934000e+03, new int[] { 0, 9, 3 });
            p.AddCoeff(3.8222906278440600e+03, new int[] { 0, 9, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c2564e50-1de3-4e27-8b8d-eab174c9754d"));
            OrthonormalPolynomials[570] = p;
            p.AddCoeff(-4.4855712597622800e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(4.4855712597622800e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-5.2331664697226600e+00, new int[] { 0, 0, 4 });
            p.AddCoeff(2.4670641928692500e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(-2.4670641928692500e+02, new int[] { 0, 2, 2 });
            p.AddCoeff(2.8782415583474600e+02, new int[] { 0, 2, 4 });
            p.AddCoeff(-2.1381223004866800e+02, new int[] { 0, 4, 0 });
            p.AddCoeff(2.1381223004866800e+03, new int[] { 0, 4, 2 });
            p.AddCoeff(-2.4944760172344700e+03, new int[] { 0, 4, 4 });
            p.AddCoeff(6.4143669014600500e+02, new int[] { 0, 6, 0 });
            p.AddCoeff(-6.4143669014600500e+03, new int[] { 0, 6, 2 });
            p.AddCoeff(7.4834280517033600e+03, new int[] { 0, 6, 4 });
            p.AddCoeff(-7.7888740946301000e+02, new int[] { 0, 8, 0 });
            p.AddCoeff(7.7888740946301000e+03, new int[] { 0, 8, 2 });
            p.AddCoeff(-9.0870197770684300e+03, new int[] { 0, 8, 4 });
            p.AddCoeff(3.2886357288438100e+02, new int[] { 0, 10, 0 });
            p.AddCoeff(-3.2886357288438100e+03, new int[] { 0, 10, 2 });
            p.AddCoeff(3.8367416836511000e+03, new int[] { 0, 10, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("1a855da8-1885-4514-b40e-a0c323f9d495"));
            OrthonormalPolynomials[571] = p;
            p.AddCoeff(1.8215977151856400e+01, new int[] { 0, 1, 1 });
            p.AddCoeff(-3.0359961919760700e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(-3.9467950495688900e+02, new int[] { 0, 3, 1 });
            p.AddCoeff(6.5779917492814900e+02, new int[] { 0, 3, 3 });
            p.AddCoeff(2.3680770297413400e+03, new int[] { 0, 5, 1 });
            p.AddCoeff(-3.9467950495688900e+03, new int[] { 0, 5, 3 });
            p.AddCoeff(-5.7510442150861000e+03, new int[] { 0, 7, 1 });
            p.AddCoeff(9.5850736918101700e+03, new int[] { 0, 7, 3 });
            p.AddCoeff(6.0705466714797800e+03, new int[] { 0, 9, 1 });
            p.AddCoeff(-1.0117577785799600e+04, new int[] { 0, 9, 3 });
            p.AddCoeff(-2.3178450927468200e+03, new int[] { 0, 11, 1 });
            p.AddCoeff(3.8630751545780400e+03, new int[] { 0, 11, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("dcb3f4e2-6364-4e2e-9663-9293e0477fd3"));
            OrthonormalPolynomials[572] = p;
            p.AddCoeff(-4.4585335662774400e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.3375600698832300e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(3.4776561816964000e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(-1.0432968545089200e+02, new int[] { 0, 2, 2 });
            p.AddCoeff(-4.3470702271204900e+02, new int[] { 0, 4, 0 });
            p.AddCoeff(1.3041210681361500e+03, new int[] { 0, 4, 2 });
            p.AddCoeff(1.9706718362946300e+03, new int[] { 0, 6, 0 });
            p.AddCoeff(-5.9120155088838900e+03, new int[] { 0, 6, 2 });
            p.AddCoeff(-4.0117248095997700e+03, new int[] { 0, 8, 0 });
            p.AddCoeff(1.2035174428799300e+04, new int[] { 0, 8, 2 });
            p.AddCoeff(3.7442764889597800e+03, new int[] { 0, 10, 0 });
            p.AddCoeff(-1.1232829466879400e+04, new int[] { 0, 10, 2 });
            p.AddCoeff(-1.3048236249405400e+03, new int[] { 0, 12, 0 });
            p.AddCoeff(3.9144708748216000e+03, new int[] { 0, 12, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d9c9b01f-a21a-4809-8e56-5cc587293643"));
            OrthonormalPolynomials[573] = p;
            p.AddCoeff(9.3315307495746500e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-2.7994592248724000e+02, new int[] { 0, 3, 1 });
            p.AddCoeff(2.3795403411415300e+03, new int[] { 0, 5, 1 });
            p.AddCoeff(-8.6116698060360400e+03, new int[] { 0, 7, 1 });
            p.AddCoeff(1.5070422160563000e+04, new int[] { 0, 9, 1 });
            p.AddCoeff(-1.2604353079743700e+04, new int[] { 0, 11, 1 });
            p.AddCoeff(4.0398567563281000e+03, new int[] { 0, 13, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e8864d99-10cf-4b63-813f-b9d45b5dc2f9"));
            OrthonormalPolynomials[574] = p;
            p.AddCoeff(-3.9882405547065600e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(4.1876525824418900e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(-7.1190093901512200e+02, new int[] { 0, 4, 0 });
            p.AddCoeff(4.5087059470957700e+03, new int[] { 0, 6, 0 });
            p.AddCoeff(-1.3526117841287300e+04, new int[] { 0, 8, 0 });
            p.AddCoeff(2.0740047356640600e+04, new int[] { 0, 10, 0 });
            p.AddCoeff(-1.5712157088364100e+04, new int[] { 0, 12, 0 });
            p.AddCoeff(4.6618488064376900e+03, new int[] { 0, 14, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("440a646d-a66f-4c22-8794-8f0d58c3301c"));
            OrthonormalPolynomials[575] = p;
            p.AddCoeff(9.3315307495746500e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-2.7994592248724000e+02, new int[] { 1, 0, 3 });
            p.AddCoeff(2.3795403411415300e+03, new int[] { 1, 0, 5 });
            p.AddCoeff(-8.6116698060360400e+03, new int[] { 1, 0, 7 });
            p.AddCoeff(1.5070422160563000e+04, new int[] { 1, 0, 9 });
            p.AddCoeff(-1.2604353079743700e+04, new int[] { 1, 0, 11 });
            p.AddCoeff(4.0398567563281000e+03, new int[] { 1, 0, 13 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("2fc26e3e-18f2-4f4c-88a6-8f9f63a8e527"));
            OrthonormalPolynomials[576] = p;
            p.AddCoeff(1.1963500960993100e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-9.3315307495746500e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(1.1664413436968300e+03, new int[] { 1, 1, 4 });
            p.AddCoeff(-5.2878674247589700e+03, new int[] { 1, 1, 6 });
            p.AddCoeff(1.0764587257545100e+04, new int[] { 1, 1, 8 });
            p.AddCoeff(-1.0046948107042000e+04, new int[] { 1, 1, 10 });
            p.AddCoeff(3.5012091888176700e+03, new int[] { 1, 1, 12 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("1906b5a7-5e17-4bd1-b95d-c0978d0877e8"));
            OrthonormalPolynomials[577] = p;
            p.AddCoeff(8.8884867156627400e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-1.9258387883935900e+02, new int[] { 1, 0, 3 });
            p.AddCoeff(1.1555032730361600e+03, new int[] { 1, 0, 5 });
            p.AddCoeff(-2.8062222345163800e+03, new int[] { 1, 0, 7 });
            p.AddCoeff(2.9621234697672900e+03, new int[] { 1, 0, 9 });
            p.AddCoeff(-1.1309925975475100e+03, new int[] { 1, 0, 11 });
            p.AddCoeff(-2.6665460146988200e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(5.7775163651807800e+02, new int[] { 1, 2, 3 });
            p.AddCoeff(-3.4665098191084700e+03, new int[] { 1, 2, 5 });
            p.AddCoeff(8.4186667035491400e+03, new int[] { 1, 2, 7 });
            p.AddCoeff(-8.8863704093018700e+03, new int[] { 1, 2, 9 });
            p.AddCoeff(3.3929777926425300e+03, new int[] { 1, 2, 11 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("78e4e740-696a-4369-b5ae-d8ab700a5866"));
            OrthonormalPolynomials[578] = p;
            p.AddCoeff(2.7407293110638800e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-1.5074011210851400e+02, new int[] { 1, 1, 2 });
            p.AddCoeff(1.3064143049404500e+03, new int[] { 1, 1, 4 });
            p.AddCoeff(-3.9192429148213500e+03, new int[] { 1, 1, 6 });
            p.AddCoeff(4.7590806822830700e+03, new int[] { 1, 1, 8 });
            p.AddCoeff(-2.0093896214084100e+03, new int[] { 1, 1, 10 });
            p.AddCoeff(-4.5678821851064700e+00, new int[] { 1, 3, 0 });
            p.AddCoeff(2.5123352018085600e+02, new int[] { 1, 3, 2 });
            p.AddCoeff(-2.1773571749007500e+03, new int[] { 1, 3, 4 });
            p.AddCoeff(6.5320715247022600e+03, new int[] { 1, 3, 6 });
            p.AddCoeff(-7.9318011371384500e+03, new int[] { 1, 3, 8 });
            p.AddCoeff(3.3489827023473500e+03, new int[] { 1, 3, 10 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("0a414f26-e82d-40ca-a969-9b8ff244521c"));
            OrthonormalPolynomials[579] = p;
            p.AddCoeff(7.3900187608663900e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-1.0838694182604000e+02, new int[] { 1, 0, 3 });
            p.AddCoeff(4.2270907312155800e+02, new int[] { 1, 0, 5 });
            p.AddCoeff(-6.0387010445936800e+02, new int[] { 1, 0, 7 });
            p.AddCoeff(2.8516088266136800e+02, new int[] { 1, 0, 9 });
            p.AddCoeff(-7.3900187608663900e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(1.0838694182604000e+03, new int[] { 1, 2, 3 });
            p.AddCoeff(-4.2270907312155800e+03, new int[] { 1, 2, 5 });
            p.AddCoeff(6.0387010445936800e+03, new int[] { 1, 2, 7 });
            p.AddCoeff(-2.8516088266136800e+03, new int[] { 1, 2, 9 });
            p.AddCoeff(8.6216885543441200e+01, new int[] { 1, 4, 1 });
            p.AddCoeff(-1.2645143213038000e+03, new int[] { 1, 4, 3 });
            p.AddCoeff(4.9316058530848400e+03, new int[] { 1, 4, 5 });
            p.AddCoeff(-7.0451512186926300e+03, new int[] { 1, 4, 7 });
            p.AddCoeff(3.3268769643826300e+03, new int[] { 1, 4, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b263aeba-e56c-4266-a9a0-e50b9e2c12a0"));
            OrthonormalPolynomials[580] = p;
            p.AddCoeff(4.2933449549966900e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-1.5456041837988100e+02, new int[] { 1, 1, 2 });
            p.AddCoeff(8.5008230108934600e+02, new int[] { 1, 1, 4 });
            p.AddCoeff(-1.4734759885548700e+03, new int[] { 1, 1, 6 });
            p.AddCoeff(7.8936213672582100e+02, new int[] { 1, 1, 8 });
            p.AddCoeff(-2.0035609789984600e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(7.2128195243944500e+02, new int[] { 1, 3, 2 });
            p.AddCoeff(-3.9670507384169500e+03, new int[] { 1, 3, 4 });
            p.AddCoeff(6.8762212799227100e+03, new int[] { 1, 3, 6 });
            p.AddCoeff(-3.6836899713871600e+03, new int[] { 1, 3, 8 });
            p.AddCoeff(1.8032048810986100e+01, new int[] { 1, 5, 0 });
            p.AddCoeff(-6.4915375719550000e+02, new int[] { 1, 5, 2 });
            p.AddCoeff(3.5703456645752500e+03, new int[] { 1, 5, 4 });
            p.AddCoeff(-6.1885991519304400e+03, new int[] { 1, 5, 6 });
            p.AddCoeff(3.3153209742484500e+03, new int[] { 1, 5, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("9dbec04d-2968-45b5-9c50-6154f89294fb"));
            OrthonormalPolynomials[581] = p;
            p.AddCoeff(5.8456259587602400e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-5.2610633628842200e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(1.1574339398345300e+02, new int[] { 1, 0, 5 });
            p.AddCoeff(-7.1650672465946900e+01, new int[] { 1, 0, 7 });
            p.AddCoeff(-1.2275814513396500e+02, new int[] { 1, 2, 1 });
            p.AddCoeff(1.1048233062056900e+03, new int[] { 1, 2, 3 });
            p.AddCoeff(-2.4306112736525100e+03, new int[] { 1, 2, 5 });
            p.AddCoeff(1.5046641217848900e+03, new int[] { 1, 2, 7 });
            p.AddCoeff(3.6827443540189500e+02, new int[] { 1, 4, 1 });
            p.AddCoeff(-3.3144699186170600e+03, new int[] { 1, 4, 3 });
            p.AddCoeff(7.2918338209575200e+03, new int[] { 1, 4, 5 });
            p.AddCoeff(-4.5139923653546600e+03, new int[] { 1, 4, 7 });
            p.AddCoeff(-2.7006791929472300e+02, new int[] { 1, 6, 1 });
            p.AddCoeff(2.4306112736525100e+03, new int[] { 1, 6, 3 });
            p.AddCoeff(-5.3473448020355200e+03, new int[] { 1, 6, 5 });
            p.AddCoeff(3.3102610679267500e+03, new int[] { 1, 6, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("fd1417a2-b5cf-4755-bb8f-c982268d7b30"));
            OrthonormalPolynomials[582] = p;
            p.AddCoeff(5.8456259587602400e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-1.2275814513396500e+02, new int[] { 1, 1, 2 });
            p.AddCoeff(3.6827443540189500e+02, new int[] { 1, 1, 4 });
            p.AddCoeff(-2.7006791929472300e+02, new int[] { 1, 1, 6 });
            p.AddCoeff(-5.2610633628842200e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(1.1048233062056900e+03, new int[] { 1, 3, 2 });
            p.AddCoeff(-3.3144699186170600e+03, new int[] { 1, 3, 4 });
            p.AddCoeff(2.4306112736525100e+03, new int[] { 1, 3, 6 });
            p.AddCoeff(1.1574339398345300e+02, new int[] { 1, 5, 0 });
            p.AddCoeff(-2.4306112736525100e+03, new int[] { 1, 5, 2 });
            p.AddCoeff(7.2918338209575200e+03, new int[] { 1, 5, 4 });
            p.AddCoeff(-5.3473448020355200e+03, new int[] { 1, 5, 6 });
            p.AddCoeff(-7.1650672465946900e+01, new int[] { 1, 7, 0 });
            p.AddCoeff(1.5046641217848900e+03, new int[] { 1, 7, 2 });
            p.AddCoeff(-4.5139923653546600e+03, new int[] { 1, 7, 4 });
            p.AddCoeff(3.3102610679267500e+03, new int[] { 1, 7, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d012334e-f12e-47d5-9c91-5a98fc770e50"));
            OrthonormalPolynomials[583] = p;
            p.AddCoeff(4.2933449549966900e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-2.0035609789984600e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(1.8032048810986100e+01, new int[] { 1, 0, 5 });
            p.AddCoeff(-1.5456041837988100e+02, new int[] { 1, 2, 1 });
            p.AddCoeff(7.2128195243944500e+02, new int[] { 1, 2, 3 });
            p.AddCoeff(-6.4915375719550000e+02, new int[] { 1, 2, 5 });
            p.AddCoeff(8.5008230108934600e+02, new int[] { 1, 4, 1 });
            p.AddCoeff(-3.9670507384169500e+03, new int[] { 1, 4, 3 });
            p.AddCoeff(3.5703456645752500e+03, new int[] { 1, 4, 5 });
            p.AddCoeff(-1.4734759885548700e+03, new int[] { 1, 6, 1 });
            p.AddCoeff(6.8762212799227100e+03, new int[] { 1, 6, 3 });
            p.AddCoeff(-6.1885991519304400e+03, new int[] { 1, 6, 5 });
            p.AddCoeff(7.8936213672582100e+02, new int[] { 1, 8, 1 });
            p.AddCoeff(-3.6836899713871600e+03, new int[] { 1, 8, 3 });
            p.AddCoeff(3.3153209742484500e+03, new int[] { 1, 8, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c8847b12-9086-402f-8620-6eb3be9c3a61"));
            OrthonormalPolynomials[584] = p;
            p.AddCoeff(7.3900187608663900e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-7.3900187608663900e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(8.6216885543441200e+01, new int[] { 1, 1, 4 });
            p.AddCoeff(-1.0838694182604000e+02, new int[] { 1, 3, 0 });
            p.AddCoeff(1.0838694182604000e+03, new int[] { 1, 3, 2 });
            p.AddCoeff(-1.2645143213038000e+03, new int[] { 1, 3, 4 });
            p.AddCoeff(4.2270907312155800e+02, new int[] { 1, 5, 0 });
            p.AddCoeff(-4.2270907312155800e+03, new int[] { 1, 5, 2 });
            p.AddCoeff(4.9316058530848400e+03, new int[] { 1, 5, 4 });
            p.AddCoeff(-6.0387010445936800e+02, new int[] { 1, 7, 0 });
            p.AddCoeff(6.0387010445936800e+03, new int[] { 1, 7, 2 });
            p.AddCoeff(-7.0451512186926300e+03, new int[] { 1, 7, 4 });
            p.AddCoeff(2.8516088266136800e+02, new int[] { 1, 9, 0 });
            p.AddCoeff(-2.8516088266136800e+03, new int[] { 1, 9, 2 });
            p.AddCoeff(3.3268769643826300e+03, new int[] { 1, 9, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("2b197b57-30be-4886-ae67-9d3bed8bdfa6"));
            OrthonormalPolynomials[585] = p;
            p.AddCoeff(2.7407293110638800e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-4.5678821851064700e+00, new int[] { 1, 0, 3 });
            p.AddCoeff(-1.5074011210851400e+02, new int[] { 1, 2, 1 });
            p.AddCoeff(2.5123352018085600e+02, new int[] { 1, 2, 3 });
            p.AddCoeff(1.3064143049404500e+03, new int[] { 1, 4, 1 });
            p.AddCoeff(-2.1773571749007500e+03, new int[] { 1, 4, 3 });
            p.AddCoeff(-3.9192429148213500e+03, new int[] { 1, 6, 1 });
            p.AddCoeff(6.5320715247022600e+03, new int[] { 1, 6, 3 });
            p.AddCoeff(4.7590806822830700e+03, new int[] { 1, 8, 1 });
            p.AddCoeff(-7.9318011371384500e+03, new int[] { 1, 8, 3 });
            p.AddCoeff(-2.0093896214084100e+03, new int[] { 1, 10, 1 });
            p.AddCoeff(3.3489827023473500e+03, new int[] { 1, 10, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b40f5be7-96b5-4a6c-ae96-42430f0abaad"));
            OrthonormalPolynomials[586] = p;
            p.AddCoeff(8.8884867156627400e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-2.6665460146988200e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(-1.9258387883935900e+02, new int[] { 1, 3, 0 });
            p.AddCoeff(5.7775163651807800e+02, new int[] { 1, 3, 2 });
            p.AddCoeff(1.1555032730361600e+03, new int[] { 1, 5, 0 });
            p.AddCoeff(-3.4665098191084700e+03, new int[] { 1, 5, 2 });
            p.AddCoeff(-2.8062222345163800e+03, new int[] { 1, 7, 0 });
            p.AddCoeff(8.4186667035491400e+03, new int[] { 1, 7, 2 });
            p.AddCoeff(2.9621234697672900e+03, new int[] { 1, 9, 0 });
            p.AddCoeff(-8.8863704093018700e+03, new int[] { 1, 9, 2 });
            p.AddCoeff(-1.1309925975475100e+03, new int[] { 1, 11, 0 });
            p.AddCoeff(3.3929777926425300e+03, new int[] { 1, 11, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b967a6ed-8dcc-45e2-80be-654ea05ca9cc"));
            OrthonormalPolynomials[587] = p;
            p.AddCoeff(1.1963500960993100e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-9.3315307495746500e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(1.1664413436968300e+03, new int[] { 1, 4, 1 });
            p.AddCoeff(-5.2878674247589700e+03, new int[] { 1, 6, 1 });
            p.AddCoeff(1.0764587257545100e+04, new int[] { 1, 8, 1 });
            p.AddCoeff(-1.0046948107042000e+04, new int[] { 1, 10, 1 });
            p.AddCoeff(3.5012091888176700e+03, new int[] { 1, 12, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("36daf47e-09f2-4872-b3ce-e6e18c56c959"));
            OrthonormalPolynomials[588] = p;
            p.AddCoeff(9.3315307495746500e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-2.7994592248724000e+02, new int[] { 1, 3, 0 });
            p.AddCoeff(2.3795403411415300e+03, new int[] { 1, 5, 0 });
            p.AddCoeff(-8.6116698060360400e+03, new int[] { 1, 7, 0 });
            p.AddCoeff(1.5070422160563000e+04, new int[] { 1, 9, 0 });
            p.AddCoeff(-1.2604353079743700e+04, new int[] { 1, 11, 0 });
            p.AddCoeff(4.0398567563281000e+03, new int[] { 1, 13, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("96620cc5-0c88-4f92-8cfc-44c2bb1ca41e"));
            OrthonormalPolynomials[589] = p;
            p.AddCoeff(-4.4585335662774400e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(3.4776561816964000e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(-4.3470702271204900e+02, new int[] { 0, 0, 4 });
            p.AddCoeff(1.9706718362946300e+03, new int[] { 0, 0, 6 });
            p.AddCoeff(-4.0117248095997700e+03, new int[] { 0, 0, 8 });
            p.AddCoeff(3.7442764889597800e+03, new int[] { 0, 0, 10 });
            p.AddCoeff(-1.3048236249405400e+03, new int[] { 0, 0, 12 });
            p.AddCoeff(1.3375600698832300e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-1.0432968545089200e+02, new int[] { 2, 0, 2 });
            p.AddCoeff(1.3041210681361500e+03, new int[] { 2, 0, 4 });
            p.AddCoeff(-5.9120155088838900e+03, new int[] { 2, 0, 6 });
            p.AddCoeff(1.2035174428799300e+04, new int[] { 2, 0, 8 });
            p.AddCoeff(-1.1232829466879400e+04, new int[] { 2, 0, 10 });
            p.AddCoeff(3.9144708748216000e+03, new int[] { 2, 0, 12 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("bea06f59-b6d0-46bc-af1d-e26966ccca61"));
            OrthonormalPolynomials[590] = p;
            p.AddCoeff(8.8884867156627400e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-1.9258387883935900e+02, new int[] { 0, 1, 3 });
            p.AddCoeff(1.1555032730361600e+03, new int[] { 0, 1, 5 });
            p.AddCoeff(-2.8062222345163800e+03, new int[] { 0, 1, 7 });
            p.AddCoeff(2.9621234697672900e+03, new int[] { 0, 1, 9 });
            p.AddCoeff(-1.1309925975475100e+03, new int[] { 0, 1, 11 });
            p.AddCoeff(-2.6665460146988200e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(5.7775163651807800e+02, new int[] { 2, 1, 3 });
            p.AddCoeff(-3.4665098191084700e+03, new int[] { 2, 1, 5 });
            p.AddCoeff(8.4186667035491400e+03, new int[] { 2, 1, 7 });
            p.AddCoeff(-8.8863704093018700e+03, new int[] { 2, 1, 9 });
            p.AddCoeff(3.3929777926425300e+03, new int[] { 2, 1, 11 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a5dd5208-0314-493f-a9f3-1cffec389d01"));
            OrthonormalPolynomials[591] = p;
            p.AddCoeff(-4.9839680664025300e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(2.7411824365213900e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(-2.3756914449852100e+02, new int[] { 0, 0, 4 });
            p.AddCoeff(7.1270743349556200e+02, new int[] { 0, 0, 6 });
            p.AddCoeff(-8.6543045495889600e+02, new int[] { 0, 0, 8 });
            p.AddCoeff(3.6540396987153400e+02, new int[] { 0, 0, 10 });
            p.AddCoeff(1.4951904199207600e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-8.2235473095641700e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(7.1270743349556200e+02, new int[] { 0, 2, 4 });
            p.AddCoeff(-2.1381223004866800e+03, new int[] { 0, 2, 6 });
            p.AddCoeff(2.5962913648766900e+03, new int[] { 0, 2, 8 });
            p.AddCoeff(-1.0962119096146000e+03, new int[] { 0, 2, 10 });
            p.AddCoeff(1.4951904199207600e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-8.2235473095641700e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(7.1270743349556200e+02, new int[] { 2, 0, 4 });
            p.AddCoeff(-2.1381223004866800e+03, new int[] { 2, 0, 6 });
            p.AddCoeff(2.5962913648766900e+03, new int[] { 2, 0, 8 });
            p.AddCoeff(-1.0962119096146000e+03, new int[] { 2, 0, 10 });
            p.AddCoeff(-4.4855712597622800e+00, new int[] { 2, 2, 0 });
            p.AddCoeff(2.4670641928692500e+02, new int[] { 2, 2, 2 });
            p.AddCoeff(-2.1381223004866800e+03, new int[] { 2, 2, 4 });
            p.AddCoeff(6.4143669014600500e+03, new int[] { 2, 2, 6 });
            p.AddCoeff(-7.7888740946301000e+03, new int[] { 2, 2, 8 });
            p.AddCoeff(3.2886357288438100e+03, new int[] { 2, 2, 10 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a4c6d991-a3f9-4625-b46d-c227997a9651"));
            OrthonormalPolynomials[592] = p;
            p.AddCoeff(1.6827812978247900e+01, new int[] { 0, 1, 1 });
            p.AddCoeff(-2.4680792368097000e+02, new int[] { 0, 1, 3 });
            p.AddCoeff(9.6255090235578100e+02, new int[] { 0, 1, 5 });
            p.AddCoeff(-1.3750727176511200e+03, new int[] { 0, 1, 7 });
            p.AddCoeff(6.4933989444636000e+02, new int[] { 0, 1, 9 });
            p.AddCoeff(-2.8046354963746500e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(4.1134653946828300e+02, new int[] { 0, 3, 3 });
            p.AddCoeff(-1.6042515039263000e+03, new int[] { 0, 3, 5 });
            p.AddCoeff(2.2917878627518600e+03, new int[] { 0, 3, 7 });
            p.AddCoeff(-1.0822331574106000e+03, new int[] { 0, 3, 9 });
            p.AddCoeff(-5.0483438934743800e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(7.4042377104290900e+02, new int[] { 2, 1, 3 });
            p.AddCoeff(-2.8876527070673400e+03, new int[] { 2, 1, 5 });
            p.AddCoeff(4.1252181529533500e+03, new int[] { 2, 1, 7 });
            p.AddCoeff(-1.9480196833390800e+03, new int[] { 2, 1, 9 });
            p.AddCoeff(8.4139064891239600e+01, new int[] { 2, 3, 1 });
            p.AddCoeff(-1.2340396184048500e+03, new int[] { 2, 3, 3 });
            p.AddCoeff(4.8127545117789100e+03, new int[] { 2, 3, 5 });
            p.AddCoeff(-6.8753635882555800e+03, new int[] { 2, 3, 7 });
            p.AddCoeff(3.2466994722318000e+03, new int[] { 2, 3, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("91715b63-b90b-4adb-b8b3-ccfb4b9286a1"));
            OrthonormalPolynomials[593] = p;
            p.AddCoeff(-5.0135467715791900e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.8048768377685100e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(-9.9268226077267900e+01, new int[] { 0, 0, 4 });
            p.AddCoeff(1.7206492520059800e+02, new int[] { 0, 0, 6 });
            p.AddCoeff(-9.2177638500320200e+01, new int[] { 0, 0, 8 });
            p.AddCoeff(5.0135467715791900e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-1.8048768377685100e+02, new int[] { 0, 2, 2 });
            p.AddCoeff(9.9268226077267900e+02, new int[] { 0, 2, 4 });
            p.AddCoeff(-1.7206492520059800e+03, new int[] { 0, 2, 6 });
            p.AddCoeff(9.2177638500320200e+02, new int[] { 0, 2, 8 });
            p.AddCoeff(-5.8491379001757200e+00, new int[] { 0, 4, 0 });
            p.AddCoeff(2.1056896440632600e+02, new int[] { 0, 4, 2 });
            p.AddCoeff(-1.1581293042347900e+03, new int[] { 0, 4, 4 });
            p.AddCoeff(2.0074241273403100e+03, new int[] { 0, 4, 6 });
            p.AddCoeff(-1.0754057825037400e+03, new int[] { 0, 4, 8 });
            p.AddCoeff(1.5040640314737600e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-5.4146305133055200e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(2.9780467823180400e+02, new int[] { 2, 0, 4 });
            p.AddCoeff(-5.1619477560179300e+02, new int[] { 2, 0, 6 });
            p.AddCoeff(2.7653291550096100e+02, new int[] { 2, 0, 8 });
            p.AddCoeff(-1.5040640314737600e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(5.4146305133055200e+02, new int[] { 2, 2, 2 });
            p.AddCoeff(-2.9780467823180400e+03, new int[] { 2, 2, 4 });
            p.AddCoeff(5.1619477560179300e+03, new int[] { 2, 2, 6 });
            p.AddCoeff(-2.7653291550096100e+03, new int[] { 2, 2, 8 });
            p.AddCoeff(1.7547413700527200e+01, new int[] { 2, 4, 0 });
            p.AddCoeff(-6.3170689321897700e+02, new int[] { 2, 4, 2 });
            p.AddCoeff(3.4743879127043800e+03, new int[] { 2, 4, 4 });
            p.AddCoeff(-6.0222723820209200e+03, new int[] { 2, 4, 6 });
            p.AddCoeff(3.2262173475112000e+03, new int[] { 2, 4, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("428a5978-f4ea-4490-a50a-541402abeace"));
            OrthonormalPolynomials[594] = p;
            p.AddCoeff(2.0825782043134200e+01, new int[] { 0, 1, 1 });
            p.AddCoeff(-1.8743203838820800e+02, new int[] { 0, 1, 3 });
            p.AddCoeff(4.1235048445405700e+02, new int[] { 0, 1, 5 });
            p.AddCoeff(-2.5526458561441600e+02, new int[] { 0, 1, 7 });
            p.AddCoeff(-9.7186982867959500e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(8.7468284581163500e+02, new int[] { 0, 3, 3 });
            p.AddCoeff(-1.9243022607856000e+03, new int[] { 0, 3, 5 });
            p.AddCoeff(1.1912347328672700e+03, new int[] { 0, 3, 7 });
            p.AddCoeff(8.7468284581163500e+01, new int[] { 0, 5, 1 });
            p.AddCoeff(-7.8721456123047200e+02, new int[] { 0, 5, 3 });
            p.AddCoeff(1.7318720347070400e+03, new int[] { 0, 5, 5 });
            p.AddCoeff(-1.0721112595805500e+03, new int[] { 0, 5, 7 });
            p.AddCoeff(-6.2477346129402500e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(5.6229611516462300e+02, new int[] { 2, 1, 3 });
            p.AddCoeff(-1.2370514533621700e+03, new int[] { 2, 1, 5 });
            p.AddCoeff(7.6579375684324800e+02, new int[] { 2, 1, 7 });
            p.AddCoeff(2.9156094860387800e+02, new int[] { 2, 3, 1 });
            p.AddCoeff(-2.6240485374349100e+03, new int[] { 2, 3, 3 });
            p.AddCoeff(5.7729067823567900e+03, new int[] { 2, 3, 5 });
            p.AddCoeff(-3.5737041986018200e+03, new int[] { 2, 3, 7 });
            p.AddCoeff(-2.6240485374349100e+02, new int[] { 2, 5, 1 });
            p.AddCoeff(2.3616436836914200e+03, new int[] { 2, 5, 3 });
            p.AddCoeff(-5.1956161041211100e+03, new int[] { 2, 5, 5 });
            p.AddCoeff(3.2163337787416400e+03, new int[] { 2, 5, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c7d7c91a-3728-4f28-9aca-2b588c2356ed"));
            OrthonormalPolynomials[595] = p;
            p.AddCoeff(-5.0182628884508000e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.0538352065746700e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(-3.1615056197240000e+01, new int[] { 0, 0, 4 });
            p.AddCoeff(2.3184374544642700e+01, new int[] { 0, 0, 6 });
            p.AddCoeff(1.0538352065746700e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(-2.2130539338068000e+02, new int[] { 0, 2, 2 });
            p.AddCoeff(6.6391618014204100e+02, new int[] { 0, 2, 4 });
            p.AddCoeff(-4.8687186543749600e+02, new int[] { 0, 2, 6 });
            p.AddCoeff(-3.1615056197240000e+01, new int[] { 0, 4, 0 });
            p.AddCoeff(6.6391618014204100e+02, new int[] { 0, 4, 2 });
            p.AddCoeff(-1.9917485404261200e+03, new int[] { 0, 4, 4 });
            p.AddCoeff(1.4606155963124900e+03, new int[] { 0, 4, 6 });
            p.AddCoeff(2.3184374544642700e+01, new int[] { 0, 6, 0 });
            p.AddCoeff(-4.8687186543749600e+02, new int[] { 0, 6, 2 });
            p.AddCoeff(1.4606155963124900e+03, new int[] { 0, 6, 4 });
            p.AddCoeff(-1.0711181039624900e+03, new int[] { 0, 6, 6 });
            p.AddCoeff(1.5054788665352400e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-3.1615056197240000e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(9.4845168591720100e+01, new int[] { 2, 0, 4 });
            p.AddCoeff(-6.9553123633928100e+01, new int[] { 2, 0, 6 });
            p.AddCoeff(-3.1615056197240000e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(6.6391618014204100e+02, new int[] { 2, 2, 2 });
            p.AddCoeff(-1.9917485404261200e+03, new int[] { 2, 2, 4 });
            p.AddCoeff(1.4606155963124900e+03, new int[] { 2, 2, 6 });
            p.AddCoeff(9.4845168591720100e+01, new int[] { 2, 4, 0 });
            p.AddCoeff(-1.9917485404261200e+03, new int[] { 2, 4, 2 });
            p.AddCoeff(5.9752456212783700e+03, new int[] { 2, 4, 4 });
            p.AddCoeff(-4.3818467889374700e+03, new int[] { 2, 4, 6 });
            p.AddCoeff(-6.9553123633928100e+01, new int[] { 2, 6, 0 });
            p.AddCoeff(1.4606155963124900e+03, new int[] { 2, 6, 2 });
            p.AddCoeff(-4.3818467889374700e+03, new int[] { 2, 6, 4 });
            p.AddCoeff(3.2133543118874800e+03, new int[] { 2, 6, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("9efa3a1d-6027-4f39-b600-576a9022dca8"));
            OrthonormalPolynomials[596] = p;
            p.AddCoeff(2.0825782043134200e+01, new int[] { 0, 1, 1 });
            p.AddCoeff(-9.7186982867959500e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(8.7468284581163500e+01, new int[] { 0, 1, 5 });
            p.AddCoeff(-1.8743203838820800e+02, new int[] { 0, 3, 1 });
            p.AddCoeff(8.7468284581163500e+02, new int[] { 0, 3, 3 });
            p.AddCoeff(-7.8721456123047200e+02, new int[] { 0, 3, 5 });
            p.AddCoeff(4.1235048445405700e+02, new int[] { 0, 5, 1 });
            p.AddCoeff(-1.9243022607856000e+03, new int[] { 0, 5, 3 });
            p.AddCoeff(1.7318720347070400e+03, new int[] { 0, 5, 5 });
            p.AddCoeff(-2.5526458561441600e+02, new int[] { 0, 7, 1 });
            p.AddCoeff(1.1912347328672700e+03, new int[] { 0, 7, 3 });
            p.AddCoeff(-1.0721112595805500e+03, new int[] { 0, 7, 5 });
            p.AddCoeff(-6.2477346129402500e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(2.9156094860387800e+02, new int[] { 2, 1, 3 });
            p.AddCoeff(-2.6240485374349100e+02, new int[] { 2, 1, 5 });
            p.AddCoeff(5.6229611516462300e+02, new int[] { 2, 3, 1 });
            p.AddCoeff(-2.6240485374349100e+03, new int[] { 2, 3, 3 });
            p.AddCoeff(2.3616436836914200e+03, new int[] { 2, 3, 5 });
            p.AddCoeff(-1.2370514533621700e+03, new int[] { 2, 5, 1 });
            p.AddCoeff(5.7729067823567900e+03, new int[] { 2, 5, 3 });
            p.AddCoeff(-5.1956161041211100e+03, new int[] { 2, 5, 5 });
            p.AddCoeff(7.6579375684324800e+02, new int[] { 2, 7, 1 });
            p.AddCoeff(-3.5737041986018200e+03, new int[] { 2, 7, 3 });
            p.AddCoeff(3.2163337787416400e+03, new int[] { 2, 7, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("744a5855-575c-4bf0-8893-759c6dedceca"));
            OrthonormalPolynomials[597] = p;
            p.AddCoeff(-5.0135467715791900e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(5.0135467715791900e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-5.8491379001757200e+00, new int[] { 0, 0, 4 });
            p.AddCoeff(1.8048768377685100e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(-1.8048768377685100e+02, new int[] { 0, 2, 2 });
            p.AddCoeff(2.1056896440632600e+02, new int[] { 0, 2, 4 });
            p.AddCoeff(-9.9268226077267900e+01, new int[] { 0, 4, 0 });
            p.AddCoeff(9.9268226077267900e+02, new int[] { 0, 4, 2 });
            p.AddCoeff(-1.1581293042347900e+03, new int[] { 0, 4, 4 });
            p.AddCoeff(1.7206492520059800e+02, new int[] { 0, 6, 0 });
            p.AddCoeff(-1.7206492520059800e+03, new int[] { 0, 6, 2 });
            p.AddCoeff(2.0074241273403100e+03, new int[] { 0, 6, 4 });
            p.AddCoeff(-9.2177638500320200e+01, new int[] { 0, 8, 0 });
            p.AddCoeff(9.2177638500320200e+02, new int[] { 0, 8, 2 });
            p.AddCoeff(-1.0754057825037400e+03, new int[] { 0, 8, 4 });
            p.AddCoeff(1.5040640314737600e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-1.5040640314737600e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(1.7547413700527200e+01, new int[] { 2, 0, 4 });
            p.AddCoeff(-5.4146305133055200e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(5.4146305133055200e+02, new int[] { 2, 2, 2 });
            p.AddCoeff(-6.3170689321897700e+02, new int[] { 2, 2, 4 });
            p.AddCoeff(2.9780467823180400e+02, new int[] { 2, 4, 0 });
            p.AddCoeff(-2.9780467823180400e+03, new int[] { 2, 4, 2 });
            p.AddCoeff(3.4743879127043800e+03, new int[] { 2, 4, 4 });
            p.AddCoeff(-5.1619477560179300e+02, new int[] { 2, 6, 0 });
            p.AddCoeff(5.1619477560179300e+03, new int[] { 2, 6, 2 });
            p.AddCoeff(-6.0222723820209200e+03, new int[] { 2, 6, 4 });
            p.AddCoeff(2.7653291550096100e+02, new int[] { 2, 8, 0 });
            p.AddCoeff(-2.7653291550096100e+03, new int[] { 2, 8, 2 });
            p.AddCoeff(3.2262173475112000e+03, new int[] { 2, 8, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("28de035c-397f-4779-a5c6-40a78b4f4079"));
            OrthonormalPolynomials[598] = p;
            p.AddCoeff(1.6827812978247900e+01, new int[] { 0, 1, 1 });
            p.AddCoeff(-2.8046354963746500e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(-2.4680792368097000e+02, new int[] { 0, 3, 1 });
            p.AddCoeff(4.1134653946828300e+02, new int[] { 0, 3, 3 });
            p.AddCoeff(9.6255090235578100e+02, new int[] { 0, 5, 1 });
            p.AddCoeff(-1.6042515039263000e+03, new int[] { 0, 5, 3 });
            p.AddCoeff(-1.3750727176511200e+03, new int[] { 0, 7, 1 });
            p.AddCoeff(2.2917878627518600e+03, new int[] { 0, 7, 3 });
            p.AddCoeff(6.4933989444636000e+02, new int[] { 0, 9, 1 });
            p.AddCoeff(-1.0822331574106000e+03, new int[] { 0, 9, 3 });
            p.AddCoeff(-5.0483438934743800e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(8.4139064891239600e+01, new int[] { 2, 1, 3 });
            p.AddCoeff(7.4042377104290900e+02, new int[] { 2, 3, 1 });
            p.AddCoeff(-1.2340396184048500e+03, new int[] { 2, 3, 3 });
            p.AddCoeff(-2.8876527070673400e+03, new int[] { 2, 5, 1 });
            p.AddCoeff(4.8127545117789100e+03, new int[] { 2, 5, 3 });
            p.AddCoeff(4.1252181529533500e+03, new int[] { 2, 7, 1 });
            p.AddCoeff(-6.8753635882555800e+03, new int[] { 2, 7, 3 });
            p.AddCoeff(-1.9480196833390800e+03, new int[] { 2, 9, 1 });
            p.AddCoeff(3.2466994722318000e+03, new int[] { 2, 9, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("da164955-167e-4b7f-b02f-e1624bdd03e9"));
            OrthonormalPolynomials[599] = p;
            p.AddCoeff(-4.9839680664025300e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.4951904199207600e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(2.7411824365213900e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(-8.2235473095641700e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(-2.3756914449852100e+02, new int[] { 0, 4, 0 });
            p.AddCoeff(7.1270743349556200e+02, new int[] { 0, 4, 2 });
            p.AddCoeff(7.1270743349556200e+02, new int[] { 0, 6, 0 });
            p.AddCoeff(-2.1381223004866800e+03, new int[] { 0, 6, 2 });
            p.AddCoeff(-8.6543045495889600e+02, new int[] { 0, 8, 0 });
            p.AddCoeff(2.5962913648766900e+03, new int[] { 0, 8, 2 });
            p.AddCoeff(3.6540396987153400e+02, new int[] { 0, 10, 0 });
            p.AddCoeff(-1.0962119096146000e+03, new int[] { 0, 10, 2 });
            p.AddCoeff(1.4951904199207600e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-4.4855712597622800e+00, new int[] { 2, 0, 2 });
            p.AddCoeff(-8.2235473095641700e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(2.4670641928692500e+02, new int[] { 2, 2, 2 });
            p.AddCoeff(7.1270743349556200e+02, new int[] { 2, 4, 0 });
            p.AddCoeff(-2.1381223004866800e+03, new int[] { 2, 4, 2 });
            p.AddCoeff(-2.1381223004866800e+03, new int[] { 2, 6, 0 });
            p.AddCoeff(6.4143669014600500e+03, new int[] { 2, 6, 2 });
            p.AddCoeff(2.5962913648766900e+03, new int[] { 2, 8, 0 });
            p.AddCoeff(-7.7888740946301000e+03, new int[] { 2, 8, 2 });
            p.AddCoeff(-1.0962119096146000e+03, new int[] { 2, 10, 0 });
            p.AddCoeff(3.2886357288438100e+03, new int[] { 2, 10, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("aae95cb4-196a-418c-92f0-dc581ce01d0d"));
            OrthonormalPolynomials[600] = p;
            p.AddCoeff(8.8884867156627400e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-1.9258387883935900e+02, new int[] { 0, 3, 1 });
            p.AddCoeff(1.1555032730361600e+03, new int[] { 0, 5, 1 });
            p.AddCoeff(-2.8062222345163800e+03, new int[] { 0, 7, 1 });
            p.AddCoeff(2.9621234697672900e+03, new int[] { 0, 9, 1 });
            p.AddCoeff(-1.1309925975475100e+03, new int[] { 0, 11, 1 });
            p.AddCoeff(-2.6665460146988200e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(5.7775163651807800e+02, new int[] { 2, 3, 1 });
            p.AddCoeff(-3.4665098191084700e+03, new int[] { 2, 5, 1 });
            p.AddCoeff(8.4186667035491400e+03, new int[] { 2, 7, 1 });
            p.AddCoeff(-8.8863704093018700e+03, new int[] { 2, 9, 1 });
            p.AddCoeff(3.3929777926425300e+03, new int[] { 2, 11, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("2e0bbd85-58ee-4c23-8404-7a3584b799c2"));
            OrthonormalPolynomials[601] = p;
            p.AddCoeff(-4.4585335662774400e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(3.4776561816964000e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(-4.3470702271204900e+02, new int[] { 0, 4, 0 });
            p.AddCoeff(1.9706718362946300e+03, new int[] { 0, 6, 0 });
            p.AddCoeff(-4.0117248095997700e+03, new int[] { 0, 8, 0 });
            p.AddCoeff(3.7442764889597800e+03, new int[] { 0, 10, 0 });
            p.AddCoeff(-1.3048236249405400e+03, new int[] { 0, 12, 0 });
            p.AddCoeff(1.3375600698832300e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-1.0432968545089200e+02, new int[] { 2, 2, 0 });
            p.AddCoeff(1.3041210681361500e+03, new int[] { 2, 4, 0 });
            p.AddCoeff(-5.9120155088838900e+03, new int[] { 2, 6, 0 });
            p.AddCoeff(1.2035174428799300e+04, new int[] { 2, 8, 0 });
            p.AddCoeff(-1.1232829466879400e+04, new int[] { 2, 10, 0 });
            p.AddCoeff(3.9144708748216000e+03, new int[] { 2, 12, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d6ef907c-dd67-4892-9516-fee6940054e2"));
            OrthonormalPolynomials[602] = p;
            p.AddCoeff(1.8215977151856400e+01, new int[] { 1, 0, 1 });
            p.AddCoeff(-3.9467950495688900e+02, new int[] { 1, 0, 3 });
            p.AddCoeff(2.3680770297413400e+03, new int[] { 1, 0, 5 });
            p.AddCoeff(-5.7510442150861000e+03, new int[] { 1, 0, 7 });
            p.AddCoeff(6.0705466714797800e+03, new int[] { 1, 0, 9 });
            p.AddCoeff(-2.3178450927468200e+03, new int[] { 1, 0, 11 });
            p.AddCoeff(-3.0359961919760700e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(6.5779917492814900e+02, new int[] { 3, 0, 3 });
            p.AddCoeff(-3.9467950495688900e+03, new int[] { 3, 0, 5 });
            p.AddCoeff(9.5850736918101700e+03, new int[] { 3, 0, 7 });
            p.AddCoeff(-1.0117577785799600e+04, new int[] { 3, 0, 9 });
            p.AddCoeff(3.8630751545780400e+03, new int[] { 3, 0, 11 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("77cbc982-2c71-4490-835b-6e7567f0cc8b"));
            OrthonormalPolynomials[603] = p;
            p.AddCoeff(2.7407293110638800e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-1.5074011210851400e+02, new int[] { 1, 1, 2 });
            p.AddCoeff(1.3064143049404500e+03, new int[] { 1, 1, 4 });
            p.AddCoeff(-3.9192429148213500e+03, new int[] { 1, 1, 6 });
            p.AddCoeff(4.7590806822830700e+03, new int[] { 1, 1, 8 });
            p.AddCoeff(-2.0093896214084100e+03, new int[] { 1, 1, 10 });
            p.AddCoeff(-4.5678821851064700e+00, new int[] { 3, 1, 0 });
            p.AddCoeff(2.5123352018085600e+02, new int[] { 3, 1, 2 });
            p.AddCoeff(-2.1773571749007500e+03, new int[] { 3, 1, 4 });
            p.AddCoeff(6.5320715247022600e+03, new int[] { 3, 1, 6 });
            p.AddCoeff(-7.9318011371384500e+03, new int[] { 3, 1, 8 });
            p.AddCoeff(3.3489827023473500e+03, new int[] { 3, 1, 10 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c15538dd-e999-4c64-a6f6-9f7ce24c7e5a"));
            OrthonormalPolynomials[604] = p;
            p.AddCoeff(1.6827812978247900e+01, new int[] { 1, 0, 1 });
            p.AddCoeff(-2.4680792368097000e+02, new int[] { 1, 0, 3 });
            p.AddCoeff(9.6255090235578100e+02, new int[] { 1, 0, 5 });
            p.AddCoeff(-1.3750727176511200e+03, new int[] { 1, 0, 7 });
            p.AddCoeff(6.4933989444636000e+02, new int[] { 1, 0, 9 });
            p.AddCoeff(-5.0483438934743800e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(7.4042377104290900e+02, new int[] { 1, 2, 3 });
            p.AddCoeff(-2.8876527070673400e+03, new int[] { 1, 2, 5 });
            p.AddCoeff(4.1252181529533500e+03, new int[] { 1, 2, 7 });
            p.AddCoeff(-1.9480196833390800e+03, new int[] { 1, 2, 9 });
            p.AddCoeff(-2.8046354963746500e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(4.1134653946828300e+02, new int[] { 3, 0, 3 });
            p.AddCoeff(-1.6042515039263000e+03, new int[] { 3, 0, 5 });
            p.AddCoeff(2.2917878627518600e+03, new int[] { 3, 0, 7 });
            p.AddCoeff(-1.0822331574106000e+03, new int[] { 3, 0, 9 });
            p.AddCoeff(8.4139064891239600e+01, new int[] { 3, 2, 1 });
            p.AddCoeff(-1.2340396184048500e+03, new int[] { 3, 2, 3 });
            p.AddCoeff(4.8127545117789100e+03, new int[] { 3, 2, 5 });
            p.AddCoeff(-6.8753635882555800e+03, new int[] { 3, 2, 7 });
            p.AddCoeff(3.2466994722318000e+03, new int[] { 3, 2, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("96f0772a-c191-4ad7-99ca-d940e4f99924"));
            OrthonormalPolynomials[605] = p;
            p.AddCoeff(6.2779535781903700e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-2.2600632881485300e+02, new int[] { 1, 1, 2 });
            p.AddCoeff(1.2430348084816900e+03, new int[] { 1, 1, 4 });
            p.AddCoeff(-2.1545936680349400e+03, new int[] { 1, 1, 6 });
            p.AddCoeff(1.1542466078758600e+03, new int[] { 1, 1, 8 });
            p.AddCoeff(-1.0463255963650600e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(3.7667721469142200e+02, new int[] { 1, 3, 2 });
            p.AddCoeff(-2.0717246808028200e+03, new int[] { 1, 3, 4 });
            p.AddCoeff(3.5909894467248900e+03, new int[] { 1, 3, 6 });
            p.AddCoeff(-1.9237443464597600e+03, new int[] { 1, 3, 8 });
            p.AddCoeff(-1.0463255963650600e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(3.7667721469142200e+02, new int[] { 3, 1, 2 });
            p.AddCoeff(-2.0717246808028200e+03, new int[] { 3, 1, 4 });
            p.AddCoeff(3.5909894467248900e+03, new int[] { 3, 1, 6 });
            p.AddCoeff(-1.9237443464597600e+03, new int[] { 3, 1, 8 });
            p.AddCoeff(1.7438759939417700e+01, new int[] { 3, 3, 0 });
            p.AddCoeff(-6.2779535781903700e+02, new int[] { 3, 3, 2 });
            p.AddCoeff(3.4528744680047100e+03, new int[] { 3, 3, 4 });
            p.AddCoeff(-5.9849824112081600e+03, new int[] { 3, 3, 6 });
            p.AddCoeff(3.2062405774329400e+03, new int[] { 3, 3, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f23fca07-8bd7-42f3-9bed-560f1e1a2625"));
            OrthonormalPolynomials[606] = p;
            p.AddCoeff(1.3373389672997100e+01, new int[] { 1, 0, 1 });
            p.AddCoeff(-1.2036050705697300e+02, new int[] { 1, 0, 3 });
            p.AddCoeff(2.6479311552534200e+02, new int[] { 1, 0, 5 });
            p.AddCoeff(-1.6391954770616400e+02, new int[] { 1, 0, 7 });
            p.AddCoeff(-1.3373389672997100e+02, new int[] { 1, 2, 1 });
            p.AddCoeff(1.2036050705697300e+03, new int[] { 1, 2, 3 });
            p.AddCoeff(-2.6479311552534200e+03, new int[] { 1, 2, 5 });
            p.AddCoeff(1.6391954770616400e+03, new int[] { 1, 2, 7 });
            p.AddCoeff(1.5602287951829900e+02, new int[] { 1, 4, 1 });
            p.AddCoeff(-1.4042059156646900e+03, new int[] { 1, 4, 3 });
            p.AddCoeff(3.0892530144623200e+03, new int[] { 1, 4, 5 });
            p.AddCoeff(-1.9123947232385800e+03, new int[] { 1, 4, 7 });
            p.AddCoeff(-2.2288982788328400e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(2.0060084509495600e+02, new int[] { 3, 0, 3 });
            p.AddCoeff(-4.4132185920890300e+02, new int[] { 3, 0, 5 });
            p.AddCoeff(2.7319924617694000e+02, new int[] { 3, 0, 7 });
            p.AddCoeff(2.2288982788328400e+02, new int[] { 3, 2, 1 });
            p.AddCoeff(-2.0060084509495600e+03, new int[] { 3, 2, 3 });
            p.AddCoeff(4.4132185920890300e+03, new int[] { 3, 2, 5 });
            p.AddCoeff(-2.7319924617694000e+03, new int[] { 3, 2, 7 });
            p.AddCoeff(-2.6003813253049800e+02, new int[] { 3, 4, 1 });
            p.AddCoeff(2.3403431927744800e+03, new int[] { 3, 4, 3 });
            p.AddCoeff(-5.1487550241038700e+03, new int[] { 3, 4, 5 });
            p.AddCoeff(3.1873245387309600e+03, new int[] { 3, 4, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c4b87ca7-f3c5-4ee1-baed-aefc6e012412"));
            OrthonormalPolynomials[607] = p;
            p.AddCoeff(9.8313826118541900e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-2.0645903484893800e+02, new int[] { 1, 1, 2 });
            p.AddCoeff(6.1937710454681400e+02, new int[] { 1, 1, 4 });
            p.AddCoeff(-4.5420987666766400e+02, new int[] { 1, 1, 6 });
            p.AddCoeff(-4.5879785521986200e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(9.6347549596171100e+02, new int[] { 1, 3, 2 });
            p.AddCoeff(-2.8904264878851300e+03, new int[] { 1, 3, 4 });
            p.AddCoeff(2.1196460911157600e+03, new int[] { 1, 3, 6 });
            p.AddCoeff(4.1291806969787600e+01, new int[] { 1, 5, 0 });
            p.AddCoeff(-8.6712794636554000e+02, new int[] { 1, 5, 2 });
            p.AddCoeff(2.6013838390966200e+03, new int[] { 1, 5, 4 });
            p.AddCoeff(-1.9076814820041900e+03, new int[] { 1, 5, 6 });
            p.AddCoeff(-1.6385637686423700e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(3.4409839141489700e+02, new int[] { 3, 1, 2 });
            p.AddCoeff(-1.0322951742446900e+03, new int[] { 3, 1, 4 });
            p.AddCoeff(7.5701646111277300e+02, new int[] { 3, 1, 6 });
            p.AddCoeff(7.6466309203310400e+01, new int[] { 3, 3, 0 });
            p.AddCoeff(-1.6057924932695200e+03, new int[] { 3, 3, 2 });
            p.AddCoeff(4.8173774798085600e+03, new int[] { 3, 3, 4 });
            p.AddCoeff(-3.5327434851929400e+03, new int[] { 3, 3, 6 });
            p.AddCoeff(-6.8819678282979400e+01, new int[] { 3, 5, 0 });
            p.AddCoeff(1.4452132439425700e+03, new int[] { 3, 5, 2 });
            p.AddCoeff(-4.3356397318277000e+03, new int[] { 3, 5, 4 });
            p.AddCoeff(3.1794691366736500e+03, new int[] { 3, 5, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d893d8dc-1bd8-40a5-b570-09213be90f48"));
            OrthonormalPolynomials[608] = p;
            p.AddCoeff(9.8313826118541900e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-4.5879785521986200e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(4.1291806969787600e+01, new int[] { 1, 0, 5 });
            p.AddCoeff(-2.0645903484893800e+02, new int[] { 1, 2, 1 });
            p.AddCoeff(9.6347549596171100e+02, new int[] { 1, 2, 3 });
            p.AddCoeff(-8.6712794636554000e+02, new int[] { 1, 2, 5 });
            p.AddCoeff(6.1937710454681400e+02, new int[] { 1, 4, 1 });
            p.AddCoeff(-2.8904264878851300e+03, new int[] { 1, 4, 3 });
            p.AddCoeff(2.6013838390966200e+03, new int[] { 1, 4, 5 });
            p.AddCoeff(-4.5420987666766400e+02, new int[] { 1, 6, 1 });
            p.AddCoeff(2.1196460911157600e+03, new int[] { 1, 6, 3 });
            p.AddCoeff(-1.9076814820041900e+03, new int[] { 1, 6, 5 });
            p.AddCoeff(-1.6385637686423700e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(7.6466309203310400e+01, new int[] { 3, 0, 3 });
            p.AddCoeff(-6.8819678282979400e+01, new int[] { 3, 0, 5 });
            p.AddCoeff(3.4409839141489700e+02, new int[] { 3, 2, 1 });
            p.AddCoeff(-1.6057924932695200e+03, new int[] { 3, 2, 3 });
            p.AddCoeff(1.4452132439425700e+03, new int[] { 3, 2, 5 });
            p.AddCoeff(-1.0322951742446900e+03, new int[] { 3, 4, 1 });
            p.AddCoeff(4.8173774798085600e+03, new int[] { 3, 4, 3 });
            p.AddCoeff(-4.3356397318277000e+03, new int[] { 3, 4, 5 });
            p.AddCoeff(7.5701646111277300e+02, new int[] { 3, 6, 1 });
            p.AddCoeff(-3.5327434851929400e+03, new int[] { 3, 6, 3 });
            p.AddCoeff(3.1794691366736500e+03, new int[] { 3, 6, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5c1107dd-f38a-4a72-9137-1cd87a9d15ac"));
            OrthonormalPolynomials[609] = p;
            p.AddCoeff(1.3373389672997100e+01, new int[] { 1, 1, 0 });
            p.AddCoeff(-1.3373389672997100e+02, new int[] { 1, 1, 2 });
            p.AddCoeff(1.5602287951829900e+02, new int[] { 1, 1, 4 });
            p.AddCoeff(-1.2036050705697300e+02, new int[] { 1, 3, 0 });
            p.AddCoeff(1.2036050705697300e+03, new int[] { 1, 3, 2 });
            p.AddCoeff(-1.4042059156646900e+03, new int[] { 1, 3, 4 });
            p.AddCoeff(2.6479311552534200e+02, new int[] { 1, 5, 0 });
            p.AddCoeff(-2.6479311552534200e+03, new int[] { 1, 5, 2 });
            p.AddCoeff(3.0892530144623200e+03, new int[] { 1, 5, 4 });
            p.AddCoeff(-1.6391954770616400e+02, new int[] { 1, 7, 0 });
            p.AddCoeff(1.6391954770616400e+03, new int[] { 1, 7, 2 });
            p.AddCoeff(-1.9123947232385800e+03, new int[] { 1, 7, 4 });
            p.AddCoeff(-2.2288982788328400e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(2.2288982788328400e+02, new int[] { 3, 1, 2 });
            p.AddCoeff(-2.6003813253049800e+02, new int[] { 3, 1, 4 });
            p.AddCoeff(2.0060084509495600e+02, new int[] { 3, 3, 0 });
            p.AddCoeff(-2.0060084509495600e+03, new int[] { 3, 3, 2 });
            p.AddCoeff(2.3403431927744800e+03, new int[] { 3, 3, 4 });
            p.AddCoeff(-4.4132185920890300e+02, new int[] { 3, 5, 0 });
            p.AddCoeff(4.4132185920890300e+03, new int[] { 3, 5, 2 });
            p.AddCoeff(-5.1487550241038700e+03, new int[] { 3, 5, 4 });
            p.AddCoeff(2.7319924617694000e+02, new int[] { 3, 7, 0 });
            p.AddCoeff(-2.7319924617694000e+03, new int[] { 3, 7, 2 });
            p.AddCoeff(3.1873245387309600e+03, new int[] { 3, 7, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("152e18ee-faeb-412c-837b-4853eccaad82"));
            OrthonormalPolynomials[610] = p;
            p.AddCoeff(6.2779535781903700e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-1.0463255963650600e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(-2.2600632881485300e+02, new int[] { 1, 2, 1 });
            p.AddCoeff(3.7667721469142200e+02, new int[] { 1, 2, 3 });
            p.AddCoeff(1.2430348084816900e+03, new int[] { 1, 4, 1 });
            p.AddCoeff(-2.0717246808028200e+03, new int[] { 1, 4, 3 });
            p.AddCoeff(-2.1545936680349400e+03, new int[] { 1, 6, 1 });
            p.AddCoeff(3.5909894467248900e+03, new int[] { 1, 6, 3 });
            p.AddCoeff(1.1542466078758600e+03, new int[] { 1, 8, 1 });
            p.AddCoeff(-1.9237443464597600e+03, new int[] { 1, 8, 3 });
            p.AddCoeff(-1.0463255963650600e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(1.7438759939417700e+01, new int[] { 3, 0, 3 });
            p.AddCoeff(3.7667721469142200e+02, new int[] { 3, 2, 1 });
            p.AddCoeff(-6.2779535781903700e+02, new int[] { 3, 2, 3 });
            p.AddCoeff(-2.0717246808028200e+03, new int[] { 3, 4, 1 });
            p.AddCoeff(3.4528744680047100e+03, new int[] { 3, 4, 3 });
            p.AddCoeff(3.5909894467248900e+03, new int[] { 3, 6, 1 });
            p.AddCoeff(-5.9849824112081600e+03, new int[] { 3, 6, 3 });
            p.AddCoeff(-1.9237443464597600e+03, new int[] { 3, 8, 1 });
            p.AddCoeff(3.2062405774329400e+03, new int[] { 3, 8, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("643a1f60-5e0a-4b47-b8c9-7e8fa58ad341"));
            OrthonormalPolynomials[611] = p;
            p.AddCoeff(1.6827812978247900e+01, new int[] { 1, 1, 0 });
            p.AddCoeff(-5.0483438934743800e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(-2.4680792368097000e+02, new int[] { 1, 3, 0 });
            p.AddCoeff(7.4042377104290900e+02, new int[] { 1, 3, 2 });
            p.AddCoeff(9.6255090235578100e+02, new int[] { 1, 5, 0 });
            p.AddCoeff(-2.8876527070673400e+03, new int[] { 1, 5, 2 });
            p.AddCoeff(-1.3750727176511200e+03, new int[] { 1, 7, 0 });
            p.AddCoeff(4.1252181529533500e+03, new int[] { 1, 7, 2 });
            p.AddCoeff(6.4933989444636000e+02, new int[] { 1, 9, 0 });
            p.AddCoeff(-1.9480196833390800e+03, new int[] { 1, 9, 2 });
            p.AddCoeff(-2.8046354963746500e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(8.4139064891239600e+01, new int[] { 3, 1, 2 });
            p.AddCoeff(4.1134653946828300e+02, new int[] { 3, 3, 0 });
            p.AddCoeff(-1.2340396184048500e+03, new int[] { 3, 3, 2 });
            p.AddCoeff(-1.6042515039263000e+03, new int[] { 3, 5, 0 });
            p.AddCoeff(4.8127545117789100e+03, new int[] { 3, 5, 2 });
            p.AddCoeff(2.2917878627518600e+03, new int[] { 3, 7, 0 });
            p.AddCoeff(-6.8753635882555800e+03, new int[] { 3, 7, 2 });
            p.AddCoeff(-1.0822331574106000e+03, new int[] { 3, 9, 0 });
            p.AddCoeff(3.2466994722318000e+03, new int[] { 3, 9, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f5e32c10-bba2-4c6f-8dbc-f59a492d6a4c"));
            OrthonormalPolynomials[612] = p;
            p.AddCoeff(2.7407293110638800e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-1.5074011210851400e+02, new int[] { 1, 2, 1 });
            p.AddCoeff(1.3064143049404500e+03, new int[] { 1, 4, 1 });
            p.AddCoeff(-3.9192429148213500e+03, new int[] { 1, 6, 1 });
            p.AddCoeff(4.7590806822830700e+03, new int[] { 1, 8, 1 });
            p.AddCoeff(-2.0093896214084100e+03, new int[] { 1, 10, 1 });
            p.AddCoeff(-4.5678821851064700e+00, new int[] { 3, 0, 1 });
            p.AddCoeff(2.5123352018085600e+02, new int[] { 3, 2, 1 });
            p.AddCoeff(-2.1773571749007500e+03, new int[] { 3, 4, 1 });
            p.AddCoeff(6.5320715247022600e+03, new int[] { 3, 6, 1 });
            p.AddCoeff(-7.9318011371384500e+03, new int[] { 3, 8, 1 });
            p.AddCoeff(3.3489827023473500e+03, new int[] { 3, 10, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("580d4498-9721-4517-aa68-ca3e68a4f04c"));
            OrthonormalPolynomials[613] = p;
            p.AddCoeff(1.8215977151856400e+01, new int[] { 1, 1, 0 });
            p.AddCoeff(-3.9467950495688900e+02, new int[] { 1, 3, 0 });
            p.AddCoeff(2.3680770297413400e+03, new int[] { 1, 5, 0 });
            p.AddCoeff(-5.7510442150861000e+03, new int[] { 1, 7, 0 });
            p.AddCoeff(6.0705466714797800e+03, new int[] { 1, 9, 0 });
            p.AddCoeff(-2.3178450927468200e+03, new int[] { 1, 11, 0 });
            p.AddCoeff(-3.0359961919760700e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(6.5779917492814900e+02, new int[] { 3, 3, 0 });
            p.AddCoeff(-3.9467950495688900e+03, new int[] { 3, 5, 0 });
            p.AddCoeff(9.5850736918101700e+03, new int[] { 3, 7, 0 });
            p.AddCoeff(-1.0117577785799600e+04, new int[] { 3, 9, 0 });
            p.AddCoeff(3.8630751545780400e+03, new int[] { 3, 11, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("9113f0ca-6e96-4f53-9e86-e01349d45fbc"));
            OrthonormalPolynomials[614] = p;
            p.AddCoeff(-4.4855712597622800e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(2.4670641928692500e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(-2.1381223004866800e+02, new int[] { 0, 0, 4 });
            p.AddCoeff(6.4143669014600500e+02, new int[] { 0, 0, 6 });
            p.AddCoeff(-7.7888740946301000e+02, new int[] { 0, 0, 8 });
            p.AddCoeff(3.2886357288438100e+02, new int[] { 0, 0, 10 });
            p.AddCoeff(4.4855712597622800e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-2.4670641928692500e+02, new int[] { 2, 0, 2 });
            p.AddCoeff(2.1381223004866800e+03, new int[] { 2, 0, 4 });
            p.AddCoeff(-6.4143669014600500e+03, new int[] { 2, 0, 6 });
            p.AddCoeff(7.7888740946301000e+03, new int[] { 2, 0, 8 });
            p.AddCoeff(-3.2886357288438100e+03, new int[] { 2, 0, 10 });
            p.AddCoeff(-5.2331664697226600e+00, new int[] { 4, 0, 0 });
            p.AddCoeff(2.8782415583474600e+02, new int[] { 4, 0, 2 });
            p.AddCoeff(-2.4944760172344700e+03, new int[] { 4, 0, 4 });
            p.AddCoeff(7.4834280517033600e+03, new int[] { 4, 0, 6 });
            p.AddCoeff(-9.0870197770684300e+03, new int[] { 4, 0, 8 });
            p.AddCoeff(3.8367416836511000e+03, new int[] { 4, 0, 10 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("58c70202-8636-4bc0-8911-71ff83decf43"));
            OrthonormalPolynomials[615] = p;
            p.AddCoeff(7.3900187608663900e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-1.0838694182604000e+02, new int[] { 0, 1, 3 });
            p.AddCoeff(4.2270907312155800e+02, new int[] { 0, 1, 5 });
            p.AddCoeff(-6.0387010445936800e+02, new int[] { 0, 1, 7 });
            p.AddCoeff(2.8516088266136800e+02, new int[] { 0, 1, 9 });
            p.AddCoeff(-7.3900187608663900e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(1.0838694182604000e+03, new int[] { 2, 1, 3 });
            p.AddCoeff(-4.2270907312155800e+03, new int[] { 2, 1, 5 });
            p.AddCoeff(6.0387010445936800e+03, new int[] { 2, 1, 7 });
            p.AddCoeff(-2.8516088266136800e+03, new int[] { 2, 1, 9 });
            p.AddCoeff(8.6216885543441200e+01, new int[] { 4, 1, 1 });
            p.AddCoeff(-1.2645143213038000e+03, new int[] { 4, 1, 3 });
            p.AddCoeff(4.9316058530848400e+03, new int[] { 4, 1, 5 });
            p.AddCoeff(-7.0451512186926300e+03, new int[] { 4, 1, 7 });
            p.AddCoeff(3.3268769643826300e+03, new int[] { 4, 1, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f356eb11-a96e-4939-af41-812f9ab2e5ef"));
            OrthonormalPolynomials[616] = p;
            p.AddCoeff(-5.0135467715791900e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.8048768377685100e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(-9.9268226077267900e+01, new int[] { 0, 0, 4 });
            p.AddCoeff(1.7206492520059800e+02, new int[] { 0, 0, 6 });
            p.AddCoeff(-9.2177638500320200e+01, new int[] { 0, 0, 8 });
            p.AddCoeff(1.5040640314737600e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-5.4146305133055200e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(2.9780467823180400e+02, new int[] { 0, 2, 4 });
            p.AddCoeff(-5.1619477560179300e+02, new int[] { 0, 2, 6 });
            p.AddCoeff(2.7653291550096100e+02, new int[] { 0, 2, 8 });
            p.AddCoeff(5.0135467715791900e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-1.8048768377685100e+02, new int[] { 2, 0, 2 });
            p.AddCoeff(9.9268226077267900e+02, new int[] { 2, 0, 4 });
            p.AddCoeff(-1.7206492520059800e+03, new int[] { 2, 0, 6 });
            p.AddCoeff(9.2177638500320200e+02, new int[] { 2, 0, 8 });
            p.AddCoeff(-1.5040640314737600e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(5.4146305133055200e+02, new int[] { 2, 2, 2 });
            p.AddCoeff(-2.9780467823180400e+03, new int[] { 2, 2, 4 });
            p.AddCoeff(5.1619477560179300e+03, new int[] { 2, 2, 6 });
            p.AddCoeff(-2.7653291550096100e+03, new int[] { 2, 2, 8 });
            p.AddCoeff(-5.8491379001757200e+00, new int[] { 4, 0, 0 });
            p.AddCoeff(2.1056896440632600e+02, new int[] { 4, 0, 2 });
            p.AddCoeff(-1.1581293042347900e+03, new int[] { 4, 0, 4 });
            p.AddCoeff(2.0074241273403100e+03, new int[] { 4, 0, 6 });
            p.AddCoeff(-1.0754057825037400e+03, new int[] { 4, 0, 8 });
            p.AddCoeff(1.7547413700527200e+01, new int[] { 4, 2, 0 });
            p.AddCoeff(-6.3170689321897700e+02, new int[] { 4, 2, 2 });
            p.AddCoeff(3.4743879127043800e+03, new int[] { 4, 2, 4 });
            p.AddCoeff(-6.0222723820209200e+03, new int[] { 4, 2, 6 });
            p.AddCoeff(3.2262173475112000e+03, new int[] { 4, 2, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("75649b34-5022-41b7-aa08-1a1460ba2e54"));
            OrthonormalPolynomials[617] = p;
            p.AddCoeff(1.3373389672997100e+01, new int[] { 0, 1, 1 });
            p.AddCoeff(-1.2036050705697300e+02, new int[] { 0, 1, 3 });
            p.AddCoeff(2.6479311552534200e+02, new int[] { 0, 1, 5 });
            p.AddCoeff(-1.6391954770616400e+02, new int[] { 0, 1, 7 });
            p.AddCoeff(-2.2288982788328400e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(2.0060084509495600e+02, new int[] { 0, 3, 3 });
            p.AddCoeff(-4.4132185920890300e+02, new int[] { 0, 3, 5 });
            p.AddCoeff(2.7319924617694000e+02, new int[] { 0, 3, 7 });
            p.AddCoeff(-1.3373389672997100e+02, new int[] { 2, 1, 1 });
            p.AddCoeff(1.2036050705697300e+03, new int[] { 2, 1, 3 });
            p.AddCoeff(-2.6479311552534200e+03, new int[] { 2, 1, 5 });
            p.AddCoeff(1.6391954770616400e+03, new int[] { 2, 1, 7 });
            p.AddCoeff(2.2288982788328400e+02, new int[] { 2, 3, 1 });
            p.AddCoeff(-2.0060084509495600e+03, new int[] { 2, 3, 3 });
            p.AddCoeff(4.4132185920890300e+03, new int[] { 2, 3, 5 });
            p.AddCoeff(-2.7319924617694000e+03, new int[] { 2, 3, 7 });
            p.AddCoeff(1.5602287951829900e+02, new int[] { 4, 1, 1 });
            p.AddCoeff(-1.4042059156646900e+03, new int[] { 4, 1, 3 });
            p.AddCoeff(3.0892530144623200e+03, new int[] { 4, 1, 5 });
            p.AddCoeff(-1.9123947232385800e+03, new int[] { 4, 1, 7 });
            p.AddCoeff(-2.6003813253049800e+02, new int[] { 4, 3, 1 });
            p.AddCoeff(2.3403431927744800e+03, new int[] { 4, 3, 3 });
            p.AddCoeff(-5.1487550241038700e+03, new int[] { 4, 3, 5 });
            p.AddCoeff(3.1873245387309600e+03, new int[] { 4, 3, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("ce458a26-e16b-4acd-ba7b-f90712cde766"));
            OrthonormalPolynomials[618] = p;
            p.AddCoeff(-5.0417551342897400e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.0587685782008400e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(-3.1763057346025300e+01, new int[] { 0, 0, 4 });
            p.AddCoeff(2.3292908720418600e+01, new int[] { 0, 0, 6 });
            p.AddCoeff(5.0417551342897400e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-1.0587685782008400e+02, new int[] { 0, 2, 2 });
            p.AddCoeff(3.1763057346025300e+02, new int[] { 0, 2, 4 });
            p.AddCoeff(-2.3292908720418600e+02, new int[] { 0, 2, 6 });
            p.AddCoeff(-5.8820476566713600e+00, new int[] { 0, 4, 0 });
            p.AddCoeff(1.2352300079009900e+02, new int[] { 0, 4, 2 });
            p.AddCoeff(-3.7056900237029600e+02, new int[] { 0, 4, 4 });
            p.AddCoeff(2.7175060173821700e+02, new int[] { 0, 4, 6 });
            p.AddCoeff(5.0417551342897400e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-1.0587685782008400e+02, new int[] { 2, 0, 2 });
            p.AddCoeff(3.1763057346025300e+02, new int[] { 2, 0, 4 });
            p.AddCoeff(-2.3292908720418600e+02, new int[] { 2, 0, 6 });
            p.AddCoeff(-5.0417551342897400e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(1.0587685782008400e+03, new int[] { 2, 2, 2 });
            p.AddCoeff(-3.1763057346025300e+03, new int[] { 2, 2, 4 });
            p.AddCoeff(2.3292908720418600e+03, new int[] { 2, 2, 6 });
            p.AddCoeff(5.8820476566713600e+01, new int[] { 2, 4, 0 });
            p.AddCoeff(-1.2352300079009900e+03, new int[] { 2, 4, 2 });
            p.AddCoeff(3.7056900237029600e+03, new int[] { 2, 4, 4 });
            p.AddCoeff(-2.7175060173821700e+03, new int[] { 2, 4, 6 });
            p.AddCoeff(-5.8820476566713600e+00, new int[] { 4, 0, 0 });
            p.AddCoeff(1.2352300079009900e+02, new int[] { 4, 0, 2 });
            p.AddCoeff(-3.7056900237029600e+02, new int[] { 4, 0, 4 });
            p.AddCoeff(2.7175060173821700e+02, new int[] { 4, 0, 6 });
            p.AddCoeff(5.8820476566713600e+01, new int[] { 4, 2, 0 });
            p.AddCoeff(-1.2352300079009900e+03, new int[] { 4, 2, 2 });
            p.AddCoeff(3.7056900237029600e+03, new int[] { 4, 2, 4 });
            p.AddCoeff(-2.7175060173821700e+03, new int[] { 4, 2, 6 });
            p.AddCoeff(-6.8623889327832500e+01, new int[] { 4, 4, 0 });
            p.AddCoeff(1.4411016758844800e+03, new int[] { 4, 4, 2 });
            p.AddCoeff(-4.3233050276534500e+03, new int[] { 4, 4, 4 });
            p.AddCoeff(3.1704236869458600e+03, new int[] { 4, 4, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f946fc7a-387f-44b0-a4f7-1889a8a58082"));
            OrthonormalPolynomials[619] = p;
            p.AddCoeff(1.5381644092705500e+01, new int[] { 0, 1, 1 });
            p.AddCoeff(-7.1781005765958900e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(6.4602905189363000e+01, new int[] { 0, 1, 5 });
            p.AddCoeff(-7.1781005765958900e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(3.3497802690780800e+02, new int[] { 0, 3, 3 });
            p.AddCoeff(-3.0148022421702700e+02, new int[] { 0, 3, 5 });
            p.AddCoeff(6.4602905189363000e+01, new int[] { 0, 5, 1 });
            p.AddCoeff(-3.0148022421702700e+02, new int[] { 0, 5, 3 });
            p.AddCoeff(2.7133220179532400e+02, new int[] { 0, 5, 5 });
            p.AddCoeff(-1.5381644092705500e+02, new int[] { 2, 1, 1 });
            p.AddCoeff(7.1781005765958900e+02, new int[] { 2, 1, 3 });
            p.AddCoeff(-6.4602905189363000e+02, new int[] { 2, 1, 5 });
            p.AddCoeff(7.1781005765958900e+02, new int[] { 2, 3, 1 });
            p.AddCoeff(-3.3497802690780800e+03, new int[] { 2, 3, 3 });
            p.AddCoeff(3.0148022421702700e+03, new int[] { 2, 3, 5 });
            p.AddCoeff(-6.4602905189363000e+02, new int[] { 2, 5, 1 });
            p.AddCoeff(3.0148022421702700e+03, new int[] { 2, 5, 3 });
            p.AddCoeff(-2.7133220179532400e+03, new int[] { 2, 5, 5 });
            p.AddCoeff(1.7945251441489700e+02, new int[] { 4, 1, 1 });
            p.AddCoeff(-8.3744506726952000e+02, new int[] { 4, 1, 3 });
            p.AddCoeff(7.5370056054256800e+02, new int[] { 4, 1, 5 });
            p.AddCoeff(-8.3744506726952000e+02, new int[] { 4, 3, 1 });
            p.AddCoeff(3.9080769805910900e+03, new int[] { 4, 3, 3 });
            p.AddCoeff(-3.5172692825319800e+03, new int[] { 4, 3, 5 });
            p.AddCoeff(7.5370056054256800e+02, new int[] { 4, 5, 1 });
            p.AddCoeff(-3.5172692825319800e+03, new int[] { 4, 5, 3 });
            p.AddCoeff(3.1655423542787900e+03, new int[] { 4, 5, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d1fa0389-019c-4a45-80fb-5437dc81eb40"));
            OrthonormalPolynomials[620] = p;
            p.AddCoeff(-5.0417551342897400e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(5.0417551342897400e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-5.8820476566713600e+00, new int[] { 0, 0, 4 });
            p.AddCoeff(1.0587685782008400e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(-1.0587685782008400e+02, new int[] { 0, 2, 2 });
            p.AddCoeff(1.2352300079009900e+02, new int[] { 0, 2, 4 });
            p.AddCoeff(-3.1763057346025300e+01, new int[] { 0, 4, 0 });
            p.AddCoeff(3.1763057346025300e+02, new int[] { 0, 4, 2 });
            p.AddCoeff(-3.7056900237029600e+02, new int[] { 0, 4, 4 });
            p.AddCoeff(2.3292908720418600e+01, new int[] { 0, 6, 0 });
            p.AddCoeff(-2.3292908720418600e+02, new int[] { 0, 6, 2 });
            p.AddCoeff(2.7175060173821700e+02, new int[] { 0, 6, 4 });
            p.AddCoeff(5.0417551342897400e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-5.0417551342897400e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(5.8820476566713600e+01, new int[] { 2, 0, 4 });
            p.AddCoeff(-1.0587685782008400e+02, new int[] { 2, 2, 0 });
            p.AddCoeff(1.0587685782008400e+03, new int[] { 2, 2, 2 });
            p.AddCoeff(-1.2352300079009900e+03, new int[] { 2, 2, 4 });
            p.AddCoeff(3.1763057346025300e+02, new int[] { 2, 4, 0 });
            p.AddCoeff(-3.1763057346025300e+03, new int[] { 2, 4, 2 });
            p.AddCoeff(3.7056900237029600e+03, new int[] { 2, 4, 4 });
            p.AddCoeff(-2.3292908720418600e+02, new int[] { 2, 6, 0 });
            p.AddCoeff(2.3292908720418600e+03, new int[] { 2, 6, 2 });
            p.AddCoeff(-2.7175060173821700e+03, new int[] { 2, 6, 4 });
            p.AddCoeff(-5.8820476566713600e+00, new int[] { 4, 0, 0 });
            p.AddCoeff(5.8820476566713600e+01, new int[] { 4, 0, 2 });
            p.AddCoeff(-6.8623889327832500e+01, new int[] { 4, 0, 4 });
            p.AddCoeff(1.2352300079009900e+02, new int[] { 4, 2, 0 });
            p.AddCoeff(-1.2352300079009900e+03, new int[] { 4, 2, 2 });
            p.AddCoeff(1.4411016758844800e+03, new int[] { 4, 2, 4 });
            p.AddCoeff(-3.7056900237029600e+02, new int[] { 4, 4, 0 });
            p.AddCoeff(3.7056900237029600e+03, new int[] { 4, 4, 2 });
            p.AddCoeff(-4.3233050276534500e+03, new int[] { 4, 4, 4 });
            p.AddCoeff(2.7175060173821700e+02, new int[] { 4, 6, 0 });
            p.AddCoeff(-2.7175060173821700e+03, new int[] { 4, 6, 2 });
            p.AddCoeff(3.1704236869458600e+03, new int[] { 4, 6, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("9377096d-ce7a-478e-b2a7-c928119c0791"));
            OrthonormalPolynomials[621] = p;
            p.AddCoeff(1.3373389672997100e+01, new int[] { 0, 1, 1 });
            p.AddCoeff(-2.2288982788328400e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(-1.2036050705697300e+02, new int[] { 0, 3, 1 });
            p.AddCoeff(2.0060084509495600e+02, new int[] { 0, 3, 3 });
            p.AddCoeff(2.6479311552534200e+02, new int[] { 0, 5, 1 });
            p.AddCoeff(-4.4132185920890300e+02, new int[] { 0, 5, 3 });
            p.AddCoeff(-1.6391954770616400e+02, new int[] { 0, 7, 1 });
            p.AddCoeff(2.7319924617694000e+02, new int[] { 0, 7, 3 });
            p.AddCoeff(-1.3373389672997100e+02, new int[] { 2, 1, 1 });
            p.AddCoeff(2.2288982788328400e+02, new int[] { 2, 1, 3 });
            p.AddCoeff(1.2036050705697300e+03, new int[] { 2, 3, 1 });
            p.AddCoeff(-2.0060084509495600e+03, new int[] { 2, 3, 3 });
            p.AddCoeff(-2.6479311552534200e+03, new int[] { 2, 5, 1 });
            p.AddCoeff(4.4132185920890300e+03, new int[] { 2, 5, 3 });
            p.AddCoeff(1.6391954770616400e+03, new int[] { 2, 7, 1 });
            p.AddCoeff(-2.7319924617694000e+03, new int[] { 2, 7, 3 });
            p.AddCoeff(1.5602287951829900e+02, new int[] { 4, 1, 1 });
            p.AddCoeff(-2.6003813253049800e+02, new int[] { 4, 1, 3 });
            p.AddCoeff(-1.4042059156646900e+03, new int[] { 4, 3, 1 });
            p.AddCoeff(2.3403431927744800e+03, new int[] { 4, 3, 3 });
            p.AddCoeff(3.0892530144623200e+03, new int[] { 4, 5, 1 });
            p.AddCoeff(-5.1487550241038700e+03, new int[] { 4, 5, 3 });
            p.AddCoeff(-1.9123947232385800e+03, new int[] { 4, 7, 1 });
            p.AddCoeff(3.1873245387309600e+03, new int[] { 4, 7, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e1e86a3f-e9b8-447d-ab80-de99a9cc0a4c"));
            OrthonormalPolynomials[622] = p;
            p.AddCoeff(-5.0135467715791900e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.5040640314737600e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(1.8048768377685100e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(-5.4146305133055200e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(-9.9268226077267900e+01, new int[] { 0, 4, 0 });
            p.AddCoeff(2.9780467823180400e+02, new int[] { 0, 4, 2 });
            p.AddCoeff(1.7206492520059800e+02, new int[] { 0, 6, 0 });
            p.AddCoeff(-5.1619477560179300e+02, new int[] { 0, 6, 2 });
            p.AddCoeff(-9.2177638500320200e+01, new int[] { 0, 8, 0 });
            p.AddCoeff(2.7653291550096100e+02, new int[] { 0, 8, 2 });
            p.AddCoeff(5.0135467715791900e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-1.5040640314737600e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(-1.8048768377685100e+02, new int[] { 2, 2, 0 });
            p.AddCoeff(5.4146305133055200e+02, new int[] { 2, 2, 2 });
            p.AddCoeff(9.9268226077267900e+02, new int[] { 2, 4, 0 });
            p.AddCoeff(-2.9780467823180400e+03, new int[] { 2, 4, 2 });
            p.AddCoeff(-1.7206492520059800e+03, new int[] { 2, 6, 0 });
            p.AddCoeff(5.1619477560179300e+03, new int[] { 2, 6, 2 });
            p.AddCoeff(9.2177638500320200e+02, new int[] { 2, 8, 0 });
            p.AddCoeff(-2.7653291550096100e+03, new int[] { 2, 8, 2 });
            p.AddCoeff(-5.8491379001757200e+00, new int[] { 4, 0, 0 });
            p.AddCoeff(1.7547413700527200e+01, new int[] { 4, 0, 2 });
            p.AddCoeff(2.1056896440632600e+02, new int[] { 4, 2, 0 });
            p.AddCoeff(-6.3170689321897700e+02, new int[] { 4, 2, 2 });
            p.AddCoeff(-1.1581293042347900e+03, new int[] { 4, 4, 0 });
            p.AddCoeff(3.4743879127043800e+03, new int[] { 4, 4, 2 });
            p.AddCoeff(2.0074241273403100e+03, new int[] { 4, 6, 0 });
            p.AddCoeff(-6.0222723820209200e+03, new int[] { 4, 6, 2 });
            p.AddCoeff(-1.0754057825037400e+03, new int[] { 4, 8, 0 });
            p.AddCoeff(3.2262173475112000e+03, new int[] { 4, 8, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("59f3a503-31e3-49a8-ba63-f0d7cfb0cbe8"));
            OrthonormalPolynomials[623] = p;
            p.AddCoeff(7.3900187608663900e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-1.0838694182604000e+02, new int[] { 0, 3, 1 });
            p.AddCoeff(4.2270907312155800e+02, new int[] { 0, 5, 1 });
            p.AddCoeff(-6.0387010445936800e+02, new int[] { 0, 7, 1 });
            p.AddCoeff(2.8516088266136800e+02, new int[] { 0, 9, 1 });
            p.AddCoeff(-7.3900187608663900e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(1.0838694182604000e+03, new int[] { 2, 3, 1 });
            p.AddCoeff(-4.2270907312155800e+03, new int[] { 2, 5, 1 });
            p.AddCoeff(6.0387010445936800e+03, new int[] { 2, 7, 1 });
            p.AddCoeff(-2.8516088266136800e+03, new int[] { 2, 9, 1 });
            p.AddCoeff(8.6216885543441200e+01, new int[] { 4, 1, 1 });
            p.AddCoeff(-1.2645143213038000e+03, new int[] { 4, 3, 1 });
            p.AddCoeff(4.9316058530848400e+03, new int[] { 4, 5, 1 });
            p.AddCoeff(-7.0451512186926300e+03, new int[] { 4, 7, 1 });
            p.AddCoeff(3.3268769643826300e+03, new int[] { 4, 9, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("81e8a9e0-29ed-48a7-af5e-68c5662ef7e6"));
            OrthonormalPolynomials[624] = p;
            p.AddCoeff(-4.4855712597622800e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(2.4670641928692500e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(-2.1381223004866800e+02, new int[] { 0, 4, 0 });
            p.AddCoeff(6.4143669014600500e+02, new int[] { 0, 6, 0 });
            p.AddCoeff(-7.7888740946301000e+02, new int[] { 0, 8, 0 });
            p.AddCoeff(3.2886357288438100e+02, new int[] { 0, 10, 0 });
            p.AddCoeff(4.4855712597622800e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-2.4670641928692500e+02, new int[] { 2, 2, 0 });
            p.AddCoeff(2.1381223004866800e+03, new int[] { 2, 4, 0 });
            p.AddCoeff(-6.4143669014600500e+03, new int[] { 2, 6, 0 });
            p.AddCoeff(7.7888740946301000e+03, new int[] { 2, 8, 0 });
            p.AddCoeff(-3.2886357288438100e+03, new int[] { 2, 10, 0 });
            p.AddCoeff(-5.2331664697226600e+00, new int[] { 4, 0, 0 });
            p.AddCoeff(2.8782415583474600e+02, new int[] { 4, 2, 0 });
            p.AddCoeff(-2.4944760172344700e+03, new int[] { 4, 4, 0 });
            p.AddCoeff(7.4834280517033600e+03, new int[] { 4, 6, 0 });
            p.AddCoeff(-9.0870197770684300e+03, new int[] { 4, 8, 0 });
            p.AddCoeff(3.8367416836511000e+03, new int[] { 4, 10, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("dfd016af-42f7-4f4a-82e1-3c8a8dae851e"));
            OrthonormalPolynomials[625] = p;
            p.AddCoeff(2.3584680961604700e+01, new int[] { 1, 0, 1 });
            p.AddCoeff(-3.4590865410353500e+02, new int[] { 1, 0, 3 });
            p.AddCoeff(1.3490437510037900e+03, new int[] { 1, 0, 5 });
            p.AddCoeff(-1.9272053585768400e+03, new int[] { 1, 0, 7 });
            p.AddCoeff(9.1006919710573000e+02, new int[] { 1, 0, 9 });
            p.AddCoeff(-1.1006184448748800e+02, new int[] { 3, 0, 1 });
            p.AddCoeff(1.6142403858165000e+03, new int[] { 3, 0, 3 });
            p.AddCoeff(-6.2955375046843400e+03, new int[] { 3, 0, 5 });
            p.AddCoeff(8.9936250066919200e+03, new int[] { 3, 0, 7 });
            p.AddCoeff(-4.2469895864934000e+03, new int[] { 3, 0, 9 });
            p.AddCoeff(9.9055660038739600e+01, new int[] { 5, 0, 1 });
            p.AddCoeff(-1.4528163472348500e+03, new int[] { 5, 0, 3 });
            p.AddCoeff(5.6659837542159100e+03, new int[] { 5, 0, 5 });
            p.AddCoeff(-8.0942625060227200e+03, new int[] { 5, 0, 7 });
            p.AddCoeff(3.8222906278440600e+03, new int[] { 5, 0, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("6fc5fb70-db0a-4b3e-ba2d-884bb7677650"));
            OrthonormalPolynomials[626] = p;
            p.AddCoeff(4.2933449549966900e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-1.5456041837988100e+02, new int[] { 1, 1, 2 });
            p.AddCoeff(8.5008230108934600e+02, new int[] { 1, 1, 4 });
            p.AddCoeff(-1.4734759885548700e+03, new int[] { 1, 1, 6 });
            p.AddCoeff(7.8936213672582100e+02, new int[] { 1, 1, 8 });
            p.AddCoeff(-2.0035609789984600e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(7.2128195243944500e+02, new int[] { 3, 1, 2 });
            p.AddCoeff(-3.9670507384169500e+03, new int[] { 3, 1, 4 });
            p.AddCoeff(6.8762212799227100e+03, new int[] { 3, 1, 6 });
            p.AddCoeff(-3.6836899713871600e+03, new int[] { 3, 1, 8 });
            p.AddCoeff(1.8032048810986100e+01, new int[] { 5, 1, 0 });
            p.AddCoeff(-6.4915375719550000e+02, new int[] { 5, 1, 2 });
            p.AddCoeff(3.5703456645752500e+03, new int[] { 5, 1, 4 });
            p.AddCoeff(-6.1885991519304400e+03, new int[] { 5, 1, 6 });
            p.AddCoeff(3.3153209742484500e+03, new int[] { 5, 1, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("8d001965-f8ec-4c66-a555-fe8362afd1f3"));
            OrthonormalPolynomials[627] = p;
            p.AddCoeff(2.0825782043134200e+01, new int[] { 1, 0, 1 });
            p.AddCoeff(-1.8743203838820800e+02, new int[] { 1, 0, 3 });
            p.AddCoeff(4.1235048445405700e+02, new int[] { 1, 0, 5 });
            p.AddCoeff(-2.5526458561441600e+02, new int[] { 1, 0, 7 });
            p.AddCoeff(-6.2477346129402500e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(5.6229611516462300e+02, new int[] { 1, 2, 3 });
            p.AddCoeff(-1.2370514533621700e+03, new int[] { 1, 2, 5 });
            p.AddCoeff(7.6579375684324800e+02, new int[] { 1, 2, 7 });
            p.AddCoeff(-9.7186982867959500e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(8.7468284581163500e+02, new int[] { 3, 0, 3 });
            p.AddCoeff(-1.9243022607856000e+03, new int[] { 3, 0, 5 });
            p.AddCoeff(1.1912347328672700e+03, new int[] { 3, 0, 7 });
            p.AddCoeff(2.9156094860387800e+02, new int[] { 3, 2, 1 });
            p.AddCoeff(-2.6240485374349100e+03, new int[] { 3, 2, 3 });
            p.AddCoeff(5.7729067823567900e+03, new int[] { 3, 2, 5 });
            p.AddCoeff(-3.5737041986018200e+03, new int[] { 3, 2, 7 });
            p.AddCoeff(8.7468284581163500e+01, new int[] { 5, 0, 1 });
            p.AddCoeff(-7.8721456123047200e+02, new int[] { 5, 0, 3 });
            p.AddCoeff(1.7318720347070400e+03, new int[] { 5, 0, 5 });
            p.AddCoeff(-1.0721112595805500e+03, new int[] { 5, 0, 7 });
            p.AddCoeff(-2.6240485374349100e+02, new int[] { 5, 2, 1 });
            p.AddCoeff(2.3616436836914200e+03, new int[] { 5, 2, 3 });
            p.AddCoeff(-5.1956161041211100e+03, new int[] { 5, 2, 5 });
            p.AddCoeff(3.2163337787416400e+03, new int[] { 5, 2, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("8ed29e14-45f4-434e-82ee-b06405c2ad15"));
            OrthonormalPolynomials[628] = p;
            p.AddCoeff(9.8313826118541900e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-2.0645903484893800e+02, new int[] { 1, 1, 2 });
            p.AddCoeff(6.1937710454681400e+02, new int[] { 1, 1, 4 });
            p.AddCoeff(-4.5420987666766400e+02, new int[] { 1, 1, 6 });
            p.AddCoeff(-1.6385637686423700e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(3.4409839141489700e+02, new int[] { 1, 3, 2 });
            p.AddCoeff(-1.0322951742446900e+03, new int[] { 1, 3, 4 });
            p.AddCoeff(7.5701646111277300e+02, new int[] { 1, 3, 6 });
            p.AddCoeff(-4.5879785521986200e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(9.6347549596171100e+02, new int[] { 3, 1, 2 });
            p.AddCoeff(-2.8904264878851300e+03, new int[] { 3, 1, 4 });
            p.AddCoeff(2.1196460911157600e+03, new int[] { 3, 1, 6 });
            p.AddCoeff(7.6466309203310400e+01, new int[] { 3, 3, 0 });
            p.AddCoeff(-1.6057924932695200e+03, new int[] { 3, 3, 2 });
            p.AddCoeff(4.8173774798085600e+03, new int[] { 3, 3, 4 });
            p.AddCoeff(-3.5327434851929400e+03, new int[] { 3, 3, 6 });
            p.AddCoeff(4.1291806969787600e+01, new int[] { 5, 1, 0 });
            p.AddCoeff(-8.6712794636554000e+02, new int[] { 5, 1, 2 });
            p.AddCoeff(2.6013838390966200e+03, new int[] { 5, 1, 4 });
            p.AddCoeff(-1.9076814820041900e+03, new int[] { 5, 1, 6 });
            p.AddCoeff(-6.8819678282979400e+01, new int[] { 5, 3, 0 });
            p.AddCoeff(1.4452132439425700e+03, new int[] { 5, 3, 2 });
            p.AddCoeff(-4.3356397318277000e+03, new int[] { 5, 3, 4 });
            p.AddCoeff(3.1794691366736500e+03, new int[] { 5, 3, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("19d775e5-2681-4f19-96eb-f425aa4d15d7"));
            OrthonormalPolynomials[629] = p;
            p.AddCoeff(1.5381644092705500e+01, new int[] { 1, 0, 1 });
            p.AddCoeff(-7.1781005765958900e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(6.4602905189363000e+01, new int[] { 1, 0, 5 });
            p.AddCoeff(-1.5381644092705500e+02, new int[] { 1, 2, 1 });
            p.AddCoeff(7.1781005765958900e+02, new int[] { 1, 2, 3 });
            p.AddCoeff(-6.4602905189363000e+02, new int[] { 1, 2, 5 });
            p.AddCoeff(1.7945251441489700e+02, new int[] { 1, 4, 1 });
            p.AddCoeff(-8.3744506726952000e+02, new int[] { 1, 4, 3 });
            p.AddCoeff(7.5370056054256800e+02, new int[] { 1, 4, 5 });
            p.AddCoeff(-7.1781005765958900e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(3.3497802690780800e+02, new int[] { 3, 0, 3 });
            p.AddCoeff(-3.0148022421702700e+02, new int[] { 3, 0, 5 });
            p.AddCoeff(7.1781005765958900e+02, new int[] { 3, 2, 1 });
            p.AddCoeff(-3.3497802690780800e+03, new int[] { 3, 2, 3 });
            p.AddCoeff(3.0148022421702700e+03, new int[] { 3, 2, 5 });
            p.AddCoeff(-8.3744506726952000e+02, new int[] { 3, 4, 1 });
            p.AddCoeff(3.9080769805910900e+03, new int[] { 3, 4, 3 });
            p.AddCoeff(-3.5172692825319800e+03, new int[] { 3, 4, 5 });
            p.AddCoeff(6.4602905189363000e+01, new int[] { 5, 0, 1 });
            p.AddCoeff(-3.0148022421702700e+02, new int[] { 5, 0, 3 });
            p.AddCoeff(2.7133220179532400e+02, new int[] { 5, 0, 5 });
            p.AddCoeff(-6.4602905189363000e+02, new int[] { 5, 2, 1 });
            p.AddCoeff(3.0148022421702700e+03, new int[] { 5, 2, 3 });
            p.AddCoeff(-2.7133220179532400e+03, new int[] { 5, 2, 5 });
            p.AddCoeff(7.5370056054256800e+02, new int[] { 5, 4, 1 });
            p.AddCoeff(-3.5172692825319800e+03, new int[] { 5, 4, 3 });
            p.AddCoeff(3.1655423542787900e+03, new int[] { 5, 4, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("bd1d05e1-0c4c-4079-a5c3-07c5917d368b"));
            OrthonormalPolynomials[630] = p;
            p.AddCoeff(1.5381644092705500e+01, new int[] { 1, 1, 0 });
            p.AddCoeff(-1.5381644092705500e+02, new int[] { 1, 1, 2 });
            p.AddCoeff(1.7945251441489700e+02, new int[] { 1, 1, 4 });
            p.AddCoeff(-7.1781005765958900e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(7.1781005765958900e+02, new int[] { 1, 3, 2 });
            p.AddCoeff(-8.3744506726952000e+02, new int[] { 1, 3, 4 });
            p.AddCoeff(6.4602905189363000e+01, new int[] { 1, 5, 0 });
            p.AddCoeff(-6.4602905189363000e+02, new int[] { 1, 5, 2 });
            p.AddCoeff(7.5370056054256800e+02, new int[] { 1, 5, 4 });
            p.AddCoeff(-7.1781005765958900e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(7.1781005765958900e+02, new int[] { 3, 1, 2 });
            p.AddCoeff(-8.3744506726952000e+02, new int[] { 3, 1, 4 });
            p.AddCoeff(3.3497802690780800e+02, new int[] { 3, 3, 0 });
            p.AddCoeff(-3.3497802690780800e+03, new int[] { 3, 3, 2 });
            p.AddCoeff(3.9080769805910900e+03, new int[] { 3, 3, 4 });
            p.AddCoeff(-3.0148022421702700e+02, new int[] { 3, 5, 0 });
            p.AddCoeff(3.0148022421702700e+03, new int[] { 3, 5, 2 });
            p.AddCoeff(-3.5172692825319800e+03, new int[] { 3, 5, 4 });
            p.AddCoeff(6.4602905189363000e+01, new int[] { 5, 1, 0 });
            p.AddCoeff(-6.4602905189363000e+02, new int[] { 5, 1, 2 });
            p.AddCoeff(7.5370056054256800e+02, new int[] { 5, 1, 4 });
            p.AddCoeff(-3.0148022421702700e+02, new int[] { 5, 3, 0 });
            p.AddCoeff(3.0148022421702700e+03, new int[] { 5, 3, 2 });
            p.AddCoeff(-3.5172692825319800e+03, new int[] { 5, 3, 4 });
            p.AddCoeff(2.7133220179532400e+02, new int[] { 5, 5, 0 });
            p.AddCoeff(-2.7133220179532400e+03, new int[] { 5, 5, 2 });
            p.AddCoeff(3.1655423542787900e+03, new int[] { 5, 5, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("6584d54a-9166-4b15-a681-c9a1f3b1b11a"));
            OrthonormalPolynomials[631] = p;
            p.AddCoeff(9.8313826118541900e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-1.6385637686423700e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(-2.0645903484893800e+02, new int[] { 1, 2, 1 });
            p.AddCoeff(3.4409839141489700e+02, new int[] { 1, 2, 3 });
            p.AddCoeff(6.1937710454681400e+02, new int[] { 1, 4, 1 });
            p.AddCoeff(-1.0322951742446900e+03, new int[] { 1, 4, 3 });
            p.AddCoeff(-4.5420987666766400e+02, new int[] { 1, 6, 1 });
            p.AddCoeff(7.5701646111277300e+02, new int[] { 1, 6, 3 });
            p.AddCoeff(-4.5879785521986200e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(7.6466309203310400e+01, new int[] { 3, 0, 3 });
            p.AddCoeff(9.6347549596171100e+02, new int[] { 3, 2, 1 });
            p.AddCoeff(-1.6057924932695200e+03, new int[] { 3, 2, 3 });
            p.AddCoeff(-2.8904264878851300e+03, new int[] { 3, 4, 1 });
            p.AddCoeff(4.8173774798085600e+03, new int[] { 3, 4, 3 });
            p.AddCoeff(2.1196460911157600e+03, new int[] { 3, 6, 1 });
            p.AddCoeff(-3.5327434851929400e+03, new int[] { 3, 6, 3 });
            p.AddCoeff(4.1291806969787600e+01, new int[] { 5, 0, 1 });
            p.AddCoeff(-6.8819678282979400e+01, new int[] { 5, 0, 3 });
            p.AddCoeff(-8.6712794636554000e+02, new int[] { 5, 2, 1 });
            p.AddCoeff(1.4452132439425700e+03, new int[] { 5, 2, 3 });
            p.AddCoeff(2.6013838390966200e+03, new int[] { 5, 4, 1 });
            p.AddCoeff(-4.3356397318277000e+03, new int[] { 5, 4, 3 });
            p.AddCoeff(-1.9076814820041900e+03, new int[] { 5, 6, 1 });
            p.AddCoeff(3.1794691366736500e+03, new int[] { 5, 6, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("1a36abd5-d9c7-4cdb-88aa-7742a2c8ed03"));
            OrthonormalPolynomials[632] = p;
            p.AddCoeff(2.0825782043134200e+01, new int[] { 1, 1, 0 });
            p.AddCoeff(-6.2477346129402500e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(-1.8743203838820800e+02, new int[] { 1, 3, 0 });
            p.AddCoeff(5.6229611516462300e+02, new int[] { 1, 3, 2 });
            p.AddCoeff(4.1235048445405700e+02, new int[] { 1, 5, 0 });
            p.AddCoeff(-1.2370514533621700e+03, new int[] { 1, 5, 2 });
            p.AddCoeff(-2.5526458561441600e+02, new int[] { 1, 7, 0 });
            p.AddCoeff(7.6579375684324800e+02, new int[] { 1, 7, 2 });
            p.AddCoeff(-9.7186982867959500e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(2.9156094860387800e+02, new int[] { 3, 1, 2 });
            p.AddCoeff(8.7468284581163500e+02, new int[] { 3, 3, 0 });
            p.AddCoeff(-2.6240485374349100e+03, new int[] { 3, 3, 2 });
            p.AddCoeff(-1.9243022607856000e+03, new int[] { 3, 5, 0 });
            p.AddCoeff(5.7729067823567900e+03, new int[] { 3, 5, 2 });
            p.AddCoeff(1.1912347328672700e+03, new int[] { 3, 7, 0 });
            p.AddCoeff(-3.5737041986018200e+03, new int[] { 3, 7, 2 });
            p.AddCoeff(8.7468284581163500e+01, new int[] { 5, 1, 0 });
            p.AddCoeff(-2.6240485374349100e+02, new int[] { 5, 1, 2 });
            p.AddCoeff(-7.8721456123047200e+02, new int[] { 5, 3, 0 });
            p.AddCoeff(2.3616436836914200e+03, new int[] { 5, 3, 2 });
            p.AddCoeff(1.7318720347070400e+03, new int[] { 5, 5, 0 });
            p.AddCoeff(-5.1956161041211100e+03, new int[] { 5, 5, 2 });
            p.AddCoeff(-1.0721112595805500e+03, new int[] { 5, 7, 0 });
            p.AddCoeff(3.2163337787416400e+03, new int[] { 5, 7, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("952a3c50-5ff4-4296-adec-4717dd95063a"));
            OrthonormalPolynomials[633] = p;
            p.AddCoeff(4.2933449549966900e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-1.5456041837988100e+02, new int[] { 1, 2, 1 });
            p.AddCoeff(8.5008230108934600e+02, new int[] { 1, 4, 1 });
            p.AddCoeff(-1.4734759885548700e+03, new int[] { 1, 6, 1 });
            p.AddCoeff(7.8936213672582100e+02, new int[] { 1, 8, 1 });
            p.AddCoeff(-2.0035609789984600e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(7.2128195243944500e+02, new int[] { 3, 2, 1 });
            p.AddCoeff(-3.9670507384169500e+03, new int[] { 3, 4, 1 });
            p.AddCoeff(6.8762212799227100e+03, new int[] { 3, 6, 1 });
            p.AddCoeff(-3.6836899713871600e+03, new int[] { 3, 8, 1 });
            p.AddCoeff(1.8032048810986100e+01, new int[] { 5, 0, 1 });
            p.AddCoeff(-6.4915375719550000e+02, new int[] { 5, 2, 1 });
            p.AddCoeff(3.5703456645752500e+03, new int[] { 5, 4, 1 });
            p.AddCoeff(-6.1885991519304400e+03, new int[] { 5, 6, 1 });
            p.AddCoeff(3.3153209742484500e+03, new int[] { 5, 8, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("8f4522fc-6beb-4587-8f7d-7c8d039cb4d1"));
            OrthonormalPolynomials[634] = p;
            p.AddCoeff(2.3584680961604700e+01, new int[] { 1, 1, 0 });
            p.AddCoeff(-3.4590865410353500e+02, new int[] { 1, 3, 0 });
            p.AddCoeff(1.3490437510037900e+03, new int[] { 1, 5, 0 });
            p.AddCoeff(-1.9272053585768400e+03, new int[] { 1, 7, 0 });
            p.AddCoeff(9.1006919710573000e+02, new int[] { 1, 9, 0 });
            p.AddCoeff(-1.1006184448748800e+02, new int[] { 3, 1, 0 });
            p.AddCoeff(1.6142403858165000e+03, new int[] { 3, 3, 0 });
            p.AddCoeff(-6.2955375046843400e+03, new int[] { 3, 5, 0 });
            p.AddCoeff(8.9936250066919200e+03, new int[] { 3, 7, 0 });
            p.AddCoeff(-4.2469895864934000e+03, new int[] { 3, 9, 0 });
            p.AddCoeff(9.9055660038739600e+01, new int[] { 5, 1, 0 });
            p.AddCoeff(-1.4528163472348500e+03, new int[] { 5, 3, 0 });
            p.AddCoeff(5.6659837542159100e+03, new int[] { 5, 5, 0 });
            p.AddCoeff(-8.0942625060227200e+03, new int[] { 5, 7, 0 });
            p.AddCoeff(3.8222906278440600e+03, new int[] { 5, 9, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("4b954ce5-e882-46a2-a97b-de968da90589"));
            OrthonormalPolynomials[635] = p;
            p.AddCoeff(-4.4911673672912700e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.6168202522248600e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(-8.8925113872367100e+01, new int[] { 0, 0, 4 });
            p.AddCoeff(1.5413686404543600e+02, new int[] { 0, 0, 6 });
            p.AddCoeff(-8.2573320024340800e+01, new int[] { 0, 0, 8 });
            p.AddCoeff(9.4314514713116600e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-3.3953225296722000e+02, new int[] { 2, 0, 2 });
            p.AddCoeff(1.8674273913197100e+03, new int[] { 2, 0, 4 });
            p.AddCoeff(-3.2368741449541600e+03, new int[] { 2, 0, 6 });
            p.AddCoeff(1.7340397205111600e+03, new int[] { 2, 0, 8 });
            p.AddCoeff(-2.8294354413935000e+01, new int[] { 4, 0, 0 });
            p.AddCoeff(1.0185967589016600e+03, new int[] { 4, 0, 2 });
            p.AddCoeff(-5.6022821739591300e+03, new int[] { 4, 0, 4 });
            p.AddCoeff(9.7106224348624800e+03, new int[] { 4, 0, 6 });
            p.AddCoeff(-5.2021191615334600e+03, new int[] { 4, 0, 8 });
            p.AddCoeff(2.0749193236885600e+01, new int[] { 6, 0, 0 });
            p.AddCoeff(-7.4697095652788300e+02, new int[] { 6, 0, 2 });
            p.AddCoeff(4.1083402609033600e+03, new int[] { 6, 0, 4 });
            p.AddCoeff(-7.1211231188991500e+03, new int[] { 6, 0, 6 });
            p.AddCoeff(3.8148873851245400e+03, new int[] { 6, 0, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("fdcb8eac-858c-4c5b-8099-6b57f30c67ff"));
            OrthonormalPolynomials[636] = p;
            p.AddCoeff(5.8456259587602400e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-5.2610633628842200e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(1.1574339398345300e+02, new int[] { 0, 1, 5 });
            p.AddCoeff(-7.1650672465946900e+01, new int[] { 0, 1, 7 });
            p.AddCoeff(-1.2275814513396500e+02, new int[] { 2, 1, 1 });
            p.AddCoeff(1.1048233062056900e+03, new int[] { 2, 1, 3 });
            p.AddCoeff(-2.4306112736525100e+03, new int[] { 2, 1, 5 });
            p.AddCoeff(1.5046641217848900e+03, new int[] { 2, 1, 7 });
            p.AddCoeff(3.6827443540189500e+02, new int[] { 4, 1, 1 });
            p.AddCoeff(-3.3144699186170600e+03, new int[] { 4, 1, 3 });
            p.AddCoeff(7.2918338209575200e+03, new int[] { 4, 1, 5 });
            p.AddCoeff(-4.5139923653546600e+03, new int[] { 4, 1, 7 });
            p.AddCoeff(-2.7006791929472300e+02, new int[] { 6, 1, 1 });
            p.AddCoeff(2.4306112736525100e+03, new int[] { 6, 1, 3 });
            p.AddCoeff(-5.3473448020355200e+03, new int[] { 6, 1, 5 });
            p.AddCoeff(3.3102610679267500e+03, new int[] { 6, 1, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("75294572-902f-4e9f-a3db-c81e7b0adb1d"));
            OrthonormalPolynomials[637] = p;
            p.AddCoeff(-5.0182628884508000e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.0538352065746700e+01, new int[] { 0, 0, 2 });
            p.AddCoeff(-3.1615056197240000e+01, new int[] { 0, 0, 4 });
            p.AddCoeff(2.3184374544642700e+01, new int[] { 0, 0, 6 });
            p.AddCoeff(1.5054788665352400e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-3.1615056197240000e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(9.4845168591720100e+01, new int[] { 0, 2, 4 });
            p.AddCoeff(-6.9553123633928100e+01, new int[] { 0, 2, 6 });
            p.AddCoeff(1.0538352065746700e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(-2.2130539338068000e+02, new int[] { 2, 0, 2 });
            p.AddCoeff(6.6391618014204100e+02, new int[] { 2, 0, 4 });
            p.AddCoeff(-4.8687186543749600e+02, new int[] { 2, 0, 6 });
            p.AddCoeff(-3.1615056197240000e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(6.6391618014204100e+02, new int[] { 2, 2, 2 });
            p.AddCoeff(-1.9917485404261200e+03, new int[] { 2, 2, 4 });
            p.AddCoeff(1.4606155963124900e+03, new int[] { 2, 2, 6 });
            p.AddCoeff(-3.1615056197240000e+01, new int[] { 4, 0, 0 });
            p.AddCoeff(6.6391618014204100e+02, new int[] { 4, 0, 2 });
            p.AddCoeff(-1.9917485404261200e+03, new int[] { 4, 0, 4 });
            p.AddCoeff(1.4606155963124900e+03, new int[] { 4, 0, 6 });
            p.AddCoeff(9.4845168591720100e+01, new int[] { 4, 2, 0 });
            p.AddCoeff(-1.9917485404261200e+03, new int[] { 4, 2, 2 });
            p.AddCoeff(5.9752456212783700e+03, new int[] { 4, 2, 4 });
            p.AddCoeff(-4.3818467889374700e+03, new int[] { 4, 2, 6 });
            p.AddCoeff(2.3184374544642700e+01, new int[] { 6, 0, 0 });
            p.AddCoeff(-4.8687186543749600e+02, new int[] { 6, 0, 2 });
            p.AddCoeff(1.4606155963124900e+03, new int[] { 6, 0, 4 });
            p.AddCoeff(-1.0711181039624900e+03, new int[] { 6, 0, 6 });
            p.AddCoeff(-6.9553123633928100e+01, new int[] { 6, 2, 0 });
            p.AddCoeff(1.4606155963124900e+03, new int[] { 6, 2, 2 });
            p.AddCoeff(-4.3818467889374700e+03, new int[] { 6, 2, 4 });
            p.AddCoeff(3.2133543118874800e+03, new int[] { 6, 2, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7f20d37e-725a-44ae-aee5-5889cc9cf5d1"));
            OrthonormalPolynomials[638] = p;
            p.AddCoeff(9.8313826118541900e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-4.5879785521986200e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(4.1291806969787600e+01, new int[] { 0, 1, 5 });
            p.AddCoeff(-1.6385637686423700e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(7.6466309203310400e+01, new int[] { 0, 3, 3 });
            p.AddCoeff(-6.8819678282979400e+01, new int[] { 0, 3, 5 });
            p.AddCoeff(-2.0645903484893800e+02, new int[] { 2, 1, 1 });
            p.AddCoeff(9.6347549596171100e+02, new int[] { 2, 1, 3 });
            p.AddCoeff(-8.6712794636554000e+02, new int[] { 2, 1, 5 });
            p.AddCoeff(3.4409839141489700e+02, new int[] { 2, 3, 1 });
            p.AddCoeff(-1.6057924932695200e+03, new int[] { 2, 3, 3 });
            p.AddCoeff(1.4452132439425700e+03, new int[] { 2, 3, 5 });
            p.AddCoeff(6.1937710454681400e+02, new int[] { 4, 1, 1 });
            p.AddCoeff(-2.8904264878851300e+03, new int[] { 4, 1, 3 });
            p.AddCoeff(2.6013838390966200e+03, new int[] { 4, 1, 5 });
            p.AddCoeff(-1.0322951742446900e+03, new int[] { 4, 3, 1 });
            p.AddCoeff(4.8173774798085600e+03, new int[] { 4, 3, 3 });
            p.AddCoeff(-4.3356397318277000e+03, new int[] { 4, 3, 5 });
            p.AddCoeff(-4.5420987666766400e+02, new int[] { 6, 1, 1 });
            p.AddCoeff(2.1196460911157600e+03, new int[] { 6, 1, 3 });
            p.AddCoeff(-1.9076814820041900e+03, new int[] { 6, 1, 5 });
            p.AddCoeff(7.5701646111277300e+02, new int[] { 6, 3, 1 });
            p.AddCoeff(-3.5327434851929400e+03, new int[] { 6, 3, 3 });
            p.AddCoeff(3.1794691366736500e+03, new int[] { 6, 3, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("44cfacec-14a7-4980-8d6e-90fc4a32cc54"));
            OrthonormalPolynomials[639] = p;
            p.AddCoeff(-5.0417551342897400e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(5.0417551342897400e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-5.8820476566713600e+00, new int[] { 0, 0, 4 });
            p.AddCoeff(5.0417551342897400e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-5.0417551342897400e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(5.8820476566713600e+01, new int[] { 0, 2, 4 });
            p.AddCoeff(-5.8820476566713600e+00, new int[] { 0, 4, 0 });
            p.AddCoeff(5.8820476566713600e+01, new int[] { 0, 4, 2 });
            p.AddCoeff(-6.8623889327832500e+01, new int[] { 0, 4, 4 });
            p.AddCoeff(1.0587685782008400e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(-1.0587685782008400e+02, new int[] { 2, 0, 2 });
            p.AddCoeff(1.2352300079009900e+02, new int[] { 2, 0, 4 });
            p.AddCoeff(-1.0587685782008400e+02, new int[] { 2, 2, 0 });
            p.AddCoeff(1.0587685782008400e+03, new int[] { 2, 2, 2 });
            p.AddCoeff(-1.2352300079009900e+03, new int[] { 2, 2, 4 });
            p.AddCoeff(1.2352300079009900e+02, new int[] { 2, 4, 0 });
            p.AddCoeff(-1.2352300079009900e+03, new int[] { 2, 4, 2 });
            p.AddCoeff(1.4411016758844800e+03, new int[] { 2, 4, 4 });
            p.AddCoeff(-3.1763057346025300e+01, new int[] { 4, 0, 0 });
            p.AddCoeff(3.1763057346025300e+02, new int[] { 4, 0, 2 });
            p.AddCoeff(-3.7056900237029600e+02, new int[] { 4, 0, 4 });
            p.AddCoeff(3.1763057346025300e+02, new int[] { 4, 2, 0 });
            p.AddCoeff(-3.1763057346025300e+03, new int[] { 4, 2, 2 });
            p.AddCoeff(3.7056900237029600e+03, new int[] { 4, 2, 4 });
            p.AddCoeff(-3.7056900237029600e+02, new int[] { 4, 4, 0 });
            p.AddCoeff(3.7056900237029600e+03, new int[] { 4, 4, 2 });
            p.AddCoeff(-4.3233050276534500e+03, new int[] { 4, 4, 4 });
            p.AddCoeff(2.3292908720418600e+01, new int[] { 6, 0, 0 });
            p.AddCoeff(-2.3292908720418600e+02, new int[] { 6, 0, 2 });
            p.AddCoeff(2.7175060173821700e+02, new int[] { 6, 0, 4 });
            p.AddCoeff(-2.3292908720418600e+02, new int[] { 6, 2, 0 });
            p.AddCoeff(2.3292908720418600e+03, new int[] { 6, 2, 2 });
            p.AddCoeff(-2.7175060173821700e+03, new int[] { 6, 2, 4 });
            p.AddCoeff(2.7175060173821700e+02, new int[] { 6, 4, 0 });
            p.AddCoeff(-2.7175060173821700e+03, new int[] { 6, 4, 2 });
            p.AddCoeff(3.1704236869458600e+03, new int[] { 6, 4, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("93fd2451-163d-4a3e-af56-0973b45d2232"));
            OrthonormalPolynomials[640] = p;
            p.AddCoeff(9.8313826118541900e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-1.6385637686423700e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(-4.5879785521986200e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(7.6466309203310400e+01, new int[] { 0, 3, 3 });
            p.AddCoeff(4.1291806969787600e+01, new int[] { 0, 5, 1 });
            p.AddCoeff(-6.8819678282979400e+01, new int[] { 0, 5, 3 });
            p.AddCoeff(-2.0645903484893800e+02, new int[] { 2, 1, 1 });
            p.AddCoeff(3.4409839141489700e+02, new int[] { 2, 1, 3 });
            p.AddCoeff(9.6347549596171100e+02, new int[] { 2, 3, 1 });
            p.AddCoeff(-1.6057924932695200e+03, new int[] { 2, 3, 3 });
            p.AddCoeff(-8.6712794636554000e+02, new int[] { 2, 5, 1 });
            p.AddCoeff(1.4452132439425700e+03, new int[] { 2, 5, 3 });
            p.AddCoeff(6.1937710454681400e+02, new int[] { 4, 1, 1 });
            p.AddCoeff(-1.0322951742446900e+03, new int[] { 4, 1, 3 });
            p.AddCoeff(-2.8904264878851300e+03, new int[] { 4, 3, 1 });
            p.AddCoeff(4.8173774798085600e+03, new int[] { 4, 3, 3 });
            p.AddCoeff(2.6013838390966200e+03, new int[] { 4, 5, 1 });
            p.AddCoeff(-4.3356397318277000e+03, new int[] { 4, 5, 3 });
            p.AddCoeff(-4.5420987666766400e+02, new int[] { 6, 1, 1 });
            p.AddCoeff(7.5701646111277300e+02, new int[] { 6, 1, 3 });
            p.AddCoeff(2.1196460911157600e+03, new int[] { 6, 3, 1 });
            p.AddCoeff(-3.5327434851929400e+03, new int[] { 6, 3, 3 });
            p.AddCoeff(-1.9076814820041900e+03, new int[] { 6, 5, 1 });
            p.AddCoeff(3.1794691366736500e+03, new int[] { 6, 5, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("21712316-69e9-4ad1-a306-6fd0fa4cef93"));
            OrthonormalPolynomials[641] = p;
            p.AddCoeff(-5.0182628884508000e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.5054788665352400e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(1.0538352065746700e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(-3.1615056197240000e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(-3.1615056197240000e+01, new int[] { 0, 4, 0 });
            p.AddCoeff(9.4845168591720100e+01, new int[] { 0, 4, 2 });
            p.AddCoeff(2.3184374544642700e+01, new int[] { 0, 6, 0 });
            p.AddCoeff(-6.9553123633928100e+01, new int[] { 0, 6, 2 });
            p.AddCoeff(1.0538352065746700e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(-3.1615056197240000e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(-2.2130539338068000e+02, new int[] { 2, 2, 0 });
            p.AddCoeff(6.6391618014204100e+02, new int[] { 2, 2, 2 });
            p.AddCoeff(6.6391618014204100e+02, new int[] { 2, 4, 0 });
            p.AddCoeff(-1.9917485404261200e+03, new int[] { 2, 4, 2 });
            p.AddCoeff(-4.8687186543749600e+02, new int[] { 2, 6, 0 });
            p.AddCoeff(1.4606155963124900e+03, new int[] { 2, 6, 2 });
            p.AddCoeff(-3.1615056197240000e+01, new int[] { 4, 0, 0 });
            p.AddCoeff(9.4845168591720100e+01, new int[] { 4, 0, 2 });
            p.AddCoeff(6.6391618014204100e+02, new int[] { 4, 2, 0 });
            p.AddCoeff(-1.9917485404261200e+03, new int[] { 4, 2, 2 });
            p.AddCoeff(-1.9917485404261200e+03, new int[] { 4, 4, 0 });
            p.AddCoeff(5.9752456212783700e+03, new int[] { 4, 4, 2 });
            p.AddCoeff(1.4606155963124900e+03, new int[] { 4, 6, 0 });
            p.AddCoeff(-4.3818467889374700e+03, new int[] { 4, 6, 2 });
            p.AddCoeff(2.3184374544642700e+01, new int[] { 6, 0, 0 });
            p.AddCoeff(-6.9553123633928100e+01, new int[] { 6, 0, 2 });
            p.AddCoeff(-4.8687186543749600e+02, new int[] { 6, 2, 0 });
            p.AddCoeff(1.4606155963124900e+03, new int[] { 6, 2, 2 });
            p.AddCoeff(1.4606155963124900e+03, new int[] { 6, 4, 0 });
            p.AddCoeff(-4.3818467889374700e+03, new int[] { 6, 4, 2 });
            p.AddCoeff(-1.0711181039624900e+03, new int[] { 6, 6, 0 });
            p.AddCoeff(3.2133543118874800e+03, new int[] { 6, 6, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b671768f-b810-4c27-b37b-11a6f5f0109f"));
            OrthonormalPolynomials[642] = p;
            p.AddCoeff(5.8456259587602400e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-5.2610633628842200e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(1.1574339398345300e+02, new int[] { 0, 5, 1 });
            p.AddCoeff(-7.1650672465946900e+01, new int[] { 0, 7, 1 });
            p.AddCoeff(-1.2275814513396500e+02, new int[] { 2, 1, 1 });
            p.AddCoeff(1.1048233062056900e+03, new int[] { 2, 3, 1 });
            p.AddCoeff(-2.4306112736525100e+03, new int[] { 2, 5, 1 });
            p.AddCoeff(1.5046641217848900e+03, new int[] { 2, 7, 1 });
            p.AddCoeff(3.6827443540189500e+02, new int[] { 4, 1, 1 });
            p.AddCoeff(-3.3144699186170600e+03, new int[] { 4, 3, 1 });
            p.AddCoeff(7.2918338209575200e+03, new int[] { 4, 5, 1 });
            p.AddCoeff(-4.5139923653546600e+03, new int[] { 4, 7, 1 });
            p.AddCoeff(-2.7006791929472300e+02, new int[] { 6, 1, 1 });
            p.AddCoeff(2.4306112736525100e+03, new int[] { 6, 3, 1 });
            p.AddCoeff(-5.3473448020355200e+03, new int[] { 6, 5, 1 });
            p.AddCoeff(3.3102610679267500e+03, new int[] { 6, 7, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7900b4c6-e056-491b-b274-fe940ed534c8"));
            OrthonormalPolynomials[643] = p;
            p.AddCoeff(-4.4911673672912700e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.6168202522248600e+01, new int[] { 0, 2, 0 });
            p.AddCoeff(-8.8925113872367100e+01, new int[] { 0, 4, 0 });
            p.AddCoeff(1.5413686404543600e+02, new int[] { 0, 6, 0 });
            p.AddCoeff(-8.2573320024340800e+01, new int[] { 0, 8, 0 });
            p.AddCoeff(9.4314514713116600e+00, new int[] { 2, 0, 0 });
            p.AddCoeff(-3.3953225296722000e+02, new int[] { 2, 2, 0 });
            p.AddCoeff(1.8674273913197100e+03, new int[] { 2, 4, 0 });
            p.AddCoeff(-3.2368741449541600e+03, new int[] { 2, 6, 0 });
            p.AddCoeff(1.7340397205111600e+03, new int[] { 2, 8, 0 });
            p.AddCoeff(-2.8294354413935000e+01, new int[] { 4, 0, 0 });
            p.AddCoeff(1.0185967589016600e+03, new int[] { 4, 2, 0 });
            p.AddCoeff(-5.6022821739591300e+03, new int[] { 4, 4, 0 });
            p.AddCoeff(9.7106224348624800e+03, new int[] { 4, 6, 0 });
            p.AddCoeff(-5.2021191615334600e+03, new int[] { 4, 8, 0 });
            p.AddCoeff(2.0749193236885600e+01, new int[] { 6, 0, 0 });
            p.AddCoeff(-7.4697095652788300e+02, new int[] { 6, 2, 0 });
            p.AddCoeff(4.1083402609033600e+03, new int[] { 6, 4, 0 });
            p.AddCoeff(-7.1211231188991500e+03, new int[] { 6, 6, 0 });
            p.AddCoeff(3.8148873851245400e+03, new int[] { 6, 8, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c7d17ced-7d0c-46ee-be60-70b5195f3513"));
            OrthonormalPolynomials[644] = p;
            p.AddCoeff(2.5377123250591500e+01, new int[] { 1, 0, 1 });
            p.AddCoeff(-2.2839410925532400e+02, new int[] { 1, 0, 3 });
            p.AddCoeff(5.0246704036171200e+02, new int[] { 1, 0, 5 });
            p.AddCoeff(-3.1105102498582200e+02, new int[] { 1, 0, 7 });
            p.AddCoeff(-2.2839410925532400e+02, new int[] { 3, 0, 1 });
            p.AddCoeff(2.0555469832979100e+03, new int[] { 3, 0, 3 });
            p.AddCoeff(-4.5222033632554100e+03, new int[] { 3, 0, 5 });
            p.AddCoeff(2.7994592248724000e+03, new int[] { 3, 0, 7 });
            p.AddCoeff(5.0246704036171200e+02, new int[] { 5, 0, 1 });
            p.AddCoeff(-4.5222033632554100e+03, new int[] { 5, 0, 3 });
            p.AddCoeff(9.9488473991619000e+03, new int[] { 5, 0, 5 });
            p.AddCoeff(-6.1588102947192700e+03, new int[] { 5, 0, 7 });
            p.AddCoeff(-3.1105102498582200e+02, new int[] { 7, 0, 1 });
            p.AddCoeff(2.7994592248724000e+03, new int[] { 7, 0, 3 });
            p.AddCoeff(-6.1588102947192700e+03, new int[] { 7, 0, 5 });
            p.AddCoeff(3.8125968491119300e+03, new int[] { 7, 0, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("dd2bba95-09cc-4409-b01d-f7c4d7f1e05a"));
            OrthonormalPolynomials[645] = p;
            p.AddCoeff(5.8456259587602400e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-1.2275814513396500e+02, new int[] { 1, 1, 2 });
            p.AddCoeff(3.6827443540189500e+02, new int[] { 1, 1, 4 });
            p.AddCoeff(-2.7006791929472300e+02, new int[] { 1, 1, 6 });
            p.AddCoeff(-5.2610633628842200e+01, new int[] { 3, 1, 0 });
            p.AddCoeff(1.1048233062056900e+03, new int[] { 3, 1, 2 });
            p.AddCoeff(-3.3144699186170600e+03, new int[] { 3, 1, 4 });
            p.AddCoeff(2.4306112736525100e+03, new int[] { 3, 1, 6 });
            p.AddCoeff(1.1574339398345300e+02, new int[] { 5, 1, 0 });
            p.AddCoeff(-2.4306112736525100e+03, new int[] { 5, 1, 2 });
            p.AddCoeff(7.2918338209575200e+03, new int[] { 5, 1, 4 });
            p.AddCoeff(-5.3473448020355200e+03, new int[] { 5, 1, 6 });
            p.AddCoeff(-7.1650672465946900e+01, new int[] { 7, 1, 0 });
            p.AddCoeff(1.5046641217848900e+03, new int[] { 7, 1, 2 });
            p.AddCoeff(-4.5139923653546600e+03, new int[] { 7, 1, 4 });
            p.AddCoeff(3.3102610679267500e+03, new int[] { 7, 1, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7a4f687e-d1b2-499e-8194-a804f1ef2087"));
            OrthonormalPolynomials[646] = p;
            p.AddCoeff(2.0825782043134200e+01, new int[] { 1, 0, 1 });
            p.AddCoeff(-9.7186982867959500e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(8.7468284581163500e+01, new int[] { 1, 0, 5 });
            p.AddCoeff(-6.2477346129402500e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(2.9156094860387800e+02, new int[] { 1, 2, 3 });
            p.AddCoeff(-2.6240485374349100e+02, new int[] { 1, 2, 5 });
            p.AddCoeff(-1.8743203838820800e+02, new int[] { 3, 0, 1 });
            p.AddCoeff(8.7468284581163500e+02, new int[] { 3, 0, 3 });
            p.AddCoeff(-7.8721456123047200e+02, new int[] { 3, 0, 5 });
            p.AddCoeff(5.6229611516462300e+02, new int[] { 3, 2, 1 });
            p.AddCoeff(-2.6240485374349100e+03, new int[] { 3, 2, 3 });
            p.AddCoeff(2.3616436836914200e+03, new int[] { 3, 2, 5 });
            p.AddCoeff(4.1235048445405700e+02, new int[] { 5, 0, 1 });
            p.AddCoeff(-1.9243022607856000e+03, new int[] { 5, 0, 3 });
            p.AddCoeff(1.7318720347070400e+03, new int[] { 5, 0, 5 });
            p.AddCoeff(-1.2370514533621700e+03, new int[] { 5, 2, 1 });
            p.AddCoeff(5.7729067823567900e+03, new int[] { 5, 2, 3 });
            p.AddCoeff(-5.1956161041211100e+03, new int[] { 5, 2, 5 });
            p.AddCoeff(-2.5526458561441600e+02, new int[] { 7, 0, 1 });
            p.AddCoeff(1.1912347328672700e+03, new int[] { 7, 0, 3 });
            p.AddCoeff(-1.0721112595805500e+03, new int[] { 7, 0, 5 });
            p.AddCoeff(7.6579375684324800e+02, new int[] { 7, 2, 1 });
            p.AddCoeff(-3.5737041986018200e+03, new int[] { 7, 2, 3 });
            p.AddCoeff(3.2163337787416400e+03, new int[] { 7, 2, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b3c0654f-a749-4021-a306-410a0e4712e3"));
            OrthonormalPolynomials[647] = p;
            p.AddCoeff(1.3373389672997100e+01, new int[] { 1, 1, 0 });
            p.AddCoeff(-1.3373389672997100e+02, new int[] { 1, 1, 2 });
            p.AddCoeff(1.5602287951829900e+02, new int[] { 1, 1, 4 });
            p.AddCoeff(-2.2288982788328400e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(2.2288982788328400e+02, new int[] { 1, 3, 2 });
            p.AddCoeff(-2.6003813253049800e+02, new int[] { 1, 3, 4 });
            p.AddCoeff(-1.2036050705697300e+02, new int[] { 3, 1, 0 });
            p.AddCoeff(1.2036050705697300e+03, new int[] { 3, 1, 2 });
            p.AddCoeff(-1.4042059156646900e+03, new int[] { 3, 1, 4 });
            p.AddCoeff(2.0060084509495600e+02, new int[] { 3, 3, 0 });
            p.AddCoeff(-2.0060084509495600e+03, new int[] { 3, 3, 2 });
            p.AddCoeff(2.3403431927744800e+03, new int[] { 3, 3, 4 });
            p.AddCoeff(2.6479311552534200e+02, new int[] { 5, 1, 0 });
            p.AddCoeff(-2.6479311552534200e+03, new int[] { 5, 1, 2 });
            p.AddCoeff(3.0892530144623200e+03, new int[] { 5, 1, 4 });
            p.AddCoeff(-4.4132185920890300e+02, new int[] { 5, 3, 0 });
            p.AddCoeff(4.4132185920890300e+03, new int[] { 5, 3, 2 });
            p.AddCoeff(-5.1487550241038700e+03, new int[] { 5, 3, 4 });
            p.AddCoeff(-1.6391954770616400e+02, new int[] { 7, 1, 0 });
            p.AddCoeff(1.6391954770616400e+03, new int[] { 7, 1, 2 });
            p.AddCoeff(-1.9123947232385800e+03, new int[] { 7, 1, 4 });
            p.AddCoeff(2.7319924617694000e+02, new int[] { 7, 3, 0 });
            p.AddCoeff(-2.7319924617694000e+03, new int[] { 7, 3, 2 });
            p.AddCoeff(3.1873245387309600e+03, new int[] { 7, 3, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f538a117-10fa-4350-8f78-d5812a82c4e6"));
            OrthonormalPolynomials[648] = p;
            p.AddCoeff(1.3373389672997100e+01, new int[] { 1, 0, 1 });
            p.AddCoeff(-2.2288982788328400e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(-1.3373389672997100e+02, new int[] { 1, 2, 1 });
            p.AddCoeff(2.2288982788328400e+02, new int[] { 1, 2, 3 });
            p.AddCoeff(1.5602287951829900e+02, new int[] { 1, 4, 1 });
            p.AddCoeff(-2.6003813253049800e+02, new int[] { 1, 4, 3 });
            p.AddCoeff(-1.2036050705697300e+02, new int[] { 3, 0, 1 });
            p.AddCoeff(2.0060084509495600e+02, new int[] { 3, 0, 3 });
            p.AddCoeff(1.2036050705697300e+03, new int[] { 3, 2, 1 });
            p.AddCoeff(-2.0060084509495600e+03, new int[] { 3, 2, 3 });
            p.AddCoeff(-1.4042059156646900e+03, new int[] { 3, 4, 1 });
            p.AddCoeff(2.3403431927744800e+03, new int[] { 3, 4, 3 });
            p.AddCoeff(2.6479311552534200e+02, new int[] { 5, 0, 1 });
            p.AddCoeff(-4.4132185920890300e+02, new int[] { 5, 0, 3 });
            p.AddCoeff(-2.6479311552534200e+03, new int[] { 5, 2, 1 });
            p.AddCoeff(4.4132185920890300e+03, new int[] { 5, 2, 3 });
            p.AddCoeff(3.0892530144623200e+03, new int[] { 5, 4, 1 });
            p.AddCoeff(-5.1487550241038700e+03, new int[] { 5, 4, 3 });
            p.AddCoeff(-1.6391954770616400e+02, new int[] { 7, 0, 1 });
            p.AddCoeff(2.7319924617694000e+02, new int[] { 7, 0, 3 });
            p.AddCoeff(1.6391954770616400e+03, new int[] { 7, 2, 1 });
            p.AddCoeff(-2.7319924617694000e+03, new int[] { 7, 2, 3 });
            p.AddCoeff(-1.9123947232385800e+03, new int[] { 7, 4, 1 });
            p.AddCoeff(3.1873245387309600e+03, new int[] { 7, 4, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("484b4eaf-1fda-44fe-b440-153536fcdb2c"));
            OrthonormalPolynomials[649] = p;
            p.AddCoeff(2.0825782043134200e+01, new int[] { 1, 1, 0 });
            p.AddCoeff(-6.2477346129402500e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(-9.7186982867959500e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(2.9156094860387800e+02, new int[] { 1, 3, 2 });
            p.AddCoeff(8.7468284581163500e+01, new int[] { 1, 5, 0 });
            p.AddCoeff(-2.6240485374349100e+02, new int[] { 1, 5, 2 });
            p.AddCoeff(-1.8743203838820800e+02, new int[] { 3, 1, 0 });
            p.AddCoeff(5.6229611516462300e+02, new int[] { 3, 1, 2 });
            p.AddCoeff(8.7468284581163500e+02, new int[] { 3, 3, 0 });
            p.AddCoeff(-2.6240485374349100e+03, new int[] { 3, 3, 2 });
            p.AddCoeff(-7.8721456123047200e+02, new int[] { 3, 5, 0 });
            p.AddCoeff(2.3616436836914200e+03, new int[] { 3, 5, 2 });
            p.AddCoeff(4.1235048445405700e+02, new int[] { 5, 1, 0 });
            p.AddCoeff(-1.2370514533621700e+03, new int[] { 5, 1, 2 });
            p.AddCoeff(-1.9243022607856000e+03, new int[] { 5, 3, 0 });
            p.AddCoeff(5.7729067823567900e+03, new int[] { 5, 3, 2 });
            p.AddCoeff(1.7318720347070400e+03, new int[] { 5, 5, 0 });
            p.AddCoeff(-5.1956161041211100e+03, new int[] { 5, 5, 2 });
            p.AddCoeff(-2.5526458561441600e+02, new int[] { 7, 1, 0 });
            p.AddCoeff(7.6579375684324800e+02, new int[] { 7, 1, 2 });
            p.AddCoeff(1.1912347328672700e+03, new int[] { 7, 3, 0 });
            p.AddCoeff(-3.5737041986018200e+03, new int[] { 7, 3, 2 });
            p.AddCoeff(-1.0721112595805500e+03, new int[] { 7, 5, 0 });
            p.AddCoeff(3.2163337787416400e+03, new int[] { 7, 5, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5bf6c6d1-b972-4d7f-805d-670411826c4d"));
            OrthonormalPolynomials[650] = p;
            p.AddCoeff(5.8456259587602400e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-1.2275814513396500e+02, new int[] { 1, 2, 1 });
            p.AddCoeff(3.6827443540189500e+02, new int[] { 1, 4, 1 });
            p.AddCoeff(-2.7006791929472300e+02, new int[] { 1, 6, 1 });
            p.AddCoeff(-5.2610633628842200e+01, new int[] { 3, 0, 1 });
            p.AddCoeff(1.1048233062056900e+03, new int[] { 3, 2, 1 });
            p.AddCoeff(-3.3144699186170600e+03, new int[] { 3, 4, 1 });
            p.AddCoeff(2.4306112736525100e+03, new int[] { 3, 6, 1 });
            p.AddCoeff(1.1574339398345300e+02, new int[] { 5, 0, 1 });
            p.AddCoeff(-2.4306112736525100e+03, new int[] { 5, 2, 1 });
            p.AddCoeff(7.2918338209575200e+03, new int[] { 5, 4, 1 });
            p.AddCoeff(-5.3473448020355200e+03, new int[] { 5, 6, 1 });
            p.AddCoeff(-7.1650672465946900e+01, new int[] { 7, 0, 1 });
            p.AddCoeff(1.5046641217848900e+03, new int[] { 7, 2, 1 });
            p.AddCoeff(-4.5139923653546600e+03, new int[] { 7, 4, 1 });
            p.AddCoeff(3.3102610679267500e+03, new int[] { 7, 6, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("fe2f1870-abf4-44a3-9f75-a010a7ef4383"));
            OrthonormalPolynomials[651] = p;
            p.AddCoeff(2.5377123250591500e+01, new int[] { 1, 1, 0 });
            p.AddCoeff(-2.2839410925532400e+02, new int[] { 1, 3, 0 });
            p.AddCoeff(5.0246704036171200e+02, new int[] { 1, 5, 0 });
            p.AddCoeff(-3.1105102498582200e+02, new int[] { 1, 7, 0 });
            p.AddCoeff(-2.2839410925532400e+02, new int[] { 3, 1, 0 });
            p.AddCoeff(2.0555469832979100e+03, new int[] { 3, 3, 0 });
            p.AddCoeff(-4.5222033632554100e+03, new int[] { 3, 5, 0 });
            p.AddCoeff(2.7994592248724000e+03, new int[] { 3, 7, 0 });
            p.AddCoeff(5.0246704036171200e+02, new int[] { 5, 1, 0 });
            p.AddCoeff(-4.5222033632554100e+03, new int[] { 5, 3, 0 });
            p.AddCoeff(9.9488473991619000e+03, new int[] { 5, 5, 0 });
            p.AddCoeff(-6.1588102947192700e+03, new int[] { 5, 7, 0 });
            p.AddCoeff(-3.1105102498582200e+02, new int[] { 7, 1, 0 });
            p.AddCoeff(2.7994592248724000e+03, new int[] { 7, 3, 0 });
            p.AddCoeff(-6.1588102947192700e+03, new int[] { 7, 5, 0 });
            p.AddCoeff(3.8125968491119300e+03, new int[] { 7, 7, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("85381387-ac84-4eb3-97cb-d84c5a76d561"));
            OrthonormalPolynomials[652] = p;
            p.AddCoeff(-4.4911673672912700e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(9.4314514713116600e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-2.8294354413935000e+01, new int[] { 0, 0, 4 });
            p.AddCoeff(2.0749193236885600e+01, new int[] { 0, 0, 6 });
            p.AddCoeff(1.6168202522248600e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(-3.3953225296722000e+02, new int[] { 2, 0, 2 });
            p.AddCoeff(1.0185967589016600e+03, new int[] { 2, 0, 4 });
            p.AddCoeff(-7.4697095652788300e+02, new int[] { 2, 0, 6 });
            p.AddCoeff(-8.8925113872367100e+01, new int[] { 4, 0, 0 });
            p.AddCoeff(1.8674273913197100e+03, new int[] { 4, 0, 2 });
            p.AddCoeff(-5.6022821739591300e+03, new int[] { 4, 0, 4 });
            p.AddCoeff(4.1083402609033600e+03, new int[] { 4, 0, 6 });
            p.AddCoeff(1.5413686404543600e+02, new int[] { 6, 0, 0 });
            p.AddCoeff(-3.2368741449541600e+03, new int[] { 6, 0, 2 });
            p.AddCoeff(9.7106224348624800e+03, new int[] { 6, 0, 4 });
            p.AddCoeff(-7.1211231188991500e+03, new int[] { 6, 0, 6 });
            p.AddCoeff(-8.2573320024340800e+01, new int[] { 8, 0, 0 });
            p.AddCoeff(1.7340397205111600e+03, new int[] { 8, 0, 2 });
            p.AddCoeff(-5.2021191615334600e+03, new int[] { 8, 0, 4 });
            p.AddCoeff(3.8148873851245400e+03, new int[] { 8, 0, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("ea309a7a-7870-444f-aa02-6d5168f53674"));
            OrthonormalPolynomials[653] = p;
            p.AddCoeff(4.2933449549966900e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-2.0035609789984600e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(1.8032048810986100e+01, new int[] { 0, 1, 5 });
            p.AddCoeff(-1.5456041837988100e+02, new int[] { 2, 1, 1 });
            p.AddCoeff(7.2128195243944500e+02, new int[] { 2, 1, 3 });
            p.AddCoeff(-6.4915375719550000e+02, new int[] { 2, 1, 5 });
            p.AddCoeff(8.5008230108934600e+02, new int[] { 4, 1, 1 });
            p.AddCoeff(-3.9670507384169500e+03, new int[] { 4, 1, 3 });
            p.AddCoeff(3.5703456645752500e+03, new int[] { 4, 1, 5 });
            p.AddCoeff(-1.4734759885548700e+03, new int[] { 6, 1, 1 });
            p.AddCoeff(6.8762212799227100e+03, new int[] { 6, 1, 3 });
            p.AddCoeff(-6.1885991519304400e+03, new int[] { 6, 1, 5 });
            p.AddCoeff(7.8936213672582100e+02, new int[] { 8, 1, 1 });
            p.AddCoeff(-3.6836899713871600e+03, new int[] { 8, 1, 3 });
            p.AddCoeff(3.3153209742484500e+03, new int[] { 8, 1, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("edf1b8cd-6f5e-44cd-bcbb-d872c3510ccf"));
            OrthonormalPolynomials[654] = p;
            p.AddCoeff(-5.0135467715791900e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(5.0135467715791900e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-5.8491379001757200e+00, new int[] { 0, 0, 4 });
            p.AddCoeff(1.5040640314737600e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-1.5040640314737600e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(1.7547413700527200e+01, new int[] { 0, 2, 4 });
            p.AddCoeff(1.8048768377685100e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(-1.8048768377685100e+02, new int[] { 2, 0, 2 });
            p.AddCoeff(2.1056896440632600e+02, new int[] { 2, 0, 4 });
            p.AddCoeff(-5.4146305133055200e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(5.4146305133055200e+02, new int[] { 2, 2, 2 });
            p.AddCoeff(-6.3170689321897700e+02, new int[] { 2, 2, 4 });
            p.AddCoeff(-9.9268226077267900e+01, new int[] { 4, 0, 0 });
            p.AddCoeff(9.9268226077267900e+02, new int[] { 4, 0, 2 });
            p.AddCoeff(-1.1581293042347900e+03, new int[] { 4, 0, 4 });
            p.AddCoeff(2.9780467823180400e+02, new int[] { 4, 2, 0 });
            p.AddCoeff(-2.9780467823180400e+03, new int[] { 4, 2, 2 });
            p.AddCoeff(3.4743879127043800e+03, new int[] { 4, 2, 4 });
            p.AddCoeff(1.7206492520059800e+02, new int[] { 6, 0, 0 });
            p.AddCoeff(-1.7206492520059800e+03, new int[] { 6, 0, 2 });
            p.AddCoeff(2.0074241273403100e+03, new int[] { 6, 0, 4 });
            p.AddCoeff(-5.1619477560179300e+02, new int[] { 6, 2, 0 });
            p.AddCoeff(5.1619477560179300e+03, new int[] { 6, 2, 2 });
            p.AddCoeff(-6.0222723820209200e+03, new int[] { 6, 2, 4 });
            p.AddCoeff(-9.2177638500320200e+01, new int[] { 8, 0, 0 });
            p.AddCoeff(9.2177638500320200e+02, new int[] { 8, 0, 2 });
            p.AddCoeff(-1.0754057825037400e+03, new int[] { 8, 0, 4 });
            p.AddCoeff(2.7653291550096100e+02, new int[] { 8, 2, 0 });
            p.AddCoeff(-2.7653291550096100e+03, new int[] { 8, 2, 2 });
            p.AddCoeff(3.2262173475112000e+03, new int[] { 8, 2, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("217ca0af-c7fa-4f55-9bb4-215c12242253"));
            OrthonormalPolynomials[655] = p;
            p.AddCoeff(6.2779535781903700e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-1.0463255963650600e+01, new int[] { 0, 1, 3 });
            p.AddCoeff(-1.0463255963650600e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(1.7438759939417700e+01, new int[] { 0, 3, 3 });
            p.AddCoeff(-2.2600632881485300e+02, new int[] { 2, 1, 1 });
            p.AddCoeff(3.7667721469142200e+02, new int[] { 2, 1, 3 });
            p.AddCoeff(3.7667721469142200e+02, new int[] { 2, 3, 1 });
            p.AddCoeff(-6.2779535781903700e+02, new int[] { 2, 3, 3 });
            p.AddCoeff(1.2430348084816900e+03, new int[] { 4, 1, 1 });
            p.AddCoeff(-2.0717246808028200e+03, new int[] { 4, 1, 3 });
            p.AddCoeff(-2.0717246808028200e+03, new int[] { 4, 3, 1 });
            p.AddCoeff(3.4528744680047100e+03, new int[] { 4, 3, 3 });
            p.AddCoeff(-2.1545936680349400e+03, new int[] { 6, 1, 1 });
            p.AddCoeff(3.5909894467248900e+03, new int[] { 6, 1, 3 });
            p.AddCoeff(3.5909894467248900e+03, new int[] { 6, 3, 1 });
            p.AddCoeff(-5.9849824112081600e+03, new int[] { 6, 3, 3 });
            p.AddCoeff(1.1542466078758600e+03, new int[] { 8, 1, 1 });
            p.AddCoeff(-1.9237443464597600e+03, new int[] { 8, 1, 3 });
            p.AddCoeff(-1.9237443464597600e+03, new int[] { 8, 3, 1 });
            p.AddCoeff(3.2062405774329400e+03, new int[] { 8, 3, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("15b8ffe8-5006-4dbf-912f-c50b6b6ece38"));
            OrthonormalPolynomials[656] = p;
            p.AddCoeff(-5.0135467715791900e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.5040640314737600e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(5.0135467715791900e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-1.5040640314737600e+01, new int[] { 0, 2, 2 });
            p.AddCoeff(-5.8491379001757200e+00, new int[] { 0, 4, 0 });
            p.AddCoeff(1.7547413700527200e+01, new int[] { 0, 4, 2 });
            p.AddCoeff(1.8048768377685100e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(-5.4146305133055200e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(-1.8048768377685100e+02, new int[] { 2, 2, 0 });
            p.AddCoeff(5.4146305133055200e+02, new int[] { 2, 2, 2 });
            p.AddCoeff(2.1056896440632600e+02, new int[] { 2, 4, 0 });
            p.AddCoeff(-6.3170689321897700e+02, new int[] { 2, 4, 2 });
            p.AddCoeff(-9.9268226077267900e+01, new int[] { 4, 0, 0 });
            p.AddCoeff(2.9780467823180400e+02, new int[] { 4, 0, 2 });
            p.AddCoeff(9.9268226077267900e+02, new int[] { 4, 2, 0 });
            p.AddCoeff(-2.9780467823180400e+03, new int[] { 4, 2, 2 });
            p.AddCoeff(-1.1581293042347900e+03, new int[] { 4, 4, 0 });
            p.AddCoeff(3.4743879127043800e+03, new int[] { 4, 4, 2 });
            p.AddCoeff(1.7206492520059800e+02, new int[] { 6, 0, 0 });
            p.AddCoeff(-5.1619477560179300e+02, new int[] { 6, 0, 2 });
            p.AddCoeff(-1.7206492520059800e+03, new int[] { 6, 2, 0 });
            p.AddCoeff(5.1619477560179300e+03, new int[] { 6, 2, 2 });
            p.AddCoeff(2.0074241273403100e+03, new int[] { 6, 4, 0 });
            p.AddCoeff(-6.0222723820209200e+03, new int[] { 6, 4, 2 });
            p.AddCoeff(-9.2177638500320200e+01, new int[] { 8, 0, 0 });
            p.AddCoeff(2.7653291550096100e+02, new int[] { 8, 0, 2 });
            p.AddCoeff(9.2177638500320200e+02, new int[] { 8, 2, 0 });
            p.AddCoeff(-2.7653291550096100e+03, new int[] { 8, 2, 2 });
            p.AddCoeff(-1.0754057825037400e+03, new int[] { 8, 4, 0 });
            p.AddCoeff(3.2262173475112000e+03, new int[] { 8, 4, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b584c8be-e20b-4fff-b114-c71e06e72a26"));
            OrthonormalPolynomials[657] = p;
            p.AddCoeff(4.2933449549966900e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-2.0035609789984600e+01, new int[] { 0, 3, 1 });
            p.AddCoeff(1.8032048810986100e+01, new int[] { 0, 5, 1 });
            p.AddCoeff(-1.5456041837988100e+02, new int[] { 2, 1, 1 });
            p.AddCoeff(7.2128195243944500e+02, new int[] { 2, 3, 1 });
            p.AddCoeff(-6.4915375719550000e+02, new int[] { 2, 5, 1 });
            p.AddCoeff(8.5008230108934600e+02, new int[] { 4, 1, 1 });
            p.AddCoeff(-3.9670507384169500e+03, new int[] { 4, 3, 1 });
            p.AddCoeff(3.5703456645752500e+03, new int[] { 4, 5, 1 });
            p.AddCoeff(-1.4734759885548700e+03, new int[] { 6, 1, 1 });
            p.AddCoeff(6.8762212799227100e+03, new int[] { 6, 3, 1 });
            p.AddCoeff(-6.1885991519304400e+03, new int[] { 6, 5, 1 });
            p.AddCoeff(7.8936213672582100e+02, new int[] { 8, 1, 1 });
            p.AddCoeff(-3.6836899713871600e+03, new int[] { 8, 3, 1 });
            p.AddCoeff(3.3153209742484500e+03, new int[] { 8, 5, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("06fc3a6f-5d64-4f4c-bca7-73b060383924"));
            OrthonormalPolynomials[658] = p;
            p.AddCoeff(-4.4911673672912700e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(9.4314514713116600e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-2.8294354413935000e+01, new int[] { 0, 4, 0 });
            p.AddCoeff(2.0749193236885600e+01, new int[] { 0, 6, 0 });
            p.AddCoeff(1.6168202522248600e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(-3.3953225296722000e+02, new int[] { 2, 2, 0 });
            p.AddCoeff(1.0185967589016600e+03, new int[] { 2, 4, 0 });
            p.AddCoeff(-7.4697095652788300e+02, new int[] { 2, 6, 0 });
            p.AddCoeff(-8.8925113872367100e+01, new int[] { 4, 0, 0 });
            p.AddCoeff(1.8674273913197100e+03, new int[] { 4, 2, 0 });
            p.AddCoeff(-5.6022821739591300e+03, new int[] { 4, 4, 0 });
            p.AddCoeff(4.1083402609033600e+03, new int[] { 4, 6, 0 });
            p.AddCoeff(1.5413686404543600e+02, new int[] { 6, 0, 0 });
            p.AddCoeff(-3.2368741449541600e+03, new int[] { 6, 2, 0 });
            p.AddCoeff(9.7106224348624800e+03, new int[] { 6, 4, 0 });
            p.AddCoeff(-7.1211231188991500e+03, new int[] { 6, 6, 0 });
            p.AddCoeff(-8.2573320024340800e+01, new int[] { 8, 0, 0 });
            p.AddCoeff(1.7340397205111600e+03, new int[] { 8, 2, 0 });
            p.AddCoeff(-5.2021191615334600e+03, new int[] { 8, 4, 0 });
            p.AddCoeff(3.8148873851245400e+03, new int[] { 8, 6, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("01111ee0-98fc-4eb6-a265-f3a579b235ea"));
            OrthonormalPolynomials[659] = p;
            p.AddCoeff(2.3584680961604700e+01, new int[] { 1, 0, 1 });
            p.AddCoeff(-1.1006184448748800e+02, new int[] { 1, 0, 3 });
            p.AddCoeff(9.9055660038739600e+01, new int[] { 1, 0, 5 });
            p.AddCoeff(-3.4590865410353500e+02, new int[] { 3, 0, 1 });
            p.AddCoeff(1.6142403858165000e+03, new int[] { 3, 0, 3 });
            p.AddCoeff(-1.4528163472348500e+03, new int[] { 3, 0, 5 });
            p.AddCoeff(1.3490437510037900e+03, new int[] { 5, 0, 1 });
            p.AddCoeff(-6.2955375046843400e+03, new int[] { 5, 0, 3 });
            p.AddCoeff(5.6659837542159100e+03, new int[] { 5, 0, 5 });
            p.AddCoeff(-1.9272053585768400e+03, new int[] { 7, 0, 1 });
            p.AddCoeff(8.9936250066919200e+03, new int[] { 7, 0, 3 });
            p.AddCoeff(-8.0942625060227200e+03, new int[] { 7, 0, 5 });
            p.AddCoeff(9.1006919710573000e+02, new int[] { 9, 0, 1 });
            p.AddCoeff(-4.2469895864934000e+03, new int[] { 9, 0, 3 });
            p.AddCoeff(3.8222906278440600e+03, new int[] { 9, 0, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("aa7c7bf5-5d83-45ea-865b-c6e7a115d13c"));
            OrthonormalPolynomials[660] = p;
            p.AddCoeff(7.3900187608663900e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-7.3900187608663900e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(8.6216885543441200e+01, new int[] { 1, 1, 4 });
            p.AddCoeff(-1.0838694182604000e+02, new int[] { 3, 1, 0 });
            p.AddCoeff(1.0838694182604000e+03, new int[] { 3, 1, 2 });
            p.AddCoeff(-1.2645143213038000e+03, new int[] { 3, 1, 4 });
            p.AddCoeff(4.2270907312155800e+02, new int[] { 5, 1, 0 });
            p.AddCoeff(-4.2270907312155800e+03, new int[] { 5, 1, 2 });
            p.AddCoeff(4.9316058530848400e+03, new int[] { 5, 1, 4 });
            p.AddCoeff(-6.0387010445936800e+02, new int[] { 7, 1, 0 });
            p.AddCoeff(6.0387010445936800e+03, new int[] { 7, 1, 2 });
            p.AddCoeff(-7.0451512186926300e+03, new int[] { 7, 1, 4 });
            p.AddCoeff(2.8516088266136800e+02, new int[] { 9, 1, 0 });
            p.AddCoeff(-2.8516088266136800e+03, new int[] { 9, 1, 2 });
            p.AddCoeff(3.3268769643826300e+03, new int[] { 9, 1, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("dcf7e77a-4c1d-45e4-854e-93313da483b9"));
            OrthonormalPolynomials[661] = p;
            p.AddCoeff(1.6827812978247900e+01, new int[] { 1, 0, 1 });
            p.AddCoeff(-2.8046354963746500e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(-5.0483438934743800e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(8.4139064891239600e+01, new int[] { 1, 2, 3 });
            p.AddCoeff(-2.4680792368097000e+02, new int[] { 3, 0, 1 });
            p.AddCoeff(4.1134653946828300e+02, new int[] { 3, 0, 3 });
            p.AddCoeff(7.4042377104290900e+02, new int[] { 3, 2, 1 });
            p.AddCoeff(-1.2340396184048500e+03, new int[] { 3, 2, 3 });
            p.AddCoeff(9.6255090235578100e+02, new int[] { 5, 0, 1 });
            p.AddCoeff(-1.6042515039263000e+03, new int[] { 5, 0, 3 });
            p.AddCoeff(-2.8876527070673400e+03, new int[] { 5, 2, 1 });
            p.AddCoeff(4.8127545117789100e+03, new int[] { 5, 2, 3 });
            p.AddCoeff(-1.3750727176511200e+03, new int[] { 7, 0, 1 });
            p.AddCoeff(2.2917878627518600e+03, new int[] { 7, 0, 3 });
            p.AddCoeff(4.1252181529533500e+03, new int[] { 7, 2, 1 });
            p.AddCoeff(-6.8753635882555800e+03, new int[] { 7, 2, 3 });
            p.AddCoeff(6.4933989444636000e+02, new int[] { 9, 0, 1 });
            p.AddCoeff(-1.0822331574106000e+03, new int[] { 9, 0, 3 });
            p.AddCoeff(-1.9480196833390800e+03, new int[] { 9, 2, 1 });
            p.AddCoeff(3.2466994722318000e+03, new int[] { 9, 2, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f9258f14-d401-4b79-b962-00e3a28f8550"));
            OrthonormalPolynomials[662] = p;
            p.AddCoeff(1.6827812978247900e+01, new int[] { 1, 1, 0 });
            p.AddCoeff(-5.0483438934743800e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(-2.8046354963746500e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(8.4139064891239600e+01, new int[] { 1, 3, 2 });
            p.AddCoeff(-2.4680792368097000e+02, new int[] { 3, 1, 0 });
            p.AddCoeff(7.4042377104290900e+02, new int[] { 3, 1, 2 });
            p.AddCoeff(4.1134653946828300e+02, new int[] { 3, 3, 0 });
            p.AddCoeff(-1.2340396184048500e+03, new int[] { 3, 3, 2 });
            p.AddCoeff(9.6255090235578100e+02, new int[] { 5, 1, 0 });
            p.AddCoeff(-2.8876527070673400e+03, new int[] { 5, 1, 2 });
            p.AddCoeff(-1.6042515039263000e+03, new int[] { 5, 3, 0 });
            p.AddCoeff(4.8127545117789100e+03, new int[] { 5, 3, 2 });
            p.AddCoeff(-1.3750727176511200e+03, new int[] { 7, 1, 0 });
            p.AddCoeff(4.1252181529533500e+03, new int[] { 7, 1, 2 });
            p.AddCoeff(2.2917878627518600e+03, new int[] { 7, 3, 0 });
            p.AddCoeff(-6.8753635882555800e+03, new int[] { 7, 3, 2 });
            p.AddCoeff(6.4933989444636000e+02, new int[] { 9, 1, 0 });
            p.AddCoeff(-1.9480196833390800e+03, new int[] { 9, 1, 2 });
            p.AddCoeff(-1.0822331574106000e+03, new int[] { 9, 3, 0 });
            p.AddCoeff(3.2466994722318000e+03, new int[] { 9, 3, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c2de9977-8a31-47fa-b855-6e861ed47160"));
            OrthonormalPolynomials[663] = p;
            p.AddCoeff(7.3900187608663900e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-7.3900187608663900e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(8.6216885543441200e+01, new int[] { 1, 4, 1 });
            p.AddCoeff(-1.0838694182604000e+02, new int[] { 3, 0, 1 });
            p.AddCoeff(1.0838694182604000e+03, new int[] { 3, 2, 1 });
            p.AddCoeff(-1.2645143213038000e+03, new int[] { 3, 4, 1 });
            p.AddCoeff(4.2270907312155800e+02, new int[] { 5, 0, 1 });
            p.AddCoeff(-4.2270907312155800e+03, new int[] { 5, 2, 1 });
            p.AddCoeff(4.9316058530848400e+03, new int[] { 5, 4, 1 });
            p.AddCoeff(-6.0387010445936800e+02, new int[] { 7, 0, 1 });
            p.AddCoeff(6.0387010445936800e+03, new int[] { 7, 2, 1 });
            p.AddCoeff(-7.0451512186926300e+03, new int[] { 7, 4, 1 });
            p.AddCoeff(2.8516088266136800e+02, new int[] { 9, 0, 1 });
            p.AddCoeff(-2.8516088266136800e+03, new int[] { 9, 2, 1 });
            p.AddCoeff(3.3268769643826300e+03, new int[] { 9, 4, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("079c9d6c-3b58-43de-b157-51f1131e3902"));
            OrthonormalPolynomials[664] = p;
            p.AddCoeff(2.3584680961604700e+01, new int[] { 1, 1, 0 });
            p.AddCoeff(-1.1006184448748800e+02, new int[] { 1, 3, 0 });
            p.AddCoeff(9.9055660038739600e+01, new int[] { 1, 5, 0 });
            p.AddCoeff(-3.4590865410353500e+02, new int[] { 3, 1, 0 });
            p.AddCoeff(1.6142403858165000e+03, new int[] { 3, 3, 0 });
            p.AddCoeff(-1.4528163472348500e+03, new int[] { 3, 5, 0 });
            p.AddCoeff(1.3490437510037900e+03, new int[] { 5, 1, 0 });
            p.AddCoeff(-6.2955375046843400e+03, new int[] { 5, 3, 0 });
            p.AddCoeff(5.6659837542159100e+03, new int[] { 5, 5, 0 });
            p.AddCoeff(-1.9272053585768400e+03, new int[] { 7, 1, 0 });
            p.AddCoeff(8.9936250066919200e+03, new int[] { 7, 3, 0 });
            p.AddCoeff(-8.0942625060227200e+03, new int[] { 7, 5, 0 });
            p.AddCoeff(9.1006919710573000e+02, new int[] { 9, 1, 0 });
            p.AddCoeff(-4.2469895864934000e+03, new int[] { 9, 3, 0 });
            p.AddCoeff(3.8222906278440600e+03, new int[] { 9, 5, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e692d58c-7a85-43b2-b393-664d3f3bb432"));
            OrthonormalPolynomials[665] = p;
            p.AddCoeff(-4.4855712597622800e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(4.4855712597622800e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(-5.2331664697226600e+00, new int[] { 0, 0, 4 });
            p.AddCoeff(2.4670641928692500e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(-2.4670641928692500e+02, new int[] { 2, 0, 2 });
            p.AddCoeff(2.8782415583474600e+02, new int[] { 2, 0, 4 });
            p.AddCoeff(-2.1381223004866800e+02, new int[] { 4, 0, 0 });
            p.AddCoeff(2.1381223004866800e+03, new int[] { 4, 0, 2 });
            p.AddCoeff(-2.4944760172344700e+03, new int[] { 4, 0, 4 });
            p.AddCoeff(6.4143669014600500e+02, new int[] { 6, 0, 0 });
            p.AddCoeff(-6.4143669014600500e+03, new int[] { 6, 0, 2 });
            p.AddCoeff(7.4834280517033600e+03, new int[] { 6, 0, 4 });
            p.AddCoeff(-7.7888740946301000e+02, new int[] { 8, 0, 0 });
            p.AddCoeff(7.7888740946301000e+03, new int[] { 8, 0, 2 });
            p.AddCoeff(-9.0870197770684300e+03, new int[] { 8, 0, 4 });
            p.AddCoeff(3.2886357288438100e+02, new int[] { 10, 0, 0 });
            p.AddCoeff(-3.2886357288438100e+03, new int[] { 10, 0, 2 });
            p.AddCoeff(3.8367416836511000e+03, new int[] { 10, 0, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c783fbc2-1ab9-45be-9533-7a2010821446"));
            OrthonormalPolynomials[666] = p;
            p.AddCoeff(2.7407293110638800e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-4.5678821851064700e+00, new int[] { 0, 1, 3 });
            p.AddCoeff(-1.5074011210851400e+02, new int[] { 2, 1, 1 });
            p.AddCoeff(2.5123352018085600e+02, new int[] { 2, 1, 3 });
            p.AddCoeff(1.3064143049404500e+03, new int[] { 4, 1, 1 });
            p.AddCoeff(-2.1773571749007500e+03, new int[] { 4, 1, 3 });
            p.AddCoeff(-3.9192429148213500e+03, new int[] { 6, 1, 1 });
            p.AddCoeff(6.5320715247022600e+03, new int[] { 6, 1, 3 });
            p.AddCoeff(4.7590806822830700e+03, new int[] { 8, 1, 1 });
            p.AddCoeff(-7.9318011371384500e+03, new int[] { 8, 1, 3 });
            p.AddCoeff(-2.0093896214084100e+03, new int[] { 10, 1, 1 });
            p.AddCoeff(3.3489827023473500e+03, new int[] { 10, 1, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b627aea9-0eb8-4eec-b6a9-e1d5101120dd"));
            OrthonormalPolynomials[667] = p;
            p.AddCoeff(-4.9839680664025300e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.4951904199207600e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(1.4951904199207600e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-4.4855712597622800e+00, new int[] { 0, 2, 2 });
            p.AddCoeff(2.7411824365213900e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(-8.2235473095641700e+01, new int[] { 2, 0, 2 });
            p.AddCoeff(-8.2235473095641700e+01, new int[] { 2, 2, 0 });
            p.AddCoeff(2.4670641928692500e+02, new int[] { 2, 2, 2 });
            p.AddCoeff(-2.3756914449852100e+02, new int[] { 4, 0, 0 });
            p.AddCoeff(7.1270743349556200e+02, new int[] { 4, 0, 2 });
            p.AddCoeff(7.1270743349556200e+02, new int[] { 4, 2, 0 });
            p.AddCoeff(-2.1381223004866800e+03, new int[] { 4, 2, 2 });
            p.AddCoeff(7.1270743349556200e+02, new int[] { 6, 0, 0 });
            p.AddCoeff(-2.1381223004866800e+03, new int[] { 6, 0, 2 });
            p.AddCoeff(-2.1381223004866800e+03, new int[] { 6, 2, 0 });
            p.AddCoeff(6.4143669014600500e+03, new int[] { 6, 2, 2 });
            p.AddCoeff(-8.6543045495889600e+02, new int[] { 8, 0, 0 });
            p.AddCoeff(2.5962913648766900e+03, new int[] { 8, 0, 2 });
            p.AddCoeff(2.5962913648766900e+03, new int[] { 8, 2, 0 });
            p.AddCoeff(-7.7888740946301000e+03, new int[] { 8, 2, 2 });
            p.AddCoeff(3.6540396987153400e+02, new int[] { 10, 0, 0 });
            p.AddCoeff(-1.0962119096146000e+03, new int[] { 10, 0, 2 });
            p.AddCoeff(-1.0962119096146000e+03, new int[] { 10, 2, 0 });
            p.AddCoeff(3.2886357288438100e+03, new int[] { 10, 2, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("73e7296d-c499-4770-8f9f-f1493cb2562c"));
            OrthonormalPolynomials[668] = p;
            p.AddCoeff(2.7407293110638800e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-4.5678821851064700e+00, new int[] { 0, 3, 1 });
            p.AddCoeff(-1.5074011210851400e+02, new int[] { 2, 1, 1 });
            p.AddCoeff(2.5123352018085600e+02, new int[] { 2, 3, 1 });
            p.AddCoeff(1.3064143049404500e+03, new int[] { 4, 1, 1 });
            p.AddCoeff(-2.1773571749007500e+03, new int[] { 4, 3, 1 });
            p.AddCoeff(-3.9192429148213500e+03, new int[] { 6, 1, 1 });
            p.AddCoeff(6.5320715247022600e+03, new int[] { 6, 3, 1 });
            p.AddCoeff(4.7590806822830700e+03, new int[] { 8, 1, 1 });
            p.AddCoeff(-7.9318011371384500e+03, new int[] { 8, 3, 1 });
            p.AddCoeff(-2.0093896214084100e+03, new int[] { 10, 1, 1 });
            p.AddCoeff(3.3489827023473500e+03, new int[] { 10, 3, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("60d77be3-8f72-4643-bfe9-93ccb0da38de"));
            OrthonormalPolynomials[669] = p;
            p.AddCoeff(-4.4855712597622800e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(4.4855712597622800e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(-5.2331664697226600e+00, new int[] { 0, 4, 0 });
            p.AddCoeff(2.4670641928692500e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(-2.4670641928692500e+02, new int[] { 2, 2, 0 });
            p.AddCoeff(2.8782415583474600e+02, new int[] { 2, 4, 0 });
            p.AddCoeff(-2.1381223004866800e+02, new int[] { 4, 0, 0 });
            p.AddCoeff(2.1381223004866800e+03, new int[] { 4, 2, 0 });
            p.AddCoeff(-2.4944760172344700e+03, new int[] { 4, 4, 0 });
            p.AddCoeff(6.4143669014600500e+02, new int[] { 6, 0, 0 });
            p.AddCoeff(-6.4143669014600500e+03, new int[] { 6, 2, 0 });
            p.AddCoeff(7.4834280517033600e+03, new int[] { 6, 4, 0 });
            p.AddCoeff(-7.7888740946301000e+02, new int[] { 8, 0, 0 });
            p.AddCoeff(7.7888740946301000e+03, new int[] { 8, 2, 0 });
            p.AddCoeff(-9.0870197770684300e+03, new int[] { 8, 4, 0 });
            p.AddCoeff(3.2886357288438100e+02, new int[] { 10, 0, 0 });
            p.AddCoeff(-3.2886357288438100e+03, new int[] { 10, 2, 0 });
            p.AddCoeff(3.8367416836511000e+03, new int[] { 10, 4, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("8d78a3c6-7e31-44a7-82bc-f0fd3445a23e"));
            OrthonormalPolynomials[670] = p;
            p.AddCoeff(1.8215977151856400e+01, new int[] { 1, 0, 1 });
            p.AddCoeff(-3.0359961919760700e+01, new int[] { 1, 0, 3 });
            p.AddCoeff(-3.9467950495688900e+02, new int[] { 3, 0, 1 });
            p.AddCoeff(6.5779917492814900e+02, new int[] { 3, 0, 3 });
            p.AddCoeff(2.3680770297413400e+03, new int[] { 5, 0, 1 });
            p.AddCoeff(-3.9467950495688900e+03, new int[] { 5, 0, 3 });
            p.AddCoeff(-5.7510442150861000e+03, new int[] { 7, 0, 1 });
            p.AddCoeff(9.5850736918101700e+03, new int[] { 7, 0, 3 });
            p.AddCoeff(6.0705466714797800e+03, new int[] { 9, 0, 1 });
            p.AddCoeff(-1.0117577785799600e+04, new int[] { 9, 0, 3 });
            p.AddCoeff(-2.3178450927468200e+03, new int[] { 11, 0, 1 });
            p.AddCoeff(3.8630751545780400e+03, new int[] { 11, 0, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("edd96b5e-b031-4a7e-8bd2-37d0a98c39b6"));
            OrthonormalPolynomials[671] = p;
            p.AddCoeff(8.8884867156627400e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-2.6665460146988200e+01, new int[] { 1, 1, 2 });
            p.AddCoeff(-1.9258387883935900e+02, new int[] { 3, 1, 0 });
            p.AddCoeff(5.7775163651807800e+02, new int[] { 3, 1, 2 });
            p.AddCoeff(1.1555032730361600e+03, new int[] { 5, 1, 0 });
            p.AddCoeff(-3.4665098191084700e+03, new int[] { 5, 1, 2 });
            p.AddCoeff(-2.8062222345163800e+03, new int[] { 7, 1, 0 });
            p.AddCoeff(8.4186667035491400e+03, new int[] { 7, 1, 2 });
            p.AddCoeff(2.9621234697672900e+03, new int[] { 9, 1, 0 });
            p.AddCoeff(-8.8863704093018700e+03, new int[] { 9, 1, 2 });
            p.AddCoeff(-1.1309925975475100e+03, new int[] { 11, 1, 0 });
            p.AddCoeff(3.3929777926425300e+03, new int[] { 11, 1, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("6743bd57-5adb-4539-8821-01265566e869"));
            OrthonormalPolynomials[672] = p;
            p.AddCoeff(8.8884867156627400e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-2.6665460146988200e+01, new int[] { 1, 2, 1 });
            p.AddCoeff(-1.9258387883935900e+02, new int[] { 3, 0, 1 });
            p.AddCoeff(5.7775163651807800e+02, new int[] { 3, 2, 1 });
            p.AddCoeff(1.1555032730361600e+03, new int[] { 5, 0, 1 });
            p.AddCoeff(-3.4665098191084700e+03, new int[] { 5, 2, 1 });
            p.AddCoeff(-2.8062222345163800e+03, new int[] { 7, 0, 1 });
            p.AddCoeff(8.4186667035491400e+03, new int[] { 7, 2, 1 });
            p.AddCoeff(2.9621234697672900e+03, new int[] { 9, 0, 1 });
            p.AddCoeff(-8.8863704093018700e+03, new int[] { 9, 2, 1 });
            p.AddCoeff(-1.1309925975475100e+03, new int[] { 11, 0, 1 });
            p.AddCoeff(3.3929777926425300e+03, new int[] { 11, 2, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("1a0330f6-6964-4dae-afe0-e80e0cf30a38"));
            OrthonormalPolynomials[673] = p;
            p.AddCoeff(1.8215977151856400e+01, new int[] { 1, 1, 0 });
            p.AddCoeff(-3.0359961919760700e+01, new int[] { 1, 3, 0 });
            p.AddCoeff(-3.9467950495688900e+02, new int[] { 3, 1, 0 });
            p.AddCoeff(6.5779917492814900e+02, new int[] { 3, 3, 0 });
            p.AddCoeff(2.3680770297413400e+03, new int[] { 5, 1, 0 });
            p.AddCoeff(-3.9467950495688900e+03, new int[] { 5, 3, 0 });
            p.AddCoeff(-5.7510442150861000e+03, new int[] { 7, 1, 0 });
            p.AddCoeff(9.5850736918101700e+03, new int[] { 7, 3, 0 });
            p.AddCoeff(6.0705466714797800e+03, new int[] { 9, 1, 0 });
            p.AddCoeff(-1.0117577785799600e+04, new int[] { 9, 3, 0 });
            p.AddCoeff(-2.3178450927468200e+03, new int[] { 11, 1, 0 });
            p.AddCoeff(3.8630751545780400e+03, new int[] { 11, 3, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c6a1ec97-7371-4078-be80-14235210ec4d"));
            OrthonormalPolynomials[674] = p;
            p.AddCoeff(-4.4585335662774400e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.3375600698832300e+00, new int[] { 0, 0, 2 });
            p.AddCoeff(3.4776561816964000e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(-1.0432968545089200e+02, new int[] { 2, 0, 2 });
            p.AddCoeff(-4.3470702271204900e+02, new int[] { 4, 0, 0 });
            p.AddCoeff(1.3041210681361500e+03, new int[] { 4, 0, 2 });
            p.AddCoeff(1.9706718362946300e+03, new int[] { 6, 0, 0 });
            p.AddCoeff(-5.9120155088838900e+03, new int[] { 6, 0, 2 });
            p.AddCoeff(-4.0117248095997700e+03, new int[] { 8, 0, 0 });
            p.AddCoeff(1.2035174428799300e+04, new int[] { 8, 0, 2 });
            p.AddCoeff(3.7442764889597800e+03, new int[] { 10, 0, 0 });
            p.AddCoeff(-1.1232829466879400e+04, new int[] { 10, 0, 2 });
            p.AddCoeff(-1.3048236249405400e+03, new int[] { 12, 0, 0 });
            p.AddCoeff(3.9144708748216000e+03, new int[] { 12, 0, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("66ca584a-7391-4e44-ad5c-2e4e49cadc29"));
            OrthonormalPolynomials[675] = p;
            p.AddCoeff(1.1963500960993100e+00, new int[] { 0, 1, 1 });
            p.AddCoeff(-9.3315307495746500e+01, new int[] { 2, 1, 1 });
            p.AddCoeff(1.1664413436968300e+03, new int[] { 4, 1, 1 });
            p.AddCoeff(-5.2878674247589700e+03, new int[] { 6, 1, 1 });
            p.AddCoeff(1.0764587257545100e+04, new int[] { 8, 1, 1 });
            p.AddCoeff(-1.0046948107042000e+04, new int[] { 10, 1, 1 });
            p.AddCoeff(3.5012091888176700e+03, new int[] { 12, 1, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f5da5012-f127-4b0e-96f4-511fa004074b"));
            OrthonormalPolynomials[676] = p;
            p.AddCoeff(-4.4585335662774400e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(1.3375600698832300e+00, new int[] { 0, 2, 0 });
            p.AddCoeff(3.4776561816964000e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(-1.0432968545089200e+02, new int[] { 2, 2, 0 });
            p.AddCoeff(-4.3470702271204900e+02, new int[] { 4, 0, 0 });
            p.AddCoeff(1.3041210681361500e+03, new int[] { 4, 2, 0 });
            p.AddCoeff(1.9706718362946300e+03, new int[] { 6, 0, 0 });
            p.AddCoeff(-5.9120155088838900e+03, new int[] { 6, 2, 0 });
            p.AddCoeff(-4.0117248095997700e+03, new int[] { 8, 0, 0 });
            p.AddCoeff(1.2035174428799300e+04, new int[] { 8, 2, 0 });
            p.AddCoeff(3.7442764889597800e+03, new int[] { 10, 0, 0 });
            p.AddCoeff(-1.1232829466879400e+04, new int[] { 10, 2, 0 });
            p.AddCoeff(-1.3048236249405400e+03, new int[] { 12, 0, 0 });
            p.AddCoeff(3.9144708748216000e+03, new int[] { 12, 2, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("43513bd6-2072-471f-9c5a-94ade8effcfd"));
            OrthonormalPolynomials[677] = p;
            p.AddCoeff(9.3315307495746500e+00, new int[] { 1, 0, 1 });
            p.AddCoeff(-2.7994592248724000e+02, new int[] { 3, 0, 1 });
            p.AddCoeff(2.3795403411415300e+03, new int[] { 5, 0, 1 });
            p.AddCoeff(-8.6116698060360400e+03, new int[] { 7, 0, 1 });
            p.AddCoeff(1.5070422160563000e+04, new int[] { 9, 0, 1 });
            p.AddCoeff(-1.2604353079743700e+04, new int[] { 11, 0, 1 });
            p.AddCoeff(4.0398567563281000e+03, new int[] { 13, 0, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("6394e57d-2bd8-42b5-bf02-15e80f7da11b"));
            OrthonormalPolynomials[678] = p;
            p.AddCoeff(9.3315307495746500e+00, new int[] { 1, 1, 0 });
            p.AddCoeff(-2.7994592248724000e+02, new int[] { 3, 1, 0 });
            p.AddCoeff(2.3795403411415300e+03, new int[] { 5, 1, 0 });
            p.AddCoeff(-8.6116698060360400e+03, new int[] { 7, 1, 0 });
            p.AddCoeff(1.5070422160563000e+04, new int[] { 9, 1, 0 });
            p.AddCoeff(-1.2604353079743700e+04, new int[] { 11, 1, 0 });
            p.AddCoeff(4.0398567563281000e+03, new int[] { 13, 1, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("1ba6bd05-8ba3-44a0-b3dc-e5a81a351d4d"));
            OrthonormalPolynomials[679] = p;
            p.AddCoeff(-3.9882405547065600e-01, new int[] { 0, 0, 0 });
            p.AddCoeff(4.1876525824418900e+01, new int[] { 2, 0, 0 });
            p.AddCoeff(-7.1190093901512200e+02, new int[] { 4, 0, 0 });
            p.AddCoeff(4.5087059470957700e+03, new int[] { 6, 0, 0 });
            p.AddCoeff(-1.3526117841287300e+04, new int[] { 8, 0, 0 });
            p.AddCoeff(2.0740047356640600e+04, new int[] { 10, 0, 0 });
            p.AddCoeff(-1.5712157088364100e+04, new int[] { 12, 0, 0 });
            p.AddCoeff(4.6618488064376900e+03, new int[] { 14, 0, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7eebd232-1a2b-4699-95f3-6cb41880c883"));
            OrthonormalPolynomials[680] = p;
            p.AddCoeff(-6.1852100426350100e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(2.4534666502452200e+02, new int[] { 0, 0, 3 });
            p.AddCoeff(-2.7969519812795600e+03, new int[] { 0, 0, 5 });
            p.AddCoeff(1.3984759906397800e+04, new int[] { 0, 0, 7 });
            p.AddCoeff(-3.5738830871905400e+04, new int[] { 0, 0, 9 });
            p.AddCoeff(4.8734769370780000e+04, new int[] { 0, 0, 11 });
            p.AddCoeff(-3.3739455718232300e+04, new int[] { 0, 0, 13 });
            p.AddCoeff(9.3185163412260600e+03, new int[] { 0, 0, 15 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b511107e-8301-4cad-946a-01ed2e960f90"));
            OrthonormalPolynomials[681] = p;
            p.AddCoeff(-6.9078352735584400e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(7.2532270372363600e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-1.2330485963301800e+03, new int[] { 0, 1, 4 });
            p.AddCoeff(7.8093077767578100e+03, new int[] { 0, 1, 6 });
            p.AddCoeff(-2.3427923330273500e+04, new int[] { 0, 1, 8 });
            p.AddCoeff(3.5922815773086000e+04, new int[] { 0, 1, 10 });
            p.AddCoeff(-2.7214254373550000e+04, new int[] { 0, 1, 12 });
            p.AddCoeff(8.0745589899543900e+03, new int[] { 0, 1, 14 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("69930745-8c49-42a9-92f6-b4ccaf542084"));
            OrthonormalPolynomials[682] = p;
            p.AddCoeff(-6.0234771979541500e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(1.8070431593862500e+02, new int[] { 0, 0, 3 });
            p.AddCoeff(-1.5359866854783100e+03, new int[] { 0, 0, 5 });
            p.AddCoeff(5.5588089569691200e+03, new int[] { 0, 0, 7 });
            p.AddCoeff(-9.7279156746959700e+03, new int[] { 0, 0, 9 });
            p.AddCoeff(8.1360749279275500e+03, new int[] { 0, 0, 11 });
            p.AddCoeff(-2.6077163230536900e+03, new int[] { 0, 0, 13 });
            p.AddCoeff(1.8070431593862500e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-5.4211294781587400e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(4.6079600564349200e+03, new int[] { 0, 2, 5 });
            p.AddCoeff(-1.6676426870907300e+04, new int[] { 0, 2, 7 });
            p.AddCoeff(2.9183747024087800e+04, new int[] { 0, 2, 9 });
            p.AddCoeff(-2.4408224783782500e+04, new int[] { 0, 2, 11 });
            p.AddCoeff(7.8231489691611000e+03, new int[] { 0, 2, 13 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f3a4da64-c655-41d2-adca-eeb8786c3a94"));
            OrthonormalPolynomials[683] = p;
            p.AddCoeff(-1.5826224176235000e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(1.2344454857463300e+02, new int[] { 0, 1, 2 });
            p.AddCoeff(-1.5430568571829100e+03, new int[] { 0, 1, 4 });
            p.AddCoeff(6.9951910858958600e+03, new int[] { 0, 1, 6 });
            p.AddCoeff(-1.4240210424859400e+04, new int[] { 0, 1, 8 });
            p.AddCoeff(1.3290863063202200e+04, new int[] { 0, 1, 10 });
            p.AddCoeff(-4.6316644008128600e+03, new int[] { 0, 1, 12 });
            p.AddCoeff(2.6377040293725000e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(-2.0574091429105500e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(2.5717614286381800e+03, new int[] { 0, 3, 4 });
            p.AddCoeff(-1.1658651809826400e+04, new int[] { 0, 3, 6 });
            p.AddCoeff(2.3733684041432400e+04, new int[] { 0, 3, 8 });
            p.AddCoeff(-2.2151438438670200e+04, new int[] { 0, 3, 10 });
            p.AddCoeff(7.7194406680214600e+03, new int[] { 0, 3, 12 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("984dd415-3dd1-4375-b78a-d6c4f2fd7729"));
            OrthonormalPolynomials[684] = p;
            p.AddCoeff(-5.1637441534121500e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(1.1188112332393000e+02, new int[] { 0, 0, 3 });
            p.AddCoeff(-6.7128673994357900e+02, new int[] { 0, 0, 5 });
            p.AddCoeff(1.6302677970058400e+03, new int[] { 0, 0, 7 });
            p.AddCoeff(-1.7208382301728200e+03, new int[] { 0, 0, 9 });
            p.AddCoeff(6.5704732424780600e+02, new int[] { 0, 0, 11 });
            p.AddCoeff(5.1637441534121500e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-1.1188112332393000e+03, new int[] { 0, 2, 3 });
            p.AddCoeff(6.7128673994357900e+03, new int[] { 0, 2, 5 });
            p.AddCoeff(-1.6302677970058400e+04, new int[] { 0, 2, 7 });
            p.AddCoeff(1.7208382301728200e+04, new int[] { 0, 2, 9 });
            p.AddCoeff(-6.5704732424780600e+03, new int[] { 0, 2, 11 });
            p.AddCoeff(-6.0243681789808400e+01, new int[] { 0, 4, 1 });
            p.AddCoeff(1.3052797721125200e+03, new int[] { 0, 4, 3 });
            p.AddCoeff(-7.8316786326750600e+03, new int[] { 0, 4, 5 });
            p.AddCoeff(1.9019790965068100e+04, new int[] { 0, 4, 7 });
            p.AddCoeff(-2.0076446018683000e+04, new int[] { 0, 4, 9 });
            p.AddCoeff(7.6655521162244200e+03, new int[] { 0, 4, 11 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("4fb4ae91-96d0-4917-8ee6-f450f007087b"));
            OrthonormalPolynomials[685] = p;
            p.AddCoeff(-2.4794928065055500e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(1.3637210435780500e+02, new int[] { 0, 1, 2 });
            p.AddCoeff(-1.1818915711009800e+03, new int[] { 0, 1, 4 });
            p.AddCoeff(3.5456747133029300e+03, new int[] { 0, 1, 6 });
            p.AddCoeff(-4.3054621518678400e+03, new int[] { 0, 1, 8 });
            p.AddCoeff(1.8178617974553100e+03, new int[] { 0, 1, 10 });
            p.AddCoeff(1.1570966430359200e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-6.3640315366975700e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(5.5154939984712300e+03, new int[] { 0, 3, 4 });
            p.AddCoeff(-1.6546481995413700e+04, new int[] { 0, 3, 6 });
            p.AddCoeff(2.0092156708716600e+04, new int[] { 0, 3, 8 });
            p.AddCoeff(-8.4833550547914600e+03, new int[] { 0, 3, 10 });
            p.AddCoeff(-1.0413869787323300e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(5.7276283830278100e+02, new int[] { 0, 5, 2 });
            p.AddCoeff(-4.9639445986241100e+03, new int[] { 0, 5, 4 });
            p.AddCoeff(1.4891833795872300e+04, new int[] { 0, 5, 6 });
            p.AddCoeff(-1.8082941037844900e+04, new int[] { 0, 5, 8 });
            p.AddCoeff(7.6350195493123200e+03, new int[] { 0, 5, 10 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("6ac12509-7054-4f1d-b6ac-aa4696366a69"));
            OrthonormalPolynomials[686] = p;
            p.AddCoeff(-4.2732085527534600e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(6.2673725440384000e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-2.4442752921749800e+02, new int[] { 0, 0, 5 });
            p.AddCoeff(3.4918218459642500e+02, new int[] { 0, 0, 7 });
            p.AddCoeff(-1.6489158717053400e+02, new int[] { 0, 0, 9 });
            p.AddCoeff(8.9737379607822600e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-1.3161482342480600e+03, new int[] { 0, 2, 3 });
            p.AddCoeff(5.1329781135674500e+03, new int[] { 0, 2, 5 });
            p.AddCoeff(-7.3328258765249300e+03, new int[] { 0, 2, 7 });
            p.AddCoeff(3.4627233305812300e+03, new int[] { 0, 2, 9 });
            p.AddCoeff(-2.6921213882346800e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(3.9484447027441900e+03, new int[] { 0, 4, 3 });
            p.AddCoeff(-1.5398934340702400e+04, new int[] { 0, 4, 5 });
            p.AddCoeff(2.1998477629574800e+04, new int[] { 0, 4, 7 });
            p.AddCoeff(-1.0388169991743600e+04, new int[] { 0, 4, 9 });
            p.AddCoeff(1.9742223513721000e+02, new int[] { 0, 6, 1 });
            p.AddCoeff(-2.8955261153457400e+03, new int[] { 0, 6, 3 });
            p.AddCoeff(1.1292551849848400e+04, new int[] { 0, 6, 5 });
            p.AddCoeff(-1.6132216928354800e+04, new int[] { 0, 6, 7 });
            p.AddCoeff(7.6179913272786700e+03, new int[] { 0, 6, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("1ea22d8d-fdfa-4c9a-88a5-59af1a16bf0d"));
            OrthonormalPolynomials[687] = p;
            p.AddCoeff(-3.3770013411936900e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(1.2157204828297300e+02, new int[] { 0, 1, 2 });
            p.AddCoeff(-6.6864626555635000e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(1.1589868602976700e+03, new int[] { 0, 1, 6 });
            p.AddCoeff(-6.2088581801661100e+02, new int[] { 0, 1, 8 });
            p.AddCoeff(3.0393012070743200e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-1.0941484345467500e+03, new int[] { 0, 3, 2 });
            p.AddCoeff(6.0178163900071500e+03, new int[] { 0, 3, 4 });
            p.AddCoeff(-1.0430881742679100e+04, new int[] { 0, 3, 6 });
            p.AddCoeff(5.5879723621494900e+03, new int[] { 0, 3, 8 });
            p.AddCoeff(-6.6864626555635000e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(2.4071265560028600e+03, new int[] { 0, 5, 2 });
            p.AddCoeff(-1.3239196058015700e+04, new int[] { 0, 5, 4 });
            p.AddCoeff(2.2947939833893900e+04, new int[] { 0, 5, 6 });
            p.AddCoeff(-1.2293539196728900e+04, new int[] { 0, 5, 8 });
            p.AddCoeff(4.1392387867774100e+01, new int[] { 0, 7, 0 });
            p.AddCoeff(-1.4901259632398700e+03, new int[] { 0, 7, 2 });
            p.AddCoeff(8.1956927978192600e+03, new int[] { 0, 7, 4 });
            p.AddCoeff(-1.4205867516220100e+04, new int[] { 0, 7, 6 });
            p.AddCoeff(7.6102861694036100e+03, new int[] { 0, 7, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("9889a767-81da-4c86-8c57-c60e3c3f9343"));
            OrthonormalPolynomials[688] = p;
            p.AddCoeff(-3.3770013411936900e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(3.0393012070743200e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-6.6864626555635000e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(4.1392387867774100e+01, new int[] { 0, 0, 7 });
            p.AddCoeff(1.2157204828297300e+02, new int[] { 0, 2, 1 });
            p.AddCoeff(-1.0941484345467500e+03, new int[] { 0, 2, 3 });
            p.AddCoeff(2.4071265560028600e+03, new int[] { 0, 2, 5 });
            p.AddCoeff(-1.4901259632398700e+03, new int[] { 0, 2, 7 });
            p.AddCoeff(-6.6864626555635000e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(6.0178163900071500e+03, new int[] { 0, 4, 3 });
            p.AddCoeff(-1.3239196058015700e+04, new int[] { 0, 4, 5 });
            p.AddCoeff(8.1956927978192600e+03, new int[] { 0, 4, 7 });
            p.AddCoeff(1.1589868602976700e+03, new int[] { 0, 6, 1 });
            p.AddCoeff(-1.0430881742679100e+04, new int[] { 0, 6, 3 });
            p.AddCoeff(2.2947939833893900e+04, new int[] { 0, 6, 5 });
            p.AddCoeff(-1.4205867516220100e+04, new int[] { 0, 6, 7 });
            p.AddCoeff(-6.2088581801661100e+02, new int[] { 0, 8, 1 });
            p.AddCoeff(5.5879723621494900e+03, new int[] { 0, 8, 3 });
            p.AddCoeff(-1.2293539196728900e+04, new int[] { 0, 8, 5 });
            p.AddCoeff(7.6102861694036100e+03, new int[] { 0, 8, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c1e4a442-5843-418e-a85b-f2796c0a88ba"));
            OrthonormalPolynomials[689] = p;
            p.AddCoeff(-4.2732085527534600e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(8.9737379607822600e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-2.6921213882346800e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(1.9742223513721000e+02, new int[] { 0, 1, 6 });
            p.AddCoeff(6.2673725440384000e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-1.3161482342480600e+03, new int[] { 0, 3, 2 });
            p.AddCoeff(3.9484447027441900e+03, new int[] { 0, 3, 4 });
            p.AddCoeff(-2.8955261153457400e+03, new int[] { 0, 3, 6 });
            p.AddCoeff(-2.4442752921749800e+02, new int[] { 0, 5, 0 });
            p.AddCoeff(5.1329781135674500e+03, new int[] { 0, 5, 2 });
            p.AddCoeff(-1.5398934340702400e+04, new int[] { 0, 5, 4 });
            p.AddCoeff(1.1292551849848400e+04, new int[] { 0, 5, 6 });
            p.AddCoeff(3.4918218459642500e+02, new int[] { 0, 7, 0 });
            p.AddCoeff(-7.3328258765249300e+03, new int[] { 0, 7, 2 });
            p.AddCoeff(2.1998477629574800e+04, new int[] { 0, 7, 4 });
            p.AddCoeff(-1.6132216928354800e+04, new int[] { 0, 7, 6 });
            p.AddCoeff(-1.6489158717053400e+02, new int[] { 0, 9, 0 });
            p.AddCoeff(3.4627233305812300e+03, new int[] { 0, 9, 2 });
            p.AddCoeff(-1.0388169991743600e+04, new int[] { 0, 9, 4 });
            p.AddCoeff(7.6179913272786700e+03, new int[] { 0, 9, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a4b4c4f3-ce7d-4111-b0c8-355fe6c9aa57"));
            OrthonormalPolynomials[690] = p;
            p.AddCoeff(-2.4794928065055500e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(1.1570966430359200e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-1.0413869787323300e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(1.3637210435780500e+02, new int[] { 0, 2, 1 });
            p.AddCoeff(-6.3640315366975700e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(5.7276283830278100e+02, new int[] { 0, 2, 5 });
            p.AddCoeff(-1.1818915711009800e+03, new int[] { 0, 4, 1 });
            p.AddCoeff(5.5154939984712300e+03, new int[] { 0, 4, 3 });
            p.AddCoeff(-4.9639445986241100e+03, new int[] { 0, 4, 5 });
            p.AddCoeff(3.5456747133029300e+03, new int[] { 0, 6, 1 });
            p.AddCoeff(-1.6546481995413700e+04, new int[] { 0, 6, 3 });
            p.AddCoeff(1.4891833795872300e+04, new int[] { 0, 6, 5 });
            p.AddCoeff(-4.3054621518678400e+03, new int[] { 0, 8, 1 });
            p.AddCoeff(2.0092156708716600e+04, new int[] { 0, 8, 3 });
            p.AddCoeff(-1.8082941037844900e+04, new int[] { 0, 8, 5 });
            p.AddCoeff(1.8178617974553100e+03, new int[] { 0, 10, 1 });
            p.AddCoeff(-8.4833550547914600e+03, new int[] { 0, 10, 3 });
            p.AddCoeff(7.6350195493123200e+03, new int[] { 0, 10, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("579f0279-0700-410f-808e-34d9365d215d"));
            OrthonormalPolynomials[691] = p;
            p.AddCoeff(-5.1637441534121500e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(5.1637441534121500e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-6.0243681789808400e+01, new int[] { 0, 1, 4 });
            p.AddCoeff(1.1188112332393000e+02, new int[] { 0, 3, 0 });
            p.AddCoeff(-1.1188112332393000e+03, new int[] { 0, 3, 2 });
            p.AddCoeff(1.3052797721125200e+03, new int[] { 0, 3, 4 });
            p.AddCoeff(-6.7128673994357900e+02, new int[] { 0, 5, 0 });
            p.AddCoeff(6.7128673994357900e+03, new int[] { 0, 5, 2 });
            p.AddCoeff(-7.8316786326750600e+03, new int[] { 0, 5, 4 });
            p.AddCoeff(1.6302677970058400e+03, new int[] { 0, 7, 0 });
            p.AddCoeff(-1.6302677970058400e+04, new int[] { 0, 7, 2 });
            p.AddCoeff(1.9019790965068100e+04, new int[] { 0, 7, 4 });
            p.AddCoeff(-1.7208382301728200e+03, new int[] { 0, 9, 0 });
            p.AddCoeff(1.7208382301728200e+04, new int[] { 0, 9, 2 });
            p.AddCoeff(-2.0076446018683000e+04, new int[] { 0, 9, 4 });
            p.AddCoeff(6.5704732424780600e+02, new int[] { 0, 11, 0 });
            p.AddCoeff(-6.5704732424780600e+03, new int[] { 0, 11, 2 });
            p.AddCoeff(7.6655521162244200e+03, new int[] { 0, 11, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("8c6b45fa-ae48-47ba-8ffc-0d5cdbef565a"));
            OrthonormalPolynomials[692] = p;
            p.AddCoeff(-1.5826224176235000e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(2.6377040293725000e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(1.2344454857463300e+02, new int[] { 0, 2, 1 });
            p.AddCoeff(-2.0574091429105500e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(-1.5430568571829100e+03, new int[] { 0, 4, 1 });
            p.AddCoeff(2.5717614286381800e+03, new int[] { 0, 4, 3 });
            p.AddCoeff(6.9951910858958600e+03, new int[] { 0, 6, 1 });
            p.AddCoeff(-1.1658651809826400e+04, new int[] { 0, 6, 3 });
            p.AddCoeff(-1.4240210424859400e+04, new int[] { 0, 8, 1 });
            p.AddCoeff(2.3733684041432400e+04, new int[] { 0, 8, 3 });
            p.AddCoeff(1.3290863063202200e+04, new int[] { 0, 10, 1 });
            p.AddCoeff(-2.2151438438670200e+04, new int[] { 0, 10, 3 });
            p.AddCoeff(-4.6316644008128600e+03, new int[] { 0, 12, 1 });
            p.AddCoeff(7.7194406680214600e+03, new int[] { 0, 12, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3a112b80-34bc-4eda-998e-5270a29e2a52"));
            OrthonormalPolynomials[693] = p;
            p.AddCoeff(-6.0234771979541500e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(1.8070431593862500e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(1.8070431593862500e+02, new int[] { 0, 3, 0 });
            p.AddCoeff(-5.4211294781587400e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(-1.5359866854783100e+03, new int[] { 0, 5, 0 });
            p.AddCoeff(4.6079600564349200e+03, new int[] { 0, 5, 2 });
            p.AddCoeff(5.5588089569691200e+03, new int[] { 0, 7, 0 });
            p.AddCoeff(-1.6676426870907300e+04, new int[] { 0, 7, 2 });
            p.AddCoeff(-9.7279156746959700e+03, new int[] { 0, 9, 0 });
            p.AddCoeff(2.9183747024087800e+04, new int[] { 0, 9, 2 });
            p.AddCoeff(8.1360749279275500e+03, new int[] { 0, 11, 0 });
            p.AddCoeff(-2.4408224783782500e+04, new int[] { 0, 11, 2 });
            p.AddCoeff(-2.6077163230536900e+03, new int[] { 0, 13, 0 });
            p.AddCoeff(7.8231489691611000e+03, new int[] { 0, 13, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c2fd0e62-2b5d-43f3-8459-110c235c833a"));
            OrthonormalPolynomials[694] = p;
            p.AddCoeff(-6.9078352735584400e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(7.2532270372363600e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-1.2330485963301800e+03, new int[] { 0, 4, 1 });
            p.AddCoeff(7.8093077767578100e+03, new int[] { 0, 6, 1 });
            p.AddCoeff(-2.3427923330273500e+04, new int[] { 0, 8, 1 });
            p.AddCoeff(3.5922815773086000e+04, new int[] { 0, 10, 1 });
            p.AddCoeff(-2.7214254373550000e+04, new int[] { 0, 12, 1 });
            p.AddCoeff(8.0745589899543900e+03, new int[] { 0, 14, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5699fbd0-dbfe-422a-b6b2-da4791ba0529"));
            OrthonormalPolynomials[695] = p;
            p.AddCoeff(-6.1852100426350100e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(2.4534666502452200e+02, new int[] { 0, 3, 0 });
            p.AddCoeff(-2.7969519812795600e+03, new int[] { 0, 5, 0 });
            p.AddCoeff(1.3984759906397800e+04, new int[] { 0, 7, 0 });
            p.AddCoeff(-3.5738830871905400e+04, new int[] { 0, 9, 0 });
            p.AddCoeff(4.8734769370780000e+04, new int[] { 0, 11, 0 });
            p.AddCoeff(-3.3739455718232300e+04, new int[] { 0, 13, 0 });
            p.AddCoeff(9.3185163412260600e+03, new int[] { 0, 15, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("bb4b13a7-ddf0-429d-ad0f-02dae73cfc73"));
            OrthonormalPolynomials[696] = p;
            p.AddCoeff(-6.9078352735584400e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(7.2532270372363600e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-1.2330485963301800e+03, new int[] { 1, 0, 4 });
            p.AddCoeff(7.8093077767578100e+03, new int[] { 1, 0, 6 });
            p.AddCoeff(-2.3427923330273500e+04, new int[] { 1, 0, 8 });
            p.AddCoeff(3.5922815773086000e+04, new int[] { 1, 0, 10 });
            p.AddCoeff(-2.7214254373550000e+04, new int[] { 1, 0, 12 });
            p.AddCoeff(8.0745589899543900e+03, new int[] { 1, 0, 14 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("6b10e409-cc79-4b62-aca4-2023480dad92"));
            OrthonormalPolynomials[697] = p;
            p.AddCoeff(1.6162685370654500e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-4.8488056111963600e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(4.1214847695169000e+03, new int[] { 1, 1, 5 });
            p.AddCoeff(-1.4915849642061200e+04, new int[] { 1, 1, 7 });
            p.AddCoeff(2.6102736873607000e+04, new int[] { 1, 1, 9 });
            p.AddCoeff(-2.1831379930653200e+04, new int[] { 1, 1, 11 });
            p.AddCoeff(6.9972371572606500e+03, new int[] { 1, 1, 13 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("06b6873e-f16a-4b6b-958f-48f1de1f033a"));
            OrthonormalPolynomials[698] = p;
            p.AddCoeff(-7.7224066640437800e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(6.0234771979541500e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-7.5293464974426600e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(3.4133037455073500e+03, new int[] { 1, 0, 6 });
            p.AddCoeff(-6.9485111962113800e+03, new int[] { 1, 0, 8 });
            p.AddCoeff(6.4852771164639400e+03, new int[] { 1, 0, 10 });
            p.AddCoeff(-2.2600208133132000e+03, new int[] { 1, 0, 12 });
            p.AddCoeff(2.3167219992131300e+00, new int[] { 1, 2, 0 });
            p.AddCoeff(-1.8070431593862500e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(2.2588039492328100e+03, new int[] { 1, 2, 4 });
            p.AddCoeff(-1.0239911236522100e+04, new int[] { 1, 2, 6 });
            p.AddCoeff(2.0845533588634200e+04, new int[] { 1, 2, 8 });
            p.AddCoeff(-1.9455831349391900e+04, new int[] { 1, 2, 10 });
            p.AddCoeff(6.7800624399395900e+03, new int[] { 1, 2, 12 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("fead6b44-767f-4af8-ae1f-efc526d1e8c3"));
            OrthonormalPolynomials[699] = p;
            p.AddCoeff(3.1550997936529100e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-6.8360495529146400e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(4.1016297317487800e+03, new int[] { 1, 1, 5 });
            p.AddCoeff(-9.9611007771041800e+03, new int[] { 1, 1, 7 });
            p.AddCoeff(1.0514495264721100e+04, new int[] { 1, 1, 9 });
            p.AddCoeff(-4.0146254647116900e+03, new int[] { 1, 1, 11 });
            p.AddCoeff(-5.2584996560881800e+01, new int[] { 1, 3, 1 });
            p.AddCoeff(1.1393415921524400e+03, new int[] { 1, 3, 3 });
            p.AddCoeff(-6.8360495529146400e+03, new int[] { 1, 3, 5 });
            p.AddCoeff(1.6601834628507000e+04, new int[] { 1, 3, 7 });
            p.AddCoeff(-1.7524158774535100e+04, new int[] { 1, 3, 9 });
            p.AddCoeff(6.6910424411861400e+03, new int[] { 1, 3, 11 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3b287074-a740-4179-b6f9-39654fd3a3c5"));
            OrthonormalPolynomials[700] = p;
            p.AddCoeff(-7.7692373228789900e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(4.2730805275834500e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-3.7033364572389900e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(1.1110009371717000e+03, new int[] { 1, 0, 6 });
            p.AddCoeff(-1.3490725665656300e+03, new int[] { 1, 0, 8 });
            p.AddCoeff(5.6960841699437900e+02, new int[] { 1, 0, 10 });
            p.AddCoeff(7.7692373228789900e+00, new int[] { 1, 2, 0 });
            p.AddCoeff(-4.2730805275834500e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(3.7033364572389900e+03, new int[] { 1, 2, 4 });
            p.AddCoeff(-1.1110009371717000e+04, new int[] { 1, 2, 6 });
            p.AddCoeff(1.3490725665656300e+04, new int[] { 1, 2, 8 });
            p.AddCoeff(-5.6960841699437900e+03, new int[] { 1, 2, 10 });
            p.AddCoeff(-9.0641102100254900e+00, new int[] { 1, 4, 0 });
            p.AddCoeff(4.9852606155140400e+02, new int[] { 1, 4, 2 });
            p.AddCoeff(-4.3205592001121300e+03, new int[] { 1, 4, 4 });
            p.AddCoeff(1.2961677600336500e+04, new int[] { 1, 4, 6 });
            p.AddCoeff(-1.5739179943265700e+04, new int[] { 1, 4, 8 });
            p.AddCoeff(6.6454315316010800e+03, new int[] { 1, 4, 10 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("cc773ec6-683b-4906-91a0-f489f5ded108"));
            OrthonormalPolynomials[701] = p;
            p.AddCoeff(4.0849865705801600e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-5.9913136368509000e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(2.3366123183718500e+03, new int[] { 1, 1, 5 });
            p.AddCoeff(-3.3380175976740700e+03, new int[] { 1, 1, 7 });
            p.AddCoeff(1.5762860877905300e+03, new int[] { 1, 1, 9 });
            p.AddCoeff(-1.9063270662707400e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(2.7959463638637500e+03, new int[] { 1, 3, 3 });
            p.AddCoeff(-1.0904190819068600e+04, new int[] { 1, 3, 5 });
            p.AddCoeff(1.5577415455812300e+04, new int[] { 1, 3, 7 });
            p.AddCoeff(-7.3560017430224900e+03, new int[] { 1, 3, 9 });
            p.AddCoeff(1.7156943596436700e+02, new int[] { 1, 5, 1 });
            p.AddCoeff(-2.5163517274773800e+03, new int[] { 1, 5, 3 });
            p.AddCoeff(9.8137717371617800e+03, new int[] { 1, 5, 5 });
            p.AddCoeff(-1.4019673910231100e+04, new int[] { 1, 5, 7 });
            p.AddCoeff(6.6204015687202500e+03, new int[] { 1, 5, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("52156fcb-ed90-4328-9807-1ee4b8b318b8"));
            OrthonormalPolynomials[702] = p;
            p.AddCoeff(-7.7789300654438300e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(2.8004148235597800e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-1.5402281529578800e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(2.6697287984603200e+02, new int[] { 1, 0, 6 });
            p.AddCoeff(-1.4302118563180300e+02, new int[] { 1, 0, 8 });
            p.AddCoeff(1.6335753137432000e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-5.8808711294755300e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(3.2344791212115400e+03, new int[] { 1, 2, 4 });
            p.AddCoeff(-5.6064304767666800e+03, new int[] { 1, 2, 6 });
            p.AddCoeff(3.0034448982678600e+03, new int[] { 1, 2, 8 });
            p.AddCoeff(-4.9007259412296100e+01, new int[] { 1, 4, 0 });
            p.AddCoeff(1.7642613388426600e+03, new int[] { 1, 4, 2 });
            p.AddCoeff(-9.7034373636346300e+03, new int[] { 1, 4, 4 });
            p.AddCoeff(1.6819291430300000e+04, new int[] { 1, 4, 6 });
            p.AddCoeff(-9.0103346948035700e+03, new int[] { 1, 4, 8 });
            p.AddCoeff(3.5938656902350500e+01, new int[] { 1, 6, 0 });
            p.AddCoeff(-1.2937916484846200e+03, new int[] { 1, 6, 2 });
            p.AddCoeff(7.1158540666654000e+03, new int[] { 1, 6, 4 });
            p.AddCoeff(-1.2334147048886700e+04, new int[] { 1, 6, 6 });
            p.AddCoeff(6.6075787761892800e+03, new int[] { 1, 6, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3ab0a3c3-0166-4875-94e8-b79b9018b776"));
            OrthonormalPolynomials[703] = p;
            p.AddCoeff(4.3954466819961800e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-3.9559020137965700e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(8.7029844303524500e+02, new int[] { 1, 1, 5 });
            p.AddCoeff(-5.3875617902181800e+02, new int[] { 1, 1, 7 });
            p.AddCoeff(-3.9559020137965700e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(3.5603118124169100e+03, new int[] { 1, 3, 3 });
            p.AddCoeff(-7.8326859873172000e+03, new int[] { 1, 3, 5 });
            p.AddCoeff(4.8488056111963600e+03, new int[] { 1, 3, 7 });
            p.AddCoeff(8.7029844303524500e+02, new int[] { 1, 5, 1 });
            p.AddCoeff(-7.8326859873172000e+03, new int[] { 1, 5, 3 });
            p.AddCoeff(1.7231909172097800e+04, new int[] { 1, 5, 5 });
            p.AddCoeff(-1.0667372344632000e+04, new int[] { 1, 5, 7 });
            p.AddCoeff(-5.3875617902181800e+02, new int[] { 1, 7, 1 });
            p.AddCoeff(4.8488056111963600e+03, new int[] { 1, 7, 3 });
            p.AddCoeff(-1.0667372344632000e+04, new int[] { 1, 7, 5 });
            p.AddCoeff(6.6036114514388600e+03, new int[] { 1, 7, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("ad59a451-52c6-4516-9905-6c1dd7bafd38"));
            OrthonormalPolynomials[704] = p;
            p.AddCoeff(-7.7789300654438300e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(1.6335753137432000e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-4.9007259412296100e+01, new int[] { 1, 0, 4 });
            p.AddCoeff(3.5938656902350500e+01, new int[] { 1, 0, 6 });
            p.AddCoeff(2.8004148235597800e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-5.8808711294755300e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(1.7642613388426600e+03, new int[] { 1, 2, 4 });
            p.AddCoeff(-1.2937916484846200e+03, new int[] { 1, 2, 6 });
            p.AddCoeff(-1.5402281529578800e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(3.2344791212115400e+03, new int[] { 1, 4, 2 });
            p.AddCoeff(-9.7034373636346300e+03, new int[] { 1, 4, 4 });
            p.AddCoeff(7.1158540666654000e+03, new int[] { 1, 4, 6 });
            p.AddCoeff(2.6697287984603200e+02, new int[] { 1, 6, 0 });
            p.AddCoeff(-5.6064304767666800e+03, new int[] { 1, 6, 2 });
            p.AddCoeff(1.6819291430300000e+04, new int[] { 1, 6, 4 });
            p.AddCoeff(-1.2334147048886700e+04, new int[] { 1, 6, 6 });
            p.AddCoeff(-1.4302118563180300e+02, new int[] { 1, 8, 0 });
            p.AddCoeff(3.0034448982678600e+03, new int[] { 1, 8, 2 });
            p.AddCoeff(-9.0103346948035700e+03, new int[] { 1, 8, 4 });
            p.AddCoeff(6.6075787761892800e+03, new int[] { 1, 8, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("eef6e21c-9d8e-40a5-a9d5-300ec115d11e"));
            OrthonormalPolynomials[705] = p;
            p.AddCoeff(4.0849865705801600e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-1.9063270662707400e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(1.7156943596436700e+02, new int[] { 1, 1, 5 });
            p.AddCoeff(-5.9913136368509000e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(2.7959463638637500e+03, new int[] { 1, 3, 3 });
            p.AddCoeff(-2.5163517274773800e+03, new int[] { 1, 3, 5 });
            p.AddCoeff(2.3366123183718500e+03, new int[] { 1, 5, 1 });
            p.AddCoeff(-1.0904190819068600e+04, new int[] { 1, 5, 3 });
            p.AddCoeff(9.8137717371617800e+03, new int[] { 1, 5, 5 });
            p.AddCoeff(-3.3380175976740700e+03, new int[] { 1, 7, 1 });
            p.AddCoeff(1.5577415455812300e+04, new int[] { 1, 7, 3 });
            p.AddCoeff(-1.4019673910231100e+04, new int[] { 1, 7, 5 });
            p.AddCoeff(1.5762860877905300e+03, new int[] { 1, 9, 1 });
            p.AddCoeff(-7.3560017430224900e+03, new int[] { 1, 9, 3 });
            p.AddCoeff(6.6204015687202500e+03, new int[] { 1, 9, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("44f23c9c-6406-4407-b1ed-c78169d8d385"));
            OrthonormalPolynomials[706] = p;
            p.AddCoeff(-7.7692373228789900e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(7.7692373228789900e+00, new int[] { 1, 0, 2 });
            p.AddCoeff(-9.0641102100254900e+00, new int[] { 1, 0, 4 });
            p.AddCoeff(4.2730805275834500e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-4.2730805275834500e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(4.9852606155140400e+02, new int[] { 1, 2, 4 });
            p.AddCoeff(-3.7033364572389900e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(3.7033364572389900e+03, new int[] { 1, 4, 2 });
            p.AddCoeff(-4.3205592001121300e+03, new int[] { 1, 4, 4 });
            p.AddCoeff(1.1110009371717000e+03, new int[] { 1, 6, 0 });
            p.AddCoeff(-1.1110009371717000e+04, new int[] { 1, 6, 2 });
            p.AddCoeff(1.2961677600336500e+04, new int[] { 1, 6, 4 });
            p.AddCoeff(-1.3490725665656300e+03, new int[] { 1, 8, 0 });
            p.AddCoeff(1.3490725665656300e+04, new int[] { 1, 8, 2 });
            p.AddCoeff(-1.5739179943265700e+04, new int[] { 1, 8, 4 });
            p.AddCoeff(5.6960841699437900e+02, new int[] { 1, 10, 0 });
            p.AddCoeff(-5.6960841699437900e+03, new int[] { 1, 10, 2 });
            p.AddCoeff(6.6454315316010800e+03, new int[] { 1, 10, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c7644948-b695-4b8c-adcc-d0d7a6eb6b51"));
            OrthonormalPolynomials[707] = p;
            p.AddCoeff(3.1550997936529100e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-5.2584996560881800e+01, new int[] { 1, 1, 3 });
            p.AddCoeff(-6.8360495529146400e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(1.1393415921524400e+03, new int[] { 1, 3, 3 });
            p.AddCoeff(4.1016297317487800e+03, new int[] { 1, 5, 1 });
            p.AddCoeff(-6.8360495529146400e+03, new int[] { 1, 5, 3 });
            p.AddCoeff(-9.9611007771041800e+03, new int[] { 1, 7, 1 });
            p.AddCoeff(1.6601834628507000e+04, new int[] { 1, 7, 3 });
            p.AddCoeff(1.0514495264721100e+04, new int[] { 1, 9, 1 });
            p.AddCoeff(-1.7524158774535100e+04, new int[] { 1, 9, 3 });
            p.AddCoeff(-4.0146254647116900e+03, new int[] { 1, 11, 1 });
            p.AddCoeff(6.6910424411861400e+03, new int[] { 1, 11, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("1426ef39-582e-4cae-8877-c113228881e2"));
            OrthonormalPolynomials[708] = p;
            p.AddCoeff(-7.7224066640437800e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(2.3167219992131300e+00, new int[] { 1, 0, 2 });
            p.AddCoeff(6.0234771979541500e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-1.8070431593862500e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(-7.5293464974426600e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(2.2588039492328100e+03, new int[] { 1, 4, 2 });
            p.AddCoeff(3.4133037455073500e+03, new int[] { 1, 6, 0 });
            p.AddCoeff(-1.0239911236522100e+04, new int[] { 1, 6, 2 });
            p.AddCoeff(-6.9485111962113800e+03, new int[] { 1, 8, 0 });
            p.AddCoeff(2.0845533588634200e+04, new int[] { 1, 8, 2 });
            p.AddCoeff(6.4852771164639400e+03, new int[] { 1, 10, 0 });
            p.AddCoeff(-1.9455831349391900e+04, new int[] { 1, 10, 2 });
            p.AddCoeff(-2.2600208133132000e+03, new int[] { 1, 12, 0 });
            p.AddCoeff(6.7800624399395900e+03, new int[] { 1, 12, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("63930613-d429-4df2-aa82-b968d0073cb2"));
            OrthonormalPolynomials[709] = p;
            p.AddCoeff(1.6162685370654500e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-4.8488056111963600e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(4.1214847695169000e+03, new int[] { 1, 5, 1 });
            p.AddCoeff(-1.4915849642061200e+04, new int[] { 1, 7, 1 });
            p.AddCoeff(2.6102736873607000e+04, new int[] { 1, 9, 1 });
            p.AddCoeff(-2.1831379930653200e+04, new int[] { 1, 11, 1 });
            p.AddCoeff(6.9972371572606500e+03, new int[] { 1, 13, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("92b72f31-d72d-4cd8-b9a1-4659f8d3213b"));
            OrthonormalPolynomials[710] = p;
            p.AddCoeff(-6.9078352735584400e-01, new int[] { 1, 0, 0 });
            p.AddCoeff(7.2532270372363600e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-1.2330485963301800e+03, new int[] { 1, 4, 0 });
            p.AddCoeff(7.8093077767578100e+03, new int[] { 1, 6, 0 });
            p.AddCoeff(-2.3427923330273500e+04, new int[] { 1, 8, 0 });
            p.AddCoeff(3.5922815773086000e+04, new int[] { 1, 10, 0 });
            p.AddCoeff(-2.7214254373550000e+04, new int[] { 1, 12, 0 });
            p.AddCoeff(8.0745589899543900e+03, new int[] { 1, 14, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("9bb5b475-c7b6-43e7-830f-a3b12344a005"));
            OrthonormalPolynomials[711] = p;
            p.AddCoeff(-6.0234771979541500e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(1.8070431593862500e+02, new int[] { 0, 0, 3 });
            p.AddCoeff(-1.5359866854783100e+03, new int[] { 0, 0, 5 });
            p.AddCoeff(5.5588089569691200e+03, new int[] { 0, 0, 7 });
            p.AddCoeff(-9.7279156746959700e+03, new int[] { 0, 0, 9 });
            p.AddCoeff(8.1360749279275500e+03, new int[] { 0, 0, 11 });
            p.AddCoeff(-2.6077163230536900e+03, new int[] { 0, 0, 13 });
            p.AddCoeff(1.8070431593862500e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-5.4211294781587400e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(4.6079600564349200e+03, new int[] { 2, 0, 5 });
            p.AddCoeff(-1.6676426870907300e+04, new int[] { 2, 0, 7 });
            p.AddCoeff(2.9183747024087800e+04, new int[] { 2, 0, 9 });
            p.AddCoeff(-2.4408224783782500e+04, new int[] { 2, 0, 11 });
            p.AddCoeff(7.8231489691611000e+03, new int[] { 2, 0, 13 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("75ab7330-9d5a-484c-b500-e81394e8e68c"));
            OrthonormalPolynomials[712] = p;
            p.AddCoeff(-7.7224066640437800e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(6.0234771979541500e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-7.5293464974426600e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(3.4133037455073500e+03, new int[] { 0, 1, 6 });
            p.AddCoeff(-6.9485111962113800e+03, new int[] { 0, 1, 8 });
            p.AddCoeff(6.4852771164639400e+03, new int[] { 0, 1, 10 });
            p.AddCoeff(-2.2600208133132000e+03, new int[] { 0, 1, 12 });
            p.AddCoeff(2.3167219992131300e+00, new int[] { 2, 1, 0 });
            p.AddCoeff(-1.8070431593862500e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(2.2588039492328100e+03, new int[] { 2, 1, 4 });
            p.AddCoeff(-1.0239911236522100e+04, new int[] { 2, 1, 6 });
            p.AddCoeff(2.0845533588634200e+04, new int[] { 2, 1, 8 });
            p.AddCoeff(-1.9455831349391900e+04, new int[] { 2, 1, 10 });
            p.AddCoeff(6.7800624399395900e+03, new int[] { 2, 1, 12 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("eaa86998-1ae7-4520-a820-7484b1aa4ba0"));
            OrthonormalPolynomials[713] = p;
            p.AddCoeff(-5.7374935037912700e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(1.2431235924881100e+02, new int[] { 0, 0, 3 });
            p.AddCoeff(-7.4587415549286600e+02, new int[] { 0, 0, 5 });
            p.AddCoeff(1.8114086633398200e+03, new int[] { 0, 0, 7 });
            p.AddCoeff(-1.9120424779698100e+03, new int[] { 0, 0, 9 });
            p.AddCoeff(7.3005258249756200e+02, new int[] { 0, 0, 11 });
            p.AddCoeff(1.7212480511373800e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-3.7293707774643300e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(2.2376224664786000e+03, new int[] { 0, 2, 5 });
            p.AddCoeff(-5.4342259900194500e+03, new int[] { 0, 2, 7 });
            p.AddCoeff(5.7361274339094200e+03, new int[] { 0, 2, 9 });
            p.AddCoeff(-2.1901577474926900e+03, new int[] { 0, 2, 11 });
            p.AddCoeff(1.7212480511373800e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-3.7293707774643300e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(2.2376224664786000e+03, new int[] { 2, 0, 5 });
            p.AddCoeff(-5.4342259900194500e+03, new int[] { 2, 0, 7 });
            p.AddCoeff(5.7361274339094200e+03, new int[] { 2, 0, 9 });
            p.AddCoeff(-2.1901577474926900e+03, new int[] { 2, 0, 11 });
            p.AddCoeff(-5.1637441534121500e+01, new int[] { 2, 2, 1 });
            p.AddCoeff(1.1188112332393000e+03, new int[] { 2, 2, 3 });
            p.AddCoeff(-6.7128673994357900e+03, new int[] { 2, 2, 5 });
            p.AddCoeff(1.6302677970058400e+04, new int[] { 2, 2, 7 });
            p.AddCoeff(-1.7208382301728200e+04, new int[] { 2, 2, 9 });
            p.AddCoeff(6.5704732424780600e+03, new int[] { 2, 2, 11 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("14101eee-e089-4ee7-9929-ccc2827dc555"));
            OrthonormalPolynomials[714] = p;
            p.AddCoeff(-1.7691331630354800e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(9.7302323966951700e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-8.4328680771358100e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(2.5298604231407400e+03, new int[] { 0, 1, 6 });
            p.AddCoeff(-3.0719733709566200e+03, new int[] { 0, 1, 8 });
            p.AddCoeff(1.2970554232927900e+03, new int[] { 0, 1, 10 });
            p.AddCoeff(2.9485552717258100e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(-1.6217053994491900e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(1.4054780128559700e+03, new int[] { 0, 3, 4 });
            p.AddCoeff(-4.2164340385679100e+03, new int[] { 0, 3, 6 });
            p.AddCoeff(5.1199556182610300e+03, new int[] { 0, 3, 8 });
            p.AddCoeff(-2.1617590388213200e+03, new int[] { 0, 3, 10 });
            p.AddCoeff(5.3073994891064500e+00, new int[] { 2, 1, 0 });
            p.AddCoeff(-2.9190697190085500e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(2.5298604231407400e+03, new int[] { 2, 1, 4 });
            p.AddCoeff(-7.5895812694222300e+03, new int[] { 2, 1, 6 });
            p.AddCoeff(9.2159201128698200e+03, new int[] { 2, 1, 8 });
            p.AddCoeff(-3.8911662698783800e+03, new int[] { 2, 1, 10 });
            p.AddCoeff(-8.8456658151774200e+00, new int[] { 2, 3, 0 });
            p.AddCoeff(4.8651161983475800e+02, new int[] { 2, 3, 2 });
            p.AddCoeff(-4.2164340385679100e+03, new int[] { 2, 3, 4 });
            p.AddCoeff(1.2649302115703700e+04, new int[] { 2, 3, 6 });
            p.AddCoeff(-1.5359866854783100e+04, new int[] { 2, 3, 8 });
            p.AddCoeff(6.4852771164639400e+03, new int[] { 2, 3, 10 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e50d648d-602d-4494-801e-21d906b0d00e"));
            OrthonormalPolynomials[715] = p;
            p.AddCoeff(-4.7702365981659800e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(6.9963470106434400e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-2.7285753341509400e+02, new int[] { 0, 0, 5 });
            p.AddCoeff(3.8979647630727700e+02, new int[] { 0, 0, 7 });
            p.AddCoeff(-1.8407055825621400e+02, new int[] { 0, 0, 9 });
            p.AddCoeff(4.7702365981659800e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-6.9963470106434400e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(2.7285753341509400e+03, new int[] { 0, 2, 5 });
            p.AddCoeff(-3.8979647630727700e+03, new int[] { 0, 2, 7 });
            p.AddCoeff(1.8407055825621400e+03, new int[] { 0, 2, 9 });
            p.AddCoeff(-5.5652760311936400e+01, new int[] { 0, 4, 1 });
            p.AddCoeff(8.1624048457506800e+02, new int[] { 0, 4, 3 });
            p.AddCoeff(-3.1833378898427600e+03, new int[] { 0, 4, 5 });
            p.AddCoeff(4.5476255569182400e+03, new int[] { 0, 4, 7 });
            p.AddCoeff(-2.1474898463225100e+03, new int[] { 0, 4, 9 });
            p.AddCoeff(1.4310709794497900e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-2.0989041031930300e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(8.1857260024528200e+02, new int[] { 2, 0, 5 });
            p.AddCoeff(-1.1693894289218300e+03, new int[] { 2, 0, 7 });
            p.AddCoeff(5.5221167476864300e+02, new int[] { 2, 0, 9 });
            p.AddCoeff(-1.4310709794497900e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(2.0989041031930300e+03, new int[] { 2, 2, 3 });
            p.AddCoeff(-8.1857260024528200e+03, new int[] { 2, 2, 5 });
            p.AddCoeff(1.1693894289218300e+04, new int[] { 2, 2, 7 });
            p.AddCoeff(-5.5221167476864300e+03, new int[] { 2, 2, 9 });
            p.AddCoeff(1.6695828093580900e+02, new int[] { 2, 4, 1 });
            p.AddCoeff(-2.4487214537252000e+03, new int[] { 2, 4, 3 });
            p.AddCoeff(9.5500136695282900e+03, new int[] { 2, 4, 5 });
            p.AddCoeff(-1.3642876670754700e+04, new int[] { 2, 4, 7 });
            p.AddCoeff(6.4424695389674900e+03, new int[] { 2, 4, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b38a6bb4-590f-46d4-9d37-871fd71714f2"));
            OrthonormalPolynomials[716] = p;
            p.AddCoeff(-2.7713422517043000e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(9.9768321061354800e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-5.4872576583745200e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(9.5112466078491600e+02, new int[] { 0, 1, 6 });
            p.AddCoeff(-5.0953106827763400e+02, new int[] { 0, 1, 8 });
            p.AddCoeff(1.2932930507953400e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-4.6558549828632300e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(2.5607202405747700e+03, new int[] { 0, 3, 4 });
            p.AddCoeff(-4.4385817503296100e+03, new int[] { 0, 3, 6 });
            p.AddCoeff(2.3778116519622900e+03, new int[] { 0, 3, 8 });
            p.AddCoeff(-1.1639637457158100e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(4.1902694845769000e+02, new int[] { 0, 5, 2 });
            p.AddCoeff(-2.3046482165173000e+03, new int[] { 0, 5, 4 });
            p.AddCoeff(3.9947235752966500e+03, new int[] { 0, 5, 6 });
            p.AddCoeff(-2.1400304867660600e+03, new int[] { 0, 5, 8 });
            p.AddCoeff(8.3140267551129000e+00, new int[] { 2, 1, 0 });
            p.AddCoeff(-2.9930496318406400e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(1.6461772975123500e+03, new int[] { 2, 1, 4 });
            p.AddCoeff(-2.8533739823547500e+03, new int[] { 2, 1, 6 });
            p.AddCoeff(1.5285932048329000e+03, new int[] { 2, 1, 8 });
            p.AddCoeff(-3.8798791523860200e+01, new int[] { 2, 3, 0 });
            p.AddCoeff(1.3967564948589700e+03, new int[] { 2, 3, 2 });
            p.AddCoeff(-7.6821607217243200e+03, new int[] { 2, 3, 4 });
            p.AddCoeff(1.3315745250988800e+04, new int[] { 2, 3, 6 });
            p.AddCoeff(-7.1334349558868700e+03, new int[] { 2, 3, 8 });
            p.AddCoeff(3.4918912371474200e+01, new int[] { 2, 5, 0 });
            p.AddCoeff(-1.2570808453730700e+03, new int[] { 2, 5, 2 });
            p.AddCoeff(6.9139446495518900e+03, new int[] { 2, 5, 4 });
            p.AddCoeff(-1.1984170725889900e+04, new int[] { 2, 5, 6 });
            p.AddCoeff(6.4200914602982000e+03, new int[] { 2, 5, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("211132c1-2741-4e6a-b42b-06f13716f202"));
            OrthonormalPolynomials[717] = p;
            p.AddCoeff(-3.7733353310726900e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(3.3960017979654200e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-7.4712039555239300e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(4.6250310200862400e+01, new int[] { 0, 0, 7 });
            p.AddCoeff(7.9240041952526500e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-7.1316037757273800e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(1.5689528306600200e+03, new int[] { 0, 2, 5 });
            p.AddCoeff(-9.7125651421811000e+02, new int[] { 0, 2, 7 });
            p.AddCoeff(-2.3772012585757900e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(2.1394811327182200e+03, new int[] { 0, 4, 3 });
            p.AddCoeff(-4.7068584919800700e+03, new int[] { 0, 4, 5 });
            p.AddCoeff(2.9137695426543300e+03, new int[] { 0, 4, 7 });
            p.AddCoeff(1.7432809229555800e+02, new int[] { 0, 6, 1 });
            p.AddCoeff(-1.5689528306600200e+03, new int[] { 0, 6, 3 });
            p.AddCoeff(3.4516962274520500e+03, new int[] { 0, 6, 5 });
            p.AddCoeff(-2.1367643312798400e+03, new int[] { 0, 6, 7 });
            p.AddCoeff(1.1320005993218100e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-1.0188005393896300e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(2.2413611866571800e+02, new int[] { 2, 0, 5 });
            p.AddCoeff(-1.3875093060258700e+02, new int[] { 2, 0, 7 });
            p.AddCoeff(-2.3772012585757900e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(2.1394811327182200e+03, new int[] { 2, 2, 3 });
            p.AddCoeff(-4.7068584919800700e+03, new int[] { 2, 2, 5 });
            p.AddCoeff(2.9137695426543300e+03, new int[] { 2, 2, 7 });
            p.AddCoeff(7.1316037757273800e+02, new int[] { 2, 4, 1 });
            p.AddCoeff(-6.4184433981546500e+03, new int[] { 2, 4, 3 });
            p.AddCoeff(1.4120575475940200e+04, new int[] { 2, 4, 5 });
            p.AddCoeff(-8.7413086279629900e+03, new int[] { 2, 4, 7 });
            p.AddCoeff(-5.2298427688667500e+02, new int[] { 2, 6, 1 });
            p.AddCoeff(4.7068584919800700e+03, new int[] { 2, 6, 3 });
            p.AddCoeff(-1.0355088682356200e+04, new int[] { 2, 6, 5 });
            p.AddCoeff(6.4102929938395300e+03, new int[] { 2, 6, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d914d177-4345-4f10-8bfe-cad724a64188"));
            OrthonormalPolynomials[718] = p;
            p.AddCoeff(-3.7733353310726900e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(7.9240041952526500e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-2.3772012585757900e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(1.7432809229555800e+02, new int[] { 0, 1, 6 });
            p.AddCoeff(3.3960017979654200e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-7.1316037757273800e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(2.1394811327182200e+03, new int[] { 0, 3, 4 });
            p.AddCoeff(-1.5689528306600200e+03, new int[] { 0, 3, 6 });
            p.AddCoeff(-7.4712039555239300e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(1.5689528306600200e+03, new int[] { 0, 5, 2 });
            p.AddCoeff(-4.7068584919800700e+03, new int[] { 0, 5, 4 });
            p.AddCoeff(3.4516962274520500e+03, new int[] { 0, 5, 6 });
            p.AddCoeff(4.6250310200862400e+01, new int[] { 0, 7, 0 });
            p.AddCoeff(-9.7125651421811000e+02, new int[] { 0, 7, 2 });
            p.AddCoeff(2.9137695426543300e+03, new int[] { 0, 7, 4 });
            p.AddCoeff(-2.1367643312798400e+03, new int[] { 0, 7, 6 });
            p.AddCoeff(1.1320005993218100e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-2.3772012585757900e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(7.1316037757273800e+02, new int[] { 2, 1, 4 });
            p.AddCoeff(-5.2298427688667500e+02, new int[] { 2, 1, 6 });
            p.AddCoeff(-1.0188005393896300e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(2.1394811327182200e+03, new int[] { 2, 3, 2 });
            p.AddCoeff(-6.4184433981546500e+03, new int[] { 2, 3, 4 });
            p.AddCoeff(4.7068584919800700e+03, new int[] { 2, 3, 6 });
            p.AddCoeff(2.2413611866571800e+02, new int[] { 2, 5, 0 });
            p.AddCoeff(-4.7068584919800700e+03, new int[] { 2, 5, 2 });
            p.AddCoeff(1.4120575475940200e+04, new int[] { 2, 5, 4 });
            p.AddCoeff(-1.0355088682356200e+04, new int[] { 2, 5, 6 });
            p.AddCoeff(-1.3875093060258700e+02, new int[] { 2, 7, 0 });
            p.AddCoeff(2.9137695426543300e+03, new int[] { 2, 7, 2 });
            p.AddCoeff(-8.7413086279629900e+03, new int[] { 2, 7, 4 });
            p.AddCoeff(6.4102929938395300e+03, new int[] { 2, 7, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a2d4d347-44d0-49bd-85ad-3da17e610a1c"));
            OrthonormalPolynomials[719] = p;
            p.AddCoeff(-2.7713422517043000e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(1.2932930507953400e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-1.1639637457158100e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(9.9768321061354800e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-4.6558549828632300e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(4.1902694845769000e+02, new int[] { 0, 2, 5 });
            p.AddCoeff(-5.4872576583745200e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(2.5607202405747700e+03, new int[] { 0, 4, 3 });
            p.AddCoeff(-2.3046482165173000e+03, new int[] { 0, 4, 5 });
            p.AddCoeff(9.5112466078491600e+02, new int[] { 0, 6, 1 });
            p.AddCoeff(-4.4385817503296100e+03, new int[] { 0, 6, 3 });
            p.AddCoeff(3.9947235752966500e+03, new int[] { 0, 6, 5 });
            p.AddCoeff(-5.0953106827763400e+02, new int[] { 0, 8, 1 });
            p.AddCoeff(2.3778116519622900e+03, new int[] { 0, 8, 3 });
            p.AddCoeff(-2.1400304867660600e+03, new int[] { 0, 8, 5 });
            p.AddCoeff(8.3140267551129000e+00, new int[] { 2, 0, 1 });
            p.AddCoeff(-3.8798791523860200e+01, new int[] { 2, 0, 3 });
            p.AddCoeff(3.4918912371474200e+01, new int[] { 2, 0, 5 });
            p.AddCoeff(-2.9930496318406400e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(1.3967564948589700e+03, new int[] { 2, 2, 3 });
            p.AddCoeff(-1.2570808453730700e+03, new int[] { 2, 2, 5 });
            p.AddCoeff(1.6461772975123500e+03, new int[] { 2, 4, 1 });
            p.AddCoeff(-7.6821607217243200e+03, new int[] { 2, 4, 3 });
            p.AddCoeff(6.9139446495518900e+03, new int[] { 2, 4, 5 });
            p.AddCoeff(-2.8533739823547500e+03, new int[] { 2, 6, 1 });
            p.AddCoeff(1.3315745250988800e+04, new int[] { 2, 6, 3 });
            p.AddCoeff(-1.1984170725889900e+04, new int[] { 2, 6, 5 });
            p.AddCoeff(1.5285932048329000e+03, new int[] { 2, 8, 1 });
            p.AddCoeff(-7.1334349558868700e+03, new int[] { 2, 8, 3 });
            p.AddCoeff(6.4200914602982000e+03, new int[] { 2, 8, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e30f0e7d-c190-466a-9db1-42d9d1a5efdc"));
            OrthonormalPolynomials[720] = p;
            p.AddCoeff(-4.7702365981659800e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(4.7702365981659800e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-5.5652760311936400e+01, new int[] { 0, 1, 4 });
            p.AddCoeff(6.9963470106434400e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-6.9963470106434400e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(8.1624048457506800e+02, new int[] { 0, 3, 4 });
            p.AddCoeff(-2.7285753341509400e+02, new int[] { 0, 5, 0 });
            p.AddCoeff(2.7285753341509400e+03, new int[] { 0, 5, 2 });
            p.AddCoeff(-3.1833378898427600e+03, new int[] { 0, 5, 4 });
            p.AddCoeff(3.8979647630727700e+02, new int[] { 0, 7, 0 });
            p.AddCoeff(-3.8979647630727700e+03, new int[] { 0, 7, 2 });
            p.AddCoeff(4.5476255569182400e+03, new int[] { 0, 7, 4 });
            p.AddCoeff(-1.8407055825621400e+02, new int[] { 0, 9, 0 });
            p.AddCoeff(1.8407055825621400e+03, new int[] { 0, 9, 2 });
            p.AddCoeff(-2.1474898463225100e+03, new int[] { 0, 9, 4 });
            p.AddCoeff(1.4310709794497900e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-1.4310709794497900e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(1.6695828093580900e+02, new int[] { 2, 1, 4 });
            p.AddCoeff(-2.0989041031930300e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(2.0989041031930300e+03, new int[] { 2, 3, 2 });
            p.AddCoeff(-2.4487214537252000e+03, new int[] { 2, 3, 4 });
            p.AddCoeff(8.1857260024528200e+02, new int[] { 2, 5, 0 });
            p.AddCoeff(-8.1857260024528200e+03, new int[] { 2, 5, 2 });
            p.AddCoeff(9.5500136695282900e+03, new int[] { 2, 5, 4 });
            p.AddCoeff(-1.1693894289218300e+03, new int[] { 2, 7, 0 });
            p.AddCoeff(1.1693894289218300e+04, new int[] { 2, 7, 2 });
            p.AddCoeff(-1.3642876670754700e+04, new int[] { 2, 7, 4 });
            p.AddCoeff(5.5221167476864300e+02, new int[] { 2, 9, 0 });
            p.AddCoeff(-5.5221167476864300e+03, new int[] { 2, 9, 2 });
            p.AddCoeff(6.4424695389674900e+03, new int[] { 2, 9, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("36e52e84-4e70-4415-ad10-6b74d29c6721"));
            OrthonormalPolynomials[721] = p;
            p.AddCoeff(-1.7691331630354800e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(2.9485552717258100e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(9.7302323966951700e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-1.6217053994491900e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(-8.4328680771358100e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(1.4054780128559700e+03, new int[] { 0, 4, 3 });
            p.AddCoeff(2.5298604231407400e+03, new int[] { 0, 6, 1 });
            p.AddCoeff(-4.2164340385679100e+03, new int[] { 0, 6, 3 });
            p.AddCoeff(-3.0719733709566200e+03, new int[] { 0, 8, 1 });
            p.AddCoeff(5.1199556182610300e+03, new int[] { 0, 8, 3 });
            p.AddCoeff(1.2970554232927900e+03, new int[] { 0, 10, 1 });
            p.AddCoeff(-2.1617590388213200e+03, new int[] { 0, 10, 3 });
            p.AddCoeff(5.3073994891064500e+00, new int[] { 2, 0, 1 });
            p.AddCoeff(-8.8456658151774200e+00, new int[] { 2, 0, 3 });
            p.AddCoeff(-2.9190697190085500e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(4.8651161983475800e+02, new int[] { 2, 2, 3 });
            p.AddCoeff(2.5298604231407400e+03, new int[] { 2, 4, 1 });
            p.AddCoeff(-4.2164340385679100e+03, new int[] { 2, 4, 3 });
            p.AddCoeff(-7.5895812694222300e+03, new int[] { 2, 6, 1 });
            p.AddCoeff(1.2649302115703700e+04, new int[] { 2, 6, 3 });
            p.AddCoeff(9.2159201128698200e+03, new int[] { 2, 8, 1 });
            p.AddCoeff(-1.5359866854783100e+04, new int[] { 2, 8, 3 });
            p.AddCoeff(-3.8911662698783800e+03, new int[] { 2, 10, 1 });
            p.AddCoeff(6.4852771164639400e+03, new int[] { 2, 10, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("472b8a98-01d8-4c1b-a9c2-239114d3424c"));
            OrthonormalPolynomials[722] = p;
            p.AddCoeff(-5.7374935037912700e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(1.7212480511373800e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(1.2431235924881100e+02, new int[] { 0, 3, 0 });
            p.AddCoeff(-3.7293707774643300e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(-7.4587415549286600e+02, new int[] { 0, 5, 0 });
            p.AddCoeff(2.2376224664786000e+03, new int[] { 0, 5, 2 });
            p.AddCoeff(1.8114086633398200e+03, new int[] { 0, 7, 0 });
            p.AddCoeff(-5.4342259900194500e+03, new int[] { 0, 7, 2 });
            p.AddCoeff(-1.9120424779698100e+03, new int[] { 0, 9, 0 });
            p.AddCoeff(5.7361274339094200e+03, new int[] { 0, 9, 2 });
            p.AddCoeff(7.3005258249756200e+02, new int[] { 0, 11, 0 });
            p.AddCoeff(-2.1901577474926900e+03, new int[] { 0, 11, 2 });
            p.AddCoeff(1.7212480511373800e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-5.1637441534121500e+01, new int[] { 2, 1, 2 });
            p.AddCoeff(-3.7293707774643300e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(1.1188112332393000e+03, new int[] { 2, 3, 2 });
            p.AddCoeff(2.2376224664786000e+03, new int[] { 2, 5, 0 });
            p.AddCoeff(-6.7128673994357900e+03, new int[] { 2, 5, 2 });
            p.AddCoeff(-5.4342259900194500e+03, new int[] { 2, 7, 0 });
            p.AddCoeff(1.6302677970058400e+04, new int[] { 2, 7, 2 });
            p.AddCoeff(5.7361274339094200e+03, new int[] { 2, 9, 0 });
            p.AddCoeff(-1.7208382301728200e+04, new int[] { 2, 9, 2 });
            p.AddCoeff(-2.1901577474926900e+03, new int[] { 2, 11, 0 });
            p.AddCoeff(6.5704732424780600e+03, new int[] { 2, 11, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("14993090-40ae-4997-81b5-62307fd351bb"));
            OrthonormalPolynomials[723] = p;
            p.AddCoeff(-7.7224066640437800e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(6.0234771979541500e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-7.5293464974426600e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(3.4133037455073500e+03, new int[] { 0, 6, 1 });
            p.AddCoeff(-6.9485111962113800e+03, new int[] { 0, 8, 1 });
            p.AddCoeff(6.4852771164639400e+03, new int[] { 0, 10, 1 });
            p.AddCoeff(-2.2600208133132000e+03, new int[] { 0, 12, 1 });
            p.AddCoeff(2.3167219992131300e+00, new int[] { 2, 0, 1 });
            p.AddCoeff(-1.8070431593862500e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(2.2588039492328100e+03, new int[] { 2, 4, 1 });
            p.AddCoeff(-1.0239911236522100e+04, new int[] { 2, 6, 1 });
            p.AddCoeff(2.0845533588634200e+04, new int[] { 2, 8, 1 });
            p.AddCoeff(-1.9455831349391900e+04, new int[] { 2, 10, 1 });
            p.AddCoeff(6.7800624399395900e+03, new int[] { 2, 12, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("bbc6f461-a93a-4ae6-9a45-f05d086478de"));
            OrthonormalPolynomials[724] = p;
            p.AddCoeff(-6.0234771979541500e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(1.8070431593862500e+02, new int[] { 0, 3, 0 });
            p.AddCoeff(-1.5359866854783100e+03, new int[] { 0, 5, 0 });
            p.AddCoeff(5.5588089569691200e+03, new int[] { 0, 7, 0 });
            p.AddCoeff(-9.7279156746959700e+03, new int[] { 0, 9, 0 });
            p.AddCoeff(8.1360749279275500e+03, new int[] { 0, 11, 0 });
            p.AddCoeff(-2.6077163230536900e+03, new int[] { 0, 13, 0 });
            p.AddCoeff(1.8070431593862500e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-5.4211294781587400e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(4.6079600564349200e+03, new int[] { 2, 5, 0 });
            p.AddCoeff(-1.6676426870907300e+04, new int[] { 2, 7, 0 });
            p.AddCoeff(2.9183747024087800e+04, new int[] { 2, 9, 0 });
            p.AddCoeff(-2.4408224783782500e+04, new int[] { 2, 11, 0 });
            p.AddCoeff(7.8231489691611000e+03, new int[] { 2, 13, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("76799c9b-b611-4a09-88a4-bcf5d982b1c4"));
            OrthonormalPolynomials[725] = p;
            p.AddCoeff(-1.5826224176235000e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(1.2344454857463300e+02, new int[] { 1, 0, 2 });
            p.AddCoeff(-1.5430568571829100e+03, new int[] { 1, 0, 4 });
            p.AddCoeff(6.9951910858958600e+03, new int[] { 1, 0, 6 });
            p.AddCoeff(-1.4240210424859400e+04, new int[] { 1, 0, 8 });
            p.AddCoeff(1.3290863063202200e+04, new int[] { 1, 0, 10 });
            p.AddCoeff(-4.6316644008128600e+03, new int[] { 1, 0, 12 });
            p.AddCoeff(2.6377040293725000e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(-2.0574091429105500e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(2.5717614286381800e+03, new int[] { 3, 0, 4 });
            p.AddCoeff(-1.1658651809826400e+04, new int[] { 3, 0, 6 });
            p.AddCoeff(2.3733684041432400e+04, new int[] { 3, 0, 8 });
            p.AddCoeff(-2.2151438438670200e+04, new int[] { 3, 0, 10 });
            p.AddCoeff(7.7194406680214600e+03, new int[] { 3, 0, 12 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("bec7e70a-4cfa-477d-a92a-f1b56bbf09f0"));
            OrthonormalPolynomials[726] = p;
            p.AddCoeff(3.1550997936529100e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-6.8360495529146400e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(4.1016297317487800e+03, new int[] { 1, 1, 5 });
            p.AddCoeff(-9.9611007771041800e+03, new int[] { 1, 1, 7 });
            p.AddCoeff(1.0514495264721100e+04, new int[] { 1, 1, 9 });
            p.AddCoeff(-4.0146254647116900e+03, new int[] { 1, 1, 11 });
            p.AddCoeff(-5.2584996560881800e+01, new int[] { 3, 1, 1 });
            p.AddCoeff(1.1393415921524400e+03, new int[] { 3, 1, 3 });
            p.AddCoeff(-6.8360495529146400e+03, new int[] { 3, 1, 5 });
            p.AddCoeff(1.6601834628507000e+04, new int[] { 3, 1, 7 });
            p.AddCoeff(-1.7524158774535100e+04, new int[] { 3, 1, 9 });
            p.AddCoeff(6.6910424411861400e+03, new int[] { 3, 1, 11 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("8465b7dd-1afa-4f3e-ab8b-289a46c0dceb"));
            OrthonormalPolynomials[727] = p;
            p.AddCoeff(-1.7691331630354800e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(9.7302323966951700e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-8.4328680771358100e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(2.5298604231407400e+03, new int[] { 1, 0, 6 });
            p.AddCoeff(-3.0719733709566200e+03, new int[] { 1, 0, 8 });
            p.AddCoeff(1.2970554232927900e+03, new int[] { 1, 0, 10 });
            p.AddCoeff(5.3073994891064500e+00, new int[] { 1, 2, 0 });
            p.AddCoeff(-2.9190697190085500e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(2.5298604231407400e+03, new int[] { 1, 2, 4 });
            p.AddCoeff(-7.5895812694222300e+03, new int[] { 1, 2, 6 });
            p.AddCoeff(9.2159201128698200e+03, new int[] { 1, 2, 8 });
            p.AddCoeff(-3.8911662698783800e+03, new int[] { 1, 2, 10 });
            p.AddCoeff(2.9485552717258100e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(-1.6217053994491900e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(1.4054780128559700e+03, new int[] { 3, 0, 4 });
            p.AddCoeff(-4.2164340385679100e+03, new int[] { 3, 0, 6 });
            p.AddCoeff(5.1199556182610300e+03, new int[] { 3, 0, 8 });
            p.AddCoeff(-2.1617590388213200e+03, new int[] { 3, 0, 10 });
            p.AddCoeff(-8.8456658151774200e+00, new int[] { 3, 2, 0 });
            p.AddCoeff(4.8651161983475800e+02, new int[] { 3, 2, 2 });
            p.AddCoeff(-4.2164340385679100e+03, new int[] { 3, 2, 4 });
            p.AddCoeff(1.2649302115703700e+04, new int[] { 3, 2, 6 });
            p.AddCoeff(-1.5359866854783100e+04, new int[] { 3, 2, 8 });
            p.AddCoeff(6.4852771164639400e+03, new int[] { 3, 2, 10 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("0ab176e7-f932-45dd-9fce-6bb90f3e2d5c"));
            OrthonormalPolynomials[728] = p;
            p.AddCoeff(5.9732810492636400e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-8.7608122055866800e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(3.4167167601788000e+03, new int[] { 1, 1, 5 });
            p.AddCoeff(-4.8810239431125800e+03, new int[] { 1, 1, 7 });
            p.AddCoeff(2.3049279731364900e+03, new int[] { 1, 1, 9 });
            p.AddCoeff(-9.9554684154394000e+01, new int[] { 1, 3, 1 });
            p.AddCoeff(1.4601353675977800e+03, new int[] { 1, 3, 3 });
            p.AddCoeff(-5.6945279336313400e+03, new int[] { 1, 3, 5 });
            p.AddCoeff(8.1350399051876300e+03, new int[] { 1, 3, 7 });
            p.AddCoeff(-3.8415466218941600e+03, new int[] { 1, 3, 9 });
            p.AddCoeff(-9.9554684154394000e+01, new int[] { 3, 1, 1 });
            p.AddCoeff(1.4601353675977800e+03, new int[] { 3, 1, 3 });
            p.AddCoeff(-5.6945279336313400e+03, new int[] { 3, 1, 5 });
            p.AddCoeff(8.1350399051876300e+03, new int[] { 3, 1, 7 });
            p.AddCoeff(-3.8415466218941600e+03, new int[] { 3, 1, 9 });
            p.AddCoeff(1.6592447359065700e+02, new int[] { 3, 3, 1 });
            p.AddCoeff(-2.4335589459963000e+03, new int[] { 3, 3, 3 });
            p.AddCoeff(9.4908798893855700e+03, new int[] { 3, 3, 5 });
            p.AddCoeff(-1.3558399841979400e+04, new int[] { 3, 3, 7 });
            p.AddCoeff(6.4025777031569300e+03, new int[] { 3, 3, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("5fab9df1-480f-4962-b62d-88ccb83e3504"));
            OrthonormalPolynomials[729] = p;
            p.AddCoeff(-1.7796325618178400e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(6.4066772225442100e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-3.5236724723993200e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(6.1076989521588200e+02, new int[] { 1, 0, 6 });
            p.AddCoeff(-3.2719815815136500e+02, new int[] { 1, 0, 8 });
            p.AddCoeff(1.7796325618178400e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-6.4066772225442100e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(3.5236724723993200e+03, new int[] { 1, 2, 4 });
            p.AddCoeff(-6.1076989521588200e+03, new int[] { 1, 2, 6 });
            p.AddCoeff(3.2719815815136500e+03, new int[] { 1, 2, 8 });
            p.AddCoeff(-2.0762379887874800e+01, new int[] { 1, 4, 0 });
            p.AddCoeff(7.4744567596349100e+02, new int[] { 1, 4, 2 });
            p.AddCoeff(-4.1109512177992000e+03, new int[] { 1, 4, 4 });
            p.AddCoeff(7.1256487775186200e+03, new int[] { 1, 4, 6 });
            p.AddCoeff(-3.8173118450992500e+03, new int[] { 1, 4, 8 });
            p.AddCoeff(2.9660542696963900e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(-1.0677795370907000e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(5.8727874539988600e+02, new int[] { 3, 0, 4 });
            p.AddCoeff(-1.0179498253598000e+03, new int[] { 3, 0, 6 });
            p.AddCoeff(5.4533026358560900e+02, new int[] { 3, 0, 8 });
            p.AddCoeff(-2.9660542696963900e+01, new int[] { 3, 2, 0 });
            p.AddCoeff(1.0677795370907000e+03, new int[] { 3, 2, 2 });
            p.AddCoeff(-5.8727874539988600e+03, new int[] { 3, 2, 4 });
            p.AddCoeff(1.0179498253598000e+04, new int[] { 3, 2, 6 });
            p.AddCoeff(-5.4533026358560900e+03, new int[] { 3, 2, 8 });
            p.AddCoeff(3.4603966479791300e+01, new int[] { 3, 4, 0 });
            p.AddCoeff(-1.2457427932724900e+03, new int[] { 3, 4, 2 });
            p.AddCoeff(6.8515853629986700e+03, new int[] { 3, 4, 4 });
            p.AddCoeff(-1.1876081295864400e+04, new int[] { 3, 4, 6 });
            p.AddCoeff(6.3621864084987700e+03, new int[] { 3, 4, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("4d1414b5-70c8-4478-b0ad-bd4c0fdaff73"));
            OrthonormalPolynomials[730] = p;
            p.AddCoeff(7.3924192867575100e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-6.6531773580817600e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(1.4636990187779900e+03, new int[] { 1, 1, 5 });
            p.AddCoeff(-9.0609939257684900e+02, new int[] { 1, 1, 7 });
            p.AddCoeff(-3.4497956671535100e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(3.1048161004381600e+03, new int[] { 1, 3, 3 });
            p.AddCoeff(-6.8305954209639400e+03, new int[] { 1, 3, 5 });
            p.AddCoeff(4.2284638320253000e+03, new int[] { 1, 3, 7 });
            p.AddCoeff(3.1048161004381600e+02, new int[] { 1, 5, 1 });
            p.AddCoeff(-2.7943344903943400e+03, new int[] { 1, 5, 3 });
            p.AddCoeff(6.1475358788675500e+03, new int[] { 1, 5, 5 });
            p.AddCoeff(-3.8056174488227700e+03, new int[] { 1, 5, 7 });
            p.AddCoeff(-1.2320698811262500e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(1.1088628930136300e+03, new int[] { 3, 1, 3 });
            p.AddCoeff(-2.4394983646299800e+03, new int[] { 3, 1, 5 });
            p.AddCoeff(1.5101656542947500e+03, new int[] { 3, 1, 7 });
            p.AddCoeff(5.7496594452558400e+02, new int[] { 3, 3, 1 });
            p.AddCoeff(-5.1746935007302600e+03, new int[] { 3, 3, 3 });
            p.AddCoeff(1.1384325701606600e+04, new int[] { 3, 3, 5 });
            p.AddCoeff(-7.0474397200421600e+03, new int[] { 3, 3, 7 });
            p.AddCoeff(-5.1746935007302600e+02, new int[] { 3, 5, 1 });
            p.AddCoeff(4.6572241506572300e+03, new int[] { 3, 5, 3 });
            p.AddCoeff(-1.0245893131445900e+04, new int[] { 3, 5, 5 });
            p.AddCoeff(6.3426957480379500e+03, new int[] { 3, 5, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("0447f2ca-edb2-4dd9-be61-b409d97ef320"));
            OrthonormalPolynomials[731] = p;
            p.AddCoeff(-1.7813066172385700e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(3.7407438962010000e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-1.1222231688603000e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(8.2296365716421900e+01, new int[] { 1, 0, 6 });
            p.AddCoeff(3.7407438962010000e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-7.8555621820220900e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(2.3566686546066300e+03, new int[] { 1, 2, 4 });
            p.AddCoeff(-1.7282236800448600e+03, new int[] { 1, 2, 6 });
            p.AddCoeff(-1.1222231688603000e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(2.3566686546066300e+03, new int[] { 1, 4, 2 });
            p.AddCoeff(-7.0700059638198800e+03, new int[] { 1, 4, 4 });
            p.AddCoeff(5.1846710401345800e+03, new int[] { 1, 4, 6 });
            p.AddCoeff(8.2296365716421900e+01, new int[] { 1, 6, 0 });
            p.AddCoeff(-1.7282236800448600e+03, new int[] { 1, 6, 2 });
            p.AddCoeff(5.1846710401345800e+03, new int[] { 1, 6, 4 });
            p.AddCoeff(-3.8020920960986900e+03, new int[] { 1, 6, 6 });
            p.AddCoeff(2.9688443620642800e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(-6.2345731603350000e+01, new int[] { 3, 0, 2 });
            p.AddCoeff(1.8703719481005000e+02, new int[] { 3, 0, 4 });
            p.AddCoeff(-1.3716060952737000e+02, new int[] { 3, 0, 6 });
            p.AddCoeff(-6.2345731603350000e+01, new int[] { 3, 2, 0 });
            p.AddCoeff(1.3092603636703500e+03, new int[] { 3, 2, 2 });
            p.AddCoeff(-3.9277810910110500e+03, new int[] { 3, 2, 4 });
            p.AddCoeff(2.8803728000747700e+03, new int[] { 3, 2, 6 });
            p.AddCoeff(1.8703719481005000e+02, new int[] { 3, 4, 0 });
            p.AddCoeff(-3.9277810910110500e+03, new int[] { 3, 4, 2 });
            p.AddCoeff(1.1783343273033100e+04, new int[] { 3, 4, 4 });
            p.AddCoeff(-8.6411184002243000e+03, new int[] { 3, 4, 6 });
            p.AddCoeff(-1.3716060952737000e+02, new int[] { 3, 6, 0 });
            p.AddCoeff(2.8803728000747700e+03, new int[] { 3, 6, 2 });
            p.AddCoeff(-8.6411184002243000e+03, new int[] { 3, 6, 4 });
            p.AddCoeff(6.3368201601644900e+03, new int[] { 3, 6, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b19c25db-8683-405b-add9-d78c5e43e6ee"));
            OrthonormalPolynomials[732] = p;
            p.AddCoeff(7.3924192867575100e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-3.4497956671535100e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(3.1048161004381600e+02, new int[] { 1, 1, 5 });
            p.AddCoeff(-6.6531773580817600e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(3.1048161004381600e+03, new int[] { 1, 3, 3 });
            p.AddCoeff(-2.7943344903943400e+03, new int[] { 1, 3, 5 });
            p.AddCoeff(1.4636990187779900e+03, new int[] { 1, 5, 1 });
            p.AddCoeff(-6.8305954209639400e+03, new int[] { 1, 5, 3 });
            p.AddCoeff(6.1475358788675500e+03, new int[] { 1, 5, 5 });
            p.AddCoeff(-9.0609939257684900e+02, new int[] { 1, 7, 1 });
            p.AddCoeff(4.2284638320253000e+03, new int[] { 1, 7, 3 });
            p.AddCoeff(-3.8056174488227700e+03, new int[] { 1, 7, 5 });
            p.AddCoeff(-1.2320698811262500e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(5.7496594452558400e+02, new int[] { 3, 1, 3 });
            p.AddCoeff(-5.1746935007302600e+02, new int[] { 3, 1, 5 });
            p.AddCoeff(1.1088628930136300e+03, new int[] { 3, 3, 1 });
            p.AddCoeff(-5.1746935007302600e+03, new int[] { 3, 3, 3 });
            p.AddCoeff(4.6572241506572300e+03, new int[] { 3, 3, 5 });
            p.AddCoeff(-2.4394983646299800e+03, new int[] { 3, 5, 1 });
            p.AddCoeff(1.1384325701606600e+04, new int[] { 3, 5, 3 });
            p.AddCoeff(-1.0245893131445900e+04, new int[] { 3, 5, 5 });
            p.AddCoeff(1.5101656542947500e+03, new int[] { 3, 7, 1 });
            p.AddCoeff(-7.0474397200421600e+03, new int[] { 3, 7, 3 });
            p.AddCoeff(6.3426957480379500e+03, new int[] { 3, 7, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c907bd1f-e495-43f5-97c8-ef7738a90d01"));
            OrthonormalPolynomials[733] = p;
            p.AddCoeff(-1.7796325618178400e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(1.7796325618178400e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-2.0762379887874800e+01, new int[] { 1, 0, 4 });
            p.AddCoeff(6.4066772225442100e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-6.4066772225442100e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(7.4744567596349100e+02, new int[] { 1, 2, 4 });
            p.AddCoeff(-3.5236724723993200e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(3.5236724723993200e+03, new int[] { 1, 4, 2 });
            p.AddCoeff(-4.1109512177992000e+03, new int[] { 1, 4, 4 });
            p.AddCoeff(6.1076989521588200e+02, new int[] { 1, 6, 0 });
            p.AddCoeff(-6.1076989521588200e+03, new int[] { 1, 6, 2 });
            p.AddCoeff(7.1256487775186200e+03, new int[] { 1, 6, 4 });
            p.AddCoeff(-3.2719815815136500e+02, new int[] { 1, 8, 0 });
            p.AddCoeff(3.2719815815136500e+03, new int[] { 1, 8, 2 });
            p.AddCoeff(-3.8173118450992500e+03, new int[] { 1, 8, 4 });
            p.AddCoeff(2.9660542696963900e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(-2.9660542696963900e+01, new int[] { 3, 0, 2 });
            p.AddCoeff(3.4603966479791300e+01, new int[] { 3, 0, 4 });
            p.AddCoeff(-1.0677795370907000e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(1.0677795370907000e+03, new int[] { 3, 2, 2 });
            p.AddCoeff(-1.2457427932724900e+03, new int[] { 3, 2, 4 });
            p.AddCoeff(5.8727874539988600e+02, new int[] { 3, 4, 0 });
            p.AddCoeff(-5.8727874539988600e+03, new int[] { 3, 4, 2 });
            p.AddCoeff(6.8515853629986700e+03, new int[] { 3, 4, 4 });
            p.AddCoeff(-1.0179498253598000e+03, new int[] { 3, 6, 0 });
            p.AddCoeff(1.0179498253598000e+04, new int[] { 3, 6, 2 });
            p.AddCoeff(-1.1876081295864400e+04, new int[] { 3, 6, 4 });
            p.AddCoeff(5.4533026358560900e+02, new int[] { 3, 8, 0 });
            p.AddCoeff(-5.4533026358560900e+03, new int[] { 3, 8, 2 });
            p.AddCoeff(6.3621864084987700e+03, new int[] { 3, 8, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("30a33e3d-16e0-49c1-b4fa-90c09f32f7c6"));
            OrthonormalPolynomials[734] = p;
            p.AddCoeff(5.9732810492636400e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-9.9554684154394000e+01, new int[] { 1, 1, 3 });
            p.AddCoeff(-8.7608122055866800e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(1.4601353675977800e+03, new int[] { 1, 3, 3 });
            p.AddCoeff(3.4167167601788000e+03, new int[] { 1, 5, 1 });
            p.AddCoeff(-5.6945279336313400e+03, new int[] { 1, 5, 3 });
            p.AddCoeff(-4.8810239431125800e+03, new int[] { 1, 7, 1 });
            p.AddCoeff(8.1350399051876300e+03, new int[] { 1, 7, 3 });
            p.AddCoeff(2.3049279731364900e+03, new int[] { 1, 9, 1 });
            p.AddCoeff(-3.8415466218941600e+03, new int[] { 1, 9, 3 });
            p.AddCoeff(-9.9554684154394000e+01, new int[] { 3, 1, 1 });
            p.AddCoeff(1.6592447359065700e+02, new int[] { 3, 1, 3 });
            p.AddCoeff(1.4601353675977800e+03, new int[] { 3, 3, 1 });
            p.AddCoeff(-2.4335589459963000e+03, new int[] { 3, 3, 3 });
            p.AddCoeff(-5.6945279336313400e+03, new int[] { 3, 5, 1 });
            p.AddCoeff(9.4908798893855700e+03, new int[] { 3, 5, 3 });
            p.AddCoeff(8.1350399051876300e+03, new int[] { 3, 7, 1 });
            p.AddCoeff(-1.3558399841979400e+04, new int[] { 3, 7, 3 });
            p.AddCoeff(-3.8415466218941600e+03, new int[] { 3, 9, 1 });
            p.AddCoeff(6.4025777031569300e+03, new int[] { 3, 9, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3dbc3b46-0436-4004-bdbc-9073385b5048"));
            OrthonormalPolynomials[735] = p;
            p.AddCoeff(-1.7691331630354800e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(5.3073994891064500e+00, new int[] { 1, 0, 2 });
            p.AddCoeff(9.7302323966951700e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-2.9190697190085500e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(-8.4328680771358100e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(2.5298604231407400e+03, new int[] { 1, 4, 2 });
            p.AddCoeff(2.5298604231407400e+03, new int[] { 1, 6, 0 });
            p.AddCoeff(-7.5895812694222300e+03, new int[] { 1, 6, 2 });
            p.AddCoeff(-3.0719733709566200e+03, new int[] { 1, 8, 0 });
            p.AddCoeff(9.2159201128698200e+03, new int[] { 1, 8, 2 });
            p.AddCoeff(1.2970554232927900e+03, new int[] { 1, 10, 0 });
            p.AddCoeff(-3.8911662698783800e+03, new int[] { 1, 10, 2 });
            p.AddCoeff(2.9485552717258100e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(-8.8456658151774200e+00, new int[] { 3, 0, 2 });
            p.AddCoeff(-1.6217053994491900e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(4.8651161983475800e+02, new int[] { 3, 2, 2 });
            p.AddCoeff(1.4054780128559700e+03, new int[] { 3, 4, 0 });
            p.AddCoeff(-4.2164340385679100e+03, new int[] { 3, 4, 2 });
            p.AddCoeff(-4.2164340385679100e+03, new int[] { 3, 6, 0 });
            p.AddCoeff(1.2649302115703700e+04, new int[] { 3, 6, 2 });
            p.AddCoeff(5.1199556182610300e+03, new int[] { 3, 8, 0 });
            p.AddCoeff(-1.5359866854783100e+04, new int[] { 3, 8, 2 });
            p.AddCoeff(-2.1617590388213200e+03, new int[] { 3, 10, 0 });
            p.AddCoeff(6.4852771164639400e+03, new int[] { 3, 10, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7b6c527c-c983-4d23-bb16-43c7cb5946e8"));
            OrthonormalPolynomials[736] = p;
            p.AddCoeff(3.1550997936529100e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-6.8360495529146400e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(4.1016297317487800e+03, new int[] { 1, 5, 1 });
            p.AddCoeff(-9.9611007771041800e+03, new int[] { 1, 7, 1 });
            p.AddCoeff(1.0514495264721100e+04, new int[] { 1, 9, 1 });
            p.AddCoeff(-4.0146254647116900e+03, new int[] { 1, 11, 1 });
            p.AddCoeff(-5.2584996560881800e+01, new int[] { 3, 1, 1 });
            p.AddCoeff(1.1393415921524400e+03, new int[] { 3, 3, 1 });
            p.AddCoeff(-6.8360495529146400e+03, new int[] { 3, 5, 1 });
            p.AddCoeff(1.6601834628507000e+04, new int[] { 3, 7, 1 });
            p.AddCoeff(-1.7524158774535100e+04, new int[] { 3, 9, 1 });
            p.AddCoeff(6.6910424411861400e+03, new int[] { 3, 11, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a6f23b87-a490-4e02-8d9a-22e6b4c8c6f1"));
            OrthonormalPolynomials[737] = p;
            p.AddCoeff(-1.5826224176235000e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(1.2344454857463300e+02, new int[] { 1, 2, 0 });
            p.AddCoeff(-1.5430568571829100e+03, new int[] { 1, 4, 0 });
            p.AddCoeff(6.9951910858958600e+03, new int[] { 1, 6, 0 });
            p.AddCoeff(-1.4240210424859400e+04, new int[] { 1, 8, 0 });
            p.AddCoeff(1.3290863063202200e+04, new int[] { 1, 10, 0 });
            p.AddCoeff(-4.6316644008128600e+03, new int[] { 1, 12, 0 });
            p.AddCoeff(2.6377040293725000e+00, new int[] { 3, 0, 0 });
            p.AddCoeff(-2.0574091429105500e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(2.5717614286381800e+03, new int[] { 3, 4, 0 });
            p.AddCoeff(-1.1658651809826400e+04, new int[] { 3, 6, 0 });
            p.AddCoeff(2.3733684041432400e+04, new int[] { 3, 8, 0 });
            p.AddCoeff(-2.2151438438670200e+04, new int[] { 3, 10, 0 });
            p.AddCoeff(7.7194406680214600e+03, new int[] { 3, 12, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a109c30f-49d6-4512-82ee-73a70640157e"));
            OrthonormalPolynomials[738] = p;
            p.AddCoeff(-5.1637441534121500e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(1.1188112332393000e+02, new int[] { 0, 0, 3 });
            p.AddCoeff(-6.7128673994357900e+02, new int[] { 0, 0, 5 });
            p.AddCoeff(1.6302677970058400e+03, new int[] { 0, 0, 7 });
            p.AddCoeff(-1.7208382301728200e+03, new int[] { 0, 0, 9 });
            p.AddCoeff(6.5704732424780600e+02, new int[] { 0, 0, 11 });
            p.AddCoeff(5.1637441534121500e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-1.1188112332393000e+03, new int[] { 2, 0, 3 });
            p.AddCoeff(6.7128673994357900e+03, new int[] { 2, 0, 5 });
            p.AddCoeff(-1.6302677970058400e+04, new int[] { 2, 0, 7 });
            p.AddCoeff(1.7208382301728200e+04, new int[] { 2, 0, 9 });
            p.AddCoeff(-6.5704732424780600e+03, new int[] { 2, 0, 11 });
            p.AddCoeff(-6.0243681789808400e+01, new int[] { 4, 0, 1 });
            p.AddCoeff(1.3052797721125200e+03, new int[] { 4, 0, 3 });
            p.AddCoeff(-7.8316786326750600e+03, new int[] { 4, 0, 5 });
            p.AddCoeff(1.9019790965068100e+04, new int[] { 4, 0, 7 });
            p.AddCoeff(-2.0076446018683000e+04, new int[] { 4, 0, 9 });
            p.AddCoeff(7.6655521162244200e+03, new int[] { 4, 0, 11 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("07054184-18b0-4581-89b6-f67ae46f72f8"));
            OrthonormalPolynomials[739] = p;
            p.AddCoeff(-7.7692373228789900e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(4.2730805275834500e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-3.7033364572389900e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(1.1110009371717000e+03, new int[] { 0, 1, 6 });
            p.AddCoeff(-1.3490725665656300e+03, new int[] { 0, 1, 8 });
            p.AddCoeff(5.6960841699437900e+02, new int[] { 0, 1, 10 });
            p.AddCoeff(7.7692373228789900e+00, new int[] { 2, 1, 0 });
            p.AddCoeff(-4.2730805275834500e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(3.7033364572389900e+03, new int[] { 2, 1, 4 });
            p.AddCoeff(-1.1110009371717000e+04, new int[] { 2, 1, 6 });
            p.AddCoeff(1.3490725665656300e+04, new int[] { 2, 1, 8 });
            p.AddCoeff(-5.6960841699437900e+03, new int[] { 2, 1, 10 });
            p.AddCoeff(-9.0641102100254900e+00, new int[] { 4, 1, 0 });
            p.AddCoeff(4.9852606155140400e+02, new int[] { 4, 1, 2 });
            p.AddCoeff(-4.3205592001121300e+03, new int[] { 4, 1, 4 });
            p.AddCoeff(1.2961677600336500e+04, new int[] { 4, 1, 6 });
            p.AddCoeff(-1.5739179943265700e+04, new int[] { 4, 1, 8 });
            p.AddCoeff(6.6454315316010800e+03, new int[] { 4, 1, 10 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("6702eb48-374c-4506-8d54-2ef0ac271f6b"));
            OrthonormalPolynomials[740] = p;
            p.AddCoeff(-4.7702365981659800e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(6.9963470106434400e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-2.7285753341509400e+02, new int[] { 0, 0, 5 });
            p.AddCoeff(3.8979647630727700e+02, new int[] { 0, 0, 7 });
            p.AddCoeff(-1.8407055825621400e+02, new int[] { 0, 0, 9 });
            p.AddCoeff(1.4310709794497900e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-2.0989041031930300e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(8.1857260024528200e+02, new int[] { 0, 2, 5 });
            p.AddCoeff(-1.1693894289218300e+03, new int[] { 0, 2, 7 });
            p.AddCoeff(5.5221167476864300e+02, new int[] { 0, 2, 9 });
            p.AddCoeff(4.7702365981659800e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-6.9963470106434400e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(2.7285753341509400e+03, new int[] { 2, 0, 5 });
            p.AddCoeff(-3.8979647630727700e+03, new int[] { 2, 0, 7 });
            p.AddCoeff(1.8407055825621400e+03, new int[] { 2, 0, 9 });
            p.AddCoeff(-1.4310709794497900e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(2.0989041031930300e+03, new int[] { 2, 2, 3 });
            p.AddCoeff(-8.1857260024528200e+03, new int[] { 2, 2, 5 });
            p.AddCoeff(1.1693894289218300e+04, new int[] { 2, 2, 7 });
            p.AddCoeff(-5.5221167476864300e+03, new int[] { 2, 2, 9 });
            p.AddCoeff(-5.5652760311936400e+01, new int[] { 4, 0, 1 });
            p.AddCoeff(8.1624048457506800e+02, new int[] { 4, 0, 3 });
            p.AddCoeff(-3.1833378898427600e+03, new int[] { 4, 0, 5 });
            p.AddCoeff(4.5476255569182400e+03, new int[] { 4, 0, 7 });
            p.AddCoeff(-2.1474898463225100e+03, new int[] { 4, 0, 9 });
            p.AddCoeff(1.6695828093580900e+02, new int[] { 4, 2, 1 });
            p.AddCoeff(-2.4487214537252000e+03, new int[] { 4, 2, 3 });
            p.AddCoeff(9.5500136695282900e+03, new int[] { 4, 2, 5 });
            p.AddCoeff(-1.3642876670754700e+04, new int[] { 4, 2, 7 });
            p.AddCoeff(6.4424695389674900e+03, new int[] { 4, 2, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("dc03b05d-2349-4fdd-b439-a791931f2829"));
            OrthonormalPolynomials[741] = p;
            p.AddCoeff(-1.7796325618178400e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(6.4066772225442100e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-3.5236724723993200e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(6.1076989521588200e+02, new int[] { 0, 1, 6 });
            p.AddCoeff(-3.2719815815136500e+02, new int[] { 0, 1, 8 });
            p.AddCoeff(2.9660542696963900e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(-1.0677795370907000e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(5.8727874539988600e+02, new int[] { 0, 3, 4 });
            p.AddCoeff(-1.0179498253598000e+03, new int[] { 0, 3, 6 });
            p.AddCoeff(5.4533026358560900e+02, new int[] { 0, 3, 8 });
            p.AddCoeff(1.7796325618178400e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-6.4066772225442100e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(3.5236724723993200e+03, new int[] { 2, 1, 4 });
            p.AddCoeff(-6.1076989521588200e+03, new int[] { 2, 1, 6 });
            p.AddCoeff(3.2719815815136500e+03, new int[] { 2, 1, 8 });
            p.AddCoeff(-2.9660542696963900e+01, new int[] { 2, 3, 0 });
            p.AddCoeff(1.0677795370907000e+03, new int[] { 2, 3, 2 });
            p.AddCoeff(-5.8727874539988600e+03, new int[] { 2, 3, 4 });
            p.AddCoeff(1.0179498253598000e+04, new int[] { 2, 3, 6 });
            p.AddCoeff(-5.4533026358560900e+03, new int[] { 2, 3, 8 });
            p.AddCoeff(-2.0762379887874800e+01, new int[] { 4, 1, 0 });
            p.AddCoeff(7.4744567596349100e+02, new int[] { 4, 1, 2 });
            p.AddCoeff(-4.1109512177992000e+03, new int[] { 4, 1, 4 });
            p.AddCoeff(7.1256487775186200e+03, new int[] { 4, 1, 6 });
            p.AddCoeff(-3.8173118450992500e+03, new int[] { 4, 1, 8 });
            p.AddCoeff(3.4603966479791300e+01, new int[] { 4, 3, 0 });
            p.AddCoeff(-1.2457427932724900e+03, new int[] { 4, 3, 2 });
            p.AddCoeff(6.8515853629986700e+03, new int[] { 4, 3, 4 });
            p.AddCoeff(-1.1876081295864400e+04, new int[] { 4, 3, 6 });
            p.AddCoeff(6.3621864084987700e+03, new int[] { 4, 3, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a0d3d723-f236-4229-af2f-c2588dd09d60"));
            OrthonormalPolynomials[742] = p;
            p.AddCoeff(-3.7909996350760400e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(3.4118996715684400e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-7.5061792774505600e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(4.6466824098503500e+01, new int[] { 0, 0, 7 });
            p.AddCoeff(3.7909996350760400e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-3.4118996715684400e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(7.5061792774505600e+02, new int[] { 0, 2, 5 });
            p.AddCoeff(-4.6466824098503500e+02, new int[] { 0, 2, 7 });
            p.AddCoeff(-4.4228329075887100e+01, new int[] { 0, 4, 1 });
            p.AddCoeff(3.9805496168298400e+02, new int[] { 0, 4, 3 });
            p.AddCoeff(-8.7572091570256500e+02, new int[] { 0, 4, 5 });
            p.AddCoeff(5.4211294781587400e+02, new int[] { 0, 4, 7 });
            p.AddCoeff(3.7909996350760400e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-3.4118996715684400e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(7.5061792774505600e+02, new int[] { 2, 0, 5 });
            p.AddCoeff(-4.6466824098503500e+02, new int[] { 2, 0, 7 });
            p.AddCoeff(-3.7909996350760400e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(3.4118996715684400e+03, new int[] { 2, 2, 3 });
            p.AddCoeff(-7.5061792774505600e+03, new int[] { 2, 2, 5 });
            p.AddCoeff(4.6466824098503500e+03, new int[] { 2, 2, 7 });
            p.AddCoeff(4.4228329075887100e+02, new int[] { 2, 4, 1 });
            p.AddCoeff(-3.9805496168298400e+03, new int[] { 2, 4, 3 });
            p.AddCoeff(8.7572091570256500e+03, new int[] { 2, 4, 5 });
            p.AddCoeff(-5.4211294781587400e+03, new int[] { 2, 4, 7 });
            p.AddCoeff(-4.4228329075887100e+01, new int[] { 4, 0, 1 });
            p.AddCoeff(3.9805496168298400e+02, new int[] { 4, 0, 3 });
            p.AddCoeff(-8.7572091570256500e+02, new int[] { 4, 0, 5 });
            p.AddCoeff(5.4211294781587400e+02, new int[] { 4, 0, 7 });
            p.AddCoeff(4.4228329075887100e+02, new int[] { 4, 2, 1 });
            p.AddCoeff(-3.9805496168298400e+03, new int[] { 4, 2, 3 });
            p.AddCoeff(8.7572091570256500e+03, new int[] { 4, 2, 5 });
            p.AddCoeff(-5.4211294781587400e+03, new int[] { 4, 2, 7 });
            p.AddCoeff(-5.1599717255201600e+02, new int[] { 4, 4, 1 });
            p.AddCoeff(4.6439745529681500e+03, new int[] { 4, 4, 3 });
            p.AddCoeff(-1.0216744016529900e+04, new int[] { 4, 4, 5 });
            p.AddCoeff(6.3246510578518300e+03, new int[] { 4, 4, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("bc78428a-8332-463a-b170-800fea269c7f"));
            OrthonormalPolynomials[743] = p;
            p.AddCoeff(-2.7869350108811700e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(5.8525635228504500e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-1.7557690568551400e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(1.2875639750271000e+02, new int[] { 0, 1, 6 });
            p.AddCoeff(1.3005696717445500e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-2.7311963106635400e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(8.1935889319906300e+02, new int[] { 0, 3, 4 });
            p.AddCoeff(-6.0086318834598000e+02, new int[] { 0, 3, 6 });
            p.AddCoeff(-1.1705127045700900e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(2.4580766795971900e+02, new int[] { 0, 5, 2 });
            p.AddCoeff(-7.3742300387915700e+02, new int[] { 0, 5, 4 });
            p.AddCoeff(5.4077686951138200e+02, new int[] { 0, 5, 6 });
            p.AddCoeff(2.7869350108811700e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-5.8525635228504500e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(1.7557690568551400e+03, new int[] { 2, 1, 4 });
            p.AddCoeff(-1.2875639750271000e+03, new int[] { 2, 1, 6 });
            p.AddCoeff(-1.3005696717445500e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(2.7311963106635400e+03, new int[] { 2, 3, 2 });
            p.AddCoeff(-8.1935889319906300e+03, new int[] { 2, 3, 4 });
            p.AddCoeff(6.0086318834598000e+03, new int[] { 2, 3, 6 });
            p.AddCoeff(1.1705127045700900e+02, new int[] { 2, 5, 0 });
            p.AddCoeff(-2.4580766795971900e+03, new int[] { 2, 5, 2 });
            p.AddCoeff(7.3742300387915700e+03, new int[] { 2, 5, 4 });
            p.AddCoeff(-5.4077686951138200e+03, new int[] { 2, 5, 6 });
            p.AddCoeff(-3.2514241793613600e+01, new int[] { 4, 1, 0 });
            p.AddCoeff(6.8279907766588600e+02, new int[] { 4, 1, 2 });
            p.AddCoeff(-2.0483972329976600e+03, new int[] { 4, 1, 4 });
            p.AddCoeff(1.5021579708649500e+03, new int[] { 4, 1, 6 });
            p.AddCoeff(1.5173312837019700e+02, new int[] { 4, 3, 0 });
            p.AddCoeff(-3.1863956957741400e+03, new int[] { 4, 3, 2 });
            p.AddCoeff(9.5591870873224100e+03, new int[] { 4, 3, 4 });
            p.AddCoeff(-7.0100705307031000e+03, new int[] { 4, 3, 6 });
            p.AddCoeff(-1.3655981553317700e+02, new int[] { 4, 5, 0 });
            p.AddCoeff(2.8677561261967200e+03, new int[] { 4, 5, 2 });
            p.AddCoeff(-8.6032683785901700e+03, new int[] { 4, 5, 4 });
            p.AddCoeff(6.3090634776327900e+03, new int[] { 4, 5, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("9e0237dc-99f6-4ddb-bc90-66374ac9fff9"));
            OrthonormalPolynomials[744] = p;
            p.AddCoeff(-2.7869350108811700e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(1.3005696717445500e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-1.1705127045700900e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(5.8525635228504500e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-2.7311963106635400e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(2.4580766795971900e+02, new int[] { 0, 2, 5 });
            p.AddCoeff(-1.7557690568551400e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(8.1935889319906300e+02, new int[] { 0, 4, 3 });
            p.AddCoeff(-7.3742300387915700e+02, new int[] { 0, 4, 5 });
            p.AddCoeff(1.2875639750271000e+02, new int[] { 0, 6, 1 });
            p.AddCoeff(-6.0086318834598000e+02, new int[] { 0, 6, 3 });
            p.AddCoeff(5.4077686951138200e+02, new int[] { 0, 6, 5 });
            p.AddCoeff(2.7869350108811700e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-1.3005696717445500e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(1.1705127045700900e+02, new int[] { 2, 0, 5 });
            p.AddCoeff(-5.8525635228504500e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(2.7311963106635400e+03, new int[] { 2, 2, 3 });
            p.AddCoeff(-2.4580766795971900e+03, new int[] { 2, 2, 5 });
            p.AddCoeff(1.7557690568551400e+03, new int[] { 2, 4, 1 });
            p.AddCoeff(-8.1935889319906300e+03, new int[] { 2, 4, 3 });
            p.AddCoeff(7.3742300387915700e+03, new int[] { 2, 4, 5 });
            p.AddCoeff(-1.2875639750271000e+03, new int[] { 2, 6, 1 });
            p.AddCoeff(6.0086318834598000e+03, new int[] { 2, 6, 3 });
            p.AddCoeff(-5.4077686951138200e+03, new int[] { 2, 6, 5 });
            p.AddCoeff(-3.2514241793613600e+01, new int[] { 4, 0, 1 });
            p.AddCoeff(1.5173312837019700e+02, new int[] { 4, 0, 3 });
            p.AddCoeff(-1.3655981553317700e+02, new int[] { 4, 0, 5 });
            p.AddCoeff(6.8279907766588600e+02, new int[] { 4, 2, 1 });
            p.AddCoeff(-3.1863956957741400e+03, new int[] { 4, 2, 3 });
            p.AddCoeff(2.8677561261967200e+03, new int[] { 4, 2, 5 });
            p.AddCoeff(-2.0483972329976600e+03, new int[] { 4, 4, 1 });
            p.AddCoeff(9.5591870873224100e+03, new int[] { 4, 4, 3 });
            p.AddCoeff(-8.6032683785901700e+03, new int[] { 4, 4, 5 });
            p.AddCoeff(1.5021579708649500e+03, new int[] { 4, 6, 1 });
            p.AddCoeff(-7.0100705307031000e+03, new int[] { 4, 6, 3 });
            p.AddCoeff(6.3090634776327900e+03, new int[] { 4, 6, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7cebbff5-757e-4d9e-bf8d-3f1706cd3124"));
            OrthonormalPolynomials[745] = p;
            p.AddCoeff(-3.7909996350760400e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(3.7909996350760400e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-4.4228329075887100e+01, new int[] { 0, 1, 4 });
            p.AddCoeff(3.4118996715684400e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-3.4118996715684400e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(3.9805496168298400e+02, new int[] { 0, 3, 4 });
            p.AddCoeff(-7.5061792774505600e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(7.5061792774505600e+02, new int[] { 0, 5, 2 });
            p.AddCoeff(-8.7572091570256500e+02, new int[] { 0, 5, 4 });
            p.AddCoeff(4.6466824098503500e+01, new int[] { 0, 7, 0 });
            p.AddCoeff(-4.6466824098503500e+02, new int[] { 0, 7, 2 });
            p.AddCoeff(5.4211294781587400e+02, new int[] { 0, 7, 4 });
            p.AddCoeff(3.7909996350760400e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-3.7909996350760400e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(4.4228329075887100e+02, new int[] { 2, 1, 4 });
            p.AddCoeff(-3.4118996715684400e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(3.4118996715684400e+03, new int[] { 2, 3, 2 });
            p.AddCoeff(-3.9805496168298400e+03, new int[] { 2, 3, 4 });
            p.AddCoeff(7.5061792774505600e+02, new int[] { 2, 5, 0 });
            p.AddCoeff(-7.5061792774505600e+03, new int[] { 2, 5, 2 });
            p.AddCoeff(8.7572091570256500e+03, new int[] { 2, 5, 4 });
            p.AddCoeff(-4.6466824098503500e+02, new int[] { 2, 7, 0 });
            p.AddCoeff(4.6466824098503500e+03, new int[] { 2, 7, 2 });
            p.AddCoeff(-5.4211294781587400e+03, new int[] { 2, 7, 4 });
            p.AddCoeff(-4.4228329075887100e+01, new int[] { 4, 1, 0 });
            p.AddCoeff(4.4228329075887100e+02, new int[] { 4, 1, 2 });
            p.AddCoeff(-5.1599717255201600e+02, new int[] { 4, 1, 4 });
            p.AddCoeff(3.9805496168298400e+02, new int[] { 4, 3, 0 });
            p.AddCoeff(-3.9805496168298400e+03, new int[] { 4, 3, 2 });
            p.AddCoeff(4.6439745529681500e+03, new int[] { 4, 3, 4 });
            p.AddCoeff(-8.7572091570256500e+02, new int[] { 4, 5, 0 });
            p.AddCoeff(8.7572091570256500e+03, new int[] { 4, 5, 2 });
            p.AddCoeff(-1.0216744016529900e+04, new int[] { 4, 5, 4 });
            p.AddCoeff(5.4211294781587400e+02, new int[] { 4, 7, 0 });
            p.AddCoeff(-5.4211294781587400e+03, new int[] { 4, 7, 2 });
            p.AddCoeff(6.3246510578518300e+03, new int[] { 4, 7, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("55693cf2-d4c8-4a92-af92-6052cdcfad54"));
            OrthonormalPolynomials[746] = p;
            p.AddCoeff(-1.7796325618178400e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(2.9660542696963900e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(6.4066772225442100e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-1.0677795370907000e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(-3.5236724723993200e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(5.8727874539988600e+02, new int[] { 0, 4, 3 });
            p.AddCoeff(6.1076989521588200e+02, new int[] { 0, 6, 1 });
            p.AddCoeff(-1.0179498253598000e+03, new int[] { 0, 6, 3 });
            p.AddCoeff(-3.2719815815136500e+02, new int[] { 0, 8, 1 });
            p.AddCoeff(5.4533026358560900e+02, new int[] { 0, 8, 3 });
            p.AddCoeff(1.7796325618178400e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-2.9660542696963900e+01, new int[] { 2, 0, 3 });
            p.AddCoeff(-6.4066772225442100e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(1.0677795370907000e+03, new int[] { 2, 2, 3 });
            p.AddCoeff(3.5236724723993200e+03, new int[] { 2, 4, 1 });
            p.AddCoeff(-5.8727874539988600e+03, new int[] { 2, 4, 3 });
            p.AddCoeff(-6.1076989521588200e+03, new int[] { 2, 6, 1 });
            p.AddCoeff(1.0179498253598000e+04, new int[] { 2, 6, 3 });
            p.AddCoeff(3.2719815815136500e+03, new int[] { 2, 8, 1 });
            p.AddCoeff(-5.4533026358560900e+03, new int[] { 2, 8, 3 });
            p.AddCoeff(-2.0762379887874800e+01, new int[] { 4, 0, 1 });
            p.AddCoeff(3.4603966479791300e+01, new int[] { 4, 0, 3 });
            p.AddCoeff(7.4744567596349100e+02, new int[] { 4, 2, 1 });
            p.AddCoeff(-1.2457427932724900e+03, new int[] { 4, 2, 3 });
            p.AddCoeff(-4.1109512177992000e+03, new int[] { 4, 4, 1 });
            p.AddCoeff(6.8515853629986700e+03, new int[] { 4, 4, 3 });
            p.AddCoeff(7.1256487775186200e+03, new int[] { 4, 6, 1 });
            p.AddCoeff(-1.1876081295864400e+04, new int[] { 4, 6, 3 });
            p.AddCoeff(-3.8173118450992500e+03, new int[] { 4, 8, 1 });
            p.AddCoeff(6.3621864084987700e+03, new int[] { 4, 8, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("68ae0b64-617a-4377-b544-5819643f37d0"));
            OrthonormalPolynomials[747] = p;
            p.AddCoeff(-4.7702365981659800e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(1.4310709794497900e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(6.9963470106434400e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-2.0989041031930300e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(-2.7285753341509400e+02, new int[] { 0, 5, 0 });
            p.AddCoeff(8.1857260024528200e+02, new int[] { 0, 5, 2 });
            p.AddCoeff(3.8979647630727700e+02, new int[] { 0, 7, 0 });
            p.AddCoeff(-1.1693894289218300e+03, new int[] { 0, 7, 2 });
            p.AddCoeff(-1.8407055825621400e+02, new int[] { 0, 9, 0 });
            p.AddCoeff(5.5221167476864300e+02, new int[] { 0, 9, 2 });
            p.AddCoeff(4.7702365981659800e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-1.4310709794497900e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(-6.9963470106434400e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(2.0989041031930300e+03, new int[] { 2, 3, 2 });
            p.AddCoeff(2.7285753341509400e+03, new int[] { 2, 5, 0 });
            p.AddCoeff(-8.1857260024528200e+03, new int[] { 2, 5, 2 });
            p.AddCoeff(-3.8979647630727700e+03, new int[] { 2, 7, 0 });
            p.AddCoeff(1.1693894289218300e+04, new int[] { 2, 7, 2 });
            p.AddCoeff(1.8407055825621400e+03, new int[] { 2, 9, 0 });
            p.AddCoeff(-5.5221167476864300e+03, new int[] { 2, 9, 2 });
            p.AddCoeff(-5.5652760311936400e+01, new int[] { 4, 1, 0 });
            p.AddCoeff(1.6695828093580900e+02, new int[] { 4, 1, 2 });
            p.AddCoeff(8.1624048457506800e+02, new int[] { 4, 3, 0 });
            p.AddCoeff(-2.4487214537252000e+03, new int[] { 4, 3, 2 });
            p.AddCoeff(-3.1833378898427600e+03, new int[] { 4, 5, 0 });
            p.AddCoeff(9.5500136695282900e+03, new int[] { 4, 5, 2 });
            p.AddCoeff(4.5476255569182400e+03, new int[] { 4, 7, 0 });
            p.AddCoeff(-1.3642876670754700e+04, new int[] { 4, 7, 2 });
            p.AddCoeff(-2.1474898463225100e+03, new int[] { 4, 9, 0 });
            p.AddCoeff(6.4424695389674900e+03, new int[] { 4, 9, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b5638dad-9715-44d0-94a0-f674e4dd935a"));
            OrthonormalPolynomials[748] = p;
            p.AddCoeff(-7.7692373228789900e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(4.2730805275834500e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-3.7033364572389900e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(1.1110009371717000e+03, new int[] { 0, 6, 1 });
            p.AddCoeff(-1.3490725665656300e+03, new int[] { 0, 8, 1 });
            p.AddCoeff(5.6960841699437900e+02, new int[] { 0, 10, 1 });
            p.AddCoeff(7.7692373228789900e+00, new int[] { 2, 0, 1 });
            p.AddCoeff(-4.2730805275834500e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(3.7033364572389900e+03, new int[] { 2, 4, 1 });
            p.AddCoeff(-1.1110009371717000e+04, new int[] { 2, 6, 1 });
            p.AddCoeff(1.3490725665656300e+04, new int[] { 2, 8, 1 });
            p.AddCoeff(-5.6960841699437900e+03, new int[] { 2, 10, 1 });
            p.AddCoeff(-9.0641102100254900e+00, new int[] { 4, 0, 1 });
            p.AddCoeff(4.9852606155140400e+02, new int[] { 4, 2, 1 });
            p.AddCoeff(-4.3205592001121300e+03, new int[] { 4, 4, 1 });
            p.AddCoeff(1.2961677600336500e+04, new int[] { 4, 6, 1 });
            p.AddCoeff(-1.5739179943265700e+04, new int[] { 4, 8, 1 });
            p.AddCoeff(6.6454315316010800e+03, new int[] { 4, 10, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c2f2cbe0-097e-4c8e-8313-6719963781f0"));
            OrthonormalPolynomials[749] = p;
            p.AddCoeff(-5.1637441534121500e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(1.1188112332393000e+02, new int[] { 0, 3, 0 });
            p.AddCoeff(-6.7128673994357900e+02, new int[] { 0, 5, 0 });
            p.AddCoeff(1.6302677970058400e+03, new int[] { 0, 7, 0 });
            p.AddCoeff(-1.7208382301728200e+03, new int[] { 0, 9, 0 });
            p.AddCoeff(6.5704732424780600e+02, new int[] { 0, 11, 0 });
            p.AddCoeff(5.1637441534121500e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-1.1188112332393000e+03, new int[] { 2, 3, 0 });
            p.AddCoeff(6.7128673994357900e+03, new int[] { 2, 5, 0 });
            p.AddCoeff(-1.6302677970058400e+04, new int[] { 2, 7, 0 });
            p.AddCoeff(1.7208382301728200e+04, new int[] { 2, 9, 0 });
            p.AddCoeff(-6.5704732424780600e+03, new int[] { 2, 11, 0 });
            p.AddCoeff(-6.0243681789808400e+01, new int[] { 4, 1, 0 });
            p.AddCoeff(1.3052797721125200e+03, new int[] { 4, 3, 0 });
            p.AddCoeff(-7.8316786326750600e+03, new int[] { 4, 5, 0 });
            p.AddCoeff(1.9019790965068100e+04, new int[] { 4, 7, 0 });
            p.AddCoeff(-2.0076446018683000e+04, new int[] { 4, 9, 0 });
            p.AddCoeff(7.6655521162244200e+03, new int[] { 4, 11, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("084d3a5a-3026-4d56-aa55-7902166a184a"));
            OrthonormalPolynomials[750] = p;
            p.AddCoeff(-2.4794928065055500e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(1.3637210435780500e+02, new int[] { 1, 0, 2 });
            p.AddCoeff(-1.1818915711009800e+03, new int[] { 1, 0, 4 });
            p.AddCoeff(3.5456747133029300e+03, new int[] { 1, 0, 6 });
            p.AddCoeff(-4.3054621518678400e+03, new int[] { 1, 0, 8 });
            p.AddCoeff(1.8178617974553100e+03, new int[] { 1, 0, 10 });
            p.AddCoeff(1.1570966430359200e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-6.3640315366975700e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(5.5154939984712300e+03, new int[] { 3, 0, 4 });
            p.AddCoeff(-1.6546481995413700e+04, new int[] { 3, 0, 6 });
            p.AddCoeff(2.0092156708716600e+04, new int[] { 3, 0, 8 });
            p.AddCoeff(-8.4833550547914600e+03, new int[] { 3, 0, 10 });
            p.AddCoeff(-1.0413869787323300e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(5.7276283830278100e+02, new int[] { 5, 0, 2 });
            p.AddCoeff(-4.9639445986241100e+03, new int[] { 5, 0, 4 });
            p.AddCoeff(1.4891833795872300e+04, new int[] { 5, 0, 6 });
            p.AddCoeff(-1.8082941037844900e+04, new int[] { 5, 0, 8 });
            p.AddCoeff(7.6350195493123200e+03, new int[] { 5, 0, 10 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d7544135-3f77-4c95-90e5-6bdc78adaf46"));
            OrthonormalPolynomials[751] = p;
            p.AddCoeff(4.0849865705801600e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-5.9913136368509000e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(2.3366123183718500e+03, new int[] { 1, 1, 5 });
            p.AddCoeff(-3.3380175976740700e+03, new int[] { 1, 1, 7 });
            p.AddCoeff(1.5762860877905300e+03, new int[] { 1, 1, 9 });
            p.AddCoeff(-1.9063270662707400e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(2.7959463638637500e+03, new int[] { 3, 1, 3 });
            p.AddCoeff(-1.0904190819068600e+04, new int[] { 3, 1, 5 });
            p.AddCoeff(1.5577415455812300e+04, new int[] { 3, 1, 7 });
            p.AddCoeff(-7.3560017430224900e+03, new int[] { 3, 1, 9 });
            p.AddCoeff(1.7156943596436700e+02, new int[] { 5, 1, 1 });
            p.AddCoeff(-2.5163517274773800e+03, new int[] { 5, 1, 3 });
            p.AddCoeff(9.8137717371617800e+03, new int[] { 5, 1, 5 });
            p.AddCoeff(-1.4019673910231100e+04, new int[] { 5, 1, 7 });
            p.AddCoeff(6.6204015687202500e+03, new int[] { 5, 1, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7f4f2bb1-5f7b-4f32-a3ca-d001a2bea747"));
            OrthonormalPolynomials[752] = p;
            p.AddCoeff(-2.7713422517043000e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(9.9768321061354800e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-5.4872576583745200e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(9.5112466078491600e+02, new int[] { 1, 0, 6 });
            p.AddCoeff(-5.0953106827763400e+02, new int[] { 1, 0, 8 });
            p.AddCoeff(8.3140267551129000e+00, new int[] { 1, 2, 0 });
            p.AddCoeff(-2.9930496318406400e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(1.6461772975123500e+03, new int[] { 1, 2, 4 });
            p.AddCoeff(-2.8533739823547500e+03, new int[] { 1, 2, 6 });
            p.AddCoeff(1.5285932048329000e+03, new int[] { 1, 2, 8 });
            p.AddCoeff(1.2932930507953400e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-4.6558549828632300e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(2.5607202405747700e+03, new int[] { 3, 0, 4 });
            p.AddCoeff(-4.4385817503296100e+03, new int[] { 3, 0, 6 });
            p.AddCoeff(2.3778116519622900e+03, new int[] { 3, 0, 8 });
            p.AddCoeff(-3.8798791523860200e+01, new int[] { 3, 2, 0 });
            p.AddCoeff(1.3967564948589700e+03, new int[] { 3, 2, 2 });
            p.AddCoeff(-7.6821607217243200e+03, new int[] { 3, 2, 4 });
            p.AddCoeff(1.3315745250988800e+04, new int[] { 3, 2, 6 });
            p.AddCoeff(-7.1334349558868700e+03, new int[] { 3, 2, 8 });
            p.AddCoeff(-1.1639637457158100e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(4.1902694845769000e+02, new int[] { 5, 0, 2 });
            p.AddCoeff(-2.3046482165173000e+03, new int[] { 5, 0, 4 });
            p.AddCoeff(3.9947235752966500e+03, new int[] { 5, 0, 6 });
            p.AddCoeff(-2.1400304867660600e+03, new int[] { 5, 0, 8 });
            p.AddCoeff(3.4918912371474200e+01, new int[] { 5, 2, 0 });
            p.AddCoeff(-1.2570808453730700e+03, new int[] { 5, 2, 2 });
            p.AddCoeff(6.9139446495518900e+03, new int[] { 5, 2, 4 });
            p.AddCoeff(-1.1984170725889900e+04, new int[] { 5, 2, 6 });
            p.AddCoeff(6.4200914602982000e+03, new int[] { 5, 2, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7ed60a09-56ba-439d-9b3e-15b421bc285c"));
            OrthonormalPolynomials[753] = p;
            p.AddCoeff(7.3924192867575100e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-6.6531773580817600e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(1.4636990187779900e+03, new int[] { 1, 1, 5 });
            p.AddCoeff(-9.0609939257684900e+02, new int[] { 1, 1, 7 });
            p.AddCoeff(-1.2320698811262500e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(1.1088628930136300e+03, new int[] { 1, 3, 3 });
            p.AddCoeff(-2.4394983646299800e+03, new int[] { 1, 3, 5 });
            p.AddCoeff(1.5101656542947500e+03, new int[] { 1, 3, 7 });
            p.AddCoeff(-3.4497956671535100e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(3.1048161004381600e+03, new int[] { 3, 1, 3 });
            p.AddCoeff(-6.8305954209639400e+03, new int[] { 3, 1, 5 });
            p.AddCoeff(4.2284638320253000e+03, new int[] { 3, 1, 7 });
            p.AddCoeff(5.7496594452558400e+02, new int[] { 3, 3, 1 });
            p.AddCoeff(-5.1746935007302600e+03, new int[] { 3, 3, 3 });
            p.AddCoeff(1.1384325701606600e+04, new int[] { 3, 3, 5 });
            p.AddCoeff(-7.0474397200421600e+03, new int[] { 3, 3, 7 });
            p.AddCoeff(3.1048161004381600e+02, new int[] { 5, 1, 1 });
            p.AddCoeff(-2.7943344903943400e+03, new int[] { 5, 1, 3 });
            p.AddCoeff(6.1475358788675500e+03, new int[] { 5, 1, 5 });
            p.AddCoeff(-3.8056174488227700e+03, new int[] { 5, 1, 7 });
            p.AddCoeff(-5.1746935007302600e+02, new int[] { 5, 3, 1 });
            p.AddCoeff(4.6572241506572300e+03, new int[] { 5, 3, 3 });
            p.AddCoeff(-1.0245893131445900e+04, new int[] { 5, 3, 5 });
            p.AddCoeff(6.3426957480379500e+03, new int[] { 5, 3, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f52e5486-ab08-4b78-b3ab-86a558183015"));
            OrthonormalPolynomials[754] = p;
            p.AddCoeff(-2.7869350108811700e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(5.8525635228504500e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-1.7557690568551400e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(1.2875639750271000e+02, new int[] { 1, 0, 6 });
            p.AddCoeff(2.7869350108811700e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-5.8525635228504500e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(1.7557690568551400e+03, new int[] { 1, 2, 4 });
            p.AddCoeff(-1.2875639750271000e+03, new int[] { 1, 2, 6 });
            p.AddCoeff(-3.2514241793613600e+01, new int[] { 1, 4, 0 });
            p.AddCoeff(6.8279907766588600e+02, new int[] { 1, 4, 2 });
            p.AddCoeff(-2.0483972329976600e+03, new int[] { 1, 4, 4 });
            p.AddCoeff(1.5021579708649500e+03, new int[] { 1, 4, 6 });
            p.AddCoeff(1.3005696717445500e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-2.7311963106635400e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(8.1935889319906300e+02, new int[] { 3, 0, 4 });
            p.AddCoeff(-6.0086318834598000e+02, new int[] { 3, 0, 6 });
            p.AddCoeff(-1.3005696717445500e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(2.7311963106635400e+03, new int[] { 3, 2, 2 });
            p.AddCoeff(-8.1935889319906300e+03, new int[] { 3, 2, 4 });
            p.AddCoeff(6.0086318834598000e+03, new int[] { 3, 2, 6 });
            p.AddCoeff(1.5173312837019700e+02, new int[] { 3, 4, 0 });
            p.AddCoeff(-3.1863956957741400e+03, new int[] { 3, 4, 2 });
            p.AddCoeff(9.5591870873224100e+03, new int[] { 3, 4, 4 });
            p.AddCoeff(-7.0100705307031000e+03, new int[] { 3, 4, 6 });
            p.AddCoeff(-1.1705127045700900e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(2.4580766795971900e+02, new int[] { 5, 0, 2 });
            p.AddCoeff(-7.3742300387915700e+02, new int[] { 5, 0, 4 });
            p.AddCoeff(5.4077686951138200e+02, new int[] { 5, 0, 6 });
            p.AddCoeff(1.1705127045700900e+02, new int[] { 5, 2, 0 });
            p.AddCoeff(-2.4580766795971900e+03, new int[] { 5, 2, 2 });
            p.AddCoeff(7.3742300387915700e+03, new int[] { 5, 2, 4 });
            p.AddCoeff(-5.4077686951138200e+03, new int[] { 5, 2, 6 });
            p.AddCoeff(-1.3655981553317700e+02, new int[] { 5, 4, 0 });
            p.AddCoeff(2.8677561261967200e+03, new int[] { 5, 4, 2 });
            p.AddCoeff(-8.6032683785901700e+03, new int[] { 5, 4, 4 });
            p.AddCoeff(6.3090634776327900e+03, new int[] { 5, 4, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("dc894727-46fb-46ef-bb92-96f8c557965c"));
            OrthonormalPolynomials[755] = p;
            p.AddCoeff(8.5025236857150800e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-3.9678443866670400e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(3.5710599480003300e+02, new int[] { 1, 1, 5 });
            p.AddCoeff(-3.9678443866670400e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(1.8516607137779500e+03, new int[] { 1, 3, 3 });
            p.AddCoeff(-1.6664946424001600e+03, new int[] { 1, 3, 5 });
            p.AddCoeff(3.5710599480003300e+02, new int[] { 1, 5, 1 });
            p.AddCoeff(-1.6664946424001600e+03, new int[] { 1, 5, 3 });
            p.AddCoeff(1.4998451781601400e+03, new int[] { 1, 5, 5 });
            p.AddCoeff(-3.9678443866670400e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(1.8516607137779500e+03, new int[] { 3, 1, 3 });
            p.AddCoeff(-1.6664946424001600e+03, new int[] { 3, 1, 5 });
            p.AddCoeff(1.8516607137779500e+03, new int[] { 3, 3, 1 });
            p.AddCoeff(-8.6410833309637700e+03, new int[] { 3, 3, 3 });
            p.AddCoeff(7.7769749978673900e+03, new int[] { 3, 3, 5 });
            p.AddCoeff(-1.6664946424001600e+03, new int[] { 3, 5, 1 });
            p.AddCoeff(7.7769749978673900e+03, new int[] { 3, 5, 3 });
            p.AddCoeff(-6.9992774980806500e+03, new int[] { 3, 5, 5 });
            p.AddCoeff(3.5710599480003300e+02, new int[] { 5, 1, 1 });
            p.AddCoeff(-1.6664946424001600e+03, new int[] { 5, 1, 3 });
            p.AddCoeff(1.4998451781601400e+03, new int[] { 5, 1, 5 });
            p.AddCoeff(-1.6664946424001600e+03, new int[] { 5, 3, 1 });
            p.AddCoeff(7.7769749978673900e+03, new int[] { 5, 3, 3 });
            p.AddCoeff(-6.9992774980806500e+03, new int[] { 5, 3, 5 });
            p.AddCoeff(1.4998451781601400e+03, new int[] { 5, 5, 1 });
            p.AddCoeff(-6.9992774980806500e+03, new int[] { 5, 5, 3 });
            p.AddCoeff(6.2993497482725900e+03, new int[] { 5, 5, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("21e651f3-e482-4d86-b6e8-de92d83b22c0"));
            OrthonormalPolynomials[756] = p;
            p.AddCoeff(-2.7869350108811700e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(2.7869350108811700e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-3.2514241793613600e+01, new int[] { 1, 0, 4 });
            p.AddCoeff(5.8525635228504500e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-5.8525635228504500e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(6.8279907766588600e+02, new int[] { 1, 2, 4 });
            p.AddCoeff(-1.7557690568551400e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(1.7557690568551400e+03, new int[] { 1, 4, 2 });
            p.AddCoeff(-2.0483972329976600e+03, new int[] { 1, 4, 4 });
            p.AddCoeff(1.2875639750271000e+02, new int[] { 1, 6, 0 });
            p.AddCoeff(-1.2875639750271000e+03, new int[] { 1, 6, 2 });
            p.AddCoeff(1.5021579708649500e+03, new int[] { 1, 6, 4 });
            p.AddCoeff(1.3005696717445500e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-1.3005696717445500e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(1.5173312837019700e+02, new int[] { 3, 0, 4 });
            p.AddCoeff(-2.7311963106635400e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(2.7311963106635400e+03, new int[] { 3, 2, 2 });
            p.AddCoeff(-3.1863956957741400e+03, new int[] { 3, 2, 4 });
            p.AddCoeff(8.1935889319906300e+02, new int[] { 3, 4, 0 });
            p.AddCoeff(-8.1935889319906300e+03, new int[] { 3, 4, 2 });
            p.AddCoeff(9.5591870873224100e+03, new int[] { 3, 4, 4 });
            p.AddCoeff(-6.0086318834598000e+02, new int[] { 3, 6, 0 });
            p.AddCoeff(6.0086318834598000e+03, new int[] { 3, 6, 2 });
            p.AddCoeff(-7.0100705307031000e+03, new int[] { 3, 6, 4 });
            p.AddCoeff(-1.1705127045700900e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(1.1705127045700900e+02, new int[] { 5, 0, 2 });
            p.AddCoeff(-1.3655981553317700e+02, new int[] { 5, 0, 4 });
            p.AddCoeff(2.4580766795971900e+02, new int[] { 5, 2, 0 });
            p.AddCoeff(-2.4580766795971900e+03, new int[] { 5, 2, 2 });
            p.AddCoeff(2.8677561261967200e+03, new int[] { 5, 2, 4 });
            p.AddCoeff(-7.3742300387915700e+02, new int[] { 5, 4, 0 });
            p.AddCoeff(7.3742300387915700e+03, new int[] { 5, 4, 2 });
            p.AddCoeff(-8.6032683785901700e+03, new int[] { 5, 4, 4 });
            p.AddCoeff(5.4077686951138200e+02, new int[] { 5, 6, 0 });
            p.AddCoeff(-5.4077686951138200e+03, new int[] { 5, 6, 2 });
            p.AddCoeff(6.3090634776327900e+03, new int[] { 5, 6, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("eef7b02e-8188-40f8-aa19-6284f6e59ae3"));
            OrthonormalPolynomials[757] = p;
            p.AddCoeff(7.3924192867575100e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-1.2320698811262500e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(-6.6531773580817600e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(1.1088628930136300e+03, new int[] { 1, 3, 3 });
            p.AddCoeff(1.4636990187779900e+03, new int[] { 1, 5, 1 });
            p.AddCoeff(-2.4394983646299800e+03, new int[] { 1, 5, 3 });
            p.AddCoeff(-9.0609939257684900e+02, new int[] { 1, 7, 1 });
            p.AddCoeff(1.5101656542947500e+03, new int[] { 1, 7, 3 });
            p.AddCoeff(-3.4497956671535100e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(5.7496594452558400e+02, new int[] { 3, 1, 3 });
            p.AddCoeff(3.1048161004381600e+03, new int[] { 3, 3, 1 });
            p.AddCoeff(-5.1746935007302600e+03, new int[] { 3, 3, 3 });
            p.AddCoeff(-6.8305954209639400e+03, new int[] { 3, 5, 1 });
            p.AddCoeff(1.1384325701606600e+04, new int[] { 3, 5, 3 });
            p.AddCoeff(4.2284638320253000e+03, new int[] { 3, 7, 1 });
            p.AddCoeff(-7.0474397200421600e+03, new int[] { 3, 7, 3 });
            p.AddCoeff(3.1048161004381600e+02, new int[] { 5, 1, 1 });
            p.AddCoeff(-5.1746935007302600e+02, new int[] { 5, 1, 3 });
            p.AddCoeff(-2.7943344903943400e+03, new int[] { 5, 3, 1 });
            p.AddCoeff(4.6572241506572300e+03, new int[] { 5, 3, 3 });
            p.AddCoeff(6.1475358788675500e+03, new int[] { 5, 5, 1 });
            p.AddCoeff(-1.0245893131445900e+04, new int[] { 5, 5, 3 });
            p.AddCoeff(-3.8056174488227700e+03, new int[] { 5, 7, 1 });
            p.AddCoeff(6.3426957480379500e+03, new int[] { 5, 7, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a1d41b35-a5de-4965-ab11-8b2a3dfe5960"));
            OrthonormalPolynomials[758] = p;
            p.AddCoeff(-2.7713422517043000e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(8.3140267551129000e+00, new int[] { 1, 0, 2 });
            p.AddCoeff(9.9768321061354800e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-2.9930496318406400e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(-5.4872576583745200e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(1.6461772975123500e+03, new int[] { 1, 4, 2 });
            p.AddCoeff(9.5112466078491600e+02, new int[] { 1, 6, 0 });
            p.AddCoeff(-2.8533739823547500e+03, new int[] { 1, 6, 2 });
            p.AddCoeff(-5.0953106827763400e+02, new int[] { 1, 8, 0 });
            p.AddCoeff(1.5285932048329000e+03, new int[] { 1, 8, 2 });
            p.AddCoeff(1.2932930507953400e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-3.8798791523860200e+01, new int[] { 3, 0, 2 });
            p.AddCoeff(-4.6558549828632300e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(1.3967564948589700e+03, new int[] { 3, 2, 2 });
            p.AddCoeff(2.5607202405747700e+03, new int[] { 3, 4, 0 });
            p.AddCoeff(-7.6821607217243200e+03, new int[] { 3, 4, 2 });
            p.AddCoeff(-4.4385817503296100e+03, new int[] { 3, 6, 0 });
            p.AddCoeff(1.3315745250988800e+04, new int[] { 3, 6, 2 });
            p.AddCoeff(2.3778116519622900e+03, new int[] { 3, 8, 0 });
            p.AddCoeff(-7.1334349558868700e+03, new int[] { 3, 8, 2 });
            p.AddCoeff(-1.1639637457158100e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(3.4918912371474200e+01, new int[] { 5, 0, 2 });
            p.AddCoeff(4.1902694845769000e+02, new int[] { 5, 2, 0 });
            p.AddCoeff(-1.2570808453730700e+03, new int[] { 5, 2, 2 });
            p.AddCoeff(-2.3046482165173000e+03, new int[] { 5, 4, 0 });
            p.AddCoeff(6.9139446495518900e+03, new int[] { 5, 4, 2 });
            p.AddCoeff(3.9947235752966500e+03, new int[] { 5, 6, 0 });
            p.AddCoeff(-1.1984170725889900e+04, new int[] { 5, 6, 2 });
            p.AddCoeff(-2.1400304867660600e+03, new int[] { 5, 8, 0 });
            p.AddCoeff(6.4200914602982000e+03, new int[] { 5, 8, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b277f64c-db05-4e31-9077-d719a5714fd9"));
            OrthonormalPolynomials[759] = p;
            p.AddCoeff(4.0849865705801600e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-5.9913136368509000e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(2.3366123183718500e+03, new int[] { 1, 5, 1 });
            p.AddCoeff(-3.3380175976740700e+03, new int[] { 1, 7, 1 });
            p.AddCoeff(1.5762860877905300e+03, new int[] { 1, 9, 1 });
            p.AddCoeff(-1.9063270662707400e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(2.7959463638637500e+03, new int[] { 3, 3, 1 });
            p.AddCoeff(-1.0904190819068600e+04, new int[] { 3, 5, 1 });
            p.AddCoeff(1.5577415455812300e+04, new int[] { 3, 7, 1 });
            p.AddCoeff(-7.3560017430224900e+03, new int[] { 3, 9, 1 });
            p.AddCoeff(1.7156943596436700e+02, new int[] { 5, 1, 1 });
            p.AddCoeff(-2.5163517274773800e+03, new int[] { 5, 3, 1 });
            p.AddCoeff(9.8137717371617800e+03, new int[] { 5, 5, 1 });
            p.AddCoeff(-1.4019673910231100e+04, new int[] { 5, 7, 1 });
            p.AddCoeff(6.6204015687202500e+03, new int[] { 5, 9, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("0fbfca94-d3cb-4a2b-80b0-5ca89808239e"));
            OrthonormalPolynomials[760] = p;
            p.AddCoeff(-2.4794928065055500e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(1.3637210435780500e+02, new int[] { 1, 2, 0 });
            p.AddCoeff(-1.1818915711009800e+03, new int[] { 1, 4, 0 });
            p.AddCoeff(3.5456747133029300e+03, new int[] { 1, 6, 0 });
            p.AddCoeff(-4.3054621518678400e+03, new int[] { 1, 8, 0 });
            p.AddCoeff(1.8178617974553100e+03, new int[] { 1, 10, 0 });
            p.AddCoeff(1.1570966430359200e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-6.3640315366975700e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(5.5154939984712300e+03, new int[] { 3, 4, 0 });
            p.AddCoeff(-1.6546481995413700e+04, new int[] { 3, 6, 0 });
            p.AddCoeff(2.0092156708716600e+04, new int[] { 3, 8, 0 });
            p.AddCoeff(-8.4833550547914600e+03, new int[] { 3, 10, 0 });
            p.AddCoeff(-1.0413869787323300e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(5.7276283830278100e+02, new int[] { 5, 2, 0 });
            p.AddCoeff(-4.9639445986241100e+03, new int[] { 5, 4, 0 });
            p.AddCoeff(1.4891833795872300e+04, new int[] { 5, 6, 0 });
            p.AddCoeff(-1.8082941037844900e+04, new int[] { 5, 8, 0 });
            p.AddCoeff(7.6350195493123200e+03, new int[] { 5, 10, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("ead54c0e-e0e4-4950-83a1-7b13bce9c0e0"));
            OrthonormalPolynomials[761] = p;
            p.AddCoeff(-4.2732085527534600e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(6.2673725440384000e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-2.4442752921749800e+02, new int[] { 0, 0, 5 });
            p.AddCoeff(3.4918218459642500e+02, new int[] { 0, 0, 7 });
            p.AddCoeff(-1.6489158717053400e+02, new int[] { 0, 0, 9 });
            p.AddCoeff(8.9737379607822600e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-1.3161482342480600e+03, new int[] { 2, 0, 3 });
            p.AddCoeff(5.1329781135674500e+03, new int[] { 2, 0, 5 });
            p.AddCoeff(-7.3328258765249300e+03, new int[] { 2, 0, 7 });
            p.AddCoeff(3.4627233305812300e+03, new int[] { 2, 0, 9 });
            p.AddCoeff(-2.6921213882346800e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(3.9484447027441900e+03, new int[] { 4, 0, 3 });
            p.AddCoeff(-1.5398934340702400e+04, new int[] { 4, 0, 5 });
            p.AddCoeff(2.1998477629574800e+04, new int[] { 4, 0, 7 });
            p.AddCoeff(-1.0388169991743600e+04, new int[] { 4, 0, 9 });
            p.AddCoeff(1.9742223513721000e+02, new int[] { 6, 0, 1 });
            p.AddCoeff(-2.8955261153457400e+03, new int[] { 6, 0, 3 });
            p.AddCoeff(1.1292551849848400e+04, new int[] { 6, 0, 5 });
            p.AddCoeff(-1.6132216928354800e+04, new int[] { 6, 0, 7 });
            p.AddCoeff(7.6179913272786700e+03, new int[] { 6, 0, 9 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("42c96878-d9ba-46b3-a99e-fc6ab716d8cf"));
            OrthonormalPolynomials[762] = p;
            p.AddCoeff(-7.7789300654438300e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(2.8004148235597800e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-1.5402281529578800e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(2.6697287984603200e+02, new int[] { 0, 1, 6 });
            p.AddCoeff(-1.4302118563180300e+02, new int[] { 0, 1, 8 });
            p.AddCoeff(1.6335753137432000e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-5.8808711294755300e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(3.2344791212115400e+03, new int[] { 2, 1, 4 });
            p.AddCoeff(-5.6064304767666800e+03, new int[] { 2, 1, 6 });
            p.AddCoeff(3.0034448982678600e+03, new int[] { 2, 1, 8 });
            p.AddCoeff(-4.9007259412296100e+01, new int[] { 4, 1, 0 });
            p.AddCoeff(1.7642613388426600e+03, new int[] { 4, 1, 2 });
            p.AddCoeff(-9.7034373636346300e+03, new int[] { 4, 1, 4 });
            p.AddCoeff(1.6819291430300000e+04, new int[] { 4, 1, 6 });
            p.AddCoeff(-9.0103346948035700e+03, new int[] { 4, 1, 8 });
            p.AddCoeff(3.5938656902350500e+01, new int[] { 6, 1, 0 });
            p.AddCoeff(-1.2937916484846200e+03, new int[] { 6, 1, 2 });
            p.AddCoeff(7.1158540666654000e+03, new int[] { 6, 1, 4 });
            p.AddCoeff(-1.2334147048886700e+04, new int[] { 6, 1, 6 });
            p.AddCoeff(6.6075787761892800e+03, new int[] { 6, 1, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("d67e5a8e-7e1c-4970-ba69-19c8a0a12f19"));
            OrthonormalPolynomials[763] = p;
            p.AddCoeff(-3.7733353310726900e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(3.3960017979654200e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-7.4712039555239300e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(4.6250310200862400e+01, new int[] { 0, 0, 7 });
            p.AddCoeff(1.1320005993218100e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-1.0188005393896300e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(2.2413611866571800e+02, new int[] { 0, 2, 5 });
            p.AddCoeff(-1.3875093060258700e+02, new int[] { 0, 2, 7 });
            p.AddCoeff(7.9240041952526500e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-7.1316037757273800e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(1.5689528306600200e+03, new int[] { 2, 0, 5 });
            p.AddCoeff(-9.7125651421811000e+02, new int[] { 2, 0, 7 });
            p.AddCoeff(-2.3772012585757900e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(2.1394811327182200e+03, new int[] { 2, 2, 3 });
            p.AddCoeff(-4.7068584919800700e+03, new int[] { 2, 2, 5 });
            p.AddCoeff(2.9137695426543300e+03, new int[] { 2, 2, 7 });
            p.AddCoeff(-2.3772012585757900e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(2.1394811327182200e+03, new int[] { 4, 0, 3 });
            p.AddCoeff(-4.7068584919800700e+03, new int[] { 4, 0, 5 });
            p.AddCoeff(2.9137695426543300e+03, new int[] { 4, 0, 7 });
            p.AddCoeff(7.1316037757273800e+02, new int[] { 4, 2, 1 });
            p.AddCoeff(-6.4184433981546500e+03, new int[] { 4, 2, 3 });
            p.AddCoeff(1.4120575475940200e+04, new int[] { 4, 2, 5 });
            p.AddCoeff(-8.7413086279629900e+03, new int[] { 4, 2, 7 });
            p.AddCoeff(1.7432809229555800e+02, new int[] { 6, 0, 1 });
            p.AddCoeff(-1.5689528306600200e+03, new int[] { 6, 0, 3 });
            p.AddCoeff(3.4516962274520500e+03, new int[] { 6, 0, 5 });
            p.AddCoeff(-2.1367643312798400e+03, new int[] { 6, 0, 7 });
            p.AddCoeff(-5.2298427688667500e+02, new int[] { 6, 2, 1 });
            p.AddCoeff(4.7068584919800700e+03, new int[] { 6, 2, 3 });
            p.AddCoeff(-1.0355088682356200e+04, new int[] { 6, 2, 5 });
            p.AddCoeff(6.4102929938395300e+03, new int[] { 6, 2, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3788c752-ce74-4fd5-be8a-7d2b795bd1d3"));
            OrthonormalPolynomials[764] = p;
            p.AddCoeff(-1.7813066172385700e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(3.7407438962010000e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-1.1222231688603000e+02, new int[] { 0, 1, 4 });
            p.AddCoeff(8.2296365716421900e+01, new int[] { 0, 1, 6 });
            p.AddCoeff(2.9688443620642800e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(-6.2345731603350000e+01, new int[] { 0, 3, 2 });
            p.AddCoeff(1.8703719481005000e+02, new int[] { 0, 3, 4 });
            p.AddCoeff(-1.3716060952737000e+02, new int[] { 0, 3, 6 });
            p.AddCoeff(3.7407438962010000e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-7.8555621820220900e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(2.3566686546066300e+03, new int[] { 2, 1, 4 });
            p.AddCoeff(-1.7282236800448600e+03, new int[] { 2, 1, 6 });
            p.AddCoeff(-6.2345731603350000e+01, new int[] { 2, 3, 0 });
            p.AddCoeff(1.3092603636703500e+03, new int[] { 2, 3, 2 });
            p.AddCoeff(-3.9277810910110500e+03, new int[] { 2, 3, 4 });
            p.AddCoeff(2.8803728000747700e+03, new int[] { 2, 3, 6 });
            p.AddCoeff(-1.1222231688603000e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(2.3566686546066300e+03, new int[] { 4, 1, 2 });
            p.AddCoeff(-7.0700059638198800e+03, new int[] { 4, 1, 4 });
            p.AddCoeff(5.1846710401345800e+03, new int[] { 4, 1, 6 });
            p.AddCoeff(1.8703719481005000e+02, new int[] { 4, 3, 0 });
            p.AddCoeff(-3.9277810910110500e+03, new int[] { 4, 3, 2 });
            p.AddCoeff(1.1783343273033100e+04, new int[] { 4, 3, 4 });
            p.AddCoeff(-8.6411184002243000e+03, new int[] { 4, 3, 6 });
            p.AddCoeff(8.2296365716421900e+01, new int[] { 6, 1, 0 });
            p.AddCoeff(-1.7282236800448600e+03, new int[] { 6, 1, 2 });
            p.AddCoeff(5.1846710401345800e+03, new int[] { 6, 1, 4 });
            p.AddCoeff(-3.8020920960986900e+03, new int[] { 6, 1, 6 });
            p.AddCoeff(-1.3716060952737000e+02, new int[] { 6, 3, 0 });
            p.AddCoeff(2.8803728000747700e+03, new int[] { 6, 3, 2 });
            p.AddCoeff(-8.6411184002243000e+03, new int[] { 6, 3, 4 });
            p.AddCoeff(6.3368201601644900e+03, new int[] { 6, 3, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("2d51e622-d35a-4513-aacb-366dabbf1c01"));
            OrthonormalPolynomials[765] = p;
            p.AddCoeff(-2.7869350108811700e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(1.3005696717445500e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-1.1705127045700900e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(2.7869350108811700e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-1.3005696717445500e+02, new int[] { 0, 2, 3 });
            p.AddCoeff(1.1705127045700900e+02, new int[] { 0, 2, 5 });
            p.AddCoeff(-3.2514241793613600e+01, new int[] { 0, 4, 1 });
            p.AddCoeff(1.5173312837019700e+02, new int[] { 0, 4, 3 });
            p.AddCoeff(-1.3655981553317700e+02, new int[] { 0, 4, 5 });
            p.AddCoeff(5.8525635228504500e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-2.7311963106635400e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(2.4580766795971900e+02, new int[] { 2, 0, 5 });
            p.AddCoeff(-5.8525635228504500e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(2.7311963106635400e+03, new int[] { 2, 2, 3 });
            p.AddCoeff(-2.4580766795971900e+03, new int[] { 2, 2, 5 });
            p.AddCoeff(6.8279907766588600e+02, new int[] { 2, 4, 1 });
            p.AddCoeff(-3.1863956957741400e+03, new int[] { 2, 4, 3 });
            p.AddCoeff(2.8677561261967200e+03, new int[] { 2, 4, 5 });
            p.AddCoeff(-1.7557690568551400e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(8.1935889319906300e+02, new int[] { 4, 0, 3 });
            p.AddCoeff(-7.3742300387915700e+02, new int[] { 4, 0, 5 });
            p.AddCoeff(1.7557690568551400e+03, new int[] { 4, 2, 1 });
            p.AddCoeff(-8.1935889319906300e+03, new int[] { 4, 2, 3 });
            p.AddCoeff(7.3742300387915700e+03, new int[] { 4, 2, 5 });
            p.AddCoeff(-2.0483972329976600e+03, new int[] { 4, 4, 1 });
            p.AddCoeff(9.5591870873224100e+03, new int[] { 4, 4, 3 });
            p.AddCoeff(-8.6032683785901700e+03, new int[] { 4, 4, 5 });
            p.AddCoeff(1.2875639750271000e+02, new int[] { 6, 0, 1 });
            p.AddCoeff(-6.0086318834598000e+02, new int[] { 6, 0, 3 });
            p.AddCoeff(5.4077686951138200e+02, new int[] { 6, 0, 5 });
            p.AddCoeff(-1.2875639750271000e+03, new int[] { 6, 2, 1 });
            p.AddCoeff(6.0086318834598000e+03, new int[] { 6, 2, 3 });
            p.AddCoeff(-5.4077686951138200e+03, new int[] { 6, 2, 5 });
            p.AddCoeff(1.5021579708649500e+03, new int[] { 6, 4, 1 });
            p.AddCoeff(-7.0100705307031000e+03, new int[] { 6, 4, 3 });
            p.AddCoeff(6.3090634776327900e+03, new int[] { 6, 4, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c0b55117-30f5-4c18-89a4-25c6725a494e"));
            OrthonormalPolynomials[766] = p;
            p.AddCoeff(-2.7869350108811700e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(2.7869350108811700e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-3.2514241793613600e+01, new int[] { 0, 1, 4 });
            p.AddCoeff(1.3005696717445500e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-1.3005696717445500e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(1.5173312837019700e+02, new int[] { 0, 3, 4 });
            p.AddCoeff(-1.1705127045700900e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(1.1705127045700900e+02, new int[] { 0, 5, 2 });
            p.AddCoeff(-1.3655981553317700e+02, new int[] { 0, 5, 4 });
            p.AddCoeff(5.8525635228504500e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-5.8525635228504500e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(6.8279907766588600e+02, new int[] { 2, 1, 4 });
            p.AddCoeff(-2.7311963106635400e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(2.7311963106635400e+03, new int[] { 2, 3, 2 });
            p.AddCoeff(-3.1863956957741400e+03, new int[] { 2, 3, 4 });
            p.AddCoeff(2.4580766795971900e+02, new int[] { 2, 5, 0 });
            p.AddCoeff(-2.4580766795971900e+03, new int[] { 2, 5, 2 });
            p.AddCoeff(2.8677561261967200e+03, new int[] { 2, 5, 4 });
            p.AddCoeff(-1.7557690568551400e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(1.7557690568551400e+03, new int[] { 4, 1, 2 });
            p.AddCoeff(-2.0483972329976600e+03, new int[] { 4, 1, 4 });
            p.AddCoeff(8.1935889319906300e+02, new int[] { 4, 3, 0 });
            p.AddCoeff(-8.1935889319906300e+03, new int[] { 4, 3, 2 });
            p.AddCoeff(9.5591870873224100e+03, new int[] { 4, 3, 4 });
            p.AddCoeff(-7.3742300387915700e+02, new int[] { 4, 5, 0 });
            p.AddCoeff(7.3742300387915700e+03, new int[] { 4, 5, 2 });
            p.AddCoeff(-8.6032683785901700e+03, new int[] { 4, 5, 4 });
            p.AddCoeff(1.2875639750271000e+02, new int[] { 6, 1, 0 });
            p.AddCoeff(-1.2875639750271000e+03, new int[] { 6, 1, 2 });
            p.AddCoeff(1.5021579708649500e+03, new int[] { 6, 1, 4 });
            p.AddCoeff(-6.0086318834598000e+02, new int[] { 6, 3, 0 });
            p.AddCoeff(6.0086318834598000e+03, new int[] { 6, 3, 2 });
            p.AddCoeff(-7.0100705307031000e+03, new int[] { 6, 3, 4 });
            p.AddCoeff(5.4077686951138200e+02, new int[] { 6, 5, 0 });
            p.AddCoeff(-5.4077686951138200e+03, new int[] { 6, 5, 2 });
            p.AddCoeff(6.3090634776327900e+03, new int[] { 6, 5, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("fcbf705e-421a-4527-8569-a2eed40ec2e0"));
            OrthonormalPolynomials[767] = p;
            p.AddCoeff(-1.7813066172385700e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(2.9688443620642800e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(3.7407438962010000e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-6.2345731603350000e+01, new int[] { 0, 2, 3 });
            p.AddCoeff(-1.1222231688603000e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(1.8703719481005000e+02, new int[] { 0, 4, 3 });
            p.AddCoeff(8.2296365716421900e+01, new int[] { 0, 6, 1 });
            p.AddCoeff(-1.3716060952737000e+02, new int[] { 0, 6, 3 });
            p.AddCoeff(3.7407438962010000e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-6.2345731603350000e+01, new int[] { 2, 0, 3 });
            p.AddCoeff(-7.8555621820220900e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(1.3092603636703500e+03, new int[] { 2, 2, 3 });
            p.AddCoeff(2.3566686546066300e+03, new int[] { 2, 4, 1 });
            p.AddCoeff(-3.9277810910110500e+03, new int[] { 2, 4, 3 });
            p.AddCoeff(-1.7282236800448600e+03, new int[] { 2, 6, 1 });
            p.AddCoeff(2.8803728000747700e+03, new int[] { 2, 6, 3 });
            p.AddCoeff(-1.1222231688603000e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(1.8703719481005000e+02, new int[] { 4, 0, 3 });
            p.AddCoeff(2.3566686546066300e+03, new int[] { 4, 2, 1 });
            p.AddCoeff(-3.9277810910110500e+03, new int[] { 4, 2, 3 });
            p.AddCoeff(-7.0700059638198800e+03, new int[] { 4, 4, 1 });
            p.AddCoeff(1.1783343273033100e+04, new int[] { 4, 4, 3 });
            p.AddCoeff(5.1846710401345800e+03, new int[] { 4, 6, 1 });
            p.AddCoeff(-8.6411184002243000e+03, new int[] { 4, 6, 3 });
            p.AddCoeff(8.2296365716421900e+01, new int[] { 6, 0, 1 });
            p.AddCoeff(-1.3716060952737000e+02, new int[] { 6, 0, 3 });
            p.AddCoeff(-1.7282236800448600e+03, new int[] { 6, 2, 1 });
            p.AddCoeff(2.8803728000747700e+03, new int[] { 6, 2, 3 });
            p.AddCoeff(5.1846710401345800e+03, new int[] { 6, 4, 1 });
            p.AddCoeff(-8.6411184002243000e+03, new int[] { 6, 4, 3 });
            p.AddCoeff(-3.8020920960986900e+03, new int[] { 6, 6, 1 });
            p.AddCoeff(6.3368201601644900e+03, new int[] { 6, 6, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("6a41588f-77c1-427c-815b-21381421c6c0"));
            OrthonormalPolynomials[768] = p;
            p.AddCoeff(-3.7733353310726900e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(1.1320005993218100e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(3.3960017979654200e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-1.0188005393896300e+02, new int[] { 0, 3, 2 });
            p.AddCoeff(-7.4712039555239300e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(2.2413611866571800e+02, new int[] { 0, 5, 2 });
            p.AddCoeff(4.6250310200862400e+01, new int[] { 0, 7, 0 });
            p.AddCoeff(-1.3875093060258700e+02, new int[] { 0, 7, 2 });
            p.AddCoeff(7.9240041952526500e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-2.3772012585757900e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(-7.1316037757273800e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(2.1394811327182200e+03, new int[] { 2, 3, 2 });
            p.AddCoeff(1.5689528306600200e+03, new int[] { 2, 5, 0 });
            p.AddCoeff(-4.7068584919800700e+03, new int[] { 2, 5, 2 });
            p.AddCoeff(-9.7125651421811000e+02, new int[] { 2, 7, 0 });
            p.AddCoeff(2.9137695426543300e+03, new int[] { 2, 7, 2 });
            p.AddCoeff(-2.3772012585757900e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(7.1316037757273800e+02, new int[] { 4, 1, 2 });
            p.AddCoeff(2.1394811327182200e+03, new int[] { 4, 3, 0 });
            p.AddCoeff(-6.4184433981546500e+03, new int[] { 4, 3, 2 });
            p.AddCoeff(-4.7068584919800700e+03, new int[] { 4, 5, 0 });
            p.AddCoeff(1.4120575475940200e+04, new int[] { 4, 5, 2 });
            p.AddCoeff(2.9137695426543300e+03, new int[] { 4, 7, 0 });
            p.AddCoeff(-8.7413086279629900e+03, new int[] { 4, 7, 2 });
            p.AddCoeff(1.7432809229555800e+02, new int[] { 6, 1, 0 });
            p.AddCoeff(-5.2298427688667500e+02, new int[] { 6, 1, 2 });
            p.AddCoeff(-1.5689528306600200e+03, new int[] { 6, 3, 0 });
            p.AddCoeff(4.7068584919800700e+03, new int[] { 6, 3, 2 });
            p.AddCoeff(3.4516962274520500e+03, new int[] { 6, 5, 0 });
            p.AddCoeff(-1.0355088682356200e+04, new int[] { 6, 5, 2 });
            p.AddCoeff(-2.1367643312798400e+03, new int[] { 6, 7, 0 });
            p.AddCoeff(6.4102929938395300e+03, new int[] { 6, 7, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("0e3f2122-b991-4dbf-9aac-ca44557e1812"));
            OrthonormalPolynomials[769] = p;
            p.AddCoeff(-7.7789300654438300e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(2.8004148235597800e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-1.5402281529578800e+02, new int[] { 0, 4, 1 });
            p.AddCoeff(2.6697287984603200e+02, new int[] { 0, 6, 1 });
            p.AddCoeff(-1.4302118563180300e+02, new int[] { 0, 8, 1 });
            p.AddCoeff(1.6335753137432000e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-5.8808711294755300e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(3.2344791212115400e+03, new int[] { 2, 4, 1 });
            p.AddCoeff(-5.6064304767666800e+03, new int[] { 2, 6, 1 });
            p.AddCoeff(3.0034448982678600e+03, new int[] { 2, 8, 1 });
            p.AddCoeff(-4.9007259412296100e+01, new int[] { 4, 0, 1 });
            p.AddCoeff(1.7642613388426600e+03, new int[] { 4, 2, 1 });
            p.AddCoeff(-9.7034373636346300e+03, new int[] { 4, 4, 1 });
            p.AddCoeff(1.6819291430300000e+04, new int[] { 4, 6, 1 });
            p.AddCoeff(-9.0103346948035700e+03, new int[] { 4, 8, 1 });
            p.AddCoeff(3.5938656902350500e+01, new int[] { 6, 0, 1 });
            p.AddCoeff(-1.2937916484846200e+03, new int[] { 6, 2, 1 });
            p.AddCoeff(7.1158540666654000e+03, new int[] { 6, 4, 1 });
            p.AddCoeff(-1.2334147048886700e+04, new int[] { 6, 6, 1 });
            p.AddCoeff(6.6075787761892800e+03, new int[] { 6, 8, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("0eea5dc8-ff14-4476-b230-be496b5954a0"));
            OrthonormalPolynomials[770] = p;
            p.AddCoeff(-4.2732085527534600e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(6.2673725440384000e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-2.4442752921749800e+02, new int[] { 0, 5, 0 });
            p.AddCoeff(3.4918218459642500e+02, new int[] { 0, 7, 0 });
            p.AddCoeff(-1.6489158717053400e+02, new int[] { 0, 9, 0 });
            p.AddCoeff(8.9737379607822600e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-1.3161482342480600e+03, new int[] { 2, 3, 0 });
            p.AddCoeff(5.1329781135674500e+03, new int[] { 2, 5, 0 });
            p.AddCoeff(-7.3328258765249300e+03, new int[] { 2, 7, 0 });
            p.AddCoeff(3.4627233305812300e+03, new int[] { 2, 9, 0 });
            p.AddCoeff(-2.6921213882346800e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(3.9484447027441900e+03, new int[] { 4, 3, 0 });
            p.AddCoeff(-1.5398934340702400e+04, new int[] { 4, 5, 0 });
            p.AddCoeff(2.1998477629574800e+04, new int[] { 4, 7, 0 });
            p.AddCoeff(-1.0388169991743600e+04, new int[] { 4, 9, 0 });
            p.AddCoeff(1.9742223513721000e+02, new int[] { 6, 1, 0 });
            p.AddCoeff(-2.8955261153457400e+03, new int[] { 6, 3, 0 });
            p.AddCoeff(1.1292551849848400e+04, new int[] { 6, 5, 0 });
            p.AddCoeff(-1.6132216928354800e+04, new int[] { 6, 7, 0 });
            p.AddCoeff(7.6179913272786700e+03, new int[] { 6, 9, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("fc5f43fc-211c-4fb7-a4d9-e8b429534374"));
            OrthonormalPolynomials[771] = p;
            p.AddCoeff(-3.3770013411936900e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(1.2157204828297300e+02, new int[] { 1, 0, 2 });
            p.AddCoeff(-6.6864626555635000e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(1.1589868602976700e+03, new int[] { 1, 0, 6 });
            p.AddCoeff(-6.2088581801661100e+02, new int[] { 1, 0, 8 });
            p.AddCoeff(3.0393012070743200e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-1.0941484345467500e+03, new int[] { 3, 0, 2 });
            p.AddCoeff(6.0178163900071500e+03, new int[] { 3, 0, 4 });
            p.AddCoeff(-1.0430881742679100e+04, new int[] { 3, 0, 6 });
            p.AddCoeff(5.5879723621494900e+03, new int[] { 3, 0, 8 });
            p.AddCoeff(-6.6864626555635000e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(2.4071265560028600e+03, new int[] { 5, 0, 2 });
            p.AddCoeff(-1.3239196058015700e+04, new int[] { 5, 0, 4 });
            p.AddCoeff(2.2947939833893900e+04, new int[] { 5, 0, 6 });
            p.AddCoeff(-1.2293539196728900e+04, new int[] { 5, 0, 8 });
            p.AddCoeff(4.1392387867774100e+01, new int[] { 7, 0, 0 });
            p.AddCoeff(-1.4901259632398700e+03, new int[] { 7, 0, 2 });
            p.AddCoeff(8.1956927978192600e+03, new int[] { 7, 0, 4 });
            p.AddCoeff(-1.4205867516220100e+04, new int[] { 7, 0, 6 });
            p.AddCoeff(7.6102861694036100e+03, new int[] { 7, 0, 8 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("2bbc3947-7795-450f-942a-8ef72dd5bef4"));
            OrthonormalPolynomials[772] = p;
            p.AddCoeff(4.3954466819961800e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-3.9559020137965700e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(8.7029844303524500e+02, new int[] { 1, 1, 5 });
            p.AddCoeff(-5.3875617902181800e+02, new int[] { 1, 1, 7 });
            p.AddCoeff(-3.9559020137965700e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(3.5603118124169100e+03, new int[] { 3, 1, 3 });
            p.AddCoeff(-7.8326859873172000e+03, new int[] { 3, 1, 5 });
            p.AddCoeff(4.8488056111963600e+03, new int[] { 3, 1, 7 });
            p.AddCoeff(8.7029844303524500e+02, new int[] { 5, 1, 1 });
            p.AddCoeff(-7.8326859873172000e+03, new int[] { 5, 1, 3 });
            p.AddCoeff(1.7231909172097800e+04, new int[] { 5, 1, 5 });
            p.AddCoeff(-1.0667372344632000e+04, new int[] { 5, 1, 7 });
            p.AddCoeff(-5.3875617902181800e+02, new int[] { 7, 1, 1 });
            p.AddCoeff(4.8488056111963600e+03, new int[] { 7, 1, 3 });
            p.AddCoeff(-1.0667372344632000e+04, new int[] { 7, 1, 5 });
            p.AddCoeff(6.6036114514388600e+03, new int[] { 7, 1, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a1bcc9af-1fbc-4ba3-9fac-56071ae683e6"));
            OrthonormalPolynomials[773] = p;
            p.AddCoeff(-3.7733353310726900e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(7.9240041952526500e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-2.3772012585757900e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(1.7432809229555800e+02, new int[] { 1, 0, 6 });
            p.AddCoeff(1.1320005993218100e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-2.3772012585757900e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(7.1316037757273800e+02, new int[] { 1, 2, 4 });
            p.AddCoeff(-5.2298427688667500e+02, new int[] { 1, 2, 6 });
            p.AddCoeff(3.3960017979654200e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-7.1316037757273800e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(2.1394811327182200e+03, new int[] { 3, 0, 4 });
            p.AddCoeff(-1.5689528306600200e+03, new int[] { 3, 0, 6 });
            p.AddCoeff(-1.0188005393896300e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(2.1394811327182200e+03, new int[] { 3, 2, 2 });
            p.AddCoeff(-6.4184433981546500e+03, new int[] { 3, 2, 4 });
            p.AddCoeff(4.7068584919800700e+03, new int[] { 3, 2, 6 });
            p.AddCoeff(-7.4712039555239300e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(1.5689528306600200e+03, new int[] { 5, 0, 2 });
            p.AddCoeff(-4.7068584919800700e+03, new int[] { 5, 0, 4 });
            p.AddCoeff(3.4516962274520500e+03, new int[] { 5, 0, 6 });
            p.AddCoeff(2.2413611866571800e+02, new int[] { 5, 2, 0 });
            p.AddCoeff(-4.7068584919800700e+03, new int[] { 5, 2, 2 });
            p.AddCoeff(1.4120575475940200e+04, new int[] { 5, 2, 4 });
            p.AddCoeff(-1.0355088682356200e+04, new int[] { 5, 2, 6 });
            p.AddCoeff(4.6250310200862400e+01, new int[] { 7, 0, 0 });
            p.AddCoeff(-9.7125651421811000e+02, new int[] { 7, 0, 2 });
            p.AddCoeff(2.9137695426543300e+03, new int[] { 7, 0, 4 });
            p.AddCoeff(-2.1367643312798400e+03, new int[] { 7, 0, 6 });
            p.AddCoeff(-1.3875093060258700e+02, new int[] { 7, 2, 0 });
            p.AddCoeff(2.9137695426543300e+03, new int[] { 7, 2, 2 });
            p.AddCoeff(-8.7413086279629900e+03, new int[] { 7, 2, 4 });
            p.AddCoeff(6.4102929938395300e+03, new int[] { 7, 2, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("97cffba2-e059-4604-a543-bc54d25cd273"));
            OrthonormalPolynomials[774] = p;
            p.AddCoeff(7.3924192867575100e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-3.4497956671535100e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(3.1048161004381600e+02, new int[] { 1, 1, 5 });
            p.AddCoeff(-1.2320698811262500e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(5.7496594452558400e+02, new int[] { 1, 3, 3 });
            p.AddCoeff(-5.1746935007302600e+02, new int[] { 1, 3, 5 });
            p.AddCoeff(-6.6531773580817600e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(3.1048161004381600e+03, new int[] { 3, 1, 3 });
            p.AddCoeff(-2.7943344903943400e+03, new int[] { 3, 1, 5 });
            p.AddCoeff(1.1088628930136300e+03, new int[] { 3, 3, 1 });
            p.AddCoeff(-5.1746935007302600e+03, new int[] { 3, 3, 3 });
            p.AddCoeff(4.6572241506572300e+03, new int[] { 3, 3, 5 });
            p.AddCoeff(1.4636990187779900e+03, new int[] { 5, 1, 1 });
            p.AddCoeff(-6.8305954209639400e+03, new int[] { 5, 1, 3 });
            p.AddCoeff(6.1475358788675500e+03, new int[] { 5, 1, 5 });
            p.AddCoeff(-2.4394983646299800e+03, new int[] { 5, 3, 1 });
            p.AddCoeff(1.1384325701606600e+04, new int[] { 5, 3, 3 });
            p.AddCoeff(-1.0245893131445900e+04, new int[] { 5, 3, 5 });
            p.AddCoeff(-9.0609939257684900e+02, new int[] { 7, 1, 1 });
            p.AddCoeff(4.2284638320253000e+03, new int[] { 7, 1, 3 });
            p.AddCoeff(-3.8056174488227700e+03, new int[] { 7, 1, 5 });
            p.AddCoeff(1.5101656542947500e+03, new int[] { 7, 3, 1 });
            p.AddCoeff(-7.0474397200421600e+03, new int[] { 7, 3, 3 });
            p.AddCoeff(6.3426957480379500e+03, new int[] { 7, 3, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("37e547c5-d680-452f-98aa-07b039981680"));
            OrthonormalPolynomials[775] = p;
            p.AddCoeff(-3.7909996350760400e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(3.7909996350760400e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-4.4228329075887100e+01, new int[] { 1, 0, 4 });
            p.AddCoeff(3.7909996350760400e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-3.7909996350760400e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(4.4228329075887100e+02, new int[] { 1, 2, 4 });
            p.AddCoeff(-4.4228329075887100e+01, new int[] { 1, 4, 0 });
            p.AddCoeff(4.4228329075887100e+02, new int[] { 1, 4, 2 });
            p.AddCoeff(-5.1599717255201600e+02, new int[] { 1, 4, 4 });
            p.AddCoeff(3.4118996715684400e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-3.4118996715684400e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(3.9805496168298400e+02, new int[] { 3, 0, 4 });
            p.AddCoeff(-3.4118996715684400e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(3.4118996715684400e+03, new int[] { 3, 2, 2 });
            p.AddCoeff(-3.9805496168298400e+03, new int[] { 3, 2, 4 });
            p.AddCoeff(3.9805496168298400e+02, new int[] { 3, 4, 0 });
            p.AddCoeff(-3.9805496168298400e+03, new int[] { 3, 4, 2 });
            p.AddCoeff(4.6439745529681500e+03, new int[] { 3, 4, 4 });
            p.AddCoeff(-7.5061792774505600e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(7.5061792774505600e+02, new int[] { 5, 0, 2 });
            p.AddCoeff(-8.7572091570256500e+02, new int[] { 5, 0, 4 });
            p.AddCoeff(7.5061792774505600e+02, new int[] { 5, 2, 0 });
            p.AddCoeff(-7.5061792774505600e+03, new int[] { 5, 2, 2 });
            p.AddCoeff(8.7572091570256500e+03, new int[] { 5, 2, 4 });
            p.AddCoeff(-8.7572091570256500e+02, new int[] { 5, 4, 0 });
            p.AddCoeff(8.7572091570256500e+03, new int[] { 5, 4, 2 });
            p.AddCoeff(-1.0216744016529900e+04, new int[] { 5, 4, 4 });
            p.AddCoeff(4.6466824098503500e+01, new int[] { 7, 0, 0 });
            p.AddCoeff(-4.6466824098503500e+02, new int[] { 7, 0, 2 });
            p.AddCoeff(5.4211294781587400e+02, new int[] { 7, 0, 4 });
            p.AddCoeff(-4.6466824098503500e+02, new int[] { 7, 2, 0 });
            p.AddCoeff(4.6466824098503500e+03, new int[] { 7, 2, 2 });
            p.AddCoeff(-5.4211294781587400e+03, new int[] { 7, 2, 4 });
            p.AddCoeff(5.4211294781587400e+02, new int[] { 7, 4, 0 });
            p.AddCoeff(-5.4211294781587400e+03, new int[] { 7, 4, 2 });
            p.AddCoeff(6.3246510578518300e+03, new int[] { 7, 4, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("83c6fc1e-ff9b-4e97-8617-916e825f6899"));
            OrthonormalPolynomials[776] = p;
            p.AddCoeff(7.3924192867575100e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-1.2320698811262500e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(-3.4497956671535100e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(5.7496594452558400e+02, new int[] { 1, 3, 3 });
            p.AddCoeff(3.1048161004381600e+02, new int[] { 1, 5, 1 });
            p.AddCoeff(-5.1746935007302600e+02, new int[] { 1, 5, 3 });
            p.AddCoeff(-6.6531773580817600e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(1.1088628930136300e+03, new int[] { 3, 1, 3 });
            p.AddCoeff(3.1048161004381600e+03, new int[] { 3, 3, 1 });
            p.AddCoeff(-5.1746935007302600e+03, new int[] { 3, 3, 3 });
            p.AddCoeff(-2.7943344903943400e+03, new int[] { 3, 5, 1 });
            p.AddCoeff(4.6572241506572300e+03, new int[] { 3, 5, 3 });
            p.AddCoeff(1.4636990187779900e+03, new int[] { 5, 1, 1 });
            p.AddCoeff(-2.4394983646299800e+03, new int[] { 5, 1, 3 });
            p.AddCoeff(-6.8305954209639400e+03, new int[] { 5, 3, 1 });
            p.AddCoeff(1.1384325701606600e+04, new int[] { 5, 3, 3 });
            p.AddCoeff(6.1475358788675500e+03, new int[] { 5, 5, 1 });
            p.AddCoeff(-1.0245893131445900e+04, new int[] { 5, 5, 3 });
            p.AddCoeff(-9.0609939257684900e+02, new int[] { 7, 1, 1 });
            p.AddCoeff(1.5101656542947500e+03, new int[] { 7, 1, 3 });
            p.AddCoeff(4.2284638320253000e+03, new int[] { 7, 3, 1 });
            p.AddCoeff(-7.0474397200421600e+03, new int[] { 7, 3, 3 });
            p.AddCoeff(-3.8056174488227700e+03, new int[] { 7, 5, 1 });
            p.AddCoeff(6.3426957480379500e+03, new int[] { 7, 5, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("23678a57-1461-4297-b336-316f6490841b"));
            OrthonormalPolynomials[777] = p;
            p.AddCoeff(-3.7733353310726900e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(1.1320005993218100e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(7.9240041952526500e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-2.3772012585757900e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(-2.3772012585757900e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(7.1316037757273800e+02, new int[] { 1, 4, 2 });
            p.AddCoeff(1.7432809229555800e+02, new int[] { 1, 6, 0 });
            p.AddCoeff(-5.2298427688667500e+02, new int[] { 1, 6, 2 });
            p.AddCoeff(3.3960017979654200e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-1.0188005393896300e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(-7.1316037757273800e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(2.1394811327182200e+03, new int[] { 3, 2, 2 });
            p.AddCoeff(2.1394811327182200e+03, new int[] { 3, 4, 0 });
            p.AddCoeff(-6.4184433981546500e+03, new int[] { 3, 4, 2 });
            p.AddCoeff(-1.5689528306600200e+03, new int[] { 3, 6, 0 });
            p.AddCoeff(4.7068584919800700e+03, new int[] { 3, 6, 2 });
            p.AddCoeff(-7.4712039555239300e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(2.2413611866571800e+02, new int[] { 5, 0, 2 });
            p.AddCoeff(1.5689528306600200e+03, new int[] { 5, 2, 0 });
            p.AddCoeff(-4.7068584919800700e+03, new int[] { 5, 2, 2 });
            p.AddCoeff(-4.7068584919800700e+03, new int[] { 5, 4, 0 });
            p.AddCoeff(1.4120575475940200e+04, new int[] { 5, 4, 2 });
            p.AddCoeff(3.4516962274520500e+03, new int[] { 5, 6, 0 });
            p.AddCoeff(-1.0355088682356200e+04, new int[] { 5, 6, 2 });
            p.AddCoeff(4.6250310200862400e+01, new int[] { 7, 0, 0 });
            p.AddCoeff(-1.3875093060258700e+02, new int[] { 7, 0, 2 });
            p.AddCoeff(-9.7125651421811000e+02, new int[] { 7, 2, 0 });
            p.AddCoeff(2.9137695426543300e+03, new int[] { 7, 2, 2 });
            p.AddCoeff(2.9137695426543300e+03, new int[] { 7, 4, 0 });
            p.AddCoeff(-8.7413086279629900e+03, new int[] { 7, 4, 2 });
            p.AddCoeff(-2.1367643312798400e+03, new int[] { 7, 6, 0 });
            p.AddCoeff(6.4102929938395300e+03, new int[] { 7, 6, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("6b45d646-da5d-4a48-a287-60e18e00c41a"));
            OrthonormalPolynomials[778] = p;
            p.AddCoeff(4.3954466819961800e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-3.9559020137965700e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(8.7029844303524500e+02, new int[] { 1, 5, 1 });
            p.AddCoeff(-5.3875617902181800e+02, new int[] { 1, 7, 1 });
            p.AddCoeff(-3.9559020137965700e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(3.5603118124169100e+03, new int[] { 3, 3, 1 });
            p.AddCoeff(-7.8326859873172000e+03, new int[] { 3, 5, 1 });
            p.AddCoeff(4.8488056111963600e+03, new int[] { 3, 7, 1 });
            p.AddCoeff(8.7029844303524500e+02, new int[] { 5, 1, 1 });
            p.AddCoeff(-7.8326859873172000e+03, new int[] { 5, 3, 1 });
            p.AddCoeff(1.7231909172097800e+04, new int[] { 5, 5, 1 });
            p.AddCoeff(-1.0667372344632000e+04, new int[] { 5, 7, 1 });
            p.AddCoeff(-5.3875617902181800e+02, new int[] { 7, 1, 1 });
            p.AddCoeff(4.8488056111963600e+03, new int[] { 7, 3, 1 });
            p.AddCoeff(-1.0667372344632000e+04, new int[] { 7, 5, 1 });
            p.AddCoeff(6.6036114514388600e+03, new int[] { 7, 7, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("242d40d7-3241-465c-8a64-f961b5eb4c48"));
            OrthonormalPolynomials[779] = p;
            p.AddCoeff(-3.3770013411936900e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(1.2157204828297300e+02, new int[] { 1, 2, 0 });
            p.AddCoeff(-6.6864626555635000e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(1.1589868602976700e+03, new int[] { 1, 6, 0 });
            p.AddCoeff(-6.2088581801661100e+02, new int[] { 1, 8, 0 });
            p.AddCoeff(3.0393012070743200e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-1.0941484345467500e+03, new int[] { 3, 2, 0 });
            p.AddCoeff(6.0178163900071500e+03, new int[] { 3, 4, 0 });
            p.AddCoeff(-1.0430881742679100e+04, new int[] { 3, 6, 0 });
            p.AddCoeff(5.5879723621494900e+03, new int[] { 3, 8, 0 });
            p.AddCoeff(-6.6864626555635000e+01, new int[] { 5, 0, 0 });
            p.AddCoeff(2.4071265560028600e+03, new int[] { 5, 2, 0 });
            p.AddCoeff(-1.3239196058015700e+04, new int[] { 5, 4, 0 });
            p.AddCoeff(2.2947939833893900e+04, new int[] { 5, 6, 0 });
            p.AddCoeff(-1.2293539196728900e+04, new int[] { 5, 8, 0 });
            p.AddCoeff(4.1392387867774100e+01, new int[] { 7, 0, 0 });
            p.AddCoeff(-1.4901259632398700e+03, new int[] { 7, 2, 0 });
            p.AddCoeff(8.1956927978192600e+03, new int[] { 7, 4, 0 });
            p.AddCoeff(-1.4205867516220100e+04, new int[] { 7, 6, 0 });
            p.AddCoeff(7.6102861694036100e+03, new int[] { 7, 8, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("da9b2b3b-62be-42da-83bf-407f52b5dd25"));
            OrthonormalPolynomials[780] = p;
            p.AddCoeff(-3.3770013411936900e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(3.0393012070743200e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-6.6864626555635000e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(4.1392387867774100e+01, new int[] { 0, 0, 7 });
            p.AddCoeff(1.2157204828297300e+02, new int[] { 2, 0, 1 });
            p.AddCoeff(-1.0941484345467500e+03, new int[] { 2, 0, 3 });
            p.AddCoeff(2.4071265560028600e+03, new int[] { 2, 0, 5 });
            p.AddCoeff(-1.4901259632398700e+03, new int[] { 2, 0, 7 });
            p.AddCoeff(-6.6864626555635000e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(6.0178163900071500e+03, new int[] { 4, 0, 3 });
            p.AddCoeff(-1.3239196058015700e+04, new int[] { 4, 0, 5 });
            p.AddCoeff(8.1956927978192600e+03, new int[] { 4, 0, 7 });
            p.AddCoeff(1.1589868602976700e+03, new int[] { 6, 0, 1 });
            p.AddCoeff(-1.0430881742679100e+04, new int[] { 6, 0, 3 });
            p.AddCoeff(2.2947939833893900e+04, new int[] { 6, 0, 5 });
            p.AddCoeff(-1.4205867516220100e+04, new int[] { 6, 0, 7 });
            p.AddCoeff(-6.2088581801661100e+02, new int[] { 8, 0, 1 });
            p.AddCoeff(5.5879723621494900e+03, new int[] { 8, 0, 3 });
            p.AddCoeff(-1.2293539196728900e+04, new int[] { 8, 0, 5 });
            p.AddCoeff(7.6102861694036100e+03, new int[] { 8, 0, 7 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("b09745ef-279f-40a0-ba4f-22c4aa2607e6"));
            OrthonormalPolynomials[781] = p;
            p.AddCoeff(-7.7789300654438300e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(1.6335753137432000e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-4.9007259412296100e+01, new int[] { 0, 1, 4 });
            p.AddCoeff(3.5938656902350500e+01, new int[] { 0, 1, 6 });
            p.AddCoeff(2.8004148235597800e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-5.8808711294755300e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(1.7642613388426600e+03, new int[] { 2, 1, 4 });
            p.AddCoeff(-1.2937916484846200e+03, new int[] { 2, 1, 6 });
            p.AddCoeff(-1.5402281529578800e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(3.2344791212115400e+03, new int[] { 4, 1, 2 });
            p.AddCoeff(-9.7034373636346300e+03, new int[] { 4, 1, 4 });
            p.AddCoeff(7.1158540666654000e+03, new int[] { 4, 1, 6 });
            p.AddCoeff(2.6697287984603200e+02, new int[] { 6, 1, 0 });
            p.AddCoeff(-5.6064304767666800e+03, new int[] { 6, 1, 2 });
            p.AddCoeff(1.6819291430300000e+04, new int[] { 6, 1, 4 });
            p.AddCoeff(-1.2334147048886700e+04, new int[] { 6, 1, 6 });
            p.AddCoeff(-1.4302118563180300e+02, new int[] { 8, 1, 0 });
            p.AddCoeff(3.0034448982678600e+03, new int[] { 8, 1, 2 });
            p.AddCoeff(-9.0103346948035700e+03, new int[] { 8, 1, 4 });
            p.AddCoeff(6.6075787761892800e+03, new int[] { 8, 1, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("39413393-f193-4f6e-886b-5ea371bca0a8"));
            OrthonormalPolynomials[782] = p;
            p.AddCoeff(-2.7713422517043000e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(1.2932930507953400e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-1.1639637457158100e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(8.3140267551129000e+00, new int[] { 0, 2, 1 });
            p.AddCoeff(-3.8798791523860200e+01, new int[] { 0, 2, 3 });
            p.AddCoeff(3.4918912371474200e+01, new int[] { 0, 2, 5 });
            p.AddCoeff(9.9768321061354800e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-4.6558549828632300e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(4.1902694845769000e+02, new int[] { 2, 0, 5 });
            p.AddCoeff(-2.9930496318406400e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(1.3967564948589700e+03, new int[] { 2, 2, 3 });
            p.AddCoeff(-1.2570808453730700e+03, new int[] { 2, 2, 5 });
            p.AddCoeff(-5.4872576583745200e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(2.5607202405747700e+03, new int[] { 4, 0, 3 });
            p.AddCoeff(-2.3046482165173000e+03, new int[] { 4, 0, 5 });
            p.AddCoeff(1.6461772975123500e+03, new int[] { 4, 2, 1 });
            p.AddCoeff(-7.6821607217243200e+03, new int[] { 4, 2, 3 });
            p.AddCoeff(6.9139446495518900e+03, new int[] { 4, 2, 5 });
            p.AddCoeff(9.5112466078491600e+02, new int[] { 6, 0, 1 });
            p.AddCoeff(-4.4385817503296100e+03, new int[] { 6, 0, 3 });
            p.AddCoeff(3.9947235752966500e+03, new int[] { 6, 0, 5 });
            p.AddCoeff(-2.8533739823547500e+03, new int[] { 6, 2, 1 });
            p.AddCoeff(1.3315745250988800e+04, new int[] { 6, 2, 3 });
            p.AddCoeff(-1.1984170725889900e+04, new int[] { 6, 2, 5 });
            p.AddCoeff(-5.0953106827763400e+02, new int[] { 8, 0, 1 });
            p.AddCoeff(2.3778116519622900e+03, new int[] { 8, 0, 3 });
            p.AddCoeff(-2.1400304867660600e+03, new int[] { 8, 0, 5 });
            p.AddCoeff(1.5285932048329000e+03, new int[] { 8, 2, 1 });
            p.AddCoeff(-7.1334349558868700e+03, new int[] { 8, 2, 3 });
            p.AddCoeff(6.4200914602982000e+03, new int[] { 8, 2, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("0458d688-5a7c-4a48-860d-e86801fb09f6"));
            OrthonormalPolynomials[783] = p;
            p.AddCoeff(-1.7796325618178400e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(1.7796325618178400e+01, new int[] { 0, 1, 2 });
            p.AddCoeff(-2.0762379887874800e+01, new int[] { 0, 1, 4 });
            p.AddCoeff(2.9660542696963900e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(-2.9660542696963900e+01, new int[] { 0, 3, 2 });
            p.AddCoeff(3.4603966479791300e+01, new int[] { 0, 3, 4 });
            p.AddCoeff(6.4066772225442100e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-6.4066772225442100e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(7.4744567596349100e+02, new int[] { 2, 1, 4 });
            p.AddCoeff(-1.0677795370907000e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(1.0677795370907000e+03, new int[] { 2, 3, 2 });
            p.AddCoeff(-1.2457427932724900e+03, new int[] { 2, 3, 4 });
            p.AddCoeff(-3.5236724723993200e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(3.5236724723993200e+03, new int[] { 4, 1, 2 });
            p.AddCoeff(-4.1109512177992000e+03, new int[] { 4, 1, 4 });
            p.AddCoeff(5.8727874539988600e+02, new int[] { 4, 3, 0 });
            p.AddCoeff(-5.8727874539988600e+03, new int[] { 4, 3, 2 });
            p.AddCoeff(6.8515853629986700e+03, new int[] { 4, 3, 4 });
            p.AddCoeff(6.1076989521588200e+02, new int[] { 6, 1, 0 });
            p.AddCoeff(-6.1076989521588200e+03, new int[] { 6, 1, 2 });
            p.AddCoeff(7.1256487775186200e+03, new int[] { 6, 1, 4 });
            p.AddCoeff(-1.0179498253598000e+03, new int[] { 6, 3, 0 });
            p.AddCoeff(1.0179498253598000e+04, new int[] { 6, 3, 2 });
            p.AddCoeff(-1.1876081295864400e+04, new int[] { 6, 3, 4 });
            p.AddCoeff(-3.2719815815136500e+02, new int[] { 8, 1, 0 });
            p.AddCoeff(3.2719815815136500e+03, new int[] { 8, 1, 2 });
            p.AddCoeff(-3.8173118450992500e+03, new int[] { 8, 1, 4 });
            p.AddCoeff(5.4533026358560900e+02, new int[] { 8, 3, 0 });
            p.AddCoeff(-5.4533026358560900e+03, new int[] { 8, 3, 2 });
            p.AddCoeff(6.3621864084987700e+03, new int[] { 8, 3, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("2e4c18b9-c911-4f75-a5ec-b3c2b8875c7e"));
            OrthonormalPolynomials[784] = p;
            p.AddCoeff(-1.7796325618178400e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(2.9660542696963900e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(1.7796325618178400e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-2.9660542696963900e+01, new int[] { 0, 2, 3 });
            p.AddCoeff(-2.0762379887874800e+01, new int[] { 0, 4, 1 });
            p.AddCoeff(3.4603966479791300e+01, new int[] { 0, 4, 3 });
            p.AddCoeff(6.4066772225442100e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-1.0677795370907000e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(-6.4066772225442100e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(1.0677795370907000e+03, new int[] { 2, 2, 3 });
            p.AddCoeff(7.4744567596349100e+02, new int[] { 2, 4, 1 });
            p.AddCoeff(-1.2457427932724900e+03, new int[] { 2, 4, 3 });
            p.AddCoeff(-3.5236724723993200e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(5.8727874539988600e+02, new int[] { 4, 0, 3 });
            p.AddCoeff(3.5236724723993200e+03, new int[] { 4, 2, 1 });
            p.AddCoeff(-5.8727874539988600e+03, new int[] { 4, 2, 3 });
            p.AddCoeff(-4.1109512177992000e+03, new int[] { 4, 4, 1 });
            p.AddCoeff(6.8515853629986700e+03, new int[] { 4, 4, 3 });
            p.AddCoeff(6.1076989521588200e+02, new int[] { 6, 0, 1 });
            p.AddCoeff(-1.0179498253598000e+03, new int[] { 6, 0, 3 });
            p.AddCoeff(-6.1076989521588200e+03, new int[] { 6, 2, 1 });
            p.AddCoeff(1.0179498253598000e+04, new int[] { 6, 2, 3 });
            p.AddCoeff(7.1256487775186200e+03, new int[] { 6, 4, 1 });
            p.AddCoeff(-1.1876081295864400e+04, new int[] { 6, 4, 3 });
            p.AddCoeff(-3.2719815815136500e+02, new int[] { 8, 0, 1 });
            p.AddCoeff(5.4533026358560900e+02, new int[] { 8, 0, 3 });
            p.AddCoeff(3.2719815815136500e+03, new int[] { 8, 2, 1 });
            p.AddCoeff(-5.4533026358560900e+03, new int[] { 8, 2, 3 });
            p.AddCoeff(-3.8173118450992500e+03, new int[] { 8, 4, 1 });
            p.AddCoeff(6.3621864084987700e+03, new int[] { 8, 4, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("1f2a036a-9b44-45c6-9fde-519dd6f74b32"));
            OrthonormalPolynomials[785] = p;
            p.AddCoeff(-2.7713422517043000e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(8.3140267551129000e+00, new int[] { 0, 1, 2 });
            p.AddCoeff(1.2932930507953400e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-3.8798791523860200e+01, new int[] { 0, 3, 2 });
            p.AddCoeff(-1.1639637457158100e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(3.4918912371474200e+01, new int[] { 0, 5, 2 });
            p.AddCoeff(9.9768321061354800e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-2.9930496318406400e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(-4.6558549828632300e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(1.3967564948589700e+03, new int[] { 2, 3, 2 });
            p.AddCoeff(4.1902694845769000e+02, new int[] { 2, 5, 0 });
            p.AddCoeff(-1.2570808453730700e+03, new int[] { 2, 5, 2 });
            p.AddCoeff(-5.4872576583745200e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(1.6461772975123500e+03, new int[] { 4, 1, 2 });
            p.AddCoeff(2.5607202405747700e+03, new int[] { 4, 3, 0 });
            p.AddCoeff(-7.6821607217243200e+03, new int[] { 4, 3, 2 });
            p.AddCoeff(-2.3046482165173000e+03, new int[] { 4, 5, 0 });
            p.AddCoeff(6.9139446495518900e+03, new int[] { 4, 5, 2 });
            p.AddCoeff(9.5112466078491600e+02, new int[] { 6, 1, 0 });
            p.AddCoeff(-2.8533739823547500e+03, new int[] { 6, 1, 2 });
            p.AddCoeff(-4.4385817503296100e+03, new int[] { 6, 3, 0 });
            p.AddCoeff(1.3315745250988800e+04, new int[] { 6, 3, 2 });
            p.AddCoeff(3.9947235752966500e+03, new int[] { 6, 5, 0 });
            p.AddCoeff(-1.1984170725889900e+04, new int[] { 6, 5, 2 });
            p.AddCoeff(-5.0953106827763400e+02, new int[] { 8, 1, 0 });
            p.AddCoeff(1.5285932048329000e+03, new int[] { 8, 1, 2 });
            p.AddCoeff(2.3778116519622900e+03, new int[] { 8, 3, 0 });
            p.AddCoeff(-7.1334349558868700e+03, new int[] { 8, 3, 2 });
            p.AddCoeff(-2.1400304867660600e+03, new int[] { 8, 5, 0 });
            p.AddCoeff(6.4200914602982000e+03, new int[] { 8, 5, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a40d3f14-9400-492f-9601-e27b50915c42"));
            OrthonormalPolynomials[786] = p;
            p.AddCoeff(-7.7789300654438300e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(1.6335753137432000e+01, new int[] { 0, 2, 1 });
            p.AddCoeff(-4.9007259412296100e+01, new int[] { 0, 4, 1 });
            p.AddCoeff(3.5938656902350500e+01, new int[] { 0, 6, 1 });
            p.AddCoeff(2.8004148235597800e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-5.8808711294755300e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(1.7642613388426600e+03, new int[] { 2, 4, 1 });
            p.AddCoeff(-1.2937916484846200e+03, new int[] { 2, 6, 1 });
            p.AddCoeff(-1.5402281529578800e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(3.2344791212115400e+03, new int[] { 4, 2, 1 });
            p.AddCoeff(-9.7034373636346300e+03, new int[] { 4, 4, 1 });
            p.AddCoeff(7.1158540666654000e+03, new int[] { 4, 6, 1 });
            p.AddCoeff(2.6697287984603200e+02, new int[] { 6, 0, 1 });
            p.AddCoeff(-5.6064304767666800e+03, new int[] { 6, 2, 1 });
            p.AddCoeff(1.6819291430300000e+04, new int[] { 6, 4, 1 });
            p.AddCoeff(-1.2334147048886700e+04, new int[] { 6, 6, 1 });
            p.AddCoeff(-1.4302118563180300e+02, new int[] { 8, 0, 1 });
            p.AddCoeff(3.0034448982678600e+03, new int[] { 8, 2, 1 });
            p.AddCoeff(-9.0103346948035700e+03, new int[] { 8, 4, 1 });
            p.AddCoeff(6.6075787761892800e+03, new int[] { 8, 6, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("94123145-23fe-48c6-82b6-f4bcb6949b32"));
            OrthonormalPolynomials[787] = p;
            p.AddCoeff(-3.3770013411936900e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(3.0393012070743200e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-6.6864626555635000e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(4.1392387867774100e+01, new int[] { 0, 7, 0 });
            p.AddCoeff(1.2157204828297300e+02, new int[] { 2, 1, 0 });
            p.AddCoeff(-1.0941484345467500e+03, new int[] { 2, 3, 0 });
            p.AddCoeff(2.4071265560028600e+03, new int[] { 2, 5, 0 });
            p.AddCoeff(-1.4901259632398700e+03, new int[] { 2, 7, 0 });
            p.AddCoeff(-6.6864626555635000e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(6.0178163900071500e+03, new int[] { 4, 3, 0 });
            p.AddCoeff(-1.3239196058015700e+04, new int[] { 4, 5, 0 });
            p.AddCoeff(8.1956927978192600e+03, new int[] { 4, 7, 0 });
            p.AddCoeff(1.1589868602976700e+03, new int[] { 6, 1, 0 });
            p.AddCoeff(-1.0430881742679100e+04, new int[] { 6, 3, 0 });
            p.AddCoeff(2.2947939833893900e+04, new int[] { 6, 5, 0 });
            p.AddCoeff(-1.4205867516220100e+04, new int[] { 6, 7, 0 });
            p.AddCoeff(-6.2088581801661100e+02, new int[] { 8, 1, 0 });
            p.AddCoeff(5.5879723621494900e+03, new int[] { 8, 3, 0 });
            p.AddCoeff(-1.2293539196728900e+04, new int[] { 8, 5, 0 });
            p.AddCoeff(7.6102861694036100e+03, new int[] { 8, 7, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("3ffa71db-e82d-4f9b-9221-8bb04a9282b7"));
            OrthonormalPolynomials[788] = p;
            p.AddCoeff(-4.2732085527534600e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(8.9737379607822600e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-2.6921213882346800e+02, new int[] { 1, 0, 4 });
            p.AddCoeff(1.9742223513721000e+02, new int[] { 1, 0, 6 });
            p.AddCoeff(6.2673725440384000e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-1.3161482342480600e+03, new int[] { 3, 0, 2 });
            p.AddCoeff(3.9484447027441900e+03, new int[] { 3, 0, 4 });
            p.AddCoeff(-2.8955261153457400e+03, new int[] { 3, 0, 6 });
            p.AddCoeff(-2.4442752921749800e+02, new int[] { 5, 0, 0 });
            p.AddCoeff(5.1329781135674500e+03, new int[] { 5, 0, 2 });
            p.AddCoeff(-1.5398934340702400e+04, new int[] { 5, 0, 4 });
            p.AddCoeff(1.1292551849848400e+04, new int[] { 5, 0, 6 });
            p.AddCoeff(3.4918218459642500e+02, new int[] { 7, 0, 0 });
            p.AddCoeff(-7.3328258765249300e+03, new int[] { 7, 0, 2 });
            p.AddCoeff(2.1998477629574800e+04, new int[] { 7, 0, 4 });
            p.AddCoeff(-1.6132216928354800e+04, new int[] { 7, 0, 6 });
            p.AddCoeff(-1.6489158717053400e+02, new int[] { 9, 0, 0 });
            p.AddCoeff(3.4627233305812300e+03, new int[] { 9, 0, 2 });
            p.AddCoeff(-1.0388169991743600e+04, new int[] { 9, 0, 4 });
            p.AddCoeff(7.6179913272786700e+03, new int[] { 9, 0, 6 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("0d7445c7-b20d-4d9f-b5e2-cd6529ddca64"));
            OrthonormalPolynomials[789] = p;
            p.AddCoeff(4.0849865705801600e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-1.9063270662707400e+02, new int[] { 1, 1, 3 });
            p.AddCoeff(1.7156943596436700e+02, new int[] { 1, 1, 5 });
            p.AddCoeff(-5.9913136368509000e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(2.7959463638637500e+03, new int[] { 3, 1, 3 });
            p.AddCoeff(-2.5163517274773800e+03, new int[] { 3, 1, 5 });
            p.AddCoeff(2.3366123183718500e+03, new int[] { 5, 1, 1 });
            p.AddCoeff(-1.0904190819068600e+04, new int[] { 5, 1, 3 });
            p.AddCoeff(9.8137717371617800e+03, new int[] { 5, 1, 5 });
            p.AddCoeff(-3.3380175976740700e+03, new int[] { 7, 1, 1 });
            p.AddCoeff(1.5577415455812300e+04, new int[] { 7, 1, 3 });
            p.AddCoeff(-1.4019673910231100e+04, new int[] { 7, 1, 5 });
            p.AddCoeff(1.5762860877905300e+03, new int[] { 9, 1, 1 });
            p.AddCoeff(-7.3560017430224900e+03, new int[] { 9, 1, 3 });
            p.AddCoeff(6.6204015687202500e+03, new int[] { 9, 1, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("626bed20-a8c2-4014-8ad8-858b3fbeff4b"));
            OrthonormalPolynomials[790] = p;
            p.AddCoeff(-4.7702365981659800e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(4.7702365981659800e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-5.5652760311936400e+01, new int[] { 1, 0, 4 });
            p.AddCoeff(1.4310709794497900e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-1.4310709794497900e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(1.6695828093580900e+02, new int[] { 1, 2, 4 });
            p.AddCoeff(6.9963470106434400e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-6.9963470106434400e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(8.1624048457506800e+02, new int[] { 3, 0, 4 });
            p.AddCoeff(-2.0989041031930300e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(2.0989041031930300e+03, new int[] { 3, 2, 2 });
            p.AddCoeff(-2.4487214537252000e+03, new int[] { 3, 2, 4 });
            p.AddCoeff(-2.7285753341509400e+02, new int[] { 5, 0, 0 });
            p.AddCoeff(2.7285753341509400e+03, new int[] { 5, 0, 2 });
            p.AddCoeff(-3.1833378898427600e+03, new int[] { 5, 0, 4 });
            p.AddCoeff(8.1857260024528200e+02, new int[] { 5, 2, 0 });
            p.AddCoeff(-8.1857260024528200e+03, new int[] { 5, 2, 2 });
            p.AddCoeff(9.5500136695282900e+03, new int[] { 5, 2, 4 });
            p.AddCoeff(3.8979647630727700e+02, new int[] { 7, 0, 0 });
            p.AddCoeff(-3.8979647630727700e+03, new int[] { 7, 0, 2 });
            p.AddCoeff(4.5476255569182400e+03, new int[] { 7, 0, 4 });
            p.AddCoeff(-1.1693894289218300e+03, new int[] { 7, 2, 0 });
            p.AddCoeff(1.1693894289218300e+04, new int[] { 7, 2, 2 });
            p.AddCoeff(-1.3642876670754700e+04, new int[] { 7, 2, 4 });
            p.AddCoeff(-1.8407055825621400e+02, new int[] { 9, 0, 0 });
            p.AddCoeff(1.8407055825621400e+03, new int[] { 9, 0, 2 });
            p.AddCoeff(-2.1474898463225100e+03, new int[] { 9, 0, 4 });
            p.AddCoeff(5.5221167476864300e+02, new int[] { 9, 2, 0 });
            p.AddCoeff(-5.5221167476864300e+03, new int[] { 9, 2, 2 });
            p.AddCoeff(6.4424695389674900e+03, new int[] { 9, 2, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("95f34bfb-39e3-4219-93da-50c2c48afb5e"));
            OrthonormalPolynomials[791] = p;
            p.AddCoeff(5.9732810492636400e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-9.9554684154394000e+01, new int[] { 1, 1, 3 });
            p.AddCoeff(-9.9554684154394000e+01, new int[] { 1, 3, 1 });
            p.AddCoeff(1.6592447359065700e+02, new int[] { 1, 3, 3 });
            p.AddCoeff(-8.7608122055866800e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(1.4601353675977800e+03, new int[] { 3, 1, 3 });
            p.AddCoeff(1.4601353675977800e+03, new int[] { 3, 3, 1 });
            p.AddCoeff(-2.4335589459963000e+03, new int[] { 3, 3, 3 });
            p.AddCoeff(3.4167167601788000e+03, new int[] { 5, 1, 1 });
            p.AddCoeff(-5.6945279336313400e+03, new int[] { 5, 1, 3 });
            p.AddCoeff(-5.6945279336313400e+03, new int[] { 5, 3, 1 });
            p.AddCoeff(9.4908798893855700e+03, new int[] { 5, 3, 3 });
            p.AddCoeff(-4.8810239431125800e+03, new int[] { 7, 1, 1 });
            p.AddCoeff(8.1350399051876300e+03, new int[] { 7, 1, 3 });
            p.AddCoeff(8.1350399051876300e+03, new int[] { 7, 3, 1 });
            p.AddCoeff(-1.3558399841979400e+04, new int[] { 7, 3, 3 });
            p.AddCoeff(2.3049279731364900e+03, new int[] { 9, 1, 1 });
            p.AddCoeff(-3.8415466218941600e+03, new int[] { 9, 1, 3 });
            p.AddCoeff(-3.8415466218941600e+03, new int[] { 9, 3, 1 });
            p.AddCoeff(6.4025777031569300e+03, new int[] { 9, 3, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c0f0a90b-f5f9-4504-a990-b339f7fe5869"));
            OrthonormalPolynomials[792] = p;
            p.AddCoeff(-4.7702365981659800e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(1.4310709794497900e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(4.7702365981659800e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-1.4310709794497900e+02, new int[] { 1, 2, 2 });
            p.AddCoeff(-5.5652760311936400e+01, new int[] { 1, 4, 0 });
            p.AddCoeff(1.6695828093580900e+02, new int[] { 1, 4, 2 });
            p.AddCoeff(6.9963470106434400e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-2.0989041031930300e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(-6.9963470106434400e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(2.0989041031930300e+03, new int[] { 3, 2, 2 });
            p.AddCoeff(8.1624048457506800e+02, new int[] { 3, 4, 0 });
            p.AddCoeff(-2.4487214537252000e+03, new int[] { 3, 4, 2 });
            p.AddCoeff(-2.7285753341509400e+02, new int[] { 5, 0, 0 });
            p.AddCoeff(8.1857260024528200e+02, new int[] { 5, 0, 2 });
            p.AddCoeff(2.7285753341509400e+03, new int[] { 5, 2, 0 });
            p.AddCoeff(-8.1857260024528200e+03, new int[] { 5, 2, 2 });
            p.AddCoeff(-3.1833378898427600e+03, new int[] { 5, 4, 0 });
            p.AddCoeff(9.5500136695282900e+03, new int[] { 5, 4, 2 });
            p.AddCoeff(3.8979647630727700e+02, new int[] { 7, 0, 0 });
            p.AddCoeff(-1.1693894289218300e+03, new int[] { 7, 0, 2 });
            p.AddCoeff(-3.8979647630727700e+03, new int[] { 7, 2, 0 });
            p.AddCoeff(1.1693894289218300e+04, new int[] { 7, 2, 2 });
            p.AddCoeff(4.5476255569182400e+03, new int[] { 7, 4, 0 });
            p.AddCoeff(-1.3642876670754700e+04, new int[] { 7, 4, 2 });
            p.AddCoeff(-1.8407055825621400e+02, new int[] { 9, 0, 0 });
            p.AddCoeff(5.5221167476864300e+02, new int[] { 9, 0, 2 });
            p.AddCoeff(1.8407055825621400e+03, new int[] { 9, 2, 0 });
            p.AddCoeff(-5.5221167476864300e+03, new int[] { 9, 2, 2 });
            p.AddCoeff(-2.1474898463225100e+03, new int[] { 9, 4, 0 });
            p.AddCoeff(6.4424695389674900e+03, new int[] { 9, 4, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e640f1e5-7dac-446b-bf63-1777bd7f7642"));
            OrthonormalPolynomials[793] = p;
            p.AddCoeff(4.0849865705801600e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-1.9063270662707400e+02, new int[] { 1, 3, 1 });
            p.AddCoeff(1.7156943596436700e+02, new int[] { 1, 5, 1 });
            p.AddCoeff(-5.9913136368509000e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(2.7959463638637500e+03, new int[] { 3, 3, 1 });
            p.AddCoeff(-2.5163517274773800e+03, new int[] { 3, 5, 1 });
            p.AddCoeff(2.3366123183718500e+03, new int[] { 5, 1, 1 });
            p.AddCoeff(-1.0904190819068600e+04, new int[] { 5, 3, 1 });
            p.AddCoeff(9.8137717371617800e+03, new int[] { 5, 5, 1 });
            p.AddCoeff(-3.3380175976740700e+03, new int[] { 7, 1, 1 });
            p.AddCoeff(1.5577415455812300e+04, new int[] { 7, 3, 1 });
            p.AddCoeff(-1.4019673910231100e+04, new int[] { 7, 5, 1 });
            p.AddCoeff(1.5762860877905300e+03, new int[] { 9, 1, 1 });
            p.AddCoeff(-7.3560017430224900e+03, new int[] { 9, 3, 1 });
            p.AddCoeff(6.6204015687202500e+03, new int[] { 9, 5, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("09245222-3f91-4bcd-b343-7527eb2fffc4"));
            OrthonormalPolynomials[794] = p;
            p.AddCoeff(-4.2732085527534600e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(8.9737379607822600e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-2.6921213882346800e+02, new int[] { 1, 4, 0 });
            p.AddCoeff(1.9742223513721000e+02, new int[] { 1, 6, 0 });
            p.AddCoeff(6.2673725440384000e+01, new int[] { 3, 0, 0 });
            p.AddCoeff(-1.3161482342480600e+03, new int[] { 3, 2, 0 });
            p.AddCoeff(3.9484447027441900e+03, new int[] { 3, 4, 0 });
            p.AddCoeff(-2.8955261153457400e+03, new int[] { 3, 6, 0 });
            p.AddCoeff(-2.4442752921749800e+02, new int[] { 5, 0, 0 });
            p.AddCoeff(5.1329781135674500e+03, new int[] { 5, 2, 0 });
            p.AddCoeff(-1.5398934340702400e+04, new int[] { 5, 4, 0 });
            p.AddCoeff(1.1292551849848400e+04, new int[] { 5, 6, 0 });
            p.AddCoeff(3.4918218459642500e+02, new int[] { 7, 0, 0 });
            p.AddCoeff(-7.3328258765249300e+03, new int[] { 7, 2, 0 });
            p.AddCoeff(2.1998477629574800e+04, new int[] { 7, 4, 0 });
            p.AddCoeff(-1.6132216928354800e+04, new int[] { 7, 6, 0 });
            p.AddCoeff(-1.6489158717053400e+02, new int[] { 9, 0, 0 });
            p.AddCoeff(3.4627233305812300e+03, new int[] { 9, 2, 0 });
            p.AddCoeff(-1.0388169991743600e+04, new int[] { 9, 4, 0 });
            p.AddCoeff(7.6179913272786700e+03, new int[] { 9, 6, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("aae883f5-f06c-4a36-8c10-32a2dab91876"));
            OrthonormalPolynomials[795] = p;
            p.AddCoeff(-2.4794928065055500e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(1.1570966430359200e+01, new int[] { 0, 0, 3 });
            p.AddCoeff(-1.0413869787323300e+01, new int[] { 0, 0, 5 });
            p.AddCoeff(1.3637210435780500e+02, new int[] { 2, 0, 1 });
            p.AddCoeff(-6.3640315366975700e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(5.7276283830278100e+02, new int[] { 2, 0, 5 });
            p.AddCoeff(-1.1818915711009800e+03, new int[] { 4, 0, 1 });
            p.AddCoeff(5.5154939984712300e+03, new int[] { 4, 0, 3 });
            p.AddCoeff(-4.9639445986241100e+03, new int[] { 4, 0, 5 });
            p.AddCoeff(3.5456747133029300e+03, new int[] { 6, 0, 1 });
            p.AddCoeff(-1.6546481995413700e+04, new int[] { 6, 0, 3 });
            p.AddCoeff(1.4891833795872300e+04, new int[] { 6, 0, 5 });
            p.AddCoeff(-4.3054621518678400e+03, new int[] { 8, 0, 1 });
            p.AddCoeff(2.0092156708716600e+04, new int[] { 8, 0, 3 });
            p.AddCoeff(-1.8082941037844900e+04, new int[] { 8, 0, 5 });
            p.AddCoeff(1.8178617974553100e+03, new int[] { 10, 0, 1 });
            p.AddCoeff(-8.4833550547914600e+03, new int[] { 10, 0, 3 });
            p.AddCoeff(7.6350195493123200e+03, new int[] { 10, 0, 5 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("73be24a5-a724-42a2-9ece-7f5f869211b1"));
            OrthonormalPolynomials[796] = p;
            p.AddCoeff(-7.7692373228789900e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(7.7692373228789900e+00, new int[] { 0, 1, 2 });
            p.AddCoeff(-9.0641102100254900e+00, new int[] { 0, 1, 4 });
            p.AddCoeff(4.2730805275834500e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-4.2730805275834500e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(4.9852606155140400e+02, new int[] { 2, 1, 4 });
            p.AddCoeff(-3.7033364572389900e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(3.7033364572389900e+03, new int[] { 4, 1, 2 });
            p.AddCoeff(-4.3205592001121300e+03, new int[] { 4, 1, 4 });
            p.AddCoeff(1.1110009371717000e+03, new int[] { 6, 1, 0 });
            p.AddCoeff(-1.1110009371717000e+04, new int[] { 6, 1, 2 });
            p.AddCoeff(1.2961677600336500e+04, new int[] { 6, 1, 4 });
            p.AddCoeff(-1.3490725665656300e+03, new int[] { 8, 1, 0 });
            p.AddCoeff(1.3490725665656300e+04, new int[] { 8, 1, 2 });
            p.AddCoeff(-1.5739179943265700e+04, new int[] { 8, 1, 4 });
            p.AddCoeff(5.6960841699437900e+02, new int[] { 10, 1, 0 });
            p.AddCoeff(-5.6960841699437900e+03, new int[] { 10, 1, 2 });
            p.AddCoeff(6.6454315316010800e+03, new int[] { 10, 1, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("bac53749-bbda-4f01-9ebe-a4229483e1e7"));
            OrthonormalPolynomials[797] = p;
            p.AddCoeff(-1.7691331630354800e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(2.9485552717258100e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(5.3073994891064500e+00, new int[] { 0, 2, 1 });
            p.AddCoeff(-8.8456658151774200e+00, new int[] { 0, 2, 3 });
            p.AddCoeff(9.7302323966951700e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-1.6217053994491900e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(-2.9190697190085500e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(4.8651161983475800e+02, new int[] { 2, 2, 3 });
            p.AddCoeff(-8.4328680771358100e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(1.4054780128559700e+03, new int[] { 4, 0, 3 });
            p.AddCoeff(2.5298604231407400e+03, new int[] { 4, 2, 1 });
            p.AddCoeff(-4.2164340385679100e+03, new int[] { 4, 2, 3 });
            p.AddCoeff(2.5298604231407400e+03, new int[] { 6, 0, 1 });
            p.AddCoeff(-4.2164340385679100e+03, new int[] { 6, 0, 3 });
            p.AddCoeff(-7.5895812694222300e+03, new int[] { 6, 2, 1 });
            p.AddCoeff(1.2649302115703700e+04, new int[] { 6, 2, 3 });
            p.AddCoeff(-3.0719733709566200e+03, new int[] { 8, 0, 1 });
            p.AddCoeff(5.1199556182610300e+03, new int[] { 8, 0, 3 });
            p.AddCoeff(9.2159201128698200e+03, new int[] { 8, 2, 1 });
            p.AddCoeff(-1.5359866854783100e+04, new int[] { 8, 2, 3 });
            p.AddCoeff(1.2970554232927900e+03, new int[] { 10, 0, 1 });
            p.AddCoeff(-2.1617590388213200e+03, new int[] { 10, 0, 3 });
            p.AddCoeff(-3.8911662698783800e+03, new int[] { 10, 2, 1 });
            p.AddCoeff(6.4852771164639400e+03, new int[] { 10, 2, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("c4355044-611b-470b-81de-bdedd1a17f43"));
            OrthonormalPolynomials[798] = p;
            p.AddCoeff(-1.7691331630354800e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(5.3073994891064500e+00, new int[] { 0, 1, 2 });
            p.AddCoeff(2.9485552717258100e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(-8.8456658151774200e+00, new int[] { 0, 3, 2 });
            p.AddCoeff(9.7302323966951700e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-2.9190697190085500e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(-1.6217053994491900e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(4.8651161983475800e+02, new int[] { 2, 3, 2 });
            p.AddCoeff(-8.4328680771358100e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(2.5298604231407400e+03, new int[] { 4, 1, 2 });
            p.AddCoeff(1.4054780128559700e+03, new int[] { 4, 3, 0 });
            p.AddCoeff(-4.2164340385679100e+03, new int[] { 4, 3, 2 });
            p.AddCoeff(2.5298604231407400e+03, new int[] { 6, 1, 0 });
            p.AddCoeff(-7.5895812694222300e+03, new int[] { 6, 1, 2 });
            p.AddCoeff(-4.2164340385679100e+03, new int[] { 6, 3, 0 });
            p.AddCoeff(1.2649302115703700e+04, new int[] { 6, 3, 2 });
            p.AddCoeff(-3.0719733709566200e+03, new int[] { 8, 1, 0 });
            p.AddCoeff(9.2159201128698200e+03, new int[] { 8, 1, 2 });
            p.AddCoeff(5.1199556182610300e+03, new int[] { 8, 3, 0 });
            p.AddCoeff(-1.5359866854783100e+04, new int[] { 8, 3, 2 });
            p.AddCoeff(1.2970554232927900e+03, new int[] { 10, 1, 0 });
            p.AddCoeff(-3.8911662698783800e+03, new int[] { 10, 1, 2 });
            p.AddCoeff(-2.1617590388213200e+03, new int[] { 10, 3, 0 });
            p.AddCoeff(6.4852771164639400e+03, new int[] { 10, 3, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("4eb743e1-bdfc-485c-ae40-7e381276e78d"));
            OrthonormalPolynomials[799] = p;
            p.AddCoeff(-7.7692373228789900e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(7.7692373228789900e+00, new int[] { 0, 2, 1 });
            p.AddCoeff(-9.0641102100254900e+00, new int[] { 0, 4, 1 });
            p.AddCoeff(4.2730805275834500e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-4.2730805275834500e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(4.9852606155140400e+02, new int[] { 2, 4, 1 });
            p.AddCoeff(-3.7033364572389900e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(3.7033364572389900e+03, new int[] { 4, 2, 1 });
            p.AddCoeff(-4.3205592001121300e+03, new int[] { 4, 4, 1 });
            p.AddCoeff(1.1110009371717000e+03, new int[] { 6, 0, 1 });
            p.AddCoeff(-1.1110009371717000e+04, new int[] { 6, 2, 1 });
            p.AddCoeff(1.2961677600336500e+04, new int[] { 6, 4, 1 });
            p.AddCoeff(-1.3490725665656300e+03, new int[] { 8, 0, 1 });
            p.AddCoeff(1.3490725665656300e+04, new int[] { 8, 2, 1 });
            p.AddCoeff(-1.5739179943265700e+04, new int[] { 8, 4, 1 });
            p.AddCoeff(5.6960841699437900e+02, new int[] { 10, 0, 1 });
            p.AddCoeff(-5.6960841699437900e+03, new int[] { 10, 2, 1 });
            p.AddCoeff(6.6454315316010800e+03, new int[] { 10, 4, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("a2f8c405-cc37-4e71-9e1f-e5cad164fe20"));
            OrthonormalPolynomials[800] = p;
            p.AddCoeff(-2.4794928065055500e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(1.1570966430359200e+01, new int[] { 0, 3, 0 });
            p.AddCoeff(-1.0413869787323300e+01, new int[] { 0, 5, 0 });
            p.AddCoeff(1.3637210435780500e+02, new int[] { 2, 1, 0 });
            p.AddCoeff(-6.3640315366975700e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(5.7276283830278100e+02, new int[] { 2, 5, 0 });
            p.AddCoeff(-1.1818915711009800e+03, new int[] { 4, 1, 0 });
            p.AddCoeff(5.5154939984712300e+03, new int[] { 4, 3, 0 });
            p.AddCoeff(-4.9639445986241100e+03, new int[] { 4, 5, 0 });
            p.AddCoeff(3.5456747133029300e+03, new int[] { 6, 1, 0 });
            p.AddCoeff(-1.6546481995413700e+04, new int[] { 6, 3, 0 });
            p.AddCoeff(1.4891833795872300e+04, new int[] { 6, 5, 0 });
            p.AddCoeff(-4.3054621518678400e+03, new int[] { 8, 1, 0 });
            p.AddCoeff(2.0092156708716600e+04, new int[] { 8, 3, 0 });
            p.AddCoeff(-1.8082941037844900e+04, new int[] { 8, 5, 0 });
            p.AddCoeff(1.8178617974553100e+03, new int[] { 10, 1, 0 });
            p.AddCoeff(-8.4833550547914600e+03, new int[] { 10, 3, 0 });
            p.AddCoeff(7.6350195493123200e+03, new int[] { 10, 5, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("798c4749-ede8-4280-bfca-7bd0e2a11064"));
            OrthonormalPolynomials[801] = p;
            p.AddCoeff(-5.1637441534121500e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(5.1637441534121500e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(-6.0243681789808400e+01, new int[] { 1, 0, 4 });
            p.AddCoeff(1.1188112332393000e+02, new int[] { 3, 0, 0 });
            p.AddCoeff(-1.1188112332393000e+03, new int[] { 3, 0, 2 });
            p.AddCoeff(1.3052797721125200e+03, new int[] { 3, 0, 4 });
            p.AddCoeff(-6.7128673994357900e+02, new int[] { 5, 0, 0 });
            p.AddCoeff(6.7128673994357900e+03, new int[] { 5, 0, 2 });
            p.AddCoeff(-7.8316786326750600e+03, new int[] { 5, 0, 4 });
            p.AddCoeff(1.6302677970058400e+03, new int[] { 7, 0, 0 });
            p.AddCoeff(-1.6302677970058400e+04, new int[] { 7, 0, 2 });
            p.AddCoeff(1.9019790965068100e+04, new int[] { 7, 0, 4 });
            p.AddCoeff(-1.7208382301728200e+03, new int[] { 9, 0, 0 });
            p.AddCoeff(1.7208382301728200e+04, new int[] { 9, 0, 2 });
            p.AddCoeff(-2.0076446018683000e+04, new int[] { 9, 0, 4 });
            p.AddCoeff(6.5704732424780600e+02, new int[] { 11, 0, 0 });
            p.AddCoeff(-6.5704732424780600e+03, new int[] { 11, 0, 2 });
            p.AddCoeff(7.6655521162244200e+03, new int[] { 11, 0, 4 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("f95b2f80-e5dc-4992-a025-ec60d0fa323d"));
            OrthonormalPolynomials[802] = p;
            p.AddCoeff(3.1550997936529100e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-5.2584996560881800e+01, new int[] { 1, 1, 3 });
            p.AddCoeff(-6.8360495529146400e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(1.1393415921524400e+03, new int[] { 3, 1, 3 });
            p.AddCoeff(4.1016297317487800e+03, new int[] { 5, 1, 1 });
            p.AddCoeff(-6.8360495529146400e+03, new int[] { 5, 1, 3 });
            p.AddCoeff(-9.9611007771041800e+03, new int[] { 7, 1, 1 });
            p.AddCoeff(1.6601834628507000e+04, new int[] { 7, 1, 3 });
            p.AddCoeff(1.0514495264721100e+04, new int[] { 9, 1, 1 });
            p.AddCoeff(-1.7524158774535100e+04, new int[] { 9, 1, 3 });
            p.AddCoeff(-4.0146254647116900e+03, new int[] { 11, 1, 1 });
            p.AddCoeff(6.6910424411861400e+03, new int[] { 11, 1, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7a79be0f-9607-4278-a353-4122cb83626d"));
            OrthonormalPolynomials[803] = p;
            p.AddCoeff(-5.7374935037912700e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(1.7212480511373800e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(1.7212480511373800e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-5.1637441534121500e+01, new int[] { 1, 2, 2 });
            p.AddCoeff(1.2431235924881100e+02, new int[] { 3, 0, 0 });
            p.AddCoeff(-3.7293707774643300e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(-3.7293707774643300e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(1.1188112332393000e+03, new int[] { 3, 2, 2 });
            p.AddCoeff(-7.4587415549286600e+02, new int[] { 5, 0, 0 });
            p.AddCoeff(2.2376224664786000e+03, new int[] { 5, 0, 2 });
            p.AddCoeff(2.2376224664786000e+03, new int[] { 5, 2, 0 });
            p.AddCoeff(-6.7128673994357900e+03, new int[] { 5, 2, 2 });
            p.AddCoeff(1.8114086633398200e+03, new int[] { 7, 0, 0 });
            p.AddCoeff(-5.4342259900194500e+03, new int[] { 7, 0, 2 });
            p.AddCoeff(-5.4342259900194500e+03, new int[] { 7, 2, 0 });
            p.AddCoeff(1.6302677970058400e+04, new int[] { 7, 2, 2 });
            p.AddCoeff(-1.9120424779698100e+03, new int[] { 9, 0, 0 });
            p.AddCoeff(5.7361274339094200e+03, new int[] { 9, 0, 2 });
            p.AddCoeff(5.7361274339094200e+03, new int[] { 9, 2, 0 });
            p.AddCoeff(-1.7208382301728200e+04, new int[] { 9, 2, 2 });
            p.AddCoeff(7.3005258249756200e+02, new int[] { 11, 0, 0 });
            p.AddCoeff(-2.1901577474926900e+03, new int[] { 11, 0, 2 });
            p.AddCoeff(-2.1901577474926900e+03, new int[] { 11, 2, 0 });
            p.AddCoeff(6.5704732424780600e+03, new int[] { 11, 2, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("9ca07729-e67e-40ad-9325-0ef53c9fb1f5"));
            OrthonormalPolynomials[804] = p;
            p.AddCoeff(3.1550997936529100e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-5.2584996560881800e+01, new int[] { 1, 3, 1 });
            p.AddCoeff(-6.8360495529146400e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(1.1393415921524400e+03, new int[] { 3, 3, 1 });
            p.AddCoeff(4.1016297317487800e+03, new int[] { 5, 1, 1 });
            p.AddCoeff(-6.8360495529146400e+03, new int[] { 5, 3, 1 });
            p.AddCoeff(-9.9611007771041800e+03, new int[] { 7, 1, 1 });
            p.AddCoeff(1.6601834628507000e+04, new int[] { 7, 3, 1 });
            p.AddCoeff(1.0514495264721100e+04, new int[] { 9, 1, 1 });
            p.AddCoeff(-1.7524158774535100e+04, new int[] { 9, 3, 1 });
            p.AddCoeff(-4.0146254647116900e+03, new int[] { 11, 1, 1 });
            p.AddCoeff(6.6910424411861400e+03, new int[] { 11, 3, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("807f6dd7-a70d-4342-989f-79b81e585d66"));
            OrthonormalPolynomials[805] = p;
            p.AddCoeff(-5.1637441534121500e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(5.1637441534121500e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(-6.0243681789808400e+01, new int[] { 1, 4, 0 });
            p.AddCoeff(1.1188112332393000e+02, new int[] { 3, 0, 0 });
            p.AddCoeff(-1.1188112332393000e+03, new int[] { 3, 2, 0 });
            p.AddCoeff(1.3052797721125200e+03, new int[] { 3, 4, 0 });
            p.AddCoeff(-6.7128673994357900e+02, new int[] { 5, 0, 0 });
            p.AddCoeff(6.7128673994357900e+03, new int[] { 5, 2, 0 });
            p.AddCoeff(-7.8316786326750600e+03, new int[] { 5, 4, 0 });
            p.AddCoeff(1.6302677970058400e+03, new int[] { 7, 0, 0 });
            p.AddCoeff(-1.6302677970058400e+04, new int[] { 7, 2, 0 });
            p.AddCoeff(1.9019790965068100e+04, new int[] { 7, 4, 0 });
            p.AddCoeff(-1.7208382301728200e+03, new int[] { 9, 0, 0 });
            p.AddCoeff(1.7208382301728200e+04, new int[] { 9, 2, 0 });
            p.AddCoeff(-2.0076446018683000e+04, new int[] { 9, 4, 0 });
            p.AddCoeff(6.5704732424780600e+02, new int[] { 11, 0, 0 });
            p.AddCoeff(-6.5704732424780600e+03, new int[] { 11, 2, 0 });
            p.AddCoeff(7.6655521162244200e+03, new int[] { 11, 4, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("7e79a8cb-b010-4b1c-b290-e326a44718ef"));
            OrthonormalPolynomials[806] = p;
            p.AddCoeff(-1.5826224176235000e+00, new int[] { 0, 0, 1 });
            p.AddCoeff(2.6377040293725000e+00, new int[] { 0, 0, 3 });
            p.AddCoeff(1.2344454857463300e+02, new int[] { 2, 0, 1 });
            p.AddCoeff(-2.0574091429105500e+02, new int[] { 2, 0, 3 });
            p.AddCoeff(-1.5430568571829100e+03, new int[] { 4, 0, 1 });
            p.AddCoeff(2.5717614286381800e+03, new int[] { 4, 0, 3 });
            p.AddCoeff(6.9951910858958600e+03, new int[] { 6, 0, 1 });
            p.AddCoeff(-1.1658651809826400e+04, new int[] { 6, 0, 3 });
            p.AddCoeff(-1.4240210424859400e+04, new int[] { 8, 0, 1 });
            p.AddCoeff(2.3733684041432400e+04, new int[] { 8, 0, 3 });
            p.AddCoeff(1.3290863063202200e+04, new int[] { 10, 0, 1 });
            p.AddCoeff(-2.2151438438670200e+04, new int[] { 10, 0, 3 });
            p.AddCoeff(-4.6316644008128600e+03, new int[] { 12, 0, 1 });
            p.AddCoeff(7.7194406680214600e+03, new int[] { 12, 0, 3 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("1491b728-ec61-4eae-86b6-d2e723799a71"));
            OrthonormalPolynomials[807] = p;
            p.AddCoeff(-7.7224066640437800e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(2.3167219992131300e+00, new int[] { 0, 1, 2 });
            p.AddCoeff(6.0234771979541500e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-1.8070431593862500e+02, new int[] { 2, 1, 2 });
            p.AddCoeff(-7.5293464974426600e+02, new int[] { 4, 1, 0 });
            p.AddCoeff(2.2588039492328100e+03, new int[] { 4, 1, 2 });
            p.AddCoeff(3.4133037455073500e+03, new int[] { 6, 1, 0 });
            p.AddCoeff(-1.0239911236522100e+04, new int[] { 6, 1, 2 });
            p.AddCoeff(-6.9485111962113800e+03, new int[] { 8, 1, 0 });
            p.AddCoeff(2.0845533588634200e+04, new int[] { 8, 1, 2 });
            p.AddCoeff(6.4852771164639400e+03, new int[] { 10, 1, 0 });
            p.AddCoeff(-1.9455831349391900e+04, new int[] { 10, 1, 2 });
            p.AddCoeff(-2.2600208133132000e+03, new int[] { 12, 1, 0 });
            p.AddCoeff(6.7800624399395900e+03, new int[] { 12, 1, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("77af773a-c6bf-49c0-9039-64861bffe47a"));
            OrthonormalPolynomials[808] = p;
            p.AddCoeff(-7.7224066640437800e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(2.3167219992131300e+00, new int[] { 0, 2, 1 });
            p.AddCoeff(6.0234771979541500e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-1.8070431593862500e+02, new int[] { 2, 2, 1 });
            p.AddCoeff(-7.5293464974426600e+02, new int[] { 4, 0, 1 });
            p.AddCoeff(2.2588039492328100e+03, new int[] { 4, 2, 1 });
            p.AddCoeff(3.4133037455073500e+03, new int[] { 6, 0, 1 });
            p.AddCoeff(-1.0239911236522100e+04, new int[] { 6, 2, 1 });
            p.AddCoeff(-6.9485111962113800e+03, new int[] { 8, 0, 1 });
            p.AddCoeff(2.0845533588634200e+04, new int[] { 8, 2, 1 });
            p.AddCoeff(6.4852771164639400e+03, new int[] { 10, 0, 1 });
            p.AddCoeff(-1.9455831349391900e+04, new int[] { 10, 2, 1 });
            p.AddCoeff(-2.2600208133132000e+03, new int[] { 12, 0, 1 });
            p.AddCoeff(6.7800624399395900e+03, new int[] { 12, 2, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("af5d008d-014f-45ed-9098-9b210f8aca92"));
            OrthonormalPolynomials[809] = p;
            p.AddCoeff(-1.5826224176235000e+00, new int[] { 0, 1, 0 });
            p.AddCoeff(2.6377040293725000e+00, new int[] { 0, 3, 0 });
            p.AddCoeff(1.2344454857463300e+02, new int[] { 2, 1, 0 });
            p.AddCoeff(-2.0574091429105500e+02, new int[] { 2, 3, 0 });
            p.AddCoeff(-1.5430568571829100e+03, new int[] { 4, 1, 0 });
            p.AddCoeff(2.5717614286381800e+03, new int[] { 4, 3, 0 });
            p.AddCoeff(6.9951910858958600e+03, new int[] { 6, 1, 0 });
            p.AddCoeff(-1.1658651809826400e+04, new int[] { 6, 3, 0 });
            p.AddCoeff(-1.4240210424859400e+04, new int[] { 8, 1, 0 });
            p.AddCoeff(2.3733684041432400e+04, new int[] { 8, 3, 0 });
            p.AddCoeff(1.3290863063202200e+04, new int[] { 10, 1, 0 });
            p.AddCoeff(-2.2151438438670200e+04, new int[] { 10, 3, 0 });
            p.AddCoeff(-4.6316644008128600e+03, new int[] { 12, 1, 0 });
            p.AddCoeff(7.7194406680214600e+03, new int[] { 12, 3, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("0301f2d1-bdcd-495a-876b-67975e6df7e1"));
            OrthonormalPolynomials[810] = p;
            p.AddCoeff(-6.0234771979541500e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(1.8070431593862500e+01, new int[] { 1, 0, 2 });
            p.AddCoeff(1.8070431593862500e+02, new int[] { 3, 0, 0 });
            p.AddCoeff(-5.4211294781587400e+02, new int[] { 3, 0, 2 });
            p.AddCoeff(-1.5359866854783100e+03, new int[] { 5, 0, 0 });
            p.AddCoeff(4.6079600564349200e+03, new int[] { 5, 0, 2 });
            p.AddCoeff(5.5588089569691200e+03, new int[] { 7, 0, 0 });
            p.AddCoeff(-1.6676426870907300e+04, new int[] { 7, 0, 2 });
            p.AddCoeff(-9.7279156746959700e+03, new int[] { 9, 0, 0 });
            p.AddCoeff(2.9183747024087800e+04, new int[] { 9, 0, 2 });
            p.AddCoeff(8.1360749279275500e+03, new int[] { 11, 0, 0 });
            p.AddCoeff(-2.4408224783782500e+04, new int[] { 11, 0, 2 });
            p.AddCoeff(-2.6077163230536900e+03, new int[] { 13, 0, 0 });
            p.AddCoeff(7.8231489691611000e+03, new int[] { 13, 0, 2 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("e20b9145-5882-4a50-bfe5-c8f6ec8848a6"));
            OrthonormalPolynomials[811] = p;
            p.AddCoeff(1.6162685370654500e+01, new int[] { 1, 1, 1 });
            p.AddCoeff(-4.8488056111963600e+02, new int[] { 3, 1, 1 });
            p.AddCoeff(4.1214847695169000e+03, new int[] { 5, 1, 1 });
            p.AddCoeff(-1.4915849642061200e+04, new int[] { 7, 1, 1 });
            p.AddCoeff(2.6102736873607000e+04, new int[] { 9, 1, 1 });
            p.AddCoeff(-2.1831379930653200e+04, new int[] { 11, 1, 1 });
            p.AddCoeff(6.9972371572606500e+03, new int[] { 13, 1, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("969cde8e-463a-4d6f-834b-35cbf0f4c6c9"));
            OrthonormalPolynomials[812] = p;
            p.AddCoeff(-6.0234771979541500e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(1.8070431593862500e+01, new int[] { 1, 2, 0 });
            p.AddCoeff(1.8070431593862500e+02, new int[] { 3, 0, 0 });
            p.AddCoeff(-5.4211294781587400e+02, new int[] { 3, 2, 0 });
            p.AddCoeff(-1.5359866854783100e+03, new int[] { 5, 0, 0 });
            p.AddCoeff(4.6079600564349200e+03, new int[] { 5, 2, 0 });
            p.AddCoeff(5.5588089569691200e+03, new int[] { 7, 0, 0 });
            p.AddCoeff(-1.6676426870907300e+04, new int[] { 7, 2, 0 });
            p.AddCoeff(-9.7279156746959700e+03, new int[] { 9, 0, 0 });
            p.AddCoeff(2.9183747024087800e+04, new int[] { 9, 2, 0 });
            p.AddCoeff(8.1360749279275500e+03, new int[] { 11, 0, 0 });
            p.AddCoeff(-2.4408224783782500e+04, new int[] { 11, 2, 0 });
            p.AddCoeff(-2.6077163230536900e+03, new int[] { 13, 0, 0 });
            p.AddCoeff(7.8231489691611000e+03, new int[] { 13, 2, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("743cc35f-b500-458e-ba9d-18293395edba"));
            OrthonormalPolynomials[813] = p;
            p.AddCoeff(-6.9078352735584400e-01, new int[] { 0, 0, 1 });
            p.AddCoeff(7.2532270372363600e+01, new int[] { 2, 0, 1 });
            p.AddCoeff(-1.2330485963301800e+03, new int[] { 4, 0, 1 });
            p.AddCoeff(7.8093077767578100e+03, new int[] { 6, 0, 1 });
            p.AddCoeff(-2.3427923330273500e+04, new int[] { 8, 0, 1 });
            p.AddCoeff(3.5922815773086000e+04, new int[] { 10, 0, 1 });
            p.AddCoeff(-2.7214254373550000e+04, new int[] { 12, 0, 1 });
            p.AddCoeff(8.0745589899543900e+03, new int[] { 14, 0, 1 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("218c3c92-ad90-45d8-a1d3-d2d88e534b99"));
            OrthonormalPolynomials[814] = p;
            p.AddCoeff(-6.9078352735584400e-01, new int[] { 0, 1, 0 });
            p.AddCoeff(7.2532270372363600e+01, new int[] { 2, 1, 0 });
            p.AddCoeff(-1.2330485963301800e+03, new int[] { 4, 1, 0 });
            p.AddCoeff(7.8093077767578100e+03, new int[] { 6, 1, 0 });
            p.AddCoeff(-2.3427923330273500e+04, new int[] { 8, 1, 0 });
            p.AddCoeff(3.5922815773086000e+04, new int[] { 10, 1, 0 });
            p.AddCoeff(-2.7214254373550000e+04, new int[] { 12, 1, 0 });
            p.AddCoeff(8.0745589899543900e+03, new int[] { 14, 1, 0 });

            //----------------------------------------------------------------------------------------
            p = new Polynomial(new Guid("59f1f3a0-f3d5-48ec-a374-bd9a60bb9112"));
            OrthonormalPolynomials[815] = p;
            p.AddCoeff(-6.1852100426350100e+00, new int[] { 1, 0, 0 });
            p.AddCoeff(2.4534666502452200e+02, new int[] { 3, 0, 0 });
            p.AddCoeff(-2.7969519812795600e+03, new int[] { 5, 0, 0 });
            p.AddCoeff(1.3984759906397800e+04, new int[] { 7, 0, 0 });
            p.AddCoeff(-3.5738830871905400e+04, new int[] { 9, 0, 0 });
            p.AddCoeff(4.8734769370780000e+04, new int[] { 11, 0, 0 });
            p.AddCoeff(-3.3739455718232300e+04, new int[] { 13, 0, 0 });
            p.AddCoeff(9.3185163412260600e+03, new int[] { 15, 0, 0 });

            //----------------------------------------------------------------------------------------
            #endregion POLY_DEF
#pragma warning restore 612
        }

        /// <summary>
        /// transforms some vertices (<paramref name="EdgeVertices"/>) from the local 2D-coordinate system of either
        /// the top, bottom, left, right, front or back edge (see <see cref="Edge"/>) to the local 
        /// coordinate system of the cube;
        /// </summary>
        /// <param name="EdgeIndex">0, 1, 2, 3, 4 or 5; <see cref="Edge"/></param>
        /// <param name="EdgeVertices">input;</param>
        /// <param name="VolumeVertices">output;</param>
        public override void TransformFaceCoordinates(int EdgeIndex, MultidimensionalArray EdgeVertices, MultidimensionalArray VolumeVertices) {
            if (EdgeVertices.Dimension != 2)
                throw new ArgumentException("dimension of EdgeVertices must be 2.", "EdgeVertices");
            if (VolumeVertices.Dimension != 2)
                throw new ArgumentException("dimension of VolumeVertices must be 2.", "VolumeVertices");
            if (VolumeVertices.GetLength(1) != 3)
                throw new ArgumentException("wrong spatial dimension of output", "VolumeVertices");
            if (EdgeVertices.GetLength(1) != 2)
                throw new ArgumentException("wrong spatial dimension of input", "EdgeVertices");
            if (EdgeVertices.GetLength(0) != VolumeVertices.GetLength(0))
                throw new ArgumentException("mismatch in number of vertices between input and output.", "EdgeVertices,VolumeVertices");

            int L = EdgeVertices.GetLength(0);

            switch (EdgeIndex) {
                case (int)Edge.Front:
                    for (int i = 0; i < L; i++) {
                        VolumeVertices[i, 0] = EdgeVertices[i, 0];
                        VolumeVertices[i, 1] = EdgeVertices[i, 1];
                        VolumeVertices[i, 2] = 1.0;
                    }
                    break;

                case (int)Edge.Back:
                    for (int i = 0; i < L; i++) {
                        VolumeVertices[i, 0] = EdgeVertices[i, 0];
                        VolumeVertices[i, 1] = EdgeVertices[i, 1];
                        VolumeVertices[i, 2] = -1.0;
                    }
                    break;


                case (int)Edge.Left:
                    for (int i = 0; i < L; i++) {
                        VolumeVertices[i, 0] = -1.0;
                        VolumeVertices[i, 1] = EdgeVertices[i, 1];
                        VolumeVertices[i, 2] = EdgeVertices[i, 0];
                    }
                    break;

                case (int)Edge.Right:
                    for (int i = 0; i < L; i++) {
                        VolumeVertices[i, 0] = 1.0;
                        VolumeVertices[i, 1] = EdgeVertices[i, 1];
                        VolumeVertices[i, 2] = EdgeVertices[i, 0];
                    }
                    break;

                case (int)Edge.Top:
                    for (int i = 0; i < L; i++) {
                        VolumeVertices[i, 0] = EdgeVertices[i, 0];
                        VolumeVertices[i, 1] = 1.0;
                        VolumeVertices[i, 2] = EdgeVertices[i, 1];
                    }
                    break;

                case (int)Edge.Bottom:
                    for (int i = 0; i < L; i++) {
                        VolumeVertices[i, 0] = EdgeVertices[i, 0];
                        VolumeVertices[i, 1] = -1.0;
                        VolumeVertices[i, 2] = EdgeVertices[i, 1];
                    }
                    break;


                default:
                    throw new ArgumentException("EdgeIndex out of range");
            }
        }

        /// <summary>
        /// partitions this cube into 8 sub-cubes of equal size;
        /// </summary>
        /// <returns></returns>
        public override AffineTrafo[] GetSubdivision() {
            AffineTrafo[] ret = new AffineTrafo[8];

            for (int i = 0; i < 8; i++) {
                ret[i] = new AffineTrafo(3);
                ret[i].Matrix = MultidimensionalArray.Create(3,3); ret[i].Matrix.AccEye(1.0);
                ret[i].Matrix.Scale(0.5);

                ret[i].Affine = Vertices.ExtractSubArrayShallow(i, -1).To1DArray();
                BLAS.dscal(3, 0.5, ret[i].Affine, 1);
            }

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
            if (pt.Length != 3)
                throw new ArgumentException("wrong spatial dimension.", "pt");

            if ((pt[0] < -1.0 - tolerance) || (pt[0] > 1.0 + tolerance)
                || (pt[1] < -1.0 - tolerance) || (pt[1] > 1.0 + tolerance)
                || (pt[2] < -1.0 - tolerance) || (pt[2] > 1.0 + tolerance))
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

            Nodes = MultidimensionalArray.Create(px * px * px, 3);
            Type = new int[Nodes.GetLength(0)];

            var Nodes1D = GenericBlas.Linspace(-1, 1, px);
            int cnt = 0;
            for (int i = 0; i < px; i++) {
                int i_edge = (i == 0 || i == px - 1) ? 1 : 0;

                for (int j = 0; j < px; j++) {
                    int j_edge = (j == 0 || j == px - 1) ? 1 : 0;

                    for (int k = 0; k < px; k++) {
                        int k_edge = (k == 0 || k == px - 1) ? 1 : 0;

                        Nodes[cnt, 0] = Nodes1D[i];
                        Nodes[cnt, 1] = Nodes1D[j];
                        Nodes[cnt, 2] = Nodes1D[k];

                        Type[cnt] = i_edge + j_edge + k_edge;

                        cnt++;
                    }
                }
            }
        }

        /// <summary>
        /// <see cref="RefElement.GetInterpolationNodes_NonLin"/>
        /// </summary>
        override protected void GetInterpolationNodes_NonLin(CellType Type, out NodeSet InterpolationNodes, out PolynomialList InterpolationPolynomials, out int[] NodeType, out int[] EntityIndex) {
            switch (Type) {
                case CellType.Cube_8: {
                        base.SelectNodalPolynomials(2, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                        return;
                    }
                case CellType.Cube_20: {
                        //MultidimensionalArray _InterpolationNodes;
                        //Polynomial[] _InterpolationPolynomials;
                        //int[] _NodeType;
                        //int[] _EntityIndex;

                        base.SelectNodalPolynomials(3, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex, NodeTypeFilter: new int[] { 0, 1 },
                            ModalBasisSelector: delegate(Polynomial p) {

                            for (int i = 0; i < p.Coeff.Length; i++) {
                                int p1 = p.Exponents[i, 0];
                                int p2 = p.Exponents[i, 1];
                                int p3 = p.Exponents[i, 2];

                                if (p1 > 2 || p2 > 2 || p3 > 2)
                                    return false;

                                if (p1 == 2 && p2 == 2)
                                    return false;
                                if (p2 == 2 && p3 == 2)
                                    return false;
                                if (p1 == 2 && p3 == 2)
                                    return false;
                            }

                            return true;
                        });
                        Debug.Assert(NodeType.Length == 20);

                        //// should be 27 Nodes, so we have to drop seven:
                        //// it will be the volume node in the center of the cell and all face nodes

                        //int _K = _InterpolationNodes.GetLength(0);
                        //Debug.Assert(_K == 27);
                        //Debug.Assert(_NodeType[0] == 0 && _NodeType[1] == 1 && _NodeType[2] == 1 && _NodeType[3] == 1 && _NodeType[4] == 1 && _NodeType[5] == 1 && _NodeType[6] == 1, "first 7 nodes should be the volume node and 6 face nodes");
                        //int D = _InterpolationNodes.GetLength(1);
                        //Debug.Assert(D == 3, "spatial dimension is expected to be 3.");


                        //int K = 20;
                        //int offset = _K - K;
                        //InterpolationNodes = MultidimensionalArray.Create(K, D);
                        //InterpolationNodes.Set(_InterpolationNodes.ExtractSubArrayShallow(new int[] { offset, 0 }, new int[] { _K - 1, D - 1 }));
                        //InterpolationPolynomials = new Polynomial[K];
                        //Array.Copy(_InterpolationPolynomials, offset, InterpolationPolynomials, 0, K);
                        //NodeType = new int[K];
                        //Array.Copy(_NodeType, offset, NodeType, 0, K);

                        //EntityIndex = new int[K];
                        //Array.Copy(_EntityIndex, offset, EntityIndex, 0, K);
                        return;
                    }
                case CellType.Cube_27: {
                        base.SelectNodalPolynomials(3, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                        return;
                    }
                case CellType.Cube_64: {
                        base.SelectNodalPolynomials(4, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                        return;
                    }
                case CellType.Cube_125: {
                        base.SelectNodalPolynomials(5, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                        return;
                    }
                case CellType.Cube_216: {
                        base.SelectNodalPolynomials(6, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                        return;
                    }
                //case CellType.Cube_343: {
                //        base.SelectNodalPolynomials(7, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                //        return;
                //    }
                //case CellType.Cube_512: {
                //        base.SelectNodalPolynomials(8, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                //        return;
                //    }
                //case CellType.Cube_729: {
                //        base.SelectNodalPolynomials(9, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                //        return;
                //    }
                //case CellType.Cube_1000: {
                //        base.SelectNodalPolynomials(10, out InterpolationNodes, out InterpolationPolynomials, out NodeType, out EntityIndex);
                //        return;
                //    }
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
                SwitchNode(ref permutationArray[3], ref permutationArray[4]);
                SwitchNode(ref permutationArray[5], ref permutationArray[6]);
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
            ForeignName = "Hexagon";
            ForeignTypeConstant = 0;
            if (conv == ExchangeFormats.Gmsh) {
                if (Type == CellType.Cube_Linear) {
                    ForeignTypeConstant = 5;
                } else if (Type == CellType.Cube_27) {
                    ForeignTypeConstant = 12;
                } else if (Type == CellType.Cube_20) {
                    ForeignTypeConstant = 27;
                } else if (Type == CellType.Cube_64) {
                    ForeignTypeConstant = 92;
                } else if (Type == CellType.Cube_125) {
                    ForeignTypeConstant = 93;
                } else {
                    throw new NotSupportedException("Wrong minor cell type");
                }
            } else if (conv == ExchangeFormats.CGNS) {
                if (Type == 0) {
                    ForeignTypeConstant = 17;
                } else if (Type == CellType.Cube_20) {
                    ForeignTypeConstant = 18;
                } else if (Type == CellType.Cube_27) {
                    ForeignTypeConstant = 19;
                } else {
                    throw new NotSupportedException("Wrong minor cell type");
                }
            } else if (conv == ExchangeFormats.GambitNeutral) {
                ForeignName = "Brick";
                if (Type == CellType.Cube_Linear) {
                    ForeignTypeConstant = 4;
                } else {
                    throw new NotSupportedException("Wrong minor cell type");
                }
            } else {
                throw new NotSupportedException("Wrong foreign convention type");
            }
        }
    }
}
