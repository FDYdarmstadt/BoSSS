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
using System.Linq;
using BoSSS.Solution.LevelSetTools.FastMarcher;
using System.Text;
using System.Threading.Tasks;
using System.Diagnostics;

namespace BoSSS.Solution.LevelSetTools.FastMarching.LocalMarcher {

    class Node : IMarchingNode {

        public LinkedList<Node> neighbors;
        double x_global;
        double y_global;
        double[] pos_global;
        double phi;

        public Node(double X, double Y, double Phi) {
            X_local = X;
            Y_local = Y;
            Pos_local = new double[] { X, Y };
            phi = Phi;
            neighbors = new LinkedList<Node>();
        }

        public void SetGlobalPosition(double X, double Y) {
            x_global = X;
            y_global = Y;
            pos_global = new double[] { X, Y };
        }

        public double X_local {
            get;
        }

        public double Y_local {
            get;
        }

        public double X_global {
            get {
                return x_global;
            }
        }

        public double Y_global {
            get {
                return y_global;
            }
        }


        public Node(double[] Pos, double Phi) {
            Pos_local = Pos;
            phi = Phi;
            neighbors = new LinkedList<Node>();
        }

        public void SetGlobalPosition(double[] Pos) {
            pos_global = Pos;
        }

        public double[] Pos_local {
            get;
        }

        public double[] Pos_global {
            get {
                return pos_global;
            }
        }


        public double Phi {
            set { 
                phi = value;
            }
            get {
                return phi;
            }
        }

        #region INode

        public IMarchingNode[] Neighbors {
            get {
                return neighbors.ToArray();
            }
        }

        //Calculate tempPhi with Eikonalapproximation 
        public void CalculateValue() {

            int D = pos_global.Length;
            Debug.Assert(D == Pos_local.Length);

            if (D == 2)
                CalculateValue_2D();

            if (D == 3)
                CalculateValue_3D();
            
        }


        public void CalculateValue_2D() {

            double phi_l = double.MaxValue;
            double phi_r = double.MaxValue;
            double phi_b = double.MaxValue;
            double phi_t = double.MaxValue;
            double h_l = double.MaxValue;
            double h_r = double.MaxValue;
            double h_b = double.MaxValue;
            double h_t = double.MaxValue;

            //Exctract u_... and h_... from neighbors
            foreach (Node neighbor in neighbors) {

                if (neighbor.Pos_global[0] > Pos_global[0]) {
                    phi_r = Math.Abs(neighbor.Phi);
                    h_r = neighbor.Pos_global[0] - Pos_global[0];
                }
                if (neighbor.Pos_global[0] < Pos_global[0]) {
                    phi_l = Math.Abs(neighbor.Phi);
                    h_l = Pos_global[0] - neighbor.Pos_global[0];
                }
                if (neighbor.Pos_global[1] > Pos_global[1]) {
                    phi_t = Math.Abs(neighbor.Phi);
                    h_t = neighbor.Pos_global[1] - Pos_global[1];
                }
                if (neighbor.Pos_global[1] < Pos_global[1]) {
                    phi_b = Math.Abs(neighbor.Phi);
                    h_b = Pos_global[1] - neighbor.Pos_global[1];
                }
            }
            //Caculate Phi so that abs(grad(phi)) = 1
            phi = Eikonal.approximate2D(phi_l, phi_r, phi_b, phi_t, h_l, h_r, h_b, h_t, 1);
        }


        public void CalculateValue_3D() {

            double phi_xm = double.MaxValue;    // x-minus (left neighbor)
            double phi_xp = double.MaxValue;    // x-minus (right neighbor)
            double phi_ym = double.MaxValue;    // y-minus (bottom neighbor)
            double phi_yp = double.MaxValue;    // y-minus (top neighbor)
            double phi_zm = double.MaxValue;    // z-minus (rear neighbor)
            double phi_zp = double.MaxValue;    // z-minus (front neighbor)
            double h_xm = double.MaxValue;
            double h_xp = double.MaxValue;
            double h_ym = double.MaxValue;
            double h_yp = double.MaxValue;
            double h_zm = double.MaxValue;
            double h_zp = double.MaxValue;

            //Exctract u_... and h_... from neighbors
            foreach (Node neighbor in neighbors) {

                if (neighbor.Pos_global[0] > Pos_global[0]) {
                    phi_xp = Math.Abs(neighbor.Phi);
                    h_xp = neighbor.Pos_global[0] - Pos_global[0];
                }
                if (neighbor.Pos_global[0] < Pos_global[0]) {
                    phi_xm = Math.Abs(neighbor.Phi);
                    h_xm = Pos_global[0] - neighbor.Pos_global[0];
                }
                if (neighbor.Pos_global[1] > Pos_global[1]) {
                    phi_yp = Math.Abs(neighbor.Phi);
                    h_yp = neighbor.Pos_global[1] - Pos_global[1];
                }
                if (neighbor.Pos_global[1] < Pos_global[1]) {
                    phi_ym = Math.Abs(neighbor.Phi);
                    h_ym = Pos_global[1] - neighbor.Pos_global[1];
                }
                if (neighbor.Pos_global[2] > Pos_global[2]) {
                    phi_zp = Math.Abs(neighbor.Phi);
                    h_zp = neighbor.Pos_global[2] - Pos_global[2];
                }
                if (neighbor.Pos_global[2] < Pos_global[2]) {
                    phi_zm = Math.Abs(neighbor.Phi);
                    h_zm = Pos_global[2] - neighbor.Pos_global[2];
                }
            }
            //Caculate Phi so that abs(grad(phi)) = 1
            phi = Eikonal.approximate3D(phi_xm, phi_xp, phi_ym, phi_yp, phi_zm, phi_zp, h_xm, h_xp, h_ym, h_yp, h_zm, h_zp, 1);
        }


        //QueueID for Heap
        public int QueueID {
            get;
            set;
        }

        public double Value {
            get {
                return phi;
            }
        }

        public void Accept() {
        }

        public MarchingNodeStatus StatusTag {
            get;
            set;
        }

        #endregion
    }
}
