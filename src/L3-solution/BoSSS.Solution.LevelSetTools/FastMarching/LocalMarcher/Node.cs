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

namespace BoSSS.Solution.LevelSetTools.FastMarching.LocalMarcher {

    class Node : IMarchingNode {

        public LinkedList<Node> neighbors;
        double x_global;
        double y_global;
        double phi;

        public Node(double X, double Y, double Phi) {
            X_local = X;
            Y_local = Y;
            phi = Phi;
            neighbors = new LinkedList<Node>();
        }

        public void SetGlobalPosition(double X, double Y) {
            x_global = X;
            y_global = Y;
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

                if (neighbor.X_global > X_global) {
                    phi_r = Math.Abs(neighbor.Phi);
                    h_r = neighbor.X_global - X_global;
                }
                if (neighbor.X_global < X_global) {
                    phi_l = Math.Abs(neighbor.Phi);
                    h_l = X_global - neighbor.X_global;
                }
                if (neighbor.Y_global > Y_global) {
                    phi_t = Math.Abs(neighbor.Phi);
                    h_t = neighbor.Y_global - Y_global;
                }
                if (neighbor.Y_global < Y_global) {
                    phi_b = Math.Abs(neighbor.Phi);
                    h_b = Y_global - neighbor.Y_global;
                }
            }
            //Caculate Phi so that abs(grad(phi)) = 1
            phi = Eikonal.approximate(phi_l, phi_r, phi_b, phi_t, h_l, h_r, h_b, h_t, 1);
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

        #endregion
    }
}
