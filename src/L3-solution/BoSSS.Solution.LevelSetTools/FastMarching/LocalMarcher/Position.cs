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
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools.FastMarching.LocalMarcher {

    struct Position {

        //public double x;
        //public double y;
        public double[] pos;

        //public Position(double X, double Y) {
        //    x = X;
        //    y = Y;
        //    pos = new double[] { X, Y };
        //}

        public Position(double[] Pos) {
            pos = Pos;
        }

        public Position(Node node) {
            //x = node.X_local;
            //y = node.Y_local;
            pos = node.Pos_local;
        }
    }

    class PositionComparer : IEqualityComparer<Position> {
        public bool Equals(Position A, Position B) {

            Debug.Assert(A.pos.Length == B.pos.Length);

            bool equal = true;
            for (int d = 0; d < A.pos.Length; d++) {
                if (Math.Abs(A.pos[d] - B.pos[d]) >= 1e-8)
                    return false;
            }

            return equal;
        }

        public int GetHashCode(Position A) {
            int temp = (int)(A.pos[0] * 1e10);
            return temp.GetHashCode();
        }
    }
}
