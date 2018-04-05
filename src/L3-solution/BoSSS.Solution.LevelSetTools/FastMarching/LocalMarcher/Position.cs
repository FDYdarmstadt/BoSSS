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
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.LevelSetTools.FastMarching.LocalMarcher {

    struct Position {

        public double x;
        public double y;

        public Position(double X, double Y) {
            x = X;
            y = Y;
        }
        public Position(Node node) {
            x = node.X_local;
            y = node.Y_local;
        }
    }

    class PositionComparer : IEqualityComparer<Position> {
        public bool Equals(Position A, Position B) {

            if ((Math.Abs(A.x - B.x)) < 1e-8) {
                if (Math.Abs(A.y - B.y) < 1e-8) {
                    return true;
                }
            }

            return false;
        }

        public int GetHashCode(Position A) {
            int temp = (int)(A.x * 1e10);
            return temp.GetHashCode();
        }
    }
}
