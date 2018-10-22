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
using System.Runtime.InteropServices;
using ilPSP;
using ilPSP.Utils;

namespace BoSSS.Platform.LinAlg {

    /*

    /// <summary>
    /// A 2D Vector2D, or Stage 1 Tensor in 3D Space;
    /// </summary>
    [StructLayout(LayoutKind.Sequential)]
    public struct Vector3D {

        

        

       

        
       
        


        

        
        

        

        /// <summary>
        /// Computes the transformation for a arbitrary vector from the
        /// standard coordinate system into a coordinate system with the axes
        /// \f$  \vec{e}_1 = \vec{n}\f$ ,
        /// \f$  \vec{e}_2 = (-n_2, n_1, 0)^T\f$  and
        /// \f$  \vec{e}_3 = \vec{e}_1 \times \vec{e}_2\f$ ,
        /// where \f$  \vec{n} = (n_1, n_2, n_3)^T\f$ 
        /// is given by <paramref name="edgeNormal"/>.
        /// </summary>
        /// <param name="edgeNormal">
        /// The normal of an edge
        /// </param>
        /// <returns>
        /// An orthonormal 3x3 matrix that, when applied to a vector,
        /// transforms this vector into the above-mentioned coordinate system.
        /// </returns>
        private MultidimensionalArray GetTransformationToEdgeCoordinates(Vector3D edgeNormal) {
            Vector3D t1 = new Vector3D(-edgeNormal[1], edgeNormal[0], 0.0);
            Vector3D t2 = edgeNormal.CrossProduct(t1);
            MultidimensionalArray trafo = MultidimensionalArray.Create(3, 3);
            for (int i = 0; i < 3; i++) {
                trafo[0, i] = edgeNormal[i];
                trafo[1, i] = t1[i];
                trafo[2, i] = t2[i];
            }

            return trafo;
        }

        /// <summary>
        /// Transforms this vector into a new vector in a coordinate system whose
        /// first axis is aligned with the given <paramref name="edgeNormal"/>.
        /// </summary>
        /// <param name="edgeNormal">
        /// The normal of an edge
        /// </param>
        /// <returns>
        /// This vector in a coordinate system with the axes
        /// \f$  \vec{e}_1 = \vec{n}\f$ ,
        /// \f$  \vec{e}_2 = (-n_2, n_1, 0)^T\f$  and
        /// \f$  \vec{e}_3 = \vec{e}_1 \times \vec{e}_2\f$ ,
        /// where \f$  \vec{n} = (n_1, n_2, n_3)^T\f$ 
        /// is given by <paramref name="edgeNormal"/>
        /// </returns>
        public Vector3D ToEdgeCoordinates(Vector3D edgeNormal) {
            double[] transformedVector = new double[3];
            GetTransformationToEdgeCoordinates(edgeNormal).gemv(
                1.0, (double[])this, 0.0, transformedVector);
            return new Vector3D(transformedVector);
        }

        /// <summary>
        /// Transforms this vector from a coordinate system as defined by
        /// <see cref="ToEdgeCoordinates"/> into a new vector in the standard
        /// coordinate system
        /// </summary>
        /// <param name="edgeNormal">
        /// The normal of an edge
        /// </param>
        /// <returns>
        /// This vector in the standard coordinate system
        /// </returns>
        public Vector3D FromEdgeCoordinates(Vector3D edgeNormal) {
            double[] transformedVector = new double[3];
            GetTransformationToEdgeCoordinates(edgeNormal).Transpose().gemv(
                1.0, (double[])this, 0.0, transformedVector);
            return new Vector3D(transformedVector);
        }
    }

    */
}
