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
using System.Diagnostics;
using ilPSP.LinSolvers;
using ilPSP.Tracing;

namespace NSE_SIMPLE {

    /// <summary>
    /// Help class for matrix assemblies
    /// </summary>
    public abstract class SIMPLEMatrixAssembly {

        bool m_IsConstant;
        bool m_OnlyAffine;

        /// <summary>
        /// Ctor
        /// </summary>
        /// <param name="IsConstant">
        /// If true, this matrix assembly is constant.
        /// </param>
        /// <param name="OnlyAffine">
        /// If true, this matrix assembly has got only an affine part.
        /// </param>
        /// <param name="MaxUsePerIterMatrix">
        /// For assemblies which are not constant the maximum number
        /// the matrix can be called via <see cref="AssemblyMatrix"/> in one SIMPLE iteration.
        /// </param>
        /// <param name="MaxUsePerIterAffine">
        /// For assemblies which are not constant the maximum number
        /// the affine part can be called via <see cref="AssemblyAffine"/> in one SIMPLE iteration.
        /// </param>
        protected SIMPLEMatrixAssembly(bool IsConstant, bool OnlyAffine, int MaxUsePerIterMatrix = 1, int MaxUsePerIterAffine = 1) {
            m_IsConstant = IsConstant;
            m_OnlyAffine = OnlyAffine;

            m_MaxUsePerIterMatrix = MaxUsePerIterMatrix;
            m_MaxUsePerIterAffine = MaxUsePerIterAffine;
        }

        bool IsInitialized = false;

        /// <summary>
        /// This method has to be called by the constructor of any derived matrix assembly.
        /// </summary>
        protected void Initialize() {
            using (new FuncTrace()) {
                if (IsInitialized) {
                    throw new ApplicationException("This method can only be called once during the lifetime of this object.");
                } else {
                    m_Affine = ComputeAffine();

                    if (!m_OnlyAffine) {
                        m_Matrix = ComputeMatrix();

                        if (m_Affine != null)
                            Debug.Assert(m_Matrix.RowPartitioning.LocalLength == m_Affine.Length, "Mismatch of dimensions");
                    }

                    IsInitialized = true;
                }
            }
        }

        /// <summary>
        /// Implement this method to define the matrix of this assembly.
        /// </summary>
        /// <returns></returns>
        protected abstract MsrMatrix ComputeMatrix();

        /// <summary>
        /// Implement this method to define the affine part of this assembly.
        /// </summary>
        /// <returns></returns>
        protected abstract double[] ComputeAffine();

        // UseCount for matrix of assemblies which are not constant
        int m_UseCntMatrix = 0;
        int m_MaxUsePerIterMatrix = 1;

        MsrMatrix m_Matrix;

        /// <summary>
        /// The matrix of this assembly.
        /// </summary>
        public MsrMatrix AssemblyMatrix {
            get {

                if (m_OnlyAffine)
                    throw new ApplicationException("This matrix assembly has got only an affine part.");

                if (IsInitialized) {
                    if (!m_IsConstant)
                        m_UseCntMatrix++;

                    if (m_IsConstant || (m_UseCntMatrix <= m_MaxUsePerIterMatrix))
                        return m_Matrix;
                    else
                        throw new ApplicationException("Non-constant AssemblyMatrix has not been updated since last usage!");
                } else
                    throw new ApplicationException("Matrix assembly has not been initialized!");
            }
        }

        // UseCount for affine part of assemblies which are not constant
        int m_UseCntAffine = 0;
        int m_MaxUsePerIterAffine = 1;

        double[] m_Affine;

        /// <summary>
        /// The affine part of this assembly.
        /// </summary>
        public double[] AssemblyAffine {
            get {
                if (IsInitialized) {
                    if (!m_IsConstant)
                        m_UseCntAffine++;

                    if (m_IsConstant || (m_UseCntAffine <= m_MaxUsePerIterAffine))
                        return m_Affine;
                    else
                        throw new ApplicationException("Non-constant AssemblyAffine has not been updated since last usage!");
                } else
                    throw new ApplicationException("Matrix assembly has not been initialized!");
            }
        }

        /// <summary>
        /// True, if this matrix assembly is constant.
        /// </summary>
        public bool IsConstant {
            get {
                return m_IsConstant;
            }
            set {
                m_IsConstant = value;
            }
        }

        /// <summary>
        /// Number of entries which are stored
        /// by this processor.
        /// </summary>
        public int LocalLength {
            get {
                return m_Affine.Length;
            }
        }

        /// <summary>
        /// <see cref="Partition.i0"/>
        /// </summary>
        public int i0 {
            get {
                return m_Matrix.RowPartitioning.i0;
            }
        }

        /// <summary>
        /// Update this assembly.
        /// </summary>
        public void Update() {
            using (new FuncTrace()) {
                if (!m_IsConstant) {
                    if (!m_OnlyAffine) {
                        m_Matrix = ComputeMatrix();
                        m_UseCntMatrix = 0;
                    }
                    m_Affine = ComputeAffine();
                    m_UseCntAffine = 0;
                } else {
                    throw new ApplicationException("This matrix assembly is constant and hence cannot be updated.");
                }
            }
        }
    }
}
