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
using BoSSS.Foundation;
using BoSSS.Solution.NSECommon;
using ilPSP;
using ilPSP.LinSolvers;
using ilPSP.Tracing;
using BoSSS.Foundation.Grid.Classic;

namespace NSE_SIMPLE {

    /// <summary>
    /// Utility class which helps implementing operators for the SIMPLE algorithm.
    /// </summary>
    public abstract class SIMPLEOperator {

        /// <summary>
        /// Na was wohl?
        /// </summary>
        protected GridData GridData {
            get {
                return (GridData)m_RowMapping.GridDat;
            }
        }

        
        UnsetteledCoordinateMapping m_RowMapping;
        UnsetteledCoordinateMapping m_ColMapping;

        protected SinglePhaseField[] m_ParameterFields;        

        bool m_IsConstant;

        bool m_OnlyAffine;

        /// <summary>
        /// If true, this operator has got only an affine part.
        /// </summary>
        public bool OnlyAffine {
            get {
                return m_OnlyAffine;
            }
        }

        bool m_OnlyBoundaryEdges;        

        /// <summary>
        /// Ctor
        /// </summary>        
        /// <param name="RowMapping"></param>
        /// <param name="ColMapping"></param>
        /// <param name="ParameterFields">
        /// Paramter fields for SpatialOperator - can be null.
        /// </param>
        /// <param name="SolverConf">
        /// Solver configuration.
        /// </param>      
        /// <param name="IsConstant">
        /// If true, the operator is constant over time and SIMPLE iterations. 
        /// That is the case for most operatores.
        /// Only operators like e.g. the <see cref="LinearizedConvection"/> are not constant.
        /// </param>
        /// <param name="ArgumentIndex">
        /// Argument index of dependent variable, e.g. spatial component u_i, for SpatialOperator
        /// - this parameter is optional.        
        /// </param>
        /// <param name="SpatialDirection">
        /// Spatial direction index of independent variable, e.g. diff(p, x_i), for SpatialOperator
        /// - this parameter is optional.        
        /// </param>
        /// <param name="OnlyAffine">
        /// If true, only the affine part of the operator is calculated.
        /// </param>
        /// <param name="OnlyBoundaryEdges">
        /// If true, the integration for calculating the affine part is carried out
        /// only on the boundary edges. Can be set to true for most operators.
        /// </param>
        /// <param name="MaxUsePerIterMatrix">
        /// For operators which are not constant the maximum number
        /// the matrix can be called via <see cref="OperatorMatrix"/> in one SIMPLE iteration.
        /// </param>
        /// <param name="MaxUsePerIterAffine">
        /// For operators which are not constant the maximum number
        /// the affine part can be called via <see cref="OperatorAffine"/> in one SIMPLE iteration.
        /// </param>
        protected SIMPLEOperator(UnsetteledCoordinateMapping RowMapping, UnsetteledCoordinateMapping ColMapping, SinglePhaseField[] ParameterFields,
            SolverConfiguration SolverConf, bool IsConstant, int ArgumentIndex = -1, int SpatialDirection = -1,
            bool OnlyAffine = false, bool OnlyBoundaryEdges = true, int MaxUsePerIterMatrix = 1, int MaxUsePerIterAffine = 1) {
                
                m_RowMapping = RowMapping;
                m_ColMapping = ColMapping;

                m_ParameterFields = ParameterFields;

                m_IsConstant = IsConstant;
                m_OnlyAffine = OnlyAffine;
                m_OnlyBoundaryEdges = OnlyBoundaryEdges;

                m_MaxUsePerIterMatrix = MaxUsePerIterMatrix;
                m_MaxUsePerIterAffine = MaxUsePerIterAffine;

                // Get SpatialOperator
                m_SpatialOperator = GetSpatialOperator(SolverConf, ArgumentIndex, SpatialDirection);

                // Initialize and compute matrix of this operator
                if (!m_OnlyAffine) {
                    m_OperatorMatrix = new MsrMatrix(m_RowMapping, m_ColMapping);
                    ComputeMatrix();
                }                

                // Initialize and compute affine part of this operator
                m_OperatorAffine = new double[m_RowMapping.LocalLength];
                ComputeAffine();
        }

        SpatialOperator m_SpatialOperator;        
        
        /// <summary>
        /// Implement this method to define the <see cref="SpatialOperator"/> of any derived <see cref="SIMPLEOperator"/>.
        /// </summary>
        /// <param name="SolverConf">
        /// Solver configuration.
        /// </param>
        /// <param name="ArgumentIndex">
        /// E.g. Spatial component index of dependent variable, e.g. u_i in momentum covenction operator.
        /// </param>
        /// <param name="SpatialDirection">
        /// Spatial direction index of independent variable, e.g. diff(p, x_i).
        /// </param>        
        /// <returns>
        /// SpatialOperator.
        /// </returns>
        protected abstract SpatialOperator GetSpatialOperator(SolverConfiguration SolverConf, int ArgumentIndex, int SpatialDirection);        

        // UseCount for matrix of operators which are not constant
        int m_UseCntMatrix = 0;
        int m_MaxUsePerIterMatrix = 1;

        MsrMatrix m_OperatorMatrix;

        /// <summary>
        /// The matrix of this operator.
        /// </summary>
        public MsrMatrix OperatorMatrix {
            get {
                if (!m_OnlyAffine) { 
                   
                    if (!m_IsConstant)
                        m_UseCntMatrix++;

                    if (m_IsConstant || (m_UseCntMatrix <= m_MaxUsePerIterMatrix))
                        return m_OperatorMatrix;
                    else
                        throw new ApplicationException("Non-constant matrix has not been updated since last usage!");
                } else {
                    throw new ApplicationException("This operator has got only an affine part.");
                }
            }
        }

        // UseCount for affine part of operators which are not constant
        int m_UseCntAffine = 0;
        int m_MaxUsePerIterAffine = 1;

        double[] m_OperatorAffine;        

        /// <summary>
        /// The affine part of this operator.
        /// </summary>
        public double[] OperatorAffine {
            get {
                if (!m_IsConstant)
                    m_UseCntAffine++;

                if (m_IsConstant || (m_UseCntAffine <= m_MaxUsePerIterAffine))
                    return m_OperatorAffine;
                else
                    throw new ApplicationException("Non-constant affine part has not been updated since last usage!");
            }
        }

        /// <summary>
        /// Number of entries which are stored
        /// by this processor.
        /// </summary>
        public int LocalLength {
            get {
                return m_OperatorAffine.Length;
            }
        }

        /// <summary>
        /// Distribution of matrix rows over MPI processors
        /// of <see cref="OperatorMatrix"/>.
        /// </summary>
        public IPartitioning RowPartition {
            get {
                return m_OperatorMatrix.RowPartitioning;
            }
        }

        void ComputeMatrix() {
            using (new FuncTrace()) {
                
                m_OperatorMatrix.Clear();

                m_SpatialOperator.ComputeMatrixEx<MsrMatrix, IList<double>>(m_ColMapping,
                    m_ParameterFields,
                    m_RowMapping,
                    m_OperatorMatrix,
                    null);
            }
        }

        void ComputeAffine() {
            using (new FuncTrace()) {

                Array.Clear(m_OperatorAffine, 0, m_OperatorAffine.Length);

                m_SpatialOperator.ComputeAffine(m_ColMapping,
                    m_ParameterFields,
                    m_RowMapping,
                    m_OperatorAffine,
                    m_OnlyBoundaryEdges);
            }
        }
        
        /// <summary>
        /// Update this operator, i.e. compute matrix and affine part.
        /// </summary>
        /// <param name="UpdateOnlyAffine">
        /// Optional parameter for special operators, where only the affine part has 
        /// to be updated and the matrix is constant (e.g. <see cref="IP2_PressureCorrectionOperator"/>).
        /// </param>
        public void Update(bool UpdateOnlyAffine = false) {
            if (!m_IsConstant) {

                if ((!m_OnlyAffine) && (!UpdateOnlyAffine)) {
                    ComputeMatrix();
                    m_UseCntMatrix = 0;
                }

                ComputeAffine();
                m_UseCntAffine = 0;
            } else {
                throw new ApplicationException("This operator is constant and hence cannot be updated.");
            }
        }
    }
}