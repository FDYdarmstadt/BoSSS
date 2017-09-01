/*
 *
 * Copyright (c) 2010, Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics)
 *
 * This file is part of the BoSSS software. 
 * The software (source code or binaries compiled from the source code) may not
 * be copied, compiled ore executed, partly or as a whole, without an explicit 
 * written permission from the Fachgebiet fuer Stroemungsdynamik (chair of fluid dynamics), TU Darmstadt.
 *
 */
using System;
using System.Collections.Generic;
using System.Text;

using BoSSS.Foundation;
using BoSSS.Foundation.Quadrature.NonLin;
using BoSSS.Foundation.Quadrature.Linear;
using BoSSS.Platform;
using ilPSP.LinSolvers;
using BoSSS.Foundation.Comm;



namespace BoSSS.Solution.Solvers {
    
    /*
    /// <summary>
    /// 
    /// </summary>
    public class Evaluator {


        /// <summary>
        /// 
        /// </summary>
        /// <param name="context"></param>
        /// <param name="eqs"></param>
        public Evaluator(Context context, ICollection<Equation> eqs)
            : this(context, eqs, new CoordinateMapping(context, Equation.CollectDependentFields(eqs))) {


        }


        /// <summary>
        /// as usual
        /// </summary>
        Context m_Context;


        /// <summary>
        /// 
        /// </summary>
        /// <param name="context">as usual</param>
        /// <param name="rightHandSide"></param>
        /// <param name="Mapping"></param>
        public Evaluator(Context context, ICollection<Equation> rightHandSide, CoordinateMapping Mapping) {

            m_Context = context;

            bool useNonLin = false, useLin = false;
            foreach (Equation r in rightHandSide) {
                if (!r.IsLinear())
                    useNonLin = true;

                if (r.IsPartlyLinear())
                    useLin = true;
            }

            foreach (Field f in Mapping.Fields) {
                m_DependentFields.Add(f);
            }

            m_Mapping = Mapping;



            // create matrix for linear components
            // -----------------------------------


            if (useLin) {
                m_RhsAffineOffset = new double[m_Mapping.NUpdate];
                Partition RhsMatrixPart = new Partition(m_Mapping.NUpdate);
                m_RhsMatrix = new MsrMatrix(RhsMatrixPart, (int)m_Mapping.GlobalCount);

                LECQuadratureEdge mxtbuilder2 = new LECQuadratureEdge(m_Context, rightHandSide, m_RhsMatrix, m_RhsAffineOffset,m_Mapping,m_Mapping);
                mxtbuilder2.Execute();

                LECVolumeQuadrature mtxBuilder = new LECVolumeQuadrature(m_Context, rightHandSide, m_RhsMatrix, m_RhsAffineOffset,m_Mapping,m_Mapping);
                mtxBuilder.Execute();
            }


            // create quartature for nonlinear components
            // ------------------------------------------

            if (useNonLin) {

                ICollection<Field> flds2 = Equation.CollectDependentFields(rightHandSide);
                foreach (Field _f in flds2)
                    if (!m_DependentFields.Contains(_f))
                        m_DependentFields.Add(_f);


                m_RhsNonlinearEdge = new NECQuadratureEdge(m_Context, m_Mapping, rightHandSide);
                m_RhsNonlinearVolume = new NECQuadratureVolume(m_Context, m_Mapping, rightHandSide);
            }


            // create transceiver
            // ------------------
            m_TRX = new Transceiver(m_Context.CommMaster, m_DependentFields);

        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="output"></param>
        /// <param name="comunicate"></param>
        public void Evaluate(double[] output, bool comunicate) {
            int th = m_Context.IOMaster.tracer.EnterFunction("Evaluate");


            if (output.Length < m_Mapping.NUpdate)
                throw new ArgumentOutOfRangeException("output vector to short.");

            if (comunicate) {
                m_TRX.TransceiveStartImReturn();
                m_TRX.TransceiveFinish();
            }
            Array.Clear(output, 0, output.Length);


            if (m_RhsMatrix != null) {

                m_RhsMatrix.Gemv<CoordinateMapping, double[]>(1.0, m_Mapping, 0.0, output);
                BLAS.daxpy(m_Mapping.NUpdate, 1.0, m_RhsAffineOffset, 1, output, 1);

            }

            if (m_RhsNonlinearEdge != null) {
                m_RhsNonlinearEdge.m_Output = output;
                m_RhsNonlinearEdge.Execute();

                m_RhsNonlinearVolume.m_Output = output;
                m_RhsNonlinearVolume.Execute();
            }

            m_Context.IOMaster.tracer.LeaveFunction(th);
        }


        /// <summary>
        /// 
        /// </summary>
        ICollection<Field> m_DependentFields = new List<Field>();



        /// <summary>
        /// <see cref="Mapping"/>
        /// </summary>
        CoordinateMapping m_Mapping;

        /// <summary>
        /// coordinate mapping which is used to order the output of <see cref="Evaluate"/>;
        /// </summary>
        public CoordinateMapping Mapping {
            get { return m_Mapping; }
        }



        /// <summary>
        /// Tranceiver for the fields within <see cref="Mapping"/>
        /// </summary>
        Transceiver m_TRX;

        /// <summary>
        /// if the right-hand-side is present and contains linear components, this is their matrix;
        /// otherwise, this member is null;
        /// </summary>
        MsrMatrix m_RhsMatrix;


        /// <summary>
        /// if the right-hand-side is present and contains linear components, this is their affine offset vector;
        /// otherwise, this member is null;
        /// </summary>
        double[] m_RhsAffineOffset;


        /// <summary>
        /// if the right-hand-side is present and contains nonlinear components, 
        /// this is the corresponding edge quadrature;
        /// otherwise, this member is null;
        /// </summary>
        NECQuadratureEdge m_RhsNonlinearEdge;

        /// <summary>
        /// if the right-hand-side is present and contains nonlinear components, 
        /// this is the corresponding volume quadrature;
        /// otherwise, this member is null;
        /// </summary>
        NECQuadratureVolume m_RhsNonlinearVolume;


    }
     */
}
