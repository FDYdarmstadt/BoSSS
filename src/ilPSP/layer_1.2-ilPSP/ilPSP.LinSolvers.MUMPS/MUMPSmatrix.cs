/* 
 * Copyright (C) 2010, Florian Kummer, Technische Universitaet Darmstadt, Fachgebiet fuer Stroemungsmechanik
 *
 * Use, modification and distribution is subject to the Boost Software
 * License, Version 1.0. (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *  
 * Authors: Florian Kummer
 * 
 */
using System;
using System.Collections.Generic;
using System.Text;

namespace ilPSP.LinSolvers.MUMPS {
    
    /// <summary>
    /// MUMPS matrix class
    /// </summary>
    class MUMPSmatrix : ISparseMatrix {

        #region ISparseMatrix Members

        public double GetDiagElement(int row) {
            throw new NotImplementedException();
        }

        public void SetDiagElement(int row, double val) {
            throw new NotImplementedException();
        }

        public Partition RowPartition {
            get {
                throw new NotImplementedException();
            }
        }

        public int NoOfRows {
            get {
                throw new NotImplementedException();
            }
        }

        public int NoOfCols {
            get {
                throw new NotImplementedException();
            }
        }

        public void SpMV<VectorType1, VectorType2>(double alpha, VectorType1 a, double beta, VectorType2 acc)
            where VectorType1 : IList<double>
            where VectorType2 : IList<double> {
            throw new NotImplementedException();
        }

        #endregion
    }
}
