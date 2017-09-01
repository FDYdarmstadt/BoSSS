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
using System.Globalization;
using System.IO;
using System.Linq;
using BoSSS.Foundation;
using BoSSS.Foundation.IO;
using BoSSS.Solution.Control;
using BoSSS.Solution.Utils;

namespace BoSSS.Solution {

    /// <summary>
    /// Logging of 'Residuals'
    /// </summary>
    public class ResidualLogger {

        /// <summary>
        /// Type of residual - used for formatting output.
        /// </summary>
        enum Residualtype {
            L2Norm,
            L2Error
        }

        int m_MPIrank;

        IDatabaseDriver m_IOMaster;

        Guid SessionGuid;

        Dictionary<string, double> m_Residuals = new Dictionary<string, double>();

        /// <summary>
        /// Ctor.
        /// </summary>
        public ResidualLogger(int MPIrank, IDatabaseDriver IOMaster, Guid _SessionGuid) {
            m_MPIrank = MPIrank;
            m_IOMaster = IOMaster;
            SessionGuid = _SessionGuid;
        }

        /// <summary>
        /// The latest residuals values which are logged, i.e. the last line of the table.
        /// - key: residual name.
        /// - value: the latest value logged under the respective name. 
        /// </summary>
        public Dictionary<string, double> Residuals {
            get {
                return m_Residuals;
            }
        }

        /// <summary>
        /// Create dictionary key and check uniqueness.
        /// </summary>
        /// <param name="ResType"></param>
        /// <param name="FieldIdentification"></param>
        /// <param name="ResidualKey"></param>
        /// <returns></returns>
        string GetDictionaryKey(Residualtype ResType, string FieldIdentification, string ResidualKey) {
            string key;

            if (ResidualKey == null)
                key = ResType.ToString() + " " + FieldIdentification;
            else
                key = ResType.ToString() + " " + ResidualKey;

            CheckDictionaryEntry(key);

            return key;
        }

        /// <summary>
        /// Some checks before writing the residual to the dictionary.
        /// </summary>
        /// <param name="KeyToCheck"></param> 
        void CheckDictionaryEntry(string KeyToCheck) {

            if (!WriteResidualsCalled) {
                // Checks for very first iteration
                if (m_Residuals.ContainsKey(KeyToCheck))
                    throw new ApplicationException("The residual " + KeyToCheck + " already exists. Names of residuals must be unique!");
            } else {
                // Checks for all other iterations
                if (!m_Residuals.ContainsKey(KeyToCheck))
                    throw new ApplicationException("The residual " + KeyToCheck + " is not defined."
                        + " All residuals have to be define prior to the first call of NextIteration() and NextTimestep(Write)!");

                if (m_Residuals[KeyToCheck] != 0.0)
                    throw new ApplicationException("The residual " + KeyToCheck + " has been calculated earlier in this iteration."
                        + " Every residual can be calculated only once per iteration!");
            }
        }

        bool m_WriteResidualsToConsole = true;

        /// <summary>
        /// If true, residuals are written to console - default is true.
        /// </summary>
        public bool WriteResidualsToConsole {
            get {
                return m_WriteResidualsToConsole;
            }
            set {
                if (!WriteResidualsCalled)
                    m_WriteResidualsToConsole = value;
                else
                    throw new NotSupportedException("Value can only be changed prior to the first call NextIteration() and NextTimestep(Write)");
            }
        }

        bool m_WriteResidualsToTextFile = true;

        /// <summary>
        /// If true, residuals are written to text file - default is true.
        /// </summary>
        public bool WriteResidualsToTextFile {
            get {
                return m_WriteResidualsToTextFile;
            }
            set {
                if (!WriteResidualsCalled)
                    m_WriteResidualsToTextFile = value;
                else
                    throw new NotSupportedException("Value can only be changed prior to the first call NextIteration() and NextTimestep(Write)");
            }
        }


        /// <summary>
        /// Name of the text file where the values are logged.
        /// </summary>
        public string TextFileFileName {
            get {
                return m_TextFileFileName;
            }
            set {
                if (!WriteResidualsCalled)
                    m_TextFileFileName = value;
                else
                    throw new NotSupportedException("Value can only be changed prior to the first call NextIteration() and NextTimestep(Write)");
            }
        }

        string m_TextFileFileName = "residual-L2";

        /// <summary>
        /// Difference of current and previous time step of <paramref name="f"/>
        /// in the L2Norm.
        /// </summary>
        /// <typeparam name="F"></typeparam>
        /// <param name="f"></param>
        /// <param name="ResidualKey">
        /// Optional parameter for key of residual. If null, the name of
        /// the field <paramref name="f"/> (cf.
        /// <see cref="DGField.Identification"/> will be used as key.
        /// </param>
        public void ComputeResidual<F>(ScalarFieldHistory<F> f, string ResidualKey = null) where F : DGField {
            ComputeDifference(f[1], f[0], ResidualKey);
        }

        /// <summary>
        /// Difference of current and previous time step for each component
        /// of <paramref name="f"/> in the L2Norm.
        /// </summary>
        /// <typeparam name="F"></typeparam>
        /// <param name="f"></param>
        /// <param name="ResidualKey">
        /// Optional parameter for key of residual. If null, identification of
        /// components of <paramref name="f"/> will be used as key.
        /// </param>
        public void ComputeResidual<F>(VectorFieldHistory<F> f, string[] ResidualKey = null) where F : DGField {

            if ((ResidualKey != null) && (ResidualKey.Length != f.Current.Dim))
                throw new ArgumentException("Number of residual keys must be equal to dimension of VectorField!");

            for (int comp = 0; comp < f.Current.Dim; comp++) {
                if (ResidualKey == null)
                    ComputeDifference(f[1][comp], f[0][comp]);
                else
                    ComputeDifference(f[1][comp], f[0][comp], ResidualKey[comp]);
            }
        }

        /// <summary>
        /// Difference between <paramref name="f1"/> and <paramref name="f0"/>
        /// in the L2Norm.
        /// </summary>
        /// <param name="f1"></param>
        /// <param name="f0"></param>
        /// <param name="ResidualKey">
        /// Optional parameter for key of residual.
        /// If null, the name of <paramref name="f1"/> will be used as key.
        /// </param>
        public void ComputeDifference(DGField f1, DGField f0, string ResidualKey = null) {

            if (f1.Basis.Degree != f0.Basis.Degree)
                throw new ArgumentException("Fields must have equal degrees!");

            string key = GetDictionaryKey(Residualtype.L2Norm, f1.Identification, ResidualKey);

            SinglePhaseField diff = new SinglePhaseField(f1.Basis);
            diff.Acc(1.0, f1);
            diff.Acc(-1.0, f0);
            double res = diff.L2Norm();

            m_Residuals[key] = res;
        }

        /// <summary>
        /// Difference between <paramref name="f1"/> and <paramref name="f0"/> for each component in the L2Norm.
        /// </summary>
        /// <typeparam name="F"></typeparam>
        /// <param name="f1"></param>
        /// <param name="f0"></param>
        /// <param name="ResidualKey">
        /// Optional parameter for key of residual.
        /// If null, identification of components of <paramref name="f1"/> will be used as key.
        /// </param>
        public void ComputeDifference<F>(VectorField<F> f1, VectorField<F> f0, string[] ResidualKey = null) where F : DGField {
            if (f1.Dim != f0.Dim)
                throw new ArgumentException("Fields must have same dimensions!");

            if ((ResidualKey != null) && (ResidualKey.Length != f1.Dim))
                throw new ArgumentException("Number of residual keys must be equal to dimension of VectorField!");

            for (int comp = 0; comp < f1.Dim; comp++) {
                if (ResidualKey == null)
                    ComputeDifference(f1[comp], f0[comp]);
                else
                    ComputeDifference(f1[comp], f0[comp], ResidualKey[comp]);
            }
        }

        /// <summary>
        /// L2Norm of given Field <paramref name="f"/>.
        /// </summary>
        /// <param name="f"></param>
        /// <param name="ResidualKey">
        /// Optional parameter for key of residual.
        /// If null, the name of <paramref name="f"/> will be used as key.
        /// </param>
        public void ComputeL2Norm(DGField f, string ResidualKey = null) {

            string key = GetDictionaryKey(Residualtype.L2Norm, f.Identification, ResidualKey);

            double L2Norm = f.L2Norm();
            m_Residuals[key] = L2Norm;
        }

        /// <summary>
        /// Logs a 'customly computed' value (<paramref name="val"/>).
        /// </summary>
        public void CustomValue(double val, string ResidualKey) {
            string key = GetDictionaryKey(Residualtype.L2Norm, null, ResidualKey);

            m_Residuals[key] = val;
        }


        /// <summary>
        /// L2Norm of each component of given VectorField <paramref name="f"/>.
        /// </summary>
        /// <typeparam name="F"></typeparam>
        /// <param name="f"></param>
        /// <param name="ResidualKey">
        /// Optional parameter for key of residual.
        /// If null, identification of components of <paramref name="f"/> will be used as key.
        /// </param>
        public void ComputeL2Norm<F>(VectorField<F> f, string[] ResidualKey = null) where F : DGField {

            if ((ResidualKey != null) && (ResidualKey.Length != f.Dim))
                throw new ArgumentException("Number of residual keys must be equal to dimension of VectorField!");

            for (int comp = 0; comp < f.Dim; comp++) {
                if (ResidualKey == null)
                    ComputeL2Norm(f[comp]);
                else
                    ComputeL2Norm(f[comp], ResidualKey[comp]);
            }
        }

        /// <summary>
        /// L2Error with respect to <paramref name="AnalyticSolution"/>.
        /// </summary>
        /// <param name="f"></param>
        /// <param name="AnalyticSolution">
        /// If null, nothing will be done.
        /// </param>
        /// <param name="QuadOrder"></param>
        /// <param name="ResidualKey">
        /// Optional parameter for key of residual.
        /// If null, the name of <paramref name="f"/> will be used as key.
        /// </param>
        public void ComputeL2Error(DGField f, Func<double[], double> AnalyticSolution, int QuadOrder, string ResidualKey = null) {
            if (AnalyticSolution != null) {
                string key = GetDictionaryKey(Residualtype.L2Error, f.Identification, ResidualKey);

                double error = f.L2Error(AnalyticSolution.Vectorize(), QuadOrder);
                m_Residuals[key] = error;
            }
        }

        /// <summary>
        /// L2Error for each component of <paramref name="f"/> with respect to <paramref name="AnalyticSolution"/>.
        /// </summary>
        /// <typeparam name="F"></typeparam>
        /// <param name="f"></param>
        /// <param name="AnalyticSolution">
        /// If null, nothing will be done.
        /// </param>
        /// <param name="QuadOrder"></param>
        /// <param name="ResidualKey">
        /// Optional parameter for key of residual.
        /// If null, identification of components of <paramref name="f"/> will be used as key.
        /// </param>
        public void ComputeL2Error<F>(VectorField<F> f, Func<double[], double>[] AnalyticSolution, int QuadOrder, string[] ResidualKey = null) where F : DGField {
            if (AnalyticSolution != null) {
                if (f.Dim != AnalyticSolution.Length)
                    throw new ArgumentException("Mismatch of dimensions!");

                if ((ResidualKey != null) && (f.Dim != ResidualKey.Length))
                    throw new ArgumentException("Number of residual keys must be equal to dimension of VectorField!");

                for (int comp = 0; comp < f.Dim; comp++) {
                    if (ResidualKey == null)
                        ComputeL2Error(f[comp], AnalyticSolution[comp], QuadOrder);
                    else
                        ComputeL2Error(f[comp], AnalyticSolution[comp], QuadOrder, ResidualKey[comp]);
                }
            }
        }

        int m_IterationCounter = 0;
        int m_LineCounter = 0;

        int m_TimeStep = 1;

        /// <summary>
        /// Hacker access to the timestep number; in order to make the timstep correct for restarted sessions.
        /// </summary>
        public int TimeStep {
            get {
                return m_TimeStep;
            }
            set {
                m_TimeStep = value;
            }
        }


        /// <summary>
        /// Hacker access to the iteration counter: some methods like GMRES only compute an 
        /// residual every e.g. 100-th iteration, so the <see cref="NextIteration"/>-method is insufficient in this case.
        /// </summary>
        public int IterationCounter {
            get {
                return m_IterationCounter;
            }
            set {
                m_IterationCounter = value;
            }
        }



        /// <summary>
        /// Write Residuals of this iteration and setup ResidualLogger for next iteration.
        /// All residuals will be set to zero!!!
        /// </summary>
        /// <param name="WriteTimeStepNo">
        /// If true, time step number is written additionally to iteration number.
        /// </param>
        public void NextIteration(bool WriteTimeStepNo) {
            m_IterationCounter++;
            m_LineCounter++;
            WriteResiduals(WriteTimeStepNo);
            ClearResiduals();
        }

        /// <summary>
        /// Setup ResidualLogger for next time step.
        /// All residuals will be set to zero!!!
        /// </summary>
        /// <param name="Write">
        /// If true, residuals of this time step are written.
        /// </param>
        public void NextTimestep(bool Write) {
            m_IterationCounter = 0;
            m_TimeStep++;
            if (Write) {
                WriteResiduals(true);
                m_LineCounter++;
            }
            ClearResiduals();
        }

        /// <summary>
        /// Set all residuals to zero.
        /// </summary>
        void ClearResiduals() {
            foreach (var key in m_Residuals.Keys.ToList())
                m_Residuals[key] = 0.0;
        }

        // False until first call of WriteResiduals(), i.e. 
        // first call of NextIteration() or NextTimestep(Write).
        // Used for writing header of text file and locking 
        // dictionary for residuals.
        bool WriteResidualsCalled = false;

        /// <summary>
        /// Writes all residuals in <see cref="m_Residuals"/>
        /// to text file and console - as configured by
        /// <see cref="WriteResidualsToTextFile"/> and <see cref="WriteResidualsToConsole"/>.
        /// </summary>
        void WriteResiduals(bool WriteTimeStepNo) {

            if (!WriteResidualsCalled) {
                if (m_WriteResidualsToConsole)
                    WriteHeader(WriteTimeStepNo, Console.Out);

                if (m_WriteResidualsToTextFile && m_MPIrank == 0) {
                    if (SessionGuid == Guid.Empty) {
                        m_LogResiduals = new StreamWriter(TextFileFileName + ".txt");
                    } else {
                        m_LogResiduals = m_IOMaster.FsDriver.GetNewLog(TextFileFileName, SessionGuid);
                    }

                    WriteHeader(WriteTimeStepNo, this.m_LogResiduals);
                }
                WriteResidualsCalled = true;
            }

            if (m_WriteResidualsToConsole)
                WriteResiduals(WriteTimeStepNo, Console.Out);
            if (m_WriteResidualsToTextFile)
                WriteResiduals(WriteTimeStepNo, this.m_LogResiduals);
        }

        void WriteHeader(bool WriteTimeStepNo, TextWriter sw) {
            if (m_MPIrank == 0) {
                sw.Write("#Line,");
                if (WriteTimeStepNo)
                    sw.Write("#Time,");
                sw.Write("#Iter");
                foreach (var key in m_Residuals.Keys) {
                    sw.Write("\t " + key);
                }
                sw.WriteLine();
                sw.Flush();
            }
        }

        void WriteResiduals(bool WriteTimeStepNo, TextWriter sw) {
            if (m_MPIrank == 0) {
                sw.Write(m_LineCounter);
                sw.Write(",");
                if (WriteTimeStepNo)
                    sw.Write(m_TimeStep + ",");
                sw.Write(m_IterationCounter);

                foreach (var value in m_Residuals.Values) {
                    sw.Write("\t " + value.ToString("E", NumberFormatInfo.InvariantInfo));
                }
                sw.WriteLine();
                sw.Flush();
            }
        }

        TextWriter m_LogResiduals;

        /// <summary>
        /// Close text writer.
        /// </summary>
        public void Close() {
            if ((m_MPIrank == 0) && (m_WriteResidualsToTextFile) && WriteResidualsCalled) {
                m_LogResiduals.Close();
            }
        }
    }
}
