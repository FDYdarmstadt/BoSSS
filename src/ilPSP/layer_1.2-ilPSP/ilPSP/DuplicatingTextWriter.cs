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
using System.IO;
using System.Diagnostics;

namespace ilPSP {

    /// <summary>
    /// Duplicates the output of one text writer to multiple text writers;
    /// This can e.g. be used to create a copy
    /// of stdout and stderr streams, or to suppress the stdout on processes with MPI rank unequal to 0;
    /// see <see cref="ilPSP.Environment.StdOut"/>, <see cref="ilPSP.Environment.StdErr"/>.
    /// </summary>
    public class DuplicatingTextWriter : TextWriter, IDisposable {

        /// <summary>
        /// Constructor.
        /// </summary>
        /// <param name="tw0"></param>
        /// <param name="FlushPerio"></param>
        /// <param name="SurpStrem0"></param>
        public DuplicatingTextWriter(TextWriter tw0, uint FlushPerio, bool SurpStrem0) {
            Writer0 = tw0;
            m_flushPeriod = FlushPerio;
            this.surpressStream0 = SurpStrem0;
        }

        TextWriter Writer0;


        HashSet<TextWriter> m__WriterS = new HashSet<TextWriter>(new FuncEqualityComparer<TextWriter>((a, b) => object.ReferenceEquals(a, b)));

        /// <summary>
        /// the text writers which are used for duplication of the console out
        /// </summary>
        /// <remarks>
        /// equality is defined by <see cref="Object.ReferenceEquals"/>
        /// </remarks>
        public ISet<TextWriter> WriterS {
            get {
                return m__WriterS;
            }
        }

  



        ///// <summary>
        ///// Adds a <see cref="StreamWriter"/>, writing to <paramref name="s"/>, to <see cref="WriterS"/>.
        ///// When <paramref name="s"/> is closed, the writer is removed.
        ///// </summary>
        //public void AddStream(Stream s) {
        //    WriterS.Add(new StreamWriter(s));
        //}
        
        
        uint m_flushPeriod;

        /// <summary>
        /// if true, the output to the 0-th stream resp. text-writer is suppressed
        /// </summary>
        public bool surpressStream0 {
            get;
            set;
        }

        private int cnt = 0;

        private void Flush1(int l) {
            if(!surpressStream0) {
                try {
                    Writer0.Flush();
                } catch(ObjectDisposedException) {

                }

            }

            cnt += l;
            if(cnt >= m_flushPeriod) {
                ExecuteOnWriter(tw => tw.Flush());
                cnt = 0;
            }
        }

        private void ExecuteOnWriter(Action<TextWriter> a) {
            try {
                a(Writer0);
            } catch(ObjectDisposedException) {

            }


            List<TextWriter> WritersToRemove = null;
            foreach(TextWriter tw in m__WriterS) {
                try {
                    a(tw);
                } catch(ObjectDisposedException) {
                    if(WritersToRemove == null)
                        WritersToRemove = new List<TextWriter>();
                    WritersToRemove.Add(tw);
                }
            }
            if(WritersToRemove != null) {
                foreach(TextWriter tw in WritersToRemove)
                    m__WriterS.Remove(tw);
            }
        }

        /// <summary>
        /// Writes a string on both writers, followed by a new line.
        /// </summary>
        /// <param name="s"></param>
        public override void WriteLine(string s) {
            if (s == null)
                s = "";
            ExecuteOnWriter(Writer1 => Writer1.WriteLine(s));
            Flush1(s.Length);
        }

        /// <summary>
        /// Writes a strong both writers.
        /// </summary>
        /// <param name="s"></param>
        public override void Write(string s) {
            if (s == null)
                s = "";
            ExecuteOnWriter(Writer1 => Writer1.Write(s));
            Flush1(s.Length);
        }

        /// <summary>
        /// Writes a double on both writers, followed by new line
        /// </summary>
        /// <param name="value"></param>
        public override void WriteLine(double value) {
            ExecuteOnWriter(Writer1 => Writer1.WriteLine(value));
            Flush1(10);
        }

        /// <summary>
        /// Writes an int on both writers, followed by new line
        /// </summary>
        /// <param name="value"></param>
        public override void WriteLine(int value) {
            ExecuteOnWriter(Writer1 => Writer1.WriteLine(value));
            Flush1(10);
        }

        /// <summary>
        /// Writes a bool on both writers, followed by new line
        /// </summary>
        /// <param name="value"></param>
        public override void WriteLine(bool value) {
            ExecuteOnWriter(Writer1 => Writer1.WriteLine(value));
            Flush1(10);
        }

        /// <summary>
        /// Writes a long on both writers, followed by new line
        /// </summary>
        /// <param name="value"></param>
        public override void WriteLine(long value) {
            ExecuteOnWriter(Writer1 => Writer1.WriteLine(value));
            Flush1(10);
        }

        /// <summary>
        /// Writes a float on both writers, followed by new line
        /// </summary>
        /// <param name="value"></param>
        public override void WriteLine(float value) {
            ExecuteOnWriter(Writer1 => Writer1.WriteLine(value));
            Flush1(10);
        }

        /// <summary>
        /// Writes a double on both writers
        /// </summary>
        /// <param name="value"></param>
        public override void Write(double value) {
            ExecuteOnWriter(Writer1 => Writer1.Write(value));
            Flush1(10);
        }

        /// <summary>
        /// Writes an int on both writers
        /// </summary>
        /// <param name="value"></param>
        public override void Write(int value) {
            ExecuteOnWriter(Writer1 => Writer1.Write(value));
            Flush1(10);
        }

        /// <summary>
        /// Writes a bool on both writers
        /// </summary>
        /// <param name="value"></param>
        public override void Write(bool value) {
            ExecuteOnWriter(Writer1 => Writer1.Write(value));
            Flush1(10);
        }

        /// <summary>
        /// Writes a long on both writers
        /// </summary>
        /// <param name="value"></param>
        public override void Write(long value) {
            ExecuteOnWriter(Writer1 => Writer1.Write(value));
            Flush1(10);
        }

        /// <summary>
        /// Writes a float on both writers
        /// </summary>
        /// <param name="value"></param>
        public override void Write(float value) {
            ExecuteOnWriter(Writer1 => Writer1.Write(value));
            Flush1(10);
        }

        /// <summary>
        /// Writes a char on both writers
        /// </summary>
        /// <param name="value"></param>
        public override void Write(char value) {
            ExecuteOnWriter(Writer1 => Writer1.Write(value));
            Flush1(1);
        }

        /// <summary>
        /// The encoding of <see cref="Writer0"/>
        /// </summary>
        public override Encoding Encoding {
            get {
                return Writer0.Encoding;
            }
        }

        /// <summary>
        /// flushing the writer and underlying streams
        /// </summary>
        public override void Flush() {
            this.Flush1((int)this.m_flushPeriod * 2);
        }

        #region IDisposable Members

        private bool disposed = false;

        /// <summary>
        /// tries to close all opened streams
        /// </summary>
        new public void Dispose() {
            if(disposed)
                return;
            try {
                Writer0.Flush();
                Writer0.Close();
            } catch(Exception) {
            }

            try {
                ExecuteOnWriter(delegate (TextWriter Writer1) {
                    Writer1.Flush();
                    Writer1.Close();
                });
            } catch(Exception) {
            }

            disposed = true;
            base.Dispose();
        }

        /// <summary>
        /// Calls <see cref="Dispose"/>
        /// </summary>
        ~DuplicatingTextWriter() {
            this.Dispose();
        }
        #endregion
    }
}
