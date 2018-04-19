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
using System.Collections;
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using System.Reflection;

namespace BoSSS.Solution {

    /// <summary>
    /// Base class for all types of histories,
    /// i.e. scalar histories (constant in space)
    /// and scalar- and vector-fields histories (not constant in space).
    /// </summary>
    /// <typeparam name="T"></typeparam>
    abstract public class HistoryBase<T> {

        /// <summary>
        /// ctor
        /// </summary>
        internal HistoryBase(T curr) {
            if (!(curr.GetType().IsValueType) && (curr == null))
                throw new ArgumentNullException();

            m_Current = curr;
        }
        

        /// <summary>
        /// <see cref="PushCount"/>
        /// </summary>
        protected int m_PushCount = 0;

        /// <summary>
        /// counts every call to <see cref="Push"/>
        /// </summary>
        public int PushCount {
            get {
                return m_PushCount;
            }
        }

        /// <summary>
        /// either a double for time-dependent scalars, which are constant in space (e.g. ambient pressure in Low-Mach solver)
        /// or a <see cref="BoSSS.Foundation.DGField"/> or a <see cref="BoSSS.Foundation.VectorField{T}"/>.
        /// </summary>
        protected T m_Current;

        /// <summary>
        /// see <see cref="HistoryLength"/>
        /// </summary>
        protected int m_HistoryLength = 0;

        /// <summary>
        /// the number of timesteps that is stored <em>in addition</em> to the current one
        /// </summary>
        public int HistoryLength {
            get { return m_HistoryLength; }
        }

        /// <summary>
        /// increases the number of previous timesteps 
        /// that should be stored in addition to the current one (see <see cref="m_Current"/> resp. <see cref="VectorFieldHistory{T}.Current"/> resp. <see cref="ScalarFieldHistory{T}.Current"/>);
        /// If the allocated length is already greater than <paramref name="val"/>,
        /// this method has no effect.
        /// </summary>
        /// <param name="val">
        /// new history length
        /// </param>
        public void IncreaseHistoryLength(int val) {
            m_HistoryLength = Math.Max(val, m_HistoryLength);
        }

        /// <summary>
        /// returns the minimum of the two numbers:
        /// <list type="bullet">
        ///   <item>number of <see cref="Push"/> - operations carried out on this object</item>
        ///   <item><see cref="m_HistoryLength"/></item>
        /// </list>
        /// If this value is, e.g. 2, this means that timesteps 1, 0, -1 are available;
        /// </summary>
        /// <returns></returns>
        public int GetPopulatedLength() {
            return History.Count;
        }

        /// <summary>
        /// container for previous values: entries are either of type <see cref="DGField"/> or <see cref="VectorField{t}"/>
        /// </summary>
        protected ArrayList History = new ArrayList();

        /// <summary>
        /// pushes a new variable set onto the top of the stack
        /// </summary>
        /// <remarks>
        /// An example of the Push-operation, for <see cref="HistoryLength"/>==3:
        /// <code>
        /// 
        ///    Timestep Index ->  1  0 -1 -2 
        ///   -----------------------------------
        ///    Before Push(..): | a  b  c  d 
        ///                     | 
        ///    After Push(..):  | a  a  b  c
        /// </code>
        /// Because of this behavior, any difference norm (sometimes incorrectly called 'residual')
        /// should be computed prior to the <see cref="Push"/>-operation.
        /// </remarks>
        abstract public void Push();
    }

    /// <summary>
    /// Common functionality of scalar- and vector field history.
    /// </summary>
    abstract public class HistoryBaseReferenceType : HistoryBase<object> {

        /// <summary>
        /// <see cref="HistoryBase{T}"/>
        /// </summary>
        /// <param name="curr"></param>
        public HistoryBaseReferenceType(object curr)
            : base(curr) {
        }

        /// <summary>
        /// copying
        /// </summary>
        abstract protected void CopyObj(object src, object dest);

        /// <summary>
        /// cloning
        /// </summary>
        abstract protected object CloneObj(object src);

        /// <summary>
        /// <see cref="HistoryBase{T}.Push"/>
        /// </summary>
        public override void Push() {
            if (History.Count >= m_HistoryLength) {
                if (History.Count > 0) {
                    object recycled = History[History.Count - 1];
                    History.RemoveAt(History.Count - 1);
                    History.Insert(0, recycled);

                    CopyObj(m_Current, recycled);
                }
            } else {
                object neu = CloneObj(m_Current);
                History.Insert(0, neu);
            }
            m_PushCount++;
        }
    }

    /// <summary>
    /// History of Scalar Fields
    /// </summary>
    public class ScalarFieldHistory<T> : HistoryBaseReferenceType where T : DGField {

        /// <summary>
        /// constructor
        /// </summary>
        /// <param name="f">
        /// Storage for the current timestep, i.e. <see cref="Current"/>;
        /// This object will be cloned to populate the history.
        /// </param>
        public ScalarFieldHistory(T f)
            : base(f) {
        }

        /// <summary>
        /// the current value
        /// </summary>
        public T Current {
            get {
                return (T)m_Current;
            }
        }

        /// <summary>
        /// gets access to previous (and actual) timesteps
        /// </summary>
        /// <param name="i">
        /// either 1 (equal to <see cref="Current"/>) or a non-positive value:
        /// 0 is the initial value of the current timestep,
        /// -1 is the previous timestep, and so on...
        /// </param>
        /// <returns></returns>
        public T this[int i] {
            get {
                if (i > 1)
                    throw new IndexOutOfRangeException();

                if (i == 1)
                    return (T)m_Current;
                else
                    return (T)History[-i];
            }
        }

        /// <summary> hot stuff </summary>
        protected override void CopyObj(object src, object dest) {
            ((T)dest).CopyFrom((T)src);
        }

        /// <summary> hot stuff </summary>
        protected override object CloneObj(object src) {
            return ((T)src).Clone();
        }
    }

    /// <summary>
    /// History of Vector Fields
    /// </summary>
    public class VectorFieldHistory<T> : HistoryBaseReferenceType where T : DGField {

        /// <summary>
        /// constructor
        /// </summary>
        /// <param name="f">
        /// Storage for the current timestep, i.e. <see cref="Current"/>;
        /// This object will be cloned to populate the history.
        /// </param>
        public VectorFieldHistory(VectorField<T> f)
            : base(f) {
        }

        /// <summary>
        /// the current value (equal to this[1])
        /// </summary>
        public VectorField<T> Current {
            get {
                return (VectorField<T>)m_Current;
            }
        }

        /// <summary>
        /// gets access to previous (and actual) timesteps
        /// </summary>
        /// <param name="i">
        /// either 1 (equal to <see cref="Current"/>) or a non-positive value:
        /// 0 is the initial value of the current timestep,
        /// -1 is the previous timestep, and so on...
        /// </param>
        /// <returns></returns>
        public VectorField<T> this[int i] {
            get {
                if (i > 1)
                    throw new IndexOutOfRangeException();

                if (i == 1)
                    return (VectorField<T>)m_Current;
                else
                    return (VectorField<T>)History[-i];
            }
        }

        /// <summary> hot stuff </summary>
        protected override void CopyObj(object src, object dest) {
            ((VectorField<T>)dest).CopyFrom((VectorField<T>)src);
        }

        /// <summary> hot stuff </summary>
        protected override object CloneObj(object src) {
            return ((VectorField<T>)src).Clone();
        }
    }

    /// <summary>
    /// History of scalars, which are constant in space.
    /// </summary>
    public class ScalarHistory : HistoryBase<double> {

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="Scalar">
        /// Storage for the current timestep, i.e. <see cref="Current"/>;
        /// This object will be cloned to populate the history
        /// </param>
        public ScalarHistory(double Scalar)
            : base(Scalar) {
        }

        /// <summary>
        /// the current value (equal to this[1])
        /// </summary>
        public double Current {
            get {
                return m_Current;
            }
            set {
                m_Current = value;
            }
        }

        /// <summary>
        /// gets access to previous (and actual) timesteps
        /// </summary>
        /// <param name="i">
        /// either 1 (equal to <see cref="Current"/>) or a non-positive value:
        /// 0 is the initial value of the current timestep,
        /// -1 is the previous timestep, and so on...
        /// </param>
        /// <returns></returns>
        public double this[int i] {
            get {
                if (i > 1)
                    throw new IndexOutOfRangeException();

                if (i == 1)
                    return m_Current;
                else
                    return (double)History[-i];
            }
        }

        /// <summary>
        /// <see cref="HistoryBase{T}.Push"/>
        /// </summary>
        public override void Push() {
            if (History.Count >= m_HistoryLength) {
                if (History.Count > 0) {                    
                    History.RemoveAt(History.Count - 1);
                    History.Insert(0, m_Current);                    
                }
            } else {                
                History.Insert(0, m_Current);
            }
            m_PushCount++;
        }
    }
}
