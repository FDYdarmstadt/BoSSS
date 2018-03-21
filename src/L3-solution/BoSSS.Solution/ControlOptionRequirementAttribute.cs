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

namespace BoSSS.Solution.Control {

    /// <summary>
    /// Base class for checks on control files attributes
    /// </summary>
    [AttributeUsage(AttributeTargets.Field, AllowMultiple = false, Inherited = false)]
    public abstract class ControlOptionRequirementAttribute : Attribute {

        /// <summary>
        /// Override this method by implementing a check on the <paramref name="obj"/>
        /// which is the current value of the property <paramref name="propertyName"/>
        /// </summary>
        /// <param name="propertyName">name used in error message</param>
        /// <param name="obj"></param>
        /// <returns>
        /// - null, if everything is ok with <paramref name="obj"/>
        /// - An error message, otherwise.
        /// </returns>
        abstract internal string Verify(string propertyName, object obj);
    }

    /// <summary>
    /// For control file entries that must not be null
    /// </summary>
    public class NotNullAttribute : ControlOptionRequirementAttribute {

        internal override string Verify(string propertyName, object obj) {
            if (obj == null) {
                return String.Format(
                    "Control file property '{0}' must not be null", propertyName);
            } else {
                return null;
            }
        }
    }

    /// <summary>
    /// Abstract class for defining a 'BoundAttribute'. Implemented by
    /// <see cref="LowerBoundAttribute"/> and <see cref="UpperBoundAttribute"/>.
    /// </summary>
    abstract public class BoundAttribute : ControlOptionRequirementAttribute {

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="bound">
        /// Value of the bound.
        /// </param>
        public BoundAttribute(double bound) {
            m_Bound = Convert.ToDouble(bound);
        }

        /// <summary>
        /// value of the bound
        /// </summary>
        protected double m_Bound;
    }

    public class InclusiveLowerBoundAttribute : BoundAttribute {

        public InclusiveLowerBoundAttribute(double lowerBound)
            : base(lowerBound) {
        }

        internal override string Verify(string propertyName, object obj) {
            double value = Convert.ToDouble(obj);

            if (value < m_Bound) {
                return "For the property '" + propertyName
                    + "' a value of " + value + " is found in the control-file, but by definition it has got a lower bound (inclusive) of "
                    + m_Bound + ".";
            } else {
                return null;
            }
        }
    }

    /// <summary>
    /// Attribute for lower bounds.
    /// </summary>    
    public class ExclusiveLowerBoundAttribute : BoundAttribute {

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="lowerBound">
        /// Value of the lower bound.
        /// </param>
        public ExclusiveLowerBoundAttribute(double lowerBound)
            : base(lowerBound) {
        }

        internal override string Verify(string propertyName, object obj) {
            double value = Convert.ToDouble(obj);

            if (value <= m_Bound) {
                return "For the property '" + propertyName
                    + "' a value of " + value + " is found in the control-file, but by definition it has got a lower bound (exclusive) of "
                    + m_Bound + ".";
            } else {
                return null;
            }
        }
    }

    public class InclusiveUpperBoundAttribute : BoundAttribute {

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="upperBound">
        /// Value of the upper bound.
        /// </param>
        public InclusiveUpperBoundAttribute(double upperBound)
            : base(upperBound) {
        }

        internal override string Verify(string propertyName, object obj) {
            double val = Convert.ToDouble(obj);

            if (val > m_Bound) {
                return "For the property '" + propertyName
                    + "' a value of " + val + " is found in the control-file, but by definition it has got a upper bound (inclusive) of "
                    + m_Bound + ".";
            } else {
                return null;
            }
        }
    }

    public class ExclusiveUpperBoundAttribute : BoundAttribute {

        public ExclusiveUpperBoundAttribute(double upperBound)
            : base(upperBound) {
        }

        internal override string Verify(string propertyName, object obj) {
            double val = Convert.ToDouble(obj);

            if (val >= m_Bound) {
                return "For the property '" + propertyName
                    + "' a value of " + val + " is found in the control-file, but by definition it has got a upper bound (exclusive) of "
                    + m_Bound + ".";
            } else {
                return null;
            }
        }
    }
}