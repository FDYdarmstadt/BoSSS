using BoSSS.Foundation;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.Utils {
    
    /// <summary>
    /// Implements a source term that may be space- and time-, but not solution-dependent.
    /// </summary>
    public class ConstantSource : IVolumeForm, ISupportsJacobianComponent {

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="func">
        /// space- and time-dependent function
        /// </param>
        public ConstantSource(Func<double[], double, double> func) {
            this.m_func = func;
        }

        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="value">
        /// constant value
        /// </param>
        public ConstantSource(double value) : this((X,time) => value) {
        }

        Func<double[], double, double> m_func;

        /// <summary>
        /// not in use, returning null
        /// </summary>
        public virtual IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

        /// <summary>
        /// to be implemented by user 
        /// </summary>
        virtual public IList<string> ArgumentOrdering {
            get {
                return new string[0];
            }
        }


        
       
        /// <summary>
        /// Active terms are <see cref="TermActivationFlags.UxV"/> and
        /// <see cref="TermActivationFlags.V"/>
        /// </summary>
        virtual public TermActivationFlags VolTerms {
            get {
                return  TermActivationFlags.V;
            }
        }

        /// <summary>
        /// translates the source term into <see cref="IVolumeForm.VolumeForm"/>
        /// </summary>
        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            return m_func(cpv.Xglobal, cpv.time) * V;
        }

        /// <summary>
        /// Constant component - returns this object itself.
        /// </summary>
        virtual public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }
    }
}
