using BoSSS.Foundation;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.Utils
{
    /// <summary>
    /// Simplified implementation of <see cref="IEdgeForm"/> and <see cref="IVolumeForm"/> for scalar equations,
    /// i.e. equations with only one argument (<see cref="ArgumentOrdering"/>).
    /// </summary>
    abstract public class ScalarEquationForm : IEdgeForm, IVolumeForm
    {
        /// <summary>
        /// ctor.
        /// </summary>
        public ScalarEquationForm(string VarName) {
            m_ArgumentOrdering = new string[] { VarName };
        }


        /// <summary>
        /// set to <see cref="TermActivationFlags.AllOn"/>
        /// </summary>
        virtual public TermActivationFlags BoundaryEdgeTerms {
            get {
                return TermActivationFlags.AllOn;
            }
        }

        /// <summary>
        /// set to <see cref="TermActivationFlags.AllOn"/>
        /// </summary>
        public TermActivationFlags InnerEdgeTerms {
            get {
                return TermActivationFlags.AllOn;
            }
        }

        string[] m_ArgumentOrdering;

        /// <summary>
        /// Single variable handed to constructor.
        /// </summary>
        public IList<string> ArgumentOrdering {
            get {
                return m_ArgumentOrdering;
            }
        }

        /// <summary>
        /// set to null
        /// </summary>
        virtual public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

        /// <summary>
        /// set to <see cref="TermActivationFlags.AllOn"/>
        /// </summary>
        public TermActivationFlags VolTerms {
            get {
                return TermActivationFlags.AllOn;
            }
        }
        /// <summary>
        /// redirected to <see cref="BoundaryEdgeForm(ref CommonParamsBnd, double, double[], double, double[])"/>
        /// </summary>
        public double BoundaryEdgeForm(ref CommonParamsBnd inp, double[] _uA, double[,] _Grad_uA, double _vA, double[] _Grad_vA) {
            return BoundaryEdgeForm(ref inp, _uA[0], _Grad_uA.GetRow(0), _vA, _Grad_vA);
        }


        /// <summary>
        ///  redirected to <see cref="InnerEdgeForm(ref CommonParams, double, double, double[], double[], double, double, double[], double[])"/>
        /// </summary>
        public double InnerEdgeForm(ref CommonParams inp, double[] _uIN, double[] _uOUT, double[,] _Grad_uIN, double[,] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT) {
            return InnerEdgeForm(ref inp, _uIN[0], _uOUT[0], _Grad_uIN.GetRow(0), _Grad_uOUT.GetRow(0), _vIN, _vOUT, _Grad_vIN, _Grad_vOUT);
        }

        /// <summary>
        /// redirected to <see cref="VolumeForm(ref CommonParamsVol, double, double[], double, double[])"/>
        /// </summary>
        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            return VolumeForm(ref cpv, U[0], GradU.GetRow(0), V, GradV);
        }

        /// <summary>
        /// override to implement the BoundaryEdgeForm
        /// </summary>
        abstract public double BoundaryEdgeForm(ref CommonParamsBnd inp, double _uA, double[] _Grad_uA, double _vA, double[] _Grad_vA);


        /// <summary>
        /// override to implement the InnerEdgeForm
        /// </summary>
        abstract public double InnerEdgeForm(ref CommonParams inp, double _uIN, double _uOUT, double[] _Grad_uIN, double[] _Grad_uOUT, double _vIN, double _vOUT, double[] _Grad_vIN, double[] _Grad_vOUT);


        /// <summary>
        /// override to implement the volume form
        /// </summary>
        abstract public double VolumeForm(ref CommonParamsVol cpv, double U, double[] GradU, double V, double[] GradV);


    }
}
