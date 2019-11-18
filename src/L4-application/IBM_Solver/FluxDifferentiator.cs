using BoSSS.Foundation;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace IBM_Solver {
    class VolumeFormDifferentiator : IVolumeForm {

        public VolumeFormDifferentiator(IVolumeForm vf) {
            m_VolForm = vf;
            m_eps = Math.Sqrt(BLAS.MachineEps);

        }

        double m_eps;

        IVolumeForm m_VolForm;

        public TermActivationFlags VolTerms => m_VolForm.VolTerms;

        public IList<string> ArgumentOrdering => m_VolForm.ArgumentOrdering;
        
        string[] m_ParameterOrdering;

        public IList<string> ParameterOrdering {
            get {
                return m_ParameterOrdering;
            }
        }

        public bool IgnoreVectorizedImplementation => false;

        public double VolumeForm(ref CommonParamsVol cpv, double[] U, double[,] GradU, double V, double[] GradV) {
            double ret = 0.0;
            int GAMMA = m_VolForm.ArgumentOrdering.Count;


            double eps = m_eps;
            double[] Utmp = new double[GAMMA];
            Utmp.SetSubVector(cpv.Parameters, 0, GAMMA);
            double[] Udelta = Utmp.CloneAs();

            int NoOfParams = m_VolForm.ParameterOrdering != null ? m_VolForm.ParameterOrdering.Count : 0;
            double[] OrgParams = new double[NoOfParams];
            OrgParams.SetSubVector(cpv.Parameters, GAMMA, NoOfParams);
            CommonParamsVol clonedParams = cpv;
            Debug.Assert(object.ReferenceEquals(cpv, clonedParams) == false);
            clonedParams.Parameters = OrgParams;


            //double[] dU = new double[GAMMA];
            for(int iVar = 0; iVar < GAMMA; iVar++) { // loop over trial variables
                if (U[iVar] != 0.0) {

                    double delta = Math.Abs(Utmp[iVar]) * eps;
                    Udelta[iVar] += delta;

                    double f0 = m_VolForm.VolumeForm(ref clonedParams, Utmp, null, V, GradV);
                    double f1 = m_VolForm.VolumeForm(ref clonedParams, Udelta, null, V, GradV);

                    Udelta[iVar] = Utmp[iVar];

                    double dU_iVar = (f1 - f0) / delta;
                    ret += dU_iVar * U[iVar];
                }
            }
            //ret += GenericBlas.InnerProd(dU, U);




            return ret;
        }
    }
}
