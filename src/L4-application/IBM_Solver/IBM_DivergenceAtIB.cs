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
using BoSSS.Foundation.XDG;
using BoSSS.Solution.NSECommon;
using System.Diagnostics;
using ilPSP.Utils;
using BoSSS.Platform;
using ilPSP;
using BoSSS.Foundation;

namespace BoSSS.Application.IBM_Solver {
    /// <summary>
    /// velocity jump penalty for the divergence operator, on the level set
    /// </summary>
    public class IBM_DivergenceAtIB : ILevelSetForm, ISupportsJacobianComponent {

        LevelSetTracker m_LsTrk;

        public IBM_DivergenceAtIB(int _D, LevelSetTracker lsTrk,
            double vorZeichen, Func<double[], double, ParticleParameters> getParticleParams) {
            this.D = _D;
            this.m_LsTrk = lsTrk;
            this.m_getParticleParams = getParticleParams;
        }

        int D;

        double pRadius;

        /// <summary>
        /// Describes: 0: velX, 1: velY, 2:rotVel,3:particleradius
        /// </summary>
        Func<double[], double, ParticleParameters> m_getParticleParams;

        /// <summary>
        /// the penalty flux
        /// </summary>
        static double DirichletFlux(double UxN_in, double UxN_out) {
            return (UxN_in - UxN_out);
        }

        public double InnerEdgeForm(ref CommonParams inp, double[] U_Neg, double[] U_Pos, double[,] Grad_uA, double[,] Grad_uB, double v_Neg, double v_Pos, double[] Grad_vA, double[] Grad_vB) {

            //double uAxN = GenericBlas.InnerProd(U_Neg, cp.Normal);
            //BoSSS.Foundation.CommonParams inp = cp;
            //var parameters_P = m_getParticleParams(inp.X, inp.time);
            //double[] uLevSet = new double[] { parameters_P.PointVelocity[0], parameters_P.PointVelocity[1] };
            ////double wLevSet = parameters_P[2];
            ////pRadius = parameters_P[3];

            //double[] _uLevSet = new double[D];

            //_uLevSet[0] = uLevSet[0]; //+pRadius*wLevSet*-cp.Normal[1];
            //_uLevSet[1] = uLevSet[1]; //+pRadius * wLevSet * cp.Normal[0];

            //double uBxN = GenericBlas.InnerProd(_uLevSet, cp.Normal);

            //// transform from species B to A: we call this the "A-fictitious" value
            //double uAxN_fict;
            //uAxN_fict = uBxN;

            //double FlxNeg = -DirichletFlux(uAxN, uAxN_fict); // flux on A-side
            ////double FlxPos = 0;

            //return FlxNeg * v_Neg;

            Vector fluidVelocity = new Vector(U_Neg);
            Vector pVelocity = new Vector(new double[D]);
            double[] X;
            switch (D) {
                case 2:
                X = new double[] { inp.X.x, inp.X.y };
                pVelocity = new Vector(new double[]{m_getParticleParams(X, 0.0).PointVelocity[0],
                        m_getParticleParams(X, 0.0).PointVelocity[1] });
                break;
                case 3:
                X = new double[] { inp.X.x, inp.X.y, inp.X.z };
                pVelocity = m_getParticleParams(X, 0.0).PointVelocity;
                break;
                default:
                throw new NotImplementedException();
            }

            return (pVelocity - fluidVelocity) * inp.Normal * v_Neg;
        }



        /*
        public override void PrimalVar_LevelSetFlux(out double FlxNeg, out double FlxPos,
            ref CommonParams cp,
            double[] U_Neg, double[] U_Pos) {
            FlxNeg = 0;
            FlxPos = 0;
        }

        public override void FluxPotential(out double G, double[] U) {
            G = 0;
        }

        public override void Nu(out double NuNeg, out double NuPos, ref CommonParams cp) {
            NuNeg = 1.0;
            NuPos = 1.0;
        }
        */

        public IList<string> ArgumentOrdering {
            get {
                return VariableNames.VelocityVector(this.D);
            }
        }

        public IList<string> ParameterOrdering {
            get {
                return null;
            }
        }

        public int LevelSetIndex {
            get {
                return 0;
            }
        }

        public SpeciesId PositiveSpecies {
            get { return this.m_LsTrk.GetSpeciesId("B"); }
        }

        public SpeciesId NegativeSpecies {
            get { return this.m_LsTrk.GetSpeciesId("A"); }
        }

        public TermActivationFlags LevelSetTerms {
            get {
                return TermActivationFlags.V | TermActivationFlags.UxV;
            }
        }

        /// <summary>
        /// Linear component - returns this object itself.
        /// </summary>
        virtual public IEquationComponent[] GetJacobianComponents(int SpatialDimension) {
            return new IEquationComponent[] { this };
        }

    }
}
