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
using BoSSS.Solution.NSECommon;
using BoSSS.Foundation.XDG;
using BoSSS.Platform;
using System.Collections;
using ilPSP;
using System.Diagnostics;
using BoSSS.Foundation;

namespace BoSSS.Solution.XNSECommon.Operator.Viscosity {

    public class ViscosityInSpeciesBulk_GradUTerm : BoSSS.Solution.NSECommon.swipViscosity_Term1, ISpeciesFilter {

        public ViscosityInSpeciesBulk_GradUTerm(double penalty, double sw, IncompressibleMultiphaseBoundaryCondMap bcMap, string spcName, SpeciesId spcId, int _d, int _D,
            double _muA, double _muB, double _betaS = 0.0)
            : base(penalty, _d, _D, bcMap, NSECommon.ViscosityOption.ConstantViscosity, constantViscosityValue: double.NegativeInfinity) {

            base.m_alpha = sw;
            this.m_bcMap = bcMap;

            m_spcId = spcId;
            switch (spcName) {
                case "A": currentMu = _muA; complementMu = _muB; break;
                case "B": currentMu = _muB; complementMu = _muA; break;
                default: throw new ArgumentException("Unknown species.");
            }

            double muFactor = Math.Max(currentMu, complementMu) / currentMu;
            base.m_penalty_base = penalty * muFactor;

            int D = base.m_D;
            base.velFunction = D.ForLoop(d => this.m_bcMap.bndFunction[VariableNames.Velocity_d(d) + "#" + spcName]);

            m_beta = _betaS;
        }

        SpeciesId m_spcId;

        public SpeciesId validSpeciesId {
            get { return m_spcId; }
        }

        double currentMu;
        double complementMu;

        IncompressibleMultiphaseBoundaryCondMap m_bcMap;


        protected override double Viscosity(double[] Parameters) {
            return currentMu;
        }


    }

    public class ViscosityInSpeciesBulk_GradUtranspTerm : BoSSS.Solution.NSECommon.swipViscosity_Term2, ISpeciesFilter {

        public ViscosityInSpeciesBulk_GradUtranspTerm(double penalty, double sw, IncompressibleMultiphaseBoundaryCondMap bcMap, string spcName, SpeciesId spcId, int _d, int _D,
            double _muA, double _muB, double _betaS = 0.0)
            : base(penalty, _d, _D, bcMap, NSECommon.ViscosityOption.ConstantViscosity, constantViscosityValue: double.NegativeInfinity) {

            base.m_alpha = sw;
            this.m_bcMap = bcMap;

            m_spcId = spcId;
            switch (spcName) {
                case "A": currentMu = _muA; complementMu = _muB; break;
                case "B": currentMu = _muB; complementMu = _muA; break;
                default: throw new ArgumentException("Unknown species.");
            }

            double muFactor = Math.Max(currentMu, complementMu) / currentMu;
            base.m_penalty_base = penalty * muFactor;

            int D = base.m_D;
            base.velFunction = D.ForLoop(d => this.m_bcMap.bndFunction[VariableNames.Velocity_d(d) + "#" + spcName]);

            m_beta = _betaS;
        }

        SpeciesId m_spcId;

        public SpeciesId validSpeciesId {
            get { return m_spcId; }
        }

        double currentMu;
        double complementMu;

        IncompressibleMultiphaseBoundaryCondMap m_bcMap;


        protected override double Viscosity(double[] Parameters) {
            return currentMu;
        }

    }

    public class DimensionlessViscosityInSpeciesBulk_GradUTerm : BoSSS.Solution.NSECommon.swipViscosity_Term1, ISpeciesFilter {

        public DimensionlessViscosityInSpeciesBulk_GradUTerm(double penalty, double sw, IncompressibleMultiphaseBoundaryCondMap bcMap, string spcName, SpeciesId spcId, int _d, int _D,
            double _reynoldsA, double _reynoldsB)
            : base(penalty, _d, _D, bcMap, NSECommon.ViscosityOption.ConstantViscosityDimensionless, reynolds: 0.0) {

            base.m_alpha = sw;
            this.m_bcMap = bcMap;

            m_spcId = spcId;
            switch (spcName) {
                case "A": base.m_reynolds = _reynoldsA; break;
                case "B": base.m_reynolds = _reynoldsB; break;
                default: throw new ArgumentException("Unknown species.");
            }

            int D = base.m_D;
            base.velFunction = D.ForLoop(d => this.m_bcMap.bndFunction[VariableNames.Velocity_d(d) + "#" + spcName]);
        }

        SpeciesId m_spcId;

        public SpeciesId validSpeciesId {
            get { return m_spcId; }
        }

        IncompressibleMultiphaseBoundaryCondMap m_bcMap;

    }

    public class DimensionlessViscosityInSpeciesBulk_GradUtranspTerm : BoSSS.Solution.NSECommon.swipViscosity_Term2, ISpeciesFilter {

        public DimensionlessViscosityInSpeciesBulk_GradUtranspTerm(double penalty, double sw, IncompressibleMultiphaseBoundaryCondMap bcMap, string spcName, SpeciesId spcId, int _d, int _D,
            double _reynoldsA, double _reynoldsB)
            : base(penalty, _d, _D, bcMap, NSECommon.ViscosityOption.ConstantViscosityDimensionless, reynolds: 0.0) {

            base.m_alpha = sw;
            this.m_bcMap = bcMap;

            m_spcId = spcId;
            switch (spcName) {
                case "A": base.m_reynolds = _reynoldsA; break;
                case "B": base.m_reynolds = _reynoldsB; break;
                default: throw new ArgumentException("Unknown species.");
            }

            int D = base.m_D;
            base.velFunction = D.ForLoop(d => this.m_bcMap.bndFunction[VariableNames.Velocity_d(d) + "#" + spcName]);
        }

        SpeciesId m_spcId;

        public SpeciesId validSpeciesId {
            get { return m_spcId; }
        }

        IncompressibleMultiphaseBoundaryCondMap m_bcMap;

    }
}
