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

using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Platform.LinAlg;
using CNS.Convection;
using CNS.EquationSystem;
using CNS.IBM;
using CNS.MaterialProperty;
using CNS.Solution;
using ilPSP.Utils;
using System;

namespace CNS.Tests.MovingIBMTests {

    class MovingIBMIsentropicVortexTests {


        public static IBMControl MovingFrameIBMIsentropicVortex(string dbPath = null, int dgDegree = 3, int noOfCellsPerDirection = 20, double initialLevelSetPosition = -0.9, double agglomerationThreshold = 0.2) {
            IBMControl c = new IBMControl();

            double advectionVelocity = 1.0;

            c.DbPath = dbPath;
            c.savetodb = dbPath != null;
            c.saveperiod = 1;

            c.ProjectName = "Moving IBM Isentropic vortex";
            c.Tags.Add("Isentropic vortex");
            c.Tags.Add("Moving IBM");

            c.DomainType = DomainTypes.MovingImmersedBoundary;
            c.ActiveOperators = Operators.Convection;
            c.ConvectiveFluxType = ConvectiveFluxTypes.MovingFrameRusanov;
            c.EquationOfState = IdealGas.Air;
            c.MachNumber = 1.0 / Math.Sqrt(c.EquationOfState.HeatCapacityRatio);

            c.TimeSteppingScheme = TimeSteppingSchemes.Explicit;
            c.TimesteppingStrategy = TimesteppingStrategies.MovingFrameFlux;
            c.ExplicitScheme = ExplicitSchemes.RungeKutta;
            c.ExplicitOrder = 1;

            c.AddVariable(Variables.Density, dgDegree);
            c.AddVariable(Variables.Momentum.xComponent, dgDegree);
            c.AddVariable(Variables.Momentum.yComponent, dgDegree);
            c.AddVariable(Variables.Energy, dgDegree);
            c.AddVariable(IBMVariables.LevelSet, 1);

            c.GridFunc = delegate {
                double[] nodes = GenericBlas.Linspace(-10.0, 10.0, noOfCellsPerDirection + 1);
                var grid = Grid2D.Cartesian2DGrid(nodes, nodes, periodicX: true, periodicY: false);
                grid.EdgeTagNames.Add(1, "adiabaticSlipWall");
                grid.DefineEdgeTags(X => 1);
                return grid;
            };

            c.MomentFittingVariant = XQuadFactoryHelper.MomentFittingVariants.Classic;
            c.SurfaceHMF_ProjectNodesToLevelSet = false;
            c.SurfaceHMF_RestrictNodes = true;
            c.SurfaceHMF_UseGaussNodes = false;
            c.VolumeHMF_NodeCountSafetyFactor = 3.0;
            c.VolumeHMF_RestrictNodes = true;
            c.VolumeHMF_UseGaussNodes = false;
            c.LevelSetQuadratureOrder = 10;
            c.LevelSetBoundaryTag = "supersonicInlet";

            c.AgglomerationThreshold = agglomerationThreshold;
            c.SaveAgglomerationPairs = true;

            double gamma = c.EquationOfState.HeatCapacityRatio;
            Func<double[], double, double> x = (X, t) => X[0] - advectionVelocity * t;
            Func<double[], double, double> r = (X, t) => Math.Sqrt(x(X, t) * x(X, t) + X[1] * X[1]);
            Func<double[], double, double> phi = (X, t) => Math.Atan2(X[1], x(X, t));
            Func<double[], double, double> rho = (X, t) => Math.Pow(
                1.0 - 0.5 * (gamma - 1.0) / gamma * Math.Exp(1.0 - r(X, t) * r(X, t)),
                1.0 / (gamma - 1.0));
            Func<double[], double, double> p = (X, t) => Math.Pow(rho(X, t), gamma);
            Func<double[], double, double> uAbs = (X, t) => r(X, t) * Math.Exp(0.5 * (1.0 - r(X, t) * r(X, t)));
            Func<double[], double, double> u = (X, t) => advectionVelocity - Math.Sin(phi(X, t)) * uAbs(X, t);
            Func<double[], double, double> v = (X, t) => Math.Cos(phi(X, t)) * uAbs(X, t);

            c.InitialValues_Evaluators.Add(Variables.Density, X => rho(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.xComponent, X => u(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Velocity.yComponent, X => v(X, 0.0));
            c.InitialValues_Evaluators.Add(Variables.Pressure, X => p(X, 0.0));

            double amplitude = 0.3;
            c.LevelSetFunction = delegate (double[] X, double time) {
                double newLevelSetPosition = initialLevelSetPosition + amplitude * Math.Sin(10.0 * time);
                return X[1] - newLevelSetPosition;
            };
            c.LevelSetVelocity = (X, t) => new Vector3D(
                0.0,
                10.0 * amplitude * Math.Cos(10.0 * t),
                0.0);

            c.AddBoundaryCondition("adiabaticSlipWall");
            c.AddBoundaryCondition("supersonicInlet", Variables.Density, rho);
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity[0], u);
            c.AddBoundaryCondition("supersonicInlet", Variables.Velocity[1], v);
            c.AddBoundaryCondition("supersonicInlet", Variables.Pressure, p);

            c.dtMin = 0.0;
            c.dtMax = 1.0;
            c.CFLFraction = 0.1;
            c.Endtime = 10;
            c.NoOfTimesteps = 1000;

            return c;
        }

    }
}
