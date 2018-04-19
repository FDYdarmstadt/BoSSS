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

using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using System.Linq;

namespace CNS.Tests.ViscousShockProfile {

    class ViscousShockProfile : Program<CNSControl> {

        //static void Main(string[] args) {
        //    Application<CNSControl>._Main(args, false, null, () => new ViscousShockProfile());
        //}
        
        protected override void SetInitial() {

            //base.SetInitial();
            WorkingSet.ProjectInitialValues(SpeciesMap, base.Control.InitialValues_Evaluators);
            string viscosityLaw = Control.ViscosityLaw.ToString();

            if (true) {
                DGField mpiRank = new SinglePhaseField(new Basis(GridData, 0), "rank");
                m_IOFields.Add(mpiRank);

                for (int j = 0; j < GridData.Cells.NoOfLocalUpdatedCells; j++) {
                    mpiRank.SetMeanValue(j, DatabaseDriver.MyRank);
                }
            }

            // exakte Lösung für 1D-Testfall
            int p = Control.MomentumDegree;
            CellQuadratureScheme scheme = new CellQuadratureScheme(true);
            int order = WorkingSet.Momentum[0].Basis.Degree * 2 + 2;
            int noOfNodesPerCell = scheme.Compile(GridData, order).First().Rule.NoOfNodes;
            int offset = Grid.CellPartitioning.i0 * noOfNodesPerCell;
            long K = Grid.NumberOfCells;

            string pathToData;
            string initial;
            if (Grid.Description.Contains("Unsteady")) {
                pathToData = @"P:\ViscousTerms\NS_exact\data\M4_" + viscosityLaw + "_unsteady\\";
                initial = "1";
            } else {
                pathToData = @"P:\ViscousTerms\NS_exact\data\M4_" + viscosityLaw + "_steady\\";
                initial = "0";
            }

            ScalarFunction func = new ScalarFunction(
                (MultidimensionalArray arrIn, MultidimensionalArray arrOut) => {

                    string filename = "m" + initial + "_" + K.ToString() + "_" + Control.MomentumDegree.ToString() + ".txt";
                    double[] output;
                    using (System.IO.StreamReader sr = new System.IO.StreamReader(pathToData + filename)) {
                        output = sr.ReadToEnd().
                            Split(',').
                            Select((str) => double.Parse(str, System.Globalization.CultureInfo.InvariantCulture.NumberFormat)).
                            Skip(offset).
                            Take(Grid.CellPartitioning.LocalLength * noOfNodesPerCell).
                            ToArray();
                    };
                    output.CopyTo(arrOut.Storage, 0);

                }
                );
            WorkingSet.Momentum[0].ProjectField(1.0, func, scheme);

            p = Control.EnergyDegree;
            scheme = new CellQuadratureScheme(true);
            order = WorkingSet.Energy.Basis.Degree * 2 + 2;
            noOfNodesPerCell = scheme.Compile(GridData, order).First().Rule.NoOfNodes;
            offset = Grid.CellPartitioning.i0 * noOfNodesPerCell;

            func = new ScalarFunction(
                (MultidimensionalArray arrIn, MultidimensionalArray arrOut) => {
                    string filename = "rhoE" + initial + "_" + K.ToString() + "_" + p.ToString() + ".txt";
                    double[] output;
                    using (System.IO.StreamReader sr = new System.IO.StreamReader(pathToData + filename)) {
                        output = sr.ReadToEnd().
                            Split(',').
                            Select((str) => double.Parse(str, System.Globalization.CultureInfo.InvariantCulture.NumberFormat)).
                            Skip(offset).
                            Take(Grid.CellPartitioning.LocalLength * noOfNodesPerCell).
                            ToArray();
                    };
                    output.CopyTo(arrOut.Storage, 0);

                }
                );
            WorkingSet.Energy.ProjectField(func);


            p = Control.DensityDegree;
            scheme = new CellQuadratureScheme(true);
            order = WorkingSet.Density.Basis.Degree * 2 + 2;
            noOfNodesPerCell = scheme.Compile(GridData, order).First().Rule.NoOfNodes;
            offset = Grid.CellPartitioning.i0 * noOfNodesPerCell;

            func = new ScalarFunction(
                (MultidimensionalArray arrIn, MultidimensionalArray arrOut) => {
                    string filename = "rho" + initial + "_" + K.ToString() + "_" + p.ToString() + ".txt";
                    double[] output;
                    using (System.IO.StreamReader sr = new System.IO.StreamReader(pathToData + filename)) {
                        output = sr.ReadToEnd().
                            Split(',').
                            Select((str) => double.Parse(str, System.Globalization.CultureInfo.InvariantCulture.NumberFormat)).
                            Skip(offset).
                            Take(Grid.CellPartitioning.LocalLength * noOfNodesPerCell).
                            ToArray();
                    };
                    output.CopyTo(arrOut.Storage, 0);

                }
                );

            WorkingSet.Density.ProjectField(func);
            WorkingSet.UpdateDerivedVariables(this, CellMask.GetFullMask(GridData));

            //string guidString = "00000000-0000-0000-0000-000000000000";

            //switch (K) {
            //    case 32:
            //        guidString = "6725b9fc-d072-44cd-ae72-196e8a692735";
            //        break;
            //    case 64:
            //        guidString = "cb714aec-4b4a-4dae-9e70-32861daa195f";
            //        break;
            //    case 128:
            //        guidString = "678dc736-bbf6-4a1e-86eb-4f80b2354748";
            //        break;
            //}


            //Guid tsGuid = new Guid(guidString);
            //var db = this.GetDatabase();
            //string fieldName = "rho";
            //ITimestepInfo tsi = db.Controller.DBDriver.LoadTimestepInfo(tsGuid, null, db);
            //previousRho = (SinglePhaseField)db.Controller.DBDriver.LoadFields(tsi, GridData, new[] { fieldName }).Single().CloneAs();
            //diffRho = WorkingSet.Density.CloneAs();
            //diffRho.Clear();


            //fieldName = "rhoE";
            //previousRhoE = (SinglePhaseField)db.Controller.DBDriver.LoadFields(tsi, GridData, new[] { fieldName }).Single().CloneAs();
            //diffRhoE = WorkingSet.Energy.CloneAs();
            //diffRhoE.Clear();

            //fieldName = "m0";
            //previousM0 = (SinglePhaseField)db.Controller.DBDriver.LoadFields(tsi, GridData, new[] { fieldName }).Single().CloneAs();
            //diffM0 = WorkingSet.Momentum[0].CloneAs();
            //diffM0.Clear();

            //diffRho.Identification = "diffRho";
            //diffM0.Identification = "diffM0";
            //diffRhoE.Identification = "diffRhoE";
            //previousRho.Identification = "Rho_analy";
            //previousM0.Identification = "M0_analy";
            //previousRhoE.Identification = "RhoE_analy";
            //m_IOFields.Add(diffM0);
            //m_IOFields.Add(diffRho);
            //m_IOFields.Add(diffRhoE);
            //m_IOFields.Add(previousM0);
            //m_IOFields.Add(previousRho);
            //m_IOFields.Add(previousRhoE);


        }


    }
}
