using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Quadrature;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Tecplot;
using ilPSP;
using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Runtime.CompilerServices;
using System.Text;
using System.Threading.Tasks;

namespace ZwoLevelSetSolver.ZLSinSituPostProcessing
{
    /// <summary>
    /// Post-processing specific to <see cref="ConvergenceTest"/>
    /// </summary>
    [Serializable]
    public class ConvergenceTestEngergieLogging : ZLSinSituPostProcessingModule
    {

        public ConvergenceTestEngergieLogging()
        {

        }
        /// <summary>
        /// 
        /// </summary>
        public const string LogfileName = "BenchmarkQuantities_ConvergenceTest_ZLS";

        /// <summary>
        /// hard-coded name for the ConvergenceTest_ZLS
        /// </summary>
        protected override string LogFileName => LogfileName;

        /// <summary>
        /// Header for the ConvergenceTest_ZLS log
        /// </summary>
        protected override void WriteHeader(TextWriter textWriter)
        {
            string header = String.Format("{0}\t{1}\t{2}", "#timestep", "kineticEnergy_L2", "elasticEnergy_L2", "Sum");
            textWriter.WriteLine(header);
            Log.Flush();
        }

        /// <summary>
        /// compute and log
        /// </summary>
        protected override void PerformTimestepPostProcessing(int TimestepNo, double phystime)
        {
            int degree = 2;
            int j0 = 0;
            //MultidimensionalArray displacementXGradientResult = MultidimensionalArray.Create(GridData.Grid.NumberOfCells, 1, 2);
            //MultidimensionalArray displacementYGradientResult = MultidimensionalArray.Create(GridData.Grid.NumberOfCells, 1, 2);
            MultidimensionalArray elasticTensor = ElasticTensor(2, 1);

            Basis basis = new Basis(CurrentVel[0].GridDat.Grid, degree);
            Basis basisPressure = new Basis(CurrentPressure.GridDat.Grid, degree - 1);
            CellMask cellMask = CellMask.GetFullMask(CurrentVel[0].GridDat);
            //NodeSet nodes = new NodeSet(GridData.Grid.GetRefElement(0), new double[] { 0, 0 });
            //nodes.LockForever();

            DGField elasticEnergyField = new SinglePhaseField(basis, "ElasticEnergy");
            DGField kineticEnergyX = new SinglePhaseField(basis, "KineticEnergyX");
            DGField kineticEnergyY = new SinglePhaseField(basis, "KineticEnergyY");
            DGField kineticEnergyField = new SinglePhaseField(basis, "KineticEnergy");
            //Gradient U
            DGField dU1dx = new SinglePhaseField(basis, "GradientU1x");
            DGField dU1dy = new SinglePhaseField(basis, "GradientU1y");
            DGField dU2dx = new SinglePhaseField(basis, "GradientU2x");
            DGField dU2dy = new SinglePhaseField(basis, "GradientU2y");
            //Strain tensor
            DGField epsilon11 = new SinglePhaseField(basis, "epsilon11");
            DGField epsilon12 = new SinglePhaseField(basis, "epsilon12");
            DGField epsilon21 = new SinglePhaseField(basis, "epsilon21");
            DGField epsilon22 = new SinglePhaseField(basis, "epsilon22");
            //Stress tesnor
            DGField sigma11 = new SinglePhaseField(basis, "sigma11");
            DGField sigma12 = new SinglePhaseField(basis, "sigma12");
            DGField sigma21 = new SinglePhaseField(basis, "sigma21");
            DGField sigma22 = new SinglePhaseField(basis, "sigma22");

            List<DGField> velocity = new List<DGField>();
            List<DGField> displacement = new List<DGField>();
            DGField pressure = new SinglePhaseField(basisPressure, "Pressure");
            velocity.Add(new SinglePhaseField(basis, "VelocityX"));
            velocity.Add(new SinglePhaseField(basis, "VelocityY"));
            displacement.Add(new SinglePhaseField(basis, "DisplacementX"));
            displacement.Add(new SinglePhaseField(basis, "DisplacementY"));

            VectorField<DGField> gradientU1 = new VectorField<DGField>(dU1dx, dU1dy);
            VectorField<DGField> gradientU2 = new VectorField<DGField>(dU2dx, dU2dy);

            velocity[0].ProjectPow(1, CurrentVel[0], 1);
            velocity[1].ProjectPow(1, CurrentVel[1], 1);
            displacement[0].ProjectPow(1, CurrentDisplacement[0], 1);
            displacement[1].ProjectPow(1, CurrentDisplacement[1], 1);
            pressure.ProjectPow(1, CurrentPressure, 1);
            gradientU1.Gradient(1.0, displacement[0]);
            gradientU2.Gradient(1.0, displacement[1]);
            epsilon11.ProjectPow(1, dU1dx, 1);
            epsilon12.ProjectPow(0.5, dU1dy + dU2dx, 1);
            epsilon21.ProjectPow(0.5, dU1dy + dU2dx, 1);
            epsilon22.ProjectPow(1, dU2dy, 1);
            sigma11.ProjectPow(1, epsilon11 + 2 * epsilon11 + 2 * epsilon22, 1);
            sigma12.ProjectPow(1, epsilon12, 1);
            sigma21.ProjectPow(1, epsilon21, 1);
            sigma22.ProjectPow(1, epsilon22 + 2 * epsilon11 + 2 * epsilon22, 1);
            //sigma11.ProjectPow(1, 2 * (epsilon11 + epsilon22) + 2 * epsilon11, 1);
            //sigma12.ProjectPow(1, epsilon12, 1);
            //sigma21.ProjectPow(1, epsilon21, 1);
            //sigma22.ProjectPow(1, 2 * (epsilon11 + epsilon22) + 2 * epsilon22, 1);

            ElasticEnergyField(elasticTensor, epsilon11, epsilon12, epsilon21, epsilon22, sigma11, sigma12, sigma21, sigma22, elasticEnergyField);
            //double elasticEnergy = elasticEnergyField.GetIntegral(cellMask);
            double elasticEnergy = elasticEnergyField.L2Norm();

            kineticEnergyX.ProjectPow(0.5, CurrentVel[0], 2);
            kineticEnergyY.ProjectPow(0.5, CurrentVel[1], 2);
            kineticEnergyField.AccLaidBack(1, kineticEnergyX + kineticEnergyY);
            double kineticEnergy = kineticEnergyField.L2Norm();
            //double kineticEnergy = kineticEnergyField.GetIntegral(cellMask);

            List<DGField> energy = new List<DGField>();
            energy.Add(kineticEnergyField);
            energy.Add(elasticEnergyField);
            Tecplot plt1 = new Tecplot(GridData, true, false, (uint)2);
            plt1.PlotFields("Energy", phystime, energy);

            List<DGField> gradient = new List<DGField>();
            gradient.Add(dU1dx);
            gradient.Add(dU1dy);
            gradient.Add(dU2dx);
            gradient.Add(dU2dy);
            Tecplot plt2 = new Tecplot(GridData, true, false, (uint)2);
            plt2.PlotFields("Gradient", phystime, gradient);
            //displacement[0].EvaluateGradient(j0, displacement[0].GridDat.Grid.NumberOfCells, nodes, displacementXGradientResult);
            //displacement[1].EvaluateGradient(j0, displacement[1].GridDat.Grid.NumberOfCells, nodes, displacementYGradientResult);
            //MultidimensionalArray strainTensor = StrainTensor(displacementXGradientResult, displacementYGradientResult, nodes);
            //MultidimensionalArray elasticEnergy = ElasticEnergy(elasticTensor, strainTensor, nodes);
            //elasticEneryField.ProjectNodal(1, , nodes);
            //double elasticEnergy_L2 = elasticEneryField.L2Norm();

            string line = String.Format("{0}\t{1}\t{2}\t{3}\t{4}", TimestepNo, phystime, kineticEnergy, elasticEnergy, kineticEnergy + elasticEnergy);

            Log.WriteLine(line);
            Log.Flush();
        }
        private void ElasticEnergyField(MultidimensionalArray elasticTensor, DGField epsilon11, DGField epsilon12, DGField epsilon21, DGField epsilon22, DGField sigma11, DGField sigma12, DGField sigma21, DGField sigma22, DGField elasticEnergyField)
        {
            DGField energy1 = new SinglePhaseField(elasticEnergyField.Basis, "energy1");
            DGField energy2 = new SinglePhaseField(elasticEnergyField.Basis, "energy2");
            DGField energy3 = new SinglePhaseField(elasticEnergyField.Basis, "energy3");
            DGField energy4 = new SinglePhaseField(elasticEnergyField.Basis, "energy4");
            //temp.ProjectProduct(0.5 * elasticTensor[0, 0, 0, 0], epsilon11, epsilon11);
            //elasticEnergyField.ProjectPow(1, elasticEnergyField + temp, 1);
            //temp.ProjectProduct(0.5 * elasticTensor[1, 0, 0, 0], epsilon21, epsilon11);
            //elasticEnergyField.ProjectPow(1, elasticEnergyField + temp, 1);
            //temp.ProjectProduct(0.5 * elasticTensor[0, 1, 0, 0], epsilon12, epsilon11);
            //elasticEnergyField.ProjectPow(1, elasticEnergyField + temp, 1);
            //temp.ProjectProduct(0.5 * elasticTensor[0, 0, 1, 0], epsilon11, epsilon21);
            //elasticEnergyField.ProjectPow(1, elasticEnergyField + temp, 1);
            //temp.ProjectProduct(0.5 * elasticTensor[0, 0, 0, 1], epsilon11, epsilon12);
            //elasticEnergyField.ProjectPow(1, elasticEnergyField + temp, 1);
            //temp.ProjectProduct(0.5 * elasticTensor[1, 1, 0, 0], epsilon22, epsilon11);
            //elasticEnergyField.ProjectPow(1, elasticEnergyField + temp, 1);
            //temp.ProjectProduct(0.5 * elasticTensor[1, 0, 1, 0], epsilon21, epsilon21);
            //elasticEnergyField.ProjectPow(1, elasticEnergyField + temp, 1);
            //temp.ProjectProduct(0.5 * elasticTensor[1, 0, 0, 1], epsilon21, epsilon12);
            //elasticEnergyField.ProjectPow(1, elasticEnergyField + temp, 1);
            //temp.ProjectProduct(0.5 * elasticTensor[0, 1, 1, 0], epsilon12, epsilon21);
            //elasticEnergyField.ProjectPow(1, elasticEnergyField + temp, 1);
            //temp.ProjectProduct(0.5 * elasticTensor[0, 1, 0, 1], epsilon12, epsilon12);
            //elasticEnergyField.ProjectPow(1, elasticEnergyField + temp, 1);
            //temp.ProjectProduct(0.5 * elasticTensor[0, 0, 1, 1], epsilon11, epsilon22);
            //elasticEnergyField.ProjectPow(1, elasticEnergyField + temp, 1);
            //temp.ProjectProduct(0.5 * elasticTensor[1, 1, 1, 0], epsilon22, epsilon21);
            //elasticEnergyField.ProjectPow(1, elasticEnergyField + temp, 1);
            //temp.ProjectProduct(0.5 * elasticTensor[1, 1, 0, 1], epsilon22, epsilon12);
            //elasticEnergyField.ProjectPow(1, elasticEnergyField + temp, 1);
            //temp.ProjectProduct(0.5 * elasticTensor[1, 0, 1, 1], epsilon21, epsilon22);
            //elasticEnergyField.ProjectPow(1, elasticEnergyField + temp, 1);
            //temp.ProjectProduct(0.5 * elasticTensor[0, 1, 1, 1], epsilon12, epsilon22);
            //elasticEnergyField.ProjectPow(1, elasticEnergyField + temp, 1);
            //temp.ProjectProduct(0.5 * elasticTensor[1, 1, 1, 1], epsilon22, epsilon22);
            //elasticEnergyField.ProjectPow(1, elasticEnergyField + temp, 1);
            energy1.ProjectProduct(1, sigma11, epsilon11);          
            energy2.ProjectProduct(1, sigma12, epsilon12);
            energy3.ProjectProduct(1, sigma21, epsilon21);
            energy4.ProjectProduct(1, sigma22, epsilon22);
            elasticEnergyField.ProjectPow(0.5, energy1 + energy2 + energy3 + energy4, 1);
        }

        private MultidimensionalArray ElasticEnergy(MultidimensionalArray elasticTensor, MultidimensionalArray strainTensor, NodeSet nodes)
        {
            MultidimensionalArray elasticEnergy = MultidimensionalArray.Create(GridData.Grid.NumberOfCells, nodes.NoOfNodes);
            int dimension = this.SolverMain.GridData.SpatialDimension;
            for(int a = 0; a < strainTensor.Lengths[0]; a++)
            {
                for(int b = 0; b < strainTensor.Lengths[1]; b++)
                {
                    for (int i = 0; i < dimension; i++)
                    {
                        for (int j = 0; j < dimension; j++)
                        {
                            for (int k = 0; k < dimension; k++)
                            {
                                for (int l = 0; l < dimension; l++)
                                {
                                    double elasticConst = elasticTensor[i, j, k, l];
                                    double strain1 = strainTensor[a, b, k, l];
                                    double strain2 = strainTensor[a, b, i, j];
                                    elasticEnergy[a, b] += 0.5 * elasticConst * strain1 * strain2;
                                }
                            }
                        }
                    }
                }
            }
            return elasticEnergy;
        }

        private MultidimensionalArray StrainTensor(MultidimensionalArray displacementXGradient, MultidimensionalArray displacementYGradient, NodeSet nodes)
        {
            int dimension = this.SolverMain.GridData.SpatialDimension;
            MultidimensionalArray strainTensor = MultidimensionalArray.Create(GridData.Grid.NumberOfCells, nodes.NoOfNodes, dimension, dimension);
            for (int a = 0; a < strainTensor.Lengths[0]; a++)
            {
                for (int b = 0; b < strainTensor.Lengths[1]; b++)
                {
                    if (dimension == 2)
                    {
                        strainTensor[a, b, 0, 0] = 2 * displacementXGradient[a, b, 0];
                        strainTensor[a, b, 0, 1] = displacementXGradient[a, b, 1] + displacementYGradient[a, b, 0];
                        strainTensor[a, b, 1, 0] = displacementXGradient[a, b, 1] + displacementYGradient[a, b, 0];
                        strainTensor[a, b, 1, 1] = 2 * displacementYGradient[a, b, 1];
                    }
                    else
                    {
                        throw new NotImplementedException();
                    }
                }
            }
            return strainTensor;
        }

        private MultidimensionalArray ElasticTensor(double lame1, double lame2)
        {
            int dimension = this.SolverMain.GridData.SpatialDimension;
            MultidimensionalArray elasticTensor = MultidimensionalArray.Create(dimension, dimension, dimension, dimension);
            for (int i = 0; i < dimension; i++)
            {
                for (int j = 0; j < dimension; j++)
                {
                    for (int k = 0; k < dimension; k++)
                    {
                        for (int l = 0; l < dimension; l++)
                        {
                            elasticTensor[i, j, k, l] = lame1* KroneckerDelta(i,j,k,l)+lame2*(KroneckerDelta(i,k,j,l)+KroneckerDelta(i,l,j,k));
                        }
                    }
                }
            }
            return elasticTensor;
        }
        
        private int KroneckerDelta(int index1, int index2, int index3, int index4)
        {
            if(index1 == index2 && index3 ==index4)
            {
                return 1;
            }
            else
            {
                return 0;
            }
        }
    }
}
