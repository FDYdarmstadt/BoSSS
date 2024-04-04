using ApplicationWithIDT;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.Utils;
using ilPSP.LinSolvers;
using ilPSP.Utils;
using MathNet.Numerics.Providers.LinearAlgebra;
using NUnit.Framework;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BUIDT.Tests {
    public static class UtilityTests {

        [Test]
        public static void ChangeOfPOrderTest() {

            BoSSS.Solution.Application.InitMPI(num_threads: 1);
            #region initialize a field and program object
            var p = new BUIDTMain();
            var C = BUIDTHardCodedControl.AccShock(
                    dbPath: null,
                    MaxIterations: 100,
                    dgDegree: 1,
                    LSDregree: 3,
                    numOfCellsX: 2,
                    numOfCellsY: 1,
                    OptiNumOfCellsX: 1,
                    OptiNumOfCellsY: 1,
                    linearization: Linearization.FD,
                    agg: 0.0,
                    ImmediatePlotPeriod: 1,
                    optiLevelSetType: OptiLevelSetType.SplineLevelSet,
                    getLevelSet: GetLevelSet.FromFunction,
                    applyReInit: true
                    );
            // ### Grid ###
            double xMin = -1.0;
            double xMax = 1.0;
            double yMin = 0.0;
            double yMax = 1.0;
            //Initiate a 2x1 grid
            C.GridFunc = delegate {
                double[] xNodes;
                xNodes = GenericBlas.Linspace(xMin, xMax, 3);

                double[] yNodes = GenericBlas.Linspace(yMin, yMax, 2);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);
                return grid;
            };

            // Initial Value
            C.LevelSetTwoInitialValue = x => 0 * x[0] + 0.5;
            C.LevelSetOneInitialValue = x => 0 * x[0] + 0.5;
            C.InitialValueFunctionsPerSpecies = new Dictionary<string, Func<double[], double>>();

            var a = 1;
            var b = 1;

            var polToTest = delegate (double[] X) {
                if(X[0] >= 0 && X[1] >= 0) {
                    return X[0] + a;
                } else {
                    return b * X[0];
                }
            };

            var projectionPolToTest = delegate (double[] X, string spc) {
                if(X[0] >= 0 && X[1] >= 0) { //cutCell
                    if(spc == "L") {
                        return a + 0.25;
                    } else {
                        return a + 0.75;
                    }
                } else {
                    return - 0.5;
                }
            };
            var projectionExact = delegate (double[] X) {
                if(X[0] >= 0.5){
                    return a + 0.75;
                } else if(X[0]<0.5 && X[0]>=0){
                    return a + 0.25;
                } else {
                    return -0.5;
                }
            };



            C.InitialValueFunctionsPerSpecies.Add("L", polToTest);
            C.InitialValueFunctionsPerSpecies.Add("R", polToTest);



            p.Init(C);
            p.InitializeEverything();
            #endregion

            p.Concentration.GetSpeciesShadowField("L").ProjectField(X => polToTest(X));
            p.Concentration.GetSpeciesShadowField("R").ProjectField(X => polToTest(X));
            p.LevelSet.ProjectField(X => X[0] - 0.5);
            p.LsTrk.UpdateTracker(0.0);

            p.ChangeTOBasisofDegree(0);
            var errorP0Projection = p.Concentration.L2Error(projectionExact);

            //Changing the Basis to deg 3 shouldn't change the values
            p.ChangeTOBasisofDegree(3);
            var errorP3 = p.Concentration.L2Error(projectionExact);

            p.Concentration.GetSpeciesShadowField("L").ProjectField(X => polToTest(X));
            p.Concentration.GetSpeciesShadowField("R").ProjectField(X => polToTest(X));

            //Changing the Basis to deg 3 shouldn't change anything as it is of deg 3
            p.ChangeTOBasisofDegree(3);

            //this shouldn't also change the function
            p.ChangeTOBasisofDegree(1);

            var errorP1 = p.Concentration.L2Error(polToTest);

            Assert.IsTrue(errorP0Projection < 1e-14 && errorP1 < 1e-14 && errorP3 < 1e-14);

        }
        /// <summary>
        /// Test for the Person Sensor
        /// </summary>
        [Test]
        public static void PerssonSensorTest() {

            BoSSS.Solution.Application.InitMPI(num_threads: 1);
            var p = new BUIDTMain();
            var C = BUIDTHardCodedControl.AccShock(
                    dbPath: null,
                    MaxIterations: 100,
                    dgDegree: 1,
                    LSDregree:3,
                    numOfCellsX: 2,
                    numOfCellsY: 1,
                    OptiNumOfCellsX: 1,
                    OptiNumOfCellsY: 1,
                    linearization: Linearization.FD,
                    agg: 0.0,
                    ImmediatePlotPeriod: 1,
                    optiLevelSetType: OptiLevelSetType.SplineLevelSet,
                    getLevelSet: GetLevelSet.FromFunction,
                    applyReInit: true
                    
                    );
            // ### Grid ###
            double xMin = -1.0;
            double xMax = 1.0;
            double yMin = 0.0;
            double yMax = 1.0;
            //Initiate a 2x1 grid
            C.GridFunc = delegate {
                double[] xNodes;
                xNodes = GenericBlas.Linspace(xMin, xMax, 3);

                double[] yNodes = GenericBlas.Linspace(yMin, yMax, 2);
                var grid = Grid2D.Cartesian2DGrid(xNodes, yNodes, periodicX: false, periodicY: false);
                return grid;
            };

            // Initial Value
            C.LevelSetTwoInitialValue = x =>  0*x[0] + 0.5;
            C.LevelSetOneInitialValue = x => 0 * x[0] + 0.5;
            C.InitialValueFunctionsPerSpecies = new Dictionary<string, Func<double[], double>>();

            var a = -0.75;
            var b = 1;

            var polToTest = delegate (double[] X) {
                if(X[0] >= 0 && X[1] >= 0) {
                    return X[0] +a;
                } else {
                    return b * X[0];
                }
            };

            var projectionPolToTest = delegate (double[] X,string spc) {
                if(X[0] >= 0 && X[1] >= 0) { //cutCell
                    if(spc=="L") {
                        return a +0.25;
                    } else {
                        return a + 0.75;
                    }
                } else {
                    return b * 0.5;
                }
            };

            var PerssonSensor = delegate (double[] X, string spc) {
                if(X[0] >= 0 && X[1] >= 0) { //cutCell
                    if(spc == "L") {
                        return Math.Log10(Math.Sqrt(0.25 / (1 + 6 * a + 12 *a * a)) );
                  } else {
                        return 0;
                    }
                } else {
                    return Math.Log10(0.5);
                }
            };



            C.InitialValueFunctionsPerSpecies.Add("L", polToTest);
            C.InitialValueFunctionsPerSpecies.Add("R", polToTest);



            p.Init(C);
            p.InitializeEverything();

            p.Concentration.GetSpeciesShadowField("L").ProjectField(X => polToTest(X));
            p.Concentration.GetSpeciesShadowField("R").ProjectField(X => polToTest(X));
            p.LevelSet.ProjectField(X => X[0] - 0.5);
            p.LsTrk.UpdateTracker(0.0);

            var l2Norm = p.Concentration.L2NormAllSpecies();


            var Projection = p.Concentration.CloneAs();
            Projection.Identification = "Projection";
            Projection.Clear();
            Projection.GetSpeciesShadowField("L").ProjectField(x => projectionPolToTest(x,"L"));
            Projection.GetSpeciesShadowField("R").ProjectField(x => projectionPolToTest(x, "R"));

            var PerssonSensorAnal = p.personField.CloneAs();
            PerssonSensorAnal.Identification = "PerssonSensor";
            PerssonSensorAnal.Clear();
            PerssonSensorAnal.GetSpeciesShadowField("L").ProjectField(x => PerssonSensor(x, "L"));
            PerssonSensorAnal.GetSpeciesShadowField("R").ProjectField(x => PerssonSensor(x, "R"));


            p.GetPerssonSensor(false);
            var fieldPMinus1 = BUIDTMain.GetPMinus1Projection(p.LsTrk, p.Concentration,2);

            var diffAlg = p.Concentration.CloneAs();
            diffAlg.Identification = "diffALg";
            diffAlg.AccLaidBack(-1.0, fieldPMinus1);

            BitArray bA= new BitArray(2);
            bA[1] = true;
            CellMask cm = new CellMask(p.GridData,bA);
            var l2Cell2A =p.Concentration.L2NormSpecies("L", cm);
            var l2Cell2B = p.Concentration.L2NormSpecies("R", cm);

            var l2Cell2A_proj = diffAlg.L2NormSpecies("L", cm);
            var l2Cell2B_proj = diffAlg.L2NormSpecies("R", cm);

            var personCell2A = Math.Log10(l2Cell2A_proj / l2Cell2A);
            var personCell2B = Math.Log10(l2Cell2B_proj / l2Cell2B);


            var diff = PerssonSensorAnal.CloneAs();
            diff.Acc(-1.0, p.personField);
            diff.Identification = "diff";
            var diffnorm = diff.L2NormAllSpecies();

            var tp = new Tecplot(p.GridData, 5);
            tp.PlotFields("PersonTest", 0.0, new DGField[] { diffAlg,fieldPMinus1, Projection,p.Concentration, PerssonSensorAnal, p.personField, diff,p.LevelSet });

            Assert.IsTrue(diff.L2NormAllSpecies()<1e-15);

        }
    }
}
