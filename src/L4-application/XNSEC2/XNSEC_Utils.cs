using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Tecplot;
using BoSSS.Solution.XNSECommon;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using BoSSS.Solution.Utils;
using BoSSS.Solution.NSECommon;
using BoSSS.Foundation.Grid;
using System.Globalization;

namespace BoSSS.Application.XNSEC {

    public partial class XNSEC {
        private void PrintConfiguration() {
            XNSEC_OperatorConfiguration config = new XNSEC_OperatorConfiguration(this.Control);

            PropertyInfo[] properties = typeof(XNSEC_OperatorConfiguration).GetProperties();
            foreach(PropertyInfo property in properties) {
                if(property.PropertyType == typeof(bool)) {
                    bool s = (bool)property.GetValue(config);
                    Console.WriteLine("     {0,-30}:{1,3}", property.Name, "[" + (s == true ? "x" : " ") + "]");
                }
            }
        }
        private DGField FindField(DGField[] Fields, string identification) {
            DGField field = null;
            foreach(DGField f in Fields) {
                if(f.Identification == identification) {
                    field = f;
                }
            }
            if(field == null)
                throw new ArgumentOutOfRangeException("Field with name" + identification + "not found");
            return field;
            //var ret = ((DGField)((field as XDGField).GetSpeciesShadowField("A"))); // take species A
            //return ret;
        }

        /// <summary>
        /// Instantiation of fields according to the domain variable names in the spatial operator.
        /// </summary>
        private IEnumerable<DGField> InstantiateErrorFields() {
            var DomNames = this.Operator.DomainVar;
            var ret = new DGField[DomNames.Count];
            for(int i = 0; i < DomNames.Count; i++) {
                string Name = DomNames[i];

                var fopts = this.Control.FieldOptions.Where(kv => kv.Key.WildcardMatch(Name)).SingleOrDefault().Value;
                if(fopts.Degree < 0) {
                    throw new ApplicationException($"Missing specification of DG degree for field {Name} in control object.");
                }

                var xb = new XDGBasis(this.LsTrk, fopts.Degree);
                ret[i] = new XDGField(xb, "Err_" + Name);
            }
            return ret;
        }

        /// <summary>
        /// Instantiation of fields according to the domain variable names in the spatial operator.
        /// </summary>
        private IEnumerable<DGField> InstantiateAnalyticalSolFields() {
            var DomNames = this.Operator.DomainVar;
            var ret = new DGField[DomNames.Count];
            for(int i = 0; i < DomNames.Count; i++) {
                string Name = DomNames[i];

                var fopts = this.Control.FieldOptions.Where(kv => kv.Key.WildcardMatch(Name)).SingleOrDefault().Value;
                if(fopts.Degree < 0) {
                    throw new ApplicationException($"Missing specification of DG degree for field {Name} in control object.");
                }

                var xb = new XDGBasis(this.LsTrk, fopts.Degree);
                ret[i] = new XDGField(xb, Name + "_an");
            }
            return ret;
        }

        /// <summary>
        /// Calculates the error in case that an analytical solution is avaliable.
        /// </summary>
        ///
        public void CalcErrors() {
            var allFields = this.m_RegisteredFields.ToArray();

            var DomNames = this.Operator.DomainVar;
            foreach(var name in DomNames) {
                var field = FindField(allFields, name);
                var anField = FindField(allFields, name + "_an");
                var errorField = FindField(allFields, "Err_" + name);
                errorField.Clear();
                errorField.CoordinateVector.AccV(1.0, anField.CoordinateVector);
                errorField.CoordinateVector.AccV(-1.0, field.CoordinateVector);

                //if(name == BoSSS.Solution.NSECommon.VariableNames.Pressure) {
                //    double meanpressure = errorField.GetMeanValueTotal(null);
                //    errorField.AccConstant(-1.0* meanpressure);
                //}

                //string[] species = new string[] { "A", "B" };
                //foreach(var sp in species) {
                //errorField.GetSpeciesShadowField(sp).Acc(1.0,anfi)
                //        }

                if (field.Identification == Solution.NSECommon.VariableNames.Pressure)
                    //base.QueryHandler.ValueQuery(errorField.Identification,
                    //    errorField.L2ErrorNoMean(((Func<double[], double>)(X => 0.0)).Vectorize(), errorField.Basis.Degree * 2),
                    //    true);
                    base.QueryHandler.ValueQuery(errorField.Identification,                        
                     errorField.L2ErrorNoMean(((Func<double[], double>)(X => 0.0)).Vectorize(), errorField.Basis.Degree * 2),
                     true);
                else
                    base.QueryHandler.ValueQuery(errorField.Identification, errorField.L2Norm(), true);

            }
        }

        private void PlotNewtonIterationsHack(double t, IReadOnlyDictionary<string, DGField> DomainVarFields, IReadOnlyDictionary<string, DGField> ParameterVarFields) {
            if(Control.PlotNewtonIterations)
                Tecplot.PlotFields(base.m_RegisteredFields, "NewtonSolverStep" + hack_TimestepIndex, 0.0, 3);
            hack_TimestepIndex++;
        }

        #region NusseltNumber

        NusseltNumber NusseltNum;
        TextWriter LogNusselt;

        /// <summary>
        /// Calculate Nusselt number for Low-Mach number flows.
        /// </summary>
        /// <param name="Timestep"></param>
        /// <param name="GridDat"></param>
        /// <param name="Temperature"></param>
        /// <param name="SolverConf"></param>
        public double[] CalculateNusselt(int Timestep, IGridData GridDat, SinglePhaseField Temperature, XNSEC_Control SolverConf) {

            // Initialize calculation of Nusselt number.
            if(NusseltNum == null) {
                NusseltNum = new NusseltNumber(GridDat, Temperature, EoS_A, SolverConf.EdgeTagsNusselt);

                if((base.MPIRank == 0) && (CurrentSessionInfo.ID != Guid.Empty)) {
                    LogNusselt = base.DatabaseDriver.FsDriver.GetNewLog("Nusselt", CurrentSessionInfo.ID);
                    LogNusselt.Write("Timestep");
                    for(int bc = 0; bc < NusseltNum.Nusselt.Length; bc++)
                        LogNusselt.Write("\t" + SolverConf.EdgeTagsNusselt[bc]);
                    LogNusselt.WriteLine();
                }
            }

            //// Calculate Nusselt number 
            //Dictionary<string, double[]> Nusselts = new Dictionary<string, double[]>();            
            //for(int i = 0; i < SolverConf.EdgeTagsNusselt.Length; i++) {
            //    Nusselts.Add(SolverConf.EdgeTagsNusselt[i], NusseltNum.Nusselt);
            //}

            // Calculate Nusselt number
                NusseltNum.CalculateNusseltNumber();

            // Write result to text file
            if((base.MPIRank == 0) && (LogNusselt != null)) {
                LogNusselt.Write(Timestep.ToString());
                for(int bc = 0; bc < NusseltNum.Nusselt.Length; bc++)
                    LogNusselt.Write("\t" + NusseltNum.Nusselt[bc].ToString("0.0000000000E+00", NumberFormatInfo.InvariantInfo));
                LogNusselt.WriteLine();
                LogNusselt.Flush();
            }
            double[] nusselts = new double[] { NusseltNum.Nusselt[0], NusseltNum.Nusselt[1], NusseltNum.Nusselt[2] };
            return nusselts;
        }
        #endregion NusseltNumber
        public static void DeleteOldPlotFiles() {
            if(ilPSP.Environment.MPIEnv.MPI_Rank == 0) {
                var dir = new DirectoryInfo(Directory.GetCurrentDirectory());
                Console.Write("rm");
                foreach(var pltFile in dir.GetFiles("*.plt").Concat(dir.GetFiles("*.curve"))) {
                    Console.Write(" " + pltFile.Name);
                    pltFile.Delete();
                }
                Console.WriteLine(";");
            }
        }

        public static void DeleteOldTextFiles() {
            if (ilPSP.Environment.MPIEnv.MPI_Rank == 0) {
                var dir = new DirectoryInfo(Directory.GetCurrentDirectory());
                Console.Write("rm");
                foreach (var pltFile in dir.GetFiles("*.txt")) {
                    Console.Write(" " + pltFile.Name);
                    pltFile.Delete();
                }
                Console.WriteLine(";");
            }
        }

    }
}