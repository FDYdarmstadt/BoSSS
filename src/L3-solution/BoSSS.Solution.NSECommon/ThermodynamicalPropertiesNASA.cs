using ilPSP.Tracing;
using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Solution.NSECommon {
    /// <summary>
    /// A class which provides methods for calculating Thermodynamical properties (cp, deltaHf, S) for individual species.
    /// Based on the NASA-polynomials used in the GRI-Mech kinetic mechanism. 
    /// http://combustion.berkeley.edu/gri_mech/version30/files30/thermo30.dat
    /// </summary>
    public class ThermodynamicalProperties {
        static private NumberFormatInfo provider;
        static ThermodynamicalProperties() {

            string AllData = Resource1.thermo30_dat;

            using (MemoryStream stream = new MemoryStream(Encoding.UTF8.GetBytes(AllData))) {
                using (StreamReader sr = new StreamReader(stream)) {


                    // Read number of lines 

                    int lineCount = getTotalNumberOfLines(sr);

                    provider = new NumberFormatInfo();
                    provider.NumberDecimalSeparator = ".";
                    provider.NumberGroupSeparator = ",";

                    sr.DiscardBufferedData();
                    sr.BaseStream.Seek(0, System.IO.SeekOrigin.Begin); // Go back to beginning of the file 

                    sr.ReadLine();
                    var data = sr.ReadLine();
                    var split = data.Split(' ');

                    int i = 0;
                    foreach (string s in split) {
                        if (!(s == "")) {
                            TemperatureLimits[i] = Convert.ToDouble(s, provider);
                            i++;
                        }
                    }


                    CreateAtomicWeigthsDictionary(); //

                    //Loop over all species and store species name and coefficients in the dictionary
                    int no_of_species = (lineCount - 6) / 4;//  Total number of species present in thermo.data
                    for (int n = 1; n < no_of_species + 1; n++) {
                        int linenumber = 4 * (n - 1) + 7;
                        string name = getname(sr, linenumber - 1);   //Get the name of the species             
                        double[][] coefficients = getCoefficients(sr, linenumber);   //Get all coefficients from this species
                        coefficientsDict.Add(name, coefficients);

                        double MW = calculateMW(name);
                        molecularWeightDict.Add(name, MW);
                    }

                }
            }
        }



        static int getTotalNumberOfLines(StreamReader sr) {
            int lines = 0;
            string line = "dummy";
            while (line != "END") {
                line = sr.ReadLine();
                lines++;
            }
            return lines;
        }
        static void CreateAtomicWeigthsDictionary() {
            // Calculate molecular weight and store in dictionary
            char[] a = new char[] { 'H', 'C', 'N', 'O', 'A' };
            double[] atomicWeigths = new double[] { 1, 12, 14, 16, 40 };
            for (int j = 0; j < atomicWeigths.Length; j++) {
                atomicWeigthsDictionary.Add(a[j], atomicWeigths[j]);
            }
        }



        static Dictionary<char, double> atomicWeigthsDictionary = new Dictionary<char, double>();
        /// <summary>
        /// Calculates molecular weight based on the molecular formula of the species
        /// Example: Molecular weight of CH4 is calculated as: AtomicWeigth(C)*1 + AtomicWeigth(H)*4
        /// </summary>
        /// <param name="name"></param>
        /// <returns></returns>
        static public double calculateMW(string name) {
            char[] charName = name.ToCharArray();
            int nameLength = name.Length;

            double MW = 0;
            for (int i = 0; i < charName.Length; i++) {
                char c = charName[i];
                char charnext = (i == nameLength - 1) ? '1' : charName[i + 1]; // Gets next char. If c is the last element, next char is empty
                if (charnext == '(' || c == 'R') // Ugly bug solution for two special cases
                    break;

                if (Char.IsLetter(c)) {
                    double AW = atomicWeigthsDictionary[c];
                    int NumberOfAtoms;
                    if (Char.IsNumber(charnext)) { // if the next char is a number, it has to be >2 and is the number of atoms
                        NumberOfAtoms = int.Parse(charnext.ToString());
                    } else if ((Char.IsLetter(charnext) || charnext == ' ')) { // if the next char is another letter, or nothing, the number of atoms is one.
                        NumberOfAtoms = 1;
                    } else {
                        throw new Exception("??????");
                    }
                    MW = MW + AW * NumberOfAtoms;
                }
            }
            return MW;
        }


        /// <summary>
        /// Get name of the species starting in line = linenumber
        /// </summary>
        /// <param name="sr"></param>
        /// <param name="linenumber"></param>
        /// <returns></returns>
        static public string getname(StreamReader sr, int linenumber) {
            string name = "";
            rewind(sr, linenumber);
            char[] charline = sr.ReadLine().ToCharArray();
            for (int i = 0; i < charline.Length; i++) {
                char c = charline[i];
                if (!(c == ' ')) { // if c is not whitespace, is part of the name
                    name = name + c;
                } else {
                    break;
                }
            }


            return name;
        }

        /// <summary>
        /// Gets the coefficients for the species starting in a given line
        /// </summary>
        /// <param name="sr"></param>
        /// <param name="line_number"></param>
        /// <returns></returns>
        static double[][] getCoefficients(StreamReader sr, int line_number) {
            rewind(sr, line_number); // now we are positioned in the line number here specified
            string line1 = sr.ReadLine();
            string line2 = sr.ReadLine();
            string line3 = sr.ReadLine();

            double[] allCoefficients = new double[14];

            int charindex = 0;
            for (int k = 0; k < 5; k++) { // loop over all 14 coefficients, stored in 3 lines. 5 coef in line1, 5 coef in line2, and 4 in line 3            
                allCoefficients[k] = Convert.ToDouble(line1.Substring(charindex, 15), provider);
                charindex = charindex + 15;
            }
            charindex = 0;
            for (int k = 5; k < 10; k++) { // loop over all 14 coefficients, stored in 3 lines. 5 coef in line1, 5 coef in line2, and 4 in line 3       
                allCoefficients[k] = Convert.ToDouble(line2.Substring(charindex, 15), provider);
                charindex = charindex + 15;
            }
            charindex = 0;
            for (int k = 10; k < 14; k++) { // loop over all 14 coefficients, stored in 3 lines. 5 coef in line1, 5 coef in line2, and 4 in line 3          
                allCoefficients[k] = Convert.ToDouble(line3.Substring(charindex, 15), provider);
                charindex = charindex + 15;
            }

            double[][] coefficients = new double[2][];
            coefficients[0] = allCoefficients.Take(7).ToArray();// low temperature coefficients
            coefficients[1] = allCoefficients.Skip(7).Take(7).ToArray(); // High temperature coefficients

            return coefficients;
        }


        /// <summary>
        /// moves the pointer of the reader to a defined line
        /// </summary>
        /// <param name="sr"></param>
        /// <param name="line"></param>
        static void rewind(StreamReader sr, int line) {
            sr.DiscardBufferedData();
            sr.BaseStream.Seek(0, System.IO.SeekOrigin.Begin); // Go back to beginning of the file 
            int actualLine = 0;
            while (actualLine < line - 1) {
                sr.ReadLine();
                actualLine++;
            }
        }

        /// <summary>
        ///  lower, middle and upper temperature limit
        ///  /// </summary>
        static private double[] TemperatureLimits = new double[3];

        /// <summary>
        /// for a given component name two double arrays with coefficients
        /// </summary>
        static private Dictionary<string, double[][]> coefficientsDict = new Dictionary<string, double[][]>();

        /// <summary>
        /// for a given component name two double arrays with coefficients
        /// </summary>
        static private Dictionary<string, double> molecularWeightDict = new Dictionary<string, double>();


        /// <summary>
        /// Returns the cp value for a given temperature <see cref="T"/>
        /// Results in Joule/mol K
        /// </summary>
        /// <param name="name">Name of species</param>
        /// <param name="T"> Temperature </param>
        /// <returns></returns>
        public double getCp(string name, double T) {
            //using (var tr = new FuncTrace()) {
            double[] coefficients;
            if (/*T >= TemperatureLimits[0] - 200.0 &&*/ T < TemperatureLimits[1]) { // Lower range, with a threshold 5.0 K
                coefficients = coefficientsDict[name][1];
            } else if (T >= TemperatureLimits[1] && T <= TemperatureLimits[2] + 5.0) { // Higher range
                coefficients = coefficientsDict[name][0];
            } else {
                throw new ArgumentOutOfRangeException("Temperature for calculation of cp is out of bounds. The used temperature is" + T);
            }

            //Force temperature
            if (T < TemperatureLimits[0])
                T = TemperatureLimits[0];

            double R = 8.314; // KJ /Kmol K
            double MW = molecularWeightDict[name]; // Kg/Kmol
            double cp;
            //using (new BlockTrace("Heat capacity calculation", tr)) {
            cp = (
            coefficients[0] +
            coefficients[1] * T +
            coefficients[2] * T * T +
            coefficients[3] * T * T * T +
            coefficients[4] * T * T * T * T) * R / MW;
            //}
            return cp; // KJ/(kg K
            //}

        }

        /// <summary>
        /// Returns the entropy 
        /// Results in Joule/mol K
        /// </summary>
        /// <param name="name">Name of species</param>
        /// <param name="T"> Temperature </param>
        /// <returns></returns>
        public double getS(string name, double T) {

            double[] coefficients;
            if (T >= TemperatureLimits[0] && T < TemperatureLimits[1]) { // Lower range
                coefficients = coefficientsDict[name][1];
            } else if (T >= TemperatureLimits[1] && T <= TemperatureLimits[2]) { // Higher range
                coefficients = coefficientsDict[name][0];
            } else {
                throw new ArgumentOutOfRangeException("Temperature for calculation of S is out of bounds");
            }
            double R = 8.314; // J /mol K

            double cp = (
                coefficients[0] * Math.Log(T) +
                coefficients[1] * T +
                coefficients[2] * T * T / 2.0 +
                coefficients[3] * T * T * T / 3.0 +
                coefficients[4] * T * T * T * T / 4 +
                coefficients[6]) * R;
            return cp;
        }




        /// <summary>
        /// Input: MassFractions and name of the components
        /// Output: The heat capacity of the mixture
        /// </summary>
        /// <returns></returns>
        public double Calculate_Cp_Mixture(double[] MassFractions, string[] Names, double Temperature) {
            double cpMixture = 0.0;
            for (int i = 0; i < MassFractions.Length; i++) {
                double cp_i = getCp(Names[i], Temperature);
                cpMixture = cpMixture + cp_i * MassFractions[i];
            }
            return cpMixture ;
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="name"></param>
        /// <returns></returns>
        public double getMW(string name) {
            double MW = molecularWeightDict[name];
            return MW;
        }

    }

}
