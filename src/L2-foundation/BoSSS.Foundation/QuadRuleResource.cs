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
using System.Diagnostics;
using System.IO;
using System.Linq;
using System.Text;
using BoSSS.Foundation.Quadrature;
using ilPSP;
using ilPSP.Utils;
using Newtonsoft.Json;
using Newtonsoft.Json.Bson;

namespace BoSSS.Foundation.Grid.RefElements {

    /// <summary>
    /// This class reads quadrature rules, stored in text files
    /// and converts them to a binary format (base64-encoded string)
    /// which can be embedded in the assembly as a resource.
    /// This approach is somewhat cleaner than the traditional one used in BoSSS 
    /// (large C# files with auto-generated code)
    /// since this, at some point ckacks the limits of compiler, JIT and other tools.
    /// </summary>
    static class QuadRuleResource {

        /// <summary>
        /// Imports a quadrature rule from a text file.
        /// </summary>
        static Tuple<int, double[,], double[]> ReadFromTextFile(string FileName) {
            using(StreamReader txt = new StreamReader(FileName)) {
                int Order;
                string header = txt.ReadLine();
                if(!header.StartsWith("order"))
                    throw new IOException();
                Order = int.Parse(header.Replace("order", ""));

                MultidimensionalArray NodesAndWeights = IMatrixExtensions.LoadFromStream(txt);

                int K = NodesAndWeights.NoOfRows;
                int D = NodesAndWeights.NoOfCols - 1;
                double[,] Nodes = NodesAndWeights.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { K - 1, D - 1 }).To2DArray();
                double[] Weigths = NodesAndWeights.ExtractSubArrayShallow(new int[] { 0, D }, new int[] { K - 1, D - 1 }).To1DArray();

                return new Tuple<int, double[,], double[]>(Order, Nodes, Weigths);
            }
        }

        /// <summary>
        /// Decodes quadrature rules from base64 - encoded binary Data.
        /// </summary>
        static public IEnumerable<Tuple<int, double[,], double[]>> DecodeFromBase64(string s) {
            byte[] BinData = Convert.FromBase64String(s);
            MemoryStream ms = new MemoryStream(BinData);
            BsonReader reader = new BsonReader(ms);
            JsonSerializer serializer = new JsonSerializer();

            Container o = serializer.Deserialize<Container>(reader);
            //var qr = (Tuple<int, double[,], double[]>[])serializer.Deserialize(reader, typeof(Tuple<int, double[,], double[]>[]));
            //R.AddRange(qr);

            return o.Quadrules;
        }

        /// <summary>
        /// Exports quadrature rules to human-readable text files.
        /// </summary>
        static void ExportRules(string dirName, string filename, IEnumerable<QuadRule> qrs) {
            int i = 0;
            foreach(var qr in qrs) {
                i++;
                int D = qr.SpatialDim;
                using(StreamWriter txt = new StreamWriter(Path.Combine("..", "..", dirName, filename + i + ".txt"))) {
                    txt.WriteLine("order" + qr.OrderOfPrecision);

                    MultidimensionalArray NodesAndWeigths = MultidimensionalArray.Create(qr.NoOfNodes, qr.SpatialDim + 1);

                    NodesAndWeigths.ExtractSubArrayShallow(new int[] { 0, 0 }, new int[] { qr.NoOfNodes - 1, D - 1 }).Set(qr.Nodes);
                    NodesAndWeigths.ExtractSubArrayShallow(new int[] { 0, D }, new int[] { qr.NoOfNodes - 1, D - 1 }).Set(qr.Weights);

                    NodesAndWeigths.SaveToStream(txt);
                }
            }
        }


        /// <summary>
        /// Seems to be necessary, otherwise Json Deserializer refuses to work.
        /// </summary>
        [Serializable]
        class Container {
            public Tuple<int, double[,], double[]>[] Quadrules;
        }
       
        
        public static void Main(string[] args) {
            //var qr = DecodeFromBase64(Resource.TetraQuadRules_bin);

            //PackQuadrulesBinary("TetraQuadRules");
            //PackQuadrulesBinary("TriangleQuadRules", Triangle.Instance.m_QuadRules.Select(qr => new Tuple<int,double[,], double[]> (qr.OrderOfPrecision,qr.Nodes.To2DArray(),qr.Weights.To1DArray()) ));
            //PackQuadrulesBinary("LineQuadRules");

            //WriteRules("TetraQuadRules", "keast", Tetra.Instance.m_QuadRules);
            //WriteRules("TriangleQuadRules", "Tri", Triangle.Instance.m_QuadRules);
            //WriteRules("LineQuadRules", "Line", Line.Instance.m_QuadRules);

            //var triImported = Import("TriangleQuadRules");
            //compare(Triangle.Instance.m_QuadRules, triImported);
        }


        private static void compare(IEnumerable<QuadRule> A, IEnumerable<Tuple<int, double[,], double[]>> B) {
            var R = new List<Tuple<QuadRule, QuadRule>>();

            foreach(var qrA in A) {
                var qrB = B.Single(_qrB => qrA.NoOfNodes == _qrB.Item2.GetLength(0) && qrA.OrderOfPrecision == _qrB.Item1);

                var Nds = MultidimensionalArray.Create(qrA.Nodes.Lengths);
                Nds.Set(qrA.Nodes);
                Nds.Acc2DArray(-1.0, qrB.Item2);
                double ERRnds = Nds.L2Norm();

                var wgt = qrA.Weights.CloneAs();
                wgt.AccVector(-1.0, qrB.Item3);
                double ERRwgt = wgt.L2Norm();

                Console.WriteLine("order {0} rule error after text in/out, nodes {1}, weights {2};", qrA.OrderOfPrecision, ERRnds, ERRwgt);
            }

        }

        private static void PackQuadrulesBinary(string qrName, IEnumerable<Tuple<int, double[,], double[]>> qrList) {
            MemoryStream ms = new MemoryStream();

            BsonWriter writer = new BsonWriter(ms);
            JsonSerializer serializer = new JsonSerializer();


            serializer.Serialize(writer, new Container() { Quadrules = qrList.ToArray() });
            writer.Flush();

            ms.Flush();
            Debug.Assert(ms.Position == ms.Length);
            byte[] _Data = ms.GetBuffer();
            Debug.Assert(_Data.Length >= ms.Position);

            byte[] BinData = new byte[ms.Length];
            Array.Copy(_Data, 0, BinData, 0, ms.Position);
            Console.WriteLine("Raw bindata: " + (BinData.Length / 1.0e6));

            string EncData = Convert.ToBase64String(BinData);
            //Console.WriteLine(EncData);
            Console.WriteLine("Encoded base64: " + (EncData.Length / 1.0e6));

            // test
            DecodeFromBase64(EncData);

            using(var txt = new StreamWriter(qrName + "_bin.txt")) {
                txt.Write(EncData);
            }
        }

        /// <summary>
        /// Imports quadrature rules from textfile to mem.
        /// </summary>
        private static List<Tuple<int, double[,], double[]>> Import(string qrName) {
            string[] Files = Directory.GetFiles(Path.Combine("..", "..", qrName), "*.txt");
            var qrList = new List<Tuple<int, double[,], double[]>>();
            foreach(string File in Files) {
                Console.WriteLine("reading " + File + "...");
                var Rule = ReadFromTextFile(File);
                qrList.Add(Rule);
            }
            return qrList;
        }
    }
}
