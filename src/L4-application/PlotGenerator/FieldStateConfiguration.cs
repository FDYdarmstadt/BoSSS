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

using BoSSS.Foundation.IO;
using System;
using System.Collections.Generic;
using System.IO;
using System.Runtime.Serialization;
using System.Xml;
using System.Xml.Serialization;

namespace BoSSS.PlotGenerator {

    /// <summary>
    /// Class representing a configuration of field states. Can be used as an
    /// exchange format for settings for a list of fields
    /// </summary>
    [DataContract]
    public class FieldStateConfiguration {

        /// <summary>
        /// Supported export formats
        /// </summary>
        [DataContract]
        public enum ExportFormats {

            /// <summary>
            /// Create binary Tecplot file
            /// </summary>
            [EnumMember]
            TecPlot = 0,

            /// <summary>
            /// The CGNS format
            /// </summary>
            [EnumMember]
            Cgns = 1,

            /// <summary>
            /// CGNS via Hdf5
            /// </summary>
            [EnumMember]
            Hdf5 = 2,

            /// <summary>
            /// generic ASCII files: see
            /// <see cref="BoSSS.Solution.ASCIIExport.CSVExportDriver"/>
            /// </summary>
            [EnumMember]
            CSV = 3,

            /// <summary>
            /// VisIts curve file: see
            /// <see cref="BoSSS.Solution.ASCIIExport.CurveExportDriver"/>
            /// </summary>
            [EnumMember]
            Curve = 4,

            /// <summary>
            /// Binary stream as it is for example used by Matlab; see
            /// <see cref="BoSSS.Solution.ASCIIExport.BinaryStreamExportDriver"/>
            /// </summary>
            [EnumMember]
            BinaryStream = 5
        }

        /// <summary>
        /// Supported reconstruction types for the originally discontinuous
        /// data
        /// </summary>
        [DataContract]
        public enum ReconstructionTypes {

            /// <summary>
            /// No reconstruction, export as is
            /// </summary>
            [EnumMember]
            None = 0,

            /// <summary>
            /// Build the node-wise average of all values
            /// </summary>
            [EnumMember]
            Average = 1
        }

        /// <summary>
        /// Supported ghost levels
        /// </summary>
        [DataContract]
        public enum GhostLevels {

            /// <summary>
            /// No overlap, export as is
            /// </summary>
            [EnumMember]
            Zero = 0,

            /// <summary>
            /// Mesh overlap of on unit (element)
            /// </summary>
            [EnumMember]
            One = 1
        }

        /// <summary>
        /// The paths to the databases
        /// </summary>
        [DataMember]
        public string[] BasePaths = null;

        /// <summary>
        /// The selected session
        /// </summary>
        [DataMember]
        public Guid SessionGuid = Guid.Empty;

        /// <summary>
        /// The fields selected for export
        /// </summary>
        [DataMember]
        public List<string> FieldNames = new List<string>();

        /// <summary>
        /// The time-step numbers of the time-steps
        /// </summary>
        [DataMember]
        public List<TimestepNumber> TimeSteps = new List<TimestepNumber>();

        /// <summary>
        /// The selected value for super sampling
        /// </summary>
        [DataMember]
        public int SuperSampling = 0;

        /// <summary>
        /// The selected reconstruction type
        /// </summary>
        [DataMember]
        public ReconstructionTypes ReconstructionType;

        /// <summary>
        /// The selected ghost level
        /// </summary>
        [DataMember]
        public GhostLevels GhostLevel;

        /// <summary>
        /// The number of processes that should be used for the export
        /// </summary>
        [DataMember]
        public uint NumberOfProcesses = 1;

        /// <summary>
        /// The selected export format
        /// </summary>
        [DataMember]
        public ExportFormats ExportFormat;

        /// <summary>
        /// Uses a <see cref="DataContractSerializer"/> to serialize the given
        /// <paramref name="config"/>
        /// </summary>
        /// <param name="fileName">The name of the file to be created</param>
        /// <param name="config">The configuration to be serialized</param>
        public static void Serialize(string fileName, FieldStateConfiguration config) {
            DataContractSerializer serializer =
                new DataContractSerializer(typeof(FieldStateConfiguration));
            using (XmlWriter writer = XmlWriter.Create(fileName)) {
                serializer.WriteObject(writer, config);
            }
        }

        /// <summary>
        /// Uses a <see cref="DataContractSerializer"/> to de-serialize the
        /// configuration stored in <paramref name="fileName"/>
        /// </summary>
        /// <param name="fileName">
        /// The name of the file containing the serialized configuration
        /// </param>
        /// <returns>
        /// The deserialized configuration or null if the deserialization
        /// failed
        /// </returns>
        public static FieldStateConfiguration Deserialize(string fileName) {
            DataContractSerializer serializer =
                new DataContractSerializer(typeof(FieldStateConfiguration));
            using (XmlReader reader = XmlReader.Create(fileName)) {
                return (FieldStateConfiguration)serializer.ReadObject(reader);
            }
        }

        /// <summary>
        /// Uses a <see cref="XmlSerializer"/> to de-serialize the content of a
        /// file with the name <paramref name="fileName"/>.
        /// </summary>
        /// <param name="fileName"></param>
        public static FieldStateConfiguration LoadConfig(string fileName) {
            if (File.Exists(fileName)) {
                return FieldStateConfiguration.Deserialize(fileName);
            } else {
                throw new ApplicationException("Could not find \"" + fileName + "\"");
            }
        }

        /// <summary>
        /// Basically compares this object to <paramref name="other"/> by
        /// comparing all public attributes for equality
        /// </summary>
        /// <param name="other">An object to be compared to</param>
        /// <returns>
        /// True, if all public attributes are equal; false otherwise
        /// </returns>
        public bool Equals(FieldStateConfiguration other) {
            int i;
            if (this.BasePaths.GetLength(0) != other.BasePaths.GetLength(0)) {
                return false;
            }

            for (i = 0; i < this.BasePaths.GetLength(0); i++) {
                if (this.BasePaths[i].Equals(other.BasePaths[i]) == false) {
                    return false;
                }
            }

            if (this.SessionGuid.Equals(other.SessionGuid) == false) {
                return false;
            }

            if (this.FieldNames.Count != other.FieldNames.Count) {
                return false;
            }

            for (i = 0; i < this.FieldNames.Count; i++) {
                if (this.FieldNames[i].Equals(other.FieldNames[i]) == false) {
                    return false;
                }
            }

            if (this.TimeSteps.Count != other.TimeSteps.Count) {
                return false;
            }

            for (i = 0; i < this.TimeSteps.Count; i++) {
                if (this.TimeSteps[i].Equals(other.TimeSteps[i]) == false) {
                    return false;
                }
            }

            if (this.ExportFormat.Equals(other.ExportFormat) == false) {
                return false;
            }

            if (this.SuperSampling != other.SuperSampling) {
                return false;
            }

            if (this.ReconstructionType.Equals(other.ReconstructionType) == false) {
                return false;
            }

            if (this.GhostLevel.Equals(other.GhostLevel) == false) {
                return false;
            }

            if (this.NumberOfProcesses.Equals(other.NumberOfProcesses) == false) {
                return false;
            }

            return true;
        }

        /// <summary>
        /// Clears all data in public attributes
        /// </summary>
        /// <param name="defaultReconstructionType">
        /// The initially selected reconstruction type
        /// </param>
        /// <param name="defaultExportFormat">
        /// The initially selected export format
        /// </param>
        /// <param name="defaultGhostLevel">
        /// </param>
        public void Clear(ReconstructionTypes defaultReconstructionType, GhostLevels defaultGhostLevel, ExportFormats defaultExportFormat) {
            FieldNames.Clear();
            TimeSteps.Clear();
            ReconstructionType = defaultReconstructionType;
            GhostLevel = defaultGhostLevel;
            ExportFormat = defaultExportFormat;
            BasePaths = null;
            SessionGuid = Guid.Empty;
            SuperSampling = 0;
            NumberOfProcesses = 1;
        }
    }
}
