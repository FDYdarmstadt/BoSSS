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
using System.Reflection;
using System.Text.RegularExpressions;
using ilPSP.Utils;
using BoSSS.Solution.Control;

namespace BoSSS.Solution {

    /// <summary>
    /// If a member in a subclass of <see cref="Application"/>
    /// is decorated by this attribute, which is only allowed for members of certain types
    /// (in particular <see cref="BoSSS.Foundation.SinglePhaseField"/>, <see cref="BoSSS.Foundation.XDG.LevelSet"/>, <see cref="BoSSS.Foundation.XDG.XDGField"/>),
    /// it will be initialized automatically during the construction of the object.
    /// </summary>
    [AttributeUsage(AttributeTargets.Field, AllowMultiple = false, Inherited = false)]
    public class InstantiateFromControlFileAttribute : Attribute {


        /// <summary>
        /// ctor
        /// </summary>
        /// <param name="FieldIdentification">
        /// Specifies the name that the created DG field will have after instantiation (member <see cref="BoSSS.Foundation.DGField.Identification"/>);
        /// If null or empty, the identification of the field will be set equal to the member name.
        /// </param>
        /// <param name="DGDegreeImpliedBy">
        /// This can be used to ensure that two or more fields in the code have equal DG polynomial Degree:
        /// if not null, the DG polynomial degree of the field is determined by another field with this name;
        /// E.g. if this argument is set to 'Psi', then the DG degree will be determined
        /// by the control file entry for field 'Psi', independent of what is specified for <paramref name="FieldIdentification"/>.
        /// </param>
        /// <param name="ioListOpt"></param>
        public InstantiateFromControlFileAttribute(string FieldIdentification = null, string DGDegreeImpliedBy = null, IOListOption ioListOpt = IOListOption.ControlFileDetermined) {
            m_ControlFileName = DGDegreeImpliedBy;
            if (m_ControlFileName != null && m_ControlFileName == "") m_ControlFileName = null;

            m_FieldIdentification = FieldIdentification;
            if (m_FieldIdentification != null && m_FieldIdentification == "") m_FieldIdentification = null;

            m_ioListOpt = ioListOpt;
            m_IsScalarField = true;            
        }


        string m_ControlFileName = null;
        string m_FieldIdentification = null;
        
        internal IOListOption m_ioListOpt = IOListOption.ControlFileDetermined;


        internal string GetControlFileName(FieldInfo fi) {
            string cName = (this.m_ControlFileName == null || this.m_ControlFileName == "") ? GetInCodeIdentification(fi) : this.m_ControlFileName;
            return cName;
        }

        internal string GetInCodeIdentification(FieldInfo fi) {
            string iName = (this.m_FieldIdentification == null || this.m_FieldIdentification == "") ? fi.Name : this.m_FieldIdentification;
            return iName;
        }
        
        internal bool m_IsVectorField = false;
        internal bool m_IsScalarField = false;
        
        /// <summary>
        /// ctor, especially for the creation of vector fields
        /// </summary>
        /// <param name="DGDegreeImpliedBy">
        /// This can be used to ensure that two or more fields in the code have equal DG polynomial Degree:
        /// if not null, the DG polynomial degree of the field is determined by another field with this name;
        /// E.g. if this argument is set to 'Psi', then the DG degree will be determined
        /// by the control file entry for field 'Psi', independent of what is specified for <paramref name="FieldIdentifications"/>.
        /// </param>
        /// <param name="DegreesMustBeEqual">
        /// An exception is thrown if the degrees of all members in <paramref name="DegreesMustBeEqual"/>
        /// are not specified equal in the control file.
        /// </param>
        /// <param name="FieldIdentifications">
        /// Names for the components of the Vector field (see <see cref="BoSSS.Foundation.DGField.Identification"/>);
        /// </param>
        /// <param name="CapDimension">
        /// Only the first <em>D</em> entries of <paramref name="FieldIdentifications"/> and <paramref name="DGDegreeImpliedBy"/>
        /// will be considered, where <em>D</em> denotes the spatial dimension of the loaded grid; e.g.
        /// if <paramref name="FieldIdentifications"/> = {"u", "v", "w" }, but <em>D</em>=2, only {"v", "w" }
        /// will be taken.
        /// </param>
        /// <param name="ioListOpt"></param>
        public InstantiateFromControlFileAttribute(string[] FieldIdentifications, string[] DGDegreeImpliedBy = null, bool DegreesMustBeEqual = true, bool CapDimension = true, IOListOption ioListOpt = IOListOption.ControlFileDetermined) {
            m_ioListOpt = ioListOpt;
            m_DegreesMustBeEqual = DegreesMustBeEqual;
            m_CapDimension = CapDimension;
            m_IsVectorField = true;
            m_ControlFileNames = DGDegreeImpliedBy;
            m_FieldIdentifications = FieldIdentifications;            
        }

        internal bool m_DegreesMustBeEqual;
        bool m_CapDimension;

        string[] m_ControlFileNames = null;
        string[] m_FieldIdentifications = null;

        internal string[] GetControlFileNames(FieldInfo fi, int D) {
            if (m_ControlFileNames != null) {
                if (m_CapDimension) {
                    if (m_ControlFileNames.Length < D)
                        throw new ApplicationException("error initializing vector field '" + fi.Name + "': mismatch in spatial dimension.");
                    return ArrayTools.GetSubVector(m_ControlFileNames, 0, D);
                } else {
                    if (m_ControlFileNames.Length != D)
                        throw new ApplicationException("error initializing vector field '" + fi.Name + "': mismatch in spatial dimension.");
                    return m_ControlFileNames;
                }
            } else {
                return GetInCodeIdentifications(fi, D);
            }
        }

        internal string[] GetInCodeIdentifications(FieldInfo fi, int D) {
            if (m_FieldIdentifications != null) {
                if (m_CapDimension) {
                    if (m_FieldIdentifications.Length < D)
                        throw new ApplicationException("Length of field identification vector is to small; (Field '" + fi.Name + "')");
                    return ArrayTools.GetSubVector(m_FieldIdentifications, 0, D);
                } else {
                    if (m_FieldIdentifications.Length != D)
                        throw new ApplicationException("Mismatch in Length of field identification vector; (Field '" + fi.Name + "')");
                    return m_FieldIdentifications;
                }
            } else {
                throw new ApplicationException("For vector fields, the in-code identification is mandatory; (Field '" + fi.Name + "')");
            }
        }
    }
}


