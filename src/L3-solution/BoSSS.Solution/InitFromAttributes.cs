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
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Foundation.XDG;
using BoSSS.Solution.Control;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Reflection;
using System.Text.RegularExpressions;

namespace BoSSS.Solution {

    /// <summary>
    /// a totally useless thing
    /// </summary>
    public static class InitFromAttributes {

        /// <summary>
        /// Instantiates those fields that are marked by an <see cref="InstantiateFromControlFileAttribute"/> - attribute;
        /// </summary>
        public static void CreateFieldsAuto(object _this, IGridData context,
            IDictionary<string, FieldOpts> FieldOptions, XQuadFactoryHelper.MomentFittingVariants cutCellQuadType,
            ICollection<DGField> IOFields, ICollection<DGField> RegisteredFields) {
            FieldInfo[] fields = _this.GetType().GetFields(BindingFlags.Instance | BindingFlags.NonPublic | BindingFlags.Public | BindingFlags.FlattenHierarchy | BindingFlags.Static);

            FieldInfo LevelSetTracker = null;

            List<FieldInfo> LevelSets = new List<FieldInfo>();
            List<FieldInfo> OtherFields = new List<FieldInfo>();

            #region DetermineTypeOfFields
            // Determine the Type of the Fields
            foreach (FieldInfo f in fields) {
                object[] atts = f.GetCustomAttributes(typeof(InstantiateFromControlFileAttribute), true);

                if (atts.Length != 0) {
                    if (atts.Length != 1)
                        throw new ApplicationException("should not happen");

                    Type HistoryType, VectorType, ComponentType;
                    GetTypes(f, out HistoryType, out VectorType, out ComponentType);

                    if (typeof(LevelSet).IsAssignableFrom(ComponentType))
                        LevelSets.Add(f);
                    else
                        OtherFields.Add(f);
                }

                atts = f.GetCustomAttributes(typeof(LevelSetTrackerAttribute), true);

                if (atts.Length != 0) {
                    if (atts.Length != 1)
                        throw new ApplicationException("should not happen");

                    if (f.FieldType != typeof(LevelSetTracker))
                        throw new ApplicationException("stupid.");

                    if (LevelSetTracker != null) {
                        throw new ApplicationException("only one level set tracker is allowed");
                    }
                    LevelSetTracker = f;

                }

            }
            #endregion

            #region CreateLevelSets
            // create Level Sets ..
            // =====================

            // must be initialized before the XDG fields ...

            FieldInfo[] LevelSetsSorted = new FieldInfo[LevelSets.Count];

            if (LevelSets.Count > 4)
                throw new ApplicationException("at maximum, 4 different level sets are supported.");

            if (LevelSets.Count > 1) {
                // more than 1 level set -> level set index is required.


                for (int i = 0; i < LevelSets.Count; i++) {
                    FieldInfo fi = LevelSets[i];
                    object[] atts2 = fi.GetCustomAttributes(typeof(LevelSetIndexAttribute), true);
                    if (atts2.Length > 1)
                        throw new ApplicationException("should not happen (AttributeUsage.AllowMultiple = false for 'LevelSetIndexAttribute')");
                    if (atts2.Length < 1)
                        throw new ApplicationException("missing 'LevelSetIndexAttribute' for level set '"
                            + fi.Name + "' in class '" + fi.ReflectedType.Name + "'; this contains " + LevelSets.Count
                            + " level sets, so level set indices in the range 0 (including) to " + LevelSets.Count
                            + " (excluding) are required.");

                    LevelSetIndexAttribute lsia = (LevelSetIndexAttribute)atts2[0];

                    if (lsia.m_idx < 0 || lsia.m_idx >= LevelSets.Count) {
                        throw new ApplicationException("error with level set index for level set '"
                            + fi.Name + "': index is specified as " + lsia.m_idx + ", but allowed range is "
                            + "0 (including) to " + LevelSets.Count + " (excluding).");
                    }

                    if (LevelSetsSorted[lsia.m_idx] != null) {
                        throw new ApplicationException("error with level set index " + lsia.m_idx + ": "
                            + " assigned for at least two level sets, i.e. for '" + fi.Name
                            + "' and for '" + LevelSetsSorted[lsia.m_idx].Name + "' -- each index is allowed only once.");
                    }

                }

            } else if (LevelSets.Count == 1) {
                LevelSetsSorted[0] = LevelSets[0];
            }


            LevelSet[] LevelSetsInstances = new LevelSet[LevelSets.Count];
            for (int i = 0; i < LevelSetsSorted.Length; i++) {
                InstantiateField(LevelSetsSorted[i], _this, IOFields, RegisteredFields, FieldOptions, context, null);

                object o = LevelSetsSorted[i].GetValue(_this);

                if (o is ScalarFieldHistory<LevelSet>)
                    LevelSetsInstances[i] = ((ScalarFieldHistory<LevelSet>)LevelSetsSorted[i].GetValue(_this)).Current;
                else
                    LevelSetsInstances[i] = (LevelSet)LevelSetsSorted[i].GetValue(_this);
            }

            #endregion

            #region CreareLevelSetTracker
            // create level set tracker
            // ========================


            LevelSetTracker lsTrk = null;

            foreach (DGField f in RegisteredFields) {
                if (f is XDGField) {
                    LevelSetTracker _lsTrk = ((XDGField)f).Basis.Tracker;

                    if (lsTrk == null) {
                        lsTrk = _lsTrk;
                    } else {
                        if (!object.ReferenceEquals(lsTrk, _lsTrk))
                            throw new NotSupportedException();
                    }
                }
            }

            if (lsTrk == null) {

                if (LevelSets.Count > 0) {
                    if (LevelSetTracker == null)
                        throw new ApplicationException("missing Level Set Tracker");
                    LevelSetTrackerAttribute att = (LevelSetTrackerAttribute)(LevelSetTracker.GetCustomAttributes(typeof(LevelSetTrackerAttribute), true)[0]);
                    

                    switch (LevelSetsSorted.Length) {
                        case 1:
                            lsTrk = new LevelSetTracker(context, cutCellQuadType, att.m_NearCellWidth, (string[])(att.GetSpeciesTable(1)), LevelSetsInstances[0]);
                            break;
                        case 2:
                            lsTrk = new LevelSetTracker(context, cutCellQuadType, att.m_NearCellWidth, (string[,])(att.GetSpeciesTable(2)), LevelSetsInstances[0], LevelSetsInstances[1]);
                            break;
                        case 3:
                            lsTrk = new LevelSetTracker(context, cutCellQuadType, att.m_NearCellWidth, (string[,,])(att.GetSpeciesTable(3)), LevelSetsInstances[0], LevelSetsInstances[1], LevelSetsInstances[2]);
                            break;
                        case 4:
                            lsTrk = new LevelSetTracker(context, cutCellQuadType, att.m_NearCellWidth, (string[,,,])(att.GetSpeciesTable(4)), LevelSetsInstances[0], LevelSetsInstances[1], LevelSetsInstances[2], LevelSetsInstances[3]);
                            break;

                    }

                    LevelSetTracker.SetValue(_this, lsTrk);
                }
            }

            #endregion


            // create other Fields, e.g. single phase, XDG, ...
            // ================================================

            foreach (FieldInfo f in OtherFields) {
                InstantiateField(f, _this, IOFields, RegisteredFields, FieldOptions, context, lsTrk);
            }
        }


        static void InstantiateField(FieldInfo f, object o, ICollection<DGField> IOFields, ICollection<DGField> RegisteredFields,
            //AppControl ctrl, 
            IDictionary<string, FieldOpts> FieldOptions,
            GridData ctx, LevelSetTracker lstrk) {
            // get attribute
            // =============
            InstantiateFromControlFileAttribute at;
            {
                object[] atts = f.GetCustomAttributes(typeof(InstantiateFromControlFileAttribute), true);
                at = (InstantiateFromControlFileAttribute)atts[0];
            }

            // instantiate
            // ===========
            var member_value = InstantiateFromAttribute(f, at, f.FieldType, IOFields, RegisteredFields, FieldOptions, ctx, lstrk);


            // set value
            // =========
            f.SetValue(o, member_value);
        }

        /// <summary>
        /// still a hack...
        /// </summary>
        static object InstantiateFromAttribute(FieldInfo f, InstantiateFromControlFileAttribute at, Type type,
            ICollection<DGField> IOFields, ICollection<DGField> RegisteredFields,
            //AppControl ctrl,
            IDictionary<string, FieldOpts> FieldOptions,
            GridData ctx, LevelSetTracker lstrk) {

            // create instance
            // ===============

            object member_value = null;

            Type HistoryType = null;
            Type VectorType = null;
            Type ComponentType = null;

            GetTypes(f, out HistoryType, out VectorType, out ComponentType);


            if (VectorType != null) {

                // vector field branch
                // +++++++++++++++++++

                if (at.m_IsScalarField)
                    throw new ApplicationException("illegal use of 'InstantiateFromControlFileAttribute' (with Scalar Declaration) on a Vector class");


                int D = ctx.SpatialDimension;
                string[] cName = at.GetControlFileNames(f, D);
                string[] iName = at.GetInCodeIdentifications(f, D);


                // determine DG polynomial degree of basis
                int[] Deg = new int[D];
                for (int d = 0; d < D; d++)
                    Deg[d] = GetDegree(cName[d], iName[d], FieldOptions);

                if (at.m_DegreesMustBeEqual) {
                    int deg0 = Deg[0];
                    for (int d = 1; d < D; d++) {
                        if (Deg[d] != deg0) {
                            StringWriter errMsg = new StringWriter();
                            errMsg.Write("DG Polynomial degree of fields {");
                            for (int dd = 0; dd < D; dd++) {
                                errMsg.Write(cName[dd]);
                                if (dd < D - 1)
                                    errMsg.Write(", ");
                            }
                            errMsg.Write("} must be equal, but found {");
                            for (int dd = 0; dd < D; dd++) {
                                errMsg.Write(Deg[dd]);
                                if (dd < D - 1)
                                    errMsg.Write(", ");
                            }
                            errMsg.Write("} in control file.");


                            throw new ApplicationException(errMsg.ToString());
                        }
                    }
                }

                // create instance: components 
                SinglePhaseField[] fld = new SinglePhaseField[D];
                XDGField[] xfld = new XDGField[D];
                DGField[] _fld = new DGField[D];
                for (int d = 0; d < D; d++) {

                    if (ComponentType == typeof(SinglePhaseField)) {
                        fld[d] = new SinglePhaseField(new Basis(ctx, Deg[d]), iName[d]);
                        _fld[d] = fld[d];

                    } else if (ComponentType == typeof(XDGField)) {
                        xfld[d] = new XDGField(new XDGBasis(lstrk, Deg[d]), iName[d]);
                        _fld[d] = xfld[d];
                        fld = null;
                    } else {
                        throw new NotSupportedException("unknown type.");
                    }
                    RegisteredFields.Add(_fld[d]);
                }

                // create instance: Vector-Field container
                var ci = VectorType.GetConstructor(new Type[] { ComponentType.MakeArrayType() });
                member_value = ci.Invoke(new object[] { (fld != null) ? ((object)fld) : ((object)xfld) });
                //member_value = ci.Invoke( new object[] { ((object)fld)  });

                // io
                for (int d = 0; d < D; d++)
                    AddToIO(iName[d], _fld[d], FieldOptions, IOFields, at);

            } else {
                // scalar field branch
                // +++++++++++++++++++

                if (at.m_IsVectorField)
                    throw new ApplicationException("illegal use of 'InstantiateFromControlFileAttribute' (with Vector Declaration) on a Non-Vector class");

                // identification 
                string cName = at.GetControlFileName(f);
                string iName = at.GetInCodeIdentification(f);

                // create basis
                int Deg = GetDegree(cName, iName, FieldOptions);
                Basis b = new Basis(ctx, Deg);

                // create instance
                DGField fld = null;
                if (ComponentType == typeof(SinglePhaseField))
                    fld = new SinglePhaseField(b, iName);
                else if (ComponentType == typeof(LevelSet))
                    fld = new LevelSet(b, iName);
                else if (ComponentType == typeof(XDGField))
                    fld = new XDGField(new XDGBasis(lstrk, Deg), iName);
                else
                    throw new NotImplementedException();

                RegisteredFields.Add(fld);
                AddToIO(iName, fld, FieldOptions, IOFields, at);

                member_value = fld;
            }

            // History, if desired
            if (HistoryType != null) {
                member_value = (HistoryType.GetConstructors()[0]).Invoke(new object[] { member_value });
            }

            return member_value;
        }

        private static void GetTypes(FieldInfo f, out Type HistoryType, out Type VectorType, out Type ComponentType) {
            HistoryType = null;
            VectorType = null;
            ComponentType = null;

            Type t = f.FieldType;

            if (t.IsGenericType && t.GetGenericTypeDefinition() == typeof(VectorField<DGField>).GetGenericTypeDefinition()) {
                // vector field init
                VectorType = t;
                ComponentType = t.GetGenericArguments()[0];
            } else if (t.IsGenericType && t.GetGenericTypeDefinition() == typeof(VectorFieldHistory<DGField>).GetGenericTypeDefinition()) {
                // vector field history
                HistoryType = t;
                ComponentType = t.GetGenericArguments()[0];
                VectorType = (typeof(VectorField<DGField>)).GetGenericTypeDefinition().MakeGenericType(ComponentType);
            } else if (t.IsGenericType && t.GetGenericTypeDefinition() == typeof(ScalarFieldHistory<DGField>).GetGenericTypeDefinition()) {
                // scalar field history
                HistoryType = t;
                ComponentType = t.GetGenericArguments()[0];
            } else {
                // scalar field
                ComponentType = t;
            }

            if (!typeof(DGField).IsAssignableFrom(ComponentType)) {
                throw new ApplicationException("illegal use of '" + typeof(InstantiateFromControlFileAttribute).Name + "' on the field "
                    + " '" + f.Name + "':  got type '" + ComponentType.Name + "', but expecting some subclass of '" + typeof(DGField).Name + "'.");
            }

            if (ComponentType.IsAbstract || ComponentType.IsGenericType) {
                throw new ApplicationException("illegal use of '" + typeof(InstantiateFromControlFileAttribute).Name + "' on the field "
                    + " '" + f.Name + "': the type is abstract or generic ('" + ComponentType.Name + "'), therefore it is not possible "
                    + "to determine what kind of object should be created.");
            }
        }

        /// <summary>
        /// used by <see cref="CreateFieldsAuto"/>
        /// </summary>
        static private void AddToIO(string iName, DGField fld,
            //AppControl ctrl,
            IDictionary<string, FieldOpts> FieldOptions,
            ICollection<DGField> IOFields, InstantiateFromControlFileAttribute at) {

            FieldOpts fopts = FieldOptions.Where(kv => WildcardToRegex(kv.Key).IsMatch(fld.Identification)).SingleOrDefault().Value;

            if (fopts != null) {
                if (at.m_ioListOpt == IOListOption.Always && fopts.SaveToDB == FieldOpts.SaveToDBOpt.FALSE)
                    throw new ApplicationException("IO for field '" + fld.Identification + "' cannot be turned OFF, i.e. 'SaveToDB==false' is illegal.");
                if (at.m_ioListOpt == IOListOption.Never && fopts.SaveToDB == FieldOpts.SaveToDBOpt.TRUE)
                    throw new ApplicationException("IO for field '" + fld.Identification + "' cannot be turned ON, i.e. 'SaveToDB==true' is illegal");

                if (at.m_ioListOpt == IOListOption.Always || fopts.SaveToDB == FieldOpts.SaveToDBOpt.TRUE)
                    IOFields.Add(fld);

            } else {
                if (at.m_ioListOpt == IOListOption.Always)
                    IOFields.Add(fld);
            }

        }

        private static Regex WildcardToRegex(string pattern) {
            return new Regex("^" + Regex.Escape(pattern).
            Replace("\\*", ".*").
            Replace("\\?", ".") + "$");
        }



        /// <summary>
        /// used by <see cref="CreateFieldsAuto"/>
        /// </summary>
        static private int GetDegree(string cName, string iName,
            //AppControl ctrl
            IDictionary<string, FieldOpts> FieldOptions) //
        {
            
           var ops = FieldOptions.Where(kv => WildcardToRegex(kv.Key).IsMatch(cName));

            int Deg = -1;
            {
                if (cName == iName) {


                    if (ops.Count() != 1)
                        throw new ApplicationException("missing DG polynomial degree specification for field '" + cName + "' in control file;");

                    Deg = (int)(ops.First().Value.Degree);

                    if (Deg < 0)
                        throw new ApplicationException("missing DG polynomial degree specification for field '" + cName + "' in control file;");
                } else {
                    // polynomial degree of field 'iName' is implied by polynomial degree of field 'cName' in ctrl file

                    if (ops.Count() != 1)
                        throw new ApplicationException("missing DG polynomial degree specification for field '" + cName +
                            "', which implies DG polynomial degree of field '" + iName + "', in control file;");

                    Deg = (int)(ops.First().Value.Degree);

                    if (Deg < 0)
                        throw new ApplicationException("missing DG polynomial degree specification for field '" + cName +
                            "', which implies DG polynomial degree of field '" + iName + "', in control file;");

                    if (FieldOptions.ContainsKey(iName) && (FieldOptions[iName].Degree >= 0 && FieldOptions[iName].Degree != Deg)) {
                        throw new ApplicationException("Control file error: specification of DG polynomial degree for field '"
                            + iName + "' is implied, by definition, from field '" + cName + "', therefore its degree may not be specified differently.");
                    }
                }
            }
            return Deg;
        }
    }
}
