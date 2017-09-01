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

using ilPSP;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Linq;
using System.Reflection;
using System.Text.RegularExpressions;
using System.Xml.Linq;
using System.Xml.XPath;
using BoSSS.Platform;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Extension methods for <see cref="MemberInfo"/> and
    /// <see cref="ParameterInfo"/> concerning the extraction of XML comments
    /// at runtime.
    /// </summary>
    /// <remarks>
    /// Adapted from 'Reading XML Documentation at Run-Time' by Bradly Smith
    /// (see http://www.brad-smith.info/blog/archives/220)
    /// </remarks>
    public static class XmlDocumentationExtensions {

        /// <summary>
        /// Retrieves the XML documentation for the given
        /// <paramref name="member"/> from the corresponding XML file.
        /// </summary>
        /// <param name="member">
        /// The <see cref="MemberInfo"/> object in question
        /// </param>
        /// <returns>
        /// A plain string containing the textual content of the XML comments
        /// for <paramref name="member"/>.
        /// </returns>
        public static string GetDocumentation(this MemberInfo member) {
            AssemblyName assemblyName = member.DeclaringType.Module.Assembly.GetName();
            XDocument xml = XDocument.Load(assemblyName.Name + ".xml");

            XElement element = xml.XPathEvaluate(String.Format(
                    "/doc/members/member[@name='{0}']/summary",
                    GetMemberElementName(member))).
                Cast<IEnumerable>().Cast<XElement>().FirstOrDefault(); // Yikes!
            if (element == null) {
                return "";
            }
            
            return element.DescendantNodes().Select(n => n.ToHumanReadableString()).CleanAggregate();
        }

        /// <summary>
        /// Retrieves the XML documentation for the given
        /// <paramref name="parameter"/> from the corresponding XML file.
        /// </summary>
        /// <param name="parameter">
        /// The <see cref="ParameterInfo"/> object in question
        /// </param>
        /// <returns>
        /// A plain string containing the textual content of the XML comments
        /// for <paramref name="parameter"/>.
        /// </returns>
        public static string GetDocumentation(this ParameterInfo parameter) {
            AssemblyName assemblyName = parameter.Member.Module.Assembly.GetName();
            XDocument xml = XDocument.Load(assemblyName.Name + ".xml");

            object result;
            if (parameter.IsRetval || String.IsNullOrEmpty(parameter.Name)) {
                result = xml.XPathEvaluate(String.Format(
                    "/doc/members/member[@name='{0}']/returns",
                    GetMemberElementName(parameter.Member)));
            } else {
                result = xml.XPathEvaluate(String.Format(
                    "/doc/members/member[@name='{0}']/param[@name='{1}']",
                    GetMemberElementName(parameter.Member),
                    parameter.Name));
            }

            XElement element = result.Cast<IEnumerable>().Cast<XElement>().FirstOrDefault(); // Yikes!

            if (element == null) {
                return "";
            }
            return element.DescendantNodes().Select(n => n.ToHumanReadableString()).CleanAggregate();
        }

        /// <summary>
        /// Joins all strings in <paramref name="dirtyStrings"/> while removing
        /// all sorts of unnecessary whitespace.
        /// </summary>
        /// <param name="dirtyStrings"></param>
        /// <returns></returns>
        private static string CleanAggregate(this IEnumerable<string> dirtyStrings) {
            if (dirtyStrings.IsNullOrEmpty()) {
                return "";
            }

            return dirtyStrings.Aggregate((s, t) => s + " " + t).WithoutExcessiveWhiteSpace().Trim();
        }

        /// <summary>
        /// Replaces each multi-character whitespace by a single space.
        /// </summary>
        /// <param name="originalString"></param>
        /// <returns></returns>
        private static string WithoutExcessiveWhiteSpace(this string originalString) {
            return Regex.Replace(originalString, @"\s+", " ");
        }

        /// <summary>
        /// Creates a human-readable version of a
        /// <paramref name="documentationNode"/>. Basically, this means that
        /// textual nodes won't change, while tags like 'paramref' and 'see'
        /// will be transformed into something that humans can understand.
        /// </summary>
        /// <param name="documentationNode">
        /// The node to be transformed.
        /// </param>
        /// <returns>
        /// A plain string representation of <paramref name="documentationNode"/>
        /// </returns>
        private static string ToHumanReadableString(this XNode documentationNode) {
            if (documentationNode is XElement) {
                XElement xelement = (XElement)documentationNode;
                switch (xelement.Name.ToString()) {
                    case "paramref":
                        return xelement.Attribute("name").Value;

                    case "see":
                        return xelement.Attribute("cref").Value;

                    default:
                        //Unknown element, ignore!
                        return "";
                }
            } else {
                return documentationNode.ToString();
            }
        }

        /// <summary>
        /// Returns the expected name for a member element in the XML
        /// documentation file. The general format is, using the example of a
        /// method,'M:Namespace.Class.Method'.
        /// </summary>
        /// <param name="member">The reflected member.</param>
        /// <returns>The name of the member element.</returns>
        private static string GetMemberElementName(MemberInfo member) {
            string memberName;
            if (member is Type) {
                memberName = member.As<Type>().FullName;
            } else {
                // member belongs to a Type
                memberName = member.DeclaringType.FullName + "." + member.Name;
            }

            char prefixCode;
            switch (member.MemberType) {
                case MemberTypes.Constructor:
                    // XML documentation uses slightly different constructor names
                    memberName = memberName.Replace(".ctor", "#ctor");
                    goto case MemberTypes.Method;
                case MemberTypes.Method:
                    prefixCode = 'M';

                    // parameters are listed according to their type, not their name
                    string paramTypesList = String.Join(
                        ",",
                        ((MethodBase)member).GetParameters()
                            .Cast<ParameterInfo>()
                            .Select(x => x.ParameterType.GetGenericTypeName(true)
                        ).ToArray()
                    );
                    if (!String.IsNullOrEmpty(paramTypesList))
                        memberName += "(" + paramTypesList + ")";
                    break;

                case MemberTypes.Event:
                    prefixCode = 'E';
                    break;

                case MemberTypes.Field:
                    prefixCode = 'F';
                    break;

                case MemberTypes.NestedType:
                    // XML documentation uses slightly different nested type names
                    memberName = memberName.Replace('+', '.');
                    goto case MemberTypes.TypeInfo;
                case MemberTypes.TypeInfo:
                    prefixCode = 'T';
                    break;

                case MemberTypes.Property:
                    prefixCode = 'P';
                    break;

                default:
                    throw new ArgumentException("Unknown member type", "member");
            }

            return String.Format("{0}:{1}", prefixCode, memberName);
        }
    }
}