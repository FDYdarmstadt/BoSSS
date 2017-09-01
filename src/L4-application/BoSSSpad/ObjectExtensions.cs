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
using System.IO;
using System.Linq;
using System.Reflection;
using System.Runtime.CompilerServices;
using System.Text;
using System.Windows.Forms;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Extension methods for <see cref="object"/>.
    /// </summary>
    public static class ObjectExtensions {

        /// <summary>
        /// Gives a detailed summary of the structure of an object.
        /// </summary>
        /// <param name="obj">The object to be described</param>
        /// <returns>
        /// A string representation of the outline of <paramref name="obj"/>.
        /// </returns>
        /// <seealso cref="Mono.CSharp.InteractiveBase.Describe"/>
        public static string Describe(this object obj) {
            return Mono.CSharp.InteractiveBase.Describe(obj);
        }

        /// <summary>
        /// Prints a summary of the public information about an object.
        /// </summary>
        /// <param name="obj">The object in question</param>
        /// <returns>
        /// A string representation of the publicly available information
        /// about an object.
        /// </returns>
        public static string Summary(this object obj) {
            StringBuilder output = new StringBuilder();

            foreach (var field in obj.GetType().GetFields().Where(f => f.IsPublic)) {
                output.AppendLine(field.Name + ": " + field.GetValue(obj));
            }

            foreach (var property in obj.GetType().GetProperties().Where(p => p.CanRead)) {
                if (property.GetIndexParameters().Count() > 0) {
                    // Indexed property, retrieving the value without
                    // specifying an index would cause an exception
                    continue;
                }

                output.AppendLine(property.Name + ": " + property.GetValue(obj, null));
            }

            return output.ToString();
        }

        /// <summary>
        /// Prints all available actions (i.e, methods and extension methods)
        /// for the given <paramref name="obj"/>.
        /// </summary>
        /// <param name="obj">
        /// The object whose actions should be described
        /// </param>
        public static void Actions(this object obj) {
            var methods = GetCallableMethods(obj, false);

            Console.WriteLine("You can invoke the following methods (more actions may exist):");
            foreach (var group in methods.GroupBy(m => m.Name).OrderBy(g => g.First().Name)) {
                foreach (MethodInfo method in group) {
                    IEnumerable<ParameterInfo> parameters = method.GetParameters();

                    if (method.IsDefined(typeof(ExtensionAttribute), false)) {
                        // Skip the first parameter of extension methods since
                        // they are not part of the call signature
                        parameters = parameters.Skip(1);
                    }

                    Console.WriteLine("- " + method.Name + parameters.ToParameterString());
                }
            }
        }

        /// <summary>
        /// Retrieves textual information about all public methods with name
        /// <paramref name="methodName"/> that are defined for
        /// typeof(<paramref name="obj"/>) from the corresponding XML
        /// documentation file.
        /// </summary>
        /// <param name="obj">
        /// The object whose method(s) should be described.
        /// </param>
        /// <param name="methodName">
        /// The name of the method(s) to be described in detail.
        /// </param>
        public static string Describe(this object obj, string methodName) {
            var methods = GetCallableMethods(obj, false).Where(m => m.Name == methodName);

            using (StringWriter s = new StringWriter()) {
                if (methods.Count() == 0) {
                    s.WriteLine("Could not find any method with name '{0}'", methodName);
                    return s.ToString();
                }

                bool moreThanOne = methods.Count() > 1;
                if (moreThanOne) {
                    s.WriteLine("Action '{0}' has {1} overloads: ", methodName, methods.Count());
                }

                int n = 1;
                foreach (MethodInfo method in methods) {
                    IEnumerable<ParameterInfo> parameters = method.GetParameters();

                    if (method.IsDefined(typeof(ExtensionAttribute), false)) {
                        // Skip the first parameter of extension methods since
                        // they are not part of the call signature
                        parameters = parameters.Skip(1);
                    }

                    s.WriteLine();
                    if (moreThanOne) {
                        s.Write(n + ") ");
                    }

                    if (method.ReturnType != typeof(void)) {
                        s.Write(method.ReturnType.GetGenericTypeName() + " ");
                    }
                    s.WriteLine(method.Name + parameters.ToParameterString() + ":");

                    s.WriteLine("   * Summary: " + method.GetDocumentation());

                    foreach (ParameterInfo parameter in parameters) {
                        s.WriteLine("   * " + parameter.Name + ": " + parameter.GetDocumentation());
                    }

                    if (method.ReturnParameter != null) {
                        s.WriteLine("   * Returns: " + method.ReturnParameter.GetDocumentation());
                    }
                    n++;
                }

                return s.ToString();
            }
        }

        /// <summary>
        /// Converts a list of parameters to a textual representation in the
        /// form (TypeOne paramOneName, TypeTwo paramTwoname, ...) 
        /// </summary>
        /// <param name="parameters">
        /// A list of parameters (typically, of a single method).
        /// </param>
        /// <returns>
        /// A string of the form
        /// (TypeOne paramOneName, TypeTwo paramTwoname, ...)
        /// </returns>
        private static string ToParameterString(this IEnumerable<ParameterInfo> parameters) {
            if (parameters.Count() == 0) {
                return "()";
            }

            var paramNamePairs = parameters.Select(
                p => p.ParameterType.GetGenericTypeName() + " " + p.Name);
            return "(" + paramNamePairs.Aggregate((s, t) => s + ", " + t) + ")";
        }

        /// <summary>
        /// Retrieves 
        /// </summary>
        /// <param name="obj"></param>
        /// <param name="includeFrameworkMethods"></param>
        /// <returns></returns>
        private static IEnumerable<MethodInfo> GetCallableMethods(object obj, bool includeFrameworkMethods) {
            Type type = obj.GetType();

            IEnumerable<MethodInfo> methods = type.GetMethods().Where(m => m.IsPublic && !m.IsSpecialName);
            IEnumerable<Assembly> assemblies = AppDomain.CurrentDomain.GetAssemblies();
            if (!includeFrameworkMethods) {
                string[] includedNamespaces = new[] { "ilPSP", "BoSSS" };

                methods = methods.Where(m => m.DeclaringType.Namespace.StartsWithAny(includedNamespaces));
                assemblies = assemblies.Where(a => a.FullName.StartsWithAny(includedNamespaces));
            }

            // Include extension methods
            return methods.Concat(assemblies.SelectMany(a => type.GetExtensionMethods(a)));
        }

        /// <summary>
        /// Saves a textual representation of the given <paramref name="obj"/>
        /// to the clipboard.
        /// </summary>
        /// <param name="obj">
        /// The object to be written to <see cref="Clipboard"/>.
        /// </param>
        public static void ToClipboard(this object obj) {
            Clipboard.SetText(obj.ToString());
        }

        /// <summary>
        /// Determines whether the string <paramref name="s"/> starts with any
        /// of the strings in <paramref name="potentialStartStrings"/>
        /// </summary>
        /// <param name="s">
        /// The string in question
        /// </param>
        /// <param name="potentialStartStrings">
        /// A list of string <paramref name="s"/> might be starting with.
        /// </param>
        /// <returns>
        /// True, if <paramref name="s"/> starts with any of the strings in
        /// <paramref name="potentialStartStrings"/>.
        /// </returns>
        private static bool StartsWithAny(this string s, IEnumerable<string> potentialStartStrings) {
            return potentialStartStrings.Any(p => s.StartsWith(p));
        }
    }
}
