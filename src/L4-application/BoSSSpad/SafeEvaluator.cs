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

using Mono.CSharp;
using System;
using System.Reflection;
using System.Threading.Tasks;
using System.Linq;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Wrapper around <see cref="Mono.CSharp.Evaluator"/> which allows to
    /// handle crashes gracefully
    /// </summary>
    public class SafeEvaluator {

        /// <summary>
        /// Function to create a new evaluator
        /// </summary>
        private readonly Func<Evaluator> evaluatorFactory;

        /// <summary>
        /// The current evaluator (which is used until its broken)
        /// </summary>
        private Evaluator evaluator;

        /// <summary>
        /// Constructs a new evaluator using
        /// <paramref name="evaluatorFactory"/>
        /// </summary>
        /// <param name="evaluatorFactory"></param>
        public SafeEvaluator(Func<Evaluator> evaluatorFactory) {
            this.evaluatorFactory = evaluatorFactory;
            this.evaluator = evaluatorFactory();
        }

        /// <summary>
        /// Returns the assembly produced by the latest call to <see cref="Evaluator.Evaluate(string, out object, out bool)"/>.
        /// </summary>
        public Assembly LatestAssembly {
            get {
                // dirty hack using Reflection in order to access private members of the evaluator
                var members = evaluator.GetType().GetFields(BindingFlags.NonPublic | BindingFlags.FlattenHierarchy | BindingFlags.Static | BindingFlags.Instance);
                FieldInfo module_fi =  members.Single(member => member.Name.Contains("module"));
                ModuleContainer mc = (ModuleContainer) module_fi.GetValue(evaluator);
                return mc.Builder.Assembly;
            }
        }

        /// <summary>
        /// Direct access to the current evaluator. Can be used for unsafe
        /// evaluation
        /// </summary>
        public Evaluator Evaluator {
            get {
                // Makes sure evaluator is run again when old evaluator has
                // been destroyed
                if (evaluator == null) {
                    evaluator = evaluatorFactory();
                }
                return evaluator;
            }
        }

        /// <summary>
        /// Safe version of
        /// <see cref="Evaluator.Compile(string)"/>
        /// </summary>
        /// <param name="input"></param>
        /// <param name="result"></param>
        /// <param name="timeout"></param>
        /// <returns></returns>
        public bool TryCompile(string input, out CompiledMethod result, int timeout) {
            return SafeEvaluate(
                () => Evaluator.Compile(input),
                timeout,
                out result);
        }

        /// <summary>
        /// Safe version of
        /// <see cref="Evaluator.Evaluate(string)"/>
        /// </summary>
        /// <param name="input"></param>
        /// <param name="result"></param>
        /// <param name="timeout"></param>
        /// <returns></returns>
        public bool TryEvaluate(string input, object result, int timeout) {
            return SafeEvaluate(
                () => Evaluator.Evaluate(input),
                timeout,
                out result);
        }

        /// <summary>
        /// Safe version of
        /// <see cref="Evaluator.Evaluate(string, out object, out bool)"/>
        /// </summary>
        /// <param name="input"></param>
        /// <param name="result"></param>
        /// <param name="result_set"></param>
        /// <param name="timeout"></param>
        /// <returns></returns>
        public bool TryEvaluate(string input, out object result, out bool result_set, int timeout) {
            string dummy;
            object myResult = null;
            bool myResultSet = false;
            bool completed = SafeEvaluate(
                () => Evaluator.Evaluate(input, out myResult, out myResultSet),
                timeout,
                out dummy);
            result = myResult;
            result_set = myResultSet;
            return completed;
        }

        /// <summary>
        /// Safe version of
        /// <see cref="Evaluator.GetCompletions(string, out string)"/>
        /// </summary>
        /// <param name="input"></param>
        /// <param name="completions"></param>
        /// <param name="prefix"></param>
        /// <param name="timeout"></param>
        /// <returns></returns>
        public bool TryGetCompletions(string input, out string[] completions, out string prefix, int timeout) {
            // only complete current line -- much more doesn't work anyway -- to prevent crashes; sometimes works even better
            string[] inputSplit = input.Split(new string[] { "\n", "\r" }, StringSplitOptions.RemoveEmptyEntries);
            if (inputSplit.Length > 0) {
                input = inputSplit[inputSplit.Length - 1];
            } else {
                completions = new string[0];
                prefix = "";
                return true;
            }

            string myPrefix = null;
            bool completed = SafeEvaluate(
                () => Evaluator.GetCompletions(input, out myPrefix),
                timeout,
                out completions);
            prefix = myPrefix;

            return completed;
        }

        /// <summary>
        /// Safe version of <see cref="Evaluator.Run(string)"/>
        /// </summary>
        /// <param name="statement"></param>
        /// <param name="timeout"></param>
        /// <param name="result"></param>
        /// <returns></returns>
        public bool TryRun(string statement, int timeout, out bool result) {
            return SafeEvaluate(
                () => Evaluator.Run(statement),
                timeout,
                out result);
        }

        /// <summary>
        /// Safely evaluates <paramref name="action"/> by wrapping a task
        /// around it which is canceled after if it has not completed after a
        /// given <paramref name="timeout"/>.
        /// </summary>
        /// <typeparam name="T"></typeparam>
        /// <param name="action"></param>
        /// <param name="timeout"></param>
        /// <param name="result"></param>
        /// <returns></returns>
        private bool SafeEvaluate<T>(Func<T> action, int timeout, out T result) {
            var task = Task.Factory.StartNew(action);

            if (Task.WaitAny(new[] { task }, TimeSpan.FromMilliseconds(timeout)) < 0) {
                // timeout
                Evaluator.Interrupt();
                result = default(T);

                // Get rid of old evaluator; it will not respond to anything
                evaluator = null;
                return false;
            } else if (task.Exception != null) {
                // Exception
                Evaluator.Interrupt();
                throw task.Exception;
            }

            result = task.Result;
            return true;
        }
    }
}
