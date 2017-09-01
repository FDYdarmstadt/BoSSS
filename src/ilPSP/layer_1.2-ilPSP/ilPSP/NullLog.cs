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

using log4net;
using log4net.Core;
using System;

namespace ilPSP {

    /// <summary>
    /// A log that logs nothing
    /// </summary>
    public class NullLog : ILog {

        /// <summary>
        /// 
        /// </summary>
        public static NullLog Instance = new NullLog();

        /// <summary>
        /// Does nothing.
        /// </summary>
        private NullLog() {
        }

        #region ILog Members

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="message"></param>
        /// <param name="exception"></param>
        public void Debug(object message, Exception exception) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="message"></param>
        public void Debug(object message) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="provider"></param>
        /// <param name="format"></param>
        /// <param name="args"></param>
        public void DebugFormat(IFormatProvider provider, string format, params object[] args) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="format"></param>
        /// <param name="arg0"></param>
        /// <param name="arg1"></param>
        /// <param name="arg2"></param>
        public void DebugFormat(string format, object arg0, object arg1, object arg2) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="format"></param>
        /// <param name="arg0"></param>
        /// <param name="arg1"></param>
        public void DebugFormat(string format, object arg0, object arg1) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="format"></param>
        /// <param name="arg0"></param>
        public void DebugFormat(string format, object arg0) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="format"></param>
        /// <param name="args"></param>
        public void DebugFormat(string format, params object[] args) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="message"></param>
        /// <param name="exception"></param>
        public void Error(object message, Exception exception) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="message"></param>
        public void Error(object message) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="provider"></param>
        /// <param name="format"></param>
        /// <param name="args"></param>
        public void ErrorFormat(IFormatProvider provider, string format, params object[] args) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="format"></param>
        /// <param name="arg0"></param>
        /// <param name="arg1"></param>
        /// <param name="arg2"></param>
        public void ErrorFormat(string format, object arg0, object arg1, object arg2) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="format"></param>
        /// <param name="arg0"></param>
        /// <param name="arg1"></param>
        public void ErrorFormat(string format, object arg0, object arg1) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="format"></param>
        /// <param name="arg0"></param>
        public void ErrorFormat(string format, object arg0) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="format"></param>
        /// <param name="args"></param>
        public void ErrorFormat(string format, params object[] args) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="message"></param>
        /// <param name="exception"></param>
        public void Fatal(object message, Exception exception) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="message"></param>
        public void Fatal(object message) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="provider"></param>
        /// <param name="format"></param>
        /// <param name="args"></param>
        public void FatalFormat(IFormatProvider provider, string format, params object[] args) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="format"></param>
        /// <param name="arg0"></param>
        /// <param name="arg1"></param>
        /// <param name="arg2"></param>
        public void FatalFormat(string format, object arg0, object arg1, object arg2) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="format"></param>
        /// <param name="arg0"></param>
        /// <param name="arg1"></param>
        public void FatalFormat(string format, object arg0, object arg1) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="format"></param>
        /// <param name="arg0"></param>
        public void FatalFormat(string format, object arg0) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="format"></param>
        /// <param name="args"></param>
        public void FatalFormat(string format, params object[] args) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="message"></param>
        /// <param name="exception"></param>
        public void Info(object message, Exception exception) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="message"></param>
        public void Info(object message) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="provider"></param>
        /// <param name="format"></param>
        /// <param name="args"></param>
        public void InfoFormat(IFormatProvider provider, string format, params object[] args) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="format"></param>
        /// <param name="arg0"></param>
        /// <param name="arg1"></param>
        /// <param name="arg2"></param>
        public void InfoFormat(string format, object arg0, object arg1, object arg2) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="format"></param>
        /// <param name="arg0"></param>
        /// <param name="arg1"></param>
        public void InfoFormat(string format, object arg0, object arg1) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="format"></param>
        /// <param name="arg0"></param>
        public void InfoFormat(string format, object arg0) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="format"></param>
        /// <param name="args"></param>
        public void InfoFormat(string format, params object[] args) {
        }

        /// <summary>
        /// Returns true
        /// </summary>
        public bool IsDebugEnabled {
            get {
                return true;
            }
        }

        /// <summary>
        /// Returns true
        /// </summary>
        public bool IsErrorEnabled {
            get {
                return true;
            }
        }

        /// <summary>
        /// Returns true
        /// </summary>
        public bool IsFatalEnabled {
            get {
                return true;
            }
        }

        /// <summary>
        /// Returns true
        /// </summary>
        public bool IsInfoEnabled {
            get {
                return true;
            }
        }

        /// <summary>
        /// Returns true
        /// </summary>
        public bool IsWarnEnabled {
            get {
                return true;
            }
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="message"></param>
        /// <param name="exception"></param>
        public void Warn(object message, Exception exception) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="message"></param>
        public void Warn(object message) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="provider"></param>
        /// <param name="format"></param>
        /// <param name="args"></param>
        public void WarnFormat(IFormatProvider provider, string format, params object[] args) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="format"></param>
        /// <param name="arg0"></param>
        /// <param name="arg1"></param>
        /// <param name="arg2"></param>
        public void WarnFormat(string format, object arg0, object arg1, object arg2) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="format"></param>
        /// <param name="arg0"></param>
        /// <param name="arg1"></param>
        public void WarnFormat(string format, object arg0, object arg1) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="format"></param>
        /// <param name="arg0"></param>
        public void WarnFormat(string format, object arg0) {
        }

        /// <summary>
        /// Does nothing.
        /// </summary>
        /// <param name="format"></param>
        /// <param name="args"></param>
        public void WarnFormat(string format, params object[] args) {
        }

        #endregion

        #region ILoggerWrapper Members

        /// <summary>
        /// Not implemented
        /// </summary>
        public ILogger Logger {
            get {
                throw new NotImplementedException();
            }
        }

        #endregion
    }
}
