using System;

namespace bcl {


    /// <summary>
    /// thrown if something is weird with the users input
    /// </summary>
    internal class UserInputException : ApplicationException {

        /// <summary>
        /// default constructor
        /// </summary>
        /// <param name="msg">
        /// not allowed to be empty
        /// </param>
        public UserInputException(string msg) 
            : base(msg) {
            if (msg == null || msg.Length <= 0)
                throw new ArgumentException("msg cannot be empty;", "msg");
        }

        /// <summary>
        /// <see cref="ApplicationException"/>
        /// </summary>
        /// <param name="msg"><see cref="ApplicationException"/></param>
        /// <param name="inner"><see cref="ApplicationException"/></param>
        public UserInputException(string msg, Exception inner) : base(msg, inner) { }
    }

    /// <summary>
    /// thrown if something is wrong with the BoSSS enviroment (e.g. missing folders,...);
    /// </summary>
    internal class EnviromentException : ApplicationException {

        /// <summary>
        /// 
        /// </summary>
        /// <param name="msg"></param>
        public EnviromentException(string msg)
            : base(msg) {
            if (msg == null || msg.Length <= 0)
                throw new ArgumentException("msg cannot be empty;", "msg");
        }

        /// <summary>
        /// <see cref="ApplicationException"/>
        /// </summary>
        /// <param name="msg"><see cref="ApplicationException"/></param>
        /// <param name="inner"><see cref="ApplicationException"/></param>
        public EnviromentException(string msg, Exception inner)
            : base(msg, inner) {
            if (msg == null || msg.Length <= 0)
                throw new ArgumentException("msg cannot be empty;", "msg");
        }

    }
}
