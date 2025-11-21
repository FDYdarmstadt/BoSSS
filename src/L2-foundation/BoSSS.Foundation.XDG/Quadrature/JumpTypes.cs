namespace BoSSS.Foundation.XDG.Quadrature {
   
    /// <summary>
    /// Different type of jumps assumed.
    /// </summary>
    public enum JumpTypes {

        /// <summary>
        /// Function is taken as is; jump is already contained in the function
        /// </summary>
        Implicit,

        /// <summary>
        /// Function is assumed zero for negative level set
        /// </summary>
        Heaviside,

        /// <summary>
        /// Function is assumed zero for positive level set
        /// </summary>
        OneMinusHeaviside,

        /// <summary>
        /// Function is multiplied by -1 for negative level and 1 for positive
        /// level set.
        /// </summary>
        Sign
    }
}
