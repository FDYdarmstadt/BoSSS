
using ApplicationWithIDT;
using System;


namespace BUIDT
{
    public class BUIDTControl : IDTControl
    {
        public BUIDTControl()
        {
        }
        public override Type GetSolverType() {
            return typeof(BUIDTMain);
        }
        /// <summary>
        /// controls smoothness of smoothed upwind flux
        /// </summary>
        public double s_alpha { get; set; } = 10;
        /// <summary>
        /// on-off switch for smoothed upwind flux
        /// </summary>
        public bool is_nf_smth { get; set; } = true;
        /// <summary>
        /// switch for initial guess: true means, that Initial Guess is obtained from p0 projection of exact solution
        /// </summary>
        public bool UseP0ProjectionAsInitialGuess { get; set; }
        /// <summary>
        /// SLSPointPath for initial Guess of Spline Level Set
        /// </summary>
        public string SLSPointPath { get; set; }
        
    }
    

}