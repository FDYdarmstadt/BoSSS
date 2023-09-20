using System;
using ApplicationWithIDT;

namespace SAIDT
{
    public class SAIDTControl : IDTControl
    {
        /// <summary>
        /// defines 
        /// </summary>
        public SAIDTControl()
        {
        }
        public double LeftValue { get; set; }
        public double RightValue { get; set; }
        public double ShockPos { get; set; }
        public Func<double, double> FlowFunc { get; set; }
        public bool UseP0ProjectionAsInitialValue { get; set; }
        public override Type GetSolverType() {
            return typeof(SAIDTMain);
        }
    }
}