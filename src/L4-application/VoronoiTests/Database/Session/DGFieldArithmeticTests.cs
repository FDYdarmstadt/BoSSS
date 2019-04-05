using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using BoSSS.Solution.Statistic;
using BoSSS.Foundation.IO;
using NUnit.Framework;

namespace VoronoiTests.Database.Session
{
    class DGFieldArithmeticTests : SessionTest
    {
        delegate void DGFieldComparisonNonEmb_ComputeError(IEnumerable<string> fieldNames,
                ITimestepInfo[] timesteps,
                out double[] gridResolutions,
                out Dictionary<string, int[]> DOFs,
                out Dictionary<string, double[]> errors,
                out Guid[] guids);

        class ErrorInfo
        {
            public double[] gridResolution;
            public Dictionary<string, int[]> DOFs;
            public Dictionary<string, double[]> Error;
            public Guid[] guids;
        }

        static ErrorInfo ComputeErrors(DGFieldComparisonNonEmb_ComputeError errorEvaluator)
        {
            string fieldName = "T";
            ISessionInfo[] sessions = GetSessions();
            ITimestepInfo[] timesteps = sessions.Select(s => s.Timesteps.Last()).ToArray();
            errorEvaluator(new[] { fieldName },
                timesteps,
                out double[] gridResolution,
                out Dictionary<string, int[]> DOFs,
                out Dictionary<string, double[]> Error,
                out Guid[] guids);
            ErrorInfo FieldComparisonError = new ErrorInfo{
                gridResolution = gridResolution,
                DOFs = DOFs,
                Error = Error,
                guids = guids};
            return FieldComparisonError;
        }

        public override void Run()
        {
            L2Error();
        }

        [Test]
        public static void L2Error()
        {
            ErrorInfo computedInfo = ComputeErrors(DGFieldComparisonNonEmb.ComputeErrors_L2);
        }

        [Test]
        public static void H1Error()
        {
            ErrorInfo computedInfo = ComputeErrors(DGFieldComparisonNonEmb.ComputeErrors_H1);
        }

        [Test]
        public static void L2noMeanError()
        {
            ErrorInfo computedInfo = ComputeErrors(DGFieldComparisonNonEmb.ComputeErrors_L2noMean);
        }
    }
}
