using ilPSP;
using ilPSP.Tracing;
using ilPSP.Utils;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;

namespace BoSSS.Solution.XdgTimestepping {

    /// <summary>
    /// Core functionality of adaptive timestepping
    /// Calling <see cref="Compute"/> spans a tree of nested timesteps.
    /// </summary>
    public class TimeLevel {

        public TimeLevel(XdgTimestepping __owner, double __dt, StateAtTime State0) {
            m_owner = __owner;
            this.TimeLevels[0] = State0;
            m_physTime = State0.time;
            m_dt = __dt;
        }


        private TimeLevel(XdgTimestepping __owner, double __physTime, double __dt) {
            m_owner = __owner;
            m_physTime = __physTime;
            m_dt = __dt;
        }



        TimeLevel m_Parrent;
        TimeLevel m_Child;


        XdgTimestepping m_owner;
        double m_physTime;
        double m_dt;

        public int iLevel {
            get {
                if(m_Parrent == null)
                    return 1;
                else
                    return m_Parrent.iLevel + 1;
            }
        }

        public int SubdivExponent {
            get {
                return 2;
            }
        }

        /// <summary>
        /// Number of single steps in this level
        /// </summary>
        public int NoOfTimesteps {
            get {
                int N = iLevel;

                int a = 1;
                for(int n = 1; n < N; n++) {
                    a *= SubdivExponent;
                }

                if(iLevel == 1)
                    Debug.Assert(a == 1);

                return a;
            }
        }

        public double dtSub {
            get {
                return m_dt / NoOfTimesteps;
            }
        }

        void Subdivide() {
            if(m_Child != null)
                throw new NotSupportedException();

            m_Child = new TimeLevel(m_owner, m_physTime, m_dt) {
                m_Parrent = this
            };

            m_Child.TimeLevels[0] = this.TimeLevels[0];
        }

        StateAtTime[] m_TimeLevels;

        StateAtTime[] TimeLevels {
            get {
                if(m_TimeLevels == null)
                    m_TimeLevels = new StateAtTime[NoOfTimesteps + 1];
                return m_TimeLevels;
            }
        }

        double[] m_ScaledDeltas;

        double[] ScaledDeltas {
            get {
                if(m_ScaledDeltas == null) {
                    m_ScaledDeltas = new double[NoOfTimesteps + 1];
                    m_ScaledDeltas.SetAll(-1.0);
                }
                return m_ScaledDeltas;
            }
        }




        StateAtTime LatestSol {
            get {
                return TimeLevels[ComputedSteps - 1];
            }
        }

        /// <summary>
        /// Number of (sub-) timesteps already performed in this level.
        /// Rem.: the computation is finished when <see cref="ComputedSteps"/> equals <see cref="NoOfTimesteps"/>.
        /// </summary>
        int ComputedSteps {
            get {
                int maxIdx = -1;
                for(int i = 0; i < TimeLevels.Length; i++) {
                    if(TimeLevels[i] != null)
                        maxIdx = i;
                }
                Debug.Assert(maxIdx >= 0);
                return maxIdx;
            }
        }



        int iParrent(int iThis) {
            int ip = iThis / this.SubdivExponent;
            return ip;
        }


        /// <summary>
        /// Next range of timesteps to compute
        /// </summary>
        void FindRange(out int i0, out int iE) {
            if(iLevel == 1) {
                Assert.IsTrue(m_Parrent == null, "internal error, first level must have index 1.");
                i0 = 0;
                iE = 1;
                return;
            }


            i0 = ComputedSteps;
            int i0_Parrent = i0 / SubdivExponent;
            //if(i0 % SubdivExponent != 0)
            //    throw new ApplicationException("error in algorithm");

            int iE_Parrent = m_Parrent.NextAvailable(i0_Parrent);

            iE = iE_Parrent * SubdivExponent;
        }

        int NextAvailable(int i0) {
            if(i0 < 0)
                throw new ArgumentOutOfRangeException();
            if(i0 >= NoOfTimesteps)
                throw new ArgumentOutOfRangeException();

            for(int i = i0 + 1; i < TimeLevels.Length; i++)
                if(TimeLevels[i] != null)
                    return i;

            if(m_Parrent == null) {
                Debug.Assert(iLevel == 1);
                Debug.Assert(NoOfTimesteps == 1);
                return 1;
            } else {
                int i0_Parrent = i0 / SubdivExponent;
                int iE_Parrent = m_Parrent.NextAvailable(i0_Parrent);
                return iE_Parrent * SubdivExponent;
            }

            //throw new NotSupportedException("unable to find next timestep.");
        }

        (int lv, StateAtTime b) GetParrentSol(int iThisLevel) {
            if(iThisLevel % SubdivExponent != 0)
                throw new ApplicationException("solution not available on parent level");
            var r = m_Parrent.TimeLevels[iThisLevel / SubdivExponent];
            if(r == null)
                return m_Parrent.GetParrentSol(iThisLevel / SubdivExponent);
            else
                return (m_Parrent.iLevel, r);
        }

        (int iLv, StateAtTime sol) GetBestSol(int iThisLevel) {
            StateAtTime ret = null;
            int iLv = -1;

            // search downward
            if(m_Child != null) {
                (iLv, ret) = m_Child.GetBestSolDownward(iThisLevel * SubdivExponent);
                if(ret != null) {
                    return (iLv, ret);
                }
            }

            // if nothing found yet, try current level next.
            ret = this.TimeLevels[iThisLevel];
            if(ret != null)
                return (this.iLevel, ret);

            // if everything fails, search upward
            (iLv, ret) = GetBestSolUpward(iThisLevel);

            return (iLv, ret);
        }

        (int iLv, StateAtTime sol) GetBestSolUpward(int iThisLevel) {
            var ret = this.TimeLevels[iThisLevel];
            if(ret != null)
                return (this.iLevel, ret);

            if(m_Parrent == null)
                return (-1, null);

            if(iThisLevel % SubdivExponent != 0)
                // the index cannot be translated to top level
                return (-1, null);

            return m_Parrent.GetBestSolUpward(iThisLevel / SubdivExponent);
        }

        (int iLv, StateAtTime sol) GetBestSolDownward(int iThisLevel) {
            StateAtTime ret = null;
            int iLv = -1;
            if(m_Child != null)
                (iLv, ret) = m_Child.GetBestSolDownward(iThisLevel * SubdivExponent);
            if(ret != null)
                return (iLv, ret);

            return (this.iLevel, TimeLevels[iThisLevel]);
        }






        /// <summary>
        /// Recursive solver call
        /// </summary>
        public void Compute() {
            using (new FuncTrace()) {

                if (iLevel == 1) {

                    FindRange(out int i0, out int _);
                    Debug.Assert(i0 == 0);
                    (int StartSolLevel, var StartSol) = GetBestSol(i0);
                    StartSol.Apply(this.m_owner);
                    Assert.IsTrue(DoubleExtensions.ApproxEqual(StartSol.time, m_physTime + dtSub * i0), "Wrong time level of initial value.");

                    Console.Write($"Level {iLevel}, rst_lv {StartSolLevel} dt_sub = {dtSub}: ");
                    Console.Write("0");
                    m_owner.TimesteppingBase.Solve(m_physTime, dtSub);
                    TimeLevels[i0 + 1] = StateAtTime.Obtain(m_owner, m_physTime + dtSub);
                    Assert.IsTrue(DoubleExtensions.ApproxEqual(TimeLevels[i0 + 1].time, m_physTime + dtSub), "wrong time level of solution, Level 1");
                    Console.WriteLine(".");

                    //this.Subdivide();
                    //m_Child.Compute();

                } else {

                    while (ComputedSteps < NoOfTimesteps) {
                        // still work to do on this level                        

                        FindRange(out int i0, out int iE);

                        // initial value:
                        (int StartSolLevel, var StartSol) = GetBestSol(i0);
                        Assert.IsTrue(StartSol.time.ApproxEqual(m_physTime + dtSub * i0), "Wrong time level of initial value.");
                        StartSol.Apply(this.m_owner);

                        // coarser solution on top time level to compare with:
                        (int CoarseLevel, var CoarseSol) = GetParrentSol(iE);
                        Assert.IsTrue(CoarseSol.time.ApproxEqual(m_physTime + dtSub * iE), "Wrong time level of coarse solution");


                        if (CoarseLevel < this.iLevel - 1)
                            return; // finish top level first

                        // perform sub-timesteps
                        Console.Write($"Level {iLevel}, rst_lv {StartSolLevel} dt_sub = {dtSub}: ");
                        for (int i = i0; i < iE; i++) {
                            double physTime = m_physTime + dtSub * i;
                            Console.Write(i);
                            m_owner.TimesteppingBase.Solve(physTime, dtSub);
                            TimeLevels[i + 1] = StateAtTime.Obtain(m_owner, physTime + dtSub);
                            Assert.IsTrue(TimeLevels[i + 1].time.ApproxEqual(physTime + dtSub), "wrong time level of solution");
                            Console.Write(".");
                        }
                        //Console.WriteLine();

                        // improved solution
                        var FineSol = TimeLevels[iE];
                        Debug.Assert(CoarseSol.time.ApproxEqual(FineSol.time));

                        // comparison:
                        double Delta = GenericBlas.L2Dist(CoarseSol.SolutionState, FineSol.SolutionState);
                        //double Delta = m_owner.ComputeL2Dist(CoarseSol.u0, FineSol.SolutionState);
                        double ScaledDelta = Delta * (1 / dtSub) * iE;
                        ScaledDeltas[iE] = ScaledDelta;
                        double RedFactor = maxParrentScaledDelta() / ScaledDelta;
                        Console.WriteLine("   Delta vs Lv" + CoarseLevel + ": " + ScaledDelta + " fact: " + RedFactor);


                        // subdivision:
                        if (iLevel <= 10 && (iLevel <= 2 || RedFactor > 2)) {
                            if (m_Child == null)
                                this.Subdivide();
                            m_Child.Compute();
                        }
                    }
                }
            }
        }


        double maxParrentScaledDelta() {
            double[] deltas = m_Parrent.ScaledDeltas.Where(delta => delta > 0).ToArray();
            if(deltas.Length <= 0)
                return 1e100;
            else
                return deltas.Max();
        }
    }


}
