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
using System.Linq;
using System.Text;
using BoSSS.Foundation;
using BoSSS.Platform;
using BoSSS.Platform.Utils;
using ilPSP;
using ilPSP.Utils;
using MPI.Wrappers;
using System.Diagnostics;
using MathNet.Numerics.Algorithms.LinearAlgebra;
using System.Numerics;

namespace BoSSS.Solution.Timestepping {

    /// <summary>
    /// ROCK4 adaptive timestepping
    /// </summary>
    public partial class ROCK4 : ITimeStepper {

        /// <summary>
        /// ctor.
        /// </summary>
        public ROCK4(SpatialOperator op, CoordinateVector V, CoordinateMapping ParamFields) {
            this.CurrentState = V;
            this.OpEv = op.GetEvaluatorEx(this.Mapping, 
                ParamFields != null ? ParamFields.Fields : new DGField[0], 
                this.Mapping);
            this.nfevals = 0;
            this.nrejected = 0;
        }

        IEvaluatorNonLin OpEv;


        /// <summary>
        /// current simulation time
        /// </summary>
        public double Time {
            get;
            private set;
        }

        /// <summary>
        /// set <see cref="Time"/> to <paramref name="NewTime"/>
        /// </summary>
        public void ResetTime(double NewTime, int timestepNumber) {
            this.Time = NewTime;
        }

        int cnt = 1;

        /// <summary>
        /// performs adaptive timesteps until 
        /// step-size <paramref name="h"/> is reached.
        /// </summary>
        /// <param name="h"></param>
        public double Perform(double h) {

            //throw new NotImplementedException();

            
            double hp = this.Timestep;

            double[] Phi = this.CurrentState.ToArray();

            int mdeg = 1;
            this.cnt++;
            if (mdeg > 152) {
                mdeg=152;
            }
            mdeg = Math.Max(mdeg, 5) - 4;
            int mr, mz;
            mdegr(ref mdeg, out mz, out mr);

            //Console.WriteLine("\nmdeg={0}  \n", mdeg);

            this.rfstep(mdeg, mr, mz, h, ref Phi);


            this.Timestep = h;
            this.PrevTimestep = hp;
            this.Time += h;

            this.CurrentState.CopyFrom(Phi, 0);
            return h;
             
        }

        double rtol = 1e-3;

        /// <summary>
        /// relative tolerance used in error estimate
        /// </summary>
        public double Reltol {
            get {
                return rtol;
            }
            set {
                if (value <= 0)
                    throw new ArgumentException("must be strictly positive.");
                this.rtol = value;
            }
        }


        double atol = 1e-6;

        /// <summary>
        /// absolute tolerance used in error estimate
        /// </summary>
        public double Abstol {
            get {
                return atol;
            }
            set {
                if (value <= 0)
                    throw new ArgumentException("must be strictly positive.");
                this.atol = value;
            }
        }

        /// <summary>
        /// Accumulated number of DG spatial operator evaluation.
        /// </summary>
        public int nfevals {
            get;
            private set;
        }

        

        private void Rock4_Driver<U, W>(double time, double alpha, U Phi, double beta, W Fode)
            where U : IList<double>
            where W : IList<double> {
            //Globals;
            //Phi=reshape(Phi,Np,K);
            //% SPZ=reshape(SPZ,Np,K);
            //[RhsHem] = Rock4_Rhs(Phi,SPZ);
            //RhsHem=RhsHem(:);
            this.CurrentState.Clear();
            this.CurrentState.Acc(1.0, Phi);
            this.OpEv.Evaluate(-1.0*alpha, beta, Fode);
            Fode.CheckForNanOrInfV(true, true, true);
            this.nfevals++;
        }

        /// <summary>
        /// <paramref name="k"/>-step Arnoldi 
        /// </summary>
        /// <param name="Phi"></param>
        /// <param name="k">number of Arnoldi steps requested</param>
        /// <param name="kact">actual Arnoldi steps taken</param>
        /// <param name="lambda"></param>
        /// <returns></returns>
        private double Calculate_Eigennvalue<U>(U Phi, int k, out int kact, double lambda = 1.0e-4)
            where U : IList<double> {
            int n = Phi.Count;
            if (n != this.Mapping.LocalLength)
                throw new ArgumentException();

            
            //bool reorth = false;
            double[][] V = new double[k + 1][];
            double[,] H = new double[k + 1, k];
            double Eigen;

            {
                double[] v0 = new double[n];
                Random R = new Random();
                for (int i = 0; i < n; i++)
                    v0[i] = R.NextDouble();

                double nrm2 = v0.L2Norm();
                if (nrm2 == 0.0)
                    throw new ArithmeticException("arnoldi: input v0 is a zero vector");

                v0.ScaleV(1.0/nrm2);
                V[0] = v0;
            }

            //double[] vA = new double[n];
            //double[] vB = new double[n];

            double tol = n*double.Epsilon;
            //double time = 0;

            double[] Fode_at_phi = new double[n];
            Rock4_Driver(this.Time, 1.0, Phi, 0.0, Fode_at_phi);
            double[] Utmp = new double[n];

            for (int j = 0; j < k; j++) {

                // vh= (Rock4_Driver(time,Phi(:)+lamda*V(:,j))-Rock4_Driver(time,Phi(:)))/lamda;
                double[] vh = new double[n];
                Utmp.SumV(0.0, 1.0, Phi, lambda, V[j]);
                Rock4_Driver(this.Time, 1.0, Utmp, 1.0/lambda, vh);
                vh.AccV(-1.0/lambda, Fode_at_phi);

                double nrmvh =  vh.MPI_L2Norm();

                //  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
                //   by MGS
                for (int i = 0; i < j; i++) {

                    double hij = ParallelBlas.ddot(n, V[i], 1, vh, 1, csMPI.Raw._COMM.WORLD); //   ( V(:,i) )' * vh;
                    //vh = vh - hij*V(:,i);:
                    vh.AccV(-hij, V[i]);
                    H[i, j] = hij;
                }

                //if(reorth == true) {
                //    for(int i= 0; i < j; i++) {
                //        double tmp = ( V(:,i) )' * vh;
                //       vh = vh - tmp*V(:,i);
                //    H(i,j)=H(i,j)+tmp;
                //}

                //%   <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

                double nrmvhOrtho = vh.MPI_L2Norm();
                H[j+1, j]= nrmvhOrtho;
                vh.ScaleV(1.0/nrmvhOrtho);

                V[j+1] = vh;
                if (nrmvhOrtho <= tol*nrmvh) {
                    Console.WriteLine("Arnoldi termination at step: " + j);
                    kact = j;

                    H = H.GetSubMatrix(0, k, 0, k);//  (1:k,1:k);
                    //Eigen = max(abs(eig(H)));
                    //return Eigen;
                    
                    //MathNet.Numerics.Algorithms.LinearAlgebra.Mkl.MklLinearAlgebraProvider

                    Eigen = EigenScheisse(H);
                    return Eigen;
                }
            }



            kact = k;
            H = H.GetSubMatrix(0, k, 0, k);
            Eigen = EigenScheisse(H);
            return Eigen;
        }

        /// <summary>
        /// Maximum of the absolute value of all Eigenvalues of <paramref name="H"/>.
        /// </summary>
        static public double EigenScheisse(double[,] H) {
            var linalg = new ManagedLinearAlgebraProvider();

            double Eigen;
            int N = H.GetLength(0);
            double[] Matrix = H.Resize(false);
            double[] EigenVect = new double[Matrix.Length];
            double[] diagM = new double[Matrix.Length];
            Complex[] EigenValues = new Complex[N];
            linalg.EigenDecomp(false, N, Matrix, EigenVect, EigenValues, diagM);
            Eigen = EigenValues.Select(ev => Complex.Abs(ev)).Max();
            return Eigen;
        }

        /// <summary>
        /// Find the optimal degree.
        /// </summary>
        /// <param name="mdeg"></param>
        /// <param name="mp1">
        /// pointer which select the degree in <see cref="ms"/>[i]\1,2,..
        /// such that <paramref name="mdeg"/> is less or equal than <see cref="ms"/>[i]
        /// </param>
        /// <param name="mp2">
        /// pointer which gives the corresponding position
        /// of a_1 in the data <see cref="recf"/> for the selected degree.
        /// </param>
        private static void mdegr(ref int mdeg, out int mp1, out int mp2) {

            // -------- Find the degree.--------     
            mp1 = 0;
            mp2 = 0;
            int i = 0;

            while (i < 50) {
                if ((ms[i]/mdeg) >= 1) {
                    mdeg=ms[i];
                    mp1 = i;
                    return;
                } else {
                    mp2 = mp2 + ms[i]*2-1;
                }

                i++;
            }
        }

     

        /// <summary>
        /// counts rejected timesteps
        /// </summary>
        public int nrejected {
            get;
            private set;
        }

        /// <summary>
        /// Initial step size guess, usually between 1d-4 and 1d-6.
        /// </summary>
        private double Timestep = 1.0e-5;


        private double PrevTimestep;

        private const double UROUND = 1.0e-16;

        /// <summary>
        /// cyclic counter which triggers Eigenvalue evaluation
        /// </summary>
        int nrho = 0;

        /// <summary>
        /// maximum Eigenvalue of spatial operator
        /// </summary>
        double eigmax = double.NaN;


        const int EigenvaluePeriod = 10;

        /// <summary>
        /// error in previous timestep, required for step size strategy 'with memory'.
        /// </summary>
        double errp = 0.0;


        
        /// <summary>
        /// number of accepted steps
        /// </summary>
        int NoAccepted = 0;

        /// <summary>
        /// number of rejected steps
        /// </summary>
        int NoRejected = 0;

        /// <summary>
        /// 
        /// </summary>
        int NoOfSteps = 0;


        double facmax=5.0;


        /// <summary>
        /// performs one timestep
        /// </summary>
        /// <param name="dtMax">upper limit for the timestep</param>
        /// <returns>the timestep which was actually taken</returns>
        public double PerformAdaptive(double dtMax) {
            if( dtMax <= UROUND)
                throw new ArgumentException("maximum timestep must be strictly positive");
            
            if(NoOfSteps <= 0)
                this.PrevTimestep = this.Timestep;


            double[] PhiOld = this.CurrentState.ToArray(); // backup of current state
            double[] Phi;
            double h = Math.Min(this.Timestep, dtMax);
            double hp = this.PrevTimestep;
            bool reject = false;

            int nrej = 0; // number of subsequent rejects
            do {
                Phi = PhiOld.CloneAs();

                // Approximate Initial number of steps
                // ------------------------------------

                if (nrho == 0) {
                    int kact;
                    eigmax = this.Calculate_Eigennvalue(Phi, 20, out kact);
                }
                
                int mdeg = (int)Math.Round(Math.Sqrt((3.00+h*eigmax)/0.3530)+1.0);
                //  Correct initial degree info.
                if (mdeg > 152) {
                    h = 0.80*(152.0.Pow2()*0.3530-3.0)/eigmax;
                    mdeg=152;
                    h = Math.Min(h, dtMax);
                }
                mdeg = Math.Max(mdeg, 5) - 4;
                int mr, mz;
                mdegr(ref mdeg, out mz, out mr);

                // Perform integration
                // -------------------

                double err = rfstep(mdeg, mr, mz, h, ref Phi);
                NoOfSteps++;


                // error control procedure
                // -----------------------

                // (seems messy...)

                double fac = (1.0/err).Pow(0.250);
                if (errp != 0.0 && !reject) {
                    double facp = errp.Pow(0.25)*fac.Pow2()*(h/hp);
                    fac = Math.Min(fac, facp);
                }

                if (reject) {
                    facmax = 1.0;
                }

                fac = Math.Min(facmax, Math.Max(0.1, 0.8*fac));
                double hnew = h*fac;


                if (err < 1.0) {
                    // accepted step
                    NoAccepted++;
                    nrho++;
                    nrho = nrho%EigenvaluePeriod;
                    facmax = 2.0;

                    if (reject) {
                        hnew = Math.Min(Math.Min(hnew, h), dtMax);
                        reject = false;
                        nrej = 0;
                    }
                    hp = h;
                    h = hnew;


                } else {
                    // rejected step
                    NoRejected++;
                    nrej++;
                    reject = true;
                    h = 0.8*hnew;
                    if (NoOfSteps <= 1)
                        h = 0.1*h;

                    if (this.nrho != 0)
                        this.nrho = 0;
                    else
                        this.nrho = 1;

                    if (nrej >= 10)
                        h = 1.0E-5;

                    Array.Copy(PhiOld, Phi, Phi.Length);
                }

            } while (reject); // after successful timestep, we return 


            this.Timestep = h;
            this.PrevTimestep = hp;
            this.Time += h;

            this.CurrentState.CopyFrom(Phi, 0);


            return hp;
        }


        /// <summary>
        /// Solution at <see cref="Time"/>+<paramref name="h"/> by an explicit (<paramref name="mdeg"/>+4)-stages formula.
        /// </summary>
        /// <param name="mdeg">number of stages minus 4.</param>
        /// <param name="mr">
        /// pointer into <see cref="recf"/>
        /// </param>
        /// <param name="mz">
        /// pointer into the first index of
        /// <see cref="fpa"/>,  <see cref="fpb"/>,  <see cref="fpbe"/>;
        /// </param>
        /// <param name="h">
        /// timestep
        /// </param>
        /// <param name="Phi">
        /// input/output:
        /// on entry: initial value;
        /// on exit: new value
        /// </param>
        private double rfstep(int mdeg, int mr, int mz, double h, ref double[] Phi) {

            // mem alloc
            // ---------------------
            int neqn = Phi.Length;
            double temp1, temp2, temp3, temp4, temp5, ci1, ci2, ci3;
            double[] Phi_jm2 = new double[neqn];
            double[] Phi_jm1 = new double[neqn];
            double[] Phi_jm3 = new double[neqn];
            double[] Phi_jm4 = new double[neqn];
            double[] fnt = new double[neqn];
            
            // First stage
            // ------------------------
            temp1 = h*recf[mr];
            ci1 = this.Time + temp1;
            ci2 = this.Time + temp1;
            ci3 = this.Time;

            Phi_jm2.CopyEntries(Phi);                        // ROCK4 paper, eq. 3.4, Zeile 1
            Rock4_Driver(this.Time, 1.0, Phi_jm2, 0.0, fnt); // ROCK4 paper, eq. 3.4, Zeile 2
            Phi_jm1.SumV(0.0, 1.0, Phi, temp1, fnt);         // ROCK4 paper, eq. 3.4, Zeile 2


            if (mdeg < 2) {
                Phi.CopyEntries(Phi_jm1);
            }

            // Stages for j=2..mdeg.-------- (ROCK4 paper, eq. 3.4)
            // --------------------------------------------------------
            for (int j = 1; j < mdeg; j++) {
                temp1 = h*recf[mr+2*(j-2+1)+1];
                temp3 = -recf[mr+2*(j-2+1)+2];
                temp2 = 1.0 - temp3;

                //Phi_new = Rock4_Driver(ci1,Phi_jm1);
                Rock4_Driver(ci1, 1.0, Phi_jm1, 0.0, Phi);

                ci1 = temp1 + temp2*ci2 + temp3*ci3;
                //Phi_new=[Phi_new,Phi_jm1,Phi_jm2]*[temp1;temp2;temp3];
                Phi.SumV(temp1, temp2, Phi_jm1, temp3, Phi_jm2);


                // -------- Shift the value "y" for the next stage.--------

                if (j < (mdeg - 1)) {
                    var buf = Phi_jm2;
                    Phi_jm2 = Phi_jm1;
                    Phi_jm1 = Phi;
                    Phi = buf;
                }

                ci3=ci2;
                ci2=ci1;
            }

            // The finishing procedure (4-stage method).
            // -----------------------------------------------------------

            // -------- Stage 1.--------
            temp1 = h*fpa[mz, 0];
            Rock4_Driver(ci1, 1.0, Phi, 0.0, Phi_jm1); //  Phi_jm1 = Rock4_Driver(ci1,Phi_new);
            Phi_jm3.SumV(0.0, 1.0, Phi, temp1, Phi_jm1); // Phi_jm3 = Phi_new + temp1*Phi_jm1;

            // -------- Stage 2.--------
            ci2 = ci1 + temp1;
            temp1 = h*fpa[mz, 1];
            temp2 = h*fpa[mz, 2];
            Rock4_Driver(ci2, 1.0, Phi_jm3, 0.0, Phi_jm2); //Phi_jm2 = Rock4_Driver(ci2,Phi_jm3);
            Phi_jm4.SumV(0.0, 1.0, Phi, temp1, Phi_jm1, temp2, Phi_jm2);  //Phi_jm4 = Phi_new+temp1*Phi_jm1+temp2*Phi_jm2;

            // -------- Stage 3.--------
            ci2 = ci1 + temp1 + temp2;
            temp1 = h*fpa[mz, 3];
            temp2 = h*fpa[mz, 4];
            temp3 = h*fpa[mz, 5];
            Rock4_Driver(ci2, 1.0, Phi_jm4, 0.0, Phi_jm3);  //Phi_jm3 = Rock4_Driver(ci2,Phi_jm4);
            fnt.SumV(0.0, 1.0, Phi, temp1, Phi_jm1, temp2, Phi_jm2, temp3, Phi_jm3); //fnt = Phi_new + temp1*Phi_jm1 + temp2*Phi_jm2 + temp3*Phi_jm3;

            // -------- Stage 4.--------
            ci2 = ci1 + temp1 + temp2 + temp3;
            temp1=h*fpb[mz, 0];
            temp2=h*fpb[mz, 1];
            temp3=h*fpb[mz, 2];
            temp4=h*fpb[mz, 3];
            Rock4_Driver(ci2, 1.0, fnt, 0.0, Phi_jm4);  //Phi_jm4 = Rock4_Driver(ci2,fnt);
            Phi.SumV(1.0, temp1, Phi_jm1, temp2, Phi_jm2, temp3, Phi_jm3, temp4, Phi_jm4);  //Phi_new = [Phi_new,Phi_jm1,Phi_jm2,Phi_jm3,Phi_jm4]*[1;temp1;temp2;temp3;temp4];

            // Error evaluation (embedded method of order 3)
            // ---------------------------------------------------------
            temp1 = h*fpbe[mz, 0] - temp1;
            temp2 = h*fpbe[mz, 1] - temp2;
            temp3 = h*fpbe[mz, 2] - temp3;
            temp4 = h*fpbe[mz, 3] - temp4;
            temp5 = h*fpbe[mz, 4];
            Rock4_Driver(Time + h, 1.0, Phi, 0.0, fnt);  //fnt = Rock4_Driver(time+dt,Phi_new);
            //Phi_jm2.SumV(0.0, temp1, Phi_jm1, temp2, Phi_jm2, temp3, Phi_jm3, temp4, Phi_jm4, temp5, fnt);  //Phi_jm2 = [Phi_jm1,Phi_jm2,Phi_jm3,Phi_jm4,fnt]*[temp1;temp2;temp3;temp4;temp5];

            double err;
            {
                err = 0;
                double ato = this.atol;
                double rto = this.rtol;
                for(int j = 1; j < neqn; j++) {
                    ci1 = Math.Abs(Phi[j])*rto;
                    err = err + ((temp1*Phi_jm1[j] + temp2*Phi_jm2[j] + temp3*Phi_jm3[j] + temp4*Phi_jm4[j] + temp5*fnt[j])/(ato+ci1)).Pow2();
        
                }
                err = err.MPISum();
                err = Math.Sqrt(err / ((double)this.Mapping.GlobalCount));
            }

            return err;
        }

        /// <summary>
        /// mapping for the DG fields on which this time integrator acts on.
        /// </summary>
        public CoordinateMapping Mapping {
            get { return CurrentState.Mapping; }
        }

        /// <summary>
        /// DG coordinate vector of the DG fields on which this time integrator acts on.
        /// </summary>
        public CoordinateVector CurrentState {
            get;
            private set;
        }

        public TimeInformation TimeInfo {
            get;
            protected set;
        }

        public void UpdateTimeInfo(TimeInformation timeInfo) {
            this.TimeInfo = timeInfo;
        }
    }
}
