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
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Numerics;
using MathNet.Numerics.IntegralTransforms;
using MathNet.Numerics.IntegralTransforms.Algorithms;
using MathNet.Numerics.Interpolation.Algorithms;
using ilPSP;
using ilPSP.Utils;
using BoSSS.Platform;
using BoSSS.Foundation;
using BoSSS.Foundation.Grid;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Solution.Statistic;

namespace BoSSS.Solution.LevelSetTools.FourierLevelSet {


    /// <summary>
    /// options for the movement of the center for PolarFourierLevSet
    /// </summary>
    public enum CenterMovement {

        /// <summary>
        /// no movement of the center point
        /// </summary>
        None,

        /// <summary>
        /// computation of the geometric center from the interface
        /// </summary>
        Reconstructed,

        /// <summary>
        /// evaluation of the velocity at the cneter point
        /// </summary>
        VelocityAtCenter,

    }

    /// <summary>
    /// option for the computation of the level set
    /// </summary>
    public enum PolarLevSetForm {

        /// <summary>
        /// quadratic LevelSet = x^2 + y^2 - r^2  
        /// </summary>
        quadratic,

        /// <summary>
        /// signed distance leveSet = sqrt(x^2 + y^2) - r  
        /// </summary>
        signedDist
    }


    /// <summary>
    /// An explicit Interface description based on control points, which form a Fourier Series in the phi-r space.
    /// This description allows the movement of the reference center point
    /// </summary>
    public class PolarFourierLevSet : FourierLevSetBase {

        /// <summary>
        /// material interface points in polar coordinates
        /// </summary>
        public MultidimensionalArray interfaceP_polar;

        /// <summary>
        /// center point of the polar coordinate system
        /// </summary>
        MultidimensionalArray center;

        /// <summary>
        /// See <see cref="CenterMovement"/>
        /// </summary>
        CenterMovement CenterMove;

        /// <summary>
        /// See <see cref="PolarLevSetForm"/>
        /// </summary>
        PolarLevSetForm LevSetForm;

        /// <summary>
        /// switch if the curvature computation extended from interface
        /// </summary>
        bool curvComp_extended;

        /// <summary>
        /// 
        /// </summary>
        public bool massConsv_correction = false;

        /// <summary>
        /// 
        /// </summary>
        public double consvPolyArea;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="ctrl"></param>
        public PolarFourierLevSet(FourierLevSetControl ctrl) : base(ctrl) {
            this.DomainSize = 2 * Math.PI;

            this.center = new MultidimensionalArray(2);
            this.center.Allocate(1, 2);
            this.center[0, 0] = ctrl.center[0];
            this.center[0, 1] = ctrl.center[1];

            this.CenterMove = ctrl.centerMove;

            this.LevSetForm = ctrl.PLevSetForm;

            this.curvComp_extended = ctrl.curvComp_extended;

            // set material polar coordinates
            interfaceP_polar = new MultidimensionalArray(2);
            interfaceP_polar.Allocate(numFp, 2);
            for (int sp = 0; sp < numFp; sp++) {
                interfaceP_polar[sp, 0] = FourierP[sp];
                interfaceP_polar[sp, 1] = current_samplP[sp];
            }

            setProjectingFourierModes();

            // set the material sample points (cartesian)
            current_interfaceP = new MultidimensionalArray(2);
            current_interfaceP.Allocate(numFp, 2);
            for (int sp = 0; sp < numFp; sp++) {
                current_interfaceP[sp, 0] = center[0, 0] + (interfaceP_polar[sp, 1] * Math.Cos(interfaceP_polar[sp, 0]));
                current_interfaceP[sp, 1] = center[0, 1] + (interfaceP_polar[sp, 1] * Math.Sin(interfaceP_polar[sp, 0]));
            }
            setInterfaceLength(current_interfaceP);
            InterfaceResolution = interfaceLength / numFp;

        }

        /// <summary>
        /// set the projecting Fourier modes for the Level set and curvature evaluation
        /// </summary>
        protected override void setProjectingFourierModes() {

            double r_min = interfaceP_polar.ExtractSubArrayShallow(new int[] { -1, 1 }).To1DArray().Min();
            cf_end = (int)Math.Ceiling(2.0 * Math.PI * r_min / (2.0 * h_min));

        }


        /// <summary>
        /// computes the length of the samnple points polygon 
        /// </summary>
        /// <returns></returns>
        protected override void setInterfaceLength(MultidimensionalArray interP) {

            interfaceLength = 0.0;
            int numP = interP.Lengths[0];
            for (int sp = 0; sp < numP - 1; sp++) {
                interfaceLength += Math.Sqrt((interP[sp + 1, 0] - interP[sp, 0]).Pow2()
                    + (interP[sp + 1, 1] - interP[sp, 1]).Pow2());
            }
            interfaceLength += Math.Sqrt((interP[0, 0] - interP[numP - 1, 0]).Pow2()
                                            + (interP[0, 1] - interP[numP - 1, 1]).Pow2());
        }


        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public double getFourierArea() {
            return Math.PI * (DFT_coeff[0].Real / (double)numFp).Pow2();
        }

        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public double getPolygonalArea(MultidimensionalArray polygon = null) {

            if (polygon == null) {
                // polygonal description of the interface
                int numP = current_interfaceP.Lengths[0];
                polygon = MultidimensionalArray.Create(numP + 1, 2);
                for (int sp = 0; sp < numP; sp++) {
                    polygon[sp, 0] = current_interfaceP[sp, 0];
                    polygon[sp, 1] = current_interfaceP[sp, 1];
                }
                polygon[numP, 0] = polygon[0, 0];
                polygon[numP, 1] = polygon[0, 1];
            }

            double area = 0.0;
            int numPp = polygon.Lengths[0] - 1;
            for (int p = 0; p < numPp; p++) {
                area += polygon[p, 0] * polygon[p + 1, 1] - polygon[p + 1, 0] * polygon[p, 1];
            }
            area *= 1.0 / 2.0;

            return area;
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="polygon"></param>
        /// <returns></returns>
        public double getPolygonalLength(MultidimensionalArray polygon = null) {

            if (polygon == null) {
                // polygonal description of the interface
                int numP = current_interfaceP.Lengths[0];
                polygon = MultidimensionalArray.Create(numP + 1, 2);
                for (int sp = 0; sp < numP; sp++) {
                    polygon[sp, 0] = current_interfaceP[sp, 0];
                    polygon[sp, 1] = current_interfaceP[sp, 1];
                }
                polygon[numP, 0] = polygon[0, 0];
                polygon[numP, 1] = polygon[0, 1];
            }

            double length = 0.0;
            int numPp = polygon.Lengths[0] - 1;
            for (int p = 0; p < numPp; p++) {
                length += Math.Sqrt((polygon[p+1, 0] - polygon[p, 0]).Pow2() + (polygon[p+1, 1] - polygon[p, 1]).Pow2());
            }
            return length;
        }


        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public override double[] GetFLSproperty() {

            switch (FourierEvolve) {
                case Fourier_Evolution.MaterialPoints:
                    int numIp = current_interfaceP.Lengths[0];
                    double[] mpoints = new double[2 * numIp ];
                    for (int sp = 0; sp < numIp; sp++) {
                        mpoints[sp * 2] = current_interfaceP[sp, 0];
                        mpoints[(sp * 2) + 1] = current_interfaceP[sp, 1];
                    }
                    return mpoints;
                case Fourier_Evolution.FourierPoints:
                    double[] center_samplP = new double[numFp + 2];
                    center_samplP[0] = center[0, 0];
                    center_samplP[1] = center[0, 1];
                    for (int i = 0; i < numFp; i++) {
                        center_samplP[i + 2] = current_samplP[i];
                    }
                    return center_samplP;
                case Fourier_Evolution.FourierModes: {
                        double[] center_DFT = new double[2 * numFp + 2];
                        center_DFT[0] = center[0, 0];
                        center_DFT[1] = center[0, 1];
                        for (int sp = 0; sp < numFp; sp++) {
                            center_DFT[sp * 2 + 2] = DFT_coeff[sp].Real;
                            center_DFT[sp * 2 + 3] = DFT_coeff[sp].Imaginary;
                        }
                        return center_DFT;
                    }
                default:
                    throw new ArgumentException();
            }
        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="velocity"></param>
        /// <returns></returns>
        public override double[] ComputeChangerate(double dt, ConventionalDGField[] velocity, double[] current_FLSprop) {

            setMaterialInterfacePoints(current_FLSprop);

            GridData grdat = (GridData)(velocity[0].GridDat);
            FieldEvaluation fEval = new FieldEvaluation(grdat);

            // Movement of the material interface points
            MultidimensionalArray VelocityAtSamplePointsXY = MultidimensionalArray.Create(current_interfaceP.Lengths);

            int outP = fEval.Evaluate(1, velocity, current_interfaceP, 0, VelocityAtSamplePointsXY);

            if (outP != 0)
                throw new Exception("points outside the grid for fieldevaluation");

            int numIp = current_interfaceP.Lengths[0];


            // change rate for the material points is the velocity at the points
            if (FourierEvolve == Fourier_Evolution.MaterialPoints) {
                double[] velAtP = new double[2 * numIp];
                for (int sp = 0; sp < numIp; sp++) {
                    velAtP[sp * 2] = VelocityAtSamplePointsXY[sp, 0];
                    velAtP[(sp * 2) + 1] = VelocityAtSamplePointsXY[sp, 1];
                }
                return velAtP;
            }

            // compute an infinitesimal change of sample points at the Fourier points/ change of Fourier modes
            MultidimensionalArray interfaceP_evo = current_interfaceP.CloneAs();
            double dt_infin = dt * 1e-3;
            interfaceP_evo.Acc(dt_infin, VelocityAtSamplePointsXY);

            // Movement of the center point
            MultidimensionalArray center_evo = center.CloneAs();
            MultidimensionalArray VelocityAtCenter = center.CloneAs();
            switch (CenterMove) {
                case CenterMovement.None: {
                        VelocityAtCenter.Clear();
                        break;
                    }
                case CenterMovement.Reconstructed: {
                        center_evo = GetGeometricCenter(interfaceP_evo);
                        //Console.WriteLine("center_evo = ({0}, {1}) / center = ({2}, {3})", center_evo[0, 0], center_evo[0, 1], center[0, 0], center[0, 1]);
                        MultidimensionalArray center_change = center_evo.CloneAs();
                        center_change.Acc(-1.0, center);
                        center_change.Scale(1.0 / dt_infin);
                        VelocityAtCenter[0, 0] = center_change[0, 0];
                        VelocityAtCenter[0, 1] = center_change[0, 1];
                        break;
                    }
                case CenterMovement.VelocityAtCenter: {
                        outP = fEval.Evaluate(1, velocity, center, 0, VelocityAtCenter);
                        if (outP != 0)
                            throw new Exception("center point outside the grid for fieldevaluation");
                        center_evo.Acc(dt_infin, VelocityAtCenter);
                        break;
                    }
            }
            //Console.WriteLine("Velocity at Center point = ({0}, {1})", VelocityAtCenter[0, 0], VelocityAtCenter[0, 1]);

            // transform to polar coordiantes
            MultidimensionalArray interP_evo_polar = interfaceP_polar.CloneAs();
            for (int sp = 0; sp < numIp; sp++) {
                double x_c = interfaceP_evo[sp, 0] - center_evo[0, 0];
                double y_c = interfaceP_evo[sp, 1] - center_evo[0, 1];

                double theta = Math.Atan2(y_c, x_c);
                if (theta < 0) {
                    theta = Math.PI * 2 + theta;
                };
                interP_evo_polar[sp, 0] = theta;
                interP_evo_polar[sp, 1] = Math.Sqrt(x_c.Pow2() + y_c.Pow2());
            }

            RearrangeOntoPeriodicDomain(interP_evo_polar);
            double[] samplP_change = current_samplP.CloneAs();
            InterpolateOntoFourierPoints(interP_evo_polar, samplP_change);

            if (FourierEvolve == Fourier_Evolution.FourierPoints) {
                samplP_change.AccV(-1.0, current_samplP);
                samplP_change.ScaleV(1.0 / dt_infin);
                double[] FPchange = new double[numFp + 2];
                FPchange[0] = VelocityAtCenter[0, 0];
                FPchange[1] = VelocityAtCenter[0, 1];
                for (int p = 0; p < numFp; p++) {
                    FPchange[p + 2] = samplP_change[p];
                }
                return FPchange;
            } else
            if (FourierEvolve == Fourier_Evolution.FourierModes) {
                Complex[] samplP_complex = new Complex[numFp];
                for (int sp = 0; sp < numFp; sp++) {
                    samplP_complex[sp] = (Complex)samplP_change[sp];
                }
                Complex[] DFTchange = DFT.NaiveForward(samplP_complex, FourierOptions.Matlab);
                double[] DFTchange_double = new double[2 * numFp + 2];
                DFTchange_double[0] = VelocityAtCenter[0, 0];
                DFTchange_double[1] = VelocityAtCenter[0, 1];
                for (int sp = 0; sp < numFp; sp++) {
                    DFTchange_double[2 + (sp * 2)] = (DFTchange[sp].Real - DFT_coeff[sp].Real) / dt_infin;
                    DFTchange_double[2 + (sp * 2) + 1] = (DFTchange[sp].Imaginary - DFT_coeff[sp].Imaginary) / dt_infin;
                }
                return DFTchange_double;
            } else
                throw new ArgumentException();

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="FLSproperty"></param>
        internal void setMaterialInterfacePoints(double[] FLSproperty) {

            //center = MultidimensionalArray.CreateWrapper(FLSproperty.GetSubVector(0, 2), new int[] { 1, 2 });
            switch (FourierEvolve) {
                case Fourier_Evolution.MaterialPoints:
                    int numSp = FLSproperty.Length / 2;
                    current_interfaceP = MultidimensionalArray.CreateWrapper(FLSproperty, new int[] { numSp, 2 });
                    break;
                case Fourier_Evolution.FourierPoints:
                    center = MultidimensionalArray.CreateWrapper(FLSproperty.GetSubVector(0, 2), new int[] { 1, 2 });
                    current_samplP = FLSproperty.GetSubVector(2, numFp);
                    for (int p = 0; p < numFp; p++) {
                        current_interfaceP[p, 0] = current_samplP[p] * Math.Cos(FourierP[p]) + center[0, 0];
                        current_interfaceP[p, 1] = current_samplP[p] * Math.Sin(FourierP[p]) + center[0, 1];
                    }
                    break;
                case Fourier_Evolution.FourierModes:
                    center = MultidimensionalArray.CreateWrapper(FLSproperty.GetSubVector(0, 2), new int[] { 1, 2 });
                    for (int sp = 0; sp < numFp; sp++) {
                        DFT_coeff[sp] = new Complex(FLSproperty[2 + (sp * 2)], FLSproperty[2 + (sp * 2) + 1]);
                    }
                    SmoothSamplePoints();
                    for (int sp = 0; sp < numFp; sp++) {
                        current_interfaceP[sp, 0] = center[0, 0] + current_samplP[sp] * Math.Cos(FourierP[sp]);
                        current_interfaceP[sp, 1] = center[0, 1] + current_samplP[sp] * Math.Sin(FourierP[sp]);
                    }
                    break;
            }

        }


        /// <summary>
        /// 
        /// </summary>
        /// <param name="current_FLSproperty"></param>
        public override void EvolveFourierLS(ref double[] current_FLSproperty) {

            switch (FourierEvolve) {
                case Fourier_Evolution.MaterialPoints:
                    int numSp = current_FLSproperty.Length / 2;
                    current_interfaceP = MultidimensionalArray.CreateWrapper(current_FLSproperty, new int[] { numSp, 2 });
                    //MultidimensionalArray prescribed_center = MultidimensionalArray.CreateWrapper(current_FLSproperty.GetSubVector(0, 2), new int[] { 1, 2 });

                    //Console.WriteLine("center_VelAtC = ({0}, {1})", prescribed_center[0, 0], prescribed_center[0, 1]);


                    SetRefCenter();

                    //bool reseed = Reseeding();
                    //if (reseed == true) {
                    //    Console.WriteLine("reseeding done");
                    //}

                    //RelocateOnInterface();

                    //MultidimensionalArray center_geom = GetGeometricCenter(current_interfaceP);
                    //Console.WriteLine("center_geometric = ({0}, {1})", center_geom[0, 0], center_geom[0, 1]);



                    int numP = current_interfaceP.Lengths[0];
                    current_FLSproperty = new double[numP * 2];
                    for (int sp = 0; sp < numP; sp++) {
                        current_interfaceP[sp, 0] = current_samplP[sp] * Math.Cos(FourierP[sp]) + center[0, 0];
                        current_interfaceP[sp, 1] = current_samplP[sp] * Math.Sin(FourierP[sp]) + center[0, 1];
                        current_FLSproperty[sp * 2] = current_interfaceP[sp, 0];
                        current_FLSproperty[(sp * 2) + 1] = current_interfaceP[sp, 1];
                    }

                    break;
                case Fourier_Evolution.FourierPoints:
                    current_samplP = current_FLSproperty.GetSubVector(2, numFp);

                    setDFTcoeff();
                    SmoothSamplePoints();

                    for (int sp = 0; sp < numFp; sp++) {
                        current_FLSproperty[sp + 2] = current_samplP[sp];
                    }

                    // set material points on Fourier points
                    center = MultidimensionalArray.CreateWrapper(current_FLSproperty.GetSubVector(0, 2), new int[] { 1, 2 });
                    //Console.WriteLine("center_prescribed = ({0}, {1})", center[0, 0], center[0, 1]);
                    for (int p = 0; p < numFp; p++) {
                        current_interfaceP[p, 0] = current_samplP[p] * Math.Cos(FourierP[p]) + center[0, 0];
                        current_interfaceP[p, 1] = current_samplP[p] * Math.Sin(FourierP[p]) + center[0, 1];
                    }
                    //MultidimensionalArray center_interP = GetCenter(current_interfaceP);
                    //double center_diff = Math.Sqrt((center_new[0, 0] - center_interP[0, 0]).Pow2() + (center_new[0, 1] - center_interP[0, 1]).Pow2());
                    //if (center_diff > 1e-9)
                    //    Console.WriteLine("Warning: significant change of the center point after smoothing: center_diff 0 {0}", center_diff);
                    //MultidimensionalArray center_geomFP = GetGeometricCenter(current_interfaceP);
                    //Console.WriteLine("center_geometric = ({0}, {1})", center_geomFP[0, 0], center_geomFP[0, 1]);


                    break;
                case Fourier_Evolution.FourierModes:
                    for (int sp = 0; sp < numFp; sp++) {
                        DFT_coeff[sp] = new Complex(current_FLSproperty[2 + (sp * 2)], current_FLSproperty[2 + (sp * 2) + 1]);
                    }
                    SmoothSamplePoints();

                    // relocate Material points in kartesian coordinates on the interface
                    center = MultidimensionalArray.CreateWrapper(current_FLSproperty.GetSubVector(0, 2), new int[] { 1, 2 });
                    for (int sp = 0; sp < current_interfaceP.Lengths[0]; sp++) {
                        current_interfaceP[sp, 0] = center[0, 0] + current_samplP[sp] * Math.Cos(FourierP[sp]);
                        current_interfaceP[sp, 1] = center[0, 1] + current_samplP[sp] * Math.Sin(FourierP[sp]);
                    }


                    break;
                default:
                    throw new ArgumentException();
            }

            //int kmax = 1;
            //for (int c = 2; c < cf_end; c++) {
            //    if (DFT_coeff[c].Imaginary > DFT_coeff[kmax].Imaginary)
            //        kmax = c;
            //}
            //Console.WriteLine("max Fourier mode imaginary part: k = {0}: {1}", kmax, DFT_coeff[kmax].Imaginary);

        }

        /// <summary>
        /// 
        /// </summary>
        /// <param name="interp_cartesian"></param>
        /// <returns></returns>
        internal MultidimensionalArray GetGeometricCenter(MultidimensionalArray interP_cartesian) {

            MultidimensionalArray cntr = center.CloneAs();
            int numSp = interP_cartesian.Lengths[0];

            // polygonal description of the interface
            MultidimensionalArray interPolygon = MultidimensionalArray.Create(numSp + 1, 2);
            for (int sp = 0; sp < numSp; sp++) {
                interPolygon[sp, 0] = interP_cartesian[sp, 0];
                interPolygon[sp, 1] = interP_cartesian[sp, 1];
            }
            interPolygon[numSp, 0] = interPolygon[0, 0];
            interPolygon[numSp, 1] = interPolygon[0, 1];

            //MultidimensionalArray center_rec = MultidimensionalArray.Create(1, 2);
            cntr.Clear();
            for (int p = 0; p < numSp; p++) {
                double a = interPolygon[p, 0] * interPolygon[p + 1, 1] - interPolygon[p + 1, 0] * interPolygon[p, 1];
                cntr[0, 0] += (interPolygon[p, 0] + interPolygon[p + 1, 0]) * a;
                cntr[0, 1] += (interPolygon[p, 1] + interPolygon[p + 1, 1]) * a;
            }
            double area = getPolygonalArea(interPolygon);
            cntr[0, 0] /= 6 * area;
            cntr[0, 1] /= 6 * area;

            return cntr;

        }


        /// <summary>
        /// sets a new refernce center point and computes the polar coordinates according to the new center point
        /// </summary>
        internal void SetRefCenter() {

            // set new reference center point
            switch (CenterMove) {
                case CenterMovement.None: {
                        break;
                    }
                case CenterMovement.Reconstructed: {
                        center = GetGeometricCenter(current_interfaceP);
                        break;
                    }
                case CenterMovement.VelocityAtCenter:
                    throw new NotSupportedException("Velocity At center is not supported for material points");
                default:
                    throw new NotImplementedException("unknown center movement");
            }


            // transform to polar coordiantes
            //interfaceP_polar = MultidimensionalArray.Create(current_interfaceP.Lengths);
            for (int sp = 0; sp < current_interfaceP.Lengths[0]; sp++) {
                double x_c = current_interfaceP[sp, 0] - center[0, 0];
                double y_c = current_interfaceP[sp, 1] - center[0, 1];

                double theta = Math.Atan2(y_c, x_c);
                if (theta < 0) {
                    theta = Math.PI * 2 + theta;
                };
                interfaceP_polar[sp, 0] = theta;
                interfaceP_polar[sp, 1] = Math.Sqrt(x_c.Pow2() + y_c.Pow2());
            }

            // mass conservation correction step
            if (massConsv_correction) {
                double area = getPolygonalArea();
                double length = getPolygonalLength();

                double cmc = (consvPolyArea - area) / length;

                for (int sp = 0; sp < current_interfaceP.Lengths[0]; sp++) {
                    interfaceP_polar[sp, 1] += cmc; 
                }
            }

            setProjectingFourierModes();

            RearrangeOntoPeriodicDomain(interfaceP_polar);

            InterpolateOntoFourierPoints(interfaceP_polar, current_samplP);

            setDFTcoeff();
            SmoothSamplePoints();

        }


        /// <summary>
        /// 
        /// </summary>
        internal void RelocateOnInterface() {

            int numIp = current_interfaceP.Lengths[0];

            // transform to polar coordiantes
            interfaceP_polar = MultidimensionalArray.Create(new int[] { numIp, 2});
            for (int sp = 0; sp < numIp; sp++) {
                double x_c = current_interfaceP[sp, 0] - center[0, 0];
                double y_c = current_interfaceP[sp, 1] - center[0, 1];

                double theta = Math.Atan2(y_c, x_c);
                if (theta < 0) {
                    theta = Math.PI * 2 + theta;
                };
                interfaceP_polar[sp, 0] = theta;
                interfaceP_polar[sp, 1] = Math.Sqrt(x_c.Pow2() + y_c.Pow2());
            }


            Complex cmplx_res;
            for (int sp = 0; sp < numIp; sp++) {
                interfaceP_polar[sp, 1] = DFT_coeff[0].Real;
                if (cf_end == numFp / 2)
                    interfaceP_polar[sp, 1] += Complex.Multiply(DFT_coeff[numFp / 2],
                                                           Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, interfaceP_polar[sp, 0] * (2 * Math.PI / DomainSize) * (numFp / 2)))
                                                           ).Real;

                for (int cf = 1; cf < cf_end; cf++) {
                    cmplx_res = DFT_coeff[cf] * Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, interfaceP_polar[sp, 0] * (2 * Math.PI / DomainSize) * cf));
                    interfaceP_polar[sp, 1] += 2.0 * cmplx_res.Real;
                }
                interfaceP_polar[sp, 1] /= (double)numFp;
            }

            // relocate Material points in kartesian coordinates on the interface
            for (int sp = 0; sp < current_interfaceP.Lengths[0]; sp++) {
                current_interfaceP[sp, 0] = center[0,0] + interfaceP_polar[sp, 1] * Math.Cos(interfaceP_polar[sp, 0]);
                current_interfaceP[sp, 1] = center[0,1] + interfaceP_polar[sp, 1] * Math.Sin(interfaceP_polar[sp, 0]);
            }

        }


        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public override ScalarFunction PhiEvaluation(int mode) {

            return (ScalarFunction)delegate (MultidimensionalArray nodes, MultidimensionalArray results) {
                //signal is real valued - i.e. for numSp sample points, numSp/2 coefficients are complex conjugated, thus they don't have to be evaluated individually
                double theta;
                double radius;
                //Complex RadiusComplex;
                Complex cmplx_res;
                for (int nd = 0; nd < nodes.GetLength(0); nd++) { // loop over nodes ...
                    double x = nodes[nd, 0] - center[0, 0];
                    double y = nodes[nd, 1] - center[0, 1];


                    theta = Math.Atan2(y, x);
                    if (theta < 0) {
                        theta = Math.PI * 2 + theta;
                    };


                    //r(theta) for  k=0 and k= numSp/2
                    radius = DFT_coeff[0].Real;
                    //RadiusComplex = DFT_coeff[0];
                    if (cf_end == numFp / 2)
                        radius = Complex.Multiply(DFT_coeff[numFp / 2],
                                                               Complex.Exp(Complex.ImaginaryOne * theta * (numFp / 2))
                                                               ).Real;
                    //radius += cmplx_res.Real;
                    //RadiusComplex += cmplx_res;
                    for (int cf = 1; cf < cf_end; cf++) {
                        cmplx_res = Complex.Multiply(DFT_coeff[cf], Complex.Exp(
                                                                        Complex.ImaginaryOne * theta * cf
                                                                        )
                                                     );
                        radius += 2.0 * cmplx_res.Real;  //-f(x) for  k=1... numSp/2-1 , values for k=numSp/2+1...k= numSp-1 are complex conjugated, thus the same real part
                        //RadiusComplex += 2 * cmplx_res;
                    }
                    radius /= (double)numFp;
                    //RadiusComplex /= (double)numFp;

                    //if (radius - RadiusComplex.Real >= 1e-14) throw new ArithmeticException();

                    //radius = 0.0;
                    //for (int cf = 0; cf < numFp / 2; cf++) {
                    //    if (DFT_coeff[cf].Magnitude > significantMode) {
                    //        cmplx_res = Complex.Multiply(DFT_coeff[cf], Complex.Exp(Complex.ImaginaryOne * theta * cf));
                    //        if (cf == 0 || cf == numFp / 2) {
                    //            radius += cmplx_res.Real;
                    //        } else {
                    //            radius += 2.0 * cmplx_res.Real;
                    //        }
                    //    }
                    //}
                    //radius /= (double)numFp;


                    switch (LevSetForm) {
                        case PolarLevSetForm.quadratic: {
                                results[nd] = x.Pow2() + y.Pow2() - radius.Pow2();
                                break;
                            }
                        case PolarLevSetForm.signedDist: {
                                results[nd] = Math.Sqrt(x.Pow2() + y.Pow2()) - radius;
                                break;
                            }
                        default:
                            throw new NotImplementedException("unknown polar levelset description");
                    }

                }
            };
        }


        /// <summary>
        /// 
        /// </summary>
        /// <returns></returns>
        public override ScalarFunction CurvEvaluation(int mode) {

            //Complex[] invDFT_coeff_curv = invDFT_coeff.CloneAs();
            //for (int sp = 0; sp < numFp; sp++) {
            //    Complex cF_x;
            //    double Fr_t = 0.0;
            //    Complex cF_xx;
            //    double Fr_tt = 0.0;

            //    // F_x
            //    if (cf_end == numFp / 2) {
            //        cF_x = Complex.ImaginaryOne * DFT_coeff[numFp / 2] * Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, FourierP[sp] * (numFp / 2)));
            //        Fr_t = (double)(numFp / 2) * cF_x.Real / (double)numFp;
            //    }

            //    for (int cf = 1; cf < cf_end; cf++) {
            //        cF_x = Complex.ImaginaryOne * DFT_coeff[cf] * Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, FourierP[sp] * cf));
            //        Fr_t += cf * 2.0 * cF_x.Real / (double)numFp;
            //    }


            //    //for (int cf = 1; cf < numFp/2; cf++) {
            //    //    if (DFT_coeff[cf].Magnitude > significantMode) {
            //    //        cF_x = Complex.ImaginaryOne * DFT_coeff[cf] * Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, FourierP[sp] * cf));
            //    //        if (cf == numFp / 2) {
            //    //            Fr_t += cf * cF_x.Real / (double)numFp;
            //    //        } else {
            //    //            Fr_t += cf * 2.0 * cF_x.Real / (double)numFp;
            //    //        }
            //    //    }
            //    //}

            //    // F_xx
            //    if (cf_end == numFp / 2) {
            //        cF_xx = -DFT_coeff[numFp / 2] * Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, FourierP[sp] * (numFp / 2)));
            //        Fr_tt = Math.Pow((double)(numFp / 2), 2) * cF_xx.Real / (double)numFp;
            //    }

            //    for (int cf = 1; cf < cf_end; cf++) {
            //        cF_xx = -DFT_coeff[cf] * Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, FourierP[sp] * cf));
            //        Fr_tt += Math.Pow(cf, 2) * 2.0 * cF_xx.Real / (double)numFp;
            //    }

            //    //for (int cf = 1; cf < numFp/2; cf++) {
            //    //    if (DFT_coeff[cf].Magnitude > significantMode) {
            //    //        cF_xx = -DFT_coeff[cf] * Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, FourierP[sp] * cf));
            //    //        if (cf == numFp / 2) {
            //    //            Fr_tt += Math.Pow(cf, 2) * cF_xx.Real / (double)numFp;
            //    //        } else {
            //    //            Fr_tt += Math.Pow(cf, 2) * 2.0 * cF_xx.Real / (double)numFp;
            //    //        }
            //    //    }
            //    //}


            //    // kappa = (radius^2 + 2*Fr_t^2 - radius*Fr_tt) / (r^2 + Fr_t^2)^3/2
            //    invDFT_coeff_curv[sp] = (Complex)((current_samplP[sp].Pow2() + 2.0 * Fr_t.Pow2() - current_samplP[sp] * Fr_tt) / (Math.Pow(current_samplP[sp].Pow2() + Fr_t.Pow2(), 3.0 / 2.0)));
            //}

            //Complex[] DFT_coeff_curv = DFT.NaiveForward(invDFT_coeff_curv, FourierOptions.Matlab);


            return (ScalarFunction)delegate (MultidimensionalArray nodes, MultidimensionalArray results) {
                //signal is real valued - i.e. for numSp sample points, numSp/2 coefficients are complex conjugated, thus they don't have to be evaluated individually

                //Complex cmplx_curv;
                for (int nd = 0; nd < nodes.GetLength(0); nd++) { // loop over nodes ...
                    double x = nodes[nd, 0] - center[0, 0];
                    double y = nodes[nd, 1] - center[0, 1];

                    double theta = Math.Atan2(y, x);
                    if (theta < 0) {
                        theta = Math.PI * 2 + theta;
                    };

                    //constant values
                    //results[nd] = DFT_coeff_curv[0].Real;
                    //if (cf_end == numFp / 2)
                    //    results[nd] += Complex.Multiply(DFT_coeff_curv[numFp / 2],
                    //                                           Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, theta * (numFp / 2)))
                    //                                           ).Real;

                    //for (int cf = 1; cf < cf_end; cf++) {
                    //    cmplx_curv = DFT_coeff_curv[cf] * Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, theta * cf));
                    //    results[nd] += 2.0 * cmplx_curv.Real;
                    //}

                    //results[nd] /= (double)numFp;

                    //results[nd] = 0.0;
                    //for (int cf = 0; cf < numFp / 2; cf++) {
                    //    if (DFT_coeff_curv[cf].Magnitude > significantMode) {
                    //        cmplx_curv = Complex.Multiply(DFT_coeff_curv[cf], Complex.Exp(Complex.ImaginaryOne * theta * cf));
                    //        if (cf == 0 || cf == numFp / 2) {
                    //            results[nd] += cmplx_curv.Real;
                    //        } else {
                    //            results[nd] += 2.0 * cmplx_curv.Real;
                    //        }
                    //    }
                    //}
                    //results[nd] /= (double)numFp;


                    double radius = 0.0;
                    if (curvComp_extended) {
                        //Complex RadiusComplex;
                        Complex cmplx_res;

                        //r(theta) for  k=0 and k= numSp/2
                        radius = DFT_coeff[0].Real;
                        //RadiusComplex = DFT_coeff[0];
                        if (cf_end == numFp / 2)
                            radius += Complex.Multiply(DFT_coeff[numFp / 2], Complex.Exp(Complex.ImaginaryOne * theta * (numFp / 2))).Real;
                        //radius += cmplx_res.Real;
                        //RadiusComplex += cmplx_res;
                        for (int cf = 1; cf < cf_end; cf++) {
                            cmplx_res = Complex.Multiply(DFT_coeff[cf], Complex.Exp(Complex.ImaginaryOne * theta * cf));
                            radius += 2.0 * cmplx_res.Real;  //-f(x) for  k=1... numSp/2-1 , values for k=numSp/2+1...k= numSp-1 are complex conjugated, thus the same real part
                            //RadiusComplex += 2 * cmplx_res;
                        }
                        radius /= (double)numFp;
                        //RadiusComplex /= (double)numFp;

                        //if (radius - RadiusComplex.Real >= 1e-14)
                        //    throw new ArithmeticException();
                    } else {
                        radius = x.Pow2() + y.Pow2();
                    }


                    Complex cF_x;
                    double Fr_t = 0.0;
                    Complex cF_xx;
                    double Fr_tt = 0.0;

                    // F_x
                    if (cf_end == numFp / 2) {
                        cF_x = Complex.ImaginaryOne * DFT_coeff[numFp / 2] * Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, theta * (numFp / 2)));
                        Fr_t = (double)(numFp / 2) * cF_x.Real / (double)numFp;
                    }

                    for (int cf = 1; cf < cf_end; cf++) {
                        cF_x = Complex.ImaginaryOne * DFT_coeff[cf] * Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, theta * cf));
                        Fr_t += cf * 2.0 * cF_x.Real / (double)numFp;
                    }

                    // F_xx
                    if (cf_end == numFp / 2) {
                        cF_xx = -DFT_coeff[numFp / 2] * Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, theta * (numFp / 2)));
                        Fr_tt = Math.Pow((double)(numFp / 2), 2) * cF_xx.Real / (double)numFp;
                    }

                    for (int cf = 1; cf < cf_end; cf++) {
                        cF_xx = -DFT_coeff[cf] * Complex.Exp(Complex.Multiply(Complex.ImaginaryOne, theta * cf));
                        Fr_tt += Math.Pow(cf, 2) * 2.0 * cF_xx.Real / (double)numFp;
                    }

                    //kappa = (radius ^ 2 + 2 * Fr_t ^ 2 - radius * Fr_tt) / (r ^ 2 + Fr_t ^ 2) ^ 3 / 2
                    results[nd] = (radius.Pow2() + 2.0 * Fr_t.Pow2() - radius * Fr_tt) / (Math.Pow(radius.Pow2() + Fr_t.Pow2(), 3.0 / 2.0));

                }
            };
        }



        /// <summary>
        /// returns the interpolated sample points as a string for the logfile entry
        /// </summary>
        /// <returns></returns>
        protected override string samplP_LogFormat() {
            string samplP_DB = "";
            // center
            samplP_DB = "\t" + center[0, 0].ToString() + "\t" + center[0, 1].ToString();
            // sample points
            for (int sp = 0; sp < numFp; sp++) {
                samplP_DB = samplP_DB + "\t" + current_samplP[sp].ToString();
            }

            return samplP_DB;
        }

        /// <summary>
        /// returns the material interface points as a string for the logfile entry
        /// </summary>
        /// <returns></returns>
        protected override string interfaceP_LogFormat() {
            string interP_DB = "";
            for (int ip = 0; ip < current_interfaceP.Lengths[0]; ip++) {
                interP_DB = interP_DB + "\t" + current_interfaceP[ip, 0].ToString();
                interP_DB = interP_DB + "\t" + current_interfaceP[ip, 1].ToString();
            }

            return interP_DB;
        }

        /// <summary>
        /// returns necessary values for restarting the Fourier Level set
        /// </summary>
        /// <returns></returns>
        public override double[] getRestartInfo() {

            double[] center_samplP = new double[numFp + 2];
            // center
            center_samplP[0] = center[0, 0];
            center_samplP[1] = center[0, 1];
            // sample points
            for (int i = 2; i < numFp + 2; i++) {
                center_samplP[i] = current_samplP[i - 2];
            }
            return center_samplP;
        }


    }

}
