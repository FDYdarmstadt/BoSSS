﻿using BoSSS.Foundation;
using BoSSS.Foundation.XDG;
using ilPSP.LinSolvers;
using System;
using System.Collections.Generic;
using ilPSP.Utils;
using ilPSP.LinSolvers.HYPRE;
using Newtonsoft.Json.Linq;
using ilPSP;
using BoSSS.Foundation.Grid.Classic;
using BoSSS.Platform;

namespace ApplicationWithIDT {


    public interface IOptProb {
        public void EvalObjective<V>(V obj_out, XDGField[] ConservativeFields) where V : IList<double>;
        public void EvalConstraint<V>(V cons_out, XDGField[] ConservativeFields) where V : IList<double>;
        public void EvalConsAndObj<V>(V obj_out, V cons_out, XDGField[] ConservativeFields) where V : IList<double>;
        public int GetObjLength(XDGField[] ConservativeFields);
        public MsrMatrix GetObjJacobian(XDGField[] ConservativeFields, LinearizationHint lHint);
        public MsrMatrix GetConsJacobian(XDGField[] ConservativeFields, LinearizationHint lHint);
        public XDGField[] CreateObjField(XDGField[] ConservativeFields);
        public (MsrMatrix,MsrMatrix) GetJacobians(XDGField[] ConservativeFields, LinearizationHint lHint);
        public XSpatialOperatorMk2 GetConsOperator();
    }


    public abstract class OptProbBase:IOptProb {

        public OptProbBase(XSpatialOperatorMk2 Op_cons, XSpatialOperatorMk2 Op_obj) {
            this.Op_cons = Op_cons;
            this.Op_obj = Op_obj;
        }

        /// <summary>
        /// This is the Spatial Operator from the PDE (Advection, Burgers, Euler,...)
        /// </summary>
        public XSpatialOperatorMk2 Op_cons { get; set; }

        /// <summary>
        /// This is the Spatial Operator defining the obj f = || R(U,\phi) ||, can be the same as Op_cons but can be different
        /// </summary>
        public XSpatialOperatorMk2 Op_obj { get; set; }
       
        /// <summary>
        /// Defines if Constraint can be evaluated from using the return from the objective
        /// </summary>
        public bool is_GetConsFromObj = false;

        public abstract void EvalObjective<V>(V obj_out, XDGField[] ConservativeFields) where V : IList<double>;
        public void EvalConstraint<V>(V cons_out, XDGField[] ConservativeFields) where V : IList<double> {
            var cons_map = new CoordinateMapping(ConservativeFields);
            var LsTrk = ConservativeFields[0].Basis.Tracker;
            cons_out.ScaleV(0.0);

            var Eval_r = Op_cons.GetEvaluatorEx(LsTrk, ConservativeFields, null, cons_map);
            Eval_r.Evaluate(1.0, 0.0, cons_out);
        }

        public abstract void GetConsFromObj<V>(V obj, V cons_out, XDGField[] ConservativeFields) where V : IList<double>;

        public abstract MsrMatrix GetConsFromObj(MsrMatrix obj,  XDGField[] ConservativeFields);
        /// <summary>
        /// Evaluates Constraint and Objective (without taking the norm)
        /// </summary>
        /// <typeparam name="V"></typeparam>
        /// <returns></returns>
        public void EvalConsAndObj<V>(V obj_out, V cons_out, XDGField[] ConservativeFields) where V : IList<double> {
            EvalObjective(obj_out, ConservativeFields);
            if(is_GetConsFromObj) {
                // Eval Objective 
                GetConsFromObj(obj_out, cons_out,ConservativeFields); 
            } else {
                EvalConstraint(cons_out, ConservativeFields);
            }
        }

        public abstract int GetObjLength(XDGField[] ConservativeFields);

        public XSpatialOperatorMk2 GetConsOperator() {
            return Op_cons;
        }

        public virtual MsrMatrix GetObjJacobian(XDGField[] ConservativeFields, LinearizationHint lHint) {
            var objFields = CreateObjField(ConservativeFields);
            var obj_map = new CoordinateMapping(objFields);


            int Length_obj = (int)obj_map.TotalLength;
            int length_con = (int)new CoordinateMapping(ConservativeFields).TotalLength;
            var obj_vec = new double[Length_obj];
            MsrMatrix Jobj = new MsrMatrix(Length_obj, length_con, 1, 1);
            switch(lHint) {
                case LinearizationHint.FDJacobi:

                var R_JacobianBuilder = Op_obj.GetFDJacobianBuilder(ConservativeFields, null, obj_map);
                R_JacobianBuilder.ComputeMatrix(Jobj, obj_vec);
                break;
                case LinearizationHint.GetJacobiOperator:

                var R_JacobiOperator = Op_obj.GetJacobiOperator(SpatialDimension: 2);
                var R_MatrixBuilder = R_JacobiOperator.GetMatrixBuilder(new CoordinateMapping(ConservativeFields), R_JacobiOperator.InvokeParameterFactory(ConservativeFields), obj_map);
                R_MatrixBuilder.ComputeMatrix(Jobj, obj_vec);
                break;
                case LinearizationHint.AdHoc:


                var r_JacobiOperator = Op_obj.GetJacobiOperator(SpatialDimension: 2);
                var R_MatrixBuilderAdhoc = Op_obj.GetMatrixBuilder(new CoordinateMapping(ConservativeFields), r_JacobiOperator.InvokeParameterFactory(ConservativeFields), obj_map);
                R_MatrixBuilderAdhoc.ComputeMatrix(Jobj, obj_vec);
                break;
                default:
                throw new ArgumentException("this should not happen");
            }
            return Jobj;

        }
        public MsrMatrix GetConsJacobian(XDGField[] ConservativeFields, LinearizationHint lHint) {
            var cons_map = new CoordinateMapping(ConservativeFields);
            var r = new double[cons_map.TotalLength];
            MsrMatrix Jobj = new MsrMatrix((int)cons_map.TotalLength, (int)cons_map.TotalLength, 1, 1);
            switch(lHint) {
                case LinearizationHint.FDJacobi:

                var JacobianBuilder = Op_cons.GetFDJacobianBuilder(ConservativeFields, null, cons_map);
                JacobianBuilder.ComputeMatrix(Jobj, r);
                break;
                case LinearizationHint.GetJacobiOperator:

                var JacobiOperator = Op_cons.GetJacobiOperator(SpatialDimension: 2);
                var MatrixBuilder = JacobiOperator.GetMatrixBuilder(cons_map, JacobiOperator.InvokeParameterFactory(ConservativeFields), cons_map);
                MatrixBuilder.ComputeMatrix(Jobj, r);
                break;
                case LinearizationHint.AdHoc:


                var r_JacobiOperator = Op_cons.GetJacobiOperator(SpatialDimension: 2);
                var MatrixBuilderAdhoc = Op_cons.GetMatrixBuilder(cons_map, r_JacobiOperator.InvokeParameterFactory(ConservativeFields), cons_map);
                MatrixBuilderAdhoc.ComputeMatrix(Jobj, r);
                break;
                default:
                throw new ArgumentException("this should not happen");
            }
            return Jobj;
        }
        public (MsrMatrix, MsrMatrix) GetJacobians(XDGField[] ConservativeFields, LinearizationHint lHint) {
            MsrMatrix Jobj = GetObjJacobian(ConservativeFields, lHint);
            MsrMatrix Jcons;
            if(is_GetConsFromObj) {
                Jcons = GetConsFromObj(Jobj, ConservativeFields);
            } else {
                Jcons = GetConsJacobian(ConservativeFields,lHint);
            }
            return (Jobj, Jcons);
        }

        public abstract XDGField[] CreateObjField(XDGField[] ConservativeFields);
    }
    public class SFFullEnRes : OptProbBase {

        int pDiff;
        public SFFullEnRes(XSpatialOperatorMk2 Op_cons,int pDiff)
                            : base(Op_cons,Op_cons) {
            this.pDiff = pDiff;
            is_GetConsFromObj = true;

        }

        public override XDGField[] CreateObjField(XDGField[] ConservativeFields) {
            var LsTrk = ConservativeFields[0].Basis.Tracker;
            var objFields = new XDGField[ConservativeFields.Length];
            for(int i = 0; i < ConservativeFields.Length; i++) {
                var enriched_basis = new XDGBasis(LsTrk, ConservativeFields[i].Basis.Degree + pDiff);
                objFields[i] = new XDGField(enriched_basis, "enres_" + ConservativeFields[i].Identification);
                //objFields[i].AccLaidBack(1, ConservativeFields[i]);
            }
            return objFields;
        }

        public override void EvalObjective<V>(V obj_out, XDGField[] ConservativeFields) {
            var objFields = CreateObjField(ConservativeFields);
            var obj_map = new CoordinateMapping(objFields);
            var LsTrk = ConservativeFields[0].Basis.Tracker;
            obj_out.ScaleV(0.0);
            var Eval_R = Op_cons.GetEvaluatorEx(LsTrk, ConservativeFields, null, obj_map);
            Eval_R.Evaluate(1.0, 0.0, obj_out);
            //if(CurrentAgglo > 0) {
            //    MultiphaseAgglomerator.ManipulateMatrixAndRHS(default(MsrMatrix), obj_out, obj_map, cons_map);
            //}
        }

        public override void GetConsFromObj<V>(V obj, V cons_out, XDGField[] ConservativeFields) {
            cons_out.ScaleV(0.0);
            var objFields = CreateObjField(ConservativeFields);
            var cons_map = new CoordinateMapping(ConservativeFields);
            var obj_map = new CoordinateMapping(objFields);
            var LsTrk = ConservativeFields[0].Basis.Tracker;

            int MaxModeR = obj_map.MaxTotalNoOfCoordinatesPerCell / ConservativeFields.Length / LsTrk.TotalNoOfSpecies;
            int MaxModer = cons_map.MaxTotalNoOfCoordinatesPerCell / ConservativeFields.Length / LsTrk.TotalNoOfSpecies;


            if(obj.Count != obj_map.TotalLength ||cons_out.Count != cons_map.TotalLength) {
                throw new Exception("enRes or res array has wrong length");
            }
            //double MaxMode2_double = (double) MaxMode2;
            //var Jr_copy = new MsrMatrix(length_r, length_r + length_l);
            int iField;
            int jCell;
            int nMode;

            for(int iRow = 0; iRow < (int)cons_map.TotalLength; iRow++) {
                cons_map.LocalFieldCoordinateIndex(iRow, out iField, out jCell, out nMode);
                double rMode = (double)nMode;
                int fac = (int)Math.Floor(rMode / MaxModer);
                int row = obj_map.LocalUniqueCoordinateIndex(iField, jCell, nMode + fac * (MaxModeR - MaxModer));
                cons_out[iRow] = obj[row];
            }
        }

        public override int GetObjLength(XDGField[] ConservativeFields) {
            var objFields = CreateObjField(ConservativeFields);
            return (int) new CoordinateVector(objFields).Count;
        }
        public override MsrMatrix GetConsFromObj(MsrMatrix Jobj, XDGField[] ConservativeFields) {
            
            var cons_map = new CoordinateMapping(ConservativeFields);
            var objFields = CreateObjField(ConservativeFields);
            var obj_map = new CoordinateMapping(objFields);
            var LsTrk = ConservativeFields[0].Basis.Tracker;

            int MaxModeR = obj_map.MaxTotalNoOfCoordinatesPerCell / ConservativeFields.Length / LsTrk.TotalNoOfSpecies;
            int MaxModer = cons_map.MaxTotalNoOfCoordinatesPerCell / ConservativeFields.Length / LsTrk.TotalNoOfSpecies;

            //double MaxMode2_double = (double) MaxMode2;
            //var Jr_copy = new MsrMatrix(length_r, length_r + length_l);
            MsrMatrix Jcons = new MsrMatrix((int)cons_map.TotalLength, (int)cons_map.TotalLength, 1, 1);
            int iField = 0;
            int jCell = 0;
            int nMode = 0;

            for(int iRow = 0; iRow < cons_map.TotalLength; iRow++) {
                cons_map.LocalFieldCoordinateIndex(iRow, out iField, out jCell, out nMode);
                double rMode = (double)nMode;

                //int row = obj_f_map.LocalUniqueCoordinateIndex(iField, jCell, nMode);

                int fac = (int)Math.Floor(rMode / MaxModer);
                int row = obj_map.LocalUniqueCoordinateIndex(iField, jCell, nMode + fac * (MaxModeR - MaxModer));
                //int row = obj_f_map.LocalUniqueCoordinateIndex(iField, jCell, nMode );
                //Console.WriteLine(iRow + ", " +row);

                Jcons.SetRow(iRow, Jobj.GetRowShallow(row).CloneAs());
            }
            return Jcons;
        }

    }
    public class SFNearBandEnRes : SFFullEnRes {
        public SFNearBandEnRes(XSpatialOperatorMk2 Op_cons, int pDiff) : base(Op_cons, pDiff) {
            is_GetConsFromObj = false;
        }
        public override void EvalObjective<V>(V obj_out, XDGField[] ConservativeFields) {
            base.EvalObjective(obj_out, ConservativeFields);
            var tracker = ConservativeFields[0].Basis.Tracker;
            var objFields = CreateObjField(ConservativeFields);
            var obj_map = new CoordinateMapping(objFields);
            var NB = tracker.Regions.GetNearFieldMask(1);
            //very ugly, because of contains method
            for(int iRow = 0; iRow < (int)obj_map.TotalLength; iRow++) {
                obj_map.LocalFieldCoordinateIndex(iRow, out int iField, out int jCell, out int nMode);
                if(NB.Contains(jCell)) {
                    continue; // do nothing if cell is in near-band
                } else {
                    obj_out[iRow] = 0; //set unneeded entries to zero
                }
            }
        }
        public override MsrMatrix GetObjJacobian(XDGField[] ConservativeFields, LinearizationHint lHint) {
            var Jobj = base.GetObjJacobian(ConservativeFields, lHint);
            var tracker = ConservativeFields[0].Basis.Tracker;
            var objFields = CreateObjField(ConservativeFields);
            var obj_map = new CoordinateMapping(objFields);
            var NB = tracker.Regions.GetNearFieldMask(1);
            //very ugly, because of contains method
            for(int iRow = 0; iRow < (int)obj_map.TotalLength; iRow++) {
                obj_map.LocalFieldCoordinateIndex(iRow, out int iField, out int jCell, out int nMode);
                if(NB.Contains(jCell)) {
                    continue; // do nothing if cell is in near-band
                } else {
                    Jobj.ClearRow(iRow); //set unneeded entries to zero
                }
            }
            return Jobj;
        }
    }
    public class SFCutCellEnRes : SFFullEnRes {
        public SFCutCellEnRes(XSpatialOperatorMk2 Op_cons, int pDiff) : base(Op_cons, pDiff) {
            is_GetConsFromObj = false;
        }
        public override void EvalObjective<V>(V obj_out, XDGField[] ConservativeFields) {
            base.EvalObjective(obj_out, ConservativeFields);
            var tracker = ConservativeFields[0].Basis.Tracker;
            var objFields = CreateObjField(ConservativeFields);
            var obj_map = new CoordinateMapping(objFields);
            var NB = tracker.Regions.GetCutCellMask();
            //very ugly, because of contains method
            for(int iRow = 0; iRow < (int)obj_map.TotalLength; iRow++) {
                obj_map.LocalFieldCoordinateIndex(iRow, out int iField, out int jCell, out int nMode);
                if(NB.Contains(jCell)) {
                    continue; // do nothing if cell is in near-band
                } else {
                    obj_out[iRow] = 0; //set unneeded entries to zero
                }
            }
        }
        public override MsrMatrix GetObjJacobian(XDGField[] ConservativeFields, LinearizationHint lHint) {
            var Jobj = base.GetObjJacobian(ConservativeFields, lHint);
            var tracker = ConservativeFields[0].Basis.Tracker;
            var objFields = CreateObjField(ConservativeFields);
            var obj_map = new CoordinateMapping(objFields);
            var NB = tracker.Regions.GetCutCellMask();
            //very ugly, because of contains method
            for(int iRow = 0; iRow < (int)obj_map.TotalLength; iRow++) {
                obj_map.LocalFieldCoordinateIndex(iRow, out int iField, out int jCell, out int nMode);
                if(NB.Contains(jCell)) {
                    continue; // do nothing if cell is in near-band
                } else {
                    Jobj.ClearRow(iRow); //set unneeded entries to zero
                }
            }
            return Jobj;
        }
    }

    public class SFRankineHugoniotBase : OptProbBase {
        int pDiff;
        public SFRankineHugoniotBase(XSpatialOperatorMk2 Op_cons, XSpatialOperatorMk2 Op_obj, int pDiff)
                            : base(Op_cons, Op_obj) {
            is_GetConsFromObj = false;
            this.pDiff = pDiff;
        }
        public override XDGField[] CreateObjField(XDGField[] ConservativeFields) {
            var LsTrk = ConservativeFields[0].Basis.Tracker;
            var objFields = new XDGField[ConservativeFields.Length];
            for(int i = 0; i < ConservativeFields.Length; i++) {
                var basis = new XDGBasis(LsTrk, ConservativeFields[i].Basis.Degree + pDiff);
                objFields[i] = new XDGField(basis, "RH_" + ConservativeFields[i].Identification);
                //objFields[i].AccLaidBack(1, ConservativeFields[i]);
            }
            return objFields;
        }

        public override void EvalObjective<V>(V obj_out, XDGField[] ConservativeFields) {
            var objFields = CreateObjField(ConservativeFields);
            var obj_map = new CoordinateMapping(objFields);
            var LsTrk = ConservativeFields[0].Basis.Tracker;
            obj_out.ScaleV(0.0);
            var Eval_R = Op_obj.GetEvaluatorEx(LsTrk, ConservativeFields, null, obj_map);
            Eval_R.Evaluate(1.0, 0.0, obj_out);
        }

        public override void GetConsFromObj<V>(V obj, V cons_out, XDGField[] ConservativeFields) {
            throw new NotImplementedException();
        }

        public override MsrMatrix GetConsFromObj(MsrMatrix obj, XDGField[] ConservativeFields) {
            throw new NotImplementedException();
        }

        public override int GetObjLength(XDGField[] ConservativeFields) {
            var objFields = CreateObjField(ConservativeFields);
            return new CoordinateVector(objFields).Count;
        }
    }
  
}
