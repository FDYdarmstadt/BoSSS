﻿using BoSSS.Foundation.IO;
using BoSSS.Solution.Statistic;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;

namespace BoSSS.Application.BoSSSpad {


    /// <summary>
    /// Workflow management.
    /// </summary>
    public partial class WorkflowMgm {

        TemporalConvergence m_tConvergence;


        /// <summary>
        /// Assistant to add temporal convergence data to the session table; 
        /// Not updated automatically, call <see cref="TemporalConvergence.Update"/> in order to re-evaluate errors
        /// </summary>
        public TemporalConvergence tConvergence {
            get {
                if(m_tConvergence == null)
                    m_tConvergence = new TemporalConvergence(this);
                return m_tConvergence;
            }
        }


        /// <summary>
        /// Automatic evaluation of spatial errors for all sessions in the current project 
        /// </summary>
        public class TemporalConvergence {


            WorkflowMgm owner;

            internal TemporalConvergence(WorkflowMgm __owner) {
                owner = __owner;
                //FieldsNormTypes = new Dictionary<string, NormType>();
            }

            /*
            /// <summary>
            /// Indicates for which fields the error should be computed in which norm
            /// </summary>
            public Dictionary<string, NormType> FieldsNormTypes {
                get;
                private set;
            }
            */

            static int GetGridSizeKey(ISessionInfo sess) {
                var Degs = sess.KeysAndQueries.Where(kv => kv.Key.StartsWith("Grid:NoOfCells")).OrderBy(kv => kv.Key);
                int a = Degs.Select(kv => Convert.ToInt32(kv.Value)).ToArray().First();
                return a;
            }

            /// <summary>
            /// Updates all columns related to convergence plots
            /// </summary>
            public void Update() {
                // Get all sessions which are successfully terminated
                // ==================================================
                var SuccSessions = owner.Sessions.Where(sess => sess.SuccessfulTermination == true).ToArray();

                // Group the sessions according to p Degree and grid size
                // =================================================
                System.Func<int, int, bool> eqFunc = (A, B) => ArrayTools.Equals(A,B);
                var comp = eqFunc.ToEqualityComparer();
                var SessionGroups = SuccSessions.GroupBy(GetGridSizeKey, comp).ToArray();

                // Temporal convergence for each session group
                // ==========================================

                // intermediate result storage
                // 1st key: Field name
                // 2nd key: session name
                // value: error norm
                var Errors = new Dictionary<string, Dictionary<Guid, double>>();


                foreach (IEnumerable<ISessionInfo> temporalSeries in SessionGroups) {
                    if(temporalSeries.Count() <= 1)
                        continue;

                    ITimestepInfo[] tsiS = temporalSeries.Select(sess => sess.Timesteps.Last()).ToArray();

                    // check that the termination time is the same
                    double t_end = tsiS[0].PhysicalTime;
                    foreach(var timeStep in tsiS) {
                        if (Math.Abs(timeStep.PhysicalTime - t_end) > 1e-6) throw new ApplicationException("temporal convergence is only valid for same timestep");
                    }

                    // find DG field identifications which are present in _all_ timesteps
                    var commonFieldIds = new HashSet<string>();
                    foreach(var fi in tsiS[0].FieldInitializers) {
                        string id = fi.Identification;

                        bool containedInOthers = true;
                        foreach(var tsi in tsiS.Skip(1)) {
                            if(tsi.FieldInitializers.Where(fii => fii.Identification == id).Count() <= 0)
                                containedInOthers = false;
                        }

                        if(containedInOthers) {
                            commonFieldIds.Add(id);
                        }
                    }

                    // check for which fields which norm should be computed
                    string[] Fields_L2emb;//, Fields_L2aprx, Fields_H1aprx;
                    {
                        //string[] Fields_L2_embedded_tmp = this.FieldsNormTypes.Where(kv => kv.Value == NormType.L2_embedded).Select(kv => kv.Key).ToArray();
                        //Fields_L2emb = commonFieldIds.Intersect(Fields_L2_embedded_tmp).ToArray();
                        Fields_L2emb = commonFieldIds.ToArray();

                        //string[] Fields_L2_approximate_tmp = this.FieldsNormTypes.Where(kv => kv.Value == NormType.L2_approximate).Select(kv => kv.Key).ToArray();
                        //Fields_L2aprx = commonFieldIds.Intersect(Fields_L2_approximate_tmp).ToArray();

                        //string[] Fields_H1_approximate_tmp = this.FieldsNormTypes.Where(kv => kv.Value == NormType.H1_approximate).Select(kv => kv.Key).ToArray();
                        //Fields_H1aprx = commonFieldIds.Intersect(Fields_H1_approximate_tmp).ToArray();
                    }

                    // compute L2-errors
                    //DGFieldComparison.ComputeErrors(Fields_L2emb, tsiS, out double[] hS_L2emb, out var DOFs_L2emb, out var ERRs_L2emb, out var tsiIdS_L2emb);
                    DGFieldComparisonEmbedded.ComputeErrorsTemporal(Fields_L2emb, tsiS, out double[] hS_L2emb, out var DOFs_L2emb, out var ERRs_L2emb, out var tsiIdS_L2emb);


                    // record errors
                    foreach (var id in Fields_L2emb) {

                        Dictionary<Guid, double> err_id;
                        if(!Errors.TryGetValue(id, out err_id)) {
                            err_id = new Dictionary<Guid, double>();
                            Errors.Add(id, err_id);
                        }

                        for(int iGrd = 0; iGrd < hS_L2emb.Length; iGrd++) {

                            ITimestepInfo tsi = tsiS.Single(t => t.ID == tsiIdS_L2emb[iGrd]);
                            ISessionInfo sess = tsi.Session;

                            err_id.Add(sess.ID, ERRs_L2emb[id][iGrd]);

                        }
                    }
                }

                // Set L2 error columns in session table
                // =====================================
                foreach(string fieldName in Errors.Keys) {
                    string colName = "L2Error_Temporal_" + fieldName;

                    if(owner.AdditionalSessionTableColums.ContainsKey(colName))
                        owner.AdditionalSessionTableColums.Remove(colName);

                    var ErrorsCol = Errors[fieldName];

                    owner.AdditionalSessionTableColums.Add(colName, delegate (ISessionInfo s) {
                        object ret = 0.0;

                        if(ErrorsCol.ContainsKey(s.ID))
                            ret = ErrorsCol[s.ID];
                        
                        return ret;
                    });

                }
            }
        }
    }
}
