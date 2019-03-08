using BoSSS.Foundation.IO;
using BoSSS.Solution.Statistic;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Encodes in which norms convergence is tested (<see cref="WorkflowMgm.SpatialConvergence"/>)
    /// </summary>
    public enum NormType {

        /// <summary>
        /// Norm by <see cref="DGFieldComparison.ComputeErrors"/>; very accurate, but requires geometrically embedded meshes
        /// </summary>
        L2_embedded,

        L2_approximate,

        H1_approximate
    }


    /// <summary>
    /// Workflow management.
    /// </summary>
    public partial class WorkflowMgm {

        SpatialConvergence m_hConvergence;


        /// <summary>
        /// Assistant to add spatial convergence data to the session table; 
        /// Not updated automatically, call <see cref="SpatialConvergence.Update"/> in order to re-evaluate errors
        /// </summary>
        public SpatialConvergence hConvergence {
            get {
                if(m_hConvergence == null)
                    m_hConvergence = new SpatialConvergence(this);
                return m_hConvergence;
            }
        }


        /// <summary>
        /// Automatic evaluation of spatial errors for all sessions in the current project 
        /// </summary>
        public class SpatialConvergence {


            WorkflowMgm owner;

            internal SpatialConvergence(WorkflowMgm __owner) {
                owner = __owner;
                FieldsNormTypes = new Dictionary<string, NormType>();
            }

            /// <summary>
            /// Indicates for which fields the error should be computed in which norm
            /// </summary>
            public Dictionary<string, NormType> FieldsNormTypes {
                get;
                private set;
            }


            
            static int[] GetDGDegreeKey(ISessionInfo sess) {
                var Degs = sess.KeysAndQueries.Where(kv => kv.Key.StartsWith("DGdegree:")).OrderBy(kv => kv.Key);
                int[] a = Degs.Select(kv => Convert.ToInt32(kv.Value)).ToArray();
                return a;
            }

            /// <summary>
            /// Updates all columns related to convergence plots
            /// </summary>
            public void Update() {
                // Get all sessions which are successfully terminated
                // ==================================================
                var SuccSessions = owner.Sessions.Where(sess => sess.SuccessfulTermination == true).ToArray();

                // Group the sessions according to polynomial degree 
                // =================================================
                System.Func<int[],int[],bool> eqFunc = (A,B) => ArrayTools.AreEqual(A,B);
                var comp = eqFunc.ToEqualityComparer();
                var SessionGroups = SuccSessions.GroupBy(GetDGDegreeKey, comp).ToArray();


                // Spatial convergence for each session group
                // ==========================================

                // intermediate result storage
                // 1st key: Field name
                // 2nd key: session name
                // value: error norm
                var Errors = new Dictionary<string, Dictionary<Guid, double>>();


                foreach (IEnumerable<ISessionInfo> spatialSeries in SessionGroups) {
                    if(spatialSeries.Count() <= 1)
                        continue;

                    ITimestepInfo[] tsiS = spatialSeries.Select(sess => sess.Timesteps.Last()).ToArray();

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

                    string[] fieldIds = commonFieldIds.Intersect()

                    // compute L2-errors
                    DGFieldComparison.ComputeErrors(fieldIds, tsiS, out double[] hS, out var DOFs, out var ERRs, out var tsiIdS);





                    // record errors
                    foreach(var id in fieldIds) {

                        Dictionary<Guid, double> err_id;
                        if(!Errors.TryGetValue(id, out err_id)) {
                            err_id = new Dictionary<Guid, double>();
                            Errors.Add(id, err_id);
                        }

                        for(int iGrd = 0; iGrd < hS.Length; iGrd++) {

                            ITimestepInfo tsi = tsiS.Single(t => t.ID == tsiIdS[iGrd]);
                            ISessionInfo sess = tsi.Session;

                            err_id.Add(sess.ID, ERRs[id][iGrd]);

                        }
                    }
                }




                // Set L2 error columns in session table
                // =====================================
                foreach(string fieldName in Errors.Keys) {
                    string colName = "L2Error_" + fieldName;

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
