using BoSSS.Foundation.IO;
using BoSSS.Solution.Statistic;
using ilPSP;
using ilPSP.Utils;
using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace BoSSS.Application.BoSSSpad
{




    /// <summary>
    /// Workflow management.
    /// </summary>
    public partial class WorkflowMgm
    {

        ConditionNumberReader m_condReader;


        /// <summary>
        /// Assistant to add condition number data to the session table; 
        /// for this CondLogger should be already added as Insitu post processing tool
        /// Not updated automatically, call <see cref="ConditionNumberReader.Update"/> in order to re-evaluate errors
        /// </summary>
        public ConditionNumberReader condReader
        {
            get
            {
                if (m_condReader == null)
                    m_condReader = new ConditionNumberReader(this);
                return m_condReader;
            }
        }


        /// <summary>
        /// Automatic condition number data for all sessions in the current project 
        /// </summary>
        public class ConditionNumberReader
        {

            WorkflowMgm owner;

            internal ConditionNumberReader(WorkflowMgm __owner)
            {
                owner = __owner;
            }


            /// <summary>
            /// Returns the all column names associated with condition number tables
            /// </summary>
            /// <param name="sessions"></param>
            public string[] GetColumnNames(ISessionInfo[] sessions = null){
                string[] allColumnNames = new string[] { "" };

                if (sessions is null)
                    sessions = owner.Sessions.Where(sess => sess.SuccessfulTermination == true).ToArray();


                int numberSessions = sessions.Count();
                for (int j = 0; j < numberSessions; j++){
                    ISessionInfo currentSession = sessions.Pick(j);


                    Console.WriteLine("Session: {0}", currentSession.ID);
                    string path = @currentSession.Database.Path + "\\sessions\\" + currentSession.ID + "\\CondNumbers.txt";
                    string[] lines;
                    string header;

                    try
                    { //reading
                        lines = File.ReadAllLines(path);
                        header = lines[0];
                    } catch
                    {
                        continue;
                    }

                    //Get column names
                    var columnsNames = header.Split(new string[] { "\t" }, StringSplitOptions.RemoveEmptyEntries);
                    if (!allColumnNames.SequenceEqual(columnsNames))
                    {
                        allColumnNames = columnsNames;
                        allColumnNames.ForEach(c => Console.Write(c + ", "));
                        Console.WriteLine("");
                    }
                }
                return allColumnNames;
            }


            /// <summary>
            /// Updates all columns related to convergence plots
            /// </summary>
            /// <param name="sessions">list of sessions for the study (if null, only the successfully terminated simulations)</param>
            /// <param name="operation">which operation to perform in the list (if null, Enumerable.Max)</param>
            /// <param name="columnNames">which colums to be extracted (if null, all the possible columns)</param> 
            /// <param name="marker">string marker in the session table (columnname + marker)</param>
            /// <remarks>Be aware that column names can vary for different dimensions.</remarks>
            public void Update(ISessionInfo[] sessions = null, string[] columnNames = null, Func<List<double>, double> operation = null, string marker = "") {
                // Get all sessions which are successfully terminated
                // ==================================================
                 sessions ??= owner.Sessions.Where(sess => sess.SuccessfulTermination == true).ToArray();

                // Check if any specific column is provided. (Be aware that column names can vary for different dimensions.)
                columnNames ??= GetColumnNames(sessions);

                // If no operation is provided, use Enumerable.Max
                operation ??= Enumerable.Max;

                // Get the table for the sessions
                var condTable = sessions.CheckForCondLogging();

                // Set columns in session table
                // =====================================
                foreach (string column in columnNames)
                {
                    string colName = column + marker;

                    if (owner.AdditionalSessionTableColums.ContainsKey(colName))
                        owner.AdditionalSessionTableColums.Remove(colName);


                    owner.AdditionalSessionTableColums.Add(colName, delegate (ISessionInfo s) {
                        object ret = 0.0;

                        if (condTable.ContainsKey(s.ID))
                            ret = operation(condTable[s.ID][colName]);

                        return ret;
                    });

                }
            }
        }
    }
}
