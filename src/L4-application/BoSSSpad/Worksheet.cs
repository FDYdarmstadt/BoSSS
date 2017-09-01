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
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Drawing;
using System.IO;
using System.Threading;
using System.Windows.Forms;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// GUI logic for the Worksheet-DBE.
    /// </summary>
    public class Worksheet : Form {

        /// <summary>
        /// font for rendering the commands entered by the user
        /// </summary>
        internal Font CommandFont;

        /// <summary>
        /// font for rendering the answer of the REPL
        /// </summary>
        internal Font ResultFont;

        /// <summary>
        /// due to a resize, the line-wrap may change the text in a box;
        /// setting this to true suppresses the action fired by the text-changed event.
        /// </summary>
        internal bool BlockTextChanged = false;

        /// <summary>
        /// if true, the <see cref="MyclientSizeChanged"/> event handler has no effect. 
        /// </summary>
        internal bool BlockResizeClient = false;

        /// <summary>
        /// the document
        /// </summary>
        WorksheetEntry CommandsListHead;

        internal WorksheetEntry CurrentCommand;


        internal Panel ClientRectPanel;

        internal Panel DocumentPanel;


        private System.Windows.Forms.MenuStrip menuStrip1;
        private System.Windows.Forms.ToolStripMenuItem MenuItem_file;
        private System.Windows.Forms.ToolStripMenuItem MenuItem_new;
        private System.Windows.Forms.ToolStripMenuItem MenuItem_open;
        private System.Windows.Forms.ToolStripMenuItem MenuItem_save;
        private System.Windows.Forms.ToolStripMenuItem MenuItem_saveAs;
        private System.Windows.Forms.ToolStripMenuItem MenuItem_edit;
        private System.Windows.Forms.ToolStripMenuItem MenuItem_commands;
        private System.Windows.Forms.ToolStripMenuItem MenuItem_insertAfter;
        private System.Windows.Forms.ToolStripMenuItem MenuItem_delete;
        private System.Windows.Forms.ToolStripMenuItem MenuItem_clearOutput;
        private System.Windows.Forms.ToolStripMenuItem MenuItem_clearAllOutput;
        private System.Windows.Forms.ToolStripMenuItem MenuItem_execute;
        private System.Windows.Forms.ToolStripMenuItem MenuItem_executeFromHere;
        private System.Windows.Forms.ToolStripMenuItem MenuItem_executeUntilHere;
        private System.Windows.Forms.ToolStripMenuItem MenuItem_executeWorksheet;
        private System.Windows.Forms.ToolStripMenuItem MenuItem_interruptCurrent;
        private System.Windows.Forms.ToolStripMenuItem MenuItem_unqueuePending;
        private System.Windows.Forms.ToolStripMenuItem MenuItem_autocomplete;

        public Worksheet(string fileToOpen = null) {

            this.components = new System.ComponentModel.Container();


            // Menu
            // ====
            {
                this.menuStrip1 = new System.Windows.Forms.MenuStrip();
                this.MenuItem_file = new System.Windows.Forms.ToolStripMenuItem();
                this.MenuItem_new = new System.Windows.Forms.ToolStripMenuItem();
                this.MenuItem_open = new System.Windows.Forms.ToolStripMenuItem();
                this.MenuItem_save = new System.Windows.Forms.ToolStripMenuItem();
                this.MenuItem_saveAs = new System.Windows.Forms.ToolStripMenuItem();
                this.MenuItem_edit = new System.Windows.Forms.ToolStripMenuItem();
                this.MenuItem_commands = new System.Windows.Forms.ToolStripMenuItem();
                this.MenuItem_insertAfter = new System.Windows.Forms.ToolStripMenuItem();
                this.MenuItem_delete = new System.Windows.Forms.ToolStripMenuItem();
                this.MenuItem_clearOutput = new System.Windows.Forms.ToolStripMenuItem();
                this.MenuItem_clearAllOutput = new System.Windows.Forms.ToolStripMenuItem();
                this.MenuItem_execute = new System.Windows.Forms.ToolStripMenuItem();
                this.MenuItem_executeFromHere = new System.Windows.Forms.ToolStripMenuItem();
                this.MenuItem_executeUntilHere = new System.Windows.Forms.ToolStripMenuItem();
                this.MenuItem_executeWorksheet = new System.Windows.Forms.ToolStripMenuItem();
                this.MenuItem_interruptCurrent = new System.Windows.Forms.ToolStripMenuItem();
                this.MenuItem_unqueuePending = new System.Windows.Forms.ToolStripMenuItem();
                this.MenuItem_autocomplete = new System.Windows.Forms.ToolStripMenuItem();

                this.menuStrip1.SuspendLayout();
                this.SuspendLayout();
                // 
                // menuStrip1
                // 
                this.menuStrip1.Items.AddRange(new System.Windows.Forms.ToolStripItem[] {
                    this.MenuItem_file,
                    this.MenuItem_edit,
                    this.MenuItem_commands });
                //this.MenuItem_texLink});
                this.menuStrip1.Location = new System.Drawing.Point(0, 0);
                this.menuStrip1.Name = "menuStrip1";
                this.menuStrip1.Size = new System.Drawing.Size(590, 24);
                this.menuStrip1.TabIndex = 0;
                this.menuStrip1.Text = "menuStrip1";
                // 
                // MenuItem_file
                // 
                this.MenuItem_file.DropDownItems.AddRange(new System.Windows.Forms.ToolStripItem[] {
                    this.MenuItem_new,
                    this.MenuItem_open,
                    this.MenuItem_save,
                    this.MenuItem_saveAs});
                this.MenuItem_file.Name = "MenuItem_file";
                this.MenuItem_file.Size = new System.Drawing.Size(37, 20);
                this.MenuItem_file.Text = "File";
                // 
                // MenuItem_new
                // 
                this.MenuItem_new.Name = "MenuItem_new";
                this.MenuItem_new.ShortcutKeys = ((System.Windows.Forms.Keys)((System.Windows.Forms.Keys.Control | System.Windows.Forms.Keys.N)));
                this.MenuItem_new.Size = new System.Drawing.Size(152, 22);
                this.MenuItem_new.Text = "New";
                this.MenuItem_new.Click += new System.EventHandler(this.MenuItem_new_Click);
                // 
                // MenuItem_open
                // 
                this.MenuItem_open.Name = "MenuItem_open";
                this.MenuItem_open.ShortcutKeys = ((System.Windows.Forms.Keys)((System.Windows.Forms.Keys.Control | System.Windows.Forms.Keys.O)));
                this.MenuItem_open.Size = new System.Drawing.Size(152, 22);
                this.MenuItem_open.Text = "Open";
                this.MenuItem_open.Click += new System.EventHandler(this.MenuItem_open_Click);
                // 
                // MenuItem_save
                // 
                this.MenuItem_save.Name = "MenuItem_save";
                this.MenuItem_save.ShortcutKeys = ((System.Windows.Forms.Keys)((System.Windows.Forms.Keys.Control | System.Windows.Forms.Keys.S)));
                this.MenuItem_save.Size = new System.Drawing.Size(152, 22);
                this.MenuItem_save.Text = "Save";
                this.MenuItem_save.Click += new System.EventHandler(this.MenuItem_save_Click);
                // 
                // saveAsToolStripMenuItem
                // 
                this.MenuItem_saveAs.Name = "saveAsToolStripMenuItem";
                this.MenuItem_saveAs.Size = new System.Drawing.Size(152, 22);
                this.MenuItem_saveAs.Text = "Save As...";
                this.MenuItem_saveAs.Click += new System.EventHandler(this.MenuItem_saveAs_Click);
                // 
                // MenuItem_edit
                // 
                this.MenuItem_edit.Name = "MenuItem_edit";
                this.MenuItem_edit.Size = new System.Drawing.Size(39, 20);
                this.MenuItem_edit.Text = "Edit";
                this.MenuItem_edit.DropDownItems.AddRange(new ToolStripItem[] {
                    this.MenuItem_delete,
                    this.MenuItem_clearOutput,
                    this.MenuItem_clearAllOutput});
                // 
                // MenuItem_commands
                // 
                this.MenuItem_commands.DropDownItems.AddRange(new ToolStripItem[] {
                    this.MenuItem_insertAfter,
                    this.MenuItem_execute,
                    this.MenuItem_executeFromHere,
                    this.MenuItem_executeUntilHere,
                    this.MenuItem_executeWorksheet,
                    this.MenuItem_interruptCurrent,
                    this.MenuItem_unqueuePending,
                    this.MenuItem_autocomplete});
                this.MenuItem_commands.Name = "MenuItem_commands";
                this.MenuItem_commands.Size = new System.Drawing.Size(81, 20);
                this.MenuItem_commands.Text = "Commands";

                // 
                // MenuItem_insertAfter
                // 
                this.MenuItem_insertAfter.Name = "MenuItem_insertAfter";
                this.MenuItem_insertAfter.ShortcutKeys = ((System.Windows.Forms.Keys)((System.Windows.Forms.Keys.Control | System.Windows.Forms.Keys.I)));
                this.MenuItem_insertAfter.Size = new System.Drawing.Size(173, 22);
                this.MenuItem_insertAfter.Text = "Insert after";
                this.MenuItem_insertAfter.Click += new System.EventHandler(this.MenuItem_insertAfter_Click);
                // 
                // MenuItem_delete
                // 
                this.MenuItem_delete.Name = "MenuItem_delete";
                this.MenuItem_delete.ShortcutKeys = ((System.Windows.Forms.Keys)((System.Windows.Forms.Keys.Control | System.Windows.Forms.Keys.Delete)));
                this.MenuItem_delete.Size = new System.Drawing.Size(173, 22);
                this.MenuItem_delete.Text = "Delete command";
                this.MenuItem_delete.Click += new System.EventHandler(this.MenuItem_delete_Click);
                // 
                // MenuItem_clearOutput
                // 
                this.MenuItem_clearOutput.Name = "MenuItem_clear";
                this.MenuItem_clearOutput.Size = new System.Drawing.Size(173, 22);
                this.MenuItem_clearOutput.Text = "Delete output of current cell";
                this.MenuItem_clearOutput.Click += new System.EventHandler(this.MenuItem_clearOutput_Click);
                // 
                // MenuItem_clearAllOutput
                // 
                this.MenuItem_clearAllOutput.Name = "MenuItem_clearAll";
                this.MenuItem_clearAllOutput.Size = new System.Drawing.Size(173, 22);
                this.MenuItem_clearAllOutput.Text = "Delete output of all cells";
                this.MenuItem_clearAllOutput.Click += new System.EventHandler(this.MenuItem_clearAllOutput_Click);
                // 
                // MenuItem_execute
                // 
                this.MenuItem_execute.Name = "MenuItem_execute";
                this.MenuItem_execute.Size = new System.Drawing.Size(173, 22);
                this.MenuItem_execute.Text = "Execute selected entry";
                this.MenuItem_execute.Click += new System.EventHandler(this.MenuItem_execute_Click);
                // 
                // MenuItem_executeFromHere
                // 
                this.MenuItem_executeFromHere.Name = "MenuItem_executeFromHere";
                this.MenuItem_executeFromHere.Size = new System.Drawing.Size(173, 22);
                this.MenuItem_executeFromHere.ShortcutKeys = ((System.Windows.Forms.Keys)((System.Windows.Forms.Keys.Control | System.Windows.Forms.Keys.F5)));
                this.MenuItem_executeFromHere.Text = "Execute from here";
                this.MenuItem_executeFromHere.Click += new System.EventHandler(this.MenuItem_executeFromHere_Click);
                // 
                // MenuItem_executeUntilHere
                // 
                this.MenuItem_executeUntilHere.Name = "MenuItem_executeUntilHere";
                this.MenuItem_executeUntilHere.Size = new System.Drawing.Size(173, 22);
                this.MenuItem_executeUntilHere.Text = "Execute until here";
                this.MenuItem_executeUntilHere.Click += new System.EventHandler(this.MenuItem_executeUntilHere_Click);
                // 
                // MenuItem_executeWorksheet
                // 
                this.MenuItem_executeWorksheet.Name = "MenuItem_executeWorksheet";
                this.MenuItem_executeWorksheet.Size = new System.Drawing.Size(173, 22);
                this.MenuItem_executeWorksheet.ShortcutKeys = System.Windows.Forms.Keys.F5;
                this.MenuItem_executeWorksheet.Text = "Execute entire worksheet";
                this.MenuItem_executeWorksheet.Click += new System.EventHandler(this.MenuItem_executeWorksheet_Click);
                // 
                // MenuItem_interruptCurrent
                // 
                this.MenuItem_interruptCurrent.Name = "MenuItem_executeWorksheet";
                this.MenuItem_interruptCurrent.Size = new System.Drawing.Size(173, 22);
                this.MenuItem_interruptCurrent.Text = "Interrupt current command";
                this.MenuItem_interruptCurrent.Click += new System.EventHandler(this.MenuItem_interruptCurrent_Click);
                // 
                // MenuItem_unqueuePending
                // 
                this.MenuItem_unqueuePending.Name = "MenuItem_unqueuePending";
                this.MenuItem_unqueuePending.Size = new System.Drawing.Size(173, 22);
                this.MenuItem_unqueuePending.Text = "De-Queue all pending commands";
                this.MenuItem_unqueuePending.Click += new System.EventHandler(this.MenuItem_unqueuePending_Click);
                // 
                // MenuItem_unqueuePending
                // 
                this.MenuItem_autocomplete.Name = "MenuItem_autocomplete";
                this.MenuItem_autocomplete.Size = new System.Drawing.Size(173, 22);
                this.MenuItem_autocomplete.Text = "Auto-Completion";
                this.MenuItem_autocomplete.CheckOnClick = true;
                this.MenuItem_autocomplete.CheckState = CheckState.Checked;
                // 
                // 
                // Form1
                // 
                this.AutoScaleDimensions = new System.Drawing.SizeF(6F, 13F);
                this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
                this.ClientSize = new System.Drawing.Size(590, 462);
                this.Controls.Add(this.menuStrip1);
                this.MainMenuStrip = this.menuStrip1;
                this.Name = "Form1";
                this.Text = "Form1";
                this.menuStrip1.ResumeLayout(false);
                this.menuStrip1.PerformLayout();
                this.ResumeLayout(false);
                this.PerformLayout();
            }

            // "document" region
            // =================
            {
                this.MinimumSize = new Size(320, 200);
                this.Size = new Size(900, 600);

                this.ClientRectPanel = new Panel();
                //this.ClientRectPanel = new PanelNoScrollOnFocus();
                this.Controls.Add(ClientRectPanel);
                this.ClientRectPanel.Top = this.menuStrip1.Bottom + 1;
                this.ClientRectPanel.Left = 0;
                this.ClientRectPanel.Width = this.ClientSize.Width - 1;
                this.ClientRectPanel.Height = this.ClientSize.Height - 2 - this.ClientRectPanel.Top;
                this.ClientRectPanel.AutoScroll = true;
                //this.ClientRectPanel.PaddingChanged   += ClientRectPanel_Region;
                this.ClientRectPanel.Show();
                //this.ClientRectPanel.Scroll += ClientRectPanel_scroll;

                //this.ClientRectPanel.LocationChanged += ClientRectPanel_Region;

                this.DocumentPanel = new Panel();
                this.DocumentPanel.Top = 0; // this.menuStrip1.Bottom + 1;
                this.DocumentPanel.Left = 0;
                DocumentPanel.Width = this.ClientSize.Width - 1;

                this.ClientRectPanel.Controls.Add(this.DocumentPanel);
                this.DocumentPanel.Show();


                this.CommandFont = new Font("Courier New", 12.0F, FontStyle.Bold, GraphicsUnit.Point, ((byte)(0))); //new Font("Courier New Bold", 15);
                this.ResultFont = new Font("Courier New", 12.0F, FontStyle.Bold, GraphicsUnit.Point, ((byte)(0))); //new Font("Courier", 10);

                this.ClientSizeChanged += MyclientSizeChanged;
                this.FormClosing += new System.Windows.Forms.FormClosingEventHandler(this.FormClosingEvent);
            }


            // misc
            // ====

            {
                this.AutosaveTimer = new System.Windows.Forms.Timer(this.components);
                this.AutosaveTimer.Interval = 1000 * 60; // jede minute
                this.AutosaveTimer.Tick += new System.EventHandler(this.AutosaveTimer_Tick);
                this.AutosaveTimer.Start();
            }

            {
                // Displays a SaveFileDialog so the user can save the Image
                // assigned to Button2.
                saveFileDialog1 = new SaveFileDialog();
                saveFileDialog1.Filter = "BoSSSPad Worksheet|*.bws|TeX files (*.tex)|*.tex|All files (*.*)|*.*";
                saveFileDialog1.Title = "Save a BoSSSPad Worksheet File";
                saveFileDialog1.OverwritePrompt = true;
                saveFileDialog1.AddExtension = true;

                openFileDialog1 = new OpenFileDialog();
                openFileDialog1.CheckFileExists = true;
                openFileDialog1.Multiselect = false;
                openFileDialog1.Title = "Open a BoSSSPad Worksheet File";
                openFileDialog1.Filter = "Supported files|*.bws;*.tex|BoSSSPad Worksheet (*.bws)|*.bws|TeX files (*.tex)|*.tex|All files (*.*)|*.*";
                openFileDialog1.AddExtension = true;

                selectTexFileDialog = new OpenFileDialog();
                selectTexFileDialog.CheckFileExists = true;
                selectTexFileDialog.Multiselect = false;
                selectTexFileDialog.Title = "Select LaTeX/TeX file for TeXlink";
                selectTexFileDialog.Filter = "LaTeX/TeX - files|*.tex";
                selectTexFileDialog.AddExtension = true;
            }

            if (fileToOpen == null) {
                this.New();
            } else {
                this.Open(fileToOpen);
                if (this.currentDocument == null) {
                    // some error opening the document
                    this.New();
                }
            }
        }

        /*
        class PanelNoScrollOnFocus : Panel {
            protected override System.Drawing.Point ScrollToControl(Control activeControl) {
                return DisplayRectangle.Location;
            }
        }

        public void ScrollInfo(ScrollEventArgs e) {
            var oldCol = Console.ForegroundColor;
            Console.ForegroundColor = ConsoleColor.Magenta;
            Console.WriteLine((e != null ? ("scrolling: " + e.OldValue + " -> " + e.NewValue) : ("        "))
                + "  VertScroll " + ClientRectPanel.VerticalScroll.Value  + "  from " + ClientRectPanel.VerticalScroll.Minimum + "-" + ClientRectPanel.VerticalScroll.Maximum);
            Console.ForegroundColor = oldCol;
        }

        void ClientRectPanel_scroll(object sender, ScrollEventArgs e) {
            ScrollInfo(e);
        }
        */

        /// <summary>
        /// Shall auto-complete be tried?
        /// </summary>
        public bool AutoCompleteEnabled {
            get {
                return this.MenuItem_autocomplete.Checked;
            }
        }


        /// <summary>
        /// Enables/disables all menu items which alter the document.
        /// </summary>
        public void EnableOrDisableDocModifiers(bool state) {
            MenuItem_insertAfter.Enabled = state;
            MenuItem_delete.Enabled = state;
            MenuItem_clearOutput.Enabled = state;
            MenuItem_clearAllOutput.Enabled = state;
            MenuItem_execute.Enabled = state;
            MenuItem_executeFromHere.Enabled = state;
            MenuItem_executeUntilHere.Enabled = state;
            MenuItem_executeWorksheet.Enabled = state;
        }

        /// <summary>
        /// Required, e.g. for <see cref="AutosaveTimer"/>.
        /// </summary>
        internal System.ComponentModel.IContainer components = null;

        /// <summary>
        /// Triggers the autosave.
        /// </summary>
        System.Windows.Forms.Timer AutosaveTimer;

        /// <summary>
        /// na was wohl?
        /// </summary>
        SaveFileDialog saveFileDialog1;

        /// <summary>
        /// na was wohl?
        /// </summary>
        OpenFileDialog openFileDialog1;

        /// <summary>
        /// texlink-file-open
        /// </summary>
        OpenFileDialog selectTexFileDialog;

        string m_DocumentPath = null;

        /// <summary>
        /// the filename of the current document
        /// </summary>
        string DocumentPath {
            get {
                return m_DocumentPath;
            }
            set {
                this.m_DocumentPath = value;
            }
        }

        string FileName {
            get {
                if (DocumentPath == null) {
                    if (this.Altered) {
                        return "Untitled*";
                    } else {
                        return "Untitled";
                    }
                } else {
                    var n = Path.GetFileName(DocumentPath);
                    if (this.Altered) {
                        return n + "*";
                    } else {
                        return n;
                    }
                }
            }
        }

        private void FormClosingEvent(object sender, FormClosingEventArgs e) {
            if (!this.CheckAltered()) {
                e.Cancel = true;
            }
        }

        ///// <summary>
        ///// File selected for TeXlink
        ///// </summary>
        //string TeXlinkFilename;

        /// <summary>
        /// the Document
        /// </summary>
        Document currentDocument;

        /// <summary>
        /// C#-Commands which should be executed.
        /// </summary>
        internal ConcurrentQueue<WorksheetEntry> m_CommandQueue = new ConcurrentQueue<WorksheetEntry>();

        /// <summary>
        /// Task which runs all commands in <see cref="m_CommandQueue"/>.
        /// </summary>
        internal System.Threading.Thread m_ExecutorOfCommandQueue;

        /// <summary>
        /// Setting this to false should terminate <see cref="m_ExecutorOfCommandQueue"/> gracefully;
        /// </summary>
        internal volatile bool m_ExecutorOfCommandQueue_RegularTermination;

        /// <summary>
        /// 
        /// </summary>
        volatile bool m_ExecutorOfCommandQueue_WorkInProgress;

        /// <summary>
        /// Main-method of task/thread <see cref="m_ExecutorOfCommandQueue"/>;
        /// </summary>
        private void ExecutorOfCommandQueue() {
            System.Threading.Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.InvariantCulture;
            m_ExecutorOfCommandQueue_WorkInProgress = false;

            try {

                while (m_ExecutorOfCommandQueue_RegularTermination) {
                    WorksheetEntry nextCommand;
                    if (!m_CommandQueue.IsEmpty && m_CommandQueue.TryDequeue(out nextCommand)) {
                        m_ExecutorOfCommandQueue_WorkInProgress = true;

                        Document.Tuple nextCommand2 = (Document.Tuple)this.Invoke(new Func<Document.Tuple>(delegate () {
                            Document.Tuple ret;
                            if (nextCommand.Deleted) {
                                ret = null;
                            } else {
                                int idx = nextCommand.Index;
                                nextCommand.GlowAnimationTimer.Start();
                                ret = this.currentDocument.CommandAndResult[idx];
                            }
                            return ret;
                        }));

                        if (nextCommand2 != null) {
                            nextCommand2.Evaluate();

                            this.Invoke(new Action(delegate () {
                                nextCommand.EvaluateCommand_Finished();
                            }));
                        }
                    }
                    m_ExecutorOfCommandQueue_WorkInProgress = false;
                    Thread.Sleep(300);
                }
            } catch (ThreadAbortException tae) {
                Console.WriteLine("Worker thread has been aborted!");
                Console.WriteLine(tae.Message);
                Console.WriteLine(tae.StackTrace);
            }
        }

        void UnqueueCommands() {
            while (!m_CommandQueue.IsEmpty) {
                WorksheetEntry dropCommand;
                if (m_CommandQueue.TryDequeue(out dropCommand)) {
                    dropCommand.SetReadOnlyState(false, false);
                }
            }
        }

        /// <summary>
        /// 
        /// </summary>
        void ResetCommandQueue() {
            UnqueueCommands();

            if (m_ExecutorOfCommandQueue != null) {
                // stop thread

                // try regular termination
                m_ExecutorOfCommandQueue_RegularTermination = false;
                Thread.Sleep(800);

                if (m_ExecutorOfCommandQueue.IsAlive) {
                    // hardcore

                    Thread.Sleep(5000);
                    if (m_ExecutorOfCommandQueue.IsAlive) {
                        m_ExecutorOfCommandQueue.Abort();
                    }

                    //m_ExecutorOfCommandQueue.Join();
                }
                m_ExecutorOfCommandQueue = null;

            }

            if (this.CommandsListHead != null)
                this.CommandsListHead.SetReadOnlyState(false, true);

            m_ExecutorOfCommandQueue = new System.Threading.Thread(ExecutorOfCommandQueue);
            m_ExecutorOfCommandQueue.Priority = ThreadPriority.Normal;
            m_ExecutorOfCommandQueue_RegularTermination = true;
            m_ExecutorOfCommandQueue_WorkInProgress = false;
            m_ExecutorOfCommandQueue.Start();
        }



        /// <summary>
        /// triggers a resize of the width of all worksheet controls
        /// </summary>
        void MyclientSizeChanged(object sender, EventArgs e) {
            BlockTextChanged = true;
            this.ClientRectPanel.Width = this.ClientSize.Width - 1;
            DocumentPanel.Width = this.ClientRectPanel.ClientSize.Width - 1;

            this.ClientRectPanel.Height = this.ClientSize.Height - 2 - this.ClientRectPanel.Top;

            if (!BlockResizeClient) {
                CommandsListHead.ResizeAll();
            }
            BlockTextChanged = false;
        }

        private void MenuItem_new_Click(object sender, EventArgs e) {
            New();
        }

        private bool altered;

        internal bool Altered {
            get {
                return altered;
            }
            set {
                altered = value;
                this.Text = "BoSSSPad (" + this.FileName + ")";
            }
        }

        /// <summary>
        /// returns true if the current document can be discarded.
        /// </summary>
        /// <returns></returns>
        internal bool CheckAltered() {
            if (currentDocument == null)
                return true;

            if (Altered) {
                DialogResult r6 = MessageBox.Show(this,
                    "Do you want to save the changes to " + this.FileName + "?",
                    "BoSSS Worksheet",
                    MessageBoxButtons.YesNoCancel,
                    MessageBoxIcon.Question,
                    MessageBoxDefaultButton.Button1);

                switch (r6) {
                    case System.Windows.Forms.DialogResult.Yes:
                        return this.Save();
                    case System.Windows.Forms.DialogResult.No:
                        return true;
                    case System.Windows.Forms.DialogResult.Cancel:
                        return false;
                    default:
                        throw new NotSupportedException();
                }
            } else {
                return true;
            }
        }

        /// <summary>
        /// creates a new document
        /// </summary>
        private void New() {
            if (CheckAltered()) {
                if (this.CommandsListHead != null) {
                    this.CommandsListHead.Destroy();
                    this.CurrentCommand = null;
                }

                this.currentDocument = new Document();
                this.currentDocument.CommandAndResult.Add(new Document.Tuple {
                    Command = "restart"
                });
                this.currentDocument.CommandAndResult.Add(new Document.Tuple());

                this.CommandsListHead = new WorksheetEntry(this, this.currentDocument);
                this.CommandsListHead.EvaluateAllCommands();
                this.CurrentCommand = this.CommandsListHead;
                this.DocumentPath = null;

                this.altered = false;

                ResetCommandQueue();
            }
        }

        private bool Save() {
            if (this.DocumentPath == null) {
                return SaveAs();
            } else if (DocumentPath.ToLowerInvariant().EndsWith(".tex")) {
                var result = MessageBox.Show(
                "This will transfer the Worksheet content to TeX-file '" + Path.GetFileName(this.DocumentPath) + "'." + Environment.NewLine
                + "The content of the tex-file will be altered!" + Environment.NewLine
                + "Ensure the file is saved in the TeX-Editor before clicking 'OK'!" + Environment.NewLine
                + "Afterwards, reload the file.",
                "Transfer from Worksheet to TeX-file...",
                MessageBoxButtons.OKCancel,
                MessageBoxIcon.Exclamation);

                if (result == DialogResult.OK) {
#if !DEBUG
                    try {
#endif
                        LatexIO.UpdateTexFile(this.DocumentPath, this.currentDocument);
#if !DEBUG
                    } catch (Exception exc) {
                        MessageBox.Show(
                            exc.GetType().Name + ":" + Environment.NewLine + exc.Message,
                            "TeXlink Error.",
                            MessageBoxButtons.OK,
                            MessageBoxIcon.Error);
                        return false;
                    }
#endif

                    return true;
                } else {
                    return false;
                }
            } else {

                try {
                    this.currentDocument.Serialize(this.DocumentPath);
                    this.Altered = false;
                } catch (Exception e) {
                    MessageBox.Show(this,
                        e.GetType().Name + ":\n" + e.Message,
                        "Error saving file",
                        MessageBoxButtons.OK,
                        MessageBoxIcon.Error);
                    return false;
                }
            }

            return true;
        }

        private void AutosaveTimer_Tick(object sender, EventArgs ea) {
            //var OldLoc = this.DocumentPanel.Location;
            //OldLoc.Y -= 15;
            //this.DocumentPanel.Location = OldLoc;
            //Console.WriteLine("DocPanel loc {0} {1} ", this.DocumentPanel.Location.X, this.DocumentPanel.Location.Y);
            //return;

            string AutosaveName = this.DocumentPath;
            if (AutosaveName == null) {
                AutosaveName = "Unnamedfile.bws";
            }
            AutosaveName = AutosaveName + "~";

            try {
                if (this.Altered) {
                    this.currentDocument.Serialize(AutosaveName);
                }
            } catch (Exception) {
                // silent Fail?


                //MessageBox.Show(this,
                //    e.GetType().Name + ":\n" + e.Message,
                //    "Error saving file",
                //    MessageBoxButtons.OK,
                //    MessageBoxIcon.Error);
            }
        }

        private bool SaveAs() {

            if (this.DocumentPath != null && this.m_DocumentPath.Length > 0) {
                try {
                    saveFileDialog1.InitialDirectory = Path.GetDirectoryName(this.DocumentPath);
                    saveFileDialog1.FileName = Path.GetFileName(this.DocumentPath);
                } catch (ArgumentException) {

                }
            }
            DialogResult r6 = saveFileDialog1.ShowDialog(this);

            if (r6 == System.Windows.Forms.DialogResult.Cancel)
                return false;

            this.DocumentPath = saveFileDialog1.FileName;
            return this.Save();
        }

        private void OpenDialog() {
            if (CheckAltered()) {
                DialogResult r6 = openFileDialog1.ShowDialog();

                if (r6 == System.Windows.Forms.DialogResult.Cancel) {
                    return;
                }

                Open(openFileDialog1.FileName);
            }
        }

        private void Open(string fileName) {
            if (this.CommandsListHead != null) {
                this.CommandsListHead.Destroy();
                this.CurrentCommand = null;
            }

            ResetCommandQueue();

            try {
                if (fileName.ToLowerInvariant().EndsWith(".tex")) {
                    // try to load as LaTeX
                    // ++++++++++++++++++++
                    List<string> dummy;
                    this.DocumentPath = fileName;
                    LatexIO.SplitTexFile(fileName, out dummy, out this.currentDocument);

                    //this.TeXlinkFilename = fileName;
                    //this.MenuItem_Bws2Tex.Enabled = true;
                } else {
                    // try to load as LaTeX
                    // ++++++++++++++++++++

                    this.currentDocument = Document.Deserialize(fileName);
                    this.Altered = false;
                    this.m_DocumentPath = fileName;
                }
            } catch (Exception e) {
                MessageBox.Show(this,
                    e.GetType().Name + ":\n" + e.Message,
                    "Error opening file",
                    MessageBoxButtons.OK,
                    MessageBoxIcon.Error);
                return;
            }

            if (this.currentDocument.CommandAndResult.Count <= 0) {
                this.currentDocument.CommandAndResult.Add(new Document.Tuple {
                    Command = "\"Empty document - lets try restart?\""
                });
            }

            this.DocumentPath = this.m_DocumentPath;

            this.CommandsListHead = new WorksheetEntry(this, this.currentDocument);
            this.CurrentCommand = this.CommandsListHead;
            //this.CommandsListHead.ResizeAll();

            string workingDirectory = Path.GetDirectoryName(fileName);
            if (workingDirectory != null && workingDirectory.Length > 0 && Directory.Exists(workingDirectory))
                Directory.SetCurrentDirectory(workingDirectory);
        }


        private void MenuItem_open_Click(object sender, EventArgs e) {
            this.OpenDialog();
        }

        private void MenuItem_save_Click(object sender, EventArgs e) {
            this.Save();
        }

        private void MenuItem_saveAs_Click(object sender, EventArgs e) {
            this.SaveAs();
        }

        private void MenuItem_insertAfter_Click(object sender, EventArgs e) {
            this.CurrentCommand.InsertAfter(null);
        }

        private void MenuItem_delete_Click(object sender, EventArgs e) {
            this.CurrentCommand.Delete();
        }

        private void MenuItem_clearOutput_Click(object sender, EventArgs e) {
            this.CurrentCommand.ClearOutput();
        }

        private void MenuItem_clearAllOutput_Click(object sender, EventArgs e) {
            this.CommandsListHead.ClearOutputRecursive();
        }

        private void MenuItem_execute_Click(object sender, EventArgs e) {
            this.CurrentCommand.EvaluateCommand();
        }

        bool CheckReset() {
            if (this.m_ExecutorOfCommandQueue_WorkInProgress || !this.m_CommandQueue.IsEmpty) {
                var res = MessageBox.Show("Commands are currently being executed. Should they be canceled?", "Interrupt execution?", MessageBoxButtons.OKCancel, MessageBoxIcon.Exclamation);

                if (res == DialogResult.Cancel)
                    return false;

                this.ResetCommandQueue();
            }

            return true;
        }


        private void MenuItem_executeFromHere_Click(object sender, EventArgs e) {
            if (!CheckReset())
                return;
            this.CurrentCommand.EvaluateAllCommands();
        }

        private void MenuItem_executeUntilHere_Click(object sender, EventArgs e) {
            if (!CheckReset())
                return;
            this.CommandsListHead.EvaluateCommandsUntil(this.CurrentCommand);
        }

        private void MenuItem_executeWorksheet_Click(object sender, EventArgs e) {
            if (!CheckReset())
                return;
            this.CommandsListHead.EvaluateAllCommands();
        }

        private void MenuItem_interruptCurrent_Click(object sender, EventArgs e) {
            //this.CommandsListHead.EvaluateAllCommands();

            this.ResetCommandQueue();
        }

        private void MenuItem_unqueuePending_Click(object sender, EventArgs e) {
            //this.CommandsListHead.EvaluateAllCommands();

            this.UnqueueCommands();
        }

        private void InitializeComponent() {
            this.SuspendLayout();
            // 
            // Worksheet
            // 
            this.AutoScroll = true;
            this.ClientSize = new System.Drawing.Size(784, 762);
            this.Name = "Worksheet";
            this.ResumeLayout(false);

        }

        /// <summary>
        /// Workaround for wrong word-wrap on start-up of the application
        /// </summary>
        /// <param name="sender"></param>
        /// <param name="e"></param>
        internal static void OnShown(object sender, EventArgs e) {
            Worksheet ws = (Worksheet)sender;
            ws.CommandsListHead.ResizeAll();
        }



    }
}
