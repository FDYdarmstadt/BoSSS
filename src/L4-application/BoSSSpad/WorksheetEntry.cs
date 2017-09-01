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

using ilPSP;
using FastColoredTextBoxNS;
using System;
using System.Collections;
using System.Collections.Generic;
using System.Diagnostics;
using System.Drawing;
using System.Linq;
using System.Text.RegularExpressions;
using System.Threading.Tasks;
using System.Windows.Forms;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Represents a single entry in the worksheet which consists of a command
    /// and the corresponding results (box). Handles the rendering of the
    /// affected part of the worksheet, the manipulation of the
    /// <see cref="Document"/>, and events that occur in worksheet controls
    /// </summary>
    class WorksheetEntry {

        Worksheet m_owner;
        WorksheetEntry m_prev = null;
        WorksheetEntry m_next = null;

        /// <summary>
        /// first item of the linked list.
        /// </summary>
        WorksheetEntry Head {
            get {
                if (m_prev == null)
                    return this;
                else
                    return m_prev.Head;
            }
        }

        int m_index = 0;

        private readonly int m_defaultResultCharHeight;

        /// <summary>
        /// index into the document list, <see cref="Document.CommandAndResult"/>
        /// </summary>
        public int Index {
            get {
                return m_index;
            }
            private set {
                if (value != m_index) {
                    m_index = value;
                    if (m_next != null) {
                        m_next.Index = value + 1;
                    }
                }
            }
        }

        Document m_doc;

        /// <summary>
        /// creates a worksheet logic from a <see cref="Document"/>
        /// </summary>
        public WorksheetEntry(Worksheet o, Document doc)
            : this(o, null, doc, 0) {
        }

        /// <summary>
        /// Creates a worksheet logic from a <see cref="Document"/>; recursively;
        /// </summary>
        private WorksheetEntry(Worksheet o, WorksheetEntry p, Document doc, int i)
            : this(o, doc, p, null, doc.CommandAndResult[i].Command, doc.CommandAndResult[i].InterpreterTextOutput) // perform recursion
        {
            Debug.Assert(this.Index == i);
            if (this.Index < doc.CommandAndResult.Count - 1) {
                new WorksheetEntry(o, this, doc, i + 1);
            }
        }

        /// <summary>
        /// Animates color changes during the execution of commands
        /// </summary>
        public System.Windows.Forms.Timer GlowAnimationTimer;

        /// <summary>
        /// private constructor
        /// </summary>
        private WorksheetEntry(Worksheet o, Document doc, WorksheetEntry p, WorksheetEntry n, string commandStr, string resultStr) {
            m_doc = doc;
            m_owner = o;
            if (p != null) {
                this.m_prev = p;
                p.m_next = this;
            }
            if (n != null) {
                this.m_next = n;
                n.m_prev = this;
            }
            if (m_prev != null) {
                this.Index = m_prev.Index + 1;
            }

            m_owner.DocumentPanel.SuspendLayout();
            m_owner.BlockResizeClient = true;

            {
                this.GlowAnimationTimer = new System.Windows.Forms.Timer(this.m_owner.components);
                this.GlowAnimationTimer.Interval = 40; // 25 fps
                this.GlowAnimationTimer.Tick += new System.EventHandler(this.GlowAnimationTimer_Tick);
            }

            {
                Command = new FastColoredTextBox();
                m_owner.DocumentPanel.Controls.Add(Command);

                Command.AutoCompleteBracketsList = new char[] { '(', ')', '{', '}', '[', ']', '\"', '\"', '\'', '\'' };
                Command.AutoIndentExistingLines = false;
                Command.AutoScrollMinSize = new System.Drawing.Size(284, 285);
                Command.BackBrush = null;
                Command.BorderStyle = System.Windows.Forms.BorderStyle.Fixed3D;

                //Command.CharHeight = 15;
                //Command.CharWidth = 7;
                Command.Cursor = System.Windows.Forms.Cursors.IBeam;
                Command.DelayedEventsInterval = 200;
                Command.DelayedTextChangedInterval = 500;
                Command.DisabledColor = System.Drawing.Color.FromArgb(((int)(((byte)(100)))), ((int)(((byte)(180)))), ((int)(((byte)(180)))), ((int)(((byte)(180)))));
                Font fa = m_owner.CommandFont;
                Font fb = new Font("Courier New", 12.0F, FontStyle.Bold, GraphicsUnit.Point, ((byte)(0)));
                Command.Font = fb;

                Command.WordWrapMode = WordWrapMode.CharWrapControlWidth;
                //Command.WordWrapMode = WordWrapMode.CharWrapPreferredWidth;
                Command.WordWrap = true;

                Command.ImeMode = System.Windows.Forms.ImeMode.Off;
                Command.IsReplaceMode = false;
                Command.Location = new System.Drawing.Point(0, 24);
                Command.Name = "fctb";
                Command.Paddings = new System.Windows.Forms.Padding(0);
                Command.PreferredLineWidth = 80;
                Command.ReservedCountOfLineNumberChars = 3;
                Command.SelectionColor = System.Drawing.Color.FromArgb(((int)(((byte)(60)))), ((int)(((byte)(0)))), ((int)(((byte)(0)))), ((int)(((byte)(255)))));
                //Command.ServiceColors = ((FastColoredTextBoxNS.ServiceColors)(resources.GetObject("fctb.ServiceColors")));
                Command.Size = new System.Drawing.Size(346, 311);
                Command.TabIndex = 3;
                Command.Text = commandStr != null ? commandStr : "";
                Command.Zoom = 100;
                Command.SelectionChangedDelayed += new System.EventHandler(Command_SelectionChangedDelayed);
                Command.AutoIndentNeeded += this.Command_AutoIndentNeeded;
                Command.LineRemoved += this.Command_LineRemoved;
                Command.LineInserted += this.Command_LineInserted;
                Command.TextChanged += this.Command_TextChanged;
                Command.ZoomChanged += this.ZoomChanged;
                Command.SelectionChanged += this.SelectionChanged;

                Command.Language = Language.CSharp;
                Command.OnSyntaxHighlight(new TextChangedEventArgs(Command.Range));

                Command.GotFocus += Command_GotFocus;
                Command.KeyDown += this.KeyDown;

                // Map redo action to Ctrl+Y (as in most editors)
                Command.HotkeysMapping.Add(Keys.Control | Keys.Y, FCTBAction.Redo);
            }

            {
                AutoCompleteBox = new AutocompleteMenu(this.Command);
                AutoCompleteBox.MinFragmentLength = 2;

                AutoCompleteBox.Items.MaximumSize = new System.Drawing.Size(200, 300);
                AutoCompleteBox.Items.Width = 200;
                
                AutoCompleteBox.Items.SetAutocompleteItems(
                    new DynamicAutoCompleteList(m_owner));
            }

            {
                m_owner.BlockTextChanged = true;
                Result = new FastColoredTextBox();
                m_owner.DocumentPanel.Controls.Add(Result);

                Result.BorderStyle = System.Windows.Forms.BorderStyle.Fixed3D;
                Result.Cursor = System.Windows.Forms.Cursors.IBeam;
                Result.DisabledColor = System.Drawing.Color.FromArgb(((int)(((byte)(100)))), ((int)(((byte)(180)))), ((int)(((byte)(180)))), ((int)(((byte)(180)))));
                Result.Font = new Font("Courier New", 12.0F, FontStyle.Bold, GraphicsUnit.Point, ((byte)(0)));

                Result.WordWrapMode = WordWrapMode.CharWrapControlWidth;
                Result.WordWrap = true;

                Result.Paddings = new System.Windows.Forms.Padding(0);
                Result.ShowLineNumbers = false;
                Result.SelectionColor = System.Drawing.Color.FromArgb(((int)(((byte)(60)))), ((int)(((byte)(0)))), ((int)(((byte)(0)))), ((int)(((byte)(255)))));

                Result.Zoom = 100;
                Result.ZoomChanged += this.ZoomChanged;

                Result.MouseWheel += Scroll;
                Result.TextChanged += Result_TextChanged;
                Result.GotFocus += Result_GotFocus;
                Result.KeyDown += KeyDown;
                Result.SelectionChanged += this.SelectionChanged;


                Result.BackColor = Color.White;
                Result.BackBrush = null;
                Result.ForeColor = Color.DarkBlue;
                Result.BorderStyle = BorderStyle.None;

                
                m_defaultResultCharHeight = Result.CharHeight;
                Result.Text = resultStr != null ? resultStr : "";
                // This seems to be the only way to tell FastColoredTextBox to
                // be of zero height if empty
                Result.CharHeight = Result.Text.IsNullOrEmpty() ? 0 : m_defaultResultCharHeight;

                Result.ReadOnly = true;
                m_owner.BlockTextChanged = false;
            }

            {
                GraphicsResult = new PictureBox();
                m_owner.DocumentPanel.Controls.Add(GraphicsResult);

                GraphicsResult.SizeMode = PictureBoxSizeMode.Zoom;

                GraphicsResult.Visible = false; // only activated in special cases.
            }

            m_owner.DocumentPanel.ResumeLayout();
            //Command.Show();
            //Result.Show();

            Resize();

            m_owner.BlockResizeClient = false;
        }

        private void ZoomChanged(object sender, EventArgs e) {
            if (!Command_ZoomChanged_locked) {
                FastColoredTextBox box = (FastColoredTextBox)sender;
                this.Head.ChangeZoomRecursive(box.Zoom);
                this.Head.ResizeAll();
            }
        }

        private void SelectionChanged(object sender, EventArgs e) {
            Debug.Assert(object.ReferenceEquals(this.Command, sender) || object.ReferenceEquals(this.Result, sender));
            FastColoredTextBox fcbt = (FastColoredTextBox)sender;
            /*
            //Console.WriteLine("SelectionChanged: from ({1},{2})--({3},{4})", object.ReferenceEquals(this.Command, sender), 
            //    this.Command.Selection.Start.iLine, this.Command.Selection.Start.iChar, this.Command.Selection.End.iLine, this.Command.Selection.End.iChar);


            int y0Cursor, yECursor;
            {
                int CursorLine = fcbt.Selection.Start.iLine;
                var lineInfos = fcbt.LineInfos;
                int LineNoWithWraps = 0;
                for (int iLi = 0; iLi < lineInfos.Count && iLi < CursorLine; iLi++) {
                    var li = lineInfos[iLi];
                    LineNoWithWraps += li.WordWrapStringsCount;
                    //Console.WriteLine("l{0}: wraps:{1} ", iLi, li.WordWrapStringsCount);
                }

                int CursorCol = fcbt.Selection.Start.iChar;
                for (int i = 0; i < lineInfos[CursorLine].WordWrapStringsCount - 1; i++) {
                    if (CursorCol > lineInfos[CursorLine].CutOffPositions[i]) {
                        LineNoWithWraps++;
                    } else {
                        break;
                    }
                }

                y0Cursor = fcbt.CharHeight * LineNoWithWraps + fcbt.Top;
                yECursor = fcbt.CharHeight + y0Cursor;

            }

            int y0Doc = -m_owner.DocumentPanel.Location.Y;
            int yEDoc = m_owner.ClientRectPanel.Height + y0Doc;

            //int WindowHeigh = m_owner.ClientRectPanel.Height;
            //Console.WriteLine("Cursor: {0}--{1}     Window: {2}--{3}", y0Cursor, yECursor, y0Doc, yEDoc);

            if ((yECursor - y0Cursor) < (yEDoc - y0Doc)) {
                // only scroll windows which are large enough 

                if (yECursor > yEDoc) {
                    int delta = yECursor - yEDoc;
                    var Loc = m_owner.DocumentPanel.Location;
                    Loc.Y -= delta;
                    m_owner.DocumentPanel.Location = Loc;
                } else if (y0Cursor < y0Doc) {
                    int delta = y0Cursor - y0Doc;
                    var Loc = m_owner.DocumentPanel.Location;
                    Loc.Y -= delta;
                    Loc.Y = Math.Min(0, Loc.Y);
                    m_owner.DocumentPanel.Location = Loc;
                    m_owner.ClientRectPanel.VerticalScroll.Value = Loc.Y;
                }



                //m_owner.ScrollInfo(null);
            } else {
                //var oldCol = Console.ForegroundColor;
                //Console.ForegroundColor = ConsoleColor.Cyan;
                //Console.WriteLine("scrolling ommited: window zu klein");
                //Console.ForegroundColor = oldCol;
            }



            //if (this.Command.Selection.Start.iLine == this.Command.Selection.End.iLine)
            //    CursorLine = this.Command.Selection.Start.iLine;
            //else if(this.Command.Selection.End.iLine < this.Command.Selection.Start.iLine)
            //    CursorLine = this.Command.Selection.


            /*
           int Y0 = Command.Top;
           int H = Command.Height;

           int Scroll = -m_owner.DocumentPanel.Location.Y;

           int WindowHeigh = m_owner.ClientRectPanel.Height;
           int DocHeight = m_owner.DocumentPanel.Height;

           Console.WriteLine("{0}-{1},  {2}-{3}, {4}", Y0, H, Scroll, WindowHeigh, DocHeight);
           */

        }

        bool Command_ZoomChanged_locked = false;

        void ChangeZoomRecursive(int newZoom) {
            this.Command_ZoomChanged_locked = true;
            this.Command.Zoom = newZoom;
            this.Result.Zoom = newZoom;
            this.Command_ZoomChanged_locked = false;
            if (m_next != null)
                m_next.ChangeZoomRecursive(newZoom);
        }

        private void Command_LineInserted(object sender, LineInsertedEventArgs e) {
            this.ResizeAll();
        }


        private void Command_LineRemoved(object sender, LineRemovedEventArgs e) {
            this.ResizeAll();
        }

        MarkerStyle SameWordsStyle = new MarkerStyle(new SolidBrush(Color.FromArgb(40, Color.Gray)));

        /// <summary>
        /// I have no real idea what this does.
        /// </summary>
        private void Command_SelectionChangedDelayed(object sender, EventArgs e) {
            Command.VisibleRange.ClearStyle(SameWordsStyle);
            if (!Command.Selection.IsEmpty)
                return;//user selected diapason

            //get fragment around caret
            var fragment = Command.Selection.GetFragment(@"\w");
            string text = fragment.Text;
            if (text.Length == 0)
                return;
            //highlight same words
            var ranges = Command.VisibleRange.GetRanges("\\b" + text + "\\b").ToArray();
            if (ranges.Length > 1)
                foreach (var r in ranges)
                    r.SetStyle(SameWordsStyle);
        }

        /// <summary>
        /// I have no real idea what this does.
        /// </summary>
        private void Command_AutoIndentNeeded(object sender, AutoIndentEventArgs args) {
            //block {}
            if (Regex.IsMatch(args.LineText, @"^[^""']*\{.*\}[^""']*$"))
                return;
            //start of block {}
            if (Regex.IsMatch(args.LineText, @"^[^""']*\{")) {
                args.ShiftNextLines = args.TabLength;
                return;
            }
            //end of block {}
            if (Regex.IsMatch(args.LineText, @"}[^""']*$")) {
                args.Shift = -args.TabLength;
                args.ShiftNextLines = -args.TabLength;
                return;
            }
            //label
            if (Regex.IsMatch(args.LineText, @"^\s*\w+\s*:\s*($|//)") &&
                !Regex.IsMatch(args.LineText, @"^\s*default\s*:")) {
                args.Shift = -args.TabLength;
                return;
            }
            //some statements: case, default
            if (Regex.IsMatch(args.LineText, @"^\s*(case|default)\b.*:\s*($|//)")) {
                args.Shift = -args.TabLength / 2;
                return;
            }
            //is unclosed operator in previous line ?
            if (Regex.IsMatch(args.PrevLineText, @"^\s*(if|for|foreach|while|[\}\s]*else)\b[^{]*$"))
                if (!Regex.IsMatch(args.PrevLineText, @"(;\s*$)|(;\s*//)"))//operator is unclosed
                {
                    args.Shift = args.TabLength;
                    return;
                }
        }

        private int CurrentLine {
            get {
                if (this.CurrentBox is TextBox) {
                    return ((TextBox)CurrentBox).GetLineFromCharIndex(((TextBox)CurrentBox).SelectionStart);
                }

                if (this.CurrentBox is FastColoredTextBox) {
                    return ((FastColoredTextBox)CurrentBox).Selection.Start.iLine;
                }


                throw new NotImplementedException();
            }
        }

        private int LastLine {
            get {
                if (this.CurrentBox is TextBox) {
                    return ((TextBox)CurrentBox).GetLineFromCharIndex(Math.Max(((TextBox)CurrentBox).Text.Length - 1, 0));
                }

                if (this.CurrentBox is FastColoredTextBox) {
                    return ((FastColoredTextBox)CurrentBox).LinesCount - 1;
                }


                throw new NotImplementedException();
            }
        }

        /// <summary>
        /// removes the controls from the GUI
        /// </summary>
        public void Destroy() {
            m_owner.DocumentPanel.Controls.Remove(Command);
            m_owner.DocumentPanel.Controls.Remove(Result);
            m_owner.DocumentPanel.Controls.Remove(GraphicsResult);
            if (m_next != null)
                m_next.Destroy();
        }

        /// <summary>
        /// determines the y-Position where the rendering of this command starts
        /// </summary>
        int GetYOffset() {
            if (m_prev == null) {
                return 0;
            } else {
                if (m_prev.GraphicsResult.Visible)
                    return m_prev.GraphicsResult.Bottom;
                else
                    return m_prev.Result.Bottom;
            }
        }

        int LineNoOffset {
            get {
                if (m_prev == null)
                    return 0;
                else
                    return m_prev.Command.LinesCount + m_prev.LineNoOffset;
            }
        }

        /// <summary>
        /// adopts the size of the controls
        /// </summary>
        /// <returns>
        /// true, if anything changed (and subsequent controls have to be updated)
        /// </returns>
        bool Resize() {
            bool heightChanged = false;

            int Width = this.m_owner.ClientRectPanel.ClientSize.Width - 1;
            int Y0 = this.GetYOffset();

            //
            // Command prompt
            //
            Command.Left = 0;
            Command.Width = Width - Command.Left;
            Command.Top = Y0;
            Command.NeedRecalc(true, true);  // Force recalculation of line-count due to word-wrap
            int newCommandHeight = ComputeHeight(Command) + 5;
            if (newCommandHeight != Command.Height) {
                Command.Height = newCommandHeight;
                heightChanged = true;
            }
            this.Command.LineNumberStartValue = (uint)this.LineNoOffset;

            //
            // Result prompt
            //
            Result.Width = Width;
            Result.Top = Command.Bottom;
            Result.Left = 0;
            Result.NeedRecalc(true, true);  // Force recalculation of line-count due to word-wrap
            int newResultHeigh = ComputeHeight(Result);
            if (newResultHeigh != Result.Height) {
                Result.Height = newResultHeigh;
                heightChanged = true;
            }

            if (m_owner.DocumentPanel.Height <= this.Result.Bottom)
                m_owner.DocumentPanel.Height = this.Result.Bottom + 50;

            
            //
            // the optional GraphicsResult
            //
            if (GraphicsResult.Visible) {
                GraphicsResult.Top = Result.Bottom;
                GraphicsResult.Left = 0;

                GraphicsResult.Width = Width;

                int h1 = GraphicsResult.Image.Height;
                int h2 = (int) Math.Round( ((double)(GraphicsResult.Image.Height*GraphicsResult.Width)) / ((double)GraphicsResult.Image.Width));
                GraphicsResult.Height = Math.Min(h1, h2);
                
                if (m_owner.DocumentPanel.Height <= this.GraphicsResult.Bottom)
                    m_owner.DocumentPanel.Height = this.GraphicsResult.Bottom + 50;
            }

            return heightChanged;
        }

        /*
        private int ComputeHeight(TextBoxBase textBox) {
            Size sz = new Size(textBox.ClientSize.Width, int.MaxValue);
            TextFormatFlags flags = TextFormatFlags.WordBreak;
            int padding = 3;
            int borders = textBox.Height - textBox.ClientSize.Height;
            sz = TextRenderer.MeasureText(textBox.Text, textBox.Font, sz, flags);
            int h = sz.Height + borders + padding;
            return h;
        }
        */

        private int ComputeHeight(FastColoredTextBox textBox) {
            var lineInfos = textBox.LineInfos;
            int VisibleNoOfLines = 0;
            foreach (var li in lineInfos) {
                VisibleNoOfLines += li.WordWrapStringsCount;
            }
            int h = textBox.CharHeight * VisibleNoOfLines;
            return h;
        }

        private void Command_GotFocus(object sender, EventArgs e) {
            m_owner.CurrentCommand = this;
        }

        private void Result_GotFocus(object sender, EventArgs e) {
            m_owner.CurrentCommand = this;
        }

        private void KeyDown(object sender, KeyEventArgs e) {
            if (!this.AutoCompleteBox.Visible) {
                switch (e.KeyValue) {
                    case (int)ConsoleKey.Enter:
                        if (!e.Shift) {
                            HandleReturn();
                            e.Handled = true;
                            e.SuppressKeyPress = true;

                        }
                        break;

                    case (int)ConsoleKey.Tab:
                        bool handled = HandleTab();
                        e.Handled = handled;
                        e.SuppressKeyPress = handled;
                        break;

                    case (int)ConsoleKey.UpArrow:
                        e.Handled = HandleUp();
                        break;

                    case (int)ConsoleKey.DownArrow:
                        e.Handled = HandleDown();
                        break;

                    case (int)ConsoleKey.PageUp:
                        e.Handled = HandlePageUp();
                        break;

                    case (int)ConsoleKey.PageDown:
                        e.Handled = HandlePageDown();
                        break;
                }
            }
        }

        /// <summary>
        /// handles the RETURN/ENTER key
        /// </summary>
        private void HandleReturn() {
            // return pressed: start evaluation
            EvaluateCommand();

            // if this is last command: add one
            if (m_next == null) {
                InsertAfter(null);
            }

            // move input focus to next command
            m_next.Command.Focus();
        }

        /// <summary>
        /// Handles the TAB key by showing an auto-completion box
        /// </summary>
        private bool HandleTab() {
            if (this.Command.SelectionLength > 0) {
                return false;
            }

            if (ReadEvalPrintLoop.eval == null) {
                MessageBox.Show("Evaluator not initialized. Try using 'restart' first");
            }

            this.AutoCompleteBox.Show(true);
            return true;
        }

        /// <summary>
        /// Divides numbers and words: "123AND456" -> "123 AND 456"
        /// Or "i=2" -> "i = 2"
        /// </summary>
        class InsertSpaceSnippet : AutocompleteItem {
            string pattern;

            public InsertSpaceSnippet(string pattern)
                : base("") {
                this.pattern = pattern;
            }

            public InsertSpaceSnippet()
                : this(@"^(\d+)([a-zA-Z_]+)(\d*)$") {
            }

            public override CompareResult Compare(string fragmentText) {
                if (Regex.IsMatch(fragmentText, pattern)) {
                    Text = InsertSpaces(fragmentText);
                    if (Text != fragmentText)
                        return CompareResult.Visible;
                }
                return CompareResult.Hidden;
            }

            public string InsertSpaces(string fragment) {
                var m = Regex.Match(fragment, pattern);
                if (m == null)
                    return fragment;
                if (m.Groups[1].Value == "" && m.Groups[3].Value == "")
                    return fragment;
                return (m.Groups[1].Value + " " + m.Groups[2].Value + " " + m.Groups[3].Value).Trim();
            }

            public override string ToolTipTitle {
                get {
                    return Text;
                }
            }
        }

        /// <summary>
        /// Handles the UP key which navigates within the current text box or,
        /// if the first line has been reached, focuses the previous worksheet
        /// element (i.e., <see cref="m_prev"/>).
        /// </summary>
        /// <returns>
        /// True, if the key press has been handled; false otherwise
        /// </returns>
        private bool HandleUp() {
            if (CurrentLine != 0 || PreviousBox == null) {
                return false;
            }

            if (PreviousBox is TextBox) {
                TextBox _PreviousBox = PreviousBox as TextBox;

                _PreviousBox.SelectionStart = _PreviousBox.GetFirstCharIndexFromLine(
                    _PreviousBox.GetLineFromCharIndex(Math.Max(_PreviousBox.Text.Length - 1, 0)));
                _PreviousBox.SelectionLength = 0;
            }
            PreviousBox.Focus();
           
            return true;
        }

        /// <summary>
        /// Handles the DOWN key which navigates within the current text box or,
        /// if the last line has been reached, focuses the previous worksheet
        /// element (i.e., <see cref="m_next"/>).
        /// </summary>
        /// <returns>
        /// True, if the key press has been handled; false otherwise
        /// </returns>
        private bool HandleDown() {
            if ((LastLine > 0 && CurrentLine != LastLine) || NextBox == null) {
                return false;
            }

            if (NextBox is TextBox) {
                ((TextBox)NextBox).SelectionStart = 0;
                ((TextBox)NextBox).SelectionLength = 0;
            }
            if (NextBox is FastColoredTextBox) {
                ((FastColoredTextBox)NextBox).GoHome();
            }

            NextBox.Focus();
            
            return true;
        }

        private Control CurrentBox {
            get {
                if (Command.Focused) {
                    return Command;
                } else {
                    return Result;
                }
            }
        }

        private Control NextBox {
            get {
                // Skip empty result boxes
                if (Command.Focused && !Result.Text.IsNullOrEmpty()) {
                    return Result;
                } else {
                    return (m_next == null ? null : m_next.Command);
                }
            }
        }

        private Control PreviousBox {
            get {
                if (Command.Focused) {
                    if (m_prev == null) {
                        return null;
                    }

                    // Skip empty result boxes
                    if (m_prev.Result.Text.IsNullOrEmpty()) {
                        return m_prev.Command;
                    } else {
                        return m_prev.Result;
                    }
                } else {
                    return Command;
                }
            }
        }

        private bool HandlePageUp() {
            Point oldPos = m_owner.ClientRectPanel.AutoScrollPosition;
            m_owner.ClientRectPanel.AutoScrollPosition = new Point(
                oldPos.X,
                Math.Abs(oldPos.Y) - m_owner.ClientRectPanel.VerticalScroll.LargeChange);
            return true;
        }

        private bool HandlePageDown() {
            Point oldPos = m_owner.ClientRectPanel.AutoScrollPosition;
            m_owner.ClientRectPanel.AutoScrollPosition = new Point(
                oldPos.X,
                Math.Abs(oldPos.Y) + m_owner.ClientRectPanel.VerticalScroll.LargeChange);
            return true;
        }

        /// <summary>
        /// inserts an object into the linked list, after this object;
        /// </summary>
        /// <param name="cmdString"></param>
        /// <returns>the newly created object</returns>
        public WorksheetEntry InsertAfter(string cmdString) {
            //if (m_EvaluationInProgress)
            //    return null;

            WorksheetEntry newEntry;
            if (m_next != null) {
                // insert
                var old_next = this.m_next;
                newEntry = new WorksheetEntry(m_owner, this.m_doc, this, this.m_next, cmdString, null);
                {
                    int currentZoom = this.Command.Zoom;
                    newEntry.Command_ZoomChanged_locked = true;
                    newEntry.Command.Zoom = currentZoom;
                    newEntry.Result.Zoom = currentZoom;
                    newEntry.Command_ZoomChanged_locked = false;
                }
                Debug.Assert(object.ReferenceEquals(m_next, newEntry));
                Debug.Assert(object.ReferenceEquals(newEntry.m_prev, this));
                Debug.Assert(object.ReferenceEquals(old_next.m_prev, newEntry));
                Debug.Assert(object.ReferenceEquals(newEntry.m_next, old_next));
                m_doc.CommandAndResult.Insert(newEntry.Index, new Document.Tuple());
            } else {
                // just add at the end
                newEntry = new WorksheetEntry(m_owner, this.m_doc, this, null, cmdString, null);
                Debug.Assert(object.ReferenceEquals(m_next, newEntry));
                Debug.Assert(object.ReferenceEquals(newEntry.m_prev, this));
                m_doc.CommandAndResult.Add(new Document.Tuple());
                Debug.Assert(newEntry.Index == m_doc.CommandAndResult.Count - 1);
            }

            this.ResizeAll();
            m_owner.Altered = true;
            return newEntry;
        }

        public bool Deleted {
            get;
            private set;
        }

        public WorksheetEntry Delete() {
            //if (m_EvaluationInProgress)
            //    return null;

            if (m_prev == null && m_next == null) {
                // Last statement cannot be deleted
                Console.Beep();
                return this;
            }

            m_doc.CommandAndResult.RemoveAt(Index);
            this.m_owner.Altered = true;

            WorksheetEntry newFocus = null;
            if (m_prev != null) {
                m_prev.m_next = this.m_next;
                newFocus = m_prev;
            }
            if (m_next != null) {
                m_next.m_prev = this.m_prev;
                newFocus = m_next;
                m_next.Index--;
            }
            Deleted = true;
            
            //m_owner.DocumentPanel.Controls.Remove(Prompt);
            m_owner.DocumentPanel.Controls.Remove(Command);
            m_owner.DocumentPanel.Controls.Remove(Result);
            m_owner.DocumentPanel.Controls.Remove(GraphicsResult);

            this.ResizeAll();
            newFocus.Command.Focus();
            m_owner.ClientRectPanel.ScrollControlIntoView(newFocus.Command);

            return newFocus;
        }

        /// <summary>
        /// Deletes the content of the result box
        /// </summary>
        public void ClearOutput() {
            //if (m_EvaluationInProgress)
            //    return;

            if (!this.Result.Text.IsNullOrEmpty()) {
                this.Result.Clear();
                //this.m_owner.Altered = true;

                if (this.GraphicsResult.Visible) {
                    this.GraphicsResult.Visible = false;
                    this.ResizeAll();
                }
            }
        }

        /// <summary>
        /// Recursively deletes the content of the result box for this entry
        /// and all following entries.
        /// </summary>
        public void ClearOutputRecursive() {
            //if (m_EvaluationInProgress)
            //    return;

            if (!this.Result.Text.IsNullOrEmpty()) {
                this.Result.Clear();
                //this.m_owner.Altered = true;
            }

            if (m_next != null) {
                m_next.ClearOutputRecursive();
            }
        }

        /// <summary>
        /// evaluates all commands that follow in the list
        /// </summary>
        public void EvaluateAllCommands() {
            EvaluateCommandsUntil(null);
        }

        /// <summary>
        /// evaluates all commands that follow in the list until
        /// <paramref name="lastEntry"/> is reached.
        /// </summary>
        public void EvaluateCommandsUntil(WorksheetEntry lastEntry) {

            List<WorksheetEntry> EntriesToEvaluate = new List<WorksheetEntry>();
            for (WorksheetEntry next = this; next != null && next != lastEntry; next = next.m_next) {
                //EntriesToEvaluate.Add(next);
                next.EvaluateCommand();
            }

            /*
            var asyncExe = new Task(delegate () {
                foreach (var entry in EntriesToEvaluate) {
                    m_owner.Invoke(new Action(delegate () {
                        entry.Command.Focus();
                        entry.EvaluateCommand();
                    }));
                    while (m_EvaluationInProgress) {
                        System.Threading.Thread.Sleep(100);
                    }
                }
            });
            
            asyncExe.Start();
            */
        }

        ///// <summary>
        ///// While one command is evaluated, this blocks the evaluation of all other commands as well as deletion, etc.
        ///// </summary>
        //volatile static bool m_EvaluationInProgress = false;

        private void GlowAnimationTimer_Tick(object sender, EventArgs ea) {
            var BaseColor = Color.Red;
            //var BaseColor = Color.BurlyWood;
            //var BaseColor = Color.SpringGreen;

            double period = (((double)DateTime.Now.Millisecond) / 1000.0)*Math.PI * 2;
            double modulation = (Math.Sin(period) + 1) * (0.8 / 2.0) + 0.2;

            int r = (int)Math.Round(((double)BaseColor.R) * modulation);
            r = Math.Max(0, Math.Min(byte.MaxValue, r));
            int g = (int)Math.Round(((double)BaseColor.G) * modulation);
            g = Math.Max(0, Math.Min(byte.MaxValue, g));
            int b = (int)Math.Round(((double)BaseColor.B) * modulation);
            b = Math.Max(0, Math.Min(byte.MaxValue, b));

            var myColor = Color.FromArgb(r, g, b);

            this.Command.BackColor = myColor;
        }



        /// <summary>
        /// evaluates this specific command
        /// </summary>
        public void EvaluateCommand() {
            if (this.Command.ReadOnly == true)
                // command is currently being executed resp. queued;
                return;

            int idx = this.Index;
            this.m_doc.CommandAndResult[idx].Command = this.Command.Text;
            this.Result.Text = "";
            if (this.GraphicsResult.Visible)
                this.GraphicsResult.Visible = false;
            this.m_owner.Altered = true;
            this.Command.ReadOnly = true;
            this.Command.BackColor = Color.Salmon;

            m_owner.m_CommandQueue.Enqueue(this); // put the command into the queue; the background thread which waits for commands will execute it.
                        

            /*
            if (m_EvaluationInProgress)
                return;
            m_EvaluationInProgress = true;
            m_owner.EnableOrDisableDocModifiers(false);
            int idx = this.Index;
            this.m_doc.CommandAndResult[idx].Command = this.Command.Text;
            this.Result.Text = "";
            if(this.GraphicsResult.Visible)
                this.GraphicsResult.Visible = false;
            this.m_owner.Altered = true;
            this.Head.SetReadOnlyStateRecursive(true);
            this.GlowAnimationTimer.Start();

            var asyncExe = new Task(delegate () {
                System.Threading.Thread.CurrentThread.CurrentCulture = System.Globalization.CultureInfo.InvariantCulture;
                this.m_doc.CommandAndResult[idx].Evaluate();
                EvaluateCommand_Finished();
            });
            //asyncExe.ContinueWith(EvaluateCommand_Finished);
            asyncExe.Start();

            /*
            this.Result.Text = this.m_doc.CommandAndResult[idx].InterpreterTextOutput;

            // special interpretation
            if (this.m_doc.CommandAndResult[idx].Result is System.Drawing.Image) {
                this.GraphicsResult.Image = (Image)this.m_doc.CommandAndResult[idx].Result;
                this.GraphicsResult.Visible = true;

                this.ResizeAll();
            } else {
                if (this.GraphicsResult.Visible) {
                    this.GraphicsResult.Visible = false;
                    this.ResizeAll();
                }
            }
            */
        }

        
        /// <summary>
        /// Sets the `ReadOnly` property for all <see cref="Command"/>'s.
        /// </summary>
        public void SetReadOnlyState(bool newState, bool Recursive) {
            if (newState == false && this.GlowAnimationTimer.Enabled)
                this.GlowAnimationTimer.Stop();
            this.Command.ReadOnly = newState;
            this.Command.BackColor = newState ?  Color.LightGray : Color.White;
            if (Recursive && this.m_next != null)
                this.m_next.SetReadOnlyState(newState, Recursive);
        }
        




        /// <summary>
        /// Called after the evaluation of the command in the C#-interpreter has been finished.
        /// </summary>
        /// <remarks>
        /// This method is 'triggered' by the background thread, but executed in the user-interface thread.
        /// </remarks>
        internal void EvaluateCommand_Finished() {

            
            int idx = this.Index;
            this.GlowAnimationTimer.Stop();
            this.Result.Text = this.m_doc.CommandAndResult[idx].InterpreterTextOutput;
            this.Command.ReadOnly = false;
            this.Command.BackColor = Color.White;
                        
            // special interpretation
            if (this.m_doc.CommandAndResult[idx].Result is System.Drawing.Image) {
                this.GraphicsResult.Image = (Image)this.m_doc.CommandAndResult[idx].Result;
                this.GraphicsResult.Visible = true;

                this.ResizeAll();
            } else {
                if (this.GraphicsResult.Visible) {
                    this.GraphicsResult.Visible = false;
                    this.ResizeAll();
                }
            }


            /*
            this.m_owner.Invoke(new Action(delegate () {

                int idx = this.Index;
                this.GlowAnimationTimer.Stop();
                this.Result.Text = this.m_doc.CommandAndResult[idx].InterpreterTextOutput;
                this.Head.SetReadOnlyStateRecursive(false);


                // special interpretation
                if (this.m_doc.CommandAndResult[idx].Result is System.Drawing.Image) {
                    this.GraphicsResult.Image = (Image)this.m_doc.CommandAndResult[idx].Result;
                    this.GraphicsResult.Visible = true;

                    this.ResizeAll();
                } else {
                    if (this.GraphicsResult.Visible) {
                        this.GraphicsResult.Visible = false;
                        this.ResizeAll();
                    }
                }

                m_EvaluationInProgress = false;
                m_owner.EnableOrDisableDocModifiers(true);
            }));
            */
        }




        private void Command_TextChanged(object sender, EventArgs e) {
            int idx = this.Index;
            this.m_doc.CommandAndResult[idx].Command = this.Command.Text;
            this.m_doc.CommandAndResult[idx].InterpreterTextOutput = this.Result.Text;
            this.m_owner.Altered = true;
        }

        private void Result_TextChanged(object sender, EventArgs e) {
            if (!Result.Text.IsNullOrEmpty()) {
                Result.CharHeight = m_defaultResultCharHeight;
            }
           
            if (!m_owner.BlockTextChanged) {
                bool hc = this.Resize();
                if (hc)
                    this.ResizeAll();
            }

            m_owner.Altered = true;
        }

        /// <summary>
        /// adopts the size of all controls to the document panel <see cref="Worksheet.DocumentPanel"/>
        /// </summary>
        public void ResizeAll() {
            m_owner.DocumentPanel.SuspendLayout();
            this.ResizeAllRecursive();
            m_owner.DocumentPanel.ResumeLayout(false);
            m_owner.DocumentPanel.PerformLayout();
        }

        /// <summary>
        /// adopts the size of all controls to the document panel <see cref="Worksheet.DocumentPanel"/>
        /// </summary>
        void ResizeAllRecursive() {
            this.Resize();
            if (m_next == null) {
                m_owner.DocumentPanel.Height = this.Result.Bottom + 50;
            } else {
                m_next.ResizeAllRecursive();
            }
        }

        FastColoredTextBox Command;

        AutocompleteMenu AutoCompleteBox;

        FastColoredTextBox Result;

        PictureBox GraphicsResult; 

        private void Scroll(object sender, MouseEventArgs e) {
            Point oldPos = m_owner.ClientRectPanel.AutoScrollPosition;
            m_owner.ClientRectPanel.AutoScrollPosition = new Point(
                oldPos.X,
                Math.Abs(oldPos.Y) - e.Delta);
        }

        /// <summary>
        /// Handles dynamic creation of auto complete items for
        /// <see cref="AutoCompleteBox"/>
        /// </summary>
        internal class DynamicAutoCompleteList : IEnumerable<AutocompleteItem> {

            private Worksheet worksheet;

            public DynamicAutoCompleteList(Worksheet worksheet) {
                this.worksheet = worksheet;
            }

            public IEnumerator<AutocompleteItem> GetEnumerator() {
                return BuildList().GetEnumerator();
            }

            IEnumerator IEnumerable.GetEnumerator() {
                return GetEnumerator();
            }

            private IEnumerable<AutocompleteItem> BuildList() {
                if(!worksheet.AutoCompleteEnabled)
                    yield break;

                if (ReadEvalPrintLoop.eval == null) {
                    yield break;
                }

                FastColoredTextBox command = worksheet.CurrentCommand.Command;
                int pos = command.SelectionStart;

                string textToBeCompleted = command.Text.Substring(0, pos);
                string rest = "";
                if (command.Text.Length > pos) {
                    rest = command.Text.Substring(pos, command.Text.Length - pos);
                }

                int timeout = 1000;
                //timeout *= 3600 * 24;

                string[] completions;
                string originalPrefix;
                bool completed = ReadEvalPrintLoop.eval.TryGetCompletions(
                    textToBeCompleted, out completions, out originalPrefix, timeout);

                if (!completed) {
                    MessageBox.Show(
                        "The requested auto-complete operations timed out, most probably"
                        + " due to a bug in Mono.Csharp which provides the completions."
                        + " Unfortunately, this is not recoverable and you will have to"
                        + " restart the application. Don't forget to save your work!");
                    if (worksheet.CheckAltered()) {
                        worksheet.Close();
                    }
                    yield break;
                }

                if (completions != null) {
                    foreach (var item in completions) {
                        yield return new MethodAutocompleteItem(originalPrefix + item);
                    }
                }

                yield return new InsertSpaceSnippet();
                yield return new InsertSpaceSnippet(@"^(\w+)([=<>!:]+)(\w+)$");
            }
        }
    }
}
