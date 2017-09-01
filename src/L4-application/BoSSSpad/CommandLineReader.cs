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
using System.Linq;
using System.Text;
using System.IO;
using System.Threading;
using System.Collections.Generic;
using BoSSS.Platform;
using ilPSP;

namespace BoSSS.Application.BoSSSpad {

    /// <summary>
    /// Defines a command line reader that provides functionality for
    /// auto-completion when the Tab key is pressed.
    /// </summary>
    /// <remarks>
    /// Adapted from Mono.Terminal.LineEditor, see
    /// https://github.com/mono/mono/blob/master/mcs/tools/csharp/getline.cs
    /// </remarks>
    public class CommandLineReader : IDisposable {

        /// <summary>
        /// Invoked when the user requests auto-completion using the tab
        /// character. Each handler should add an entry to
        /// <see cref="AutoCompleteEventArgs.Completions"/> for the completions
        /// it found.
        /// </summary>
        public EventHandler<AutoCompleteEventArgs> AutoCompleteEvent;

        /// <summary>
        /// Defines arguments for <see cref="AutoCompleteEvent"/>.
        /// </summary>
        public class AutoCompleteEventArgs : EventArgs {

            /// <summary>
            /// The text entered so far.
            /// </summary>
            public readonly string Text;

            /// <summary>
            /// A list of completions. Any handler for
            /// <see cref="AutoCompleteEvent"/> should add his results to this
            /// array.
            /// </summary>
            public readonly List<CompletionList> Completions = new List<CompletionList>();

            /// <summary>
            /// The text to be completed
            /// </summary>
            /// <param name="text"></param>
            public AutoCompleteEventArgs(string text) {
                Text = text;
            }
        }

        /// <summary>
        /// The text being edited
        /// </summary>
        private StringBuilder text;

        /// <summary>
        /// The current cursor position, i.e. index into <see cref="text"/>.
        /// For an index into <see cref="renderedText"/>, use
        /// <see cref="GetRenderedTextPosition"/>
        /// </summary>
        private int textPosition;

        /// <summary>
        /// The text as it is rendered (replaces (char)1 with ^A on display for
        /// example).
        /// </summary>
        private StringBuilder renderedText;

        /// <summary>
        /// The prompt shown to the user at the beginning of each input line
        /// </summary>
        private string prompt;

        /// <summary>
        /// The row where we started displaying data
        /// </summary>
        private int homeRow;

        /// <summary>
        /// The maximum length that has been displayed on the screen
        /// </summary>
        private int maxRenderedLength;

        /// <summary>
        /// If we are done editing, this breaks the interactive loop
        /// </summary>
        private bool done = false;

        /// <summary>
        /// The object that tracks the history of issued commands
        /// </summary>
        private History history;

        /// <summary>
        /// Handlers for special keys.
        /// </summary>
        private readonly KeyHandler[] keyHandlers;

        /// <summary>
        /// Creates a new reader with the given <paramref name="name"/>.
        /// </summary>
        /// <param name="name">
        /// The name of the application for which this reader is used.
        /// </param>
        /// <param name="historyStoragePath">
        /// The path to the directory where the history of commands should be
        /// stored.
        /// </param>
        /// <param name="maxHistoryLength">
        /// The maximum number of history entries
        /// </param>
        /// <param name="historyIgnoredCommands">
        /// Commands that are ignored in the history, e.g. exit commands that
        /// don't need to be recorded
        /// </param>
        public CommandLineReader(string name, string historyStoragePath = null, int maxHistoryLength = 300, string[] historyIgnoredCommands = null) {
            keyHandlers = new KeyHandler[] {
				new KeyHandler(ConsoleKey.Home, CmdMoveToStartOfLine),
				new KeyHandler(ConsoleKey.End, CmdMoveToEndOfLine),
				new KeyHandler(ConsoleKey.LeftArrow, CmdMoveLeft),
				new KeyHandler(ConsoleKey.RightArrow, CmdMoveRight),
				new KeyHandler(ConsoleKey.UpArrow, CmdHistoryPrevious),
				new KeyHandler(ConsoleKey.DownArrow, CmdHistoryNext),
				new KeyHandler(ConsoleKey.Enter, CmdDone),
				new KeyHandler(ConsoleKey.Backspace, CmdBackspace),
				new KeyHandler(ConsoleKey.Delete, CmdDeleteChar),
				new KeyHandler(ConsoleKey.Tab, CmdAutoComplete),
                new KeyHandler(ConsoleKey.PageUp, CmdHistorySearchBackwards),
                new KeyHandler(ConsoleKey.PageDown, CmdHistorySearchForwards)
			};

            renderedText = new StringBuilder();
            text = new StringBuilder();

            history = new History(name, historyStoragePath, maxHistoryLength, historyIgnoredCommands);
        }

        /// <summary>
        /// Starts reading a command from the command line
        /// </summary>
        /// <param name="prompt">
        /// The prompt to be shown at the beginning of the line. Cannot be
        /// deleted by the user and has no effect on the result.
        /// </param>
        /// <param name="initialCommand">
        /// The command to be displayed initially. Can be deleted by the user.
        /// </param>
        /// <returns>
        /// The command entered by the user after the enter key has been
        /// pressed
        /// </returns>
        public string ReadCommand(string prompt, string initialCommand) {
            done = false;
            history.MoveToNextFreeEntry();
            maxRenderedLength = 0;

            this.prompt = prompt;
            InitText(initialCommand);
            history.AppendEntry(initialCommand);

            do {
                try {
                    EditLoop();
                } catch (ThreadAbortException) {
                    Thread.ResetAbort();
                    Console.WriteLine();
                    SetPrompt(prompt);
                    SetText("");
                }
            } while (!done);
            Console.WriteLine();

            if (text == null) {
                history.Dispose();
                return null;
            }

            string result = text.ToString();
            if (result != "") {
                history.UpdateLastEntry(result);
            } else {
                history.RemoveLastEntry();
            }

            history.Save();

            return result;
        }

        /// <summary>
        /// Updates <see cref="homeRow"/>
        /// </summary>
        /// <param name="screenpos"></param>
        private void UpdateHomeRow(int screenpos) {
            int lines = 1 + (screenpos / Console.WindowWidth);
            homeRow = Console.CursorTop - (lines - 1);
            if (homeRow < 0) {
                homeRow = 0;
            }
        }

        /// <summary>
        /// Writes <see cref="prompt"/> and <see cref="renderedText"/> while
        /// performing wrapping if required.
        /// </summary>
        private void Render() {
            Console.Write(prompt);
            Console.Write(renderedText);

            int max = Math.Max(renderedText.Length + prompt.Length, maxRenderedLength);

            for (int i = renderedText.Length + prompt.Length; i < maxRenderedLength; i++) {
                Console.Write(' ');
            }
            maxRenderedLength = prompt.Length + renderedText.Length;

            // Write one more to ensure that we always wrap around properly if we are at the
            // end of a line.
            Console.Write(' ');

            UpdateHomeRow(max + 1);
        }

        /// <summary>
        /// Similar to <see cref="Render"/>
        /// </summary>
        /// <param name="pos"></param>
        private void RenderFrom(int pos) {
            int rpos = GetRenderedTextPosition(pos);
            int i;

            for (i = rpos; i < renderedText.Length; i++) {
                Console.Write(renderedText[i]);
            }

            if ((prompt.Length + renderedText.Length) > maxRenderedLength) {
                maxRenderedLength = prompt.Length + renderedText.Length;
            } else {
                int max_extra = maxRenderedLength - prompt.Length;
                for (; i < max_extra; i++) {
                    Console.Write(' ');
                }
            }
        }

        /// <summary>
        /// Updates <see cref="renderedText"/> given the current value of
        /// <see cref="text"/>.
        /// </summary>
        private void UpdateRenderedText() {
            renderedText.Length = 0;

            for (int i = 0; i < text.Length; i++) {
                int c = (int)text[i];
                if (c < 26) {
                    if (c == '\t') {
                        renderedText.Append("    ");
                    } else {
                        renderedText.Append('^');
                        renderedText.Append((char)(c + (int)'A' - 1));
                    }
                } else {
                    renderedText.Append((char)c);
                }
            }
        }

        /// <summary>
        /// Computes the position in <see cref="renderedText"/> for a given
        /// position in <see cref="textPosition"/>.
        /// </summary>
        /// <param name="textPosition"></param>
        /// <returns></returns>
        private int GetRenderedTextPosition(int textPosition) {
            int p = 0;

            for (int i = 0; i < textPosition; i++) {
                char c = text[i];

                if (c < 26) {
                    if (c == (int)ConsoleKey.Tab) {
                        p += 4;
                    } else {
                        p += 2;
                    }
                } else {
                    p++;
                }
            }

            return p;
        }

        /// <summary>
        /// Calculates the screen position for a given
        /// <see cref="textPosition"/>
        /// </summary>
        /// <param name="textPosition">
        /// A position within <see cref="text"/>.
        /// </param>
        /// <returns>
        /// The position on the screen corresponding to
        /// <see cref="textPosition"/> within <see cref="text"/>.
        /// </returns>
        private int TextToScreenPos(int textPosition) {
            return prompt.Length + GetRenderedTextPosition(textPosition);
        }

        /// <summary>
        /// The number of lines on the screen.
        /// </summary>
        private int LineCount {
            get {
                return (prompt.Length + renderedText.Length) / Console.WindowWidth;
            }
        }

        /// <summary>
        /// Forcibly moves the cursor to the specified position on the console.
        /// </summary>
        /// <param name="newpos"></param>
        private void MoveCursor(int newpos) {
            textPosition = newpos;

            int actual_pos = prompt.Length + GetRenderedTextPosition(textPosition);
            int row = homeRow + (actual_pos / Console.WindowWidth);
            int col = actual_pos % Console.WindowWidth;

            if (row >= Console.BufferHeight) {
                row = Console.BufferHeight - 1;
            }

            Console.SetCursorPosition(col, row);
        }

        /// <summary>
        /// Moves the cursor position, but only if it has changed.
        /// </summary>
        /// <param name="newpos"></param>
        private void MoveCursorIfChanged(int newpos) {
            if (textPosition == newpos)
                return;

            MoveCursor(newpos);
        }

        /// <summary>
        /// Prints a new character to the console.
        /// </summary>
        /// <param name="c"></param>
        private void InsertChar(char c) {
            int prev_lines = LineCount;
            text = text.Insert(textPosition, c);
            UpdateRenderedText();
            if (prev_lines != LineCount) {
                Console.SetCursorPosition(0, homeRow);
                Render();
                MoveCursor(++textPosition);
            } else {
                RenderFrom(textPosition);
                MoveCursor(++textPosition);
                UpdateHomeRow(TextToScreenPos(textPosition));
            }
        }

        /// <summary>
        /// Moves the cursor to <paramref name="textPosition"/> and renders
        /// <see cref="text"/> starting from <paramref name="textPosition"/>.
        /// </summary>
        /// <param name="textPosition"></param>
        private void MoveToAndRenderAfter(int textPosition) {
            MoveCursor(textPosition);
            RenderFrom(textPosition);
            MoveCursor(textPosition);
        }

        #region Commands

        /// <summary>
        /// Finishes the current line
        /// </summary>
        private void CmdDone() {
            done = true;
        }

        /// <summary>
        /// Handles and auto-complete request according to the result of the
        /// <see cref="AutoCompleteEvent"/> (also, see
        /// <see cref="CompletionList"/>)
        /// </summary>
        private void CmdAutoComplete() {
            if (AutoCompleteEvent == null) {
                InsertChar('\t');
                return;
            }

            AutoCompleteEventArgs eventArgs = new AutoCompleteEventArgs(
                text.ToString().Substring(0, textPosition));
            AutoCompleteEvent(this, eventArgs);

            string[] completions = eventArgs.Completions.
                SelectMany(c => c.Completions).
                Distinct().
                OrderBy(s => s).
                ToArray();

            int noOfCompletions = completions.Length;
            if (noOfCompletions == 0) {
                return;
            }

            if (completions.Length == 1) {
                InsertTextAtCursor(completions[0]);
            } else {
                string commonPrefix = eventArgs.Completions.
                    Select(c => c.CommonPrefix).
                    ElementAtMin(s => s.Length);

                int last = -1;

                for (int p = 0; p < completions[0].Length; p++) {
                    char c = completions[0][p];

                    for (int i = 1; i < noOfCompletions; i++) {
                        if (completions[i].Length < p) {
                            goto mismatch;
                        }

                        if (completions[i][p] != c) {
                            goto mismatch;
                        }
                    }
                    last = p;
                }

            mismatch:
                if (last != -1) {
                    InsertTextAtCursor(completions[0].Substring(0, last + 1));
                }
                Console.WriteLine();
                foreach (string s in completions) {
                    Console.Write(commonPrefix);
                    Console.Write(s);
                    Console.WriteLine();
                }
                Console.WriteLine();
                Render();
                MoveCursor(textPosition);
            }
        }

        /// <summary>
        /// Move the cursor to the beginning of the line
        /// </summary>
        private void CmdMoveToStartOfLine() {
            MoveCursorIfChanged(0);
        }

        /// <summary>
        /// Moves the cursor the end of the line.
        /// </summary>
        private void CmdMoveToEndOfLine() {
            MoveCursorIfChanged(text.Length);
        }

        /// <summary>
        /// Moves the cursor one position to the left
        /// </summary>
        private void CmdMoveLeft() {
            if (textPosition == 0) {
                return;
            }

            MoveCursorIfChanged(textPosition - 1);
        }

        /// <summary>
        /// Moves the cursor one position to the right.
        /// </summary>
        private void CmdMoveRight() {
            if (textPosition == text.Length)
                return;

            MoveCursorIfChanged(textPosition + 1);
        }

        /// <summary>
        /// Delete the last character.
        /// </summary>
        private void CmdBackspace() {
            if (textPosition == 0)
                return;

            text.Remove(--textPosition, 1);
            UpdateRenderedText();
            MoveToAndRenderAfter(textPosition);
        }

        /// <summary>
        /// Deletes the character at the current position.
        /// </summary>
        private void CmdDeleteChar() {
            if (text.Length == 0 || textPosition == text.Length) {
                return;
            }

            text.Remove(textPosition, 1);
            UpdateRenderedText();
            MoveToAndRenderAfter(textPosition);
        }

        /// <summary>
        /// Navigates backwards in the command history.
        /// </summary>
        private void CmdHistoryPrevious() {
            if (!history.PreviousEntryExists)
                return;

            HistoryUpdateLine();

            SetText(history.MoveToPrevious());
        }

        /// <summary>
        /// Navigates forwards in the command history.
        /// </summary>
        private void CmdHistoryNext() {
            if (!history.NextEntryExists)
                return;

            history.UpdateCurrentEntry(text.ToString());
            SetText(history.MoveToNext());

        }

        /// <summary>
        /// Searches backwards through command history with current. 
        /// Known from Linux: history-search-backward 
        /// </summary>
        private void CmdHistorySearchBackwards() {
            int oldTextPosition = textPosition;
            string seachString = text.ToString().Substring(0, textPosition);
            string completeEntry = history.SearchEntryBackwards(seachString);
            if (completeEntry == null)
                return;

            SetText(completeEntry);
            MoveCursor(oldTextPosition);
        }

        /// <summary>
        /// Searches forwards through command history with current text.
        /// Known from Linux: history-search-forward
        /// </summary>
        private void CmdHistorySearchForwards() {
            int oldTextPosition = textPosition;
            string seachString = text.ToString().Substring(0, textPosition);
            string completeEntry = history.SearchEntryForwards(seachString);
            if (completeEntry == null)
                return;

            SetText(completeEntry);
            MoveCursor(oldTextPosition);
        }

        #endregion

        /// <summary>
        /// Adds the current line to the history if needed
        /// </summary>
        private void HistoryUpdateLine() {
            history.UpdateCurrentEntry(text.ToString());
        }

        /// <summary>
        /// Inserts the given <paramref name="newText"/> at the current
        /// position of the cursor
        /// </summary>
        /// <param name="newText">
        /// The text to be inserted.
        /// </param>
        private void InsertTextAtCursor(string newText) {
            int prev_lines = LineCount;
            text.Insert(textPosition, newText);
            UpdateRenderedText();
            if (prev_lines != LineCount) {
                Console.SetCursorPosition(0, homeRow);
                Render();
                textPosition += newText.Length;
                MoveCursor(textPosition);
            } else {
                RenderFrom(textPosition);
                textPosition += newText.Length;
                MoveCursor(textPosition);
                UpdateHomeRow(TextToScreenPos(textPosition));
            }
        }

        /// <summary>
        /// Upon each key press: Checks if there is a handler in
        /// <see cref="keyHandlers"/> that is responsible for this particular
        /// key and calls it if needed.
        /// </summary>
        private void EditLoop() {
            ConsoleKeyInfo cki;

            while (!done) {
                cki = Console.ReadKey(true);

                bool handled = false;
                foreach (KeyHandler handler in keyHandlers) {
                    ConsoleKeyInfo t = handler.KeyInfo;

                    if (t.Key == cki.Key) {
                        handled = true;
                        handler.Action();
                        break;
                    }
                }

                if (handled) {
                    continue;
                }

                // No special key, just insert
                if (cki.KeyChar != (char)0) {
                    InsertChar(cki.KeyChar);
                }
            }
        }

        /// <summary>
        /// Initializes <see cref="text"/> and <see cref="renderedText"/> with
        /// <paramref name="initial"/>.
        /// </summary>
        /// <param name="initial">
        /// The initial text.
        /// </param>
        private void InitText(string initial) {
            text = new StringBuilder(initial);
            UpdateRenderedText();
            textPosition = text.Length;
            Render();
            MoveCursor(textPosition);
        }

        /// <summary>
        /// Overwrites the content of the current line by
        /// <paramref name="newText"/>.
        /// </summary>
        /// <param name="newText"></param>
        private void SetText(string newText) {
            Console.SetCursorPosition(0, homeRow);
            InitText(newText);
        }

        /// <summary>
        /// Writes the prompt <paramref name="newPrompt"/> to the console.
        /// </summary>
        /// <param name="newPrompt"></param>
        private void SetPrompt(string newPrompt) {
            prompt = newPrompt;
            Console.SetCursorPosition(0, homeRow);
            Render();
            MoveCursor(textPosition);
        }

        /// <summary>
        /// Saves the history
        /// </summary>
        public void SaveHistory() {
            if (history != null) {
                history.Save();
            }
        }

        #region IDisposable Members

        /// <summary>
        /// <see cref="SaveHistory"/>
        /// </summary>
        public void Dispose() {
            SaveHistory();
        }

        #endregion

        /// <summary>
        /// Encapsulates a set of auto-completion options for a given prefix.
        /// </summary>
        public struct CompletionList {

            /// <summary>
            /// The common prefix of all completions, i.e. the part of the
            /// completion that can already be applied since it is valid for
            /// all available completions.
            /// </summary>
            public readonly string CommonPrefix;

            /// <summary>
            /// A set of potential completions for <see cref="CommonPrefix"/>.
            /// If this is given by
            /// <list type="bullet">
            ///     <item>
            ///     an array with a single entry, the string should be the text
            ///     to be <b>appended</b> to the existing text. For example,
            ///     if the entered text is "T", the result for a completion of
            ///     "ToString" should be "oString", <b>not</b> "ToString"
            ///     </item>
            ///     <item>
            ///     an array with multiple entries, the result should be the
            ///     <b>full</b> text
            ///     </item>
            /// </list>
            /// </summary>
            public readonly string[] Completions;

            /// <summary>
            /// Constructs a new set of auto-completion options.
            /// </summary>
            /// <param name="commonPrefix">
            /// <see cref="CommonPrefix"/>
            /// </param>
            /// <param name="completions">
            /// See <see cref="Completions"/>
            /// </param>
            public CompletionList(string commonPrefix, string[] completions) {
                if (completions == null) {
                    throw new ArgumentException("Completions must not be null", "completions");
                }

                CommonPrefix = commonPrefix;
                Completions = completions;
            }
        }

        /// <summary>
        /// Combines a <see cref="ConsoleKey"/> with an <see cref="Action"/> in
        /// order to define a handler that can deal with special keys.
        /// </summary>
        private struct KeyHandler {

            /// <summary>
            /// Information about the key that has been pressed.
            /// </summary>
            public readonly ConsoleKeyInfo KeyInfo;

            /// <summary>
            /// The action that should be performed upon this key press.
            /// </summary>
            public readonly Action Action;

            /// <summary>
            /// Creates a new handler
            /// </summary>
            /// <param name="key"></param>
            /// <param name="action"></param>
            public KeyHandler(ConsoleKey key, Action action) {
                KeyInfo = new ConsoleKeyInfo((char)0, key, false, false, false);
                Action = action;
            }
        }

        /// <summary>
        /// Stores the history of commands, even across multiple sessions by
        /// storing it to a file. Emulates the bash-like behavior, where edits
        /// done to the history are recorded
        /// </summary>
        private class History : IDisposable {

            /// <summary>
            /// The entries of the rotating history. Since the history can
            /// rotate, the history does not necessarily start and [0] but
            /// rather at [<see cref="firstEntry"/>].
            /// </summary>
            private string[] history;

            /// <summary>
            /// Points to the next free entry in <see cref="history"/>.
            /// </summary>
            private int nextFreeEntry;

            /// <summary>
            /// Points to the first valid entry in <see cref="history"/>.
            /// </summary>
            private int firstEntry;

            /// <summary>
            /// The current position while navigating through the history.
            /// </summary>
            private int currentEntry;

            /// <summary>
            /// The current number of entries in the history.
            /// </summary>
            private int length;

            /// <summary>
            /// The full path to the history file.
            /// </summary>
            private string historyFilePath;

            /// <summary>
            /// Custom commands that are not included in the history
            /// </summary>
            private string[] ignoredCommands;

            /// <summary>
            /// Stores the latest used searchKey for historySearchBackward/Forward  
            /// </summary>
            private string lastSearchKey;

            /// <summary>
            /// Creates a new history object
            /// </summary>
            /// <param name="appName">
            /// The name of the app for which the history is stored.
            /// </param>
            /// <param name="storageDirectory">
            /// If storage directory is not null, the history will be read from
            /// and written to
            /// <paramref name="storageDirectory"/>/<paramref name="appName"/>.history
            /// </param>
            /// <param name="maxHistorySize">
            /// The maximum number of history entries.
            /// </param>
            /// <param name="historyIgnoredCommands">
            /// Commands that are ignored in the history, e.g. exit commands that
            /// don't need to be recorded
            /// </param>
            public History(string appName, string storageDirectory, int maxHistorySize, string[] historyIgnoredCommands) {
                if (maxHistorySize < 1) {
                    throw new ArgumentException("Invalid history size", "maxHistorySize");
                }

                this.ignoredCommands = historyIgnoredCommands ?? new string[0];
                lastSearchKey = "";

                history = new string[maxHistorySize];

                if (storageDirectory != null) {
                    if (!Directory.Exists(storageDirectory)) {
                        try {
                            Directory.CreateDirectory(storageDirectory);
                        } catch {
                            throw new IOException(String.Format(
                                "Directory '{0}' does not exist and could not be created",
                                storageDirectory));
                        }
                    }

                    historyFilePath = Path.Combine(storageDirectory, appName) + ".history";
                }

                if (File.Exists(historyFilePath)) {
                    using (StreamReader sr = File.OpenText(historyFilePath)) {
                        string line;

                        while ((line = sr.ReadLine()) != null) {
                            if (line != "")
                                AppendEntry(line);
                        }
                    }
                }
            }

            /// <summary>
            /// Appends a value to the history
            /// </summary>
            /// <param name="s">
            /// The value to be appended
            /// </param>
            public void AppendEntry(string s) {
                if (ignoredCommands.Contains(s)) {
                    return;
                }

                history[nextFreeEntry] = s;
                nextFreeEntry = (nextFreeEntry + 1) % history.Length;

                if (nextFreeEntry == firstEntry) {
                    firstEntry = (firstEntry + 1 % history.Length);
                }

                if (length != history.Length) {
                    length++;
                }
            }

            /// <summary>
            /// Updates the current cursor location with the string, to support
            /// editing of history items.
            /// </summary>
            /// <param name="s"></param>
            /// <remarks>
            /// For the current line to participate, an <see cref="AppendEntry"/>
            /// must be called before
            /// </remarks>
            public void UpdateCurrentEntry(string s) {
                history[currentEntry] = s;
            }

            /// <summary>
            /// Removes the last valid entry from the history.
            /// </summary>
            public void RemoveLastEntry() {
                nextFreeEntry = nextFreeEntry - 1;
                if (nextFreeEntry < 0) {
                    nextFreeEntry = history.Length - 1;
                }
            }

            /// <summary>
            /// Updates the last valid entry of the history.
            /// </summary>
            /// <param name="s"></param>
            public void UpdateLastEntry(string s) {
                int t = nextFreeEntry - 1;
                if (t < 0)
                    t = history.Length - 1;

                history[t] = s;
            }

            /// <summary>
            /// Returns true if an entry previous to <see cref="currentEntry"/>
            /// exists.
            /// </summary>
            public bool PreviousEntryExists {
                get {
                    if (length == 0) {
                        return false;
                    }

                    int next = currentEntry - 1;
                    if (next < 0) {
                        next = length - 1;
                    }

                    if (next == nextFreeEntry) {
                        return false;
                    }

                    return true;
                }
            }

            /// <summary>
            /// Returns true if an entry subsequent to
            /// <see cref="currentEntry"/> exists.
            /// </summary>
            public bool NextEntryExists {
                get {
                    if (length == 0) {
                        return false;
                    }

                    int next = (currentEntry + 1) % history.Length;
                    if (next == nextFreeEntry) {
                        return false;
                    }

                    return true;
                }
            }

            /// <summary>
            /// Move the current history position to the previous history entry
            /// and returns the corresponding data.
            /// </summary>
            /// <returns>
            /// If there is sufficient data in the history: A string with the
            /// contents of the previous history entry. Otherwise, null is
            /// returned
            /// </returns>
            public string MoveToPrevious() {
                if (!PreviousEntryExists)
                    return null;

                currentEntry--;
                if (currentEntry < 0)
                    currentEntry = history.Length - 1;

                return history[currentEntry];
            }

            /// <summary>
            /// Move the current history position to the next history entry
            /// and returns the corresponding data.
            /// </summary>
            /// <returns>
            /// If there is sufficient data in the history: A string with the
            /// contents of the next history entry. Otherwise, null is
            /// returned
            /// </returns>
            public string MoveToNext() {
                if (!NextEntryExists)
                    return null;

                currentEntry = (currentEntry + 1) % history.Length;
                return history[currentEntry];
            }

            /// <summary>
            /// Sets the current history position to the next free entry.
            /// </summary>
            public void MoveToNextFreeEntry() {
                if (nextFreeEntry == firstEntry)
                    return;

                currentEntry = nextFreeEntry;
            }

            /// <summary>
            /// Searches for the last entry in the history which starts
            /// with the given searchKey. Search starts from currentEntry.
            /// </summary>
            /// <param name="searchKey"></param>
            /// <returns>
            /// If a command is found, which starts with searchKey: A string  
            /// with the complete command. Otherwise, null is returned.
            /// </returns>
            public string SearchEntryBackwards(string searchKey) {
                if (!lastSearchKey.Equals(searchKey))
                    currentEntry = nextFreeEntry - 1;
                string validEntry = null;
                int i = 0;
                while (currentEntry != firstEntry) {
                    currentEntry--;
                    i++;
                    if (currentEntry < 0)
                        currentEntry = history.Length - 1;
                    if (history[currentEntry].Length > searchKey.Length) {
                        if (history[currentEntry].Substring(0, searchKey.Length).Equals(searchKey)) {
                            validEntry = history[currentEntry];
                            lastSearchKey = searchKey;
                            break;
                        }
                    } 
                }

                return validEntry;
            }

            /// <summary>
            /// Searches for the next entry in the history which starts
            /// with the given searchKey. Search starts from currentEntry.
            /// </summary>
            /// <param name="searchKey"></param>
            /// <returns>
            /// If a command is found, which starts with searchKey: A string  
            /// with the complete command. Otherwise, null is returned.
            /// </returns>
            public string SearchEntryForwards(string searchKey) {
                if (!lastSearchKey.Equals(searchKey))
                    currentEntry = nextFreeEntry - 1;
                string validEntry = null;
                while (currentEntry != nextFreeEntry - 1) {
                    currentEntry++;
                    if (currentEntry > history.Length - 1)
                        currentEntry = 0;
                    if (history[currentEntry].Length > searchKey.Length) {
                        if (history[currentEntry].Substring(0, searchKey.Length).Equals(searchKey)) {
                            validEntry = history[currentEntry];
                            lastSearchKey = searchKey;
                            break;
                        }
                    }
                }

                return validEntry;
            }

            /// <summary>
            /// Saves the history.
            /// </summary>
            public void Save() {
                if (historyFilePath == null) {
                    return;
                }

                using (StreamWriter sw = File.CreateText(historyFilePath)) {
                    int start = (length == history.Length) ? nextFreeEntry : firstEntry;
                    for (int i = start; i < start + length; i++) {
                        int p = i % history.Length;
                        sw.WriteLine(history[p]);
                    }
                }
            }

            #region IDisposable Members

            /// <summary>
            /// Saves the history.
            /// </summary>
            public void Dispose() {
                Save();
            }

            #endregion
        }
    }
}
