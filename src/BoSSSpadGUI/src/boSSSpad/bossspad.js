
import {InteractiveList} from './interactiveList';
import {Editor} from './editor';
import boSSSRuntime from './bosssInterface';
import Split from "split.js";
import {Range as monacoRange} from 'monaco-editor';
import {BoxWithMenu} from './commandBoxes.js'
import Status from './status.js'
require('./css/bossspad.css');

export class BoSSSpad{
    constructor(element){
        this.editor;
        this.locked = false;
        this.userInput;
        this.createHtmlSkeleton(element);
        this.status = Status;
        this.boSSS = new Editor(this.editor, this.status);
        this.userGUI = new InteractiveList(this.userInput, this.status);
        
        //Start loading the different parts of BoSSSpad
        var that = this;
        this.boSSS.load()
        .then(
            function() {
                that.userGUI.load(that.boSSS);
                that.split = Split([that.editor, that.userInput], {
                sizes: [60, 40],
                minSize: 180,
                gutterSize: 4,
                onDrag: that.update.bind(that)
                });
            }
        ).then(
            this.register.bind(this)
        ).then(
            this.update.bind(this)
        );
    }

    status(){
        return this.status;
    } 

    createHtmlSkeleton(element){
        this.userInput = document.createElement("DIV");
        this.userInput.className="userInput";
        this.editor = document.createElement("DIV");
        this.editor.className="editor";
        element.appendChild(this.editor);
        element.appendChild(this.userInput);
    }

    register(){
        this.boSSS.registerContextMenu(
            this.addNewRunCommand.bind(this), 
            'runCommand', 
            'Insert Run Box',
            [
                monaco.KeyMod.CtrlCmd | monaco.KeyCode.KEY_F,
            ]
        );
        this.boSSS.registerContextMenu(
            this.addNewCommentCommand.bind(this), 
            'runCommentCommand', 
            'Insert Comment Box',
            [
                monaco.KeyMod.CtrlCmd | monaco.KeyCode.KEY_K,
            ]
        );
        this.boSSS.registerContextMenu(
            this.removeCommand.bind(this), 
            'removeCommand', 
            'Remove Box',
            [
                monaco.KeyMod.chord(monaco.KeyMod.CtrlCmd | monaco.KeyCode.KEY_D)
            ]
        );
        this.boSSS.onDidChangeModelContent( this.deleteHandler.bind(this));
        this.boSSS.onDidScrollChange(this.userGUI.updateScroll.bind(this.userGUI));
        this.boSSS.registerLanguage_BoSSS(boSSSRuntime.provideAutoComplete.bind(boSSSRuntime));
        window.addEventListener("resize", this.update.bind(this));
        this.status.registerFunctionOnLock(this.boSSS.setReadonly.bind(this.boSSS));
        this.status.registerFunctionOnUnlock(this.boSSS.unsetReadonly.bind(this.boSSS));
    }

    addNewRunCommand(ed){
        if(!this.status.isLocked()){
        var range = ed.getSelection();
        range.startColumn = 1;
        var runCommand = this.userGUI.addNewRunCommand(range);
        }
    }

    addNewCommentCommand(ed){
        if(!this.status.isLocked()){
        var range = ed.getSelection();
        range.startColumn = 1;
        var runCommand = this.userGUI.addNewCommentCommand(range);
        }
    }

    deleteHandler(IModelContentChange){
        this.status.changed = true;
        function constructRangeFromText( text, startLineNumber){
            var formatedText = text.split('\n');
            var lineNumber = formatedText.length - 1;
            var endColumn = formatedText[formatedText.length - 1 ].length;
            return new monacoRange(startLineNumber, 1, startLineNumber + lineNumber, endColumn + 1);
        };
        
        var length = IModelContentChange.changes[0].rangeLength; 
        //If standard delete or if delete with enter-key
        if(  length > 0 ){
            var oldRange = IModelContentChange.changes[0].range;
            var text = IModelContentChange.changes[0].text;
            var newRange = constructRangeFromText(text, oldRange.startLineNumber);
            this.userGUI.deleteCommandSection(oldRange, newRange);
        }
        this.userGUI.update();
    }

    removeCommand(ed){
        if(!this.status.isLocked()){
            var range = ed.getSelection();
            this.userGUI.removeCommand(range);
        }
    }

    update(){
        this.boSSS.update();
        this.userGUI.updateBoxesInternal();
    }

    //Command actions
    //===========================================================================
    async executeFromHere(){
        if(!this.status.isLocked()){
            this.status.lock();
            var lineNumber = this.boSSS.getSelection().selectionStartLineNumber;
            var boxArray = this.userGUI.getAllBoxesFromLine(lineNumber);
            await BoxWithMenu.runBoxes(boxArray);
            this.status.unlock();
        }
    }

    async executeEntireWorksheet(){
        if(!this.status.isLocked()){
            this.status.lock();
            var boxArray = this.userGUI.getAllBoxes();
            await BoxWithMenu.runBoxes(boxArray);
            this.status.unlock();
        }
    }

    async executeUntilHere(){
        if(!this.status.isLocked()){
            this.status.lock();
            var lineNumber = this.boSSS.getSelection().selectionStartLineNumber;
            var boxArray = this.userGUI.getAllBoxesUntilLine(lineNumber);
            await BoxWithMenu.runBoxes(boxArray);
            this.status.unlock();
        }
    }

    async interruptCurrentCommand(){
        this.deQueueAllPendingCommands();
        await boSSSRuntime.forceAbort();
    }

    deQueueAllPendingCommands(){
        BoxWithMenu.deQueue(); 
    }

    //File actions
    //===========================================================================
    resetFile(){
        this.userGUI.reset();
        this.boSSS.reset();
        this.status.changed = false;
        this.status.locked = false;
        BoxWithMenu.reset();
    }

    openFile(path){
        var that = this;
        boSSSRuntime.load(path).then( function(result){
            that.resetFile();
            that.userGUI.setCommandBoxValues(result);
            that.status.changed = false;
        });
    }

    async saveFile(myPath){
        this.status.changed = false;
        var value = this.userGUI.getCommandBoxValues();
        var data = {
            path : myPath,
            commands : value.commands,
            results : value.results
        }
        var result = await boSSSRuntime.save(data);
        return result;
    }

    hasChanged(){
        return this.status.changed;
    }
}

