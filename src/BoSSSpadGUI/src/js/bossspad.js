
import {InteractiveList} from './interactiveList';
import {Editor} from './editor';
import boSSSRuntime from './bosssInterface';
import Split from "split.js";
require('../css/bossspad.css');

export class BoSSSpad{
  constructor(element){
    this.editor;
    this.locked = false;
    this.userInput;
    this.createHtmlSkeleton(element);
    this.status = new (this.status());
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
    return class {
      constructor(){
        this.locked = false;
      }

      isLocked(){
        return this.locked;
      }
      toggleLock(){
        this.locked = !this.locked;
      }
    }   
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
    this.boSSS.registerContextMenu(this.addNewRunCommand.bind(this), 'runCommand', 'Insert Run Box');
    this.boSSS.registerContextMenu(this.addNewCommentCommand.bind(this), 'runCommentCommand', 'Insert Comment Box');
    this.boSSS.registerContextMenu(this.removeCommand.bind(this), 'removeCommand', 'Remove');
    this.boSSS.onDidScrollChange(this.userGUI.update.bind(this.userGUI));
    window.addEventListener("resize", this.update.bind(this));
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

  removeCommand(ed){
    if(!this.status.isLocked()){
      var range = ed.getSelection();
      this.userGUI.removeCommand(range);
    }
  }

  update(){
    this.boSSS.update();
  }

  resetFile(){
    this.userGUI.reset();
    this.boSSS.reset();
  }

  openFile(path){
    var that = this;
    boSSSRuntime.load(path).then( function(result){
      that.resetFile();
      that.userGUI.setCommandBoxValues(result);
    });
    
  }

  saveFile(myPath){
    var value = this.userGUI.getCommandBoxValues();
    var data = {
      path : myPath,
      commands : value.commands,
      results : value.results
    }
    boSSSRuntime.save(data);
  }

}

