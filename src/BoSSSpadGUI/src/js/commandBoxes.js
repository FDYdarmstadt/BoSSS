
//var pandoc = require('node-pandoc');
var markdown = require('markdown').markdown;
var mathjax = require('mathjax-electron');
import boSSSRuntime from './bosssInterface';


class BoxWithMenu{
    constructor(element, parentBox){
      this.parentBox = parentBox;
      this.div = document.createElement("DIV");
      this.div.className = "runBox";
      element.appendChild(this.div);
      
      var UL = document.createElement("UL");
      UL.className ="runUL";
      this.div.appendChild(UL);
  
      this.readoutLI = document.createElement("LI");
      this.readoutLI.className = "readoutLI";
      var buttonsLI = document.createElement("LI");
      buttonsLI.className = "runButtonLI";
      UL.appendChild(this.readoutLI);
      UL.appendChild(buttonsLI);
      this.buildButtons(buttonsLI);
  
      this.registerButtons(parentBox);
      this.registerSizeChange();
      this.readOutInterval;
      
      this.decorationClass = 'runDecoration';
      this.overviewRulerColor = 'rgba(196, 223, 231, 0.5)';
      this.readoutLIIsSmall = true;
    }

    registerSizeChange(){
        var that = this;
        this.readoutLI.addEventListener("click", () => that.changeSize());
    }  

    changeSize(){
        if(this.readoutLIIsSmall){
            this.div.classList.add("enlargedRunBox");
            this.div.style.zIndex = 10;
            this.div.style.backgroundColor = "#f2e6ff";
            this.readoutLI.style.overflow = "hidden";
            var height = this.readoutLI.scrollHeight + 29;
            var parentHeight = this.parentBox.getHeight();
            if(height > parentHeight){
                this.div.style.height = height  + "px";
            }
        }else{
            this.div.classList.remove("enlargedRunBox");
            this.div.style.zIndex = 1;
            this.readoutLI.style.overflow = "auto";
            this.div.style.height = "";
            this.div.style.backgroundColor = "";
        }
        this.readoutLIIsSmall = !this.readoutLIIsSmall;        
    }
  
    //------------------------------abstract funcs ------------------------------------
    runButtonSymbol(){
      throw "Not Implemented";
      return buttonNode;
    }
  
    runUntilHereButtonSymbol(){
      throw "Not Implemented";
      return buttonNode;
    }
  
    async createReadoutContent( readoutNode, value){
      throw "Not Implemented";
      return;
    }

    setValue(editorValue, ReadOutValue){
      throw "Not Implemented";
      return;
    }
    //----------------------------------------------------------------------------------
  
    buildButtons(root){
      this.runButton = document.createElement("SPAN");
      this.runButton.className = "button";
      this.runButton.classList.toggle("bigButton");
      var runSymbol = this.runButtonSymbol();
      this.runButton.appendChild(runSymbol);
      root.appendChild(this.runButton);
  
      this.runUntilHereButton = document.createElement("SPAN");
      this.runUntilHereButton.className = "button";
      this.runUntilHereButton.classList.toggle("bigButton");
      var runAllSymbol = this.runUntilHereButtonSymbol();
      this.runUntilHereButton.appendChild(runAllSymbol);
      root.appendChild(this.runUntilHereButton);
  
      this.deleteButton = document.createElement("Span");
      this.deleteButton.className = "button";
      this.deleteButton.style ="float:right";
      var deleteSymbol = document.createTextNode("x");
      this.deleteButton.appendChild(deleteSymbol);
      root.appendChild(this.deleteButton);
    }
  
    run(){
      var that = this;
      //Start Blinking
      this.readoutLI.classList.toggle("blinkingreadoutLI");
      this.readOutInterval = setInterval( function() {
        that.readoutLI.classList.toggle("blinkingreadoutLI");
      }, 1200);
  
      //Start calculating
      var calculating = new Promise(async function(resolve, reject){
        var value = that.parentBox.parentList.getSelectionValue(that.parentBox);
        
        //Testdelay
        await that.createReadoutContent(that.readoutLI, value);
        
        //Stop Blinking
        clearInterval(that.readOutInterval, value);
        setTimeout(() => {that.readoutLI.classList.remove("blinkingreadoutLI");}, 150);
        resolve();
      });
      return calculating;
    }

    toggleWaitingForWork(){
        this.readoutLI.classList.toggle("blinkingreadoutLI");
    }
  
    async runUntilHere(){
      var boxSubarray = this.parentBox.parentList.getAllBoxesUntil(this.parentBox, this.constructor.prototype.constructor);
      //Change Backgroundcolor
      for(var i = 1; i < boxSubarray.length; ++i){
        await boxSubarray[i].boxContent.readoutLI.classList.toggle("blinkingreadoutLI");
      }
      //Start calculation in respective order
      for(var i = 0; i < boxSubarray.length; ++i){
        await boxSubarray[i].boxContent.run();
      }
    }
  
    delete(){
      this.parentBox.parentList.deleteBoxByBox(this.parentBox);
    }
  
    registerButtons(){
      var that = this;
      var status = this.parentBox.parentList.status; 
      
      //Run button
      this.runButton.onclick = async function(){
        if (!status.isLocked()){
          status.toggleLock();
          await that.run();
          status.toggleLock();
        }      
      };
      
      //Run until here button
      this.runUntilHereButton.onclick = async function(){
        if (!status.isLocked()){
          status.toggleLock();
          await that.runUntilHere();
          status.toggleLock();
        }
      };
      
      //Delete button
      this.deleteButton.onclick = async function(){
        if (!status.isLocked()){
          status.toggleLock();
          await that.delete();
          status.toggleLock();
        }
      };
    }
}

export class RunBox extends BoxWithMenu {
    constructor(element, parentBox){
        super(element, parentBox);
        var hasReadoutContent = false;
        this.img = null;
    }

    runButtonSymbol(){
        var runSymbol = document.createTextNode("Run");
        return runSymbol;
    }

    runUntilHereButtonSymbol(){
        var runSymbol = document.createTextNode("...Run");
        return runSymbol;
    }

    setValue( ReadOutValue){
      //ReadoutContent
      var text = document.createTextNode(ReadOutValue);
      if(this.hasReadoutContent){
        this.readoutLI.replaceChild(text, readoutNode.firstChild);
      }
      else{
        this.hasReadoutContent = true;
        this.readoutLI.appendChild(text);
      }
    }

    async createReadoutContent( readoutNode, value){
        var result = await boSSSRuntime.runCommand(value);
        
        //Write readout into HTML Element
        var text = document.createTextNode(result.Item1);
        if(this.hasReadoutContent){
          readoutNode.replaceChild(text, readoutNode.firstChild);
        }
        else{
          this.hasReadoutContent = true;
          readoutNode.appendChild(text);
        }
        if(this.img != null){
          readoutNode.removeChild(this.img);
          this.img = null;
        }
        if (result.Item2 != null){
           
          this.img = new Image();
          this.img.src = 'data:image/gif;base64,'+ result.Item2;
          readoutNode.appendChild(this.img);
        }
    }
}

export class CommentBox extends BoxWithMenu{
  constructor(element, parentBox){
    super(element, parentBox);
    this.decorationClass = 'commentDecoration';
    this.overviewRulerColor = 'rgba(196, 231, 196, 0.5)';
  }

  runButtonSymbol(){
    var runSymbol = document.createTextNode("Render");
    return runSymbol;
  }

  runUntilHereButtonSymbol(){
    var runSymbol = document.createTextNode("...Render");
    return runSymbol;
  }

  async createReadoutContent( readoutNode, value){
    /*
    var callback = function(err, result) {
      readoutNode.innerHTML = result;
      mathjax.typesetMath(readoutNode.innerHTML);
    }
    var args = '-f latex -t html';
    pandoc(value, args, callback);
    */
    readoutNode.innerHTML = markdown.toHTML(value);
    mathjax.typesetMath(readoutNode);
  }
}