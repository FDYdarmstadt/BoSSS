import {BoxWithMenu} from './commandBoxes.js'
import boSSSRuntime from './bosssInterface';

export class RunBox extends BoxWithMenu {
    constructor(element, parentBox){
        super(element, parentBox);
        var hasReadoutContent = false;
        this.img = null;
        this.createRunBoxButtons();
        this.registerRunBoxButtons();
        this.registerSizeChange();
        this.ErrorIsOnDisplay = false;
    }

    registerSizeChange(){
        var that = this;
        this.readoutLI.addEventListener("click", () => that.parentBox.toggleEnlarged());
    } 

    update(){
        if(!this.parentBox.IsSmall){
            this.adaptHeightToContent();
        }
    }

    enlarge(){
        this.div.classList.add("enlargedRunBox");
        this.div.style.zIndex = BoxWithMenu.zIndex;
        BoxWithMenu.zIndex = BoxWithMenu.zIndex + 1;

        this.readoutLI.style.overflow = "hidden";
        this.adaptHeightToContent();
    }

    adaptHeightToContent(){
        var innerHeight = this.readoutLI.scrollHeight;
        var height = innerHeight + 27;
        var parentHeight = this.parentBox.getHeight();
        if(height > parentHeight){
            this.div.style.height = height  + "px";
        }

    }

    reduce(){
        this.div.classList.remove("enlargedRunBox");
        this.div.style.zIndex = 1;
        this.readoutLI.style.overflow = "auto";
        this.div.style.height = "";
        this.div.style.backgroundColor = "";
    }

    createRunBoxButtons(){
        this.lastErrorButton = document.createElement("SPAN");
        this.lastErrorButton.className = "button";
        this.lastErrorButton.classList.toggle("bigButton");
        var lastErrorSymbol = this.lastErrorSymbol();
        this.lastErrorButton.appendChild(lastErrorSymbol);
        this.buttonsLI.appendChild(this.lastErrorButton);
    }

    registerRunBoxButtons(){
        var status = this.parentBox.parentList.status; 
        var that = this;
        this.lastErrorButton.onclick = function(){
            if (!status.isLocked()){
                status.lock();
                that.toggleDisplayError(); 
                status.unlock();
            }      
        }
    }

    lastErrorSymbol(){
        var runSymbol = document.createTextNode("Error");
        return runSymbol;
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

    toggleDisplayError(){
        if(!this.ErrorIsOnDisplay){
            this.display(this.error);
        }else{
            this.display(this.result);
        }
        this.ErrorIsOnDisplay = !this.ErrorIsOnDisplay;
    }

    display(boSSSobject){
        if(boSSSobject != null){
            var readoutNode = this.readoutLI;
            var text = document.createTextNode(boSSSobject.Item1);
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
            if (boSSSobject.Item2 != null){
               
              this.img = new Image();
              this.img.src = 'data:image/gif;base64,'+ boSSSobject.Item2;
              readoutNode.appendChild(this.img);
            }
        }
        this.update();
    }

    async createReadoutContent(command){
        this.result = await boSSSRuntime.runCommand(command);
        this.error = await boSSSRuntime.runCommand("LastError");
        this.ErrorIsOnDisplay = false;
        this.display(this.result);
    }

    getResultString(){
        if(this.result != null){
            return this.result.Item1
        }else{
            return "";
        }
    }
}