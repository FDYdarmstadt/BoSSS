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
    }

    registerSizeChange(){
        var that = this;
        this.readoutLI.addEventListener("click", () => that.parentBox.toggleEnlarged());
    }  

    enlarge(){
        this.div.classList.add("enlargedRunBox");
        this.div.style.zIndex = BoxWithMenu.zIndex;
        BoxWithMenu.zIndex = BoxWithMenu.zIndex + 1;

        this.readoutLI.style.overflow = "hidden";

        var height = this.readoutLI.scrollHeight + 29;
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
        var that = this;
        var status = this.parentBox.parentList.status; 
        this.lastErrorButton.onclick = async function(){
            if (!status.isLocked()){
              status.toggleLock();
              await that.executeCommand("LastError");
              status.toggleLock();
            }      
        };
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

    async createReadoutContent( readoutNode, command){
        var result = await boSSSRuntime.runCommand(command);
        
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