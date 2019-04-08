import { runInDebugContext } from "vm";

export class BoxWithMenu{
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
        this.buttonsLI = document.createElement("LI");
        this.buttonsLI.className = "runButtonLI";
        UL.appendChild(this.readoutLI);
        UL.appendChild(this.buttonsLI);
        this.buildButtons(this.buttonsLI);
    
        this.registerButtons(parentBox);
        
        this.readOutInterval;
        
        this.decorationClass = 'runDecoration';
        this.overviewRulerColor = 'rgba(196, 223, 231, 0.5)';
        
        BoxWithMenu.zIndex = 2;
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
        var value = this.parentBox.parentList.getSelectionValue(that.parentBox);
        return this.executeCommand(value);
    }

    executeCommand(value){
        var that = this;
        //Start Blinking
        this.readoutLI.classList.toggle("blinkingreadoutLI");
        this.readOutInterval = setInterval( function() {
            that.readoutLI.classList.toggle("blinkingreadoutLI");
        }, 1200);
    
        //Start calculating
        var calculating = new Promise(async function(resolve, reject){
            
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