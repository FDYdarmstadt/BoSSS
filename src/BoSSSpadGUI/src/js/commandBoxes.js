import { runInDebugContext } from "vm";
import Status from './status.js'

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
        this.IsSelectedToRun = true;

        if(BoxWithMenu.zIndex == null){
            BoxWithMenu.zIndex = 2;
        }
        if(BoxWithMenu.BoxesWaitingForWork == null){
            BoxWithMenu.BoxesWaitingForWork = [];
        }
        if(BoxWithMenu.abort == null){
            BoxWithMenu.abort = false;
        }
    }

    static reset(){
        BoxWithMenu.abort = false;
        BoxWithMenu.zIndex = 2;
        BoxWithMenu.BoxesWaitingForWork = [];
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
  
    async createReadoutContent( value){
        throw "Not Implemented";
        return;
    }

    setValue(editorValue, ReadOutValue){
        throw "Not Implemented";
        return;
    }

    update(){
        throw "Not Implemented";
        return;
    };
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

    registerButtons(){
        var that = this;
        var status = this.parentBox.parentList.status; 
        
        //Run button
        this.runButton.onclick = async function(){
            if (!status.isLocked()){
                status.lock();
                await that.run();
                status.unlock();
            }      
        };
      
        //Run until here button
        this.runUntilHereButton.onclick = async function(){
            if (!status.isLocked()){
                status.lock();
                await that.runUntilHere();
                status.unlock();
            }
        };
        
        //Delete button
        this.deleteButton.onclick = async function(){
            if (!status.isLocked()){
                status.lock();
                await that.delete();
                status.unlock();
            }
        };
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
            that.readoutLI.classList.toggle("blinkingreadoutLI");;
        }, 1200);
    
        //Start calculating
        var calculating = new Promise(async function(resolve, reject){
            //Testdelay
            try{
                await that.createReadoutContent( value);
            }
            catch(e){
                console.log(e);
            }
            //Stop Blinking
            clearInterval(that.readOutInterval, value);
            setTimeout(() => {that.readoutLI.classList.remove("blinkingreadoutLI");}, 150);
            resolve();
        });
        return calculating;
    }

    static AddToBoxesWaitingForWorkList(box){
        box.readoutLI.classList.toggle("blinkingreadoutLI");
        var workBoxes = BoxWithMenu.BoxesWaitingForWork;
        workBoxes.push(box);
    }

    static RemoveFromBoxesWaitingForWorkList(box){
        var workBoxes = BoxWithMenu.BoxesWaitingForWork;
        var index = workBoxes.indexOf(box);
        if(index > -1){
            box.readoutLI.classList.toggle("blinkingreadoutLI");
            workBoxes.splice(index,1);
        }
    }

    static async runBoxes(boxArray){
        
        for(var i = 1; i < boxArray.length; ++i){
            var box = boxArray[i].boxContent;
            if(box.IsSelectedToRun){
                BoxWithMenu.AddToBoxesWaitingForWorkList(box);
            }
        }
        //Start calculation in respective order
        for(var i = 0; i < boxArray.length; ++i){
            if(BoxWithMenu.abort == true){
                BoxWithMenu.abort = false;
                break; 
            }
            var box = boxArray[i].boxContent;

            if(box.IsSelectedToRun){
                await box.run();
                BoxWithMenu.RemoveFromBoxesWaitingForWorkList(box);
            }
        }
    }
  
    async runUntilHere(){
        var boxSubarray = this.parentBox.parentList.getAllBoxesUntil(this.parentBox, this.constructor.prototype.constructor);
        //Change Backgroundcolor
        BoxWithMenu.runBoxes(boxSubarray);
    }

    static deQueue(){
        var status = Status;
        if(status.isLocked() == true){
            BoxWithMenu.abort = true;
            BoxWithMenu.removeBlinkingReadout();
        }
    }

    static removeBlinkingReadout(){
        var workBoxes = BoxWithMenu.BoxesWaitingForWork;
        for(var i = workBoxes.length - 1; i > -1; --i ){
            BoxWithMenu.RemoveFromBoxesWaitingForWorkList(workBoxes[i]);
        }
    }
  
    delete(){
        this.parentBox.parentList.deleteBoxByBox(this.parentBox);
    }
}