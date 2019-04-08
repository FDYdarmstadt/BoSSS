import {Selection} from './selection.js';

export class Box{
    constructor(range, BoxType, parentList){
        this.LI = document.createElement("LI");
        this.LI.classList.add("userInputRegion");
        this.createDivisionsForSelectAndBox();
        this.bottomHeight;
        this.id;
        this.range = range;
        this.BoxType = BoxType; 
        this.parentList = parentList;
        this.boxContent = new BoxType(this.box, this);
        this.Selection = new Selection(this.selectionDiv, this.boxContent);
        this.IsSmall = true;
    }

    update(){
        this.boxContent.update();
    }

    createDivisionsForSelectAndBox(){
        this.selectionDiv = document.createElement("DIV");
        this.selectionDiv.classList.add("selection");
        this.box = document.createElement("DIV");
        this.box.classList.add("box");
        this.LI.appendChild(this.selectionDiv);
        this.LI.appendChild(this.box);
    }

    toggleEnlarged(){
        if(this.IsSmall){
            this.boxContent.enlarge();
            this.Selection.addHighlight();
        }else{
            this.boxContent.reduce();
            this.Selection.removeHighlight();
        }
        this.IsSmall = ! this.IsSmall;
    }
  
    getDomNode(){
      return this.LI;
    }

    setHeight(height){
      this.LI.style.height = height +"px"; 
    }

    getHeight(){
        //Convert string to int
        var height = this.LI.style.height;
        height = height.substring(0, height.length - 2);
        var heightInt = parseInt(height);
        return heightInt;
    }
    
    setRange(range){
      this.range = range;
    }
}