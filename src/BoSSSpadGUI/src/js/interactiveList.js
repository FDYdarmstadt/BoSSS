import {RunBox, CommentBox} from './commandBoxes.js'
import * as monaco from 'monaco-editor';

export class InteractiveList{
    constructor(element, status){
      this.element = element;
      this.UL;
      this.boxes = [];
      this.editor;
      this.status = status;
    }
    
    load(editor){
      var that = this;
      var isLoading = new Promise(function(resolve, reject){
        that.UL = document.createElement("UL");
        that.UL.className = "interactiveListUL";
        that.UL.style.position = "absolute";
        that.element.appendChild(that.UL);        
        that.editor = editor;
        resolve();
      }); 
      return isLoading;
    }
  
    update(){
      if(this.boxes.length > 0){
        this.updateRange();   
        this.updateBoxes();
      }
    }

    updateScroll(){
      if(this.boxes.length > 0){
        this.updateBoxes();
      }
    }
  
    updateBoxes(){
      if(this.boxes.length > 0){
        //Reread position and height from monaco decorations
        var heightInLines = this.boxes[0].range.endLineNumber;
        this.boxes[0].setHeight(heightInLines * 19 );
        
        for(var i = 1; i < this.boxes.length; ++i){
          heightInLines = this.boxes[i].range.endLineNumber - this.boxes[i - 1].range.endLineNumber;
          this.boxes[i].setHeight(heightInLines * 19 );
        }
    
        //Update positions of boxes
        this.UL.style.top = -this.editor.getOffset() + "px";
        //console.log(this.boxes);
      }
    }
  
    updateRange(){
      for(var i = 0; i < this.boxes.length; ++i){
        var range = this.editor.getDecorationRange(this.boxes[i].id);
        this.boxes[i].range = range;
      }
    }
  
    addNewRunCommand(range){
      var newBox = this.insertBox( range, RunBox);
      this.updateBoxes();
      return newBox;
    }
  
    addNewCommentCommand(range){
      var newBox = this.insertBox( range, CommentBox);
      this.updateBoxes();
      return newBox;
    }
  
    getSelectionValue(box){
      //update range
      var range = this.editor.getDecorationRange(box.id);
      box.range = range;
      //get value
      var value = this.editor.getValueInRange(box.range);
      return value;
    }
   
    getAllBoxes(){
        return this.boxes;
    }

    getAllBoxesUntil(box, boxType){
      var check = function(someBox){
        return someBox.BoxType === boxType && someBox.range.startLineNumber <= box.range.startLineNumber;
      }
      return this.boxes.filter(check);
    }
  
    deleteBoxByBox(box){
      var index = this.boxes.indexOf(box);
      this.deleteBox(index);
      this.update();
    }
  
    removeCommand(range){
      this.updateRange();
      this.removeBox(range);
      this.update();
    }

    deleteCommandSection(oldRange, newRange){

      var oldBox = this.findBox(oldRange.startLineNumber);
      //Check if Range is contained in oldBox, if so do nothing
      var rangeContainedInBox = false;
      if(oldBox != null)
      {
        var rangeContainedInBox = monaco.Range.containsRange(oldBox.range, oldRange);
      }
      if(rangeContainedInBox === false)
      {
        this.removeCommand(newRange);
        if(oldBox != null)
        {
          oldBox.range.endLineNumber = newRange.endLineNumber;
        }
      }
      this.updateBoxes();  
    }

    findBox(startLineNumber){
      var range = new monaco.Range(startLineNumber, 1, startLineNumber, 1);
      for(var i = 0; i < this.boxes.length; ++i){
        //End if out of Range
        if(range.endLineNumber < this.boxes[i].startLineNumber)
          return null;
  
        var intersection = this.boxes[i].range.intersectRanges(range);
        if (intersection != null){
          return this.boxes[i];  
        }
      }
    }
  
    removeBox( range){
      //find all intersections, then cut intersections into new Boxes
      for(var i = 0; i < this.boxes.length; ++i){
        //End if out of Range
        if(range.endLineNumber < this.boxes[i].startLineNumber)
          return;
  
        var intersection = this.boxes[i].range.intersectRanges(range);
        if (intersection != null){
          var range1 = null;
          var range2 = null;
          
          //Box in Range
          if(intersection.equalsRange(this.boxes[i].range)){
            this.deleteBox(i);
            --i;
            continue;
          }
          
          var BoxType = this.boxes[i].BoxType;
          //Range in Box
          if(intersection.equalsRange(range)){
            //Cut into 2 new Boxes
            if(intersection.startLineNumber - 1 >= this.boxes[i].range.startLineNumber ){
              var endColumn = this.editor.getLineLastNonWhitespaceColumn(intersection.startLineNumber - 1);
              range1 = new monaco.Range(this.boxes[i].range.startLineNumber,1, intersection.startLineNumber - 1 ,endColumn);
            }
            if(this.boxes[i].range.endLineNumber >= intersection.endLineNumber + 1){
              var endColumn = this.editor.getLineLastNonWhitespaceColumn(this.boxes[i].range.endLineNumber);
              range2 = new monaco.Range(intersection.endLineNumber + 1, 1, this.boxes[i].range.endLineNumber, endColumn);
            }
              
            this.deleteBox(i);
            if(range1 != null )
              this.insertBox(range1, BoxType);
            if(range2 != null)
              this.insertBox(range2, BoxType);
          }
  
          //Range part of Box
          else{
            //Resize Box
            if(intersection.startLineNumber != this.boxes[i].range.startLineNumber){
              var endColumn = this.editor.getLineLastNonWhitespaceColumn(intersection.startLineNumber - 1);
              range1 = new monaco.Range(this.boxes[i].range.startLineNumber,1, intersection.startLineNumber - 1,endColumn);
            } 
            else{
              var endColumn = this.editor.getLineLastNonWhitespaceColumn(this.boxes[i].range.endLineNumber);
              range1 = new monaco.Range(intersection.endLineNumber + 1, 1, this.boxes[i].range.endLineNumber, endColumn);
            }
            this.deleteBox(i);
            this.insertBox(range1, BoxType);
            continue;
          }
        }
      }
    }
  
    insertBox( range, BoxType){
      //change range, so that it holds full endline
      range.endColumn = this.editor.getLineLastNonWhitespaceColumn(range.endLineNumber);
      range.endColumn = range.endColumn + 1;

      var newBox;
      //If first box 
      if(this.boxes.length === 0){
        newBox = new Box( range, BoxType, this);
        this.boxes.push(newBox);
        //Handle Dom Stuff
        this.UL.appendChild(newBox.getDomNode());
      }else{
        //Update Range Information
        this.updateRange();
        //find box that is on bottom of new box, if there is none return null
        //while doing that find all ranges that intersect with this range
        var bottomBoxAdress = null;
        var adressOfBoxesToChange = [];  
        for(var i = 0; i < this.boxes.length; ++i){
          var intersection = this.boxes[i].range.intersectRanges(range);
          if (intersection != null ){
            if(this.boxes[i].BoxType === BoxType){
              range = range.plusRange(this.boxes[i].range);
            }
            if(range.endLineNumber < this.boxes[i].range.startLineNumber){
              bottomBoxAdress = i;
              break;
            }
            adressOfBoxesToChange.push(i);
          }
          if(range.endLineNumber < this.boxes[i].range.startLineNumber){
            bottomBoxAdress = i;
            break;
          }
        }
        //If there aren't any intersections, add new box at respective position
        //Else Extend first box in adressOfBoxesToChange and delete the remaining boxes in adressOfBoxesToChange   
        newBox = new Box( range, BoxType, this);
        this.deleteBoxIndiceRange(adressOfBoxesToChange);
        if (bottomBoxAdress === null){
          this.UL.appendChild(newBox.getDomNode());
          this.boxes.push(newBox);
        }else{
          bottomBoxAdress -= adressOfBoxesToChange.length;
          this.UL.insertBefore(newBox.getDomNode(), this.boxes[bottomBoxAdress].getDomNode()); 
          this.boxes.splice(bottomBoxAdress,0,newBox);
        }
        
      }
  
      //Add Box to Editor and give it its ID.
      var id = this.editor.addDecoration(range, newBox.boxContent.decorationClass, newBox.boxContent.overviewRulerColor);
      newBox.id = id;
      return newBox;
    }

    appendBox(BoxType){
      if(this.boxes.length === 0){
        var range = new monaco.Range (1,0,1,0);
        this.editor.setValue(range, "\n\n");
      }else{
        var range = new monaco.Range( 
          this.boxes[this.boxes.length - 1].range.endLineNumber + 2,
          0,
          this.boxes[this.boxes.length - 1].range.endLineNumber + 2,
          0
        );
        this.editor.setValue(range, "\n\n");
      }
      var newBox = this.insertBox(range, BoxType);
      return newBox;
    }
  
    //indice Array must be sorted
    deleteBoxIndiceRange(indiceArray){
      for(var i = 0; i < indiceArray.length; ++i){
        this.deleteBox(indiceArray[i]-i);
      }
    }
    
    deleteBox(indice){
      //Remove from editor
      this.editor.removeDecoration(this.boxes[indice].id);
  
      //Remove from DOM
      this.UL.removeChild(this.boxes[indice].getDomNode());
  
      //Remove from boxArray
      this.boxes.splice(indice , 1);  
    }

    getCommandBoxValues(){
      var myCommands = [];
      var myResults = [];
      for(var i = 0; i < this.boxes.length; ++i){
        var box = this.boxes[i];
        if(box.BoxType  === RunBox ){
          myCommands.push(this.getSelectionValue(box));
          try{
            myResults.push(box.boxContent.readoutLI.firstChild.innerHTML);
          }catch(err){
            console.log("Did not get result: " + err );
            myResults.push("");
          }
        }
      }
      return{
        commands : myCommands,
        results : myResults
      }
    }

    setCommandBoxValues(data){
      var myCommands = data.Item1;
      var myResults = data.Item2;
      for(var i = 0; i < myCommands.length; ++i){
        var commandBox = this.appendBox(RunBox);
        commandBox.boxContent.setValue(myResults[i]);
        this.editor.setValue(commandBox.range, myCommands[i]);
      }
    }

    reset(){
      for(var i = this.boxes.length - 1; i >= 0; --i){
        this.deleteBox(i);
      }
    }
}
  
class Box{
    constructor(range, BoxType, parentList){
      this.LI = document.createElement("LI");
      this.LI.style="position:relative";
      this.bottomHeight;
      this.id;
      this.range = range;
      this.BoxType = BoxType; 
      this.parentList = parentList;
      this.boxContent = new BoxType(this.LI, this);
    }
  
    getDomNode(){
      return this.LI;
    }

    setHeight(height){
      this.LI.style.height = height +"px"; 
    }
    
    setRange(range){
      this.range = range;
    }


}
  
