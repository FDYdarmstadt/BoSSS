export class Selection{
    constructor(element, commandBox){
        this.element = element;
        this.commandBox = commandBox;
        this.addCheckBox();
    }

    addCheckBox(){
        var checkBox = document.createElement("INPUT");
        checkBox.setAttribute("type", "checkbox");
        checkBox.checked = true;
        this.element.appendChild(checkBox);
        this.registerCheckbox(checkBox);
    }

    registerCheckbox(checkBox){
        var that = this;
        checkBox.onclick = function(){
            if(checkBox.checked == true){
                that.commandBox.IsSelectedToRun = true;
            }else{
                that.commandBox.IsSelectedToRun = false;
            }
        }
    }

    addHighlight(){
        this.element.classList.add("enlargedSelection");
    }

    removeHighlight(){
        this.element.classList.remove("enlargedSelection");
    }
}