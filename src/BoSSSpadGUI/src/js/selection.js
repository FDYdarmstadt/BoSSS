export class Selection{
    constructor(element){
        this.element = element;
        var x = document.createElement("INPUT");
        x.setAttribute("type", "checkbox");
        element.appendChild(x);
    }

    addHighlight(){
        this.element.classList.add("enlargedSelection");
    }

    removeHighlight(){
        this.element.classList.remove("enlargedSelection");
    }
}