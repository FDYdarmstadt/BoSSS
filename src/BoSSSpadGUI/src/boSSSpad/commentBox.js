import {BoxWithMenu} from './commandBoxes.js'
//var pandoc = require('node-pandoc');
var markdown = require('markdown').markdown;
var mathjax = require('mathjax-electron');

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

    update(){
        //Nothing to update manually
    };
  }