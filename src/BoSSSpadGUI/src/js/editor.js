
import * as monaco from 'monaco-editor';


export class Editor{
    constructor(element, status){
      this.element = element;
      this.status = status;
    }
  
    load(){
      var that = this;
      var loading = new Promise(function(resolve, reject){
        var myValue = [
              'restart'
            ].join('\n');
        attachEditor(that.element, myValue, 'csharp')
        .then(
          function(editor){
            that.monaco = editor;
            resolve("Monaco loaded");
          }
        )
      });
      return loading;
    }
  
    registerContextMenu(func, myId, myName){
      //Add Run Template
      this.monaco.addAction({
        // An unique identifier of the contributed action.
        id: myId,
  
        // A label of the action that will be presented to the user.
        label: myName,
  
        // A precondition for this action.
        precondition: null,
  
        // A rule to evaluate on top of the precondition in order to dispatch the keybindings.
        keybindingContext: null,
  
        contextMenuGroupId: 'navigation',
  
        contextMenuOrder: 1.5,
  
        // Method that will be executed when the action is triggered.
        // @param editor The editor instance is passed in as a convinience
        run: func
      });
    }
  
    onDidScrollChange( func){
      this.monaco.onDidScrollChange(func);
    }
  
    addDecoration(myRange, decorationClassName, overviewRulerColor){
  
      var newDecorationId = this.monaco.deltaDecorations([], [
          {
              range: myRange,
              options: {
                  isWholeLine: true,
                  className: decorationClassName,
            overviewRuler: {
              color: overviewRulerColor
            } 
              }
          }
      ]);
      return newDecorationId;
    }
    
    removeDecoration(id){
      this.monaco.deltaDecorations([id], []);
    }
  
    getDecorationRange(id){
      return this.monaco.getModel().getDecorationRange(id);
    }
  
    getOffset(){
      return this.monaco.getScrollTop();
    }
  
    getValueInRange(range){
      return this.monaco.getModel().getValueInRange(range);
    }
  
    update(){
      this.monaco.layout();
    }
}

self.MonacoEnvironment = {
	getWorkerUrl: function (moduleId, label) {
		if (label === 'json') {
			return './json.worker.bundle.js';
		}
		if (label === 'css') {
			return './css.worker.bundle.js';
		}
		if (label === 'html') {
			return './html.worker.bundle.js';
		}
		if (label === 'typescript' || label === 'javascript') {
			return './ts.worker.bundle.js';
		}
		return './editor.worker.bundle.js';
	}
}

function attachEditor( element, myValue, myLanguage, myLineNumbers){
  var promise = new Promise(function(resolve, reject) {
    var editor =  monaco.editor.create(element, {
      value: myValue,
      language: myLanguage,
      lineNumbers: myLineNumbers,
      minimap: {
        enabled: false
      }
    });
    resolve(editor);
  });
  return promise;
}