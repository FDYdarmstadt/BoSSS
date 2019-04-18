//Js singleton, approch from: 
//https://stackoverflow.com/questions/1479319/simplest-cleanest-way-to-implement-singleton-in-javascript

class Status {
    constructor(){
        this.locked = false;
        this.changed = false;
        this.constructor.prototype.instance = this;
        this.functionsOnLock = [];
        this.functionsOnUnlock = [];
    }

    isLocked(){
        console.log(this.locked);
        return this.locked;
    }

    lock(){
        this.locked = true;
        this.OnLock();
    }
    
    OnLock(){
        console.log("Locking down...");
        for(var i = 0; i < this.functionsOnLock.length; ++i){
            this.functionsOnLock[i]();
        }
    }

    unlock(){
        this.locked = false;
        this.OnUnlock();
    }

    OnUnlock(){
        console.log("Unlocking...");
        for(var i = 0; i < this.functionsOnUnlock.length; ++i){
            this.functionsOnUnlock[i]();
        }
    }

    registerFunctionOnLock(func){
        this.functionsOnLock.push(func);
    }

    registerFunctionOnUnlock(func){
        this.functionsOnUnlock.push(func);
    }


}

var SingletonFactory = (function() {

    var instance; 

    return{

        getInstance: function(){
            if(instance == null){
                instance = new Status();
                instance.constructor = null;
            }
            return instance;
        }
    }   
})();   

const instance = SingletonFactory.getInstance();

export default instance;