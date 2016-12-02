var M = 10;	//Used in Streamlines	
	var L = 10; //Used in Streamlines	
var Ht = 0.5;

var Streamline = function(){
	var fwd = [];
	var bwd = [];
	var origin = new Point();
	this.fwd = fwd;
	this.bwd = bwd;
	this.origin = origin;
	for(var i=0;i<M+L-1;i++){
		fwd.push(new Point());
		bwd.push(new Point());
	}
};

/*
inline Point &operator[](int m) {
    if (m == 0) 
      return origin;
    else if (m>0)
      return fwd[m-1];
    else
      return bwd[-m-1];
    
  }

*/
Streamline.prototype.getIth = function(m){
	if(m === 0)
		return this.origin;
	else if(m>0)
		return this.fwd[m-1];
	else
		return this.bwd[-m-1];

};

Streamline.prototype.getFwd = function(){
	return this.fwd;
};

Streamline.prototype.setFwd = function(f,k){
	 this.fwd[k] = f;
};

Streamline.prototype.getOrigin = function(){
	return this.origin;
};

Streamline.prototype.setOrigin = function(o){
	 this.origin = o;
};

Streamline.prototype.getBwd = function(){
	return this.bwd;
};
Streamline.prototype.setBwd = function(b,k){
	 this.bwd[k] = b;
};