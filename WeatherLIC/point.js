var Point = function(x, y, z, h) {
  this.x = x || 0;
  this.y = y || 0;
  this.z = z || 0;
  this.h = h || 1;
};

Point.prototype.getX = function(){
	return this.x;
}

Point.prototype.getY = function(){
	return this.y;
}
Point.prototype.getZ = function(){
	return this.z;
}


Point.prototype.setX = function(x){
	this.x = x;
}

Point.prototype.setY = function(y){
	this.y = y;
}

