/*  Global class for simulating the movement of particle through a 1km wind grid

    credit: All the credit for this work goes to: https://github.com/cambecc for creating the repo:
      https://github.com/cambecc/earth. The majority of this code is directly take nfrom there, since its awesome.

    This class takes a canvas element and an array of data (1km GFS from http://www.emc.ncep.noaa.gov/index.php?branch=GFS)
    and then uses a mercator (forward/reverse) projection to correctly map wind vectors in "map space".

    The "start" method takes the bounds of the map at its current extent and starts the whole gridding,
    interpolation and animation process.
*/



var Windy = function( params ){



	var VELOCITY_SCALE = 0.011;             // scale for wind velocity (completely arbitrary--this value looks nice)
	var INTENSITY_SCALE_STEP = 10;            // step size of particle intensity color scale
	var MAX_WIND_INTENSITY = 40;              // wind velocity at which particle intensity is maximum (m/s)
	var MAX_PARTICLE_AGE = 100;                // max number of frames a particle is drawn before regeneration
	var PARTICLE_LINE_WIDTH = 2;              // line width of a drawn particle
	var PARTICLE_MULTIPLIER = 1/30;              // particle count scalar (completely arbitrary--this values looks nice)
	var PARTICLE_REDUCTION = 0.75;            // reduce particle count to this much of normal for mobile devices
	var FRAME_RATE = 20;                      // desired milliseconds per frame
	var BOUNDARY = 0.45;

	var NULL_WIND_VECTOR = [NaN, NaN, null];  // singleton for no wind in the form: [u, v, magnitude]
	var TRANSPARENT_BLACK = [255, 0, 0, 0];

	var MAX_VECTOR_MAG = 0;						//MAx vector magnitude
	var τ = 2 * Math.PI;
	var H = Math.pow(10, -5.2);
	var MIN_NUM_HITS = 12;
	var NUM_COLORS = 4;



  	var ctx = document.createElement('canvas').getContext('2d');
  	//Generate Noise of dim
	function genNoise(width, height){
		var img, imgData, i, to;
		img = ctx.createImageData(width, height);
		imgData = img.data;
		for (i = 0, to = width * height * 4; i < to; i += 4) {
			imgData[i] = imgData[i + 1] = imgData[i + 2] = Math.random() >= 0.5 ? 255 : 0;
			imgData[i + 3] = 255;
		}
		return img;
	}
	
	var point = [];
	point.push(new Point(100,100,100,1));
	point.push(new Point(150,0,0,1));
	point.push(new Point(255,0,0,1));
	point.push(new Point(0,0,255,1));


  // interpolation for vectors like wind (u,v,m)
  var bilinearInterpolateVector = function(x, y, g00, g10, g01, g11) {
      var rx = (1 - x);
      var ry = (1 - y);
      var a = rx * ry,  b = x * ry,  c = rx * y,  d = x * y;
      var u = g00[0] * a + g10[0] * b + g01[0] * c + g11[0] * d;
      var v = g00[1] * a + g10[1] * b + g01[1] * c + g11[1] * d;
      return [u, v, Math.sqrt(u * u + v * v)];
  };


	var createWindBuilder = function(uComp, vComp) {
	  var uData = uComp.data, vData = vComp.data;
	  return {
	      header: uComp.header,
	      //recipe: recipeFor("wind-" + uComp.header.surface1Value),
	      data: function(i) {
	          return [uData[i], vData[i]];
	      },
	      //interpolate: bilinearInterpolateVector
	  }
	};

	function zeros(dimensions) {
	    var array = [];
	    for (var i = 0; i < dimensions[0]; ++i) {
	        array.push(dimensions.length == 1 ? 0 : zeros(dimensions.slice(1)));
	    }
	    return array;
	}

	var createBuilder = function(data) {
		var uComp = null, vComp = null, scalar = null;

		data.forEach(function(record) {

			switch (record.header.parameterCategory + "," + record.header.parameterNumber) {
			  case "2,2": uComp = record; break;
			  case "2,3": vComp = record; break;
			  default:
			    scalar = record;
			}
		});
		//console.log(uComp.data[65159],vComp.data[65159], scalar);
		return createWindBuilder(uComp, vComp);
	};

  	var buildGrid = function(data, callback) {
		var builder = createBuilder(data);

		var header = builder.header;
		var λ0 = header.lo1, φ0 = header.la1;  // the grid's origin (e.g., 0.0E, 90.0N)
		var Δλ = header.dx, Δφ = header.dy;    // distance between grid points (e.g., 2.5 deg lon, 2.5 deg lat)
		var ni = header.nx, nj = header.ny;    // number of grid points W-E and N-S (e.g., 144 x 73)
		var date = new Date(header.refTime);
		date.setHours(date.getHours() + header.forecastTime);

		//display date
		document.getElementById('time').innerHTML = date; 

		var hitData = zeros([nj, ni+1]);
		var texData = zeros([nj, ni+1]);
		var Idata = zeros([nj, ni+1]);

		var imgNoise = genNoise(nj,ni+1);
		// Scan mode 0 assumed. Longitude increases from λ0, and latitude decreases from φ0.
		// http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_table3-4.shtml
		var grid = [];
		var p = 0;
		var isContinuous = Math.floor(ni * Δλ) >= 360;

		var rows = nj, cols = ni;

		for (var j = 0; j < nj; j++) {
			var row = [];
			for (var i = 0; i < ni; i++, p++) {
		  		row[i] = builder.data(p);
		  		var vx = builder.data(p)[0];
		  		var vy = builder.data(p)[1];
		  		var mag = Math.sqrt(vx * vx + vy * vy)
		  		if(mag > MAX_VECTOR_MAG)
		  			MAX_VECTOR_MAG = mag;
		  		texData[j][i] = imgNoise.data[j*nj + i];
			}
			if (isContinuous) {
	  			// For wrapped grids, duplicate first column as last column to simplify interpolation logic
	  			row.push(row[0]);
			}
  			grid[j] = row;
		}


		function gridData(i,j){
			return grid[i][j];
		};

		function getVector(p){
			var i,j,bool;
			var arr = getPixel(p);
			i = arr[0],j=arr[1],bool = arr[2];
			
		    if (bool) {
		    	//console.log(gridData(i,j));
	    		return gridData(i,j);
		    }
		    //console.log("here");
		    return([0,0]);
		};

		function getPixel(p){
			var i_ = parseInt(p.getX() + rows)%rows;
		    var j_ = parseInt(p.getY() + cols)%cols;
		    //console.log(i_,(p.getY() + cols)%cols, p.getY(), cols);
		    if (i_ >= 0 && i_ < rows && j_ >=0 && j_ < cols)
		      return [i_,j_,1];
		    return [i_,j_,0];
		};
      /*
      function interpolate(λ, φ) {
          var i = floorMod(λ - λ0, 360) / Δλ;  // calculate longitude index in wrapped range [0, 360)
          var j = (φ0 - φ) / Δφ;                 // calculate latitude index in direction +90 to -90

          var fi = Math.floor(i), ci = fi + 1;
          var fj = Math.floor(j), cj = fj + 1;

          var row;
          if ((row = grid[fj])) {
              var g00 = row[fi];
              var g10 = row[ci];
              if (isValue(g00) && isValue(g10) && (row = grid[cj])) {
                  var g01 = row[fi];
                  var g11 = row[ci];
                  if (isValue(g01) && isValue(g11)) {
                      // All four points found, so interpolate the value.
                      return builder.interpolate(i - fi, j - fj, g00, g10, g01, g11);
                  }
              }
          }
          return null;
      }
      */
      	function iszero(v){
      		//console.log(v.length);
      		if(v[0]*v[0] + v[1]*v[1] === 0)
      			return 1;
      		else
      			return 0;
      	};

      	function unit(v){
      		var size = Math.sqrt(v[0]*v[0] + v[1]*v[1]);
      		return [v[0]/size, v[1]/size];
      	};

      	function RK(p,h,s){
			var f_ = new Point();
			
			//var k1  = new Vector(),k2  = new Vector(), k3  = new Vector(),k4  = new Vector();
			var v = getVector(p);
			if (!iszero(v))
		    	v = unit(v);

		    var x1 = v[0]*h*0.5, y1 = v[1]*h*0.5;
		    var p1 = new Point(p.getX() + x1, p.getY() + y1);
		    v = getVector(p1);
		    if (!iszero(v))
		    	v = unit(v);

		    var x2 = v[0]*h*0.5, y2 = v[1]*h*0.5;
		    var p2 = new Point(p.getX() + x2, p.getY() + y2);
		    v = getVector(p2);
		    if (!iszero(v))
		    	v = unit(v);

		    var x3 = v[0]*h*0.5, y3 = v[1]*h*0.5;
		    var p3 = new Point(p.getX() + x3, p.getY() + y3);
		    v = getVector(p3);
		    if (!iszero(v))
		    	v = unit(v);


		    var x4 = v[0]*h*0.5, y4 = v[1]*h*0.5;

		    f_.setX(p.getX() + x1/6 + x2/3 + x3/3 + x4/6);
		    f_.setY(p.getY() + y1/6 + y2/3 + y3/3 + y4/6);
		    
		    /*
			p += k1/6 + k2/3 + k3/3 + k4/6;
			*/
			return f_;
		};

		function getT(p){
			/*
			inline double getT(Point &p)
			  {
			    int i,j;

			    if (getPixel(p,i,j)) 
			      return texdata[i][j];
			    return 0;
			  }
			*/
			var arr = getPixel(p);
			if(arr[2] == true)
				return texData[arr[0]][arr[1]];
			return 0;

		};


		function validpt(p){
			
			/*
				inline int validpt(Point &p) {
				    int i,j;

				    if (getPixel(p,i,j))
				      return 1;
				    return 0;
				  }
			*/
			var i,j;
			var arr = getPixel(p);
			if(arr[2] == true)
				return 1;
			return 0;
		};

		function ComputeI(s, nsum){
			/*
						double T,k,I;
					    int i,j;

					    T=0;
					    numvalid = 0;
					    
					    for(i=-L; i<= L; i++) {
					      if (validpt(s[i])) {
							T += getT(s[i]);
							numvalid++;
					      }
					    }
					    if (getPixel(s[0],i,j)) {
					      k = 1./numvalid;
					      Idata[i][j] += I = T*k;
					      hitdata[i][j]++;
					      return I;
					    }
					    return 0;
			*/
			var T = 0,k,I;
			var i,j;
			nsum = 0;

			for(i = -L;i<=L;i++){
				if(validpt(s.getIth(i))){
					T += getT(s.getIth(i));
					nsum++;
				}
			}

			var arr = getPixel(s.getIth(0));
			if(arr[2] == true){
				k = 1/nsum;
				Idata[arr[0]][arr[1]] += I = T*k;
				hitData[arr[0]][arr[1]]++;
				return [I,nsum];
			}

			return [0,nsum];
		};


		function ComputeIFwd(s,I,m,numvalid){
			/*inline double ComputeIFwd(StreamLine &s, double &I, int m, int &numvalid)
			  {
			    int i,j;
			    double k;

			    if (getPixel(s[m],i,j)) {
			      if (validpt(s[m+L]))
				numvalid++;
			      if (validpt(s[m-1-L]))
				numvalid--;
			      k = 1./numvalid;
			      Idata[i][j] += I += k*(getT(s[m+L]) - getT(s[m-1-L]));
			      hitdata[i][j]++;
			      return I;
			    }
			    return 0;
			  }*/
			  var k;
			  var arr = getPixel(s.getIth(m));
			  if(arr[2] == true){
				if(validpt(s.getIth(m + L)))
					numvalid++;
				if(validpt(s.getIth(m - 1 - L)))
					numvalid--;
				k = 1/numvalid;
				Idata[arr[0]][arr[1]] += I += k*(getT(s.getIth(m + L)) - getT(s.getIth(m - 1 - L)));
				hitData[arr[0]][arr[1]]++;
				return [I, numvalid];

			  }
			  return [0,numvalid];
		}

		function ComputeIBwd(s,I,m,numvalid){

			var k;
			  var arr = getPixel(s.getIth(m));
			  if(arr[2] == true){
				if(validpt(s.getIth(m - L)))
					numvalid++;
				if(validpt(s.getIth(m + 1 + L)))
					numvalid--;
				k = 1/numvalid;
				Idata[arr[0]][arr[1]] += I += k*(getT(s.getIth(m - L)) - getT(s.getIth(m + 1 + L)));
				hitData[arr[0]][arr[1]]++;
				return [I, numvalid];

			  }
			  return [0,numvalid];
		}
		function Normalize(){
			 for (var i=0; i<rows; i++) 
      			for (var j=0; j<cols; j++) {
					Idata[i][j] /= hitData[i][j];
     			}

     		//console.log(Idata[1][1],Idata[2][2]);
		}

		function docolors(){
			/*
				void
				docolors() {
				  table[0] = Point(100,100,100);
				  table[1] = Point(150,0,0);
				  table[2] = Point(255,0,0);
				  table[3] = Point(0,0,255);
				}
			*/
		};

		callback( {
			date: date,
			hitData : hitData,
			texData :texData,
			Idata :Idata,
			gridData:gridData,
			getVector : getVector,
			RK:RK,
			getPixel : getPixel,
			ComputeI : ComputeI,
			ComputeIFwd: ComputeIFwd,
			ComputeIBwd:ComputeIBwd,
			Normalize:Normalize,
			rows: rows,
			colms : cols
          //interpolate: interpolate
		});
	};

	

  	var lic = function(grid,callback){
  		
  		var I0;
  		var I;
  		var m;
  		var nsum=0, tmpsum=0;
  		var gridRows = grid.texData.length, gridColms = grid.texData[0].length;

  		for(var i =0;i<gridRows;i++){
  			for(var j =0;j<gridColms;j++){
  				if(grid.hitData[i][j] < MIN_NUM_HITS)
  				{
  					//Generate StreamLine
  					var s = new Streamline();
	  				var b,f;
					b = f = new Point(i+1, j+0.5, 0, 1);
					s.setOrigin(b);
					for(var k =0;k<M+L-1;k++){
						f = grid.RK(f, Ht, s);
						s.setFwd(f,k);
						b = grid.RK(b, -Ht, s);
						s.setBwd(b,k);
					} 

					//Compute I
					/*
						I = I0 = grid->ComputeI(s,nsum);
						tmpsum = nsum;
					*/
					var arr =  grid.ComputeI(s,nsum);
					I = I0 = arr[0];
					nsum = arr[1];
					tmpsum = nsum;

					//Compute Fwd
					for (m=1; m < M; m++) 
					{
						var arr = grid.ComputeIFwd(s,I,m,tmpsum);
						I = arr[0];
						tmpsum = arr[1];
					}
					I = I0;
					tmpsum = nsum;

					//Compute Bwd
					for (m=1; m < M; m++){
	  					var arr = grid.ComputeIBwd(s,I,-m,tmpsum);
	  					I = arr[0];
						tmpsum = arr[1];
					}
				}
  			}
  		}

  		grid.Normalize();

  		callback( {
			grid:grid,
		});


  	};

  /**
   * @returns {Boolean} true if the specified value is not null and not undefined.
   */
  var isValue = function(x) {
      return x !== null && x !== undefined;
  }

  /**
   * @returns {Number} returns remainder of floored division, i.e., floor(a / n). Useful for consistent modulo
   *          of negative numbers. See http://en.wikipedia.org/wiki/Modulo_operation.
   */
  var floorMod = function(a, n) {
      return a - n * Math.floor(a / n);
  }

  /**
   * @returns {Number} the value x clamped to the range [low, high].
   */
  var clamp = function(x, range) {
      return Math.max(range[0], Math.min(x, range[1]));
  }

  /**
   * @returns {Boolean} true if agent is probably a mobile device. Don't really care if this is accurate.
   */
  var isMobile = function() {
      return (/android|blackberry|iemobile|ipad|iphone|ipod|opera mini|webos/i).test(navigator.userAgent);
  }

  /**
   * Calculate distortion of the wind vector caused by the shape of the projection at point (x, y). The wind
   * vector is modified in place and returned by this function.
   */
  var distort = function(projection, λ, φ, x, y, scale, wind, windy) {
      var u = wind[0] * scale;
      var v = wind[1] * scale;
      var d = distortion(projection, λ, φ, x, y, windy);

      // Scale distortion vectors by u and v, then add.
      wind[0] = d[0] * u + d[2] * v;
      wind[1] = d[1] * u + d[3] * v;
      return wind;
  };

  var distortion = function(projection, λ, φ, x, y, windy) {
      var τ = 2 * Math.PI;
      var H = Math.pow(10, -5.2);
      var hλ = λ < 0 ? H : -H;
      var hφ = φ < 0 ? H : -H;

      var pλ = project(φ, λ + hλ,windy);
      var pφ = project(φ + hφ, λ, windy);

      // Meridian scale factor (see Snyder, equation 4-3), where R = 1. This handles issue where length of 1º λ
      // changes depending on φ. Without this, there is a pinching effect at the poles.
      var k = Math.cos(φ / 360 * τ);
      return [
          (pλ[0] - x) / hλ / k,
          (pλ[1] - y) / hλ / k,
          (pφ[0] - x) / hφ,
          (pφ[1] - y) / hφ
      ];
  };



  var createField = function(columns, bounds, callback) {

      /**
       * @returns {Array} wind vector [u, v, magnitude] at the point (x, y), or [NaN, NaN, null] if wind
       *          is undefined at that point.
       */
      function field(x, y) {
          var column = columns[Math.round(x)];
          return column && column[Math.round(y)] || NULL_WIND_VECTOR;
      }

      // Frees the massive "columns" array for GC. Without this, the array is leaked (in Chrome) each time a new
      // field is interpolated because the field closure's context is leaked, for reasons that defy explanation.
      field.release = function() {
          columns = [];
      };

      field.randomize = function(o) {  // UNDONE: this method is terrible
          var x, y;
          var safetyNet = 0;
          do {
              x = Math.round(Math.floor(Math.random() * bounds.width) + bounds.x);
              y = Math.round(Math.floor(Math.random() * bounds.height) + bounds.y)
          } while (field(x, y)[2] === null && safetyNet++ < 30);
          o.x = x;
          o.y = y;
          return o;
      };

      //field.overlay = mask.imageData;
      //return field;
      callback( bounds, field );
  };

  var buildBounds = function( bounds, width, height ) {
      var upperLeft = bounds[0];
      var lowerRight = bounds[1];
      var x = Math.round(upperLeft[0]); //Math.max(Math.floor(upperLeft[0], 0), 0);
      var y = Math.max(Math.floor(upperLeft[1], 0), 0);
      var xMax = Math.min(Math.ceil(lowerRight[0], width), width - 1);
      var yMax = Math.min(Math.ceil(lowerRight[1], height), height - 1);
      return {x: x, y: y, xMax: width, yMax: yMax, width: width, height: height};
  };

  var deg2rad = function( deg ){
    return (deg / 180) * Math.PI;
  };

  var rad2deg = function( ang ){
    return ang / (Math.PI/180.0);
  };

  //Inverse Mapping Equation(Points on flat(x,y) -> lamda, sigma)
  var invert = function(x, y, windy){
    var mapLonDelta = windy.east - windy.west;
    var worldMapRadius = windy.width / rad2deg(mapLonDelta) * 360/(2 * Math.PI);
    var mapOffsetY = ( worldMapRadius / 2 * Math.log( (1 + Math.sin(windy.south) ) / (1 - Math.sin(windy.south))  ));
    var equatorY = windy.height + mapOffsetY;
    var a = (equatorY-y)/worldMapRadius;

    var lat = 180/Math.PI * (2 * Math.atan(Math.exp(a)) - Math.PI/2);
    var lon = rad2deg(windy.west) + x / windy.width * rad2deg(mapLonDelta);
    return [lon, lat];
  };

  var mercY = function( lat ) {
    return Math.log( Math.tan( lat / 2 + Math.PI / 4 ) );
  };


  var project = function( lat, lon, windy) { // both in radians, use deg2rad if neccessary
    var ymin = mercY(windy.south);
    var ymax = mercY(windy.north);
    var xFactor = windy.width / ( windy.east - windy.west );
    var yFactor = windy.height / ( ymax - ymin );

    var y = mercY( deg2rad(lat) );
    var x = (deg2rad(lon) - windy.west) * xFactor;
    var y = (ymax - y) * yFactor; // y points south
    return [x, y];
  };




  var interpolateField = function( grid, bounds, extent, callback ) {

    var projection = {};
    var velocityScale = VELOCITY_SCALE;

    var columns = [];
    var x = bounds.x;

    function interpolateColumn(x) {
        var column = [];
        for (var y = bounds.y; y <= bounds.yMax; y += 2) {
                var coord = invert( x, y, extent );
                if (coord) {
                    var λ = coord[0], φ = coord[1];
                    if (isFinite(λ)) {
                        var wind = grid.interpolate(λ, φ);
                        if (wind) {
                            wind = distort(projection, λ, φ, x, y, velocityScale, wind, extent);
                            column[y+1] = column[y] = wind;

                        }
                    }
                }
        }
        columns[x+1] = columns[x] = column;
    }

    (function batchInterpolate() {
                var start = Date.now();
                while (x < bounds.width) {
                    interpolateColumn(x);
                    x += 2;
                    if ((Date.now() - start) > 1000) { //MAX_TASK_TIME) {
                        setTimeout(batchInterpolate, 25);
                        return;
                    }
                }
          createField(columns, bounds, callback);
    })();
  };

/*
	var animate = function(bounds, field) {

	function asColorStyle(r, g, b, a) {
	    return "rgba(" + 243 + ", " + 243 + ", " + 238 + ", " + a + ")";
	}

	function hexToR(h) {return parseInt((cutHex(h)).substring(0,2),16)}
	function hexToG(h) {return parseInt((cutHex(h)).substring(2,4),16)}
	function hexToB(h) {return parseInt((cutHex(h)).substring(4,6),16)}
	function cutHex(h) {return (h.charAt(0)=="#") ? h.substring(1,7):h}

	function windIntensityColorScale(step, maxWind) {

	    var result = [
	      // blue to red
	      //"rgba(" + hexToR('#178be7') + ", " + hexToG('#178be7') + ", " + hexToB('#178be7') + ", " + 0.5 + ")",
	      //"rgba(" + hexToR('#8888bd') + ", " + hexToG('#8888bd') + ", " + hexToB('#8888bd') + ", " + 0.5 + ")",
	      //"rgba(" + hexToR('#b28499') + ", " + hexToG('#b28499') + ", " + hexToB('#b28499') + ", " + 0.5 + ")",
	      //"rgba(" + hexToR('#cc7e78') + ", " + hexToG('#cc7e78') + ", " + hexToB('#cc7e78') + ", " + 0.5 + ")",
	      //"rgba(" + hexToR('#de765b') + ", " + hexToG('#de765b') + ", " + hexToB('#de765b') + ", " + 0.5 + ")",
	      //"rgba(" + hexToR('#ec6c42') + ", " + hexToG('#ec6c42') + ", " + hexToB('#ec6c42') + ", " + 0.5 + ")",
	      //"rgba(" + hexToR('#f55f2c') + ", " + hexToG('#f55f2c') + ", " + hexToB('#f55f2c') + ", " + 0.5 + ")",
	      //"rgba(" + hexToR('#fb4f17') + ", " + hexToG('#fb4f17') + ", " + hexToB('#fb4f17') + ", " + 0.5 + ")",
	      //"rgba(" + hexToR('#fe3705') + ", " + hexToG('#fe3705') + ", " + hexToB('#fe3705') + ", " + 0.5 + ")",
	      //"rgba(" + hexToR('#ff0000') + ", " + hexToG('#ff0000') + ", " + hexToB('#ff0000') + ", " + 0.5 + ")"
	      //
	      "rgba(" + hexToR('#00ffff') + ", " + hexToG('#00ffff') + ", " + hexToB('#00ffff') + ", " + 0.5 + ")",
	      "rgba(" + hexToR('#64f0ff') + ", " + hexToG('#64f0ff') + ", " + hexToB('#64f0ff') + ", " + 0.5 + ")",
	      "rgba(" + hexToR('#87e1ff') + ", " + hexToG('#87e1ff') + ", " + hexToB('#87e1ff') + ", " + 0.5 + ")",
	      "rgba(" + hexToR('#a0d0ff') + ", " + hexToG('#a0d0ff') + ", " + hexToB('#a0d0ff') + ", " + 0.5 + ")",
	      "rgba(" + hexToR('#b5c0ff') + ", " + hexToG('#b5c0ff') + ", " + hexToB('#b5c0ff') + ", " + 0.5 + ")",
	      "rgba(" + hexToR('#c6adff') + ", " + hexToG('#c6adff') + ", " + hexToB('#c6adff') + ", " + 0.5 + ")",
	      "rgba(" + hexToR('#d49bff') + ", " + hexToG('#d49bff') + ", " + hexToB('#d49bff') + ", " + 0.5 + ")",
	      "rgba(" + hexToR('#e185ff') + ", " + hexToG('#e185ff') + ", " + hexToB('#e185ff') + ", " + 0.5 + ")",
	      "rgba(" + hexToR('#ec6dff') + ", " + hexToG('#ec6dff') + ", " + hexToB('#ec6dff') + ", " + 0.5 + ")",
	      "rgba(" + hexToR('#ff1edb') + ", " + hexToG('#ff1edb') + ", " + hexToB('#ff1edb') + ", " + 0.5 + ")"
	    ]
	    
	    //var result = [];
	    //for (var j = 225; j >= 100; j = j - step) {
	    //  result.push(asColorStyle(j, j, j, 1));
	    //}
	    
	    result.indexFor = function(m) {  // map wind speed to a style
	        return Math.floor(Math.min(m, maxWind) / maxWind * (result.length - 1));
	    };
	    return result;
	}

	var colorStyles = windIntensityColorScale(INTENSITY_SCALE_STEP, MAX_WIND_INTENSITY);
	var buckets = colorStyles.map(function() { return []; });

	var particleCount = Math.round(bounds.width * bounds.height * PARTICLE_MULTIPLIER);
	if (isMobile()) {
	  particleCount *= PARTICLE_REDUCTION;
	}

	var fadeFillStyle = "rgba(0, 0, 0, 0.97)";

	var particles = [];
	for (var i = 0; i < particleCount; i++) {
	    particles.push(field.randomize({age: Math.floor(Math.random() * MAX_PARTICLE_AGE) + 0}));
	}

	function evolve() {
	    buckets.forEach(function(bucket) { bucket.length = 0; });
	    particles.forEach(function(particle) {
	        if (particle.age > MAX_PARTICLE_AGE) {
	            field.randomize(particle).age = 0;
	        }
	        var x = particle.x;
	        var y = particle.y;
	        var v = field(x, y);  // vector at current position
	        var m = v[2];
	        if (m === null) {
	            particle.age = MAX_PARTICLE_AGE;  // particle has escaped the grid, never to return...
	        }
	        else {
	            var xt = x + v[0];
	            var yt = y + v[1];
	            if (field(xt, yt)[2] !== null) {
	                // Path from (x,y) to (xt,yt) is visible, so add this particle to the appropriate draw bucket.
	                particle.xt = xt;
	                particle.yt = yt;
	                buckets[colorStyles.indexFor(m)].push(particle);
	            }
	             else {
	                // Particle isn't visible, but it still moves through the field.
	                particle.x = xt;
	                particle.y = yt;
	            }
	        }
	        particle.age += 1;
	    });
	}

	var g = params.canvas.getContext("2d");
	g.lineWidth = PARTICLE_LINE_WIDTH;
	g.fillStyle = fadeFillStyle;

	function draw() {
	    // Fade existing particle trails.
	    var prev = g.globalCompositeOperation;
	    g.globalCompositeOperation = "destination-in";
	    g.fillRect(bounds.x, bounds.y, bounds.width, bounds.height);
	    g.globalCompositeOperation = prev;

	    // Draw new particle trails.
	    buckets.forEach(function(bucket, i) {
	        if (bucket.length > 0) {
	            g.beginPath();
	            g.strokeStyle = colorStyles[i];
	            bucket.forEach(function(particle) {
	                g.moveTo(particle.x, particle.y);
	                g.lineTo(particle.xt, particle.yt);
	                particle.x = particle.xt;
	                particle.y = particle.yt;
	            });
	            g.stroke();
	        }
	    });
	}


	(function frame() {
	    try {
	        windy.timer = setTimeout(function() {
	          requestAnimationFrame(frame);
	          evolve();
	          draw();
	        }, 1000 / FRAME_RATE);
	    }
	    catch (e) {
	        console.error(e);
	    }
	})();

	}
*/

	var drawField = function(mapBounds, bounds, grid){

		//console.log(bounds);
		var img;
		var ctx = params.canvas.getContext('2d');
		img = ctx.getImageData(0,0,grid.rows, grid.colms);
		

		for(var i =0;i<grid.rows;i++)
			for(var j =0;j<grid.colms;j++){

				var scale = grid.Idata[i][j]/255.0;
				var gr = grid.gridData(i,j);
				var mag = Math.sqrt(gr[0]*gr[0] + gr[1]*gr[1]);
				var magind = mag/MAX_VECTOR_MAG;
				if(i ===0 && j=== 0) console.log(magind);
				var kn = (NUM_COLORS * (magind))>>1;
				//if(i === 0 && j === 1) console.log(kn);
				var r,g,b;

				var p1 = point[kn];
				//console.log(p1);
				//var p2 = point[kn+1];
				var p2;
				if(kn == 3)
					p2 = point[kn-1];
				else
					p2 = point[kn+1];

				r = p1.getX()*(1-(NUM_COLORS * magind)+kn)*scale + p2.getX()*((NUM_COLORS * magind)-kn)*scale;
				g = p1.getY()*(1-(NUM_COLORS * magind)+kn)*scale + p2.getY()*((NUM_COLORS * magind)-kn)*scale;
				b = p1.getZ()*(1-(NUM_COLORS * magind)+kn)*scale + p2.getZ()*((NUM_COLORS * magind)-kn)*scale;

				img.data[(i * grid.colms) + j] =  r;//Math.random() >= 0.5 ? 255 : 0;//(grid.Idata[i][j]/255.0) * 40;
				img.data[(i * grid.colms) + j + 1] =  g;//Math.random() >= 0.5 ? 255 : 0;//(grid.Idata[i][j]/255.0)* 40;
				img.data[(i * grid.colms) + j + 2] =  b;//Math.random() >= 0.5 ? 255 : 0;//(grid.Idata[i][j]/255.0)* 40;
				img.data[(i * grid.colms) + j + 3] =  255;
			}


		console.log(grid.Idata[100][1]/255.0, grid.Idata[111][2]/255.0, grid.Idata[123][1]/255.0, grid.Idata[10][1]/255.0);
		
		//console.log(imgData[100 * grid.rows+1], imgData[111* grid.rows+1], imgData[123* grid.rows+1], imgData[10* grid.rows+1]);


		//var g = params.canvas.getContext("2d");
		//g.drawImage(genNoise(grid.rows, grid.colms), 0,0);
		
		/*Put the scaled image on the canvas*/
		//ctx.putImageData(img, 0,0);
		var newCanvas = $("<canvas>")
		    .attr("width", grid.colms)
		    .attr("height",  grid.rows)[0];

		newCanvas.getContext("2d").putImageData(img, 0, 0);
		//Scale the grid to the whole
		ctx.scale( 1.8 * bounds[1][0]/grid.colms, 2 * bounds[1][1]/grid.rows);
		ctx.drawImage(newCanvas, 0, 0);




		//ctx = params.canvas.getContext("2d");
		//ctx.globalCompositeOperation = "destination-in";
		//ctx.fillStyle = "green";
		//ctx.fillRect(bounds[0][0], bounds[0][1], bounds[1][0], bounds[1][1]);

		// save canvas image as data url (png format by default)
      	//var dataURL = params.canvas.toDataURL();

		// set canvasImg image src to dataURL
		// so it can be saved as an image
      	//document.getElementById('canvasImg').src = dataURL;


	};

	var start = function( bounds, width, height, extent ){

		var mapBounds = {
			south: deg2rad(extent[0][1]),
			north: deg2rad(extent[1][1]),
			east: deg2rad(extent[1][0]),
			west: deg2rad(extent[0][0]),
			width: width,
			height: height
		};
		//console.log(mapBounds);
		stop();
		// build grid, hitdata, texdata, idata
		buildGrid( params.data, function(grid){
			/*
			interpolateField( grid, buildBounds( bounds, width, height), mapBounds, function( bounds, field ){
				//console.log("interpolate field");
				// animate the canvas with random points
				console.log("Bounds: ",bounds.x, bounds.y, bounds.width, bounds.height);
				windy.field = field;
				animate( bounds, field );
			});
			*/
			//console.log(grid);
			lic(grid, function(){
				//Draw Field
				drawField(mapBounds, bounds ,grid);
			});

		});
	};

	var stop = function(){
		//if (windy.field) windy.field.release();
		//if (windy.timer) clearTimeout(windy.timer)
	};


	var windy = {
		params: params,
		start: start,
		stop: stop
	};

	return windy;
}



// shim layer with setTimeout fallback
window.requestAnimationFrame = (function(){
  return  window.requestAnimationFrame       ||
          window.webkitRequestAnimationFrame ||
          window.mozRequestAnimationFrame    ||
          window.oRequestAnimationFrame ||
          window.msRequestAnimationFrame ||
          function( callback ){
            window.setTimeout(callback, 1000 / 20);
          };
})();
