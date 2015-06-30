/******************\
|     Deciduous    |
|     Cover Art    |
| @author Anthony  |
| @version 0.3     |
| @date 2015/06/27 |
| @edit 2015/06/30 |
\******************/

var DeciduousCoverArt = (function() {
  /**********
   * config */
  var DRAW_HELPERS = {
      coordSystems: false,
      axes: false,
      contourFunction: false,
      generatorPoints: false,
      boundingBox: false,
      booleanColors: true,
      colorThresh: 0.19
  };

  var cDIMS = [500, 500]; //size in pixels
  var mDIMS = [2, 2]; //size in mathematical units
  var mBOUND_BOX = [
    [-1, 1],
    [1, -1]
  ]; //top left and bottom right corner in mathematical coords

  var BG = {
      numRects: 4,
      theta: -20*Math.PI/180,
      maxOpacity: 0.07,
      colorIntensities: [
          [0, 1, 3, 1],
          [1, 2, 3, 0]
      ]
  }; //the background

  var A = 0.4, C = -0.3; //parameters of the contour function
  var points = [
    [1.21, 0.19], [-1.21, 0.19],
    [0.49, -0.26], [-0.49, -0.26],
    [1.0, -0.3], [-1.0, -0.3],
    [0, -0.58],
    [1.15, -0.69], [-1.15, -0.69],
    [0.54, -0.74], [-0.54, -0.74],
    [0, -1.03],
    [0.65, -1.10], [-0.65, -1.10]
  ]; //these points spawn leaves

  /*************
   * constants */
  var cORIGIN = cDIMS.map(function(a) {
    return a/2;
  }); //the canvas point that corresponds to the mathematical (0,0)
  var mORIGIN = [0, 0]; //mathematical origin

  //some helper variables for the background
  BG.rW = (
      Math.abs(
          Math.sin(BG.theta)*cDIMS[1]
      )+Math.abs(Math.cos(BG.theta)*cDIMS[0])
  )/BG.numRects;
  BG.rL = Math.abs(
      cDIMS[1]/Math.cos(BG.theta)
  ) + Math.abs(BG.rW*Math.tan(BG.theta));
  BG.offset = cDIMS[1]*Math.tan(BG.theta);

  /*********************
   * working variables */
  var canvas, ctx;
  var contFunc;

  /******************
   * work functions */
  function initDeciduousCoverArt() {
    //assign working variables
    canvas = $s('#canvas');
    canvas.width = cDIMS[0];
    canvas.height = cDIMS[1];
    ctx = canvas.getContext('2d');

    contFunc = getContourFunction([A, C]);
    $s('#color-thresh').value = DRAW_HELPERS.colorThresh;

    //event listeners
    canvas.addEventListener('click', function(e) {
        var pos = getCanvMousePos(e);
        var mPos = cToM(pos);
        if (Math.abs(mPos[0]) < 0.02) {
            mPos[0] = 0;
            points.push(mPos);
        } else {
            points.push(mPos);
            points.push([-mPos[0], mPos[1]]);
        }
        drawCoverArt();
    });

    $s('#toggle-color-btn').addEventListener('click', function(e) {
        DRAW_HELPERS.booleanColors = !DRAW_HELPERS.booleanColors;
        drawCoverArt();
    });

    $s('#color-thresh-btn').addEventListener('click', function() {
        DRAW_HELPERS.colorThresh = parseFloat($s('#color-thresh').value);
        drawCoverArt();
    });

    $s('#save-btn').addEventListener('click', promptSaveCanvas);

    //do the work
    drawCoverArt();
  }

  function drawCoverArt() {
      //clean slate
      clearCanvas();

      //draw the background
      ctx.translate(BG.offset, 0);
      ctx.rotate(BG.theta);
      for (var ai = 0; ai < BG.numRects; ai++) {
          var opcty = BG.maxOpacity*(BG.colorIntensities[0][ai]+1)/BG.numRects;
          ctx.fillStyle = 'rgba(0,0,0,'+opcty+')';
          ctx.fillRect(
              ai*BG.rW,
              -ai*BG.rW*Math.tan(BG.theta),
              BG.rW, BG.rL
          );
      }
      ctx.rotate(-BG.theta);
      ctx.translate(-BG.offset, 0);

      ctx.rotate(-BG.theta);
      for (var ai = 0; ai < BG.numRects; ai++) {
          var opcty = BG.maxOpacity*(BG.colorIntensities[1][ai]+1)/BG.numRects;
          ctx.fillStyle = 'rgba(0,0,0,'+opcty+')';
          ctx.fillRect(
              ai*BG.rW,
              (ai+1)*BG.rW*Math.tan(BG.theta),
              BG.rW, BG.rL
          );
      }
      ctx.rotate(BG.theta);

      //draw the distance-based colors
      var s = +new Date();
      paintDistances();
      console.log(((+new Date()) - s) + 'ms');

      //draw the axes
      if (DRAW_HELPERS.axes) {
          drawAxes('#CCC', 1);
      }

      //outlines the viewport with a box
      if (DRAW_HELPERS.boundingBox) {
          var topLeft = mToC(mBOUND_BOX[0]);
          var botRight = mToC(mBOUND_BOX[1]);
          ctx.strokeStyle = 'black';
          ctx.strokeRect(
            topLeft[0], topLeft[1],
            botRight[0]-topLeft[0], botRight[1] - topLeft[1]
          );
      }

      //draw the generator points
      for (var ai = 0; ai < points.length; ai++) {
        var p = points[ai];
        if (DRAW_HELPERS.generatorPoints) {
            drawPoint(mToC(p), 1, 'red');
        }

        if (DRAW_HELPERS.coordSystems) {
            //draw the closest normal
            var closestNormal = contFunc.normalThruP(p);
            if (closestNormal(0) === Infinity) {
                drawLine(
                    mToC([0, -mDIMS[1]/2]),
                    mToC([0, mDIMS[1]/2]),
                    'orange', 1
                );
            }
            drawLineFromFunc(closestNormal, 'orange', 1);

            //draw the tan line
            var tanFunc = contFunc.tanThruP(p);
            drawLineFromFunc(tanFunc, 'orange', 1);
        }
      }

      //draw F(X)
      if (DRAW_HELPERS.contourFunction) {
          plot(
              contFunc.F,
              mORIGIN[0]-mDIMS[0]/2, mORIGIN[0]+mDIMS[0]/2,
              250, 'red'
          );
      }
  }

  function paintDistances() {
    var coordSystems = points.map(function(g) {
        return contFunc.getLeafVectors(g);
    });
    var distFuncs = points.map(function(g) {
        return getDistFunc(g);
    });
    var farthest = 0;
    var allTheNthClosests = [];
    for (var y = 0; y < canvas.height; y++) {
        var thisRowsClosests = [];
        for (var x = 0; x < canvas.width; x++) {
        	var distances = [];
            var mCoords = cToM([x, y]);
        	for (var ai = 0; ai < points.length; ai++) {
        		distances.push(
                    [ai, distFuncs[ai](mCoords)]
                );
        	}
        	distances.sort(function(a,b) { return a[1]-b[1]; });
        	var nthClosest = [
                distances[0][0],
                distances[0][1] - 0.05*distances[1][1]
            ];
        	thisRowsClosests.push(nthClosest);
        	if (nthClosest[1] > farthest) farthest = nthClosest[1];
        }
        allTheNthClosests.push(thisRowsClosests);
    }

    var currImageData = ctx.getImageData(0, 0, canvas.width, canvas.height);
    for (var y = 0; y < canvas.height; y++) {
    	for (var x = 0; x < canvas.width; x++) {
    		var b = [0, 0, 0, 0];

            //black and white
            if (DRAW_HELPERS.booleanColors) {
                var frac = allTheNthClosests[y][x][1]/farthest;
                if (frac < DRAW_HELPERS.colorThresh) {
                    var globalCoords = cToM([x, y]);
                    var transformedCoords = getCoordsIn(
                        cToM([x, y]),
                        coordSystems[allTheNthClosests[y][x][0]]
                    );
                    color = getLeafColor(
                        transformedCoords, globalCoords,
                        coordSystems[allTheNthClosests[y][x][0]],
                        allTheNthClosests[y][x][1],
                        farthest
                    );
                } else {
                    continue; //skip this pixel
                }
            } else { //cool red to blue kinda thing
                color = getCoolColor(
                    allTheNthClosests[y][x][1],
                    [0, farthest]
                );
            }

    		var idx = 4*(y*canvas.width + x);
    		currImageData.data[idx+0] = color[0];
    		currImageData.data[idx+1] = color[1];
    		currImageData.data[idx+2] = color[2];
    		currImageData.data[idx+3] = color[3];
    	}
    }
    ctx.putImageData(currImageData, 0, 0);
  }

  /* getLeafColor(localCoords, globalCoords, coordSystem, dist, farthest)
   * Given all the defining information about a point in mathematical space,
   * this function returns its color.
   */
  function getLeafColor(
      localCoords, globalCoords, coordSystem, dist, farthest
  ) {
      var frac = dist/farthest;
      var color1 = HSVtoRGB(
          numMap(
              (globalCoords[0]<mORIGIN[0]?1:-1)*localCoords[0],
              [-0.4375, 0.4375], [0,0.8]
          ), 0.7, 1
      ); //rainbow gradient
      color1 = getGradient(
          color1,
          [255, 255, 255],
          (frac > DRAW_HELPERS.colorThresh - 0.01 ? numMap(
              dist/farthest, [
                  DRAW_HELPERS.colorThresh - 0.01, DRAW_HELPERS.colorThresh
              ], [0.9, 0.1]
          ) : 0.9)
      ); //white highlights
      if (localCoords[1] > 0) {
          return getGradient(
              color1,
              [0, 0, 0],
              (frac > DRAW_HELPERS.colorThresh - 0.02 ? numMap(
                  dist/farthest, [
                      DRAW_HELPERS.colorThresh - 0.02, DRAW_HELPERS.colorThresh
                  ], [1, 0.85]
              ) : 1)
          ).concat([255]); //black shadows
      } else return color1.concat([255]);
  }

  /* getDistFunc(g)
   * Given a generator point, this function returns the distance function that
   * defines how far away points are from this generator.
   */
  function getDistFunc(g) {
      var coordSystem = contFunc.getLeafVectors(g);
      return function(p) {
          var transformedCoords = getCoordsIn(p, coordSystem);
          var x = Math.abs(transformedCoords[0]);
          var y = Math.abs(transformedCoords[1]);
          var params = transformedCoords[1] > 0 ? [
              //spikeyContr, roundyContr, ...
              0.3, tightNumMap(y, [0, 0.4], [0.7, 2]), 0.9, 0.4, 2.8, 2.0
          ] : [
              0.3, tightNumMap(y, [0, 0.4], [0.7, 2]), 0.9, 0.4, 2.8, 2.0
          ];
          var spikey = Math.sqrt(
              Math.pow(x, params[2]) + Math.pow(y, params[3])
          ); //pointy
          var roundy = Math.sqrt(
              Math.pow(x, params[4]) + Math.pow(y, params[5])
          ); //rounded
          var dist = params[0]*spikey + params[1]*roundy;
          return Math.pow(dist, 0.85);
      };
  }

  /***********
   * objects */
  function getContourFunction(params) {
    return {
      /* F(x)
       * This function defines the shape of the leaf group. It'll be referred to
       * as the "contour" function.
       */
      F: function(x) {
        return params[0]*Math.pow(x, 2) + params[1];
      },

      /* dFdx(x)
       * This is the derivative of F(x).
       */
      dFdx: function(x) {
        return params[0]*x;
      },

      /* normalsThruP(p)
       * Returns all the lines through point p that are normal to F(x).
       */
      normalsThruP: function(p) {
        var k = (1/(2*params[0]))+(-p[1])+(params[1]);
        //potential x-coords of intersections between the normal and F(x)
        var x0s = solveCubic(params[0], 0, k, -p[0]/(2*params[0]));

        //create the normal functions from the x0s
        var funcs = [];
        for (var ai = 0; ai < x0s.length; ai++) {
          funcs.push((function(x0, CF) {
            return function(x) {
              if (x === 'get-intersection-x') {
                return x0; //special case for efficiency later
              } else if (x0 === 0) {
                  return Infinity;
              } else {
                return (-1/(2*params[0]*x0))*(x-x0) + CF(x0);
              }
            }
          })(x0s[ai], this.F));
        }

        return funcs;
      },

      /* normalThruP(p)
       * Returns the normal through p with the closest intersection with F.
       */
      normalThruP: function(p) {
        var normals = this.normalsThruP(p);

        //we need to find the closest such intersection
        var x0s = normals.map(function(normal) {
          return normal('get-intersection-x');
        });
        var x0 = -1, minDist = Infinity;
        for (var ai = 0; ai < x0s.length; ai++) {
          var cand = x0s[ai];
          var dist = getDist(p, [cand, this.F(cand)]);
          if (dist < minDist) {
            x0 = cand;
            minDist = dist;
          }
        }

        //return the normal that corresponded to the closest intersection
        for (var ai = 0; ai < normals.length; ai++) {
          if (x0 === normals[ai]('get-intersection-x')) {
            return normals[ai];
          }
        }
      },

      /* tanThruP(p)
       * Returns the perpendicular line to normalThruP(p) that goes through p.
       */
      tanThruP: function(p) {
          var theNormThruP = this.normalThruP(p);
          var x0 = theNormThruP('get-intersection-x');
          var m = 2*params[0]*x0;
          return function(x) {
              return m*(x - p[0]) + p[1];
          };
      },

      /* getLeafVectors(p)
       * Returns the coordinate system that defines the distance function of
       * each generator point.
       */
      getLeafVectors: function(p) {
          var tang = this.tanThruP(p);
          var norm = this.normalThruP(p);
          var xAxis = normalize([1, tang(1)-tang(0)]);
          var yAxis = norm(0) === Infinity ? [0, 1] : normalize(
              [1, norm(1)-norm(0)]
          );
          return {
              x: xAxis,
              y: p[0] > 0 ? yAxis : [-yAxis[0], -yAxis[1]],
              origin: p
          };
      }
    };
  }

  /********************
   * helper functions */
  function promptSaveCanvas() {
      var downloadUrl = canvas.toDataURL('image/png').replace('image/png', 'image/octet-stream');
      if (downloadUrl.length < 1024*1024 || confirm(
          'This is a relatively large image, so your browser may crash. Continue?'
      )) {
          window.location = downloadUrl;
      }
  }

  function getJuliaColorFromCoord(realConst, imConst, initX, initY) {
    var colorMult = 2;
    var maxIterations = 400;
    var palette = [
        [255,255,255], [254,254,254], [255, 255, 0], [0, 255, 0], [0, 0, 255],
        [255, 0, 255], [255, 0, 0], [255,255,255], [254,254,254]
    ];
    var color = [0, 0, 0];

    ///////////////////////////////////////////
    //continue with the escape time algorithm//
    var x_ = initX, y_ = initY;
    var iteration = 0;
    var xsq = initX*initX, ysq = initY*initY;
    var val = xsq + ysq;

    while (val <= 4 && iteration < maxIterations) {
        y_ = x_*y_;
        y_ += y_; //times 2
        y_ += imConst;
        x_ = (xsq - ysq) + realConst;

        xsq = x_*x_;
        ysq = y_*y_;
        val = xsq + ysq;
        iteration += 1;
    }

    if (iteration != maxIterations) { //if it didn't survive all the iterations, it has a color
        var mu = iteration - (Math.log(Math.log(val))); //fractional stopping iteration
        mu = palette.length * (mu/maxIterations); //spread the colors out
        if (mu > palette.length) mu = palette.length; //not too big
        else if (mu < 0) mu = 0; //not too small
        mu *= colorMult; //allow the color scheme to repeat
        var intPartMu = Math.floor(mu); //integer part of mu
        var bucket = intPartMu%palette.length; //base the color on the integer part of mu
        var nextBucket = (bucket+1)%palette.length; //the color right after that
        var percentToNextBucket = mu - intPartMu; //the fractional part of mu

        color = getGradient(
            palette[bucket], palette[nextBucket],
            1-percentToNextBucket
        ); //invert it because small fractions means more of the first color
    }

    return [color[0], color[1], color[2], 255];
  }

  /* getGradient(c1, c2, percent)
   * Returns an RGB color that's a percentage of the first and the second
   */
  function getGradient(c1, c2, percent) {
  	var ret = [0, 0, 0];

  	ret[0] = Math.floor((percent * c1[0]) + ((1 - percent) * c2[0]))%256;
  	ret[1] = Math.floor((percent * c1[1]) + ((1 - percent) * c2[1]))%256;
  	ret[2] = Math.floor((percent * c1[2]) + ((1 - percent) * c2[2]))%256;

  	return ret;
  }

  function plot(fn, left, right, n, color) {
      for (var x = left; x < right; x += (right - left)/n) {
        var p = mToC([x, fn(x)]);
        drawPoint(p, 1, color);
      }
  }

  function drawAxes(color, thickness) {
    drawLine(
      [0, cORIGIN[1]], [cDIMS[0], cORIGIN[1]], color, thickness
    ); //x axis
    drawLine(
      [cORIGIN[0], 0], [cORIGIN[0], cDIMS[1]], color, thickness
    ); //y axis
  }

  function clearCanvas() {
    ctx.fillStyle = 'white';
    ctx.fillRect(0, 0, canvas.width, canvas.height);
  }

  function drawPoint(pos, r, color) {
    ctx.fillStyle = color || 'rgba(255, 255, 255, 0.3)';
    ctx.beginPath();
    ctx.arc(pos[0], pos[1], r, 0, 2*Math.PI, true);
    ctx.closePath();
    ctx.fill();
  }

  function drawLine(start, end, color, thickness) {
    ctx.strokeStyle = color || 'rgba(0, 0, 0, 1)';
    ctx.beginPath();
    ctx.moveTo(start[0], start[1]);
    ctx.lineTo(end[0], end[1]);
    ctx.lineWidth = thickness || 3;
    ctx.stroke();
  }

  function drawLineFromFunc(fn, color, thickness) {
      //draw the tan line
      var startX = mORIGIN[0]-mDIMS[0]/2;
      var startY = fn(startX);
      var endX = mORIGIN[0]+mDIMS[0]/2;
      var endY = fn(endX);
      drawLine(
        mToC([startX, startY]), mToC([endX, endY]),
        color || 'orange', thickness || 1
      );
  }

  /* mToC(p)
   * This function converts mathematical coordinates to canvas coordinates.
   */
  function mToC(p) {
    return [
      (p[0] - mORIGIN[0])*(cDIMS[0]/mDIMS[0]) + cORIGIN[0],
      -(p[1] - mORIGIN[1])*(cDIMS[1]/mDIMS[1]) + cORIGIN[1],
    ];
  }

  /* cToM(p)
   * This function converts canvas coordinates to mathematical coordinates.
   */
  function cToM(p) {
    return [
      (p[0] - cORIGIN[0])*(mDIMS[0]/cDIMS[0]) + mORIGIN[0],
      -(p[1] - cORIGIN[1])*(mDIMS[1]/cDIMS[1]) + mORIGIN[1],
    ];
  }

  /* getCoordsIn(p, coordSystem)
   * Returns the coordinates of p in the 2d system described by coordSystem.
   */
  function getCoordsIn(p, coordSystem) {
      var shiftedP = [
          p[0] - coordSystem.origin[0], p[1] - coordSystem.origin[1]
      ];
      return [
          getProjOn(shiftedP, coordSystem.x),
          getProjOn(shiftedP, coordSystem.y)
      ];
  }

  /* getProjOn(a, b)
   * Returns the projection of vector a on b, where a and b are 2d vectors.
   */
  function getProjOn(a, b) {
    return (a[0]*b[0] + a[1]*b[1])/getMag(b);
  }

  function getCanvMousePos(e) {
    var rect = canvas.getBoundingClientRect();
    return [e.clientX-rect.left, e.clientY-rect.top];
  }

  function $s(id) { //for convenience
    if (id.charAt(0) !== '#') return false;
    return document.getElementById(id.substring(1));
  }

  function getRandInt(low, high) { //output is in [low, high)
    return Math.floor(low + Math.random()*(high-low));
  }

  function normalize(a) {
      var mag = getMag(a);
      return [a[0]/mag, a[1]/mag];
  }

  function getDist(a, b) {
    return getMag([
        a[0] - b[0],
        a[1] - b[1]
    ]);
  }

  function getMag(a) {
      return Math.sqrt(a[0]*a[0] + a[1]*a[1]);
  }

  function round(n, places) {
    var mult = Math.pow(10, places);
    return Math.round(mult*n)/mult;
  }

  /* FROM http://stackoverflow.com/questions/7706339/
                 grayscale-to-red-green-blue-matlab-jet-color-scale */
  function getCoolColor(n, range) {
  	var raw = [1.0, 1.0, 1.0]; //white

  	if (n < range[0]) n = range[0];
  	if (n > range[1]) n = range[1];
  	var dn = range[1] - range[0];

  	if (n < (range[0] + 0.25 * dn)) {
  		raw[0] = 0;
  		raw[1] = 4 * (n - range[0]) / dn;
  	} else if (n < (range[0] + 0.5 * dn)) {
  		raw[0] = 0;
  		raw[2] = 1 + 4 * (range[0] + 0.25 * dn - n) / dn;
  	} else if (n < (range[0] + 0.75 * dn)) {
  		raw[0] = 4 * (n - range[0] - 0.5 * dn) / dn;
  		raw[2] = 0;
  	} else {
  		raw[1] = 1 + 4 * (range[0] + 0.75 * dn - n) / dn;
  		raw[2] = 0;
  	}

  	var color = [
          tightNumMap(raw[0], [0, 1], [0, 255]),
          tightNumMap(raw[1], [0, 1], [0, 255]),
          tightNumMap(raw[2], [0, 1], [0, 255]),
          255 //alpha
      ];
  	return color;
  }

  /* http://stackoverflow.com/questions/17242144/
            javascript-convert-hsb-hsv-color-to-rgb-accurately */
  function HSVtoRGB(h, s, v) {
    var r, g, b, i, f, p, q, t;
    i = Math.floor(h * 6);
    f = h * 6 - i;
    p = v * (1 - s);
    q = v * (1 - f * s);
    t = v * (1 - (1 - f) * s);
    switch (i % 6) {
        case 0: r = v, g = t, b = p; break;
        case 1: r = q, g = v, b = p; break;
        case 2: r = p, g = v, b = t; break;
        case 3: r = p, g = q, b = v; break;
        case 4: r = t, g = p, b = v; break;
        case 5: r = v, g = p, b = q; break;
    }
    return [
        Math.floor(r * 255),
        Math.floor(g * 255),
        Math.floor(b * 255)
    ];
  }

  /* FROM http://stackoverflow.com/questions/27176423/
                 function-to-solve-cubic-equation-analytically */
  function solveCubic(a, b, c, d) {
    function cuberoot(x) {
      var y = Math.pow(Math.abs(x), 1/3);
      return x < 0 ? -y : y;
    }

    if (Math.abs(a) < 1e-8) { // Quadratic case, ax^2+bx+c=0
      a = b; b = c; c = d;
      if (Math.abs(a) < 1e-8) { // Linear case, ax+b=0
        a = b; b = c;
        if (Math.abs(a) < 1e-8) return [];// Degenerate case

        return [-b/a];
      }

      var D = b*b - 4*a*c;
      if (Math.abs(D) < 1e-8) [-b/(2*a)];
      else if (D > 0) {
        return [(-b+Math.sqrt(D))/(2*a), (-b-Math.sqrt(D))/(2*a)];
      }

      return [];
    }

    // Convert to depressed cubic t^3+pt+q = 0 (subst x = t - b/3a)
    var p = (3*a*c - b*b)/(3*a*a);
    var q = (2*b*b*b - 9*a*b*c + 27*a*a*d)/(27*a*a*a);
    var roots;

    if (Math.abs(p) < 1e-8) { // p = 0 -> t^3 = -q -> t = -q^1/3
      roots = [cuberoot(-q)];
    } else if (Math.abs(q) < 1e-8) { // q = 0 -> t^3 + pt = 0 -> t(t^2+p)=0
      roots = [0].concat(p < 0 ? [Math.sqrt(-p), -Math.sqrt(-p)] : []);
    } else {
      var D = q*q/4 + p*p*p/27;
      if (Math.abs(D) < 1e-8) {       // D = 0 -> two roots
        roots = [-1.5*q/p, 3*q/p];
      } else if (D > 0) {             // Only one real root
        var u = cuberoot(-q/2 - Math.sqrt(D));
        roots = [u - p/(3*u)];
      } else {                        // D < 0, three roots, but needs to use complex numbers/trigonometric solution
        var u = 2*Math.sqrt(-p/3);
        var t = Math.acos(3*q/p/u)/3;  // D < 0 implies p < 0 and acos argument in [-1..1]
        var k = 2*Math.PI/3;
        roots = [u*Math.cos(t), u*Math.cos(t-k), u*Math.cos(t-2*k)];
      }
    }

    // Convert back from depressed cubic
    for (var i = 0; i < roots.length; i++) roots[i] -= b/(3*a);

    return roots;
  }

  /* tightNumMap(n, d, r)
   * Returns numMap(n, d, r), bounding the output to r.
   */
  function tightNumMap(n, d, r) { //enforces boundaries
  	var raw = numMap(n, d, r);
  	if (raw < r[0]) return r[0];
  	else if (raw > r[1]) return r[1];
  	else return raw;
  }

  /* numMap(n, d, r)
   * Given an n in [d[0], d[1]], this function returns a linearly related
   * number in [r[0], r[1]]
   */
  function numMap(n, d, r) {
  	var Rd = d[1]-d[0];
  	var Rr = r[1]-r[0];
  	return (Rr/Rd)*(n - d[0]) + r[0];
  }

  return {
    init: initDeciduousCoverArt,
    getCoordsIn: getCoordsIn,
    getPoints: function() {
        return points;
    }
  };
})();

window.addEventListener('load', DeciduousCoverArt.init);
