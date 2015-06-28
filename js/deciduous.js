/******************\
|     Deciduous    |
|     Cover Art    |
| @author Anthony  |
| @version 0.1     |
| @date 2015/06/27 |
| @edit 2015/06/27 |
\******************/

var DeciduousCoverArt = (function() {
  /**********
   * config */
  var DRAW_HELPERS = {
      coordSystems: true,
      axes: true,
      contourFunction: true,
      boundingBox: true
  };
  var cDIMS = [500, 500]; //size in pixels
  var mDIMS = [3, 3]; //size in mathematical units
  var mBOUND_BOX = [
    [-1, 1],
    [1, -1]
  ]; //top left and bottom right corner in mathematical coords
  var A = 0.75, C = -0.5; //parameters of the contour function
  var points = [
    [0.5, -0.65],
    [-0.5, -0.65]
  ]; //these points spawn leaves

  /*************
   * constants */
  var cORIGIN = cDIMS.map(function(a) {
    return a/2;
  }); //the canvas point that corresponds to the mathematical (0,0)
  var mORIGIN = [0, 0]; //mathematical origin

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

    //clean slate
    clearCanvas();

    //draw the distance estimations
    var leftX = mORIGIN[0]-mDIMS[0]/2;
    var rightX = leftX + mDIMS[0];
    var leftY = mORIGIN[1]-mDIMS[1]/2;
    var rightY = leftY + mDIMS[1];
    var detail = 250;
    for (var ai = 0; ai < points.length; ai++) {
        var p = points[ai];
        var distFunc = getDistFunc(p);
        var maxDist = Math.max(
            Math.max(distFunc([leftX, leftY]), distFunc([leftX, rightY])),
            Math.max(distFunc([rightX, leftY]), distFunc([rightX, rightY]))
        );
        for (var x = leftX; x < rightX; x+=mDIMS[0]/detail) {
            for (var y = leftY; y < rightY; y+=mDIMS[1]/detail) {
                var dist = distFunc([x, y]);
                var col = numMap(dist, [0, maxDist], [0, 255]);
                col = Math.floor(col);
                var colString = 'rgba('+col+','+col+','+col+', 1)';
                if (col < 37) {
                    drawPoint(mToC([x, y]), 1, 'black');
                }
            }
        }
    }

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
      drawPoint(mToC(p), 4, '#5b7');

      if (DRAW_HELPERS.coordSystems) {
          //draw the closest normal
          plot(
              contFunc.normalThruP(p),
              mORIGIN[0]-mDIMS[0]/2, mORIGIN[0]+mDIMS[0]/2,
              1000, 'orange'
          );

          //draw the tan line
          plot(
              contFunc.tanThruP(p),
              mORIGIN[0]-mDIMS[0]/2, mORIGIN[0]+mDIMS[0]/2,
              1000, 'orange'
          );
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

  /* getDistFunc(g)
   * Given a generator point, this function returns the distance function that
   * defines how far away points are from this generator.
   */
  function getDistFunc(g) {
      var coordSystem = contFunc.getLeafVectors(g);
      return function(p) {
          var shiftedP = [p[0] - g[0], p[1] - g[1]];
          var transformedCoords = [
              getProjOn(shiftedP, coordSystem.x),
              getProjOn(shiftedP, coordSystem.y)
          ];
          var x = Math.abs(transformedCoords[0]);
          var y = Math.abs(transformedCoords[1]);
          var dist = Math.sqrt(Math.pow(x, 4) + Math.pow(y, 2));
          return Math.pow(dist, 0.5);
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
          //warning: can't handle vertical lines
          var xAxis = normalize([1, tang(1)-tang(0)]);
          var yAxis = normalize([1, norm(1)-norm(0)]);
          return {
              x: xAxis,
              y: yAxis
          };
      },

      /* D(p, q)
       * Returns the leafy-distance from generator g to point p.
       */
      D: function(g, p) {
        return false;
      }
    };
  }

  /********************
   * helper functions */
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

  /* mToC(p)
   * This function converts mathematical coordinates to canvas coordinates.
   */
  function mToC(p) {
    return [
      p[0]*(cDIMS[0]/mDIMS[0]) + cORIGIN[0],
      -p[1]*(cDIMS[1]/mDIMS[1]) + cORIGIN[1],
    ];
  }

  /* getProjOn(a, b)
   * Returns the projection of vector a on b, where a and b are 2d vectors.
   */
  function getProjOn(a, b) {
    return (a[0]*b[0] + a[1]*b[1])/getMag(b);
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
    init: initDeciduousCoverArt
  };
})();

window.addEventListener('load', DeciduousCoverArt.init);
