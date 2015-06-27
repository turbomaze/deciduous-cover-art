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
  var cDIMS = [600, 600]; //size in pixels
  var mDIMS = [4, 4]; //size in mathematical units
  var A = 0.75, C = -0.5; //parameters of the contour function
  var points = [
    [0.5, 0.6]
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

    //draw the axes
    drawAxes('#CCC', 1);

    //draw the generator points
    for (var ai = 0; ai < points.length; ai++) {
      var p = points[ai];
      drawPoint(mToC(p), 4, '#5b7');

      //draw the normals through p
      var normsThruP = contFunc.normalsThruP(p);
      for (var bi = 0; bi < normsThruP.length; bi++) {
        var normThruP = normsThruP[bi];
        console.log(normThruP('get-intersection-x'));
        var n = 1000;
        var leftmostX = mORIGIN[0]-mDIMS[0]/2;
        for (var x = leftmostX; x < leftmostX+mDIMS[0]; x += mDIMS[0]/n) {
          var y = normThruP(x);
          var p = mToC([x, y]);
          drawPoint(p, 1, 'orange');
        }
      }
    }

    //draw F(X)
    var n = 250;
    var leftmostX = mORIGIN[0]-mDIMS[0]/2;
    for (var x = leftmostX; x < leftmostX+mDIMS[0]; x += mDIMS[0]/n) {
      var y = contFunc.F(x);
      var p = mToC([x, y]);
      drawPoint(p, 1, 'red');
    }
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
        var k = 2*params[0]*p[1] - params[1] + 1;
        //potential x-coords of intersections between the normal and F(x)
        var x0s = solveCubic(2*Math.pow(params[0], 2), 0, -k, -p[0]);

        //create the normal functions from the x0s
        var funcs = [];
        for (var ai = 0; ai < x0s.length; ai++) {
          funcs.push((function(x0) {
            return function(x) {
              if (x === 'get-intersection-x') {
                return x0; //special case for efficiency later
              } else {
                return (-1/(2*params[0]*x0))*(x-x0) + params[0]*x0*x0 + params[1];
              }
            }
          })(x0s[ai]));
        }

        return funcs;
      },

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
    var bMag = Math.sqrt(b[0]*b[0] + b[1]*b[1]);
    return (a[0]*b[0] + a[1]*b[1])/bMag;
  }

  function $s(id) { //for convenience
    if (id.charAt(0) !== '#') return false;
    return document.getElementById(id.substring(1));
  }

  function getRandInt(low, high) { //output is in [low, high)
    return Math.floor(low + Math.random()*(high-low));
  }

  function getDist(a, b) {
    return Math.sqrt(Math.pow(a[0]-b[0]) + Math.pow(a[1]-b[1]));
  }

  function round(n, places) {
    var mult = Math.pow(10, places);
    return Math.round(mult*n)/mult;
  }

  /* FROM http://stackoverflow.com/questions/27176423/
                 function-to-solve-cubic-equation-analytically */
  function cuberoot(x) {
    var y = Math.pow(Math.abs(x), 1/3);
    return x < 0 ? -y : y;
  }

  /* FROM http://stackoverflow.com/questions/27176423/
                 function-to-solve-cubic-equation-analytically */
  function solveCubic(a, b, c, d) {
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
