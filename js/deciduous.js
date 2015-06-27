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
  var mDIMS = [2, 2]; //size in mathematical units
  var points = [
    [0, 0]
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

  /******************
   * work functions */
  function initDeciduousCoverArt() {
    //assign working variables
    canvas = $s('#canvas');
    canvas.width = cDIMS[0];
    canvas.height = cDIMS[1];
    ctx = canvas.getContext('2d');

    //clean slate
    clearCanvas();

    //draw the axes
    drawAxes('#CCC', 1);

    //draw the generator points
    for (var ai = 0; ai < points.length; ai++) {
      var p = points[ai];
      drawPoint(mToC(p), 4, '#5b7');
    }

    //draw F(X)
    var n = 250;
    var leftmostX = mORIGIN[0]-mDIMS[0]/2;
    for (var x = leftmostX; x < leftmostX+mDIMS[0]; x += mDIMS[0]/n) {
      var y = F(x);
      var p = mToC([x, y]);
      drawPoint(p, 1, 'red');
    }
  }

  /* F(x)
   * This function defines the shape of the leaf group.
   */
  function F(x) {
    return 0.75*Math.pow(x, 2) - 0.5;
  }

  /***********
   * objects */

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

  function $s(id) { //for convenience
    if (id.charAt(0) !== '#') return false;
    return document.getElementById(id.substring(1));
  }

  function getRandInt(low, high) { //output is in [low, high)
    return Math.floor(low + Math.random()*(high-low));
  }

  function round(n, places) {
    var mult = Math.pow(10, places);
    return Math.round(mult*n)/mult;
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
