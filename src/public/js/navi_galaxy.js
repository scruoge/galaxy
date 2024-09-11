function distance(p1, p2) {
	var dx = p1.x - p2.x,
		dy = p1.y - p2.y,
		dz = p1.z - p2.z;

	return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

function dist2d2(p1, p2) {
	var dx = p1.x - p2.x,
		dy = p1.y - p2.y;

	return dx * dx + dy * dy;
}

function dist2d(p1, p2) {
	return Math.sqrt(dist2d2(p1, p2));
}

function getCrossPoints(p1, r1, p2, r2) {
	var d = dist2d(p1, p2);
	var b = (r2 * r2 - r1 * r1 + d * d) / (2 * d);
	var a = d - b;
	var h = Math.sqrt(r1 * r1 - a * a);
	var c0 = {
		x: p1.x + a / d * (p2.x - p1.x),
		y: p1.y + a / d * (p2.y - p1.y)
	};
	var c1 = {
			x: c0.x + h / d * (p2.y - p1.y),
			y: c0.y - h / d * (p2.x - p1.x)
		},
		c2 = {
			x: c0.x - h / d * (p2.y - p1.y),
			y: c0.y + h / d * (p2.x - p1.x)
		};

	return [c1, c2];
}

function delegate(o, f) {
    var a = [];
    var l = arguments.length;
    for (var i = 2; i < l; i++) a[i - 2] = arguments[i];
    return function () {
        var aP = Array.from(arguments);
        return f.apply(o, aP.concat(a));
    }
}

/*
// tweening
var Tween = {
    tweens: [],
    init: function () {
        requestAnimationFrame(delegate(this, this.onEnterFrame));
        return this;
    },
    camTween: function (cam, func, begin, end, dur, cb, efCb) {
        efCb = efCb || false;
        var _obj = {
            type: 'cam',
            cam: cam,
            func: func,
            begin: begin,
            end: end,
            dur: dur,
            cb: cb,
            efCb: efCb,
            started: null,
            delta: [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]],
            finished: false
        };

        for (var i = 0; i < 4; i++) {
            for (var j = 0; j < 4; j++) {
                _obj.delta[i][j] = _obj.end[i][j] - _obj.begin[i][j];
            }
        }
        this.tweens.push(_obj);
    },
    viewTween: function (cam, func, begin, end, dur, cb, efCb) {
        efCb = efCb || false;
        var _obj = {
            type: 'view',
            cam: cam,
            func: func,
            begin: begin,
            end: end,
            dur: dur,
            cb: cb,
            efCb: efCb,
            started: null,
            delta: {
                start: {
                    a: end.start.a - begin.start.a,
                    r: end.start.r - begin.start.r
                },
                center: {
                    a: end.center.a - begin.center.a,
                    r: end.center.r - begin.center.r,
                    z: end.center.z - begin.center.z
                },
                end: {
                    a: end.end.a - begin.end.a,
                    r: end.end.r - begin.end.r
                }
            },
            finished: false
        };
        if (Math.abs(_obj.delta.center.a) > Math.PI) {
            _obj.delta.center.a -= Math.sign(_obj.delta.center.a) * Math.PI * 2;
            _obj.delta.start.a -= Math.sign(_obj.delta.start.a) * Math.PI * 2;
            _obj.delta.end.a -= Math.sign(_obj.delta.end.a) * Math.PI * 2;
        }

        this.tweens.push(_obj);
    },
    tween: function (obj, prop, func, begin, end, dur, cb, efCb) {
        efCb = efCb || false;
        var _obj = {
            type: 'obj',
            obj: obj,
            prop: prop,
            func: func,
            begin: begin,
            end: end,
            dur: dur,
            cb: cb,
            efCb: efCb,
            started: null,
            delta: end - begin,
            finished: false
        };

        this.tweens.push(_obj);
    },
    onEnterFrame: function (ts) {
        var t = Math.round(ts);
        var tw = [];
        for (var i = 0; i < this.tweens.length; i++) {
            this.updateTween(this.tweens[i], t);
            if (this.tweens[i].efCb) {
                this.tweens[i].efCb();
            }
            if (!this.tweens[i].finished) {
                tw.push(this.tweens[i]);
            } else {
                this.tweens[i].cb();
            }
        }
        this.tweens = tw;
        requestAnimationFrame(delegate(this, this.onEnterFrame));
    },
    updateTween: function (tween, t) {
        if (!tween.started) {
            tween.started = t;
        }

        switch (tween.type) {
            case 'cam':
                for (var i = 0; i < 4; i++) {
                    for (var j = 0; j < 4; j++) {
                        tween.cam.tm[i][j] = tween.func(t - tween.started, tween.begin[i][j], tween.delta[i][j], tween.dur);
                    }
                }

                tween.cam.tm = normalizeMatrix(tween.cam.tm);
                tween.cam.isValid(false);
                break;
            case 'view':
                tween.cam.setView({
                    start: {
                        a: tween.func(t - tween.started, tween.begin.start.a, tween.delta.start.a, tween.dur),
                        r: tween.func(t - tween.started, tween.begin.start.r, tween.delta.start.r, tween.dur)
                    },
                    center: {
                        a: tween.func(t - tween.started, tween.begin.center.a, tween.delta.center.a, tween.dur),
                        r: tween.func(t - tween.started, tween.begin.center.r, tween.delta.center.r, tween.dur),
                        z: tween.func(t - tween.started, tween.begin.center.z, tween.delta.center.z, tween.dur)
                    },
                    end: {
                        a: tween.func(t - tween.started, tween.begin.end.a, tween.delta.end.a, tween.dur),
                        r: tween.func(t - tween.started, tween.begin.end.r, tween.delta.end.r, tween.dur)
                    }
                });

                tween.cam.isValid(false);
                break;
            case 'obj':
                tween.obj[tween.prop] = tween.func(t - tween.started, tween.begin, tween.delta, tween.dur);
                break;
        }

        tween.finished = (t >= (tween.started + tween.dur));
    },
    easeIn: function (t, b, c, d) {
        t = Math.min(t, d);
        return c * (t /= d) * t + b;
    },
    easeOut: function (t, b, c, d) {
        t = Math.min(t, d);
        return -c * (t /= d) * (t - 2) + b;
    },
    easeInOut: function (t, b, c, d) {
        t = Math.min(t, d);
        if ((t /= d * .5) < 1) {
            return c / 2 * t * t + b;
        }
        return -c * .5 * ((--t) * (t - 2) - 1) + b;
    }
}.init();
*/

/*
function levelRadius(level) {
    var base = 1.1;
    var mul = 833.349487788449;
    var add = -1766.35863684255;
    var pow = 8;
    return Math.pow(base, level + pow) * mul + add;
}
*/
/*
function levelHeight(level) {
    var base = 1.1;
    var mul = 148.756660081479;
    var add = -100.760423159995;
    var pow = 8;
    return -(Math.pow(base, level + pow) * mul + add);
}
*/

function round(number, digits) {
    digits = isNaN(digits) ? 2 : digits;
    var e = Math.pow(10, Math.round(digits));
    return Math.round(number * e) / e;
}

/*
function galaxyGuides() {
    var arr = [];

    var step = 32;

    for (var i = 1; i <= 12; i++) {
        if (i > 0) {
            var r = levelRadius(i);
            var da = .1 / i;
            var a = 0;
            do {
                arr.push({pos: {x: r * Math.cos(a), y: r * Math.sin(a), z: 0}});
                a += da;
            } while (a < 2 * Math.PI);
        }
        a = Math.PI / 12 * i;
        for (var j = -5; j > -120; --j) {
            arr.push({pos: {x: step * j * Math.cos(a), y: step * j * Math.sin(a), z: 0}});
            arr.push({pos: {x: step * -j * Math.cos(a), y: step * -j * Math.sin(a), z: 0}});
        }
    }

    return arr;
}
*/

// function compareObjZ(a, b) {
// 	return a.pos.z - b.pos.z;
// }
function drawObject(context, obj) {
    var src = starTypes[obj.cls];

    var dw = src.canvas.width / src.scale * obj.pos.xScale * .3;
	if (dw < 3) {
		return;
	}
	var dh = src.canvas.width / src.scale * obj.pos.yScale * .3;
	context.globalAlpha = 1 / Math.max(1, obj.pos.xScale * .2);
/*
    if (obj.cls === 0) {
        dw = Math.min(4, dw);
        dh = Math.min(4, dh);
    } else if (obj.cls !== 6) {
        dw /= Math.max(1, 2 * Math.log(obj.pos.xScale));
        // noinspection JSSuspiciousNameCombination
		dh /= Math.max(1, 2 * Math.log(obj.pos.yScale));
    }
*/

    context.drawImage(src.canvas, obj.pos.x - dw * .5, obj.pos.y - dh * .5, dw, dh);
}

var objectProjections = [];
/*
var sectorProjections = [];
*/

function concatSubArr(c, arr) {
    if (!c) {
        return [];
    }

    return arr[arr.length - c].concat(concatSubArr(c - 1, arr));
}

var ts = Date.now() * 1E-3;

function calcSystemProjections(cam) {
    if (currentView !== 'system') {
        return;
	}

    cam.isValid(false);
    /*
        orbital speed = SQRT(G * M / R), where G = 6.67192E-11 (gravitational constant), M - center mass, R - orbit radius.
        Mass in kg, radius in meters, speed in meters/second
        angular speed in degree/second = 360 / (2 * R * PI / V), where R - orbit radius in meters, PI ~3.14, V - orbital speed.
        all in one formula: 360 / (2 * PI * SQRT(G * M * R^3) / (G * M)) - degrees per second.
    */

    objectProjections = [];
    var o_ae = 149599300000;
    var o_mult = 0.22352941;
    var o_base = 1.7;
    var cm_mult = 1055203174765460;
    var cm_base = 1.42152718630078;
    system.object.chars.size = 102; // stub;
    var cm = cm_mult * Math.exp(system.object.chars.size * Math.log(cm_base));
    var G = 6.67192E-11;
    var Gm = 1E-9;
    var s_Gm = 1E-7;
    var i, j;

    var s_mm = 1E6;
    var s_km = .05;
    var s_mult = 63.2455;
    var s_base = 1.32879;

    /* Star */
    var cAPos = {
        x: 0,
        y: 0,
        z: 0
    };
    var cPos = cam.projection(cAPos);
    if (cPos.z > 0 && cPos.inView) {
        objectProjections.push({pos: cPos, cls: 6});
    }

    var tm = (ts - Date.now() * 1E-3) * 86400 * 10;

    for (i = 0; i < system.children.length; i++) {
        // noinspection JSUnresolvedVariable
        var orbit = system.children[i].coordinates.distance;
        var r = o_ae * o_mult * Math.exp(orbit * Math.log(o_base));
        var av = 360 / (2 * Math.PI * Math.sqrt(G * cm * r * r * r) / (G * cm));
        // noinspection JSUnresolvedVariable
        var a = ((system.children[i].coordinates.angle + tm * av) % 360) * Math.PI / 180;
        cAPos = {
            x: Gm * r * Math.cos(a),
            y: Gm * r * Math.sin(a),
            z: 0
        };
        system.children[i].orbit = r;
        cPos = cam.projection(cAPos);
        if (cPos.z > 0 && cPos.inView) {
            objectProjections.push({pos: cPos, cls: 5});
        }
        if (typeof(system.children[i].children) !== 'undefined') {
            var s_cm = cm_mult * Math.exp(system.children[i].chars['size'] * Math.log(cm_base));
            var s_ok = system.children[i].chars['size'] * s_km;
            for (j = 0; j < system.children[i].children.length; j++) {
                // noinspection JSUnresolvedVariable
                var s_orbit = system.children[i].children[j].coordinates.distance;
                var s_r = s_ok * s_mm * s_mult * Math.exp(s_orbit * Math.log(s_base));
                var s_av = 360 / (2 * Math.PI * Math.sqrt(G * s_cm * s_r * s_r * s_r) / (G * s_cm));
                // noinspection JSUnresolvedVariable
                var s_a = ((system.children[i].children[j].coordinates.angle + tm * s_av) % 360) * Math.PI / 180;
                var s_cAPos = {
                    x: cAPos.x + s_Gm * s_r * Math.cos(s_a),
                    y: cAPos.y + s_Gm * s_r * Math.sin(s_a),
                    z: 0
                };
                // noinspection JSUndefinedPropertyAssignment
                system.children[i].children[j].orbit = {c: cAPos, r: s_r};
                cPos = cam.projection(s_cAPos);
                if (cPos.z > 0 && cPos.inView) {
                    objectProjections.push({pos: cPos, cls: 2});
                }
            }
        }
    }

    setTimeout(delegate(this, arguments.callee, cam), 100);
}

var Hyperbola = (function () {
	function Hyperbola(a, b, O, beta) {
		this.a = a;
		this.b = b;
		this.O = O;
		this.beta = beta;
		this.alpha = Math.abs(Math.atan2(this.b, this.a) - Math.atan2(this.b, -this.a));
		this.c = Math.sqrt(a * a + b * b);
	}

	Hyperbola.prototype.x = function (y, s) {
		return sign(s) * this.a / this.b * Math.sqrt(this.b * this.b + y * y);
	};

	Hyperbola.prototype.y = function (x, s) {
		return sign(s) * this.b / this.a * Math.sqrt(x * x - this.a * this.a);
	};

	return Hyperbola;
})();

function sign(s) {
	if (s < 0)
		return -1;

	return 1;
}

var Mat2d = {
	det: function(A) {
		return A[0]*A[4]*A[8] + A[1]*A[5]*A[6] + A[2]*A[3]*A[7] - A[2]*A[4]*A[6] - A[1]*A[3]*A[8] - A[0]*A[5]*A[7];
	},

	invMatrix: function (A) {
		var d = this.det(A);
		return [
			(A[4]*A[8]-A[5]*A[7])/d,
			-(A[1]*A[8]-A[2]*A[7])/d,
			(A[1]*A[5]-A[2]*A[4])/d,
			-(A[3]*A[8]-A[5]*A[6])/d,
			(A[0]*A[8]-A[2]*A[6])/d,
			-(A[0]*A[5]-A[2]*A[3])/d,
			(A[3]*A[7]-A[4]*A[6])/d,
			-(A[0]*A[7]-A[1]*A[6])/d,
			(A[0]*A[4]-A[1]*A[3])/d
		];
	},

	mulVector: function (v, m) {
		return [
			v[0] * m[0] + v[1] * m[1] + v[2] * m[2],
			v[0] * m[3] + v[1] * m[4] + v[2] * m[5],
			v[0] * m[6] + v[1] * m[7] + v[2] * m[8]
		];
	}
};

function calcProjections(cam) {
    if (currentView === 'system') {
        calcSystemProjections(cam);
        return;
    }
    var i, j;

    var pos;
/*
    var guidesProjections = [];
    for (i = 0; i < guides.length; ++i) {
        pos = cam.projection(guides[i].pos);
        if (pos.z > 0 && pos.inView) {
            guidesProjections.push({pos: pos, cls: 0});
        }
    }
*/

    var starProjections = [];
    for (i = 0; i < starTypes.length; i++) {
        starProjections.push([]);
    }

    if (showBH) {
		var invCamTm = inverseMatrix(cam.tm);
		var camPoint = mulVectorMatrix([0, 0, 0, 1], invCamTm);
		for (i = 0; i < stars.length; i++) {
			var vz = normalizeVector(crossVectorVector([stars[i].pos.x, stars[i].pos.y, stars[i].pos.z, 1], camPoint));
			var vx = normalizeVector(crossVectorVector([stars[i].pos.x, stars[i].pos.y, stars[i].pos.z, 1], vz));
			var vy = normalizeVector(crossVectorVector(vz, vx));
			var vMatrix = [vx, vy, vz, [0, 0, 0, 1]];
			var invVMatrix = normalizeMatrix(inverseMatrix(vMatrix));
			var vA = mulVectorMatrix(camPoint, invVMatrix);
			var vB = mulVectorMatrix([stars[i].pos.x, stars[i].pos.y, stars[i].pos.z, 1], invVMatrix);
			var vC = mulVectorMatrix([0, 0, 0, 1], invVMatrix);
			var A = {
					x: vA[0],
					y: vA[1]
				},
				B = {
					x: vB[0],
					y: vB[1]
				},
				C = {
					x: vC[0],
					y: vC[1]
				},
				r1A = dist2d(A, C),
				r1B = dist2d(B, C),
				hypA = CR,
				r2A = 2 * hypA + r1A,
				r2B = 2 * hypA + r1B,
				O1 = getCrossPoints(A, r2A, B, r2B);
			var h, mat, imat, s;
			for (j = 0; j < 2; ++j) {
				var hypC = dist2d(C, O1[j]) / 2;
				if ((hypC - hypA) < CR) {
					continue;
				}
				var hypB = Math.sqrt(hypC * hypC - hypA * hypA);
				h = new Hyperbola(
					hypA,
					hypB,
					{
						x: (C.x + .5 * (O1[j].x - C.x)),
						y: (C.y + .5 * (O1[j].y - C.y))
					},
					Math.atan2(O1[j].y - C.y, O1[j].x - C.x)
				);
				s = 1;
				mat = [
					Math.cos(h.beta), -Math.sin(h.beta), h.O.x,
					Math.sin(h.beta), Math.cos(h.beta), h.O.y,
					0, 0, 1];
				imat = Mat2d.invMatrix(mat);
				var _vB = [B.x, B.y, 1];
				var _vBm1 = Mat2d.mulVector(_vB, imat);
				if (Math.abs(_vBm1[0] - h.x(_vBm1[1], s)) > E) {
					s = -1;
				}
				_vBm1[0] = h.x(_vBm1[1], -s);
				_vB = Mat2d.mulVector(_vBm1, mat);
				var _vB3d = mulVectorMatrix([_vB[0], _vB[1], 0, 1], vMatrix);

				var vStarPosG = {
					x: _vB3d[0],
					y: _vB3d[1],
					z: _vB3d[2]
				};
				pos = cam.projection(vStarPosG);
				if (pos.z > 0 && pos.inView) {
					starProjections[stars[i].cls].push({star: i, pos: pos, cls: stars[i].cls});
				}
			}
		}
	} else {
		for (i = 0; i < stars.length; i++) {
			pos = cam.projection(stars[i].pos);
			if (pos.z > 0 && pos.inView) {
				starProjections[stars[i].cls].push({star: i, pos: pos, cls: stars[i].cls});
			}
		}
	}

//    objectProjections = guidesProjections.concat(concatSubArr(starProjections.length, starProjections));
    objectProjections = concatSubArr(starProjections.length, starProjections);

/*
    sectorProjections = [];

    for (i = 0; i < sectors.length; i++) {
        var inView = sectors[i].focus.length;
        var points = [];
        for (var j = 0; j < sectors[i].focus.length; j++) {
            sectors[i].focus[j].pos.z = 0;
            pos = cam.projection(sectors[i].focus[j].pos);
            if (pos.z <= 0 || !pos.inView) {
                inView--;
            }
            points.push(pos);
        }
        if (inView) {
            sectors[i].projection = points;
            sectorProjections.push({sector: sectors[i], proj: points});
        }
    }
*/
}

function renderTriangle(ctx, bitmap, vs, sx0, sy0, sx1, sy1, sx2, sy2)
{
	var x0 = vs[0].x, y0 = vs[0].y,
		x1 = vs[1].x, y1 = vs[1].y,
		x2 = vs[2].x, y2 = vs[2].y;
	ctx.beginPath();
	ctx.moveTo(x0, y0);
	ctx.lineTo(x1, y1);
	ctx.lineTo(x2, y2);
	ctx.closePath();
	ctx.clip();

	// Textured triangle transformation code originally by Thatcher Ulrich
	// TODO: figure out if drawImage goes faster if we specify the rectangle that bounds the source coords.
	// TODO: this is far from perfect - due to perspective corrected texture mapping issues see:
	//       http://tulrich.com/geekstuff/canvas/perspective.html

	// collapse terms
	var denom = 1.0 / (sx0 * (sy2 - sy1) - sx1 * sy2 + sx2 * sy1 + (sx1 - sx2) * sy0);
	// calculate context transformation matrix
	var m11 = - (sy0 * (x2 - x1) - sy1 * x2 + sy2 * x1 + (sy1 - sy2) * x0) * denom,
		m12 = (sy1 * y2 + sy0 * (y1 - y2) - sy2 * y1 + (sy2 - sy1) * y0) * denom,
		m21 = (sx0 * (x2 - x1) - sx1 * x2 + sx2 * x1 + (sx1 - sx2) * x0) * denom,
		m22 = - (sx1 * y2 + sx0 * (y1 - y2) - sx2 * y1 + (sx2 - sx1) * y0) * denom,
		dx = (sx0 * (sy2 * x1 - sy1 * x2) + sy0 * (sx1 * x2 - sx2 * x1) + (sx2 * sy1 - sx1 * sy2) * x0) * denom,
		dy = (sx0 * (sy2 * y1 - sy1 * y2) + sy0 * (sx1 * y2 - sx2 * y1) + (sx2 * sy1 - sx1 * sy2) * y0) * denom;

	ctx.transform(m11, m12, m21, m22, dx, dy);

	// Draw the whole texture image. Transform and clip will map it onto the correct output polygon.
	ctx.drawImage(bitmap, 0, 0);
}

function calcMidPosFraction(a, b, i, f) {
	return {
		x: a.x + (b.x - a.x) * i / f,
		y: a.y + (b.y - a.y) * i / f,
		z: a.z + (b.z - a.z) * i / f
	};
}

function inflatePolygon (vertices, pixels) {
	pixels = pixels || 0.5;
	var inflatedVertices = new Array(vertices.length);
	for (var i1 = 0; i1 < vertices.length; i1++) {
		inflatedVertices[i1] = {x: vertices[i1].x, y: vertices[i1].y};
	}

	for (var i = 0, j = vertices.length, k, x1, y1, x2, y2; i < j; i++) {
		k = (i < j - 1) ? (i + 1) : 0;
		x1 = inflatedVertices[i].x;
		y1 = inflatedVertices[i].y;
		x2 = inflatedVertices[k].x;
		y2 = inflatedVertices[k].y;
		var x = x2 - x1, y = y2 - y1,
			det = x * x + y * y, idet;

		if (det === 0) det = 1E-6;

		idet = pixels / Math.sqrt(det);

		x *= idet;
		y *= idet;

		inflatedVertices[i].x -= x;
		inflatedVertices[i].y -= y;
		inflatedVertices[k].x += x;
		inflatedVertices[k].y += y;
	}
	return inflatedVertices;
}

function renderRectangle(texture, nw, ne, se, sw, fractions, ctx, cam) {
	fractions = fractions || 10;
	var shift = .1;
	for (var i = 0; i < fractions; i++) {
		var twi = calcMidPosFraction(nw, ne, i, fractions),
			tei = calcMidPosFraction(nw, ne, i + 1, fractions),
			bwi = calcMidPosFraction(sw, se, i, fractions),
			bei = calcMidPosFraction(sw, se, i + 1, fractions);
		for (var j = 0; j < fractions; j++) {
			var nwj = calcMidPosFraction(twi, bwi, j, fractions),
				nej = calcMidPosFraction(tei, bei, j, fractions),
				swj = calcMidPosFraction(twi, bwi, j + 1, fractions),
				sej = calcMidPosFraction(tei, bei, j + 1, fractions),
				_nw = cam.projection(nwj),
				_ne = cam.projection(nej),
				_sw = cam.projection(swj),
				_se = cam.projection(sej);
			if (_sw.inView && _sw.z > 0 && _nw.inView && _nw.z > 0 && _ne.inView && _ne.z > 0) {
				var tx = texture.naturalWidth / fractions,
					ty = texture.naturalHeight / fractions;
				var vs = [_sw, _nw, _ne],
					tx0 = i * tx,
					ty0 = (j + 1) * ty,
					tx1 = i * tx,
					ty1 = j * ty,
					tx2 = (i + 1) * tx,
					ty2 = j * ty;
				ctx.save();
				ctx.globalAlpha = Math.max(0, Math.min(1, ((_sw.z + _ne.z) - 1000) * 2.5E-4));
				renderTriangle(ctx, texture, inflatePolygon(vs, shift), tx0, ty0, tx1, ty1, tx2, ty2);
				ctx.restore();
			}
			if (_ne.inView && _ne.z > 0 && _se.inView && _se.z > 0 && _sw.inView && _sw.z > 0) {
				vs = [_ne, _se, _sw];
				tx0 = (i + 1) * tx;
				ty0 = j * ty;
				tx1 = (i + 1) * tx;
				ty1 = (j + 1) * ty;
				tx2 = i * tx;
				ty2 = (j + 1) * ty;
				ctx.save();
				ctx.globalAlpha = Math.max(0, Math.min(1, ((_sw.z + _ne.z) - 1000) * 2.5E-4));
				renderTriangle(ctx, texture, inflatePolygon(vs, shift), tx0, ty0, tx1, ty1, tx2, ty2);
				ctx.restore();
			}
		}
	}
}

function draw(ctx, cam) {
    var x0 = cam.limits.x * .5;
    var y0 = cam.limits.y * .5;

	ctx.setTransform(1, 0, 0, 1, x0, y0);
	ctx.globalCompositeOperation = 'screen';

	var bgImg = MyCache.getImage('img/galaxy_bg1.png');
// 	var bgImg = {complete: true, naturalWidth: 600, naturalHeight: 600};
	if (bgImg && bgImg.complete) {
	    var mul = 3.6,
			ww = mul * bgImg.naturalWidth,
			hh = mul * bgImg.naturalHeight;
	    var cen = cam.projection({x: 0, y: 0, z: 0});
		if (cen.z > 0) {
		    renderRectangle(bgImg,
			    {x: -(ww >> 1), y: -(hh >> 1), z: 0},
			    {x: (ww >> 1), y: -(hh >> 1), z: 0},
			    {x: (ww >> 1), y: (hh >> 1), z: 0},
			    {x: -(ww >> 1), y: (hh >> 1), z: 0},
			    20, ctx, cam);
			ctx.globalAlpha = 1;
		}
	}

    for (var i = 0; i < objectProjections.length; i++) {
        drawObject(ctx, objectProjections[i]);
    }
}

/*
function galaxySectors() {
    var sn = [
        'A', 'B', 'C', 'D', 'E', 'F',
        'G', 'H', 'I', 'J', 'K', 'L',
        'M', 'N', 'O', 'P', 'Q', 'R',
        'S', 'T', 'U', 'V', 'W', 'X'
    ];

    var arr = [];

    var da = Math.PI / 72;
    for (var i = 0; i < 12; i++) {
        for (var a = 0; a < 24; a++) {
            var sector = {
                name: sn[a] + (i + 1),
                neighbours: {},
                start: {r: 0, a: 0},
                end: {r: 0, a: 0},
                center: {x: 0, y: 0, z: 0, r: 0, a: 0},
                focus: [],
                mesh: {meridians: [], parallels: []}
            };
            sector.neighbours.left = (!a ? sn[23] : sn[a - 1]) + (i + 1);
            sector.neighbours.right = (a === 23 ? sn[0] : sn[a + 1]) + (i + 1);
            sector.neighbours.up = i === 11 ? null : sn[a] + (i + 2);
            sector.neighbours.down = !i ? null : sn[a] + (i);

            sector.start.r = levelRadius(i);
            sector.end.r = levelRadius(i + 1);
            sector.start.a = Math.PI * a / 12;
            sector.end.a = Math.PI * (a + 1) / 12;

            var pCnt = 5;
            for (var j = 0; j <= pCnt; j++) {
                var _a = sector.start.a;
                var parallel = [];
                while (_a <= sector.end.a) {
                    if (!j) {
                        sector.focus.push({
                            pos: {
                                x: sector.start.r * Math.cos(_a),
                                y: sector.start.r * Math.sin(_a),
                                z: 0
                            }
                        });

                        sector.mesh.meridians.push([
                            {
                                x: sector.start.r * Math.cos(_a),
                                y: sector.start.r * Math.sin(_a),
                                z: 0
                            },
                            {
                                x: sector.end.r * Math.cos(_a),
                                y: sector.end.r * Math.sin(_a),
                                z: 0
                            }
                        ]);
                    }
                    var r = sector.start.r + (sector.end.r - sector.start.r) * j / pCnt;
                    parallel.push({
                        x: r * Math.cos(_a),
                        y: r * Math.sin(_a),
                        z: 0
                    });
                    _a += da;
                }
                sector.mesh.parallels.push(parallel);
            }
            _a = sector.end.a;
            while (_a >= sector.start.a - 0.001) {
                sector.focus.push({
                    pos: {
                        x: sector.end.r * Math.cos(_a),
                        y: sector.end.r * Math.sin(_a),
                        z: 0
                    }
                });
                _a -= da;
            }

            sector.center.a = Math.PI * (a + .5) / 12;
            sector.center.r = levelRadius(i + .5);
            sector.center.x = sector.center.r * Math.cos(sector.center.a);
            sector.center.y = sector.center.r * Math.sin(sector.center.a);
            sector.center.z = levelHeight(i + 1);
            arr.push(sector);
        }
    }

    return arr;
}
*/

/*
function sectorHitTest(ctx, cam) {
    var sector = cam.zoomSector();
    if (sector) {
        return sector;
    }

    var x0 = cam.limits.x * .5;
    var y0 = cam.limits.y * .5;

    var lx = mousePos.x;
    var ly = mousePos.y;

    ctx.setTransform(1, 0, 0, 1, x0, y0);
    ctx.fillStyle = '#999';
    ctx.textAlign = 'left';
    ctx.fillText('MOUSE: ' + lx + ', ' + ly, -x0 + 37, -y0 + 52);

    for (var i = 0; i < sectorProjections.length; i++) {
        ctx.beginPath();
        var t = [];
        for (var j = 0; j < sectorProjections[i].proj.length; j++) {
            if (j === 0) {
                ctx.moveTo(sectorProjections[i].proj[j].x, sectorProjections[i].proj[j].y);
            } else {
                ctx.lineTo(sectorProjections[i].proj[j].x, sectorProjections[i].proj[j].y);
            }
            t.push({x: sectorProjections[i].proj[j].x, y: sectorProjections[i].proj[j].y});
        }
        ctx.closePath();
        if (ctx.isPointInPath(lx, ly)) {
            return sectorProjections[i].sector;
        }
    }

    return null;
}
*/

var Router = (function () {
    function Router(start) {
        this.start = start;
        this.mapData = {
            show: false,
            complete: false,
            map: [],
            iter: 0,
            status: 'idle',
            data: null,
            index: null
        };
        this.routesData = {
            complete: false,
            dist: null,
            parent: null,
            mark: null,
            iter: 0
        };
	    this.jumpsData = {
		    dist: null,
            jumps: null,
		    parent: null,
		    mark: null,
			iter: 0
	    };
        this.stage = '';
        this.progress = 0;
        this.step = 30;
        this.iterTimer = setTimeout(delegate(this, iteration), 10);
    }

    function iteration() {
    	clearTimeout(this.iterTimer);
    	this.iterTimer = 0;
        if (!this.mapData.complete) {
            if (this.mapData.status === 'loaded') {
                this.stage = 'INITIALIZATION';
            } else {
                this.stage = 'LOADING ROUTES MAP';
            }
            this.progress = calcMapIter.call(this);
        } else if (!this.routesData.complete) {
            this.stage = 'ROUTING';
            this.progress = calcRoutesIter.call(this, this.start);
        } else {
            return;
        }
        this.iterTimer = setTimeout(delegate(this, iteration), 10);
    }

	function findMapDataIndex(a) {
    	return a.i === this.t;
	}

	function addTunnel(mData, tunnel) {
		var isource = mData.index[tunnel.from];
		var idest = mData.index[tunnel.to];

		if (mData.map[isource].findIndex(findMapDataIndex, {t: idest}) < 0)
			mData.map[mData.index[tunnel.from]].push({
				i: mData.index[tunnel.to],
				d: tunnel.dist
			});
		if (mData.map[idest].findIndex(findMapDataIndex, {t: isource}) < 0)
			mData.map[mData.index[tunnel.to]].push({
				i: mData.index[tunnel.from],
				d: tunnel.dist
			});
	}

	function calcMapIter() {
        if (this.mapData.status === 'idle') {
            this.mapData.status = 'loading';
            var that = this;
            this.mapData.xhr = jQuery.ajax({
				url: 'resources/tunnels.json',
				dataType: 'json',
				jsonp: false,
				crossDomain: true,
				xhrFields: {
					withCredentials: true
				}
			}).done(
				function (data) {
					that.mapData.data = data;
					that.mapData.status = 'loaded';
				}
			);
        }

        if (this.mapData.status === 'loaded') {
            var i;

            if (this.mapData.iter === 0) {
                this.mapData.map = [];
                this.mapData.index = {};
                for (i = 0; i < stars.length; i++) {
                    this.mapData.map[i] = [];
                    this.mapData.index[stars[i].id] = i;
                    stars[i].i = i;
                }
            }

            for (var _i = 0; _i < this.step; _i++) {
				i = this.mapData.iter++;
				if (i < this.mapData.data.length) {
					addTunnel(this.mapData, this.mapData.data[i]);
				} else {
					this.mapData.data = [];
					this.mapData.iter = this.mapData.data.length;
                    this.mapData.complete = true;
                    return 1;
                }
            }

            return (i / (this.mapData.data.length - 1));
        }

        return 0;
    }

    function cmpJumps(v1, v2, len) {
        if (this.jumpsData.jumps[v1] + 1 < this.jumpsData.jumps[v2]) {
	        return -1;
        }
        else if (this.jumpsData.jumps[v1] + 1 === this.jumpsData.jumps[v2]) {
	        if (this.jumpsData.dist[v1] + len < this.jumpsData.dist[v2])
		        return -1;
        }
        else {
	        return 1;
        }
    }

    function calcRoutesIter() {
        var inf = 1000000000,
            i, j;

        if (this.routesData.iter === 0) {
            this.routesData.dist = [];
            this.routesData.parent = [];
            this.routesData.mark = [];
            this.jumpsData.dist = [];
            this.jumpsData.jumps = [];
            this.jumpsData.parent = [];
            this.jumpsData.mark = [];

            for (i = 0; i < this.mapData.map.length; i++) {
                this.routesData.dist[i] = inf;
	            this.jumpsData.dist[i] = inf;
                this.jumpsData.jumps[i] = inf;
            }

            this.routesData.dist[this.start] = 0;
	        this.jumpsData.dist[this.start] = 0;
	        this.jumpsData.jumps[this.start] = 0;
			this.jumpsData.iter = this.routesData.iter;
        }

        for (var _i = 0; _i < this.step; _i++) {
            i = Math.min(this.routesData.iter++, this.mapData.map.length - 1);
            var v = -1, vj = -1;
            for (j = 0; j < this.mapData.map.length; j++) {
                if (!this.routesData.mark[j] && (v === -1 || this.routesData.dist[j] < this.routesData.dist[v])) {
                    v = j;
                }
	            if (!this.jumpsData.mark[j] && (vj === -1 || this.jumpsData.jumps[j] < this.jumpsData.jumps[vj])) {
		            vj = j;
	            }
            }
            var to, len;
            if (v === -1) {
                this.routesData.iter = this.mapData.map.length;
            } else {
                this.routesData.mark[v] = true;
                for (j = 0; j < this.mapData.map[v].length; j++) {
                    to = this.mapData.map[v][j].i;
                    len = this.mapData.map[v][j].d;
                    if (this.routesData.dist[v] + len < this.routesData.dist[to]) {
                        this.routesData.dist[to] = this.routesData.dist[v] + len;
                        this.routesData.parent[to] = v;
                    }
                }
            }
	        if (vj === -1) {
                this.jumpsData.iter = this.mapData.map.length;
			} else {
		        this.jumpsData.mark[vj] = true;
		        for (j = 0; j < this.mapData.map[vj].length; j++) {
			        to = this.mapData.map[vj][j].i;
			        len = this.mapData.map[vj][j].d;
			        if (cmpJumps.call(this, vj, to, len) < 0) {
				        this.jumpsData.jumps[to] = this.jumpsData.jumps[vj] + 1;
				        this.jumpsData.dist[to] = this.jumpsData.dist[vj] + len;
				        this.jumpsData.parent[to] = vj;
			        }
		        }
	        }
            if ((this.routesData.iter === this.mapData.map.length)
                && (this.jumpsData.iter === this.mapData.map.length)) {
                this.routesData.complete = true;
				jQuery(window).trigger('routerUpdated');
                return 1;
            }
        }

        return (i / (this.mapData.map.length - 1));
    }

    function equalRoutes(a, b) {
    	if (a.path.length !== b.path.length || a.length !== b.length) {
    		return false;
		}

		for (var i = 0; i < a.path.length; ++i) {
    		if (a.path[i] !== b.path[i]) {
    			return false;
			}
		}

		return true;
	}

    Router.prototype.getRoutes = function (end) {
        var routes = [],
            path = [],
            dist = [],
            fail = false,
			v;

		for (v = end.i; v !== this.start; v = this.routesData.parent[v]) {
			if (typeof v === 'undefined') {
				fail = true;
				break;
			}
			path.unshift(v);
			dist.unshift(this.routesData.dist[v]);
		}
		path.unshift(this.start);
		dist.unshift(0);

		if (!fail) {
			routes.push({path: path.slice(), dist: dist.slice(), length: this.routesData.dist[end.i]});
		}

		fail = false;
	    path = [];
		dist = [];

	    for (v = end.i; v !== this.start; v = this.jumpsData.parent[v]) {
		    if (typeof v === 'undefined') {
				fail = true;
				break;
		    }
		    path.unshift(v);
		    dist.unshift(this.jumpsData.dist[v]);
	    }
	    path.unshift(this.start);
	    dist.unshift(0);

		if (!fail) {
			routes.push({path: path, dist: dist, length: this.jumpsData.dist[end.i]});
		}

		var i = 0, j;
		while (i < routes.length - 1) {
			j = i + 1;
			while (j < routes.length) {
				if (equalRoutes(routes[i], routes[j])) {
					routes.splice(j, 1);
				} else {
					++j;
				}
			}
			++i;
		}

        return routes;
    };

    Router.prototype.progressFn = function () {
        return {action: this.stage, progress: this.progress};
    };

    Router.prototype.setLocation = function(star) {
    	this.start = star;
		if (this.mapData.status !== 'loading') {
			var that = this;
			this.mapData.xhr = jQuery.ajax({
				url: 'mock/star_tunnels.json',
				dataType: 'json',
				jsonp: false,
				crossDomain: true,
				xhrFields: {
					withCredentials: true
				}
			}).done(
				function (data) {
					jQuery.each(data, function(i, v) {
						addTunnel(that.mapData, v);
					});
					that.routesData.iter = 0;
					that.routesData.complete = false;

					if (!that.iterTimer) {
						that.iterTimer = setTimeout(delegate(that, iteration), 10);
					}
				}
			);
		}
	};

    return Router;
})();

function starHitTest(ctx, cam) {
    var x0 = cam.limits.x * .5;
    var y0 = cam.limits.y * .5;

    var lx = mousePos.x;
    var ly = mousePos.y;

    ctx.setTransform(1, 0, 0, 1, x0, y0);
    ctx.fillStyle = '#999';
    ctx.textAlign = 'left';
    ctx.fillText('MOUSE: ' + lx + ', ' + ly, -x0 + 37, -y0 + 52);

    for (var i = 0; i < objectProjections.length; i++) {
        var obj = objectProjections[i];
        if (obj.pos.z > 500) {
        	continue;
		}
        if (obj.cls > 0 && obj.cls < 6) {
            var r = 10;
            ctx.beginPath();
            ctx.arc(obj.pos.x, obj.pos.y, r, 0, Math.PI * 2);

            if (ctx.isPointInPath(lx, ly)) {
                return obj;
            }
        }
    }

    return null;
}

var focusArc = {ACW: true, last: 0};

function drawStarFocus(ctx, cam) {
    if (currentView !== 'galaxy') {
        return;
    }

    var starProjection = starHitTest(ctx, cam);
    if (starProjection === null) {
		jQuery(window).trigger({type: 'starHover', star: -1});
        return;
    }
	jQuery(window).trigger({type: 'starHover', star: starProjection.star});

    var x0 = cam.limits.x * .5;
    var y0 = cam.limits.y * .5;

    ctx.setTransform(1, 0, 0, 1, x0, y0);
    ctx.globalCompositeOperation = 'source-over';

    ctx.strokeStyle = '#809dff';
    ctx.lineWidth = 2;
    ctx.beginPath();
    var t = Date.now();
    var k1 = 1000;
    var k2 = k1 * .6;
    var a1 = 2 * Math.PI * ((t % k1) / k1);
    var a2 = 2 * Math.PI * ((t % k2) / k2);
    if (a2 > a1 && (a2 - a1) < .1) {
        if ((t - focusArc.last) > k2) {
            focusArc.last = t;
            focusArc.ACW = !focusArc.ACW;
        }
    }
    ctx.arc(starProjection.pos.x, starProjection.pos.y, 18, a1, a2, focusArc.ACW);
    ctx.stroke();
}

/*
function drawFocus(ctx, cam) {
    if (currentView !== 'galaxy') {
        return;
    }

    var sector = sectorHitTest(ctx, cam);
    if (sector === null) {
        return;
    }

    var x0 = cam.limits.x * .5;
    var y0 = cam.limits.y * .5;

    ctx.setTransform(1, 0, 0, 1, x0, y0);
    ctx.globalCompositeOperation = 'source-over';

    ctx.strokeStyle = 'rgba(0, 255, 0, .5)';
    ctx.lineWidth = 1;
    ctx.fillStyle = 'rgba(0, 255, 0, .5)';

    ctx.beginPath();
    for (var i = 0; i < sector.projection.length; i++) {
        var pos = sector.projection[i];
        if (i > 0) {
            ctx.lineTo(pos.x, pos.y);
        } else {
            ctx.moveTo(pos.x, pos.y);
        }
    }
    ctx.closePath();
    ctx.fill();
    ctx.stroke();
}
*/
/*
function drawSectorMesh(ctx, cam) {
    if (currentView !== 'sector' || (!cam.sectorMesh.length)) {
        return;
    }

    var x0 = cam.limits.x * .5;
    var y0 = cam.limits.y * .5;

    ctx.setTransform(1, 0, 0, 1, x0, y0);
    ctx.globalCompositeOperation = 'source-over';

    ctx.strokeStyle = 'rgba(0, 255, 0, .5)';
    ctx.lineWidth = 1;
    ctx.fillStyle = 'rgba(0, 255, 0, .5)';

    var pos, i, j, k;

    ctx.beginPath();
    for (k = 0; k < cam.sectorMesh.length; k++) {
        for (i = 0; i < cam.sectorMesh[k].meridians.length; i++) {
            pos = cam.projection(cam.sectorMesh[k].meridians[i][0]);
            ctx.moveTo(pos.x, pos.y);
            pos = cam.projection(cam.sectorMesh[k].meridians[i][1]);
            ctx.lineTo(pos.x, pos.y);
        }
        for (i = 0; i < cam.sectorMesh[k].parallels.length; i++) {
            for (j = 0; j < cam.sectorMesh[k].parallels[i].length; j++) {
                pos = cam.projection(cam.sectorMesh[k].parallels[i][j]);
                if (!j) {
                    ctx.moveTo(pos.x, pos.y);
                } else {
                    ctx.lineTo(pos.x, pos.y);
                }
            }
        }
    }
    ctx.stroke();
}
*/
function drawOrbits(ctx, cam) {
	if (currentView !== 'system') {
		return;
	}

    var Gm = 1E-9;
    var s_Gm = 1E-7;
    var i, j;

    var x0 = cam.limits.x * .5;
    var y0 = cam.limits.y * .5;

    ctx.setTransform(1, 0, 0, 1, x0, y0);
    ctx.globalCompositeOperation = 'source-over';
    ctx.lineWidth = 1;
    ctx.strokeStyle = '#666666';
    ctx.beginPath();
    for (i = 0; i < system.children.length; i++) {
        var clPos = null;
        var csPos = null;
        if (typeof(system.children[i].orbit) !== 'undefined') {
            for (var ad = 0; ad < 360; ad++) {
                var ar = ad * Math.PI / 180;
                var cAPos = {
                    x: Gm * system.children[i].orbit * Math.cos(ar),
                    y: Gm * system.children[i].orbit * Math.sin(ar),
                    z: 0
                };
                var cPos = cam.projection(cAPos);
                if (!(cPos.z > 0 && cPos.inView)) {
                    clPos = null;
                    cPos = null;
                } else {
                    if (ad === 0) {
                        csPos = cPos;
                    }
                    if (clPos !== null) {
                        ctx.moveTo(clPos.x, clPos.y);
                        ctx.lineTo(cPos.x, cPos.y);
                    }
                    clPos = cPos;
                }
            }
            if (csPos !== null && cPos !== null) {
                ctx.moveTo(cPos.x, cPos.y);
                ctx.lineTo(csPos.x, csPos.y);
            }
        }
        if (typeof(system.children[i].children) !== 'undefined') {
            for (j = 0; j < system.children[i].children.length; j++) {
                var s_clPos = null;
                var s_csPos = null;
                if (typeof(system.children[i].children[j].orbit) !== 'undefined') {
                    for (var s_ad = 0; s_ad < 360; s_ad++) {
                        var s_ar = s_ad * Math.PI / 180;
                        var s_cAPos = {
                            x: system.children[i].children[j].orbit.c.x + s_Gm * system.children[i].children[j].orbit.r * Math.cos(s_ar),
                            y: system.children[i].children[j].orbit.c.y + s_Gm * system.children[i].children[j].orbit.r * Math.sin(s_ar),
                            z: 0
                        };
                        var s_cPos = cam.projection(s_cAPos);
                        if (!(s_cPos.z > 0 && s_cPos.inView)) {
                            s_clPos = null;
                            s_cPos = null;
                        } else {
                            if (s_ad === 0) {
                                s_csPos = s_cPos;
                            }
                            if (s_clPos !== null) {
                                ctx.moveTo(s_clPos.x, s_clPos.y);
                                ctx.lineTo(s_cPos.x, s_cPos.y);
                            }
                            s_clPos = s_cPos;
                        }
                    }
                    if (s_csPos !== null && s_cPos !== null) {
                        ctx.moveTo(s_cPos.x, s_cPos.y);
                        ctx.lineTo(s_csPos.x, s_csPos.y);
                    }
                }
            }
        }
    }
    ctx.stroke();
}

function drawTunnels(ctx, cam) {
	if (currentView !== 'galaxy') {
		return;
	}

	if (!(router && router.mapData.show && router.mapData.complete)) {
		return;
	}

	var x0 = cam.limits.x * .5;
	var y0 = cam.limits.y * .5;

	ctx.setTransform(1, 0, 0, 1, x0, y0);
	ctx.globalCompositeOperation = 'screen';
	ctx.globalAlpha = 0.5;
	ctx.lineWidth = 1;
	ctx.strokeStyle = '#999999';

	ctx.beginPath();
	for (var i = 0; i < router.mapData.map.length; i++) {
	    for (var j = 0; j < router.mapData.map[i].length; j++) {
	        if (router.mapData.map[i][j].i < i) {
	            continue;
			}

			var pos1 = cam.projection(stars[i].pos),
	            pos2 = cam.projection(stars[router.mapData.map[i][j].i].pos);
	        if ((pos1.z > 0 && pos1.inView) || (pos2.z > 0 && pos2.inView)) {
	            ctx.moveTo(pos1.x, pos1.y);
	            ctx.lineTo(pos2.x, pos2.y);
			}
		}
	}
	ctx.stroke();
}

function drawRoutes(ctx, cam) {
    if (currentView === 'system') {
        return;
    }
    if (routes === null) {
        return;
    }

    var x0 = cam.limits.x * .5;
    var y0 = cam.limits.y * .5;

    ctx.setTransform(1, 0, 0, 1, x0, y0);
    ctx.globalCompositeOperation = 'source-over';

    ctx.lineWidth = 2;

    var target = pLocation.star;
    /* todo: optimize next nested loops */
    for (var i = 0; i < routes.length; i++) {
	    ctx.strokeStyle = routeStrokeStyles[i % 6];
        ctx.beginPath();
        for (var j = 0; j < routes[i].path.length - 1; j++) {
            var nodes = [];
            for (var k = 0; k < objectProjections.length; k++) {
                target = routes[i].path[j + 1];
                if (objectProjections[k].star === routes[i].path[j]) {
                    nodes.push(k);
                    if (nodes.length > 1) {
                        break;
                    }
                }
                if (objectProjections[k].star === routes[i].path[j + 1]) {
                    nodes.push(k);
                    if (nodes.length > 1) {
                        break;
                    }
                }
            }
            if (nodes.length > 1) {
                ctx.moveTo(objectProjections[nodes[0]].pos.x, objectProjections[nodes[0]].pos.y);
                ctx.lineTo(objectProjections[nodes[1]].pos.x, objectProjections[nodes[1]].pos.y);
            }
        }
        ctx.stroke();
    }
}

function mulMatrix(a, b) {
    return [
        [
            a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0] + a[0][3] * b[3][0],
            a[0][0] * b[0][1] + a[0][1] * b[1][1] + a[0][2] * b[2][1] + a[0][3] * b[3][1],
            a[0][0] * b[0][2] + a[0][1] * b[1][2] + a[0][2] * b[2][2] + a[0][3] * b[3][2],
            a[0][0] * b[0][3] + a[0][1] * b[1][3] + a[0][2] * b[2][3] + a[0][3] * b[3][3]
        ],
        [
            a[1][0] * b[0][0] + a[1][1] * b[1][0] + a[1][2] * b[2][0] + a[1][3] * b[3][0],
            a[1][0] * b[0][1] + a[1][1] * b[1][1] + a[1][2] * b[2][1] + a[1][3] * b[3][1],
            a[1][0] * b[0][2] + a[1][1] * b[1][2] + a[1][2] * b[2][2] + a[1][3] * b[3][2],
            a[1][0] * b[0][3] + a[1][1] * b[1][3] + a[1][2] * b[2][3] + a[1][3] * b[3][3]
        ],
        [
            a[2][0] * b[0][0] + a[2][1] * b[1][0] + a[2][2] * b[2][0] + a[2][3] * b[3][0],
            a[2][0] * b[0][1] + a[2][1] * b[1][1] + a[2][2] * b[2][1] + a[2][3] * b[3][1],
            a[2][0] * b[0][2] + a[2][1] * b[1][2] + a[2][2] * b[2][2] + a[2][3] * b[3][2],
            a[2][0] * b[0][3] + a[2][1] * b[1][3] + a[2][2] * b[2][3] + a[2][3] * b[3][3]
        ],
        [
            a[3][0] * b[0][0] + a[3][1] * b[1][0] + a[3][2] * b[2][0] + a[3][3] * b[3][0],
            a[3][0] * b[0][1] + a[3][1] * b[1][1] + a[3][2] * b[2][1] + a[3][3] * b[3][1],
            a[3][0] * b[0][2] + a[3][1] * b[1][2] + a[3][2] * b[2][2] + a[3][3] * b[3][2],
            a[3][0] * b[0][3] + a[3][1] * b[1][3] + a[3][2] * b[2][3] + a[3][3] * b[3][3]
        ]
    ];
}

/*
function copyMatrix(m) {
    return [
        [m[0][0], m[0][1], m[0][2], m[0][3]],
        [m[1][0], m[1][1], m[1][2], m[1][3]],
        [m[2][0], m[2][1], m[2][2], m[2][3]],
        [m[3][0], m[3][1], m[3][2], m[3][3]]
    ];
}
*/
function mulVectorMatrix(v, m) {
    return [
        v[0] * m[0][0] + v[1] * m[1][0] + v[2] * m[2][0] + v[3] * m[3][0],
        v[0] * m[0][1] + v[1] * m[1][1] + v[2] * m[2][1] + v[3] * m[3][1],
        v[0] * m[0][2] + v[1] * m[1][2] + v[2] * m[2][2] + v[3] * m[3][2],
        v[0] * m[0][3] + v[1] * m[1][3] + v[2] * m[2][3] + v[3] * m[3][3]
    ];
}

function normalizeMatrix(m) {
    var m1 = [
        [0, 0, 0, m[0][3]],
        [0, 0, 0, m[1][3]],
        [0, 0, 0, m[2][3]],
        [m[3][0], m[3][1], m[3][2], m[3][3]]
    ];

    for (var i = 0; i < 3; i++) {
        var l = Math.sqrt(m[i][0] * m[i][0] + m[i][1] * m[i][1] + m[i][2] * m[i][2]);
        m1[i][0] = m[i][0] / l;
        m1[i][1] = m[i][1] / l;
        m1[i][2] = m[i][2] / l;
    }

    return m1;
}

function normalizeVector(v) {
	var v1 = [0, 0, 0, v[3]],
		l = Math.sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);

	v1[0] = v[0] / l;
	v1[1] = v[1] / l;
	v1[2] = v[2] / l;

	return v1;
}

function inverseMatrix(a) {
	var d = a[0][0] * a[1][1] * a[2][2] * a[3][3]  +  a[0][0] * a[1][2] * a[2][3] * a[3][1]  +  a[0][0] * a[1][3] * a[2][1] * a[3][2]
		+ a[0][1] * a[1][0] * a[2][3] * a[3][2]  +  a[0][1] * a[1][2] * a[2][0] * a[3][3]  +  a[0][1] * a[1][3] * a[2][2] * a[3][0]
		+ a[0][2] * a[1][0] * a[2][1] * a[3][3]  +  a[0][2] * a[1][1] * a[2][3] * a[3][0]  +  a[0][2] * a[1][3] * a[2][0] * a[3][1]
		+ a[0][3] * a[1][0] * a[2][2] * a[3][1]  +  a[0][3] * a[1][1] * a[2][0] * a[3][2]  +  a[0][3] * a[1][2] * a[2][1] * a[3][0]
		- a[0][0] * a[1][1] * a[2][3] * a[3][2]  -  a[0][0] * a[1][2] * a[2][1] * a[3][3]  -  a[0][0] * a[1][3] * a[2][2] * a[3][1]
		- a[0][1] * a[1][0] * a[2][2] * a[3][3]  -  a[0][1] * a[1][2] * a[2][3] * a[3][0]  -  a[0][1] * a[1][3] * a[2][0] * a[3][2]
		- a[0][2] * a[1][0] * a[2][3] * a[3][1]  -  a[0][2] * a[1][1] * a[2][0] * a[3][3]  -  a[0][2] * a[1][3] * a[2][1] * a[3][0]
		- a[0][3] * a[1][0] * a[2][1] * a[3][2]  -  a[0][3] * a[1][1] * a[2][2] * a[3][0]  -  a[0][3] * a[1][2] * a[2][0] * a[3][1];
	if (!d) {
		console.log('Determinant == 0!');
	}

	return [
		[
			(a[1][1] * a[2][2] * a[3][3]  +  a[1][2] * a[2][3] * a[3][1]  +  a[1][3] * a[2][1] * a[3][2]  -  a[1][1] * a[2][3] * a[3][2]  -  a[1][2] * a[2][1] * a[3][3]  -  a[1][3] * a[2][2] * a[3][1]) / d,
			(a[0][1] * a[2][3] * a[3][2]  +  a[0][2] * a[2][1] * a[3][3]  +  a[0][3] * a[2][2] * a[3][1]  -  a[0][1] * a[2][2] * a[3][3]  -  a[0][2] * a[2][3] * a[3][1]  -  a[0][3] * a[2][1] * a[3][2]) / d,
			(a[0][1] * a[1][2] * a[3][3]  +  a[0][2] * a[1][3] * a[3][1]  +  a[0][3] * a[1][1] * a[3][2]  -  a[0][1] * a[1][3] * a[3][2]  -  a[0][2] * a[1][1] * a[3][3]  -  a[0][3] * a[1][2] * a[3][1]) / d,
			(a[0][1] * a[1][3] * a[2][2]  +  a[0][2] * a[1][1] * a[2][3]  +  a[0][3] * a[1][2] * a[2][1]  -  a[0][1] * a[1][2] * a[2][3]  -  a[0][2] * a[1][3] * a[2][1]  -  a[0][3] * a[1][1] * a[2][2]) / d
		],
		[
			(a[1][0] * a[2][3] * a[3][2]  +  a[1][2] * a[2][0] * a[3][3]  +  a[1][3] * a[2][2] * a[3][0]  -  a[1][0] * a[2][2] * a[3][3]  -  a[1][2] * a[2][3] * a[3][0]  -  a[1][3] * a[2][0] * a[3][2]) / d,
			(a[0][0] * a[2][2] * a[3][3]  +  a[0][2] * a[2][3] * a[3][0]  +  a[0][3] * a[2][0] * a[3][2]  -  a[0][0] * a[2][3] * a[3][2]  -  a[0][2] * a[2][0] * a[3][3]  -  a[0][3] * a[2][2] * a[3][0]) / d,
			(a[0][0] * a[1][3] * a[3][2]  +  a[0][2] * a[1][0] * a[3][3]  +  a[0][3] * a[1][2] * a[3][0]  -  a[0][0] * a[1][2] * a[3][3]  -  a[0][2] * a[1][3] * a[3][0]  -  a[0][3] * a[1][0] * a[3][2]) / d,
			(a[0][0] * a[1][2] * a[2][3]  +  a[0][2] * a[1][3] * a[2][0]  +  a[0][3] * a[1][0] * a[2][2]  -  a[0][0] * a[1][3] * a[2][2]  -  a[0][2] * a[1][0] * a[2][3]  -  a[0][3] * a[1][2] * a[2][0]) / d
		],
		[
			(a[1][0] * a[2][1] * a[3][3]  +  a[1][1] * a[2][3] * a[3][0]  +  a[1][3] * a[2][0] * a[3][1]  -  a[1][0] * a[2][3] * a[3][1]  -  a[1][1] * a[2][0] * a[3][3]  -  a[1][3] * a[2][1] * a[3][0]) / d,
			(a[0][0] * a[2][3] * a[3][1]  +  a[0][1] * a[2][0] * a[3][3]  +  a[0][3] * a[2][1] * a[3][0]  -  a[0][0] * a[2][1] * a[3][3]  -  a[0][1] * a[2][3] * a[3][0]  -  a[0][3] * a[2][0] * a[3][1]) / d,
			(a[0][0] * a[1][1] * a[3][3]  +  a[0][1] * a[1][3] * a[3][0]  +  a[0][3] * a[1][0] * a[3][1]  -  a[0][0] * a[1][3] * a[3][1]  -  a[0][1] * a[1][0] * a[3][3]  -  a[0][3] * a[1][1] * a[3][0]) / d,
			(a[0][0] * a[1][3] * a[2][1]  +  a[0][1] * a[1][0] * a[2][3]  +  a[0][3] * a[1][1] * a[2][0]  -  a[0][0] * a[1][1] * a[2][3]  -  a[0][1] * a[1][3] * a[2][0]  -  a[0][3] * a[1][0] * a[2][1]) / d
		],
		[
			(a[1][0] * a[2][2] * a[3][1]  +  a[1][1] * a[2][0] * a[3][2]  +  a[1][2] * a[2][1] * a[3][0]  -  a[1][0] * a[2][1] * a[3][2]  -  a[1][1] * a[2][2] * a[3][0]  -  a[1][2] * a[2][0] * a[3][1]) / d,
			(a[0][0] * a[2][1] * a[3][2]  +  a[0][1] * a[2][2] * a[3][0]  +  a[0][2] * a[2][0] * a[3][1]  -  a[0][0] * a[2][2] * a[3][1]  -  a[0][1] * a[2][0] * a[3][2]  -  a[0][2] * a[2][1] * a[3][0]) / d,
			(a[0][0] * a[1][2] * a[3][1]  +  a[0][1] * a[1][0] * a[3][2]  +  a[0][2] * a[1][1] * a[3][0]  -  a[0][0] * a[1][1] * a[3][2]  -  a[0][1] * a[1][2] * a[3][0]  -  a[0][2] * a[1][0] * a[3][1]) / d,
			(a[0][0] * a[1][1] * a[2][2]  +  a[0][1] * a[1][2] * a[2][0]  +  a[0][2] * a[1][0] * a[2][1]  -  a[0][0] * a[1][2] * a[2][1]  -  a[0][1] * a[1][0] * a[2][2]  -  a[0][2] * a[1][1] * a[2][0]) / d
		]
	];
}

function crossVectorVector(v1, v2) {
	//a × b = {aybz - azby; azbx - axbz; axby - aybx}
	return [
		v1[1]*v2[2] - v1[2]*v2[1],
		v1[2]*v2[0] - v1[0]*v2[2],
		v1[0]*v2[1] - v1[1]*v2[0],
		0
	]
}

waitImagesLoading();

var cProps = (function () {
    var o = {};
    var p = 576;
    o.hw = p / 9 * 16;
    o.hh = p;
    o.hl = Math.floor((document.body.clientWidth - o.hw) * .5) + 'px';
    o.ht = Math.floor((document.body.clientHeight - o.hh) * .5) + 'px';
    o.w = o.hw - 262;
    o.h = o.hh - 58;
    o.dl = 236;
    o.dt = 36;
    o.l = (parseInt(o.hl) + o.dl) + 'px';
    o.t = (parseInt(o.ht) + o.dt) + 'px';

    return o;
})();

/*
function sectorProjection(pos) {
    var pos3d = {x: pos.x, y: pos.y, z: pos.z * this.modulus};

    var p = polar(pos3d);

    var ad = Math.abs((p.a - this.view.center.a + Math.PI * 2) % (Math.PI * 2));
    if (ad > 1.6 && ad < 4.7) {
        return {
            x: cProps.w << 1,
            y: cProps.h << 1,
            z: -1,
            xScale: 0,
            yScale: 0,
            inView: false
        };
    }

    var pos2d = this._old_projection(pos3d);

    var sx = (-1 + (this.limits.x - 30) / (Math.sin(this.view.end.a - this.view.start.a) * this.view.center.r * pos2d.xScale)) * (1 - this.modulus);
    var sy = (-1 + (this.limits.y - 30) / ((this.view.end.r - this.view.start.r) * pos2d.yScale)) * (1 - this.modulus);

    var x2d = (1 + sx) * (pos2d.x - (pos2d.x - Math.sin(-this.view.center.a + p.a) * this.view.center.r * pos2d.xScale) * (1 - this.modulus));
    var y2d = (1 + sy) * (pos2d.y - (pos2d.y - (this.view.center.r - p.r) * pos2d.yScale) * (1 - this.modulus));

    return {
        x: x2d,
        y: y2d,
        z: pos2d.z,
        xScale: pos2d.xScale,
        yScale: pos2d.yScale,
        inView: (Math.abs(x2d) < (this.limits.x * .999)) && (Math.abs(y2d) < (this.limits.y * .999))
    };
}
*/
/*
function sectorSetView(view) {
    this.view.start.a = view.start.a;
    this.view.start.r = view.start.r;
    this.view.center.a = view.center.a;
    this.view.center.r = view.center.r;
    this.view.center.z = view.center.z;
    this.view.end.a = view.end.a;
    this.view.end.r = view.end.r;

    var scA = -view.center.a - Math.PI * .5;
    var sinA = Math.sin(scA);
    var cosA = Math.cos(scA);

    this.tm = mulMatrix(
        [
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [-Math.cos(view.center.a) * view.center.r, -Math.sin(view.center.a) * view.center.r, -view.center.z, 1]
        ],
        [
            [cosA, sinA, 0, 0],
            [-sinA, cosA, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
        ]
    );
}
*/
/*
function copyView(view) {
    return {
        start: {
            a: view.start.a,
            r: view.start.r
        },
        center: {
            a: view.center.a,
            r: view.center.r,
            z: view.center.z
        },
        end: {
            a: view.end.a,
            r: view.end.r
        }
    };
}
*/
/*
function pushSectorMesh(sector) {
    this.sectorMesh.push(sector.mesh);
}
*/
/*
function shiftSectorMesh() {
    this.sectorMesh.shift();
}
*/
/*
function setSectorName(name) {
    this.sectorName = name;
}
*/
/*
function sectorCam(fromCam) {
    var cam = new Cam();
    cam.tm = fromCam.tm;
    cam.pos.fov = fromCam.pos.fov;
    cam.view = {start: {a: 0, r: 0}, center: {a: 0, r: 0, z: 0}, end: {a: 0, r: 0}};
    cam.setView = sectorSetView;
    cam.sectorMesh = [];
    cam.pushSectorMesh = pushSectorMesh;
    cam.shiftSectorMesh = shiftSectorMesh;
    cam.sectorName = '';
    cam.setSectorName = setSectorName;
    cam.update();
    cam.modulus = 1;
    cam._old_projection = cam.projection;
    cam.projection = sectorProjection;
    return cam;
}
*/
var Cam = (function () {
    function Cam(x, y, z, rx, ry, rz, fov) {
        this.pos = {
            sx: 0, sy: 0, sz: 0,
            rx: 0, ry: 0, rz: 0,
            fov: 0
        };
        this.depth = 0;
        this.limits = {
            x: cProps.w,
            y: cProps.h
        };
        this.tm = null;
        this.valid = false;
/*
        this.zoomToSector = null;
*/

        this.pos.fov = fov;
        this.depth = (cProps.h >> 1) / Math.tan(this.pos.fov * Math.PI / 360);

        this.initTM(x, y, z, rx, ry, rz);
    }

    Cam.prototype.initTM = function (x, y, z, rx, ry, rz) {
        var deg = Math.PI / 180;

        var cosx = Math.cos(rx * deg);
        var sinx = Math.sin(rx * deg);
        var cosy = Math.cos(ry * deg);
        var siny = Math.sin(ry * deg);
        var cosz = Math.cos(rz * deg);
        var sinz = Math.sin(rz * deg);

        var T = [
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [-x, -y, -z, 1]
        ];
        var X = [
            [1, 0, 0, 0],
            [0, cosx, sinx, 0],
            [0, -sinx, cosx, 0],
            [0, 0, 0, 1]
        ];
        var Y = [
            [cosy, 0, -siny, 0],
            [0, 1, 0, 0],
            [siny, 0, cosy, 0],
            [0, 0, 0, 1]
        ];
        var Z = [
            [cosz, sinz, 0, 0],
            [-sinz, cosz, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
        ];

        this.tm = normalizeMatrix(mulMatrix(mulMatrix(mulMatrix(T, Z), Y), X));
        this.isValid(false);
    };

    Cam.prototype.isValid = function (valid) {
        var ret = this.valid;

        if (typeof valid === 'undefined') {
            return ret;
        }

        this.valid = valid;

        return ret;
    };

    Cam.prototype.applyTransform = function () {
        var transform = [];
        var deg = Math.PI / 180;
		var E = 1E-7;

        if ((Math.abs(this.pos.sx) + Math.abs(this.pos.sy) + Math.abs(this.pos.sz)) > E) {
            transform.push(
                [
                    [1, 0, 0, 0],
                    [0, 1, 0, 0],
                    [0, 0, 1, 0],
                    [-this.pos.sx, -this.pos.sy, -this.pos.sz, 1]
                ]
            );

            this.pos.sx = this.pos.sy = this.pos.sz = 0;
        }

        if (Math.abs(this.pos.rz) > .001) {
            var cosz = Math.cos((this.pos.rz) * deg);
            var sinz = Math.sin((this.pos.rz) * deg);

            transform.push(
                [
                    [cosz, sinz, 0, 0],
                    [-sinz, cosz, 0, 0],
                    [0, 0, 1, 0],
                    [0, 0, 0, 1]
                ]
            );

            this.pos.rz = 0;
        }


        if (Math.abs(this.pos.ry) > .001) {
            var cosy = Math.cos((this.pos.ry) * deg);
            var siny = Math.sin((this.pos.ry) * deg);

            transform.push(
                [
                    [cosy, 0, -siny, 0],
                    [0, 1, 0, 0],
                    [siny, 0, cosy, 0],
                    [0, 0, 0, 1]
                ]
            );

            this.pos.ry = 0;
        }

        if (Math.abs(this.pos.rx) > .001) {
            var cosx = Math.cos((this.pos.rx) * deg);
            var sinx = Math.sin((this.pos.rx) * deg);

            transform.push(
                [
                    [1, 0, 0, 0],
                    [0, cosx, sinx, 0],
                    [0, -sinx, cosx, 0],
                    [0, 0, 0, 1]
                ]
            );

            this.pos.rx = 0;
        }

        if (transform.length) {
            var tm = this.tm;

            for (var i = 0; i < transform.length; i++) {
                tm = mulMatrix(tm, transform[i]);
            }

            this.tm = normalizeMatrix(tm);

            this.isValid(false);
        }
    };

    Cam.prototype.update = function () {
        this.applyTransform();

        var d = this.depth;
        this.depth = (cProps.h >> 1) / Math.tan(this.pos.fov * Math.PI / 360);
        if (d !== this.depth) {
            this.isValid(false);
        }
    };

    Cam.prototype.projection = function (pos) {
        var v3d = mulVectorMatrix([pos.x, pos.y, pos.z, 1], this.tm);

        if (Math.abs(v3d[2]) < .001)
            return {
                x: 0,
                y: 0,
                z: -1,
                xScale: 0,
                yScale: 0,
                inView: false
            };

        var scale = this.depth / (v3d[2]);
        var x2d = (v3d[0]) * scale;
        var y2d = (v3d[1]) * scale;

        // noinspection JSSuspiciousNameCombination
		return {
            x: x2d,
            y: y2d,
            z: v3d[2],
            xScale: scale,
            yScale: scale,
            inView: (Math.abs(x2d) < (this.limits.x * .999)) && (Math.abs(y2d) < (this.limits.y * .999))
        };
    };

/*
    Cam.prototype.zoomSector = function (sector) {
        var ret = this.zoomToSector;

        if (typeof sector === 'undefined') {
            return ret;
        }

        this.zoomToSector = sector;

        return ret;
    };
*/

    return Cam;
})();

//var guides = galaxyGuides();
var stars = null;

var routeStrokeStyles = [
	'#67ffe3',
	'#ff67e3',
	'#ffe367',
	'#e3ff67',
	'#e367ff',
	'#67e3ff'
];
var routes = null;
//var sectors = galaxySectors();
/*
var sectorMap = sectors.reduce(function (acc, val) {
    acc[val.name] = val;
    return acc;
}, {});
*/
var currentView;

var G = 6.67E-11;
var PK = 3.0856776E+16;
var C = 3E+8;
var CM = 2E+30 * 4.5E+6;
var CR = 2 * G * CM / (C * C) / PK;
//CR = 10;
var E = 1E-3;
var showBH = false;

// background
var back = document.createElement('canvas');
back.width = cProps.w;
back.height = cProps.h;
back.style.left = cProps.l;
back.style.top = cProps.t;
document.body.appendChild(back);

//foreground
var front = document.createElement('canvas');
front.width = cProps.w;
front.height = cProps.h;
front.style.left = cProps.l;
front.style.top = cProps.t;
document.body.appendChild(front);

//focus and location
var active = document.createElement('canvas');
active.width = cProps.w;
active.height = cProps.h;
active.style.left = cProps.l;
active.style.top = cProps.t;
document.body.appendChild(active);

// hud
var clip = document.createElement('canvas');
clip.width = cProps.hw;
clip.height = cProps.hh;
clip.style.left = cProps.hl;
clip.style.top = cProps.ht;
document.body.appendChild(clip);

var fContext = front.getContext('2d');
var aContext = active.getContext('2d');

// var btn1 = document.createElement('button');
// btn1.style.position = 'absolute';
// btn1.style.left = '0px';
// btn1.style.top = '0px';
// btn1.innerHTML = '<i>To Galaxy View</i>';
// btn1.addEventListener('click', function () {
//     hud.sectorCtrls(false);
//     switchView('galaxy');
//     currentCam = galaxyCam;
// /*
//     currentCam.zoomSector(null);
// */
// });
// document.body.appendChild(btn1);
//
// var btn2 = document.createElement('button');
// btn2.style.position = 'absolute';
// btn2.style.left = '0px';
// btn2.style.top = '22px';
// btn2.innerHTML = '<i>To System View</i>';
// btn2.addEventListener('click', function () {
//     hud.sectorCtrls(false);
//     switchView('system');
//     currentCam = systemCam;
// });
// document.body.appendChild(btn2);

// Mouse events and settings
var mousePos = {
    x: 0,
    y: 0,
    buttons: 0,
    update: function (ev, ctx) {
        if (ev) {
            this.buttons = ev.buttons;
            if (ctx) {
                this.x = ev.clientX - ctx.offsetLeft - cProps.dl;
                this.y = ev.clientY - ctx.offsetTop - cProps.dt;
            }
        }
    }
};
var lDragStarted = null;
/*
var sectorClick = function () {
    if (!(currentView === 'galaxy')) {
        return;
    }

    var sfTime = 2000;
    var seTime = 1500;
    var sector = sectorHitTest(aContext, currentCam);
    if (sector) {
        currentCam.zoomSector(sector);
        var scA = -sector.center.a - Math.PI * .5;
        var cosA = Math.cos(scA);
        var sinA = Math.sin(scA);
        var scm = mulMatrix(
            [
                [1, 0, 0, 0],
                [0, 1, 0, 0],
                [0, 0, 1, 0],
                [-sector.center.x, -sector.center.y, -sector.center.z, 1]
            ],
            [
                [cosA, sinA, 0, 0],
                [-sinA, cosA, 0, 0],
                [0, 0, 1, 0],
                [0, 0, 0, 1]
            ]
        );

        Tween.camTween(currentCam, Tween.easeInOut, copyMatrix(currentCam.tm), scm, sfTime,
            function () {
                console.log('Tween finished!');
                currentCam = sectorCam(currentCam);
                currentCam.setView(sector);
                currentCam.pushSectorMesh(sector);
                currentCam.setSectorName(sector.name);
                switchView('sector');
                Tween.tween(currentCam, 'modulus', Tween.easeInOut, 1, 0, seTime,
                    function () {
                        console.log('Sector modulus tween finished!');
                        hud.sectorCtrls(true);
                    },
                    function () {
                        currentCam.isValid(false);
                    });
            });
    }
};
*/
/*
clip.addEventListener('dblclick', sectorClick);
*/
function hudClickHandler() {
/*
    if (currentView === 'sector') {
		var button = hud.btnHitTest();
		if (button) {
			var seTime = 1200;
			var sector;

			switch (button.data) {
				case 'up':
				case 'down':
					sector = sectorMap[sectorMap[currentCam.sectorName].neighbours[button.data]];
					if (sector) {
						hud.sectorCtrls(false);
						currentCam.pushSectorMesh(sector);
						currentCam.setSectorName(sector.name);
						Tween.viewTween(currentCam, Tween.easeInOut, copyView(currentCam.view), copyView(sector), seTime,
							function () {
								console.log('Intersector tween finished!');
								currentCam.setView(sector);
								currentCam.shiftSectorMesh();
								hud.sectorCtrls(true);
							}
						);
					}
					break;
				case 'left':
				case 'right':
					sector = sectorMap[sectorMap[currentCam.sectorName].neighbours[button.data]];
					if (sector) {
						hud.sectorCtrls(false);
						currentCam.pushSectorMesh(sector);
						currentCam.setSectorName(sector.name);
						Tween.viewTween(currentCam, Tween.easeInOut, copyView(currentCam.view), copyView(sector), seTime,
							function () {
								console.log('Intersector tween finished!');
								currentCam.setView(sector);
								currentCam.shiftSectorMesh();
								hud.sectorCtrls(true);
							}
						);
					}
					break;
			}
		}
	}
*/
	var starProjection = starHitTest(aContext, currentCam);
	if (starProjection) {
		var star = stars[starProjection.star];
		if (router) {
			routes = router.getRoutes(star);
			console.log(routes);
			pLocation.target = starProjection.star;
			jQuery(window).trigger('routesUpdate');
		}
	}
}

var TouchHandler = (function () {
	var tapTimeMs = 200;

	var handlers = [];

	function findTouch(t) {
		return t.id === this.identifier;
	}

	function findHandler(h) {
		return h.$target.get(0) === this.target;
	}

	function onTouchStart(ev) {
		for (var i = 0; i < ev.changedTouches.length; i++) {
			var touch = {
				id: ev.changedTouches[i].identifier,
				clientX: ev.changedTouches[i].clientX,
				clientY: ev.changedTouches[i].clientY,
				start: new Date().getTime(),
				track: []
			};
			var touchi = this.touches.findIndex(findTouch, ev.changedTouches[i]);
			if (touchi >= 0) {
				this.touches[touchi] = touch;
			} else {
				this.touches.push(touch);
			}
		}
		if (this.touches.length === 1) {
			this.$target.trigger(
				$.Event('mousedown', {
					target: this.$target.get(0),
					clientX: touch.clientX,
					clientY: touch.clientY,
					button: 0
				})
			);
		} else {
			this.$target.trigger(
				$.Event('mouseup', {
					target: this.$target.get(0),
					clientX: this.touches[0].clientX,
					clientY: this.touches[0].clientY,
					button: 0
				})
			);
		}
		ev.preventDefault();
		ev.stopPropagation();
	}

	function onTouchMove(ev) {
		var tc = this.touches.length;
		for (var i = 0; i < ev.changedTouches.length; i++) {
			var touch = null,
				touchi = this.touches.findIndex(findTouch, ev.changedTouches[i]);
			if (touchi >= 0) {
				touch = this.touches[touchi];
			}
			if (touch) {
				touch.track.push({x: touch.clientX, y: touch.clientY});
				touch.clientX = ev.changedTouches[i].clientX;
				touch.clientY = ev.changedTouches[i].clientY;
			}
		}
		if (touch) {
			ev.preventDefault();
			ev.stopPropagation();
		}
		if (tc === 1 && touch) {
			this.$target.trigger(
				$.Event('mousemove', {
					target: this.$target.get(0),
					clientX: touch.clientX,
					clientY: touch.clientY,
					button: 0
				})
			);
		} else if (tc === 2) {
			var p1 = {x: this.touches[0].clientX, y: this.touches[0].clientY},
				p2 = {x: this.touches[1].clientX, y: this.touches[1].clientY},
				pp1, pp2;
			pp1 = p1;
			pp2 = p2;
			if (this.touches[0].track.length) {
				pp1 = this.touches[0].track[0];
			}
			if (this.touches[1].track.length) {
				pp2 = this.touches[1].track[0];
			}
			var pdx = (p1.x - p2.x),
				pdy = (p1.y - p2.y),
				ppdx = (pp1.x - pp2.x),
				ppdy = (pp1.y - pp2.y);
			var pd = Math.sqrt(pdx*pdx + pdy*pdy),
				ppd = Math.sqrt(ppdx*ppdx + ppdy*ppdy);

			var w = (ppd - pd) * .5;
			this.$target.trigger(
				$.Event('wheel', {
					target: this.$target.get(0),
					deltaX: w,
					deltaY: w
				})
			);
		}
	}

	function onTouchEnd(ev) {
		var ct = new Date().getTime();
		var tc = this.touches.length;
		for (var i = 0; i < ev.changedTouches.length; i++) {
			var touch = null,
				touchi = this.touches.findIndex(findTouch, ev.changedTouches[i]);
			if (touchi >= 0) {
				touch = this.touches[touchi];
			}
			this.touches.splice(touchi, 1);
		}
		if (touch) {
			ev.preventDefault();
			ev.stopPropagation();
		}
		if (touch && tc === 1) {
			this.$target.trigger(
				$.Event('mouseup', {
					target: this.$target.get(0),
					clientX: touch.clientX,
					clientY: touch.clientY,
					button: 0
				})
			);
			if ((ct - touch.start) < tapTimeMs) {
				this.$target.trigger(
					$.Event('click', {
						target: this.$target.get(0),
						clientX: touch.clientX,
						clientY: touch.clientY,
						button: 0
					})
				);
			}
		}
	}

	function TouchHandler (target) {
		this.$target = $(target);
		this.touches = [];
		target.addEventListener('touchstart', this.onStart = delegate(this, onTouchStart), {passive: false});
		window.addEventListener('touchmove', this.onMove = delegate(this, onTouchMove), {passive: false});
		window.addEventListener('touchend', this.onEnd = delegate(this, onTouchEnd), {passive: false});
	}

	return {
		handle: function (target) {
			handlers.push(new TouchHandler(target));
		},
		unhandle: function (target) {
			var iHandler = handlers.findIndex(findHandler, {target: target});
			if (iHandler >= 0) {
				var handler = handlers[iHandler];
				target.removeEventListener('touchstart', handler.onStart);
				window.removeEventListener('touchmove', handler.onMove);
				window.removeEventListener('touchend', handler.onEnd);
				handlers.splice(iHandler, 1);
			}
		}
	};
})();

TouchHandler.handle(clip);

$(clip).on('click', hudClickHandler)
	.on('mouseenter', function (ev) {
		if (ev.originalEvent)
			ev = ev.originalEvent;
		mousePos.update(ev, this);
	})
	.on('mousedown', function (ev) {
		if (ev.originalEvent)
			ev = ev.originalEvent;
		mousePos.update(ev, this);
		if (ev.button === 0) {
			lDragStarted = {x: mousePos.x, y: mousePos.y};
		}
	})
	.on('mouseup', function (ev) {
		if (ev.originalEvent)
			ev = ev.originalEvent;
		mousePos.update(ev, this);
		if (ev.button === 0) {
			lDragStarted = null;
		}
	})
	.on('mousemove', function (ev) {
		if (ev.originalEvent)
			ev = ev.originalEvent;
		mousePos.update(ev, this);

		var dx, dy, multiplier;
		if (lDragStarted) {

			dx = mousePos.x - lDragStarted.x;
			dy = mousePos.y - lDragStarted.y;

			lDragStarted = {x: mousePos.x, y: mousePos.y};

			multiplier = currentCam.tm[3][2] * Math.sin(currentCam.pos.fov * Math.PI / 360) / (cProps.h >> 1);

			currentCam.pos.sx -= dx * multiplier;
			currentCam.pos.sy -= dy * multiplier;

		}
	})
	.on('wheel', function (ev) {
		if (ev.originalEvent)
			ev = ev.originalEvent;
		switch (currentView) {
			case 'sector':
				return;
			case 'galaxy':
				currentCam.pos.sz = currentCam.pos.sz - ev.deltaY;
				ev.preventDefault();
				return;
			case 'system':
				currentCam.pos.sz = currentCam.pos.sz - ev.deltaY * 5;
				ev.preventDefault();
				return;
		}
	});

var background = {
    complete: false,
    w: 0,
    h: 0,
    pattern: null,
    pCtx: null,
    bgImg: null,
    bRoll: 2000,
    bCtx: null,
    x: 0,
    y: 0,
    r: Math.PI / 6,
    cx: 0,
    cy: 0,
    init: function (backCanvas) {
        this.bCtx = backCanvas.getContext('2d');
    },
    update: function () {
        if (!this.complete) {
            this.hudImg = MyCache.getImage("img/hud.png");
            if (!this.hudImg.complete) {
                return;
            }
            this.bgImg = MyCache.getImage("img/stars/star_bg2.jpg");
            if (!this.bgImg.complete) {
                return;
            }
            this.w = this.bgImg.naturalWidth;
            this.h = this.bgImg.naturalHeight;
            this.pattern = document.createElement('canvas');
            this.pattern.width = (this.bRoll << 1);
            this.pattern.height = (this.bRoll << 1);
            this.pCtx = this.pattern.getContext('2d');
            this.pCtx.fillStyle = this.pCtx.createPattern(this.bgImg, 'repeat');

            this.pCtx.setTransform(1, 0, 0, 1, 0, 0);
            this.pCtx.fillRect(0, 0, this.bRoll << 1, this.bRoll << 1);

            this.complete = true;
        }

        var fovmul = 45 / currentCam.pos.fov;
        var mul = cProps.h / 45;

        this.x = (this.x + (currentCam.pos.ry) * mul);
        this.y = (this.y - (currentCam.pos.rx) * mul);
        this.r = (this.r + (currentCam.pos.rz) * Math.PI / 180) % (2 * Math.PI);

        var sinr = Math.sin(-this.r);
        var cosr = Math.cos(-this.r);

        this.cx = (this.cx + this.x * cosr - this.y * sinr) % this.w;
        this.cy = (this.cy + this.x * sinr + this.y * cosr) % this.h;

        this.x = 0;
        this.y = 0;

        this.bCtx.setTransform(1, 0, 0, 1, (cProps.w >> 1), (cProps.h >> 1));
        this.bCtx.rotate(this.r);

        this.bCtx.drawImage(this.pattern,
            0, 0, this.pattern.width, this.pattern.height,
            fovmul * (this.cx - this.bRoll), fovmul * (this.cy - this.bRoll),
            fovmul * this.pattern.width, fovmul * this.pattern.height
        );
    }
};

background.init(back);

var Button = (function () {
    function Button(x, y, width, height, image, trans) {
        trans = trans || [1, 0, 0, 1, 0, 0];
        this.x = x;
        this.y = y;
        this.w = width;
        this.h = height;
        this.image = image;
        this.trans = trans;
        this.data = null;
        this.enabled = true;
    }

    Button.prototype.draw = function (ctx, state) {
        var sxm = 0;
        switch (state) {
            case 'hover':
                sxm = 1;
                break;
            case 'down':
                sxm = 2;
                break;
            case 'disabled':
                sxm = 3;
                break;
            case 'normal':
            default:
                sxm = 0;
                state = 'normal';
        }

        var matrix = this.trans.slice(0);
        matrix[4] += this.x;
        matrix[5] += this.y;
        ctx.setTransform.apply(ctx, matrix);
        ctx.drawImage(this.image.img, sxm * this.w, 0, this.w, this.h, 0, 0, this.w, this.h);
    };

    Button.prototype.setData = function (data) {
        this.data = data;

        return this;
    };

    Button.prototype.enable = function (state) {
        this.enabled = !!state;
    };

    return Button;
})();

var hud = {
    complete: false,
    hudImg: null,
    hudCtx: null,
    btnRightImg: {img: null},
    btnDownImg: {img: null},
    btnLeftImg: {img: null},
    btnUpImg: {img: null},
    valid: false,
    showBtns: false,
    btnSpace: 3,
    btnSizeX: 50,
    btnSizeY: 50,
    buttons: [],
    init: function (hudCanvas) {
        this.hudCtx = hudCanvas.getContext('2d');
        this.buttons.push(new Button(
            cProps.dl + cProps.w - this.btnSizeX - this.btnSpace,
            cProps.dt + (cProps.h >> 1) - (this.btnSizeY >> 1),
            this.btnSizeX,
            this.btnSizeY,
            this.btnRightImg).setData('right'));
        this.buttons.push(new Button(
            cProps.dl + (cProps.w >> 1) - (this.btnSizeX >> 1),
            cProps.dt + cProps.h - this.btnSizeY - this.btnSpace,
            this.btnSizeX,
            this.btnSizeY,
            this.btnDownImg).setData('down'));
        this.buttons.push(new Button(
            cProps.dl + this.btnSpace,
            cProps.dt + (cProps.h >> 1) - (this.btnSizeY >> 1),
            this.btnSizeX,
            this.btnSizeY,
            this.btnLeftImg).setData('left'));
        this.buttons.push(new Button(
            cProps.dl + (cProps.w >> 1) - (this.btnSizeX >> 1),
            cProps.dt + this.btnSpace,
            this.btnSizeX,
            this.btnSizeY,
            this.btnUpImg).setData('up'));
    },
    update: function () {
        if (!this.complete) {
            this.hudImg = MyCache.getImage("img/hud.png");
            this.btnRightImg.img = MyCache.getImage("img/btn_right.png");
            this.btnDownImg.img = MyCache.getImage("img/btn_down.png");
            this.btnLeftImg.img = MyCache.getImage("img/btn_left.png");
            this.btnUpImg.img = MyCache.getImage("img/btn_up.png");
            if (!(this.hudImg.complete
                    && this.btnRightImg.img.complete
                    && this.btnDownImg.img.complete
                    && this.btnLeftImg.img.complete
                    && this.btnUpImg.img.complete)
            ) {
                return;
            }
            this.complete = true;
        }

        this.hudCtx.setTransform(1, 0, 0, 1, 0, 0);
        this.hudCtx.clearRect(0, 0, cProps.hw, cProps.hh);
        this.hudCtx.drawImage(this.hudImg, 0, 0, cProps.hw, cProps.hh);

        // if (this.showBtns) {
        //     var btn = this.btnHitTest();
        //     var sector = sectorMap[currentCam.sectorName];
        //     for (var i = 0; i < this.buttons.length; i++) {
        //         this.buttons[i].enable(sector.neighbours[this.buttons[i].data] !== null);
        //         var state = 'normal';
        //         if (!this.buttons[i].enabled) {
        //             state = 'disabled';
        //         } else if (btn === this.buttons[i]) {
        //             state = (mousePos.buttons & 1) === 1 ? 'down' : 'hover';
        //         }
        //         this.buttons[i].draw(this.hudCtx, state);
        //     }
        // }

        this.isValid(true);
    },
    isValid: function (valid) {
        var ret = this.valid;

        if (typeof valid === 'undefined') {
            return ret;
        }

        this.valid = valid;

        return ret;
    }
/*
    ,
    sectorCtrls: function (show) {
        this.showBtns = !!show;
        this.isValid(false);
    },
	btnHitTest: function () {
		if (this.showBtns) {
			for (var i = 0; i < this.buttons.length; i++) {
				if ((mousePos.x + cProps.dl) > this.buttons[i].x
					&& (mousePos.x + cProps.dl) < (this.buttons[i].x + this.buttons[i].w)
					&& (mousePos.y + cProps.dt) > this.buttons[i].y
					&& (mousePos.y + cProps.dt) < (this.buttons[i].y + this.buttons[i].h)
				) {
					return this.buttons[i];
				}
			}
		}

		return null;
	}
*/
};

hud.init(clip);

function switchView(view) {
    switch (view) {
        case 'system':
        case 'sector':
        case 'galaxy':
            if (currentView !== view) {
                currentView = view;
                hud.isValid(false);
            }
    }
}

var starTypes = [];

function prepareStarTypes() {
    starTypes = [];
    var img = MyCache.getImage('img/stars/small.png');
    if (!img.complete) {
        return;
    }

    var sss = [
        {type: 'render', fs: ['rgba(68, 85, 51, 255)', 'rgba(68, 85, 51, 0)'], size: 2, gsize: 4, scale: 1},
        {type: 'image', tx: 256, ty: 128, tw: 128, th: 128, scale: 3}, // Class G
        {type: 'image', tx: 128, ty: 128, tw: 128, th: 128, scale: 2}, // Red giant
        {type: 'image', tx: 0, ty: 0, tw: 128, th: 128, scale: 4}, // White dwarf
        {type: 'image', tx: 256, ty: 0, tw: 128, th: 128, scale: 4}, // Yellow dwarf
        {type: 'image', tx: 0, ty: 128, tw: 128, th: 128, scale: 2}, // Blue giant
        {type: 'image', tx: 0, ty: 0, tw: 128, th: 128, scale: .2} // Galaxy core
    ];

    for (var i = 0; i < sss.length; i++) {
        var context;
        var obj = {
            canvas: document.createElement('canvas'),
            scale: sss[i].scale
        };

        if (sss[i].type === 'image') {
            obj.canvas.width = sss[i].tw;
            obj.canvas.height = sss[i].th;
            context = obj.canvas.getContext('2d');
            context.setTransform(1, 0, 0, 1, 0, 0);
            context.drawImage(img,
                sss[i].tx, sss[i].ty, sss[i].tw, sss[i].th,
                0, 0, sss[i].tw, sss[i].th);
        } else {
            var size = sss[i].size * obj.scale;
            var gs = sss[i].gsize * obj.scale;

            obj.canvas.width = gs * 2;
            obj.canvas.height = gs * 2;
            context = obj.canvas.getContext('2d');
            context.setTransform(1, 0, 0, 1, 0, 0);
            var grad = context.createRadialGradient(gs, gs, size, gs, gs, gs);
            grad.addColorStop(0, sss[i].fs[0]);
            grad.addColorStop(1, sss[i].fs[1]);
            context.fillStyle = grad;
            context.fillRect(0, 0, gs * 2, gs * 2);
        }

        starTypes.push(obj);
    }
}

var pLocation = {
	data: null,
	type: null,
    star: null,
	target: null,
	timer: null,
	src: null,
	total: null,
	start: null
};

function drawSelection(ctx, cam, obj, ss, lw, ccw) {
    var lPos = cam.projection(obj.pos);

    if (!(lPos.z > 0 && lPos.inView)) {
        return;
    }
    ccw = (!!ccw) ? -1 : 1;

    var x0 = cam.limits.x * .5;
    var y0 = cam.limits.y * .5;

    ctx.globalCompositeOperation = 'source-over';

    var a = ccw * 2 * Math.PI * ((Date.now() % 3000) / 3000),
        r = 24 * Math.max(.4, Math.min(1, lPos.xScale)),
        r1 = 24 * Math.max(.4, Math.min(1, lPos.xScale)),
        r2 = 30 * Math.max(.4, Math.min(1, lPos.xScale));

    ctx.setTransform(Math.cos(a), Math.sin(a), -Math.sin(a), Math.cos(a), x0 + lPos.x, y0 + lPos.y);

    ctx.strokeStyle = ss;
    ctx.lineWidth = lw;

    ctx.beginPath();
    ctx.arc(0, 0, r, 0, 2 * Math.PI);
    ctx.moveTo(r1, 0);
    ctx.lineTo(r2, 0);
    ctx.moveTo(0, r1);
    ctx.lineTo(0, r2);
    ctx.moveTo(-r1, 0);
    ctx.lineTo(-r2, 0);
    ctx.moveTo(0, -r1);
    ctx.lineTo(0, -r2);
    ctx.stroke();
}

function drawTransition(ctx, cam, src, dst, time, rTime) {
	var tPos = {
		x: src.pos.x + (dst.pos.x - src.pos.x) * time,
		y: src.pos.y + (dst.pos.y - src.pos.y) * time,
		z: src.pos.z + (dst.pos.z - src.pos.z) * time
	};
	var lPos = cam.projection(tPos),
		bPos = cam.projection(src.pos),
		ePos = cam.projection(dst.pos);

	if (!(lPos.z > 0 && lPos.inView)) {
		return;
	}

	var x0 = cam.limits.x * .5,
		y0 = cam.limits.y * .5,
		dr = 8 * Math.max(.4, Math.min(1, lPos.xScale));

	ctx.save();
	ctx.globalCompositeOperation = 'source-over';
	ctx.setTransform(1, 0, 0, 1, x0 + lPos.x, y0 + lPos.y);
	ctx.strokeStyle = '#ffffff';
	if (bPos.z > 0 && bPos.inView && ePos.z > 0 && ePos.inView) {
		ctx.lineWidth = .5;
		ctx.setLineDash([10, 4, 2, 4]);
		ctx.beginPath();
		ctx.moveTo(bPos.x - lPos.x, bPos.y - lPos.y);
		ctx.lineTo(ePos.x - lPos.x, ePos.y - lPos.y);
		ctx.stroke();
		ctx.setLineDash([]);
	}
	for (var i = 0; i < 3; i++) {
		let p = (Date.now() % 400) / 400;
		var r = (dr * i + dr * p) % (3 * dr);
		ctx.lineWidth = i < 2 ? 1 : (1 - p);
		ctx.beginPath();
		ctx.moveTo(r, 0);
		ctx.arc(0, 0, r, 0, 2 * Math.PI);
		ctx.stroke();
	}
	ctx.fillStyle = '#ffffff';
	ctx.textAlign = 'left';
	ctx.fillText(rTime, 2.16 * dr, -2.16 * dr);
	ctx.restore();
}

function drawLocation(ctx, cam) {
    if (currentView === 'system' || pLocation.star === null) {
        return;
    }

	if (pLocation.type === 3) {
		var duration = (Date.now() - pLocation.start) / pLocation.total;
		if (duration > 1) {
			drawSelection(ctx, cam, stars[pLocation.star], '#ffffff', 2);
		}
		else {
			var src = pLocation.src === null ? pLocation.star : pLocation.src,
				rTime = round(pLocation.total * (1 - duration) * 1E-3, 0).toString();
			if (typeof String.prototype.padStart !== 'undefined') {
				rTime = rTime.padStart(3, '0');
			}
			drawTransition(ctx, cam, stars[src], stars[pLocation.star], duration, rTime);
		}
	} else if (pLocation.type === 1) {
		drawSelection(ctx, cam, stars[pLocation.star], '#ffffff', 2);
	}

    if (pLocation.target !== null) {
		drawSelection(ctx, cam, stars[pLocation.target], '#5dff7e', 2, 1);
	}
}

function drawGizmo(ctx, cam) {
	var l = 40;
	var gizmo = mulMatrix([
		[1, 0, 0, 0],
		[0, 1, 0, 0],
		[0, 0, 1, 0],
		[0, 0, 0, 1]
	], cam.tm);

	var x0 = 40 + cProps.w * -.5;
	var y0 = 70 + cProps.h * -.5;

	ctx.globalAlpha = 1;
	ctx.strokeStyle = '#f00';
	ctx.beginPath();
	ctx.moveTo(x0, y0);
	ctx.lineTo(x0 + l * gizmo[0][0], y0 + l * gizmo[0][1]);
	ctx.stroke();

	ctx.strokeStyle = '#0f0';
	ctx.beginPath();
	ctx.moveTo(x0, y0);
	ctx.lineTo(x0 + l * gizmo[1][0], y0 + l * gizmo[1][1]);
	ctx.stroke();

	ctx.strokeStyle = '#00f';
	ctx.beginPath();
	ctx.moveTo(x0, y0);
	ctx.lineTo(x0 + l * gizmo[2][0], y0 + l * gizmo[2][1]);
	ctx.stroke();
}

/*
function polar(pos) {
    var a = Math.atan2(pos.y, pos.x);
    return {
        r: Math.sqrt(pos.x * pos.x + pos.y * pos.y),
        a: a >= 0 ? a : a + 2 * Math.PI
    };
}
*/

var Preloader = (function () {
    function Preloader() {
        this.progress = 0;
        this.action = '';
        this.progressFn = false;
        this.deactivateTimeout = false;
    }

    Preloader.prototype.draw = function (ctx) {
        var w = 200;
        var h = 20;

        ctx.setTransform(1, 0, 0, 1, cProps.w >> 1, cProps.h >> 1);
        ctx.globalCompositeOperation = 'source-over';
        ctx.strokeStyle = '#dddddd';
        ctx.lineWidth = 2;
        ctx.fillStyle = '#5e59dd';

        ctx.strokeRect(-(w >> 1), -(h >> 1), w, h);
        ctx.fillRect(2 - (w >> 1), 2 - (h >> 1), (w - 4), h - 4);
        ctx.fillStyle = '#a0a0a0';
        ctx.fillRect(2 - (w >> 1), 2 - (h >> 1), (w - 4) * this.progress, h - 4);

        ctx.fillStyle = '#ffffff';
        ctx.textAlign = 'center';
        ctx.fillText((this.action.length ? this.action + ' ' : '') + Math.round(this.progress * 100) + '%', 0, (h >> 1) - 5);
    };

    Preloader.prototype.update = function (ctx) {
        this.action = '';
        this.progress = 0;

        if (this.progressFn) {
            var a = this.progressFn();
            this.action = a.action;
            this.progress = a.progress;
            this.draw(ctx);
            if (this.progress >= 1 && !this.deactivateTimeout) {
                this.deactivateTimeout = setTimeout(delegate(this, this.deactivate), 500);
            } else if (this.progress < 1 && this.deactivateTimeout) {
                clearTimeout(this.deactivateTimeout);
                this.deactivateTimeout = false;
            }
        }
    };

    Preloader.prototype.activate = function (progressFn) {
        this.progressFn = progressFn;
    };

    Preloader.prototype.deactivate = function () {
    	this.deactivateTimeout = false;
        this.progressFn = false;
    };

    return Preloader;
}());


switchView('galaxy');

var frameCounter = 0;
var frameTimer = 0;
var FPS = 0;
var beginTime = Date.now();

function main(data) {
    // noinspection JSUnresolvedVariable
    if (typeof data.stars !== 'undefined') {
        // noinspection JSUnresolvedVariable
        stars = data.stars;
		initControls();
	}

    var dt = Date.now() - beginTime;
    beginTime = Date.now();
    frameTimer += dt;
    frameCounter++;
    if (frameTimer > 1000) {
        FPS = frameCounter * 1000 / frameTimer;
        frameCounter = 0;
        frameTimer = 0;
    }

    background.update();
    hud.update();
    currentCam.update();
    if (starTypes.length) {
        if (!currentCam.isValid()) {
            fContext.setTransform(1, 0, 0, 1, 0, 0);
            fContext.clearRect(0, 0, front.width, front.height);
            calcProjections(currentCam);
            drawOrbits(fContext, currentCam);
            draw(fContext, currentCam);
            drawTunnels(fContext, currentCam);
            drawGizmo(fContext, currentCam);
            currentCam.isValid(true);
        }
        aContext.setTransform(1, 0, 0, 1, 0, 0);
        aContext.clearRect(0, 0, front.width, front.height);
        aContext.fillStyle = '#999';
        aContext.textAlign = 'left';
        aContext.fillText("FPS: " + round(FPS, 2), 37, 32);
        aContext.fillText("FOV, DEP: " + round(currentCam.pos.fov) + ', ' + round(currentCam.depth), 37, 42);
    } else {
        prepareStarTypes();
    }
/*
    drawFocus(aContext, currentCam);
    drawSectorMesh(aContext, currentCam);
*/
    drawLocation(aContext, currentCam);
    drawStarFocus(aContext, currentCam);
    drawRoutes(aContext, currentCam);
    preloader.update(aContext);

    requestAnimationFrame(main);
}

function updateRouterData(star) {
	router.setLocation(star);
	preloader.activate(delegate(router, router.progressFn));
}

function onLocationLoad(data) {
	pLocation.data = data.place;
	jQuery(window).trigger({
		type: 'locationLoaded',
		location: data.place
	});
	var current = pLocation.star;
	var curType = pLocation.type;
	var i, location, start;
	var place = data.place;
	pLocation.type = place.type;
	var reloadTime = 60000;
	switch (place.type) {
		case 3: // in motion
			start = place.source;
			while (start !== null && start.type !== 2) {
				start = start.parent || null;
			}
			if (start !== null && start.type === 2) {
				for (i = 0; i < stars.length; ++i) {
					if (stars[i].id === start.uid) {
						pLocation.src = i;
						break;
					}
				}
			}
			pLocation.start = Date.now() - (place.total - place.duration) * 1000;
			pLocation.total = place.total * 1000;
			reloadTime = Math.max(500, Math.min(reloadTime, place.duration * 1000));
			location = place.destination;
			break;
		case 1: // orbit
			location = place.params;
			break;
		case 2: // in dock
		default:
			location = null;
	}

	while (location !== null && location.type !== 2) {
		location = location.parent || null;
	}

	if (location !== null && location.type === 2) {
		for (i = 0; i < stars.length; i++) {
			if (stars[i].id === location.uid) {
				pLocation.star = i;
				break;
			}
		}

		if (!router) {
			router = new Router(pLocation.star);
			preloader.activate(delegate(router, router.progressFn));
		} else if ((curType === 3 && pLocation.type === 1 && pLocation.src !== pLocation.star)
			|| (pLocation.type === 1 && current !== pLocation.star)) {
			console.log('Location changed. Rebuilding routes');
			updateRouterData(pLocation.star);
		} else {
			console.log('Location unchanged.');
		}
	}

	pLocation.timer = setTimeout(reloadLocation, reloadTime);
}

var system = null;

function onSystemLoad(data) {
    if (typeof data.children !== 'undefined') {
        system = data;
    }
}

var galaxyCam = new Cam(1400, 2666, -2000, -55, 12, 20, 45);
//var galaxyCam = new Cam(0, 0, -5000, 0, 0, 0, 45);
//var systemCam = new Cam(1000, 5000, -3000, -55, 15, 0, 45);
var currentCam = galaxyCam;
var router = null;
var preloader = new Preloader();

function reloadLocation() {
	if (pLocation.timer) {
		clearTimeout(pLocation.timer);
		pLocation.timer = null;
	}

	jQuery.ajax({
		url: 'mock/whereami.json',
		dataType: 'json',
		jsonp: false,
		crossDomain: true,
		xhrFields: {
			withCredentials: true
		}
	}).done(onLocationLoad);
}

jQuery.ajax({
    url: 'resources/galaxy.json',
    dataType: 'json',
    jsonp: false,
    crossDomain: true,
    xhrFields: {
        withCredentials: true
    }
}).done(main).then(reloadLocation);

jQuery.getJSON('resources/system.json').done(onSystemLoad);

function onRouterUpdate() {
	currentCam.isValid(false);
	if (pLocation.target) {
		routes = router.getRoutes(stars[pLocation.target]);
		console.log(routes);
		jQuery(window).trigger('routesUpdate');
	}
}

function initControls() {
    console.log('initControls');

	var $body = jQuery('body');

	var $focusInfo = jQuery('<div/>').addClass('focusInfo');
	var $locationInfo = jQuery('<div/>').addClass('locInfo');
	var $nextHopInfo = jQuery('<div/>').addClass('hopInfo');

	jQuery('<div/>')
		.addClass('display')
		.css({
			left: (parseInt(cProps.hl) + 30) + 'px',
			top: (parseInt(cProps.ht) + 39) + 'px'
		})
		.append($locationInfo)
		.append($focusInfo)
		.append($nextHopInfo)
		.appendTo($body);

	function updateLocationInfo(ev) {
		if (typeof ev.location === 'undefined') {
			$locationInfo.empty();
			return;
		}

		var i, t = ['Location is not defined'];
		if (ev.location.type === 3) {
			var start = ev.location.source;
			while (start !== null && start.type !== 2) {
				start = start.parent || null;
			}
			if (start !== null && start.type === 2) {
				for (i = 0; i < stars.length; ++i) {
					if (stars[i].id === start.uid) {
						start = i;
						break;
					}
				}
			}
			var finish = ev.location.destination;
			while (finish !== null && finish.type !== 2) {
				finish = finish.parent || null;
			}
			if (finish !== null && finish.type === 2) {
				for (i = 0; i < stars.length; ++i) {
					if (stars[i].id === finish.uid) {
						finish = i;
						break;
					}
				}
			}

			if (start !== null && start === finish) {
				t = [
					'<strong>Transferring:</strong>&nbsp;',
					'<strong>System:</strong>&nbsp;' + stars[start].name
				];
			} else if (start !== null && finish !== null) {
				t = [
					'<strong>Transferring:</strong>&nbsp;',
					'<strong>From:</strong>&nbsp;' + stars[start].name,
					'<strong>To:</strong>&nbsp;' + stars[finish].name
				];
			}
		} else if (ev.location.type === 1) {
			var loc = ev.location.params;
			if (typeof loc.uid === 'undefined') loc.uid = ev.location.uid;
			while (loc !== null && loc.type !== 2) {
				loc = loc.parent || null;
			}
			if (loc !== null && loc.type === 2) {
				for (i = 0; i < stars.length; ++i) {
					if (stars[i].id === loc.uid) {
						loc = i;
						break;
					}
				}
			}
			if (loc !== null) {
				t = ['<strong>System:</strong>&nbsp;' + stars[loc].name];
			}
		}

		if (router) {
			t.push('<strong>Available:</strong>&nbsp;' + router.mapData.map.filter(function(i){return i.length > 0}).length);
		}

		$locationInfo.html(t.join('<br/>'));
	}

	function updateFocusInfo(ev) {
		if (ev.star === $focusInfo.data('star')) {
			return;
		}

		$focusInfo.data('star', ev.star);

		if (ev.star < 0 || typeof stars[ev.star] === 'undefined') {
			$focusInfo.empty();
			return;
		}

		var t = ['<strong>Target:</strong>&nbsp;' + stars[ev.star].name];

		if (pLocation.star !== null) {
			t.push('<strong>Distance:</strong>&nbsp;'
				+ round(distance(stars[pLocation.star].pos, stars[ev.star].pos), 2)
				+ '&nbsp;pc');

			if (router) {
				if (router.mapData.map[ev.star].length) {
					t.push('(available)');
				} else {
					t.push('(not available)');
				}
			}
		}

		$focusInfo.html(t.join('<br/>'));
	}

	function updateNextHopInfo() {
		$nextHopInfo.empty();

		if (pLocation.data === null || pLocation.target === null) {
			return;
		}

		if (pLocation.data.type === 1) {
			var i, loc = pLocation.data.params;
			if (typeof loc.uid === 'undefined') loc.uid = pLocation.data.uid;
			if (loc !== null && loc.type === 8) {
				while (loc !== null && loc.type !== 2) {
					loc = loc.parent || null;
				}
				if (loc !== null && loc.type === 2) {
					for (i = 0; i < stars.length; ++i) {
						if (stars[i].id === loc.uid) {
							loc = i;
							break;
						}
					}
				}
				if (loc !== null && routes.length) {
					var e, j, k, next = [];
					for (i = 0; i < routes.length; ++i) {
						e = false;
						for (j = 1; j < routes[i].path.length; ++j) {
							if (routes[i].path[j] !== loc) {
								for (k = 0; k < next.length; ++k) {
									if (next[k] === routes[i].path[j]) {
										break;
									}
								}
								if (k >= next.length) {
									next.push(routes[i].path[j]);
								}
								break;
							}
						}
					}
					if (next.length) {
						jQuery(next).each(function (i, v) {
							$nextHopInfo.append(
								jQuery('<button/>')
									.addClass('btnJump')
									.append('Transfer:&nbsp;' + stars[v].name)
									.bind('click', function () {
										jQuery.ajax({
											url: 'mock/dummy.json',
											dataType: 'html',
											jsonp: false,
											crossDomain: true,
											xhrFields: {
												withCredentials: true
											}
										}).done(reloadLocation);
									}));
						});
						return;
					}
				}
			}
		}

		$nextHopInfo.html('Hyperjump not available');
	}

	jQuery(window).on('starHover', updateFocusInfo);
	jQuery(window).on('locationLoaded', updateLocationInfo);
	jQuery(window).on('locationLoaded', updateNextHopInfo);
	jQuery(window).on('routesUpdate', updateNextHopInfo);

    var $starDD = jQuery('<input/>')
        .addClass('loc')
        .attr({
            type: 'text'
        });

    var $pathDisplay = jQuery('<div/>').addClass('routeText');

    jQuery('<div/>')
        .addClass('controls')
        .css({
            left: (parseInt(cProps.hl) + 30) + 'px',
            top: (parseInt(cProps.ht) + 274) + 'px'
        })
        .append([
            jQuery('<div>')
                .append([
                    $starDD,
                    jQuery('<button/>')
                        .addClass('btnQ')
                        .append('&raquo;')
                        .bind('click', function () {
                            var i = $starDD.attr('data-value');
                            if (i && stars[i]) {
                                var star = stars[i];
                                if (router) {
	                                routes = router.getRoutes(star);
	                                console.log(routes);
									pLocation.target = i;
	                                jQuery(window).trigger('routesUpdate');
                                }
                            }
                        })
                ]),
            jQuery('<div>').addClass('routeTextScroller').append($pathDisplay),
			jQuery('<input/>').addClass('showTun').attr({id: 'showTunCkBox', type: 'checkbox'}).on('change',
				function() {
					if (router) {
						router.mapData.show = jQuery(this).is(':checked');
						currentCam.isValid(false);
					}
				}
			),
			jQuery('<label/>').addClass('showTun').attr({for: 'showTunCkBox'}).text('hypertunnels')
		])
        .appendTo($body);

    $starDD.autocomplete({
        minLength: 2,
        source: function (req, res) {
            var found = [];
            var re = new RegExp(jQuery.ui.autocomplete.escapeRegex(req.term), 'i');
            for (var i = 0, n = stars.length; i < n; ++i) {
                if ((typeof stars[i].name !== 'undefined') && stars[i].name.match(re)) {
                    found.push({i: stars[i].i});
                }
                if (found.length >= 20) {
                    break;
                }
            }
            res(found);
        },
        focus: function (event, ui) {
            $starDD.val(stars[ui.item.i].name);
            return false;
        },
        select: function (event, ui) {
            $starDD
                .attr('data-value', ui.item.i)
                .val(stars[ui.item.i].name);
            return false;
        }
    }).autocomplete('instance')._renderItem = function (ul, item) {
        return jQuery('<li/>')
            .append(stars[item.i].name)
            .appendTo(ul);
    };

    function updateRoutes() {
		var target = null;
		var tt = 'Routes not found';
		if (routes.length > 0) {
			var ttt = routes.length % 10;
			if (routes.length > 10 && routes.length < 20) {
				ttt += 10;
			}
			switch (ttt) {
				case 1:
					tt = 'Found ' + routes.length + ' route';
					break;
				default:
					tt = 'Found ' + routes.length + ' routes';
			}
		}
        var text = ['<span style="font-weight: bold;">' + tt + '</span>'];
        for (var j = 0; j < routes.length; ++j) {
            text.push('<span style="font-weight: bold; color: ' + routeStrokeStyles[j % 6] + ';"> Route #' + (j + 1)
                + ', jumps: ' + (routes[j].path.length - 1)
                + ', distance: ' + round(routes[j].length, 2) + 'pc</span>'
                + '<span style="color: ' + routeStrokeStyles[j % 6] + ';">');
            var lastDist = 0;
            for (var k = 1; k < routes[j].path.length; ++k) {
                text.push('  ' + stars[routes[j].path[k]].name + ' [' + round(routes[j].dist[k] - lastDist, 2) + ']');
                target = routes[j].path[k];
                lastDist = routes[j].dist[k];
            }
            text.push('</span>');
        }
        $pathDisplay.html(text.join('<br/>'));
        if (pLocation.target) {
            $starDD.attr('data-value', pLocation.target).val(stars[pLocation.target].name);
        }
    }

    jQuery(window).on('routesUpdate', updateRoutes);
    jQuery(window).on('routerUpdated', onRouterUpdate)
}

jQuery('<button>')
	.text('Update position')
	.css({
		position: 'absolute',
		left: '0px',
		top: '0px'
	})
	.on('click', reloadLocation)
	.appendTo(jQuery('body'));
