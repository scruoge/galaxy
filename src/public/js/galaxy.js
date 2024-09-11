(function () {
	Gaussian = {
		next: 0,
		hasNext: false,
		getNext: function () {
			// polar method of G. E. P. Box, M. E. Muller, and G. Marsaglia
			// See Knuth, ACP, Section 3.4.1 Algorithm C.
			if (this.hasNext) {
				this.hasNext = false;

				return this.next;
			}

			let v1, v2, s, m;

			do
			{
				v1 = Math.random() * 2 - 1;  // between -1 and 1
				v2 = Math.random() * 2 - 1;  // between -1 and 1
				s = v1 * v1 + v2 * v2;
			}
			while (s >= 1 || s <= 1E-7);

			m = Math.sqrt(-2 * Math.log(s) / s);
			this.next = v2 * m;
			this.hasNext = true;

			return v1 * m;
		}
	};
})();

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

function mulVectorMatrix(v, m) {
	return [
		v[0] * m[0][0] + v[1] * m[1][0] + v[2] * m[2][0] + v[3] * m[3][0],
		v[0] * m[0][1] + v[1] * m[1][1] + v[2] * m[2][1] + v[3] * m[3][1],
		v[0] * m[0][2] + v[1] * m[1][2] + v[2] * m[2][2] + v[3] * m[3][2],
		v[0] * m[0][3] + v[1] * m[1][3] + v[2] * m[2][3] + v[3] * m[3][3]
	];
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

	if (l > 0) {
		v1[0] = v[0] / l;
		v1[1] = v[1] / l;
		v1[2] = v[2] / l;
	}

	return v1;
}

function rotateVector(v, dir, angle) {
	let sinD = Math.sin(dir),
		cosD = Math.cos(dir),
		sinA = Math.sin(angle),
		cosA = Math.cos(angle);
	let vz = [cosD, sinD, 0, 0];
	let vx = normalizeVector(crossVectorVector(v, vz));
	let vy = normalizeVector(crossVectorVector(vz, vx));
	let vMatrix = [vx, vy, vz, [0, 0, 0, 1]];
	let invVMatrix = normalizeMatrix(inverseMatrix(vMatrix));
	let n = mulVectorMatrix(v, invVMatrix);
	n = mulVectorMatrix(n, [
		[cosA, -sinA, 0, 0],
		[sinA, cosA, 0, 0],
		[0, 0, 1, 0],
		[0, 0, 0, 1]]);

	return mulVectorMatrix(n, vMatrix);
}

function ellipseGalaxy2(stars, time) {
	let arr = [
		{
			cls: 5,
			pos: {x: 0, y: 0, z: 0},
			c: {sA: 0, sE: 0, sI: 0, sID: 0, sOP: 0, sCP: 0}
		}
	];
	let r = 1600 << 1;
	let tt = -Math.log(100);
	let n = (typeof stars === 'number');
	let nStars = n ? stars : stars.length;
	let minSA = r * .02;
	let maxSA = minSA + r * .5;
	for (let i = 1; i < nStars; ++i) {
		let sA, // Orbit big radius
			sE, // Eccentricity
			sI, // Inclination
			sID, // Inclination directive
			sPhi, // Precession per period
			sP, // Orbital period
			sIOT, // Initial orbital time (0..1)
			sIP; // Initial precession (angle)

		let tP, // temporary value
			cR, // Current distance
			M,  // Mean anomaly
			E,  // Eccentric anomaly
			tA; // True anomaly

		if (n) {
			let kSA = Math.abs(Gaussian.getNext()),
				kSB = Math.abs(Gaussian.getNext());
			sA = minSA + 0.25 * r * kSA; // random
			//tP = Math.PI + (1 - kSA *.25) * (4 - Math.PI);
			//sE = Math.sqrt(Math.max(0, (tP - Math.PI)/(tP + Math.PI)));
			sE = Math.random() * 0.3
			sPhi = -1 / Math.exp(1/6 * Math.log(sA * (1 - sE * sE)));
			sP = Math.exp(1.5 * Math.log(sA));
			sIP = Math.round(Math.random()) * Math.PI;
			sIOT = Math.random();
			sI = Math.PI * .25 * kSB / Math.exp(4 * Math.log(1 + kSA));
			sID = Math.PI * Math.random() * Math.PI * 2;
		} else {
			sA = stars[i].c.sA;
			sE = stars[i].c.sE;
			sPhi = stars[i].c.sPhi;
			sP = stars[i].c.sP;
			sIP = stars[i].c.sIP;
			sIOT = stars[i].c.sIOT;
			sI = stars[i].c.sI;
			sID = stars[i].c.sID;
		}

		let periods = (time) / sP + sIOT;
		M = (periods - Math.floor(periods)) * 2 * Math.PI; // Mean anomaly
		E = M;
		let dE = 0;
		do {
			E = M + dE;
			dE = sE * Math.sin(E);
		} while (Math.abs(M - E + dE) > 1E-8);

		tA = 2 * Math.atan(Math.sqrt((1 + sE) / (1 - sE)) * Math.tan(E / 2)); // true anomaly

		cR = sA * (1 - sE * Math.cos(E));

		let sCP = periods * sPhi + sIP // Current precession
		while (sCP > 2 * Math.PI) {
			sCP -= 2 * Math.PI;
		}

		let nS = rotateVector([Math.cos(tA) * cR, Math.sin(tA) * cR, 0, 1], sID, sI);

		let rS = mulVectorMatrix(nS, [
			[Math.cos(sCP), -Math.sin(sCP), 0, 0],
			[Math.sin(sCP), Math.cos(sCP), 0, 0],
			[0, 0, 1, 0],
			[0, 0, 0, 1]]);

		if (n) {
			let obj = {
				cls: Math.round(Math.random() * 4),
				pos: {x: rS[0], y: rS[1], z: rS[2]},
				c: {sA: sA, sE: sE, sI: sI, sID: sID, sIOT: sIOT, sPhi: sPhi, sIP: sIP, sP: sP}
			};
			arr.push(obj);
		} else {
			stars[i].pos = {x: rS[0], y: rS[1], z: rS[2]};
		}
	}

	if (n) {
		return arr;
	} else {
		return stars;
	}
}

function galaxyGuides() {
	var arr = [];

	var step = mult * 8;

	for (var i = 0; i <= 12; i++) {
		if (i > 0) {
			var r = i * step * 10;
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

function projection(pos, cam) {
	var v3d = mulVectorMatrix([pos.x, pos.y, pos.z, 1], cam.tm);

	var scale = cam.depth / (v3d[2]);
	var x2d = (v3d[0]) * scale;
	var y2d = (v3d[1]) * scale;

	return {
		x: x2d,
		y: y2d,
		z: v3d[2],
		xScale: scale,
		yScale: scale,
		inView: (Math.abs(x2d) < (cam.limits.x * .999)) && (Math.abs(y2d) < (cam.limits.y * .999))
	};
}

function compareObjZ(a, b) {
	return a.pos.z - b.pos.z;
}

function drawObject(context, obj) {
	let m = starScale,
		src = starTypes[obj.cls],
		dw = m * src.canvas.width * obj.pos.xScale / src.scale,
		dh = m * src.canvas.height * obj.pos.yScale / src.scale;

	if (obj.cls !== 5) {
		dw = Math.max(4, dw / Math.max(1, 2 * Math.log(obj.pos.xScale)));
		dh = Math.max(4, dh / Math.max(1, 2 * Math.log(obj.pos.yScale)));
	}

	context.drawImage(src.canvas, obj.pos.x - dw * .5, obj.pos.y - dh * .5, dw, dh);
}

function drawOrbits (context, cam) {
	let i, a, pos;
	var x0 = cam.limits.x *.5;
	var y0 = cam.limits.y *.5;

	context.setTransform(1, 0, 0, 1, x0, y0);
	context.globalCompositeOperation = 'screen';
	context.lineWidth = 0.2;
	context.strokeStyle = '#ffffff';

	for (i = 1; i < stars.length; i++) {
		if (stars[i].c.sA > maxOrbit) {
			continue;
		}
		context.beginPath();
		let started = false;
		for (a = 0; a < Math.PI * 2; a += Math.PI * .01) {
			let nR = stars[i].c.sA * (1 - stars[i].c.sE * stars[i].c.sE) / (1 + stars[i].c.sE * Math.cos(a));
			let vS = rotateVector([Math.cos(a) * nR, Math.sin(a) * nR, 0, 1], stars[i].c.sID, stars[i].c.sI);

			let rS = mulVectorMatrix(vS, [
				[Math.cos(stars[i].c.sCP), -Math.sin(stars[i].c.sCP), 0, 0],
				[Math.sin(stars[i].c.sCP), Math.cos(stars[i].c.sCP), 0, 0],
				[0, 0, 1, 0],
				[0, 0, 0, 1]]);

			pos = projection({x: rS[0], y: rS[1], z: rS[2]}, cam);
			context[started ? 'lineTo' : 'moveTo'](pos.x, pos.y);
			started = true;
		}
		context.closePath();
		context.stroke();
	}
}

function draw(context, cam) {
	var i;
	var x0 = cam.limits.x *.5;
	var y0 = cam.limits.y *.5;

	var pos;

	var aPos = [];

	for (i = 0; i < stars.length; i++) {
		pos = projection(stars[i].pos, cam);
		if (pos.z > 0 && pos.inView) {
			aPos.push({pos: pos, cls: stars[i].cls});
		}
	}

	aPos.sort(compareObjZ);

	context.setTransform(1, 0, 0, 1, x0, y0);
	context.globalCompositeOperation = 'screen';

	for (i = 0; i < aPos.length; i++) {
		drawObject(context, aPos[i]);
	}
}

var cProps = {w: 0, h: 0};
var frameCounter = 0;
var frameTimer = 0;
var FPSCounter = 0;
var beginTime = Date.now();
var pathTimer = 0;

var mult = 4;
var cam = {
	pos: {
		rx: 0,
		ry: 0,
		rz: 0,

		fov: 0
	},
	depth: 0,
	limits: {
		x: cProps.w,
		y: cProps.h
	},
	tm: null,
	init: function(x, y, z, rx, ry, rz, fov) {
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

		this.limits = {
			x: cProps.w,
			y: cProps.h
		};

		this.tm = mulMatrix(mulMatrix(mulMatrix(T, Z), Y), X);

		this.pos.fov = fov;
		this.depth = (cProps.h >> 1) / Math.tan(this.pos.fov * deg / 2)
	},
	applyTransform: function () {
		var transform = [];
		var deg = Math.PI / 180;

		var s = this.shift.update();

		if (Math.abs(s.x + s.y + s.z) > .001) {
			transform.push(
				[
					[1, 0, 0, 0],
					[0, 1, 0, 0],
					[0, 0, 1, 0],
					[-s.x, -s.y, -s.z, 1]
				]
			);
		}

		if (Math.abs(this.pos.rz + this.rotation.rz) > .001) {
			var cosz = Math.cos((this.pos.rz + this.rotation.rz) * deg);
			var sinz = Math.sin((this.pos.rz + this.rotation.rz) * deg);

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


		if (Math.abs(this.pos.ry + this.rotation.ry) > .001) {
			var cosy = Math.cos((this.pos.ry + this.rotation.ry) * deg);
			var siny = Math.sin((this.pos.ry + this.rotation.ry) * deg);

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

		if (Math.abs(this.pos.rx + this.rotation.rx) > .001) {
			var cosx = Math.cos((this.pos.rx + this.rotation.rx) * deg);
			var sinx = Math.sin((this.pos.rx + this.rotation.rx) * deg);

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

			this.tm = tm;
		}
	},
	update: function () {
		this.applyTransform();

		this.depth = (cProps.h >> 1) / Math.tan(this.pos.fov * Math.PI / 360);
	},
	accel: function(axis, dir) {
		this.shift.accel(axis, dir)
	},
	shift: {
		sx: 0,
		sy: 0,
		sz: 0,
		tx: 0,
		ty: 0,
		tz: 0,
		a: 5,
		maxS: 500,
		update: function () {
			this.sx = (this.sx == this.tx)
				? this.sx
				: ((this.sx > this.tx)
					? Math.max(this.tx, this.sx - this.a)
					: Math.min(this.tx, this.sx + this.a));
			this.sy = (this.sy == this.ty)
				? this.sy
				: ((this.sy > this.ty)
					? Math.max(this.ty, this.sy - this.a)
					: Math.min(this.ty, this.sy + this.a));
			this.sz = (this.sz == this.tz)
				? this.sz
				: ((this.sz > this.tz)
					? Math.max(this.tz, this.sz - this.a)
					: Math.min(this.tz, this.sz + this.a));

			return {x: this.sx, y: this.sy, z: this.sz};
		},
		accel: function (axis, dir) {
			var _var;

			switch (axis) {
				case 'x':
				case 'y':
				case 'z':
					_var = 't' + axis;
					break;
				default:
					return;
			}

			switch (true) {
				case (dir > 0):
					this[_var] = this.maxS;
					break;
				case (dir < 0):
					this[_var] = -this.maxS;
					break;
				case (dir == 0):
					this[_var] = 0;
					break;
			}
		}
	},
	rotate: function (axis, dir) {
		this.rotation.rotate(axis, dir);
	},
	rotation: {
		rS: 1,
		rx: 0,
		ry: 0,
		rz: 0,
		rotate: function(axis, dir) {
			var _var;

			switch (axis) {
				case 'x':
				case 'y':
				case 'z':
					_var = 'r' + axis;
					break;
				default:
					return;
			}

			switch (true) {
				case (dir > 0):
					this[_var] = this.rS;
					break;
				case (dir < 0):
					this[_var] = -this.rS;
					break;
				case (dir == 0):
					this[_var] = 0;
					break;
			}
		}
	}
};

var timeScale = .005;
var timeZero = 3.3E5;

var next = 0;
var stars = 1.5E4;
var orbits = false;
var maxOrbit = 100;
var starScale = .25;

if (!(typeof module !== 'undefined' && module.exports)) {
	waitImagesLoading();

	cProps = {
		w: document.body.clientWidth,
		h: document.body.clientHeight
	};

// background
	var back = document.createElement('canvas');
	back.width = cProps.w;
	back.height = cProps.h;
	document.body.appendChild(back);

//foreground
	var front = document.createElement('canvas');
	front.width = cProps.w;
	front.height = cProps.h;
	document.body.appendChild(front);
	var fContext = front.getContext('2d');

// Mouse events and settings
	var dragStarted = null;

	front.addEventListener('mousedown', function (ev) {
		dragStarted = {x: ev.screenX, y: ev.screenY};
	});
	document.addEventListener('mouseup', function () {
		dragStarted = null;
	});
	document.addEventListener('mousemove', function (ev) {
		if (!dragStarted) {
			return;
		}

		var dx = ev.screenX - dragStarted.x;
		var dy = ev.screenY - dragStarted.y;

		dragStarted = {x: ev.screenX, y: ev.screenY};

		var multiplier = cam.pos.fov / cProps.h;

		cam.pos.rx = (cam.pos.rx + dy * multiplier) % 360;
		cam.pos.ry = (cam.pos.ry - dx * multiplier) % 360;
	});
	document.addEventListener('keydown', function (ev) {
		switch (ev.code) {
			case 'KeyA':
				cam.accel('x', -1);
				break;
			case 'KeyD':
				cam.accel('x', 1);
				break;
			case 'KeyW':
				cam.accel('z', 1);
				break;
			case 'KeyS':
				cam.accel('z', -1);
				break;
			case 'KeyZ':
				cam.rotate('z', 1);
				break;
			case 'KeyX':
				cam.rotate('z', -1);
				break;
			case 'KeyE':
				cam.accel('y', -1);
				break;
			case 'KeyC':
				cam.accel('y', 1);
				break;
		}
	});

	document.addEventListener('keyup', function (ev) {
		switch (ev.code) {
			case 'KeyZ':
			case 'KeyX':
				cam.rotate('z', 0);
				break;
			case 'KeyW':
			case 'KeyS':
				cam.accel('z', 0);
				break;
			case 'KeyE':
			case 'KeyC':
				cam.accel('y', 0);
				break;
			case 'KeyA':
			case 'KeyD':
				cam.accel('x', 0);
				break;
		}
	});

	front.addEventListener('wheel', function (ev) {
		cam.pos.fov = Math.min(90, Math.max(0.1, cam.pos.fov + cam.pos.fov * ev.deltaY * .0002));
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

			var fovmul = 45 / cam.pos.fov;
			var mul = cProps.h / 45;

			this.x = (this.x + (cam.pos.ry + cam.rotation.ry) * mul);
			this.y = (this.y - (cam.pos.rx + cam.rotation.rx) * mul);
			this.r = (this.r + (cam.pos.rz + cam.rotation.rz) * Math.PI / 180) % (2 * Math.PI);

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
}

var starTypes = [];

function prepareStarTypes () {
	starTypes = [];
	var img = MyCache.getImage('img/stars/small.png');
	if (!img.complete) {
		return;
	}

	var sss = [
		{type: 'image', tx: 256, ty: 128, tw: 128, th: 128, scale: 3}, // Class G
		{type: 'image', tx: 128, ty: 128, tw: 128, th: 128, scale: 2}, // Red giant
		{type: 'image', tx: 0, ty: 0, tw: 128, th: 128, scale: 4}, // White dwarf
		{type: 'image', tx: 256, ty: 0, tw: 128, th: 128, scale: 4}, // Yellow dwarf
		{type: 'image', tx: 0, ty: 128, tw: 128, th: 128, scale: 2}, // Blue giant
		{type: 'image', tx: 0, ty: 0, tw: 128, th: 128, scale: .2}, // Galaxy core
		{type: 'render', fs: ['rgba(68, 85, 51, 255)','rgba(68, 85, 51, 0)'], size: 2, gsize: 4, scale: 1}
	];

	for (var i = 0; i < sss.length; i++) {
		var context;
		var obj = {
			canvas: document.createElement('canvas'),
			scale: sss[i].scale
		};

		if (sss[i].type == 'image') {
			obj.canvas.width = sss[i].tw;
			obj.canvas.height = sss[i].th;
			context = obj.canvas.getContext('2d');
			context.setTransform(1, 0, 0, 1, 0, 0);
			context.drawImage(img,
				sss[i].tx, sss[i].ty, sss[i].tw, sss[i].th,
				0, 0,  sss[i].tw, sss[i].th);
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

var startTime = Date.now();

function main(data) {
	// if (typeof data.children != 'undefined')
	// {
	// 	stars = [{
	// 		pos: {
	// 			x: 0,
	// 			y: 0,
	// 			z: 0
	// 		},
	// 		cls: 5
	// 	}];
	// 	for (var i = data.children.length - 1; i >= 0; --i) {
	// 		stars.push({
	// 			pos: {
	// 				x: data.children[i].coordinates.x,
	// 				y: data.children[i].coordinates.y,
	// 				z: data.children[i].coordinates.z
	// 			},
	// 			cls: data.children[i].genus - 1
	// 		});
	// 	}
	// }
	let galaxyTime = timeZero + (Date.now() - startTime) * timeScale;
	stars = ellipseGalaxy2(stars, galaxyTime);

	if (typeof module !== 'undefined' && module.exports) {
		process.stdout.write('{"stars":[\n');
		for (let i = stars.length - 1; i >= 0; --i) {
			process.stdout.write(JSON.stringify({id: i, cls: stars[i].cls, pos: stars[i].pos}));
			process.stdout.write(i > 0 ? ',\n' : '\n');
		}
		process.stdout.write(']}\n');
		return;
	}

	var dt = Date.now() - beginTime;
	beginTime = Date.now();

	// расчет FPS
	frameCounter = frameCounter + 1;
	if (frameCounter > 49) {
		FPSCounter = frameCounter * 1000 / frameTimer;
		frameCounter = 0;
		frameTimer = 0;
	}
	frameTimer = frameTimer + dt;
	pathTimer = pathTimer + dt;

	if (starTypes.length) {
		background.update();
		cam.update();
		fContext.setTransform(1, 0, 0, 1, 0, 0);
		fContext.clearRect(0, 0, front.width, front.height);
		fContext.fillStyle = '#999';
		fContext.fillText("FPS: " + FPSCounter.toFixed(1), 0, 10);
		fContext.fillText("FOV, DEP: " + cam.pos.fov + ', ' + cam.depth, 0, 20);
		fContext.fillText("TIME: " + galaxyTime, 0, 30);
		if (orbits) drawOrbits(fContext, cam);
		draw(fContext, cam);
	} else {
		prepareStarTypes();
	}

	requestAnimationFrame(main);
}

//cam.init(0, 3000, -1000, -80, 0, -15, 45);
cam.init(0, -20000, -8000, 67.5, 15, 6, 3);

//jQuery.getJSON('stars.json').done(main);
main();

if (!(typeof module !== 'undefined' && module.exports)) {
	document.getElementById('exportBtn').addEventListener('click', function () {
		let i, cStars = {stars: []};
		for (i = 0; i < stars.length; ++i) {
			cStars.stars.push({id: i, dzc: 0, cls: stars[i].cls, pos: stars[i].pos});
		}

		let w = window.open();

		w.document.write(JSON.stringify(cStars));
		w.document.close();
	});
}
