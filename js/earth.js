window.gl = null;
World = {
	gl : null,
	Canvas : null,
	idCounter : 0,
	BASE_LEVEL : 6,
	EARTH_RADIUS : 6378137,
	MAX_LATITUDE : 85.05112877980659,
	MAX_PROJECTED_COORD : 20037508.3427892,
	CURRENT_LEVEL : -1,
	OLD_POSITION : null,
	BELOW_LEVEL : 1
};
World.Enum = {
	FULL_IN : "FULL_IN",
	FULL_OUT : "FULL_OUT",
	IN_OUT : "IN_OUT",
	NOKIA_TILED_MAP : "NOKIA_TILED_MAP",
	Google_TILED_MAP : "Google_TILED_MAP",
	OSM_TILED_MAP : "OSM_TILED_MAP",
	BLENDED_TILED_MAP : "BLENDED_TILED_MAP"
};
World.Math = {
	ONE_RADIAN_EQUAL_DEGREE : 57.29577951308232,
	ONE_DEGREE_EQUAL_RADIAN : 0.017453292519943295,
	LEFT_TOP : "LEFT_TOP",
	RIGHT_TOP : "RIGHT_TOP",
	LEFT_BOTTOM : "LEFT_BOTTOM",
	RIGHT_BOTTOM : "RIGHT_BOTTOM"
};
World.Math.isZero = function (value) {
	if (Math.abs(value) < 0.000001) {
		return true
	} else {
		return false
	}
};
World.Math.numerationSystemTo10 = function (numSys, strNum) {
	var sum = 0;
	for (var i = 0; i < strNum.length; i++) {
		var level = strNum.length - 1 - i;
		var key = parseInt(strNum[i]);
		sum += key * Math.pow(numSys, level)
	}
	return sum
};
World.Math.numerationSystemFrom10 = function (numSys, num) {
	var tempResultArray = [];
	var quotient = Math.floor(num / numSys);
	var remainder = num % numSys;
	tempResultArray.push(remainder);
	while (quotient != 0) {
		num = quotient;
		quotient = Math.floor(num / numSys);
		remainder = num % numSys;
		tempResultArray.push(remainder)
	}
	tempResultArray.reverse();
	var strResult = tempResultArray.join("");
	return strResult
};
World.Math.numerationSystemChange = function (numSysFrom, numSysTo, strNumFrom) {
	var temp10 = this.numerationSystemTo10(numSysFrom, strNumFrom);
	var strResult = this.numerationSystemFrom10(numSysTo, temp10);
	return strResult
};
World.Math.getTriangleArea = function (v1, v2, v3) {
	var v1Copy = v1.getCopy();
	var v2Copy = v2.getCopy();
	var v3Copy = v3.getCopy();
	var direction = v3Copy.minus(v2Copy);
	var line = new World.Line(v2Copy, direction);
	var h = this.getLengthFromVerticeToLine(v1Copy, line);
	var w = this.getLengthFromVerticeToVertice(v2Copy, v3Copy);
	var area = 0.5 * w * h;
	return area
};
World.Math.getLengthFromVerticeToVertice = function (vertice1, vertice2) {
	var vertice1Copy = vertice1.getCopy();
	var vertice2Copy = vertice2.getCopy();
	var length2 = Math.pow(vertice1Copy.x - vertice2Copy.x, 2) + Math.pow(vertice1Copy.y - vertice2Copy.y, 2) + Math.pow(vertice1Copy.z - vertice2Copy.z, 2);
	var length = Math.sqrt(length2);
	return length
};
World.Math.getLengthFromVerticeToLine = function (vertice, line) {
	var verticeCopy = vertice.getCopy();
	var lineCopy = line.getCopy();
	var x0 = verticeCopy.x;
	var y0 = verticeCopy.y;
	var z0 = verticeCopy.z;
	var verticeOnLine = lineCopy.vertice;
	var x1 = verticeOnLine.x;
	var y1 = verticeOnLine.y;
	var z1 = verticeOnLine.z;
	var lineVector = lineCopy.vector;
	lineVector.normalize();
	var a = lineVector.x;
	var b = lineVector.y;
	var c = lineVector.z;
	var A = (y0 - y1) * c - b * (z0 - z1);
	var B = (z0 - z1) * a - c * (x0 - x1);
	var C = (x0 - x1) * b - a * (y0 - y1);
	var result = Math.sqrt(A * A + B * B + C * C);
	return result
};
World.Math.getLengthFromVerticeToPlan = function (vertice, plan) {
	var verticeCopy = vertice.getCopy();
	var planCopy = plan.getCopy();
	var x = verticeCopy.x;
	var y = verticeCopy.y;
	var z = verticeCopy.z;
	var A = planCopy.A;
	var B = planCopy.B;
	var C = planCopy.C;
	var D = planCopy.D;
	var numerator = Math.abs(A * x + B * y + C * z + D);
	var denominator = Math.sqrt(A * A + B * B + C * C);
	var length = numerator / denominator;
	return length
};
World.Math.getVerticeVerticalIntersectPointWidthPlan = function (vertice, plan) {
	var verticeCopy = vertice.getCopy();
	var planCopy = plan.getCopy();
	var x0 = verticeCopy.x;
	var y0 = verticeCopy.y;
	var z0 = verticeCopy.z;
	var normalVector = new World.Vector(planCopy.A, planCopy.B, planCopy.C);
	normalVector.normalize();
	var a = normalVector.x;
	var b = normalVector.y;
	var c = normalVector.z;
	var d = planCopy.D * a / planCopy.A;
	var k =  - (a * x0 + b * y0 + c * z0 + d);
	var x = k * a + x0;
	var y = k * b + y0;
	var z = k * c + z0;
	var intersectVertice = new World.Vertice(x, y, z);
	return intersectVertice
};
World.Math.getIntersectPointByLineAdPlan = function (line, plan) {
	var lineCopy = line.getCopy();
	var planCopy = plan.getCopy();
	lineCopy.vector.normalize();
	var A = planCopy.A;
	var B = planCopy.B;
	var C = planCopy.C;
	var D = planCopy.D;
	var x0 = lineCopy.vertice.x;
	var y0 = lineCopy.vertice.y;
	var z0 = lineCopy.vertice.z;
	var a = lineCopy.vector.x;
	var b = lineCopy.vector.y;
	var c = lineCopy.vector.z;
	var k =  - (A * x0 + B * y0 + C * z0 + D) / (A * a + B * b + C * c);
	var x = k * a + x0;
	var y = k * b + y0;
	var z = k * c + z0;
	var intersectVertice = new World.Vertice(x, y, z);
	return intersectVertice
};
World.Math.getLineIntersectPointWidthEarth = function (line) {
	var lineCopy = line.getCopy();
	vertice = lineCopy.vertice;
	var direction = lineCopy.vector;
	direction.normalize();
	var r = World.EARTH_RADIUS;
	result = [];
	direction.normalize();
	var a = direction.x;
	var b = direction.y;
	var c = direction.z;
	var x0 = vertice.x;
	var y0 = vertice.y;
	var z0 = vertice.z;
	var a2 = a * a;
	var b2 = b * b;
	var c2 = c * c;
	var r2 = r * r;
	var ay0 = a * y0;
	var az0 = a * z0;
	var bx0 = b * x0;
	var bz0 = b * z0;
	var cx0 = c * x0;
	var cy0 = c * y0;
	var deltaA = ay0 * bx0 + az0 * cx0 + bz0 * cy0;
	var deltaB = ay0 * ay0 + az0 * az0 + bx0 * bx0 + bz0 * bz0 + cx0 * cx0 + cy0 * cy0;
	var deltaC = a2 + b2 + c2;
	var delta = 8 * deltaA - 4 * deltaB + 4 * r2 * deltaC;
	var t = a * x0 + b * y0 + c * z0;
	var A = a2 + b2 + c2;
	if (delta < 0) {
		result = []
	} else if (delta == 0) {
		var k = -t / A;
		var x = k * a + x0;
		var y = k * b + y0;
		var z = k * c + z0;
		var p = new World.Vertice(x, y, z);
		result.push(p)
	} else if (delta > 0) {
		var sqrtDelta = Math.sqrt(delta);
		var k1 = (-2 * t + sqrtDelta) / (2 * A);
		var x1 = k1 * a + x0;
		var y1 = k1 * b + y0;
		var z1 = k1 * c + z0;
		var p1 = new World.Vertice(x1, y1, z1);
		result.push(p1);
		var k2 = (-2 * t - sqrtDelta) / (2 * A);
		var x2 = k2 * a + x0;
		var y2 = k2 * b + y0;
		var z2 = k2 * c + z0;
		var p2 = new World.Vertice(x2, y2, z2);
		result.push(p2)
	}
	return result
};
World.Math.getCrossPlaneByLine = function (vertice, direction) {
	var verticeCopy = vertice.getCopy();
	var directionCopy = direction.getCopy();
	directionCopy.normalize();
	var a = directionCopy.x;
	var b = directionCopy.y;
	var c = directionCopy.z;
	var x0 = verticeCopy.x;
	var y0 = verticeCopy.y;
	var z0 = verticeCopy.z;
	var d =  - (a * x0 + b * y0 + c * z0);
	var plan = new World.Plan(a, b, c, d);
	return plan
};
World.Math.isVerticeInSquare = function (vertice, squareVertices) {
	var verticeCopy = vertice.getCopy();
	var squareVerticesCopy = [];
	for (var i = 0; i < squareVertices.length; i++) {
		squareVerticesCopy[i] = squareVertices[i].getCopy()
	}
	var isIn = false;
	var s = 0.0;
	for (var i = 0; i < squareVerticesCopy.length; i++) {
		var v1 = squareVerticesCopy[i];
		var v2 = squareVerticesCopy[(i + 1) % squareVerticesCopy.length];
		s += this.getTriangleArea(verticeCopy, v1, v2)
	}
	var p0 = squareVerticesCopy[0];
	var p1 = squareVerticesCopy[1];
	var p2 = squareVerticesCopy[2];
	var p3 = squareVerticesCopy[3];
	var area = this.getTriangleArea(p0, p1, p2) + this.getTriangleArea(p0, p3, p2);
	var radio = s / area;
	if (radio >= 0.99 && radio <= 1.01) {
		isIn = true
	} else {
		isIn = false
	}
	return isIn
};
World.Math.isVectorInSquarePyramidVectors = function (vertice, vector, squarePyramidVectors) {
	var verticeCopy = vertice.getCopy();
	var vectorCopy = vector.getCopy();
	var squarePyramidVectorsCopy = [];
	for (var i = 0; i < squarePyramidVectors.length; i++) {
		squarePyramidVectorsCopy[i] = squarePyramidVectors[i].getCopy()
	}
	var length = 10000000;
	vectorCopy.normalize();
	var centerVector = new World.Vector(0, 0, 0);
	for (var i = 0; i < squarePyramidVectorsCopy.length; i++) {
		var vectori = squarePyramidVectorsCopy[i];
		vectori.normalize();
		centerVector.x += vectori.x;
		centerVector.y += vectori.y;
		centerVector.z += vectori.z
	}
	centerVector.normalize();
	centerVector.setLength(length);
	var endVertice = verticeCopy.plus(centerVector);
	var plan = this.getCrossPlaneByLine(endVertice, centerVector);
	var squareVertices = [];
	for (var i = 0; i < squarePyramidVectorsCopy.length; i++) {
		var vectori = squarePyramidVectorsCopy[i];
		var linei = new World.Line(verticeCopy, vectori);
		var intersectVertice = this.getIntersectPointByLineAdPlan(linei, plan);
		squareVertices.push(intersectVertice)
	}
	var isIn = this.isVerticeInSquare(endVertice, squareVertices);
	return isIn
};
World.Math.is2SquarePyramidsOverlap = function (vertice, squarePyramidVectors1, squarePyramidVectors2) {
	var verticeCopy = vertice.getCopy();
	var squarePyramidVectors1Copy = [];
	var squarePyramidVectors2Copy = [];
	for (var i = 0; i < squarePyramidVectors1.length; i++) {
		squarePyramidVectors1Copy[i] = squarePyramidVectors1[i].getCopy()
	}
	for (var i = 0; i < squarePyramidVectors2.length; i++) {
		squarePyramidVectors2Copy[i] = squarePyramidVectors2[i].getCopy()
	}
	var isOverlap = false;
	for (var i = 0; i < squarePyramidVectors1Copy.length; i++) {
		var vectori = squarePyramidVectors1Copy[i];
		isOverlap = this.isVectorInSquarePyramidVectors(verticeCopy, vectori, squarePyramidVectors2Copy);
		if (isOverlap) {
			isOverlap = true;
			return isOverlap
		}
	}
	for (var i = 0; i < squarePyramidVectors2Copy.length; i++) {
		var vectori = squarePyramidVectors2Copy[i];
		isOverlap = this.isVectorInSquarePyramidVectors(verticeCopy, vectori, squarePyramidVectors1Copy);
		if (isOverlap) {
			isOverlap = true;
			return isOverlap
		}
	}
	return isOverlap
};
World.Math.getPositionRelativeToCanvasCenter = function (absoluteX, absoluteY) {
	var w = World.Canvas.width;
	var h = World.Canvas.height;
	var relativeX = absoluteX - w / 2;
	var relativeY = h / 2 - absoluteY;
	var result = [relativeX, relativeY];
	return result
};
World.Math.getLengthFromCamera2EarthSurface = function (level) {
	var length;
	if (level == 0) {
		length = 11110948
	} else {
		length = 7820683 / Math.pow(2, level)
	}
	return length
};
World.Math.geographicToCartesianCoord = function (lon, lat) {
	var radianLon = this.degreeToRadian(lon);
	var radianLat = this.degreeToRadian(lat);
	var sin1 = Math.sin(radianLon);
	var cos1 = Math.cos(radianLon);
	var sin2 = Math.sin(radianLat);
	var cos2 = Math.cos(radianLat);
	var x = World.EARTH_RADIUS * sin1 * cos2;
	var y = World.EARTH_RADIUS * sin2;
	var z = World.EARTH_RADIUS * cos1 * cos2;
	var p = new World.Vertice(x, y, z);
	return p
};
World.Math.cartesianCoordToGeographic = function (vertice) {
	var verticeCopy = vertice.getCopy();
	var x = verticeCopy.x;
	var y = verticeCopy.y;
	var z = verticeCopy.z;
	var sin2 = y / World.EARTH_RADIUS;
	if (sin2 > 1) {
		sin2 = 2
	} else if (sin2 < -1) {
		sin2 = -1
	}
	var radianLat = Math.asin(sin2);
	var cos2 = Math.cos(radianLat);
	var sin1 = x / (World.EARTH_RADIUS * cos2);
	if (sin1 > 1) {
		sin1 = 1
	} else if (sin1 < -1) {
		sin1 = -1
	}
	var cos1 = z / (World.EARTH_RADIUS * cos2);
	if (cos1 > 1) {
		cos1 = 1
	} else if (cos1 < -1) {
		cos1 = -1
	}
	var radianLog = Math.asin(sin1);
	if (sin1 >= 0) {
		if (cos1 >= 0) {
			radianLog = radianLog
		} else {
			radianLog = Math.PI - radianLog
		}
	} else {
		if (cos1 >= 0) {
			radianLog = radianLog
		} else {
			radianLog = -radianLog - Math.PI
		}
	}
	var degreeLat = World.Math.radianToDegree(radianLat);
	var degreeLog = World.Math.radianToDegree(radianLog);
	var result = [degreeLog, degreeLat];
	return result
};
World.Math.getTileGridByParent = function (parentLevel, parentRow, parentColumn, position) {
	var level = parentLevel + 1;
	var row = -1;
	var column = -1;
	if (position == this.LEFT_TOP) {
		row = 2 * parentRow;
		column = 2 * parentColumn
	} else if (position == this.RIGHT_TOP) {
		row = 2 * parentRow;
		column = 2 * parentColumn + 1
	} else if (position == this.LEFT_BOTTOM) {
		row = 2 * parentRow + 1;
		column = 2 * parentColumn
	} else if (position == this.RIGHT_BOTTOM) {
		row = 2 * parentRow + 1;
		column = 2 * parentColumn + 1
	}
	var info = {
		level : level,
		row : row,
		column : column
	};
	return info
};
World.Math.getTileGridByGeo = function (lon, lat, level) {
	var info = {
		level : level,
		row : null,
		column : null
	};
	var coordWebMercator = this.degreeGeographicToWebMercator(lon, lat);
	var x = coordWebMercator[0];
	var y = coordWebMercator[1];
	var horX = x + World.MAX_PROJECTED_COORD;
	var verY = World.MAX_PROJECTED_COORD - y;
	var size = World.MAX_PROJECTED_COORD / Math.pow(2, level - 1);
	info.row = Math.floor(verY / size);
	info.column = Math.floor(horX / size);
	return info
};
World.Math.degreeToRadian = function (degree) {
	var radian = degree * this.ONE_DEGREE_EQUAL_RADIAN;
	return radian
};
World.Math.radianToDegree = function (radian) {
	var degree = radian * this.ONE_RADIAN_EQUAL_DEGREE;
	return degree
};
World.Math.webMercatorXToRadianLog = function (x) {
	var radianLog = x / World.EARTH_RADIUS;
	return radianLog
};
World.Math.webMercatorXToDegreeLog = function (x) {
	var radianLog = this.webMercatorXToRadianLog(x);
	var degreeLog = this.radianToDegree(radianLog);
	return degreeLog
};
World.Math.webMercatorYToRadianLat = function (y) {
	var a = y / World.EARTH_RADIUS;
	var b = Math.pow(Math.E, a);
	var c = Math.atan(b);
	var radianLat = 2 * c - Math.PI / 2;
	return radianLat
};
World.Math.webMercatorYToDegreeLat = function (y) {
	var radianLat = this.webMercatorYToRadianLat(y);
	var degreeLat = this.radianToDegree(radianLat);
	return degreeLat
};
World.Math.webMercatorToRadianGeographic = function (x, y) {
	var radianLog = this.webMercatorXToRadianLog(x);
	var radianLat = this.webMercatorYToRadianLat(y);
	var radianResult = [radianLog, radianLat];
	return radianResult
};
World.Math.webMercatorToDegreeGeographic = function (x, y) {
	var degreeLog = this.webMercatorXToDegreeLog(x);
	var degreeLat = this.webMercatorYToDegreeLat(y);
	var degreeResult = [degreeLog, degreeLat];
	return degreeResult
};
World.Math.radianLogToWebMercatorX = function (radianLog) {
	var x = World.EARTH_RADIUS * radianLog;
	return x
};
World.Math.degreeLogToWebMercatorX = function (degreeLog) {
	var radianLog = this.degreeToRadian(degreeLog);
	var x = this.radianLogToWebMercatorX(radianLog);
	return x
};
World.Math.radianLatToWebMercatorY = function (radianLat) {
	var a = Math.PI / 4 + radianLat / 2;
	var b = Math.tan(a);
	var c = Math.log(b);
	var y = World.EARTH_RADIUS * c;
	return y
};
World.Math.degreeLatToWebMercatorY = function (degreeLat) {
	var radianLat = this.degreeToRadian(degreeLat);
	var y = this.radianLatToWebMercatorY(radianLat);
	return y
};
World.Math.radianGeographicToWebMercator = function (radianLog, radianLat) {
	var x = this.radianLogToWebMercatorX(radianLog);
	var y = this.radianLatToWebMercatorY(radianLat);
	var result = [x, y];
	return result
};
World.Math.degreeGeographicToWebMercator = function (degreeLog, degreeLat) {
	var x = this.degreeLogToWebMercatorX(degreeLog);
	var y = this.degreeLatToWebMercatorY(degreeLat);
	var result = [x, y];
	return result
};
World.Math.getTileWebMercatorEnvelopeByGrid = function (level, row, column) {
	var k = World.MAX_PROJECTED_COORD;
	var size = 2 * k / Math.pow(2, level);
	var minX = -k + column * size;
	var maxX = minX + size;
	var maxY = k - row * size;
	var minY = maxY - size;
	var Eproj = {
		"minX" : minX,
		"minY" : minY,
		"maxX" : maxX,
		"maxY" : maxY
	};
	return Eproj
};
World.Math.getTileGeographicEnvelopByGrid = function (level, row, column) {
	var Eproj = this.getTileWebMercatorEnvelopeByGrid(level, row, column);
	var pMin = this.webMercatorToDegreeGeographic(Eproj.minX, Eproj.minY);
	var pMax = this.webMercatorToDegreeGeographic(Eproj.maxX, Eproj.maxY);
	var Egeo = {
		"minLon" : pMin[0],
		"minLat" : pMin[1],
		"maxLon" : pMax[0],
		"maxLat" : pMax[1]
	};
	return Egeo
};
World.Math.getTileCartesianEnvelopByGrid = function (level, row, column) {
	var Egeo = this.getTileGeographicEnvelopByGrid(level, row, column);
	var minLon = Egeo.minLon;
	var minLat = Egeo.minLat;
	var maxLon = Egeo.maxLon;
	var maxLat = Egeo.maxLat;
	var pLeftBottom = this.geographicToCartesianCoord(minLon, minLat);
	var pLeftTop = this.geographicToCartesianCoord(minLon, maxLat);
	var pRightTop = this.geographicToCartesianCoord(maxLon, maxLat);
	var pRightBottom = this.geographicToCartesianCoord(maxLon, minLat);
	var Ecar = {
		"pLeftBottom" : pLeftBottom,
		"pLeftTop" : pLeftTop,
		"pRightTop" : pRightTop,
		"pRightBottom" : pRightBottom,
		"minLon" : minLon,
		"minLat" : minLat,
		"maxLon" : maxLon,
		"maxLat" : maxLat
	};
	return Ecar
};
World.Math.getGeographicTileCenter = function (level, row, column) {
	var lonlatTileCenter;
	var Egeo = this.getTileGeographicEnvelopByGrid(level, row, column);
	var minLon = Egeo.minLon;
	var minLat = Egeo.minLat;
	var maxLon = Egeo.maxLon;
	var maxLat = Egeo.maxLat;
	var centerLon = (minLon + maxLon) / 2;
	var centerLat = (minLat + maxLat) / 2;
	var lonlatTileCenter = [centerLon, centerLat];
	return lonlatTileCenter
};
World.Math.getCartesianTileCenter = function (level, row, column) {
	var lonLat = this.getGeographicTileCenter(level, row, column);
	var vertice = this.geographicToCartesianCoord(lonLat[0], lonLat[1]);
	return vertice
};
World.Math.calculateNormals = function (vs, ind) {
	var x = 0;
	var y = 1;
	var z = 2;
	var ns = [];
	for (var i = 0; i < vs.length; i = i + 3) {
		ns[i + x] = 0.0;
		ns[i + y] = 0.0;
		ns[i + z] = 0.0
	}
	for (var i = 0; i < ind.length; i = i + 3) {
		var v1 = [];
		var v2 = [];
		var normal = [];
		v1[x] = vs[3 * ind[i + 2] + x] - vs[3 * ind[i + 1] + x];
		v1[y] = vs[3 * ind[i + 2] + y] - vs[3 * ind[i + 1] + y];
		v1[z] = vs[3 * ind[i + 2] + z] - vs[3 * ind[i + 1] + z];
		v2[x] = vs[3 * ind[i] + x] - vs[3 * ind[i + 1] + x];
		v2[y] = vs[3 * ind[i] + y] - vs[3 * ind[i + 1] + y];
		v2[z] = vs[3 * ind[i] + z] - vs[3 * ind[i + 1] + z];
		normal[x] = v1[y] * v2[z] - v1[z] * v2[y];
		normal[y] = v1[z] * v2[x] - v1[x] * v2[z];
		normal[z] = v1[x] * v2[y] - v1[y] * v2[x];
		for (j = 0; j < 3; j++) {
			ns[3 * ind[i + j] + x] = ns[3 * ind[i + j] + x] + normal[x];
			ns[3 * ind[i + j] + y] = ns[3 * ind[i + j] + y] + normal[y];
			ns[3 * ind[i + j] + z] = ns[3 * ind[i + j] + z] + normal[z]
		}
	}
	for (var i = 0; i < vs.length; i = i + 3) {
		var nn = [];
		nn[x] = ns[i + x];
		nn[y] = ns[i + y];
		nn[z] = ns[i + z];
		var len = Math.sqrt((nn[x] * nn[x]) + (nn[y] * nn[y]) + (nn[z] * nn[z]));
		if (len == 0)
			len = 1.0;
		nn[x] = nn[x] / len;
		nn[y] = nn[y] / len;
		nn[z] = nn[z] / len;
		ns[i + x] = nn[x];
		ns[i + y] = nn[y];
		ns[i + z] = nn[z]
	}
	return ns
};
window.requestAnimationFrame = window.requestAnimationFrame || window.mozRequestAnimationFrame || window.webkitRequestAnimationFrame || window.msRequestAnimationFrame || window.oRequestAnimationFrame || function (callback) {
	setTimeout(callback, 1000 / 60)
};
World.ShaderContent = {
	SIMPLE_SHADER : {
		VS_CONTENT : "attribute vec3 aVertexPosition;\nattribute vec2 aTextureCoord;\nvarying vec2 vTextureCoord;\nuniform mat4 uMVMatrix;\nuniform mat4 uPMatrix;\nvoid main()\n{\ngl_Position = uPMatrix * uMVMatrix * vec4(aVertexPosition,1.0);\nvTextureCoord = aTextureCoord;\n}",
		FS_CONTENT : "#ifdef GL_ES\nprecision highp float;\n#endif\nuniform bool uUseTexture;\nuniform float uShininess;\nuniform vec3 uLightDirection;\nuniform vec4 uLightAmbient;\nuniform vec4 uLightDiffuse;\nuniform vec4 uLightSpecular;\nvarying vec2 vTextureCoord;\nuniform sampler2D uSampler;\nvoid main()\n{\ngl_FragColor = texture2D(uSampler, vec2(vTextureCoord.s, vTextureCoord.t));\n}"
	}
};
World.WebGLRenderer = function (canvas, vertexShaderText, fragmentShaderText) {
	window.renderer = this;
	this.scene = null;
	this.camera = null;
	this.bAutoRefresh = false;
	function initWebGL(canvas) {
		try {
			var contextList = ["webgl", "experimental-webgl", "webkit-3d", "moz-webgl"];
			for (var i = 0; i < contextList.length; i++) {
				var g = canvas.getContext(contextList[i], {
						antialias : true
					});
				if (g) {
					window.gl = g;
					World.gl = g;
					World.Canvas = canvas;
					break
				}
			}
		} catch (e) {}

	}
	function getShader(gl, shaderType, shaderText) {
		if (!shaderText) {
			return null
		}
		var shader = null;
		if (shaderType == "VERTEX_SHADER") {
			shader = gl.createShader(gl.VERTEX_SHADER)
		} else if (shaderType == "FRAGMENT_SHADER") {
			shader = gl.createShader(gl.FRAGMENT_SHADER)
		} else {
			return null
		}
		gl.shaderSource(shader, shaderText);
		gl.compileShader(shader);
		if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) {
			alert(gl.getShaderInfoLog(shader));
			console.error(gl.getShaderInfoLog(shader));
			gl.deleteShader(shader);
			return null
		}
		return shader
	}
	function initShaders(vertexShaderText, fragmentShaderText) {
		var vertexShader = getShader(World.gl, "VERTEX_SHADER", vertexShaderText);
		var fragmentShader = getShader(World.gl, "FRAGMENT_SHADER", fragmentShaderText);
		var shaderProgram = gl.createProgram();
		gl.attachShader(shaderProgram, vertexShader);
		gl.attachShader(shaderProgram, fragmentShader);
		gl.linkProgram(shaderProgram);
		if (!gl.getProgramParameter(shaderProgram, gl.LINK_STATUS)) {
			console.error("Could not link program!");
			gl.deleteProgram(shaderProgram);
			gl.deleteShader(vertexShader);
			gl.deleteShader(fragmentShader);
			return
		}
		gl.useProgram(shaderProgram);
		gl.shaderProgram = shaderProgram;
		gl.shaderProgram.aVertexPosition = gl.getAttribLocation(gl.shaderProgram, "aVertexPosition");
		gl.shaderProgram.aTextureCoord = gl.getAttribLocation(gl.shaderProgram, "aTextureCoord");
		gl.shaderProgram.uMVMatrix = gl.getUniformLocation(gl.shaderProgram, "uMVMatrix");
		gl.shaderProgram.uPMatrix = gl.getUniformLocation(gl.shaderProgram, "uPMatrix");
		gl.shaderProgram.uSampler = gl.getUniformLocation(gl.shaderProgram, "uSampler");
		var squareArray = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
		var squareMatrix = new Float32Array(squareArray);
		gl.uniformMatrix4fv(gl.shaderProgram.uMVMatrix, false, squareMatrix)
	}
	initWebGL(canvas);
	if (!window.gl) {
		alert("浏览器不支持WebGL或将WebGL禁用!");
		console.debug("浏览器不支持WebGL或将WebGL禁用!");
		return
	}
	initShaders(vertexShaderText, fragmentShaderText);
	gl.clearColor(0.25, 0.25, 0.25, 1.0);
	gl.disable(gl.DEPTH_TEST);
	gl.depthFunc(gl.LEQUAL);
	gl.enable(gl.CULL_FACE);
	gl.frontFace(gl.CCW);
	gl.cullFace(gl.BACK)
};
World.WebGLRenderer.prototype = {
	constructor : World.WebGLRenderer,
	render : function (scene, camera) {
		gl.viewport(0, 0, World.Canvas.width, World.Canvas.height);
		gl.clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT);
		camera.viewMatrix = null;
		camera.viewMatrix = camera.getViewMatrix();
		for (var i = 0; i < scene.objectList.length; i++) {
			var obj = scene.objectList[i];
			if (obj) {
				if (obj.isvisible) {
					obj.draw(camera)
				}
			}
		}
	},
	bindScene : function (scene) {
		this.scene = scene
	},
	bindCamera : function (camera) {
		this.camera = camera
	},
	tick : function () {
		var renderer = window.renderer;
		if (renderer.scene && renderer.camera) {
			renderer.render(renderer.scene, renderer.camera)
		}
		if (renderer.bAutoRefresh) {
			window.requestAnimationFrame(renderer.tick)
		}
	},
	setIfAutoRefresh : function (bAuto) {
		this.bAutoRefresh = bAuto;
		if (this.bAutoRefresh) {
			this.tick()
		}
	}
};
World.Scene = function () {
	this.objectList = []
};
World.Scene.prototype = {
	constructor : World.Scene,
	findObjById : function (objId) {
		for (var i = 0; i < this.objectList.length; i++) {
			var obj = this.objectList[i];
			if (obj.id == objId) {
				obj.index = i;
				return obj
			}
		}
		return null
	},
	findTile : function (level, row, column) {
		for (var i = 0; i < this.objectList.length; i++) {
			var obj = this.objectList[i];
			if (obj instanceof World.Tile) {
				if (obj.level == level && obj.row == row && obj.column == column) {
					return obj
				}
			}
		}
		return null
	},
	add : function (obj) {
		if (obj) {
			if (obj instanceof World.Object3D || obj instanceof World.Object3DComponents) {
				if (this.findObjById(obj.id) != null) {
					console.debug("obj已经存在于Scene中，无法将其再次加入！");
					return
				} else {
					this.objectList.push(obj);
					obj.scene = this
				}
			}
		}
	},
	remove : function (obj) {
		if (obj) {
			var result = this.findObjById(obj.id);
			if (result == null) {
				console.debug("obj不存在于Scene中，所以无法将其从中删除！");
				return false
			}
			obj.destroy();
			this.objectList.splice(result.index, 1);
			obj = null;
			return true
		} else {
			return false
		}
	},
	clear : function () {
		for (var i = 0; i < this.objectList.length; i++) {
			var obj = this.objectList[i];
			obj.destroy()
		}
		this.objectList = []
	}
};
World.Vertice = function (x, y, z) {
	this.x = x || 0;
	this.y = y || 0;
	this.z = z || 0
};
World.Vertice.prototype = {
	constructor : World.Vertice,
	minus : function (otherVertice) {
		var x = this.x - otherVertice.x;
		var y = this.y - otherVertice.y;
		var z = this.z - otherVertice.z;
		return new World.Vector(x, y, z)
	},
	plus : function (otherVector) {
		var x = this.x + otherVector.x;
		var y = this.y + otherVector.y;
		var z = this.z + otherVector.z;
		return new World.Vertice(x, y, z)
	},
	getVector : function () {
		return new World.Vector(this.x, this.y, this.z)
	},
	getArray : function () {
		return [this.x, this.y, this.z]
	},
	getCopy : function () {
		return new World.Vertice(this.x, this.y, this.z)
	},
	getOpposite : function () {
		return new World.Vertice(-this.x, -this.y, -this.z)
	}
};
World.Vector = function (x, y, z) {
	this.x = x || 0;
	this.y = y || 0;
	this.z = z || 0
};
World.Vector.prototype = {
	constructor : World.Vector,
	getVertice : function () {
		return new World.Vertice(this.x, this.y, this.z)
	},
	getArray : function () {
		return [this.x, this.y, this.z]
	},
	getCopy : function () {
		return new World.Vector(this.x, this.y, this.z)
	},
	getOpposite : function () {
		return new World.Vector(-this.x, -this.y, -this.z)
	},
	getLength : function () {
		return Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z)
	},
	normalize : function () {
		var length = this.getLength();
		if (!World.Math.isZero(length)) {
			this.x /= length;
			this.y /= length;
			this.z /= length
		} else {
			this.x = 0;
			this.y = 0;
			this.z = 0
		}
		return this
	},
	setLength : function (length) {
		this.normalize();
		this.x *= length;
		this.y *= length;
		this.z *= length;
		return this
	},
	getRandomVerticalVector : function () {
		var result;
		var length = this.getLength();
		if (length == 0) {
			result = new World.Vector(0, 0, 0)
		} else {
			var x2,
			y2,
			z2;
			if (this.x != 0) {
				y2 = 1;
				z2 = 0;
				x2 = -this.y / this.x
			} else if (this.y != 0) {
				z2 = 1;
				x2 = 0;
				y2 = -this.z / this.y
			} else if (this.z != 0) {
				x2 = 1;
				y2 = 0;
				z2 = -this.x / this.z
			}
			result = new World.Vector(x2, y2, z2);
			result.normalize()
		}
		return result
	},
	cross : function (other) {
		var x = this.y * other.z - this.z * other.y;
		var y = this.z * other.x - this.x * other.z;
		var z = this.x * other.y - this.y * other.x;
		return new World.Vector(x, y, z)
	},
	dot : function (other) {
		var result = this.x * other.x + this.y * other.y + this.z * other.z;
		return result
	},
	rotateVertice : function (vertice, radian) {
		var mat = new World.Matrix();
		mat.worldRotateByVector(radian, this);
		var point = [vertice.x, vertice.y, vertice.z, 1];
		var newPoint = mat.multiplyColumn(point);
		var newVertice = new World.Vertice(newPoint[0], newPoint[1], newPoint[2]);
		return newVertice
	},
	rotateVector : function (other, radian) {
		var vertice = other.getVertice();
		var newVertice = this.rotateVertice(vertice, radian);
		var newVector = newVertice.getVector();
		return newVector
	}
};
World.Ray = function (position, direction) {
	this.vertice = position.getCopy();
	this.vector = direction.getCopy();
	this.vector.normalize()
};
World.Ray.prototype.constructor = World.Ray;
World.Ray.prototype.setVertice = function (position) {
	this.vertice = position.getCopy();
	return this
};
World.Ray.prototype.setVector = function (direction) {
	this.vector = direction.getCopy();
	this.vector.normalize();
	return this
};
World.Ray.prototype.getCopy = function () {
	var rayCopy = new World.Ray(this.vertice, this.vector);
	return rayCopy
};
World.Line = function (position, direction) {
	this.vertice = position.getCopy();
	this.vector = direction.getCopy();
	this.vector.normalize()
};
World.Line.prototype.constructor = World.Line;
World.Line.prototype.setVertice = function (position) {
	this.vertice = position.getCopy();
	return this
};
World.Line.prototype.setVector = function (direction) {
	this.vector = direction.getCopy();
	this.vector.normalize();
	return this
};
World.Line.prototype.getCopy = function () {
	var lineCopy = new World.Line(this.vertice, this.vector);
	return lineCopy
};
World.Plan = function (A, B, C, D) {
	this.A = A;
	this.B = B;
	this.C = C;
	this.D = D
};
World.Plan.prototype.constructor = World.Plan;
World.Plan.prototype.getCopy = function () {
	var planCopy = new World.Plan(this.A, this.B, this.C, this.D);
	return planCopy
};
World.Matrix = function (m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34, m41, m42, m43, m44) {
	this.elements = new Float32Array(16);
	this.setElements((m11 == undefined ? 1 : m11), (m12 == undefined ? 0 : m12), (m13 == undefined ? 0 : m13), (m14 == undefined ? 0 : m14), (m21 == undefined ? 0 : m21), (m22 == undefined ? 1 : m22), (m23 == undefined ? 0 : m23), (m24 == undefined ? 0 : m24), (m31 == undefined ? 0 : m31), (m32 == undefined ? 0 : m32), (m33 == undefined ? 1 : m33), (m34 == undefined ? 0 : m34), (m41 == undefined ? 0 : m41), (m42 == undefined ? 0 : m42), (m43 == undefined ? 0 : m43), (m44 == undefined ? 1 : m44))
};
World.Matrix.prototype = {
	constructor : World.Matrix,
	setElements : function (m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34, m41, m42, m43, m44) {
		var values = this.elements;
		values[0] = m11;
		values[4] = m12;
		values[8] = m13;
		values[12] = m14;
		values[1] = m21;
		values[5] = m22;
		values[9] = m23;
		values[13] = m24;
		values[2] = m31;
		values[6] = m32;
		values[10] = m33;
		values[14] = m34;
		values[3] = m41;
		values[7] = m42;
		values[11] = m43;
		values[15] = m44
	},
	setColumnX : function (x, y, z) {
		this.elements[0] = x;
		this.elements[1] = y;
		this.elements[2] = z
	},
	getColumnX : function () {
		return new World.Vertice(this.elements[0], this.elements[1], this.elements[2])
	},
	setColumnY : function (x, y, z) {
		this.elements[4] = x;
		this.elements[5] = y;
		this.elements[6] = z
	},
	getColumnY : function () {
		return new World.Vertice(this.elements[4], this.elements[5], this.elements[6])
	},
	setColumnZ : function (x, y, z) {
		this.elements[8] = x;
		this.elements[9] = y;
		this.elements[10] = z
	},
	getColumnZ : function () {
		return new World.Vertice(this.elements[8], this.elements[9], this.elements[10])
	},
	setColumnTrans : function (x, y, z) {
		this.elements[12] = x;
		this.elements[13] = y;
		this.elements[14] = z
	},
	getColumnTrans : function () {
		return new World.Vertice(this.elements[12], this.elements[13], this.elements[14])
	},
	setLastRowDefault : function () {
		this.elements[3] = 0;
		this.elements[7] = 0;
		this.elements[11] = 0;
		this.elements[15] = 1
	},
	transpose : function () {
		var result = this.getTransposeMatrix();
		this.setMatrixByOther(result)
	},
	getTransposeMatrix : function () {
		var result = new World.Matrix();
		result.elements[0] = this.elements[0];
		result.elements[4] = this.elements[1];
		result.elements[8] = this.elements[2];
		result.elements[12] = this.elements[3];
		result.elements[1] = this.elements[4];
		result.elements[5] = this.elements[5];
		result.elements[9] = this.elements[6];
		result.elements[13] = this.elements[7];
		result.elements[2] = this.elements[8];
		result.elements[6] = this.elements[9];
		result.elements[10] = this.elements[10];
		result.elements[14] = this.elements[11];
		result.elements[3] = this.elements[12];
		result.elements[7] = this.elements[13];
		result.elements[11] = this.elements[14];
		result.elements[15] = this.elements[15];
		return result
	},
	inverse : function () {
		var result = this.getInverseMatrix();
		this.setMatrixByOther(result)
	},
	getInverseMatrix : function () {
		var a = this.elements;
		var result = new World.Matrix();
		var b = result.elements;
		var c = a[0],
		d = a[1],
		e = a[2],
		g = a[3],
		f = a[4],
		h = a[5],
		i = a[6],
		j = a[7],
		k = a[8],
		l = a[9],
		n = a[10],
		o = a[11],
		m = a[12],
		p = a[13],
		r = a[14],
		s = a[15];
		var A = c * h - d * f;
		var B = c * i - e * f;
		t = c * j - g * f;
		u = d * i - e * h;
		v = d * j - g * h;
		w = e * j - g * i;
		x = k * p - l * m;
		y = k * r - n * m;
		z = k * s - o * m;
		C = l * r - n * p;
		D = l * s - o * p;
		E = n * s - o * r;
		q = A * E - B * D + t * C + u * z - v * y + w * x;
		if (!q)
			return null;
		q = 1 / q;
		b[0] = (h * E - i * D + j * C) * q;
		b[1] = (-d * E + e * D - g * C) * q;
		b[2] = (p * w - r * v + s * u) * q;
		b[3] = (-l * w + n * v - o * u) * q;
		b[4] = (-f * E + i * z - j * y) * q;
		b[5] = (c * E - e * z + g * y) * q;
		b[6] = (-m * w + r * t - s * B) * q;
		b[7] = (k * w - n * t + o * B) * q;
		b[8] = (f * D - h * z + j * x) * q;
		b[9] = (-c * D + d * z - g * x) * q;
		b[10] = (m * v - p * t + s * A) * q;
		b[11] = (-k * v + l * t - o * A) * q;
		b[12] = (-f * C + h * y - i * x) * q;
		b[13] = (c * C - d * y + e * x) * q;
		b[14] = (-m * u + p * B - r * A) * q;
		b[15] = (k * u - l * B + n * A) * q;
		return result
	},
	setMatrixByOther : function (otherMatrix) {
		for (var i = 0; i < otherMatrix.elements.length; i++) {
			this.elements[i] = otherMatrix.elements[i]
		}
	},
	setSquareMatrix : function () {
		this.setElements(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1)
	},
	copy : function () {
		var clone = new World.Matrix(this.elements[0], this.elements[4], this.elements[8], this.elements[12], this.elements[1], this.elements[5], this.elements[9], this.elements[13], this.elements[2], this.elements[6], this.elements[10], this.elements[14], this.elements[3], this.elements[7], this.elements[11], this.elements[15]);
		return clone
	},
	multiplyMatrix : function (otherMatrix) {
		var values1 = this.elements;
		var values2 = otherMatrix.elements;
		var m11 = values1[0] * values2[0] + values1[4] * values2[1] + values1[8] * values2[2] + values1[12] * values2[3];
		var m12 = values1[0] * values2[4] + values1[4] * values2[5] + values1[8] * values2[6] + values1[12] * values2[7];
		var m13 = values1[0] * values2[8] + values1[4] * values2[9] + values1[8] * values2[10] + values1[12] * values2[11];
		var m14 = values1[0] * values2[12] + values1[4] * values2[13] + values1[8] * values2[14] + values1[12] * values2[15];
		var m21 = values1[1] * values2[0] + values1[5] * values2[1] + values1[9] * values2[2] + values1[13] * values2[3];
		var m22 = values1[1] * values2[4] + values1[5] * values2[5] + values1[9] * values2[6] + values1[13] * values2[7];
		var m23 = values1[1] * values2[8] + values1[5] * values2[9] + values1[9] * values2[10] + values1[13] * values2[11];
		var m24 = values1[1] * values2[12] + values1[5] * values2[13] + values1[9] * values2[14] + values1[13] * values2[15];
		var m31 = values1[2] * values2[0] + values1[6] * values2[1] + values1[10] * values2[2] + values1[14] * values2[3];
		var m32 = values1[2] * values2[4] + values1[6] * values2[5] + values1[10] * values2[6] + values1[14] * values2[7];
		var m33 = values1[2] * values2[8] + values1[6] * values2[9] + values1[10] * values2[10] + values1[14] * values2[11];
		var m34 = values1[2] * values2[12] + values1[6] * values2[13] + values1[10] * values2[14] + values1[14] * values2[15];
		var m41 = values1[3] * values2[0] + values1[7] * values2[1] + values1[11] * values2[2] + values1[15] * values2[3];
		var m42 = values1[3] * values2[4] + values1[7] * values2[5] + values1[11] * values2[6] + values1[15] * values2[7];
		var m43 = values1[3] * values2[8] + values1[7] * values2[9] + values1[11] * values2[10] + values1[15] * values2[11];
		var m44 = values1[3] * values2[12] + values1[7] * values2[13] + values1[11] * values2[14] + values1[15] * values2[15];
		return new World.Matrix(m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34, m41, m42, m43, m44)
	},
	multiplyColumn : function (column) {
		var values1 = this.elements;
		var values2 = column;
		var m11 = values1[0] * values2[0] + values1[4] * values2[1] + values1[8] * values2[2] + values1[12] * values2[3];
		var m21 = values1[1] * values2[0] + values1[5] * values2[1] + values1[9] * values2[2] + values1[13] * values2[3];
		var m31 = values1[2] * values2[0] + values1[6] * values2[1] + values1[10] * values2[2] + values1[14] * values2[3];
		var m41 = values1[3] * values2[0] + values1[7] * values2[1] + values1[11] * values2[2] + values1[15] * values2[3];
		var result = [m11, m21, m31, m41];
		return result
	},
	checkZero : function () {
		for (var i = 0; i < this.elements.length; i++) {
			if (World.Math.isZero(this.elements[i])) {
				this.elements[i] = 0
			}
		}
	},
	worldTranslate : function (x, y, z) {
		this.elements[12] += x;
		this.elements[13] += y;
		this.elements[14] += z
	},
	worldRotateX : function (radian) {
		var c = Math.cos(radian);
		var s = Math.sin(radian);
		var m = new World.Matrix(1, 0, 0, 0, 0, c, -s, 0, 0, s, c, 0, 0, 0, 0, 1);
		var result = m.multiplyMatrix(this);
		result.checkZero();
		this.setMatrixByOther(result)
	},
	worldRotateY : function (radian) {
		var c = Math.cos(radian);
		var s = Math.sin(radian);
		var m = new World.Matrix(c, 0, s, 0, 0, 1, 0, 0, -s, 0, c, 0, 0, 0, 0, 1);
		var result = m.multiplyMatrix(this);
		result.checkZero();
		this.setMatrixByOther(result)
	},
	worldRotateZ : function (radian) {
		var c = Math.cos(radian);
		var s = Math.sin(radian);
		var m = new World.Matrix(c, -s, 0, 0, s, c, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1);
		var result = m.multiplyMatrix(this);
		result.checkZero();
		this.setMatrixByOther(result)
	},
	worldRotateByVector : function (radian, vector) {
		var x = vector.x;
		var y = vector.y;
		var z = vector.z;
		var length,
		s,
		c;
		var xx,
		yy,
		zz,
		xy,
		yz,
		zx,
		xs,
		ys,
		zs,
		one_c;
		s = Math.sin(radian);
		c = Math.cos(radian);
		length = Math.sqrt(x * x + y * y + z * z);
		x /= length;
		y /= length;
		z /= length;
		xx = x * x;
		yy = y * y;
		zz = z * z;
		xy = x * y;
		yz = y * z;
		zx = z * x;
		xs = x * s;
		ys = y * s;
		zs = z * s;
		one_c = 1.0 - c;
		var m11 = (one_c * xx) + c;
		var m12 = (one_c * xy) - zs;
		var m13 = (one_c * zx) + ys;
		var m14 = 0.0;
		var m21 = (one_c * xy) + zs;
		var m22 = (one_c * yy) + c;
		var m23 = (one_c * yz) - xs;
		var m24 = 0.0;
		var m31 = (one_c * zx) - ys;
		var m32 = (one_c * yz) + xs;
		var m33 = (one_c * zz) + c;
		var m34 = 0.0;
		var m41 = 0.0;
		var m42 = 0.0;
		var m43 = 0.0;
		var m44 = 1.0;
		var mat = new World.Matrix(m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34, m41, m42, m43, m44);
		var result = mat.multiplyMatrix(this);
		result.checkZero();
		this.setMatrixByOther(result)
	},
	localRotateX : function (radian) {
		var transX = this.elements[12];
		var transY = this.elements[13];
		var transZ = this.elements[14];
		this.worldTranslate(-transX, -transY, -transZ);
		var columnX = this.getColumnX();
		this.worldRotateByVector(radian, columnX);
		this.worldTranslate(transX, transY, transZ)
	},
	localRotateY : function (radian) {
		var transX = this.elements[12];
		var transY = this.elements[13];
		var transZ = this.elements[14];
		this.worldTranslate(-transX, -transY, -transZ);
		var columnY = this.getColumnY();
		this.worldRotateByVector(radian, columnY);
		this.worldTranslate(transX, transY, transZ)
	},
	localRotateZ : function (radian) {
		var transX = this.elements[12];
		var transY = this.elements[13];
		var transZ = this.elements[14];
		this.worldTranslate(-transX, -transY, -transZ);
		var columnZ = this.getColumnZ();
		this.worldRotateByVector(radian, columnZ);
		this.worldTranslate(transX, transY, transZ)
	}
};
World.TextureMaterial = function (url) {
	if (url) {
		this.setImageUrl(url)
	}
};
World.TextureMaterial.prototype.setImageUrl = function (url) {
	function handleLoadedTexture(texture) {
		gl.bindTexture(gl.TEXTURE_2D, texture);
		gl.pixelStorei(gl.UNPACK_FLIP_Y_WEBGL, true);
		gl.texImage2D(gl.TEXTURE_2D, 0, gl.RGBA, gl.RGBA, gl.UNSIGNED_BYTE, texture.image);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR_MIPMAP_NEAREST);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR_MIPMAP_NEAREST);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
		gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
		gl.generateMipmap(gl.TEXTURE_2D);
		gl.bindTexture(gl.TEXTURE_2D, null);
		texture.loaded = true
	}
	this.texture = gl.createTexture();
	this.texture.image = new Image();
	this.texture.image.crossOrigin = 'anonymous';
	this.texture.loaded = false;
	var texture = this.texture;
	this.texture.image.onload = function () {
		handleLoadedTexture(texture)
	};
	this.texture.image.src = url
};
World.Object3D = function (args) {
	this.matrix = new World.Matrix();
	this.id = ++World.idCounter;
	this.scene = null;
	this.vertices = [];
	this.vertexBuffer = null;
	this.indices = [];
	this.indexBuffer = null;
	this.textureCoords = [];
	this.textureCoordBuffer = null;
	this.material = new World.TextureMaterial();
	this.isvisible = true;
	if (args) {
		if (args.material) {
			this.material = args.material
		}
	}
	this.init(args)
};
World.Object3D.prototype = {
	constructor : World.Object3D,
	init : function (params) {},
	initBuffers : function (infos) {
		if (infos) {
			this.vertices = infos.vertices || [];
			this.indices = infos.indices || [];
			this.textureCoords = infos.textureCoords || [];
			if (this.vertices.length > 0 && this.indices.length > 0) {
				this.vertexBuffer = gl.createBuffer();
				gl.bindBuffer(gl.ARRAY_BUFFER, this.vertexBuffer);
				gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(this.vertices), gl.STATIC_DRAW);
				this.indexBuffer = gl.createBuffer();
				gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.indexBuffer);
				gl.bufferData(gl.ELEMENT_ARRAY_BUFFER, new Uint16Array(this.indices), gl.STATIC_DRAW)
			}
			if (this.material instanceof World.TextureMaterial) {
				if (this.textureCoords.length > 0) {
					this.textureCoordBuffer = gl.createBuffer();
					gl.bindBuffer(gl.ARRAY_BUFFER, this.textureCoordBuffer);
					gl.bufferData(gl.ARRAY_BUFFER, new Float32Array(this.textureCoords), gl.STATIC_DRAW)
				}
			}
		}
		gl.bindBuffer(gl.ARRAY_BUFFER, null);
		gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, null)
	},
	setShaderMatrix : function (camera) {
		camera.viewMatrix = (camera.viewMatrix instanceof World.Matrix) ? camera.viewMatrix : camera.getViewMatrix();
		var mvMatrix = camera.viewMatrix.multiplyMatrix(this.matrix);
		gl.uniformMatrix4fv(gl.shaderProgram.uMVMatrix, false, mvMatrix.elements);
		gl.uniformMatrix4fv(gl.shaderProgram.uPMatrix, false, camera.projMatrix.elements)
	},
	draw : function (camera) {
		if (this.material instanceof World.TextureMaterial) {
			if (this.isvisible && this.material.texture.loaded) {
				gl.enableVertexAttribArray(gl.shaderProgram.aTextureCoord);
				gl.bindBuffer(gl.ARRAY_BUFFER, this.textureCoordBuffer);
				gl.vertexAttribPointer(gl.shaderProgram.aTextureCoord, 2, gl.FLOAT, false, 0, 0);
				gl.activeTexture(gl.TEXTURE0);
				gl.bindTexture(gl.TEXTURE_2D, this.material.texture);
				gl.uniform1i(gl.shaderProgram.uSampler, 0);
				this.setShaderMatrix(camera);
				gl.enableVertexAttribArray(gl.shaderProgram.aVertexPosition);
				gl.bindBuffer(gl.ARRAY_BUFFER, this.vertexBuffer);
				gl.vertexAttribPointer(gl.shaderProgram.aVertexPosition, 3, gl.FLOAT, false, 0, 0);
				gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, this.indexBuffer);
				gl.drawElements(gl.TRIANGLES, this.indices.length, gl.UNSIGNED_SHORT, 0);
				gl.bindBuffer(gl.ARRAY_BUFFER, null);
				gl.bindBuffer(gl.ELEMENT_ARRAY_BUFFER, null);
				gl.bindTexture(gl.TEXTURE_2D, null)
			}
		}
	},
	destroy : function () {
		this.scene = null;
		if (this.vertexBuffer) {
			gl.deleteBuffer(this.vertexBuffer)
		}
		if (this.indexBuffer) {
			gl.deleteBuffer(this.indexBuffer)
		}
		if (this.textureCoordBuffer) {
			gl.deleteBuffer(this.textureCoordBuffer)
		}
		if (this.material instanceof World.TextureMaterial) {
			if (this.material.texture) {
				gl.deleteTexture(this.material.texture)
			}
		}
		this.vertexBuffer = null;
		this.indexBuffer = null;
		this.textureCoordBuffer = null;
		if (this.material instanceof World.TextureMaterial) {
			if (this.material.texture) {
				this.material.texture = null
			}
		}
	},
	getPosition : function () {
		var position = this.matrix.getColumnTrans();
		return position
	},
	setPosition : function (x, y, z) {
		this.matrix.setColumnTrans(x, y, z)
	},
	worldTranslate : function (x, y, z) {
		this.matrix.worldTranslate(x, y, z)
	},
	worldRotateX : function (radian) {
		this.matrix.worldRotateX(radian)
	},
	worldRotateY : function (radian) {
		this.matrix.worldRotateY(radian)
	},
	worldRotateZ : function (radian) {
		this.matrix.worldRotateZ(radian)
	},
	worldRotateByVector : function (radian, vector) {
		this.matrix.worldRotateByVector(radian, vector)
	},
	localRotateX : function (radian) {
		this.matrix.localRotateX(radian)
	},
	localRotateY : function (radian) {
		this.matrix.localRotateY(radian)
	},
	localRotateZ : function (radian) {
		this.matrix.localRotateZ(radian)
	},
	getXAxisDirection : function () {
		var columnX = this.matrix.getColumnX();
		var directionX = columnX.getVector();
		directionX.normalize();
		return directionX
	},
	getYAxisDirection : function () {
		var columnY = this.matrix.getColumnY();
		var directionY = columnY.getVector();
		directionY.normalize();
		return directionY
	},
	getZAxisDirection : function () {
		var columnZ = this.matrix.getColumnZ();
		var directionZ = columnZ.getVector();
		directionZ.normalize();
		return directionZ
	},
	lookAt : function (vector) {},
	getEnvelop : function () {}

};
World.Object3DComponents = function (params) {
	this.id = ++World.idCounter;
	this.isvisible = true;
	this.scene = null;
	this.objectList = []
};
World.Object3DComponents.prototype = {
	constructor : World.Object3DComponents,
	draw : function (camera) {
		if (camera instanceof World.PerspectiveCamera) {
			for (var i = 0; i < this.objectList.length; i++) {
				var obj = this.objectList[i];
				if (obj) {
					if (obj.isvisible) {
						obj.draw(camera)
					}
				}
			}
		}
	},
	findObjById : function (objId) {
		for (var i = 0; i < this.objectList.length; i++) {
			var obj = this.objectList[i];
			if (obj.id == objId) {
				obj.index = i;
				return obj
			}
		}
		return null
	},
	add : function (obj) {
		if (obj) {
			if (obj instanceof World.Object3D || obj instanceof World.Object3DComponents) {
				if (this.findObjById(obj.id) != null) {
					console.debug("obj已经存在于Object3DComponents中，无法将其再次加入！");
					return
				} else {
					this.objectList.push(obj);
					obj.scene = this
				}
			}
		}
	},
	destroy : function () {
		this.scene = null;
		this.clear()
	},
	remove : function (obj) {
		if (obj) {
			var result = this.findObjById(obj.id);
			if (result == null) {
				console.debug("obj不存在于Object3DComponents中，所以无法将其从中删除！");
				return false
			}
			obj.destroy();
			this.objectList.splice(result.index, 1);
			obj = null;
			return true
		} else {
			return false
		}
	},
	clear : function () {
		for (var i = 0; i < this.objectList.length; i++) {
			var obj = this.objectList[i];
			obj.destroy()
		}
		this.objectList = []
	}
};
World.Tile = function (args) {
	if (args) {
		this.level = 0;
		this.row = 0;
		this.column = 0;
		this.url = args.url;
		this.subTiledLayer = null;
		World.Object3D.apply(this, arguments)
	}
};
World.Tile.prototype = new World.Object3D();
World.Tile.prototype.constructor = World.Tile;
World.Tile.prototype.init = function (args) {
	if (args) {
		if (args.level != undefined && args.row != undefined && args.column != undefined) {
			this.level = args.level;
			this.row = args.row;
			this.column = args.column
		}
	}
	this.setTileInfoByGrid(this.level, this.row, this.column);
	var minCoord = World.Math.degreeGeographicToWebMercator(this.minLon, this.minLat);
	var maxCoord = World.Math.degreeGeographicToWebMercator(this.maxLon, this.maxLat);
	var minX = minCoord[0] - 1;
	var minY = minCoord[1] - 1;
	var maxX = maxCoord[0] + 1;
	var maxY = maxCoord[1] + 1;
	var vertices = [];
	var indices = [];
	var textureCoords = [];
	var deltaX = (maxX - minX) / this.segment;
	var deltaY = (maxY - minY) / this.segment;
	var deltaTextureCoord = 1.0 / this.segment;
	for (var i = 0; i < this.segment; i++) {
		for (var j = 0; j < this.segment; j++) {
			var xLeft = minX + deltaX * j;
			var xRight = xLeft + deltaX;
			var yBottom = minY + deltaY * i;
			var yTop = yBottom + deltaY;
			var pLB = World.Math.webMercatorToDegreeGeographic(xLeft, yBottom);
			var pRT = World.Math.webMercatorToDegreeGeographic(xRight, yTop);
			var logLeft = pLB[0];
			var logRight = pRT[0];
			var latBottom = pLB[1];
			var latTop = pRT[1];
			var tLeft = deltaTextureCoord * j;
			var tRight = tLeft + deltaTextureCoord;
			var tBottom = deltaTextureCoord * i;
			var tTop = tBottom + deltaTextureCoord;
			var idx = (this.segment * i + j) * 4;
			var idx0 = idx + 0;
			var pLeftTop = World.Math.geographicToCartesianCoord(logLeft, latTop).getArray();
			var tLeftTop = [tLeft, tTop];
			vertices = vertices.concat(pLeftTop);
			textureCoords = textureCoords.concat(tLeftTop);
			var idx1 = idx + 1;
			var pLeftBottom = World.Math.geographicToCartesianCoord(logLeft, latBottom).getArray();
			var tLeftBottom = [tLeft, tBottom];
			vertices = vertices.concat(pLeftBottom);
			textureCoords = textureCoords.concat(tLeftBottom);
			var idx2 = idx + 2;
			var pRightBottom = World.Math.geographicToCartesianCoord(logRight, latBottom).getArray();
			var tRightBottom = [tRight, tBottom];
			vertices = vertices.concat(pRightBottom);
			textureCoords = textureCoords.concat(tRightBottom);
			var idx3 = idx + 3;
			var pRightTop = World.Math.geographicToCartesianCoord(logRight, latTop).getArray();
			var tRightTop = [tRight, tTop];
			vertices = vertices.concat(pRightTop);
			textureCoords = textureCoords.concat(tRightTop);
			indices = indices.concat([idx0, idx1, idx2]);
			indices = indices.concat([idx2, idx3, idx0])
		}
	}
	var infos = {
		vertices : vertices,
		indices : indices,
		textureCoords : textureCoords
	};
	this.initBuffers(infos)
};
World.Tile.prototype.setTileInfoByGrid = function (level, row, column) {
	if (level != undefined && row != undefined && column != undefined) {
		var Egeo = World.Math.getTileGeographicEnvelopByGrid(level, row, column);
		this.minLon = Egeo.minLon;
		this.minLat = Egeo.minLat;
		this.maxLon = Egeo.maxLon;
		this.maxLat = Egeo.maxLat;
		if (this.level < World.BASE_LEVEL) {
			var changeLevel = World.BASE_LEVEL - this.level;
			this.segment = Math.pow(2, changeLevel)
		} else {
			this.segment = 1
		}
		this.material = new World.TextureMaterial(this.url)
	}
};
World.Tile.prototype.destroy = function () {
	World.Object3D.prototype.destroy.apply(this, arguments);
	this.subTiledLayer = null
};
World.SubTiledLayer = function (args) {
	World.Object3DComponents.apply(this, arguments);
	this.level = -1;
	this.tiledLayer = null;
	if (args) {
		if (args.level != undefined) {
			this.level = args.level
		}
	}
};
World.SubTiledLayer.prototype = new World.Object3DComponents();
World.SubTiledLayer.prototype.constructor = World.SubTiledLayer;
World.SubTiledLayer.prototype.add = function (tile) {
	if (tile instanceof World.Tile) {
		if (tile.level == this.level) {
			World.Object3DComponents.prototype.add.apply(this, arguments);
			tile.subTiledLayer = this
		}
	}
};
World.SubTiledLayer.prototype.getImageUrl = function (level, row, column) {
	var url = "";
	if (this.tiledLayer) {
		url = this.tiledLayer.getImageUrl(level, row, column)
	}
	return url
};
World.SubTiledLayer.prototype.destroy = function () {
	World.Object3DComponents.prototype.destroy.apply(this, arguments);
	this.tiledLayer = null
};
World.SubTiledLayer.prototype.checkIfLoaded = function () {
	for (var i = 0; i < this.objectList.length; i++) {
		var tile = this.objectList[i];
		if (tile) {
			var isTileLoaded = tile.material.texture.loaded;
			if (!isTileLoaded) {
				return false
			}
		}
	}
	return true
};
World.TiledLayer = function (args) {
	World.Object3DComponents.apply(this, arguments)
};
World.TiledLayer.prototype = new World.Object3DComponents();
World.TiledLayer.prototype.constructor = World.TiledLayer;
World.TiledLayer.prototype.add = function (subTiledLayer) {
	if (subTiledLayer instanceof World.SubTiledLayer) {
		World.Object3DComponents.prototype.add.apply(this, arguments);
		subTiledLayer.tiledLayer = this
	}
};
World.TiledLayer.prototype.draw = function (camera) {
	World.Object3DComponents.prototype.draw.apply(this, arguments)
};
World.TiledLayer.prototype.getImageUrl = function (level, row, column) {
	return ""
};
World.GoogleTiledLayer = function (args) {
	World.TiledLayer.apply(this, arguments)
};
World.GoogleTiledLayer.prototype = new World.TiledLayer();
World.GoogleTiledLayer.prototype.constructor = World.GoogleTiledLayer;
World.GoogleTiledLayer.prototype.getImageUrl = function (level, row, column) {
	var sum = level + row + column;
	var idx = 1 + sum % 3;
	var url = "http://mt" + idx + ".google.cn/vt/lyrs=m@212000000&hl=zh-CN&gl=CN&src=app&x=" + column + "&y=" + row + "&z=" + level + "&s=Galil";
	return url
};
World.OsmTiledLayer = function (args) {
	World.TiledLayer.apply(this, arguments)
};
World.OsmTiledLayer.prototype = new World.TiledLayer();
World.OsmTiledLayer.prototype.constructor = World.OsmTiledLayer;
World.OsmTiledLayer.prototype.getImageUrl = function (level, row, column) {
	var sum = level + row + column;
	var idx = 1 + sum % 4;
	var url = "http://otile" + idx + ".mqcdn.com/tiles/1.0.0/osm/" + level + "/" + column + "/" + row + ".jpg";
	return url
};
World.NokiaTiledLayer = function (args) {
	World.TiledLayer.apply(this, arguments)
};
World.NokiaTiledLayer.prototype = new World.TiledLayer();
World.NokiaTiledLayer.prototype.constructor = World.NokiaTiledLayer;
World.NokiaTiledLayer.prototype.getImageUrl = function (level, row, column) {
	var sum = level + row + column;
	var idx = 1 + sum % 4;
	var url = "http://" + idx + ".maps.nlp.nokia.com.cn/maptile/2.1/maptile/newest/satellite.day/" + level + "/" + column + "/" + row + "/256/jpg?token=BIl5zlMQF2fUaQLYWcxE&app_id=tGvvOZHNwj1-4guzYwIU";
	return url
};
World.BingTiledLayer = function (args) {
	World.TiledLayer.apply(this, arguments)
};
World.BingTiledLayer.prototype = new World.TiledLayer();
World.BingTiledLayer.prototype.constructor = World.BingTiledLayer;
World.BingTiledLayer.prototype.getImageUrl = function (level, row, column) {
	var url = "";
	var tileX = column;
	var tileY = row;
	var strTileX2 = World.Math.numerationSystemFrom10(2, tileX);
	var strTileY2 = World.Math.numerationSystemFrom10(2, tileY);
	var delta = strTileX2.length - strTileY2.length;
	if (delta > 0) {
		for (var i = 0; i < delta; i++) {
			strTileY2 = '0' + strTileY2
		}
	} else if (delta < 0) {
		delta = -delta;
		for (var i = 0; i < delta; i++) {
			strTileX2 = '0' + strTileX2
		}
	}
	var strMerge2 = "";
	for (var i = 0; i < strTileY2.length; i++) {
		var charY = strTileY2[i];
		var charX = strTileX2[i];
		strMerge2 += charY + charX
	}
	var strMerge4 = World.Math.numerationSystemChange(2, 4, strMerge2);
	if (strMerge4.length < level) {
		var delta = level - strMerge4.length;
		for (var i = 0; i < delta; i++) {
			strMerge4 = '0' + strMerge4
		}
	}
	var sum = level + row + column;
	var serverIdx = sum % 8;
	url = "http://ecn.t" + serverIdx + ".tiles.virtualearth.net/tiles/h" + strMerge4 + ".jpeg?g=1239&mkt=en-us";
	return url
};
World.TiandituTiledLayer = function (args) {
	World.TiledLayer.apply(this, arguments)
};
World.TiandituTiledLayer.prototype = new World.TiledLayer();
World.TiandituTiledLayer.prototype.constructor = World.TiandituTiledLayer;
World.TiandituTiledLayer.prototype.getImageUrl = function (level, row, column) {
	var url = "";
	var sum = level + row + column;
	var serverIdx = sum % 8;
	url = "http://t" + serverIdx + ".tianditu.com/DataServer?T=vec_w&x=" + column + "&y=" + row + "&l=" + level;
	return url
};
World.SosoTiledLayer = function (args) {
	World.TiledLayer.apply(this, arguments)
};
World.SosoTiledLayer.prototype = new World.TiledLayer();
World.SosoTiledLayer.prototype.constructor = World.SosoTiledLayer;
World.SosoTiledLayer.prototype.getImageUrl = function (level, row, column) {
	var url = "";
	var tileCount = Math.pow(2, level);
	var a = column;
	var b = tileCount - row - 1;
	var A = Math.floor(a / 16);
	var B = Math.floor(b / 16);
	var sum = level + row + column;
	var serverIdx = sum % 4;
	var sateUrl = "http://p" + serverIdx + ".map.soso.com/sateTiles/" + level + "/" + A + "/" + B + "/" + a + "_" + b + ".jpg";
	url = sateUrl;
	return url
};
World.BlendTiledLayer = function (args) {
	World.TiledLayer.apply(this, arguments)
};
World.BlendTiledLayer.prototype = new World.TiledLayer();
World.BlendTiledLayer.prototype.constructor = World.BlendTiledLayer;
World.BlendTiledLayer.prototype.getImageUrl = function (level, row, column) {
	var array = [World.NokiaTiledLayer, World.GoogleTiledLayer, World.OsmTiledLayer];
	var sum = level + row + column;
	var idx = sum % 3;
	var url = array[idx].prototype.getImageUrl.apply(this, arguments);
	return url
};
World.ArcGISTiledLayer = function (args) {
	World.TiledLayer.apply(this, arguments);
	this.service = "";
	if (args) {
		if (args.url) {
			this.service = args.url
		}
	}
};
World.ArcGISTiledLayer.prototype = new World.TiledLayer();
World.ArcGISTiledLayer.prototype.constructor = World.ArcGISTiledLayer;
World.ArcGISTiledLayer.prototype.getImageUrl = function (level, row, column) {
	var url = "proxy.jsp?" + this.service + "/tile/" + level + "/" + row + "/" + column;
	return url
};
World.AutonaviTiledLayer = function (args) {
	World.TiledLayer.apply(this, arguments)
};
World.AutonaviTiledLayer.prototype = new World.TiledLayer();
World.AutonaviTiledLayer.prototype.constructor = World.AutonaviTiledLayer;
World.AutonaviTiledLayer.prototype.getImageUrl = function (level, row, column) {
	var sum = level + row + column;
	var serverIdx = 1 + sum % 4;
	var url = "proxy.jsp?http://webrd0" + serverIdx + ".is.autonavi.com/appmaptile?x=" + column + "&y=" + row + "&z=" + level + "&lang=zh_cn&size=1&scale=1&style=8";
	return url
};
World.TestTiledLayer = function (args) {
	World.TiledLayer.apply(this, arguments)
};
World.TestTiledLayer.prototype = new World.TiledLayer();
World.TestTiledLayer.prototype.constructor = World.TestTiledLayer;
World.TestTiledLayer.prototype.getImageUrl = function (level, row, column) {
	var url = "";
	return url
};
World.PerspectiveCamera = function (fov, aspect, near, far) {
	World.Object3D.apply(this, null);
	this.fov = fov || 90;
	this.aspect = aspect || 1;
	this.near = near || 0.1;
	this.far = far || 100;
	this.projMatrix = new World.Matrix();
	this.setPerspectiveMatrix(this.fov, this.aspect, this.near, this.far)
};
World.PerspectiveCamera.prototype = new World.Object3D();
World.PerspectiveCamera.prototype.constructor = World.PerspectiveCamera;
World.PerspectiveCamera.prototype.setPerspectiveMatrix = function (fov, aspect, near, far) {
	this.fov = fov || this.fov;
	this.aspect = aspect || this.aspect;
	this.near = near || this.near;
	this.far = far || this.far;
	var mat = [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1];
	var f = this.fov * Math.PI / 360;
	var a = this.far - this.near;
	var e = Math.cos(f) / Math.sin(f);
	mat[0] = e / this.aspect;
	mat[5] = e;
	mat[10] =  - (this.far + this.near) / a;
	mat[11] = -1;
	mat[14] = -2 * this.near * this.far / a;
	mat[15] = 0;
	this.projMatrix.setElements(mat[0], mat[1], mat[2], mat[3], mat[4], mat[5], mat[6], mat[7], mat[8], mat[9], mat[10], mat[11], mat[12], mat[13], mat[14], mat[15])
};
World.PerspectiveCamera.prototype.getLightDirection = function () {
	var dirVertice = this.matrix.getColumnZ();
	var direction = new World.Vector(-dirVertice.x, -dirVertice.y, -dirVertice.z);
	direction.normalize();
	return direction
};
World.PerspectiveCamera.prototype.getProjViewMatrix = function () {
	var viewMatrix = this.getViewMatrix();
	var projViewMatrix = this.projMatrix.multiplyMatrix(viewMatrix);
	return projViewMatrix
};
World.PerspectiveCamera.prototype.convertVerticeFromWorldToNDC = function (verticeInWorld, projViewMatrix) {
	if (!(projViewMatrix instanceof World.Matrix)) {
		projViewMatrix = this.getProjViewMatrix()
	}
	var columnWorld = [verticeInWorld.x, verticeInWorld.y, verticeInWorld.z, 1];
	var columnProject = projViewMatrix.multiplyColumn(columnWorld);
	var w = columnProject[3];
	var columnNDC = [];
	columnNDC[0] = columnProject[0] / w;
	columnNDC[1] = columnProject[1] / w;
	columnNDC[2] = columnProject[2] / w;
	columnNDC[3] = 1;
	var verticeInNDC = new World.Vertice(columnNDC[0], columnNDC[1], columnNDC[2]);
	return verticeInNDC
};
World.PerspectiveCamera.prototype.converPoint3DToPoint2D = function (verticeInWorld) {
	var W = World.Canvas.width;
	var H = World.Canvas.height;
	var verticeInNDC = this.convertVerticeFromWorldToNDC(verticeInWorld);
	var x = (1 + verticeInNDC.x) * W / 2;
	var y = (1 - verticeInNDC.y) * H / 2;
	var z = verticeInNDC.z;
	var result = new World.Vertice(x, y, z);
	return result
};
World.PerspectiveCamera.prototype.convertVerticeFromCameraToWorld = function (verticeInCamera) {
	var verticeInCameraCopy = verticeInCamera.getCopy();
	var viewMatrix = this.getViewMatrix();
	var inverseMatrix = viewMatrix.getInverseMatrix();
	var column = [verticeInCameraCopy.x, verticeInCameraCopy.y, verticeInCameraCopy.z, 1];
	var column2 = inverseMatrix.multiplyColumn(column);
	var verticeInWorld = new World.Vertice(column2[0], column2[1], column2[2]);
	return verticeInWorld
};
World.PerspectiveCamera.prototype.convertVectorFromCameraToWorld = function (vectorInCamera) {
	var vectorInCameraCopy = vectorInCamera.getCopy();
	var verticeInCamera = vectorInCameraCopy.getVertice();
	var verticeInWorld = this.convertVerticeFromCameraToWorld(verticeInCamera);
	var originInWorld = this.getPosition();
	var vectorInWorld = verticeInWorld.minus(originInWorld);
	vectorInWorld.normalize();
	return vectorInWorld
};
World.PerspectiveCamera.prototype.setFov = function () {
	this.setPerspectiveMatrix(fov, this.aspect, this.near, this.far)
};
World.PerspectiveCamera.prototype.setAspect = function (aspect) {
	this.setPerspectiveMatrix(this.fov, aspect, this.near, this.far)
};
World.PerspectiveCamera.prototype.setNear = function (near) {
	this.setPerspectiveMatrix(this.fov, this.aspect, near, this.far)
};
World.PerspectiveCamera.prototype.setFar = function (far) {
	this.setPerspectiveMatrix(this.fov, this.aspect, this.near, far)
};
World.PerspectiveCamera.prototype.getViewMatrix = function () {
	var columnTrans = this.matrix.getColumnTrans();
	var transX = columnTrans.x;
	var transY = columnTrans.y;
	var transZ = columnTrans.z;
	var mat1 = new World.Matrix();
	mat1.setMatrixByOther(this.matrix);
	mat1.transpose();
	mat1.setColumnTrans(0, 0, 0);
	mat1.setLastRowDefault();
	var mat2 = new World.Matrix();
	mat2.setColumnTrans(-transX, -transY, -transZ);
	var viewMatrix = mat1.multiplyMatrix(mat2);
	return viewMatrix
};
World.PerspectiveCamera.prototype.look = function (cameraPnt, targetPnt, upDirection) {
	var cameraPntCopy = cameraPnt.getCopy();
	var targetPntCopy = targetPnt.getCopy();
	upDirection = upDirection || new World.Vector(0, 1, 0);
	var up = upDirection.getCopy();
	var transX = cameraPntCopy.x;
	var transY = cameraPntCopy.y;
	var transZ = cameraPntCopy.z;
	var zAxis = new World.Vector(cameraPntCopy.x - targetPntCopy.x, cameraPntCopy.y - targetPntCopy.y, cameraPntCopy.z - targetPntCopy.z).normalize();
	var xAxis = up.cross(zAxis).normalize();
	var yAxis = zAxis.cross(xAxis).normalize();
	this.matrix.setColumnX(xAxis.x, xAxis.y, xAxis.z);
	this.matrix.setColumnY(yAxis.x, yAxis.y, yAxis.z);
	this.matrix.setColumnZ(zAxis.x, zAxis.y, zAxis.z);
	this.matrix.setColumnTrans(transX, transY, transZ);
	this.matrix.setLastRowDefault();
	var deltaX = cameraPntCopy.x - targetPntCopy.x;
	var deltaY = cameraPntCopy.y - targetPntCopy.y;
	var deltaZ = cameraPntCopy.z - targetPntCopy.z;
	var far = Math.sqrt(deltaX * deltaX + deltaY * deltaY + deltaZ * deltaZ);
	this.setFar(far)
};
World.PerspectiveCamera.prototype.lookAt = function (targetPnt, upDirection) {
	var targetPntCopy = targetPnt.getCopy();
	var upDirectionCopy = upDirection.getCopy();
	var position = this.getPosition();
	this.look(position, targetPntCopy, upDirectionCopy)
};
World.PerspectiveCamera.prototype.getPickDirection = function (absoluteX, absoluteY) {
	var relativeCoords = World.Math.getPositionRelativeToCanvasCenter(absoluteX, absoluteY);
	var relativeX = relativeCoords[0];
	var relativeY = relativeCoords[1];
	var ndcX = relativeX / (World.Canvas.width / 2);
	var ndcY = relativeY / (World.Canvas.height / 2);
	var ndcZ = 0.5;
	var ndcW = 1;
	var columnNDC = [ndcX, ndcY, ndcZ, ndcW];
	var inverseProj = this.projMatrix.getInverseMatrix();
	var columnCameraTemp = inverseProj.multiplyColumn(columnNDC);
	var cameraX = columnCameraTemp[0] / columnCameraTemp[3];
	var cameraY = columnCameraTemp[1] / columnCameraTemp[3];
	var cameraZ = columnCameraTemp[2] / columnCameraTemp[3];
	var cameraW = 1;
	var columnCamera = [cameraX, cameraY, cameraZ, cameraW];
	var viewMatrix = this.getViewMatrix();
	var inverseView = viewMatrix.getInverseMatrix();
	var columnWorld = inverseView.multiplyColumn(columnCamera);
	var verticeInWorld = new World.Vertice(columnWorld[0], columnWorld[1], columnWorld[2]);
	var cameraPositon = this.getPosition();
	var pickDirection = verticeInWorld.minus(cameraPositon);
	pickDirection.normalize();
	return pickDirection
};
World.PerspectiveCamera.prototype.getPickCartesianCoordInEarth = function (absoluteX, absoluteY) {
	var result = [];
	var pickDirection = this.getPickDirection(absoluteX, absoluteY);
	var cameraVertice = this.getPosition();
	var line = new World.Line(cameraVertice, pickDirection);
	var pickVertices = World.Math.getLineIntersectPointWidthEarth(line);
	if (pickVertices.length == 0) {
		result = []
	} else if (pickVertices.length == 1) {
		var pickVertice = pickVertices[0];
		result = [pickVertice]
	} else if (pickVertices.length == 2) {
		var pickVerticeA = pickVertices[0];
		var pickVerticeB = pickVertices[1];
		var lengthA = World.Math.getLengthFromVerticeToVertice(cameraVertice, pickVerticeA);
		var lengthB = World.Math.getLengthFromVerticeToVertice(cameraVertice, pickVerticeB);
		if (lengthA <= lengthB) {
			result = [pickVerticeA, pickVerticeB]
		} else {
			result = [pickVerticeB, pickVerticeA]
		}
	}
	return result
};
World.PerspectiveCamera.prototype.getLengthFromCameraToOrigin = function () {
	var origin = new World.Vertice(0, 0, 0);
	var cameraPosition = this.getPosition();
	var length = World.Math.getLengthFromVerticeToVertice(origin, cameraPosition);
	return length
};
World.PerspectiveCamera.prototype.getPlanX0Z = function () {
	var position = this.getPosition();
	var direction = this.getLightDirection();
	var plan = World.Math.getCrossPlaneByLine(position, direction);
	return plan
};
World.PerspectiveCamera.prototype.getEarthInOutInfo = function () {
	var result = {
		status : "",
		ltPickResult : [],
		rtPickResult : [],
		rbPickResult : [],
		lbPickResult : [],
		lcPickResult : [],
		tcPickResult : [],
		rcPickResult : [],
		bcPickResult : []
	};
	var WIDTH = World.Canvas.width;
	var HEIGHT = World.Canvas.height;
	result.ltPickResult = this.getPickCartesianCoordInEarth(0, 0);
	result.rtPickResult = this.getPickCartesianCoordInEarth(WIDTH, 0);
	result.rbPickResult = this.getPickCartesianCoordInEarth(WIDTH, HEIGHT);
	result.lbPickResult = this.getPickCartesianCoordInEarth(0, HEIGHT);
	result.lcPickResult = this.getPickCartesianCoordInEarth(0, HEIGHT / 2.0);
	result.tcPickResult = this.getPickCartesianCoordInEarth(WIDTH / 2.0, 0);
	result.rcPickResult = this.getPickCartesianCoordInEarth(WIDTH, HEIGHT / 2.0);
	result.bcPickResult = this.getPickCartesianCoordInEarth(WIDTH / 2.0, HEIGHT);
	var isFullOut = result.ltPickResult.length > 0 && result.rtPickResult.length > 0 && result.rbPickResult.length > 0 && result.lbPickResult.length > 0;
	var isFullIn = result.lcPickResult.length == 0 && result.tcPickResult.length == 0 && result.rcPickResult.length == 0 && result.bcPickResult.length == 0;
	if (isFullOut && !isFullIn) {
		result.status = World.Enum.FULL_OUT
	} else if (!isFullOut && isFullIn) {
		result.status = World.Enum.FULL_IN
	} else {
		result.status = World.Enum.IN_OUT
	}
	return result
};
World.PerspectiveCamera.prototype.setLevel = function (level) {
	var length = World.Math.getLengthFromCamera2EarthSurface(level) + World.EARTH_RADIUS;
	var origin = new World.Vertice(0, 0, 0);
	var oldPosition = this.getPosition();
	var vector = oldPosition.minus(origin);
	var newPosition = null;
	if (vector.getLength() != 0) {
		vector.normalize();
		vector.setLength(length);
		newPosition = vector.getVertice()
	} else {
		newPosition = new World.Vertice(0, 0, length)
	}
	this.look(newPosition, origin);
	World.CURRENT_LEVEL = level
};
World.PerspectiveCamera.prototype.getEarthTangencyCircleRadius = function () {
	var R = World.EARTH_RADIUS;
	var L = this.getLengthFromCameraToOrigin();
	var length = R / L * Math.sqrt(L * L - R * R);
	return length
};
World.PerspectiveCamera.prototype.getLengthFromCameraToTangencyCircle = function () {
	var R = World.EARTH_RADIUS;
	var L = this.getLengthFromCameraToOrigin();
	var distance = L - R * R / L;
	return distance
};
World.PerspectiveCamera.prototype.getIntersectPointByLightDirectionAndTangencyCircle = function () {
	var lengthToTangencyCircle = this.getLengthFromCameraToTangencyCircle();
	var lightDirection = this.getLightDirection();
	lightDirection.setLength(lengthToTangencyCircle);
	var cameraPosition = this.getPosition();
	var intersectPoint = cameraPosition.plus(lightDirection);
	return intersectPoint
};
World.PerspectiveCamera.prototype.isTileVisible = function (level, row, column, geoExtent, projViewMatrix) {
	if (level == 0 && row == 0 && column == 0) {
		return true
	}
	if (!geoExtent) {
		geoExtent = this.getCurrentGeoExtent()
	}
	if (!(projViewMatrix instanceof World.Matrix)) {
		if (geoExtent.projViewMatrix instanceof World.Matrix) {
			projViewMatrix = geoExtent.projViewMatrix
		}
	}
	var earthLonLeft = geoExtent.lonLeft;
	var earthLonRight = geoExtent.lonRight;
	var earthLatBottom = geoExtent.latBottom;
	var earthLatTop = geoExtent.latTop;
	var Egeo = World.Math.getTileGeographicEnvelopByGrid(level, row, column);
	var tileMinLon = Egeo.minLon;
	var tileMaxLon = Egeo.maxLon;
	var tileMinLat = Egeo.minLat;
	var tileMaxLat = Egeo.maxLat;
	var isLatVisible = this.isLatVisible(earthLatBottom, earthLatTop, tileMinLat, tileMaxLat);
	if (!isLatVisible) {
		return false
	}
	if (earthLatBottom != -90 && earthLatTop != 90) {
		var isLonVisible = this.isLonVisible(earthLonLeft, earthLonRight, tileMinLon, tileMaxLon);
		if (!isLonVisible) {
			return false
		}
	}
	if (geoExtent.status != World.Enum.FULL_IN) {
		var pLeftBottom = World.Math.geographicToCartesianCoord(tileMinLon, tileMinLat);
		var pLeftTop = World.Math.geographicToCartesianCoord(tileMinLon, tileMaxLat);
		var pRightTop = World.Math.geographicToCartesianCoord(tileMaxLon, tileMaxLat);
		var pRightBottom = World.Math.geographicToCartesianCoord(tileMaxLon, tileMinLat);
		var vertices = [pLeftBottom, pLeftTop, pRightTop, pRightBottom];
		var ndcVertices = [];
		var outCanvasCountArray = [0, 0, 0, 0];
		if (!(projViewMatrix instanceof World.Matrix)) {
			projViewMatrix = this.getProjViewMatrix()
		}
		for (var i = 0; i < vertices.length; i++) {
			var verticeInWorld = vertices[i];
			var verticeNDC = this.convertVerticeFromWorldToNDC(verticeInWorld, projViewMatrix);
			ndcVertices.push(verticeNDC);
			var ndcX = verticeNDC.x;
			var ndcY = verticeNDC.y;
			var ndcZ = verticeNDC.z;
			if (ndcX < -1) {
				outCanvasCountArray[0]++
			} else if (ndcY > 1) {
				outCanvasCountArray[1]++
			} else if (ndcX > 1) {
				outCanvasCountArray[2]++
			} else if (ndcY < -1) {
				outCanvasCountArray[3]++
			}
		}
		var countLeft = outCanvasCountArray[0];
		var countTop = outCanvasCountArray[1];
		var countRight = outCanvasCountArray[2];
		var countBottom = outCanvasCountArray[3];
		var isFullOut = countLeft == 4 || countTop == 4 || countRight == 4 || countBottom == 4;
		if (isFullOut) {
			return false
		}
	}
	return true
};
World.PerspectiveCamera.prototype.isLonVisible = function (earthLonLeft, earthLonRight, tileMinLon, tileMaxLon) {
	var isvisible = false;
	var isTileLonPositive = tileMaxLon > 0;
	if (earthLonLeft >= 0 && earthLonRight >= 0) {
		if (isTileLonPositive) {
			isvisible = this._isNumberIntersect(earthLonLeft, earthLonRight, tileMinLon, tileMaxLon)
		} else {
			isvisible = false
		}
	} else if (earthLonLeft <= 0 && earthLonRight <= 0) {
		if (isTileLonPositive) {
			isvisible = false
		} else {
			isvisible = this._isNumberIntersect(earthLonLeft, earthLonRight, tileMinLon, tileMaxLon)
		}
	} else if (earthLonLeft < 0 && earthLonRight > 0) {
		if (isTileLonPositive) {
			earthLonLeft = 0;
			isvisible = this._isNumberIntersect(earthLonLeft, earthLonRight, tileMinLon, tileMaxLon)
		} else {
			earthLonRight = 0;
			isvisible = this._isNumberIntersect(earthLonLeft, earthLonRight, tileMinLon, tileMaxLon)
		}
	} else if (earthLonLeft > 0 && earthLonRight < 0) {
		if (isTileLonPositive) {
			earthLonRight = 180;
			isvisible = this._isNumberIntersect(earthLonLeft, earthLonRight, tileMinLon, tileMaxLon)
		} else {
			earthLonLeft = -180;
			isvisible = this._isNumberIntersect(earthLonLeft, earthLonRight, tileMinLon, tileMaxLon)
		}
	}
	return isvisible
};
World.PerspectiveCamera.prototype.isLatVisible = function (earthLatBottom, earthLatTop, tileMinLat, tileMaxLat) {
	var isvisible = this._isNumberIntersect(earthLatBottom, earthLatTop, tileMinLat, tileMaxLat);
	return isvisible
};
World.PerspectiveCamera.prototype._isNumberIntersect = function (from1, to1, from2, to2) {
	var min1 = Math.min(from1, to1);
	var max1 = Math.max(from1, to1);
	var min2 = Math.min(from2, to2);
	var max2 = Math.max(from2, to2);
	var isNotIntersect = max1 <= min2 || max2 <= min1;
	return !isNotIntersect
};
World.PerspectiveCamera.prototype.getCurrentGeoExtent = function (earthInOutInfo) {
	var result = {
		status : "",
		lonLeft : null,
		lonRight : null,
		latBottom : null,
		latTop : null,
		lonCenter : null,
		latCenter : null,
		earthInOutInfo : null,
		projViewMatrix : null
	};
	if (!earthInOutInfo) {
		earthInOutInfo = this.getEarthInOutInfo()
	}
	result.earthInOutInfo = earthInOutInfo;
	result.status = earthInOutInfo.status;
	if (earthInOutInfo.status == World.Enum.FULL_OUT) {
		var ltGeo = World.Math.cartesianCoordToGeographic(earthInOutInfo.ltPickResult[0]);
		var rtGeo = World.Math.cartesianCoordToGeographic(earthInOutInfo.rtPickResult[0]);
		var rbGeo = World.Math.cartesianCoordToGeographic(earthInOutInfo.rbPickResult[0]);
		var lbGeo = World.Math.cartesianCoordToGeographic(earthInOutInfo.lbPickResult[0]);
		var lcGeo = World.Math.cartesianCoordToGeographic(earthInOutInfo.lcPickResult[0]);
		var tcGeo = World.Math.cartesianCoordToGeographic(earthInOutInfo.tcPickResult[0]);
		var rcGeo = World.Math.cartesianCoordToGeographic(earthInOutInfo.rcPickResult[0]);
		var bcGeo = World.Math.cartesianCoordToGeographic(earthInOutInfo.bcPickResult[0]);
		var isPoleIn = false;
		var isNorth;
		var minLatValue = Math.min(ltGeo[1], rtGeo[1], rbGeo[1], lbGeo[1], lcGeo[1], tcGeo[1], rcGeo[1], bcGeo[1]);
		var maxLatValue = Math.max(ltGeo[1], rtGeo[1], rbGeo[1], lbGeo[1], lcGeo[1], tcGeo[1], rcGeo[1], bcGeo[1]);
		if (maxLatValue < 0 || minLatValue > 0) {
			var geoPole = null;
			if (minLatValue > 0) {
				geoPole = [0, 90];
				isNorth = true
			} else if (maxLatValue < 0) {
				geoPole = [0, -90];
				isNorth = false
			}
			if (geoPole) {
				var poleInWorld = World.Math.geographicToCartesianCoord(geoPole[0], geoPole[1]);
				result.projViewMatrix = this.getProjViewMatrix();
				var ndcPole = this.convertVerticeFromWorldToNDC(poleInWorld, result.projViewMatrix);
				if (ndcPole.x >= -1 && ndcPole.x <= 1 && ndcPole.y >= -1 && ndcPole.y <= 1) {
					isPoleIn = true
				} else {
					isPoleIn = false
				}
			}
		}
		if (isPoleIn) {
			if (isNorth) {
				result.latTop = 90;
				result.latBottom = minLatValue
			} else {
				result.latBottom = -90;
				result.latTop = maxLatValue
			}
		} else {
			result.latBottom = minLatValue;
			result.latTop = maxLatValue
		}
		var minLeftValue = Math.min(ltGeo[0], lbGeo[0], lcGeo[0]);
		var maxLeftValue = Math.max(ltGeo[0], lbGeo[0], lcGeo[0]);
		if (maxLeftValue <= 0 || minLeftValue >= 0) {
			result.lonLeft = minLeftValue
		} else {
			var averLeftValue = (Math.abs(ltGeo[0]) + Math.abs(lbGeo[0]) + Math.abs(lcGeo[0])) / 3.0;
			var a = Math.abs(averLeftValue);
			var b = Math.abs(180 - averLeftValue);
			if (a <= b) {
				result.lonLeft = minLeftValue
			} else {
				result.lonLeft = maxLeftValue
			}
		}
		var minRightValue = Math.min(rtGeo[0], rbGeo[0], rcGeo[0]);
		var maxRightValue = Math.max(rtGeo[0], rbGeo[0], rcGeo[0]);
		if (maxRightValue <= 0 || minRightValue >= 0) {
			result.lonRight = maxRightValue
		} else {
			var averRightValue = (Math.abs(rtGeo[0]) + Math.abs(rbGeo[0]) + Math.abs(rcGeo[0])) / 3.0;
			var a = Math.abs(averRightValue);
			var b = Math.abs(180 - averRightValue);
			if (a <= b) {
				result.lonRight = maxRightValue
			} else {
				result.lonRight = minRightValue
			}
		}
	} else {
		var centerX = World.Canvas.width / 2.0;
		var centerY = World.Canvas.height / 2.0;
		var centerPickResult = this.getPickCartesianCoordInEarth(centerX, centerY);
		var centerVertice = centerPickResult[0];
		var ctGeo = World.Math.cartesianCoordToGeographic(centerVertice);
		var ctLon = ctGeo[0];
		var ctLat = ctGeo[1];
		var L = this.getLengthFromCameraToOrigin();
		var R = World.EARTH_RADIUS;
		var degree = World.Math.radianToDegree(Math.acos(R / L));
		result.lonLeft = ctLon - degree;
		result.lonRight = ctLon + degree;
		result.latBottom = ctLat - degree;
		result.latTop = ctLat + degree
	}
	if (result.lonLeft > 0) {
		if (result.lonLeft > 180) {
			result.lonLeft = result.lonLeft - 360
		}
	} else if (result.lonLeft < 0) {
		if (result.lonLeft < -180) {
			result.lonLeft = result.lonLeft + 360
		}
	}
	if (result.lonRight > 0) {
		if (result.lonRight > 180) {
			result.lonRight = result.lonRight - 360
		}
	} else if (result.lonRight < 0) {
		if (result.lonRight < -180) {
			result.lonRight = result.lonRight + 360
		}
	}
	if (result.latBottom > 0) {
		result.latBottom = Math.min(result.latBottom, 90)
	} else if (result.latBottom < 0) {
		result.latBottom = Math.max(result.latBottom, -90)
	}
	if (result.latTop > 0) {
		result.latTop = Math.min(result.latTop, 90)
	} else if (result.latTop < 0) {
		result.latTop = Math.max(result.latTop, -90)
	}
	result.latBottom = Math.min(result.latBottom, result.latTop);
	result.latTop = Math.max(result.latBottom, result.latTop);
	return result
};
World.PerspectiveCamera.prototype.getVisibleTileIndex = function (level, geoExtent, projViewMatrix) {
	if (!geoExtent) {
		geoExtent = this.getCurrentGeoExtent()
	}
	if (!(projViewMatrix instanceof World.Matrix)) {
		projViewMatrix = this.getProjViewMatrix()
	}
	var sumTilesLevel = [[]];
	if (level == 0) {
		var tile = {
			index : [0, 0, 0],
			parentIndex : null
		};
		sumTilesLevel = [[tile]]
	} else if (level > 0) {
		var sumTilesLevel_1 = this.getVisibleTileIndex(level - 1, geoExtent, projViewMatrix);
		sumTilesLevel = sumTilesLevel_1;
		sumTilesLevel.push([]);
		var singleTilesLevel_1 = sumTilesLevel[level - 1];
		for (var i = 0; i < singleTilesLevel_1.length; i++) {
			var tile = singleTilesLevel_1[i];
			var childTileLT = World.Math.getTileGridByParent(tile.index[0], tile.index[1], tile.index[2], World.Math.LEFT_TOP);
			var childTileRT = World.Math.getTileGridByParent(tile.index[0], tile.index[1], tile.index[2], World.Math.RIGHT_TOP);
			var childTileLB = World.Math.getTileGridByParent(tile.index[0], tile.index[1], tile.index[2], World.Math.LEFT_BOTTOM);
			var childTileRB = World.Math.getTileGridByParent(tile.index[0], tile.index[1], tile.index[2], World.Math.RIGHT_BOTTOM);
			var childTiles = [childTileLT, childTileRT, childTileLB, childTileRB];
			for (var j = 0; j < childTiles.length; j++) {
				var childTile = childTiles[j];
				var isvisible = this.isTileVisible(childTile.level, childTile.row, childTile.column, geoExtent, projViewMatrix);
				if (isvisible) {
					var childTile = {
						index : [childTile.level, childTile.row, childTile.column],
						parentIndex : [tile.index[0], tile.index[1], tile.index[2]]
					};
					sumTilesLevel[level].push(childTile)
				}
			}
		}
	}
	sumTilesLevel.geoExtent = geoExtent;
	sumTilesLevel.projViewMatrix = projViewMatrix;
	return sumTilesLevel
};
World.Globe = function (canvas, args) {
	var vertexShaderContent = World.ShaderContent.SIMPLE_SHADER.VS_CONTENT;
	var fragmentShaderContent = World.ShaderContent.SIMPLE_SHADER.FS_CONTENT;
	this.renderer = null;
	this.scene = null;
	this.camera = null;
	this.tiledLayer = null;
	this.renderer = new World.WebGLRenderer(canvas, vertexShaderContent, fragmentShaderContent);
	this.scene = new World.Scene();
	var radio = canvas.width / canvas.height;
	this.camera = new World.PerspectiveCamera(30, radio, 1.0, 20000000.0);
	this.cameraMoveInfo = {
		horDegree : null,
		verDegree : null
	};
	this.setLevel(0);
	this.renderer.bindScene(this.scene);
	this.renderer.bindCamera(this.camera);
	this.renderer.setIfAutoRefresh(true)
};
World.Globe.prototype = {
	constructor : World.Globe,
	setTiledLayer : function (tiledLayer) {
		if (tiledLayer instanceof World.TiledLayer) {
			if (this.tiledLayer) {
				var b = this.scene.remove(this.tiledLayer);
				if (!b) {
					console.debug("LINE:2887,this.scene.remove(this.tiledLayer)失败")
				}
				this.scene.tiledLayer = null
			}
			this.tiledLayer = tiledLayer;
			this.scene.add(this.tiledLayer);
			this.refresh(true)
		}
	},
	setLevel : function (level) {
		level = (level >= 0) ? level : 0;
		level = level > 18 ? 18 : level;
		if (level != World.CURRENT_LEVEL) {
			if (this.camera) {
				this.camera.setLevel(level);
				var geoExtent = this.camera.getCurrentGeoExtent();
				this.calculateCameraMoveInfo(geoExtent);
				this.showBelowLevelSubTiledLayer(true, geoExtent);
				this.refresh(false, geoExtent)
			}
		}
	},
	calculateCameraMoveInfo : function (geoExtent) {
		geoExtent = geoExtent ? geoExtent : this.camera.getCurrentGeoExtent();
		this.cameraMoveInfo = {
			horDegree : null,
			verDegree : null
		};
		if (geoExtent.status != World.Enum.FULL_OUT) {
			this.cameraMoveInfo.horDegree = World.Canvas.width / World.Canvas.height * this.camera.fov;
			this.cameraMoveInfo.verDegree = this.camera.fov
		} else {
			var D = 2 * World.EARTH_RADIUS;
			var earthInOutInfo = geoExtent.earthInOutInfo;
			var lcVertice = earthInOutInfo.lcPickResult[0];
			var rcVertice = earthInOutInfo.rcPickResult[0];
			var horLength = World.Math.getLengthFromVerticeToVertice(lcVertice, rcVertice);
			var horRadian = 2 * Math.asin(horLength / D);
			this.cameraMoveInfo.horDegree = World.Math.radianToDegree(horRadian);
			var tcVertice = earthInOutInfo.tcPickResult[0];
			var bcVertice = earthInOutInfo.bcPickResult[0];
			var verLength = World.Math.getLengthFromVerticeToVertice(tcVertice, bcVertice);
			var verRadian = 2 * Math.asin(verLength / D);
			this.cameraMoveInfo.verDegree = World.Math.radianToDegree(verRadian)
		}
	},
	showBelowLevelSubTiledLayer : function (isCalculateVisubalTiles, geoExtent) {
		if (this.tiledLayer) {
			var subTiledLayer = this.tiledLayer.objectList[World.BELOW_LEVEL];
			subTiledLayer.isvisible = true;
			if (isCalculateVisubalTiles) {
				var visualLevels = this.camera.getVisibleTileIndex(1, geoExtent);
				var visualTiles = visualLevels[1];
				var existTiles = subTiledLayer.objectList;
				for (var i = 0; i < existTiles.length; i++) {
					var existTile = existTiles[i];
					var row = existTile.row;
					var column = existTile.column;
					for (var j = 0; j < visualTiles.length; j++) {
						var visualTileInfo = visualTiles[j];
						if (existTile.row == visualTileInfo.index[1] && existTile.column == visualTileInfo.index[2]) {
							existTile.isvisible = true;
							break
						} else {
							existTile.isvisible = false
						}
					}
				}
			}
		}
	},
	refresh : function (isForceRefresh, geoExtent) {
		if (this.tiledLayer instanceof World.TiledLayer) {
			if (this.scene && this.camera) {
				geoExtent = geoExtent ? geoExtent : this.camera.getCurrentGeoExtent();
				var level = World.CURRENT_LEVEL + 3;
				var idxLevelBelow = World.BELOW_LEVEL;
				var idxLastLevel_1 = level - 1;
				var idxLastLevel = level;
				if (!isForceRefresh && this.tiledLayer.objectList.length == (level + 1)) {
					var position = this.camera.getPosition();
					var isCameraChanged = true;
					if (World.OLD_POSITION) {
						var offset = 1;
						var deltaX = Math.abs(World.OLD_POSITION.x - position.x);
						var deltaY = Math.abs(World.OLD_POSITION.y - position.y);
						var deltaZ = Math.abs(World.OLD_POSITION.z - position.z);
						isCameraChanged = !(deltaX <= offset && deltaY <= offset && deltaZ <= offset)
					} else {
						isCameraChanged = true
					}
					World.OLD_POSITION = position;
					if (!isCameraChanged) {
						var isShowBelow = true;
						var subTiledLayerLastLevel = this.tiledLayer.objectList[idxLastLevel];
						if (subTiledLayerLastLevel.checkIfLoaded()) {
							for (var i = 0; i < idxLastLevel; i++) {
								var subTiledLayerI = this.tiledLayer.objectList[i];
								subTiledLayerI.isvisible = false
							}
							isShowBelow = false
						}
						if (isShowBelow) {
							this.showBelowLevelSubTiledLayer(true, geoExtent)
						}
						return
					}
				}
				var visualLevels = this.camera.getVisibleTileIndex(level, geoExtent);
				var subLayerCount = this.tiledLayer.objectList.length;
				if (subLayerCount < visualLevels.length) {
					var delta = visualLevels.length - subLayerCount;
					for (var i = 0; i < delta; i++) {
						var args = {
							level : i + subLayerCount
						};
						var subTiledLayer = new World.SubTiledLayer(args);
						this.tiledLayer.add(subTiledLayer);
						if (args.level == 1) {
							World.Canvas.style.cursor = "wait";
							for (var m = 0; m <= 1; m++) {
								for (var n = 0; n <= 1; n++) {
									var args1mn = {
										level : 1,
										row : m,
										column : n,
										url : ""
									};
									args1mn.url = this.tiledLayer.getImageUrl(args1mn.level, args1mn.row, args1mn.column);
									var tile1mn = new World.Tile(args1mn);
									subTiledLayer.add(tile1mn)
								}
							}
							World.Canvas.style.cursor = "default"
						}
					}
				} else if (subLayerCount > visualLevels.length) {
					var delta = subLayerCount - visualLevels.length;
					for (var i = 0; i < delta; i++) {
						var subTiledLayer = this.tiledLayer.objectList[visualLevels.length];
						if (subTiledLayer.level != 1) {
							var b = this.tiledLayer.remove(subTiledLayer);
							if (!b) {
								console.debug("LINE:2925,this.tiledLayer.remove(subTiledLayer)失败")
							}
						}
					}
				}
				function checkTileExist(tileArray, lev, row, col) {
					var result = {
						isExist : false,
						index : -1
					};
					for (var m = 0; m < tileArray.length; m++) {
						var tileInfo = tileArray[m];
						var index = tileInfo.index;
						if (lev == index[0] && row == index[1] && col == index[2]) {
							result.isExist = true;
							result.index = m
						}
					}
					return result
				}
				subLayerCount = this.tiledLayer.objectList.length;
				for (var i = 0; i < subLayerCount; i++) {
					var subTiledLayer = this.tiledLayer.objectList[i];
					if (i == World.BELOW_LEVEL || i == subLayerCount - 1 || i == subLayerCount - 2) {
						subTiledLayer.isvisible = true
					} else {
						subTiledLayer.isvisible = false
					}
					if (i == 0) {
						continue
					}
					var needDeletedTiles = [];
					for (var j = 0; j < subTiledLayer.objectList.length; j++) {
						var tile = subTiledLayer.objectList[j];
						var checkResult = checkTileExist(visualLevels[i], tile.level, tile.row, tile.column);
						var isExist = checkResult.isExist;
						if (isExist) {
							visualLevels[i].splice(checkResult.index, 1);
							tile.isvisible = subTiledLayer.isvisible
						} else {
							if (i == 1) {
								tile.isvisible = false
							} else {
								needDeletedTiles.push(tile)
							}
						}
					}
					while (needDeletedTiles.length > 0) {
						var b = subTiledLayer.remove(needDeletedTiles[0]);
						needDeletedTiles.splice(0, 1);
						if (!b) {
							console.debug("LINE:2969,subTiledLayer.remove(needDeletedTiles[0])失败")
						}
					}
					if (subTiledLayer.isvisible) {
						for (var k = 0; k < visualLevels[i].length; k++) {
							var tileInfo = visualLevels[i][k];
							var index = tileInfo.index;
							var args = {
								level : index[0],
								row : index[1],
								column : index[2],
								url : ""
							};
							args.url = this.tiledLayer.getImageUrl(args.level, args.row, args.column);
							var tile = new World.Tile(args);
							subTiledLayer.add(tile)
						}
					}
				}
			}
		}
	}
};
var bMouseDown = false;
var previousX = -1;
var previousY = -1;
function onMouseDown(evt) {
	bMouseDown = true;
	previousX = evt.layerX || evt.offsetX;
	previousY = evt.layerY || evt.offsetY;
	canvas.addEventListener("mousemove", onMouseMove);
	canvas.style.cursor = "move"
}
function onMouseMove(evt) {
	var currentX = evt.layerX || evt.offsetX;
	var currentY = evt.layerY || evt.offsetY;
	globe.showBelowLevelSubTiledLayer();
	if (bMouseDown) {
		onRotateMouseMove(currentX, currentY)
	}
	previousX = currentX;
	previousY = currentY;
	canvas.style.cursor = "move"
}
function onRotateMouseMove(currentX, currentY) {
	if (previousX > 0 && previousY > 0) {
		var changeX = currentX - previousX;
		var changeY = currentY - previousY;
		var changeHorAngle = changeX / World.Canvas.width * globe.cameraMoveInfo.horDegree;
		var changeVerAngle = changeY / World.Canvas.height * globe.cameraMoveInfo.verDegree;
		var rotateRadianHor = -changeHorAngle * Math.PI / 180;
		globe.camera.worldRotateY(rotateRadianHor);
		var lightDir = globe.camera.getLightDirection();
		var plumbVector = getPlumbVector(lightDir, false);
		var rotateRadianVer = -changeVerAngle * Math.PI / 180;
		globe.camera.worldRotateByVector(rotateRadianVer, plumbVector)
	}
}
function onMouseWheel(evt) {
	if (evt.wheelDelta) {
		var delta = evt.wheelDelta;
		var deltaLevel = parseInt(delta / 120);
		var newLevel = World.CURRENT_LEVEL + deltaLevel;
		globe.setLevel(newLevel)
	} else if (evt.detail) {
		var delta = evt.detail;
		var deltaLevel = -parseInt(delta / 3);
		var newLevel = World.CURRENT_LEVEL + deltaLevel;
		globe.setLevel(newLevel)
	}
}
function onMouseUp(evt) {
	bMouseDown = false;
	previousX = -1;
	previousY = -1;
	canvas.removeEventListener("mousemove", onMouseMove);
	canvas.style.cursor = "default"
}
function onDbClick(evt) {
	var absoluteX = evt.layerX || evt.offsetX;
	var absoluteY = evt.layerY || evt.offsetY;
	console.log("absoluteX:" + absoluteX + "," + "absoluteY:" + absoluteY);
	var pickResult = globe.camera.getPickCartesianCoordInEarth(absoluteX, absoluteY);
	if (pickResult.length == 0) {
		console.log("没有单击到地球");
		globe.setLevel(World.CURRENT_LEVEL + 1)
	} else if (pickResult.length >= 1) {
		var pickVertice = pickResult[0];
		var origin = new World.Vertice(0, 0, 0);
		globe.camera.look(pickVertice, origin);
		globe.setLevel(World.CURRENT_LEVEL + 1)
	}
}
function getPlumbVector(direction, bLeft) {
	direction.y = 0;
	direction.normalize();
	var plumbVector = new World.Vector(-direction.z, 0, direction.x);
	plumbVector.normalize();
	return plumbVector
}
var canvas, globe;
function initLayout() {
	canvas.width = document.body.clientWidth;
	canvas.height = document.body.clientHeight;
	if (globe) {
		globe.camera.setAspect(canvas.width / canvas.height);
		globe.calculateCameraMoveInfo();
		globe.refresh(true)
	}
}
function initEvents() {
	window.onresize = initLayout;
	canvas.addEventListener("mousedown", onMouseDown);
	canvas.addEventListener("mouseup", onMouseUp);
	canvas.addEventListener("mousewheel", onMouseWheel);
	canvas.addEventListener("DOMMouseScroll", onMouseWheel);
	canvas.addEventListener("dblclick", onDbClick)
}
function startWebGL() {
	globe = new World.Globe(canvas);
	var mapSelector = document.getElementById("mapSelector");
	mapSelector.onchange = changeTiledLayer;
	changeTiledLayer();
	setTimeout(window.refresh, 500)
}
function changeTiledLayer() {
	var mapSelector = document.getElementById("mapSelector");
	var newTiledLayer = null;
	if (mapSelector.value == "google") {
		newTiledLayer = new World.GoogleTiledLayer()
	} else if (mapSelector.value == "osm") {
		newTiledLayer = new World.OsmTiledLayer()
	} else if (mapSelector.value == "nokia") {
		newTiledLayer = new World.NokiaTiledLayer()
	} else if (mapSelector.value == "bing") {
		newTiledLayer = new World.BingTiledLayer()
	} else if (mapSelector.value == "arcgisonline") {
		var args = {
			url : "http://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer"
		};
		newTiledLayer = new World.ArcGISTiledLayer(args)
	} else if (mapSelector.value == "tianditu") {
		newTiledLayer = new World.TiandituTiledLayer()
	} else if (mapSelector.value == "soso") {
		newTiledLayer = new World.SosoTiledLayer()
	} else if (mapSelector.value == "gistech") {
		var args = {
			url : "http://map.geoq.cn/ArcGIS/rest/services/ChinaOnlineCommunity/MapServer"
		};
		newTiledLayer = new World.ArcGISTiledLayer(args)
	} else if (mapSelector.value == "autonavi") {
		newTiledLayer = new World.AutonaviTiledLayer()
	} else if (mapSelector.value == "blend") {
		newTiledLayer = new World.BlendTiledLayer()
	} else if (mapSelector.value == "test") {
		newTiledLayer = new World.TestTiledLayer()
	}
	if (newTiledLayer) {
		globe.setTiledLayer(newTiledLayer)
	}
}
function initAll() {
	canvas = document.getElementById("canvasId");
	initLayout();
	initEvents();
	startWebGL()
}
function refresh() {
	setTimeout(window.refresh, 750);
	globe.refresh()
}
window.onload = initAll;
