### Point-in-Polygon

*The algorithm telling if a point is inside or outside of the give convex polygon.*

**Steps:**
	
	1) Choose a point inside of a polygon: 
		
		centroid of a tiangle made of any 3 vertices

	2) Order the vertices by their angle relative to the chosen point:
		
		```atan2((vertex - center).x, (vertex - center).y)```

	3) Using logarithmic-complexity search find a wedge the located point belongs to:
		
		```map<float, Vector3f> vertices;  vertices.upper_bound(angle);```

	4) Determine if the point is in internal part of the wedge or not.