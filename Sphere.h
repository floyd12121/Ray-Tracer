/*
 * Sphere.h
 *
 *  Created on: Feb 9, 2014
 *      Author: Matthew Reeping
 */

#ifndef SPHERE_H_
#define SPHERE_H_

#include "MyVector.h"
#include <math.h>

class Sphere {

private:
	float radius;
	float Kd;
	float Ks;
	float e;
	float T;
	float ir;
	MyVector* center;
	float* color;

	/*
	 * calculates the discrimnate for intersection equation of
	 * a sphere (stuff under the square root).
	 */
	float discriminate(MyVector* d, MyVector* ec);

public:
	Sphere(float x, float y, float z, float r, float* color, float Kd,
			float Ks, float e, float T, float ir);
	virtual ~Sphere();

	/**
	 * Calculates the intersection of the ray with "this" sphere object
	 */
	float calc_intersection(MyVector* d, MyVector* eye, float hither);

	// accessors
	MyVector* get_center();
	float get_radius();
	float* get_color();
	float get_Ks();
	float get_Kd();
	float get_Shine();
	float get_IndexofRefraction();
	float get_T();
};

#endif /* SPHERE_H_ */
