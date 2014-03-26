/*
 * Sphere.cpp
 *
 *  Created on: Feb 9, 2014
 *      Author: Matthew Reeping
 *
 *      Sphere object used in ray tracer. Holds position
 *      in world space, color conponents (RGB, diffuse and
 *      specular coeff, e for specular), index of refraction,
 *      and T.
 */

#include "Sphere.h"

Sphere::Sphere(float x, float y, float z, float r, float* color,
				float Kd, float Ks, float e, float T, float ir){
	float* vect = new float[3];
	vect[0]=x;
	vect[1] = y;
	vect[2] = z;
	center = new MyVector(3, vect);
	radius = r;
	this->Kd = Kd;
	this->Ks = Ks;
	this->T = T;
	this->e = e;
	this->ir = ir;
	this->color = new float[3];
	this->color[0] = color[0];
	this->color[1] = color[1];
	this->color[2] = color[2];
	delete[] vect;
}

Sphere::~Sphere() {
	delete center;
	delete color;
}

// accessors
MyVector* Sphere::get_center(){ return center; }
float Sphere::get_radius(){ return radius; }
float* Sphere::get_color(){ return color; }
float Sphere::get_Ks(){ return Ks; }
float Sphere::get_Kd(){ return Kd; }
float Sphere::get_Shine(){ return e; }
float Sphere::get_IndexofRefraction(){ return ir; }
float Sphere::get_T(){ return T; }

/*
 * calculate 't' to determine intersection of ray and sphere.
 * d - directional vector
 * eye - origin of view vector
 * hither - closest object can be
 */
float Sphere::calc_intersection(MyVector* d, MyVector* eye, float hither){
	MyVector* ec = eye->subtract(center);
	float disc = discriminate(d, ec);
	float t=0;
	if(disc < 0){ // if negative then no intersection
		delete ec;
		return -10;
	}

	// make calculations to solve for t
	else {
		disc = sqrt(disc);
		float closest = 10000000;
		d->scalar_multiply(-1);
		float dec = d->dot_product(ec);
		float dd = d->dot_product(d);
		t += dec;
		t += disc;
		t = t/dd;
		if(t<closest && t>hither){
			closest = t;
		}
		t=0;
		t += dec;
		t -= disc;
		t = t/dd;
		if(t<closest && t>hither){
			closest = t;
		}
		d->scalar_multiply(-1);
		delete ec;
		return closest;		// return min t
	}
}

/*
 * calculate the discriminate for intersection equation of sphere and ray.
 * if negative, there is no intersection.
 */
float Sphere::discriminate(MyVector* d, MyVector* ec){
	float dec = d->dot_product(ec);
	dec = pow(dec, 2);
	float ecr = ec->dot_product(ec);
	ecr -= pow(radius, 2);
	float ddecr = d->dot_product(d);
	ddecr *= ecr;
	return dec - ddecr;
}

