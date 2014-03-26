/*
 * Image.h
 *
 *  Created on: Feb 7, 2014
 *      Author: Matthew Reeping
 */

#ifndef IMAGE_H_
#define IMAGE_H_

#include "MyVector.h"

class Image {

private:

	float distance; // distance from eye to image plane
	int* res; // x then y
	float angle; // Field of view
	MyVector* background; // r g b
	MyVector* up; // x,y,z
	MyVector* at; // Where the camera is pointing to (x,y,z)
	float hither;
	MyVector* from; // Camera position(x,y,z)
	float* material; // Object material (r g b, Kd, Ks, e, T, Index_of_Refraction)

public:
	Image();
	virtual ~Image();

	// setters
	void set_resolution(int x, int y);
	void set_back_color(float* vect);
	void set_material(float red, float green, float blue, float Kd, float Ks, float e, float T, float ior);
	void set_from(float* vect);
	void set_at(float* vect);
	void set_up(float* vect);
	void set_view(float view);
	void set_hither(float hith);
	void set_distance(float dis);

	// accessors
	int get_resx();
	int get_resy();
	float get_view();
	float* get_background();
	MyVector* get_at();
	MyVector* get_up();
	float get_hither();
	MyVector* get_from();
	float* get_material();
	float get_distance();


};

#endif /* IMAGE_H_ */
