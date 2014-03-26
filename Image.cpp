/*
 * Image.cpp
 *
 *  Created on: Feb 7, 2014
 *      Author: Matthew Reeping
 *
 *      Object to hold all pertinent image
 *      data: background color, from, at,
 *      and eye vectors, etc.
 */

#include "Image.h"

Image::Image() {
	res = new int[2];
	background = new MyVector(3);
	up = new MyVector(3);
	from = new MyVector(3);
	material = new float[8];
	at = new MyVector(3);
	hither = 0;
	angle = 0;
	distance = 1; // Unit vector distance
}

Image::~Image() {
	delete[] res;
	delete[] material;
	delete from;
	delete up;
	delete at;
	delete background;
}

void Image::set_resolution(int x, int y){
	res[0] = x;
	res[1] = y;
}

void Image::set_hither(float hith){
	hither = hith;
}

void Image::set_distance(float dis){
	distance = dis;
}

void Image::set_view(float view){
	angle = view;
}

/*
 * red - red color component
 * green - green color component
 * blue - blue color component
 * Kd - diffuse coefficient
 * Ks - specular coefficient
 * e - exponent for Blinn-Phong
 * T - transmittance
 * ior - index of refraction
 */
void Image::set_material(float red, float green ,float blue, float Kd,
				float Ks, float e, float T, float ior){
	material[0] = red;
	material[1] = green;
	material[2] = blue;
	material[3] = Kd;
	material[4] = Ks;
	material[5] = e;
	material[6] = T;
	material[7] = ior;
}

void Image::set_at(float* vect){
	at->set_vector(vect);
}

void Image::set_from(float* vect){
	from->set_vector(vect);
}

void Image::set_up(float* vect){
	up->set_vector(vect);
}

void Image::set_back_color(float* vect){
	background->set_vector(vect);
}

int Image::get_resy(){
	return res[1];
}

int Image::get_resx(){
	return res[0];
}

float Image::get_view(){
	return angle;
}

float* Image::get_background(){
	return background->get_point();
}

MyVector* Image::get_at(){
	return at;
}

MyVector* Image::get_from(){
	return from;
}

MyVector* Image::get_up(){
	return up;
}

float* Image::get_material(){
	return material;
}

float Image::get_hither(){
	return hither;
}

float Image::get_distance(){
	return distance;
}
