/*
 * MyVector.cpp
 *
 *  Created on: Feb 9, 2014
 *      Author: Matthew Reeping
 *
 *      Simple vector class geared more toward use in
 *      ray tracer, that hold arrays of floats.
 *      Does contain ability to be multidimensional
 *      and calculate dot and cross(only 3-dimensional) products.
 */

#include "MyVector.h"
#include <math.h>

/**
 * Sets the dimension and initializes vector to all 0.
 */
MyVector::MyVector(int dimension){
	this->dimension = dimension;
	vector = new float[dimension];
	int i;
	for(i=0; i<dimension; i++){
			vector[i] = 0;
	}
	this->intensity = new float[3];
	for(i=0; i<3; i++){
		this->intensity[i] = 1;
	}
}

/*
 * Sets dimension and initializes vector as copy of vect
 */
MyVector::MyVector(int dimension, float* vect) {
	vector = new float[dimension];
	int i;
	for(i=0; i<dimension; i++){
		vector[i] = vect [i];
	}
	this->dimension = dimension;
	this->intensity = new float[3];
	for(i=0; i<3; i++){
		this->intensity[i] = 1;
	}
}

/*
 *
 * Sets dimension and initializes vector as copy of vect
 * intensity as copy of intensity
 *
 */
MyVector::MyVector(int dimension, float* vect, float* intensity){
	vector = new float[dimension];
	int i;
	for(i=0; i<dimension; i++){
		vector[i] = vect [i];
	}
	this->dimension = dimension;
	this->intensity = new float[3];
	for(i=0; i<3; i++){
		this->intensity[i] = intensity [i];
	}
}

/*
 * copy constructor
 */
MyVector::MyVector(MyVector* copy){
	dimension = copy->dimension;
	int i;
	vector = new float[dimension];
	for(i=0; i<dimension; i++){
		vector[i] = copy->vector[i];
	}
	this->intensity = new float[3];
	for(i=0; i<3; i++){
		this->intensity[i] = copy->intensity [i];
	}

}

/*
 * destructor
 */
MyVector::~MyVector(){
	delete [] vector;
	delete [] intensity;
}

void MyVector::set_vector(float* vect){
	int i;
	for(i=0; i<dimension; i++){
		vector[i] = vect[i];
	}
}

// accessors
int MyVector::get_dimension(){ return dimension; }
float* MyVector::get_point(){ return vector; }
float* MyVector::get_intensity(){ return intensity; }

/*
 * adds two vectors together
 */
MyVector* MyVector::add(MyVector* b){
	MyVector* temp = new MyVector(dimension);
	int i;
	for(i=0; i<dimension; i++){
		temp->vector[i] = vector[i] + b->vector[i];
	}
	return temp;
}

/*
 * subtracts 'b' from vector
 */
MyVector* MyVector::subtract(MyVector* b){
	MyVector* temp = new MyVector(dimension);
	int i;
	for(i=0; i<dimension; i++){
		temp->vector[i] = vector[i] - b->vector[i];
	}
	return temp;
}

/*
 * scales vector by some constant
 */
void MyVector::scalar_multiply(float scale){
	int i;
	for(i=0; i<dimension; i++){
		vector[i] *= scale;
	}
}

/*
 * gives vector unit length
 */
float MyVector::normalize(){
	float temp=0;
	int i;
	for(i=0; i<dimension; i++){
		temp += pow(vector[i], 2);
	}
	return sqrt(temp);
}

/*
 * Only done for 3 dimensional vectors
 */
MyVector* MyVector::cross_product(MyVector* b){
	MyVector* temp = new MyVector(3);
	temp->vector[0] = vector[1]*b->vector[2] - vector[2]*b->vector[1]; // aybz - azby
	temp->vector[1] = vector[2]*b->vector[0] - vector[0]*b->vector[2];
	temp->vector[2] = vector[0]*b->vector[1] - vector[1]*b->vector[0];
	return temp;
}

/*
 * dot product of vector
 */
float MyVector::dot_product(MyVector* b){
	float temp=0;
	int i;
	for(i=0; i<dimension; i++){
		temp += vector[i]*b->vector[i];
	}
	return temp;
}
