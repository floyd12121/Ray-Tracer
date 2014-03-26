/*
 * MyVector.h
 *
 *  Created on: Feb 9, 2014
 *      Author: Matthew Reeping
 */

#ifndef MYVECTOR_H_
#define MYVECTOR_H_

class MyVector {
private:
	float* vector;
	int dimension;
	float* intensity;

public:
	// initializers and destructor
	MyVector(int dimension, float* vect);
	MyVector(int dimension);
	MyVector(MyVector* copy); // copy constructor
	virtual ~MyVector();
	MyVector(int dimension, float* vect, float* intensity);

	// accessors
	int get_dimension();
	float* get_point();
	float* get_intensity();

	/*
	 * Set the vector. vect must be size of dimension
	 */
	void set_vector(float* vect);

	/*
	 * calculates the dot product with vector 'b'
	 */
	float dot_product(MyVector* b);

	/*
	 * calculates the cross product with vector 'b'
	 */
	MyVector* cross_product(MyVector* b);

	/*
	 * scales vector by some constant "scale"
	 */
	void scalar_multiply(float scale);

	/*
	 * add with vector 'b'
	 */
	MyVector* add(MyVector* b);

	/*
	 * subtract 'b' from vector
	 */
	MyVector* subtract(MyVector* b);

	/*
	 * return unit length of vector
	 */
	float normalize();

};

#endif /* MYVECTOR_H_ */
