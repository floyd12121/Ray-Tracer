/*
 * main.cpp
 *
 *  Created on: Feb 3, 2014
 *  Last edited: March 25, 2014
 *      Author: Matthew Reeping
 *
 *		More details can be found in about.txt.
 *      Description: Ray Tracing program. This is a simple ray tracer which
 *    	calculates Sphere intersection within a given world space. program
 *    	calculates world space, viewing plane, etc. from information
 *    	parsed in the given .nff file. It calculates sphere intersection,
 *   	shadowing, and object reflection and refraction... The output is
 *    	given in binary data contained in a .ppm file, which then can be
 *    	converted into any image file extension (.jpeg, .png, .gif, etc.).
 *    	Images can only be in square resolutions (ex: 512 x 512 pixels)
 *
 *    	Project is built from a custom Makefile. Just enter make into command
 *    	prompt to build project, which will give an executable called trace.
 */

#include <iostream>
#include <fstream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include "Image.h"
#include "Sphere.h"
using namespace std;

#define PI 3.14159265	// Used to convert deg into rad and visa versa


void parse(char* filename, Image* image, vector<Sphere*>* s_list, vector<MyVector*>* l_list);
float* get_nextfloats(char* ptr, float* temp, int num_want);
float calc_tan_deg(float theta, float distance);
MyVector* pixel_position(MyVector* u, MyVector* v, MyVector* eye, MyVector* wd_norm,
		float* pix_ptr, int i, int j);
float* calc_pixcolor(MyVector* light, MyVector* normal, MyVector* H, float Ks,
		float Kd, float* color, float e, float* intens);
Sphere* closest_Sphere(vector<Sphere*>* list, MyVector* from, MyVector* to,
		float bias, float& closest);
MyVector* get_reflection(MyVector* v, MyVector* n);
void raycolor( MyVector* view_vector, MyVector* from,  vector<Sphere*>* slist, vector<MyVector*>* lights,
		float bias,  float* color, float* background, float contribution, int steps);



int main(int argc, char* argv[]){
	vector<Sphere*>* sptr = new vector<Sphere*>(); // list of spheres
	vector<MyVector*>* lights = new vector<MyVector*>(); // list of lights
	Image* im = new Image();		// Holds image details

	// parse file and create objects
	parse(argv[1], im, sptr, lights);

	int i, j;

	MyVector* w_norm = im->get_from()->subtract(im->get_at()); // Temporarily the d vector
	im->set_distance(w_norm->normalize());

	// Pixel location calculations
	float*** pix_ptr = new float**[im->get_resy()];
	for(j=0; j<im->get_resy(); j++){
		pix_ptr[j] = new float*[im->get_resx()];
		for(i=0; i<im->get_resx(); i++){
			pix_ptr[j][i] = new float[2];
		}
	}
	float HEIGHT = 2*calc_tan_deg(im->get_view()/2, im->get_distance()); // divide by 2 for u
	float WIDTH = HEIGHT; // divide by two for r
	for(j=0; j<im->get_resy(); j++){
		for(i=0; i<im->get_resx(); i++){ // pixel scales will be in (x, y) order
			pix_ptr[j][i][0] = (WIDTH*(i+0.5)/im->get_resx()) - WIDTH/2; // Us
			pix_ptr[j][i][1] = (HEIGHT*(j+0.5)/im->get_resy()) - HEIGHT/2; // Vs
		}
	}

	unsigned char pixels[im->get_resy()][im->get_resx()][3]; // 3D array holding pixel colors

	// Normalize w vector
	w_norm->scalar_multiply(1/w_norm->normalize());

	// Calculate normal vectors for u and v
	MyVector* u_norm = im->get_up()->cross_product(w_norm);
	u_norm->scalar_multiply(1/u_norm->normalize());
	MyVector* v_norm = w_norm->cross_product(u_norm);
	v_norm->scalar_multiply(-1/v_norm->normalize()); // Had to negate to get correct output, flipped vertically
	MyVector* wd_norm = new MyVector(w_norm);
	wd_norm->scalar_multiply(im->get_distance());

	// Loop for image creation
	for(j=0; j<im->get_resy(); j++){
		for(i=0; i<im->get_resx(); i++){
			MyVector* s = pixel_position(u_norm, v_norm, im->get_from(), wd_norm, pix_ptr[j][i], i, j);
			MyVector* temp_d = s->subtract(im->get_from()); // get's vector d = (s-e)

			float color[3] = {0, 0, 0};		// Set color as background
			int k=0;
			raycolor(temp_d, im->get_from(), sptr, lights, im->get_hither(),
					color, im->get_background(), 1, 0);

			// Set pixel color
			for(k=0; k<3; k++){
				float temp =color[k]*255;
				if(temp>255)
					temp = 255;
				else if(temp<0)
					temp =0;
				pixels[j][i][k] = temp;
			}
			delete s; // pixel_position() returns a new instance of object
			delete temp_d;
		}
	}

	// Printing binary color code to .ppm file, P6, of size (x, y)
	FILE *f = fopen("trace.ppm","wb");
	fprintf(f, "P6\n%d %d\n%d\n", im->get_resx(), im->get_resy(), 255);
	fwrite(pixels, 1, im->get_resx()*im->get_resy()*3, f);
	fclose(f);

	// Delete the array's of pixel locations: 3D array
	for(j=0; j<im->get_resy(); j++){
		for(i=0; i<im->get_resx(); i++){
			delete[] pix_ptr[j][i];
		}
	}
	for(j=0; j<im->get_resy(); j++){
		delete[] pix_ptr[j];
	}
	delete[] pix_ptr;
	delete wd_norm;
	delete w_norm;
	delete u_norm;
	delete v_norm;
	delete sptr;
	delete lights;

	return 0;		// exit program with no issues
}

/*
 * parses a .nff file and store data in its corresponding objects.
 * Does not use polygons, cones, or polygon patch primitives.
 * filename - input .nff file name to be parsed
 * image - data object for certain types of data (background color, vectors, etc.)
 * s_list - empty list for spheres
 * l_list - empty list for light components
 */
void parse(char* filename, Image* image, vector<Sphere*>* s_list, vector<MyVector*>* l_list){
	ifstream nfile(filename, ios::in);
	string line;
	if(nfile.is_open()){ 							   	// opens the file to be parsed
		while ( getline(nfile, line)){				    // accesses a single line from the file
		  char * cline = new char[line.length()+1];
		  strcpy(cline, line.c_str());
		  char* ptr = strtok(cline, " ");

		  // Conditionals for certain types of input data
		  if(!strcmp(ptr, "b")){ 			// Background color data
		    float temp[3];
		    get_nextfloats(ptr, temp, 3);
		    image->set_back_color(temp);
		  }
		  else if(!strcmp(ptr, "from")){ 	// View Data
		    float temp[3];
		    get_nextfloats(ptr, temp, 3);
		    image->set_from(temp);
		  }
		  else if(!strcmp(ptr, "at")){
		    float temp[3];
		    get_nextfloats(ptr, temp, 3);
		    image->set_at(temp);
		  }
		  else if(!strcmp(ptr, "angle")){
		    float temp[1];
		    get_nextfloats(ptr, temp, 1);
		    image->set_view(*temp);
		  }
		  else if(!strcmp(ptr, "hither")){
		    float temp;
		    get_nextfloats(ptr, &temp, 1);
		    image->set_hither(temp);
		  }
		  else if(!strcmp(ptr, "up")){
		    float temp[3];
		    get_nextfloats(ptr, temp, 3);
		    image->set_up(temp);
		  }
		  else if(!strcmp(ptr, "resolution")){
		    float temp[2];
		    get_nextfloats(ptr, temp, 2);
		    image->set_resolution((int)temp[0], (int)temp[1]);
		  }
		  else if(!strcmp(ptr, "l")){ 		// Positional Light Data
			  float temp[6];
			  get_nextfloats(ptr, temp, 6);
			  float templ[3] = {temp[3], temp[4], temp[5]};
			  l_list->push_back(new MyVector(3, temp, templ));
		  }
		  else if(!strcmp(ptr, "f")){ 		// Fill color and shading
		    float temp[3];
		    get_nextfloats(ptr, temp, 8);
		    image->set_material(temp[0], temp[1], temp[2], temp[3],
		    					temp[4], temp[5], temp[6], temp[7]);
		  }
		  else if(!strcmp(ptr, "s")){ 		// Sphere data
		    float temp [4];
		    get_nextfloats(ptr, temp, 4);
		    float* mat = image->get_material();
		    Sphere* s = new Sphere(temp[0], temp[1], temp[2], temp[3],
		    					mat, mat[3], mat[4], mat[5], mat[6], mat[7]);
		    s_list->push_back(s);
		  }
		  else if(!strcmp(ptr, "p")){ 		// Polygon data
		    // skip over certain number of lines
		    // Ex: p 4, skip over the next 4 lines in .nff file
		    ptr = strtok(NULL, " ");
		    int skip = *ptr -'0';
		    int i;
		    for(i=0; i<skip; i++){
		      getline(nfile, line);
		    }
		  }
		  else { /*Do nothing */}
		  delete [] cline; // clear the memory of cline for new allocation

		}
		nfile.close();
	}
	else{
		cerr << "There is an error opening " << filename << endl;
		nfile.close();
		exit(EXIT_FAILURE);

	}
}

/**
 * Parses the next desired number of floats from a line in the ifstream.
 * ptr - current buffer of line in ifstream to be parsed.
 * temp - array of floats to hold new values
 * num_want - number of float values desired from ptr
 */
float* get_nextfloats(char* ptr, float* temp, int num_want){
	ptr = strtok(NULL, " "); // Access first value
	int i;
	for(i=0; i<num_want; i++){
		if(ptr!=NULL)
			temp[i] = atof(ptr);
		else
			temp[i] = 1;
		ptr = strtok(NULL, " ");
	}
	return temp;
}

/*
 * calculate the tangent of theta (in degrees) and
 * return value of distance times tan(theta) in radians.
 * helps with scaling image to proper size.
 */
float calc_tan_deg(float theta, float distance){
	float temp = tan(theta * PI /180);
	return temp*distance;
}

/*
 * Caculates the s vector for pixel position in camera view
 */
MyVector* pixel_position(MyVector* u, MyVector* v, MyVector* eye, MyVector* wd_norm,
							float* temp, int i, int j){
	MyVector* s = new MyVector(eye);
	MyVector* temp_s;
	MyVector* temp_u = new MyVector(u);
	MyVector* temp_v = new MyVector(v);
	// Compute position for each pixel
	temp_u->scalar_multiply(temp[0]); // scaling u vector to pixel location
	temp_v->scalar_multiply(temp[1]); // scaling v vector to pixel location
	// Find a better way to calculate using vectors and not having to delete after every operation
	temp_s = s->add(temp_u);
	delete s;
	s = temp_s->add(temp_v);
	delete temp_s;
	temp_s = s->subtract(wd_norm);
	delete s;
	delete temp_u;
	delete temp_v;
	return temp_s;
}

/*
 * Calculate color based on diffuse and specular components.
 * light - light position
 * normal - normal to surface who's surface color is being calculated
 * H - half-way vector between view and reflection
 * Ks - specular coefficient
 * Kd - diffuse coefficient
 * color - RGB array of color
 * e - specular exponent
 * intens - light intensity in RGB
 */
float* calc_pixcolor(MyVector* light, MyVector* normal, MyVector* H, float Ks, float Kd,
						float* color, float e, float* intens){
	float* temp_color = new float[3];
	float diffuse = max((float)0, normal->dot_product(light));
	float specular = pow(max((float)0, H->dot_product(normal)), e);

	// calculate each color
	temp_color[0] = (Kd*color[0] + Ks*specular) * diffuse * intens[0];	// red
	temp_color[1] = (Kd*color[1] + Ks*specular) * diffuse * intens[1];	// blue
	temp_color[2] = (Kd*color[2] + Ks*specular) * diffuse * intens[2];	// green

	return temp_color;
}

/**
 * Determines the closest sphere from the given list of spheres
 * based on an origin and a destination.
 * Precondition: Closest must be set as a very large number for intersection to work properly
 * list - list of spheres to test intersection
 * from - origin vector
 * to - destination vector
 * bias - smallest value t can be
 * closest - the current closest sphere's distance
 */
Sphere* closest_Sphere(vector<Sphere*>* list, MyVector* from, MyVector* to, float bias, float& closest){
	int k;
	float t=0;
	Sphere* closest_s = NULL;
	for(k=0; k<list->size(); k++){
		Sphere* temp_s = list->at(k); // Current sphere being tested for intersection
		t = temp_s->calc_intersection(to, from, bias);
		if( t< closest && t>bias){
			closest = t;
			closest_s = temp_s;
		}
	}
	return closest_s;			// return closest sphere
}

/*
 * If contribution is < 1/255 or steps >= 5 then end recursion
 * v - view vector in opposite direction
 * p - point on the sphere from which the reflection is coming from
 * n - normal to the surface reflecting
 * bias - the bias used to keep object from intersecting themselves
 */
MyVector* get_reflection(MyVector* v, MyVector* n){
	// Calculate reflection vector
	MyVector* temp_r = new MyVector(n);
	temp_r->scalar_multiply(2*v->dot_product(n));
	MyVector* reflect = temp_r->subtract(v);
	delete temp_r; // free memory of temp_r
	return reflect;
}

/*
 * calculate the refracted ray.
 * v - incoming vector
 * n - normal to surface of object
 * nt - index of refraction of object entering
 * nv - index of refraction of object exiting
 * in - true if going into object, otherwise false
 */
MyVector* get_refraction(MyVector* v, MyVector* n, float nt, float nv, bool in){
	// Calculate ncos(theta t)
	float sqrtstuff = 1-pow(n->dot_product(v),2); // (1 - (n dot v)^2)
	sqrtstuff = sqrtstuff*(pow(nv,2)/pow(nt,2)); // above * (nv^2)/(nt^2)
	sqrtstuff = 1 - sqrtstuff; // sqr(1- above)
	if(sqrtstuff < 0)
		return NULL;
	MyVector* ncos = new MyVector(n);

	if(in)
		ncos->scalar_multiply(-sqrt(sqrtstuff));		// entering object
	else
		ncos->scalar_multiply(sqrt(sqrtstuff));			// exiting object

	// Calculate msin(theta t)
	MyVector* temp_m = new MyVector(n);
	temp_m->scalar_multiply(n->dot_product(v)); //(n dot v)
	MyVector* msin = temp_m->subtract(v); // (n(n dot v) - v)
	msin->scalar_multiply(nv/nt);

	MyVector* t = msin->subtract(ncos); // t = msin - ncos

	delete temp_m;
	delete ncos;
	delete msin;
	return t;
}

/*
 * Recursive function to calculate the given ray's color
 * view_vector - view vector
 * from - origin of view vector
 * slist - list of spheres in world
 * lights - list of lights in world
 * bias - lowest values for intersection
 * color - current ray color
 * background - background color
 * contribution - the ratio the given material can contribute to ray's total color (min = 1/255)
 * steps - number of steps deep recursion is (max = 5)
 */
void raycolor(MyVector* view_vector, MyVector* from,  vector<Sphere*>* slist, vector<MyVector*>* lights,
		float bias, float* color, float* background, float contribution, int steps){

	if(contribution<(1/255) || steps>5){
		return;
	}
	float closest = 10000000; // sets the closest object in scene, initialized far away
	Sphere* closest_s = closest_Sphere(slist, from, view_vector, bias, closest);
	int k;
	float* new_color;
	if(closest_s != NULL){
		// Calculate normal to intersection
		MyVector* temp_p = new MyVector(view_vector);
		temp_p->scalar_multiply(closest);
		MyVector* p = from->add(temp_p); // Point of intersection
		delete temp_p; // Free memory of temp_p
		MyVector* normal = p->subtract(closest_s->get_center()); // Normal to intersection
		normal->scalar_multiply(1/normal->normalize());
		if(closest_s->get_radius()<0){
			normal->scalar_multiply(-1);
		}

		// Calculate color based on each light
		for(k=0; k<lights->size(); k++){
			MyVector* light = lights->at(k)->subtract(p);
			light->scalar_multiply(1/light->normalize());
			int kk;
			int intersect = 0;

			// Shadow ray for object intersection
			for(kk=0; kk<slist->size(); kk++){
				float t2 = slist->at(kk)->calc_intersection(light, p, bias);
				if(t2>=bias && t2 <= 1){ // If within this range, then shadow
					intersect = 10;
				}
			}
			if(intersect==0){	 // No object intersection of shadow ray
				MyVector* H = light->subtract(view_vector);
				H->scalar_multiply(1/H->normalize()); // H vector for specular reflection
				new_color = calc_pixcolor(light, normal, H, closest_s->get_Ks(), closest_s->get_Kd(),
								closest_s->get_color(), closest_s->get_Shine(), light->get_intensity());
				for(kk=0; kk<3; kk++){
					 // add new color based on each light
					color[kk] += (new_color[kk]/sqrt((double)lights->size()))*contribution;
				}
				delete H;
				delete[] new_color;
				delete light;
			}
		}

		steps++;
		view_vector->scalar_multiply(-1); // view coming on from opposite direction

		// if reflective calculate reflection ray color
		if(closest_s->get_Ks()>0){
			MyVector* reflect = get_reflection(view_vector, normal);
			raycolor(reflect, p, slist, lights, bias, color, background,
						contribution*closest_s->get_Ks(), steps);
			delete reflect;
		}

		// if refractive calculate refraction ray color
		if(closest_s->get_T()>0){
			float air_ior = 1;
			MyVector* refraction;
			if(view_vector->dot_product(normal)>0) // air to sphere
				refraction = get_refraction(view_vector, normal, closest_s->get_IndexofRefraction(),
								air_ior, false);
			else // Sphere to air refraction
				refraction = get_refraction(view_vector, normal, air_ior,
								closest_s->get_IndexofRefraction(), true);
			if(refraction != NULL){
				raycolor(refraction, p, slist, lights, bias, color, background,
							contribution*closest_s->get_T(), steps);
				delete refraction;
			}
		}
		view_vector->scalar_multiply(-1); // flipping view back
		delete p;
		delete normal;
	}
	else{ // If no intersection then set to background color
		for(k=0; k<3; k++){
			color[k] += background[k]*contribution; // add new color based on each light
		}
	}
}
