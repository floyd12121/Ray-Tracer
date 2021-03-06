Ray-Tracer
==========

Project: Ray Tracer
Last Modified: March 25, 2014
Author: Matthew Reeping

This is a simple ray tracer program that implements sphere primitives with
reflection, refraction, sphere intersection, lighting and shadows. 
Image data created using ".nff" files and output files with a ".ppm" extension.

Description:
	As stated above this ray tracer is simple and implements only the sphere
	primitives along with other visual feature previously stated. 
	If future time permits, polygons and cones will be implemented. This program parses 
	data from ".nff" files and places data into corresponding data structures for
	further calculations of ray-sphere intersection and pixel colors (RGB) involving
	shadowing, reflection, refraction, and Blinn-Phong effects. After pixel colors are
	calculated, pixels are written to a ".ppm" file as binary data. In order to 
	view the image, the ".ppm" file must be converted to an image file of choice 
	(.jpeg, .png, .gif, etc.). 
	
	This program uses a custom makefile. Follow instructions below on how to run properly.
	
		To run program:
			1.	run command "make" (without quotes) in command prompt to build project.
			2.	then run executable "trace" (without quotes) with .nff file "filename.nff"
			(with quotes).
			
				Ex: make
					trace "balls2.nff"
					
		After project is run, trace.ppm will be the output file, unless change within source code.
		
	Example Code: 
		There is an example .nff and an image (.png) file which demonstrates reflection, lighting
		and shadows. To run:
				trace "balls3.nff"
				
		This may take a few minutes to run. Image file is balls3.png.
		
		Another example demonstrates refraction, reflection, lighting, shadows, and viewing
		from within a larger sphere:
				trace "test.nff"
		
		Image file is test.png.
		
	Known Bugs:
		-	When using refraction, colors are not always 100% correct. Don't believe to be overflow
			because program handles that case, still under debugging process. 
