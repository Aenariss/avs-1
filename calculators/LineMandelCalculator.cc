/**
 * @file LineMandelCalculator.cc
 * @author Vojtech Fiala <xfiala61@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over lines
 * @date DATE
 */
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <mm_malloc.h>


#include <stdlib.h>


#include "LineMandelCalculator.h"
#define BLOCK_SIZE 64

LineMandelCalculator::LineMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "LineMandelCalculator")
{	
	data = (int *) _mm_malloc(height * width * sizeof(int), BLOCK_SIZE);
	reals_precounted = (float *) _mm_malloc(width * sizeof(float), BLOCK_SIZE);
	imags_precounted = (float *) _mm_malloc(half * sizeof(float), BLOCK_SIZE);
	reals = (float *) _mm_malloc(width * sizeof(float), BLOCK_SIZE);
	imags = (float *) _mm_malloc(width * sizeof(float), BLOCK_SIZE);

	// init values
	#pragma omp simd
	for (int i = 0; i < width*height; i++)
		data[i] = 0;
	
	// count only once here, negligent impact on performance but better to count in initialization imo
	for (int i = 0; i < half; i++)
		imags_precounted[i] = y_start + i * dy;

	// count only once here, will be used on many places later
	for (int i = 0; i < width; i++)
		reals_precounted[i] = x_start + i * dx;

}

LineMandelCalculator::~LineMandelCalculator() {
	_mm_free(data);
	_mm_free(reals);
	_mm_free(imags);
	_mm_free(reals_precounted);
	_mm_free(imags_precounted);
	imags = NULL;
	reals_precounted = NULL;
	data = NULL;
	reals = NULL;
	imags_precounted = NULL;
}


int * LineMandelCalculator::calculateMandelbrot () {

	int *pdata = data;
	int times = 1;
	const int mid = height*width/2;
	for (int i = 0; i < half; i++) {

		float y = imags_precounted[i];	// y value is the same for all values in line
		// because im going through the line many times (cuz of limit), i need to fill base values for each value in the line 
		// x is the real precounted base, y is the precounted y
		for (int z = 0; z < width; z++) {
			reals[z] = reals_precounted[z];
			imags[z] = y;
		}
			
		for (int j = 0; j < limit; j++) {
			
			#pragma omp simd
			for (int k = 0; k < width; k++) {
				float x = reals_precounted[k];	// dont count it again, its always the same for each value in line
				float zReal = reals[k];	// real value in line
				float zImag = imags[k];	// img value
				float r2 = zReal * zReal;
				float i2 = zImag * zImag;
				if (r2 + i2 < 4.0f) {	// count new values only if it makes sense, no need to calculate if the result value is not increasing anymore
					pdata[i*width + k]++;	// max achievable value is the limit, cuz it can only be run on each address ``limit`` times.
					reals[k] = r2 - i2 + x;
					imags[k] = 2.0f * zReal * zImag + y;
				}
			}
		}
	}

	// matrix is symmetric, no need to count the other half, just fill it based on the already calculated one
	int z = 0;
	int ct = mid;
	for (int i=0; i < mid; i++) {
       pdata[ct++] = pdata[mid-width*times+z++];
	   if (!(ct % width)) {
			times++;
			z = 0;
	   }
    }
	
	return data;
}
