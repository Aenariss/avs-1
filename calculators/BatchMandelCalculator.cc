/**
 * @file BatchMandelCalculator.cc
 * @author Vojtech Fiala <xfiala61@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over small batches
 * @date DATE
 */

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <mm_malloc.h>

#include <stdlib.h>
#include <stdexcept>

#include "BatchMandelCalculator.h"
#define BLOCK_SIZE 64

BatchMandelCalculator::BatchMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "BatchMandelCalculator")
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


BatchMandelCalculator::~BatchMandelCalculator() {
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


int * BatchMandelCalculator::calculateMandelbrot () {
	int *pdata = data;
	float *preals_precounted = reals_precounted;
	float *pimags_precounted = imags_precounted;
	float *preals = reals;
	float *pimags = imags;
	int times = 1;
	const int mid = height*width/2;

	const int block_divide = (int) width/BLOCK_SIZE;

	for (int i = 0; i < half; i++) {

		float y = pimags_precounted[i];	// y value is the same for all values in line
		// because im going through the line many times (cuz of limit), i need to fill base values for each value in the line 
		// x is the real precounted base, y is the precounted y
		for (int z = 0; z < width; z++) {
			reals[z] = preals_precounted[z];
			imags[z] = y;
		}
		
		for (int q = 0; q < block_divide; q++) {
			for (int j = 0; j < limit; j++) {
				// divide into blocks sized 64 cuz of max size
				int q_blocksize = q*BLOCK_SIZE;
				#pragma omp simd
				for (int k = 0; k < BLOCK_SIZE; k++) {
					int index = k+q_blocksize;
					float x = preals_precounted[index];	// dont count it again, its always the same for each value in line
					float zReal = preals[index];	// real value in line
					float zImag = pimags[index];	// img value
					float r2 = zReal * zReal;
					float i2 = zImag * zImag;
					if (r2 + i2 < 4.0f) {	// count new values only if it makes sense, no need to calculate if the result value is not increasing anymore
						pdata[i*width + index]++;	// max achievable value is the limit, cuz it can only be run on each address ``limit`` times.
						preals[index] = r2 - i2 + x;
						pimags[index] = 2.0f * zReal * zImag + y;
					}
				}
			}
		}
	}

	// matrix is symmetric, no need to count the other half, just fill it based on the already calculated one
	int z = 0;
	int ct = mid;
	for (int i=mid; i < mid*2; i++) {
       pdata[i] = pdata[mid-width*times+z++];
	   if (!(++ct % width)) {
			times++;
			z = 0;
	   }
    }
	
	return data;
}