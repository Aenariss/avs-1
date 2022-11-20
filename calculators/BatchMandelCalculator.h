/**
 * @file BatchMandelCalculator.h
 * @author Vojtech Fiala <xfiala61@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over small batches
 * @date 18.11.2022
 */
#ifndef BATCHMANDELCALCULATOR_H
#define BATCHMANDELCALCULATOR_H

#include <BaseMandelCalculator.h>

class BatchMandelCalculator : public BaseMandelCalculator
{
public:
    BatchMandelCalculator(unsigned matrixBaseSize, unsigned limit);
    ~BatchMandelCalculator();
    int * calculateMandelbrot();

private:
    const int half = (int) height/2;
    int *data;
    float *reals;
    float *imags;
    float *reals_precounted;
    float *imags_precounted;
};

#endif