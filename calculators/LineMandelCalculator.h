/**
 * @file LineMandelCalculator.h
 * @author Vojtech Fiala <xfiala61@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over lines
 * @date 18.11.2022
 */

#include <BaseMandelCalculator.h>

class LineMandelCalculator : public BaseMandelCalculator
{
public:
    LineMandelCalculator(unsigned matrixBaseSize, unsigned limit);
    ~LineMandelCalculator();
    int *calculateMandelbrot();

private:
    const int half = (int) height/2;
    int *data;
    float *reals;
    float *imags;
    float *reals_precounted;
    float *imags_precounted;
};