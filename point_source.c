#include <math.h>

#include "point_source.h"

float eval_point_source(float t, float nu){

	// Time shift to ensure causality
	float t0 = 6.0 / (M_PI * nu * sqrt(2.0));
	float ts = (t-t0);
	float v = -1*M_PI*M_PI*nu*nu*ts*ts;
	float w = (1+2*v)*exp(v);

	return w;
}