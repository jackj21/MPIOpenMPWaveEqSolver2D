/*
* Authors:
*   Russell J. Hewett (rhewett@vt.edu)
*
* Not licensed for external distribution.
*/

#ifndef __ISING_H__
#define __ISING_H__

#include "array_2d.h"

typedef struct WaveSimParameters_tag {

    float T;
    int N;
    int Mx;
    int My;
    float alpha;

} WaveSimParameters;

int error_norm(Array2D_f* u1, Array2D_f* u2, float* err);

int evaluate_standing_wave(Array2D_f* u, unsigned int Mx, unsigned int My, float dx, float dy, float t);

int wave_timestep(Array2D_f* u_prev, Array2D_f* u_curr, Array2D_f* u_next, float dt, float dx);

int standing_wave_simulation(int nt, unsigned int N, unsigned int Mx, unsigned int My, float alpha);

int standing_wave_simulation_nsteps(Array2D_f* u_prev, Array2D_f* u_curr, Array2D_f* u_next, float dt, float dx, int nt);


#endif
