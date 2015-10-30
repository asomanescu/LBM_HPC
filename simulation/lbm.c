/*
** code to implement a d2q9-bgk lattice boltzmann scheme.
** 'd2' inidates a 2-dimensional grid, and
** 'q9' indicates 9 velocities per grid cell.
** 'bgk' refers to the bhatnagar-gross-krook collision step.
**
** the 'speeds' in each cell are numbered as follows:
**
** 6 2 5
**  \|/
** 3-0-1
**  /|\
** 7 4 8
**
** a 2d grid:
**
**           cols
**       --- --- ---
**      | d | e | f |
** rows  --- --- ---
**      | a | b | c |
**       --- --- ---
**
** 'unwrapped' in row major order to give a 1d array:
**
**  --- --- --- --- --- ---
** | a | b | c | d | e | f |
**  --- --- --- --- --- ---
**
** grid indicies are:
**
**          ny
**          ^       cols(jj)
**          |  ----- ----- -----
**          | | ... | ... | etc |
**          |  ----- ----- -----
** rows(ii) | | 1,0 | 1,1 | 1,2 |
**          |  ----- ----- -----
**          | | 0,0 | 0,1 | 0,2 |
**          |  ----- ----- -----
**          ----------------------> nx
**
** note the names of the input parameter and obstacle files
** are passed on the command line, e.g.:
**
**   d2q9-bgk.exe input.params obstacles.dat
**
** be sure to adjust the grid dimensions in the parameter file
** if you choose a different obstacle file.
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <omp.h>

#include "lbm.h"

/*
** main program:
** initialise, timestep loop, finalise
*/

void swap_cells(float **c, float **t) {
    float *tmp = *c;
    *c = *t;
    *t = tmp;
}
int main(int argc, char* argv[])
{
    char * final_state_file = NULL;
    char * av_vels_file = NULL;
    char * param_file = NULL;

    accel_area_t accel_area;

    param_t  params;              /* struct to hold parameter values */
    unsigned char*     obstacles = NULL;    /* grid indicating which cells are blocked */
    float*  av_vels   = NULL;    /* a record of the av. velocity computed for each timestep */
    float *cells;
    float *tmp_cells;
    unsigned total_cells;
    unsigned    i;                    /*  generic counter */
    struct timeval timstr;        /* structure to hold elapsed time */
    struct rusage ru;             /* structure to hold CPU time--system and user */
    double tic,toc;               /* doubleing point numbers to calculate elapsed wallclock time */
    double usrtim;                /* doubleing point number to record elapsed user CPU time */
    double systim;                /* doubleing point number to record elapsed system CPU time */

    int heights[4];
    parse_args(argc, argv, &final_state_file, &av_vels_file, &param_file);

    initialise(param_file, &accel_area, &params, &cells, &tmp_cells, &obstacles, &av_vels, &total_cells, heights);
    const float w1 = params.density * params.accel / 9.0;
    const float w2 = params.density * params.accel / 36.0;
    const unsigned total_num = params.nx * params.ny;
    const float omega_min_1 = (1 - params.omega);
    const float omega_magic = params.omega/omega_min_1;
    const float w0 = omega_magic*4.0/9.0;    /* weighting factor */
    const float w[] = {omega_magic/9.0, omega_magic/36.0};    /* weighting factor */
    const float one = 1.0;
    const float f_constants[] = {one, w0, w[0], w[1], omega_min_1, w1, w2 };
    const unsigned u_constants[] = {total_num, accel_area.col_or_row, accel_area.idx, heights[0], heights[1], params.nx, params.ny};
    float tot_u = 0.0;

    gettimeofday(&timstr,NULL);
    tic=timstr.tv_sec+(timstr.tv_usec/1000000.0);
    for (i = 0; i <= params.max_iters; i++)
    {
        // lbm_context
        tot_u = 0.0;
    #pragma omp parallel reduction(+:tot_u) proc_bind(master)
    {
        unsigned local_vars[7]; //x_e,x_w,y_n,y_s;
        unsigned ii;
        float tot_u_thread = 0.0;
        float local_density = 0.0;
        float u[5];
        float t[NSPEEDS];
        // unsigned ii, ri, rj;
        #pragma omp for schedule(auto)
        for (ii = *(u_constants + 3)**(u_constants + 5); ii < *(u_constants + 4)**(u_constants + 5); ii++) 
        {                
            *(local_vars + 1) = ii /(*(u_constants + 5));
            *(local_vars + 2) = ii%(*(u_constants + 5));
            if (*(local_vars + 1) == 0) {
                *(local_vars + 3) = *(u_constants + 6) - 1;
                *(local_vars + 4) = *(local_vars + 1) + 1;
            } else if (*(local_vars + 1) == *(u_constants + 6) - 1) {
                *(local_vars + 4) = 0;
                *(local_vars + 3) = *(local_vars + 1) - 1;
            } else {
                *(local_vars + 4) = *(local_vars + 1) + 1;
                *(local_vars + 3) = *(local_vars + 1) - 1;
            }
            if (*(local_vars + 2) == 0 ) {
                *(local_vars + 5) = *(u_constants + 5) - 1;
                *(local_vars + 6) = *(local_vars + 2) + 1;
            } else if ( *(local_vars + 2) == *(u_constants + 5) - 1) {
                *(local_vars + 6) = 0;
                *(local_vars + 5) = *(local_vars + 2) - 1;
            } else {
                *(local_vars + 6) = *(local_vars + 2) + 1;
                *(local_vars + 5) = *(local_vars + 2) - 1;
            }
            if (*(cells + ii) == -1) {
                *(tmp_cells + ii) = *(cells + ii);
                *(tmp_cells + 1**u_constants + ii) = *(cells + *u_constants*3 + *(local_vars + 1)**(u_constants + 5)  + *(local_vars + 6));
                *(tmp_cells + 2**u_constants + ii) = *(cells + *u_constants*4 + *(local_vars + 4)**(u_constants + 5) + *(local_vars + 2));
                *(tmp_cells + 3**u_constants + ii) = *(cells + *u_constants   + *(local_vars + 1)**(u_constants + 5)  + *(local_vars + 5));
                *(tmp_cells + 4**u_constants + ii) = *(cells + *u_constants*2 + *(local_vars + 3)**(u_constants + 5) + *(local_vars + 2));
                *(tmp_cells + 5**u_constants + ii) = *(cells + *u_constants*7 + *(local_vars + 4)**(u_constants + 5) + *(local_vars + 6));
                *(tmp_cells + 6**u_constants + ii) = *(cells + *u_constants*8 + *(local_vars + 4)**(u_constants + 5) + *(local_vars + 5));
                *(tmp_cells + 7**u_constants + ii) = *(cells + *u_constants*5 + *(local_vars + 3)**(u_constants + 5) + *(local_vars + 5));
                *(tmp_cells + 8**u_constants + ii) = *(cells + *u_constants*6 + *(local_vars + 3)**(u_constants + 5) + *(local_vars + 6));
            } else {
                *(t) = *(f_constants + 4)**(cells + ii);
                *(t + 1) = *(f_constants + 4)**(cells + *u_constants   + *(local_vars + 1)**(u_constants + 5)  + *(local_vars + 5));
                *(t + 2) = *(f_constants + 4)**(cells + *u_constants*2 + *(local_vars + 3)**(u_constants + 5) + *(local_vars + 2));
                *(t + 3) = *(f_constants + 4)**(cells + *u_constants*3 + *(local_vars + 1)**(u_constants + 5)  + *(local_vars + 6));
                *(t + 4) = *(f_constants + 4)**(cells + *u_constants*4 + *(local_vars + 4)**(u_constants + 5) + *(local_vars + 2));
                *(t + 5) = *(f_constants + 4)**(cells + *u_constants*5 + *(local_vars + 3)**(u_constants + 5) + *(local_vars + 5));
                *(t + 6) = *(f_constants + 4)**(cells + *u_constants*6 + *(local_vars + 3)**(u_constants + 5) + *(local_vars + 6));
                *(t + 7) = *(f_constants + 4)**(cells + *u_constants*7 + *(local_vars + 4)**(u_constants + 5) + *(local_vars + 6));
                *(t + 8) = *(f_constants + 4)**(cells + *u_constants*8 + *(local_vars + 4)**(u_constants + 5) + *(local_vars + 5));
                
                local_density = *(t) + *(t + 1) + *(t + 2) + *(t + 3) + *(t + 4) + *(t + 5) + *(t + 6) + *(t + 7) + *(t + 8);
                /* compute x velocity component */
                *(u + 1) = (*(t + 1) + *(t + 5) + *(t + 8) - (*(t + 3) + *(t + 6) + *(t + 7)))/local_density;
                /* compute y velocity component */
                *(u + 2) = (*(t + 2) + *(t + 5) + *(t + 6) - (*(t + 4) + *(t + 7) + *(t + 8)))/local_density;
                /* velocity squared */
                *u = *f_constants - 1.5*(*(u + 1) * *(u + 1) + *(u + 2) * *(u + 2));                                                                                                                        
                *(u + 3) = *(u + 1) + *(u + 2);
                *(u + 4) = -*(u + 1) + *(u + 2);                                                          
                *(tmp_cells + ii) =  *t + *(f_constants + 1) * local_density**u;
                *(tmp_cells + *u_constants + ii) = *(t + 1) + *(f_constants + 2) * local_density * (*u + (3.0**(u + 1))
                        + 4.5**(u + 1)**(u + 1));
                *(tmp_cells + 2**u_constants + ii) = *(t + 2)+ *(f_constants + 2) * local_density * (*u + (3.0**(u + 2))
                        + 4.5**(u + 2)**(u + 2));
                *(tmp_cells + 3**u_constants + ii) = *(t + 3)+ *(f_constants + 2) * local_density * (*u - (3.0**(u + 1))
                        + 4.5**(u + 1)**(u + 1));
                *(tmp_cells + 4**u_constants + ii) = *(t + 4)+ *(f_constants + 2) * local_density * (*u - (3.0**(u + 2))
                        + 4.5**(u + 2)**(u + 2));
                *(tmp_cells + 5**u_constants + ii) = *(t + 5)+ *(f_constants + 3) * local_density * (*u + (3.0**(u + 3))
                        + 4.5**(u + 3)**(u + 3));
                *(tmp_cells + 6**u_constants + ii) = *(t + 6)+ *(f_constants + 3) * local_density * (*u + (3.0**(u + 4))
                        + 4.5**(u + 4)**(u + 4));
                *(tmp_cells + 7**u_constants + ii) = *(t + 7)+ *(f_constants + 3) * local_density * (*u - (3.0**(u + 3))
                        + 4.5**(u + 3)**(u + 3));
                *(tmp_cells + 8**u_constants + ii) = *(t + 8)+ *(f_constants + 3) * local_density * (*u - (3.0**(u + 4))
                        + 4.5**(u + 4)**(u + 4));
                tot_u_thread = tot_u_thread + sqrt(*f_constants - *u);
                // tot_u = tot_u + 0.81*sqrt(one - *u);  
                if (!*(u_constants + 1) && *(local_vars + 2) == *(u_constants + 2)) {
                    /* if the cell is not occupied and
                    ** we don't send a density negative */
                    if (
                    (*(tmp_cells + 4**u_constants + ii) - *(f_constants + 5)) > 0.0 &&
                    (*(tmp_cells + 7**u_constants + ii) - *(f_constants + 6)) > 0.0 &&
                    (*(tmp_cells + 8**u_constants + ii) - *(f_constants + 6)) > 0.0 )
                    {
                         // increase 'north-side' densities 
                        *(tmp_cells + 2**u_constants + ii) += *(f_constants + 5);
                        *(tmp_cells + 5**u_constants + ii) += *(f_constants + 6);
                        *(tmp_cells + 6**u_constants + ii) += *(f_constants + 6);
                        /* decrease 'south-side' densities */
                        *(tmp_cells + 4**u_constants + ii) -= *(f_constants + 5);
                        *(tmp_cells + 7**u_constants + ii) -= *(f_constants + 6);
                        *(tmp_cells + 8**u_constants + ii) -= *(f_constants + 6);
                    }
                } else if ( *(u_constants + 1) && *(local_vars + 1) == *(u_constants + 2)) {
                    if (
                    (*(tmp_cells + 3**u_constants + ii) - *(f_constants + 5)) > 0.0 &&
                    (*(tmp_cells + 6**u_constants + ii) - *(f_constants + 6)) > 0.0 &&
                    (*(tmp_cells + 7**u_constants + ii) - *(f_constants + 6)) > 0.0 )
                    {
                        /* increase 'east-side' densities */
                        *(tmp_cells + 1**u_constants + ii) += *(f_constants + 5);
                        *(tmp_cells + 5**u_constants + ii) += *(f_constants + 6);
                        *(tmp_cells + 8**u_constants + ii) += *(f_constants + 6);
                        /* decrease 'west-side' densities */
                        *(tmp_cells + 3**u_constants + ii) -= *(f_constants + 5);
                        *(tmp_cells + 6**u_constants + ii) -= *(f_constants + 6);
                        *(tmp_cells + 7**u_constants + ii) -= *(f_constants + 6);
                    }
                }
            }
        }
            tot_u = tot_u + tot_u_thread;
        }

        *(av_vels + i) = tot_u;
        swap_cells(&cells, &tmp_cells);


    }

    gettimeofday(&timstr,NULL);
    toc=timstr.tv_sec+(timstr.tv_usec/1000000.0);
    getrusage(RUSAGE_SELF, &ru);
    timstr=ru.ru_utime;
    usrtim=timstr.tv_sec+(timstr.tv_usec/1000000.0);
    timstr=ru.ru_stime;
    systim=timstr.tv_sec+(timstr.tv_usec/1000000.0);

    printf("==done==\n");
    printf("Reynolds number:\t\t%.12E\n", calc_reynolds(params,cells,obstacles));
    printf("Elapsed time:\t\t\t%.6f (s)\n", toc-tic);
    printf("Elapsed user CPU time:\t\t%.6f (s)\n", usrtim);
    printf("Elapsed system CPU time:\t%.6f (s)\n", systim);

    write_values(final_state_file, av_vels_file, params, cells, obstacles, av_vels, total_cells);
    finalise(&cells, &tmp_cells, &obstacles, &av_vels);
    return EXIT_SUCCESS;
}

void write_values(const char * final_state_file, const char * av_vels_file,
    const param_t params, float* cells, unsigned char* obstacles, float* av_vels, unsigned total_cells)
{
    FILE* fp;                     /* file pointer */
    int ii,kk;                 /* generic counters */
    const float c_sq = 1.0/3.0;  /* sq. of speed of sound */
    const float my_constant = sqrt(2.0/3.0);
    float local_density;         /* per grid cell sum of densities */
    float pressure;              /* fluid pressure in grid cell */
    float u_x;                   /* x-component of velocity in grid cell */
    float u_y;                   /* y-component of velocity in grid cell */
    float u;                     /* norm--root of summed squares--of u_x and u_y */

    fp = fopen(final_state_file, "w");

    if (fp == NULL)
    {
        DIE("could not open file output file");
    }

    for (ii = 0; ii< params.ny * params.nx; ii++) {
        int obs_ii;
        if (cells[ii] == -1) obs_ii = 1;
        else obs_ii = 0;
        if (obstacles[ii]) {
            u_x = u_y = u = 0.0;
            pressure = params.density * c_sq;
        } else {
            local_density = 0.0;

            for (kk = 0; kk < NSPEEDS; kk++)
            {
                local_density += cells[kk*params.nx*params.ny+ii];
            }

            /* compute x velocity component */
            u_x = (cells[1*params.nx*params.ny+ii] +
                    cells[5*params.nx*params.ny+ii] +
                    cells[8*params.nx*params.ny+ii]
                - (cells[3*params.nx*params.ny+ii] +
                    cells[6*params.nx*params.ny+ii] +
                    cells[7*params.nx*params.ny+ii]))
                / local_density;

            /* compute y velocity component */
            u_y = (cells[2*params.nx*params.ny+ii] +
                    cells[5*params.nx*params.ny+ii] +
                    cells[6*params.nx*params.ny+ii]
                - (cells[4*params.nx*params.ny+ii] +
                    cells[7*params.nx*params.ny+ii] +
                    cells[8*params.nx*params.ny+ii]))
                / local_density;

            /* compute norm of velocity */
            u = sqrt((u_x * u_x) + (u_y * u_y));

            /* compute pressure */
            pressure = local_density * c_sq;
        }
        fprintf(fp,"%d %d %.12E %.12E %.12E %.12E %u\n",
            ii%params.nx,ii/params.nx,u_x,u_y,u,pressure, obstacles[ii]);
    }

    fclose(fp);

    fp = fopen(av_vels_file, "w");
    if (fp == NULL)
    {
        DIE("could not open file output file");
    }

    for (ii = 1; ii <= params.max_iters; ii++)
    {
        fprintf(fp,"%d:\t%.12E\n", ii-1, my_constant * av_vels[ii]/(float) total_cells);
    }

    fclose(fp);
}

float calc_reynolds(const param_t params, float* cells, unsigned char* obstacles)
{
    const float viscosity = 1.0 / 6.0 * (2.0 / params.omega - 1.0);

    return av_velocity(params,cells,obstacles) * params.reynolds_dim / viscosity;
}

float total_density(const param_t params, speed_t* cells)
{
    int ii,jj,kk;        /* generic counters */
    float total = 0.0;  /* accumulator */

    for (ii = 0; ii < params.ny; ii++)
    {
        for (jj = 0; jj < params.ny; jj++)
        {
            for (kk = 0; kk < NSPEEDS; kk++)
            {
                total += cells[ii*params.nx + jj].speeds[kk];
            }
        }
    }

    return total;
}



