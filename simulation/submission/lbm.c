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
    const float w0 = 4.0/9.0;    /* weighting factor */
    const float w[] = {1.0/9.0, 1.0/36.0};    /* weighting factor */
    const float w1 = params.density * params.accel * w[0];
    const float w2 = params.density * params.accel * w[1];
    const unsigned total_num = params.nx * params.ny;
    const float one = 1.0;
    float tot_u = 0.0;

    gettimeofday(&timstr,NULL);
    tic=timstr.tv_sec+(timstr.tv_usec/1000000.0);
    for (i = 0; i <= params.max_iters; i++)
    {
        tot_u = 0.0;
    #pragma omp parallel reduction(+:tot_u) proc_bind(master)
    {
        unsigned x_e,x_w,y_n,y_s;
        float tot_u_thread = 0.0;
        float local_density = 0.0;
        float t[NSPEEDS];
        float u[5];
        unsigned ii, ri, rj;
        if(*(heights + 3) == 1) {
            #pragma omp for schedule(auto)
                for (ii = 0; ii< params.nx; ii++) {
                    ri = params.ny - 1;
                    rj = ii;
                    y_n = 0;
                    // y_s = ri - 1;
                    if (rj == 0 ) {
                        x_w = rj + params.nx - 1;
                        x_e = rj + 1;
                    } else if ( rj == params.nx - 1) {
                        x_e = 0;
                        x_w = rj - 1;
                    } else {
                        x_e = rj + 1;
                        x_w = rj - 1;
                    }
                    *(t + 2) = *(cells + total_num*2 + y_s*params.nx + rj);
                    *(t + 5) = *(cells + total_num*5 + y_s*params.nx + x_w);
                    *(t + 6) = *(cells + total_num*6 + y_s*params.nx + x_e);
                    *(tmp_cells + 4*total_num + ii) = *(t + 2);
                    *(tmp_cells + 7*total_num + ii) = *(t + 5);
                    *(tmp_cells + 8*total_num + ii) = *(t + 6);
                }
        }
        if (*(heights + 2) == 1)
        {
            #pragma omp for schedule(auto)
                for (ii = params.nx*(params.ny - 1); ii< total_num; ii++)
                {

                    ri = 0;
                    rj = ii;
                    y_s = ri + params.ny - 1;
                    if (rj == 0 ) {
                        x_w = rj + params.nx - 1;
                        x_e = rj + 1;
                    } else if ( rj == params.nx - 1) {
                        x_e = 0;
                        x_w = rj - 1;
                    } else {
                        x_e = rj + 1;
                        x_w = rj - 1;
                    }
                    *(t + 4) = *(cells + total_num*4 + y_n*params.nx + rj);
                    *(t + 7) = *(cells + total_num*7 + y_n*params.nx + x_e);
                    *(t + 8) = *(cells + total_num*8 + y_n*params.nx + x_w);
                    *(tmp_cells + 2*total_num + ii) = *(t + 4);
                    *(tmp_cells + 5*total_num + ii) = *(t + 7);
                    *(tmp_cells + 6*total_num + ii) = *(t + 8);
                }
        }
        #pragma omp for schedule(auto)
            for (ii = *heights*params.nx; ii < *(heights + 1)*params.nx; ii++) 
            {                
                ri = ii/params.nx;
                rj = ii%params.nx;
                if (ri == 0) {
                    y_s = ri + params.ny - 1;
                    y_n = ri + 1;
                } else if (ri == params.ny - 1) {
                    y_n = 0;
                    y_s = ri - 1;
                } else {
                    y_n = ri + 1;
                    y_s = ri - 1;
                }
                if (rj == 0 ) {
                    x_w = rj + params.nx - 1;
                    x_e = rj + 1;
                } else if ( rj == params.nx - 1) {
                    x_e = 0;
                    x_w = rj - 1;
                } else {
                    x_e = rj + 1;
                    x_w = rj - 1;
                }
                *t = *(cells + ii);
                *(t + 1) = *(cells + total_num   + ri*params.nx  + x_w);
                *(t + 2) = *(cells + total_num*2 + y_s*params.nx + rj);
                *(t + 3) = *(cells + total_num*3 + ri*params.nx  + x_e);
                *(t + 4) = *(cells + total_num*4 + y_n*params.nx + rj);
                *(t + 5) = *(cells + total_num*5 + y_s*params.nx + x_w);
                *(t + 6) = *(cells + total_num*6 + y_s*params.nx + x_e);
                *(t + 7) = *(cells + total_num*7 + y_n*params.nx + x_e);
                *(t + 8) = *(cells + total_num*8 + y_n*params.nx + x_w);
                if (obstacles[ii]) {
                    *(tmp_cells + 1*total_num + ii) = *(t + 3);
                    *(tmp_cells + 2*total_num + ii) = *(t + 4);
                    *(tmp_cells + 3*total_num + ii) = *(t + 1);
                    *(tmp_cells + 4*total_num + ii) = *(t + 2);
                    *(tmp_cells + 5*total_num + ii) = *(t + 7);
                    *(tmp_cells + 6*total_num + ii) = *(t + 8);
                    *(tmp_cells + 7*total_num + ii) = *(t + 5);
                    *(tmp_cells + 8*total_num + ii) = *(t + 6);
                } else {
                    
		    local_density = *t + *(t + 1) + *(t + 2) + *(t + 3) + *(t + 4) + *(t + 5) + *(t + 6) + *(t + 7) + *(t + 8);                    
                                                                                

                    /* compute x velocity component */
                    *(u + 1) = (*(t + 1) + *(t + 5) + *(t + 8) - (*(t + 3) + *(t + 6) + *(t + 7)))/local_density;

                    /* compute y velocity component */
                    *(u + 2) = (*(t + 2) + *(t + 5) + *(t + 6) - (*(t + 4) + *(t + 7) + *(t + 8)))/local_density;

                    /* velocity squared */
                    *u = *(u + 1) * *(u + 1) + *(u + 2) * *(u + 2);                                                                                                                        
                    *(u + 3) = *(u + 1) + *(u + 2);
                    *(u + 4) = -*(u + 1) + *(u + 2);
                    *(tmp_cells + ii) = (*t + params.omega * ((w0 * local_density - local_density*(2.0**u) / (3.0)) - *t));
                    *(tmp_cells + total_num + ii) = (*(t + 1) + params.omega * ((*w * local_density * (one + (3.0**(u + 1))
                            +(4.5**(u + 1)**(u + 1) - 1.5**u))) - *(t + 1)));
                    *(tmp_cells + 2*total_num + ii) = (*(t + 2) + params.omega * ((*w * local_density * (one + (3.0**(u + 2))
                            +(4.5**(u + 2)**(u + 2) - 1.5**u))) - *(t + 2)));
                    *(tmp_cells + 3*total_num + ii) = (*(t + 3) + params.omega * ((*w * local_density * (one - (3.0**(u + 1))
                            +(4.5**(u + 1)**(u + 1) - 1.5**u))) - *(t + 3)));
                    *(tmp_cells + 4*total_num + ii) = (*(t + 4) + params.omega * ((*w * local_density * (one - (3.0**(u + 2))
                            +(4.5**(u + 2)**(u + 2) - 1.5**u))) - *(t + 4)));
                    *(tmp_cells + 5*total_num + ii) = (*(t + 5) + params.omega * ((*(w + 1) * local_density * (one + (3.0**(u + 3))
                            +(4.5**(u + 3)**(u + 3) - 1.5**u))) - *(t + 5)));
                    *(tmp_cells + 6*total_num + ii) = (*(t + 6) + params.omega * ((*(w + 1) * local_density * (one + (3.0**(u + 4))
                            +(4.5**(u + 4)**(u + 4) - 1.5**u))) - *(t + 6)));
                    *(tmp_cells + 7*total_num + ii) = (*(t + 7) + params.omega * ((*(w + 1) * local_density * (one - (3.0**(u + 3))
                            +(4.5**(u + 3)**(u + 3) - 1.5**u))) - *(t + 7)));
                    *(tmp_cells + 8*total_num + ii) = (*(t + 8) + params.omega * ((*(w + 1) * local_density * (one - (3.0**(u + 4))
                            +(4.5**(u + 4)**(u + 4) - 1.5**u))) - *(t + 8)));
                    tot_u_thread = tot_u_thread + sqrt(*u);  
                    if (accel_area.col_or_row == ACCEL_COLUMN && rj == accel_area.idx) {
                        /* if the cell is not occupied and
                        ** we don't send a density negative */
                        if (
                        (*(tmp_cells + 4*total_num + ii) - w1) > 0.0 &&
                        (*(tmp_cells + 7*total_num + ii) - w2) > 0.0 &&
                        (*(tmp_cells + 8*total_num + ii) - w2) > 0.0 )
                        {
                             // increase 'north-side' densities 
                            *(tmp_cells + 2*total_num + ii) += w1;
                            *(tmp_cells + 5*total_num + ii) += w2;
                            *(tmp_cells + 6*total_num + ii) += w2;
                            /* decrease 'south-side' densities */
                            *(tmp_cells + 4*total_num + ii) -= w1;
                            *(tmp_cells + 7*total_num + ii) -= w2;
                            *(tmp_cells + 8*total_num + ii) -= w2;
                        }
                    } else if ( accel_area.col_or_row == ACCEL_ROW && ri == accel_area.idx) {
                        if (
                        (*(tmp_cells + 3*total_num + ii) - w1) > 0.0 &&
                        (*(tmp_cells + 6*total_num + ii) - w2) > 0.0 &&
                        (*(tmp_cells + 7*total_num + ii) - w2) > 0.0 )
                        {
                            /* increase 'east-side' densities */
                            *(tmp_cells + 1*total_num + ii) += w1;
                            *(tmp_cells + 5*total_num + ii) += w2;
                            *(tmp_cells + 8*total_num + ii) += w2;
                            /* decrease 'west-side' densities */
                            *(tmp_cells + 3*total_num + ii) -= w1;
                            *(tmp_cells + 6*total_num + ii) -= w2;
                            *(tmp_cells + 7*total_num + ii) -= w2;
                        }
                    }
                }
            }
        tot_u = tot_u + tot_u_thread;

    }

    *(av_vels +i) = tot_u;
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
            ii%params.nx,ii/params.nx,u_x,u_y,u,pressure,obstacles[ii]);
    }

    fclose(fp);

    fp = fopen(av_vels_file, "w");
    if (fp == NULL)
    {
        DIE("could not open file output file");
    }

    for (ii = 1; ii <= params.max_iters; ii++)
    {
        fprintf(fp,"%d:\t%.12E\n", ii-1, av_vels[ii]/(float) total_cells);
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

