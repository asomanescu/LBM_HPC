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
int main(int argc, char* argv[])
{
    char * final_state_file = NULL;
    char * av_vels_file = NULL;
    char * param_file = NULL;

    accel_area_t accel_area;

    param_t  params;              /* struct to hold parameter values */
    float* cells     = NULL;    /* grid containing fluid densities */
    float* tmp_cells = NULL;    /* scratch space */
    int*     obstacles = NULL;    /* grid indicating which cells are blocked */
    float*  av_vels   = NULL;    /* a record of the av. velocity computed for each timestep */
    int total_cells;
    int    i;                    /*  generic counter */
    struct timeval timstr;        /* structure to hold elapsed time */
    struct rusage ru;             /* structure to hold CPU time--system and user */
    double tic,toc;               /* doubleing point numbers to calculate elapsed wallclock time */
    double usrtim;                /* doubleing point number to record elapsed user CPU time */
    double systim;                /* doubleing point number to record elapsed system CPU time */

    parse_args(argc, argv, &final_state_file, &av_vels_file, &param_file);

    initialise(param_file, &accel_area, &params, &cells, &tmp_cells, &obstacles, &av_vels, &total_cells);
    float* sp[2];
    sp[0] = cells;
    sp[1] = tmp_cells;
    const float w0 = 4.0/9.0;    /* weighting factor */
    const float w[] = {1.0/9.0, 1.0/36.0};    /* weighting factor */
    const int total_num = params.nx * params.ny;
    const float one = 1.0;
    float w1,w2;  /* weighting factors */
    // const float w2 = 1.0/36.0;   /* weighting factor */
           /* directional velocities */
    float tot_u = 0.0;
            /* determine indices of axis-direction neighbours
            ** respecting periodic boundary conditions (wrap around) */
    /* loop over the cells in the grid
    ** NB the collision step is called after
    ** the propagate step and so values of interest
    ** are in the scratch-space grid */
    int ii,kk,jj,ri,rj;                 /* generic counters */
    float u_sq;                  /* squared velocity */
    float local_density;         /* sum of densities in a particular cell */
    float t[NSPEEDS];
    float u[NSPEEDS]; 
    float d_equ;        /* equilibrium densities */
    /* iterate for max_iters timesteps */
    gettimeofday(&timstr,NULL);
    tic=timstr.tv_sec+(timstr.tv_usec/1000000.0);
// const float c_sq = 1.0/3.0;  /* square of speed of sound */

    for (i = 0; i < params.max_iters; i++)
    {

        /* compute weighting factors */
        tot_u = 0.0;
        w1 = params.density * params.accel / 9.0;
        w2 = params.density * params.accel / 36.0;

        if (accel_area.col_or_row == ACCEL_COLUMN)
        {
            // jj = accel_area.idx;
            // for (ii = 0; ii < params.ny; ii++)
            // {
            //     /* if the cell is not occupied and
            //     ** we don't send a density negative */
            //     if (!obstacles[ii*params.nx + jj] &&
            //     (cells[ii*params.nx + jj].speeds[4] - w1) > 0.0 &&
            //     (cells[ii*params.nx + jj].speeds[7] - w2) > 0.0 &&
            //     (cells[ii*params.nx + jj].speeds[8] - w2) > 0.0 )
            //     {
            //          // increase 'north-side' densities 
            //         cells[ii*params.nx + jj].speeds[2] += w1;
            //         cells[ii*params.nx + jj].speeds[5] += w2;
            //         cells[ii*params.nx + jj].speeds[6] += w2;
            //         /* decrease 'south-side' densities */
            //         cells[ii*params.nx + jj].speeds[4] -= w1;
            //         cells[ii*params.nx + jj].speeds[7] -= w2;
            //         cells[ii*params.nx + jj].speeds[8] -= w2;
            //     }
            // }
        }
        else
        {
            // Know where to start from
            ii = accel_area.idx * params.nx;

#pragma omp parallel for private(jj) schedule(guided)
            for (jj = ii; jj < (ii+params.nx); jj++) {
                /* if the cell is not occupied and
                ** we don't send a density negative */
                if (!obstacles[jj] &&
                (*(sp[i%2] + 3*total_num + jj) - w1) > 0.0 &&
                (*(sp[i%2] + 6*total_num + jj) - w2) > 0.0 &&
                (*(sp[i%2] + 7*total_num + jj) - w2) > 0.0 )
                {
                    /* increase 'east-side' densities */
                    *(sp[i%2] + 1*total_num + jj) += w1;
                    *(sp[i%2] + 5*total_num + jj) += w2;
                    *(sp[i%2] + 8*total_num + jj) += w2;
                    /* decrease 'west-side' densities */
                    *(sp[i%2] + 3*total_num + jj) -= w1;
                    *(sp[i%2] + 6*total_num + jj) -= w2;
                    *(sp[i%2] + 7*total_num + jj) -= w2;
                }
            }
        }

#pragma omp parallel private(ii, kk, ri, rj, u_sq, local_density, t , u, d_equ)
    {

#pragma omp for reduction(+:tot_u) schedule(guided)
    for (ii = 0; ii < total_num; ii++) {
        int x_e,x_w,y_n,y_s;  /* indices of neighbouring cells */
        // printf("%d .... %d \n", omp_get_thread_num(), ii);
        ri = ii/params.nx;
        rj = ii%params.nx;
        y_n = ri + 1;
        x_e = rj + 1;
        y_s = ri - 1;
        x_w = rj - 1;
        if (ri == 0) {
            y_s += params.ny;
        } else if (ri == params.ny - 1) {
            y_n = 0;
        }
        if (rj == 0 ) {
            x_w += params.nx;
        } else if ( rj == params.nx - 1) {
            x_e = 0;
        }
        *(t + 0) = *(sp[i%2] + ii);
        *(t + 1) = *(sp[i%2] + total_num   + ri*params.nx  + x_w);
        *(t + 2) = *(sp[i%2] + total_num*2 + y_s*params.nx + rj);
        *(t + 3) = *(sp[i%2] + total_num*3 + ri*params.nx  + x_e);
        *(t + 4) = *(sp[i%2] + total_num*4 + y_n*params.nx + rj);
        *(t + 5) = *(sp[i%2] + total_num*5 + y_s*params.nx + x_w);
        *(t + 6) = *(sp[i%2] + total_num*6 + y_s*params.nx + x_e);
        *(t + 7) = *(sp[i%2] + total_num*7 + y_n*params.nx + x_e);
        *(t + 8) = *(sp[i%2] + total_num*8 + y_n*params.nx + x_w);
        if (obstacles[ii]) {
            *(sp[(i+1)%2] + 1*total_num + ii) = *(t + 3);
            *(sp[(i+1)%2] + 2*total_num + ii) = *(t + 4);
            *(sp[(i+1)%2] + 3*total_num + ii) = *(t + 1);
            *(sp[(i+1)%2] + 4*total_num + ii) = *(t + 2);
            *(sp[(i+1)%2] + 5*total_num + ii) = *(t + 7);
            *(sp[(i+1)%2] + 6*total_num + ii) = *(t + 8);
            *(sp[(i+1)%2] + 7*total_num + ii) = *(t + 5);
            *(sp[(i+1)%2] + 8*total_num + ii) = *(t + 6);
        } else {
            
            local_density = 0.0;

            for (kk = 0; kk < NSPEEDS; kk++)
            {
                local_density += *(t + kk);
            }

            /* compute x velocity component */
            u[1] = (t[1] + t[5] + t[8] - (t[3] + t[6] + t[7])) / local_density;

            /* compute y velocity component */
            u[2] = (t[2] + t[5] + t[6] - (t[4] + t[7] + t[8])) / local_density;

            /* velocity squared */
            u_sq = u[1] * u[1] + u[2] * u[2];

            /* directional velocity components */
            // u[1] =   u_x;        /* east */
            // u[2] =         u_y;  /* north */
            u[3] = - u[1];        /* west */
            u[4] =        - u[2];  /* south */
            u[5] =   u[1] + u[2];  /* north-east */
            u[6] = - u[1] + u[2];  /* north-west */
            u[7] = - u[1] - u[2];  /* south-west */
            u[8] =   u[1] - u[2];  /* south-east */

            d_equ = w0 * local_density * (one - (3.0*u_sq) / 2.0);
            /* relaxation step */
            *(sp[(i+1)%2] + ii) = (t[0] + params.omega * (d_equ - t[0]));
            for (kk = 1; kk < NSPEEDS; kk++)
            {
                d_equ = w[(kk-1)/4] * local_density * (one + (3.0*u[kk])/ 1.0
                    +(9.0*u[kk]*u[kk]) / 2.0
                    - (3.0* u_sq )/ 2.0);
                *(sp[(i+1)%2] + kk*total_num + ii) = (t[kk] + params.omega * (d_equ - t[kk]));
            }
            tot_u = tot_u + sqrt(u_sq);  
        }

    }
}
    av_vels[i] = tot_u / (float)total_cells;
        // accelerate_flow(params,accel_area,sp[(i)%2],obstacles);
    // propagate(params,cells,tmp_cells);
    // rebound(params,cells,tmp_cells,obstacles);
        // av_vels[i] = collision(params,sp[(i)%2],sp[(i+1)%2],obstacles, total_cells);
        // av_vels[ii] = timestep(params, accel_area, sp[ii%2], sp[(ii+1)%2], obstacles, total_cells);
        // av_vels[ii] = av_velocity(params, sp[(ii+1)%2].spd, obstacles);

        #ifdef DEBUG
            printf("==timestep: %d==\n", ii);
            printf("av velocity: %.12E\n", av_vels[ii]);
            printf("tot density: %.12E\n", total_density(params, cells));
        #endif
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

    write_values(final_state_file, av_vels_file, params, cells, obstacles, av_vels);
    finalise(&cells, &tmp_cells, &obstacles, &av_vels);

    return EXIT_SUCCESS;
}

void write_values(const char * final_state_file, const char * av_vels_file,
    const param_t params, float* cells, int* obstacles, float* av_vels)
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
        fprintf(fp,"%d %d %.12E %.12E %.12E %.12E %d\n",
            ii%params.nx,ii/params.nx,u_x,u_y,u,pressure,obstacles[ii]);
    }
    // for (ii = 0; ii < params.ny; ii++)
    // {
    //     for (jj = 0; jj < params.nx; jj++)
    //     {
    //         /* an occupied cell */
    //         if (obstacles[ii*params.nx + jj])
    //         {
    //             u_x = u_y = u = 0.0;
    //             pressure = params.density * c_sq;
    //         }
    //         /* no obstacle */
    //         else
    //         {
    //             local_density = 0.0;

    //             for (kk = 0; kk < NSPEEDS; kk++)
    //             {
    //                 local_density += cells[ii*params.nx + jj].speeds[kk];
    //             }

    //             /* compute x velocity component */
    //             u_x = (cells[ii*params.nx + jj].speeds[1] +
    //                     cells[ii*params.nx + jj].speeds[5] +
    //                     cells[ii*params.nx + jj].speeds[8]
    //                 - (cells[ii*params.nx + jj].speeds[3] +
    //                     cells[ii*params.nx + jj].speeds[6] +
    //                     cells[ii*params.nx + jj].speeds[7]))
    //                 / local_density;

    //             /* compute y velocity component */
    //             u_y = (cells[ii*params.nx + jj].speeds[2] +
    //                     cells[ii*params.nx + jj].speeds[5] +
    //                     cells[ii*params.nx + jj].speeds[6]
    //                 - (cells[ii*params.nx + jj].speeds[4] +
    //                     cells[ii*params.nx + jj].speeds[7] +
    //                     cells[ii*params.nx + jj].speeds[8]))
    //                 / local_density;

    //             /* compute norm of velocity */
    //             u = sqrt((u_x * u_x) + (u_y * u_y));

    //             /* compute pressure */
    //             pressure = local_density * c_sq;
    //         }

    //         /* write to file */
    //         fprintf(fp,"%d %d %.12E %.12E %.12E %.12E %d\n",
    //             jj,ii,u_x,u_y,u,pressure,obstacles[ii*params.nx + jj]);
    //     }
    // }

    fclose(fp);

    fp = fopen(av_vels_file, "w");
    if (fp == NULL)
    {
        DIE("could not open file output file");
    }

    for (ii = 0; ii < params.max_iters; ii++)
    {
        fprintf(fp,"%d:\t%.12E\n", ii, av_vels[ii]);
    }

    fclose(fp);
}

float calc_reynolds(const param_t params, float* cells, int* obstacles)
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

