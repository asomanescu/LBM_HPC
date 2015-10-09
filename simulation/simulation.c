/* Functions pertinent to the outer simulation steps */

#include <math.h>
#include <stdio.h>

#include "lbm.h"

void timestep(const param_t params, const accel_area_t accel_area,
    speed_t* cells, speed_t* tmp_cells, int* obstacles)
{
    accelerate_flow(params,accel_area,cells,obstacles);
    // propagate(params,cells,tmp_cells);
    // rebound(params,cells,tmp_cells,obstacles);
    collision(params,cells,tmp_cells,obstacles);
}

void accelerate_flow(const param_t params, const accel_area_t accel_area,
    speed_t* cells, int* obstacles)
{
    int ii,jj;     /* generic counters */
    float w1,w2;  /* weighting factors */

    /* compute weighting factors */
    w1 = params.density * params.accel / 9.0;
    w2 = params.density * params.accel / 36.0;

    if (accel_area.col_or_row == ACCEL_COLUMN)
    {
        jj = accel_area.idx;
        for (ii = 0; ii < params.ny; ii++)
        {
            /* if the cell is not occupied and
            ** we don't send a density negative */
            if (!obstacles[ii*params.nx + jj] &&
            (cells[ii*params.nx + jj].speeds[4] - w1) > 0.0 &&
            (cells[ii*params.nx + jj].speeds[7] - w2) > 0.0 &&
            (cells[ii*params.nx + jj].speeds[8] - w2) > 0.0 )
            {
                /* increase 'north-side' densities */
                cells[ii*params.nx + jj].speeds[2] += w1;
                cells[ii*params.nx + jj].speeds[5] += w2;
                cells[ii*params.nx + jj].speeds[6] += w2;
                /* decrease 'south-side' densities */
                cells[ii*params.nx + jj].speeds[4] -= w1;
                cells[ii*params.nx + jj].speeds[7] -= w2;
                cells[ii*params.nx + jj].speeds[8] -= w2;
            }
        }
    }
    else
    {
        // Know where to start from
        ii = accel_area.idx * params.nx;
        for (jj = ii; jj < (ii+params.nx); jj++)
        {
            /* if the cell is not occupied and
            ** we don't send a density negative */
            if (!obstacles[jj] &&
            (cells[jj].speeds[3] - w1) > 0.0 &&
            (cells[jj].speeds[6] - w2) > 0.0 &&
            (cells[jj].speeds[7] - w2) > 0.0 )
            {
                /* increase 'east-side' densities */
                cells[jj].speeds[1] += w1;
                cells[jj].speeds[5] += w2;
                cells[jj].speeds[8] += w2;
                /* decrease 'west-side' densities */
                cells[jj].speeds[3] -= w1;
                cells[jj].speeds[6] -= w2;
                cells[jj].speeds[7] -= w2;
            }
        }
    }
}

void propagate(const param_t params, speed_t* cells, speed_t* tmp_cells)
{
    int ii,jj;            /* generic counters */

    /* loop over _all_ cells */
    for (ii = 0; ii < params.ny; ii++)
    {
        for (jj = 0; jj < params.nx; jj++)
        {
            int x_e,x_w,y_n,y_s;  /* indices of neighbouring cells */
            /* determine indices of axis-direction neighbours
            ** respecting periodic boundary conditions (wrap around) */
            y_n = (ii + 1) % params.ny;
            x_e = (jj + 1) % params.nx;
            y_s = (ii == 0) ? (ii + params.ny - 1) : (ii - 1);
            x_w = (jj == 0) ? (jj + params.nx - 1) : (jj - 1);
            /* propagate densities to neighbouring cells, following
            ** appropriate directions of travel and writing into
            ** scratch space grid */
            tmp_cells[ii *params.nx + jj].speeds[0]  = cells[ii*params.nx + jj].speeds[0]; /* central cell, */
                                                     /* no movement   */
            tmp_cells[ii *params.nx + x_e].speeds[1] = cells[ii*params.nx + jj].speeds[1]; /* east */
            tmp_cells[y_n*params.nx + jj].speeds[2]  = cells[ii*params.nx + jj].speeds[2]; /* north */
            tmp_cells[ii *params.nx + x_w].speeds[3] = cells[ii*params.nx + jj].speeds[3]; /* west */
            tmp_cells[y_s*params.nx + jj].speeds[4]  = cells[ii*params.nx + jj].speeds[4]; /* south */
            tmp_cells[y_n*params.nx + x_e].speeds[5] = cells[ii*params.nx + jj].speeds[5]; /* north-east */
            tmp_cells[y_n*params.nx + x_w].speeds[6] = cells[ii*params.nx + jj].speeds[6]; /* north-west */
            tmp_cells[y_s*params.nx + x_w].speeds[7] = cells[ii*params.nx + jj].speeds[7]; /* south-west */
            tmp_cells[y_s*params.nx + x_e].speeds[8] = cells[ii*params.nx + jj].speeds[8]; /* south-east */
        }
    }
}

void rebound(const param_t params, speed_t* cells, speed_t* tmp_cells, int* obstacles)
{
    int ii,jj;  /* generic counters */

    /* loop over the cells in the grid */
    for (ii = 0; ii < params.ny; ii++)
    {
        for (jj = 0; jj < params.nx; jj++)
        {
            /* if the cell contains an obstacle */
            if (obstacles[ii*params.nx + jj])
            {
                /* called after propagate, so taking values from scratch space
                ** mirroring, and writing into main grid */
                cells[ii*params.nx + jj].speeds[1] = tmp_cells[ii*params.nx + jj].speeds[3];
                cells[ii*params.nx + jj].speeds[2] = tmp_cells[ii*params.nx + jj].speeds[4];
                cells[ii*params.nx + jj].speeds[3] = tmp_cells[ii*params.nx + jj].speeds[1];
                cells[ii*params.nx + jj].speeds[4] = tmp_cells[ii*params.nx + jj].speeds[2];
                cells[ii*params.nx + jj].speeds[5] = tmp_cells[ii*params.nx + jj].speeds[7];
                cells[ii*params.nx + jj].speeds[6] = tmp_cells[ii*params.nx + jj].speeds[8];
                cells[ii*params.nx + jj].speeds[7] = tmp_cells[ii*params.nx + jj].speeds[5];
                cells[ii*params.nx + jj].speeds[8] = tmp_cells[ii*params.nx + jj].speeds[6];
            }
        }
    }
}

void collision(const param_t params, speed_t* cells, speed_t* tmp_cells, int* obstacles)
{
    int ii,iii,jj,kk;                 /* generic counters */
    // const float c_sq = 1.0/3.0;  /* square of speed of sound */
    const float w0 = 4.0/9.0;    /* weighting factor */
    const float w[] = {1.0/9.0, 1.0/36.0};    /* weighting factor */
    const float one = 1.0;
    // const float w2 = 1.0/36.0;   /* weighting factor */

    float u_x,u_y;               /* av. velocities in x and y directions */
    float u_sq;                  /* squared velocity */
    float local_density;         /* sum of densities in a particular cell */
    float t[NSPEEDS];
    float u[NSPEEDS];            /* directional velocities */
    float d_equ;        /* equilibrium densities */

    int x_e,x_w,y_n,y_s;  /* indices of neighbouring cells */
            /* determine indices of axis-direction neighbours
            ** respecting periodic boundary conditions (wrap around) */
    /* loop over the cells in the grid
    ** NB the collision step is called after
    ** the propagate step and so values of interest
    ** are in the scratch-space grid */
    for (iii = 0; iii < params.ny; iii++)
    {
        for (jj = 0; jj< params.nx; jj++) {
            ii = iii*params.nx + jj;
            y_n = iii + 1;
            x_e = jj + 1;
            y_s = iii - 1;
            x_w = jj - 1;
            if (iii == 0) {
                y_s += params.ny;
            } else if (iii == params.ny - 1) {
                y_n = 0;
            }
            if (jj == 0 ) {
                x_w += params.nx;
            } else if ( jj == params.nx - 1) {
                x_e = 0;
            }
            // y_n = (iii + 1) % params.ny;
            // x_e = (jj + 1) % params.nx;
            // y_s = (iii == 0) ? (iii + params.ny - 1) : (iii - 1);
            // x_w = (jj == 0) ? (jj + params.nx - 1) : (jj - 1);
            t[0] = cells[ii].speeds[0];
            t[1] = cells[iii*params.nx + x_w].speeds[1];
            t[2] = cells[y_s*params.nx + jj].speeds[2];
            t[3] = cells[iii*params.nx + x_e].speeds[3];
            t[4] = cells[y_n*params.nx + jj].speeds[4];
            t[5] = cells[y_s*params.nx + x_w].speeds[5];
            t[6] = cells[y_s*params.nx + x_e].speeds[6];
            t[7] = cells[y_n*params.nx + x_e].speeds[7];
            t[8] = cells[y_n*params.nx + x_w].speeds[8];
            /* don't consider occupied cells */
            if (obstacles[ii]) {
                tmp_cells[ii].speeds[1] = t[3];
                tmp_cells[ii].speeds[2] = t[4];
                tmp_cells[ii].speeds[3] = t[1];
                tmp_cells[ii].speeds[4] = t[2];
                tmp_cells[ii].speeds[5] = t[7];
                tmp_cells[ii].speeds[6] = t[8];
                tmp_cells[ii].speeds[7] = t[5];
                tmp_cells[ii].speeds[8] = t[6];
            } else {
                /* compute local density total */
                local_density = 0.0;

                for (kk = 0; kk < NSPEEDS; kk++)
                {
                    local_density += t[kk];
                }

                /* compute x velocity component */
                u_x = (t[1] +
                        t[5] +
                        t[8]
                    - (t[3] +
                        t[6] +
                        t[7]))
                    / local_density;

                /* compute y velocity component */
                u_y = (t[2] +
                        t[5] +
                        t[6]
                    - (t[4] +
                        t[7] +
                        t[8]))
                    / local_density;

                /* velocity squared */
                u_sq = u_x * u_x + u_y * u_y;

                /* directional velocity components */
                u[1] =   u_x;        /* east */
                u[2] =         u_y;  /* north */
                u[3] = - u_x;        /* west */
                u[4] =       - u_y;  /* south */
                u[5] =   u_x + u_y;  /* north-east */
                u[6] = - u_x + u_y;  /* north-west */
                u[7] = - u_x - u_y;  /* south-west */
                u[8] =   u_x - u_y;  /* south-east */

                d_equ = w0 * local_density * (one - (3.0*u_sq) / 2.0);
                /* relaxation step */
                tmp_cells[ii].speeds[0] = 
                    (t[0] + params.omega * 
                        (d_equ - t[0]));
                for (kk = 1; kk < NSPEEDS; kk++)
                {
                    d_equ = w[(kk-1)/4] * local_density * (one + (3.0*u[kk])/ 1.0
                        +(9.0*u[kk]*u[kk]) / 2.0
                        - (3.0* u_sq )/ 2.0);
                    tmp_cells[ii].speeds[kk] =
                        (t[kk] + params.omega * 
                        (d_equ - t[kk]));
                }
            }
        }
    }
}

float av_velocity(const param_t params, speed_t* cells, int* obstacles)
{
    int    ii,kk;       /* generic counters */
    int    tot_cells = 0;  /* no. of cells used in calculation */
    float tot_u;          /* accumulated magnitudes of velocity for each cell */

    float local_density;  /* total density in cell */
    float u_x;            /* x-component of velocity for current cell */
    float u_y;            /* y-component of velocity for current cell */

    /* initialise */
    tot_u = 0.0;

    /* loop over all non-blocked cells */
    for (ii = 0; ii < params.ny * params.nx; ii++)
    {
        /* ignore occupied cells */
        if (!obstacles[ii])
        {
            /* local density total */
            local_density = 0.0;

            for (kk = 0; kk < NSPEEDS; kk++)
            {
                local_density += cells[ii].speeds[kk];
            }

            /* x-component of velocity */
            u_x = (cells[ii].speeds[1] +
                    cells[ii].speeds[5] +
                    cells[ii].speeds[8]
                - (cells[ii].speeds[3] +
                    cells[ii].speeds[6] +
                    cells[ii].speeds[7])) /
                local_density;

            /* compute y velocity component */
            u_y = (cells[ii].speeds[2] +
                    cells[ii].speeds[5] +
                    cells[ii].speeds[6]
                - (cells[ii].speeds[4] +
                    cells[ii].speeds[7] +
                    cells[ii].speeds[8])) /
                local_density;

            /* accumulate the norm of x- and y- velocity components */
            tot_u += sqrt(u_x*u_x + u_y*u_y);
            /* increase counter of inspected cells */
            ++tot_cells;
        }
    }

    return tot_u / (float)tot_cells;
}

