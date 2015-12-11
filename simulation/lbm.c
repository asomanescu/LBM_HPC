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

#include "lbm.h"

#define MASTER 0


void swap_cells(speed_t **c, speed_t **t) {
  speed_t *tmp = *c;
  *c = *t;
  *t = tmp;
}

int main(int argc, char* argv[]) {
  int rank;               /* 'rank' of process among it's cohort */ 
  int size;               /* size of cohort, i.e. num processes started */
  int flag;               /* for checking whether MPI_Init() has been called */
  int strlen;             /* length of a character array */
  enum bool {FALSE,TRUE}; /* enumerated type: false = 0, true = 1 */  
  char hostname[MPI_MAX_PROCESSOR_NAME];  /* character array to hold hostname running process */

  // speed_t t;
  // t.speeds[0] = 0;
  // t.speeds[1] = 0;
  // t.speeds[2] = 0;
  // t.speeds[3] = 0;
  // t.speeds[4] = 0;
  // t.speeds[5] = 0;
  // t.speeds[6] = 0;
  // t.speeds[7] = 0;
  // t.speeds[8] = 0;

  MPI_Status status;
  MPI_Datatype MPI_SPEED, MPI_PARAMS, MPI_ACCEL_AREA;

  MPI_Datatype types[1] = {MPI_FLOAT}, params_types[7] = {MPI_INT, MPI_INT, MPI_INT, MPI_INT, MPI_FLOAT, MPI_FLOAT, MPI_FLOAT};
  MPI_Datatype accel_types[2] = {MPI_INT, MPI_INT};

  int blocklen[1] = {9}, params_len[7] = {1,1,1,1,1,1,1}, accel_len[2] = {1, 1}; 

  /* initialise our MPI environment */
  MPI_Init( &argc, &argv );

  /* check whether the initialisation was successful */
  MPI_Initialized(&flag);
  if ( flag != TRUE ) {
    MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
  }
  
  MPI_Aint disp[1] = {0}, params_disp[7],accel_disp[2];
  
  params_disp[0] = 0;
  params_disp[1] = sizeof(int);
  params_disp[2] = params_disp[1] + sizeof(int);
  params_disp[3] = params_disp[2] + sizeof(int);
  params_disp[4] = params_disp[3] + sizeof(int);
  params_disp[5] = params_disp[4] + sizeof(float);
  params_disp[6] = params_disp[5] + sizeof(float);

  accel_disp[0] = 0;
  accel_disp[1] = sizeof(int);


  // printf("%d", params_disp[rank]);
  MPI_Type_create_struct(1, blocklen, disp, types, &MPI_SPEED);
  MPI_Type_commit(&MPI_SPEED);
  MPI_Type_create_struct(7, params_len, params_disp, params_types, &MPI_PARAMS);
  MPI_Type_commit(&MPI_PARAMS);
  MPI_Type_create_struct(2, accel_len, accel_disp, accel_types, &MPI_ACCEL_AREA);
  MPI_Type_commit(&MPI_ACCEL_AREA);

  MPI_Comm_size( MPI_COMM_WORLD, &size );
  
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );

  char * final_state_file = NULL;
  char * av_vels_file = NULL;
  char * param_file = NULL;

  accel_area_t accel_area;

  param_t  params;              /* struct to hold parameter values */
  speed_t* cells     = NULL;    /* grid containing fluid densities */
  speed_t* local_cells = NULL;
  speed_t* tmp_cells = NULL;    /* scratch space */
  int*     obstacles = NULL;    /* grid indicating which cells are blocked */
  float*  av_vels   = NULL;    /* a record of the av. velocity computed for each timestep */
  unsigned total_cells;

  // int    ii;                    /*  generic counter */
  struct timeval timstr;        /* structure to hold elapsed time */
  struct rusage ru;             /* structure to hold CPU time--system and user */
  double tic,toc;               /* floating point numbers to calculate elapsed wallclock time */
  double usrtim;                /* floating point number to record elapsed user CPU time */
  double systim;                /* floating point number to record elapsed system CPU time */
  unsigned heights[2];

  // MPI_Type_struct(1, count, offsets, types, &myData);
  // MPI_Type_commit(&myData);


  if (rank == MASTER) {
    parse_args(argc, argv, &final_state_file, &av_vels_file, &param_file);
    initialise(param_file, &accel_area, &params, &cells, &tmp_cells, &obstacles, &av_vels, &total_cells, heights);
  }
  MPI_Bcast(&params, 1, MPI_PARAMS, MASTER, MPI_COMM_WORLD);
  MPI_Bcast(&accel_area, 1, MPI_ACCEL_AREA, MASTER, MPI_COMM_WORLD);
  MPI_Bcast(heights, 2, MPI_INT, MASTER, MPI_COMM_WORLD);
  unsigned bb_height = heights[1] - heights[0];
  unsigned gsz = params.nx*( bb_height/size );
  unsigned gsz_count[size];
  unsigned gsz_disp[size];
  int i;
  // if (rank == MASTER) {
  //   printf("%d %d %d %d\n", heights[0], heights[1], heights[1]-heights[0], params.nx*(heights[1]-heights[0])/size);
  // }
  for( i = 0; i < size; i++) {
    gsz_count[i] = gsz;
    if ( i < bb_height % size ) {
      gsz_count[i] += params.nx;
    }
    gsz_disp[i] = (i != 0) ? (gsz_disp[i - 1] + gsz_count[i - 1]) : heights[0]*params.nx;
  }
  unsigned local_gsz = gsz_count[rank];
  // printf("%d \n", local_gsz);
  local_cells = (speed_t*) malloc(sizeof(speed_t)*local_gsz);
  // free(cells);
  speed_t* local_tmp_cells = (speed_t *) malloc(sizeof(speed_t)*local_gsz);
  MPI_Scatterv(cells, gsz_count, gsz_disp, MPI_SPEED, local_cells, local_gsz, MPI_SPEED, MASTER, MPI_COMM_WORLD);
  // MPI_Scatterv(cells, gsz_count, gsz_disp, MPI_SPEED, local_tmp_cells, local_gsz, MPI_SPEED, MASTER, MPI_COMM_WORLD);
  speed_t* hallo_1 = (speed_t *)malloc(sizeof(speed_t)*params.nx);
  speed_t* hallo_2 = (speed_t *)malloc(sizeof(speed_t)*params.nx);
  float *local_av_vels = (float *)malloc(sizeof(float)*params.max_iters);
  float tot_u;
  unsigned x_e,x_w,y_n,y_s;
  float local_density = 0.0;
  float t[NSPEEDS];
  float u[5];
  const unsigned rank_b = ( rank != 0 ) ? (rank - 1) : (size - 1);
  const unsigned rank_t = ( rank != size - 1) ? ( rank + 1 ) : 0;
  const float w0 = 4.0/9.0;    /* weighting factor */
  const float w[] = {1.0/9.0, 1.0/36.0};    /* weighting factor */
  const float w1 = params.density * params.accel * w[0];
  const float w2 = params.density * params.accel * w[1];
  unsigned iii, ii, ri, rj;
  const unsigned nx = params.nx;
  const float omega = params.omega;
  const float minus_omega = 1.0f - omega;
  accel_area.idx = (accel_area.col_or_row == 0) ? accel_area.idx * nx - gsz_disp[rank] : accel_area.idx;




//   // 0.044444 0.011111 0.002778
// //      iterate for max_iters timesteps 
  MPI_Request recv_up_req, send_up_req, recv_down_req, send_down_req;
  MPI_Status recv_up_status, recv_down_status, send_up_status, send_down_status;
  if ( rank == MASTER ) {
    gettimeofday(&timstr,NULL);
    tic=timstr.tv_sec+(timstr.tv_usec/1000000.0);
  }


  for (iii = 0; iii < params.max_iters; iii++)
  {
    tot_u = 0.0f;
    MPI_Isend(local_cells, nx, MPI_SPEED, (size + rank - 1)%size, iii, MPI_COMM_WORLD, &send_down_req);
    MPI_Isend(&local_cells[local_gsz - nx], nx, MPI_SPEED, (size + rank + 1)%size, iii, MPI_COMM_WORLD, &send_up_req);
    MPI_Irecv(hallo_2, nx, MPI_SPEED, (size + rank - 1)%size, iii, MPI_COMM_WORLD, &recv_down_req);
    MPI_Irecv(hallo_1, nx, MPI_SPEED, (size + rank + 1)%size, iii, MPI_COMM_WORLD, &recv_up_req);
    // if ( rank % 2 == 0 ) {
    //   // MPI_Ssend(local_cells, params.nx, MPI_SPEED, (size + rank - 1)%size, iii, MPI_COMM_WORLD);
    //   // MPI_Recv(hallo_2, params.nx, MPI_SPEED, (size + rank - 1)%size, iii, MPI_COMM_WORLD, &status);
    //   // MPI_Ssend(&local_cells[local_gsz - params.nx], params.nx, MPI_SPEED, (size + rank + 1)%size, iii, MPI_COMM_WORLD);
    //   // MPI_Recv(hallo_1, params.nx, MPI_SPEED, (size + rank + 1)%size, iii, MPI_COMM_WORLD, &status);
    // } else {
    //   MPI_Irecv(hallo_1, nx, MPI_SPEED, (size + rank + 1)%size, iii, MPI_COMM_WORLD, &req);
    //   MPI_Wait(&req, &status);
    //   // MPI_Recv(hallo_1, params.nx, MPI_SPEED, (size + rank + 1)%size, iii, MPI_COMM_WORLD, &status);
    //   // MPI_Ssend(&local_cells[local_gsz - params.nx], params.nx, MPI_SPEED, (size + rank + 1)%size, iii, MPI_COMM_WORLD);
    //   // MPI_Recv(hallo_2, params.nx, MPI_SPEED, (size + rank - 1)%size, iii, MPI_COMM_WORLD, &status);
    //   // MPI_Ssend(local_cells, params.nx, MPI_SPEED, (size + rank - 1)%size, iii, MPI_COMM_WORLD);
    // }
  //   // printf("Iter %d on process %d\n", iii, rank);
    for(ii = nx; ii < local_gsz - nx; ii++) {
      ri = ii/nx;
      ri = ri * nx;
      // rj = ii%nx;
      rj = ii - ri;
      y_s = (ri - nx);
      y_n = ( ri + nx );
      x_w = (rj != 0) ? (rj - 1) : (nx - 1);
      x_e = (rj != nx - 1) ? (rj + 1) : 0;

      float t[NSPEEDS]; //local_cells[ii];

      *t = local_cells[ii].speeds[0];
      *(t + 1) = local_cells[ri + x_w].speeds[1];
      *(t + 2) = local_cells[y_s + rj].speeds[2];
      *(t + 3) = local_cells[ri + x_e].speeds[3];
      *(t + 4) = local_cells[y_n + rj].speeds[4];
      *(t + 5) = local_cells[y_s + x_w].speeds[5];
      *(t + 6) = local_cells[y_s + x_e].speeds[6];
      *(t + 7) = local_cells[y_n + x_e].speeds[7];
      *(t + 8) = local_cells[y_n + x_w].speeds[8];

      // if ( iii == 0 && rj == 0) {
      //   printf("%d %d %d %f \n", ri, gsz_disp[rank] + ri, accel_area.idx*nx, t[0]);
      // }
      if(t[0] != -1) {
        local_density = t[0] + t[1] + t[2] + t[3] + t[4] + t[5] + t[6] + t[7] + t[8];

        u[1] = (t[1] + t[5] + t[8] - (t[3] + t[6] + t[7])) / local_density;

        u[2] = (t[2] + t[5] + t[6] - (t[4] + t[7] + t[8])) / local_density;

        u[0] = u[1]*u[1] + u[2]*u[2];

        *(u + 3) = *(u + 1) + *(u + 2);
        *(u + 4) = -*(u + 1) + *(u + 2);

        tot_u += sqrt(*u);

        *u = 1.0f - 1.5f**u;
        local_density = 0.4444444f*local_density*omega;
        // *u = local_density**u;
        t[0] = t[0] * minus_omega + local_density**u; 
        local_density *= 0.25f;
        t[1] = t[1] * minus_omega + local_density * (u[0] + 4.5f * u[1] * u[1] + 3.0f * u[1]);
        t[2] = t[2] * minus_omega + local_density * (u[0] + 4.5f * u[2] * u[2] + 3.0f * u[2]);
        t[3] = t[3] * minus_omega + local_density * (u[0] + 4.5f * u[1] * u[1] - 3.0f * u[1]);
        t[4] = t[4] * minus_omega + local_density * (u[0] + 4.5f * u[2] * u[2] - 3.0f * u[2]);
        local_density *= 0.25f;
        t[5] = t[5] * minus_omega + local_density * (u[0] + 4.5f * u[3] * u[3] + 3.0f * u[3]);
        t[6] = t[6] * minus_omega + local_density * (u[0] + 4.5f * u[4] * u[4] + 3.0f * u[4]);
        t[7] = t[7] * minus_omega + local_density * (u[0] + 4.5f * u[3] * u[3] - 3.0f * u[3]);
        t[8] = t[8] * minus_omega + local_density * (u[0] + 4.5f * u[4] * u[4] - 3.0f * u[4]);
        if (accel_area.col_or_row == 1 && rj == accel_area.idx) {
          if (
          (*(t + 4) - w1) > 0.0 &&
          (*(t + 7) - w2) > 0.0 &&
          (*(t + 8) - w2) > 0.0 )
          {
             // increase 'north-side' densities 
            *(t + 2) += w1;
            *(t + 5) += w2;
            *(t + 6) += w2;
            /* decrease 'south-side' densities */
            *(t + 4) -= w1;
            *(t + 7) -= w2;
            *(t + 8) -= w2;
          }
        } else if ( accel_area.col_or_row == 0 && (ri == accel_area.idx)) {
          if (
          (*(t + 3) - w1) > 0.0 &&
          (*(t + 6) - w2) > 0.0 &&
          (*(t + 7) - w2) > 0.0 )
          {
            /* increase 'east-side' densities */
            *(t + 1) += w1;
            *(t + 5) += w2;
            *(t + 8) += w2;
            /* decrease 'west-side' densities */
            *(t + 3) -= w1;
            *(t + 6) -= w2;
            *(t + 7) -= w2;
          }
        }
      }
      // speed_t t_to_send;
      local_tmp_cells[ii].speeds[0] = t[0];
      local_tmp_cells[ii].speeds[1] = (t[0] == -1) ? t[3] : t[1];
      local_tmp_cells[ii].speeds[2] = (t[0] == -1) ? t[4] : t[2];
      local_tmp_cells[ii].speeds[3] = (t[0] == -1) ? t[1] : t[3];
      local_tmp_cells[ii].speeds[4] = (t[0] == -1) ? t[2] : t[4];
      local_tmp_cells[ii].speeds[5] = (t[0] == -1) ? t[7] : t[5];
      local_tmp_cells[ii].speeds[6] = (t[0] == -1) ? t[8] : t[6];
      local_tmp_cells[ii].speeds[7] = (t[0] == -1) ? t[5] : t[7];
      local_tmp_cells[ii].speeds[8] = (t[0] == -1) ? t[6] : t[8];
      // local_tmp_cells[ii] = t_to_send;
    }

    MPI_Wait(&send_up_req, &send_up_status);
    MPI_Wait(&recv_up_req, &recv_up_status);
    for(ii = local_gsz - nx; ii < local_gsz; ii++) {
      ri = local_gsz - nx;
      // ri = ri * nx;
      // rj = ii%nx;
      rj = ii - ri;
      y_s = (ri - nx);
      // y_n = (ri != local_gsz - nx) ? ( ri + nx ) : -1;
      x_w = (rj != 0) ? (rj - 1) : (nx - 1);
      x_e = (rj != nx - 1) ? (rj + 1) : 0;

      float t[NSPEEDS]; //local_cells[ii];

      *t = local_cells[ii].speeds[0];
      *(t + 1) = local_cells[ri + x_w].speeds[1];
      *(t + 2) = local_cells[y_s + rj].speeds[2];
      *(t + 3) = local_cells[ri + x_e].speeds[3];
      *(t + 4) = hallo_1[rj].speeds[4];
      *(t + 5) = local_cells[y_s + x_w].speeds[5];
      *(t + 6) = local_cells[y_s + x_e].speeds[6];
      *(t + 7) = hallo_1[x_e].speeds[7];
      *(t + 8) = hallo_1[x_w].speeds[8];

      // if ( iii == 0 && rj == 0) {
      //   printf("%d %d %d %f \n", ri, gsz_disp[rank] + ri, accel_area.idx*nx, t[0]);
      // }
      if(t[0] != -1) {
        local_density = t[0] + t[1] + t[2] + t[3] + t[4] + t[5] + t[6] + t[7] + t[8];

        u[1] = (t[1] + t[5] + t[8] - (t[3] + t[6] + t[7])) / local_density;

        u[2] = (t[2] + t[5] + t[6] - (t[4] + t[7] + t[8])) / local_density;

        u[0] = u[1]*u[1] + u[2]*u[2];

        *(u + 3) = *(u + 1) + *(u + 2);
        *(u + 4) = -*(u + 1) + *(u + 2);

        tot_u += sqrt(*u);

        *u = 1.0f - 1.5f**u;
        local_density = 0.4444444f*local_density*omega;
        // *u = local_density**u;
        t[0] = t[0] * minus_omega + local_density**u; 
        local_density *= 0.25f;
        t[1] = t[1] * minus_omega + local_density * (u[0] + 4.5f * u[1] * u[1] + 3.0f * u[1]);
        t[2] = t[2] * minus_omega + local_density * (u[0] + 4.5f * u[2] * u[2] + 3.0f * u[2]);
        t[3] = t[3] * minus_omega + local_density * (u[0] + 4.5f * u[1] * u[1] - 3.0f * u[1]);
        t[4] = t[4] * minus_omega + local_density * (u[0] + 4.5f * u[2] * u[2] - 3.0f * u[2]);
        local_density *= 0.25f;
        t[5] = t[5] * minus_omega + local_density * (u[0] + 4.5f * u[3] * u[3] + 3.0f * u[3]);
        t[6] = t[6] * minus_omega + local_density * (u[0] + 4.5f * u[4] * u[4] + 3.0f * u[4]);
        t[7] = t[7] * minus_omega + local_density * (u[0] + 4.5f * u[3] * u[3] - 3.0f * u[3]);
        t[8] = t[8] * minus_omega + local_density * (u[0] + 4.5f * u[4] * u[4] - 3.0f * u[4]);
        if (accel_area.col_or_row == 1 && rj == accel_area.idx) {
          if (
          (*(t + 4) - w1) > 0.0 &&
          (*(t + 7) - w2) > 0.0 &&
          (*(t + 8) - w2) > 0.0 )
          {
             // increase 'north-side' densities 
            *(t + 2) += w1;
            *(t + 5) += w2;
            *(t + 6) += w2;
            /* decrease 'south-side' densities */
            *(t + 4) -= w1;
            *(t + 7) -= w2;
            *(t + 8) -= w2;
          }
        } else if ( accel_area.col_or_row == 0 && (ri == accel_area.idx)) {
          if (
          (*(t + 3) - w1) > 0.0 &&
          (*(t + 6) - w2) > 0.0 &&
          (*(t + 7) - w2) > 0.0 )
          {
            /* increase 'east-side' densities */
            *(t + 1) += w1;
            *(t + 5) += w2;
            *(t + 8) += w2;
            /* decrease 'west-side' densities */
            *(t + 3) -= w1;
            *(t + 6) -= w2;
            *(t + 7) -= w2;
          }
        }
      }
      // speed_t t_to_send;
      local_tmp_cells[ii].speeds[0] = t[0];
      local_tmp_cells[ii].speeds[1] = (t[0] == -1) ? t[3] : t[1];
      local_tmp_cells[ii].speeds[2] = (t[0] == -1) ? t[4] : t[2];
      local_tmp_cells[ii].speeds[3] = (t[0] == -1) ? t[1] : t[3];
      local_tmp_cells[ii].speeds[4] = (t[0] == -1) ? t[2] : t[4];
      local_tmp_cells[ii].speeds[5] = (t[0] == -1) ? t[7] : t[5];
      local_tmp_cells[ii].speeds[6] = (t[0] == -1) ? t[8] : t[6];
      local_tmp_cells[ii].speeds[7] = (t[0] == -1) ? t[5] : t[7];
      local_tmp_cells[ii].speeds[8] = (t[0] == -1) ? t[6] : t[8];
      // local_tmp_cells[ii] = t_to_send;
    }
    // printf("%d finished?", MPI_Test(&req, &iii, &status));
    // MPI_Wait(&req, &status);
    // printf("wait finished\n");
    MPI_Wait(&send_down_req, &send_down_status);
    MPI_Wait(&recv_down_req, &recv_down_status);
    for(ii = 0; ii< nx; ii++) {
      ri = 0;
      rj = ii - ri;
      // y_s = (ri != 0) ? (ri - nx) : -1;
      y_n = ( ri + nx );
      x_w = (rj != 0) ? (rj - 1) : (nx - 1);
      x_e = (rj != nx - 1) ? (rj + 1) : 0;

      float t[NSPEEDS]; //local_cells[ii];

      *t = local_cells[ii].speeds[0];
      *(t + 1) = local_cells[ri + x_w].speeds[1];
      *(t + 2) = hallo_2[rj].speeds[2];
      *(t + 3) = local_cells[ri + x_e].speeds[3];
      *(t + 4) = local_cells[y_n + rj].speeds[4];
      *(t + 5) = hallo_2[x_w].speeds[5];
      *(t + 6) = hallo_2[x_e].speeds[6];
      *(t + 7) = local_cells[y_n + x_e].speeds[7];
      *(t + 8) = local_cells[y_n + x_w].speeds[8];

      // if ( iii == 0 && rj == 0) {
      //   printf("%d %d %d %f \n", ri, gsz_disp[rank] + ri, accel_area.idx*nx, t[0]);
      // }
      if(t[0] != -1) {
        local_density = t[0] + t[1] + t[2] + t[3] + t[4] + t[5] + t[6] + t[7] + t[8];

        u[1] = (t[1] + t[5] + t[8] - (t[3] + t[6] + t[7])) / local_density;

        u[2] = (t[2] + t[5] + t[6] - (t[4] + t[7] + t[8])) / local_density;

        u[0] = u[1]*u[1] + u[2]*u[2];

        *(u + 3) = *(u + 1) + *(u + 2);
        *(u + 4) = -*(u + 1) + *(u + 2);

        tot_u += sqrt(*u);

        *u = 1.0f - 1.5f**u;
        local_density = 0.4444444f*local_density*omega;
        // *u = local_density**u;
        t[0] = t[0] * minus_omega + local_density**u; 
        local_density *= 0.25f;
        t[1] = t[1] * minus_omega + local_density * (u[0] + 4.5f * u[1] * u[1] + 3.0f * u[1]);
        t[2] = t[2] * minus_omega + local_density * (u[0] + 4.5f * u[2] * u[2] + 3.0f * u[2]);
        t[3] = t[3] * minus_omega + local_density * (u[0] + 4.5f * u[1] * u[1] - 3.0f * u[1]);
        t[4] = t[4] * minus_omega + local_density * (u[0] + 4.5f * u[2] * u[2] - 3.0f * u[2]);
        local_density *= 0.25f;
        t[5] = t[5] * minus_omega + local_density * (u[0] + 4.5f * u[3] * u[3] + 3.0f * u[3]);
        t[6] = t[6] * minus_omega + local_density * (u[0] + 4.5f * u[4] * u[4] + 3.0f * u[4]);
        t[7] = t[7] * minus_omega + local_density * (u[0] + 4.5f * u[3] * u[3] - 3.0f * u[3]);
        t[8] = t[8] * minus_omega + local_density * (u[0] + 4.5f * u[4] * u[4] - 3.0f * u[4]);
        if (accel_area.col_or_row == 1 && rj == accel_area.idx) {
          if (
          (*(t + 4) - w1) > 0.0 &&
          (*(t + 7) - w2) > 0.0 &&
          (*(t + 8) - w2) > 0.0 )
          {
             // increase 'north-side' densities 
            *(t + 2) += w1;
            *(t + 5) += w2;
            *(t + 6) += w2;
            /* decrease 'south-side' densities */
            *(t + 4) -= w1;
            *(t + 7) -= w2;
            *(t + 8) -= w2;
          }
        } else if ( accel_area.col_or_row == 0 && (ri == accel_area.idx)) {
          if (
          (*(t + 3) - w1) > 0.0 &&
          (*(t + 6) - w2) > 0.0 &&
          (*(t + 7) - w2) > 0.0 )
          {
            /* increase 'east-side' densities */
            *(t + 1) += w1;
            *(t + 5) += w2;
            *(t + 8) += w2;
            /* decrease 'west-side' densities */
            *(t + 3) -= w1;
            *(t + 6) -= w2;
            *(t + 7) -= w2;
          }
        }
      }
      // speed_t t_to_send;
      local_tmp_cells[ii].speeds[0] = t[0];
      local_tmp_cells[ii].speeds[1] = (t[0] == -1) ? t[3] : t[1];
      local_tmp_cells[ii].speeds[2] = (t[0] == -1) ? t[4] : t[2];
      local_tmp_cells[ii].speeds[3] = (t[0] == -1) ? t[1] : t[3];
      local_tmp_cells[ii].speeds[4] = (t[0] == -1) ? t[2] : t[4];
      local_tmp_cells[ii].speeds[5] = (t[0] == -1) ? t[7] : t[5];
      local_tmp_cells[ii].speeds[6] = (t[0] == -1) ? t[8] : t[6];
      local_tmp_cells[ii].speeds[7] = (t[0] == -1) ? t[5] : t[7];
      local_tmp_cells[ii].speeds[8] = (t[0] == -1) ? t[6] : t[8];
      // local_tmp_cells[ii] = t_to_send;
    }
  //     // timestep(params, accel_area, cells, tmp_cells, obstacles);
    local_av_vels[iii] = tot_u; //av_velocity(params, cells, obstacles);
    swap_cells(&local_cells, &local_tmp_cells);
  }

  if ( rank == MASTER ) {
    gettimeofday(&timstr,NULL);
    toc=timstr.tv_sec+(timstr.tv_usec/1000000.0);
    getrusage(RUSAGE_SELF, &ru);
    timstr=ru.ru_utime;
    usrtim=timstr.tv_sec+(timstr.tv_usec/1000000.0);
    timstr=ru.ru_stime;
    systim=timstr.tv_sec+(timstr.tv_usec/1000000.0);
  }
  MPI_Reduce(local_av_vels, av_vels, params.max_iters, MPI_FLOAT, MPI_SUM, MASTER, MPI_COMM_WORLD);
  if( rank == MASTER) {
    printf("==done==\n");
    printf("Reynolds number:\t\t%.12E\n", calc_reynolds(params, av_vels[params.max_iters - 1], total_cells));
    printf("Elapsed time:\t\t\t%.6f (s)\n", toc-tic);
    printf("Elapsed user CPU time:\t\t%.6f (s)\n", usrtim);
    printf("Elapsed system CPU time:\t%.6f (s)\n", systim);
  }
  // int i;
  MPI_Gatherv(local_cells, local_gsz, MPI_SPEED, cells, gsz_count, gsz_disp, MPI_SPEED, MASTER, MPI_COMM_WORLD);
  if (rank == MASTER) {
    write_values(final_state_file, av_vels_file, params, cells, obstacles, av_vels, total_cells);
    finalise(&cells, &tmp_cells, &obstacles, &av_vels);    
  }
  // MPI_Scatter(checks, 1, MPI_INT, &check, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
  // free(checks);
  // printf("%d id ==== changed %d\n", rank, check);

  /* determine the hostname */
  // MPI_Get_processor_name(hostname,&strlen);


  // printf("Hello, world; from host %s: process %d of %d\n", hostname, rank, size);

  /* finialise the MPI enviroment */
  MPI_Finalize();

  /* and exit the program */
  return EXIT_SUCCESS;
}

/*
** main program:
** initialise, timestep loop, finalise
*/
// int main(int argc, char* argv[])
// {
//     char * final_state_file = NULL;
//     char * av_vels_file = NULL;
//     char * param_file = NULL;

//     accel_area_t accel_area;

//     param_t  params;              /* struct to hold parameter values */
//     speed_t* cells     = NULL;    /* grid containing fluid densities */
//     speed_t* tmp_cells = NULL;    /* scratch space */
//     int*     obstacles = NULL;    /* grid indicating which cells are blocked */
//     float*  av_vels   = NULL;    /* a record of the av. velocity computed for each timestep */

//     int    ii;                    /*  generic counter */
//     struct timeval timstr;        /* structure to hold elapsed time */
//     struct rusage ru;             /* structure to hold CPU time--system and user */
//     double tic,toc;               /* floating point numbers to calculate elapsed wallclock time */
//     double usrtim;                /* floating point number to record elapsed user CPU time */
//     double systim;                /* floating point number to record elapsed system CPU time */

//     parse_args(argc, argv, &final_state_file, &av_vels_file, &param_file);

//     initialise(param_file, &accel_area, &params, &cells, &tmp_cells, &obstacles, &av_vels);

//      iterate for max_iters timesteps 
//     gettimeofday(&timstr,NULL);
//     tic=timstr.tv_sec+(timstr.tv_usec/1000000.0);

//     for (ii = 0; ii < params.max_iters; ii++)
//     {
//         timestep(params, accel_area, cells, tmp_cells, obstacles);
//         av_vels[ii] = av_velocity(params, cells, obstacles);

//         #ifdef DEBUG
//             printf("==timestep: %d==\n", ii);
//             printf("av velocity: %.12E\n", av_vels[ii]);
//             printf("tot density: %.12E\n", total_density(params, cells));
//         #endif
//     }

//     gettimeofday(&timstr,NULL);
//     toc=timstr.tv_sec+(timstr.tv_usec/1000000.0);
//     getrusage(RUSAGE_SELF, &ru);
//     timstr=ru.ru_utime;
//     usrtim=timstr.tv_sec+(timstr.tv_usec/1000000.0);
//     timstr=ru.ru_stime;
//     systim=timstr.tv_sec+(timstr.tv_usec/1000000.0);

//     printf("==done==\n");
//     printf("Reynolds number:\t\t%.12E\n", calc_reynolds(params,cells,obstacles));
//     printf("Elapsed time:\t\t\t%.6f (s)\n", toc-tic);
//     printf("Elapsed user CPU time:\t\t%.6f (s)\n", usrtim);
//     printf("Elapsed system CPU time:\t%.6f (s)\n", systim);

//     write_values(final_state_file, av_vels_file, params, cells, obstacles, av_vels);
//     finalise(&cells, &tmp_cells, &obstacles, &av_vels);

//     return EXIT_SUCCESS;
// }

void write_values(const char * final_state_file, const char * av_vels_file,
    const param_t params, speed_t* cells, int* obstacles, float* av_vels, unsigned total_cells)
{
    FILE* fp;                     /* file pointer */
    int ii,jj,kk;                 /* generic counters */
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

    for (ii = 0; ii < params.ny; ii++)
    {
        for (jj = 0; jj < params.nx; jj++)
        {
            /* an occupied cell */
            if (cells[ii*params.nx + jj].speeds[0] == -1)
            {
                u_x = u_y = u = 0.0;
                pressure = params.density * c_sq;
            }
            /* no obstacle */
            else
            {
                local_density = 0.0;

                for (kk = 0; kk < NSPEEDS; kk++)
                {
                    local_density += cells[ii*params.nx + jj].speeds[kk];
                }

                /* compute x velocity component */
                u_x = (cells[ii*params.nx + jj].speeds[1] +
                        cells[ii*params.nx + jj].speeds[5] +
                        cells[ii*params.nx + jj].speeds[8]
                    - (cells[ii*params.nx + jj].speeds[3] +
                        cells[ii*params.nx + jj].speeds[6] +
                        cells[ii*params.nx + jj].speeds[7]))
                    / local_density;

                /* compute y velocity component */
                u_y = (cells[ii*params.nx + jj].speeds[2] +
                        cells[ii*params.nx + jj].speeds[5] +
                        cells[ii*params.nx + jj].speeds[6]
                    - (cells[ii*params.nx + jj].speeds[4] +
                        cells[ii*params.nx + jj].speeds[7] +
                        cells[ii*params.nx + jj].speeds[8]))
                    / local_density;

                /* compute norm of velocity */
                u = sqrt((u_x * u_x) + (u_y * u_y));

                /* compute pressure */
                pressure = local_density * c_sq;
            }

            /* write to file */
            fprintf(fp,"%d %d %.12E %.12E %.12E %.12E %d\n",
                jj,ii,u_x,u_y,u,pressure,obstacles[ii*params.nx + jj]);
        }
    }

    fclose(fp);

    fp = fopen(av_vels_file, "w");
    if (fp == NULL)
    {
        DIE("could not open file output file");
    }

    for (ii = 0; ii < params.max_iters; ii++)
    {
        fprintf(fp,"%d:\t%.12E\n", ii, av_vels[ii]/total_cells);
    }

    fclose(fp);
}

float calc_reynolds(const param_t params, float av_vel, unsigned total_cells)
{
    const float viscosity = 1.0 / 6.0 * (2.0 / params.omega - 1.0);

    return (av_vel / (float) total_cells) * params.reynolds_dim / viscosity;
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
