#ifndef LBM_HDR_FILE
#define LBM_HDR_FILE

#define NSPEEDS         9

/* Size of box in imaginary 'units */
#define BOX_X_SIZE (100.0)
#define BOX_Y_SIZE (100.0)

/* struct to hold the parameter values */
typedef struct {
    int nx;            /* no. of cells in x-direction */
    int ny;            /* no. of cells in y-direction */
    int max_iters;      /* no. of iterations */
    int reynolds_dim;  /* dimension for Reynolds number */
    float density;       /* density per link */
    float accel;         /* density redistribution */
    float omega;         /* relaxation parameter */
} param_t;

/* obstacle positions */
typedef struct {
    float obs_x_min;
    float obs_x_max;
    float obs_y_min;
    float obs_y_max;
} obstacle_t;

/* struct to hold the 'speed' values */
typedef struct {
    float speeds[NSPEEDS];
} speed_t;

typedef struct {
    float *spd;
} speed_t2;

typedef enum { ACCEL_ROW, ACCEL_COLUMN } accel_e;
typedef struct {
    accel_e col_or_row;
    int idx;
} accel_area_t;

/* Parse command line arguments to get filenames */
void parse_args (int argc, char* argv[],
    char** final_state_file, char** av_vels_file, char** param_file);

void initialise(const char* paramfile, accel_area_t * accel_area,
    param_t* params, float** cells_ptr, float** tmp_cells_ptr,
    unsigned** obstacles_ptr, float** av_vels_ptr, unsigned* total_cells);

void write_values(const char * final_state_file, const char * av_vels_file,
    const param_t params, float* cells, unsigned* obstacles, float* av_vels, unsigned total_cells);

void finalise(float** cells_ptr, float** tmp_cells_ptr,
    unsigned** obstacles_ptr, float** av_vels_ptr);

float timestep(const param_t params, const accel_area_t accel_area,
    float* cells, float* tmp_cells, int* obstacles, int total_cells);
void accelerate_flow(const param_t params, const accel_area_t accel_area,
    float* cells, int* obstacles);
void propagate(const param_t params, speed_t* cells, speed_t* tmp_cells);
void rebound(const param_t params, speed_t* cells, speed_t* tmp_cells, int* obstacles);
float collision(const param_t params, float* cells, float* tmp_cells, int* obstacles, int total_cells);

/* Sum all the densities in the grid.
** The total should remain constant from one timestep to the next. */
float total_density(const param_t params, speed_t* cells);

/* compute average velocity */
float av_velocity(const param_t params, float* cells, unsigned* obstacles);

/* calculate Reynolds number */
float calc_reynolds(const param_t params, float* cells, unsigned* obstacles);

/* Exit, printing out formatted string */
#define DIE(...) exit_with_error(__LINE__, __FILE__, __VA_ARGS__)
void exit_with_error(int line, const char* filename, const char* format, ...)
__attribute__ ((format (printf, 3, 4)));

#endif
