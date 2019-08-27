/*
 * Custom hydro solver as part of Bert's course!
 *
 * Going to attempt to write some C, lol.
 *
 * This a second order (space and time) finite volume solver.
 */

#include <stdio.h>
#include <stdlib.h>

#include "riemann/minmax.h"
#include "riemann/riemann_exact.h"

/* Constants for setting up the state of the system */
#define NUMBER_OF_CELLS 100
#define BOX_SIZE 2.f
#define CELL_VOLUME (BOX_SIZE / NUMBER_OF_CELLS)
#define CELL_MOMENTUM 0.f
#define GAS_GAMMA 1.6666f

/* Left state */
#define LEFT_MASS 0.01f
#define LEFT_ENERGY 0.015f

/* Right state */
#define RIGHT_MASS 0.00125f
#define RIGHT_ENERGY 0.0015f

/* Timestep */
#define TIME_STEP 0.0001f
#define END_TIME 0.2f

/* Slope limiter parameters (cell-wide) */
#define SLOPE_LIMITER_BETA 0.5f

/* Structs for use in the rest of the code */
struct cell {
  /* Cell volume */
  float volume;

  /* Conserved quantities */
  float mass;
  float momentum;
  float energy;

  /* Primitives */
  float density;
  float velocity;
  float pressure;

  /* Primitive gradients */
  float density_gradient;
  float velocity_gradient;
  float pressure_gradient;

  /* Pointer to our friends */
  struct cell *right_neighbour;
  struct cell *left_neighbour;
};

/* Sets up the state of a cell given a left or right state */
void setup_state_single_cell(struct cell *current_cell,
                             struct cell *right_neighbour,
                             struct cell *left_neighbour, int state) {
  /* All cells have the same volume */
  current_cell->volume = CELL_VOLUME;
  current_cell->momentum = CELL_MOMENTUM;

  if (state >= 0) {
    /* Right state */
    current_cell->mass = RIGHT_MASS;
    current_cell->energy = RIGHT_ENERGY;
  } else {
    /* Left state */
    current_cell->mass = LEFT_MASS;
    current_cell->energy = LEFT_ENERGY;
  }

  /* Initialise the primitives for good measure */
  current_cell->density = 0.f;
  current_cell->velocity = 0.f;
  current_cell->pressure = 0.f;

  /* Link the neighbour */
  current_cell->right_neighbour = right_neighbour;
  current_cell->left_neighbour = left_neighbour;
}

/* Converts all conserved variables in a given cell into primitive
 * variables. */
void convert_conserved_to_primitive_single_cell(struct cell *cell) {
  cell->density = cell->mass / cell->volume;
  cell->velocity = cell->momentum / cell->mass;

  /* Total Energy */
  const float energy_density = cell->energy / cell->volume;
  /* Need to take off kinetic energy to recover thermal energy */
  const float kinetic_energy =
      0.5f * cell->density * cell->density * cell->velocity;

  cell->pressure = (GAS_GAMMA - 1.f) * (energy_density - kinetic_energy);
}

/* Converts all primitive variables to conserved variables */
void convert_primitive_to_conserved_single_cell(struct cell *cell) {
  cell->mass = cell->density * cell->volume;
  cell->momentum = cell->velocity * cell->mass;

  /* Thermal energy */
  const float thermal_energy = cell->pressure / (GAS_GAMMA - 1.f);
  /* Need to add on kinetic energy to recover total energy */
  const float kinetic_energy =
      0.5f * cell->density * cell->density * cell->velocity;

  cell->energy = cell->volume * (thermal_energy + kinetic_energy);
}

/* Converts all cells variables from conserved to primitive */
void convert_conserved_to_primitive_all(struct cell *cells) {
  for (int i = 0; i < NUMBER_OF_CELLS; i++) {
    convert_conserved_to_primitive_single_cell(&cells[i]);
  }
}

/* Converts all cells variables from primitive to conserved */
void convert_primitive_to_conserved_all(struct cell *cells) {
  for (int i = 0; i < NUMBER_OF_CELLS; i++) {
    convert_primitive_to_conserved_single_cell(&cells[i]);
  }
}

/* Sets up the contents of all of the cells with periodic boundary
 * conditions */
void setup_all_cells(struct cell *cells) {
  for (int i = 0; i < NUMBER_OF_CELLS; i++) {
    const int right_neighbouring_cell = (i + 1) % NUMBER_OF_CELLS;
    const int left_neighbouring_cell =
        (i + NUMBER_OF_CELLS - 1) % NUMBER_OF_CELLS;
    const int state = i - NUMBER_OF_CELLS / 2;

    setup_state_single_cell(&cells[i], &cells[right_neighbouring_cell],
                            &cells[left_neighbouring_cell], state);
  }

  /* Now need to convert the conserved quantities to primitives before starting
   */
  convert_conserved_to_primitive_all(cells);
}

/* Do the slope limiting step based on neighbour cells */
void slope_limit_cell(struct cell *left_cell, struct cell *cell,
                      struct cell *right_cell) {

  /* Limit density */
  const float density_max =
      max(left_cell->density, right_cell->density) - cell->density;
  const float density_min =
      min(left_cell->density, right_cell->density) - cell->density;
  const float density_gradient_dx = cell->density_gradient * CELL_VOLUME * 0.5f;
  const float density_pred_max =
      max(cell->density + density_gradient_dx,
          cell->density_gradient - density_gradient_dx) -
      cell->density;
  const float density_pred_min =
      min(cell->density + density_gradient_dx,
          cell->density_gradient - density_gradient_dx) -
      cell->density;

  const float density_max_ratio =
      density_pred_max > 0.f ? density_max / density_pred_max : FLT_MAX;
  const float density_min_ratio =
      density_pred_min > 0.f ? density_min / density_pred_min : FLT_MAX;

  const float alpha_density =
      min(1.f, SLOPE_LIMITER_BETA * min(density_max_ratio, density_min_ratio));

  /* Now we can _actually_ limit our cell */
  cell->density_gradient *= alpha_density;

  /* Limit velocity */
  const float velocity_max =
      max(left_cell->velocity, right_cell->velocity) - cell->velocity;
  const float velocity_min =
      min(left_cell->velocity, right_cell->velocity) - cell->velocity;
  const float velocity_gradient_dx =
      cell->velocity_gradient * CELL_VOLUME * 0.5f;
  const float velocity_pred_max =
      max(cell->velocity + velocity_gradient_dx,
          cell->velocity_gradient - velocity_gradient_dx) -
      cell->velocity;
  const float velocity_pred_min =
      min(cell->velocity + velocity_gradient_dx,
          cell->velocity_gradient - velocity_gradient_dx) -
      cell->velocity;

  const float velocity_max_ratio =
      velocity_pred_max > 0.f ? velocity_max / velocity_pred_max : FLT_MAX;
  const float velocity_min_ratio =
      velocity_pred_min > 0.f ? velocity_min / velocity_pred_min : FLT_MAX;

  const float alpha_velocity = min(
      1.f, SLOPE_LIMITER_BETA * min(velocity_max_ratio, velocity_min_ratio));

  /* Now we can _actually_ limit our cell */
  cell->velocity_gradient *= alpha_velocity;

  /* Limit pressure */
  const float pressure_max =
      max(left_cell->pressure, right_cell->pressure) - cell->pressure;
  const float pressure_min =
      min(left_cell->pressure, right_cell->pressure) - cell->pressure;
  const float pressure_gradient_dx =
      cell->pressure_gradient * CELL_VOLUME * 0.5f;
  const float pressure_pred_max =
      max(cell->pressure + pressure_gradient_dx,
          cell->pressure_gradient - pressure_gradient_dx) -
      cell->pressure;
  const float pressure_pred_min =
      min(cell->pressure + pressure_gradient_dx,
          cell->pressure_gradient - pressure_gradient_dx) -
      cell->pressure;

  const float pressure_max_ratio =
      pressure_pred_max > 0.f ? pressure_max / pressure_pred_max : FLT_MAX;
  const float pressure_min_ratio =
      pressure_pred_min > 0.f ? pressure_min / pressure_pred_min : FLT_MAX;

  const float alpha_pressure = min(
      1.f, SLOPE_LIMITER_BETA * min(pressure_max_ratio, pressure_min_ratio));

  /* Now we can _actually_ limit our cell */
  cell->pressure_gradient *= alpha_pressure;
}

/* Generates a drifted state based on the time-step and cell gradients */
void generate_drifted_state(struct cell *cell, float *state) {
  state[0] = cell->density - 0.5f * TIME_STEP *
                                 (cell->density * cell->velocity_gradient +
                                  cell->velocity * cell->density_gradient);
  state[1] = cell->velocity - 0.5f * TIME_STEP *
                                  (cell->velocity * cell->velocity_gradient +
                                   cell->pressure_gradient / cell->density);
  state[2] = 0.f; /* 1D */
  state[3] = 0.f; /* 1D */
  state[4] = cell->pressure -
             0.5f * TIME_STEP *
                 (GAS_GAMMA * cell->pressure * cell->velocity_gradient +
                  cell->velocity * cell->pressure_gradient);
}

/* Perform one full hydro step on all cells */
void step(struct cell *cells) {
  /* Constant ratio for all cells */
  const float gamma_ratio = GAS_GAMMA / (GAS_GAMMA - 1.f);

  /* Convert all variables from primitive to conserved */
  convert_primitive_to_conserved_all(cells);

  for (int i = 0; i < NUMBER_OF_CELLS; i++) {
    struct cell *current_cell = &cells[i];
    struct cell *right_neighbour_cell = cells[i].right_neighbour;
    struct cell *left_neighbour_cell = cells[i].left_neighbour;

    /* Compute gradients */
    const float grad_fac = 0.5 / CELL_VOLUME; /* Will change in ND! */
    const float density_gradient = grad_fac * (right_neighbour_cell->density -
                                               left_neighbour_cell->density);
    const float velocity_gradient = grad_fac * (right_neighbour_cell->velocity -
                                                left_neighbour_cell->velocity);
    const float pressure_gradient = grad_fac * (right_neighbour_cell->pressure -
                                                left_neighbour_cell->velocity);

    /* Whack 'em in */
    current_cell->density_gradient = density_gradient;
    current_cell->velocity_gradient = velocity_gradient;
    current_cell->pressure_gradient = pressure_gradient;
  }

  /* Now perform the slope limiting step */
  for (int i = 0; i < NUMBER_OF_CELLS; i++) {
    struct cell *current_cell = &cells[i];
    struct cell *right_neighbour_cell = cells[i].right_neighbour;
    struct cell *left_neighbour_cell = cells[i].left_neighbour;

    slope_limit_cell(left_neighbour_cell, current_cell, right_neighbour_cell);
  }

  /* We have now calculated the gradients for each of the cells, so we
   * can now do the actual Riemann problem */

  for (int i = 0; i < NUMBER_OF_CELLS; i++) {
    struct cell *current_cell = &cells[i];
    /* Always go from left to right */
    struct cell *neighbour_cell = cells[i].right_neighbour;

    /* Predict half step forward */

    float left_state[5];
    generate_drifted_state(current_cell, &left_state[0]);

    float right_state[5];
    generate_drifted_state(neighbour_cell, &right_state[0]);

    float output_state[5];

    /* Needed for the SWIFT riemann solver; the unit vector for the face*/
    float n_unit[3];
    n_unit[0] = 1.f;
    n_unit[1] = 0.f;
    n_unit[2] = 0.f;

    /* Call our favourite riemann solver */
    riemann_solver_solve(&left_state[0], &right_state[0], &output_state[0],
                         &n_unit[0]);

    /* Solution values */
    float density_S = output_state[0];
    float velocity_S = output_state[1];
    float pressure_S = output_state[4];

    /* Now we can use these solutions to generate the fluxes */
    const float kinetic_energy = 0.5f * density_S * velocity_S * velocity_S;

    const float flux_mass = density_S * velocity_S;
    const float flux_momentum = 2.f * kinetic_energy + pressure_S;
    const float flux_energy =
        ((pressure_S * gamma_ratio) + kinetic_energy) * velocity_S;

    /* Now we can do the flux exchange; note here surface area is taken to be
     * 1.0 for all particles */
    current_cell->mass -= flux_mass * TIME_STEP;
    current_cell->momentum -= flux_momentum * TIME_STEP;
    current_cell->energy -= flux_energy * TIME_STEP;

    /* In one ear and right out the other */
    neighbour_cell->mass += flux_mass * TIME_STEP;
    neighbour_cell->momentum += flux_momentum * TIME_STEP;
    neighbour_cell->energy += flux_energy * TIME_STEP;

    /* We're done! */
  }

  /* Convert all variables from conserved to primitive */
  convert_conserved_to_primitive_all(cells);
}

/* Prints the current state of all cells */
void print_state_of_cells(struct cell *cells) {
  /* Print a nice header so people can actually read the data */
  printf("# Volume, Mass, Momentum, Energy, Density, Velocity, Pressure, "
         "Density Gradient, Velocity Gradient, Pressure Gradient\n");
  for (int i = 0; i < NUMBER_OF_CELLS; i++) {
    printf("%e, %e, %e, %e, %e, %e, %e, %e, %e, %e\n", cells[i].volume,
           cells[i].mass, cells[i].momentum, cells[i].energy, cells[i].density,
           cells[i].velocity, cells[i].pressure, cells[i].density_gradient,
           cells[i].velocity_gradient, cells[i].pressure_gradient);
  }
}

int main(int argc, char *argv[]) {
  /* Welcome! */
  printf("# Welcome to your dodgy hydro code!\n");
  printf("# We are going to integrate a SodShock over some time!\n");

  printf("# Setting up initial conditions...\n");
  struct cell *cells = malloc(NUMBER_OF_CELLS * sizeof(struct cell));
  setup_all_cells(cells);

  // printf("# Printing initial state of all cells:\n");
  // print_state_of_cells(cells);

  /* Now we can actually do the iteration! */
  printf("# Running iterations. Wish me luck!\n");
  const int number_of_steps = END_TIME / TIME_STEP;
  float current_time = 0.f;

  for (int current_step = 0; current_step < number_of_steps; current_step++) {
    step(cells);
    current_time += TIME_STEP;
  }

  printf("# Printing final state of all cells:\n");
  print_state_of_cells(cells);

  printf("# Reached time t=%e.\n", current_time);
  printf("# Bye!\n");

  return 0;
}
