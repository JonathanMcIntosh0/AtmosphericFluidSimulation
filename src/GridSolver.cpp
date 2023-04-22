#define _USE_MATH_DEFINES
#include <iostream>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <utility>

#include "GridSolver.h"

#define TOP_TEMP (-1)
#define BOT_TEMP 1
#define EPSILON 0.0001f

using std::vector;
using std::swap;
using std::deque;
using std::shared_ptr;
using std::make_shared;

float nv = 1.0f, kappa = 0.1f; // diffuse rate. nv is for velocity field and kapp is for heat transfer.
float dt = 0.01f; // simulation step size
float simTime;
float dx = 1; // size of a cell in simulation
//float dz = 1; // distance between vertical layers in simulation
float velScale = 10.0f; // scaling for display
float buoyancy = 0.5f; // buoyancy force coefficient
float planetary_rotation = 2.0f; // Rotation = 2 for videos (1 for external heat vids, then 2 for last one)
float external_heat_factor = 1.0f; // .1, .5, 1.
int iterations = 20; // number of iterations for iterative solvers (Gauss–Seidel, Conjugate Gradients, etc.)

/*
 * Quantities in the grid.
 * Each has 2 entries (e.g., vx[0] and vy[1] for storing the current/next state)
 * But you had better use pointers like cur_vx, next_vx, etc. defined below.
 * Note that rows are horizontal and columns are vertical.
 * e.g., cur_vx[3][4] is the x component of the velocity at position (4, 3).
 * Also note that the first rows/columns and last rows/columns are for the boundary.
 */
// velocity
GridArray<double> vx [2] { GridArray<double>(), GridArray<double>() };
GridArray<double> vy [2] { GridArray<double>(), GridArray<double>() };
GridArray<double> vz [2] { GridArray<double>(), GridArray<double>() };

// temperature
GridArray<double> tp [2] { GridArray<double>(), GridArray<double>() };

// pointers for the corresponding quantities in the current and the next states
GridArray<double> *cur_vx = &vx[0], *next_vx = &vx[1];
GridArray<double> *cur_vy = &vy[0], *next_vy = &vy[1];
GridArray<double> *cur_vz = &vz[0], *next_vz = &vz[1];
GridArray<double> *cur_tp = &tp[0], *next_tp = &tp[1];

// heat source
GridArray<double> tpSrc;
//double tpAvg[lats];

/*
 * Flow type:
 *  zonal: along latitude lines (East-West)
 *	meridional: along longitude lines (North-South)
 *	vertical: vertical flow
 *	other: quantities like temperature, density, etc.
 */
enum FlowType { zonal, meridional, vertical, temperature, other };
void setBnd(GridArray<double> &a, FlowType flowType);

//Vertex gridVertices[1][1][1];
Array3D<Vertex, longs + 1, lats + 1, layers + 1> gridVertices;
std::vector<std::vector<GLuint>> gridIndices_long(longs + 1, std::vector<GLuint>(6 * lats * layers));
std::vector<std::vector<GLuint>> gridIndices_lat(lats + 1, std::vector<GLuint>(6 * longs * layers));
std::vector<std::vector<GLuint>> gridIndices_vert(layers + 1, std::vector<GLuint>(6 * longs * lats));

Array3D<Vertex[2], longs, lats, layers> velVertices;
std::vector<std::vector<GLuint>> velIndices_long(longs, std::vector<GLuint>(2 * lats * layers));
std::vector<std::vector<GLuint>> velIndices_lat(lats, std::vector<GLuint>(2 * longs * layers));
std::vector<std::vector<GLuint>> velIndices_vert(layers, std::vector<GLuint>(2 * longs * lats));

void initGrid() {
    vx[0].fill(0); vx[1].fill(0);
    vy[0].fill(0); vy[1].fill(0);
    vz[0].fill(0); vz[1].fill(0);
    tp[0].fill(0); tp[1].fill(0);
    
    tpSrc.fill(0);
    
    cur_vx = &vx[0];
    cur_vy = &vy[0];
    cur_vz = &vz[0];
    cur_tp = &tp[0];
    next_vx = &vx[1];
    next_vy = &vy[1];
    next_vz = &vz[1];
    next_tp = &tp[1];

    //Set tpAvg (linear function from -1 to 1)
//    tpAvg[0] = -1.0;
//    double incr_size = 2.0 / (lats - 1);
//    for (size_t j = 1; j < lats; ++j) {
//        tpAvg[j] = tpAvg[j - 1] + incr_size;
//    }

    //Set tpSrc
//    tpSrc[longs/2][lats/2][1] = 50;
    for (size_t j = 1; j <= lats; ++j) {
        double ext_heat = -cos(2 * (j - .5) * M_PI / lats);
        for (size_t l = 1; l <= layers; ++l) {
            for (size_t i = 1; i <= longs; ++i) {
                tpSrc.at(i, j, l) = external_heat_factor * ext_heat;
            }
        }
    }

    // Init temperature buffers
    Array3D<GLuint, longs + 1, lats + 1, layers + 1> indices;
    GLuint idx = 0;
    for (size_t i = 0; i <= longs; ++i) {
        for (size_t j = 0; j <= lats; ++j) {
            for (size_t l = 0; l <= layers; ++l) {
                // TODO possibly deal with flip of coords??
                gridVertices.at(i, j, l) = Vertex((float) i * cellSize, (float) j * cellSize, (float) l * cellSize, 0.f,
                                               0.f, 0.f, 1.f);
                indices.at(i, j, l) = idx++;
            }
        }
    }

    for (size_t i = 0; i <= longs; ++i) {
        for (size_t j = 0; j <= lats; ++j) {
            for (size_t l = 0; l <= layers; ++l) {
                size_t tmp;
                if (j != lats && l != layers) {
                    tmp = 6 * (j * layers + l);
                    gridIndices_long[i][tmp++] = indices.at(i, j, l);
                    gridIndices_long[i][tmp++] = indices.at(i, j + 1, l);
                    gridIndices_long[i][tmp++] = indices.at(i, j + 1, l + 1);
                    gridIndices_long[i][tmp++] = indices.at(i, j, l);
                    gridIndices_long[i][tmp++] = indices.at(i, j, l + 1);
                    gridIndices_long[i][tmp++] = indices.at(i, j + 1, l + 1);
                }
                if (i != longs && l != layers) {
                    tmp = 6 * (i * layers + l);
                    gridIndices_lat[j][tmp++] = indices.at(i, j, l);
                    gridIndices_lat[j][tmp++] = indices.at(i + 1, j, l);
                    gridIndices_lat[j][tmp++] = indices.at(i + 1, j, l + 1);
                    gridIndices_lat[j][tmp++] = indices.at(i, j, l);
                    gridIndices_lat[j][tmp++] = indices.at(i, j, l + 1);
                    gridIndices_lat[j][tmp++] = indices.at(i + 1, j, l + 1);
                }
                if (i != longs && j != lats) {
                    tmp = 6 * (i * lats + j);
                    gridIndices_vert[l][tmp++] = indices.at(i, j, l);
                    gridIndices_vert[l][tmp++] = indices.at(i + 1, j, l);
                    gridIndices_vert[l][tmp++] = indices.at(i + 1, j + 1, l);
                    gridIndices_vert[l][tmp++] = indices.at(i, j, l);
                    gridIndices_vert[l][tmp++] = indices.at(i, j + 1, l);
                    gridIndices_vert[l][tmp++] = indices.at(i + 1, j + 1, l);
                }
            }
        }
    }

    // Init velocity buffers
//    GLuint indicesVel[longs][lats][layers][2];
    idx = 0;
    for (size_t i = 0; i < longs; ++i) {
        for (size_t j = 0; j < lats; ++j) {
            for (size_t l = 0; l < layers; ++l) {
                velVertices.at(i, j, l)[0] = Vertex(
                        (float) (i + 0.5) * cellSize,
                        (float) (j + 0.5) * cellSize,
                        (float) (l + 0.5) * cellSize,
                        0.f, 1.f, 0.f, 1.f);
                velVertices.at(i, j, l)[1] = Vertex(
                        (float) (i + 0.5) * cellSize,
                        (float) (j + 0.5) * cellSize,
                        (float) (l + 0.5) * cellSize,
                        0.5f, .5f, 0.5f, 1.f);

                // (lats - 1) * layers + layers - 1 = lats * layers - 1
                size_t tmp = 2 * (j * layers + l);
                velIndices_long[i][tmp] = idx;
                velIndices_long[i][tmp + 1] = idx + 1;

                tmp = 2 * (i * layers + l);
                velIndices_lat[j][tmp] = idx;
                velIndices_lat[j][tmp + 1] = idx + 1;

                tmp = 2 * (i * lats + j);
                velIndices_vert[l][tmp] = idx;
                velIndices_vert[l][tmp + 1] = idx + 1;

                idx += 2;

//                indicesVel[i][j][l][0] = idx++;
//                indicesVel[i][j][l][1] = idx++;
            }
        }
    }

//    for (size_t i = 0; i < longs; ++i) {
//        for (size_t j = 0; j < lats; ++j) {
//            for (size_t l = 0; l < layers; ++l) {
//                size_t tmp = 2 * j * l;
//                velIndices_long[i][tmp] = indicesVel[i][j][l][0];
//                velIndices_long[i][tmp + 1] = indicesVel[i][j][l][1];
//
//                tmp = 2 * i * l;
//                velIndices_lat[j][tmp] = indicesVel[i][j][l][0];
//                velIndices_lat[j][tmp + 1] = indicesVel[i][j][l][1];
//
//                tmp = 2 * i * j;
//                velIndices_vert[l][tmp] = indicesVel[i][j][l][0];
//                velIndices_vert[l][tmp + 1] = indicesVel[i][j][l][1];
//            }
//        }
//    }

	updateGrid();
	simTime = 0;
}

/*
 * Unit square bilinear-interpolation (https://en.wikipedia.org/wiki/Bilinear_interpolation#On_the_unit_square).
 * You need this for interpolating the value in the cell (e.g., in backtracing advect).
 * Parameters:
 *	x, y:
 *	    position of the point whose value you want to interpolate
 *	v00:
 *	    value at (0, 0)
 *	v01:
 *	    value at (0, 1)
 *	v10:
 *	    value at (1, 0)
 *	v11:
 *	    value at (1, 1)
 *  return:
 *	interpolated value
 */
double bilinearInterpolate(double x, double y, double v00, double v01, double v10, double v11) {
	// DONE: Do a unit square bilinear interpolation for a point at (x, y).
	return v00 * (1 - x) * (1 - y) + v10 * x * (1 - y) + v01 * (1 - x) * y + v11 * x * y;
}

/*
 * Set up the boundary.
 * Parameters:
 *	a:
 *	    2D array whose boundary need to be set
 *	flowType:
 *	    type of flow: horizontal flow (the interpolated quantity will goes to zero at the vertical boundaries),
 *	    vertical flow (the interpolated quantity will goes to zero at the horizontal boundaries),
 *	    other (e.g., temperature, density, which are not really flow... only continuity need to be guaranteed.)
 */
void setBnd(GridArray<double> &a, FlowType flowType) {
	// DONE: Set up the boundary according to the flow type.

    // TOP and BOTTOM boundaries
    for (int i = 1; i <= longs; ++i)
        for (int j = 1; j <= lats; ++j)
        {
            a.at(i, j, 0) = flowType == temperature ? 2*(BOT_TEMP) - a.at(i, j, 1) :
                         flowType == vertical ? -a.at(i, j, 1) :
                         a.at(i, j, 1);
//                         flowType == other ? a.at(i, j, 1) :
//                         -a.at(i, j, 1); // No slip condition
            a.at(i, j, layers + 1) = flowType == temperature ? 2*(TOP_TEMP) - a.at(i,j,layers) :
                                  flowType == vertical ? -a.at(i,j,layers) :
                                  a.at(i,j,layers);
        }
    // NORTH and SOUTH boundaries
    for (int i = 1; i <= longs; ++i) {
        for (int l = 1; l <= layers; ++l)
        {
            a.at(i,0,l) = flowType == meridional ? -a.at(i,1,l) : a.at(i,1,l);
            a.at(i, lats + 1, l) = flowType == meridional ? -a.at(i,lats,l) : a.at(i,lats,l);
        }
    }

    // NORTH and SOUTH corner zonal lines
    for (int i = 1; i <= longs; ++i) {
        a.at(i, 0, 0) = 0.5*(a.at(i, 1, 0) + a.at(i, 0, 1));
        a.at(i, 0, layers + 1) = 0.5 * (a.at(i, 1, layers + 1) + a.at(i, 0, layers));
        a.at(i, lats + 1, 0) = 0.5 * (a.at(i, lats, 0) + a.at(i, lats + 1, 1));
        a.at(i, lats + 1, layers + 1) = 0.5 * (a.at(i, lats, layers + 1) + a.at(i, lats + 1, layers));
    }

    // EAST and WEST boundaries TODO check that other boundaries with wrap-around part is fine
    for (int j = 0; j <= lats + 1; ++j)
        for (int l = 0; l <= layers + 1; ++l)
        {
//            a[0][i][j] = flowType == other ? a[1][i][j] : a[longs][i][j];
//            a[longs + 1][i][j] = flowType == other ? a[longs][i][j] : a[1][i][j];
            a.at(0, j, l) = a.at(longs, j, l);
            a.at(longs + 1, j, l) = a.at(1, j, l);
        }
}

/*
 * Add source to the field.
 * Parameters:
 *	src:
 *	    2D array containing the source
 *	a:
 *	    target 2D array storing the quantity to modify
 */
void addSource(GridArray<double> src, GridArray<double> &a) {
	// DONE: Add the source from array *src* to the target array *a*: e.g., a = a + dt * src
    for (int i = 0; i < longs + 2; ++i) {
        for (int j = 0; j < lats + 2; ++j) {
            for (int l = 0; l < layers + 2; ++l) {
                a.at(i,j,l) += dt * src.at(i,j,l);
            }
        }
    }
}

/*
 * Compute the diffusion part.
 * Parameters:
 *	a0:
 *	    2D array storing the quantities in current state
 *	a1:
 *	    2D array storing the quantities in next state
 *	nv:
 *	    diffusion rate (nv/kappa in the equations)
 *	flowType:
 *	    flow type
 */
void diffuse(GridArray<double> a0, GridArray<double> &a1, double nv, FlowType flowType) {
	// DONE: diffusion
	// Compute the diffusion part and update array *a1*.
	// Use dt and dx for step size and cell size.
	// Do a implicit solve for stability.
	// Use a Gauss Seidel solve (or better, Conjugate Gradients).
	// Use *iterations* as the number of iterations.
	// Call setBnd to fix the boundary.
    bool hasConverged;
    double max_diff;

    double a = dt * nv / (dx * dx);
    for (int k = 0; k < iterations; k++) {
        hasConverged = true;
        max_diff = 0;
        for (int i = 1; i <= longs; i++) {
            for (int j = 1; j <= lats; j++) {
                for (int l = 1; l <= layers; l++) {
                    double b = (a0.at(i,j,l) +
                                a * (a1.at(i-1,j,l) + a1.at(i+1,j,l)
                                     + a1.at(i,j-1,l) + a1.at(i,j+1,l)
                                     + a1.at(i,j,l-1) + a1.at(i,j,l+1))) / (1 + 6*a);
                    double diff = abs(b - a1.at(i,j,l));
                    if (diff > EPSILON) {
                        hasConverged = false;
                        max_diff = fmax(max_diff, diff);
                    }

                    a1.at(i,j,l) = b;
                }
            }
        }
        setBnd(a1, flowType);
        if (hasConverged) break;
    }
    if (!hasConverged) std::cout << "DIFFUSE: NO CONVERGENCE: DIFF = " << max_diff << std::endl;


}

/*
 * Compute the advection part.
 * Parameters:
 *	a0:
 *	    2D array storing the quantities in current state
 *	a1:
 *	    2D array storing the quantities in next state
 *	vx, vy:
 *	    2D arrays storing the velocity field
 *	flowType:
 *	    flow type
 */
void advect(GridArray<double> a0, GridArray<double> &a1,
            GridArray<double> vx, GridArray<double> vy, GridArray<double> vz,
            FlowType flowType) {
	// DONE: advection
	// Compute the advection part and update array *a1*.
	// Use dt and dx for step size and cell size.
	// Do a linear (or better, higher order or adaptive) backtrace for each center of the cells.
	// Compute the quantity in the previous state using the bilinear interpolation.
	// Call setBnd to fix the boundary.

    float dt0 = dt / dx;
    for (int i = 1; i <= longs; i++) {
        for (int j = 1; j <= lats; j++) {
            for (int l = 1; l <= layers; l++) {
                double x = i - dt0 * vx.at(i,j,l);
                double y = j - dt0 * vy.at(i,j,l);
                double z = l - dt0 * vz.at(i,j,l);

                // Deal with wrap around for EAST-WEST (0 = longs)
                if (x < .5)
                    x += longs;
                else if (x > longs + .5)
                    x -= longs;
                int i0 = (int) x;
                int i1 = i0 + 1;

                if (y < 0.5) y = 0.5;
                else if (y > lats + 0.5) y = lats + 0.5;
                int j0 = (int) y;
                int j1 = j0 + 1;

                if (z < 0.5) z = 0.5;
                else if (z > layers + 0.5) z = layers + 0.5;
                int l0 = (int) z;
                int l1 = l0 + 1;

                z -= l0; // z is now between 0 and 1 (0 at l0, 1 at l1)
                double v00 = (1 - z) * a0.at(i0, j0, l0) + z * a0.at(i0, j0, l1);
                double v01 = (1 - z) * a0.at(i0, j1, l0) + z * a0.at(i0, j1, l1);
                double v10 = (1 - z) * a0.at(i1, j0, l0) + z * a0.at(i1, j0, l1);
                double v11 = (1 - z) * a0.at(i1, j1, l0) + z * a0.at(i1, j1, l1);

                a1.at(i,j,l) = bilinearInterpolate(x - i0, y - j0, v00, v01, v10, v11);
            }
        }
    }
    setBnd(a1, flowType);
}

GridArray<double> s_p;
Array3D<double, longs, lats, layers> s_div;

/*
 * Projection for the mass conservation.
 * Parameter:
 *	vx, vy:
 *	    the velocity field to be fixed
 */
void project(GridArray<double> &vx, GridArray<double> &vy, GridArray<double> &vz) {
	// DONE: projection
	// Do a Poisson Solve to get a divergence free velocity field.
	// Use a Gauss Seidel solve (or better, Conjugate Gradients).
	// Use *iterations* as the number of iterations.
	// Call setBnd to fix the boundary.
    for (int i = 1; i <= longs; i++)
        for (int j = 1; j <= lats; j++) {
            for (int l = 1; l <= layers; l++) {
                s_div.at(i-1, j-1, l-1) =
                        -0.5*dx*(
                                vx.at(i+1,j,l) - vx.at(i-1,j,l)
                                + vy.at(i,j+1,l) - vy.at(i,j-1,l)
                                + vz.at(i,j,l+1) - vz.at(i,j,l-1));
                s_p.at(i,j,l) = 0;
            }
        }
//    setBnd(s_div, other);
    setBnd(s_p, other);

    bool hasConverged;
    double max_diff;
    for (int k = 0; k < iterations; k++) {
        hasConverged = true;
        max_diff = 0;
        for (int i = 1; i <= longs; i++) {
            for (int j = 1; j <= lats; j++) {
                for (int l = 1; l <= layers; l++) {
                    double p = (s_div.at(i-1, j-1, l-1)
                            + s_p.at(i-1,j,l) + s_p.at(i+1,j,l)
                            + s_p.at(i,j-1,l) + s_p.at(i,j+1,l)
                            + s_p.at(i,j,l-1) + s_p.at(i,j,l+1))/6;
                    double diff = abs(p - s_p.at(i,j,l));
                    if (diff > EPSILON) {
                        hasConverged = false;
                        max_diff = fmax(max_diff, diff);
                    }
                    s_p.at(i,j,l) = p;
                }
            }
        }
        setBnd(s_p, other);
        if (hasConverged) break;
    }
    if (!hasConverged) std::cout << "PROJECT: NO CONVERGENCE: DIFF = " << max_diff << std::endl;

    for (int i = 1; i <= longs; i++) {
        for (int j = 1; j <= lats; j++) {
            for (int l = 1; l <= layers; l++) {
                vx.at(i,j,l) -= 0.5*(s_p.at(i + 1, j, l) - s_p.at(i - 1, j, l))/dx;
                vy.at(i,j,l) -= 0.5*(s_p.at(i, j + 1, l) - s_p.at(i, j - 1, l))/dx;
                vz.at(i,j,l) -= 0.5*(s_p.at(i, j, l + 1) - s_p.at(i, j, l - 1))/dx;
            }
        }
    }
    setBnd(vx, zonal); setBnd(vy, meridional); setBnd(vz, vertical);
}

/*
 * Get the reference (average) temperature of the grid.
 * Parameters:
 *	temperature:
 *	    2D array storing the temperatures
 */
double getReferenceTemperature(GridArray<double> tp) {
	// DONE: Sum up array *temperature* and compute the average.
    // Note: assuming we don't include the boundaries in the average...
    double sum = 0;
    for (int i = 1; i <= longs; i++) {
        for (int j = 1; j <= lats; j++) {
            for (int l = 1; l <= layers; l++) {
                sum += tp.at(i, j, l);
            }
        }
    }
	return sum / (longs * lats * layers);
}

/*
 * Apply the buoyancy force due to the temperature difference.
 * Parameters:
 *	a:
 *	    2D array storing the velocity
 *	tp:
 *	    2D array storing the temperature
 *	beta:
 *	    buoyancy force coefficient
 *	flowType:
 *	    flow type
 */
void applyTemperatureForce(GridArray<double> &a, GridArray<double> tp, double beta, FlowType flowType) {
	// DONE: buoyancy forces
	// Apply the buoyancy force and update array *a*.
	// For more details, see Foster and Metaxas [1997] Equation 2.

    if (flowType != vertical) return;
    double t_ref = getReferenceTemperature(tp);
//    double t_ref = 0;
    for (int j = 1; j <= lats; j++) {
        for (int i = 1; i <= longs; i++) {
            for (int l = 1; l <= layers; l++) {
                a.at(i,j,l) += dt * beta * (tp.at(i,j,l) - t_ref);
//                a[i][j] += dt * beta * (tp[i][j] - tpAvg[j - 1]);
//                a[i][j] += dt * beta * (tp[i][j]);
            }
        }
    }
    setBnd(a, flowType);
}

void applyCoriolisForce(GridArray<double> a0, GridArray<double> &a1, GridArray<double> vx, GridArray<double> vy,
                        double omega, FlowType flowType) {
    // Apply the coriolis force and update array *a*.
    // F = (fv, -fu, ~0) where f = 2*omega*sin(lat)

    if (flowType != zonal && flowType != meridional) return;
    for (int j = 1; j <= lats; j++) {
        double f = 2 * omega * cos((j - .5) * M_PI / lats);
        for (int i = 1; i <= longs; i++) {
            for (int l = 1; l <= layers; l++) {
                a1.at(i,j,l) = a0.at(i,j,l) + dt * f * ((flowType == zonal) ? vy.at(i,j,l) : -vx.at(i,j,l));
            }
        }
    }
    setBnd(a1, flowType);
}

/*
 * One stimulation step.
 */
void step() {
	// DONE: step the simulation
	// Compute and update the velocity and temperature (in cur_vx,cur_vy, cur_tp) in the next state based on the current state.
	// You need to apply source, diffuse, advect forces for temperature and velocity fields.
	// For velocity field, you need to do the projection to get a divergence free field.
	// You also need to apply the buoyancy force to the velocity field.
	// Don't forget to swap pointers (e.g., cur_vx and next_vx, etc.)!
	// Have a look at GDC03.
	// Change the paramters (dt, dx, etc.) in the setting panel to check whether your solver is stable!


    // Vel_step
//    addSource(vxSrc, cur_vx); addSource(vySrc, cur_vy);

    applyTemperatureForce(*cur_vz, *cur_tp, buoyancy, vertical);
    applyCoriolisForce(*cur_vx, *next_vx, *cur_vx, *cur_vy, planetary_rotation, zonal);
    applyCoriolisForce(*cur_vy, *next_vy, *cur_vx, *cur_vy, planetary_rotation, meridional);
    swap(cur_vx, next_vx); swap(cur_vy, next_vy);

    diffuse(*cur_vx, *next_vx, nv, zonal);
    diffuse(*cur_vy, *next_vy, nv, meridional);
    diffuse(*cur_vz, *next_vz, nv, vertical);
    swap(cur_vx, next_vx); swap(cur_vy, next_vy); swap(cur_vz, next_vz);

    project(*cur_vx, *cur_vy, *cur_vz); // Project to improve result of advect
    advect(*cur_vx, *next_vx, *cur_vx, *cur_vy, *cur_vz, zonal);
    advect(*cur_vy, *next_vy, *cur_vx, *cur_vy, *cur_vz, meridional);
    advect(*cur_vz, *next_vz, *cur_vx, *cur_vy, *cur_vz, vertical);
    project(*next_vx, *next_vy, *next_vz);
    swap(cur_vx, next_vx); swap(cur_vy, next_vy); swap(cur_vz, next_vz);

    // Temp_step
    addSource(tpSrc, *cur_tp);
    diffuse(*cur_tp, *next_tp, kappa, temperature);
    swap(cur_tp, next_tp);
    advect(*cur_tp, *next_tp, *cur_vx, *cur_vy, *cur_vz, temperature);
//    advect(cur_tp, next_tp, next_vx, next_vy, other);
    swap(cur_tp, next_tp);

	// Please DO NOT change the following
	simTime += dt;
	//updateGrid();
//	updateFilament();
//	memset(vxSrc, 0, sizeof(vxSrc));
//	memset(vySrc, 0, sizeof(vySrc));

//    memset(tpSrc, 0, sizeof(tpSrc));
}

void updateGrid() {
    for (size_t l = 0; l <= layers; ++l)
        for (size_t i = 0; i <= longs; ++i)
            for (size_t j = 0; j <= lats; ++j) {
                double v00 = .5*(cur_tp->at(i,j,l) + cur_tp->at(i, j, l + 1));
                double v01 = .5*(cur_tp->at(i, j + 1, l) + cur_tp->at(i, j + 1, l + 1));
                double v10 = .5*(cur_tp->at(i + 1, j, l) + cur_tp->at(i + 1, j, l + 1));
                double v11 = .5*(cur_tp->at(i + 1, j + 1, l) + cur_tp->at(i + 1, j + 1, l + 1));
                // Weird interpolation because the quad is rendered as 2 triangles?
                // Maybe this can help: https://jcgt.org/published/0011/03/04/paper.pdf
                double t = bilinearInterpolate(0.5, 0.5, v00, v01, v10, v11);
                gridVertices.at(i, j, l).r = 0;
                gridVertices.at(i, j, l).g = 0;
                gridVertices.at(i, j, l).b = 0;
                if (t > 0) {
                    gridVertices.at(i, j, l).r = (GLfloat)std::min(1.0, t);
                } else if (t < 0) {
                    gridVertices.at(i, j, l).b = (GLfloat)std::min(1.0, -t);
                }
                gridVertices.at(i, j, l).a = (abs(t) < .1f) ? 0 : abs(t);

            }

    for (size_t l = 0; l < layers; ++l)
        for (size_t i = 0; i < longs; ++i)
            for (size_t j = 0; j < lats; ++j) {
                // TODO deal with flip of coordinates
                velVertices.at(i, j, l)[1].x = velVertices.at(i, j, l)[0].x + cur_vx->at(i + 1, j + 1, l + 1) * velScale;
                velVertices.at(i, j, l)[1].y = velVertices.at(i, j, l)[0].y + cur_vy->at(i + 1, j + 1, l + 1) * velScale;
                velVertices.at(i, j, l)[1].z = velVertices.at(i, j, l)[0].z + cur_vz->at(i + 1, j + 1, l + 1) * velScale;

//                if (abs(cur_vx[i + 1][j + 1][l + 1]) < EPSILON
//                    && abs(cur_vy[i + 1][j + 1][l + 1]) < EPSILON
//                    && abs(cur_vz[i + 1][j + 1][l + 1]) < EPSILON) {
//                    velVertices[i][j][l][0].a = 0;
//                    velVertices[i][j][l][1].a = 0;
//                } else {
//                    velVertices[i][j][l][0].a = 1;
//                    velVertices[i][j][l][1].a = 1;
//                }

            }
}

//double dist(const Vertex &a, const Vertex &b) {
//    double x = a.x - b.x;
//    double y = a.y - b.y;
//    return std::sqrt(x*x + y*y);
//}

//void addFilament(double x, double y) {
//	filaments.push_back(make_shared<vector<Vertex>>());
//	for (size_t i = 0; i <= longs; ++i) {
//		filaments.back()->push_back(Vertex(x, (float)i * cellSize, 0.f, 1.0f, 1.0f, 1.0f, 1.0f));
//	}
//	ages.push_back(0);
//
//	filaments.push_back(make_shared<vector<Vertex>>());
//	for (size_t i = 0; i <= lats; ++i) {
//		filaments.back()->push_back(Vertex((float)i * cellSize, y, 0.f, 1.0f, 1.0f, 1.0f, 1.0f));
//	}
//	ages.push_back(0);
//}

//void updateFilament() {
//	while (!ages.empty() && ages.front() > maxAge) {
//		ages.pop_front();
//		filaments.pop_front();
//	}
//	for (size_t i = 0; i < filaments.size(); ++i) {
//		ages[i] += dt;
//		shared_ptr<vector<Vertex>> tmp = make_shared<vector<Vertex>>();
//		tmp->push_back(filaments[i]->front());
//		for (size_t j = 1; j < filaments[i]->size(); ++j) {
//			const Vertex &a = tmp->back();
//			const Vertex &b = filaments[i]->at(j);
//			if (dist(a, b) > refinementThreshold) {
//				tmp->push_back(Vertex((a.x + b.x)/2, (a.y + b.y)/2, 0.f, b.r, b.g, b.b, b.a));
//			}
//			tmp->push_back(b);
//		}
//		filaments[i] = tmp;
//		for (Vertex &v: *filaments[i]) {
//			double x = v.x/cellSize + 0.5;
//			double y = v.y/cellSize + 0.5;
//			size_t j0 = (size_t)x;
//			size_t i0 = (size_t)y;
//			size_t i1 = i0 + 1, j1 = j0 + 1;
//			x -= j0;
//			y -= i0;
//			double v00, v01, v10, v11;
//			double vx, vy;
//			v00 = cur_vx[i0][j0];
//			v01 = cur_vx[i1][j0];
//			v10 = cur_vx[i0][j1];
//			v11 = cur_vx[i1][j1];
//			vx = bilinearInterpolate(x, y, v00, v01, v10, v11);
//			v00 = cur_vy[i0][j0];
//			v01 = cur_vy[i1][j0];
//			v10 = cur_vy[i0][j1];
//			v11 = cur_vy[i1][j1];
//			vy = bilinearInterpolate(x, y, v00, v01, v10, v11);
//			v.x += (float)vx * dt * cellSize;
//			v.y += (float)vy * dt * cellSize;
//			v.x = std::max(0.f, v.x);
//			v.x = std::min((float)frameWidth, v.x);
//			v.y = std::max(0.f, v.y);
//			v.y = std::min((float)frameHeight, v.y);
//			if (ages[i] > maxAge / 2)
//			v.a = 2 - (float)ages[i] / (maxAge / 2);
//		}
//	}
//}
