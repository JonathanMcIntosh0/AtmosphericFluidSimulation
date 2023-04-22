#ifndef GRID_SOLVER_H
#define GRID_SOLVER_H

#define BIG_GRID

#include "GLHeaders.h"

#include <vector>
#include <deque>
#include <memory>

#include "Vertex.h"
#include "Array3D.h"

const int res = 5;
#ifndef BIG_GRID
const int cols = 25, rows = 20;
const int cellSize = 40; // cell size in display
#else
const int lats = 18 * res, longs = 36 * res, layers = 4 * res;
const int cellSize = 40 / res; // cell size in display
#endif
const int frameWidth = longs * cellSize, frameHeight = layers * cellSize, frameDepth = lats * cellSize;

template<class T>
using GridArray = Array3D<T, longs + 2, lats + 2, layers + 2>;

extern float nv, kappa, buoyancy, dt, velScale;
extern float simTime;
extern float dx;
//extern float refinementThreshold;
extern int iterations;

//extern double tpAvg[lats];
extern GridArray<double> tpSrc;
//extern double vxSrc[rows + 2][cols + 2][layers + 2];
//extern double vySrc[rows + 2][cols + 2][layers + 2];

extern GridArray<double> *cur_tp;
//extern double (*cur_tp)[lats + 2][layers + 2];

//extern Vertex gridVertices[longs + 1][lats + 1][layers + 1];
extern Array3D<Vertex, longs + 1, lats + 1, layers + 1> gridVertices;
extern std::vector<std::vector<GLuint>> gridIndices_long;
extern std::vector<std::vector<GLuint>> gridIndices_lat;
extern std::vector<std::vector<GLuint>> gridIndices_vert;

//extern Vertex velVertices[longs][lats][layers][2];
extern Array3D<Vertex[2], longs, lats, layers> velVertices;
extern std::vector<std::vector<GLuint>> velIndices_long;
extern std::vector<std::vector<GLuint>> velIndices_lat;
extern std::vector<std::vector<GLuint>> velIndices_vert;

void initGrid();
double getReferenceTemperature(GridArray<double> tp);
void step();
void updateGrid();
//void addFilament(double x, double y);
//void updateFilament();

#endif
