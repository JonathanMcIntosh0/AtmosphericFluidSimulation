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
const int N_CELLS = ((longs + 2) * (lats + 2) * (layers + 2));

typedef Array3D<double, longs + 2, lats + 2, layers + 2> GridArray;

extern float nv, kappa, buoyancy, dt, velScale;
extern float simTime;
extern float dx;
//extern float refinementThreshold;
extern int iterations;

//extern double tpAvg[lats];
extern GridArray tpSrc;
//extern double vxSrc[rows + 2][cols + 2][layers + 2];
//extern double vySrc[rows + 2][cols + 2][layers + 2];

extern GridArray *cur_tp;
//extern double (*cur_tp)[lats + 2][layers + 2];

//extern std::deque<std::shared_ptr<std::vector<Vertex>>> filaments;
//extern std::deque<double> ages;
//extern float maxAge;

extern Vertex gridVertices[1][1][1];
//extern Vertex gridVertices[longs + 1][lats + 1][layers + 1];
extern std::vector<std::vector<GLuint>> gridIndices_long;
extern std::vector<std::vector<GLuint>> gridIndices_lat;
extern std::vector<std::vector<GLuint>> gridIndices_vert;

extern Vertex velVertices[1][1][1][2];
//extern Vertex velVertices[longs][lats][layers][2];
extern std::vector<std::vector<GLuint>> velIndices_long;
extern std::vector<std::vector<GLuint>> velIndices_lat;
extern std::vector<std::vector<GLuint>> velIndices_vert;

void initGrid();
double getReferenceTemperature(double tp[longs + 2][lats + 2][layers + 2]);
void step();
void updateGrid();
//void addFilament(double x, double y);
//void updateFilament();

#endif
