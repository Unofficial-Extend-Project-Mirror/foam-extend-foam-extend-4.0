/*
 * Copyright 2001, Regents of the University of Minnesota
 *
 * mgridgen.c
 *
 * This file contains the top level routines for the sparse hierarchical
 * clustering algorithm.
 *
 * George Irene
 */

#include "mgridgen.h"
#include "rand48.h"


/*************************************************************************
* This function is the entry point for the SHCluster() routine.
**************************************************************************/
void MGridGen(int nvtxs, idxtype *xadj, realtype *vvol, realtype *vsurf,
              idxtype *adjncy, realtype *adjwgt, int minsize, int maxsize,
              int *options, int *nmoves, int *nparts, idxtype *part)
{
  GraphType graph;
  CtrlType ctrl;

  srand(4321);
  srand48(7654321L);

  /*------------------------------------------------------------
   * Set up the various control structures
   *------------------------------------------------------------*/
  ctrl.CType = options[OPTION_CTYPE];
  ctrl.RType = options[OPTION_RTYPE];
  ctrl.dbglvl = options[OPTION_DBGLVL];
  ctrl.dim = options[OPTION_DIM];
  ctrl.minsize = minsize;
  ctrl.maxsize = maxsize;
  ctrl.nparts = -1;

  /*------------------------------------------------------------
   * Set up the graph
   *------------------------------------------------------------*/
  SetUpGraph(&graph, nvtxs, xadj, vvol, vsurf, adjncy, adjwgt);

  CreateGrid(&ctrl, &graph);

  *nparts = ctrl.nparts;
  icopy(nvtxs, graph.where, part);
  *nmoves = graph.nmoves;

  FreeGraph(&graph);
}


/*************************************************************************
* This function is the entry point for performing refinement
**************************************************************************/
void MGridGenRefine(int nvtxs, idxtype *xadj, realtype *vvol, realtype *vsurf,
                    idxtype *adjncy, idxtype *fusedinfo, realtype *adjwgt,
                    int minsize, int maxsize, int *options, int *nmoves,
                    int *nparts, idxtype *part)
{
  int i;
  GraphType graph;
  CtrlType ctrl;

  srand(4321);
  srand48(7654321L);

  /*------------------------------------------------------------
   * Set up the various control structures
   *------------------------------------------------------------*/
  ctrl.CType = options[OPTION_CTYPE];
  ctrl.RType = options[OPTION_RTYPE];
  ctrl.dbglvl = options[OPTION_DBGLVL];
  ctrl.dim = options[OPTION_DIM];
  ctrl.minsize = minsize;
  ctrl.maxsize = maxsize;
  ctrl.nparts = -1;

  /*------------------------------------------------------------
   * Set up the graph
   *------------------------------------------------------------*/
  SetUpGraph(&graph, nvtxs, xadj, vvol, vsurf, adjncy, adjwgt);
  graph.cmap = NULL;

  graph.where = idxmalloc(graph.nvtxs, "graph.where");
  for (i=0; i<graph.nvtxs; i++)
    graph.where[i] = fusedinfo[i];

  RefineKWayOnce(&ctrl, &graph, 10);

  *nparts = ctrl.nparts;
  icopy(nvtxs, graph.where, part);
  *nmoves = graph.nmoves;

  FreeGraph(&graph);
}


/*************************************************************************
* This function creates the coarse grid
**************************************************************************/
void CreateGrid(CtrlType *ctrl, GraphType *graph)
{
  GraphType *cgraph;

  cgraph = Coarsen(ctrl, graph);

  RefineKWay(ctrl, graph, cgraph, 10);
}
