/* ----------------------------------------------------------------------------
  JAM_AXI_CYLIN_RMS

    Wrapper for second moment calculator designed to interface with Python.

    INPUTS
      xp : projected x' [pc]
      yp : projected y' [pc]
      nxy : number of x' and y' values given
      incl : inclination [radians]
      lum_area : projected luminous MGE area
      lum_sigma : projected luminous MGE width
      lum_q : projected luminous MGE flattening
      pot_area : projected potential MGE area
      pot_sigma : projected potential MGE sigma
      pot_q : projected potential MGE flattening
      beta : velocity anisotropy (1-vz^2/vr^2)
      nrad : number of radial bins in interpolation grid
      nang : number of angular bins in interpolation grid
      rxx : array to hold the xx second moments calculated
      ryy : array to hold the yy second moments calculated
      rzz : array to hold the zz second moments calculated
      rxy : array to hold the xy second moments calculated
      rxz : array to hold the xz second moments calculated
      ryz : array to hold the yz second moments calculated
---------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "jam.h"
#include "../mge/mge.h"


void jam_axi_pot(double *r, double *z, int nrz, double incl, \
double *pot_area, double *pot_sigma, double *pot_q, int pot_total, \
int nrad, int nang, double *potval) {

    struct multigaussexp pot;
    double* mu;
    int i;

    // put potential MGE components into structure
    pot.area = pot_area;
    pot.sigma = pot_sigma;
    pot.q = pot_q;
    pot.ntotal = pot_total;

    // calculate xx moments and put into results array
    mu = jam_axi_pot_mmt(r, z, nrz, incl, &pot, nrad, nang);
    for (i=0; i<nrz; i++) {
        potval[i] = mu[i];
    }

    // free memory
    free( mu );

    return;
}
