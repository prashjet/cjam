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


void jam_axi_cylin_rms(double *xp, double *yp, int nxy, double incl, \
double *lum_area, double *lum_sigma, double *lum_q, int lum_total, \
double *pot_area, double *pot_sigma, double *pot_q, int pot_total, \
double *beta, int nrad, int nang, double *rrr, double *rff, double *rzz) {

    struct multigaussexp lum, pot;
    double** mu;
    int i, check;

    // put luminous MGE components into structure
    lum.area = lum_area;
    lum.sigma = lum_sigma;
    lum.q = lum_q;
    lum.ntotal = lum_total;

    // put potential MGE components into structure
    pot.area = pot_area;
    pot.sigma = pot_sigma;
    pot.q = pot_q;
    pot.ntotal = pot_total;

    // check for any non-zero beta or non-unity flattening
    check = 0;
    for (i=0; i<lum.ntotal; i++) if (beta[i]!=0.) check++;
    for (i=0; i<lum.ntotal; i++) if (lum_q[i]!=1.) check++;
    for (i=0; i<pot.ntotal; i++) if (pot_q[i]!=1.) check++;

    // calculate xx moments and put into results array
    mu = jam_axi_cylin_rms_mmt(xp, yp, nxy, incl, &lum, &pot, beta, nrad, nang, check);
    for (i=0; i<nxy; i++) {
        rrr[i] = mu[0][i];
        rff[i] = mu[1][i];
        rzz[i] = mu[2][i];
    }

    // free memory
    for ( i = 0; i < 3; i++ ) \
        free( mu[i] );
    free( mu );

    return;
}
