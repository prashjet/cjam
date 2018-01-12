/* ----------------------------------------------------------------------------
  JAM_AXI_CYLIN_RMS_MGEINT

    Integrand for the MGE integral required for second moment calculation.

    INPUTS
      u      : integration variable
      params : function parameters passed as a structure

    NOTES
      * Based on janis2_jeans_mge_integrand IDL code by Michele Cappellari.

  Mark den Brok
  Laura L Watkins [lauralwatkins@gmail.com]

  This code is released under a BSD 2-clause license.
  If you use this code for your research, please cite:
  Watkins et al. 2013, MNRAS, 436, 2598
  "Discrete dynamical models of omega Centauri"
  http://adsabs.harvard.edu/abs/2013MNRAS.436.2598W
---------------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../mge/mge.h"
#include "jam.h"


double jam_axi_cylin_rms_mgeint( double u, void *params ) {

    struct params_potint *p;
    double e2u2p, f, expnt, sum;
    int j, k;

    double u2 = u * u;

    p = params;

    sum = 0.;
    for ( j = 0; j < p->pot->ntotal; j++ ) { //mass gaussians

        e2u2p = u2 * p->e2p[j];
        expnt = - 0.5 * u2 / p->s2p[j] * ( p->r2 + p->z2 / ( 1. - e2u2p ) );
        sum += p->pot->q[j] * p->pot->area[j] * pow( 2. * M_PI , 1.5 ) \
                * p->s2p[j] / sqrt( 1. - e2u2p ) \
                * exp( expnt );

    }

    return - sqrt(2./M_PI) * G * sum;

}
