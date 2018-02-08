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

    struct params_cylinrmsint *p;
    double e2u2p, c, d, f, expnt, expnt2, sum;
    int j, k;

    double u2 = u * u;

    p = params;

    sum = 0.;
    for ( j = 0; j < p->pot->ntotal; j++ ) { //mass gaussians

        e2u2p = u2 * p->e2p[j];
        expnt = - 0.5 * u2 / p->s2p[j] * ( p->r2 + p->z2 / ( 1. - e2u2p ) );

        for ( k = 0; k < p->lum->ntotal; k++ ) { // luminous gaussians

            c = p->e2p[j] - p->s2q2l[k] / p->s2p[j];

            expnt2 = - 0.5 / p->s2l[k] * ( p->r2 + p->z2 / p->q2l[k] );
            // expnt2 for taking nu(0,0) to local density nu(R,z)

            switch ( p->vv ) {
                case 1: // v2zz
                    f = p->s2q2l[k];
                    break;
                case 2: // v2rr
                    f = p->kani[k] * p->s2q2l[k];
                    break;
                case 3: // v2ff
                    d = 1. - p->kani[k] * p->q2l[k] \
                        - (( 1. - p->kani[k] ) * c + p->e2p[j] * p->kani[k])*u2;
                    f = d * p->r2 + p->kani[k] * p->s2q2l[k];
                    break;
                default:
                    printf( "No integral selected.  Options: 1=v2zz, " );
                    printf( "2=v2rr, 3=v2ff.\n" );
                    f = 1.;
                    break;
            }

            sum += p->lum->area[k] * exp( expnt2 ) \
                * p->pot->q[j] * p->pot->area[j] * u2 \
                / ( 1. - c * u2 ) / sqrt( 1. - e2u2p ) \
                * f * exp( expnt );

        }
    }

    return 4. * M_PI * G * sum;

}
