/* ----------------------------------------------------------------------------
  JAM_AXI_CYLIN_RMS_WMMT

    Calculates weighted second moment.

    INPUTS
      xp    : projected x' [pc]
      yp    : projected y' [pc]
      nxy   : number of x' and y' values given
      incl  : inclination [radians]
      lum   : projected luminous MGE
      pot   : projected potential MGE
      beta  : velocity anisotropy (1 - vz^2 / vr^2)
      vv    : velocity integral selector (1=xx, 2=yy, 3=zz, 4=xy, 5=xz, 6=yz)

    NOTES
    * Based on janis2_weighted_second_moment_squared IDL code by Michele
      Cappellari.

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
#include <gsl/gsl_integration.h>
#include "jam.h"
#include "../mge/mge.h"
#include "../tools/tools.h"


double *jam_axi_pot_wmmt( double *r, double *z, int nrz, double incl, \
        struct multigaussexp *pot ) {

    struct params_potint p;
    struct multigaussexp ipot;
    double *s2p, *e2p;
    double result, error, *sb_mu2;
    int i;

    // convert from projected MGEs to intrinsic MGEs
    ipot = mge_deproject( pot, incl );

    // mge component combinations
    s2p = (double *) malloc( pot->ntotal * sizeof( double ) );
    e2p = (double *) malloc( pot->ntotal * sizeof( double ) );

    for ( i = 0; i < ipot.ntotal; i++ ) {
        s2p[i] = pow( ipot.sigma[i], 2. );
        e2p[i] = 1. - pow( ipot.q[i], 2. );
    }

    // parameters for the integrand function
    p.pot = &ipot;
    p.s2p = s2p;
    p.e2p = e2p;

    // perform integration

    gsl_integration_workspace *w = gsl_integration_workspace_alloc( 1000 );
    gsl_function F;
    F.function = &jam_axi_pot_mgeint;

    // set up interpolation grid arrays
    sb_mu2 = (double *) malloc( nrz * sizeof( double ) );

    // calculate vz dispersion
    for ( i = 0; i < nrz; i++ ) {
        p.r2 = r[i] * r[i];
        p.z2 = z[i] * z[i];
        F.params = &p;
        gsl_integration_qag( &F, 0., 1., 0., 1e-5, 1000, 6, w, \
            &result, &error );
        sb_mu2[i] = result;
    }

    gsl_integration_workspace_free( w );

    free( s2p );
    free( e2p );

    return sb_mu2;

}
