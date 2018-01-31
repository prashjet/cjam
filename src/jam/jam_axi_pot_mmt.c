/* ----------------------------------------------------------------------------
  JAM_AXI_CYLIN_RMS_MMT

    Calculates second moment.

    INPUTS
      xp    : projected x' [pc]
      yp    : projected y' [pc]
      nxy   : number of x' and y' values given
      incl  : inclination [radians]
      lum   : projected luminous MGE
      pot   : projected potential MGE
      beta  : velocity anisotropy (1 - vz^2 / vr^2)
      nrad  : number of radial bins in interpolation grid
      nang  : number of angular bins in interpolation grid
      vv    : velocity integral selector (1=xx, 2=yy, 3=zz, 4=xy, 5=xz, 6=yz)

    NOTES
      * Based on janis2_second_moment IDL code by Michele Cappellari.
      * This version does not implement PDF convolution.

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
#include "jam.h"
#include "../mge/mge.h"
#include "../tools/tools.h"
#include "../interp/interp.h"


double* jam_axi_pot_mmt( double *r, double *z, int nrz, double incl, \
        struct multigaussexp *pot, \
        int nrad, int nang) {

    int i, j, k, l, npol;
    double qmed, *rell, *rinp, *einp, step, rmax, *lograd, *rad, *ang, *angvec;
    double *wm2, *mu, *mui, *rpol, *zpol, **mupol;
    struct multigaussexp ipot;

    // convert from projected MGEs to intrinsic MGEs
    ipot = mge_deproject( pot, incl );

    // array to store second moment
    mu = (double *) malloc( nrz * sizeof( double ) );

    // skip the interpolation when computing just a few points
    if ( nrad * nang > nrz ) {

        // potential
        wm2 = jam_axi_pot_wmmt( r, z, nrz, incl, pot );

        for ( i = 0; i < nrz; i++ ) {
            mu[i] = wm2[i];
        }

        free( wm2 );

        return mu;
    }


    // ---------------------------------


    // elliptical radius of input (r,z)
    qmed = mge_qmed( &ipot, maximum( r, nrz ) );
    rell = (double *) malloc( nrz * sizeof( double ) );
    for ( i = 0; i < nrz; i++ ) \
        rell[i] = sqrt( r[i] * r[i] + z[i] * z[i] / qmed / qmed );

    // set interpolation grid parameters
    step = minimum( rell, nrz );
    if ( step <= 0.001 ) step = 0.001;          // minimum radius of 0.001 pc
    npol = nrad * nang;

    // make linear grid in log of elliptical radius
    rmax = maximum( rell, nrz );
    lograd = range( log( step * 0.99 ), log( rmax * 1.01 ), nrad, False );
    rad = (double *) malloc( nrad * sizeof( double ) );
    for ( i = 0; i < nrad; i++ ) rad[i] = exp( lograd[i] );

    // make linear grid in eccentric anomaly
    ang = range( -M_PI, -M_PI / 2., nang, False );
    angvec = range( -M_PI, M_PI, 4 * nang - 3, False );

    // convert grid to cartesians
    rpol = (double *) malloc( npol * sizeof( double ) );
    zpol = (double *) malloc( npol * sizeof( double ) );
    for ( i = 0; i < nrad; i++ ) {
        for ( j = 0; j < nang; j++ ) {
            rpol[i*nang+j] = rad[i] * cos( ang[j] );
            zpol[i*nang+j] = rad[i] * sin( ang[j] ) * qmed;
        }
    }

    // set up interpolation grid arrays
      mupol = (double **) malloc( nrad * sizeof( double * ) );
      for ( j = 0; j < nrad; j++ ) {
          mupol[j] = (double *) malloc( ( 4 * nang - 3 ) * sizeof( double ) );
      }

    // ---------------------------------


    // weighted second moment on polar grid
    wm2 = jam_axi_pot_wmmt( rpol, zpol, npol, incl, pot );

    // model velocity on the polar grid
    for ( j = 0; j < nrad; j++ ) {
        for ( k = 0; k < nang; k++ ) {

            l = j * nang + k;
            mupol[j][k] = wm2[l];
            mupol[j][2*nang-2-k] = mupol[j][k];
            mupol[j][2*nang-2+k] = mupol[j][k];
            mupol[j][4*nang-4-k] = mupol[j][k];

        }
    }

    // ---------------------------------


    // elliptical radius and eccentric anomaly of inputs
    rinp = (double *) malloc( nrz * sizeof( double ) );
    einp = (double *) malloc( nrz * sizeof( double ) );
    for ( i = 0; i < nrz; i++ ) {
        rinp[i] = sqrt( pow( r[i], 2. ) + pow( z[i] / qmed, 2. ) );
        einp[i] = atan2( z[i] / qmed, r[i] );
    }

    mui = interp2dpol( mupol, rad, angvec, rinp, einp, nrad, 4*nang-3, nrz );
    for ( i = 0; i < nrz; i++ ) \
        mu[i] = mui[i];
    free( mui );

    free( rell );
    free( rad );
    free( ang );
    free( angvec );
    free( rpol );
    free( zpol );
    for ( j = 0; j < nrad; j++ ) free( mupol[j] );
    free( mupol );
    free( wm2 );
    free( rinp );
    free( einp );

    return mu;

}
