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


double** jam_axi_cylin_rms_mmt( double *xp, double *yp, int nxy, double incl, \
        struct multigaussexp *lum, struct multigaussexp *pot, double *beta, \
        int nrad, int nang, int check ) {

    int i, j, k, l, npol;
    double qmed, *rell, *r, *e, step, rmax, *lograd, *rad, *ang, *angvec;
    double **wm2, *surf, **mu, *mui, *xpol, *ypol, ***mupol;

    // array to store second moment
    mu = (double **) malloc( 3 * sizeof( double * ) );
    for ( i = 0; i < 3; i++ ) \
        mu[i] = (double *) malloc( nxy * sizeof( double ) );

    // skip the interpolation when computing just a few points
    if ( nrad * nang > nxy ) {

        // weighted second moment
        wm2 = jam_axi_cylin_rms_wmmt( xp, yp, nxy, incl, lum, pot, beta, check );

        // surface brightness
        surf = mge_surf( lum, xp, yp, nxy );

        // second moment
        for ( i = 0; i < 3; i++ ) {
            for ( j = 0; j < nxy; j++ ) {
                mu[i][j] = wm2[i][j] / surf[j];
            }
        }

        for ( i = 0; i < 3; i++ ) \
            free( wm2[i] );
        free( wm2 );
        free( surf );

        return mu;
    }


    // ---------------------------------


    // elliptical radius of input (x,y)
    qmed = mge_qmed( lum, maximum( xp, nxy ) );
    rell = (double *) malloc( nxy * sizeof( double ) );
    for ( i = 0; i < nxy; i++ ) \
        rell[i] = sqrt( xp[i] * xp[i] + yp[i] * yp[i] / qmed / qmed );

    // set interpolation grid parameters
    step = minimum( rell, nxy );
    if ( step <= 0.001 ) step = 0.001;          // minimum radius of 0.001 pc
    npol = nrad * nang;

    // make linear grid in log of elliptical radius
    rmax = maximum( rell, nxy );
    lograd = range( log( step * 0.99 ), log( rmax * 1.01 ), nrad, False );
    rad = (double *) malloc( nrad * sizeof( double ) );
    for ( i = 0; i < nrad; i++ ) rad[i] = exp( lograd[i] );

    // make linear grid in eccentric anomaly
    ang = range( -M_PI, -M_PI / 2., nang, False );
    angvec = range( -M_PI, M_PI, 4 * nang - 3, False );

    // convert grid to cartesians
    xpol = (double *) malloc( npol * sizeof( double ) );
    ypol = (double *) malloc( npol * sizeof( double ) );
    for ( i = 0; i < nrad; i++ ) {
        for ( j = 0; j < nang; j++ ) {
            xpol[i*nang+j] = rad[i] * cos( ang[j] );
            ypol[i*nang+j] = rad[i] * sin( ang[j] ) * qmed;
        }
    }

    // set up interpolation grid arrays
    mupol = (double ***) malloc( 3 * sizeof( double ** ) );
    for ( i = 0; i < 3; i++ ) {
        mupol[i] = (double **) malloc( nrad * sizeof( double * ) );
        for ( j = 0; j < nrad; j++ ) {
            mupol[i][j] = (double *) malloc( ( 4 * nang - 3 ) * sizeof( double ) );
        }
    }

    // ---------------------------------


    // weighted second moment on polar grid
    wm2 = jam_axi_cylin_rms_wmmt( xpol, ypol, npol, incl, lum, pot, beta, check );

    // surface brightness on polar grid
    surf = mge_surf( lum, xpol, ypol, npol );

    // model velocity on the polar grid
    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < nrad; j++ ) {
            for ( k = 0; k < nang; k++ ) {

                l = j * nang + k;
                mupol[i][j][k] = wm2[i][l] / surf[l];
                mupol[i][j][2*nang-2-k] = mupol[i][j][k];
                mupol[i][j][2*nang-2+k] = mupol[i][j][k];
                mupol[i][j][4*nang-4-k] = mupol[i][j][k];

            }
        }
    }

    // ---------------------------------


    // elliptical radius and eccentric anomaly of inputs
    r = (double *) malloc( nxy * sizeof( double ) );
    e = (double *) malloc( nxy * sizeof( double ) );
    for ( i = 0; i < nxy; i++ ) {
        r[i] = sqrt( pow( xp[i], 2. ) + pow( yp[i] / qmed, 2. ) );
        e[i] = atan2( yp[i] / qmed, xp[i] );
    }

    mui = interp2dpol( mupol[0], rad, angvec, r, e, nrad, 4*nang-3, nxy );
    for ( i = 0; i < nxy; i++ ) \
        mu[0][i] = mui[i];
    free( mui );
    mui = interp2dpol( mupol[1], rad, angvec, r, e, nrad, 4*nang-3, nxy );
    for ( i = 0; i < nxy; i++ ) \
        mu[1][i] = mui[i];
    free( mui );
    mui = interp2dpol( mupol[2], rad, angvec, r, e, nrad, 4*nang-3, nxy );
    for ( i = 0; i < nxy; i++ ) \
        mu[2][i] = mui[i];
    free( mui );

    free( rell );
    free( rad );
    free( ang );
    free( angvec );
    free( xpol );
    free( ypol );
    for ( i = 0; i < 3; i++ ) {
        for ( j = 0; j < nrad; j++ ) free( mupol[i][j] );
        free( mupol[i] );
    }
    free( mupol );
    for ( i = 0; i < 3; i++ ) \
        free( wm2[i] );
    free( wm2 );
    free( surf );
    free( r );
    free( e );

    return mu;

}
