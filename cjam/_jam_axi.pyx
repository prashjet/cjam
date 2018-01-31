# Python wrappers for JAM functions in C.

import numpy as np
from astropy import table, units as u
cimport cython_jam


def axi_vel(xp, yp, incl, lum_area, lum_sigma, lum_q, pot_area, pot_sigma, pot_q,
    beta, kappa, nrad=30, nang=7):

    # get array lengths needed for C
    nxy = len(xp)
    lum_total = len(lum_area)
    pot_total = len(pot_area)

    # set array types for C
    cdef double [:] c_xp
    cdef double [:] c_yp
    cdef double [:] c_lum_area
    cdef double [:] c_lum_sigma
    cdef double [:] c_lum_q
    cdef double [:] c_pot_area
    cdef double [:] c_pot_sigma
    cdef double [:] c_pot_q
    cdef double [:] c_beta
    cdef double [:] c_kappa
    cdef double [:] c_vx
    cdef double [:] c_vy
    cdef double [:] c_vz

    # set C arrays to be views into the input arrays
    c_xp = np.array(xp, dtype=np.double, copy=False)
    c_yp = np.array(yp, dtype=np.double, copy=False)
    c_lum_area = np.array(lum_area, dtype=np.double, copy=False)
    c_lum_sigma = np.array(lum_sigma, dtype=np.double, copy=False)
    c_lum_q = np.array(lum_q, dtype=np.double, copy=False)
    c_pot_area = np.array(pot_area, dtype=np.double, copy=False)
    c_pot_sigma = np.array(pot_sigma, dtype=np.double, copy=False)
    c_pot_q = np.array(pot_q, dtype=np.double, copy=False)
    c_beta = np.array(beta, dtype=np.double, copy=False)
    c_kappa = np.array(kappa, dtype=np.double, copy=False)

    # create empty arrays to store the result and create C views into them
    vx = np.zeros(nxy)
    c_vx = vx
    vy = np.zeros(nxy)
    c_vy = vy
    vz = np.zeros(nxy)
    c_vz = vz

    # now call the JAM code
    cython_jam.jam_axi_vel(&c_xp[0], &c_yp[0], nxy, incl,
        &c_lum_area[0], &c_lum_sigma[0], &c_lum_q[0], lum_total,
        &c_pot_area[0], &c_pot_sigma[0], &c_pot_q[0], pot_total,
        &c_beta[0], &c_kappa[0], nrad, nang, &c_vx[0], &c_vy[0], &c_vz[0])

    return vx, vy, vz



def axi_rms(xp, yp, incl, lum_area, lum_sigma, lum_q, pot_area, pot_sigma, pot_q,
    beta, nrad=30, nang=7):

    # get array lengths needed for C
    nxy = len(xp)
    lum_total = len(lum_area)
    pot_total = len(pot_area)

    # set array types for C
    cdef double [:] c_xp
    cdef double [:] c_yp
    cdef double [:] c_lum_area
    cdef double [:] c_lum_sigma
    cdef double [:] c_lum_q
    cdef double [:] c_pot_area
    cdef double [:] c_pot_sigma
    cdef double [:] c_pot_q
    cdef double [:] c_beta
    cdef double [:] c_rxx
    cdef double [:] c_ryy
    cdef double [:] c_rzz
    cdef double [:] c_rxy
    cdef double [:] c_rxz
    cdef double [:] c_ryz

    # set C arrays to be views into the input arrays
    c_xp = np.array(xp, dtype=np.double, copy=False)
    c_yp = np.array(yp, dtype=np.double, copy=False)
    c_lum_area = np.array(lum_area, dtype=np.double, copy=False)
    c_lum_sigma = np.array(lum_sigma, dtype=np.double, copy=False)
    c_lum_q = np.array(lum_q, dtype=np.double, copy=False)
    c_pot_area = np.array(pot_area, dtype=np.double, copy=False)
    c_pot_sigma = np.array(pot_sigma, dtype=np.double, copy=False)
    c_pot_q = np.array(pot_q, dtype=np.double, copy=False)
    c_beta = np.array(beta, dtype=np.double, copy=False)

    # create empty arrays to store the results and create C views into them
    rxx = np.empty(nxy)
    c_rxx = rxx
    ryy = np.empty(nxy)
    c_ryy = ryy
    rzz = np.empty(nxy)
    c_rzz = rzz
    rxy = np.empty(nxy)
    c_rxy = rxy
    rxz = np.empty(nxy)
    c_rxz = rxz
    ryz = np.empty(nxy)
    c_ryz = ryz

    # now call the JAM code
    cython_jam.jam_axi_rms(&c_xp[0], &c_yp[0], nxy, incl,
        &c_lum_area[0], &c_lum_sigma[0], &c_lum_q[0], lum_total,
        &c_pot_area[0], &c_pot_sigma[0], &c_pot_q[0], pot_total,
        &c_beta[0], nrad, nang, &c_rxx[0], &c_ryy[0], &c_rzz[0], &c_rxy[0],
        &c_rxz[0], &c_ryz[0])

    return rxx, ryy, rzz, rxy, rxz, ryz



def axisymmetric(xp, yp, tracer_mge, potential_mge, distance, beta=0, kappa=0,
    nscale=1, mscale=1, incl=np.pi/2.*u.rad, mbh=0*u.Msun, rbh=0*u.arcsec,
    nrad=30, nang=7):

    # make sure anisotropy and rotation arrays are the correct length
    beta = np.ones(len(tracer_mge))*beta
    kappa = np.ones(len(tracer_mge))*kappa

    # copy MGEs so that changes we make here aren't propagated
    tracer_copy = tracer_mge.copy()
    potential_copy = potential_mge.copy()

    # adjust tracer MGE by Nscale
    tracer_copy["i"] *= nscale

    # adjust potential MGE by M/L
    potential_copy["i_mass"] = potential_copy["i"] * mscale

    # add BH to potential gaussian
    if mbh>0 and rbh>0:
        potential_copy.add_row([0, mbh/2/np.pi/(rbh*distance/u.rad).to("pc")**2,
            rbh, 1])
        potential_copy.sort("n")

    # calculate first moments
    vx, vy, vz = axi_vel(\
        (xp*distance/u.rad).to("pc").value,
        (yp*distance/u.rad).to("pc").value,
        incl.to("rad").value,
        (tracer_copy["i"]).to("Lsun/pc**2").value,
        (tracer_copy["s"]*distance/u.rad).to("pc").value,
        tracer_copy["q"],
        potential_copy["i_mass"].to("Msun/pc**2").value,
        (potential_copy["s"]*distance/u.rad).to("pc").value,
        potential_copy["q"],
        beta,
        kappa)

    # calculate second moments
    rxx, ryy, rzz, rxy, rxz, ryz = axi_rms(\
        (xp*distance/u.rad).to("pc").value,
        (yp*distance/u.rad).to("pc").value,
        incl.to("rad").value,
        (tracer_copy["i"]).to("Lsun/pc**2").value,
        (tracer_copy["s"]*distance/u.rad).to("pc").value,
        tracer_copy["q"],
        potential_copy["i"].to("Msun/pc**2").value,
        (potential_copy["s"]*distance/u.rad).to("pc").value,
        potential_copy["q"],
        beta)

    # put results into astropy table, also convert PMs to mas/yr
    kms2masyr = (u.km/u.s*u.rad/distance).to("mas/yr")
    moments = table.QTable()
    moments["vx"] = vx*kms2masyr
    moments["vy"] = vy*kms2masyr
    moments["vz"] = vz*u.km/u.s
    moments["v2xx"] = rxx*kms2masyr**2
    moments["v2yy"] = ryy*kms2masyr**2
    moments["v2zz"] = rzz*(u.km/u.s)**2
    moments["v2xy"] = rxy*kms2masyr**2
    moments["v2xz"] = rxz*kms2masyr*u.km/u.s
    moments["v2yz"] = ryz*kms2masyr*u.km/u.s

    return moments



def axi_cylin_rms(r, z, incl, lum_area, lum_sigma, lum_q, pot_area, pot_sigma, pot_q,
    beta, nrad, nang):

    # get array lengths needed for C
    nrz = len(r)
    lum_total = len(lum_area)
    pot_total = len(pot_area)

    # set array types for C
    cdef double [:] c_r
    cdef double [:] c_z
    cdef double [:] c_lum_area
    cdef double [:] c_lum_sigma
    cdef double [:] c_lum_q
    cdef double [:] c_pot_area
    cdef double [:] c_pot_sigma
    cdef double [:] c_pot_q
    cdef double [:] c_beta
    cdef double [:] c_rrr
    cdef double [:] c_rff
    cdef double [:] c_rzz

    # set C arrays to be views into the input arrays
    c_r = np.array(r, dtype=np.double, copy=False)
    c_z = np.array(z, dtype=np.double, copy=False)
    c_lum_area = np.array(lum_area, dtype=np.double, copy=False)
    c_lum_sigma = np.array(lum_sigma, dtype=np.double, copy=False)
    c_lum_q = np.array(lum_q, dtype=np.double, copy=False)
    c_pot_area = np.array(pot_area, dtype=np.double, copy=False)
    c_pot_sigma = np.array(pot_sigma, dtype=np.double, copy=False)
    c_pot_q = np.array(pot_q, dtype=np.double, copy=False)
    c_beta = np.array(beta, dtype=np.double, copy=False)

    # create empty arrays to store the results and create C views into them
    rrr = np.empty(nrz)
    c_rrr = rrr
    rff = np.empty(nrz)
    c_rff = rff
    rzz = np.empty(nrz)
    c_rzz = rzz

    # now call the JAM code
    cython_jam.jam_axi_cylin_rms(&c_r[0], &c_z[0], nrz, incl,
        &c_lum_area[0], &c_lum_sigma[0], &c_lum_q[0], lum_total,
        &c_pot_area[0], &c_pot_sigma[0], &c_pot_q[0], pot_total,
        &c_beta[0], nrad, nang, &c_rrr[0], &c_rff[0], &c_rzz[0])

    return rrr, rff, rzz



def axisymmetric_cylin(r, z, tracer_mge, potential_mge, distance, beta=0, kappa=0,
    nscale=1, mscale=u.Msun/u.Lsun, incl=np.pi/2.*u.rad, mbh=0*u.Msun, rbh=0*u.arcsec,
    nrad=30, nang=7):

    # make sure anisotropy and rotation arrays are the correct length
    beta = np.ones(len(tracer_mge))*beta
    kappa = np.ones(len(tracer_mge))*kappa

    # copy MGEs so that changes we make here aren't propagated
    tracer_copy = tracer_mge.copy()
    potential_copy = potential_mge.copy()

    # adjust tracer MGE by Nscale
    tracer_copy["i"] *= nscale

    # adjust potential MGE by M/L
    potential_copy["i_mass"] = potential_copy["i"] * mscale

    # add BH to potential gaussian
    if mbh>0 and rbh>0:
        potential_copy.add_row([0, mbh/2/np.pi/(rbh*distance/u.rad).to("pc")**2,
            rbh, 1])
        potential_copy.sort("n")

    # calculate second moments
    rrr, rff, rzz = axi_cylin_rms(\
        r.to("pc").value,
        z.to("pc").value,
        incl.to("rad").value,
        (tracer_copy["i"]).to("Lsun/pc**2").value,
        (tracer_copy["s"]*distance/u.rad).to("pc").value,
        tracer_copy["q"],
        potential_copy["i_mass"].to("Msun/pc**2").value,
        (potential_copy["s"]*distance/u.rad).to("pc").value,
        potential_copy["q"],
        beta,
        nrad,
        nang
        )

    # put results into astropy table
    moments = table.QTable()
    moments["v2rr"] = rrr*(u.km/u.s)**2
    moments["v2ff"] = rff*(u.km/u.s)**2
    moments["v2zz"] = rzz*(u.km/u.s)**2

    return moments

def axi_pot(r, z, incl, pot_area, pot_sigma, pot_q, nrad, nang):

    # get array lengths needed for C
    nrz = len(r)
    pot_total = len(pot_area)

    # set array types for C
    cdef double [:] c_r
    cdef double [:] c_z
    cdef double [:] c_pot_area
    cdef double [:] c_pot_sigma
    cdef double [:] c_pot_q
    cdef double [:] c_pot

    # set C arrays to be views into the input arrays
    c_r = np.array(r, dtype=np.double, copy=False)
    c_z = np.array(z, dtype=np.double, copy=False)
    c_pot_area = np.array(pot_area, dtype=np.double, copy=False)
    c_pot_sigma = np.array(pot_sigma, dtype=np.double, copy=False)
    c_pot_q = np.array(pot_q, dtype=np.double, copy=False)

    # create empty arrays to store the results and create C views into them
    pot = np.empty(nrz)
    c_pot = pot

    # now call the JAM code
    cython_jam.jam_axi_pot(&c_r[0], &c_z[0], nrz, incl,
        &c_pot_area[0], &c_pot_sigma[0], &c_pot_q[0], pot_total,
        nrad, nang, &c_pot[0])

    return pot



def axisymmetric_pot(r, z, potential_mge, distance,
    mscale=u.Msun/u.Lsun, incl=np.pi/2.*u.rad, mbh=0*u.Msun, rbh=0*u.arcsec,
    nrad=30, nang=7):

    # copy MGEs so that changes we make here aren't propagated
    potential_copy = potential_mge.copy()

    # adjust potential MGE by M/L
    potential_copy["i_mass"] = potential_copy["i"] * mscale

    # add BH to potential gaussian
    if mbh>0 and rbh>0:
        potential_copy.add_row([0, mbh/2/np.pi/(rbh*distance/u.rad).to("pc")**2,
            rbh, 1])
        potential_copy.sort("n")

    # calculate potential
    pot = axi_pot(\
        r.to("pc").value,
        z.to("pc").value,
        incl.to("rad").value,
        potential_copy["i_mass"].to("Msun/pc**2").value,
        (potential_copy["s"]*distance/u.rad).to("pc").value,
        potential_copy["q"],
        nrad,
        nang
        )

    return pot*(u.km/u.s)**2
