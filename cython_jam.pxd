# Cython declarations of C functions.

cdef extern from "../src/jam/jam.h":

    void jam_axi_rms(double *xp, double *yp, int nxy, double incl, \
        double *lum_area, double *lum_sigma, double *lum_q, int lum_total, \
        double *pot_area, double *pot_sigma, double *pot_q, int pot_total, \
        double *beta, int nrad, int nang, double *rxx, double *ryy, \
        double *rzz, double *rxy, double *rxz, double *ryz)

    void jam_axi_rms_los(double *xp, double *yp, int nxy, double incl, \
        double *lum_area, double *lum_sigma, double *lum_q, int lum_total, \
        double *pot_area, double *pot_sigma, double *pot_q, int pot_total, \
        double *beta, int nrad, int nang, double *rzz)

    void jam_axi_vel(double *xp, double *yp, int nxy, double incl, \
        double *lum_area, double *lum_sigma, double *lum_q, int lum_total, \
        double *pot_area, double *pot_sigma, double *pot_q, int pot_total, \
        double *beta, double *kappa, int nrad, int nang, double *vx, \
        double *vy, double *vz)

    void jam_axi_cylin_rms(double *r, double *z, int nrz, double incl, \
        double *lum_area, double *lum_sigma, double *lum_q, int lum_total, \
        double *pot_area, double *pot_sigma, double *pot_q, int pot_total, \
        double *beta, int nrad, int nang, double *rrr, double *rff, \
        double *rzz);

    void jam_axi_pot(double *r, double *z, int nrz, double incl, \
        double *pot_area, double *pot_sigma, double *pot_q, int pot_total, \
        int nrad, int nang, double *potval);
