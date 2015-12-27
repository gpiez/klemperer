#ifndef INTEGRATOR_H
#define INTEGRATOR_H

#include "body.h"

#include <QThread>

class Integrator: public QThread {

public:
    Integrator(unsigned n, double ecc);
    ~Integrator();
    void correct_step(const double old_pos[][NDIM], const double old_vel[][NDIM],
                      const double old_acc[][NDIM], const double old_jerk[][NDIM],
                      double dt);
    void evolve(double& t, bool init_out, bool x_flag);
    void evolve_step(double& t, double dt, double& epot, double& coll_time);
    void get_acc_jerk_pot_coll(double& epot,
                               double& coll_time);
    void get_snapshot();
    void predict_step(double dt);
    void put_snapshot(double t);
    bool read_options(int argc, char* argv[], bool& i_flag, bool& x_flag);
    void run() override;
//    void start();
//    void stop();

    unsigned getN() const;

    double getPos(unsigned, unsigned) const;
    double getMass(unsigned) const;
    double getRadius(unsigned) const;

    double getDt_t() const;

    double getEkin() const;

    double getEpot() const;

    qlonglong getNSteps() const;

    double getEtot() const;

    void setDt_param(double value);

private:
    double  dt_param;     // control parameter to determine time step size
    qlonglong nStepsOut;          // time interval between output of snapshots
    qlonglong nSteps;
    double  dt_tot;         // duration of the integration
    double  dt_t;

    double ekin;
    double epot;                        // potential energy of the n-body system
    double etot;
    unsigned n;

    double (*pos)[NDIM];
    double (*vel)[NDIM];
    double (*acc)[NDIM];
    double (*jerk)[NDIM];
    double (*mass);
    double (*radius);
    void write_diagnostics(double &einit, bool init_flag, bool x_flag);
};

#endif // INTEGRATOR_H
