#include "integrator.h"
#include <iostream>
#include <cmath>    // to include sqrt(), etc.
#include <cstdlib>  // for atoi() and atof()
#include <unistd.h> // for getopt()

// http://arxiv.org/pdf/0708.0738v4.pdf
// Time-stamp: <2002-01-18 21:51:36 piet>
//================================================================
//                                                                |
//           /__----__                         ........            |
//       .           \                     ....:        :.          |
//       :                 _\|/_         ..:                         |
//       :                   /|\         :                     _\|/_  |
//  ___   ___                  _____                      ___    /|\   |
//  /     |   \    /\ \     / |   |   \  / |        /\    |   \         |
//  |   __ |___/   /  \ \   /  |   |    \/  |       /  \   |___/         |
//   |    | |  \   /____\ \ /   |   |    /   |      /____\  |   \     \/  |
//     \___| |   \ /      \ V    |   |   /    |____ /      \ |___/     |   |
//                                                                      /   |
//              :                       _/|     :..                    |/    |
//                :..               ____/           :....          ..         |
/*   o   //          :.    _\|/_     /                   :........:           |
*  O  `//\                 /|\                                               |
*  |     /\                                                                  |
*=============================================================================
*
*  nbody_sh1.C:  an N-body integrator with a shared but variable time step
*                (the same for all particles but changing in time), using
*                the Hermite integration scheme.
*
*                ref.: Hut, P., Makino, J. & McMillan, S., 1995,
*                      Astrophysical Journal Letters 443, L93-L96.
*
*  note: in this first version, all functions are included in one file,
*        without any use of a special library or header files.
*_____________________________________________________________________________
*
*  usage: nbody_sh1 [-h (for help)] [-d step_size_control_parameter]
*                   [-e diagnostics_interval] [-o output_interval]
*                   [-t total_duration] [-i (start output at t = 0)]
*                   [-x (extra debugging diagnostics)]
*
*         "step_size_control_parameter" is a coefficient determining the
*            the size of the shared but variable time step for all particles
*
*         "diagnostics_interval" is the time between output of diagnostics,
*            in the form of kinetic, potential, and total energy; with the
*            -x option, a dump of the internal particle data is made as well
*
*         "output_interval" is the time between successive snapshot outputs
*
*         "total_duration" is the integration time, until the program stops
*
*         Input and output are written from the standard i/o streams.  Since
*         all options have sensible defaults, the simplest way to run the code
*         is by only specifying the i/o files for the N-body snapshots:
*
*            nbody_sh1 < data.in > data.out
*
*         The diagnostics information will then appear on the screen.
*         To capture the diagnostics information in a file, capture the
*         standard error stream as follows:
*
*            (nbody_sh1 < data.in > data.out) >& data.err
*
*  Note: if any of the times specified in the -e, -o, or -t options are not an
*        an integer multiple of "step", output will occur slightly later than
*        predicted, after a full time step has been taken.  And even if they
*        are integer multiples, round-off error may induce one extra step.
*_____________________________________________________________________________
*
*  External data format:
*
*     The program expects input of a single snapshot of an N-body system,
*     in the following format: the number of particles in the snapshot n;
*     the time t; mass mi, position ri and velocity vi for each particle i,
*     with position and velocity given through their three Cartesian
*     coordinates, divided over separate lines as follows:
*
*                      n
*                      t
*                      m1 r1_x r1_y r1_z v1_x v1_y v1_z
*                      m2 r2_x r2_y r2_z v2_x v2_y v2_z
*                      ...
*                      mn rn_x rn_y rn_z vn_x vn_y vn_z
*
*     Output of each snapshot is written according to the same format.
*
*  Internal data format:
*
*     The data for an N-body system is stored internally as a 1-dimensional
*     array for the masses, and 2-dimensional arrays for the positions,
*     velocities, accelerations and jerks of all particles.
*_____________________________________________________________________________
*
*    version 1:  Jan 2002   Piet Hut, Jun Makino
*_____________________________________________________________________________
*/

/*-----------------------------------------------------------------------------
*  read_options  --  reads the command line options, and implements them.
*
*  note: when the help option -h is invoked, the return value is set to false,
*        to prevent further execution of the main program; similarly, if an
*        unknown option is used, the return value is set to false.
*-----------------------------------------------------------------------------
*/
Integrator::Integrator(unsigned nSat, double centralMass)
{
    int c = (centralMass > 0.0);
    n = nSat + c;

    pos = new double[n][NDIM];
    vel = new double[n][NDIM];
    acc = new double[n][NDIM];
    jerk = new double[n][NDIM];
    mass = new double[n];
    radius = new double[n];

    dt_param = 1/32768.0;     // control parameter to determine time step size
    nStepsOut = 1000;          // time interval between output of snapshots
    dt_tot = 1e300;         // duration of the integration
    dt_t = 0;
    nSteps = 0;

    for(unsigned i=0; i<n; i++) mass[i] = 1;
    if (c) mass[0] = centralMass;

    double mass_sum = 0;
    for(unsigned i=0; i<n; i++) {
        radius[i] = 10*pow(mass[i], 1.0/9.0);
        mass_sum += mass[i];
    }

    constexpr double WIDTH = 100.0;
    double phi = M_PI*2/nSat;
    for(unsigned i=c; i<n; i++) {
        pos[i][0] = cos(phi*i) * WIDTH;
        pos[i][1] = sin(phi*i) * WIDTH;
        pos[i][2] = rand()/(RAND_MAX+1.0) - 0.5;
    }

    double epot, coll;
    get_acc_jerk_pot_coll(epot, coll);

    for(unsigned i=c; i<n; i++) {
        double SPEED = sqrt(epot/-2.0/n*2.0/mass[i]);
        vel[i][0] = sin(phi*i) * SPEED;
        vel[i][1] = -cos(phi*i) * SPEED;
        vel[i][2] = 0;
    }
}

Integrator::~Integrator() {
    delete [] pos;
    delete [] vel;
    delete [] acc;
    delete [] jerk;
    delete mass;
    delete radius;
}

bool Integrator::read_options(int argc, char* argv[], bool& i_flag, bool& x_flag) {
    int c;
    while ((c = getopt(argc, argv, "hd:e:o:t:ix")) != -1)
        switch(c) {
        case 'h': std::cerr << "usage: " << argv[0]
                           << " [-h (for help)]"
                           << " [-d step_size_control_parameter]\n"
                           << "         [-e diagnostics_interval]"
                           << " [-o output_interval]\n"
                           << "         [-t total_duration]"
                           << " [-i (start output at t = 0)]\n"
                           << "         [-x (extra debugging diagnostics)]"
                           << std::endl;
            return false;         // execution should stop after help
        case 'd': dt_param = atof(optarg);
            break;
        case 'i': i_flag = true;
            break;
        case 't': dt_tot = atof(optarg);
            break;
        case 'x': x_flag = true;
            break;
        case '?': std::cerr << "usage: " << argv[0]
                           << " [-h (for help)]"
                           << " [-d step_size_control_parameter]\n"
                           << "         [-e diagnostics_interval]"
                           << " [-o output_interval]\n"
                           << "         [-t total_duration]"
                           << " [-i (start output at t = 0)]\n"
                           << "         [-x (extra debugging diagnostics)]"
                           << std::endl;
            return false;        // execution should stop after error
        }

    return true;                         // ready to continue program execution
}

/*-----------------------------------------------------------------------------
*  get_snapshot  --  reads a single snapshot from the input stream cin.
*
*  note: in this implementation, only the particle data are read in, and it
*        is left to the main program to first read particle number and time
*-----------------------------------------------------------------------------
*/

void Integrator::get_snapshot(/*double mass[], double pos[][NDIM], double vel[][NDIM]*/) {
    for (unsigned i = 0; i < n ; i++) {
        std::cin >> mass[i];                       // mass of particle i
        for (unsigned k = 0; k < NDIM; k++)
            std::cin >> pos[i][k];                 // position of particle i
        for (unsigned k = 0; k < NDIM; k++)
            std::cin >> vel[i][k];                 // velocity of particle i
    } }

/*-----------------------------------------------------------------------------
*  put_snapshot  --  writes a single snapshot on the output stream cout.
*
*  note: unlike get_snapshot(), put_snapshot handles particle number and time
*-----------------------------------------------------------------------------
*/

void Integrator::put_snapshot(/*const double mass[], const double pos[][NDIM],
                              const double vel[][NDIM], */double t) {
    // TODO put data in to temporary buffer
}

/*-----------------------------------------------------------------------------
*  write_diagnostics  --  writes diagnostics on the error stream cerr:
*                         current time; number of integration steps so far;
*                         kinetic, potential, and total energy; absolute and
*                         relative energy errors since the start of the run.
*                         If x_flag (x for eXtra data) is true, all internal
*                         data are dumped for each particle (mass, position,
*                         velocity, acceleration, and jerk).
*
*  note: the kinetic energy is calculated here, while the potential energy is
*        calculated in the function get_acc_jerk_pot_coll().
*-----------------------------------------------------------------------------
*/

void Integrator::write_diagnostics(double& einit, bool init_flag,
                                   bool x_flag) {
    double eKinSum = 0;                       // kinetic energy of the n-body system
    for (unsigned i = 0; i < n ; i++)
        for (unsigned k = 0; k < NDIM ; k++)
            eKinSum += 0.5 * mass[i] * vel[i][k] * vel[i][k];

    ekin = eKinSum;
    etot = ekin + epot;             // total energy of the n-body system

    if (init_flag)                       // at first pass, pass the initial
        einit = etot;                    // energy back to the calling function

//    std::cerr << "at time t = " << t << " , after " << nsteps
//         << " steps :\n  E_kin = " << ekin
//         << " , E_pot = " << epot
//         << " , E_tot = " << etot << std::endl;
//    std::cerr << "                "
//         << "absolute energy error: E_tot - E_init = "
//         << etot - einit << std::endl;
//    std::cerr << "                "
//         << "relative energy error: (E_tot - E_init) / E_init = "
//         << (etot - einit) / einit << std::endl;

    if (x_flag) {
        std::cerr << "  for debugging purposes, here is the internal data "
             << "representation:\n";
        for (unsigned i = 0; i < n ; i++) {
            std::cerr << "    internal data for particle " << i+1 << " : " << std::endl;
            std::cerr << "      ";
            std::cerr << mass[i];
            for (unsigned k = 0; k < NDIM; k++)
                std::cerr << ' ' << pos[i][k];
            for (unsigned k = 0; k < NDIM; k++)
                std::cerr << ' ' << vel[i][k];
            for (unsigned k = 0; k < NDIM; k++)
                std::cerr << ' ' << acc[i][k];
            for (unsigned k = 0; k < NDIM; k++)
                std::cerr << ' ' << jerk[i][k];
            std::cerr << std::endl; } } }

/*-----------------------------------------------------------------------------
*  evolve  --  integrates an N-body system, for a total duration dt_tot.
*              Snapshots are sent to the standard output stream once every
*              time interval dt_out.  Diagnostics are sent to the standard
*              error stream once every time interval dt_dia.
*
*  note: the integration time step, shared by all particles at any given time,
*        is variable.  Before each integration step we use coll_time (short
*        for collision time, an estimate of the time scale for any significant
*        change in configuration to happen), multiplying it by dt_param (the
*        accuracy parameter governing the size of dt in units of coll_time),
*        to obtain the new time step size.
*
*  Before moving any particles, we start with an initial diagnostics output
*  and snapshot output if desired.  In order to write the diagnostics, we
*  first have to calculate the potential energy, with get_acc_jerk_pot_coll().
*  That function also calculates accelerations, jerks, and an estimate for the
*  collision time scale, all of which are needed before we can enter the main
*  integration loop below.
*       In the main loop, we take as many integration time steps as needed to
*  reach the next output time, do the output required, and continue taking
*  integration steps and invoking output this way until the final time is
*  reached, which triggers a `break' to jump out of the infinite loop set up
*  with `while(true)'.
*-----------------------------------------------------------------------------
*/

void Integrator::evolve(double& t, bool init_out, bool x_flag) {
    std::cerr << "Starting a Hermite integration for a " << n
         << "-body system,\n  from time t = " << t
         << " with time step control parameter dt_param = " << dt_param
         << "  until time " << t + dt_tot
         << std::endl;

//    double (* acc)[NDIM] = new double[n][NDIM];          // accelerations and jerks
//    double (* jerk)[NDIM] = new double[n][NDIM];         // for all particles

    int nsteps = 0;               // number of integration time steps completed
    double einit;                   // initial total energy of the system

    double coll_time;                   // collision (close encounter) time scale
    get_acc_jerk_pot_coll(epot, coll_time);
    write_diagnostics(einit, true, x_flag);
    put_snapshot(t);

    qlonglong outSteps = nSteps + nStepsOut;           // next time for snapshot output
    double t_end = t + dt_tot;           // final time, to finish the integration

    while (!isInterruptionRequested()) {
        while (nSteps < outSteps  && t < t_end) {
            double dt = dt_param * coll_time;
            evolve_step(t, dt, epot, coll_time);
            write_diagnostics(einit, true, x_flag);
            nSteps++; }
        if (nSteps >= outSteps) {
            put_snapshot(t);
            outSteps += nSteps; }
    }
}

/*-----------------------------------------------------------------------------
*  evolve_step  --  takes one integration step for an N-body system, using the
*                   Hermite algorithm.
*-----------------------------------------------------------------------------
*/

void Integrator::evolve_step(/*const double mass[], double pos[][NDIM], double vel[][NDIM],
                             double acc[][NDIM], double jerk[][NDIM], */double& t,
                             double dt, double& epot, double& coll_time) {
    double (* old_pos)[NDIM] = new double[n][NDIM];
    double (* old_vel)[NDIM] = new double[n][NDIM];
    double (* old_acc)[NDIM] = new double[n][NDIM];
    double (* old_jerk)[NDIM] = new double[n][NDIM];

    for (unsigned i = 0; i < n ; i++)
        for (unsigned k = 0; k < NDIM ; k++) {
            old_pos[i][k] = pos[i][k];
            old_vel[i][k] = vel[i][k];
            old_acc[i][k] = acc[i][k];
            old_jerk[i][k] = jerk[i][k]; }

    predict_step(dt);
    get_acc_jerk_pot_coll(epot, coll_time);
    correct_step(old_pos, old_vel, old_acc, old_jerk,
                 dt);
    t += dt;

    delete[] old_pos;
    delete[] old_vel;
    delete[] old_acc;
    delete[] old_jerk; }

/*-----------------------------------------------------------------------------
*  predict_step  --  takes the first approximation of one Hermite integration
*                    step, advancing the positions and velocities through a
*                    Taylor series development up to the order of the jerks.
*-----------------------------------------------------------------------------
*/

void Integrator::predict_step(double dt) {
    for (unsigned i = 0; i < n ; i++)
        for (unsigned k = 0; k < NDIM ; k++) {
            pos[i][k] += vel[i][k]*dt + acc[i][k]*dt*dt/2
                         + jerk[i][k]*dt*dt*dt/6;
            vel[i][k] += acc[i][k]*dt + jerk[i][k]*dt*dt/2; } }

/*-----------------------------------------------------------------------------
*  correct_step  --  takes one iteration to improve the new values of position
*                    and velocities, effectively by using a higher-order
*                    Taylor series constructed from the terms up to jerk at
*                    the beginning and the end of the time step.
*-----------------------------------------------------------------------------
*/

void Integrator::correct_step(const double old_pos[][NDIM], const double old_vel[][NDIM],
                              const double old_acc[][NDIM], const double old_jerk[][NDIM],
                              double dt) {
    for (unsigned i = 0; i < n ; i++)
        for (unsigned k = 0; k < NDIM ; k++) {
            vel[i][k] = old_vel[i][k] + (old_acc[i][k] + acc[i][k])*dt/2
                        + (old_jerk[i][k] - jerk[i][k])*dt*dt/12;
            pos[i][k] = old_pos[i][k] + (old_vel[i][k] + vel[i][k])*dt/2
                        + (old_acc[i][k] - acc[i][k])*dt*dt/12; } }

/*-----------------------------------------------------------------------------
*  get_acc_jerk_pot_coll  --  calculates accelerations and jerks, and as side
*                             effects also calculates potential energy and
*                             the time scale coll_time for significant changes
*                             in local configurations to occur.
*                                                  __                     __
*                                                 |          -->  -->       |
*               M                           M     |           r  . v        |
*   -->          j    -->       -->          j    | -->        ji   ji -->  |
*    a   ==  --------  r    ;    j   ==  -------- |  v   - 3 ---------  r   |
*     ji     |-->  |3   ji        ji     |-->  |3 |   ji      |-->  |2   ji |
*            | r   |                     | r   |  |           | r   |       |
*            |  ji |                     |  ji |  |__         |  ji |     __|
*
*  note: it would be cleaner to calculate potential energy and collision time
*        in a separate function.  However, the current function is by far the
*        most time consuming part of the whole program, with a double loop
*        over all particles that is executed every time step.  Splitting off
*        some of the work to another function would significantly increase
*        the total computer time (by an amount close to a factor two).
*
*  We determine the values of all four quantities of interest by walking
*  through the system in a double {i,j} loop.  The first three, acceleration,
*  jerk, and potential energy, are calculated by adding successive terms;
*  the last, the estimate for the collision time, is found by determining the
*  minimum value over all particle pairs and over the two choices of collision
*  time, position/velocity and sqrt(position/acceleration), where position and
*  velocity indicate their relative values between the two particles, while
*  acceleration indicates their pairwise acceleration.  At the start, the
*  first three quantities are set to zero, to prepare for accumulation, while
*  the last one is set to a very large number, to prepare for minimization.
*       The integration loops only over half of the pairs, with j > i, since
*  the contributions to the acceleration and jerk of particle j on particle i
*  is the same as those of particle i on particle j, apart from a minus sign
*  and a different mass factor.
*-----------------------------------------------------------------------------
*/

void Integrator::get_acc_jerk_pot_coll(double& epot,
                                       double& coll_time) {
    for (unsigned i = 0; i < n ; i++)
        for (unsigned k = 0; k < NDIM ; k++)
            acc[i][k] = jerk[i][k] = 0;
    epot = 0;
    const double VERY_LARGE_NUMBER = 1e300;
    double coll_time_q = VERY_LARGE_NUMBER;      // collision time to 4th power
    double coll_est_q;                           // collision time scale estimate
    // to 4th power (quartic)
    for (unsigned i = 0; i < n ; i++) {
        for (unsigned j = i+1; j < n ; j++) {           // rji[] is the vector from
            double rji[NDIM];                        // particle i to particle j
            double vji[NDIM];                        // vji[] = d rji[] / d t
            for (unsigned k = 0; k < NDIM ; k++) {
                rji[k] = pos[j][k] - pos[i][k];
                vji[k] = vel[j][k] - vel[i][k]; }
            double r2 = 0;                           // | rji |^2
            double v2 = 0;                           // | vji |^2
            double rv_r2 = 0;                        // ( rij . vij ) / | rji |^2
            for (unsigned k = 0; k < NDIM ; k++) {
                r2 += rji[k] * rji[k];
                v2 += vji[k] * vji[k];
                rv_r2 += rji[k] * vji[k]; }
            rv_r2 /= r2;
            double r = sqrt(r2);                     // | rji |
            double r3 = r * r2;                      // | rji |^3

// add the {i,j} contribution to the total potential energy for the system:

            epot -= mass[i] * mass[j] / r;

// add the {j (i)} contribution to the {i (j)} values of acceleration and jerk:

            double da[NDIM];                            // main terms in pairwise
            double dj[NDIM];                            // acceleration and jerk
            for (unsigned k = 0; k < NDIM ; k++) {
                da[k] = rji[k] / r3;                           // see equations
                dj[k] = (vji[k] - 3 * rv_r2 * rji[k]) / r3;    // in the header
            }
            for (unsigned k = 0; k < NDIM ; k++) {
                acc[i][k] += mass[j] * da[k];                 // using symmetry
                acc[j][k] -= mass[i] * da[k];                 // find pairwise
                jerk[i][k] += mass[j] * dj[k];                // acceleration
                jerk[j][k] -= mass[i] * dj[k];                // and jerk
            }

// first collision time estimate, based on unaccelerated linear motion:

            coll_est_q = (r2*r2) / (v2*v2);
            if (coll_time_q > coll_est_q)
                coll_time_q = coll_est_q;

// second collision time estimate, based on free fall:

            double da2 = 0;                                  // da2 becomes the
            for (unsigned k = 0; k < NDIM ; k++)                // square of the
                da2 += da[k] * da[k];                      // pair-wise accel-
            double mij = mass[i] + mass[j];                // eration between
            da2 *= mij * mij;                              // particles i and j

            coll_est_q = r2/da2;
            if (coll_time_q > coll_est_q)
                coll_time_q = coll_est_q; } }                                             // from q for quartic back
    coll_time = sqrt(sqrt(coll_time_q));            // to linear collision time
}

void Integrator::run() {
    evolve(dt_t, false, false);
}

//void Integrator::start()
//{
//    running = true;
//    QThread::start();
//}

//void Integrator::stop()
//{
//    running = false;
//}

unsigned Integrator::getN() const
{
    return n;
}

double Integrator::getPos(unsigned i, unsigned j) const
{
    return pos[i][j];
}

double Integrator::getMass(unsigned i) const
{
    return mass[i];
}

double Integrator::getRadius(unsigned i) const
{
    return radius[i];
}

double Integrator::getDt_t() const
{
    return dt_t;
}

double Integrator::getEkin() const
{
    return ekin;
}

double Integrator::getEpot() const
{
    return epot;
}

qlonglong Integrator::getNSteps() const
{
    return nSteps;
}

double Integrator::getEtot() const
{
    return etot;
}

void Integrator::setDt_param(double value)
{
    dt_param = value;
}

