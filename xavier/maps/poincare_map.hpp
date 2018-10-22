/*************************************************************************
POMATO: POincare MAp TOpology: Extraction and visualization of the
        topological structure of the circular restricted three-body 
        problem for orbital mechanics applications. 

Authors: Wayne Schlei and Xavier Tricoche 

Copyright (c) 2013-2018, Purdue University 

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
**************************************************************************/


#ifndef __POINCARE_MAP_HPP__
#define __POINCARE_MAP_HPP__

#include <vector>
#include <exception>
#include <iostream>
#include <fstream>
#include <assert.h>
#include <stdexcept>
#include <map>
#include <iomanip>

#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <math/bounding_box.hpp>
#include <util/timer.hpp>
#include <maps/mapNd.hpp>
#include <maps/mapExceptions.hpp>
#include <maps/metric.hpp>
#include <maps/section.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace {
template<typename V>
struct PassiveTracker {
    typedef V  vec_type;
    void initialize(const vec_type&) const {}
    std::vector<double> operator()(const double&, const vec_type&) const
    {
        return std::vector<double>(1,0);
    }
    void mark_crossing() {}
};
}

namespace xavier {

//Note:  Ideally, we should improve the error-handling in this class to be
//consistent throughout, but we are simply trying to pass out the initial state
//that causes a mapping failure at some point downstream.
// FUTURE WORK:  Build a single exception (with different error codes or enums) that
// passes out detailed information to the user.

template<typename RHS, typename ODESolver, typename Section>
class poincare_map : public mapNd<typename RHS::value_type, Section::local_dimension> {
public:

    static const int N = Section::local_dimension;
    
    typedef RHS                                                     rhs_type;
    typedef typename RHS::value_type                                value_type;
    typedef typename RHS::state_type                                state_type;
    typedef typename RHS::xstate_type                               xstate_type;
    typedef Section                                                 section_type;
    typedef typename Section::gvec_type                             gvec_type;
    typedef typename Section::gmat_type                             gmat_type;
    typedef typename Section::lvec_type                             lvec_type;
    typedef typename Section::lmat_type                             lmat_type;
    typedef typename Section::lbox_type                             lbox_type;
    typedef ODESolver                                               solver_type;
    typedef typename solver_type::step                              step_type;
    typedef return_map_info<value_type, N>                          return_type;
    typedef state_info<value_type, RHS::dimension>                  return_state;
    typedef metric<value_type, N>                                   metric_type;
    typedef typename mapNd<value_type, N>::bvec_type                bvec_type;
    typedef poincare_map<RHS, ODESolver, Section>                   self_type;
    typedef FailedIntegration<lvec_type,gvec_type>                  IntegError;
    
    poincare_map(const rhs_type& rhs, const section_type& section, int periodic=0)
        : _rhs(rhs), _Ham(rhs.desired_hamiltonian()), _throwHamError(false), _section(section), _periodic(periodic), _prec(1.0e-8), _max_iter(1000000),
          _verbose(false) {}
          
    poincare_map(const poincare_map& other)
        : _rhs(other._rhs), _Ham(other._Ham), _section(other._section)
    {
        _throwHamError = other._throwHamError;
        _periodic = other._periodic;
        _prec     = other._prec;
        _max_iter = other._max_iter;
        _verbose  = other._verbose;
    }
    
    bool verbose() const
    {
        return _verbose;
    }
    
    bool& verbose()
    {
        return _verbose;
    }
    
    /// Setting the integration precision
    void setPrecision(const double& prec)
    {
        _prec = prec;
    }
    
    const double& precision() const
    {
        return _prec;
    }
    
    size_t& max_iterations()
    {
        return _max_iter;
    }
    
    const size_t& max_iterations() const
    {
        return _max_iter;
    }
    
    
    const section_type& section() const
    {
        return _section;
    }
    
    /// Setting a new section
    void setSection(const section_type& s)
    {
        _section = s;
    }
    
    const rhs_type& rhs() const
    {
        return _rhs;
    }
    
    /// Setting a new right-hand side (EOMs)
    void setRHS(const rhs_type& r)
    {
        _rhs = r;
        setHamiltonian( r.desired_hamiltonian() );
    }
    
    /// Setting a new section and rhs
    void setNewComponents( const rhs_type& r, const section_type& s)
    {
        setRHS( r );
        setSection( s );
    }
    
    /// Return the Hamiltonian value stored in poincare_map object
    const value_type& hamiltonian() const
    {
        return _Ham;
    }
    
    /// Set the Hamiltonian value
    void setHamiltonian( const double& ham)
    {
        _Ham = ham;
    }
    
    
    /// Enable or disable poincare_map from throwing a runtime_error when Hamiltonian is violated
    void stopAtHamiltonianError(const bool ch)
    {
        _throwHamError = ch;
    }
    
    /// Is Hamiltonian violation error being thrown
    bool isCheckingHamiltonian() const
    {
        return _throwHamError;
    }
    
    
    const lbox_type& bounds() const
    {
        return _section.local_bounds();
    }
    
    bvec_type periodic() const
    {
        return bvec_type(false);
    }
    
    self_type* clone() const
    {
        return new self_type(*this);
    }
    
    /// Transform a gvec_type to a state type
    state_type getState(const gvec_type& x) const
    {
        state_type x0(0.0);
        const int rhsDIM = RHS::dimension;
        for(int i=0; i<rhsDIM; i++) {
            x0[i] = x[i];
        }
        return x0;
    }
    
    /// Transform a state_type to a gvec_type (with identity STM)
    gvec_type getXState(const state_type& x) const
    {
        gvec_type x0(0.0);
        const int rhsDIM = RHS::dimension;
        for(int i=0; i<rhsDIM; i++) {
            x0[i] = x[i];
            x0[ rhsDIM + rhsDIM*i + i  ] = 1.0; //STM identity
        }
        return x0;
    }
    
    // Mapping Functions ---------------------------------------------------------------------------
    /// Simplistic map function for given number of iterates, returning only last crossing map state
    lvec_type map(const lvec_type& in, int niter) const
    {
        return_type r = map_complete(in, niter);
        return r.x;
    }
    
    /// Execute Poincare map for niter and return std::vector of returns on section
    void map(const lvec_type& in, std::vector<lvec_type>& out, int niter) const
    {
        std::vector<return_type> tmp;
        PassiveTracker<gvec_type> passive;
        return_map<PassiveTracker<gvec_type> >(in, tmp, niter, passive);
        if ((int) tmp.size() < std::abs(niter)) {
            std::ostringstream os;
            os << tmp.size() << " iterations instead of " << niter << " at " << in << '\n';
            std::cerr << os.str();
            MapUndefined m_u;
            m_u.where = in;
            throw m_u;
        }
        out.resize(tmp.size());
        for (int i=0 ; i<(int)tmp.size() ; ++i) {
            out[i] = tmp[i].x;
        }
    }
    
    /// Execute Poincare map for niter and return std::vector of state_info objects
    void map(const lvec_type& in, std::vector<return_state>& out, int niter) const
    {
        PassiveTracker<gvec_type> passive;
        return_map<PassiveTracker<gvec_type> >(in, out, niter, passive);
        if ((int) out.size() < std::abs(niter)) {
            std::ostringstream os;
            os << out.size() << " iterations instead of " << niter << " at " << in << '\n';
            std::cerr << os.str();
            MapUndefined m_u;
            m_u.where = in;
            throw m_u;
        }
    }
    
    
    /// Execute a Poincare map with a custom state returning map conditions
    void map(const state_type& in, std::vector<return_type>& out, int niter) const
    {
        PassiveTracker<gvec_type> passive;
        return_map<PassiveTracker<gvec_type> >(in, out, niter, passive);
        if ((int) out.size() < std::abs(niter)) {
            std::ostringstream os;
            os << out.size() << " iterations instead of " << niter << " at " << in << '\n';
            std::cerr << os.str();
            MapUndefined m_u;
            gvec_type gv(0);
            for(int i=0; i<rhs_type::dimension; i++) {
                gv[i] = in[i];
            }
            m_u.where = _section.project(gv).first;
            throw m_u;
        }
    }
    
    /// Execute a Poincare map with a custom state returning state info with STM at crossings
    void map(const state_type& in, std::vector<return_state>& out, int niter) const
    {
        PassiveTracker<gvec_type> passive;
        std::vector<return_state> temp;
        return_map<PassiveTracker<gvec_type> >(in, out, temp, niter, passive);
        if ((int) out.size() < std::abs(niter)) {
            std::ostringstream os;
            os << out.size() << " iterations instead of " << niter << " at " << in << '\n';
            std::cerr << os.str();
            MapUndefined m_u;
            gvec_type gv(0);
            for(int i=0; i<rhs_type::dimension; i++) {
                gv[i] = in[i];
            }
            m_u.where = _section.project(gv).first;
            throw m_u;
        }
    }
    
    
    void map_complete(const lvec_type& in, std::vector<return_type>& out, int niter) const
    {
        PassiveTracker<gvec_type> passive;
        return_map<PassiveTracker<gvec_type> >(in, out, niter, passive);
        if ((int) out.size() < std::abs(niter)) {
            MapUndefined m_u;
            m_u.where = in;
            throw m_u;
        }
    }
    
    void map_complete(const lvec_type& in, std::vector<return_state>& stateout, int niter) const
    {
        PassiveTracker<gvec_type> passive;
        return_map<PassiveTracker<gvec_type> >(in, stateout, niter, passive);
        if ((int)stateout.size() < std::abs(niter)) {
            MapUndefined m_u;
            m_u.where = in;
            throw m_u;
        }
    }
    
    return_type map_complete(const lvec_type& in, int niter) const
    {
        std::vector<return_type> out;
        PassiveTracker<gvec_type> passive;
        return_map<PassiveTracker<gvec_type> >(in, out, niter, passive);
        if ((int) out.size() < std::abs(niter)) {
            MapUndefined m_u;
            m_u.where = in;
            throw m_u;
        }
        return out.back();
    }
    
    /// Mapping function with Tracking System - returns only return_type at last crossing
    template<typename Tracker>
    return_type map_and_track_complete(const lvec_type& in, int niter, Tracker& tracker) const
    {
        std::vector<return_type> out;
        return_map<Tracker>(in, out, niter, tracker);
        if ((int) out.size() < std::abs(niter)) {
            MapUndefined m_u;
            m_u.where = in;
            throw m_u;
        }
        return out.back();
    }
    
    template<typename Tracker>
    void map_and_track_complete(const lvec_type& in, std::vector<return_type>& out, int niter, Tracker& tracker) const
    {
        return_map<Tracker>(in, out, niter, tracker);
        if ((int) out.size() < std::abs(niter)) {
            MapUndefined m_u;
            m_u.where = in;
            throw m_u;
        }
    }
    
    template<typename Tracker>
    void map_and_track_complete(const lvec_type& in, std::vector<lvec_type>& out, int niter, Tracker& tracker) const
    {
        std::vector<return_type> tmp;
        return_map<Tracker>(in,tmp,niter,tracker);
        if ((int)tmp.size() < std::abs(niter)) {
            std::ostringstream os;
            os << tmp.size() << " iterations instead of " << niter << " at " << in << '\n';
            std::cerr << os.str();
            MapUndefined m_u;
            m_u.where = in;
            throw m_u;
        }
        out.resize(tmp.size());
        for (int i=0 ; i<(int)tmp.size() ; ++i) {
            out[i] = tmp[i].x;
        }
    }
    
    template<typename Tracker>
    void map_and_track_complete(const lvec_type& x0, std::vector<return_state>& itersOut,
                                std::vector<return_state>& internalOut, int niter, Tracker& tracker) const
    {
        gvec_type y0 = _section.unproject(x0);
        state_type ys(0.0);
        for(int i=0; i<rhs_type::dimension; i++) {
            ys[i] = y0[i];
        }
        return_map<Tracker>(ys, itersOut, internalOut, niter, tracker, 0.0, true);
        if ((int) itersOut.size() < std::abs(niter)) {
            std::ostringstream os;
            os << itersOut.size() << " iterations instead of " << niter << " at " << x0 << '\n';
            std::cerr << os.str();
            MapUndefined m_u;
            m_u.where = x0;
            throw m_u;
        }
    }
    
    /// Return state_f and STM for a specified time given a full (N+N*N) state
    return_state integrate_state(const gvec_type& y0, double t) const
    {
        return_state out;
        try {
            _integrate_state(y0,out,t);
        } catch(IntegError e) {
            throw e;
        } catch(...) {
            throw("Unkown Error");
        }
        return out;
    }
    
    /// Return state_f and STM for a specified time given a state
    return_state integrate_state(const state_type& x0, double t) const
    {
        return_state out;
        //Convert to a full state by using identity for STM
        gvec_type y0(0.0);
        const int rhsDIM = RHS::dimension;
        for(int i=0; i<rhsDIM; i++) {
            y0[i] = x0[i];
            y0[ rhsDIM + rhsDIM*i + i  ] = 1.0; //STM identity
        }
        try {
            _integrate_state(y0,out,t);
        } catch(IntegError e) {
            throw e;
        } catch(...) {
            throw("Unkown Error");
        }
        return out;
    }
    
    /// Return state_f and STM for a specified time given on-map states
    return_state flow_map(const lvec_type& in, double t) const
    {
        return_state out;
        try {
            _flow_map(in, out, t);
        } catch(IntegError e) {
            throw e;
        } catch(...) {
            throw("Unknown Error");
        }
        return out;
    }
    
    /// Return all states and STMs for a specified map state and given number of iterates
    void flow_map(const lvec_type& x0, std::vector<return_state>& itersOut,
                  std::vector<return_state>& internalOut, int niter) const
    {
        gvec_type y0 = _section.unproject(x0);
        state_type ys(0.0);
        for(int i=0; i<rhs_type::dimension; i++) {
            ys[i] = y0[i];
        }
        PassiveTracker<gvec_type> passive;
        return_map<PassiveTracker<gvec_type> >(ys, itersOut, internalOut, niter, passive, 0.0, true);
        if ((int) itersOut.size() < std::abs(niter)) {
            std::ostringstream os;
            os << itersOut.size() << " iterations instead of " << niter << " at " << x0 << '\n';
            std::cerr << os.str();
            MapUndefined m_u;
            m_u.where = x0;
            throw m_u;
        }
    }
    
    /// Return all states and STMs for a specified state and given number of iterates
    void flow_map(const state_type& x0, std::vector<return_state>& itersOut,
                  std::vector<return_state>& internalOut, int niter) const
    {
        PassiveTracker<gvec_type> passive;
        return_map<PassiveTracker<gvec_type> >(x0, itersOut, internalOut, niter, passive, 0.0, true);
        if ((int) itersOut.size() < std::abs(niter)) {
            std::ostringstream os;
            os << itersOut.size() << " iterations instead of " << niter << " at " << x0 << '\n';
            std::cerr << os.str();
            MapUndefined m_u;
            gvec_type gv(0);
            for(int i=0; i<rhs_type::dimension; i++) {
                gv[i] = x0[i];
            }
            m_u.where = _section.project(gv).first;
            throw m_u;
        }
    }
    
    /// Return all states and STMs for a specified on-map state and number of iterates
    void flow_map(const lvec_type& in, std::vector<return_state>& out, int niter) const
    {
        PassiveTracker<gvec_type> passive;
        //turn on "dense" output for total time-history
        return_map<PassiveTracker<gvec_type> >(in, out, niter, passive, 0.0, true);
        if ((int)out.size() < std::abs(niter)) {
            MapUndefined m_u;
            m_u.where = in;
            throw m_u;
        }
    }
    
    /// Return all states and Map iterates for a specified on-map state and number of iters
    void flow_map(const lvec_type& in, std::vector<return_type>& itersOut,
                  std::vector<return_state>& internalOut, int niter) const
    {
        PassiveTracker<gvec_type> passive;
        //turn on "dense output for total time-history
        return_map<PassiveTracker<gvec_type> >(in, itersOut, internalOut, niter, passive, 0.0, true);
        if ((int)itersOut.size() < std::abs(niter) ) {
            MapUndefined m_u;
            m_u.where = in;
            throw m_u;
        }
    }
    
    /// Return all states and Map iterates for a specified on-map state and time of flight
    void flow_map(const lvec_type& in, const double t, std::vector<return_type>& itersOut,
                  std::vector<return_state>& internalOut) const
    {
        PassiveTracker<gvec_type> passive;
        //turn on "dense output for total time-history
        return_map<PassiveTracker<gvec_type> >(in, t, itersOut, internalOut, passive, 0.0, true);
        //Check if it reached appropriate time
        if ((t>0 && internalOut.back().t<t) || (t<0 && internalOut.back().t>t)) {
            MapUndefined m_u;
            m_u.where = in;
            throw m_u;
        }
    }
    
    /** Mapping function that resets the STM after each crossing
     *  - construct total through multiplication
     *  - Note: not particularly consistent as crossings won't be at equal time steps
     */
    void map_STMperReturn(const lvec_type& in, std::vector<return_state>& stateout, int niter) const
    {
        PassiveTracker<gvec_type> passive;
        return_map_STM<PassiveTracker<gvec_type> >(in, stateout, niter, passive);
        if ((int)stateout.size() < std::abs(niter)) {
            MapUndefined m_u;
            m_u.where = in;
            throw m_u;
        }
    }
    
    /** Run the map to compute better approximations of the linearSTM's by limiting
     *   the maximum time step between propagation nodes.  The STM is then returned
     *   at each crossing as the multiplication of STMs.
     *   -Note: Usually use this for fixed points!
     *   -Note: For p=3, tCrossings has 3 values at p=1,p=2,p=3.  Do NOT store p=0 (t=0).
     */
    void map_LinearMonodromy(const lvec_type& in, const double& dtMax,
                             const std::vector<double>& tCrossings,
                             std::vector<return_state>& stateout) const
    {
        return_monodromy(in, dtMax, tCrossings, stateout);
    }
    
    
    /// Evaluate a given number of hits (assuming nHits>=maxIters) to a section region
    void mapRegion(const lvec_type& in, std::vector<return_type>& out,
                   int nHits, int maxIters, const lbox_type& regionBox) const
    {
        PassiveTracker<gvec_type> passive;
        compute_hits<PassiveTracker<gvec_type> >(in,out,nHits,maxIters,regionBox,passive);
        if ((int)out.size() < std::abs(nHits)) {
            std::ostringstream os;
            os << "Insufficient hits : " << out.size() << " instead of "
               << nHits << " at " << in << "\n";
            if(_verbose) {
                std::cerr << os.str();
            }
            MapUndefined m_u;
            m_u.where = in;
            throw m_u;
        }
        
    }
    
    
    /// Approx winding number with weighted periodic subtraction (periodic Poincare map only!)
    double winding_number(const lvec_type& x0, const std::vector<lvec_type>& orbit) const
    {
        //Note: This does NOT WORK for systems like CR3BP
        double w = 0;
        double number = 0;
        for (size_t i=0 ; i<orbit.size() ; ++i) {
            const lvec_type& x = orbit[i];
            double d = metric_type::periodic_subtraction(x0[_periodic], x[_periodic], bounds().size()[_periodic]);
            if (d == 0) {
                return (x[_periodic] - x0[_periodic])/(float)(i+1);
            } else {
                w += 1./d; // weight is inverse propertional to periodic distance to x0
                number += (x[_periodic] - x0[_periodic])/(float)(i+1)/d;
            }
        }
        return number/w;
    }
    
private:
    bool stepHamiltonianCheck(const step_type& step) const;
    template<typename Tracker>
    void return_map(const lvec_type& seed, std::vector<return_type>& out,
                    int niter, Tracker& tracker, double h=0) const;
    template<typename Tracker>
    void return_map(const lvec_type& seed, std::vector<return_state>& out,
                    int niter, Tracker& tracker, double h=0, bool dense=false) const;
    template<typename Tracker>
    void return_map(const lvec_type& seed, std::vector<return_type>& itersOut,
                    std::vector<return_state>& internalOut, int niter,
                    Tracker& tracker, double h=0, bool dense=false) const;
    template<typename Tracker>
    void return_map(const state_type& seed, std::vector<return_type>& out,
                    int niter, Tracker& tracker, double h=0) const;
    template<typename Tracker>
    void return_map(const state_type& seed, std::vector<return_state>& itersOut,
                    std::vector<return_state>& internalOut, int niter,
                    Tracker& tracker, double h=0, bool dense=false) const;
    template<typename Tracker>
    void return_map(const lvec_type& seed, const double& t,
                    std::vector<return_type>& itersOut,
                    std::vector<return_state>& internalOut,
                    Tracker& tracker, double h=0, bool dense=false) const;
    ///Note: This is not a great function!
    template<typename Tracker>
    void return_map_STM(const lvec_type& seed, std::vector<return_state>& out,
                        int niter, Tracker& tracker, double h=0, bool dense=false) const;
    ///Computing a better linear STM approximation limited by dtMax
    void return_monodromy(const lvec_type& seed, const double& dtMax,
                          const std::vector<double>& tVals, std::vector<return_state>& out
                         ) const;
    template<typename Tracker>
    void compute_hits(const lvec_type& seed, std::vector<return_type>& hitsOut,
                      int nHits, int maxIters, const lbox_type& regionBB, Tracker& tracker, double h=0) const;
    void _flow_map(const lvec_type& seed, return_state& out, double t) const;
    
    void _integrate_state(const gvec_type& y0, return_state& out, double t) const;
    
    rhs_type        _rhs;
    value_type      _Ham;
    bool            _throwHamError;
    section_type    _section;
    int             _periodic;
    double          _prec;
    size_t          _max_iter;
    mutable bool    _verbose;
};


namespace {

template<typename Section>
class intersect_stop {
public:
    typedef Section                                   section_type;
    typedef typename section_type::gvec_type          vec_type;
    typedef typename section_type::gmat_type          mat_type;
    
private:
    const section_type&                _sec;
    bool                               _stopped;
    vec_type                           _inter;
    double                             _inter_time;
    double                             _sign;
    
public:
    intersect_stop(const section_type& sec, bool fwd=true) :
        _sec(sec), _stopped(false), _sign(fwd? 1 : -1)
    {
    }
    
    static double secant_method(double v0, double t0, double v1, double t1)
    {
        return (v1*t0 - v0*t1)/(v1-v0);
    }
    
    template<typename STEP>
    bool operator()(const STEP& s, bool verbose=false)
    {
        _stopped = false;
        double t0 = std::min(s.t0(), s.t1());
        double t1 = std::max(s.t0(), s.t1());
        double f0 = _sign*_sec.distance(s.y(t0));
        double f1 = _sign*_sec.distance(s.y(t1));
        
        std::ostringstream os;
        
        if (verbose) {
            os << "poincare: stop: f0=" << f0 << ", f1=" << f1 << ", t0=" << t0 << ", t1=" << t1 << ", dt=" << t1-t0 << '\n';
            std::cerr << os.str();
            os.str("");
            os.clear();
        }
        
        //If appropriate sign change occurs, compute the crossing
        bool fwd = (_sign > 0);
        if ((fwd && f0 < 0 && f1 > 0) || (!fwd && f0 > 0 && f1 < 0)) {
            double dt, t, df, f;
            dt = t1 - t0;
            df = f1 - f0;
            
            // apply secant method
            f = 1; // dummy;
            while (fabs(f) > 1e-8) {
                t = secant_method(f0, t0, f1, t1);
                f = _sign*_sec.distance(s.y(t));
                
                if (f > 0) {
                    t1 = t;
                    f1 = f;
                } else {
                    t0 = t;
                    f0 = f;
                }
                
                double _dt = fabs(t1-t0);
                double _df = fabs(f1-f0);
                if (_dt >= dt || _df >= df) {
                    break;
                }
                dt = _dt;
                df = _df;
                
                if (verbose) {
                    os << "poincare: stop: " << std::setprecision(20) << "t=" << t << ", f=" << f << ", dt=" << dt << ", df=" << df << "\n";
                    std::cerr << os.str();
                    os.str("");
                    os.clear();
                }
            }
            
            _inter = s.y(t);
            _inter_time = t;
            _stopped = true;
        }
        
        return _stopped;
    }
    
    bool did_stop() const
    {
        return _stopped;
    }
    
    const vec_type& where() const
    {
        return _inter;
    }
    
    double when() const
    {
        return _inter_time;
    }
};

} // anonymous namespace

// --------------------------------------------------------------------------
///Check for Hamiltonian violation (useful for detecting failure due to singularities)
template<typename RHS, typename ODESolver, typename Section>
bool poincare_map<RHS, ODESolver, Section>::
stepHamiltonianCheck(const typename ODESolver::step& step) const
{
    value_type ham1 = _rhs.hamiltonian(step.t1(), step.y1());
    if (fabs(ham1-_Ham) > 100.*_prec) {
        //Step induced too much error in the Hamiltonian (which should be constant)
        std::cerr << " Hamiltonian evaluates to " << ham1 << " when it should be " << _Ham
                  << " (dH = " << fabs(ham1-_Ham) << ")\n";
        return false;
    }
    //Hamiltonian is ok
    return true;
}


///Return Map -> sends back vector of return_map_info
template<typename RHS, typename ODESolver, typename Section>
template<typename Tracker>
void poincare_map<RHS, ODESolver, Section>::
return_map(const lvec_type& seed, std::vector<return_type>& out, int niter, Tracker& tracker, double h) const
{
    typedef ODESolver                solver_type;
    typename solver_type::result     res;
    typename solver_type::step       step;
    
    out.clear();
    out.reserve(std::abs(niter));
    
    solver_type intg;
    // set initial condition of ODE solver
    gvec_type x0 = _section.unproject(seed);
    tracker.initialize(x0);
    intg.set_init_cond(x0, 0.);
    // sign of tmax determines direction
    intg.set_t_max(niter < 0 ? -1e8 : 1e8);
    intg.set_init_step(1e-5);
    intg.set_precision(_prec);
    intersect_stop<section_type> stop(_section, niter>0);
    nvis::timer _t;
    
    for (size_t counter=0 ; counter<_max_iter && (int)out.size() < std::abs(niter) ; ++counter) {
    
        if (_verbose) {
            std::ostringstream os;
            os << "poincare at " << intg.y << " / ";
            if ((int) out.size()) {
                os << out.back().x;
            } else {
                os << seed;
            }
            os << " after " << counter << " steps and " << _t.elapsed() << " seconds\n";
            std::cerr << os.str();
        }
        
        try {
            res = intg.do_step(_rhs, step);
        } catch (MapSingularityError<state_type>& msError) {
            if(_verbose) {
                std::cerr << "Caught Singularity in RHS at " << msError.where << std::endl;
            }
            throw msError;
        } catch (...) {
            if(_verbose) {
                std::cerr << "caught exception in return_map" << std::endl;
            }
            break;
        }
        //Check for Hamiltonian violation (useful for detecting failure due to singularities)
        if (_throwHamError && !stepHamiltonianCheck(step) ) {
            if(_verbose) {
                std::cerr << "constant Hamiltonian violated in return_map" << std::endl;
            }
            value_type ham1 = _rhs.hamiltonian(step.t1(), step.y1());
            if(_verbose) std::cerr << " Hamiltonian evaluates to " << ham1 << " when it should be " << _Ham
                                       << " (dH = " << fabs(ham1-_Ham) << ", prec=" << _prec << ")\n";
                                       
            break;
        }
        
        if (_verbose) {
            std::ostringstream os;
            os << "poincare: step completed\n";
            std::cerr << os.str();
        }
        
        if (res != solver_type::OK) {
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: integration ended with result: " << (int)res << '\n';
                std::cerr << os.str();
            }
            break;
        }
        if (nvis::norm(step.y0() - step.y1())<1.0e-8) {
            std::ostringstream os;
            os << "poincare: step size undeflow\n";
            std::cerr << os.str();
            break;
        }
        
        stop(step, _verbose);
        if (_verbose) {
            std::ostringstream os;
            os << "poincare: stop check completed\n";
            std::cerr << os.str();
        }
        
        if (stop.did_stop()) {
            if (_verbose) {
                std::ostringstream os;
                os << "**** poincare: stop detected ****\n";
                std::cerr << os.str();
            }
            out.push_back(return_type());
            return_type& r = out.back();
            std::pair<lvec_type, lmat_type> tmp = _section.project(stop.where());
            r.x = tmp.first;
            r.J = tmp.second;
            r.delta_theta = tracker(stop.when(),stop.where());
            tracker.mark_crossing();
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: intersection position computed\n";
                std::cerr << os.str();
            }
            r.t = stop.when();
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: intersection time stored\n";
                std::cerr << os.str();
            }
            
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: stop information stored\n";
                std::cerr << os.str();
            }
        }
        
        tracker(step.t1(),step.y1());
        
        if (_verbose) {
            std::ostringstream os;
            os << "poincare: Jacobian integration completed\n";
            std::cerr << os.str();
        }
    }
}

///Return Map -> sends back vector of state_info
template<typename RHS, typename ODESolver, typename Section>
template<typename Tracker>
void poincare_map<RHS, ODESolver, Section>::
return_map(const lvec_type& seed, std::vector<return_state>& out, int niter, Tracker& tracker,
           double h, bool dense) const
{
    typedef ODESolver                solver_type;
    typename solver_type::result     res;
    typename solver_type::step       step;
    
    out.clear();
    out.reserve(std::abs(niter));
    if (dense) {
        out.reserve(std::abs(niter)*1000);
    }
    size_t iterCount = 0;
    
    solver_type intg;
    // set initial condition of ODE solver
    gvec_type x0 = _section.unproject(seed);
    tracker.initialize(x0);
    intg.set_init_cond(x0, 0.);
    // sign of tmax determines direction
    intg.set_t_max(niter < 0 ? -1e8 : 1e8);
    intg.set_init_step(1e-5);
    intg.set_precision(_prec);
    intersect_stop<section_type> stop(_section, niter>0);
    nvis::timer _t;
    
    for (size_t counter=0 ; counter<_max_iter && (int)iterCount < std::abs(niter) ; ++counter) {
    
        if (_verbose) {
            std::ostringstream os;
            os << "poincare at " << intg.y << " / ";
            if (out.size()) {
                os << out.back().x;
            } else {
                os << seed;
            }
            os << " after " << counter << " steps and " << _t.elapsed() << " seconds\n";
            std::cerr << os.str();
        }
        
        try {
            res = intg.do_step(_rhs, step);
        } catch (MapSingularityError<state_type>& msError) {
            if(_verbose) {
                std::cerr << "Caught Singularity in RHS at " << msError.where << std::endl;
            }
            throw msError;
        } catch (...) {
            if(_verbose) {
                std::cerr << "caught exception in return_map" << std::endl;
            }
            break;
        }
        //Check for Hamiltonian violation (useful for detecting failure due to singularities)
        if ( _throwHamError && !stepHamiltonianCheck(step) ) {
            if(_verbose) {
                std::cerr << "constant Hamiltonian violated in return_map" << std::endl;
            }
            value_type ham1 = _rhs.hamiltonian(step.t1(), step.y1());
            if(_verbose) std::cerr << " Hamiltonian evaluates to " << ham1 << " when it should be " << _Ham
                                       << " (dH = " << fabs(ham1-_Ham) << ", prec=" << _prec << ")\n";
            break;
        }
        if (_verbose) {
            std::ostringstream os;
            os << "poincare: step completed\n";
            std::cerr << os.str();
        }
        
        if (res != solver_type::OK) {
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: integration ended with result: " << (int)res << '\n';
                std::cerr << os.str();
            }
            break;
        }
        if (nvis::norm(step.y0() - step.y1())<1.0e-8) {
            std::ostringstream os;
            os << "poincare: step size undeflow\n";
            std::cerr << os.str();
            break;
        }
        
        //Check for a stopping condition
        stop(step, _verbose);
        if (_verbose) {
            std::ostringstream os;
            os << "poincare: stop check completed\n";
            std::cerr << os.str();
        }
        
        if (stop.did_stop()) {
            if (_verbose) {
                std::ostringstream os;
                os << "**** poincare: stop detected ****\n";
                std::cerr << os.str();
            }
            out.push_back(return_state());
            return_state& r = out.back();
            r.setState(stop.where());
            tracker.mark_crossing();
            iterCount++;
            
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: intersection position computed\n";
                std::cerr << os.str();
            }
            r.t = stop.when();
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: intersection time stored\n";
                std::cerr << os.str();
            }
            
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: stop information stored\n";
                std::cerr << os.str();
            }
        }
        
        //Store info if dense output is on
        if (dense) {
            //Store every step
            out.push_back(return_state());
            return_state& r = out.back();
            r.setState(step.y1()); //Store the step
            r.t = step.t1();
        }
        tracker(step.t1(),step.y1());
        
        if (_verbose) {
            std::ostringstream os;
            os << "poincare: Jacobian integration completed\n";
            std::cerr << os.str();
        }
    }
}

///Return Map -> sends back iterates in itersOut and dense output vector of state_info
template<typename RHS, typename ODESolver, typename Section>
template<typename Tracker>
void poincare_map<RHS, ODESolver, Section>::
return_map(const lvec_type& seed, std::vector<return_type>& itersOut,
           std::vector<return_state>& internalOut, int niter, Tracker& tracker,
           double h, bool dense) const
{
    typedef ODESolver                solver_type;
    typename solver_type::result     res;
    typename solver_type::step              step;
    
    itersOut.clear();
    itersOut.reserve(std::abs(niter));
    internalOut.clear();
    internalOut.reserve(std::abs(niter));
    if (dense) {
        internalOut.reserve(std::abs(niter)*1000);
    }
    size_t iterCount = 0;
    
    solver_type intg;
    // set initial condition of ODE solver
    gvec_type x0 = _section.unproject(seed);
    tracker.initialize(x0);
    intg.set_init_cond(x0, 0.);
    // sign of tmax determines direction
    intg.set_t_max(niter < 0 ? -1e8 : 1e8);
    intg.set_init_step(1e-5);
    intg.set_precision(_prec);
    intersect_stop<section_type> stop(_section, niter>0);
    nvis::timer _t;
    
    for (size_t counter=0 ; counter<_max_iter && (int)iterCount < std::abs(niter) ; ++counter) {
    
        if (_verbose) {
            std::ostringstream os;
            os << "poincare at " << intg.y << " / ";
            if (itersOut.size()) {
                os << itersOut.back().x;
            } else {
                os << seed;
            }
            os << " after " << counter << " steps and " << _t.elapsed() << " seconds\n";
            std::cerr << os.str();
        }
        
        try {
            res = intg.do_step(_rhs, step);
        } catch (MapSingularityError<state_type>& msError) {
            if(_verbose) {
                std::cerr << "Caught Singularity in RHS at " << msError.where << std::endl;
            }
            throw msError;
        } catch (...) {
            if(_verbose) {
                std::cerr << "caught exception in return_map" << std::endl;
            }
            break;
        }
        //Check for Hamiltonian violation (useful for detecting failure due to singularities)
        if ( _throwHamError && !stepHamiltonianCheck(step) ) {
            if(_verbose) {
                std::cerr << "constant Hamiltonian violated in return_map" << std::endl;
            }
            value_type ham1 = _rhs.hamiltonian(step.t1(), step.y1());
            if(_verbose) std::cerr << " Hamiltonian evaluates to " << ham1 << " when it should be " << _Ham
                                       << " (dH = " << fabs(ham1-_Ham) << ", prec=" << _prec << ")\n";
            break;
        }
        if (_verbose) {
            std::ostringstream os;
            os << "poincare: step completed\n";
            std::cerr << os.str();
        }
        
        if (res != solver_type::OK) {
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: integration ended with result: " << (int)res << '\n';
                std::cerr << os.str();
            }
            break;
        }
        if (nvis::norm(step.y0() - step.y1())<1.0e-8) {
            std::ostringstream os;
            os << "poincare: step size undeflow\n";
            std::cerr << os.str();
            break;
        }
        
        //Check for a stopping condition
        stop(step, _verbose);
        if (_verbose) {
            std::ostringstream os;
            os << "poincare: stop check completed\n";
            std::cerr << os.str();
        }
        
        if (stop.did_stop()) {
            if (_verbose) {
                std::ostringstream os;
                os << "**** poincare: stop detected ****\n";
                std::cerr << os.str();
            }
            //Return info
            itersOut.push_back(return_type());
            return_type& r = itersOut.back();
            std::pair<lvec_type, lmat_type> tmp = _section.project(stop.where());
            r.x = tmp.first;
            r.J = tmp.second;
            r.delta_theta = tracker(stop.when(),stop.where());
            //state_info
            internalOut.push_back(return_state());
            return_state& rs = internalOut.back();
            rs.setState(stop.where());
            tracker.mark_crossing();
            iterCount++;
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: intersection position computed\n";
                std::cerr << os.str();
            }
            r.t = stop.when();
            rs.t = stop.when();
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: intersection time stored\n";
                std::cerr << os.str();
            }
            
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: stop information stored\n";
                std::cerr << os.str();
            }
        }
        
        //Store info if dense output is on
        if (dense) {
            //Store every step
            internalOut.push_back(return_state());
            return_state& rs = internalOut.back();
            rs.setState(step.y1()); //Store the step
            rs.t = step.t1();
        }
        tracker(step.t1(),step.y1());
        
        if (_verbose) {
            std::ostringstream os;
            os << "poincare: Jacobian integration completed\n";
            std::cerr << os.str();
        }
    }
}

///Return Map -> sends back iterates in itersOut and dense output vector of state_info
template<typename RHS, typename ODESolver, typename Section>
template<typename Tracker>
void poincare_map<RHS, ODESolver, Section>::
return_map(const lvec_type& seed, const double& t, std::vector<return_type>& itersOut,
           std::vector<return_state>& internalOut, Tracker& tracker,
           double h, bool dense) const
{
    typedef ODESolver                solver_type;
    typename solver_type::result     res;
    typename solver_type::step       step;
    
    int niter = (int)sign(t)*5;
    itersOut.clear();
    itersOut.reserve(std::abs(niter));
    internalOut.clear();
    internalOut.reserve(std::abs(niter));
    if (dense) {
        internalOut.reserve(std::abs(niter)*1000);
    }
    size_t iterCount = 0;
    
    solver_type intg;
    // set initial condition of ODE solver
    gvec_type x0 = _section.unproject(seed);
    tracker.initialize(x0);
    intg.set_init_cond(x0, 0.);
    // sign of tmax determines direction
    intg.set_t_max(niter < 0 ? -1e8 : 1e8);
    intg.set_init_step(1e-5);
    intg.set_precision(_prec);
    intersect_stop<section_type> stop(_section, niter>0);
    nvis::timer _t;
    
    for (size_t counter=0 ; counter<_max_iter ; ++counter) {
    
        if (_verbose) {
            std::ostringstream os;
            os << "poincare at " << intg.y << " / ";
            if (itersOut.size()) {
                os << itersOut.back().x;
            } else {
                os << seed;
            }
            os << " after " << counter << " steps and " << _t.elapsed() << " seconds\n";
            std::cerr << os.str();
        }
        
        try {
            res = intg.do_step(_rhs, step);
        } catch (MapSingularityError<state_type>& msError) {
            if(_verbose) {
                std::cerr << "Caught Singularity in RHS at " << msError.where << std::endl;
            }
            throw msError;
        } catch (...) {
            if(_verbose) {
                std::cerr << "caught exception in return_map" << std::endl;
            }
            break;
        }
        //Check for Hamiltonian violation (useful for detecting failure due to singularities)
        if ( _throwHamError && !stepHamiltonianCheck(step) ) {
            if(_verbose) {
                std::cerr << "constant Hamiltonian violated in return_map" << std::endl;
            }
            value_type ham1 = _rhs.hamiltonian(step.t1(), step.y1());
            if(_verbose) std::cerr << " Hamiltonian evaluates to " << ham1 << " when it should be " << _Ham
                                       << " (dH = " << fabs(ham1-_Ham) << ", prec=" << _prec << ")\n";
            break;
        }
        if (_verbose) {
            std::ostringstream os;
            os << "poincare: step completed\n";
            std::cerr << os.str();
        }
        
        if (res != solver_type::OK) {
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: integration ended with result: " << (int)res << '\n';
                std::cerr << os.str();
            }
            break;
        }
        if (nvis::norm(step.y0() - step.y1())<1.0e-8) {
            std::ostringstream os;
            os << "poincare: step size undeflow\n";
            std::cerr << os.str();
            break;
        }
        
        //Check for a stopping condition
        stop(step, _verbose);
        if (_verbose) {
            std::ostringstream os;
            os << "poincare: stop check completed\n";
            std::cerr << os.str();
        }
        
        if (stop.did_stop()) {
            if (_verbose) {
                std::ostringstream os;
                os << "**** poincare: stop detected ****\n";
                std::cerr << os.str();
            }
            //Return info
            itersOut.push_back(return_type());
            return_type& r = itersOut.back();
            std::pair<lvec_type, lmat_type> tmp = _section.project(stop.where());
            r.x = tmp.first;
            r.J = tmp.second;
            r.delta_theta = tracker(stop.when(),stop.where());
            //state_info
            internalOut.push_back(return_state());
            return_state& rs = internalOut.back();
            rs.setState(stop.where());
            tracker.mark_crossing();
            iterCount++;
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: intersection position computed\n";
                std::cerr << os.str();
            }
            r.t = stop.when();
            rs.t = stop.when();
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: intersection time stored\n";
                std::cerr << os.str();
            }
            
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: stop information stored\n";
                std::cerr << os.str();
            }
        }
        
        //Store info if dense output is on
        if (dense) {
            //Store every step
            internalOut.push_back(return_state());
            return_state& rs = internalOut.back();
            rs.setState(step.y1()); //Store the step
            rs.t = step.t1();
        }
        tracker(step.t1(),step.y1());
        
        //Check if we are done with the provided time
        if ((t>0 && step.t1()>=t) || (t<0 && step.t1()<=t)) {
            break;
        }
    }
    
    if (_verbose) {
        std::ostringstream os;
        os << "poincare: Jacobian integration completed\n";
        std::cerr << os.str();
    }
}


///Return Map -> sends back vector of return_map_info but input is a state_type
template<typename RHS, typename ODESolver, typename Section>
template<typename Tracker>
void poincare_map<RHS, ODESolver, Section>::
return_map(const state_type& seed, std::vector<return_type>& out, int niter, Tracker& tracker, double h) const
{
    typedef ODESolver                solver_type;
    typename solver_type::result     res;
    typename solver_type::step              step;
    
    out.clear();
    out.reserve(std::abs(niter));
    
    solver_type intg;
    
    // set initial condition of ODE solver
    gvec_type x0(0.0);
    const int rhsDIM = RHS::dimension;
    for(int i=0; i<rhsDIM; i++) {
        x0[i] = seed[i];
        x0[ rhsDIM + rhsDIM*i + i  ] = 1.0; //STM identity
    }
    
    tracker.initialize(x0);
    intg.set_init_cond(x0, 0.);
    // sign of tmax determines direction
    intg.set_t_max(niter < 0 ? -1e8 : 1e8);
    intg.set_init_step(1e-5);
    intg.set_precision(_prec);
    intersect_stop<section_type> stop(_section, niter>0);
    nvis::timer _t;
    
    for (size_t counter=0 ; counter<_max_iter && (int)out.size() < std::abs(niter) ; ++counter) {
    
        if (_verbose) {
            std::ostringstream os;
            os << "poincare at " << intg.y << " / ";
            if (out.size()) {
                os << out.back().x;
            } else {
                os << seed;
            }
            os << " after " << counter << " steps and " << _t.elapsed() << " seconds\n";
            std::cerr << os.str();
        }
        
        try {
            res = intg.do_step(_rhs, step);
        } catch (MapSingularityError<state_type>& msError) {
            if(_verbose) {
                std::cerr << "Caught Singularity in RHS at " << msError.where << std::endl;
            }
            throw msError;
        } catch (...) {
            if(_verbose) {
                std::cerr << "caught exception in return_map" << std::endl;
            }
            break;
        }
        //Check for Hamiltonian violation (useful for detecting failure due to singularities)
        if (_throwHamError && !stepHamiltonianCheck(step) ) {
            if(_verbose) {
                std::cerr << "constant Hamiltonian violated in return_map" << std::endl;
            }
            value_type ham1 = _rhs.hamiltonian(step.t1(), step.y1());
            if(_verbose) std::cerr << " Hamiltonian evaluates to " << ham1 << " when it should be " << _Ham
                                       << " (dH = " << fabs(ham1-_Ham) << ", prec=" << _prec << ")\n";
            break;
        }
        if (_verbose) {
            std::ostringstream os;
            os << "poincare: step completed\n";
            std::cerr << os.str();
        }
        
        if (res != solver_type::OK) {
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: integration ended with result: " << (int)res << '\n';
                std::cerr << os.str();
            }
            break;
        }
        if (nvis::norm(step.y0() - step.y1())<1.0e-8) {
            std::ostringstream os;
            os << "poincare: step size undeflow\n";
            std::cerr << os.str();
            break;
        }
        
        stop(step, _verbose);
        if (_verbose) {
            std::ostringstream os;
            os << "poincare: stop check completed\n";
            std::cerr << os.str();
        }
        
        if (stop.did_stop()) {
            if (_verbose) {
                std::ostringstream os;
                os << "**** poincare: stop detected ****\n";
                std::cerr << os.str();
            }
            out.push_back(return_type());
            return_type& r = out.back();
            std::pair<lvec_type, lmat_type> tmp = _section.project(stop.where());
            r.x = tmp.first;
            r.J = tmp.second;
            r.delta_theta = tracker(stop.when(),stop.where());
            tracker.mark_crossing();
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: intersection position computed\n";
                std::cerr << os.str();
            }
            r.t = stop.when();
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: intersection time stored\n";
                std::cerr << os.str();
            }
            
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: stop information stored\n";
                std::cerr << os.str();
            }
        }
        
        tracker(step.t1(),step.y1());
        
        if (_verbose) {
            std::ostringstream os;
            os << "poincare: Jacobian integration completed\n";
            std::cerr << os.str();
        }
    }
}

///Return Map -> sends back vector of return_map_info but input is a state_type
template<typename RHS, typename ODESolver, typename Section>
template<typename Tracker>
void poincare_map<RHS, ODESolver, Section>::
return_map(const state_type& seed, std::vector<return_state>& out,
           std::vector<return_state>& internalOut,
           int niter, Tracker& tracker,
           double h, bool dense) const
{
    typedef ODESolver                solver_type;
    typename solver_type::result     res;
    typename solver_type::step       step;
    
    out.clear();
    out.reserve(std::abs(niter));
    internalOut.clear();
    internalOut.reserve(std::abs(niter));
    if (dense) {
        internalOut.reserve(std::abs(niter)*1000);
    }
    size_t iterCount = 0;
    
    solver_type intg;
    
    // set initial condition of ODE solver
    gvec_type x0(0.0);
    const int rhsDIM = RHS::dimension;
    for(int i=0; i<rhsDIM; i++) {
        x0[i] = seed[i];
        x0[ rhsDIM + rhsDIM*i + i  ] = 1.0; //STM identity
    }
    
    tracker.initialize(x0);
    intg.set_init_cond(x0, 0.);
    // sign of tmax determines direction
    intg.set_t_max(niter < 0 ? -1e8 : 1e8);
    intg.set_init_step(1e-5);
    intg.set_precision(_prec);
    intersect_stop<section_type> stop(_section, niter>0);
    nvis::timer _t;
    
    for (size_t counter=0 ; counter<_max_iter && (int)iterCount < std::abs(niter) ; ++counter) {
    
        if (_verbose) {
            std::ostringstream os;
            os << "poincare at " << intg.y << " / ";
            if (out.size()) {
                os << out.back().x;
            } else {
                os << seed;
            }
            os << " after " << counter << " steps and " << _t.elapsed() << " seconds\n";
            std::cerr << os.str();
        }
        
        try {
            res = intg.do_step(_rhs, step);
        } catch (MapSingularityError<state_type>& msError) {
            if(_verbose) {
                std::cerr << "Caught Singularity in RHS at " << msError.where << std::endl;
            }
            throw msError;
        } catch (...) {
            if(_verbose) {
                std::cerr << "caught exception in return_map" << std::endl;
            }
            break;
        }
        //Check for Hamiltonian violation (useful for detecting failure due to singularities)
        //  - Note: This may be really bad for this function!
        if (_throwHamError && !stepHamiltonianCheck(step) ) {
            if(_verbose) {
                std::cerr << "constant Hamiltonian violated in return_map" << std::endl;
            }
            value_type ham1 = _rhs.hamiltonian(step.t1(), step.y1());
            if(_verbose) std::cerr << " Hamiltonian evaluates to " << ham1 << " when it should be " << _Ham
                                       << " (dH = " << fabs(ham1-_Ham) << ", prec=" << _prec << ")\n";
            break;
        }
        if (_verbose) {
            std::ostringstream os;
            os << "poincare: step completed\n";
            std::cerr << os.str();
        }
        
        if (res != solver_type::OK) {
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: integration ended with result: " << (int)res << '\n';
                std::cerr << os.str();
            }
            break;
        }
        if (nvis::norm(step.y0() - step.y1())<1.0e-8) {
            std::ostringstream os;
            os << "poincare: step size undeflow\n";
            std::cerr << os.str();
            break;
        }
        
        stop(step, _verbose);
        if (_verbose) {
            std::ostringstream os;
            os << "poincare: stop check completed\n";
            std::cerr << os.str();
        }
        
        if (stop.did_stop()) {
            if (_verbose) {
                std::ostringstream os;
                os << "**** poincare: stop detected ****\n";
                std::cerr << os.str();
            }
            out.push_back(return_state());
            return_state& r = out.back();
            r.setState(stop.where());
            
            //Store the internal point
            internalOut.push_back(return_state());
            return_state& rs = internalOut.back();
            rs.setState(stop.where());
            
            tracker.mark_crossing();
            iterCount++;
            
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: intersection position computed\n";
                std::cerr << os.str();
            }
            rs.t = stop.when();
            r.t = stop.when();
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: intersection time stored\n";
                std::cerr << os.str();
            }
            
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: stop information stored\n";
                std::cerr << os.str();
            }
        }
        
        //Store info if dense output is on
        if (dense) {
            //Store every step
            internalOut.push_back(return_state());
            return_state& rs = internalOut.back();
            rs.setState(step.y1()); //Store the step
            rs.t = step.t1();
        }
        tracker(step.t1(),step.y1());
        
        if (_verbose) {
            std::ostringstream os;
            os << "poincare: Jacobian integration completed\n";
            std::cerr << os.str();
        }
    }
}


/**Return Map -> sends back vector of state_info
 *  Note: This version resets the STM at every crossing in hopes the
 *  resulting STM occurs less error over the whole orbit.
 */
template<typename RHS, typename ODESolver, typename Section>
template<typename Tracker>
void poincare_map<RHS, ODESolver, Section>::
return_map_STM(const lvec_type& seed, std::vector<return_state>& out, int niter, Tracker& tracker,
               double h, bool dense) const
{
    typedef ODESolver                solver_type;
    typename solver_type::result     res;
    typename solver_type::step       step;
    
    static const int spaceDim = RHS::dimension;
    out.clear();
    out.reserve(std::abs(niter));
    if (dense) {
        out.reserve(std::abs(niter)*1000);
    }
    size_t iterCount = 0;
    
    solver_type intg;
    // set initial condition of ODE solver
    gvec_type x0 = _section.unproject(seed);
    tracker.initialize(x0);
    intg.set_init_cond(x0, 0.);
    // sign of tmax determines direction
    intg.set_t_max(niter < 0 ? -1e8 : 1e8);
    intg.set_init_step(1e-5);
    intg.set_precision(_prec);
    intersect_stop<section_type> stop(_section, niter>0);
    nvis::timer _t;
    
    for (size_t counter=0 ; counter<_max_iter && iterCount < std::abs(niter) ; ++counter) {
    
        if (_verbose) {
            std::ostringstream os;
            os << "poincare at " << intg.y << " / ";
            if (out.size()) {
                os << out.back().x;
            } else {
                os << seed;
            }
            os << " after " << counter << " steps and " << _t.elapsed() << " seconds\n";
            std::cerr << os.str();
        }
        
        try {
            res = intg.do_step(_rhs, step);
        } catch (MapSingularityError<state_type>& msError) {
            if(_verbose) {
                std::cerr << "Caught Singularity in RHS at " << msError.where << std::endl;
            }
            throw msError;
        } catch (...) {
            if(_verbose) {
                std::cerr << "caught exception in return_map" << std::endl;
            }
            break;
        }
        //Check for Hamiltonian violation (useful for detecting failure due to singularities)
        if ( _throwHamError && !stepHamiltonianCheck(step) ) {
            if(_verbose) {
                std::cerr << "constant Hamiltonian violated in return_map" << std::endl;
            }
            value_type ham1 = _rhs.hamiltonian(step.t1(), step.y1());
            if(_verbose) std::cerr << " Hamiltonian evaluates to " << ham1 << " when it should be " << _Ham
                                       << " (dH = " << fabs(ham1-_Ham) << ", prec=" << _prec << ")\n";
            break;
        }
        if (_verbose) {
            std::ostringstream os;
            os << "poincare: step completed\n";
            std::cerr << os.str();
        }
        
        if (res != solver_type::OK) {
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: integration ended with result: " << (int)res << '\n';
                std::cerr << os.str();
            }
            break;
        }
        if (nvis::norm(step.y0() - step.y1())<1.0e-8) {
            std::ostringstream os;
            os << "poincare: step size undeflow\n";
            std::cerr << os.str();
            break;
        }
        
        //Check for a stopping condition
        stop(step, _verbose);
        if (_verbose) {
            std::ostringstream os;
            os << "poincare: stop check completed\n";
            std::cerr << os.str();
        }
        
        if (stop.did_stop()) {
            if (_verbose) {
                std::ostringstream os;
                os << "**** poincare: stop detected ****\n";
                std::cerr << os.str();
            }
            out.push_back(return_state());
            return_state& r = out.back();
            r.setState(stop.where());
            
            tracker.mark_crossing();
            iterCount++;
            
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: intersection position computed\n";
                std::cerr << os.str();
            }
            r.t = stop.when();
            if (_verbose) {
http://www.pandora.com/inactive
                std::ostringstream os;
                os << "poincare: intersection time stored\n";
                std::cerr << os.str();
            }
            
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: stop information stored\n";
                std::cerr << os.str();
            }
            //Reset the linear STM to identity
            gvec_type ny = x0; //Identity STM
            gvec_type& yCurrent = intg.y;
            for(int i=0; i<spaceDim; i++) {
                ny[i] = yCurrent[i];
            }
            intg.y = ny;
            
        }
        
        //Store info if dense output is on
        if (dense) {
            //Store every step
            out.push_back(return_state());
            return_state& r = out.back();
            r.setState(step.y1()); //Store the step
            r.t = step.t1();
        }
        tracker(step.t1(),step.y1());
        
        if (_verbose) {
            std::ostringstream os;
            os << "poincare: Jacobian integration completed\n";
            std::cerr << os.str();
        }
    }
}

template<typename RHS, typename ODESolver, typename Section>
void poincare_map<RHS, ODESolver, Section>::
return_monodromy(const lvec_type& seed, const double& dtMax,
                 const std::vector<double>& tVals, std::vector<return_state>& out
                ) const
{
    typedef typename return_state::mat_type  MatType;
    typename solver_type::step       step;
    const static int N = RHS::dimension;
    
    out.clear();
    
    //Compute the time steps
    int niter = (int) tVals.size();
    std::vector<double> times(1,0.0);
    std::vector<int> crossingIdx;
    double tTotal = 0.0;
    for(int i=0; i<niter; i++) {
        int numDT = floor( tVals[i] / dtMax );
        //Add the partial time steps
        if(numDT>0) {
            for(int j=1; j<=numDT; j++) {
                times.push_back( tTotal + j*dtMax );
            }
        }
        tTotal = tVals[i];
        //Add the final time
        times.push_back( tTotal );
        crossingIdx.push_back( (int)times.size() - 1 );
    }
    
    //Run the map for each time step, reseting the STM to identity each time
    int numProp = crossingIdx.back();
    std::vector<return_state> dtStates;
    gvec_type y = _section.unproject(seed);
    for(int i=0; i<numProp; i++) {
        double dt = times[i+1] - times[i];
        //Get state, time for prop
        dtStates.push_back( integrate_state(y,dt) );
        //Set the state based on results
        y = dtStates.back().getState();
        //Reset STM
        for (int k=N ; k<((N+1)*N); k++) {
            y[k] = 0.0;
        }
        for (int k=0 ; k<N ; ++k) {
            y[N*(k+1)+k] = 1.0;
        }
    }
    //Compute the STM's for each crossing by multiplication
    MatType phi;
    phi.identity();
    std::vector<int>::iterator cid = crossingIdx.begin();
    int cross = 0;
    for(int i=0; i<numProp; i++) {
        MatType& J = dtStates[i].J;
        phi = J * phi; //STM(t2,0) = STM(t2,t1)*STM(t1,0)
        //Store information if on crossing index
        if(i == (*cid)) {
            out.push_back(dtStates[i]);
            out[cross].t = tVals[cross];
            out[cross].J = phi; //Accumulated STM
            //Increment counters
            ++cid;
            ++cross;
        }
    }
    
}



template<typename RHS, typename ODESolver, typename Section>
void poincare_map<RHS, ODESolver, Section>::
_flow_map(const lvec_type& seed, return_state& out, double t) const
{
    typedef ODESolver                         solver_type;
    typename solver_type::result              res;
    typename solver_type::step                step;
    
    solver_type intg;
    // set initial condition of ODE solver
    intg.set_init_cond(_section.unproject(seed), 0.);
    // sign of tmax determines direction
    intg.set_t_max(t>0 ? 1.0e+8 : -1.0e+8);
    intg.set_init_step(1e-5);
    intg.set_precision(_prec);
    nvis::timer _t;
    
    size_t counter = 0;
    while (counter<_max_iter) {
        ++counter;
        if (_verbose) {
            std::ostringstream os;
            os << "flow_map at " << intg.y << " / ";
            os << " after " << counter << " steps and " << _t.elapsed() << " seconds\n";
            std::cerr << os.str();
        }
        
        try {
            res = intg.do_step(_rhs, step);
        } catch (MapSingularityError<state_type>& msError) {
            if(_verbose) {
                std::cerr << "Caught Singularity in RHS at " << msError.where << std::endl;
            }
            throw msError;
        } catch (...) {
            if(_verbose) {
                std::cerr << "caught exception in return_map" << std::endl;
            }
            IntegError fInteg(std::string("Integrator Error"));
            fInteg.where = seed;
            fInteg.y = step.y0();
            fInteg.t = step.t0();
            throw fInteg;
            break;
        }
        //Check for Hamiltonian violation (useful for detecting failure due to singularities)
        if ( _throwHamError && !stepHamiltonianCheck(step) ) {
            if(_verbose) {
                std::cerr << "constant Hamiltonian violated in _flow_map" << std::endl;
            }
            value_type ham1 = _rhs.hamiltonian(step.t1(), step.y1());
            if(_verbose) std::cerr << " Hamiltonian evaluates to " << ham1 << " when it should be " << _Ham
                                       << " (dH = " << fabs(ham1-_Ham) << ", prec=" << _prec << ")\n";
            IntegError fInteg(std::string("Hamiltonian Violation Detected"));
            fInteg.where = seed;
            fInteg.y = step.y0();
            fInteg.t = step.t0();
            throw fInteg;
            break;
        }
        if (_verbose) {
            std::ostringstream os;
            os << "flow_map: step completed\n";
            std::cerr << os.str();
        }
        
        if (res != solver_type::OK) {
            if (_verbose) {
                std::ostringstream os;
                os << "flow_map: integration ended with result: " << (int)res << '\n';
                std::cerr << os.str();
            }
            break;
        }
        if (nvis::norm(step.y0() - step.y1())<1.0e-8) {
            std::ostringstream os;
            os << "flow_map: step size undeflow\n"; //Throw error?
            std::cerr << os.str();
            break;
        }
        
        if ((t>0 && step.t1()>=t) || (t<0 && step.t1()<=t)) {
            //std::pair<lvec_type, lmat_type> tmp = _section.project(step.y(t));
            //out.x = tmp.first;
            //out.J = tmp.second;
            out.setState(step.y(t));
            out.t = t;
            out.nsteps = counter;
            
            if (_verbose) {
                std::ostringstream os;
                os << "flow_map: integration completed\n";
                std::cerr << os.str();
            }
            return;
        }
    }
}

template<typename RHS, typename ODESolver, typename Section>
void poincare_map<RHS, ODESolver, Section>::
_integrate_state(const gvec_type& y0, return_state& out, double t) const
{
    typedef ODESolver                         solver_type;
    typename solver_type::result         res;
    typename solver_type::step                 step;
    
    solver_type intg;
    // set initial condition of ODE solver
    intg.set_init_cond(y0, 0.);
    // sign of tmax determines direction
    intg.set_t_max(t>0 ? 1.0e+8 : -1.0e+8);
    intg.set_init_step(1e-5);
    intg.set_precision(_prec);
    nvis::timer _t;
    
    size_t counter = 0;
    while (counter<_max_iter) {
        ++counter;
        if (_verbose) {
            std::ostringstream os;
            os << "integrate_state at " << intg.y << " / ";
            os << " after " << counter << " steps and " << _t.elapsed() << " seconds\n";
            std::cerr << os.str();
        }
        
        try {
            res = intg.do_step(_rhs, step);
        } catch (MapSingularityError<state_type>& msError) {
            if(_verbose) {
                std::cerr << "Caught Singularity in RHS at " << msError.where << std::endl;
            }
            throw msError;
        } catch (...) {
            if(_verbose) {
                std::cerr << "caught exception in integrate_state" << std::endl;
            }
            IntegError fInteg(std::string("Integrator Error"));
            fInteg.where = _section.project(y0).first;
            fInteg.y = step.y0();
            fInteg.t = step.t0();
            throw fInteg;
            break;
        }
        //Check for Hamiltonian violation (useful for detecting failure due to singularities)
        if (_throwHamError && !stepHamiltonianCheck(step) ) {
            if(_verbose) {
                std::cerr << "constant Hamiltonian violated in integrate_state" << std::endl;
            }
            value_type ham1 = _rhs.hamiltonian(step.t1(), step.y1());
            if(_verbose) std::cerr << " Hamiltonian evaluates to " << ham1 << " when it should be " << _Ham
                                       << " (dH = " << fabs(ham1-_Ham) << ", prec=" << _prec << ")\n";
            IntegError fInteg(std::string("Hamiltonian Violation Detected"));
            fInteg.where = _section.project(y0).first;
            fInteg.y = step.y0();
            fInteg.t = step.t0();
            throw fInteg;
            break;
        }
        if (_verbose) {
            std::ostringstream os;
            os << "integrate_state: step completed\n";
            std::cerr << os.str();
        }
        
        if (res != solver_type::OK) {
            if (_verbose) {
                std::ostringstream os;
                os << "integrate_state: integration ended with result: " << (int)res << '\n';
                std::cerr << os.str();
            }
            break;
        }
        if (nvis::norm(step.y0() - step.y1())<1.0e-8) {
            std::ostringstream os;
            os << "integrate_state: step size undeflow\n"; //Throw Error?
            std::cerr << os.str();
            break;
        }
        
        if ((t>0 && step.t1()>=t) || (t<0 && step.t1()<=t)) {
            //std::pair<lvec_type, lmat_type> tmp = _section.project(step.y(t));
            //out.x = tmp.first;
            //out.J = tmp.second;
            out.setState(step.y(t));
            out.t = t;
            out.nsteps = (unsigned int)counter;
            
            if (_verbose) {
                std::ostringstream os;
                os << "integrate_state: integration completed\n";
                std::cerr << os.str();
            }
            return;
        }
    }
}


///Compute Hits to a Region ends back vector of return_map_info
template<typename RHS, typename ODESolver, typename Section>
template<typename Tracker>
void poincare_map<RHS, ODESolver, Section>::
compute_hits(const lvec_type& seed, std::vector<return_type>& out, int nHits, int maxIters, const lbox_type& regionBB, Tracker& tracker, double h) const
{
    typedef ODESolver                solver_type;
    typename solver_type::result     res;
    typename solver_type::step       step;
    
    //Declare some sizes
    int hitsFound = 0;
    int numStops = 0;
    out.clear();
    out.reserve(std::abs(nHits));
    
    solver_type intg;
    // set initial condition of ODE solver
    gvec_type x0 = _section.unproject(seed);
    tracker.initialize(x0);
    intg.set_init_cond(x0, 0.);
    // sign of tmax determines direction
    intg.set_t_max(nHits < 0 ? -1e8 : 1e8);
    intg.set_init_step(1e-5);
    intg.set_precision(_prec);
    intersect_stop<section_type> stop(_section, nHits>0);
    nvis::timer _t;
    
    // Loop through hits & iters until maxIters or filled hits
    for (size_t counter=0 ; counter < _max_iter && numStops < std::abs(maxIters) && hitsFound < std::abs(nHits) ; ++counter) {
    
        if (_verbose) {
            std::ostringstream os;
            os << "poincare at " << intg.y << " / ";
            if (out.size()) {
                os << out.back().x;
            } else {
                os << seed;
            }
            os << " after " << counter << " steps and " << _t.elapsed() << " seconds\n";
            std::cerr << os.str();
        }
        
        try {
            res = intg.do_step(_rhs, step);
        } catch (MapSingularityError<state_type>& msError) {
            if(_verbose) {
                std::cerr << "Caught Singularity in RHS at " << msError.where << std::endl;
            }
            throw msError;
        } catch (...) {
            if(_verbose) {
                std::cerr << "caught exception in return_map" << std::endl;
            }
            break;
        }
        //Check for Hamiltonian violation (useful for detecting failure due to singularities)
        if (_throwHamError && !stepHamiltonianCheck(step) ) {
            if(_verbose) {
                std::cerr << "constant Hamiltonian violated in return_map" << std::endl;
            }
            value_type ham1 = _rhs.hamiltonian(step.t1(), step.y1());
            if(_verbose) std::cerr << " Hamiltonian evaluates to " << ham1 << " when it should be " << _Ham
                                       << " (dH = " << fabs(ham1-_Ham) << ", prec=" << _prec << ")\n";
            break;
        }
        
        if (_verbose) {
            std::ostringstream os;
            os << "poincare: step completed\n";
            std::cerr << os.str();
        }
        
        if (res != solver_type::OK) {
            if (_verbose) {
                std::ostringstream os;
                os << "poincare: integration ended with result: " << (int)res << '\n';
                std::cerr << os.str();
            }
            break;
        }
        if (nvis::norm(step.y0() - step.y1())<1.0e-8) {
            std::ostringstream os;
            os << "poincare: step size undeflow\n";
            std::cerr << os.str();
            break;
        }
        
        stop(step, _verbose);
        if (_verbose) {
            std::ostringstream os;
            os << "poincare: stop check completed\n";
            std::cerr << os.str();
        }
        
        // Section crossing detected
        if (stop.did_stop()) {
            if (_verbose) {
                std::ostringstream os;
                os << "**** poincare: stop detected ****\n";
                std::cerr << os.str();
            }
            //The stop information
            std::pair<lvec_type, lmat_type> tmp = _section.project(stop.where());
            numStops++;
            //We can only store if it returns to the specified region
            if( regionBB.inside( tmp.first ) ) {
                //Storing
                hitsFound++;
                out.push_back(return_type());
                return_type& r = out.back();
                r.x = tmp.first;
                r.J = tmp.second;
                r.delta_theta = tracker(stop.when(),stop.where());
                tracker.mark_crossing();
                if (_verbose) {
                    std::ostringstream os;
                    os << "poincare: intersection position computed\n";
                    std::cerr << os.str();
                }
                r.t = stop.when();
                if (_verbose) {
                    std::ostringstream os;
                    os << "poincare: intersection time stored\n";
                    std::cerr << os.str();
                }
                
                if (_verbose) {
                    std::ostringstream os;
                    os << "poincare: stop information stored\n";
                    std::cerr << os.str();
                }
            }
        }
        
        if (_verbose) {
            std::ostringstream os;
            os << "poincare: Jacobian integration completed\n";
            std::cerr << os.str();
        }
        
    }
}

} // namespace xavier

#endif















