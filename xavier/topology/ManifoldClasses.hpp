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


/* Manifold Classes for computing manifolds on the map
 *  - Various classes for use with England method (and modified method)
 * Author: Wayne Schlei and Xavier Tricoche (Purdue University)
 */
#ifndef MANIFOLD_CLASSES_HPP
#define MANIFOLD_CLASSES_HPP

#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <math/fixed_vector.hpp>
#include <math/bezier_spline.hpp>
#include <math/bounding_box.hpp>
#include <maps/metric.hpp>
#include <maps/fixpoints.hpp>
#include <maps/map_analysis.hpp>
#include <topology/EdgeRotationFailure.hpp>

namespace xavier {
//struct fixpoint;

/// Storing all fixed points associated with a particular periodic orbit
struct fp_chain {
    typedef xavier::fixpoint        fp_type;
    
    fp_chain(const std::vector<fp_type>& _fps, const metric_type& _metric)
        : fixed_points(_fps), _metric(_metric) {}
        
    template<typename Iterator>
    fp_chain(Iterator begin, Iterator end, const metric_type& _metric)
        : fixed_points(begin, end), _metric(_metric) {}
        
    unsigned int period() const
    {
        return (unsigned int) fixed_points.size();
    }
    
    bool saddle() const
    {
        return fixed_points.size() && fixed_points[0].saddle;
    }
    
    const fp_type& operator[](unsigned int i) const
    {
        if (!period()) {
            throw std::runtime_error("empty chain");
        }
        
        unsigned int j = (i % period());
        return fixed_points[j];
    }
    
    unsigned int closest(const nvis::vec2& x) const
    {
        std::vector<double> dist(period());
        for (int i = 0 ; i < (int) dist.size() ; ++i) {
            dist[i] = _metric.distance(x, fixed_points[i].pos);
        }
        return std::distance(dist.begin(), std::min_element(dist.begin(), dist.end()));
    }
    unsigned int closest(const fp_type& fp) const
    {
        return closest(fp.pos);
    }
    
    double distance(const nvis::vec2& x) const
    {
        return _metric.distance(x, fixed_points[closest(x)].pos);
    }
    double distance(const fp_type& fp) const
    {
        return distance(fp.pos);
    }
    
    double distance(const fp_chain& c) const
    {
        if (c.period() != period()) {
            return std::numeric_limits<double>::max();
        }
        
        std::vector<double> dist(period());
        for (int i = 0 ; i < (int)dist.size() ; ++i) {
            dist[i] = distance(c[i]);
        }
        return *std::min_element(dist.begin(), dist.end());
    }
    
    size_t size() const
    {
        return fixed_points.size();
    }
    
    std::vector<xavier::fixpoint>        fixed_points;
    metric_type                          _metric;
};

struct Neighbor {
    Neighbor() : chain_id(-1), p_id(-1) {}
    Neighbor(int cid, int pid) : chain_id(cid), p_id(pid) {}
    Neighbor(const Neighbor& n) : chain_id(n.chain_id), p_id(n.p_id) {}
    int chain_id;
    int p_id;
};

struct Connection {
    Connection() : neigh(), d(std::numeric_limits<double>::max()) {}
    Connection(const Neighbor& n, double dist) : neigh(n), d(dist) {}
    Neighbor neigh;
    double d;
};

struct Lt_Connection {
    bool operator()(const Connection& c0, const Connection& c1) const
    {
        return (c0.d < c1.d);
    }
};

typedef std::vector<nvis::vec2>            manifold_type;
typedef std::pair<int, int>                fp_index_type;



/// Data Structure for storing manifold infomration
/// - Crude version originally used for testing, but is now deprecated.
struct Separatrix {
    typedef EdgeRotationFailure<nvis::vec2> MapDiscont;
    fp_index_type         start, end;
    manifold_type         manifold;
    std::set<int>         breakIDs;
    std::list<MapDiscont> discontList;
    int                   lastSeed;
    bool                  forward;
    double                length;
};



/// Enum for manifold type
enum Perturbation {
    UNSTABLE_PLUS,
    UNSTABLE_MINUS,
    STABLE_PLUS,
    STABLE_MINUS
};


//Later:  Need to implement something more like this for Separatrix in
//order to construct designs!
/*template <typename VEC>
struct ManifoldPoint {
public:
  typedef std::pair<int,int> SeedSegID;
  VEC       fixedPointSource;
  VEC       point;
  SeedSegID id; //(-1,0) indicates starting segment
  double    tau;
  int       period;
}

template <typename VEC>
class Manifold {
public:
  typedef std::vector<ManifoldPoint<VEC> > History;
  History         points;
  Perturbation    pert;
}
typedef Manifold<nvis::vec2>   manifold_type; */

/// Data structure for 1D Global Manifolds by Continuation settings
struct ManifoldSettings {
    //Defaults
    ManifoldSettings() : eps(1e-7), manualStep(false),
        sdelta_max(1.e-5), max_si(1e10),
        alpha_max(0.3),  liberal_alpha_max(M_PI/2.),
        delta_alpha_max(1.e-3), delta_min(1.0e-4),
        max_arc_length(30), enableMaxArclength(true),
        maxSeedArclength(0.5), maxTreeDepth(8),
        maxNumSeg(8000), maxNumSeedSeg(500), enableMaxSeg(false),
        delta_tau_min(5.e-6), delta_max(0.3) {}
        
    //Tune these parameters per system
    /// eps = used to find initial step
    double eps;
    /// Option of whether or not to use a manual initial step or a function of lambda_max
    bool manualStep;
    /// sdelta_max = Maximum step distance for first seeding point from fixed-points
    double sdelta_max; //Maximum initial allowed step
    /// Maximum allowed stability index magnitude (useful for cropping computation)
    double max_si; //Originally off at 1e10
    /// alpha_max = Maximum angle between 3 successive points (rad)
    double alpha_max;
    /// liberal_alpha_max = A broader range version
    double liberal_alpha_max;
    /// delta_alpha_max = Control sampling rate of segments (delta*alpha=0 means straight line)
    double delta_alpha_max;
    /// delta_min = Minimum distance between manifold points
    double delta_min;
    /// max_arc_length = Maximum arc length to be achieved by the method (stopping condition)
    double max_arc_length;
    /// Enable the maximum arc length stopping condition
    bool enableMaxArclength;
    /// Maximum seeding arc-length (arc that the seeding segments make)
    double maxSeedArclength;
    /// Maximum depth level to achieve in invariant manifold tree
    int maxTreeDepth;
    /// Maximum number of segments
    int maxNumSeg;
    /// Maximum number of seeding segments
    int maxNumSeedSeg;
    /// Enable the maximum number of segments and seeding segments constraints
    bool enableMaxSeg;
    /// Minimum change in tau allowed between steps
    double delta_tau_min;
    /// Maximum distance between downstream points (Disable if testing changes to section sep)
    double delta_max;
    /// Bounds for viable segments (display and counting arc length)
    nvis::bbox2 bounds;
    
    ///Write settings to file
    bool write(const char* filename) const
    {
        FILE* f = fopen(filename, "w");
        if(!f) {
            std::cerr << " Bad filename for WRITING ManifoldSettings...\n";
            return false;
        }
        
        //Header
        fprintf(f,"# Invariant Manifold Advection Settings File for STHManifold algorithm\n\n");
        
        fprintf(f,"-----------------------------------------\n");
        fprintf(f," PARAMETERS:\n");
        fprintf(f,"-----------------------------------------\n");
        fprintf(f,"EPS= %.15f\n",eps);
        fprintf(f,"STEP_MAX= %.15f\n",sdelta_max);
        fprintf(f,"MANUAL_STEP= %d\n",(manualStep)? 1 : 0);
        fprintf(f,"MAX_SI= %.15f\n",max_si);
        fprintf(f,"ENABLE_MAX_ARC_LENGTH= %d\n",(enableMaxArclength)? 1 : 0);
        fprintf(f,"ENABLE_MAX_SEG_COUNT= %d\n",(enableMaxSeg)? 1 : 0);
        fprintf(f,"BOUNDS= %.15f %.15f %.15f %.15f\n",
                bounds.min()[0],bounds.min()[1],
                bounds.max()[0],bounds.max()[1]);
                
                
        fprintf(f,"\n\n");
        fprintf(f,"-----------------------------------------\n");
        fprintf(f," HEURISTICS:\n");
        fprintf(f,"-----------------------------------------\n");
        fprintf(f,"ALPHA_MAX= %.15f\n",alpha_max);
        fprintf(f,"LIBERAL_ALPHA= %.15f\n",liberal_alpha_max);
        fprintf(f,"DELTA_MIN= %.15f\n",delta_min);
        fprintf(f,"DELTA_MAX= %.15f\n",delta_max);
        fprintf(f,"DELTA_ALPHA_MAX= %.15f\n",delta_alpha_max);
        fprintf(f,"DTAU_MIN= %.15f\n",delta_tau_min);
        fprintf(f,"MAX_ARC_LENGTH= %.15f\n",max_arc_length);
        fprintf(f,"MAX_TREE_DEPTH= %d\n",maxTreeDepth);
        fprintf(f,"MAX_SEED_ARC_LENGTH= %.15f\n",maxSeedArclength);
        fprintf(f,"MAX_NUM_SEED_SEGS= %d\n",maxNumSeedSeg);
        fprintf(f,"MAX_NUM_SEGS= %d\n",maxNumSeg);
        fprintf(f,"\n");
        
        fclose(f);
        return true;
    }
    
    /// Read from file (NOTE: Fixed-format read!)
    bool read(const char* filename)
    {
        FILE* f = fopen(filename, "r");
        if(!f) {
            std::cerr << " Bad filename for READING ManifoldSettings...\n";
            return false;
        }
        //Ideally, use regex, but no time!
        
        char buffer[80];
        int tempBool = 0;
        //Headers
        for(int i=0; i<5; i++) {
            fgets(buffer, 80, f);
        }
        fscanf(f, "%s %lf", buffer, &eps);
        fscanf(f, "%s %lf", buffer, &sdelta_max);
        fscanf(f, "%s %d", buffer, &tempBool);
        manualStep = (tempBool==1);
        fscanf(f, "%s %lf", buffer, &max_si);
        fscanf(f, "%s %d", buffer, &tempBool);
        enableMaxArclength = (tempBool==1);
        fscanf(f, "%s %d", buffer, &tempBool);
        enableMaxSeg = (tempBool==1);
        fscanf( f, "%s %lf %lf %lf %lf",buffer,
                &(bounds.min()[0]),&(bounds.min()[1]),
                &(bounds.max()[0]),&(bounds.max()[1])  );
        //Heuristics
        for(int i=0; i<3; i++) {
            fscanf(f,"%s",buffer);    //Read separator lines
        }
        fscanf(f, "%s %lf", buffer, &alpha_max);
        fscanf(f, "%s %lf", buffer, &liberal_alpha_max);
        fscanf(f, "%s %lf", buffer, &delta_min);
        fscanf(f, "%s %lf", buffer, &delta_max);
        fscanf(f, "%s %lf", buffer, &delta_alpha_max);
        fscanf(f, "%s %lf", buffer, &delta_tau_min);
        fscanf(f, "%s %lf", buffer, &max_arc_length);
        fscanf(f, "%s %d", buffer, &maxTreeDepth);
        fscanf(f, "%s %lf", buffer, &maxSeedArclength);
        fscanf(f, "%s %d", buffer, &maxNumSeedSeg);
        fscanf(f, "%s %d", buffer, &maxNumSeg);
        
        fclose(f);
        return true;
    }
    
};

struct ManifoldProgress {
    ManifoldProgress() : length(0), seedingLength(0) {}
    
    double length;
    double seedingLength;
    double tau;
    unsigned int segment;
};


///Structure indicating manifold stopping criteria (for EnglandManifold)
/// NOTE:  This does NOT apply to STHManifold advection (i.e., curve refinement)
template<typename SECTION>
class ManifoldStopper {
    typedef typename SECTION::lbox_type         bbox_type;
    typedef typename SECTION::lvec_type         lvec_type;
    
public :
    ManifoldStopper(const std::vector<fp_chain>& all_chains, const xavier::metric_type& metric_,
                    const bbox_type bounds, const ManifoldSettings& settings)
        : _chains(all_chains), _metric(metric_), _bounds(bounds), _eps(settings.eps), _maxlength(settings.max_arc_length),
          _left_dist(settings.sdelta_max), _left(false), _periodOne(false), _homoclinicDist(0.0)
    {
        //The old default _left_distance
        //_left_dist = 0.3*min_dist(all_chains, _metric); //Large value for a non-periodic system
    }
    ManifoldStopper(const ManifoldStopper& other) :
        _chains(other._chains),
        _metric(other._metric),
        _bounds(other._bounds),
        _eps(other._eps),
        _maxlength(other._maxlength),
        _left_dist(other._left_dist),
        _left(other._left),
        _periodOne(other._periodOne),
        _homoclinicDist(other._homoclinicDist)
    {}
    
    ///Initiate the propagation
    void start(const fp_index_type fpi)
    {
        _fpi = fpi;
        _periodOne = (_chains[fpi.first].period() == 1 ? true : false);
        _homoclinicDist = 1000.0;
        _start = _chains[fpi.first][fpi.second].pos;
        if(!_periodOne) {
            int p = _chains[fpi.first].period();
            for(int i=0; i<p; i++) {
                //Compute minimum distance
                double d = _metric.distance(_chains[fpi.first][i].pos, _start);
                if (d<_homoclinicDist) {
                    _homoclinicDist = d;
                }
            }
            //Use the half distance to indicate approach
        }
        _end_chain = _end_fp = -1;
        _min_dist_id = fp_index_type(-1, -1);
        _min_dist = std::numeric_limits<double>::max();
        _approaching = false;
        _counter = 0;
        _last = _start;
        _length = 0;
        _left = false;
    }
    
    ///Checking function -> Evaluates stopping criteria
    bool operator()(const lvec_type& x)
    {
        ++_counter;
        /*std::cerr << "In STOP: Step " << _counter << " \n";
        std::cerr << "             current value x = " << x << "\n";
        std::cerr << "             last point = " << _last << " \n";
        std::cerr << "             current length = " << _length << "\n";
        std::cerr << "             disp = " << _metric.displacement(_last,x) << " with norm = " << _metric.distance(_last,x) << "\n";
        std::cerr << "             new length = " << _length + _metric.distance(_last,x) << "\n";
        */
        _length += _metric.distance(_last, x);
        //Arc has not sufficiently departed fixed-point
        if (!_left && _metric.distance(x, _start) < _left_dist) {
            if (_counter > 15) {
                std::cerr << "unable to leave start. giving up.\n";
                throw std::runtime_error("wrong type suspected");
            }
            return false;
        } else if (!_left) {
            std::cerr << "left start\n";
            _left = true;
        }
        
        //Stop if left visible domain
        if (!_bounds.inside(x)) {
            std::cerr << "   segment left bounds\n";
            return true;
        }
        
        //Stop on maxArcLength (optional)
        if (_length > _maxlength) {
            std::cerr << "   segment has completed length criterion with length [" << _length
                      << "] and maxArcLength [" << _maxlength << "]\n";
            return true;
        }
        
        //Set new _last after arcLength check
        _last = x;
        
        //Check if we are approaching a fixed_point (but different branch for p=1 and p>1)
        if (!_periodOne) {
            //Check to see if we are nearing another fixed_point in chain (p>1)
            std::map<double, fp_index_type> dist;
            for (int i=0 ; i<(int)_chains.size() ; ++i) {
                for (int j=0 ; j<(int)_chains[i].size() ; ++j) {
                    double d = nvis::norm(x-_chains[i][j].pos);
                    dist[d] = fp_index_type(i, j);
                }
            }
            std::map<double, fp_index_type>::iterator it = dist.begin();
            double d = it->first;
            fp_index_type fp_index = it->second;
            //lvec_type pos = _chains[it->second.first][it->second.second].pos;
            if (!_approaching &&
                    (fp_index.first != _fpi.first || fp_index.second != _fpi.second)) {
                if (d < 0.5*_metric.distance(x, _start)) {
                    _approaching = true;
                }
                return false;
            }
            if (_approaching && d < 0.002*_length) {
                _end_chain = dist.begin()->second.first;
                _end_fp = dist.begin()->second.second;
                std::cerr << "   segment approaching Period-n fixed point\n";
                return true;
            }
        } else {
            //Period 1 saddle
            double p1dist = _metric.distance(x, _start);
            if (!_left && p1dist > 50*_eps) {
                _left = true;
            }
            if (!_left) {
                return false;
            }
            if (_left && !_approaching && p1dist < 10*_eps) {
                _approaching = true;
                _nb_backtracked = 0;
            }
            if (_approaching) {
                if (p1dist<_min_dist) {
                    _min_dist = p1dist;
                    _nb_backtracked = 0;
                } else {
                    ++_nb_backtracked;
                }
                if (_nb_backtracked > 3) { //Note: backtracking doesn't really apply to STHManifold
                    std::cerr << "too many backtrack steps. stop.\n";
                    return true;
                }
            }
            if (_approaching && _metric.distance(x, _start) < _eps) {
                std::cerr << "   segment has returned to period-1 fixed point\n";
                return true;
            }
            return false;
        }
        
        //Could we match stable/unstable arcs?
        // - better stopping or stablization of tangles?
        
        return false;
    }
    
    ///Check to see if we are nearing another fixed_point in the same chain (p>1) to stop manifold
    bool homoclinicTest(const lvec_type& x)
    {
        std::map<double, fp_index_type> dist;
        ///Checking ALL chains for an approach...
        /*for (int i=0 ; i<(int)_chains.size() ; ++i) {
            for (int j=0 ; j<(int)_chains[i].size() ; ++j) {
                double d = nvis::norm(x-_chains[i][j].pos);
                dist[d] = fp_index_type(i, j);
            }
        }*/
        ///Checking fixed points within the current chain(orbit)
        int i=_fpi.first; //Orbit ID
        for(int j=0; j<(int)_chains[i].size(); j++) {
            if(j==_fpi.second) {
                continue;
            }
            double d = _metric.distance(x,_chains[i][j].pos);
            dist[d] = fp_index_type(i,j);
        }
        
        std::map<double, fp_index_type>::iterator it = dist.begin();
        double d = it->first;
        //fp_index_type fp_index = it->second;
        //lvec_type pos = _chains[it->second.first][it->second.second].pos;
        if (!_approaching ) {
            if (d < _homoclinicDist) {
                _approaching = true;
            }
            return false;
        }
        //Use the half distance to indicate approach
        if (_approaching && d < _homoclinicDist/2.0) {
            _end_chain = dist.begin()->second.first;
            _end_fp = dist.begin()->second.second;
            std::cerr << " ManCurve: Manifold approaching Homoclinic connection.  Stopping...\n";
            return true;
        }
        
        return false;
    }
    
    
    
    fp_index_type connection() const
    {
        return std::make_pair(_end_chain, _end_fp);
    }
    
private :
    const std::vector<fp_chain>&         _chains;
    xavier::metric_type                  _metric;
    bbox_type                            _bounds;
    double                               _eps;
    double                               _maxlength;
    double                               _left_dist;
    lvec_type                            _start;
    fp_index_type                        _fpi;
    int                                  _end_chain;
    int                                  _end_fp;
    fp_index_type                        _min_dist_id;
    double                               _min_dist;
    bool                                 _approaching;
    int                                  _counter;
    lvec_type                            _last;
    double                               _length;
    bool                                 _left;
    bool                                 _periodOne;
    double                               _homoclinicDist;
    int                                  _nb_backtracked;
};



} //end xavier


#endif