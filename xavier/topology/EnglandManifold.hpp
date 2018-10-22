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


//////////////////////////////////////////////////////////////////////////////
// EnglandManifold.hpp
// Author:  Xavier Tricoche & Wayne Schlei
// Date:  2/21/2013
// Purpose:  Utilize the ManBVP algorithm for computing manifolds on a
// Poincare section (from Computing 1D Global Manifolds By Continuation by
// England et al.).
//
// Update notice:
// The original version of this has been modified to make algorithm
// parameters more tractable for CR3BP.  Also, the stopping criteria was
// modified to include leaving a region of interest.
//
// Future Work:
// 1)This propagates manifolds on a 2D surface of section; for work with higher
// dimensional sections, be sure to conform to a template that specifies dimensions.
// 2)This stops when you leave a region of interest, but it might be nice to
// check if a manifold will come back (or add ability to add segements that return).
///////////////////////////////////////////////////////////////////////////////

#ifndef ENGLAND_MANIFOLD_HPP
#define ENGLAND_MANIFOLD_HPP

#include <vector>
#include <list>
#include <math/fixed_vector.hpp>
#include <math/bezier_spline.hpp>
#include <math/bounding_box.hpp>
#include <maps/metric.hpp>
#include <maps/fixpoints.hpp>
#include <maps/map_analysis.hpp>
#include <topology/ManifoldClasses.hpp>

using namespace xavier;

namespace topology {

//////////////////////////////////////////////////////////////////////////////
//
// England et al.'s 1D manifold construction technique
// (SIAM J. Appl. Dynamical Systems, 2005)
//
//////////////////////////////////////////////////////////////////////////////


/// Boundary Value Problem solution step
template<typename MAP,typename VEC>
inline VEC
BVP_Step(const MAP& map, int period, const metric_type& _metric,
         const VEC& start, const VEC& end,                        // seeding segment
         const VEC& prev, const VEC& cur,                         // end segment
         const ManifoldSettings& settings,                        // quality control
         //TurningPoint& tp,
         double& tau, const map_analysis_param& param)            // where are we along seeding segment
{
    double delta_min = settings.delta_min;
    double alpha_max = settings.alpha_max;
    double delta_alpha_max = settings.delta_alpha_max;
    
    double delta, alpha;
    double t = 1;
    VEC x, y, span;
    span = _metric.displacement(start, end);
    
    bool verbose = param.verbose;
    //verbose = false;
    
    // VEC q = (1 - tau) * start + tau * end;
    VEC q = _metric.modulo(start + tau*span);
    
    if (verbose)
        std::cerr << "BVP called at " << q << " (" << tau << ") on segment "
                  << start << " - " << end << ", prev = " << prev
                  << ", cur = " << cur << '\n';
                  
    int iterationCount = 0;
    while (t > tau) {
        x = _metric.modulo(start + t*span);
        if (verbose) {
            std::cerr << "currently trying at " << x << ". ";
        }
        bool integrationOK = true;
        VEC _y;
        try {
            _y = map.map(x, period);
        } catch(...) {
            std::cerr << " BVP_Step:  Map integration failure. Skipping a step.\n";
            integrationOK = false;
        }
        iterationCount++;
        if (!integrationOK) {
            if (verbose) {
                std::cerr << "FAILED\n";
            }
            t = 0.5*(tau+t);
        }
        if (verbose) {
            std::cerr << "map(" << x << ") = " << _y;
        }
        y = _metric.modulo(_y);
        if (verbose) {
            std::cerr << " (=" << y << ")\n";
            std::cerr << "\tat " << y << " for t=" << t << '\n';
        }
        alpha = _metric.angle(prev, cur, y);
        delta = _metric.distance(cur, y);
        
        //Run a turning point check for the computed alpha if no turning point detected
        //if (!tp.detected) tp.checkAlpha(alpha);
        
        if (verbose) {
            std::cerr << "delta = " << delta << " (" << delta_min << "), alpha = " << alpha
                      << ", delta*alpha = " << delta* alpha << " (" << delta_alpha_max << "): ";
        }
        if (fabs(delta) < delta_min ||
                (fabs(alpha) < alpha_max && fabs(delta*alpha) < delta_alpha_max)) {
            if (verbose) {
                std::cerr << "PASSED\n";
            }
            break;
        } else {
            if (verbose) {
                std::cerr << "FAILED\n";
            }
            //t = (tau + t)/2.0; //Average the current t and tau for a better solution
            t = (9.0*tau + t)/10.0; //Attempt to accelerate convergence
            //Small steps are usually required in stiff problems (like CR3BP), so use weighted average to last guess
        }
        if (verbose) {
            std::cerr << "delta = " << delta << ", alpha = " << alpha << '\n';
        }
    }
    
    std::cerr << " In BVP_Step:  This solution took " << iterationCount << " map integrations.\n";
    tau = t;
    return y;
}

/// The Man BVP algorithm for computing 1D manifolds with continuation
template<typename MAP, typename STOP>
inline void ManBVP(manifold_type& manifold, std::vector<int>& breakIDs, ManifoldProgress& progress,
                   const MAP& map, STOP& stop, int period, const metric_type& _metric,
                   const xavier::fixpoint& saddle, const Perturbation& mType,
                   const ManifoldSettings& settings, const map_analysis_param& param)
{
    typedef typename MAP::section_type::lvec_type lvec_type;
    unsigned int i0, i1;
    double tau;
    double length;
    //Propagation direction (stable or unstable)
    bool fwd = (mType == UNSTABLE_PLUS || mType == UNSTABLE_MINUS) ? true : false;
    //Plus or minus perturbation
    bool p_or_m = (mType == STABLE_PLUS || mType == UNSTABLE_PLUS) ? true : false;
    
    //Find the initial step
    if (progress.length == 0) {
        if (param.verbose) {
            manifold.clear();
        }
        manifold.push_back(saddle.pos);
        
        if (param.verbose) {
            std::cerr << " ManBVP:  Computing phi1 (intial guess phi1=phi0+eps*v0)\n";
        }
        // move along eigenvector until p-step aligns within alpha_max
        // with the eigenvector
        lvec_type evec = (fwd ? saddle.evec[1] : saddle.evec[0]);
        if (param.verbose) std::cerr << "        Perturbation: type = " << (fwd ? "unstable" : "stable")
                                         << ", period = " << period << "\n";
        if (param.verbose) std::cerr << "                    :  evec = " << evec << ",  step = " << settings.eps
                                         << " with (" << (p_or_m ? "+":"-") << ")\n";
        lvec_type p = saddle.pos;
        bool already_aligned = false;
        int tempCount = 0;
        while (true) {
            p += settings.eps * (p_or_m ? 1 : -1) * evec;
            if (param.verbose) {
                std::cerr << "  Step " << tempCount << " (eps = " << (tempCount+1)*settings.eps << "):  p = " << p;
            }
            lvec_type map_of_p = map.map(p,period);
            if (param.verbose) {
                std::cerr << " map(p) = " << map_of_p << "\n";
            }
            lvec_type dir = _metric.displacement(p, map_of_p);
            //        std::cerr << "     at [u1-u0] = " << dir << " norm = " << nvis::norm(dir) << " ( where delta_min = " << settings.delta_min << ") \n";
            //        std::cerr << "     (test : u1-u0 = " << map_of_p-p << " and _metric.displacement(u0,u1) = " << dir << ")\n";
            dir /= nvis::norm(dir);
            lvec_type p0u0 = _metric.displacement(saddle.pos,p);
            //        std::cerr << "     (test : u0-p0 = " << p-saddle.pos << " and _metric.displacement(p0,u0) = " << p0u0 << ")\n";
            p0u0 /= nvis::norm(p0u0);
            double cosalpha = nvis::inner( p0u0 , dir );
            //        std::cerr << "     at cosAngle = " << cosalpha << " so alpha = " << acos(cosalpha) << " (alpha_max = " << settings.alpha_max << ") \n";
            if (!already_aligned && cosalpha > cos(settings.alpha_max)) {
                already_aligned = true;
                break; //break if aligned
            } else if (already_aligned) {
                break;
            }
            tempCount++;
        }
        manifold.push_back(p);
        length = _metric.distance(manifold[0], manifold[1]);
        if (length > 2.0*settings.sdelta_max) { //Max distance allowed (50x eps was in original code)
            manifold[1] = settings.sdelta_max * (p_or_m ? 1 : -1)*evec;
            progress.length = settings.sdelta_max;
            std::cerr << " ManBVP:  Initial step larger than allowed first step (2*sdelta_max).\n";
            return;
        }
        i0 = 0;
        i1 = 1;
        tau = 0;
    } else {
        length = progress.length;
        i0 = progress.segment;
        i1 = i0 + 1;
        tau = progress.tau;
    }
    
    //Turning point detector
    //TurningPoint tp(settings.liberal_alpha_max, settings.delta_tau_min, settings.update_min);
    
    //Sovle for a manifold step using BVP_Step
    int count = 0;
    while (!stop(manifold.back())) {
    
        const lvec_type& q0 = manifold[manifold.size()-2];
        const lvec_type& q1 = manifold.back();
        const lvec_type& p0 = manifold[i0];
        const lvec_type& p1 = manifold[i1];
        
        lvec_type next =
            BVP_Step<MAP,lvec_type>(map, period, _metric, p0, p1, q0, q1, settings, tau, param);
            
        //Observe output
        count++;
        double dl = _metric.distance(manifold.back(), next);
        std::cerr << " In ManBVP: Step " << count << " computed point " << next << "\n";
        std::cerr << "       i0 = " << i0 << "  i1 = " << i1 << "\n";
        std::cerr << "       phi0 = " << p0 << "  phi1 = " << p1 << "\n";
        std::cerr << "       tau = " << tau << "\n";
        std::cerr << "       u(0) = " << p0* tau+(1-tau)*p1 << " u(1) = " << next << "\n";
        std::cerr << "       current ArcLength = " << length << "\n";
        std::cerr << "       last point = " << manifold.back() << "\n";
        std::cerr << "       distance = " << dl << "\n";
        std::cerr << "       new ArcLength = " << length+dl << "\n";
        
        length += dl;
        
        //if (length > 3*_metric.diameter()/(double)fabs(period)) {
        //    std::cerr << "reached length " << length << " exceeds prescribed upper bound\n";
        //    break;
        //} //Periodic domain only -> metric.diameter() is zero for infinite domain
        manifold.push_back(next);
        
        
        /*
        // Run a check to see if we need to start a new segment (turning_point detected)
        tp.check(tau,dl);
        // Start a new segment by indicating break if startNew is attained or a stop while detected
        //  two new suitable points after singularity
        if (tp.startNew || (tp.detected && stop(manifold.back())) ) {
                //Indicate break point index in manifold set
                breakIDs.push_back((int)manifold.size()-1);
                //Move through singularity to next valid pair of points
                double tc1 = tp.tau_c;
                double tc2 = tp.tau_c;
        
                lvec_type psi0(0),psi1(0);
                bool integrationOK = false;
                int k=2, phi0ID = i0, phi1ID = i1;
                int advance = 0;
                while (!integrationOK) {
                        tc1 = tp.tau_c + (double) k*tp.delta_tau_min;
                        tc2 = tp.tau_c + (double) (k+1)*tp.delta_tau_min;
                        lvec_type a0 = tc1*manifold[phi0ID] + (1-tc1)*manifold[phi1ID];
                        lvec_type b0 = tc2*manifold[phi0ID] + (1-tc2)*manifold[phi1ID];
                        if (tc1 > 1.0) {
                          tc1 -= 1.0;
                          tc2 -= 1.0;
                          advance = 1;//both
                          phi0ID++; phi1ID++;
                          a0 = tc1*manifold[phi0ID] + (1-tc1)*manifold[phi1ID];
                          b0 = tc2*manifold[phi0ID] + (1-tc2)*manifold[phi1ID];
                        } else if (tc2 > 1.0) {
                          tc2 -= 1.0;
                          advance = 2; //Just tau_c_2
                          b0 = tc2*manifold[phi0ID+1] + (1-tc2)*manifold[phi1ID+1];
                        }
        
                        try {
                           psi0 = map.map(a0, period);
                           psi1 = map.map(b0, period);
                           integrationOK = true;
                        } catch(...) {
                           //Note: Integration usually fails in singularity
                           if (param.verbose) std::cerr << " ManBVP:  Working through Singularity\n";
                           integrationOK = false;
                        }
                        k++;
                }
                //Add two more points to manifold set to start new q0,q1
                manifold.push_back(psi0);
                manifold.push_back(psi1);
        
                //Increment working segment if necessary
                tau = tc2; //next guess
                if (advance > 0) {
                        ++i0;
                        ++i1;
                }
        
                //Reset turning point detector
                tp.reset(i0);
        }*/
        
        // End of working segment, move to next linear segment
        if (tau == 1) {
            ++i0;
            ++i1;
            tau = 0;
            //tp.reset(i0);
        }
        
        
    }
    
    std::cerr << "stop criterion tested true at " << manifold.back() << std::endl;
    
    progress.length = length;
    progress.segment = i0;
    progress.tau = tau;
}
} //end xavier


#endif
