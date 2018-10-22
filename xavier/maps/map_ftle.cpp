#include <iostream>
#include <iomanip>
#include <vector>
#include <complex>
#include <sstream>
#include <math.h>
#include <iostream>
#include <list>

#include <boost/format.hpp>
#include <boost/limits.hpp>

// chris
#include <math/fixed_vector.hpp>
#include <math/fixed_matrix.hpp>
#include <util/wall_timer.hpp>

// xavier
#include <data/grid.hpp>
#include <data/raster_data.hpp>
#include <util/wall_timer.hpp>
#include <image/nrrd_wrapper.hpp>
#include "standard_map.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace nvis;

double eps;
int resx, resy, maxp, maxit, it_step;
char* outs, *file, *ts;
double minx, maxx, miny, maxy, _K;
double _minx, _maxx, _miny, _maxy;
double hx, hy;

void initialize(int argc, char* argv[])
{
    hestOpt* hopt = NULL;
    hestParm* hparm;
    airArray* mop;
    char* me;
    
    mop = airMopNew();
    me = argv[0];
    hparm = hestParmNew();
    airMopAdd(mop, hparm, AIR_CAST(airMopper, hestParmFree), airMopAlways);
    hparm->elideSingleOtherType = AIR_TRUE;
    hestOptAdd(&hopt, "k",      "K",                        airTypeDouble,  0, 1, &_K,          "0.5",      "K parameter");
    hestOptAdd(&hopt, "o",      "output",               airTypeString,  1, 1, &outs,        NULL,       "output name");
    hestOptAdd(&hopt, "rx",     "resolution x",     airTypeInt,         0, 1, &resx,        "1024",     "sampling resolution in X");
    hestOptAdd(&hopt, "ry",     "resolution y",     airTypeInt,         0, 1, &resy,        "1024",     "sampling resolution in Y");
    hestOptAdd(&hopt, "maxi",       "max iterations",       airTypeInt,         0, 1, &maxit,       "10",           "max number of map iterations");
    hestOptAdd(&hopt, "s",          "iterations step",  airTypeInt,         0, 1, &it_step, "0",            "iteration step size");
    hestOptAdd(&hopt, "minx",       "min x coord",          airTypeDouble,  0, 1, &minx,        "-10000",   "min x in bounding box");
    hestOptAdd(&hopt, "maxx",       "max x coord",          airTypeDouble,  0, 1, &maxx,        "10000",        "max x in bounding box");
    hestOptAdd(&hopt, "miny",       "min y coord",          airTypeDouble,  0, 1, &miny,        "-10000",   "min y in bounding box");
    hestOptAdd(&hopt, "maxy",       "max y coord",          airTypeDouble,  0, 1, &maxy,        "10000",        "max y in bounding box");
    
    hestParseOrDie(hopt, argc - 1, (const char**)argv + 1, hparm,
                   (const char*)me, "Compute FTLE value of standard map after a given number of iterations. Intermediate steps are saved to disk if requested.",
                   AIR_TRUE, AIR_TRUE, AIR_TRUE);
}

inline double lmax(int n, const std::vector< nvis::vec2 >& pos, const xavier::metric<double, 2>& metric)
{
    static double hx = (maxx - minx) / (double)resx;
    static double hy = (maxy - miny) / (double)resy;
    
    unsigned int i = n % resx;
    unsigned int j = n / resx;
    if (i == 0 || i == resx - 1 || j == 0 || j == resy - 1) {
        return -1.0;
    }
    
    // look for valid neighboring values to compute derivatives
    nvis::vec2 Jx = 1. / (2.*hx) * metric.displacement(pos[n-1], pos[n+1]);
    nvis::vec2 Jy = 1. / (2.*hy) * metric.displacement(pos[n-resx], pos[n+resx]);
    double a = nvis::inner(Jx, Jx);
    double b = nvis::inner(Jx, Jy);
    double c = nvis::inner(Jy, Jy);
    
    double lmaj = 0.5 * (a * a + 2 * b * b + c * c + (a + c) * sqrt((a - c) * (a - c) + 4 * b * b));
    
    if (lmaj < 1 && !(n%1000)) {
        std::ostringstream os;
        os << "Jx = " << Jx << ", Jy = " << Jy << std::endl;
        os << "pos[" << n+1 << "] = " << pos[n+1] << ", pos[" << n-1 << "] = " << pos[n-1] << std::endl;
        std::cerr << os.str();
    }
    
    return lmaj;
}

using namespace xavier;

double __mod(double a, double b)
{
    return a >= 0 ? fmod(a, b) : b + fmod(a, b);
}

typedef grid<double, 3>                                 grid_type;
typedef raster_data<nvis::vec3, double, 3>      field_type;
typedef nvis::ivec3                                     ivec_type;
typedef xavier::standard_map                            map_type;

int main(int argc, char* argv[])
{
    initialize(argc, argv);
    
#if _OPENMP
    std::cout << omp_get_max_threads() << " threads available\n";
#endif
    
    std::cout << "parameters: max period = " << maxp
              << ", max iterations = " << maxit << ", resolution = " << resx << " x " << resy
              << ", eps = " << eps << '\n';
              
    map_type map(_K);
    
    nvis::bbox2 _bounds(nvis::vec2(0,0), nvis::vec2(1,1));
    
    std::cerr << "bounding box = " << _bounds << '\n';
    
    bool per[2] = {true, true};
    xavier::metric<double, 2> metric2d(_bounds, per);
    
    _minx = _bounds.min()[0];
    _miny = _bounds.min()[1];
    _maxx = _bounds.max()[0];
    _maxy = _bounds.max()[1];
    
    minx = std::max(_minx, minx);
    maxx = std::min(_maxx, maxx);
    miny = std::max(_miny, miny);
    maxy = std::min(_maxy, maxy);
    
    assert(minx < maxx && miny < maxy);
    
    hx = (maxx - minx) / (double)resx;
    hy = (maxy - miny) / (double)resy;
    
    unsigned int blocked = 0;
    float last_pct = 0.;
    
    std::cout << "minx = " << minx << ", miny = " << miny << ", maxx = " << maxx << ", maxy = " << maxy
              << ", hx = " << hx << ", hy = " << hy << std::endl;
    std::cout << "it_step = " << it_step << std::endl;
    
    nvis::timer timer;
    
    std::cout << "initializing coordinates\n";
    std::vector<nvis::vec2> pos_fwd(resx*resy), pos_bwd(resx*resy);
    for (int n = 0 ; n < resx*resy ; ++n) {
        unsigned int i = n % resx;
        unsigned int j = n / resx;
        nvis::vec2 x(minx + hx*(double)i, miny + hy*(double)j);
        pos_fwd[n] = pos_bwd[n] = x;
    }
    
    if (it_step == 0) {
        it_step = maxit;
    }
    
    for (int iteration = it_step ; iteration <= maxit ; iteration += it_step) {
        std::cout << "Computing map..." << std::endl;
        std::cout << "\niteration " << iteration << " from " << maxit << '\n';
        
        // do that in parallel
        #pragma omp parallel
        {
            int thread_id = 0;
            #pragma omp for schedule(dynamic,1)
            for (int n = 0 ; n < resx*resy ; ++n) {
#if _OPENMP
                thread_id = omp_get_thread_num();
#endif
                if (thread_id == 0) {
                    float pct = 100 * n / (resx * resy);
                    std::ostringstream os;
                    os << (int)pct << "% forward integration completed      \r" << std::flush;
                    std::cout << os.str();
                }
                
                const map_type* clone = map.clone();
                
                nvis::vec2 prev = pos_fwd[n];
                
                try {
                    pos_fwd[n] = clone->map(pos_fwd[n], it_step);
                } catch (...) {}
            }
        }
        std::cerr << '\n';
        
        #pragma omp parallel
        {
            int thread_id = 0;
            #pragma omp for schedule(dynamic,1)
            for (int n = 0 ; n < resx*resy ; ++n) {
#if _OPENMP
                thread_id = omp_get_thread_num();
#endif
                if (thread_id == 0) {
                    float pct = 100 * n / (resx * resy);
                    std::cout << (int)pct << "% backward integration completed      \r";
                }
                
                const map_type* clone = map.clone();
                
                unsigned int i = n % resx;
                unsigned int j = n / resx;
                nvis::vec2 x(minx + hx*(double)i, miny + hy*(double)j);
                
                try {
                    pos_bwd[n] = clone->map(pos_bwd[n], -it_step);
                } catch (...) {}
            }
        }
        std::cerr << '\n';
        
        std::cout << "computing FTLE...\n";
        float* _ftle = (float*)calloc(3 * resx * resy, sizeof(float));
        float* _lmax = (float*)calloc(3 * resx * resy, sizeof(float));
        #pragma omp parallel
        {
            int thread_id = 0;
            #pragma omp for schedule(dynamic,1)
            for (int n = 0 ; n < resx*resy ; ++n) {
#if _OPENMP
                thread_id = omp_get_thread_num();
#endif
                double lmax_f = 0., lmax_b = 0.;
                double ftle_f = 0., ftle_b = 0.;
                double lm = lmax(n, pos_fwd, metric2d);
                
                lmax_f = log(std::max(1., lm));
                ftle_f = lmax_f;
                
                lm = lmax(n, pos_bwd, metric2d);
                lmax_b = log(std::max(1., lm));
                ftle_b = lmax_b;
                
                if (!std::isinf(ftle_f) && !std::isnan(ftle_f)) {
                    _ftle[3*n] = ftle_f;
                }
                if (!std::isinf(ftle_b) && !std::isnan(ftle_b)) {
                    // std::cerr << ftle_b << " in backward\n";
                    _ftle[3*n+2] = ftle_b;
                }
            }
        }
        
        std::cerr << timer.elapsed() << " seconds needed for the computation\n";
        
        std::ostringstream os;
        os << outs << "-ftle-p=" << iteration << "_from_" << maxit
           << "-res=" << resx << "x" << resy
           << ".nrrd";
        Nrrd* nout = nrrdNew();
        size_t  size[]  = {3, resx, resy};
        double  spc[]   = {airNaN(), hx, hy};
        int     kind[]  = {nrrdKindUnknown, nrrdKindSpace, nrrdKindSpace};
        int     center[] = {nrrdCenterUnknown, nrrdCenterCell, nrrdCenterCell};
        if (nrrdWrap_nva(nout, _ftle, nrrdTypeFloat, 3, size)) {
            std::cout << "ERROR while wrapping data: " << biffGetDone(NRRD)
                      << std::endl;
            if (_ftle) {
                delete[] _ftle;
            } else {
                nrrdNuke(nout);
            }
            exit(-1);
        }
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoKind, kind);
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoCenter, center);
        nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, spc);
        if (nrrdSave(os.str().c_str(), nout, NULL)) {
            std::cout << "ERROR while exporting file: " << biffGetDone(NRRD)
                      << std::endl;
            exit(-1);
        }
        nrrdNuke(nout);
        std::cout << "exported " << os.str() << std::endl;
    }
    
    return 0;
}






























