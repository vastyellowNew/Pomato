#include "standard_map.hpp"
#include <math/fixed_vector.hpp>
#include <iostream>
#include <teem/nrrd.h>
#include <sstream>

int main(int argc, char* argv[])
{
    if (argc != 5) {
        std::cerr << "USAGE: " << argv[0] << " <size> <#iter> <K> <base name>\n";
        exit(-1);
    }
    
    int size = atoi(argv[1]);
    int niter = atoi(argv[2]);
    float k = atof(argv[3]);
    nvis::vec2 step(1./(float)size, 1./(float)size);
    float* data = (float*)calloc(size*size, sizeof(float));
    float* sf = (float*)calloc(size*size, sizeof(float));
    
    srand48(time(0));
#pragma openmp parallel for
    for (int i=0 ; i<size*size ; ++i) {
        nvis::vec2 c(i%size, i/size);
        nvis::vec2 x = c*step;
        xavier::standard_map stdmap(k);
        std::vector<nvis::vec2> orbit;
        stdmap.map(x, orbit, niter);
        
        sf[i] = (orbit.back()[0]-x[0])/(float)niter;
        data[i] = stdmap.winding_number(x, orbit);
    }
    
    std::ostringstream os;
    os << argv[4] << "-avg_wn.nrrd";
    
    size_t sz[] = {size, size};
    double spc[] = {step[0], step[1]};
    Nrrd* nout = nrrdNew();
    nrrdWrap_nva(nout, data, nrrdTypeFloat, 2, sz);
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, spc);
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoMin, 0);
    nrrdSave(os.str().c_str(), nout, NULL);
    nrrdNuke(nout);
    
    os.clear();
    os.str("");
    os << argv[4] << "-direct_wn.nrrd";
    nout = nrrdNew();
    nrrdWrap_nva(nout, sf, nrrdTypeFloat, 2, sz);
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoSpacing, spc);
    nrrdAxisInfoSet_nva(nout, nrrdAxisInfoMin, 0);
    nrrdSave(os.str().c_str(), nout, NULL);
    nrrdNuke(nout);
    
    return 0;
}