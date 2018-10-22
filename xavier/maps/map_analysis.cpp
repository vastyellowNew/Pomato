/* map_analysis functions
   Author:  Wayne Schlei & Xavier Tricoche
   Date:  7/1/2013
   Implementation source for map_analysis.hpp (functions without templates)
*/
#include <algorithm>
#include <boost/rational.hpp>
//API
#include <math/rational.hpp>
#include <maps/map_analysis.hpp>
//Avizo
#ifdef QT_CLEAN_NAMESPACE
#include <QString>
#endif

namespace xavier {

double average_distance(const std::vector<nvis::vec2>& steps, int p,
                        const metric_type& _metric)
{
    if ((int) steps.size() < p) {
        return std::numeric_limits<double>::max();
    }
    double d = 0;
    int count = 0;
    for (int i=0 ; i<(int)steps.size() ; ++i) {
        if (i+p<(int)steps.size()) {
            d += _metric.distance(steps[i], steps[i+p]);
            ++count;
        }
        if (i>=p) {
            d += _metric.distance(steps[i-p], steps[i]);
            ++count;
        }
    }
    return d/(double)count;
}

//Best periods based on the average distance -> Works best for a periodic domain
int best_period(const std::vector<vec_type>& steps, int maxp, const metric_type& m)
{
    std::map<double, int> p2d;
    for (int p=1 ; p<=maxp ; ++p) {
        p2d[average_distance(steps, p, m)] = p;
    }
    return p2d.begin()->second;
}

//Computing best approximate toroidal period based on average distance of all points
// -> Works best for a periodic domain
void best_periods(std::vector<int>& periods, const std::vector<vec_type>& steps,
                  int maxp, int nper, const metric_type& m)
{
    typedef std::map<double, int>  map_type;
    map_type p2d;
    for (int p=1 ; p<=maxp ; ++p) {
        double d = average_distance(steps, p, m);
        p2d[d] = p;
        std::cerr << "period " << p << " -> average distance = " << d << '\n';
    }
    size_t counter=0;
    map_type::const_iterator it;
    periods.clear();
    std::ostringstream os;
    os << "best_periods: ";
    for (it=p2d.begin() ; it!=p2d.end() && (int)counter<nper; ++it, ++counter) {
        os << it->second << ", ";
        periods.push_back(it->second);
    }
    os << '\n';
    std::cerr << os.str();
}




} //end xavier
