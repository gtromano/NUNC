#ifndef ___NUNCSTEP_H___
#define ___NUNCSTEP_H___

#include <boost/circular_buffer.hpp>
#include <vector>
#include <tuple>
#include <list>


// Information that gets passed at each time step
typedef struct {
  long long t;
  long long changepoint;
  double Q;
  std::vector<double> zVals;
} Info;

Info nuncLongMem (const boost::circular_buffer<double>&, Info, const double&, const std::vector<double>&);
//long long nuncOnline (const boost::circular_buffer<double>&, const double&, const std::vector<double>&, const long unsigned& );
Info nuncOnline (const boost::circular_buffer<double>&,  Info, const double&, const std::vector<double>&, const std::vector<long>& );
std::vector<double> quantilesRcpp (const boost::circular_buffer<double>&, const long&, const long&, const double&);
std::vector<double> initCDF (const boost::circular_buffer<double>&, std::vector<double>&);

#endif
