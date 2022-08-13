#include <iostream>
#include <algorithm>
#include <string>
#include <numeric>      // std::iota
#include <cmath>
#include "NUNC_step.h"
#include <stdint.h>
#include <set>

// just a function to return the position of the min element in a vector
template < typename Type >
int whichMin(const std::vector<Type>& v) {
  return distance(v.begin(), std::min_element(v.begin(), v.end()));
}


// a function to return the sum of all the elements in a vector
template < typename Type >
decltype(auto) elemSum(const std::vector<Type>& v) {
  return std::accumulate(v.begin(), v.end(), 0.0);
}


// Quantiles
// as from https://stackoverflow.com/questions/11964552/finding-quartiles
// TODO: possibly need to change type of quantiles since being a vector addidion is non constant...
template<typename T>
double Lerp(T v0, T v1, T t) {
  return (1 - t)*v0 + t*v1;
}

template<typename T, typename dataT>
std::vector<T> Quantile(const dataT& inData, const std::vector<T>& probs){
  if (inData.empty()){
    return std::vector<T>();
  }
  if (1 == inData.size()){
    return std::vector<T>(1, inData[0]);
  }
  auto data = inData;
  std::sort(data.begin(), data.end());
  std::vector<T> quantiles;
  for (size_t i = 0; i < probs.size(); ++i){
    T poi = Lerp<T>(-0.5, data.size() - 0.5, probs[i]);
    size_t left = std::max(int64_t(std::floor(poi)), int64_t(0));
    size_t right = std::min(int64_t(std::ceil(poi)), int64_t(data.size() - 1));
    T datLeft = data.at(left);
    T datRight = data.at(right);
    T quantile = Lerp<T>(datLeft, datRight, poi - left);
    quantiles.push_back(quantile);
  }
  
  return quantiles;
}


template<typename T>
auto computeQuantiles(const T& inData, const int& K) {
  std::vector<double> probs(K);
  double step = 1.0 / (K + 1);
  for (size_t i = 0; i < K; ++i) {
    probs[i] = (i + 1) * step;
  }
  
  return Quantile(inData, probs);
}

// This function will give quantiles like in NP-PELT, rather than the evenly spaced quantiles
template<typename T>
auto computeQuantilesNPPELT(const T& inData, const int& K, const int& w, const double& c) {
  std::vector<double> probs(K);
  for( size_t i = 0; i < K; ++i) {
    probs[i] = 1 / (1 +  (2 * (double)w + 1) * exp( (c / K) * (2 * i - 1)) );
  }
  
  return Quantile(inData, probs);
}
// can't find transform reduce...
// double empDist (const auto& from, const auto& to, const double& q) {
//   auto len = distance(from, to);
//   //std::cout << len << std::endl;
//   if (len == 0) {len = 1;}
//   
//   double result = std::transform_reduce(from, to, from, 0.0,
//                                         std::plus<>(), [&q](auto x1, auto x2) {
//                                           return (double)(x < q) + 0.5 * (double)(fabs(x - q) < 0.0000001);
//                                         });
//   
//   return result / (double)len;
// }

template<typename T>
double empDist (const T& from, const T& to, const double& q) {
  auto len = distance(from, to);
  //std::cout << len << std::endl;
  if (len == 0) {len = 1;}
  
  double total = 0;
  std::for_each(from, to, [&q, &total](auto& x){
    total += (double)(x < q) + 0.5 * (double)(fabs(x - q) < 0.0000001);
  });
  
  return total / (double)len;
}


// computes the cdf cost (in the paper as L)
template<typename T>
double cdfCost (const double& cdfP, const T& segLen) {
  if ((cdfP == 0) | (cdfP == 1)) {
    return 0;
  }
  auto conj = 1 - cdfP;
  return segLen * (cdfP * log(cdfP) - conj * log(conj));
}

template<typename T1, typename T2>
std::vector<double> updateCDF (std::vector<double> zVals, const std::vector<double>& qVec, const T1& data, const T2& currLen) {
  
  // new point cdf
  std::vector<double> cdfPoints(zVals.size());
  transform(qVec.begin(), qVec.end(), cdfPoints.begin(), [&data] (auto& q) {
    return empDist(data.begin(), data.begin() + 1, q); // this is just the last point that it's going to leave the buffer!
  });
  
  
  std::transform(zVals.begin(), zVals.end(), cdfPoints.begin(), zVals.begin(), [&currLen](const auto& z, const auto& newP) {
    return ((currLen * z) + newP) / (currLen + 1);
  });
  
  return zVals;
}


template<typename T, typename IntT>
std::vector<double> computeSegCost (const std::vector<double>& qVec, const T& vFrom, const T& vTo, const IntT& segLen) {
  std::vector<double> segCosts(qVec.size());
  transform(qVec.begin(), qVec.end(), segCosts.begin(), [&vFrom, &vTo, &segLen] (auto& q) {
    auto cdfP =  empDist(vFrom, vTo, q);
    return cdfCost(cdfP, segLen);
  });
  
  return segCosts;
}


Info nuncLongMem (const boost::circular_buffer<double>& data, Info I, const double& threshold, const std::vector<double>& fixedQuantiles) {
  if (fixedQuantiles.size() == I.zVals.size()) {
    
    I.t++; // incrementing the t
    auto w = data.size();
    auto l = I.t - w; // last point index
    
    // for (auto i : data) {std::cout << " " << i;}
    // std::cout << std::endl;
    
    std::vector<double> TEMP(I.zVals.size()); // this is a vector needed for computing the costs
    
    // computing the current window cdf and the relative cost
    std::vector<double> currSegCDF(fixedQuantiles.size());
    transform(fixedQuantiles.begin(), fixedQuantiles.end(), currSegCDF.begin(), [&data] (auto& q) {
      return empDist(data.begin(), data.end(), q);
    });
    transform(currSegCDF.begin(), currSegCDF.end(), TEMP.begin(), [&w] (auto& cdfP) {
      return cdfCost(cdfP, w);
    });
    auto currSegCost = elemSum(TEMP);
    
    
    // computing the before the window cost (longrun cost)  
    transform(I.zVals.begin(), I.zVals.end(), TEMP.begin(), [&l] (auto& cdfP) {
      return cdfCost(cdfP, (double)l);
    });
    auto longRunCost = elemSum(TEMP);
    
    // computing the cost of the total data as a mixture of longrun CDF and window segment
    transform(I.zVals.begin(), I.zVals.end(), currSegCDF.begin(), TEMP.begin(),
              [&w, &l, &I](auto& z, auto& curr) {
                auto cdfP = (((double)l / I.t) * z) + (((double)w / I.t) * curr);
                return cdfCost(cdfP, I.t);
              });
    auto totCost = elemSum(TEMP);
    
    
    // stopping condition
    I.Q = (2 * (longRunCost + currSegCost - totCost));
    
    if (I.Q >= threshold) {
      I.changepoint = l - 1;
    }
    
    // updating the cdf
    I.zVals = updateCDF(std::move(I.zVals), fixedQuantiles, data, l - 1);
  } else {
    // std::cout << "Longrun cdf (zVals) and quantiles not of the same dimention." << std::endl;
  }
  
  return I;
}




Info nuncOnline (const boost::circular_buffer<double>& data, Info I, const double& threshold, const std::vector<double>& quant, const std::vector<long>& grid) {

  //std::cout << std::endl;
  
  I.t++; // incrementing the t

  auto w = data.size();  // window size

  auto fullWinCost = elemSum(computeSegCost(quant, data.begin(), data.end(), (double)w));

  double Q = -INFINITY;
  long cp;

  //auto maxChecks = std::min(checks, w-1);
  
  //auto low = std::max((long unsigned) 1, (w / 2) - (maxChecks / 2));
  //auto high = std::min((w / 2) + (maxChecks / 2), w-1);
  
  //for (size_t index=1; index<=maxChecks; index++) {
  // for (size_t index=low; index<=high; index++) {
  for (auto index : grid) {
    double total = 0.0;
    std::for_each(quant.begin(), quant.end(), [&data, &index, &total] (auto& q) {
      total += cdfCost(empDist(data.begin(), data.begin() + index, q), distance(data.begin(), data.begin() + index)) +
        cdfCost(empDist(data.begin() + index, data.end(), q), distance(data.begin() + index, data.end()));
    });
    auto cost = (2 * (total - fullWinCost));
    //std::cout << (cost) << ", ";
    if (cost > Q) {
      Q = cost;
      cp = index;
    }

  }

  //std::cout << std::endl;

  I.Q = Q;
  if (Q >= threshold) {
    //std::cout << "condition running "<< std::endl;
    I.changepoint = (I.t - w) + cp;
  }
  
  return I;
}



std::vector<double> quantilesRcpp (const boost::circular_buffer<double>& data, const long& K, const long& w, const double& c) {
  return(computeQuantilesNPPELT(data, K, w, c));
}

std::vector<double> initCDF (const boost::circular_buffer<double>& data, std::vector<double>& quantiles) {
  std::vector<double> CDF(quantiles.size());
  transform(quantiles.begin(), quantiles.end(), CDF.begin(), [&data] (auto& q) {
    return empDist(data.begin(), data.end(), q);
  });
  return CDF;
}
