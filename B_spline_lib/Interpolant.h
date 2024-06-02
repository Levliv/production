#ifndef INTERPOLATION_PROJECT_INTERPOLANT_H
#define INTERPOLATION_PROJECT_INTERPOLANT_H

#include <vector>
#include <functional>
typedef long double datatype;

class Interpolant
{
public:
  explicit Interpolant(datatype*, datatype*, int, std::function<datatype (datatype)>);
  datatype* BuildInterpolant(datatype*, int, int);

private:
  datatype* grid;
  int gridSize;
  datatype* functionValues;
  std::function<datatype (datatype)> interpolationFunction;
  std::pair<int, int> FindClosed(datatype);
};
#endif //INTERPOLATION_PROJECT_INTERPOLANT_H
