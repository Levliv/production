#include <valarray>
#include <functional>
#include <iostream>
#include "Interpolant.h"

Interpolant::Interpolant(datatype* grid, datatype* functionValues, int gridSize, std::function<datatype (datatype)> interpolationFunction)
{
  this->grid = grid;
  this->gridSize = gridSize;
  this->functionValues = functionValues;
  this->interpolationFunction = interpolationFunction;
}

datatype* Interpolant::BuildInterpolant(datatype* actualGrid, int actualGridSize)
{
  auto* result = new datatype [actualGridSize];
  std::pair<int, int> point;
  //#pragma omp parallel for num_threads(1)
  for (int i = 0; i < actualGridSize; ++i) {
    point = FindClosed(actualGrid[i]);
    result[i] = functionValues[point.first] +
        (functionValues[point.second] - functionValues[point.first]) *
        (interpolationFunction(actualGrid[i]) - interpolationFunction(grid[point.first])) / (interpolationFunction(grid[point.second]) - interpolationFunction(grid[point.first]));
  }
  return result;
}

std::pair<int, int> Interpolant::FindClosed(datatype point) {
  for (int i = 0; i < gridSize-1; ++i) {
    if (point < grid[i]) {
      return {i-1, i};
    }
  }
  return {gridSize-2, gridSize-1};
}
