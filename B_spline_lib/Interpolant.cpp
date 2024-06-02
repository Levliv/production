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

datatype* Interpolant::BuildInterpolant(datatype* actualGrid, int actualGridSize, int numThreads)
{
  auto* result = new datatype [actualGridSize];
  std::pair<int, int> point;
  #pragma omp parallel for private (point) num_threads(numThreads)
  for (ptrdiff_t i = 0; i < actualGridSize; ++i) {
    point = FindClosed(actualGrid[i]);
    result[i] = functionValues[point.first] +
                (functionValues[point.second] - functionValues[point.first]) *
                (interpolationFunction(actualGrid[i]) - interpolationFunction(grid[point.first])) /
                (interpolationFunction(grid[point.second]) - interpolationFunction(grid[point.first]));
  }
  return result;
}

std::pair<int, int> Interpolant::FindClosed(datatype point) {
  int leftIndex = 0;
  int rightIndex = gridSize - 1;
  int pivot;
  while (leftIndex <= rightIndex) {
    pivot = (int)((leftIndex + rightIndex) / 2);
    if (point >= grid[pivot]) {
      leftIndex = pivot + 1;
    }
    else {
      rightIndex = pivot - 1;
    }
  }
  return {leftIndex - 1, leftIndex};
}
