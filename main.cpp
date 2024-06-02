#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include "Interpolant.h"

typedef long double datatype;

using namespace std;

datatype func(datatype point){
    return atan(point);
}

datatype quasiLinerSplineFunc(datatype point){
    return tanh(point);
}

datatype linearSplineFunc(datatype point){
    return point;
}

datatype getFuncDerivative(datatype point)
{
  datatype offset = 0.000001;
  return (func(point + offset) - func(point - offset)) / (2 * offset);
}

vector<datatype> Grid_generator(datatype left_border, // левый конец промежутка
                                datatype right_border, // правый конец промежутка
                                int number_of_approx_segments, // Количество промежутков аппроксимации
                                int spline_order // степень сплайна (порядок - 1)
                                ){
    datatype offset = 0.0001;
    vector<datatype> grid(number_of_approx_segments + 1 + 2 * spline_order);
    // Основная часть отрезка
    for (int i = spline_order; i < number_of_approx_segments + 1 + spline_order; ++i) {
        grid[i] = left_border + (i - spline_order) * (right_border - left_border) / number_of_approx_segments;
    }
    // Левая граница отрезка
    for (int i = 0; i < spline_order; ++i) {
        grid[i] = grid[spline_order] - (spline_order - i) * offset;
    }
    // Правая граница отрезка
    for (int i = number_of_approx_segments + 1 + spline_order; i < number_of_approx_segments + 1 + 2 * spline_order; ++i) {
        grid[i] = grid[number_of_approx_segments + spline_order] + (i - number_of_approx_segments - spline_order) * offset;
    }
    return grid;
}

datatype* GridGenerator(datatype left_border, // левый конец промежутка
                                datatype right_border, // правый конец промежутка
                                int number_of_approx_segments, // Количество промежутков аппроксимации
                                int spline_order // степень сплайна (порядок - 1)
                                ){
  datatype offset = 0.0001;
  auto* grid = new datatype[number_of_approx_segments + 1 + 2 * spline_order];
  // Основная часть отрезка
  for (int i = spline_order; i < number_of_approx_segments + 1 + spline_order; ++i) {
    grid[i] = left_border + (i - spline_order) * (right_border - left_border) / number_of_approx_segments;
  }
  // Левая граница отрезка
  for (int i = 0; i < spline_order; ++i) {
    grid[i] = grid[spline_order] - (spline_order - i) * offset;
  }
  // Правая граница отрезка
  for (int i = number_of_approx_segments + 1 + spline_order; i < number_of_approx_segments + 1 + 2 * spline_order; ++i) {
    grid[i] = grid[number_of_approx_segments + spline_order] + (i - number_of_approx_segments - spline_order) * offset;
  }
  return grid;
}

void testTimePerThread(datatype leftBorder, datatype rightBorder, int spDegree, int multiplier, int numberRuns, int numberOfThreadsUsed){
  int minNumberOfSegments = 100000;
  int maxNumberOfSegments = 300000;
  int pointsSegmentStep = (maxNumberOfSegments - minNumberOfSegments) / 20;
  std::ofstream out;
  string fileName = "time_test_threads_" + to_string(numberOfThreadsUsed) + ".txt";
  out.open(fileName);
  if (!out.is_open()) {
    cerr << "File was not opened!" << std::endl;
  }
  for (int numberOfSegments = minNumberOfSegments; numberOfSegments <= maxNumberOfSegments; numberOfSegments += pointsSegmentStep) {
    auto *grid = new datatype[numberOfSegments + 2 * spDegree + 1];
    auto *functionValues = new datatype[numberOfSegments + 2 * spDegree + 1];

    datatype offset = 0.0001;
    // Основная часть отрезка
    for (int i = spDegree; i < numberOfSegments + 1 + spDegree; ++i) {
      grid[i] = leftBorder + (i - spDegree) * (rightBorder - leftBorder) / numberOfSegments;
    }
    // Левая граница отрезка
    for (int i = 0; i < spDegree; ++i) {
      grid[i] = grid[spDegree] - (spDegree - i) * offset;
    }
    // Правая граница отрезка
    for (int i = numberOfSegments + 1 + spDegree; i < numberOfSegments + 1 + 2 * spDegree; ++i) {
      grid[i] = grid[numberOfSegments + spDegree] + (i - numberOfSegments - spDegree) * offset;
    }

    for (int i = 0; i < numberOfSegments + 2 * spDegree + 1; ++i) {
      functionValues[i] = func(grid[i]);
    }

    auto *interpolant = new Interpolant(grid, functionValues, numberOfSegments + 2 * spDegree + 1, quasiLinerSplineFunc);

    datatype step;
    auto *shallowGrid = new datatype[numberOfSegments * multiplier + 1];
    for (int i = 0; i < numberOfSegments * multiplier + 1; ++i) {
      step = (grid[i / multiplier + spDegree + 1] - grid[i / multiplier + spDegree]) / multiplier;
      shallowGrid[i] = grid[i / multiplier + spDegree] + i % multiplier * step;
    }
    auto *ags = new datatype[numberOfSegments * multiplier + 1];
    auto *timeData = new double[numberRuns];
    double sum = 0;
    double sumDistSquare = 0;
    for (int i = 0; i < numberRuns; ++i) {
      clock_t start = clock();
      ags = interpolant->BuildInterpolant(shallowGrid, numberOfSegments * multiplier + 1, numberOfThreadsUsed);
      clock_t end = clock();
      timeData[i] = (double) (end - start) / CLOCKS_PER_SEC;
    }
    for (int i = 0; i < numberRuns; ++i) {
      sum += timeData[i];
    }

    double average = sum / numberRuns;

    for (int i = 0; i < numberRuns; ++i) {
      sumDistSquare += (timeData[i] - average) * (timeData[i] - average);
    }
    double mse = sqrt(sumDistSquare / numberRuns);

    out << numberOfSegments << " " << average << " " << mse << std::endl;
  }
  out.close();
  std::cout << "File with time test results for " << numberOfThreadsUsed << " threads successfully saved in build directory!" << std::endl;
}

void testTime(datatype leftBorder, datatype rightBorder, int spDegree, int multiplier, int numberRuns){
  testTimePerThread(leftBorder, rightBorder, spDegree, multiplier, numberRuns, 1);
  testTimePerThread(leftBorder, rightBorder, spDegree, multiplier, numberRuns, 2);
  testTimePerThread(leftBorder, rightBorder, spDegree, multiplier, numberRuns, 4);
}

void test_error(datatype leftBorder, datatype rightBorder, int numberOfSegments, int spDegree, int multiplier)
{
  auto* grid = new datatype[numberOfSegments + 2 * spDegree + 1];
  auto* functionValues = new datatype[numberOfSegments + 2 * spDegree + 1];

  datatype offset = 0.0001;
  // Основная часть отрезка
  for (int i = spDegree; i < numberOfSegments + 1 + spDegree; ++i) {
    grid[i] = leftBorder + (i - spDegree) * (rightBorder - leftBorder) / numberOfSegments;
  }
  // Левая граница отрезка
  for (int i = 0; i < spDegree; ++i) {
    grid[i] = grid[spDegree] - (spDegree - i) * offset;
  }
  // Правая граница отрезка
  for (int i = numberOfSegments + 1 + spDegree; i < numberOfSegments + 1 + 2 * spDegree; ++i) {
    grid[i] = grid[numberOfSegments + spDegree] + (i - numberOfSegments - spDegree) * offset;
  }

  for (int i = 0; i < numberOfSegments + 2 * spDegree + 1; ++i) {
    functionValues[i] = func(grid[i]);
  }

  auto* interpolant = new Interpolant(grid, functionValues, numberOfSegments + 2 * spDegree + 1, quasiLinerSplineFunc);
  auto* polynomial_interpolant = new Interpolant(grid, functionValues, numberOfSegments + 2 * spDegree + 1,
                                                 linearSplineFunc);

  datatype step;
  auto* shallowGrid = new datatype [numberOfSegments * multiplier + 1];
  for (int i = 0; i < numberOfSegments * multiplier + 1; ++i) {
    step = (grid[i/multiplier + spDegree + 1] - grid[i/multiplier + spDegree]) / multiplier;
    shallowGrid[i] = grid[i / multiplier + spDegree] + i % multiplier * step;
  }
  auto* ags = new datatype [numberOfSegments * multiplier + 1];
  auto* polynomial_ags = new datatype [numberOfSegments * multiplier + 1];
  clock_t start = clock();
  ags = interpolant->BuildInterpolant(shallowGrid, numberOfSegments * multiplier + 1, 1);
  clock_t end = clock();
  double seconds = (double)(end - start) / CLOCKS_PER_SEC;
  std::cout << "Consumed time: " << seconds << endl;

  start = clock();
  polynomial_ags = polynomial_interpolant->BuildInterpolant(shallowGrid, numberOfSegments * multiplier + 1, 1);
  end = clock();
  seconds = (double)(end - start) / CLOCKS_PER_SEC;
  std::cout << "Consumed time: " << seconds << endl;

  datatype currentError = abs(func(shallowGrid[0]) - ags[0]);
  datatype maxError = currentError;

  for (int i = 0; i < numberOfSegments * multiplier + 1; ++i)
  {
    currentError = std::abs(ags[i] - func(shallowGrid[i]));
    if (currentError > maxError)
    {
      maxError = currentError;
    }
  }
  std::cout << "Quasi linear spline max error: " << maxError << std::endl;

  currentError = abs(func(shallowGrid[0]) - polynomial_ags[0]);
  maxError = currentError;

  for (int i = 0; i < numberOfSegments * multiplier + 1; ++i)
  {
    currentError = std::abs(polynomial_ags[i] - func(shallowGrid[i]));
    if (currentError > maxError)
    {
      maxError = currentError;
    }
  }
  std::cout << "Linear spline max error: " << maxError << std::endl;
}

void adaptiveGrid(datatype left_border, // левый конец промежутка
                  datatype right_border, // правый конец промежутка
                  int number_of_approx_segments, // Количество промежутков аппроксимации
                  int spline_order // степень сплайна (порядок - 1))
){
  auto* derivativeGrid = new datatype[number_of_approx_segments + 2 * spline_order + 1];
  derivativeGrid = GridGenerator(left_border, right_border, number_of_approx_segments, spline_order);
  auto* derivativeValues = new datatype[number_of_approx_segments + 2 * spline_order + 1];
  for (int i = 0; i < number_of_approx_segments + 2 * spline_order + 1; ++i) {
    derivativeValues[i] = getFuncDerivative(derivativeGrid[i]);
  }
  auto* interpolant = new Interpolant(derivativeGrid, derivativeValues,
                                      number_of_approx_segments + 2 * spline_order + 1, linearSplineFunc);
  
}

int main() {
  testTime(0.1, 6000, 1, 10, 1);

  return 0;
}







void B_spline_creator(){
  int N = 5; // количество промежутков апроксимации
  int m = 1; // степень сплайна (порядок - 1)
  datatype a = 0.1; // левый конец промежутка
  datatype b = 0.6; // правый конец промежутка
  vector<datatype> x(N + 1 + 2 * m);
  vector<datatype> y(N + 1 + 2 * m);
  x = Grid_generator(a, b, N, m);
  int multiplier = 20;
  vector<datatype> shallowGrid((N + 1 + m) * multiplier);
  for (int i = m; i < N + m; ++i){
    for (int j = 0; j <= multiplier; ++j){
      int idx = i * multiplier + j;
      datatype step = (x[i+1] - x[i]) / multiplier;
      shallowGrid[idx] = x[i]  + step * j;
    }
  }
  std::ofstream out;
  out.open("B_spline_data.txt");
  if (!out.is_open())
  {
    cout << "Error! File is not opened";
  }
  // Нарисуем первый Bphi сплайн
  for (int i = m; i < m + 1; ++i){
    for (int j = 0; j <= multiplier; ++j){
      int idx = i * multiplier + j;
      //cout << shallowGrid[idx] << endl;
      out << shallowGrid[idx] << " " << (quasiLinerSplineFunc(shallowGrid[idx]) - quasiLinerSplineFunc(x[i])) / (quasiLinerSplineFunc(
          x[i + 1]) -
                                                                                                                 quasiLinerSplineFunc(x[i])) << endl;
    }
  }
  //vector<datatype> weight();
  out.close();
}
