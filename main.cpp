#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>

typedef long double datatype;

using namespace std;

long double func(datatype point){
    return point*sin(10*point);
}

long double spline_func(datatype point){
    return pow(exp(point), 3);
}

long double B_spline_func(datatype point){
    return pow(point, 3);
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

void spline_create(){
    int N = 20; // количество промежутков апроксимации
    int m = 1; // степень сплайна (порядок - 1)
    // расширенное множество узлов N+1+2*m
    datatype a = 0.1; // левый конец промежутка
    datatype b = 3; // правый конец промежутка
    vector<datatype> grid(N + 1 + 2 * m);
    vector<datatype> func_values(N + 1 + m);
    int grid_multiplier = 5;
    vector<datatype> small_grid((N + 1 + m) * grid_multiplier);
    vector<datatype> spline_values((N + 1 + m) * grid_multiplier);
    vector<datatype> b_spline_values((N + 1 + m) * grid_multiplier);
    datatype max_error = 0.;
    datatype current_error = 0.;
    datatype b_spline_max_error = 0.;
    datatype b_spline_current_error = 0.;

    // Generating a non-uniform grid with multiple knots
    grid = Grid_generator(a, b, N, m);

    // Counting func values on grid knots
    for (int i = m; i < N + 1 + m; ++i)
        func_values[i] = func(grid[i]);

    std::ofstream out;
    out.open("data.txt");
    if (!out.is_open())
    {
        cout << "Error! File is not opened";
    }

    //Counting Spline values for shallower grid
    for (int i = m; i < N + m; ++i){
        for (int j = 0; j <= grid_multiplier; ++j){
            int idx = i * grid_multiplier + j;
            datatype step = (grid[i+1] - grid[i]) / grid_multiplier;
            small_grid[idx] = grid[i]  + step * j;
            spline_values[idx] = func_values[i] + (func_values[i+1] - func_values[i]) * (spline_func(small_grid[idx]) -
                                                                                         spline_func(grid[i])) / (spline_func(grid[i+1]) - spline_func(grid[i]));
            b_spline_values[idx] = func_values[i] + (func_values[i+1] - func_values[i]) * (B_spline_func(small_grid[idx]) -
                    B_spline_func(grid[i])) / (B_spline_func(grid[i+1]) - B_spline_func(grid[i]));
        }
    }

    // Counting an error for spline func on this shallower grid
    for (int i = m*grid_multiplier; i <= (N + m)*grid_multiplier; ++i){
        current_error = abs(func(small_grid[i]) - spline_values[i]);
        b_spline_current_error = abs(func(small_grid[i]) - b_spline_values[i]);
        out << small_grid[i] << " " << func(small_grid[i]) << " " << spline_values[i] << " " << b_spline_values[i] << std::endl;
        if (max_error < current_error)
            max_error = current_error;
        if (b_spline_max_error < b_spline_current_error)
            b_spline_max_error = b_spline_current_error;
    }
    out.close();
    cout << "Bphi spline max error on shallow grid: " << max_error << endl;
    cout << "B spline max error on shallow grid: " << b_spline_max_error << endl;
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
            cout << shallowGrid[idx] << endl;
            out << shallowGrid[idx] << " " << (spline_func(shallowGrid[idx]) - spline_func(x[i])) / (spline_func(x[i+1]) - spline_func(x[i])) << endl;
        }
    }
    //vector<datatype> weight();
    out.close();
}

int main() {
    spline_create();
    //B_spline_creator();
    return 0;
}
