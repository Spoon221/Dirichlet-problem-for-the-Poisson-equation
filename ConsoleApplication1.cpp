#include <iostream>
#include <fstream>
#include <cmath>
#include <chrono>
#include <omp.h>

using namespace std;
using namespace std::chrono;

const int width = 512; // сетка
double h = 1.0 / (width - 1); // шаг сетки
double diff = 1; // разница между предыдущей и текущей итерацией
int iter = 0; // счетчик итераций
double** last_matrix; // объявление указателя на указатель предыдущей сетки

double f(double x, double y);
void JacobiMethod();
void GaussSeidelMethod();
void updateLastMatrix(double** last_matrix, double** current_matrix, int width);
void JacobiMethodParallel();
void writeResultsToFile(double** last_matrix, int width);

int main()
{
	setlocale(LC_ALL, "Russian");

	//JacobiMethod();
	//GaussSeidelMethod();
	//upperRelaxation();
	JacobiMethodParallel();
	//Dsxbcktybz2024

	return 0;
}

// вычисление правой части
double f(double x, double y)
{
	return 4 + 2 * x * x - 2 * x + 2 * y * y - 2 * y;
}

void JacobiMethod()
{
	// Цикл создания и заполнения нулями матрицы last_matrix
	last_matrix = new double* [width];
	for (int i = 0; i < width; i++)
	{
		last_matrix[i] = new double[width];
		fill(last_matrix[i], last_matrix[i] + width, 0);
	}

	// начало отсчета времени
	high_resolution_clock::time_point start_time = high_resolution_clock::now();

	// правило остановки цикла
	while (diff > 0.001)
	{
		double** current_matrix = new double* [width]; /* указатель на указатель массива на double,
														  указывающий на первый элемент массива указателей */

														  // выделение памяти под каждый элемент массива указателей 
		for (int i = 0; i < width; i++)
			current_matrix[i] = new double[width];

		double diff_max = 0; // хранилище разницы между правильным результатом и текущим значением

		// цикл по индексам оси Х в матрице
		for (int xi = 0; xi < width; xi++)
		{
			// вложенный цикл по индексам оси У в матрице
			for (int yi = 0; yi < width; yi++)
			{
				// вычисление текущих значений Х и У на основе индексов и шага сетки
				double x = xi * h;
				double y = yi * h;

				// условия проверка на граничные условия по краям сетки и установка значений в current_matrix
				if (xi == 0 || xi == width - 1)
					current_matrix[xi][yi] = y * y - y + 1;
				else if (yi == 0 || yi == width - 1)
					current_matrix[xi][yi] = x * x - x + 1;

				// в случае, когда не на краях сетки, вычисляется новое значение в точке с помощью метода итерации
				else
				{
					// вычисление по методу Якоби
					current_matrix[xi][yi] = 0.25 * (last_matrix[xi - 1][yi] + last_matrix[xi + 1][yi] + last_matrix[xi][yi - 1] + last_matrix[xi][yi + 1] - h * h * f(x, y));
					double right_result = (x * x - x + 1) * (y * y - y + 1); // вычисление правильного результата в текущей точке
					double df = abs(right_result - current_matrix[xi][yi]); /* вычисление разницы между правильным
																			   результатом и текущим значением в точке */

																			   // обновление diff_max на максимальное значение разницы в текущей итерации по всем точкам
					if (df > diff_max)
						diff_max = df;
				}
			}
		}

		cout << iter++ << " " << diff_max << "\n"; // вывод значений переменной в консоль
		diff = diff_max; // присвоение diff значение diff_max

		// копирование значений из current_matrix в last_matrix
		for (int i = 0; i < width; i++)
		{
			// перебирание всех элементов двумерных массивов и копирует значения из одного массива в другой
			for (int j = 0; j < width; j++)
				last_matrix[i][j] = current_matrix[i][j];
		}

		// освобождения памяти, выделенной под каждый вложенный массив
		for (int i = 0; i < width; i++)
			delete[] current_matrix[i];
		// удаление одномерного массива
		delete[] current_matrix;
	}

	// остановка времени и вывод результатов
	high_resolution_clock::time_point end_time = high_resolution_clock::now();
	duration<double> run_time = duration_cast<duration<double>>(end_time - start_time);
	cout << "Время выполнения метода Якоби: " << run_time.count() << " сек\n";

	// запись результатов в txt файл
	writeResultsToFile(last_matrix, width);

	// освобождение памяти под двумерный массив last_matrix
	for (int i = 0; i < width; i++)
		delete[] last_matrix[i];
	// удаление одномерного массива
	delete[] last_matrix;
}

void GaussSeidelMethod()
{
	// Цикл создания и заполнения нулями матрицы last_matrix
	last_matrix = new double* [width];
	for (int i = 0; i < width; i++)
	{
		last_matrix[i] = new double[width];
		fill(last_matrix[i], last_matrix[i] + width, 0);
	}

	// начало отсчета времени
	high_resolution_clock::time_point start_time = high_resolution_clock::now();

	// правило остановки цикла
	while (diff > 0.001)
	{
		double diff_max = 0; // хранилище разницы между правильным результатом и текущим значением

		// цикл по индексам оси Х в матрице
		for (int xi = 0; xi < width; xi++)
		{
			// вложенный цикл по индексам оси У в матрице
			for (int yi = 0; yi < width; yi++)
			{
				// вычисление текущих значений Х и У на основе индексов и шага сетки
				double x = xi * h;
				double y = yi * h;

				// условия проверка на граничные условия по краям сетки и установка значений в last_matrix
				if (xi == 0 || xi == width - 1)
					last_matrix[xi][yi] = y * y - y + 1;
				else if (yi == 0 || yi == width - 1)
					last_matrix[xi][yi] = x * x - x + 1;
				// в случае, когда не на краях сетки, вычисляется новое значение в точке с помощью метода итерации
				else
				{
					// 
					double old_value = last_matrix[xi][yi]; /* сохраняем старое значение элемента матрицы
															   last_matrix перед обновлением этого элемента */
															   // вычисление по методу Гаусса-Зейделя 
					last_matrix[xi][yi] = 0.25 * (last_matrix[xi - 1][yi] + last_matrix[xi + 1][yi] + last_matrix[xi][yi - 1] + last_matrix[xi][yi + 1] - h * h * f(x, y));
					double right_result = (x * x - x + 1) * (y * y - y + 1); // вычисление правильного результата в текущей точке
					double df = abs(right_result - last_matrix[xi][yi]); /* вычисление разницы между правильным
																			результатом и текущим значением в точке */

																			// обновление diff_max на максимальное значение разницы в текущей итерации по всем точкам
					if (df > diff_max)
						diff_max = df;

					last_matrix[xi][yi] = old_value + 1.5 * (last_matrix[xi][yi] - old_value); /* шаг линейной интерполяции
																								  для обновления значения элемента
																								  матрицы last_matrix[xi][yi] */
				}
			}
		}

		cout << iter++ << " " << diff_max << "\n"; // вывод значений переменной в консоль
		diff = diff_max; // присвоение diff значение diff_max
	}

	// остановка времени и вывод результатов
	high_resolution_clock::time_point end_time = high_resolution_clock::now();
	duration<double> run_time = duration_cast<duration<double>>(end_time - start_time);
	cout << "Время выполнения метода Гаусса-Зейделя: " << run_time.count() << " сек\n";

	// запись результатов в txt файл
	writeResultsToFile(last_matrix, width);

	// освобождение памяти под двумерный массив last_matrix
	for (int i = 0; i < width; i++)
		delete[] last_matrix[i];
	// удаление одномерного массива
	delete[] last_matrix;
}

void upperRelaxation()
{
	double omega = 1.99; // омега

	// Цикл создания и заполнения нулями матрицы last_matrix
	last_matrix = new double* [width];
	for (int i = 0; i < width; i++)
	{
		last_matrix[i] = new double[width];
		fill(last_matrix[i], last_matrix[i] + width, 0);
	}

	// начало отсчета времени
	high_resolution_clock::time_point start_time = high_resolution_clock::now();

	// правило остановки цикла
	while (diff > 0.001)
	{
		double diff_max = 0; // хранилище разницы между правильным результатом и текущим значением

		// цикл по индексам оси Х в матрице
		for (int xi = 0; xi < width; xi++)
		{
			// вложенный цикл по индексам оси У в матрице
			for (int yi = 0; yi < width; yi++)
			{
				// вычисление текущих значений Х и У на основе индексов и шага сетки
				double x = xi * h;
				double y = yi * h;

				// проверка условия на четность Х, У и номера итераций
				if ((xi + yi) % 2 == iter % 2)
				{
					// условия проверка на граничные условия по краям сетки и установка значений в last_matrix
					if (xi == 0 || xi == width - 1)
						last_matrix[xi][yi] = y * y - y + 1;
					else if (yi == 0 || yi == width - 1)
						last_matrix[xi][yi] = x * x - x + 1;
					// в случае, когда не на краях сетки, вычисляется новое значение в точке с помощью метода итерации
					else
					{
						// вычисление по методу верхней релаксации 
						double new_value = 0.25 * (last_matrix[xi - 1][yi] + last_matrix[xi + 1][yi] + last_matrix[xi][yi - 1] + last_matrix[xi][yi + 1] - h * h * f(x, y));
						double right_result = (x * x - x + 1) * (y * y - y + 1);
						double df = abs(right_result - new_value);

						// обновление diff_max на максимальное значение разницы в текущей итерации по всем точкам
						if (df > diff_max)
							diff_max = df;

						//обновление значения элемента матрицы `last_matrix[xi][yi]` с использованием метода
						last_matrix[xi][yi] = last_matrix[xi][yi] + omega * (new_value - last_matrix[xi][yi]);
					}
				}
			}
		}

		cout << iter++ << " " << diff_max << "\n"; // вывод значений переменной в консоль
		diff = diff_max; // присвоение diff значение diff_max
	}

	// остановка времени и вывод результатов
	high_resolution_clock::time_point end_time = high_resolution_clock::now();
	duration<double> run_time = duration_cast<duration<double>>(end_time - start_time);
	cout << "Время выполнения метода верхней релаксации: " << run_time.count() << " сек\n";

	// запись результатов в txt файл
	writeResultsToFile(last_matrix, width);
}

void JacobiMethodParallel()
{
	// Цикл создания и заполнения нулями матрицы last_matrix
	last_matrix = new double* [width];
	for (int i = 0; i < width; i++)
	{
		last_matrix[i] = new double[width];
		fill(last_matrix[i], last_matrix[i] + width, 0);
	}

	// начало отсчета времени
	high_resolution_clock::time_point start_time = high_resolution_clock::now();

	// правило остановки цикла
	while (diff > 0.001)
	{
		double** current_matrix = new double* [width]; /* указатель на указатель массива на double,
														  указывающий на первый элемент массива указателей */
														  // выделение памяти под каждый элемент массива указателей 
		for (int i = 0; i < width; i++)
			current_matrix[i] = new double[width];

		double diff_max = 0; // хранилище разницы между правильным результатом и текущим значением

#pragma omp parallel for collapse(2)
		// цикл по индексам оси Х в матрице
		for (int xi = 0; xi < width; xi++)
		{
			// вложенный цикл по индексам оси У в матрице
			for (int yi = 0; yi < width; yi++)
			{
				// вычисление текущих значений Х и У на основе индексов и шага сетки
				double x = xi * h;
				double y = yi * h;

				// условия проверка на граничные условия по краям сетки и установка значений в current_matrix
				if (xi == 0 || xi == width - 1)
					current_matrix[xi][yi] = y * y - y + 1;
				else if (yi == 0 || yi == width - 1)
					current_matrix[xi][yi] = x * x - x + 1;

				// в случае, когда не на краях сетки, вычисляется новое значение в точке с помощью метода итерации
				else
				{
					// вычисление по методу Якоби
					current_matrix[xi][yi] = 0.25 * (last_matrix[xi - 1][yi] + last_matrix[xi + 1][yi] + last_matrix[xi][yi - 1] + last_matrix[xi][yi + 1] - h * h * f(x, y));
					double right_result = (x * x - x + 1) * (y * y - y + 1); // вычисление правильного результата в текущей точке
					double df = abs(right_result - current_matrix[xi][yi]); /* вычисление разницы между правильным
																			   результатом и текущим значением в точке */

																			   // обновление diff_max на максимальное значение разницы в текущей итерации по всем точкам
					if (df > diff_max)
						diff_max = df;
				}
			}
		}


		if (iter % 1000 == 0)
		{
			diff = diff_max;
			cout << iter << " " << diff << endl;
		}
		iter++;
#pragma omp parallel for
		// копирование значений из current_matrix в last_matrix
		for (int i = 0; i < width; i++)
		{
			// перебирание всех элементов двумерных массивов и копирует значения из одного массива в другой
			for (int j = 0; j < width; j++)
				last_matrix[i][j] = current_matrix[i][j];
		}
#pragma omp parallel for 
		// освобождения памяти, выделенной под каждый вложенный массив
		for (int i = 0; i < width; i++)
			delete[] current_matrix[i];
		// удаление одномерного массива
		delete[] current_matrix;
	}

	// остановка времени и вывод результатов
	high_resolution_clock::time_point end_time = high_resolution_clock::now();
	duration<double> run_time = duration_cast<duration<double>>(end_time - start_time);
	cout << "Время выполнения метода Якоби: " << run_time.count() << " сек\n";

	// запись результатов в txt файл
	writeResultsToFile(last_matrix, width);

	// освобождение памяти под двумерный массив last_matrix
	for (int i = 0; i < width; i++)
		delete[] last_matrix[i];
	// удаление одномерного массива
	delete[] last_matrix;
}

void writeResultsToFile(double** last_matrix, int width)
{
	ofstream out("data.txt");
	for (int i = 0; i < width; i++)
	{
		for (int j = 0; j < width; j++)
			out << i << " " << j << " " << last_matrix[i][j] << "\n";
	}
}