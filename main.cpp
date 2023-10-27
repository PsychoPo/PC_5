#include <iostream>
#include <omp.h>
#include <vector>
#include <fstream>
#include <cmath>
#include <tbb/parallel_for.h>
#include <tbb/info.h>
#include <tbb/tbb.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_reduce.h>
#include <iomanip>
#include "Header.h"
#include "Header1.h"

using namespace std;
using namespace tbb;

//______________________НАЧАЛО СУММЫ______________________//
inline void sum_posled(double** arr1, double** arr2, int size1, int size2) {

	double** arr_res = new double* [size1];
	for (int i = 0; i < size1; i++)
		arr_res[i] = new double[size2];

	for (size_t i = 0; i < size1; i++) {
		for (size_t j = 0; j < size2; j++) {
			arr_res[i][j] = arr1[i][j] + arr2[i][j];
		}
	}

	for (int i = 0; i < size1; i++)
		delete[] arr_res[i];
	delete[] arr_res;
}

//Параллельное сложение openmp parallel for
inline void sum_par_omp(double** arr1, double** arr2, int size1, int size2) {

	double** arr_res = new double* [size1];
	for (int i = 0; i < size1; i++)
		arr_res[i] = new double[size2];

#pragma omp parallel for
	for (long i = 0; i < size1; i++) {
		for (long j = 0; j < size2; j++) {
			arr_res[i][j] = arr1[i][j] + arr2[i][j];
		}
	}

	for (int i = 0; i < size1; i++)
		delete[] arr_res[i];
	delete[] arr_res;
}

//Параллельное сложение openmp sections
inline void sum_par_omp_sections(double** arr1, double** arr2, int size1, int size2) {

	int num_thread = omp_get_max_threads();

	double** arr_res = new double* [size1];

	for (int i = 0; i < size1; i++)
		arr_res[i] = new double[size2];

	int p0 = 0,
		p1 = size1 / num_thread,
		p2 = 2 * size1 / num_thread,
		p3 = 3 * size1 / num_thread,
		p4 = size1;

#pragma omp parallel sections
	{
#pragma omp section
		{
			for (int i = p0; i < p1; i++) {
				for (int j = 0; j < size1; j++) {
					arr_res[i][j] = arr1[i][j] + arr2[i][j];
				}
			}
		}
#pragma omp section
		{
			for (int i = p1; i < p2; i++) {
				for (int j = 0; j < size1; j++) {
					arr_res[i][j] = arr1[i][j] + arr2[i][j];
				}
			}
		}
#pragma omp section
		{
			if (num_thread > 2)
				for (int i = p2; i < p3; i++)
					for (int j = 0; j < size1; j++) {
						arr_res[i][j] = arr1[i][j] + arr2[i][j];
					}
		}
#pragma omp section
		{
			if (num_thread > 3)
				for (int i = p3; i < p4; i++)
					for (int j = 0; j < size1; j++) {
						arr_res[i][j] = arr1[i][j] + arr2[i][j];
					}
		}
	}

	for (int i = 0; i < size1; i++)
		delete[] arr_res[i];
	delete[] arr_res;

}

//Параллельное сложение oneTBB task_group
inline void sum_par_tbb_task_group(double** arr1, double** arr2, int size1, int size2) {

	int num_thread = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);

	double** arr_res = new double* [size1];
	for (int i = 0; i < size1; i++)
		arr_res[i] = new double[size2];

	int p0 = 0,
		p1 = size1 / num_thread,
		p2 = 2 * size1 / num_thread,
		p3 = 3 * size1 / num_thread,
		p4 = size1;

tbb:task_group g;

	g.run([&] {	for (long i = p0; i < p1; i++) {
		for (int j = 0; j < size1; j++) {
			arr_res[i][j] = arr1[i][j] + arr2[i][j];
		}
	}; });

	g.run([&] {	for (long i = p1; i < p2; i++) {
		for (int j = 0; j < size1; j++) {
			arr_res[i][j] = arr1[i][j] + arr2[i][j];
		}
	}; });
	if (num_thread > 2)
		g.run([&] {	for (long i = p2; i < p3; i++) {
		for (int j = 0; j < size1; j++) {
			arr_res[i][j] = arr1[i][j] + arr2[i][j];
		}
	}; });
	if (num_thread > 3)
		g.run([&] {	for (long i = p3; i < p4; i++) {
		for (int j = 0; j < size1; j++) {
			arr_res[i][j] = arr1[i][j] + arr2[i][j];
		}
	}; });


	g.wait();

	for (int i = 0; i < size1; i++)
		delete[] arr_res[i];
	delete[] arr_res;
}

//Параллельное сложение oneTBB parallel_for лямбда выражение
inline void sum_par_tbb_par_for_lambda(double** arr1, double** arr2, int size1, int size2) {

	double** arr_res = new double* [size1];
	for (int i = 0; i < size1; i++)
		arr_res[i] = new double[size2];

	tbb::parallel_for(tbb::blocked_range2d<double>(0, size1, 0, size2), [&](tbb::blocked_range2d<double> r)
		{
			for (size_t i = r.rows().begin(); i < r.rows().end(); ++i) {
				for (size_t j = r.cols().begin(); j < r.cols().end(); ++j) {
					arr_res[i][j] = arr1[i][j] + arr2[i][j];
				}
			}
		});

	for (int i = 0; i < size1; i++)
		delete[] arr_res[i];
	delete[] arr_res;
}

//Сложение oneTBB parallel_for на основе класса
class ArraySummer {

	double** arr1_;
	double** arr2_;
	double** arr_res_;

public:

	ArraySummer(double** arr1, double** arr2, double** arr_res) :arr1_(arr1), arr2_(arr2), arr_res_(arr_res) {}
	void operator()(const tbb::blocked_range2d<size_t>& r) const {
		for (size_t i = r.rows().begin(); i < r.rows().end(); i++) {
			for (size_t j = r.cols().begin(); j < r.cols().end(); j++) {
				arr_res_[i][j] = arr1_[i][j] + arr2_[i][j];
			}
		}

	}
};

inline void sum_par_for_class(double** arr1, double** arr2, int size1, int size2) {

	double** arr_res = new double* [size1];
	for (int i = 0; i < size1; i++)
		arr_res[i] = new double[size2];

	tbb::parallel_for(tbb::blocked_range2d<size_t>(0, size1, 0, size2), ArraySummer(arr1, arr2, arr_res));

	for (int i = 0; i < size1; i++)
		delete[] arr_res[i];
	delete[] arr_res;
}

//____________________________КОНЕЦ СУММЫ________________________//


//____________________________НАЧАЛО УМНОЖЕНИЯ___________________//
// 
//Последовательное перемножение
inline void mul_posled(double** arr1, double** arr2, int size1, int size2) {

	double** arr_res = new double* [size1];
	for (int i = 0; i < size1; i++)
		arr_res[i] = new double[size2];

	for (size_t i = 0; i < size1; i++) {
		for (size_t j = 0; j < size2; j++) {
			arr_res[i][j] = arr1[i][j] * arr2[i][j];
		}
	}

	for (int i = 0; i < size1; i++)
		delete[] arr_res[i];
	delete[] arr_res;

}

//Параллельное перемножение openmp parallel for
inline void mul_par_omp(double** arr1, double** arr2, int size1, int size2) {

	double** arr_res = new double* [size1];
	for (int i = 0; i < size1; i++)
		arr_res[i] = new double[size2];

#pragma omp parallel for
	for (long i = 0; i < size1; i++) {
		for (long j = 0; j < size2; j++) {
			arr_res[i][j] = arr1[i][j] * arr2[i][j];
		}
	}

	for (int i = 0; i < size1; i++)
		delete[] arr_res[i];
	delete[] arr_res;
}

//Параллельное перемножение openmp sections
inline void mul_par_omp_sections(double** arr1, double** arr2, int size1, int size2) {

	int num_thread = omp_get_max_threads();

	double** arr_res = new double* [size1];
	for (int i = 0; i < size1; i++)
		arr_res[i] = new double[size2];

	int p0 = 0,
		p1 = size1 / num_thread,
		p2 = 2 * size1 / num_thread,
		p3 = 3 * size1 / num_thread,
		p4 = size1;

#pragma omp parallel sections
	{
#pragma omp section
		{
			for (int i = p0; i < p1; i++) {
				for (int j = 0; j < size1; j++) {
					arr_res[i][j] = arr1[i][j] * arr2[i][j];
				}
			}
		}
#pragma omp section
		{
			for (int i = p1; i < p2; i++) {
				for (int j = 0; j < size1; j++) {
					arr_res[i][j] = arr1[i][j] * arr2[i][j];
				}
			}
		}
#pragma omp section
		{
			if (num_thread > 2)
				for (int i = p2; i < p3; i++)
					for (int j = 0; j < size1; j++) {
						arr_res[i][j] = arr1[i][j] * arr2[i][j];
					}
		}
#pragma omp section
		{
			if (num_thread > 3)
				for (int i = p3; i < p4; i++)
					for (int j = 0; j < size1; j++) {
						arr_res[i][j] = arr1[i][j] * arr2[i][j];
					}
		}
	}

	for (int i = 0; i < size1; i++)
		delete[] arr_res[i];
	delete[] arr_res;
}

//Параллельное перемножение oneTBB task_group
inline void mul_par_tbb_task_group(double** arr1, double** arr2, int size1, int size2) {

	int num_thread = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);

	double** arr_res = new double* [size1];
	for (int i = 0; i < size1; i++)
		arr_res[i] = new double[size2];

	int p0 = 0,
		p1 = size1 / num_thread,
		p2 = 2 * size1 / num_thread,
		p3 = 3 * size1 / num_thread,
		p4 = size1;

tbb:task_group g;

	g.run([&] {	for (long i = p0; i < p1; i++) {
		for (int j = 0; j < size1; j++) {
			arr_res[i][j] = arr1[i][j] * arr2[i][j];
		}
	}; });

	g.run([&] {	for (long i = p1; i < p2; i++) {
		for (int j = 0; j < size1; j++) {
			arr_res[i][j] = arr1[i][j] * arr2[i][j];
		}
	}; });
	if (num_thread > 2)
		g.run([&] {	for (long i = p2; i < p3; i++) {
		for (int j = 0; j < size1; j++) {
			arr_res[i][j] = arr1[i][j] * arr2[i][j];
		}
	}; });
	if (num_thread > 3)
		g.run([&] {	for (long i = p3; i < p4; i++) {
		for (int j = 0; j < size1; j++) {
			arr_res[i][j] = arr1[i][j] * arr2[i][j];
		}
	}; });

	g.wait();

	for (int i = 0; i < size1; i++)
		delete[] arr_res[i];
	delete[] arr_res;
}

//Параллельное перемножение oneTBB parallel_for лямбда выражение
inline void mul_tbb_par_for_lambda(double** arr1, double** arr2, int size1, int size2) {

	double** arr_res = new double* [size1];
	for (int i = 0; i < size1; i++)
		arr_res[i] = new double[size2];

	tbb::parallel_for(tbb::blocked_range2d<double>(0, size1, 0, size2), [&](tbb::blocked_range2d<double> r)
		{
			for (size_t i = r.rows().begin(); i < r.rows().end(); ++i) {
				for (size_t j = r.cols().begin(); j < r.cols().end(); ++j) {
					arr_res[i][j] = arr1[i][j] * arr2[i][j];
				}
			}
		});


	for (int i = 0; i < size1; i++)
		delete[] arr_res[i];
	delete[] arr_res;
}

//Параллельное перемножение oneTBB parallel_for на основе класса
class ArrayMul {

	double** arr1_;
	double** arr2_;
	double** arr_res_;

public:
	ArrayMul(double** arr1, double** arr2, double** arr_res) :arr1_(arr1), arr2_(arr2), arr_res_(arr_res) {}
	void operator()(const tbb::blocked_range2d<size_t>& r) const {
		for (size_t i = r.rows().begin(); i < r.rows().end(); i++) {
			for (size_t j = r.cols().begin(); j < r.cols().end(); j++) {
				arr_res_[i][j] = arr1_[i][j] * arr2_[i][j];
			}
		}

	}
};

inline void mul_par_for_class(double** arr1, double** arr2, int size1, int size2) {

	double** arr_res = new double* [size1];
	for (int i = 0; i < size1; i++)
		arr_res[i] = new double[size2];

	tbb::parallel_for(tbb::blocked_range2d<size_t>(0, size1, 0, size2), ArrayMul(arr1, arr2, arr_res));

	for (int i = 0; i < size1; i++)
		delete[] arr_res[i];
	delete[] arr_res;
}

//______________________________КОНЕЦ УМНОЖЕНИЯ________________//
//______________________________НАЧАЛО СУММЫ ЭЛЕМЕНТОВ_________//

//Последовательное сложение всех элементов
inline void sum_all_posled(double** arr1, double** arr2, int size1, int size2) {

	double sum = 0;
	for (size_t i = 0; i < size1; i++) {
		for (size_t j = 0; j < size2; j++) {
			sum += arr1[i][j] + arr2[i][j];
		}
	}
}

//Параллельное сложение всех элементов openmp parallel for
inline void sum_all_par_omp(double** arr1, double** arr2, int size1, int size2) {

	double sum = 0;
#pragma omp parallel for reduction(+:sum)
	for (long i = 0; i < size1; i++) {
		for (long j = 0; j < size2; j++) {
			sum += arr1[i][j] + arr2[i][j];
		}
	}

}

//Параллельное сложение всех элементов openmp sections
inline void sum_all_par_omp_sections(double** arr1, double** arr2, int size1, int size2) {

	int num_thread = omp_get_max_threads();

	double s1 = 0, s2 = 0, s3 = 0, s4 = 0, sum = 0;
	int p0 = 0,
		p1 = size1 / num_thread,
		p2 = 2 * size1 / num_thread,
		p3 = 3 * size1 / num_thread,
		p4 = size1;

#pragma omp parallel sections
	{
#pragma omp section
		{
			for (int i = p0; i < p1; i++) {
				for (int j = 0; j < size1; j++) {
					s1 += arr1[i][j] + arr2[i][j];
				}
			}
		}
#pragma omp section
		{
			for (int i = p1; i < p2; i++) {
				for (int j = 0; j < size1; j++) {
					s2 += arr1[i][j] + arr2[i][j];
				}
			}
		}
#pragma omp section
		{
			if (num_thread > 2)
				for (int i = p2; i < p3; i++)
					for (int j = 0; j < size1; j++) {
						s3 += arr1[i][j] + arr2[i][j];
					}
		}
#pragma omp section
		{
			if (num_thread > 3)
				for (int i = p3; i < p4; i++)
					for (int j = 0; j < size1; j++) {
						s4 += arr1[i][j] + arr2[i][j];
					}
		}
	}
	sum = s1 + s2 + s3 + s4;
}

//Параллельное сложение всех элементов oneTBB task_group
inline void sum_all_par_tbb_task_group(double** arr1, double** arr2, int size1, int size2) {

	int num_thread = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);

	double s1 = 0, s2 = 0, s3 = 0, s4 = 0, sum = 0;
	int p0 = 0,
		p1 = size1 / num_thread,
		p2 = 2 * size1 / num_thread,
		p3 = 3 * size1 / num_thread,
		p4 = size1;

tbb:task_group g;
	g.run([&] {	for (long i = p0; i < p1; i++) {
		for (int j = 0; j < size1; j++) {
			s1 += arr1[i][j] + arr2[i][j];
		}
	}; });

	g.run([&] {	for (long i = p1; i < p2; i++) {
		for (int j = 0; j < size1; j++) {
			s2 += arr1[i][j] + arr2[i][j];
		}
	}; });
	if (num_thread > 2)
		g.run([&] {	for (long i = p2; i < p3; i++) {
		for (int j = 0; j < size1; j++) {
			s3 += arr1[i][j] + arr2[i][j];
		}
	}; });
	if (num_thread > 3)
		g.run([&] {	for (long i = p3; i < p4; i++) {
		for (int j = 0; j < size1; j++) {
			s4 += arr1[i][j] + arr2[i][j];
		}
	}; });

	g.wait();

	sum = s1 + s2 + s3 + s4;
}

//Параллельное сложение всех элементов oneTBB parallel_reduce лямбда выражение
inline void sum_all_tbb_par_reduce_lambda(double** arr1, double** arr2, int size1, int size2) {

	auto total = tbb::parallel_reduce(tbb::blocked_range2d<size_t>(0, size1, 0, size2), 0.0, [&](tbb::blocked_range2d<size_t> r, double running_total)
		{
			for (size_t i = r.rows().begin(); i < r.rows().end(); ++i) {
				for (size_t j = r.cols().begin(); j < r.cols().end(); ++j) {
					running_total += arr1[i][j] + arr2[i][j];
				}
			}
			return running_total;
		}, std::plus<double>());
}

//Параллельное сложение всех элементов oneTBB parallel_reduce на основе класс
class reduce_par {
public:

	double sum;
	void operator()(const tbb::blocked_range2d<size_t>& r) {
		double sum_local = sum;
		double** arr1loc = arr1_;
		double** arr2loc = arr2_;

		for (size_t i = r.rows().begin(); i < r.rows().end(); i++) {
			for (size_t j = r.cols().begin(); j < r.cols().end(); j++) {
				sum_local += arr1loc[i][j] + arr2loc[i][j];
			}
		}
		sum = sum_local;
	}
	reduce_par(reduce_par& r, tbb::split) : sum(0.0), arr1_(r.arr1_), arr2_(r.arr2_) {}

	void join(const reduce_par& r) { sum += r.sum; }

	reduce_par(double** arr1, double** arr2) : sum(0.0), arr1_(arr1), arr2_(arr2) {}

private:
	double** arr1_;
	double** arr2_;
};

inline void sum_all_par_reduce_class(double** arr1, double** arr2, int size1, int size2) {

	reduce_par r(arr1, arr2);
	parallel_reduce(tbb::blocked_range2d<size_t>(0, size1, 0, size2), r);
}

//__________________КОНЕЦ СУММЫ ЭЛЕМЕНТОВ___________//
//__________________НАЧАЛО СУММЫ МАКС И МИН ЭЛЕМЕНТА_______//

//Последовательное нахождение min/max
inline void min_max_posled(double** arr1, double** arr2, int size1, int size2) {

	double MaxValue = -DBL_MAX;
	double MinValue = DBL_MAX;

	for (size_t i = 0; i < size1; i++) {
		for (size_t j = 0; j < size2; j++) {
			if (arr1[i][j] > MaxValue) MaxValue = arr1[i][j];
			if (arr2[i][j] > MaxValue) MaxValue = arr2[i][j];

			if (arr1[i][j] < MinValue) MinValue = arr1[i][j];
			if (arr2[i][j] < MinValue) MinValue = arr2[i][j];
		}
	}
}

//Параллельное нахождение min/max openmp parallel for
inline void min_max_par_omp(double** arr1, double** arr2, int size1, int size2) {

	double MaxValue = -DBL_MAX;
	double MinValue = DBL_MAX;

#pragma omp parallel for reduction(max: MaxValue) reduction(min: MinValue)
	for (int i = 0; i < size1; i++) {
		for (int j = 0; j < size2; j++) {
			if (arr1[i][j] > MaxValue) MaxValue = arr1[i][j];
			if (arr2[i][j] > MaxValue) MaxValue = arr2[i][j];

			if (arr1[i][j] < MinValue) MinValue = arr1[i][j];
			if (arr2[i][j] < MinValue) MinValue = arr2[i][j];
		}
	}
}

//Параллельное нахождение min/max openmp sections
inline void min_max_par_omp_sections(double** arr1, double** arr2, int size1, int size2) {

	int num_thread = omp_get_max_threads();

	int p0 = 0,
		p1 = size1 / num_thread,
		p2 = 2 * size1 / num_thread,
		p3 = 3 * size1 / num_thread,
		p4 = size1;
	double MaxValue = -DBL_MAX;
	double MinValue = DBL_MAX;

#pragma omp parallel sections
	{
#pragma omp section
		{
			for (int i = p0; i < p1; i++) {
				for (int j = 0; j < size1; j++) {
					if (arr1[i][j] > MaxValue) MaxValue = arr1[i][j];
					if (arr2[i][j] > MaxValue) MaxValue = arr2[i][j];

					if (arr1[i][j] < MinValue) MinValue = arr1[i][j];
					if (arr2[i][j] < MinValue) MinValue = arr2[i][j];
				}
			}
		}
#pragma omp section
		{
			for (int i = p1; i < p2; i++) {
				for (int j = 0; j < size1; j++) {
					if (arr1[i][j] > MaxValue) MaxValue = arr1[i][j];
					if (arr2[i][j] > MaxValue) MaxValue = arr2[i][j];

					if (arr1[i][j] < MinValue) MinValue = arr1[i][j];
					if (arr2[i][j] < MinValue) MinValue = arr2[i][j];
				}
			}
		}
#pragma omp section
		{
			if (num_thread > 2)
				for (int i = p2; i < p3; i++)
					for (int j = 0; j < size1; j++) {
						if (arr1[i][j] > MaxValue) MaxValue = arr1[i][j];
						if (arr2[i][j] > MaxValue) MaxValue = arr2[i][j];

						if (arr1[i][j] < MinValue) MinValue = arr1[i][j];
						if (arr2[i][j] < MinValue) MinValue = arr2[i][j];
					}
		}
#pragma omp section
		{
			if (num_thread > 3)
				for (int i = p3; i < p4; i++)
					for (int j = 0; j < size1; j++) {
						if (arr1[i][j] > MaxValue) MaxValue = arr1[i][j];
						if (arr2[i][j] > MaxValue) MaxValue = arr2[i][j];

						if (arr1[i][j] < MinValue) MinValue = arr1[i][j];
						if (arr2[i][j] < MinValue) MinValue = arr2[i][j];
					}
		}
	}
}

//Параллельное нахождение min/max oneTBB task_group
inline void min_max_par_tbb_task_group(double** arr1, double** arr2, int size1, int size2) {

	int num_thread = tbb::global_control::active_value(tbb::global_control::max_allowed_parallelism);
	int p0 = 0,
		p1 = size1 / num_thread,
		p2 = 2 * size1 / num_thread,
		p3 = 3 * size1 / num_thread,
		p4 = size1;
	double MaxValue = -DBL_MAX;
	double MinValue = DBL_MAX;
tbb:task_group g;

	g.run([&] {	for (long i = p0; i < p1; i++) {
		for (int j = 0; j < size1; j++) {
			if (arr1[i][j] > MaxValue) MaxValue = arr1[i][j];
			if (arr2[i][j] > MaxValue) MaxValue = arr2[i][j];

			if (arr1[i][j] < MinValue) MinValue = arr1[i][j];
			if (arr2[i][j] < MinValue) MinValue = arr2[i][j];
		}
	}; });

	g.run([&] {	for (long i = p1; i < p2; i++) {
		for (int j = 0; j < size1; j++) {
			if (arr1[i][j] > MaxValue) MaxValue = arr1[i][j];
			if (arr2[i][j] > MaxValue) MaxValue = arr2[i][j];

			if (arr1[i][j] < MinValue) MinValue = arr1[i][j];
			if (arr2[i][j] < MinValue) MinValue = arr2[i][j];
		}
	}; });
	if (num_thread > 2)
		g.run([&] {	for (long i = p2; i < p3; i++) {
		for (int j = 0; j < size1; j++) {
			if (arr1[i][j] > MaxValue) MaxValue = arr1[i][j];
			if (arr2[i][j] > MaxValue) MaxValue = arr2[i][j];

			if (arr1[i][j] < MinValue) MinValue = arr1[i][j];
			if (arr2[i][j] < MinValue) MinValue = arr2[i][j];
		}
	}; });
	if (num_thread > 3)
		g.run([&] {	for (long i = p3; i < p4; i++) {
		for (int j = 0; j < size1; j++) {
			if (arr1[i][j] > MaxValue) MaxValue = arr1[i][j];
			if (arr2[i][j] > MaxValue) MaxValue = arr2[i][j];

			if (arr1[i][j] < MinValue) MinValue = arr1[i][j];
			if (arr2[i][j] < MinValue) MinValue = arr2[i][j];
		}
	}; });

	g.wait();
}

//Параллельное нахождение min/max oneTBB parallel_for лямбда выражение

struct minmax_st
{
	double v_min;
	double v_max;
	minmax_st() { v_min = DBL_MAX; v_max = -DBL_MAX; };
	minmax_st(double val) :v_min(val), v_max(val) {};
};

struct minmax_join {
	minmax_st operator()(const minmax_st& _Left, const minmax_st&
		_Right) const {
		minmax_st tmp(0);
		tmp.v_max = max(_Left.v_max, _Right.v_max);
		tmp.v_min = min(_Left.v_min, _Right.v_min);
		return tmp;
	}
};

inline void min_max_tbb_par_for_lambda(double** arr1, double** arr2, int size1, int size2) {

	minmax_st total_mm = tbb::parallel_reduce(tbb::blocked_range2d<double>(0, size1, 0, size2), minmax_st(), [&](tbb::blocked_range2d<double> r, minmax_st running_maxmin)
		{
			for (size_t i = r.rows().begin(); i < r.rows().end(); ++i) {
				for (size_t j = r.cols().begin(); j < r.cols().end(); ++j) {
					running_maxmin.v_max = max(running_maxmin.v_max, *arr1[i]);
					running_maxmin.v_max = max(running_maxmin.v_max, *arr2[i]);

					running_maxmin.v_min = min(running_maxmin.v_min, *arr1[i]);
					running_maxmin.v_min = min(running_maxmin.v_min, *arr2[i]);
				}
			}
			return running_maxmin;
		}, minmax_join());
}

//Параллельное нахождение min/max на основе класса
class MinMaxCalc {

private:
	double** arr1_;
	double** arr2_;

public:
	double MinValue;
	double MaxValue;
	void operator()(const tbb::blocked_range2d<size_t>& r) {
		for (size_t i = r.rows().begin(); i < r.rows().end(); i++) {
			for (size_t j = r.cols().begin(); j < r.cols().end(); j++) {
				if (arr1_[i][j] > MaxValue) MaxValue = arr1_[i][j];
				if (arr2_[i][j] > MaxValue) MaxValue = arr2_[i][j];

				if (arr1_[i][j] < MinValue) MinValue = arr1_[i][j];
				if (arr2_[i][j] < MinValue) MinValue = arr2_[i][j];
			}
		}
	}
	MinMaxCalc(MinMaxCalc& x, tbb::split) :arr1_(x.arr1_), arr2_(x.arr2_), MinValue(DBL_MAX), MaxValue(-DBL_MAX) {}

	void join(const MinMaxCalc& y) {
		if (y.MinValue < MinValue) MinValue = y.MinValue;
		if (y.MaxValue > MaxValue) MaxValue = y.MaxValue;
	}

	MinMaxCalc(double** arr1, double** arr2) :
		arr1_(arr1),
		arr2_(arr2),
		MinValue(DBL_MAX),
		MaxValue(-DBL_MAX)
	{}
};


inline void min_max_par_reduce_class(double** arr1, double** arr2, int size1, int size2) {

	MinMaxCalc MMC(arr1, arr2);
	tbb::parallel_reduce(tbb::blocked_range2d<size_t>(0, size1, 0, size2), MMC);
}

//______________________КОНЕЦ СУММ МАКС И МИН ЭЛЕМЕНТА___________________//

//__________________НАЧАЛО Медианный фильтр____________________//
template<typename T>
void quickSort(T* arr, long size) {
	long i = 0;
	long j = size - 1;
	T pivot = arr[size / 2];

	do {
		while (arr[i] < pivot)
			i++;
		while (arr[j] > pivot)
			j--;

		if (i <= j) {
			T temp = arr[i];
			arr[i] = arr[j];
			arr[j] = temp;
			i++;
			j--;
		}
	} while (i < j);

	if (j > 0)
	{
		quickSort(arr, j + 1);
	}

	if (size > i)
	{
		quickSort(arr + i, size - i);
	}
}

void median_posled(int height, int width, RGBTRIPLE** rgb_input, RGBTRIPLE** rgb_output, int ksize) {
	int size = ksize * ksize;
	int RH = ksize / 2, RW = ksize / 2;


	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			int* medmas_red = new int[size];
			int* medmas_green = new int[size];
			int* medmas_blue = new int[size];
			int masind = 0;
			for (int dy = -RH; dy <= RH; dy++) {
				int ky = y + dy;
				if (ky < 0) ky = 0;
				if (ky > height - 1) ky = height - 1;

				for (int dx = -RW; dx <= RW; dx++) {
					int kx = x + dx;
					if (kx < 0) kx = 0;
					if (kx > width - 1) kx = width - 1;

					medmas_red[masind] = rgb_input[ky][kx].rgbtRed;
					medmas_blue[masind] = rgb_input[ky][kx].rgbtBlue;
					medmas_green[masind] = rgb_input[ky][kx].rgbtGreen;
					masind++;
				}
			}
			quickSort(medmas_red, size);
			rgb_output[y][x].rgbtRed = medmas_red[size / 2];
			quickSort(medmas_green, size);
			rgb_output[y][x].rgbtGreen = medmas_green[size / 2];
			quickSort(medmas_blue, size);
			rgb_output[y][x].rgbtBlue = medmas_blue[size / 2];

			delete[] medmas_blue;
			delete[] medmas_red;
			delete[] medmas_green;
		}
	}
}

void median_parallel_omp(int height, int width, RGBTRIPLE** rgb_input, RGBTRIPLE** rgb_output, int ksize) {
	int size = ksize * ksize;
	int RH = ksize / 2, RW = ksize / 2;

#pragma omp parallel for
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			int* medmas_red = new int[size];
			int* medmas_green = new int[size];
			int* medmas_blue = new int[size];
			int masind = 0;
			for (int dy = -RH; dy <= RH; dy++) {
				int ky = y + dy;
				if (ky < 0) ky = 0;
				if (ky > height - 1) ky = height - 1;
				for (int dx = -RW; dx <= RW; dx++) {
					int kx = x + dx;
					if (kx < 0) kx = 0;
					if (kx > width - 1) kx = width - 1;

					medmas_red[masind] = rgb_input[ky][kx].rgbtRed;
					medmas_blue[masind] = rgb_input[ky][kx].rgbtBlue;
					medmas_green[masind] = rgb_input[ky][kx].rgbtGreen;
					masind++;
				}
			}

			quickSort<int>(medmas_red, size);
			rgb_output[y][x].rgbtRed = medmas_red[size / 2];
			quickSort<int>(medmas_green, size);
			rgb_output[y][x].rgbtGreen = medmas_green[size / 2];
			quickSort<int>(medmas_blue, size);
			rgb_output[y][x].rgbtBlue = medmas_blue[size / 2];

			delete[] medmas_blue;
			delete[] medmas_red;
			delete[] medmas_green;
		}
	}
}

void median_parallel_tbb(int height, int width, RGBTRIPLE** rgb_input, RGBTRIPLE** rgb_output, int ksize) {
	int size = ksize * ksize;
	int RH = ksize / 2, RW = ksize / 2;

	tbb::parallel_for(
		tbb::blocked_range2d<int>(0, height, 0, width), [&](tbb::blocked_range2d<int> r) {
			for (int y = r.rows().begin(); y < r.rows().end(); y++) {
				for (int x = r.cols().begin(); x < r.cols().end(); x++) {
					int* medmas_red = new int[size];
					int* medmas_green = new int[size];
					int* medmas_blue = new int[size];
					int masind = 0;
					for (int dy = -RH; dy <= RH; dy++) {
						int ky = y + dy;
						if (ky < 0) ky = 0;
						if (ky > height - 1) ky = height - 1;
						for (int dx = -RW; dx <= RW; dx++) {
							int kx = x + dx;
							if (kx < 0) kx = 0;
							if (kx > width - 1) kx = width - 1;

							medmas_red[masind] = rgb_input[ky][kx].rgbtRed;
							medmas_blue[masind] = rgb_input[ky][kx].rgbtBlue;
							medmas_green[masind] = rgb_input[ky][kx].rgbtGreen;
							masind++;
						}
					}

					quickSort<int>(medmas_red, size);
					rgb_output[y][x].rgbtRed = medmas_red[size / 2];
					quickSort<int>(medmas_green, size);
					rgb_output[y][x].rgbtGreen = medmas_green[size / 2];
					quickSort<int>(medmas_blue, size);
					rgb_output[y][x].rgbtBlue = medmas_blue[size / 2];

					delete[] medmas_blue;
					delete[] medmas_red;
					delete[] medmas_green;
				}
			}
		}
	);
}
//__________________КОНЕЦ Медианный фильтр__________________________//
//__________________НАЧАЛО Фильтр Гаусса ___________________________//

double** gauss_matrix_posled(int size) {
	int RH = size / 2, RW = size / 2;
	double** matrix = new double* [size];
	for (int i = 0; i < size; i++)
		matrix[i] = new double[size];

	double SUM = 0;

	for (int y = -RH; y <= RH; y++) {
		for (int x = -RW; x <= RW; x++) {
			int YK = y + RH;
			int XK = x + RW;

			double CF = (1 / (2 * 3.14 * pow(0.8, 2))) * exp(-1 * (pow(x, 2) + pow(y, 2)) / (2 * pow(0.8, 2)));
			matrix[YK][XK] = CF;
			SUM += CF;
		}
	}

	for (int y = -RH; y <= RH; y++) {
		for (int x = -RW; x <= RW; x++) {
			int YK = y + RH;
			int XK = x + RW;
			matrix[YK][XK] /= SUM;
		}
	}

	return matrix;
}


void gauss_filter_posled(int height, int width, RGBTRIPLE** rgb_input, RGBTRIPLE** rgb_output, int ksize) {
	int size = ksize * ksize;
	int RH = ksize / 2, RW = ksize / 2;
	double** matrix = gauss_matrix_posled(ksize);

	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			int LinF_Value_Red = 0;
			int LinF_Value_Blue = 0;
			int LinF_Value_Green = 0;

			for (int dy = -RH; dy <= RH; dy++) {
				int ky = y + dy;
				if (ky < 0) ky = 0;
				if (ky > height - 1) ky = height - 1;

				for (int dx = -RW; dx <= RW; dx++) {
					int kx = x + dx;
					if (kx < 0) kx = 0;
					if (kx > width - 1) kx = width - 1;
					LinF_Value_Red += rgb_input[ky][kx].rgbtRed * matrix[dy + RH][dx + RW];
					LinF_Value_Blue += rgb_input[ky][kx].rgbtBlue * matrix[dy + RH][dx + RW];
					LinF_Value_Green += rgb_input[ky][kx].rgbtGreen * matrix[dy + RH][dx + RW];
				}
			}

			LinF_Value_Red = LinF_Value_Red < 0 ? 0 : LinF_Value_Red;
			LinF_Value_Red = LinF_Value_Red > 255 ? 255 : LinF_Value_Red;
			rgb_output[y][x].rgbtRed = LinF_Value_Red;

			LinF_Value_Green = LinF_Value_Green < 0 ? 0 : LinF_Value_Green;
			LinF_Value_Green = LinF_Value_Green > 255 ? 255 : LinF_Value_Green;
			rgb_output[y][x].rgbtGreen = LinF_Value_Green;

			LinF_Value_Blue = LinF_Value_Blue < 0 ? 0 : LinF_Value_Blue;
			LinF_Value_Blue = LinF_Value_Blue > 255 ? 255 : LinF_Value_Blue;
			rgb_output[y][x].rgbtBlue = LinF_Value_Blue;

		}
	}
}

// gauss omp

double** gauss_matrix_parallel_omp(int size) {
	int RH = size / 2, RW = size / 2;

	double** matrix = new double* [size];
	for (int i = 0; i < size; i++)
		matrix[i] = new double[size];

	double SUM = 0;

#pragma omp parallel for reduction(+:SUM)
	for (int y = -RH; y <= RH; y++) {
		for (int x = -RW; x <= RW; x++) {
			int YK = y + RH;
			int XK = x + RW;

			double CF = (1 / (2 * 3.14 * pow(RH, 2))) * exp(-1 * (pow(x, 2) + pow(y, 2)) / (2 * pow(RH, 2)));
			matrix[YK][XK] = CF;
			SUM += CF;
		}
	}

#pragma omp parallel for
	for (int y = -RH; y <= RH; y++) {
		for (int x = -RW; x <= RW; x++) {
			int YK = y + RH;
			int XK = x + RW;
			matrix[YK][XK] /= SUM;
		}
	}

	return matrix;
}

void gauss_filter_parallel_omp(int height, int width, RGBTRIPLE** rgb_input, RGBTRIPLE** rgb_output, int ksize) {
	int size = ksize * ksize;
	int RH = ksize / 2, RW = ksize / 2;

	double** matrix = gauss_matrix_parallel_omp(ksize);

#pragma omp parallel for
	for (int y = 0; y < height; y++) {
		for (int x = 0; x < width; x++) {
			int LinF_Value_Red = 0;
			int LinF_Value_Blue = 0;
			int LinF_Value_Green = 0;

			for (int dy = -RH; dy <= RH; dy++) {
				int ky = y + dy;
				if (ky < 0) ky = 0;
				if (ky > height - 1) ky = height - 1;

				for (int dx = -RW; dx <= RW; dx++) {
					int kx = x + dx;
					if (kx < 0) kx = 0;
					if (kx > width - 1) kx = width - 1;
					LinF_Value_Red += rgb_input[ky][kx].rgbtRed * matrix[dy + RH][dx + RW];
					LinF_Value_Blue += rgb_input[ky][kx].rgbtBlue * matrix[dy + RH][dx + RW];
					LinF_Value_Green += rgb_input[ky][kx].rgbtGreen * matrix[dy + RH][dx + RW];
				}
			}

			LinF_Value_Red = LinF_Value_Red < 0 ? 0 : LinF_Value_Red;
			LinF_Value_Red = LinF_Value_Red > 255 ? 255 : LinF_Value_Red;
			rgb_output[y][x].rgbtRed = LinF_Value_Red;

			LinF_Value_Green = LinF_Value_Green < 0 ? 0 : LinF_Value_Green;
			LinF_Value_Green = LinF_Value_Green > 255 ? 255 : LinF_Value_Green;
			rgb_output[y][x].rgbtGreen = LinF_Value_Green;

			LinF_Value_Blue = LinF_Value_Blue < 0 ? 0 : LinF_Value_Blue;
			LinF_Value_Blue = LinF_Value_Blue > 255 ? 255 : LinF_Value_Blue;
			rgb_output[y][x].rgbtBlue = LinF_Value_Blue;
		}
	}
}

//gauss_tbb

double** gauss_matrix_parallel_tbb(int size) {
	int RH = size / 2, RW = size / 2;

	double** matrix = new double* [size];
	for (int i = 0; i < size; i++)
		matrix[i] = new double[size];

	double SUM = 0;

tbb:task_group g;
	g.run([&] {	for (long Y = -RH; Y <= RH; Y++) {
		for (int X = -RW; X <= RW; X++) {
			int YK = Y + RH;
			int XK = X + RW;

			double CF = (1 / (2 * 3.14 * pow(RH, 2))) * exp(-1 * (pow(X, 2) + pow(Y, 2)) / (2 * pow(RH, 2)));
			matrix[YK][XK] = CF;
			SUM += CF;
		}
	};
		});


	g.wait();

	tbb::parallel_for(tbb::blocked_range2d<double>(-RH, RH, -RW, RW), [&](tbb::blocked_range2d<double> r) {
		for (size_t Y = r.rows().begin(); Y <= r.rows().end(); Y++) {
			for (size_t X = r.cols().begin(); X <= r.cols().end(); X++) {
				int YK = Y + RH;
				int XK = X + RW;
				matrix[YK][XK] /= SUM;
			}
		}
		});

	return matrix;
}


void gauss_filter_parallel_tbb(int height, int width, RGBTRIPLE** rgb_input, RGBTRIPLE** rgb_output, int ksize) {
	int size = ksize * ksize;
	int RH = ksize / 2, RW = ksize / 2;

	double** matrix = gauss_matrix_parallel_tbb(ksize);

	tbb::parallel_for(
		tbb::blocked_range2d<int>(0, height, 0, width), [&](tbb::blocked_range2d<int> r) {
			for (int y = r.rows().begin(); y < r.rows().end(); y++) {
				for (int x = r.cols().begin(); x < r.cols().end(); x++) {
					int LinF_Value_Red = 0;
					int LinF_Value_Blue = 0;
					int LinF_Value_Green = 0;

					for (int dy = -RH; dy <= RH; dy++) {
						int ky = y + dy;
						if (ky < 0) ky = 0;
						if (ky > height - 1) ky = height - 1;

						for (int dx = -RW; dx <= RW; dx++) {
							int kx = x + dx;
							if (kx < 0) kx = 0;
							if (kx > width - 1) kx = width - 1;
							LinF_Value_Red += rgb_input[ky][kx].rgbtRed * matrix[dy + RH][dx + RW];
							LinF_Value_Blue += rgb_input[ky][kx].rgbtBlue * matrix[dy + RH][dx + RW];
							LinF_Value_Green += rgb_input[ky][kx].rgbtGreen * matrix[dy + RH][dx + RW];
						}
					}

					LinF_Value_Red = LinF_Value_Red < 0 ? 0 : LinF_Value_Red;
					LinF_Value_Red = LinF_Value_Red > 255 ? 255 : LinF_Value_Red;
					rgb_output[y][x].rgbtRed = LinF_Value_Red;

					LinF_Value_Green = LinF_Value_Green < 0 ? 0 : LinF_Value_Green;
					LinF_Value_Green = LinF_Value_Green > 255 ? 255 : LinF_Value_Green;
					rgb_output[y][x].rgbtGreen = LinF_Value_Green;

					LinF_Value_Blue = LinF_Value_Blue < 0 ? 0 : LinF_Value_Blue;
					LinF_Value_Blue = LinF_Value_Blue > 255 ? 255 : LinF_Value_Blue;
					rgb_output[y][x].rgbtBlue = LinF_Value_Blue;
				}
			}
		}
	);
}

//__________________КОНЕЦ Фильтр Гаусса ___________________________//

void((*functions[24]))(double**, double**, int, int) = {
	sum_posled,
	sum_par_omp,
	sum_par_omp_sections,
	sum_par_tbb_task_group,
	sum_par_tbb_par_for_lambda,
	sum_par_for_class,

	mul_posled,
	mul_par_omp,
	mul_par_omp_sections,
	mul_par_tbb_task_group,
	mul_tbb_par_for_lambda,
	mul_par_for_class,

	sum_all_posled,
	sum_all_par_omp,
	sum_all_par_omp_sections,
	sum_all_par_tbb_task_group,
	sum_all_tbb_par_reduce_lambda,
	sum_all_par_reduce_class,

	min_max_posled,
	min_max_par_omp,
	min_max_par_omp_sections,
	min_max_par_tbb_task_group,
	min_max_tbb_par_for_lambda,
	min_max_par_reduce_class
};

const char* functions_name[24] = {
	"Сложение последовательно",
	"Сложение omp for",
	"Сложение omp sections",
	"Сложение tbb task group",
	"Сложение tbb for lambda",
	"Сложение tbb for class",

	"Умножение последовательно",
	"Умножение omp for",
	"Умножение omp sections",
	"Умножение tbb task group",
	"Умножение tbb for lambda",
	"Умножение tbb for class",

	"Сумма последовательно",
	"Сумма omp for",
	"Сумма omp sections",
	"Сумма tbb task group",
	"Сумма tbb for lambda",
	"Сумма tbb reduce class",

	"МинМакс последовательно",
	"МинМакс omp for",
	"МинМакс omp sections",
	"МинМакс tbb task group",
	"МинМакс tbb for lambda",
	"МинМакс tbb reduce class",
};

const char* images[3] = {
	"1280x720.bmp",
	"1600x900.bmp",
	"2580x1080.bmp",
};

const char* images_out[3] = {
	"1280x720_out.bmp",
	"1600x900_out.bmp",
	"2580x1080_out.bmp",
};

void((*filters[6]))(int, int, RGBTRIPLE**, RGBTRIPLE**, int) = {
	gauss_filter_posled, gauss_filter_parallel_omp, gauss_filter_parallel_tbb, median_posled, median_parallel_omp, median_parallel_tbb
};

const char* NAMES[6] = {
	"Гаусс последовательно",
	"Гаусс OMP",
	"Гаусс for TBB",

	"Медианный последовательно",
	"Медианный OMP",
	"Медианный for TBB"
};

void((*texture[3]))(RGBTRIPLE**&, BITMAPFILEHEADER&, BITMAPINFOHEADER&, const char*, int) = {
	texture_posled, texture_omp, texture_tbb
};

const char* texture_name[3] = {
"Texture posled",
"Texture omp",
"Texture tbb",
};

const char* images_n3[2] = {
"1600x900.bmp",
"2580x1080.bmp",
};

void sort(double* A, int s) {
	for (int i = 0; i < s; i++) {
		for (int j = 0; j < s - i - 1; j++) {
			if (A[j + 1] < A[j]) {
				swap(A[j], A[j + 1]);
			}
		}
	}
}

//Задание 1
/*
int main() {
	setlocale(LC_ALL, "en_US.UTF-8");
	ofstream n1("n1.txt");

	for (int function = 0; function <= 23; function++) {
		cout << functions_name[function] << endl;
		n1 << functions_name[function] << endl;
		for (int th = 2; th <= 4; th++) {
			cout << "Thread-" << th << ": ";
			n1 << "П-" << th << ": ";
			for (int sizeA = 3000; sizeA <= 6000; sizeA += 1000) {
				omp_set_num_threads(th);
				tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, th);

				double** A;
				A = new double* [sizeA];
				A[0] = new double[sizeA * sizeA];

				for (long i = 1; i < sizeA; i++)
					A[i] = &A[0][i * sizeA];

				for (long i = 0; i < sizeA; i++)
					for (int j = 0; j < sizeA; j++)
						A[i][j] = sin(i) + cos(j / 2.3) * 5;

				double** B;
				B = new double* [sizeA];
				B[0] = new double[sizeA * sizeA];

				for (int i = 1; i < sizeA; i++)
					B[i] = &B[0][i * sizeA];

				for (long i = 0; i < sizeA; i++)
					for (int j = 0; j < sizeA; j++)
						B[i][j] = sin(i) + cos(j / 2.3) * 3;

				double* median = new double[10]; // 100

				for (int run = 0; run < 10; run++) { // 100
					double start_time = omp_get_wtime();
					functions[function](A, B, sizeA, sizeA);
					median[run] = omp_get_wtime() - start_time;
				}

				sort(median, 10); // 100
				cout << median[4] * 1000 << endl; // 49
				n1 << median[4] * 1000 << " "; // 49

				delete A[0];
				delete[]A;

				delete B[0];
				delete[]B;
			}

			cout << endl;
			n1 << endl;

			if (function == 0 || function == 6 || function == 12 || function == 18)
				break;
		}
		cout << endl;
		n1 << endl;
	}

	return 0;
}
*/

//Задание 2
/*
int main() {
	setlocale(LC_ALL, "en_US.UTF-8");
	ofstream n2("n2.txt");
	for (int function = 0; function < 6; function++) {
		cout << NAMES[function] << endl;
		n2 << NAMES[function] << endl;
		for (int th = 2; th <= 4; th++) {
			for (int k = 5; k <= 9; k += 2) {
				cout << "Thread-" << th << " K" << k << ": ";
				n2 << "П-" << th << " K" << k << ": ";
				for (int image = 0; image < 3; image++) {
					omp_set_num_threads(th);
					tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, th);
					double* median = new double[1]; // 20

					RGBTRIPLE** rgb_input, ** rgb_output;
					BITMAPFILEHEADER header;
					BITMAPINFOHEADER bmiHeader;

					BMPRead(rgb_input, header, bmiHeader, images[image]);


					rgb_output = new RGBTRIPLE * [bmiHeader.biHeight];
					for (int i = 0; i < bmiHeader.biHeight; i++) {
						rgb_output[i] = new RGBTRIPLE[bmiHeader.biWidth];
					}


					for (int run = 0; run < 1; run++) { // 20
						double start_time = omp_get_wtime();
						filters[function](bmiHeader.biHeight, bmiHeader.biWidth, rgb_input, rgb_output, k);
						median[run] = omp_get_wtime() - start_time;
					}

					BMPWrite(rgb_output, bmiHeader.biWidth, bmiHeader.biHeight, images_out[image]);

					sort(median, 1); // 20

					cout << median[0] * 1000 << " "; // 9
					n2 << median[0] * 1000 << " "; // 9
				}
				cout << endl;
				n2 << endl;
			}

			//Res << std::endl << std::endl;

			if (function == 0 || function == 3)
				break;
			n2 << endl << endl;
		}
	}
	return 0;
}
*/

//Задание 3
int main() {
	setlocale(LC_ALL, "en_US.UTF-8");
	ofstream n3("n3.txt");
	for (int filter = 0; filter < 3; filter++) {
		cout << texture_name[filter] << endl;
		n3 << texture_name[filter] << std::endl;
		for (int th = 2; th <= 4; th++) {
			for (int image = 0; image < 2; image++) {
				for (int k = 5; k <= 9; k += 2) {
					cout << "k = " << k << endl;
					n3 << "k = " << k << std::endl;
					for (int image = 0; image < 2; image += 1) {
						RGBTRIPLE** rgb_input;
						BITMAPFILEHEADER header;
						BITMAPINFOHEADER bmiHeader;

						omp_set_num_threads(th);
						tbb::global_control global_limit(tbb::global_control::max_allowed_parallelism, th);

						double* median = new double[1]; // 10

						for (int run = 0; run < 1; run++) { // 10
							double start_time = omp_get_wtime();
							texture[filter](rgb_input, header, bmiHeader, images_n3[image], k);
							median[run] = omp_get_wtime() - start_time;
						}

						sort(median, 1); // 10

						cout << median[0] * 1000 << " "; // 4
						n3 << median[0] * 1000 << " "; // 4

						//system("pause");
					}
					cout << std::endl;
					n3 << std::endl;
				}


			}
			if (filter == 0)
				break;
			cout << std::endl << std::endl;
			n3 << std::endl << std::endl;
		}
	}
	return 0;
}