#pragma once
#include "Header1.h"

void texture_posled(RGBTRIPLE**& rgb_input, BITMAPFILEHEADER& header, BITMAPINFOHEADER& bmiHeader, const char* in, int ksize) {
	BMPRead(rgb_input, header, bmiHeader, in);
	RGBTRIPLE** rgbout_m2 = new RGBTRIPLE * [bmiHeader.biHeight];
	RGBTRIPLE** rgbout_u = new RGBTRIPLE * [bmiHeader.biHeight];
	RGBTRIPLE** rgbout_r = new RGBTRIPLE * [bmiHeader.biHeight];
	RGBTRIPLE** rgbout_e = new RGBTRIPLE * [bmiHeader.biHeight];

	for (int i = 0; i < bmiHeader.biHeight; i++) {
		rgbout_m2[i] = new RGBTRIPLE[bmiHeader.biWidth];
		rgbout_u[i] = new RGBTRIPLE[bmiHeader.biWidth];
		rgbout_r[i] = new RGBTRIPLE[bmiHeader.biWidth];
		rgbout_e[i] = new RGBTRIPLE[bmiHeader.biWidth];
	}

	int** BIMAGE = new int* [bmiHeader.biHeight];
	double** M2 = new double* [bmiHeader.biHeight];
	double** R = new double* [bmiHeader.biHeight];
	double** U = new double* [bmiHeader.biHeight];
	double** E = new double* [bmiHeader.biHeight];
	double Hist[256];

	for (int i = 0; i < bmiHeader.biHeight; i++) {
		BIMAGE[i] = new int[bmiHeader.biWidth];
		M2[i] = new double[bmiHeader.biWidth];
		R[i] = new double[bmiHeader.biWidth];
		U[i] = new double[bmiHeader.biWidth];
		E[i] = new double[bmiHeader.biWidth];
	}

	int size = ksize * ksize;
	int rh = ksize / 2;
	int rw = ksize / 2;

	for (int Y = 0; Y < bmiHeader.biHeight; Y++) {
		for (int X = 0; X < bmiHeader.biWidth; X++) {

			BIMAGE[Y][X] = rgb_input[Y][X].rgbtRed * 0.299 + rgb_input[Y][X].rgbtGreen * 0.587 + rgb_input[Y][X].rgbtBlue * 0.114;
		}
	}
	double min_m2 = DBL_MAX, min_u = DBL_MAX, min_r = DBL_MAX, min_e = DBL_MAX;
	double max_m2 = -DBL_MAX, max_u = -DBL_MAX, max_r = -DBL_MAX, max_e = -DBL_MAX;

	for (int Y = 0; Y < bmiHeader.biHeight; Y++) {
		for (int X = 0; X < bmiHeader.biWidth; X++) {

			for (int i = 0; i <= 255; i++) {
				Hist[i] = 0;
			}

			for (int DY = -rh; DY <= rh; DY++) {
				int KY = Y + DY;
				if (KY < 0) {
					KY = 0;
				}
				if (KY > bmiHeader.biHeight - 1) {
					KY = bmiHeader.biHeight - 1;
				}
				for (int DX = -rw; DX <= rw; DX++) {
					int KX = X + DX;
					if (KX < 0)
						KX = 0;
					if (KX > bmiHeader.biWidth - 1)
						KX = bmiHeader.biWidth - 1;
					Hist[BIMAGE[KY][KX]] = Hist[BIMAGE[KY][KX]] + 1;
				}
			}

			int size = (rh * 2 + 1) * (rw * 2 + 1);
			double m = 0, m2 = 0, e = 0, u = 0;

			for (int i = 0; i <= 255; i++) {
				Hist[i] = Hist[i] / size;
			}

			for (int i = 0; i <= 255; i++) {
				m += Hist[i] * i;
			}

			for (int i = 0; i <= 255; i++) {
				m2 += pow((i - m), 2) * Hist[i];
				if (Hist[i] != 0) {
					e += Hist[i] * log2(Hist[i]);
				}
				u += pow(Hist[i], 2);
			}

			M2[Y][X] = m2;
			U[Y][X] = u;
			R[Y][X] = 1 - (1 / (1 + m2));
			E[Y][X] = -1 * e;

			if (Y == 0 && X == 0)
			{
				max_m2 = m2;
				max_u = u;
				max_r = (1 - (1 / (1 + m2)));
				max_e = (-1 * e);
				min_m2 = m2;
				min_u = u;
				min_r = (1 - (1 / (1 + m2)));
				min_e = (-1 * e);
			}
			else
			{
				if (max_m2 < m2) {
					max_m2 = m2;
				}
				if (max_u < u) {
					max_u = u;
				}
				if (max_r < (1 - (1 / (1 + m2)))) {
					max_r = (1 - (1 / (1 + m2)));
				}
				if (max_e < (-1 * e)) {
					max_e = (-1 * e);
				}
				if (min_m2 > m2) {
					min_m2 = m2;
				}
				if (min_u > u) {
					min_u = u;
				}
				if (min_r > (1 - (1 / (1 + m2)))) {
					min_r = (1 - (1 / (1 + m2)));
				}
				if (min_e > (-1 * e)) {
					min_e = (-1 * e);
				}
			}
		}
	}

	double T1_m2 = (max_m2 - min_m2) * 0.2 + min_m2;
	double T2_m2 = (max_m2 - min_m2) * 0.7 + min_m2;
	double T1_u = (max_u - min_u) * 0.2 + min_u;
	double T2_u = (max_u - min_u) * 0.7 + min_u;
	double T1_r = (max_r - min_r) * 0.2 + min_r;
	double T2_r = (max_r - min_r) * 0.7 + min_r;
	double T1_e = (max_e - min_e) * 0.2 + min_e;
	double T2_e = (max_e - min_e) * 0.7 + min_e;

	std::cout << "T1_m2 = " << T1_m2 << std::endl << "T2_m2 = " << T2_m2 << std::endl <<
		"T1_u = " << T1_u << std::endl << "T2_u = " << T2_u << std::endl << "T1_r = " << T1_r << std::endl <<
		"T2_r = " << T2_r << std::endl << "T1_e = " << T1_e << std::endl << "T2_e = " << T2_e << std::endl;

	for (int Y = 0; Y < bmiHeader.biHeight; Y++) {
		for (int X = 0; X < bmiHeader.biWidth; X++) {
			if (M2[Y][X] < T1_m2) {
				rgbout_m2[Y][X].rgbtRed = 0;
				rgbout_m2[Y][X].rgbtGreen = 255;
				rgbout_m2[Y][X].rgbtBlue = 0;
			}
			else if (M2[Y][X] >= T1_m2 && M2[Y][X] <= T2_m2) {
				rgbout_m2[Y][X].rgbtRed = 255;
				rgbout_m2[Y][X].rgbtGreen = 255;
				rgbout_m2[Y][X].rgbtBlue = 0;
			}
			else {
				rgbout_m2[Y][X].rgbtRed = 255;
				rgbout_m2[Y][X].rgbtGreen = 0;
				rgbout_m2[Y][X].rgbtBlue = 0;
			}


			if (U[Y][X] < T1_u) {
				rgbout_u[Y][X].rgbtRed = 0;
				rgbout_u[Y][X].rgbtGreen = 255;
				rgbout_u[Y][X].rgbtBlue = 0;
			}
			else if (U[Y][X] >= T1_u && U[Y][X] <= T2_u) {
				rgbout_u[Y][X].rgbtRed = 255;
				rgbout_u[Y][X].rgbtGreen = 255;
				rgbout_u[Y][X].rgbtBlue = 0;
			}
			else {
				rgbout_u[Y][X].rgbtRed = 255;
				rgbout_u[Y][X].rgbtGreen = 0;
				rgbout_u[Y][X].rgbtBlue = 0;
			}


			if (R[Y][X] < T1_r) {
				rgbout_r[Y][X].rgbtRed = 0;
				rgbout_r[Y][X].rgbtGreen = 255;
				rgbout_r[Y][X].rgbtBlue = 0;
			}
			else if (R[Y][X] >= T1_r && R[Y][X] <= T2_r) {
				rgbout_r[Y][X].rgbtRed = 255;
				rgbout_r[Y][X].rgbtGreen = 255;
				rgbout_r[Y][X].rgbtBlue = 0;
			}
			else {
				rgbout_r[Y][X].rgbtRed = 255;
				rgbout_r[Y][X].rgbtGreen = 0;
				rgbout_r[Y][X].rgbtBlue = 0;
			}


			if (E[Y][X] < T1_e) {
				rgbout_e[Y][X].rgbtRed = 0;
				rgbout_e[Y][X].rgbtGreen = 255;
				rgbout_e[Y][X].rgbtBlue = 0;
			}
			else if (E[Y][X] >= T1_e && M2[Y][X] <= T2_e) {
				rgbout_e[Y][X].rgbtRed = 255;
				rgbout_e[Y][X].rgbtGreen = 255;
				rgbout_e[Y][X].rgbtBlue = 0;
			}
			else {
				rgbout_e[Y][X].rgbtRed = 255;
				rgbout_e[Y][X].rgbtGreen = 0;
				rgbout_e[Y][X].rgbtBlue = 0;
			}
		}
	}

	BMPWrite(rgbout_m2, bmiHeader.biWidth, bmiHeader.biHeight, "posl__output_m2.bmp");
	BMPWrite(rgbout_u, bmiHeader.biWidth, bmiHeader.biHeight, "posl__output_u.bmp");
	BMPWrite(rgbout_r, bmiHeader.biWidth, bmiHeader.biHeight, "posl__output_r.bmp");
	BMPWrite(rgbout_e, bmiHeader.biWidth, bmiHeader.biHeight, "posl__output_e.bmp");

	for (int i = 0; i < bmiHeader.biHeight; i++)
	{
		delete[] rgbout_m2[i];
		delete[] rgbout_u[i];
		delete[] rgbout_r[i];
		delete[] rgbout_e[i];
		delete[] BIMAGE[i];
		delete[] M2[i];
		delete[] R[i];
		delete[] U[i];
		delete[] E[i];
	}

	delete[] rgbout_m2;
	delete[] rgbout_u;
	delete[] rgbout_r;
	delete[] rgbout_e;
	delete[] BIMAGE;
	delete[] M2;
	delete[] R;
	delete[] U;
	delete[] E;
}

void texture_omp(RGBTRIPLE**& rgb_input, BITMAPFILEHEADER& header, BITMAPINFOHEADER& bmiHeader, const char* in, int ksize) {
	BMPRead(rgb_input, header, bmiHeader, in);
	RGBTRIPLE** rgbout_m2 = new RGBTRIPLE * [bmiHeader.biHeight];
	RGBTRIPLE** rgbout_u = new RGBTRIPLE * [bmiHeader.biHeight];
	RGBTRIPLE** rgbout_r = new RGBTRIPLE * [bmiHeader.biHeight];
	RGBTRIPLE** rgbout_e = new RGBTRIPLE * [bmiHeader.biHeight];

	for (int i = 0; i < bmiHeader.biHeight; i++)
	{
		rgbout_m2[i] = new RGBTRIPLE[bmiHeader.biWidth];
		rgbout_u[i] = new RGBTRIPLE[bmiHeader.biWidth];
		rgbout_r[i] = new RGBTRIPLE[bmiHeader.biWidth];
		rgbout_e[i] = new RGBTRIPLE[bmiHeader.biWidth];
	}


	int** BIMAGE = new int* [bmiHeader.biHeight];
	double** M2 = new double* [bmiHeader.biHeight];
	double** R = new double* [bmiHeader.biHeight];
	double** U = new double* [bmiHeader.biHeight];
	double** E = new double* [bmiHeader.biHeight];
	double Hist[256];

	for (int i = 0; i < bmiHeader.biHeight; i++) {
		BIMAGE[i] = new int[bmiHeader.biWidth];
		M2[i] = new double[bmiHeader.biWidth];
		R[i] = new double[bmiHeader.biWidth];
		U[i] = new double[bmiHeader.biWidth];
		E[i] = new double[bmiHeader.biWidth];
	}

	int size = ksize * ksize;
	int RH = ksize / 2;
	int RW = ksize / 2;

#pragma omp parallel for
	for (int Y = 0; Y < bmiHeader.biHeight; Y++) {
		for (int X = 0; X < bmiHeader.biWidth; X++) {

			BIMAGE[Y][X] = rgb_input[Y][X].rgbtRed * 0.299 + rgb_input[Y][X].rgbtGreen * 0.587 + rgb_input[Y][X].rgbtBlue * 0.114;
		}
	}

	double min_m2 = DBL_MAX, min_u = DBL_MAX, min_r = DBL_MAX, min_e = DBL_MAX;
	double max_m2 = -DBL_MAX, max_u = -DBL_MAX, max_r = -DBL_MAX, max_e = -DBL_MAX;

#pragma omp parallel for reduction(max: max_m2, max_u, max_r, max_e) reduction(min: min_m2, min_u, min_r, min_e)
	for (int y = 0; y < bmiHeader.biHeight; y++) {
		for (int x = 0; x < bmiHeader.biWidth; x++) {
			double Hist[256];
			for (int i = 0; i < 256; i++)
				Hist[i] = 0;

			for (int dy = -RH; dy <= RH; dy++) {
				int ky = y + dy;
				if (ky < 0) ky = 0;
				if (ky > bmiHeader.biHeight - 1) ky = bmiHeader.biHeight - 1;

				for (int dx = -RW; dx <= RW; dx++) {
					int kx = x + dx;
					if (kx < 0) kx = 0;
					if (kx > bmiHeader.biWidth - 1) kx = bmiHeader.biWidth - 1;
					Hist[BIMAGE[ky][kx]] = Hist[BIMAGE[ky][kx]] + 1;
				}
			}

			double m = 0, m2 = 0, e = 0, u = 0;

			for (int i = 0; i <= 255; i++)
				Hist[i] = Hist[i] / size;
			for (int i = 0; i <= 255; i++)
				m += Hist[i] * i;
			for (int i = 0; i <= 255; i++) {
				m2 += pow(i - m, 2) * Hist[i];
				if (Hist[i] != 0) e += Hist[i] * log2(Hist[i]);
				u += pow(Hist[i], 2);
			}
			M2[y][x] = m2;
			U[y][x] = u;
			R[y][x] = 1 - (1 / (1 + m2));
			E[y][x] = -1 * e;

			if (y == 0 && x == 0)
			{
				max_m2 = m2;
				max_u = u;
				max_r = (1 - (1 / (1 + m2)));
				max_e = (-1 * e);
				min_m2 = m2;
				min_u = u;
				min_r = (1 - (1 / (1 + m2)));
				min_e = (-1 * e);
			}
			else
			{
				if (max_m2 < m2) {
					max_m2 = m2;
				}
				if (max_u < u) {
					max_u = u;
				}
				if (max_r < (1 - (1 / (1 + m2)))) {
					max_r = (1 - (1 / (1 + m2)));
				}
				if (max_e < (-1 * e)) {
					max_e = (-1 * e);
				}
				if (min_m2 > m2) {
					min_m2 = m2;
				}
				if (min_u > u) {
					min_u = u;
				}
				if (min_r > (1 - (1 / (1 + m2)))) {
					min_r = (1 - (1 / (1 + m2)));
				}
				if (min_e > (-1 * e)) {
					min_e = (-1 * e);
				}
			}
		}
	}

	double T1_m2 = (max_m2 - min_m2) * 0.2 + min_m2;
	double T2_m2 = (max_m2 - min_m2) * 0.7 + min_m2;
	double T1_u = (max_u - min_u) * 0.2 + min_u;
	double T2_u = (max_u - min_u) * 0.7 + min_u;
	double T1_r = (max_r - min_r) * 0.2 + min_r;
	double T2_r = (max_r - min_r) * 0.7 + min_r;
	double T1_e = (max_e - min_e) * 0.2 + min_e;
	double T2_e = (max_e - min_e) * 0.7 + min_e;


	for (int Y = 0; Y < bmiHeader.biHeight; Y++) {
		for (int X = 0; X < bmiHeader.biWidth; X++) {
			if (M2[Y][X] < T1_m2) {
				rgbout_m2[Y][X].rgbtRed = 0;
				rgbout_m2[Y][X].rgbtGreen = 255;
				rgbout_m2[Y][X].rgbtBlue = 0;
			}
			else if (M2[Y][X] >= T1_m2 && M2[Y][X] <= T2_m2) {
				rgbout_m2[Y][X].rgbtRed = 255;
				rgbout_m2[Y][X].rgbtGreen = 255;
				rgbout_m2[Y][X].rgbtBlue = 0;
			}
			else {
				rgbout_m2[Y][X].rgbtRed = 255;
				rgbout_m2[Y][X].rgbtGreen = 0;
				rgbout_m2[Y][X].rgbtBlue = 0;
			}


			if (U[Y][X] < T1_u) {
				rgbout_u[Y][X].rgbtRed = 0;
				rgbout_u[Y][X].rgbtGreen = 255;
				rgbout_u[Y][X].rgbtBlue = 0;
			}
			else if (U[Y][X] >= T1_u && U[Y][X] <= T2_u) {
				rgbout_u[Y][X].rgbtRed = 255;
				rgbout_u[Y][X].rgbtGreen = 255;
				rgbout_u[Y][X].rgbtBlue = 0;
			}
			else {
				rgbout_u[Y][X].rgbtRed = 255;
				rgbout_u[Y][X].rgbtGreen = 0;
				rgbout_u[Y][X].rgbtBlue = 0;
			}


			if (R[Y][X] < T1_r) {
				rgbout_r[Y][X].rgbtRed = 0;
				rgbout_r[Y][X].rgbtGreen = 255;
				rgbout_r[Y][X].rgbtBlue = 0;
			}
			else if (R[Y][X] >= T1_r && R[Y][X] <= T2_r) {
				rgbout_r[Y][X].rgbtRed = 255;
				rgbout_r[Y][X].rgbtGreen = 255;
				rgbout_r[Y][X].rgbtBlue = 0;
			}
			else {
				rgbout_r[Y][X].rgbtRed = 255;
				rgbout_r[Y][X].rgbtGreen = 0;
				rgbout_r[Y][X].rgbtBlue = 0;
			}


			if (E[Y][X] < T1_e) {
				rgbout_e[Y][X].rgbtRed = 0;
				rgbout_e[Y][X].rgbtGreen = 255;
				rgbout_e[Y][X].rgbtBlue = 0;
			}
			else if (E[Y][X] >= T1_e && M2[Y][X] <= T2_e) {
				rgbout_e[Y][X].rgbtRed = 255;
				rgbout_e[Y][X].rgbtGreen = 255;
				rgbout_e[Y][X].rgbtBlue = 0;
			}
			else {
				rgbout_e[Y][X].rgbtRed = 255;
				rgbout_e[Y][X].rgbtGreen = 0;
				rgbout_e[Y][X].rgbtBlue = 0;
			}
		}
	}

	BMPWrite(rgbout_m2, bmiHeader.biWidth, bmiHeader.biHeight, "omp__output_m2.bmp");
	BMPWrite(rgbout_u, bmiHeader.biWidth, bmiHeader.biHeight, "omp__output_u.bmp");
	BMPWrite(rgbout_r, bmiHeader.biWidth, bmiHeader.biHeight, "omp__output_r.bmp");
	BMPWrite(rgbout_e, bmiHeader.biWidth, bmiHeader.biHeight, "omp__output_e.bmp");

#pragma omp parallel for
	for (int i = 0; i < bmiHeader.biHeight; i++)
	{
		delete[] rgbout_m2[i];
		delete[] rgbout_u[i];
		delete[] rgbout_r[i];
		delete[] rgbout_e[i];
		delete[] BIMAGE[i];
		delete[] M2[i];
		delete[] R[i];
		delete[] U[i];
		delete[] E[i];
	}

	delete[] rgbout_m2;
	delete[] rgbout_u;
	delete[] rgbout_r;
	delete[] rgbout_e;
	delete[] BIMAGE;
	delete[] M2;
	delete[] R;
	delete[] U;
	delete[] E;
}

//Параллельное нахождение min/max reduce в texture на основе класса
class MinMaxTexture {

private:
	int RH_;
	int RW_;
	int** BIMAGE_;
	double** M2_;
	double** R_;
	double** U_;
	double** E_;

public:
	double min_m2, min_u, min_r, min_e;
	double max_m2, max_u, max_r, max_e;
	void operator()(const tbb::blocked_range2d<size_t>& r) {
		for (size_t y = r.rows().begin(); y < r.rows().end(); y++) {
			for (size_t x = r.cols().begin(); x < r.cols().end(); x++) {

				double Hist[256];
				for (int i = 0; i < 256; i++)
					Hist[i] = 0;

				for (int dy = -RH_; dy <= RH_; dy++) {
					int ky = y + dy;
					if (ky < 0) ky = 0;
					if (ky > r.rows().end() - 1) ky = r.rows().end() - 1;

					for (int dx = -RW_; dx <= RW_; dx++) {
						int kx = x + dx;
						if (kx < 0) kx = 0;
						if (kx > r.cols().end() - 1) kx = r.cols().end() - 1;
						Hist[BIMAGE_[ky][kx]] = Hist[BIMAGE_[ky][kx]] + 1;
					}
				}

				int Size = (RH_ * 2 + 1) * (RW_ * 2 + 1);
				double m = 0, m2 = 0, e = 0, u = 0;

				for (int i = 0; i <= 255; i++)
					Hist[i] = Hist[i] / Size;
				for (int i = 0; i <= 255; i++)
					m += Hist[i] * i;
				for (int i = 0; i <= 255; i++) {
					m2 += pow(i - m, 2) * Hist[i];
					if (Hist[i] != 0) e += Hist[i] * log2(Hist[i]);
					u += pow(Hist[i], 2);
				}
				M2_[y][x] = m2;
				U_[y][x] = u;
				R_[y][x] = 1 - (1 / (1 + m2));
				E_[y][x] = -1 * e;

				if (y == 0 && x == 0)
				{
					max_m2 = m2;
					max_u = u;
					max_r = (1 - (1 / (1 + m2)));
					max_e = (-1 * e);
					min_m2 = m2;
					min_u = u;
					min_r = (1 - (1 / (1 + m2)));
					min_e = (-1 * e);
				}
				else
				{
					if (max_m2 < m2) {
						max_m2 = m2;
					}
					if (max_u < u) {
						max_u = u;
					}
					if (max_r < (1 - (1 / (1 + m2)))) {
						max_r = (1 - (1 / (1 + m2)));
					}
					if (max_e < (-1 * e)) {
						max_e = (-1 * e);
					}
					if (min_m2 > m2) {
						min_m2 = m2;
					}
					if (min_u > u) {
						min_u = u;
					}
					if (min_r > (1 - (1 / (1 + m2)))) {
						min_r = (1 - (1 / (1 + m2)));
					}
					if (min_e > (-1 * e)) {
						min_e = (-1 * e);
					}
				}

			}
		}
	}
	MinMaxTexture(MinMaxTexture& x, tbb::split) :RH_(x.RH_), RW_(x.RW_), BIMAGE_(x.BIMAGE_), M2_(x.M2_), R_(x.R_), U_(x.U_), E_(x.E_),
												min_m2(DBL_MAX), min_u(DBL_MAX), min_r(DBL_MAX), min_e(DBL_MAX),
												max_m2(-DBL_MAX), max_u(-DBL_MAX), max_r(-DBL_MAX), max_e(-DBL_MAX) {}

	void join(const MinMaxTexture& y) {
		if (max_m2 < y.max_m2) max_m2 = y.max_m2;
		if (max_u < y.max_u) max_u = y.max_u;
		if (max_r < (1 - (1 / (1 + y.max_m2)))) max_r = (1 - (1 / (1 + y.max_m2)));
		if (max_e < (-1 * y.max_e)) max_e = (-1 * y.max_e);
		if (min_m2 > y.min_m2) min_m2 = y.min_m2;
		if (min_u > y.min_u) min_u = y.min_u;
		if (min_r > (1 - (1 / (1 + y.min_m2)))) min_r = (1 - (1 / (1 + y.min_m2)));
		if (min_e > (-1 * y.min_e)) min_e = (-1 * y.min_e);
	}

	MinMaxTexture(int RH, int RW, int** BIMAGE, double** M2,
			      double** R, double** U, double** E) :
		RH_(RH),
		RW_(RW),
		BIMAGE_(BIMAGE),
		M2_(M2),
		R_(R),
		U_(U),
		E_(E),
		min_m2(DBL_MAX), min_u(DBL_MAX), min_r(DBL_MAX), min_e(DBL_MAX),
		max_m2(-DBL_MAX), max_u(-DBL_MAX), max_r(-DBL_MAX), max_e(-DBL_MAX)
	{}
};

void texture_tbb(RGBTRIPLE**& rgb_input, BITMAPFILEHEADER& header, BITMAPINFOHEADER& bmiHeader, const char* in, int ksize) {
	BMPRead(rgb_input, header, bmiHeader, in);
	RGBTRIPLE** rgbout_m2 = new RGBTRIPLE * [bmiHeader.biHeight];
	RGBTRIPLE** rgbout_u = new RGBTRIPLE * [bmiHeader.biHeight];
	RGBTRIPLE** rgbout_r = new RGBTRIPLE * [bmiHeader.biHeight];
	RGBTRIPLE** rgbout_e = new RGBTRIPLE * [bmiHeader.biHeight];


	int size = ksize * ksize;
	int RH = ksize / 2;
	int RW = ksize / 2;

	for (int i = 0; i < bmiHeader.biHeight; i++) {
		rgbout_m2[i] = new RGBTRIPLE[bmiHeader.biWidth];
		rgbout_u[i] = new RGBTRIPLE[bmiHeader.biWidth];
		rgbout_r[i] = new RGBTRIPLE[bmiHeader.biWidth];
		rgbout_e[i] = new RGBTRIPLE[bmiHeader.biWidth];
	}

	int imWidth = bmiHeader.biWidth;
	int imHeight = bmiHeader.biHeight;

	int** BIMAGE = new int* [bmiHeader.biHeight];
	double** M2 = new double* [bmiHeader.biHeight];
	double** R = new double* [bmiHeader.biHeight];
	double** U = new double* [bmiHeader.biHeight];
	double** E = new double* [bmiHeader.biHeight];
	double Hist[256];

	for (int i = 0; i < bmiHeader.biHeight; i++) {
		BIMAGE[i] = new int[bmiHeader.biWidth];
		M2[i] = new double[bmiHeader.biWidth];
		R[i] = new double[bmiHeader.biWidth];
		U[i] = new double[bmiHeader.biWidth];
		E[i] = new double[bmiHeader.biWidth];
	}

	tbb::parallel_for(tbb::blocked_range2d<double>(0, bmiHeader.biHeight, 0, bmiHeader.biWidth), [&](tbb::blocked_range2d<double> r) {
		for (int Y = r.rows().begin(); Y < r.rows().end(); Y++) {
			for (int X = r.cols().begin(); X < r.cols().end(); X++) {

				BIMAGE[Y][X] = rgb_input[Y][X].rgbtRed * 0.299 + rgb_input[Y][X].rgbtGreen * 0.587 + rgb_input[Y][X].rgbtBlue * 0.114;
			}
		}
		});

	MinMaxTexture MMT(RH, RW, BIMAGE, M2, R, U, E);
	tbb::parallel_reduce(tbb::blocked_range2d<size_t>(0, imHeight, 0, imWidth), MMT);

	double min_m2 = MMT.min_m2, min_u = MMT.min_u, min_r = MMT.min_r, min_e = MMT.min_e;
	double max_m2 = MMT.max_m2, max_u = MMT.max_u, max_r = MMT.max_r, max_e = MMT.max_e;

	double T1_m2 = (max_m2 - min_m2) * 0.2 + min_m2;
	double T2_m2 = (max_m2 - min_m2) * 0.7 + min_m2;
	double T1_u = (max_u - min_u) * 0.2 + min_u;
	double T2_u = (max_u - min_u) * 0.7 + min_u;
	double T1_r = (max_r - min_r) * 0.2 + min_r;
	double T2_r = (max_r - min_r) * 0.7 + min_r;
	double T1_e = (max_e - min_e) * 0.2 + min_e;
	double T2_e = (max_e - min_e) * 0.7 + min_e;


	for (int Y = 0; Y < bmiHeader.biHeight; Y++) {
		for (int X = 0; X < bmiHeader.biWidth; X++) {
			if (M2[Y][X] < T1_m2) {
				rgbout_m2[Y][X].rgbtRed = 0;
				rgbout_m2[Y][X].rgbtGreen = 255;
				rgbout_m2[Y][X].rgbtBlue = 0;
			}
			else if (M2[Y][X] >= T1_m2 && M2[Y][X] <= T2_m2) {
				rgbout_m2[Y][X].rgbtRed = 255;
				rgbout_m2[Y][X].rgbtGreen = 255;
				rgbout_m2[Y][X].rgbtBlue = 0;
			}
			else {
				rgbout_m2[Y][X].rgbtRed = 255;
				rgbout_m2[Y][X].rgbtGreen = 0;
				rgbout_m2[Y][X].rgbtBlue = 0;
			}


			if (U[Y][X] < T1_u) {
				rgbout_u[Y][X].rgbtRed = 0;
				rgbout_u[Y][X].rgbtGreen = 255;
				rgbout_u[Y][X].rgbtBlue = 0;
			}
			else if (U[Y][X] >= T1_u && U[Y][X] <= T2_u) {
				rgbout_u[Y][X].rgbtRed = 255;
				rgbout_u[Y][X].rgbtGreen = 255;
				rgbout_u[Y][X].rgbtBlue = 0;
			}
			else {
				rgbout_u[Y][X].rgbtRed = 255;
				rgbout_u[Y][X].rgbtGreen = 0;
				rgbout_u[Y][X].rgbtBlue = 0;
			}


			if (R[Y][X] < T1_r) {
				rgbout_r[Y][X].rgbtRed = 0;
				rgbout_r[Y][X].rgbtGreen = 255;
				rgbout_r[Y][X].rgbtBlue = 0;
			}
			else if (R[Y][X] >= T1_r && R[Y][X] <= T2_r) {
				rgbout_r[Y][X].rgbtRed = 255;
				rgbout_r[Y][X].rgbtGreen = 255;
				rgbout_r[Y][X].rgbtBlue = 0;
			}
			else {
				rgbout_r[Y][X].rgbtRed = 255;
				rgbout_r[Y][X].rgbtGreen = 0;
				rgbout_r[Y][X].rgbtBlue = 0;
			}


			if (E[Y][X] < T1_e) {
				rgbout_e[Y][X].rgbtRed = 0;
				rgbout_e[Y][X].rgbtGreen = 255;
				rgbout_e[Y][X].rgbtBlue = 0;
			}
			else if (E[Y][X] >= T1_e && M2[Y][X] <= T2_e) {
				rgbout_e[Y][X].rgbtRed = 255;
				rgbout_e[Y][X].rgbtGreen = 255;
				rgbout_e[Y][X].rgbtBlue = 0;
			}
			else {
				rgbout_e[Y][X].rgbtRed = 255;
				rgbout_e[Y][X].rgbtGreen = 0;
				rgbout_e[Y][X].rgbtBlue = 0;
			}
		}
	}

	BMPWrite(rgbout_m2, bmiHeader.biWidth, bmiHeader.biHeight, "tbb__output_m2.bmp");
	BMPWrite(rgbout_u, bmiHeader.biWidth, bmiHeader.biHeight, "tbb__output_u.bmp");
	BMPWrite(rgbout_r, bmiHeader.biWidth, bmiHeader.biHeight, "tbb__output_r.bmp");
	BMPWrite(rgbout_e, bmiHeader.biWidth, bmiHeader.biHeight, "tbb__output_e.bmp");

	tbb::parallel_for(tbb::blocked_range<int>(0, bmiHeader.biHeight), [&](tbb::blocked_range<int> r) {
		for (int i = r.begin(); i < r.end(); i++)
		{
			delete[] rgbout_m2[i];
			delete[] rgbout_u[i];
			delete[] rgbout_r[i];
			delete[] rgbout_e[i];
			delete[] BIMAGE[i];
			delete[] M2[i];
			delete[] R[i];
			delete[] U[i];
			delete[] E[i];
		}
		});

	delete[] rgbout_m2;
	delete[] rgbout_u;
	delete[] rgbout_r;
	delete[] rgbout_e;
	delete[] BIMAGE;
	delete[] M2;
	delete[] R;
	delete[] U;
	delete[] E;
}
