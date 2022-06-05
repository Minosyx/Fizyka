#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define G 6.67408e-11
#define RZ 6371.008e3
#define RK 1737.064e3
#define vII 11.19e3
#define PI 3.141592654
#define N 3

float* createArray(int size) {
	float* tab;
	tab = (float*)malloc(sizeof(float) * size);
	return tab;
}

float distance(float* src, const int s1, const int s2, int n) {
	int i;
	float sqSum = 0;

	for (i = 0; i < n; i++) {
		sqSum += pow(src[s1 + i] - src[s2 + i], 2);
	}

	return sqrt(sqSum);
}

void freeArray(float* tab) {
	free(tab);
}

void force(float* dest, float* src, const int s1, const int s2, float m1, float m2, int n)
{
	int i;
	float Ftmp = -G * m1 * m2 / pow(distance(src, s1, s2, n), 3);
	for (i = 0; i < n; i++)
		dest[i] = Ftmp * (src[s1 + i] - src[s2 + i]);
}

void euler(float* r, float* f, float dt, int n)
{
	for (int i = 0; i < n; i++)
		r[i] += f[i] * dt;
}

float energyK(float* r, const int s, float m, int n) {
	int k;
	float sumv = 0;
	for (k = 0; k < n; k++) {
		sumv += r[s + k] * r[s + k];
	}
	return m * sumv / 2;
}

float energyP(float* r, float m1, float m2, const int s1, const int s2, int n) {
	float dist = distance(r, s1, s2, n);
	return -G * m1 * m2 / dist;
}


void obliczf(float* f, float* r, float t, const float MZ, const float MK, const float MR, int n, const int dim) {
	float* tmp1, * tmp2;
	int j, k;

	tmp1 = createArray(dim);
	tmp2 = createArray(dim);

	for (j = 0; j < n; j += 2 * dim) {
		for (k = 0; k < dim; k++)
			f[j + k] = r[dim + k + j];
	}

	force(tmp1, r, 0, 2 * dim, MZ, MK, dim);
	force(tmp2, r, 0, 2 * 2 * dim, MZ, MR, dim);

	for (k = dim, j = 0; j < dim; k++, j++) {
		f[k] = (tmp1[j] + tmp2[j]) / MZ;
	}

	force(tmp1, r, 2 * dim, 0, MK, MZ, dim);
	force(tmp2, r, 2 * dim, 2 * 2 * dim, MK, MR, dim);

	for (k = 2 * dim + dim, j = 0; j < dim; k++, j++) {
		f[k] = (tmp1[j] + tmp2[j]) / MK;
	}

	force(tmp1, r, 2 * 2 * dim, 0, MR, MZ, dim);
	force(tmp2, r, 2 * 2 * dim, 2 * dim, MR, MK, dim);

	for (k = 2 * 2 * dim + dim, j = 0; j < dim; k++, j++) {
		f[k] = (tmp1[j] + tmp2[j]) / MR;
	}
	freeArray(tmp1);
	freeArray(tmp2);
}

void rkf45(float* r, float* f, float* r4, float* r5, float t, float h, int n, const float MZ, const float MK, const float MR, const int dim)
{
	float* k1, *k2, *k3, *k4, *k5, *k6, *rTmp;

	k1 = createArray(n);
	k2 = createArray(n);
	k3 = createArray(n);
	k4 = createArray(n);
	k5 = createArray(n);
	k6 = createArray(n);
	rTmp = createArray(n);

	obliczf(f, r, t, MZ, MK, MR, n, dim);
	for (int i = 0; i < n; i++)
	{
		k1[i] = f[i] * h;            // k1
		rTmp[i] = r[i] + k1[i] / 4;
	}

	obliczf(f, rTmp, t + h / 4, MZ, MK, MR, n, dim);
	for (int i = 0; i < n; i++)
	{
		k2[i] = f[i] * h;            // k2
		rTmp[i] = r[i] + (k1[i] * 3 + k2[i] * 9) / 32;
	}

	obliczf(f, rTmp, t + h * 3 / 8, MZ, MK, MR, n, dim);
	for (int i = 0; i < n; i++)
	{
		k3[i] = f[i] * h;            // k3
		rTmp[i] = r[i] + (k1[i] * 1932 - k2[i] * 7200 + k3[i] * 7296) / 2197;
	}

	obliczf(f, rTmp, t + h * 12 / 13, MZ, MK, MR, n, dim);
	for (int i = 0; i < n; i++)
	{
		k4[i] = f[i] * h;            // k4
		rTmp[i] = r[i] + k1[i] * 439 / 216 - k2[i] * 8 + k3[i] * 3680 / 513 - k4[i] * 845 / 4104;
	}

	obliczf(f, rTmp, t + h, MZ, MK, MR, n, dim);
	for (int i = 0; i < n; i++)
	{
		k5[i] = f[i] * h;            // k5
		r4[i] = r[i] + k1[i] * 25 / 216 + k3[i] * 1408 / 2565 + k4[i] * 2197 / 4104 - k5[i] / 5;      // RK4
		rTmp[i] = r[i] - k1[i] * 8 / 27 + k2[i] * 2 - k3[i] * 3544 / 2565 + k4[i] * 1859 / 4104 - k5[i] * 11 / 40;
	}

	obliczf(f, rTmp, t + h / 2, MZ, MK, MR, n, dim);
	for (int i = 0; i < n; i++)
	{
		k6[i] = f[i] * h;           // k6
		r5[i] = r[i] + k1[i] * 16 / 135 + k3[i] * 6656 / 12825 + k4[i] * 28561 / 56430 - k5[i] * 9 / 50 + k6[i] * 2 / 55;  //RK5
	}

	freeArray(k6);
	freeArray(k5);
	freeArray(k4);
	freeArray(k3);
	freeArray(k2);
	freeArray(k1);
	freeArray(rTmp);
}

void calculate(float* r, const float MZ, const float MK, const float MR, const int n, const int dim) {
	int i = 0, k = 0;
	float t0 = 0;
	float t = t0;
	float dt = 60;
	float distRZ = 0, distRK = 0, distZK = 0;
	float hmin = 0.01, hmax = 1000;
	float *f, *r4, *r5;
	float Ec, EKZ, EKR, EKK, EPZK, EPKR, EPZR;
	float epsmin = 0.1, epsmax = 10, eps;
	int counter;

	FILE* wyniki = fopen("./wyniki.txt", "w");

	fprintf(wyniki, "#	t		Polozenie Ziemi			Polozenie Ksiezyca		Polozenie rakiety		Predkosc ziemi			Predkosc Ksiezyca		Predkosc rakiety		Ec		dt\n");
	if (dim < 2) {
		printf("Wymiar musi być conajmniej równy 2");
		return;
	}

	f = createArray(n);
	r4 = createArray(n);
	r5 = createArray(n);

	i = 0;
	while (1) {
		int j;

		distZK = distance(r, 0, 2 * dim, dim);
		distRK = distance(r, 4, 2 * 2 * dim, dim);
		distRZ = distance(r, 0, 2 * 2 * dim, dim);

		if (distRK <= RK) {
			printf("\nRakieta trafila w ksiezyc");
			fprintf(wyniki, "\nRakieta trafila w ksiezyc");
			break;
		}
		else if (distRK > RK && distRZ > distZK) {
			printf("\nRakieta nie trafila w ksiezyc");
			fprintf(wyniki, "\nRakieta nie trafila w ksiezyc");
			break;
		}

		EKZ = energyK(r, 2, MZ, dim);
		EKK = energyK(r, 2 * dim + dim, MK, dim);
		EKR = energyK(r, 2 * 2 * dim + dim, MR, dim);
		EPZK = energyP(r, MZ, MK, 0, 2 * dim, dim);
		EPZR = energyP(r, MZ, MR, 0, 2 * 2 * dim, dim);
		EPKR = energyP(r, MK, MR, 2 * dim, 2 * 2 * dim, dim);
		Ec = EKZ + EKK + EKR + EPZK + EPZR + EPKR;

		fprintf(wyniki, "%d\t%lf	", ++i, t);
		for (j = 0; j < n; j++) {
			fprintf(wyniki, "%e\t", r[j]);
		}
		fprintf(wyniki, "%e\t", Ec);
		fprintf(wyniki, "%e\t", dt);
		fprintf(wyniki, "\n");

		counter = 0;
		do {
			rkf45(r, f, r4, r5, t, dt, n, MZ, MK, MR, dim);
			eps = 0;
			for (j = 0; j < n; j++) {
				eps += fabs((double)r5[j] - r4[j]);
			}
			eps = sqrt(eps);
			if (eps < epsmin) dt *= 2;
			else if (eps > epsmax) dt /= 2;
			counter++;
		} while ((eps > epsmax || eps < epsmin) && counter < 100);

		if (dt < hmin) dt = hmin;
		else if (dt > hmax) dt = hmax;

		memcpy(r, r4, n * sizeof(float));

		t += dt;
	}
	freeArray(f);
	freeArray(r4);
	freeArray(r5);
	fclose(wyniki);
}

int main()
{
	int dim = 2;		// size
	const float MZ = 5.97219e24;	// kg
	const float MR = 1e4;			// kg
	const float MK = 7.347673e22;	// kg
	float* r;
	int n = N * 2 * dim;
	int alfa;
	float rad;

	while (1) {
		printf("Podaj kat: ");
		scanf("%d", &alfa);
		if (alfa >= 0 && alfa <= 360) break;
		else printf("Podaj odpowiednia wartosc kata!\n");
	}
	rad = (alfa * PI) / 180;

	r = createArray(n);

	//Ziemia
	r[0] = 0;			// x
	r[1] = 0;			// y
	r[2] = 0;			// vx
	r[3] = 0;			// vy
	//Ksiezyc
	r[4] = 405696e3;	// x
	r[5] = 0;			// y
	r[6] = 0;			// vx
	r[7] = 0.968e3;		// vy
	//rakieta
	r[8] = RZ * cos(rad);	// x
	r[9] = RZ * sin(rad);	// y
	r[10] = vII * cos(rad);	// vx
	r[11] = vII * sin(rad);	// vy

	calculate(r, MZ, MK, MR, n, dim);

	freeArray(r);
}
