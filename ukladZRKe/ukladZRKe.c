#include <stdio.h>
#include <stdlib.h>
#include <math.h>
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

void calculate(float* r, float* f, const float MZ, const float MK, const float MR, const int n, const int dim, const float repsR, const float vepsR) {
	int i = 0, k = 0;
	float t0 = 0;
	float t = t0;
	float dt = 0;
	float tr, tv;
	float distRZ = 0, distRK = 0, distZK = 0;
	float absF, dabsF;
	float dSum = 0;
	float* tmp1, * tmp2;
	float dtmin = 0.1, dtmax = 60;
	float* Aold;
	float Ec, EKZ, EKR, EKK, EPZK, EPKR, EPZR;

	FILE* wyniki = fopen("./wyniki.txt", "w");

	fprintf(wyniki, "#	t		Polozenie Ziemi			Polozenie Ksiezyca		Polozenie rakiety		Predkosc ziemi			Predkosc Ksiezyca		Predkosc rakiety		Ec		tr		tv		dt\n");
	if (dim < 2) {
		printf("Wymiar musi być conajmniej równy 2");
		return;
	}

	tmp1 = createArray(dim);
	tmp2 = createArray(dim);
	Aold = createArray(dim);

	i = 0;
	while (1) {
		int j;

		distZK = distance(r, 0, 2*dim, dim);
		distRK = distance(r, 4, 2*2*dim, dim);
		distRZ = distance(r, 0, 2*2*dim, dim);

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

		for (j = 0; j < dim; j++) {
			Aold[j] = f[2 * 2 * dim + dim + j];
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

		absF = 0;
		for (j = 0; j < dim; j++) {
			absF += pow((double)f[2 * 2 * dim + dim + j] * MR, 2);
		}
		absF = sqrt(absF);
		
		dabsF = 0;
		if (t != 0) {
			for (j = 0; j < dim; j++) {
				dabsF += pow(((double)f[2 * 2 * dim + dim + j] - Aold[j]) * MR, 2);
			}
			dabsF = sqrt(dabsF);
		}
		
		tr = sqrt((2 * MR * repsR) / absF);
		if (t != 0)
			tv = sqrt((2 * MR * dt * vepsR) / dabsF);
		else tv = 0;

		if (tr > tv && tv != 0) dt = tv;
		else dt = tr;

		if (dt > dtmax) dt = dtmax;
		else if (dt < dtmin) dt = dtmin;

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
		fprintf(wyniki, "%e\t%e\t%e", tr, tv, dt);
		fprintf(wyniki, "\n");

		euler(r, f, dt, n);

		for (k = 0; k < N; k++) {
			for (j = 0; j < dim; j++)
				f[2 * k * dim + j] = r[2 * k * dim + dim + j];
		}
		t += dt;
	}

	freeArray(tmp1);
	freeArray(tmp2);
	fclose(wyniki);
}

int main()
{
	int dim = 2;		// size
	const float MZ = 5.97219e24;	// kg
	const float MR = 1e4;			// kg
	const float MK = 7.347673e22;	// kg
	float* r, *f;
	int n = N * 2 * dim;
	const float repsR = 10;		// m
	const float vepsR = 0.01;	// m/s
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
	f = createArray(n);

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

	f[0] = r[2];
	f[1] = r[3];
	f[4] = r[6];
	f[5] = r[7];
	f[8] = r[10];
	f[9] = r[11];

	calculate(r, f, MZ, MK, MR, n, dim, repsR, vepsR);

	freeArray(r);
	freeArray(f);
}
