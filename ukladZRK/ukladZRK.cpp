#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define G 6.67408e-11
#define RZ 6371.008e3
#define RK 1737.064e3
#define vII 11.19e3
#define PI 3.141592654

float* createArray(int size) {
	float* tab;
	tab = (float*)malloc(sizeof(float) * size);
	return tab;
}

float distance(float* r1, float* r2, int n) {
	int i;
	float sqSum = 0;

	for (i = 0; i < n; i++) {
		sqSum += pow(r1[i] - r2[i], 2);
	}

	return sqrt(sqSum);
}

void freeArray(float* tab) {
	free(tab);
}

void force(float *dest, float* r1, float* r2, float m1, float m2, int n)
{
	int i;
	float Ftmp = -G * m1 * m2 / pow(distance(r1, r2, n), 3);
	for (i = 0; i < n; i++)
		dest[i] = Ftmp * (r1[i] - r2[i]);
}

void euler(float* r, float* f, float dt, int n)
{
	for (int i = 0; i < n; i++)
		r[i] += f[i] * dt;
}

void calculate(float* rZ, float* rK, float* rR, float* vZ, float* vK, float* vR, float* FZ, float* FK, float* FR, const float MZ, const float MK, const float MR, int n, const float repsR, const float vepsR) {
	int i = 0, k = 0;
	int dim = 3 * 2 * n;
	float t0 = 0;
	float t = t0;
	float dt = 0;
	float tr, tv;
	float distRZ = 0, distRK = 0, distZK = 0;
	float absF, dabsF;
	float distSum = 0;
	float* r, * f, * tmp1, * tmp2;
	float dtmin = 0.1, dtmax = 10;

	FILE* wyniki = fopen("./wyniki.txt", "w");

	fprintf(wyniki, "#	t		Polozenie Ziemi			Polozenie Ksiezyca		Polozenie rakiety		Predkosc ziemi			Predkosc Ksiezyca		Predkosc rakiety		tr		tv		dt\n");
	if (n < 2) {
		printf("Wymiar musi być conajmniej równy 2");
		return;
	}

	r = createArray(dim);
	f = createArray(dim);
	tmp1 = createArray(n);
	tmp2 = createArray(n);

	force(tmp1, rZ, rK, MZ, MK, n);
	force(tmp2, rZ, rR, MZ, MR, n);
	for (k = 0; k < n; k++) {
		FZ[k] = tmp1[k] + tmp2[k];
	}
	force(tmp1, rK, rZ, MK, MZ, n);
	force(tmp2, rK, rR, MK, MR, n);
	for (k = 0; k < n; k++) {
		FK[k] = tmp1[k] + tmp2[k];
	}
	force(tmp1, rR, rZ, MR, MZ, n);
	force(tmp2, rR, rK, MR, MK, n);
	for (k = 0; k < n; k++) {
		FR[k] = tmp1[k] + tmp2[k];
	}

	for (k = 0; k < n; k++) {
		r[k] = rZ[k];
		f[k] = vZ[k];
	}
	for (k = n; k < 2*n; k++) {
		r[k] = vZ[k - n];
		f[k] = FZ[k] / MZ;
	}
	for (k = 2 * n; k < 3 * n; k++) {
		r[k] = rK[k - 2 * n];
		f[k] = vK[k - 2 * n];
	}
	for (k = 3 * n; k < 4 * n; k++) {
		r[k] = vK[k - 3 * n];
		f[k] = FK[k] / MK;
	}
	for (k = 4 * n; k < 5 * n; k++) {
		r[k] = rR[k - 4 * n];
		f[k] = vR[k - 4 * n];
	}
	for (k = 5 * n; k < 6 * n; k++) {
		r[k] = vR[k - 5 * n];
		f[k] = FR[k] / MR;
	}


	while (1) {
		int j;

		distZK = distance(rZ, rK, n);
		distRK = distance(rR, rK, n);
		distRZ = distance(rR, rZ, n);

		if (distRK <= RK) break;
		else if (distRK > RK && distRZ > distZK) break;

		absF = 0;
		for (j = 0; j < n; j++) {
			absF += pow(FR[j], 2);
		}
		absF = sqrt(absF);

		for (j = 0; j < n; j++) {
			distSum += pow(rR[j], 2);
		}
		distSum = sqrt(distSum);
		dabsF = 2 * G * MZ * MR / pow(distSum, 3) * sqrt(2 * G * MZ / distSum);

		tr = sqrt((2 * MR * repsR) / absF);
		if (t != 0)
			tv = sqrt((2 * MR * vepsR) / dabsF);
		else tv = 0;

		if (tr > tv && tv != 0) dt = tv;
		else if (tr < tv || (tr > tv && tv == 0)) dt = tr;

		if (dt > dtmax) dt = dtmax;
		else if (dt < dtmin) dt = dtmin;

		fprintf(wyniki, "%d	%lf	", ++i, t);
		for (j = 0; j < n; j++) fprintf(wyniki, "%e	", rZ[j]);
		for (j = 0; j < n; j++) fprintf(wyniki, "%e	", rK[j]);
		for (j = 0; j < n; j++) fprintf(wyniki, "%e	", rR[j]);
		for (j = 0; j < n; j++) fprintf(wyniki, "%e	", vZ[j]);
		for (j = 0; j < n; j++) fprintf(wyniki, "%e	", vK[j]);
		for (j = 0; j < n; j++) fprintf(wyniki, "%e	", vR[j]);
		fprintf(wyniki, "%e	%e	%e", tr, tv, dt);
		fprintf(wyniki, "\n");

		euler(r, f, dt, dim);

		for (k = 0; k < n; k++) {
			rZ[k] = r[k];
			vZ[k] = f[k];
		}
		for (k = 2 * n; k < 3 * n; k++) {
			rK[k - 2 * n] = r[k];
			vK[k - 2 * n] = f[k];
		}
		for (k = 4 * n; k < 5 * n; k++) {
			rR[k - 4 * n] = r[k];
			vR[k - 4 * n] = f[k];
		}

		force(tmp1, rZ, rK, MZ, MK, n);
		force(tmp2, rZ, rR, MZ, MR, n);
		for (k = 0; k < n; k++) {
			FZ[k] = tmp1[k] + tmp2[k];
		}
		force(tmp1, rK, rZ, MK, MZ, n);
		force(tmp2, rK, rR, MK, MR, n);
		for (k = 0; k < n; k++) {
			FK[k] = tmp1[k] + tmp2[k];
		}
		force(tmp1, rR, rZ, MR, MZ, n);
		force(tmp2, rR, rK, MR, MK, n);
		for (k = 0; k < n; k++) {
			FR[k] = tmp1[k] + tmp2[k];
		}

		t += dt;
	}

	freeArray(tmp1);
	freeArray(tmp2);
	fclose(wyniki);
}

int main()
{
	int n = 2;		// size
	const float MZ = 5.97219e24;
	const float MR = 1e4;
	const float MK = 7.347673e22;
	float* rZ, *rR, *rK;
	float* FZ, *FR, *FK;
	float* vZ, *vR, *vK;
	const float repsR = 10;
	const float vepsR = 0.1;
	int alfa;
	float rad;

	while (1) {
		printf("Podaj kat: ");
		scanf("%d", &alfa);
		if (alfa >= 0 && alfa <= 180) break;
		else printf("Podaj odpowiednia wartosc kata!\n");
	}
	rad = (alfa * PI) / 180;

	rZ = createArray(n);
	FZ = createArray(n);
	vZ = createArray(n);
	rR = createArray(n);
	FR = createArray(n);
	vR = createArray(n);
	rK = createArray(n);
	vK = createArray(n);
	FK = createArray(n);

	rZ[0] = 0;
	rZ[1] = 0;
	vZ[0] = 0;
	vZ[1] = 0;

	rK[0] = 405696e3;
	rK[1] = 0;
	vK[0] = 0;
	vK[1] = 0.968e3;

	rR[0] = RZ * cos(rad);
	rR[1] = RZ * sin(rad);
	vR[0] = vII * cos(rad);
	vR[1] = vII * sin(rad);

	calculate(rZ, rK, rR, vZ, vK, vR, FZ, FK, FR, MZ, MK, MR, n, repsR, vepsR);

	freeArray(rZ);
	freeArray(rR);
	freeArray(FZ);
	freeArray(FR);
	freeArray(vZ);
	freeArray(vR);
	freeArray(FK);
	freeArray(vK);
	freeArray(rK);
}
