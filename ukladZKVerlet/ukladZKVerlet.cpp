#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define G 6.67408e-11

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

void force(float* F, float* r1, float* r2, float m1, float m2, int n)
{
	int i;
	float Ftmp = -G * m1 * m2 / pow(distance(r1, r2, n), 3);
	for (i = 0; i < n; i++)
		F[i] = Ftmp * (r1[i] - r2[i]);
}

void calculate(float* rZ, float* rK, float* vZ, float* vK, float* FZ, float* FK, const float MZ, const float MK, int n, float dt, float tN) {
	int i = 0;
	int l;
	float t0 = 0;
	float t = t0;
	float* rZPlus;
	float* rZMinus;
	float* rKPlus;
	float* rKMinus;

	rZPlus = createArray(n);
	rZMinus = createArray(n);
	rKPlus = createArray(n);
	rKMinus = createArray(n);

	FILE* wyniki = fopen("./wyniki.txt", "w");

	if (n == 2)
		fprintf(wyniki, "#	t		xZ		yZ		xK		yK		vxZ		vyZ		vxK		vyK\n");
	else if (n == 3)
		fprintf(wyniki, "#	t		xZ		yZ		zZ		xK		yK		zK		vxZ		vyZ		vzZ		vxK		vyK		vzK\n");
	else if (n < 2) {
		printf("Wymiar musi być conajmniej równy 2");
		return;
	}

	force(FZ, rZ, rK, MZ, MK, n);
	force(FK, rK, rZ, MK, MZ, n);

	for (l = 0; l < n; l++) {
		rZPlus[l] = rZ[l] + vZ[l] * dt + FZ[l] / (2 * MZ) * dt * dt;
		rKPlus[l] = rK[l] + vK[l] * dt + FK[l] / (2 * MK) * dt * dt;
	}

	while (t < tN + dt) {
		int j;

		fprintf(wyniki, "%d	%lf	", ++i, t);
		for (j = 0; j < n; j++) fprintf(wyniki, "%e	", rZ[j]);
		for (j = 0; j < n; j++) fprintf(wyniki, "%e	", rK[j]);
		for (j = 0; j < n; j++) fprintf(wyniki, "%e	", vZ[j]);
		for (j = 0; j < n; j++) fprintf(wyniki, "%e	", vK[j]);
		fprintf(wyniki, "\n");

		for (j = 0; j < n; j++) {
			rZMinus[j] = rZ[j];
			rKMinus[j] = rK[j];
			rZ[j] = rZPlus[j];
			rK[j] = rKPlus[j];

		}
		force(FZ, rZ, rK, MZ, MK, n);
		force(FK, rK, rZ, MK, MZ, n);
		for (j = 0; j < n; j++) {
			rZPlus[j] = -rZMinus[j] + 2 * rZ[j] + FZ[j] / MZ * dt * dt;
			rKPlus[j] = -rKMinus[j] + 2 * rK[j] + FK[j] / MK * dt * dt;
			vZ[j] = (rZPlus[j] - rZMinus[j]) / (2 * dt);
			vK[j] = (rKPlus[j] - rKMinus[j]) / (2 * dt);
		}
		t += dt;
	}

	freeArray(rZPlus);
	freeArray(rKPlus);
	freeArray(rZMinus);
	freeArray(rKMinus);
	fclose(wyniki);
}

int main()
{
	int n = 2;		// size
	const float MZ = 5.97219e24;
	const float MK = 7.347673e22;
	float* rZ;
	float* FZ;
	float* vZ;
	float* rK;
	float* FK;
	float* vK;
	float dt = 0.1 * 3600;
	float tN = 30 * 24 * 3600;

	rZ = createArray(n);
	FZ = createArray(n);
	vZ = createArray(n);
	rK = createArray(n);
	FK = createArray(n);
	vK = createArray(n);

	rZ[0] = 0;
	rZ[1] = 0;
	vZ[0] = 0;
	vZ[1] = 0;

	rK[0] = 405696e3;
	rK[1] = 0;
	vK[0] = 0;
	vK[1] = 0.968e+3;

	calculate(rZ, rK, vZ, vK, FZ, FK, MZ, MK, n, dt, tN);

	freeArray(rZ);
	freeArray(rK);
	freeArray(FZ);
	freeArray(FK);
	freeArray(vZ);
	freeArray(vK);
}
