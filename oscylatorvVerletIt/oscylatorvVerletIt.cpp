#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float* createArray(int n) {
	float* tab;
	tab = (float*)malloc(sizeof(float) * n);
	return tab;
}

void freeArray(float* tab) {
	free(tab);
}

void force(float* F, float* r, float* v, float k, float b, int n)
{
	int i;
	for (i = 0; i < n; i++)
		F[i] = -k * r[i] - b * v[i];
}

float energy(float* r, float* v, float k, int n, float m) {
	int i;
	float sv = 0, sr = 0;
	for (i = 0; i < n; i++) {
		sv += v[i] * v[i];
		sr += r[i] * r[i];
	}
	return m * sv / 2 + k * sr / 2;
}

void calculate(float* r, float* v, float* F, int n, float m, float dt, float tN, float k, float b) {
	int i = 0;
	float t0 = 0;
	float t = t0;
	float Ec, Wop = 0;
	float thr = 1e-7;
	float vDiff = 2*thr;
	int counter = 0;
	int maxCount = 100;

	float* Fold;
	float* rOld;
	float* vTemp;
	float* vOld;

	Fold = createArray(n);
	rOld = createArray(n);
	vTemp = createArray(n);
	vOld = createArray(n);


	FILE* wyniki = fopen("./wyniki.txt", "w");

	fprintf(wyniki, "#	t		x		y		vx		vy		Ec		Ec-Wop\n");

	force(F, r, v, k, b, n);

	while (t < tN + dt) {
		int j;
		i++;

		Ec = energy(r, v, k, n, m);

		fprintf(wyniki, "%d	%f	", i, t);
		for (j = 0; j < n; j++) fprintf(wyniki, "%e	", r[j]);
		for (j = 0; j < n; j++) fprintf(wyniki, "%e	", v[j]);
		fprintf(wyniki, "%e	%e\n", Ec, Ec - Wop);

		for (j = 0; j < n; j++) {
			rOld[j] = r[j];
			r[j] += v[j] * dt + F[j] / (2 * m) * dt * dt;
			Wop += -b * v[j] * (r[j] - rOld[j]);
			Fold[j] = F[j];
			vTemp[j] = v[j] + F[j] / m * dt;
		}

		do{
			float sumDiff = 0, sumLen = 0;
			float vLength;
			force(F, r, vTemp, k, b, n);
			counter++;

			for (j = 0; j < n; j++) {
				vOld[j] = vTemp[j];
				vTemp[j] = v[j] + (F[j] + Fold[j]) / (2 * m) * dt;
				sumDiff += pow(vTemp[j] - vOld[j], 2);
				sumLen += pow(vTemp[j], 2);
			}

			vDiff = sqrt(sumDiff);
			vLength = sqrt(sumLen);

			if (vLength > thr)
				vDiff /= vLength;
		} while (vDiff > thr && counter < maxCount);

		for (j = 0; j < n; j++) {
			v[j] = vTemp[j];
		}

		t += dt;

	}
	freeArray(Fold);
	freeArray(rOld);
	freeArray(vTemp);
	freeArray(vOld);
	fclose(wyniki);
}


int main()
{
	int n = 2;		// size
	float k = 2;		// kg/s^2
	float m = 0.5;		// kg
	float b = 0.1;		// kg/s
	float* r;
	float* F;
	float* v;
	float dt = 0.01;
	float tN = 20;

	r = createArray(n);
	F = createArray(n);
	v = createArray(n);

	r[0] = 0.5;
	r[1] = 0.75;
	v[0] = 0;
	v[1] = 0;

	calculate(r, v, F, n, m, dt, tN, k, b);

	freeArray(r);
	freeArray(F);
	freeArray(v);
}
