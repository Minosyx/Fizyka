#include <stdio.h>
#include <stdlib.h>

float* createArray(int n) {
	float* tab;
	tab = (float*)malloc(sizeof(float) * n);
	return tab;
}

void freeArray(float* tab) {
	free(tab);
}

void force(float* F, float* r, float k, int n)
{
	int i;
	for (i = 0; i < n; i++)
		F[i] = -k * r[i];
}

float energy(float* r, float* v, float k, int n, float m) {
	int i;
	float Ec, sv = 0, sr = 0;
	for (i = 0; i < n; i++) {
		sv += v[i] * v[i];
		sr += r[i] * r[i];
	}
	return Ec = m * sv / 2 + k * sr / 2;
}

void calculate(float* r, float* v, float* F, int n, float m, float dt, float tN, float k) {
	int i = 0;
	int l;
	float t0 = 0;
	float t = t0;
	float Ec;
	float* rPlus;
	float* rMinus;

	rPlus = createArray(n);
	rMinus = createArray(n);

	FILE* wyniki = fopen("./wyniki.txt", "w");

	fprintf(wyniki, "#	t		x		y		vx		vy		Ec\n");

	force(F, r, k, n);

	for (l = 0; l < n; l++) {
		rPlus[l] = r[l] + v[l] * dt + F[l] / (2 * m) * dt * dt;
	}

	while (t < tN + dt) {
		int j;
		i++;

		Ec = energy(r, v, k, n, m);
		
		fprintf(wyniki, "%d	%f	", i, t);
		for (j = 0; j < n; j++) fprintf(wyniki, "%e	", r[j]);
		for (j = 0; j < n; j++) fprintf(wyniki, "%e	", v[j]);
		fprintf(wyniki, "%e\n", Ec);

		for (j = 0; j < n; j++) {
			rMinus[j] = r[j];
			r[j] = rPlus[j];
		}
        force(F, r, k, n);
        for (j = 0; j < n; j++) {
            rPlus[j] = -rMinus[j] + 2 * r[j] + F[j] / m * dt * dt;
			v[j] = (rPlus[j] - rMinus[j]) / (2 * dt);
        }
		t += dt;
	}

	freeArray(rPlus);
	freeArray(rMinus);
	fclose(wyniki);
}

int main()
{
	int n = 2;		// n
	float k = 2;		// kg/s^2
	float m = 0.5;		// kg
	float* r;
	float* F;
	float* v;
	float dt = 0.01;
	float tN = 10;

	r = createArray(n);
	F = createArray(n);
	v = createArray(n);

	r[0] = 0.5;
	r[1] = 0.75;
	v[0] = 0;
	v[1] = 0;

	calculate(r, v, F, n, m, dt, tN, k);

	freeArray(r);
	freeArray(F);
	freeArray(v);
}
