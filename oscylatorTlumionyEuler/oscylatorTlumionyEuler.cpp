#include <stdio.h>
#include <stdlib.h>

float* createArray(int size) {
	float* tab;
	tab = (float*)malloc(sizeof(float) * size);
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
	float* rOld;

	rOld = createArray(n);

	FILE* wyniki = fopen("./wyniki.txt", "w");

	fprintf(wyniki, "#	t		x		y		vx		vy		Ec		Ec-Wop\n");

	while (t < tN + dt) {
		int j;

		i++;

		force(F, r, v, k, b, n);
		Ec = energy(r, v, k, n, m);

		fprintf(wyniki, "%d	%lf	", i, t);
		for (j = 0; j < n; j++) fprintf(wyniki, "%e	", r[j]);
		for (j = 0; j < n; j++) fprintf(wyniki, "%e	", v[j]);
		fprintf(wyniki, "%e	%e\n", Ec, Ec-Wop);

		for (j = 0; j < n; j++) {
			rOld[j] = r[j];
			r[j] += v[j] * dt;
			Wop += -b * v[j] * (r[j] - rOld[j]);
			v[j] += (F[j] / m) * dt;
		}

		t += dt;
	}

	freeArray(rOld);
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
