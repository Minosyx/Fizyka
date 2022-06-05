#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define G 6.67408e-11
#define RZ 6371.008e3 //promien Ziemi
#define vII 11.19e3 // druga predkosc kosmiczna
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

void force(float* F, float* r1, float* r2, float m1, float m2, int n)
{
	int i;
	float Ftmp = -G * m1 * m2 / pow(distance(r1, r2, n), 3);
	for (i = 0; i < n; i++)
		F[i] = Ftmp * (r1[i] - r2[i]);
}

void calculate(float* rZ, float* rR, float* vZ, float* vR, float* FZ, float* FR, const float MZ, const float MR, int n, const float repsR, const float vepsR) {
	int i = 0;
	float t0 = 0;
	float t = t0;
	float dt = 0;
	float tr, tv;
	float distR = 0;
	float absF, dabsF, absF2;
	float distSum = 0;

	FILE* wyniki = fopen("./wyniki.txt", "w");

	float* Fold;
	Fold = createArray(n);

	if (n == 2)
		fprintf(wyniki, "#	t		xZ		yZ		xR		yR		vxZ		vyZ		vxR		vyR		tr		tv		dt\n");
	else if (n == 3)
		fprintf(wyniki, "#	t		xZ		yZ		zZ		xR		yR		zR		vxZ		vyZ		vzZ		vxR		vyR		vzR		tr		tv		dt\n");
	else if (n < 2) {
		printf("Wymiar musi być conajmniej równy 2");
		return;
	}

	while (distR <= 405696e3) {
		int j;

		force(FZ, rZ, rR, MZ, MR, n);
		force(FR, rR, rZ, MR, MZ, n);

		distR = distance(rR, rZ, n);

		// szacowanie kroku czasowego

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

		float tempt = 0;
		if (t != 0) {
			int sumS = 0;
			for (j = 0; j < n; j++) {
				sumS += pow((FR[j] - Fold[j]), 2);
			}
			sumS = sqrt(sumS);

			tempt = sqrt((2 * MR * dt * vepsR) / sumS);
		}

		tr = sqrt((2 * MR * repsR) / absF);
		if (t != 0)
			tv = tempt;//sqrt((2 * MR * vepsR) / dabsF);
		else tv = 0;

		printf("%e : %e\n", tv, tempt);

		if (tr > tv && tv != 0) dt = tv;
		else dt = tr;

		// koniec szacowania

		fprintf(wyniki, "%d	%lf	", ++i, t);
		for (j = 0; j < n; j++) fprintf(wyniki, "%e	", rZ[j]);
		for (j = 0; j < n; j++) fprintf(wyniki, "%e	", rR[j]);
		for (j = 0; j < n; j++) fprintf(wyniki, "%e	", vZ[j]);
		for (j = 0; j < n; j++) fprintf(wyniki, "%e	", vR[j]);
		fprintf(wyniki, "%e	%e	%e", tr, tv, dt);
		fprintf(wyniki, "\n");

		for (j = 0; j < n; j++) {
			rZ[j] += vZ[j] * dt;
			rR[j] += vR[j] * dt;
			vZ[j] += (FZ[j] / MZ) * dt;
			vR[j] += (FR[j] / MR) * dt;
			Fold[j] = FR[j];
		}
		t += dt;
	}
	freeArray(Fold);
	fclose(wyniki);
}

int main()
{
	int n = 2;		// size
	const float MZ = 5.97219e24;
	const float MR = 1e4;
	float* rZ;
	float* FZ;
	float* vZ;
	float* rR;
	float* FR;
	float* vR;
	const float repsR = 10; // epsilon rR
	const float vepsR = 0.01; // epsilon vR
	int alfa;
	float rad;

	while (1) {
		printf("Podaj kat: ");
		scanf("%d", &alfa);
		if (alfa >= 0 && alfa <= 360) break;
		else printf("Podaj odpowiednia wartosc kata!\n");
	}
	rad = (alfa * PI) / 180;

	rZ = createArray(n);
	FZ = createArray(n);
	vZ = createArray(n);
	rR = createArray(n);
	FR = createArray(n);
	vR = createArray(n);

	rZ[0] = 0;
	rZ[1] = 0;
	vZ[0] = 0;
	vZ[1] = 0;

	rR[0] =	RZ * cos(rad);
	rR[1] = RZ * sin(rad);
	vR[0] = vII * cos(rad);
	vR[1] = vII * sin(rad);

	calculate(rZ, rR, vZ, vR, FZ, FR, MZ, MR, n, repsR, vepsR);

	freeArray(rZ);
	freeArray(rR);
	freeArray(FZ);
	freeArray(FR);
	freeArray(vZ);
	freeArray(vR);
}
