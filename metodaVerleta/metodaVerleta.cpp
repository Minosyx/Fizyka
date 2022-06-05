#include <stdio.h>
#include <stdlib.h>
#define g 9.81


float* createArray(int size) {
	float* tab;
	tab = (float*)malloc(sizeof(float) * size);
	return tab;
}

void freeArray(float* tab) {
	free(tab);
}

void calculate(float* r, float* v, float* F, int size, float m, float dt, float tN) {
	int i = 0;
	int k;
	float t0 = 0;
	float t = t0;
	float Ec, Ep;
	float* rPlus;
	float* rMinus;

	rPlus = createArray(size);
	rMinus = createArray(size);

	FILE* wyniki = fopen("./wyniki.txt", "w");

	if (size == 2)	fprintf(wyniki, "#	t		x		y		vx		vy		Ec\n");
	else if (size == 3) fprintf(wyniki, "#	t		x		y		z		vx		vy		vz		Ec\n");
	else if (size == 1) fprintf(wyniki, "#	t		x		vx		Ec\n");

	for (k = 0; k < size; k++) {
		rPlus[k] = r[k] + v[k] * dt + F[k] / (2 * m) * dt * dt;
	}

	while (t < tN + dt || ((r[size - 1] > 0) && F[size - 1]) != 0) {
		int j;
		float sumav2 = 0;

		i++;

		for (j = 0; j < size; j++) {	//zliczanie v^2
			sumav2 += v[j] * v[j];
		}
		if (F[size - 1] != 0) Ep = m * g * r[size - 1];
		else Ep = 0;
		Ec = Ep + (m / 2 * sumav2);

		fprintf(wyniki, "%d	%f	", i, t);
		for (j = 0; j < size; j++) fprintf(wyniki, "%e	", r[j]);
		for (j = 0; j < size; j++) fprintf(wyniki, "%e	", v[j]);
		fprintf(wyniki, "%e\n", Ec);

		for (j = 0; j < size; j++) {
			rMinus[j] = r[j];
			r[j] = rPlus[j];

			rPlus[j] = -rMinus[j] + 2 * r[j] + F[j] / m * dt * dt;
			v[j] = (rPlus[j] - rMinus[j]) / (2 * dt);
		}
		t += dt;
	}

	fclose(wyniki);
    freeArray(rPlus);
    freeArray(rMinus);
}


int main()
{
	int size, i, choice = 0;
	float* r;
	float* v;
	float* F;
	float m, tN, dt;

	printf("Podaj ilosc rozpatrywanych wymiarow: ");	scanf("%d", &size);
	r = createArray(size);
	v = createArray(size);
	F = createArray(size);
	printf("Podaj mase ciala (kg): ");	scanf("%f", &m);
	printf("Podaj ograniczenie czasowe ruchu (jezeli wysokosc = 0 jest warunkiem konca obliczen wpisz 0) (s): ");	scanf("%f", &tN);
	printf("Podaj wartosc kroku czasowego dt (s): ");	scanf("%f", &dt);

	printf("\nPodaj wartosc skladowych polozenia (m)\n");

	for (i = 0; i < size; i++) {
		printf("%d. skladowa: ", i + 1);	scanf("%f", &r[i]);
	}

	printf("\nPodaj wartosc skladowych predkosci (m/s)\n");

	for (i = 0; i < size; i++) {
		printf("%d. skladowa: ", i + 1);	scanf("%f", &v[i]);
	}

	printf("\nWybierz wartosc skladowych sily\n");

	for (i = 0; i < size; i++) {
		while (choice != 1 || choice != 2) {
			printf("%d. skladowa: (1. 0 2. -mg)\n", i + 1);
			printf("Podaj swoj wybor: "); scanf("%d", &choice);
			switch (choice) {
			case 1:
				F[i] = 0;
				break;
			case 2:
				F[i] = -m * g;
				break;
			default:
				printf("Wybrano bledna opcje!!!\n");
			}
			if (choice == 1 || choice == 2) break;
		}
	}

	calculate(r, v, F, size, m, dt, tN);

	freeArray(r);
	freeArray(v);
	freeArray(F);
}
