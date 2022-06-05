#include <stdio.h>
#include <stdlib.h>
#define g 9.81


double* createArray(int size) {
	double* tab;
	tab = (double*)malloc(sizeof(double)*size);
	return tab;
}

void freeArray(double* tab) {
	free(tab);
}

void calculate(double* r, double* v, double* F, int size, double m, double dt, double tN) {
	int i = 0;
	double t0 = 0;
	double t = t0;
	double Ec, Ep;

	FILE* wyniki = fopen("./wyniki.txt", "w");

	if (size == 2)	fprintf(wyniki, "#	t		x		y		vx		vy		Ec\n");
	else if (size == 3) fprintf(wyniki, "#	t		x		y		z		vx		vy		vz		Ec\n");
	else if (size == 1) fprintf(wyniki, "#	t		x		vx		Ec\n");

	while (t < tN + dt || (r[size - 1] > 0)) {
		int j;
		double sumav2 = 0;

		i++;

		for (j = 0; j < size; j++) {	//zliczanie v^2
			sumav2 += v[j] * v[j];
		}
		Ep = m * g * r[size - 1];
		Ec = Ep + (m / 2 * sumav2);

		fprintf(wyniki, "%d	%lf	", i, t);
		for (j = 0; j < size; j++) fprintf(wyniki, "%e	", r[j]);
		for (j = 0; j < size; j++) fprintf(wyniki, "%e	", v[j]);
		fprintf(wyniki, "%e\n", Ec);

		for (j = 0; j < size; j++) {
			r[j] += v[j] * dt;
			v[j] += (F[j] / m) * dt;
		}
		t += dt;
	}

	fclose(wyniki);
}


int main()
{

	/*Pierwsze z trzech zadan nalezy rozpatrywac w 3D, gdyz energia potencjalna okreslona jest w programie na stale jako m * g * h, gdzie h to z w 3D, y w 2D, x w 1D.
	Rzut ukosny mozna rozpatrywac w 2D i 3D, a spadek swobodny w 1D, 2D badz 3D.*/
	int size, i;
	double* r;
	double* v;
	double* F;
	double m, tN, dt;

	printf("Podaj ilosc rozpatrywanych wymiarow: ");	scanf("%d", &size);
	r = createArray(size);
	v = createArray(size);
	F = createArray(size);
	printf("Podaj mase ciala (kg): ");	scanf("%lf", &m);
	printf("Podaj ograniczenie czasowe ruchu (jezeli wysokosc = 0 jest warunkiem konca obliczen wpisz 0) (s): ");	scanf("%lf", &tN);
	printf("Podaj wartosc kroku czasowego dt (s): ");	scanf("%lf", &dt);

	printf("\nPodaj wartosc skladowych polozenia (m)\n");

	for (i = 0; i < size; i++) {
		printf("%d. skladowa: ", i + 1);	scanf("%lf", &r[i]);
	}
	
	printf("\nPodaj wartosc skladowych predkosci (m/s)\n");

	for (i = 0; i < size; i++) {
		printf("%d. skladowa: ", i + 1);	scanf("%lf", &v[i]);
	}

	printf("\nPodaj wartosc skladowych sily (za g przyjmij 9.81)(N)\n");

	for (i = 0; i < size; i++) {
		printf("%d. skladowa: ", i + 1);	scanf("%lf", &F[i]);
	}

	calculate(r, v, F, size, m, dt, tN);

	freeArray(r);
	freeArray(v);
	freeArray(F);
}
