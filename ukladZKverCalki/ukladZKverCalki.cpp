#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define G 6.67408e-11

double* createArray(int size) {
	double* tab;
	tab = (double*)malloc(sizeof(double) * size);
	return tab;
}

double distance(double* r1, double* r2, int n) {
	int i;
	double sqSum = 0;

	for (i = 0; i < n; i++) {
		sqSum += pow(r1[i] - r2[i], 2);
	}

	return sqrt(sqSum);
}

void freeArray(double* tab) {
	free(tab);
}

double Ek(double m, double* v, int n) {
	int i;
	double vsum = 0;

	for (i = 0; i < n; i++)
		vsum += v[i] * v[i];

	return m * vsum / 2;
}

double Ep(double m1, double m2, double* r1, double* r2, int n) {
	double dist = distance(r1, r2, n);

	return -G * m1 * m2 / dist;
}

double arealVel(double* r, double* v, int n) {
	double outc = 0;
	if (n == 2) {
		outc = sqrt(pow(r[0] * v[1] - r[1] * v[0], 2));
	}
	else if (n == 3) {
		outc = sqrt(pow(r[0] * v[1] - r[1] * v[0], 2) + pow(r[1] * v[2] - r[2] * v[1], 2) + pow(r[2] * v[0] - r[0] * v[2], 2));
	}
	return (0.5) * outc;
}

void momentum(double m, double* v, double* p, int n) {
	int i;
	for (i = 0; i < n; i++) {
		p[i] = m * v[i];
	}
}

void force(double* F, double* r1, double* r2, double m1, double m2, int n)
{
	int i;
	double Ftmp = -G * m1 * m2 / pow(distance(r1, r2, n), 3);
	for (i = 0; i < n; i++)
		F[i] = Ftmp * (r1[i] - r2[i]);
}

void calculate(double* rZ, double* rK, double* vZ, double* vK, double* FZ, double* FK, const double MZ, const double MK, int n, double dt, double tN) {
	int i = 0;
	int l;
	double t0 = 0;
	double t = t0;
	double* rZPlus;
	double* rZMinus;
	double* rKPlus;
	double* rKMinus;
	double* elipEl;
	double* Pc, * PZ, * PK;
	double elipSum = 0, elipRMS = 0;
	double Ec, aV;

	double time = 0;
	double a, b, c, e;

	double xmax = 0;
	double yXmax = 0;

	double xmin = 0;
	double yXmin = 0;

	double ymax = 0;
	double xYmax = 0;

	double ymin = 0;
	double xYmin = 0;

	int check = 0;

	double* x;
	double* y;

	rZPlus = createArray(n);
	rZMinus = createArray(n);
	rKPlus = createArray(n);
	rKMinus = createArray(n);
	Pc = createArray(n);
	PZ = createArray(n);
	PK = createArray(n);

	x = createArray(int((tN / dt)) + 1);
	y = createArray(int((tN / dt)) + 1);
	elipEl = createArray(int((tN / dt)) + 1);

	FILE* wyniki = fopen("./wyniki.txt", "w");

	FILE* kepler = fopen("./Kepler.txt", "w");

	FILE* calki = fopen("./calki.txt", "w");

	if (n == 2)
		fprintf(wyniki, "#	t		xZ		yZ		xK		yK		vxZ		vyZ		vxK		vyK\n");
	else if (n == 3)
		fprintf(wyniki, "#	t		xZ		yZ		zZ		xK		yK		zK		vxZ		vyZ		vzZ		vxK		vyK		vzK\n");
	else if (n < 2) {
		printf("Wymiar musi być conajmniej równy 2");
		return;
	}

	fprintf(calki, "t		Ec		Pc[x]		Pc[y]		aV\n");

	force(FZ, rZ, rK, MZ, MK, n);
	force(FK, rK, rZ, MK, MZ, n);

	for (l = 0; l < n; l++) {
		rZPlus[l] = rZ[l] + vZ[l] * dt + FZ[l] / (2 * MZ) * dt * dt;
		rKPlus[l] = rK[l] + vK[l] * dt + FK[l] / (2 * MK) * dt * dt;
	}

	while (t < tN + dt) {
		int j;

		x[i] = rK[0] - rZ[0];
		y[i] = rK[1] - rZ[1];

		// przygotowanie pod wyznaczanie parametrów elipsy
		if (x[i] > xmax) {
			xmax = x[i];
			yXmax = y[i];
		}
		else if (x[i] < xmin) {
			xmin = x[i];
			yXmin = y[i];
		}
		if (y[i] > ymax) {
			ymax = y[i];
			xYmax = x[i];
		}
		else if (y[i] < ymin) {
			ymin = y[i];
			xYmin = x[i];
			check++;
		}

		if (check > 0 && y[i] <= 0) {
			time = t;
		}
		// koniec

		// wyznaczanie energii, pędu i prędkości polowej
		Ec = Ek(MZ, vZ, n) + Ek(MK, vK, n) + Ep(MZ, MK, rZ, rK, n);
		momentum(MZ, vZ, PZ, n);
		momentum(MK, vK, PK, n);
		for (j = 0; j < n; j++) {
			Pc[j] = PZ[j] + PK[j];
		}
		aV = arealVel(rK, vK, n);
		fprintf(calki, "%lf	%e	", t, Ec);
		for (j = 0; j < n; j++) fprintf(calki, "%e	", Pc[j]);
		fprintf(calki, "%e\n", aV);
		// koniec wyznaczania 

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
	// Wyznaczanie parametrów elipsy
	fprintf(kepler, "Skrajne punkty:\n");
	fprintf(kepler, "	Prawa: %e km, %e km\n", xmax / 1000, yXmax / 1000);
	fprintf(kepler, "	Lewa: %e km, %e km\n", xmin / 1000, yXmin / 1000);
	fprintf(kepler, "	Gora: %e km, %e km\n", xYmax / 1000, ymax / 1000);
	fprintf(kepler, "	Dol: %e km, %e km\n", xYmin / 1000, ymin / 1000);

	if (xmax > 0 && xmin < 0 && ymax > 0 && ymin < 0) {
		a = (xmax - xmin) / 2;
		b = (ymax - ymin) / 2;
		c = (xYmax + xYmin) / 2;
		e = c / a;

		fprintf(kepler, "Parametry elipsy:\n");
		fprintf(kepler, "	a = %lf km (384400 km)\n", a / 1000);
		fprintf(kepler, "	b = %lf km\n", b / 1000);
		fprintf(kepler, "	c = %lf km\n", c / 1000);
		fprintf(kepler, "	e = %lf (0.0554)\n", e);

		for (i = 0; i <= int(tN / dt); i++) {
			elipEl[i] = pow((x[i] - c), 2) / pow(a, 2) + pow(y[i], 2) / pow(b, 2);
			elipSum += elipEl[i];
		}
		elipSum /= (tN / dt);
		for (i = 0; i <= int(tN / dt); i++) {
			elipRMS += pow((elipEl[i] - elipSum), 2);
		}
		elipRMS = sqrt(elipRMS / int(tN / dt));
		fprintf(kepler, "	x^2/a^2 + y^2/b^2 ?= 1 : %.7lf +- %.7lf\n", elipSum, elipRMS);
	}
	fprintf(kepler, "Okres obiegu >= %lf d (27.321661 d)\n", time / 24 / 3600);
	// koniec wyznaczania parametrów elipsy

	freeArray(x);
	freeArray(y);
	freeArray(elipEl);
	freeArray(Pc);
	freeArray(PZ);
	freeArray(PK);
	freeArray(rZPlus);
	freeArray(rKPlus);
	freeArray(rZMinus);
	freeArray(rKMinus);
	fclose(wyniki);
	fclose(kepler);
	fclose(calki);
}

int main()
{
	int n = 2;		// size
	const double MZ = 5.97219e24;
	const double MK = 7.347673e22;
	double* rZ;
	double* FZ;
	double* vZ;
	double* rK;
	double* FK;
	double* vK;
	double dt = 0.1 * 3600;
	double tN = (double)30 * 24 * 3600;

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
