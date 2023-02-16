#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <string>
#include <Eigen/Dense>

// Integration parameters
const double TMAX = 30;
const double h = 0.001;

// Derived integration parameters
const double h2 = h/2;
const double h6 = h/6;

// Defined constants of the problem
const double beta = 1.0/50.0;
const double g = 1.0;
const double mu = 1000.0; // mu = 1.0, 10.0, 100.0
const double BASIC_REPR_NUMBER = 2.5;

// Initial conditions
const double S0 = 50;
const double R0 = 0;

// Derived constants of the equations
// We fix landa to keep the basic reproduction number constant for all the simulations
const double landa = BASIC_REPR_NUMBER*g/S0*mu/beta;
const double betap = landa/mu*beta;

// Derived initial conditions
const double I0 = mu/landa*0.01;

// Functions
// S -> 0, I -> 1, R -> 2, P -> 3
inline double fS(const Eigen::Vector3d& y) {
	return -betap*y[1]*y[0];
}

inline double fI(const Eigen::Vector3d& y) {
	return (betap*y[0] - g)*y[1];
}

inline double fR(const Eigen::Vector3d& y) {
	return g*y[1];
}

inline double get_P(const Eigen::Vector3d& y) {
	return landa/mu*y[1];
}

int main() {
	unsigned i, n;

	double t = .0;
	double tWrite = .1, tLastWrite = .0;

	Eigen::Vector3d y(S0, I0, R0);
	Eigen::Vector3d aux;
	std::vector< Eigen::Vector3d > k(4);

	std::vector<double (*)(const Eigen::Vector3d&)> f(3);
	f[0] = fS; f[1] = fI; f[2] = fR;

	std::ofstream fout("SIRP_fastslow_approx_" + std::to_string(mu) + ".txt");

	// We print the header and the initial condition
	fout << "#t\tS\tI\tR\tP\n" << t << ' ' << y.transpose() << ' ' << get_P(y) << '\n';

	while (t < TMAX) {
		// We get k¹
		for (i = 0; i < 3; i++) k[0][i] = f[i](y);
		
		// We get k², k³
		for (n = 1; n < 3; n++) {
			aux = y + h2*k[n-1];
			for (i = 0; i < 3; i++) k[n][i] = f[i](aux);
		}

		// We get k⁴
		aux = y + h*k[2];
		for (i = 0; i < 3; i++) k[3][i] = f[i](aux);

		// We update to the new position
		y += h6*(k[0] + 2*k[1] + 2*k[2] + k[3]);

		// We update the time
		t += h;

		// We write things if we have to
		tLastWrite += h;
		if (tLastWrite > tWrite) {
			tLastWrite -= tWrite;
			fout << t << ' ' << y.transpose() << ' ' << get_P(y) << '\n';
		}
	}

	fout.close();
}
