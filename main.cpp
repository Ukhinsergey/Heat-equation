#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
const double tx = 0.01; // =h
const double ty = 0.000001; // t
const double time0 = 1;
const double x0 = 1;
const int Mx = (int) (x0/tx);
const int My = (int) (time0/ty);


double u(double t, double x) {
	return 3 * exp(-M_PI*M_PI * t) * sin(2 * M_PI * x) - t * x + t;
}

double u0(double x) {
	return 3 * sin(2*M_PI * x);
}

double 	fi1(double t) {
	return t;
}

double fi2(double t) {
	return 0;
}

double psi1(double t) {
	return t;
}

double psi2(double t) {
	return t;
}


double f(double t, double x) {
	return 1 - x;
}

const double alpha = 0.5;


void explisit(int mode,int leftborder,int rightborder, ofstream &fout) {
	double *lay1 = new double [(int) Mx + 1];
	double *lay2 = new double [(int) My + 1];
	for(int i = 0; i <= Mx; ++i) {
		lay1[i] = u0(i * tx);
	}
	for(int i = 1; i <= My; ++i){
		for (int j = 1; j <= Mx - 1; ++j) {
			lay2[j] = lay1[j] + ty * alpha/tx * alpha/tx * (lay1[j+1] - 2 * lay1[j] + lay1[j - 1]) + ty * f(i * ty, j * tx);
		}  
		if (leftborder == 1) {
			lay2[0] = fi1(ty * i);
		} else if (leftborder == 2) {
			lay2[0] = (4 * lay2[1] - lay2[2] - psi1(ty * i) * 2 * tx) / 3;
		} else {
			cout << "error leftborder" << endl;
			return ;
		}
		if (rightborder == 1) {
			lay2[Mx] = fi2(ty * i);
		} else if (rightborder == 2) {
			lay2[Mx] = ( 2 * tx * psi2(ty * i) + 4 * lay2[Mx - 1] - lay2[Mx - 2]);
		}
		double *tmp = lay2;
		lay2 = lay1;
		lay1 = tmp;
		if  (mode == 2) {
			for (int j = 0; j <= Mx; ++j) {
				fout << j * tx << ' ' << lay1[j] << endl;
			}
			fout << endl << endl;
		}
	}
	if (mode == 1) {
		for (int j = 0; j <= Mx; ++j) {
			cout << j * tx << ' ' << lay1[j] << ' ' << u(time0, j*tx) << endl;
		}
	}
	delete[] lay1;
	delete[] lay2;
	return;
}

void implicit(int mode, int leftborder, int rightborder, ofstream &out) {



}




int main() {
	cout << "enter mode: 1 - test, 2 - visual" << endl;
	int mode ;
	cin >> mode;
	cout << "enter method 1 - explisit, 2 - implicit" << endl;
	int method ;
	cin  >> method;
	int leftborder, rightborder;
	cout << "type of leftborder: ";
	cin >> leftborder;
	cout << "type of rightborder: ";
	cin >> rightborder;

	ofstream fout("out.txt");

	if (method == 1) {
		explisit(mode, leftborder, rightborder, fout);
	} else if (method == 2) {
		implicit(mode, leftborder, rightborder, fout);
	} else {
		cout << "error method" << endl ;
	}
	cout << "hi" ;
	fflush(stdout);
	fout.close();
	cout << endl << M_PI;
	return 0;
}