#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;
const double tx = 0.0005; // =h
const double ty = 0.00025; // tau
const double time0 = 1;
const double x0 = 1;
const int Mx = (int) (x0/tx);
const int My = (int) (time0/ty);
const int steps = 100;
const int timefile = My/steps;
double maxy = 0;
double miny = 10000;
double u(double t, double x) {
	//return exp(5.0 - M_PI * M_PI * t) * sin(x * M_PI);
	return M_PI * x * tan(t) * t + 2 * cos(10 * t) * sin(5 * x);
}

double u0(double x) {
	//return exp(5.0) * sin(x * M_PI);
	return  2 * sin(5 * x);
}

double 	fi1(double t) {
	//return 0;
	return 0; // 1

}

double fi2(double t) {
	//return 0;
	return M_PI * tan(t) * t + 2 * cos(10 * t) * sin(5);
}

double psi1(double t) {
	//return t;
	return M_PI * tan(t) * t + 10 * cos(10 * t);
}

double psi2(double t) {
	//return 	-exp(5.0 - M_PI * M_PI * t) * M_PI;
	return M_PI * tan(t) * t + 10 * cos(10 * t) * cos(5);
}


double f(double t, double x) {
	//return 0;
	return M_PI * x * tan(t) + M_PI * x * t /cos(t) / cos(t) + sin(5 * x) * (50 * cos( 10 * t) - 20 * sin(10 * t));
}

const double alpha = 1;


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
			lay2[Mx] = ( 2 * tx * psi2(ty * i) + 4 * lay2[Mx - 1] - lay2[Mx - 2]) /3;
		}
		double *tmp = lay2;
		lay2 = lay1;
		lay1 = tmp;
		if  ((i % timefile == 0) && mode == 2 ) {
			for (int j = 0; j <= Mx; ++j) {
				fout << j * tx << ' ' << lay1[j] << endl;
			}
			fout << endl << endl;
		}
		for(int j = 0; j <=Mx; ++j) {
			if (abs(lay1[j]) > maxy) {
				maxy = abs(lay1[j]);
			}
			if ( lay1[j] < miny) {
				miny = lay1[j];
			}
		}
	}
	if (mode == 1) {
		double abspogr = 0;
		for(int i = 0 ; i <= Mx; ++i) {
			if( abs(u(time0, i * tx) - lay1[i]) > abspogr) {
				abspogr = abs(u(time0, i * tx) - lay1[i]);
			}
		}
		double maxrealfunk = 0;
		for(int i = 0 ; i <= Mx; ++i) {
			if (abs(u(time0, i * tx)) > maxrealfunk) {
				maxrealfunk = abs(u(time0, i * tx));
			}
		}
		double relativepogr = abspogr / maxrealfunk;
		double pogr2 = 0;
		for(int i = 0 ; i <= Mx; ++i) {
			pogr2 += (u(time0, i * tx) - lay1[i]) * (u(time0, i * tx) - lay1[i]);
		}
		pogr2 *= tx;
		pogr2 = sqrt(pogr2);
		double pogr3 = pogr2/maxrealfunk;
		cout << "abspogr : " << abspogr << endl;
		cout << "relativepogr : " << relativepogr << endl;
		cout << "evklid abs : " << pogr2 << endl;
		cout << "evklid relative : " << pogr3 << endl;
	}
	delete[] lay1;
	delete[] lay2;
	return;
}

void implicit(int mode, int leftborder, int rightborder, ofstream &fout) {
	double ksi1;
	double mu1;
	double ksi2;
	double mu2;
	double bi = alpha * alpha * ty / tx / tx;
	double ci = (1 + 2 * bi);
	double ai = bi;

	double alpha1[Mx + 1];
	double beta[Mx +1];
	double lay[Mx + 1];
	for(int i = 0 ; i <= Mx; ++i) {
		lay[i] = u0(i * tx);
	}
	if  (mode == 2 ) {
		for (int j = 0; j <= Mx; ++j) {
			fout << j * tx << ' ' << lay[j] << endl;
		}
		fout << endl << endl;
	}
	for(int i = 1; i <= My; ++i) {
		if (leftborder == 1) {
			ksi1 = 0;
			mu1 = fi1( ty * i);
		} else if (leftborder == 2) {
			ksi1 = 1;
			mu1 = -tx * psi1( ty * i);
		} else {
			cout << "error leftborder" << endl;
			return;
		}
		if (rightborder == 1) {
			ksi2 = 0;
			mu2 = fi2 (ty * i);
		} else if (rightborder == 2) {
			ksi2 = 1;
			mu2 = tx * psi2(ty * i);
		} else {
			cout << "error rightborder" << endl;
			return;
		}
		alpha1[1] = ksi1;
		beta[1] = mu1;
		for (int j = 2 ; j <= Mx; ++j) {
			alpha1[j] = bi/(ci - ai*alpha1[j-1]);
			//double fi = ty * f(i * ty, j * tx) + lay[j - 1];
			beta[j] = (ty * f(i * ty, (j - 1) * tx) + lay[j - 1] + ai * beta[j-1]) / (ci - ai * alpha1[j - 1]); 	
		}
		if (rightborder == 1) {
			lay[Mx] = fi2(ty * i);
		} else {
			lay[Mx] = (mu2 + ksi2 * beta[Mx]) / (1 - ksi2 * alpha1[Mx]); 
		}
		for (int j = Mx - 1; j != -1 ; --j) {
			lay[j] = alpha1[j + 1] * lay[j + 1] + beta[j + 1];
		}


		if  ((i % timefile == 0) && mode == 2 ) {
			for (int j = 0; j <= Mx; ++j) {
				fout << j * tx << ' ' << lay[j] << endl;
			}
			fout << endl << endl;
		}
		for(int j = 0; j <=Mx; ++j) {
			if (abs(lay[j]) > maxy) {
				maxy = abs(lay[j]);
			}
			if ( lay[j] < miny) {
				miny = lay[j];
			}
		}
		
	}

	if (mode == 1) {
		double abspogr = 0;
		for(int i = 0 ; i <= Mx; ++i) {
			if( abs(u(time0, i * tx) - lay[i]) > abspogr) {
				abspogr = abs(u(time0, i * tx) - lay[i]);
			}
		}
		double maxrealfunk = 0;
		for(int i = 0 ; i <= Mx; ++i) {
			if (abs(u(time0, i * tx)) > maxrealfunk) {
				maxrealfunk = abs(u(time0, i * tx));
			}
		}
		double relativepogr = abspogr / maxrealfunk;
		double pogr2 = 0;
		for(int i = 0 ; i <= Mx; ++i) {
			pogr2 += (u(time0, i * tx) - lay[i]) * (u(time0, i * tx) - lay[i]);
		}
		pogr2 *= tx;
		pogr2 = sqrt(pogr2);
		double pogr3 = pogr2/maxrealfunk;
		cout << "abspogr : " << abspogr << endl;
		cout << "relativepogr : " << relativepogr << endl;
		cout << "evklid abs : " << pogr2 << endl;
		cout << "evklid relative : " << pogr3 << endl;
	}
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
		return 1;
	}
	ofstream omain("main.gn"), oplotter("plotter.gn");
	omain << "set xrange [0:1]\nset yrange [" << (int) miny - 1<<": " <<  (int) maxy + 1 << "]\niter = 0\nload\"plotter.gn\"";
	oplotter << "iter = iter + 1\nplot \"out.txt\" i iter u 1:2 w l lt 6 notitle\npause 0.1\nif (iter < " << steps << ") reread\n";
	fout.close();
	return 0;
}