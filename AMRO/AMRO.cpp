#include "stdafx.h"
using namespace std;

//sdf
//int veloX(Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out);
//int veloY(Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out);
//int veloZ(Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out);

int veloXsub(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out);//derivative of e[kx,ky,kz]
int veloYsub(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out);
int veloZsub(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out);

int veloX(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *tempvel, Ipp64f *out);//derivative of eigen[{e[kx,ky,kz],delta},{delta,e[kx+Pi/a,ky+Pi/a,kz]}]
int veloY(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *tempvel, Ipp64f *out);
int veloZ(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *tempvel, Ipp64f *out);

int func(Ipp64f *params, Ipp64f * kz, Ipp64f * kx, Ipp64f * ky, int length, Ipp64f *temp, Ipp64f *out);

int veloXsubD(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out);//derivative of e[kx,ky,kz]
int veloYsubD(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out);
int veloZsubD(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out);

int veloXD(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *tempvel, Ipp64f *out);//derivative of eigen[{e[kx,ky,kz],delta},{delta,e[kx+Pi/a,ky+Pi/a,kz]}]
int veloYD(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *tempvel, Ipp64f *out);
int veloZD(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *tempvel, Ipp64f *out);

int funcD(Ipp64f *params, Ipp64f * kz, Ipp64f * kx, Ipp64f * ky, int length, Ipp64f *temp, Ipp64f *out);

int fx(Ipp64f * field, Ipp64f *vy, Ipp64f *vz, int length, Ipp64f *temp, Ipp64f *out);
int fy(Ipp64f * field, Ipp64f *vx, Ipp64f *vz, int length, Ipp64f *temp, Ipp64f *out);
int fz(Ipp64f * field, Ipp64f *vx, Ipp64f *vy, int length, Ipp64f *temp, Ipp64f *out);


//int taufun(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, int length, Ipp64f *temp, Ipp64f *out);
int taufun(Ipp64f *params, Ipp64f minDos, Ipp64f maxDos, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out, Ipp64f *ones);
int _tmain(int argc, _TCHAR* argv[])
{
	std::clock_t startT;
	startT = std::clock();
	//Ipp64f thetas[numthetas] = {0,3.1415926 };
	const int numthetas = 1;
	Ipp64f thetas[1] = { 1.48353 };

	
	//const int numthetas = 21;
	//Ipp64f thetas[numthetas] = { 0., 0.0872665, 0.174533, 0.261799, 0.349066, 0.436332, 0.523599, \
0.610865, 0.698132, 0.785398, 0.872665, 0.959931, 1.0472, 1.13446, \
1.22173, 1.309, 1.39626, 1.48353, 1.5708, 1.65806, 1.74533 };
	Ipp64f *condout = new Ipp64f[numthetas];
	//Ipp64f*condout = new Ipp64f[46];
	//Ipp64f tau = .5;
	//Ipp64f final = 10* params[0];//time final?
	
	Ipp64f field45 = 7.91209; // 45 tesla in appropriate units

	Ipp64f *params = new Ipp64f[12]; //last parameter is phi
	DataExtractor extractor("params.dat");
	DataExtractor extractdat("data.dat");
	params = extractor.getDataArray();
	Ipp64f final = 8 * params[1 - 1];
	long steps = 800;//number of time steps?
	Ipp64f h = final / steps;
	FindFermi Fermi( params);
	int nPoints = Fermi.nPoints;
	cout << nPoints << endl;


	
	
	Ipp64f *starts = new Ipp64f[nPoints * 3];
	Fermi.ReturnStart(starts);
	
	Ipp64f *output = new Ipp64f[nPoints*steps * 4]; //stores evolution of orbit around Fermi surface
	Ipp64f *times = new Ipp64f[steps * nPoints]; //time steps

	//std::clock_t startT;
	Ipp64f duration;

	Ipp64f *field = new Ipp64f[3];
	Ipp64f *vzStorage = new Ipp64f[steps*nPoints];
	Ipp64f *vz0Storage = new Ipp64f[nPoints];
	Ipp64f *DOS = new Ipp64f[nPoints];
	Ipp64f *ones = new Ipp64f[nPoints];//for inverting
	ippsSet_64f(1, ones, nPoints);
	Ipp64f *taus = new Ipp64f[nPoints];//phi dependent taus

	Ipp64f *data= new Ipp64f[numthetas];
	data = extractdat.getDataArray();
	Ipp64f residual=0;
	Ipp64f *vx = new Ipp64f[nPoints];
	Ipp64f *vy = new Ipp64f[nPoints];
	Ipp64f *vz = new Ipp64f[nPoints];

	Ipp64f *argx = new Ipp64f[nPoints];
	Ipp64f *argy = new Ipp64f[nPoints];
	Ipp64f *argz = new Ipp64f[nPoints];

	Ipp64f *tempx = new Ipp64f[16 * nPoints];
	Ipp64f *tempy = new Ipp64f[16 * nPoints];
	Ipp64f *tempz = new Ipp64f[16 * nPoints];



	Ipp64f *k1x = new Ipp64f[nPoints];
	Ipp64f *k1y = new Ipp64f[nPoints];
	Ipp64f *k1z = new Ipp64f[nPoints];
	Ipp64f *k2x = new Ipp64f[nPoints];
	Ipp64f *k2y = new Ipp64f[nPoints];
	Ipp64f *k2z = new Ipp64f[nPoints];
	Ipp64f *k3x = new Ipp64f[nPoints];
	Ipp64f *k3y = new Ipp64f[nPoints];
	Ipp64f *k3z = new Ipp64f[nPoints];
	Ipp64f *k4x = new Ipp64f[nPoints];
	Ipp64f *k4y = new Ipp64f[nPoints];
	Ipp64f *k4z = new Ipp64f[nPoints];//what are these? why do we need k1-4?
	Ipp64f *exptau = new Ipp64f[steps * nPoints];
	Ipp64f total = 0;
	Ipp64f *tempfunc = new Ipp64f[16*nPoints];
	//Ipp64f tau = params[1-1];
	//startT = std::clock();
	Ipp64f *minDos = new Ipp64f;
	Ipp64f *maxDos = new Ipp64f;
	
	for (int j = 0; j < nPoints; j++) { //initialize Fermi surface to starting grid, only needs to be done once
		output[0 * nPoints + j] = starts[j * 3 + 0];
		output[1 * nPoints + j] = starts[j * 3 + 1];
		output[2 * nPoints + j] = starts[j * 3 + 2];
	}
	
	ippsCopy_64f(&output[nPoints * (0)], argx, nPoints);//initial velocities for DOS calc;
	ippsCopy_64f(&output[nPoints * (1)], argy, nPoints);
	ippsCopy_64f(&output[nPoints * (2)], argz, nPoints);
	/*
	veloX(params, argx, argy, argz, nPoints, tempx,  tempfunc, vx); //velocities for DOS are stored in vx, vy, and vz buffers.
	veloY(params, argx, argy, argz, nPoints, tempy,  tempfunc, vy);
	veloZ(params, argx, argy, argz, nPoints, tempz,  tempfunc, vz);
	*/

	veloXD(params, argx, argy, argz, nPoints, tempx, tempfunc, vx); //velocities for DOS are stored in vx, vy, and vz buffers.
	veloYD(params, argx, argy, argz, nPoints, tempy, tempfunc, vy);
	veloZD(params, argx, argy, argz, nPoints, tempz, tempfunc, vz);


	ippsSqr_64f_I(vx, nPoints);//in-place square of velocities
	ippsSqr_64f_I(vy, nPoints);
	ippsSqr_64f_I(vz, nPoints);

	ippsAdd_64f(vx, vy, tempx, nPoints);//add all square velocities
	ippsAdd_64f_I(vz, tempx, nPoints);
	ippsSqrt_64f_I(tempx, nPoints);//square root
	ippsDiv_64f(tempx, ones, DOS, nPoints);
	//ippsMin_64f(DOS, nPoints, Dos0);
	ippsMax_64f(DOS, nPoints, maxDos);
	ippsMin_64f(DOS, nPoints, minDos);
	for (int th = 0; th < numthetas; th++) {

		//for (int p = 0; p < steps; p++) { //re-initialize times SLOW STEP CREATE TEMP VARIABLE
			//times[p] = -p;//why -p?
		//	ippsSet_64f(-p, &times[nPoints * p],nPoints);
		//}
		ippsSet_64f(-1, times, steps*nPoints);
		ippsMulC_64f_I(h, times, steps*nPoints);//time stamps

		//field[0] = field45 * sin(thetas[th]) * cos(params[9]);  //set field(phi == 0?)
		//field[1] = field45 * sin(thetas[th]) * sin(params[9]);
		field[0] = field45 * sin(thetas[th]) * cos(params[10]);  //for testing DOS tau
		field[1] = field45 * sin(thetas[th]) * sin(params[10]);
		field[2] = field45 * cos(thetas[th]);
		//for (int z = 0; z < nPoints; z++) {
		//	if (th == 0 ) {
		//		cout << setprecision(20) << times[nPoints*2 + z] << " " << endl;
		//	}
		//}


		for (int i = 1; i < steps; i++) {

			ippsCopy_64f(&output[nPoints * (3 * (i - 1) + 0)], argx, nPoints);//copy arguments for k1;
			ippsCopy_64f(&output[nPoints * (3 * (i - 1) + 1)], argy, nPoints);
			ippsCopy_64f(&output[nPoints * (3 * (i - 1) + 2)], argz, nPoints);
			/*
			veloX(params, argx, argy, argz, nPoints, tempx, tempfunc, vx); //calculate velocities;
			veloY(params, argx, argy, argz, nPoints, tempy, tempfunc, vy);
			veloZ(params, argx, argy, argz, nPoints, tempz, tempfunc, vz);
			*/

			veloXD(params, argx, argy, argz, nPoints, tempx, tempfunc, vx); //calculate velocities;
			veloYD(params, argx, argy, argz, nPoints, tempy, tempfunc, vy);
			veloZD(params, argx, argy, argz, nPoints, tempz, tempfunc, vz);


			ippsCopy_64f(vz, &vzStorage[nPoints * (i - 1)], nPoints);//store vz for conductivity later

			taufun(params, *minDos, *maxDos, argx, argy, argz, nPoints, tempx, taus, ones);// calculate k dependent tau
			ippsDiv_64f_I(taus, &times[nPoints * (i - 1)], nPoints);
			//ippsDivC_64f_I(tau, &times[nPoints * (i - 1)], nPoints);
			//ippsExp_64f_I(&times[nPoints * (i-1)], nPoints);
			//ippsMulC_64f_I((1E-12) * h, &times[nPoints * (i-1)], nPoints);

			//taufun(params, argx, argy, nPoints, tempx, taus);
			//for (int z = 0; z < nPoints; z++) {
			//	if (th == 0 && i == 30) {
			//		cout << setprecision(20) << taus[z] << " " << endl;
			//	}
			//}
			/*for (int z = 0; z < nPoints; z++) {
			if (th == 1 && i == 1) {
			cout << setprecision(20)<< argx[z] << " " << endl;
			}
			}*/


			fx(field, vy, vz, nPoints, tempx, k1x); //calculate evolution in k and store in k1
			fy(field, vx, vz, nPoints, tempy, k1y);
			fz(field, vx, vy, nPoints, tempz, k1z);

			ippsMulC_64f(k1x, h / 2, tempx, nPoints); //prep evolved k step for k2
			ippsMulC_64f(k1y, h / 2, tempy, nPoints);
			ippsMulC_64f(k1z, h / 2, tempz, nPoints);
			ippsAdd_64f(tempx, &output[nPoints * (3 * (i - 1) + 0)], argx, nPoints); //add step to previous k point, load into arguments for k2;
			ippsAdd_64f(tempy, &output[nPoints * (3 * (i - 1) + 1)], argy, nPoints);
			ippsAdd_64f(tempz, &output[nPoints * (3 * (i - 1) + 2)], argz, nPoints);
			/*
			veloX(params, argx, argy, argz, nPoints, tempx, tempfunc, vx); //calculate velocities;
			veloY(params, argx, argy, argz, nPoints, tempy, tempfunc, vy);
			veloZ(params, argx, argy, argz, nPoints, tempz, tempfunc, vz);
			*/
			veloXD(params, argx, argy, argz, nPoints, tempx, tempfunc, vx); //calculate velocities;
			veloYD(params, argx, argy, argz, nPoints, tempy, tempfunc, vy);
			veloZD(params, argx, argy, argz, nPoints, tempz, tempfunc, vz);
			
			fx(field, vy, vz, nPoints, tempx, k2x); //calculate evolution in k and store in k2
			fy(field, vx, vz, nPoints, tempy, k2y);
			fz(field, vx, vy, nPoints, tempz, k2z);

			ippsMulC_64f(k2x, h / 2, tempx, nPoints); //prep evolved k step for k3
			ippsMulC_64f(k2y, h / 2, tempy, nPoints);
			ippsMulC_64f(k2z, h / 2, tempz, nPoints);
			ippsAdd_64f(tempx, &output[nPoints * (3 * (i - 1) + 0)], argx, nPoints); //add step to previous k point, load into arguments for k3;
			ippsAdd_64f(tempy, &output[nPoints * (3 * (i - 1) + 1)], argy, nPoints);
			ippsAdd_64f(tempz, &output[nPoints * (3 * (i - 1) + 2)], argz, nPoints);
			/*
			veloX(params, argx, argy, argz, nPoints, tempx, tempfunc, vx); //calculate velocities;
			veloY(params, argx, argy, argz, nPoints, tempy, tempfunc, vy);
			veloZ(params, argx, argy, argz, nPoints, tempz, tempfunc, vz);
			*/
			veloXD(params, argx, argy, argz, nPoints, tempx, tempfunc, vx); //calculate velocities;
			veloYD(params, argx, argy, argz, nPoints, tempy, tempfunc, vy);
			veloZD(params, argx, argy, argz, nPoints, tempz, tempfunc, vz);
			
			fx(field, vy, vz, nPoints, tempx, k3x); //calculate evolution in k and store in k3
			fy(field, vx, vz, nPoints, tempy, k3y);
			fz(field, vx, vy, nPoints, tempz, k3z);

			ippsMulC_64f(k3x, h, tempx, nPoints); //prep evolved k step for k4
			ippsMulC_64f(k3y, h, tempy, nPoints);
			ippsMulC_64f(k3z, h, tempz, nPoints);
			ippsAdd_64f(tempx, &output[nPoints * (3 * (i - 1) + 0)], argx, nPoints); //add step to previous k point, load into arguments for k4;
			ippsAdd_64f(tempy, &output[nPoints * (3 * (i - 1) + 1)], argy, nPoints);
			ippsAdd_64f(tempz, &output[nPoints * (3 * (i - 1) + 2)], argz, nPoints);
			/*
			veloX(params, argx, argy, argz, nPoints, tempx, tempfunc, vx); //calculate velocities;
			veloY(params, argx, argy, argz, nPoints, tempy, tempfunc, vy);
			veloZ(params, argx, argy, argz, nPoints, tempz, tempfunc, vz);
			*/
			veloXD(params, argx, argy, argz, nPoints, tempx, tempfunc, vx); //calculate velocities;
			veloYD(params, argx, argy, argz, nPoints, tempy, tempfunc, vy);
			veloZD(params, argx, argy, argz, nPoints, tempz, tempfunc, vz);
			
			fx(field, vy, vz, nPoints, tempx, k4x); //calculate evolution in k and store in k4
			fy(field, vx, vz, nPoints, tempy, k4y);
			fz(field, vx, vy, nPoints, tempz, k4z);

			ippsMulC_64f_I(2, k2x, nPoints); //scale k2
			ippsMulC_64f_I(2, k2y, nPoints);
			ippsMulC_64f_I(2, k2z, nPoints);
			ippsMulC_64f_I(2, k3x, nPoints); //scale k3
			ippsMulC_64f_I(2, k3y, nPoints);
			ippsMulC_64f_I(2, k3z, nPoints);

			ippsAdd_64f(k1x, k2x, tempx, nPoints); //add k1 + k2 to temp
			ippsAdd_64f(k1y, k2y, tempy, nPoints);
			ippsAdd_64f(k1z, k2z, tempz, nPoints);

			ippsAdd_64f_I(k3x, tempx, nPoints); //add in k3
			ippsAdd_64f_I(k3y, tempy, nPoints);
			ippsAdd_64f_I(k3z, tempz, nPoints);

			ippsAdd_64f_I(k4x, tempx, nPoints); //add in k4
			ippsAdd_64f_I(k4y, tempy, nPoints);
			ippsAdd_64f_I(k4z, tempz, nPoints);

			ippsMulC_64f_I(h / 6, tempx, nPoints); //scale the entire sum
			ippsMulC_64f_I(h / 6, tempy, nPoints); //scale the entire sum
			ippsMulC_64f_I(h / 6, tempz, nPoints); //scale the entire sum

			ippsAdd_64f(&output[nPoints * (3 * (i - 1) + 0)], tempx, &output[nPoints * (3 * i + 0)], nPoints); //add sum to previous output and store
			ippsAdd_64f(&output[nPoints * (3 * (i - 1) + 1)], tempy, &output[nPoints * (3 * i + 1)], nPoints);
			ippsAdd_64f(&output[nPoints * (3 * (i - 1) + 2)], tempz, &output[nPoints * (3 * i + 2)], nPoints);
		}


		ippsCopy_64f(&output[nPoints * (3 * (steps - 1) + 0)], argx, nPoints);//get velocity for last point
		ippsCopy_64f(&output[nPoints * (3 * (steps - 1) + 1)], argy, nPoints);
		ippsCopy_64f(&output[nPoints * (3 * (steps - 1) + 2)], argz, nPoints);
		
		//veloZ(params, argx, argy, argz, nPoints, tempz, tempfunc, vz);
		veloZD(params, argx, argy, argz, nPoints, tempz, tempfunc, vz);
		ippsCopy_64f(vz, &vzStorage[nPoints * (steps - 1)], nPoints);

		taufun(params, *minDos, *maxDos, argx, argy, argz, nPoints, tempx, taus, ones);// calculate k dependent tau for last point
		ippsDiv_64f_I(taus, &times[nPoints * (steps - 1)], nPoints);
		//ippsDivC_64f_I(tau, &times[nPoints * (steps - 1)], nPoints);
		//ippsExp_64f_I(&times[nPoints * (steps - 1)], nPoints);
		//ippsMulC_64f_I((1E-12)*h, &times[nPoints * (steps - 1)], nPoints);

	/*	ippsCopy_64f(&output[nPoints * (0)], argx, nPoints);//initial velocities for DOS calc;
		ippsCopy_64f(&output[nPoints * (1)], argy, nPoints);
		ippsCopy_64f(&output[nPoints * (2)], argz, nPoints);
		veloX(params, argx, argy, argz, nPoints, tempx, vx); //velocities for DOS are stored in vx, vy, and vz buffers.
		veloY(params, argx, argy, argz, nPoints, tempy, vy);
		veloZ(params, argx, argy, argz, nPoints, tempz, vz);

		ippsSqr_64f_I(vx, nPoints);//in-place square of velocities
		ippsSqr_64f_I(vy, nPoints);
		ippsSqr_64f_I(vz, nPoints);

		ippsAdd_64f(vx, vy, tempx, nPoints);//add all square velocities
		ippsAdd_64f_I(vz, tempx, nPoints);
		ippsSqrt_64f_I(tempx, nPoints);//square root
		ippsDiv_64f(tempx, ones, DOS, nPoints);
		*/

		//need to change!!
		  //ippsDiv_64f_I(taus, times, steps);//exponential stuff, negative tau is taken care of in time
		  //ippsExp_64f_I(times, steps);
		  //ippsMulC_64f_I((1E-12)*h, times, steps);

		ippsCopy_64f(&vzStorage[0], vz0Storage, nPoints);//save initial velocity before exp

		for (int i = 0; i < steps; i++) {

			if (i == 0) { ippsSet_64f(1, &exptau[nPoints * (i)], nPoints); }//start time to be 0 
			else if (i > 1) {
				ippsAdd_64f_I(&times[(i - 2) * nPoints], &times[(i - 1)* nPoints], nPoints);
				ippsExp_64f(&times[nPoints * (i - 1)], &exptau[nPoints * (i)], nPoints);
			}
			else {
				ippsExp_64f(&times[nPoints * (i - 1)], &exptau[nPoints * (i)], nPoints);
			}//integration of (1/tau)

		//	ippsExp_64f(&times[nPoints * (i)], &exptau[nPoints * (i)], nPoints);
			ippsMulC_64f_I((1E-12) * h, &exptau[nPoints * (i)], nPoints);
			ippsMul_64f_I(&exptau[nPoints * (i)], &vzStorage[i * nPoints], nPoints); //multiply velocities by exp time factor

		}

		for (int i = 0; i < (steps - 1); i++) {
			ippsAdd_64f_I(&vzStorage[i*nPoints], &vzStorage[(i + 1)*nPoints], nPoints); //add all and accumulate in last vector
		}

		ippsMul_64f_I(DOS, &vzStorage[(steps - 1)*nPoints], nPoints);
		ippsMul_64f_I(vz0Storage, &vzStorage[(steps - 1)*nPoints], nPoints);//multiply by initial velocities

		ippsSum_64f(&vzStorage[(steps - 1)*nPoints], nPoints, &total);//sum all elements of velocity vector

		condout[th] = total;
	}


	ippsDivC_64f(condout, condout[0], &tempx[2 * numthetas], numthetas);//normalize conductivity
	ippsDiv_64f(&tempx[2 * numthetas], ones, tempx, numthetas);// 1/conductivity to get resistivity 
	ippsSub_64f(tempx, data, &tempx[numthetas], numthetas);// (cal[theta]-dat[theta])
	ippsMul_64f_I(&tempx[numthetas], &tempx[numthetas], numthetas);// (cal[theta]-dat[theta])^2
	ippsSum_64f(&tempx[numthetas], numthetas, &residual);//sum( (cal[theta]-dat[theta])^2)

	duration = (std::clock() - startT) / (Ipp64f)CLOCKS_PER_SEC;
	cout<<residual<<endl;
	cout << "time: " << duration << endl;
	
	ofstream fout;
	fout.open("conductivity.dat");
	fout.precision(20);
	
	//fout << residual << endl;
	
	for (int i = 0; i < numthetas; i++) {

		fout << thetas[i] << '\t' << condout[i] << endl;
		//fout <<  condout[i] << endl;
		//cout << thetas[i] << "\t" << condout[i] << endl;
	}

	/*or (int j = 0; j < 10; ++j) {
		fout << params[j]<<"  ";
	}
	fout << endl;*/
	fout.close();
	//fout.open("Fermitraject.dat");
	//fout.precision(15);
	//cout << nPoints << endl;
	
	for (int j = 0; j < nPoints; ++j) {
		fout.open("Fermitraject"+std::to_string(j)+".dat");
		fout.precision(15);
		for (int i = 1; i < steps; i++) {

			fout << output[nPoints * (3 * (i - 1) + 0) + j] << '\t' << output[nPoints * (3 * (i - 1) + 1) +j ] << '\t' << output[nPoints * (3 * (i - 1) + 2)+j] << endl;
			//cout << thetas[i] << "\t" << condout[i] << endl;
		}
		fout.close();
	}
	
	//for (int i = 0; i < steps; ++i) {
	//	cout << times[i] << " " << endl;
	//}
	cout << "end"<<endl;

	/*ofstream fout2;
	fout2.open("vztimes.dat");
	fout2.precision(15);

	for (int i = 0; i < steps; i++) {

	fout2 << times[i] << endl;
	}

	fout2.close();*/

	/*char a;
	cin >> a;


	delete[] output;


	while (true);*/
	return 0;
}

int func(Ipp64f *params, Ipp64f * kz, Ipp64f * kx, Ipp64f * ky, int length, Ipp64f *temp, Ipp64f *out) {

	ippsMulC_64f(kx, 3.747665940, temp, length);

	vdCos(length, temp, &temp[1 * length]); // cos cos
	ippsMulC_64f(ky, 3.747665940, temp, length);

	vdCos(length, temp, &temp[2 * length]); // cos sin
	ippsMulC_64f(kx, 3.747665940 / 2, temp, length);

	vdCos(length, temp, &temp[3 * length]); // cos cos/2
	ippsMulC_64f(ky, 3.747665940 / 2, temp, length);

	vdCos(length, temp, &temp[4 * length]); // cos sin/2
	ippsMulC_64f(kx, 3.747665940 * 2, temp, length);

	vdCos(length, temp, &temp[5 * length]); // cos 2 cos
	ippsMulC_64f(ky, 3.747665940 * 2, temp, length);

	vdCos(length, temp, &temp[6 * length]); // cos 2 sin
	ippsMulC_64f(kz, 0.5*13.2, &temp[8 * length], length); //kzc / 2
	vdCos(length, &temp[8 * length], &temp[7 * length]); // cos kzc/2
	vdCos(length, &temp[8 * length], &temp[9 * length]); // cos kzc ***made same as above, cos kzc/2

	ippsAdd_64f(&temp[5 * length], &temp[6 * length], temp, length);// param 5
	ippsMulC_64f(temp, -35164.83516*params[5 - 1], out, length);

	ippsMul_64f(&temp[1 * length], &temp[2 * length], temp, length);// param 4
	ippsMulC_64f_I(-35164.83516 * 2 * params[4 - 1], temp, length);
	ippsAdd_64f_I(temp, out, length);

	ippsAdd_64f(&temp[1 * length], &temp[2 * length], temp, length);// param 3
	ippsMulC_64f_I(-35164.83516 * params[3 - 1], temp, length);
	ippsAdd_64f_I(temp, out, length);
	//ippsAddC_64f_I(35164.83516 / 2 * params[2 - 1], out, length);// param 2
	ippsSub_64f(&temp[2 * length], &temp[1 * length], temp, length);// param 6
	ippsSqr_64f_I(temp, length); //square
	ippsMul_64f_I(&temp[3 * length], temp, length); // mult by cos cos/2
	ippsMul_64f_I(&temp[4 * length], temp, length); // mult by cos sin/2
	ippsMul_64f_I(&temp[7 * length], temp, length); // mult by cos  kz/2
	ippsMulC_64f_I(-35164.83516 * params[6 - 1], temp, length);
	ippsMulC_64f_I(-35164.83516 *params[7 - 1], &temp[9 * length], length);
	ippsAdd_64f_I(temp, out, length);
	ippsAdd_64f_I(&temp[9 * length], out, length);

	return 0;
}

int veloXsub(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out)
{
	ippsMulC_64f(kx, 3.74767, temp, length); //term for sin(kx), param3
	vdSin(length, temp, &temp[1 * length]);
	ippsMulC_64f(&temp[1 * length], params[3 - 1] * 11.4215, out, length);

	ippsMulC_64f(kx, 3.74767, temp, length); //term for sin(kx)cos(ky), param4
	vdSin(length, temp, &temp[1 * length]);
	ippsMulC_64f(ky, 3.74767, temp, length);
	vdCos(length, temp, &temp[2 * length]);
	ippsMul_64f_I(&temp[1 * length], &temp[2 * length], length);
	ippsMulC_64f(&temp[2 * length], 22.8429*params[4 - 1], temp, length);
	ippsAdd_64f_I(temp, out, length);

	ippsMulC_64f(kx, 2 * 3.74767, temp, length); //term for sin(2 kx), param5
	vdSin(length, temp, &temp[1 * length]);
	ippsMulC_64f(&temp[1 * length], params[5 - 1] * 11.4215 * 2, temp, length);
	ippsAdd_64f_I(temp, out, length);

	ippsMulC_64f(kx, 3.74767, temp, length); //term for long complicated kz term, param6
	vdSin(length, temp, &temp[1 * length]); // sin kx
	vdCos(length, temp, &temp[2 * length]); // cos kx
	ippsMulC_64f(ky, 3.74767, temp, length);
	vdSin(length, temp, &temp[3 * length]); // sin ky
	vdCos(length, temp, &temp[4 * length]); // cos ky
	ippsMulC_64f(kx, 3.74767 / 2, temp, length);
	vdSin(length, temp, &temp[5 * length]); // sin kx/2
	vdCos(length, temp, &temp[6 * length]); // cos kx/2
	ippsMulC_64f(ky, 3.74767 / 2, temp, length);
	vdSin(length, temp, &temp[7 * length]); // sin ky/2
	vdCos(length, temp, &temp[8 * length]); // cos ky/2
	ippsMulC_64f(kz, 6.6, temp, length);//kz*c/2(c=13.2)
	vdCos(length, temp, &temp[9 * length]); // cos kz/2
	ippsMul_64f(&temp[6 * length], &temp[1 * length], &temp[10 * length], length);//cos kx/2 * sin kx
	ippsMulC_64f_I(65893 * 4, &temp[10 * length], length); // mult by constant
	ippsSub_64f(&temp[4 * length], &temp[2 * length], &temp[11 * length], length);//cos kx - cos ky
	ippsMulC_64f(&temp[11 * length], 65893, &temp[12 * length], length);//mult by constant
	ippsMul_64f(&temp[5 * length], &temp[12 * length], &temp[13 * length], length);// mult by sin kx/2
	ippsAdd_64f(&temp[13 * length], &temp[10 * length], temp, length);//add those two together
	ippsMul_64f_I(&temp[9 * length], temp, length);//mult by cos kz/2
	ippsMul_64f_I(&temp[11 * length], temp, length);//mult by cos kx - cos ky
	ippsMul_64f_I(&temp[8 * length], temp, length);//mult by cos ky/2
	ippsMulC_64f_I(params[6 - 1] * 0.0000866667, temp, length);
	//ippsMulC_64f_I(params[6 - 1], temp, length);
	ippsAdd_64f_I(temp, out, length);

	return 0;
}


int veloYsub(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out) {
	ippsMulC_64f(ky, 3.74767, temp, length); //term for sin(ky), param3
	vdSin(length, temp, &temp[1 * length]);
	ippsMulC_64f(&temp[1 * length], params[3 - 1] * 11.4215, out, length);

	ippsMulC_64f(kx, 3.74767, temp, length); //term for cos(kx)sin(ky), param4
	vdCos(length, temp, &temp[1 * length]);
	ippsMulC_64f(ky, 3.74767, temp, length);
	vdSin(length, temp, &temp[2 * length]);
	ippsMul_64f_I(&temp[1 * length], &temp[2 * length], length);
	ippsMulC_64f(&temp[2 * length], 22.8429*params[4 - 1], temp, length);
	ippsAdd_64f_I(temp, out, length);

	ippsMulC_64f(ky, 2 * 3.74767, temp, length); //term for sin(2 ky), param5
	vdSin(length, temp, &temp[1 * length]);
	ippsMulC_64f(&temp[1 * length], params[5 - 1] * 11.4215 * 2, temp, length);
	ippsAdd_64f_I(temp, out, length);

	ippsMulC_64f(kx, 3.74767, temp, length); //term for long complicated kz term, param6
	vdSin(length, temp, &temp[1 * length]); // sin kx
	vdCos(length, temp, &temp[2 * length]); // cos kx
	ippsMulC_64f(ky, 3.74767, temp, length);
	vdSin(length, temp, &temp[3 * length]); // sin ky
	vdCos(length, temp, &temp[4 * length]); // cos ky
	ippsMulC_64f(kx, 3.74767 / 2, temp, length);
	vdSin(length, temp, &temp[5 * length]); // sin kx/2
	vdCos(length, temp, &temp[6 * length]); // cos kx/2
	ippsMulC_64f(ky, 3.74767 / 2, temp, length);
	vdSin(length, temp, &temp[7 * length]); // sin ky/2
	vdCos(length, temp, &temp[8 * length]); // cos ky/2
	ippsMulC_64f(kz, 6.6, temp, length);//kz*c/2(c=13.2)
	vdCos(length, temp, &temp[9 * length]); // cos kz/2
	ippsMul_64f(&temp[8 * length], &temp[3 * length], &temp[10 * length], length);//cos ky/2 * sin ky
	ippsMulC_64f_I(65893 * 4, &temp[10 * length], length); // mult by constant
	ippsSub_64f(&temp[4 * length], &temp[2 * length], &temp[11 * length], length);//cos kx - cos ky CHECK SIGN
	ippsMulC_64f(&temp[11 * length], 65893, &temp[12 * length], length);//mult by constant
	ippsMul_64f(&temp[7 * length], &temp[12 * length], &temp[13 * length], length);// mult by sin kx/2
	ippsSub_64f(&temp[10 * length], &temp[13 * length], temp, length);//add those two together (neg sign)	
	ippsMul_64f_I(&temp[9 * length], temp, length);//mult by cos kz/2
	ippsMul_64f_I(&temp[11 * length], temp, length);//mult by cos kx - cos ky
	ippsMul_64f_I(&temp[6 * length], temp, length);//mult by cos kx/2	
	ippsMulC_64f_I(params[6 - 1] * 0.0000866667, temp, length);
	//ippsMulC_64f_I(params[6 - 1], temp, length);
	ippsAdd_64f_I(temp, out, length);
	return 0;
	
	
}

int veloZsub(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out) {
	ippsMulC_64f(kx, 3.74767, temp, length); //term for long complicated kz term, param6
	vdCos(length, temp, &temp[2 * length]); // cos kx
	ippsMulC_64f(ky, 3.74767, temp, length);
	vdCos(length, temp, &temp[4 * length]); // cos ky
	ippsMulC_64f(kx, 3.74767 / 2, temp, length);
	vdCos(length, temp, &temp[6 * length]); // cos kx/2
	ippsMulC_64f(ky, 3.74767 / 2, temp, length);
	vdCos(length, temp, &temp[8 * length]); // cos ky/2
	ippsMulC_64f(kz, 6.6, temp, length);
	vdSin(length, temp, &temp[9 * length]); // sin kzc/2
	ippsMulC_64f(kz, 6.6, temp, length); // changed to kz c/2
	vdSin(length, temp, &temp[13 * length]);// sin kzc **fixed to kz c/2

	ippsSub_64f(&temp[4 * length], &temp[2 * length], &temp[11 * length], length);//cos kx - cos ky
	ippsSqr_64f_I(&temp[11 * length], length);//square it
	ippsMul_64f_I(&temp[9 * length], &temp[11 * length], length);// times sin kz/2
	ippsMul_64f_I(&temp[8 * length], &temp[11 * length], length);// times cos ky/2
	ippsMul_64f_I(&temp[6 * length], &temp[11 * length], length);// times cos ky/2
	ippsMulC_64f(&temp[11 * length], params[6 - 1] * 10.0571*2, out, length);
	ippsMulC_64f(&temp[13 * length], params[7 - 1] * 10.0571*2, &temp[12 * length], length);//h7 term ***removed factor of two when switching to cos kz from cos kz/2
	ippsAdd_64f_I(&temp[12 * length], out, length);


	return 0;

	
}

int veloX(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f* tempvel, Ipp64f *out) {
	veloXsub(params, kx, ky, kz, length, temp, tempvel);// (9.1*10^-31) 10^-20/(10^-12)*1/\[HBar]*D[e(kx,ky,kz),kx]
	func(params, kz, kx, ky, length, temp, &tempvel[1 * length]);//e(kx,ky,kz)
	ippsAddC_64f(kx, 3.1415926535 / 3.747665940, &tempvel[2 * length], length);//kx+Pi/a
	ippsAddC_64f(ky, 3.1415926535 / 3.747665940, &tempvel[3 * length], length);// ky+Pi/b
	veloXsub(params, &tempvel[2 * length], &tempvel[3 * length], kz, length, temp, &tempvel[4 * length]); //(9.1*10^-31) 10^-20/(10^-12)*1/\[HBar]*D[e(kx+Pi/a, ky+Pi/b, kz), kx]
	func(params, kz, &tempvel[2 * length], &tempvel[3 * length], length, temp, &tempvel[5 * length]);//e(kx+Pi/a,ky+Pi/b,kz)
	ippsSub_64f(&tempvel[1 * length], &tempvel[5 * length], &tempvel[6 * length], length);//(e(kx+Pi/a,ky+Pi/b,kz)-e(kx,ky,kz))
	ippsSub_64f(tempvel, &tempvel[4 * length], &tempvel[7 * length], length);//(D[e(kx+Pi/2a, ky+Pi/2b, kz), kx]- D[e(kx,ky,kz),kx])
	ippsMul_64f(&tempvel[6 * length], &tempvel[6 * length], &tempvel[8 * length], length);//(e(kx+Pi/a,ky+Pi/b,kz)-e(kx,ky,kz))^2
	ippsAddC_64f_I(params[10 - 1] * params[10 - 1] * (-35164.83516)*(-35164.83516), &tempvel[8 * length], length);//4*delta^2+(e-ePi)^2
	ippsSqrt_64f_I(&tempvel[8 * length], length);//Sqrt(4*delta^2+(e-ePi)^2)
	ippsMul_64f(&tempvel[6 * length], &tempvel[7 * length], &tempvel[9 * length], length);//(e(kx+Pi/2a,ky+Pi/2b,kz)-e(kx,ky,kz))*(D[e(kx+Pi/2a, ky+Pi/2b, kz), kx]- D[e(kx,ky,kz),kx])
	ippsDiv_64f_I(&tempvel[8 * length], &tempvel[9 * length], length);//(-)(-)/sqrt()
	//ippsAdd_64f_I(tempvel,&tempvel[9 * length],length);
	ippsSub_64f_I(&tempvel[9 * length], tempvel, length);// D[e(kx, ky, kz), kx]- (-)(-) / sqrt()**********(Sub -> func2)(Add->func1)
	ippsAdd_64f_I(&tempvel[4 * length], tempvel, length); //D[e(kx, ky, kz), kx] - (-)(-) / sqrt()+/D[e(kx+Pi/2a, ky+Pi/2b, kz), kx]
	ippsMulC_64f(tempvel, 0.5, out, length);//0.5*()
	return 0;
}
int veloY(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f* tempvel, Ipp64f *out) {
	veloYsub(params, kx, ky, kz, length, temp, tempvel);// D[e(kx,ky,kz),ky]
	func(params, kz, kx, ky, length, temp, &tempvel[1 * length]);//e(kx,ky,kz)
	ippsAddC_64f(kx, 3.1415926535 / 3.747665940, &tempvel[2 * length], length);//kx+Pi/a
	ippsAddC_64f(ky, 3.1415926535 / 3.747665940, &tempvel[3 * length], length);// ky+Pi/b
	veloYsub(params, &tempvel[2 * length], &tempvel[3 * length], kz, length, temp, &tempvel[4 * length]); //D[e(kx+Pi/2a, ky+Pi/2b, kz), ky]
	func(params, kz, &tempvel[2 * length], &tempvel[3 * length], length, temp, &tempvel[5 * length]);//e(kx+Pi/2a,ky+Pi/2b,kz)
	ippsSub_64f(&tempvel[1 * length], &tempvel[5 * length], &tempvel[6 * length], length);//(e(kx+Pi/2a,ky+Pi/2b,kz)-e(kx,ky,kz))
	ippsSub_64f(tempvel, &tempvel[4 * length], &tempvel[7 * length], length);//(D[e(kx+Pi/2a, ky+Pi/2b, kz), ky]- D[e(kx,ky,kz),ky])
	ippsMul_64f(&tempvel[6 * length], &tempvel[6 * length], &tempvel[8 * length], length);//(e(kx+Pi/2a,ky+Pi/2b,kz)-e(kx,ky,kz))^2
	ippsAddC_64f_I( params[10 - 1] * params[10 - 1] * (-35164.83516)*(-35164.83516), &tempvel[8 * length], length);//4*delta^2+(e-ePi)^2
	ippsSqrt_64f_I(&tempvel[8 * length], length);//Sqrt(4*delta^2+(e-ePi)^2)
	ippsMul_64f(&tempvel[6 * length], &tempvel[7 * length], &tempvel[9 * length], length);//(e(kx+Pi/2a,ky+Pi/2b,kz)-e(kx,ky,kz))*(D[e(kx+Pi/2a, ky+Pi/2b, kz), ky]- D[e(kx,ky,kz),ky])
	ippsDiv_64f_I(&tempvel[8 * length], &tempvel[9 * length], length);//(-)(-)/sqrt()
	ippsSub_64f_I(&tempvel[9 * length], tempvel, length);// D[e(kx, ky, kz), ky]- (-)(-) / sqrt()**********(Sub -> func2)(Add->func1)
	ippsAdd_64f_I(&tempvel[4 * length], tempvel, length); //D[e(kx, ky, kz), ky] - (-)(-) / sqrt()+/D[e(kx+Pi/2a, ky+Pi/2b, kz), ky]
	ippsMulC_64f(tempvel, 0.5, out, length);//0.5*()
	return 0;
}
int veloZ(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f* tempvel, Ipp64f *out) {
	veloZsub(params, kx, ky, kz, length, temp, tempvel);// D[e(kx,ky,kz),ky]
	func(params, kz, kx, ky, length, temp, &tempvel[1 * length]);//e(kx,ky,kz)
	ippsAddC_64f(kx, 3.1415926535 / 3.747665940, &tempvel[2 * length], length);//kx+Pi/a
	ippsAddC_64f(ky, 3.1415926535 / 3.747665940, &tempvel[3 * length], length);// ky+Pi/b
	veloZsub(params, &tempvel[2 * length], &tempvel[3 * length], kz, length, temp, &tempvel[4 * length]); //D[e(kx+Pi/2a, ky+Pi/2b, kz), ky]
	func(params, kz, &tempvel[2 * length], &tempvel[3 * length], length, temp, &tempvel[5 * length]);//e(kx+Pi/2a,ky+Pi/2b,kz)
	ippsSub_64f(&tempvel[1 * length], &tempvel[5 * length], &tempvel[6 * length], length);//(e(kx+Pi/2a,ky+Pi/2b,kz)-e(kx,ky,kz))
	ippsSub_64f(tempvel, &tempvel[4 * length], &tempvel[7 * length], length);//(D[e(kx+Pi/2a, ky+Pi/2b, kz), ky]- D[e(kx,ky,kz),ky])
	ippsMul_64f(&tempvel[6 * length], &tempvel[6 * length], &tempvel[8 * length], length);//(e(kx+Pi/2a,ky+Pi/2b,kz)-e(kx,ky,kz))^2
	ippsAddC_64f_I(params[10 - 1] * params[10 - 1] * (-35164.83516)*(-35164.83516), &tempvel[8 * length], length);//4*delta^2+(e-ePi)^2
	ippsSqrt_64f_I(&tempvel[8 * length], length);//Sqrt(4*delta^2+(e-ePi)^2)
	ippsMul_64f(&tempvel[6 * length], &tempvel[7 * length], &tempvel[9 * length], length);//(e(kx+Pi/2a,ky+Pi/2b,kz)-e(kx,ky,kz))*(D[e(kx+Pi/2a, ky+Pi/2b, kz), ky]- D[e(kx,ky,kz),ky])
	ippsDiv_64f_I(&tempvel[8 * length], &tempvel[9 * length], length);//(-)(-)/sqrt()
	ippsSub_64f_I(&tempvel[9 * length], tempvel, length);// D[e(kx, ky, kz), ky]- (-)(-) / sqrt()**********(Sub -> func2)(Add->func1)
	ippsAdd_64f_I(&tempvel[4 * length], tempvel, length); //D[e(kx, ky, kz), ky] - (-)(-) / sqrt()+/D[e(kx+Pi/2a, ky+Pi/2b, kz), ky]
	ippsMulC_64f(tempvel, 0.5, out, length);//0.5*()
	return 0;
}


int funcD(Ipp64f *params, Ipp64f * kz, Ipp64f * kx, Ipp64f * ky, int length, Ipp64f *temp, Ipp64f *out) {

	ippsMulC_64f(kx, 3.747665940, temp, length);

	vdCos(length, temp, &temp[1 * length]); // cos cos
	ippsMulC_64f(ky, 3.747665940, temp, length);

	vdCos(length, temp, &temp[2 * length]); // cos sin
	ippsMulC_64f(kx, 3.747665940 / 2, temp, length);

	vdCos(length, temp, &temp[3 * length]); // cos cos/2
	ippsMulC_64f(ky, 3.747665940 / 2, temp, length);

	vdCos(length, temp, &temp[4 * length]); // cos sin/2
	ippsMulC_64f(kx, 3.747665940 * 2, temp, length);

	vdCos(length, temp, &temp[5 * length]); // cos 2 cos
	ippsMulC_64f(ky, 3.747665940 * 2, temp, length);

	vdCos(length, temp, &temp[6 * length]); // cos 2 sin
	ippsMulC_64f(kz, 0.5*13.2, &temp[8 * length], length); //kzc / 2
	vdCos(length, &temp[8 * length], &temp[7 * length]); // cos kzc/2
	vdCos(length, &temp[8 * length], &temp[9 * length]); // cos kzc ***made same as above, cos kzc/2

	ippsAdd_64f(&temp[5 * length], &temp[6 * length], temp, length);// param 5
	ippsMulC_64f(temp, -35164.83516*params[5 - 1], out, length);

	ippsMul_64f(&temp[1 * length], &temp[2 * length], temp, length);// param 4
	ippsMulC_64f_I(-35164.83516 * 2 * params[4 - 1], temp, length);
	ippsAdd_64f_I(temp, out, length);

	ippsAdd_64f(&temp[1 * length], &temp[2 * length], temp, length);// param 3
	ippsMulC_64f_I(-35164.83516 * params[3 - 1], temp, length);
	ippsAdd_64f_I(temp, out, length);
	//ippsAddC_64f_I(35164.83516 / 2 * params[2 - 1], out, length);// param 2
	ippsSub_64f(&temp[2 * length], &temp[1 * length], temp, length);// param 6
	ippsSqr_64f_I(temp, length); //square
	ippsMul_64f_I(&temp[3 * length], temp, length); // mult by cos cos/2
	ippsMul_64f_I(&temp[4 * length], temp, length); // mult by cos sin/2
	ippsMul_64f_I(&temp[7 * length], temp, length); // mult by cos  kz/2
	ippsMulC_64f_I(-35164.83516 * 0, temp, length);
	ippsMulC_64f_I(-35164.83516 * 0, &temp[9 * length], length);
	ippsAdd_64f_I(temp, out, length);
	ippsAdd_64f_I(&temp[9 * length], out, length);

	return 0;
}

int veloXsubD(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out)
{
	ippsMulC_64f(kx, 3.74767, temp, length); //term for sin(kx), param3
	vdSin(length, temp, &temp[1 * length]);
	ippsMulC_64f(&temp[1 * length], params[3 - 1] * 11.4215, out, length);

	ippsMulC_64f(kx, 3.74767, temp, length); //term for sin(kx)cos(ky), param4
	vdSin(length, temp, &temp[1 * length]);
	ippsMulC_64f(ky, 3.74767, temp, length);
	vdCos(length, temp, &temp[2 * length]);
	ippsMul_64f_I(&temp[1 * length], &temp[2 * length], length);
	ippsMulC_64f(&temp[2 * length], 22.8429*params[4 - 1], temp, length);
	ippsAdd_64f_I(temp, out, length);

	ippsMulC_64f(kx, 2 * 3.74767, temp, length); //term for sin(2 kx), param5
	vdSin(length, temp, &temp[1 * length]);
	ippsMulC_64f(&temp[1 * length], params[5 - 1] * 11.4215 * 2, temp, length);
	ippsAdd_64f_I(temp, out, length);

	ippsMulC_64f(kx, 3.74767, temp, length); //term for long complicated kz term, param6
	vdSin(length, temp, &temp[1 * length]); // sin kx
	vdCos(length, temp, &temp[2 * length]); // cos kx
	ippsMulC_64f(ky, 3.74767, temp, length);
	vdSin(length, temp, &temp[3 * length]); // sin ky
	vdCos(length, temp, &temp[4 * length]); // cos ky
	ippsMulC_64f(kx, 3.74767 / 2, temp, length);
	vdSin(length, temp, &temp[5 * length]); // sin kx/2
	vdCos(length, temp, &temp[6 * length]); // cos kx/2
	ippsMulC_64f(ky, 3.74767 / 2, temp, length);
	vdSin(length, temp, &temp[7 * length]); // sin ky/2
	vdCos(length, temp, &temp[8 * length]); // cos ky/2
	ippsMulC_64f(kz, 6.6, temp, length);//kz*c/2(c=13.2)
	vdCos(length, temp, &temp[9 * length]); // cos kz/2
	ippsMul_64f(&temp[6 * length], &temp[1 * length], &temp[10 * length], length);//cos kx/2 * sin kx
	ippsMulC_64f_I(65893 * 4, &temp[10 * length], length); // mult by constant
	ippsSub_64f(&temp[4 * length], &temp[2 * length], &temp[11 * length], length);//cos kx - cos ky
	ippsMulC_64f(&temp[11 * length], 65893, &temp[12 * length], length);//mult by constant
	ippsMul_64f(&temp[5 * length], &temp[12 * length], &temp[13 * length], length);// mult by sin kx/2
	ippsAdd_64f(&temp[13 * length], &temp[10 * length], temp, length);//add those two together
	ippsMul_64f_I(&temp[9 * length], temp, length);//mult by cos kz/2
	ippsMul_64f_I(&temp[11 * length], temp, length);//mult by cos kx - cos ky
	ippsMul_64f_I(&temp[8 * length], temp, length);//mult by cos ky/2
	ippsMulC_64f_I(params[6 - 1] * 0.0000866667*0, temp, length);
	//ippsMulC_64f_I(params[6 - 1], temp, length);
	ippsAdd_64f_I(temp, out, length);

	return 0;
}


int veloYsubD(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out) {
	ippsMulC_64f(ky, 3.74767, temp, length); //term for sin(ky), param3
	vdSin(length, temp, &temp[1 * length]);
	ippsMulC_64f(&temp[1 * length], params[3 - 1] * 11.4215, out, length);

	ippsMulC_64f(kx, 3.74767, temp, length); //term for cos(kx)sin(ky), param4
	vdCos(length, temp, &temp[1 * length]);
	ippsMulC_64f(ky, 3.74767, temp, length);
	vdSin(length, temp, &temp[2 * length]);
	ippsMul_64f_I(&temp[1 * length], &temp[2 * length], length);
	ippsMulC_64f(&temp[2 * length], 22.8429*params[4 - 1], temp, length);
	ippsAdd_64f_I(temp, out, length);

	ippsMulC_64f(ky, 2 * 3.74767, temp, length); //term for sin(2 ky), param5
	vdSin(length, temp, &temp[1 * length]);
	ippsMulC_64f(&temp[1 * length], params[5 - 1] * 11.4215 * 2, temp, length);
	ippsAdd_64f_I(temp, out, length);

	ippsMulC_64f(kx, 3.74767, temp, length); //term for long complicated kz term, param6
	vdSin(length, temp, &temp[1 * length]); // sin kx
	vdCos(length, temp, &temp[2 * length]); // cos kx
	ippsMulC_64f(ky, 3.74767, temp, length);
	vdSin(length, temp, &temp[3 * length]); // sin ky
	vdCos(length, temp, &temp[4 * length]); // cos ky
	ippsMulC_64f(kx, 3.74767 / 2, temp, length);
	vdSin(length, temp, &temp[5 * length]); // sin kx/2
	vdCos(length, temp, &temp[6 * length]); // cos kx/2
	ippsMulC_64f(ky, 3.74767 / 2, temp, length);
	vdSin(length, temp, &temp[7 * length]); // sin ky/2
	vdCos(length, temp, &temp[8 * length]); // cos ky/2
	ippsMulC_64f(kz, 6.6, temp, length);//kz*c/2(c=13.2)
	vdCos(length, temp, &temp[9 * length]); // cos kz/2
	ippsMul_64f(&temp[8 * length], &temp[3 * length], &temp[10 * length], length);//cos ky/2 * sin ky
	ippsMulC_64f_I(65893 * 4, &temp[10 * length], length); // mult by constant
	ippsSub_64f(&temp[4 * length], &temp[2 * length], &temp[11 * length], length);//cos kx - cos ky CHECK SIGN
	ippsMulC_64f(&temp[11 * length], 65893, &temp[12 * length], length);//mult by constant
	ippsMul_64f(&temp[7 * length], &temp[12 * length], &temp[13 * length], length);// mult by sin kx/2
	ippsSub_64f(&temp[10 * length], &temp[13 * length], temp, length);//add those two together (neg sign)	
	ippsMul_64f_I(&temp[9 * length], temp, length);//mult by cos kz/2
	ippsMul_64f_I(&temp[11 * length], temp, length);//mult by cos kx - cos ky
	ippsMul_64f_I(&temp[6 * length], temp, length);//mult by cos kx/2	
	ippsMulC_64f_I(params[6 - 1] * 0.0000866667*0, temp, length);
	//ippsMulC_64f_I(params[6 - 1], temp, length);
	ippsAdd_64f_I(temp, out, length);
	return 0;


}

int veloZsubD(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out) {
	ippsMulC_64f(kx, 3.74767, temp, length); //term for long complicated kz term, param6
	vdCos(length, temp, &temp[2 * length]); // cos kx
	ippsMulC_64f(ky, 3.74767, temp, length);
	vdCos(length, temp, &temp[4 * length]); // cos ky
	ippsMulC_64f(kx, 3.74767 / 2, temp, length);
	vdCos(length, temp, &temp[6 * length]); // cos kx/2
	ippsMulC_64f(ky, 3.74767 / 2, temp, length);
	vdCos(length, temp, &temp[8 * length]); // cos ky/2
	ippsMulC_64f(kz, 6.6, temp, length);
	vdSin(length, temp, &temp[9 * length]); // sin kzc/2
	ippsMulC_64f(kz, 6.6, temp, length); // changed to kz c/2
	vdSin(length, temp, &temp[13 * length]);// sin kzc **fixed to kz c/2

	ippsSub_64f(&temp[4 * length], &temp[2 * length], &temp[11 * length], length);//cos kx - cos ky
	ippsSqr_64f_I(&temp[11 * length], length);//square it
	ippsMul_64f_I(&temp[9 * length], &temp[11 * length], length);// times sin kz/2
	ippsMul_64f_I(&temp[8 * length], &temp[11 * length], length);// times cos ky/2
	ippsMul_64f_I(&temp[6 * length], &temp[11 * length], length);// times cos ky/2
	ippsMulC_64f(&temp[11 * length], params[6 - 1] * 10.0571 * 2, out, length);
	ippsMulC_64f(&temp[13 * length], params[7 - 1] * 10.0571 * 2, &temp[12 * length], length);//h7 term ***removed factor of two when switching to cos kz from cos kz/2
	ippsAdd_64f_I(&temp[12 * length], out, length);


	return 0;


}

int veloXD(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f* tempvel, Ipp64f *out) {
	veloXsubD(params, kx, ky, kz, length, temp, tempvel);// (9.1*10^-31) 10^-20/(10^-12)*1/\[HBar]*D[e(kx,ky,kz),kx]
	funcD(params, kz, kx, ky, length, temp, &tempvel[1 * length]);//e(kx,ky,kz)
	ippsAddC_64f(kx, 3.1415926535 / 3.747665940, &tempvel[2 * length], length);//kx+Pi/a
	ippsAddC_64f(ky, 3.1415926535 / 3.747665940, &tempvel[3 * length], length);// ky+Pi/b
	veloXsubD(params, &tempvel[2 * length], &tempvel[3 * length], kz, length, temp, &tempvel[4 * length]); //(9.1*10^-31) 10^-20/(10^-12)*1/\[HBar]*D[e(kx+Pi/a, ky+Pi/b, kz), kx]
	funcD(params, kz, &tempvel[2 * length], &tempvel[3 * length], length, temp, &tempvel[5 * length]);//e(kx+Pi/a,ky+Pi/b,kz)
	ippsSub_64f(&tempvel[1 * length], &tempvel[5 * length], &tempvel[6 * length], length);//(e(kx+Pi/a,ky+Pi/b,kz)-e(kx,ky,kz))
	ippsSub_64f(tempvel, &tempvel[4 * length], &tempvel[7 * length], length);//(D[e(kx+Pi/2a, ky+Pi/2b, kz), kx]- D[e(kx,ky,kz),kx])
	ippsMul_64f(&tempvel[6 * length], &tempvel[6 * length], &tempvel[8 * length], length);//(e(kx+Pi/a,ky+Pi/b,kz)-e(kx,ky,kz))^2
	ippsAddC_64f_I(params[10 - 1] * params[10 - 1] * (-35164.83516)*(-35164.83516), &tempvel[8 * length], length);//4*delta^2+(e-ePi)^2
	ippsSqrt_64f_I(&tempvel[8 * length], length);//Sqrt(4*delta^2+(e-ePi)^2)
	ippsMul_64f(&tempvel[6 * length], &tempvel[7 * length], &tempvel[9 * length], length);//(e(kx+Pi/2a,ky+Pi/2b,kz)-e(kx,ky,kz))*(D[e(kx+Pi/2a, ky+Pi/2b, kz), kx]- D[e(kx,ky,kz),kx])
	ippsDiv_64f_I(&tempvel[8 * length], &tempvel[9 * length], length);//(-)(-)/sqrt()
	//ippsAdd_64f_I(tempvel,&tempvel[9 * length],length);
	ippsAdd_64f_I(&tempvel[9 * length], tempvel, length);// D[e(kx, ky, kz), kx]- (-)(-) / sqrt()**********(Sub -> func2)(Add->func1)
	ippsAdd_64f_I(&tempvel[4 * length], tempvel, length); //D[e(kx, ky, kz), kx] - (-)(-) / sqrt()+/D[e(kx+Pi/2a, ky+Pi/2b, kz), kx]
	ippsMulC_64f(tempvel, 0.5, out, length);//0.5*()

	ippsMulC_64f(kx, 3.74767, temp, length); //term for long complicated kz term, param6
	vdSin(length, temp, &temp[1 * length]); // sin kx
	vdCos(length, temp, &temp[2 * length]); // cos kx
	ippsMulC_64f(ky, 3.74767, temp, length);
	vdSin(length, temp, &temp[3 * length]); // sin ky
	vdCos(length, temp, &temp[4 * length]); // cos ky
	ippsMulC_64f(kx, 3.74767 / 2, temp, length);
	vdSin(length, temp, &temp[5 * length]); // sin kx/2
	vdCos(length, temp, &temp[6 * length]); // cos kx/2
	ippsMulC_64f(ky, 3.74767 / 2, temp, length);
	vdSin(length, temp, &temp[7 * length]); // sin ky/2
	vdCos(length, temp, &temp[8 * length]); // cos ky/2
	ippsMulC_64f(kz, 6.6, temp, length);//kz*c/2(c=13.2)
	vdCos(length, temp, &temp[9 * length]); // cos kz/2
	ippsMul_64f(&temp[6 * length], &temp[1 * length], &temp[10 * length], length);//cos kx/2 * sin kx
	ippsMulC_64f_I(65893 * 4, &temp[10 * length], length); // mult by constant
	ippsSub_64f(&temp[4 * length], &temp[2 * length], &temp[11 * length], length);//cos kx - cos ky
	ippsMulC_64f(&temp[11 * length], 65893, &temp[12 * length], length);//mult by constant
	ippsMul_64f(&temp[5 * length], &temp[12 * length], &temp[13 * length], length);// mult by sin kx/2
	ippsAdd_64f(&temp[13 * length], &temp[10 * length], temp, length);//add those two together
	ippsMul_64f_I(&temp[9 * length], temp, length);//mult by cos kz/2
	ippsMul_64f_I(&temp[11 * length], temp, length);//mult by cos kx - cos ky
	ippsMul_64f_I(&temp[8 * length], temp, length);//mult by cos ky/2
	ippsMulC_64f_I(params[6 - 1] * 0.0000866667 , temp, length);

	ippsAdd_64f_I(temp,out,length);

	return 0;
}
int veloYD(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f* tempvel, Ipp64f *out) {
	veloYsubD(params, kx, ky, kz, length, temp, tempvel);// D[e(kx,ky,kz),ky]
	funcD(params, kz, kx, ky, length, temp, &tempvel[1 * length]);//e(kx,ky,kz)
	ippsAddC_64f(kx, 3.1415926535 / 3.747665940, &tempvel[2 * length], length);//kx+Pi/a
	ippsAddC_64f(ky, 3.1415926535 / 3.747665940, &tempvel[3 * length], length);// ky+Pi/b
	veloYsubD(params, &tempvel[2 * length], &tempvel[3 * length], kz, length, temp, &tempvel[4 * length]); //D[e(kx+Pi/2a, ky+Pi/2b, kz), ky]
	funcD(params, kz, &tempvel[2 * length], &tempvel[3 * length], length, temp, &tempvel[5 * length]);//e(kx+Pi/2a,ky+Pi/2b,kz)
	ippsSub_64f(&tempvel[1 * length], &tempvel[5 * length], &tempvel[6 * length], length);//(e(kx+Pi/2a,ky+Pi/2b,kz)-e(kx,ky,kz))
	ippsSub_64f(tempvel, &tempvel[4 * length], &tempvel[7 * length], length);//(D[e(kx+Pi/2a, ky+Pi/2b, kz), ky]- D[e(kx,ky,kz),ky])
	ippsMul_64f(&tempvel[6 * length], &tempvel[6 * length], &tempvel[8 * length], length);//(e(kx+Pi/2a,ky+Pi/2b,kz)-e(kx,ky,kz))^2
	ippsAddC_64f_I(params[10 - 1] * params[10 - 1] * (-35164.83516)*(-35164.83516), &tempvel[8 * length], length);//4*delta^2+(e-ePi)^2
	ippsSqrt_64f_I(&tempvel[8 * length], length);//Sqrt(4*delta^2+(e-ePi)^2)
	ippsMul_64f(&tempvel[6 * length], &tempvel[7 * length], &tempvel[9 * length], length);//(e(kx+Pi/2a,ky+Pi/2b,kz)-e(kx,ky,kz))*(D[e(kx+Pi/2a, ky+Pi/2b, kz), ky]- D[e(kx,ky,kz),ky])
	ippsDiv_64f_I(&tempvel[8 * length], &tempvel[9 * length], length);//(-)(-)/sqrt()
	ippsAdd_64f_I(&tempvel[9 * length], tempvel, length);// D[e(kx, ky, kz), ky]- (-)(-) / sqrt()**********(Sub -> func2)(Add->func1)
	ippsAdd_64f_I(&tempvel[4 * length], tempvel, length); //D[e(kx, ky, kz), ky] - (-)(-) / sqrt()+/D[e(kx+Pi/2a, ky+Pi/2b, kz), ky]
	ippsMulC_64f(tempvel, 0.5, out, length);//0.5*()

	ippsMulC_64f(kx, 3.74767, temp, length); //term for long complicated kz term, param6
	vdSin(length, temp, &temp[1 * length]); // sin kx
	vdCos(length, temp, &temp[2 * length]); // cos kx
	ippsMulC_64f(ky, 3.74767, temp, length);
	vdSin(length, temp, &temp[3 * length]); // sin ky
	vdCos(length, temp, &temp[4 * length]); // cos ky
	ippsMulC_64f(kx, 3.74767 / 2, temp, length);
	vdSin(length, temp, &temp[5 * length]); // sin kx/2
	vdCos(length, temp, &temp[6 * length]); // cos kx/2
	ippsMulC_64f(ky, 3.74767 / 2, temp, length);
	vdSin(length, temp, &temp[7 * length]); // sin ky/2
	vdCos(length, temp, &temp[8 * length]); // cos ky/2
	ippsMulC_64f(kz, 6.6, temp, length);//kz*c/2(c=13.2)
	vdCos(length, temp, &temp[9 * length]); // cos kz/2
	ippsMul_64f(&temp[8 * length], &temp[3 * length], &temp[10 * length], length);//cos ky/2 * sin ky
	ippsMulC_64f_I(65893 * 4, &temp[10 * length], length); // mult by constant
	ippsSub_64f(&temp[4 * length], &temp[2 * length], &temp[11 * length], length);//cos kx - cos ky CHECK SIGN
	ippsMulC_64f(&temp[11 * length], 65893, &temp[12 * length], length);//mult by constant
	ippsMul_64f(&temp[7 * length], &temp[12 * length], &temp[13 * length], length);// mult by sin kx/2
	ippsSub_64f(&temp[10 * length], &temp[13 * length], temp, length);//add those two together (neg sign)	
	ippsMul_64f_I(&temp[9 * length], temp, length);//mult by cos kz/2
	ippsMul_64f_I(&temp[11 * length], temp, length);//mult by cos kx - cos ky
	ippsMul_64f_I(&temp[6 * length], temp, length);//mult by cos kx/2	
	ippsMulC_64f_I(params[6 - 1] * 0.0000866667 , temp, length);

	ippsAdd_64f_I(temp, out, length);


	return 0;
}
int veloZD(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f* tempvel, Ipp64f *out) {
	veloZsubD(params, kx, ky, kz, length, temp, out);// D[e(kx,ky,kz),ky]
	
	return 0;
}

int fx(Ipp64f * field, Ipp64f *vy, Ipp64f *vz, int length, Ipp64f *temp, Ipp64f *out) {
	/*return -1 / (11538.5) * (vy * field[2] - vz*field[1]);*/
	ippsMulC_64f(vy, field[2], temp, length);
	ippsMulC_64f(vz, field[1], &temp[length], length);
	ippsMulC_64f_I(-1, &temp[length], length);
	ippsAdd_64f_I(&temp[length], temp, length);
	ippsMulC_64f(temp, 1 / (11538.5), out, length);//+1 to run back in time

	return 0;

}
int fy(Ipp64f * field, Ipp64f *vx, Ipp64f *vz, int length, Ipp64f *temp, Ipp64f *out) {
	//return -1 / (11538.5) * (vz * field[0] - vx*field[2]);
	ippsMulC_64f(vz, field[0], temp, length);
	ippsMulC_64f(vx, field[2], &temp[length], length);
	ippsMulC_64f_I(-1, &temp[length], length);
	ippsAdd_64f_I(&temp[length], temp, length);
	ippsMulC_64f(temp, 1 / (11538.5), out, length);//+1 to run back in time

	return 0;
}
int fz(Ipp64f * field, Ipp64f *vx, Ipp64f *vy, int length, Ipp64f *temp, Ipp64f *out) {
	//return -1 / (11538.5) * (vx * field[1] - vy*field[0]);
	ippsMulC_64f(vx, field[1], temp, length);
	ippsMulC_64f(vy, field[0], &temp[length], length);
	ippsMulC_64f_I(-1, &temp[length], length);
	ippsAdd_64f_I(&temp[length], temp, length);
	ippsMulC_64f(temp, 1 / (11538.5), out, length);//+1 to run back in time

	return 0;
}
/*
int taufun(Ipp64f *params, Ipp64f *kx, Ipp64f *ky, int length, Ipp64f *temp, Ipp64f *out, Ipp64f *ones) {
	ippsDiv_64f(kx, ky, temp, length);
	//cout << temp[0] << endl;
	ippsAtan_64f_A50(temp, &temp[length], length);
	vdSin(length, &temp[length], &temp[2 * length]);//sin(arctan(ky/kx))
	ippsSqr_64f(&temp[2 * length], &temp[3 * length], length);
	vdCos(length, &temp[length], &temp[4 * length]);//cos(arctan(ky/kx))
	ippsSqr_64f(&temp[4 * length], &temp[5 * length], length);
	ippsSub_64f(&temp[3 * length], &temp[5 * length], &temp[6 * length], length);//sin(arctan(ky/kx))^2-cos(arctan(ky/kx))^2
	ippsSqr_64f_I(&temp[6 * length], length);//(sin(arctan(ky/kx))^2-cos(arctan(ky/kx))^2)
	//cout << temp[2 * length] << endl;
	ippsMul_64f(&temp[2 * length], &temp[4 * length], &temp[7 * length], length);
	ippsMulC_64f_I(2, &temp[7 * length], length);//2*sin(arctan(ky / kx))cos(arctan(ky/kx))
	ippsSqr_64f_I(&temp[7 * length], length);//4*sin(arctan(ky / kx))^2cos(arctan(ky/kx))^2
	ippsSub_64f_I(&temp[7 * length], &temp[6 * length], length);// sin(arctan(ky / kx)) ^ 2 - cos(arctan(ky / kx)) ^ 2-4*sin(arctan(ky / kx))^2cos(arctan(ky/kx))^2
	ippsMulC_64f_I(-( params[7]), &temp[6 * length], length);
	ippsAddC_64f(&temp[6 * length], params[0] , out, length);
	//cout << temp[6 * length] << endl;





	return 0;
}*/
int taufun(Ipp64f *params, Ipp64f minDos, Ipp64f maxDos, Ipp64f *kx, Ipp64f *ky, Ipp64f *kz, int length, Ipp64f *temp, Ipp64f *out, Ipp64f *ones) {
	//Ipp64f *minDos = new Ipp64f;
	/*veloX(params, kx, ky, kz, length, &temp[3 * length], temp); //velocities for DOS are stored in vx, vy, and vz buffers.
	veloY(params, kx, ky, kz, length, &temp[3 * length], &temp[length]);
	veloZ(params, kx, ky, kz, length, &temp[3 * length], &temp[2 * length]);
	*/
	veloX(params, kx, ky, kz, length, &temp[3 * length], &temp[4 * length], temp); //calculate velocities;
	veloY(params, kx, ky, kz, length, &temp[3 * length], &temp[4 * length], &temp[length]);
	veloZ(params, kx, ky, kz, length, &temp[3 * length], &temp[4 * length], &temp[2 * length]);
	
	ippsSqr_64f_I(temp, length);//in-place square of velocities
	ippsSqr_64f_I(&temp[length], length);
	ippsSqr_64f_I(&temp[2 * length], length);

	ippsAdd_64f(temp, &temp[length], &temp[3 * length], length);//add all square velocities
	ippsAdd_64f_I(&temp[2 * length], &temp[3 * length], length);
	ippsSqrt_64f_I(&temp[3 * length], length);//square root
	//ippsMulC_64f_I(1 / temp[3 * length], &temp[3 * length], length);//density of state
	//ippsMulC_64f_I(params[8 - 1], &temp[3 * length], length);
	//ippsAddC_64f(&temp[3 * length], params[1 - 1], out, length);
	ippsDiv_64f(&temp[3 * length], ones, &temp[4 * length], length);//density of state
	ippsAddC_64f_I(-maxDos, &temp[4 * length], length);

	ippsMulC_64f_I((1 / params[8 - 1] - 1 / params[1 - 1]) / (maxDos - minDos), &temp[4 * length], length);
	//	ippsMulC_64f_I(1 / params[1 - 1], &temp[4 * length], length);
		//ippsAddC_64f_I(1/params[8 - 1]+params[10 - 1], &temp[4 * length], length);
	ippsAddC_64f_I(1 / params[8 - 1], &temp[4 * length], length);

	ippsDiv_64f(kx, ky, temp, length);
	ippsAtan_64f_A50(temp, &temp[length], length);
	vdSin(length, &temp[length], &temp[2 * length]);//sin(arctan(ky/kx))
	ippsSqr_64f_I(&temp[2 * length], length);
	vdCos(length, &temp[length], &temp[3 * length]);//cos(arctan(ky/kx))
	ippsSqr_64f_I(&temp[3 * length], length);
	ippsSub_64f_I(&temp[2 * length], &temp[3 * length], length);//sin(arctan(ky/kx))^2-cos(arctan(ky/kx))^2
	ippsMul_64f_I(&temp[3 * length], &temp[3 * length], length); //(sin(arctan(ky / kx)) ^ 2 - cos(arctan(ky / kx)) ^ 2)^2
	ippsMulC_64f_I(params[9 - 1], &temp[3 * length], length);
	ippsAdd_64f_I(&temp[3 * length], &temp[4 * length], length);
	//ippsAddC_64f_I(params[10 - 1], &temp[4 * length], length);

	ippsDiv_64f(&temp[4 * length], ones, out, length);
	//cout << *minDos << endl;
	//delete minDos;






	return 0;
}