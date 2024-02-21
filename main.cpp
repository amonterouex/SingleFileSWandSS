#include <iostream>
#include <iomanip>
#include <fstream>
#include <stdio.h> 
#include <unistd.h>
#include <cmath>
#include <chrono>
#include <pthread.h>
#include <omp.h>
#include <complex>

#include "quasi1d.hpp"
#include "iofunctions.hpp"

int main() 
{
    InputData indata;
    if(ReadInputFile("input.dat", &indata)){
        std::cout << "Exiting...\n" << std::endl;
        std::cout << "\n Press enter to exit \n" << std::endl;
        std::cin.get();
        return 1;
    }
    if(CheckData(&indata)){
        std::cout << "Exiting...\n" << std::endl;
        std::cout << "\n Press enter to exit \n" << std::endl;
        std::cin.get();
        return 1;
    }

    //Discretization parameter
    int nsize = indata.nsize;

    //Lo primero es ver si el usuario ha dado presion o temperatura
    Quasi1d f = Quasi1d(nsize, indata.eps, indata.r0, 1.0/indata.temperature);
    if (indata.density>0.0){
        f.SetDensity(indata.density);
        std::cout << std::fixed << std::setprecision(9) << "The pressure of the system has been set to " << f.bp << std::endl;
        std::cout << "The compressiblity factor is " << f.bp/indata.density << std::endl;
        indata.bp = f.bp;
    }
    else{
        f.SetPressure(indata.bp);
        double density = f.Density(indata.bp);
        std::cout << "The density of the system has been set to " << density << std::endl;
        std::cout << "The compressiblity factor is " << f.bp/density << std::endl;
    }
    std::ofstream ofile;

    //Compute equation of state
    if (indata.comValue == 1){
        std::vector<double> pvalues (indata.npoints, 0.0); //vector storing pressure values in the interval
        std::vector<double> dvalues (indata.npoints, 0.0); //vector storing density values for the corresponding pressures in pvalues
        std::vector<double> uvalues (indata.npoints, 0.0); //vector storing density values for the corresponding pressures in pvalues

    //Compute step size of the pressure interval
    double pstep = (indata.xmax-indata.xmin)/(indata.npoints-1);

    //Compute density for each of the pressure values in the interval
    for(int pidx=0;pidx<indata.npoints;pidx++){
        pvalues[pidx] = indata.xmin + pstep * (double)pidx;
        dvalues[pidx] = f.Density(pvalues[pidx]);
        uvalues[pidx] = f.InternalEnergy(pvalues[pidx]);
    }
    //Write results

    ofile.open ("equation_of_state.txt");
    ofile << "#Parameters of the system: eps=" << f.eps << ", r0 = " << f.r0 << ", temperature:" << 1.0/f.beta << std::endl;
    ofile << "bp density Z uex" << std::endl;
    for(int i=0;i<pvalues.size();i++){
            //if (i<=4)  std::cout << std::fixed << std::setprecision(6)<<  "\t" <<dvalues[i] << "\t" << pvalues[i]/dvalues[i] << std::endl;
            ofile << pvalues[i] << "\t" << dvalues[i] << "\t" << pvalues[i]/dvalues[i] <<  "\t" << uvalues[i] <<std::endl;           
    }
    ofile.close();

    }
    else if (indata.comValue==2){    //Compute Gr

        //vectors to store results
        double xi;
        double rstep = (indata.xmax-indata.xmin)/(indata.npoints-1);
        std::vector<double> rvalues(indata.npoints,0.0);
        std::vector<double> Fvalues(indata.npoints,0.0);
        int nmax=0;

        std::vector<double> Grtemp;
        std::vector<double> Grtemp2;
        std::vector<double> gvalues11(indata.npoints,0.0);
        std::vector<double> gvalues12(indata.npoints,0.0);
        std::vector<double> gvalues22(indata.npoints,0.0);
        std::vector<double> gvalues13(indata.npoints,0.0);
    
        auto begin = std::chrono::high_resolution_clock::now();
        double kmax = 0.0;

        std::cout << "Computing..." << std::endl;



        #pragma omp parallel for num_threads(16) schedule(dynamic)
        for(int k=0;k<indata.npoints;k++){
            Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> Qtt;
            Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> Qttm;
            Qtt.resize(f.nsize,f.nsize);
            Qttm.resize(f.nsize,f.nsize);
            xi = indata.xmin + k*rstep;
            rvalues[k] = xi;
            if(xi < 3.0*sqrt(1-indata.eps*indata.eps)){
            //if(0){
            //Siempre como inversa de laplace
                Grtemp = f.Gx(xi, f.bp);
                Fvalues[k] = Grtemp[0];
                gvalues11[k] = Grtemp[1];
                gvalues13[k] = Grtemp[3];
                gvalues12[k] = Grtemp[2];
                gvalues22[k] = Grtemp[4];  
            }
            else{
                Grtemp = f.GxInv(xi, f.bp, 300,40,Qtt, Qttm);
                Fvalues[k] = Grtemp[0];
                gvalues11[k] = Grtemp[1];
                gvalues13[k] = Grtemp[3];
                gvalues12[k] = Grtemp[2];
                gvalues22[k] = Grtemp[4];  
            }         

            //Progress bar
            #pragma omp critical
            ProgressBar(double(k)/(double(indata.npoints)-1),70,&kmax);
            //std::cout << k << std::endl;

        }
        std::cout << std::endl;
        std::cout << "Done!\n" << std::endl;

        //Write results to file
        std::cout << "******************\nWriting results to file " << "output_Gx.txt" << std::endl;
        std::cout << "Formatting:\n x \t G(x)\n\n";
        std::ofstream ofile;
        ofile.open("output_Gx.txt");
        for(int k=0;k<(int)indata.npoints;k++){
            ofile << std::fixed<< std::setprecision(15) << rvalues[k] << "\t" << Fvalues[k] << " " << gvalues11[k] << " " << gvalues12[k]<< " " << gvalues13[k]<< " " << gvalues22[k] << "\n";
        }
        ofile.close();
        std::cout << "******************" << std::endl;
    }
    else if (indata.comValue==3){ //Compute density profile
    //Export density profile
        ofile.open ("density_profile.txt");
        for(int i=0;i<f.xvalues.size();i++){
            ofile << f.y(i) << "\t" << f.xvalues[i] << std::endl;           
        }
        ofile.close();
    }
    else if (indata.comValue==4){ //Compute Zeno and Seno line
    double tmin = 0.01;
    double tmax = 0.5885;
    int tnum = 200;
    double tstep = (tmax-tmin)/tnum;
    for (int ti=0; ti<tnum+1; ti++){
        double tc = tmin + tstep*ti;
        Quasi1d f = Quasi1d(nsize, indata.eps, indata.r0, 1.0/tc);

        std::cout << tc << "\t" << f.Zeno() << "\t" << f.Seno() << std::endl;
    }
    }
    else if (indata.comValue==5){ //Compute Sq
        //vectors to store results
        double xi;
        double rstep = (indata.xmax-indata.xmin)/(indata.npoints-1);
        std::vector<double> qvalues(indata.npoints,0.0);
        std::vector<double> Sq(indata.npoints,0.0);
        int nmax=0;

        double kmax = 0.0;

        std::cout << "Computing..." << std::endl;

        //Fuera del bucle for porque no voy a paralelizar
        Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> Qtt;
        Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> Qttm;
        Qtt.resize(f.nsize,f.nsize);
        Qttm.resize(f.nsize,f.nsize);
        for(int k=0;k<indata.npoints;k++){
            xi = indata.xmin + k*rstep;
            qvalues[k] = xi;
            std::complex<double> Gs1 = f.Gs(std::complex<double>(0.0,xi), f.bp, Qtt, Qttm)[0];
            std::complex<double> Gs2 = f.Gs(std::complex<double>(0.0,-xi), f.bp, Qtt, Qttm)[0];
            Sq[k] = std::real(1.0+f.density*(Gs1+Gs2));

            //Progress bar
            ProgressBar(double(k)/(double(indata.npoints)-1),70,&kmax);
            //std::cout << k << std::endl;

        }
        std::cout << std::endl;
        std::cout << "Done!\n" << std::endl;

        //Write results to file
        std::cout << "******************\nWriting results to file " << "output_Sq.txt" << std::endl;
        std::cout << "Formatting:\n x \t G(x)\n\n";
        std::ofstream ofile;
        ofile.open("output_Sq.txt");
        for(int k=0;k<(int)indata.npoints;k++){
            ofile << std::fixed<< std::setprecision(15) << qvalues[k] << "\t" << Sq[k] << "\n";
        }
        ofile.close();
        std::cout << "******************" << std::endl;
    }
        
/*
    std::string outfiles="y";

    double density;
    Eigen::Matrix<double, Eigen::Dynamic, 1> eigvec;

    std::vector<double> pvalues (pnum, 0.0); //vector storing pressure values in the interval
    std::vector<double> dvalues (pnum, 0.0); //vector storing density values for the corresponding pressures in pvalues

    //Compute step size of the pressure interval
    double pstep = (pmax-pmin)/(pnum-1);

    //Compute density for each of the pressure values in the interval
    for(int pidx=0;pidx<pnum;pidx++){
        pvalues[pidx] = pmin + pstep * (double)pidx;
        computeDensity(eps, pvalues[pidx], 1.0/temperature, nsize, r0, &dvalues[pidx], &eigvec);
    }

    //Write results
    std::cout << "\n************************" << std::endl;
    std::cout << "Epsilon =   " << eps << std::endl;
    std::cout << "Pressure =  [" << pmin << ", " << pmax << "]" << std::endl;
    std::cout << "************************" << std::endl;
    if (outfiles == "n"){
        std::cout << "No output files produced" << std::endl;
    }
    if (outfiles == "y"){
        std::cout << "Writing file equation_of_state.txt..." << std::endl;
        std::cout << "   ...exporting in a column-like format without headers:" << std::endl;
        std::cout << "\t (density) \t (Z)" << std::endl;
        double dy = (eps/(nsize-1));
        std::ofstream ofile;
        ofile.open ("equation_of_state.txt");
        for(int i=0;i<pvalues.size();i++){
            if (i<=4)  std::cout << std::fixed << std::setprecision(6)<<  "\t" <<dvalues[i] << "\t" << pvalues[i]/dvalues[i] << std::endl;
            ofile << dvalues[i] << "\t" << pvalues[i]/dvalues[i] << std::endl;           
        }
        ofile.close();
        std::cout << "    ...done." << std::endl;
        
    }
    std::cout << "************************" << std::endl;
*/
    return 0;
    
}