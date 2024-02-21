#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif 

#ifndef QUASI1D_H_INCLUDED
#define QUASI1D_H_INCLUDED

//Efficient way of computing binomial coefficients
double binomialCoefficients(int n, int k) {
   double C[k+1];
   memset(C, 0, sizeof(C));
   C[0] = 1;
   for (int i = 1; i <= n; i++) {
      for (int j = std::min(i, k); j > 0; j--)
         C[j] = C[j] + C[j-1];
   }
   return C[k];
}

//Class Quasi1d
class Quasi1d {       
  public:             
    int nsize; //size of the discretization (number of species).
    double eps; //excess width of the pore
    double bp; //pressure
    double beta; //beta
    double r0; //corona size
    double amin; //min distance between two particles
    double density;

    double Aval; //A coming from the eigenvalue equation.
    Eigen::Matrix<double, Eigen::Dynamic, 1> xvalues; //vector of molar fractions.

    //int nprocs = omp_get_num_procs();

    //Constructor
    Quasi1d(int nsize0, double eps0, double r00, double beta0){
        nsize = nsize0;
        eps = eps0;
        r0 = r00;
        beta = beta0;
        amin = sqrt(1-eps*eps);

        //Internal stuff
        am.resize(nsize, nsize);
        M.resize(nsize, nsize);
        UpdateMatrixA();

        //Useful matrices for later on
        //Qt.resize(nsize,nsize);
        //Qtm.resize(nsize,nsize);
        
    }

    //Set pressure of the system and update xvalues and Lmean.
    int SetPressure(double bp0){
        //Create matrix M
        //std::cout << "Setting pressure " << bp0 << std::endl;
        for (int i=0; i<nsize; i++){
            for(int j=0; j<nsize; j++){
                M(i,j) = omega(bp0, i, j);
                //std::cout << M(i,j) << std::endl;
            }
        }

        //Solve eigenvalue equation
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>> es;
        es.compute(M);
        xvalues = es.eigenvectors().col(nsize-1);
        Aval = sqrt(1.0/es.eigenvalues()(nsize-1));

        

        //Prepare values
        double sum = 0;
        for(int i = 0; i<nsize; i++){
            xvalues(i) = xvalues(i)*xvalues(i);
            sum += xvalues(i);
        }
        xvalues /= sum;

        //check normalization
        //double norm = 0.0;
        //for(int i = 0; i<nsize; i++){
        //    norm = norm + xvalues(i);
        //}
        //std::cout << "Normalization: " << norm << std::endl;

        //Update internal bp
        this->bp = bp0;
        this->density = Density(bp);
        return 0;

    }

    //If instead of pressure the user wishes to set the density, use this function.
    void SetDensity(double rho){
        //Bolzano method to compute the pressure associated with that density
        double a = 0.000001;
        double b=170;
  
        double c = a;
        while ((b-a) >= 0.00000001){
 
            // Find middle point
            c = (a+b)/2.0;
           
            // Decide the side to repeat the steps
            if ((Density(c)-rho)*(Density(a)-rho) < 0)   b = c;
            else  a = c;
        }

        if (abs(c-170.0)<=0.0001){
            std::cout << "WARNING: required density set too high, numerical errors might occur. Setting new density to " << Density(c) << std::endl;
        }

        //Set the corresponding pressure
        SetPressure(c);


    }

    // Compute density of the system at a certain bp
    double Density(double bp0){
        if (bp0 != bp)  SetPressure(bp0);
        double nsum = 0;
        //Igual puedo optimizar este bucle y hacer la mitad del bucle
        for (int i=0; i<nsize; i++){
            for(int j=0; j<nsize; j++){
                nsum += sqrt(xvalues(i)*xvalues(j))*omegad(bp,i,j)*Aval*Aval;
            }
        }
        return -beta/nsum;
    }

    double Seno(){
        Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> Qtt;
        Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> Qttm;
        Qtt.resize(nsize,nsize);
        Qttm.resize(nsize,nsize);

        std::complex<double> s1(0.0000001,0.0);
        std::complex<double> s2(0.001,0.0);
        //double s1=0.001;
        //double s2=0.01;
        

        //Bolzano method
        double a = 0.0000001;
        double b=10.0;
  
        double c = a;
        while ((b-a) >= 0.00000001){
 
            // Find middle point
            c = (a+b)/2.0;
            //compute slope
            std::complex<double> slopec = (s2*Gs(s2, c, Qtt, Qttm)[0]- s1*Gs(s1, c, Qtt, Qttm)[0])/(s2-s1);
            std::complex<double> slopea = (s2*Gs(s2, a, Qtt, Qttm)[0]- s1*Gs(s1, a, Qtt, Qttm)[0])/(s2-s1);

            //std::cout << c << "\t" <<density << "\t" <<1.0/beta << "\t" << slopea << "\t" << slopec  << std::endl;

            // Decide the side to repeat the steps
            if (std::real((slopec-0.0)*(slopea-0.0)) < 0)   b = c;
            else  a = c;
        }

        if (abs(c-170.0)<=0.0001){
            std::cout << "WARNING: required density set too high, numerical errors might occur. Setting new density to " << Density(c) << std::endl;
        }
        //std::cout << "Zeno value is " << c << std::endl;
        return Density(c);

    }

    double Zeno(){
        double a = 0.000001;
        double b=5.0;
  
        double c = a;
        while ((b-a) >= 0.00000001){
 
            // Find middle point
            c = (a+b)/2.0;

            //
           
            // Decide the side to repeat the steps
            if ((c/Density(c)-1.0)*(a/Density(a)-1.0) < 0)   b = c;
            else  a = c;
        }

        if (abs(c-170.0)<=0.0001){
            std::cout << "WARNING: required density set too high, numerical errors might occur. Setting new density to " << Density(c) << std::endl;
        }
        //std::cout << "Zeno value is " << c << std::endl;
        return c;

        //Set the corresponding pressure
        //SetPressure(c);
    }

    // Compute G(s)
    template<class T>
    std::vector<T> Gs(T s, double bp0, Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> &Qtt,Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> &Qttm){
        if (bp0 != bp && bp0 >1e-10)  SetPressure(bp0);
        std::vector<T> gtotal = std::vector<T>(5,0.0);

        for (int i=0; i<nsize; i++){
            for(int j=0; j<nsize; j++){
                Qtt(i,j) = Aval*Aval*omega<T>(s+bp0,i,j);
                Qttm(i,j) = -Qtt(i,j);
                if(i==j) Qttm(i,j) = 1.0+Qttm(i,j);

            }
        }

        //Invert matrix and all that
        Qtt = Qtt*(Qttm.inverse());
        T gmean = 0;
        double gmeanRE = 0;
        double gmeanIM = 0;

        double density = Density(bp);
        T dummy = 0.0;
        //#pragma omp parallel for reduction(+:gmeanRE,gmeanIM) num_threads(nprocs) private(dummy)
        for (int i=0; i<nsize; i++){
            for(int j=i; j<nsize; j++){
                dummy = sqrt(xvalues(i)*xvalues(j))*Qtt(i,j);
                if (i!=j){
                    gmeanRE += std::real(2.0*dummy);
                    gmeanIM += std::imag(2.0*dummy);
                }
                else{
                    gmeanRE += std::real(dummy);
                    gmeanIM += std::imag(dummy);
                }
                //std::cout << s << " "<< i << " " << j << " " << gmeanRE << " " << gmeanIM << std::endl;
                //Get partials G(r)
                if(i==0 && j==0)          gtotal[1] = Qtt(i,j)/sqrt(xvalues(i)*xvalues(j))/density;
                else if(i==0 && j==(nsize+1)/2) gtotal[2] = Qtt(i,j)/sqrt(xvalues(i)*xvalues(j))/density;
                else if(i==0 && j==nsize-1)    gtotal[3] = Qtt(i,j)/sqrt(xvalues(i)*xvalues(j))/density;
                else if(i==(nsize+1)/2 && j==(nsize+1)/2) gtotal[4] = Qtt(i,j)/sqrt(xvalues(i)*xvalues(j))/density;
            }
        }
        if constexpr (std::is_same_v<T, std::complex<double>>){
            gtotal[0] = std::complex<double>(gmeanRE/density, gmeanIM/density);
        }
        else{
            gtotal[0] = gmeanRE/density;
        }
        return gtotal;
    }

    //G(x) calculated from numerical Laplace inverse of G(s)
    std::vector<double> GxInv(double t, double bp0, int ntr, int meuler, Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> &Qtt,Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> &Qttm){
        double aa=30.0;
        double hh = M_PI/t;

        double uu = exp(aa/2.0)/t;
        double xx = aa/(2.0*t);
        //Sum[]
        std::vector<std::complex<double>> Gstemp = Gs<std::complex<double>>(xx,bp0,Qtt, Qttm);
        std::complex<double> suma = Gstemp[0]/2.0;
        std::complex<double> suma1 = Gstemp[1]/2.0;
        std::complex<double> suma2 = Gstemp[2]/2.0;
        std::complex<double> suma3 = Gstemp[3]/2.0;
        std::complex<double> suma4 = Gstemp[4]/2.0;

        std::complex<double> dummy;
        for (int i=1;i<ntr+1;i++){
            Gstemp = Gs<std::complex<double>>(std::complex<double>(xx, i*hh),bp0,Qtt, Qttm);
            //std::cout << std::complex<double>(xx, i*hh) << " " << bp0 << " " << Gstemp[0] << " " << Gstemp[1]<< " " << Gstemp[2]<< std::endl;
            dummy = std::complex<double>(pow(-1,i),0);
            suma += dummy * Gstemp[0];
            suma1 += dummy * Gstemp[1];
            suma2 += dummy * Gstemp[2];
            suma3 += dummy * Gstemp[3];
            suma4 += dummy * Gstemp[4];
        }

        //#su[]
        std::vector<std::complex<double>> su(meuler+2,std::complex<double>(0,0));
        std::vector<std::complex<double>> su1(meuler+2,std::complex<double>(0,0));
        std::vector<std::complex<double>> su2(meuler+2,std::complex<double>(0,0));
        std::vector<std::complex<double>> su3(meuler+2,std::complex<double>(0,0));
        std::vector<std::complex<double>> su4(meuler+2,std::complex<double>(0,0));
        su[0] = suma;
        su1[0] = suma1;
        su2[0] = suma2;
        su3[0] = suma3;
        su4[0] = suma4;

        //#avgsu
        //std::complex<double> avgsu(0,0);
        std::complex<double> avgsu(0,0);
        std::complex<double> avgsu1(0,0);
        std::complex<double> avgsu2(0,0);
        std::complex<double> avgsu3(0,0);
        std::complex<double> avgsu4(0,0);
        double bn;
        for(int i=0;i<meuler+1;i++){
            double nn = ntr+i+1.0;
            dummy = pow(-1,nn);
            Gstemp = Gs<std::complex<double>>(std::complex<double>(xx,nn*hh),bp0,Qtt, Qttm);
            su[i+1] = su[i] + dummy * Gstemp[0];
            su1[i+1] = su1[i] + dummy * Gstemp[1];
            su2[i+1] = su2[i] + dummy * Gstemp[2];
            su3[i+1] = su3[i] + dummy * Gstemp[3];
            su4[i+1] = su4[i] + dummy * Gstemp[4];
            //Lo de abajo estaba originalmente en otro bucle
            bn = binomialCoefficients(meuler,(i+1)-1);
            //avgsu += bn*su[i];
            avgsu += bn*su[i+1];
            avgsu1 += bn*su1[i+1];
            avgsu2 += bn*su2[i+1];
            avgsu3 += bn*su3[i+1];
            avgsu4 += bn*su4[i+1];

        }
        double pp = pow(2.0,meuler);
        //std::complex<double> fun = uu*avgsu/(pp);
        std::vector<double> fun = std::vector<double>(5,0.0);
        fun[0] = real(uu*avgsu/(pp));
        fun[1] = real(uu*avgsu1/(pp));
        fun[2] = real(uu*avgsu2/(pp));
        fun[3] = real(uu*avgsu3/(pp));
        fun[4] = real(uu*avgsu4/(pp));

        return fun;
        //return 0.0;

    }

    //Get internal energy
    double InternalEnergy(double bp0){
        if (bp0 != bp)  SetPressure(bp0);
        double nsum=0.0;
        for (int i=0; i<nsize; i++){
            for(int j=0; j<nsize; j++){
                double dy = y(i)-y(j);
                nsum += sqrt(xvalues(i)*xvalues(j))*exp(-bp*sqrt(1.2*1.2-dy*dy));

                //nsum += sqrt(xvalues(i)*xvalues(j))*exp(-bp*b(i,j)); // esta linea solo es la buena
            }
        }
        return -1.0 + Aval*Aval*nsum/(bp);
    }

    //Analytical Gx
    double f(int n, double x, double b){
        if(x-b>0){
            return exp(-bp*x)*pow(x-b,n-1)/tgamma(n); //tgamma(n+1)=factorial(n)
        }
        else{
            return 0.0;
        }
    }

    double Gijx(double x, int i, int j){
 
        double factor = Aval*Aval*exp(beta)/(density*sqrt(xvalues[i]*xvalues[j]));
        double v = 1-exp(-beta);
        double t1 = f(1,x,a(i,j))-v*f(1,x,b(i,j));

        double t2 = 0.0;
        double t3 = 0.0;
        for (int k =0; k<nsize; k++){
            t2 += f(2,x,a(i,k)+a(k,j))-v*f(2,x,a(i,k)+b(k,j))-v*f(2,x,b(i,k)+a(k,j))+v*v*f(2,x,b(i,k)+b(k,j));
            //for (int l =0; l<nsize; l++){
             //   t3 += f(3,x,a(i,k)+a(k,l)+a(l,j)) -v*f(3,x,a(i,k)+a(k,l)+b(l,j))-v*f(3,x,a(i,k)+b(k,l)+a(l,j))-v*f(3,x,b(i,k)+a(k,l)+a(l,j)) +v*v*f(3,x,a(i,k)+b(k,l)+b(l,j))+v*v*f(3,x,b(i,k)+a(k,l)+b(l,j))+v*v*f(3,x,b(i,k)+b(k,l)+a(l,j))-v*v*v*f(3,x,b(i,k)+b(k,l)+b(l,j));
            //}
        }
        t2 *= Aval*Aval*exp(beta);
        //t3 *= pow(Aval,4)*exp(2*beta);

        return factor*(t1+t2);

    }

        // Compute  G(x)
    std::vector<double> Gx(double x, double bp0){
        if (bp0 != bp)  SetPressure(bp0);
        std::vector<double> gtotal = std::vector<double>(5,0.0);

        double gtemp = 0.0;
        double gmean = 0.0;
        for (int ix=0; ix<nsize; ix++){
            for(int jx=ix; jx<nsize; jx++){
                gtemp = Gijx(x,ix,jx);
                if (ix!=jx){
                    gmean += 2.0*xvalues(ix)*xvalues(jx)*gtemp;
                }
                else{
                    gmean += xvalues(ix)*xvalues(jx)*gtemp;
                } 
                //Get partials G(r)
                if(ix==0 && jx==0)          gtotal[1] = gtemp;
                else if(ix==0 && jx==(nsize+1)/2) gtotal[2] = gtemp;
                else if(ix==0 && jx==nsize-1)    gtotal[3] = gtemp;
                else if(ix==(nsize+1)/2 && jx==(nsize+1)/2) gtotal[4] = gtemp;
            }
        }
        gtotal[0] = gmean;
        return gtotal;
    }

   //private:

    Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> am;
    Eigen::Matrix<double,Eigen::Dynamic, Eigen::Dynamic> M;
    //Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> Qt;
    //Eigen::Matrix<std::complex<double>,Eigen::Dynamic, Eigen::Dynamic> Qtm;

    //Longitudinal separation at contact of two disks
    double a(int i, int j){
        double dy = y(i)-y(j);
        return sqrt(1.0-dy*dy);
    }

    double b(int i, int j){
        double dy = y(i)-y(j);
        
        
        return sqrt(r0*r0-dy*dy);
    }

    template<class T>
    T omega(T s, int i, int j){
        return exp(beta)*(exp(-a(i,j)*s)-(1.0-exp(-beta))*exp(-b(i,j)*s))/s;
    }

    template<class T>
    T omegad(T s, int i, int j){
        double aa = a(i,j);
        double bb = b(i,j);
        return exp(-(aa+bb)*s)*(-exp(beta+bb*s)*(1+aa*s)+exp(aa*s)*(-1.0+exp(beta))*(1+bb*s))*beta/(s*s);
    }

    double y(int i){
        return -eps/2.0 + i*eps/(nsize-1);
    }

    int UpdateMatrixA(){
        double temp;
        for (int i=0; i<nsize; i++){
            for(int j=0; j<nsize; j++){
                temp = (y(i)-y(j))*(y(i)-y(j));
                am(i,j) = sqrt(1.0-temp);
            }
        }
        return 0;
    }  

    
};

#endif
