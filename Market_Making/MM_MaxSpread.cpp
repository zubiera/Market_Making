//#include <iostream>
//#include <fstream>
//#include <cmath>
//#include <nag.h>
//#include <nag_stdlib.h>
//#include <nage04.h>
//#include <nagx02.h>
//#include <vector>
//#include <iostream>
//#include <string>
//#include <functional>
//#include <sstream>
//#include <algorithm>
//using namespace std;
//
///** \brief Function to wrap around NAG's e04ucc routine
// * \details Use a c++11 function wrapper to call NAG routine.
// * \param x On entry a vector containing the initial guess at the solution. On exit the location of the minimum point \f$\underline{x}^*\f$.
// * \param a A matrix (in row-column format) containing the linear set of constraints.
// * \param bl A vector containing the lower bound for {\bf all} constraints
// * \param bu A vector containing the upper bound for {\bf all} constraint
// * \param minValue On exit the value of the \f$F(\underline{x}^*)\f$ at the minimum point \f$\underline{x}^*\f$
// * \param F On entry, a function with the definition
// * ~~~~~~~~~~~~~~~{.c}
// * double funct(const double* xc){ ... }
// * ~~~~~~~~~~~~~~~
// * or a class with the member function
// * ~~~~~~~~~~~~~~~{.c}
// * double operator()(const double* xc) const { ... }
// * ~~~~~~~~~~~~~~~
// * \param tol Defines the accuracy of the function \f$F\f$.
// * \param MAXITER Maximium number of iterations to perform at each stage.
// * \return An integer defining the success of the routine. A non-zero number indicates a fail.
// */
//template<class T>
//int findMinNonLinOpt(vector<double> &x,vector<double> a,vector<double> bl,vector<double> bu,double &minValue,T &f,double tol=1.e-8,int MAXITER=100)
//{
//    // some NAG stuff
//    // comm.p is used to store a pointer to the object
//    Integer exit_status = 0, i, j, n, nclin, ncnlin, tda, totalvars;
//    Nag_Comm    comm;
//    // NAG's way of dealing with failure messages
//    NagError    fail;
//    Nag_E04_Opt options;
//    INIT_FAIL(fail);
//    std::vector<double> objgrd(x.size());
//
//    n=x.size();
//    nclin=(a.size()/n);
//    ncnlin=0;
//    totalvars = n + nclin + ncnlin;
//    tda=n;
//
//    if( bl.size()!=uint(totalvars) || bu.size()!=uint(totalvars) )
//    {
//        cout << " Inconsistent bounds on x ";
//        cout << " bl.size() = " << bl.size() << endl;
//        cout << " bu.size() = " << bu.size() << endl;
//        cout << " x.size() = " << x.size() << endl;
//        return 1;
//    }
//
//    if( a.size()!= uint(n*nclin) )
//    {
//        cout << " Wrong number of columns in matrix a ";
//        cout << " a.size() = " << a.size() << endl;
//        cout << " expected size = " << n*nclin << endl;
//        return 1;
//    }
//
//
//    int maxIter=100;
//    //int maxIter=100;
//    //  #pragma GCC diagnostic ignored "-Wunused-parameter"
//    // set up a wrapper to call the function
//    auto lambdaWrapper=[](Integer n, const double x[], double *objf,
//                          double objgrd[], Nag_Comm *comm)
//    {
//        /* objfun */
//        if (comm->flag == 0 || comm->flag == 2)
//        {
//            T* fx=reinterpret_cast<T*> (comm->p);
//            // compiler should check here that this function exists in the supplied object
//            *objf=fx->operator()(x);
//        }
//    };
//    void (*nagCallf)(Integer n, const double x[], double *objf,double objgrd[], Nag_Comm *comm)=lambdaWrapper;
//
//    // link your input function/object to the NAG comm
//    // note that f MUST be an object -- it cannot just be a standard function
//    comm.p=&f;
//
//    /* nag_opt_init (e04xxc).
//     * Initialization function for option setting
//     */
//    nag_opt_init(&options);
//    options.obj_deriv = Nag_FALSE;
//    options.con_deriv = Nag_FALSE;
//    options.max_iter = MAXITER;
//    options.minor_max_iter = MAXITER;
//    options.optim_tol = tol;
//    options.f_prec = tol;
//    // run the algorithm with no monitoring
//    nag_opt_nlp(n, nclin, ncnlin, a.data(), tda, bl.data(), bu.data(), nagCallf, NULLFN, x.data(), &minValue,objgrd.data(), &options, &comm, &fail);
//    // check for errors
//    if (fail.code != NE_NOERROR) {
//        printf("Error from nag_opt_nlp (e04ucc).\n%s\n", fail.message);
//        return 1;
//    }
//    // return 0 on success
//    return 0;
//}
//
//
////Including code for making M vector
//
//class MVector
//{
//    unsigned int N;
//    double* v;
//public:
//    explicit MVector():N(0),v(nullptr){}
//    explicit MVector(unsigned int n):N(n),v(new double[n]){}
//    explicit MVector(unsigned int n,double x):N(n),v(new double[n]){for(unsigned int i=0;i<N;i++)v[i]=x;}
//    ~MVector(){if(v!=nullptr) delete [] v;}
//    // normal copies
//    MVector(const MVector& X):N(X.N),v(new double[X.N]){std::copy(X.v,X.v+X.N,v);}
//    MVector& operator=(const MVector &X){
//        if(v!=nullptr && N!=X.N){delete [] v;}
//        if(N!=X.N){v = new double [X.N];}
//        N = X.N;
//        std::copy(X.v,X.v+X.N,v);
//        return *this;
//    }
//    // c++11 initialiser and move construct, move=
//    explicit MVector(std::initializer_list<double> list):N(list.size()),v(new double[list.size()]){auto  vPos=v;for(auto lPos=list.begin();lPos!=list.end();lPos++,vPos++)*vPos=*lPos;}
//    MVector(MVector&& X):N(X.N),v(X.v){X.N=0;X.v=nullptr;}
//    MVector& operator=(MVector&& X){
//        if(v!=nullptr){
//            if(v==X.v){return *this;}
//            delete [] v;
//        }
//        N=X.N;v=X.v;X.N=0;X.v=nullptr;
//        return *this;
//    }
//
//    double& operator[](int index){return v[index];}
//    double operator[](int index) const {return v[index];}
//    unsigned int size() const {return N;}
//    // resize vector
//    void resize(int n){
//        if(v!=nullptr) delete [] v;
//        v=new double[n];
//        N=n;
//    }
//    void resize(int n,double x){
//        resize(n);
//        double *vPos=v;
//        for(unsigned int i=0;i<N;i++)*vPos++=x;
//    }
//    // return array
//    double* returnArray(){return v;}
//    const double* returnArray() const {return v;}
//};
//
//std::ostream& operator<<(std::ostream &output,const MVector &X)
//{
//    //output << "( ";
//    unsigned int morethanone=std::min((unsigned int)(1),X.size());
//    for(unsigned int i=0;i<morethanone;i++)
//        output << X[i];
//    for(unsigned int i=morethanone;i<X.size();i++)
//        output << " , " << X[i];
//    //output << " )";
//    return output;
//}
//
//MVector operator*(const double& lhs,const MVector &rhs)
//{
//    MVector temp(rhs);
//    for(unsigned int i=0;i<rhs.size();i++)
//        temp[i]*=lhs;
//    return temp;
//}
//
//MVector operator*(const double& lhs,MVector&& rhs)
//{
//    MVector temp(std::move(rhs));
//    double *tempPos=temp.returnArray();
//    for(unsigned int i=0;i<temp.size();i++)
//        (*tempPos++)*=lhs;
//    return temp;
//}
//
//MVector operator*(const MVector& lhs,const double &rhs)
//{
//    MVector temp(lhs);
//    for(unsigned int i=0;i<lhs.size();i++)
//        temp[i]*=rhs;
//    return temp;
//}
//
//MVector operator/(const MVector& lhs,const double &rhs)
//{
//    MVector temp(lhs);
//    for(unsigned int i=0;i<lhs.size();i++)
//        temp[i]/=rhs;
//    return temp;
//}
//
//MVector operator+(const MVector& lhs,const MVector &rhs)
//{
//    MVector temp(lhs);
//    double *tempPos=temp.returnArray();
//    const double *rPos=rhs.returnArray();
//    for(unsigned int i=0;i<lhs.size();i++)
//        (*tempPos++)+=(*rPos++);
//    return temp;
//}
//
//MVector operator+(const MVector& lhs,MVector&& rhs)
//{
//    MVector temp(std::move(rhs));
//    double *tempPos=temp.returnArray();
//    const double *lPos=lhs.returnArray();
//    for(unsigned int i=0;i<lhs.size();i++)
//        (*tempPos++)+=(*lPos++);
//    return temp;
//}
//
//MVector operator-(const MVector& lhs,const MVector &rhs)
//{
//    MVector temp(lhs);
//    for(unsigned int i=0;i<lhs.size();i++)
//        temp[i]-=rhs[i];
//    return temp;
//}
//
//MVector operator-(const MVector& lhs,MVector&& rhs)
//{
//    MVector temp(std::move(rhs));
//    double *tempPos=temp.returnArray();
//    const double *lPos=lhs.returnArray();
//    for(unsigned int i=0;i<lhs.size();i++)
//    {
//        *tempPos*=-1;
//        (*tempPos++)+=(*lPos++);
//    }
//    return temp;
//}
//
//
//// n is number of steps
//// t_0 = a
//// t_n = b
//// t_i = a + ih where h=(b-a)/n
//// x_0 = alpha
//// this will solve to find x_n
//template<class SCALAR,class VECTOR,class FUNCTION>
//VECTOR RK4MethodTemplate(int n,SCALAR a,SCALAR b,
//                         VECTOR alpha,const FUNCTION &f)
//{
//    // local variables
//    SCALAR h,t;
//    VECTOR x;
//    // intialise values
//    h=(b-a)/(SCALAR)(n);
//    t=a;
//    x=alpha;
//    // implement Euler's method
//    for(int i=0;i<n;i++)
//    {
//        t = a + i*h; // update value of t to t_i
//        VECTOR k1 = h*f(x,t); // update x to x_{i+1}
//        VECTOR k2 = h*f(x+0.5*k1,t+0.5*h); // update x to x_{i+1}
//        VECTOR k3 = h*f(x+0.5*k2,t+0.5*h); // update x to x_{i+1}
//        VECTOR k4 = h*f(x+k3,t+h); // update x to x_{i+1}
//        x = x + 1./6.*(k1+k4+(2.*(k2+k3)));
//    }
//    return x; // returns x_n
//}
//
//
//
//
//int main()
//{
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//    // *********parameters from geuant 2017 Optimal Market Making Page 31 section 6******
//
//    int qMax     = 4;  ///divided 40M by delta = 10M so we have delta=1 and qmax=4
//    double sigma = 5.38e-6;     //2.15e-5;
//    double A     = 9.10e-4; //1.06e-3;
//    double Delta = 50e6;   //10e6;
//    double k     = 1.79e4/Delta; //5.47e3/Delta;
//    double gamma = 6e-5;
//    double xi    = gamma; // this is the sepcial case in the 2017 geuant paper
//    double b     = 0.*2.4248e-08*Delta*Delta; //2.4248e-08;     // spread for HY Index (no info in paper abuut penalty)
//    double T     = 7200; //corresponds to two hours
//    double max_spread = 10000.;
//
//
//    //Modelling paramters.
//    int     n = 100;  //number of time observations
//    double dT = T/n;  //time step
//
//
//    //store the value of the omega and delta at each time step and q value
//    vector<MVector> omega(n+1, MVector(2*qMax+1)), delta_b(n+1,MVector(2*qMax+1)), delta_a(n+1,MVector(2*qMax+1));
//
//    //no we initalise the solution at t=T for delta_a
//    for (int q=-qMax; q<=qMax; q++)
//    {
//        // from the terminal conditions of (3.9)
//        omega[n][qMax + q] =    -b*q*q;  //exp(-k*b*q*q); //penalty funcion is q
//        // from formula for delta_a assuming Big Delta=1
//        if(q>-qMax)
//            delta_a[n][qMax + q] = (1./gamma)*log(1 + gamma/k) + (omega[n][qMax +q]-omega[n][qMax +q-1]);
//        else
//            delta_a[n][qMax + q]=1.;
//    }
//
//    for (int q = qMax; q>=-qMax; q--)
//    {
//        // from the terminal conditions
//        omega[n][qMax+q] = -b*q*q; //exp(-k*b*q*q); //penalty funcion is q
//        // from formula for delta from (4.7) in the 2017 paper (assuming that Delta=1)
//        if(q<qMax)
//            delta_b[n][qMax + q] = (1./gamma)*log(1 + gamma/k) + (omega[n][qMax +q]-omega[n][qMax + q + 1]);
//        else
//            delta_b[n][qMax + q]=1.;
//
//    }
//
//    //Here we use a numerical integration scheme to find the value at T_i given T_{i+1}
//    for(int i=n-1; i>=0; i--)
//    {
//        // some of the constants in gueant 2017 eqn (3.9)
//        double alpha = 1./2.*gamma*sigma*sigma*Delta*Delta;
//        double beta  = A/gamma;
//        // constants in delta_a and delta_b
//        double c = 1./gamma*log(1 + gamma/k);
//
//        //Now we can solve gueant 2017 eqn (3.9) but calculating delta at each time step s to that (3.9) can be a function of delta.
//        // with initial condition  w(q,T_i) = omega(q,T_i)
//        // so that omega(q,T_{i-1}) = w(q,T_{i-1})
//
//        omega[i] = RK4MethodTemplate(100,i*dT,(i-1)*dT,omega[i+1],
//                                     [&]
//                                     (const  MVector &w,double t)
//                                     {
//                                         MVector F(2*qMax+1);
//
//                                         double db0 = c + w[0] - w[1];
//                                         double daQ = c + w[2*qMax] - w[2*qMax-1];
//
//                                         //This is the value when inventory is at -Q(no sell quotes)
//                                         F[0] = (alpha*qMax*qMax) - beta*exp(-k*db0)*(1-exp(-gamma*(db0 - (w[0]-w[1]))));
//
//                                         for(int q=-qMax+1;q<qMax;q++)
//                                         {
//                                             double delt_a, delt_b;
//                                             delt_a = c + w[qMax+q] - w[qMax+q-1];
//                                             delt_b = c + w[qMax+q] - w[qMax+q+1];
//                                             //delt_a = min(10000.,delt_a); // constrain
//                                             //delt_b = min(10000.,delt_b); // constrain
//
//
//                                             //The code below will find optimim choices of delta_a and delta_b within max spread
//                                             if(delt_a+delt_b < max_spread)
//                                                 //The normal calculation
//                                                 F[qMax + q] = (alpha*q*q) -beta*exp(-k*delt_b)*(1-exp(-gamma*(delt_b - (w[qMax+q]-w[qMax+q+1])))) - beta*exp(-k*delt_a)*(1-exp(-gamma*(delt_a - (w[qMax+q]-w[qMax+q-1]))));
//                                             else
//                                             {
//                                                 //define the function Func to minimise
//                                                 auto Func = [alpha,beta,gamma,q,qMax,k,w](const double *x){
//                                                     return -( (alpha*q*q) -beta*exp(-k*x[0])*(1-exp(-gamma*(x[0] - (w[qMax+q]-w[qMax+q+1])))) - beta*exp(-k*x[1])*(1-exp(-gamma*(x[1] - (w[qMax+q] -w[qMax+q-1])))) );
//                                                 };
//                                                 // a vector for the initial guess
//                                                 vector<double> x = {1.,0.};
//                                                 // store the min value
//                                                 double f;
//                                                 // run the minimisation
//                                                 findMinNonLinOpt( x , {1.,1.}, {-max_spread,-max_spread,-2.*max_spread} , {max_spread,max_spread,max_spread} , f, Func, 1.e-8,100);
//                                               F[qMax + q] = -f;
////                                               F[qMax + q] = (alpha*q*q) -beta*exp(-k*delt_b)*(1-exp(-gamma*(delt_b - (w[qMax+q]-w[qMax+q+1])))) - beta*exp(-k*delt_a)*(1-exp(-gamma*(delt_a - (w[qMax+q]-w[qMax+q-1]))));
//                                             //Need to add something here to store and overwrite the  delta_a and delta_b
//                                             }
//                                         }
//                                         //The value when inventory is at Q. (no buy quotes)
//                                         F[2*qMax] = (alpha*qMax*qMax) - beta*exp(-k*daQ)*(1-exp(-gamma*(daQ - (w[2*qMax]-w[2*qMax-1]))));
//                                         return F;
//                                     }
//                                     );
//
//
//        // Now using delta it is possible to caculate the value of otpimal ask proice
//        //because we have delta's formula
//        for(int q=0;q<=2*qMax-1;q++)
//            //delta_b[i][q] = min(10000.,omega[i][q] - omega[i][q+1] + 1/xi*log(1+xi/k));
//        delta_b[i][q] = omega[i][q] - omega[i][q+1] + 1/xi*log(1+xi/k);
//        delta_b[i][2*qMax] = 1.; //Is this needed???
//
//
//        for(int q=1;q<=2*qMax;q++)
//           // delta_a[i][q] = min(10000.,omega[i][q] - omega[i][q-1] + 1/xi*log(1+xi/k));
//        delta_a[i][q] = omega[i][q] - omega[i][q-1] + 1/xi*log(1+xi/k);
//        delta_a[i][0] = 1.; //Is this needed???
//
//
//
//    }
//
//
//
//
//    // Below is to plot the results for
//    const char *path_b="/Users/amitarfan1/Documents/Phd/3yr/Code/Output/delta_b_bound_10000.csv";
//    //and output results to file
//    ofstream output_b(path_b);
//    output_b << "Time_Step";
//    for(int ch = -4; ch <= 4; ++ch)
//    {
//        output_b << ",q="<<ch;
//    }
//    output_b << endl;
//    for( int i=0; i<=n; i++)
//    {
//        //The below is for frmatting the file for matplotlib
//        output_b << i*dT <<  " , " << delta_b[i]/Delta << endl;
//    }
//
//    const char *path_a="/Users/amitarfan1/Documents/Phd/3yr/Code/Output/delta_a_bound_10000.csv";
//    //and output results to file
//    ofstream output_a(path_a);
//    output_a << "Time_Step";
//    for(int ch = -4; ch <= 4; ++ch)
//    {
//        output_a << ",q="<<ch;
//    }
//    output_a << endl;
//    for( int i=0; i<=n; i++)
//    {
//
//        //The below is for frmatting the file for matplotlib
//        output_a << i*dT <<  " , " << delta_a[i]/Delta << endl;
//    }
//
//
//
//}

#include <iostream>
#include <fstream>
#include <cmath>
#include <nag.h>
#include <nag_stdlib.h>
#include <nage04.h>
#include <nagx02.h>
#include <vector>
#include <iostream>
#include <string>
#include <functional>
#include <sstream>
#include <algorithm>

using namespace std;

/** \brief Function to wrap around NAG's e04ucc routine
 * \details Use a c++11 function wrapper to call NAG routine.
 * \param x On entry a vector containing the initial guess at the solution. On exit the location of the minimum point \f$\underline{x}^*\f$.
 * \param a A matrix (in row-column format) containing the linear set of constraints.
 * \param bl A vector containing the lower bound for {\bf all} constraints
 * \param bu A vector containing the upper bound for {\bf all} constraint
 * \param minValue On exit the value of the \f$F(\underline{x}^*)\f$ at the minimum point \f$\underline{x}^*\f$
 * \param F On entry, a function with the definition
 * ~~~~~~~~~~~~~~~{.c}
 * double funct(const double* xc){ ... }
 * ~~~~~~~~~~~~~~~
 * or a class with the member function
 * ~~~~~~~~~~~~~~~{.c}
 * double operator()(const double* xc) const { ... }
 * ~~~~~~~~~~~~~~~ my_funcs
 * \param tol Defines the accuracy of the function \f$F\f$.
 * \param MAXITER Maximium number of iterations to perform at each stage.
 * \return An integer defining the success of the routine. A non-zero number indicates a fail.
 */
template<class T>
int findMinNonLinOpt(vector<double> &x,vector<double> a,vector<double> bl,vector<double> bu,double &minValue,T &f,double tol,int MAXITER)
{
    // some NAG stuff
    // comm.p is used to store a pointer to the object
    Integer exit_status = 0, i, j, n, nclin, ncnlin, tda, totalvars;
    Nag_Comm    comm;
    // NAG's way of dealing with failure messages
    NagError    fail;
    Nag_E04_Opt options;
    INIT_FAIL(fail);
    std::vector<double> objgrd(x.size());
    
    n=x.size();
    nclin=(a.size()/n);
    ncnlin=0;
    totalvars = n + nclin + ncnlin;
    tda=n;
    
    if( bl.size()!=uint(totalvars) || bu.size()!=uint(totalvars) )
    {
        cout << " Inconsistent bounds on x ";
        cout << " bl.size() = " << bl.size() << endl;
        cout << " bu.size() = " << bu.size() << endl;
        cout << " x.size() = " << x.size() << endl;
        return 1;
    }
    
    if( a.size()!= uint(n*nclin) )
    {
        cout << " Wrong number of columns in matrix a ";
        cout << " a.size() = " << a.size() << endl;
        cout << " expected size = " << n*nclin << endl;
        return 1;
    }
    
    
    int maxIter=100;
    //int maxIter=100;
    //  #pragma GCC diagnostic ignored "-Wunused-parameter"
    // set up a wrapper to call the function
    auto lambdaWrapper=[](Integer n, const double x[], double *objf,
                          double objgrd[], Nag_Comm *comm)
    {
        /* objfun */
        if (comm->flag == 0 || comm->flag == 2)
        {
            T* fx=reinterpret_cast<T*> (comm->p);
            // compiler should check here that this function exists in the supplied object
            *objf=fx->operator()(x);
        }
    };
    void (*nagCallf)(Integer n, const double x[], double *objf,double objgrd[], Nag_Comm *comm)=lambdaWrapper;
    
    // link your input function/object to the NAG comm
    // note that f MUST be an object -- it cannot just be a standard function
    comm.p=&f;
    
    /* nag_opt_init (e04xxc).
     * Initialization function for option setting
     */
    nag_opt_init(&options);
    options.obj_deriv = Nag_FALSE;
    options.con_deriv = Nag_FALSE;
    options.max_iter = MAXITER;
    options.minor_max_iter = MAXITER;
    options.optim_tol = tol;
    options.f_prec = tol;
    options.list=Nag_FALSE;
    options.print_level = Nag_NoPrint;
    options.print_deriv = Nag_D_NoPrint;
    // run the algorithm with no monitoring
    nag_opt_nlp(n, nclin, ncnlin, a.data(), tda, bl.data(), bu.data(), nagCallf, NULLFN, x.data(), &minValue,objgrd.data(), &options, &comm, &fail);
    // check for errors
    if (fail.code != NE_NOERROR) {
        printf("Error from nag_opt_nlp (e04ucc).\n%s\n", fail.message);
        return 1;
    }
    // return 0 on success
    return 0;
}


/** \brief Function to wrap around NAG's e04ucc routine
 * \details Use a c++11 function wrapper to call NAG routine.
 * \param x On entry a vector containing the initial guess at the solution. On exit the location of the minimum point \f$\underline{x}^*\f$.
 * \param a A matrix (in row-column format) containing the linear set of constraints.
 * \param bl A vector containing the lower bound for {\bf all} constraints
 * \param bu A vector containing the upper bound for {\bf all} constraint
 * \param minValue On exit the value of the \f$F(\underline{x}^*)\f$ at the minimum point \f$\underline{x}^*\f$
 * \param F On entry, a function with the definition
 * ~~~~~~~~~~~~~~~{.c}
 * double funct(const double* xc){ ... }
 * ~~~~~~~~~~~~~~~
 * or a class with the member function
 * ~~~~~~~~~~~~~~~{.c}
 * double operator()(const double* xc) const { ... }
 * ~~~~~~~~~~~~~~~
 * \param tol Defines the accuracy of the function \f$F\f$.
 * \param MAXITER Maximium number of iterations to perform at each stage.
 * \return An integer defining the success of the routine. A non-zero number indicates a fail.
 */
template<class FUNCTION,class JACOBIAN>
int findMinNonLinOpt(vector<double> &x,vector<double> a,vector<double> bl,vector<double> bu,double &minValue,FUNCTION &f,JACOBIAN &J,double tol,int MAXITER)
{
    // some NAG stuff
    // comm.p is used to store a pointer to the object
    Integer exit_status = 0, i, j, n, nclin, ncnlin, tda, totalvars;
    Nag_Comm    comm;
    // NAG's way of dealing with failure messages
    NagError    fail;
    Nag_E04_Opt options;
    INIT_FAIL(fail);
    std::vector<double> objgrd(x.size());
    
    n=x.size();
    nclin=(a.size()/n);
    ncnlin=0;
    totalvars = n + nclin + ncnlin;
    tda=n;
    
    if( bl.size()!=uint(totalvars) || bu.size()!=uint(totalvars) )
    {
        cout << " Inconsistent bounds on x ";
        cout << " bl.size() = " << bl.size() << endl;
        cout << " bu.size() = " << bu.size() << endl;
        cout << " x.size() = " << x.size() << endl;
        return 1;
    }
    
    if( a.size()!= uint(n*nclin) )
    {
        cout << " Wrong number of columns in matrix a ";
        cout << " a.size() = " << a.size() << endl;
        cout << " expected size = " << n*nclin << endl;
        return 1;
    }
    
    
    int maxIter=100;
    //int maxIter=100;
    //  #pragma GCC diagnostic ignored "-Wunused-parameter"
    // set up a wrapper to call the function
    auto lambdaWrapper=[](Integer n, const double x[], double *objf,
                          double objgrd[], Nag_Comm *comm)
    {
        vector<void*>* my_funcs = reinterpret_cast< vector<void*>* > (comm->p);
        /* objfun */
        if (comm->flag == 0 || comm->flag == 2)
        {
            FUNCTION* fx=reinterpret_cast<FUNCTION*> ((*my_funcs)[0]);
            // compiler should check here that this function exists in the supplied object
            *objf=fx->operator()(x);
        }
        if ( comm->flag == 2)
        {
            JACOBIAN* gx=reinterpret_cast<JACOBIAN*> ((*my_funcs)[1]);
            // compiler should check here that this function exists in the supplied object
            gx->operator()(x,objgrd);
        }
    };
    void (*nagCallf)(Integer n, const double x[], double *objf,double objgrd[], Nag_Comm *comm)=lambdaWrapper;
    
    // link your input function/object to the NAG comm
    // note that f MUST be an object -- it cannot just be a standard function
    vector<void*> my_funcs = {&f,&J};
    comm.p=&my_funcs;
    
    /* nag_opt_init (e04xxc).
     * Initialization function for option setting
     */
    nag_opt_init(&options);
    //     options.obj_deriv = Nag_FALSE;
    options.step_limit = 100;
    options.con_deriv = Nag_FALSE;
    options.max_iter = MAXITER;
    options.minor_max_iter = MAXITER;
    options.optim_tol = tol;
    options.f_prec = tol;
    options.list=Nag_FALSE;
    options.print_level = Nag_NoPrint;
    options.print_deriv = Nag_D_NoPrint;
    // run the algorithm with no monitoring
    nag_opt_nlp(n, nclin, ncnlin, a.data(), tda, bl.data(), bu.data(), nagCallf, NULLFN, x.data(), &minValue,objgrd.data(), &options, &comm, &fail);
    // check for errors
    if (fail.code != NE_NOERROR) {
        printf("Error from nag_opt_nlp (e04ucc).\n%s\n", fail.message);
        return 1;
    }
    // return 0 on success
    return 0;
}


//Including code for making M vector

class MVector
{
    unsigned int N;
    double* v;
public:
    explicit MVector():N(0),v(nullptr){}
    explicit MVector(unsigned int n):N(n),v(new double[n]){}
    explicit MVector(unsigned int n,double x):N(n),v(new double[n]){for(unsigned int i=0;i<N;i++)v[i]=x;}
    ~MVector(){if(v!=nullptr) delete [] v;}
    // normal copies
    MVector(const MVector& X):N(X.N),v(new double[X.N]){std::copy(X.v,X.v+X.N,v);}
    MVector& operator=(const MVector &X){
        if(v!=nullptr && N!=X.N){delete [] v;}
        if(N!=X.N){v = new double [X.N];}
        N = X.N;
        std::copy(X.v,X.v+X.N,v);
        return *this;
    }
    // c++11 initialiser and move construct, move=
    explicit MVector(std::initializer_list<double> list):N(list.size()),v(new double[list.size()]){auto  vPos=v;for(auto lPos=list.begin();lPos!=list.end();lPos++,vPos++)*vPos=*lPos;}
    MVector(MVector&& X):N(X.N),v(X.v){X.N=0;X.v=nullptr;}
    MVector& operator=(MVector&& X){
        if(v!=nullptr){
            if(v==X.v){return *this;}
            delete [] v;
        }
        N=X.N;v=X.v;X.N=0;X.v=nullptr;
        return *this;
    }
    
    double& operator[](int index){return v[index];}
    double operator[](int index) const {return v[index];}
    unsigned int size() const {return N;}
    // resize vector
    void resize(int n){
        if(v!=nullptr) delete [] v;
        v=new double[n];
        N=n;
    }
    void resize(int n,double x){
        resize(n);
        double *vPos=v;
        for(unsigned int i=0;i<N;i++)*vPos++=x;
    }
    // return array
    double* returnArray(){return v;}
    const double* returnArray() const {return v;}
};

std::ostream& operator<<(std::ostream &output,const MVector &X)
{
    //output << "( ";
    unsigned int morethanone=std::min((unsigned int)(1),X.size());
    for(unsigned int i=0;i<morethanone;i++)
        output << X[i];
    for(unsigned int i=morethanone;i<X.size();i++)
        output << " , " << X[i];
    //output << " )";
    return output;
}

MVector operator*(const double& lhs,const MVector &rhs)
{
    MVector temp(rhs);
    for(unsigned int i=0;i<rhs.size();i++)
        temp[i]*=lhs;
    return temp;
}

MVector operator*(const double& lhs,MVector&& rhs)
{
    MVector temp(std::move(rhs));
    double *tempPos=temp.returnArray();
    for(unsigned int i=0;i<temp.size();i++)
        (*tempPos++)*=lhs;
    return temp;
}

MVector operator*(const MVector& lhs,const double &rhs)
{
    MVector temp(lhs);
    for(unsigned int i=0;i<lhs.size();i++)
        temp[i]*=rhs;
    return temp;
}

MVector operator/(const MVector& lhs,const double &rhs)
{
    MVector temp(lhs);
    for(unsigned int i=0;i<lhs.size();i++)
        temp[i]/=rhs;
    return temp;
}

MVector operator+(const MVector& lhs,const MVector &rhs)
{
    MVector temp(lhs);
    double *tempPos=temp.returnArray();
    const double *rPos=rhs.returnArray();
    for(unsigned int i=0;i<lhs.size();i++)
        (*tempPos++)+=(*rPos++);
    return temp;
}

MVector operator+(const MVector& lhs,MVector&& rhs)
{
    MVector temp(std::move(rhs));
    double *tempPos=temp.returnArray();
    const double *lPos=lhs.returnArray();
    for(unsigned int i=0;i<lhs.size();i++)
        (*tempPos++)+=(*lPos++);
    return temp;
}

MVector operator-(const MVector& lhs,const MVector &rhs)
{
    MVector temp(lhs);
    for(unsigned int i=0;i<lhs.size();i++)
        temp[i]-=rhs[i];
    return temp;
}

MVector operator-(const MVector& lhs,MVector&& rhs)
{
    MVector temp(std::move(rhs));
    double *tempPos=temp.returnArray();
    const double *lPos=lhs.returnArray();
    for(unsigned int i=0;i<lhs.size();i++)
    {
        *tempPos*=-1;
        (*tempPos++)+=(*lPos++);
    }
    return temp;
}


// n is number of steps
// t_0 = a
// t_n = b
// t_i = a + ih where h=(b-a)/n
// x_0 = alpha
// this will solve to find x_n
template<class SCALAR,class VECTOR,class FUNCTION>
VECTOR RK4MethodTemplate(int n,SCALAR a,SCALAR b,
                         VECTOR alpha,const FUNCTION &f)
{
    // local variables
    SCALAR h,t;
    VECTOR x;
    // intialise values
    h=(b-a)/(SCALAR)(n);
    t=a;
    x=alpha;
    // implement Euler's method
    for(int i=0;i<n;i++)
    {
        t = a + i*h; // update value of t to t_i
        VECTOR k1 = h*f(x,t); // update x to x_{i+1}
        VECTOR k2 = h*f(x+0.5*k1,t+0.5*h); // update x to x_{i+1}
        VECTOR k3 = h*f(x+0.5*k2,t+0.5*h); // update x to x_{i+1}
        VECTOR k4 = h*f(x+k3,t+h); // update x to x_{i+1}
        x = x + 1./6.*(k1+k4+(2.*(k2+k3)));
    }
    return x; // returns x_n
}




int main(int argc, char *argv[])
{
    
    // *********parameters from geuant 2017 Optimal Market Making Page 31 section 6******
    
    int qMax     = 4;  ///divided 40M by delta = 10M so we have delta=1 and qmax=4
    double sigma = 5.38e-6;     //2.15e-5;
    double A     = 9.10e-4; //1.06e-3;
    double Delta = 50e6;   //10e6;
    double k     = 1.79e4/Delta; //5.47e3/Delta;
    double gamma = 6e-5;
    double xi    = gamma; // this is the sepcial case in the 2017 geuant paper
    double b     = 0.;     // spread for HY Index (no info in paper abuut penalty)
    double T     = 7200.; //corresponds to two hours
    double max_spread = std::stod(argv[1]);
    
    
    //Modelling paramters.
    int     n = 100;  //number of time observations
    double dT = T/n;  //time step
    
    
    //store the value of the omega and delta at each time step and q value
    vector<MVector> omega(n+1, MVector(2*qMax+1)), delta_b(n+1,MVector(2*qMax+1)), delta_a(n+1,MVector(2*qMax+1));
    
    //no we initalise the solution at t=T for delta_a
    for (int q=-qMax; q<=qMax; q++)
    {
        // from the terminal conditions of (3.9)
        omega[n][qMax + q] =    -b*q*q;  //exp(-k*b*q*q); //penalty funcion is q
        // from formula for delta_a assuming Big Delta=1
        if(q>-qMax)
            delta_a[n][qMax + q] = (1./gamma)*log(1 + gamma/k) + (omega[n][qMax +q]-omega[n][qMax +q-1]);
        else
            delta_a[n][qMax + q]=1.;
    }
    
    for (int q = qMax; q>=-qMax; q--)
    {
        // from the terminal conditions
        omega[n][qMax+q] = -b*q*q; //exp(-k*b*q*q); //penalty funcion is q
        // from formula for delta from (4.7) in the 2017 paper (assuming that Delta=1)
        if(q<qMax)
            delta_b[n][qMax + q] = (1./gamma)*log(1 + gamma/k) + (omega[n][qMax +q]-omega[n][qMax + q + 1]);
        else
            delta_b[n][qMax + q]=1.;
        
    }
    
    //Here we use a numerical integration scheme to find the value at T_i given T_{i+1}
    for(int i=n-1; i>=0; i--)
    {
        // some of the constants in gueant 2017 eqn (3.9)
        double alpha = 1./2.*gamma*sigma*sigma*Delta*Delta;
        double beta  = A/gamma;
        // constants in delta_a and delta_b
        double c = 1./gamma*log(1 + gamma/k);
        
        vector<double> bl ={-max_spread,-max_spread,0.} , bu ={max_spread,max_spread,max_spread};
        //Now we can solve gueant 2017 eqn (3.9) but calculating delta at each time step s to that (3.9) can be a function of delta.
        // with initial condition  w(q,T_i) = omega(q,T_i)
        // so that omega(q,T_{i-1}) = w(q,T_{i-1})
        
        omega[i] = RK4MethodTemplate(100,i*dT,(i-1)*dT,omega[i+1],
                                     [&]
                                     (const  MVector &w,double t)
                                     {
                                         MVector F(2*qMax+1);
                                         
                                         double db0 = c + w[0] - w[1];
                                         delta_b[i][0] = db0;
                                         F[0] = (alpha*qMax*qMax) - beta*exp(-k*db0)*(1-exp(-gamma*(db0 - (w[0]-w[1]))));
                                         for(int q=-qMax+1;q<qMax;q++)
                                         {
                                             double delt_a, delt_b;
                                             // JUST TO SHOW THIS WORKS, PUT WRONG GUESS IN HERE
                                             delt_a = c;// + w[qMax+q] - w[qMax+q-1];
                                             delt_b = c;// + w[qMax+q] - w[qMax+q+1];
                                             //define the function Func to minimise
                                             auto Func = [alpha,beta,gamma,q,qMax,k,w](const double *x){
                                                 return ( (alpha*q*q)
                                                         - beta*exp(-k*x[0])*(1-exp(-gamma*(x[0] - (w[qMax+q]-w[qMax+q+1]))))
                                                         - beta*exp(-k*x[1])*(1-exp(-gamma*(x[1] - (w[qMax+q]-w[qMax+q-1])))) );
                                             };
                                             auto Jac = [beta,gamma,q,qMax,k,w](const double *x,double *J){
                                                 J[0] =  k*beta*exp(-k*x[0])*(1-exp(-gamma*(x[0] - (w[qMax+q]-w[qMax+q+1]))))-gamma*beta*exp(-k*x[0])*exp(-gamma*(x[0] - (w[qMax+q]-w[qMax+q+1])));
                                                 J[1] =  k*beta*exp(-k*x[1])*(1-exp(-gamma*(x[1] - (w[qMax+q]-w[qMax+q-1])))) -gamma*beta*exp(-k*x[1])*exp(-gamma*(x[1] - (w[qMax+q]-w[qMax+q-1])));
                                             };
                                             std::vector<double> x = { delt_b , delt_a};
                                             // store the min value
                                             double f;
                                             
                                             
                                             // The initial Guess
                                             x = { 1000. , 1000.};
                                             
                                             
                                             // run the minimisation
                                             findMinNonLinOpt( x , {1.,1.}, bl , bu , f, Func,Jac, 1.e-6,100);
                                             delta_b[i][qMax+q] = x[0];
                                             delta_a[i][qMax+q] = x[1];
                                             F[qMax + q] = f;
                                         }
                                         double daQ = c + w[2*qMax] - w[2*qMax-1];
                                         delta_a[i][2*qMax] = daQ;
                                         F[2*qMax] = (alpha*qMax*qMax) - beta*exp(-k*daQ)*(1-exp(-gamma*(daQ - (w[2*qMax]-w[2*qMax-1]))));
                                         return F;
                                     }
                                     );
        
    }
    
    
    
    // Below is to plot the results for
        std::string path_b = "/Users/amitarfan1/Documents/Phd/3yr/Code/Output/delta_b_bound_"+std::string(argv[1])+".csv";
    //const char *path_b="/Users/amitarfan1/Documents/Phd/3yr/Code/Output/delta_b_bound_15000.csv";
    //and output results to file
    ofstream output_b(path_b);
    output_b << "Time_Step";
    for(int ch = -4; ch <= 4; ++ch)
    {
        output_b << ",q="<<ch;
    }
    output_b << endl;
    for( int i=0; i<=n; i++)
    {
        //The below is for frmatting the file for matplotlib
        output_b << i*dT <<  " , " << delta_b[i]/Delta <<  endl;
        //output_b << i*dT <<  " , " << omega[i] << endl;
    }
    std::string path_a = "/Users/amitarfan1/Documents/Phd/3yr/Code/Output/delta_a_bound_"+std::string(argv[1])+".csv";
    //std::cout << path_a << std::endl;
    
    //and output results to file
    ofstream output_a(path_a);
    output_a << "Time_Step";
    for(int ch = -4; ch <= 4; ++ch)
    {
        output_a << ",q="<<ch;
    }
    output_a << endl;
    for( int i=0; i<=n; i++)
    {
        //The below is for frmatting the file for matplotlib
        output_a << i*dT <<  " , " << delta_a[i]/Delta << endl;
    }
    
    
    
}
