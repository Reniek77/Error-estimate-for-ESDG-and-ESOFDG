#include<bits/stdc++.h>

using namespace std;

/************************************************************************** 
PART 1: Basic Definition
***************************************************************************/

const int MAXX=20202;
const int Max_Round=6;
const double CFL=0.1;
const double PI=acos(-1.0);
const double X_LEFT  = 0.0;
const double X_RIGHT = 2.0*PI;
const double DOMAIN_LENGTH = X_RIGHT - X_LEFT;

double theta_in_LF=0.5;

using State = std::vector<double>[MAXX];

int k,N,hat_type;
int OF_YES;
double h,tau,now_time,T_end,real_time,sigma0;
std::vector<double> mesh_point, cell_h;

class Matrix{
    private:
        double matrix[7][7];
    public:
        Matrix(){}
        void init(){
            for(int i=0;i<7;++i)
            for(int j=0;j<7;++j)
            matrix[i][j]=0.0;
        }
        double& operator()(int i,int j){return matrix[i][j];}
        const double& operator()(int i,int j)const{return matrix[i][j];}
};

class Numerical_Integral_and_Differential{
    private:
        int Number_of_Points;
        vector<pair<double,double> >Gauss_Labatto_Points;
        Matrix Difference_Matrix;
        Matrix Lagrange_Derivate;
    public:
        Numerical_Integral_and_Differential(){}

        int Get_Number_of_Points() const {return Number_of_Points;}

        void Init(int Degree_of_Poly){
            Number_of_Points=Degree_of_Poly+1;
            Difference_Matrix.init();
            Gauss_Labatto_Points.resize(Number_of_Points+3,std::make_pair(0.0,0.0));
            Lagrange_Derivate.init();
        }

        void Calculate_Gauss_Labatto_Points(){
            int Poin=Number_of_Points;
            if(Poin==2){
                Gauss_Labatto_Points[0]={-1.00000,1.0000000};
                Gauss_Labatto_Points[1]={1.000000000,1.0000000};
            }
            if(Poin==3){
                Gauss_Labatto_Points[0]={-1.00000,0.33333333333333333333333333};
                Gauss_Labatto_Points[1]={0.000000,1.3333333333333333333333333};
                Gauss_Labatto_Points[2]={1.0000000,0.33333333333333333333333};
            }
            if(Poin==4){
                Gauss_Labatto_Points[0]={-1.0000000,1.0/6.0};
                Gauss_Labatto_Points[1]={-0.447213595499958,5.0/6.0};
                Gauss_Labatto_Points[2]={0.447213595499958,5.0/6.0};
                Gauss_Labatto_Points[3]={1.0000000,1.0/6.0};
            }
            if(Poin==5){
                Gauss_Labatto_Points[0]={-1.0000000,0.1};
                Gauss_Labatto_Points[1]={-0.654653670707977,0.544444444444444};
                Gauss_Labatto_Points[2]={0.0,0.711111111111111};
                Gauss_Labatto_Points[3]={0.654653670707977,0.544444444444444};
                Gauss_Labatto_Points[4]={1.0000000,0.1};
            }
            return;
        }

        void Calculate_Lagrange_Derivative(){
            for(register int i=0;i<=k;++i)
            {
                for(register int j=0;j<=k;++j)
                if(i==j)Lagrange_Derivate(i,j)=1.0;
                else{
                    Lagrange_Derivate(i,j)=1.0;
                    for(register int l=0;l<=k;++l)
                    if(l==j||l==i)continue;
                    else Lagrange_Derivate(i,j)=(Lagrange_Derivate(i,j)*(Gauss_Labatto_Points[j].first-Gauss_Labatto_Points[l].first)/(Gauss_Labatto_Points[i].first-Gauss_Labatto_Points[l].first));
                }
            }
            return;
        }

        void Get_Difference_Matrix(){
            for(register int i=0;i<=k;++i)// i-th point
            for(register int j=0;j<=k;++j)// j-th 
            {
                double sum=0.0;
                if(i==j)
                {
                    for(register int l=0;l<=k;++l)
                    if(l==j)continue;
                    else sum+=(1.0/(Gauss_Labatto_Points[i].first-Gauss_Labatto_Points[l].first));
                    Difference_Matrix(i,j)=sum;
                }
                else
                {
                    sum=1.0;
                    double fac=1.0/(Gauss_Labatto_Points[j].first-Gauss_Labatto_Points[i].first);
                    for(register int l=0;l<=k;++l)
                    if(l==j||l==i)continue;
                    else sum*=(Gauss_Labatto_Points[i].first-Gauss_Labatto_Points[l].first)/(Gauss_Labatto_Points[j].first-Gauss_Labatto_Points[l].first);

                    Difference_Matrix(i,j)=fac*sum;
                }
            }
            return;
        }

        vector<pair<double,double> > get_GL(){return Gauss_Labatto_Points;}
        const vector<pair<double,double> >& get_GL()const{return Gauss_Labatto_Points;}
        Matrix get_DF(){return Difference_Matrix;}
        const Matrix get_DF()const{return Difference_Matrix;}
        Matrix get_LD(){return Lagrange_Derivate;}
        const Matrix get_LD()const{return Lagrange_Derivate;}
}NID;

std::vector<double>uh[MAXX],u1[MAXX],u2[MAXX],u3[MAXX];

double norm_h[Max_Round<<1],trun[Max_Round<<1],norm_L2[Max_Round<<1],norm_Linf[Max_Round<<1];

/************************************************************************** 
PART 1.5: Mesh Generation
***************************************************************************/

inline void Generate_Mesh(){

    h = DOMAIN_LENGTH / N;
    mesh_point.resize(N + 1);
    cell_h.resize(N + 2);

    static std::mt19937 randd(20020702u);
    static std::uniform_real_distribution <double> dist(-0.2, 0.2);

    mesh_point[0] = X_LEFT;

    for(int i=1; i<=N; i+=2){
        double delta = dist(randd);
        cell_h[i] = h * (1.0 + delta);
        cell_h[i+1] = h * (1.0 - delta);
        
        mesh_point[i] = mesh_point[i-1] + cell_h[i];
        mesh_point[i+1] = mesh_point[i] + cell_h[i+1];
    }

    mesh_point[N] = X_RIGHT;
    return;
}

inline double Min_Cell_h(){

    double res=1e100;
    for(register int i=1;i<=N;++i)res=min(res,cell_h[i]);
    return res;
}

/************************************************************************** 
PART 2: Real Information
***************************************************************************/

double entropy_var(double u){
    // return 0.9*u + 0.1*exp(u);
    return 2.0*u + exp(u);
}

double entropy_var_prime(double u){
    // return 0.9 + 0.1*exp(u);
    return 2.0 + exp(u);
}

double f(double u){
    return u*u*0.5;
    // return u;
}

double f_prime(double u){
    return u;
    // return 1.0;
}

double f_second(double u){
    return 1.0;
    // return 0.0;
}
double entropy_potential(double u){
    // return 0.15*u*u*u+ (0.05*u*u - 0.1*u + 0.1)*exp(u);
    return u*u*u/3.0+ (0.5*u*u - u + 1.0)*exp(u);
}

double f_s(double uL, double uR){
    if(fabs(uR-uL) < 1.0e-12){
        return f(0.5*(uL+uR));
    }

    return (entropy_potential(uR)-entropy_potential(uL))
         / (entropy_var(uR)-entropy_var(uL));
}

double f_hat_lf(double uL, double uR){
    // const double alpha = 5.0; 
    const double alpha=max(fabs(f_prime(uL)),fabs(f_prime(uR)));
    return 0.5 * (f(uL) + f(uR)) - theta_in_LF * alpha * (uR - uL);
}

double f_hat(double uL, double uR){
    if(hat_type==1)return f_hat_lf(uL,uR);
    if(hat_type==2)return f_hat_lf(uL,uR);
    return 0.0;
}

double init_solution(double x){
    return exp(cos(x))*sin(x)+sin(x)*sin(x);
}

double init_solution_prime(double x){
    return (-sin(x)*sin(x)*exp(cos(x))+cos(x)*exp(cos(x)))+2.0*cos(x)*sin(x);
}

double true_solution(double x, double t){
    // return 0.25+sin(x);
    // Periodic reduction to [X_LEFT, X_RIGHT).
    x = X_LEFT + fmod(x - X_LEFT, DOMAIN_LENGTH);
    if(x < X_LEFT) x += DOMAIN_LENGTH;
    if(x >= X_RIGHT) x -= DOMAIN_LENGTH;

    // Characteristic formula:
    //
    //      u(x,t)=u0(xi),
    //      x = xi + t*f'(u0(xi))    modulo DOMAIN_LENGTH.
    //
    // Before shock formation, the lifted map
    //      Phi(xi)=xi+t*f'(u0(xi))
    // satisfies Phi(xi+L)=Phi(xi)+L.

    double phi0 = X_LEFT + t * f_prime(init_solution(X_LEFT));

    double target = x + ceil((phi0 - x) / DOMAIN_LENGTH) * DOMAIN_LENGTH;

    double left = X_LEFT, right = X_RIGHT;

    for(int it = 0; it < 100; ++it){
        double mid = 0.5 * (left + right);
        double phi_mid = mid + t * f_prime(init_solution(mid));

        if(phi_mid < target) left = mid;
        else right = mid;
    }

    double xi = 0.5 * (left + right);
    return init_solution(xi);
}

inline double estimate_shock_time(){

    double min_val = 1e100;
    const int M = 20000;

    for(int m = 0; m < M; ++m){
        double x = X_LEFT + (m + 0.5) * DOMAIN_LENGTH / M;

        double u0  = init_solution(x);
        double u0x = init_solution_prime(x);

        double val = f_second(u0) * u0x;

        min_val = min(min_val, val);
    }

    if(min_val >= 0.0) return 1e100;

    return -1.0 / min_val;
}

double estimate_max_speed(){

    double alpha = 0.0;

    for(int i=1;i<=N;++i)
    {
        for(int j=0;j<=k;++j)
        alpha = max(alpha, fabs(f_prime(uh[i][j])));
    }

    return alpha;
}

/*********************************************************************************
Part 3: Calculate
*********************************************************************************/

inline double Cell_Average(int i, const State& u){

    double numerator = 0.0;
    double denominator = 0.0;

    auto gl = NID.get_GL();

    for(int l=0;l<=k;++l){
        double weight = gl[l].second;
        numerator += weight*u[i][l];
        denominator += weight;
    }

    return numerator/denominator;
}

inline double Derivative_At_Endpoint(
    int i,
    int endpoint,
    const State& u
){
    double derivative = 0.0;
    Matrix D = NID.get_DF();

    for(int l=0;l<=k;++l){
        derivative += D(endpoint,l)*u[i][l];
    }

    return 2.0*derivative/cell_h[i];
}

inline double Damping_Coefficient(int i, const State& u){

    int left_cell  = (i==1 ? N : i-1);
    int right_cell = (i==N ? 1 : i+1);

    // Jumps of u_h.
    double jump_u_left =
        u[i][0] - u[left_cell][k];

    double jump_u_right =
        u[right_cell][0] - u[i][k];

    // Jumps of the first derivative.
    double jump_ux_left =
        Derivative_At_Endpoint(i,0,u)
        - Derivative_At_Endpoint(left_cell,k,u);

    double jump_ux_right =
        Derivative_At_Endpoint(right_cell,0,u)
        - Derivative_At_Endpoint(i,k,u);

    double sigma_squared =
        0.5*(
            jump_u_left*jump_u_left
            + jump_u_right*jump_u_right
        )
        + 0.25*cell_h[i]*cell_h[i]*(
            jump_ux_left*jump_ux_left
            + jump_ux_right*jump_ux_right
        );
    
    double sigma_now=sqrt(sigma_squared);
    sigma0=max(sigma0,sigma_now);
    return sigma_now;
}

inline double compute_tau(){
    double hmin = Min_Cell_h();
    double max_speed = estimate_max_speed();

    double CFL_nw=CFL / ((max(1.0, max_speed)*2.0*theta_in_LF)+sigma0);

    double dt=0.0;
    if(k <= 2){
        dt = CFL_nw * hmin;
    }else{
        dt = CFL_nw * pow(hmin, (k + 1.0) / 3.0);
    }

    return dt;
}

inline void init(){

    Generate_Mesh();    now_time=0.0;
    for(register int i=0;i<=N+1;++i)
    uh[i].resize(k+5),u1[i].resize(k+5),u2[i].resize(k+5),u3[i].resize(k+5);

    NID.Init(k);NID.Calculate_Gauss_Labatto_Points();NID.Get_Difference_Matrix();NID.Calculate_Lagrange_Derivative();

    for(register int i=1;i<=N;++i)
    {
        double left=mesh_point[i-1];
        double right=mesh_point[i];
        double hi=cell_h[i];
        for(register int j=0;j<=k;++j)
        {
            double nw_x=NID.get_GL()[j].first*hi*0.5+(left+right)*0.5;
            uh[i][j]=init_solution(nw_x);
        }
    }

    tau=compute_tau();
    return;
}

inline double RHS(int i,int j,const State& u){

    double sum=0.0;

    for(register int l=0;l<=k;++l)
    sum-=2.0*NID.get_DF()(j,l)*f_s(u[i][j],u[i][l]);

    double tmp=0.0;
    if(j==0){
        tmp=f(u[i][j])-f_hat(u[(i-1==0?N:i-1)][k],u[i][0]);
        tmp*=(-1.0/NID.get_GL()[j].second);
    }
    if(j==k){
        tmp=f(u[i][j])-f_hat(u[i][k],u[(i+1==N+1?1:i+1)][0]);
        tmp*=(1.0/NID.get_GL()[j].second);
    }

    sum+=tmp;

    double left=mesh_point[i-1];
    double right=mesh_point[i];
    double hi=cell_h[i];
    double nw_x=NID.get_GL()[j].first*hi*0.5+(left+right)*0.5;

    double ans=(2.0/cell_h[i])*(sum);
    if(OF_YES==1){
        double sigma = Damping_Coefficient(i,u);
        double Qu = u[i][j] - Cell_Average(i,u);

        ans -= sigma*Qu;
    }

    return ans;
}

inline void SSP(){
    
    for(;now_time<T_end;){
        if(now_time+tau<=T_end)now_time+=tau;
        else {tau=T_end-now_time;now_time=T_end;}

        sigma0=0.0;
        // Stage 1

        real_time=now_time-tau;
        for(int i=1;i<=N;++i){
            for(int j=0;j<=k;++j)u1[i][j]=uh[i][j]+tau*RHS(i,j,uh);
        }

        // Stage 2
   
        real_time=now_time;
        for(int i=1;i<=N;++i){
            for(int j=0;j<=k;++j)u2[i][j]=3.0/4.0*uh[i][j]+1.0/4.0*(u1[i][j]+tau*RHS(i,j,u1));
        }

        // Stage 3

        real_time=now_time-0.5*tau;
        for(int i=1;i<=N;++i){
            for(int j=0;j<=k;++j)u3[i][j]=1.0/3.0*uh[i][j]+2.0/3.0*(u2[i][j]+tau*RHS(i,j,u2));
        }

        for(int i=1;i<=N;++i)for(int j=0;j<=k;++j)uh[i][j]=u3[i][j];

        tau=compute_tau();
    }

    return;
}

inline void Error_Estimate(int round){

    for(register int i=1;i<=N;++i)
    {
        double left=mesh_point[i-1];
        double right=mesh_point[i];
        double hi=cell_h[i];
        for(register int j=0;j<=k;++j)
        {
            double nw_x=NID.get_GL()[j].first*hi*0.5+(left+right)*0.5;
            double tmp=fabs(uh[i][j]-true_solution(nw_x,now_time));
            norm_h[round]+=tmp*tmp*NID.get_GL()[j].second*hi*0.5;
        }
    }

    norm_h[round]=sqrt(norm_h[round]);
    // cout<<norm_h[round]<<" ";
    // if(round>1){
    //     cout<<log(norm_h[round-1]/norm_h[round])/log(2)<<endl;
    // }
    // else cout<<endl;

    return;
}


inline vector<pair<double,double>> Gauss_Legendre_16(){
    return {
        {-0.9894009349916499, 0.0271524594117541},
        {-0.9445750230732326, 0.0622535239386479},
        {-0.8656312023878318, 0.0951585116824928},
        {-0.7554044083550030, 0.1246289712555339},
        {-0.6178762444026438, 0.1495959888165767},
        {-0.4580167776572274, 0.1691565193950025},
        {-0.2816035507792589, 0.1826034150449236},
        {-0.0950125098376374, 0.1894506104550685},
        { 0.0950125098376374, 0.1894506104550685},
        { 0.2816035507792589, 0.1826034150449236},
        { 0.4580167776572274, 0.1691565193950025},
        { 0.6178762444026438, 0.1495959888165767},
        { 0.7554044083550030, 0.1246289712555339},
        { 0.8656312023878318, 0.0951585116824928},
        { 0.9445750230732326, 0.0622535239386479},
        { 0.9894009349916499, 0.0271524594117541}
    };
}

inline double Lagrange_Basis_Value(int j, double xi){
    double val = 1.0;
    const auto& gl = NID.get_GL();

    for(int l = 0; l <= k; ++l){
        if(l == j) continue;
        val *= (xi - gl[l].first) / (gl[j].first - gl[l].first);
    }

    return val;
}

inline double Evaluate_Uh_In_Cell(int i, double xi){
    double val = 0.0;

    for(int j = 0; j <= k; ++j){
        val += uh[i][j] * Lagrange_Basis_Value(j, xi);
    }

    return val;
}

inline void Error_Estimate_L2(int round){

    norm_L2[round] = 0.0;

    vector<pair<double,double>> quad = Gauss_Legendre_16();

    for(int i = 1; i <= N; ++i){
        double left  = mesh_point[i-1];
        double right = mesh_point[i];
        double hi    = cell_h[i];

        for(auto &qw : quad){
            double xi = qw.first;
            double w  = qw.second;

            double x = 0.5 * (left + right) + 0.5 * hi * xi;

            double uh_val = Evaluate_Uh_In_Cell(i, xi);
            double u_val  = true_solution(x, now_time);

            double diff = uh_val - u_val;

            norm_L2[round] += diff * diff * w * 0.5 * hi;
        }
    }

    norm_L2[round] = sqrt(norm_L2[round]);
    return;
}

int main(){

    // cerr << "Estimated shock time = " << estimate_shock_time() << endl;

    cout<<"Shock time: "<<estimate_shock_time()<<endl;

    cout<<"Please Input k and T_end: ";
    cin>>k>>T_end;
    cout<<"Please decide f_hat type(1:God 2:LF): ";
    cin>>hat_type;
    cout<<"Please decide whether use OF: ";
    cin>>OF_YES;

    N=16;

    for(int round=1;round<=Max_Round;++round,N<<=1)
    {
        // cout<<"ssssss"<<endl;
        init();
        SSP();
        Error_Estimate(round);
        Error_Estimate_L2(round);
        
        cout << N << "  " << norm_h[round];

        if(round > 1){
        cout << "  "
        << log(norm_h[round-1]/norm_h[round])/log(2.0);
        }
        cout<<endl;
    }

    return 0;
}