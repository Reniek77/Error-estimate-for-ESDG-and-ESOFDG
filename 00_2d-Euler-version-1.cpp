#include<bits/stdc++.h>
using namespace std;

/***************************************************************************************************************************************
Part 1: Basic Definition
****************************************************************************************************************************************/

const double PI=acos(-1.0);
const double GAS_GAMMA=1.4;
const double VORTEX_EPSILON=5.0/(2.0*PI);
const double VORTEX_ALPHA=0.5;
const double VORTEX_CENTER_X=10.0;
const double VORTEX_CENTER_Y=10.0;
const double DOMAIN_LENGTH=20.0;

using Matrix = vector<vector<double> >;

inline Matrix zeros(int nrow, int ncol){
    return Matrix(nrow, vector<double>(ncol, 0.0));
}
inline Matrix identity(int n) {
    Matrix I(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) I[i][i] = 1.0;
    return I;
}
inline Matrix transpose(const Matrix& A) {
    int nrow = (int)A.size();
    int ncol = (int)A[0].size();
    Matrix AT(ncol, vector<double>(nrow, 0.0));
    for (int i = 0; i < nrow; ++i)
        for (int j = 0; j < ncol; ++j)
            AT[j][i] = A[i][j];
    return AT;
}
inline Matrix matmul(const Matrix& A, const Matrix& B) {
    int nrow = (int)A.size();
    int nmid = (int)A[0].size();
    int ncol = (int)B[0].size();
    Matrix C(nrow, vector<double>(ncol, 0.0));
    for (int i = 0; i < nrow; ++i)
        for (int k = 0; k < nmid; ++k) {
            double aik = A[i][k];
            if (aik == 0.0) continue;
            for (int j = 0; j < ncol; ++j)
                C[i][j] += aik * B[k][j];
        }
    return C;
}
inline Matrix matmul_diag_right(const Matrix& A, const vector<double>& d) {
    int nrow = (int)A.size();
    int ncol = (int)A[0].size();
    assert((int)d.size() == ncol);
    Matrix C = A;
    for (int i = 0; i < nrow; ++i)
        for (int j = 0; j < ncol; ++j)
            C[i][j] *= d[j];
    return C;
}
inline Matrix matmul_diag_left(const vector<double>& d, const Matrix& A) {
    int nrow = (int)A.size();
    int ncol = (int)A[0].size();
    assert((int)d.size() == nrow);
    Matrix C = A;
    for (int i = 0; i < nrow; ++i)
        for (int j = 0; j < ncol; ++j)
            C[i][j] *= d[i];
    return C;
}
inline Matrix invert(const Matrix& A_in) {
    int n = (int)A_in.size();
    Matrix A = A_in;
    Matrix Inv = identity(n);

    for (int col = 0; col < n; ++col) {
        // pivot
        int pivot = col;
        double max_abs = std::fabs(A[col][col]);
        for (int i = col+1; i < n; ++i) {
            double val = std::fabs(A[i][col]);
            if (val > max_abs) {
                max_abs = val;
                pivot = i;
            }
        }
        assert(max_abs > 1e-14); // assume non-singular

        if (pivot != col) {
            std::swap(A[pivot], A[col]);
            std::swap(Inv[pivot], Inv[col]);
        }

        double diag = A[col][col];
        double inv_diag = 1.0 / diag;

        for (int j = 0; j < n; ++j) {
            A[col][j] *= inv_diag;
            Inv[col][j] *= inv_diag;
        }

        for (int i = 0; i < n; ++i) {
            if (i == col) continue;
            double factor = A[i][col];
            if (factor == 0.0) continue;
            for (int j = 0; j < n; ++j) {
                A[i][j] -= factor * A[col][j];
                Inv[i][j] -= factor * Inv[col][j];
            }
        }
    }
    return Inv;
}

using Vector=std::array<double,4>;

inline Vector operator+(const Vector& a,const Vector& b){
    Vector c{};
    for(int s=0;s<4;++s)c[s]=a[s]+b[s];
    return c;
}
inline Vector operator-(const Vector& a,const Vector& b){
    Vector c{};
    for(int s=0;s<4;++s)c[s]=a[s]-b[s];
    return c;
}
inline Vector operator*(double a,const Vector& b){
    Vector c{};
    for(int s=0;s<4;++s)c[s]=a*b[s];
    return c;
}
inline Vector operator*(const Vector& a,double b){
    return b*a;
}
inline Vector operator/(const Vector& a,double b){
    return (1.0/b)*a;
}
inline Vector& operator+=(Vector& a,const Vector& b){
    for(int s=0;s<4;++s)a[s]+=b[s];
    return a;
}
inline Vector& operator-=(Vector& a,const Vector& b){
    for(int s=0;s<4;++s)a[s]-=b[s];
    return a;
}

inline Vector Matrix_Vector_Product(const Matrix& A,const Vector& x){
    assert(A.size()==4&&A[0].size()==4);
    Vector b{};
    for(int i=0;i<4;++i)
        for(int j=0;j<4;++j)b[i]+=A[i][j]*x[j];
    return b;
}

inline double Average_Under_Log(double a,double b){
    assert(a>0.0&&b>0.0);
    double f=a/b-1.0;
    double u=f*f/((f+2.0)*(f+2.0));
    if(fabs(f)<1e-2){
        return (a+b)/(2.0+2.0/3.0*u+2.0/5.0*u*u);
    }
    return (a-b)/log(a/b);
}

class State_of_Fluid{
private:
    double density=0.0;
    double velocity_x=0.0;
    double velocity_y=0.0;
    double pressure=0.0;
    double momentum_x=0.0;
    double momentum_y=0.0;
    double energy=0.0;
    double enthalpy=0.0;
    Vector cal_data{};

public:
    State_of_Fluid()=default;

    void Init(double rho,double u,double v,double p){
        if(!(rho>0.0&&p>0.0)){
            throw runtime_error("Non-positive density or pressure in State_of_Fluid::Init().");
        }
        density=rho;
        velocity_x=u;
        velocity_y=v;
        pressure=p;
        momentum_x=rho*u;
        momentum_y=rho*v;
        energy=0.5*rho*(u*u+v*v)+p/(GAS_GAMMA-1.0);
        enthalpy=(energy+p)/rho;
        cal_data={density,momentum_x,momentum_y,energy};
    }

    void Set_Conserved(const Vector& U){
        cal_data=U;
        Synchronize();
    }

    void Synchronize(){
        density=cal_data[0];
        momentum_x=cal_data[1];
        momentum_y=cal_data[2];
        energy=cal_data[3];
        if(!(density>0.0)){
            throw runtime_error("Non-positive density in State_of_Fluid::Synchronize().");
        }
        velocity_x=momentum_x/density;
        velocity_y=momentum_y/density;
        pressure=(GAS_GAMMA-1.0)
            *(energy-0.5*density*(velocity_x*velocity_x+velocity_y*velocity_y));
        if(!(pressure>0.0)){
            throw runtime_error("Non-positive pressure in State_of_Fluid::Synchronize().");
        }
        enthalpy=(energy+pressure)/density;
    }

    double get_density()const{return density;}
    double get_velocity_x()const{return velocity_x;}
    double get_velocity_y()const{return velocity_y;}
    double get_pressure()const{return pressure;}
    double get_momentum_x()const{return momentum_x;}
    double get_momentum_y()const{return momentum_y;}
    double get_energy()const{return energy;}
    double get_enthalpy()const{return enthalpy;}
    const Vector& get_cal_data()const{return cal_data;}
};

inline Vector Physical_Flux_1(const State_of_Fluid& U){
    return {
        U.get_momentum_x(),
        U.get_density()*U.get_velocity_x()*U.get_velocity_x()+U.get_pressure(),
        U.get_density()*U.get_velocity_x()*U.get_velocity_y(),
        (U.get_energy()+U.get_pressure())*U.get_velocity_x()
    };
}

inline Vector Physical_Flux_2(const State_of_Fluid& U){
    return {
        U.get_momentum_y(),
        U.get_density()*U.get_velocity_x()*U.get_velocity_y(),
        U.get_density()*U.get_velocity_y()*U.get_velocity_y()+U.get_pressure(),
        (U.get_energy()+U.get_pressure())*U.get_velocity_y()
    };
}

using PairD=std::pair<double,double>;
inline PairD operator+(const PairD& a,const PairD& b){
    return {a.first+b.first,a.second+b.second};
}
inline PairD operator-(const PairD& a,const PairD& b){
    return {a.first-b.first,a.second-b.second};
}
inline PairD operator*(double s,const PairD& p){
    return {s*p.first,s*p.second};
}
inline PairD operator*(const PairD& p,double s){
    return {s*p.first,s*p.second};
}

static inline uint64_t EdgeKey(int a,int b){
    if(a>b)std::swap(a, b);
    return(uint64_t)(uint32_t)a<<32|(uint32_t)b;
}
inline void Strip_CR(string& line){
    if(!line.empty()&&line.back()=='\r')line.pop_back();
    return;
}

/***************************************************************************************************************************************
Part 2: Variable
****************************************************************************************************************************************/

const double CFL=0.1;
const string MESH_DIRECTORY="";
const bool USE_DAMPING=false; // false: ESDG; true: OFESDG

double now_time,stage_time;
double tau,hmin,hmax,T_end;

int N_GLP,k_poly;
std::vector<PairD>Gauss_1D;
std::vector<std::pair<PairD,double> >Gauss_2D;

std::vector<vector<State_of_Fluid> >u_h,u_1,u_2,u_3;

Matrix D_1,D_2;
std::vector<double>B_1,B_2;// in fact: diag matrix

// M-orthogonal P^k projection followed by evaluation at the three vertices.
// The rows correspond to v0=(0,0), v1=(0,1), v2=(1,0).
Matrix Projected_Value_At_Vertex;
Matrix Projected_Dr_At_Vertex;
Matrix Projected_Ds_At_Vertex;

int Nnode,Nele;
std::vector<std::array<int,3> >NeiEdge;  // Id in Neighbour
std::vector<std::array<int,3> >FaceFlip; // 0 same, 1 reversed, -1 unset

int is_egde[37];// -1: no 0:Left 1:Right 2:Down
int from_2d_to_1d[37];
int from_1d_to_2d[3][37];

struct Half_Edge_Info{
    int elem;     // element id
    int ledge;    // local edge index (0..2)
    int a,b;     // directed endpoints a -> b for this half-edge in that element
};

struct Node{
    int Id=0;
    std::pair<double,double>Coordinate={0.0,0.0};
};

struct Element{
    int Id=0;
    int Vertical[3]={0,0,0};
    int Neighbor[3]={-1,-1,-1};
    Matrix Jaco;
    double det_Jaco=0.0;
    double h_K=0.0;
};

std::vector<Node>Nod;
std::vector<Element>Ele;

std::vector<Vector>u_ave;
std::vector<double>damping_coef;
std::vector<std::array<Vector,3> >projected_vertex_value;
std::vector<std::array<Vector,3> >projected_vertex_dx;
std::vector<std::array<Vector,3> >projected_vertex_dy;
double sigma0=0.0;

/***************************************************************************************************************************************
Part 3: Mesh Read and Pre
****************************************************************************************************************************************/

inline bool Read_MSH2(const string& filename) {

    ifstream fin(filename.c_str());
    if (!fin) {
        cout << "Cannot open mesh file: " << filename << "\n";
        return false;
    }
    string line;

    Nnode=0;Nele=0;
    Nod.clear();Ele.clear();

    // Step 1: Read Nodes
    bool found_nodes=false;
    while(getline(fin,line)){
        Strip_CR(line);
        if(line=="$Nodes"){
            found_nodes=true;
            break;
        }
    }
    if(!found_nodes||!getline(fin,line)){
        cout << "No $Nodes block in " << filename << "\n";
        return false;
    }

    Strip_CR(line);
    int numNodes=stoi(line);
    Nod.assign(numNodes+1,Node());

    for(int i=0;i<numNodes;++i){
        if(!getline(fin,line))return false;
        Strip_CR(line);
        istringstream sin(line);

        int id;double x,y,z;
        sin>>id>>x>>y>>z;
        if(!sin){
            cout << "Invalid node line in " << filename << "\n";
            return false;
        }

        if(id>=(int)Nod.size())Nod.resize(id+1);
        Nod[id].Id=id;Nod[id].Coordinate={x,y};
        Nnode=max(Nnode,id);
    }

    // Step 2: Read Elements line by line, so non-triangle elements can be skipped safely.
    bool found_elements=false;
    while(getline(fin,line)){
        Strip_CR(line);
        if(line=="$Elements"){
            found_elements=true;
            break;
        }
    }
    if(!found_elements||!getline(fin,line)){
        cout << "No $Elements block in " << filename << "\n";
        return false;
    }

    Strip_CR(line);
    int numElements=stoi(line);
    Ele.assign(1,Element()); // The physical elements are numbered from 1.

    for(int i=0;i<numElements;++i){
        if(!getline(fin,line))return false;
        Strip_CR(line);
        istringstream sin(line);

        int elmId,elmType,numTags;
        sin>>elmId>>elmType>>numTags;
        if(!sin){
            cout << "Invalid element line in " << filename << "\n";
            return false;
        }

        for(int t=0;t<numTags;++t){
            int tag;sin>>tag;
        }

        if(elmType==2){
            Element K;
            K.Id=elmId;
            sin>>K.Vertical[0]>>K.Vertical[1]>>K.Vertical[2];
            if(!sin){
                cout << "Invalid triangle line in " << filename << "\n";
                return false;
            }
            Ele.push_back(K);
        }
    }

    Nele=(int)Ele.size()-1;

    return true;
}

inline void Scale_Mesh(double factor){
    for(int id=1;id<=Nnode;++id){
        Nod[id].Coordinate.first*=factor;
        Nod[id].Coordinate.second*=factor;
    }
}

inline void FixOrientation(){
    for(int id=1; id<=Nele; ++id)
    {
        int v0 = Ele[id].Vertical[0];
        int v1 = Ele[id].Vertical[1];
        int v2 = Ele[id].Vertical[2];

        auto p0 = Nod[v0].Coordinate;
        auto p1 = Nod[v1].Coordinate;
        auto p2 = Nod[v2].Coordinate;

        double det =(p2.first-p0.first)*(p1.second-p0.second)-(p1.first-p0.first)*(p2.second-p0.second);
        if(det < 0)std::swap(Ele[id].Vertical[1], Ele[id].Vertical[2]);
    }
    return;
}

inline void Build_Adjacency() {
    NeiEdge.assign(Nele+1,std::array<int,3>{-1,-1,-1});
    FaceFlip.assign(Nele+1,std::array<int,3>{-1,-1,-1});

    for(int id=1;id<=Nele;++id){
        for(int e=0;e<3;++e){
            Ele[id].Neighbor[e]=-1;
        }
    }

    unordered_map<uint64_t, Half_Edge_Info> mp;
    mp.reserve((size_t)(3ULL * Nele));
    mp.max_load_factor(0.7f);

    for (int id=1; id<=Nele; ++id){
        int v0=Ele[id].Vertical[0];
        int v1=Ele[id].Vertical[1];
        int v2=Ele[id].Vertical[2];

        int ea[3] = {v0, v1, v0};
        int eb[3] = {v1, v2, v2};// left skew down

        for(int e=0;e<3;++e){
            uint64_t key=EdgeKey(ea[e], eb[e]);
            auto it=mp.find(key);
            if(it==mp.end())mp.emplace(key, Half_Edge_Info{id,e, ea[e], eb[e]});
            else{
                int id2 = it->second.elem;
                int e2  = it->second.ledge;
                int a2  = it->second.a;
                int b2  = it->second.b;

                Ele[id].Neighbor[e]   = id2;  NeiEdge[id][e]   = e2;
                Ele[id2].Neighbor[e2] = id;   NeiEdge[id2][e2] = e;

                int a = ea[e], b = eb[e];
                int flip=-1;
                if (a == b2 && b == a2) flip = 1;        // opposite direction
                else if (a == a2 && b == b2) flip = 0;   // same direction

                if(flip==-1){
                    cout<<"Invalid orientation when pairing an internal edge.\n";
                    exit(1);
                }

                FaceFlip[id][e]   = flip;
                FaceFlip[id2][e2] = flip;

                mp.erase(it);
            }
        }
    }
    return;
}

inline void Get_Local_Edge_Vertices(int id,int e,int& a,int& b){
    int v0=Ele[id].Vertical[0];
    int v1=Ele[id].Vertical[1];
    int v2=Ele[id].Vertical[2];

    if(e==0){a=v0;b=v1;}
    else if(e==1){a=v1;b=v2;}
    else {a=v0;b=v2;}
}

/*********************************************************************************************************************************
Part 4: Data on Reference Triangles
**********************************************************************************************************************************/

inline PairD Get_Phy_Coordinate(int id, int poin){
    double x_ref = Gauss_2D[poin].first.first;
    double y_ref = Gauss_2D[poin].first.second;
    PairD v0 = Nod[Ele[id].Vertical[0]].Coordinate;
    PairD v1 = Nod[Ele[id].Vertical[1]].Coordinate;
    PairD v2 = Nod[Ele[id].Vertical[2]].Coordinate;
    double x_phy = (v2.first  - v0.first ) * x_ref + (v1.first  - v0.first ) * y_ref + v0.first;
    double y_phy = (v2.second - v0.second) * x_ref + (v1.second - v0.second) * y_ref + v0.second;
    return {x_phy, y_phy};
}

inline void Get_Quadrature_Rules_On_Reference_Triangles(int k){
    double ob111[5];double ob21[5];
    std::pair<double,double>vert[5];
    vert[0]={0.0,0.0};vert[1]={0.0,1.0};vert[2]={1.0,0.0};

    Gauss_1D.resize(k+1);

    if(k==1)
    {
        Gauss_1D[0]={-1.0/sqrt(3),1.0};Gauss_1D[1]={1.0/sqrt(3),1.0};

        N_GLP=6;
        Gauss_2D.resize(N_GLP);
        ob111[0]=0.0;ob111[1]=0.5-sqrt(3.0)/6.0;ob111[2]=1.0-ob111[0]-ob111[1];
        std::sort(ob111,ob111+3);
        int temp=0;
        do
        {
            Gauss_2D[temp]=make_pair(ob111[0]*vert[0]+ob111[1]*vert[1]+ob111[2]*vert[2],1.0/12.0);
            temp++;
        }while(next_permutation(ob111,ob111+3));
    }
    else if(k==2)
    {
        Gauss_1D[0]={-sqrt(0.6),5.0/9.0};Gauss_1D[1]={0.0,8.0/9.0};Gauss_1D[2]={sqrt(0.6),5.0/9.0};

        N_GLP=10;
        Gauss_2D.resize(N_GLP);
        ob111[0]=0.0;ob111[1]=0.5-sqrt(15.0)/10.0;ob111[2]=1.0-ob111[0]-ob111[1];
        std::sort(ob111,ob111+3);
        int temp=0;
        do
        {
            Gauss_2D[temp]=make_pair(ob111[0]*vert[0]+ob111[1]*vert[1]+ob111[2]*vert[2],1.0/48.0);
            temp++;
        }while(next_permutation(ob111,ob111+3));

        ob21[0]=0.5,ob21[1]=0.5,ob21[2]=0.0;
        std::sort(ob21,ob21+3);
        do
        {
            Gauss_2D[temp]=make_pair(ob21[0]*vert[0]+ob21[1]*vert[1]+ob21[2]*vert[2],1.0/20.0);
            temp++;
        }while(next_permutation(ob21,ob21+3));

        Gauss_2D[temp]={{1.0/3.0,1.0/3.0},9.0/40.0};
    }
    else if(k==3)
    {
        Gauss_1D[0] = {-sqrt(3.0/7.0 + 2.0/7.0*sqrt(6.0/5.0)), (18.0 - sqrt(30.0)) / 36.0};Gauss_1D[1] = {-sqrt(3.0/7.0 - 2.0/7.0*sqrt(6.0/5.0)), (18.0 + sqrt(30.0)) / 36.0};
        Gauss_1D[2] = { sqrt(3.0/7.0 - 2.0/7.0*sqrt(6.0/5.0)), (18.0 + sqrt(30.0)) / 36.0};Gauss_1D[3] = { sqrt(3.0/7.0 + 2.0/7.0*sqrt(6.0/5.0)), (18.0 - sqrt(30.0)) / 36.0};

        N_GLP=18;
        Gauss_2D.resize(N_GLP);
        ob111[0]=0.0;ob111[1]=0.330009478207572;ob111[2]=1.0-ob111[0]-ob111[1];
        std::sort(ob111,ob111+3);
        int temp=0;
        do
        {
            Gauss_2D[temp]=make_pair(ob111[0]*vert[0]+ob111[1]*vert[1]+ob111[2]*vert[2], 0.0202282703414950);
            temp++;
        }while(next_permutation(ob111,ob111+3));

        ob111[0]=0.0;ob111[1]=0.0694318442029737;ob111[2]=1.0-ob111[0]-ob111[1];
        std::sort(ob111,ob111+3);
        do
        {
            Gauss_2D[temp]=make_pair(ob111[0]*vert[0]+ob111[1]*vert[1]+ob111[2]*vert[2], 0.00754950743628280);
            temp++;
        }while(next_permutation(ob111,ob111+3));

        ob111[0]=0.1870738791912771;ob111[1]=0.5841571139756569;ob111[2]=1.0-ob111[0]-ob111[1];
        std::sort(ob111,ob111+3);
        do
        {
            Gauss_2D[temp]=make_pair(ob111[0]*vert[0]+ob111[1]*vert[1]+ob111[2]*vert[2], 0.0555555555555556);
            temp++;
        }while(next_permutation(ob111,ob111+3));
    }
    else if(k==4)
    {
        Gauss_1D[0] = {-0.9061798459386640, 0.2369268850561891};Gauss_1D[1] = {-0.5384693101056831, 0.4786286704993665};
        Gauss_1D[2] = { 0.0,                0.5688888888888889};
        Gauss_1D[3] = { 0.5384693101056831, 0.4786286704993665};Gauss_1D[4] = { 0.9061798459386640, 0.2369268850561891};

        N_GLP=22;
        Gauss_2D.resize(N_GLP);
        ob111[0]=0.0;ob111[1]=0.230765344947159;ob111[2]=1.0-ob111[0]-ob111[1];
        std::sort(ob111,ob111+3);
        int temp=0;
        do
        {
            Gauss_2D[temp]=make_pair(ob111[0]*vert[0]+ob111[1]*vert[1]+ob111[2]*vert[2], 0.0102652298402145);
            temp++;
        }while(next_permutation(ob111,ob111+3));

        ob111[0]=0.0;ob111[1]=0.046910077030668;ob111[2]=1.0-ob111[0]-ob111[1];
        std::sort(ob111,ob111+3);
        do
        {
            Gauss_2D[temp]=make_pair(ob111[0]*vert[0]+ob111[1]*vert[1]+ob111[2]*vert[2], 0.00330065754050081);
            temp++;
        }while(next_permutation(ob111,ob111+3));

        ob21[0]= 0.4384239524408185,ob21[1]= 0.4384239524408185,ob21[2]=1.0-ob21[0]-ob21[1];
        std::sort(ob21,ob21+3);
        do
        {
            Gauss_2D[temp]=make_pair(ob21[0]*vert[0]+ob21[1]*vert[1]+ob21[2]*vert[2], 0.0623683661448868);
            temp++;
        }while(next_permutation(ob21,ob21+3));

        ob21[0]= 0.1394337314154536,ob21[1]= 0.1394337314154536,ob21[2]=1.0-ob21[0]-ob21[1];
        std::sort(ob21,ob21+3);
        do
        {
            Gauss_2D[temp]=make_pair(ob21[0]*vert[0]+ob21[1]*vert[1]+ob21[2]*vert[2], 0.0527146648104222);
            temp++;
        }while(next_permutation(ob21,ob21+3));

        ob21[0]= 0.5,ob21[1]= 0.5,ob21[2]=1.0-ob21[0]-ob21[1];
        std::sort(ob21,ob21+3);
        do
        {
            Gauss_2D[temp]=make_pair(ob21[0]*vert[0]+ob21[1]*vert[1]+ob21[2]*vert[2], 0.00926854241697489);
            temp++;
        }while(next_permutation(ob21,ob21+3));

        Gauss_2D[temp]={{1.0/3.0,1.0/3.0}, 0.0455499555988567};
    }

    return;
}

inline void Get_Boundary_Matrix(int k){

    B_1.resize(N_GLP);B_2.resize(N_GLP);

    memset(is_egde,-1,sizeof(is_egde));

    for(int id=0;id<N_GLP;++id)
    {
        double x_ref=Gauss_2D[id].first.first;
        double y_ref=Gauss_2D[id].first.second;

        double tau_1=0.0;
        double tau_2=0.0;

        const double eps=1e-12;

        //Left: [0,1] \to [-1,1], y \to 2y-1
        if(fabs(x_ref)<=eps){
            double s=y_ref;//in [0,1]
            double w=2.0*s-1.0;// in [-1,1]

            int cor = 0;
            double best = 1e300;
            for(int j=0;j<=k;++j){
                double d = fabs(w - Gauss_1D[j].first);
                if(d < best){ best=d; cor=j; }
            }

            tau_1=-Gauss_1D[cor].second*0.5;
            tau_2=0.0;

            is_egde[id]=0;
            from_2d_to_1d[id]=cor;
            from_1d_to_2d[0][cor]=id;
        }

        //Down [0,1] \to [-1,1]
        else if(fabs(y_ref)<=eps){
            double s=x_ref;//in [0,1]
            double w=2.0*s-1.0;// in [-1,1]

            int cor = 0;
            double best = 1e300;
            for(int j=0;j<=k;++j){
                double d = fabs(w - Gauss_1D[j].first);
                if(d < best){ best=d; cor=j; }
            }

            tau_1=0.0;
            tau_2=-Gauss_1D[cor].second*0.5;

            is_egde[id]=2;
            from_2d_to_1d[id]=cor;
            from_1d_to_2d[2][cor]=id;
        }

        //Skew
        else if(fabs(x_ref+y_ref-1.0)<eps){
            double s=x_ref;//in [0,1]
            double w=2.0*s-1.0;// in [-1,1]

            int cor = 0;
            double best = 1e300;
            for(int j=0;j<=k;++j){
                double d = fabs(w - Gauss_1D[j].first);
                if(d < best){ best=d; cor=j; }
            }

            tau_1=Gauss_1D[cor].second*0.5;
            tau_2=Gauss_1D[cor].second*0.5;

            is_egde[id]=1;
            from_2d_to_1d[id]=cor;
            from_1d_to_2d[1][cor]=id;
        }

        B_1[id]=tau_1;
        B_2[id]=tau_2;
    }

    return;
}

inline double Binom_Real(double a, int m) {
    // C(a,0) = 1
    if (m == 0) return 1.0;
    double num = std::tgamma(a + 1.0);
    double den = std::tgamma(m + 1.0) * std::tgamma(a - m + 1.0);
    return num / den;
}

inline double JacobiP(int n, double alpha, double beta, double x) {
    if (n == 0) return 1.0;
    if (n == 1) {
        return 0.5*((alpha+beta+2.0)*x+(alpha-beta));
    }

    double sum = 0.0;
    double xm1_over_2 = 0.5 * (x - 1.0);
    double xp1_over_2 = 0.5 * (x + 1.0);

    for (int m = 0; m <= n; ++m) {
        double c1 = Binom_Real(n + alpha, m);
        double c2 = Binom_Real(n + beta,  n - m);
        double term = c1 * c2
                      * std::pow(xm1_over_2, n - m)
                      * std::pow(xp1_over_2, m);
        sum += term;
    }
    return sum;
}

inline double JacobiP_Deriv(int n, double alpha, double beta, double x) {
    if (n == 0) return 0.0;
    double factor = 0.5 * (n + alpha + beta + 1.0);
    // P_{n-1}^{(alpha+1,beta+1)}(x)
    double Pnm1 = JacobiP(n - 1, alpha + 1.0, beta + 1.0, x);
    return factor * Pnm1;
}

inline void DubinerBasis(double x, double y, int k, std::vector<double>& phi) {
    const double eps = 1e-14;
    double eta = 2.0 * y - 1.0;   // in [-1,1]
    double xi;
    if (1.0 - y > eps) {
        xi = 2.0 * x / (1.0 - y) - 1.0;  // in [-1,1]
    } else {
        // near y=1 -> the top vertex; xi not really used
        xi = -1.0;
    }

    int n_basis = (k + 1) * (k + 2) / 2;
    phi.assign(n_basis, 0.0);

    int idx = 0;
    for (int i = 0; i <= k; ++i) {
        double one_minus_y_pow_i = std::pow(1.0 - y, i);
        double Pi_xi = JacobiP(i, 0.0, 0.0, xi);

        for (int j = 0; j <= k - i; ++j) {
            double Pj_eta = JacobiP(j, 2.0 * i + 1.0, 0.0, eta);
            phi[idx] = Pi_xi * one_minus_y_pow_i * Pj_eta;
            ++idx;
        }
    }
    return;
}

inline void DubinerBasisDxDy(double x, double y, int k,std::vector<double>& phi_x,std::vector<double>& phi_y) {
    // Compute d/dx and d/dy of Dubiner basis at (x,y).
    // Output phi_x, phi_y same ordering as DubinerBasis.
    const double eps = 1e-14;
    double eta = 2.0 * y - 1.0;
    double xi;
    bool near_top = (1.0 - y <= eps);

    if (!near_top) {
        xi = 2.0 * x / (1.0 - y) - 1.0;
    } else {
        // near top vertex, we treat as xi=-1 in a limiting sense
        xi = -1.0;
    }

    int n_basis = (k + 1) * (k + 2) / 2;
    phi_x.assign(n_basis, 0.0);
    phi_y.assign(n_basis, 0.0);

    // Precompute geometric derivatives if not at top vertex
    double xi_x = 0.0, xi_y = 0.0;
    double eta_y = 2.0;  // eta = 2y - 1
    if (!near_top) {
        double one_minus_y = 1.0 - y;
        xi_x = 2.0 / one_minus_y;
        xi_y = 2.0 * x / (one_minus_y * one_minus_y);
    }

    int idx = 0;
    for (int i = 0; i <= k; ++i) {
        // g(y) = (1 - y)^i
        double one_minus_y = 1.0 - y;
        double g = std::pow(one_minus_y, i);

        // g_y
        double g_y = 0.0;
        if (i > 0) {
            g_y = -i * std::pow(one_minus_y, i - 1);
        }

        // A(xi) = P_i^{(0,0)}(xi), A_xi = d/dxi P_i^{(0,0)}(xi)
        double A = JacobiP(i, 0.0, 0.0, xi);
        double A_xi = JacobiP_Deriv(i, 0.0, 0.0, xi);

        // For fixed i, we will reuse P_j^{(2i+1,0)}(eta) and its derivative.
        for (int j = 0; j <= k - i; ++j) {
            // B(eta) = P_j^{(2i+1,0)}(eta), B_eta = d/deta
            double B = JacobiP(j, 2.0 * i + 1.0, 0.0, eta);
            double B_eta = JacobiP_Deriv(j, 2.0 * i + 1.0, 0.0, eta);

            // ----- d/dx -----
            double dphi_dx = 0.0;
            if (!near_top) {
                // phi_x = A_x g B,  A_x = A_xi * xi_x
                double A_x = A_xi * xi_x;
                dphi_dx = A_x * g * B;
            } else {
                // near top vertex, contribution tends to 0 for i>0,
                // and for i=0 it is harmless (single point)
                dphi_dx = 0.0;
            }

            // ----- d/dy -----
            double dphi_dy = 0.0;
            if (!near_top) {
                // A_y = A_xi * xi_y
                double A_y = A_xi * xi_y;
                // B_y = B_eta * eta_y (eta_y=2)
                double B_y = B_eta * eta_y;

                // phi_y = A_y*g*B + A*g_y*B + A*g*B_y
                dphi_dy = A_y * g * B
                          + A * g_y * B
                          + A * g * B_y;
            } else {
                // y ~ 1, one_minus_y ~ 0, most modes vanish; we can set 0 safely
                dphi_dy = 0.0;
            }

            phi_x[idx] = dphi_dx;
            phi_y[idx] = dphi_dy;
            ++idx;
        }
    }
    return;
}

inline void Build_Diff_Matrix()
{
    int k = k_poly;
    int N = N_GLP;  // number of quadrature points (nodes)
    int n_basis = (k + 1) * (k + 2) / 2;

    Matrix V(N, vector<double>(n_basis, 0.0));
    Matrix Vx(N, vector<double>(n_basis, 0.0));
    Matrix Vy(N, vector<double>(n_basis, 0.0));

    std::vector<double> phi, phi_x, phi_y;

    for (int j = 0; j < N; ++j) {
        double x = Gauss_2D[j].first.first;
        double y = Gauss_2D[j].first.second;

        DubinerBasis(x, y, k, phi);
        DubinerBasisDxDy(x, y, k, phi_x, phi_y);

        for (int l = 0; l < n_basis; ++l) {
            V[j][l]  = phi[l];
            Vx[j][l] = phi_x[l];
            Vy[j][l] = phi_y[l];
        }
    }

    vector<double> w(N), w_inv(N);
    for (int j = 0; j < N; ++j) {
        w[j]     = Gauss_2D[j].second;
        w_inv[j] = 1.0 / w[j];
    }

    Matrix G(n_basis, vector<double>(n_basis, 0.0));
    for (int a = 0; a < n_basis; ++a) {
        for (int b = 0; b < n_basis; ++b) {
            double sum = 0.0;
            for (int j = 0; j < N; ++j) {
                sum += V[j][a] * w[j] * V[j][b];
            }
            G[a][b] = sum;
        }
    }

    Matrix R = zeros(n_basis, n_basis);
    for (int i = 0; i < n_basis; ++i) {
        for (int j = i; j < n_basis; ++j) {
            double sum = G[i][j];
            for (int k2 = 0; k2 < i; ++k2)
                sum -= R[k2][i] * R[k2][j];
            if (i == j) {
                R[i][i] = std::sqrt(sum);
            } else {
                R[i][j] = sum / R[i][i];
            }
        }
    }
    Matrix R_inv = invert(R);

    Matrix T = R_inv;

    Matrix V_ortho  = matmul(V,  T);
    Matrix Vx_ortho = matmul(Vx, T);
    Matrix Vy_ortho = matmul(Vy, T);

    V  = V_ortho;
    Vx = Vx_ortho;
    Vy = Vy_ortho;

    Matrix VVt(N, vector<double>(N, 0.0));
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            double sum = 0.0;
            for (int l = 0; l < n_basis; ++l)
                sum += V[i][l] * V[j][l];
            VVt[i][j] = sum;
        }

    Matrix A(N, vector<double>(N, 0.0));
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            double val = VVt[i][j];
            if (i == j) val += w_inv[i]; 
            A[i][j] = val;
        }
    Matrix Cmat(N, vector<double>(N, 0.0));
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            double vvM = VVt[i][j] * w[j];
            Cmat[i][j] = (i == j ? 1.0 : 0.0) - vvM;
        }

    Matrix VxVtM(N, vector<double>(N, 0.0));
    Matrix VyVtM(N, vector<double>(N, 0.0));

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            double sumx = 0.0;
            double sumy = 0.0;
            for (int l = 0; l < n_basis; ++l) {
                sumx += Vx[i][l] * V[j][l];
                sumy += Vy[i][l] * V[j][l];
            }
            VxVtM[i][j] = sumx * w[j];
            VyVtM[i][j] = sumy * w[j];
        }

    Matrix AB1(N, vector<double>(N, 0.0));
    Matrix AB2(N, vector<double>(N, 0.0));

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            AB1[i][j] = A[i][j] * B_1[j];
            AB2[i][j] = A[i][j] * B_2[j];
        }

    Matrix AB1C(N, vector<double>(N, 0.0));
    Matrix AB2C(N, vector<double>(N, 0.0));

    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            double sum1 = 0.0;
            double sum2 = 0.0;
            for (int l = 0; l < N; ++l) {
                sum1 += AB1[i][l] * Cmat[l][j];
                sum2 += AB2[i][l] * Cmat[l][j];
            }
            AB1C[i][j] = sum1;
            AB2C[i][j] = sum2;
        }

    // 10) D1, D2
    D_1 = zeros(N, N);
    D_2 = zeros(N, N);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) {
            D_1[i][j] = 0.5 * AB1C[i][j] + VxVtM[i][j];
            D_2[i][j] = 0.5 * AB2C[i][j] + VyVtM[i][j];
        }
    return;
}

inline double Integer_Power(double x,int p){
    double ans=1.0;
    for(int i=0;i<p;++i)ans*=x;
    return ans;
}

inline void Build_Damping_Projection_Operators(){
    // The OF coefficient uses traces of the M-orthogonal P^k projection of
    // the nodal vector.  A Bernstein basis on the reference triangle keeps
    // the small Gram matrix well conditioned and is regular at all vertices.
    vector<array<int,3> >modes;
    for(int a=0;a<=k_poly;++a){
        for(int b=0;b<=k_poly-a;++b){
            int c=k_poly-a-b;
            modes.push_back({a,b,c});
        }
    }

    int n_basis=(int)modes.size();
    Matrix V=zeros(N_GLP,n_basis);
    for(int j=0;j<N_GLP;++j){
        double r=Gauss_2D[j].first.first;
        double s=Gauss_2D[j].first.second;
        double lambda[3]={1.0-r-s,s,r};
        for(int m=0;m<n_basis;++m){
            int a=modes[m][0];
            int b=modes[m][1];
            int c=modes[m][2];
            double factor=tgamma(k_poly+1.0)
                /(tgamma(a+1.0)*tgamma(b+1.0)*tgamma(c+1.0));
            V[j][m]=factor*Integer_Power(lambda[0],a)
                *Integer_Power(lambda[1],b)*Integer_Power(lambda[2],c);
        }
    }

    Matrix Gram=zeros(n_basis,n_basis);
    for(int a=0;a<n_basis;++a){
        for(int b=0;b<n_basis;++b){
            for(int j=0;j<N_GLP;++j){
                Gram[a][b]+=Gauss_2D[j].second*V[j][a]*V[j][b];
            }
        }
    }
    Matrix Gram_Inverse=invert(Gram);

    // Coefficient_Operator maps nodal values to the coefficients of their
    // M-orthogonal P^k projection: c=(V^T M V)^{-1}V^T M u.
    Matrix Coefficient_Operator=zeros(n_basis,N_GLP);
    for(int m=0;m<n_basis;++m){
        for(int j=0;j<N_GLP;++j){
            for(int q=0;q<n_basis;++q){
                Coefficient_Operator[m][j]
                    +=Gram_Inverse[m][q]*V[j][q]*Gauss_2D[j].second;
            }
        }
    }

    const PairD reference_vertex[3]={{0.0,0.0},{0.0,1.0},{1.0,0.0}};
    Projected_Value_At_Vertex=zeros(3,N_GLP);
    Projected_Dr_At_Vertex=zeros(3,N_GLP);
    Projected_Ds_At_Vertex=zeros(3,N_GLP);

    for(int vertex=0;vertex<3;++vertex){
        double r=reference_vertex[vertex].first;
        double s=reference_vertex[vertex].second;
        double lambda[3]={1.0-r-s,s,r};

        for(int m=0;m<n_basis;++m){
            int a=modes[m][0];
            int b=modes[m][1];
            int c=modes[m][2];
            double factor=tgamma(k_poly+1.0)
                /(tgamma(a+1.0)*tgamma(b+1.0)*tgamma(c+1.0));

            double value=factor*Integer_Power(lambda[0],a)
                *Integer_Power(lambda[1],b)*Integer_Power(lambda[2],c);
            double common_a=(a==0?0.0:a*Integer_Power(lambda[0],a-1)
                *Integer_Power(lambda[1],b)*Integer_Power(lambda[2],c));
            double common_b=(b==0?0.0:b*Integer_Power(lambda[0],a)
                *Integer_Power(lambda[1],b-1)*Integer_Power(lambda[2],c));
            double common_c=(c==0?0.0:c*Integer_Power(lambda[0],a)
                *Integer_Power(lambda[1],b)*Integer_Power(lambda[2],c-1));
            double dr=factor*(-common_a+common_c);
            double ds=factor*(-common_a+common_b);

            for(int j=0;j<N_GLP;++j){
                Projected_Value_At_Vertex[vertex][j]
                    +=value*Coefficient_Operator[m][j];
                Projected_Dr_At_Vertex[vertex][j]
                    +=dr*Coefficient_Operator[m][j];
                Projected_Ds_At_Vertex[vertex][j]
                    +=ds*Coefficient_Operator[m][j];
            }
        }
    }
    return;
}

/*********************************************************************************************************************************
Part 5: Flux and Functions
**********************************************************************************************************************************/

struct Exact_Euler_Data{
    State_of_Fluid state;
    Vector derivative_x{};
    Vector derivative_y{};
};

inline Exact_Euler_Data Exact_Solution_Data(PairD Coordinate,double t){
    const double xi=Coordinate.first-t-VORTEX_CENTER_X;
    const double eta=Coordinate.second-t-VORTEX_CENTER_Y;
    const double r_squared=xi*xi+eta*eta;

    const double phi=VORTEX_EPSILON
        *exp(VORTEX_ALPHA*(1.0-r_squared));
    const double phi_x=-2.0*VORTEX_ALPHA*xi*phi;
    const double phi_y=-2.0*VORTEX_ALPHA*eta*phi;

    const double temperature=1.0-(GAS_GAMMA-1.0)/(2.0*GAS_GAMMA)*phi*phi;
    if(!(temperature>0.0)){
        throw runtime_error("The isentropic-vortex temperature is non-positive.");
    }

    const double temperature_x=-(GAS_GAMMA-1.0)/GAS_GAMMA*phi*phi_x;
    const double temperature_y=-(GAS_GAMMA-1.0)/GAS_GAMMA*phi*phi_y;

    const double rho=pow(temperature,1.0/(GAS_GAMMA-1.0));
    const double p=pow(temperature,GAS_GAMMA/(GAS_GAMMA-1.0));
    const double velocity_x=1.0-eta*phi;
    const double velocity_y=1.0+xi*phi;

    const double rho_x=rho/(GAS_GAMMA-1.0)*temperature_x/temperature;
    const double rho_y=rho/(GAS_GAMMA-1.0)*temperature_y/temperature;
    const double p_x=GAS_GAMMA/(GAS_GAMMA-1.0)*p*temperature_x/temperature;
    const double p_y=GAS_GAMMA/(GAS_GAMMA-1.0)*p*temperature_y/temperature;

    const double velocity_x_x=-eta*phi_x;
    const double velocity_x_y=-phi-eta*phi_y;
    const double velocity_y_x= phi+xi*phi_x;
    const double velocity_y_y= xi*phi_y;

    Exact_Euler_Data exact;
    exact.state.Init(rho,velocity_x,velocity_y,p);

    exact.derivative_x[0]=rho_x;
    exact.derivative_x[1]=rho_x*velocity_x+rho*velocity_x_x;
    exact.derivative_x[2]=rho_x*velocity_y+rho*velocity_y_x;
    exact.derivative_x[3]=p_x/(GAS_GAMMA-1.0)
        +0.5*rho_x*(velocity_x*velocity_x+velocity_y*velocity_y)
        +rho*(velocity_x*velocity_x_x+velocity_y*velocity_y_x);

    exact.derivative_y[0]=rho_y;
    exact.derivative_y[1]=rho_y*velocity_x+rho*velocity_x_y;
    exact.derivative_y[2]=rho_y*velocity_y+rho*velocity_y_y;
    exact.derivative_y[3]=p_y/(GAS_GAMMA-1.0)
        +0.5*rho_y*(velocity_x*velocity_x+velocity_y*velocity_y)
        +rho*(velocity_x*velocity_x_y+velocity_y*velocity_y_y);

    return exact;
}

inline State_of_Fluid Init_Solution(PairD Coordinate){
    return Exact_Solution_Data(Coordinate,0.0).state;
}

inline State_of_Fluid Exact_Solution(PairD Coordinate,double t){
    return Exact_Solution_Data(Coordinate,t).state;
}

inline Vector EC_FLUX_1(const State_of_Fluid& uL,const State_of_Fluid& uR){
    const double rhoL=uL.get_density();
    const double rhoR=uR.get_density();
    const double uL_x=uL.get_velocity_x();
    const double uR_x=uR.get_velocity_x();
    const double uL_y=uL.get_velocity_y();
    const double uR_y=uR.get_velocity_y();
    const double betaL=rhoL/(2.0*uL.get_pressure());
    const double betaR=rhoR/(2.0*uR.get_pressure());

    const double rho_bar=0.5*(rhoL+rhoR);
    const double u_bar=0.5*(uL_x+uR_x);
    const double v_bar=0.5*(uL_y+uR_y);
    const double beta_bar=0.5*(betaL+betaR);
    const double rho_log=Average_Under_Log(rhoL,rhoR);
    const double beta_log=Average_Under_Log(betaL,betaR);
    const double p_tilde=rho_bar/(2.0*beta_bar);
    const double p_log_tilde=rho_log/(2.0*beta_log);
    const double velocity_square_bar=0.5*(
        uL_x*uL_x+uR_x*uR_x+uL_y*uL_y+uR_y*uR_y
    );
    const double energy_tilde=0.5*rho_log
        *(2.0*(u_bar*u_bar+v_bar*v_bar)-velocity_square_bar)+p_tilde;

    return {
        rho_log*u_bar,
        rho_log*u_bar*u_bar+p_tilde,
        rho_log*u_bar*v_bar,
        u_bar*(p_log_tilde/(GAS_GAMMA-1.0)+energy_tilde)
    };
}

inline Vector EC_FLUX_2(const State_of_Fluid& uL,const State_of_Fluid& uR){
    const double rhoL=uL.get_density();
    const double rhoR=uR.get_density();
    const double uL_x=uL.get_velocity_x();
    const double uR_x=uR.get_velocity_x();
    const double uL_y=uL.get_velocity_y();
    const double uR_y=uR.get_velocity_y();
    const double betaL=rhoL/(2.0*uL.get_pressure());
    const double betaR=rhoR/(2.0*uR.get_pressure());

    const double rho_bar=0.5*(rhoL+rhoR);
    const double u_bar=0.5*(uL_x+uR_x);
    const double v_bar=0.5*(uL_y+uR_y);
    const double beta_bar=0.5*(betaL+betaR);
    const double rho_log=Average_Under_Log(rhoL,rhoR);
    const double beta_log=Average_Under_Log(betaL,betaR);
    const double p_tilde=rho_bar/(2.0*beta_bar);
    const double p_log_tilde=rho_log/(2.0*beta_log);
    const double velocity_square_bar=0.5*(
        uL_x*uL_x+uR_x*uR_x+uL_y*uL_y+uR_y*uR_y
    );
    const double energy_tilde=0.5*rho_log
        *(2.0*(u_bar*u_bar+v_bar*v_bar)-velocity_square_bar)+p_tilde;

    return {
        rho_log*v_bar,
        rho_log*u_bar*v_bar,
        rho_log*v_bar*v_bar+p_tilde,
        v_bar*(p_log_tilde/(GAS_GAMMA-1.0)+energy_tilde)
    };
}

inline State_of_Fluid get_u_out(
    int id,
    int poin,
    const std::vector<std::vector<State_of_Fluid>>& v
){
    int e = is_egde[poin];
    int neib   = Ele[id].Neighbor[e];
    int neib_e = NeiEdge[id][e];

    if(neib<=0){
        return Exact_Solution(Get_Phy_Coordinate(id,poin),stage_time);
    }
    int idx = from_2d_to_1d[poin]; // 0..k_poly
    int idx2 = idx;

    int flip = FaceFlip[id][e];
    if (flip == 1) idx2 = k_poly - idx;  // only reverse when needed

    int neib_p = from_1d_to_2d[neib_e][idx2];

    return v[neib][neib_p];
}

inline double Local_LF_Alpha(
    const State_of_Fluid& uL,
    const State_of_Fluid& uR,
    double nx,
    double ny
){
    const double unL=uL.get_velocity_x()*nx+uL.get_velocity_y()*ny;
    const double unR=uR.get_velocity_x()*nx+uR.get_velocity_y()*ny;
    const double pL=uL.get_pressure();
    const double pR=uR.get_pressure();
    const double cL=sqrt(GAS_GAMMA*uL.get_pressure()/uL.get_density());
    const double cR=sqrt(GAS_GAMMA*uR.get_pressure()/uR.get_density());

    // Two-rarefaction approximation applied to the normal Riemann problem.
    const double pressure_exponent=
        (GAS_GAMMA-1.0)/(2.0*GAS_GAMMA);
    const double numerator=cL+cR
        +0.5*(GAS_GAMMA-1.0)*(unL-unR);
    const double denominator=
        cL/pow(pL,pressure_exponent)
        +cR/pow(pR,pressure_exponent);

    // A non-positive numerator corresponds to the vacuum case; p_star=0
    // then gives q_L=q_R=1 in the formulas below.
    double p_star=0.0;
    if(numerator>0.0){
        p_star=pow(
            numerator/denominator,
            2.0*GAS_GAMMA/(GAS_GAMMA-1.0)
        );
    }

    auto q=[&](double p){
        if(p_star<=p)return 1.0;
        return sqrt(
            1.0+(GAS_GAMMA+1.0)/(2.0*GAS_GAMMA)
                *(p_star/p-1.0)
        );
    };

    const double lambda_left=unL-cL*q(pL);
    const double lambda_right=unR+cR*q(pR);
    return max(fabs(lambda_left),fabs(lambda_right));
}

inline Vector ES_FLUX(
    int id,
    int poin,
    const vector<vector<State_of_Fluid>>& v
)
{
    double t1 = B_1[poin];
    double t2 = B_2[poin];

    const double xr = Ele[id].Jaco[0][0];
    const double xs = Ele[id].Jaco[0][1];
    const double yr = Ele[id].Jaco[1][0];
    const double ys = Ele[id].Jaco[1][1];

    double Nx = ys*t1 - yr*t2;
    double Ny = -xs*t1 + xr*t2;

    double Len = sqrt(Nx*Nx+Ny*Ny);

    double nx = Nx/Len;
    double ny = Ny/Len;

    const State_of_Fluid& uL=v[id][poin];
    State_of_Fluid uR=get_u_out(id,poin,v);

    Vector fnL=nx*Physical_Flux_1(uL)+ny*Physical_Flux_2(uL);
    Vector fnR=nx*Physical_Flux_1(uR)+ny*Physical_Flux_2(uR);

    double alpha=Local_LF_Alpha(uL,uR,nx,ny);

    return 0.5*(fnL+fnR)
        -0.5*alpha*(uR.get_cal_data()-uL.get_cal_data());
}

/*********************************************************************************************************************************
Part 6: Evolution
**********************************************************************************************************************************/

inline void Get_Cell_Average(
    const std::vector<std::vector<State_of_Fluid>>& u
){

    for(int id=1;id<=Nele;++id)
    {
        Vector num{};
        double den=0.0;
        for(int j=0;j<N_GLP;++j){
            double w = Gauss_2D[j].second;
            num+=w*u[id][j].get_cal_data();
            den += w;
        }
        u_ave[id]=num/den;
    }
    return;
}

inline void Get_Mesh_Size(){
    hmin=1e100;hmax=0.0;

    for(int id=1;id<=Nele;++id){
        PairD p[3];
        for(int j=0;j<3;++j)p[j]=Nod[Ele[id].Vertical[j]].Coordinate;

        double hK=0.0;
        for(int j=0;j<3;++j){
            PairD d=p[(j+1)%3]-p[j];
            hK=max(hK,sqrt(d.first*d.first+d.second*d.second));
        }

        Ele[id].h_K=hK;
        hmin=min(hmin,hK);
        hmax=max(hmax,hK);
    }
    return;
}

inline double estimate_lambda0(
    const std::vector<std::vector<State_of_Fluid>>& v
){
    double lambda0=0.0;

    for(int id=1;id<=Nele;++id){
        for(int j=0;j<N_GLP;++j){
            const State_of_Fluid& state=v[id][j];
            double sound_speed=sqrt(
                GAS_GAMMA*state.get_pressure()/state.get_density()
            );
            // rho(A_1)+rho(A_2)
            // = |w_1|+c+|w_2|+c for the two-dimensional Euler equations.
            double speed=fabs(state.get_velocity_x())
                +fabs(state.get_velocity_y())+2.0*sound_speed;
            lambda0=max(lambda0,speed);
        }
    }
    return lambda0;
}

inline double compute_tau(){
    double lambda0=estimate_lambda0(u_h);
    double CFL_nw=CFL/(lambda0+sigma0);

    if(k_poly<=2)return CFL_nw*hmin;
    return CFL_nw*pow(hmin,(k_poly+1.0)/3.0);
}

inline void Get_Projected_Vertex_Data(
    const std::vector<std::vector<State_of_Fluid>>& v
){
    for(int id=1;id<=Nele;++id){
        const Matrix& Jaco=Ele[id].Jaco;
        const double det_Jaco=Ele[id].det_Jaco;

        const double xr=Jaco[0][0];
        const double xs=Jaco[0][1];
        const double yr=Jaco[1][0];
        const double ys=Jaco[1][1];

        for(int vertex=0;vertex<3;++vertex){
            Vector value{};
            Vector derivative_r{};
            Vector derivative_s{};

            for(int j=0;j<N_GLP;++j){
                const Vector& U=v[id][j].get_cal_data();
                value+=Projected_Value_At_Vertex[vertex][j]*U;
                derivative_r+=Projected_Dr_At_Vertex[vertex][j]*U;
                derivative_s+=Projected_Ds_At_Vertex[vertex][j]*U;
            }

            projected_vertex_value[id][vertex]=value;
            projected_vertex_dx[id][vertex]
                =(ys*derivative_r-yr*derivative_s)/det_Jaco;
            projected_vertex_dy[id][vertex]
                =(-xs*derivative_r+xr*derivative_s)/det_Jaco;
        }
    }
    return;
}

inline PairD Get_Edge_Unit_Normal(int id,int edge){
    int a,b;
    Get_Local_Edge_Vertices(id,edge,a,b);
    PairD p=Nod[a].Coordinate;
    PairD q=Nod[b].Coordinate;

    PairD midpoint=0.5*(p+q);
    PairD centroid=(
        Nod[Ele[id].Vertical[0]].Coordinate
        +Nod[Ele[id].Vertical[1]].Coordinate
        +Nod[Ele[id].Vertical[2]].Coordinate
    )*(1.0/3.0);

    PairD tangent=q-p;
    PairD normal={tangent.second,-tangent.first};
    PairD to_centroid=centroid-midpoint;
    if(normal.first*to_centroid.first+normal.second*to_centroid.second>0.0){
        normal=-1.0*normal;
    }

    double length=sqrt(normal.first*normal.first+normal.second*normal.second);
    return normal*(1.0/length);
}

inline Matrix Euler_Left_Characteristic_Matrix(
    const Vector& UL,
    const Vector& UR,
    double nx,
    double ny
){
    State_of_Fluid state_L,state_R;
state_L.Set_Conserved(UL);
state_R.Set_Conserved(UR);

const double sqrt_rho_L=sqrt(state_L.get_density());
const double sqrt_rho_R=sqrt(state_R.get_density());
const double denominator=sqrt_rho_L+sqrt_rho_R;

const double w1=(
    sqrt_rho_L*state_L.get_velocity_x()
    +sqrt_rho_R*state_R.get_velocity_x()
)/denominator;

const double w2=(
    sqrt_rho_L*state_L.get_velocity_y()
    +sqrt_rho_R*state_R.get_velocity_y()
)/denominator;

const double enthalpy=(
    sqrt_rho_L*state_L.get_enthalpy()
    +sqrt_rho_R*state_R.get_enthalpy()
)/denominator;

const double c=sqrt(
    (GAS_GAMMA-1.0)
    *(enthalpy-0.5*(w1*w1+w2*w2))
);
    const double gamma_minus_one=GAS_GAMMA-1.0;
    const double b2=0.5*gamma_minus_one*(w1*w1+w2*w2);
    const double wn=w1*nx+w2*ny;
    const double wt=-w1*ny+w2*nx;
    const double factor=gamma_minus_one/c;

    Matrix L=zeros(4,4);

    L[0][0]=0.5*(b2+wn*c);
    L[0][1]=-0.5*(gamma_minus_one*w1+nx*c);
    L[0][2]=-0.5*(gamma_minus_one*w2+ny*c);
    L[0][3]=0.5*gamma_minus_one;

    L[1][0]=wt*c;
    L[1][1]=ny*c;
    L[1][2]=-nx*c;

    L[2][0]=c*c-b2;
    L[2][1]=gamma_minus_one*w1;
    L[2][2]=gamma_minus_one*w2;
    L[2][3]=-gamma_minus_one;

    L[3][0]=0.5*(b2-wn*c);
    L[3][1]=-0.5*(gamma_minus_one*w1-nx*c);
    L[3][2]=-0.5*(gamma_minus_one*w2-ny*c);
    L[3][3]=0.5*gamma_minus_one;

    for(int i=0;i<4;++i)
        for(int j=0;j<4;++j)
            L[i][j]*=factor;

    return L;
}

inline void Get_Sigma(
    const std::vector<std::vector<State_of_Fluid>>& v
){
    sigma0=0.0;

    if(!USE_DAMPING){
        fill(damping_coef.begin(),damping_coef.end(),0.0);
        return;
    }

    Get_Cell_Average(v);
    Get_Projected_Vertex_Data(v);

    // Local edge endpoints in the same direction used by FaceFlip.
    const int edge_vertex[3][2]={{0,1},{1,2},{0,2}};
    const double number_of_faces=3.0;

    for(int id=1;id<=Nele;++id){
        Vector sigma_squared{};

        for(int e=0;e<3;++e){
            int neighbor=Ele[id].Neighbor[e];
            PairD normal=Get_Edge_Unit_Normal(id,e);

            for(int endpoint=0;endpoint<2;++endpoint){
                int vertex=edge_vertex[e][endpoint];

                Vector outside_value{};
                Vector outside_dx{};
                Vector outside_dy{};

                if(neighbor>0){
                    int neighbor_edge=NeiEdge[id][e];
                    int flip=FaceFlip[id][e];
                    int neighbor_endpoint=(flip==1?1-endpoint:endpoint);
                    int neighbor_vertex=
                        edge_vertex[neighbor_edge][neighbor_endpoint];

                    outside_value=projected_vertex_value[neighbor][neighbor_vertex];
                    outside_dx=projected_vertex_dx[neighbor][neighbor_vertex];
                    outside_dy=projected_vertex_dy[neighbor][neighbor_vertex];
                }
                else{
                    PairD physical_vertex=
                        Nod[Ele[id].Vertical[vertex]].Coordinate;
                    Exact_Euler_Data exact=
                        Exact_Solution_Data(physical_vertex,stage_time);
                    outside_value=exact.state.get_cal_data();
                    outside_dx=exact.derivative_x;
                    outside_dy=exact.derivative_y;
                }

                Vector jump_value=
                    projected_vertex_value[id][vertex]-outside_value;
                Vector jump_dx=projected_vertex_dx[id][vertex]-outside_dx;
                Vector jump_dy=projected_vertex_dy[id][vertex]-outside_dy;

                Matrix L=Euler_Left_Characteristic_Matrix(
                    projected_vertex_value[id][vertex],
                    outside_value,
                    normal.first,
                    normal.second
                );
                Vector characteristic_jump=Matrix_Vector_Product(L,jump_value);
                Vector characteristic_jump_dx=Matrix_Vector_Product(L,jump_dx);
                Vector characteristic_jump_dy=Matrix_Vector_Product(L,jump_dy);

                for(int s=0;s<4;++s){
                    sigma_squared[s]
                        +=characteristic_jump[s]*characteristic_jump[s]
                            /number_of_faces;
                    sigma_squared[s]
                        +=0.5*Ele[id].h_K*Ele[id].h_K
                            *(characteristic_jump_dx[s]*characteristic_jump_dx[s]
                            +characteristic_jump_dy[s]*characteristic_jump_dy[s])
                            /number_of_faces;
                }
            }
        }

        damping_coef[id]=0.0;
        for(int s=0;s<4;++s){
            damping_coef[id]=max(
                damping_coef[id],sqrt(max(0.0,sigma_squared[s]))
            );
        }
        sigma0=max(sigma0,damping_coef[id]);
    }
    return;
}

inline Vector RHS(
    int id,
    int poin,
    const std::vector<vector<State_of_Fluid>>& v
){

    const Matrix Jaco=Ele[id].Jaco;
    const double det_Jaco=Ele[id].det_Jaco;

    Vector sum{};

    for(int i=0;i<N_GLP;++i){
        double coefficient_x=
            (Jaco[1][1]*D_1[poin][i]-Jaco[1][0]*D_2[poin][i])
            /det_Jaco;
        double coefficient_y=
            (-Jaco[0][1]*D_1[poin][i]+Jaco[0][0]*D_2[poin][i])
            /det_Jaco;
        sum-=2.0*(
            coefficient_x*EC_FLUX_1(v[id][poin],v[id][i])
            +coefficient_y*EC_FLUX_2(v[id][poin],v[id][i])
        );
    }

    Vector tmp{};

    double t1 = B_1[poin];
    double t2 = B_2[poin];

    const double xr = Jaco[0][0];
    const double xs = Jaco[0][1];
    const double yr = Jaco[1][0];
    const double ys = Jaco[1][1];

    double Nx = ( ys * t1 - yr * t2);
    double Ny = (-xs * t1 + xr * t2);

    double Len=sqrt(Nx*Nx+Ny*Ny);

    tmp+=Nx*Physical_Flux_1(v[id][poin])
        +Ny*Physical_Flux_2(v[id][poin]);

    if(is_egde[poin]!=-1){
        tmp -=Len*ES_FLUX(id,poin,v);
    }

    tmp=tmp/(Gauss_2D[poin].second*det_Jaco);

    sum+=tmp;

    if(USE_DAMPING){
        sum-=damping_coef[id]*(v[id][poin].get_cal_data()-u_ave[id]);
    }

    return sum;
}

inline void SSP(){

    now_time=0.0;

    while(now_time<T_end-1e-14){
        stage_time=now_time;
        Get_Sigma(u_h);
        tau=compute_tau();
        if(now_time+tau>T_end)tau=T_end-now_time;

        // ---------------- stage 1 ----------------
        stage_time=now_time;
        for(int i=1;i<=Nele;++i){
            for(int j=0;j<N_GLP;++j){
                u_1[i][j].Set_Conserved(
                    u_h[i][j].get_cal_data()+tau*RHS(i,j,u_h)
                );
            }
        }

        // ---------------- stage 2 ----------------
        stage_time=now_time+tau;
        Get_Sigma(u_1);
        for(int i=1;i<=Nele;++i){
            for(int j=0;j<N_GLP;++j){
                u_2[i][j].Set_Conserved(
                    3.0/4.0*u_h[i][j].get_cal_data()
                    +1.0/4.0*(u_1[i][j].get_cal_data()+tau*RHS(i,j,u_1))
                );
            }
        }

        // ---------------- stage 3 ----------------
        stage_time=now_time+0.5*tau;
        Get_Sigma(u_2);
        for(int i=1;i<=Nele;++i){
            for(int j=0;j<N_GLP;++j){
                u_3[i][j].Set_Conserved(
                    1.0/3.0*u_h[i][j].get_cal_data()
                    +2.0/3.0*(u_2[i][j].get_cal_data()+tau*RHS(i,j,u_2))
                );
            }
        }

        for(int i=1;i<=Nele;++i){
            for(int j=0;j<N_GLP;++j){
                u_h[i][j] = u_3[i][j];
            }
        }
        now_time+=tau;
    }

    return;
}

inline void Init_Data_Cal(){
    for(int id=1;id<=Nele;++id)
    {
        Ele[id].Jaco.resize(2, vector<double>(2));

        PairD v0=Nod[Ele[id].Vertical[0]].Coordinate;
        PairD v1=Nod[Ele[id].Vertical[1]].Coordinate;
        PairD v2=Nod[Ele[id].Vertical[2]].Coordinate;

        for(int j=0;j<N_GLP;++j)
        {
            double x_ref=Gauss_2D[j].first.first;
            double y_ref=Gauss_2D[j].first.second;

            double x_phy=(v2.first-v0.first)*x_ref+(v1.first-v0.first)*y_ref+v0.first;
            double y_phy=(v2.second-v0.second)*x_ref+(v1.second-v0.second)*y_ref+v0.second;
            u_h[id][j]=Init_Solution({x_phy, y_phy});
        }

        const double x0 = v0.first;
        const double y0 = v0.second;

        const double x1 = v1.first;
        const double y1 = v1.second;

        const double x2 = v2.first;
        const double y2 = v2.second;

        const double xr = x2 - x0;
        const double xs = x1 - x0;
        const double yr = y2 - y0;
        const double ys = y1 - y0;

        Ele[id].Jaco[0][0] = xr;  // x_r
        Ele[id].Jaco[0][1] = xs;  // x_s
        Ele[id].Jaco[1][0] = yr;  // y_r
        Ele[id].Jaco[1][1] = ys;  // y_s

        Ele[id].det_Jaco = xr * ys - xs * yr;
    }
    return;
}

inline void Input_Parameter(){
    cout<<"Please Input k and T_end: ";
    cin>>k_poly>>T_end;

    if(k_poly<1||k_poly>4){
        cout<<"This program supports k=1,2,3,4.\n";
        exit(1);
    }

    if(!(T_end>0.0)){
        cout<<"T_end must be positive.\n";
        exit(1);
    }

    return;
}

inline void Init_Reference_Data(){
    Get_Quadrature_Rules_On_Reference_Triangles(k_poly);

    int control_max=max(N_GLP,(k_poly+1)*(k_poly+2)/2);

    D_1.assign(control_max, std::vector<double>(control_max, 0.0));
    D_2.assign(control_max, std::vector<double>(control_max, 0.0));

    B_1.assign(control_max, 0.0);
    B_2.assign(control_max, 0.0);

    for (int i = 0; i < 37; ++i) {
        is_egde[i] = -1;
        from_2d_to_1d[i] = -1;
    }
    for (int e = 0; e < 3; ++e)
    for (int i = 0; i < 37; ++i)
    from_1d_to_2d[e][i] = -1;

    Get_Boundary_Matrix(k_poly);
    Build_Diff_Matrix();
    Build_Damping_Projection_Operators();
    return;
}

inline void Init_Mesh_Data(){
    int control_max=max(N_GLP,(k_poly+1)*(k_poly+2)/2);

    u_h.assign(Nele+1,std::vector<State_of_Fluid>(control_max));
    u_1.assign(Nele+1,std::vector<State_of_Fluid>(control_max));
    u_2.assign(Nele+1,std::vector<State_of_Fluid>(control_max));
    u_3.assign(Nele+1,std::vector<State_of_Fluid>(control_max));

    for (int id = 1; id <= Nele; ++id) {
        for (int e = 0; e < 3; ++e) {
            Ele[id].Neighbor[e] = -1;
        }
        Ele[id].det_Jaco = 0.0;
        Ele[id].h_K = 0.0;
        Ele[id].Jaco = Matrix(2, std::vector<double>(2, 0.0));
    }

    u_ave.assign(Nele+1,Vector{});
    damping_coef.assign(Nele+1,0.0);
    projected_vertex_value.assign(Nele+1,std::array<Vector,3>{});
    projected_vertex_dx.assign(Nele+1,std::array<Vector,3>{});
    projected_vertex_dy.assign(Nele+1,std::array<Vector,3>{});
    sigma0=0.0;
    return;
}

/*********************************************************************************************************************************
Part 7: Output
**********************************************************************************************************************************/

inline bool Element_Is_Inside_Box(
    int id,
    double xmin,
    double xmax,
    double ymin,
    double ymax
){
    const double tolerance=1e-12*DOMAIN_LENGTH;

    for(int vertex=0;vertex<3;++vertex){
        const PairD p=Nod[Ele[id].Vertical[vertex]].Coordinate;

        if(p.first<xmin-tolerance||p.first>xmax+tolerance
            ||p.second<ymin-tolerance||p.second>ymax+tolerance){
            return false;
        }
    }

    return true;
}

inline std::array<double,3> Discrete_h_Errors(
    const std::vector<std::vector<State_of_Fluid>>& u,
    double t_exact
){
    std::array<double,3> error_square={0.0,0.0,0.0};

    for(int id=1;id<=Nele;++id){
        const bool inside_5_15=
            Element_Is_Inside_Box(id,5.0,15.0,5.0,15.0);
        const bool inside_2p5_17p5=
            Element_Is_Inside_Box(id,2.5,17.5,2.5,17.5);

        for(int j=0;j<N_GLP;++j){
            PairD p=Get_Phy_Coordinate(id,j);
            State_of_Fluid exact=Exact_Solution(p,t_exact);
            Vector error=exact.get_cal_data()-u[id][j].get_cal_data();

            double local_error_square=0.0;
            for(int s=0;s<4;++s){
                local_error_square+=error[s]*error[s];
            }

            local_error_square*=
                Ele[id].det_Jaco*Gauss_2D[j].second;

            error_square[0]+=local_error_square;
            if(inside_5_15)error_square[1]+=local_error_square;
            if(inside_2p5_17p5)error_square[2]+=local_error_square;
        }
    }

    for(double& value:error_square)value=sqrt(value);
    return error_square;
}

/*********************************************************************************************************************************
Part 8: Main
**********************************************************************************************************************************/

int main(){
    const vector<pair<string,double> >mesh_data={
        {"square_h0p25.msh",      0.25},
        {"square_h0p125.msh",     0.125},
        {"square_h0p0625.msh",    0.0625},
        {"square_h0p03125.msh",   0.03125},
        {"square_h0p015625.msh",  0.015625},
        {"square_h0p0078125.msh", 0.0078125},
        {"square_h0p00390625.msh",0.00390625},
        {"square_h0p001953125.msh",0.001953125}
    };

    Input_Parameter();
    Init_Reference_Data();

    cout
        <<"20/h  full-domain error  order  "
        <<"[5,15]^2 error  order  "
        <<"[2.5,17.5]^2 error  order\n";

    std::array<double,3> previous_errors={0.0,0.0,0.0};
    double previous_h=0.0;

    for(size_t round=0;round<mesh_data.size();++round){
        string mesh_file=MESH_DIRECTORY+mesh_data[round].first;
        double nominal_h=mesh_data[round].second;

        if(!Read_MSH2(mesh_file))return 1;

        // Reuse the meshes of [0,1]^2 by mapping them onto [0,20]^2.
        Scale_Mesh(DOMAIN_LENGTH);
        FixOrientation();
        Init_Mesh_Data();
        Build_Adjacency();
        Init_Data_Cal();
        Get_Mesh_Size();

        SSP();

        std::array<double,3> errors=Discrete_h_Errors(u_h,T_end);
        double physical_h=DOMAIN_LENGTH*nominal_h;

        cout<<setw(4)<<llround(DOMAIN_LENGTH/physical_h);

        for(int region=0;region<3;++region){
            cout<<"  "<<scientific<<setprecision(6)<<errors[region]<<"  ";

            if(round==0)cout<<"  -- ";
            else{
                double order=
                    log(previous_errors[region]/errors[region])
                    /log(previous_h/physical_h);
                cout<<fixed<<setprecision(3)<<setw(5)<<order;
            }
        }
        cout<<"\n";

        previous_errors=errors;
        previous_h=physical_h;
    }

    return 0;
}
