#ifndef FLUIDSIM_H
#define FLUIDSIM_H
#define N 20
#define IX(i, j) (i  + (N+2)*j)

class FluidSimulation {
public:
    int n;
    float* vel_x, * vel_y;
    float* vel_x_prev, * vel_y_prev;
    float* dens, * dens_prev;

    FluidSimulation()
    {
        int size = (N + 2) * (N + 2);
        
        n = N;
        dens = new float[size];
        dens_prev = new float[size];

        vel_x = new float[size];
        vel_y = new float[size];
        vel_x_prev = new float[size];
        vel_y_prev = new float[size];    
    }

    void set_bnd(int n, int b, float* x);
    void add_source(int n, float* x, float* s, float dt);
    void diffuse(int n, int b, float* x, float* x0, float diff, float dt);
    void advect(int n, int b, float* d, float* d0, float* vel_x, float* vel_y, float dt);
    void project(int n, float* vel_x, float* vel_y, float* p, float* div);
    void vel_step(int n, float* vel_x, float* vel_y, float* vel_x_prev, float* vel_y_prev, float visc, float dt);
    void dens_step(int n, float* dens, float* dens_prev, float* vel_x, float* vel_y, float diff, float dt);
    

};

#endif // FLUIDSIM_H
