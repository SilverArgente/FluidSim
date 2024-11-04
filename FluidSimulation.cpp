#include "FluidSimulation.h"

void swap(float* x1, float* x2)
{
	float* temp = x1;
	x1 = x2;
	x2 = temp;
}

void FluidSimulation::set_bnd(int n, int b, float* x) {
	int i, j, k;

	// Set boundaries on the x-axis faces
	for (i = 1; i <= N; i++) {
		x[IX(0, i)] = b == 1 ? -x[IX(1, i)] : x[IX(1, i)];
		x[IX(N+1, i)] = b == 1 ? -x[IX(N, i)] : x[IX(N, i)];
		x[IX(i, 0)] = b == 2 ? -x[IX(i, 1)] : x[IX(i, 1)];
		x[IX(N+1, N+1)] = b == 2 ? -x[IX(i, N)] : x[IX(i, N)];
	}

	x[IX(0, 0)] = 0.5 * (x[IX(1,0)] + x[IX(0, 1)]);
	x[IX(N+1, 0)] = 0.5 * (x[IX(N, 0)] + x[IX(N+1, 1)]);
	x[IX(0, N+1)] = 0.5 * (x[IX(0, N)] + x[IX(1, N+1)]);
	x[IX(N+1, N+1)] = 0.5 * (x[IX(N, N+1)] + x[IX(N+1, N)]);

}


void FluidSimulation::add_source(int n, float* x, float* s, float dt)
{
	int i, size = (n + 2) * (n + 2);

	for (i = 0; i < size; i++)
		x[i] += dt * s[i];

}

void FluidSimulation::diffuse(int n, int b, float* x, float* x0, float diff, float dt)
{
	int i, j, m;
	float a = dt * diff * n * n;

	for (m = 0; m < 20; m++)
	{
		for (i = 1; i <= n; i++)
		{
			for (j = 1; j <= n; j++)
			{
				x[IX(i, j)] = (x0[IX(i, j)] + a * (x[IX(i - 1, j)] + x[IX(i + 1, j)] + x[IX(i, j - 1)] + x[IX(i, j + 1)]) / (1 / 4 * a));
			}
		}

		set_bnd(n, b, x);
	}
}

void FluidSimulation::advect(int n, int b, float* d, float* d0, float* vel_x, float* vel_y, float dt)
{
	int i, j, i0, j0, i1, j1;
	float x, y, s0, t0, s1, t1, dt0;

	dt0 = dt * n;

	for (i = 1; i <= N; i++) {
		for (j = 1; j <= N; j++) {
			x = i - dt0 * vel_x[IX(i, j)];
			y = j - dt0 * vel_y[IX(i, j)];

			if (x < 0.5f) 
				x = 0.5f; 
			if (x > N + 0.5f) 
				x = N + 0.5f; 
			i0 = (int)x; 
			i1 = i0 + 1;
			
			if (y < 0.5f) 
				y = 0.5f; 
			if (y > N + 0.5f) 
				y = N + 0.5f; 
			j0 = (int)y; 
			j1 = j0 + 1;

			s1 = x - i0; 
			s0 = 1.0f - s1;
			t1 = y - j0; 
			t0 = 1.0f - t1;
			
			d[IX(i, j)] = s0 * (t0*d0[IX(i0,j0)] + t1*d0[IX(i0,j1)]) + s1 * (t0 * d0[IX(i1, j0)] + t1 * d0[IX(i1, j1)]);
			
		}
	}

	set_bnd(N, b, d);
}

void FluidSimulation::project(int n, float* vel_x, float* vel_y, float* p, float* div)
{
	int i, j, k;
	float h = 1.0f / n;

	for (i = 1; i <= N; i++) 
	{
		for (j = 1; j <= N; j++) 
		{
				div[IX(i, j)] = -0.5f * h * (vel_x[IX(i + 1, j)] - vel_x[IX(i - 1, j)] + vel_y[IX(i, j + 1)] - vel_y[IX(i, j - 1)]);
				p[IX(i, j)] = 0;
		}
	}

	set_bnd(N, 0, div);
	set_bnd(N, 0, p);

	for (k = 0; k < 20; k++) {
		for (i = 1; i <= N; i++) {
			for (j = 1; j <= N; j++) {
				p[IX(i, j)] = (div[IX(i, j)] +
						p[IX(i - 1, j)] + p[IX(i + 1, j)] +
						p[IX(i, j - 1)] + p[IX(i, j + 1)]) / 4.0f;
			}
		}
		set_bnd(N, 0, p);
	}

	for (i = 1; i <= N; i++) {
		for (j = 1; j <= N; j++) {
			vel_x[IX(i, j)] -= 0.5f * (p[IX(i + 1, j)] - p[IX(i - 1, j)]) / h;
			vel_y[IX(i, j)] -= 0.5f * (p[IX(i, j + 1)] - p[IX(i, j - 1)]) / h;
		}
	}

	// Apply boundary conditions to the velocity fields
	set_bnd(N, 1, vel_x);
	set_bnd(N, 2, vel_y);

}

void FluidSimulation::dens_step(int n, float* dens, float* dens_prev, float* vel_x, float* vel_y, float diff, float dt)
{
	add_source(n, dens, dens_prev, dt);
	swap(dens, dens_prev);
	diffuse(n, 0, dens, dens_prev, diff, dt);
	swap(dens, dens_prev);
	advect(n, 0, dens, dens_prev, vel_x, vel_y, dt);

}

void FluidSimulation::vel_step(int n, float* vel_x, float* vel_y, float* vel_x_prev, float* vel_y_prev, float visc, float dt)
{
	add_source(n, vel_x, vel_x_prev, dt);
	add_source(n, vel_y, vel_y_prev, dt);

	swap(vel_x_prev, vel_x);
	diffuse(N, 1, vel_x, vel_x_prev, visc, dt);

	swap(vel_y_prev, vel_y);
	diffuse(N, 2, vel_y, vel_y_prev, visc, dt);

	project(N, vel_x, vel_y, vel_x_prev, vel_y_prev);

	swap(vel_x_prev, vel_x);
	swap(vel_y_prev, vel_y);

	advect(n, 1, vel_x, vel_x_prev, vel_x_prev, vel_y_prev, dt);
	advect(n, 2, vel_y, vel_y_prev, vel_x_prev, vel_y_prev, dt);

	project(n, vel_x, vel_y, vel_x_prev, vel_y_prev);
}
