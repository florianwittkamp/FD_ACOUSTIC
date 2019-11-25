// FD_2D_DX4_DT2 2-D acoustic Finite-Difference modelling
//
// GNU General Public License v3.0
//
// Author: Nicolas Salmieri
//
// Finite-Difference acoustic seismic wave simulation
//
// Discretization of the first-order acoustic wave equation

#include "numpylike.h"
#include "bmp.h"

int main(){
	Np np;

	// Discretization
	auto c_1=20;    // Number of grid points per dominant wavelength
	auto c_2=0.5;   // CFL-Number
	auto nx=200;    // Number of grid points in X
	auto ny=200;    // Number of grid points in Y
	auto T=1;     // Total propagation time

	// Source Signal
	auto f0= 5.0;     // Center frequency Ricker-wavelet
	auto q0= 1.0;     // Maximum amplitude Ricker-Wavelet
	std::vector<int> sp {100,100}; // Source position (in grid points) y,x

	// Velocity and density
	auto modell_v = 3000*np.ones({ny,nx});
	auto rho=2.2*np.ones({ny,nx});

	// Calculate first Lame-Paramter
	auto l=rho * modell_v * modell_v;

	auto cmin=min(modell_v.flatten());  // Lowest P-wave velocity
	auto cmax=max(modell_v.flatten());  // Highest P-wave velocity
	auto fmax=2.0*f0;                     // Maximum frequency
	auto dx=cmin/(fmax*c_1);             // Spatial discretization (in m)
	auto dy=dx;                         // Spatial discretization (in m)
	auto dt=dx/(cmax)*c_2;               // Temporal discretization (in s)
	auto lampda_min=cmin/fmax;          // Smallest wavelength

	// Output model parameter:
	print("Model size: x:",dx*nx,"in m, y:",dy*ny,"in m");
	print("Temporal discretization: ",dt," s");
	print("Spatial discretization: ",dx," m");
	print("Number of gridpoints per minimum wavelength: ",lampda_min/dx);

	// Create space and time vector

	auto x=np.arange(0,dx*nx,dx); // Space vector in X
	auto y=np.arange(0,dy*ny,dy); // Space vector in Y
	auto t=np.arange(0,T,dt);     // Time vector
	int nt=np.size(t)[0];         // Number of time steps




	// Source signal - Ricker-wavelet
	auto tau=np.pi*f0*(t-1.5/f0);
	auto q=q0*(1.0-2.0*tau.pow(2.0))*np.exp(0-tau.pow(2.0));

	// Calculation of some coefficients
	//	auto i_dx=1.0/(dx);
	//	auto i_dy=1.0/(dy);
	auto c1=9.0/(8.0*dx);
	auto c2=1.0/(24.0*dx);
	auto c3=9.0/(8.0*dy);
	auto c4=1.0/(24.0*dy);
	//	auto c5=1.0/np.power(dx,3);
	//	auto c6=1.0/np.power(dy,3);
	//	auto c7=1.0/np.power(dx,2);
	//	auto c8=1.0/np.power(dy,2);
	//	auto c9=np.power(dt,3)/24.0;

	print("Starting time stepping...");

	// Init wavefields
	auto vx=np.zeros({ny,nx});
	auto vy=np.zeros({ny,nx});
	auto p=np.zeros({ny,nx});
	for (auto& n : range(0,nt)) {
		// Update velocity
		for (auto& kx : range(5,nx-4)) {
			for (auto& ky : range(5,ny-4)) {
				auto p_x=c1*(p[{ky,kx+1}]-p[{ky,kx}])-c2*(p[{ky,kx+2}]-p[{ky,kx-1}]);
				auto p_y=c3*(p[{ky+1,kx}]-p[{ky,kx}])-c4*(p[{ky+2,kx}]-p[{ky-1,kx}]);

				vx[{ky,kx}]=vx[{ky,kx}]-dt/rho[{ky,kx}]*p_x;
				vy[{ky,kx}]=vy[{ky,kx}]-dt/rho[{ky,kx}]*p_y;
			}
		}

		// Inject source wavelet
		p[sp]=p[sp]+q[n];

		// Update pressure
		for (auto& kx : range(5,nx-4)) {
			for (auto& ky : range(5,ny-4)) {
				auto vx_x=c1*(vx[{ky,kx}]-vx[{ky,kx-1}])-c2*(vx[{ky,kx+1}]-vx[{ky,kx-2}]);
				auto vy_y=c3*(vy[{ky,kx}]-vy[{ky-1,kx}])-c4*(vy[{ky+1,kx}]-vy[{ky-2,kx}]);

				p[{ky,kx}]=p[{ky,kx}]-l[{ky,kx}]*dt*(vx_x+vy_y);
			}
		}
	}
	print("Finished time stepping!");

	auto aze=p.flatten();
	auto mi=min(aze);
	auto ma=max(aze);
	BMP bmp(nx,ny,false);
	for (auto& y : range(1,ny+1)) {
		/*y is inverted in the file format compared to cartesian*/
		for (auto& x : range(0,nx)) {
			uint32_t val=254*(p[{y,x}]-mi)/(ma-mi);
			bmp.set_pixel(x,ny-y,val,val,val,1);
		}
	}
	bmp.write("out.bmp");

	return 0;
}
