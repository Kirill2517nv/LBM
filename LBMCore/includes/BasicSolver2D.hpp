#pragma once
#include <vector>

namespace Engine {
	class BasicSolver2D {

	public:
		BasicSolver2D(int Nx, int Ny, double T, int numspec);
		BasicSolver2D(const BasicSolver2D& other) = delete;
		BasicSolver2D(BasicSolver2D&& other) = delete;
		virtual void SaveVTKFile(int tStep);
		virtual ~BasicSolver2D();

		virtual void LBM_Step();
		virtual void set_initial_conditions() = 0;
		virtual void set_border_conditions() = 0;
		virtual void movement_step() = 0;
		virtual void calculate_moments() = 0;
		virtual void calculate_pressure() = 0;
		virtual void set_Phi() = 0;
		virtual void calculate_force() = 0;
		virtual void collision_step() = 0;

		virtual void check_rho();

		const int get_Nx() { return mNx; }
		const int get_Ny() { return mNy; }
		const double get_max_rho() { return max_rho; }
		const double get_min_rho() { return min_rho; }
		std::vector<std::vector<std::vector<double>>> get_rhomulticomponent() { return rhomulticomponent; }

	protected:
		virtual void eq_func(double rho, double ux, double uy, double* f_eq);
		double k = 1.;
		double R = 8.31446; // [J/(mol*K)]
		double max_rho = 0;
		double min_rho = 100;
		double A2 = -0.580;  // -0.8080 -0.1381 Coefficient for calculating forces Peng-Robinson
		double tau = 1.0; // relaxation time
		double temperature; // Temperature of fluid [K]
		int mNx; // Number of nodes along the x-axis
		int mNy; // Number of nodes along the y-axis
		int number_of_species; // Number of species on fluid
		std::vector<double> gamma; // Coefficient ??????????
		std::vector<double> omega; // Acentric factor of each component
		std::vector<double> critical_temperatures; // Critical temperatures of each component [K]
		std::vector<double> critical_pressures; // Critical pressures of each component [MPa]
		std::vector<double> critical_rho; // Critical densities of each component [kg/m^3]
		std::vector<double> molarmass; // Molar masses of each component [kg/mol]
		std::vector<std::vector<double>> gamma_multiply_rho; // need for calculate multicomponent forces
		std::vector<std::vector<double>> rho; // Total density of fluid [kg/m^3]
		std::vector<std::vector<double>> effrho; // need for calculate forces
		std::vector<std::vector<double>> sqr_effrho;
		std::vector<std::vector<double>> ux; // Projection of the total fluid velocity on the x-axis at each node
		std::vector<std::vector<double>> uy; // Projection of the total fluid velocity on the y-axis at each node
		std::vector<std::vector<double>> dux_force; // Projection of the total force on the x-axis at each node
		std::vector<std::vector<double>> duy_force; // Projection of the total force on the y-axis at each node
		std::vector<std::vector<std::vector<double>>> ux_spec; // Projection of the component of fluid velocity on the x-axis at each node
		std::vector<std::vector<std::vector<double>>> uy_spec; // Projection of the component of fluid velocity on the y-axis at each node
		std::vector<std::vector<std::vector<double>>> dux_force_spec; // Projection of the component force on the x-axis at each node
		std::vector<std::vector<std::vector<double>>> duy_force_spec; // Projection of the component force on the y-axis at each node
		std::vector<std::vector<std::vector<double>>> rhomulticomponent; // Component density[]
		std::vector<std::vector<std::vector<std::vector<double>>>> fmulticomponent; // Function of distribution
		std::vector<std::vector<double>> pressure; // Pressure
		std::vector<double> Ai; // coefficient in PR EoS
		std::vector<double> Bi; // coefficient in PR EoS
		const int dx[9] = { 0,   1, 0, -1,  0,   1, -1, -1,  1 };
		const int dy[9] = { 0,   0, 1,  0, -1,   1,  1, -1, -1 };


	private:

	};
}