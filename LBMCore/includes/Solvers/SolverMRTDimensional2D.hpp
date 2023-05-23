#pragma once
#include "BasicSolver2D.hpp"

namespace Engine {
	class SolverMRTDimensional2D : public Engine::BasicSolver2D{
	public:
		SolverMRTDimensional2D(int Nx, int Ny, double T, int numspec);
		virtual void SolverMRTDimensional2D::set_initial_conditions() override;
		virtual void SolverMRTDimensional2D::set_border_conditions() override;
		virtual void SolverMRTDimensional2D::movement_step() override;
		virtual void SolverMRTDimensional2D::calculate_moments() override;
		virtual void SolverMRTDimensional2D::calculate_pressure() override;
		virtual void SolverMRTDimensional2D::set_Phi() override;
		virtual void SolverMRTDimensional2D::calculate_force() override;
		virtual void SolverMRTDimensional2D::collision_step() override;
	private:
		virtual void eq_func(double rho, double ux, double uy, double* f_eq) override;
		double moments[9];
		double h = 1.0e-6;
		double delta_t = 1.0e-9;
		double wettability1 = 1.0;
		double b0 = 0.07780669;
		double a0 = 0.4572793;
		double s2 = 1.64;
		double s3 = 1.54;
		double s5 = 1.9;
		double s7 = 1.9;

		//double s2 = 1.0 / tau;
		//double s3 = 1.0 / tau;
		//double s5 = 1.0 / tau;
		//double s7 = 1.0 / tau;
		double s8 = 1.0 / tau;
		double s9 = 1.0 / tau;
		double c1 = -2.;
		double alfa2 = -8.;
		double alfa3 = 4.;
		double gamma1 = 2. / 3.;
		double gamma2 = 18.;
		double gamma3 = 2. / 3.;
		double gamma4 = -18.;
		std::vector<double> s_i;
		std::vector<std::vector<double>> k_ij;
		std::vector<std::vector<std::vector<double>>> qx_spec;
		std::vector<std::vector<std::vector<double>>> qy_spec;
		std::vector<std::vector<std::vector<double>>> pxx_spec;
		std::vector<std::vector<std::vector<double>>> pxy_spec;
		std::vector<std::vector<std::vector<double>>> e_spec;
		std::vector<std::vector<std::vector<double>>> epsilon_spec;
		double M[9][9] =
		{
			{1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0 },
			{-4.0, -1.0, -1.0, -1.0, -1.0,  2.0,  2.0,  2.0,  2.0},
			{4.0, -2.0, -2.0, -2.0, -2.0,  1.0,  1.0,  1.0,  1.0},
			{0.0,  1.0,  0.0, -1.0,  0.0,  1.0, -1.0, -1.0,  1.0},
			{0.0, -2.0,  0.0,  2.0,  0.0,  1.0, -1.0, -1.0,  1.0},
			{0.0,  0.0,  1.0,  0.0, -1.0,  1.0,  1.0, -1.0, -1.0},
			{0.0,  0.0, -2.0,  0.0,  2.0,  1.0,  1.0, -1.0, -1.0},
			{0.0,  1.0, -1.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0},
			{0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0},
		};

		double M_1[9][9] =
		{
			{1.0 / 9.0 , -1.0 / 9.0 ,  1.0 / 9.0 ,    0.0     ,     0.0    ,    0.0    ,     0.0    ,    0.0 ,  0.0 },
			{1.0 / 9.0 , -1.0 / 36.0, -1.0 / 18.0,  1.0 / 6.0 , -1.0 / 6.0 ,    0.0    ,     0.0    ,    0.25,  0.0 },
			{1.0 / 9.0 , -1.0 / 36.0, -1.0 / 18.0,    0.0     ,     0.0    , 1.0 / 6.0 , -1.0 / 6.0 ,   -0.25,  0.0 },
			{1.0 / 9.0 , -1.0 / 36.0, -1.0 / 18.0, -1.0 / 6.0 ,  1.0 / 6.0 ,    0.0    ,     0.0    ,    0.25,  0.0 },
			{1.0 / 9.0 , -1.0 / 36.0, -1.0 / 18.0,    0.0     ,     0.0    , -1.0 / 6.0,  1.0 / 6.0 ,   -0.25,  0.0 },
			{1.0 / 9.0 ,  1.0 / 18.0,  1.0 / 36.0,  1.0 / 6.0 ,  1.0 / 12.0,  1.0 / 6.0,  1.0 / 12.0,    0.0 ,  0.25},
			{1.0 / 9.0 ,  1.0 / 18.0,  1.0 / 36.0, -1.0 / 6.0 , -1.0 / 12.0,  1.0 / 6.0,  1.0 / 12.0,    0.0 , -0.25},
			{1.0 / 9.0 ,  1.0 / 18.0,  1.0 / 36.0, -1.0 / 6.0 , -1.0 / 12.0, -1.0 / 6.0, -1.0 / 12.0,    0.0 ,  0.25},
			{1.0 / 9.0 ,  1.0 / 18.0,  1.0 / 36.0,  1.0 / 6.0 ,  1.0 / 12.0, -1.0 / 6.0, -1.0 / 12.0,    0.0 , -0.25},
		};
	};
}