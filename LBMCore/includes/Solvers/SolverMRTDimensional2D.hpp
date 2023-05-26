#pragma once
#include "BasicSolver2D.hpp"

namespace Engine {
	class SolverMRTDimensional2D : public Engine::BasicSolver2D{
	public:
		SolverMRTDimensional2D() = default;
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
		void eq_func_dimensional_moment(double rho, double ux, double uy, double* moment_eq);
		double h = 1.0e-6;
		double delta_t = 1.0e-9;
		double wettability1 = 1.0;
		double b0 = 0.07780669;
		double a0 = 0.4572793;
		//std::vector<double> s_k = { 0, 1.63, 1.14, 1., 1.92, 1., 1.92, 1 / tau, 1 / tau };
		std::vector<double> s_k = { 0, 1 / tau, 1 / tau, 1., 1 / tau, 1., 1 / tau, 1 / tau, 1 / tau };

		std::vector<double> s_i;
		std::vector<std::vector<double>> k_ij;

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