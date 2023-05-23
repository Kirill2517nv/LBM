#include "Solvers/SolverMRTDimensional2D.hpp"

namespace Engine {
	SolverMRTDimensional2D::SolverMRTDimensional2D(int Nx, int Ny, double T, int numspec) :
		BasicSolver2D(Nx, Ny, T, numspec),
		s_i(numspec),
		qx_spec(numspec, std::vector<std::vector<double>>(Nx + 2, std::vector<double>(Ny + 2))),
		qy_spec(numspec, std::vector<std::vector<double>>(Nx + 2, std::vector<double>(Ny + 2))),
		pxx_spec(numspec, std::vector<std::vector<double>>(Nx + 2, std::vector<double>(Ny + 2))),
		pxy_spec(numspec, std::vector<std::vector<double>>(Nx + 2, std::vector<double>(Ny + 2))),
		e_spec(numspec, std::vector<std::vector<double>>(Nx + 2, std::vector<double>(Ny + 2))),
		epsilon_spec(numspec, std::vector<std::vector<double>>(Nx + 2, std::vector<double>(Ny + 2)))
	{
		double f_eq[9];
		for (int numspec = 0; numspec < number_of_species; numspec++)
		{
			for (int i = 1; i <= mNx; i++)
			{
				for (int j = 1; j <= mNy; j++)
				{
					eq_func(rhomulticomponent[numspec][i][j], ux_spec[numspec][i][j], uy_spec[numspec][i][j], f_eq);
					for (int k = 0; k < 9; k++)
					{
						fmulticomponent[numspec][k][i][j] = f_eq[k];
					}
				}
			}
		}
		Bi = { b0 * R * critical_temperatures[0] / (critical_pressures[0] * molarmass[0]),
			 b0 * R * critical_temperatures[1] / (critical_pressures[1] * molarmass[1]) };
		Ai = { R * critical_temperatures[0] / (sqrt(critical_pressures[0]) * molarmass[0]),
			 R * critical_temperatures[1] / (sqrt(critical_pressures[1]) * molarmass[1]) };
		k_ij = { {0., 0.03}, {0.03, 0.} };


		set_initial_conditions();
	}


	void SolverMRTDimensional2D::set_initial_conditions()
	{
		double rho0up = 5; // 2.45
		double rho0down = 1; // 0.048
		double rho1up = 2.5; //0.056
		double rho1down = 0.5; //0.044
		double radius = 5;
		//droplet
	for (int i = 0; i <= mNx + 1; i++)
	{
		for (int j = 0; j <= mNy + 1; j++)
		{
			rhomulticomponent[0][i][j] = rho0down;
			rhomulticomponent[1][i][j] = rho1down;
			if (i <= 3)
			{
				rhomulticomponent[0][i][j] = rho0up;
				rhomulticomponent[1][i][j] = rho1up;
			}
		}
	}
	}


	void SolverMRTDimensional2D::set_border_conditions()
	{
		for (int numspec = 0; numspec < number_of_species; numspec++)
		{
			fmulticomponent[numspec][5][0][0] = fmulticomponent[numspec][7][1][1];
			fmulticomponent[numspec][6][mNx + 1][0] = fmulticomponent[numspec][8][mNx][1];
			fmulticomponent[numspec][7][mNx + 1][mNy + 1] = fmulticomponent[numspec][5][mNx][mNy];
			fmulticomponent[numspec][8][0][mNy + 1] = fmulticomponent[numspec][6][1][mNy];
			for (int i = 1; i < mNy + 1; i++)
			{
				fmulticomponent[numspec][1][0][i] = fmulticomponent[numspec][3][1][i];
				fmulticomponent[numspec][3][mNx + 1][i] = fmulticomponent[numspec][1][mNx][i];
				fmulticomponent[numspec][5][0][i - 1] = fmulticomponent[numspec][7][1][i];
				fmulticomponent[numspec][8][0][i + 1] = fmulticomponent[numspec][6][1][i];
				fmulticomponent[numspec][6][mNx + 1][i - 1] = fmulticomponent[numspec][8][mNx][i];
				fmulticomponent[numspec][7][mNx + 1][i + 1] = fmulticomponent[numspec][5][mNx][i];
			}
			for (int i = 1; i < mNx + 1; i++)
			{
				fmulticomponent[numspec][2][i][0] = fmulticomponent[numspec][4][i][1];
				fmulticomponent[numspec][4][i][mNy + 1] = fmulticomponent[numspec][2][i][mNy];
				fmulticomponent[numspec][5][i - 1][0] = fmulticomponent[numspec][7][i][1];
				fmulticomponent[numspec][6][i + 1][0] = fmulticomponent[numspec][8][i][1];
				fmulticomponent[numspec][7][i + 1][mNy + 1] = fmulticomponent[numspec][5][i][mNy];
				fmulticomponent[numspec][8][i - 1][mNy + 1] = fmulticomponent[numspec][6][i][mNy];
			}
		}
	}


	void SolverMRTDimensional2D::movement_step()
	{
		for (int numspec = 0; numspec < number_of_species; numspec++)
		{
			std::vector<std::vector<double>> f_temp(mNx + 2, std::vector<double>(mNy + 2));
			for (int k = 1; k < 9; k++)
			{
				for (int i = 1; i <= mNx; i++)
				{
					for (int j = 1; j <= mNy; j++)
					{
						int i_adj = i - dx[k];
						int j_adj = j - dy[k];
						f_temp[i][j] = fmulticomponent[numspec][k][i_adj][j_adj];

					}
				}
				fmulticomponent[numspec][k].swap(f_temp);
			}
		}
	}


	void SolverMRTDimensional2D::calculate_moments()
	{
		for (int i = 1; i <= mNx; i++) {
			for (int j = 1; j <= mNy; j++) {
				rho[i][j] = 0;
				gamma_multiply_rho[i][j] = 0;
				ux[i][j] = 0;
				uy[i][j] = 0;
				for (int numspec = 0; numspec < number_of_species; numspec++)
				{
					e_spec[numspec][i][j] = 0;
					epsilon_spec[numspec][i][j] = 0;
					qx_spec[numspec][i][j] = 0;
					qy_spec[numspec][i][j] = 0;
					pxx_spec[numspec][i][j] = 0;
					pxy_spec[numspec][i][j] = 0;
					rhomulticomponent[numspec][i][j] = 0;
					ux_spec[numspec][i][j] = 0;
					uy_spec[numspec][i][j] = 0;
					for (int k = 0; k < 9; k++)
					{
						rhomulticomponent[numspec][i][j] += fmulticomponent[numspec][k][i][j];
						e_spec[numspec][i][j] += M[1][k] * fmulticomponent[numspec][k][i][j];
						epsilon_spec[numspec][i][j] += M[2][k] * fmulticomponent[numspec][k][i][j];
						ux_spec[numspec][i][j] += M[3][k] * fmulticomponent[numspec][k][i][j];
						qx_spec[numspec][i][j] += M[4][k] * fmulticomponent[numspec][k][i][j];
						uy_spec[numspec][i][j] += M[5][k] * fmulticomponent[numspec][k][i][j];
						qy_spec[numspec][i][j] += M[6][k] * fmulticomponent[numspec][k][i][j];
						pxx_spec[numspec][i][j] += M[7][k] * fmulticomponent[numspec][k][i][j];
						pxy_spec[numspec][i][j] += M[8][k] * fmulticomponent[numspec][k][i][j];
					}
					gamma_multiply_rho[i][j] += gamma[numspec] * rhomulticomponent[numspec][i][j];
					rho[i][j] += rhomulticomponent[numspec][i][j];
					ux_spec[numspec][i][j] *= h / delta_t;
					uy_spec[numspec][i][j] *= h / delta_t;
					ux[i][j] += ux_spec[numspec][i][j];
					uy[i][j] += uy_spec[numspec][i][j];
					ux_spec[numspec][i][j] /= rhomulticomponent[numspec][i][j];
					uy_spec[numspec][i][j] /= rhomulticomponent[numspec][i][j];

				}
				ux[i][j] /= rho[i][j];
				uy[i][j] /= rho[i][j];
			}
		}
	}


	double a(const double& temperature, double omega) // need for PengRobinson EOS
	{
		double m = omega <= 0.49 ? 0.37464 + 1.54226 * omega - 0.269922 * omega * omega :
			0.379642 + 1.48503 * omega - 0.164423 * omega * omega + 0.016666 * omega * omega * omega; //Classic
		/*double m = omega <= 0.49 ? 0.382144 + 1.476905 * omega - 0.134488 * omega * omega :
			0.379642 + 1.48503 * omega - 0.164423 * omega * omega + 0.016666 * omega * omega * omega;*/  // Kalashnikov
		double a = pow((1 + m * (1 - sqrt(temperature))), 2);
		return a;
	}


	void SolverMRTDimensional2D::calculate_pressure()
	{
		for (int i = 1; i <= mNx; i++) {
			for (int j = 1; j <= mNy; j++)
			{
				double D = 0;
				double B = 0;
				double A = 0;
				double S = 0;
				for (int numspec = 0; numspec < number_of_species; numspec++)
				{
					D += rhomulticomponent[numspec][i][j] / molarmass[numspec];
					B += rhomulticomponent[numspec][i][j] * Bi[numspec];
					S += s_i[numspec] * Bi[numspec] * rhomulticomponent[numspec][i][j];
					for (int numspec2 = 0; numspec2 < number_of_species; numspec2++)
					{
						A += rhomulticomponent[numspec][i][j] * rhomulticomponent[numspec2][i][j] * Ai[numspec] * Ai[numspec2] *
							sqrt(a(temperature / critical_temperatures[numspec], omega[numspec]) *
								a(temperature / critical_temperatures[numspec2], omega[numspec2])
							) * (1 - k_ij[numspec][numspec2]);
					}

				}
				D *= R;
				A *= a0;
				pressure[i][j] = D * temperature / (1 + S - B) - A / ((1 + S) * (1 + S) + 2 * B * (1 + S) - B * B);
			}
		}
	}


	void SolverMRTDimensional2D::set_Phi()
	{
		for (int i = 1; i <= mNx; i++) {
			for (int j = 1; j <= mNy; j++) {
				sqr_effrho[i][j] = rho[i][j] * h * h / (3.0 * delta_t * delta_t) - k * pressure[i][j];
				effrho[i][j] = sqrt(sqr_effrho[i][j]);
			}
		}
		effrho[0][mNy + 1] = wettability1 * effrho[1][mNy];
		effrho[mNx + 1][mNy + 1] = wettability1 * effrho[mNx][mNy];
		effrho[0][0] = wettability1 * effrho[1][1];
		effrho[mNx + 1][0] = wettability1 * effrho[mNx][1];
		sqr_effrho[0][0] = effrho[0][0] * effrho[0][0];
		sqr_effrho[0][mNy + 1] = effrho[0][mNy + 1] * effrho[0][mNy + 1];
		sqr_effrho[mNx + 1][0] = effrho[mNx + 1][0] * effrho[mNx + 1][0];
		sqr_effrho[mNx + 1][mNy + 1] = effrho[mNx + 1][mNy + 1] * effrho[mNx + 1][mNy + 1];
		for (int i = 0; i <= mNx + 1; i++) {
			effrho[i][mNy + 1] = wettability1 * effrho[i][mNy];
			effrho[i][0] = wettability1 * effrho[i][1];
			sqr_effrho[i][mNy + 1] = effrho[i][mNy + 1] * effrho[i][mNy + 1];
			sqr_effrho[i][0] = effrho[i][0] * effrho[i][0];
		}
		for (int i = 1; i <= mNy; i++) {
			effrho[0][i] = effrho[1][i];
			effrho[mNx + 1][i] = effrho[mNx][i];
			sqr_effrho[0][i] = effrho[0][i] * effrho[0][i];
			sqr_effrho[mNx + 1][i] = effrho[mNx + 1][i] * effrho[mNx + 1][i];
		}
	}


	void SolverMRTDimensional2D::calculate_force()
	{
		for (int i = 1; i <= mNx; i++) {
			for (int j = 1; j <= mNy; j++) {
				double force_x = 0;
				force_x = 2.0 / 3 / h * delta_t * ((1 - 2 * A2) * effrho[i][j] * (effrho[i + 1][j] - effrho[i - 1][j] +
					0.25 * (effrho[i + 1][j + 1] + effrho[i + 1][j - 1] - effrho[i - 1][j + 1] - effrho[i - 1][j - 1])) +
					A2 * (sqr_effrho[i + 1][j] - sqr_effrho[i - 1][j] + 0.25 *
						(sqr_effrho[i + 1][j + 1] + sqr_effrho[i + 1][j - 1] - sqr_effrho[i - 1][j + 1] - sqr_effrho[i - 1][j - 1])));
				dux_force[i][j] = force_x / rho[i][j];
				double force_y = 0;
				force_y = 2.0 / 3 / h * delta_t * ((1 - 2 * A2) * effrho[i][j] * (effrho[i][j + 1] - effrho[i][j - 1] +
					0.25 * (effrho[i + 1][j + 1] + effrho[i - 1][j + 1] - effrho[i + 1][j - 1] - effrho[i - 1][j - 1])) +
					A2 * (sqr_effrho[i][j + 1] - sqr_effrho[i][j - 1] + 0.25 *
						(sqr_effrho[i + 1][j + 1] + sqr_effrho[i - 1][j + 1] - sqr_effrho[i + 1][j - 1] - sqr_effrho[i - 1][j - 1])));
				duy_force[i][j] = force_y / rho[i][j];
				for (int numspec = 0; numspec < number_of_species; numspec++) {
					dux_force_spec[numspec][i][j] = force_x * gamma[numspec] / gamma_multiply_rho[i][j];
					duy_force_spec[numspec][i][j] = force_y * gamma[numspec] / gamma_multiply_rho[i][j];
				}

			}
		}
	}


	void SolverMRTDimensional2D::collision_step()
	{
		for (int numspec = 0; numspec < number_of_species; numspec++)
		{
			double f_eq_spec[9], f_eq_spec_n[9];
			for (int i = 1; i <= mNx; i++)
			{
				for (int j = 1; j <= mNy; j++)
				{
					e_spec[numspec][i][j] -= s2 * (e_spec[numspec][i][j]
						- 1. / 4. * alfa2 * rhomulticomponent[numspec][i][j]
						- 1. / 6. * gamma2 * rhomulticomponent[numspec][i][j] * rhomulticomponent[numspec][i][j] *
						(ux_spec[numspec][i][j] * ux_spec[numspec][i][j] + uy_spec[numspec][i][j] * uy_spec[numspec][i][j]));
					epsilon_spec[numspec][i][j] -= s3 * (epsilon_spec[numspec][i][j]
						- 1. / 4. * alfa3 * rhomulticomponent[numspec][i][j]
						- 1. / 6. * gamma4 * rhomulticomponent[numspec][i][j] * rhomulticomponent[numspec][i][j] *
						(ux_spec[numspec][i][j] * ux_spec[numspec][i][j] + uy_spec[numspec][i][j] * uy_spec[numspec][i][j]));
					qx_spec[numspec][i][j] -= s5 * (qx_spec[numspec][i][j]
						- 0.5 * c1 * rhomulticomponent[numspec][i][j] * ux_spec[numspec][i][j]);
					qy_spec[numspec][i][j] -= s7 * (qy_spec[numspec][i][j]
						- 0.5 * c1 * rhomulticomponent[numspec][i][j] * uy_spec[numspec][i][j]);
					pxx_spec[numspec][i][j] -= s8 * (pxx_spec[numspec][i][j]
						- 0.5 * gamma1 * rhomulticomponent[numspec][i][j] * rhomulticomponent[numspec][i][j] *
						(ux_spec[numspec][i][j] * ux_spec[numspec][i][j] - uy_spec[numspec][i][j] * uy_spec[numspec][i][j]));
					pxy_spec[numspec][i][j] -= s9 * (pxy_spec[numspec][i][j]
						- 0.5 * gamma3 * rhomulticomponent[numspec][i][j] * rhomulticomponent[numspec][i][j] *
						ux_spec[numspec][i][j] * uy_spec[numspec][i][j]);
					moments[0] = rhomulticomponent[numspec][i][j];
					moments[1] = e_spec[numspec][i][j];
					moments[2] = epsilon_spec[numspec][i][j];
					moments[3] = rhomulticomponent[numspec][i][j] * ux_spec[numspec][i][j];
					moments[4] = qx_spec[numspec][i][j];
					moments[5] = rhomulticomponent[numspec][i][j] * uy_spec[numspec][i][j];
					moments[6] = qy_spec[numspec][i][j];
					moments[7] = pxx_spec[numspec][i][j];
					moments[8] = pxy_spec[numspec][i][j];
					eq_func(rhomulticomponent[numspec][i][j], ux_spec[numspec][i][j], uy_spec[numspec][i][j], f_eq_spec);
					eq_func(rhomulticomponent[numspec][i][j], ux_spec[numspec][i][j] + dux_force_spec[numspec][i][j],
						uy_spec[numspec][i][j] + duy_force_spec[numspec][i][j], f_eq_spec_n);
					for (int k = 0; k < 9; k++)
					{
						fmulticomponent[numspec][k][i][j] = 0;
						for (int z = 0; z < 9; z++) {
							fmulticomponent[numspec][k][i][j] += M_1[k][z] * moments[z];
						}
						fmulticomponent[numspec][k][i][j] += f_eq_spec_n[k] - f_eq_spec[k];
					}
				}
			}
		}
	}


	void SolverMRTDimensional2D::eq_func(double rho, double ux, double uy, double* f_eq)
	{
		double du2 = 1.0 - 1.5 * (ux * ux * delta_t / h * delta_t / h + uy * uy * delta_t / h * delta_t / h);
		f_eq[0] = 4.0 / 9.0 * rho * du2;
		f_eq[1] = rho / 9.0 * (du2 + ux * delta_t / h * (3.0 + 4.5 * ux * delta_t / h));
		f_eq[2] = rho / 9.0 * (du2 + uy * delta_t / h * (3.0 + 4.5 * uy * delta_t / h));
		f_eq[3] = rho / 9.0 * (du2 - ux * delta_t / h * (3.0 - 4.5 * ux * delta_t / h));
		f_eq[4] = rho / 9.0 * (du2 - uy * delta_t / h * (3.0 - 4.5 * uy * delta_t / h));
		f_eq[5] = rho / 36.0 * (du2 + (ux * delta_t / h + uy * delta_t / h) * (3.0 + 4.5 * (ux * delta_t / h + uy * delta_t / h)));
		f_eq[6] = rho / 36.0 * (du2 + (uy * delta_t / h - ux * delta_t / h) * (3.0 + 4.5 * (uy * delta_t / h - ux * delta_t / h)));
		f_eq[7] = rho / 36.0 * (du2 - (ux * delta_t / h + uy * delta_t / h) * (3.0 - 4.5 * (ux * delta_t / h + uy * delta_t / h)));
		f_eq[8] = rho / 36.0 * (du2 + (ux * delta_t / h - uy * delta_t / h) * (3.0 + 4.5 * (ux * delta_t / h - uy * delta_t / h)));
	}
}