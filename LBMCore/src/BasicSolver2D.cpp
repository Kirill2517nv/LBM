#include "BasicSolver2D.hpp"
#include "Log.hpp"
#include <iostream>
#include <sstream>
#include <fstream>

Engine::BasicSolver2D::BasicSolver2D(int Nx, int Ny, double T, int numspec):
    mNx(Nx), mNy(Ny), temperature(T), number_of_species(numspec)
{
}

void Engine::BasicSolver2D::SaveVTKFile(int tStep)
{
    std::stringstream fname;
    fname << "VTK/adv_";
    if (tStep < 10) fname << "0";
    if (tStep < 100) fname << "0";
    if (tStep < 1000) fname << "0";
    if (tStep < 10000) fname << "0";
    if (tStep < 100000) fname << "0";
    if (tStep < 1000000) fname << "0";
    if (tStep < 10000000) fname << "0";
    fname << tStep << ".vtk";
    std::ofstream vtk_file(fname.str().c_str());
    vtk_file << "# vtk DataFile Version 3.0\n";
    vtk_file << "Immiscible displacement\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET RECTILINEAR_GRID\nDIMENSIONS " << mNx << " " << mNy << " 1\n";
    vtk_file << "X_COORDINATES " << mNx << " double\n";
    for (int i = 1; i <= mNx; i++) vtk_file << i << " ";
    vtk_file << std::endl;
    vtk_file << "Y_COORDINATES " << mNy << " double\n";
    for (int i = 1; i <= mNy; i++) vtk_file << i << " ";
    vtk_file << std::endl;
    vtk_file << "Z_COORDINATES 1 double\n0\n";
    vtk_file << "POINT_DATA " << mNx * mNy << std::endl;
    /*  vtk_file << "SCALARS delta_rho double 1\n";
      vtk_file << "LOOKUP_TABLE default\n";
      for (int j = 1; j <= mNy; j++)
          for (int i = 1; i <= mNx; i++)
              if (mask[i][j] == 0) vtk_file << rho_spec[0][i][j] - rho_spec[1][i][j] << " ";
              else vtk_file << -2.0 << " ";
      vtk_file << endl;*/
    vtk_file << "SCALARS rho double 1\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for (int j = 1; j <= mNy; j++)
        for (int i = 1; i <= mNx; i++) vtk_file << rho[i][j] << " ";
    vtk_file << std::endl;
    for (int numspec = 0; numspec < number_of_species; numspec++)
    {
        vtk_file << "SCALARS rho" << numspec << " double 1\n";
        vtk_file << "LOOKUP_TABLE default\n";
        for (int j = 1; j <= mNy; j++)
            for (int i = 1; i <= mNx; i++) vtk_file << rhomulticomponent[numspec][i][j] << " ";
        vtk_file << std::endl;
    }
    vtk_file << "SCALARS pressure double 1\n";
    vtk_file << "LOOKUP_TABLE default\n";
    for (int j = 1; j <= mNy; j++)
        for (int i = 1; i <= mNx; i++) vtk_file << pressure[i][j] << " ";
    vtk_file << std::endl;
    vtk_file << "VECTORS uflow double\n";
    for (int j = 1; j <= mNy; j++)
        for (int i = 1; i <= mNx; i++) vtk_file << ux[i][j] << "  " << uy[i][j] << "  0.0" << " ";
    vtk_file << std::endl;

    for (int numspec = 0; numspec < number_of_species; numspec++)
    {
        vtk_file << "VECTORS uflow" << numspec << " double\n";
        for (int j = 1; j <= mNy; j++)
            for (int i = 1; i <= mNx; i++) vtk_file << ux_spec[numspec][i][j] << "  " << uy_spec[numspec][i][j] << "  0.0" << " ";
        vtk_file << std::endl;
    }
    vtk_file.close();

    LOG_INFO("File {0} written", fname.str());
}

Engine::BasicSolver2D::~BasicSolver2D()
{
}

void Engine::BasicSolver2D::LBM_Step()
{
	void set_border_conditions();
	void movement_step();
	void calculate_moments();
	void calculate_pressure();
	void set_Phi();
	void calculate_force();
	void collision_step();
}

void Engine::BasicSolver2D::eq_func(double rho, double ux, double uy, double* f_eq)
{
    double du2 = 1.0 - 1.5 * (ux * ux + uy * uy);
    f_eq[0] = 4.0 / 9.0 * rho * du2;
    f_eq[1] = rho / 9.0 * (du2 + ux * (3.0 + 4.5 * ux));
    f_eq[2] = rho / 9.0 * (du2 + uy * (3.0 + 4.5 * uy));
    f_eq[3] = rho / 9.0 * (du2 - ux * (3.0 - 4.5 * ux));
    f_eq[4] = rho / 9.0 * (du2 - uy * (3.0 - 4.5 * uy));
    f_eq[5] = rho / 36.0 * (du2 + (ux + uy) * (3.0 + 4.5 * (ux + uy)));
    f_eq[6] = rho / 36.0 * (du2 + (uy - ux) * (3.0 + 4.5 * (uy - ux)));
    f_eq[7] = rho / 36.0 * (du2 - (ux + uy) * (3.0 - 4.5 * (ux + uy)));
    f_eq[8] = rho / 36.0 * (du2 + (ux - uy) * (3.0 + 4.5 * (ux - uy)));
}
