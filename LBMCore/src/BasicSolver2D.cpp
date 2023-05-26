#include "BasicSolver2D.hpp"
#include "Log.hpp"
#include <iostream>
#include <sstream>
#include <fstream>

Engine::BasicSolver2D::BasicSolver2D(int Nx, int Ny, double T, int numspec):
    mNx(Nx), 
    mNy(Ny), 
    temperature(T), 
    number_of_species(numspec),
    rhomulticomponent(numspec, std::vector<std::vector<double>>(Nx + 2, std::vector<double>(Ny + 2))),
    critical_temperatures(std::vector<double>(numspec)),
    critical_rho(std::vector<double>(numspec)),
    critical_pressures(std::vector<double>(numspec)),
    molarmass(std::vector<double>(numspec)),
    omega(std::vector<double>(numspec)),
    gamma(std::vector<double>(numspec)),
    ux_spec(numspec, std::vector<std::vector<double>>(Nx + 2, std::vector<double>(Ny + 2))),
    uy_spec(numspec, std::vector<std::vector<double>>(Nx + 2, std::vector<double>(Ny + 2))),
    gamma_multiply_rho(std::vector<std::vector<double>>(Nx + 2, std::vector<double>(Ny + 2))),
    rho(std::vector<std::vector<double>>(Nx + 2, std::vector<double>(Ny + 2))),
    pressure(std::vector<std::vector<double>>(Nx + 2, std::vector<double>(Ny + 2))),
    ux(std::vector<std::vector<double>>(Nx + 2, std::vector<double>(Ny + 2))),
    uy(std::vector<std::vector<double>>(Nx + 2, std::vector<double>(Ny + 2))),
    dux_force(std::vector<std::vector<double>>(Nx + 2, std::vector<double>(Ny + 2))),
    duy_force(std::vector<std::vector<double>>(Nx + 2, std::vector<double>(Ny + 2))),
    dux_force_spec(numspec, std::vector<std::vector<double>>(Nx + 2, std::vector<double>(Ny + 2))),
    duy_force_spec(numspec, std::vector<std::vector<double>>(Nx + 2, std::vector<double>(Ny + 2))),
    effrho(std::vector<std::vector<double>>(Nx + 2, std::vector<double>(Ny + 2))),
    sqr_effrho(std::vector<std::vector<double>>(Nx + 2, std::vector<double>(Ny + 2))),
    fmulticomponent(numspec, std::vector<std::vector<std::vector<double>>>
        (9, std::vector<std::vector<double>>(Nx + 2, std::vector<double>(Ny + 2))))
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
    set_border_conditions();
    movement_step();
    calculate_moments();
    calculate_pressure();
    set_Phi();
    calculate_force();
    collision_step();
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

void Engine::BasicSolver2D::check_rho()
{
    double full_rho = 0;
    double full_rho_nC5 = 0;
    double full_rho_C1 = 0;
    std::vector<double> rhol(number_of_species, 0.);
    std::vector<double> rhov(number_of_species, 0.);
    std::vector<double> ful_rhol(number_of_species, 0.);
    std::vector<double> ful_rhov(number_of_species, 0.);
    int S_rhol = 0;
    int S_rhov = 0;
    double max_rho_temp = 0;
    double min_rho_temp = 1000;
    int max_rho_x = 0;
    int max_rho_y = 0;
    int min_rho_x = 0;
    int min_rho_y = 0;
    double r_liquid = 0;
    for (int numspec = 0; numspec < number_of_species; numspec++) {
        for (int i = 1; i <= mNx; i++)
        {
            for (int j = 1; j <= mNy; j++)
            {
                if (numspec == 0) {
                    full_rho += rho[i][j];
                    full_rho_nC5 += rhomulticomponent[numspec][i][j];
                    full_rho_C1 += rhomulticomponent[1][i][j];
                    if (rho[i][j] < min_rho_temp)
                    {
                        min_rho = rho[i][j];
                        min_rho_temp = rho[i][j];
                    }
                    if (rho[i][j] > max_rho_temp)
                    {
                        max_rho = rho[i][j];
                        max_rho_temp = rho[i][j];
                    }
                    if (max_rho_temp == rho[i][j])
                    {
                        max_rho_x = i;
                        max_rho_y = j;
                    }
                    if (min_rho_temp == rho[i][j])
                    {
                        min_rho_x = i;
                        min_rho_y = j;
                    }
                }
            }
        }
    }
    /*for (int i = 1; i <= mNx; i++)
    {
        for (int j = 1; j <= mNy; j++)
        {
            if (rho[i][j] > (max_rho + min_rho) / 2)
            {
                ful_rhol[0] += rhomulticomponent[0][i][j];
                ful_rhol[1] += rhomulticomponent[1][i][j];
                S_rhol++;
            }
            if (rho[i][j] <= (max_rho + min_rho) / 2)
            {
                ful_rhov[0] += rhomulticomponent[0][i][j];
                ful_rhov[1] += rhomulticomponent[1][i][j];
                S_rhov++;
            }
        }
    }*/
   /* mol_fraction[0] = (ful_rhov[0] * S_rhov / molarmass[0] + ful_rhov[1] * S_rhov / molarmass[1]) /
        (ful_rhov[0] * S_rhov / molarmass[0] + ful_rhov[1] * S_rhov / molarmass[1] +
            ful_rhol[0] * S_rhol / molarmass[0] + ful_rhol[1] * S_rhol / molarmass[1]);
    mol_fraction[1] = 1. - mol_fraction[0];
    fraction[0] = full_rho_nC5 / molarmass[0] / (full_rho_nC5 / molarmass[0] + full_rho_C1 / molarmass[1]);
    fraction[1] = 1 - fraction[0];
    fraction_vap[0] = ful_rhov[0] * S_rhov / molarmass[0] /
        (ful_rhov[0] * S_rhov / molarmass[0] + ful_rhov[1] * S_rhov / molarmass[1]);
    fraction_vap[1] = 1. - fraction_vap[0];
    fraction_liq[0] = ful_rhol[0] * S_rhol / molarmass[0] /
        (ful_rhol[0] * S_rhol / molarmass[0] + ful_rhol[1] * S_rhol / molarmass[1]);
    fraction_liq[1] = 1. - fraction_liq[0];*/
    std::cout << "full rho = " << full_rho << std::endl;
    std::cout << "max rho = " << max_rho << "\t" << "min rho = " << min_rho << std::endl;
    std::cout << "max rho(x) = " << max_rho_x << "\t" << "max rho(y) = " << max_rho_y << std::endl;
    std::cout << "min rho(x) = " << min_rho_x << "\t" << "min rho(y) = " << min_rho_y << std::endl;
}
