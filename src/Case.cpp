#include "Case.hpp"
#include "Enums.hpp"
#include "Communication.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <vector>

namespace filesystem = std::filesystem;

#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredGrid.h>
#include <vtkStructuredGridWriter.h>
#include <vtkTuple.h>

#define IS_DETAIL_LOG (0)
#if IS_DETAIL_LOG
#include <memory>
    int simIter = 10;
#endif

int g_rank = 0; //for debug use
Case::Case(std::string file_name, int argn, char **args) {
    // Read input parameters
    const int MAX_LINE_LENGTH = 1024;
    std::ifstream file(file_name);
    double nu = 0;      /* viscosity   */
    double UI = 0;      /* velocity x-direction */
    double VI = 0;      /* velocity y-direction */
    double PI = 0;      /* pressure */
    double GX = 0;      /* gravitation x-direction */
    double GY = 0;      /* gravitation y-direction */
    double xlength = 0; /* length of the domain x-dir.*/
    double ylength = 0; /* length of the domain y-dir.*/
    double dt = 0;      /* time step */
    int imax = 0;       /* number of cells x-direction*/
    int jmax = 0;       /* number of cells y-direction*/
    double gamma = 0;   /* uppwind differencing factor*/
    double omg = 0;     /* relaxation factor */
    double tau = 0;     /* safety factor for time step*/
    int itermax = 0;    /* max. number of iterations for pressure per time step */
    double eps =0;     /* accuracy bound for pressure*/

    std::string geom_name = _geom_name;          /*geometry file name for the problem*/
    double UIN = 0;                     /* inflow velocity x-direction */
    double VIN = 0;                     /* inflow velocity y-direction */
    int wallnum = 0;                    /* wall num */
    std::map<int, double> wall_vel; /* wall velocities */

    bool energy_on = false;
    std::string energy_eq = "off";
    double TI = 0;                       /* initial temperature */
    double Pr = 0;                       /* Prandtl number (Pr = nu / alpha), here directly given in .dat file*/
    double beta = 0;                     /* thermal expansion coefficient*/
    std::map<int, double> wall_temp; /* wall temperatures */

    int iproc = 1;       /* number of processes per x   */
    int jproc = 1;       /* number of processes per y   */
    bool parallel_on = false;

    if (file.is_open()) {

        std::string var;
        while (!file.eof() && file.good()) {
            file >> var;
            if (var[0] == '#') { /* ignore comment line*/
                file.ignore(MAX_LINE_LENGTH, '\n');
            } else {
                if (var == "iproc") file >> iproc; // TODO: check invalid user input e.g iproc not integer
                if (var == "jproc") file >> jproc;
                if (var == "xlength") file >> xlength;
                if (var == "ylength") file >> ylength;
                if (var == "nu") file >> nu;
                if (var == "t_end") file >> _t_end;
                if (var == "dt") file >> dt;
                if (var == "omg") file >> omg;
                if (var == "eps") file >> eps;
                if (var == "tau") file >> tau;
                if (var == "gamma") file >> gamma;
                if (var == "dt_value") file >> _output_freq;
                if (var == "UI") file >> UI;
                if (var == "VI") file >> VI;
                if (var == "GX") file >> GX;
                if (var == "GY") file >> GY;
                if (var == "PI") file >> PI;
                if (var == "itermax") file >> itermax;
                if (var == "imax") file >> imax;
                if (var == "jmax") file >> jmax;
                if (var == "geo_file") file >> geom_name;
                if (var == "UIN") file >> UIN;
                if (var == "VIN") file >> VIN;
                if (var == "energy_eq") file >> energy_eq;
                if (var == "TI") file >> TI;
                if (var == "Pr") file >> Pr;
                if (var == "alpha") {
                    file >> Pr;
                    Pr = nu / Pr;
                } // read in alpha, convert to Pr
                if (var == "beta") file >> beta;
                if (var == "num_of_walls" || var == "num_walls") file >> wallnum;
                for (int i = 0; i < wallnum; ++i) {
                    int wallIdx = i + 3;
                    std::string str = "wall_vel_" + std::to_string(wallIdx);
                    std::string temp_str = "wall_temp_" + std::to_string(wallIdx);
                    if (var == str) {
                        double tmpVelocity = 0.;
                        file >> tmpVelocity;
                        wall_vel.insert(std::pair<int, double>(wallIdx, tmpVelocity));
                    }
                    if (var == temp_str) {
                        double tmpTemperature = 0.;
                        file >> tmpTemperature;
                        wall_temp.insert(std::pair<int, double>(wallIdx, tmpTemperature));
                    }
                }
            }
        }
    }
    file.close();
    //MPI_initialize or communitate::initialize
    // Parameter safety check
    const double zeroEpilon = 1e-05;
    if (Pr < zeroEpilon) {
        char buffer[1024];
        snprintf(buffer, 1024, "Pr number = %f invalid.Reset to 1.0\n", Pr);
        std::cerr << buffer;
        Pr = 1.0;
    }

    _geom_name = geom_name;
    std::cout << "geom_name = " << geom_name << "\n";
    if (_geom_name.compare("NONE") == 0) {
        wall_vel.insert(std::pair<int, double>(LidDrivenCavity::moving_wall_id, LidDrivenCavity::wall_velocity));
    }

    // Check if need energy equation
    if (energy_eq.compare("on") == 0) {
        energy_on = true;
    }
    _energy_On = energy_on;
    _rank = 0; // init to 0;
    _left_neighbour_rank = -1;
    _right_neighbour_rank = -1;
    _top_neighbour_rank = -1;
    _bottom_neighbour_rank = -1;

    // Check if need parallel
    if (iproc != 1 || jproc != 1) {
        parallel_on = true;
    }
     _parallel_On = parallel_on;
     _jproc = jproc;
     _iproc = iproc;

    // Set file names for geometry file and output directory
    set_file_names(file_name);

    if (_parallel_On)
    {
        int rank = 0;
        int num_proc = 0;
        Communication::init_parallel(argn,args,rank,num_proc);
        _rank = rank; 
        g_rank= rank;
        assert(num_proc == _jproc * _iproc);
    }

    // Build up the domain
    Domain domain;
    domain.dx = xlength / (double)imax;
    domain.dy = ylength / (double)jmax;
    domain.domain_size_x = imax;
    domain.domain_size_y = jmax;

    build_domain(domain, imax, jmax);

    _grid = Grid(_geom_name, domain);
    _field = Fields(_grid, nu, dt, tau, _grid.domain().size_x, _grid.domain().size_y, UI, VI, PI, GX, GY, energy_on, TI,
                    Pr, beta);

    _discretization = Discretization(domain.dx, domain.dy, gamma);
    _pressure_solver = std::make_unique<SOR>(omg);
    _max_iter = itermax;
    _tolerance = eps;

    // Construct boundaries
    if (_geom_name.compare("NONE") == 0) { // LidDrivenCavity
        if (not _grid.moving_wall_cells().empty()) {
            _boundaries.push_back(
                std::make_unique<MovingWallBoundary>(_grid.moving_wall_cells(), LidDrivenCavity::wall_velocity));
        }
        if (not _grid.fixed_wall_cells().empty()) {
            _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.fixed_wall_cells()));
        }
    } else {

        if (not _grid.fixed_wall_cells().empty()) {
            if (energy_on)
                _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.fixed_wall_cells(), wall_temp));
            else
                _boundaries.push_back(std::make_unique<FixedWallBoundary>(_grid.fixed_wall_cells()));
        }

        if (not _grid.inflow_cells().empty()) {
            _boundaries.push_back(std::make_unique<InFlowBoundary>(_grid.inflow_cells(), UIN, VIN));
        }

        if (not _grid.outflow_cells().empty()) {
            _boundaries.push_back(std::make_unique<OutFlowBoundary>(_grid.outflow_cells()));
        }
    }
}

Case::~Case()
{
    if (_parallel_On)
    {
        Communication::finalize();
    }
}
   
void Case::set_file_names(std::string file_name) {
    std::string temp_dir;
    bool case_name_flag = true;
    bool prefix_flag = false;

    for (int i = file_name.size() - 1; i > -1; --i) {
        if (file_name[i] == '/') {
            case_name_flag = false;
            prefix_flag = true;
        }
        if (case_name_flag) {
            _case_name.push_back(file_name[i]);
        }
        if (prefix_flag) {
            _prefix.push_back(file_name[i]);
        }
    }

    for (int i = file_name.size() - _case_name.size() - 1; i > -1; --i) {
        temp_dir.push_back(file_name[i]);
    }

    std::reverse(_case_name.begin(), _case_name.end());
    std::reverse(_prefix.begin(), _prefix.end());
    std::reverse(temp_dir.begin(), temp_dir.end());

    _case_name.erase(_case_name.size() - 4);
    _dict_name = temp_dir;
    _dict_name.append(_case_name);
    _dict_name.append("_Output");

    if (_geom_name.compare("NONE") != 0) {
        _geom_name = _prefix + _geom_name;
    }

    // Create output directory
    filesystem::path folder(_dict_name);
    try {
        filesystem::create_directory(folder);
    } catch (const std::exception &e) {
        std::cerr << "Output directory could not be created." << std::endl;
        std::cerr << "Make sure that you have write permissions to the "
                     "corresponding location"
                  << std::endl;
    }
}

/**
 * This function is the main simulation loop. In the simulation loop, following steps are required
 * - Calculate and apply boundary conditions for all the boundaries in _boundaries container
 *   using apply() member function of Boundary class
 * - Calculate fluxes (F and G) using calculate_fluxes() member function of Fields class.
 *   Flux consists of diffusion and convection part, which are located in Discretization class
 * - Calculate right-hand-side of PPE using calculate_rs() member function of Fields class
 * - Iterate the pressure poisson equation until the residual becomes smaller than the desired tolerance
 *   or the maximum number of the iterations are performed using solve() member function of PressureSolver class
 * - Calculate the velocities u and v using calculate_velocities() member function of Fields class
 * - Calculat the maximal timestep size for the next iteration using calculate_dt() member function of Fields class
 * - Write vtk files using output_vtk() function
 *
 * Please note that some classes such as PressureSolver, Boundary are abstract classes which means they only provide the
 * interface. No member functions should be defined in abstract classes. You need to define functions in inherited
 * classes such as MovingWallBoundary class.
 *
 * For information about the classes and functions, you can check the header files.
 */
void Case::simulate() {

    double t = 0.0;
    double dt = _field.dt();
    int timestep = 0;
    double output_counter = 0.0;

#if IS_DETAIL_LOG
    while (timestep < simIter) {
#else
    while (t < _t_end) {
#endif
        /*****
         apply boundary
        ******/
        for (int i = 0; i < _boundaries.size(); ++i) {
            _boundaries[i]->apply(_field);
        }

        /*****
         update field
        ******/
        _field.calculate_temperature(_grid);
        if (_parallel_On){
            Communication::communicate(
                _grid,
                _field.T_matrix(), 
                _left_neighbour_rank,
                _right_neighbour_rank,
                _bottom_neighbour_rank,
                _top_neighbour_rank);
        }
 
        _field.calculate_fluxes(_grid);
        if (_parallel_On){
            Communication::communicate(
                _grid,
                _field.F_matrix(),
                _left_neighbour_rank,
                _right_neighbour_rank,
                _bottom_neighbour_rank,
                _top_neighbour_rank);
            Communication::communicate(
                _grid,
                _field.G_matrix(),
                _left_neighbour_rank,
                _right_neighbour_rank,
                _bottom_neighbour_rank,
                _top_neighbour_rank);
        }

        _field.calculate_rs(_grid);

        // loops here a number of times
        int it = 0;
        double res = _tolerance + 1.0; // init res to be greater than tolerance to enter the loop
        while (res > _tolerance && it++ < _max_iter) {
            res = _pressure_solver->solve(_field, _grid, _boundaries,_parallel_On);
            if (_parallel_On){
                Communication::communicate(
                    _grid,
                    _field.p_matrix(),
                    _left_neighbour_rank,
                    _right_neighbour_rank,
                    _bottom_neighbour_rank,
                    _top_neighbour_rank);
                res = Communication::reduce_sum(_rank,res); //can be more accuracte by considering global res
                int totalGridSize = Communication::reduce_sum(_rank,_grid.fluid_cells().size());
                res = res / totalGridSize;
                res = std::sqrt(res);
            }
        }
        if (it >= _max_iter) {
            if ((_parallel_On && _rank == 0) || !_parallel_On){
             std::cerr << "Pressure Solver fails to converge at timestep " << timestep << "!\n";
            }
        }
 
        _field.calculate_velocities(_grid);
        if (_parallel_On){
                Communication::communicate(
                    _grid,
                    _field.u_matrix(),
                    _left_neighbour_rank,
                    _right_neighbour_rank,
                    _bottom_neighbour_rank,
                    _top_neighbour_rank);
                Communication::communicate(
                    _grid,
                    _field.v_matrix(),
                    _left_neighbour_rank,
                    _right_neighbour_rank,
                    _bottom_neighbour_rank,
                    _top_neighbour_rank);
        }
        /*****
        increment time
        ******/
        timestep++;
        t += dt;
   
        /*****
        intermediate output field
        ******/
        if (t > _output_freq * output_counter) {
        //if (1){
            output_vtk(timestep,_rank);
            output_counter = output_counter + 1;
        }

        /*****
         Compute Courant number
        ******/
        double CourantNum = fabs(_field.u(_grid.imax() / 2, _grid.jmax() / 2)) * dt / _grid.dx();

        /*****
         Console logging
        ******/
        if ((_parallel_On && _rank == 0) || !_parallel_On){
            char buffer[1024];
            snprintf(buffer, 1024, "step = %d, t = %.3f, p.solver res = %.3e, CNum = %.3e\n", timestep, t, res,
                     CourantNum);
            std::cout << buffer;

#if IS_DETAIL_LOG
            std::FILE* detaillog = std::fopen("./detaillog.log", "a");
            std::fprintf(detaillog, "%s",buffer);

            int pointToObserveI =50;
            int pointToObserveJ =_grid.jmax() / 2;
            snprintf(buffer, 1024, "T = %.3e,%.3e,%.3e,%.3e,%.3e\n", 
            _field.T(pointToObserveI,pointToObserveJ), _field.T(pointToObserveI-1, pointToObserveJ), _field.T(pointToObserveI+1, pointToObserveJ), 
            _field.T(pointToObserveI, pointToObserveJ - 1), _field.T(pointToObserveI, pointToObserveJ + 1));
            std::cout << buffer;
            std::fprintf(detaillog, "%s",buffer);

            snprintf(buffer, 1024, "f = %.3e,%.3e,%.3e,%.3e,%.3e\n", 
            _field.f(pointToObserveI,pointToObserveJ), _field.f(pointToObserveI-1, pointToObserveJ), _field.f(pointToObserveI+1, pointToObserveJ), 
            _field.f(pointToObserveI, pointToObserveJ - 1), _field.f(pointToObserveI, pointToObserveJ + 1));
            std::cout << buffer;
            std::fprintf(detaillog, "%s",buffer);

            snprintf(buffer, 1024, "g = %.3e,%.3e,%.3e,%.3e,%.3e\n", 
            _field.g(pointToObserveI,pointToObserveJ), _field.g(pointToObserveI-1, pointToObserveJ), _field.g(pointToObserveI+1, pointToObserveJ), 
            _field.g(pointToObserveI, pointToObserveJ - 1), _field.g(pointToObserveI, pointToObserveJ + 1));
            std::cout << buffer;
            std::fprintf(detaillog, "%s",buffer);

            snprintf(buffer, 1024, "rs = %.3e,%.3e,%.3e,%.3e,%.3e\n", 
            _field.rs(pointToObserveI,pointToObserveJ), _field.rs(pointToObserveI-1, pointToObserveJ), _field.rs(pointToObserveI+1, pointToObserveJ), 
            _field.rs(pointToObserveI, pointToObserveJ - 1), _field.rs(pointToObserveI, pointToObserveJ + 1));
            std::cout << buffer;
            std::fprintf(detaillog, "%s",buffer);

            snprintf(buffer, 1024, "P = %.3e,%.3e,%.3e,%.3e,%.3e\n", 
            _field.p(pointToObserveI,pointToObserveJ), _field.p(pointToObserveI-1, pointToObserveJ), _field.p(pointToObserveI+1, pointToObserveJ), 
            _field.p(pointToObserveI, pointToObserveJ - 1), _field.p(pointToObserveI, pointToObserveJ + 1));
            std::cout << buffer;
            std::fprintf(detaillog,"%s", buffer);
            
            snprintf(buffer, 1024, "u = %.3e,%.3e,%.3e,%.3e,%.3e\n", 
            _field.u(pointToObserveI,pointToObserveJ), _field.u(pointToObserveI-1, pointToObserveJ), _field.u(pointToObserveI+1, pointToObserveJ), 
            _field.u(pointToObserveI, pointToObserveJ - 1), _field.u(pointToObserveI, pointToObserveJ + 1));
            std::cout << buffer;
            std::fprintf(detaillog, "%s",buffer);

            snprintf(buffer, 1024, "v = %.3e,%.3e,%.3e,%.3e,%.3e\n", 
            _field.v(pointToObserveI,pointToObserveJ), _field.v(pointToObserveI-1, pointToObserveJ), _field.v(pointToObserveI+1, pointToObserveJ), 
            _field.v(pointToObserveI, pointToObserveJ - 1), _field.v(pointToObserveI, pointToObserveJ + 1));
            std::cout << buffer;
            std::fprintf(detaillog, "%s",buffer);

            std::fclose(detaillog);
#endif
        }

        /*****
         Compute time step for next iteration
        ******/
        dt = _field.calculate_dt(_grid);
        if (_parallel_On){
            dt = Communication::reduce_min(_rank,dt);
        }
    }

    /*****
     Output the fields in the end
    ******/
    output_vtk(timestep,_rank);
    output_counter++;
}

void Case::output_vtk(int timestep, int my_rank) {
    // Create a new structured grid
    vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();

    // Create grid
    vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

    double dx = _grid.dx();
    double dy = _grid.dy();

    double x = _grid.domain().imin * dx;
    double y = _grid.domain().jmin * dy;

    { y += dy; }
    { x += dx; }

    double z = 0;
    for (int col = 0; col < _grid.domain().size_y + 1; col++) {
        x = _grid.domain().imin * dx;
        { x += dx; }
        for (int row = 0; row < _grid.domain().size_x + 1; row++) {
            points->InsertNextPoint(x, y, z);
            x += dx;
        }
        y += dy;
    }

    // Specify the dimensions of the grid
    structuredGrid->SetDimensions(_grid.domain().size_x + 1, _grid.domain().size_y + 1, 1);
    structuredGrid->SetPoints(points);

    // Pressure Array
    vtkDoubleArray *Pressure = vtkDoubleArray::New();
    Pressure->SetName("pressure");
    Pressure->SetNumberOfComponents(1);

    // Velocity Array
    vtkDoubleArray *Velocity = vtkDoubleArray::New();
    Velocity->SetName("velocity");
    Velocity->SetNumberOfComponents(3);

    // Temperature Array
    vtkDoubleArray *Temperature = vtkDoubleArray::New();
    Temperature->SetName("temperature");
    Temperature->SetNumberOfComponents(1);

    // Geometry Array
    vtkDoubleArray *Geometry = vtkDoubleArray::New();
    Geometry->SetName("geometry");
    Geometry->SetNumberOfComponents(1);

    // Print pressure and temperature from bottom to top
    for (int j = 1; j < _grid.domain().size_y + 1; j++) {
        for (int i = 1; i < _grid.domain().size_x + 1; i++) {
            double pressure = _field.p(i, j);
            Pressure->InsertNextTuple(&pressure);
        }
    }
    if (_field.energy_on()) {
        for (int j = 1; j < _grid.domain().size_y + 1; j++) {
            for (int i = 1; i < _grid.domain().size_x + 1; i++) {
                double temperature = _field.T(i, j);
                Temperature->InsertNextTuple(&temperature);
            }
        }
    }

    // Temp Velocity
    float vel[3];
    vel[2] = 0; // Set z component to 0

    // Print Velocity from bottom to top
    for (int j = 0; j < _grid.domain().size_y + 1; j++) {
        for (int i = 0; i < _grid.domain().size_x + 1; i++) {
            vel[0] = (_field.u(i, j) + _field.u(i, j + 1)) * 0.5;
            vel[1] = (_field.v(i, j) + _field.v(i + 1, j)) * 0.5;
            Velocity->InsertNextTuple(vel);
        }
    }

    // Print Geometry from bottom to top
    for (int j = 1; j < _grid.domain().size_y + 1; j++) {
        for (int i = 1; i < _grid.domain().size_x + 1; i++) {
            double geometryData = _grid.cell(i, j).wall_id();
            Geometry->InsertNextTuple(&geometryData);
        }
    }

    // Add Pressure to Structured Grid
    structuredGrid->GetCellData()->AddArray(Pressure);

    // Add Velocity to Structured Grid
    structuredGrid->GetPointData()->AddArray(Velocity);

    // Add Geometry to Structured Grid
    structuredGrid->GetCellData()->AddArray(Geometry);

    // Add Temperature to Structured Grid
    if (_field.energy_on()) {
        structuredGrid->GetCellData()->AddArray(Temperature);
    }

#if IS_DETAIL_LOG
    // F Array
    vtkDoubleArray *farray = vtkDoubleArray::New();
    farray->SetName("F");
    farray->SetNumberOfComponents(1);

    // G Array
    vtkDoubleArray *garray = vtkDoubleArray::New();
    garray->SetName("G");
    garray->SetNumberOfComponents(1);

    // rs Array
    vtkDoubleArray *rsarray = vtkDoubleArray::New();
    rsarray->SetName("rs");
    rsarray->SetNumberOfComponents(1);

    for (int j = 1; j < _grid.domain().size_y + 1; j++) {
        for (int i = 1; i < _grid.domain().size_x + 1; i++) {
            double f = _field.f(i, j);
            farray->InsertNextTuple(&f);

            double g = _field.g(i, j);
            garray->InsertNextTuple(&g);
        }
    }

    for (int j = 1; j < _grid.domain().size_y + 1; j++) {
        for (int i = 1; i < _grid.domain().size_x + 1; i++) {
            double rs = _field.rs(i, j);
            rsarray->InsertNextTuple(&rs);
        }
    }
    structuredGrid->GetCellData()->AddArray(farray);
    structuredGrid->GetCellData()->AddArray(garray);
    structuredGrid->GetCellData()->AddArray(rsarray);
#endif

    // Write Grid
    vtkSmartPointer<vtkStructuredGridWriter> writer = vtkSmartPointer<vtkStructuredGridWriter>::New();

    // Create Filename
    std::string outputname =
        _dict_name + '/' + _case_name + "_" + std::to_string(my_rank) + "." + std::to_string(timestep) + ".vtk";

    writer->SetFileName(outputname.c_str());
    writer->SetInputData(structuredGrid);
    writer->Write();
}

void Case::build_domain(Domain &domain, int imax_domain, int jmax_domain) {

    if (_parallel_On)
    {
        /* eg: iproc = 4, jproc = 3
            (0,2)   (1,2)   (2,2)   (3,2)  ==  proc 2   proc 5   proc 8   proc 11
            (0,1)   (1,1)   (2,1)   (3,1)  ==  proc 1   proc 4   proc 7   proc 10
            (0,0)   (1,0)   (2,0)   (3,0)  ==  proc 0   proc 3   proc 6   proc 9
        */
        // if rank 0, compute decomposition boundaries like imin imax for each other processes
        // if rank 0, send the information out
        
        if (_rank == 0) //master
        {
            #if IS_DETAIL_LOG
                std::FILE* domainlog = std::fopen("./Domain.log", "w");
            #endif 
            /*
            * Compute partition size according to number of processes assigned
            */
            int subsize_x = 0;
            int subsize_y = 0;
            
            if(imax_domain % _iproc != 0){
                subsize_x = imax_domain / _iproc + 1; //ceiling
            }else{
                subsize_x = imax_domain / _iproc;
            }
            if(jmax_domain % _jproc != 0){
                subsize_y = jmax_domain / _jproc + 1; //ceiling
            }else{
                subsize_y = jmax_domain / _jproc;
            }
            /*
            * Compute domain min max according to partition size and communicate
            */
            #if IS_DETAIL_LOG
            char buffer[1024];
            snprintf(buffer, 1024, "domain subsize_x,y = %d,%d\n", subsize_x,subsize_y);
            std::cout << buffer;
            std::fprintf(domainlog, "%s",buffer);
            #endif 

            for (int i = 0; i < _iproc; ++i){
                 for (int j = 0; j < _jproc; ++j){
                    int i_domain_min = i * subsize_x;
                    int i_domain_max = std::min((i+1) * subsize_x,imax_domain);
                    int j_domain_min = j * subsize_y;
                    int j_domain_max = std::min((j+1) * subsize_y,jmax_domain);

                    int left_neighbour_rank = -1;
                    if (i > 0){
                        left_neighbour_rank = (i - 1) * _jproc + j;
                    }
                    int right_neighbour_rank = -1;
                    if (i < _iproc-1){
                        right_neighbour_rank = (i + 1) * _jproc + j;
                    }
                    int bottom_neighbour_rank= -1;
                    if (j > 0){
                        bottom_neighbour_rank = (i) * _jproc + j - 1;
                    }
                    int top_neighbour_rank = -1;
                     if (j < _jproc-1){
                        top_neighbour_rank = (i) * _jproc + j + 1;
                    }
                   
                    #if IS_DETAIL_LOG
                    char buffer[1024];
                    snprintf(buffer, 1024, "partition(%d,%d): min_i= %d, max_i= %d, min_j = %d, max_j = %d\n", 
                                i,j,i_domain_min,i_domain_max,j_domain_min,j_domain_max);
                    std::cout << buffer;
                    std::fprintf(domainlog, "%s",buffer);

                    snprintf(buffer, 1024, "partition(%d,%d): LeftNeighbourRank= %d, RightNeighbourRank= %d, BotNeighbourRank = %d, TopNeighbourRank = %d\n", 
                                i,j,left_neighbour_rank,right_neighbour_rank,bottom_neighbour_rank,top_neighbour_rank);
                    std::cout << buffer;
                    std::fprintf(domainlog, "%s",buffer);
                    #endif 

                    if(i == 0 && j == 0){ // master own domain
                        domain.imin = i_domain_min;
                        domain.jmin = j_domain_min;
                        domain.imax = i_domain_max + 2;
                        domain.jmax = j_domain_max + 2;
                        domain.size_x = i_domain_max - i_domain_min;
                        domain.size_y = j_domain_max - j_domain_min;

                        _left_neighbour_rank = left_neighbour_rank;
                        _right_neighbour_rank = right_neighbour_rank;
                        _bottom_neighbour_rank = bottom_neighbour_rank;
                        _top_neighbour_rank = top_neighbour_rank;
                    }
                    else // send slave domain
                    {
                        //communicate
                        int their_rank = i * _jproc + j;
                        std::cout << "Send Domain info to rank " << their_rank << "\n";
                        Communication::communicateDomainInfo(_rank,their_rank,i_domain_min,i_domain_max,j_domain_min,j_domain_max);
                        Communication::communicateNeighbourInfo(_rank,their_rank,left_neighbour_rank,right_neighbour_rank,bottom_neighbour_rank,top_neighbour_rank);
                    }
                }
            }
            std::cout << "rank " << _rank << " Successful Domain Communication.\n";
             #if IS_DETAIL_LOG
                std::fclose(domainlog);
            #endif 
        }
        else //slave
        {
            int i_domain_min = 0;
            int i_domain_max = 0;
            int j_domain_min = 0;
            int j_domain_max = 0;
            int their_rank = 0;
            Communication::communicateDomainInfo(_rank,their_rank,i_domain_min,i_domain_max,j_domain_min,j_domain_max);
            domain.imin = i_domain_min;
            domain.jmin = j_domain_min;
            domain.imax = i_domain_max + 2;
            domain.jmax = j_domain_max + 2;
            domain.size_x = i_domain_max - i_domain_min;
            domain.size_y = j_domain_max - j_domain_min;

            int left_neighbour_rank = -1;
            int right_neighbour_rank = -1;
            int bottom_neighbour_rank= -1;
            int top_neighbour_rank = -1;
            Communication::communicateNeighbourInfo(_rank,their_rank,left_neighbour_rank,right_neighbour_rank,bottom_neighbour_rank,top_neighbour_rank);
            _left_neighbour_rank = left_neighbour_rank;
            _right_neighbour_rank = right_neighbour_rank;
            _bottom_neighbour_rank = bottom_neighbour_rank;
            _top_neighbour_rank = top_neighbour_rank;
             std::cout << "rank " << _rank << " Successful Domain Communication.\n";
        }
    }
    else // not parallel
    {
        domain.imin = 0;
        domain.jmin = 0;
        domain.imax = imax_domain + 2;
        domain.jmax = jmax_domain + 2;
        domain.size_x = imax_domain;
        domain.size_y = jmax_domain;
    }
}
