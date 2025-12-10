#include "left_ventricle.hpp"
#include <iostream>


int main(int argc, char* argv[]){
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);

  LV model = LV("../mesh/prova_mesh.msh", 2);
  model.setup();
  model.assemble_system();


  model.solve_newton();
  
}
