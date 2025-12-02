#include "left_ventricle.hpp"

void
LV::setup()
{
  pcout << "===============================================" << std::endl;

  // Create the mesh.
  {
    pcout << "Initializing the mesh" << std::endl;

    Triangulation<dim> mesh_serial;

    GridIn<dim> grid_in;
    grid_in.attach_triangulation(mesh_serial);

    std::ifstream grid_in_file(mesh_file_name);
    grid_in.read_msh(grid_in_file);

    GridTools::partition_triangulation(mpi_size, mesh_serial);
    const auto construction_data = TriangulationDescription::Utilities::
      create_description_from_triangulation(mesh_serial, MPI_COMM_WORLD);
    mesh.create_triangulation(construction_data);

    pcout << "  Number of elements = " << mesh.n_global_active_cells()
          << std::endl;
  }
  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the finite element space. This is the same as in serial codes.
  {
    pcout << "Initializing the finite element space" << std::endl;

    fe = std::make_unique<FE_SimplexP<dim>>(r);

    pcout << "  Degree                     = " << fe->degree << std::endl;
    pcout << "  DoFs per cell              = " << fe->dofs_per_cell
          << std::endl;

    quadrature = std::make_unique<QGaussSimplex<dim>>(fe->degree + 1);
    quadrature_face = std::make_unique<QGaussSimplex<dim - 1>>(fe->degree + 1);

    pcout << "  Quadrature points per cell = " << quadrature->size()
          << std::endl;
    pcout << "  Quadrature points per face = " << quadrature_face->size()
          << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the DoF handler.
  {
    pcout << "Initializing the DoF handler" << std::endl;

    dof_handler.reinit(mesh);
    dof_handler.distribute_dofs(*fe);

    // We retrieve the set of locally owned DoFs, which will be useful when
    // initializing linear algebra classes.
    locally_owned_dofs = dof_handler.locally_owned_dofs();
    DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

    pcout << "  Number of DoFs = " << dof_handler.n_dofs() << std::endl;
  }

  pcout << "-----------------------------------------------" << std::endl;

  // Initialize the linear system.
  {
    pcout << "Initializing the linear system" << std::endl;

    pcout << "  Initializing the sparsity pattern" << std::endl;

    // To initialize the sparsity pattern, we use Trilinos' class, that manages
    // some of the inter-process communication.
    TrilinosWrappers::SparsityPattern sparsity(locally_owned_dofs,
                                               MPI_COMM_WORLD);
    DoFTools::make_sparsity_pattern(dof_handler, sparsity);
    sparsity.compress();

    pcout << "  Initializing the Jacobian matrix" << std::endl;
    jacobian_matrix.reinit(sparsity);


    pcout << "  Initializing the system right-hand side" << std::endl;
    system_rhs.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    pcout << "  Initializing the solution vector" << std::endl;
    solution_owned.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    delta_owned.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    
    solution.reinit(locally_owned_dofs,locally_relevant_dofs, MPI_COMM_WORLD);
  }
}

void
LV::assemble()
{
  pcout << "===============================================" << std::endl;

  pcout << "  Assembling the linear system" << std::endl;

  const unsigned int dofs_per_cell = fe->dofs_per_cell;
  const unsigned int n_q           = quadrature->size();
  const unsigned int n_q_face      = quadrature_face->size();

  FEValues<dim> fe_values(*fe,
                          *quadrature,
                          update_values | update_gradients |
                            update_quadrature_points | update_JxW_values);
  
  FEFaceValues<dim> fe_face_values(*fe,
                                   *quadrature_face,
                                   update_values | update_normal_vectors |
                                     update_JxW_values);

  FullMatrix<double> cell_matrix(dofs_per_cell, dofs_per_cell);
  Vector<double>     cell_rhs(dofs_per_cell);

  std::vector<types::global_dof_index> dof_indices(dofs_per_cell);

  jacobian_matrix = 0.0;
  system_rhs      = 0.0;

//necesary for Robin
std::vector<double> solution_loc(n_q_face);
std::vector<Tensor<1, dim>> solution_gradient_loc(n_q_face);

//forse una cosa del genere per la pressione e/o F????



  for (const auto &cell : dof_handler.active_cell_iterators())
    {
      // If current cell is not owned locally, we skip it.
      if (!cell->is_locally_owned())
        continue;

      // On all other cells (which are owned by current process), we perform the
      // assembly as usual.

      fe_values.reinit(cell);

      cell_matrix = 0.0;
      cell_rhs    = 0.0;
      

      
      cell->get_dof_indices(dof_indices);

      for (unsigned int q = 0; q < n_q; ++q)
        {
          input_data data;
          //riempire data con i valori corretti
          
          Tensor<2,dim> F=;//TODO calcolare il gradiente della mappatura in quel punto

          Tensor<2,dim> P=compute_P(data);
          Tensor<4,dim> dP_dF=compute_dP_dF(data);
          

          for (unsigned int i = 0; i < dofs_per_cell; ++i){
            cell_rhs(i)+=scalar_product(P,fe_values.shape_grad(i,q))*fe_values.JxW(q);
          }
        }

          if(cell->at_boundary()){
            for (unsigned int f=0;f< cell->n_faces();++f){
              if (cell->face(f)->at_boundary() && cell->face(f)->boundary_id()== 3) {
                fe_face_values.reinit(cell,f);
                  
                  for (unsigned int q=0; q<n_q_face;++q){

                    //da capire la dipendenza rispetto a q


                    Tensor <2,dim> Fh=F(q);
                    
                    //Tensor<2,dim> H=det(Fh)*transpose(inverse(Fh)); 

                    double pressure=pressure(q) ; //TODO definire la pressione in quel punto
                    
                    for (unsigned int i=0;i<dofs_per_cell;++i){
                        cell_rhs[i] += pressure*scalar_product(H*fe_face_values.normal_vector(q),fe_face_values.shape_value(i,q))*fe_face_values.JxW(q);
                  
                      //compute derivative and calculate the cell_matrix[i,j]=d cell_rhs[i]/d u_j * (-1)
                      }
                  }
              }


      //ROBIN TERM
      //mi serve sapere uh sulla faccia e onn l'ho mai trovato nei laboratori quindi boh
            
      /*
      if (cell->face(f)->at_boundary() && cell->face(f)->boundary_id()== 1) {
        fe_face_values.reinit(cell,f);
        

        
        std::vector<double> solution_loc(n_q_face);
        //std::vector<Tensor<1, dim>> solution_gradient_loc(n_q_face);
        fe_face_values.get_function_values(solution, solution_loc);
        //fe_face_values.get_function_gradients(solution, solution_gradient_loc);



        for (unsigned int q=0; q<n_q_face;++q){
          
        
        for (unsigned int i=0;i<dofs_per_cell;++i){
          cell_rhs(i) += alpha * solution_loc[q] *
                               fe_face_values.shape_value(i, q) *
                               fe_face_values.JxW(q);
        }
        
        for (unsigned int i=0;i<dofs_per_cell;++i){
          for (unsigned int j=0;j<dofs_per_cell;++j){
            cell_matrix(i,j)+= alpha*fe_face_values.shape_value(i,q)*
                              fe_face_values.shape_value(j,q)*
                              fe_face_values.JxW(q);
        
      }
    }   
    */



  }
}





      cell->get_dof_indices(dof_indices);

      system_matrix.add(dof_indices, cell_matrix);
      system_rhs.add(dof_indices, cell_rhs);


}
  system_matrix.compress(VectorOperation::add);
  system_rhs.compress(VectorOperation::add);

  // Dirichlet Boundary conditions.
  {
    std::map<types::global_dof_index, double> boundary_values;

    std::map<types::boundary_id, const Function<dim> *> boundary_functions;


    boundary_functions[2] = &function_g;


    VectorTools::interpolate_boundary_values(dof_handler,
                                             boundary_functions,
                                             boundary_values);

    MatrixTools::apply_boundary_values(
      boundary_values, system_matrix, solution, system_rhs, true);
  }
}
void
LV::solve_linear_system(){

   SolverControl solver_control(1000, 1e-6 * system_rhs.l2_norm());

  SolverGMRES<TrilinosWrappers::MPI::Vector> solver(solver_control);
  TrilinosWrappers::PreconditionSSOR      preconditioner;
  preconditioner.initialize(
    jacobian_matrix, TrilinosWrappers::PreconditionSSOR::AdditionalData(1.0));

  solver.solve(jacobian_matrix, delta_owned, system_rhs, preconditioner);
  pcout << "  " << solver_control.last_step() << " GMRES iterations" << std::endl;
}

}

void
LV::solve_newton()
{
  const unsigned int n_max_iters        = 1000;
  const double       residual_tolerance = 1e-6;

  unsigned int n_iter        = 0;
  double       residual_norm = residual_tolerance + 1;

  // We apply the boundary conditions to the initial guess (which is stored in
  // solution_owned and solution).
  {
    IndexSet dirichlet_dofs = DoFTools::extract_boundary_dofs(dof_handler);
    dirichlet_dofs          = dirichlet_dofs & dof_handler.locally_owned_dofs();

    function_g.set_time(time);

    TrilinosWrappers::MPI::Vector vector_dirichlet(solution_owned);
    VectorTools::interpolate(dof_handler, function_g, vector_dirichlet);

    for (const auto &idx : dirichlet_dofs)
      solution_owned[idx] = vector_dirichlet[idx];

    solution_owned.compress(VectorOperation::insert);
    solution = solution_owned;
  }

  while (n_iter < n_max_iters && residual_norm > residual_tolerance)
    {
      assemble_system();
      residual_norm = system_rhs.l2_norm();

      pcout << "  Newton iteration " << n_iter << "/" << n_max_iters
            << " - ||r|| = " << std::scientific << std::setprecision(6)
            << residual_norm << std::flush;

      // We actually solve the system only if the residual is larger than the
      // tolerance.
      if (residual_norm > residual_tolerance)
        {
          solve_linear_system();

          solution_owned += delta_owned;
          solution = solution_owned;
        }
      else
        {
          pcout << " < tolerance" << std::endl;
        }

      ++n_iter;
    }
}





void
LV::output() const
{
  pcout << "===============================================" << std::endl;

  IndexSet locally_relevant_dofs;
  DoFTools::extract_locally_relevant_dofs(dof_handler, locally_relevant_dofs);

  // To correctly export the solution, each process needs to know the solution
  // DoFs it owns, and the ones corresponding to elements adjacent to the ones
  // it owns (the locally relevant DoFs, or ghosts). We create a vector to store
  // them.
  TrilinosWrappers::MPI::Vector solution_ghost(locally_owned_dofs,
                                               locally_relevant_dofs,
                                               MPI_COMM_WORLD);

  // This performs the necessary communication so that the locally relevant DoFs
  // are received from other processes and stored inside solution_ghost.
  solution_ghost = solution;

  // Then, we build and fill the DataOut class as usual.
  DataOut<dim> data_out;
  data_out.add_data_vector(dof_handler, solution_ghost, "solution");

  // We also add a vector to represent the parallel partitioning of the mesh.
  std::vector<unsigned int> partition_int(mesh.n_active_cells());
  GridTools::get_subdomain_association(mesh, partition_int);
  const Vector<double> partitioning(partition_int.begin(), partition_int.end());
  data_out.add_data_vector(partitioning, "partitioning");

  data_out.build_patches();

  const std::filesystem::path mesh_path(mesh_file_name);
  const std::string output_file_name = "output-" + mesh_path.stem().string();

  // Finally, we need to write in a format that supports parallel output. This
  // can be achieved in multiple ways (e.g. XDMF/H5). We choose VTU/PVTU files,
  // because the interface is nice and it is quite robust.
  data_out.write_vtu_with_pvtu_record("./",
                                      output_file_name,
                                      0,
                                      MPI_COMM_WORLD);

  pcout << "Output written to " << output_file_name << std::endl;

  pcout << "===============================================" << std::endl;
}


  Tensor<2,dim> compute_P(input_data data) const
  {
    Tensor<2,dim> C=transpose(data.F) *data.F;

    Tensor<2,dim> E=0.5*(C- unit_symmetric_tensor<dim>());

    double E_ff= (E*data.f0)*data.f0;
    double E_ss= (E*data.s0)*data.s0;
    double E_nn= (E*data.n0)*data.n0;

//aggiungere altre direzioni se si aggiungono sopra
  double Q=data.b_f*E_ff*E_ff+
            data.b_s*E_ss*E_ss+
            data.b_n*E_nn*E_nn;


    double W= 0.5*(std::exp(Q)-1);

    double I4_f= data.f0* C* data.f0;

    Tensor<2,dim> P_att= Ta*((data.F*data.f0)*data.f0)/
                        std::sqrt(I4_f);


    Tensor<2,dim> P_pass= ;//TODO capire come si fa la derivata di W rispetto a F

    //si potrebbe fare con la regola della catena ma non sarebbe scalabile in alcun modo

  }

Tensor<4, LV::dim> LV::compute_dP_dF(const Vector<double> &d,const Point<dim> &p) const
{
    Tensor<4, dim> tangent;
    //TODO 
    //tensore tangente dP/dF (chissà)
    return tangent;
}



