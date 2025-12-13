SetFactory("OpenCASCADE");

// Crea il cubo di lato 1
Box(1) = {0, 0, 0, 1, 1, 1};

// Sincronizza la geometria
Synchronize();

// Crea Physical Surface per tutte le 6 facce
Physical Surface("Dirichlet") = {1,2,3};
Physical Surface("Neumann") = {4,5,6};

// Crea Physical Volume per il dominio
Physical Volume("Domain") = {1};
