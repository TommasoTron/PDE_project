SetFactory("OpenCASCADE"); //Necessary for doing the intern subraction

//Given parameters from: "A transmurally heterogeneous orthotropic activation model for ventricular contraction and its numerical validation"

//Ventricle geometry with fiber orientation

//Mesh size
ms = 0.5;

//Focal length 
d = 2.91;
xi_in = 0.6;       // inter ellipsoidale
xi_out = 1.02;     // outer ellipsoidale
z_cut = 1.19;      // base

//Calculate the 4 semi axes of the ellipsoids (2 inner and 2 outer)
R_in_xy = d * Sinh(xi_in);
R_in_z  = d * Cosh(xi_in);

R_out_xy = d * Sinh(xi_out);
R_out_z  = d * Cosh(xi_out);

//Create the outer ellipsoid
v_out = newv; // Create a new volume
Sphere(v_out) = {0, 0, 0, 1}; // Unitary sphere in the origin with tag = newv
Dilate {{0, 0, 0}, {R_out_xy, R_out_xy, R_out_z}} { Volume{v_out}; } 

//Create the inner ellipsoid 
v_in = newv; 
Sphere(v_in) = {0, 0, 0, 1};
Dilate {{0, 0, 0}, {R_in_xy, R_in_xy, R_in_z}} { Volume{v_in}; }

//Create the ventricle by subtracting the inner from the outer ellipsoid
v_wall = newv;
BooleanDifference(v_wall) = { Volume{v_out}; Delete; }{ Volume{v_in}; Delete; };

// Cut the ventricle at the base (z = z_cut)
Box(1000) = {-50, -50, z_cut, 100, 100, 100};

// With the same concept we cancel the upper part of the ventricle
BooleanDifference{ Volume{v_wall}; Delete; }{ Volume{1000}; Delete; }


// --- MODIFIED SECTION: Precise Selection ---

eps = 1e-2; // Slightly larger tolerance for safety

// 1. Base: Selection by Z-height (Standard)
s_base() = Surface In BoundingBox{-100, -100, z_cut-eps, 100, 100, z_cut+eps};

// 2. Endo: Selection by probe point at equator (X = R_in_xy)
// We look for the surface passing exactly through (R_in_xy, 0, 0)
s_endo() = Surface In BoundingBox{R_in_xy-eps, -eps, -eps, R_in_xy+eps, eps, eps};

// 3. Epi: Selection by probe point at equator (X = R_out_xy)
// We look for the surface passing exactly through (R_out_xy, 0, 0)
s_epi() = Surface In BoundingBox{R_out_xy-eps, -eps, -eps, R_out_xy+eps, eps, eps};

// -------------------------------------------

Physical Surface("BASE", 100) = { 2 }; 
Physical Surface("ENDO_IN", 200) = { 3 };
Physical Surface("EPI_OUT", 300)  = { 1 };
Physical Volume("VOLUME", 400) = { v_wall };

//Mesh size specification
Mesh.CharacteristicLengthMin = ms;
Mesh.CharacteristicLengthMax = ms;

// Base (ID 2) -> ROSSO
Color Red { Surface{ 2 }; }

// Endocardio (ID 3) -> VERDE
Color Green { Surface{ 3 }; }

// Epicardio (ID 1) -> VIOLA
Color Purple { Surface{ 1 }; }

// Volume -> Giallo
Color Yellow { Volume{ v_wall }; }