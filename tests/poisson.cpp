#include "mfem.hpp"
#include <fstream>
#include <iostream>

using namespace std;
using namespace mfem;

// Device sample runs:
//               ex1 -pa -d cuda
//               ex1 -fa -d cuda
//               ex1 -pa -d raja-cuda
//             * ex1 -pa -d raja-hip
//               ex1 -pa -d occa-cuda
//               ex1 -pa -d raja-omp
//               ex1 -pa -d occa-omp
//               ex1 -pa -d ceed-cpu
// 
//               ex1 -pa -d ceed-cpu-something-blocked

// note: partial assembly only supported on tensor product elements (?)

int main(int argc, char *argv[])
{


   // 1. Parse command-line options.
   //const char *mesh_file = "../data/star.mesh";
   const char *mesh_file = "/Users/sam/Dropbox/meshes/rocker_arm_hex.vtk"; // TODO: DELETE
   int order = 1;
   bool pa = false;
   bool fa = false;
   const char *device_config = "cpu";

   int num_refinements = 0;

   OptionsParser args(argc, argv);
   args.AddOption(&mesh_file, "-m", "--mesh",
                  "Mesh file to use.");
   args.AddOption(&order, "-o", "--order",
                  "Finite element order (polynomial degree) or -1 for"
                  " isoparametric space.");
   args.AddOption(&pa, "-pa", "--partial-assembly", "-no-pa",
                  "--no-partial-assembly", "Enable Partial Assembly.");
   args.AddOption(&fa, "-fa", "--full-assembly", "-no-fa",
                  "--no-full-assembly", "Enable Full Assembly.");
   args.AddOption(&num_refinements, "-ref", "--refinement", "how many refinements to apply to the mesh");
   args.AddOption(&device_config, "-d", "--device",
                  "Device configuration string, see Device::Configure().");
#ifdef MFEM_USE_CEED
   bool algebraic_ceed = false;
   args.AddOption(&algebraic_ceed, "-a", "--algebraic", "-no-a", "--no-algebraic",
                  "Use algebraic Ceed solver");
#endif
   args.Parse();
   if (!args.Good())
   {
      args.PrintUsage(cout);
      return 1;
   }
   args.PrintOptions(cout);

   // 2. Enable hardware devices such as GPUs, and programming models such as
   //    CUDA, OCCA, RAJA and OpenMP based on command line options.
   Device device(device_config);
   device.Print();

   // 3. Read the mesh from the given mesh file. We can handle triangular,
   //    quadrilateral, tetrahedral, hexahedral, surface and volume meshes with
   //    the same code.
   Mesh mesh(mesh_file, 1, 1);
   int dim = mesh.Dimension();

   // 4. Refine the mesh to increase the resolution. In this example we do
   //    'ref_levels' of uniform refinement. We choose 'ref_levels' to be the
   //    largest number that gives a final mesh with no more than 50,000
   //    elements.
   for (int l = 0; l < num_refinements; l++)
   {
      mesh.UniformRefinement();
   }

   // 5. Define a finite element space on the mesh. Here we use continuous
   //    Lagrange finite elements of the specified order. If order < 1, we
   //    instead use an isoparametric/isogeometric space.
   FiniteElementCollection *fec = new H1_FECollection(order, dim);
   FiniteElementSpace fespace(&mesh, fec);
   cout << "Number of finite element unknowns: " << fespace.GetTrueVSize() << endl;

   LinearForm b(&fespace);
   ConstantCoefficient one(1.0);
   b.AddDomainIntegrator(new DomainLFIntegrator(one));
   if (pa) { 
     b.UseFastAssembly(true);
   }
   b.Assemble();

   // 8. Define the solution vector x as a finite element grid function
   //    corresponding to fespace. Initialize x with initial guess of zero,
   //    which satisfies the boundary conditions.
   GridFunction x(&fespace);
   x.Randomize();

   // 9. Set up the bilinear form a(.,.) on the finite element space
   //    corresponding to the Laplacian operator -Delta, by adding the Diffusion
   //    domain integrator.
   BilinearForm a(&fespace);
   a.AddDomainIntegrator(new MassIntegrator(one));
   a.AddDomainIntegrator(new DiffusionIntegrator(one));

   if (pa) { a.SetAssemblyLevel(AssemblyLevel::PARTIAL); }
   if (fa) { a.SetAssemblyLevel(AssemblyLevel::FULL); }

   a.Assemble();
   a.Finalize();

   Vector R(x.Size());

   for (int i = 0; i < 5; i++) {
     tic();
     a.Mult(x, R);
     R += b;
     double elapsed = toc();
     cout << elapsed * 1000 << "ms" << endl;
   }

   delete fec;

   return 0;
}