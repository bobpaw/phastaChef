#include <cassert>
#include <random>
#include <iostream>
#include <chrono>

#include <apfMesh2.h>
#include <apfShape.h>
#include <PCU.h>
#include <lionPrint.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <phInput.h>

namespace pc {
  using Shocks = std::vector<std::vector<apf::MeshEntity*>>;

  // Function we are testing:
  Shocks labelShocksSerial(const ph::Input& in, apf::Mesh2* m);
}

void print_usage(std::ostream& str, char* argv0) {
  str << "USAGE: " << argv0 << " <modelFileName> <meshFileName>" << std::endl;
}

int main (int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();

  lion_set_verbosity(1);

  if (argc != 3) {
    print_usage(std::cerr, argv[0]);
    return -1;
  }

  apf::Mesh2* m = nullptr;

  ph::Input in;
  in.modelFileName = std::string(argv[1]);
  in.meshFileName = std::string(argv[2]);

  // Init gmi.
  gmi_register_mesh();

  // Load mesh.
  m = apf::loadMdsMesh(in.modelFileName.c_str(), in.meshFileName.c_str());

  // Enable the 1ms literal.
  using namespace std::literals;

  auto clock_start = std::chrono::system_clock::now();
  pc::Shocks shocks = pc::labelShocksSerial(in, m);
  auto clock_end = std::chrono::system_clock::now();

	std::cout << "Labeled " << shocks.size() << " distinct shocks." << std::endl;

  std::cout << "Labeling: " << (clock_end - clock_start) / 1ms << "ms." << std::endl;

  // Write post-label VTK files.
  apf::writeVtkFiles("td_label.vtk", m);

  // Write MDS mesh.
  m->writeNative("labeled.smb");

  // Cleanup.
  m->destroyNative();
  apf::destroyMesh(m);

  PCU_Comm_Free();
  MPI_Finalize();

  return 0;
}
