#include <cassert>
#include <random>
#include <iostream>
#include <chrono>

#include <apfMesh2.h>
#include <apfShape.h>
#include <PCU.h>
#include <lionPrint.h>
#include <phstream.h>
#include <gmi_mesh.h>
#include <apfMDS.h>

#include <chef.h>
#include <phInput.h>

#include "pcShockParam.h"
#include "ShockFunc.h"

namespace pc {
  // Function we are testing.
  void defragmentShocksSerial(const ph::Input& in, apf::Mesh2* m);
}

void print_usage(std::ostream& str, char* argv0) {
	str << "USAGE: " << argv0
			<< " <modelFileName> <meshFileName> <A> <B> <C> <D> <tol>"
			<< std::endl;
}

int mockShockParam(apf::Mesh* m, pc::_test::ShockFunc* sf) {
	assert(m);
	assert(sf);

  std::mt19937 mt(std::random_device{}());
  std::uniform_int_distribution<int> dist(0, 100);

	apf::Field* Shock_Param =
			apf::createField(m, "Shock_Param", apf::SCALAR, apf::getConstant(3));

  int n = 0;

	apf::MeshIterator* it = m->begin(3);
  for (apf::MeshEntity* e = m->iterate(it); e; e = m->iterate(it)) {
    if ((*sf)(m, e) && dist(mt) < 30) {
      apf::setScalar(Shock_Param, e, 0, pc::ShockParam::SHOCK);
      ++n;
    } else {
      apf::setScalar(Shock_Param, e, 0, pc::ShockParam::NONE);
    }
  }
  m->end(it);

  return n;
}

int main (int argc, char* argv[]) {
  // MPI/PCU init?
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();

  // lion verbosity
  lion_set_verbosity(1);

  if (argc != 8) {
    print_usage(std::cerr, argv[0]);
    return -1;
  }

  apf::Mesh2* m = nullptr;
  pc::_test::ShockFunc* sf = nullptr;

  ph::Input in;

  try {
    in.modelFileName = std::string(argv[1]);
    in.meshFileName = std::string(argv[2]);
    double a = std::atof(argv[3]), b = std::atof(argv[4]);
    double c = std::atof(argv[5]), d = std::atof(argv[6]);
    double tol = std::atof(argv[7]);

    // FIXME: Replace tol with tolerance based on a fraction of mesh size.
    sf = new pc::_test::PlanarShockFunc(a, b, c, d, tol);
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    print_usage(std::cerr, argv[0]);
    return -1;
  }

  // Init gmi.
  gmi_register_mesh();

  // Load mesh.
  m = apf::loadMdsMesh(in.modelFileName.c_str(), in.meshFileName.c_str());

  apf::Vector3 c_min, c_max;
  apf::MeshIterator* it = m->begin(0);
  for (apf::MeshEntity* e = m->iterate(it); e; e = m->iterate(it)) {
    apf::Vector3 lc = apf::getLinearCentroid(m, e);
    c_min.x() = c_min.x() < lc.x() ? c_min.x() : lc.x();
    c_min.y() = c_min.y() < lc.y() ? c_min.y() : lc.y();
    c_min.z() = c_min.z() < lc.z() ? c_min.z() : lc.z();
    c_max.x() = c_max.x() > lc.x() ? c_max.x() : lc.x();
    c_max.y() = c_max.y() > lc.y() ? c_max.y() : lc.y();
    c_max.z() = c_max.z() > lc.z() ? c_max.z() : lc.z();
  }
  m->end(it);

  apf::Vector3 m_size = c_max - c_min;

  // Enable the 1ms literal.
  using namespace std::literals;

  // Print mesh statistics.
  std::cout << "Mesh size: (" << m_size.x() << ", " << m_size.y() << ", " << m_size.z() << ")." << std::endl;
  std::cout << "Mesh elements (D=3): " << m->count(3) << std::endl;

  auto clock_start = std::chrono::system_clock::now();
  int n = mockShockParam(m, sf);
  auto clock_end = std::chrono::system_clock::now();

  // Print Shock_Param statistics.
  // - number of elements.
  std::cout << "Assigned " << n << " shock elements: " << (clock_end - clock_start) / 1ms << "ms." << std::endl;

  // Write pre-defrag VTK file.
  apf::writeVtkFiles("defrag_pre.vtk", m);

  // Defrag step
  clock_start = std::chrono::system_clock::now();
  pc::defragmentShocksSerial(in, m);
  clock_end = std::chrono::system_clock::now();
  std::cout << "Defragmentation: " << (clock_end - clock_start) / 1ms << "ms." << std::endl;
  
  // Write post-defrag VTK files.
  apf::writeVtkFiles("defrag_post.vtk", m);

  // Write Simmetrix file.

  // Cleanup
  m->destroyNative();
  apf::destroyMesh(m);

  PCU_Comm_Free();
  MPI_Finalize();

  return 0;
}
