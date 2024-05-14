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
#include <phInput.h>
#include <maSize.h>

#include "pcShockParam.h"
#include "ShockFunc.h"

namespace pc {
  // Function we are testing.
  void defragmentShocksSerial(const ph::Input& in, apf::Mesh2* m);
}

void print_usage(std::ostream& str, char* argv0) {
  str << "USAGE: " << argv0 << " <modelFileName> <meshFileName>"
      << " type:arg,arg,... <tol_alpha> <fragmentation>"
      << std::endl;
}

int mockShockParam(apf::Mesh* m, pc::_test::ShockFunc* sf) {
  assert(m);
  assert(sf);

  apf::Field* Shock_Param =
      apf::createField(m, "Shock_Param", apf::SCALAR, apf::getConstant(3));

  int n = 0;

  apf::MeshIterator* it = m->begin(3);
  for (apf::MeshEntity* e = m->iterate(it); e; e = m->iterate(it)) {
    if ((*sf)(m, e)) {
      apf::setScalar(Shock_Param, e, 0, pc::ShockParam::SHOCK);
      ++n;
    } else {
      apf::setScalar(Shock_Param, e, 0, pc::ShockParam::NONE);
    }
  }
  m->end(it);

  return n;
}

// Compute difference between two time values.
template <typename T>
double my_tdiff(T start, T end) {
  return (end - start) / std::chrono::duration<double>(1.0);
}

/**
 * Remove a random sample of shock elements.
 *
 * @param percent The percentage of shock elements to remove.
 */
int fragmentShockParam(apf::Mesh* m, int percent) {
  assert(m);
  assert(0 <= percent && percent <= 100);

  apf::Field* Shock_Param = m->findField("Shock_Param");

  std::vector<apf::MeshEntity*> shock_elements;

  apf::MeshIterator* it = m->begin(3);
  for (apf::MeshEntity* e = m->iterate(it); e; e = m->iterate(it)) {
    if (apf::getScalar(Shock_Param, e, 0) == pc::ShockParam::SHOCK) {
      shock_elements.push_back(e);
    }
  }
  m->end(it);

  std::mt19937 mt(std::random_device{}());
  std::uniform_int_distribution<int> index;
  int n = shock_elements.size() * percent / 100;

  for (int i = 0; i < n; ++i) {
    int id = index(mt) % shock_elements.size();
    apf::setScalar(Shock_Param, shock_elements[id], 0, pc::ShockParam::NONE);
    shock_elements.erase(shock_elements.begin() + id);
  }

  return n;
}

double calculateTolerance(double alpha, apf::Mesh2* m) {
  return alpha * ma::getAverageEdgeLength(m);
}

int main (int argc, char* argv[]) {
  // MPI/PCU init?
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();

  // lion verbosity
  lion_set_verbosity(1);

  if (argc != 6) {
    print_usage(std::cerr, argv[0]);
    return -1;
  }

  apf::Mesh2* m = nullptr;
  pc::_test::ShockFunc* sf = nullptr;

  ph::Input in;
  in.modelFileName = std::string(argv[1]);
  in.meshFileName = std::string(argv[2]);

  // Init gmi.
  gmi_register_mesh();

  // Load mesh.
  m = apf::loadMdsMesh(in.modelFileName.c_str(), in.meshFileName.c_str());

  int fragmentation;

  try {
    double tol_alpha = std::atof(argv[4]);
    double tol = calculateTolerance(tol_alpha, m);
    fragmentation = std::atoi(argv[5]);

    std::cout << "Tolerance: " << tol << std::endl;

    sf = pc::_test::ShockFunc::makeFromString(argv[3], tol);
		if (sf == nullptr) throw std::invalid_argument("ShockFunc argument");
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    print_usage(std::cerr, argv[0]);
    return -1;
  }

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

  // Print mesh statistics.
  std::cout << "Mesh bounding box size: (" << m_size.x() << ", " << m_size.y() << ", " << m_size.z() << ")." << std::endl;
  std::cout << "Mesh elements (D=3): " << m->count(3) << std::endl;

  auto clock_start = std::chrono::steady_clock::now();
  int n = mockShockParam(m, sf);
  auto clock_end = std::chrono::steady_clock::now();
  std::cout << "Assigned " << n << " shock elements: " << my_tdiff(clock_start, clock_end) << " seconds." << std::endl;

  // Write full shock VTK file.
  apf::writeVtkFiles("td_shock.vtk", m);

  // Fragmentation step
  clock_start = std::chrono::steady_clock::now();
  n = fragmentShockParam(m, fragmentation);
  clock_end = std::chrono::steady_clock::now();
  std::cout << "Fragmented (removed) " << n << " shock elements: " << my_tdiff(clock_start, clock_end) << " seconds." << std::endl;

  // Write pre-defrag VTK file.
  apf::writeVtkFiles("td_frag.vtk", m);

  // Defragmentation step.
  clock_start = std::chrono::steady_clock::now();
  pc::defragmentShocksSerial(in, m);
  clock_end = std::chrono::steady_clock::now();
  std::cout << "Defragmentation: " << my_tdiff(clock_start, clock_end) << " seconds." << std::endl;
  
  // Write post-defrag VTK files.
  apf::writeVtkFiles("td_defrag.vtk", m);

  // Write mesh file.
  m->writeNative("defragmented.smb");

  // Cleanup
  m->destroyNative();
  apf::destroyMesh(m);

  PCU_Comm_Free();
  MPI_Finalize();

  return 0;
}
