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
#include <ma.h>
#include <maSize.h>

#include "pcShockParam.h"
#include "ShockFunc.h"

namespace pc {
  // Function we are testing.
  void defragmentShocksSerial(const ph::Input& in, apf::Mesh2* m);
  using Shocks = std::vector<std::vector<apf::MeshEntity*>>;
	Shocks labelShocksSerial(const ph::Input& in, apf::Mesh2* m);
  void denoiseShocksSerial(const ph::Input& in, apf::Mesh2* m, Shocks& shocks, int minsize = 10);
}

void print_usage(std::ostream& str, char* argv0) {
  str << "USAGE: " << argv0 << " <modelFileName> <meshFileName>"
      << " type:arg,arg,... type2:arg,arg,... <tol_alpha> <fragmentation>"
      << std::endl;
}

int mockShockParam(apf::Mesh* m, pc::_test::ShockFunc* sf) {
  assert(m);
  assert(sf);

  apf::Field* Shock_Param = m->findField("Shock_Param");
  if (Shock_Param == nullptr) {
    Shock_Param = apf::createField(m, "Shock_Param", apf::SCALAR, apf::getConstant(3));
    apf::zeroField(Shock_Param);
  }

  int n = 0;

  apf::MeshIterator* it = m->begin(3);
  for (apf::MeshEntity* e = m->iterate(it); e; e = m->iterate(it)) {
    if ((*sf)(m, e)) {
      apf::setScalar(Shock_Param, e, 0, pc::ShockParam::SHOCK);
      ++n;
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

void addNoise(apf::Mesh* m) {
  std::mt19937 mt(std::random_device{}());
  std::uniform_int_distribution<int> chance(0, 1000);

  apf::Field* Shock_Param = m->findField("Shock_Param");

  int max_noise = 10;

  apf::MeshIterator* it = m->begin(3);
  for (apf::MeshEntity* e = m->iterate(it); e && max_noise > 0; e = m->iterate(it)) {
    if (chance(mt) < 5) {
      apf::setScalar(Shock_Param, e, 0, pc::ShockParam::SHOCK);
      --max_noise;
    }
  }
  m->end(it);
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

  if (argc != 7) {
    print_usage(std::cerr, argv[0]);
    return -1;
  }

  apf::Mesh2* m = nullptr;
  pc::_test::ShockFunc* sf1 = nullptr, *sf2 = nullptr;

  ph::Input in;
  in.modelFileName = std::string(argv[1]);
  in.meshFileName = std::string(argv[2]);

  // Init gmi.
  gmi_register_mesh();

  // Load mesh.
  m = apf::loadMdsMesh(in.modelFileName.c_str(), in.meshFileName.c_str());
  ma::runUniformRefinement(m);
  m->verify();

  int fragmentation;

  try {
    double tol_alpha = std::atof(argv[5]);
    double tol = calculateTolerance(tol_alpha, m);
    fragmentation = std::atoi(argv[6]);

    std::cout << "Tolerance: " << tol << std::endl;

    sf1 = pc::_test::ShockFunc::makeFromString(argv[3], 4);
    sf2 = pc::_test::ShockFunc::makeFromString(argv[4], 25);
		if (sf1 == nullptr) throw std::invalid_argument("ShockFunc argument");
		if (sf2 == nullptr) throw std::invalid_argument("ShockFunc argument");
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
  int n = mockShockParam(m, sf1);
  n += mockShockParam(m, sf2);
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

  addNoise(m);

  apf::writeVtkFiles("td_noise.vtk", m);

  // Defragmentation step.
  clock_start = std::chrono::steady_clock::now();
  pc::defragmentShocksSerial(in, m);
  clock_end = std::chrono::steady_clock::now();
  std::cout << "Defragmentation: " << my_tdiff(clock_start, clock_end) << " seconds." << std::endl;

  // Write post-defrag VTK files.
  apf::writeVtkFiles("td_defrag.vtk", m);

  pc::Shocks shocks = pc::labelShocksSerial(in, m);

  apf::writeVtkFiles("td_label.vtk", m);

  pc::denoiseShocksSerial(in, m, shocks, 30);

  apf::writeVtkFiles("td_denoise.vtk", m);

  // Cleanup
  m->destroyNative();
  apf::destroyMesh(m);

  PCU_Comm_Free();
  MPI_Finalize();

  return 0;
}
