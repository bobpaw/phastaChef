#include <cassert>
#include <cstring>
#include <random>
#include <iostream>
#include <fstream>
#include <chrono>
#include <vector>

#include <apfMesh2.h>
#include <apfShape.h>
#include <apfNumbering.h>
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
}

int mockShockParam(apf::Mesh* m, pc::_test::ShockFunc* sf) {
  assert(m);
  assert(sf);

  apf::Field* Shock_Param = m->findField("Shock_Param");
  if (Shock_Param == nullptr) {
    Shock_Param = apf::createField(m, "Shock_Param", apf::SCALAR,
      apf::getConstant(3));
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

double calculateTolerance(double alpha, apf::Mesh2* m) {
  return alpha * ma::getAverageEdgeLength(m);
}

class TestDefragArgs {
public:
  static void print_usage(std::ostream& str, char* argv0) {
    str << "USAGE: " << argv0 << " -g <modelFileName> -m <meshFileName>"
        << "-a <tol_alpha> -f <fragmentation> -s type:args,... [-d] [-r]"
        << std::endl << std::endl;
    str << "Multiple shock functions may be specified." << std::endl;
    str << "-d\tAssign element IDs and dump info to 'mesh_elms.txt'."
        << std::endl;
    str << "-r\tDo initial uniform refinement." << std::endl;
  }

  TestDefragArgs() {}

  void parse(int argc, char* argv[]) {
    bool g_given = false, m_given = false, a_given = false, f_given = false,
      s_given = false;
    for (int i = 1; i < argc; ++i) {
      if (strcmp(argv[i], "-g") == 0) {
        ++i;
        if (i < argc) {
          model_ = argv[i];
          g_given = true;
        } else {
          throw std::invalid_argument("Expected argument after '-g'.");
        }
      } else if (strcmp(argv[i], "-m") == 0) {
        ++i;
        if (i < argc) {
          mesh_ = argv[i];
          m_given = true;
        } else {
          throw std::invalid_argument("Expected argument after '-m'.");
        }
      } else if (strcmp(argv[i], "-a") == 0) {
        ++i;
        if (i < argc) {
          alpha_ = std::atof(argv[i]);
          a_given = true;
        } else {
          throw std::invalid_argument("Expected argument after '-a'.");
        }
      } else if (strcmp(argv[i], "-f") == 0) {
        ++i;
        if (i < argc) {
          fragmentation_ = std::atoi(argv[i]);
          f_given = true;
        } else {
          throw std::invalid_argument("Expected argument after '-f'.");
        }
      } else if (strcmp(argv[i], "-s") == 0) {
        ++i;
        if (i < argc) {
          shockfuncs_.push_back(argv[i]);
          s_given = true;
        } else {
          throw std::invalid_argument("Expected argument after '-s'.");
        }
      } else if (strcmp(argv[i], "-d") == 0) {
        dump_ = true;
      } else if (strcmp(argv[i], "-r") == 0) {
        refine_ = true;
      } else {
        throw std::invalid_argument(std::string("Unknown argument '") + argv[i]
          + "'.");
      }
    }

    if (!g_given) throw std::invalid_argument("Missing required argument -g.");
    if (!m_given) throw std::invalid_argument("Missing required argument -m.");
    if (!a_given) throw std::invalid_argument("Missing required argument -a.");
    if (!f_given) throw std::invalid_argument("Missing required argument -f.");
    if (!s_given) throw std::invalid_argument("Missing required argument -s.");
  }

  int fragmentation(void) const noexcept { return fragmentation_; }
  double alpha(void) const noexcept { return alpha_; }
  bool dump(void) const noexcept { return dump_; }
  bool refine(void) const noexcept { return refine_; }
  const std::string& mesh(void) const noexcept { return mesh_; }
  const std::string& model(void) const noexcept { return model_; }
  const std::vector<std::string>& shockfuncs(void) const noexcept {
    return shockfuncs_;
  }
  const std::string& shockfunc(int i) const { return shockfuncs_[i]; }

private:
  int fragmentation_{0};
  double alpha_{0.0};
  bool dump_{false}, refine_{false};
  std::string mesh_, model_;
  std::vector<std::string> shockfuncs_;
};

int main (int argc, char* argv[]) {
  // MPI/PCU init?
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();

  // lion verbosity
  lion_set_verbosity(1);

  TestDefragArgs args;

  try {
    args.parse(argc, argv);
  } catch (const std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    args.print_usage(std::cerr, argv[0]);
    return -1;
  }

  apf::Mesh2* m = nullptr;

  ph::Input in;
  in.modelFileName = args.model();
  in.meshFileName = args.mesh();

  // Init gmi.
  gmi_register_mesh();

  // Load mesh.
  m = apf::loadMdsMesh(in.modelFileName.c_str(), in.meshFileName.c_str());

  if (args.refine()) {
    ma::runUniformRefinement(m);
    m->verify();
  }

  if (args.dump()) {
    apf::Numbering *elm_id = apf::numberOwnedDimension(m, "elm_id", 3);
    apf::writeVtkFiles("mesh.vtk", m);
    std::ofstream outfile("mesh_elms.txt");
    outfile<<"id,lc,handle"<<std::endl;
    apf::MeshIterator* it = m->begin(3);
    for (apf::MeshEntity* e = m->iterate(it); e; e = m->iterate(it)) {
      int id = apf::getNumber(elm_id, e, 0, 0);
      apf::Vector3 lc = apf::getLinearCentroid(m,e);
      outfile << id << ',' << lc << ',' << e << std::endl;
    }
    m->end(it);
  }

  double tol = calculateTolerance(args.alpha(), m);
  std::cout << "Tolerance: " << tol << std::endl;

  std::vector<pc::_test::ShockFunc*> shockfuncs;
  try {
    for (int i = 0; i < args.shockfuncs().size(); ++i) {
      pc::_test::ShockFunc* sf = pc::_test::ShockFunc::makeFromString(
        args.shockfunc(i), tol);
      shockfuncs.push_back(sf);
    }
  } catch (std::exception& e) {
    std::cerr << "ERROR: " << e.what() << std::endl;
    args.print_usage(std::cerr, argv[0]);
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
  int n = 0;
  for (int i = 0; i < shockfuncs.size(); ++i) {
    n += mockShockParam(m, shockfuncs[i]);
  }
  auto clock_end = std::chrono::steady_clock::now();
  std::cout << "Assigned " << n << " shock elements: " << my_tdiff(clock_start, clock_end) << " seconds." << std::endl;

  // Write full shock VTK file.
  apf::writeVtkFiles("td_shock.vtk", m);

  // Fragmentation step
  clock_start = std::chrono::steady_clock::now();
  n = fragmentShockParam(m, args.fragmentation());
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
