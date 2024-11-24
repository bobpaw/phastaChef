#include <exception>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <apf.h>
#include <apfSIM.h>
#include <apfShape.h>
#include <gmi.h>
#include <gmi_sim.h>
#include <lionPrint.h>
#include <PCU.h>
#include <pcu_util.h>

#include <SimPartitionedMesh.h>

#include "../pcShockParam.h"

namespace {
  std::vector<apf::Vector3> readCSVPoints(std::string csvFile);
}

int main(int argc, char* argv[]) {
  // Initialize parallelism.
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  // Enable SCOREC library output.
  lion_set_verbosity(1);
  // Check arguments.
  if (argc != 4) {
    std::cout << "USAGE: " << argv[0] << " <sim_model.smd> <sim_mesh.sms> "
      "<shock_lcs.csv>" << std::endl;
    PCU_Comm_Free();
    MPI_Finalize();
    return 1;
  }
  // Initialize Simmetrix.
  gmi_register_sim();
  gmi_sim_start();
  // Read geometry.
  gmi_model* gmodel = gmi_load(argv[1]);
  // Read mesh.
  pParMesh sim_mesh = PM_load(argv[2], gmi_export_sim(gmodel), nullptr);
  // Create mesh wrapper.
  apf::Mesh2* apf_mesh = apf::createMesh(sim_mesh);
  apf::printStats(apf_mesh);
  // Create Shock_Param field.
  apf::Field* Shock_Param = apf::createField(apf_mesh, "Shock_Param",
    apf::SCALAR, apf::getConstant(3));
  // Read CSV file.
  std::vector<apf::Vector3> shock_pts = readCSVPoints(argv[3]);
  // Map CSV data to apf::MeshEntity*.
  double tol = 10e-3;
  apf::MeshIterator* it = apf_mesh->begin(3);
  for (apf::MeshEntity* e; e = apf_mesh->iterate(it);) {
    apf::Vector3 lc = apf::getLinearCentroid(apf_mesh, e);
    lc = lc * 1000;
    double sp = pc::ShockParam::NONE;
    for (size_t i = 0; i < shock_pts.size(); ++i) {
      apf::Vector3 d = shock_pts[i] - lc;
      if (d * d < tol * tol) {
        sp = pc::ShockParam::SHOCK;
        break;
      }
    }
    apf::setScalar(Shock_Param, e, 0, sp);
  }
  // Write VTK file.
  apf::writeVtkFiles("load_cev.vtk", apf_mesh);
  // Cleanup.
  apf::destroyMesh(apf_mesh);
  M_release(sim_mesh);
  gmi_destroy(gmodel);
  gmi_sim_stop();
  PCU_Comm_Free();
  MPI_Finalize();
  return 0;
}

namespace {
  std::vector<std::string> tokenize(std::string line) {
    std::vector<std::string> tokens;
    for (size_t n = 0, end = line.find(','); n < line.size(); n = end + 1) {
      tokens.push_back(line.substr(n, end - n));
      end = line.find(',', n);
      if (end == std::string::npos) break;
    }
    return tokens;
  }

  class CSVException : public std::exception {
  public:
    CSVException(int line, std::string m) {
      msg = "line " + std::to_string(line) + ": " + m;
    }
    const char* what(void) const noexcept {
      return msg.c_str();
    }
  private:
    std::string msg;
  };

  std::vector<apf::Vector3> readCSVPoints(std::string csvFile) {
    std::vector<apf::Vector3> points;
    std::ifstream str(csvFile);
    std::string line;
    std::getline(str, line); // skip header
    for (int lineno = 2; std::getline(str, line); ++lineno) {
      // Remove newline from line.
      line.erase(line.find_first_of("\r\n"));
      std::vector<std::string> fields = tokenize(line);
      if (fields.size() != 3)
        throw CSVException(lineno, "wrong number of fields");
      try {
        apf::Vector3 pt(std::stod(fields[0]), std::stod(fields[1]),
          std::stod(fields[2]));
        points.push_back(pt);
      } catch (const std::exception& err) {
        throw CSVException(lineno, err.what());
      }
    }
    return points;
  }
} // namespace
