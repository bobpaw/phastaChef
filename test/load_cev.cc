#include <algorithm>
#include <exception>
#include <fstream>
#include <iostream>
#include <queue>
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
#include "csv.h"

namespace pc {
  void processShocksSerial(ph::Input& in, apf::Mesh2* m);
}

int main(int argc, char* argv[]) {
  // Initialize parallelism.
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  // Enable SCOREC library output.
  lion_set_verbosity(1);
  // Check arguments.
  if (argc != 6) {
    std::cout << "USAGE: " << argv[0] << " <nat_model.x_t> <sim_model.smd> "
      "<sim_mesh.sms> <shock_lcs.csv> <out_mesh.sms>" << std::endl;
    PCU_Comm_Free();
    MPI_Finalize();
    return 1;
  }
  // Initialize Simmetrix.
  gmi_register_sim();
  gmi_sim_start();
  Sim_readLicenseFile(0);
  // Read geometry.
  gmi_model* gmodel = gmi_sim_load(argv[1], argv[2]);
  // Read mesh.
  pParMesh sim_mesh = PM_load(argv[3], gmi_export_sim(gmodel), nullptr);
  // Create mesh wrapper.
  apf::Mesh2* apf_mesh = apf::createMesh(sim_mesh);
  apf::printStats(apf_mesh);
  // Create Shock_Param field.
  apf::Field* Shock_Param = apf::createField(apf_mesh, "Shock_Param",
    apf::SCALAR, apf::getConstant(3));
  // Read CSV file.
  std::vector<apf::Vector3> shock_pts = test::readCSVPoints(argv[4]);
  // Map CSV data to apf::MeshEntity*.
  size_t hits = 0;
  double tol = 10e-3;
  apf::MeshIterator* it = apf_mesh->begin(3);
  for (apf::MeshEntity* e; e = apf_mesh->iterate(it);) {
    apf::Vector3 lc = apf::getLinearCentroid(apf_mesh, e);
    lc = lc * 1000;
    double sp = pc::ShockParam::NONE;
    for (size_t i = 0; i < shock_pts.size(); ++i) {
      apf::Vector3 d = shock_pts[i] - lc;
      if (d * d < tol * tol) {
        ++hits;
        sp = pc::ShockParam::SHOCK;
        break;
      }
    }
    apf::setScalar(Shock_Param, e, 0, sp);
  }
  apf_mesh->end(it);
  std::cout << "Matched " << hits << '/' << shock_pts.size() << " centroids"
    << std::endl;
  // Write field.
  std::vector<test::FieldRow> rows;
  it = apf_mesh->begin(3);
  for (apf::MeshEntity* e; e = apf_mesh->iterate(it);) {
    double sp = apf::getScalar(Shock_Param, e, 0);
    apf::Vector3 lc = apf::getLinearCentroid(apf_mesh, e);
    rows.emplace_back(sp, lc);
  }
  apf_mesh->end(it);
  test::writeCSVField("Shock_Param.csv", "Shock_Param", rows);
  // Write VTK file.
  apf::writeVtkFiles("load_cev.vtk", apf_mesh);
  // Cleanup.
  apf::destroyMesh(apf_mesh);
  M_release(sim_mesh);
  gmi_destroy(gmodel);
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  PCU_Comm_Free();
  MPI_Finalize();
  return 0;
}
