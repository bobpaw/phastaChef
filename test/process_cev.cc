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
#include <phInput.h>

#include <SimAdvMeshing.h>
#include <SimPartitionedMesh.h>

#include "csv.h"

namespace pc {
  void processShocksSerial(ph::Input& in, apf::Mesh2* m);
}

int main(int argc, char* argv[]) {
  // Initialize parallelism.
  SimPartitionedMesh_start(&argc, &argv); // MPI_Init()
  PCU_Comm_Init();
  // Check arguments.
  if (argc != 4) {
    std::cout << "USAGE: " << argv[0] << " <nat_model.x_t> <sim_model.smd> "
      "<sim_mesh.sms>" << std::endl;
    PCU_Comm_Free();
    MPI_Finalize();
    return 1;
  }
  // Enable SCOREC library output.
  lion_set_verbosity(1);
  Sim_logOn("sim.log");
  // Initialize Simmetrix.
  MS_init();
  SimAdvMeshing_start();
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
  // Load Shock_Param.
  apf::Field* Shock_Param = apf::createField(apf_mesh, "Shock_Param",
    apf::SCALAR, apf::getConstant(3));
  std::vector<test::FieldRow> rows =
    test::readCSVField("Shock_Param.csv", "Shock_Param");
  std::cout << "CSV file contained " << rows.size() << " rows for "
    << apf_mesh->count(3) << " elements." << std::endl;
  size_t hits = 0, i = 0;
  constexpr double tol = 1e-3;
  apf::MeshIterator* it = apf_mesh->begin(3);
  for (apf::MeshEntity* e; e = apf_mesh->iterate(it); ++i) {
    PCU_DEBUG_ASSERT(i < rows.size());
    // Confirm LC.
    apf::Vector3 lc = apf::getLinearCentroid(apf_mesh, e);
    apf::Vector3 dist = lc - rows[i].lc;
    if (dist * dist < tol * tol) {
      ++hits;
      apf::setScalar(Shock_Param, e, 0, rows[i].val);
    }
  }
  apf_mesh->end(it);
  std::cout << "CSV mapped " << hits << '/' << rows.size() << std::endl;
  // Make dummy ph::Input.
  ph::Input in;
  // Process shock.
  pc::processShocksSerial(in, apf_mesh);
  // Cleanup.
  apf::destroyMesh(apf_mesh);
  M_release(sim_mesh);
  gmi_destroy(gmodel);
  SimAdvMeshing_stop();
  gmi_sim_stop();
  Sim_unregisterAllKeys();
  MS_exit();
  Sim_logOff();
  PCU_Comm_Free();
  SimPartitionedMesh_stop(); // MPI_Finalize()
  return 0;
}
