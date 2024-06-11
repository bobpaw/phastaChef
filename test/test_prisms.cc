#include <iostream>
#include <PCU.h>
#include <lionPrint.h>
#include <gmi_null.h>
#include <gmi_mesh.h>
#include <apfMDS.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfConvert.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <pcu_util.h>

namespace pc {
  // Private test target.
  void rayTrace(apf::Mesh* m, apf::MeshEntity* start, apf::MeshEntity* end, std::function<void(apf::Mesh* m, apf::MeshEntity* e)>);
}

apf::Mesh2* createMesh(void) {
  gmi_register_null();
  gmi_register_mesh();
  gmi_model* g = gmi_load(".null");
  apf::Mesh2* m = apf::makeEmptyMdsMesh(g, 3, false);

  const int prisms = 4 * 2, tets = prisms * 3, planeverts = 6,
    verts = planeverts * 3;

  double coordsZ0[planeverts * 2] = {
    0, 0, 1, 0, 0, 1, 1, 1, 2, 0, 2, 1
  };
  double coords[verts * 3];
  for (int i = 0; i < planeverts; ++i) {
    coords[i * 3] = coordsZ0[i * 2]; // x component
    coords[i * 3 + 1] = coordsZ0[i * 2 + 1]; // y component
    coords[i * 3 + 2] = 0;
    coords[(i + planeverts) * 3] = coordsZ0[i * 2]; // x component
    coords[(i + planeverts) * 3 + 1] = coordsZ0[i * 2 + 1]; // y component
    coords[(i + planeverts) * 3 + 2] = 1;
    coords[(i + planeverts * 2) * 3] = coordsZ0[i * 2]; // x component
    coords[(i + planeverts * 2) * 3 + 1] = coordsZ0[i * 2 + 1]; // y component
    coords[(i + planeverts * 2) * 3 + 2] = -1;
  }
  apf::Gid conn[tets * 4] = {
    0,6,7,8,/**/ 0,1,8,7,/**/ 0,1,2,8, // Prism 1
    1,3,8,9,/**/ 1,3,2,8,/**/ 1,8,7,9, // Prism 2
    1,7,10,9,/**/ 1,4,9,10,/**/ 1,4,3,9, // Prism 3
    4,5,9,11,/**/ 4,5,3,9,/**/ 4,9,10,11, // Prism 4
    12,0,1,2,/**/ 12,13,2,1,/**/ 12,13,14,2, // Prism 5
    13,15,2,3,/**/ 13,15,14,2,/**/ 13,2,1,3, // Prism 6
    13,1,4,3,/**/ 13,16,3,4,/**/ 13,16,15,3, // Prism 7
    16,17,3,5,/**/ 16,17,15,3,/**/ 16,3,4,5 // Prism 8
  };

  apf::GlobalToVert outmap;
  apf::NewElements elms = apf::assemble(m, conn, tets, apf::Mesh::TET, outmap);
  apf::deriveMdsModel(m);
  apf::finalise(m, outmap);
  apf::setCoords(m, coords, verts, outmap);
  m->verify();

  apf::numberOwnedDimension(m, "tet id", 3);

  apf::MeshIterator* it = m->begin(0);
  for (apf::MeshEntity* e = m->iterate(it); e; e = m->iterate(it)) {
    apf::Vector3 pt;
    m->getPoint(e, 0, pt);
    std::cout << "Vert " << e << ": " << pt << std::endl;
  }
  m->end(it);

  return m;
}

std::vector<apf::MeshEntity*> tetsById;

void testRay(apf::Mesh* m, const char* fieldname, int start_id, int end_id, const std::vector<int>& steps) {
  apf::Field* ray = apf::createField(m, fieldname, apf::SCALAR, apf::getConstant(3));
  apf::zeroField(ray);
  int step = 0;
  pc::rayTrace(m, tetsById[start_id], tetsById[end_id], [ray, &step, &steps](apf::Mesh* m, apf::MeshEntity* e) {
    apf::setScalar(ray, e, 0, step + 1);
    if (step < steps.size()) {
      PCU_ALWAYS_ASSERT(e == tetsById[steps[step]]);
    } else {
      fprintf(stderr, "Error: too many steps. Infinite loop?\n");
      exit(EXIT_FAILURE);
    }
    ++step;
  });
  if (step != steps.size()) {
    std::cerr << "ERROR: Ray " << fieldname << " stopped at step " << step << std::endl;
    exit(1);
  }
}

void testRays(apf::Mesh2* m) {
  apf::Numbering* tet_id = m->findNumbering("tet id");

  tetsById.resize(m->count(3));
  apf::MeshIterator* it = m->begin(3);
  for (apf::MeshEntity* e = m->iterate(it); e; e = m->iterate(it)) {
    int id = apf::getNumber(tet_id, e, 0, 0);
    tetsById[id] = e;
  }
  m->end(it);

  testRay(m, "ray1", 10, 8, {10, 8});
  testRay(m, "ray2", 10, 0, {10, 8, 3, 5, 1, 0});
}

int main (int argc, char* argv[]) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  lion_set_verbosity(1);
  apf::Mesh2* m = createMesh();
  apf::writeVtkFiles("test_prisms_pre.vtk", m, 3);
  testRays(m);
  apf::writeVtkFiles("test_prisms_post.vtk", m, 3);
  m->destroyNative();
  apf::destroyMesh(m);
  PCU_Comm_Free();
  MPI_Finalize();
}
