#include <functional>
#include <iostream>
#include <PCU.h>
#include <lionPrint.h>
#include <gmi_null.h>
#include <apfMDS.h>
#include <apf.h>
#include <apfMesh2.h>
#include <apfConvert.h>
#include <apfNumbering.h>
#include <apfShape.h>
#include <pcu_util.h>

namespace pc {
  // Private test target.
  void rayTrace(apf::Mesh* m, apf::MeshEntity* start, apf::MeshEntity* end, std::function<bool(apf::Mesh* m, apf::MeshEntity* e)>);
}

apf::Mesh2* createMesh(void) {
	gmi_register_null();
	gmi_model* g = gmi_load(".null");
	apf::Mesh2* m = apf::makeEmptyMdsMesh(g, 3, false);

  const int tets = 6, verts = 8;

	apf::Gid conn[tets * 4] = {0, 2, 3, 1, /**/ 0, 2, 4, 3, /**/ 1, 5, 2, 3,
  4, 2, 0, 6, /**/ 4, 0, 3, 7, /**/ 4, 0, 7, 6};
	double coords[verts * 3] = {0, 1, 1,
		0, -1, 1,
		-1, 0, 0,
		1, 0, 0,
		0, 1, 0,
		0, -1, 0,
    -1, 1.5, 0.5,
    1, 1.5, 0.5};

	apf::GlobalToVert outmap;
	apf::NewElements elms = apf::assemble(m, conn, tets, apf::Mesh::TET, outmap);
	apf::deriveMdsModel(m);
	apf::setCoords(m, coords, verts, outmap);
	apf::finalise(m, outmap);
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

namespace pc {
  int getFaceVertsCCW(apf::Mesh* m, apf::MeshEntity* e, int face, apf::MeshEntity** verts);
  int edgeIndicator(apf::Mesh* m, const apf::Vector3& src, const apf::Vector3& ray,
                apf::MeshEntity* v1, apf::MeshEntity* v2);
}

void testReorder(apf::Mesh* m) {
  apf::MeshIterator* it = m->begin(3);
  for (apf::MeshEntity* tet = m->iterate(it); tet; tet = m->iterate(it)) {
    apf::Vector3 lc = apf::getLinearCentroid(m, tet);
    apf::Downward faces;
    int face_ct = m->getDownward(tet, 2, faces);
    for (int i = 0; i < face_ct; ++i) {
      apf::Vector3 face_ray = apf::getLinearCentroid(m, faces[i]) - lc;
      apf::Downward verts;
      int vert_ct = pc::getFaceVertsCCW(m, tet, i, verts); //m->getDownward(faces[i], 0, verts);
      verts[vert_ct] = verts[0];
      for (int j = 0; j < vert_ct; ++j) {
        int ind = pc::edgeIndicator(m, lc, face_ray, verts[j], verts[j + 1]);
        if (ind != 1) {
          fprintf(stderr, "ERROR: Tet %p face %d (%p,%p,%p) "
                  "does not point outward!\n", tet, i, verts[0], verts[1], verts[2]);
          exit(EXIT_FAILURE);
        }
      }
    }
  }
  m->end(it);
}

std::vector<apf::MeshEntity*> tetsById;

void testRay(apf::Mesh* m, const char* fieldname, int start_id, int end_id, const std::vector<int>& steps) {
  std::cout << "Testing ray " << fieldname << '.' << std::endl;
  apf::Field* ray = apf::createField(m, fieldname, apf::SCALAR, apf::getConstant(3));
  apf::zeroField(ray);
  int step = 0;
  pc::rayTrace(m, tetsById[start_id], tetsById[end_id], [ray, &step, &steps](apf::Mesh* m, apf::MeshEntity* e) {
    apf::setScalar(ray, e, 0, step + 1);
    if (step < steps.size()) {
      if (e != tetsById[steps[step]]) {
        std::cerr << "ERROR: Expected to be at element "
                  << tetsById[steps[step]] << " but actually at element "
                  << e << std::endl;
        exit(EXIT_FAILURE);
      }
    } else {
      std::cerr << "ERROR: Too many steps. Infinite loop?" << std::endl;
      exit(EXIT_FAILURE);
    }
    ++step;
    return true;
  });
  if (step != steps.size()) {
    std::cerr << "ERROR: Ray " << fieldname << " stopped at step " << step << std::endl;
    exit(EXIT_FAILURE);
  }
}

void testRays(apf::Mesh2* m) {
	apf::Numbering* tet_id = m->findNumbering("tet id");

	tetsById.resize(m->count(3));
	apf::MeshIterator* it = m->begin(3);
	for (apf::MeshEntity* e = m->iterate(it); e; e = m->iterate(it)) {
		int id = apf::getNumber(tet_id, e, 0, 0);
    tetsById[id] = e;
    std::cout << "Tet " << id << ": " << e << ": " <<
      apf::getLinearCentroid(m, e) << std::endl;
	}
	m->end(it);

  testRay(m, "ray1", 2, 5, {2, 0, 1, 5});
  testRay(m, "ray2", 2, 0, {2, 0});
  testRay(m, "ray3", 0, 1, {0, 1});
  testRay(m, "ray4", 1, 3, {1, 3});
  testRay(m, "ray5", 1, 4, {1, 4});
  testRay(m, "ray6", 1, 0, {1, 0});
  testRay(m, "ray7", 3, 1, {3, 1});
  testRay(m, "ray8", 4, 1, {4, 1});
  testRay(m, "ray9", 0, 1, {0, 1});
}

int main (int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	PCU_Comm_Init();
	lion_set_verbosity(1);
	apf::Mesh2* m = createMesh();

  apf::writeVtkFiles("test_raytracing_pre.vtk", m, 3);

  testReorder(m);
	testRays(m);

	apf::writeVtkFiles("test_raytracing_tets.vtk", m, 3);
	gmi_write_dmg(m->getModel(), "test_raytracing_tets.dmg");
	m->writeNative("test_raytracing_tets.smb");

	m->destroyNative();
	apf::destroyMesh(m);
	PCU_Comm_Free();
	MPI_Finalize();
}
