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
  void rayTrace(apf::Mesh2* m, apf::MeshEntity* start, apf::MeshEntity* end, std::function<void(apf::Mesh2* m, apf::MeshEntity* e)>);
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

	return m;
}

void testRay(apf::Mesh2* m) {
	apf::Numbering* tet_id = m->findNumbering("tet id");

	apf::MeshEntity* tet1, *tet2;
	apf::MeshIterator* it = m->begin(3);
	for (apf::MeshEntity* e = m->iterate(it); e; e = m->iterate(it)) {
		int id = apf::getNumber(tet_id, e, 0, 0);
		if (id == 2) tet1 = e;
		if (id == 5) tet2 = e;
	}
	m->end(it);
	
	apf::Field* ray = apf::createField(m, "ray", apf::SCALAR, apf::getConstant(3));
  apf::zeroField(ray);

  int step = 1;
  pc::rayTrace(m, tet1, tet2, [tet_id, ray, &step](apf::Mesh2* m, apf::MeshEntity* e) {
    apf::setScalar(ray, e, 0, step);
    switch (step) {
    case 1:
      PCU_ALWAYS_ASSERT(apf::getNumber(tet_id, e, 0, 0) == 2);
      break;
    case 2:
      PCU_ALWAYS_ASSERT(apf::getNumber(tet_id, e, 0, 0) == 0);
      break;
    case 3:
      PCU_ALWAYS_ASSERT(apf::getNumber(tet_id, e, 0, 0) == 1);
      break;
    case 4:
      PCU_ALWAYS_ASSERT(apf::getNumber(tet_id, e, 0, 0) == 5);
      break;
    default:
      fprintf(stderr, "Error: too many steps. Infinite loop?\n");
      exit(EXIT_FAILURE);
    }
    ++step;
  });
}

int main (int argc, char* argv[]) {
	MPI_Init(&argc, &argv);
	PCU_Comm_Init();
	lion_set_verbosity(1);
	apf::Mesh2* m = createMesh();

	testRay(m);

	apf::writeVtkFiles("tets.vtk", m, 3);
	apf::writeVtkFiles("tris.vtk", m, 2);
	gmi_write_dmg(m->getModel(), "tets.dmg");
	m->writeNative("tets.smb");

	m->destroyNative();
	apf::destroyMesh(m);
	PCU_Comm_Free();
	MPI_Finalize();
}
