#include "chefPhasta.h"
#include <PCU.h>
#include <chef.h>
#include <phasta.h>
#include <phstream.h>
#include <sam.h>
#include <apfMDS.h>
#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>

namespace {
  void freeMesh(apf::Mesh* m) {
    m->destroyNative();
    apf::destroyMesh(m);
  }

  void mychdir(const char* path) {
    const int fail = chdir(path);
    if(fail) {
      fprintf(stderr, "ERROR failed to change to %s dir... exiting\n", path);
      exit(1);
    }
  }

  void mychdir(int step) {
    std::stringstream path;
    path << step;
    string s = path.str();
    const int fail = chdir(s.c_str());
    if(fail) {
      fprintf(stderr, "ERROR failed to change to %d dir... exiting\n", step);
      exit(1);
    }
  }

  std::string makeMeshName(int step) {
    std::stringstream meshname;
    meshname  << "bz2:" << "t" << step << "p" << PCU_Comm_Peers() << "_.smb";
    return meshname.str();
  }

  std::string makeRestartName() {
    std::stringstream restartname;
    restartname << PCU_Comm_Peers() << "-procs_case/restart";
    return restartname.str();
  }

  std::string prefixCwd(std::string name) {
    char cwd[4096] = "\0";
    getcwd(cwd,4096);
    std::stringstream s;
    s << cwd << "/" << name;
    return s.str();
  }

  apf::Field* getField(apf::Mesh* m) {
    /* if the value of the fldIdx'th index from the fldName
     * field is greater than fldLimit then multiply the current
     * isotropic mesh size at the vertex by szFactor */
    const unsigned fldIdx = 5;
    const double fldLimit = 1e-6;
    const double szFactor = 0.5;
    const char* fldName = "errors";
    return sam::errorThreshold(m,fldName,fldIdx,fldLimit,szFactor);
  }

  static FILE* openfile_read(ph::Input&, const char* path) {
    return fopen(path, "r");
  }

  static FILE* openstream_read(ph::Input& in, const char* path) {
    std::string fname(path);
    std::string restartStr("restart");
    FILE* f = NULL;
    if( fname.find(restartStr) != std::string::npos )
      f = openRStreamRead(in.rs);
    else {
      fprintf(stderr,
        "ERROR %s type of stream %s is unknown... exiting\n",
        __func__, fname.c_str());
      exit(1);
    }
    return f;
  }

  void setupChef(ph::Input& ctrl, int step) {
    //don't split or tetrahedronize
    ctrl.splitFactor = 1;
    ctrl.tetrahedronize = 0;
    ctrl.timeStepNumber = step;
    ctrl.solutionMigration = 1;
    if(step>1) {
      if(!PCU_Comm_Self()) {
        fprintf(stderr, "STATUS error based adapt %d\n", step);
        fprintf(stderr, "STATUS ctrl.attributeFileName %s step %d\n",
            ctrl.attributeFileName.c_str(), step);
      }
      ctrl.adaptStrategy = 1; //error field adapt
      ctrl.adaptFlag = 1;
    }
    ctrl.outMeshFileName = makeMeshName(step);
    ctrl.restartFileName = makeRestartName();
  }
}
int main(int argc, char** argv) {
  MPI_Init(&argc, &argv);
  PCU_Comm_Init();
  PCU_Protect();
  if( argc != 6 ) {
    if(!PCU_Comm_Self())
      fprintf(stderr, "Usage: %s <maxTimeStep> <ramdisk path> <solver.inp> <input.config> <adapt.inp>\n",argv[0]);
    exit(EXIT_FAILURE);
  }
  int maxStep = atoi(argv[1]);
  const char* ramdisk = argv[2];
  const char* solverinp = argv[3];
  const char* inputcfg = argv[4];
  const char* adaptinp = argv[5];
  gmi_model* g = 0;
  apf::Mesh2* m = 0;

  mychdir(ramdisk);
  ph::Input ctrl;
  ctrl.load(adaptinp);
  ctrl.outMeshFileName = makeMeshName(0);
  chefPhasta::initModelers();
  chef::cook(g,m,ctrl);
  freeMesh(m); m = NULL;
  phSolver::Input inp(solverinp,inputcfg);
  int step = 0;
  do {
    ctrl.meshFileName = makeMeshName(step);
    step = phasta(inp);
    assert(step >= 0);
    if(!PCU_Comm_Self())
      fprintf(stderr, "STATUS ran to step %d\n", step);
    if( step >= maxStep )
      break;
    setupChef(ctrl,step);
    chef::readAndAttachFields(ctrl,m);
    apf::Field* szFld = getField(m);
    assert(szFld);
    chef::adapt(m,szFld);
    apf::destroyField(szFld);
    chef::balanceAndReorder(ctrl,m);
    chef::preprocess(m,ctrl);
    freeMesh(m); m = NULL;
    mychdir(step);
  } while( step < maxStep );
  chefPhasta::finalizeModelers();
  PCU_Comm_Free();
  MPI_Finalize();
}