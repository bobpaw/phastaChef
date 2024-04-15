#include "pcAdapter.h"
#include "pcShock.h"
#include "pcError.h"
#include "pcUpdateMesh.h"
#include "pcSmooth.h"
#include "pcWriteFiles.h"
#include <SimUtil.h>
#include <SimPartitionedMesh.h>
#include <SimDiscrete.h>
#include "SimMeshMove.h"
#include "SimMeshTools.h"
#include <SimAdvMeshing.h>
#include "apfSIM.h"
#include "gmi_sim.h"
#include <PCU.h>
#include <cassert>
#include <phastaChef.h>
#include <maStats.h>
#include <apfShape.h>
#include <math.h>
#include <fstream>
#include <limits>
#include <iomanip>
#include <set>
#include <map>
#include <queue>
#include <limits>
#include <algorithm>
#include <random>
#include <spr.h>
#include <ctime>

// // for shock surface ID
// #include <armadillo>
// #include <unordered_set>
#include <vector>

// typedef unordered_set<int>::iterator nbd_iter;
// typedef unordered_set<int>& set_ref;

extern void MSA_setBLSnapping(pMSAdapt, int onoff);

namespace pc {

     std::vector<std::string> ParseDeliminatedString(std::string& FullLine, std::string& delim){
          std::vector<std::string> output;
          size_t pos=0;
          std::string token;

          while ((pos = FullLine.find(delim)) != std::string::npos) {
               token = FullLine.substr(0, pos);
               output.push_back(token);
               FullLine.erase(0, pos + delim.length());
          }
          output.push_back(token);


          return output;          
     }

     std::set<apf::MeshEntity*> TraceSurf(std::map<apf::MeshEntity*, std::pair<std::set<apf::MeshEntity*>,int> >& LocGraph,
                                        std::set<apf::MeshEntity*>& StartSet){
          
          std::set<apf::MeshEntity*> output=StartSet;
          for (auto itr1=StartSet.begin(); itr1!=StartSet.end(); itr1++){
               for (auto nei=LocGraph[*itr1].first.begin(); nei!=LocGraph[*itr1].first.end(); nei++){
                    output.insert(*nei); 
               }     
          }
     return output;
     }

  bool vertexIsInCylinder(apf::MeshEntity* v) {
    double x_min = 0.0;
    double x_max = 2.01;
    double r_max = 0.07;//outer radius

    apf::Vector3 xyz = apf::Vector3(0.0,0.0,0.0);
    m->getPoint(v,0,xyz);
    double xyz_r = sqrt(xyz[1]*xyz[1] + xyz[2]*xyz[2]);
    if (xyz[0] > x_min && xyz[0] < x_max && xyz_r < r_max)
      return true;
    return false;
  }


  apf::Field* convertField(apf::Mesh* m,
    const char* inFieldname,
    const char* outFieldname) {
    apf::Field* inf = m->findField(inFieldname);
    assert(inf);
    int size = apf::countComponents(inf);
    apf::Field* outf = m->findField(outFieldname);
    if (outf)
      apf::destroyField(outf);
    outf = apf::createPackedField(m, outFieldname, size);
    apf::NewArray<double> inVal(size);
    apf::NewArray<double> outVal(size);
    apf::MeshEntity* vtx;
    apf::MeshIterator* it = m->begin(0);
    while ((vtx = m->iterate(it))) {
      apf::getComponents(inf, vtx, 0, &inVal[0]);
      for (int i = 0; i < size; i++){
        outVal[i] = inVal[i];
      }
      apf::setComponents(outf,vtx, 0, &outVal[0]);
    }
    m->end(it);
    apf::destroyField(inf);
    return outf;
  }

  apf::Field* convertVtxFieldToElm(apf::Mesh* m,
                  const char* vFieldname,
                  const char* eFieldname) {
    apf::Field* vf = m->findField(vFieldname);
    assert(vf);
    int size = apf::countComponents(vf);
    apf::Field* ef = m->findField(eFieldname);
    if (ef) apf::destroyField(ef);
    ef = apf::createPackedField(m, eFieldname, size, apf::getConstant(m->getDimension()));
    apf::NewArray<double> vVal(size);
    apf::NewArray<double> eVal(size);
    apf::MeshEntity* elm;
    apf::MeshIterator* it = m->begin(m->getDimension());
    while ((elm = m->iterate(it))) {
      for (int i = 0; i < size; i++) eVal[i] = 0.0;
      apf::Downward vtx;
      int nbv = m->getDownward(elm, 0, vtx);
      for (int j = 0; j < nbv; j++){
        apf::getComponents(vf, vtx[j], 0, &vVal[0]);
        for (int i = 0; i < size; i++){
          eVal[i] += vVal[i]/(double)nbv;
        }
      }
      apf::setComponents(ef, elm, 0, &eVal[0]);
    }
    m->end(it);
    apf::destroyField(vf);
    return ef;
  }

  void attachMeshSizeField(apf::Mesh2*& m, ph::Input& in, phSolver::Input& inp) {
    /* create a field to store mesh size */
    if(m->findField("sizes")) apf::destroyField(m->findField("sizes"));
    apf::Field* sizes = apf::createSIMFieldOn(m, "sizes", apf::VECTOR);
    /* switch between VMS error mesh size and initial mesh size */
    if((string)inp.GetValue("Error Estimation Option") != "False") {
      pc::attachVMSSizeField(m, in, inp);
    }
    else {
      if(m->findField("frames")) apf::destroyField(m->findField("frames"));
      apf::Field* frames = apf::createSIMFieldOn(m, "frames", apf::MATRIX);
      ph::attachSIMSizeField(m, sizes, frames);
    }
// prescribe mesh size field for the projectile case
// this is hardcoded, please comment out this call for other usage
//    pc::prescribe_proj_mesh_size(m, sizes, in.rbParamData[0]);
  }

  int getNumOfMappedFields(apf::Mesh2*& m) {
    /* initially, we have 7 fields: pressure, velocity, temperature,
       time der of pressure, time der of velocity, time der of temperature,
       ,mesh velocity and 1 optional field: time resource bound factor field */
    int numOfMappedFields;
    if (m->findField("ctcn_elm")) numOfMappedFields = 8;
    else numOfMappedFields = 7;
    return numOfMappedFields;
  }

  /* remove all fields except for solution, time
     derivative of solution, mesh velocity and keep
     certain field if corresponding option is on */
  void removeOtherFields(apf::Mesh2*& m, phSolver::Input& inp) {
    int index = 0;
    int numOfPackFields = 4;
    while (m->countFields() > numOfPackFields) {
      apf::Field* f = m->getField(index);
      if ( f == m->findField("solution") ||
           f == m->findField("time derivative of solution") ||
           f == m->findField("mesh_vel") ||
           f == m->findField("ctcn_elm") ) {
        index++;
        continue;
      }
      m->removeField(f);
      apf::destroyField(f);
    }
    m->verify();
  }

  int getSimFields(apf::Mesh2*& m, int simFlag, pField* sim_flds, phSolver::Input& inp) {
    int num_flds = 0;
    if (m->findField("solution")) {
      num_flds += 3;
      sim_flds[0] = apf::getSIMField(chef::extractField(m,"solution","pressure",1,apf::SCALAR,simFlag));
      sim_flds[1] = apf::getSIMField(chef::extractField(m,"solution","velocity",2,apf::VECTOR,simFlag));
      sim_flds[2] = apf::getSIMField(chef::extractField(m,"solution","temperature",5,apf::SCALAR,simFlag));
      apf::destroyField(m->findField("solution"));
    }

    if (m->findField("time derivative of solution")) {
      num_flds += 3;
      sim_flds[3] = apf::getSIMField(chef::extractField(m,"time derivative of solution","der_pressure",1,apf::SCALAR,simFlag));
      sim_flds[4] = apf::getSIMField(chef::extractField(m,"time derivative of solution","der_velocity",2,apf::VECTOR,simFlag));
      sim_flds[5] = apf::getSIMField(chef::extractField(m,"time derivative of solution","der_temperature",5,apf::SCALAR,simFlag));
      apf::destroyField(m->findField("time derivative of solution"));
    }

    if (m->findField("mesh_vel")) {
      num_flds += 1;
      sim_flds[6] = apf::getSIMField(chef::extractField(m,"mesh_vel","mesh_vel_sim",1,apf::VECTOR,simFlag));
      apf::destroyField(m->findField("mesh_vel"));
    }

    if (m->findField("ctcn_elm")) {
      num_flds += 1;
      sim_flds[7] = apf::getSIMField(chef::extractField(m,"ctcn_elm","ctcn_elm_sim",1,apf::SCALAR,simFlag));
      apf::destroyField(m->findField("ctcn_elm"));
    }

    return num_flds;
  }

  /* unpacked solution into serveral fields,
     put these field explicitly into pPList */
  pPList getSimFieldList(ph::Input& in, apf::Mesh2*& m){
    /* load input file for solver */
    phSolver::Input inp("solver.inp", "input.config");
    int num_flds = getNumOfMappedFields(m);
    removeOtherFields(m,inp);
    pField* sim_flds = new pField[num_flds];
    getSimFields(m, in.simmetrixMesh, sim_flds, inp);
    pPList sim_fld_lst = PList_new();
    for (int i = 0; i < num_flds; i++) {
      PList_append(sim_fld_lst, sim_flds[i]);
    }
    assert(num_flds == PList_size(sim_fld_lst));
    delete [] sim_flds;
    return sim_fld_lst;
  }

  void transferSimFields(apf::Mesh2*& m) {
    if (m->findField("pressure")) // assume we had solution before
      chef::combineField(m,"solution","pressure","velocity","temperature");
    if (m->findField("der_pressure")) // assume we had time derivative of solution before
      chef::combineField(m,"time derivative of solution","der_pressure","der_velocity","der_temperature");
    if (m->findField("mesh_vel_sim"))
      convertField(m, "mesh_vel_sim", "mesh_vel");
    if (m->findField("ctcn_elm_sim"))
      convertVtxFieldToElm(m, "ctcn_elm_sim", "err_tri_f");
    // destroy mesh size field
    if(m->findField("sizes"))  apf::destroyField(m->findField("sizes"));
    if(m->findField("frames")) apf::destroyField(m->findField("frames"));
  }

  void attachCurrentSizeField(apf::Mesh2*& m) {
    int  nsd = m->getDimension();
    if(m->findField("cur_size")) apf::destroyField(m->findField("cur_size"));
    apf::Field* cur_size = apf::createField(m, "cur_size", apf::SCALAR, apf::getConstant(nsd));
    // loop over non-BL elements
    apf::MeshEntity* e;
    apf::MeshIterator* eit = m->begin(nsd);
    while ((e = m->iterate(eit))) {
      pRegion meshRegion = reinterpret_cast<pRegion>(e);
      if (EN_isBLEntity(meshRegion)) continue;
      // set mesh size field
      double h = 0.0;
      if (m->getType(e) == apf::Mesh::TET)
        h = apf::computeShortestHeightInTet(m,e) * sqrt(3.0);
      else
        h = pc::getShortestEdgeLength(m,e);
      apf::setScalar(cur_size, e, 0, h);
    }
    m->end(eit);

    // get sim model
    apf::MeshSIM* sim_m = dynamic_cast<apf::MeshSIM*>(m);
    pParMesh sim_pm = sim_m->getMesh();
    pMesh pm = PM_mesh(sim_pm,0);

    gmi_model* gmiModel = sim_m->getModel();
    pGModel model = gmi_export_sim(gmiModel);

    // loop over model faces
    pGFace modelFace;
    pFace meshFace;
    pFace blFace;
    pRegion blRegion;
    FIter fIter;
    pEntity seed;
    pPList growthRegions = PList_new();
    pPList growthFaces = PList_new();
    GFIter gfIter = GM_faceIter(model);
    while((modelFace=GFIter_next(gfIter))){
      // loop over mesh faces on model face
      fIter = M_classifiedFaceIter(pm, modelFace, 1);
      while((meshFace = FIter_next(fIter))){
        apf::MeshEntity* apf_f = reinterpret_cast<apf::MeshEntity*>(meshFace);
        // check if BL base
        if (BL_isBaseEntity(meshFace, modelFace)) {
          // loop over BL regions and layers
          for(int faceSide = 0; faceSide < 2; faceSide++){
            int hasSeed = BL_stackSeedEntity(meshFace, modelFace, faceSide, NULL, &seed);
            if (hasSeed) {
              BL_growthRegionsAndLayerFaces((pRegion)seed, growthRegions, growthFaces, Layer_Entity);
              if (PList_size(growthRegions) >= (PList_size(growthFaces)-1)*3) { // tet
                for(int i = 0; i < PList_size(growthFaces); i++) {
                  blFace = (pFace)PList_item(growthFaces,i);
                  apf::MeshEntity* apf_f = reinterpret_cast<apf::MeshEntity*>(blFace);
                  double h = apf::computeShortestHeightInTri(m,apf_f) * sqrt(2.0);
                  for(int j = 0; j < 3; j++) {
                    if (i*3+j == PList_size(growthRegions)) break;
                    blRegion = (pRegion)PList_item(growthRegions,i*3+j);
                    // set mesh size field
                    apf::MeshEntity* apf_r = reinterpret_cast<apf::MeshEntity*>(blRegion);
                    apf::setScalar(cur_size, apf_r, 0, h);
                  }
                }
              }
              else if (PList_size(growthRegions) >= (PList_size(growthFaces)-1)) { // wedge
                for(int i = 0; i < PList_size(growthFaces); i++) {
                  if (i == PList_size(growthRegions)) break;
                  blFace = (pFace)PList_item(growthFaces,i);
                  apf::MeshEntity* apf_f = reinterpret_cast<apf::MeshEntity*>(blFace);
                  double h = apf::computeShortestHeightInTri(m,apf_f) * sqrt(2.0);
                  blRegion = (pRegion)PList_item(growthRegions,i);
                  // set mesh size field
                  apf::MeshEntity* apf_r = reinterpret_cast<apf::MeshEntity*>(blRegion);
                  apf::setScalar(cur_size, apf_r, 0, h);
                }
              }
            }
            else if (hasSeed < 0) {
              printf("not support blending BL mesh or miss some info!\n");
              exit(0);
            }
          }
        }
      }
      FIter_delete(fIter);
    }
    GFIter_delete(gfIter);
  }

  void initializeCtCn(apf::Mesh2*& m) {
    apf::Field* ctcn = apf::createSIMFieldOn(m, "ctcn_elm", apf::SCALAR);
    apf::MeshEntity* v;
    apf::MeshIterator* vit = m->begin(0);
    while ((v = m->iterate(vit))) {
      apf::setScalar(ctcn,v,0,1.0);
    }
    m->end(vit);
  }

  void applyMaxSizeBound(apf::Mesh2*& m, apf::Field* sizes, ph::Input& in) {
    //int barrelTag = 69;  // 2mm case
    //int domainTag = 113; // 2mm case
    //int projctTag = 108; // 2mm case
    //int domainTag = 1; //hollow barrel case
    //int a_barrelTag = 285; //back of barrel
    //int b_barrelTag = 377; //front of barrel


    // int domainTag = 171; //diamond (unfitted)
    //  int domainTag = 195; //capsule (unfitted)
    // int domainTag = 152; //capsule (m2_fit)
    // int domainTag = 153; //cpasule (m_fit_long)
    int domainTag = 92; //bullet stationary cube
    // int domainTag = 254;
    //int a_barrelTag = 185;
    //int b_barrelTag = 1;
    //int domainTag = 128;//2d dw case
    
    //apf::ModelEntity* bme = m->findModelEntity(3, barrelTag);
    
    apf::ModelEntity* dme = m->findModelEntity(3, domainTag); //for barrel problems
    
    //apf::ModelEntity* pme = m->findModelEntity(3, projctTag);

    //apf::ModelEntity* a_bme = m->findModelEntity(3,a_barrelTag); // for hb problems
    //apf::ModelEntity* b_bme = m->findModelEntity(3,b_barrelTag);//

    apf::Vector3 v_mag = apf::Vector3(0.0,0.0,0.0);
    //double dmeSizeUpperBound = 0.1; //upper bound for farfield

    apf::MeshEntity* v;
    apf::MeshIterator* vit = m->begin(0);
    while ((v = m->iterate(vit))) {
      apf::getVector(sizes,v,0,v_mag);
      for (int i = 0; i < 3; i++) {
        apf::ModelEntity* me = m->toModel(v);
        //if ((m->isInClosureOf(me, b_bme) || m->isInClosureOf(me, a_bme) ) && v_mag[i] > 0.016) v_mag[i] = 0.016;
        // if (m->isInClosureOf(me, dme) && v_mag[i] > in.simSizeUpperBound) v_mag[i] = in.simSizeUpperBound;
        if (v_mag[i] > in.simSizeUpperBound) v_mag[i] = in.simSizeUpperBound;
        //if (m->isInClosureOf(me, pme) && v_mag[i] > 0.004) v_mag[i] = 0.004;
        /*if(vertexIsInCylinder(v) && v_mag[i] > in.simSizeUpperBound)
          v_mag[i] = in.simSizeUpperBound; //set in adapt.inp*/
      }
      apf::setVector(sizes,v,0,v_mag);
    }
    m->end(vit);
  }

  void syncMeshSize(apf::Mesh2*& m, apf::Field* sizes) {
    PCU_Comm_Begin();
    apf::Copies remotes;
    apf::Vector3 v_mag = apf::Vector3(0.0,0.0,0.0);
    apf::MeshEntity* v;
    apf::MeshIterator* vit = m->begin(0);
    while ((v = m->iterate(vit))) {
      apf::getVector(sizes,v,0,v_mag);
      if(m->isShared(v)) {
        m->getRemotes(v, remotes);
        APF_ITERATE(apf::Copies, remotes, rit) {
          PCU_COMM_PACK(rit->first, rit->second);
          PCU_Comm_Pack(rit->first, &(v_mag[0]), sizeof(double));
        }
      }
    }
    m->end(vit);

    PCU_Comm_Send();
    while (PCU_Comm_Receive()) {
      apf::MeshEntity* rv;
      PCU_COMM_UNPACK(rv);
      double rv_mag;
      PCU_Comm_Unpack(&(rv_mag), sizeof(double));
      apf::getVector(sizes,v,0,v_mag);
      if(rv_mag < v_mag[0]) { // smaller wins
        v_mag = apf::Vector3(rv_mag,rv_mag,rv_mag);
        apf::setVector(sizes,rv,0,v_mag);
      }
    }
  }

  void syncMeshSize(apf::Mesh2*& m){
    apf::Field* sizes = m->findField("sizes");
    assert(sizes);

    if(!PCU_Comm_Self())
      std::cerr << "Begin sync mesh size" << std::endl;
    apf::synchronize(sizes);
  }

  void setupSimImprover(pVolumeMeshImprover vmi, pPList sim_fld_lst) {
    VolumeMeshImprover_setModifyBL(vmi, 1);
    VolumeMeshImprover_setShapeMetric(vmi, ShapeMetricType_VolLenRatio, 0.3);
    VolumeMeshImprover_setSmoothType(vmi, 1); // 0:Laplacian-based; 1:Gradient-based

    /* set fields to be mapped */
    if (PList_size(sim_fld_lst))
      VolumeMeshImprover_setMapFields(vmi, sim_fld_lst);
  }

  void setupSimAdapter(pMSAdapt adapter, ph::Input& in, apf::Mesh2*& m, pPList& sim_fld_lst) {
    MSA_setAdaptBL(adapter, 1);
    MSA_setExposedBLBehavior(adapter,BL_DisallowExposed);
    MSA_setBLSnapping(adapter, 0); // currently needed for parametric model
    MSA_setBLMinLayerAspectRatio(adapter, 0.0); // needed in parallel

    // double grad_rate = 1.0;
    // MSA_setSizeGradation(adapter, 0, grad_rate);
	

    /* attach mesh size field */
    phSolver::Input inp("solver.inp", "input.config");
    attachMeshSizeField(m, in, inp);
    apf::Field* sizes = m->findField("sizes");
    assert(sizes);
	
    /* initial ctcn field */
    pc::initializeCtCn(m);

    /* apply upper bound */
    pc::applyMaxSizeBound(m, sizes, in);

    /* attach current size field for verification */
    pc::attachCurrentSizeField(m);

    /* apply max number of element */
    // double cn = pc::applyMaxNumberElement(m, sizes, in);

    /* scale mesh if reach time resource bound */
    // pc::applyMaxTimeResource(m, sizes, in, inp);

    /* apply upper bound */
    // pc::applyMaxSizeBound(m, sizes, in);

    /* add mesh smooth/gradation function here */
    // pc::addSmoother(m, in.gradingFactor);

    pc::syncMeshSize(m);

    apf::Field* position = m->findField("motion_coords");
    apf::Vector3 pos;
    apf::Field* P_Filt = m->findField("P_Filt");
    apf::Field* PG_avg = m->findField("PG_avg");
    apf::Vector3 PG;
    apf::Field* shock_param = m->findField("Shock Param");
    apf::Field* surf_id = m->findField("Shock_ID");
    apf::Field* shock_vert = m->findField("shock_vert");
    apf::Field* plan_field = m->findField("plan_vert");
    apf::Field* shock_line = m->findField("shock_line");
    apf::Field* shock_line_marker = m->findField("shock_line_marker");

    apf::Field* aniso_size = apf::createFieldOn(m,"aniso_size",apf::MATRIX);
    apf::Field* shk_id = apf::createField(m,"shk_id",apf::SCALAR,apf::getConstant(3));

    apf::Matrix3x3 an_size;

    apf::Vector3 v_mag = apf::Vector3(0.0,0.0,0.0);
    apf::MeshEntity* v;
    apf::MeshIterator* vit = m->begin(0);
    while ((v = m->iterate(vit))) {
      //set isotropic/anisotropic size field at vertices

      double shock = apf::getScalar(shock_vert,v,0);
      double planarity = apf::getScalar(plan_field,v,0);
      double scale_factor = in.sizeScaleFactor;
 
      apf::getVector(sizes,v,0,v_mag);
      pVertex meshVertex = reinterpret_cast<pVertex>(v);
      bool aniso = in.anisotropicShockAdaptation;
      double aspect_ratio = in.anisotropicShockAR;

      apf::Vector3 pt;
      m->getPoint(v,0,pt);
      double iso_marker = 1;
      if(in.extendShocks){
        iso_marker = apf::getScalar(shock_line_marker,v,0);
      }

      if(aniso && shock > 1.0) {
        printf("Doing mesh intersection.\n");
        // Intersecting shock elements. This can replace the next else-if
        // conditional. Previously they were distinguishing by iso_marker, we
        // will distinguish by unique adjacent shock id count. (3. can we replace it?)
        // Get set of unique adjacent shock Ids.
        // For each ID, get an associated shock element.
        // - Only pick one element per shock ID.
        // Loop over shock ID element set.
        // - Get PG for that shock element. <- 2a. can we get this from
        //                                         elements? Yes. Field P_Filt
        // - Later, conditionally calculate average PG from regions for each shock ID.
        // - Also need shock_line (t1). <- 2b. can we get this from elements?
        // - Math. Compute t2 like below on elements.
        //   - Get metric intersection of anisoSizes.
        //     - Common diagonalization. <-- Don't worry about this yet.

        // Get adjacent shock elements.
        apf::Adjacent adj;
        m->getAdjacent(v, 3, adj);
        std::vector<apf::MeshEntity*> shocks;
        for (int i = 0; i < adj.size(); ++i) {
          if (apf::getScalar(shock_param, adj[i], 0)) {
            shocks.push_back(adj[i]);
          }
        }

        // Filter to get one element per shock id.
        std::set<double> shock_ids;
        std::vector<apf::Vector3> p_filts;

        for (int i = 0; i < shocks.size(); ++i) {
          auto result = shock_ids.insert(std::abs(apf::getScalar(surf_id, shocks[i], 0)));
          if (result.second) {
            apf::Vector3 pressure;
            apf::getVector(P_Filt, shocks[i], 0, pressure);
            p_filts.push_back(pressure);
          }
        }

        // Only works for 2-shock case.
        // For more shocks, if in the same plane as two current shocks, use
        // the circle. If on a new plane, use a sphere.
        assert(p_filts.size() == 2);

        // Shocks should not be parallel. If they were, they would have been
        // identified as the same shock.
        assert(abs(p_filts[0] * p_filts[1]) < 0.9);

        apf::Vector3 t1 = apf::cross(p_filts[0], p_filts[1]).normalize();
        apf::Vector3 inplane1 = p_filts[0].normalize();
        apf::Vector3 inplane2 = apf::cross(inplane1, t1);

        inplane1 = inplane1 * v_mag[0] / aspect_ratio * scale_factor;
        inplane2 = inplane2 * v_mag[0] / aspect_ratio * scale_factor;
        t1 = t1 * v_mag[0] * scale_factor;

        double anisoSize[3][3];
        inplane1.toArray(anisoSize[0]);
        inplane2.toArray(anisoSize[1]);
        t1.toArray(anisoSize[2]);

        MSA_setAnisoVertexSize(adapter, meshVertex, anisoSize);
      }
      else if(aniso && shock && iso_marker < 0.4){
        //this section is unverified
        //vertices which get 2 short component, 1 long (intersecting shocks)
        apf::getComponents(PG_avg,v,0,&PG[0]);
        apf::Vector3 n = PG.normalize();
        apf::Vector3 t1;
        apf::getVector(shock_line,v,0,t1);
        apf::Vector3 t2 = apf::cross(n,t1).normalize();

        n = n * v_mag[0] / aspect_ratio * scale_factor;
        t1 = t1 * v_mag[0] * scale_factor;
        t2 = t2 * v_mag[0] / aspect_ratio * scale_factor;

        double anisoSize[3][3] = {{n[0],n[1],n[2]},{t1[0],t1[1],t1[2]},{t2[0],t2[1],t2[2]}}; 

        MSA_setAnisoVertexSize(adapter, meshVertex, anisoSize);
      }
      else if(aniso && shock){
        // vertices which get 1 short component, 2 long
        // rewrite to use apf::Vector calculations as opposed to by hand
        apf::getVector(PG_avg,v,0,PG);
        apf::Vector3 n1 = PG.normalize();

        apf::Vector3 e = {1,0,0};
        if(n1.y() < n1.x()){
          if(n1.z() < n1.y()){
            e.x() = 0;
            e.z() = 1;
          }
          else{
            e.x() = 0;
            e.y() = 1;
          } 
        }
        else if (n1.z() < n1.x()){
          e.x() = 0;
          e.z() = 1;
        }

        apf::Vector3 t1 = apf::cross(n1, e).normalize();
        apf::Vector3 t2 = apf::cross(n1, t1).normalize();

        n1 = n1 * v_mag[0] / aspect_ratio * scale_factor;
        t1 = t1 * v_mag[0] * scale_factor;
        t2 = t2 * v_mag[0] * scale_factor;

        apf::Vector3 Vecs[3] = {n1, t1, t2};
        an_size = apf::Matrix3x3(Vecs);

        apf::setMatrix(aniso_size,v,0,an_size);

        double anisoSize[3][3];
        an_size.toArray(anisoSize);

        double n_val = v_mag[0] / aspect_ratio * scale_factor;
        double t_val = v_mag[0] * scale_factor;

        v_mag[0] = n_val; v_mag[1] = t_val; v_mag[2] = t_val;
        MSA_setAnisoVertexSize(adapter, meshVertex, anisoSize);
      }
      else {
        // isotropic refinement based on VMS error
        an_size = apf::Matrix3x3(v_mag[0],0,0,0,v_mag[0],0,0,0,v_mag[0]);
        apf::setMatrix(aniso_size,v,0,an_size);
        MSA_setVertexSize(adapter, meshVertex, v_mag[0]*scale_factor);
      }
    } //end iterate over vertices
    m->end(vit);

    /* use current size field */
    if(!PCU_Comm_Self())
      printf("Start mesh adapt of setting size field\n");

    /* write error and mesh size */
    pc::writeSequence(m, in.timeStepNumber, "error_mesh_size_");

    /* set fields to be mapped */
    PList_clear(sim_fld_lst);
    if (in.solutionMigration) {
      sim_fld_lst = getSimFieldList(in, m);
      MSA_setMapFields(adapter, sim_fld_lst);
    }
  }

  void runMeshAdapter(ph::Input& in, apf::Mesh2*& m, apf::Field*& orgSF, int step) {
    /* use the size field of the mesh before mesh motion */
    apf::Field* szFld = orgSF;

    if(in.simmetrixMesh == 1) {
      if (in.writeSimLog)
        Sim_logOn("sim_mesh_adaptation.log");
      pProgress progress = Progress_new();
      Progress_setDefaultCallback(progress);

      apf::MeshSIM* sim_m = dynamic_cast<apf::MeshSIM*>(m);
      pParMesh sim_pm = sim_m->getMesh();
      pMesh pm = PM_mesh(sim_pm,0);

      gmi_model* gmiModel = sim_m->getModel();
      pGModel model = gmi_export_sim(gmiModel);

      // declaration
      VIter vIter;
      pVertex meshVertex;

      pProgress progress_tmp = Progress_new();
      Progress_setDefaultCallback(progress_tmp);
      apf::MeshSIM* apf_msim = dynamic_cast<apf::MeshSIM*>(m);
      pParMesh ppm = apf_msim->getMesh();
  
      
      /*Begin shock detection code: Written by Isaac Tam, 2020*/ 
      
      PM_write(ppm, "ShockInd.sms",progress_tmp);
 
      apf::Field* S_Ind = m->findField("Shock_Ind");
      apf::Field* P_Filt = m->findField("P_Filt");
      apf::Field* vms_err = m->findField("VMS_error");\

      apf::Field* sol = m->findField("solution");
      apf::Field* td_sol = m->findField("time derivative of solution");
     

      //mesh resampling hacks
      // bool mesh_resample = false; // resample mesh from solution transfer in paraview
      int mesh_resample = in.manualMeshResample;
      if (mesh_resample && step == 1){ // just the initial step
        std::cout << "Enter manual solution migration" << std::endl;

        // read in the resampled data
        std::ifstream in_stream("outpoints.csv");
        if (!in_stream.good())
        {
          std::cerr << "can't find the resampled data" << std::endl;
          exit(1);
        }
        std::string debugstring1;
        in_stream >> debugstring1; // remove header
        std::string Lin_Parse;
        std::string delimin = ",";

        // std::map<int, std::vector<double>> resample_map; // vert_ID, solution vec
        std::vector<std::vector<double>> resample_map;

        std::cout << "Begin loop file" << std::endl;
        while (in_stream >> Lin_Parse)
        {
          std::vector<std::string> working_vec = ParseDeliminatedString(Lin_Parse, delimin);
          // populate map

          if (PCU_Comm_Self() == std::stoi(working_vec[5]))
          {
            std::vector<double> dropin;
            dropin.push_back(std::stod(working_vec[0]));
            dropin.push_back(std::stod(working_vec[1]));
            dropin.push_back(std::stod(working_vec[2]));
            dropin.push_back(std::stod(working_vec[3]));
            dropin.push_back(std::stod(working_vec[4]));

            // resample_map[std::stoi(working_vec[6])] = dropin;

            resample_map.push_back(dropin);
            // if (!PCU_Comm_Self()) std::cout << working_vec[5] <<" " <<working_vec[6]<< " " << dropin[0] <<std::endl;
          }
        }
        std::cout << "Begin iterate mesh to change solution" << std::endl;

        // solution can be changed here
        apf::MeshEntity *vert_tmp;
        apf::MeshIterator *v_iter = m->begin(0);
        while ((vert_tmp = m->iterate(v_iter)))
        {
          // grab solution
          apf::NewArray<double> solu_tmp(apf::countComponents(sol));
          apf::getComponents(sol, vert_tmp, 0, &solu_tmp[0]);
          pEntity ent1 = reinterpret_cast<pEntity>(vert_tmp);
          int v_ID = EN_id(ent1);

          std::vector<double> insert_vec = resample_map[v_ID];

          // std::cout <<  insert_vec[0] << std::endl;
          // std::cout <<  insert_vec[1] << std::endl;
          // std::cout <<  insert_vec[2] << std::endl;
          // std::cout <<  insert_vec[3] << std::endl;
          // std::cout <<  insert_vec[4] << std::endl;

          solu_tmp[0] = insert_vec[0];
          solu_tmp[1] = insert_vec[1];
          solu_tmp[2] = insert_vec[2];
          solu_tmp[3] = insert_vec[3];
          solu_tmp[4] = insert_vec[4];

          apf::setComponents(sol, vert_tmp, 0, &solu_tmp[0]);
        }
        pc::writeSequence(m, in.timeStepNumber, "resampled_prior_");
      }// end mesh resampling hacks

      /* Get smooth pressure gradient field, chef level shock detector*/
     
      apf::Field* PG_avg= apf::createFieldOn(m,"PG_avg",apf::VECTOR);
      apf::Field* shk_det= apf::createFieldOn(m,"shk_det",apf::SCALAR);
      apf::Field* num_elms= apf::createFieldOn(m,"num_elms",apf::SCALAR);
 
      apf::NewArray<double> holder(apf::countComponents(P_Filt));
      apf::Vector3 PG_add = apf::Vector3(0.0,0.0,0.0);

      //parallel implementation for pressure gradient field - Steven Spreizer, 2023
      PCU_Comm_Begin();
 
      apf::MeshEntity* v_tmp;
      apf::MeshIterator* v_itr = m->begin(0);
      while ((v_tmp = m->iterate(v_itr))) {
        //get adjacent nodes
        PG_add[0]=0.0; PG_add[1]=0.0; PG_add[2]=0.0;
        apf::Adjacent Adja;
        m->getAdjacent(v_tmp,3,Adja);
        int num_elm = Adja.getSize(); //number of elements for current sum
          
        //loop for contributions from neighboring elements
        for (size_t i=0; i<num_elm; i++){
          apf::getComponents(P_Filt, Adja[i], 0, &holder[0]);
          PG_add[0]=PG_add[0]+holder[0];
          PG_add[1]=PG_add[1]+holder[1];
          PG_add[2]=PG_add[2]+holder[2];
        }

        //set value
        apf::setVector(PG_avg, v_tmp, 0, PG_add); 
        apf::setScalar(num_elms,v_tmp,0, Adja.getSize());

        //pack up data for communications
        if(!m->isOwned(v_tmp)){
          apf::Copies remotes;
          m->getRemotes(v_tmp,remotes);
          int owningPart = m->getOwner(v_tmp);

          PCU_COMM_PACK(owningPart,remotes[owningPart]); //send entity
          PCU_COMM_PACK(owningPart,num_elm); //send int 
          PCU_COMM_PACK(owningPart,PG_add); //send apf::Vector3
        }
      } //end PG accumulation and comm packaging iteration
      m->end(v_itr);

      PCU_Comm_Send();

      apf::MeshEntity* ent;
      int received_num_elm;
      apf::Vector3 received_PG;
      apf::Vector3 current_PG;
      int current_num_elm;

      while(PCU_Comm_Receive()){
        //receive comms 
        PCU_COMM_UNPACK(ent);
        PCU_COMM_UNPACK(received_num_elm);
        PCU_COMM_UNPACK(received_PG);

        if(!m->isOwned(ent)){
          std::cout << "ERROR: Data sent to non-owner entity" << std::endl;
          std::exit(1);
        }

        apf::getComponents(PG_avg,ent,0,&current_PG[0]);
        current_num_elm = apf::getScalar(num_elms,ent,0);
        int total_elms = current_num_elm + received_num_elm;
        //compute accumulation from remotes
        current_PG[0] = current_PG[0] + received_PG[0];
        current_PG[1] = current_PG[1] + received_PG[1];
        current_PG[2] = current_PG[2] + received_PG[2];

        //update current fields (on owner)
        apf::setComponents(PG_avg, ent, 0, &current_PG[0]); 
        apf::setScalar(num_elms,ent,0, total_elms);
      } //end receive comms

      //perform averaging
      int total_elms;
      v_itr = m->begin(0);
      while((v_tmp = m->iterate(v_itr))){
        if(m->isOwned(v_tmp)){
          total_elms = apf::getScalar(num_elms,v_tmp,0);
          apf::getComponents(PG_avg,v_tmp,0,&current_PG[0]);
          current_PG[0] /= (double)total_elms;
          current_PG[1] /= (double)total_elms;
          current_PG[2] /= (double)total_elms;
          apf::setComponents(PG_avg, v_tmp, 0, &current_PG[0]); 
        }
      } //end averaging iteration
      m->end(v_itr);

      //synchronize data on all processes
      apf::synchronize(PG_avg);
      apf::destroyField(num_elms);    

      //iterate over vertices and calculate normal mach number (Lovely & Haimes), store in shk_det
      apf::Vector3 PG_tmp = apf::Vector3(0.0,0.0,0.0);
      v_itr = m->begin(0);
      while((v_tmp = m->iterate(v_itr))){
        double loc_det=0.0;
        apf::NewArray<double> sol_tmp(apf::countComponents(sol));
        apf::NewArray<double> td_sol_tmp(apf::countComponents(td_sol));

        apf::getComponents(sol, v_tmp, 0, &sol_tmp[0]);
        apf::getComponents(td_sol, v_tmp, 0, &td_sol_tmp[0]);
        apf::getComponents(PG_avg, v_tmp, 0, &PG_tmp[0]);

        loc_det= sol_tmp[1]*PG_tmp[0] + sol_tmp[2]*PG_tmp[1] + sol_tmp[3]*PG_tmp[2];//term 2
        loc_det= loc_det + td_sol_tmp[0];//term 1
        loc_det= loc_det/ ( sqrt(1.4*287*sol_tmp[4])  );//speed of sound
        loc_det= loc_det/ ( sqrt(PG_tmp[0]*PG_tmp[0] + PG_tmp[1]*PG_tmp[1] + PG_tmp[2]*PG_tmp[2]) );
        // P Grad Mag
        
        apf::setScalar(shk_det, v_tmp, 0, loc_det);
      } //end iterate setting normal mach number
      m->end(v_itr);

      /* Loop through elements and output IDs that contain a shock after filtering */      
      int nsd = 3;
      apf::MeshEntity* elm;

      apf::NewArray<double> Shock_Ind(apf::countComponents(S_Ind));
      apf::NewArray<double> P_filter(apf::countComponents(P_Filt));
      apf::NewArray<double> VMS_err(apf::countComponents(vms_err));

      //default values for filtering
      //pressure gradient threshold
      double P_thres_max   = 10000000000000000000000.0; // accept the highest values
      double P_thres_min   = 0.0; // accept lowest values
      //VMS error threshold
      double VMS_thres_max = 10000000000000000000000.0;
      double VMS_thres_min = 0.0;
      
      //read from "Shock.inp"
      std::ifstream in_str("Shock.inp");
      if (!in_str.good()){
        if (!PCU_Comm_Self()){
          std::cout << "Can't open Shock.inp to read. Using defaults.\n" 
                    << std::endl;
        }
      }
      else{
        std::string parse;
        while (in_str >> parse){
          if (parse == "P_thres_max"){
            in_str >> P_thres_max;
          }
          else if (parse == "P_thres_min"){
            in_str >> P_thres_min;
          }
          else if (parse == "VMS_thres_max"){
            in_str >> VMS_thres_max;
          }
          else if (parse == "VMS_thres_min"){
            in_str >> VMS_thres_min;
          }
        }
      }

      //Setup for shock parameter recording
      std::vector<int> ShkIDs;
      std::vector<apf::Vector3 > ShkLocs;
      std::vector<double> ShkFilt;
      
      //1 or 0 indicator of shock in element
      apf::Field* Shock_Param = apf::createField(m, "Shock_Param", apf::SCALAR, apf::getConstant(nsd));
      //shock system ID 
      apf::Field* Shock_IDs   = apf::createField(m, "Shock_ID", apf::SCALAR, apf::getConstant(nsd));
      //planarity indicator
      apf::Field* plan_val = apf::createField(m, "planarity", apf::SCALAR, apf::getConstant(nsd));
      
      //iterate over elements to set shock parameter
      int num_shock_elms = 0;
      apf::MeshIterator* it = m->begin(nsd);
      while ((elm = m->iterate(it))) {
        apf::getComponents(S_Ind, elm, 0, &Shock_Ind[0]);//phasta det calc
        apf::getComponents(P_Filt, elm, 0, &P_filter[0]);
        apf::getComponents(vms_err, elm, 0, &VMS_err[0]);

        apf::setScalar(Shock_Param,elm,0,0.0);
        apf::setScalar(Shock_IDs,elm,0,0.0);
        apf::setScalar(plan_val,elm,0,-1.0);

        //get momentum vms err
        double moment_err = 0.0;
        moment_err = sqrt(VMS_err[1]*VMS_err[1]
                          +VMS_err[2]*VMS_err[2]
                          +VMS_err[3]*VMS_err[3]);      
        
        //actual filtering
        if (P_thres_max > P_filter[3] && P_filter[3] > P_thres_min){//pressure filter, between max and min
          if ( VMS_thres_max > moment_err && moment_err > VMS_thres_min) {//similar VMS_filter
            apf::Adjacent Adja;
            m->getAdjacent(elm,0,Adja);
            //find iso surface elms
            bool loc_gt = false; bool loc_lt=false; double loc_check=0.0;
            //loop adjacent vertices to see if element contains a normal 
            //mach number of 1 => element contains shock
            for (size_t i=0; i<Adja.getSize(); i++){
              apf::getComponents(shk_det, Adja[i], 0, &loc_check);
              if (loc_check >= 1.0){
                loc_gt=true;
              }
              if (loc_check <= 1.0){
                loc_lt=true;
              }
            }

            if (loc_gt && loc_lt){ // iso surface elements from chef    
              /* option for geometric filtering as well:
              if desired add if condtiional to only look at elements in filter          
              apf::Vector3 x3p = apf::Vector3(0.0,0.0,0.0);
              m->getPoint(Adja[0],0,x3p);//for geom filter
              double x3p_r= sqrt(x3p[1]*x3p[1] + x3p[2]*x3p[2]);
              !(x3p[0] < 2.01 && x3p_r < 0.07)){ // geom filtering
              */
              pEntity ent = reinterpret_cast<pEntity>(elm);
              int rID = EN_id(ent);
              
              apf::Vector3 CellCent = apf::getLinearCentroid(m, elm);
              ShkLocs.push_back(CellCent);

              ShkIDs.push_back(rID);
              ShkFilt.push_back(P_filter[3]);
              apf::setScalar(Shock_Param, elm, 0, 1.0);
              if (num_shock_elms <= 10) printf("Setting Shock_Param\n");
              // apf::setScalar(Shock_IDs, elm, 0, rID);
              apf::setScalar(Shock_IDs, elm, 0, 0);
              ++num_shock_elms;
            }
          }
        }
      } //end iterate over elements to set shock parameter
      m->end(it);


      //option to output shock information to file
      bool shock_to_file = true;
      if(shock_to_file){
        std::clock_t shock_to_file_c = std::clock();
        std::string ShkFile ="ShockElms-"+ std::to_string(PCU_Comm_Self()) + ".txt";
        std::ofstream ShockOut(ShkFile);
        ShockOut << std::setprecision(std::numeric_limits<double>::digits10 + 1);
        if(!ShockOut.good()){
            std::cerr << "Can't open "<< ShkFile <<" to write.\n" <<std::endl;
            exit(1);
        }
        double tmpO1; double tmpO2; double tmpO3;
        for (unsigned int i=0; i<ShkIDs.size();i++){
            ShockOut <<ShkIDs[i] << ",";
            tmpO1=ShkLocs[i][0]; tmpO2=ShkLocs[i][1]; tmpO3=ShkLocs[i][2];
            ShockOut << tmpO1 << "," << tmpO2 <<"," <<tmpO3 << "," << ShkFilt[i] <<std::endl;
        }
        std::clock_t shock_to_file_c_end = std::clock();
        std::cout << "Writing ShockElms took " << 1000*(shock_to_file_c_end - shock_to_file_c) / CLOCKS_PER_SEC << "ms" << std::endl;
      }

      PCU_Add_Int(num_shock_elms);
      if(!PCU_Comm_Self()){
        std::cout<< "Shock Elms count: " << num_shock_elms <<std::endl;     
      }

      // if(true){
      if(in.shockIDSegment){
        std::cerr << "Entering shock system ID and segmentation code " << std::endl;
        std::clock_t ssid_c = std::clock();
        /* Shock system ID and segmentation using Armadillo */
        if(!PCU_Comm_Self()){ //proceed in serial 
          std::vector<sElm> all_shock_elements;
          double spacing = 0.1; //closeness parameter
          
          for(int i = 0; i < PCU_Comm_Peers();i++){
            //loop over each file and add to working vector
            ReadData(all_shock_elements,i,all_shock_elements.size());
          }
          GetNbs(all_shock_elements,spacing);

          std::vector<int> SurfSizes;
          std::cout << "Beginning SurfaceSorter" << std::endl;
          SurfSizes = SurfaceSorter(all_shock_elements);
          std::cout << "SurfaceSorter complete" << std::endl;
 
          int threshold = 10; //threshold num elements to consider shock system as just noise
          for (int i = 1; i < SurfSizes.size(); i++)
          { // check surface sizes
            if(SurfSizes[i] < threshold){
              for(int j = 0; j < all_shock_elements.size(); j++){
                if(all_shock_elements[j].SurfID == i){
                  all_shock_elements[j].SurfID = 0;
                }
              }
            }
            else{
              std::cout << i << " " << SurfSizes[i] << std::endl;
            }
          }
        
          //segmentation procedure
          CalcPlanarity(all_shock_elements);
          WriteData(all_shock_elements);

          std::clock_t ssid_c_endw = std::clock();
          std::cout << "Writing to files took " << 1000 * (ssid_c_endw - ssid_c) / CLOCKS_PER_SEC << "ms." << std::endl;
        }

        //add code to now read back in file, compare element centroid position with positions on current process
        //if element is from that process, read the shock surface ID and write to the field data

        PCU_Barrier();

        std::string file = "processed_shock_elements.txt";
        in_str = std::ifstream(file);
        if(!in_str){
          std::cerr << "Could not open " << file << " to read." << std::endl;
        }

        std::string parse;
        string delimiter = ",";
        int id_iter = 0;

        if(!PCU_Comm_Self()){
          std::cout << "Looping shock elements to add ID to field data" << std::endl;
        }

        std::clock_t ssid_read_c = std::clock();

        while (in_str >> parse)
        {
          std::vector<std::string> work_vec;

          std::string s = parse;
          size_t pos = 0;
          std::string token;
          while ((pos = s.find(delimiter)) != std::string::npos)
          {
            token = s.substr(0, pos);
            work_vec.push_back(token);
            s.erase(0, pos + delimiter.length());
          }
          work_vec.push_back(s); // work vec has the x,y,z coords

          sElm drop;
          double x_pos = stod(work_vec[1]);
          double y_pos = stod(work_vec[2]);
          double z_pos = stod(work_vec[3]);
          int surf_id = stoi(work_vec[4]);
          double plan_ind = stod(work_vec[5]);

          //iterate elements to see if there's a match
          double eps = 1e-7; //tolerance for checking if same point
          apf::MeshEntity* e;
          apf::Vector3 pt; 
          apf::MeshIterator* itr = m->begin(3);
          while((e = m->iterate(itr))){
            pt = apf::getLinearCentroid(m,e);
            bool match = abs(pt[0] - x_pos) < eps && abs(pt[1] - y_pos) < eps 
                      && abs(pt[2] - z_pos) < eps;
            if(match){
              apf::setScalar(Shock_IDs,e,0,surf_id);
              apf::setScalar(plan_val,e,0,plan_ind);

              //undo shock detection of elements that are just noise
              if(surf_id == 0){
                // apf::setScalar(Shock_Param,e,0,0);
              }
              break;
            }
          }
        }

        std::clock_t ssid_read_c_end = std::clock();

        cout << "["<< PCU_Comm_Self() << "]\tDone reading and adding to field data" << endl;
        std::cout << "That took " << 1000 * (ssid_read_c_end - ssid_read_c) / CLOCKS_PER_SEC << "ms." << std::endl;

        std::cerr << "Exiting shock system ID and segmentation code " << std::endl;
        std::cout << "SSID/Seg took " << 1000 * (ssid_read_c_end - ssid_c) / CLOCKS_PER_SEC << "ms total." << std::endl;
      }
      PCU_Barrier();

      if(in.extendShocks){
        //routine to extend the shock detection result
        extendShocks(m,in);

        /* 
        recalculate planarity after extension
        */
        
        // Make this conditional upon newly found shock elements.
        double spacing = 0.1;
        calcPlanarity(m,spacing,in);

        /*
        create second PG field for special handling of 2 short-1 long dimension adaptation
        where interactions are found (high planarity)
        */
        interactionHandling(m,in);
      }


      /* Create vertex level shock indicator for aniso adaptation */
      // Mark vertices that have an adjacent shock containing element 
      PCU_Comm_Begin();
      apf::Field* shock_vert = apf::createFieldOn(m,"shock_vert",apf::SCALAR);
      v_itr = m->begin(0);
      while((v_tmp = m->iterate(v_itr))){
        double v_shock = 0;
        apf::Adjacent elements;
        m->getAdjacent(v_tmp,3,elements);
        std::set<double> adjacentIds;
        for(int i = 0; i < elements.size(); ++i){
          if(apf::getScalar(Shock_Param,elements[i],0)){
            double id = std::abs(apf::getScalar(Shock_IDs, elements[i], 0));
            adjacentIds.insert(id);
          }
        }
        apf::setScalar(shock_vert,v_tmp,0,double(adjacentIds.size()));

        if(!m->isOwned(v_tmp)){
          apf::Copies remotes;
          m->getRemotes(v_tmp,remotes);
          int owningPart = m->getOwner(v_tmp);

          PCU_COMM_PACK(owningPart,remotes[owningPart]); //send entity
          PCU_COMM_PACK(owningPart,v_shock); //send int 
        }
      }
      m->end(v_itr);

      PCU_Comm_Send();

      double v_shock_rec;
      while(PCU_Comm_Receive()){
        PCU_COMM_UNPACK(ent);
        PCU_COMM_UNPACK(v_shock_rec);

        if(!m->isOwned(ent)){
          std::cerr << "ERROR: Data sent to non-owner entity" << std::endl;
          std::exit(1);
        }

        double v_shock_curr = apf::getScalar(shock_vert,ent,0);
        v_shock_curr = v_shock_curr + v_shock_rec;

        apf::setScalar(shock_vert,ent,0,v_shock_curr);
      }

      apf::synchronize(shock_vert);

      /* Vertex level planarity indicator */

      PCU_Comm_Begin();
      /* not implemented correctly need to fix
      need to average only neighbor shock elements not all elements */
      apf::Field* plan_vert = apf::createFieldOn(m,"plan_vert",apf::SCALAR);
      v_itr = m->begin(0);
      while((v_tmp = m->iterate(v_itr))){
        double v_plan = 0;
        apf::Adjacent elements;
        m->getAdjacent(v_tmp,3,elements);
        int count = 0;
        for(int i = 0; i < elements.size(); ++i){
          if(apf::getScalar(plan_val,elements[i],0) >= 0){
            count++;
            v_plan += apf::getScalar(plan_val,elements[i],0);
          }
        }
        v_plan /= count;
        apf::setScalar(plan_vert,v_tmp,0,v_plan);

        if(!m->isOwned(v_tmp)){
          apf::Copies remotes;
          m->getRemotes(v_tmp,remotes);
          int owningPart = m->getOwner(v_tmp);

          PCU_COMM_PACK(owningPart,remotes[owningPart]); //send entity
          PCU_COMM_PACK(owningPart,v_plan); //send int 
        }
      }
      m->end(v_itr);

      PCU_Comm_Send();

      double v_plan_rec;
      while(PCU_Comm_Receive()){
        PCU_COMM_UNPACK(ent);
        PCU_COMM_UNPACK(v_plan_rec);

        if(!m->isOwned(ent)){
          std::cerr << "ERROR: Data sent to non-owner entity" << std::endl;
          std::exit(1);
        }

        double v_plan_curr = apf::getScalar(plan_vert,ent,0);
        v_plan_curr = std::min(v_plan_curr,v_plan_rec);

        apf::setScalar(plan_vert,ent,0,v_plan_curr);
      }

      apf::synchronize(plan_vert);
      /* End vertex level shock indicator */
      
      // create the Simmetrix adapter *
      if(!PCU_Comm_Self())
        printf("Start mesh adapt\n");
      pMSAdapt adapter = MSA_new(sim_pm, 1);
      pPList sim_fld_lst = PList_new();
      setupSimAdapter(adapter, in, m, sim_fld_lst);

      // run the adapter *
      if(!PCU_Comm_Self())
        printf("do real mesh adapt\n");
      if(step > 1)
        MSA_adapt(adapter, progress);
      // MSA_adapt(adapter, progress); //temporary to test solution migration
      MSA_delete(adapter);


      if(step > 1) {
      // create Simmetrix improver *
      pVolumeMeshImprover vmi = VolumeMeshImprover_new(sim_pm);
      setupSimImprover(vmi, sim_fld_lst);

      // run the improver *
      VolumeMeshImprover_execute(vmi, progress);
      VolumeMeshImprover_delete(vmi);

      PList_clear(sim_fld_lst);
      PList_delete(sim_fld_lst);

      // load balance *
      pc::balanceEqualWeights(sim_pm, progress);
      }


      // write mesh *
      if(!PCU_Comm_Self())
        printf("write mesh after mesh adaptation\n");
      writeSIMMesh(sim_pm, in.timeStepNumber, "sim_mesh_");
      
      
      Progress_delete(progress);

      // transfer data back to apf *
      if (in.solutionMigration)
        transferSimFields(m); 
    }
    else {
      assert(szFld);
      apf::synchronize(szFld);
      apf::synchronize(m->getCoordinateField());
      /* do SCOREC mesh adaptation */
      chef::adapt(m,szFld,in);
      chef::balance(in,m);
    }
    m->verify();
  }

  /**
   * @brief Generate averaged vertex pressure gradient field from region field.
   * @param m APF mesh.
   * @input-field "P_Filt" VECTOR
   * @working-field "pc_spf_num_rgns"
   * @output-field "PG_avg" VECTOR
   * @return A pointer to the "PG_avg" field.
   */
  apf::Field* smoothP_Filt(apf::Mesh2* m) {
    // Get input field.
    apf::Field* P_Filt = m->findField("P_Filt");

    // Create working field.
    apf::Field* num_rgns = apf::createFieldOn(m, "pc_spf_num_rgns", apf::SCALAR);

    // Create output field.
    apf::Field* PG_avg = apf::createFieldOn(m, "PG_avg", apf::VECTOR);

    // P_Filt is currently a 8-component field, but this may change.
    apf::NewArray<double> p_filt_val(apf::countComponents(P_Filt));

    // Iterate through every vertex.
    apf::MeshIterator* it = m->begin(0);
    for (apf::MeshEntity* v = m->iterate(it); v; v = m->iterate(it)) {
      apf::Vector3 pg_sum;
      apf::Adjacent adja;
      m->getAdjacent(v, 3, adja);

      // Sum P_Filt for owned adjacent regions.
      int adjacent_owned = 0;
      for (int i = 0; i < adja.size(); ++i) {
        if (m->isOwned(adja[i])) {
          apf::getComponents(P_Filt, adja[i], 0, p_filt_val.begin());

          // Although P_Filt has 8 components, we only need the first 3.
          pg_sum += apf::Vector3(p_filt_val[0], p_filt_val[1], p_filt_val[2]);
          ++adjacent_owned;
        }
      }

      apf::setVector(PG_avg, v, 0, pg_sum);
      apf::setScalar(num_rgns, v, 0, adjacent_owned);
    }
    m->end(it);

    // Accumulate sums along partition boundary.
    apf::accumulate(PG_avg);
    apf::accumulate(num_rgns);

    // Calculate averages.
    it = m->begin(0);
    for (apf::MeshEntity* v = m->iterate(it); v; v = m->iterate(it)) {
      apf::Vector3 pg_sum;
      apf::getVector(PG_avg, v, 0, pg_sum);
      double count = apf::getScalar(num_rgns, v, 0);

      apf::setVector(PG_avg, v, 0, pg_sum / count);
    }
    m->end(it);

    // No need to re-synchronize because all part boundaries have the same field
    // data.

    // Clean up intermediate fields.
    apf::destroyField(num_rgns);

    return PG_avg;
  }

  /**
   * @brief Calculate normal mach number for mesh m.
   * 
   * @param m APF Mesh.
   * @input-field "PG_avg"
   * @input-field "solution"
   * @input-field "time derivative of solution"
   * @output-field "shk_det" (SCALAR)
   * @return A pointer to the shk_det Field.
   */
  apf::Field* calcNormalMach(apf::Mesh* m) {
    apf::Field* PG_avg = m->findField("PG_avg");
    apf::Field* sol = m->findField("solution");
    apf::Field* td_sol = m->findField("time derivative of solution");
    apf::Field* shk_det = apf::createFieldOn(m, "shk_det", apf::SCALAR);

    apf::Vector3 pg_avg;
    apf::NewArray<double> sol_tmp(apf::countComponents(sol));
    apf::NewArray<double> td_sol_tmp(apf::countComponents(td_sol));

    apf::MeshIterator* it = m->begin(0);
    for (apf::MeshEntity* v = m->iterate(it); v; v = m->iterate(it)) {
      apf::getVector(PG_avg, v, 0, pg_avg);
      apf::getComponents(sol, v, 0, sol_tmp.begin());
      apf::getComponents(td_sol, v, 0, td_sol_tmp.begin());

      apf::Vector3 sol_v(sol_tmp.begin() + 1);

      double a = sqrt(1.4*287 * sol_tmp[4]); // = speed of sound = sqrt(gamma*R)
      double mach_n = (td_sol_tmp[0]/a + sol_v * pg_avg)/pg_avg.getLength();

      apf::setScalar(shk_det, v, 0, mach_n);
    }
    m->end(it);

    return shk_det;
  }

  /**
   * @brief ShockFilter is a class that handles reading of Shock.inp files and
   * filtering data based on those inputs.
   */
	class ShockFilter {
	public:
		static ShockFilter ReadFile(const char* ShockInp);

    /**
     * @brief Check if a pressure value is within input bounds.
     * 
     * @param P A pressure value to test.
     * @return true if P is in bounds.
     */
    bool checkPressure(double P) const noexcept {
      return P_min <= P && P <= P_max;
    }

    /**
     * @brief Check if a VMS error value is within input bounds.
     * 
     * @param VMS A VMS error to test.
     * @return true if VMS is in bounds.
     */
    bool checkVMS(double VMS) const noexcept {
      return VMS_min <= VMS && VMS <= VMS_max;
    }

	private:
		ShockFilter(double pmin, double pmax, double vmin, double vmax):
				P_min(pmin), P_max(pmax), VMS_min(vmin), VMS_max(vmax) {}

		ShockFilter(): ShockFilter(0, HUGE_VAL, 0, HUGE_VAL) {}

		// Pressure thresholds.
		double P_min, P_max;

		// VMS error thresholds.
		double VMS_min, VMS_max;
	};

	/**
   * @brief Read Shock.inp file.
   * 
   * @param ShockInp The filename to parse. Usually "Shock.inp".
   * @return ShockFilter A newly constructed ShockFilter.
   */
  ShockFilter ShockFilter::ReadFile(const char* ShockInp) {
    ShockFilter sf;

		std::ifstream in_str("Shock.inp");
		if (in_str.good()) {
			std::string word;
			while (in_str >> word) {
				if (word == "P_thres_max") {
					in_str >> sf.P_max;
				} else if (word == "P_thres_min") {
					in_str >> sf.P_min;
				} else if (word == "VMS_thres_max") {
					in_str >> sf.VMS_max;
				} else if (word == "VMS_thres_min") {
					in_str >> sf.VMS_min;
				}
			}
		} else {
			if (!PCU_Comm_Self()) {
				std::cerr << "ERROR: failed to open " << ShockInp
									<< "Falling back to defaults." << std::endl;
			}
		}

		return sf;
	}

  /**
   * @brief Mark elements as shock based on normal mach = 1 in adjacent edge.
   * 
   * 1-mach line is detected in edges by IVT then propogated to all adjacent
   * elements.
   * 
   * @param m An apf::Mesh.
   * @param shk_det The input field with normal mach numbers.
   * @param Shock_Param The output field for shock indicators.
   */
  void locateShockEdges(apf::Mesh* m, apf::Field* shk_det, apf::Field* Shock_Param) {
    // Iterate over edges.
    apf::MeshIterator* it = m->begin(1);
    for (apf::MeshEntity* e = m->iterate(it); e; e = m->iterate(it)) {
      apf::Downward down;
      m->getDownward(e, 0, down);
      double M0 = apf::getScalar(shk_det, down[0], 0);
      double M1 = apf::getScalar(shk_det, down[1], 0);

      double target_mach_n = 1;
      // Intermediate Value Theorem. Shock line of 1 on an edge indicates shock
      // in the surrounding elements (regions).
      /* FIXME: Some potential for numerical error (subtracting close values)
       * but the sign should remain correct. */
      if ((M0 - target_mach_n) * (M1 - target_mach_n) < 0) {
        apf::Adjacent adja;
        m->getAdjacent(e, 3, adja);
        for (size_t i = 0; i < adja.size(); ++i) {
          apf::setScalar(Shock_Param, adja[i], 0, 1);
        }
      }
    }
    m->end(it);

    // Propogate Shock_Param along mesh partition boundary.
    apf::Sharing* sharing = apf::getSharing(m);
    apf::sharedReduction(Shock_Param, sharing, false,
                         apf::ReductionMax<double>());
  }

  /**
   * @brief Filter Shock_Param indicator via pressure residual and VMS error
   * magnitude.
   * 
   * @param m An apf::Mesh.
   * @param Shock_Param The Shock_Param field to filter.
   * @param filt ShockFilter inputs.
   * @input-field "P_Filt"
   * @input-field "VMS_error"
   */
  void filterShockParam(apf::Mesh* m, apf::Field* Shock_Param, const ShockFilter& filt) {
    apf::Field* P_Filt = m->findField("P_Filt");
    apf::Field* VMS_error = m->findField("VMS_error");

    // P_Filt and VMS_error are 8*real in the Fortran code.
    apf::NewArray<double> p(apf::countComponents(P_Filt));
    apf::NewArray<double> vms(apf::countComponents(VMS_error));

    apf::MeshIterator* it = m->begin(3);
    for (apf::MeshEntity* e = m->iterate(it); e; m->iterate(it)) {
      apf::Vector3 vms_v(vms.begin() + 1);
      if (!(filt.checkPressure(p[3]) && filt.checkVMS(vms_v.getLength()))) {
        apf::setScalar(Shock_Param, e, 0, 0);
      }
    }
    m->end(it);

    // We only operate on individual elements. Given field data is synchronized,
    // partition boundaries stay synchronized.
  }

  /**
   * @brief Detect shocks on a mesh using normal mach number with pressure gradient filtering.
   * 
   * Lovely and Haimes.
   * 
   * @param in PHASTA input config.
   * @param m An apf::Mesh.
   * @input-field "P_Filt"
   * @working-field "shk_det"
   * @output-field "PG_avg"
   * @output-field "Shock_Param"
  // FIXME: * @output-field "Shock_ID"
  // FIXME: * @output-field "plan_val"
   */
  void detectShocksSerial(const ph::Input& in, apf::Mesh2* m) {
    apf::Field* PG_avg = smoothP_Filt(m);

    // Calculate normal mach number.
    apf::Field* shk_det = calcNormalMach(m);

    // Shock_Param is a scalar indicator (1 or 0) field on 3D elements.
    apf::Field* Shock_Param = apf::createField(m, "Shock_Param", apf::SCALAR,
                                               apf::getConstant(3));

    // Find shock line using IVT for edges.
    locateShockEdges(m, shk_det, Shock_Param);

    // Filter Shock_Param using pressure/vms error.
    filterShockParam(m, Shock_Param, ShockFilter::ReadFile("Shock.inp"));

    // Destroy intermediate fields.
    apf::destroyField(shk_det);
  }

  using BFScheck = std::function<bool(apf::Mesh* m, apf::MeshEntity* c, apf::MeshEntity* e)>;
  using BFSaction = std::function<void(apf::Mesh* m, apf::MeshEntity* e, int component, int distance)>;

  void serialBFS(apf::Mesh* m, int dimension, BFScheck check, BFSaction action) {
    int bridgeDimension = 0;
    if (dimension == 0) bridgeDimension = 1;

    std::set<apf::MeshEntity*> visited;

    int component = 0;
    apf::MeshIterator *it = m->begin(dimension);
    for (apf::MeshEntity* e = m->iterate(it); e; e = m->iterate(it)) {
      int distance = 0;
      std::queue<apf::MeshEntity*> Q, Q_next;

      Q.push(e);
      visited.insert(e);

      // BFS loop.
      for (; !Q.empty(); Q.pop()) {
        if (check(m, e, Q.front())) {
          action(m, e, component, distance);

          apf::Adjacent neighbors;
          apf::getBridgeAdjacent(m, Q.front(), bridgeDimension, dimension, neighbors);
          for (size_t i = 0; i < neighbors.size(); ++i) {
            if (visited.find(neighbors[i]) == visited.end()) {
              Q_next.push(neighbors[i]);
            }
          }

          if (Q.empty()) {
            Q = Q_next;
            ++distance;
          }
        }
      }

      ++component;
    }
  }

  /**
   * @brief Connect fragmented shock elements. Breadth-first search stopping at
   * after 3 neighbors.
   * 
   * @param m APF Mesh.
   * @input-field "Shock_Param"
   * @output-field "Shock_Param"
   */
  void defragmentShocksSerial(apf::Mesh2* m) {
    // for each element:
    // - if no adjacent neighbors are shock, do BFS.
    // - if majority of
    // spacing = 0.1
  }

  /**
   * @brief Connect shock components.
   * 
   * @param in PHASTA input.
   * @param m APF mesh.
   * @input-field "Shock_Param" Shock indicator field. 1 or 0.
   * @working-field "pc_cs_bfs"
   * @output-field "Shock_Param"
   */
  void labelShocksSerial(const ph::Input& in, apf::Mesh2* m) {
    // Do BFS

    // Get input field.
    apf::Field* Shock_Param = m->findField("Shock_Param");

    // Idea: each value is a (PCU_Rank, Component ID, distance) triplet.
    apf::Field* pc_cs_bfs = apf::createField(m, "pc_cs_bfs", apf::VECTOR,
                                             apf::getConstant(3));
    
    // Initialize BFS field.
    apf::MeshIterator* it = m->begin(3);
    for (apf::MeshEntity* e = m->iterate(it); e; e = m->iterate(it)) {
      apf::Vector3 bfs_trip(-1, -1, std::numeric_limits<double>::infinity());
      apf::setVector(pc_cs_bfs, e, 0, bfs_trip);
    }
    m->end(it);

    serialBFS(m, 3, [Shock_Param, pc_cs_bfs](apf::Mesh* m, apf::MeshEntity* c, apf::MeshEntity* e) {
      apf::Vector3 v;
      apf::getVector(pc_cs_bfs, e, 0, v);
      return apf::getScalar(Shock_Param, e, 0) == 1 && v.y() < 0;
    }, [pc_cs_bfs](apf::Mesh* m, apf::MeshEntity* e, int c, int d) {
      apf::Vector3 v;
      apf::getVector(pc_cs_bfs, e, 0, v);
      v.x() = PCU_Comm_Self();
      v.y() = c;
      v.z() = d;
      apf::setVector(pc_cs_bfs, e, 0, v);
    });

    // FIXME: Process part boundary (repeated triplet max.)
  }

  void denoiseShocksSerial(const ph::Input& in, apf::Mesh2* m);
  void segmentShocksSerial(const ph::Input& in, apf::Mesh2* m);
  void setupSizeField(const ph::Input& in, apf::Mesh2* m, apf::Field* szFld);

  /**
   * @brief Run mesh adaptation.
   * 
   * @param in PHASTA input config.
   * @param m APF mesh.
   * @param szFld Size field.
   * @param step Current execution step.
   * @input-field "P_Filt"
   * @working-field ""
   */
  void runMeshAdapterSerial(ph::Input& in, apf::Mesh2* m, apf::Field* szFld, int step) {
    PCU_ALWAYS_ASSERT(PCU_Comm_Peers() == 1);
    m->verify();

    // 0. Detect shock using normal mach number.
    detectShocksSerial(in, m);

    // 1. Gaps and fragments.
    // 1.a. Breadth-first search radius 3 or 5 neighbors.
    defragmentShocksSerial(m);

    // 2. Component connection (breadth first search).
    // 2.a. OK to connect distinct shock systems.
    labelShocksSerial(in, m);

    // 3. Filter components to remove noise (systems with <10 elements).
    denoiseShocksSerial(in, m);

    // 4. Divide shock systems with planarity.
    if (in.shockIDSegment) {
      segmentShocksSerial(in, m);
    }

    // 5. Extension and intersection.


    // 6. Size-field.
    setupSizeField(in, m, szFld);
    if (in.simmetrixMesh == 1) {
      setupSimAdapter(...);
      // run adapter, improver, etc.
    } else {
      chef::adapt(m, szFld, in);
      chef::balance(in, m);
    }

    m->verify();
  }
}
