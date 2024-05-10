#include <iostream>
#include <fstream>
#include <queue>
#include <vector>
#include <unordered_set>

#include <apf.h>
#include <apfMesh2.h>
#include <apfShape.h>
#include <PCU.h>
#include <pcu_util.h>

namespace pc {
  class ShockParam {
  public:
    static constexpr double NONE = 0;
    static constexpr double SHOCK = 1;
    static constexpr double FRAGMENT = 2;

  private:
    // Static only class has private constructor.
    ShockParam();
  };

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
          apf::setScalar(Shock_Param, adja[i], 0, ShockParam::SHOCK);
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
        apf::setScalar(Shock_Param, e, 0, ShockParam::NONE);
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

    // Initalize Shock_Param to NONE (ensure it is a complete field).
    apf::MeshIterator *it = m->begin(3);
    for (apf::MeshEntity *e = m->iterate(it); e; e = m->iterate(it)) {
      apf::setScalar(Shock_Param, e, 0, ShockParam::NONE);
    }
    m->end(it);

    // Find shock line using IVT for edges.
    locateShockEdges(m, shk_det, Shock_Param);

    // Filter Shock_Param using pressure/vms error.
    filterShockParam(m, Shock_Param, ShockFilter::ReadFile("Shock.inp"));

    // Destroy intermediate fields.
    apf::destroyField(shk_det);
  }

	struct BFSarg {
		apf::Mesh* m;
		apf::MeshEntity* c; // First element for this BFS loop.
    apf::MeshEntity* e; // Current element.
		int component, distance;
	};

  // Return value for BFScheck.
  enum class BFSresult { ACT, CONTINUE, BREAK };

  using BFScheck = std::function<BFSresult(const BFSarg& arg)>;
  using BFSaction = std::function<void(const BFSarg& arg)>;

  /**
   * Do a breadth-first search on mesh `m` for all elements of dimension `dim`.
   *
   * @param m The APF mesh to search through.
   * @param dim The dimenion to search through.
   * @param check A callable argument that determines what to do for a given
   *              entity.
   * @param action The action to perform on entities where `check()` returns
   *               BFSresult::ACT.
   */
	void serialBFS(apf::Mesh* m, int dim, BFScheck check, BFSaction action) {
    // Bridge through points by default.
    int bridgeDim = 0;
    if (dim == 0) bridgeDim = 1;

    std::unordered_set<apf::MeshEntity*> visited(m->count(dim));
    BFSarg arg = (BFSarg){m, nullptr, nullptr, 0, 0};

    bool bfs_done = false;
    int component = 0;
    apf::MeshIterator *it = m->begin(dim);
    for (apf::MeshEntity* e = m->iterate(it); e && !bfs_done; e = m->iterate(it)) {
      int distance = 0;
      std::queue<apf::MeshEntity*> Q, Q_next;

      arg.c = e;

      // BFS loop.
      Q.push(e);
      while (!Q.empty()) {
        // Mark element as visited.
        visited.insert(Q.front());

        // Setup BFSarg.
        arg.e = Q.front();
        arg.component = component;
        arg.distance = distance;

        BFSresult br = check(arg);
        if (br == BFSresult::BREAK) {
          bfs_done = true;
          break;
        } else if (br == BFSresult::ACT) {
          action(arg);

          apf::Adjacent neighbors;
          apf::getBridgeAdjacent(m, Q.front(), bridgeDim, dim, neighbors);
          for (size_t i = 0; i < neighbors.size(); ++i) {
            if (visited.find(neighbors[i]) == visited.end()) {
              Q_next.push(neighbors[i]);
            }
          }
        }

        Q.pop();
        if (Q.empty()) {
          Q.swap(Q_next);
          ++distance;
        }
      }

      ++component;
    }
    m->end(it);
  }

  /**
   * Do a breadth-first search on mesh `m` starting from `e` through dimension
   * `dim`.
   *
   * @param m An APF mesh.
   * @param dim The dimension to search through.
   * @param e The starting element.
   * @param check A callable argument that determines what to do for a given
   *              entity.
   * @param action The action to perform on entities where `check()` returns
                   BFSresult::ACT.
   */
  void serialBFS(apf::Mesh* m, int dim, apf::MeshEntity* e, BFScheck check,
                 BFSaction action) {
    int bridgeDim = 0;
    if (dim == 0) bridgeDim = 1;

    std::set<apf::MeshEntity*> visited;
    BFSarg arg = (BFSarg){m, e, e, 0, 0};

    int component = 0, distance = 0;
		std::queue<apf::MeshEntity*> Q, Q_next;

		// BFS loop.
		for (Q.push(e); !Q.empty();) {
			// Mark element as visited.
			visited.insert(Q.front());

			// Setup BFSarg.
			arg.e = Q.front();
			arg.distance = distance;

			BFSresult br = check(arg);
			if (br == BFSresult::BREAK) {
				break;
			} else if (br == BFSresult::ACT) {
				action(arg);

				apf::Adjacent neighbors;
				apf::getBridgeAdjacent(m, Q.front(), bridgeDim, dim, neighbors);
				for (size_t i = 0; i < neighbors.size(); ++i) {
					if (visited.find(neighbors[i]) == visited.end()) {
						Q_next.push(neighbors[i]);
					}
				}
      }

      Q.pop();
      if (Q.empty()) {
        Q.swap(Q_next);
        ++distance;
      }
    }
  }

  /**
   * Get the max radius of an entity -- that is, the distance between the \ref
   * apf::getLinearCentroid "linear centroid" and the farthest point from the
   * centroid.
   *
   * @param m An APF mesh on which `e` resides.
   * @param e The entity to measure.
   */
  double getMaxRadius(apf::Mesh* m, apf::MeshEntity* e) {
    apf::Vector3 lc = apf::getLinearCentroid(m, e);

    double radius = 0;

    // Get points.
    apf::Downward points;
    int points_size = m->getDownward(e, 0, points);
    for (int i = 0; i < points_size; ++i) {
      apf::Vector3 p;
      m->getPoint(points[i], 0, p);
      double r = (p - lc).getLength();
      radius = r > radius ? r : radius;
    }

    return radius;
  }

  /**
   * @brief Fuzzily defragment shock elements. Mark elements neighboring shock
   * elements as shock.
   * 
   * @param m APF Mesh.
   * @input-field "Shock_Param"
   * @output-field "Shock_Param"
   */
  void defragmentShocksSerial(const ph::Input& in, apf::Mesh2* m) {
    apf::Field* Shock_Param = m->findField("Shock_Param");

    apf::MeshIterator* it = m->begin(3);
    for (apf::MeshEntity* e = m->iterate(it); e; e = m->iterate(it)) {
      if (apf::getScalar(Shock_Param, e, 0) == ShockParam::SHOCK) {
        double dist_max = getMaxRadius(m, e); // Max linear search distance.
        serialBFS(m, 3, e, [dist_max](const BFSarg& arg) -> BFSresult {
          if (arg.distance > 1) return BFSresult::BREAK;
          apf::Vector3 e_lc = apf::getLinearCentroid(arg.m, arg.e);
          apf::Vector3 c_lc = apf::getLinearCentroid(arg.m, arg.c);
          apf::Vector3 dist =  e_lc - c_lc;
          // Only process elements within spacing distance.
          return dist * dist < dist_max * dist_max ? BFSresult::ACT : BFSresult::CONTINUE;
        }, [Shock_Param](const BFSarg& arg){
          if (apf::getScalar(Shock_Param, arg.e, 0) == ShockParam::NONE) {
            apf::setScalar(Shock_Param, arg.e, 0, ShockParam::FRAGMENT);
          }
        });
      }
    }
  }

  using Shocks = std::vector<std::vector<apf::MeshEntity*>>;

  /**
   * @brief Connect shock components.
   * 
   * @param in PHASTA input.
   * @param m APF mesh.
   * @input-field "Shock_Param" Shock indicator field. 1 or 0.
   * @working-field "pc_cs_bfs"
   * @output-field "Shock_Param"
   * @return An array of labeled shocks.
   */
	Shocks labelShocksSerial(const ph::Input& in, apf::Mesh2* m) {
		// Get input field.
		apf::Field* Shock_Param = m->findField("Shock_Param");

		// Idea: each value is a (PCU_Rank, Component ID, distance) triplet.
		apf::Field* pc_cs_bfs =
				apf::createField(m, "pc_cs_bfs", apf::VECTOR, apf::getConstant(3));

		// Initialize BFS field.
		apf::MeshIterator* it = m->begin(3);
		for (apf::MeshEntity* e = m->iterate(it); e; e = m->iterate(it)) {
			apf::Vector3 bfs_trip(-1, -1, std::numeric_limits<double>::infinity());
			apf::setVector(pc_cs_bfs, e, 0, bfs_trip);
		}
		m->end(it);

		Shocks shocks;

		serialBFS(
				m, 3,
				[Shock_Param, pc_cs_bfs](const BFSarg& arg) -> BFSresult {
					apf::Vector3 v;
					apf::getVector(pc_cs_bfs, arg.e, 0, v);
					return apf::getScalar(Shock_Param, arg.e, 0) == ShockParam::SHOCK &&
								 v.y() < 0 ? BFSresult::ACT : BFSresult::CONTINUE;
				},
				[pc_cs_bfs, &shocks](const BFSarg& arg) {
					apf::Vector3 v;
					apf::getVector(pc_cs_bfs, arg.e, 0, v);
					v.x() = PCU_Comm_Self();
					v.y() = arg.component;
					v.z() = arg.distance;
					apf::setVector(pc_cs_bfs, arg.e, 0, v);

					if (arg.component >= shocks.size()) {
						shocks.push_back(std::vector<apf::MeshEntity*>());
						shocks[arg.component].push_back(arg.e);
					}
				});

		// FIXME: Process part boundary (repeated triplet max.)

		return shocks;
	}

	/**
   * @brief Remove components with < minsize elements.
   * 
   * @param in PHASTA input config.
   * @param m APF mesh.
   * @param shocks A collection of shocks which are each collections of entities.
   * @param minsize The minimum component size to retain.
   */
  void denoiseShocksSerial(const ph::Input& in, apf::Mesh2* m, Shocks& shocks, int minsize = 10) {
    apf::Field* Shock_Param = m->findField("Shock_Param");

    for (auto it = shocks.begin(); it != shocks.end();) {
      if (it->size() < minsize) {
        for (auto it2 = it->begin(); it2 != it->end(); ++it2) {
          apf::setScalar(Shock_Param, *it2, 0, ShockParam::NONE);
        }
        it = shocks.erase(it);
      } else {
        ++it;
      }
    }
  }

  /**
   * @brief Use planarity to segment shocks.
   * 
   * @param in PHASTA input config.
   * @param m APF mesh.
   */
  void segmentShocksSerial(const ph::Input& in, apf::Mesh2* m, Shocks& shocks) {}

  /**
   * @brief Setup Simmetrix/APF size field depending on input.
   * 
   * @param in PHASTA input config.
   * @param m APF mesh.
   * @param szFld Size field to setup.
   */
  void setupSizeField(const ph::Input& in, apf::Mesh2* m, apf::Field* szFld) {}

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
    defragmentShocksSerial(in, m);

    // 2. Component connection (breadth first search).
    // 2.a. OK to connect distinct shock systems.
    Shocks shocks = labelShocksSerial(in, m);

    // 3. Filter components to remove noise (systems with <10 elements).
    denoiseShocksSerial(in, m, shocks, 10);

    // 4. Divide shock systems with planarity.
    if (in.shockIDSegment) {
      segmentShocksSerial(in, m, shocks);
    }

    // 5. Extension and intersection.


    // 6. Size-field.
    setupSizeField(in, m, szFld);
    if (in.simmetrixMesh == 1) {
      // setupSimAdapter(...);
      // run adapter, improver, etc.
    } else {
      chef::adapt(m, szFld, in);
      chef::balance(in, m);
    }

    m->verify();
  }
} // namespace pc
