#ifndef PC_TEST_SHOCKFUNC_H
#define PC_TEST_SHOCKFUNC_H

#include <apfMesh.h>

namespace pc {
  namespace _test {

    /**
     * @brief A ShockFunc is a test fixture that returns true if the element is
     * detected as shock and false otherwise.
     * 
     */
    class ShockFunc {
    public:
      bool operator()(apf::Mesh* m, apf::MeshEntity* e) {
        return eval(m, e);
      }
    
    private:
      virtual bool eval(apf::Mesh* m, apf::MeshEntity* e) = 0;
    };

    /**
     * @brief A 3D Planar shock function |Ax+By+Cz+D|<=tol.
     * 
     */
    class PlanarShockFunc : public ShockFunc {
    public:
      PlanarShockFunc(double A, double B, double C, double D, double tol) : 
      a(A), b(B), c(C), d(D), tol_(tol) {}

    private:
      virtual bool eval(apf::Mesh* m, apf::MeshEntity* e) {
        apf::Vector3 lc = apf::getLinearCentroid(m, e);
        return std::abs(a * lc.x() + b * lc.y() + c * lc.z() + d) <= tol_;
      }
      double a, b, c, d, tol_;
    };
  }
}

#endif // PC_TEST_SHOCK_FUNC
