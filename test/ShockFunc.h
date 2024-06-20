#ifndef PC_TEST_SHOCKFUNC_H
#define PC_TEST_SHOCKFUNC_H

#include <string>

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

			static ShockFunc* makeFromString(std::string s, double tol);
    
    private:
      virtual bool eval(apf::Mesh* m, apf::MeshEntity* e) = 0;
    };

    /**
     * @brief A 3D Planar shock function |Ax+By+Cz-D|<=tol.
     * 
     */
    class PlanarShockFunc : public ShockFunc {
    public:
      PlanarShockFunc(double A, double B, double C, double D, double tol) : 
      a(A), b(B), c(C), d(D), tol_(tol) {}

    private:
      virtual bool eval(apf::Mesh* m, apf::MeshEntity* e) {
        apf::Vector3 lc = apf::getLinearCentroid(m, e);
        return std::abs(a * lc.x() + b * lc.y() + c * lc.z() - d) <= tol_;
      }
      double a, b, c, d, tol_;
    };

		/**
		 * @brief A 3D Spherical shock function |(x-A)^2+(y-B)^2+(z-C)^2-R^2|<=tol.
     */
		class SphericalShockFunc : public ShockFunc {
		public:
			SphericalShockFunc(double A, double B, double C, double R, double tol) :
			a(A), b(B), c(C), r(R), tol_(tol) {}

		private:
			virtual bool eval(apf::Mesh* m, apf::MeshEntity* e) {
				apf::Vector3 lc = apf::getLinearCentroid(m, e);
				double i = lc.x() - a, j = lc.y() - b, k = lc.z() - c;
				return std::abs(i * i + j * j + k * k - r * r) <= tol_;
			}
			double a, b, c, r, tol_;
		};

		/**
		 * @brief A 3D cylindrical (x-aligned) shock function |(y-B)^2+(z-C)^2-R^2|<=tol.
		 */
		class CylindricalXShockFunc : public ShockFunc {
		public:
			CylindricalXShockFunc(double B, double C, double R, double tol) :
			b(B), c(C), r(R), tol_(tol) {}

		private:
			virtual bool eval(apf::Mesh* m, apf::MeshEntity* e) {
				apf::Vector3 lc = apf::getLinearCentroid(m, e);
				double j = lc.y() - b, k = lc.z() - c;
				return std::abs(j * j + k * k - r * r) <= tol_;
			}
			double b, c, r, tol_;
		};

		/**
		 * @brief A 3D cylindrical (z-aligned) shock function |(x-A)^2+(y-B)^2-R^2|<=tol.
		 */
		class CylindricalZShockFunc : public ShockFunc {
		public:
			CylindricalZShockFunc(double A, double B, double R, double tol) :
			a(A), b(B), r(R), tol_(tol) {}

		private:
			virtual bool eval(apf::Mesh* m, apf::MeshEntity* e) {
				apf::Vector3 lc = apf::getLinearCentroid(m, e);
				double i = lc.x() - a, j = lc.y() - b;
				return std::abs(i * i + j * j - r * r) <= tol_;
			}
			double a, b, r, tol_;
		};

		ShockFunc* ShockFunc::makeFromString(std::string s, double tol) {
			std::string::size_type col = s.find(':');
			if (col == std::string::npos) {
				throw std::invalid_argument("malformed ShockFunc descriptor");
			}

			std::string s_type = s.substr(0, col);

			++col;
			std::vector<double> args;
			while (col < s.size()) {
				std::string::size_type comma = s.find(',', col);
				args.push_back(std::stof(s.substr(col, comma - col)));
				if (comma == std::string::npos) break;
				col = comma + 1;
			}

			if (s_type == "plane") {
				if (args.size() != 4) {
					throw std::invalid_argument("invalid ShockFunc arguments");
				}
				return new PlanarShockFunc(args[0], args[1], args[2], args[3], tol);
			} else if (s_type == "sphere") {
				if (args.size() != 4) {
					throw std::invalid_argument("invalid ShockFunc arguments");
				}
				return new SphericalShockFunc(args[0], args[1], args[2], args[3], tol);
			} else if (s_type == "cyl_x") {
				if (args.size() != 3) {
					throw std::invalid_argument("invalid ShockFunc arguments");
				}
				return new CylindricalXShockFunc(args[0], args[1], args[2], tol);
			} else if (s_type == "cyl_z") {
				if (args.size() != 3) {
					throw std::invalid_argument("invalid ShockFunc arguments");
				}
				return new CylindricalZShockFunc(args[0], args[1], args[2], tol);
			} else {
				throw std::invalid_argument("invalid ShockFunc type");
			}
		}
  }
}

#endif // PC_TEST_SHOCK_FUNC
