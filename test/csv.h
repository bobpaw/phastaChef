#include <vector>
#include <string>
#include <apfVector.h>

namespace test {
  std::vector<apf::Vector3> readCSVPoints(std::string csvFile);

  struct FieldRow {
    double val;
    apf::Vector3 lc;
    FieldRow(double v, const apf::Vector3& c): val(v), lc(c) {}
  };
  std::vector<FieldRow> readCSVField(std::string csvFile,
    std::string fieldName);
  void writeCSVField(std::string csvFile, std::string fieldName,
    const std::vector<FieldRow>& rows);
}
