#include <exception>
#include <fstream>

#include "csv.h"

namespace test {
  std::vector<std::string> tokenize(std::string line) {
    std::vector<std::string> tokens;
    for (size_t n = 0, end; n < line.size(); n = end + 1) {
      end = line.find(',', n);
      tokens.push_back(line.substr(n, end - n));
      if (end == std::string::npos) break;
    }
    // Account for trailing comma.
    if (line.size() && line.back() == ',') tokens.push_back("");
    return tokens;
  }

  class CSVException : public std::exception {
  public:
    CSVException(int line, std::string m) {
      msg = "line " + std::to_string(line) + ": " + m;
    }
    const char* what(void) const noexcept {
      return msg.c_str();
    }
  private:
    std::string msg;
  };

  std::vector<apf::Vector3> readCSVPoints(std::string csvFile) {
    std::vector<apf::Vector3> points;
    std::ifstream str(csvFile);
    std::string line;
    std::getline(str, line); // skip header
    for (int lineno = 2; std::getline(str, line); ++lineno) {
      // Remove newline from line.
      size_t newline = line.find_first_of("\r\n");
      if (newline != std::string::npos) line.erase(newline);
      std::vector<std::string> fields = tokenize(line);
      if (fields.size() != 3)
        throw CSVException(lineno, "wrong number of fields");
      try {
        apf::Vector3 pt(std::stod(fields[0]), std::stod(fields[1]),
          std::stod(fields[2]));
        points.push_back(pt);
      } catch (const std::exception& err) {
        throw CSVException(lineno, err.what());
      }
    }
    return points;
  }

  std::vector<FieldRow> readCSVField(std::string csvFile,
    std::string fieldName) {
    std::vector<FieldRow> rows;
    std::ifstream str(csvFile);
    std::string line;
    std::getline(str, line);
    std::vector<std::string> header = tokenize(line);
    if (header.size() != 4 || header[0] != fieldName ||
      header[1] != "x" || header[2] != "y" || header[3] != "z") {
      throw CSVException(0, "invalid header");
    }
    for (int lineno = 2; std::getline(str, line); ++lineno) {
      // Remove newline from line.
      size_t newline = line.find_first_of("\r\n");
      if (newline != std::string::npos) line.erase(newline);
      std::vector<std::string> fields = tokenize(line);
      if (fields.size() != 4)
        throw CSVException(lineno, "wrong number of fields");
      try {
        apf::Vector3 pt(std::stod(fields[1]), std::stod(fields[2]),
          std::stod(fields[3]));
        rows.emplace_back(std::stod(fields[0]), pt);
      } catch (const std::exception& err) {
        throw CSVException(lineno, err.what());
      }
    }
    return rows;
  }

  void writeCSVField(std::string csvFile, std::string fieldName,
    const std::vector<FieldRow>& rows) {
    std::ofstream file(csvFile);
    file << fieldName << ",x,y,z\n";
    for (size_t i = 0; i < rows.size(); ++i) {
      file << rows[i].val << ',' << rows[i].lc.x() << ','
        << rows[i].lc.y() << ',' << rows[i].lc.z() << '\n';
    }
  }
} // namespace test
