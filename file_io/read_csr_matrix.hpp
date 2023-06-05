#pragma once

#include <boost/lexical_cast.hpp>
#include <memory>
#include <string>
#include "data_structures/csr_matrix.hpp"

std::unique_ptr<CSRMatrix>
readCSRMatrix(const std::string& row_file, const std::string& col_file, const std::string& val_file);

template <typename T>
std::vector<T>
readArray(const std::string& file)
{
  std::vector<T> values;
  std::string    string_input;
  std::ifstream  filein(file);

  if (filein.is_open()) {
    filein >> string_input;
    const int nvals = boost::lexical_cast<int>(string_input);
    values.reserve(nvals);

    for (int iline = 0; iline < nvals; ++iline) {
      filein >> string_input;
      values.push_back(boost::lexical_cast<T>(string_input));
    }
  }

  return values;
}
