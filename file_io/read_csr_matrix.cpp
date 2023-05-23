#include "read_csr_matrix.hpp"

#include <fstream>



std::unique_ptr<CSRMatrix>
readCSRMatrix(const std::string& row_file, const std::string& col_file, const std::string& val_file)
{
  std::ifstream rows(row_file);
  std::string   string_input;

  const std::vector<int>     row_index = readArray<int>(row_file);
  const std::vector<int>     col_index = readArray<int>(col_file);
  const std::vector<real_wp> vals      = readArray<real_wp>(val_file);

//  std::cout << "reading csr matrix\n";
//  std::cout << "row_index size: " << row_index.size() << "\n";
//  std::cout << "col_index size: " << col_index.size() << "\n";
//  std::cout << "vals size: " << vals.size() << "\n";

  const int nrows = row_index.size()-1;
  const int nnz   = vals.size();

  auto csr_mat = std::make_unique<CSRMatrix>(nrows, nrows, nnz);

  {
    auto col_id = write_access_host(csr_mat->colIdx());
    auto val = write_access_host(csr_mat->values());
    for (int i = 0; i < col_id.size(); ++i){
      col_id[i] = col_index[i];
      val[i] = vals[i];
    }
  }

  {
    auto row_id = write_access_host(csr_mat->rowIdx());
    for (int i = 0; i < row_id.size(); ++i){
      row_id[i] = row_index[i];
    }
  }

  return csr_mat;
}
