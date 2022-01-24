#include "initial_conditions_from_file.hpp"

#include "time_integrable_rhs.hpp"
#include "file_io/vtk_solution_reader.hpp"


void
InitialConditionsFromFile::set(const TimeIntegrableRHS& rhs, SolutionState& state) const
{
  reader_.setSolution(rhs.interiorIndices(), state);
}
