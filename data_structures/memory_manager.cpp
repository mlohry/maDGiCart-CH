#include "memory_manager.hpp"


template <typename BaseType>
MemoryManager<BaseType>::~MemoryManager()
{
    if (!host_data_.empty()) {
      std::cout << "Warning: " << this->objectName() << " has " << host_data_.size()
                << " arrays still allocated on destruction.\n";

#ifndef NDEBUG
      std::cout << "look at my stack trace of size " << malloc_stack_trace_.size() << "\n";
  for (const auto& trace : malloc_stack_trace_){
    std::cout << trace.second << std::endl;
  }
#endif
    }
    else {

    }


}

// explicit instantiations to make sure they're compiled by the device compiler
template class MemoryManager<float>;
template class MemoryManager<double>;
template class MemoryManager<int>;
