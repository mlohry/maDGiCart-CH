#pragma once

#include <iostream>
#include <string>
#include "typedefs.hpp"

struct MemoryHandle {
  int_t       index      = -1;
  int_t       generation = -1;
  std::string name       = "UNNAMED";
};

inline std::ostream& operator<<(std::ostream& os, const MemoryHandle& handle)
{
  os << "index: " << handle.index << " name: " << handle.name
     << " generation: " << handle.generation;
  return os;
}


/*
 *
 * See
 * https://blog.molecular-matters.com/2013/05/17/adventures-in-data-oriented-design-part-3b-internal-references/
 *
 * "The idea is quite simple: instead of only using an index,
 * a handle also stores the generation in which the index was
 * created. The generation is simply a monotonically increasing
 * counter that gets incremented each time a data item is
 *  deleted. The generation is stored both inside the handle,
 * and for each data item. Whenever we want to access data
 * using a handle, the index’ generation and the data item’s
 *  generation need to match."
 *
 *
 * Generation could be a global iteration id
 *
 *
 * namespace backend
{
  VertexBufferHandle CreateVertexBuffer(unsigned int vertexCount, unsigned int stride, const void*
initialData); VertexBufferHandle CreateDynamicVertexBuffer(unsigned int vertexCount, unsigned int
stride); void DestroyVertexBuffer(VertexBufferHandle handle);

  IndexBufferHandle CreateIndexBuffer(unsigned int indexCount, IndexBuffer::Format::Enum format,
const void* initialData); void DestroyIndexBuffer(IndexBufferHandle handle);
}
 */
