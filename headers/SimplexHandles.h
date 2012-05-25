#ifndef SIMPLEXHANDLES_H
#define SIMPLEXHANDLES_H

namespace SimplexMesh {


/** Handle for referring to a vertex */
class VertexHandle 
{

private:
  int m_idx;
  explicit VertexHandle(int idx) : m_idx(idx) {}
  int idx() const { return m_idx; }

public:

  explicit VertexHandle() : m_idx(-1) {}

  friend class VertexIterator;
  friend class VertexVertexIterator;
  friend class VertexEdgeIterator;
  friend class EdgeVertexIterator;
  friend class VertexFaceIterator;
  
  friend class VertexIterator;
  friend class VertexEdgeIterator;friend class EdgeVertexIterator;

  friend class SimplicialComplex;
  
  template<class T> friend class VertexProperty;

  bool operator== (const VertexHandle& rhs) const { return m_idx == rhs.m_idx; }
  bool operator!= (const VertexHandle& rhs) const { return m_idx != rhs.m_idx; }
  bool operator< (const VertexHandle& rhs) const { return m_idx < rhs.m_idx; }

  bool isValid() const { return m_idx >= 0; }
  static VertexHandle invalid() { return VertexHandle(-1); }
};

/** Handle for referring to an edge */
class EdgeHandle 
{

private:
  int m_idx;
  explicit EdgeHandle(int idx) : m_idx(idx) {} 
  int idx() const { return m_idx; }

public:
  
  explicit EdgeHandle() : m_idx(-1) {}

  friend class EdgeIterator;
  friend class VertexEdgeIterator;
  friend class EdgeVertexIterator;
  friend class FaceEdgeIterator;
  friend class EdgeFaceIterator;
  friend class FaceVertexIterator;

  friend class EdgeIterator;
   friend class VertexEdgeIterator; friend class EdgeVertexIterator;
   friend class EdgeFaceIterator; friend class FaceEdgeIterator;

  friend class SimplicialComplex;

  template<class T> friend class EdgeProperty;

  bool operator== (const EdgeHandle& rhs) const { return m_idx == rhs.m_idx; }
  bool operator!= (const EdgeHandle& rhs) const { return m_idx != rhs.m_idx; }
  bool operator< (const EdgeHandle& rhs) const { return m_idx < rhs.m_idx; }

  bool isValid() const { return m_idx >= 0; }
  
  static EdgeHandle invalid() { return EdgeHandle(-1); }
};

/** Handle for referring to a face */
class FaceHandle
{
private:
  int m_idx;
  explicit FaceHandle(int idx) : m_idx(idx) {}
  int idx() const { return m_idx; }

public:

  explicit FaceHandle() : m_idx(-1) {}

  friend class FaceIterator;
  friend class FaceEdgeIterator;
  friend class EdgeFaceIterator;
  friend class FaceVertexIterator;
  friend class VertexFaceIterator;
  friend class SimplicialComplex;

  friend class FaceIterator;
  friend class FaceEdgeIterator;
  friend class EdgeFaceIterator;
  friend class FaceTetIterator;
  friend class TetFaceIterator;

  template<class T> friend class FaceProperty;

  bool operator== (const FaceHandle& rhs) const { return m_idx == rhs.m_idx; }
  bool operator!= (const FaceHandle& rhs) const { return m_idx != rhs.m_idx; }
  bool operator< (const FaceHandle& rhs) const { return m_idx < rhs.m_idx; }

  bool isValid() const { return m_idx >= 0; }
  
  static FaceHandle invalid() { return FaceHandle(-1); }
};

/** Handle for referring to a tet */
class TetHandle
{

private:
  int m_idx;
  explicit TetHandle(int idx) : m_idx(idx) {}
  int idx() const { return m_idx; }

public:

  explicit TetHandle() : m_idx(-1) {}

  friend class TetIterator;
  friend class SimplicialComplex;
  
  friend class TetIterator;
  friend class FaceTetIterator; friend class TetFaceIterator;
  
  template<class T> friend class TetProperty;

  bool operator== (const TetHandle& rhs) const { return m_idx == rhs.m_idx; }
  bool operator!= (const TetHandle& rhs) const { return m_idx != rhs.m_idx; }
  bool operator< (const TetHandle& rhs) const { return m_idx < rhs.m_idx; }

  bool isValid() const { return m_idx >= 0; }

  static TetHandle invalid() { return TetHandle(-1); }
};



} // namespace SimplexMesh

#endif // SIMPLEXHANDLES_H
