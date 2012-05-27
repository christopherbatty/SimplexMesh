#ifndef SIMPLEXITERATORS_H
#define SIMPLEXITERATORS_H

#include "SimplicialComplex.h"
#include "SimplexHandles.h"

#include <set>

namespace SimplexMesh 
{

//////////////////////////////////////////////////////////////////////////
//Basic iterators

class VertexIterator {
public:
   VertexIterator(const SimplicialComplex& obj);
   void advance();
   bool done() const;
   VertexHandle current() const;

private:
   int m_idx;
   const SimplicialComplex& m_obj;

};

class EdgeIterator {
public:
   EdgeIterator(const SimplicialComplex& obj);
   void advance();
   bool done() const;
   EdgeHandle current() const;

private:
   int m_idx;
   const SimplicialComplex& m_obj;

};


class FaceIterator {
public:
   FaceIterator(const SimplicialComplex& obj);
   void advance();
   bool done() const;
   FaceHandle current() const;

private:
   int m_idx;
   const SimplicialComplex& m_obj;

};

class TetIterator {
public:
   TetIterator(const SimplicialComplex& obj);
   void advance();
   bool done() const;
   TetHandle current() const;

private:
   int m_idx;
   const SimplicialComplex& m_obj;

};



//////////////////////////////////////////////////////////////////////////
//Adjacency iterators

class VertexEdgeIterator {
public:
   VertexEdgeIterator(const SimplicialComplex& obj, const VertexHandle& vh);
   void advance();
   bool done() const;
   EdgeHandle current() const;

private:
   int m_idx;
   VertexHandle m_vh;
   const SimplicialComplex& m_obj;

};

class EdgeVertexIterator {
public:
   EdgeVertexIterator(const SimplicialComplex& obj, const EdgeHandle& eh, bool ordered);
   void advance();
   bool done() const;
   VertexHandle current() const;

private:
   bool m_ordered;
   int m_idx;
   EdgeHandle m_eh;
   const SimplicialComplex& m_obj;

};

class EdgeFaceIterator {
public:
   EdgeFaceIterator(const SimplicialComplex& obj, const EdgeHandle& eh);
   void advance();
   bool done() const;
   FaceHandle current() const;

private:
   int m_idx;
   EdgeHandle m_eh;
   const SimplicialComplex& m_obj;

};

class FaceEdgeIterator {
public:
   FaceEdgeIterator(const SimplicialComplex& obj, const FaceHandle& eh, bool ordered);
   void advance();
   bool done() const;
   EdgeHandle current() const;

private:
   bool m_ordered;
   int m_idx;
   EdgeHandle m_cur;
   FaceHandle m_fh;
   const SimplicialComplex& m_obj;

};


class FaceTetIterator {
public:
   FaceTetIterator(const SimplicialComplex& obj, const FaceHandle& eh);
   void advance();
   bool done() const;
   TetHandle current() const;

private:
   int m_idx;
   FaceHandle m_fh;
   const SimplicialComplex& m_obj;

};

class TetFaceIterator {
public:
   TetFaceIterator(const SimplicialComplex& obj, const TetHandle& eh);
   void advance();
   bool done() const;
   FaceHandle current() const;

private:
   int m_idx;
   TetHandle m_th;
   const SimplicialComplex& m_obj;

};

//////////////////////////////////////////////////////////////////////////
//Trickier adjacency iterators

class VertexFaceIterator {
public:
   VertexFaceIterator(const SimplicialComplex& obj, const VertexHandle& vh);
   void advance();
   bool done() const;
   FaceHandle current() const;

private:

   const SimplicialComplex& m_obj;
   std::set<FaceHandle> m_faces;
   std::set<FaceHandle>::iterator m_fiter;

};

class FaceVertexIterator {
public:
   FaceVertexIterator(const SimplicialComplex& obj, const FaceHandle& fh, bool ordered);
   void advance();
   bool done() const;
   VertexHandle current() const;

private:
   FaceHandle m_fh;
   FaceEdgeIterator m_feit;
   const SimplicialComplex& m_obj;

};

class VertexTetIterator {
public:
   VertexTetIterator(const SimplicialComplex& obj, const VertexHandle& vh);
   void advance();
   bool done() const;
   TetHandle current() const;

private:

   const SimplicialComplex& m_obj;
   std::set<TetHandle> m_tets;
   std::set<TetHandle>::iterator m_titer;

};

class TetVertexIterator {
public:
   TetVertexIterator(const SimplicialComplex& obj, const TetHandle& th);
   void advance();
   bool done() const;
   VertexHandle current() const;

private:
   const SimplicialComplex& m_obj;
   std::set<VertexHandle> m_verts;
   std::set<VertexHandle>::iterator m_viter;

};


class EdgeTetIterator {
public:
   EdgeTetIterator (const SimplicialComplex& obj, const EdgeHandle& eh);
   void advance();
   bool done() const;
   TetHandle current() const;

private:

   const SimplicialComplex& m_obj;
   std::set<TetHandle> m_tets;
   std::set<TetHandle>::iterator m_titer;

};

class TetEdgeIterator {
public:
   TetEdgeIterator(const SimplicialComplex& obj, const TetHandle& th);
   void advance();
   bool done() const;
   EdgeHandle current() const;

private:
   const SimplicialComplex& m_obj;
   std::set<EdgeHandle> m_edges;
   std::set<EdgeHandle>::iterator m_eiter;

};


}

#endif