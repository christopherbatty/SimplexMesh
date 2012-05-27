#include "SimplexIterators.h"

namespace SimplexMesh {

//VertexIterator

VertexIterator::VertexIterator(const SimplicialComplex& obj): m_obj(obj), m_idx(-1) {
   advance();
}

void VertexIterator::advance() {
   do {
      m_idx++;
   } while(m_idx < (int)m_obj.numVertexSlots() && !m_obj.m_V[m_idx]);
}

bool VertexIterator::done() const {
   return m_idx >= (int)m_obj.numVertexSlots();
}

VertexHandle VertexIterator::current() const {
   return VertexHandle(m_idx);
}

//EdgeIterator

EdgeIterator::EdgeIterator(const SimplicialComplex& obj): m_obj(obj), m_idx(-1) {
   advance();
}

void EdgeIterator::advance() {
   do {
      m_idx++;
   } while(m_idx < (int)m_obj.numEdgeSlots() && m_obj.m_EV.getNumEntriesInRow(m_idx) == 0);
}

bool EdgeIterator::done() const {
   return m_idx >= (int)m_obj.numEdgeSlots();
}

EdgeHandle EdgeIterator::current() const {
   return EdgeHandle(m_idx);
}

//FaceIterator

FaceIterator::FaceIterator(const SimplicialComplex& obj): m_obj(obj), m_idx(-1) {
   advance();
}

void FaceIterator::advance() {
   do {
      m_idx++;
   } while(m_idx < (int)m_obj.numFaceSlots() && m_obj.m_FE.getNumEntriesInRow(m_idx) == 0);
}

bool FaceIterator::done() const {
   return m_idx >= (int)m_obj.numFaceSlots();
}

FaceHandle FaceIterator::current() const {
   return FaceHandle(m_idx);
}


//TetIterator

TetIterator::TetIterator(const SimplicialComplex& obj): m_obj(obj), m_idx(-1) {
   advance();
}

void TetIterator::advance() {
   do {
      m_idx++;
   } while(m_idx < (int)m_obj.numTetSlots() && m_obj.m_TF.getNumEntriesInRow(m_idx) == 0);
}

bool TetIterator::done() const {
   return m_idx >= (int)m_obj.numTetSlots();
}

TetHandle TetIterator::current() const {
   return TetHandle(m_idx);
}


//////////////////////////////////////////////////////////////////////////
//Adjacency iterators


//VertexEdgeIterator
VertexEdgeIterator::VertexEdgeIterator(const SimplicialComplex& obj, const VertexHandle& vh): m_obj(obj), m_idx(0), m_vh(vh) {
}

void VertexEdgeIterator::advance() {
   m_idx++;
}

bool VertexEdgeIterator::done() const {
   return m_idx >= (int)m_obj.vertexIncidentEdgeCount(m_vh);
}

EdgeHandle VertexEdgeIterator::current() const {
   return m_idx >= (int)m_obj.vertexIncidentEdgeCount(m_vh) ?
      EdgeHandle::invalid() : EdgeHandle(m_obj.m_VE.getColByIndex(m_vh.idx(), m_idx));
}

//EdgeVertexIterator
EdgeVertexIterator::EdgeVertexIterator(const SimplicialComplex& obj, const EdgeHandle& eh, bool ordered): m_obj(obj), m_idx(0), m_eh(eh), m_ordered(ordered) {
}

void EdgeVertexIterator::advance() {
   m_idx++;
}

bool EdgeVertexIterator::done() const {
   return m_idx >= 2;
}

VertexHandle EdgeVertexIterator::current() const {
   if(m_ordered) {
      if(m_idx == 0)
         return m_obj.fromVertex(m_eh);
      else if(m_idx == 1)
         return m_obj.toVertex(m_eh);
      else return VertexHandle::invalid();
   }
   else {
      return m_idx >= 2 ?
         VertexHandle::invalid() : VertexHandle(m_obj.m_EV.getColByIndex(m_eh.idx(), m_idx));
   }
}


//EdgeFaceIterator
EdgeFaceIterator::EdgeFaceIterator(const SimplicialComplex& obj, const EdgeHandle& eh): m_obj(obj), m_idx(0), m_eh(eh) {
}

void EdgeFaceIterator::advance() {
   m_idx++;
}

bool EdgeFaceIterator::done() const {
   return m_idx >= (int)m_obj.edgeIncidentFaceCount(m_eh);
}

FaceHandle EdgeFaceIterator::current() const {
   return m_idx >= (int)m_obj.edgeIncidentFaceCount(m_eh) ?
      FaceHandle::invalid() : FaceHandle(m_obj.m_EF.getColByIndex(m_eh.idx(), m_idx));
}

//FaceEdgeIterator
FaceEdgeIterator::FaceEdgeIterator(const SimplicialComplex& obj, const FaceHandle& fh, bool ordered): m_obj(obj), m_idx(0), m_fh(fh), m_ordered(ordered) {
   
   if(m_ordered) { //start at an arbitrary edge
      m_cur = EdgeHandle(m_obj.m_FE.getColByIndex(m_fh.idx(), m_idx));
   }

}

void FaceEdgeIterator::advance() {
   m_idx++;
   if(m_ordered)
      m_cur = m_obj.nextEdge(m_fh, m_cur);
}

bool FaceEdgeIterator::done() const {
   return m_idx >= 3;
}

EdgeHandle FaceEdgeIterator::current() const {
   if(m_ordered) {
      return m_idx >= 3? EdgeHandle::invalid() : m_cur;
   }
   else {
      return m_idx >= 3 ? EdgeHandle::invalid() : EdgeHandle(m_obj.m_FE.getColByIndex(m_fh.idx(), m_idx));
   }
}


//FaceTetIterator
FaceTetIterator::FaceTetIterator(const SimplicialComplex& obj, const FaceHandle& fh): m_obj(obj), m_idx(0), m_fh(fh) {
}

void FaceTetIterator::advance() {
   m_idx++;
}

bool FaceTetIterator::done() const {
   return m_idx >= (int)m_obj.m_FT.getNumEntriesInRow(m_fh.idx());
}

TetHandle FaceTetIterator::current() const {
   return m_idx >= (int)m_obj.m_FT.getNumEntriesInRow(m_fh.idx()) ?
      TetHandle::invalid():
   TetHandle(m_obj.m_FT.getColByIndex(m_fh.idx(), m_idx));
}

//TetFaceIterator
TetFaceIterator::TetFaceIterator(const SimplicialComplex& obj, const TetHandle& th): m_obj(obj), m_idx(0), m_th(th) {
}

void TetFaceIterator::advance() {
   m_idx++;
}

bool TetFaceIterator::done() const {
   return m_idx >= 4;
}

FaceHandle TetFaceIterator::current() const {
   return m_idx >= 4 ?
      FaceHandle::invalid():
   FaceHandle(m_obj.m_TF.getColByIndex(m_th.idx(), m_idx));
}

//VertexFaceIterator
VertexFaceIterator::VertexFaceIterator(const SimplicialComplex& obj, const VertexHandle& vh) : m_obj(obj) {
  
   //build the set of (unique) faces
   for(VertexEdgeIterator veit(m_obj, vh); !veit.done(); veit.advance())
      for(EdgeFaceIterator efit(m_obj, veit.current()); !efit.done(); efit.advance())
         m_faces.insert(efit.current());

   //...and set up an iterator
   m_fiter = m_faces.begin();
}


void VertexFaceIterator::advance() {
   if(!done())
      ++m_fiter;
}

bool VertexFaceIterator::done() const {
   return m_fiter == m_faces.end();
}

FaceHandle VertexFaceIterator::current() const {
   if(!done())
      return *m_fiter;
   else
      return FaceHandle::invalid();
}


//FaceVertexIterator
FaceVertexIterator::FaceVertexIterator(const SimplicialComplex& obj, const FaceHandle& fh, bool ordered) : 
   m_obj(obj), m_feit(obj, fh, ordered), m_fh(fh) 
{

}

void FaceVertexIterator::advance() {
   m_feit.advance();
}

bool FaceVertexIterator::done() const {
   return m_feit.done();
}

VertexHandle FaceVertexIterator::current() const {
   EdgeHandle curEdge = m_feit.current();
   int direction = m_obj.getRelativeOrientation(m_fh, curEdge);
   return direction > 0? m_obj.fromVertex(curEdge) : m_obj.toVertex(curEdge);
}

//VertexTetIterator
VertexTetIterator::VertexTetIterator(const SimplicialComplex& obj, const VertexHandle& vh) : m_obj(obj) {

   //build the set of (unique) tets
   for(VertexEdgeIterator veit(m_obj, vh); !veit.done(); veit.advance())
      for(EdgeFaceIterator efit(m_obj, veit.current()); !efit.done(); efit.advance())
         for(FaceTetIterator ftit(m_obj, efit.current()); !ftit.done(); ftit.advance())
            m_tets.insert(ftit.current());  

   //...and set up an iterator
   m_titer = m_tets.begin();
}


void VertexTetIterator::advance() {
   if(!done())
      ++m_titer;
}

bool VertexTetIterator::done() const {
   return m_titer == m_tets.end();
}

TetHandle VertexTetIterator::current() const {
   if(!done())
      return *m_titer;
   else
      return TetHandle::invalid();
}


//TetVertexIterator
TetVertexIterator::TetVertexIterator(const SimplicialComplex& obj, const TetHandle& th) : 
m_obj(obj)
{
   //build the set of (unique) vertices
   for(TetFaceIterator tfit(m_obj, th); !tfit.done(); tfit.advance())
      for(FaceEdgeIterator feit(m_obj, tfit.current(), false); !feit.done(); feit.advance())
         for(EdgeVertexIterator evit(m_obj, feit.current(), false); !evit.done(); evit.advance())
            m_verts.insert(evit.current());
   
   //...and set up an iterator
   m_viter = m_verts.begin();

}

void TetVertexIterator::advance() {
   if(!done())
      ++m_viter;
}

bool TetVertexIterator::done() const {
   return m_viter == m_verts.end();
}

VertexHandle TetVertexIterator::current() const {
   if(!done())
      return *m_viter;
   else
      return VertexHandle::invalid();
}


//EdgeTetIterator
EdgeTetIterator::EdgeTetIterator(const SimplicialComplex& obj, const EdgeHandle& eh) : m_obj(obj) {

   //build the set of (unique) tets
   for(EdgeFaceIterator efit(m_obj, eh); !efit.done(); efit.advance())
      for(FaceTetIterator ftit(m_obj, efit.current()); !ftit.done(); ftit.advance())
         m_tets.insert(ftit.current());  

   //...and set up an iterator
   m_titer = m_tets.begin();
}


void EdgeTetIterator::advance() {
   if(!done())
      ++m_titer;
}

bool EdgeTetIterator::done() const {
   return m_titer == m_tets.end();
}

TetHandle EdgeTetIterator::current() const {
   if(!done())
      return *m_titer;
   else
      return TetHandle::invalid();
}


//TetEdgeIterator
TetEdgeIterator::TetEdgeIterator(const SimplicialComplex& obj, const TetHandle& th) : 
m_obj(obj)
{
   //build the set of (unique) vertices
   for(TetFaceIterator tfit(m_obj, th); !tfit.done(); tfit.advance())
      for(FaceEdgeIterator feit(m_obj, tfit.current(), false); !feit.done(); feit.advance())
         m_edges.insert(feit.current());

   //...and set up an iterator
   m_eiter = m_edges.begin();

}

void TetEdgeIterator::advance() {
   if(!done())
      ++m_eiter;
}

bool TetEdgeIterator::done() const {
   return m_eiter == m_edges.end();
}

EdgeHandle TetEdgeIterator::current() const {
   if(!done())
      return *m_eiter;
   else
      return EdgeHandle::invalid();
}

}