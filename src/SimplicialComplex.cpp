#include "SimplicialComplex.h"
#include "SimplexProperty.h"
#include "SimplexIterators.h"

#include <iostream>
#include <utility>
#include <set>
#include <map>
#include <queue>

namespace SimplexMesh {

   int SimplicialComplex::numVerts() const { return m_nVerts;}

   int SimplicialComplex::numEdges() const { return m_nEdges;}

   int SimplicialComplex::numFaces() const { return m_nFaces;}

   int SimplicialComplex::numTets() const { return m_nTets;}
   
   SimplicialComplex::SimplicialComplex()
   {

      m_nVerts = 0;
      m_nEdges = 0;
      m_nFaces = 0;
      m_nTets = 0;

   }


   bool SimplicialComplex::vertexExists(const VertexHandle& vertex) const {
      if(vertex.idx() < 0 || vertex.idx() >= (int)m_V.size())
         return false;
      else
         return m_V[vertex.idx()];
   }

   bool SimplicialComplex::edgeExists(const EdgeHandle& edge) const {
      if(edge.idx() < 0 || edge.idx() >= (int)m_EV.getNumRows())
         return false;
      else
         return m_EV.getNumEntriesInRow(edge.idx()) > 0;

   }

   bool SimplicialComplex::faceExists(const FaceHandle& face) const {
      if(face.idx() < 0 || face.idx() >= (int)m_FE.getNumRows())
         return false;
      else
         return m_FE.getNumEntriesInRow(face.idx()) > 0;
   }

   bool SimplicialComplex::tetExists(const TetHandle& tet) const {
      if(tet.idx() < 0 || tet.idx() >= (int)m_TF.getNumRows())
         return false;
      else
         return m_TF.getNumEntriesInRow(tet.idx()) > 0;
   }

   
   int SimplicialComplex::getRelativeOrientation(const TetHandle& th, const FaceHandle& fh) const {
      if(th.idx() < 0 || th.idx() >= (int)m_TF.getNumRows() || fh.idx() < 0 || fh.idx() >= (int)m_TF.getNumCols())
         return 0;
      return m_TF.get(th.idx(), fh.idx());
   }

   int SimplicialComplex::getRelativeOrientation(const FaceHandle& fh, const EdgeHandle& eh) const {
      if(fh.idx() < 0 || fh.idx() >= (int)m_FE.getNumRows() || eh.idx() < 0 || eh.idx() >= (int)m_FE.getNumCols())
         return 0;
      return m_FE.get(fh.idx(), eh.idx());
   }

   int SimplicialComplex::getRelativeOrientation(const EdgeHandle& eh, const VertexHandle& vh) const {
      if(eh.idx() < 0 || eh.idx() >= (int)m_EV.getNumRows() || vh.idx() < 0 || vh.idx() >= (int)m_EV.getNumCols())
         return 0;
      return m_EV.get(eh.idx(), vh.idx());
   }

   
   bool SimplicialComplex::isIncident(const VertexHandle& vh, const EdgeHandle& eh) const {
     return m_EV.get(eh.idx(), vh.idx()) != 0;
   }

   bool SimplicialComplex::isIncident(const EdgeHandle& eh, const FaceHandle& fh) const {
     return m_FE.get(fh.idx(), eh.idx()) != 0;
   }

   bool SimplicialComplex::isIncident(const FaceHandle& fh, const TetHandle& th) const {
     return m_TF.get(th.idx(), fh.idx()) != 0;
   }

   
   VertexHandle SimplicialComplex::addVertex()
   {

      int new_index;
      if(m_deadVerts.size() == 0) {
         //create new vertex in incidence matrices
         m_EV.addCols(1);
         m_VE.addRows(1);

         //add the vertex to the data/property array
         m_V.push_back(true);
         for(unsigned int i = 0; i < m_vertProperties.size(); ++i) m_vertProperties[i]->resize(m_V.size());

         new_index = m_V.size()-1;
      }
      else {
         new_index = m_deadVerts.back();
         m_deadVerts.pop_back();

         assert(!m_V[new_index]);

         m_V[new_index] = true;
      }

      m_nVerts += 1;

      return VertexHandle(new_index);
   }


   EdgeHandle SimplicialComplex::addEdge(const VertexHandle& v0, const VertexHandle& v1)
   {

      //Cheap safety checks
      if(!vertexExists(v0) || !vertexExists(v1))
         return EdgeHandle::invalid();
      if(v0.idx() == v1.idx())
         return EdgeHandle::invalid();

      
      if(m_safetyChecks) {

         //full duplication: check if any edge joining these vertices already exists
         for(unsigned int i = 0; i < m_VE.getNumEntriesInRow(v0.idx()); ++i) {
            int edgeID = m_VE.getColByIndex(v0.idx(), i);
            for(unsigned int j = 0; j < m_EV.getNumEntriesInRow(edgeID); ++j){
               int vertID = m_EV.getColByIndex(edgeID, j);
               if(vertID == v1.idx()) //this creates a duplicate edge
                  return EdgeHandle::invalid();
            }
         }
      }
      

      //get the next free edge index, or add space
      int new_index;
      if(m_deadEdges.size() == 0) {
         //make room for the new edge
         m_FE.addCols(1);
         m_EF.addRows(1);

         m_EV.addRows(1);
         m_VE.addCols(1);

         new_index = m_EV.getNumRows()-1;

         //add property slots for the new edge
         for(unsigned int i = 0; i < m_edgeProperties.size(); ++i) m_edgeProperties[i]->resize(m_EV.getNumRows());

         assert(m_EV.getNumRows() == m_VE.getNumCols());
         assert(m_EV.getNumCols() == m_VE.getNumRows());
         assert(m_FE.getNumRows() == m_EF.getNumCols());
         assert(m_FE.getNumCols() == m_EF.getNumRows());
      }
      else {
         //grab the first dead edge off the pile
         new_index = m_deadEdges.back();
         m_deadEdges.pop_back();
      }

      //build new edge connectivity
      //choose indices explicitly, to indicate ordering
      m_EV.setByIndex(new_index, 0, v0.idx(), -1);
      m_EV.setByIndex(new_index, 1, v1.idx(), +1);
      
      m_VE.set(v0.idx(), new_index, -1);
      m_VE.set(v1.idx(), new_index, 1);

      //adjust edge count
      m_nEdges += 1;

      return EdgeHandle(new_index);
   }

   FaceHandle SimplicialComplex::addFace(const EdgeHandle& e0, 
      const EdgeHandle& e1, 
      const EdgeHandle& e2)
   {

      //cheap safety checks
      if(!edgeExists(e0) || !edgeExists(e1) || !edgeExists(e2)) //make sure all edges exists
         return FaceHandle::invalid();
      
      if(e0 == e1 || e1 == e2 || e0 == e2) //prevent degenerate faces
         return FaceHandle::invalid();

      
      if(m_safetyChecks) {

         //check if a duplicate face already exists (including partial matches and reversed orientations)
         for(unsigned int i = 0; i < m_EF.getNumEntriesInRow(e0.idx()); ++i) {
            int faceInd = m_EF.getColByIndex(e0.idx(), i);
            for(unsigned int j = 0; j < m_FE.getNumEntriesInRow(faceInd); ++j) {
               int edgeInd = m_FE.getColByIndex(faceInd, j);
               if(edgeInd == e1.idx() || edgeInd == e2.idx())
                  return FaceHandle::invalid();
            }
         }
         for(unsigned int i = 0; i < m_EF.getNumEntriesInRow(e1.idx()); ++i) {
            int faceInd = m_EF.getColByIndex(e1.idx(), i);
            for(unsigned int j = 0; j < m_FE.getNumEntriesInRow(faceInd); ++j) {
               int edgeInd = m_FE.getColByIndex(faceInd, j);
               if(edgeInd == e0.idx() || edgeInd == e2.idx())
                  return FaceHandle::invalid();
            }
         }
         for(unsigned int i = 0; i < m_EF.getNumEntriesInRow(e2.idx()); ++i) {
            int faceInd = m_EF.getColByIndex(e2.idx(), i);
            for(unsigned int j = 0; j < m_FE.getNumEntriesInRow(faceInd); ++j) {
               int edgeInd = m_FE.getColByIndex(faceInd, j);
               if(edgeInd == e0.idx() || edgeInd == e1.idx())
                  return FaceHandle::invalid();
            }
         }
        
         //check that the composing edges actually share the same 3 vertices, and use them twice each
         std::map<int,int> vertList;
         vertList[fromVertex(e0).idx()]++;
         vertList[toVertex(e0).idx()]++;
         vertList[fromVertex(e1).idx()]++;
         vertList[toVertex(e1).idx()]++;
         vertList[fromVertex(e2).idx()]++;
         vertList[toVertex(e2).idx()]++;
         if(vertList.size() != 3) 
            return FaceHandle::invalid();
         for(std::map<int,int>::iterator iter = vertList.begin(); iter != vertList.end(); ++iter) {
            if(iter->second != 2) {
               return FaceHandle::invalid();
            }
         }
      }

      //get the next free face or add one
      int new_index;
      if(m_deadFaces.size() == 0) {
         //create space for new face
         m_TF.addCols(1);
         m_FT.addRows(1);

         m_FE.addRows(1);
         m_EF.addCols(1);

         new_index = m_FE.getNumRows()-1;

         //allocate space for properties associated to the edge
         for(unsigned int i = 0; i < m_faceProperties.size(); ++i) m_faceProperties[i]->resize(m_FE.getNumRows());

         assert(m_FE.getNumRows() == m_EF.getNumCols());
         assert(m_FE.getNumCols() == m_EF.getNumRows());
         assert(m_FT.getNumRows() == m_TF.getNumCols());
         assert(m_FT.getNumCols() == m_TF.getNumRows());
      }
      else {
         //grab the next empty face off the pile
         new_index = m_deadFaces.back();
         m_deadFaces.pop_back();
      }

      //Signs are chosen to follow the ordering of edges provided as input,
      //so we flip to ensure edge vertices connect properly.

      //If the head of the first edge doesn't match either of the second edge's vertices, we must flip it,
      //since we want it oriented towards the second edge.
      bool flip0 = (toVertex(e0) != fromVertex(e1) && toVertex(e0) != toVertex(e1));

      //Now determine the shared vertex between edges 0 and 1.
      //Then flip edge1 if the shared_vertex is at the head; it should be at the tail.
      VertexHandle shared_vert0_hnd = flip0? fromVertex(e0) : toVertex(e0);
      bool flip1 = (shared_vert0_hnd != fromVertex(e1));

      //Determine shared vertex between edge 1 and 2.
      //Then flip edge2 if the shared_vertex is at the head.
      VertexHandle shared_vert1_hnd = flip1 ? fromVertex(e1) : toVertex(e1);
      bool flip2 = (shared_vert1_hnd != fromVertex(e2));

      //build face connectivity
      //get a unique ordering by arbitrarily choosing smallest index to go first
      int smallest = std::min(std::min(e0.idx(), e1.idx()), e2.idx());
      
      //add'em in requested ordering
      m_FE.setByIndex(new_index, 0, e0.idx(), flip0?-1:1);
      m_FE.setByIndex(new_index, 1, e1.idx(), flip1?-1:1);
      m_FE.setByIndex(new_index, 2, e2.idx(), flip2?-1:1);

      //cycle to get the smallest one first, for consistency
      while(m_FE.getColByIndex((unsigned int)new_index, (unsigned int)0) != (unsigned int)smallest)
        m_FE.cycleRow(new_index);

      //build the other one the usual way
      m_EF.set(e0.idx(), new_index, flip0?-1:1);
      m_EF.set(e1.idx(), new_index, flip1?-1:1);
      m_EF.set(e2.idx(), new_index, flip2?-1:1);

      m_nFaces += 1;

      return FaceHandle(new_index);
   }

   TetHandle SimplicialComplex::addTet(const FaceHandle& f0,
      const FaceHandle& f1,
      const FaceHandle& f2,
      const FaceHandle& f3, bool flip_face0)
   {

      //ensure all the faces actually exist
      if(!faceExists(f0) || !faceExists(f1) || !faceExists(f2) || !faceExists(f3))
         return TetHandle::invalid();

      //prevent degenerate tets
      if(f0 == f1 || f0 == f2 || f0 == f3 || f1 == f2 || f1 == f3 || f2 == f3)
         return TetHandle::invalid();

      if(m_safetyChecks) {
         
         //check for tets that share 2 (or possibly more) of the same faces, since this isn't really valid when embedded in 3D.
         //[technically perhaps sharing 3 faces should be the no-no]
         FaceHandle faceList[4] = {f0, f1, f2, f3};
         for(int i = 0; i < 4; ++i) {
            FaceHandle curF = faceList[i];
            for(unsigned int j = 0; j < m_FT.getNumEntriesInRow(curF.idx()); ++j) {
               unsigned tetId = m_FT.getColByIndex(curF.idx(), j);
               for(unsigned int k = 0; k < m_TF.getNumEntriesInRow(tetId); ++k) {
                  int otherF = m_TF.getColByIndex(tetId, k);
                  if(otherF == f0.idx() || otherF == f1.idx() || otherF == f2.idx() || otherF == f3.idx())
                     return TetHandle::invalid();
               }
            }
         }
         

         //Check that these faces share the same 6 edges, twice each
         FaceHandle faces[4] = {f0, f1, f2, f3};
         std::map<int,int> edgeList;
         for(int i = 0; i < 4; ++i) {
            for(unsigned int j = 0; j < m_FE.getNumEntriesInRow(faces[i].idx()); ++j) {
               int edgeInd = m_FE.getColByIndex(faces[i].idx(), j);
               edgeList[edgeInd]++;
            }
         }
         
         if(edgeList.size() != 6) 
            return TetHandle::invalid();
         for(std::map<int,int>::iterator iter = edgeList.begin(); iter != edgeList.end(); ++iter) {
            if(iter->second != 2) {
               return TetHandle::invalid();
            }
         }
      }

      //get the next free tet or add one
      int new_index;
      if(m_deadTets.size() == 0) {
         //create new tet
         m_TF.addRows(1);
         m_FT.addCols(1);

         new_index = m_TF.getNumRows()-1;

         //allocate space for the properties
         for(unsigned int i = 0; i < m_tetProperties.size(); ++i) m_tetProperties[i]->resize(m_TF.getNumRows());

         assert(m_FT.getNumRows() == m_TF.getNumCols());
         assert(m_FT.getNumCols() == m_TF.getNumRows());
      }
      else {
         //grab the next unused tet
         new_index = m_deadTets.back();
         m_deadTets.pop_back();
      }

      //Need to figure out signs to be consistent with the choice of the first face

      //Determine the shared edge between two adjacent faces. Flip the 2nd so its direction is opposed (signs differ) after the flips.
      EdgeHandle shared_edge0 = getSharedEdge(f0, f1);
      bool flip1 = flip_face0 ? m_FE.get(f0.idx(), shared_edge0.idx()) == m_FE.get(f1.idx(), shared_edge0.idx()) : 
         m_FE.get(f0.idx(), shared_edge0.idx()) != m_FE.get(f1.idx(), shared_edge0.idx());

      EdgeHandle shared_edge1 = getSharedEdge(f0, f2);
      bool flip2 = flip_face0 ? m_FE.get(f0.idx(), shared_edge1.idx()) == m_FE.get(f2.idx(), shared_edge1.idx()) : 
         m_FE.get(f0.idx(), shared_edge1.idx()) != m_FE.get(f2.idx(), shared_edge1.idx());

      EdgeHandle shared_edge2 = getSharedEdge(f0, f3);
      bool flip3 = flip_face0 ? m_FE.get(f0.idx(), shared_edge2.idx()) == m_FE.get(f3.idx(), shared_edge2.idx()) : 
         m_FE.get(f0.idx(), shared_edge2.idx()) != m_FE.get(f3.idx(), shared_edge1.idx());

      //build tet connectivity
      m_TF.set(new_index, f0.idx(), flip_face0?1:-1);
      m_TF.set(new_index, f1.idx(), flip1?1:-1);
      m_TF.set(new_index, f2.idx(), flip2?1:-1);
      m_TF.set(new_index, f3.idx(), flip3?1:-1);

      m_FT.set(f0.idx(), new_index, flip_face0?1:-1);
      m_FT.set(f1.idx(), new_index, flip1?1:-1);
      m_FT.set(f2.idx(), new_index, flip2?1:-1);
      m_FT.set(f3.idx(), new_index, flip3?1:-1);

      m_nTets += 1;

      //invalidate the relevant cached neighbour data

      return TetHandle(new_index);
   }

   FaceHandle SimplicialComplex::addFace(const VertexHandle& v0, const VertexHandle& v1, const VertexHandle& v2)
   {
      //This alternate method for adding a face requires 
      //searching for the edges containing the desired vertices, and creating them if they don't exist.
      //Likely to be slower, albeit more convenient.
      //ensure the vertices exist, and aren't duplicated
      if(!vertexExists(v0) || !vertexExists(v1) || !vertexExists(v2))
         return FaceHandle::invalid();

      //Careful about making the ordering match the desired input ordering of the vertices...

      //Find the edges we need, or create them.
      EdgeHandle e01 = getEdge(v0,v1);
      if(!e01.isValid()) e01 = addEdge(v0,v1);

      EdgeHandle e02 = getEdge(v2,v0);
      if(!e02.isValid()) e02 = addEdge(v2,v0);

      EdgeHandle e12 = getEdge(v1,v2);
      if(!e12.isValid()) e12 = addEdge(v1,v2);

      return addFace(e01, e12, e02);
   }

   TetHandle SimplicialComplex::addTet(const VertexHandle& v0, const VertexHandle& v1, const VertexHandle& v2, const VertexHandle& v3) {
      if(!v0.isValid() || !v1.isValid() || !v2.isValid() || !v3.isValid())
         return TetHandle::invalid();

      EdgeHandle e0 = getEdge(v0,v1);
      if(!e0.isValid()) e0 = addEdge(v0,v1);
      EdgeHandle e1 = getEdge(v0,v2);
      if(!e1.isValid()) e1 = addEdge(v0,v2);
      EdgeHandle e2 = getEdge(v0,v3);
      if(!e2.isValid()) e2 = addEdge(v0,v3);
      EdgeHandle e3 = getEdge(v1,v2);
      if(!e3.isValid()) e3 = addEdge(v1,v2);
      EdgeHandle e4 = getEdge(v1,v3);
      if(!e4.isValid()) e4 = addEdge(v1,v3);
      EdgeHandle e5 = getEdge(v2,v3);
      if(!e5.isValid()) e5 = addEdge(v2,v3);

      FaceHandle f0 = getFace(e0,e2,e4);
      if(!f0.isValid()) f0 = addFace(e0,e2,e4);
      FaceHandle f1 = getFace(e3,e4,e5);
      if(!f1.isValid()) f1 = addFace(e3,e4,e5);
      FaceHandle f2 = getFace(e0,e1,e3);
      if(!f2.isValid()) f2 = addFace(e0,e1,e3); 
      FaceHandle f3 = getFace(e1,e2,e5);
      if(!f3.isValid()) f3 = addFace(e1,e2,e5);

      return addTet(f0,f1,f2,f3);
   }



   bool SimplicialComplex::deleteVertex(const VertexHandle& vertex)
   {
      if(!vertexExists(vertex))
         return false;

      //test for safety - don't perform the delete if the simplex is not orphaned, or we get inconsistency.
      if(m_VE.getNumEntriesInRow(vertex.idx()) != 0)
         return false;

      //set the vertex to inactive
      m_V[vertex.idx()] = false;
      m_deadVerts.push_back(vertex.idx());

      //adjust the vertex count
      m_nVerts -= 1;

      return true;
   }

   bool SimplicialComplex::deleteEdge(const EdgeHandle& edge, bool recurse)
   {
      if(!edgeExists(edge))
         return false;

      //test for safety - don't perform delete if not orphaned, or it creates inconsistency.
      if(m_EF.getNumEntriesInRow(edge.idx()) != 0)
         return false;

      //determine the corresponding vertices
      for(unsigned int i = 0; i < m_EV.getNumEntriesInRow(edge.idx()); ++i) {

         //delete the edge entry in the transpose
         int col = m_EV.getColByIndex(edge.idx(),i);
         m_VE.remove(col, edge.idx()); 

         //delete the composing vertices if desired
         if(recurse)
            deleteVertex(VertexHandle(col));       
      }

      //...and delete the row
      m_EV.zeroRow(edge.idx());
      m_deadEdges.push_back(edge.idx());

      //adjust the edge count
      m_nEdges -= 1;

      return true;
   }

   bool SimplicialComplex::deleteFace(const FaceHandle& face, bool recurse)
   {
      if(!faceExists(face))
         return false;

      //test for safety - don't perform delete if not orphaned, or we'll have inconsistencies.
      if(m_FT.getNumEntriesInRow(face.idx()) != 0)
         return false;

      //determine the corresponding edges
      for(unsigned int i = 0; i < m_FE.getNumEntriesInRow(face.idx()); ++i) {

         //remove face entry from the transpose
         int col = m_FE.getColByIndex(face.idx(), i);
         m_EF.remove(col, face.idx());

         //delete the composing edges
         if(recurse)
            deleteEdge(EdgeHandle(col), recurse);
      }

      //...and delete the row
      m_FE.zeroRow(face.idx());
      m_deadFaces.push_back(face.idx());

      //adjust the face count
      m_nFaces -= 1;

      return true;
   }

   bool SimplicialComplex::deleteTet(const TetHandle& tet, bool recurse) {
      if(!tetExists(tet))
         return false;

      //There are no higher dimensional simplices in 3D, so deleting the tet cannot introduce inconsistencies
      //as it can in the lower cases.

      //determine the corresponding faces
      for(unsigned int i = 0; i < m_TF.getNumEntriesInRow(tet.idx()); ++i) {

         //clear the tet entries in the transpose
         int col = m_TF.getColByIndex(tet.idx(), i);
         m_FT.remove(col, tet.idx());

         //delete composing faces if desired
         if(recurse)
            deleteFace(FaceHandle(col), recurse);
      }

      //...and delete the row
      m_TF.zeroRow(tet.idx());
      m_deadTets.push_back(tet.idx());

      //adjust the tet count
      m_nTets -= 1;

      //invalidate cached relationships

      return true;

   }


   
   VertexHandle SimplicialComplex::getVertex(const EdgeHandle& edge, int index) const {
      assert(edgeExists(edge));
      assert(index >= 0 && index <= 1);
      return VertexHandle(m_EV.getColByIndex(edge.idx(), index));
   }

   EdgeHandle SimplicialComplex::getEdge(const FaceHandle& face, int index) const {
      assert(faceExists(face));
      assert(index >= 0 && index <= 2);
      return EdgeHandle(m_FE.getColByIndex(face.idx(), index));
   }

   FaceHandle SimplicialComplex::getFace(const TetHandle& tet, int index) const {
      assert(tetExists(tet));
      assert(index >= 0 && index <= 3);
      return FaceHandle(m_TF.getColByIndex(tet.idx(), index));
   }

   EdgeHandle SimplicialComplex::getEdge(const VertexHandle& v0, const VertexHandle& v1) const {
      //returns the appropriate edge if it exists, ignoring orientation
      if(!vertexExists(v0) || !vertexExists(v1))
         return EdgeHandle::invalid();

      for(VertexEdgeIterator veit(*this, v0); !veit.done(); veit.advance()) {
         EdgeHandle curEdge = veit.current();
         if(fromVertex(curEdge) == v1 || toVertex(curEdge) == v1)
            return curEdge;
      }

      return EdgeHandle::invalid();
   }

   FaceHandle SimplicialComplex::getFace(const EdgeHandle& e0, const EdgeHandle& e1, const EdgeHandle& e2) const {
      //returns the appropriate face if it exists, ignoring orientation
      if(!edgeExists(e0) || !edgeExists(e1) || !edgeExists(e2)) 
         return FaceHandle::invalid();

      for(EdgeFaceIterator efit(*this, e0); !efit.done(); efit.advance()) {
         FaceHandle curFace = efit.current();
         bool foundE1 = false, foundE2 = false;
         for(FaceEdgeIterator feit(*this, curFace, false); !feit.done(); feit.advance()) {
            EdgeHandle curEdge = feit.current();
            if(curEdge == e1) foundE1 = true;
            if(curEdge == e2) foundE2 = true;
         }
         if(foundE1 && foundE2)
            return curFace;
      }

      return FaceHandle::invalid();
   }

   TetHandle SimplicialComplex::getTet(const FaceHandle& f0, const FaceHandle& f1, const FaceHandle& f2, const FaceHandle& f3) const {
      //returns the appropriate tet if it exists, ignoring orientation
      if(!faceExists(f0) || !faceExists(f1) || !faceExists(f2) || !faceExists(f3)) 
         return TetHandle::invalid();

      for(FaceTetIterator ftit(*this, f0); !ftit.done(); ftit.advance()) {
         TetHandle curTet = ftit.current();
         bool foundf1 = false, foundf2 = false, foundf3 = false;
         for(TetFaceIterator tfit(*this, curTet); !tfit.done(); tfit.advance()) {
            FaceHandle curFace = tfit.current();
            if(curFace == f1) foundf1 = true;
            if(curFace == f2) foundf2 = true;
            if(curFace == f3) foundf3 = true;
         }
         if(foundf1 && foundf2 && foundf3)
            return curTet;
      }

      return TetHandle::invalid();
   }

   
   EdgeHandle SimplicialComplex::getSharedEdge(const FaceHandle& f0, const FaceHandle& f1) const {
      assert(faceExists(f0) && faceExists(f1));
      assert(m_FE.getNumEntriesInRow(f0.idx()) == 3);
      assert(m_FE.getNumEntriesInRow(f1.idx()) == 3);

      //Iterate over all pairs, stop when we hit the match.
      for(unsigned int ind0 = 0; ind0 < 3; ++ind0) for(unsigned ind1 = 0; ind1 < 3; ++ind1) {
         int col0 = m_FE.getColByIndex(f0.idx(), ind0);
         int col1 = m_FE.getColByIndex(f1.idx(), ind1);
         if(col0 == col1)
            return EdgeHandle(col0);
      }

      return EdgeHandle(-1);
   }

   
   VertexHandle SimplicialComplex::fromVertex(const EdgeHandle& eh) const
   {
      assert(eh.isValid());
      assert(m_EV.getNumEntriesInRow(eh.idx()) == 2);

      //1st vertex is the from vertex, by design
      return VertexHandle(m_EV.getColByIndex(eh.idx(), 0));
   }

   VertexHandle SimplicialComplex::toVertex(const EdgeHandle& eh) const
   {
      assert(eh.isValid());
      assert(m_EV.getNumEntriesInRow(eh.idx()) == 2);

      //2nd vertex is the to vertex, by design
      return VertexHandle(m_EV.getColByIndex(eh.idx(), 1));
   }

   FaceHandle SimplicialComplex::frontFace(const EdgeHandle& eh) const {
      assert(eh.isValid());
      assert(m_EF.getNumEntriesInRow(eh.idx()) <= 2);

      int edgeIdx = eh.idx(); 

      if(m_EF.getNumCols() == 0)
         return FaceHandle::invalid();
      else if(m_EF.getNumCols() == 1)
         return (m_EF.getValueByIndex(edgeIdx, 0) == 1 ? FaceHandle(m_EF.getColByIndex(edgeIdx, 0)) : FaceHandle::invalid());
      else
         return (m_EF.getValueByIndex(edgeIdx, 0) == 1 ? FaceHandle(m_EF.getColByIndex(edgeIdx, 0)) : FaceHandle(m_EF.getColByIndex(edgeIdx, 1)));
   }

   FaceHandle SimplicialComplex::backFace(const EdgeHandle& eh) const {
      assert(eh.isValid());
      assert(m_EF.getNumEntriesInRow(eh.idx()) <= 2);

      int edgeIdx = eh.idx();
      if(m_EF.getNumCols() == 0)
         return FaceHandle::invalid();
      else if(m_EF.getNumCols() == 1)
         return (m_EF.getValueByIndex(edgeIdx, 0) == -1 ? FaceHandle(m_EF.getColByIndex(edgeIdx, 0)) : FaceHandle::invalid());
      else
         return (m_EF.getValueByIndex(edgeIdx, 0) == -1 ? FaceHandle(m_EF.getColByIndex(edgeIdx, 0)) : FaceHandle(m_EF.getColByIndex(edgeIdx, 1)));
   }

   TetHandle SimplicialComplex::frontTet(const FaceHandle& fh) const {
      assert(fh.isValid());
      assert(m_FT.getNumEntriesInRow(fh.idx()) <= 2);

      int faceIdx = fh.idx();
      if(m_FT.getNumCols() == 0)
         return TetHandle::invalid();
      else if(m_FT.getNumCols() == 1)
         return (m_FT.getValueByIndex(faceIdx, 0) == 1 ? TetHandle(m_FT.getColByIndex(faceIdx, 0)) : TetHandle::invalid());
      else
         return (m_FT.getValueByIndex(faceIdx, 0) == 1 ? TetHandle(m_FT.getColByIndex(faceIdx, 0)) : TetHandle(m_FT.getColByIndex(faceIdx, 1)));
   }

   TetHandle SimplicialComplex::backTet(const FaceHandle& fh) const {
      assert(fh.isValid());
      assert(m_FT.getNumEntriesInRow(fh.idx()) <= 2);

      int faceIdx = fh.idx();
      if(m_FT.getNumCols() == 0)
         return TetHandle::invalid();
      else if(m_FT.getNumCols() == 1)
         return (m_FT.getValueByIndex(faceIdx, 0) == -1 ? TetHandle(m_FT.getColByIndex(faceIdx, 0)) : TetHandle::invalid());
      else
         return (m_FT.getValueByIndex(faceIdx, 0) == -1 ? TetHandle(m_FT.getColByIndex(faceIdx, 0)) : TetHandle(m_FT.getColByIndex(faceIdx, 1)));
   }

  
   VertexHandle SimplicialComplex::nextVertex(const VertexHandle& curVertex) const {
      assert(vertexExists(curVertex));

      int idx = curVertex.idx();
      do {
         ++idx;
         idx %= numVertexSlots();
      } while(!m_V[idx]);

      return VertexHandle(idx);
   }
   
   VertexHandle SimplicialComplex::prevVertex(const VertexHandle& curVertex) const {
      assert(vertexExists(curVertex));
   
      int idx = curVertex.idx();
      do {
         idx+=numVertexSlots()-1;
         idx %= numVertexSlots();
      } while(!m_V[idx]);

      return VertexHandle(idx);
   }

   EdgeHandle SimplicialComplex::nextEdge(const EdgeHandle& curEdge) const {
      assert(edgeExists(curEdge));

      int idx = curEdge.idx();
      do {
         ++idx;
         idx %= numEdgeSlots();
      }while (m_EV.getNumEntriesInRow(idx) == 0);

      return EdgeHandle(idx);
   }

   EdgeHandle SimplicialComplex::prevEdge(const EdgeHandle& curEdge) const {
      assert(edgeExists(curEdge));

      int idx = curEdge.idx();
      do {
         idx += numEdgeSlots()-1;
         idx %= numEdgeSlots();
      }while (m_EV.getNumEntriesInRow(idx) == 0);

      return EdgeHandle(idx);
   }

   FaceHandle SimplicialComplex::nextFace(const FaceHandle& curFace) const {
      assert(faceExists(curFace));

      int idx = curFace.idx();
      do {
         ++idx;
         idx %= numFaceSlots();
      }while (m_FE.getNumEntriesInRow(idx) == 0);

      return FaceHandle(idx);
   }

   FaceHandle SimplicialComplex::prevFace(const FaceHandle& curFace) const {
      assert(faceExists(curFace));

      int idx = curFace.idx();
      do {
         idx += numFaceSlots()-1;
         idx %= numFaceSlots();
      } while (m_FE.getNumEntriesInRow(idx) == 0);

      return FaceHandle(idx);
   }

   TetHandle SimplicialComplex::nextTet(const TetHandle& curTet) const {
      assert(tetExists(curTet));

      int idx = curTet.idx();
      do {
         ++idx;
         idx %= numTetSlots();
      }while (m_TF.getNumEntriesInRow(idx) == 0);

      return TetHandle(idx);
   }

   TetHandle SimplicialComplex::prevTet(const TetHandle& curTet) const {
      assert(tetExists(curTet));

      int idx = curTet.idx();
      do {
         idx += numTetSlots()-1;
         idx %= numTetSlots();
      }while (m_TF.getNumEntriesInRow(idx) == 0);

      return TetHandle(idx);
   }

   //---------------------------------------

   VertexHandle SimplicialComplex::nextVertex(const EdgeHandle& edge, const VertexHandle& curVert) const {
      assert(edgeExists(edge));
      assert(vertexExists(curVert));
      assert(m_EV.getNumEntriesInRow(edge.idx()) == 3);
      assert(isIncident(curVert, edge));

      //just return whichever vert is the opposite
      int edgeIdx = edge.idx();
      int col0 = m_EV.getColByIndex(edgeIdx, 0);
      int vertIdx = (col0 == curVert.idx() ? m_EV.getColByIndex(edgeIdx, 1) : col0);
      
      return VertexHandle(vertIdx);
   }

   VertexHandle SimplicialComplex::prevVertex(const EdgeHandle& edge, const VertexHandle& curVert) const {
      return nextVertex(edge, curVert);
   }

   
   EdgeHandle SimplicialComplex::nextEdge(const FaceHandle& face, const EdgeHandle& curEdge) const {
      assert(faceExists(face));
      assert(edgeExists(curEdge));
      assert(m_FE.getNumEntriesInRow(face.idx()) == 3);
      assert(isIncident(curEdge, face));

      //take advantage of known fixed ordering of indices
      int faceIdx = face.idx();
      int colIdx = curEdge.idx();
      int col0 = m_FE.getColByIndex(faceIdx, 0);
      int col1 = m_FE.getColByIndex(faceIdx, 1);
      if(col0 == colIdx)
         return EdgeHandle(col1);
      else if(col1 == colIdx)
         return EdgeHandle(m_FE.getColByIndex(faceIdx, 2));
      else 
         return EdgeHandle(col0);
   }

   EdgeHandle SimplicialComplex::prevEdge(const FaceHandle& face, const EdgeHandle& curEdge) const {
      assert(faceExists(face));
      assert(edgeExists(curEdge));
      assert(m_FE.getNumEntriesInRow(face.idx()) == 3);

      //take advantage of known fixed ordering of indices
      int faceIdx = face.idx();
      int colIdx = curEdge.idx();
      int col0 = m_FE.getColByIndex(faceIdx, 0);
      int col2 = m_FE.getColByIndex(faceIdx, 2);
      if(col0 == colIdx)
         return EdgeHandle(col2);
      else if(col2 == colIdx)
         return EdgeHandle(col0);
      else 
         return EdgeHandle(m_FE.getColByIndex(faceIdx, 1));
   }

   FaceHandle SimplicialComplex::nextFace(const TetHandle& tet, const FaceHandle& curFace) const {
      assert(tetExists(tet));
      assert(faceExists(curFace));
      assert(m_TF.getNumEntriesInRow(tet.idx()) == 4);
      assert(isIncident(curFace, tet));

      int tetIdx = tet.idx();

      //find the current face
      int colIdx = 0;
      for(; colIdx < 4; ++colIdx) {
         if(curFace.idx() == (int)m_TF.getColByIndex(tetIdx, colIdx)) {
            break;
         }
      }
      int nextIdx = (colIdx+1)%4;
      return FaceHandle(m_TF.getColByIndex(tetIdx, nextIdx));
   }

   FaceHandle SimplicialComplex::prevFace(const TetHandle& tet, const FaceHandle& curFace) const {
      assert(tetExists(tet));
      assert(faceExists(curFace));
      assert(m_TF.getNumEntriesInRow(tet.idx()) == 4);
      assert(isIncident(curFace, tet));

      int tetIdx = tet.idx();

      //find the current face
      int colIdx = 0;
      for(; colIdx < 4; ++colIdx) {
         if(curFace.idx() == (int)m_TF.getColByIndex(tetIdx, colIdx)) {
            break;
         }
      }
      int nextIdx = (colIdx+3)%4;
      return FaceHandle(m_TF.getColByIndex(tetIdx, nextIdx));
   }

   //--------------------------------

   
   VertexHandle SimplicialComplex::collapseEdge(const EdgeHandle& eh, const VertexHandle& vertToRemove) {

      VertexHandle fromV = fromVertex(eh);
      VertexHandle toV = toVertex(eh);
      assert(fromV == vertToRemove || toV == vertToRemove); //make sure the vertex selected for deletion is actually used by the edge

      //grab the vertex that we are keeping
      int vertToKeep = fromV == vertToRemove? toV.idx() : fromV.idx();

      //determine adjacent faces to the collapsing edge
      int faceCount = m_EF.getNumEntriesInRow(eh.idx());

      //look at faces on left and right of collapsing edge
      std::vector<int> facesToDelete(faceCount);
      for(int f = 0; f < faceCount; ++f) {
         int face_idx = m_EF.getColByIndex(eh.idx(), f);

         //visit the face's other edges
         std::set<int> neighbourEdges;
         for(unsigned int i = 0; i < m_FE.getNumEntriesInRow(face_idx); ++i) {
            int edge_idx = m_FE.getColByIndex(face_idx, i);
            if(edge_idx != eh.idx()) {

               //walk over the other faces that this edge belongs to, checking for a shared edge
               for(unsigned int j = 0; j < m_EF.getNumEntriesInRow(edge_idx); ++j) {
                  int curFace = m_EF.getColByIndex(edge_idx, j);

                  //if we see a shared edge, then the collapsing edge merges two faces
                  //and that's unacceptable.
                  if(curFace != face_idx) {
                     for(unsigned int k = 0; k < m_FE.getNumEntriesInRow(curFace); ++k) {
                        int curEdge = m_FE.getColByIndex(curFace, k);

                        if(neighbourEdges.find(curEdge) != neighbourEdges.end()) {
                           return VertexHandle(-1); //there's a shared face/edge, don't collapse
                        }
                        else {
                           neighbourEdges.insert(curEdge);
                        }
                     }
                  }
               }
            }
         }

         facesToDelete[f] = face_idx;
      }


      //delete the faces and then edge, leaving a hole to be stitched
      for(unsigned int i = 0; i < facesToDelete.size(); ++i) {
         bool success = deleteFace(FaceHandle(facesToDelete[i]), false);
         assert(success);
      }
      bool success = deleteEdge(eh, false);
      assert(success);

      //determine all existing edges using the vertex being eliminated
      std::vector< std::pair<int,int> > edgeIndices;
      for(unsigned int e = 0; e < m_VE.getNumEntriesInRow(vertToRemove.idx()); ++e) {
         unsigned int edgeInd = m_VE.getColByIndex(vertToRemove.idx(), e);
         int sign = m_VE.getValueByIndex(vertToRemove.idx(), e);
         edgeIndices.push_back(std::make_pair(edgeInd,sign));
      }

      //relabel all the edges' to-be-deleted endpoints to the vertex being kept.
      //doing it "in place" like this rather than using safer atomic add/deletes
      //ensures that the original data on the modified edges gets maintained.
      for(unsigned int i = 0; i < edgeIndices.size(); ++i) {
         unsigned int edgeInd = edgeIndices[i].first;
         int vertSign = edgeIndices[i].second;

         m_VE.remove(vertToRemove.idx(), edgeInd);
         m_EV.remove(edgeInd, vertToRemove.idx());

         m_VE.set(vertToKeep, edgeInd, vertSign);
         m_EV.set(edgeInd, vertToKeep, vertSign);
      }


      //now we have some edges that are duplicates, possibly pointing in opposite directions
      //identify duplicate edges for deletion.
      std::vector< std::pair<int,int> > duplicateEdges; //list of (edge,edge) pairs that are duplicates
      std::map<int,int> vertEdgeMap;//for each vertex in the set of neighbours, the first edge we hit that uses it.
      for(unsigned int e = 0; e < m_VE.getNumEntriesInRow(vertToKeep); ++e) {
         int edgeInd = m_VE.getColByIndex(vertToKeep,e);
         int fromV = fromVertex(EdgeHandle(edgeInd)).idx();
         int toV = toVertex(EdgeHandle(edgeInd)).idx();
         int otherVert = fromV == vertToKeep? toV : fromV;

         //check if an edge with this end vertex has been seen yet
         //if not, add to the set; if so, log it as a duplicate
         std::map<int,int>::iterator iter = vertEdgeMap.find(otherVert);
         if(iter != vertEdgeMap.end())
            duplicateEdges.push_back(std::make_pair(edgeInd,(*iter).second));
         else
            vertEdgeMap[otherVert] = edgeInd;
      }


      //now replace the duplicate edges with their partner everywhere they're used (relabelling again)
      for(unsigned int i = 0; i < duplicateEdges.size(); ++i) {
         EdgeHandle e0(duplicateEdges[i].first);
         EdgeHandle e1(duplicateEdges[i].second);

         //if the edges pointed in opposite directions, their directions in faces need to be swapped during the relabelling
         int flipSign = fromVertex(e0) != fromVertex(e1) ? -1 : 1;

         //let's choose to remove e1 arbitrarily.

         //collect all faces that use this edge
         std::vector< std::pair<int,int> > faceIndices;
         for(unsigned int f = 0; f < m_EF.getNumEntriesInRow(e1.idx()); ++f) {
            unsigned int faceInd = m_EF.getColByIndex(e1.idx(), f);
            int sign = m_EF.getValueByIndex(e1.idx(), f);
            faceIndices.push_back(std::make_pair(faceInd,sign));
         }

         //relabel all the faces' to-be-deleted edge to the other duplicate edge being kept.
         //doing it "in place" like this rather than using safer atomic add/deletes
         //ensures that the original data on the retained edges gets maintained.
         for(unsigned int i = 0; i < faceIndices.size(); ++i) {
            unsigned int faceInd = faceIndices[i].first;
            int edgeSign = faceIndices[i].second;

            int newSign = flipSign*edgeSign;
            m_EF.remove(e1.idx(), faceInd);
            m_FE.remove(faceInd, e1.idx());

            m_EF.set(e0.idx(), faceInd, newSign);
            m_FE.set(faceInd, e0.idx(), newSign);
         }

         //finally, delete the orphaned edge
         bool success = deleteEdge(e1, false);
         assert(success);
      }

      success = deleteVertex(vertToRemove);
      assert(success);

      return VertexHandle(vertToKeep);
   }

   //Two silly utility functions
   VertexHandle getSharedVertexFromEdgePair(const SimplicialComplex& obj, const EdgeHandle& e0, const EdgeHandle& e1) {
      VertexHandle v0 = obj.fromVertex(e0), v1 = obj.toVertex(e0),
         v2 = obj.fromVertex(e1), v3 = obj.toVertex(e1);

      if(v0 == v2 || v0 == v3) return v0;
      else if(v1 == v2 || v1 == v3) return v1;
      else return VertexHandle::invalid();
   }

   EdgeHandle getEdgeFromVertexPair( const SimplicialComplex& obj, const VertexHandle& v0, const VertexHandle& v1 ) {
      for(VertexEdgeIterator ve_it(obj, v0); !ve_it.done(); ve_it.advance()) {
         EdgeHandle eh = ve_it.current();
         if(obj.fromVertex(eh) == v1 || obj.toVertex(eh) == v1)
            return eh;
      }
      return EdgeHandle::invalid();
   }

   VertexHandle SimplicialComplex::splitEdge(const EdgeHandle& splitEdge, std::vector<FaceHandle>& newFaces) {

      EdgeFaceIterator ef_iter(*this, splitEdge);

      newFaces.clear();

      //get the edge's vertices
      VertexHandle from_vh, to_vh;
      from_vh = fromVertex(splitEdge);
      to_vh = toVertex(splitEdge);

      //add a new midpoint vertex
      VertexHandle newVert = addVertex();

      //add the two new edges that are part of the original split edge
      EdgeHandle e_0 = addEdge(from_vh, newVert);
      EdgeHandle e_1 = addEdge(to_vh, newVert);

      //now iterate over the existing faces, splitting them in two appropriately
      std::vector<FaceHandle> facesToDelete;
      EdgeFaceIterator ef_it(*this, splitEdge);
      for(;!ef_it.done(); ef_it.advance()) {
         FaceHandle fh = ef_it.current();

         //store this for deletion later.
         facesToDelete.push_back(fh);

         //find the other vertex in the first face
         FaceVertexIterator fv_it(*this, fh);
         while((fv_it.current() == from_vh) || (fv_it.current() == to_vh)) fv_it.advance();
         VertexHandle other_vh = fv_it.current();

         //create the new edge that splits this face
         EdgeHandle e_faceSplit = addEdge(other_vh, newVert);

         //build the two new faces
         FaceEdgeIterator fe_it(*this, fh);
         for(;!fe_it.done(); fe_it.advance()) {
            EdgeHandle cur = fe_it.current();

            //for each edge that isn't the splitEdge...      
            if(cur == splitEdge) continue;

            //accumulate a list of the new edges for the associated new face.
            //it is done in this manner to ensure new orientation is consistent with the old
            std::vector<EdgeHandle> edgeList;

            //determine which half of the splitEdge to use for this new face
            EdgeHandle halfEdge = fromVertex(cur) == from_vh || toVertex(cur) == from_vh? e_0 : e_1;

            FaceEdgeIterator fe_it2(*this, fh);
            for(;!fe_it2.done(); fe_it2.advance()) {
               EdgeHandle cur2 = fe_it2.current();
               if(cur2 == cur) edgeList.push_back(cur); //the current edge
               else if(cur2 == splitEdge) edgeList.push_back(halfEdge); //half of the original split edge
               else edgeList.push_back(e_faceSplit); //the new edge that splits the face
            }
            FaceHandle newFace = addFace(edgeList[0], edgeList[1], edgeList[2]);
            newFaces.push_back(newFace);
         }

      }

      //Delete the previous faces and edge. Done as a post-process so as not to mess up the iterators.
      for(unsigned int i = 0; i < facesToDelete.size(); ++i)
         deleteFace(facesToDelete[i], false);
      deleteEdge(splitEdge, false);

      //Return a handle to the vertex we created
      return newVert;
   }

   EdgeHandle SimplicialComplex::flipEdge(const EdgeHandle& eh) {
      assert(edgeExists(eh));

      VertexHandle from_vh, to_vh;
      from_vh = fromVertex(eh);
      to_vh = toVertex(eh);

      assert(from_vh != to_vh);

      //Just use the first two faces we hit.
      EdgeFaceIterator ef_it(*this, eh);
      if(ef_it.done()) return EdgeHandle::invalid(); //if there is no face at all, can't flip
      const FaceHandle fh = ef_it.current();
      ef_it.advance();
      if(ef_it.done()) return EdgeHandle::invalid(); //there is no 2nd face, we also can't flip
      const FaceHandle fh2 = ef_it.current();
      ef_it.advance();
      assert(fh != fh2); //we should never hit the same face
      assert(ef_it.done()); //this edge should have no more faces. A non-manifold edge flip doesn't make sense.

      //find the 3rd vertex in the first face
      FaceVertexIterator fv_it(*this, fh);
      while((fv_it.current() == from_vh) || (fv_it.current() == to_vh)) fv_it.advance();
      VertexHandle f1_vh = fv_it.current();

      //find the 3rd vertex in the second face
      FaceVertexIterator fv_it2(*this, fh2);
      while((fv_it2.current() == from_vh) || (fv_it2.current() == to_vh)) fv_it2.advance();
      VertexHandle f2_vh = fv_it2.current();

      assert(f1_vh != f2_vh);

      assert(from_vh != f1_vh);
      assert(from_vh != f2_vh);
      assert(to_vh != f1_vh);
      assert(to_vh != f2_vh);

      //check for an edge already matching this description... and don't do the flip.
      EdgeHandle edge = getEdgeFromVertexPair(*this, f1_vh, f2_vh);
      if(edge.isValid()) return EdgeHandle::invalid();

      EdgeHandle newEdge = addEdge(f1_vh, f2_vh);

      //grab all the current edges in proper order, starting from the shared face
      EdgeHandle e0 = nextEdge(fh, eh);
      EdgeHandle e1 = nextEdge(fh, e0);

      EdgeHandle e2 = nextEdge(fh2, eh);
      EdgeHandle e3 = nextEdge(fh2, e2);

      assert(e0 != e1);
      assert(e0 != e2);
      assert(e0 != e3);
      assert(e0 != eh);
      assert(e0 != newEdge);

      assert(e1 != e2);
      assert(e1 != e3);
      assert(e1 != eh);
      assert(e1 != newEdge);

      assert(e2 != e3);
      assert(e2 != eh);
      assert(e2 != newEdge);

      assert(e3 != eh);
      assert(e3 != newEdge);

      assert(eh != newEdge);

      //flip the edges of the second face so it matches the first face 
      //(if the original faces didn't sync, it doesn't matter, since there's no way to flip and maintain consistency)
      VertexHandle sharedV = getSharedVertexFromEdgePair(*this, e1, e2);
      if(!sharedV.isValid())
         std::swap(e2,e3);

      //add in the new faces
      addFace(e1, e2, newEdge);
      addFace(e3, e0, newEdge);

      //Delete the old patch
      bool success = deleteFace(fh, false);
      assert(success);
      success = deleteFace(fh2, false);
      assert(success);
      success = deleteEdge(eh, false);
      assert(success);

      return newEdge;
   }

   
   
   void SimplicialComplex::registerVertexProperty(SimplexPropertyBase* prop) { 
      m_vertProperties.push_back(prop); 
      prop->resize(numVertexSlots());
   }
   void SimplicialComplex::removeVertexProperty(SimplexPropertyBase* prop) { 
      std::vector<SimplexPropertyBase*>::iterator it = std::find(m_vertProperties.begin(), m_vertProperties.end(), prop);
      if(it != m_vertProperties.end()) 
         m_vertProperties.erase(it); 
   }

   
   void SimplicialComplex::registerEdgeProperty(SimplexPropertyBase* prop) { 
      m_edgeProperties.push_back(prop); 
      prop->resize(numEdgeSlots());
   }
   void SimplicialComplex::removeEdgeProperty(SimplexPropertyBase* prop) { 
      std::vector<SimplexPropertyBase*>::iterator it = std::find(m_edgeProperties.begin(), m_edgeProperties.end(), prop);
      if(it != m_edgeProperties.end()) 
         m_edgeProperties.erase(it); 
   }

   
   void SimplicialComplex::registerFaceProperty(SimplexPropertyBase* prop) { 
      m_faceProperties.push_back(prop); 
      prop->resize(numFaceSlots());
   }
   void SimplicialComplex::removeFaceProperty(SimplexPropertyBase* prop) { 
      std::vector<SimplexPropertyBase*>::iterator it = std::find(m_faceProperties.begin(), m_faceProperties.end(), prop);
      if(it != m_faceProperties.end()) 
         m_faceProperties.erase(it); 
   }

   
   void SimplicialComplex::registerTetProperty(SimplexPropertyBase* prop) { 
      m_tetProperties.push_back(prop); 
      prop->resize(numTetSlots());
   }
   void SimplicialComplex::removeTetProperty(SimplexPropertyBase* prop) { 
      std::vector<SimplexPropertyBase*>::iterator it = std::find(m_tetProperties.begin(), m_tetProperties.end(), prop);
      if(it != m_tetProperties.end()) 
         m_tetProperties.erase(it); 
   }


   
   bool SimplicialComplex::isOnBoundary(const FaceHandle& fh) const {
      int numTets = faceIncidentTetCount(fh);

      //in 3D, it is boundary if it has only one tet
      //in 2D is is never on the boundary (is this true?)
      //in 1D it cannot exist
      return numTets == 1;
   }

   bool SimplicialComplex::isOnBoundary(const EdgeHandle& eh) const {

      //Ignoring mixed dimensional points here...

      //it's boundary if at least one of the tets is a boundary tet
      bool partOfAnyTets = false;
      for(EdgeFaceIterator efit(*this, eh); !efit.done(); efit.advance()) {
         FaceHandle fh = efit.current();
         int tets = faceIncidentTetCount(fh);
         if(tets > 0)
            partOfAnyTets = true;
         if(tets == 1)
            return true; //it's 
      }

      if(partOfAnyTets) //all the tets are connected, there's no boundary.
         return false;

      //if it's part of a face, it's on the boundary if it has only 1 face
      int numFaces = edgeIncidentFaceCount(eh);
      if(numFaces > 0)
         return numFaces == 1;
      else //if it's just an edge, it's never part of the boundary
         return false;

   }

   bool SimplicialComplex::isOnBoundary(const VertexHandle& vh) const {

      //ignore mixed dimension points!

      //3D
      //if it's part of a tet, it's on the boundary if there are any faces without 2 tets
      bool partOfAnyTets = false;
      for(VertexEdgeIterator veit(*this, vh); !veit.done(); veit.advance()) {
         EdgeHandle eh = veit.current();
         for(EdgeFaceIterator efit(*this, eh); !efit.done(); efit.advance()) {
            FaceHandle fh = efit.current();
            int tets = faceIncidentTetCount(fh);
            if(tets > 0)
               partOfAnyTets = true;
            if(tets == 1)
               return true;
         }
      }

      if(partOfAnyTets)
         return false;


      //2D
      //if it's part of a face, it's on the boundary if, the faces don't form a closed loop
      bool partOfAnyFace = false;
      for(VertexEdgeIterator veit(*this, vh); !veit.done(); veit.advance()) {
         EdgeHandle eh = veit.current();
         int faces = edgeIncidentFaceCount(eh);
         if(faces > 0)
            partOfAnyFace = true;
         if(faces == 1)
            return true;
      }
   
      if(partOfAnyFace)
         return false;

      //1D
      int numEdges = vertexIncidentEdgeCount(vh);
      if(numEdges > 0) //it's on the boundary if there's only 1 outgoing edge
         return numEdges == 1;
      else  //otherwise it's an isolated vertex, which we define as not being on the boundary
         return false;
   }

   
   bool SimplicialComplex::isManifold(const FaceHandle& fh) const {
      //in a 3D scenario, it's manifold if it belongs to one or two tets
      //in a 2D scenario, it is guaranteed to be manifold
      int numTets = faceIncidentTetCount(fh);
      return numTets < 3;
   }

   bool SimplicialComplex::isManifold(const EdgeHandle& eh) const {

      //check if we're part of any tetrahedra, looking for stray faces or non-manifold faces
      bool partOfAnyTets = false;
      bool freeFace = false;

      FaceHandle boundaryFace;
      std::set<FaceHandle> faceSet;

      for(EdgeFaceIterator efit(*this, eh); !efit.done(); efit.advance()) {
         FaceHandle fh = efit.current();
         faceSet.insert(fh);

         int tets = faceIncidentTetCount(fh);
         if(tets == 0) {
            freeFace = true;
         }
         else if(tets == 1 || tets == 2) {
           partOfAnyTets = true;
           if(tets == 1) 
              boundaryFace = fh; //save the boundary face for later
         }
         else if (tets > 2) { //an edge is non-manifold if it has a non-manifold face
           return false;
         }
      }

      if(partOfAnyTets) {
         if(freeFace) {
            return false;
         }
         else { 
            //if we have only tets, then the test is whether we can walk around the edge, via tet-face connections
            //and in doing so visit all the faces.  Because of the check above, guaranteed no non-manifold faces,
            //so there are definitely no branches in our walk.
            
            std::set<FaceHandle> unvisitedSet = faceSet; 
           
            //start from one boundary if it exists, otherwise start anywhere, since it should be a closed loop
            FaceHandle startFace = boundaryFace.isValid() ? boundaryFace : *unvisitedSet.begin(); 
            unvisitedSet.erase(startFace);

            //try to walk around the edge, via tets, either to the other boundary or back to the start again, 
            //and see if we visit all the faces. if not, we know there are (at least) two disconnected components.
            FaceHandle prevFace = startFace;
            
            TetHandle prevTet = TetHandle::invalid();
            do {
               
               //Get the next tet connected to this face.
               FaceTetIterator ftit(*this, prevFace);
               while(ftit.current() == prevTet && !ftit.done())
                  ftit.advance();
               if(ftit.done()) break; //ran out of tets, dead end

               //find the next face of the tet that is part of the faceSet (and is not the preceding one)
               TetHandle curTet = ftit.current(); 
               TetFaceIterator tfit(*this, curTet);
               while(!tfit.done() && (tfit.current() == prevFace || faceSet.find(tfit.current()) == faceSet.end()))
                  tfit.advance();

               assert(!tfit.done());
               
               //walk to the next face, and mark this one visited
               prevFace = tfit.current(); 
               prevTet = curTet;
               unvisitedSet.erase(prevFace);

            } while(prevFace != startFace && unvisitedSet.size() > 0);

            //if our cycle around the edges's faces left no face unvisited, then we're manifold
            return unvisitedSet.size() == 0;
         }
      }
      else {
         //dimension 2/1/0
         int numFaces = edgeIncidentFaceCount(eh);
         return numFaces < 3;  
      }
     
   }

   bool SimplicialComplex::isManifold(const VertexHandle& vh) const {

      //dimension 3
      //the surrounding boundary tets form a set that is fully connected via faces
      //and there are no stray faces or edges
      bool partOfAnyTets = false;
      bool freeFace = false;
      bool freeEdge = false;
      std::set<FaceHandle> boundaryFaces;
      for(VertexEdgeIterator veit(*this, vh); !veit.done(); veit.advance()) {
         EdgeHandle eh = veit.current();

         //check for solo edges
         if(edgeIncidentFaceCount(eh) == 0) {
            freeEdge = true;
         }

         for(EdgeFaceIterator efit(*this, eh); !efit.done(); efit.advance()) { //this will visit some faces multiple times, no biggie.
            FaceHandle fh = efit.current();
            int tets = faceIncidentTetCount(fh);
            if(tets > 0) { //we have encountered at least one tet, so we know we're in 3D mode.
               partOfAnyTets = true;
            }
            if(tets == 0) { //check for solo faces
               freeFace = true;
            }
         }
      }
      
      if(partOfAnyTets) {
         if(freeFace || freeEdge) {
            return false;
         }
         else {
            //test if all faces touching the vertex can be reached by starting at one faces and only walking across tets

           std::set<FaceHandle> boundaryFaces;

            //collect all the faces, and boundaryfaces.
            std::set<FaceHandle> faceSet;
            for(VertexEdgeIterator veit(*this, vh); !veit.done(); veit.advance()) {
               EdgeHandle eh = veit.current();
               for(EdgeFaceIterator efit(*this, eh); !efit.done(); efit.advance()) {
                  FaceHandle fh = efit.current();
                  faceSet.insert(fh);
                  if(faceIncidentTetCount(fh) == 1) {
                    boundaryFaces.insert(fh);
                  }
               }
            }

            std::set<FaceHandle> unvisitedFaces = faceSet;
            //now do an exhaustive search across faces to see if we can reach them all.
            std::queue<FaceHandle> facesToVisit;
            facesToVisit.push(*(faceSet.begin())); //push the first tet, whatever it may be
            while(facesToVisit.size() > 0) {
               FaceHandle curFace = facesToVisit.front();
               facesToVisit.pop();
               if(unvisitedFaces.find(curFace) == unvisitedFaces.end()) continue; //already dealt with this one, so skip it
               unvisitedFaces.erase(curFace);

               //look at all neighbouring faces for unvisited ones
               for(FaceTetIterator ftit(*this, curFace); !ftit.done(); ftit.advance()) {
                  TetHandle curTet = ftit.current();
                     for(TetFaceIterator tfit(*this, curTet); !tfit.done(); tfit.advance()) {
                        FaceHandle nbrFace = tfit.current();
                     if(faceSet.find(nbrFace) == faceSet.end() || nbrFace == curFace) continue; //skip irrelevant faces
                     
                     //if we found an unvisited face in the ring, push it
                     if(unvisitedFaces.find(nbrFace) != unvisitedFaces.end())
                        facesToVisit.push(nbrFace);
                  }
               }
            }
            bool allFacesReachable = (unvisitedFaces.size() == 0);

            //Also need to check for the case where there are tets all the way around, EXCEPT
            //for two spots that ony meet at a vertex, (or a single edge?). right?
            
            
            std::set<FaceHandle> unvisitedSet = boundaryFaces;
            if(boundaryFaces.size() > 0) {
               //try to walk around the boundary, via faces/edges, back to the start again. We know it has to be a loop.
               //see if we visit all the edges
               EdgeHandle prevEdge_ = EdgeHandle::invalid();
               FaceHandle startFace = *boundaryFaces.begin();
               FaceHandle prevFace =  startFace;
               do {
                  //get the next edge that is connected to our relevant vertex
                  FaceEdgeIterator feit(*this, prevFace, false);
                  while((feit.current() == prevEdge_ || !isIncident(vh, feit.current())) && !feit.done() )
                     feit.advance();

                  EdgeHandle curEdge = feit.current();

                  //count how many boundary faces this edge is connected to. if more than 2, we're in non-manifold scenario
                  int boundaryFaceCount = 0;
                  for(EdgeFaceIterator efit(*this, curEdge); !efit.done(); efit.advance()) {
                     if(boundaryFaces.find(efit.current()) != boundaryFaces.end())
                        ++boundaryFaceCount;
                  }
                  
                  if(boundaryFaceCount > 2) //this is non-manifold connection between two tets, via an edge
                     return false;                     
                  assert(boundaryFaceCount == 2);

                  //find the next face that is part of the faceSet, but is not the preceding one
                  EdgeFaceIterator efit(*this, curEdge);
                  while(!efit.done() && (efit.current() == prevFace || boundaryFaces.find(efit.current()) == boundaryFaces.end()))
                     efit.advance();

                  assert(!feit.done());
                  
                  //move to next edge
                  prevFace = efit.current(); 
                  prevEdge_ = curEdge;
                  unvisitedSet.erase(prevFace);
               } while(prevFace != startFace && unvisitedSet.size() > 0);
            }
            //if our linear cycle around the vertex's boundary faces left no boundary face unvisited, then we're manifold
            bool boundaryConnected = (unvisitedSet.size() == 0);

            return allFacesReachable && boundaryConnected; //if we successfully walked over all the faces via tets, then it has to be manifold

         }
      }

      //dimension 2
      //need to check if all adjacent faces form a *single* closed loop (or half-loop)
      bool partOfAnyFaces = false;
      std::set<EdgeHandle> edgeSet;
      EdgeHandle boundaryEdge;
      for(VertexEdgeIterator veit(*this, vh); !veit.done(); veit.advance()) {
         EdgeHandle eh = veit.current();
         edgeSet.insert(eh);

         int faces = edgeIncidentFaceCount(eh);
         
         if(faces > 0) //we know we're in 2D
            partOfAnyFaces = true;
         
         if(faces == 1) //we have found one end to the potential series of faces around the vertex
            boundaryEdge = eh;

         if(faces == 3) //if there is a non-manifold edge, then we know the vertex is non-manifold too.
            return false;
         
      }

      if(partOfAnyFaces && freeEdge) { //edge is part of at least one face, so do the dimension 2 case
         return false;
      }
      else if(partOfAnyFaces) { //test if we have a single closed loop
         
         std::set<EdgeHandle> unvisitedSet = edgeSet;
         
         EdgeHandle startEdge = boundaryEdge.isValid() ? boundaryEdge : *unvisitedSet.begin();
         unvisitedSet.erase(startEdge);

         //try to walk around the vertex, via faces, either to a dead end or back to the start again, and see if we visit all the edges
         EdgeHandle prevEdge_ = startEdge;
         FaceHandle prevFace =  FaceHandle::invalid();
         do {
            EdgeFaceIterator efit(*this, prevEdge_);
            while(efit.current() == prevFace && !efit.done())
               efit.advance();
            if(efit.done()) break; //ran out of faces, dead end

            FaceHandle curFace = efit.current(); //find the next edge that is part of the edgeSet, but not the preceding one
            FaceEdgeIterator feit(*this, curFace, false);
            while(!feit.done() && (feit.current() == prevEdge_ || edgeSet.find(feit.current()) == edgeSet.end()))
               feit.advance();
            
            assert(!feit.done());
            prevEdge_ = feit.current(); //move to next edge
            prevFace = curFace;

            unvisitedSet.erase(prevEdge_);
            
         } while(prevEdge_ != startEdge && unvisitedSet.size() > 0);
         
         //if our cycle around the vertex's edges left no edge unvisited, then we're manifold
         return unvisitedSet.size() == 0;
      }
      
      //dimension 1/0
      int numEdges = vertexIncidentEdgeCount(vh);
      return numEdges < 3;
   }




} //namespace SimplexMesh