#ifndef SIMPLICIALCOMPLEX_H
#define SIMPLICIALCOMPLEX_H

#include <algorithm>

#include "SimplexHandles.h"
#include "IncidenceMatrix.h"

namespace SimplexMesh {

   class SimplexPropertyBase;

   //Determine what level of duplication is allowed
   enum DuplicateSimplexMode {
      Arbitrary, //Any kind of duplicate allowed
      Relaxed,   //Allow simplices that differ in orientation or which share some sub-simplices (e.g. two faces that have 2 of the same 3 edges)
      None       //No duplication at all
   };

   // An object that represents a collection of vertices, edges, faces and tets
   // with associated connectivity information.
   class SimplicialComplex
   {
   public:

      SimplicialComplex();

      //possibly do some additional safety checks when adding simplices
      void setSafeMode(bool safe) { m_safetyChecks = safe; }
      void setDuplicateMode(DuplicateSimplexMode duplication) { m_allowDuplicates = duplication; }

      int numVerts() const;
      int numEdges() const;
      int numFaces() const;
      int numTets() const;

      //Addition: note that the resulting orientation is dependent on the order of parameters
      VertexHandle addVertex();
      EdgeHandle addEdge(const VertexHandle& v0, const VertexHandle& v1); //ordered from v0 to v1
      FaceHandle addFace(const EdgeHandle& e0, const EdgeHandle& e1, const EdgeHandle& e2); //in order of given edges
      TetHandle addTet(const FaceHandle& f0, const FaceHandle& f1, //choose in/out orientation based on face0, possibly flipped
                       const FaceHandle& f2, const FaceHandle& f3,
                       bool flip_face0 = false);

      //auxiliary addition functions for convenience. These will be slower.
      FaceHandle addFace(const VertexHandle& v0, const VertexHandle& v1, const VertexHandle& v2);
      TetHandle addTet(const VertexHandle& v0, const VertexHandle& v1, const VertexHandle& v2, const VertexHandle& v3);

      //deletion: recurse=true will recursively delete its composing sub-simplices if they are not used in any other simplices
      bool deleteVertex(const VertexHandle& vertex);
      bool deleteEdge(const EdgeHandle& edge, bool recurse);
      bool deleteFace(const FaceHandle& face, bool recurse);
      bool deleteTet(const TetHandle& tet, bool recurse);

      //existence: check if simplices still exist
      bool vertexExists(const VertexHandle& vertex) const;
      bool edgeExists(const EdgeHandle& edge) const;
      bool faceExists(const FaceHandle& face) const;
      bool tetExists(const TetHandle& tet) const;
    
      //Exploit fixed ordering to provide simpler access to sub-elements of a simplex
      //These are quick
      VertexHandle getVertex(const EdgeHandle& eh, int index) const;
      EdgeHandle getEdge(const FaceHandle& fh, int index) const;
      FaceHandle getFace(const TetHandle& th, int index) const; //Is there an important/inherent ordering here? should we uniquify by sorting

      //get functions - grab a simplex by its constitutive simplices - careful, slow!
      EdgeHandle getEdge(const VertexHandle& v0, const VertexHandle& v1) const;
      FaceHandle getFace(const EdgeHandle& e0, const EdgeHandle& e1, const EdgeHandle& e2) const;
      TetHandle getTet(const FaceHandle& f0, const FaceHandle& f1, const FaceHandle& f2, const FaceHandle& f3) const;

      //determine orientation (+1 or -1)
      int getRelativeOrientation(const TetHandle& th, const FaceHandle& fh) const;
      int getRelativeOrientation(const FaceHandle& fh, const EdgeHandle& eh) const;
      int getRelativeOrientation(const EdgeHandle& eh, const VertexHandle& vh) const;

      //incidence counts for simplices of neighbouring dimensions
      int vertexIncidentEdgeCount(const VertexHandle& v) const { return m_VE.getNumEntriesInRow(v.idx()); }
      int edgeIncidentFaceCount(const EdgeHandle& e) const { return m_EF.getNumEntriesInRow(e.idx()); }
      int faceIncidentTetCount(const FaceHandle& f) const { return m_FT.getNumEntriesInRow(f.idx()); }

      //manifoldness tests - complicated...
      bool isManifold(const VertexHandle& v) const;
      bool isManifold(const EdgeHandle& e) const;
      bool isManifold(const FaceHandle& f) const;

      //boundary tests
      bool isOnBoundary(const VertexHandle& v) const;
      bool isOnBoundary(const EdgeHandle& e) const;
      bool isOnBoundary(const FaceHandle& f) const;

      //incidence tests
      bool isIncident(const VertexHandle& vh, const EdgeHandle& eh) const;
      bool isIncident(const EdgeHandle& eh, const FaceHandle& fh) const;
      bool isIncident(const FaceHandle& fh, const TetHandle& th) const;


      //Pairwise relationships / traversal
      //----------------------------------
      //Edge/Vert
      VertexHandle fromVertex(const EdgeHandle& eh) const;
      VertexHandle toVertex(const EdgeHandle& eh) const;

      //Face/Edge - assumes consistent orientation and 2D
      FaceHandle frontFace(const EdgeHandle& fh) const;
      FaceHandle backFace(const EdgeHandle& fh) const;

      //Tet/Face - assumes consistent orientation
      TetHandle frontTet(const FaceHandle& fh) const;
      TetHandle backTet(const FaceHandle& fh) const;

      //Simple global traversal functions (terser than iterators, safer/slower when editing meshes)
      //------------------------------------
      //Verts
      VertexHandle nextVertex(const VertexHandle& curVertex) const;
      VertexHandle prevVertex(const VertexHandle& curVertex) const;
      
      //Edges
      EdgeHandle nextEdge(const EdgeHandle& curEdge) const; 
      EdgeHandle prevEdge(const EdgeHandle& curEdge) const;
      
      //Faces
      FaceHandle nextFace(const FaceHandle& curFace) const;
      FaceHandle prevFace(const FaceHandle& curFace) const;

      //Tets
      TetHandle nextTet(const TetHandle& curTet) const;
      TetHandle prevTet(const TetHandle& curTet) const;

      //Simple local traversal functions (stateless, terser than iterators, safer/slower when editing meshes)
      //---------------------------------
      //Vertex traversal
      VertexHandle nextVertex(const EdgeHandle&edge, const VertexHandle& curVertex) const;
      VertexHandle prevVertex(const EdgeHandle&edge, const VertexHandle& curVertex) const;
      VertexHandle nextVertex(const FaceHandle&face, const VertexHandle& curVertex) const;
      VertexHandle prevVertex(const FaceHandle&face, const VertexHandle& curVertex) const;
      
      //Edge traversal
      EdgeHandle nextEdge(const VertexHandle& vert, const EdgeHandle& curEdge) const;
      EdgeHandle prevEdge(const VertexHandle& vert, const EdgeHandle& curEdge) const;
      EdgeHandle nextEdge(const FaceHandle& face, const EdgeHandle& curEdge) const;
      EdgeHandle prevEdge(const FaceHandle& face, const EdgeHandle& curEdge) const;
      
      //Face traversal
      FaceHandle nextFace(const TetHandle& tet, const FaceHandle& curFace) const;
      FaceHandle prevFace(const TetHandle& tet, const FaceHandle& curFace) const;
      FaceHandle nextFace(const EdgeHandle& edge, const FaceHandle& curFace) const;
      FaceHandle prevFace(const EdgeHandle& edge, const FaceHandle& curFace) const;

      //Tet traversal
      TetHandle nextTet(const FaceHandle& face, const TetHandle& curTet) const;
      TetHandle prevTet(const FaceHandle& face, const TetHandle& curTet) const;

      //--------------------------------
   private:

      VertexHandle collapseEdge(const EdgeHandle& eh, const VertexHandle& vertToRemove);

      //Friendship relations
      //////////////////////////////////////////////////////////////////////////

      //Basic iterators
      friend class VertexIterator; friend class EdgeIterator; 
      friend class FaceIterator; friend class TetIterator;

      //Neighbour iterators
      friend class VertexEdgeIterator; friend class EdgeVertexIterator;
      friend class EdgeFaceIterator; friend class FaceEdgeIterator;
      friend class FaceTetIterator; friend class TetFaceIterator;

      //allow properties to access object internals (for registering/unregistering themselves) 
      template<class T> friend class VertexProperty;
      template<class T> friend class EdgeProperty;
      template<class T> friend class FaceProperty;
      template<class T> friend class TetProperty;

      //Internal functions
      //////////////////////////////////////////////////////////////////////////

      EdgeHandle getSharedEdge(const FaceHandle& f0, const FaceHandle& f1) const;
      FaceHandle getSharedFace(const TetHandle &t0, const TetHandle& t1) const;

      //Functions for registering/unregistering properties associated to simplex elements
      void registerVertexProperty(SimplexPropertyBase* prop);
      void removeVertexProperty(SimplexPropertyBase* prop);

      void registerEdgeProperty(SimplexPropertyBase* prop);
      void removeEdgeProperty(SimplexPropertyBase* prop);

      void registerFaceProperty(SimplexPropertyBase* prop);
      void removeFaceProperty(SimplexPropertyBase* prop);

      void registerTetProperty(SimplexPropertyBase* prop);
      void removeTetProperty(SimplexPropertyBase* prop);

      //The number of spaces currently allocated for each simplex type. (Note this is different
      //from the number of active simplices of each type.)
      unsigned int numVertexSlots() const {return m_V.size();}
      unsigned int numEdgeSlots() const {return m_EV.getNumRows();}
      unsigned int numFaceSlots() const {return m_FE.getNumRows();}
      unsigned int numTetSlots() const {return m_TF.getNumRows();}

      //Core Data
      //////////////////////////////////////////////////////////////////////////

      //Simplex counts
      int m_nVerts, m_nEdges, m_nFaces, m_nTets;

      //Fundamental mesh data (incidence matrix format)
      IncidenceMatrix m_TF;  ///< tet-to-face relations
      IncidenceMatrix m_FE;  ///< face-to-edge relations
      IncidenceMatrix m_EV;  ///< edge-to-vert relations
      std::vector<bool> m_V; ///< vertex existence, to support isolated vertices

      //Transposes, needed for efficient deletion/traversal/etc
      IncidenceMatrix m_FT; ///< face-to-tet relations
      IncidenceMatrix m_EF; ///< edge-to-face relations
      IncidenceMatrix m_VE; ///< vert-to-edge relations

      //Pools of empty rows/columns in the above matrices, to efficiently add
      //data in previously deleted slots.
      std::vector<unsigned int> m_deadVerts, m_deadEdges, m_deadFaces, m_deadTets;

      //Lists of simplex properties. These are just pointers so that memory can be managed by the SimplexMesh for adding/deleting.
      //But the data "lives" wherever it has been created.
      std::vector<SimplexPropertyBase*> m_vertProperties;
      std::vector<SimplexPropertyBase*> m_edgeProperties;
      std::vector<SimplexPropertyBase*> m_faceProperties;
      std::vector<SimplexPropertyBase*> m_tetProperties;

      //Option flags
      bool m_safetyChecks;
      DuplicateSimplexMode m_allowDuplicates;

   };



} // namespace SimplexMesh

#include "SimplexProperty.h"
#include "SimplexIterators.h"
#include "SimplexUtil.h"

#endif // SIMPLICIALCOMPLEX_H
