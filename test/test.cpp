#include "SimplicialComplex.h"

#include <iostream>

using namespace SimplexMesh;

bool test_constructTetAndIterateSimplices();
bool test_constructSimplicesFromVerts();
bool test_edgeDuplication();
bool test_faceDuplication();
bool test_faceCreationValid();
bool test_tetDuplication();
bool test_tetCreationValid();
bool test_vertexVertexIterator();

typedef bool (*test_func)();

const int test_count = 8;
test_func tests[] = {test_constructSimplicesFromVerts,
                     test_constructTetAndIterateSimplices,
                     test_edgeDuplication,
                     test_faceDuplication,
                     test_faceCreationValid,
                     test_tetDuplication,
                     test_tetCreationValid,
                     test_vertexVertexIterator};


void main() {
   int success_count = 0;
   for(int i = 0; i < test_count; ++i) {
      if( (*tests[i])() )
         ++success_count;
   }

   std::cout << "Tests passed: " << success_count << " of " << test_count << std::endl;
   
}

bool test_constructTetAndIterateSimplices() {
   
   //construct one tetrahedron, from the ground up.
   //Iterate through its elements.
   int failures = 0;

   SimplicialComplex mesh;
   
   VertexHandle v0 = mesh.addVertex();
   VertexHandle v1 = mesh.addVertex();
   VertexHandle v2 = mesh.addVertex();
   VertexHandle v3 = mesh.addVertex();   
   
   EdgeHandle e0 = mesh.addEdge(v0, v1);
   EdgeHandle e1 = mesh.addEdge(v0, v2);
   EdgeHandle e2 = mesh.addEdge(v0, v3);
   EdgeHandle e3 = mesh.addEdge(v1, v2);
   EdgeHandle e4 = mesh.addEdge(v1, v3);
   EdgeHandle e5 = mesh.addEdge(v2, v3);
   
   FaceHandle f0 = mesh.addFace(e0, e1, e3);
   FaceHandle f1 = mesh.addFace(e3, e4, e5);
   FaceHandle f2 = mesh.addFace(e0, e2, e4);
   FaceHandle f3 = mesh.addFace(e1, e2, e5);

   TetHandle t0 = mesh.addTet(f0,f1,f2,f3);

   //Assign data.
   VertexProperty<int> vertexID(mesh);
   vertexID[v0] = 0;
   vertexID[v1] = 1;
   vertexID[v2] = 2;
   vertexID[v3] = 3;

   EdgeProperty<int> edgeID(mesh);
   edgeID[e0] = 0;
   edgeID[e1] = 1;
   edgeID[e2] = 2;
   edgeID[e3] = 3;
   edgeID[e4] = 4;
   edgeID[e5] = 5;

   FaceProperty<int> faceID(mesh);
   faceID[f0] = 0;
   faceID[f1] = 1;
   faceID[f2] = 2;
   faceID[f3] = 3;

   TetProperty<int> tetID(mesh);
   tetID[t0] = 0;
   
   int i = 0;
   for(VertexIterator it(mesh); !it.done(); it.advance()) {
      VertexHandle cur = it.current();
      if(vertexID[cur] != i) ++failures;
      ++i;
   }
   
   i = 0;
   for(EdgeIterator it(mesh); !it.done(); it.advance()) {
      EdgeHandle cur = it.current();
      if(edgeID[cur] != i) ++failures;
      ++i;
   }
   
   i = 0;
   for(FaceIterator it(mesh); !it.done(); it.advance()) {
      FaceHandle cur = it.current();
      if(faceID[cur] != i) ++failures;
      ++i;
   }
   
   i = 0;
   for(TetIterator it(mesh); !it.done(); it.advance()) {
      TetHandle cur = it.current();
      if(tetID[cur] != i) ++failures;
      ++i;
   }
   
   return failures == 0;
}

bool test_constructSimplicesFromVerts() {
   
   //Try adding faces directly from vertices
   SimplicialComplex mesh;

   VertexHandle v0 = mesh.addVertex();
   VertexHandle v1 = mesh.addVertex();
   VertexHandle v2 = mesh.addVertex();
   VertexHandle v3 = mesh.addVertex();   

   mesh.addFace(v0, v1, v2);
   mesh.addFace(v0, v1, v3);
   mesh.addFace(v0, v2, v3);
   mesh.addFace(v1, v2, v3);

   //Try adding a tet directly from vertices
   SimplicialComplex mesh2;

   v0 = mesh2.addVertex();
   v1 = mesh2.addVertex();
   v2 = mesh2.addVertex();
   v3 = mesh2.addVertex();   

   TetHandle th = mesh.addTet(v0, v1, v2, v3);

   return true;
}

bool test_edgeDuplication() {
   SimplicialComplex mesh;
   
   mesh.setSafeMode(true);
   
   VertexHandle v0 = mesh.addVertex();
   VertexHandle v1 = mesh.addVertex();
   VertexHandle v2 = mesh.addVertex();
  
   EdgeHandle edge0 = mesh.addEdge(v0, v1);
   EdgeHandle edge1 = mesh.addEdge(v0, v1);
   
   return edge0.isValid() && !edge1.isValid();
}

bool test_faceDuplication() {
   SimplicialComplex mesh;

   mesh.setSafeMode(true);

   VertexHandle v0 = mesh.addVertex();
   VertexHandle v1 = mesh.addVertex();
   VertexHandle v2 = mesh.addVertex();

   FaceHandle face0 = mesh.addFace(v0, v1, v2);
   FaceHandle face1 = mesh.addFace(v0, v1, v2);
   FaceHandle face2 = mesh.addFace(v0, v2, v1);

   return face0.isValid() && !face1.isValid() && !face1.isValid();
}


bool test_faceCreationValid() {
   SimplicialComplex mesh;

   mesh.setSafeMode(true);

   VertexHandle v0 = mesh.addVertex();
   VertexHandle v1 = mesh.addVertex();
   VertexHandle v2 = mesh.addVertex();
   VertexHandle v3 = mesh.addVertex();

   EdgeHandle e0 = mesh.addEdge(v0, v1);
   EdgeHandle e1 = mesh.addEdge(v1, v2);
   EdgeHandle e2 = mesh.addEdge(v0, v2);
   EdgeHandle e3 = mesh.addEdge(v0, v3);

   FaceHandle face0 = mesh.addFace(e0, e1, e2); //good face
   FaceHandle face1 = mesh.addFace(e0, e1, e3); //bad face

   return face0.isValid() && !face1.isValid();
}

bool test_tetDuplication() {
   SimplicialComplex mesh;

   mesh.setSafeMode(true);

   VertexHandle v0 = mesh.addVertex();
   VertexHandle v1 = mesh.addVertex();
   VertexHandle v2 = mesh.addVertex();
   VertexHandle v3 = mesh.addVertex();

   TetHandle tet0 = mesh.addTet(v0, v1, v2, v3);
   TetHandle tet1 = mesh.addTet(v0, v1, v2, v3);
   TetHandle tet2 = mesh.addTet(v0, v1, v3, v2);
   
   return tet0.isValid() && !tet1.isValid() && !tet2.isValid();
}

bool test_tetCreationValid() {
   SimplicialComplex mesh;

   mesh.setSafeMode(true);

   VertexHandle v0 = mesh.addVertex();
   VertexHandle v1 = mesh.addVertex();
   VertexHandle v2 = mesh.addVertex();
   VertexHandle v3 = mesh.addVertex();
   VertexHandle v4 = mesh.addVertex();

   FaceHandle f0 = mesh.addFace(v0,v1,v2);
   FaceHandle f1 = mesh.addFace(v0,v2,v3);
   FaceHandle f2 = mesh.addFace(v0,v1,v3);
   FaceHandle f3 = mesh.addFace(v1,v2,v3);
   FaceHandle f4 = mesh.addFace(v0, v2, v4);

   TetHandle t0 = mesh.addTet(f0,f1,f2,f3);
   TetHandle t1 = mesh.addTet(f0,f1,f2,f4);
   TetHandle t2 = mesh.addTet(f0,f1,f3,f4);

   return t0.isValid() && !t1.isValid() && !t2.isValid();
}

bool test_vertexVertexIterator() {
    SimplicialComplex mesh;

    mesh.setSafeMode(true);

    VertexHandle v0 = mesh.addVertex();
    VertexHandle v1 = mesh.addVertex();
    VertexHandle v2 = mesh.addVertex();
    VertexHandle v3 = mesh.addVertex();
    VertexHandle v4 = mesh.addVertex();

    FaceHandle f0 = mesh.addFace(v0,v1,v2);
    FaceHandle f1 = mesh.addFace(v0,v2,v3);
    FaceHandle f2 = mesh.addFace(v0,v1,v3);
    FaceHandle f3 = mesh.addFace(v1,v2,v3);
    FaceHandle f4 = mesh.addFace(v0,v2,v4);

    int count = 0;
    for(VertexVertexIterator vit(mesh, v0); !vit.done(); vit.advance()) {
        VertexHandle vh = vit.current();
        ++count;
    }
    if(count != 4) 
        return false;

    count = 0;
    for(VertexVertexIterator vit(mesh, v4); !vit.done(); vit.advance()) {
        VertexHandle vh = vit.current();
        ++count;
    }
    if(count != 2) 
        return false;

    count = 0;
    for(VertexVertexIterator vit(mesh, v3); !vit.done(); vit.advance()) {
        VertexHandle vh = vit.current();
        ++count;
    }
    if(count != 3) 
        return false;
    
    return true;
}