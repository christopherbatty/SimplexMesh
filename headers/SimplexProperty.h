#ifndef SIMPLEXPROPERTY_H
#define SIMPLEXPROPERTY_H

#include "SimplicialComplex.h"

namespace SimplexMesh {

//Base class, so TopologicalObject can store pointers to TopObjProperties of different types in a single list.
class SimplexPropertyBase {

public:

  SimplexPropertyBase(SimplicialComplex* obj) : m_obj(obj) {}
  virtual ~SimplexPropertyBase() {}

protected:

  //Only the TopologicalObject that owns this object is allowed to manipulate its size, 
  //in order to match the number of the given simplex type. (eg. #of edge properties should equal #of edges)
  virtual size_t size() const = 0;
  virtual void resize(size_t n) = 0;

  //The simplex mesh this property is associated with.
  SimplicialComplex* m_obj;

  friend class SimplicialComplex;
};


//The templated object property that stores the data. This is a base class that should not be used, since it doesn't get registered with
//a particular simplex type (edge, face, etc.), and cannot be resized.
template <class T>
class SimplexProperty : public SimplexPropertyBase {

public:
  SimplexProperty(SimplicialComplex* obj, size_t n) : SimplexPropertyBase(obj), m_data(n) {}

  virtual ~SimplexProperty() {}

  void assign(const T& data_value) { for(unsigned int i = 0; i < m_data.size(); ++i) m_data[i] = data_value; }
  
protected: 
  
  size_t size() const { return m_data.size(); }
  void resize(size_t n) { m_data.resize(n); }

  std::vector<T> m_data;
  
};

//The properties associated to simplices of different dimensions. They are instantiated with a pointer to the the TopologicalObject that
//"owns" them, and registered with that object, so they can be automatically resized to match the mesh when necessary.
template <class T>
class VertexProperty : public SimplexProperty<T> {

public:

  VertexProperty(SimplicialComplex* obj) : SimplexProperty<T>(obj,obj->numVertexSlots()) {
    this->m_obj->registerVertexProperty(this);
  }

  explicit VertexProperty(const VertexProperty& prop) : SimplexProperty<T>(prop.m_obj,prop.m_obj->numVertexSlots()) {
    this->m_obj->registerVertexProperty(this);
    this->m_data = prop.m_data;
  }

  ~VertexProperty() {
    this->m_obj->removeVertexProperty(this);
  }

  VertexProperty& operator=(const VertexProperty& other) {

    if(this != &other) {
        this->m_obj->removeVertexProperty(this);

        this->m_obj = other.m_obj;
        this->m_data = other.m_data;

        this->m_obj->registerVertexProperty(this);
    }

    return *this;
  }

  T& operator[] (const VertexHandle& h) { 
    assert(h.idx() >= 0 && h.idx() < (int)m_data.size());
    return m_data[h.idx()]; 
  }

  T const& operator[] (const VertexHandle& h) const { 
    assert(h.idx() >= 0 && h.idx() < (int)m_data.size()); 
    return m_data[h.idx()]; 
  }

};


template <class T>
class EdgeProperty : public SimplexProperty<T> {
public:
  
  EdgeProperty(SimplicialComplex* obj) : SimplexProperty<T>(obj,obj->numEdgeSlots()) {
    this->m_obj->registerEdgeProperty(this);
  }

  explicit EdgeProperty(const EdgeProperty& prop) : SimplexProperty<T>(prop.m_obj,prop.m_obj->numEdgeSlots()) {
    this->m_obj->registerEdgeProperty(this);
    this->m_data = prop.m_data;
  }

  ~EdgeProperty() {
    this->m_obj->removeEdgeProperty(this);
  }

  EdgeProperty& operator=(const EdgeProperty& other) {

    if(this != &other) {
      this->m_obj->removeEdgeProperty(this);

      this->m_obj = other.m_obj;
      this->m_data = other.m_data;

      this->m_obj->registerEdgeProperty(this);
    }

    return *this;
  }


  T& operator[] (const EdgeHandle& h) { 
    assert(h.idx() >= 0 && h.idx() < (int)m_data.size());
    return m_data[h.idx()]; 
  }

  T const& operator[] (const EdgeHandle& h) const { 
    assert(h.idx() >= 0 && h.idx() < (int)m_data.size()); 
    return m_data[h.idx()]; 
  }
};

template <class T>
class FaceProperty : public SimplexProperty<T> {
public:
  FaceProperty(SimplicialComplex* obj) : SimplexProperty<T>(obj,obj->numFaceSlots()) {
    this->m_obj->registerFaceProperty(this);
  }

  explicit FaceProperty(const FaceProperty& prop) : SimplexProperty<T>(prop.m_obj,prop.m_obj->numFaceSlots()) {
    this->m_obj->registerFaceProperty(this);
    this->m_data = prop.m_data;
  }

  ~FaceProperty() {
    this->m_obj->removeFaceProperty(this);
  }

  FaceProperty& operator=(const FaceProperty& other) {

    if(this != &other) {
      this->m_obj->removeFaceProperty(this);

      this->m_obj = other.m_obj;
      this->m_data = other.m_data;

      this->m_obj->registerFaceProperty(this);
    }

    return *this;
  }

  T& operator[] (const FaceHandle& h) { 
    assert(h.idx() >= 0 && h.idx() < (int)m_data.size());
    return m_data[h.idx()]; 
  }

  T const& operator[] (const FaceHandle& h) const { 
    assert(h.idx() >= 0 && h.idx() < (int)m_data.size()); 
    return m_data[h.idx()]; 
  }

};


template <class T>
class TetProperty : public SimplexProperty<T> {
public:
  TetProperty(SimplicialComplex* obj) : SimplexProperty<T>(obj, obj->numTetSlots()) {
    this->m_obj->registerTetProperty(this);
  }

  explicit TetProperty(const TetProperty& prop) : SimplexProperty<T>(prop.m_obj,prop.m_obj->numTetSlots()) {
    this->m_obj->registerTetProperty(this);
    this->m_data = prop.m_data;
  }

  ~TetProperty() {
    this->m_obj->removeTetProperty(this);
  }

  TetProperty& operator=(const TetProperty& other) {

    if(this != &other) {
      this->m_obj->removeTetProperty(this);

      this->m_obj = other.m_obj;
      this->m_data = other.m_data;

      this->m_obj->registerTetProperty(this);
    }

    return *this;
  }

  T& operator[] (const TetHandle& h) { 
    assert(h.idx() >= 0 && h.idx() < (int)m_data.size());
    return m_data[h.idx()]; 
  }

  T const& operator[] (const TetHandle& h) const { 
    assert(h.idx() >= 0 && h.idx() < (int)m_data.size()); 
    return m_data[h.idx()]; 
  }
  
};

}

#endif