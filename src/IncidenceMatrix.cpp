#include "IncidenceMatrix.h"

namespace SimplexMesh {

int signum(int val) {
  return (val >= 0 ? 1 : -1);

}

IncidenceMatrix::IncidenceMatrix() : 
   n_rows(0), n_cols(0), m_indices(0)
{
}

IncidenceMatrix::IncidenceMatrix(unsigned int rows, unsigned int cols) : 
   n_rows(rows), n_cols(cols), 
   m_indices(rows, std::vector<int>())
{
}

unsigned int IncidenceMatrix::getNumEntriesInRow(unsigned int row) const {
 assert(row < n_rows);
 return m_indices[row].size();
}

int IncidenceMatrix::getValueByIndex(unsigned int i, unsigned int index_in_row) const {
 assert(i < n_rows);
 assert(index_in_row < m_indices[i].size());

 return signum(m_indices[i][index_in_row]);
}

unsigned int IncidenceMatrix::getColByIndex(unsigned int i, unsigned int index_in_row) const {
  assert(i < n_rows);
  assert(index_in_row < m_indices[i].size());

  return abs(m_indices[i][index_in_row]) - 1;
}

void IncidenceMatrix::cycleRow(unsigned int i) {
  int t = m_indices[i][0];
  int row_len = m_indices[i].size();
  for(int j = 0; j < row_len-1; ++j)
    m_indices[i][j] = m_indices[i][j+1];
  m_indices[i][row_len-1] = t;
}

void IncidenceMatrix::setByIndex(unsigned int i, unsigned int index_in_row, unsigned int col, int value) {
   assert(value == 1 || value == -1);
   if(index_in_row >= m_indices[i].size()) m_indices[i].resize(index_in_row+1);
  
   m_indices[i][index_in_row] = (col+1)*value;
}

void IncidenceMatrix::set(unsigned int i, unsigned int j, int new_val) {
   assert(i < n_rows && j < n_cols);
   if(new_val == 0) {
     remove(i,j);
     return;
   }

   assert(new_val == 1 || new_val == -1);
  
   int colShift = j+1;
   bool found = false;
   for(unsigned int cur = 0; cur < m_indices[i].size(); ++cur) {
     if(abs(m_indices[i][cur]) == (int)colShift ) {
        m_indices[i][cur] = signum(new_val)*colShift;
        found = true;
        break;
     }
   }
   if(!found)
      m_indices[i].push_back(signum(new_val)*colShift);
   
}

int IncidenceMatrix::get(unsigned int i, unsigned int j) const {
   assert(i < n_rows && j < n_cols);
  
   int colShift = j+1;
   for(unsigned int k=0; k<m_indices[i].size(); ++k){
      if(abs(m_indices[i][k])==(int)colShift){
         return signum(m_indices[i][k]);
      }
   }
   return 0;
}

void IncidenceMatrix::remove(unsigned int i, unsigned int j) {
   assert(i<n_rows && j < n_cols);

   int colShift = j+1;
   for(unsigned int k=0; k<m_indices[i].size(); ++k){
      if(abs(m_indices[i][k])==colShift){
         m_indices[i].erase(m_indices[i].begin()+k);
         return;
      }
   }
}

bool IncidenceMatrix::exists(unsigned int i, unsigned int j) const {
   assert(i<n_rows && j < n_cols);
  
   int colShift = j+1;
   for(unsigned int k=0; k<m_indices[i].size(); ++k){
      if(abs(m_indices[i][k])==colShift){
         return true;
      }
   }
   return false;
}

void IncidenceMatrix::addRows(unsigned int rows) {
   n_rows += rows;
   m_indices.resize(n_rows);
}

void IncidenceMatrix::addCols(unsigned int cols) {
   n_cols += cols;
}

void IncidenceMatrix::zeroRow( unsigned int i )
{
   m_indices[i].clear();
}

void IncidenceMatrix::printMatrix() const 
{
   printf("Dimensions (%d,%d):\n", n_rows, n_cols);
   for(unsigned int row = 0; row < m_indices.size(); ++row) {
      printf("%d: ", row);
      for(unsigned int i = 0; i < m_indices[row].size(); ++i)
        printf("%c%d ", m_indices[row][i] > 0?'+':'-', abs(m_indices[row][i])-1);
      printf("\n");
   }
}

void IncidenceMatrix::zeroAll()
{
  for(unsigned int i = 0; i < getNumRows(); ++i)
    zeroRow(i);
}

}