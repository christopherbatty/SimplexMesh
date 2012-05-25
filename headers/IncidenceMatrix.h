#ifndef INCIDENCEMATRIX_H
#define INCIDENCEMATRIX_H

#include <cassert>
#include <vector>

namespace SimplexMesh {

  //A simple std::vector-based sparse compressed row incidence matrix 
  //to store the topology of our simplex mesh structure.
  //It needs to be resize-able in order to add/delete simplices.
  class IncidenceMatrix {

  private:

    //Matrix dimensions
    unsigned int n_rows, n_cols;

    //For each row, a list of all column indices, in oriented order where appropriate. 
    //The sign indicates whether the value is intended to be: +1 or -1.
    //NOTE: We shift all column indices up by 1, so the zero'th column is enabled to have a sign!
    std::vector< std::vector<int> > m_indices; 

  public:
    IncidenceMatrix();
    IncidenceMatrix(unsigned int rows, unsigned int cols);

    //Matrix dimensions
    unsigned int getNumRows() const {return n_rows;}
    unsigned int getNumCols() const {return n_cols;}
    void addRows(unsigned int rows);
    void addCols(unsigned int cols);

    //Regular accessors
    void set(unsigned int i, unsigned int j, int new_val);
    int  get(unsigned int i, unsigned int j) const;
    bool exists(unsigned int i, unsigned int j) const;
    void remove(unsigned int i, unsigned int j);
    void zeroRow(unsigned int i);
    void zeroAll();

    void cycleRow(unsigned int index); //permute the row by shifting them all over by 1

    //Constant-time access within the row, by row-index rather than column number
    unsigned int getNumEntriesInRow(unsigned int row) const;
    unsigned int getColByIndex(unsigned int i, unsigned int index_in_row) const;
    int getValueByIndex(unsigned int i, unsigned int index_in_row) const;
    void setByIndex(unsigned int i, unsigned int index_in_row, unsigned int col, int val);

    //Debugging
    void printMatrix() const;
  };

} // namespace SimplexMesh


#endif //INCIDENCEMATRIX_H