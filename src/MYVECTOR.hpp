#ifndef MYVECTOR_HPP
#define MYVECTOR_HPP

class myVector{
    
private:
    int size;
    double *data;
    
public:
    myVector(int size);
    ~myVector();

    int Dimension() const;


    //USer-defined operators
    void Print();
    void Initialize(double a);
    void Initialize(double *v);
};


class myMatrix{
    
private:
    int nrows;
    int ncols;
    double **data;
    
public:
    myMatrix(int nrows, int ncols);
    ~myMatrix();

    int Dimension() const;
    
    //User-defined operators
    void Print();
    void Initialize(double a);
    void Initialize(double **data);
};


#endif
