//
// Created by jonat on 2023-04-20.
//

#ifndef COMP559FINALPROJ_ARRAY3D_H
#define COMP559FINALPROJ_ARRAY3D_H

#include <vector>
#include <iostream>

template <class T, size_t X, size_t Y, size_t Z>
class Array3D
{
public:
    const int size_X = X;
    const int size_Y = Y;
    const int size_Z = Z;
    typedef typename T type;

    Array3D() : buffer(size_X*size_Y*size_Z)
    {
    }

    inline type& at(unsigned int x, unsigned int y, unsigned int z)
    {
        return buffer[x*size_YZ + y*size_Z + z];
    }

    inline const type& at(unsigned int x, unsigned int y, unsigned int z) const
    {
        return buffer[x*size_YZ + y*size_Z + z];
    }

    inline void fill(type val) { std::fill(buffer.begin(), buffer.end(), val); }
    inline type* data() { return &buffer[0]; }

    std::vector<T> buffer;
private:
    const int size_YZ = Y * Z;
};

#endif //COMP559FINALPROJ_ARRAY3D_H
