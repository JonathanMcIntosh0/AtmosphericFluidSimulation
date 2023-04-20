//
// Created by jonat on 2023-04-20.
//

#ifndef COMP559FINALPROJ_ARRAY3D_H
#define COMP559FINALPROJ_ARRAY3D_H

#include <vector>

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

private:
    const int size_YZ = Y * Z;
    std::vector<T> buffer;
};

#endif //COMP559FINALPROJ_ARRAY3D_H
