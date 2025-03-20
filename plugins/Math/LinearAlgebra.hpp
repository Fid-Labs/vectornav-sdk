// The MIT License (MIT)
// 
//  VectorNav Software Development Kit (v0.14.2)
// Copyright (c) 2024 VectorNav Technologies, LLC
// 
// 
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
// 
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
// THE SOFTWARE.

#ifndef LINEARALGEBRA_HPP
#define LINEARALGEBRA_HPP

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <initializer_list>
#include <iostream>
#include <cmath>
#include <array>
#include <limits>
#include "TemplateLibrary/Matrix.hpp"

#include "Debug.hpp"

namespace VN
{

namespace LinAlg
{

// Vector math
template <typename T, typename S>
T dot(const Matrix<3, 1, T>& lhs, const Matrix<3, 1, S>& rhs)
{
    return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2];
}

template <typename T, typename S>
T dot(const Matrix<4, 1, T>& lhs, const Matrix<4, 1, S>& rhs)
{
    return lhs[0] * rhs[0] + lhs[1] * rhs[1] + lhs[2] * rhs[2] + lhs[3] * rhs[3];
}

template <size_t m, typename T, typename S>
T dot(const Matrix<m, 1, T>& lhs, const Matrix<m, 1, S>& rhs)
{
    T sum = 0;
    for (size_t i = 0; i < m; i++) { sum += lhs[i] * rhs[i]; }
    return sum;
}

template <typename T>
Matrix<3, 1, T> cross(const Matrix<3, 1, T>& lhs, const Matrix<3, 1, T>& rhs)
{
    Matrix<3, 1, T> retMatrix;
    retMatrix[0] = lhs(1) * rhs(2) - lhs(2) * rhs(1);
    retMatrix[1] = lhs(2) * rhs(0) - lhs(0) * rhs(2);
    retMatrix[2] = lhs(0) * rhs(1) - lhs(1) * rhs(0);
    return retMatrix;
}

// If non-invertable it will return an identity matrix
template <size_t n, typename T>
Matrix<n, n, T> inverse(const Matrix<n, n, T>& mat)
{
    Matrix<n, n, T> ac = mat;
    Matrix<n, n, T> inverse;

    size_t i, j, iPass, imx, icol, irow;
    T det = 1;
    T temp, pivot, factor;

    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++) { inverse(i, j) = 0; }
        inverse(i, i) = 1;
    }

    for (iPass = 0; iPass < n; iPass++)
    {
        imx = iPass;
        for (irow = iPass; irow < n; irow++)
        {
            if (fabs(ac(irow, iPass)) > fabs(ac(imx, iPass))) { imx = irow; }
        }

        if (imx != iPass)
        {
            for (icol = 0; icol < n; icol++)
            {
                temp = inverse(iPass, icol);
                inverse(iPass, icol) = inverse(imx, icol);
                inverse(imx, icol) = temp;
                if (icol >= iPass)
                {
                    temp = ac(iPass, icol);
                    ac(iPass, icol) = ac(imx, icol);
                    ac(imx, icol) = temp;
                }
            }
        }

        pivot = ac(iPass, iPass);
        det = det * pivot;

        if (std::fabs(det) < std::numeric_limits<T>::epsilon())
        {
            return Matrix<n, n, T>::identity();  // Matrix is not invertible
        }

        for (icol = 0; icol < n; icol++)
        {
            inverse(iPass, icol) = inverse(iPass, icol) / pivot;
            if (icol >= iPass) { ac(iPass, icol) = ac(iPass, icol) / pivot; }
        }

        for (irow = 0; irow < n; irow++)
        {
            if (irow != iPass)
            {
                factor = ac(irow, iPass);
                for (icol = 0; icol < n; icol++)
                {
                    if (irow != iPass)
                    {
                        inverse(irow, icol) -= factor * inverse(iPass, icol);
                        ac(irow, icol) -= factor * ac(iPass, icol);
                    }
                }
            }
        }
    }

    return inverse;  // Matrix inversion successful
}

template <typename T>
Matrix<2, 2, T> inverse(const Matrix<2, 2, T>& mat)
{
    Matrix<2, 2, T> nm({mat[3], -mat[1], -mat[2], mat[0]});

    T det = mat[0] * mat[3] - mat[1] * mat[2];

    if (std::fabs(det) < 10 * std::numeric_limits<T>::epsilon()) { return Matrix<2, 2, T>::identity(); }

    T invDet = 1 / det;
    return nm * invDet;
}

template <typename T>
Matrix<3, 3, T> inverse(const Matrix<3, 3, T>& mat)
{
    Matrix<3, 3, T> nm;
    T det = mat(0, 0) * (mat(1, 1) * mat(2, 2) - mat(2, 1) * mat(1, 2)) - mat(0, 1) * (mat(1, 0) * mat(2, 2) - mat(1, 2) * mat(2, 0)) +
            mat(0, 2) * (mat(1, 0) * mat(2, 1) - mat(1, 1) * mat(2, 0));

    if (std::fabs(det) < std::numeric_limits<T>::epsilon()) { return Matrix<3, 3, T>::identity(); }

    T invDet = 1 / det;

    nm(0, 0) = (mat(1, 1) * mat(2, 2) - mat(2, 1) * mat(1, 2)) * invDet;
    nm(0, 1) = (mat(0, 2) * mat(2, 1) - mat(0, 1) * mat(2, 2)) * invDet;
    nm(0, 2) = (mat(0, 1) * mat(1, 2) - mat(0, 2) * mat(1, 1)) * invDet;
    nm(1, 0) = (mat(1, 2) * mat(2, 0) - mat(1, 0) * mat(2, 2)) * invDet;
    nm(1, 1) = (mat(0, 0) * mat(2, 2) - mat(0, 2) * mat(2, 0)) * invDet;
    nm(1, 2) = (mat(1, 0) * mat(0, 2) - mat(0, 0) * mat(1, 2)) * invDet;
    nm(2, 0) = (mat(1, 0) * mat(2, 1) - mat(2, 0) * mat(1, 1)) * invDet;
    nm(2, 1) = (mat(2, 0) * mat(0, 1) - mat(0, 0) * mat(2, 1)) * invDet;
    nm(2, 2) = (mat(0, 0) * mat(1, 1) - mat(1, 0) * mat(0, 1)) * invDet;

    return nm;
}

template <size_t m, size_t n, typename T>
Matrix<n, m, T> transpose(const Matrix<m, n, T>& mat)
{
    Matrix<n, m, T> nm;

    for (size_t row = 0; row < m; row++)
    {
        for (size_t col = 0; col < n; col++) { nm[row * n + col] = mat[col * m + row]; }
    }

    return nm;
}

template <size_t m, typename T>
T norm(const Matrix<m, 1, T>& mat)
{
    return std::sqrt(dot(mat, mat));
}

template <size_t m, typename T>
Matrix<m, 1, T> normalize(Matrix<m, 1, T>& mat)
{
    return mat / norm(mat);
}

}  // namespace LinAlg
}  // namespace VN

#endif  // LINEARALGEBRA_HPP
