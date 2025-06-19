/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2025 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "Matrix.H"
#include "MatrixSpace.H"

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Form, class Type>
void Foam::Matrix<Form, Type>::allocate()
{
    if (mRows_ > 0 && nCols_ > 0)
    {
        v_ = new Type[size()];
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(const label m, const label n)
:
    mRows_(m),
    nCols_(n),
    v_(nullptr)
{
    if (mRows_ < 0 || nCols_ < 0)
    {
        FatalErrorInFunction
            << "Incorrect m, n " << mRows_ << ", " << nCols_
            << abort(FatalError);
    }

    allocate();
}

template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(const label m, const label n, const zero)
:
    mRows_(m),
    nCols_(n),
    v_(nullptr)
{
    if (mRows_ < 0 || nCols_ < 0)
    {
        FatalErrorInFunction
            << "Incorrect m, n " << mRows_ << ", " << nCols_
            << abort(FatalError);
    }

    allocate();

    if (v_)
    {
        std::fill(v_, v_ + size(), Zero);
    }
}

template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(const label m, const label n, const Type& s)
:
    mRows_(m),
    nCols_(n),
    v_(nullptr)
{
    if (mRows_ < 0 || nCols_ < 0)
    {
        FatalErrorInFunction
            << "Incorrect m, n " << mRows_ << ", " << nCols_
            << abort(FatalError);
    }

    allocate();

    if (v_)
    {
        std::fill(v_, v_ + size(), s);
    }
}

template<class Form, class Type>
template<class InputIterator>
Foam::Matrix<Form, Type>::Matrix
(
    const label m,
    const label n,
    InputIterator first,
    InputIterator last
)
:
    mRows_(m),
    nCols_(n),
    v_(nullptr)
{
    if (mRows_ < 0 || nCols_ < 0)
    {
        FatalErrorInFunction
            << "Incorrect m, n " << mRows_ << ", " << nCols_
            << abort(FatalError);
    }

    if (std::distance(first, last) != mRows_ * nCols_)
    {
        FatalErrorInFunction
            << "Number of values provided " << std::distance(first, last)
            << " is not the same as the number of matrix elements "
            << mRows_ * nCols_
            << abort(FatalError);
    }

    allocate();

    if (v_)
    {
        std::copy(first, last, v_);
    }
}

template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix
(
    const label m,
    const label n,
    std::initializer_list<Type> lst
)
:
    Matrix(m, n, lst.begin(), lst.end())
{}

template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix
(
    std::initializer_list<std::initializer_list<Type>> lstLst
)
:
    mRows_(lstLst.size()),
    nCols_(lstLst.begin()->size())
{
    allocate();

    label rowi = 0;
    for (const auto& rowIter : lstLst)
    {
        if (rowIter.size() != nCols_)
        {
            FatalErrorInFunction
                << "Number of columns in row " << rowi
                << " is not equal to " << nCols_
                << abort(FatalError);
        }

        std::copy(rowIter.begin(), rowIter.end(), v_ + rowi * nCols_);
        rowi++;
    }
}

template<class Form, class Type>
Foam::Matrix<Form, Type>::Matrix(const Matrix<Form, Type>& M)
:
    mRows_(M.mRows_),
    nCols_(M.nCols_),
    v_(nullptr)
{
    if (M.v_)
    {
        allocate();
        std::copy(M.v_, M.v_ + size(), v_);
    }
}

template<class Form, class Type>
template<class Form2>
Foam::Matrix<Form, Type>::Matrix(const Matrix<Form2, Type>& M)
:
    mRows_(M.m()),
    nCols_(M.n()),
    v_(nullptr)
{
    if (M.v())
    {
        allocate();
        std::copy(M.v(), M.v() + size(), v_);
    }
}

// * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * * //

template<class Form, class Type>
Foam::Matrix<Form, Type>::~Matrix()
{
    delete[] v_;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Form, class Type>
void Foam::Matrix<Form, Type>::clear()
{
    delete[] v_;
    v_ = nullptr;
    mRows_ = nCols_ = 0;
}

template<class Form, class Type>
void Foam::Matrix<Form, Type>::transfer(Matrix<Form, Type>& M)
{
    clear();

    mRows_ = M.mRows_;
    nCols_ = M.nCols_;
    v_ = M.v_;
    M.v_ = nullptr;
    M.mRows_ = M.nCols_ = 0;
}

template<class Form, class Type>
void Foam::Matrix<Form, Type>::setSize(const label m, const label n)
{
    Matrix<Form, Type> newMatrix(m, n, Zero);
    label minM = std::min(m, mRows_);
    label minN = std::min(n, nCols_);

    for (label i = 0; i < minM; ++i)
    {
        for (label j = 0; j < minN; ++j)
        {
            newMatrix(i, j) = (*this)(i, j);
        }
    }

    transfer(newMatrix);
}

template<class Form, class Type>
Form Foam::Matrix<Form, Type>::T() const
{
    Form At(n(), m());
    for (label i = 0; i < m(); ++i)
    {
        for (label j = 0; j < n(); ++j)
        {
            At(j, i) = (*this)(i, j);
        }
    }
    return At;
}

// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator=(const Matrix<Form, Type>& M)
{
    if (this == &M)
    {
        FatalErrorInFunction
            << "Attempted assignment to self"
            << abort(FatalError);
    }

    if (mRows_ != M.mRows_ || nCols_ != M.nCols_)
    {
        clear();
        mRows_ = M.mRows_;
        nCols_ = M.nCols_;
        allocate();
    }

    std::copy(M.v_, M.v_ + size(), v_);
}

template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator=(const Type& s)
{
    std::fill(v_, v_ + size(), s);
}

template<class Form, class Type>
void Foam::Matrix<Form, Type>::operator=(const zero)
{
    std::fill(v_, v_ + size(), Zero);
}

// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

template<class Form, class Type>
const Type& Foam::max(const Matrix<Form, Type>& M)
{
    const label mn = M.size();
    if (mn)
    {
        const Type* Mv = M.v();
        return *std::max_element(Mv, Mv + mn);
    }
    FatalErrorInFunction
        << "Matrix is empty"
        << abort(FatalError);
}
