#! /usr/bin/env python

import itertools

'''
An extremely inefficient and not well guarded matrix class that enables us to
avoid numpy
'''

class Matrix(object):
    '''
    A very simple matrix implementation that supports basic matrix algebra
    '''

    def __init__(self,m):
        '''
        Initialize a matrix with list `m`. If `m` is a list of floats, make the
        matrix a single column matrix with the contents of `m`. Otherwise, `m`
        must be a list of lists of floats representing a matrix.
        '''
        if type(m) is not list:
            raise ValueError, 'Matrices must be intialized with a list or list of lists'
        if type(m[0]) is not list:
            self.m = [[v] for v in m]
        else:
            self.m = m

    @classmethod
    def diag(cls,m):
        '''
        Return a symmetric diagonal matrix with the contents of `m`.
        '''
        return Matrix(
            [ [0.]*i + [m[i]] + [0.]*(len(m)-i-1) for i in range(len(m)) ])

    def T(self):
        '''
        Return a new matrix which is the transposition of this matrix
        '''
        t = [list() for _ in range(len(self.m[0]))]
        for i in range(len(self.m)):
            for j in range(len(self.m[0])):
                t[j].append(self.m[i][j])
        return Matrix(t)

    def _minor(self,i,j):
        return [row[:j] + row[j+1:] for row in (self.m[:i]+self.m[i+1:])]

    def minor(self,i,j):
        '''
        Return a matrix which does not contain the row `i` and the column `j`.
        '''
        return Matrix(self._minor(i,j))

    def _determinant_helper(self,m):
        #base case for 2x2 matrix
        if len(m) == 2:
            return m[0][0]*m[1][1]-m[0][1]*m[1][0]

        det = 0
        for c in range(len(m)):
            det += ((-1)**c)*m[0][c]*self._determinant_helper(self._minor(0,c))
        return det
    def determinant(self):
        '''
        Compute the determinant of the matrix.
        '''
        return self._determinant_helper(self.m)

    def inverse(self):
        '''
        Return the inverse of the matrix. Matrix must be invertable.
        '''
        det = self.determinant()
        if det == 0.0:
            raise ValueError, 'Matrix is not invertable.'
        #special case for 2x2 matrix:
        if len(self.m) == 2:
            return Matrix([[self.m[1][1]/det, -1*self.m[0][1]/det],
                    [-1*self.m[1][0]/det, self.m[0][0]/det]])

        #find matrix of cofactors
        cofactors = []
        for r in range(len(self.m)):
            cofactorRow = []
            for c in range(len(self.m)):
                m_minor = self._minor(r,c)
                cofactorRow.append(((-1)**(r+c)) * self._determinant_helper(m_minor))
            cofactors.append(cofactorRow)
        cofactors = transpose(cofactors)
        for r in range(len(cofactors)):
            for c in range(len(cofactors)):
                cofactors[r][c] = cofactors[r][c]/det
        return Matrix(cofactors)

    def dot(self,m):
        '''
        Return a new matrix which is the product of this matrix and `m`
        '''
        zip_b = zip(*m.m)
        return Matrix([[sum(ele_a*ele_b for ele_a, ele_b in zip(row_a, col_b)) 
                 for col_b in zip_b] for row_a in self.m])

    def __str__(self):
        return '[' + '\n'.join(str(r)+ ',' for r in self.m) + ']'

