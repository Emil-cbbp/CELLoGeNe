""" Script used to generate M and inversed M for given order and saves it to a file. """

import os
import numpy as np

# Provide a PATH to where the matrices should be stored, e.g.:
PATH = os.getcwd() +'/Matrices/'
existingMatrices = sorted(os.listdir(PATH))

iMatrices = []
matrices = []

def generateInvMatrix(order):
	if(order < len(iMatrices)):
		return
	if(order == 0):
		iMatrices.append([[1]])
		return
	generateInvMatrix(order - 1)
	
	size = 1 << order
	existingMatrices = sorted(os.listdir(PATH))
	if('InverseMatrix_' + str(order) + '.matrix' not in existingMatrices):
		if('InverseMatrix_' + str(order - 1) + '.matrix' in existingMatrices):
			with open(PATH + 'InverseMatrix_' + str(order - 1) + '.matrix') as M:
				prev = np.loadtxt(M)
		else:
			prev = iMatrices[order - 1]
		
		matrix = [[0 for x in range(size)] for y in range(size)]
		m = 0
		size /= 2
		size = int(size)
		for i, ii in zip(range(size), range(size, 2 * size)):
			for j, jj in zip(range(i + 1), range(size, 2 * size + 1)):
				m = prev[i][j]
				matrix[i ][j ] =  m
				matrix[ii][j ] = -m
				matrix[ii][jj] =  m
		np.savetxt(PATH + 'InverseMatrix_' + str(order) + '.matrix', np.array(matrix).astype(np.byte), fmt='%d', delimiter='\t')
		existingMatrices = sorted(os.listdir(PATH))
		iMatrices.append(matrix)
	
	return iMatrices[-1]
	
	
def generateMatrix(order):
	if(order < len(matrices)):
		return
	if(order == 0):
		matrices.append([[1]])
		return
	
	print('Order',order)
	generateMatrix(order - 1)
	
	size = 1 << order
	print('size', int(size), 'order', order)
	existingMatrices = sorted(os.listdir(PATH))
	if('Matrix_' + str(order) + '.matrix' not in existingMatrices):
		if('Matrix_' + str(order - 1) + '.matrix' in existingMatrices):
			with open(PATH + 'Matrix_' + str(order - 1) + '.matrix') as M:
				prev = np.loadtxt(M)
		else:
			prev = matrices[order - 1]
		matrix = [[0 for x in range(size)] for y in range(size)]
		m = 0
		
		size /= 2
		size = int(size)
		for i, ii in zip(range(size), range(size, 2 * size)):
			for j, jj in zip(range(i + 1), range(size, 2 * size + 1)):
				m = prev[i][j]
				matrix[i ][j ] = m
				matrix[ii][j ] = m
				matrix[ii][jj] = m
		np.savetxt(PATH + 'Matrix_' + str(order) + '.matrix', np.array(matrix).astype(np.byte), fmt='%d', delimiter=' ')
		existingMatrices = sorted(os.listdir(PATH))
		matrices.append(matrix)
	return matrices[-1]
