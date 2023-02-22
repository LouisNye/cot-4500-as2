import numpy as np
#Alg Q1
def nevilles_method(x, y, datapoint):
    #acquire length for iteration
    sizeofx = len(x)
    P = [[0] * sizeofx for num in range(sizeofx)]
    #filling in P
    for num in range(sizeofx):
        P[num][0] = y[num]
    #computation loop
    for i in range(1, sizeofx):
        for j in range(sizeofx - i):
            P[j][i] = ((datapoint - x[i + j]) * P[j][i - 1] - (datapoint - x[j]) * P[j + 1][i - 1]) / (x[j] - x[i + j])
    return P[0][sizeofx - 1]

#1
#setting up vars for passing as argument
x = [3.6, 3.8, 3.9]
y = [1.675, 1.436, 1.318]
datapoint = 3.7
#print
answer1 = nevilles_method(x, y, datapoint)
print(answer1)
print()


#Alg Q2
def forward_method(x, fx):
    #init length and list to append to with results, divtable to be used in Q3
    sizeofx = len(x)
    divtable = np.zeros((sizeofx, sizeofx))
    list = []
    for i in range(sizeofx):
        divtable[i][0] = fx[i]
    for i in range(1, sizeofx):
        for j in range(1, i + 1):
            divtable[i][j] = (divtable[i][j - 1] - divtable[i - 1][j - 1]) / (x[i] - x[i - j])
            if i == j:
                list.append(divtable[i][j])
    #print inside func def, loop for print to match expected output format
    for i in range(len(list)):
        print(list[i]) #values dont match?
    return divtable

#2
x = [7.2, 7.4, 7.5, 7.6]
fx = [23.5492, 25.3913, 26.8224, 27.4589]
#store divided difference table for Q3
divtable = forward_method(x, fx) 
print()


#Alg Q3
def approximate_result(divtable, x, targetval, yknot):
    #init reocurring vars
    recx = 1
    recpx = yknot
    #compute
    for i in range(1, len(divtable)):
        poly = divtable[i][i]
        recx *= (targetval - x[i - 1])
        op = poly * recx
        recpx += op
    #printing inside func def
    print(recpx)

#3
#print in func
approximate_result(divtable, x, 7.3, fx[0])
print()


#Alg Q4
def hermite_interpolation(x, fx, fp_of_x):
    #initialize length
    sizeofx = (2 * len(x))
    #init zero matrix to plug in givens
    m = np.zeros((sizeofx, sizeofx + 1))
    for i in range(sizeofx):
        m[i][0] = x[i // 2]
        m[i][1] = fx[i // 2]
    for i in range(1, sizeofx, 2):
        m[i][2] = fp_of_x[i // 2]
    #compute
    for i in range(2, sizeofx):
        for j in range(2, i + 2):
            if m[i][j] != 0.:
                continue           
            m[i][j] = (m[i][j - 1] - m[i - 1][j - 1]) / (m[i][0] - m[i - j + 1][0])
    #drop column
    m = np.delete(m, sizeofx, 1)
    #print options to match expected output format
    np.set_printoptions(precision=7, suppress=False, linewidth=100)
    #printing inside func def
    print(np.matrix(m))

#4
#init givens as args
x = [3.6, 3.8, 3.9]
fx = [1.675, 1.436, 1.318]
fp_of_x = [-1.195, -1.188, -1.182]
#passing args, print in func def
hermite_interpolation(x, fx, fp_of_x)
print()


#Alg Q5
def cubic_spline_interpolation(x, fx):
    #init length
    sizeofx = len(x)
    #init diff & zeros matrix (matrix A) to plug in
    diffx = np.diff(x)
    A = np.zeros((sizeofx, sizeofx))
    A[0, 0] = 1
    A[sizeofx-1, sizeofx-1] = 1
    #compute for parts a,b, c and print
    for i in range(1, sizeofx-1):
        A[i, i-1] = diffx[i-1]
        A[i, i] = 2 * (diffx[i-1] + diffx[i])
        A[i, i+1] = diffx[i]
    #print matrix A
    print(A)
    #compute vector b
    b = np.zeros(sizeofx)
    for i in range(1, sizeofx-1):
        b[i] = (3/diffx[i]) * (fx[i + 1]- fx[i]) - (3/diffx[i-1]) * (fx[i] - fx[i - 1])
    #print vector b
    print(b)
    #compute and print vector x
    x = np.linalg.solve(A, b)
    print(x)

#5
x = np.array([2, 5, 8, 10])
fx = np.array([3, 5, 7, 9])
#printing inside func def, passing args
cubic_spline_interpolation(x, fx)