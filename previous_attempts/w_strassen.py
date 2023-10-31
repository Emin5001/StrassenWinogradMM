import numpy as np
from numpy.testing import assert_array_equal
from pypapi import papi_high
from pypapi import events as papi_events

papi_high.start_counters([
    papi_events.PAPI_FP_OPS,
    papi_events.PAPI_TOT_CYC,
])


def winograd(A, B):
    n = A.shape[0]
    
    if n == 1:
        return A * B 
    elif n == 2:
        # Handle 2x2 base case
        return np.array([[A[0,0]*B[0,0] + A[0,1]*B[1,0], A[0,0]*B[0,1] + A[0,1]*B[1,1]],
                         [A[1,0]*B[0,0] + A[1,1]*B[1,0], A[1,0]*B[0,1] + A[1,1]*B[1,1]]])
    
    A11, A12, A21, A22 = A[:n//2,:n//2], A[:n//2,n//2:], A[n//2:,:n//2], A[n//2:,n//2:]
    B11, B12, B21, B22 = B[:n//2,:n//2], B[:n//2,n//2:], B[n//2:,:n//2], B[n//2:,n//2:]
    
    P1 = winograd(A11 + A22, B11 + B22)
    P2 = winograd(A21 + A22, B11)        
    P3 = winograd(A11, B12 - B22)
    P4 = winograd(A22, B21 - B11)        
    P5 = winograd(A11 + A12, B22)
    P6 = winograd(A21 - A11, B11 + B12) 
    P7 = winograd(A12 - A22, B21 + B22)
    
    C11 = P1 + P4 - P5 + P7
    C12 = P3 + P5
    C21 = P2 + P4
    C22 = P1 - P2 + P3 + P6
    
    C = np.zeros((n, n))
    C[:n//2, :n//2] = C11
    C[:n//2, n//2:] = C12
    C[n//2:, :n//2] = C21 
    C[n//2:, n//2:] = C22
    
    return C

def test_2x2():
    A = np.array([[1, 2], [3, 4]]) 
    B = np.array([[5, 6], [7, 8]])
    
    C = winograd(A, B)
    
    expected = np.array([[19, 22], [43, 50]])
    print("C is: ")
    print(C)
    print("expected is: ")
    print(expected)

def test_4x4():
    A = np.random.rand(4,4)
    B = np.random.rand(4,4)
    
    C = winograd(A, B)
    D = A @ B 
    
    print("C is: ")
    print(C)
    print("D is: ")
    print(D)
    
def test_8x8():
    A = np.random.rand(8,8)
    B = np.random.rand(8,8)
     
    C = winograd(A, B)
    D = A @ B
    
    print("C is: ")
    print(C)
    print("D is: ")
    print(D)

test_2x2()
test_4x4()
test_8x8()