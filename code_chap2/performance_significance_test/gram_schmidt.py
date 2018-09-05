from numpy import linalg as la
def gram_schmidt(X):

	Q,_=la.qr(X)

	return Q