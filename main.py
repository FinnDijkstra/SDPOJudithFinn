import math
def lp_sample(n, theta, d, npoints, filename):
    # S = [-1 + math.cos(i * (theta- math.pi)/(npoints-1) for i in range(npoints)]
    S = [math.cos(math.pi + i/(npoints-1) * (-math.pi +theta)) for i in range(npoints)]

    with open(filename, "w") as out:
        #number of constraints
        out.write("%d\n" % npoints)

        #number of blocks and sizes
        out.write("2\n")
        out.write("-%d -%d\n" % (d+1, npoints))

        # The right-hand side of the linear constraints.
        out.write("-1.0" + (npoints - 1)* " -1.0" + "\n")

        #the constraints
        for i in range(npoints):
            Plist = jacobi_List(n, d, S[i])
            for j in range(d+1):
                P = Plist[j]
                out.write("%d 1 %d %d %f\n"% ((i+1), (j+1), (j+1), P))
            out.write("%d 2 %d %d 1.0\n" % ((i+1), (i+1), (i+1)))

        #objective function
        for j in range(1, d+2):
            out.write("0 1 %d %d -1.0\n" %(j,j))









def jacobi_(n,k,u):
    alpha = (n-3)
    if k == 0:
        return 1
    if k == 1:
        return u
    leftP = jacobi_(n,k-1,u)
    rightP = jacobi_(n,k-2,u)
    P = (2*k+alpha-1)/(k+alpha)*u*leftP -(k-1)/(k+alpha)*rightP
    print(f"{P} at {k}")
    return P

def jacobi_List(n,k,u):
    Plist = [1.0,u]
    for i2 in range(2,k+1):
        alpha = (n - 3)
        leftP = Plist[i2-1]
        print(leftP)
        rightP = Plist[i2-2]
        print(rightP)
        P = (2 * k + alpha - 1) / (k + alpha) * u * Plist[i2-1] - (k - 1) / (k + alpha) * Plist[i2-2]
        Plist.append(P)
    return Plist


# def polymin_sdpa(filename, p):
#     """Write an SDPA file with an SDP whose optimal value is the minimum of p.
#
#     The polynomial p is given as a list of coefficients.  For instance, the
#     polynomial p(x) = 2x^2 - 1 is given by the list [-1, 0, 2].
#
#     This function writes to a file with name `filename` an SDP in SDPA format
#     whose optimal value is the global minimum of p.
#
#     To understand the setup, read the file `polymin.md`.
#
#     """
#     if len(p) % 2 != 1:
#         raise ValueError("polynomial must have even degree")
#
#     # Our polynomial has degree 2d.
#     d = len(p) // 2
#
#     # Generate the SDPA file.
#     with open(filename, "w") as out:
#         # Number of constraints.
#         out.write("%d\n" % (2 * d + 1))
#
#         # Number of blocks and block sizes.
#         out.write("2\n")
#         out.write("%d -2\n" % (d + 1))
#
#         # The right-hand side of the linear constraints.
#         out.write(" ".join(map(str, p)) + "\n")
#
#         # The constraint <F_1, Q> + <B, L> = p[0].
#         out.write("1 1 1 1 1.0\n")
#         out.write("1 2 1 1 1.0\n")
#         out.write("1 2 2 2 -1.0\n")
#
#         # The constraints <F_k, Q> = p[k - 1].
#         for k in range(1, len(p)):
#             for i in range(0, d + 1):
#                 # We need to have i + j == k, so j = k - i.  But we also want to
#                 # have j >= i, because we only want to give the upper-diagonal
#                 # entries.  We also need to have j <= d.
#                 j = k - i
#
#                 if i <= j <= d:
#                     out.write("%d 1 %d %d 1.0\n" % (k + 1, i + 1, j + 1))
#
#         # The objective function.
#         out.write("0 2 1 1 1.0\n")
#         out.write("0 2 2 2 -1.0\n")


def main():
    # Let's try our program with 4.0*x^4 + 2.0*x^3 - 3.0*x^2 + 2.0.  This should
    # give an SDP with optimal solution 0.688022543737423, which is really the
    # minimum of the polynomial!
    print(jacobi_List(8,10,1/2))
    jacobi_(8,10,1/2)
    lp_sample(8, math.pi/3.0, 30, 10, "test_1.sdpa")

main()


