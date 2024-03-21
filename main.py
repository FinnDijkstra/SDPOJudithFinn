##### Finn Dijkstra & Judith Capel #######
import math
def lp_sample(n, theta, d, npoints, filename):
    #S = [-1 + i*(math.cos( theta )+1)/(npoints-1) for i in range(npoints)]
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
                out.write("%d 1 %d %d %s\n"% ((i+1), (j+1), (j+1), P))
            out.write("%d 2 %d %d 1.0\n" % ((i+1), (i+1), (i+1)))

        #objective function
        for j in range(1, d+2):
            out.write("0 1 %d %d -1.0\n" %(j,j))

# def jacobi_(n,k,u):
#     alpha = (n-3)
#     if k == 0:
#         return 1
#     if k == 1:
#         return u
#     leftP = jacobi_(n,k-1,u)
#     rightP = jacobi_(n,k-2,u)
#     P = (2*k+alpha-1)/(k+alpha)*u*leftP -(k-1)/(k+alpha)*rightP
#     print(f"{P} at {k}")
#     return P

def jacobi_List(n,k,u):
    Plist = [1.0,u]
    for i2 in range(2,k+1):
        alpha = (n - 3)
        leftP = Plist[i2-1]
        rightP = Plist[i2-2]
        P = (2*i2+alpha-1)/(i2+alpha)*u*leftP -(i2-1)/(i2+alpha)*rightP
        Plist.append(P)
    return Plist

def jacobi_coef(n,d):
    alpha = (n-3)/2
    C= []
    for l in range(d+1):
        C.append((d+1)*[0])
    C[0][0] = 1
    C[1][1] =1
    for k in range(2,d+1):
        A = ((2*k)+2*alpha -1)/(k+2*alpha)
        B = -((k-1)/(k+2*alpha))
        for j in range(d+1):
            if j == 0:
                C[k][j] = B*C[k-2][j]
            else:
                C[k][j]= A*C[k-1][j-1] + B* C[k-2][j]

    return C

def lp_sos(n, theta, d, filename):
    with open(filename, "w") as out:
        #number of constraints
        out.write("%d\n" % (d+1))

        #number of blocks and sizes
        out.write("3\n")
        out.write("-%d %d %d\n" % ((d+1), d//2+1, d//2+1))

        # The right-hand side of the linear constraints.
        out.write("-1.0" + (d) * " 0.0" + "\n")


        #the constraints
        C = jacobi_coef(n,d)
        for i in range(d+1):
            for k in range(d+1):
                out.write("%d 1 %d %d %.12f\n" % ((i + 1), (k + 1), (k + 1), C[k][i]))
            for j in range(d//2+1):
                j2 = i-j
                k = i-1-j
                if j2>=0:
                    if j2 <= j:
                        out.write("%d 2 %d %d %.12f\n" % ((i + 1), (j2 + 1), (j + 1),(1)))
                        out.write("%d 3 %d %d %.12f\n" % ((i + 1), (j2 + 1), (j + 1), (math.cos(theta))))
                if k >=0:
                    if k<=j:
                        out.write("%d 2 %d %d %.12f\n" % ((i + 1), (k + 1), (j + 1), 1))
                        out.write("%d 3 %d %d %.12f\n" % ((i + 1), (k + 1), (j + 1), -1))

        #objective function
        for j in range(1, d+2):
            out.write("0 1 %d %d -1.0\n" %(j,j))



def main():
    lp_sample(24, math.pi/3.0, 19, 1000, "test_1.sdpa")
    lp_sos(24,math.pi/3.0, 19,"test_sos.sdpa")
main()


