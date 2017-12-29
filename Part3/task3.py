# !/usr/bin/python

from dolfin import *
from math import *
import numpy
import scipy.optimize as opt

# Create mesh and function space
mesh = Mesh("circle2.xml")

# Define parameters
endTime = 20
deltaTime = 0.1
alpha = 0.01
r = 0.2
R = 0.5
rho = 10

# Create subdomain for Dirichlet boundary
class DirichletBoundary(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary
    
# Define initial data
class InitialCondition(Expression):
    def eval(self, value, x):
        norm = hypot(x[0], x[1])
        circle = numpy.zeros(len(x))
        circle[0], circle[1] = x[0], x[1]
        if norm == 0:
            circle[0] = 1
        else:
            circle[0] *= R/norm
            circle[1] *= R/norm
        circle = numpy.subtract(circle, x)
        if sqrt(numpy.dot(circle, circle)) < r:
            value[0] = rho
        else:
            value[0] = 0.0

def Simulate(outFileName, outLossFileName):
    # Set up simulation
    functionSpace = FunctionSpace(mesh, "CG", 1)

    # Set up boundary condition
    uBoundary = Constant(0.0)
    uBoundary.t = 0
    boundaryCond = DirichletBC(functionSpace, uBoundary, DirichletBoundary())

    # Set initial data
    uPrev = Function(functionSpace)
    initialCondition = InitialCondition(element = functionSpace.ufl_element())
    uPrev.interpolate(initialCondition)
    uInitial = Function(functionSpace)
    uInitial.assign(uPrev)

    # Define variational problem
    u = TrialFunction(functionSpace)
    v = TestFunction(functionSpace)
    a = u*v*dx + 0.5*alpha*deltaTime*inner(grad(u), grad(v))*dx
    L = uPrev*v*dx - 0.5*deltaTime*alpha*inner(grad(uPrev), grad(v))*dx

    # Save results
    if not (outFileName is None):
        outFile = File(outFileName)
        outFile << uPrev

    # Prepare problem
    A = assemble(a)

    currentTime = deltaTime
    uk = Function(functionSpace)

    massLossFunctional = (uInitial-uk)*dx
    massLoss = []


    while currentTime <= endTime:
        # Get solution
        b = assemble(L)
        uBoundary.t = currentTime
        boundaryCond.apply(A, b)
        solve(A, uk.vector(), b)

        # Prepare next iteration
        currentTime += deltaTime
        uPrev.assign(uk)

        # Save data
        massLoss.append(assemble(massLossFunctional))
        if not (outFileName is None):
            outFile << uPrev

    if not (outLossFileName is None):
        outLossFile = open(outLossFileName, 'w')
        outLossFile.write("0,")
        outLossFile.write(",".join([str(x) for x in massLoss]))

    return massLoss
    

# Part C.1
print("Solving part C.1...")
rho = 10
R = 0.5
r = 0.2
endTime = 20
mesh = Mesh("circle3.xml")
Simulate("./solution/circle.pvd", None)
mesh = Mesh("sphere2.xml")
Simulate("./solution/sphere.pvd", None)

# Part C.2
print("Solving part C.2...")
mesh = Mesh("sphere2.xml")
R = 0.5
r = 0.2
endTime = 50
rho = 10
Simulate(None, "mass_loss_1.txt")
rho = 20
Simulate(None, "mass_loss_2.txt")
rho = 40
Simulate(None, "mass_loss_3.txt")

# Part C.3
print("Solving part C.3...")
mesh = Mesh("sphere1.xml")
endTime = 35
def CostFunction(x):
    global rho, R, r
    rho = x[0]
    R = x[1]
    r = x[2]

    massLoss = Simulate(None, None)

    keyTimes = [5, 7, 30]
    keyDoses = [10, 15, 30]
    cost = 0
    for i in range(0,3):
        keyLoss = massLoss[int(keyTimes[i]/deltaTime)]
        cost += (keyLoss - keyDoses[i])**2

    return cost

optimOptions = {'xtol': 1e-3, 'disp': True}
optBounds = [(0, None), (0, None), (0, 1)] # "L-BFGS-B"
optimal = opt.minimize(CostFunction, [20, 0.5, 0.1], method='nelder-mead', options=optimOptions).x #, bounds=optBounds).x
print("Optimial solution is: {}".format(optimal))
rho = optimal[0]
R = optimal[1]
r = optimal[2]
Simulate(None, "mass_loss_optimal.txt")

