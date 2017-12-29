#plot the mass loss

import numpy as np
import matplotlib.pyplot as pyplot

def readfile(myfile):
    datafile = open(myfile,"r")
    return [float(j) for j in datafile.readline().split(",")]
    # time = []
    # mass_loss = []
    # for row in datafile:
    #     vars = [float(j) for j in row.split(",")]
    #     time.append(vars[0])
    #     mass_loss.append(vars[1])
    # return [time, mass_loss]

result1 = readfile("mass_loss_1.txt")
result2 = readfile("mass_loss_2.txt")
result3 = readfile("mass_loss_3.txt")

deltaTime = 0.1
times = np.arange(0, 50, deltaTime)

pyplot.figure(figsize = (10,10))
pyplot.plot(times, result1, label = "p = 10")
pyplot.plot(times, result2, label = "p = 20")
pyplot.plot(times, result3, label = "p = 40")

pyplot.xlabel("Time (s)", color = "black", fontsize=14)
pyplot.ylabel("Mass loss", color = "black", fontsize=14)
pyplot.title("Mass loss", fontsize=20)
pyplot.legend()

pyplot.savefig("/home/petiaccja/UppsalaUniversity/UU_AppliedFEM/images/p3_ml_all.png")


result1 = readfile("mass_loss_optimal.txt")
times = np.arange(0, 35, deltaTime)

pyplot.figure(figsize = (10,10))
pyplot.plot(times, result1, label = "optimal")

pyplot.xlabel("Time (s)", color = "black", fontsize=14)
pyplot.ylabel("Mass loss", color = "black", fontsize=14)
pyplot.title("Mass loss", fontsize=20)
pyplot.legend()
pyplot.plot([5, 7, 30], [10, 15, 30], "ro", label="target")

pyplot.savefig("/home/petiaccja/UppsalaUniversity/UU_AppliedFEM/images/p3_ml_optimal.png")