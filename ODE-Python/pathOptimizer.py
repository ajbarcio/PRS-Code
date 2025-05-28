import numpy as np
from modules.PATHDEF import *
from scipy.optimize import minimize
from matplotlib import pyplot as plt

def main():
    path = RadiallyEndedPolynomial(1, 6)

    print(path.measure_length())
    plt.plot(path.get_neutralSurface(100))
    try:
        plt.show()
    except KeyboardInterrupt:
        return 0

if __name__ == "__main__":
    main()
