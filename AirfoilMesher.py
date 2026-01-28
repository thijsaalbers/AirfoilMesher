import sys
import numpy as np
import matplotlib.pyplot as plt

from GMSHunstructured import createUnstructuredGridGMSH 
from EllipticStructuredOgrid import createEllipticStructuredOGrid


def parseArguments():
    """
    Parse command-line arguments.
    """
    # Sanity checks for all input arguments
    if len(sys.argv) < 3:
        raise RuntimeError("\n\nNot enough input arguments provided, please look at the user guide\n\n")

    # Argument 1 is the airfoil code which should be in format NACAXXXX
    airfoilCode = sys.argv[1].upper()
    if not airfoilCode.startswith("NACA"):
        raise RuntimeError("\n\nOnly NACA airfoils are supported\n\n")

    # Argument 2 says which meshing algorithm should be used
    meshingMethod = sys.argv[2].upper()

    return airfoilCode, meshingMethod


def naca4Geometry(code, nPoints=20):
    """
    Generate NACA 4-digit airfoil coordinates.
    """
    m = int(code[0]) / 100.0
    p = int(code[1]) / 10.0
    t = int(code[2:]) / 100.0

    finalCoeff = -0.1036 # only consider a closer geometry for now, for open geometry the result is -0.1015

    # Code for the camber line of the airfoil (naca 4 digit series)
    def yc(x):
        if m == 0:
            return np.zeros_like(x)  

        if p == 0:
            return m * x

        if p == 1:
            return m * np.ones_like(x)

        y = np.where(
            x <= p,
            m/p**2 * (2*p*x - x**2),
            m/(1-p)**2 * ((1 - 2*p) + 2*p*x - x**2)
        )
        return y

    def yt(x):
        return 5 * t * (0.2969 * np.sqrt(x) - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 + finalCoeff*x**4)

    theta = np.linspace(0, np.pi, nPoints)
    xCoords = 0.5*(1 - np.cos(theta)) # Focus more points near the start and end of the airfoil

    yCamber = yc(xCoords)
    thickness = yt(xCoords)

    yUpper = yCamber + thickness
    yLower = yCamber - thickness

    return xCoords, yUpper, yLower, yCamber


def plotAirfoil(xCoords, yUpper, yLower, yCamber, airfoilCode):
    """
    Plot the airfoil and camber line.
    """
    plt.figure(figsize=(10, 6))
    plt.plot(xCoords, yUpper, color='blue', label='Airfoil shape')
    plt.scatter(xCoords, yUpper, color='blue', s = 5)
    plt.plot(xCoords, yLower, color='blue')
    plt.scatter(xCoords, yLower, color='blue', s = 5)
    plt.plot(xCoords, yCamber, color='red', label='Camber line')
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title(f"{airfoilCode} airfoil point distribution")
    plt.axis('equal')
    plt.legend()
    plt.savefig("geometryPlots/" +airfoilCode + ".png", dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":

    # Start code by parsing all command line arguments
    airfoilCode, meshingMethod = parseArguments()
    
    numberPart = airfoilCode[4:]
    if not numberPart.isdigit():
        raise RuntimeError("\n\nThe provided airfoil code is not valid as it does not consist of only numbers\n\n")

    # Define airfoil geometry
    xCoords = yUpper = yLower = yCamber = None
    if len(numberPart) == 4:
        xCoords, yUpper, yLower, yCamber = naca4Geometry(numberPart)
    elif len(numberPart) > 4:
        raise RuntimeError('\n\nSupport for airfoil geometries of series other than the NACA 4-series have not yet been implemented\n\n')

    # Plot the airfoil geometry as parsed from the provided airfoil code and save the result as a .png
    plotAirfoil(xCoords, yUpper, yLower, yCamber, airfoilCode)

    # Create a mesh for the airfoil geometry
    if meshingMethod == 'GMSH':
        createUnstructuredGridGMSH(xCoords, yUpper, yLower)
    elif meshingMethod == 'ELLIPTIC':
        createEllipticStructuredOGrid(xCoords, yUpper, yLower)
    else:
        raise RuntimeError('\n\nProvided meshing method \'' + meshingMethod + '\' is not a valid option, please look at the user guide\n\n')
