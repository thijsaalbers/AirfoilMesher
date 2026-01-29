import sys
import numpy as np
import matplotlib.pyplot as plt
import argparse

from GMSHunstructured import createUnstructuredGridGMSH 
from EllipticStructuredOgrid import createEllipticStructuredOGrid

def NACA4(code, nodeSpacing, Nc):
    """
    Generate NACA 4-digit airfoil coordinates.
    """

    # Parge the airfoil code digits
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

    # Define the number of points required for the mesh according to the total circumferential amount of points
    # note that the leading edge and trailing edge 
    nPoints = int(Nc / 2) + 1

    xCoords = None
    if nodeSpacing == "UNIFORM":
        xCoords = np.linspace(0.0, 1.0, nPoints)
    elif nodeSpacing == "COSLE":
        theta = np.linspace(0, np.pi, nPoints)
        xCoords = (1 - np.cos(theta)) # Focus more points near the leading edge of the airfoil
    elif nodeSpacing == "COSLETE":
        theta = np.linspace(0, np.pi, nPoints)
        xCoords = 0.5*(1 - np.cos(theta)) # Focus more points near both leading- and trailing edges of the airfoil
    else:
        raise RuntimeError("\n\nNot valid node spacing method in NACA geometry definition\n\n")

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
    plt.title(f"NACA{airfoilCode} airfoil point distribution")
    plt.axis('equal')
    plt.legend()
    plt.savefig("geometryPlots/NACA" +airfoilCode + ".png", dpi=300, bbox_inches="tight")
    plt.close()


if __name__ == "__main__":

    # Start by parsing input arguments using the argparse library 
    parser = argparse.ArgumentParser()

    # Airfoil code from NACA series
    parser.add_argument("--NACA", 
        type=str, 
        default="0012", 
        help="NACA airfoil code"
    )


    # Meshing method 
    parser.add_argument(
        "--meshingMethod",
        type=str.upper,
        default="GMSH",
        choices=["GMSH", "ELLIPTIC"],
        help="The method chosen for generating the mesh"
    )

    # Node spacing of airfoil
    parser.add_argument(
        "--nodeSpacing",
        type=str.upper,
        choices=["UNIFORM", "COSLE", "COSLETE"],
        default="COSLETE",
        help=(
            "Node spacing distribution along the airfoil chord:\n"
            "  UNIFORM: uniform node spacing;"
            "  COSLE: cosine clustering toward the leading edge;"
            "  COSLETE: cosine clustering toward both leading- and trailing edges;"
        )
    )

    parser.add_argument("--Nc", 
        type=int, 
        default=40, 
        help="Number of circumferential points used to define the airfoil geometry, MUST BE AN EVEN NUMBER.  "
            "Please note that more points define a more accurate airfoil however it has no direct influence on the "
            "unstructured automaticaly generated mesh resolution. For the elliptic grid generation this "
            "DOES define the amount of nodes in circumferential direction making it important for mesh generation."
    )

    parser.add_argument("--Nr", 
        type=int, 
        default=40, 
        help="Number of points in radial direction, note that this also includes the boundaries so it must be greater than 2."
            "This option is only relevant to the O-grid generation and is disregarded for automatic gmsh mesh generation"
    )

    # Read arguments
    args = parser.parse_args()
    NACAdigits=args.NACA
    meshingMethod=args.meshingMethod
    nodeSpacing=args.nodeSpacing
    Nc=args.Nc
    Nr=args.Nr

    # Sanity checks for defining airfoil geometry
    if not NACAdigits.isdigit():
        raise RuntimeError("\n\nThe provided airfoil code is not valid as it does not consist of only numbers\n\n")

    # Define airfoil geometry
    xCoords = yUpper = yLower = yCamber = None
    if len(NACAdigits) < 4:
        raise RuntimeError('\n\nNACA airfoil series need at least 4 digits\n\n')
    elif len(NACAdigits) == 4:
        xCoords, yUpper, yLower, yCamber = NACA4(NACAdigits, nodeSpacing, Nc)
    elif len(NACAdigits) > 4:
        raise RuntimeError('\n\nSupport for airfoil geometries of series other than the NACA 4-series have not yet been implemented\n\n')

    # Plot the airfoil geometry as parsed from the provided airfoil code and save the result as a .png
    # This way the user can visually confirm the correct airfoil geometry is used
    plotAirfoil(xCoords, yUpper, yLower, yCamber, NACAdigits)

    # Create a mesh for the airfoil geometry
    if meshingMethod == "GMSH":
        createUnstructuredGridGMSH(xCoords, yUpper, yLower)
    elif meshingMethod == "ELLIPTIC":
        createEllipticStructuredOGrid(xCoords, yUpper, yLower, Nr)
    else:
        raise RuntimeError('\n\nProvided meshing method \'' + meshingMethod + '\' is not a valid option, please look at the user guide\n\n')
