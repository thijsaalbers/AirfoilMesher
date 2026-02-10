import numpy as np

def NACA4(code, nodeSpacing, Nc):
    """
    Generate NACA 4-digit airfoil coordinates.
    """

    # Parge the airfoil code digits
    m = int(code[0]) / 100.0
    p = int(code[1]) / 10.0
    t = int(code[2:4]) / 100.0

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


# Work in progress
# def NACA5(code, nodeSpacing, Nc):
#     """
#     Generate NACA 5-digit airfoil coordinates.
#     """

#     # Parge the airfoil code digits
#     L = int(code[0]) * 3.0 / 20.0
#     P = int(code[1]) / 20.0
#     Q = int(code[2]) 
#     XX = int(code[3:6]) / 100.0

#     if Q==0: # Standard version


#     elif Q==1: # Reflec version


#     else:
#         raise RuntimeError("The 3rd digit of the NACA 5 digit series airfoil code must be either 0 (standard) or 1 (reflex), the current version is invalid")



#     finalCoeff = -0.1036 # only consider a closer geometry for now, for open geometry the result is -0.1015

#     # Code for the camber line of the airfoil (naca 4 digit series)
#     def yc(x):
#         if m == 0:
#             return np.zeros_like(x)  

#         if p == 0:
#             return m * x

#         if p == 1:
#             return m * np.ones_like(x)

#         y = np.where(
#             x <= p,
#             m/p**2 * (2*p*x - x**2),
#             m/(1-p)**2 * ((1 - 2*p) + 2*p*x - x**2)
#         )
#         return y

#     def yt(x):
#         return 5 * t * (0.2969 * np.sqrt(x) - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 + finalCoeff*x**4)

#     # Define the number of points required for the mesh according to the total circumferential amount of points
#     # note that the leading edge and trailing edge 
#     nPoints = int(Nc / 2) + 1

#     xCoords = None
#     if nodeSpacing == "UNIFORM":
#         xCoords = np.linspace(0.0, 1.0, nPoints)
#     elif nodeSpacing == "COSLE":
#         theta = np.linspace(0, np.pi, nPoints)
#         xCoords = (1 - np.cos(theta)) # Focus more points near the leading edge of the airfoil
#     elif nodeSpacing == "COSLETE":
#         theta = np.linspace(0, np.pi, nPoints)
#         xCoords = 0.5*(1 - np.cos(theta)) # Focus more points near both leading- and trailing edges of the airfoil
#     else:
#         raise RuntimeError("\n\nNot valid node spacing method in NACA geometry definition\n\n")

#     yCamber = yc(xCoords)
#     thickness = yt(xCoords)

#     yUpper = yCamber + thickness
#     yLower = yCamber - thickness

#     return xCoords, yUpper, yLower, yCamber