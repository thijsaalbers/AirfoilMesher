import numpy as np
import matplotlib.pyplot as plt
import gmsh

def exportToGmsh(x, y):

    gmsh.initialize()
    gmsh.model.add("grid")

    Nc, Nr = x.shape

    # Create discrete entities
    surf_tag = 1
    airfoil_tag = 2
    farfield_tag = 3

    gmsh.model.addDiscreteEntity(2, surf_tag)
    gmsh.model.addDiscreteEntity(1, airfoil_tag)
    gmsh.model.addDiscreteEntity(1, farfield_tag)

    # Add nodes
    node_tags = []
    node_coords = []

    for j in range(Nr):
        for i in range(Nc):
            node_id = j * Nc + i + 1
            node_tags.append(node_id)
            node_coords.extend([x[i, j], y[i, j], 0.0])

    gmsh.model.mesh.addNodes(2, surf_tag, node_tags, node_coords)

    # Add quad elements (surface)
    quad_type = 3
    quad_tags = []
    quad_nodes = []

    elem_id = 1

    for j in range(Nr - 1):
        for i in range(Nc):
            ip = (i + 1) % Nc

            n1 =  j    * Nc + i  + 1
            n2 =  j    * Nc + ip + 1
            n3 = (j+1) * Nc + ip + 1
            n4 = (j+1) * Nc + i  + 1

            quad_tags.append(elem_id)
            quad_nodes.extend([n1, n2, n3, n4])
            elem_id += 1

    gmsh.model.mesh.addElements(
        2, surf_tag,
        [quad_type],
        [quad_tags],
        [quad_nodes]
    )

    # Add boundary line elements
    line_type = 1

    airfoil_elems = []
    airfoil_nodes = []

    farfield_elems = []
    farfield_nodes = []

    # Airfoil (j = 0)
    for i in range(Nc):
        ip = (i + 1) % Nc

        n1 = i + 1
        n2 = ip + 1

        airfoil_elems.append(elem_id)
        airfoil_nodes.extend([n1, n2])
        elem_id += 1

    # Farfield (j = Nr-1)
    offset = (Nr - 1) * Nc
    for i in range(Nc):
        ip = (i + 1) % Nc

        n1 = offset + i  + 1
        n2 = offset + ip + 1

        farfield_elems.append(elem_id)
        farfield_nodes.extend([n1, n2])
        elem_id += 1

    gmsh.model.mesh.addElements(
        1, airfoil_tag,
        [line_type],
        [airfoil_elems],
        [airfoil_nodes]
    )

    gmsh.model.mesh.addElements(
        1, farfield_tag,
        [line_type],
        [farfield_elems],
        [farfield_nodes]
    )

    # Physical groups
    gmsh.model.addPhysicalGroup(2, [surf_tag], 1)
    gmsh.model.setPhysicalName(2, 1, "fluid")

    gmsh.model.addPhysicalGroup(1, [airfoil_tag], 2)
    gmsh.model.setPhysicalName(1, 2, "airfoil")

    gmsh.model.addPhysicalGroup(1, [farfield_tag], 3)
    gmsh.model.setPhysicalName(1, 3, "farfield")

    gmsh.fltk.run()

def createEllipticStructuredOGrid(xCoords, yUpper, yLower):

    # Grid parameters
    Rfarfield = 20.0 # Outer circle (farfield) in amount of chord lengths
    Nr = 40        # radial points
    # Circumferential points are defined by the airfoil geometry that is provided

    # Gauss-Seidel iteration parameters
    max_iter = 5000
    tolerance = 1e-6

    # Make airfoil geometry suitable for O-grid by removing repeated nodes and making it a singular CCW loop
    xUpper = xCoords[::-1]
    xLower = xCoords[1:-1] # lower points neglect the leading edge and trailing edge points as these are already included in the upper points

    x_airfoil = np.append(xUpper, xLower)
    y_airfoil = np.append(yUpper[::-1], yLower[1:-1])

    # Number of circumferential points
    Nc = len(x_airfoil) 

    thetaFarfield = np.linspace(0, 2.0*np.pi, int(Nc), endpoint = False)

    # Initialize grid arrays
    x = np.zeros((Nc, Nr))
    y = np.zeros((Nc, Nr))

    # Set inner boundary (airfoil surface)
    x[:, 0] = x_airfoil
    y[:, 0] = y_airfoil

    # Set outer boundary (far-field circle)
    x[:, -1] = Rfarfield * np.cos(thetaFarfield)+0.5 # add 0.5 to center the airfoil at x=0 to x=1 in the middle
    y[:, -1] = Rfarfield * np.sin(thetaFarfield)

    # Use linear interpolation for a proper initial guess to speed up convergence
    for j in range(1, Nr-1):
        s = j / (Nr-1)
        x[:, j] = (1-s)*x[:,0] + s*x[:, -1]
        y[:, j] = (1-s)*y[:,0] + s*y[:, -1]


    # Iterative Gauss-Seidel algorithm of solving the system of elliptical equations defined by Thomas 1974
    for it in range(max_iter):
        x_old = x.copy()
        y_old = y.copy()

        # -1 loops back to the end in Python thus the upper loop automatically enforces periodicity
        for i in range(Nc):
            ip = (i+1) % Nc
            im = (i-1) % Nc
            for j in range(1, Nr-1):
                jp = j+1
                jm = j-1

                # Central finite differences
                x_xi = 0.5*(x[ip,j] - x[im,j])
                x_eta = 0.5*(x[i,jp] - x[i,jm])
                y_xi = 0.5*(y[ip,j] - y[im,j])
                y_eta = 0.5*(y[i,jp] - y[i,jm])

                # Metric coefficients
                alpha = x_xi**2 + y_xi**2
                beta  = x_xi*x_eta + y_xi*y_eta
                gamma = x_eta**2 + y_eta**2

                # Cross derivatives
                x_xi_eta = 0.25*(x[ip,jp] - x[ip,jm] - x[im,jp] + x[im,jm])
                y_xi_eta = 0.25*(y[ip,jp] - y[ip,jm] - y[im,jp] + y[im,jm])

                # Laplace-like update with metrics (explicit Gauss-Seidel)
                denom = 2*(alpha + gamma)
                x_new = ((alpha*(x[i,jp]+x[i,jm]) + gamma*(x[ip,j]+x[im,j]) - 2*beta*x_xi_eta)/denom)
                y_new = ((alpha*(y[i,jp]+y[i,jm]) + gamma*(y[ip,j]+y[im,j]) - 2*beta*y_xi_eta)/denom)

                # Succesive over relaxation
                omega = 1.8
                x[i,j] = (1-omega)*x[i,j] + omega*x_new
                y[i,j] = (1-omega)*y[i,j] + omega*y_new

        # Convergence check
        err = np.max(np.abs(x - x_old) + np.abs(y - y_old))
        if err < tolerance:
            print(f"Gauss-Seidel iteration has converged at iteration {it}")
            break

    exportToGmsh(x,y)

    # Postprocessing by first adding the first element of x and y to the end of each array to plot the completely looped O-grid
    # x = np.vstack([x, x[0:1, :]])  # add first row at the end
    # y = np.vstack([y, y[0:1, :]])  # same for y
    # plt.figure(figsize=(20,20))
    # for j in range(Nr):
    #     plt.plot(x[:, j], y[:, j], 'k', linewidth = 0.1)  # radial lines
    # for i in range(Nc):
    #     plt.plot(x[i, :], y[i, :], 'k', linewidth = 0.1)  # circumferential lines
    # plt.axis('equal')
    # plt.savefig("grid.png", dpi=600, bbox_inches="tight")
    # plt.close()
