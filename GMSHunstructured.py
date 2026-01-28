
import gmsh


def createUnstructuredGridGMSH(xCoords, yUpper, yLower):

    # Variables determining mesh
    r = 20 # Farfield radius
    lc = 2.5  # mesh size
    LcMin = 0.002 * lc # Min mesh size near airfoil surface
    LcMinDist = 0.002 * lc # Distance of minimal mesh size near airfoil surface
    meshAlgorithm = 6 # 6 = frontal delaunay
    meshOrder = 1
    recombinationAlgorithm = 0 # 0 = less recombinations , 1 = more recombinations
    recombine = 0 # 0 = no recombination, 1 = recombine

    # Start of GMSH mesh creation
    gmsh.initialize()
    gmsh.model.add("Unstructured_Airfoil_mesh")

    # Create airfoil surface

    upperPoints = []
    lowerPoints = []

    # Leading edge
    pLE = gmsh.model.geo.addPoint(xCoords[0], yUpper[0], 0, lc)
    upperPoints.append(pLE)

    # Upper surface interior
    for xi, yi in zip(xCoords[1:], yUpper[1:]):
        upperPoints.append(
            gmsh.model.geo.addPoint(float(xi), float(yi), 0.0, lc)
        )

    # Trailing edge
    pTE = upperPoints[-1]

    # Lower surface interior (reverse)
    for xi, yi in zip(xCoords[-2:0:-1], yLower[-2:0:-1]):
        lowerPoints.append(
            gmsh.model.geo.addPoint(float(xi), float(yi), 0.0, lc)
        )

    lowerPoints.append(pLE)

    upperSpline = gmsh.model.geo.addSpline(upperPoints)
    lowerSpline = gmsh.model.geo.addSpline([pTE] + lowerPoints)

    airfoil_loop = gmsh.model.geo.addCurveLoop([
        upperSpline,
        lowerSpline
    ])

    # Create farfield
    
    cx, cy = 0.0, 0.0  # center

    # circle points (4 quarters)
    pC1 = gmsh.model.geo.addPoint(cx + r, cy, 0, lc)
    pC2 = gmsh.model.geo.addPoint(cx, cy + r, 0, lc)
    pC3 = gmsh.model.geo.addPoint(cx - r, cy, 0, lc)
    pC4 = gmsh.model.geo.addPoint(cx, cy - r, 0, lc)

    # center point (for arcs)
    center = gmsh.model.geo.addPoint(cx, cy, 0, lc)

    # create 4 quarter-circle arcs
    a1 = gmsh.model.geo.addCircleArc(pC1, center, pC2)
    a2 = gmsh.model.geo.addCircleArc(pC2, center, pC3)
    a3 = gmsh.model.geo.addCircleArc(pC3, center, pC4)
    a4 = gmsh.model.geo.addCircleArc(pC4, center, pC1)

    circle_loop = gmsh.model.geo.addCurveLoop([a1, a2, a3, a4])
    circle_surface = gmsh.model.geo.addPlaneSurface([circle_loop, airfoil_loop])  # subtract square


    # Create names
    airfoil_tag = gmsh.model.addPhysicalGroup(1, [upperSpline, lowerSpline])  
    gmsh.model.setPhysicalName(1, airfoil_tag, "airfoil")

    farfield_tag = gmsh.model.addPhysicalGroup(1, [a1, a2, a3, a4])  
    gmsh.model.setPhysicalName(1, farfield_tag, "farfield")

    # Detail mesh generation
    airfoil_curves = [upperSpline, lowerSpline]
    dist_field = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(dist_field, "EdgesList", airfoil_curves)

    threshold_field = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(threshold_field, "IField", dist_field)
    gmsh.model.mesh.field.setNumber(threshold_field, "LcMin", LcMin)
    gmsh.model.mesh.field.setNumber(threshold_field, "LcMax", lc)
    gmsh.model.mesh.field.setNumber(threshold_field, "DistMin", LcMinDist)
    gmsh.model.mesh.field.setNumber(threshold_field, "DistMax", r)  # farfield distance

    gmsh.model.mesh.field.setAsBackgroundMesh(threshold_field)

    gmsh.option.setNumber("Mesh.Algorithm", meshAlgorithm)  # 6 = Frontal-Delaunay
    gmsh.option.setNumber("Mesh.RecombinationAlgorithm", recombinationAlgorithm) # 0 = less recombinations , 1 = more recombinations
    gmsh.option.setNumber("Mesh.RecombineAll", recombine) # 0 = no recombination, 1 = recombine

    gmsh.option.setNumber("Mesh.ElementOrder", meshOrder) # mesh order, for now just set to default value of 1 ie no higher order meshin
    gmsh.option.setNumber("Mesh.HighOrderOptimize", 1)

    # Finalize mesh and open GUI
    gmsh.model.geo.synchronize()
    gmsh.model.mesh.generate(2)
    gmsh.fltk.run()