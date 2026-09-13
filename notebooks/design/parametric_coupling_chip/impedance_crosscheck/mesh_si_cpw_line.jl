# Minimal single coplanar-waveguide THROUGH line: one 2-fin lumped port at each
# end. Si substrate. Adapted (heavily stripped) from the Palace source example
# examples/cpw/mesh/mesh.jl :: generate_cpw_lumped_mesh -- one CPW instead of the
# side-by-side pair, two ports instead of four.
#
# A uniform 2-port line lets Z0 and gamma be recovered EXACTLY per frequency
# from (S11, S21) by the classic uniform-line inversion (no fitting) -- far more
# robust than fitting an open stub's coth() response, and it is the same
# measurement as ../driven_sweep/check_line_impedance.py (which got ~40 ohm).
#
# Explicit physical group tags (so si_cpw_line.json is stable):
#   1  air (3D)      2  si (3D)
#   10 PEC (trace + both ground planes)      11 farfield (absorbing box)
#   12 gap surfaces
#   20/21 port-1 fins (x=0,   lower/upper gap)
#   22/23 port-2 fins (x=XL,  lower/upper gap)
#
#   include("mesh_si_cpw_line.jl")
#   generate(; filename="mesh.msh", order=1, refinement=1, length_um=4000.0)

using Gmsh: gmsh

function generate(;
    filename::AbstractString,
    order::Integer            = 1,
    refinement::Integer       = 1,
    trace_width_um::Real      = 20.0,
    gap_width_um::Real        = 11.5,
    ground_width_um::Real     = 250.0,
    substrate_height_um::Real = 525.0,
    length_um::Real           = 4000.0,
    port_len_um::Real         = -1.0,   # fin extent along x; <0 => = gap_width
    air_pad_um::Real          = 400.0,
    verbose::Integer          = 5,
    gui::Bool                 = false,
)
    @assert order > 0 && refinement >= 0
    port_len_um < 0 && (port_len_um = gap_width_um)
    k = gmsh.model.occ

    gmsh.initialize()
    gmsh.option.setNumber("General.Verbosity", verbose)
    "cpwline" in gmsh.model.list() && (gmsh.model.setCurrent("cpwline"); gmsh.model.remove())
    gmsh.model.add("cpwline")

    W, S, GW = trace_width_um, gap_width_um, ground_width_um
    L, H, pad = length_um, substrate_height_um, air_pad_um
    PW = port_len_um

    l_trace    = 1.5 * W * (2.0^-refinement)
    l_farfield = 1.0 * H * (2.0^-refinement)

    y1 = GW; y2 = GW + S; y3 = GW + S + W; y4 = GW + 2S + W; y5 = 2GW + 2S + W

    g1 = k.addRectangle(0.0, 0.0, 0.0, L, GW)      # lower ground
    n1 = k.addRectangle(0.0, y1,  0.0, L, S)       # lower gap
    t1 = k.addRectangle(0.0, y2,  0.0, L, W)       # trace
    n2 = k.addRectangle(0.0, y3,  0.0, L, S)       # upper gap
    g2 = k.addRectangle(0.0, y4,  0.0, L, GW)      # upper ground

    p1a = k.addRectangle(0.0,      y1, 0.0, PW, S)     # port 1, lower fin
    p1b = k.addRectangle(0.0,      y3, 0.0, PW, S)     # port 1, upper fin
    p2a = k.addRectangle(L - PW,   y1, 0.0, PW, S)     # port 2, lower fin
    p2b = k.addRectangle(L - PW,   y3, 0.0, PW, S)     # port 2, upper fin

    # Substrate is flush with the metal pattern in x AND y -- no padding --
    # exactly matching Palace's own cpw_lumped example (examples/cpw/mesh/mesh.jl
    # :: generate_cpw_lumped_mesh: `addBox(0, 0, -H, L, dy, H)`, no padding
    # anywhere). Ports also sit flush at x=0/x=L, same as that reference.
    #
    # What that reference DOES pad is the surrounding AIR domain, in all three
    # directions including x: `addBox(-0.5*sep_dy, -sep_dy, -sep_dz, L+sep_dy,
    # dy+2*sep_dy, 2*sep_dz)` extends 0.5*sep_dy beyond each port face in x, not
    # just in y/z. Our domain previously did NOT pad x at all (line ended flush
    # with the domain, so the lumped port sat directly on the absorbing
    # boundary with no open space beyond it to radiate/decay into) -- that is
    # the real deviation from the reference, now fixed by padding x here too.
    substrate = k.addBox(0.0, 0.0, -H, L, y5, H)
    domain    = k.addBox(-pad, -pad, -H - pad, L + 2pad, y5 + 2pad, H + 2pad)
    _, dbound = k.getSurfaceLoops(domain)
    @assert length(dbound) == 1
    dbound = first(dbound)

    geom = filter(x -> x[1] == 2 || x[1] == 3, k.getEntities())
    _, gmap = k.fragment(geom, [])
    k.synchronize()

    si_dom  = first(last.(gmap[findfirst(x -> x == (3, substrate), geom)]))
    air_dom = last.(gmap[findfirst(x -> x == (3, domain), geom)])
    filter!(x -> x != si_dom, air_dom)
    @assert length(air_dom) == 1
    air_dom = first(air_dom)
    gmsh.model.addPhysicalGroup(3, [air_dom], 1, "air")
    gmsh.model.addPhysicalGroup(3, [si_dom], 2, "si")

    flat(ids) = last.(collect(Iterators.flatten(gmap[findall(x -> x[1] == 2 && x[2] in ids, geom)])))
    pec      = flat([g1, g2, t1])
    farfield = flat(dbound)
    gaps     = flat([n1, n2])
    f1a = last.(gmap[findfirst(x -> x == (2, p1a), geom)])
    f1b = last.(gmap[findfirst(x -> x == (2, p1b), geom)])
    f2a = last.(gmap[findfirst(x -> x == (2, p2a), geom)])
    f2b = last.(gmap[findfirst(x -> x == (2, p2b), geom)])
    setdiff!(gaps, f1a, f1b, f2a, f2b)
    setdiff!(farfield, pec, gaps, f1a, f1b, f2a, f2b)

    gmsh.model.addPhysicalGroup(2, pec, 10, "pec")
    gmsh.model.addPhysicalGroup(2, farfield, 11, "farfield")
    gmsh.model.addPhysicalGroup(2, gaps, 12, "gap")
    gmsh.model.addPhysicalGroup(2, f1a, 20, "port1a")
    gmsh.model.addPhysicalGroup(2, f1b, 21, "port1b")
    gmsh.model.addPhysicalGroup(2, f2a, 22, "port2a")
    gmsh.model.addPhysicalGroup(2, f2b, 23, "port2b")

    gmsh.option.setNumber("Mesh.MeshSizeMin", l_trace)
    gmsh.option.setNumber("Mesh.MeshSizeMax", l_farfield)
    gmsh.option.setNumber("Mesh.MeshSizeFromPoints", 0)
    gmsh.option.setNumber("Mesh.MeshSizeFromCurvature", 0)
    gmsh.option.setNumber("Mesh.MeshSizeExtendFromBoundary", 0)

    gpts = last.(filter(x -> x[1] == 0, gmsh.model.getBoundary([(2, z) for z in gaps], false, true, true)))
    gcrv = last.(filter(x -> x[1] == 1, gmsh.model.getBoundary([(2, z) for z in gaps], false, false, false)))
    gmsh.model.mesh.field.add("Distance", 1)
    gmsh.model.mesh.field.setNumbers(1, "PointsList", gpts)
    gmsh.model.mesh.field.setNumbers(1, "CurvesList", gcrv)
    gmsh.model.mesh.field.setNumbers(1, "SurfacesList", gaps)
    gmsh.model.mesh.field.setNumber(1, "Sampling", ceil(L / l_trace))
    gmsh.model.mesh.field.add("Threshold", 2)
    gmsh.model.mesh.field.setNumber(2, "InField", 1)
    gmsh.model.mesh.field.setNumber(2, "SizeMin", l_trace)
    gmsh.model.mesh.field.setNumber(2, "SizeMax", l_farfield)
    gmsh.model.mesh.field.setNumber(2, "DistMin", W)
    gmsh.model.mesh.field.setNumber(2, "DistMax", 0.9 * (H + 2pad))
    gmsh.model.mesh.field.add("Min", 101)
    gmsh.model.mesh.field.setNumbers(101, "FieldsList", [2])
    gmsh.model.mesh.field.setAsBackgroundMesh(101)

    gmsh.option.setNumber("Mesh.Algorithm", 6)
    gmsh.option.setNumber("Mesh.Algorithm3D", 10)
    for tag in gaps
        gmsh.model.mesh.setAlgorithm(2, tag, 8)
    end

    gmsh.model.mesh.generate(3)
    gmsh.model.mesh.setOrder(order)
    gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
    gmsh.option.setNumber("Mesh.Binary", 1)
    gmsh.write(filename)

    if verbose > 0
        n_el = length(gmsh.model.mesh.getElements(3)[2][1])
        println("\n[mesh] $filename  order=$order refinement=$refinement  ~$(n_el) tets")
        println("[mesh] l_trace=$(round(l_trace, digits=2))um l_far=$(round(l_farfield, digits=1))um "
                * "len=$(L)um air_pad=$(pad)um")
        println("[mesh] tags: air=1 si=2 pec=10 farfield=11 gap=12  port1=20/21 port2=22/23\n")
    end
    gui && gmsh.fltk.run()
    gmsh.finalize()
    return filename
end
