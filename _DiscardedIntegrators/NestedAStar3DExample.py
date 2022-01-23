from Contour3D.Integrand3D import Integrand3D
from Contour3D.Integrator3D import SparseGridIntegrator3D
from _DiscardedIntegrators.NestedAStar3D import NestedAStar3D


def intfunc1(x: complex, y: complex, z: complex) -> complex:
    return 1 / ((0.2 + 0.1 * x + 0.2 * y + 0.3 * z
                 + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
                 + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
                 + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5
                 - 0.1 * x * x * x * x + 0.2 * y * y * y * y - 0.3 * z * z * z * z
                 + 0.3 * x * x * x * y - 0.2 * x * x * x * z + 0.1 * y * y * y * x
                 + 0.1 * y * y * y * z - 0.2 * z * z * z * x + 0.3 * z * z * z * y
                 + 0.1 * x * x * y * z - 0.2 * x * x * y * y + 0.3 * x * x * z * z
                 - 0.3 * y * y * x * z - 0.2 * y * y * z * z + 0.1 * z * z * y * x) ** 2)


intgrand1 = Integrand3D(intfunc1, 0, 1, 0, 1, 0, 1)
intgrand1.SetMathematicaExpress("""1 / ((0.2 + 0.1 * x + 0.2 * y + 0.3 * z
                 + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
                 + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
                 + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5
                 - 0.1 * x * x * x * x + 0.2 * y * y * y * y - 0.3 * z * z * z * z
                 + 0.3 * x * x * x * y - 0.2 * x * x * x * z + 0.1 * y * y * y * x
                 + 0.1 * y * y * y * z - 0.2 * z * z * z * x + 0.3 * z * z * z * y
                 + 0.1 * x * x * y * z - 0.2 * x * x * y * y + 0.3 * x * x * z * z
                 - 0.3 * y * y * x * z - 0.2 * y * y * z * z + 0.1 * z * z * y * x) ^ 2)""")


def intfunc2(x: complex, y: complex, z: complex) -> complex:
    return 1 / ((0.2 + 0.1 * x + 0.2 * y + 0.3 * z
                 + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
                 + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
                 + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5
                 - 0.1 * x * x * x * x + 0.2 * y * y * y * y - 0.3 * z * z * z * z
                 + 0.3 * x * x * x * y - 0.2 * x * x * x * z + 0.1 * y * y * y * x
                 + 0.1 * y * y * y * z - 0.2 * z * z * z * x + 0.3 * z * z * z * y
                 + 0.1 * x * x * y * z - 0.2 * x * x * y * y + 0.3 * x * x * z * z
                 + 0.3 * y * y * x * z - 0.2 * y * y * z * z + 0.1 * z * z * y * x) ** 2)


intgrand2 = Integrand3D(intfunc2, 0, 1, 0, 1, 0, 1)
intgrand2.SetMathematicaExpress("""1 / ((0.2 + 0.1 * x + 0.2 * y + 0.3 * z
                 + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
                 + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
                 + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5
                 - 0.1 * x * x * x * x + 0.2 * y * y * y * y - 0.3 * z * z * z * z
                 + 0.3 * x * x * x * y - 0.2 * x * x * x * z + 0.1 * y * y * y * x
                 + 0.1 * y * y * y * z - 0.2 * z * z * z * x + 0.3 * z * z * z * y
                 + 0.1 * x * x * y * z - 0.2 * x * x * y * y + 0.3 * x * x * z * z
                 + 0.3 * y * y * x * z - 0.2 * y * y * z * z + 0.1 * z * z * y * x) ^ 2)""")

intgrator = SparseGridIntegrator3D()

grid3d = NestedAStar3D(5, 5, 0, intgrator)
print(grid3d.Integrate(intgrand1))
print(grid3d.GatherInfo())

print(grid3d.Integrate(intgrand2))
print(grid3d.GatherInfo())

intgrator.Finish()


"""
(* =========== Copy these to Mathematica ========== *)

Print["Original Integral is Integrate[f[x,y,z], {x,0,1}, {y,0,1}, {z,0,1}]"]
f[x_, y_, z_]:=1 / ((0.2 + 0.1 * x + 0.2 * y + 0.3 * z
                 + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
                 + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
                 + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5
                 - 0.1 * x * x * x * x + 0.2 * y * y * y * y - 0.3 * z * z * z * z
                 + 0.3 * x * x * x * y - 0.2 * x * x * x * z + 0.1 * y * y * y * x
                 + 0.1 * y * y * y * z - 0.2 * z * z * z * x + 0.3 * z * z * z * y
                 + 0.1 * x * x * y * z - 0.2 * x * x * y * y + 0.3 * x * x * z * z
                 - 0.3 * y * y * x * z - 0.2 * y * y * z * z + 0.1 * z * z * y * x) ^ 2);
g[x_, y_, z_]:=0.5 * 0.5 * 0.5 *  f[0.5 * (x + 1) + 0, 0.5 * (y + 1) + 0, 0.5 * (z + 1) + 0];

zr={-1.0,-0.5,-0.5,1.0,1.0};
zi={0.0,0.0,-0.5,-0.5,0.0};
yr={-1.0,1.0};
yi={0.0,0.0};
xr={-1.0,-0.5,-0.5,1.0,1.0};
xi={0.0,0.0,-0.5,-0.5,0.0};
Print["expecting: (-3.5785889419590213+2.0191294496573158j)"]
        res = Sum[NIntegrate[g[x, y, z], 
        {x, xr[[u]] + xi[[u]] I, xr[[u + 1]] + xi[[u + 1]] I}, 
        {y, yr[[v]] + yi[[v]] I, yr[[v + 1]] + yi[[v + 1]] I}, 
        {z, zr[[w]] + zi[[w]] I, zr[[w + 1]] + zi[[w + 1]] I}], 
        {u, 1, Length[xr] - 1}, {v, 1, Length[yr] - 1}, {w, 1, Length[zr] - 1}]
        ListLinePlot[Transpose[{xr, xi}], AxesLabel -> {"Re[x]", "Im[x]"}]
ListLinePlot[Transpose[{yr, yi}], AxesLabel -> {"Re[y]", "Im[y]"}]
ListLinePlot[Transpose[{zr, zi}], AxesLabel -> {"Re[z]", "Im[z]"}]

(* =========== Copy these to Mathematica ========== *)

Print["Original Integral is Integrate[f[x,y,z], {x,0,1}, {y,0,1}, {z,0,1}]"]
f[x_, y_, z_]:=1 / ((0.2 + 0.1 * x + 0.2 * y + 0.3 * z
                 + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z
                 + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z
                 + 0.1 * y * y * x + 0.2 * y * y * z + 0.3 * y * x * x + 0.4 * y * z * z - 0.5
                 - 0.1 * x * x * x * x + 0.2 * y * y * y * y - 0.3 * z * z * z * z
                 + 0.3 * x * x * x * y - 0.2 * x * x * x * z + 0.1 * y * y * y * x
                 + 0.1 * y * y * y * z - 0.2 * z * z * z * x + 0.3 * z * z * z * y
                 + 0.1 * x * x * y * z - 0.2 * x * x * y * y + 0.3 * x * x * z * z
                 + 0.3 * y * y * x * z - 0.2 * y * y * z * z + 0.1 * z * z * y * x) ^ 2);
g[x_, y_, z_]:=0.5 * 0.5 * 0.5 *  f[0.5 * (x + 1) + 0, 0.5 * (y + 1) + 0, 0.5 * (z + 1) + 0];

zr={-1.0,-0.5,-0.5,1.0,1.0};
zi={0.0,0.0,-0.5,-0.5,0.0};
yr={-1.0,1.0};
yi={0.0,0.0};
xr={-1.0,-0.5,-0.5,1.0,1.0};
xi={0.0,0.0,-0.5,-0.5,0.0};
Print["expecting: (-3.4175425109135347+1.8925545860438218j)"]
        res = Sum[NIntegrate[g[x, y, z], 
        {x, xr[[u]] + xi[[u]] I, xr[[u + 1]] + xi[[u + 1]] I}, 
        {y, yr[[v]] + yi[[v]] I, yr[[v + 1]] + yi[[v + 1]] I}, 
        {z, zr[[w]] + zi[[w]] I, zr[[w + 1]] + zi[[w + 1]] I}], 
        {u, 1, Length[xr] - 1}, {v, 1, Length[yr] - 1}, {w, 1, Length[zr] - 1}]
        ListLinePlot[Transpose[{xr, xi}], AxesLabel -> {"Re[x]", "Im[x]"}]
ListLinePlot[Transpose[{yr, yi}], AxesLabel -> {"Re[y]", "Im[y]"}]
ListLinePlot[Transpose[{zr, zi}], AxesLabel -> {"Re[z]", "Im[z]"}]



"""