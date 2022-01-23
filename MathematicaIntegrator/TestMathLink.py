from MathematicaIntegrator.Constants import MathLink

mathlink = MathLink()

print(mathlink.Call("""f[x_, y_, z_]:=1 / ((0.1 * x + 0.2 * y + 0.3 * z + 0.6 * x * x - 0.5 * y * y + 0.4 * x * y + 0.3 * z * z - 0.2 * x * z + 0.1 * y * z + 0.1 * x * x * x - 0.2 * y * y * y + 0.3 * z * z * z + 0.1 * y * y * x + 0.2 * y * y * z - 0.3 * y * x * x + 0.4 * y * z * z - 0.5) ^ 2);
g[x_, y_, z_]:=0.5 * 0.5 * 0.5 *  f[0.5 * (x + 1) + 0, 0.5 * (y + 1) + 0, 0.5 * (z + 1) + 0];"""))

[bDone, realPart, imagePart] = mathlink.Call("Quiet[If[res = Check[NIntegrate[g[x, y, z], {x, (-1-1 I), (1+1 I)}, {y, (-1-1 I), (1+1 I)}, {z, (-1-1 I), (1+1 I)}], False, {NIntegrate::slwcon, NIntegrate::ncvb}]; NumberQ[res], {True, Re[res], Im[res]}, {False, 0, 0}], {NIntegrate::slwcon, NIntegrate::ncvb}]")
print(bDone)
print(realPart)
print(imagePart)
mathlink.Quit()
