"""
translate from https://people.math.sc.edu/Burkardt/f_src/patterson_rule_compute/patterson_rule_compute.html
"""
import math


def dot_product(v1: list, v2: list, start: int, end: int) -> float:
    ret = 0.0
    for i in range(start, end + 1):
        ret = ret + v1[i] * v2[i]
    return ret


def idamax(n: int, start: int, dx: list, incx: int) -> int:
    """
    IDAMAX indexes the array element of maximum absolute value.
    Input, integer ( kind = 4 ) N, the number of entries in the vector.
    Input, real ( kind = 8 ) X(*), the vector to be examined.
    Input, integer ( kind = 4 ) INCX, the increment between successive entries of SX.
    Output, integer ( kind = 4 ) IDAMAX, the index of the element of SX of maximum absolute value.
    """
    if n < 1 or incx <= 0:
        return 0
    if 1 == n:
        return 1
    ret: int = 1
    if 1 == incx:
        dmax = abs(dx[start + 0])
        for i in range(1, n):
            if dmax < abs(dx[start + i]):
                ret = i + 1
                dmax = abs(dx[start + i])
    else:
        ix = 1
        dmax = abs(dx[start + 0])
        ix = ix + incx
        for i in range(2, n + 1):
            if dmax < abs(dx[start + ix - 1]):
                ret = i
                dmax = abs(dx[start + ix - 1])
            ix = ix + incx
    return ret


def recura(k: int) -> [float, float, float]:
    """
    This is an example of a user supplied subroutine to define the orthogonal polynomials.
    RECURA ( K, C, D, E ) gives the coefficients C, D and E such that:
    P(K+1,X) = ( C * X + D ) * P(K,X) + E * P(K-1,X)
    Input, integer ( kind = 4 ) K, the index.
    Output, real ( kind = 8 ) C, D, E, the recurrence relation parameters.
    """
    f: float = k + 1
    c: float = (2 * k + 1) / f
    d: float = 0.0
    e: float = - k / f
    return [c, d, e]


def dgefa(a: list, n: int) -> [int, list]:
    """
    DGEFA factors a real general matrix.
    Input/output, real ( kind = 8 ) A(LDA,N).
        On intput, the matrix to be factored.
        On output, an upper triangular matrix and the multipliers used to obtain it.
    The factorization can be written A=L*U, where L is a product of permutation and unit lower triangular matrices,
    and U is upper triangular.
    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
    Input, integer ( kind = 4 ) N, the order of the matrix A.
    Output, integer ( kind = 4 ) IPVT(N), the pivot indices.
    Output, integer ( kind = 4 ) INFO, singularity indicator.
        0, normal value.
        K, if U(K,K) == 0.  This is not an error condition for this subroutine,
        but it does indicate that DGESL or DGEDI will divide by zero if called.
        Use RCOND in DGECO for a reliable indication of singularity.
    """
    # Gaussian elimination with partial pivoting.
    info: int = 0
    ipvt = [0 for _ in range(0, n)]
    for k in range(1, n):
        # Find L = pivot index.
        l: int = idamax(n - k + 1, k - 1, a[k - 1], 1) + k - 1
        ipvt[k - 1] = l
        # Zero pivot implies this column already triangularized.
        if abs(a[l - 1][k - 1]) < 1.0e-22:
            continue
        # Interchange if necessary.
        if l != k:
            t = a[l - 1][k - 1]
            a[l - 1][k - 1] = a[k - 1][k - 1]
            a[k - 1][k - 1] = t
        # Compute multipliers.
        t = -1.0 / a[k - 1][k - 1]
        for kk in range(k, n):
            a[kk][k - 1] = a[kk][k - 1] * t
        # Row elimination with column indexing.
        for j in range(k + 1, n + 1):
            t = a[l - 1][j - 1]
            if l != k:
                a[l - 1][j - 1] = a[k - 1][j - 1]
                a[k - 1][j - 1] = t
            for kk in range(k, n):
                a[kk][j - 1] = a[kk][j - 1] + t * a[kk][j - 1]
    ipvt[n - 1] = n
    if abs(a[n - 1][n - 1]) < 1.0e-22:
        info = n
    return [info, ipvt]


def lfact(gamma: list, n: int, xnode: float, recur) -> [list, float, float]:
    """
    LFACT removes a linear factor from a polynomial expansion.

    This function removes the linear factor (X-XNODE) from the polynomial expansion:
    SUM ( 0 <= I <= N ) GAMMA(I) * P(I,X)
    to give the quotient:
    SUM ( 0 <= I <= N-1 ) DELTA(I) * P(I,X).
    and the remainder and its derivative at XNODE.

    Input, real ( kind = 8 ) GAMMA(0:N), the polynomial from which the factor is to be removed and expressed as:
    GAMMA = SUM (0 <= I <= N) GAMMA(I) * P(I,X)
    Output, real ( kind = 8 ) DELTA(0:N-1), the quotient polynomial expressed as:
    DELTA = SUM (0 <= I <= N-1) DELTA(I) * P(I,X)
    Input, integer ( kind = 4 ) N, the degree of GAMMA.
    Input, real ( kind = 8 ) XNODE, the node to be removed.
    Input, external RECUR ( ), the function which defines the orthogonal polynomials.  See EXTEND for a full description.
    Output, real ( kind = 8 ) R, the remainder from division.
    Output, real ( kind = 8 ) DR, the derivative of R with respect to XNODE.
    (-R/DR is the Newton correction).

    Note: delta is 0:n-1, gamma is 0:n
    """
    r: float = 0
    dr: float = 0
    delta = [0.0 for _ in range(0, n)]
    [ck, dk, ek] = recur(n)
    [ck1, dk1, ek1] = recur(n + 1)
    bk1 = 0
    bk2 = 0
    dbk1 = 0
    dbk2 = 0
    for k in range(n, -1, -1):
        r = gamma[k] + (dk + xnode * ck) * bk1 + ek1 * bk2
        dr = (dk + xnode * ck) * dbk1 + ek1 * dbk2 + ck * bk1
        bk2 = bk1
        bk1 = r
        dbk2 = dbk1
        dbk1 = dr
        ckm1 = 0
        dkm1 = 0
        ekm1 = 0
        if k != 0:
            [ckm1, dkm1, ekm1] = recur(k - 1)
            delta[k - 1] = r * ckm1
        ek1 = ek
        ck = ckm1
        dk = dkm1
        ek = ekm1
    return [delta, r, dr]


def newton(t: list, n: int, xnode: float, recur, idigit: int = 5) -> [int, list, float, float]:
    """
    NEWTON applies Newton's method for a single root of a polynomial.
    This function applies Newton's method to find a single root of the polynomial T:
    T = SUM (0 <= I <= N) T(I) * P(I,X)
    where P(I,X) are the orthogonal polymonials whose recurrence relation is defined by RECUR.
    The value of T is found from the remainder when T is divided by (X-XNODE).
    The derivative (of the remainder) is calculated simultaneously.
    The deflated polynomial:
    DELTA = SUM (I = 0 to N-1) DELTA(I) * P(I,X)
    is also computed.

    Input, real ( kind = 8 ) T(0:N), the polynomial whose roots define the nodes of the quadrature rule and expressed as:
    T = SUM (0 <= I <= N) T(I) * P(I,X)

    Input, integer ( kind = 4 ) N, the degree of the expansion of T.

    Input/output, real ( kind = 8 ) XNODE,
        On input, a rough estimate for the root.
        On output, this estimate has been improved using Newton's method.

    Input, external RECUR ( ), the function which defines the orthogonal polynomials.  See EXTEND for a full description.
    Input, integer ( kind = 4 ) IDIGIT, the node convergence parameter, an integer greater than 0.
        An attempt is made to calculate the nodes to the maximum accuracy possible by the machine precision available.
        IDIGIT controls the assessment procedure to take account of round-off errors
        and specifies the number of least significan decimal digits that can be ignored, that is, attributed to round-off, in the computed relative error.
        A typical value is 5.
    Output, real ( kind = 8 ) DELTA(0:N-1), the coefficients of the deflated polynomial.
    Output, real ( kind = 8 ) ERRVAL, the value of the correction.
        May be used as a measure of the root accuracy when convergence is not achieved.
    Output, integer ( kind = 4 ) IFAIL, error flag.
        * 0, Convergence OK.
        * 1, Unsatisfactory convergence after 50 iterations.

    Note: delta is 0:n-1, t is 0:n
    """
    delta = [0.0 for _ in range(0, n)]
    itera: int = 50
    tol: float = 10.0 ** (-idigit)
    while True:
        itera = itera - 1
        if itera < 0:
            return [1, delta, xnode, 0.0]
        [delta, r, dr] = lfact(t, n, xnode, recur)
        eps = - r / dr
        xnode = xnode + eps
        if tol * abs(eps) < 1.0e-22:
            break
    # Final iteration
    [delta, r, dr] = lfact(t, n, xnode, recur)
    eps = - r / dr
    xnode = xnode + eps
    return [0, delta, xnode, abs(r / dr)]


def eprod(n: int, j: int, work: list, ix: list, recur) -> [int, list]:
    """
    EPROD expands a product of two orthogonal polynomials.
    This function calculates the expansion of a product of two orthogonal polynomials:
        P(N,X) * P(J,X) = SUM (I = N-J to N+J ) COEFF(I) * P(I,X)
    where J must not exceed N.
    The orthogonal polynomials are defined by the recurrence relation calculated by RECUR.
    For proper initialization, the function must first be called with J = 0 and the required value of N.
    Subsequent calls must be in the order J = 1,2,,,,,N with the appropriate expansion being generated from previous values and returned in COEFF(*).
    The coefficients of P(N-J,X),...., P(N+J,X) are stored in the array COEFF(1),...,COEFF(2*J+1).
    Input, integer ( kind = 4 ) N, the highest polynomial degree.
        Note that after the initial call with J = 0 the value of N in this argument is ignored.
    Input, integer ( kind = 4 ) J, the current product of P(J,X) with P(N,X) to be calculated.
        Note that this function must be first called with J = 0 and the required largest N.
        Subsequent calls must be in the order J = 1,2,..,N.
    Output, real ( kind = 8 ) COEFF(2*J+1), the coefficients of the expansion.
    Workspace, real ( kind = 8 ) WORK(2*J+1,2).
    The contents of this work area must not be altered between calls by the calling program.
    Input, integer ( kind = 4 ) LW, the leading dimension of WORK.
    Input, external RECUR ( ), the function which defines the orthogonal polynomials.  See EXTEND for a full description.
    Output, integer ( kind = 4 ) IFAIL, error flag.
        * 0, Result OK.
        * 1, J exceeds N.
        * 2, J has not been called sequentially.

    Note: last is stored as ix[2], ss is stored as ix[3]
    """
    coeff = [0.0 for _ in range(0, 2 * j + 1)]
    # Initialize.
    if 0 == j:
        ix[0] = 1
        ix[1] = 2
        ix[2] = 0
        ix[3] = n
        coeff[0] = 1.0
        work[0][1] = 1.0
        return [0, coeff]
    s = ix[3]
    if s < j:
        print("EPROD - Fatal error!, S < J")
        return [1, coeff]
    if ix[2] != j - 1:
        print("EPROD - Fatal error!, J was not used sequentially.")
        return [2, coeff]
    ix[2] = j
    j2 = j + j
    [cj1, dj1, ej1] = recur(j - 1)
    if 1 == j:
        for k in range(0, j2 + 1):
            coeff[k] = 0
    else:
        for i in range(1, j2 - 2):
            coeff[i + 1] = work[i - 1][ix[0] - 1] * ej1
        coeff[0] = 0
        coeff[1] = 0
        coeff[j2 - 1] = 0.0
        coeff[j2] = 0.0
    ibot = s - j + 1
    itop = s + j - 1
    for ii in range(ibot, itop + 1):
        i = ii - ibot + 1
        [ci, di, ei] = recur(ii)
        coeff[i + 1] = coeff[i + 1] + (work[i - 1][ix[1] - 1] / ci) * cj1
        coeff[i] = coeff[i] + work[i - 1][ix[1] - 1] * (dj1 - (cj1 / ci) * di)
        coeff[i - 1] = coeff[i - 1] - (work[i - 1][ix[1] - 1] / ci) * cj1 * ei
    ii = ix[0]
    ix[0] = ix[1]
    ix[1] = ii
    for i in range(1, j2 + 2):
        work[i - 1][ix[1] - 1] = coeff[i - 1]
    return [0, coeff]


def qfact(n: int, gamma: list, recur, alpha: float, beta: float) -> [list, float, float, float, float, float, float]:
    """
    QFACT divides a polynomial by a quadratic factor.
    This function divides the polynomial:
        sum ( 0 <= I <= N ) GAMMA(I) * P(I,X)
    by the quadratic factor:
        P(2,X) - ALPHA * P(1,X) - BETA * P(0,X)
    giving the quotient:
        sum ( 0 <= I <= N-2 ) DELTA(I) * P(I,X)
    and remainder:
        A * P(1,X) + B * P(0,X)
    where P(I,X) is the orthogonal polynomial of degree I defined by RECUR.

    Input, integer ( kind = 4 ) N, the degree of GAMMA.
    Input, real ( kind = 8 ) GAMMA(0:N), a polynomial to be divided by a quadratic factor.

    Output, real ( kind = 8 ) DELTA(0:N-2), the quotient polynomial of degree N-2.

    Input, external RECUR(), the function which defines the orthogonal polynomials.  See EXTEND for a full description.
    Input, real ( kind = 8 ) ALPHA, BETA, the coefficients of the quadratic factor.

    Output, real ( kind = 8 ) A, B, the remainder coefficients.
    Output, real ( kind = 8 ) AA, the partial of A with respect to ALPHA.
    Output, real ( kind = 8 ) AB, the partial of A with respect to BETA.
    Output, real ( kind = 8 ) BA, the partial of B with respect to ALPHA.
    Output, real ( kind = 8 ) BB, the partial of B with respect to BETA.

    Note: gamma is 0:n and delta is 0:n-2
    """
    delta = [0.0 for _ in range(0, n - 1)]
    a = 0.0
    b = 0.0
    aa = 0.0
    ab = 0.0
    ba = 0.0
    bb = 0.0
    # Initialize coefficients.
    dnp2 = 0.0
    dnp1 = 0.0
    dn = 0.0
    dnm1 = 0.0
    # Partial coefficients wrt ALPHA.
    adnp2 = 0.0
    adnp1 = 0.0
    adn = 0.0
    adnm1 = 0.0
    # Partial coefficients wrt BETA.
    bdnp2 = 0.0
    bdnp1 = 0.0
    bdn = 0.0
    bdnm1 = 0.0
    # Scaling parameters.
    sn1 = 1.0
    sn2 = 1.0
    sn3 = 1.0
    sn4 = 1.0
    [c0, d0, e0] = recur(0)
    [c1, d1, e1] = recur(1)
    [c2, d2, e2] = recur(2)
    [c3, d3, e3] = recur(3)
    r0 = - c0 * e1 / c1
    r1 = - c0 * e2 / c2
    r2 = - c0 * e3 / c3
    vm1 = d0 - c0 * d1 / c1
    vm2 = d0 - c0 * d2 / c2
    w0 = - r1 * e1
    w1 = - c1 * r2 * e2 / c2
    v1 = d1 * r1 - c1 * vm2 * e2 / c2 - c1 * r1 * d1 / c1
    k = n - 2
    [ck4, dk4, ek4] = recur(k + 4)
    [ck3, dk3, ek3] = recur(k + 3)
    [ck2, dk2, ek2] = recur(k + 2)
    [ck1, dk1, ek1] = recur(k + 1)
    vlk4 = c0 / ck3
    vlk3 = c0 / ck2
    vlk2 = c0 / ck1
    rk3 = - c0 * ek4 / ck4
    rk2 = - c0 * ek3 / ck3
    vmk3 = d0 - dk3 * vlk4
    vmk2 = d0 - dk2 * vlk3
    # Extract quadratic factor and find partial derivatives
    for k in range(n - 2, -1, -1):
        [ck, dk, ek] = recur(k)
        vlk1 = c0 / ck
        rk1 = - c0 * ek2 / ck2
        vmk1 = d0 - dk1 * vlk2
        sk2 = c1 * vlk1 * vlk2 / c0
        tk2 = vlk2 * (d1 - c1 * dk2 / ck2) + c1 * vmk1 / ck1
        uk2 = d1 * vmk2 + e1 - c1 * vlk3 * ek3 / ck3 - c1 * vmk2 * dk2 / ck2 + c1 * rk1 / ck1
        vk2 = d1 * rk2 - c1 * vmk3 * ek3 / ck3 - c1 * rk2 * dk2 / ck2
        wk2 = - c1 * rk3 * ek3 / ck3
        cf1 = (alpha * vlk2 - tk2) / sn1
        cf2 = (beta + alpha * vmk2 - uk2) / sn2
        cf3 = (alpha * rk2 - vk2) / sn3
        cf4 = - wk2 / sn4
        rs = gamma[k + 2]
        d = rs + cf1 * dnm1 + cf2 * dn + cf3 * dnp1 + cf4 * dnp2
        delta[k] = d / sk2
        da = vlk2 * dnm1 / sn1 + vmk2 * dn / sn2 + rk2 * dnp1 / sn3 + cf1 * adnm1 + cf2 * adn + cf3 * adnp1 + cf4 * adnp2
        db = dn / sn2 + cf1 * bdnm1 + cf2 * bdn + cf3 * bdnp1 + cf4 * bdnp2
        # Recycle old values.
        sn4 = sn3
        sn3 = sn2
        sn2 = sn1
        sn1 = sk2
        dnp2 = dnp1
        dnp1 = dn
        dn = dnm1
        dnm1 = d
        adnp2 = adnp1
        adnp1 = adn
        adn = adnm1
        adnm1 = da
        bdnp2 = bdnp1
        bdnp1 = bdn
        bdn = bdnm1
        bdnm1 = db
        ck4 = ck3
        ck3 = ck2
        ck2 = ck1
        ck1 = ck
        dk4 = dk3
        dk3 = dk2
        dk2 = dk1
        dk1 = dk
        ek4 = ek3
        ek3 = ek2
        ek2 = ek1
        ek1 = ek
        vlk4 = vlk3
        vlk3 = vlk2
        vlk2 = vlk1
        rk3 = rk2
        rk2 = rk1
        vmk3 = vmk2
        vmk2 = vmk1
    cf1 = alpha
    cf2 = beta + alpha * vm1 - r1
    cf3 = alpha * r1 - v1
    cf4 = - w1
    cf5 = alpha * r0
    rs0 = gamma[0]
    rs1 = gamma[1]
    dnm1 = dnm1 / sn1
    dn = dn / sn2
    dnp1 = dnp1 / sn3
    dnp2 = dnp2 / sn4
    adnm1 = adnm1 / sn1
    adn = adn / sn2
    adnp1 = adnp1 / sn3
    adnp2 = adnp2 / sn4
    bdnm1 = bdnm1 / sn1
    bdn = bdn / sn2
    bdnp1 = bdnp1 / sn3
    bdnp2 = bdnp2 / sn4
    # Remainder.
    a = rs1 + cf1 * dnm1 + cf2 * dn + cf3 * dnp1 + cf4 * dnp2
    b = rs0 + beta * dnm1 + cf5 * dn - w0 * dnp1
    # Partials.
    aa = dnm1 + vm1 * dn + r1 * dnp1 + cf1 * adnm1 + cf2 * adn + cf3 * adnp1 + cf4 * adnp2
    ab = dn + cf1 * bdnm1 + cf2 * bdn + cf3 * bdnp1 + cf4 * bdnp2
    ba = r0 * dn + beta * adnm1 + cf5 * adn - w0 * adnp1
    bb = dnm1 + beta * bdnm1 + cf5 * bdn - w0 * bdnp1
    return [delta, a, b, aa, ab, ba, bb]


def bair(n: int, polin: list, a012: list, recur, idigit: int = 5) -> [int, list, float]:
    """
    BAIR seeks roots of a polynomial.

    This function carries out a generalized Bairstow root extraction for the polynomial:
        SUM(I = 0 to N)  POLIN(I) * P(I,X).
    It calculates the root as a quadratic factor:
        A2 * P(2,X) - A1 * P(1,X) - A0 * P(0,X)
    where P(I,X) is a general orthogonal polynomial of degree I.

    Input, integer ( kind = 4 ) N, the degree of input polynomial POLIN.
    Input, real ( kind = 8 ) POLIN(0:N), coefficients of the polynomial whose quadratic factor is to be found, i.e.
        POLIN = SUM(I = 0 to N) POLIN(I) * P(I,X)
    Output, real ( kind = 8 ) POLOUT(0:N-2), coefficients of the deflated polynomial of degree N-2 with the quadratic factor removed, i.e.
        POLOUT = SUM(I = 0 to N-2) POLOUT(I) * P(I,X)

    Input/output, real ( kind = 8 ) A0, A1, A2,
        On input, the estimated quadratic factors.
        On output, the estimate has been improved.

    Input, external RECUR ( ), the function which defines the orthogonal polynomials.
        See EXTEND for full description.

    Input, integer ( kind = 4 ) IDIGIT, the node convergence parameter, an integer greater than 0.
        An attempt is made to calculate the nodes to the maximum accuracy possible by the machine precision available.
        IDIGIT controls the assessment procedure to take account of round-off errors
        and specifies the number of least signific decimal digits that can be ignored (i.e. attributed to round-off)
        in the computed relative error.  A typical value is 5.

    Output, real ( kind = 8 ) ERRVAL, the mean value of the correction to the coefficients of the quadratic factor.
        May be used as a measure of the root accuracy when convergence is not achieved.

    Output, integer ( kind = 4 ) IFAIL, error flag.
        * 0, Quadratic factor found.
        * 1, Convergence not achieved after 50 iterations.

    Note: polin is 0:n, polout is 0:n-2
    """
    polout = [0 for _ in range(0, n - 1)]
    iter = 100
    # Special cases.
    if 1 == n:
        a012[0] = - polin[0]
        a012[1] = - polin[1]
        a012[2] = 0.0
        return [0, polout, 0.0]
    if 2 == n:
        a012[0] = - polin[0]
        a012[1] = - polin[1]
        a012[2] = polin[2]
        return [0, polout, 0.0]
    tol = 10.0 ** (- max(1, idigit))
    alpha = a012[1] / a012[2]
    beta = a012[0] / a012[2]
    eps = 0
    delta = 0
    while True:
        iter = iter - 1
        if iter < 0:
            print("BAIR - Warning! Iteration did not meet convergence criterion.")
            a012[0] = beta
            a012[1] = alpha
            a012[2] = 1.0
            return [1, polout, 0.5 * (abs(eps) + abs(delta))]
        [polout, a, b, aa, ab, ba, bb] = qfact(n, polin, recur, alpha, beta)
        scale = max(abs(ab), abs(bb))
        f1 = ab / scale
        f2 = bb / scale
        delta = (b * f1 - a * f2) / (aa * f2 - ba * f1)
        scale = max(abs(ba), abs(aa))
        f1 = ba / scale
        f2 = aa / scale
        eps = (a * f1 - b * f2) / (bb * f2 - ab * f1)
        alpha = alpha + delta
        beta = beta + eps
        # Test for convergence. Stop if correction is less than 1/TOL times the smallest machine relative error.
        if abs(delta) <= tol * (abs(alpha) + 1.0) and abs(eps) <= tol * (abs(beta) + 1.0):
            break
        if tol * abs(delta) < 1.0e-22 and tol * abs(eps) < 1.0e-22:
            break
    # Final iteration to tidy up result.
    [polout, a, b, aa, ab, ba, bb] = qfact(n, polin, recur, alpha, beta)
    scale = max(abs(ab), abs(bb))
    f1 = ab / scale
    f2 = bb / scale
    delta = (b * f1 - a * f2) / (aa * f2 - ba * f1)
    scale = max(abs(ba), abs(aa))
    f1 = ba / scale
    f2 = aa / scale
    eps = (a * f1 - b * f2) / (bb * f2 - ab * f1)
    alpha = alpha + delta
    beta = beta + eps
    a012[0] = beta
    a012[1] = alpha
    a012[2] = 1.0
    return [0, polout, 0.5 * (abs(eps) + abs(delta))]


def dgesl(a: list, n: int, ipvt: list, b: list, job: int):
    """
    DGESL solves a real general linear system A * X = B.

    DGESL can solve either of the systems A * X = B or A' * X = B.
    The system matrix must have been factored by DGECO or DGEFA.
    A division by zero will occur if the input factor contains a zero on the diagonal.
    Technically this indicates singularity but it is often caused by improper arguments or improper setting of LDA.
    It will not occur if the subroutines are called correctly and if DGECO has set 0.0 < RCOND or DGEFA has set INFO == 0.

    Input, real ( kind = 8 ) A(LDA,N), the output from DGECO or DGEFA.
    Input, integer ( kind = 4 ) LDA, the leading dimension of A.
    Input, integer ( kind = 4 ) N, the order of the matrix A.
    Input, integer ( kind = 4 ) IPVT(N), the pivot vector from DGECO or DGEFA.
    Input/output, real ( kind = 8 ) B(N).
        On input, the right hand side vector.
        On output, the solution vector.
    Input, integer ( kind = 4 ) JOB.
        0, solve A * X = B;
        nonzero, solve A' * X = B.
    """
    # Solve A * X = B.
    if 0 == job:
        for k in range(1, n - 2):
            l = ipvt[k - 1]
            t = b[l - 1]
            if l != k:
                b[l - 1] = b[k - 1]
                b[k - 1] = t
            for kk in range(k, n):
                b[kk] = b[kk] + t * a[kk][k]
        for k in range(n, 0, -1):
            b[k - 1] = b[k - 1] / a[k - 1][k - 1]
            t = -b[k - 1]
            for kk in range(0, k - 1):
                b[kk] = b[kk] + t * a[kk][k - 1]
    else:
        # Solve A' * X = B.
        for k in range(1, n + 1):
            alist = [a[kk][k - 1] for kk in range(0, k - 1)]
            t = dot_product(alist, b, 0, k - 2)
            b[k - 1] = (b[k - 1] - t) / a[k - 1][k - 1]
        for k in range(n - 1, 0, -1):
            alist = [a[kk][k - 1] for kk in range(k, n)]
            b[k - 1] = b[k - 1] + dot_product(alist, b, k, n - 1)
            l = ipvt[k - 1]
            if l != k:
                t = b[l - 1]
                b[l - 1] = b[k - 1]
                b[k - 1] = t


def gener(t: list, m0n: list, m: int, recur, symmet: bool, worka: list) -> [int, list]:
    """
    GENER calculates the polynomial defining the optimal new nodes.

    Given N preassigned quadrature nodes defined as the roots of the polynomial expansion:
        SUM (I = M0 to N) (TI/HI) * P(I,X)
    calculate the polynomial expansion:
        SUM (I = 0 to M) SI * P(I,X)
    whose roots are the M optimal nodes and the new expansion:
        SUM (I = M to N+M) (RI/HI) * P(I,X)
    whose roots are to the (N + M) nodes of the full extended quadrature rule.

    Input/output, real ( kind = 8 ) T(0:N);
        On input, the coefficients TI of the polynomial whose roots
            define the N preassigned nodes of the quadrature rule and expressed as:
            SUM (I = M0 to N) (TI/HI) * P(I,X)
        where HI is the integral of W(X) * P(I,X)**2 over the interval for which orthogonality
            with respect the weight W(X) is defined (moment integrals) and P(I,X) is the orthogonal polynomial of degree I.
            T(I-M0) holds the value of TI. This array should be declared
            to have at least max(N-M0+1,M+1) elements in the calling program.
        On output, the coefficients of the new orthogonal expansion whose roots are the nodes of the extended quadrature rule
            (that is the preassigned nodes plus the extended nodes).
            It is expressed as:
            SUM (I = M to N+M) (TI/HI) * P(I,X)
            where N and M have their original values. T(I-M) holds the value of TI. See input argument of T for definitions.
    Input/output, integer ( kind = 4 ) M0,
        On input, the lower limit to the expansion of T.
        On output, this is updated to correspond with the output value of T.
    Input/output, integer ( kind = 4 ) N,
        On input, the upper limit to expansion of T.
        On output, this is updated to correspond with the output value of T.
    Input, integer ( kind = 4 ) M, the number of nodes to be optimally added.
    Input, external RECUR ( ), the user supplied function which defines the orthogonal polynomials.
        Given K, CALL RECUR(K,C,D,E) gives the coefficients C,D and E such that,
        P(K+1,X) = ( C * X + D ) * P(K,X) + E * P(K-1,X)
        The parameters are defined as follows:
            K = Index
            C, D, E = parameters in the recurrence relation (functions of K).
    Input, SYMMET
        * FALSE, if no advantage is to be taken of symmetry, if any,
            about x = 0 in the interval of integration and the orthogonality  weight function.
            Note that if symmetry in fact does exist setting this parameter to zero
            will still produce correct results - only efficiency is effected.
        * TRUE, if the extended rule computations should exploit symmetry about x = 0 in the interval of integration
        and the orthogonality  weight function.
        This reduces the size of the system of linear equations determining EXT by a factor of about 2 (see WORKA).
        If symmetry does not in fact exist erroneous results will be produced.

    Output, real ( kind = 8 ) EXT(0:M),
        the coefficients of the polynomial whose roots are the  new extended nodes and expressed as:
        EXT = SUM (I = 0 to M) EXT(I) * P(I,X)
    Workspace, integer ( kind = 4 ) IWORK(M).
    Workspace, real ( kind = 8 ) WORKA(M+1,max(M+1,N+1)).
        If SYMMET = TRUE, the dimension can be reduced to M/2+1 by max(M/2+1,N/2+1).
    Input, integer ( kind = 4 ) LDA, the number of elements in the leading dimension of WORKA declared in the calling program.
    Input, real ( kind = 8 ) WORKB(2*M+1,3).
    Input, integer ( kind = 4 ) LDB, the number of elements in the leading dimension of WORKB declared in the calling program.
    Output, integer ( kind = 4 ) IFLAG, error flag.
        * 0, No error detected
        * 1, The linear system of equations defining the polynomial
        whose roots are the extended nodes became singular or very ill-conditioned.

    Note: T, EXT from 0 to N/M
    IWORK(M), WORKA(M+1,max(M+1,N+1)), WORKB(2*M+1,3) ALSO start from 0
    """
    ext = [0.0 for _ in range(0, m + 1)]
    # Look for trivial case.
    if 0 == m0n[1]:
        ext[m] = 1.0
        t[0] = 1.0
        m0n[0] = m
        m0n[1] = m
        return [0, ext]
    # General case.
    neven = (0 == (m0n[1] % 2))
    nm = m0n[1] + m
    # Form matrix.
    for s in range(0, m + 1):
        msodd = (1 == ((m + s) % 2))
        if neven and msodd and symmet:
            continue
        ixForEprod = [0, 0, 0, 0]
        workForEprod = [[0.0, 0.0] for _ in range(0, 2 * m + 1)]
        for j in range(0, s + 1):
            [_, res] = eprod(s, j, workForEprod, ixForEprod, recur)
            if ((m0n[1] + s + j) % 2) != 1 or (not symmet):
                iref = s - j
                itop = min(m0n[1], j + s)
                ibot = max(m0n[0], iref)
                total = 0.0
                for i in range(ibot, itop + 1):
                    total = total + t[i - m0n[0]] * res[i - iref]
                if not symmet:
                    worka[s][j] = total
                    worka[j][s] = total
                else:
                    if neven:
                        worka[s // 2][j // 2] = total
                        worka[j // 2][s // 2] = total
                    elif msodd:
                        worka[s // 2][j // 2] = total
                    else:
                        worka[j // 2][s // 2] = total
    neq: int = 0
    if symmet:
        neq = m // 2
    else:
        neq = m
    # Solve for expansion coefficients.
    if 0 < neq:
        [info, ipvt] = dgefa(worka, neq)
        if info != 0:
            return [1, ext]
        # bcopy = [worka[0][ib] for ib in range(neq, neq + neq)]
        dgesl(worka, neq, ipvt, worka[0], 0)
        # for ii in range(0, neq):
        #     worka[0][ii + neq] = bcopy[ii]
        # Store expansion coefficients.
        for ii in range(0, neq):
            ext[ii] = - worka[ii][neq]
        ext[neq] = 1.0
    # Calculate new T polynomial.
    if not symmet:
        # Non-symmetric case.
        for s in range(m, nm + 1):
            if s != m:
                ixForEprod = [0, 0, 0, 0]
                workForEprod = [[0.0, 0.0] for _ in range(0, 2 * m + 1)]
                for j in range(0, m + 1):
                    [_, res] = eprod(s, j, workForEprod, ixForEprod, recur)
                    iref = s - j
                    itop = min(m0n[1], j + s)
                    ibot = max(m0n[0], iref)
                    total = 0.0
                    for i in range(ibot, itop + 1):
                        # ir = i - iref
                        total = total + t[i - m0n[0]] * res[i - iref]
                    worka[m][j] = total
            dotWorka = [worka[m][ii] for ii in range(0, m + 1)]
            worka[m - 1][s - m] = dot_product(dotWorka, ext, 0, m)
        # Overwrite old values of T.
        for i in range(0, m0n[1] + 1):
            t[i] = worka[m - 1][i]
    else:
        # Symmetric case.
        for s in range(m, nm + 1):
            if ((m + m0n[1] + s) % 2) != 1:
                ixForEprod = [0, 0, 0, 0]
                workForEprod = [[0.0, 0.0] for _ in range(0, 2 * m + 1)]
                for j in range(0, m + 1):
                    [_, res] = eprod(s, j, workForEprod, ixForEprod, recur)
                    if ((m0n[1] + s + j) % 2) != 1:
                        iref = s - j
                        itop = min(m0n[1], j + s)
                        ibot = max(m0n[0], iref)
                        total = 0.0
                        for i in range(ibot, itop + 1):
                            # ir = i - iref
                            total = total + t[i - m0n[0]] * res[i - iref]
                        worka[neq][j // 2] = total
            dotWorka = [worka[neq][ii] for ii in range(0, neq + 1)]
            worka[neq - 1][(s - m) // 2] = dot_product(dotWorka, ext, 0, neq)
        # Overwrite old values of T in full unsymmetric form.
        ic = m0n[1] // 2
        miss = True
        for j in range(m0n[1], -1, -1):
            miss = not miss
            if miss:
                t[j] = 0.0
            else:
                if 0 <= (neq - 1):
                    t[j] = worka[neq - 1][ic]
                    ic = ic - 1
        # Convert EXT to full unsymmetric form.
        extCopy = [ipvt[j] if j < len(ipvt) else 0 for j in range(0, len(ext))]
        extCopy[m] = 1.0
        ic = neq - 1
        miss = False
        for j in range(m - 1, -1, -1):
            miss = not miss
            if miss:
                extCopy[j] = 0
            else:
                extCopy[j] = ext[ic]
                ic = ic - 1
        ext = extCopy.copy()
    # Scale the new T polynomial.
    pmax = -1
    for i in range(0, m0n[1] + 1):
        if pmax < abs(t[i]):
            pmax = abs(t[i])
    for i in range(0, m0n[1] + 1):
        t[i] = t[i] / pmax
    m0n[1] = nm
    m0n[0] = m
    return [0, ext]


def rsort(a: list, n: int, iflag: int):
    """
    RSORT carries out a simple ripple sort.
    Input/output, real ( kind = 8 ) A(N), the array to be sorted.
    Input, integer ( kind = 4 ) N, the number of elements to be sorted
    Input, integer ( kind = 4 ) IFLAG, determines the sort direction.
        * 0, for ascending sort;
        * 1, for descending sort.
    """
    ascend: bool = (iflag == 0)
    # Begin scans.
    for j in range(n - 1, 0, -1):
        done: bool = True
        for k in range(1, j + 1):
            k1 = k
            k2 = k
            if ascend:
                k2 = k + 1
            else:
                k1 = k + 1
            # Exchange elements.
            if a[k2 - 1] < a[k1 - 1]:
                val = a[k1 - 1]
                a[k1 - 1] = a[k2 - 1]
                a[k2 - 1] = val
                done = False
        if done:
            return


def transf(t: list, m: int, n: int, recur, iflag: int):
    """
    TRANSF scales a polynomial expansion with respect to the moments.
    This function scales the polynomial expansion:
        SUM (M to N) TI * P(I,X)
    with respect to the moments HI of the orthogonality weight function
    giving the expansion:
        H0 * SUM (M to N) (TI/HI) * P(I,X)
    or
        (1/H0) * SUM (M to N) (TI*HI) * P(I,X)
    depending on the value of IFLAG.

    Input/output, real ( kind = 8 ) T(0:N), the coefficients TI of the polynomial expansion to be scaled and expressed as:
        SUM (I = M to N) TI * P(I,X)
        T(I-M) holds the value of TI.  On output, the polynomial has been scaled.
    Input, integer ( kind = 4 ) M, the lower limit to the expansion of T.
    Input, integer ( kind = 4 ) N, the upper limit to the expansion of T.
    Input, external RECUR ( ), the function which defines the orthogonal polynomials.  See EXTEND for a full description.
    Input, integer ( kind = 4 ) IFLAG, the operation to be carried out.
        * 0, if coefficient TI is to be replaced by TI*(H0/HI).
        * 1, if coefficient TI is to be replaced by TI*(HI/H0).
    """
    h = 1.0
    ckm1 = 0
    for k in range(0, n + 1):
        [ck, dk, ek] = recur(k)
        if k != 0:
            h = - ckm1 / ck * ek * h
        if m <= k:
            if 0 == iflag:
                t[k - m] = t[k - m] / h
            else:
                t[k - m] = t[k - m] * h
        ckm1 = ck


def assign(n: int, pnodes: list, recur) -> [int, list]:
    """
    ASSIGN generates the polynomial whose roots are the preassigned nodes.

    This routine generates the initial polynomial T whose roots are the required preassigned nodes.
    It requires a user-supplied routine RECUR(K,C,D,E) which defines the recursion coefficients of the orthogonal polynomials.
    The routine takes as input the polynomial index K, and returns coefficients C, D and E such that:
        P(K+1,X) = ( C * X + D ) * P(K,X) + E * P(K-1,X)

    Input, integer ( kind = 4 ) N, the number of preassigned nodes.
    Input, real ( kind = 8 ) PNODES(N), the preassigned nodes.
    Input, external RECUR, the user-supplied function which defines the orthogonal polynomials.
    Output, real ( kind = 8 ) T(0:N), the coefficients of the polynomial
    whose roots define the preassigned nodes of the quadrature rule and expressed as:
        H0 * Sum ( 0 <= I <= N ) T(I)/HI * P(I,X)
    Output, integer ( kind = 4 ) IERR.
        * 0, No error detected
        * 1, The linear system of equations used to generate the polynomial T became singular or very ill-conditioned.
    """
    t = [0.0 for _ in range(0, n + 1)]
    a = [[0.0 for _ in range(0, n)] for _ in range(0, n)]
    b = [0.0 for _ in range(0, n)]
    # Set up the linear system of equations.
    for i in range(1, n + 1):
        x = pnodes[i - 1]
        p0 = 0.0
        p1 = 1.0
        p = p1
        for j in range(1, n + 1):
            a[i - 1][j - 1] = p
            [c0, d0, e0] = recur(j)
            p = (c0 * x + d0) * p1 + e0 * p0
            p0 = p1
            p1 = p
        b[i - 1] = p
    # Solve the linear system.
    [info, ipvt] = dgefa(a, n)
    if info != 0:
        print("ASSIGN - Fatal error!")
        print("The linear system is too ill-conditioned.")
        return [1, t]
    dgesl(a, n, ipvt, b, 0)
    # Set T.
    for i in range(0, n):
        t[i] = -b[i]
    t[n] = 1.0
    # Weight with moments.
    transf(t, 0, n, recur, 1)


def roots(a0: float, a1: float, a2: float, recur) -> [int, complex, complex]:
    """
    ROOTS calculates roots of a quadratic factor.

    This function calculates the roots corresponding to the quadratic factor
        A2 * P(2,X) - A1 * P(1,X) - A0 * P(0,X)
    where P(I,X) is a general orthogonal polynomial of degree I defined by the recurrence calculated by RECUR.

    Input, real ( kind = 8 ) A0, A1, A2, the coefficients of the quadratic factor.
    Output, real ( kind = 8 ) ZREAL1, ZIMAG1, the real and imaginary parts of root 1.
    Output, real ( kind = 8 ) ZREAL2, ZIMAG2, the real and Imaginary parts of root 2.
    Input, external RECUR ( ), the function which defines the orthogonal polynomials.  See EXTEND for full description.
    Output, integer ( kind = 4 ) INFO, error flag.
        * 0, two roots found.
        * 1, one root only (A2 = 0).
    """
    [c0, d0, e0] = recur(0)
    if abs(a2) < 1.0e-22:
        if abs(a0 + a1 * d0) < 1.0e-22:
            return [1, 0 + 0j, 0 + 0j]
        return [1, -(a0 + a1 * d0) / a1 / c0 + 0j, 0 + 0j]
    [c1, d1, e1] = recur(1)
    aa = - c0 * c1 * a2
    bb = - a2 * (c0 * d1 + d0 * c1) + c0 * a1
    cc = - d0 * d1 * a2 - e1 * a2 + a0 + a1 * d0
    z = bb * bb - 4.0 * aa * cc
    zr = math.sqrt(abs(z))
    if 0 <= z:
        zreal1 = 0.5 * (- bb - (zr if bb >= 0 else (-zr))) / aa
        return [0, zreal1 + 0j, cc / aa / zreal1 + 0j]
    zreal1 = - 0.5 * bb / aa
    zimag1 = 0.5 * zr / aa
    return [0, zreal1 + zimag1 * 1j + zreal1 - zimag1 * 1j]


def solve(ext: list, m: int, symmet: bool, recur, work1: list, work2: list, idigit: int = 5) -> [int, int, list, list, list, list]:
    """
    SOLVE calculates roots of an orthogonal polynomial expansion.
    This function calculates the roots of the orthogonal polynomial expansion:
        SUM ( 0 <= I <= M ) EXT(I) * P(I,X)

    Input, real ( kind = 8 ) EXT(0:M), the coefficients of the polynomial whose roots are required
        (nodes of the quadrature rule) and expressed as:
        SUM (I = 0 to M) EXT(I) * P(I,X)
    The recurrence relation for the orthogonal polynomials P(I,X) is defined by RECUR.

    Input, integer ( kind = 4 ) M, the upper limit to expansion EXT (polynomial degree).
    Input, logical SYMMET
        * FALSE, if no advantage can be taken of symmetry about x = 0
        in the interval of integration and the orthogonality weight function.
        * TRUE, if symmetry exists about x = 0 in the interval of integration and the orthogonality weight function.
    Input, external RECUR ( ), the user supplied function which defines the orthogonal polynomials.
        Given K, RECUR ( K, C, D, E ) gives the coefficients C,D and E such that,
        P(K+1,X) = ( C * X + D ) * P(K,X) + E * P(K-1,X)
    The parameters are defined as follows:
        K = Index
        C, D, E = parameters in the recurrence relation (functions of K).
    Input, integer ( kind = 4 ) IDIGIT, the node convergence parameter, an integer greater than 0.
        An attempt is made to calculate the nodes to the maximum accuracy possible by the machine precision available.
        IDIGIT controls the assessment procedure to take account of round-off errors
        and specifies the number of least significan decimal digits that can be ignored (i.e. attributed to round-off) in the computed relative error.
        A typical value is 5.
    Output, real ( kind = 8 ) QR(M), the real parts of the roots fo EXT.
    Output, real ( kind = 8 ) QI(M), the imaginary parts of the roots of EXT.
        (Hopefully these values are zero!).
    Output, integer ( kind = 4 ) NODES, the number of extended nodes found.
        Normally equals M but see IERR.
    Output, real ( kind = 8 ) ERR(M), a measure of the relative error in the roots.
        This may be inspected if the convergence error flag has been raised (IERR = 2)
        to decide if the roots in question are acceptable.
    (ERR(*) actually gives the mean last correction to the quadratic factor
        in the generalized Bairstow root finder (see BAIR).
    Output, integer ( kind = 4 ) ICHECK(M), root convergence flags.
        Elements 1 to NODES give information on the convergence of the roots of the polynomial EXT.
        * 0, Convergence of I th root satisfactory;
        * 1, Convergence of I th root unsatisfactory.
    Worskpace, real ( kind = 8 ) WORK(M+1,2).
    Input, integer ( kind = 4 ) LDW, the leading dimension of WORK (which must be at least M+1).
    Output, integer ( kind = 4 ) IERR
        * 0, No error detected
        * 1, Poor convergence has been detected in the calculation of the roots of EXT
            or all roots have not been found (M not equal to NODES).  See also ERR.
        * 2, Possible imaginary nodes detected.
    """
    nodes = 0
    nroot = 0
    ierr = 0
    vrt1 = 0.0000001
    icheck = [0 for _ in range(0, m)]
    err = [0.0 for _ in range(0, m)]
    qr = [0.0 for _ in range(0, m)]
    qi = [0.0 for _ in range(0, m)]
    # If M is odd, find and remove initial real root using Newton iteration.
    # Set WORK(*,1) to polynomial to be processed.
    if 1 == (m % 2):
        [ifail, delta, zr1, errval] = newton(ext, m, vrt1, recur, idigit)
        nodes = nodes + 1
        icheck[nodes - 1] = ifail
        err[nodes - 1] = errval
        qr[nodes - 1] = zr1
        qi[nodes - 1] = 0.0
        for ii in range(0, m):
            work1[ii] = delta[ii]
        nroot = m - 1
    else:
        for ii in range(0, m + 1):
            work1[ii] = ext[ii]
        nroot = m
    # Find the remaining root pairs.
    # Calculate seed approximation for quadratic factor.
    if nroot != 0:
        [c0, d0, e0] = recur(0)
        [c1, d1, e1] = recur(1)
        rt1 = vrt1
        if symmet:
            rt2 = - rt1
        else:
            rt2 = 0.0
        p1a = c0 * rt1 + d0
        p1b = c0 * rt2 + d0
        p2a = (c1 * rt1 + d1) * p1a + e1
        p2b = (c1 * rt2 + d1) * p1b + e1
        det = c0 * (rt1 - rt2)
        sa1 = (p2a - p2b) / det
        sa0 = (p1a * p2b - p1b * p2a) / det
        reset = True
        # Alternate approximation which introduces a small complex component.
        rt1 = vrt1
        rt2 = vrt1
        sfa1 = (c0 * d1 + d0 * c1) / c0 + 2.0 * c1 * rt1
        sfa0 = d0 * d1 + e1 - d0 * sfa1 - c0 * c1 * (rt1 * rt1 + rt2 * rt2)
        # IP1 points to the current deflated polynomial.
        ip1 = 1
        ip2 = 2
        works = [work1, work2]
        while True:
            worka = work1 if (1 == ip1) else work2
            workb = work1 if (1 == ip2) else work2
            if reset:
                a2 = 1.0
                a1 = sa1
                a0 = sa0
                reset = False
            a012 = [a0, a1, a2]
            [ifail, plout, errval] = bair(nroot, worka, a012, recur, idigit)
            a0 = a012[0]
            a1 = a012[1]
            a2 = a012[2]
            for kk in range(0, nroot - 1):
                workb[kk] = plout[kk]
            # On failure, try again with complex components introduced.
            if ifail != 0:
                a2 = 1.0
                a1 = sfa1
                a0 = sfa0
                reset = True
                a012 = [a0, a1, a2]
                [ifail, plout, errval] = bair(nroot, works[ip1 - 1], a012, recur, idigit)
                a0 = a012[0]
                a1 = a012[1]
                a2 = a012[2]
                for kk in range(0, nroot - 1):
                    works[ip2 - 1][kk] = plout[kk]
            # Apply Bairstow to full expansion to avoid error accumulation.
            a012 = [a0, a1, a2]
            [ifail, plout, errval] = bair(m, ext, a012, recur, idigit)
            a0 = a012[0]
            a1 = a012[1]
            a2 = a012[2]
            for kk in range(0, nroot - 1):
                workb[kk] = plout[kk]
            # Tidy up the quotient polynomial.
            [delta, qfa, qfb, qfaa, qfab, qfba, qfbb] = qfact(nroot, worka, recur, a1, a0)
            for kk in range(0, nroot - 1):
                workb[kk] = delta[kk]
            [info, z1, z2] = roots(a0, a1, a2, recur)
            # Record node information.
            # If IFAIL is nonzero, then the calculation is going to be rejected.
            nodes = nodes + 1
            icheck[nodes - 1] = ifail
            err[nodes - 1] = errval
            qr[nodes - 1] = z1.real
            qi[nodes - 1] = z1.imag
            nodes = nodes + 1
            icheck[nodes - 1] = ifail
            err[nodes - 1] = errval
            qr[nodes - 1] = z2.real
            qi[nodes - 1] = z2.imag
            nroot = nroot - 2
            # Make the deflated polynomial current.
            i = ip1
            ip1 = ip2
            ip2 = i
            worka = work1 if (1 == ip1) else work2
            workb = work1 if (1 == ip2) else work2
            # Scale the deflated polynomial.
            if nroot <= 0:
                break
            pmax = abs(worka[0])
            for ii in range(1, nroot + 1):
                if pmax < abs(worka[ii]):
                    pmax = abs(worka[ii])
            for ii in range(0, nroot + 1):
                worka[ii] = worka[ii] / pmax
    # Check for poor convergence.
    print("qr = ", qr)
    i = sum(icheck)
    if i != 0:
        print("SOLVE - Warning:")
        print("Poor convergence for some roots.")
        return [1, nodes, qr, qi, icheck, err]
    if nodes != m:
        print("SOLVE - Warning:")
        print("NODE /= M.")
        return [1, nodes, qr, qi, icheck, err]
    # Look for possible imaginary nodes.
    bcomplexnode = False
    for j in range(1, nodes + 1):
        if abs(qi[j - 1]) > 1.0e-22:
            print("node {} is complex: {} + {} i".format(j, qr[j - 1], qi[j - 1]))
            bcomplexnode = True
    if bcomplexnode:
        return [2, nodes, qr, qi, icheck, err]
    return [0, nodes, qr, qi, icheck, err]


def weight(t: list, m: int, n: int, xnode: int, recur, h0: float, nexp: int = -38) -> float:
    """
    WEIGHT calculates quadrature weights.

    This function calculates the quadrature weight associated with the node XNODE
        in the rule whose nodes are defined by the roots of polynomial T.
    The weight is calculated by dividing T by (X-XNODE) to give:
        S(X) = T(X)/(X-XNODE) = SUM (0 to N-1) G(I) * P(I,X).
        S(X) is then divided by (X-XNODE) to give the remainder R.
    The weight is finally given by H0*G(0)/R. If N = M the
    Christoffel-Darboux identity result is used to reduce extreme cancellation effects at high degree.

    Input, real ( kind = 8 ) T(0:N), the coefficients TI of the polynomial
        whose roots define the N preassigned nodes of the quadrature rule and expressed as:
        SUM (I = M to N) (TI/HI) * P(I,X)
    where HI is the integral of W(X) * P(I,X)**2 over the interval
    for which orthogonality with respect the weight W(X) is defined (moment integrals) and P(I,X) is the orthogonal polynomial of degree I.
    T(I-M) holds the value of TI.
    This array should be declared to have at least N-M+1 elements in the calling program.
    Input, integer ( kind = 4 ) M, the lower limit to the expansion of T.
    input, integer ( kind = 4 ) N, the upper limit to the expansion of T.
    Input, real ( kind = 8 ) XNODE, the node whose weight is required
    Input, external RECUR ( ), the function which defines the orthogonal polynomials.
        See EXTEND for a full description.
    Input, real ( kind = 8 ) H0, the integral of the orthogonality weight function over the interval of integration.
        Zero moment integral.  Note that P(0,X) is arbitrarily taken to be 1.0
    Input, integer ( kind = 4 ) NEXP, the largest negative decimal exponent supported on the computer. (Positive number - typical value 38).
    Weights less than approximately 10**(-NEXP) are set to zero when the Christoffel-Darboux identity is used (N = M).
    Output, real ( kind = 8 ) WT, the weight associated with XNODE.
    """
    # Use Christoffel-Darboux result.
    if m == n:
        bk1 = 0.0
        bk2 = 1.0
        dk1 = 0.0
        dk2 = 0.0
        iscale = 0
        k = 0
        [h, d0, e0] = recur(k)
        for k in range(0, n):
            [ck, dk, ek] = recur(k)
            if 0 < k:
                h = - ek * h
            bb = (ck * xnode + dk) * bk2 + ek * bk1
            dd = (ck * xnode + dk) * dk2 + ek * dk1 + ck * bk2
            bk1 = bk2
            bk2 = bb
            dk1 = dk2
            dk2 = dd
            if abs(bk2) > 1.0e-22:
                j = int(math.log(abs(bk2), 10))
                if 2 < abs(j):
                    # Scale to control overflow/underflow.
                    iscale = iscale - 2 * j
                    scale = 10.0 ** j
                    bk2 = bk2 / scale
                    bk1 = bk1 / scale
                    dk1 = dk1 / scale
                    dk2 = dk2 / scale
            if abs(h) < 1.0e-22:
                j = int(math.log(abs(h), 10))
                if 2 <= abs(j):
                    iscale = iscale + j
                    h = h / (10.0 ** j)
        wt = h0 * h / dk2 / bk1
        if abs(wt) < 1.0e-22:
            itest = int(math.log(abs(wt), 10)) + iscale
            if (-nexp) <= itest:
                wt = wt * 10.0 ** iscale
            else:
                wt = 0.0
        return wt
    # General case.
    bk2 = 0.0
    bk1 = 0.0
    rk2 = 0.0
    rk1 = 0.0
    [ck, dk, ek] = recur(n)
    [ck1, dk1, ek1] = recur(n + 1)
    h = 1.0
    iscale = 0
    for k in range(n, 0, -1):
        rs = 0.0
        if m <= k:
            rs = t[k - m] / h
            # Scale and adjust for possible overflow/underflow.
            if nexp < iscale:
                rs = 0.0
            else:
                rs = rs / (10.0 ** iscale)
        bb = rs + (dk + xnode * ck) * bk1 + ek1 * bk2
        bk2 = bk1
        bk1 = bb
        [ckm1, dkm1, ekm1] = recur(k - 1)
        if n != m:
            h = - h * ck / ek / ckm1
        bb = bb * ckm1
        wt = bb + (dkm1 + xnode * ckm1) * rk1 + ek * rk2
        rk2 = rk1
        rk1 = wt
        ck1 = ck
        dk1 = dk
        ek1 = ek
        ck = ckm1
        dk = dkm1
        ek = ekm1
        if abs(bk1) > 1.0e-22:
            j = int(math.log(abs(bk1), 10))
            # Scale to control overflow/underflow.
            if 2 < abs(j):
                iscale = iscale + j
                scale = 10.0 ** j
                bk1 = bk1 / scale
                bk2 = bk2 / scale
                rk1 = rk1 / scale
                rk2 = rk2 / scale
    wt = h0 * bb / wt
    return wt


def extend(m: int, m0n: list, t: list, recur,
           symmet: bool, start: bool, pnodes: list,
           h0: float, nexp: int, idigit: int = 5) -> [int, list, list, list]:
    """
    EXTEND extends a quadrature rule by adding new nodes.

    This function calculates the N+M node quadrature rule composed of N preassigned nodes together
        with M nodes chosen optimally to achieve algebraic degree of precision of at least N+2*M-1.
    The orthogonal system of polynomials associated with the quadrature weight is defined
        generally by the recurrence relation specified in the user supplied function RECUR.

    Input/output, integer ( kind = 4 ) N;
            On input, the number of preassigned nodes, and the upper limit to the expansion.
            On output, if the computation was successful, N is reset to N+M,
            which is the appropriate next input value in cases where this routine is called iteratively.
    Input, integer ( kind = 4 ) M, the number of nodes to be optimally added.
    Input/ouput, integer ( kind = 4 ) M0;
        On input, the lower limit to the expansion of T.
            This is ignored if START is TRUE.  If the computation is successful, the output value of M0 is reset to M,
            which is the appropriate next input value in cases where this routine is called iteratively.
        On output, the lower limit defining the new orthogonal expansion T.  (Set to M).
    Input/output, real ( kind = 8 ) T(0:max(N-M0,M));
        On input, the coefficients TI of the polynomial whose roots define the N preassigned nodes
            of the quadrature rule and expressed as:
            SUM (I = M0 to N) (TI/HI) * P(I,X)
            where HI is the integral of W(X) * P(I,X)^2 over the interval for which orthogonality
            with respect the weight W(X) is defined (moment integrals) and P(I,X) is the orthogonal polynomial of degree I.
            Element T(I-M0) holds the value of TI.
            Note that T is either,
            (1) provided explicitly,
            (2) generated automatically from the N preassigned nodes given in PNODES(*) (if START is TRUE.)
            (3) generated from a previous call to the function.
            This array should be declared to have at least max(N-M0+1,M+1) elements in the calling program.
            The service function TRANSF can be used to transform the expansion to the required input form
            if desired with the parameter IFLAG set to 1.
        On output, the coefficients TI of the new orthogonal expansion whose roots are the nodes of the extended quadrature
            (that is, the preassigned nodes plus the extended nodes) and expressed as:
            SUM (I = M to N+M) (TI/HI) * P(I,X)
            T(I-M) holds the value of TI.
            (For definitions see description of input argument T).
            This polynomial can be used as input for further extensions.
            The service function TRANSF can be used to remove the moment factors from the expansion
            if desired with the parameter IFLAG set to 0.
    Input, external RECUR ( ), the user supplied function which defines the orthogonal polynomials.
        Given K, CALL RECUR(K,C,D,E) gives the coefficients C,D and E such that,
        P(K+1,X) = ( C * X + D ) * P(K,X) + E * P(K-1,X)
        The parameters are defined as follows:
            K = Index
            C, D, E = Parameters in the recurrence relation (functions of K)
    Input, logical SYMMET =
        * FALSE, if no advantage is to be taken of symmetry, if any, about x = 0 in the interval of integration and the orthogonality weight function.
            Note that if symmetry in fact does exist setting this parameter to zero will still produce correct results - only efficiency is effected.
        * TRUE, if the extended rule computations should exploit symmetry about x = 0 in the interval of integration and the orthogonality  weight function.
            This reduces the size of the system of linear equations determining EXT by a factor of about 2 (see WORKA).
            If symmetry does not in fact exist erroneous results will be produced.
    Input, logical START,
        * TRUE, then the polynomial T is generated to have the preassigned nodes (PNODES) as its roots.
        * FALSE. then the supplied values of the coefficients of T are used directly.
    Input, real ( kind = 8 ) PNODES(N+M),
        On input, the preassigned nodes.
        On output, the nodes of the extended quadrature rule made up from the original preassigned nodes and the new optimally extended nodes.
        These values can be used in subsequent iterative use.
    Input, real ( kind = 8 ) H0, the integral of the orthogonality weight function over the interval of integration.
        Zero moment integral.
    Input, integer ( kind = 4 ) NEXP, the largest negative decimal exponent supported on the computer.
        (Positive number - typical value 38).
        Weights less than approximately 10^(-NEXP) are set to zero when the Christoffel-Darboux identity is used (N = M).
        This may be set to INT(LOG10(X1MACH(2))) where X is set to correspond to the appropriate precision in the PORT library.
    Input, integer ( kind = 4 ) IDIGIT, the node convergence parameter, an integer greater than 0.
        An attempt is made to calculate the nodes to the maximum accuracy possible by the machine precision available.
        IDIGIT controls the assessment procedure to take account of round-off errors and specifies the number of
        least significan decimal digits that can be ignored (i.e. attributed to round-off) in the computed relative error.
        A typical value is 5.
    Output, real ( kind = 8 ) WT(N+M), the quadrature weights for the extended rule associated with the nodes in PNODES.
    Output, integer ( kind = 4 ) NODES, the number of extended nodes found.
        Normally equals M, but see IFAIL.
    Output, real ( kind = 8 ) QR(M), the real parts of the extended nodes.
    Output, real ( kind = 8 ) QI(M), the imaginary parts of the extended nodes (1,..,NODES).
        Hopefully these values are zero!
    Output, real ( kind = 8 ) ERR(M), a measure of the relative error in the nodes.
        This may be inspected if the convergence error flag has been raised (IFLAG = 3) to decide if the nodes in question are acceptable.
        (ERR(*) actually gives the mean last correction to the quadratic factor in the generalized Bairstow root finder (see BAIR).
    Output, real ( kind = 8 ) EXT(M+1), the coefficients of the polynomial whose roots are the  extended nodes (QRS(*),QIS(*)) and expressed as:
        EXT = SUM (0 <= I <= M) EXT(I) * P(I,X).
    Output, integer ( kind = 4 ) IWORK(max(M,N)), node convergence flags.
        Elements 1 to NODES give information on the convergence of the roots of the polynomial EXT corresponding to each extended node.
        * 0, Convergence of I th root satisfactory;
        * 1, Convergence of I th root unsatisfactory.
    Workspace, real ( kind = 8 ) WORKA(max(M+1,N+1),max(M+1,N+1)).
        If SYMMET = TRUE, the dimensions can be reduced to max(M/2+1,N) by max(M/2+1,N+1).
    Input, integer ( kind = 4 ) LDA, the leading dimension of WORKA.
    Input, real ( kind = 8 ) WORKB(2*M+1,3).
    Input, integer ( kind = 4 ) LDB, the leading dimension of WORKB.
    Output, integer ( kind = 4 ) IFLAG  = 0, No error detected
        * 1, The linear system of equations defining the polynomial whose roots are the extended nodes became singular or very  ill-conditioned.   (FATAL).
        * 2, The linear system of equations used to generate the polynomial T when START is TRUE became singular or very ill-conditioned. (FATAL).
        * 3, Poor convergence has been detected in the calculation of the roots of EXT corresponding to the new nodes or all nodes have not been found (M not equal to NODES). See also ERR(*).
        * 4, Possible imaginary nodes detected.
        * 5, Value of N and M incompatible for SYMMET = TRUE. Both cannot be odd. (FATAL)
        * 6, Test of new quadrature rule has failed.
    """
    tol = 0.0000001
    iflag = 0
    nodes = 0
    ideg = m0n[1] + 2 * m - 1
    wt = [0.0 for _ in range(0, m + m0n[1])]
    if symmet and ((m0n[1] % 2) == 1) and ((m % 2) == 1):
        print("EXTEND - Fatal error!")
        print("N and M cannot both be odd for a symmetric rule.")
    # If required, generate the initial T polynomial corresponding to prescribed preassigned nodes.
    if start and m0n[1] != 0:
        [ierr, outt] = assign(m0n[1], pnodes, recur)
        m0n[0] = 0
        for i in range(0, m0n[1]):
            t[i] = outt[i]
        if ierr != 0:
            print("EXTEND - Fatal error!")
            print("Unable to generate the initial T polynomial for the preassigned nodes.")
    nlast = m0n[1]
    # Generate extended expansion coefficients and overwrite T.
    worka = [[0.0 for _ in range(0, max(m + 1, m0n[1] + 1))] for _ in range(0, max(m + 1, m0n[1] + 1))]
    print("before gener m0, n = ", m0n, " m =", m)
    print("before gener t = ", t)
    [ierr, ext] = gener(t, m0n, m, recur, symmet, worka)
    print("after gener m0, n = ", m0n, " m =", m)
    print("after gener t = ", t)
    print("after gener ext = ", ext)
    if ierr != 0:
        print("EXTEND - Fatal error!")
        print("gener failed")
    # Find extended nodes as roots of EXT(*).
    workb = [[0.0 for _ in range(0, 2 * m + 1)] for _ in range(0, 3)]
    print("ext = ", ext, " m = ", m)
    [ierr, nodes, qr, qi, err, iwork] = solve(ext, m, symmet, recur, workb[0], workb[1], idigit)
    print("qr = ", qr, " nodes = ", nodes, " nlast = ", nlast, " m = ", m)
    if ierr != 0:
        print("EXTEND - Fatal error!")
        print("solve failed")
    # Accumulate nodes for extended rule.
    for i in range(0, m):
        pnodes[nlast + i] = qr[i]
    # Reorder.
    rsort(pnodes, m0n[1], 1)
    # Compute weights (only for positive nodes if symmetric).
    print("pnodes = ", pnodes)
    num = m0n[1]
    if symmet:
        num = (m0n[1] + 1) // 2
    for i in range(1, num + 1):
        wt[i - 1] = weight(t, m0n[0], m0n[1], pnodes[i - 1], recur, h0, nexp)
        if symmet:
            wt[m0n[1] - i] = wt[i - 1]
    return [nodes, wt, qr, qi, ext]


def GenerateGaussianPattersonOneOrder(order: int) -> [list, list]:
    if 1 == order:
        return [[0.0], [2.0]]
    if 2 == order:
        return [[
            -0.77459666924148337704,
            0.0,
            0.77459666924148337704
        ], [
            0.555555555555555555556,
            0.888888888888888888889,
            0.555555555555555555556
        ]]
    length = (1 << order) - 1
    n = 0
    m = 3
    m0 = 0
    t = [0.0 for _ in range(0, length)]
    pnodes = [0.0 for _ in range(0, length)]
    wt = [0.0 for _ in range(0, length)]
    t[0] = 1.0
    for i in range(1, order):
        # Calculate the extension.
        print("========= Start i =", i, " === pnodes = ", pnodes, " wt = ", wt)
        h0 = 2.0
        ideg = n + 2 * m - 1
        m0n = [m0, n]
        [_, wt, _, _, _] = extend(m, m0n, t, recura, True, False, pnodes, h0, 1021, 8)
        m0 = m0n[0]
        n = m0n[1]
        m = n + 1
        print("========= End i =", i, " === pnodes = ", pnodes, " wt = ", wt)
    return [[pnodes[length - i - 1] for i in range(0, length)], [wt[length - i - 1] for i in range(0, length)]]





