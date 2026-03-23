import numpy as np

def mcs(fcn, data, u, v, prt=1, smax=None, nf=None, stop=None,
        iinit=None, local=50, gamma=None, hess=None):
    """
    Python translation of MATLAB mcs: Multi-level Coordinate Search
    """
    u = np.asarray(u, dtype=float)
    v = np.asarray(v, dtype=float)
    n = len(u)
    
    # Check box bounds
    if np.any(v < u):
        raise ValueError('Incompatible box bounds')
    if np.any(u == v):
        raise ValueError('Degenerate box bound')
    
    # Default values for parameters
    if smax is None: smax = 5 * n + 10
    if nf is None: nf = 50 * n ** 2
    if stop is None: stop = [3 * n]
    if gamma is None: gamma = np.finfo(float).eps
    if hess is None: hess = np.ones((n, n))
    if iinit is None:
        if not (np.any(np.isinf(u)) or np.any(np.isinf(v))):
            iinit = 0
        else:
            iinit = 1

    # Placeholder for later outputs
    xbest, fbest, xmin, fmi, ncloc, flag = None, None, None, None, 0, None

    # Next: initialization list and function evaluation
    x0, l, L = init_list(u, v, iinit)

    # Evaluate at these points
    f0, istar = eval_init_list(fcn, data, x0, l, L)
    
    # For now, print them out to check
    print("Initial points x0:\n", x0)
    print("Function values at x0:", f0)
    print("Index of best init point:", istar)

    # --- Next steps: main optimization loop, etc. ---

    # For now, return initial data (to check correctness)
    return x0, f0, istar

# ---- Helper functions ----

def init_list(u, v, iinit):
    """Constructs the initialization list (x0, l, L) based on iinit flag."""
    n = len(u)
    if iinit == 0:
        # Corners and midpoint
        x0 = np.column_stack((u, (u+v)/2, v))
        l = 2 * np.ones(n, dtype=int)  # Index of the midpoint in each dim
        L = 3 * np.ones(n, dtype=int)
    elif iinit == 1:
        x0 = np.zeros((n, 3))
        for i in range(n):
            if u[i] >= 0:
                x0[i, 0] = u[i]
                x0[i, 2], x0[i, 1] = subint(u[i], v[i])
                x0[i, 1] = 0.5 * (x0[i, 0] + x0[i, 2])
            elif v[i] <= 0:
                x0[i, 2] = v[i]
                x0[i, 1], x0[i, 0] = subint(v[i], u[i])
                x0[i, 1] = 0.5 * (x0[i, 0] + x0[i, 2])
            else:
                x0[i, 1] = 0
                _, x0[i, 0] = subint(0, u[i])
                _, x0[i, 2] = subint(0, v[i])
        l = 2 * np.ones(n, dtype=int)
        L = 3 * np.ones(n, dtype=int)
    elif iinit == 2:
        x0 = np.column_stack((
            (5*u + v)/6,
            (u + v)/2,
            (u + 5*v)/6
        ))
        l = 2 * np.ones(n, dtype=int)
        L = 3 * np.ones(n, dtype=int)
    else:
        raise NotImplementedError("Custom/self-defined initialization lists not supported yet.")
    return x0, l, L

def subint(a, b):
    """MATLAB: [c,d] = subint(a, b) -- returns two points between a and b."""
    # We'll use points at 1/3 and 2/3 between a and b
    c = (2*a + b) / 3
    d = (a + 2*b) / 3
    return c, d

def eval_init_list(fcn, data, x0, l, L):
    """
    Evaluate the objective at each initial point in x0.
    Return function values and index of the best.
    """
    n, m = x0.shape
    f0 = np.zeros(m)
    for j in range(m):
        f0[j] = fcn(data, x0[:, j])
    istar = np.argmin(f0)
    return f0, istar



