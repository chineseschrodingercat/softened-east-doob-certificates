"""Small-ring floating-point diagnostics; not a rigorous certificate."""
import argparse
import json
from pathlib import Path

import numpy as np


def check(n, kappa, s, candidate):
    c = 0.3
    size = 1 << n
    states = ((np.arange(size)[:, None] >> np.arange(n)) & 1).astype(float)
    pi = np.prod(np.where(states == 1, c, 1 - c), axis=1)
    rates = (kappa + np.roll(states, 1, axis=1)) * np.where(states == 0, c, 1 - c)
    ci = np.exp(-s) * rates
    flip = np.arange(size)[:, None] ^ (1 << np.arange(n))
    tilted = np.zeros((size, size))
    tilted[np.arange(size)[:, None], flip] = ci
    a0 = -rates.sum(axis=1)
    np.fill_diagonal(tilted, a0)
    symmetric = np.sqrt(pi)[:, None] * tilted / np.sqrt(pi)[None, :]
    assert np.max(np.abs(symmetric - symmetric.T)) < 1e-12
    eigenvalues, eigenvectors = np.linalg.eigh(symmetric)
    lam = eigenvalues[-1]
    vector = eigenvectors[:, -1]
    if vector.sum() < 0:
        vector = -vector
    assert np.all(vector > 0)
    h = vector / np.sqrt(pi)
    logh = np.log(h)
    singles = states.sum(axis=1)
    pairs = (states * np.roll(states, -1, axis=1)).sum(axis=1)
    triples = (states * np.roll(states, -1, axis=1) * np.roll(states, -2, axis=1)).sum(axis=1)
    U = {
        "constant": np.zeros(size),
        "pair": 0.11 * singles - 0.035 * pairs,
        "triple": 0.05 * singles - 0.025 * pairs + 0.04 * triples,
    }[candidate]
    nu = pi * np.exp(2 * U)
    nu /= nu.sum()
    b = ci * np.exp(U[flip] - U[:, None])
    q = ci * np.exp(logh[flip] - logh[:, None])
    V = a0 + b.sum(axis=1)
    gap = float((lam - nu @ V) / n)
    kl_rate_density = float(nu @ (b * np.log(b / q) - b + q).sum(axis=1) / n)
    cancellation = float(nu @ (b * np.log(b / q)).sum(axis=1))
    G = np.zeros_like(tilted)
    G[np.arange(size)[:, None], flip] = b
    np.fill_diagonal(G, -b.sum(axis=1))
    row = dict(n=n, c=c, kappa=kappa, s=s, candidate=candidate,
               spectral_gap_from_center=gap, kl_rate_density=kl_rate_density,
               identity_error=abs(gap - kl_rate_density),
               gradient_cancellation=abs(cancellation),
               stationarity_residual=float(np.max(np.abs(nu @ G))),
               perron_residual=float(np.max(np.abs(tilted @ h - lam * h)) / np.max(h)))
    assert gap > -1e-11 and kl_rate_density > -1e-11
    for key in ("identity_error", "gradient_cancellation", "stationarity_residual", "perron_residual"):
        assert row[key] < 1e-10, row
    return row


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise SystemExit("Refusing to overwrite an existing diagnostic result")
    rows = [check(n, kappa, s, candidate)
            for kappa, s in [(1.0, 0.02), (0.1, 0.01), (0.01, 0.001)]
            for n in (3, 4, 5, 6)
            for candidate in ("constant", "pair", "triple")]
    result = dict(status="PASS", method="ordinary float64 diagnostics; not strict certification",
                  numpy_version=np.__version__, cases=len(rows),
                  max_identity_error=max(row["identity_error"] for row in rows), results=rows)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2)
    print(json.dumps({key: value for key, value in result.items() if key != "results"}, indent=2))
