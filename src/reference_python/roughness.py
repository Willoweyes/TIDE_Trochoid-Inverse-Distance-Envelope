# -*- coding: utf-8 -*-
"""Shared ISO 25178 areal roughness metrics — used by ALL methods for fair comparison."""
import numpy as np


def compute_roughness(z_um):
    """z_um: array of surface heights in micrometers. Returns dict of areal params."""
    z = np.asarray(z_um, dtype=float).ravel()
    z = z[~np.isnan(z)]
    if z.size == 0:
        return dict(Sa=np.nan, Sq=np.nan, Sz=np.nan, Sp=np.nan, Sv=np.nan,
                    Ssk=np.nan, Sku=np.nan, n_points=0)
    zc = z - z.mean()
    Sq = float(np.sqrt(np.mean(zc**2)))
    return dict(
        Sa=float(np.mean(np.abs(zc))),
        Sq=Sq,
        Sz=float(z.max() - z.min()),
        Sp=float(zc.max()),
        Sv=float(-zc.min()),
        Ssk=float(np.mean(zc**3) / Sq**3) if Sq > 1e-12 else 0.0,
        Sku=float(np.mean(zc**4) / Sq**4) if Sq > 1e-12 else 0.0,
        n_points=int(z.size),
    )
