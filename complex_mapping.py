from dataclasses import dataclass, field
from typing import Self, Optional

import matplotlib as plt
import numpy as np
from scipy.special import jv, hankel1


@dataclass
class ComplexMapping:
    num: int
    numPoints: Optional[int] = 100
    cols: Optional[int] = None
    rows: Optional[int] = None
    zeta0: complex = 0+0j
    kappa0: float = 0.
    zetaMatrix:  Optional[np.ndarray] = field(default_factory=lambda: np.ndarray([]))
    kappaMatrix: Optional[np.ndarray] = field(default_factory=lambda: np.ndarray([]))

    def delta(self: Self, zetas: np.ndarray, kappas: np.ndarray) -> np.ndarray:
        n = self.num
        return abs(jv(n, kappas * zetas) * hankel1(n + 1, 1, kappas) -
                   zetas * hankel1(n, 1, kappas) * jv(n + 1, kappas * zetas))

    def dKappa_dZeta(self: Self, zeta: complex, kappa: float) -> complex:
        n = self.num
        
        JnKappaZeta =  jv(n,     kappa * zeta)
        Jn1KappaZeta = jv(n + 1, kappa * zeta)
        
        H1nKappa  = hankel1(n,     1, kappa)
        H1n1Kappa = hankel1(n + 1, 1, kappa)
        
        dJn_dZeta  = kappa * (jv(n - 1, kappa * zeta) - jv(n + 1, kappa * zeta)) / 2
        dJn1_dZeta = kappa * (jv(n,     kappa * zeta) - jv(n + 2, kappa * zeta)) / 2
        dDeltaZeta = dJn_dZeta * H1n1Kappa - (H1nKappa * Jn1KappaZeta + zeta * H1nKappa * dJn1_dZeta)
        
        dJn_dKappa  = zeta * (jv(n - 1, kappa * zeta) - jv(n + 1, kappa * zeta)) / 2
        dJn1_dKappa = zeta * (jv(n,     kappa * zeta) - jv(n + 2, kappa * zeta)) / 2
        dH1n_dKappa  = (hankel1(n - 1, 1, kappa) - hankel1(n + 1, 1, kappa)) / 2
        dH1n1_dKappa = (hankel1(n,     1, kappa) - hankel1(n + 2, 1, kappa)) / 2;

        dDeltaKappa = (dJn_dKappa * H1n1Kappa + JnKappaZeta * dH1n1_dKappa) - \
            (zeta * dH1n_dKappa * Jn1KappaZeta + zeta * H1nKappa * dJn1_dKappa)
        return -dDeltaZeta / dDeltaKappa

    def kappaSolution(self: Self, col: int) -> np.ndarray:
        kappaSol = np.zeros(self.rows)
        zetaSpan = self.zetaMatrix[col, :]
        kappaSol[0] = self.kappa0
        for i in range(self.rows - 1):
            h = zetaSpan[i+1] - zetaSpan[i]
            currentZeta = zetaSpan[i]
            currentKappa = kappaSol[i]
            k1 = h * self.dKappa_dZeta(currentZeta,         currentKappa)
            k2 = h * self.dKappa_dZeta(currentZeta + h / 2, currentKappa + k1 / 2)
            k3 = h * self.dKappa_dZeta(currentZeta + h / 2, currentKappa + k2 / 2)
            k4 = h * self.dKappa_dZeta(currentZeta + h,     currentKappa + k3)
            kappaSol[i+1] = currentKappa + (k1 + 2 * k2 + 2 * k3 + k4) / 6
        return kappaSol

    def calcKappa(self: Self) -> None:
        shape = self.cols, self.rows = self.zetaMatrix.shape
        self.kappaMatrix = np.zeros(shape)
        for col in range(self.cols):
            self.kappaMatrix[col, :] = self.kappaSolution(col)

    # TODO: other methods
