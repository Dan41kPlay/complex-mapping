from dataclasses import dataclass, field
from typing import Self, Optional
from decimal import Decimal, getcontext

import cv2
import matplotlib as plt
import numpy as np
#import skimage as sk
from scipy.special import jv, hankel1, fresnel


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
    def line(self: Self) -> None:
        ray_length = 0.0001
        numRays = 500
        t_param = np.linspace(0, 1, self.numPoints).reshape(1, -1)
        angles = np.linspace(0, 2*np.pi, numRays + 1)[:-1].reshape(-1, 1).conj()
        directions = np.exp(1j * angles)
        self.zetaMatrix = (self.zeta0 + ray_length * directions * t_param).conj().T

    def parabols(self: Self) -> None:
        a = self.zeta0.real
        b = self.zeta0.imag
        c = -0.8
        C = 0.5
        n_values = np.arange(2, 13, 1)

        num_points_curved = 1000
        num_points_linear = 1000
        t1 = np.linspace(a, c, num_points_curved)
        t2 = np.linspace(c, C, num_points_linear + 1)[1:]

        self.zetaMatrix = np.zeros((max(t1.shape) + max(t2.shape), max(n_values.shape)), dtype = complex)

        for k in range(max(n_values.shape)):
            n = n_values[k]
            zeta_curved = t1 + 1j * b * ((t1 - c) / (a - c)) ** n
            zeta_linear = t2 + 0j
            self.zetaMatrix[:, k] = np.concatenate([zeta_curved, zeta_linear])

    def powerCurveVectorized(self: Self) -> None:
        t = np.linspace(3, 5.001, self.numPoints).reshape(1, -1)
        n = np.arange(1, 13, 1).reshape(-1, 1)
        self.zetaMatrix = t - 1j * ((t - 3) / (-2)) ** n

    def kornus(self: Self) -> None:
        num_points = 500;
        t = np.linspace(-np.pi, np.pi, num_points)
        a_values = np.linspace(0.00000001, 0.0000001, 200).reshape(1, -1)
        y, x = fresnel(t)
        base_spiral = (x + 1j * y).reshape(1, -1).T
        self.zetaMatrix = (base_spiral * a_values + self.zeta0).T

    def spirals(self: Self) -> None:
        t_values = np.linspace(0, 8*np.pi, self.numPoints).reshape(1, -1)
        k_values = np.arange(0.00001, 0.006 + 0.00001/2, 0.00001).reshape(-1, 1)
        self.zetaMatrix = self.zeta0 + k_values * t_values * np.exp(1j * t_values)

    def koch(self: Self) -> None:
        it = 6
        angles = np.arange(0, 11, 1)
        L = 0.0000000001
        z = np.linspace(0, L, self.numPoints).reshape(-1, 1)
        for i in range(it):
            dz = np.diff(z, axis=0)
            p1 = z[:-1]
            p2 = p1 + dz/3
            p3 = p2 + dz/3 * np.exp(1j * np.pi / 3)
            p4 = p1 + 2*dz/3
            z_new = np.column_stack((p1,p2,p3,p4)).reshape(-1,1)
            z = np.concatenate([z_new, z[-1].reshape(-1, 1)])
        self.zetaMatrix = np.zeros((max(angles.shape), max(z.shape)), dtype = complex)
        for k in range(max(angles.shape)):
            self.zetaMatrix[k, :] = self.zeta0 + z.T * np.exp(1j * np.deg2rad(angles[k]))

    def mandelbrot(self: Self) -> None:
        res = 1000
        max_iter = 150
        x = np.linspace(-2.1, 0.6, res)
        y = np.linspace(-1.2, 1.2, res)
        X, Y = np.meshgrid(x, y)
        C = X + 1j * Y

        Z_calc = np.zeros_like(C)
        Inside = np.ones_like(C, dtype = bool)
        for n in range(max_iter):
            mask = np.abs(Z_calc) <= 2
            Z_calc[mask] = Z_calc[mask] ** 2 + C[mask]
            Inside = Inside & (np.abs(Z_calc) <= 2)

        #Под вопросом
        #boundaries = sk.measure.find_contours(np.pad(Inside, pad_width=1, mode='constant', constant_values=0), 0.9)
        #boundaries = [i[np.append((np.diff(i[:,0])!=0)|(np.diff(i[:,1])!=0), True)]-1 for i in [np.round(i).astype(int) for i in boundaries]]
        '''
        Так как skimage возвращает точки, ноходящиеся между 0 и 1, их приходится округлять и удалять лишние чтобы повторить поведение Matlab,
        возможно лучше использовать opencv (вероятно быстрее работает)
        '''
        #Под вопросом
        Inside = Inside.astype(np.uint8) * 255
        boundaries, _ = cv2.findContours(Inside, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        mainIdx = np.argmax([max(i.shape) for i in boundaries])
        mainBoundary = boundaries[mainIdx].reshape(-1, 2)
        mainBoundary = np.concatenate((mainBoundary, mainBoundary[0].reshape(-1, 2)))
        z_edge = x[mainBoundary[:, 0]] + 1j * y[mainBoundary[:, 1]]

        refIdx = np.argmax(z_edge.real)
        #Под вопросом
        reference_point = z_edge[refIdx]
        '''
        Из-за того, что массив точек mainBoundary смещён относительно матлаба (алгоритмы начинают контур из разных точек), 
        эта точка получается другой (тк массив может содержать несколько точек с максимальным real значением, и max берёт первую точку в массиве), 
        поэтому комплексные значения итоговой матрицы отличаются
        Имеет ли это значение?
Python (против часовой)
[[False False False False False False False False False False]
 [False False False False False False False False False False]
 [False False False False False False False False False False]
 [False False False False False  True  True  True  True False]
 [False False False  True  True  True  True  True  True False]
 [False False False  True  True  True  True  True  True False]
 [False False False False False  True  True  True  True False]
 [False False False False False False False False False False]
 [False False False False False False False False False False]
 [False False False False False False False False False False]]
[[5 3]
 [4 4]
 [3 4]
 [3 5]
 [4 5]
 [5 6]
 [6 6]
 [7 6]
 [8 6]
 [8 5]
 [8 4]
 [8 3]
 [7 3]
 [6 3]
 [5 3]]
 Matlab (по часовой)
   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   1   1   1   1   0
   0   0   0   1   1   1   1   1   1   0
   0   0   0   1   1   1   1   1   1   0
   0   0   0   0   0   1   1   1   1   0
   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0
   0   0   0   0   0   0   0   0   0   0

     5     4
     5     5
     4     6
     4     7
     4     8
     4     9
     5     9
     6     9
     7     9
     7     8
     7     7
     7     6
     6     5
     6     4
     5     4
        '''
        #Под вопросом
        offset = self.zeta0 - reference_point
        z_base = z_edge + offset
        scales = np.arange(0.000000001, 0.000000001+0.00000000005/2,0.00000000005)
        self.zetaMatrix = np.zeros((max(z_base.shape), max(scales.shape)), dtype = complex)
        for k in range(max(scales.shape)):
            s = scales[k]
            self.zetaMatrix[:, k] = self.zeta0 + s * (z_base - self.zeta0)
        self.zetaMatrix = self.zetaMatrix.T

    def julia(self: Self) -> None:
        res = 1000
        max_iter = 150
        C_constant = -0.123 + 0.745j
        x = np.linspace(-1.5, 1.5, res)
        y = np.linspace(-1.5, 1.5, res)
        X, Y = np.meshgrid(x, y)
        Z_calc = X + 1j * Y
        Inside = np.ones_like(Z_calc, dtype=bool)
        for n in range(max_iter):
            mask = np.abs(Z_calc) <= 2
            Z_calc[mask] = Z_calc[mask] ** 2 + C_constant
            Inside = Inside & (np.abs(Z_calc) <= 2)
        Inside = Inside.astype(np.uint8) * 255
        boundaries, _ = cv2.findContours(Inside, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        mainIdx = np.argmax([max(i.shape) for i in boundaries])
        mainBoundary = boundaries[mainIdx].reshape(-1, 2)
        mainBoundary = np.concatenate((mainBoundary, mainBoundary[0].reshape(-1, 2)))
        z_edge = x[mainBoundary[:, 0]] + 1j * y[mainBoundary[:, 1]]
        refIdx = np.argmax(z_edge.real)
        reference_point = z_edge[refIdx]
        offset = self.zeta0 - reference_point
        z_base = z_edge + offset
        scales = np.arange(0.000000000001, 0.00000000001 + 0.0000000000001 / 2, 0.0000000000001)
        self.zetaMatrix = np.zeros((max(z_base.shape), max(scales.shape)), dtype=complex)
        for k in range(max(scales.shape)):
            s = scales[k]
            self.zetaMatrix[:, k] = self.zeta0 + s * (z_base - self.zeta0)
        self.zetaMatrix = self.zetaMatrix.T

    def circles(self: Self) -> None:
        Radii = np.linspace(0.00000011, 0.000012, 200).reshape(-1, 1)
        theta = np.linspace(0, 2 * np.pi, self.numPoints).reshape(1, -1)
        Theta_grid, Radii_grid = np.meshgrid(theta, Radii)
        Centers = self.zeta0 + Radii
        Circles_adjusted = Centers + Radii_grid * (np.exp(1j * Theta_grid) - 1.) - Radii
        self.zetaMatrix = Circles_adjusted.T.conj()

