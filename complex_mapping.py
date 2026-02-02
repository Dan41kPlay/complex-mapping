from dataclasses import dataclass, field
from typing import Self, Optional
from decimal import Decimal, getcontext
from concurrent import futures

import cv2
import matplotlib.pyplot as plt
import numpy as np
import mpmath as mp
#import skimage as sk
from scipy.special import fresnel

mp.mp.dps = 24

@dataclass
class ComplexMapping:
    num: int
    numPoints: Optional[int] = 100
    cols: Optional[int] = None
    rows: Optional[int] = None
    zeta0:  Optional[mp.mpc] = mp.mpc(0, 0)
    kappa0: Optional[mp.mpf] = mp.mpf(0)
    zetaMatrix:  Optional[np.ndarray] = field(default_factory=lambda: np.ndarray([], dtype=mp.mpc))
    kappaMatrix: Optional[np.ndarray] = field(default_factory=lambda: np.ndarray([], dtype=mp.mpf))

    def delta(self: Self, zetas: np.ndarray, kappas: np.ndarray) -> np.ndarray:
        n = self.num
        '''return abs(mp.besselj(n, kappas * zetas) * mp.hankel1(n + 1, 1, kappas) -
                   zetas * mp.hankel1(n, 1, kappas) * mp.besselj(n + 1, kappas * zetas))'''
        zeta_list = zetas.flatten().tolist()
        kappa_list = kappas.flatten().tolist()
        result = []
        for z, k in zip(zeta_list, kappa_list):
            z_mp = mp.mpc(z)
            k_mp = mp.mpf(k)
            val = mp.besselj(n, k_mp * z_mp) * mp.hankel1(n + 1, k_mp) - \
                  z_mp * mp.hankel1(n, k_mp) * mp.besselj(n + 1, k_mp * z_mp)
            result.append(float(mp.fabs(val)))
        return np.array(result).reshape(zetas.shape)

    def dKappa_dZeta(self: Self, zeta: mp.mpc, kappa: mp.mpf) -> mp.mpc:
        n = self.num
        
        JnKappaZeta =  mp.besselj(n,     kappa * zeta)
        Jn1KappaZeta = mp.besselj(n + 1, kappa * zeta)
        
        H1nKappa  = mp.hankel1(n,     kappa)
        H1n1Kappa = mp.hankel1(n + 1, kappa)
        
        dJn_dZeta  = kappa * (mp.besselj(n - 1, kappa * zeta) - mp.besselj(n + 1, kappa * zeta)) / 2
        dJn1_dZeta = kappa * (mp.besselj(n,     kappa * zeta) - mp.besselj(n + 2, kappa * zeta)) / 2
        dDeltaZeta = dJn_dZeta * H1n1Kappa - (H1nKappa * Jn1KappaZeta + zeta * H1nKappa * dJn1_dZeta)
        
        dJn_dKappa  = zeta * (mp.besselj(n - 1, kappa * zeta) - mp.besselj(n + 1, kappa * zeta)) / 2
        dJn1_dKappa = zeta * (mp.besselj(n,     kappa * zeta) - mp.besselj(n + 2, kappa * zeta)) / 2
        dH1n_dKappa  = (mp.hankel1(n - 1, kappa) - mp.hankel1(n + 1, kappa)) / 2
        dH1n1_dKappa = (mp.hankel1(n,     kappa) - mp.hankel1(n + 2, kappa)) / 2;

        dDeltaKappa = (dJn_dKappa * H1n1Kappa + JnKappaZeta * dH1n1_dKappa) - \
            (zeta * dH1n_dKappa * Jn1KappaZeta + zeta * H1nKappa * dJn1_dKappa)
        return -dDeltaZeta / dDeltaKappa

    def kappaSolution(self: Self, col: int) -> np.ndarray:
        kappaSol = np.zeros(self.rows, dtype=mp.mpf)
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
            kappaSol[i+1] = currentKappa + mp.mpf(1) / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
        return kappaSol

    def calcKappa(self: Self) -> None:
        def acalcKappaCol(col):
            self.kappaMatrix[col, :] = self.kappaSolution(col)
            
        shape = self.cols, self.rows = self.zetaMatrix.shape
        self.kappaMatrix = np.zeros(shape, dtype=mp.mpf)
        executor = futures.ProcessPoolExecutor(10)
        tasks = [executor.submit(acalcKappaCol, col) for col in range(self.cols)]
        futures.wait(tasks)
            
    # Matvey
    def line(self: Self, numRays: int = 500, rayLength: mp.mpf = mp.mpf('.0001')) -> None:
        tParam = np.linspace(0, 1, self.numPoints, dtype=mp.mpf).reshape(1, -1)
        angles = np.linspace(0, 2 * mp.pi, numRays + 1)[:-1].reshape(-1, 1).conj()
        directions = np.array([mp.e**(1j * angle) for angle in angles])
        self.zetaMatrix = (self.zeta0 + rayLength * directions * tParam).conj().T

    def parabolas(self: Self) -> None:
        a = self.zeta0.real
        b = self.zeta0.imag
        c = mp.mpf('-.8')
        C = mp.mpf('.5')
        n_values = np.arange(2, 13, 1, dtype=mp.mpf)

        num_points_curved = 1000
        num_points_linear = 1000
        t1 = np.linspace(a, c, num_points_curved)
        t2 = np.linspace(c, C, num_points_linear + 1)[1:]

        self.zetaMatrix = np.zeros((max(t1.shape) + max(t2.shape), max(n_values.shape)), dtype=mp.mpc)

        for k in range(max(n_values.shape)):
            n = n_values[k]
            zeta_curved = t1 + 1j * b * ((t1 - c) / (a - c)) ** n
            zeta_linear = t2 + 0j
            self.zetaMatrix[:, k] = np.concatenate([zeta_curved, zeta_linear])

    def powerCurveVectorized(self: Self) -> None:
        t = np.linspace(3, mp.mpf('5.001'), self.numPoints).reshape(1, -1)
        n = np.arange(1, 13, 1, dtype=mp.mpf).reshape(-1, 1)
        self.zetaMatrix = t - 1j * ((t - 3) / (-2)) ** n

    def kornus(self: Self) -> None:
        num_points = 500
        t = np.linspace(-np.pi, np.pi, num_points, dtype=mp.mpf)
        a_values = np.linspace(mp.mpf('1e-8'), mp.mpf('1e-8'), 200).reshape(1, -1)
        y, x = fresnel(t)
        base_spiral = (x + 1j * y).reshape(1, -1).T
        self.zetaMatrix = (base_spiral * a_values + self.zeta0).T

    def spirals(self: Self) -> None:
        t_values = np.linspace(0, 8 * mp.pi, self.numPoints).reshape(1, -1)
        k_values = np.arange(mp.mpf('1e-5'), mp.mpf('6e-3') + mp.mpf('5e-7'), mp.mpf('1e-5')).reshape(-1, 1)
        self.zetaMatrix = self.zeta0 + k_values * t_values * np.exp(1j * t_values)

    def koch(self: Self) -> None:
        it = 6
        angles = np.arange(0, 11, 1)
        L = mp.mpf('1e-10')
        z = np.linspace(0, L, self.numPoints).reshape(-1, 1)
        for i in range(it):
            dz = np.diff(z, axis=0)
            p1 = z[:-1]
            p2 = p1 + dz / 3
            p3 = p2 + dz / 3 * np.exp(1j * np.pi / 3)
            p4 = p1 + 2 * dz / 3
            z_new = np.column_stack((p1, p2, p3, p4)).reshape(-1, 1)
            z = np.concatenate([z_new, z[-1].reshape(-1, 1)])
        self.zetaMatrix = np.zeros((max(angles.shape), max(z.shape)), dtype=mp.mpc)
        for k in range(max(angles.shape)):
            self.zetaMatrix[k, :] = self.zeta0 + z.T * np.exp(1j * np.deg2rad(angles[k]))

    def mandelbrot(self: Self) -> None:
        res = 1000
        iters = 150
        x = np.linspace(-2.1, 0.6, res)
        y = np.linspace(-1.2, 1.2, res)
        X, Y = np.meshgrid(x, y)
        C = X + 1j * Y

        Zcalc = np.zeros_like(C)
        inside = np.ones_like(C, dtype=bool)
        for n in range(iters):
            mask = np.abs(Zcalc) <= 2
            Zcalc[mask] = Zcalc[mask] ** 2 + C[mask]
            inside &= np.abs(Zcalc) <= 2

        #Под вопросом
        #boundaries = sk.measure.find_contours(np.pad(inside, pad_width=1, mode='constant', constant_values=0), 0.9)
        #boundaries = [i[np.append((np.diff(i[:,0])!=0)|(np.diff(i[:,1])!=0), True)]-1 for i in [np.round(i).astype(int) for i in boundaries]]
        '''
        Так как skimage возвращает точки, ноходящиеся между 0 и 1, их приходится округлять и удалять лишние чтобы повторить поведение Matlab,
        возможно лучше использовать opencv (вероятно быстрее работает)
        '''
        #Под вопросом
        inside = inside.astype(np.uint8) * 255
        boundaries, _ = cv2.findContours(inside, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
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
        scales = np.arange(mp.mpf('1e-9'), mp.mpf('1e-9') + mp.mpf('5e-11') / 2, mp.mpf('5e-11'))
        self.zetaMatrix = np.zeros((max(z_base.shape), max(scales.shape)), dtype=mp.mpc)
        for k in range(max(scales.shape)):
            s = scales[k]
            self.zetaMatrix[:, k] = self.zeta0 + s * (z_base - self.zeta0)
        self.zetaMatrix = self.zetaMatrix.T

    def julia(self: Self, res: int = 100, maxIter: int = 20) -> None:
        C_constant = mp.mpc(-0.123, 0.745)
        x = np.linspace(-1.5, 1.5, res, dtype=mp.mpf)
        y = np.linspace(-1.5, 1.5, res, dtype=mp.mpf)
        X, Y = np.meshgrid(x, y)
        Zcalc = X + 1j * Y
        inside = np.ones_like(Zcalc, dtype=bool)
        for n in range(maxIter):
            mask = np.abs(Zcalc) <= 2
            Zcalc[mask] = Zcalc[mask] ** 2 + C_constant
            inside &= np.abs(Zcalc) <= 2
        inside = inside.astype(np.uint8) * 255
        boundaries, _ = cv2.findContours(inside, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_NONE)
        mainIdx = np.argmax([max(i.shape) for i in boundaries])
        mainBoundary = boundaries[mainIdx].reshape(-1, 2)
        mainBoundary = np.concatenate((mainBoundary, mainBoundary[0].reshape(-1, 2)))
        z_edge = x[mainBoundary[:, 0]] + 1j * y[mainBoundary[:, 1]]
        refIdx = np.argmax(z_edge.real)
        reference_point = z_edge[refIdx]
        offset = self.zeta0 - reference_point
        z_base = z_edge + offset
        scales = np.arange(1e-12, 1e-11 + 1e-13 / 2, 1e-13, dtype=mp.mpf)
        self.zetaMatrix = np.zeros((max(z_base.shape), max(scales.shape)), dtype=mp.mpc)
        for k in range(max(scales.shape)):
            s = scales[k]
            self.zetaMatrix[:, k] = self.zeta0 + s * (z_base - self.zeta0)
        self.zetaMatrix = self.zetaMatrix.T

    def circles(self: Self) -> None:
        radii = np.linspace(0.000011, 0.0012, 10).reshape(-1, 1)
        theta = np.linspace(0, 2 * np.pi, self.numPoints).reshape(1, -1)
        theta_grid, radii_grid = np.meshgrid(theta, radii)
        centers = self.zeta0 + radii
        circles_adjusted = centers + radii_grid * (np.exp(1j * theta_grid) - 1.) - radii
        self.zetaMatrix = circles_adjusted.T.conj()

    def plot_gradient(self: Self, matrix: np.ndarray, draw_arrows: bool = False, 
                     labels: Dict[str, str] = None, start_point_label: str = '') -> None:
        """
        Plot gradient visualization.
        
        Parameters:
        -----------
        matrix : np.ndarray
            Complex matrix to plot
        draw_arrows : bool, optional
            Whether to draw gradient arrows
        labels : Dict[str, str], optional
            Dictionary with 'X' and 'Y' axis labels
        start_point_label : str, optional
            Label for the starting point
        """
        if labels is None:
            labels = {'X': '', 'Y': ''}

        matrix = matrix.astype(complex)
        _, m_cols = matrix.shape
        fig, ax = plt.subplots()
        colors = plt.cm.jet(np.linspace(0, 1, m_cols))
        
        # Plot starting point
        start_x = np.real(matrix[0, 0])
        start_y = np.imag(matrix[0, 0])
        ax.plot(start_x, start_y, 'ro', markersize=10, linewidth=2, markerfacecolor='r')
        
        if start_point_label:
            ax.text(start_x, start_y, f' {start_point_label}', 
                   verticalalignment='bottom', horizontalalignment='left',
                   fontname='Times New Roman', fontsize=14)
        
        # Plot each column
        for i in range(m_cols):
            x_data = np.real(matrix[:, i])
            y_data = np.imag(matrix[:, i])
            ax.plot(x_data, y_data, color=colors[i], linewidth=1.5)
            
            if draw_arrows and len(x_data) > 1:
                dx = np.gradient(x_data)
                dy = np.gradient(y_data)
                idx = np.round(np.linspace(0, len(x_data) - 1, 10)).astype(int)
                ax.quiver(x_data[idx], y_data[idx], dx[idx], dy[idx], 
                         color=colors[i], linewidth=1.2, scale=30)
        
        # Formatting
        ax.grid(True)
        ax.set_xlabel(labels.get('X', ''), fontsize=15)
        ax.set_ylabel(labels.get('Y', ''), fontsize=15)
        ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%.12f'))
        ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.12f'))
        
        ax.set_xlabel(labels.get('X', ''), fontsize=20)
        ax.set_ylabel(labels.get('Y', ''), fontsize=20)
        
        plt.show()
    
    def discrepancy(self: Self) -> None:
        """Plot the discrepancy function."""
        zetas = self.zetaMatrix
        kappas = self.kappaMatrix
        m, _ = zetas.shape
        
        plt.figure()
        colors = plt.cm.jet(np.linspace(0, 1, m))
        
        for i in range(m):
            zeta = zetas[i, :]
            kappa = kappas[i, :]
            delta = self.delta(zeta, kappa)
            delta_plot = [float(d) for d in delta]
            plt.plot(delta_plot, color=colors[i], linewidth=1.3)
        
        plt.grid(True)
        plt.show()


def main():
    cm = ComplexMapping(1, 3)
    
    # Set parameters
    cm.zeta0 = mp.mpc(1, 2)
    cm.kappa0 = mp.mpf('.5')
    
    cm.line()
    print('line done')
    cm.calcKappa()
    print('calc kappa done')
    #cm.julia()
    
    cm.plot_gradient(cm.zetaMatrix, draw_arrows=True, labels={'X': 'Real', 'Y': 'Imag'})
    #cm.discrepancy()

if __name__ == '__main__':
    main()
