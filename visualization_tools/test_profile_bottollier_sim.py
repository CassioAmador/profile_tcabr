# -*- coding: utf-8 -*-
"""Test simulated data with the Bottollier method"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate, optimize, interpolate
import sys
sys.path.insert(0, './../src/')

import proc_profile_bottollier as ppb


def density(r, n0=1.65e19, alpha=1.2):
    """Perfil radial de densidade."""
    return n0 * np.power((1 - (r / .18) ** 2), alpha)


def r2f_plasma(r, n0=1.65e19, alpha=1.2):
    """Perfil radial da frequência de plasma."""
    return 8.978663 * np.sqrt(density(r, n0, alpha))


def f_plasma2ne(f):
    """Converte frequência de plasma na densidade associada."""
    return (f / 8.978663) ** 2


def f_plasma2rc(f, n0=1.65e19, alpha=1.2):
    """Calcula a posicao radial associada a frequência de plasma especificado."""
#    Solução analítica para a densidade parabólica
    return np.sqrt(1 - np.power((f / 8.978663) ** 2 / (n0), 1 / alpha)) * 0.18

#    Solução genérica para qualquer formato de densidade
#    (desabilitado porque nem sempre consegue convergir)
    if np.size(f) == 1:
        fun = lambda r: 1e-9 * density(r, n0, alpha) - 1e-9 * f_plasma2ne(f)
        rc = optimize.fsolve(fun, 0.12)
    else:
        rc = np.zeros(np.size(f))
        for i in range(np.size(rc)):
            fun = lambda r: 1e-9 * \
                density(r, n0, alpha) - 1e-9 * f_plasma2ne(f[i])
            rc[i] = optimize.fsolve(fun, 0.12)
    return rc


def n_index(r, f,  n0=1.65e19, alpha=1.2):
    """Calcula o indece de refração para o plasma na posição r"""
    return np.sqrt(f**2 - r2f_plasma(r, n0, alpha)**2) / f


def phase_shift(fc, n0=1.65e19, alpha=1.2):
    phi = np.zeros(len(fc))
    for i in range(len(phi)):
        phi[i] = (4. * np.pi * fc[i] / 3e8) * integrate.quad(n_index, f_plasma2rc(fc[i],
                                                                                  n0, alpha), 0.18, args=(fc[i], n0, alpha,), epsabs=1e-14)[0] - np.pi / 2
    return phi


def v_group_inv(r, f, n0=1.65e19, alpha=1.2):
    """Calcula o inverso da velocidade de grupo para o modo O."""
    return (f * 1e-9 / np.sqrt((f * 1e-9) ** 2 - (r2f_plasma(r, n0, alpha) * 1e-9) ** 2)) / 3e8


def group_delay(f_probe, n0=1.65e19, alpha=1.2):
    """Calcula o atraso de grupo para a frequencia de sondagem."""
    rc = f_plasma2rc(f_probe, n0, alpha)
    tau = np.zeros(len(f_probe))
    for i in range(len(tau)):
        tau[i] = 2. * integrate.quad(v_group_inv, rc[i], 0.18, args=(
            f_probe[i], n0, alpha,), epsrel=1e-14, epsabs=1e-14)[0]
    return tau


if __name__ == "__main__":
    n0 = 1.65e19
    alpha = 1.2
    # frequência de sondagem experimental
    f_probe = np.linspace(16e9, np.min([35e9, r2f_plasma(0.001)]), 100)
    tau = group_delay(f_probe, n0, alpha)
    phi = phase_shift(f_probe, n0, alpha)
    
    # inicializacao linear
    # f_probe = np.append(np.linspace(1e9, f_probe[0], num=16, endpoint=False), f_probe)
    # phi = np.append(np.linspace(-np.pi/2, phi[0], num=16, endpoint=False), phi)

    # inicializacao cubica
    # phi = np.append(np.polyval([(phi[0] + np.pi) / f_probe[1] ** 3., 0, 0, -np.pi/2], np.linspace(1e9, f_probe[0], num=10, endpoint=False)), phi)
    # f_probe = np.append(np.linspace(1e9, f_probe[0], num=10, endpoint=False), f_probe)
    phi = np.append(interpolate.interp1d([0, 1e9, f_probe[0], f_probe[5]], [-np.pi / 2, -0.999 * np.pi / 2, phi[
                    0], phi[5]], kind='cubic')(np.linspace(1e9, f_probe[0], num=20, endpoint=False)), phi)
    f_probe = np.append(np.linspace(
        1e9, f_probe[0], num=20, endpoint=False), f_probe)

    # plt.plot(f_probe * 1e-9, phase_shift(f_probe, n0, alpha), 'r--')
    # plt.plot(f_probe * 1e-9, phi, '.')
    # plt.ylabel('phase [rad]')
    # plt.xlabel('probing frequency [GHz]')
    # plt.show()

    r_real = np.linspace(0, 0.18, 50)
    ne = density(r_real, n0, alpha)
    rc_BC = ppb.find_pos(f_probe * 1e-9, phi)
    ne = f_plasma2ne(f_probe)
    fig = plt.figure()
    plt.plot(rc_BC, ne * 1e-19, 'b.', label='BC')
    plt.plot(r_real[::-1], density(r_real, n0, alpha)
             * 1e-19, 'r--', label='real')
    plt.ylabel('ne [10^19 /m^3]')
    plt.xlabel('posicao radial [m]')
    plt.xlim([0, 0.18])
    plt.legend(loc=2)
    plt.show()
