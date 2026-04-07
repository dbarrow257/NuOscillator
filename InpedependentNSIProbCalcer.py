
# --------------------------------------------------------------------------------
# Caroline Fengler, 06.01.2026
# Independent NSI oscillation probability calculator for crosscheck with NuOscillator
# --------------------------------------------------------------------------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm
from matplotlib.colors import ListedColormap

# ----------------------------
# Unit conversions
# ----------------------------
GEV_TO_EV = 1.0e9
KM_TO_EVINV = 5.067730716e9   # 1 km = 5.067730716e9 eV^-1
V_COEFF = 7.63247e-14         # eV * (g/cm^3)^-1 * Ye^-1


# Flavor indices: 0=e, 1=mu, 2=tau

def pmns_matrix(theta12, theta13, theta23, delta_cp):
    """Standard PDG PMNS matrix."""
    s12, c12 = np.sin(theta12), np.cos(theta12)
    s13, c13 = np.sin(theta13), np.cos(theta13)
    s23, c23 = np.sin(theta23), np.cos(theta23)
    ed = np.exp(-1j * delta_cp)
    ep = np.exp(+1j * delta_cp)

    U = np.array([
        [ c12*c13,                          s12*c13,              s13*ed ],
        [-s12*c23 - c12*s23*s13*ep,   c12*c23 - s12*s23*s13*ep,   s23*c13 ],
        [ s12*s23 - c12*c23*s13*ep,  -c12*s23 - s12*c23*s13*ep,   c23*c13 ]
    ], dtype=complex)
    return U


def matter_nsi_matrix(eps_ee=0.0, eps_emu=0.0, eps_etau=0.0,
                      eps_mumu=0.0, eps_mutau=0.0, eps_tautau=0.0):
    """Full 3x3 NSI matrix epsilon."""
    return np.array([
        [eps_ee,            eps_emu,             eps_etau],
        [np.conj(eps_emu),  eps_mumu,            eps_mutau],
        [np.conj(eps_etau), np.conj(eps_mutau),  eps_tautau]
    ], dtype=complex)


def hamiltonian_const_density(
    E_GeV,
    rho_gcm3,
    Ye,
    dm21,
    dm31,
    theta12,
    theta13,
    theta23,
    delta_cp,
    eps=None,
    antineutrino=False
):
    """
    Return the flavor-basis Hamiltonian in eV for constant density matter.
    """
    E_eV = E_GeV * GEV_TO_EV
    U = pmns_matrix(theta12, theta13, theta23, delta_cp)

    M2 = np.diag([0.0, dm21, dm31])
    H_vac = U @ M2 @ U.conj().T / (2.0 * E_eV)

    V = V_COEFF * rho_gcm3 * Ye

    if eps is None:
        eps = np.zeros((3, 3), dtype=complex)

    M_matter = np.array([[1.0, 0.0, 0.0],
                         [0.0, 0.0, 0.0],
                         [0.0, 0.0, 0.0]], dtype=complex) + eps

    if antineutrino:
        # U -> U*, V -> -V, epsilon -> epsilon*
        H = H_vac.conj() - V * M_matter.conj()
    else:
        H = H_vac + V * M_matter

    return H


def probability_mumu_layers(
    E_GeV,
    cosz,
    layers,
    depth_km=1.5,
    rho_crust=2.7,
    **kwargs
    # Ye=0.5,
    # dm21=7.65e-5,
    # dm31=2.5e-3,
    # theta12=0.56, # np.deg2rad(33.44),
    # theta13=0.148, # np.deg2rad(8.57),
    # theta23=0.79, # np.deg2rad(49.2),
    # delta_cp=4.19, # np.deg2rad(195.0),
    # eps=None,
    # antineutrino=False
):
    """
    P(nu_mu -> nu_mu) in constant-density matter with NSI.
    """
    if cosz >= 0:
        # Down-going: only the short overburden above the detector matters
        L_km = baseline_from_cosz(cosz, depth_km=depth_km)
        H = hamiltonian_const_density(E_GeV, rho_crust, kwargs.get("Ye", 0.5),
                                      kwargs.get("dm21", 7.65e-5),
                                      kwargs.get("dm31", 2.5e-3),
                                      kwargs.get("theta12", 0.56),
                                      kwargs.get("theta13", 0.148),
                                      kwargs.get("theta23", 0.79),
                                      kwargs.get("delta_cp", 4.19),
                                      eps=kwargs.get("eps", None),
                                      antineutrino=kwargs.get("antineutrino", False))
        S = expm(-1j * H * (L_km * KM_TO_EVINV))
        return float(np.abs(S[1, 1])**2)

    # Up-going: use the shell model
    S = np.eye(3, dtype=complex)

    for (R_outer, R_inner, rho) in layers:
        L_i = path_length_in_shell(cosz, R_outer, R_inner, depth_km=depth_km)
        if L_i <= 0:
            continue

        H_i = hamiltonian_const_density(
            E_GeV, rho, kwargs.get("Ye", 0.5),
            kwargs.get("dm21", 7.42e-5),
            kwargs.get("dm31", 2.517e-3),
            kwargs.get("theta12", np.deg2rad(33.44)),
            kwargs.get("theta13", np.deg2rad(8.57)),
            kwargs.get("theta23", np.deg2rad(49.2)),
            kwargs.get("delta_cp", np.deg2rad(195.0)),
            eps=kwargs.get("eps", None),
            antineutrino=kwargs.get("antineutrino", False)
        )

        S_i = expm(-1j * H_i * (L_i * KM_TO_EVINV))
        S = S_i @ S

    return float(np.abs(S[1, 1])**2)

#########################################################################
###########################################################################

R_EARTH_KM = 6371.0

def baseline_from_cosz(cosz, depth_km=1.5, R_km=R_EARTH_KM):
    cosz = np.clip(cosz, -1.0, 1.0)
    r = R_km - depth_km
    return -r * cosz + np.sqrt(R_km**2 - r**2 * (1.0 - cosz**2))


def distance_to_radius(cosz, R, depth_km=0.0):
    """
    Distance along the neutrino direction from the detector to the point
    where the trajectory intersects a sphere of radius R.

    cosz is signed:
      cosz > 0  -> downward-going
      cosz < 0  -> upward-going
    """
    r_det = R_EARTH_KM - depth_km
    disc = R**2 - r_det**2 * (1.0 - cosz**2)
    if disc < 0:
        return 0.0
    return -r_det * cosz + np.sqrt(disc)


def path_length_in_shell(cosz, R_outer, R_inner, depth_km=0.0):
    """
    Path length through the shell [R_inner, R_outer]
    along the ray defined by cosz.
    """
    s_outer = distance_to_radius(cosz, R_outer, depth_km=depth_km)
    s_inner = distance_to_radius(cosz, R_inner, depth_km=depth_km)
    return max(0.0, s_outer - s_inner)

#########################################################################
###########################################################################

# New Colormap
newcolors = ['#FF0500','#FF1E00','#FF3800','#FF5100','#FF6B00','#FF8400','#FF9E00','#FFB700','#FFD100','#FFEA00','#F4F40A','#C1C13D','#8E8E70','#5B5BA3','#2828D6','#0000F9','#0000E0','#0000C6','#0000AD','#000093']
newcolors = newcolors[::-1]
cmap = ListedColormap(newcolors)


# Earth layers simple
layers = [
    (6371.0, 5701.0, 3.3),   # crust
    (5701.0, 3480.0, 5.0),   # mantle
    (3480.0, 1220.0, 11.3),  # outer core
    (1220.0, 0.0,    13.0),  # inner core
]

# Example NSI parameters
eps = matter_nsi_matrix(
    eps_ee=0.0,
    eps_emu=0.0,
    eps_etau=0.0,
    eps_mumu=0.05,
    eps_mutau=0.0, # + 0.00j,
    eps_tautau=0.25
)

E_vals = np.logspace(-1, 2, 600)      # GeV
cosz_vals = np.linspace(-1.0, 1.0, 601)

P = np.zeros((len(cosz_vals), len(E_vals)))

for i, cz in enumerate(cosz_vals):
    L_km = baseline_from_cosz(cz, depth_km=1.5)
    rho = 2.7
    for j, E in enumerate(E_vals):
        P[i, j] = probability_mumu_layers(
            E_GeV=E,
            cosz=cz,
            layers=layers,
            depth_km=1.5,
            eps=eps,
            antineutrino=False
        )

fig, ax = plt.subplots(figsize=(6, 5))
pcm = ax.pcolormesh(E_vals, cosz_vals, P, shading="auto", cmap=cmap)
ax.set_xscale("log")
ax.set_xlabel("Energy [GeV]")
ax.set_ylabel(r"$\cos\theta_z$")
cbar = fig.colorbar(pcm, ax=ax)
cbar.set_label("Probability")
plt.tight_layout()
plt.show()

