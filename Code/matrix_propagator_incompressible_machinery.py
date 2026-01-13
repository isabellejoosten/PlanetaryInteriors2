"""
Created on 2025/11/08
Author: Marc Rovira-Navarro
Machinery to compute tidal response using matrix propagator method for incompressible bodies.
It defines a TidalInterior that is made of TidalLayers with constant properties.
"""
from rheologies import Rheology, MaxwellRheology, LiquidRheology, ElasticRheology
from math import pi
from scipy.constants import G
from typing import Optional, List, Tuple, Callable
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import time
import numpy as np
from numpy import pi, logspace, zeros
from matplotlib.ticker import ScalarFormatter


class TidalLayer:
    """
    Represents a layer in the tidal model with constant material properties.
    """
    def __init__(self,
                 name: str,
                 thickness: Optional[float] = None,
                 inner_radius: Optional[float] = None,
                 outer_radius: Optional[float] = None,
                 inner_gravity: Optional[float] = None,
                 outer_gravity: Optional[float] = None,
                 density: float = None,
                 viscosity: Optional[float] = None,
                 shear_modulus: Optional[float] = None,
                 rheology: Optional[Rheology] = None,
                 liquid: bool = False):
        """
        Parameters
        ----------
        name : str
            Name of the layer.
        thickness : float, optional
            Layer thickness. If not provided, it is inferred from inner/outer radius.
        inner_radius : float, optional
            Radius at the inner boundary. Computed if thickness and outer_radius are given.
        outer_radius : float, optional
            Radius at the outer boundary. Computed if thickness and inner_radius are given.
        inner_gravity : float, optional
            Gravitational acceleration at the inner boundary (computed if not provided).
        outer_gravity : float, optional
            Gravitational acceleration at the outer boundary (computed if not provided).
        density : float
            Constant density of the layer (kg/m³).
        viscosity: float, optional
            Viscosity in Pa·s.
        shear_modulus : float, optional
            Elastic shear modulus in Pa.
        rheology : Rheology, optional
            Rheological model. Defaults to MaxwellRheology if not provided.
        liquid : bool, optional
            Whether the layer behaves as a liquid.
        """
        self.name = name    # layer name
        self.inner_radius= inner_radius
        self.outer_radius= outer_radius
        # compute thickness if inner and outer radius are provided
        if self.inner_radius is not None and self.outer_radius is not None:
            self.thickness = outer_radius - inner_radius
        else:
            self.thickness = thickness  # m
        self.inner_radius = inner_radius  # m
        self.outer_radius = outer_radius  # m
        self.inner_gravity = inner_gravity  # m/s^2
        self.outer_gravity = outer_gravity
        self.density = density      # kg/m^3
        self.eta=viscosity  # Pa.s
        self.mu=shear_modulus              # Pa
        self.liquid=liquid  # bool
        self.rheology=rheology
        if rheology is None:
            self.rheology = MaxwellRheology()
        else:
            self.rheology = rheology
    def __repr__(self):
        thickness_km = f"{self.thickness / 1e3:.2f}" if self.thickness is not None else "None"
        eta_str = f"{self.eta:.2e}" if self.eta is not None else "None"
        mu_str = f"{self.mu:.2e}" if self.mu is not None else "None"
        density_str = f"{self.density:.1f}" if self.density is not None else "None"
        rheology_name = type(self.rheology).__name__ if self.rheology is not None else "None"
        return (
            f"TidalLayer(name={self.name}, "
            f"thickness={thickness_km} km, "
            f"density={density_str} kg/m^3, "
            f"viscosity={eta_str} Pa.s, "
            f"shear_modulus={mu_str} Pa, "
            f"rheology={rheology_name})"
        )
    def mass_gravity(self, inner_radius:float, inner_mass:float) -> float:
        """Calculates the mass and gravity of the layer given the inner radius of the layer and the inner mass
        Parameters
        ----------
        inner_radius : float
            inner_radius of the layer (m)
        inner_mass : float
            The mass enclosed within the inner radius (kg)
        """
        self.inner_radius = inner_radius
        self.outer_radius = self.inner_radius + self.thickness
        mass = (4/3) * pi* (self.outer_radius**3 - self.inner_radius**3)* self.density
        self.mass=mass
        self.inner_gravity= G * inner_mass / (self.inner_radius**2) if inner_radius > 0 else 0
        self.outer_gravity= G * (inner_mass + mass) / (self.outer_radius**2)
        return mass
class TidalInterior:
    """Represents the interior structure for tidal computations"""
    def __init__(self, name: str, layers: List[TidalLayer]):
        self.name = name
        self.layers = layers
        self.radius = None
        self.mass = None
        # Compute mass and gravty profile if not defined
        if self.layers[-1].outer_gravity is None:
            self.compute_mass_and_gravity()
        # normalize parameters
        self.normalize_parameters()
        self.rho_avg=self.mass/(4/3*pi*self.radius**3)
        self.surf_gravity=G*self.mass/(self.radius**2)
    def total_radius(self) -> float:
        """Calculates the total radius of the planet based on its layers."""
        return sum(layer.thickness for layer in self.layers)
    def compute_mass_and_gravity(self):
        """Computes mass of each layer, total mass, and gravity."""
        inner_radius = 0
        inner_mass = 0
        for layer in self.layers:
            layer_mass = layer.mass_gravity(inner_radius, inner_mass)
            inner_mass += layer_mass
            inner_radius += layer.thickness
        self.mass = inner_mass
        self.radius = self.total_radius()
        self.gravity=self.layers[-1].outer_gravity

    def normalize_parameters(self, period: float = 24 * 3600):
        """
        Normalizes shear modulus, density, radii, thickness, viscosity, and gravity.
        Stores them as attributes on each Layer.
        Also stores the normalized gravitational constant and reference values as attributes.

        Parameters
        ----------
        period : float, optional
            Tidal forcing period in seconds. Default is 1 day (24*3600 s).
        """
        # Reference parameters
        self.mu0 = self.layers[-1].mu  # reference shear modulus (last layer)
        self.rho0 = self.layers[-1].density  # reference density (outermost layer)
        self.R0 = self.total_radius()  # total planetary radius
        self.T0 = period  # reference tidal period

        if self.mu0 is None:
            raise ValueError("Last layer must have a defined shear modulus.")
        if self.rho0 is None:
            raise ValueError("Last layer must have a defined density.")

        for i, layer in enumerate(self.layers):
            # check that layers have required properties
            if not isinstance(layer.rheology,LiquidRheology) and layer.mu is None:
                raise ValueError(f"Layer {layer.name} must have a shear modulus defined for solid rheologies.")
            if not (isinstance(layer.rheology,LiquidRheology) or isinstance(layer.rheology,ElasticRheology)) and layer.eta is None:
                raise ValueError(f"Layer {layer.name} must have a viscosity defined for viscoelastic rheologies.")
            if isinstance(layer.rheology,LiquidRheology):
                if i == 0:  # liquid core
                    layer.liquid=True
                else:
                    #another liquid layers (subsurface ocean), set shear modulus to very low value
                    layer.mu=1e-8*self.layers[-1].mu
                    layer.rheology=ElasticRheology()
            # Shear modulus
            layer.mu_norm = layer.mu / self.mu0 if layer.mu is not None else None
            # Density
            layer.density_norm = layer.density / self.rho0 if layer.density is not None else None
            # Inner and outer radius
            layer.inner_radius_norm = layer.inner_radius / self.R0
            layer.outer_radius_norm = layer.outer_radius / self.R0
            # Thickness
            layer.thickness_norm = layer.thickness / self.R0
            # Dimensionless viscosity using period
            layer.eta_norm = (
                        2 * np.pi * layer.eta / (self.mu0 * self.T0)) if layer.eta is not None else None
            # Normalized gravity
            layer.inner_gravity_norm = layer.inner_gravity * self.rho0 * self.R0 / self.mu0
            layer.outer_gravity_norm = layer.outer_gravity * self.rho0 * self.R0 / self.mu0
        # Normalized gravitational constant
        self.G_norm = G * (self.rho0 ** 2 * self.R0 ** 2 / self.mu0)
        self.radius_norm = self.radius / self.R0
    def __repr__(self):
        layer_strings = "\n  ".join(repr(layer) for layer in self.layers)
        return (
            f"Interior_Model(\n"
            f"  name={self.name!r},\n"
            f"  radius={self.radius:.3e},\n"
            f"  mass={self.mass:.3e},\n"
            f"  layers=[\n"
            f"  {layer_strings}\n"
            f"  ]\n"
            f")"
        )

class TidalResponse:
    """
    Compute the tidal response of a body, including the Love number :math:`k_2`,
    radial displacement :math:`y(r)`, and dissipation function :math:`H(r)`,
    for a given interior model and tidal frequency.
    """

    def __init__(self, interior_model: TidalInterior, omega: float):
        """
        Parameters
        ----------
        interior_model : TidalInterior
            A tidal interior model defining the layered structure of the body.
            It must provide access to layer radii, density, viscosity, shear
            modulus, gravitational parameters, and normalization constants
            (`G_norm`, `T0`).
        omega : float
            Tidal forcing frequency in rad/s.
        """
        self.interior_model = interior_model
        self.omega = omega
        self.G = interior_model.G_norm
        self.omega_norm = omega * interior_model.T0 / (2 * pi)
        self._prepare_layers()          # compute effective rigidity and initialize arrays
        self._compute_layer_matrices()  # compute layer propagation matrices and k2
        self._setup_y_H_functions()     # setup vectorized y(r) and H(r) functions

    def __repr__(self):
        return (
            f"TidalResponse(\n"
            f"  model={self.interior_model!r},\n"
            f"  omega={self.omega:.3e} rad/s,\n"
            f"  k2={self.k2.item().real:.5e} + {self.k2.item().imag:.5e}j\n"
            f"  h2={self.h2.item().real:.5e} + {self.h2.item().imag:.5e}j\n"
            f")"
        )
    # ----------------------
    # Internal setup methods
    # ----------------------
    def _prepare_layers(self):
        """Compute effective rigidity for each layer and initialize arrays."""
        for layer in self.interior_model.layers:
            layer.mu_eff_norm =layer.rheology(layer.mu_norm,layer.eta_norm, self.omega_norm)

        # Initialize arrays for vectorized computation
        self.num_layers = len(self.interior_model.layers)
        self.layer_bottom_matrices_array = np.zeros((self.num_layers, 6, 3), dtype=np.complex128)
        self.layer_inner_radius_array = np.zeros(self.num_layers)
        self.layer_outer_radius_array = np.zeros(self.num_layers)
        self.layer_mu_eff_array = np.zeros(self.num_layers, dtype=np.complex128)
        self.layer_density_array = np.zeros(self.num_layers, dtype=np.complex128)
        self.layer_inner_gravity_array = np.zeros(self.num_layers, dtype=np.complex128)

        for i, layer in enumerate(self.interior_model.layers):
            self.layer_inner_radius_array[i] = layer.inner_radius_norm
            self.layer_outer_radius_array[i] = layer.outer_radius_norm
            self.layer_mu_eff_array[i] = layer.mu_eff_norm if layer.mu_eff_norm is not None else 0.0
            self.layer_density_array[i] = layer.density_norm
            self.layer_inner_gravity_array[i] = layer.inner_gravity_norm

    def _compute_layer_matrices(self):
        """Compute propagation matrices from center to surface and obtain k2"""
        self.layer_bottom_matrices = {}
        layers = self.interior_model.layers

        # Boundary conditions at the center
        if layers[0].liquid:
            r_core = layers[0].outer_radius_norm
            g_core = layers[0].outer_gravity_norm
            rho_core = layers[0].density_norm
            layer_start_indx = 1
            IC = np.array([
                [-r_core**2 / g_core, 0, 1],
                [0, 1, 0],
                [0, 0, rho_core * g_core],
                [0, 0, 0],
                [r_core**2, 0, 0],
                [2 * r_core, 0, 4 * pi * self.G * rho_core]
            ], dtype=np.complex128)
            self.layer_bottom_matrices[layers[0]] = np.zeros((6, 3), dtype=np.complex128)
        else:
            layer_start_indx = 0
            IC = np.array([
                [1, 0, 0],
                [0, 1, 0],
                [0, 0, 1],
                [0, 0, 0],
                [0, 0, 0],
                [0, 0, 0]
            ], dtype=np.complex128)

        # Surface boundary conditions
        BC = np.array([[0], [0], [5 / self.interior_model.radius_norm]], dtype=np.complex128)
        P1 = np.array([
            [0, 0, 1, 0, 0, 0],
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 0, 1]
        ], dtype=np.complex128)
        P2 = np.array([
            [1, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 0]
        ], dtype=np.complex128)

        # Propagate through layers
        Prop_matrix = IC.copy()
        first_layer = layers[layer_start_indx]
        self.layer_bottom_matrices[first_layer] = Prop_matrix.copy()

        # Special handling for first solid layer after liquid core
        if layer_start_indx == 1:
            upper = layers[1]
            Yinv_upper = np.squeeze(self.get_Yinv(
                upper.inner_radius_norm, upper.mu_eff_norm, upper.density_norm, upper.inner_gravity_norm, 2, self.G))
            Prop_matrix = Yinv_upper @ Prop_matrix
            self.layer_bottom_matrices[upper] = Prop_matrix.copy()

        # Loop through layers 2..N-1
        for i in range(layer_start_indx, len(layers) - 1):
            lower = layers[i]
            upper = layers[i + 1]
            Y_lower = np.squeeze(self.get_Y(lower.outer_radius_norm, lower.mu_eff_norm,
                                           lower.density_norm, lower.outer_gravity_norm, 2, self.G))
            Yinv_upper = np.squeeze(self.get_Yinv(upper.inner_radius_norm, upper.mu_eff_norm,
                                                 upper.density_norm, upper.inner_gravity_norm, 2, self.G))
            Prop_matrix = Yinv_upper @ Y_lower @ Prop_matrix
            self.layer_bottom_matrices[upper] = Prop_matrix.copy()

        # Last layer to surface
        last_layer = layers[-1]
        Y_last = np.squeeze(self.get_Y(last_layer.outer_radius_norm, last_layer.mu_eff_norm,
                                      last_layer.density_norm, last_layer.outer_gravity_norm, 2, self.G))
        Prop_matrix = Y_last @ Prop_matrix

        # Solve integration constants
        C = np.linalg.inv(P1 @ Prop_matrix) @ BC
        self.C = C
        self.k2 = (P2 @ Prop_matrix @ C)[-1] - 1
        self.h2= -self.interior_model.layers[-1].outer_gravity_norm*(P2 @ Prop_matrix @ C)[0]

        # Convert bottom matrices to arrays
        for i, layer in enumerate(layers):
            self.layer_bottom_matrices_array[i] = self.layer_bottom_matrices[layer]

    def _setup_y_H_functions(self):
        """Define vectorized y(r) and H(r) functions."""
        def y_func(r_norm_array: np.ndarray) -> np.ndarray:
            N = len(r_norm_array)
            layer_indices = np.searchsorted(self.layer_outer_radius_array, r_norm_array, side='left')
            Prop_to_bottom_all = self.layer_bottom_matrices_array[layer_indices]
            mu_eff_all = self.layer_mu_eff_array[layer_indices]
            rho_all = self.layer_density_array[layer_indices]
            inner_g_all = self.layer_inner_gravity_array[layer_indices]
            inner_r_all = self.layer_inner_radius_array[layer_indices]

            mass_in_layer = (4 / 3) * np.pi * (r_norm_array ** 3 - inner_r_all ** 3) * rho_all
            g_r = inner_g_all * (inner_r_all / r_norm_array) ** 2 + self.G * mass_in_layer / r_norm_array ** 2

            Y_r = self.get_Y(r_norm_array, mu_eff_all, rho_all, g_r, 2, self.G)
            YProp = np.einsum('nij,njk->nik', Y_r, Prop_to_bottom_all)
            y_all = np.einsum('nij,jk->nik', YProp, self.C)
            return y_all[:, :, 0].T

        def H_func(r_norm_array: np.ndarray) -> np.ndarray:
            y_all = y_func(r_norm_array)
            y1 = y_all[0]
            y2 = y_all[1]
            y4 = y_all[3]

            y1_d = 1 / r_norm_array * (-2 * y1 + 6 * y2)
            fA = 4 / 3 * np.abs(r_norm_array * y1_d - y1 + 3 * y2) ** 2
            fB = 6 * np.abs(r_norm_array * y4 / self.layer_mu_eff_array[np.searchsorted(
                self.layer_outer_radius_array, r_norm_array, side='left')]) ** 2
            fC = 24 * np.abs(y2) ** 2

            mu_eff_all = self.layer_mu_eff_array[np.searchsorted(self.layer_outer_radius_array, r_norm_array, side='left')]
            H_all = mu_eff_all.imag * (fA + fB + fC)
            H_all = np.where(mu_eff_all == 0.0, 0.0, H_all)
            return H_all

        self.y = y_func
        self.H = H_func

    # ----------------------
    # Public methods
    # ----------------------
    def plot(self, axes=None, ax_big=None, label=None, scale=None):
        r_grid = np.linspace(1e-4, 1, 5000)
        y_r = self.y(r_grid)
        H_r = self.H(r_grid)

        if scale is None:
            scale = np.max(np.abs(H_r))
            exponent = int(np.floor(np.log10(scale)))
        else:
            exponent = scale

        H_scaled = H_r / 10**exponent

        ax_big.plot(
            r_grid,
            H_scaled,
            label=fr"{label} ($\times 10^{{{exponent}}}$)"
        )
        ax_big.set_ylabel(r"$H / 10^{n}$")  
        layer_boundaries = [layer.outer_radius_norm for layer in self.interior_model.layers]

        titles = [r"$U_2$", r"$V_2$", r"$R_2$", r"$S_2$", r"$\phi_2$", r"$Q_2$"]
        ylabels = [r"$y_1$", r"$y_2$", r"$y_3$", r"$y_4$", r"$y_5$", r"$y_6$"]

        for i, ax in enumerate(axes):
            ax.plot(r_grid, y_r[i, :].real, label=f"{label} real")
            ax.plot(r_grid, y_r[i, :].imag, linestyle="--", label=f"{label} imag")

            for boundary in layer_boundaries:
                ax.axvline(boundary, color="red", linestyle=":", alpha=0.7)

    def plot_H_only(self, ax=None, label=None):
        r_grid = np.linspace(1e-4, 1, 5000)
        H_r = self.H(r_grid)
        layer_boundaries = [layer.outer_radius_norm for layer in self.interior_model.layers]

        # Create axis if needed
        if ax is None:
            fig, ax = plt.subplots(figsize=(6, 6))
        else:
            fig = ax.figure

        # --- scaling logic ---
        max_val = np.max(np.abs(H_r))
        if max_val == 0:
            exponent = 0
        else:
            exponent = int(np.floor(np.log10(max_val)))

        H_scaled = H_r / 10**exponent

        # Plot scaled curve
        ax.plot(
            r_grid,
            H_scaled,
            label=fr"{label} ($\times 10^{{{exponent}}}$)" if label else None
        )

        # Layer boundaries
        for boundary in layer_boundaries:
            ax.axvline(boundary, color="red", linestyle=":", alpha=0.7)

        # Axis formatting
        ax.set_xlim(0, 1)
        ax.set_xlabel(r"$r / R$")
        ax.set_ylabel(r"$H$ (scaled)")
        ax.set_title(r"Tidal Dissipation Function $H(r)$")
        ax.grid(True, linestyle=":", alpha=0.5)

        # Disable scientific notation on axis
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.ticklabel_format(style="plain", axis="y")

        if label:
            ax.legend(title="Rheology")

        return fig, ax


    def validate_dissipation(self):
        """Compare integrated H(r) against k2 imaginary part."""
        N_grid = np.logspace(1, 6, num=10)
        difference = np.zeros_like(N_grid, dtype=float)
        layer_boundaries = [1e-4] + [layer.outer_radius_norm for layer in self.interior_model.layers]
        times = np.zeros_like(N_grid, dtype=float)

        for i, N_per_layer in enumerate(N_grid):
            start_time = time.perf_counter()
            r_grid_layers = []
            for j in range(len(layer_boundaries) - 1):
                r_start = layer_boundaries[j]
                r_end = layer_boundaries[j + 1]
                r_layer = np.linspace(r_start, r_end, int(N_per_layer), endpoint=False)
                r_grid_layers.append(r_layer)
            r_grid_layers.append(np.array([1.0]))
            r_grid = np.concatenate(r_grid_layers)

            H_r = self.H(r_grid)
            integral_H = np.trapz(H_r, r_grid)
            tidal_dissipation_k2 = -5 * self.interior_model.radius_norm / (4 * pi * self.G) * self.k2.imag
            difference[i] = (integral_H - tidal_dissipation_k2) / tidal_dissipation_k2
            times[i] = time.perf_counter() - start_time

        fig, ax1 = plt.subplots(figsize=(10, 5))
        color = 'tab:blue'
        ax1.set_xlabel("Number of radial points per layer")
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_ylabel("Relative difference |(∫H - k2_diss)/k2_diss|", color=color)
        ax1.plot(N_grid, np.abs(difference), marker='o', linestyle='-', color=color)
        ax1.tick_params(axis='y', labelcolor=color)
        ax1.grid(True, which="both", linestyle="--", alpha=0.5)

        ax2 = ax1.twinx()
        color = 'tab:red'
        ax2.set_ylabel("Computation time [s]", color=color)
        ax2.plot(N_grid, times, marker='x', linestyle='--', color=color)
        ax2.tick_params(axis='y', labelcolor=color)
        ax2.set_yscale('log')
        plt.title(f"Tidal Dissipation Convergence for {self.interior_model.name}")
        plt.tight_layout()
        #plt.show()

    # ----------------------
    # Static methods
    # ----------------------
    @staticmethod
    def get_Y(
            r: np.ndarray[float],
            mu: np.ndarray[float],
            rho: np.ndarray[float],
            g: np.ndarray[float],
            l: int = 2,
            G: float = 6.67430e-11,
    ) -> np.ndarray[complex]:
        """
        Vectorized incompressible propagation matrix for arrays of radius, shear modulus, density, and gravity.
        Define matrix propagator from Sabadini 2016 Eq (2.42)
        Parameters
        ----------
        r : np.ndarray[float]
            Radius array of shape (N,)
        mu : np.ndarray[float]
            Shear modulus array of shape (N,) or scalar
        rho : np.ndarray[float]
            Density array of shape (N,) or scalar
        g : np.ndarray[float]
            Gravity array of shape (N,) or scalar
        l : int, optional
            Spherical harmonic degree, default 2
        G : float, optional
            Gravitational constant, default 6.67430e-11

        Returns
        -------
        Y : np.ndarray[complex]
            N x 6 x 6 propagation matrix
        """
        r = np.atleast_1d(r)
        N = r.size
        mu = np.broadcast_to(mu, N)
        rho = np.broadcast_to(rho, N)
        g = np.broadcast_to(g, N)

        Y = np.zeros((N, 6, 6), dtype=np.complex128)

        # Powers
        r_l = r ** l
        r_lm1 = r_l / r  # r ** (l - 1)
        r_lm2 = r_lm1 / r  # r ** (l - 2)
        r_l1 = r_l * r  # r ** (l + 1)
        r_neg_l = r ** (-l)
        r_neg_l1 = r_neg_l / r  # r ** (-l - 1)
        r_neg_l2 = r_neg_l1 / r  # r ** (-l - 2)
        r_neg_l3 = r_neg_l2 / r  # r ** (-l - 3)

        # --- Y matrix ---
        Y[:, 0, 0] = l * r_l1 / (2 * (2 * l + 3))
        Y[:, 0, 1] = r_lm1
        Y[:, 0, 3] = (l + 1) * r_neg_l / (2 * (2 * l - 1))
        Y[:, 0, 4] = r_neg_l2

        Y[:, 1, 0] = (l + 3) * r_l1 / (2 * (2 * l + 3) * (l + 1))
        Y[:, 1, 1] = r_lm1 / l
        Y[:, 1, 3] = (2 - l) * r_neg_l / (2 * l * (2 * l - 1))
        Y[:, 1, 4] = -r_neg_l2 / (l + 1)

        Y[:, 2, 0] = (l * rho * g * r + 2 * (l ** 2 - l - 3) * mu) * r_l / (2 * (2 * l + 3))
        Y[:, 2, 1] = (rho * g * r + 2 * (l - 1) * mu) * r_lm2
        Y[:, 2, 2] = rho * r_l
        Y[:, 2, 3] = ((l + 1) * rho * g * r - 2 * (l ** 2 + 3 * l - 1) * mu) / (2 * (2 * l - 1) * r_l1)
        Y[:, 2, 4] = (rho * g * r - 2 * (l + 2) * mu) * r_neg_l3
        Y[:, 2, 5] = rho / r_l1

        Y[:, 3, 0] = l * (l + 2) * mu * r_l / ((2 * l + 3) * (l + 1))
        Y[:, 3, 1] = 2 * (l - 1) * mu * r_lm2 / l
        Y[:, 3, 3] = (l ** 2 - 1) * mu / (l * (2 * l - 1) * r_l1)
        Y[:, 3, 4] = 2 * (l + 2) * mu / ((l + 1)) * r_neg_l3

        Y[:, 4, 2] = r_l
        Y[:, 4, 5] = r_neg_l1

        Y[:, 5, 0] = 2 * np.pi * G * rho * l * r_l1 / (2 * l + 3)
        Y[:, 5, 1] = 4 * np.pi * G * rho * r_lm1
        Y[:, 5, 2] = (2 * l + 1) * r ** (l - 1)
        Y[:, 5, 3] = 2 * np.pi * G * rho * (l + 1) / (2 * l - 1) * r_neg_l
        Y[:, 5, 4] = 4 * np.pi * G * rho / r ** (l + 2)

        return Y

    @staticmethod
    def get_Yinv(
            r: np.ndarray[float],
            mu: np.ndarray[float],
            rho: np.ndarray[float],
            g: np.ndarray[float],
            l: int = 2,
            G: float = 6.67430e-11,
    ) -> np.ndarray[complex]:
        """
        Vectorized incompressible inverse propagation matrix for arrays of radius, shear modulus, density, and gravity.

        Parameters
        ----------
        r : np.ndarray[float]
            Radius array of shape (N,)
        mu : np.ndarray[float]
            Shear modulus array of shape (N,) or scalar
        rho : np.ndarray[float]
            Density array of shape (N,) or scalar
        g : np.ndarray[float]
            Gravity array of shape (N,) or scalar
        l : int, optional
            Spherical harmonic degree, default 2
        G : float, optional
            Gravitational constant, default 6.67430e-11

        Returns
        -------
        Yinv : np.ndarray[complex]
            N x 6 x 6 inverse propagation matrix
        """
        r = np.atleast_1d(r)
        N = r.size
        mu = np.broadcast_to(mu, N)
        rho = np.broadcast_to(rho, N)
        g = np.broadcast_to(g, N)

        c = 1.0 / (2 * l + 1)

        # Powers
        r_l = r ** l
        r_lm1 = r_l / r  # r ** (l - 1)
        r_l1 = r_l * r  # r ** (l + 1)
        r_l2 = r_l1 * r  # r ** (l + 2)

        # --- diagonal row scalars d_j ---
        d0 = c * (l + 1) / r_l1
        d1 = c * l * (l + 1) / (2 * (2 * l - 1) * r_lm1)
        d2 = c * (-1.0 / r_lm1)
        d3 = c * (l * r_l)
        d4 = c * (l * (l + 1) / (2 * (2 * l + 3)) * r_l2)
        d5 = c * r_l1

        Yinv = np.zeros((N, 6, 6), dtype=np.complex128)

        # -------------------------
        # Row 0
        # -------------------------
        Yinv[:, 0, 0] = d0 * (rho * g * r / mu - 2 * (l + 2))
        Yinv[:, 0, 1] = d0 * (2 * l * (l + 2))
        Yinv[:, 0, 2] = d0 * (-r / mu)
        Yinv[:, 0, 3] = d0 * (l * r / mu)
        Yinv[:, 0, 4] = d0 * (rho * r / mu)

        # -------------------------
        # Row 1
        # -------------------------
        Yinv[:, 1, 0] = d1 * (-rho * g * r / mu + 2 * (l ** 2 + 3 * l - 1) / (l + 1))
        Yinv[:, 1, 1] = d1 * (-2 * (l ** 2 - 1))
        Yinv[:, 1, 2] = d1 * (r / mu)
        Yinv[:, 1, 3] = d1 * ((2 - l) * r / mu)
        Yinv[:, 1, 4] = d1 * (-rho * r / mu)

        # -------------------------
        # Row 2
        # -------------------------
        Yinv[:, 2, 0] = d2 * (4 * np.pi * G * rho)
        Yinv[:, 2, 5] = d2 * (-1)

        # -------------------------
        # Row 3
        # -------------------------
        Yinv[:, 3, 0] = d3 * (rho * g * r / mu + 2 * (l - 1))
        Yinv[:, 3, 1] = d3 * (2 * (l ** 2 - 1))
        Yinv[:, 3, 2] = d3 * (-r / mu)
        Yinv[:, 3, 3] = d3 * (-(l + 1) * r / mu)
        Yinv[:, 3, 4] = d3 * (rho * r / mu)

        # -------------------------
        # Row 4
        # -------------------------
        Yinv[:, 4, 0] = d4 * (-rho * g * r / mu - 2 * (l ** 2 - l - 3) / l)
        Yinv[:, 4, 1] = d4 * (-2 * l * (l + 2))
        Yinv[:, 4, 2] = d4 * (r / mu)
        Yinv[:, 4, 3] = d4 * ((l + 3) * r / mu)
        Yinv[:, 4, 4] = d4 * (-rho * r / mu)

        # -------------------------
        # Row 5
        # -------------------------
        Yinv[:, 5, 0] = d5 * (4 * np.pi * G * rho * r)
        Yinv[:, 5, 4] = d5 * (2 * l + 1)
        Yinv[:, 5, 5] = d5 * (-r)
        return Yinv

    # Examples of usage
    class examples:
        # ----------------------------------------------------------
        # (2) Mantle Dominated Io
        # ----------------------------------------------------------
        @staticmethod
        def mantle_dominated_io():
            from numpy import pi
            layers = [
               TidalLayer("Core", thickness=980e3, density=5150, rheology=LiquidRheology()),
               TidalLayer("Mantle", thickness=811e3, density=3260, viscosity=2.8e15, shear_modulus=1e10),
               TidalLayer("Lithosphere", thickness=30e3, density=3260, viscosity=1e23, shear_modulus=6.5e10),
            ]
            omega = 2 * pi / (1.769 * 24 * 3600)
            model = TidalInterior("Io Mantle Dominated", layers)
            tide = TidalResponse(model, omega)
            print(tide)
            tide.plot()
            tide.validate_dissipation()
            return model, tide

        # ----------------------------------------------------------
        # (3) Asthenosphere Dominated Io
        # ----------------------------------------------------------
        @staticmethod
        def asthenosphere_dominated_io():
            from numpy import pi
            layers = [
               TidalLayer("Core", thickness=965e3, density=5150, rheology=LiquidRheology()),
               TidalLayer("Lower Mantle", thickness=626.6e3, density=3244, shear_modulus=6e10, viscosity=1e20),
               TidalLayer("Asthenosphere", thickness=200e3, density=3244, shear_modulus=7.8e5, viscosity=1e11),
               TidalLayer("Lithosphere", thickness=30e3, density=3244, shear_modulus=6.5e10, viscosity=1e23)
            ]
            omega = 2 * pi / (1.769 * 24 * 3600)
            model = TidalInterior("Io Asthenosphere Dominated", layers)
            tide = TidalResponse(model, omega)
            tide.plot()
            tide.validate_dissipation()
            print(tide)
            return model, tide

        # ----------------------------------------------------------
        # (4) Titan
        # ----------------------------------------------------------
        @staticmethod
        def titan():
            layers = [
               TidalLayer("Silicate Mantle", thickness=1.984e6, density=2650, shear_modulus=70e9, viscosity=1e20),
               TidalLayer("High Pressure Ice", thickness=241e3, density=1351, shear_modulus=4.8e9, viscosity=1e14),
               TidalLayer("Subsurface Ocean", thickness=250e3, density=1281.3, rheology=LiquidRheology()),
               TidalLayer("Convective Ice", thickness=10e3, density=931, shear_modulus=3.3e9, viscosity=1e14),
               TidalLayer("Conductive Ice", thickness=90e3, density=931, shear_modulus=3.3e9, viscosity=1e20),
            ]
            omega = 4.56e-6
            model = TidalInterior("Titan", layers)
            tide = TidalResponse(model, omega)
            tide.plot()
            print(tide)
            return model, tide

        # ----------------------------------------------------------
        # (5) Europa
        # ----------------------------------------------------------
        @staticmethod
        def europa():
            layers = [
               TidalLayer("Core", thickness=587e3, density=5500, rheology=LiquidRheology()),
               TidalLayer("Silicate Mantle", thickness=840e3, density=3500, shear_modulus=70e9, viscosity=1e20),
               TidalLayer("Ocean", thickness=114e3, density=1000, rheology=LiquidRheology()),
               TidalLayer("Convective Ice", thickness=12e3, density=920, shear_modulus=3.3e9, viscosity=1e14),
               TidalLayer("Conductive Ice", thickness=8e3, density=920, shear_modulus=3.3e9, viscosity=1e20),
            ]
            omega = 2.047e-5
            model = TidalInterior("Europa", layers)
            tide = TidalResponse(model, omega)
            tide.plot()
            print(tide)
            return model, tide

 # Run all examples if simply called
#if __name__ == "__main__":
#    TidalResponse.examples.mantle_dominated_io()
#    TidalResponse.examples.asthenosphere_dominated_io()
#    TidalResponse.examples.titan()
#    TidalResponse.examples.europa()

