"""
Created on 2025/11/28
Author: Quirijn B. van Woerkom
Implementations and validation of some basic rheologies.
"""

import numpy as np
from scipy.special import gamma


class Element:
    """Class representing a rheological element."""

    def complex_compliance(
        self,
        mu_elastic: np.ndarray[float],
        eta: np.ndarray[float],
        omega: np.ndarray[float],
        **kwargs,
    ) -> np.ndarray[float]:
        """Compute the complex compliance of this element.

        This method should be overridden by subclasses to provide
        specific element models. Make sure that this is vectorised.

        Parameters
        ----------
        mu_elastic : np.ndarray[float]
            Elastic shear modulus.
        eta : np.ndarray[float]
            Viscosity. Must have the same shape as mu_elastic.
        omega : np.ndarray[float]
            Angular frequency. Need not have the same shape as
            mu_elastic.
        **kwargs : dict
            Element-specific parameters that are not used in this model.

        Returns
        -------
        J_complex : np.ndarray[float]
            Complex compliance at the specified conditions and
            frequencies.
        """
        raise NotImplementedError(
            "Subclasses of Element must implement the "
            "complex_compliance method."
        )

class ElasticElement(Element):
    """Pure elastic Hookean solid."""
    @staticmethod
    def complex_compliance(
        mu_elastic, eta, omega, **kwargs
    ):
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.where(
                mu_elastic == 0,
                np.inf,
                1.0 / mu_elastic
            )

class MaxwellElement(Element):
    """Rheological element based on the Maxwell model."""

    @staticmethod
    def complex_compliance(
        mu_elastic: np.ndarray[float],
        eta: np.ndarray[float],
        omega: np.ndarray[float],
        **kwargs,
    ) -> np.ndarray[float]:
        """Compute the complex compliance using the Maxwell model.

        Parameters
        ----------
        mu_elastic : np.ndarray[float]
            Elastic shear modulus.
        eta : np.ndarray[float]
            Viscosity. Must have the same shape as mu_elastic.
        omega : np.ndarray[float]
            Angular frequency. Need not have the same shape as
            mu_elastic.
        **kwargs : dict
            Element-specific parameters that are not used in this model.

        Returns
        -------
        J_complex : np.ndarray[float]
            Complex compliance at the specified conditions and
            frequencies.
        """
        # Wherever mu_elastic, omega, or eta are zero, return infinity
        with np.errstate(divide="ignore", invalid="ignore"):
            J_complex = np.where(
                (mu_elastic == 0) | (omega == 0) | (eta == 0),
                np.inf,
                1 / mu_elastic - 1j / (omega * eta),
            )
        return J_complex


class AndradeElement(Element):
    """Rheological element based on the Andrade (zeta) model."""

    def __init__(self, zeta: float, n: float) -> None:
        """Initialise an Andrade element.

        Initialise an Andrade element with the specified parameters.
        For the nomenclature, we follow Bierson (2024).

        Parameters
        ----------
        zeta : float
            Andrade proportionality parameter zeta.
        n : float
            Andrade exponent n.
        """
        self.zeta = zeta
        self.n = n
        self.gamma_val = gamma(n + 1)

    def complex_compliance(
        self,
        mu_elastic: np.ndarray[float],
        eta: np.ndarray[float],
        omega: np.ndarray[float],
        **kwargs,
    ) -> np.ndarray[float]:
        """Compute the complex compliance using the Andrade model.

        Parameters
        ----------
        mu_elastic : np.ndarray[float]
            Elastic shear modulus.
        eta : np.ndarray[float]
            Viscosity. Must have the same shape as mu_elastic.
        omega : np.ndarray[float]
            Angular frequency. Need not have the same shape as
            mu_elastic.
        **kwargs : dict
            Element-specific parameters that are not used in this model.

        Returns
        -------
        J_complex : np.ndarray[float]
            Complex compliance at the specified conditions and
            frequencies.
        """
        # Wherever mu_elastic is zero, return infinity
        with np.errstate(divide="ignore", invalid="ignore"):
            J_complex = np.where(
                mu_elastic == 0,
                np.inf,
                self.gamma_val
                / mu_elastic
                * (self.zeta * eta / mu_elastic * omega * 1j) ** (-self.n),
            )
        return J_complex


class KelvinVoigtElement(Element):
    """Rheological element based on the Kelvin-Voigt model."""

    @staticmethod
    def complex_compliance(
        mu_elastic: np.ndarray[float],
        eta: np.ndarray[float],
        omega: np.ndarray[float],
        mu_KV: np.ndarray[float],
        eta_KV: np.ndarray[float],
        **kwargs,
    ) -> np.ndarray[float]:
        """Compute the complex compliance using the Kelvin-Voigt model.

        Parameters
        ----------
        mu_elastic : np.ndarray[float]
            Elastic shear modulus.
        eta : np.ndarray[float]
            Viscosity. Must have the same shape as mu_elastic.
        omega : np.ndarray[float]
            Angular frequency. Need not have the same shape as
            mu_elastic.
        mu_KV : np.ndarray[float]
            Kelvin-Voigt spring modulus. Must have the same shape as
            mu_elastic.
        eta_KV : np.ndarray[float]
            Kelvin-Voigt viscosity. Must have the same shape as
            mu_elastic.
        **kwargs : dict
            Element-specific parameters that are not used in this model.

        Returns
        -------
        J_complex : np.ndarray[float]
            Complex compliance at the specified conditions and
            frequencies.
        """
        # Wherever mu_elastic, omega, and eta are zero, return infinity
        with np.errstate(divide="ignore", invalid="ignore"):
            J_complex = np.where(
                (mu_elastic == 0) & (omega == 0) & (eta == 0),
                np.inf,
                (
                    (mu_KV - 1j * omega * eta_KV)
                    / (mu_KV**2 + (omega * eta_KV) ** 2)
                ),
            )
        return J_complex


# TODO: add Cole rheological element


class Rheology:
    """Parent class representing a rheological law for a Layer."""

    def __init__(self, elements: list[Element]) -> None:
        """Initialise a rheology.

        Parameters
        ----------
        elements : list[str]
            List of elements that this rheology is comprised of. The
            complex compliances is computed as the sum of the complex
            compliances of each element.
        """
        # Save elements
        self.elements = elements

    def __call__(
        self,
        mu_elastic: np.ndarray[float],
        eta: np.ndarray[float],
        omega: np.ndarray[float],
        **kwargs,
    ) -> np.ndarray[float]:
        """Compute the complex shear modulus.

        Parameters
        ----------
        mu_elastic : np.ndarray[float]
            Elastic shear modulus.
        eta : np.ndarray[float]
            Viscosity. Must have the same shape as mu_elastic.
        omega : np.ndarray[float]
            Angular frequency. Need not have the same shape as
            mu_elastic.
        **kwargs : dict
            Rheology-specific parameters.

        Returns
        -------
        mu_complex : np.ndarray[float]
            Complex shear modulus at the specified conditions and
            frequencies.
        """
        # Validate input shapes
        # Test whether mu_elastic and eta have the same shape
        # or whether they are scalars
        if np.isscalar(mu_elastic) or np.isscalar(eta):
            pass
        elif mu_elastic.shape != eta.shape:
            raise ValueError("mu_elastic and eta must have the same shape.")
        return 1 / self.complex_compliance(mu_elastic, eta, omega, **kwargs)

    def complex_compliance(
        self,
        mu_elastic: np.ndarray[float],
        eta: np.ndarray[float],
        omega: np.ndarray[float],
        **kwargs,
    ) -> np.ndarray[float]:
        """Implementation of the complex compliance calculation.

        This method should be overridden by subclasses to provide
        specific rheological models. Make sure that this is vectorised.

        Parameters
        ----------
        mu_elastic : np.ndarray[float]
            Elastic shear modulus.
        eta : np.ndarray[float]
            Viscosity. Must have the same shape as mu_elastic.
        omega : np.ndarray[float]
            Angular frequency. Need not have the same shape as
            mu_elastic.
        **kwargs : dict
            Rheology-specific parameters.

        Returns
        -------
        J_complex : np.ndarray[float]
            Complex compliance at the specified conditions and
            frequencies.
        """
        return sum(
            element.complex_compliance(mu_elastic, eta, omega, **kwargs)
            for element in self.elements
        )


class MaxwellRheology(Rheology):
    """Rheology based on a single Maxwell element."""

    def __init__(self) -> None:
        """Initialise a Maxwell rheology."""
        super().__init__([MaxwellElement()])


class AndradeRheology(Rheology):
    """Rheology based on a single Andrade element."""

    def __init__(self, zeta: float, n: float) -> None:
        """Initialise an Andrade rheology.

        Parameters
        ----------
        zeta : float
            Andrade proportionality parameter zeta.
        n : float
            Andrade exponent n.
        """
        super().__init__([MaxwellElement(), AndradeElement(zeta, n)])


class BurgersRheology(Rheology):
    """Rheology based on a Maxwell and Kelvin-Voigt element in series."""

    def __init__(
        self,
    ) -> None:
        """Initialise a Burgers rheology."""
        super().__init__([MaxwellElement(), KelvinVoigtElement()])


class SundbergCooperRheology(Rheology):
    """Rheology with Maxwell, Andrade and Kelvin-Voigt elements."""

    def __init__(self, zeta: float, n: float) -> None:
        """Initialise a Sundberg-Cooper rheology.

        Parameters
        ----------
        zeta : float
            Andrade proportionality parameter zeta.
        n : float
            Andrade exponent n.
        """
        super().__init__(
            [
                MaxwellElement(),
                AndradeElement(zeta, n),
                KelvinVoigtElement(),
            ]
        )

class ElasticRheology(Rheology):
    """Purely elastic rheology."""
    def __init__(self):
        super().__init__([ElasticElement()])
class LiquidRheology(Rheology):
    """Rheology representing a fluid with no shear resistance."""
    def __init__(self):
        super().__init__([])
    def __call__(
        self,
        mu_elastic=None,
        eta=None,
        omega=None,
        **kwargs,
    ):
        """Return zero complex shear modulus (no resistance to shear)."""
        return 0.0

if __name__ == "__main__":
    # Validate against Tobie et al. (2025) Fig. 3
    # Shared parameters
    n1 = 0.2
    n2 = 0.3
    zeta = 1
    omega_values = np.logspace(0, -10, 500)  # rad/s
    P_values = 2 * np.pi / omega_values / (24 * 3600)  # days

    # Ices
    mu_elastic_ice = 3.3e9  # Pa
    eta_ice = 1e14  # Pa s
    mu_KV_ice = 3.3e9  # Pa
    eta_KV_ice = 60e9  # Pa s

    # Silicates
    mu_elastic_sil = 70e9  # Pa
    eta_sil = 1e20  # Pa s
    mu_KV_sil = 70e9  # Pa
    eta_KV_sil = 7e16  # Pa s

    # Create rheologies
    maxwell = MaxwellRheology()
    andrade1 = AndradeRheology(zeta, n1)
    andrade2 = AndradeRheology(zeta, n2)
    burgers = BurgersRheology()
    sundberg_cooper1 = SundbergCooperRheology(zeta, n1)
    sundberg_cooper2 = SundbergCooperRheology(zeta, n2)

    maxwell_ice = maxwell(
        mu_elastic_ice,
        eta_ice,
        omega_values,
    )
    andrade1_ice = andrade1(
        mu_elastic_ice,
        eta_ice,
        omega_values,
    )
    andrade2_ice = andrade2(
        mu_elastic_ice,
        eta_ice,
        omega_values,
    )
    burgers_ice = burgers(
        mu_elastic_ice,
        eta_ice,
        omega_values,
        mu_KV=mu_KV_ice,
        eta_KV=eta_KV_ice,
    )
    sundberg_cooper1_ice = sundberg_cooper1(
        mu_elastic_ice,
        eta_ice,
        omega_values,
        mu_KV=mu_KV_ice,
        eta_KV=eta_KV_ice,
    )
    sundberg_cooper2_ice = sundberg_cooper2(
        mu_elastic_ice,
        eta_ice,
        omega_values,
        mu_KV=mu_KV_ice,
        eta_KV=eta_KV_ice,
    )
    maxwell_sil = maxwell(
        mu_elastic_sil,
        eta_sil,
        omega_values,
    )
    andrade1_sil = andrade1(
        mu_elastic_sil,
        eta_sil,
        omega_values,
    )
    andrade2_sil = andrade2(
        mu_elastic_sil,
        eta_sil,
        omega_values,
    )
    burgers_sil = burgers(
        mu_elastic_sil,
        eta_sil,
        omega_values,
        mu_KV=mu_KV_sil,
        eta_KV=eta_KV_sil,
    )
    sundberg_cooper1_sil = sundberg_cooper1(
        mu_elastic_sil,
        eta_sil,
        omega_values,
        mu_KV=mu_KV_sil,
        eta_KV=eta_KV_sil,
    )
    sundberg_cooper2_sil = sundberg_cooper2(
        mu_elastic_sil,
        eta_sil,
        omega_values,
        mu_KV=mu_KV_sil,
        eta_KV=eta_KV_sil,
    )
    # Plot real and imaginary parts of the complex shear modulus
    import matplotlib.pyplot as plt

    # First do ices
    plt.figure(figsize=(10, 6))
    plt.plot(
        P_values,
        maxwell_ice.real,
        label="Maxwell (Ice)",
        color="blue",
        linestyle="--",
    )
    plt.plot(
        P_values,
        andrade1_ice.real,
        label="Andrade n=0.2 (Ice)",
        color="blue",
        linestyle="-",
    )
    plt.plot(
        P_values,
        andrade2_ice.real,
        label="Andrade n=0.3 (Ice)",
        color="cyan",
        linestyle="-",
    )
    plt.plot(
        P_values,
        burgers_ice.real,
        label="Burgers (Ice)",
        color="navy",
        linestyle=":",
    )
    plt.plot(
        P_values,
        sundberg_cooper1_ice.real,
        label="Sundberg-Cooper n=0.2 (Ice)",
        color="blue",
        linestyle="-.",
    )
    plt.plot(
        P_values,
        sundberg_cooper2_ice.real,
        label="Sundberg-Cooper n=0.3 (Ice)",
        color="cyan",
        linestyle="-.",
    )
    plt.xscale("log")
    plt.xlabel("Period (days)")
    plt.ylabel("Real Part of Complex Shear Modulus (Pa)")
    plt.title("Complex Shear Modulus - Real Part")
    plt.xlim(1e-4, 1e5)
    plt.legend()
    plt.grid(True, which="both", ls="--")
    plt.show()
    # Then do silicates in a separate plot
    plt.figure(figsize=(10, 6))
    plt.plot(
        P_values,
        maxwell_sil.real,
        label="Maxwell (Silicate)",
        color="red",
        linestyle="--",
    )
    plt.plot(
        P_values,
        andrade1_sil.real,
        label="Andrade n=0.2 (Silicate)",
        color="red",
        linestyle="-",
    )
    plt.plot(
        P_values,
        andrade2_sil.real,
        label="Andrade n=0.3 (Silicate)",
        color="orange",
        linestyle="-",
    )
    plt.plot(
        P_values,
        burgers_sil.real,
        label="Burgers (Silicate)",
        color="darkred",
        linestyle=":",
    )
    plt.plot(
        P_values,
        sundberg_cooper1_sil.real,
        label="Sundberg-Cooper n=0.2 (Silicate)",
        color="red",
        linestyle="-.",
    )
    plt.plot(
        P_values,
        sundberg_cooper2_sil.real,
        label="Sundberg-Cooper n=0.3 (Silicate)",
        color="orange",
        linestyle="-.",
    )
    plt.xscale("log")
    plt.xlabel("Period (days)")
    plt.ylabel("Real Part of Complex Shear Modulus (Pa)")
    plt.title("Complex Shear Modulus - Real Part")
    plt.xlim(1e-4, 1e5)
    plt.legend()
    plt.grid(True, which="both", ls="--")
    plt.show()

    # Now do imaginary parts
    # First do ices
    plt.figure(figsize=(10, 6))
    plt.plot(
        P_values,
        maxwell_ice.imag,
        label="Maxwell (Ice)",
        color="blue",
        linestyle="--",
    )
    plt.plot(
        P_values,
        andrade1_ice.imag,
        label="Andrade n=0.2 (Ice)",
        color="blue",
        linestyle="-",
    )
    plt.plot(
        P_values,
        andrade2_ice.imag,
        label="Andrade n=0.3 (Ice)",
        color="cyan",
        linestyle="-",
    )
    plt.plot(
        P_values,
        burgers_ice.imag,
        label="Burgers (Ice)",
        color="navy",
        linestyle=":",
    )
    plt.plot(
        P_values,
        sundberg_cooper1_ice.imag,
        label="Sundberg-Cooper n=0.2 (Ice)",
        color="blue",
        linestyle="-.",
    )
    plt.plot(
        P_values,
        sundberg_cooper2_ice.imag,
        label="Sundberg-Cooper n=0.3 (Ice)",
        color="cyan",
        linestyle="-.",
    )
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Period (days)")
    plt.ylabel("Imaginary Part of Complex Shear Modulus (Pa)")
    plt.title("Complex Shear Modulus - Imaginary Part")
    plt.xlim(1e-4, 1e5)
    plt.legend()
    plt.grid(True, which="both", ls="--")
    plt.show()
    # Then do silicates in a separate plot
    plt.figure(figsize=(10, 6))
    plt.plot(
        P_values,
        maxwell_sil.imag,
        label="Maxwell (Silicate)",
        color="red",
        linestyle="--",
    )
    plt.plot(
        P_values,
        andrade1_sil.imag,
        label="Andrade n=0.2 (Silicate)",
        color="red",
        linestyle="-",
    )
    plt.plot(
        P_values,
        andrade2_sil.imag,
        label="Andrade n=0.3 (Silicate)",
        color="orange",
        linestyle="-",
    )
    plt.plot(
        P_values,
        burgers_sil.imag,
        label="Burgers (Silicate)",
        color="darkred",
        linestyle=":",
    )
    plt.plot(
        P_values,
        sundberg_cooper1_sil.imag,
        label="Sundberg-Cooper n=0.2 (Silicate)",
        color="red",
        linestyle="-.",
    )
    plt.plot(
        P_values,
        sundberg_cooper2_sil.imag,
        label="Sundberg-Cooper n=0.3 (Silicate)",
        color="orange",
        linestyle="-.",
    )
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Period (days)")
    plt.ylabel("Imaginary Part of Complex Shear Modulus (Pa)")
    plt.title("Complex Shear Modulus - Imaginary Part")
    plt.xlim(1e-4, 1e5)
    plt.legend()
    plt.grid(True, which="both", ls="--")
    plt.show()
