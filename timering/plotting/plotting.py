import matplotlib.pyplot as plt
import numpy as np


def gaussian(x, a, mu, sig):
    return a * np.exp(-((x - mu) ** 2) / (2 * sig**2))


def result_plot(result_hdu, *, title=None):
    """
    Simple function to quickly display Zn2 results
    """
    plt.style.use("dark_background")
    fig, a = plt.subplots()
    a.set_title(title)
    a.set_xlabel(r"$\nu$ (Hz)")
    a.set_ylabel(r"$Z_1^2$")
    data = result_hdu.data
    metadata = result_hdu.header
    hm = 0.5 * metadata["Zn2"]
    fwhm = np.linspace(
        metadata["nu"] - metadata["error"], metadata["nu"] + metadata["error"]
    )
    gfit = gaussian(
            data["Frequency"],
            metadata["G_Zn2"],
            metadata["G_nu"],
            metadata["G_error"],
        )
    hm_range = hm * np.ones(len(fwhm))
    a.plot(fwhm, hm_range, label="fwhm", linestyle="dotted")
    a.plot(data["Frequency"], data["Z12"], label=r"$Z_n^2$ Result", color="cyan")
    a.plot(data["Frequency"], gfit, label="Gaussian Fit", linestyle="dashed")
    plt.legend()
    return a


def spin_curve_plot(plc_hdu, *, title=None):
    hist = plc_hdu.data["Phase Rate"]
    plt.style.use("dark_background")
    fig, a = plt.subplots()
    phase = np.linspace(0, 1, num=len(hist))
    a.plot(phase, hist)
    a.set_title(title)
    a.set_xlabel(r"Phase ($\phi$)")
    a.set_ylabel("Arrival Times")
    return a
