import matplotlib.pyplot as plt
import numpy as np
import plotly.express as px
from astropy.io import fits
import pathlib


def gaussian(x, a, mu, sig):
    return a * np.exp(-((x - mu) ** 2) / (2 * sig**2))


def resultplot(result_hdu, *, title=None):
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
    a.plot(data["Frequency"], data["Z12"],
           label=r"$Z_n^2$ Result",
           color="cyan")
    a.plot(data["Frequency"], gfit, label="Gaussian Fit", linestyle="dashed")
    plt.legend()
    return fig


def spincurveplot(plc_hdu, *, title=None):
    hist = plc_hdu.data["Phase Rate"]
    plt.style.use("dark_background")
    fig, a = plt.subplots()
    phase = np.linspace(0, 1, num=len(hist))
    a.plot(phase, hist)
    a.set_title(title)
    a.set_xlabel(r"Phase ($\phi$)")
    a.set_ylabel("Arrival Times")
    return fig


def evo_plot(table, plottype):
    fig = 0
    if plottype == "Scatter":
        fig = px.scatter(
            table,
            x="TIME",
            y="NU",
            color="MISSION",
            error_y="NU_ERR",
        )
    else:
        fig = px.line(
            table,
            x="TIME",
            y="NU",
            color="MISSION",
            error_y="NU_ERR",
            markers=True,
        )
    fig.update_layout(xaxis_title="Time", yaxis_title=r"Spin Frequency (hz)")
    fig.update_xaxes(showgrid=True)
    fig.update_yaxes(showgrid=True)
    return fig


class Result:
    def __init__(self, resultfile):
        self.resultfile = pathlib.Path(resultfile).resolve()
        self.hdul = self.open_twr()
        self.phdr = self.get_pheader()

    def open_twr(self):
        return fits.open(self.resultfile)

    def get_hdu(self, index):
        hdu = self.hdul[index]
        return hdu.header, hdu

    def get_pheader(self):
        phdr = self.hdul[0].header
        return phdr

    def get_title(self, header):
        return (f"{self.phdr['MISSION']} " +
                f"{self.phdr['OBSID']} @ " +
                f"{header['NU_DT']}")

    def close_file(self):
        return self.hdul.close()


def obsid_plots(resultfile):
    """
    Plots all measurement results and phase curves
    from a specfied result file
    """
    obsid = Result(resultfile)
    plots = []
    for i in range(1, len(obsid.hdul)):
        header, table = obsid.get_hdu(i)
        title = obsid.get_title(header)
        if header["EXTNAME"] == "PHASE CURVE":
            plot = spincurveplot(table, title=title)
            plots.append(plot)
        else:
            plot = resultplot(obsid.hdul[i], title=title)
            plots.append(plot)
    return plots
