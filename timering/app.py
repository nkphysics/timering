import streamlit as st
import pandas as pd
import argparse
import pathlib
import plotly.express as px
import sqlite3
import numpy as np
from astropy.time import Time


p = argparse.ArgumentParser(description="Set local database")

p.add_argument("--database",
               "-db",
               default="",
               type=str,
               help="Local path to .twdb"
               )

pargs = p.parse_args()


def querymission(mission):
    mtable = table_in.loc[table_in["MISSION"] == mission.upper()]
    return mtable


def filter_intmode(table, intmode):
    if intmode != "all":
        table = table.loc[table["MODE"] == intmode.lower()]
    return table


def filtermax(table, column, maxbound):
    table = table.loc[table[column] <= float(maxbound)]
    return table


def filtermin(table, column, minbound):
    table = table.loc[table[column] >= float(minbound)]
    return table


if "show_df" not in st.session_state:
    st.session_state.show_df = "Off"

with st.sidebar:
    st.markdown("# Timering")
    dbpath = st.text_input("Local database path", value=pargs.database)
    dbpath = pathlib.Path(dbpath).resolve()
    plottype = st.radio("Timing Evolution Plot Type", ["Line", "Scatter"])
    show_df = st.radio(
        "Show Data Table",
        ["Off", "On"],
    )

con = sqlite3.connect(dbpath)
table_in = pd.read_sql_query("SELECT * FROM nu_results", con,
                             parse_dates=["TIME"])
table_in["TIME"] = pd.to_datetime(table_in["TIME"])
table_in = table_in.sort_values(by="TIME")
with st.sidebar:
    with st.expander("Timing Evolution Filters"):
        st.markdown("Interval Modes:")
        set1, set2 = st.columns(2)
        with set1:
            xteintmode = st.selectbox("XTE",
                                      options=["all", "full", "gti"]
                                      )
        with set2:
            nicerintmode = st.selectbox("NICER",
                                        options=["all", "full", "gti"]
                                        )
        nitable = querymission("NICER")
        nitable = filter_intmode(nitable, nicerintmode)
        xtetable = querymission("XTE")
        xtetable = filter_intmode(xtetable, xteintmode)
        table_in = pd.merge(xtetable, nitable, how="outer")
        st.markdown(r"$\nu$ Error (hz)")
        set3, set4 = st.columns(2)
        with set3:
            minerr = st.text_input("Min Error", value=0.0)
        with set4:
            maxerr = st.text_input("Max Error", value=1.0)
        table_in = filtermax(table_in, "NU_ERR", maxerr)
        table_in = filtermin(table_in, "NU_ERR", minerr)
        minzn2 = st.text_input("Min $Z_n^2$", value=30.0)
        table_in = filtermin(table_in, "ZN2", minzn2)
        st.markdown(r"$\nu$ Gaussian Fit Error (hz)")
        set5, set6 = st.columns(2)
        with set5:
            mingerr = st.text_input("Min Gaussian Error", value=-0.1)
        with set6:
            maxgerr = st.text_input("Max Gaussian Error", value=1.0)
        table_in = filtermax(table_in, "G_NU_ERR", maxgerr)
        table_in = filtermin(table_in, "G_NU_ERR", mingerr)
        st.markdown(r"Exposure (s)")
        set7, set8 = st.columns(2)
        with set7:
            minexpo = st.text_input("Min Exposure", value=0.0)
        with set8:
            maxexpo = st.text_input("Max Exposure", value=1000000)
        table_in = filtermax(table_in, "EXPOSURE", maxexpo)
        table_in = filtermin(table_in, "EXPOSURE", minexpo)
        st.markdown(r"Arrival Times (counts)")
        set9, set10 = st.columns(2)
        with set9:
            minats = st.text_input("Min ATs", value=0.0)
        with set10:
            maxats = st.text_input("Max ATs", value=1000000000)
        table_in = filtermax(table_in, "ATS", maxats)
        table_in = filtermin(table_in, "ATS", minats)
st.session_state.show_df = show_df


def tevo_plot(table, plottype):
    fig = 0
    if plottype == "Scatter":
        fig = px.scatter(
            table,
            x="TIME",
            y="NU",
            title="Timing Evolution",
            color="MISSION",
            error_y="NU_ERR",
        )
    else:
        fig = px.line(
            table,
            x="TIME",
            y="NU",
            title="Timing Evolution",
            color="MISSION",
            error_y="NU_ERR",
            markers=True,
        )
    fig.update_layout(xaxis_title="Time", yaxis_title=r"Spin Frequency (hz)")
    fig.update_xaxes(showgrid=True)
    fig.update_yaxes(showgrid=True)
    return fig


st.plotly_chart(tevo_plot(table_in, plottype))

total, tnicer, txte = st.columns(3)
total.metric("Total", len(table_in["NU"]))
tnicer.metric("NICER", len(table_in.loc[table_in["MISSION"] == "NICER"]))
txte.metric("XTE", len(table_in.loc[table_in["MISSION"] == "XTE"]))
if st.session_state.show_df == "On":
    crab_table = st.dataframe(table_in)


def boxcarfit(df, tpts=1, lpts=1, order=1):
    """
    Performs a boxcar fit by least squares.
    Requires time and corresponding spin dataframe
    Returns boxcar fitting coefficients
    """
    df = df.reset_index(drop=True)
    coeffs = {"TIME": [], "DNU": []}
    for num, _ in enumerate(df["TIME"], start=tpts):
        try:
            times = []
            spins = []
            for j in range(1, tpts):
                times.append(df["TIME"][num - j])
                spins.append(df["NU"][num - j])
            for j in range(1, lpts):
                times.append(df["TIME"][num + j])
                spins.append(df["NU"][num + j])
            times = Time(times, format="datetime").mjd
            twoorder = np.polynomial.polynomial.polyfit(times, spins, order)
            coeffs["TIME"].append(df["TIME"][num])
            coeffs["DNU"].append(twoorder[1])
        except KeyError:
            break
    return coeffs


with st.sidebar:
    with st.expander(r"$\delta \nu$ Fitting"):
        trail = st.slider("# of Trailing Boxcar Points", 2, 10, 2)
        lead = st.slider("# of Leading Boxcar Points", 2, 10, 2)
        order = st.number_input("Order of Polynomial", min_value=1,
                                max_value=3, step=1)

nuresiduals = boxcarfit(table_in, tpts=trail, lpts=lead, order=order)
nures_plot = px.line(
            x=nuresiduals["TIME"],
            y=nuresiduals["DNU"],
            title="Nu Residuals",
            markers=True,
        )
nures_plot.update_layout(xaxis_title="Time",
                         yaxis_title=r"Nu Residual")
st.plotly_chart(nures_plot, order=order)
