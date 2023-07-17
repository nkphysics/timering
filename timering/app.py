import streamlit as st
import pandas as pd
import argparse
import pathlib
import plotly.graph_objects as go
import sqlite3
import numpy as np
from astropy.time import Time
import plotting
from astroquery.heasarc import Heasarc


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


def getcrabtime():
    """
    Retrieves the CRABTIME database (Crab Pulsar Monthly
    Ephemeris from Jodrell Bank) from HEASARC
    """
    heasarc = Heasarc()
    crabtime = heasarc.query_object(object_name='PSR B0531+21', mission='crabtime')
    jbcrab = crabtime.to_pandas()
    jbcrab.rename(columns={"MJD": "TIME", "NU_ERROR": "NU_ERR"}, inplace=True)
    jbcrab["TIME"] = Time(jbcrab["TIME"], format="mjd").to_datetime()
    jbcrab.drop('RA', axis=1, inplace=True)
    jbcrab.drop('DEC', axis=1, inplace=True)
    jbcrab.drop('SEARCH_OFFSET_', axis=1, inplace=True)
    return jbcrab


if "show_df" not in st.session_state:
    st.session_state.show_df = "Off"

with st.sidebar:
    st.markdown("# Timering")
    dbpath = st.text_input("Local database path", value=pargs.database)
    dbpath = pathlib.Path(dbpath).resolve()
    plottype = st.radio("Nu Evolution Plot Type", ["Line", "Scatter"])
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
st.markdown(r"## $\nu$ Evolution")

st.plotly_chart(plotting.evo_plot(table_in, plottype))

total, tnicer, txte = st.columns(3)
total.metric("Total", len(table_in["NU"]))
tnicer.metric("NICER", len(table_in.loc[table_in["MISSION"] == "NICER"]))
txte.metric("XTE", len(table_in.loc[table_in["MISSION"] == "XTE"]))
if st.session_state.show_df == "On":
    crab_table = st.dataframe(table_in, use_container_width=True)
    st.markdown(r"""
                **NOTE:**

                Data presneted in this table rounds to $10^{-4}$ Hz by default.
                However, in most cases measurements of $\nu$ are obtained
                at resolutions of $<10^{-5}$ Hz. Selecting single cells of the
                dataframe will display the full values.
                """)


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
            time0 = Time(df["TIME"][num], format="datetime").mjd
            polycoes = np.polynomial.polynomial.polyfit(times, spins, order)
            poly = np.poly1d(polycoes)
            coeffs["TIME"].append(df["TIME"][num])
            prediction = poly(time0)
            diff = df["NU"][num] - prediction
            coeffs["DNU"].append(diff)
        except KeyError:
            break
    return coeffs


with st.sidebar:
    with st.expander(r"$\delta \nu$ Fitting from $\nu$"):
        show_nufit = st.radio("Show Plot",
                              ["Off", "On"])
        set1, set2 = st.columns(2)
        with set1:
            addcrabtime = st.radio("Include CRABTIME",
                                   ["Off", "On"])
        with set2:
            addxray = st.radio("Include X-Ray",
                               ["Off", "On"],
                               index=1)
        st.session_state.show_nufit = show_nufit
        trail = st.slider("# of Trailing Boxcar Points", 2, 10, 2)
        lead = st.slider("# of Leading Boxcar Points", 2, 10, 2)
        order = st.number_input("Order of Polynomial", min_value=1,
                                max_value=3, step=1)

if "show_nufit" not in st.session_state:
    st.session_state.show_nufit = "Off"
if st.session_state.show_nufit == "On":
    st.divider()
    st.markdown(r"## $\delta \nu$ Fitting of $\nu$ Data")
    dnuplot = go.Figure()
    if addcrabtime == "On":
        jbcrab = getcrabtime()
        radiodnu = boxcarfit(jbcrab, tpts=trail, lpts=lead, order=order)
        dnuplot.add_trace(go.Scatter(x=radiodnu["TIME"],
                                     y=radiodnu["DNU"],
                                     name='Radio'))
    if addxray == "On":
        nuresiduals = boxcarfit(table_in, tpts=trail, lpts=lead, order=order)
        dnuplot.add_trace(go.Scatter(x=nuresiduals["TIME"],
                                     y=nuresiduals["DNU"],
                                     name='X-Ray'))
    dnuplot.update_layout(xaxis_title="Time",
                          yaxis_title=r"Nu Residual")
    st.plotly_chart(dnuplot, order=order)


nif = pd.read_sql("SELECT NICER.OBSID, NICER.TWR_FILE FROM NICER " +
                  "WHERE NICER.TWR_FILE IS NOT NULL", con)
xtef = pd.read_sql("SELECT XTE.OBSID, XTE.TWR_FILE FROM XTE " +
                   "WHERE XTE.TWR_FILE IS NOT NULL", con)
resdf = pd.merge(nif, xtef, how="outer")
resdf = pd.merge(resdf, table_in, how="inner", on="OBSID")
resdf = resdf.drop_duplicates(subset=["OBSID"])
st.divider()

st.markdown("# Individual Measurements")
obsid = st.selectbox("OBSID", resdf["OBSID"])
obsid_filt = resdf.loc[resdf["OBSID"] == obsid]
obsid_filt.reset_index(drop=True, inplace=True)
rplots = plotting.obsid_plots(obsid_filt["TWR_FILE"][0])
for num, i in enumerate(rplots):
    if num % 2 == 0:
        if num <= 1:
            st.markdown("### Full")
        else:
            subindex = int(num - (num / 2.0))
            st.markdown(f"### Interval {num - subindex}")
        sluicing, phasecurve = st.columns(2)
        with sluicing:
            st.pyplot(rplots[num])
        with phasecurve:
            st.pyplot(rplots[num + 1])
with st.expander("More Info On Individual Results"):
    st.markdown("""
                Individual Results are processed in two different modes, FULL and GTI.
                It will be best to start with GTI mode

                ### GTI Mode

                GTI mode may be misleading from the name.
                Results processed in GTI mode are measurements optained from a singular exposure
                interval during an observation. This is much more applicable to missions like
                NICER that group multiple exposure intervals into a single .evt datafile.
                To explain as simply as possible, missions in orbit do not get continueous
                exposure from a source due to the missions' orbit, as well as the orbits of other
                objects that may impeded the view of a source by the mission. This is dealt with
                by not collecting, or clearing out data when the source is not in view. This
                leads to periods across an observation of for example 1200s of exposure followed
                by >2000s of nothing before another exposure. In some cases this short duration
                between exposures can be lead to noticable differences in measurements of spin
                Therefore each interval is isolated and measured on its own.

                ### Full Mode

                Full mode measurements are obtained by taking all photon arrival times from
                an observation and using them to obtain a single measurement of spin. This is
                done disregarding total observation exposure or the number of exposure intervals.
                For sources that are easy to obtain measurements from, this can lead to noisy
                measurement results. But, for sources in which it is harder to obtain timing
                measurements, the full mode can be useful in obtaining something as opposed to
                nothing.

                DO NOTE:

                Observations or .evt files with only one exposure interval will only be processed
                in full mode. This is very much so the case with XTE results since XTE .evt files
                generally come in the form of one exposure interval per .evt file.
                """)
