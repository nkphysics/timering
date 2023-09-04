import streamlit as st
import pandas as pd
import argparse
import pathlib
import plotly.graph_objects as go
import sqlite3
import numpy as np
from astropy.time import Time
import plotting
import messages
from astroquery.heasarc import Heasarc
import logging
import sys
import os

TIMERING = os.path.basename(sys.argv[0])


def parseargs():
    p = argparse.ArgumentParser(description="Set local database")

    p.add_argument("--database",
                   "-db",
                   default="",
                   type=str,
                   help="Local path to .twdb")
    p.add_argument("--debug",
                   default=False,
                   action="store_true",
                   help="Sets logging mode to debug")

    pargs = p.parse_args()
    return pargs


def querymission(table_in, mission):
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
    crabtime = heasarc.query_object(object_name='PSR B0531+21',
                                    mission='crabtime')
    jbcrab = crabtime.to_pandas()
    jbcrab.rename(columns={"MJD": "TIME", "NU_ERROR": "NU_ERR"}, inplace=True)
    jbcrab["TIME"] = Time(jbcrab["TIME"], format="mjd").to_datetime()
    jbcrab.drop('RA', axis=1, inplace=True)
    jbcrab.drop('DEC', axis=1, inplace=True)
    jbcrab.drop('SEARCH_OFFSET_', axis=1, inplace=True)
    return jbcrab


def boxcarfit(df, tpts=1, lpts=1, order=1,
              mode="Difference"):
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
            coeffs["TIME"].append(df["TIME"][num])
            if mode == "Difference":
                poly = np.poly1d(polycoes)
                prediction = poly(time0)
                diff = df["NU"][num] - prediction
                coeffs["DNU"].append(diff)
            elif mode == "Coefficient":
                coeffs["DNU"].append(polycoes[1])
        except KeyError:
            break
    return coeffs


def main(pargs: argparse.Namespace):
    level = logging.WARNING
    if pargs.debug is True:
        level = logging.DEBUG
    logging.basicConfig(stream=sys.stdout,
                        level=level)
    logger = logging.getLogger(TIMERING)

    if "show_df" not in st.session_state:
        st.session_state.show_df = "Off"

    with st.sidebar:
        st.markdown("# Timering")
        dbpath = st.text_input("Local database path", value=pargs.database)
        dbpath = pathlib.Path(dbpath).resolve()
        show_df = st.radio(
            "Show Data Table",
            ["Off", "On"],
        )

    con = sqlite3.connect(dbpath)
    src = pd.read_sql_query("SELECT Source FROM df_metadata WHERE rowid = 1",
                            con)
    st.markdown(f"# {src['Source'][0]}")
    table_in = pd.read_sql_query("SELECT * FROM nu_results", con,
                                 parse_dates=["TIME"])
    table_in = table_in.sort_values(by="TIME")
    with st.sidebar:
        with st.expander("Timing Evolution Filters"):
            set0, set00 = st.columns(2)
            with set0:
                plottype = st.radio("Nu Plot Type", ["Line", "Scatter"])
                logger.debug(f"Nu Evolution plot type set to {plottype}")
            with set00:
                datatype = st.radio("Plot Data", ["Natural", "Gaussian"])
                logger.debug(f"Nu Evolution plot shows {plottype} ZN2 data")
            st.markdown("Interval Modes:")
            set1, set2 = st.columns(2)
            with set1:
                xteintmode = st.selectbox("XTE",
                                          options=["all", "full", "gti"]
                                          )
                logger.debug(f"XTE Interval Mode set to {xteintmode}")
            with set2:
                nicerintmode = st.selectbox("NICER",
                                            options=["all", "full", "gti"]
                                            )
                logger.debug(f"NICER Interval Mode set to {nicerintmode}")
            nitable = querymission(table_in, "NICER")
            nitable = filter_intmode(nitable, nicerintmode)
            xtetable = querymission(table_in, "XTE")
            xtetable = filter_intmode(xtetable, xteintmode)
            table_in = pd.merge(xtetable, nitable, how="outer")
            st.markdown(r"$\nu$ Error (hz)")
            set3, set4 = st.columns(2)
            with set3:
                minerr = st.text_input("Min Error", value=0.0)
                logger.debug(f"Min Error set to {minerr}")
            with set4:
                maxerr = st.text_input("Max Error", value=1.0)
                logger.debug(f"Max Eror set to {maxerr}")
            table_in = filtermax(table_in, "NU_ERR", maxerr)
            table_in = filtermin(table_in, "NU_ERR", minerr)
            minzn2 = st.text_input("Min $Z_n^2$", value=30.0)
            table_in = filtermin(table_in, "ZN2", minzn2)
            st.markdown(r"$\nu$ Gaussian Fit Error (hz)")
            set5, set6 = st.columns(2)
            with set5:
                mingerr = st.text_input("Min Gaussian Error", value=-0.1)
                logger.debug(f"Min Gaussian Fit error set to {mingerr}")
            with set6:
                maxgerr = st.text_input("Max Gaussian Error", value=1.0)
                logger.debug(f"Max Gaussian error set to {maxgerr}")
            table_in = filtermax(table_in, "G_NU_ERR", maxgerr)
            table_in = filtermin(table_in, "G_NU_ERR", mingerr)
            st.markdown(r"Exposure (s)")
            set7, set8 = st.columns(2)
            with set7:
                minexpo = st.text_input("Min Exposure", value=0.0)
                logger.debug(f"Min exposure set to {minexpo}")
            with set8:
                maxexpo = st.text_input("Max Exposure", value=1000000)
                logger.debug(f"Max Exposure set to {maxexpo}")
            table_in = filtermax(table_in, "EXPOSURE", maxexpo)
            table_in = filtermin(table_in, "EXPOSURE", minexpo)
            st.markdown(r"Arrival Times (counts)")
            set9, set10 = st.columns(2)
            with set9:
                minats = st.text_input("Min ATs", value=0.0)
                logger.debug(f"Min Arrival times set to {minats}")
            with set10:
                maxats = st.text_input("Max ATs", value=1000000000)
                logger.debug(f"Max Arrival times set to {maxats}")
            table_in = filtermax(table_in, "ATS", maxats)
            table_in = filtermin(table_in, "ATS", minats)
    st.session_state.show_df = show_df
    st.markdown(r"## $\nu$ Evolution")

    st.plotly_chart(plotting.evo_plot(table_in, plottype, datatype))

    total, tnicer, txte = st.columns(3)
    total.metric("Total", len(table_in["NU"]))
    tnicer.metric("NICER", len(table_in.loc[table_in["MISSION"] == "NICER"]))
    txte.metric("XTE", len(table_in.loc[table_in["MISSION"] == "XTE"]))
    if st.session_state.show_df == "On":
        st.dataframe(table_in, use_container_width=True)
        st.markdown(messages.tableinfo())
    with st.sidebar:
        with st.expander(r"$\delta \nu$ Fitting from $\nu$"):
            set1, set2 = st.columns(2)
            with set1:
                show_nufit = st.radio("Show Plot",
                                      ["Off", "On"])
                logger.debug(f"Delta nu plot {show_nufit}")
                addcrabtime = st.radio("Include CRABTIME",
                                       ["Off", "On"])
                logger.debug(f"CRABTIME {addcrabtime} for dnu plot")
            with set2:
                rmode = st.radio("Residual Mode",
                                 ["Difference", "Coefficient"])
                logger.debug(f"Residual mode set to {rmode}")
                addxray = st.radio("Include X-Ray",
                                   ["Off", "On"],
                                   index=1)
                logger.debug(f"X-Ray {addxray} for dnu plot")
            st.session_state.show_nufit = show_nufit
            trail = st.slider("# of Trailing Boxcar Points", 2, 10, 2)
            logger.debug(f"Boxcar set to {trail} trailing pts")
            lead = st.slider("# of Leading Boxcar Points", 2, 10, 2)
            logger.debug(f"Boxcar set to {lead} leading pts")
            order = st.number_input("Order of Polynomial", min_value=1,
                                    max_value=3, step=1)
            logger.debug(f"Boxcar set to {order} polynomial")

    if "show_nufit" not in st.session_state:
        st.session_state.show_nufit = "Off"
    if st.session_state.show_nufit == "On":
        st.divider()
        st.markdown(r"## $\delta \nu$ Fitting of $\nu$ Data")
        dnuplot = go.Figure()
        if addcrabtime == "On":
            jbcrab = getcrabtime()
            radiodnu = boxcarfit(jbcrab, tpts=trail,
                                 lpts=lead, order=order,
                                 mode=rmode)
            dnuplot.add_trace(go.Scatter(x=radiodnu["TIME"],
                                         y=radiodnu["DNU"],
                                         name='Radio'))
            logger.debug("CRABTIME radio included dnu plot")
        if addxray == "On":
            nuresiduals = boxcarfit(table_in, tpts=trail,
                                    lpts=lead, order=order,
                                    mode=rmode)
            dnuplot.add_trace(go.Scatter(x=nuresiduals["TIME"],
                                         y=nuresiduals["DNU"],
                                         name='X-Ray'))
            logger.debug("X-Ray results included in dnu plot")
        dnuplot.update_layout(xaxis_title="Time",
                              yaxis_title=r"Nu Residual")
        st.plotly_chart(dnuplot, order=order)
        if addcrabtime == "On":
            st.markdown(messages.crabtime_credit())

    try:
        nif = pd.read_sql("SELECT NICER.OBSID, NICER.TWR_FILE FROM NICER " +
                          "WHERE NICER.TWR_FILE IS NOT NULL", con)
    except pd.errors.DatabaseError:
        logger.warning("No NICER Table Found")
        nif = pd.DataFrame({"OBSID": []})

    try:
        xtef = pd.read_sql("SELECT XTE.OBSID, XTE.TWR_FILE FROM XTE " +
                           "WHERE XTE.TWR_FILE IS NOT NULL", con)
    except pd.errors.DatabaseError:
        logger.warning("No XTE Table Found")
        xtef = pd.DataFrame({"OBSID": []})

    resdf = pd.merge(nif, xtef, how="outer")
    resdf = pd.merge(resdf, table_in, how="inner", on="OBSID")
    resdf = resdf.drop_duplicates(subset=["OBSID"])
    st.divider()

    st.markdown("# Individual Measurements")
    obsid = st.selectbox("OBSID", resdf["OBSID"])
    obsid_filt = resdf.loc[resdf["OBSID"] == obsid]
    obsid_filt.reset_index(drop=True, inplace=True)
    rplots = plotting.obsid_plots(obsid_filt["TWR_FILE"][0])
    for num, fig in enumerate(rplots["ZN2"]):
        if fig is not False:
            if num <= 0:
                st.markdown("### Full")
            else:
                st.markdown(f"### Interval {num}")
            sluicing, phasecurve = st.columns(2)
            with sluicing:
                st.pyplot(rplots["ZN2"][num])
            with phasecurve:
                st.pyplot(rplots["Phase"][num])
    with st.expander("More Info On Individual Results"):
        st.markdown(messages.iresult_info())

    st.divider()
    st.markdown(messages.poweredby())


if __name__ == "__main__":
    pargs = parseargs()
    main(pargs)
