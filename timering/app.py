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
import json

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
    p.add_argument("--local",
                   default=False,
                   action="store_true",
                   help="Sets the path mode to local rather than web")
    p.add_argument("--config",
                   "-cfg",
                   default=None,
                   type=str,
                   help="Config file for multiple databases")

    pargs = p.parse_args()
    return pargs


def querymission(table_in, mission):
    mtable = table_in.loc[table_in["MISSION"] == mission.upper()]
    return mtable


def filter_intmode(table, intmode):
    """
    Filters a datatable by interval mode
    """
    if intmode != "all":
        modetable = table.loc[table["MODE"] == intmode.lower()]
        if intmode == "gti":
            fulltable = table.loc[table["MODE"] == "full"]
            merged = fulltable[~fulltable.OBSID.isin(modetable.OBSID)]
            modetable = pd.concat([modetable, merged])
    else:
        modetable = table
    return modetable


def filtermax(table, column, maxbound):
    table = table.loc[table[column] <= float(maxbound)]
    return table


def filtermin(table, column, minbound):
    table = table.loc[table[column] >= float(minbound)]
    return table


def filter_obsidout(table, obsid):
    """
    Filters out the input table of a specific OBSID
    """
    table = table.loc[table["OBSID"] != obsid]
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


def parse_config(configpath):
    """
    Parses a config file for multiple databases
    """
    with open(configpath) as config:
        cfg = json.load(config)
        return cfg


def makecsv(df):
    """
    Converts a pandas dataframe to csv for download
    """
    return df.to_csv().encode('utf-8')


class Dashboard:

    def __init__(self, pargs, logger):
        self.logger = logger
        self.sources = []
        self.srcsel = ""
        self.dbpath = self.catalog_options_parsing()
        self.con = self.connect_twdb()
        self.active_src = self.get_src()
        self.obsids = self.get_obsids()
        self.nuresults = self.get_nuresults()

    def catalog_options_parsing(self):
        """
        Parses setup for timering from config or input database

        Returns: Path to selected source database
        """
        with st.sidebar:
            st.markdown("# Catalog Options")
            if pargs.config is None:
                dbpath_in = st.text_input("Local database path",
                                          value=pargs.database)
                dbpath = pathlib.Path(dbpath_in).resolve()
            else:
                self.sources = parse_config(pargs.config)
                self.srcsel = st.selectbox("Source", self.sources)
                dbpath = pathlib.Path(self.sources[self.srcsel]
                                                  ["database"]).resolve()
            return dbpath

    def connect_twdb(self):
        try:
            con = sqlite3.connect(self.dbpath)
            return con
        except sqlite3.OperationalError:
            self.logger.critical("Unable to connect to database")

    def get_src(self):
        src_query = pd.read_sql_query(("""SELECT Source FROM df_metadata
                                          WHERE rowid = 1"""),
                                      self.con)
        src = src_query['Source'][0]
        return src

    def get_obsids(self):
        obsids = pd.read_sql_query(("""SELECT DISTINCT OBSID
                                       from nu_results"""), self.con)
        return obsids["OBSID"]

    def get_nuresults(self):
        table_in = pd.read_sql_query("SELECT * FROM nu_results",
                                     self.con,
                                     parse_dates=["TIME"])
        table_in = table_in.sort_values(by="TIME")
        return table_in

    def min_max_filters(self, column: str, minval: float,
                        maxval: float):
        self.nuresults = filtermax(self.nuresults, column, maxval)
        self.nuresults = filtermin(self.nuresults, column, minval)
        return self.nuresults


def main(pargs: argparse.Namespace):
    level = logging.WARNING
    if pargs.debug is True:
        level = logging.DEBUG
    logging.basicConfig(stream=sys.stdout,
                        level=level)
    logger = logging.getLogger(TIMERING)

    if "show_df" not in st.session_state:
        st.session_state.show_df = False
    dashboard = Dashboard(pargs, logger)
    with st.sidebar:
        show_df = st.toggle("Show Data Table")
    st.markdown(f"# {dashboard.active_src}")

    try:
        alias = dashboard.sources[dashboard.srcsel]["alias"]
        st.markdown(f"**Aliases:** {alias}")
    except KeyError:
        logger.debug("No aliases for {dashboard.srcsel}")
    table_in = dashboard.nuresults
    unfiltereddf = table_in.copy()
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
            obsouts = st.multiselect("OBSID Exclusions", dashboard.obsids)
            st.session_state.obsouts = obsouts
            for i in obsouts:
                table_in = filter_obsidout(table_in, i)
            st.markdown(r"$\nu$ Error (hz)")
            set3, set4 = st.columns(2)
            with set3:
                minerr = st.text_input("Min Error", value=0.0)
                logger.debug(f"Min Error set to {minerr}")
            with set4:
                maxerr = st.text_input("Max Error", value=1.0)
                logger.debug(f"Max Eror set to {maxerr}")
            table_in = dashboard.min_max_filters("NU_ERR", minerr, maxerr)
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
            table_in = dashboard.min_max_filters("G_NU_ERR", mingerr, maxgerr)
            st.markdown(r"Exposure (s)")
            set7, set8 = st.columns(2)
            with set7:
                minexpo = st.text_input("Min Exposure", value=0.0)
                logger.debug(f"Min exposure set to {minexpo}")
            with set8:
                maxexpo = st.text_input("Max Exposure", value=1000000)
                logger.debug(f"Max Exposure set to {maxexpo}")
            table_in = dashboard.min_max_filters("EXPOSURE", minexpo, maxexpo)
            st.markdown(r"Arrival Times (counts)")
            set9, set10 = st.columns(2)
            with set9:
                minats = st.text_input("Min ATs", value=0.0)
                logger.debug(f"Min Arrival times set to {minats}")
            with set10:
                maxats = st.text_input("Max ATs", value=1000000000)
                logger.debug(f"Max Arrival times set to {maxats}")
            table_in = dashboard.min_max_filters("ATS", minats, maxats)
    table_in = table_in.sort_values(by="TIME")
    st.session_state.show_df = show_df
    st.markdown(r"## $\nu$ Evolution")

    st.plotly_chart(plotting.evo_plot(table_in, plottype, datatype))

    total, tnicer, txte = st.columns(3)
    total.metric("Total", len(table_in["NU"]))
    tnicer.metric("NICER", len(table_in.loc[table_in["MISSION"] == "NICER"]))
    txte.metric("XTE", len(table_in.loc[table_in["MISSION"] == "XTE"]))

    dall, dfilt = st.columns(2)
    with dall:
        unfilteredcsv = makecsv(unfiltereddf)
        st.download_button(label="Download Unfiltered Data",
                           data=unfilteredcsv,
                           file_name=f"{dashboard.active_src}-all.csv",
                           mime='text/csv')
    with dfilt:
        unfilteredcsv = makecsv(table_in)
        st.download_button(label="Download Filtered Data",
                           data=unfilteredcsv,
                           file_name=f"{dashboard.active_src}-filtered.csv",
                           mime='text/csv')
    if st.session_state.show_df is True:
        st.dataframe(table_in, use_container_width=True)
        st.markdown(messages.tableinfo())
    with st.sidebar:
        with st.expander(r"$\delta \nu$ Fitting from $\nu$"):
            set1, set2 = st.columns(2)
            with set1:
                show_nufit = st.toggle("Show Plot")
                logger.debug(f"Delta nu plot {show_nufit}")
                addcrabtime = st.toggle("Include CRABTIME")
                logger.debug(f"CRABTIME {addcrabtime} for dnu plot")
                addxray = st.toggle("Include X-Ray", value=True)
                logger.debug(f"X-Ray {addxray} for dnu plot")
            with set2:
                rmode = st.radio("Residual Mode",
                                 ["Difference", "Coefficient"])
                logger.debug(f"Residual mode set to {rmode}")
            st.session_state.show_nufit = show_nufit
            trail = st.slider("# of Trailing Boxcar Points", 2, 10, 2)
            logger.debug(f"Boxcar set to {trail} trailing pts")
            lead = st.slider("# of Leading Boxcar Points", 2, 10, 2)
            logger.debug(f"Boxcar set to {lead} leading pts")
            order = st.number_input("Order of Polynomial", min_value=1,
                                    max_value=3, step=1)
            logger.debug(f"Boxcar set to {order} polynomial")

    if "show_nufit" not in st.session_state:
        st.session_state.show_nufit = False
    if st.session_state.show_nufit:
        st.divider()
        st.markdown(r"## $\delta \nu$ Fitting of $\nu$ Data")
        dnuplot = go.Figure()
        if addcrabtime:
            jbcrab = getcrabtime()
            radiodnu = boxcarfit(jbcrab, tpts=trail,
                                 lpts=lead, order=order,
                                 mode=rmode)
            dnuplot.add_trace(go.Scatter(x=radiodnu["TIME"],
                                         y=radiodnu["DNU"],
                                         name='Radio'))
            logger.debug("CRABTIME radio included dnu plot")
        if addxray:
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
        if addcrabtime:
            st.markdown(messages.crabtime_credit())

    try:
        nif = pd.read_sql("SELECT NICER.OBSID, NICER.TWR_FILE FROM NICER " +
                          "WHERE NICER.TWR_FILE IS NOT NULL", dashboard.con)
    except pd.errors.DatabaseError:
        logger.warning("No NICER Table Found")
        nif = pd.DataFrame({"OBSID": []})

    try:
        xtef = pd.read_sql("SELECT XTE.OBSID, XTE.TWR_FILE FROM XTE " +
                           "WHERE XTE.TWR_FILE IS NOT NULL", dashboard.con)
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
    if pargs.local:
        rplots = plotting.obsid_plots(obsid_filt["TWR_FILE"][0])
    else:
        rplots = plotting.obsid_plots(f'data/{obsid_filt["TWR_FILE"][0]}')
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
    st.markdown("## Acknowledgments")
    st.markdown(messages.heasarc_credit())
    st.markdown(messages.poweredby())


if __name__ == "__main__":
    pargs = parseargs()
    main(pargs)
