import streamlit as st
import pandas as pd
import argparse
import pathlib
import plotly.express as px
import sqlite3


p = argparse.ArgumentParser(description="Set local database")

p.add_argument("--database",
               "-db",
               default="",
               type=str,
               help="Local path to .twdb"
               )

pargs = p.parse_args()

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

st.session_state.show_df = show_df

con = sqlite3.connect(dbpath)
table_in = pd.read_sql_query("SELECT * FROM nu_results", con)
table_in["TIME"] = pd.to_datetime(table_in["TIME"])
table_in = table_in.sort_values(by="TIME")


def querymission(mission):
    mtable = table_in.loc[table_in["MISSION"] == mission.upper()]
    return mtable


def filter_intmode(table, intmode):
    if intmode != "all":
        table = table.loc[table["MODE"] == intmode.lower()]
    return table


nitable = querymission("NICER")
nitable = filter_intmode(nitable, nicerintmode)
xtetable = querymission("XTE")
xtetable = filter_intmode(xtetable, xteintmode)
table_in = pd.merge(xtetable, nitable, how="outer")


def tevo_plot(plottype):
    fig = 0
    if plottype == "Scatter":
        fig = px.scatter(
            table_in,
            x="TIME",
            y="NU",
            title="Timing Evolution",
            color="MISSION",
            error_y="NU_ERR",
        )
    else:
        fig = px.line(
            table_in,
            x="TIME",
            y="NU",
            title="Timing Evolution",
            color="MISSION",
            error_y="NU_ERR",
            markers=True,
        )
    fig.update_layout(xaxis_title="Time", yaxis_title=r"Spin Frequency $\nu$")
    fig.update_xaxes(showgrid=True)
    fig.update_yaxes(showgrid=True)
    return fig


st.plotly_chart(tevo_plot(plottype))
if st.session_state.show_df == "On":
    crab_table = st.dataframe(table_in)
