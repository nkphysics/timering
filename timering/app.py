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
    intmode = st.selectbox("Interval mode selection:",
                           options=["all", "full", "gti"]
                           )
    show_df = st.radio(
        "Show Data Table",
        ["Off", "On"],
    )

st.session_state.show_df = show_df

con = sqlite3.connect(dbpath)
table_in = pd.read_sql_query("SELECT * FROM nu_results", con)
table_in["TIME"] = pd.to_datetime(table_in["TIME"])
table_in = table_in.sort_values(by="TIME")
if intmode != "all":
    table_in = table_in.loc[table_in["MODE"] == intmode]


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
    fig.update_xaxes(showgrid=False)
    fig.update_yaxes(showgrid=False)
    return fig


st.plotly_chart(tevo_plot(plottype))
if st.session_state.show_df == "On":
    crab_table = st.dataframe(table_in)
