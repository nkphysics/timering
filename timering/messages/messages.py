def iresult_info():
    """
    Information about individual result measurements
    """
    info = ("""
            Individual Results are processed in two different modes, FULL and GTI.

            ### GTI Mode (Interval results)

            GTI mode may be misleading from the name.
            Results processed in GTI mode are measurements optained from a singular exposure
            interval during an observation. This is much more applicable to missions like
            NICER that group multiple exposure intervals into a single .evt file.
            To explain as simply as possible, missions in orbit do not get continueous
            exposure from a source due to the missions' orbit, as well as the orbits of other
            objects that may impede the view of a source by the mission. This is dealt with
            by not collecting, or clearing out data when the source is not in view. This
            leads to periods across an observation of for example 1200s of exposure followed
            by >2000s of nothing before another exposure. In some cases this short duration
            between exposures can be lead to noticable differences in measurements of spin
            Therefore each exposure interval is isolated and measured on its own.
            In the above example you can see there are three exposure intervals,
            leading to 3 seperate measurement attempts, one for each interval.


            ### Full Mode

            Full mode measurements are obtained by taking all photon arrival times from
            an observation and using them to obtain a single measurement of spin. This is
            done disregarding total observation exposure or the number of exposure intervals.
            For sources that are easy to obtain measurements from, this can lead to noisy
            measurement results. For sources in which it is harder to obtain timing
            measurements, the full mode can be useful in obtaining something as opposed to
            nothing.

            **DO NOTE:**

            .evt files with only one exposure interval will only be processed
            in full mode. This is very much the case with XTE results since XTE .evt files
            generally come in the form of one exposure interval per .evt file.
            """)
    return info


def tableinfo():
    """
    Information for dataframe dropdown about rounding
    """
    info = (r"""
            **NOTE:**

            Data presneted in this table rounds to $10^{-4}$ Hz by default.
            However, in most cases measurements of $\nu$ are obtained
            at resolutions of $<10^{-5}$ Hz. Selecting single cells of the
            dataframe will display the full values.
            """)
    return info


def crabtime_credit():
    """
    Citation for CRABTIME and JB Radio for the Crab
    """
    info = ("""
            **Credit**

            Radio timing data for the Crab Pulsar (PSR B0531+21) comes from 
            the *Crab Pulsar Monthly Ephemeris* covering a timespan of Feb 15, 1982 
            to the present which was created Dr. Andrew Lyne and multiple 
            collaborators at Jodrell Bank Observatory in the UK, and was accessed
            via the HEASARC's CRABTIME database where the HEASARC created two new 
            parameters which were not present in the original Jodrell Bank tables, 
            the pulsar period and its first derivative.

            **References**

            Lyne, A.G., Jordan, C.A., Roberts, M.E., 
            "Jodrell Bank Crab Pulsar Timing Results, 
            Monthly Ephemeris", http://www.jb.man.ac.uk/~pulsar/crab.html, 
            University of Manchester, Jodrell Bank Observatory, Macclesfield, 
            Cheshire, SK11 9DL, UK.
            """)
    return info


def poweredby():
    """
    Message acknowledging all of the packages uses in development
    """
    info = ("""
            ### Site

            This web dashboard is powered by the following packages. 

            If you find this useful in any way please make sure to give
            credit to all the listed software packages (In no particular order).
            (Software listed without DOIs just means I could not find one. Reach out if you know of one please!)

            - Astropy 
            [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8136839.svg)](https://doi.org/10.5281/zenodo.8136839)

            - Astroquery
            [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5804082.svg)](https://doi.org/10.5281/zenodo.5804082)

            - Matplotlib
            [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8118151.svg)](https://doi.org/10.5281/zenodo.8118151)

            - Numpy
            [![Nature Paper](https://img.shields.io/badge/DOI-10.1038%2Fs41586--020--2649--2-blue)](https://doi.org/10.1038/s41586-020-2649-2)

            - Pandas
            [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3509134.svg)](https://doi.org/10.5281/zenodo.3509134)

            - Plotly
            - Sqlite3
            - Streamlit

            """)
    return info


def heasarc_credit():
    """
    Information to credit HEASARC and missions' archives
    (i.e. NICER and RXTE as of 2023-10-15)
    """
    info = ("""
            ### Research

            This research has made use of data and/or software provided by the 
            High Energy Astrophysics Science Archive Research Center (HEASARC), 
            which is a service of the Astrophysics Science Division at NASA/GSFC.
            """)   
    return info


def sidebar_footer():
    info = ("""
            If you find any of the data and/or dashboard features 
            useful in your projects or work, please consider acknowledging
            *Nicholas Kuechel* and reaching out.
            """)
    return info
