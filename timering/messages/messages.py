def iresult_info():
    """
    Information about individual result measurements
    """
    info = ("""
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

            **DO NOTE:**

            Observations or .evt files with only one exposure interval will only be processed
            in full mode. This is very much so the case with XTE results since XTE .evt files
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
