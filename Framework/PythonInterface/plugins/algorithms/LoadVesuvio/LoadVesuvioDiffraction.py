from mantid.kernel import config
import mantid.simpleapi as ms

def _exec_single_foil_state_mode(runs, spectra, SUMMED_WS, _LOGGING_):
    """
    Execution path when a single foil state is requested
    """

    if len(runs) > 1:
        raise RuntimeError("Single foil state mode does not currently support summing "
                           "multiple files")

    isis = config.getFacility("ISIS")
    inst_prefix = isis.instrument("VESUVIO").shortName()

    try:
        run_str = inst_prefix + runs[0]
    except ValueError:
        run_str = runs[0]

    all_spectra = [item for sublist in spectra for item in sublist]
    ms.LoadRaw(Filename=run_str, OutputWorkspace=SUMMED_WS, SpectrumList=all_spectra,
               EnableLogging=_LOGGING_)
    raw_group = mtd[SUMMED_WS]
    nperiods = raw_group.size()
    first_ws = raw_group[0]
    foil_out = WorkspaceFactory.create(first_ws)
    x_values = first_ws.readX(0)

    foil_map = SpectraToFoilPeriodMap(nperiods)
    for ws_index, spectrum_no in enumerate(all_spectra):
        self._set_spectra_type(spectrum_no)
        foil_out_periods, foil_thin_periods, _ = self._get_foil_periods()

        if self._diff_opt == "FoilOut":
            raw_grp_indices = foil_map.get_indices(spectrum_no, foil_out_periods)
        elif self._diff_opt == "FoilIn":
            indices_thin = foil_map.get_indices(spectrum_no, foil_thin_periods)
            indices_thick = foil_map.get_indices(spectrum_no, foil_thin_periods)
            raw_grp_indices = indices_thin + indices_thick
        elif self._diff_opt == "FoilInOut":
            raw_grp_indices = range(0, self._nperiods)
        else:
            raise RuntimeError("Unknown single foil mode: %s." % (self._diff_opt))

        dataY = foil_out.dataY(ws_index)
        dataE = foil_out.dataE(ws_index)
        for group_index in raw_grp_indices:
            dataY += raw_group[group_index].readY(ws_index)
            dataE += np.square(raw_group[group_index].readE(ws_index))
        np.sqrt(dataE, dataE)
        foil_out.setX(ws_index, x_values)

    ip_file = self.getPropertyValue(INST_PAR_PROP)
    if len(ip_file) > 0:
        self.foil_out = self._load_ip_file(self.foil_out, ip_file)

    if self._sumspectra:
        self._sum_all_spectra()

    ms.DeleteWorkspace(Workspace=SUMMED_WS)
    self._store_results()

#----------------------------------------------------------------------------------------