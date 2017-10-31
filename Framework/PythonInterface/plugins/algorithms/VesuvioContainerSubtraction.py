from mantid.api import *
from mantid.kernel import *
from mantid.simpleapi import Divide, MonteCarloAbsorption, VesuvioCalculateMS

from vesuvio.base import VesuvioBase

class VesuvioContainerSubtraction(VesuvioBase):

    def summary(self):
        return "A container subtraction routine for Vesuvio, which accounts for multiple scattering in the container."

    def category(self):
        return 'Inelastic\\Indirect\\Vesuvio'

    def PyInit(self):
        self.declareProperty(MatrixWorkspaceProperty("InputWorkspace", "",
                                                     direction=Direction.Input),
                             doc="The input sample workspace in TOF")

        self.declareProperty(MatrixWorkspaceProperty("ContainerWorkspace", "",
                                                     direction=Direction.Input),
                             doc="The input container workspace in TOF")

        self.declareProperty(MatrixWorkspaceProperty("OutputWorkspace", "",
                                                     direction=Direction.Output),
                             doc="The single scattering sample workspace")

    def _get_properties(self):
        self._sample_workspace = self.getProperty("InputWorkspace").value
        self._container_workspace = self.getProperty("ContainerWorkspace").value

    def PyExec(self):
        self._sample_container_sum = self._sample_workspace + self._container_workspace

        # Calculate absorption factors using a Monte Carlo routine
        container_absorption_factor = MonteCarloAbsorption(InputWorkspace=self._container_workspace,
                                                           OutputWorkspace="container_abs_factor",
                                                           StoreInADS=False)
        sum_absorption_factor = MonteCarloAbsorption(InputWorkspace=self._sample_container_sum,
                                                     OutputWorkspace="sum_abs_factor",
                                                     StoreInADS=False)

        # Divide original container and Sample+Container by absorption factors
        absorption_corrected_container = Divide(LHSWorkspace=self._container_workspace,
                                                RHSWorkspace=container_absorption_factor,
                                                OutputWorkspace="abs_corrected_container",
                                                StoreInADS=False)
        absorption_corrected_sum = Divide(LHSWorkspace=self._sample_workspace,
                                          RHSWorkspace=sum_absorption_factor,
                                          OutputWorkspace="abs_corrected_sum",
                                          StoreInADS=False)

        # Calculate multiple scattering workspaces from
        _, container_multiple_scattering = VesuvioCalculateMS(InputWorkspace=absorption_corrected_container,
                                                              TotalScatteringWS="container_total_scattering",
                                                              MultipleScatteringWS="container_multiple_scattering",
                                                              StoreInADS=False)
        _, sum_multiple_scattering = VesuvioCalculateMS(InputWorkspace=absorption_corrected_sum,
                                                        TotalScatteringWS="sum_total_scattering",
                                                        MultipleScatteringWS="sum_multiple_scattering",
                                                        StoreInADS=False)

        container_single_scattering = self._container_workspace - container_multiple_scattering
        sum_single_scattering = self._sample_container_sum - sum_multiple_scattering
        sample_single_scattering = self._sample_workspace - (container_single_scattering + sum_single_scattering)

        self.setProperty("OutputWorkspace", sample_single_scattering)