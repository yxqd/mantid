
class SpectraToFoilPeriodMap(object):
    """Defines the mapping between a spectrum number
    & the period index into a WorkspaceGroup for a foil state.
    2 period :: forward scattering
    3 period :: back scattering
    6 period :: back & forward scattering

    one_to_one          :: Only used in back scattering where there is a single
                           static foil
    odd_even/even_odd   :: Only used in forward scatter models when the foil
                           is/isn't infront of each detector. First bank 135-142
                           is odd_even, second (143-150) is even_odd and so on.
    """

    def __init__(self, nperiods=6):
        """Constructor. For nperiods set up the mappings"""
        if nperiods == 2:
            self._one_to_one = {1:1, 2:2}   # Kept for use in reorder method
            self._odd_even =   {1:1, 2:2}
            self._even_odd =   {1:2, 2:1}
        elif nperiods == 3:
            self._one_to_one = {1:1, 2:2, 3:3}
        elif nperiods == 6:
            self._one_to_one = {1:1, 2:2, 3:3, 4:4, 5:5, 6:6}
            self._odd_even =   {1:1, 2:3, 3:5, 4:2, 5:4, 6:6}
            self._even_odd =   {1:2, 2:4, 3:6, 4:1, 5:3, 6:5}
        else:
            raise RuntimeError("Unsupported number of periods given: " + str(nperiods) +
                               ". Supported number of periods=2,3,6")

#----------------------------------------------------------------------------------------

    def reorder(self, arr):
        """
           Orders the given array by increasing value. At the same time
           it reorders the 1:1 map to match this order
           numpy
        """
        vals = np.array(self._one_to_one.values())
        sorted_indices = arr.argsort()
        vals = vals[sorted_indices]
        arr = arr[sorted_indices]
        self._one_to_one = {}
        for index,val in enumerate(vals):
            self._one_to_one[index+1] = int(val)

        return arr

#----------------------------------------------------------------------------------------

    def get_foilout_periods(self, spectrum_no):
        """Returns a list of the foil-out periods for the given
        spectrum number. Note that these start from 1 not zero
            @param spectrum_no :: A spectrum number (1->nspectra)
            @returns A list of period numbers for foil out state
        """
        return self.get_foil_periods(spectrum_no, state=0)

#----------------------------------------------------------------------------------------

    def get_foilin_periods(self, spectrum_no):
        """Returns a list of the foil-out periods for the given
        spectrum number. Note that these start from 1 not zero
            @param spectrum_no :: A spectrum number (1->nspectra)
            @returns A list of period numbers for foil out state
        """
        return self.get_foil_periods(spectrum_no, state=1)

#----------------------------------------------------------------------------------------

    def get_foil_periods(self, spectrum_no, state):
        """Returns a list of the periods for the given
        spectrum number & foil state. Note that these start from 1 not zero
            @param spectrum_no :: A spectrum number (1->nspectra)
            @param state :: 0 = foil out, 1 = foil in.
            @returns A list of period numbers for foil out state
        """
        self._validate_spectrum_number(spectrum_no)

        foil_out = (state==0)

        if spectrum_no < 135:
            foil_periods = [1,2,3]
        elif (spectrum_no >= 135 and spectrum_no <= 142) or \
             (spectrum_no >= 151 and spectrum_no <= 158) or \
             (spectrum_no >= 167 and spectrum_no <= 174) or \
             (spectrum_no >= 183 and spectrum_no <= 190):
            foil_periods = [2,4,6] if foil_out else [1,3,5]
        else:
            foil_periods = [1,3,5] if foil_out else [2,4,6]
        return foil_periods

#----------------------------------------------------------------------------------------

    def get_indices(self, spectrum_no, foil_state_numbers):
        """
        Returns a tuple of indices that can be used to access the Workspace within
        a WorkspaceGroup that corresponds to the foil state numbers given
        @param spectrum_no :: A spectrum number (1->nspectra)
        @param foil_state_no :: A number between 1 & 6(inclusive) that defines which foil
                                state is required
        @returns A tuple of indices in a WorkspaceGroup that gives the associated Workspace
        """
        indices = []
        for state in foil_state_numbers:
            indices.append(self.get_index(spectrum_no, state))
        return tuple(indices)

#----------------------------------------------------------------------------------------

    def get_index(self, spectrum_no, foil_state_no):
        """Returns an index that can be used to access the Workspace within
        a WorkspaceGroup that corresponds to the foil state given
            @param spectrum_no :: A spectrum number (1->nspectra)
            @param foil_state_no :: A number between 1 & 6(inclusive) that defines which
                                        foil state is required
            @returns The index in a WorkspaceGroup that gives the associated Workspace
        """

        self._validate_foil_number(foil_state_no)
        self._validate_spectrum_number(spectrum_no)

        # For the back scattering banks or foil states > 6 then there is a 1:1 map
        if foil_state_no > 6 or spectrum_no < 135:
            foil_periods = self._one_to_one
        elif (spectrum_no >= 135 and spectrum_no <= 142) or \
             (spectrum_no >= 151 and spectrum_no <= 158) or \
             (spectrum_no >= 167 and spectrum_no <= 174) or \
             (spectrum_no >= 183 and spectrum_no <= 190):
             # For each alternating forward scattering bank :: foil_in = 1,3,5, foil out = 2,4,6
            foil_periods = self._odd_even
        else:
            # foil_in = 2,4,6 foil out = 1,3,5
            foil_periods = self._even_odd

        foil_period_no = foil_periods[foil_state_no]
        return foil_period_no - 1 # Minus 1 to get to WorkspaceGroup index

#----------------------------------------------------------------------------------------

    def _validate_foil_number(self, foil_number):
        if foil_number < 1 or foil_number > 6:
            raise ValueError("Invalid foil state given, expected a number between "
                             "1 and 6. number=%d" % foil_number)

#----------------------------------------------------------------------------------------

    def _validate_spectrum_number(self, spectrum_no):
        if spectrum_no < 1 or spectrum_no > 198:
            raise ValueError("Invalid spectrum given, expected a number between 3 "
                             "and 198. spectrum=%d" % spectrum_no)
