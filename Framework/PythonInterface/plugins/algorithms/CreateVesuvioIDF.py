
import mantid.simpleapi as api
from mantid.kernel import *
from mantid.api import *

import os

IP_PROP = "IPFilename"
IDF_PROP = "OutputFilename"

class CreateVesuvioIDF(PythonAlgorithm):

    _phi_data = [214.730191, 208.517121, 202.010235, 195.421267, 188.96977, 182.844226, 211.242111, 204.716592, 197.985822, 191.286993, 184.846089, 213.93161, 207.466845, 200.65919, 193.756495, 187.017631, 210.210156,
    203.405588, 196.360378, 189.356977, 206.177646, 199.068035, 191.851598, 184.82908, 208.915322, 201.835797, 194.480414, 187.177719, 204.608279, 197.205619, 189.685619, 199.975591, 192.323117, 184.737973, 202.726124,
    195.048535, 187.244444, 197.801361, 189.877995, 192.587675, 184.533139, 187.150599, 189.834159, 184.176128, 221.938415, 223.306516, 223.92399, 223.686699, 222.63518, 220.938585, 220.959129, 221.965524, 222.159315,
    221.504586, 220.11456, 218.367157, 219.706913, 220.315223, 220.08216, 219.04887, 217.177614, 218.156711, 218.344202, 217.706834, 215.701391, 216.286262, 216.061609, 215.069539, 212.987781, 213.915873, 214.095722,
    213.490117, 211.261975, 211.809309, 211.597774, 209.221666, 209.386564, 208.831155, 206.376447, 206.866572, 206.67748, 204.071389, 204.214864, 201.467663, 201.308245, 198.603367, 195.656907, 192.61408, 185.205584,
    191.595814, 198.222015, 204.873935, 211.325168, 217.375886, 187.421932, 194.101789, 200.929073, 207.65962, 214.061296, 183.145456, 189.806538, 196.745008, 203.71478, 210.44977, 185.375212, 192.350125, 199.497462,
    206.535367, 187.780971, 195.030842, 202.317551, 209.335787, 183.090513, 190.348958, 197.813852, 205.153973, 185.505069, 193.053415, 200.650003, 188.084353, 195.849832, 203.477304, 182.973707, 190.795307, 198.682133,
    185.551916, 193.588352, 188.257475, 196.399504, 191.033554, 185.461205, 188.213211, 34.554567, 26.679833, 18.558404, 10.507388, 212.725413, 204.626927, 196.398361, 188.370964, 30.461603, 24.581168, 18.615915, 12.691706,
    209.016303, 203.06463, 197.079934, 191.19273, 34.260946, 26.398807, 18.301171, 10.282116, 212.755407, 204.689201, 196.509329, 188.538513, 30.206595, 24.297026, 18.31852, 12.399372, 209.956038, 204.020647, 198.023792,
    192.09454, 30.648007, 24.733247, 18.745522, 12.814072, 209.999575, 204.039724, 198.015424, 192.059775, 34.155443, 26.180059, 17.998443, 9.934835, 213.645682, 205.570229, 197.316915, 189.21912, 30.075572, 24.190557,
    18.223872, 12.304715, 209.828899, 203.877328, 197.86838, 191.93202, 33.381438, 25.419044, 17.266051, 9.24428, 213.672724, 205.694631, 197.552196, 189.565927]


    def summary(self):
        return "Uses an IP file to create an IDF for the Vesuvio instrument at ISIS."

    def version(self):
        return 1

    def category(self):
        return "Inelastic;PythonAlgorithms"

    def PyInit(self):
        self.declareProperty(FileProperty(IP_PROP, "", action=FileAction.OptionalLoad,
                                          extensions=["par"]),
                              doc="The file path to the Vesuvio IP file.")

        self.declareProperty(FileProperty(IDF_PROP, "", Direction.Output,
                                          extensions=["xml"]),
                              doc="The file path to the desired IDF file.")

    def validateInputs(self):
        # Requires check to ensure files are valid
        # .xml on idf
        issues = dict()
        ip_filename = self.getPropertyValue(IP_PROP)
        idf_filename = self.getPropertyValue(IDF_PROP)

        return issues

    def PyExec(self):
        ip_filename = self.getPropertyValue(IP_PROP)
        idf_filename = self.getPropertyValue(IDF_PROP)
        ip_lines = self._read_ip(ip_filename)

        # Read headings
        ip_heading = ip_lines[0].split()

        # Read Monitors
        ip_monitor_one = ip_lines[1].split()
        ip_monitor_two = ip_lines[2].split()

        # Construct full IDF
        idf_string = self._construct_idf_header()
        idf_string += self._construct_idf_source_sample(ip_lines[3].split()[4])
        idf_string += self._construct_idf_monitors(ip_monitor_one, ip_monitor_two)
        idf_string += self._construct_idf_foils()
        idf_string += self._construct_idf_forward(ip_lines)
        idf_string += self._construct_idf_back(ip_lines)
        idf_string += self._construct_idf_pixel()
        idf_string += self._construct_idf_det_list()

        self._write_idf(idf_filename, idf_string)


#===========================================================================================================================================#
    def _read_ip(self, ip_file_name):
        # Read IP file
        ip_file = open(ip_file_name, 'r')
        ip_lines = ip_file.readlines()
        return ip_lines


#===========================================================================================================================================#
    def _construct_idf_header(self):
        # Define the header of the IDF
        idf_header = (
        '<?xml version="1.0" encoding="UTF-8"?>\n' +
        '<!-- For help on the notation used to specify an Instrument Definition File\n' +
        '     see http://www.mantidproject.org/IDF -->\n' +
        '<instrument xmlns="http://www.mantidproject.org/IDF/1.0"\n' +
        '            xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"\n' +
        '            xsi:schemaLocation="http://www.mantidproject.org/IDF/1.0 http://schema.mantidproject.org/IDF/1.0/IDFSchema.xsd"\n' +
        '            name="VESUVIO"\n' +
        '            valid-from="2006-08-11 00:00:00"\n' +
        '            valid-to="2100-01-31 23:59:59"\n' +
        '            last-modified="2011-07-20 00:00:21">\n' +
        '\n' +
        '<defaults>\n' +
        '   <length unit="meter"/>\n' +
        '   <angle unit="degree"/>\n' +
        '   <reference-frame>\n' +
        '      <!-- The z-axis is set parallel to and in the direction of the beam. the\n' +
        '           y-axis points up and the coordinate system is right handed. -->\n' +
        '      <along-beam axis="z"/>\n' +
        '      <pointing-up axis="y"/>\n' +
        '      <handedness val="right"/>\n' +
        '   </reference-frame>\n' +
        '</defaults>' +
        '\n')

        return idf_header


#===========================================================================================================================================#
    def _construct_idf_source_sample(self, non_monitor_L0):
        # Define the source and sample positions
        idf_src_sam = (
        '  <!--  SOURCE AND SAMPLE POSITION -->\n'+
        '\n' +
        '  <component type="moderator">\n' +
        ('    <location z="-%s" />\n' % non_monitor_L0) +
        '  </component>\n' +
        '\n' +
        '  <type name="moderator" is="Source" />\n' +
        '\n' +
        '  <component type="sample-position">\n' +
        '    <location />\n' +
        '  </component>\n' +
        '\n' +
        '  <type name="sample-position" is="SamplePos">\n' +
        '  </type>\n' +
        '\n')

        return idf_src_sam


#===========================================================================================================================================#
    def _construct_idf_monitors(self, monitor_one, monitor_two):
        # Define the Monitors
        idf_monitors = (
        '  <!-- MONITORS -->\n' +
        '\n' +
        '  <component type="monitor" idlist="monitors">\n' +
        ('    <location z="%s" name="Monitor 1"/>\n' % monitor_one[5]) +
        ('    <location z="%s" name="Monitor 2"/>\n' % monitor_two[5]) +
        '  </component>\n' +
        '\n' +
        '  <type name="monitor" is="monitor">\n' +
        '   <cuboid id="shape">\n' +
        '      <left-front-bottom-point x="0.01" y="-0.01" z="0.0"  />\n' +
        '      <left-front-top-point  x="0.01" y="-0.01" z="0.02"  />\n' +
        '      <left-back-bottom-point  x="-0.01" y="-0.01" z="0.0"  />\n' +
        '      <right-front-bottom-point  x="0.01" y="0.01" z="0.0"  />\n' +
        '   </cuboid>\n' +
        '  </type>\n' +
        '\n')

        return idf_monitors


#===========================================================================================================================================#
    def _construct_idf_foils(self):
        # Define the Foils and Foil changers
        idf_foils = (
        '  <!-- FOILS  -->\n' +
        '\n' +
        '  <!-- midpoint positions\n' +
        '       phi is simply set to put them all in the scattering plane\n' +
        '  -->\n' +
        '  <component name="foil-pos0" type="foil">\n' +
        '   <!-- position 0 -->\n' +
        '   <location t="-42" r="0.225" p="0"/>\n' +
        '   <location t="-64" r="0.225" p="0"/>\n' +
        '   <location t="35" r="0.225" p="0"/>\n' +
        '   <location t="55" r="0.225" p="0"/>\n' +
        '  </component>\n' +
        '\n' +
        '  <component name="foil-pos1" type="foil">\n' +
        '   <!-- position 1 -->\n' +
        '   <location t="-31" r="0.225" p="0"/>\n' +
        '   <location t="-50" r="0.225" p="0"/>\n' +
        '   <location t="46" r="0.225" p="0"/>\n' +
        '   <location t="67" r="0.225" p="0"/>\n' +
        '  </component>\n' +
        '\n' +
        '  <!-- The foil shape is currently set to give the correct dTheta as seen\n' +
        '       from the sample position -->\n' +
        '  <type name="foil">\n' +
        '    <cuboid id="foil-shape-shape">\n' +
        '      <left-front-bottom-point x="-0.02000121679" y="0.001" z="-0.001"/>\n' +
        '      <left-front-top-point x="-0.02000121679" y="0.001" z="0.001"/>\n' +
        '      <left-back-bottom-point x="-0.02000121679" y="-0.001" z="-0.001" />\n' +
        '      <right-front-bottom-point x="0.02000121679" y="0.001" z="-0.001" />\n' +
        '    </cuboid>\n' +
        '  </type>\n' +
        '\n' +
        '  <!-- FOIL CHANGER -->\n' +
        '  <!-- The shape & location are made up to provide inputs for the\n' +
        '       CalculateGammaBackground algorithm -->\n' +
        '  <component type="foil-changer">\n' +
        '  <location x="0" y ="0" z="0"/>\n' +
        '  </component>\n' +
        '\n' +
        '  <!-- Important quantity currently is the height. It defines the Z range\n' +
        '       for the integration of the foil positions  (-0.2,0.2)-->\n' +
        '  <type name="foil-changer">\n' +
        '    <cylinder id="foil-changer-shape">\n' +
        '    <centre-of-bottom-base x="0.0" y="-0.2" z="0.0" />\n' +
        '    <axis x="0.0" y="1" z="0" />\n' +
        '    <radius val="0.1" />\n' +
        '    <height val="0.4" />\n' +
        '    </cylinder>\n' +
        '  </type>\n' +
        '\n')
        return idf_foils


#===========================================================================================================================================#
    def _construct_idf_forward(self, ip_lines):
        # Define the forward scattering detectors
        idf_forward_det_head = (
        '  <!-- DETECTORS -->\n' +
        '\n' +
        '  <component type="forward" idlist="forward">\n' +
        '    <location />\n' +
        '  </component>\n' +
        '\n' +
        '\n' +
        '  <component type="back" idlist="back">\n' +
        '    <location />\n' +
        '  </component>\n' +
        '\n' +
        '  <type name="forward">\n' +
        '    <component type="pixel forward" >\n')

        # Define all locations of forward scattering detectors
        idf_forward_locations = ''
        for idx in range(135,  199):
            # Get corresponding detector from IP file
            ip_det = ip_lines[idx].split()
            idf_forward_locations += ('      <location t="%f" r="%s" p="%f" name="S%s"/>\n' % (float(ip_det[2]), ip_det[5], self._phi_data[idx - 3], ip_det[1]))

        idf_forward_det_end = (
        '    </component>\n' +
        '  </type>\n' +
        '\n')

        # Construct forward detector definitions string
        idf_forward_dets = idf_forward_det_head + idf_forward_locations + idf_forward_det_end
        return idf_forward_dets


#===========================================================================================================================================#
    def _construct_idf_back(self, ip_lines):
        # Define the back scattering detectors
        idf_back_head =(
        '   <type name="back">\n' +
        '     <component type="pixel back" >\n')

        # Define all locations of back scattering detectors
        idf_back_locations = ''
        for idx in range (3,135):
            # Get corresponding detector from IP file
            ip_det = ip_lines[idx].split()
            idf_back_locations += ('      <location t="%f" r="%s" p="%f" name="S%s"/>\n' % (float(ip_det[2]), ip_det[5], self._phi_data[idx - 3], ip_det[1]))

        idf_back_end = (
        '    </component>\n' +
        '  </type>\n')

        # Construct back detector definitions string
        idf_back_dets = idf_back_head + idf_back_locations + idf_back_end
        return idf_back_dets


#===========================================================================================================================================#
    def _construct_idf_pixel(self):
        # Define pixel forward/back
        idf_pix = (
        '  <type name="pixel forward" is="detector">\n' +
        '    <cuboid id="shape">\n' +
        '      <left-front-bottom-point x="0.0125" y="-0.0395" z= "0.0045" />\n' +
        '      <left-front-top-point x="0.0125" y="0.0395" z= "0.0045" />\n' +
        '      <left-back-bottom-point x="0.0125" y="-0.0395" z= "-0.0045" />\n' +
        '      <right-front-bottom-point x="-0.0125" y="-0.0395" z= "0.0045" />\n' +
        '    </cuboid>\n' +
        '    <algebra val="shape" />\n' +
        '  </type>\n' +
        '\n' +
        '  <type name="pixel back" is="detector">\n' +
        '    <cuboid id="shape" >\n' +
        '      <left-front-bottom-point x="0.0125" y="-0.04" z= "0.0045" />\n' +
        '      <left-front-top-point x="0.0125" y="0.04" z= "0.0045" />\n' +
        '      <left-back-bottom-point x="0.0125" y="-0.04" z= "-0.0045" />\n' +
        '      <right-front-bottom-point x="-0.0125" y="-0.04" z= "0.0045" />\n' +
        '    </cuboid>\n' +
        '  <algebra val="shape" />\n' +
        '  </type>\n' +
        '\n')

        return idf_pix


#===========================================================================================================================================#
    def _construct_idf_det_list(self):
        # Detector ID Lists + tail
        idf_det_id_list = (
        '  <!-- DETECTOR ID LISTS -->\n' +
        '\n' +
        '  <idlist idname="back">\n' +
        '    <id start="2101" end="2144" />\n' +
        '    <id start="2201" end="2244" />\n' +
        '    <id start="2301" end="2344" />\n' +
        '  </idlist>\n' +
        '\n' +
        '  <idlist idname="monitors">\n' +
        '    <id start="1101" end="1102" />\n' +
        '  </idlist>\n' +
        '\n' +
        '  <idlist idname="forward">\n' +
        '    <id start="3101" end="3132" />\n' +
        '    <id start="3201" end="3232" />\n' +
        '  </idlist>\n' +
        '\n' +
        '</instrument>\n')

        return idf_det_id_list


#===========================================================================================================================================#
    def _write_idf(self, idf_filename, idf_string):
        # Write IDF out to .xml file
        idf = open(idf_filename, 'w')
        idf.write(idf_string)
        idf.close()

# Register algorithm with Mantid
AlgorithmFactory.subscribe(CreateVesuvioIDF)
