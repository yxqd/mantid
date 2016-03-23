
import mantid.simpleapi as api
from mantid.kernel import *
from mantid.api import *

import os

IP_PROP = "IPFilename"
IDF_PROP = "OutputFilename"

class CreateVesuvioIDF(PythonAlgorithm):

    _ip_filename = None
    _idf_filename = None
    _ip_lines = None

    _phi_data = [-128.0874448587, -137.4504101855, -147.5076669307, -157.6896897223, -167.3819019002, -176.1271875294,
                 -131.6778173977, -141.9927254194, -152.8143229161, -163.4092563484, -173.1233208227, -125.5728641099,
                 -135.920364941, -147.2440564461, -158.7466822376, -169.5515011966, -129.3597185468, -140.9477882504,
                 -153.2771829659, -165.2706728903, -133.9553425886, -146.9007277786, -160.1136938725, -172.2910414134,
                 -126.3867483977, -139.5716272365, -153.8908842374, -167.6910709804, -131.351084127, -146.4357144159,
                 -161.9469223465, -137.6759059782, -154.7430250679, -170.9489870959, -127.7502825245, -145.7691164574,
                 -164.532494474, -134.9062407181, -156.0029934465, -144.7334002164, -168.4268599076, -158.0625805251,
                 -142.9148179861, -162.0079403863, -68.1257015311, -77.57486040664, -87.75716862714, -98.08253795803,
                 -107.9103673217, -116.765676209, -71.71883636505, -82.14998513379, -93.12701665152, -103.8819784063,
                 -113.7338283201, -65.53020065802, -75.97697601007, -87.45720162623, -99.14286950758, -110.1186688899,
                 -69.31240990386, -81.04350979932, -93.57207630107, -105.7735702416, -73.92089334184, -87.06496894932,
                 -100.5208546301, -112.9145564586, -66.21877860108, -79.57784832523, -94.16469751451, -108.2413122069,
                 -71.1838519919, -86.53046844038, -102.380546841, -77.5485466938, -94.99358022811, -111.5828177633,
                 -67.3825036465, -85.75896123864, -105.0206934853, -74.56321187153, -96.2334165996, -84.54850079755,
                 -109.0417372572, -98.28877633342, -82.38103301845, -102.3389223272, -7.286972691881, -16.74717300156,
                 -26.9628059271, -37.34310958122, -47.23527744391, -56.15102701443, -10.8610028062, -21.31478934394,
                 -32.34284398044, -43.1663899362, -53.09016806635, -4.647687993777, -15.1019959278, -26.62198853621,
                 -38.37898474685, -49.43630126273, -8.40042018502, -20.15613881487, -32.74895831707, -45.04244475045,
                 -12.98195005424, -26.17627045553, -39.72582542777, -52.23013788624, -5.238226866502, -18.61797021173,
                 -33.2925121702, -47.49621735801, -10.15991004596, -25.56695777916, -41.55059224045, -16.4877440155,
                 -34.05411064343, -50.82612175418, -6.239303832486, -24.68344275415, -44.15037262656, -13.34825884272,
                 -35.19633286048, -23.29566877746, -48.13132146988, -37.10502214547, -20.79453568439, -40.91480165425,
                 38.54093084891, 30.15352921739, 21.2028416006, 12.09091666095, -36.49183354547, -27.82571557711,
                 -18.7223644126, -9.620127243978, 37.39080432975, 30.74759797313, 23.67019659512, 16.34386944053,
                 -36.44072467876, -29.49322449279, -22.16376160952, -14.67460757133, 46.7771602103, 37.86029267043,
                 27.43596591805, 15.92533800454, -44.97943637152, -35.54610080611, -24.74581923439, -13.14843524179,
                 49.78939754585, 42.64273002831, 34.13076467976, 24.30605588768, -49.76072934472, -42.39027877539,
                 -33.65580915066, -23.65726347469, 134.0380602071, 141.1787547125, 149.3227541915, 158.3041801021,
                 -134.2400476383, -141.6213992416, -150.03697541, -159.2857909202, 136.1229316768, 145.1887786041,
                 155.3618978243, 166.1336385912, -136.2873427672, -145.5924013339, -156.0151372065, -166.9948039869,
                 144.5585845472, 151.0823954816, 157.9467472491, 164.969861814, -144.7288972613, -151.371545539,
                 -158.3264383956, -165.4029832258, 143.922577583, 152.2803073279, 161.036377481, 169.8005699682,
                 -143.6543487997, -151.9955223602, -160.7178549141, -169.4348763528]


    def summary(self):
        return "Uses an IP file to create an IDF for the Vesuvio instrument at ISIS."

    def version(self):
        return 1

    def category(self):
        return "Inelastic;PythonAlgorithms"

    def PyInit(self):
        self.declareProperty(FileProperty(IP_PROP, "", action=FileAction.OptionalLoad,
                                          extensions=["par","dat"]),
                              doc="The file path to the Vesuvio IP file.")

        self.declareProperty(FileProperty(IDF_PROP, "", Direction.Output,
                                          extensions=["xml"]),
                              doc="The file path to the desired IDF file.")

    def validateInputs(self):
        # Requires check to ensure files are valid
        # .xml on idf
        self.getProperties()
        issues = dict()

        # Validate output filename
        output_ext = self._idf_filename[self._idf_filename.find('.'):]
        if output_ext != '.xml':
            issues[IDF_PROP] = 'OutputFileName does not have .xml extension. Please supply an OutputFileName ending in .xml'

        if len(self._ip_lines) != 199:
            issues[IP_PROP] = 'IP file was not the correct length. expected length is 199 lines'

        return issues

    def PyExec(self):

        # Read headings
        ip_heading = self._ip_lines[0].split()

        # Read Monitors
        ip_monitor_one = self._ip_lines[1].split()
        ip_monitor_two = self._ip_lines[2].split()

        # Construct full IDF
        idf_string = self._construct_idf_header()
        idf_string += self._construct_idf_source_sample(self._ip_lines[3].split()[4])
        idf_string += self._construct_idf_monitors(ip_monitor_one, ip_monitor_two)
        idf_string += self._construct_idf_foils()
        idf_string += self._construct_idf_forward()
        idf_string += self._construct_idf_back()
        idf_string += self._construct_idf_pixel()
        idf_string += self._construct_idf_det_list()

        self._write_idf(self._idf_filename, idf_string)

    def getProperties(self):
        self._ip_filename = self.getPropertyValue(IP_PROP)
        self._idf_filename = self.getPropertyValue(IDF_PROP)
        self._read_ip(self._ip_filename)


#===========================================================================================================================================#
    def _read_ip(self, ip_file_name):
        # Read IP file
        ip_file = open(ip_file_name, 'r')
        self._ip_lines = ip_file.readlines()
        return self._ip_lines


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
    def _construct_idf_forward(self):
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
        '  <type name="forward">\n')

        # Define all locations of forward scattering detectors
        idf_forward_locations = ''
        for idx in range(135,  199):
            # Get corresponding detector from IP file
            ip_det = self._ip_lines[idx].split()
            idf_forward_locations += ('    <component type="pixel back" >\n')
            idf_forward_locations += ('      <location t="%f" r="%s" p="%f" name="S%s"/>\n' % (float(ip_det[2]), ip_det[5], self._phi_data[idx - 3], ip_det[1]))
            idf_forward_locations += ('      <parameter name="t0"> <value val="%f"/> </parameter>\n' % (float(ip_det[3])))
            idf_forward_locations += ('    </component>\n')

        idf_forward_det_end = (
        '  </type>\n' +
        '\n')

        # Construct forward detector definitions string
        idf_forward_dets = idf_forward_det_head + idf_forward_locations + idf_forward_det_end
        return idf_forward_dets


#===========================================================================================================================================#
    def _construct_idf_back(self):
        # Define the back scattering detectors
        idf_back_head =(
        '   <type name="back">\n')

        # Define all locations of back scattering detectors
        idf_back_locations = ''
        for idx in range (3,135):
            # Get corresponding detector from IP file
            ip_det = self._ip_lines[idx].split()
            idf_back_locations += ('    <component type="pixel back" >\n')
            idf_back_locations += ('      <location t="%f" r="%s" p="%f" name="S%s"/>\n' % (float(ip_det[2]), ip_det[5], self._phi_data[idx - 3], ip_det[1]))
            idf_back_locations += ('      <parameter name="t0"> <value val="%f"/> </parameter>\n' % (float(ip_det[3])))
            idf_back_locations += ('    </component>\n')

        idf_back_end = (
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
