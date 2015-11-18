
import mantid.simpleapi as api
from mantid.kernel import *
from mantid.api import *

import os

IP_PROP = "IPFilename"
IDF_PROP = "OutputFilename"

class CreateVesuvioIDF(PythonAlgorithm):

    _phi_data = ['-34.73019', '-28.51712', '-22.01023', '-15.42126', '-8.96977', '-2.84422', '-31.24211', '-24.71659',
                      '-17.98582', '-11.28699', '-4.84608', '-33.93161', '-27.46684', '-20.65919', '-13.75649', '-7.01763',
                      '-30.21015', '-23.40558', '-16.36037', '-9.35697', '-26.17764', '-19.06803', '-11.85159', '-4.82908',
                      '-28.91532', '-21.83579', '-14.48041', '-7.17771', '-24.60827', '-17.20561', '-9.68561', '-19.97559',
                      '-12.32311', '-4.73797', '-22.72612', '-15.04853', '-7.24444', '-17.80136', '-9.87799', '-12.58767',
                      '-4.53313', '-7.15059', '-9.83415', '-4.17612', '-41.93841', '-43.30651', '-43.92399', '-43.68669',
                      '-42.63518', '-40.93858', '-40.95912', '-41.96552', '-42.15931', '-41.50458', '-40.11456', '-38.36715',
                      '-39.70691', '-40.31522', '-40.08216', '-39.04887', '-37.17761', '-38.15671', '-38.34420', '-37.70683',
                      '-35.70139', '-36.28626', '-36.06160', '-35.06953', '-32.98778', '-33.91587', '-34.09572', '-33.49011',
                      '-31.26197', '-31.80930', '-31.59777', '-29.22166', '-29.38656', '-28.83115', '-26.37644', '-26.86657',
                      '-26.67748', '-24.07138', '-24.21486', '-21.46766', '-21.30824', '-18.60336', '-15.65690', '-12.61408',
                      '-5.20558', '-11.59581', '-18.22201', '-24.87393', '-31.32516', '-37.37588', '-7.42193', '-14.10178',
                      '-20.92907', '-27.65962', '-34.06129', '-3.14545', '-9.80653', '-16.74500', '-23.71478', '-30.44977',
                      '-5.37521', '-12.35012', '-19.49746', '-26.53536', '-7.78097', '-15.03084', '-22.31755', '-29.33578',
                      '-3.09051', '-10.34895', '-17.81385', '-25.15397', '-5.50506', '-13.05341', '-20.65000', '-8.08435',
                      '-15.84983', '-23.47730', '-2.97370', '-10.79530', '-18.68213', '-5.55191', '-13.58835', '-8.25747',
                      '-16.39950', '-11.03355', '-5.46120', '-8.21321', '34.55456', '26.67983', '18.55840', '10.50738',
                      '-32.72541', '-24.62692', '-16.39836', '-8.37096', '30.46160', '24.58116', '18.61591', '12.69170',
                      '-29.01630', '-23.06463', '-17.07993', '-11.19273', '34.26094', '26.39880', '18.30117', '10.28211',
                      '-32.75540', '-24.68920', '-16.50932', '-8.53851', '30.20659', '24.29702', '18.31852', '12.39937',
                      '-29.95603', '-24.02064', '-18.02379', '-12.09454', '30.64800', '24.73324', '18.74552', '12.81407',
                      '-29.99957', '-24.03972', '-18.01542', '-12.05977', '34.15544', '26.18005', '17.99844', '9.93483',
                      '-33.64568', '-25.57022', '-17.31691', '-9.21912', '30.07557', '24.19055', '18.22387', '12.30471',
                      '-29.82889', '-23.87732', '-17.86838', '-11.93202', '33.38143', '25.41904', '17.26605', '9.24428',
                      '-33.67272', '-25.69463', '-17.55219', '-9.56592']


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
            idf_forward_locations += ('      <location t="-%s" r="%s" p="%s" name="S%s"/>\n' % (ip_det[2], ip_det[5], self._phi_data[idx-3], ip_det[1]))

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

        # Define all locations of forward scattering detectors
        idf_back_locations = ''
        for idx in range (3,135):
            # Get corresponding detector from IP file
            ip_det = ip_lines[idx].split()
            idf_back_locations += ('      <location t="-%s" r="%s" p="%s" name="S%s"/>\n' % (ip_det[2], ip_det[5], self._phi_data[idx - 3], ip_det[1]))

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
