
import mantid.simpleapi as api
from mantid.kernel import *
from mantid.api import *

IP_PROP = "IP File"
IDF_PROP = "IDF File"

class CreateVesuvioIDF(PythonAlgorithm):

    def summary(self):
        return "Uses an IP file to create an IDF for the Vesuvio instrument at ISIS."

    def PyInit(self):
        self.declareProperty(FileProperty(IP_PROP, "", extensions["par"]),
                              doc="The file path to the Vesuvio IP file.")

        self.declareProperty(FileProperty(IDF_PROP, "", Direction.Output),
                              doc="The file path to the desired IDF file.")


    def PyExec(self):
        ip_filename = self.getPropertyValue(IP_PROP)
        idf_filename = self.getPropertyValue(IDF_PROP)
        ip_lines = _read_ip(ip_filename)

        # Read headings
        ip_heading = ip_lines[0].split()

        # Read Monitors
        ip_monitor_one = ip_lines[1].split()
        ip_monitor_two = ip_lines[2].split()

        # Construct full IDF
        idf_string = _construct_idf_header()
        idf_string += _construct_idf_source_sample(ip_lines[3].split()[4])
        idf_string += _construct_idf_monitors(ip_monitor_one, ip_monitor_two)
        idf_string += _construct_idf_foils()
        idf_string += _construct_idf_forward(ip_lines)
        idf_string += _construct_idf_back(ip_lines)
        idf_string += _construct_idf_pixel()
        idf_string += _construct_idf_det_list()

        if write_flag:
            # Write to file
            _write_idf(filename, idf_string)


#===========================================================================================================================================#
    def _read_ip(self, ip_file_name):
        # Read IP file
        ip_file = open(ip_file_name, 'r')
        ip_lines = ip_file.readlines()
        return ip_lines


#===========================================================================================================================================#
    def _construct_idf_header(self):
        # Define the header of the IDF - input=()
        idf_header = ("""<?xml version="1.0" encoding="UTF-8"?>
        <!-- For help on the notation used to specify an Instrument Definition File
             see http://www.mantidproject.org/IDF -->
        <instrument xmlns="http://www.mantidproject.org/IDF/1.0"
                    xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
                    xsi:schemaLocation="http://www.mantidproject.org/IDF/1.0 http://schema.mantidproject.org/IDF/1.0/IDFSchema.xsd"
                    name="VESUVIO"
                    valid-from="2006-08-11 00:00:00"
                    valid-to="2100-01-31 23:59:59"
                    last-modified="2011-07-20 00:00:21">

        <defaults>
            <length unit="meter"/>
            <angle unit="degree"/>
            <reference-frame>
              <!-- The z-axis is set parallel to and in the direction of the beam. the
                   y-axis points up and the coordinate system is right handed. -->
              <along-beam axis="z"/>
              <pointing-up axis="y"/>
              <handedness val="right"/>
            </reference-frame>
        </defaults>

        """)

        return idf_header


#===========================================================================================================================================#
    def _construct_idf_source_sample(self, non_monitor_L0):
        # Define the source and sample positions - input=(non-monitor['L0'])
        idf_src_sam = ("""
          <!--  SOURCE AND SAMPLE POSITION -->

          <component type="moderator">
            <location z="-%s" />
          </component>

          <type name="moderator" is="Source" />

          <component type="sample-position">
            <location />
          </component>

          <type name="sample-position" is="SamplePos">
          </type>

        """ % non_monitor_L0)

        return idf_src_sam


#===========================================================================================================================================#
    def _construct_idf_monitors(self, monitor_one, monitor_two):
        # Define the Monitors - input=(monitor1['L1'], monitor2['L1'])
        idf_monitors = ("""
          <!-- MONITORS -->

          <component type="monitor" idlist="monitors">
            <location z="%s" name="Monitor 1"/>
            <location z="%s" name="Monitor 2"/>
          </component>

          <type name="monitor" is="monitor">
           <cuboid id="shape">
              <left-front-bottom-point x="0.01" y="-0.01" z="0.0"  />
              <left-front-top-point  x="0.01" y="-0.01" z="0.02"  />
              <left-back-bottom-point  x="-0.01" y="-0.01" z="0.0"  />
              <right-front-bottom-point  x="0.01" y="0.01" z="0.0"  />
           </cuboid>
          </type>

        """ % (monitor_one[5], monitor_two[5]))

        return idf_monitors


#===========================================================================================================================================#
    def _construct_idf_foils(self):
        # Define the Foils and Foil changers - input=()
        idf_foils = ("""
          <!-- FOILS  -->

          <!-- midpoint positions
               phi is simply set to put them all in the scattering plane
          -->
          <component name="foil-pos0" type="foil">
           <!-- position 0 -->
           <location t="-42" r="0.225" p="0"/>
           <location t="-64" r="0.225" p="0"/>
           <location t="35" r="0.225" p="0"/>
           <location t="55" r="0.225" p="0"/>
          </component>

          <component name="foil-pos1" type="foil">
           <!-- position 1 -->
           <location t="-31" r="0.225" p="0"/>
           <location t="-50" r="0.225" p="0"/>
           <location t="46" r="0.225" p="0"/>
           <location t="67" r="0.225" p="0"/>
          </component>

          <!-- The foil shape is currently set to give the correct dTheta as seen
               from the sample position -->
          <type name="foil">
            <cuboid id="foil-shape-shape">
              <left-front-bottom-point x="-0.02000121679" y="0.001" z="-0.001"/>
              <left-front-top-point x="-0.02000121679" y="0.001" z="0.001"/>
              <left-back-bottom-point x="-0.02000121679" y="-0.001" z="-0.001" />
              <right-front-bottom-point x="0.02000121679" y="0.001" z="-0.001" />
            </cuboid>
          </type>

          <!-- FOIL CHANGER -->
          <!-- The shape & location are made up to provide inputs for the
               CalculateGammaBackground algorithm -->
          <component type="foil-changer">
          <location x="0" y ="0" z="0"/>
          </component>

          <!-- Important quantity currently is the height. It defines the Z range
               for the integration of the foil positions  (-0.2,0.2)-->
          <type name="foil-changer">
            <cylinder id="foil-changer-shape">
            <centre-of-bottom-base x="0.0" y="-0.2" z="0.0" />
            <axis x="0.0" y="1" z="0" />
            <radius val="0.1" />
            <height val="0.4" />
            </cylinder>
          </type>

        """)
        return idf_foils


#===========================================================================================================================================#
    def _construct_idf_forward(self, ip_lines):
        # Define the forward scattering detectors
        idf_forward_det_head = ("""
          <!-- DETECTORS -->

          <component type="forward" idlist="forward">
            <location />
          </component>


          <component type="back" idlist="back">
            <location />
          </component>

          <type name="forward">
            <component type="pixel forward" >
        """)

        # Define all locations of forward scattering detectors - inputs (ip_det['theta'], ip_det['L1'], ???)
        idf_forward_locations = ''
        for idx in range(136,  199):
            # Get corresponding detector from IP file
            ip_det = ip_lines[idx].split()
            idf_forward_locations += ('      <location t="-%s" r="%s" p="%i"/>\n' % (ip_det[2], ip_det[5], 0))

        idf_forward_det_end = ("""    </component>
          </type>

        """)

        # Construct forward detector definitions string
        idf_forward_dets = idf_forward_det_head + idf_forward_locations + idf_forward_det_end
        return idf_forward_dets


#===========================================================================================================================================#
    def _construct_idf_back(self, ip_lines):
        # Define the back scattering detectors
        idf_back_head =("""  <type name="back">
            <component type="pixel back" >
        """)

        # Define all locations of forward scattering detectors - inputs (ip_det['theta'], ip_det['L1'], ???, ip_det['det'] )
        idf_back_locations = ''
        for idx in range (4,135):
            # Get corresponding detector from IP file
            ip_det = ip_lines[idx].split()
            idf_back_locations += ('      <location t="-%s" r="%s" p="%i" name="S%s"/>\n' % (ip_det[2], ip_det[5], 0, ip_det[1]))

        idf_back_end = ("""    </component>
          </type>
        """)

        # Construct back detector definitions string
        idf_back_dets = idf_back_head + idf_back_locations + idf_back_end
        return idf_back_dets


#===========================================================================================================================================#
    def _construct_idf_pixel(self):
        # Define pixel forward/back
        idf_pix = ("""
          <type name="pixel forward" is="detector">
            <cuboid id="shape">
              <left-front-bottom-point x="0.0125" y="-0.0395" z= "0.0045" />
              <left-front-top-point x="0.0125" y="0.0395" z= "0.0045" />
              <left-back-bottom-point x="0.0125" y="-0.0395" z= "-0.0045" />
              <right-front-bottom-point x="-0.0125" y="-0.0395" z= "0.0045" />
            </cuboid>
            <algebra val="shape" />
          </type>

          <type name="pixel back" is="detector">
            <cuboid id="shape" >
              <left-front-bottom-point x="0.0125" y="-0.04" z= "0.0045" />
              <left-front-top-point x="0.0125" y="0.04" z= "0.0045" />
              <left-back-bottom-point x="0.0125" y="-0.04" z= "-0.0045" />
              <right-front-bottom-point x="-0.0125" y="-0.04" z= "0.0045" />
            </cuboid>
          <algebra val="shape" />
          </type>

        """)
        return idf_pix


#===========================================================================================================================================#
    def _construct_idf_det_list(self):
        # Detector ID Lists + tail
        idf_det_id_list = ("""
          <!-- DETECTOR ID LISTS -->

          <idlist idname="back">
            <id start="2101" end="2144" />
            <id start="2201" end="2244" />
            <id start="2301" end="2344" />
          </idlist>

          <idlist idname="monitors">
            <id start="1101" end="1102" />
          </idlist>

          <idlist idname="forward">
            <id start="3101" end="3132" />
            <id start="3201" end="3232" />
          </idlist>

        </instrument>
        """)
        return idf_det_id_list


#===========================================================================================================================================#
    def _write_idf(idf_filename, idf_string):
        # Write IDF out to .xml file
        idf = open(idf_filename, 'w')
        idf.write(idf_string)
        idf.close()

# Register algorithm with Mantid
AlgorithmFactory.subscribe(CreateVesuvioIDF)
