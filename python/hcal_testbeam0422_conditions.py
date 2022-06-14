"""Load the April 2022 Testbeam conditions for HCal

We assume that the conditions-data repository has been cloned into LDMX_BASE.

    cd ${LDMX_BASE}
    git clone https://github.com/LDMX-Software/conditions-data

You can modify this path if you have your conditions somewhere else.
**Make sure these conditions are available to the container.

    import hcal_testbeam0422_conditions
    hcal_testbeam0422_conditions.adc_pedestal.conditions_baseURL = 'file:///full/path/my/url/'

"""

from LDMX.Conditions.SimpleCSVTableProvider import SimpleCSVIntegerTableProvider, SimpleCSVDoubleTableProvider
import os

adc_pedestal = SimpleCSVDoubleTableProvider("hcal_adc_pedestal",["PEDESTAL_ADC","PEDESTAL_RMS_ADC"])
adc_pedestal.conditions_baseURL = f'file://{os.environ["LDMX_BASE"]}/conditions-data/'
adc_pedestal.entriesURL = '${LDMX_CONDITION_BASEURL}/Hcal/testbeam04-2022/pedestals/index_v1_0_0.csv'

adc_gain = SimpleCSVDoubleTableProvider("hcal_adc_gain",["gain"])
adc_gain.validForAllRows([1.2]) # 4 ADCs per PE - maxADCRange/readoutPadCapacitance/1024

tot_pedestal = SimpleCSVDoubleTableProvider("hcal_tot_pedestal",["pedestal"])
tot_pedestal.validForAllRows([1.]) # dummy value since TOT is not implemented

tot_gain = SimpleCSVDoubleTableProvider("hcal_tot_gain",["gain"])
tot_gain.validForAllRows([2.5]) # dummy value - conversion to estimated charge deposited in TOT mode

from .conditions import HcalReconConditionsProvider
HcalReconConditionsProvider(adc_pedestal, adc_gain, tot_pedestal, tot_gain)
