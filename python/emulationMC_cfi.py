'''

Emulate the L1 and UCT upgrade primitives, and put them in the event.

For MC we need to add an extra step to work around an HCAL bug.  After the
hcalDigis, we need to re-emulate the TPGs from the full readout.
This wraps the normal from-data behaviour in emulator_cfi and adds MC specific
HCAL TPG emulation step.

Authors: M. Cepeda, T. Perry, E. Friis

'''

from L1Trigger.UCT2015.emulation_cfi import *  # NOQA
from L1Trigger.Configuration.ValL1Emulator_cff import *  # NOQA

print "Using workarounds for MC RAW data content bugs"

# We use valHcalTriggerPrimitiveDigis since we emulate these from the full
# hcalDigis readout.  Otherwise it depends on the presence of the sim digis.
uctDigiStep += valHcalTriggerPrimitiveDigis
hackHCALMIPs.src = cms.InputTag("valHcalTriggerPrimitiveDigis")
# Otherwise everything is zero.
HcalTPGCoderULUT.LUTGenerationMode = True
