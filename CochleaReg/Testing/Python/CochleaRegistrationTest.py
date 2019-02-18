# run:  ~/sw/Slicer-4.10.0/Slicer --no-main-window --python-script  CochleaRegistrationTest.py
#TODO when called from test, no display. add number instead of boolean 0:internal, 1, external, 2 test


import CochleaReg

CT = CochleaReg.CochleaRegTest()
CT.testSlicerCochleaRegistration()
