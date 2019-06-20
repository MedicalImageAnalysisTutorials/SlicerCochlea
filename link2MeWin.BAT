REM Create links to github sources:
Rem /home/ibr/.config/NA-MIC/Extensions-27931
REM open cmd as admin then run: mklink dest src
set ext="SlicerCochlea"
set vr="\lib\Slicer-4.10"
set p1="C:\Users\ibr\AppData\Roaming\NA-MIC\Extensions-28257\"
set p2="\qt-scripted-modules\"
set src="D:\myGitHub\"%ext%"\4.10\"

rename %p1%%ext%%vr%%p2%"CochleaSeg.py" "CochleaSeg.py.bk"
mklink %p1%%ext%%vr%%p2%"CochleaSeg.py" %src%"CochleaSeg""\CochleaSeg.py"

rename %p1%%ext%%vr%%p2%"CochleaReg.py" "CochleaReg.py.bk"
mklink %p1%%ext%%vr%%p2%"CochleaReg.py" %src%"CochleaReg""\CochleaReg.py"

rename %p1%%ext%%vr%%p2%"VisSimCommon.py" "VisSimCommon.py.bk"
mklink %p1%%ext%%vr%%p2%"VisSimCommon.py" %src%"VisSimCommon""\VisSimCommon.py"