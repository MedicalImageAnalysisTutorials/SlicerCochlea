# this link the installed extension source file to this repository.

ewd1=$HOME"/.config/NA-MIC/Extensions-28283"
ewd2="lib/Slicer-4.11/qt-scripted-modules"
ghb=""

#Cochlea
ext="SlicerCochlea"
mv $ewd1/$ext/$ewd2/CochleaReg.py $ewd1/$ext/$ewd2/CochleaReg.py.bk
ln -s $PWD/CochleaReg/CochleaReg.py $ewd1/$ext/$ewd2/CochleaReg.py  
mv $ewd1/$ext/$ewd2/CochleaSeg.py $ewd1/$ext/$ewd2/CochleaSeg.py.bk
ln -s $PWD/CochleaSeg/CochleaSeg.py $ewd1/$ext/$ewd2/CochleaSeg.py  
mv $ewd1/$ext/$ewd2/VisSimCommon.py $ewd1/$ext/$ewd2/VisSimCommon.py.bk
ln -s $PWD/VisSimCommon/VisSimCommon.py $ewd1/$ext/$ewd2/VisSimCommon.py

