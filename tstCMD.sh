clear
ExtVer="4.10"
ExtID="28257"
Ext="SlicerCochlea"
Module="CochleaSeg"
SlicerP=" /home/"$USER"/sw/Slicer-4.10.2/Slicer "
addSettingPath=" /home/"$USER"/myGitHub/SlicerCochlea/zbuild-"$ExtVer"/AdditionalLauncherSettings.ini"
moduleSPath=" /home/"$USER"/myGitHub/SlicerCochlea/zbuild-"$ExtVer"/lib/Slicer-"$ExtVer"/qt-scripted-modules "
moduleCPath=" /home/"$USER"/myGitHub/SlicerCochlea/zbuild-"$ExtVer"/lib/Slicer-"$ExtVer"/cli-modules "
moduleLPath=" /home/"$USER"/myGitHub/SlicerCochlea/zbuild-"$ExtVer"/lib/Slicer-"$ExtVer"/qt-loadable-modules "
#tstPath=" /home/"$USER"/.config/NA-MIC/Extensions-"$ExtID"/"$Ext
tstSrcPath="'"/home/$USER/myGitHub/SlicerCochlea/$ExtVer/$Module"'"
tstBldPath="'"/home/$USER/myGitHub/SlicerCochlea/zbuild-$ExtVer/$Module"'"
echo $tstSrcPath
echo $tstBldPath

$SlicerP --no-splash --testing --launcher-additional-settings  $addSettingPath --additional-module-paths $moduleSPath  $moduleCPath $moduleLPath   --python-code "import slicer.testing; slicer.testing.runUnitTest(["$tstBldPath","$tstSrcPath"],'"$Module"')"




#D:\D\S\Slicer-4102-build\Slicer-build\Slicer.exe "--no-splash" "--testing" "--launcher-additional-settings" "D:/D/S/S-4102-E-b/SlicerCochlea-build/AdditionalLauncherSettings.ini" "--additional-module-paths" "D:/D/S/S-4102-E-b/SlicerCochlea-build/lib/Slicer-4.10/qt-scripted-modules" "D:/D/S/S-4102-E-b/SlicerCochlea-build/lib/Slicer-4.10/cli-modules" "D:/D/S/S-4102-E-b/SlicerCochlea-build/lib/Slicer-4.10/qt-loadable-modules" "--python-code" "import slicer.testing; slicer.testing.runUnitTest(['D:/D/S/S-4102-E-b/SlicerCochlea-build/CochleaSeg', 'D:/D/S/S-4102-E-b/SlicerCochlea/CochleaSeg'], 'CochleaSeg')"
