::Adjust the project path
cd ..\VS2010\LocalWaveletPCA

call Release\LocalWavletPCA.exe ..\..\Model ..\..\Model\MeanFace_NE00_Scaled.off ..\..\Model\ScaledAlignedLandmarks.txt ..\..\Example\stereo_pointcloud.off ..\..\Example\stereo_pointcloud_landmarks.txt ..\..\Example\stereo_pointcloud_fitting

pause