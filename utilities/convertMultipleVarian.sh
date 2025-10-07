mkdir nii_images
for f in fsems_*.img ; do
    
    A=$(cut -d'/' -f1 <<<$f)
    B=$(cut -d'/' -f3 <<<$f)
    C=$(cut -d'.' -f1 <<<$B)
    file="/home/MEDECINE/ricm2409/Bureau/treluc_20240423_003/nii_images/"$A"_"$C
    echo $file
    python /home/MEDECINE/ricm2409/Documents/Code/my_scripts/fdf2nii.py $f $file --scale 1
done
