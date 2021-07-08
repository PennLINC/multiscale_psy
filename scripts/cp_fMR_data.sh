NewSubjs=$(</cbica/projects/pinesParcels/data_psy/NewPsySubjs.txt)
for i in ${NewSubjs}; do
cp /cbica/projects/pncSingleFuncParcel/pncSingleFuncParcel_psycho/data/SurfaceData/CombinedData/${i} -r /cbica/projects/pinesParcels/data_psy/CombinedData
done 
