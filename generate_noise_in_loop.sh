for i in {5880..6000}
do
	echo $i
	cd noisefiles/hdfonly
	pycbc_condition_strain --gps-start-time 0 --gps-end-time 128 --sample-rate 2048 --fake-strain aLIGOZeroDetHighPower --fake-strain-seed $i --fake-strain-flow 20.0 --output-strain-file noise_$i.hdf
	cd ../../
	python hdf_to_txt.py $i
done

