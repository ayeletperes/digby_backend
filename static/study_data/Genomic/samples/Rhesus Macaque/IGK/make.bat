rm samples*.zip 
rm -rf samples/*
python ../../../../digby_backend/make_genomic_db.py "Rhesus Macaque" IGK
cd samples
rm ../samples.z*
7z a -tzip ..\samples.zip *
cd ..
#python ../../../../digby_data/python/commit_sample_file.py
