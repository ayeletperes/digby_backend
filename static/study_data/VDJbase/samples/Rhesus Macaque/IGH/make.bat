python ../../../../digby_backend/make_vdjbase_db.py Rhesus_Macaque IGH >make.log
cd samples
zip -r ../samples.zip *
cd ..
python ../../../../digby_data/python/commit_sample_file.py
