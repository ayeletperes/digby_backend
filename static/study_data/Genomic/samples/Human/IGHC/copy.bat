cp db.sqlite3 ../../../../digby_backend/static/study_data/Genomic/db/Human/IGHC/.
cp db_description.txt ../../../../digby_backend/static/study_data/Genomic/db/Human/IGHC/.
rm -rf ../../../../digby_backend/static/study_data/Genomic/samples/Human/IGHC/samples/*
cp -r samples/* ../../../../digby_backend/static/study_data/Genomic/samples/Human/IGHC/.
