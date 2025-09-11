cp db.sqlite3 ../../../../digby_backend/static/study_data/Genomic/db/Human/IGH/.
cp db_description.txt ../../../../digby_backend/static/study_data/Genomic/db/Human/IGH/.
rm -rf ../../../../digby_backend/static/study_data/Genomic/samples/Human/IGH/samples/*
cp -r samples/* ../../../../digby_backend/static/study_data/Genomic/samples/Human/IGH/.
