cp db.sqlite3 ../../../../digby_backend/static/study_data/Genomic/db/Human/IGK/.
cp db_description.txt ../../../../digby_backend/static/study_data/Genomic/db/Human/IGK/.
rm -rf ../../../../digby_backend/static/study_data/Genomic/samples/Human/IGK/samples/*
cp -r samples/* ../../../../digby_backend/static/study_data/Genomic/samples/Human/IGK/.
