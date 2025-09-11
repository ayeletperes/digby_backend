cp db.sqlite3 "../../../../digby_backend/static/study_data/Genomic/db/Rhesus Macaque/IGK/."
cp db_description.txt "../../../../digby_backend/static/study_data/Genomic/db/Rhesus Macaque/IGK/."
rm -rf "../../../../digby_backend/static/study_data/Genomic/samples/Rhesus Macaque/IGK/samples/*"
cp -r samples/* "../../../../digby_backend/static/study_data/Genomic/samples/Rhesus Macaque/IGK/."

