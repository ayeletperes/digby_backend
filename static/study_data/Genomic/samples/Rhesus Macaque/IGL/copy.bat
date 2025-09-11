cp db.sqlite3 "../../../../digby_backend/static/study_data/Genomic/db/Rhesus Macaque/IGL/."
cp db_description.txt "../../../../digby_backend/static/study_data/Genomic/db/Rhesus Macaque/IGL/."
rm -rf "../../../../digby_backend/static/study_data/Genomic/samples/Rhesus Macaque/IGL/samples/*"
cp -r samples/* "../../../../digby_backend/static/study_data/Genomic/samples/Rhesus Macaque/IGL/."

