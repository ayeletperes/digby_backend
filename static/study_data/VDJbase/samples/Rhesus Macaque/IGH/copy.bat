cp db.sqlite3 "../../../../digby_backend/static/study_data/VDJbase/db/Rhesus Macaque/IGH/."
cp db_description.txt "../../../../digby_backend/static/study_data/VDJbase/db/Rhesus Macaque/IGH/."
rm -rf "../../../../digby_backend/static/study_data/VDJbase/samples/Rhesus Macaque/IGH/samples/*"
cp -r samples/* "../../../../digby_backend/static/study_data/VDJbase/samples/Rhesus Macaque/IGH/."

