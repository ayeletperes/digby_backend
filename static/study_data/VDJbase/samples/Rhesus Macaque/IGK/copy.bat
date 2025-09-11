cp db.sqlite3 "../../../../digby_backend/static/study_data/VDJbase/db/Rhesus Macaque/IGK/."
cp db_description.txt "../../../../digby_backend/static/study_data/VDJbase/db/Rhesus Macaque/IGK/."
rm -rf "../../../../digby_backend/static/study_data/VDJbase/samples/Rhesus Macaque/IGK/samples/*"
cp -r samples/* "../../../../digby_backend/static/study_data/VDJbase/samples/Rhesus Macaque/IGK/."

