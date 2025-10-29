import os
import shutil
import zipfile

base_dir = "/scratch/users/gelnesr/afcluster"
os.chdir(base_dir)

# Step 1: zip each inner folder
for folder in os.listdir(base_dir):
    folder_path = os.path.join(base_dir, folder)
    if os.path.isdir(folder_path):
        zip_path = f"{folder}.zip"
        print(f"Zipping {folder_path} → {zip_path}")
        shutil.make_archive(folder, 'zip', folder_path)

# Step 2: zip all the generated zip files into one big zip
final_zip = "28oct25_afcluster_ge.zip"
with zipfile.ZipFile(final_zip, 'w', zipfile.ZIP_DEFLATED) as zf:
    for z in os.listdir(base_dir):
        if z.endswith(".zip") and z != final_zip:
            zf.write(z)
print(f"Created {final_zip}")