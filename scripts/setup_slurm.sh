set -e

ml gcc/12.4.0

mkdir -p $SCRATCH/tools
cd $SCRATCH/tools
#wget https://raw.githubusercontent.com/YoshitakaMo/localcolabfold/main/install_colabbatch_linux.sh
#bash install_colabbatch_linux.sh
export PATH="$SCRATCH/tools/localcolabfold/colabfold-conda/bin:$PATH"

echo "ColabFold environment prepared."

echo "Making uv enviornment"

VENV_DIR="afc"
if [ ! -d "$VENV_DIR" ]; then
    echo "Environment not found. Creating new uv environment..."
    uv venv "$VENV_DIR"
else
    echo "Environment already exists. Skipping creation."
fi

source afc/bin/activate
uv pip install numpy pandas biopython scikit-learn tqdm requests pyyaml
