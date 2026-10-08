#!/usr/bin/env bash
# =====================================================================
#  BioEmu sobre la quimera LBPCheY (modelo 0 de AlphaFold3)
#  Uso:   bash run_bioemu.sh
#  Ejecutar desde cualquier lado; las rutas son absolutas.
# =====================================================================
set -eo pipefail

FOLD=~/fold_lbpchey                         # prediccion de AlphaFold3
MODEL=$FOLD/fold_lbpchey_model_0.cif        # primer modelo
BASE=~/cursos/proteinas/quimerasLBPChey
OUT=$BASE/bioemu
ENV=bioemu
NSAMPLES=${NSAMPLES:-300}      # subelo a 1000 si tienes tiempo/VRAM
BATCH=${BATCH:-10}             # bajalo a 5 si te da CUDA out of memory

mkdir -p "$OUT"
cd "$BASE"

# ---------------------------------------------------------------------
# 1. Entorno
#    OJO: 'pip' en el PATH es ~/.local/bin/pip (python 3.12 del usuario),
#    por eso se usa siempre 'python -m pip' y se ignora el site de usuario.
# ---------------------------------------------------------------------
export PYTHONNOUSERSITE=1

# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
if ! conda env list | grep -qE "^$ENV\s"; then
  echo ">> creando entorno conda '$ENV'"
  conda create -y -n "$ENV" python=3.11
fi
conda activate "$ENV"

python -c "import bioemu, mdtraj, matplotlib" 2>/dev/null || {
  echo ">> instalando bioemu (esto tarda, baja ColabFold la primera vez)"
  python -m pip install --upgrade pip
  python -m pip install "bioemu[md]" mdtraj matplotlib   # [md]: cadenas laterales
}

nvidia-smi --query-gpu=name,memory.total --format=csv

# ---------------------------------------------------------------------
# 0. Secuencia: se extrae del modelo 0 de AlphaFold3
# ---------------------------------------------------------------------
LBPCHEY=$(python - "$MODEL" <<'EOF'
import sys
from Bio.PDB import MMCIFParser, PPBuilder
s = MMCIFParser(QUIET=True).get_structure("m", sys.argv[1])
print("".join(str(pp.get_sequence()) for pp in PPBuilder().build_peptides(s[0])))
EOF
)
echo ">> modelo: $MODEL"
echo ">> longitud LBPCheY: ${#LBPCHEY}  (debe decir 339)"
[ "${#LBPCHEY}" -eq 339 ] || { echo "ERROR: longitud incorrecta"; exit 1; }
printf ">LBPCheY_model0\n%s\n" "$LBPCHEY" > "$OUT/lbpchey_model0.fasta"

# ---------------------------------------------------------------------
# 2. Muestreo del ensamble
# ---------------------------------------------------------------------
echo ">> muestreando $NSAMPLES estructuras de LBPCheY"
python -m bioemu.sample \
    --sequence "$OUT/lbpchey_model0.fasta" \
    --num_samples "$NSAMPLES" \
    --batch_size_100 "$BATCH" \
    --output_dir "$OUT/lbpchey"

# ---------------------------------------------------------------------
# 3. Cadenas laterales + minimizacion local
#    (BioEmu entrega solo backbone)
# ---------------------------------------------------------------------
echo ">> reconstruyendo cadenas laterales"
python -m bioemu.sidechain_relax \
    --pdb-path  "$OUT/lbpchey/topology.pdb" \
    --xtc-path  "$OUT/lbpchey/samples.xtc" \
    --outpath   "$OUT/lbpchey" || \
  echo "!! sidechain_relax fallo; el analisis de backbone sigue siendo valido"

echo
echo "===================================================="
echo " Listo. Resultados en: $OUT/lbpchey"
echo " Ahora corre:  python $FOLD/analiza_ensamble.py"
echo "===================================================="
