#!/bin/bash
#
#SBATCH --job-name=grs-docker-tests
#SBATCH --output=
#SBATCH --time=00:30:00
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem-per-cpu=4G
#SBATCH --account=campus

#cnes_level2

handle_test_result () {
  test_name=${1}
  test_result=${2}
  test_time=${3}
  test_params=${4}

  if [[ ${test_result} -ne 0 ]]
  then
    echo "<testcase classname=\"${test_name}\" name=\"${test_name}\" time=\"${test_time}\">" >> ${docker_report}
    echo "<failure type=\"Failed\"> GRS launcher exited with status ${test_result}. See log for more info. </failure>" >> ${docker_report}
    echo "</testcase>" >> ${docker_report}
    export test_failed="fail"
    return
  fi

  _OUT=$(grep output_dir ${test_params} | cut -d '"' -f4)
  suffix=$(grep suffix ${test_params} | cut -d '"' -f4)
  OUT_DIR=$(eval echo ${_OUT})
  OUT_REF=${REF_DIR}/tests/docker_tests/
  mkdir -p "${OUT_REF}"

  test_ok="yes"
  nout=0
  for out in $(find "${OUT_DIR}" -type f -name "*${suffix}*.nc" 2>/dev/null); do
    nout=$((nout+1))
    ref=$(find "${OUT_REF}" -iname "$(basename ${out})" 2>/dev/null)
    if [[ -z "${ref}" ]]; then
      echo "[REF GEN] No reference for $(basename ${out}); creating from current output"
      cp "${out}" "${OUT_REF}/$(basename ${out})"
      continue
    fi
    echo "grs_compare_ref --nc ${out} --nc-ref ${ref}"
    grs_compare_ref --nc "${out}" --nc-ref "${ref}" &> "reports/ref_diff_docker_$(basename ${out} .nc).log"
    ref_status=${?}
    if [[ ${ref_status} -ne 0 ]]
    then
      echo "[REF DIFF] $(basename ${out}) differs from reference"
      test_ok="no"
    fi
  done

  if [[ ${test_ok} == "yes" ]]
  then
    echo "<testcase classname=\"${test_name}\" name=\"${test_name}\" time=\"${test_time}\"/>" >> ${docker_report}
  else
    echo "<testcase classname=\"${test_name}\" name=\"${test_name}\" time=\"${test_time}\">" >> ${docker_report}
    echo "<failure type=\"Failed\"> Output differs from reference. See reports/ref_diff_docker.log </failure>" >> ${docker_report}
    echo "</testcase>" >> ${docker_report}
    export test_failed="fail"
  fi
}

cd ${MY_JOB_WORKDIR}

# Activate conda env for grs_compare_ref (reference comparison runs outside the container)
module load conda
conda activate /work/scratch/data/${USER}/.conda/envs/grs2_tests
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}"

# singularity is a module on the compute nodes (not in the default PATH)
module load singularity/3.10.2

# Locate the Singularity image (rsynced into grs2/ by the CI job)
SIF=$(ls ${MY_JOB_WORKDIR}/grs2/*.sif 2>/dev/null | head -n 1)
if [ -z "${SIF}" ]; then
  echo "ERROR: no .sif image found in ${MY_JOB_WORKDIR}/grs2/" >&2
  exit 1
fi
echo "Using SIF image: ${SIF}"

mkdir -p reports
export docker_report="reports/docker-report.xml"
export TEST_PARAMS=${MY_JOB_WORKDIR}/grs2/grs/tests/docker_tests/test_1_inputs.yml
echo "start docker CI"
echo "<testsuite tests=\"1\">" >> ${docker_report}

# /work holds the input product, the CAMS forecasts and the scratch output directory.
# --cleanenv + PYTHONNOUSERSITE=1 prevent the host's ~/.local user-site (bind-mounted
# via $HOME) from shadowing the container's packages.
# USER and MY_JOB_WORKDIR are passed explicitly so os.path.expandvars() inside
# grs_launcher can resolve $USER/$MY_JOB_WORKDIR in the YAML config paths.
start=$(date +%s)
singularity exec --cleanenv --bind /work ${SIF} env PYTHONNOUSERSITE=1 PYTHONPATH= USER=${USER} MY_JOB_WORKDIR=${MY_JOB_WORKDIR} grs_launcher ${TEST_PARAMS}
test_status=${?}
end=$(date +%s)
handle_test_result docker_test_nominal ${test_status} $((end-start)) ${TEST_PARAMS}

echo "</testsuite>" >> ${docker_report}

if [ -z ${test_failed+x} ]
then
  exit 0
else
  exit 1
fi
