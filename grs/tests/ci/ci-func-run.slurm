#!/bin/bash
#
#SBATCH --job-name=grs-functionnal-tests
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
    echo "<testcase classname=\"${test_name}\" name=\"${test_name}\" time=\"${test_time}\">" >> ${func_report}
    echo "<failure type=\"Failed\"> GRS processing failed. See log for more info. </failure>" >> ${func_report}
    echo "</testcase>" >> ${func_report}
    export test_failed="fail"
    return
  fi

  _OUT=$(grep output_dir ${test_params} |cut -d '"' -f4)
  suffix=$(grep suffix ${test_params} |cut -d '"' -f4)
  OUT_DIR=$(eval echo ${_OUT})
  OUT_REF=${REF_DIR}/tests/functionnal_tests/
  mkdir -p "${OUT_REF}"

  test_ok="yes"
  for out in $(find "${OUT_DIR}" -type f -name "*${suffix}*.nc"); do
    ref=$(find "${OUT_REF}" -iname "$(basename ${out})")
    if [[ -z "${ref}" ]]; then
      echo "[REF GEN] No reference for $(basename ${out}); creating from current output"
      cp "${out}" "${OUT_REF}/$(basename ${out})"
      continue
    fi
    echo "grs_compare_ref --nc ${out} --nc-ref ${ref}"
    grs_compare_ref --nc "${out}" --nc-ref "${ref}" &> "reports/ref_diff_${suffix}.log"
    ref_status=${?}
    if [[ ${ref_status} -ne 0 ]]
    then
      echo "[REF DIFF] $(basename ${out}) differs from reference"
      test_ok="no"
    fi
  done

  if [[ ${test_ok} == "yes" ]]
  then
    echo "<testcase classname=\"${test_name}\" name=\"${test_name}\" time=\"${test_time}\"/>" >> ${func_report}
  else
    echo "<testcase classname=\"${test_name}\" name=\"${test_name}\" time=\"${test_time}\">" >> ${func_report}
    echo "<failure type=\"Failed\"> Output differs from reference. See ref_diff_${suffix}.log </failure>" >> ${func_report}
    echo "</testcase>" >> ${func_report}
    export test_failed="fail"
  fi
}

export COVERAGE_RCFILE=${MY_JOB_WORKDIR}/grs2/pyproject.toml
ml conda
conda activate /work/scratch/data/${USER}/.conda/envs/grs2_tests
# Prepend conda env lib so its newer libstdc++ (CXXABI_1.3.15) is found
# before the older one from the spack/mamba module load
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}"

cd ${MY_JOB_WORKDIR}

mkdir -p reports
export func_report="reports/functionnal_tests-report.xml"
export TEST_DATA_PATH=${PWD}/tests/data
echo "start CI"
echo "<testsuite tests=\"5\">" >> ${func_report}

find /work/scratch/data/$USER/OBS2CO/CI/tests/coverage/ -type f -exec rm {} +

start=$(date +%s)
coverage run --data-file=/work/scratch/data/$USER/OBS2CO/CI/tests/coverage/test_1.cov \
  grs2/exe/launcher.py ${MY_JOB_WORKDIR}/grs2/grs/tests/functionnal_tests/test_1_inputs.yml
test_status=${?}
end=$(date +%s)
handle_test_result grs_nominal_1 ${test_status} $((end-start)) ${MY_JOB_WORKDIR}/grs2/grs/tests/functionnal_tests/test_1_inputs.yml

start=$(date +%s)
coverage run --data-file=/work/scratch/data/$USER/OBS2CO/CI/tests/coverage/test_2.cov \
  grs2/exe/launcher.py ${MY_JOB_WORKDIR}/grs2/grs/tests/functionnal_tests/test_2_inputs.yml
test_status=${?}
end=$(date +%s)
handle_test_result grs_nominal_without_mask ${test_status} $((end-start)) ${MY_JOB_WORKDIR}/grs2/grs/tests/functionnal_tests/test_2_inputs.yml

start=$(date +%s)
coverage run --data-file=/work/scratch/data/$USER/OBS2CO/CI/tests/coverage/test_3.cov \
  grs2/exe/launcher.py ${MY_JOB_WORKDIR}/grs2/grs/tests/functionnal_tests/test_3_inputs.yml
test_status=${?}
end=$(date +%s)
handle_test_result grs_nominal_CAMSv2 ${test_status} $((end-start)) ${MY_JOB_WORKDIR}/grs2/grs/tests/functionnal_tests/test_3_inputs.yml

start=$(date +%s)
grs /work/scratch/data/$USER/OBS2CO/CI/test_1/S2A_MSIL1C_20231012T104951_N0509_R051_T31TCJ_20231012T143114.SAFE \
   --cams_file /work/datalake/static_aux/CAMS/OBS2CO_GRS_v1/2023/10/12/2023-10-12-cams-global-atmospheric-composition-forecasts.nc \
   -o ${MY_JOB_WORKDIR}/CI/OUT --resolution 20 --suffix _cli_run
test_status=${?}
end=$(date +%s)
#handle_test_result grs_nominal_cli ${test_status} $((end-start)) # TODO handle test from cli

echo "</testsuite>" >> ${func_report}

# combine et create the coverage report
cd /work/scratch/data/$USER/OBS2CO/CI/tests/coverage/
coverage combine --keep *.cov
coverage xml -o ${MY_JOB_WORKDIR}/reports/cov_report.xml

if [ -z ${test_failed+x} ]
then
  exit 0
else
  exit 1
fi
