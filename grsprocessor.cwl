# Based on OGC best practises (8.4.  Application)
cwlVersion: v1.2
$namespaces:
  s: https://schema.org/
  js: https://json-schema.org/
$schemas:
- http://schema.org/version/latest/schemaorg-current-http.rdf
$graph:
- class: Workflow
  id: main
  label: Antflow CWL Bencher
  doc: Agregate metrics and generate insight on CWL workflow run
  requirements:
    InlineJavascriptRequirement: {}
  inputs:
    bucket_name:
      type: string
      doc: S3 bucket holding input_file, cams_file and dem_file
    input_file:
      type: string
      doc: Path (key) of the L1C SAFE product within bucket_name
    cams_file:
      type: ["null", string]
      doc: Path (key) of the CAMS netCDF file within bucket_name
    resolution:
      type: ["null", int]
      doc: spatial resolution of the scene pixels
    no_clobber:
      type: boolean
      doc: Do not process <input_file> if <output_file> already exists.
      default: True
    dem_file:
      type: ["null", string]
      doc: Path (key) of the DEM geotiff file within bucket_name (already subset for the S2 tile)
  steps:
    retrieve_input:
      run: '#download_dir_from_s3'
      in:
        bucket_name: bucket_name
        file_path: input_file
      out:
        [data]
    retrieve_cams:
      run: '#download_file_from_s3'
      when: $(inputs.file_path !== null)
      in:
        bucket_name: bucket_name
        file_path: cams_file
      out:
        [data]
    retrieve_dem:
      run: '#download_file_from_s3'
      when: $(inputs.file_path !== null)
      in:
        bucket_name: bucket_name
        file_path: dem_file
      out:
        [data]
    bench_run:
      run: '#bench_chain'
      in:
        input_file: retrieve_input/data
        cams_file: retrieve_cams/data
        resolution: resolution
        no_clobber: no_clobber
        dem_file: retrieve_dem/data
      out:
        [nc]
  outputs:
    nc:
      type: File[]
      outputSource: bench_run/nc

- class: CommandLineTool
  id: download_file_from_s3
  doc: |
    Download a single object from bucket_name/file_path as a staged local File.
  baseCommand: download_from_s3
  arguments: ["--out", $(runtime.outdir)/staged]
  hints:
    DockerRequirement:
      dockerPull: artifactory.cnes.fr/platform-incub-docker-prod-local/cwl-processing/cwl_helpers:v1.3.3
    ResourceRequirement:
      ramMin: 200
  requirements:
    NetworkAccess:
      networkAccess: true
    ToolTimeLimit:
      timelimit: 3600
  inputs:
    bucket_name:
      type: string
      inputBinding:
        prefix: "--bucket"
        position: 1
    file_path:
      type: string
      inputBinding:
        prefix: "--path"
        position: 2
  outputs:
    data:
      type: File
      outputBinding:
        glob: $(runtime.outdir)/staged/*

- class: CommandLineTool
  id: download_dir_from_s3
  doc: |
    Download an S3 prefix (directory) from bucket_name/file_path as a staged local Directory.
    baseCommand: download_from_s3
  arguments: ["--out", $(runtime.outdir)/staged]
  hints:
    DockerRequirement:
      dockerPull: artifactory.cnes.fr/platform-incub-docker-prod-local/cwl-processing/cwl_helpers:v1.3.3
    ResourceRequirement:
      ramMin: 200
  requirements:
    NetworkAccess:
      networkAccess: true
    InlineJavascriptRequirement: {}
    ToolTimeLimit:
      timelimit: 3600
  inputs:
    bucket_name:
      type: string
      inputBinding:
        prefix: "--bucket"
        position: 1
    file_path:
      type: string
      doc: Bucket-relative prefix. A trailing "/" is required for directory downloads and is added automatically if missing.
      inputBinding:
        prefix: "--path"
        position: 2
        valueFrom: '$(self.endsWith("/") ? self : self + "/")'
  outputs:
    data:
      type: Directory
      outputBinding:
        glob: $(runtime.outdir)/staged/

- class: CommandLineTool
  id: bench_chain
  baseCommand: "grs"
  arguments: ["--odir", $(runtime.outdir)]
  hints:
    DockerRequirement:
      dockerPull: guillaumeeb/grs:2.1.9
    ResourceRequirement:
      ramMax: 10000
  requirements:
    NetworkAccess:
      networkAccess: true
    ToolTimeLimit:
      timelimit: 43200
  inputs:
    input_file:
      type: Directory
      inputBinding:
        position: 1
    cams_file:
      type: ["null", File]
      inputBinding:
        prefix: "--cams_file"
        position: 2
    resolution:
      type: ["null", int]
      inputBinding:
        prefix: "--resolution"
        position: 3
    no_clobber:
      type: boolean
      inputBinding:
        prefix: "--no_clobber"
        position: 4
    dem_file:
      type: ["null", File]
      inputBinding:
        prefix: "--dem_file"
        position: 5
  outputs:
    nc:
      type: File[]
      outputBinding:
        glob: $(runtime.outdir)/**/*.nc
