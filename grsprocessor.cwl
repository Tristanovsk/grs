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
  inputs:
    input_file:
      type: string
      doc: File to process
    cams_file:
      type: ["null", string]
      doc: File to process
    resolution
      type: ["null", integer]
      doc: spatial resolution of the scene pixels
    no_clobber
      type: boolean
      doc: Do not process <input_file> if <output_file> already exists.
      default: True
    dem_file
      type: ["null", string]
      doc: Absolute path of the DEM geotiff file (already subset for the S2 tile)   
  steps:
    bench_run:
      run: '#bench_chain'
      in:
        input_file: input_file
        cams_file: cams_file
        resolution: resolution
        no_clobber: no_clobber
        dem_file: dem_file
      out:
        [nc]
  outputs:
    nc:
      type: File[]
      outputSource: bench_run/nc
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
  inputs:
    input_file:
      type: string
      inputBinding:
        position: 1
    cams_file:
      type: ["null", string]
      inputBinding:
        prefix: "--cams_file"
        position: 2
    resolution
      type: ["null", integer]
      inputBinding:
        prefix: "--resolution"
        position: 3
    no_clobber
      type: boolean
      inputBinding:
        prefix: "--no_clobber"
        position: 4
    dem_file
      type: ["null", string]
      inputBinding:
        prefix: "--dem_file"
        position: 5
  outputs:
    nc:
      type: File[]
      outputBinding:
        glob: $(runtime.outdir)/**/*.nc


