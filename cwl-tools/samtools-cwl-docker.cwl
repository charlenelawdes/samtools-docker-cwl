cwlVersion: v1.0
class: CommandLineTool
baseCommand: samtools
id: samtools-tool

inputs:
  input_bam:
    type: File
    inputBinding:
      position: 1
    doc: "Input BAM file"
  
  output_bam:
    type: File
  
  command:
    type: string
    inputBinding:
      position: 2
    doc: "Samtools subcommand"

outputs:
  - id: output_bam
    type: File
    outputBinding:
      glob: "$(inputs.command).bam"

requirements:
  - class: InlineJavascriptRequirement

hints:
  DockerRequirement:
    dockerPull: "ghcr.io/bwbioinfo/samtools-docker-cwl:latest"

