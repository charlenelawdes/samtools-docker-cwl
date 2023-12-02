cwlVersion: v1.0
class: CommandLineTool
baseCommand: samtools
id: samtools-tool

inputs:
  command:
    type: string
    inputBinding:
      position: 2
    doc: "Samtools subcommand"

outputs:
  test_out:
    type: stdout

hints:
  DockerRequirement:
    dockerPull: "ghcr.io/bwbioinfo/samtools:latest"
