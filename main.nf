#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def versionMessage()
{
  log.info"""

  nf-humann-regroup - Version: ${workflow.manifest.version}
  """.stripIndent()
}

def helpMessage()
{
  log.info"""

nf-humann-regroup - Version: ${workflow.manifest.version}

  Mandatory arguments:
    --utility_db path   folder for the HUMAnN utility mapping database
    --input      path   path to input list of CMD studies to regroup
    --groups     list   comma-delimited list of HUMAnN mapping groups to regroup into
    --split      bool   boolean on whether to output split stratified tables [true, false]
"""
}

/**
Prints version when asked for
*/
if (params.version) {
  versionMessage()
  exit 0
}

/**
Prints help when asked for
*/

if (params.help) {
  helpMessage()
  exit 0
}



// Header log info
log.info """---------------------------------------------
nf-humann-regroup
---------------------------------------------

Analysis introspection:

"""

def summary = [:]

summary['Starting time'] = new java.util.Date()
//Environment
summary['Environment'] = ""
summary['Pipeline Name'] = 'nf-humann-regroup'
summary['Pipeline Version'] = workflow.manifest.version

summary['Config Profile'] = workflow.profile
summary['Resumed'] = workflow.resume

summary['Nextflow version'] = nextflow.version.toString() + " build " + nextflow.build.toString() + " (" + nextflow.timestamp + ")"

summary['Java version'] = System.getProperty("java.version")
summary['Java Virtual Machine'] = System.getProperty("java.vm.name") + "(" + System.getProperty("java.vm.version") + ")"

summary['Operating system'] = System.getProperty("os.name") + " " + System.getProperty("os.arch") + " v" +  System.getProperty("os.version")
summary['User name'] = System.getProperty("user.name") //User's account name

summary['Container Engine'] = workflow.containerEngine
if(workflow.containerEngine) summary['Container'] = workflow.container
summary['HUMAnN'] = params.docker_container_humann3
summary['CMD'] = params.docker_container_cmd

//General
summary['Running parameters'] = ""
summary['Studies file'] = params.input
summary['Uniref database'] = params.utility_db
summary['Output groups'] = params.groups
summary['Split stratified'] = params.split
summary['Split size'] = params.split_size

//Folders
summary['Folders'] = ""
summary['Output dir'] = params.outdir

log.info summary.collect { k,v -> "${k.padRight(27)}: $v" }.join("\n")
log.info ""




/**
  Prepare workflow introspection

  This process adds the workflow introspection (also printed at runtime) in the logs
  This is NF-CORE code.
*/

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve("workflow_summary_mqc.yaml")
    yaml_file.text  = """
    id: 'workflow-summary'
    description: "This information is collected when the pipeline is started."
    section_name: 'nf-humann-regroup Workflow Summary'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd>$v</dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}


/////// new stuff


process download_cmd_tables {
  tag "$study"
  label "medmem"
  // conda "/home/jonsan/miniforge3/envs/R"
  container params.docker_container_cmd

  publishDir {"${params.outdir}/" }, mode: 'copy', pattern: "*.{biom}"
  
  input:
  val(study)

  output:
  tuple val(study), path('*.biom'), emit: gf_biom

  script:
  study = task.ext.study ?: "${study}"
  """
  >&2 R --version 
  mkdir -p ${params.hub_cache_path}
  export EXPERIMENT_HUB_CACHE='${params.hub_cache_path}'
  >&2 echo "exported"
  download_cmd.R ${study}

  convert_cmd.py data.mtx rows.txt cols.txt ${study}.biom
  """

  stub:
  study = task.ext.study ?: "${study}"
  """
  touch ${study}_genefamilies.biom
  """
}


process regroup_humann_tables {
  tag "$study_$group"
  label "bigmem"
  // conda "/home/jonsan/miniforge3/envs/humann"
  container params.docker_container_humann3

  publishDir {"${params.outdir}/${group}" }, mode: 'copy', pattern: "*.{biom}"
  
  input:
  tuple val(group), val(study), path(gf_biom)

  output:
  tuple val(group), val(study), path('*.biom'), emit: gp_biom


  script:
  study = task.ext.study ?: "${study}"
  group = task.ext.group ?: "${group}"
  """
  safe_regroup.py \
    ${gf_biom} \
    ${group} \
    ${study}_${group}.biom \
    ${params.split_size}

  """

  stub:
  study = task.ext.study ?: "${study}"
  group = task.ext.group ?: "${group}"
  """
  touch ${study}_${group}.biom
  """
}



process split_humann_tables {
  tag "$study_$group"
  label "bigmem"
  // conda "/home/jonsan/miniforge3/envs/humann"
  container params.docker_container_humann3

  publishDir {"${params.outdir}/${group}/split" }, mode: 'copy', pattern: "*.{biom}"
  
  input:
  tuple val(group), val(study), path(biom)

  output:
  path('*_stratified.biom'), emit: stratified_biom
  path('*_unstratified.biom'), emit: unstratified_biom

  when:
  params.split

  script:
  study = task.ext.study ?: "${study}"
  group = task.ext.group ?: "${group}"
  """
  humann_split_stratified_table \
    -i ${biom} \
    -o output

  mv output/*.biom .
  """

  stub:
  study = task.ext.study ?: "${study}"
  group = task.ext.group ?: "${group}"
  """
  mkdir output
  touch output/${study}_${group}_stratified.biom
  touch output/${study}_${group}_unstratified.biom

  mv output/*.biom .
  """
}




workflow {
  studies = channel
    .fromPath(params.input)
    .splitText() { it.trim() }

  // studies.view()

  groups = channel
    .fromList(params.groups.split(',') as List)
  // groups.view()

  download_cmd_tables(studies)
  gf_bioms = download_cmd_tables.out.gf_biom

  // count samples
  // split channel to small and large
  // 
  // gf_bioms.view()
  groups
    .combine(gf_bioms) | regroup_humann_tables | split_humann_tables
    

}


