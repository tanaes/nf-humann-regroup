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


//Creates working dir
workingpath = params.outdir + "/" + params.project
workingdir = file(workingpath)
if( !workingdir.exists() ) {
  if( !workingdir.mkdirs() )  {
    exit 1, "Cannot create working directory: $workingpath"
  }
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

//Folders
summary['Folders'] = ""
summary['Output dir'] = workingpath
summary['Working dir'] = workflow.workDir
summary['Output dir'] = params.outdir
summary['Script dir'] = workflow.projectDir
summary['Lunching dir'] = workflow.launchDir

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

  container params.docker_container_cmd

  publishDir {"${params.outdir}/" }, mode: 'copy', pattern: "*.{biom}"
  
  input:
  val(study)

  output:
  tuple val(study), path('*_genefamilies.biom'), emit: gf_biom

  script:
  study = task.ext.study ?: "${study}"
  """
  export EXPERIMENT_HUB_CACHE='${params.hub_cache_path}'
  download_cmd.R ${study}

  convert_cmd.py data.mtx rows.txt cols.txt ${study}_genefamilies.biom
  """

  stub:
  study = task.ext.study ?: "${study}"
  """
  touch ${study}_genefamilies.biom
  """
}


process regroup_humann_tables {
  tag "$study"

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
  humann_config --update database_folders utility_mapping ${params.utility_db}

  humann_regroup_table \
    -i ${gf_biom} \
    -g ${group} \
    -o ${study}_${group}.biom
  """

  stub:
  study = task.ext.study ?: "${study}"
  group = task.ext.group ?: "${group}"
  """
  touch ${study}_${group}.biom
  """
}



process split_humann_tables {
  tag "$study"

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
  // gf_bioms.view()
  groups
    .combine(gf_bioms) | regroup_humann_tables | split_humann_tables
    

}


