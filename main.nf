#!/usr/bin/env nextflow


def runDir = WorkdirManager.resolveRunDir(workflow.launchDir, workflow.runName)

params.run_dir = runDir.toString()
params.input = ""
params.mode = ""
params.sampleName = ""
params.cohort = ""
params.r1 = ""
params.r2 = ""
params.patientId = ""


if (params.mode == 'single') {

    if (!params.input) {
        error "Single-sample mode requires --input <samplesheet.csv>."
    }
    sample_ch = Channel.fromPath(params.input)
        .splitCsv(header: true, sep: ",")
        .map { row ->
            def meta = [sampleName: row.sample_name, pairedEnd: row.paired]
            def r1 = row.file_path_R1
            def r2 = row.file_path_R2
            return [meta, r1, r2]
        }


} else if (params.mode == 'multi') {
    if (!params.cohort) {
        error "Parameter --cohort is required for multi-sample mode."
    }
    if (!params.input) {
        error "Multi-sample mode requires --input <samplesheet.csv>."
    }
    sample_ch = Channel.fromPath(params.input)
        .splitCsv(header: true, sep: ",")
        .map { row ->
            def meta = [sampleName: row.sample_name, pairedEnd: row.paired]
            if (row.patient_id) {
                meta.patientId = row.patient_id
            }
            def r1 = row.file_path_R1
            def r2 = row.file_path_R2
            return [meta, r1, r2]
        }

} else if (params.mode == 'multi_patient') {
    if (!params.input) {
        error "Multi-patient mode requires --input <samplesheet.csv> with patient_id column."
    }

    sample_ch = Channel.fromPath(params.input)
        .splitCsv(header: true, sep: ",")
        .map { row ->
            if (!row.patient_id) {
                error "Multi-patient CSV must contain a patient_id column."
            }
            def meta = [sampleName: row.sample_name, patientId: row.patient_id, pairedEnd: row.paired]
            def r1 = row.file_path_R1
            def r2 = row.file_path_R2
            return [meta, r1, r2]
        }

} else {
    error "Unknown params.mode value '${params.mode}'. Use 'single', 'multi', or 'multi_patient'."
}


include {fasterq_dump_DownloadAndSplit} from "./modules/download_samples.nf"
include {trimgalore_trim_reads} from "./modules/trimgalore"
include {bwa_index; bwa_align_reads} from "./modules/bwa"
include {gatk_mark_duplicates; gatk_base_recalibrator; gatk_applybqsr} from "./modules/gatk_preprocessing"
include {gatk_HaplotypeCaller; gatk_combine_gvcfs; gatk_GenotypeGVCFs; gatk_select_variants_SNPs; gatk_select_variants_INDELs; gatk_FilterVariants_SNPs; gatk_FilterVariants_INDELs} from "./modules/gatk_variant_calling"
include {gatk_extract_filtered_variants;DOWNLOAD_SNPEFF_DB;SNPEFF_ANNOTATE; clinvar_download; clinvar_annotation; gnomad_download; gnomad_annotation; vep_annotation; bcftools_filter} from "./modules/gatk_Filter&Annotation"


workflow {
    if (params.mode == 'single') {
        single_sample()
    } else if (params.mode == 'multi') {
        multi_samples()
    } else if (params.mode == 'multi_patient') {
        multi_patient_samples()
    }
}


workflow common_preprocessing {
    take:
    sample_ch

    main:
    trimmed_reads = trimgalore_trim_reads(sample_ch)
    bwa_index(params.ref)
    bwa_align_reads(trimgalore_trim_reads.out.trimmed_fastq)
    gatk_mark_duplicates(bwa_align_reads.out)
    gatk_base_recalibrator(gatk_mark_duplicates.out.dedup_bam)

    applybqsr_input = gatk_mark_duplicates.out.dedup_bam.join(gatk_base_recalibrator.out)
    gatk_applybqsr(applybqsr_input)
    gatk_HaplotypeCaller(gatk_applybqsr.out.dedup_bqsr_bam, gatk_applybqsr.out.dedup_bqsr_index)

    emit:
    gvcf = gatk_HaplotypeCaller.out.gvcf
    gvcf_index = gatk_HaplotypeCaller.out.gvcf_index
}


workflow single_sample {
    main:
    gvcf_output = common_preprocessing(sample_ch)

    gvcf_ch = gvcf_output.gvcf
    gvcf_index_ch = gvcf_output.gvcf_index

    gatk_GenotypeGVCFs(gvcf_ch, gvcf_index_ch)
    gatk_select_variants_SNPs(gatk_GenotypeGVCFs.out.vcf, gatk_GenotypeGVCFs.out.vcf_index)
    gatk_select_variants_INDELs(gatk_GenotypeGVCFs.out.vcf, gatk_GenotypeGVCFs.out.vcf_index)
    gatk_FilterVariants_SNPs(gatk_select_variants_SNPs.out.snp_vcf, gatk_select_variants_SNPs.out.snp_vcf_index)
    gatk_FilterVariants_INDELs(gatk_select_variants_INDELs.out.indel_vcf, gatk_select_variants_INDELs.out.indel_vcf_index)

    snp_ch = gatk_FilterVariants_SNPs.out.filtered_snp.join(gatk_FilterVariants_SNPs.out.filtered_snp_index)
    indel_ch = gatk_FilterVariants_INDELs.out.filtered_indels.join(gatk_FilterVariants_INDELs.out.filtered_indels_index)

    // Extract SNPs and INDELs separately so both are annotated independently
    snp_ch = gatk_FilterVariants_SNPs.out.filtered_snp.join(
        gatk_FilterVariants_SNPs.out.filtered_snp_index
    )

    indel_ch = gatk_FilterVariants_INDELs.out.filtered_indels.join(
        gatk_FilterVariants_INDELs.out.filtered_indels_index
    )


    // Merge SNP and INDEL channels
    snp_indel_ch = snp_ch.mix(indel_ch)


    // Run extraction ONCE for both file types
    extracted_ch = gatk_extract_filtered_variants(snp_indel_ch)


    // Download snpEff database once
    DOWNLOAD_SNPEFF_DB()


    // Annotate both SNP and INDEL extracted VCFs
    ann_ch = SNPEFF_ANNOTATE(
        gatk_extract_filtered_variants.out.extracted_filtered_variants,
        DOWNLOAD_SNPEFF_DB.out.snpeff_db_path
    )


    // Download ClinVar once
    clinvar_download()


    // Annotate both SNP and INDEL VCFs with ClinVar
    clin_ch = clinvar_annotation(
        SNPEFF_ANNOTATE.out.ann_vcf,
        SNPEFF_ANNOTATE.out.ann_vcf_index,
        clinvar_download.out.clinvar_vcf,
        clinvar_download.out.clinvar_vcf_index
    )


    // Download gnomAD once
    gnomad_download()


    // Annotate both SNP and INDEL VCFs with gnomAD
    gnomad_annotation(
        clinvar_annotation.out.clinvar_vcf,
        clinvar_annotation.out.clinvar_vcf_index,
        gnomad_download.out.gnomad_vcf,
        gnomad_download.out.gnomad_vcf_index
    )

    vep_annotation(
        gnomad_annotation.out.gnomad_vcf,
        gnomad_annotation.out.gnomad_vcf_index
    )

    bcftools_filter(
        vep_annotation.out.vep_vcf,
        vep_annotation.out.vep_vcf_index
    )


    //emit:
    //annotated_variants = SNPEFF_ANNOTATE.out.ann_vcf
}


workflow multi_samples {
    main:
    gvcf_output = common_preprocessing(sample_ch)

    vcf_list = gvcf_output.gvcf.map { meta, gvcf -> gvcf }.collect()
    tbi_list = gvcf_output.gvcf_index.map { meta, gvcf_index -> gvcf_index }.collect()

    cohort_metadata = [cohortName: params.cohort]
    gatk_combine_gvcfs(vcf_list.map { [cohort_metadata, it] }, tbi_list.map { [cohort_metadata, it] })

    gatk_GenotypeGVCFs(gatk_combine_gvcfs.out.vcf, gatk_combine_gvcfs.out.vcf_index)
    gatk_select_variants_SNPs(gatk_GenotypeGVCFs.out.vcf, gatk_GenotypeGVCFs.out.vcf_index)
    gatk_select_variants_INDELs(gatk_GenotypeGVCFs.out.vcf, gatk_GenotypeGVCFs.out.vcf_index)
    gatk_FilterVariants_SNPs(gatk_select_variants_SNPs.out.snp_vcf, gatk_select_variants_SNPs.out.snp_vcf_index)
    gatk_FilterVariants_INDELs(gatk_select_variants_INDELs.out.indel_vcf, gatk_select_variants_INDELs.out.indel_vcf_index)

    snp_ch = gatk_FilterVariants_SNPs.out.filtered_snp.join(
        gatk_FilterVariants_SNPs.out.filtered_snp_index
    )

    indel_ch = gatk_FilterVariants_INDELs.out.filtered_indels.join(
        gatk_FilterVariants_INDELs.out.filtered_indels_index
    )


    // Merge SNP and INDEL channels
    snp_indel_ch = snp_ch.mix(indel_ch)


    // Run extraction ONCE for both file types
    extracted_ch = gatk_extract_filtered_variants(snp_indel_ch)


    // Download snpEff database once
    DOWNLOAD_SNPEFF_DB()


    // Annotate both SNP and INDEL extracted VCFs
    ann_ch = SNPEFF_ANNOTATE(
        gatk_extract_filtered_variants.out.extracted_filtered_variants,
        DOWNLOAD_SNPEFF_DB.out.snpeff_db_path
    )


    // Download ClinVar once
    clinvar_download()


    // Annotate both SNP and INDEL VCFs with ClinVar
    clin_ch = clinvar_annotation(
        SNPEFF_ANNOTATE.out.ann_vcf,
        SNPEFF_ANNOTATE.out.ann_vcf_index,
        clinvar_download.out.clinvar_vcf,
        clinvar_download.out.clinvar_vcf_index
    )


    // Download gnomAD once
    gnomad_download()


    // Annotate both SNP and INDEL VCFs with gnomAD
    gnomad_annotation(
        clinvar_annotation.out.clinvar_vcf,
        clinvar_annotation.out.clinvar_vcf_index,
        gnomad_download.out.gnomad_vcf,
        gnomad_download.out.gnomad_vcf_index
    )

    vep_annotation(
        gnomad_annotation.out.gnomad_vcf,
        gnomad_annotation.out.gnomad_vcf_index
    )

    bcftools_filter(
        vep_annotation.out.vep_vcf,
        vep_annotation.out.vep_vcf_index
    )

    //emit:
    //annotated_variants = gnomad_annotation.out.ann_vcf
}


workflow multi_patient_samples {
    main:
    // Preprocess all samples
    gvcf_output = common_preprocessing(sample_ch)

    // Group GVCFs by patient ID
    patient_vcf_ch = gvcf_output.gvcf
        .map { meta, gvcf -> [meta.sampleName,meta.patientId, gvcf] }
        .groupTuple(by: 1)
        .map { sampleName, patientId, gvcfs ->
            [ [sampleName: sampleName, patientId: patientId], gvcfs ]
        }

    // Group GVCF indices by patient ID
    patient_tbi_ch = gvcf_output.gvcf_index
        .map { meta, gvcf_index -> [meta.sampleName,meta.patientId, gvcf_index] }
        .groupTuple(by: 1)
        .map { sampleName, patientId, gvcf_indices ->
            [ [sampleName: sampleName, patientId: patientId], gvcf_indices ]
        }

    // Combine GVCFs for each patient separately
    gatk_combine_gvcfs(patient_vcf_ch, patient_tbi_ch)

    // Genotype each patient's combined GVCF using patient metadata
    gatk_GenotypeGVCFs(gatk_combine_gvcfs.out.vcf, gatk_combine_gvcfs.out.vcf_index)
    gatk_select_variants_SNPs(gatk_GenotypeGVCFs.out.vcf, gatk_GenotypeGVCFs.out.vcf_index)
    gatk_select_variants_INDELs(gatk_GenotypeGVCFs.out.vcf, gatk_GenotypeGVCFs.out.vcf_index)
    gatk_FilterVariants_SNPs(gatk_select_variants_SNPs.out.snp_vcf, gatk_select_variants_SNPs.out.snp_vcf_index)
    gatk_FilterVariants_INDELs(gatk_select_variants_INDELs.out.indel_vcf, gatk_select_variants_INDELs.out.indel_vcf_index)


    snp_ch = gatk_FilterVariants_SNPs.out.filtered_snp.join(
        gatk_FilterVariants_SNPs.out.filtered_snp_index
    )

    indel_ch = gatk_FilterVariants_INDELs.out.filtered_indels.join(
        gatk_FilterVariants_INDELs.out.filtered_indels_index
    )


    // Merge SNP and INDEL channels
    snp_indel_ch = snp_ch.mix(indel_ch)


    // Run extraction ONCE for both file types
    extracted_ch = gatk_extract_filtered_variants(snp_indel_ch)


    // Download snpEff database once
    DOWNLOAD_SNPEFF_DB()


    // Annotate both SNP and INDEL extracted VCFs
    ann_ch = SNPEFF_ANNOTATE(
        extracted_ch.extracted_filtered_variants,
        DOWNLOAD_SNPEFF_DB.out.snpeff_db_path
    )


    // Download ClinVar once
    clinvar_download()


    // Annotate both SNP and INDEL VCFs with ClinVar
    clin_ch = clinvar_annotation(
        ann_ch.ann_vcf,
        ann_ch.ann_vcf_index,
        clinvar_download.out.clinvar_vcf,
        clinvar_download.out.clinvar_vcf_index
    )


    // Download gnomAD once
    gnomad_download()


    // Annotate both SNP and INDEL VCFs with gnomAD

    gnomad_annotation(
        clin_ch.clinvar_vcf,
        clin_ch.clinvar_vcf_index,
        gnomad_download.out.gnomad_vcf,
        gnomad_download.out.gnomad_vcf_index
    )


    vep_annotation(
        gnomad_annotation.out.gnomad_vcf,
        gnomad_annotation.out.gnomad_vcf_index
    )

    bcftools_filter(
        vep_annotation.out.vep_vcf,
        vep_annotation.out.vep_vcf_index
    )

        //emit:
        //annotated_variants = gnomad_annotation.out.ann_vcf
}
