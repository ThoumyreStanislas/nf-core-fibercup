#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    nf-core/fibercup
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/nf-core/fibercup
    Website: https://nf-co.re/fibercup
    Slack  : https://nfcore.slack.com/channels/fibercup
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

params.help = false

// Importing modules and processes

include { DENOISING_MPPCA } from "./modules/nf-scil/denoising/mppca/main.nf"
include { UTILS_EXTRACTB0 } from "../modules/nf-scil/utils/extractb0/main.nf"
include { BETCROP_FSLBETCROP } from "./modules/nf-scil/betcrop/fslbetcrop/main.nf"
include { PREPROC_N4 } from "./modules/nf-scil/preproc/n4/main.nf"
include { RECONST_DTIMETRICS } from "./modules/nf-scil/reconst/dtimetrics/main.nf"
include { RECONST_FRF } from "./modules/nf-scil/reconst/frf/main.nf"
include { RECONST_FODF } from "./modules/nf-scil/reconst/fodf/main.nf"
include { TRACKING_LOCALTRACKING } from "./modules/nf-scil/tracking/localtracking/main.nf"

workflow {

    main:

        dwi_channel = Channel.fromFilePairs("$input/**/*dwi.nii.gz", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        bval_channel = Channel.fromFilePairs("$input/**/*bval", size: 1, flat: true)
            { fetch_id(it.parent, input) }
        bvec_channel = Channel.fromFilePairs("$input/**/*bvec", size: 1, flat: true)
            { fetch_id(it.parent, input) }

        // ** Denoising ** //
        DENOISING_MPPCA(dwi_channel)

         // ** Bet ** //
        bet_channel = DENOISING_MPPCA.out.dwi
            .combine(bval_channel)
            .combine(bvec_channel)
        BETCROP_FSLBETCROP(bet_channel)

        // ** Extract b0 ** //
        b0_channel = BETCROP_FSLBETCROP.out.dwi
            .combine(bval_channel)
            .combine(bvec_channel)
        UTILS_EXTRACTB0(b0_channel)

        // ** N4 ** //
        n4_channel = BETCROP_FSLBETCROP.out.dwi
            .combine(UTILS_EXTRACTB0.out.b0)
            .combine(BETCROP_FSLBETCROP.out.mask)
        PREPROC_N4(n4_channel)

        // ** DTI ** //
        dti_channel = PREPROC_N4.out.dwi
            .combine(bval_channel)
            .combine(bvec_channel)
        RECONST_DTIMETRICS(dti_channel)

        // ** FRF ** //
        frf_channel = PREPROC_N4.out.dwi
            .combine(bval_channel)
            .combine(bvec_channel)
            .combine(b0_mask_channel)
        RECONST_FRF(frf_channel)

        // ** FODF ** //
        fodf_channel = PREPROC_N4.out.dwi
            .combine(bval_channel)
            .combine(bvec_channel)
            .combine(b0_mask_channel)
            .combine(RECONST_DTIMETRICS.out.fa)
            .combine(RECONST_DTIMETRICS.out.md)
            .combine(RECONST_FRF.out.frf)
        RECONST_FODF(fodf_channel)

        // ** Local Tracking ** //
        tracking_channel = RECONST_FODF.out.fodf
            .combine(tracking_mask_channel)
            .combine(seed_channel)
        TRACKING_LOCALTRACKING(tracking_channel)

}

