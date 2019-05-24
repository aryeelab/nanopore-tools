workflow combine_sample_sheets {
    call rbind
}

task rbind {
    Array[File] sample_sheets

    command <<<
        cat ${sep=' ' sample_sheets} | awk '!seen[$0]++' > samples.csv
    >>>
    
    output {
        File samples = "samples.csv"    
    }

}