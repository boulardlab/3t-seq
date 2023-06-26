checkpoint parse_runinfo:
    input: expand(data_folder.joinpath("{serie}_sra.csv"), serie=geo_series)
    output: data_folder.joinpath("serie_sheet.txt")
    params:
        data_folder=data_folder,
        log_folder=log_folder,
        alignments_folder=alignments_folder
    log: log_folder.joinpath("download/parse_runinfo.log")
    shell:
        """
        set -x
        for runinfo in {input}; do
            # extract serie id from filename
            SERIE=$(basename $runinfo)
            SERIE=${{SERIE%_sra.csv}}
            
            # determine experiment type by looking at the LibraryStrategy column
            ETYPE=$(cat $runinfo | cut -d ',' -f 13 | tail -n +2 | uniq)            
            ETYPE=${{ETYPE,,}}
            
            # determine library layout
            LAYOUT=$(cat $runinfo | cut -d ',' -f 16 | tail -n +2 | uniq)
            LAYOUT=${{LAYOUT,,}}
            
            echo $runinfo $SERIE $ETYPE $LAYOUT >> {log}
            echo $SERIE $ETYPE $LAYOUT >> {output}
            
        done
        """