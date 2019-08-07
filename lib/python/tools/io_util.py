
def combine_csv( input_csvs, output_csv, header_lines = 1 ):
    ''' Combine a list of CSV files specified in input_csvs into
        a single output csv in output_csv. It is assumed that all
        CSVs have the same header structure. The number of header
        lines is specified in header_lines
    '''
    with open(output_csv, "w") as out:
        for i, icsv in enumerate(input_csvs):
            with open( icsv, "r" ) as infile:
                header = []
                for h in xrange(header_lines):
                    header.append(infile.next( ))
                if i == 0:
                    for hline in header:
                        out.write(hline)
                for line in infile:
                    out.write(line)
