library(TwoSampleMR)

mrTwoSampleMR = function(in_file_exposure, in_file_outcome, out_file, harmonise_action = 1) {
    df_exposure = read_exposure_data(in_file_exposure, sep="\t")
    df_outcome = read_outcome_data(in_file_outcome, sep="\t")
    data = harmonise_data(df_exposure, df_outcome, action = harmonise_action)
    df = mr(data)
    write.table(df, file=out_file, sep="\t", row.names=FALSE, quote=FALSE)
}

args = commandArgs(trailingOnly=TRUE)
in_file_exposure = args[1]
in_file_outcome = args[2]
out_file = args[3]
harmonise_action = if (length(args) > 3) as.integer(args[4]) else 1
mrTwoSampleMR(in_file_exposure, in_file_outcome, out_file, harmonise_action)
