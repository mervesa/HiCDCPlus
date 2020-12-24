test_functions <- function() {
    checkEquals(get_chrs("Hsapiens","hg19"), paste0("chr",c(1:22,"X")))
}
